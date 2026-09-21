"""
Metadata-cube interpolation tests for the NISAR GDAL driver, driven by earthaccess.

Granules are located through NASA Earthdata (CMR) and read through the driver over
HTTPS (CloudFront) with `.netrc` credentials, or over S3 when running inside AWS
us-west-2. The DEM is the public NISAR DEM VRT.

Run with pytest (or directly with `python test_interpolation_earthaccess.py`):

    earthaccess login --strategy netrc        # or place ~/.netrc first
    pytest -v conda-build/tests/test_interpolation_earthaccess.py

Environment overrides:
    NISAR_TEST_GRANULE   granule ID stem shared by the RSLC/GCOV/GSLC triple
                         (default: a 2025-11-09 DHDH acquisition, ~1 GB each)
    NISAR_TEST_DATA_DIR  directory with local copies of the granules; when the
                         files exist there they are used instead of remote reads
    NISAR_TEST_DEM       DEM path/URL (default: public EPSG4326.vrt over HTTPS)
    NISAR_TEST_ACCESS    "external" (HTTPS, default) or "direct" (S3, in-region)
"""

import os
import re
import subprocess
import sys

import numpy as np
import pytest

try:
    import earthaccess
except ImportError:  # pragma: no cover
    earthaccess = None

from osgeo import gdal

gdal.UseExceptions()

DEFAULT_STEM = (
    "004_159_A_024_2005_DHDH_A_20251109T051646_20251109T051650_X05010_N_P_J_001"
)
DEFAULT_DEM = (
    "/vsicurl/https://nisar.asf.earthdatacloud.nasa.gov/NISAR/DEM/v1.2/"
    "EPSG4326/EPSG4326.vrt"
)
SHORT_NAMES = {
    "RSLC": "NISAR_L1_RSLC_BETA_V1",
    "GCOV": "NISAR_L2_GCOV_BETA_V1",
    "GSLC": "NISAR_L2_GSLC_BETA_V1",
}
QUANTITY = "incidenceAngle"
WINDOW_M = 10000.0  # side of the test window, metres


def _in_us_west_2(target="s3.us-west-2.amazonaws.com", threshold_ms=1.0):
    try:
        proc = subprocess.run(
            ["ping", "-c", "3", "-i", "0.2", target],
            capture_output=True, text=True, timeout=5,
        )
        m = re.search(r"(\d+\.\d+)/(\d+\.\d+)/(\d+\.\d+)/", proc.stdout)
        return proc.returncode == 0 and m is not None and float(m.group(2)) < threshold_ms
    except Exception:
        return False


def _configure_gdal(access):
    gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")
    gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    gdal.SetConfigOption("GDAL_CACHEMAX", "536870912")
    if access == "direct":
        gdal.SetConfigOption("AWS_REGION", "us-west-2")
    else:
        gdal.SetConfigOption("GDAL_ENABLE_HDF5_ROS3", "NO")
        gdal.SetConfigOption("GDAL_HTTP_NETRC", "YES")
        gdal.SetConfigOption("GDAL_HTTP_UNRESTRICTED_AUTH", "YES")
        gdal.SetConfigOption("GDAL_HTTP_TIMEOUT", "60")
        gdal.SetConfigOption("GDAL_HTTP_MAX_RETRY", "5")
        gdal.SetConfigOption("GDAL_HTTP_VERSION", "2")
        gdal.SetConfigOption("GDAL_HTTP_MULTIPLEX", "YES")
        gdal.SetConfigOption("CPL_VSIL_CURL_CHUNK_SIZE", "8388608")
        gdal.SetConfigOption("CPL_VSIL_CURL_CACHE_SIZE", "536870912")
        gdal.SetConfigOption("GDAL_INGESTED_BYTES_AT_OPEN", "32768")
        gdal.SetConfigOption("GDAL_HTTP_MERGE_CONSECUTIVE_RANGES", "YES")
        cookies = os.path.join(os.path.expanduser("~"), ".urs_cookies")
        gdal.SetConfigOption("GDAL_HTTP_COOKIEFILE", cookies)
        gdal.SetConfigOption("GDAL_HTTP_COOKIEJAR", cookies)


class Granules:
    """Resolves product -> URI for one RSLC/GCOV/GSLC acquisition."""

    def __init__(self):
        self.stem = os.environ.get("NISAR_TEST_GRANULE", DEFAULT_STEM)
        self.local_dir = os.environ.get("NISAR_TEST_DATA_DIR")
        self.dem = os.environ.get("NISAR_TEST_DEM", DEFAULT_DEM)
        self.access = os.environ.get(
            "NISAR_TEST_ACCESS", "direct" if _in_us_west_2() else "external"
        )
        _configure_gdal(self.access)
        self._uris = {}
        self._window = None

    def _filename(self, product):
        level = "L1" if product == "RSLC" else "L2"
        return f"NISAR_{level}_PR_{product}_{self.stem}.h5"

    def uri(self, product):
        if product in self._uris:
            return self._uris[product]
        fn = self._filename(product)
        local = os.path.join(self.local_dir, fn) if self.local_dir else None
        if local and os.path.exists(local):
            self._uris[product] = local
            return local
        if earthaccess is None:
            pytest.skip("earthaccess not installed and no local granule copy")
        earthaccess.login(strategy="netrc")
        results = earthaccess.search_data(
            short_name=SHORT_NAMES[product], granule_name=fn[:-3], count=1
        )
        if not results:
            pytest.skip(f"{fn} not found in CMR")
        links = results[0].data_links(access=self.access)
        if not links:
            pytest.skip(f"no {self.access} data link for {fn}")
        self._uris[product] = links[0]
        return links[0]

    def conn(self, product, path=None):
        s = f'NISAR:"{self.uri(product)}"'
        return s + ":" + path if path else s

    def window(self):
        """WINDOW_M square at the GCOV grid centre, snapped to 40 m so it tiles
        both the 20 m GCOV and the 10 x 5 m GSLC postings."""
        if self._window is None:
            ref = gdal.OpenEx(self.conn("GCOV"), gdal.OF_RASTER,
                              open_options=["FREQ=A", "POL=HHHH"])
            gt = ref.GetGeoTransform()
            cx = gt[0] + gt[1] * ref.RasterXSize / 2.0
            cy = gt[3] + gt[5] * ref.RasterYSize / 2.0
            ulx = gt[0] + 40.0 * round((cx - WINDOW_M / 2 - gt[0]) / 40.0)
            uly = gt[3] - 40.0 * round((gt[3] - cy - WINDOW_M / 2) / 40.0)
            self._window = [ulx, uly, ulx + WINDOW_M, uly - WINDOW_M]
        return self._window


@pytest.fixture(scope="module")
def granules():
    if gdal.GetDriverByName("NISAR") is None:
        pytest.skip("NISAR driver not registered (set GDAL_DRIVER_PATH)")
    return Granules()


def _open_interp(granules, product, path=None, **oo):
    opts = [f"QUANTITY={QUANTITY}", f"DEM_FILE={granules.dem}"]
    opts += [f"{k}={v}" for k, v in oo.items()]
    return gdal.OpenEx(granules.conn(product, path), gdal.OF_RASTER, open_options=opts)


def _read_window(granules, ds):
    return gdal.Translate("", ds, format="MEM", projWin=granules.window())


def test_driver_version():
    drv = gdal.GetDriverByName("NISAR")
    assert drv is not None
    assert drv.GetMetadataItem("DRIVER_VERSION").startswith("v0.7.0")


def test_quantity_requires_dem(granules):
    with pytest.raises(RuntimeError, match="DEM_FILE open option is REQUIRED"):
        gdal.OpenEx(granules.conn("GCOV"), gdal.OF_RASTER,
                    open_options=[f"QUANTITY={QUANTITY}"])


def test_gcov_explicit_cube_path_matches_reference_grid(granules):
    cube = f"/science/LSAR/GCOV/metadata/radarGrid/{QUANTITY}"
    ds = _open_interp(granules, "GCOV", cube)
    ref = gdal.OpenEx(granules.conn("GCOV"), gdal.OF_RASTER,
                      open_options=["FREQ=A", "POL=HHHH"])
    assert (ds.RasterXSize, ds.RasterYSize) == (ref.RasterXSize, ref.RasterYSize)
    assert ds.GetGeoTransform() == ref.GetGeoTransform()
    assert ds.GetSpatialRef().IsSame(ref.GetSpatialRef())
    assert ds.RasterCount == 1
    assert ds.GetRasterBand(1).DataType == gdal.GDT_Float32
    md = ds.GetMetadata()
    assert md["NISAR_PRODUCT_TYPE"] == "GCOV"
    assert md["NISAR_CUBE_PATH"] == cube
    assert md["NISAR_REFERENCE_GRID"] == "/science/LSAR/GCOV/grids/frequencyA/HHHH"
    assert md["NISAR_QUANTITY"] == QUANTITY


def test_gcov_auto_resolved_cube_equals_explicit(granules):
    cube = f"/science/LSAR/GCOV/metadata/radarGrid/{QUANTITY}"
    explicit = _read_window(granules, _open_interp(granules, "GCOV", cube)).ReadAsArray()
    auto_ds = _open_interp(granules, "GCOV")
    assert auto_ds.GetMetadataItem("NISAR_CUBE_PATH") == cube
    auto = _read_window(granules, auto_ds).ReadAsArray()
    gt = auto_ds.GetGeoTransform()
    assert explicit.shape == (round(WINDOW_M / -gt[5]), round(WINDOW_M / gt[1]))
    assert np.array_equal(explicit, auto, equal_nan=True)
    assert np.isfinite(explicit).mean() > 0.9
    assert 20.0 < np.nanmin(explicit) < np.nanmax(explicit) < 60.0


def test_gslc_reference_grid_default(granules):
    ds = _open_interp(granules, "GSLC")
    ref = gdal.OpenEx(granules.conn("GSLC"), gdal.OF_RASTER,
                      open_options=["FREQ=A", "POL=HH"])
    assert (ds.RasterXSize, ds.RasterYSize) == (ref.RasterXSize, ref.RasterYSize)
    assert ds.GetGeoTransform() == ref.GetGeoTransform()
    assert ds.GetSpatialRef().IsSame(ref.GetSpatialRef())
    md = ds.GetMetadata()
    assert md["NISAR_PRODUCT_TYPE"] == "GSLC"
    assert md["NISAR_CUBE_PATH"] == f"/science/LSAR/GSLC/metadata/radarGrid/{QUANTITY}"
    assert md["NISAR_REFERENCE_GRID"] == "/science/LSAR/GSLC/grids/frequencyA/HH"


def test_gslc_reference_grid_freq_b_pol(granules):
    ds = _open_interp(granules, "GSLC", FREQ="B", POL="HV")
    ref = gdal.OpenEx(granules.conn("GSLC"), gdal.OF_RASTER,
                      open_options=["FREQ=B", "POL=HV"])
    assert (ds.RasterXSize, ds.RasterYSize) == (ref.RasterXSize, ref.RasterYSize)
    assert ds.GetGeoTransform() == ref.GetGeoTransform()
    assert ds.GetMetadataItem("NISAR_REFERENCE_GRID") == "/science/LSAR/GSLC/grids/frequencyB/HV"


def test_gslc_matches_gcov_on_common_ground(granules):
    """Same acquisition, same DEM: GSLC (10x5 m) block-averaged to the GCOV posting
    (20 m) must reproduce the GCOV interpolation to well below the cube precision."""
    gslc = _read_window(granules, _open_interp(granules, "GSLC")).ReadAsArray()
    gcov = _read_window(granules, _open_interp(granules, "GCOV")).ReadAsArray()
    fy, fx = gslc.shape[0] // gcov.shape[0], gslc.shape[1] // gcov.shape[1]
    assert (fy * gcov.shape[0], fx * gcov.shape[1]) == gslc.shape
    gslc_coarse = gslc.reshape(gcov.shape[0], fy, gcov.shape[1], fx).mean(axis=(1, 3))
    assert np.nanmax(np.abs(gslc_coarse - gcov)) < 1e-3


def test_gslc_explicit_cube_path_equals_auto(granules):
    cube = f"/science/LSAR/GSLC/metadata/radarGrid/{QUANTITY}"
    explicit = _read_window(granules, _open_interp(granules, "GSLC", cube)).ReadAsArray()
    auto = _read_window(granules, _open_interp(granules, "GSLC")).ReadAsArray()
    assert np.array_equal(explicit, auto, equal_nan=True)


def test_missing_quantity_fails(granules):
    with pytest.raises(RuntimeError, match="Failed to open valid 3D coarse metadata cube"):
        gdal.OpenEx(granules.conn("GCOV"), gdal.OF_RASTER,
                    open_options=["QUANTITY=notACube", f"DEM_FILE={granules.dem}"])


def test_level1_rejected(granules):
    with pytest.raises(RuntimeError, match="Level-1 product RSLC is not supported yet"):
        _open_interp(granules, "RSLC")


if __name__ == "__main__":
    sys.exit(pytest.main(["-v", __file__] + sys.argv[1:]))
