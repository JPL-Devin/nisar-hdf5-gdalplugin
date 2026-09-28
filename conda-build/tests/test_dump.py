"""
NISAR_DUMP metadata-domain tests for the NISAR GDAL driver.

Exercises `-oo DUMP=YES` (h5dump-style listing exposed through the NISAR_DUMP
metadata domain) in container mode on one granule per product level:
L1 RSLC, L2 GUNW, L2 GCOV and L3 SME2.

Granules are located through NASA Earthdata (CMR) with earthaccess and read
through the driver over HTTPS with `.netrc` credentials, unless a local copy
is found in NISAR_TEST_DATA_DIR.

    pytest -v -p no:cacheprovider conda-build/tests/test_dump.py

Environment overrides:
    NISAR_TEST_DATA_DIR  directory with local NISAR_L?_PR_<PRODUCT>_*.h5 files;
                         the first match per product is used instead of CMR
    NISAR_TEST_DUMP_<PRODUCT>  explicit path/URL for one product (RSLC, GUNW,
                         GCOV, SME2), bypassing both the local dir and CMR
    NISAR_TEST_ACCESS    "external" (HTTPS, default) or "direct" (S3, in-region)
"""

import glob
import os

import pytest

try:
    import earthaccess
except ImportError:  # pragma: no cover
    earthaccess = None

from osgeo import gdal

gdal.UseExceptions()

PRODUCTS = ("RSLC", "GUNW", "GCOV", "SME2")
SHORT_NAMES = {
    "RSLC": "NISAR_L1_RSLC_BETA_V1",
    "GUNW": "NISAR_L2_GUNW_BETA_V1",
    "GCOV": "NISAR_L2_GCOV_BETA_V1",
    "SME2": "NISAR_L3_SME2_BETA_V1",
}
LEVEL = {"RSLC": "L1", "GUNW": "L2", "GCOV": "L2", "SME2": "L3"}
IDENT = "/science/LSAR/identification"


def _configure_gdal():
    gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")
    gdal.SetConfigOption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    if os.environ.get("NISAR_TEST_ACCESS", "external") == "direct":
        gdal.SetConfigOption("AWS_REGION", "us-west-2")
    else:
        gdal.SetConfigOption("GDAL_HTTP_NETRC", "YES")
        gdal.SetConfigOption("GDAL_HTTP_UNRESTRICTED_AUTH", "YES")
        gdal.SetConfigOption("GDAL_HTTP_TIMEOUT", "60")
        gdal.SetConfigOption("GDAL_HTTP_MAX_RETRY", "5")
        cookies = os.path.join(os.path.expanduser("~"), ".urs_cookies")
        gdal.SetConfigOption("GDAL_HTTP_COOKIEFILE", cookies)
        gdal.SetConfigOption("GDAL_HTTP_COOKIEJAR", cookies)


class Granules:
    def __init__(self):
        self.local_dir = os.environ.get("NISAR_TEST_DATA_DIR")
        self.access = os.environ.get("NISAR_TEST_ACCESS", "external")
        _configure_gdal()
        self._uris = {}

    def uri(self, product):
        if product in self._uris:
            return self._uris[product]
        explicit = os.environ.get(f"NISAR_TEST_DUMP_{product}")
        if explicit:
            self._uris[product] = explicit
            return explicit
        if self.local_dir:
            pattern = os.path.join(
                self.local_dir, f"NISAR_{LEVEL[product]}_PR_{product}_*.h5")
            found = sorted(glob.glob(pattern))
            if found:
                self._uris[product] = found[0]
                return found[0]
        if earthaccess is None:
            pytest.skip("earthaccess not installed and no local granule copy")
        try:
            earthaccess.login(strategy="netrc")
            results = earthaccess.search_data(
                short_name=SHORT_NAMES[product], count=1)
        except Exception as exc:  # credentials / network
            pytest.skip(f"earthaccess unavailable: {exc}")
        if not results:
            pytest.skip(f"no {product} granule found in CMR")
        links = [l for l in results[0].data_links(access=self.access)
                 if l.endswith(".h5")]
        if not links:
            pytest.skip(f"no {self.access} .h5 link for {product}")
        self._uris[product] = links[0]
        return links[0]

    def conn(self, product, path=None):
        s = f'NISAR:"{self.uri(product)}"'
        return s + ":" + path if path else s


@pytest.fixture(scope="module")
def granules():
    if gdal.GetDriverByName("NISAR") is None:
        pytest.skip("NISAR driver not available (check GDAL_DRIVER_PATH)")
    return Granules()


def _open(granules, product, **oo):
    opts = [f"{k}={v}" for k, v in oo.items()]
    return gdal.OpenEx(granules.conn(product), gdal.OF_RASTER, open_options=opts)


def _dump(ds):
    lines = ds.GetMetadata_List("NISAR_DUMP")
    assert lines, "NISAR_DUMP domain is empty"
    return lines


def _objects(lines, kind):
    """Paths of GROUP/DATASET headers in a dump."""
    prefix = f'{kind} "'
    return [l[len(prefix):l.index('"', len(prefix))]
            for l in lines if l.startswith(prefix)]


# --------------------------------------------------------------------------
# Default behaviour is untouched
# --------------------------------------------------------------------------

def test_open_options_registered():
    drv = gdal.GetDriverByName("NISAR")
    if drv is None:
        pytest.skip("NISAR driver not available")
    xml = drv.GetMetadataItem("DMD_OPENOPTIONLIST")
    for name in ("DUMP", "DUMP_ROOT", "DUMP_MODE"):
        assert f"name='{name}'" in xml


@pytest.mark.parametrize("product", PRODUCTS)
def test_default_without_dump(granules, product):
    ds = _open(granules, product)
    assert "NISAR_DUMP" not in ds.GetMetadataDomainList()
    with gdal.quiet_errors():
        assert not ds.GetMetadata_List("NISAR_DUMP")
    sub = ds.GetMetadata("SUBDATASETS")
    descs = [v for k, v in sub.items() if k.endswith("_DESC")]
    assert descs
    assert not any("not openable" in d for d in descs)
    assert not any("(string" in d for d in descs)


# --------------------------------------------------------------------------
# HEADER dump
# --------------------------------------------------------------------------

@pytest.mark.parametrize("product", PRODUCTS)
def test_header_dump(granules, product):
    ds = _open(granules, product, DUMP="YES")
    assert "NISAR_DUMP" in ds.GetMetadataDomainList()
    lines = _dump(ds)

    assert lines[0].startswith('HDF5 "') and lines[0].endswith(" HEADER {")
    assert lines[-1] == "}"

    groups = _objects(lines, "GROUP")
    datasets = _objects(lines, "DATASET")
    assert "/science/LSAR" in groups
    assert IDENT in groups
    assert f"/science/LSAR/{product}" in groups
    assert f"{IDENT}/productType" in datasets
    assert all(p.startswith("/science/LSAR") for p in groups + datasets)

    # datatype / dataspace / attribute lines are present
    assert any(l.strip().startswith("DATATYPE  H5T_STRING") for l in lines)
    assert any(l.strip().startswith("DATATYPE  H5T_IEEE_F") for l in lines)
    assert any(l.strip().startswith("DATASPACE  SCALAR") for l in lines)
    assert any(l.strip().startswith("DATASPACE  SIMPLE") for l in lines)
    assert any(l.strip().startswith('ATTRIBUTE "description"') for l in lines)

    # HEADER mode prints attribute values but no dataset DATA blocks
    idx = lines.index(f'DATASET "{IDENT}/productType" {{')
    block = lines[idx:idx + 12]
    assert not any(l.startswith("   DATA {") for l in block)
    assert any(l.startswith("      DATA {") for l in block)  # attribute value


@pytest.mark.parametrize("product", PRODUCTS)
def test_dump_lists_string_subdatasets(granules, product):
    ds = _open(granules, product, DUMP="YES")
    sub = ds.GetMetadata("SUBDATASETS")
    descs = [v for k, v in sub.items() if k.endswith("_DESC")]
    assert any(f"{IDENT}/productType (string, not openable)" in d for d in descs)
    # rasters keep their normal description
    assert any(d.endswith("(Float32)") or d.endswith("(complex, Float32)")
               for d in descs)


@pytest.mark.parametrize("product", PRODUCTS)
def test_dump_root_scoping(granules, product):
    ds = _open(granules, product, DUMP="YES", DUMP_ROOT=IDENT)
    lines = _dump(ds)
    groups = _objects(lines, "GROUP")
    datasets = _objects(lines, "DATASET")
    assert groups == [IDENT]
    assert datasets
    assert all(p.startswith(IDENT + "/") for p in datasets)
    assert "/science/LSAR" not in groups


def test_dump_root_file_root(granules):
    ds = _open(granules, "GCOV", DUMP="YES", DUMP_ROOT="/")
    groups = _objects(_dump(ds), "GROUP")
    assert groups[0] == "/"
    assert "/science" in groups
    assert "/science/LSAR" in groups


def test_dump_root_missing(granules):
    with pytest.raises(RuntimeError):
        _open(granules, "GCOV", DUMP="YES", DUMP_ROOT="/no/such/group")


def test_dump_mode_invalid(granules):
    with pytest.raises(RuntimeError):
        _open(granules, "GCOV", DUMP="YES", DUMP_MODE="BOGUS")


# --------------------------------------------------------------------------
# FULL dump
# --------------------------------------------------------------------------

def test_full_dump_values(granules):
    ds = _open(granules, "GCOV", DUMP="YES", DUMP_MODE="FULL", DUMP_ROOT=IDENT)
    lines = _dump(ds)
    assert lines[0].endswith(" FULL {")

    idx = lines.index(f'DATASET "{IDENT}/productType" {{')
    block = lines[idx:idx + 4]
    assert '   DATA { "GCOV" }' in block

    idx = lines.index(f'DATASET "{IDENT}/listOfFrequencies" {{')
    data = [l for l in lines[idx:idx + 6] if l.startswith("   DATA {")]
    assert data and '"A"' in data[0]


def test_full_dump_skips_large_arrays(granules):
    root = "/science/LSAR/GCOV/grids/frequencyA"
    ds = _open(granules, "GCOV", DUMP="YES", DUMP_MODE="FULL", DUMP_ROOT=root)
    lines = _dump(ds)
    idx = lines.index(f'DATASET "{root}/xCoordinates" {{')
    data = [l for l in lines[idx:idx + 4] if l.startswith("   DATA {")]
    assert data and "elements, not printed" in data[0]
