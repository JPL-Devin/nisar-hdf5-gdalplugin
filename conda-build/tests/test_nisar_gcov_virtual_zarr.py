# Copyright 2025, by the California Institute of Technology.
# ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
#
# SPDX-License-Identifier: Apache-2.0

"""
Offline tests for nisar_gcov_virtual_zarr.py (no Earthdata access required).

A small synthetic GCOV-like HDF5 file exercises the local-file fallback: dimension-scale
coordinates, shuffle+deflate grids with unallocated chunks, a complex calibration grid, a
rank-3 metadata cube, an LZF (non-deflate) dataset, string/compound datasets and hard and
external links.

    pytest -v conda-build/tests/test_nisar_gcov_virtual_zarr.py
"""

import argparse
import csv
import json
import os
import sys
import zipfile

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("zarr")
pytest.importorskip("xarray")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nisar_gcov_virtual_zarr as nvz  # noqa: E402

GRID = "/science/LSAR/GCOV/grids/frequencyA"
CAL = "/science/LSAR/GCOV/metadata/calibrationInformation/frequencyA/elevationAntennaPattern"
CUBE = "/science/LSAR/GCOV/metadata/radarGrid"


def _scales(f, grp, ny, nx, y0=4000000.0, x0=500000.0, dx=20.0):
    y = grp.create_dataset("yCoordinates", data=y0 - dx * np.arange(ny))
    x = grp.create_dataset("xCoordinates", data=x0 + dx * np.arange(nx))
    for d, name in ((y, "projection_y_coordinate"), (x, "projection_x_coordinate")):
        d.attrs["standard_name"] = np.bytes_(name)
        d.attrs["units"] = np.bytes_("meters")
        d.make_scale()
    return y, x


def make_granule(path):
    rng = np.random.default_rng(0)
    with h5py.File(path, "w") as f:
        ident = f.create_group("/science/LSAR/identification")
        ident["granuleId"] = np.bytes_("NISAR_L2_PR_GCOV_SYNTHETIC_001")
        ident["productType"] = np.bytes_("GCOV")
        ident["zeroDopplerStartTime"] = np.bytes_("2025-10-12T18:25:08.000000")
        g = f.create_group(GRID)
        y, x = _scales(f, g, 100, 120)
        g["projection"] = np.uint32(32611)
        g["projection"].attrs["grid_mapping_name"] = np.bytes_("universal_transverse_mercator")
        hhhh = g.create_dataset("HHHH", shape=(100, 120), dtype="f4", chunks=(32, 32),
                                compression="gzip", compression_opts=1, shuffle=True, fillvalue=np.nan)
        hhhh[0:64, 0:96] = rng.random((64, 96), dtype="f4")
        hhhh.attrs["units"] = np.bytes_("1")
        hhhh.attrs["_FillValue"] = np.float32(np.nan)
        hhhh.attrs["grid_mapping"] = np.bytes_("projection")
        mask = g.create_dataset("mask", data=rng.integers(0, 4, (100, 120), dtype="u1"), chunks=(32, 32),
                                compression="gzip", compression_opts=9, fillvalue=255)
        lzf = g.create_dataset("numberOfLooks", data=rng.random((100, 120), dtype="f4"), chunks=(32, 32),
                               compression="lzf")
        for d in (hhhh, mask, lzf):
            d.dims[0].attach_scale(y)
            d.dims[1].attach_scale(x)
        g["listOfPolarizations"] = np.array([b"HH"], dtype="S2")
        g["HHHH_alias"] = hhhh  # hard link
        g["external"] = h5py.ExternalLink("does_not_exist.h5", "/data")
        c = f.create_group(CAL)
        cy, cx = _scales(f, c, 10, 12, dx=500.0)
        hh = c.create_dataset("HH", data=(rng.random((10, 12)) + 1j * rng.random((10, 12))).astype("c8"),
                              chunks=(10, 12), compression="gzip")
        hh.dims[0].attach_scale(cy)
        hh.dims[1].attach_scale(cx)
        cube = f.create_group(CUBE)
        cube.create_dataset("heightAboveEllipsoid", data=np.linspace(-500, 1500, 5))
        cube.create_dataset("incidenceAngle", data=rng.random((5, 10, 12)), chunks=(1, 10, 12),
                            compression="gzip")
        cube.create_dataset("compoundThing", data=np.zeros(3, dtype=[("a", "i4"), ("b", "f8")]))


def run_generate(tmp_path, *extra):
    src = tmp_path / "granule.h5"
    make_granule(src)
    out = tmp_path / "out"
    rc = nvz.main(["generate", "--local-file", str(src), "--output-dir", str(out), *extra])
    gid = "NISAR_L2_PR_GCOV_SYNTHETIC_001"
    return rc, src, out, out / gid, gid


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    return run_generate(tmp_path_factory.mktemp("gcov"), "--remote-uri",
                        "https://example.invalid/NISAR/granule.h5")


def _chunk_rows(outdir):
    with open(outdir / "gcov_chunk_table.csv") as f:
        return list(csv.DictReader(f))


def test_artifacts_and_zip(generated):
    rc, _src, out, outdir, gid = generated
    assert rc == 3  # validation reports the LZF dataset as not included, see test_unsupported_recorded
    zpath = out / f"nisar_gcov_virtual_zarr_{gid}.zip"
    assert zpath.exists()
    names = set(zipfile.ZipFile(zpath).namelist())
    for n in ("gcov_chunk_table.csv", "gcov_chunk_table.md", "granule_metadata.json", "environment.json",
              "validation_report.json", "validation_report.md", "unsupported_lossy_report.json",
              "unsupported_lossy_report.md", f"nisar_gcov_{gid}.kerchunk.json"):
        assert n in names


def test_manifest_targets_remote_uri_and_allocated_chunks(generated):
    _rc, src, _out, outdir, gid = generated
    refs = json.load(open(outdir / f"nisar_gcov_{gid}.kerchunk.json"))["refs"]
    assert json.loads(refs[".zgroup"])["zarr_format"] == 2
    za = json.loads(refs["grids/frequencyA/HHHH/.zarray"])
    assert za["shape"] == [100, 120] and za["chunks"] == [32, 32] and za["dtype"] == "<f4"
    assert za["compressor"]["id"] == "zlib" and za["filters"][0]["id"] == "shuffle"
    dims = json.loads(refs["grids/frequencyA/HHHH/.zattrs"])["_ARRAY_DIMENSIONS"]
    assert dims == ["yCoordinates", "xCoordinates"]
    keys = sorted(k.rsplit("/", 1)[1] for k in refs if k.startswith("grids/frequencyA/HHHH/")
                  and k.rsplit("/", 1)[1][0].isdigit())
    with h5py.File(src) as f:
        n_alloc = f[GRID + "/HHHH"].id.get_num_chunks()
    assert n_alloc == 2 * 3 and len(keys) == n_alloc  # 2x3 of 4x4 logical chunks were written
    targets = {v[0] for v in refs.values() if isinstance(v, list)}
    assert targets == {"https://example.invalid/NISAR/granule.h5"}
    assert "grids/frequencyA/xCoordinates/.zarray" in refs


def test_chunk_table(generated):
    rows = _chunk_rows(generated[3])
    by = {r["hdf5_path"]: r for r in rows}
    hhhh = by[GRID + "/HHHH"]
    assert hhhh["logical_chunks"] == "16" and hhhh["allocated_chunks"] == "6" and hhhh["skipped_chunks"] == "10"
    # Warnings: unallocated chunks and the HHHH_alias hard link.
    assert hhhh["category"] == "principal_grid" and hhhh["manifest_status"] == "included_with_warnings"
    assert by[GRID + "/mask"]["manifest_status"] == "included"
    assert by[GRID + "/mask"]["category"] == "mask"
    assert by[CAL + "/HH"]["category"] == "calibration_grid"
    assert by[CAL + "/HH"]["manifest_status"] == "included_with_warnings"
    assert by[GRID + "/numberOfLooks"]["manifest_status"] == "skipped_unsupported"
    first_total = next(i for i, r in enumerate(rows) if r["row_type"] != "dataset")
    data = rows[:first_total]
    assert [r["hdf5_path"] for r in data] == sorted(r["hdf5_path"] for r in data)
    assert len(data) == 4
    labels = [r["hdf5_path"] for r in rows if r["row_type"] != "dataset"]
    for cat in ("principal_grid", "mask", "coordinate", "calibration_grid", "metadata"):
        assert f"SUBTOTAL {cat}" in labels
    assert "TOTAL 2-D datasets" in labels


def test_unsupported_recorded(generated):
    rep = json.load(open(generated[3] / "unsupported_lossy_report.json"))
    text = json.dumps(rep)
    for needle in ("lzf", "external", "hard", "complex", "compound", "string", "not allocated"):
        assert needle in text.lower(), needle
    val = json.load(open(generated[3] / "validation_report.json"))
    assert [f for f in val["failures"] if f["scope"] == GRID + "/numberOfLooks"]
    assert not [f for f in val["failures"] if f["scope"] == GRID + "/HHHH"]
    assert val["raw_decode"]["decoded_equals_hdf5"]


def test_xarray_reads_manifest(generated):
    import xarray as xr
    _rc, src, _out, outdir, gid = generated
    doc = json.load(open(outdir / f"nisar_gcov_{gid}.kerchunk.json"))
    fo = {k: ([str(src), v[1], v[2]] if isinstance(v, list) else v) for k, v in doc["refs"].items()}
    so = {"fo": fo, "remote_protocol": "file", "asynchronous": True}
    tree = xr.open_datatree("reference://", engine="zarr", consolidated=True, zarr_format=2,
                            storage_options=so, mask_and_scale=False, chunks=None)
    da = tree["grids/frequencyA"].ds["HHHH"]
    assert da.dims == ("yCoordinates", "xCoordinates")
    assert "xCoordinates" in da.coords and "yCoordinates" in da.coords
    with h5py.File(src) as f:
        np.testing.assert_array_equal(da.values, f[GRID + "/HHHH"][()])
        np.testing.assert_array_equal(tree[CUBE.split("GCOV/")[1]].ds["incidenceAngle"].values,
                                      f[CUBE + "/incidenceAngle"][()])


def test_retarget(generated, tmp_path):
    _rc, _src, _out, outdir, gid = generated
    new = tmp_path / "m.json"
    nvz.main(["retarget", str(outdir / f"nisar_gcov_{gid}.kerchunk.json"), "s3://bucket/granule.h5",
              "-o", str(new)])
    refs = json.load(open(new))["refs"]
    assert {v[0] for v in refs.values() if isinstance(v, list)} == {"s3://bucket/granule.h5"}


def test_retarget_rejects_unmatched_old_uri(generated, tmp_path):
    _rc, _src, _out, outdir, gid = generated
    new = tmp_path / "m.json"
    rc = nvz.main(["retarget", str(outdir / f"nisar_gcov_{gid}.kerchunk.json"), "s3://bucket/granule.h5",
                   "--old-uri", "https://example.invalid/other.h5", "-o", str(new)])
    assert rc == 2
    assert not new.exists()


def test_signed_uris_rejected(generated, tmp_path):
    _rc, src, _out, outdir, gid = generated
    signed = "https://d1.cloudfront.net/g.h5?A-userid=u&Expires=1&Signature=x&Key-Pair-Id=k"
    with pytest.raises(ValueError, match="signed"):
        nvz.main(["generate", "--local-file", str(src), "--output-dir", str(tmp_path), "--remote-uri", signed])
    with pytest.raises(ValueError, match="signed"):
        nvz.main(["retarget", str(outdir / f"nisar_gcov_{gid}.kerchunk.json"), signed, "-o",
                  str(tmp_path / "m.json")])


def test_unreferenced_projection_and_cross_group_scales(tmp_path):
    src = tmp_path / "g.h5"
    grid = "/science/LSAR/GCOV/grids/frequencyB"
    with h5py.File(src, "w") as f:
        f["/science/LSAR/identification/granuleId"] = np.bytes_("NISAR_L2_PR_GCOV_SYNTHETIC_002")
        coords = f.create_group("/science/LSAR/GCOV/metadata/coords")
        cy, cx = _scales(f, coords, 4, 6, x0=10.0)
        g = f.create_group(grid)
        ly, lx = _scales(f, g, 4, 6, x0=30.0)
        g["projection"] = np.uint32(32611)
        hh = g.create_dataset("HHHH", data=np.ones((4, 6), "f4"), chunks=(2, 3), compression="gzip")
        hh.dims[0].attach_scale(cy)
        hh.dims[1].attach_scale(cx)
        mask = g.create_dataset("mask", data=np.ones((4, 6), "u1"), chunks=(2, 3))
        mask.dims[0].attach_scale(ly)
        mask.dims[1].attach_scale(lx)
    out = tmp_path / "out"
    assert nvz.main(["generate", "--local-file", str(src), "--output-dir", str(out)]) == 0
    gid = "NISAR_L2_PR_GCOV_SYNTHETIC_002"
    refs = json.load(open(out / gid / f"nisar_gcov_{gid}.kerchunk.json"))["refs"]
    assert "grids/frequencyB/projection/.zarray" in refs
    dims = json.loads(refs["grids/frequencyB/HHHH/.zattrs"])["_ARRAY_DIMENSIONS"]
    assert dims == ["metadata__coords__yCoordinates", "metadata__coords__xCoordinates"]
    assert json.loads(refs["grids/frequencyB/mask/.zattrs"])["_ARRAY_DIMENSIONS"] == ["yCoordinates", "xCoordinates"]


def test_seeded_selection_is_deterministic(monkeypatch):
    def fake(n):
        return {"umm": {"GranuleUR": f"NISAR_L2_PR_GCOV_{n:03d}", "RelatedUrls": []}, "meta": {}}
    pool = [fake(n) for n in (5, 3, 9, 1, 7)]
    monkeypatch.setattr(nvz, "earthdata_login", lambda strategy: "netrc")
    monkeypatch.setattr(nvz, "search_candidates", lambda args: (list(reversed(pool)), {}))
    args = argparse.Namespace(local_file=None, granule_ur=None, concept_id=None, seed=nvz.DEFAULT_SEED,
                              auth="auto")
    a = nvz.select_granule(args)
    monkeypatch.setattr(nvz, "search_candidates", lambda args: (pool, {}))
    b = nvz.select_granule(args)
    assert a["granule"]["granule_ur"] == b["granule"]["granule_ur"]
    assert a["selection"]["candidate_list_sha256"] == b["selection"]["candidate_list_sha256"]
