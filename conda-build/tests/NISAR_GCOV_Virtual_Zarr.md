# NISAR GCOV virtual Zarr (Kerchunk) workflow

`nisar_gcov_virtual_zarr.py` retrieves one NISAR Level-2 GCOV granule from NASA Earthdata,
inventories every HDF5 dataset under `/science/LSAR/GCOV`, and writes a consolidated
Kerchunk (reference spec v1) JSON that exposes the granule as a Zarr v2 hierarchy. Every chunk
reference points at the granule's archive URI, so no pixel data is copied: only HDF5 metadata,
the chunk index and a few validation chunks are read with HTTP range requests.

Offline tests (synthetic HDF5, no network): `test_nisar_gcov_virtual_zarr.py`.

## Requirements

The `~/nisar-env` environment plus `h5py`, `zarr` (>= 3, reading Zarr v2), `numcodecs`, `xarray`
and `fsspec`/`aiohttp` (installed with `earthaccess`):

```bash
micromamba install -y -p ~/nisar-env -c conda-forge h5py zarr numcodecs xarray
```

GDAL is **not** used; validation decodes chunks with fsspec + numcodecs + zarr/xarray only.

## Earthdata authentication

Login is attempted in this order (`--auth auto`, the default):

1. an existing authenticated `earthaccess` session (when the module is imported from Python
   after `earthaccess.login()`);
2. environment variables: `EARTHDATA_TOKEN`, or `EARTHDATA_USERNAME` + `EARTHDATA_PASSWORD`;
3. `~/.netrc` with `machine urs.earthdata.nasa.gov login <user> password <password>`
   (mode `600`).

Force one strategy with `--auth environment|netrc|interactive`. Credentials, cookies and signed
URLs are never written to the outputs.

## Granule selection

| Mode | Command |
|---|---|
| Seeded random (default seed `20251012`, October 2025) | `python nisar_gcov_virtual_zarr.py generate` |
| Other seed / window / area | `... generate --seed 7 --temporal 2025-11-01T00:00:00Z 2025-11-30T23:59:59Z --bbox -125 32 -114 42` |
| Whole collection | `... generate --all-time --max-candidates 50000` |
| Explicit concept ID | `... generate --concept-id G4062016788-ASF` |
| Explicit granule UR | `... generate --granule-ur NISAR_L2_PR_GCOV_002_109_D_063_4005_DHDH_A_20251012T182508_20251012T182531_X05010_N_P_J_001` |
| Only print the selection | add `--select-only` |

Seeded selection: CMR candidates for `--short-name` (default `NISAR_L2_GCOV_BETA_V1`) matching
`--temporal` / `--bbox` are sorted by GranuleUR and the index is
`random.Random(seed).randrange(len(candidates))`. `granule_metadata.json` records the query,
CMR hit count, candidate count, selected index and a SHA-256 of the sorted candidate list, so a
run is reproducible as long as the candidate set is unchanged. The default window is fixed in
the past, but ASF can still reprocess or add granules; compare `candidate_list_sha256` to
confirm that a rerun saw the same candidates, or use `--concept-id` / `--granule-ur` for an
exact rerun. Searches with more than `--max-candidates` (default 5000) hits are rejected.

## Remote access and URL persistence

`--access external` (default) uses the HTTPS `GET DATA` link from CMR, e.g.
`https://nisar.asf.earthdatacloud.nasa.gov/NISAR/NISAR_L2_GCOV_BETA_V1/<UR>/<UR>.h5`.
For ASF NISAR this is a TEA endpoint: an authenticated request gets `303 See Other` to a
**CloudFront signed URL** (`*.cloudfront.net/...sds-n-cumulus-prod-nisar-products.s3.us-west-2...`
with `Expires`, `Signature`, `Key-Pair-Id`, `A-userid` query parameters, observed lifetime
3600 s) which serves byte ranges (`206`) from S3. The script probes this chain and stores the
redirect host, redacted path, parameter *names*, expiry and lifetime in
`granule_metadata.json` (`access_probe`); the signature is never stored.

`--access direct` uses the `s3://sds-n-cumulus-prod-nisar-products/...` link with temporary S3
credentials from `earthaccess` (only works in AWS `us-west-2`).

The manifest targets the **stable archive URI** (the HTTPS link, or the `s3://` link with
`--access direct`), never the signed URL. Readers therefore need an Earthdata-authenticated
fsspec session and the redirect is resolved per request. The references stay valid only while
the archive keeps that URI and the object is not reprocessed in place (CMR `revision_id`, size,
MD5 checksum, `ETag` and `Last-Modified` are recorded to detect this). If the URI changes, or
you want S3 targets instead of HTTPS, rewrite every reference without regenerating:

```bash
python nisar_gcov_virtual_zarr.py retarget output/<granule_id>/nisar_gcov_<granule_id>.kerchunk.json \
    s3://sds-n-cumulus-prod-nisar-products/NISAR_L2_GCOV_BETA_V1/<UR>/<UR>.h5 -o manifest_s3.json
```

`retarget` fails without writing anything if no chunk reference targets the old URI (the
manifest's `source_uri`, or `--old-uri`). `generate --remote-uri` and `retarget` reject URLs
whose query carries signatures or tokens (`Signature`, `Key-Pair-Id`, `Expires`,
`X-Amz-Signature`, ...), so pre-signed URLs never end up in shared manifests. Always target the
stable URI.

## Local-file fallback

```bash
python nisar_gcov_virtual_zarr.py generate --local-file /data/<UR>.h5 \
    --remote-uri https://nisar.asf.earthdatacloud.nasa.gov/NISAR/NISAR_L2_GCOV_BETA_V1/<UR>/<UR>.h5
```

Metadata, chunk offsets and validation come from the local copy (no network). References still
target `--remote-uri` (or the CMR HTTPS link when `--granule-ur`/`--concept-id` is also given;
`file://<path>` if neither is available). Validation decodes from the local copy with the same
byte offsets; add `--validate-remote` to decode over the network instead. Byte offsets are only
valid for the exact archived file, so do not use a locally re-packed copy.

## What is generated

Default output directory: `conda-build/tests/output/` (git-ignored; `--output-dir` to change).

```
output/nisar_gcov_virtual_zarr_<granule_id>.zip
output/<granule_id>/
  nisar_gcov_<granule_id>.kerchunk.json   consolidated Kerchunk v1 / Zarr v2 references (.zmetadata included)
  gcov_chunk_table.csv / .md              one row per rank-2 dataset + subtotals + total
  gcov_dataset_inventory.csv              every dataset and link under /science/LSAR/GCOV
  granule_metadata.json                   CMR metadata, selection, stable URI, access probe
  environment.json                        Python/library/HDF5 versions, script hash, git commit
  validation_report.json / .md
  unsupported_lossy_report.json / .md
  README.md                               how to open the manifest
```

Absolute paths and byte sizes are printed at the end. Exit status is `0` when validation
passes and `3` when any validation failure was recorded (artifacts are still written).

### Manifest contents

* Every **rank-2** dataset under `/science/LSAR/GCOV` is a required 2-D subdataset regardless of
  GDAL `SUBDATASETS`; the chunk table counts them and reports why any is not in the manifest.
* Supporting arrays are also emitted and reported separately from the 2-D count: dimension-scale
  coordinates (`xCoordinates`, `yCoordinates`, ...), the scalar `projection` grid-mapping
  variables (always emitted, even if no array names them in `grid_mapping`), and numeric rank-3
  metadata cubes (`--no-rank3` leaves the cubes out).
* Chunk keys come from the HDF5 chunk index (`H5Dchunk_iter` via `h5py`, falling back to
  `H5Dget_chunk_info`); unallocated chunks get no key, so readers return `fill_value`, like HDF5.
  Contiguous datasets become one chunk; compact datasets are inlined.
* HDF5 filters map to numcodecs: deflate -> `zlib`, shuffle -> `shuffle`,
  Fletcher32 -> `fletcher32`. Any other filter (LZF, SZIP, Blosc, ...), or a chunk whose filter
  mask skips a filter, makes the dataset `skipped_unsupported`.
* Attributes are copied to `.zattrs` (bytes decoded to UTF-8, NaN/Inf as JSON strings,
  `DIMENSION_LIST`/`REFERENCE_LIST` replaced by `_ARRAY_DIMENSIONS`). A dimension whose scale
  lives in another group gets the qualified name `<group>__<scale>` (for example
  `metadata__coords__xCoordinates`) and a warning, because Zarr v2 cannot link to a coordinate
  in another group. This stops a same-named local coordinate from being attached instead. Scalar identification
  values are copied into the root `.zattrs` under `nisar_identification`.

### Chunk table columns

`hdf5_path, zarr_path, category, dimensions, dtype, chunk_shape, logical_chunks,
allocated_chunks, skipped_chunks, compressed_chunk_bytes, approx_json_bytes, filter_pipeline,
manifest_status`. Rows are sorted by HDF5 path, followed by subtotals for principal grids, masks,
coordinates, calibration grids and metadata groups, and a total row (datasets, allocated chunks,
skipped chunks, JSON bytes). `approx_json_bytes` is the serialized size of the dataset's
`.zarray`, `.zattrs` and chunk entries.

### Unsupported / lossy report

Lists, per dataset: non-deflate filters, virtual datasets, external and soft links (not
followed), extra hard links, compound/string/enum/reference/variable-length datatypes, missing
(unallocated) chunks, complex arrays (`<c8`/`<c16` are zarr-python/xarray readable but not
portable to every Zarr implementation), and attributes that could not be represented. Notable
global caveat: Zarr v2 has one `fill_value`; the manifest uses the HDF5 fill value while the
NISAR `_FillValue` attribute (often NaN) is kept in `.zattrs` - open with
`mask_and_scale=False`.

### Validation

* structure: valid Zarr keys, every key under a `.zgroup`, every reference targets the stable
  URI and lies inside the object size;
* per array: shape, dtype, chunks, compressor/filters and fill value against the source HDF5
  dataset, and chunk keys equal to the allocated HDF5 chunks (same offsets and sizes);
* reader side, without GDAL or HDF5: a raw range request decoded with numcodecs, `zarr.open_group`
  on the reference filesystem, sample chunks per array compared with h5py values, and
  `xarray.open_datatree` over the whole hierarchy.

Failures are recorded in `validation_report.*`, never silently dropped.

## Opening the manifest

```python
import earthaccess, xarray as xr
earthaccess.login()
fs = earthaccess.get_fsspec_https_session()
so = {"fo": "nisar_gcov_<granule_id>.kerchunk.json", "remote_protocol": "https",
      "remote_options": {"client_kwargs": fs.client_kwargs, "asynchronous": True},
      "asynchronous": True}
tree = xr.open_datatree("reference://", engine="zarr", consolidated=True, zarr_format=2,
                        storage_options=so, mask_and_scale=False, chunks=None)
tree["grids/frequencyA"]["HHHH"][:512, :512].values
```

`asynchronous=True` is needed with zarr-python 3, which drives fsspec asynchronously.

## Regenerating with a fixed seed

```bash
cd conda-build/tests
python nisar_gcov_virtual_zarr.py generate --seed 20251012 \
    --temporal 2025-10-01T00:00:00Z 2025-10-31T23:59:59Z --select-only   # check selection
python nisar_gcov_virtual_zarr.py generate --seed 20251012 \
    --temporal 2025-10-01T00:00:00Z 2025-10-31T23:59:59Z
pytest -v -p no:cacheprovider test_nisar_gcov_virtual_zarr.py            # offline tests
```

Never commit `output/`, `.netrc`, `.urs_cookies`, signed URLs or granule files.
