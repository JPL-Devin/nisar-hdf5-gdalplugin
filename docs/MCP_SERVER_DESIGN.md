# NISAR MCP server: design

Status: **draft for review**. Nothing described here is implemented yet.

## 1. Goal and scope

Expose NISAR HDF5 products to MCP clients (Claude Desktop, Devin, IDE agents) through a small
[Model Context Protocol](https://modelcontextprotocol.io) server that drives the existing
`gdal-driver-nisar` plugin. The server adds no product logic of its own: every answer comes from
GDAL opening a `NISAR:` connection string, so the driver stays the single source of truth for
product identification, subdatasets, georeferencing, masks and cube interpolation.

Non-goals:

- No changes to the driver. It stays a pure, read-only GDAL raster plugin with no Python or
  MCP dependency (section 7).
- No write/update paths into granules, no search of the ASF/CMR catalogue (use `earthaccess` or
  CMR directly and pass the resulting URLs in), no data hosting.

## 2. Architecture

```text
MCP client ──stdio / streamable HTTP──> nisar-mcp-server (Python, `mcp` SDK 2.x, MCPServer)
                                           │  osgeo.gdal (in-process)      ← default
                                           │  gdalinfo / gdal_translate    ← subprocess fallback
                                           ▼
                                  libgdal + gdal_NISAR.so (conda env with gdal-driver-nisar)
                                           ▼
                           local .h5 │ /vsicurl/https://… │ /vsis3/bucket/…
```

- **Language / SDK:** Python ≥ 3.11 with the official `mcp` SDK **2.x** (`>=2.3,<3`), using
  `mcp.server.mcpserver.MCPServer`, which is what FastMCP is called in 2.x. Tools are plain typed
  functions; the SDK derives the JSON schemas from type hints. 2.x was chosen over 1.x to match
  `earthdata-mcp-server` and `reactor_mcp_server` (section 8), and because no 1.x code exists that
  would have to be migrated.
- **Execution mode:** the default is the in-process `osgeo.gdal` API (`gdal.Info`,
  `gdal.Translate`, `gdal.Open` + `ReadAsArray`). It avoids a process start per call and keeps
  GDAL's block cache and `/vsicurl/` metadata cache warm across calls on the same granule.
  `NISAR_MCP_EXEC=subprocess` switches to running `gdalinfo` / `gdal_translate` from the same
  conda env's `bin/` as a child process. That isolates a driver crash (a segfault in GDAL/HDF5
  kills only the child) and makes it easy to enforce a hard timeout, at the cost of re-opening the
  granule each call. Both modes must produce identical tool output; the tables in section 4 give
  both forms.
- **Driver selection:** always open with the `NISAR:` prefix (or `allowed_drivers=["NISAR"]`).
  Without it a remote file is only claimed when the URL contains `NISAR`, and driver order can
  hand local files to the generic `HDF5` driver (see the `nisar-gdal` skill).
- **Concurrency:** GDAL dataset handles are not thread-safe. Each tool call opens its own
  dataset handle (opens are cheap after the first, because the VSI cache is shared); blocking GDAL
  work runs in a worker thread (`anyio.to_thread.run_sync`) so the MCP event loop stays responsive.
  A small bounded semaphore (`NISAR_MCP_MAX_CONCURRENCY`, default 4) caps parallel remote reads.
- **Startup check:** on launch the server runs the equivalent of `gdalinfo --format NISAR` and
  refuses to start if the driver is missing. It records `DRIVER_VERSION` and reports it in the
  MCP `serverInfo` and in every `get_metadata` result.

## 3. Common conventions

### 3.1 Product references

Every tool and resource takes a `path` that names a granule: a local path, `s3://bucket/key.h5`,
`/vsis3/…`, `https://…` or `/vsicurl/https://…`, exactly as the driver accepts it. Tools that act
on one raster also take an optional `subdataset` (an HDF5 path such as
`/science/LSAR/GCOV/grids/frequencyA/HHHH`), or accept a full `NISAR:"…":/…` string copied
from `list_subdatasets`. The server builds connection strings as

```text
NISAR:"<path>"                    container
NISAR:"<path>":<hdf5-dataset>     one raster
```

and never asks the model to type an HDF5 path from memory: tool descriptions say to call
`list_subdatasets` first and to copy paths from it. That is the same rule as in the `nisar-gdal`
skill.

Optional `open_options` (a dict) map 1:1 onto `-oo KEY=VALUE`. They are validated against the
driver's `DMD_OPENOPTIONLIST` (read once at startup), so currently `INST`, `FREQ`, `POL`,
`METADATA`, `MASK`, `QUANTITY`, `DEM_FILE`, `DEM_RESAMPLING`, `DEM_NODATA_HEIGHT`, `DUMP`,
`DUMP_ROOT`, `DUMP_MODE` and `ENABLE_PAGE_BUFFERING`. Unknown keys are rejected rather than
passed through.

### 3.2 Path policy

- Local paths must resolve (after `realpath`) under one of the roots in `NISAR_MCP_ALLOWED_ROOTS`
  (colon-separated; defaults to `NISAR_TEST_DATA_DIR` if set, else the working directory).
- Remote prefixes must match `NISAR_MCP_ALLOWED_PREFIXES` (default:
  `https://nisar.asf.earthdatacloud.nasa.gov/`, `s3://sds-n-cumulus-prod-nisar-products/`,
  `s3://nisar-ops-rs-fwd/`).
- `/vsi*` chains other than `/vsicurl/` and `/vsis3/` (for example `/vsizip/`,
  `/vsistdin/` or `/vsimem/`) are rejected in input paths.
- Output files (`get_subset`, `translate`) are written only under `NISAR_MCP_OUTPUT_DIR`
  (default `$XDG_CACHE_HOME/nisar-mcp/out`), with server-generated names. The caller gets back a
  `file://` URI plus a `nisar://output/{id}` resource (section 5).

### 3.3 Errors and limits

- GDAL errors are caught (`gdal.UseExceptions()`, or the subprocess's stderr) and returned as MCP
  tool errors (`isError: true`) with the driver's message verbatim. Driver messages such as
  "DEM_FILE open option is REQUIRED when QUANTITY is specified." are already actionable.
- Every call gets a wall-clock timeout (`NISAR_MCP_TIMEOUT_S`, default 300 s). In subprocess
  mode the child is killed when it runs out. In-process, the server stops waiting and keeps the
  worker slot busy until GDAL returns, which is one reason subprocess mode exists.
- Size guards: `get_subset` with `output="array"` or `"stats"` refuses windows over
  `NISAR_MCP_MAX_PIXELS` (default 16 M pixels per band) unless the call asks for a decimated
  overview (`overview_level`). `translate` has no pixel cap but streams to disk.

## 4. Tools

All tools are read-only and idempotent (MCP annotations `readOnlyHint: true`,
`idempotentHint: true`; `openWorldHint: true` for remote paths). `get_subset` and `translate`
write only into the server's output dir.

### 4.1 `list_products`

| | |
|---|---|
| Arguments | `prefix: str` (local dir, `/vsicurl/https://…/`, `/vsis3/bucket/prefix/`, or `s3://…`), `pattern: str = "product_*.h5"`, `recursive: bool = False`, `limit: int = 500` |
| Returns | `[{path, size_bytes, product_type?, level?}]`. The product type and level are parsed from the `NISAR_<L>_PR_<TYPE>_…` file-name convention when it matches, else omitted. The tool does not open the files. |
| GDAL call | `gdal.ReadDirRecursive(prefix)` / `gdal.ReadDir(prefix)` + `fnmatch`; `gdal.VSIStatL` for sizes. Subprocess equivalent: `gdal vsi list [-R] <prefix>` (GDAL ≥ 3.11 CLI). |

Notes:

- `ReadDir` works on `/vsis3/` (ListObjectsV2) and on local directories. On `/vsicurl/` it
  works only where the server returns an HTML index or an S3-style listing. ASF's CloudFront
  `https://` endpoint does not provide one, so for HTTPS the tool documents that the caller
  should pass explicit URLs, for example from `earthaccess`.
- The default `pattern` is `product_*.h5` as specified. Operational granules are named
  `NISAR_L?_PR_<TYPE>_….h5`, so clients will usually pass `pattern="NISAR_*.h5"`. **Open
  question:** should the default become `*.h5`?
- Set `GDAL_DISABLE_READDIR_ON_OPEN=EMPTY_DIR` for the *open* tools (below) but not for this
  one.

### 4.2 `get_metadata`

| | |
|---|---|
| Arguments | `path: str`, `subdataset: str \| None`, `open_options: dict = {}`, `domains: list[str] = []` (for example `["NISAR_ORBIT"]`, `["NISAR_DUMP"]`), `include_gcps: bool = False` |
| Returns | gdalinfo JSON, trimmed: `driverShortName`, `driverLongName`, `nisar_driver_version`, `size` `[x, y]`, `bands` (type, block size, nodata, description, unit), `coordinateSystem.wkt` + `epsg` (from `GetSpatialRef().GetAuthorityCode`), `geoTransform` (L2/L3) **or** `gcp_count` + `gcp_crs` (L1, which has no GeoTransform), `cornerCoordinates`/`wgs84Extent` when they exist, `metadata` for the requested domains, `subdatasets` (container only, as in 4.3). |
| GDAL call | `gdal.Info(ds, format="json", showGCPs=include_gcps, extraMDDomains=domains, computeMinMax=False, stats=False)` on `gdal.OpenEx('NISAR:"<path>"[:<sds>]', open_options=[…])`. Subprocess: `gdalinfo -json [-nogcp] [-mdd <d>]… [-oo K=V]… 'NISAR:"<path>"[:<sds>]'`. |

Notes:

- `showGCPs=False` / `-nogcp` is the default because L1 swaths carry thousands of GCPs.
- The tool never requests `-stats`/`-approx_stats`. The driver answers those from
  `min_value`/`max_value` attributes and synthesises mean and stddev (see the skill's
  "Statistics" section), which a model could mistake for measured values. Measured statistics
  come only from `get_subset(output="stats")`.
- `domains=["NISAR_DUMP"]` works without `DUMP=YES`, which is the driver's own behaviour.
  `METADATA=ORBIT` etc. must be passed in `open_options` for the corresponding `NISAR_<GROUP>`
  domain to be populated; the tool adds it automatically when a `NISAR_<GROUP>` domain is
  requested.

### 4.3 `list_subdatasets`

| | |
|---|---|
| Arguments | `path: str`, `filter: str \| None` (substring or glob on the HDF5 path, for example `"frequencyA"`, `"*HHHH"`), `open_options: dict = {}` (for example `INST=SSAR`) |
| Returns | `[{name, hdf5_path, dims, dtype, description}]`, where `name` is the verbatim `SUBDATASET_n_NAME` (`NISAR:"…":/science/…`, ready to pass back as `subdataset`), and `dims`/`dtype` are parsed from `SUBDATASET_n_DESC` (`[21x720x748] …/incidenceAngle (Float32)`). |
| GDAL call | `gdal.OpenEx('NISAR:"<path>"').GetMetadata("SUBDATASETS")`. Subprocess: `gdalinfo -json 'NISAR:"<path>"'` → `.metadata.SUBDATASETS`. |

Notes:

- The list contains only entries the driver can open. String scalars such as
  `boundingPolygon` are not subdatasets; `get_metadata` returns them in the container's default
  domain.
- Derived views (`DERIVED_SUBDATASET:AMPLITUDE:NISAR:…`) are returned with
  `include_derived=True`, from the raster's `DERIVED_SUBDATASETS` domain. They cost one extra
  open per raster, so they are off by default.

### 4.4 `get_subset`

| | |
|---|---|
| Arguments | `path`, `subdataset`, `open_options`; exactly one of `window: {xoff, yoff, xsize, ysize}` (pixels) or `bbox: {minx, miny, maxx, maxy, crs="EPSG:4326"}`; `bands: list[int] \| None`; `overview_level: int \| None`; `output: "cog" \| "stats" \| "array" = "stats"`; `derived: "AMPLITUDE" \| "PHASE" \| "INTENSITY" \| "LOGAMPLITUDE" \| … \| None` |
| Returns | `stats`: per band `{min, max, mean, std, count_valid, count_nodata}`, computed from pixels (NaN- and mask-aware with `MASK=YES`). `cog`: `{uri, resource, size_bytes, width, height, crs, geotransform}`. `array`: a small inline array (capped at 512×512) as nested lists. |
| GDAL call | **window:** `gdal.Translate(dst, ds, srcWin=[xoff, yoff, xsize, ysize], bandList=bands, format="COG", creationOptions=["COMPRESS=ZSTD", "PREDICTOR=YES", "BIGTIFF=IF_SAFER"])` = `gdal_translate -of COG -srcwin … [-b N]… -co COMPRESS=ZSTD …`. **bbox on L2/L3:** `projWin=[minx, maxy, maxx, miny], projWinSRS=crs` = `gdal_translate -projwin … -projwin_srs <crs>`. **bbox on L1** (GCPs only, no GeoTransform): the bbox is converted to a pixel window with `gdal.Transformer(ds, None, ["METHOD=GCP_POLYNOMIAL", "DST_SRS=<crs>"])` (inverse transform of the bbox corners plus densified edges), then `srcWin` as above. The output keeps the GCPs; the result says it is radar geometry. **stats:** the same `srcWin` read as `band.ReadAsArray(xoff, yoff, xsize, ysize)` (per block for large windows) and reduced in NumPy, or the subprocess equivalent `gdal_translate -of MEM` piped to `gdalinfo -stats` on a temporary VRT. |

Notes:

- `derived` wraps the connection string as `DERIVED_SUBDATASET:<D>:NISAR:"<path>":<sds>`.
  The mask is not inherited by derived views, which is driver behaviour and is documented in
  the tool description.
- `overview_level` reads from the driver's virtual power-of-two overviews
  (`band.GetOverview(i)`, or `gdal_translate -ovr i`). It is the cheap way to get a quick look of a
  whole GCOV or GSLC.
- For large remote windows the server sets `NISAR_PREFETCH_GRID` (for example `24`) in
  `gdal.config_options(...)` around the call. For small windows it leaves the default `1`.
- `MASK` defaults to `YES` for `stats` and `array` (so fill values are excluded) and to the
  driver default (`NO`) for `cog`, which copies the mask as an internal TIFF mask when the
  caller sets it. The value used is echoed in the result.

### 4.5 `sample_at_points`

| | |
|---|---|
| Arguments | `path`, `subdataset`, `open_options`, `points: list[{lat, lon, height?}]` (≤ 10 000), `resampling: "nearest" \| "bilinear" = "nearest"`, `quantity: str \| None`, `dem_file: str \| None`, `dem_resampling: str \| None` |
| Returns | `[{lat, lon, pixel, line, value[s], inside: bool}]`, plus `{product_type, grid_type: "MAP" \| "RADAR", subdataset, quantity?}` |
| GDAL call | see the three paths below |

Lat/lon is mapped to a grid value differently by product level:

1. **L2/L3 grids (GSLC, GCOV, GUNW, GOFF, SME2) and radar-grid cubes:** the dataset has a
   GeoTransform and SRS, so `osr.CoordinateTransformation(EPSG:4326 → ds SRS)` (with
   `OAMS_TRADITIONAL_GIS_ORDER`), then the inverse GeoTransform, then
   `band.ReadAsArray(px, ln, 1, 1)`. For `bilinear` the server reads a 2×2 neighbourhood and
   interpolates. Subprocess: `gdallocationinfo -valonly -wgs84 [-r bilinear] [-b N] -oo …
   'NISAR:"<path>":<sds>'` with points on stdin.
2. **L1 swaths (RSLC, RIFG, RUNW):** the dataset has GCPs and no GeoTransform, and
   `gdallocationinfo -wgs84` fails on it (verified with GDAL 3.12.4: exit 1, no output).
   The server maps (lon, lat) to (pixel, line) with the inverse of
   `gdal.Transformer(ds, None, ["METHOD=GCP_TPS", "DST_SRS=EPSG:4326"])`. `GCP_TPS` is the
   README's recommendation for accuracy, and `GCP_POLYNOMIAL` can be selected through
   `l1_method`. The value is then read as in (1). The result is flagged as a GCP-fit
   approximation. It is a quick look, not the full geometry, as the skill's warping section
   says for L1.
3. **Cube quantities via the driver's interpolation mode** (`quantity` set, for example
   `incidenceAngle`, `losUnitVectorX`, `slantRange`): the server opens
   `NISAR:"<path>"[:<cube>]` with `open_options` `QUANTITY=<q>`, `DEM_FILE=<dem>`,
   `DEM_RESAMPLING=…`, `DEM_NODATA_HEIGHT=…` and `FREQ`/`POL`/`INST`. The driver then
   returns a `NisarInterpolatedDataset` on the product's target grid (L2/L3: map grid,
   interpolated in height at DEM terrain; L1: radar grid with the terrain height solved per
   pixel and GCPs passed through). Points are located on that grid with (1) or (2) and read
   one block at a time. All height interpolation is done by the driver; the server does no
   cube maths. The default `dem_file` comes from `NISAR_MCP_DEFAULT_DEM` (for example
   `https://nisar.asf.earthdatacloud.nasa.gov/NISAR/DEM/v1.2/EPSG4326/EPSG4326.vrt`). If none is
   configured, the server returns the driver's own "DEM_FILE open option is REQUIRED" error.

Points outside the raster get `inside: false` and `value: null`. They are never clamped to
the edge.

### 4.6 `translate`

| | |
|---|---|
| Arguments | `path`, `subdataset`, `open_options`, `format: "COG" \| "netCDF" = "COG"`, `bands`, `window`/`bbox` (as in 4.4), `derived`, `creation_options: list[str] = []`, `warp_to: str \| None` (target CRS, for example `"EPSG:4326"`), `resampling: str = "near"` |
| Returns | `{uri, resource, format, size_bytes, width, height, crs, driver_version, command}`, where `command` is the equivalent CLI line for reproducibility |
| GDAL call | **COG:** `gdal.Translate(dst, src, format="COG", creationOptions=["COMPRESS=ZSTD", "PREDICTOR=YES", "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS", *creation_options])` = `gdal_translate -of COG -co …`. **netCDF:** `gdal.Translate(dst, src, format="netCDF", creationOptions=["FORMAT=NC4", "COMPRESS=DEFLATE", *creation_options])` = `gdal_translate -of netCDF -co FORMAT=NC4 -co COMPRESS=DEFLATE …`. **With `warp_to`:** `gdal.Warp(dst, src, dstSRS=warp_to, resampleAlg=resampling, format=…, multithread=True, srcNodata="nan", dstNodata="nan" (float bands), tps=True (L1))` = `gdalwarp -t_srs … -r … [-tps] -multi -wo NUM_THREADS=ALL_CPUS …`. |

Notes:

- The netCDF output needs the `libgdal-netcdf` plugin in the same environment. That is a
  dependency of the server package, not of the driver (section 7).
- The tool description carries the resampling guidance from the skill: nearest for phase,
  complex, mask and connected-component layers; bilinear for backscatter and coherence.
- Long translations report MCP progress notifications from `gdal.TermProgress`-style callbacks
  (`callback=` on `gdal.Translate`/`gdal.Warp`).

### 4.7 Summary: tool → GDAL

| Tool | Python API (default) | Subprocess fallback |
|---|---|---|
| `list_products` | `gdal.ReadDir[Recursive]`, `gdal.VSIStatL` | `gdal vsi list [-R]` |
| `get_metadata` | `gdal.Info(ds, format="json", showGCPs=False, extraMDDomains=…)` | `gdalinfo -json -nogcp [-mdd …] [-oo …]` |
| `list_subdatasets` | `ds.GetMetadata("SUBDATASETS")` | `gdalinfo -json` → `metadata.SUBDATASETS` |
| `get_subset` (cog) | `gdal.Translate(format="COG", srcWin/projWin)` | `gdal_translate -of COG -srcwin/-projwin` |
| `get_subset` (stats) | `band.ReadAsArray(win)` + NumPy | `gdal_translate -of VRT -srcwin` + `gdalinfo -stats -json` |
| `sample_at_points` L2/L3 | `osr` transform + inverse GeoTransform + `ReadAsArray` | `gdallocationinfo -valonly -wgs84` |
| `sample_at_points` L1 | `gdal.Transformer(METHOD=GCP_TPS)` inverse + `ReadAsArray` | (none; `gdallocationinfo` needs a GeoTransform) |
| `sample_at_points` cube | as above on `OpenEx(..., open_options=[QUANTITY, DEM_FILE, …])` | `gdallocationinfo -oo QUANTITY=… -oo DEM_FILE=…` (L2/L3 only) |
| `translate` | `gdal.Translate` / `gdal.Warp` | `gdal_translate` / `gdalwarp` |

## 5. Resources

Resources expose per-product metadata that a client can attach to context without a tool call.
`{path}` is the product reference from 3.1, percent-encoded (RFC 6570 `{+path}` is not used, so
`/` and `:` inside the path are encoded).

| URI template | MIME | Content | GDAL call |
|---|---|---|---|
| `nisar://subdatasets/{path}` | `application/json` | Same payload as `list_subdatasets(path)` | `GetMetadata("SUBDATASETS")` |
| `nisar://metadata/{path}` | `application/json` | Container `get_metadata(path)`: product type, level, instrument, frequencies/polarizations, `boundingPolygon`, time span, driver version | `gdal.Info(format="json")` on `NISAR:"<path>"` |
| `nisar://metadata/{path}?subdataset={sds}` | `application/json` | Raster-level `get_metadata` | as above on `NISAR:"<path>":<sds>` |
| `nisar://dump/{path}` *(optional)* | `text/plain` | `NISAR_DUMP` HEADER walk of `/science/<INST>` | `ds.GetMetadata_List("NISAR_DUMP")` |
| `nisar://output/{id}` | `image/tiff` / `application/x-netcdf` | Files produced by `get_subset`/`translate` | — |

- Resources are listed through `resources/templates/list`, not enumerated, because the set of
  products is unbounded. After `list_products`, the server may send
  `notifications/resources/list_changed` with concrete URIs for the products it found.
- Metadata results are cached in memory by `(path, subdataset, open_options, driver_version)`
  with a TTL (`NISAR_MCP_CACHE_TTL_S`, default 3600 s). Granules are immutable once published,
  so the TTL mainly bounds memory.

## 6. Credentials

**Credentials never appear in tool arguments, resource URIs, tool results or logs.** They are
supplied only through the MCP server process environment, which the MCP client sets in its
server config. GDAL then picks them up exactly as it does for the CLI. The driver has no
credential code of its own: all remote I/O goes through GDAL VSI.

| Need | Environment variables / files (server env only) |
|---|---|
| Local test data | `NISAR_TEST_DATA_DIR` (also the default `NISAR_MCP_ALLOWED_ROOTS`) |
| Direct S3 (in-region, `us-west-2`) | `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_SESSION_TOKEN`, or `AWS_PROFILE` (+ `~/.aws/`), `AWS_REGION=us-west-2`; `AWS_NO_SIGN_REQUEST=YES` for public buckets; instance/role credentials work unchanged |
| HTTPS from ASF / Earthdata | `~/.netrc` with `machine urs.earthdata.nasa.gov login … password …` (path overridable with `NETRC`), `GDAL_HTTP_NETRC=YES` (default), `GDAL_HTTP_COOKIEFILE`/`GDAL_HTTP_COOKIEJAR` pointing to a server-private cookie file; or a bearer token through `EARTHDATA_TOKEN` → the server sets `GDAL_HTTP_HEADERS="Authorization: Bearer $EARTHDATA_TOKEN"`; or `EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` (the variables `earthdata-mcp-server` reads) → at startup the server writes a server-private netrc (mode 0600, in its cache dir) and points `GDAL_HTTP_NETRC_FILE` at it. GDAL does not read the `EARTHDATA_*` variables itself. |
| Tuning (not secret) | `GDAL_HTTP_MAX_RETRY`, `GDAL_CACHEMAX`, `GDAL_NUM_THREADS`, `NISAR_PREFETCH_GRID`, `NISAR_MAX_MEGAFETCH_BYTES` |

Example client config (Claude Desktop / Devin MCP JSON):

```json
{
  "mcpServers": {
    "nisar": {
      "command": "/opt/conda/envs/nisar/bin/nisar-mcp-server",
      "env": {
        "NISAR_TEST_DATA_DIR": "/data/nisar",
        "AWS_PROFILE": "saml-pub",
        "AWS_REGION": "us-west-2",
        "NETRC": "/home/me/.netrc",
        "NISAR_MCP_DEFAULT_DEM": "https://nisar.asf.earthdatacloud.nasa.gov/NISAR/DEM/v1.2/EPSG4326/EPSG4326.vrt"
      }
    }
  }
}
```

Rules the server enforces:

- Tool argument schemas have no credential fields. An `open_options` key or `path` containing
  `AWS_`, `GDAL_HTTP_`, `token`, `password` or a userinfo part (`https://user:pass@…`) is
  rejected.
- Per-call GDAL config is limited to an allow-list of non-secret tuning keys. Callers cannot
  call `gdal.SetConfigOption` / `--config` for `AWS_*`, `GDAL_HTTP_HEADERS`,
  `GDAL_HTTP_USERPWD` and similar.
- GDAL error text is scrubbed for signed-URL query strings (`X-Amz-Signature=…`,
  `X-Amz-Credential=…`, CloudFront `Signature=`/`Policy=`) before it is returned. `CPL_DEBUG`
  output is never forwarded to the client.
- In subprocess mode the child inherits the server env. Credentials are never put on the
  command line, so they cannot leak through `ps` or the returned `command` string.
- Streamable-HTTP transport binds to `127.0.0.1` by default. Remote exposure needs MCP auth
  (OAuth per the MCP authorization spec) in front of it, and the server's credentials are those
  of the deployment, not of the caller.

## 7. Packaging

The server ships as a **separate conda package, `nisar-mcp-server`**, on the same `nisar-forge`
channel. The driver package `gdal-driver-nisar` is unchanged: it stays one shared library plus
docs under `lib/gdalplugins/` and `share/doc/`, with no Python and no MCP dependency. That
keeps it a pure GDAL plugin that can be upstreamed into OSGeo/gdal (`UPSTREAM_CHECKLIST.md`).

```yaml
# conda-build/nisar-mcp-recipe/meta.yaml (sketch)
package:
  name: nisar-mcp-server
  version: 0.1.0
source:
  path: ../../mcp            # pure-Python package: mcp/nisar_mcp/…, mcp/pyproject.toml
build:
  noarch: python
  script: {{ PYTHON }} -m pip install . --no-deps -vv
  entry_points:
    - nisar-mcp-server = nisar_mcp.server:main
requirements:
  host: [python >=3.11, pip, hatchling]
  run:
    - python >=3.11
    - gdal-driver-nisar >=0.7.2      # brings the matching gdal / libgdal-core pins
    - gdal                           # osgeo.gdal Python bindings + CLI tools
    - libgdal-netcdf                 # translate(format="netCDF")
    - numpy
    - mcp >=2.3,<3                   # official Python SDK 2.x (conda-forge; brings httpx2, mcp-types, opentelemetry-api)
test:
  imports: [nisar_mcp]
  commands:
    - gdalinfo --formats | grep NISAR
    - nisar-mcp-server --self-test   # loads driver, lists tools, exits 0
```

- `noarch: python` (one build for all platforms). The platform-specific binary stays in
  `gdal-driver-nisar`, which already pins `gdal`/`libgdal-core`/`hdf5` to its build versions,
  so the server cannot pull in a GDAL the plugin was not built against.
- The source lives in this repository under `mcp/` (Python package `nisar_mcp`), outside
  `frmts/nisar/`, so the driver's CMake build, the in-tree GDAL build and the upstreaming
  diff never see it. Versioning is independent of the driver. The server declares the minimum
  `DRIVER_VERSION` it needs and checks it at startup.
- Install:

  ```bash
  conda create -n nisar-mcp -c nisar-forge -c conda-forge nisar-mcp-server
  conda run -n nisar-mcp nisar-mcp-server            # stdio
  conda run -n nisar-mcp nisar-mcp-server --http 127.0.0.1:8765
  ```

- Also publishable to PyPI as `nisar-mcp-server`, but PyPI cannot provide the GDAL plugin. The
  conda package is the supported route, and the PyPI install documents that it needs an env
  with `gdal-driver-nisar`.
- The `nisar-gdal` skill (`SKILL.md`) is reused as the server's MCP `instructions` string
  (trimmed) and exposed as a prompt `nisar://docs/skill`. The tool descriptions therefore carry
  the same guidance on inventory-first use, `MASK` and statistics as the CLI docs.

## 8. Composing with earthdata-mcp-server

[`datalayer/earthdata-mcp-server`](https://github.com/datalayer/earthdata-mcp-server) does
discovery and download, and this server reads inside granules, so the two complement each other.
Its tools are `search_earth_datasets`, `search_earth_datagranules` and
`download_earth_data_granules` (`mode="manifest" | "script" | "download"`). It does CMR search
through `earthaccess`, and it never opens a file.

Workflow:

1. `search_earth_datagranules(short_name=…, temporal=…, bounding_box=…)` or
   `download_earth_data_granules(..., mode="manifest")` on the earthdata server. Each manifest item
   carries `links` (earthaccess `data_links`, HTTPS by default).
2. Pass the `.h5` link as `path` to `list_subdatasets` / `get_metadata` / `get_subset` /
   `sample_at_points` / `translate` here. The data is read by range over `/vsicurl/`, with no
   download. This replaces `list_products` for HTTPS, where ASF offers no directory listing.

Integration points:

- **Same SDK major.** Both servers use `mcp` 2.x. They can run as two separate MCP servers in the
  client config, or as two toolsets on one `reactor_mcp_server` host (`/mcp?earthdata&nisar`). For
  the latter this server would also ship a `reactor.mcp.extensions` entry point. That is optional,
  because the stand-alone `nisar-mcp-server` entry point stays the primary interface.
- **Shared credentials.** One env block with `EARTHDATA_USERNAME`/`EARTHDATA_PASSWORD` serves both
  servers (section 6).
- **Link selection.** The manifest returns at most 5 links per granule, and NISAR granules can carry
  non-`.h5` links. The tool descriptions here say to pass the link that ends in `.h5`, and `path`
  validation rejects anything else, which matches the driver's own `.h5` requirement. Still to be
  checked on real NISAR manifests: whether the `.h5` link is always within the first 5.
- **No overlap.** Its `mode="download"` (whole granules to local disk) is needed only for repeated
  offline work. The local path it returns works here as well, provided it is under
  `NISAR_MCP_ALLOWED_ROOTS`.

## 9. Testing

- Unit tests (no network): a tiny synthetic NISAR-layout HDF5 built with `h5py` in a fixture
  (as `conda-build/tests/test_dump.py` does) covers argument validation, path policy, credential
  rejection, JSON shapes, and parity between the in-process and subprocess modes.
- Integration: reuse the granule discovery of `conda-build/tests/test_interpolation_earthaccess.py`
  (`NISAR_TEST_DATA_DIR` for local copies, else earthaccess + `~/.netrc` over HTTPS;
  `NISAR_TEST_ACCESS=direct` + `AWS_REGION=us-west-2` in-region). Each tool is called through the
  SDK's in-memory client session against GCOV, GSLC, GUNW and RSLC, and results are compared
  with direct `osgeo.gdal` calls on the same connection strings.
- `sample_at_points` cube mode is checked against `gdal_translate -oo QUANTITY=…` output at the
  same pixels, so it inherits the driver's interpolation tests rather than duplicating them.

## 10. Open questions

1. Default `list_products` pattern: keep `product_*.h5` or change to `NISAR_*.h5` / `*.h5`?
2. Default execution mode: in-process (faster) or subprocess (crash-isolated)?
3. Should the server live in this repository under `mcp/`, or in its own repository so this
   repository holds only the driver?
4. Should `get_subset(output="cog")` also return a small PNG quick-look as MCP image content
   for clients that render it?
