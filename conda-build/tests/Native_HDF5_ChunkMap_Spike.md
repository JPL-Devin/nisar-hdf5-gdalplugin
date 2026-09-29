# Native HDF5 chunk-map parser: feasibility spike

Scope: replace only the chunk-index discovery in `NisarRasterBand::MapChunks()` (`H5Dchunk_iter`) with a
minimal read-only parser (`conda-build/nisar-gdal-recipe/hdf5native.{h,cpp}`) that reads through GDAL VSI.
libhdf5 is unchanged and remains the default; `NISAR_NATIVE_HDF5=YES` enables the native path, and any
non-OK `NisarHDF5Native::Status` falls back to `H5Dchunk_iter`.

Environment: GDAL 3.12.4, HDF5 2.2.0 (conda-forge), x86_64 Linux, outside us-west-2 (HTTPS via Earthdata).

## Reproduce

```bash
cd conda-build/nisar-gdal-recipe
cmake -S . -B build -DNISAR_BUILD_CHUNKMAP_BENCH=ON && cmake --build build -j
# validation (exit code 1 if any dataset diverges)
build/nisar_chunkmap_bench NISAR_*.h5
# synthetic files covering format features absent from current granules
build/nisar_chunkmap_synth /tmp/synth && build/nisar_chunkmap_bench /tmp/synth/*.h5
# benchmark (local or /vsicurl/, /vsis3/)
build/nisar_chunkmap_bench --no-validate --repeat 3 --block-sizes 0,65536,4194304 \
    --bench /science/LSAR/GCOV/grids/frequencyA/HHHH /vsicurl/https://.../NISAR_L2_PR_GCOV_....h5
# in-driver
NISAR_NATIVE_HDF5=YES CPL_DEBUG=NISAR_CHUNKMAP gdalinfo -checksum NISAR:granule.h5:/science/...
```

Validation compares, per chunked 2-D/3-D dataset: dims, chunk dims, filter IDs, and every
`(offset[], filter_mask, addr, size)` record from `H5Dchunk_iter`, then the derived per-band
`(blockY, blockX, offset, length, isMissing)` maps built exactly as `MapChunks()` does (3-D: one map per band).

## Validation summary (real granules)

| product | file | superblock | HDF5 objects | chunked 2-D/3-D | identical | divergent | native fallback | index types |
|---|---|---:|---:|---:|---:|---:|---:|---|
| GCOV | NISAR_L2_PR_GCOV_004_159_A_024_2005_DHDH_A_20251109T051646_20251109T051650_X05010_N_P_J_001 | v2 | 285 | 38 | 38 | 0 | 0 | btree_v1 |
| RSLC | NISAR_L1_PR_RSLC_004_159_A_024_2005_DHDH_A_20251109T051646_20251109T051650_X05010_N_P_J_001 | v2 | 180 | 13 | 13 | 0 | 0 | btree_v1 |
| GUNW | NISAR_L2_PR_GUNW_003_005_A_018_004_4000_SH_20251017T125111_20251017T125127_20251029T125111_20251029T125131_X05010_N_P_J_001 | v2 | 217 | 29 | 29 | 0 | 0 | btree_v1 |
| GSLC | NISAR_L2_PR_GSLC_004_159_A_024_2005_DHDH_A_20251109T051646_20251109T051650_X05010_N_P_J_001 | v2 | 270 | 33 | 33 | 0 | 0 | btree_v1 |

All 113 chunked 2-D/3-D datasets are byte-identical. All current granules use superblock v2, object header v1,
symbol-table groups, layout message v3 and v1 B-tree chunk indexes; fixed/extensible array indexes were not
encountered. In-driver `gdalinfo -checksum` is identical with `NISAR_NATIVE_HDF5=YES|NO` for GCOV HHHH, GCOV
frequencyB HVHV, the 21-band GCOV `incidenceAngle` cube, RSLC HH, GSLC HH and GUNW unwrappedPhase.

## Synthetic coverage (`nisar_chunkmap_synth`)

Each file has 18 chunked datasets (filtered/unfiltered, partially written so some chunks are missing) at the
root, inside a 20 000-link group (`/big`, dense for v1.8+: fractal heap with nested indirect blocks + v2 B-tree
name index) and inside a subgroup in the middle of that group.

| file | superblock / OH | identical | divergent | native fallback | index types |
|---|---|---:|---:|---:|---|
| synth_sb0_earliest.h5 | v0 / v1, symbol tables | 18 | 0 | 0 | btree_v1 |
| synth_sb2_v18.h5 | v2 / v2, dense groups | 18 | 0 | 0 | btree_v1 |
| synth_sb3_v110.h5 | v3 / v2, dense groups, layout v4 | 17 | 1 | 0 | single_chunk, implicit, fixed_array (paged), extensible_array, btree_v2 |
| synth_sb3_latest.h5 | v3 / v2, layout v4/v5 | 6 | 0 | 12 | btree_v2, fixed_array, implicit, single_chunk, extensible_array |

Divergences / fallbacks:

- `synth_sb3_v110.h5:/earray_2d_inner_unlim` (dims 40x600, maxdims 40xUNLIMITED, chunk 8x8, extensible array
  with the unlimited dimension not first): **libhdf5 2.2.0 is wrong, native is correct.** `H5Dchunk_iter`,
  `H5Dget_chunk_info` and `H5Dget_chunk_info_by_coord` report swizzled/linearised offsets (e.g. `(0,8)`,
  `(0,24)`, ..., up to `(0,2992)` for a 600-column dataset), while `H5Dread` returns the correct values. Decoding
  the raw chunk bytes (inflate + unshuffle) confirms the native mapping (e.g. address 26209251 holds chunk
  `(8,0)`, first value 94 = expected for that chunk; libhdf5 reports it as `(0,8)`). The existing libhdf5
  `MapChunks()` path would therefore mis-map such datasets; none exist in current NISAR products.
- `synth_sb3_latest.h5`: 12 filtered datasets use data layout message **version 5** (written by HDF5 2.x with
  `libver=latest`). The spike does not implement v5; it returns `ERR_LAYOUT` and falls back cleanly.

## Benchmark: chunk-map discovery, GCOV `/science/LSAR/GCOV/grids/frequencyA/HHHH` (1225 chunks)

Rows:
- *libhdf5 cold*: `H5Fopen` (plugin VFL, 4 MiB page buffer as in the driver) + `H5Dopen` + `H5Dchunk_iter`.
- *libhdf5 iter only*: `H5Dchunk_iter` on an already-open dataset (what `MapChunks()` costs today).
- *native cold*: fresh VSI handle, superblock -> groups -> object header -> chunk index. `block=N` reads metadata
  in aligned N-byte blocks cached per `File`; `block=0` issues exact-size reads.
- *native in-driver*: native run after libhdf5 has opened the file (VSI/curl cache warm), i.e. what
  `MapChunks()` costs with `NISAR_NATIVE_HDF5=YES`.

"VSI reads/bytes" are `VSIFReadL` calls/bytes (libhdf5: counted in the VFL; native: counted in the parser).
"HTTP" columns come from `VSINetworkStatsGetAsSerializedJSON` after `VSICurlClearCache()`.

### Remote, HTTPS /vsicurl, `CPL_VSIL_CURL_CHUNK_SIZE=8388608` (repository test settings), median of 3

| path | wall ms | VSI reads | VSI bytes | HTTP GET | HTTP HEAD | HTTP bytes |
|---|---:|---:|---:|---:|---:|---:|
| libhdf5 cold: H5Fopen+H5Dopen+H5Dchunk_iter | 1330.5 | 8 | 20972088 | 3 | 1 | 16777216 |
| libhdf5 H5Dchunk_iter only (file+dataset already open) | 0.2 | 0 | 0 | 0 | 0 | 0 |
| native cold, block=0 | 1702.4 | 119 | 56032 | 3 | 1 | 16777216 |
| native cold, block=4194304 | 1365.3 | 2 | 8388608 | 3 | 1 | 16777216 |
| native in-driver (after libhdf5 H5Dopen), block=0 | 59.7 | 119 | 56032 | 0 | 0 | 0 |
| native in-driver (after libhdf5 H5Dopen), block=4194304 | 2.6 | 2 | 8388608 | 0 | 0 | 0 |

### Remote, HTTPS /vsicurl, default `CPL_VSIL_CURL_CHUNK_SIZE`, median of 3

| path | wall ms | VSI reads | VSI bytes | HTTP GET | HTTP HEAD | HTTP bytes |
|---|---:|---:|---:|---:|---:|---:|
| libhdf5 cold: H5Fopen+H5Dopen+H5Dchunk_iter | 1381.9 | 8 | 20972088 | 4 | 1 | 8388608 |
| libhdf5 H5Dchunk_iter only (file+dataset already open) | 0.2 | 0 | 0 | 0 | 0 | 0 |
| native cold, block=0 | 1595.5 | 119 | 56032 | 8 | 1 | 131072 |
| native cold, block=65536 | 1174.1 | 4 | 262144 | 5 | 1 | 262144 |
| native cold, block=4194304 | 1260.8 | 2 | 8388608 | 3 | 1 | 8388608 |
| native in-driver (after libhdf5 H5Dopen), block=0 | 0.2 | 119 | 56032 | 0 | 0 | 0 |
| native in-driver (after libhdf5 H5Dopen), block=65536 | 0.1 | 4 | 262144 | 0 | 0 | 0 |
| native in-driver (after libhdf5 H5Dopen), block=4194304 | 4.9 | 2 | 8388608 | 0 | 0 | 0 |

### Local file, median of 5 (largest-chunk-count dataset per product)

| product / dataset | path | wall ms | VSI reads | VSI bytes | chunks |
|---|---|---:|---:|---:|---:|
| RSLC /science/LSAR/RSLC/swaths/frequencyA/HH | libhdf5 cold: H5Fopen+H5Dopen+H5Dchunk_iter | 0.5 | 4 | 4194872 | 780 |
| RSLC /science/LSAR/RSLC/swaths/frequencyA/HH | libhdf5 H5Dchunk_iter only (file+dataset already open) | 0.1 | 0 | 0 | 780 |
| RSLC /science/LSAR/RSLC/swaths/frequencyA/HH | native cold, block=0 | 0.0 | 83 | 37128 | 780 |
| RSLC /science/LSAR/RSLC/swaths/frequencyA/HH | native cold, block=65536 | 0.0 | 3 | 196608 | 780 |
| RSLC /science/LSAR/RSLC/swaths/frequencyA/HH | native in-driver (after libhdf5 H5Dopen), block=0 | 0.1 | 83 | 37128 | 780 |
| RSLC /science/LSAR/RSLC/swaths/frequencyA/HH | native in-driver (after libhdf5 H5Dopen), block=65536 | 0.1 | 3 | 196608 | 780 |
| GCOV /science/LSAR/GCOV/grids/frequencyA/HHHH | libhdf5 cold: H5Fopen+H5Dopen+H5Dchunk_iter | 1.8 | 8 | 20972088 | 1225 |
| GCOV /science/LSAR/GCOV/grids/frequencyA/HHHH | libhdf5 H5Dchunk_iter only (file+dataset already open) | 0.1 | 0 | 0 | 1225 |
| GCOV /science/LSAR/GCOV/grids/frequencyA/HHHH | native cold, block=0 | 0.1 | 119 | 56032 | 1225 |
| GCOV /science/LSAR/GCOV/grids/frequencyA/HHHH | native cold, block=65536 | 0.1 | 4 | 262144 | 1225 |
| GCOV /science/LSAR/GCOV/grids/frequencyA/HHHH | native in-driver (after libhdf5 H5Dopen), block=0 | 0.1 | 119 | 56032 | 1225 |
| GCOV /science/LSAR/GCOV/grids/frequencyA/HHHH | native in-driver (after libhdf5 H5Dopen), block=65536 | 0.1 | 4 | 262144 | 1225 |
| GSLC /science/LSAR/GSLC/grids/frequencyA/HH | libhdf5 cold: H5Fopen+H5Dopen+H5Dchunk_iter | 0.5 | 5 | 4195009 | 1181 |
| GSLC /science/LSAR/GSLC/grids/frequencyA/HH | libhdf5 H5Dchunk_iter only (file+dataset already open) | 0.1 | 0 | 0 | 1181 |
| GSLC /science/LSAR/GSLC/grids/frequencyA/HH | native cold, block=0 | 0.1 | 114 | 54608 | 1181 |
| GSLC /science/LSAR/GSLC/grids/frequencyA/HH | native cold, block=65536 | 0.1 | 7 | 458752 | 1181 |
| GSLC /science/LSAR/GSLC/grids/frequencyA/HH | native in-driver (after libhdf5 H5Dopen), block=0 | 0.1 | 114 | 54608 | 1181 |
| GSLC /science/LSAR/GSLC/grids/frequencyA/HH | native in-driver (after libhdf5 H5Dopen), block=65536 | 0.1 | 7 | 458752 | 1181 |
| GUNW /science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/HH/coherenceMagnitude | libhdf5 cold: H5Fopen+H5Dopen+H5Dchunk_iter | 0.8 | 6 | 8389313 | 1156 |
| GUNW /science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/HH/coherenceMagnitude | libhdf5 H5Dchunk_iter only (file+dataset already open) | 0.1 | 0 | 0 | 1156 |
| GUNW /science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/HH/coherenceMagnitude | native cold, block=0 | 0.1 | 112 | 52920 | 1156 |
| GUNW /science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/HH/coherenceMagnitude | native cold, block=65536 | 0.1 | 5 | 327680 | 1156 |
| GUNW /science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/HH/coherenceMagnitude | native in-driver (after libhdf5 H5Dopen), block=0 | 0.1 | 112 | 52920 | 1156 |
| GUNW /science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/HH/coherenceMagnitude | native in-driver (after libhdf5 H5Dopen), block=65536 | 0.1 | 5 | 327680 | 1156 |

Note: the HTTP byte counter for the libhdf5 rows is lower than its 20 MiB of VSI reads; the cause was not
investigated (likely `/vsicurl/` read coalescing/caching not attributed by the network-stats counters). Treat the
VSI columns as the per-path request/byte counts.

## Findings

1. **Correctness is achievable with a small parser.** ~1.6 kLOC reproduces `H5Dchunk_iter` byte-for-byte on
   every chunked 2-D/3-D dataset of GCOV, RSLC, GSLC and GUNW, plus synthetic superblock v0/v2/v3, object header
   v2, dense groups, and all five layout-v4 chunk indexes. Real NISAR products only need superblock v2,
   object header v1, symbol-table groups and v1 B-trees.
2. **Metadata footprint drops ~80-400x.** For GCOV HHHH the native path needs 56 KB in 119 exact reads (or 256 KB
   in 4 reads with 64 KiB blocks) versus ~20 MB of 4 MiB page reads by libhdf5 open + `H5Dchunk_iter`.
3. **Wall time does not improve for `MapChunks()` alone.** Inside the driver the file is already open through
   libhdf5, so `H5Dchunk_iter` costs ~0.2 ms and issues no new I/O; the native path is equally fast only when its
   reads hit the warm `/vsicurl/` cache (0.1-0.3 ms at block <= 64 KiB, 2-60 ms otherwise). Cold, both paths are
   latency-bound by ~4-5 HTTP round trips (~1.2-1.6 s from outside us-west-2); native with 64 KiB blocks was
   ~15% faster than libhdf5 cold with default vsicurl chunking and indistinguishable with 8 MiB chunking.
4. **The benefit only materialises if libhdf5 is bypassed for open as well** (e.g. a metadata-only open path
   or a native chunk reader), because today's cost is dominated by `H5Fopen`/`H5Dopen` page-buffer reads, not by
   chunk iteration. That would require a fuller parser (attributes, strings, compound types, layout v5).
5. **Found a libhdf5 2.2.0 bug** in chunk-query APIs for extensible-array datasets whose unlimited dimension is
   not dimension 0 (details above). Not triggered by current NISAR products.

Recommendation: keep `NISAR_NATIVE_HDF5=NO` as default. A fuller parser is only justified if the goal becomes
eliminating libhdf5 metadata I/O at open (dataset discovery + chunk map) for remote access; for chunk mapping
alone it adds no measurable wall-time benefit.

## Appendix: per-dataset validation (real granules)
### NISAR_L1_PR_RSLC_004_159_A_024_2005_DHDH_A_20251109T051646_20251109T051650_X05010_N_P_J_001.h5

| dataset | shape | chunk | index | chunks | verdict | detail |
|---|---|---|---|---:|---|---|
| `/science/LSAR/RSLC/metadata/geolocationGrid/alongTrackUnitVectorX` | 20x86x337 | 1x86x337 | btree_v1 | 20 | IDENTICAL | 20 band map(s), 20 entries |
| `/science/LSAR/RSLC/metadata/geolocationGrid/alongTrackUnitVectorY` | 20x86x337 | 1x86x337 | btree_v1 | 20 | IDENTICAL | 20 band map(s), 20 entries |
| `/science/LSAR/RSLC/metadata/geolocationGrid/coordinateX` | 20x86x337 | 1x86x337 | btree_v1 | 20 | IDENTICAL | 20 band map(s), 20 entries |
| `/science/LSAR/RSLC/metadata/geolocationGrid/coordinateY` | 20x86x337 | 1x86x337 | btree_v1 | 20 | IDENTICAL | 20 band map(s), 20 entries |
| `/science/LSAR/RSLC/metadata/geolocationGrid/elevationAngle` | 20x86x337 | 1x86x337 | btree_v1 | 20 | IDENTICAL | 20 band map(s), 20 entries |
| `/science/LSAR/RSLC/metadata/geolocationGrid/groundTrackVelocity` | 20x86x337 | 1x86x337 | btree_v1 | 20 | IDENTICAL | 20 band map(s), 20 entries |
| `/science/LSAR/RSLC/metadata/geolocationGrid/incidenceAngle` | 20x86x337 | 1x86x337 | btree_v1 | 20 | IDENTICAL | 20 band map(s), 20 entries |
| `/science/LSAR/RSLC/metadata/geolocationGrid/losUnitVectorX` | 20x86x337 | 1x86x337 | btree_v1 | 20 | IDENTICAL | 20 band map(s), 20 entries |
| `/science/LSAR/RSLC/metadata/geolocationGrid/losUnitVectorY` | 20x86x337 | 1x86x337 | btree_v1 | 20 | IDENTICAL | 20 band map(s), 20 entries |
| `/science/LSAR/RSLC/swaths/frequencyA/HH` | 7600x26126 | 512x512 | btree_v1 | 780 | IDENTICAL | 1 band map(s), 780 entries |
| `/science/LSAR/RSLC/swaths/frequencyA/HV` | 7600x26126 | 512x512 | btree_v1 | 780 | IDENTICAL | 1 band map(s), 780 entries |
| `/science/LSAR/RSLC/swaths/frequencyB/HH` | 7600x6532 | 512x512 | btree_v1 | 195 | IDENTICAL | 1 band map(s), 195 entries |
| `/science/LSAR/RSLC/swaths/frequencyB/HV` | 7600x6532 | 512x512 | btree_v1 | 195 | IDENTICAL | 1 band map(s), 195 entries |

### NISAR_L2_PR_GCOV_004_159_A_024_2005_DHDH_A_20251109T051646_20251109T051650_X05010_N_P_J_001.h5

| dataset | shape | chunk | index | chunks | verdict | detail |
|---|---|---|---|---:|---|---|
| `/science/LSAR/GCOV/grids/frequencyA/HHHH` | 17640x17820 | 512x512 | btree_v1 | 1225 | IDENTICAL | 1 band map(s), 1225 entries |
| `/science/LSAR/GCOV/grids/frequencyA/HVHV` | 17640x17820 | 512x512 | btree_v1 | 1225 | IDENTICAL | 1 band map(s), 1225 entries |
| `/science/LSAR/GCOV/grids/frequencyA/mask` | 17640x17820 | 512x512 | btree_v1 | 1225 | IDENTICAL | 1 band map(s), 1225 entries |
| `/science/LSAR/GCOV/grids/frequencyA/numberOfLooks` | 17640x17820 | 512x512 | btree_v1 | 1225 | IDENTICAL | 1 band map(s), 1225 entries |
| `/science/LSAR/GCOV/grids/frequencyA/rtcGammaToSigmaFactor` | 17640x17820 | 512x512 | btree_v1 | 1225 | IDENTICAL | 1 band map(s), 1225 entries |
| `/science/LSAR/GCOV/grids/frequencyB/mask` | 4410x4455 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GCOV/grids/frequencyB/HVHV` | 4410x4455 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GCOV/grids/frequencyB/rtcGammaToSigmaFactor` | 4410x4455 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GCOV/grids/frequencyB/HHHH` | 4410x4455 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GCOV/grids/frequencyB/numberOfLooks` | 4410x4455 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/frequencyA/elevationAntennaPattern/HH` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/frequencyA/elevationAntennaPattern/HV` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/frequencyA/noiseEquivalentBackscatter/HH` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/frequencyA/noiseEquivalentBackscatter/HV` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/frequencyB/elevationAntennaPattern/HH` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/frequencyB/elevationAntennaPattern/HV` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/frequencyB/noiseEquivalentBackscatter/HH` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/frequencyB/noiseEquivalentBackscatter/HV` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/crosstalk/txHorizontalCrosspol` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/crosstalk/txVerticalCrosspol` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/crosstalk/rxVerticalCrosspol` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/calibrationInformation/crosstalk/rxHorizontalCrosspol` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/processingInformation/parameters/frequencyA/dopplerCentroid` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/processingInformation/parameters/frequencyB/dopplerCentroid` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/processingInformation/parameters/referenceTerrainHeight` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/processingInformation/timingCorrections/frequencyA/azimuthIonosphere` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/processingInformation/timingCorrections/frequencyA/slantRangeIonosphere` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/processingInformation/timingCorrections/frequencyB/azimuthIonosphere` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/processingInformation/timingCorrections/frequencyB/slantRangeIonosphere` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GCOV/metadata/radarGrid/alongTrackUnitVectorX` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GCOV/metadata/radarGrid/alongTrackUnitVectorY` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GCOV/metadata/radarGrid/elevationAngle` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GCOV/metadata/radarGrid/groundTrackVelocity` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GCOV/metadata/radarGrid/incidenceAngle` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GCOV/metadata/radarGrid/losUnitVectorX` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GCOV/metadata/radarGrid/losUnitVectorY` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GCOV/metadata/radarGrid/slantRange` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GCOV/metadata/radarGrid/zeroDopplerAzimuthTime` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |

### NISAR_L2_PR_GSLC_004_159_A_024_2005_DHDH_A_20251109T051646_20251109T051650_X05010_N_P_J_001.h5

| dataset | shape | chunk | index | chunks | verdict | detail |
|---|---|---|---|---:|---|---|
| `/science/LSAR/GSLC/grids/frequencyA/HH` | 70560x35640 | 512x512 | btree_v1 | 1181 | IDENTICAL | 1 band map(s), 9660 entries |
| `/science/LSAR/GSLC/grids/frequencyA/HV` | 70560x35640 | 512x512 | btree_v1 | 1181 | IDENTICAL | 1 band map(s), 9660 entries |
| `/science/LSAR/GSLC/grids/frequencyA/mask` | 70560x35640 | 512x512 | btree_v1 | 1181 | IDENTICAL | 1 band map(s), 9660 entries |
| `/science/LSAR/GSLC/grids/frequencyB/HH` | 70560x8910 | 512x512 | btree_v1 | 459 | IDENTICAL | 1 band map(s), 2484 entries |
| `/science/LSAR/GSLC/grids/frequencyB/HV` | 70560x8910 | 512x512 | btree_v1 | 459 | IDENTICAL | 1 band map(s), 2484 entries |
| `/science/LSAR/GSLC/grids/frequencyB/mask` | 70560x8910 | 512x512 | btree_v1 | 459 | IDENTICAL | 1 band map(s), 2484 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/frequencyA/elevationAntennaPattern/HH` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/frequencyA/elevationAntennaPattern/HV` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/frequencyA/noiseEquivalentBackscatter/HH` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/frequencyA/noiseEquivalentBackscatter/HV` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/frequencyB/elevationAntennaPattern/HH` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/frequencyB/elevationAntennaPattern/HV` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/frequencyB/noiseEquivalentBackscatter/HH` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/frequencyB/noiseEquivalentBackscatter/HV` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/crosstalk/txHorizontalCrosspol` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/crosstalk/txVerticalCrosspol` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/crosstalk/rxVerticalCrosspol` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/crosstalk/rxHorizontalCrosspol` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/geometry/sigma0` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/geometry/beta0` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/calibrationInformation/geometry/gamma0` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/processingInformation/parameters/frequencyA/dopplerCentroid` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/processingInformation/parameters/frequencyB/dopplerCentroid` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/processingInformation/parameters/referenceTerrainHeight` | 354x357 | 354x357 | btree_v1 | 1 | IDENTICAL | 1 band map(s), 1 entries |
| `/science/LSAR/GSLC/metadata/radarGrid/alongTrackUnitVectorX` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GSLC/metadata/radarGrid/alongTrackUnitVectorY` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GSLC/metadata/radarGrid/elevationAngle` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GSLC/metadata/radarGrid/groundTrackVelocity` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GSLC/metadata/radarGrid/incidenceAngle` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GSLC/metadata/radarGrid/losUnitVectorX` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GSLC/metadata/radarGrid/losUnitVectorY` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GSLC/metadata/radarGrid/slantRange` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |
| `/science/LSAR/GSLC/metadata/radarGrid/zeroDopplerAzimuthTime` | 21x718x725 | 1x512x512 | btree_v1 | 84 | IDENTICAL | 21 band map(s), 84 entries |

### NISAR_L2_PR_GUNW_003_005_A_018_004_4000_SH_20251017T125111_20251017T125127_20251029T125111_20251029T125131_X05010_N_P_J_001.h5

| dataset | shape | chunk | index | chunks | verdict | detail |
|---|---|---|---|---:|---|---|
| `/science/LSAR/GUNW/grids/frequencyA/pixelOffsets/HH/alongTrackOffset` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/pixelOffsets/HH/correlationSurfacePeak` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/pixelOffsets/HH/slantRangeOffset` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/pixelOffsets/mask` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/coherenceMagnitude` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/connectedComponents` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/ionospherePhaseScreen` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/ionospherePhaseScreenUncertainty` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/HH/unwrappedPhase` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/mask` | 4239x4311 | 512x512 | btree_v1 | 81 | IDENTICAL | 1 band map(s), 81 entries |
| `/science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/HH/coherenceMagnitude` | 16956x17244 | 512x512 | btree_v1 | 1156 | IDENTICAL | 1 band map(s), 1156 entries |
| `/science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/HH/wrappedInterferogram` | 16956x17244 | 512x512 | btree_v1 | 1156 | IDENTICAL | 1 band map(s), 1156 entries |
| `/science/LSAR/GUNW/grids/frequencyA/wrappedInterferogram/mask` | 16956x17244 | 512x512 | btree_v1 | 1156 | IDENTICAL | 1 band map(s), 1156 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/alongTrackUnitVectorX` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/alongTrackUnitVectorY` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/elevationAngle` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/groundTrackVelocity` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/hydrostaticTroposphericPhaseScreen` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/incidenceAngle` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/losUnitVectorX` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/losUnitVectorY` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/parallelBaseline` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/perpendicularBaseline` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/referenceSlantRange` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/referenceZeroDopplerAzimuthTime` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/secondarySlantRange` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/secondaryZeroDopplerAzimuthTime` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/slantRangeSolidEarthTidesPhase` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |
| `/science/LSAR/GUNW/metadata/radarGrid/wetTroposphericPhaseScreen` | 20x690x702 | 1x512x512 | btree_v1 | 80 | IDENTICAL | 20 band map(s), 80 entries |

