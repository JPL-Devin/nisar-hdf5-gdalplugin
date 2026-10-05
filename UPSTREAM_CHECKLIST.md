# Upstreaming the NISAR driver to GDAL: remaining human items

`frmts/nisar/` is laid out as a GDAL in-tree driver. A local GDAL branch (`nisar-driver`,
based on OSGeo/gdal `master`, 3.14.0dev) has been prepared and validated:

- `frmts/nisar/` = the sources from this repository, with `CMakeLists.gdal.txt` copied to `CMakeLists.txt`
- registration: `gdal_dependent_format(nisar ... "GDAL_USE_HDF5")` in `frmts/CMakeLists.txt`,
  `NISAR` in `frmts/drivers.ini`, `GDALRegister_NISAR()` in `gcore/gdal_frmts.h` and
  `frmts/gdalallregister.cpp` (registered before HDF5)
- `doc/source/drivers/raster/nisar.rst` + entry in `doc/source/drivers/raster/index.rst`
- `autotest/gdrivers/nisar.py` (synthetic fixture generated with h5py at test time)
- verified: `-DGDAL_ENABLE_DRIVER_NISAR=ON` (builtin) and `-DGDAL_ENABLE_DRIVER_NISAR_PLUGIN=ON`
  (`gdalplugins/gdal_NISAR.so`) both build, and `autotest/gdrivers/nisar.py` passes in both modes
  and against the standalone plugin on GDAL 3.12

None of this has been posted, pushed or proposed upstream. Everything below needs a person.

## 1. GDAL AI tool policy (blocking)

GDAL's [AI tool policy](https://github.com/OSGeo/gdal/blob/master/doc/source/community/ai_tool_policy.rst)
(linked from its `AGENTS.md`) bans "vibe-coded" contributions and agent-submitted PRs. LLMs may
only be used for autocomplete or repeated mechanical refactoring, any use must be disclosed, and
a human must be the primary author, responsible for every line.

- [ ] Decide how this policy applies to the driver itself and to the integration work prepared
      by Devin (registration, docs, autotest, the GetMetadata/CMake fixes). Most likely the
      upstream patch has to be rewritten or line-by-line reviewed and owned by a JPL developer.
- [ ] Disclose AI assistance in the PR, as the policy requires.

## 2. JPL legal sign-off on licensing

- [ ] SPDX header choice. All NISAR-authored C++ sources now carry
      `Copyright 2025 California Institute of Technology` / `U.S. Government sponsorship acknowledged.` /
      `SPDX-License-Identifier: Apache-2.0`, replacing the former "ALL RIGHTS RESERVED" /
      export-control / commercial-negotiation notice. Confirm JPL legal and the software release
      process (open-source release approval, export-control review) cover this.
- [ ] `hdf5vfl.{h,cpp}` derive from GDAL's MIT-licensed `frmts/hdf5/hdf5vfl.h`. They keep the
      original MIT notice (Denis Nadeau, Sam Gillingham, Even Rouault) and `SPDX-License-Identifier: MIT`,
      with a Caltech line for the modifications. Confirm this is acceptable.
- [ ] GDAL core code is MIT. Ask whether GDAL accepts an Apache-2.0 driver in-tree (some vendored
      third-party code in GDAL is Apache-2.0, see `LICENSE.TXT`), or whether the upstream copy must
      be relicensed to MIT.
- [ ] `autotest/gdrivers/nisar.py` and `nisar.rst` were written with a Caltech/MIT header, matching
      GDAL autotest convention. Confirm.

## 3. gdal-dev RFC / proposal

- [ ] Post to [gdal-dev](https://lists.osgeo.org/mailman/listinfo/gdal-dev) introducing the driver
      and asking whether the PSC wants an RFC or a plain PR. Include:
      scope (read-only, L1/L2/L3), dependencies (libhdf5 only; zlib), maintenance commitment
      (who maintains it, for how long), test data plan, licensing.
- [ ] Decide `frmts/nisar` (separate driver, as prepared) vs. an HDF5 subdriver
      (`frmts/hdf5/nisardataset.cpp`, like BAG / S-102 / S-104 / S-111). Arguments:
  - separate driver: own `NISAR:` prefix, own open options, independent plugin, reuses
    nothing from `HDF5Dataset` today;
  - HDF5 subdriver: matches how GDAL handles other HDF5 profiles, shares the HDF5 global
    lock, the `hdf5drivercore` deferred-plugin machinery and `hdf5vfl.h` (the copied VFL
    could be dropped), one HDF5 plugin to package.

## 4. Test data hosting

- [ ] The autotest currently builds a tiny GCOV-like fixture with h5py at run time and is
      skipped if h5py is missing. GDAL CI may not install h5py. Options: commit small `.h5`
      fixtures under `autotest/gdrivers/data/nisar/` with a `generate_test.py` (the
      `data/s102`, `data/bag` pattern), and/or host real subsetted granules in
      [OSGeo/gdal-test-datasets](https://github.com/OSGeo/gdal-test-datasets) for slow tests.
- [ ] Add L1 (RSLC/GCP), complex (GSLC derived subdatasets), mask, metadata-cube
      (QUANTITY/DEM_FILE) and `/vsicurl/` coverage. Confirm redistribution rights for any
      real NISAR data used.

## 5. CLA / commit authorship

- [ ] Confirm what GDAL requires from contributors (CLA, DCO or nothing) and whether Caltech
      needs to sign anything on behalf of JPL employees.
- [ ] Pick the commit author(s) (a JPL developer, real name and `@jpl.nasa.gov` address) and
      decide how to credit the original driver authors.

## 6. Technical items found during the in-tree build

These are not blocking for the standalone plugin, but GDAL reviewers will likely ask for them:

- [ ] Compiler warnings under GDAL's flags: about 130 (`-Weffc++`, `-Wold-style-cast`,
      `-Wreorder`, `-Wshadow`, unused functions/variables). GDAL CI builds with `-Werror`.
      Remove the `#pragma message` in `nisardataset.h`.
- [ ] `NisarDataset::GetMetadata()` returns `CSLDuplicate(...)` copies for `NISAR_GLOBAL` and
      `SUBDATASETS`. GDAL callers do not free the result, so these leak. Return owned members.
- [ ] Deferred plugin loading: GDAL's HDF5 driver splits `Identify()` into
      `hdf5drivercore.cpp` (`CORE_SOURCES`, `NO_SHARED_SYMBOL_WITH_CORE`) so a plugin is only
      loaded when needed. `NisarDataset::Identify()` calls `H5Fopen`, so it cannot be moved into
      core as is.
- [ ] `Identify()` opens local files with `H5Fopen` directly (not via VSI), so `/vsimem/` and
      other non-network virtual paths are not identified.
- [ ] Thread safety: GDAL's HDF5 driver uses a global HDF5 lock when libhdf5 is not
      thread-safe (`GDAL_ENABLE_HDF5_GLOBAL_LOCK`). The NISAR driver decompresses in parallel
      but has no such lock.
- [ ] `NISAR_EXPORT_ZARR` writes to `/tmp`; upstream will expect a user-given path or removal.
- [ ] The `#if GDAL_VERSION_*` compatibility blocks (GDAL < 3.12 / < 3.13) are not needed
      in-tree; keep them only in the standalone copy, or drop them once the plugin requires a
      recent GDAL.
- [ ] `DRIVER_VERSION` is only set by the standalone build (`NISAR_DRIVER_VERSION` defined in
      `frmts/nisar/CMakeLists.txt`); in-tree builds do not report it.
- [ ] `versionadded:: 3.14` in `nisar.rst` is a guess; set it to the release the driver
      lands in.

## 7. Behaviour change to release-note

- [ ] The driver no longer sets `GDAL_HTTP_MAX_RETRY=5` when it is unset (a driver must not
      change global configuration). Users reading remote granules should set it themselves.
      Add this to the next release's notes and the HOWTO version history.
