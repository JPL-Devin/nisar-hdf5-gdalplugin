#!/usr/bin/env python3
# Copyright 2025, by the California Institute of Technology.
# ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
#
# SPDX-License-Identifier: Apache-2.0

"""
Retrieve, inspect, virtualize and package one NISAR L2 GCOV granule.

Workflow (see NISAR_GCOV_Virtual_Zarr.md for full usage):

  1. Select a GCOV granule from NASA Earthdata (CMR) with earthaccess, either by a
     documented deterministic seed or explicitly (granule concept ID / granule UR),
     optionally narrowed by temporal range and bounding box. A local HDF5 copy can be
     used instead of remote reads.
  2. Walk /science/LSAR/GCOV with h5py, inventory every link and dataset (layout,
     chunking, filters, fill value, storage size) and enumerate the *allocated*
     HDF5 chunks through the chunk index (H5Dchunk_iter).
  3. Emit one consolidated Kerchunk (reference spec v1) JSON describing a Zarr v2
     hierarchy whose chunk references all point at the archive URI of the granule.
  4. Write gcov_chunk_table.{csv,md}, a dataset inventory, granule metadata,
     environment versions, a validation report and an unsupported/lossy report.
  5. Validate the manifest with fsspec + zarr + xarray (remote range requests decoded
     without GDAL or HDF5) against the source HDF5 metadata and sample values.
  6. Zip everything into nisar_gcov_virtual_zarr_<granule_id>.zip.

Pixel data is never downloaded to build references; only HDF5 metadata, the chunk
index and a few validation sample chunks are read over HTTP range requests.

Usage:
    python nisar_gcov_virtual_zarr.py generate [--seed N | --granule-ur UR |
        --concept-id G...] [--temporal START END] [--bbox W S E N]
        [--local-file PATH [--remote-uri URI]] [--output-dir DIR]
    python nisar_gcov_virtual_zarr.py retarget MANIFEST NEW_URI [--old-uri URI]
"""

import argparse
import base64
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import random
import re
import subprocess
import sys
import time
import urllib.parse
import zipfile

import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover
    h5py = None

try:
    import earthaccess
except ImportError:  # pragma: no cover
    earthaccess = None

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTPUT_DIR = os.path.join(HERE, "output")
DEFAULT_SHORT_NAME = "NISAR_L2_GCOV_BETA_V1"
DEFAULT_SEED = 20251012
DEFAULT_TEMPORAL = ("2025-10-01T00:00:00Z", "2025-10-31T23:59:59Z")
DEFAULT_MAX_CANDIDATES = 5000
GCOV_ROOT = "/science/LSAR/GCOV"
IDENT_ROOT = "/science/LSAR/identification"
TOOL_NAME = "nisar_gcov_virtual_zarr"
TOOL_VERSION = "1.0.0"

LAYOUTS = {0: "compact", 1: "contiguous", 2: "chunked", 3: "virtual"}
FILTER_NAMES = {
    1: "deflate", 2: "shuffle", 3: "fletcher32", 4: "szip", 5: "nbit",
    6: "scaleoffset", 305: "lzo", 307: "bzip2", 32000: "lzf", 32001: "blosc",
    32004: "lz4", 32008: "bitshuffle", 32013: "zfp", 32015: "zstd", 32017: "sz",
}
# Categories used for chunk-table subtotals.
CAT_PRINCIPAL = "principal_grid"
CAT_MASK = "mask"
CAT_COORD = "coordinate"
CAT_CALIB = "calibration_grid"
CAT_META = "metadata"
CAT_CUBE = "metadata_cube"
CAT_GRIDMAP = "grid_mapping"
CAT_OTHER = "other"
SUBTOTAL_ORDER = [CAT_PRINCIPAL, CAT_MASK, CAT_COORD, CAT_CALIB, CAT_META]
SUPPORT_ORDER = [CAT_COORD, CAT_GRIDMAP, CAT_CUBE]
COORD_NAMES = {"xCoordinates", "yCoordinates", "heightAboveEllipsoid"}
DIMSCALE_ATTRS = {"CLASS", "NAME", "REFERENCE_LIST", "DIMENSION_LIST", "_Netcdf4Dimid",
                  "_Netcdf4Coordinates", "_nc3_strict"}

STATUS_OK = "included"
STATUS_WARN = "included_with_warnings"
STATUS_SKIP = "skipped_unsupported"
STATUS_INVENTORY = "inventory_only"
STATUS_ERROR = "error"


# --------------------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------------------

def log(msg):
    print(f"[{TOOL_NAME}] {msg}", file=sys.stderr, flush=True)


def utcnow():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, bytes):
        return o.decode("utf-8", "replace")
    if isinstance(o, (set, tuple)):
        return list(o)
    raise TypeError(f"not JSON serializable: {type(o)}")


def dumps(obj, **kw):
    return json.dumps(obj, default=_json_default, allow_nan=False, **kw)


def write_json(path, obj):
    with open(path, "w") as f:
        f.write(dumps(obj, indent=2, sort_keys=False))
        f.write("\n")


def safe_id(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def redact_url(url):
    """Drop query string and fragment (signatures, tokens) from a URL."""
    p = urllib.parse.urlsplit(url)
    return urllib.parse.urlunsplit((p.scheme, p.netloc, p.path, "", ""))


def encode_float(v):
    v = float(v)
    if math.isnan(v):
        return "NaN"
    if math.isinf(v):
        return "Infinity" if v > 0 else "-Infinity"
    return v


def encode_fill_value(value, dtype):
    """Zarr v2 .zarray fill_value encoding for a numeric numpy dtype."""
    if value is None:
        return None
    kind = dtype.kind
    if kind == "f":
        return encode_float(value)
    if kind == "c":
        c = complex(value)
        return [encode_float(c.real), encode_float(c.imag)]
    if kind == "b":
        return bool(value)
    if kind in "iu":
        return int(value)
    return None


def encode_attr(value, issues, name):
    """HDF5 attribute value -> JSON-compatible value. Appends lossy notes to issues."""
    if h5py is not None and isinstance(value, (h5py.Reference, h5py.RegionReference)):
        issues.append(f"attribute {name}: object/region reference dropped")
        return None
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8")
        except UnicodeDecodeError:
            issues.append(f"attribute {name}: non-UTF-8 bytes decoded as latin-1")
            return value.decode("latin-1")
    if isinstance(value, str):
        return value
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        enc = encode_float(value)
        if isinstance(enc, str):
            issues.append(f"attribute {name}: non-finite float stored as JSON string {enc!r}")
        return enc
    if isinstance(value, (complex, np.complexfloating)):
        issues.append(f"attribute {name}: complex value stored as [real, imag]")
        return [encode_float(value.real), encode_float(value.imag)]
    if isinstance(value, np.void) and value.dtype.names:
        issues.append(f"attribute {name}: compound value stored as JSON object")
        return {k: encode_attr(value[k], issues, f"{name}.{k}") for k in value.dtype.names}
    if isinstance(value, np.ndarray):
        if value.dtype.names:
            issues.append(f"attribute {name}: compound array stored as list of JSON objects")
            return [encode_attr(v, issues, name) for v in value.reshape(-1)]
        if value.dtype.kind == "O" and value.size and h5py is not None and isinstance(
                value.reshape(-1)[0], (h5py.Reference, h5py.RegionReference)):
            issues.append(f"attribute {name}: object reference array dropped")
            return None
        sub = []
        local = []
        for v in value.reshape(-1):
            sub.append(encode_attr(v, local, name))
        for msg in sorted(set(local)):
            issues.append(msg)
        return np.array(sub, dtype=object).reshape(value.shape).tolist() if value.ndim else sub[0]
    issues.append(f"attribute {name}: unsupported type {type(value).__name__} dropped")
    return None


def dtype_class(dtype):
    if dtype.names:
        return "compound"
    if dtype.subdtype is not None:
        return "array"
    if h5py is not None:
        if h5py.check_string_dtype(dtype) is not None:
            return "string"
        if h5py.check_vlen_dtype(dtype) is not None:
            return "vlen"
        if h5py.check_ref_dtype(dtype) is not None:
            return "reference"
        if h5py.check_enum_dtype(dtype) is not None:
            return "enum"
    return {"f": "float", "i": "integer", "u": "integer", "c": "complex", "b": "bool",
            "S": "string", "U": "string", "O": "object", "V": "opaque"}.get(dtype.kind, dtype.kind)


# --------------------------------------------------------------------------------------
# Earthdata authentication and granule selection
# --------------------------------------------------------------------------------------

def earthdata_login(strategy="auto"):
    """Log in to Earthdata; returns the strategy used. Never records credentials."""
    if earthaccess is None:
        raise RuntimeError("earthaccess is not installed")
    auth = earthaccess.__auth__
    if auth is not None and auth.authenticated:
        return "existing-session"
    tried = []
    if strategy == "auto":
        order = []
        if os.environ.get("EARTHDATA_TOKEN") or (
                os.environ.get("EARTHDATA_USERNAME") and os.environ.get("EARTHDATA_PASSWORD")):
            order.append("environment")
        order.append("netrc")
    else:
        order = [strategy]
    for s in order:
        try:
            a = earthaccess.login(strategy=s)
            if a is not None and a.authenticated:
                return s
            tried.append(f"{s}: not authenticated")
        except Exception as e:  # noqa: BLE001
            tried.append(f"{s}: {e}")
    raise RuntimeError("Earthdata login failed (" + "; ".join(tried) + "). Provide "
                       "EARTHDATA_USERNAME/EARTHDATA_PASSWORD or EARTHDATA_TOKEN, or a "
                       "~/.netrc entry for urs.earthdata.nasa.gov.")


def _umm_links(granule):
    urls = granule["umm"].get("RelatedUrls", [])
    https = [u["URL"] for u in urls if u.get("Type") == "GET DATA"
             and u["URL"].startswith("https://") and u["URL"].endswith(".h5")]
    s3 = [u["URL"] for u in urls if u.get("Type") == "GET DATA VIA DIRECT ACCESS"
          and u["URL"].startswith("s3://") and u["URL"].endswith(".h5")]
    return https, s3


def _granule_record(g):
    umm = g["umm"]
    meta = g["meta"]
    https, s3 = _umm_links(g)
    h5_name = os.path.basename((https or s3 or [umm["GranuleUR"] + ".h5"])[0])
    archive = [a for a in umm.get("DataGranule", {}).get("ArchiveAndDistributionInformation", [])
               if a.get("Name") == h5_name]
    addl = {a["Name"]: a.get("Values") for a in umm.get("AdditionalAttributes", [])}
    return {
        "granule_ur": umm["GranuleUR"],
        "granule_concept_id": meta.get("concept-id"),
        "provider_id": meta.get("provider-id"),
        "revision_id": meta.get("revision-id"),
        "revision_date": meta.get("revision-date"),
        "collection": umm.get("CollectionReference"),
        "collection_concept_id": meta.get("collection-concept-id"),
        "filename": h5_name,
        "size_bytes": archive[0].get("SizeInBytes") if archive else None,
        "checksum": archive[0].get("Checksum") if archive else None,
        "temporal_extent": umm.get("TemporalExtent"),
        "spatial_extent": umm.get("SpatialExtent"),
        "orbit_calculated_spatial_domains": umm.get("OrbitCalculatedSpatialDomains"),
        "production_datetime": umm.get("DataGranule", {}).get("ProductionDateTime"),
        "pge_version": umm.get("PGEVersionClass"),
        "platforms": umm.get("Platforms"),
        "input_granules": umm.get("InputGranules"),
        "additional_attributes": addl,
        "https_links": https,
        "s3_links": s3,
    }


def search_candidates(args):
    """Return (granule list, query description) according to CLI selection options."""
    query = {"short_name": args.short_name}
    if args.concept_id:
        res = earthaccess.search_data(concept_id=args.concept_id, count=10)
        res = [g for g in res if g["meta"].get("concept-id") == args.concept_id]
        return res, {"concept_id": args.concept_id}
    if args.granule_ur:
        ur = args.granule_ur[:-3] if args.granule_ur.endswith(".h5") else args.granule_ur
        res = earthaccess.search_data(short_name=args.short_name, granule_name=ur, count=10)
        res = [g for g in res if g["umm"]["GranuleUR"] == ur]
        return res, {"short_name": args.short_name, "granule_ur": ur}
    kw = {"short_name": args.short_name}
    if not args.all_time:
        kw["temporal"] = tuple(args.temporal)
        query["temporal"] = list(args.temporal)
    if args.bbox:
        kw["bounding_box"] = tuple(args.bbox)
        query["bounding_box"] = list(args.bbox)
    q = earthaccess.DataGranules().short_name(args.short_name)
    if "temporal" in kw:
        q = q.temporal(*kw["temporal"])
    if "bounding_box" in kw:
        q = q.bounding_box(*kw["bounding_box"])
    hits = q.hits()
    query["cmr_hits"] = hits
    if hits > args.max_candidates:
        raise RuntimeError(
            f"{hits} candidate granules exceed --max-candidates {args.max_candidates}; "
            "narrow --temporal/--bbox so the candidate set (and therefore the seeded "
            "selection) stays reproducible.")
    res = earthaccess.search_data(count=-1, **kw)
    return res, query


def select_granule(args):
    """Resolve the granule to process; returns a metadata dict (no secrets)."""
    sel = {"selection": {}, "auth_strategy": None}
    if args.local_file and not (args.granule_ur or args.concept_id):
        sel["selection"] = {"method": "local-file", "local_file": os.path.abspath(args.local_file)}
        sel["granule"] = None
        return sel
    sel["auth_strategy"] = earthdata_login(args.auth)
    log(f"Earthdata login via {sel['auth_strategy']}")
    res, query = search_candidates(args)
    if not res:
        raise RuntimeError(f"no GCOV granule matched {query}")
    res = sorted(res, key=lambda g: g["umm"]["GranuleUR"])
    urs = [g["umm"]["GranuleUR"] for g in res]
    if args.concept_id or args.granule_ur:
        idx = 0
        method = "concept-id" if args.concept_id else "granule-ur"
    else:
        idx = random.Random(args.seed).randrange(len(res))
        method = "seeded-random"
    sel["selection"] = {
        "method": method,
        "seed": args.seed if method == "seeded-random" else None,
        "rng": "python random.Random(seed).randrange(len(candidates))",
        "candidate_order": "sorted by GranuleUR (ascending)",
        "query": query,
        "candidate_count": len(res),
        "candidate_list_sha256": hashlib.sha256("\n".join(urs).encode()).hexdigest(),
        "selected_index": idx,
        "local_file": os.path.abspath(args.local_file) if args.local_file else None,
    }
    sel["granule"] = _granule_record(res[idx])
    sel["_granule_obj"] = res[idx]
    return sel


def probe_https_access(url):
    """Characterise the redirect chain of an Earthdata HTTPS URL without storing signatures."""
    import requests
    out = {"stable_uri": url}
    try:
        session = earthaccess.get_requests_https_session()
        r = session.get(url, headers={"Range": "bytes=0-7"}, allow_redirects=False, timeout=60)
        out["initial_status"] = r.status_code
        loc = r.headers.get("Location")
        if r.status_code in (301, 302, 303, 307, 308) and loc:
            p = urllib.parse.urlsplit(loc)
            q = urllib.parse.parse_qs(p.query)
            out["redirect_host"] = p.netloc
            out["redirect_url_redacted"] = redact_url(loc)
            out["redirect_query_parameters"] = sorted(q)
            if "urs.earthdata.nasa.gov" in p.netloc:
                out["access_type"] = "Earthdata Login OAuth redirect (request not authenticated)"
            elif p.netloc.endswith("cloudfront.net") and "Signature" in q:
                out["access_type"] = "CloudFront signed URL (HTTP redirect from the DAAC TEA endpoint)"
            elif "X-Amz-Signature" in q:
                out["access_type"] = "S3 presigned HTTPS URL (HTTP redirect)"
            else:
                out["access_type"] = "HTTP redirect to unsigned URL"
            if "Expires" in q:
                exp = int(q["Expires"][0])
                out["signed_url_expires_utc"] = dt.datetime.fromtimestamp(
                    exp, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
                out["signed_url_lifetime_s"] = exp - int(time.time())
            elif "X-Amz-Expires" in q:
                out["signed_url_lifetime_s"] = int(q["X-Amz-Expires"][0])
            r2 = requests.get(loc, headers={"Range": "bytes=0-7"}, timeout=60)
        else:
            out["access_type"] = "direct HTTPS (no redirect)"
            r2 = r
        out["final_status"] = r2.status_code
        out["hdf5_signature_ok"] = r2.content[:8] == b"\x89HDF\r\n\x1a\n"
        cr = r2.headers.get("Content-Range", "")
        if "/" in cr:
            out["content_length"] = int(cr.split("/")[-1])
        out["final_server"] = r2.headers.get("Server")
        out["cloudfront"] = "cloudfront" in (r2.headers.get("Via", "") + r2.headers.get("X-Cache", "")).lower()
        out["etag"] = r2.headers.get("ETag")
        out["last_modified"] = r2.headers.get("Last-Modified")
    except Exception as e:  # noqa: BLE001
        out["error"] = f"{type(e).__name__}: {e}"
    return out


# --------------------------------------------------------------------------------------
# Opening the source HDF5
# --------------------------------------------------------------------------------------

class Source:
    """HDF5 handle plus what is needed to read the reference target through fsspec."""

    def __init__(self, h5, reference_uri, access, remote_protocol, remote_options,
                 local_file=None, fileobj=None):
        self.h5 = h5
        self.reference_uri = reference_uri
        self.access = access
        self.remote_protocol = remote_protocol
        self.remote_options = remote_options
        self.local_file = local_file
        self.fileobj = fileobj

    def close(self):
        try:
            self.h5.close()
        finally:
            if self.fileobj is not None:
                self.fileobj.close()


def open_source(sel, args):
    g = sel.get("granule")
    if args.local_file:
        path = os.path.abspath(args.local_file)
        h5 = h5py.File(path, "r")
        if args.remote_uri:
            uri = args.remote_uri
        elif g and g["https_links"]:
            uri = g["https_links"][0]
        else:
            uri = "file://" + path
        proto, opts = _remote_fs_options(uri, sel)
        return Source(h5, uri, "local-file", proto, opts, local_file=path)
    if args.access == "direct":
        uri = g["s3_links"][0]
        fs = earthaccess.get_s3_filesystem(results=[sel["_granule_obj"]])
        proto = "s3"
        opts = {k: v for k, v in fs.storage_options.items() if k in ("key", "secret", "token")}
    else:
        if not g["https_links"]:
            raise RuntimeError("granule has no HTTPS data link")
        uri = g["https_links"][0]
        fs = earthaccess.get_fsspec_https_session()
        proto, opts = _remote_fs_options(uri, sel)
    fobj = fs.open(uri, mode="rb", block_size=args.block_size, cache_type="blockcache")
    h5 = h5py.File(fobj, "r")
    return Source(h5, uri, args.access, proto, opts, fileobj=fobj)


def _remote_fs_options(uri, sel):
    """fsspec protocol + options for resolving references to ``uri`` (in memory only)."""
    scheme = urllib.parse.urlsplit(uri).scheme or "file"
    if scheme in ("http", "https"):
        opts = {}
        if earthaccess is not None and earthaccess.__auth__ is not None and earthaccess.__auth__.authenticated:
            fs = earthaccess.get_fsspec_https_session()
            opts = {"client_kwargs": dict(fs.client_kwargs)}
        return "https", opts
    if scheme == "s3":
        return "s3", {}
    return "file", {}


# --------------------------------------------------------------------------------------
# Inventory
# --------------------------------------------------------------------------------------

def _filters(dcpl):
    out = []
    for i in range(dcpl.get_nfilters()):
        fid, flags, cd, name = dcpl.get_filter(i)
        label = FILTER_NAMES.get(int(fid), name.decode(errors="replace") or f"filter-{fid}")
        out.append({"id": int(fid), "name": label,
                    "flags": int(flags), "cd_values": [int(c) for c in cd]})
    return out


def pipeline_str(filters):
    if not filters:
        return "none"
    return "|".join(f"{f['name']}({','.join(str(c) for c in f['cd_values'])})" for f in filters)


def _enumerate_chunks(dsid):
    chunks = []
    try:
        dsid.chunk_iter(lambda si: chunks.append((tuple(si.chunk_offset), int(si.filter_mask),
                                                  int(si.byte_offset), int(si.size))))
        return chunks, "H5Dchunk_iter"
    except Exception:  # noqa: BLE001
        chunks = []
        for i in range(dsid.get_num_chunks()):
            si = dsid.get_chunk_info(i)
            chunks.append((tuple(si.chunk_offset), int(si.filter_mask), int(si.byte_offset), int(si.size)))
        return chunks, "H5Dget_chunk_info"


def categorize(rel, rank, is_dimscale, is_gridmap):
    base = rel.rsplit("/", 1)[-1]
    if rank == 2:
        if "mask" in base.lower():
            return CAT_MASK
        if is_dimscale or base in COORD_NAMES:
            return CAT_COORD
        if rel.startswith("grids/"):
            return CAT_PRINCIPAL
        if rel.startswith("metadata/calibrationInformation/"):
            return CAT_CALIB
        return CAT_META
    if rank == 3:
        return CAT_CUBE
    if is_dimscale or base in COORD_NAMES:
        return CAT_COORD
    if is_gridmap:
        return CAT_GRIDMAP
    return CAT_OTHER


def inventory(h5, root=GCOV_ROOT):
    """Inventory every link under ``root``; returns (datasets, links, groups)."""
    grp = h5[root]
    links = []
    names = []
    grp.visit_links(names.append)
    seen_obj = {}
    datasets = {}
    groups = {root: {"path": root, "attrs": dict(grp.attrs.items())}}
    for name in sorted(names):
        path = f"{root}/{name}"
        lk = grp.get(name, getlink=True)
        if isinstance(lk, h5py.ExternalLink):
            links.append({"path": path, "type": "external", "target": f"{lk.filename}:{lk.path}"})
            continue
        if isinstance(lk, h5py.SoftLink):
            links.append({"path": path, "type": "soft", "target": lk.path})
            continue
        obj = grp[name]
        info = h5py.h5o.get_info(obj.id)
        key = (info.fileno, info.addr)
        if key in seen_obj:
            links.append({"path": path, "type": "hard-alias", "target": seen_obj[key]})
            if seen_obj[key] in datasets:
                datasets[seen_obj[key]]["aliases"].append(path)
            continue
        seen_obj[key] = path
        if isinstance(obj, h5py.Group):
            groups[path] = {"path": path, "attrs": dict(obj.attrs.items())}
            continue
        if isinstance(obj, h5py.Dataset):
            datasets[path] = describe_dataset(obj, path, root)
    for lk in links:
        if lk["type"] in ("soft", "hard-alias"):
            tgt = lk["target"]
            if tgt in datasets:
                datasets[tgt]["issues"].append(
                    {"kind": f"{lk['type']}-link", "severity": "info",
                     "detail": f"also reachable as {lk['path']} ({lk['type']} link); "
                               "manifest keeps the primary path only"})
    mark_dimscales(h5, datasets)
    return datasets, links, groups


def describe_dataset(d, path, root):
    dcpl = d.id.get_create_plist()
    layout = LAYOUTS.get(dcpl.get_layout(), str(dcpl.get_layout()))
    filters = _filters(dcpl) if layout == "chunked" else []
    shape = tuple(int(s) for s in d.shape) if d.shape is not None else None
    rec = {
        "path": path,
        "rel": path[len(root) + 1:],
        "rank": len(shape) if shape is not None else None,
        "shape": list(shape) if shape is not None else None,
        "dtype": d.dtype.str if not d.dtype.names else "compound",
        "dtype_class": dtype_class(d.dtype),
        "itemsize": d.dtype.itemsize,
        "layout": layout,
        "chunk_shape": list(d.chunks) if d.chunks else None,
        "filters": filters,
        "pipeline": pipeline_str(filters),
        "shuffle": any(f["id"] == 2 for f in filters),
        "fill_value_defined": int(dcpl.fill_value_defined()),
        "fill_value": None,
        "storage_size": int(d.id.get_storage_size()),
        "external_storage_files": int(dcpl.get_external_count()),
        "logical_chunks": None,
        "allocated_chunks": None,
        "chunk_index_api": None,
        "chunks": None,
        "contiguous_offset": None,
        "aliases": [],
        "issues": [],
        "dims": None,
        "dimscale_of": [],
        "is_dimscale": False,
        "gridmap_of": [],
        "category": None,
        "zarr_path": None,
        "status": STATUS_INVENTORY,
    }
    try:
        fv = d.fillvalue
        if d.dtype.kind in "fciub":
            rec["fill_value"] = encode_fill_value(fv, d.dtype)
        else:
            rec["fill_value"] = repr(fv)
    except Exception as e:  # noqa: BLE001
        rec["issues"].append({"kind": "fill-value", "severity": "warning", "detail": str(e)})
    attrs = {}
    for k in d.attrs.keys():
        try:
            attrs[k] = d.attrs[k]
        except Exception as e:  # noqa: BLE001
            rec["issues"].append({"kind": "attribute-read", "severity": "warning",
                                  "detail": f"attribute {k} unreadable: {e}"})
    rec["_attrs"] = attrs
    if attrs.get("CLASS") in (b"DIMENSION_SCALE", "DIMENSION_SCALE"):
        rec["is_dimscale"] = True
    if layout == "chunked":
        chunks, api = _enumerate_chunks(d.id)
        rec["chunks"] = chunks
        rec["chunk_index_api"] = api
        rec["allocated_chunks"] = len(chunks)
        rec["logical_chunks"] = math.prod(math.ceil(s / c) if s else 0 for s, c in zip(shape, d.chunks)) if shape else 1
        rec["compressed_bytes"] = sum(c[3] for c in chunks)
    elif layout == "contiguous":
        off = d.id.get_offset()
        rec["contiguous_offset"] = int(off) if off is not None else None
        rec["logical_chunks"] = 1
        rec["allocated_chunks"] = 1 if off is not None and rec["storage_size"] > 0 else 0
        rec["compressed_bytes"] = rec["storage_size"] if rec["allocated_chunks"] else 0
    elif layout == "compact":
        rec["logical_chunks"] = 1
        rec["allocated_chunks"] = 1
        rec["compressed_bytes"] = rec["storage_size"]
    else:
        rec["logical_chunks"] = 0
        rec["allocated_chunks"] = 0
        rec["compressed_bytes"] = 0
    return rec


def mark_dimscales(h5, datasets):
    for rec in datasets.values():
        dl = rec["_attrs"].get("DIMENSION_LIST")
        gm = rec["_attrs"].get("grid_mapping")
        if dl is not None and rec["rank"]:
            dims = []
            for i in range(rec["rank"]):
                try:
                    refs = dl[i]
                    tgt = h5[refs[0]].name if len(refs) else None
                except Exception:  # noqa: BLE001
                    tgt = None
                dims.append(tgt)
            rec["dims"] = dims
            for t in dims:
                if t in datasets:
                    datasets[t]["dimscale_of"].append(rec["path"])
        if gm is not None:
            gname = gm.decode() if isinstance(gm, bytes) else str(gm)
            gpath = rec["path"].rsplit("/", 1)[0] + "/" + gname
            rec["grid_mapping_path"] = gpath
            if gpath in datasets:
                datasets[gpath]["gridmap_of"].append(rec["path"])
    for rec in datasets.values():
        rec["category"] = categorize(rec["rel"], rec["rank"], rec["is_dimscale"] or bool(rec["dimscale_of"]),
                                     bool(rec["gridmap_of"]))


# --------------------------------------------------------------------------------------
# Manifest planning and generation
# --------------------------------------------------------------------------------------

def zarr_codecs(rec):
    """Map the HDF5 filter pipeline to (compressor, filters, unsupported reasons)."""
    mapped = []
    bad = []
    warn = []
    for f in rec["filters"]:
        if f["id"] == 1:
            mapped.append({"id": "zlib", "level": f["cd_values"][0] if f["cd_values"] else 1})
        elif f["id"] == 2:
            mapped.append({"id": "shuffle", "elementsize": rec["itemsize"]})
        elif f["id"] == 3:
            mapped.append({"id": "fletcher32"})
            warn.append("fletcher32 checksum mapped to numcodecs fletcher32")
        else:
            bad.append(f"HDF5 filter {f['name']} (id {f['id']}) has no portable Zarr v2 codec mapping")
    compressor = None
    if mapped and mapped[-1]["id"] == "zlib":
        compressor = mapped.pop()
    return compressor, (mapped or None), bad, warn


def plan_manifest(datasets, include_rank3=True):
    """Decide which datasets enter the manifest and why others do not."""
    for rec in datasets.values():
        reasons = []
        warns = []
        cls = rec["dtype_class"]
        wanted = rec["rank"] == 2 or (rec["rank"] == 3 and include_rank3)
        supporting = bool(rec["dimscale_of"]) or bool(rec["gridmap_of"])
        if not (wanted or supporting):
            if rec["rank"] == 3:
                rec["issues"].append({"kind": "rank-3", "severity": "info",
                                      "detail": "rank-3 array excluded (--no-rank3)"})
            continue
        if rec["layout"] == "virtual":
            reasons.append("virtual dataset (VDS) layout")
        if rec["external_storage_files"]:
            reasons.append("raw data stored in external files")
        if cls in ("compound", "string", "vlen", "reference", "opaque", "object", "array"):
            reasons.append(f"{cls} datatype is not representable as a portable Zarr v2 numeric array")
        if rec["layout"] == "chunked":
            _, _, bad, w = zarr_codecs(rec)
            reasons += bad
            warns += w
            fm = [c for c in rec["chunks"] if c[1] != 0]
            if fm:
                reasons.append(f"{len(fm)} chunks have a non-zero HDF5 filter mask (per-chunk skipped filters)")
            missing = rec["logical_chunks"] - rec["allocated_chunks"]
            if missing > 0:
                rec["issues"].append({"kind": "missing-chunks", "severity": "info",
                                      "detail": f"{missing} of {rec['logical_chunks']} logical chunks are "
                                                "not allocated; "
                                                f"readers return fill_value {rec['fill_value']!r} (same as HDF5)"})
        if rec["layout"] == "contiguous" and rec["allocated_chunks"] == 0 and rec["shape"] and math.prod(rec["shape"]):
            rec["issues"].append({"kind": "missing-chunks", "severity": "info",
                                  "detail": "contiguous storage not allocated; readers return fill_value"})
        if cls == "complex":
            warns.append("complex-valued array: Zarr v2 '<c8'/'<c16' dtype is readable by zarr-python/xarray "
                         "but not by all Zarr implementations")
        if cls == "enum":
            warns.append("HDF5 enum stored as its integer base type; enum mapping not preserved")
        if rec["aliases"]:
            warns.append("dataset has additional hard links: " + ", ".join(rec["aliases"]))
        if reasons:
            rec["status"] = STATUS_SKIP
            for r in reasons:
                rec["issues"].append({"kind": "unsupported", "severity": "error", "detail": r})
        else:
            rec["status"] = STATUS_WARN if warns else STATUS_OK
            rec["zarr_path"] = rec["rel"]
        for w in warns:
            rec["issues"].append({"kind": "lossy-or-portability", "severity": "warning", "detail": w})


def zattrs_for(rec, datasets, root):
    issues = []
    out = {}
    for k, v in rec["_attrs"].items():
        if k in DIMSCALE_ATTRS:
            continue
        enc = encode_attr(v, issues, k)
        if enc is not None:
            out[k] = enc
    dims = []
    for i in range(rec["rank"] or 0):
        tgt = rec["dims"][i] if rec["dims"] else None
        if tgt and tgt.startswith(root + "/"):
            dims.append(tgt.rsplit("/", 1)[-1])
        elif rec["rank"] == 1 and rec["dimscale_of"]:
            dims.append(rec["rel"].rsplit("/", 1)[-1])
        else:
            dims.append(f"{rec['rel'].rsplit('/', 1)[-1]}_dim{i}")
    out["_ARRAY_DIMENSIONS"] = dims
    dropped = [k for k in rec["_attrs"] if k in DIMSCALE_ATTRS]
    if dropped:
        issues.append("HDF5 dimension-scale attributes " + ", ".join(sorted(dropped)) +
                      " replaced by _ARRAY_DIMENSIONS")
    out["_hdf5_path"] = rec["path"]
    return out, issues


def group_attrs(attrs, where):
    issues = []
    out = {}
    for k, v in attrs.items():
        enc = encode_attr(v, issues, k)
        if enc is not None:
            out[k] = enc
    return out, [f"{where}: {i}" for i in issues]


def read_identification(h5):
    out = {}
    if IDENT_ROOT not in h5:
        return out
    g = h5[IDENT_ROOT]
    for k in sorted(g.keys()):
        obj = g[k]
        if isinstance(obj, h5py.Dataset) and obj.size is not None and obj.size <= 64:
            try:
                out[k] = encode_attr(obj[()], [], k)
            except Exception:  # noqa: BLE001
                pass
    return out


def build_manifest(src, datasets, groups, sel):
    """Return (refs dict, per-dataset JSON byte estimate, attribute issues by path)."""
    h5 = src.h5
    root = GCOV_ROOT
    refs = {}
    per_bytes = {}
    attr_issues = {}
    included = sorted((r for r in datasets.values() if r["zarr_path"]), key=lambda r: r["path"])
    # Group hierarchy: every ancestor of an included array needs .zgroup (+ .zattrs).
    needed = {""}
    for r in included:
        parts = r["zarr_path"].split("/")[:-1]
        for i in range(1, len(parts) + 1):
            needed.add("/".join(parts[:i]))
    root_attrs, gi = group_attrs(dict(h5.attrs.items()), "/")
    if root in groups:
        ga, gi2 = group_attrs(groups[root]["attrs"], root)
        root_attrs.update(ga)
        gi += gi2
    root_attrs["nisar_identification"] = read_identification(h5)
    g = sel.get("granule") or {}
    root_attrs["nisar_virtual_zarr"] = {
        "generator": f"{TOOL_NAME} {TOOL_VERSION}",
        "created_utc": utcnow(),
        "hdf5_root": root,
        "source_uri": src.reference_uri,
        "granule_ur": g.get("granule_ur"),
        "granule_concept_id": g.get("granule_concept_id"),
        "collection": g.get("collection"),
        "filename": g.get("filename") or os.path.basename(src.reference_uri),
        "note": "Chunk references point at source_uri; HTTPS archive URLs require Earthdata Login.",
    }
    if gi:
        attr_issues["/"] = gi
    for gp in sorted(needed):
        key = f"{gp}/" if gp else ""
        refs[f"{key}.zgroup"] = dumps({"zarr_format": 2})
        if gp:
            hp = f"{root}/{gp}"
            ga, gis = group_attrs(groups.get(hp, {"attrs": {}})["attrs"], hp)
            if gis:
                attr_issues[hp] = gis
            refs[f"{key}.zattrs"] = dumps(ga)
        else:
            refs[".zattrs"] = dumps(root_attrs)
    for r in included:
        zp = r["zarr_path"]
        shape = r["shape"]
        dtype = np.dtype(r["dtype"])
        if r["layout"] == "chunked":
            chunks = r["chunk_shape"]
            compressor, filters, _, _ = zarr_codecs(r)
        else:
            chunks = [max(1, s) for s in shape]
            compressor, filters = None, None
        zarray = {"zarr_format": 2, "shape": shape, "chunks": chunks, "dtype": dtype.str,
                  "compressor": compressor, "filters": filters,
                  "fill_value": r["fill_value"] if r["fill_value_defined"] or r["fill_value"] is not None else None,
                  "order": "C", "dimension_separator": "."}
        zattrs, ai = zattrs_for(r, datasets, root)
        if ai:
            attr_issues[r["path"]] = ai
        local = {f"{zp}/.zarray": dumps(zarray), f"{zp}/.zattrs": dumps(zattrs)}
        uri = src.reference_uri
        if r["layout"] == "chunked":
            for off, _fm, byte_off, size in r["chunks"]:
                key = ".".join(str(o // c) for o, c in zip(off, chunks)) if shape else "0"
                local[f"{zp}/{key}"] = [uri, byte_off, size]
        elif r["layout"] == "contiguous":
            if r["allocated_chunks"]:
                key = ".".join("0" for _ in shape) if shape else "0"
                local[f"{zp}/{key}"] = [uri, r["contiguous_offset"], r["storage_size"]]
        elif r["layout"] == "compact":
            arr = np.ascontiguousarray(h5[r["path"]][()], dtype=dtype)
            key = ".".join("0" for _ in shape) if shape else "0"
            local[f"{zp}/{key}"] = "base64:" + base64.b64encode(arr.tobytes()).decode()
            r["issues"].append({"kind": "inline-data", "severity": "info",
                                "detail": "compact-layout dataset stored inline (base64) in the manifest"})
        refs.update(local)
        per_bytes[r["path"]] = sum(len(dumps(k)) + len(dumps(v)) + 4 for k, v in local.items())
    # Consolidated metadata so zarr/xarray can open without listing.
    meta = {k: json.loads(v) for k, v in refs.items()
            if k.rsplit("/", 1)[-1] in (".zgroup", ".zattrs", ".zarray")}
    refs[".zmetadata"] = dumps({"zarr_consolidated_format": 1, "metadata": meta})
    return refs, per_bytes, attr_issues


def check_dim_consistency(datasets, refs):
    """Same dimension name within one Zarr group must have one length (xarray requirement)."""
    problems = []
    by_group = {}
    for r in datasets.values():
        if not r["zarr_path"]:
            continue
        za = json.loads(refs[f"{r['zarr_path']}/.zattrs"])
        grp = r["zarr_path"].rsplit("/", 1)[0] if "/" in r["zarr_path"] else ""
        for name, size in zip(za["_ARRAY_DIMENSIONS"], r["shape"]):
            prev = by_group.setdefault(grp, {}).setdefault(name, (size, r["path"]))
            if prev[0] != size:
                problems.append(f"group {grp or '/'}: dimension {name} is {prev[0]} in {prev[1]} "
                                f"but {size} in {r['path']}")
    return problems


# --------------------------------------------------------------------------------------
# Chunk table
# --------------------------------------------------------------------------------------

TABLE_COLUMNS = ["row_type", "hdf5_path", "zarr_path", "category", "rank", "dimensions", "dtype",
                 "chunk_shape", "logical_chunks", "allocated_chunks", "skipped_chunks",
                 "compressed_chunk_bytes", "approx_json_bytes", "filter_pipeline", "manifest_status",
                 "counted_as_2d"]


def _dims_str(r):
    names = []
    for i, s in enumerate(r["shape"] or []):
        t = r["dims"][i] if r["dims"] and r["dims"][i] else None
        names.append(f"{t.rsplit('/', 1)[-1]}[{s}]" if t else f"[{s}]")
    return " x ".join(names) if names else "scalar"


def _table_row(r, per_bytes):
    alloc = r["allocated_chunks"] or 0
    ref_chunks = alloc if r["zarr_path"] else 0
    return {
        "row_type": "dataset",
        "hdf5_path": r["path"],
        "zarr_path": r["zarr_path"] or "",
        "category": r["category"],
        "rank": r["rank"],
        "dimensions": _dims_str(r),
        "dtype": r["dtype"],
        "chunk_shape": "x".join(str(c) for c in r["chunk_shape"]) if r["chunk_shape"] else r["layout"],
        "logical_chunks": r["logical_chunks"],
        "allocated_chunks": alloc,
        "skipped_chunks": (r["logical_chunks"] or 0) - ref_chunks,
        "compressed_chunk_bytes": r["compressed_bytes"],
        "approx_json_bytes": per_bytes.get(r["path"], 0),
        "filter_pipeline": r["pipeline"],
        "manifest_status": r["status"],
        "counted_as_2d": "yes" if r["rank"] == 2 else "no",
    }


def _sum_rows(rows, row_type, label):
    return {
        "row_type": row_type, "hdf5_path": label, "zarr_path": "", "category": "",
        "rank": "", "dimensions": f"{len(rows)} datasets", "dtype": "", "chunk_shape": "",
        "logical_chunks": sum(r["logical_chunks"] or 0 for r in rows),
        "allocated_chunks": sum(r["allocated_chunks"] for r in rows),
        "skipped_chunks": sum(r["skipped_chunks"] for r in rows),
        "compressed_chunk_bytes": sum(r["compressed_chunk_bytes"] or 0 for r in rows),
        "approx_json_bytes": sum(r["approx_json_bytes"] for r in rows),
        "filter_pipeline": "", "manifest_status": "",
        "counted_as_2d": "",
    }


def chunk_table(datasets, per_bytes, manifest_bytes):
    two_d = [_table_row(r, per_bytes) for r in sorted(datasets.values(), key=lambda r: r["path"]) if r["rank"] == 2]
    support = [_table_row(r, per_bytes) for r in sorted(datasets.values(), key=lambda r: r["path"])
               if r["rank"] != 2 and r["zarr_path"]]
    subtotals = [_sum_rows([r for r in two_d if r["category"] == c], "subtotal", f"SUBTOTAL {c}")
                 for c in SUBTOTAL_ORDER]
    total = _sum_rows(two_d, "total", "TOTAL 2-D datasets")
    total["manifest_status"] = (f"{sum(1 for r in two_d if r['manifest_status'] in (STATUS_OK, STATUS_WARN))} "
                                f"in manifest / {len(two_d)}")
    sup_sub = [_sum_rows([r for r in support if r["category"] == c], "support_subtotal", f"SUPPORT {c}")
               for c in SUPPORT_ORDER]
    manifest = _sum_rows(two_d + support, "manifest", "MANIFEST (all arrays incl. support + metadata keys)")
    manifest["approx_json_bytes"] = manifest_bytes
    manifest["dimensions"] = f"{len(two_d) + len(support)} arrays"
    return two_d, subtotals, total, support, sup_sub, manifest


def write_chunk_table(outdir, tables, sel_summary):
    two_d, subtotals, total, support, sup_sub, manifest = tables
    csv_path = os.path.join(outdir, "gcov_chunk_table.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=TABLE_COLUMNS)
        w.writeheader()
        for r in two_d + subtotals + [total] + support + sup_sub + [manifest]:
            w.writerow(r)
    md_path = os.path.join(outdir, "gcov_chunk_table.md")
    cols = ["hdf5_path", "zarr_path", "category", "dimensions", "dtype", "chunk_shape", "logical_chunks",
            "allocated_chunks", "skipped_chunks", "compressed_chunk_bytes", "approx_json_bytes",
            "filter_pipeline", "manifest_status"]

    def md_table(rows, columns):
        lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
        for r in rows:
            vals = []
            for c in columns:
                v = r[c]
                vals.append(f"`{v}`" if c in ("hdf5_path", "zarr_path") and v and not str(v).isupper() else str(v))
            lines.append("| " + " | ".join(vals) + " |")
        return "\n".join(lines)

    sum_cols = ["hdf5_path", "dimensions", "logical_chunks", "allocated_chunks", "skipped_chunks",
                "compressed_chunk_bytes", "approx_json_bytes", "manifest_status"]
    with open(md_path, "w") as f:
        f.write(f"# NISAR GCOV chunk table\n\n{sel_summary}\n\n")
        f.write("## 2-D datasets (one row per rank-2 HDF5 dataset under `/science/LSAR/GCOV`)\n\n")
        f.write(md_table(two_d, cols) + "\n\n")
        f.write("## Subtotals and summary (2-D datasets)\n\n")
        f.write(md_table(subtotals + [total], sum_cols) + "\n\n")
        f.write("`skipped_chunks` = logical chunks with no reference in the manifest (never allocated in HDF5, "
                "or belonging to a skipped dataset). `approx_json_bytes` = serialized size of the dataset's "
                "`.zarray`, `.zattrs` and chunk-reference entries.\n\n")
        f.write("## Supporting arrays in the manifest (not counted as 2-D subdatasets)\n\n")
        f.write(md_table(support, cols) + "\n\n")
        f.write(md_table(sup_sub + [manifest], sum_cols) + "\n")
    return csv_path, md_path


def write_inventory(outdir, datasets, links):
    path = os.path.join(outdir, "gcov_dataset_inventory.csv")
    cols = ["hdf5_path", "rank", "shape", "dtype", "dtype_class", "layout", "chunk_shape",
            "allocated_chunks", "logical_chunks", "fill_value", "fill_value_defined", "filter_pipeline",
            "shuffle", "storage_size", "chunk_index_api", "category", "zarr_path", "manifest_status",
            "issues"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in sorted(datasets.values(), key=lambda r: r["path"]):
            w.writerow([r["path"], r["rank"], "x".join(map(str, r["shape"] or [])) or "scalar", r["dtype"],
                        r["dtype_class"], r["layout"],
                        "x".join(map(str, r["chunk_shape"])) if r["chunk_shape"] else "",
                        r["allocated_chunks"], r["logical_chunks"], dumps(r["fill_value"]),
                        r["fill_value_defined"], r["pipeline"], r["shuffle"], r["storage_size"],
                        r["chunk_index_api"] or "", r["category"], r["zarr_path"] or "", r["status"],
                        "; ".join(i["detail"] for i in r["issues"])])
        for lk in links:
            w.writerow([lk["path"], "", "", "", f"link:{lk['type']}", "", "", "", "", "", "", "", "", "",
                        "", "", "", STATUS_SKIP if lk["type"] != "hard-alias" else STATUS_INVENTORY,
                        f"{lk['type']} link -> {lk['target']}"])
    return path


# --------------------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------------------

def _expected_codecs(dset):
    dcpl = dset.id.get_create_plist()
    rec = {"filters": _filters(dcpl) if dcpl.get_layout() == 2 else [], "itemsize": dset.dtype.itemsize}
    comp, filt, _, _ = zarr_codecs(rec)
    return comp, filt


def _codec_ids(codec_json):
    if codec_json is None:
        return None
    if isinstance(codec_json, list):
        return [c["id"] for c in codec_json] or None
    return codec_json["id"]


def _sample_chunk_indices(rec, zchunks, n):
    grid = [math.ceil(s / c) for s, c in zip(rec["shape"], zchunks)]
    if rec["layout"] != "chunked":
        return [tuple(0 for _ in grid)], []
    alloc = sorted(tuple(o // c for o, c in zip(off, zchunks)) for off, _, _, _ in rec["chunks"])
    picks = []
    if alloc:
        for i in sorted({0, len(alloc) // 2, len(alloc) - 1})[:max(1, n)]:
            picks.append(alloc[i])
    missing = []
    if rec["logical_chunks"] > rec["allocated_chunks"]:
        aset = set(alloc)
        for idx in np.ndindex(*grid):
            if idx not in aset:
                missing.append(idx)
                break
    return picks, missing


def _region(idx, zchunks, shape):
    return tuple(slice(i * c, min((i + 1) * c, s)) for i, c, s in zip(idx, zchunks, shape))


def validate(src, refs, datasets, manifest_path, args):
    """Validate the manifest; every failure is recorded, nothing is dropped silently."""
    import fsspec
    import zarr
    report = {"started_utc": utcnow(), "reference_uri": src.reference_uri, "arrays": [],
              "structure": {}, "xarray": {}, "raw_decode": {}, "failures": []}
    fail = report["failures"]
    h5 = src.h5

    # Structural checks on the reference JSON itself.
    key_re = re.compile(r"^([A-Za-z0-9_.-]+/)*(\.zgroup|\.zattrs|\.zarray|\.zmetadata|\d+(\.\d+)*)$")
    bad_keys = [k for k in refs if not key_re.match(k)]
    wrong_target = [k for k, v in refs.items() if isinstance(v, list) and v[0] != src.reference_uri]
    oob = []
    total = os.path.getsize(src.local_file) if src.local_file else args.file_size
    if total:
        oob = [k for k, v in refs.items() if isinstance(v, list) and v[1] + v[2] > total]
    orphan = []
    for k in refs:
        parts = k.split("/")[:-1]
        for i in range(len(parts) + 1):
            gk = "/".join(parts[:i])
            if f"{gk + '/' if gk else ''}.zgroup" not in refs and not (
                    i == len(parts) and f"{gk + '/' if gk else ''}.zarray" in refs):
                orphan.append(k)
                break
    n_chunk_refs = sum(1 for v in refs.values() if isinstance(v, list))
    report["structure"] = {"keys": len(refs), "chunk_references": n_chunk_refs,
                           "inline_references": sum(1 for v in refs.values()
                                                    if isinstance(v, str) and v.startswith("base64:")),
                           "invalid_keys": bad_keys[:20], "invalid_key_count": len(bad_keys),
                           "references_not_targeting_source_uri": len(wrong_target),
                           "references_beyond_file_size": len(oob), "file_size_checked": total,
                           "keys_without_parent_group": orphan[:20], "orphan_key_count": len(orphan)}
    for label, lst in (("invalid Zarr key", bad_keys), ("reference not targeting source URI", wrong_target),
                       ("reference beyond end of file", oob), ("key without parent .zgroup", orphan)):
        if lst:
            fail.append({"scope": "manifest", "check": label, "detail": f"{len(lst)} keys, e.g. {lst[:3]}"})

    # Reader-side access through fsspec ReferenceFileSystem (no GDAL / HDF5 in this path).
    if src.local_file and not args.validate_remote:
        fo = {k: ([src.local_file, v[1], v[2]] if isinstance(v, list) else v) for k, v in refs.items()}
        so = {"fo": fo, "remote_protocol": "file"}
        report["decode_source"] = f"local copy {src.local_file} (references re-targeted in memory)"
    else:
        so = {"fo": manifest_path, "remote_protocol": src.remote_protocol,
              "remote_options": src.remote_options}
        report["decode_source"] = f"remote {src.reference_uri} via fsspec {src.remote_protocol}"
    # zarr-python 3 drives fsspec asynchronously; the reference FS and its async target must agree.
    aso = dict(so, asynchronous=True)
    if src.remote_protocol not in ("file", None):
        aso["remote_options"] = dict(so.get("remote_options") or {}, asynchronous=True)
    try:
        zg = zarr.open_group("reference://", mode="r", zarr_format=2, storage_options=aso,
                             use_consolidated=True)
    except Exception as e:  # noqa: BLE001
        zg = None
        fail.append({"scope": "manifest", "check": "zarr.open_group", "detail": f"{type(e).__name__}: {e}"})

    # Raw range request + numcodecs decode for one principal chunk, bypassing zarr too.
    try:
        import numcodecs
        prin = sorted((r for r in datasets.values() if r["category"] == CAT_PRINCIPAL and r["zarr_path"]
                       and r["layout"] == "chunked" and r["allocated_chunks"]), key=lambda r: r["path"])
        if prin:
            r = prin[0]
            za = json.loads(refs[f"{r['zarr_path']}/.zarray"])
            pre = r["zarr_path"] + "/"
            key = sorted(k for k in refs if k.startswith(pre) and k[len(pre):][0].isdigit())[0]
            target, off, size = refs[key]
            rfs = fsspec.filesystem("reference", **so)
            t0 = time.time()
            raw = rfs.cat_file(key)
            buf = raw
            if za["compressor"]:
                buf = numcodecs.get_codec(za["compressor"]).decode(buf)
            for fc in reversed(za["filters"] or []):
                buf = numcodecs.get_codec(fc).decode(buf)
            arr = np.frombuffer(bytes(buf), dtype=za["dtype"]).reshape(za["chunks"])
            idx = tuple(int(i) for i in key.rsplit("/", 1)[1].split("."))
            reg = _region(idx, za["chunks"], za["shape"])
            ref_vals = h5[r["path"]][reg]
            ok = np.array_equal(arr[tuple(slice(0, s.stop - s.start) for s in reg)], ref_vals, equal_nan=True)
            report["raw_decode"] = {"hdf5_path": r["path"], "key": key, "byte_offset": off, "byte_count": size,
                                    "fetched_bytes": len(raw), "seconds": round(time.time() - t0, 3),
                                    "decoded_equals_hdf5": bool(ok)}
            if not ok or len(raw) != size:
                fail.append({"scope": r["path"], "check": "raw range request decode",
                             "detail": "decoded bytes differ from HDF5 values"})
    except Exception as e:  # noqa: BLE001
        fail.append({"scope": "manifest", "check": "raw range request decode", "detail": f"{type(e).__name__}: {e}"})

    for rec in sorted(datasets.values(), key=lambda r: r["path"]):
        if rec["status"] not in (STATUS_OK, STATUS_WARN):
            if rec["rank"] == 2 or rec["dimscale_of"] or rec["gridmap_of"]:
                entry = {"hdf5_path": rec["path"], "status": rec["status"], "checks": {},
                         "note": "; ".join(i["detail"] for i in rec["issues"] if i["severity"] == "error")}
                report["arrays"].append(entry)
                if rec["rank"] == 2:
                    fail.append({"scope": rec["path"], "check": "included in manifest",
                                 "detail": entry["note"] or rec["status"]})
            continue
        zp = rec["zarr_path"]
        entry = {"hdf5_path": rec["path"], "zarr_path": zp, "status": rec["status"], "checks": {}}
        report["arrays"].append(entry)
        chk = entry["checks"]

        def record(name, ok, detail=None, _entry=entry, _chk=chk):
            _chk[name] = {"ok": bool(ok)} if detail is None else {"ok": bool(ok), "detail": detail}
            if not ok:
                fail.append({"scope": _entry["hdf5_path"], "check": name, "detail": detail})
        try:
            dset = h5[rec["path"]]
            za = json.loads(refs[f"{zp}/.zarray"])
            exp_shape = list(dset.shape)
            exp_chunks = list(dset.chunks) if dset.chunks else [max(1, s) for s in exp_shape]
            exp_comp, exp_filt = _expected_codecs(dset)
            record("manifest_shape", za["shape"] == exp_shape, f"{za['shape']} vs HDF5 {exp_shape}")
            record("manifest_dtype", np.dtype(za["dtype"]) == dset.dtype, f"{za['dtype']} vs HDF5 {dset.dtype.str}")
            record("manifest_chunks", za["chunks"] == exp_chunks, f"{za['chunks']} vs HDF5 {exp_chunks}")
            record("manifest_codecs", za["compressor"] == exp_comp and za["filters"] == exp_filt,
                   f"compressor={_codec_ids(za['compressor'])} filters={_codec_ids(za['filters'])} vs HDF5 "
                   f"pipeline {rec['pipeline']}")
            exp_fill = encode_fill_value(dset.fillvalue, dset.dtype)
            record("manifest_fill_value", za["fill_value"] == exp_fill, f"{za['fill_value']!r} vs HDF5 {exp_fill!r}")
            # Chunk keys must resolve to exactly the chunks of this HDF5 dataset.
            prefix = zp + "/"
            got = {k[len(prefix):]: v for k, v in refs.items()
                   if k.startswith(prefix) and "/" not in k[len(prefix):] and k[len(prefix)][0].isdigit()}
            if rec["layout"] == "chunked":
                chunks, _ = _enumerate_chunks(dset.id)
                want = {".".join(str(o // c) for o, c in zip(off, exp_chunks)): [src.reference_uri, bo, sz]
                        for off, _fm, bo, sz in chunks}
            elif rec["layout"] == "contiguous":
                off = dset.id.get_offset()
                key = ".".join("0" for _ in exp_shape) or "0"
                want = {key: [src.reference_uri, int(off), int(dset.id.get_storage_size())]} if off is not None else {}
            else:
                want = {k: v for k, v in got.items()}
            mism = sorted(set(want) ^ set(got)) + sorted(k for k in set(want) & set(got) if want[k] != got[k])
            record("chunk_keys_resolve", not mism,
                   f"{len(got)} refs vs {len(want)} allocated HDF5 chunks; {len(mism)} mismatched" +
                   (f" e.g. {mism[:3]}" if mism else ""))
            if zg is None:
                continue
            za_obj = zg[zp]
            record("zarr_shape", list(za_obj.shape) == exp_shape, f"{list(za_obj.shape)}")
            record("zarr_dtype", za_obj.dtype == dset.dtype, f"{za_obj.dtype}")
            record("zarr_chunks", list(za_obj.chunks) == exp_chunks, f"{list(za_obj.chunks)}")
            md = za_obj.metadata
            zc = md.compressor.get_config() if md.compressor is not None else None
            zf = [f.get_config() for f in md.filters] if md.filters else None
            record("zarr_codecs", _codec_ids(zc) == _codec_ids(exp_comp) and _codec_ids(zf) == _codec_ids(exp_filt),
                   f"compressor={_codec_ids(zc)} filters={_codec_ids(zf)}")
            picks, missing = _sample_chunk_indices(rec, exp_chunks, args.sample_chunks)
            samples = []
            t0 = time.time()
            for idx in picks + missing:
                reg = _region(idx, exp_chunks, exp_shape) if exp_shape else ()
                zv = za_obj[reg] if exp_shape else za_obj[()]
                hv = dset[reg] if exp_shape else dset[()]
                ok = np.array_equal(np.asarray(zv), np.asarray(hv), equal_nan=True)
                samples.append({"chunk": ".".join(map(str, idx)) or "0", "allocated": idx not in missing,
                                "equal_to_hdf5": bool(ok)})
                if not ok:
                    fail.append({"scope": rec["path"], "check": "decoded sample equals HDF5",
                                 "detail": f"chunk {idx}"})
            chk["decoded_samples"] = {"ok": all(s["equal_to_hdf5"] for s in samples), "samples": samples,
                                      "seconds": round(time.time() - t0, 3)}
        except Exception as e:  # noqa: BLE001
            record("exception", False, f"{type(e).__name__}: {e}")

    # xarray: open the whole tree, check each included array appears with the right dims.
    try:
        import xarray as xr
        t0 = time.time()
        tree = xr.open_datatree("reference://", engine="zarr", consolidated=True, zarr_format=2,
                                storage_options=aso, mask_and_scale=False, decode_times=False,
                                decode_timedelta=False, chunks=None)
        found = {}
        for node in tree.subtree:
            for name, var in node.ds.variables.items():
                found[(node.path.strip("/") + "/" + name).strip("/")] = (list(var.dims), list(var.shape))
        missing = []
        for rec in datasets.values():
            if rec["zarr_path"]:
                v = found.get(rec["zarr_path"])
                if v is None or v[1] != rec["shape"]:
                    missing.append(rec["zarr_path"])
        report["xarray"] = {"ok": not missing, "groups": len(list(tree.subtree)), "variables": len(found),
                            "arrays_missing_or_wrong_shape": missing, "seconds": round(time.time() - t0, 3)}
        prin = sorted(r["zarr_path"] for r in datasets.values() if r["category"] == CAT_PRINCIPAL and r["zarr_path"])
        if prin:
            grp, var = prin[0].rsplit("/", 1)
            da = tree[grp].ds[var]
            win = da.isel({da.dims[0]: slice(0, 64), da.dims[1]: slice(0, 64)}).values
            ref_vals = h5[f"{GCOV_ROOT}/{prin[0]}"][0:64, 0:64]
            report["xarray"]["window_read"] = {"variable": prin[0], "window": "[0:64, 0:64]",
                                               "equal_to_hdf5": bool(np.array_equal(win, ref_vals, equal_nan=True)),
                                               "coords": [c for c in da.coords]}
            if not report["xarray"]["window_read"]["equal_to_hdf5"]:
                fail.append({"scope": prin[0], "check": "xarray window read", "detail": "values differ"})
        if missing:
            fail.append({"scope": "xarray", "check": "open_datatree", "detail": f"missing/wrong shape: {missing[:5]}"})
    except Exception as e:  # noqa: BLE001
        report["xarray"] = {"ok": False, "error": f"{type(e).__name__}: {e}"}
        fail.append({"scope": "xarray", "check": "open_datatree", "detail": f"{type(e).__name__}: {e}"})
    report["finished_utc"] = utcnow()
    report["failure_count"] = len(fail)
    report["ok"] = not fail
    return report


def write_validation_md(path, rep):
    with open(path, "w") as f:
        f.write("# Validation report\n\n")
        f.write(f"- Result: **{'PASS' if rep['ok'] else 'FAIL'}** ({rep['failure_count']} failures)\n")
        f.write(f"- Reference URI: `{rep['reference_uri']}`\n- Decode source: {rep.get('decode_source')}\n")
        s = rep["structure"]
        f.write(f"- Manifest keys: {s.get('keys')}, chunk references: {s.get('chunk_references')}, "
                f"inline: {s.get('inline_references')}, invalid keys: {s.get('invalid_key_count')}, "
                f"refs not targeting source URI: {s.get('references_not_targeting_source_uri')}, "
                f"refs beyond EOF: {s.get('references_beyond_file_size')} (file size {s.get('file_size_checked')})\n")
        if rep.get("raw_decode"):
            r = rep["raw_decode"]
            f.write(f"- Raw range request (fsspec + numcodecs, no zarr/HDF5): `{r['hdf5_path']}` key `{r['key']}` "
                    f"bytes {r['byte_offset']}+{r['byte_count']} -> equal to HDF5: {r['decoded_equals_hdf5']}\n")
        x = rep.get("xarray", {})
        f.write(f"- xarray open_datatree: ok={x.get('ok')} groups={x.get('groups')} variables={x.get('variables')}"
                + (f" error={x.get('error')}" if x.get("error") else "") + "\n")
        if x.get("window_read"):
            f.write(f"- xarray window read `{x['window_read']['variable']}{x['window_read']['window']}` equal to "
                    f"HDF5: {x['window_read']['equal_to_hdf5']}\n")
        f.write("\n## Failures\n\n")
        if not rep["failures"]:
            f.write("None.\n")
        for fl in rep["failures"]:
            f.write(f"- `{fl['scope']}` {fl['check']}: {fl['detail']}\n")
        f.write("\n## Per-array checks\n\n| HDF5 path | status | checks passed | decoded samples |\n"
                "|---|---|---|---|\n")
        for a in rep["arrays"]:
            c = a["checks"]
            npass = sum(1 for v in c.values() if v["ok"])
            ds = c.get("decoded_samples", {})
            f.write(f"| `{a['hdf5_path']}` | {a['status']} | {npass}/{len(c)} | "
                    f"{len(ds.get('samples', []))} ok={ds.get('ok', '')} |\n")


def write_issue_report(outdir, datasets, links, attr_issues, dim_problems):
    rows = []
    for r in sorted(datasets.values(), key=lambda r: r["path"]):
        for i in r["issues"]:
            rows.append({"hdf5_path": r["path"], "rank": r["rank"], "dtype_class": r["dtype_class"],
                         "status": r["status"], **i})
        if r["zarr_path"] is None and r["status"] == STATUS_INVENTORY:
            why = {0: "scalar dataset", 1: "rank-1 dataset not used as a coordinate"}.get(
                r["rank"], f"rank-{r['rank']} dataset")
            if r["dtype_class"] in ("string", "compound", "vlen", "reference", "opaque"):
                why += f" with {r['dtype_class']} datatype"
            rows.append({"hdf5_path": r["path"], "rank": r["rank"], "dtype_class": r["dtype_class"],
                         "status": r["status"], "kind": "not-in-manifest", "severity": "info",
                         "detail": f"{why}; inventoried only (see gcov_dataset_inventory.csv)"})
    for lk in links:
        rows.append({"hdf5_path": lk["path"], "rank": None, "dtype_class": None,
                     "status": STATUS_SKIP if lk["type"] != "hard-alias" else STATUS_INVENTORY,
                     "kind": f"{lk['type']}-link", "severity": "warning" if lk["type"] == "external" else "info",
                     "detail": f"{lk['type']} link to {lk['target']}"})
    attr_rows = [{"hdf5_path": p, "detail": d} for p, lst in sorted(attr_issues.items()) for d in lst]
    global_notes = [
        "Zarr v2 has a single fill_value per array. The manifest uses the HDF5 dataset fill value (what HDF5 "
        "returns for unallocated chunks). NISAR also defines a CF '_FillValue' attribute (often NaN) that can "
        "differ; it is kept in .zattrs, but xarray's Zarr v2 backend derives _FillValue from .zarray "
        "fill_value, so open with mask_and_scale=False or apply the attribute yourself.",
        "HDF5 dimension scales (DIMENSION_LIST/REFERENCE_LIST object references, CLASS/NAME) are not "
        "representable in Zarr v2; they are replaced by _ARRAY_DIMENSIONS naming the referenced coordinate.",
        "Non-finite float attributes are stored as the JSON strings 'NaN', 'Infinity', '-Infinity'.",
        "HDF5 string/compound/enum/reference datatypes and scalar metadata datasets are inventoried but not "
        "exported as Zarr arrays; scalar identification values are copied to the root .zattrs "
        "'nisar_identification'.",
    ]
    rep = {"issues": rows, "attribute_issues": attr_rows, "dimension_conflicts": dim_problems,
           "global_notes": global_notes,
           "counts": {"skipped_unsupported": sum(1 for r in datasets.values() if r["status"] == STATUS_SKIP),
                      "included_with_warnings": sum(1 for r in datasets.values() if r["status"] == STATUS_WARN),
                      "inventory_only": sum(1 for r in datasets.values() if r["status"] == STATUS_INVENTORY),
                      "links_not_followed": sum(1 for lk in links if lk["type"] in ("soft", "external"))}}
    jp = os.path.join(outdir, "unsupported_lossy_report.json")
    write_json(jp, rep)
    mp = os.path.join(outdir, "unsupported_lossy_report.md")
    with open(mp, "w") as f:
        f.write("# Unsupported / lossy dataset report\n\n")
        f.write("Counts: " + ", ".join(f"{k}={v}" for k, v in rep["counts"].items()) + "\n\n## Global notes\n\n")
        for n in global_notes:
            f.write(f"- {n}\n")
        for sev, title in (("error", "Unsupported (skipped)"), ("warning", "Lossy / portability warnings"),
                           ("info", "Informational")):
            sel = [r for r in rows if r["severity"] == sev]
            f.write(f"\n## {title} ({len(sel)})\n\n")
            if not sel:
                f.write("None.\n")
                continue
            f.write("| HDF5 path | rank | dtype class | status | kind | detail |\n|---|---|---|---|---|---|\n")
            for r in sel:
                f.write(f"| `{r['hdf5_path']}` | {r['rank']} | {r['dtype_class']} | {r['status']} | "
                        f"{r['kind']} | {r['detail']} |\n")
        f.write(f"\n## Dimension conflicts ({len(dim_problems)})\n\n")
        for p in dim_problems or ["None."]:
            f.write(f"- {p}\n")
        f.write(f"\n## Attribute conversion notes ({len(attr_rows)})\n\n")
        for r in attr_rows:
            f.write(f"- `{r['hdf5_path']}`: {r['detail']}\n")
    return jp, mp


# --------------------------------------------------------------------------------------
# Environment, packaging, CLI
# --------------------------------------------------------------------------------------

def environment_info():
    pkgs = ["earthaccess", "h5py", "numpy", "fsspec", "aiohttp", "requests", "zarr", "numcodecs",
            "xarray", "s3fs", "python-cmr", "pandas", "gdal"]
    vers = {}
    for p in pkgs:
        try:
            vers[p] = importlib.metadata.version(p)
        except importlib.metadata.PackageNotFoundError:
            vers[p] = None
    try:
        from osgeo import gdal
        vers["gdal"] = gdal.__version__
    except ImportError:
        pass
    info = {"python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
            "packages": vers, "hdf5_library": h5py.version.hdf5_version if h5py else None,
            "tool": f"{TOOL_NAME} {TOOL_VERSION}"}
    with open(os.path.abspath(__file__), "rb") as f:
        info["script_sha256"] = hashlib.sha256(f.read()).hexdigest()
    try:
        info["git_commit"] = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True,
                                            text=True, timeout=10).stdout.strip() or None
    except Exception:  # noqa: BLE001
        info["git_commit"] = None
    return info


def artifact_readme(gid, manifest_name, sel, src):
    return f"""# NISAR GCOV virtual Zarr: {gid}

Files:
- `{manifest_name}` - consolidated Kerchunk (reference spec v1) JSON describing a Zarr v2 hierarchy rooted
  at HDF5 group `{GCOV_ROOT}`. Every chunk reference targets `{src.reference_uri}`.
- `gcov_chunk_table.csv` / `gcov_chunk_table.md` - one row per rank-2 dataset, subtotals and summary.
- `gcov_dataset_inventory.csv` - every HDF5 dataset/link under `{GCOV_ROOT}`.
- `granule_metadata.json` - CMR metadata, selection parameters, access URLs and redirect behaviour.
- `environment.json` - tool/library versions.
- `validation_report.json` / `validation_report.md` - validation results (failures are listed, not hidden).
- `unsupported_lossy_report.json` / `unsupported_lossy_report.md` - skipped and lossy datasets/attributes.

Open (requires Earthdata Login for HTTPS references):

```python
import earthaccess, xarray as xr
earthaccess.login()
fs = earthaccess.get_fsspec_https_session()
so = {{"fo": "{manifest_name}", "remote_protocol": "https",
      "remote_options": {{"client_kwargs": fs.client_kwargs, "asynchronous": True}},
      "asynchronous": True}}  # zarr-python 3 drives fsspec asynchronously
tree = xr.open_datatree("reference://", engine="zarr", consolidated=True, zarr_format=2,
                        storage_options=so, mask_and_scale=False, chunks=None)
print(tree["grids/frequencyA"])
```

The HTTPS archive URL redirects to a short-lived signed CloudFront URL on every request; only the
stable archive URL is stored. If the archive URL ever changes, rewrite the targets with
`python nisar_gcov_virtual_zarr.py retarget {manifest_name} <new-uri>`.
"""


def make_zip(zip_path, files):
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for p in files:
            z.write(p, arcname=os.path.basename(p))
    return zip_path


def cmd_generate(args):
    if h5py is None:
        raise SystemExit("h5py is required")
    t_start = time.time()
    sel = select_granule(args)
    g = sel.get("granule")
    if g:
        log(f"selected {g['granule_ur']} ({g['granule_concept_id']}) via {sel['selection']['method']} "
            f"from {sel['selection'].get('candidate_count')} candidates")
    if args.select_only:
        print(dumps({k: v for k, v in sel.items() if not k.startswith("_")}, indent=2))
        return 0
    access_probe = None
    if g and g["https_links"] and not (args.local_file and not args.probe_url):
        access_probe = probe_https_access(g["https_links"][0])
        log(f"access: {access_probe.get('access_type')}")
    src = open_source(sel, args)
    try:
        ident = read_identification(src.h5)
        gid = safe_id(args.granule_id or (g["granule_ur"] if g else None) or ident.get("granuleId")
                      or os.path.splitext(os.path.basename(src.reference_uri))[0])
        outdir = os.path.join(os.path.abspath(args.output_dir), gid)
        os.makedirs(outdir, exist_ok=True)
        log(f"inventory of {GCOV_ROOT} ...")
        t0 = time.time()
        datasets, links, groups = inventory(src.h5)
        plan_manifest(datasets, include_rank3=not args.no_rank3)
        log(f"{len(datasets)} datasets, {len(links)} non-primary links in {time.time() - t0:.1f}s")
        refs, per_bytes, attr_issues = build_manifest(src, datasets, groups, sel)
        dim_problems = check_dim_consistency(datasets, refs)
        manifest_name = f"nisar_gcov_{gid}.kerchunk.json"
        manifest_path = os.path.join(outdir, manifest_name)
        with open(manifest_path, "w") as f:
            f.write(dumps({"version": 1, "refs": refs}, separators=(",", ":")))
        manifest_bytes = os.path.getsize(manifest_path)
        sel_summary = (f"Granule `{gid}`; reference URI `{src.reference_uri}`; manifest `{manifest_name}` "
                       f"({manifest_bytes} bytes).")
        tables = chunk_table(datasets, per_bytes, manifest_bytes)
        csv_path, md_path = write_chunk_table(outdir, tables, sel_summary)
        inv_path = write_inventory(outdir, datasets, links)
        meta = {
            "generated_utc": utcnow(),
            "granule_id": gid,
            "reference_uri": src.reference_uri,
            "reference_uri_note": "Stable archive URI used as the target of every chunk reference. It is not a "
                                  "signed URL; HTTPS access needs Earthdata Login and is redirected per request.",
            "access_mode": src.access,
            "auth_strategy": sel.get("auth_strategy"),
            "selection": sel["selection"],
            "cmr": g,
            "access_probe": access_probe,
            "hdf5_identification": ident,
            "hdf5_root_attributes": group_attrs(dict(src.h5.attrs.items()), "/")[0],
        }
        meta_path = os.path.join(outdir, "granule_metadata.json")
        write_json(meta_path, meta)
        env_path = os.path.join(outdir, "environment.json")
        write_json(env_path, environment_info())
        if access_probe and access_probe.get("content_length") and not args.file_size:
            args.file_size = access_probe["content_length"]
        elif g and g.get("size_bytes") and not args.file_size:
            args.file_size = g["size_bytes"]
        log("validating manifest ...")
        if args.skip_validation:
            rep = {"ok": False, "skipped": True, "failure_count": 1, "structure": {}, "arrays": [],
                   "reference_uri": src.reference_uri,
                   "failures": [{"scope": "all", "check": "validation", "detail": "skipped by --skip-validation"}]}
        else:
            rep = validate(src, refs, datasets, manifest_path, args)
        val_json = os.path.join(outdir, "validation_report.json")
        write_json(val_json, rep)
        val_md = os.path.join(outdir, "validation_report.md")
        write_validation_md(val_md, rep)
        iss_json, iss_md = write_issue_report(outdir, datasets, links, attr_issues, dim_problems)
        readme = os.path.join(outdir, "README.md")
        with open(readme, "w") as f:
            f.write(artifact_readme(gid, manifest_name, sel, src))
        files = [manifest_path, csv_path, md_path, inv_path, meta_path, env_path, val_json, val_md,
                 iss_json, iss_md, readme]
        zip_path = make_zip(os.path.join(os.path.abspath(args.output_dir),
                                         f"nisar_gcov_virtual_zarr_{gid}.zip"), files)
    finally:
        src.close()
    two_d, _, total, support, _, _ = tables
    print(f"\n2-D datasets: {len(two_d)} ({total['manifest_status']}), supporting arrays: {len(support)}, "
          f"allocated chunks: {total['allocated_chunks']}, skipped chunks: {total['skipped_chunks']}")
    print(f"validation: {'PASS' if rep['ok'] else 'FAIL'} ({rep['failure_count']} failures); "
          f"elapsed {time.time() - t_start:.1f}s\n")
    print("Artifacts:")
    for p in files + [zip_path]:
        print(f"  {os.path.getsize(p):>12,d}  {os.path.abspath(p)}")
    return 0 if rep["ok"] else 3


def cmd_retarget(args):
    with open(args.manifest) as f:
        doc = json.load(f)
    refs = doc["refs"]
    old = args.old_uri
    if old is None:
        old = json.loads(refs[".zattrs"])["nisar_virtual_zarr"]["source_uri"]
    n = 0
    for k, v in refs.items():
        if isinstance(v, list) and v[0] == old:
            v[0] = args.new_uri
            n += 1
    for key in (".zattrs",):
        za = json.loads(refs[key])
        za.setdefault("nisar_virtual_zarr", {})["source_uri"] = args.new_uri
        refs[key] = dumps(za)
    zm = json.loads(refs[".zmetadata"])
    zm["metadata"][".zattrs"] = json.loads(refs[".zattrs"])
    refs[".zmetadata"] = dumps(zm)
    out = args.output or args.manifest
    with open(out, "w") as f:
        f.write(dumps(doc, separators=(",", ":")))
    print(f"retargeted {n} references from {old} to {args.new_uri} -> {os.path.abspath(out)}")
    return 0


def build_parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("generate", help="select a granule and build the virtual Zarr artifact")
    sg = g.add_argument_group("granule selection")
    sg.add_argument("--short-name", default=DEFAULT_SHORT_NAME, help=f"CMR collection (default {DEFAULT_SHORT_NAME})")
    sg.add_argument("--seed", type=int, default=DEFAULT_SEED, help=f"selection seed (default {DEFAULT_SEED})")
    sg.add_argument("--concept-id", help="explicit granule concept ID (G...-ASF)")
    sg.add_argument("--granule-ur", help="explicit granule UR (with or without .h5)")
    sg.add_argument("--temporal", nargs=2, metavar=("START", "END"), default=list(DEFAULT_TEMPORAL),
                    help="candidate temporal range (default %(default)s)")
    sg.add_argument("--all-time", action="store_true", help="do not restrict candidates by time")
    sg.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"), help="candidate bounding box")
    sg.add_argument("--max-candidates", type=int, default=DEFAULT_MAX_CANDIDATES)
    sg.add_argument("--select-only", action="store_true", help="print the selection and exit")
    ag = g.add_argument_group("access")
    ag.add_argument("--auth", choices=["auto", "environment", "netrc", "interactive"], default="auto")
    ag.add_argument("--access", choices=["external", "direct"], default="external",
                    help="external = HTTPS archive URL (default); direct = s3:// (AWS us-west-2 only)")
    ag.add_argument("--local-file", help="read a local HDF5 copy instead of the remote object")
    ag.add_argument("--remote-uri", help="reference target to record when using --local-file")
    ag.add_argument("--granule-id", help="override the granule ID used in output names")
    ag.add_argument("--block-size", type=int, default=8 * 2**20, help="fsspec block size for metadata reads")
    ag.add_argument("--probe-url", action="store_true", help="probe HTTPS redirects even with --local-file")
    og = g.add_argument_group("output / validation")
    og.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="default: conda-build/tests/output")
    og.add_argument("--no-rank3", action="store_true", help="leave rank-3 metadata cubes out of the manifest")
    og.add_argument("--sample-chunks", type=int, default=3, help="allocated chunks decoded per array")
    og.add_argument("--validate-remote", action="store_true",
                    help="with --local-file, decode samples from the remote URI instead of the local copy")
    og.add_argument("--file-size", type=int, help="object size for out-of-range reference checks")
    og.add_argument("--skip-validation", action="store_true")
    g.set_defaults(func=cmd_generate)
    r = sub.add_parser("retarget", help="rewrite the reference target URI of an existing manifest")
    r.add_argument("manifest")
    r.add_argument("new_uri")
    r.add_argument("--old-uri", help="default: nisar_virtual_zarr.source_uri from the root .zattrs")
    r.add_argument("-o", "--output", help="write to a new file instead of in place")
    r.set_defaults(func=cmd_retarget)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
