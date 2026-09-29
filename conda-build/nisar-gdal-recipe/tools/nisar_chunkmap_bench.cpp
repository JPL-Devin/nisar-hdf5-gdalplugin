// nisar_chunkmap_bench.cpp
/**************************************************************************************************************************/
/* Copyright 2025, by the California Institute of Technology.                                                             */
/* ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.                                                */
/* Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.*/
/*                                                                                                                        */
/* This software may be subject to U.S. export control laws.                                                              */
/* By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations.            */
/* User has the responsibility to obtain export licenses, or other export authority as may be required                    */
/* before exporting such information to foreign countries or providing access to foreign persons.                         */
/**************************************************************************************************************************/

// Validation + benchmark harness for the native HDF5 chunk-map parser.
//
//   nisar_chunkmap_bench [options] FILE [FILE...]
//
//   --validate           Compare native vs H5Dchunk_iter for every chunked 2-D/3-D dataset (default on).
//   --no-validate        Skip validation.
//   --bench PATH         Benchmark chunk-map discovery for this HDF5 dataset path (repeatable).
//   --bench-auto         Benchmark the chunked dataset with the most chunks in each file.
//   --repeat N           Benchmark repetitions (median wall time reported), default 3.
//   --block-sizes A,B    Native metadata read granularity(ies) in bytes, default 0 (exact reads).
//   --quiet              Only print divergences and summary tables.
//
// FILE may be a local path or any GDAL VSI path (/vsicurl/https://..., /vsis3/...).
// libhdf5 always goes through the plugin's VSI VFL with the same 4 MiB page buffer the driver uses,
// so both paths are measured at the same (VSI) layer.

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "cpl_conv.h"
#include "cpl_json.h"
#include "cpl_string.h"
#include "cpl_vsi.h"
#include "hdf5.h"

#include "../hdf5native.h"
#include "../hdf5vfl.h"

using NisarHDF5Native::ChunkRecord;
using NisarHDF5Native::DatasetInfo;
using NisarHDF5Native::Status;

namespace {

struct NetStats {
    uint64_t nGet = 0;
    uint64_t nHead = 0;
    uint64_t nBytes = 0;
};

NetStats GetNetStats()
{
    NetStats o;
    char *pszJSON = VSINetworkStatsGetAsSerializedJSON(nullptr);
    if (!pszJSON) return o;
    CPLJSONDocument oDoc;
    if (oDoc.LoadMemory(pszJSON)) {
        CPLJSONObject oRoot = oDoc.GetRoot();
        o.nGet = static_cast<uint64_t>(oRoot.GetLong("methods/GET/count", 0));
        o.nHead = static_cast<uint64_t>(oRoot.GetLong("methods/HEAD/count", 0));
        o.nBytes = static_cast<uint64_t>(oRoot.GetLong("methods/GET/downloaded_bytes", 0));
    }
    CPLFree(pszJSON);
    return o;
}

struct Measure {
    double dfMs = 0;
    uint64_t nReads = 0;     // VSIFReadL calls issued by the path under test
    uint64_t nReadBytes = 0; // bytes requested by those calls
    NetStats oNet;           // HTTP requests actually sent (remote only)
    size_t nChunks = 0;
    bool bOk = false;
};

double NowMs()
{
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

hid_t OpenHDF5(const std::string &osPath)
{
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_driver(fapl, NisarVFL::HDF5VFLGetFileDriver(), nullptr);
    H5Pset_page_buffer_size(fapl, 4 * 1024 * 1024, 0, 0);
    hid_t hFile = H5Fopen(osPath.c_str(), H5F_ACC_RDONLY, fapl);
    H5Pclose(fapl);
    return hFile;
}

struct IterCtx {
    int nRank;
    std::vector<ChunkRecord> *paoRecs;
};

int ChunkIterCB(const hsize_t *offset, unsigned filter_mask, haddr_t addr, hsize_t size, void *op_data)
{
    IterCtx *p = static_cast<IterCtx *>(op_data);
    ChunkRecord o;
    o.anOffset.assign(offset, offset + p->nRank);
    o.nFilterMask = filter_mask;
    o.nAddr = addr;
    o.nSize = size;
    p->paoRecs->push_back(std::move(o));
    return 0;
}

bool LibChunks(hid_t hDset, int nRank, std::vector<ChunkRecord> &aoRecs)
{
    aoRecs.clear();
    IterCtx ctx{nRank, &aoRecs};
    return H5Dchunk_iter(hDset, H5P_DEFAULT, ChunkIterCB, &ctx) >= 0;
}

// Same mapping as NisarRasterBand::MapChunks() for one band.
struct MapEntry {
    int nBlockX, nBlockY;
    uint64_t nOffset;
    uint64_t nLength;
    bool bIsMissing;
    bool operator==(const MapEntry &o) const
    {
        return nBlockX == o.nBlockX && nBlockY == o.nBlockY && nOffset == o.nOffset && nLength == o.nLength &&
               bIsMissing == o.bIsMissing;
    }
};

std::vector<MapEntry> BuildMap(const std::vector<ChunkRecord> &aoRecs, int nRank, int nBand, uint64_t nXSize,
                               uint64_t nYSize, int nBlockX, int nBlockY)
{
    const int nPerRow = static_cast<int>((nXSize + nBlockX - 1) / nBlockX);
    const int nPerCol = static_cast<int>((nYSize + nBlockY - 1) / nBlockY);
    std::vector<MapEntry> ao(static_cast<size_t>(nPerRow) * nPerCol);
    for (int y = 0; y < nPerCol; y++)
        for (int x = 0; x < nPerRow; x++) ao[y * nPerRow + x] = {x, y, 0, 0, true};
    for (const ChunkRecord &r : aoRecs) {
        int bx, by;
        if (nRank == 3) {
            if (r.anOffset[0] != static_cast<uint64_t>(nBand - 1)) continue;
            by = static_cast<int>(r.anOffset[1] / nBlockY);
            bx = static_cast<int>(r.anOffset[2] / nBlockX);
        }
        else {
            by = static_cast<int>(r.anOffset[0] / nBlockY);
            bx = static_cast<int>(r.anOffset[1] / nBlockX);
        }
        const int idx = by * nPerRow + bx;
        if (idx >= 0 && idx < static_cast<int>(ao.size())) {
            ao[idx].nOffset = r.nAddr;
            ao[idx].nLength = r.nSize;
            ao[idx].bIsMissing = false;
        }
    }
    return ao;
}

struct VisitCtx {
    std::vector<std::string> aosPaths;
};

herr_t VisitCB(hid_t, const char *name, const H5O_info2_t *info, void *op_data)
{
    if (info->type == H5O_TYPE_DATASET) static_cast<VisitCtx *>(op_data)->aosPaths.push_back(std::string("/") + name);
    return 0;
}

struct DsetResult {
    std::string osPath;
    std::string osShape;
    std::string osChunk;
    std::string osIndex;
    size_t nChunks = 0;
    std::string osVerdict;  // IDENTICAL / DIVERGENT / NATIVE_FALLBACK
    std::string osDetail;
};

std::string Join(const std::vector<uint64_t> &an)
{
    std::string s;
    for (size_t i = 0; i < an.size(); i++) s += (i ? "x" : "") + std::to_string(an[i]);
    return s;
}

struct FileSummary {
    std::string osFile;
    int nSB = -1;
    size_t nDatasets = 0, nChecked = 0, nIdentical = 0, nDivergent = 0, nFallback = 0, nSkipped = 0;
    std::vector<std::string> aosIndexTypes;
};

bool Validate(const std::string &osFile, bool bQuiet, FileSummary &oSum, std::vector<std::string> &aosLargest)
{
    oSum.osFile = osFile;
    hid_t hFile = OpenHDF5(osFile);
    if (hFile < 0) {
        fprintf(stderr, "ERROR: libhdf5 cannot open %s\n", osFile.c_str());
        return false;
    }
    std::unique_ptr<NisarHDF5Native::File> poNative;
    Status eOpen = NisarHDF5Native::File::Open(osFile, poNative, 0);
    oSum.nSB = poNative ? poNative->GetSuperblockVersion() : -1;
    VisitCtx vctx;
    H5Ovisit3(hFile, H5_INDEX_NAME, H5_ITER_NATIVE, VisitCB, &vctx, H5O_INFO_BASIC);
    oSum.nDatasets = vctx.aosPaths.size();

    std::vector<DsetResult> aoResults;
    size_t nMostChunks = 0;
    std::string osMost;
    for (const std::string &osPath : vctx.aosPaths) {
        hid_t hDset = H5Dopen2(hFile, osPath.c_str(), H5P_DEFAULT);
        if (hDset < 0) continue;
        hid_t hSpace = H5Dget_space(hDset);
        const int nRank = H5Sget_simple_extent_ndims(hSpace);
        std::vector<hsize_t> anDims(std::max(nRank, 0));
        if (nRank > 0) H5Sget_simple_extent_dims(hSpace, anDims.data(), nullptr);
        H5Sclose(hSpace);
        hid_t hDcpl = H5Dget_create_plist(hDset);
        const H5D_layout_t eLayout = H5Pget_layout(hDcpl);
        if ((nRank != 2 && nRank != 3) || eLayout != H5D_CHUNKED) {
            // Not a MapChunks() target; confirm the native parser rejects it cleanly.
            if (eOpen == Status::OK && (nRank == 2 || nRank == 3)) {
                DatasetInfo oInfo;
                Status e = poNative->OpenDataset(osPath, oInfo);
                if (e != Status::ERR_NOT_CHUNKED) {
                    DsetResult r;
                    r.osPath = osPath;
                    r.osVerdict = "DIVERGENT";
                    r.osDetail = std::string("non-chunked dataset: native returned ") +
                                 NisarHDF5Native::StatusToString(e) + " instead of ERR_NOT_CHUNKED";
                    aoResults.push_back(r);
                    oSum.nDivergent++;
                }
            }
            oSum.nSkipped++;
            H5Pclose(hDcpl);
            H5Dclose(hDset);
            continue;
        }
        std::vector<hsize_t> anChunk(nRank);
        H5Pget_chunk(hDcpl, nRank, anChunk.data());
        std::vector<uint16_t> anLibFilters;
        for (int i = 0; i < H5Pget_nfilters(hDcpl); i++) {
            unsigned flags;
            size_t n = 0;
            anLibFilters.push_back(static_cast<uint16_t>(H5Pget_filter2(hDcpl, i, &flags, &n, nullptr, 0, nullptr, nullptr)));
        }
        H5Pclose(hDcpl);
        oSum.nChecked++;

        DsetResult r;
        r.osPath = osPath;
        r.osShape = Join(std::vector<uint64_t>(anDims.begin(), anDims.end()));
        r.osChunk = Join(std::vector<uint64_t>(anChunk.begin(), anChunk.end()));

        std::vector<ChunkRecord> aoLib, aoNat;
        LibChunks(hDset, nRank, aoLib);
        H5Dclose(hDset);
        r.nChunks = aoLib.size();
        if (aoLib.size() > nMostChunks) { nMostChunks = aoLib.size(); osMost = osPath; }

        DatasetInfo oInfo;
        Status e = eOpen;
        if (e == Status::OK) e = poNative->OpenDataset(osPath, oInfo);
        if (e == Status::OK) e = poNative->GetChunks(oInfo, aoNat);
        r.osIndex = NisarHDF5Native::ChunkIndexToString(oInfo.eIndex);
        if (std::find(oSum.aosIndexTypes.begin(), oSum.aosIndexTypes.end(), r.osIndex) == oSum.aosIndexTypes.end())
            oSum.aosIndexTypes.push_back(r.osIndex);
        if (e != Status::OK) {
            r.osVerdict = "NATIVE_FALLBACK";
            r.osDetail = std::string(NisarHDF5Native::StatusToString(e)) + ": " +
                         (poNative ? poNative->GetLastError() : std::string("open failed"));
            oSum.nFallback++;
            aoResults.push_back(r);
            continue;
        }

        std::vector<std::string> aosDiff;
        if (oInfo.anDims != std::vector<uint64_t>(anDims.begin(), anDims.end())) aosDiff.push_back("dims");
        if (oInfo.anChunkDims != std::vector<uint64_t>(anChunk.begin(), anChunk.end())) aosDiff.push_back("chunk dims");
        if (oInfo.anFilterIds != anLibFilters) aosDiff.push_back("filter pipeline");

        // Raw H5Dchunk_iter records: same set and same iteration order.
        bool bSameOrder = aoLib.size() == aoNat.size();
        for (size_t i = 0; bSameOrder && i < aoLib.size(); i++)
            bSameOrder = aoLib[i].anOffset == aoNat[i].anOffset && aoLib[i].nAddr == aoNat[i].nAddr &&
                         aoLib[i].nSize == aoNat[i].nSize && aoLib[i].nFilterMask == aoNat[i].nFilterMask;
        if (!bSameOrder) {
            auto Less = [](const ChunkRecord &a, const ChunkRecord &b) { return a.anOffset < b.anOffset; };
            std::vector<ChunkRecord> a = aoLib, b = aoNat;
            std::sort(a.begin(), a.end(), Less);
            std::sort(b.begin(), b.end(), Less);
            size_t nMismatch = a.size() == b.size() ? 0 : std::max(a.size(), b.size());
            for (size_t i = 0; a.size() == b.size() && i < a.size(); i++)
                if (!(a[i].anOffset == b[i].anOffset && a[i].nAddr == b[i].nAddr && a[i].nSize == b[i].nSize &&
                      a[i].nFilterMask == b[i].nFilterMask))
                    nMismatch++;
            aosDiff.push_back(nMismatch ? "chunk records (" + std::to_string(nMismatch) + " of " +
                                              std::to_string(a.size()) + "/" + std::to_string(b.size()) + ")"
                                        : "chunk record order");
        }

        // MapChunks() map for every band, compared entry by entry.
        const int nBands = nRank == 3 ? static_cast<int>(anDims[0]) : 1;
        const uint64_t nY = anDims[nRank - 2], nX = anDims[nRank - 1];
        const int nBY = static_cast<int>(anChunk[nRank - 2]), nBX = static_cast<int>(anChunk[nRank - 1]);
        size_t nMapMismatch = 0, nMapEntries = 0;
        for (int b = 1; b <= nBands; b++) {
            const auto m1 = BuildMap(aoLib, nRank, b, nX, nY, nBX, nBY);
            const auto m2 = BuildMap(aoNat, nRank, b, nX, nY, nBX, nBY);
            nMapEntries += m1.size();
            for (size_t i = 0; i < m1.size(); i++)
                if (!(m1[i] == m2[i])) nMapMismatch++;
        }
        if (nMapMismatch)
            aosDiff.push_back("MapChunks map (" + std::to_string(nMapMismatch) + "/" + std::to_string(nMapEntries) + " entries)");

        if (aosDiff.empty()) {
            r.osVerdict = "IDENTICAL";
            r.osDetail = std::to_string(nBands) + " band map(s), " + std::to_string(nMapEntries) + " entries";
            oSum.nIdentical++;
        }
        else {
            r.osVerdict = "DIVERGENT";
            for (size_t i = 0; i < aosDiff.size(); i++) r.osDetail += (i ? "; " : "") + aosDiff[i];
            oSum.nDivergent++;
        }
        aoResults.push_back(r);
    }
    H5Fclose(hFile);
    if (!osMost.empty()) aosLargest.push_back(osMost);

    printf("\n## Validation: %s\n\n", CPLGetFilename(osFile.c_str()));
    printf("| dataset | shape | chunk | index | chunks | verdict | detail |\n");
    printf("|---|---|---|---|---:|---|---|\n");
    for (const DsetResult &r : aoResults) {
        if (bQuiet && r.osVerdict == "IDENTICAL") continue;
        printf("| `%s` | %s | %s | %s | %zu | %s | %s |\n", r.osPath.c_str(), r.osShape.c_str(), r.osChunk.c_str(),
               r.osIndex.c_str(), r.nChunks, r.osVerdict.c_str(), r.osDetail.c_str());
    }
    fflush(stdout);
    return true;
}

void ResetCaches()
{
    VSICurlClearCache();
    VSINetworkStatsReset();
}

// libhdf5, cold: H5Fopen + H5Dopen + H5Dchunk_iter. pIter receives the H5Dchunk_iter-only portion.
Measure BenchLib(const std::string &osFile, const std::string &osDset, Measure *pIter)
{
    Measure m;
    ResetCaches();
    uint64_t c0, b0, c1, b1, c2, b2;
    NisarVFL::HDF5VFLGetReadStats(&c0, &b0);
    const double t0 = NowMs();
    hid_t hFile = OpenHDF5(osFile);
    if (hFile < 0) return m;
    hid_t hDset = H5Dopen2(hFile, osDset.c_str(), H5P_DEFAULT);
    hid_t hSpace = H5Dget_space(hDset);
    const int nRank = H5Sget_simple_extent_ndims(hSpace);
    H5Sclose(hSpace);
    const double t1 = NowMs();
    NisarVFL::HDF5VFLGetReadStats(&c1, &b1);
    const NetStats n1 = GetNetStats();
    std::vector<ChunkRecord> aoRecs;
    m.bOk = LibChunks(hDset, nRank, aoRecs);
    const double t2 = NowMs();
    NisarVFL::HDF5VFLGetReadStats(&c2, &b2);
    const NetStats n2 = GetNetStats();
    H5Dclose(hDset);
    H5Fclose(hFile);
    m.dfMs = t2 - t0;
    m.nReads = c2 - c0;
    m.nReadBytes = b2 - b0;
    m.oNet = n2;
    m.nChunks = aoRecs.size();
    if (pIter) {
        pIter->dfMs = t2 - t1;
        pIter->nReads = c2 - c1;
        pIter->nReadBytes = b2 - b1;
        pIter->oNet.nGet = n2.nGet - n1.nGet;
        pIter->oNet.nHead = n2.nHead - n1.nHead;
        pIter->oNet.nBytes = n2.nBytes - n1.nBytes;
        pIter->nChunks = aoRecs.size();
        pIter->bOk = m.bOk;
    }
    return m;
}

// Native, cold (or warm when bResetCaches is false).
Measure BenchNative(const std::string &osFile, const std::string &osDset, size_t nBlockSize, bool bResetCaches)
{
    Measure m;
    if (bResetCaches) ResetCaches();
    const NetStats n0 = GetNetStats();
    const double t0 = NowMs();
    std::unique_ptr<NisarHDF5Native::File> poFile;
    Status e = NisarHDF5Native::File::Open(osFile, poFile, nBlockSize);
    DatasetInfo oInfo;
    std::vector<ChunkRecord> aoRecs;
    if (e == Status::OK) e = poFile->OpenDataset(osDset, oInfo);
    if (e == Status::OK) e = poFile->GetChunks(oInfo, aoRecs);
    m.dfMs = NowMs() - t0;
    const NetStats n1 = GetNetStats();
    m.bOk = e == Status::OK;
    if (poFile) {
        m.nReads = poFile->GetStats().nReadCalls;
        m.nReadBytes = poFile->GetStats().nBytesRead;
    }
    m.oNet.nGet = n1.nGet - n0.nGet;
    m.oNet.nHead = n1.nHead - n0.nHead;
    m.oNet.nBytes = n1.nBytes - n0.nBytes;
    m.nChunks = aoRecs.size();
    if (!m.bOk)
        fprintf(stderr, "native failed on %s: %s %s\n", osDset.c_str(), NisarHDF5Native::StatusToString(e),
                poFile ? poFile->GetLastError().c_str() : "");
    return m;
}

Measure Median(std::vector<Measure> v)
{
    std::sort(v.begin(), v.end(), [](const Measure &a, const Measure &b) { return a.dfMs < b.dfMs; });
    return v[v.size() / 2];
}

void PrintRow(const char *pszPath, const Measure &m)
{
    printf("| %s | %s | %.1f | %" PRIu64 " | %" PRIu64 " | %" PRIu64 " | %" PRIu64 " | %" PRIu64 " | %zu |\n", pszPath,
           m.bOk ? "ok" : "FAIL", m.dfMs, m.nReads, m.nReadBytes, m.oNet.nGet, m.oNet.nHead, m.oNet.nBytes, m.nChunks);
}

void Bench(const std::string &osFile, const std::string &osDset, int nRepeat, const std::vector<size_t> &anBlockSizes)
{
    printf("\n## Benchmark: %s `%s` (median of %d)\n\n", CPLGetFilename(osFile.c_str()), osDset.c_str(), nRepeat);
    printf("| path | status | wall ms | VSI reads | VSI bytes | HTTP GET | HTTP HEAD | HTTP bytes | chunks |\n");
    printf("|---|---|---:|---:|---:|---:|---:|---:|---:|\n");

    std::vector<Measure> aoLib, aoIter, aoWarm;
    for (int i = 0; i < nRepeat; i++) {
        Measure oIter;
        aoLib.push_back(BenchLib(osFile, osDset, &oIter));
        aoIter.push_back(oIter);
    }
    PrintRow("libhdf5 cold: H5Fopen+H5Dopen+H5Dchunk_iter", Median(aoLib));
    PrintRow("libhdf5 H5Dchunk_iter only (file+dataset already open)", Median(aoIter));
    for (size_t nBS : anBlockSizes) {
        std::vector<Measure> ao;
        for (int i = 0; i < nRepeat; i++) ao.push_back(BenchNative(osFile, osDset, nBS, true));
        const std::string osLabel = "native cold, block=" + std::to_string(nBS);
        PrintRow(osLabel.c_str(), Median(ao));
    }
    // Warm: what MapChunks() pays when the driver has already opened the file through libhdf5
    // (shared /vsicurl/ block cache is populated by H5Fopen/H5Dopen).
    for (size_t nBS : anBlockSizes) {
        std::vector<Measure> ao;
        for (int i = 0; i < nRepeat; i++) {
            ResetCaches();
            hid_t hFile = OpenHDF5(osFile);
            hid_t hDset = H5Dopen2(hFile, osDset.c_str(), H5P_DEFAULT);
            ao.push_back(BenchNative(osFile, osDset, nBS, false));
            H5Dclose(hDset);
            H5Fclose(hFile);
        }
        const std::string osLabel = "native in-driver (after libhdf5 H5Dopen), block=" + std::to_string(nBS);
        PrintRow(osLabel.c_str(), Median(ao));
    }
    fflush(stdout);
}

}  // namespace

int main(int argc, char **argv)
{
    bool bValidate = true, bQuiet = false, bBenchAuto = false;
    int nRepeat = 3;
    std::vector<size_t> anBlockSizes = {0};
    std::vector<std::string> aosFiles, aosBench;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        if (a == "--validate") bValidate = true;
        else if (a == "--no-validate") bValidate = false;
        else if (a == "--quiet") bQuiet = true;
        else if (a == "--bench-auto") bBenchAuto = true;
        else if (a == "--bench" && i + 1 < argc) aosBench.push_back(argv[++i]);
        else if (a == "--repeat" && i + 1 < argc) nRepeat = std::max(1, atoi(argv[++i]));
        else if (a == "--block-sizes" && i + 1 < argc) {
            anBlockSizes.clear();
            char **papsz = CSLTokenizeString2(argv[++i], ",", 0);
            for (char **p = papsz; p && *p; p++) anBlockSizes.push_back(static_cast<size_t>(strtoull(*p, nullptr, 10)));
            CSLDestroy(papsz);
        }
        else if (!a.empty() && a[0] == '-') { fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
        else aosFiles.push_back(a);
    }
    if (aosFiles.empty()) {
        fprintf(stderr, "usage: nisar_chunkmap_bench [--no-validate] [--quiet] [--bench PATH] [--bench-auto] "
                        "[--repeat N] [--block-sizes A,B] FILE...\n");
        return 2;
    }
    CPLSetConfigOption("CPL_VSIL_NETWORK_STATS_ENABLED", "YES");
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

    int nRet = 0;
    std::vector<FileSummary> aoSums;
    for (const std::string &osFile : aosFiles) {
        std::vector<std::string> aosLargest;
        if (bValidate) {
            FileSummary oSum;
            if (!Validate(osFile, bQuiet, oSum, aosLargest)) { nRet = 1; continue; }
            if (oSum.nDivergent) nRet = 1;
            aoSums.push_back(oSum);
        }
        else if (bBenchAuto) {
            FileSummary oSum;
            Validate(osFile, true, oSum, aosLargest);
        }
        std::vector<std::string> aosTargets = aosBench;
        if (bBenchAuto) aosTargets.insert(aosTargets.end(), aosLargest.begin(), aosLargest.end());
        for (const std::string &osDset : aosTargets) Bench(osFile, osDset, nRepeat, anBlockSizes);
    }

    if (!aoSums.empty()) {
        printf("\n## Validation summary\n\n");
        printf("| file | superblock | datasets | chunked 2-D/3-D | identical | divergent | native fallback | index types |\n");
        printf("|---|---:|---:|---:|---:|---:|---:|---|\n");
        for (const FileSummary &s : aoSums) {
            std::string osIdx;
            for (size_t i = 0; i < s.aosIndexTypes.size(); i++) osIdx += (i ? "," : "") + s.aosIndexTypes[i];
            printf("| %s | v%d | %zu | %zu | %zu | %zu | %zu | %s |\n", CPLGetFilename(s.osFile.c_str()), s.nSB,
                   s.nDatasets, s.nChecked, s.nIdentical, s.nDivergent, s.nFallback, osIdx.c_str());
        }
    }
    NisarVFL::HDF5VFLUnloadFileDriver();
    return nRet;
}
