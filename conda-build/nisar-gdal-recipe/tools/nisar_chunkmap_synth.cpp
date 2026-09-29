// nisar_chunkmap_synth.cpp
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

// Writes small synthetic HDF5 files that exercise native-parser code paths not present in
// current NISAR granules (superblock v0/v3, object header v2, dense groups, layout v4 indexes).
//
//   nisar_chunkmap_synth OUTDIR

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "hdf5.h"

namespace {

void Check(int64_t e, const char *pszWhat)
{
    if (e < 0) { fprintf(stderr, "HDF5 call failed: %s\n", pszWhat); exit(1); }
}

// Creates a chunked float32 dataset and writes every chunk whose linear index is not a multiple of nSkip.
void MakeChunked(hid_t hLoc, const char *pszName, const std::vector<hsize_t> &anDims,
                 const std::vector<hsize_t> &anMax, const std::vector<hsize_t> &anChunk, bool bDeflate,
                 int nSkip, H5D_alloc_time_t eAlloc = H5D_ALLOC_TIME_DEFAULT)
{
    const int nRank = static_cast<int>(anDims.size());
    hid_t hSpace = H5Screate_simple(nRank, anDims.data(), anMax.data());
    hid_t hDcpl = H5Pcreate(H5P_DATASET_CREATE);
    Check(H5Pset_chunk(hDcpl, nRank, anChunk.data()), "set_chunk");
    if (bDeflate) { H5Pset_shuffle(hDcpl); H5Pset_deflate(hDcpl, 1); }
    if (eAlloc != H5D_ALLOC_TIME_DEFAULT) H5Pset_alloc_time(hDcpl, eAlloc);
    hid_t hDset = H5Dcreate2(hLoc, pszName, H5T_NATIVE_FLOAT, hSpace, H5P_DEFAULT, hDcpl, H5P_DEFAULT);
    Check(hDset, pszName);

    std::vector<hsize_t> anGrid(nRank);
    size_t nChunkElems = 1, nChunks = 1;
    for (int i = 0; i < nRank; i++) {
        anGrid[i] = (anDims[i] + anChunk[i] - 1) / anChunk[i];
        nChunkElems *= anChunk[i];
        nChunks *= anGrid[i];
    }
    std::vector<float> afBuf(nChunkElems);
    hid_t hMem = H5Screate_simple(nRank, anChunk.data(), nullptr);
    std::vector<hsize_t> anScaled(nRank, 0), anStart(nRank), anCount(nRank);
    for (size_t n = 0; n < nChunks; n++) {
        if (nSkip <= 0 || n % nSkip != 0) {
            for (size_t k = 0; k < nChunkElems; k++) afBuf[k] = static_cast<float>((n * 31 + k) % 97);
            for (int i = 0; i < nRank; i++) {
                anStart[i] = anScaled[i] * anChunk[i];
                anCount[i] = std::min(anChunk[i], anDims[i] - anStart[i]);
            }
            hid_t hFSel = H5Dget_space(hDset);
            H5Sselect_hyperslab(hFSel, H5S_SELECT_SET, anStart.data(), nullptr, anCount.data(), nullptr);
            hid_t hMSel = H5Scopy(hMem);
            std::vector<hsize_t> anZero(nRank, 0);
            H5Sselect_hyperslab(hMSel, H5S_SELECT_SET, anZero.data(), nullptr, anCount.data(), nullptr);
            Check(H5Dwrite(hDset, H5T_NATIVE_FLOAT, hMSel, hFSel, H5P_DEFAULT, afBuf.data()), "write");
            H5Sclose(hMSel);
            H5Sclose(hFSel);
        }
        for (int d = nRank - 1; d >= 0; d--) {
            if (++anScaled[d] < anGrid[d]) break;
            anScaled[d] = 0;
        }
    }
    H5Sclose(hMem);
    H5Dclose(hDset);
    H5Pclose(hDcpl);
    H5Sclose(hSpace);
}

void MakeContiguous(hid_t hLoc, const char *pszName, H5D_layout_t eLayout)
{
    hsize_t anDims[2] = {8, 8};
    hid_t hSpace = H5Screate_simple(2, anDims, nullptr);
    hid_t hDcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(hDcpl, eLayout);
    hid_t hDset = H5Dcreate2(hLoc, pszName, H5T_NATIVE_FLOAT, hSpace, H5P_DEFAULT, hDcpl, H5P_DEFAULT);
    std::vector<float> af(64, 1.0f);
    H5Dwrite(hDset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, af.data());
    H5Dclose(hDset);
    H5Pclose(hDcpl);
    H5Sclose(hSpace);
}

// Group with nLinks sibling subgroups (long names so the link heap grows indirect blocks) and
// datasets inserted at the beginning, middle and end of the name order.
hid_t MakeBigGroup(hid_t hParent, const char *pszName, int nLinks, bool bDense)
{
    hid_t hGcpl = H5Pcreate(H5P_GROUP_CREATE);
    if (bDense) H5Pset_link_phase_change(hGcpl, 0, 0);
    hid_t hGroup = H5Gcreate2(hParent, pszName, H5P_DEFAULT, hGcpl, H5P_DEFAULT);
    H5Pclose(hGcpl);
    for (int i = 0; i < nLinks; i++) {
        char szName[128];
        snprintf(szName, sizeof(szName), "sibling_group_with_a_rather_long_name_to_fill_heap_%06d", i);
        hid_t h = H5Gcreate2(hGroup, szName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(h);
    }
    return hGroup;
}

void WriteFile(const std::string &osPath, H5F_libver_t eLow, H5F_libver_t eHigh, hsize_t nUserBlock = 0)
{
    hid_t hFapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_libver_bounds(hFapl, eLow, eHigh);
    hid_t hFcpl = H5Pcreate(H5P_FILE_CREATE);
    if (nUserBlock) Check(H5Pset_userblock(hFcpl, nUserBlock), "set_userblock");
    hid_t hFile = H5Fcreate(osPath.c_str(), H5F_ACC_TRUNC, hFcpl, hFapl);
    Check(hFile, osPath.c_str());
    H5Pclose(hFapl);
    H5Pclose(hFcpl);

    const bool bNewGroups = eLow >= H5F_LIBVER_V18;
    hid_t hG1 = MakeBigGroup(hFile, "big", 20000, bNewGroups);
    hid_t hG2 = H5Gopen2(hG1, "sibling_group_with_a_rather_long_name_to_fill_heap_010000", H5P_DEFAULT);
    const hsize_t U = H5S_UNLIMITED;

    for (hid_t hLoc : {hG1, hG2, hFile}) {
        MakeChunked(hLoc, "aaa_btree_2d", {1000, 1500}, {1000, 1500}, {128, 128}, true, 5);
        MakeChunked(hLoc, "zzz_cube_3d", {7, 300, 400}, {7, 300, 400}, {1, 128, 128}, true, 4);
    }
    // Fixed dims: layout v4 uses a fixed array (paged when > 1024 chunks).
    MakeChunked(hFile, "fixed_small", {600, 700}, {600, 700}, {64, 64}, false, 3);
    MakeChunked(hFile, "fixed_small_filtered", {600, 700}, {600, 700}, {64, 64}, true, 3);
    MakeChunked(hFile, "fixed_paged_filtered", {2048, 2048}, {2048, 2048}, {32, 32}, true, 7);
    MakeChunked(hFile, "fixed_paged", {2048, 2048}, {2048, 2048}, {32, 32}, false, 0);
    // One unlimited dimension: extensible array (large enough for super blocks and paged data blocks).
    MakeChunked(hFile, "earray_filtered", {6000, 64, 64}, {U, 64, 64}, {1, 64, 64}, true, 9);
    MakeChunked(hFile, "earray", {40000, 4}, {U, 4}, {1, 4}, false, 11);
    MakeChunked(hFile, "earray_2d_inner_unlim", {40, 600}, {40, U}, {8, 8}, true, 6);
    // Two unlimited dimensions: v2 B-tree.
    MakeChunked(hFile, "btree2_filtered", {900, 900}, {U, U}, {32, 32}, true, 5);
    MakeChunked(hFile, "btree2", {900, 900}, {U, U}, {32, 32}, false, 5);
    // Single chunk (filtered and not).
    MakeChunked(hFile, "single", {300, 200}, {300, 200}, {300, 200}, false, 0);
    MakeChunked(hFile, "single_filtered", {300, 200}, {300, 200}, {300, 200}, true, 0);
    // Implicit index: early allocation, no filters.
    MakeChunked(hFile, "implicit", {500, 300}, {500, 300}, {64, 64}, false, 4, H5D_ALLOC_TIME_EARLY);
    // Not chunked: must be rejected with ERR_NOT_CHUNKED.
    MakeContiguous(hFile, "contiguous", H5D_CONTIGUOUS);
    MakeContiguous(hFile, "compact", H5D_COMPACT);

    H5Gclose(hG2);
    H5Gclose(hG1);
    H5Fclose(hFile);
    printf("wrote %s\n", osPath.c_str());
}

}  // namespace

int main(int argc, char **argv)
{
    if (argc != 2) { fprintf(stderr, "usage: nisar_chunkmap_synth OUTDIR\n"); return 2; }
    const std::string osDir = argv[1];
    WriteFile(osDir + "/synth_sb0_earliest.h5", H5F_LIBVER_EARLIEST, H5F_LIBVER_LATEST);
    WriteFile(osDir + "/synth_sb2_v18.h5", H5F_LIBVER_V18, H5F_LIBVER_V18);
    WriteFile(osDir + "/synth_sb3_v110.h5", H5F_LIBVER_V110, H5F_LIBVER_LATEST);
    WriteFile(osDir + "/synth_sb3_latest.h5", H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    WriteFile(osDir + "/synth_sb0_userblock512.h5", H5F_LIBVER_EARLIEST, H5F_LIBVER_LATEST, 512);
    WriteFile(osDir + "/synth_sb3_userblock4096.h5", H5F_LIBVER_V110, H5F_LIBVER_LATEST, 4096);
    return 0;
}
