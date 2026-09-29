// hdf5native.h
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

// Minimal read-only HDF5 parser for chunk-index discovery via GDAL VSI (feasibility spike).
// Unsupported file-format features return a non-OK Status so callers can fall back to libhdf5.

#ifndef NISAR_HDF5_NATIVE_H
#define NISAR_HDF5_NATIVE_H

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "cpl_vsi.h"

namespace NisarHDF5Native {

enum class Status : int {
    OK = 0,
    ERR_IO = 1,            // VSI open/read failure or read past EOF
    ERR_SIGNATURE = 2,     // missing or wrong structure signature
    ERR_SUPERBLOCK = 3,    // no superblock / unsupported superblock version
    ERR_OHDR = 4,          // unsupported object header version or encoding
    ERR_MESSAGE = 5,       // unrecognized, shared, or unsupported message encoding
    ERR_GROUP = 6,         // unsupported group storage or link type (soft/external/huge)
    ERR_NOT_FOUND = 7,     // path component does not exist
    ERR_NOT_DATASET = 8,   // object has no layout message
    ERR_NOT_CHUNKED = 9,   // dataset layout is contiguous / compact
    ERR_LAYOUT = 10,       // unsupported layout message version or class
    ERR_INDEX = 11,        // unsupported chunk index type or index encoding
    ERR_CORRUPT = 12,      // inconsistent structure or checksum mismatch
    ERR_LIMIT = 13,        // exceeded a sanity limit (depth, size, count)
};

const char *StatusToString(Status eStatus);

enum class ChunkIndex : int {
    NONE = 0,
    BTREE_V1,
    SINGLE_CHUNK,
    IMPLICIT,
    FIXED_ARRAY,
    EXTENSIBLE_ARRAY,
    BTREE_V2,
};

const char *ChunkIndexToString(ChunkIndex eIndex);

struct IOStats {
    uint64_t nReadCalls = 0;   // VSIFReadL() calls issued by the parser
    uint64_t nBytesRead = 0;   // bytes returned by those calls
};

struct DatasetInfo {
    std::string osPath;
    uint64_t nObjectHeaderAddr = 0;
    int nObjectHeaderVersion = 0;
    int nRank = 0;
    std::vector<uint64_t> anDims;
    std::vector<uint64_t> anMaxDims;       // UINT64_MAX == unlimited
    std::vector<uint64_t> anChunkDims;     // rank entries (element size excluded)
    uint32_t nElementSize = 0;
    int nTypeClass = -1;
    int nLayoutVersion = 0;
    int nLayoutClass = -1;                 // 0 compact, 1 contiguous, 2 chunked, 3 virtual
    ChunkIndex eIndex = ChunkIndex::NONE;
    uint64_t nIndexAddr = 0;
    uint8_t nLayoutFlags = 0;              // v4 chunked layout flags
    std::vector<uint16_t> anFilterIds;

    // Index-type specific parameters (layout v4)
    uint64_t nSingleFilteredSize = 0;
    uint32_t nSingleFilterMask = 0;
    uint8_t nFarrayPageBits = 0;
    uint8_t nEarrayMaxBits = 0;
    uint8_t nEarrayIdxBlkElmts = 0;
    uint8_t nEarrayMinPtrs = 0;
    uint8_t nEarrayMinElmts = 0;
    uint8_t nEarrayPageBits = 0;

    uint64_t ChunkBytes() const;           // uncompressed chunk size in bytes
};

struct ChunkRecord {
    std::vector<uint64_t> anOffset;        // element offsets, rank entries
    uint32_t nFilterMask = 0;
    uint64_t nAddr = 0;                    // file address, relative to base
    uint64_t nSize = 0;                    // stored bytes
};

class Reader;

class File {
  public:
    ~File();

    // nBlockSize: metadata read granularity in bytes (0 = exact-range reads).
    static Status Open(const std::string &osVSIPath, std::unique_ptr<File> &poFile,
                       size_t nBlockSize);

    Status OpenDataset(const std::string &osPath, DatasetInfo &oInfo);
    Status GetChunks(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks);

    const IOStats &GetStats() const;
    int GetSuperblockVersion() const { return m_nSuperblockVersion; }

    // Last failure detail (human readable).
    const std::string &GetLastError() const { return m_osLastError; }

  private:
    File() = default;

    struct Message {
        uint16_t nType = 0;
        uint8_t nFlags = 0;
        std::vector<uint8_t> abyData;
    };

    struct LocalHeap {
        std::vector<uint8_t> abyData;
    };

    struct FractalHeap {
        uint16_t nIdLen = 0;
        uint16_t nFilterLen = 0;
        uint8_t nFlags = 0;
        uint32_t nMaxManSize = 0;
        uint16_t nTableWidth = 0;
        uint64_t nStartBlockSize = 0;
        uint64_t nMaxDirectBlockSize = 0;
        uint16_t nMaxHeapBits = 0;
        uint64_t nRootAddr = 0;
        uint16_t nRootRows = 0;
        unsigned nHeapOffSize = 0;
        unsigned nHeapLenSize = 0;
        unsigned nMaxDirectRows = 0;
    };

    struct BTree2Header {
        uint8_t nType = 0;
        uint32_t nNodeSize = 0;
        uint16_t nRecordSize = 0;
        uint16_t nDepth = 0;
        uint64_t nRootAddr = 0;
        uint16_t nRootRecords = 0;
        uint64_t nTotalRecords = 0;
        unsigned nMaxNrecSize = 0;
        std::vector<unsigned> anCumMaxNrecSize;   // per depth
    };

    Status Fail(Status eStatus, const std::string &osMsg);

    Status ReadSuperblock();
    Status ReadObjectHeader(uint64_t nAddr, std::vector<Message> &aoMsgs, int *pnVersion);
    Status LookupChild(const std::vector<Message> &aoGroupMsgs, const std::string &osName,
                       uint64_t &nChildAddr);
    Status LookupSymbolTable(uint64_t nBTreeAddr, uint64_t nHeapAddr, const std::string &osName,
                             uint64_t &nChildAddr);
    Status ReadLocalHeap(uint64_t nAddr, const LocalHeap *&poHeap);
    Status LookupCompact(const std::vector<Message> &aoGroupMsgs, const std::string &osName,
                         uint64_t &nChildAddr);
    Status LookupDense(uint64_t nHeapAddr, uint64_t nNameIndexAddr, const std::string &osName,
                       uint64_t &nChildAddr);
    Status DecodeLinkMessage(const uint8_t *pabyData, size_t nLen, std::string &osName,
                             int &nLinkType, uint64_t &nAddr);

    Status ReadFractalHeapHeader(uint64_t nAddr, FractalHeap &oHeap);
    Status ReadFractalHeapObject(const FractalHeap &oHeap, const uint8_t *pabyId,
                                 std::vector<uint8_t> &abyObj);

    Status ReadBTree2Header(uint64_t nAddr, BTree2Header &oHdr);
    Status CollectBTree2Records(const BTree2Header &oHdr, uint64_t nNodeAddr, unsigned nDepth,
                                unsigned nRecords, std::vector<std::vector<uint8_t>> &aoRecords,
                                int nDepthGuard);

    Status DecodeDataset(const std::vector<Message> &aoMsgs, DatasetInfo &oInfo);

    Status IterBTreeV1(const DatasetInfo &oInfo, uint64_t nAddr, int nExpectedLevel,
                       std::vector<ChunkRecord> &aoChunks, int nDepthGuard);
    Status IterSingle(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks);
    Status IterImplicit(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks);
    Status IterFixedArray(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks);
    Status IterExtensibleArray(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks);
    Status IterBTreeV2(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks);

    std::unique_ptr<Reader> m_poReader;
    int m_nSuperblockVersion = -1;
    unsigned m_nSizeofAddr = 8;
    unsigned m_nSizeofLen = 8;
    uint64_t m_nBaseAddr = 0;
    uint64_t m_nRootAddr = 0;
    std::string m_osLastError;
    std::map<uint64_t, LocalHeap> m_oLocalHeaps;
};

}  // namespace NisarHDF5Native

#endif  // NISAR_HDF5_NATIVE_H
