// hdf5native.cpp
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

// Structure layouts follow the "HDF5 File Format Specification Version 3.0".

#include "hdf5native.h"

#include <algorithm>
#include <cstring>

namespace NisarHDF5Native {

static constexpr uint64_t UNDEF_ADDR = ~static_cast<uint64_t>(0);
static constexpr uint64_t UNLIMITED = ~static_cast<uint64_t>(0);
static constexpr int MAX_DEPTH = 64;
static constexpr size_t MAX_OHDR_CHUNKS = 4096;
static constexpr size_t MAX_METADATA_READ = 256 * 1024 * 1024;
static constexpr size_t MAX_CACHE_BYTES = 256 * 1024 * 1024;

// HDF5 object header message types
enum : uint16_t {
    MSG_NIL = 0x00,
    MSG_DATASPACE = 0x01,
    MSG_LINK_INFO = 0x02,
    MSG_DATATYPE = 0x03,
    MSG_FILL_OLD = 0x04,
    MSG_FILL = 0x05,
    MSG_LINK = 0x06,
    MSG_EXTERNAL_FILES = 0x07,
    MSG_LAYOUT = 0x08,
    MSG_BOGUS = 0x09,
    MSG_GROUP_INFO = 0x0A,
    MSG_FILTER_PIPELINE = 0x0B,
    MSG_ATTRIBUTE = 0x0C,
    MSG_COMMENT = 0x0D,
    MSG_MTIME_OLD = 0x0E,
    MSG_SHMESG_TABLE = 0x0F,
    MSG_CONTINUATION = 0x10,
    MSG_SYMBOL_TABLE = 0x11,
    MSG_MTIME = 0x12,
    MSG_BTREE_K = 0x13,
    MSG_DRIVER_INFO = 0x14,
    MSG_ATTRIBUTE_INFO = 0x15,
    MSG_REFCOUNT = 0x16,
    MSG_FSINFO = 0x17,
    MSG_LAST_KNOWN = MSG_FSINFO,
};

static constexpr uint8_t MSG_FLAG_SHARED = 0x02;

const char *StatusToString(Status eStatus)
{
    switch (eStatus) {
        case Status::OK: return "OK";
        case Status::ERR_IO: return "ERR_IO";
        case Status::ERR_SIGNATURE: return "ERR_SIGNATURE";
        case Status::ERR_SUPERBLOCK: return "ERR_SUPERBLOCK";
        case Status::ERR_OHDR: return "ERR_OHDR";
        case Status::ERR_MESSAGE: return "ERR_MESSAGE";
        case Status::ERR_GROUP: return "ERR_GROUP";
        case Status::ERR_NOT_FOUND: return "ERR_NOT_FOUND";
        case Status::ERR_NOT_DATASET: return "ERR_NOT_DATASET";
        case Status::ERR_NOT_CHUNKED: return "ERR_NOT_CHUNKED";
        case Status::ERR_LAYOUT: return "ERR_LAYOUT";
        case Status::ERR_INDEX: return "ERR_INDEX";
        case Status::ERR_CORRUPT: return "ERR_CORRUPT";
        case Status::ERR_LIMIT: return "ERR_LIMIT";
    }
    return "ERR_UNKNOWN";
}

const char *ChunkIndexToString(ChunkIndex eIndex)
{
    switch (eIndex) {
        case ChunkIndex::NONE: return "none";
        case ChunkIndex::BTREE_V1: return "btree_v1";
        case ChunkIndex::SINGLE_CHUNK: return "single_chunk";
        case ChunkIndex::IMPLICIT: return "implicit";
        case ChunkIndex::FIXED_ARRAY: return "fixed_array";
        case ChunkIndex::EXTENSIBLE_ARRAY: return "extensible_array";
        case ChunkIndex::BTREE_V2: return "btree_v2";
    }
    return "unknown";
}

uint64_t DatasetInfo::ChunkBytes() const
{
    uint64_t n = nElementSize;
    for (uint64_t d : anChunkDims) n *= d;
    return n;
}

/************************************************************************/
/*                         Helpers                                      */
/************************************************************************/

static unsigned Log2Gen(uint64_t n)
{
    unsigned r = 0;
    while (n >>= 1) r++;
    return r;
}

// Number of bytes needed to encode values up to n (H5VM_limit_enc_size).
static unsigned LimitEncSize(uint64_t n)
{
    return Log2Gen(n) / 8 + 1;
}

static inline uint32_t Rot(uint32_t x, int k)
{
    return (x << k) | (x >> (32 - k));
}

// Bob Jenkins' lookup3 hashlittle(), byte-wise variant (H5_checksum_lookup3).
static uint32_t Lookup3(const uint8_t *k, size_t length, uint32_t initval)
{
    uint32_t a, b, c;
    a = b = c = 0xdeadbeef + static_cast<uint32_t>(length) + initval;
    while (length > 12) {
        a += k[0] + (static_cast<uint32_t>(k[1]) << 8) + (static_cast<uint32_t>(k[2]) << 16) +
             (static_cast<uint32_t>(k[3]) << 24);
        b += k[4] + (static_cast<uint32_t>(k[5]) << 8) + (static_cast<uint32_t>(k[6]) << 16) +
             (static_cast<uint32_t>(k[7]) << 24);
        c += k[8] + (static_cast<uint32_t>(k[9]) << 8) + (static_cast<uint32_t>(k[10]) << 16) +
             (static_cast<uint32_t>(k[11]) << 24);
        a -= c; a ^= Rot(c, 4);  c += b;
        b -= a; b ^= Rot(a, 6);  a += c;
        c -= b; c ^= Rot(b, 8);  b += a;
        a -= c; a ^= Rot(c, 16); c += b;
        b -= a; b ^= Rot(a, 19); a += c;
        c -= b; c ^= Rot(b, 4);  b += a;
        length -= 12;
        k += 12;
    }
    switch (length) {
        case 12: c += static_cast<uint32_t>(k[11]) << 24; [[fallthrough]];
        case 11: c += static_cast<uint32_t>(k[10]) << 16; [[fallthrough]];
        case 10: c += static_cast<uint32_t>(k[9]) << 8; [[fallthrough]];
        case 9: c += k[8]; [[fallthrough]];
        case 8: b += static_cast<uint32_t>(k[7]) << 24; [[fallthrough]];
        case 7: b += static_cast<uint32_t>(k[6]) << 16; [[fallthrough]];
        case 6: b += static_cast<uint32_t>(k[5]) << 8; [[fallthrough]];
        case 5: b += k[4]; [[fallthrough]];
        case 4: a += static_cast<uint32_t>(k[3]) << 24; [[fallthrough]];
        case 3: a += static_cast<uint32_t>(k[2]) << 16; [[fallthrough]];
        case 2: a += static_cast<uint32_t>(k[1]) << 8; [[fallthrough]];
        case 1: a += k[0]; break;
        case 0: return c;
    }
    c ^= b; c -= Rot(b, 14);
    a ^= c; a -= Rot(c, 11);
    b ^= a; b -= Rot(a, 25);
    c ^= b; c -= Rot(b, 16);
    a ^= c; a -= Rot(c, 4);
    b ^= a; b -= Rot(a, 14);
    c ^= b; c -= Rot(b, 24);
    return c;
}

// Bounds-checked little-endian decoder over a byte buffer.
class Cursor {
  public:
    Cursor(const uint8_t *p, size_t n) : m_p(p), m_n(n) {}
    explicit Cursor(const std::vector<uint8_t> &v) : m_p(v.data()), m_n(v.size()) {}

    bool Ok() const { return m_bOk; }
    size_t Pos() const { return m_nPos; }
    size_t Remaining() const { return m_bOk ? m_n - m_nPos : 0; }
    const uint8_t *Ptr() const { return m_p + m_nPos; }

    uint64_t U(size_t nBytes)
    {
        if (!m_bOk || nBytes > 8 || m_nPos + nBytes > m_n) { m_bOk = false; return 0; }
        uint64_t v = 0;
        for (size_t i = 0; i < nBytes; i++) v |= static_cast<uint64_t>(m_p[m_nPos + i]) << (8 * i);
        m_nPos += nBytes;
        return v;
    }
    uint8_t U8() { return static_cast<uint8_t>(U(1)); }
    uint16_t U16() { return static_cast<uint16_t>(U(2)); }
    uint32_t U32() { return static_cast<uint32_t>(U(4)); }
    uint64_t U64() { return U(8); }

    uint64_t Addr(unsigned nSize)
    {
        uint64_t v = U(nSize);
        uint64_t nAllOnes = nSize >= 8 ? UNDEF_ADDR : ((static_cast<uint64_t>(1) << (8 * nSize)) - 1);
        return v == nAllOnes ? UNDEF_ADDR : v;
    }

    void Skip(size_t n)
    {
        if (!m_bOk || m_nPos + n > m_n) { m_bOk = false; return; }
        m_nPos += n;
    }

    bool Sig(const char *pszSig)
    {
        if (!m_bOk || m_nPos + 4 > m_n) { m_bOk = false; return false; }
        bool b = memcmp(m_p + m_nPos, pszSig, 4) == 0;
        m_nPos += 4;
        return b;
    }

  private:
    const uint8_t *m_p;
    size_t m_n;
    size_t m_nPos = 0;
    bool m_bOk = true;
};

/************************************************************************/
/*                         Reader                                       */
/************************************************************************/

// VSI reader with an optional aligned block cache; counts every VSIFReadL().
class Reader {
  public:
    Reader(VSILFILE *fp, uint64_t nFileSize, size_t nBlockSize)
        : m_fp(fp), m_nFileSize(nFileSize), m_nBlockSize(nBlockSize) {}
    ~Reader() { if (m_fp) VSIFCloseL(m_fp); }

    uint64_t FileSize() const { return m_nFileSize; }
    const IOStats &Stats() const { return m_oStats; }

    bool Read(uint64_t nOff, size_t nLen, std::vector<uint8_t> &abyOut)
    {
        abyOut.resize(nLen);
        if (nLen == 0) return true;
        if (nLen > MAX_METADATA_READ || nOff >= m_nFileSize || nLen > m_nFileSize - nOff) return false;
        if (m_nBlockSize == 0) return RawRead(nOff, nLen, abyOut.data());

        if (m_nCacheBytes > MAX_CACHE_BYTES) { m_oBlocks.clear(); m_nCacheBytes = 0; }

        const uint64_t nFirst = nOff / m_nBlockSize;
        const uint64_t nLast = (nOff + nLen - 1) / m_nBlockSize;
        uint64_t iBlock = nFirst;
        while (iBlock <= nLast) {
            if (m_oBlocks.count(iBlock)) { iBlock++; continue; }
            uint64_t iRunEnd = iBlock;
            while (iRunEnd + 1 <= nLast && !m_oBlocks.count(iRunEnd + 1)) iRunEnd++;
            const uint64_t nRunOff = iBlock * m_nBlockSize;
            const uint64_t nRunLen = std::min<uint64_t>((iRunEnd - iBlock + 1) * m_nBlockSize, m_nFileSize - nRunOff);
            std::vector<uint8_t> abyRun(static_cast<size_t>(nRunLen));
            if (!RawRead(nRunOff, abyRun.size(), abyRun.data())) return false;
            for (uint64_t b = iBlock; b <= iRunEnd; b++) {
                const size_t nStart = static_cast<size_t>((b - iBlock) * m_nBlockSize);
                const size_t nSz = std::min<size_t>(m_nBlockSize, abyRun.size() - nStart);
                m_oBlocks[b].assign(abyRun.begin() + nStart, abyRun.begin() + nStart + nSz);
                m_nCacheBytes += nSz;
            }
            iBlock = iRunEnd + 1;
        }

        size_t nCopied = 0;
        for (uint64_t b = nFirst; b <= nLast; b++) {
            const std::vector<uint8_t> &abyBlock = m_oBlocks[b];
            const uint64_t nBlockOff = b * m_nBlockSize;
            const size_t nSrcStart = (b == nFirst) ? static_cast<size_t>(nOff - nBlockOff) : 0;
            const size_t nAvail = abyBlock.size() > nSrcStart ? abyBlock.size() - nSrcStart : 0;
            const size_t nTake = std::min(nAvail, nLen - nCopied);
            if (nTake == 0) return false;
            memcpy(abyOut.data() + nCopied, abyBlock.data() + nSrcStart, nTake);
            nCopied += nTake;
        }
        return nCopied == nLen;
    }

  private:
    bool RawRead(uint64_t nOff, size_t nLen, uint8_t *pDst)
    {
        if (VSIFSeekL(m_fp, static_cast<vsi_l_offset>(nOff), SEEK_SET) != 0) return false;
        const size_t nGot = VSIFReadL(pDst, 1, nLen, m_fp);
        m_oStats.nReadCalls++;
        m_oStats.nBytesRead += nGot;
        return nGot == nLen;
    }

    VSILFILE *m_fp;
    uint64_t m_nFileSize;
    size_t m_nBlockSize;
    std::map<uint64_t, std::vector<uint8_t>> m_oBlocks;
    size_t m_nCacheBytes = 0;
    IOStats m_oStats;
};

/************************************************************************/
/*                         File                                         */
/************************************************************************/

File::~File() = default;

const IOStats &File::GetStats() const
{
    return m_poReader->Stats();
}

Status File::Fail(Status eStatus, const std::string &osMsg)
{
    m_osLastError = osMsg;
    return eStatus;
}

#define READ_OR_FAIL(off, len, buf)                                                          \
    do {                                                                                     \
        if (!m_poReader->Read((off), (len), (buf)))                                          \
            return Fail(Status::ERR_IO, "read failed at " + std::to_string(off) + " (" +     \
                                            std::to_string(len) + " bytes)");                \
    } while (0)

Status File::Open(const std::string &osVSIPath, std::unique_ptr<File> &poFile, size_t nBlockSize)
{
    poFile.reset();
    VSILFILE *fp = VSIFOpenL(osVSIPath.c_str(), "rb");
    if (!fp) return Status::ERR_IO;
    if (VSIFSeekL(fp, 0, SEEK_END) != 0) { VSIFCloseL(fp); return Status::ERR_IO; }
    const uint64_t nFileSize = static_cast<uint64_t>(VSIFTellL(fp));

    std::unique_ptr<File> poNew(new File());
    poNew->m_poReader.reset(new Reader(fp, nFileSize, nBlockSize));
    const Status eStatus = poNew->ReadSuperblock();
    if (eStatus != Status::OK) {
        poFile = std::move(poNew);  // keep for error string / stats
        return eStatus;
    }
    poFile = std::move(poNew);
    return Status::OK;
}

Status File::ReadSuperblock()
{
    static const uint8_t abySig[8] = {0x89, 'H', 'D', 'F', '\r', '\n', 0x1a, '\n'};
    std::vector<uint8_t> aby;
    uint64_t nSBAddr = UNDEF_ADDR;
    for (uint64_t nTry = 0; nTry < m_poReader->FileSize(); nTry = (nTry == 0) ? 512 : nTry * 2) {
        const size_t nLen = static_cast<size_t>(std::min<uint64_t>(256, m_poReader->FileSize() - nTry));
        if (nLen < 8) break;
        READ_OR_FAIL(nTry, nLen, aby);
        if (memcmp(aby.data(), abySig, 8) == 0) { nSBAddr = nTry; break; }
    }
    if (nSBAddr == UNDEF_ADDR) return Fail(Status::ERR_SUPERBLOCK, "HDF5 signature not found");

    Cursor c(aby);
    c.Skip(8);
    m_nSuperblockVersion = c.U8();
    if (m_nSuperblockVersion == 0 || m_nSuperblockVersion == 1) {
        c.U8();  // free-space version
        c.U8();  // root group symbol table entry version
        c.U8();  // reserved
        c.U8();  // shared header message format version
        m_nSizeofAddr = c.U8();
        m_nSizeofLen = c.U8();
        c.U8();  // reserved
        c.U16(); // group leaf node K
        c.U16(); // group internal node K
        c.U32(); // file consistency flags
        if (m_nSuperblockVersion == 1) { c.U16(); c.U16(); }
        if (m_nSizeofAddr < 2 || m_nSizeofAddr > 8 || m_nSizeofLen < 2 || m_nSizeofLen > 8)
            return Fail(Status::ERR_SUPERBLOCK, "bad sizeof offsets/lengths");
        m_nBaseAddr = c.U(m_nSizeofAddr);
        c.Addr(m_nSizeofAddr);  // free-space info address
        c.Addr(m_nSizeofAddr);  // end of file address
        c.Addr(m_nSizeofAddr);  // driver information block address
        c.U(m_nSizeofAddr);     // root symbol table entry: link name offset
        m_nRootAddr = c.Addr(m_nSizeofAddr);
    }
    else if (m_nSuperblockVersion == 2 || m_nSuperblockVersion == 3) {
        m_nSizeofAddr = c.U8();
        m_nSizeofLen = c.U8();
        c.U8();  // file consistency flags
        if (m_nSizeofAddr < 2 || m_nSizeofAddr > 8 || m_nSizeofLen < 2 || m_nSizeofLen > 8)
            return Fail(Status::ERR_SUPERBLOCK, "bad sizeof offsets/lengths");
        m_nBaseAddr = c.U(m_nSizeofAddr);
        c.Addr(m_nSizeofAddr);  // superblock extension address
        c.Addr(m_nSizeofAddr);  // end of file address
        m_nRootAddr = c.Addr(m_nSizeofAddr);
        const size_t nChecked = c.Pos();
        const uint32_t nStored = c.U32();
        if (!c.Ok()) return Fail(Status::ERR_SUPERBLOCK, "truncated superblock");
        if (Lookup3(aby.data(), nChecked, 0) != nStored)
            return Fail(Status::ERR_CORRUPT, "superblock checksum mismatch");
    }
    else {
        return Fail(Status::ERR_SUPERBLOCK, "unsupported superblock version " +
                                                std::to_string(m_nSuperblockVersion));
    }
    if (!c.Ok() || m_nRootAddr == UNDEF_ADDR) return Fail(Status::ERR_SUPERBLOCK, "truncated superblock");
    return Status::OK;
}

/************************************************************************/
/*                         Object headers                               */
/************************************************************************/

Status File::ReadObjectHeader(uint64_t nAddr, std::vector<Message> &aoMsgs, int *pnVersion)
{
    aoMsgs.clear();
    if (nAddr == UNDEF_ADDR) return Fail(Status::ERR_OHDR, "undefined object header address");
    const uint64_t nAbs = m_nBaseAddr + nAddr;

    std::vector<uint8_t> aby;
    READ_OR_FAIL(nAbs, 16, aby);

    struct Pending { uint64_t nAddr; uint64_t nLen; };
    std::vector<Pending> aoPending;

    const bool bV2 = memcmp(aby.data(), "OHDR", 4) == 0;
    if (!bV2) {
        // Version 1: version, reserved, #messages(2), refcount(4), header size(4), pad(4)
        Cursor c(aby);
        const uint8_t nVersion = c.U8();
        if (nVersion != 1) return Fail(Status::ERR_OHDR, "unsupported object header version " + std::to_string(nVersion));
        if (pnVersion) *pnVersion = 1;
        c.U8();
        c.U16();
        c.U32();
        const uint32_t nHdrSize = c.U32();
        std::vector<uint8_t> abyChunk;
        READ_OR_FAIL(nAbs + 16, nHdrSize, abyChunk);
        size_t iChunk = 0;
        for (;;) {
            Cursor m(abyChunk);
            while (m.Remaining() >= 8) {
                const uint16_t nType = m.U16();
                const uint16_t nSize = m.U16();
                const uint8_t nFlags = m.U8();
                m.Skip(3);
                if (m.Remaining() < nSize) return Fail(Status::ERR_CORRUPT, "v1 message overruns chunk");
                Message oMsg;
                oMsg.nType = nType;
                oMsg.nFlags = nFlags;
                oMsg.abyData.assign(m.Ptr(), m.Ptr() + nSize);
                m.Skip(nSize);
                if (nType > MSG_LAST_KNOWN)
                    return Fail(Status::ERR_MESSAGE, "unrecognized message type " + std::to_string(nType));
                if (nType == MSG_CONTINUATION) {
                    Cursor cc(oMsg.abyData);
                    Pending p;
                    p.nAddr = cc.Addr(m_nSizeofAddr);
                    p.nLen = cc.U(m_nSizeofLen);
                    if (!cc.Ok() || p.nAddr == UNDEF_ADDR) return Fail(Status::ERR_CORRUPT, "bad continuation message");
                    aoPending.push_back(p);
                }
                aoMsgs.push_back(std::move(oMsg));
            }
            if (iChunk >= aoPending.size()) break;
            if (iChunk >= MAX_OHDR_CHUNKS) return Fail(Status::ERR_LIMIT, "too many object header chunks");
            READ_OR_FAIL(m_nBaseAddr + aoPending[iChunk].nAddr, static_cast<size_t>(aoPending[iChunk].nLen), abyChunk);
            iChunk++;
        }
        return Status::OK;
    }

    // Version 2
    Cursor c(aby);
    c.Skip(4);
    const uint8_t nVersion = c.U8();
    if (nVersion != 2) return Fail(Status::ERR_OHDR, "unsupported OHDR version " + std::to_string(nVersion));
    if (pnVersion) *pnVersion = 2;
    const uint8_t nFlags = c.U8();
    size_t nPrefix = 6;
    if (nFlags & 0x20) nPrefix += 16;
    if (nFlags & 0x10) nPrefix += 4;
    const size_t nSizeBytes = static_cast<size_t>(1) << (nFlags & 0x03);
    std::vector<uint8_t> abyPrefix;
    READ_OR_FAIL(nAbs, nPrefix + nSizeBytes, abyPrefix);
    Cursor cp(abyPrefix);
    cp.Skip(nPrefix);
    const uint64_t nChunk0 = cp.U(nSizeBytes);
    nPrefix += nSizeBytes;
    const bool bTrackCorder = (nFlags & 0x04) != 0;
    const size_t nMsgHdr = bTrackCorder ? 6 : 4;

    std::vector<uint8_t> abyChunk;
    READ_OR_FAIL(nAbs, static_cast<size_t>(nPrefix + nChunk0 + 4), abyChunk);
    {
        Cursor ck(abyChunk.data() + nPrefix + nChunk0, 4);
        if (Lookup3(abyChunk.data(), static_cast<size_t>(nPrefix + nChunk0), 0) != ck.U32())
            return Fail(Status::ERR_CORRUPT, "OHDR checksum mismatch");
    }
    size_t nStart = nPrefix;
    size_t nEnd = static_cast<size_t>(nPrefix + nChunk0);
    size_t iChunk = 0;
    for (;;) {
        Cursor m(abyChunk.data() + nStart, nEnd - nStart);
        while (m.Remaining() >= nMsgHdr) {
            const uint8_t nType = m.U8();
            const uint16_t nSize = m.U16();
            const uint8_t nMsgFlags = m.U8();
            if (bTrackCorder) m.U16();
            if (m.Remaining() < nSize) return Fail(Status::ERR_CORRUPT, "v2 message overruns chunk");
            Message oMsg;
            oMsg.nType = nType;
            oMsg.nFlags = nMsgFlags;
            oMsg.abyData.assign(m.Ptr(), m.Ptr() + nSize);
            m.Skip(nSize);
            if (nType > MSG_LAST_KNOWN)
                return Fail(Status::ERR_MESSAGE, "unrecognized message type " + std::to_string(nType));
            if (nType == MSG_CONTINUATION) {
                Cursor cc(oMsg.abyData);
                Pending p;
                p.nAddr = cc.Addr(m_nSizeofAddr);
                p.nLen = cc.U(m_nSizeofLen);
                if (!cc.Ok() || p.nAddr == UNDEF_ADDR || p.nLen < 8) return Fail(Status::ERR_CORRUPT, "bad continuation message");
                aoPending.push_back(p);
            }
            aoMsgs.push_back(std::move(oMsg));
        }
        if (iChunk >= aoPending.size()) break;
        if (iChunk >= MAX_OHDR_CHUNKS) return Fail(Status::ERR_LIMIT, "too many object header chunks");
        const Pending &p = aoPending[iChunk++];
        READ_OR_FAIL(m_nBaseAddr + p.nAddr, static_cast<size_t>(p.nLen), abyChunk);
        if (memcmp(abyChunk.data(), "OCHK", 4) != 0) return Fail(Status::ERR_SIGNATURE, "missing OCHK signature");
        Cursor ck(abyChunk.data() + abyChunk.size() - 4, 4);
        if (Lookup3(abyChunk.data(), abyChunk.size() - 4, 0) != ck.U32())
            return Fail(Status::ERR_CORRUPT, "OCHK checksum mismatch");
        nStart = 4;
        nEnd = abyChunk.size() - 4;
    }
    return Status::OK;
}

/************************************************************************/
/*                         Groups                                       */
/************************************************************************/

Status File::LookupChild(const std::vector<Message> &aoGroupMsgs, const std::string &osName,
                         uint64_t &nChildAddr)
{
    for (const Message &oMsg : aoGroupMsgs) {
        if (oMsg.nType == MSG_SYMBOL_TABLE) {
            if (oMsg.nFlags & MSG_FLAG_SHARED) return Fail(Status::ERR_MESSAGE, "shared symbol table message");
            Cursor c(oMsg.abyData);
            const uint64_t nBTree = c.Addr(m_nSizeofAddr);
            const uint64_t nHeap = c.Addr(m_nSizeofAddr);
            if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated symbol table message");
            return LookupSymbolTable(nBTree, nHeap, osName, nChildAddr);
        }
    }
    for (const Message &oMsg : aoGroupMsgs) {
        if (oMsg.nType == MSG_LINK_INFO) {
            if (oMsg.nFlags & MSG_FLAG_SHARED) return Fail(Status::ERR_MESSAGE, "shared link info message");
            Cursor c(oMsg.abyData);
            const uint8_t nVersion = c.U8();
            const uint8_t nFlags = c.U8();
            if (nVersion != 0) return Fail(Status::ERR_MESSAGE, "unsupported link info version");
            if (nFlags & 0x01) c.U64();
            const uint64_t nHeap = c.Addr(m_nSizeofAddr);
            const uint64_t nNameIdx = c.Addr(m_nSizeofAddr);
            if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated link info message");
            if (nHeap != UNDEF_ADDR) return LookupDense(nHeap, nNameIdx, osName, nChildAddr);
            break;
        }
    }
    return LookupCompact(aoGroupMsgs, osName, nChildAddr);
}

Status File::ReadLocalHeap(uint64_t nAddr, const LocalHeap *&poHeap)
{
    auto it = m_oLocalHeaps.find(nAddr);
    if (it != m_oLocalHeaps.end()) { poHeap = &it->second; return Status::OK; }
    std::vector<uint8_t> aby;
    const size_t nHdr = 8 + 2 * m_nSizeofLen + m_nSizeofAddr;
    READ_OR_FAIL(m_nBaseAddr + nAddr, nHdr, aby);
    Cursor c(aby);
    if (!c.Sig("HEAP")) return Fail(Status::ERR_SIGNATURE, "missing HEAP signature");
    if (c.U8() != 0) return Fail(Status::ERR_GROUP, "unsupported local heap version");
    c.Skip(3);
    const uint64_t nSegSize = c.U(m_nSizeofLen);
    c.U(m_nSizeofLen);
    const uint64_t nSegAddr = c.Addr(m_nSizeofAddr);
    if (!c.Ok() || nSegAddr == UNDEF_ADDR) return Fail(Status::ERR_CORRUPT, "bad local heap header");
    LocalHeap oHeap;
    READ_OR_FAIL(m_nBaseAddr + nSegAddr, static_cast<size_t>(nSegSize), oHeap.abyData);
    poHeap = &(m_oLocalHeaps[nAddr] = std::move(oHeap));
    return Status::OK;
}

static bool HeapString(const std::vector<uint8_t> &abyHeap, uint64_t nOff, std::string &os)
{
    if (nOff >= abyHeap.size()) return false;
    const uint8_t *p = abyHeap.data() + nOff;
    const void *pEnd = memchr(p, 0, abyHeap.size() - static_cast<size_t>(nOff));
    if (!pEnd) return false;
    os.assign(reinterpret_cast<const char *>(p), static_cast<const uint8_t *>(pEnd) - p);
    return true;
}

Status File::LookupSymbolTable(uint64_t nBTreeAddr, uint64_t nHeapAddr, const std::string &osName,
                               uint64_t &nChildAddr)
{
    const LocalHeap *poHeap = nullptr;
    Status eStatus = ReadLocalHeap(nHeapAddr, poHeap);
    if (eStatus != Status::OK) return eStatus;

    uint64_t nNode = nBTreeAddr;
    const size_t nNodeHdr = 8 + 2 * m_nSizeofAddr;
    std::vector<uint8_t> aby;
    for (int nGuard = 0; nGuard < MAX_DEPTH; nGuard++) {
        READ_OR_FAIL(m_nBaseAddr + nNode, nNodeHdr, aby);
        Cursor c(aby);
        if (!c.Sig("TREE")) return Fail(Status::ERR_SIGNATURE, "missing TREE signature");
        const uint8_t nType = c.U8();
        const uint8_t nLevel = c.U8();
        const uint16_t nEntries = c.U16();
        if (nType != 0) return Fail(Status::ERR_CORRUPT, "expected group B-tree node");
        const size_t nBody = static_cast<size_t>(nEntries) * (m_nSizeofLen + m_nSizeofAddr) + m_nSizeofLen;
        READ_OR_FAIL(m_nBaseAddr + nNode + nNodeHdr, nBody, aby);
        Cursor b(aby);
        std::vector<uint64_t> anKeys(nEntries + 1), anChildren(nEntries);
        for (unsigned i = 0; i < nEntries; i++) {
            anKeys[i] = b.U(m_nSizeofLen);
            anChildren[i] = b.Addr(m_nSizeofAddr);
        }
        anKeys[nEntries] = b.U(m_nSizeofLen);
        if (!b.Ok()) return Fail(Status::ERR_CORRUPT, "truncated group B-tree node");

        // Child i holds names in (key[i], key[i+1]].
        int iFound = -1;
        for (unsigned i = 0; i < nEntries; i++) {
            std::string osRight;
            if (!HeapString(poHeap->abyData, anKeys[i + 1], osRight))
                return Fail(Status::ERR_CORRUPT, "bad group B-tree key");
            if (osName.compare(osRight) <= 0) {
                if (i > 0) {
                    std::string osLeft;
                    if (!HeapString(poHeap->abyData, anKeys[i], osLeft))
                        return Fail(Status::ERR_CORRUPT, "bad group B-tree key");
                    if (osName.compare(osLeft) <= 0) break;
                }
                iFound = static_cast<int>(i);
                break;
            }
        }
        if (iFound < 0) return Fail(Status::ERR_NOT_FOUND, "'" + osName + "' not found");

        if (nLevel > 0) { nNode = anChildren[iFound]; continue; }

        // Symbol table node
        const uint64_t nSnod = anChildren[iFound];
        READ_OR_FAIL(m_nBaseAddr + nSnod, 8, aby);
        Cursor s(aby);
        if (!s.Sig("SNOD")) return Fail(Status::ERR_SIGNATURE, "missing SNOD signature");
        if (s.U8() != 1) return Fail(Status::ERR_GROUP, "unsupported SNOD version");
        s.U8();
        const uint16_t nSyms = s.U16();
        const size_t nEntrySize = 2 * m_nSizeofAddr + 24;
        READ_OR_FAIL(m_nBaseAddr + nSnod + 8, nSyms * nEntrySize, aby);
        Cursor e(aby);
        for (unsigned i = 0; i < nSyms; i++) {
            const uint64_t nNameOff = e.U(m_nSizeofAddr);
            const uint64_t nOH = e.Addr(m_nSizeofAddr);
            const uint32_t nCacheType = e.U32();
            e.Skip(4 + 16);
            std::string osEntry;
            if (!HeapString(poHeap->abyData, nNameOff, osEntry))
                return Fail(Status::ERR_CORRUPT, "bad symbol table entry name");
            if (osEntry == osName) {
                if (nCacheType == 2) return Fail(Status::ERR_GROUP, "soft link in symbol table");
                nChildAddr = nOH;
                return Status::OK;
            }
        }
        return Fail(Status::ERR_NOT_FOUND, "'" + osName + "' not found");
    }
    return Fail(Status::ERR_LIMIT, "group B-tree too deep");
}

Status File::DecodeLinkMessage(const uint8_t *pabyData, size_t nLen, std::string &osName,
                               int &nLinkType, uint64_t &nAddr)
{
    Cursor c(pabyData, nLen);
    if (c.U8() != 1) return Fail(Status::ERR_MESSAGE, "unsupported link message version");
    const uint8_t nFlags = c.U8();
    nLinkType = 0;
    if (nFlags & 0x08) nLinkType = c.U8();
    if (nFlags & 0x04) c.U64();
    if (nFlags & 0x10) c.U8();
    const uint64_t nNameLen = c.U(static_cast<size_t>(1) << (nFlags & 0x03));
    if (!c.Ok() || nNameLen > c.Remaining()) return Fail(Status::ERR_CORRUPT, "bad link message");
    osName.assign(reinterpret_cast<const char *>(c.Ptr()), static_cast<size_t>(nNameLen));
    c.Skip(static_cast<size_t>(nNameLen));
    nAddr = UNDEF_ADDR;
    if (nLinkType == 0) nAddr = c.Addr(m_nSizeofAddr);
    if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated link message");
    return Status::OK;
}

Status File::LookupCompact(const std::vector<Message> &aoGroupMsgs, const std::string &osName,
                           uint64_t &nChildAddr)
{
    bool bAnyLink = false;
    bool bLinkInfo = false;
    for (const Message &oMsg : aoGroupMsgs) {
        if (oMsg.nType == MSG_LINK_INFO) bLinkInfo = true;
        if (oMsg.nType != MSG_LINK) continue;
        bAnyLink = true;
        if (oMsg.nFlags & MSG_FLAG_SHARED) return Fail(Status::ERR_MESSAGE, "shared link message");
        std::string osLink;
        int nLinkType = 0;
        uint64_t nAddr = UNDEF_ADDR;
        Status eStatus = DecodeLinkMessage(oMsg.abyData.data(), oMsg.abyData.size(), osLink, nLinkType, nAddr);
        if (eStatus != Status::OK) return eStatus;
        if (osLink == osName) {
            if (nLinkType != 0) return Fail(Status::ERR_GROUP, "non-hard link '" + osName + "'");
            nChildAddr = nAddr;
            return Status::OK;
        }
    }
    if (!bAnyLink && !bLinkInfo) return Fail(Status::ERR_GROUP, "object is not a group");
    return Fail(Status::ERR_NOT_FOUND, "'" + osName + "' not found");
}

Status File::ReadFractalHeapHeader(uint64_t nAddr, FractalHeap &oHeap)
{
    const size_t nFixed = 4 + 1 + 2 + 2 + 1 + 4 + m_nSizeofLen + m_nSizeofAddr + m_nSizeofLen + m_nSizeofAddr +
                          8 * m_nSizeofLen + 2 + 2 * m_nSizeofLen + 2 + 2 + m_nSizeofAddr + 2;
    std::vector<uint8_t> aby;
    READ_OR_FAIL(m_nBaseAddr + nAddr, nFixed, aby);
    Cursor c(aby);
    if (!c.Sig("FRHP")) return Fail(Status::ERR_SIGNATURE, "missing FRHP signature");
    if (c.U8() != 0) return Fail(Status::ERR_GROUP, "unsupported fractal heap version");
    oHeap.nIdLen = c.U16();
    oHeap.nFilterLen = c.U16();
    oHeap.nFlags = c.U8();
    oHeap.nMaxManSize = c.U32();
    c.U(m_nSizeofLen);        // next huge object id
    c.Addr(m_nSizeofAddr);    // huge objects v2 B-tree
    c.U(m_nSizeofLen);        // free space in managed blocks
    c.Addr(m_nSizeofAddr);    // free-space manager
    for (int i = 0; i < 8; i++) c.U(m_nSizeofLen);
    oHeap.nTableWidth = c.U16();
    oHeap.nStartBlockSize = c.U(m_nSizeofLen);
    oHeap.nMaxDirectBlockSize = c.U(m_nSizeofLen);
    oHeap.nMaxHeapBits = c.U16();
    c.U16();                  // starting # rows in root indirect block
    oHeap.nRootAddr = c.Addr(m_nSizeofAddr);
    oHeap.nRootRows = c.U16();
    if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated fractal heap header");
    if (oHeap.nFilterLen != 0) return Fail(Status::ERR_GROUP, "filtered fractal heap");
    if (oHeap.nTableWidth == 0 || oHeap.nStartBlockSize == 0 || oHeap.nMaxDirectBlockSize < oHeap.nStartBlockSize)
        return Fail(Status::ERR_CORRUPT, "bad fractal heap doubling table");
    oHeap.nHeapOffSize = (oHeap.nMaxHeapBits + 7) / 8;
    const unsigned nMaxDirOffSize = (Log2Gen(oHeap.nMaxDirectBlockSize) + 7) / 8;
    oHeap.nHeapLenSize = std::min(nMaxDirOffSize, LimitEncSize(oHeap.nMaxManSize));
    oHeap.nMaxDirectRows = Log2Gen(oHeap.nMaxDirectBlockSize) - Log2Gen(oHeap.nStartBlockSize) + 2;
    return Status::OK;
}

Status File::ReadFractalHeapObject(const FractalHeap &oHeap, const uint8_t *pabyId,
                                   std::vector<uint8_t> &abyObj)
{
    const uint8_t nIdFlags = pabyId[0];
    if ((nIdFlags >> 6) != 0) return Fail(Status::ERR_GROUP, "unsupported heap ID version");
    const int nIdType = (nIdFlags >> 4) & 0x03;
    if (nIdType == 2) {
        // Tiny object stored inline in the heap ID
        const bool bExtended = oHeap.nIdLen - 1 > 16;
        size_t nLen, nStart;
        if (bExtended) { nLen = (((nIdFlags & 0x0F) << 8) | pabyId[1]) + 1; nStart = 2; }
        else { nLen = (nIdFlags & 0x0F) + 1; nStart = 1; }
        if (nStart + nLen > oHeap.nIdLen) return Fail(Status::ERR_CORRUPT, "bad tiny heap ID");
        abyObj.assign(pabyId + nStart, pabyId + nStart + nLen);
        return Status::OK;
    }
    if (nIdType != 0) return Fail(Status::ERR_GROUP, "huge fractal heap objects not supported");

    Cursor c(pabyId + 1, oHeap.nIdLen - 1);
    uint64_t nOff = c.U(oHeap.nHeapOffSize);
    const uint64_t nLen = c.U(oHeap.nHeapLenSize);
    if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "bad managed heap ID");
    if (oHeap.nRootAddr == UNDEF_ADDR) return Fail(Status::ERR_CORRUPT, "empty fractal heap");

    const uint64_t W = oHeap.nTableWidth;
    const uint64_t S = oHeap.nStartBlockSize;
    auto RowBlockSize = [&](unsigned r) -> uint64_t { return r == 0 ? S : S << (r - 1); };
    auto RowStart = [&](unsigned r) -> uint64_t { return r == 0 ? 0 : W * S << (r - 1); };

    uint64_t nBlockAddr = oHeap.nRootAddr;
    uint64_t nBlockHeapOff = 0;   // heap offset where the current block starts
    unsigned nRows = oHeap.nRootRows;
    std::vector<uint8_t> aby;
    for (int nGuard = 0; nRows > 0; nGuard++) {
        if (nGuard > MAX_DEPTH) return Fail(Status::ERR_LIMIT, "fractal heap too deep");
        const uint64_t nRel = nOff - nBlockHeapOff;
        unsigned r = 0;
        while (r + 1 < nRows && RowStart(r + 1) <= nRel) r++;
        const uint64_t nCol = (nRel - RowStart(r)) / RowBlockSize(r);
        if (nCol >= W) return Fail(Status::ERR_CORRUPT, "heap offset outside indirect block");
        const uint64_t nEntry = r * W + nCol;
        const size_t nIBHdr = 5 + m_nSizeofAddr + oHeap.nHeapOffSize;
        std::vector<uint8_t> abySig;
        READ_OR_FAIL(m_nBaseAddr + nBlockAddr, 4, abySig);
        if (memcmp(abySig.data(), "FHIB", 4) != 0) return Fail(Status::ERR_SIGNATURE, "missing FHIB signature");
        READ_OR_FAIL(m_nBaseAddr + nBlockAddr + nIBHdr + nEntry * m_nSizeofAddr, m_nSizeofAddr, aby);
        Cursor e(aby);
        const uint64_t nChild = e.Addr(m_nSizeofAddr);
        if (nChild == UNDEF_ADDR) return Fail(Status::ERR_CORRUPT, "heap object in unallocated block");
        nBlockHeapOff += RowStart(r) + nCol * RowBlockSize(r);
        nBlockAddr = nChild;
        if (r < oHeap.nMaxDirectRows) { nRows = 0; break; }
        nRows = Log2Gen(RowBlockSize(r)) - Log2Gen(S * W) + 1;
    }

    // Direct block
    std::vector<uint8_t> abyHdr;
    READ_OR_FAIL(m_nBaseAddr + nBlockAddr, 4, abyHdr);
    if (memcmp(abyHdr.data(), "FHDB", 4) != 0) return Fail(Status::ERR_SIGNATURE, "missing FHDB signature");
    READ_OR_FAIL(m_nBaseAddr + nBlockAddr + (nOff - nBlockHeapOff), static_cast<size_t>(nLen), abyObj);
    return Status::OK;
}

Status File::ReadBTree2Header(uint64_t nAddr, BTree2Header &oHdr)
{
    const size_t nLen = 4 + 1 + 1 + 4 + 2 + 2 + 1 + 1 + m_nSizeofAddr + 2 + m_nSizeofLen + 4;
    std::vector<uint8_t> aby;
    READ_OR_FAIL(m_nBaseAddr + nAddr, nLen, aby);
    Cursor c(aby);
    if (!c.Sig("BTHD")) return Fail(Status::ERR_SIGNATURE, "missing BTHD signature");
    if (c.U8() != 0) return Fail(Status::ERR_INDEX, "unsupported v2 B-tree version");
    oHdr.nType = c.U8();
    oHdr.nNodeSize = c.U32();
    oHdr.nRecordSize = c.U16();
    oHdr.nDepth = c.U16();
    c.U8();
    c.U8();
    oHdr.nRootAddr = c.Addr(m_nSizeofAddr);
    oHdr.nRootRecords = c.U16();
    oHdr.nTotalRecords = c.U(m_nSizeofLen);
    const size_t nChecked = c.Pos();
    const uint32_t nStored = c.U32();
    if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated BTHD");
    if (Lookup3(aby.data(), nChecked, 0) != nStored) return Fail(Status::ERR_CORRUPT, "BTHD checksum mismatch");
    if (oHdr.nRecordSize == 0 || oHdr.nNodeSize <= 10u + oHdr.nRecordSize || oHdr.nDepth > MAX_DEPTH)
        return Fail(Status::ERR_CORRUPT, "bad v2 B-tree parameters");

    // Record-count field widths (H5B2__hdr_init)
    const uint64_t nLeafMax = (oHdr.nNodeSize - 10) / oHdr.nRecordSize;
    oHdr.nMaxNrecSize = LimitEncSize(nLeafMax);
    std::vector<uint64_t> anCumMax(oHdr.nDepth + 1);
    oHdr.anCumMaxNrecSize.assign(oHdr.nDepth + 1, 0);
    anCumMax[0] = nLeafMax;
    for (unsigned d = 1; d <= oHdr.nDepth; d++) {
        const uint64_t nPtr = m_nSizeofAddr + oHdr.nMaxNrecSize + (d > 1 ? oHdr.anCumMaxNrecSize[d - 1] : 0);
        if (oHdr.nNodeSize <= 10 + nPtr) return Fail(Status::ERR_CORRUPT, "bad v2 B-tree node size");
        const uint64_t nIntMax = (oHdr.nNodeSize - (10 + nPtr)) / (oHdr.nRecordSize + nPtr);
        anCumMax[d] = (nIntMax + 1) * anCumMax[d - 1] + nIntMax;
        oHdr.anCumMaxNrecSize[d] = LimitEncSize(anCumMax[d]);
    }
    return Status::OK;
}

Status File::CollectBTree2Records(const BTree2Header &oHdr, uint64_t nNodeAddr, unsigned nDepth,
                                  unsigned nRecords, std::vector<std::vector<uint8_t>> &aoRecords,
                                  int nDepthGuard)
{
    if (nDepthGuard > MAX_DEPTH) return Fail(Status::ERR_LIMIT, "v2 B-tree too deep");
    if (nNodeAddr == UNDEF_ADDR) return Status::OK;
    std::vector<uint8_t> aby;
    READ_OR_FAIL(m_nBaseAddr + nNodeAddr, oHdr.nNodeSize, aby);
    Cursor c(aby);
    const char *pszSig = nDepth == 0 ? "BTLF" : "BTIN";
    if (!c.Sig(pszSig)) return Fail(Status::ERR_SIGNATURE, std::string("missing ") + pszSig + " signature");
    if (c.U8() != 0) return Fail(Status::ERR_INDEX, "unsupported v2 B-tree node version");
    if (c.U8() != oHdr.nType) return Fail(Status::ERR_CORRUPT, "v2 B-tree node type mismatch");
    std::vector<const uint8_t *> apRec(nRecords);
    for (unsigned i = 0; i < nRecords; i++) { apRec[i] = c.Ptr(); c.Skip(oHdr.nRecordSize); }
    if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated v2 B-tree node");
    if (nDepth == 0) {
        for (unsigned i = 0; i < nRecords; i++) aoRecords.emplace_back(apRec[i], apRec[i] + oHdr.nRecordSize);
        return Status::OK;
    }
    struct Child { uint64_t nAddr; unsigned nRecs; };
    std::vector<Child> aoChildren(nRecords + 1);
    for (unsigned i = 0; i <= nRecords; i++) {
        aoChildren[i].nAddr = c.Addr(m_nSizeofAddr);
        aoChildren[i].nRecs = static_cast<unsigned>(c.U(oHdr.nMaxNrecSize));
        if (nDepth > 1) c.U(oHdr.anCumMaxNrecSize[nDepth - 1]);
    }
    if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated v2 B-tree internal node");
    for (unsigned i = 0; i <= nRecords; i++) {
        Status eStatus = CollectBTree2Records(oHdr, aoChildren[i].nAddr, nDepth - 1, aoChildren[i].nRecs, aoRecords, nDepthGuard + 1);
        if (eStatus != Status::OK) return eStatus;
        if (i < nRecords) aoRecords.emplace_back(apRec[i], apRec[i] + oHdr.nRecordSize);
    }
    return Status::OK;
}

Status File::LookupDense(uint64_t nHeapAddr, uint64_t nNameIndexAddr, const std::string &osName,
                         uint64_t &nChildAddr)
{
    FractalHeap oHeap;
    Status eStatus = ReadFractalHeapHeader(nHeapAddr, oHeap);
    if (eStatus != Status::OK) return eStatus;
    BTree2Header oBT;
    eStatus = ReadBTree2Header(nNameIndexAddr, oBT);
    if (eStatus != Status::OK) return eStatus;
    if (oBT.nType != 5) return Fail(Status::ERR_GROUP, "unexpected name index B-tree type");
    if (oBT.nRecordSize != 4 + oHeap.nIdLen) return Fail(Status::ERR_CORRUPT, "name index record size mismatch");

    const uint32_t nHash = Lookup3(reinterpret_cast<const uint8_t *>(osName.data()), osName.size(), 0);

    // Descend by hash; keep every record with a matching hash (collisions are rare).
    std::vector<std::vector<uint8_t>> aoMatches;
    struct Frame { uint64_t nAddr; unsigned nDepth; unsigned nRecs; };
    std::vector<Frame> aoStack = {{oBT.nRootAddr, oBT.nDepth, oBT.nRootRecords}};
    std::vector<uint8_t> aby;
    size_t nVisited = 0;
    while (!aoStack.empty()) {
        const Frame f = aoStack.back();
        aoStack.pop_back();
        if (f.nAddr == UNDEF_ADDR || ++nVisited > 100000) continue;
        READ_OR_FAIL(m_nBaseAddr + f.nAddr, oBT.nNodeSize, aby);
        Cursor c(aby);
        if (!c.Sig(f.nDepth == 0 ? "BTLF" : "BTIN")) return Fail(Status::ERR_SIGNATURE, "bad v2 B-tree node signature");
        c.U8();
        c.U8();
        std::vector<uint32_t> anHash(f.nRecs);
        for (unsigned i = 0; i < f.nRecs; i++) {
            const uint8_t *p = c.Ptr();
            c.Skip(oBT.nRecordSize);
            if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated name index node");
            anHash[i] = static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8) |
                        (static_cast<uint32_t>(p[2]) << 16) | (static_cast<uint32_t>(p[3]) << 24);
            if (anHash[i] == nHash) aoMatches.emplace_back(p + 4, p + oBT.nRecordSize);
        }
        if (f.nDepth == 0) continue;
        for (unsigned i = 0; i <= f.nRecs; i++) {
            const uint64_t nAddr = c.Addr(m_nSizeofAddr);
            const unsigned nRecs = static_cast<unsigned>(c.U(oBT.nMaxNrecSize));
            if (f.nDepth > 1) c.U(oBT.anCumMaxNrecSize[f.nDepth - 1]);
            if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated name index node");
            const bool bLeftOk = (i == 0) || anHash[i - 1] <= nHash;
            const bool bRightOk = (i == f.nRecs) || anHash[i] >= nHash;
            if (bLeftOk && bRightOk) aoStack.push_back({nAddr, f.nDepth - 1u, nRecs});
        }
    }

    for (const std::vector<uint8_t> &abyId : aoMatches) {
        std::vector<uint8_t> abyLink;
        eStatus = ReadFractalHeapObject(oHeap, abyId.data(), abyLink);
        if (eStatus != Status::OK) return eStatus;
        std::string osLink;
        int nLinkType = 0;
        uint64_t nAddr = UNDEF_ADDR;
        eStatus = DecodeLinkMessage(abyLink.data(), abyLink.size(), osLink, nLinkType, nAddr);
        if (eStatus != Status::OK) return eStatus;
        if (osLink == osName) {
            if (nLinkType != 0) return Fail(Status::ERR_GROUP, "non-hard link '" + osName + "'");
            nChildAddr = nAddr;
            return Status::OK;
        }
    }
    return Fail(Status::ERR_NOT_FOUND, "'" + osName + "' not found");
}

/************************************************************************/
/*                         Datasets                                     */
/************************************************************************/

Status File::OpenDataset(const std::string &osPath, DatasetInfo &oInfo)
{
    oInfo = DatasetInfo();
    oInfo.osPath = osPath;
    uint64_t nAddr = m_nRootAddr;
    std::vector<Message> aoMsgs;
    size_t nPos = 0;
    int nVersion = 0;
    Status eStatus;
    while (nPos < osPath.size()) {
        const size_t nSlash = osPath.find('/', nPos);
        const std::string osComp = osPath.substr(nPos, nSlash == std::string::npos ? std::string::npos : nSlash - nPos);
        nPos = (nSlash == std::string::npos) ? osPath.size() : nSlash + 1;
        if (osComp.empty() || osComp == ".") continue;
        eStatus = ReadObjectHeader(nAddr, aoMsgs, &nVersion);
        if (eStatus != Status::OK) return eStatus;
        eStatus = LookupChild(aoMsgs, osComp, nAddr);
        if (eStatus != Status::OK) return eStatus;
    }
    eStatus = ReadObjectHeader(nAddr, aoMsgs, &nVersion);
    if (eStatus != Status::OK) return eStatus;
    oInfo.nObjectHeaderAddr = nAddr;
    oInfo.nObjectHeaderVersion = nVersion;
    return DecodeDataset(aoMsgs, oInfo);
}

Status File::DecodeDataset(const std::vector<Message> &aoMsgs, DatasetInfo &oInfo)
{
    const Message *poSpace = nullptr, *poType = nullptr, *poLayout = nullptr, *poPipeline = nullptr;
    for (const Message &oMsg : aoMsgs) {
        switch (oMsg.nType) {
            case MSG_DATASPACE: poSpace = &oMsg; break;
            case MSG_DATATYPE: poType = &oMsg; break;
            case MSG_LAYOUT: poLayout = &oMsg; break;
            case MSG_FILTER_PIPELINE: poPipeline = &oMsg; break;
            case MSG_EXTERNAL_FILES: return Fail(Status::ERR_LAYOUT, "external data files");
            default: break;
        }
    }
    if (!poLayout) return Fail(Status::ERR_NOT_DATASET, "object has no layout message");
    if (!poSpace || !poType) return Fail(Status::ERR_MESSAGE, "dataset lacks dataspace/datatype message");
    for (const Message *p : {poSpace, poType, poLayout, poPipeline})
        if (p && (p->nFlags & MSG_FLAG_SHARED)) return Fail(Status::ERR_MESSAGE, "shared dataset message");

    // Dataspace
    {
        Cursor c(poSpace->abyData);
        const uint8_t nVersion = c.U8();
        const uint8_t nRank = c.U8();
        const uint8_t nFlags = c.U8();
        if (nVersion == 1) { c.Skip(5); }
        else if (nVersion == 2) {
            const uint8_t nSpaceType = c.U8();
            if (nSpaceType != 1) return Fail(Status::ERR_NOT_CHUNKED, "scalar or null dataspace");
        }
        else return Fail(Status::ERR_MESSAGE, "unsupported dataspace version");
        oInfo.nRank = nRank;
        oInfo.anDims.resize(nRank);
        for (auto &d : oInfo.anDims) d = c.U(m_nSizeofLen);
        oInfo.anMaxDims = oInfo.anDims;
        if (nFlags & 0x01)
            for (auto &d : oInfo.anMaxDims) d = c.U(m_nSizeofLen);
        if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated dataspace message");
        for (auto &d : oInfo.anMaxDims)
            if (d == ((m_nSizeofLen >= 8) ? UNLIMITED : ((static_cast<uint64_t>(1) << (8 * m_nSizeofLen)) - 1))) d = UNLIMITED;
    }

    // Datatype
    {
        Cursor c(poType->abyData);
        const uint8_t nClassVer = c.U8();
        c.Skip(3);
        oInfo.nTypeClass = nClassVer & 0x0F;
        oInfo.nElementSize = c.U32();
        if (!c.Ok() || (nClassVer >> 4) < 1 || (nClassVer >> 4) > 5)
            return Fail(Status::ERR_MESSAGE, "unsupported datatype message");
    }

    // Filter pipeline
    if (poPipeline) {
        Cursor c(poPipeline->abyData);
        const uint8_t nVersion = c.U8();
        const uint8_t nFilters = c.U8();
        if (nVersion == 1) c.Skip(6);
        else if (nVersion != 2) return Fail(Status::ERR_MESSAGE, "unsupported filter pipeline version");
        for (unsigned i = 0; i < nFilters; i++) {
            const uint16_t nId = c.U16();
            uint16_t nNameLen = 0;
            if (nVersion == 1 || nId >= 256) nNameLen = c.U16();
            c.U16();  // flags
            const uint16_t nValues = c.U16();
            if (nVersion == 1) c.Skip((nNameLen + 7) & ~7u);
            else c.Skip(nNameLen);
            c.Skip(4u * nValues);
            if (nVersion == 1 && (nValues & 1)) c.Skip(4);
            oInfo.anFilterIds.push_back(nId);
        }
        if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated filter pipeline message");
    }

    // Layout
    Cursor c(poLayout->abyData);
    oInfo.nLayoutVersion = c.U8();
    if (oInfo.nLayoutVersion != 3 && oInfo.nLayoutVersion != 4)
        return Fail(Status::ERR_LAYOUT, "unsupported layout version " + std::to_string(oInfo.nLayoutVersion));
    oInfo.nLayoutClass = c.U8();
    if (oInfo.nLayoutClass == 0 || oInfo.nLayoutClass == 1) return Fail(Status::ERR_NOT_CHUNKED, "dataset is not chunked");
    if (oInfo.nLayoutClass != 2) return Fail(Status::ERR_LAYOUT, "unsupported layout class " + std::to_string(oInfo.nLayoutClass));

    std::vector<uint64_t> anLayoutDims;
    if (oInfo.nLayoutVersion == 3) {
        const uint8_t nDims = c.U8();
        oInfo.nIndexAddr = c.Addr(m_nSizeofAddr);
        anLayoutDims.resize(nDims);
        for (auto &d : anLayoutDims) d = c.U32();
        oInfo.eIndex = ChunkIndex::BTREE_V1;
    }
    else {
        oInfo.nLayoutFlags = c.U8();
        const uint8_t nDims = c.U8();
        const uint8_t nEncLen = c.U8();
        if (nEncLen < 1 || nEncLen > 8) return Fail(Status::ERR_LAYOUT, "bad layout dimension encoding");
        anLayoutDims.resize(nDims);
        for (auto &d : anLayoutDims) d = c.U(nEncLen);
        const uint8_t nIndexType = c.U8();
        switch (nIndexType) {
            case 1:
                oInfo.eIndex = ChunkIndex::SINGLE_CHUNK;
                if (oInfo.nLayoutFlags & 0x02) {
                    oInfo.nSingleFilteredSize = c.U(m_nSizeofLen);
                    oInfo.nSingleFilterMask = c.U32();
                }
                break;
            case 2: oInfo.eIndex = ChunkIndex::IMPLICIT; break;
            case 3:
                oInfo.eIndex = ChunkIndex::FIXED_ARRAY;
                oInfo.nFarrayPageBits = c.U8();
                break;
            case 4:
                oInfo.eIndex = ChunkIndex::EXTENSIBLE_ARRAY;
                oInfo.nEarrayMaxBits = c.U8();
                oInfo.nEarrayIdxBlkElmts = c.U8();
                oInfo.nEarrayMinPtrs = c.U8();
                oInfo.nEarrayMinElmts = c.U8();
                oInfo.nEarrayPageBits = c.U8();
                break;
            case 5:
                oInfo.eIndex = ChunkIndex::BTREE_V2;
                c.U32();
                c.U8();
                c.U8();
                break;
            default:
                return Fail(Status::ERR_INDEX, "unsupported chunk index type " + std::to_string(nIndexType));
        }
        oInfo.nIndexAddr = c.Addr(m_nSizeofAddr);
    }
    if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated layout message");
    if (anLayoutDims.size() != static_cast<size_t>(oInfo.nRank) + 1)
        return Fail(Status::ERR_LAYOUT, "layout dimensionality does not match dataspace rank");
    if (anLayoutDims.back() != oInfo.nElementSize)
        return Fail(Status::ERR_LAYOUT, "layout element size does not match datatype size");
    oInfo.anChunkDims.assign(anLayoutDims.begin(), anLayoutDims.end() - 1);
    for (uint64_t d : oInfo.anChunkDims)
        if (d == 0) return Fail(Status::ERR_CORRUPT, "zero chunk dimension");
    return Status::OK;
}

/************************************************************************/
/*                         Chunk indexes                                */
/************************************************************************/

Status File::GetChunks(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks)
{
    aoChunks.clear();
    if (oInfo.nLayoutClass != 2) return Fail(Status::ERR_NOT_CHUNKED, "dataset is not chunked");
    if (oInfo.nIndexAddr == UNDEF_ADDR) return Status::OK;  // no storage allocated yet
    Status eStatus = Status::OK;
    switch (oInfo.eIndex) {
        case ChunkIndex::BTREE_V1: eStatus = IterBTreeV1(oInfo, oInfo.nIndexAddr, -1, aoChunks, 0); break;
        case ChunkIndex::SINGLE_CHUNK: eStatus = IterSingle(oInfo, aoChunks); break;
        case ChunkIndex::IMPLICIT: eStatus = IterImplicit(oInfo, aoChunks); break;
        case ChunkIndex::FIXED_ARRAY: eStatus = IterFixedArray(oInfo, aoChunks); break;
        case ChunkIndex::EXTENSIBLE_ARRAY: eStatus = IterExtensibleArray(oInfo, aoChunks); break;
        case ChunkIndex::BTREE_V2: eStatus = IterBTreeV2(oInfo, aoChunks); break;
        case ChunkIndex::NONE: return Fail(Status::ERR_INDEX, "no chunk index");
    }
    if (eStatus != Status::OK) return eStatus;
    // Index entries are relative to the superblock base; callers get physical file offsets.
    for (ChunkRecord &oRec : aoChunks) {
        if (oRec.nAddr > UINT64_MAX - m_nBaseAddr) return Fail(Status::ERR_CORRUPT, "chunk address overflow");
        oRec.nAddr += m_nBaseAddr;
    }
    return Status::OK;
}

Status File::IterBTreeV1(const DatasetInfo &oInfo, uint64_t nAddr, int nExpectedLevel,
                         std::vector<ChunkRecord> &aoChunks, int nDepthGuard)
{
    if (nDepthGuard > MAX_DEPTH) return Fail(Status::ERR_LIMIT, "chunk B-tree too deep");
    const size_t nNodeHdr = 8 + 2 * m_nSizeofAddr;
    std::vector<uint8_t> aby;
    READ_OR_FAIL(m_nBaseAddr + nAddr, nNodeHdr, aby);
    Cursor h(aby);
    if (!h.Sig("TREE")) return Fail(Status::ERR_SIGNATURE, "missing TREE signature");
    const uint8_t nType = h.U8();
    const int nLevel = h.U8();
    const uint16_t nEntries = h.U16();
    if (nType != 1) return Fail(Status::ERR_CORRUPT, "expected raw-data chunk B-tree node");
    if (nExpectedLevel >= 0 && nLevel != nExpectedLevel) return Fail(Status::ERR_CORRUPT, "chunk B-tree level mismatch");

    const unsigned nKeyDims = static_cast<unsigned>(oInfo.nRank) + 1;
    const size_t nKeySize = 8 + 8 * nKeyDims;
    READ_OR_FAIL(m_nBaseAddr + nAddr + nNodeHdr, nEntries * (nKeySize + m_nSizeofAddr) + nKeySize, aby);
    Cursor c(aby);
    for (unsigned i = 0; i < nEntries; i++) {
        ChunkRecord oRec;
        oRec.nSize = c.U32();
        oRec.nFilterMask = c.U32();
        oRec.anOffset.resize(oInfo.nRank);
        for (auto &o : oRec.anOffset) o = c.U64();
        c.U64();  // element-size dimension, always 0
        const uint64_t nChild = c.Addr(m_nSizeofAddr);
        if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated chunk B-tree node");
        if (nLevel > 0) {
            Status eStatus = IterBTreeV1(oInfo, nChild, nLevel - 1, aoChunks, nDepthGuard + 1);
            if (eStatus != Status::OK) return eStatus;
        }
        else {
            oRec.nAddr = nChild;
            if (nChild != UNDEF_ADDR) aoChunks.push_back(std::move(oRec));
        }
    }
    return Status::OK;
}

Status File::IterSingle(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks)
{
    ChunkRecord oRec;
    oRec.anOffset.assign(oInfo.nRank, 0);
    oRec.nAddr = oInfo.nIndexAddr;
    if (oInfo.nLayoutFlags & 0x02) {
        oRec.nSize = oInfo.nSingleFilteredSize;
        oRec.nFilterMask = oInfo.nSingleFilterMask;
    }
    else {
        oRec.nSize = oInfo.ChunkBytes();
    }
    aoChunks.push_back(std::move(oRec));
    return Status::OK;
}

static std::vector<uint64_t> ChunkGrid(const std::vector<uint64_t> &anDims, const std::vector<uint64_t> &anChunk)
{
    std::vector<uint64_t> an(anDims.size());
    for (size_t i = 0; i < an.size(); i++)
        an[i] = anDims[i] == UNLIMITED ? UNLIMITED : (anDims[i] + anChunk[i] - 1) / anChunk[i];
    return an;
}

Status File::IterImplicit(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks)
{
    const std::vector<uint64_t> anCur = ChunkGrid(oInfo.anDims, oInfo.anChunkDims);
    const std::vector<uint64_t> anMax = ChunkGrid(oInfo.anMaxDims, oInfo.anChunkDims);
    const int nRank = oInfo.nRank;
    uint64_t nTotal = 1;
    for (int i = 0; i < nRank; i++) {
        if (anMax[i] == UNLIMITED) return Fail(Status::ERR_INDEX, "implicit index with unlimited dimension");
        nTotal *= anCur[i];
    }
    std::vector<uint64_t> anDown(nRank, 1);
    for (int i = nRank - 2; i >= 0; i--) anDown[i] = anDown[i + 1] * anMax[i + 1];
    const uint64_t nChunkBytes = oInfo.ChunkBytes();
    std::vector<uint64_t> anScaled(nRank, 0);
    for (uint64_t n = 0; n < nTotal; n++) {
        uint64_t nIdx = 0;
        for (int i = 0; i < nRank; i++) nIdx += anScaled[i] * anDown[i];
        ChunkRecord oRec;
        oRec.anOffset.resize(nRank);
        for (int i = 0; i < nRank; i++) oRec.anOffset[i] = anScaled[i] * oInfo.anChunkDims[i];
        oRec.nAddr = oInfo.nIndexAddr + nIdx * nChunkBytes;
        oRec.nSize = nChunkBytes;
        aoChunks.push_back(std::move(oRec));
        for (int d = nRank - 1; d >= 0; d--) {
            if (++anScaled[d] < anCur[d]) break;
            anScaled[d] = 0;
        }
    }
    return Status::OK;
}

// Decodes one array element (chunk address, optional size + filter mask).
static bool DecodeArrayElement(Cursor &c, unsigned nSizeofAddr, bool bFiltered, unsigned nSizeLen,
                               uint64_t nChunkBytes, ChunkRecord &oRec)
{
    oRec.nAddr = c.Addr(nSizeofAddr);
    if (bFiltered) {
        oRec.nSize = c.U(nSizeLen);
        oRec.nFilterMask = c.U32();
    }
    else {
        oRec.nSize = nChunkBytes;
        oRec.nFilterMask = 0;
    }
    return c.Ok();
}

static bool BitGet(const uint8_t *pabyBitmap, uint64_t nBit)
{
    return (pabyBitmap[nBit / 8] & (0x80 >> (nBit % 8))) != 0;
}

Status File::IterFixedArray(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks)
{
    const size_t nHdrLen = 4 + 1 + 1 + 1 + 1 + m_nSizeofLen + m_nSizeofAddr + 4;
    std::vector<uint8_t> aby;
    READ_OR_FAIL(m_nBaseAddr + oInfo.nIndexAddr, nHdrLen, aby);
    Cursor h(aby);
    if (!h.Sig("FAHD")) return Fail(Status::ERR_SIGNATURE, "missing FAHD signature");
    if (h.U8() != 0) return Fail(Status::ERR_INDEX, "unsupported fixed array version");
    const uint8_t nClient = h.U8();
    const uint8_t nEntrySize = h.U8();
    const uint8_t nPageBits = h.U8();
    const uint64_t nElmts = h.U(m_nSizeofLen);
    const uint64_t nDblkAddr = h.Addr(m_nSizeofAddr);
    const size_t nChecked = h.Pos();
    const uint32_t nStored = h.U32();
    if (!h.Ok()) return Fail(Status::ERR_CORRUPT, "truncated FAHD");
    if (Lookup3(aby.data(), nChecked, 0) != nStored) return Fail(Status::ERR_CORRUPT, "FAHD checksum mismatch");
    if (nClient > 1) return Fail(Status::ERR_INDEX, "unsupported fixed array client");
    const bool bFiltered = nClient == 1;
    if (bFiltered != !oInfo.anFilterIds.empty())
        return Fail(Status::ERR_CORRUPT, "fixed array client does not match filter pipeline");
    if (bFiltered && nEntrySize <= m_nSizeofAddr + 4) return Fail(Status::ERR_CORRUPT, "bad fixed array entry size");
    if (!bFiltered && nEntrySize != m_nSizeofAddr) return Fail(Status::ERR_CORRUPT, "bad fixed array entry size");
    const unsigned nSizeLen = bFiltered ? nEntrySize - m_nSizeofAddr - 4 : 0;

    const std::vector<uint64_t> anMax = ChunkGrid(oInfo.anMaxDims, oInfo.anChunkDims);
    uint64_t nExpected = 1;
    for (uint64_t n : anMax) {
        if (n == UNLIMITED) return Fail(Status::ERR_INDEX, "fixed array with unlimited dimension");
        nExpected *= n;
    }
    if (nExpected != nElmts) return Fail(Status::ERR_CORRUPT, "fixed array size does not match chunk grid");
    if (nDblkAddr == UNDEF_ADDR) return Status::OK;

    const uint64_t nPageElmts = static_cast<uint64_t>(1) << nPageBits;
    const bool bPaged = nElmts > nPageElmts;
    const uint64_t nPages = bPaged ? (nElmts + nPageElmts - 1) / nPageElmts : 0;
    const size_t nBitmapLen = static_cast<size_t>((nPages + 7) / 8);
    const size_t nDblkHdr = 4 + 1 + 1 + m_nSizeofAddr;

    const uint64_t nChunkBytes = oInfo.ChunkBytes();
    const int nRank = oInfo.nRank;
    std::vector<uint64_t> anScaled(nRank, 0);
    auto Emit = [&](Cursor &c) -> bool {
        ChunkRecord oRec;
        if (!DecodeArrayElement(c, m_nSizeofAddr, bFiltered, nSizeLen, nChunkBytes, oRec)) return false;
        if (oRec.nAddr != UNDEF_ADDR) {
            oRec.anOffset.resize(nRank);
            for (int i = 0; i < nRank; i++) oRec.anOffset[i] = anScaled[i] * oInfo.anChunkDims[i];
            aoChunks.push_back(std::move(oRec));
        }
        for (int d = nRank - 1; d >= 0; d--) {
            if (++anScaled[d] < anMax[d]) break;
            anScaled[d] = 0;
        }
        return true;
    };
    auto SkipElmts = [&](uint64_t n) {
        for (uint64_t k = 0; k < n; k++)
            for (int d = nRank - 1; d >= 0; d--) {
                if (++anScaled[d] < anMax[d]) break;
                anScaled[d] = 0;
            }
    };

    if (!bPaged) {
        READ_OR_FAIL(m_nBaseAddr + nDblkAddr, static_cast<size_t>(nDblkHdr + nElmts * nEntrySize + 4), aby);
        Cursor c(aby);
        if (!c.Sig("FADB")) return Fail(Status::ERR_SIGNATURE, "missing FADB signature");
        if (Lookup3(aby.data(), aby.size() - 4, 0) != Cursor(aby.data() + aby.size() - 4, 4).U32())
            return Fail(Status::ERR_CORRUPT, "FADB checksum mismatch");
        c.Skip(nDblkHdr - 4);
        for (uint64_t i = 0; i < nElmts; i++)
            if (!Emit(c)) return Fail(Status::ERR_CORRUPT, "truncated fixed array data block");
        return Status::OK;
    }

    READ_OR_FAIL(m_nBaseAddr + nDblkAddr, nDblkHdr + nBitmapLen + 4, aby);
    if (memcmp(aby.data(), "FADB", 4) != 0) return Fail(Status::ERR_SIGNATURE, "missing FADB signature");
    if (Lookup3(aby.data(), aby.size() - 4, 0) != Cursor(aby.data() + aby.size() - 4, 4).U32())
        return Fail(Status::ERR_CORRUPT, "FADB checksum mismatch");
    const std::vector<uint8_t> abyBitmap(aby.begin() + nDblkHdr, aby.begin() + nDblkHdr + nBitmapLen);
    const uint64_t nPageBytes = nPageElmts * nEntrySize + 4;
    const uint64_t nFirstPage = nDblkAddr + nDblkHdr + nBitmapLen + 4;
    for (uint64_t p = 0; p < nPages; p++) {
        const uint64_t nThis = std::min(nPageElmts, nElmts - p * nPageElmts);
        if (!BitGet(abyBitmap.data(), p)) { SkipElmts(nThis); continue; }
        READ_OR_FAIL(m_nBaseAddr + nFirstPage + p * nPageBytes, static_cast<size_t>(nThis * nEntrySize + 4), aby);
        if (Lookup3(aby.data(), aby.size() - 4, 0) != Cursor(aby.data() + aby.size() - 4, 4).U32())
            return Fail(Status::ERR_CORRUPT, "fixed array page checksum mismatch");
        Cursor c(aby);
        for (uint64_t i = 0; i < nThis; i++)
            if (!Emit(c)) return Fail(Status::ERR_CORRUPT, "truncated fixed array page");
    }
    return Status::OK;
}

Status File::IterExtensibleArray(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks)
{
    const size_t nHdrLen = 4 + 1 + 1 + 6 + 6 * m_nSizeofLen + m_nSizeofAddr + 4;
    std::vector<uint8_t> aby;
    READ_OR_FAIL(m_nBaseAddr + oInfo.nIndexAddr, nHdrLen, aby);
    Cursor h(aby);
    if (!h.Sig("EAHD")) return Fail(Status::ERR_SIGNATURE, "missing EAHD signature");
    if (h.U8() != 0) return Fail(Status::ERR_INDEX, "unsupported extensible array version");
    const uint8_t nClient = h.U8();
    const uint8_t nElmtSize = h.U8();
    const uint8_t nMaxBits = h.U8();
    const uint8_t nIdxBlkElmts = h.U8();
    const uint8_t nDblkMinElmts = h.U8();
    const uint8_t nSblkMinPtrs = h.U8();
    const uint8_t nPageBits = h.U8();
    for (int i = 0; i < 4; i++) h.U(m_nSizeofLen);
    const uint64_t nMaxIdxSet = h.U(m_nSizeofLen);
    h.U(m_nSizeofLen);
    const uint64_t nIblkAddr = h.Addr(m_nSizeofAddr);
    const size_t nChecked = h.Pos();
    const uint32_t nStored = h.U32();
    if (!h.Ok()) return Fail(Status::ERR_CORRUPT, "truncated EAHD");
    if (Lookup3(aby.data(), nChecked, 0) != nStored) return Fail(Status::ERR_CORRUPT, "EAHD checksum mismatch");
    if (nClient > 1) return Fail(Status::ERR_INDEX, "unsupported extensible array client");
    const bool bFiltered = nClient == 1;
    if (bFiltered != !oInfo.anFilterIds.empty())
        return Fail(Status::ERR_CORRUPT, "extensible array client does not match filter pipeline");
    if ((bFiltered && nElmtSize <= m_nSizeofAddr + 4) || (!bFiltered && nElmtSize != m_nSizeofAddr))
        return Fail(Status::ERR_CORRUPT, "bad extensible array element size");
    if (nDblkMinElmts == 0 || (nDblkMinElmts & (nDblkMinElmts - 1)) || nSblkMinPtrs == 0 ||
        (nSblkMinPtrs & (nSblkMinPtrs - 1)) || nMaxBits > 64)
        return Fail(Status::ERR_CORRUPT, "bad extensible array parameters");
    const unsigned nSizeLen = bFiltered ? nElmtSize - m_nSizeofAddr - 4 : 0;
    if (nIblkAddr == UNDEF_ADDR || nMaxIdxSet == 0) return Status::OK;

    // Super block table (H5EA__hdr_init)
    struct SblkInfo { uint64_t nDblks; uint64_t nDblkElmts; uint64_t nStartIdx; uint64_t nStartDblk; };
    const unsigned nSblks = 1 + (nMaxBits - Log2Gen(nDblkMinElmts));
    std::vector<SblkInfo> aoSblk(nSblks);
    uint64_t nStartIdx = 0, nStartDblk = 0;
    for (unsigned u = 0; u < nSblks; u++) {
        aoSblk[u].nDblks = static_cast<uint64_t>(1) << (u / 2);
        aoSblk[u].nDblkElmts = (static_cast<uint64_t>(1) << ((u + 1) / 2)) * nDblkMinElmts;
        aoSblk[u].nStartIdx = nStartIdx;
        aoSblk[u].nStartDblk = nStartDblk;
        nStartIdx += aoSblk[u].nDblks * aoSblk[u].nDblkElmts;
        nStartDblk += aoSblk[u].nDblks;
    }
    const unsigned nIblkSblks = 2 * Log2Gen(nSblkMinPtrs);
    const size_t nIblkDblkAddrs = 2 * (static_cast<size_t>(nSblkMinPtrs) - 1);
    if (nIblkSblks > nSblks) return Fail(Status::ERR_CORRUPT, "bad extensible array super block count");
    const size_t nIblkSblkAddrs = nSblks - nIblkSblks;
    const size_t nArrOffSize = (nMaxBits + 7) / 8;
    const uint64_t nPageElmts = static_cast<uint64_t>(1) << nPageBits;

    // Index block
    const size_t nIblkLen = 4 + 1 + 1 + m_nSizeofAddr + nIdxBlkElmts * static_cast<size_t>(nElmtSize) +
                            (nIblkDblkAddrs + nIblkSblkAddrs) * m_nSizeofAddr + 4;
    std::vector<uint8_t> abyIblk;
    READ_OR_FAIL(m_nBaseAddr + nIblkAddr, nIblkLen, abyIblk);
    Cursor ci(abyIblk);
    if (!ci.Sig("EAIB")) return Fail(Status::ERR_SIGNATURE, "missing EAIB signature");
    if (Lookup3(abyIblk.data(), abyIblk.size() - 4, 0) != Cursor(abyIblk.data() + abyIblk.size() - 4, 4).U32())
        return Fail(Status::ERR_CORRUPT, "EAIB checksum mismatch");
    ci.Skip(2 + m_nSizeofAddr);
    const uint8_t *pabyIblkElmts = ci.Ptr();
    ci.Skip(nIdxBlkElmts * static_cast<size_t>(nElmtSize));
    std::vector<uint64_t> anIblkDblk(nIblkDblkAddrs), anIblkSblk(nIblkSblkAddrs);
    for (auto &a : anIblkDblk) a = ci.Addr(m_nSizeofAddr);
    for (auto &a : anIblkSblk) a = ci.Addr(m_nSizeofAddr);
    if (!ci.Ok()) return Fail(Status::ERR_CORRUPT, "truncated EAIB");

    // Swizzled chunk grid: the unlimited dimension becomes the slowest-varying one.
    const int nRank = oInfo.nRank;
    const std::vector<uint64_t> anMax = ChunkGrid(oInfo.anMaxDims, oInfo.anChunkDims);
    int nUnlim = -1;
    for (int i = 0; i < nRank; i++)
        if (anMax[i] == UNLIMITED) { if (nUnlim < 0) nUnlim = i; else return Fail(Status::ERR_INDEX, "extensible array with >1 unlimited dims"); }
    if (nUnlim < 0) return Fail(Status::ERR_INDEX, "extensible array without unlimited dimension");
    std::vector<int> anOrder;
    anOrder.push_back(nUnlim);
    for (int i = 0; i < nRank; i++) if (i != nUnlim) anOrder.push_back(i);
    std::vector<uint64_t> anScaled(nRank, 0);  // swizzled coordinates

    const uint64_t nChunkBytes = oInfo.ChunkBytes();
    auto Advance = [&]() {
        for (int d = nRank - 1; d >= 1; d--) {
            if (++anScaled[d] < anMax[anOrder[d]]) return;
            anScaled[d] = 0;
        }
        ++anScaled[0];
    };
    uint64_t nIdx = 0;
    auto EmitFrom = [&](Cursor &c) -> bool {
        ChunkRecord oRec;
        if (!DecodeArrayElement(c, m_nSizeofAddr, bFiltered, nSizeLen, nChunkBytes, oRec)) return false;
        if (oRec.nAddr != UNDEF_ADDR) {
            oRec.anOffset.resize(nRank);
            for (int i = 0; i < nRank; i++) oRec.anOffset[anOrder[i]] = anScaled[i] * oInfo.anChunkDims[anOrder[i]];
            aoChunks.push_back(std::move(oRec));
        }
        Advance();
        nIdx++;
        return true;
    };
    auto SkipN = [&](uint64_t n) { for (uint64_t k = 0; k < n && nIdx < nMaxIdxSet; k++) { Advance(); nIdx++; } };

    // Elements stored directly in the index block
    {
        Cursor c(pabyIblkElmts, nIdxBlkElmts * static_cast<size_t>(nElmtSize));
        for (unsigned i = 0; i < nIdxBlkElmts && nIdx < nMaxIdxSet; i++)
            if (!EmitFrom(c)) return Fail(Status::ERR_CORRUPT, "truncated EAIB elements");
    }

    // Reads one data block (optionally paged) and emits up to nDblkElmts elements.
    std::vector<uint8_t> abyBlk;
    auto ReadDblk = [&](uint64_t nAddr, uint64_t nDblkElmts, const uint8_t *pabyPageInit) -> Status {
        const uint64_t nThisBlk = std::min(nDblkElmts, nMaxIdxSet - nIdx);
        if (nAddr == UNDEF_ADDR) { SkipN(nThisBlk); return Status::OK; }
        const size_t nPrefix = 4 + 1 + 1 + m_nSizeofAddr + nArrOffSize;
        if (nDblkElmts <= nPageElmts) {
            if (!m_poReader->Read(m_nBaseAddr + nAddr, static_cast<size_t>(nPrefix + nDblkElmts * nElmtSize + 4), abyBlk))
                return Fail(Status::ERR_IO, "extensible array data block read failed");
            if (memcmp(abyBlk.data(), "EADB", 4) != 0) return Fail(Status::ERR_SIGNATURE, "missing EADB signature");
            if (Lookup3(abyBlk.data(), abyBlk.size() - 4, 0) != Cursor(abyBlk.data() + abyBlk.size() - 4, 4).U32())
                return Fail(Status::ERR_CORRUPT, "EADB checksum mismatch");
            Cursor c(abyBlk.data() + nPrefix, static_cast<size_t>(nDblkElmts * nElmtSize));
            for (uint64_t i = 0; i < nThisBlk; i++)
                if (!EmitFrom(c)) return Fail(Status::ERR_CORRUPT, "truncated EADB");
            return Status::OK;
        }
        if (!pabyPageInit) return Fail(Status::ERR_INDEX, "paged extensible array data block without bitmap");
        const uint64_t nPages = nDblkElmts / nPageElmts;
        const uint64_t nPageBytes = nPageElmts * nElmtSize + 4;
        const uint64_t nFirstPage = nAddr + nPrefix + 4;
        for (uint64_t p = 0; p < nPages && nIdx < nMaxIdxSet; p++) {
            const uint64_t nThis = std::min(nPageElmts, nMaxIdxSet - nIdx);
            if (!BitGet(pabyPageInit, p)) { SkipN(nThis); continue; }
            if (!m_poReader->Read(m_nBaseAddr + nFirstPage + p * nPageBytes, static_cast<size_t>(nPageBytes), abyBlk))
                return Fail(Status::ERR_IO, "extensible array page read failed");
            if (Lookup3(abyBlk.data(), abyBlk.size() - 4, 0) != Cursor(abyBlk.data() + abyBlk.size() - 4, 4).U32())
                return Fail(Status::ERR_CORRUPT, "extensible array page checksum mismatch");
            Cursor c(abyBlk);
            for (uint64_t i = 0; i < nThis; i++)
                if (!EmitFrom(c)) return Fail(Status::ERR_CORRUPT, "truncated extensible array page");
        }
        return Status::OK;
    };

    for (unsigned s = 0; s < nSblks && nIdx < nMaxIdxSet; s++) {
        const SblkInfo &si = aoSblk[s];
        if (s < nIblkSblks) {
            for (uint64_t d = 0; d < si.nDblks && nIdx < nMaxIdxSet; d++) {
                const uint64_t nDblk = si.nStartDblk + d;
                if (nDblk >= anIblkDblk.size()) return Fail(Status::ERR_CORRUPT, "bad extensible array data block index");
                Status eStatus = ReadDblk(anIblkDblk[static_cast<size_t>(nDblk)], si.nDblkElmts, nullptr);
                if (eStatus != Status::OK) return eStatus;
            }
            continue;
        }
        const uint64_t nSblkAddr = anIblkSblk[s - nIblkSblks];
        if (nSblkAddr == UNDEF_ADDR) { SkipN(si.nDblks * si.nDblkElmts); continue; }
        const bool bPaged = si.nDblkElmts > nPageElmts;
        const size_t nPageInitSize = bPaged ? static_cast<size_t>((si.nDblkElmts / nPageElmts + 7) / 8) : 0;
        const size_t nSblkLen = 4 + 1 + 1 + m_nSizeofAddr + nArrOffSize +
                                static_cast<size_t>(si.nDblks) * (nPageInitSize + m_nSizeofAddr) + 4;
        std::vector<uint8_t> abySblk;
        READ_OR_FAIL(m_nBaseAddr + nSblkAddr, nSblkLen, abySblk);
        if (memcmp(abySblk.data(), "EASB", 4) != 0) return Fail(Status::ERR_SIGNATURE, "missing EASB signature");
        if (Lookup3(abySblk.data(), abySblk.size() - 4, 0) != Cursor(abySblk.data() + abySblk.size() - 4, 4).U32())
            return Fail(Status::ERR_CORRUPT, "EASB checksum mismatch");
        const uint8_t *pabyInit = abySblk.data() + 4 + 1 + 1 + m_nSizeofAddr + nArrOffSize;
        Cursor ca(pabyInit + si.nDblks * nPageInitSize, static_cast<size_t>(si.nDblks) * m_nSizeofAddr);
        for (uint64_t d = 0; d < si.nDblks && nIdx < nMaxIdxSet; d++) {
            const uint64_t nDblkAddr = ca.Addr(m_nSizeofAddr);
            Status eStatus = ReadDblk(nDblkAddr, si.nDblkElmts, bPaged ? pabyInit + d * nPageInitSize : nullptr);
            if (eStatus != Status::OK) return eStatus;
        }
    }
    return Status::OK;
}

Status File::IterBTreeV2(const DatasetInfo &oInfo, std::vector<ChunkRecord> &aoChunks)
{
    BTree2Header oHdr;
    Status eStatus = ReadBTree2Header(oInfo.nIndexAddr, oHdr);
    if (eStatus != Status::OK) return eStatus;
    const bool bFiltered = oHdr.nType == 11;
    if (oHdr.nType != 10 && oHdr.nType != 11) return Fail(Status::ERR_INDEX, "unexpected chunk v2 B-tree type");
    const int nRank = oInfo.nRank;
    const size_t nFixed = m_nSizeofAddr + 8 * static_cast<size_t>(nRank) + (bFiltered ? 4 : 0);
    if (oHdr.nRecordSize < nFixed + (bFiltered ? 1 : 0)) return Fail(Status::ERR_CORRUPT, "bad chunk v2 B-tree record size");
    const unsigned nSizeLen = static_cast<unsigned>(oHdr.nRecordSize - nFixed);
    if (!bFiltered && nSizeLen != 0) return Fail(Status::ERR_CORRUPT, "bad chunk v2 B-tree record size");

    std::vector<std::vector<uint8_t>> aoRecords;
    eStatus = CollectBTree2Records(oHdr, oHdr.nRootAddr, oHdr.nDepth, oHdr.nRootRecords, aoRecords, 0);
    if (eStatus != Status::OK) return eStatus;
    const uint64_t nChunkBytes = oInfo.ChunkBytes();
    for (const auto &abyRec : aoRecords) {
        Cursor c(abyRec);
        ChunkRecord oRec;
        if (!DecodeArrayElement(c, m_nSizeofAddr, bFiltered, nSizeLen, nChunkBytes, oRec))
            return Fail(Status::ERR_CORRUPT, "truncated chunk v2 B-tree record");
        oRec.anOffset.resize(nRank);
        for (int i = 0; i < nRank; i++) oRec.anOffset[i] = c.U64() * oInfo.anChunkDims[i];
        if (!c.Ok()) return Fail(Status::ERR_CORRUPT, "truncated chunk v2 B-tree record");
        if (oRec.nAddr != UNDEF_ADDR) aoChunks.push_back(std::move(oRec));
    }
    return Status::OK;
}

}  // namespace NisarHDF5Native
