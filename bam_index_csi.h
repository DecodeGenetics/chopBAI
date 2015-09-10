#ifndef INCLUDE_SEQAN_BAM_IO_BAM_INDEX_CSI_H_
#define INCLUDE_SEQAN_BAM_IO_BAM_INDEX_CSI_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Bai
// ----------------------------------------------------------------------------

struct Csi_;
typedef Tag<Csi_> Csi;

// ----------------------------------------------------------------------------
// Helper Class BaiBamIndexBinData_
// ----------------------------------------------------------------------------

// Store the information of a bin.

struct CsiBamIndexBinData_
{
    __uint64 loffset;
    String<Pair<__uint64, __uint64> > chunkBegEnds;
};

// ----------------------------------------------------------------------------
// Spec CSI BamIndex
// ----------------------------------------------------------------------------

/*!
 * @class CsiBamIndex
 * @headerfile <seqan/bam_io.h>
 * @extends BamIndex
 * @brief Access to CSI (samtools-style).
 *
 * @signature template <>
 *            class BamIndex<Csi>;
 */

/*!
 * @fn CsiBamIndex::BamIndex
 * @brief Constructor.
 *
 * @signature BamIndex::BamIndex();
 *
 * @section Remarks
 *
 * Only the default constructor is provided.
 */

template <>
class BamIndex<Csi>
{
public:
    typedef std::map<__uint32, CsiBamIndexBinData_> TBinIndex_;

    __int32 _minShift;
    __int32 _depth;
    __uint64 _unalignedCount;

    String<__uint8> _aux;
    String<TBinIndex_> _binIndices;

    BamIndex() : _minShift(14), _depth(5), _unalignedCount(maxValue<__uint64>())
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function jumpToRegion()
// ----------------------------------------------------------------------------

static inline void
_csiReg2bins(String<__uint32> & list, __uint64 beg, __uint64 end, __int32 minShift, __int32 depth)
{
    int l = 0, t = 0;
    unsigned s = minShift + depth*3;
    for (--end; l <= depth; s -= 3, t += 1<<l*3, ++l)
    {
        unsigned b = t + (beg>>s);
        unsigned e = t + (end>>s);
        for (__uint32 i = b; i <= e; ++i)
            appendValue(list, i);
    }
}

template <typename TSpec>
inline bool
jumpToRegion(FormattedFile<Bam, Input, TSpec> & bamFile,
             bool & hasAlignments,
             __int32 refId,
             __int64 pos,
             __int64 posEnd,
             BamIndex<Csi> const & index)
{
    if (!isEqual(format(bamFile), Bam()))
        return false;

    hasAlignments = false;
    if (refId < 0)
        return false;  // Cannot seek to invalid reference.
    if (static_cast<unsigned>(refId) >= length(index._binIndices))
        return false;  // Cannot seek to invalid reference.

    // ------------------------------------------------------------------------
    // Compute offset in BGZF file.
    // ------------------------------------------------------------------------
    __uint64 offset = MaxValue<__uint64>::VALUE;

    // Retrieve the candidate bin identifiers for [pos, posEnd).
    String<__uint32> candidateBins;
    _csiReg2bins(candidateBins, pos, posEnd, index._minShift, index._depth);

    // Convert candidate bins into candidate offset.
    typedef std::set<__uint64> TOffsetCandidates;
    TOffsetCandidates offsetCandidates;
    typedef typename Iterator<String<__uint32>, Rooted>::Type TCandidateIter;
    for (TCandidateIter it = begin(candidateBins, Rooted()); !atEnd(it); goNext(it))
    {
        typedef typename std::map<__uint32, CsiBamIndexBinData_>::const_iterator TMapIter;
        TMapIter mIt = index._binIndices[refId].find(*it);
        if (mIt == index._binIndices[refId].end())
            continue;  // Candidate is not in index!

        typedef typename Iterator<String<Pair<__uint64, __uint64> > const, Rooted>::Type TBegEndIter;
        for (TBegEndIter it2 = begin(mIt->second.chunkBegEnds, Rooted()); !atEnd(it2); goNext(it2))
            offsetCandidates.insert(it2->i1);
    }

    // Search through candidate offsets, find rightmost possible.
    // Note that it is not necessarily the first.
    typedef typename TOffsetCandidates::const_iterator TOffsetCandidateIter;
    BamAlignmentRecord record;
    for (TOffsetCandidateIter candIt = offsetCandidates.begin(); candIt != offsetCandidates.end(); ++candIt)
    {
        setPosition(bamFile, *candIt);

        readRecord(record, bamFile);

        // std::cerr << "record.beginPos == " << record.beginPos << "\n";
        // __int32 endPos = record.beginPos + getAlignmentLengthInRef(record);
        if (record.rID != refId)
            continue;  // Wrong contig.
        if (!hasAlignments || record.beginPos <= pos)
        {
            // Found a valid alignment.
            hasAlignments = true;
            offset = *candIt;
        }

        if (record.beginPos >= posEnd)
            break;  // Cannot find overlapping any more.
    }

    if (offset != MaxValue<__uint64>::VALUE)
        setPosition(bamFile, offset);

    // Finding no overlapping alignment is not an error, hasAlignments is false.
    return true;
}

// ----------------------------------------------------------------------------
// Function jumpToOrphans()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TNameStore, typename TNameStoreCache>
bool jumpToOrphans(FormattedFile<Bam, Input, TSpec> & bamFile,
                   bool & hasAlignments,
                   BamIndex<Csi> const & index)
{

    SEQAN_FAIL("This does not work yet!");
    return false;
}

// ----------------------------------------------------------------------------
// Function getUnalignedCount()
// ----------------------------------------------------------------------------

// TODO(bkehr): Could be combined in a template function with implementation for BamIndex<Bai>
inline __uint64
getUnalignedCount(BamIndex<Csi> const & index)
{
    return index._unalignedCount;
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template<typename TStream>
bool
open(BamIndex<Csi> & index, TStream & fin)
{
    // Read the magic number.
    CharString buffer;
    resize(buffer, 4);

    fin.read(&buffer[0], 4);
 
    if (!fin.good())
        return false;

    if (buffer != "CSI\1")
        return false;  // Magic string is wrong.
        
    // Read the minShift.
    fin.read(reinterpret_cast<char *>(&index._minShift), 4);
    if (!fin.good())
        return false;
    
    // Read the depth.
    fin.read(reinterpret_cast<char *>(&index._depth), 4);
    if (!fin.good())
        return false;
    
    // Read auxiliary data.
    __int32 lAux = 0;
    fin.read(reinterpret_cast<char *>(&lAux), 4);
    if (!fin.good())
        return false;
    
    resize(index._aux, lAux);
    
    for (int i = 0; i < lAux; ++i)
    {
        fin.read(reinterpret_cast<char *>(&index._aux[i]), 1);
        if (!fin.good())
            return false;
    }

    // Read number of references.
    __int32 nRef = 0;
    fin.read(reinterpret_cast<char *>(&nRef), 4);
    if (!fin.good())
        return false;

    resize(index._binIndices, nRef);

    // Read bin index for each reference.
    for (int i = 0; i < nRef; ++i)
    {
        // Read number of bins.
        __int32 nBin = 0;
        fin.read(reinterpret_cast<char *>(&nBin), 4);
        if (!fin.good())
            return false;

        index._binIndices[i].clear();
        CsiBamIndexBinData_ data;

        // Read loffset and chunks for each bin.
        for (int j = 0; j < nBin; ++j)
        {
            clear(data.chunkBegEnds);

            // Read distinct bin.
            __uint32 bin = 0;
            fin.read(reinterpret_cast<char *>(&bin), 4);
            if (!fin.good())
                return false;

            // Read loffset.
            fin.read(reinterpret_cast<char *>(&data.loffset), 8);
            if (!fin.good())
                return false;

            // Read number of chunks.
            __int32 nChunk = 0;
            fin.read(reinterpret_cast<char *>(&nChunk), 4);
            if (!fin.good())
                return false;

            reserve(data.chunkBegEnds, nChunk);

            // Read begin and end for each chunk.
            for (int k = 0; k < nChunk; ++k)
            {
                __uint64 chunkBeg = 0;
                __uint64 chunkEnd = 0;
                fin.read(reinterpret_cast<char *>(&chunkBeg), 8);
                fin.read(reinterpret_cast<char *>(&chunkEnd), 8);
                if (!fin.good())
                    return false;
                appendValue(data.chunkBegEnds, Pair<__uint64>(chunkBeg, chunkEnd));
            }

            // Copy bin data into index.
            index._binIndices[i][bin] = data;
        }
    }

    if (!fin.good())
        return false;

    // Read (optional) number of alignments without coordinate.
    __uint64 nNoCoord = 0;
    fin.read(reinterpret_cast<char *>(&nNoCoord), 8);
    if (!fin.good())
    {
        fin.clear();
        nNoCoord = 0;
    }
    index._unalignedCount = nNoCoord;

    return true;
}

// ----------------------------------------------------------------------------

bool
open(BamIndex<Csi> & index, char const * filename)
{
    std::ifstream iss(filename);
    if (!iss.good())
        return false;  // Could not open file.
    
    if (guessFormatFromStream(iss, BgzfFile()))
    {
        bgzf_istream fin(iss);
        if (!fin.good())
            return false;  // Could not open file.
        return open(index, fin);
    }
    else
    {
        return open(index, iss);
    }
}

// ----------------------------------------------------------------------------
// Function saveIndex()
// ----------------------------------------------------------------------------

inline bool saveIndex(BamIndex<Csi> const & index, char const * filename)
{
    typedef BamIndex<Csi> const                TBamIndex;
    typedef TBamIndex::TBinIndex_ const        TBinIndex;
    typedef TBinIndex::const_iterator          TBinIndexIter;

    // Open output file.
    std::ofstream out(filename, std::ios::binary | std::ios::out);

    // Write header.
    out.write("CSI\1", 4);
    __int32 minShift = index._minShift;
    out.write(reinterpret_cast<char *>(&minShift), 4);
    __int32 depth = index._depth;
    out.write(reinterpret_cast<char *>(&depth), 4);
    
    // Write auxiliary data if present.
    __int32 lenAux = length(index._aux);
    out.write(reinterpret_cast<char *>(&lenAux), 4);
    for (int i = 0; i < lenAux; ++i)
        out.write(reinterpret_cast<char *>(index._aux[i]), 1);

    // Write out binning index.
    __int32 numRefSeqs = length(index._binIndices);
    out.write(reinterpret_cast<char *>(&numRefSeqs), 4);
    for (int i = 0; i < numRefSeqs; ++i)
    {
        TBinIndex const & binIndex = index._binIndices[i];

        // Write out binning index.
        __int32 numBins = binIndex.size();
        out.write(reinterpret_cast<char *>(&numBins), 4);
        for (TBinIndexIter itB = binIndex.begin(), itBEnd = binIndex.end(); itB != itBEnd; ++itB)
        {
            // Write out bin id.
            out.write(reinterpret_cast<char const *>(&itB->first), 4);
            // Write out loffset.
            __uint64 loffset = itB->second.loffset;
            out.write(reinterpret_cast<char *>(&loffset), 8);
            // Write out number of chunks.
            __uint32 numChunks = length(itB->second.chunkBegEnds);
            out.write(reinterpret_cast<char *>(&numChunks), 4);
            // Write out all chunks.
            typedef Iterator<String<Pair<__uint64> > const, Rooted>::Type TChunkIter;
            for (TChunkIter itC = begin(itB->second.chunkBegEnds); !atEnd(itC); goNext(itC))
            {
                out.write(reinterpret_cast<char const *>(&itC->i1), 8);
                out.write(reinterpret_cast<char const *>(&itC->i2), 8);
            }
        }
    }

    // Write the number of unaligned reads if set.
    if (index._unalignedCount != maxValue<__uint64>())
        out.write(reinterpret_cast<char const *>(&index._unalignedCount), 8);

    return out.good();  // false on error, true on success.
}

// ----------------------------------------------------------------------------
// Function buildIndex()
// ----------------------------------------------------------------------------

inline bool
buildIndex(BamIndex<Csi> & /*index*/, char const * /*filename*/)
{
    SEQAN_FAIL("This does not work yet!");
    return false;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_BAM_INDEX_CSI_H_
