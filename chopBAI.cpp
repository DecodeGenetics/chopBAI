#include <fstream>
#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>

#include "bam_index_csi.h"

using namespace seqan;


// -----------------------------------------------------------------------------

struct ChopBaiOptions {
    // Input arguments
    CharString bamfile;
    String<CharString> regions;

    // Output options
    CharString outputPrefix;
    bool writeLinear;
    
    ChopBaiOptions() :
        outputPrefix("."), writeLinear(false)
    {}
};


// -----------------------------------------------------------------------------

struct GenomicInterval {
    size_t chrId;
    __uint32 begin;
    __uint32 end;
};


// -----------------------------------------------------------------------------
// Function setupParser()
// -----------------------------------------------------------------------------

void setupParser(ArgumentParser & parser, ChopBaiOptions & options)
{
    setShortDescription(parser, "chops a bam index file into pieces");
    
    setVersion(parser, "0.1 beta");
    setDate(parser, DATE);
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM-FILE\\fP \\fIREGION1\\fP [... \\fIREGIONn\\fP]");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM-FILE\\fP \\fIREGION-FILE\\fP");
    
    addDescription(parser, "Writes small index files for the specified regions based on an existing bai or csi file for "
                           "the input bamfile. The regions have to be specified in the format \'chr:begin-end\' and can "
                           "be listed directly on the command line separated by spaces or in a file listing one region "
                           "per line. The program writes a smaller index file for each region to the directory "
                           "\'<output prefix>/<region>/<bamfile>.[bai|csi]\'. The output directories are created if "
                           "they do not exist.");
    
    // Required arguments.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM-FILE"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "REGIONS", true));

    // Options.
    addSection(parser, "Output options");
    addOption(parser, ArgParseOption("p", "prefix", "Output prefix.", ArgParseArgument::STRING, "STR"));
    addOption(parser, ArgParseOption("l", "linear", "Include linear index of BAI in the output."));

    // Set defualt values.
    setDefaultValue(parser, "prefix", "current directory");
    setDefaultValue(parser, "linear", options.writeLinear?"true":"false");
}


// -----------------------------------------------------------------------------
// Function getOptionValues()
// -----------------------------------------------------------------------------

void getOptionValues(ChopBaiOptions & options, ArgumentParser & parser)
{
    // Get argument values.
    getArgumentValue(options.bamfile, parser, 0);
    options.regions = getArgumentValues(parser, 1);
    
    // Get option values.
    if (isSet(parser, "prefix"))
        getOptionValue(options.outputPrefix, parser, "prefix");
    if (isSet(parser, "linear"))
        options.writeLinear = true;
}


// -----------------------------------------------------------------------------
// Function parseCommandLine()
// -----------------------------------------------------------------------------

int parseCommandLine(ChopBaiOptions & options, int argc, char const ** argv)
{
    // Setup the parser.
    ArgumentParser parser(argv[0]);
    setupParser(parser, options);

    // Parse the command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return 1;

    // Collect the option values.
    getOptionValues(options, parser);

    return 0;
}


// -----------------------------------------------------------------------------
// Function readRegions()
// -----------------------------------------------------------------------------

bool readRegions(String<CharString> & regions, CharString const & regionsFile)
{
    std::ifstream stream(toCString(regionsFile));
    if (!stream.is_open())
    {
        std::cerr << "ERROR: Could not open file listing the regions: " << regionsFile << std::endl;
        return 1;
    }
    
    std::string region;
    while (getline(stream, region))
        appendValue(regions, region);

    return 0;
}

// -----------------------------------------------------------------------------
// Function fileExists()
// -----------------------------------------------------------------------------

bool fileExists(CharString & filename)
{
    FILE *file = fopen(toCString(filename), "r");
    if (file != NULL)
    {
        fclose(file);
        return true;
    }
    return false;
}


// -----------------------------------------------------------------------------
// Function findIndexFile()
// -----------------------------------------------------------------------------

bool findIndexFile(CharString & indexfile, CharString & bamfile)
{
    // Check if BAI file exists: Given FILE.bam, first look for FILE.bam.bai and then for FILE.bai
    indexfile = bamfile;
    indexfile += ".bai";
    
    if (fileExists(indexfile))
        return 0;
    
    if (suffix(bamfile, length(bamfile)-4) == ".bam")
    {
        indexfile = prefix(bamfile, length(bamfile)-4);
        indexfile += ".bai";

        if (fileExists(indexfile))
            return 0;
    }

    // Check if CSI file exists: Given FILE.bam, first look for FILE.bam.csi and then for FILE.csi
    indexfile = bamfile;
    indexfile += ".csi";
    
    if (fileExists(indexfile))
        return 0;
    
    if (suffix(bamfile, length(bamfile)-4) == ".bam")
    {
        indexfile = prefix(bamfile, length(bamfile)-4);
        indexfile += ".csi";

        if (fileExists(indexfile))
            return 0;
    }
    
    std::cerr << "ERROR: Could not find .bai or .csi file for input bam file " << bamfile << std::endl;
    return 1;
}

// -----------------------------------------------------------------------------
// Function parseInterval()
// -----------------------------------------------------------------------------

bool parseInterval(GenomicInterval & interval, CharString & region, NameStoreCache<StringSet<CharString> > & refNames)
{
    typedef Iterator<CharString, Rooted>::Type TIter;
    
    // Read the chromosome name.
    CharString chrName;
    TIter it = begin(region);
    TIter itEnd = end(region);
    for (; it != itEnd; ++it)
    {
        if (*it == ':')
        {
            chrName = prefix(region, position(it));
            break;
        }
    }
    if (it == itEnd) return 1;
    interval.chrId = nameToId(refNames, chrName);
    
    // Read the begin and end position.
    ++it;
    size_t infixBeg = position(it);
    for (; it != itEnd; ++it)
    {
        if (*it == '-')
        {
            if (!lexicalCast(interval.begin, infix(region, infixBeg, position(it)))) return 1;
            if (!lexicalCast(interval.end, suffix(region, position(it) + 1))) return 1;
            if (interval.begin != 0) --interval.begin;
            return 0;
        }
    }
    
    return 1;
}


// -----------------------------------------------------------------------------
// Function parseIntervals()
// -----------------------------------------------------------------------------

bool parseIntervals(String<GenomicInterval> & intervals, String<CharString> & regions, CharString & bamfile)
{
    typedef NameStoreCache<StringSet<CharString> > TNamesCache;
    
    // Convert chromosome name into chrId using the header of bamfile.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamfile)))
    {
        std::cerr << "ERROR: Could not open " << bamfile << std::endl;
        return 1;
    }
    BamHeader header;
    readHeader(header, bamFileIn);
    TNamesCache refNames = contigNamesCache(context(bamFileIn));
    
    if (length(regions) == 1)
    {
        // Check if it is a single region.
        GenomicInterval interval;
        if (parseInterval(interval, regions[0], refNames) == 0)
        {
            appendValue(intervals, interval);
            return 0;
        }
        
        // Read file listing the regions.
        CharString regionsFile = regions[0];
        clear(regions);
        readRegions(regions, regionsFile);
    }
    
    for(size_t i = 0; i < length(regions); ++i)
    {
        GenomicInterval interval;
        if (parseInterval(interval, regions[i], refNames) != 0)
        {
            std::cerr << "ERROR: Could not parse region " << i << ": " << regions[i] << std::endl;
            std::cerr << "       Please specify region in format chr:beg-end." << std::endl;
            return 1;
        }
        appendValue(intervals, interval);
    }
    
    return 0;
}

// -----------------------------------------------------------------------------
// Function cropInterval()
// -----------------------------------------------------------------------------

void cropInterval(BamIndex<Csi> & outcsi, BamIndex<Csi> & incsi, GenomicInterval & interval, bool )
{
    // Initialize the output bam index
    outcsi._minShift = incsi._minShift;
    outcsi._depth = incsi._depth;
    __int32 nRef = length(incsi._binIndices);
    resize(outcsi._binIndices, nRef);

    // --- Crop the region from the bin index ---
    
    String<__uint32> candidateBins;
    _csiReg2bins(candidateBins, interval.begin, interval.end, incsi._minShift, incsi._depth);
    
    typedef Iterator<String<__uint32>, Rooted>::Type TCandidateIter;
    for (TCandidateIter it = begin(candidateBins, Rooted()); !atEnd(it); ++it)
    {
        typedef std::map<__uint32, CsiBamIndexBinData_>::const_iterator TMapIter;
        TMapIter mIt = incsi._binIndices[interval.chrId].find(*it);
        if (mIt == incsi._binIndices[interval.chrId].end())
            continue;  // Candidate is not in index!

        // TODO: What about chunks.loffset?
        CsiBamIndexBinData_ chunks;
        typedef Iterator<String<Pair<__uint64, __uint64> > const, Rooted>::Type TBegEndIter;
        for (TBegEndIter it2 = begin(mIt->second.chunkBegEnds, Rooted()); !atEnd(it2); goNext(it2))
            appendValue(chunks.chunkBegEnds, *it2);
        
        if (length(chunks.chunkBegEnds) > 0)
            outcsi._binIndices[interval.chrId][*it] = chunks;
    }
}

// -----------------------------------------------------------------------------

void cropInterval(BamIndex<Bai> & outbai, BamIndex<Bai> & inbai, GenomicInterval & interval, bool writeLinear)
{
    // Initialize the output bam index
    __int32 nRef = length(inbai._linearIndices);
    resize(outbai._linearIndices, nRef);
    resize(outbai._binIndices, nRef);

    // --- Crop the region from the linear index ---

    unsigned windowIdx = interval.begin >> 14;  // Linear index consists of 16kb windows.
    unsigned windowEndIdx = (interval.end >> 14) + 1;

    __uint64 linearMinOffset = 0;
    
    if (windowIdx < length(inbai._linearIndices[interval.chrId]))
    {
        linearMinOffset = inbai._linearIndices[interval.chrId][windowIdx];
        
        __uint64 minOffset = 0u;
        for (unsigned i = interval.chrId; i < length(inbai._linearIndices); ++i)
        {
            if (!empty(inbai._linearIndices[i]))
            {
                minOffset = front(inbai._linearIndices[i]);
                break;
            }
        }
    
        // Crop the linear index if user asked for it.
        if (writeLinear)
        {
            for (unsigned i = 0; i < windowIdx; ++i)
                appendValue(outbai._linearIndices[interval.chrId], minOffset);
            for (unsigned i = windowIdx; i < _min(windowEndIdx, length(inbai._linearIndices[interval.chrId])); ++i)
                appendValue(outbai._linearIndices[interval.chrId], inbai._linearIndices[interval.chrId][i]);
        }
    }
    else  // set linearMinOffset to next non-zero entry in linear indices
    {
        if (empty(inbai._linearIndices[interval.chrId]))
        {
            for (unsigned i = interval.chrId; i < length(inbai._linearIndices); ++i)
            {
                if (!empty(inbai._linearIndices[i]))
                {
                    linearMinOffset = front(inbai._linearIndices[i]);
                    if (linearMinOffset != 0u)
                        break;
                    for (unsigned j = 1; j < length(inbai._linearIndices[i]); ++j)
                    {
                        if (inbai._linearIndices[i][j] != 0u)
                        {
                            linearMinOffset = inbai._linearIndices[i][j];
                            break;
                        }
                    }
                    if (linearMinOffset != 0u)
                        break;
                }
            }
        }
        else
        {
            linearMinOffset = back(inbai._linearIndices[interval.chrId]);
        }
    }
    
    // --- Crop the region from the bin index ---
    
    String<__uint16> candidateBins;
    _baiReg2bins(candidateBins, interval.begin, interval.end);
    
    typedef Iterator<String<__uint16>, Rooted>::Type TCandidateIter;
    for (TCandidateIter it = begin(candidateBins, Rooted()); !atEnd(it); ++it)
    {
        typedef std::map<__uint32, BaiBamIndexBinData_>::const_iterator TMapIter;
        TMapIter mIt = inbai._binIndices[interval.chrId].find(*it);
        if (mIt == inbai._binIndices[interval.chrId].end())
            continue;  // Candidate is not in index!

        BaiBamIndexBinData_ chunks;
        typedef Iterator<String<Pair<__uint64, __uint64> > const, Rooted>::Type TBegEndIter;
        for (TBegEndIter it2 = begin(mIt->second.chunkBegEnds, Rooted()); !atEnd(it2); goNext(it2))
            if (it2->i2 >= linearMinOffset)
                appendValue(chunks.chunkBegEnds, *it2);
        
        if (length(chunks.chunkBegEnds) > 0)
            outbai._binIndices[interval.chrId][*it] = chunks;
    }
}


// -----------------------------------------------------------------------------
// Function printBamIndex()                              // only for debugging
// -----------------------------------------------------------------------------

void printBamIndex(BamIndex<Csi> & csi)
{
    std::cout << "Number of bits for the minimal interval: " << csi._minShift << std::endl;
    std::cout << "Depth of the binning index: " << csi._depth << std::endl;
    std::cout << "Length of auxiliary data: " << length(csi._aux) << std::endl;
    std::cout << "Length of bin indices: " << length(csi._binIndices) << std::endl;
    std::cout << std::endl;
    
    std::cout << "Auxiliary data: " << std::endl;    
    for (unsigned i = 0; i < length(csi._aux); ++i)
        std::cout << "  Aux data " << i << ": " << csi._aux[i] << std::endl;
    std::cout << std::endl;

    for (unsigned i = 0; i < length(csi._binIndices); ++i)
    {
        if (csi._binIndices[i].size() == 0) continue;
        std::cout << "Bin index " << i << " has size " << csi._binIndices[i].size() << std::endl;
        for (BamIndex<Csi>::TBinIndex_::iterator it = csi._binIndices[i].begin(); it != csi._binIndices[i].end(); ++it)
        {
            std::cout << "  Bin " << it->first << ", loffset " << (it->second).loffset << " has " << length((it->second).chunkBegEnds) << " chunks." << std::endl;
            for (unsigned j = 0; j < length((it->second).chunkBegEnds); ++j)
                std::cout << "    chunkBeg: " << (it->second).chunkBegEnds[j].i1
                          << ", chunkEnd: " << (it->second).chunkBegEnds[j].i2
                          << ", diff: " << ((it->second).chunkBegEnds[j].i2 - (it->second).chunkBegEnds[j].i1) << std::endl;
        }
        std::cout << std::endl;
    }
}

// -----------------------------------------------------------------------------

void printBamIndex(BamIndex<Bai> & bai)
{
    std::cout << "Length of bin indices: " << length(bai._binIndices) << std::endl;
    std::cout << "Length of linear indices: " << length(bai._linearIndices) << std::endl;
    std::cout << std::endl;

    for (unsigned i = 0; i < length(bai._binIndices); ++i)
    {
        if (bai._binIndices[i].size() == 0) continue;
        std::cout << "Bin index " << i << " has size " << bai._binIndices[i].size() << std::endl;
        for (BamIndex<Bai>::TBinIndex_::iterator it = bai._binIndices[i].begin(); it != bai._binIndices[i].end(); ++it)
        {
            std::cout << "  Bin " << it->first << " has " << length((it->second).chunkBegEnds) << " chunks." << std::endl;
            for (unsigned j = 0; j < length((it->second).chunkBegEnds); ++j)
                std::cout << "    chunkBeg: " << (it->second).chunkBegEnds[j].i1
                          << ", chunkEnd: " << (it->second).chunkBegEnds[j].i2
                          << ", diff: " << ((it->second).chunkBegEnds[j].i2 - (it->second).chunkBegEnds[j].i1) << std::endl;
        }
        std::cout << std::endl;
    }

    for (unsigned i = 0; i < length(bai._linearIndices); ++i)
    {
        if (length(bai._linearIndices[i]) == 0) continue;
        std::cout << "Linear index " << i << " has size " << length(bai._linearIndices[i]) << std::endl;
        for (unsigned j = _max(0, (int)length(bai._linearIndices[i])-100); j < length(bai._linearIndices[i]); ++j)
            std::cout << "  " << j << "/" << length(bai._linearIndices[i]) << ": " << bai._linearIndices[i][j] << std::endl;
        std::cout << std::endl;
    }
}


// -----------------------------------------------------------------------------
// Function saveIndex()
// -----------------------------------------------------------------------------

bool saveIndex(BamIndex<Bai> const & index, char const * filename)
{
    // Open output file.
    std::ofstream out(filename, std::ios::binary | std::ios::out);
    if (!out.is_open())
    {
        std::cerr << "ERROR: Could not open output file: " << filename << std::endl;
        return 1;
    }

    SEQAN_ASSERT_EQ(length(index._binIndices), length(index._linearIndices));

    // Write header.
    out.write("BAI\1", 4);
    __int32 numRefSeqs = length(index._binIndices);
    out.write(reinterpret_cast<char *>(&numRefSeqs), 4);

    // Write out indices.
    typedef BamIndex<Bai> const                TBamIndex;
    typedef TBamIndex::TBinIndex_ const        TBinIndex;
    typedef TBinIndex::const_iterator          TBinIndexIter;
    typedef TBamIndex::TLinearIndex_           TLinearIndex;
    for (int i = 0; i < numRefSeqs; ++i)
    {
        TBinIndex const & binIndex = index._binIndices[i];
        TLinearIndex const & linearIndex = index._linearIndices[i];

        // Write out binning index.
        __int32 numBins = binIndex.size();
        out.write(reinterpret_cast<char *>(&numBins), 4);
        for (TBinIndexIter itB = binIndex.begin(), itBEnd = binIndex.end(); itB != itBEnd; ++itB)
        {
            // Write out bin id.
            out.write(reinterpret_cast<char const *>(&itB->first), 4);

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

        // Write out linear index.
        __int32 numIntervals = length(linearIndex);
        out.write(reinterpret_cast<char *>(&numIntervals), 4);
        typedef Iterator<String<__uint64> const, Rooted>::Type TLinearIndexIter;
        for (TLinearIndexIter it = begin(linearIndex, Rooted()); !atEnd(it); goNext(it))
            out.write(reinterpret_cast<char const *>(&*it), 8);
    }

    // Write the number of unaligned reads if set.
    if (index._unalignedCount != maxValue<__uint64>())
        out.write(reinterpret_cast<char const *>(&index._unalignedCount), 8);

    return out.good();  // false on error, true on success.
}


// -----------------------------------------------------------------------------
// Function chopIndex()
// -----------------------------------------------------------------------------

template<typename TTag>
int chopIndex(String<GenomicInterval> & intervals, CharString & indexfile, ChopBaiOptions & options, TTag)
{
    // Load the input bam index.
    BamIndex<TTag> inIndex;
    if (open(inIndex, toCString(indexfile)) != true)
    {
        std::cerr << "ERROR: Open failed on bam index file " << indexfile << std::endl;
        return 1;
    }
    
    // Crop the filename from the path of the index file.
    size_t i = length(indexfile);
    while (i > 0u)
    {
        --i;
        if (indexfile[i] == '/') break;
    }
    CharString indexfilename = suffix(indexfile, i);

    // Iterate regions.
    for (unsigned i = 0; i < length(intervals); ++i)
    {
        // Crop the region from the input bam index.
        BamIndex<TTag> outIndex;
        cropInterval(outIndex, inIndex, intervals[i], options.writeLinear);
        
        // Create output directory if not exists.
        std::stringstream outdir;
        outdir << options.outputPrefix << "/" << options.regions[i];
        mkdir(toCString(outdir.str()), 0755);
        
        std::stringstream outfile;
        outfile << outdir.str() << "/" << indexfilename;
    
        // Write the output bam index for the region.
        if (!saveIndex(outIndex, toCString(outfile.str())))
            return 1;
    }
    
    return 0;
}



// -----------------------------------------------------------------------------
// Function main()
// -----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse command line parameters.
    ChopBaiOptions options;
    if (parseCommandLine(options, argc, argv) != 0)
        return 1;

    // Parse the regions.
    String<GenomicInterval> intervals;
    if (parseIntervals(intervals, options.regions, options.bamfile) != 0)
        return 1;

    // Look for the index file given a BAM file.
    CharString indexfile;
    if (findIndexFile(indexfile, options.bamfile) != 0)
        return 1;
    
    // Chop the index file.
    if (suffix(indexfile, length(indexfile) - 3) == "bai") // Found BAI file
    {
        if (chopIndex(intervals, indexfile, options, Bai()) != 0)
            return 1;
    }
    else if (suffix(indexfile, length(indexfile) - 3) == "csi") // Found CSI file
    {
        if (chopIndex(intervals, indexfile, options, Csi()) != 0)
            return 1;
    }
    
    return 0;
}
