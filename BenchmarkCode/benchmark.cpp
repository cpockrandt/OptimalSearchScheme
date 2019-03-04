#include <benchmark/benchmark.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

#include <type_traits>

#include "bav.h"
#include "common.h"
#include "paper_optimum_schemes.h"
#include "find2_index_approx_extension.h"

using namespace seqan;

typedef Index<StringSet<String<Dna, Alloc<>>, Owner<ConcatDirect<> > >, TIndexConfig> TIndex;
typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

TIndex fm_index;
StringSet<DnaString> reads;

uint64_t no_verifications;

class OSSContextOff
{
public:
    template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks>
    bool itvConditionComp(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter, uint32_t const needleLeftPos, uint32_t const needleRightPos, uint8_t const errors, OptimalSearch<nbrBlocks> const & s, uint8_t const blockIndex)
    { return false; }
};
OSSContextOff ossContextOff;

template <unsigned occ>
class OSSContextOcc
{
public:
    template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks>
    bool itvConditionComp(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter, uint32_t const needleLeftPos, uint32_t const needleRightPos,
			  uint8_t const errors, OptimalSearch<nbrBlocks> const & s, uint8_t const blockIndex)
    {
        if (countOccurrences(iter) < occ)
        {
            no_verifications += countOccurrences(iter);
            return true;
        }
        return false;
        //return countOccurrences(iter) < occ;
    }
};
OSSContextOcc<25> ossContextOcc25;
OSSContextOcc<50> ossContextOcc50;

template <unsigned blocks>
class OSSContextIndex
{
public:
    template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks>
    bool itvConditionComp(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter, uint32_t const needleLeftPos, uint32_t const needleRightPos, uint8_t const errors, OptimalSearch<nbrBlocks> const & s, uint8_t const blockIndex)
    {
        if (blockIndex >= s.pi.size() - blocks)
        {
            no_verifications += countOccurrences(iter);
            return true;
        }    
        return false;
        //return blockIndex >= s.pi.size() - 1;
    }
};
OSSContextIndex<1> ossContextIndex1;

template <typename TOSSContext, typename TSearchScheme>
void OSS_AnyDistance(benchmark::State& state, TOSSContext & ossContext, TSearchScheme scheme)
{
    typedef HammingDistance TDistanceTag;

    TIter it(fm_index);

    uint64_t hitsNbr, uniqueHits;

    auto delegate = [&hitsNbr](TOSSContext & /*ossContext*/, auto const & it, DnaString const & /*needle*/, uint32_t const /*needleId*/, uint8_t /*errors*/, bool const /*rev*/)
    {
        ++hitsNbr;
        unsigned x = 0;
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
            x += getOccurrences(it)[i].i2;
    };
    auto delegateDirect = [&hitsNbr](TOSSContext & /*ossContext*/, Pair<uint16_t, uint32_t> const & pos, Pair<uint16_t, uint32_t> const & posEnd, DnaString const & needle, uint32_t const needleId, uint8_t const errors)
    {
        ++hitsNbr;
        //unsigned x = pos.i2;
    };

    std::vector<std::pair<TBitvector, TSupport>> empty_bitvectors;
    calcConstParameters(scheme);
    for (auto _ : state)
    {
        hitsNbr = 0;
        uniqueHits = 0;
        no_verifications = 0;
        for (unsigned i = 0; i < length(reads); ++i)
        {
            uint64_t oldHits = hitsNbr;
            find(ossContext, delegate, delegateDirect, it, empty_bitvectors, scheme, reads[i], i, TDistanceTag());
            reverseComplement(reads[i]);
            find(ossContext, delegate, delegateDirect, it, empty_bitvectors, scheme, reads[i], i, TDistanceTag());
            benchmark::DoNotOptimize(uniqueHits += oldHits != hitsNbr);
        }
        // std::cout << "Backtracking: " << ((double)((time*100)/CLOCKS_PER_SEC)/100) << " s. ";
        // std::cout       << "Hits: " << uniqueHits << " (" << hitsNbr << ")" << std::endl;
    }
    std::cout       << "Hits: " << uniqueHits << " (" << hitsNbr << ")" << ": verifications: " << no_verifications << std::endl;
}

BENCHMARK_CAPTURE(OSS_AnyDistance, 1_OSS_itv_off  , ossContextOff   , PaperOptimumSearchSchemes<1>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 1_OSS_itv_occ25, ossContextOcc25 , PaperOptimumSearchSchemes<1>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 1_OSS_itv_occ50, ossContextOcc50 , PaperOptimumSearchSchemes<1>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 1_OSS_itv_ind1 , ossContextIndex1, PaperOptimumSearchSchemes<1>::VALUE_plus_one)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(OSS_AnyDistance, 2_OSS_itv_off  , ossContextOff   , PaperOptimumSearchSchemes<2>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 2_OSS_itv_occ25, ossContextOcc25 , PaperOptimumSearchSchemes<2>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 2_OSS_itv_occ50, ossContextOcc50 , PaperOptimumSearchSchemes<2>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 2_OSS_itv_ind1 , ossContextIndex1, PaperOptimumSearchSchemes<2>::VALUE_plus_two)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(OSS_AnyDistance, 3_OSS_itv_off  , ossContextOff   , PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 3_OSS_itv_occ25, ossContextOcc25 , PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 3_OSS_itv_occ50, ossContextOcc50 , PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 3_OSS_itv_ind1 , ossContextIndex1, PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(OSS_AnyDistance, 1_BACK_itv_off  , ossContextOff   , BacktrackingSchemes<1>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 1_BACK_itv_occ25, ossContextOcc25 , BacktrackingSchemes<1>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 1_BACK_itv_occ50, ossContextOcc50 , BacktrackingSchemes<1>::VALUE)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(OSS_AnyDistance, 2_BACK_itv_off  , ossContextOff   , BacktrackingSchemes<2>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 2_BACK_itv_occ25, ossContextOcc25 , BacktrackingSchemes<2>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 2_BACK_itv_occ50, ossContextOcc50 , BacktrackingSchemes<2>::VALUE)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(OSS_AnyDistance, 3_BACK_itv_off  , ossContextOff   , BacktrackingSchemes<3>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 3_BACK_itv_occ25, ossContextOcc25 , BacktrackingSchemes<3>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 3_BACK_itv_occ50, ossContextOcc50 , BacktrackingSchemes<3>::VALUE)->Unit(benchmark::kMillisecond);

// BENCHMARK_MAIN();
int main(int argc, char** argv)
{
    // Argument parser
    ArgumentParser parser("SearchSchemes - Benchmarking");
    addDescription(parser,
        "App for creating the benchmark of Optimum Search Schemes from the paper.");

    addOption(parser, ArgParseOption("G", "genome", "Path to the indexed genome", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "genome");

    addOption(parser, ArgParseOption("R", "reads", "Path to the reads", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "reads", "fa fasta fastq");
	setRequired(parser, "reads");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    CharString indexPath, readsPath;
    getOptionValue(indexPath, parser, "genome");
    getOptionValue(readsPath, parser, "reads");

    open(fm_index, toCString(indexPath), OPEN_RDONLY);
    StringSet<CharString> ids;
    SeqFileIn seqFileIn(toCString(readsPath));
    readRecords(ids, reads, seqFileIn);

    for (unsigned i = 1; i < length(reads); ++i)
    {
        if (length(reads[i]) != length(reads[0]))
        {
            std::cerr << "ERROR: Not all reads have the same length." << std::endl;
            return 1;
        }
    }

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
