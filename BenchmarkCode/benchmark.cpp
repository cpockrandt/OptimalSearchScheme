#include <benchmark/benchmark.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

#include <sdsl/bit_vectors.hpp>

#include <type_traits>

#include "bav.h"
#include "common.h"
#include "paper_optimum_schemes.h"
#include "find2_index_approx_extension.h"

using namespace sdsl;
using namespace seqan;

typedef Index<StringSet<String<Dna, Alloc<>>, Owner<ConcatDirect<> > >, TIndexConfig> TIndex;
typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

TIndex fm_index;
StringSet<DnaString> reads;

uint64_t no_verifications;

/*class OSSContext
{
public:
    constexpr bool itv = true;
    template <typename TText, typename TIndex, typename TIndexSpec,
              size_t nbrBlocks>
    bool itvConditionComp(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                      uint32_t const needleLeftPos,
                      uint32_t const needleRightPos,
                      uint8_t const errors,
                      OptimalSearch<nbrBlocks> const & s,
                      uint8_t const blockIndex)
    {
        return(itv && blockIndex > 0/* && iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1 + 0) * 2* /);
    }
};*/

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


void BM_HammingDistance(benchmark::State& state, uint8_t const maxErrors)
{
    typedef HammingDistance TDistanceTag;

    TIter it(fm_index);

    uint64_t hitsNbr, uniqueHits;
    auto delegate = [&hitsNbr](auto const &it, DnaString const & /*read*/, unsigned const errors = 0) {
        ++hitsNbr;
        unsigned x = 0;
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
            x += getOccurrences(it)[i].i2;
    };

    for (auto _ : state)
    {
        hitsNbr = 0;
        uniqueHits = 0;
        for (unsigned i = 0; i < length(reads); ++i)
        {
            uint64_t oldHits = hitsNbr;
            _findBacktracking(it.revIter, reads[i], begin(reads[i], Standard()), 0, (signed) maxErrors, delegate,
                              TDistanceTag());
            reverseComplement(reads[i]);
            _findBacktracking(it.revIter, reads[i], begin(reads[i], Standard()), 0, (signed) maxErrors, delegate,
                              TDistanceTag());
            benchmark::DoNotOptimize(uniqueHits += oldHits != hitsNbr);
        }
        // std::cout << "Backtracking: " << ((double)((time*100)/CLOCKS_PER_SEC)/100) << " s. "
        //           << "Hits: " << uniqueHits << " (" << hitsNbr << ")" << std::endl;
        std::cout       << "Hits: " << uniqueHits << " (" << hitsNbr << ")" << std::endl;
    }
}

template <class TSearchScheme>
void BM_HammingDistance(benchmark::State& state, TSearchScheme scheme)
{
    typedef HammingDistance TDistanceTag;

    TIter it(fm_index);
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(reads[0]));

    uint64_t hitsNbr, uniqueHits;
    auto delegate = [&hitsNbr](auto const &it, DnaString const & /*read*/, unsigned const errors = 0) {
        ++hitsNbr;
        unsigned x = 0;
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
            x += getOccurrences(it)[i].i2;
    };


    for (auto _ : state)
    {
        hitsNbr = 0;
        uniqueHits = 0;
        for (unsigned i = 0; i < length(reads); ++i)
        {
            uint64_t oldHits = hitsNbr;
            _optimalSearchScheme(delegate, it, reads[i], scheme, TDistanceTag());
            reverseComplement(reads[i]);
            _optimalSearchScheme(delegate, it, reads[i], scheme, TDistanceTag());
            benchmark::DoNotOptimize(uniqueHits += oldHits != hitsNbr);
        }
        std::cout       << "Hits: " << uniqueHits << " (" << hitsNbr << ")" << std::endl;
//         std::cout << "Hits SS: " << uniqueHits << " (" << hitsNbr << ")" << std::endl; // Hits compare
        // std::cout << "Opt.-Schemes: " << ((double)((time*100)/CLOCKS_PER_SEC)/100) << " s. "
        //           << "Hits: " << uniqueHits << " (" << hitsNbr << ")" << std::endl;
    }
}

void BM_EditDistance(benchmark::State& state, uint8_t const maxErrors)
{
    typedef EditDistance TDistanceTag;

    TIter it(fm_index);

    uint64_t hitsNbr, uniqueHits;
    auto delegate = [&hitsNbr](auto const &it, DnaString const & /*read*/, unsigned const errors = 0) {
        ++hitsNbr;
        unsigned x = 0;
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
            x += getOccurrences(it)[i].i2;
    };

    for (auto _ : state)
    {
        hitsNbr = 0;
        uniqueHits = 0;
        for (unsigned i = 0; i < length(reads); ++i)
        {
            uint64_t oldHits = hitsNbr;
            _findBacktracking(it.revIter, reads[i], begin(reads[i], Standard()), 0, (signed) maxErrors, delegate,
                              TDistanceTag());
            reverseComplement(reads[i]);
            _findBacktracking(it.revIter, reads[i], begin(reads[i], Standard()), 0, (signed) maxErrors, delegate,
                              TDistanceTag());
            benchmark::DoNotOptimize(uniqueHits += oldHits != hitsNbr);
        }
        // std::cout << "Backtracking: " << ((double)((time*100)/CLOCKS_PER_SEC)/100) << " s. "
        //           << "Hits: " << uniqueHits << " (" << hitsNbr << ")" << std::endl;
    }
}

template <class TSearchScheme>
void BM_EditDistance(benchmark::State& state, TSearchScheme scheme)
{
    typedef EditDistance TDistanceTag;

    TIter it(fm_index);
    _optimalSearchSchemeComputeFixedBlocklength(scheme, length(reads[0]));

    uint64_t hitsNbr, uniqueHits;
    auto delegate = [&hitsNbr](auto const &it, DnaString const & /*read*/, unsigned const errors = 0) {
        ++hitsNbr;
        unsigned x = 0;
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
            x += getOccurrences(it)[i].i2;
    };

    for (auto _ : state)
    {
        hitsNbr = 0;
        uniqueHits = 0;
        for (unsigned i = 0; i < length(reads); ++i)
        {
            uint64_t oldHits = hitsNbr;
            _optimalSearchScheme(delegate, it, reads[i], scheme, TDistanceTag());
            reverseComplement(reads[i]);
            _optimalSearchScheme(delegate, it, reads[i], scheme, TDistanceTag());
            benchmark::DoNotOptimize(uniqueHits += oldHits != hitsNbr);
        }
        // std::cout << "Opt.-Schemes: " << ((double)((time*100)/CLOCKS_PER_SEC)/100) << " s. "
        //           << "Hits: " << uniqueHits << " (" << hitsNbr << ")" << std::endl;
    }
}

template <typename TPredictify>
void BM_010Seeds(benchmark::State& state, uint8_t const maxErrors, bool const indels_param, TPredictify & predictify)
{
    //bool const indels_param = false;

    TIter it(fm_index);

    uint64_t hitsNbr, uniqueHits;
    auto delegate = [&hitsNbr](auto const textpos) {
        ++hitsNbr;
	unsigned x = 0;
        benchmark::DoNotOptimize(x = getSeqOffset(textpos));
        //unsigned x = textpos.i2;
        //for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
        //    x += getOccurrences(it)[i].i2;
    };

    for (auto _ : state)
    {
        hitsNbr = 0;
        uniqueHits = 0;
        for (unsigned i = 0; i < length(reads); ++i)
        {
            uint64_t oldHits = hitsNbr;
            search(delegate, it, predictify, maxErrors, reads[i], indels_param);
            reverseComplement(reads[i]);
            search(delegate, it, predictify, maxErrors, reads[i], indels_param);
            benchmark::DoNotOptimize(uniqueHits += oldHits != hitsNbr);
        }
        std::cout << "Hits 01*0: " << uniqueHits << " (" << hitsNbr << ")" << std::endl; // Hits compare
    }
}

auto predictify_bidirectional = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                                    unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
{
    return false;
};

auto predictify_unidirectional = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                                     unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
{
    return pos_right == length(pattern) - 1;
};

class OSSContextOn
{
public:
    template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks>
    bool itvConditionComp(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter, uint32_t const needleLeftPos, uint32_t const needleRightPos, uint8_t const errors, OptimalSearch<nbrBlocks> const & s, uint8_t const blockIndex)
    { return true; //countOccurrences(iter) < 10;
    }
};
OSSContextOn ossContextOn;

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
//iOSSContextOcc<15> ossContextOcc15;
OSSContextOcc<50> ossContextOcc50;

template <unsigned blocks>
class OSSContextIndex
{
public:
    template <typename TText, typename TIndex, typename TIndexSpec, size_t nbrBlocks>
    bool itvConditionComp(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter, uint32_t const needleLeftPos, uint32_t const needleRightPos, uint8_t const errors, OptimalSearch<nbrBlocks> const & s, uint8_t const blockIndex)
    {
        //return needleRightPos - needleLeftPos > 95;
        //return 4*(blockIndex) < countOccurrences(iter);
        
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
OSSContextIndex<2> ossContextIndex2;
OSSContextIndex<3> ossContextIndex3;

BENCHMARK_CAPTURE(OSS_AnyDistance, 1_OSS_itv_off  , ossContextOff  , OptimalSearchSchemes<0, 1>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 1_OSS_itv_occ25, ossContextOcc25, OptimalSearchSchemes<0, 1>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 1_OSS_itv_occ50, ossContextOcc50, OptimalSearchSchemes<0, 1>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 1_OSS_itv_ind1, ossContextIndex1, OptimalSearchSchemes<0, 1>::VALUE)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(OSS_AnyDistance, 2_OSS_itv_off  , ossContextOff  , OptimalSearchSchemes<0, 2>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 2_OSS_itv_occ25, ossContextOcc25, OptimalSearchSchemes<0, 2>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 2_OSS_itv_occ50, ossContextOcc50, OptimalSearchSchemes<0, 2>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 2_OSS_itv_ind1, ossContextIndex1, OptimalSearchSchemes<0, 2>::VALUE)->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(OSS_AnyDistance, 3_TOP_itv_off  , ossContextOff  , PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 3_TOP_itv_occ25, ossContextOcc25, PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 3_TOP_itv_occ50, ossContextOcc50, PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(OSS_AnyDistance, 3_TOP_itv_ind1, ossContextIndex1, PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);

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
    CharString indexPath, readsPath, bitvectorPath;
    getOptionValue(indexPath, parser, "genome");
    getOptionValue(readsPath, parser, "reads");

    open(fm_index, toCString(indexPath), OPEN_RDONLY);
    std::cerr << "Index loaded." << std::endl;
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
    std::cerr << "Reads loaded." << std::endl;

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
