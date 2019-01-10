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


class OSSContext
{
public:
    bool itv = true;
    template <typename TText, typename TIndex, typename TIndexSpec,
              size_t nbrBlocks>
    bool itvConditionComp(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                      uint32_t const needleLeftPos,
                      uint32_t const needleRightPos,
                      uint8_t const errors,
                      OptimalSearch<nbrBlocks> const & s,
                      uint8_t const blockIndex)
    {
        return(itv && iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < (s.pi.size() - blockIndex - 1 + 0) * 2);
    }
};

template <class TSearchScheme>
void OSS_HammingDistance(benchmark::State& state, uint8_t const maxErrors, TSearchScheme scheme)
{
    typedef HammingDistance TDistanceTag;
    OSSContext ossContext;

    TIter it(fm_index);

    uint64_t hitsNbr, uniqueHits;

    auto delegate = [&hitsNbr](OSSContext & ossContext, auto const & it, DnaString const & /*needle*/, uint32_t const /*needleId*/, uint8_t /*errors*/, bool const /*rev*/)
    {
        ++hitsNbr;
        unsigned x = 0;
        for (unsigned i = 0; i < length(getOccurrences(it)); ++i)
            x += getOccurrences(it)[i].i2;
    };
    auto delegateDirect = [&hitsNbr](OSSContext & ossContext, Pair<uint16_t, uint32_t> const & pos, Pair<uint16_t, uint32_t> const & posEnd, DnaString const & needle, uint32_t const needleId, uint8_t const errors)
    {
        ++hitsNbr;
        unsigned x = pos.i2;
    };

//     auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;

    for (auto _ : state)
    {
        hitsNbr = 0;
        uniqueHits = 0;

        std::vector<std::pair<TBitvector, TSupport>> empty_bitvectors;
        calcConstParameters(scheme);

        for (unsigned i = 0; i < length(reads); ++i)
        {
            uint64_t oldHits = hitsNbr;
//             find(0, maxErrors, myOSSContext, delegate, delegateDirect, fm_index, reads[i], HammingDistance());

            find(ossContext, delegate, delegateDirect, it, empty_bitvectors, scheme, reads[i], i, TDistanceTag());

            reverseComplement(reads[i]);

            find(ossContext, delegate, delegateDirect, it, empty_bitvectors, scheme, reads[i], i, TDistanceTag());

            benchmark::DoNotOptimize(uniqueHits += oldHits != hitsNbr);
        }
        // std::cout << "Backtracking: " << ((double)((time*100)/CLOCKS_PER_SEC)/100) << " s. "
        //           << "Hits: " << uniqueHits << " (" << hitsNbr << ")" << std::endl;
    }
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
            std:cout << i << "\n";
            uint64_t oldHits = hitsNbr;
            _optimalSearchScheme(delegate, it, reads[i], scheme, TDistanceTag());
            reverseComplement(reads[i]);
            _optimalSearchScheme(delegate, it, reads[i], scheme, TDistanceTag());
            benchmark::DoNotOptimize(uniqueHits += oldHits != hitsNbr);
        }
        std::cout << "Hits SS: " << uniqueHits << " (" << hitsNbr << ")" << std::endl; // Hits compare
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


BENCHMARK_CAPTURE(OSS_HammingDistance, OSS               , (uint8_t)1, OptimalSearchSchemes<0, 1>::VALUE)->Unit(benchmark::kMillisecond);
/*
BENCHMARK_CAPTURE(BM_HammingDistance, errors_1_backtracking      , (uint8_t)1)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BM_HammingDistance, errors_1_pig               , PigeonholeOptimumSearchSchemes<1>::VALUE)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BM_HammingDistance, errors_1_oss_parts_k_plus_1, PaperOptimumSearchSchemes<1>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
BENCHMARK_CAPTURE(BM_HammingDistance, errors_1_oss_parts_k_plus_2, PaperOptimumSearchSchemes<1>::VALUE_plus_two)->Unit(benchmark::kMillisecond);*/
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_1_oss_parts_k_plus_3, PaperOptimumSearchSchemes<1>::VALUE_plus_three)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_1_top               , OptimalSearchSchemes<0, 1>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_1_010_ss            , VrolandOptimumSearchSchemes<1>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds       , errors_1_010_jan_uni       , (uint8_t)1, false, predictify_unidirectional)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds       , errors_1_010_jan_bi        , (uint8_t)1, false, predictify_bidirectional)->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_2_backtracking      , (uint8_t)2)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_2_pig               , PigeonholeOptimumSearchSchemes<2>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_2_oss_parts_k_plus_1, PaperOptimumSearchSchemes<2>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_2_oss_parts_k_plus_2, PaperOptimumSearchSchemes<2>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_2_oss_parts_k_plus_3, PaperOptimumSearchSchemes<2>::VALUE_plus_three)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_2_kuc_parts_k_plus_1, KucherovOptimumSearchSchemes<2>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_2_kuc_parts_k_plus_2, KucherovOptimumSearchSchemes<2>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_2_top               , OptimalSearchSchemes<0, 2>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_2_010_ss            , VrolandOptimumSearchSchemes<2>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds       , errors_2_010_jan_uni       , (uint8_t)2, false, predictify_unidirectional)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds       , errors_2_010_jan_bi        , (uint8_t)2, false, predictify_bidirectional)->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_3_backtracking      , (uint8_t)3)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_3_pig               , PigeonholeOptimumSearchSchemes<3>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_3_oss_parts_k_plus_1, PaperOptimumSearchSchemes<3>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_3_oss_parts_k_plus_2, PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_3_oss_parts_k_plus_3, PaperOptimumSearchSchemes<3>::VALUE_plus_three)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_3_kuc_parts_k_plus_1, KucherovOptimumSearchSchemes<3>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_3_kuc_parts_k_plus_2, KucherovOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_3_top               , OptimalSearchSchemes<0, 3>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_3_010_ss            , VrolandOptimumSearchSchemes<3>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds       , errors_3_010_jan_uni       , (uint8_t)3, false, predictify_unidirectional)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds       , errors_3_010_jan_bi        , (uint8_t)3, false, predictify_bidirectional)->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_4_backtracking      , (uint8_t)4)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_4_pig               , PigeonholeOptimumSearchSchemes<4>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_4_oss_parts_k_plus_1, PaperOptimumSearchSchemes<4>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_4_oss_parts_k_plus_2, PaperOptimumSearchSchemes<4>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_4_oss_parts_k_plus_3, PaperOptimumSearchSchemes<4>::VALUE_plus_three)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_4_kuc_parts_k_plus_1, KucherovOptimumSearchSchemes<4>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_4_kuc_parts_k_plus_2, KucherovOptimumSearchSchemes<4>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_4_top               , OptimalSearchSchemes<0, 4>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_HammingDistance, errors_4_010_ss            , VrolandOptimumSearchSchemes<4>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds       , errors_4_010_jan_uni       , (uint8_t)4, false, predictify_unidirectional)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds       , errors_4_010_jan_bi        , (uint8_t)4, false, predictify_bidirectional)->Unit(benchmark::kMillisecond);
//
//
//
// BENCHMARK_CAPTURE(BM_EditDistance, errors_1_backtracking      , (uint8_t)1)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_1_pig               , PigeonholeOptimumSearchSchemes<1>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_1_oss_parts_k_plus_1, PaperOptimumSearchSchemes<1>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_1_oss_parts_k_plus_2, PaperOptimumSearchSchemes<1>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_1_oss_parts_k_plus_3, PaperOptimumSearchSchemes<1>::VALUE_plus_three)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_1_top               , OptimalSearchSchemes<0, 1>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_1_010_ss            , VrolandOptimumSearchSchemes<1>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds    , errors_1_010_jan_uni       , (uint8_t)1, true, predictify_unidirectional)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds    , errors_1_010_jan_bi        , (uint8_t)1, true, predictify_bidirectional)->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(BM_EditDistance, errors_2_backtracking      , (uint8_t)2)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_2_pig               , PigeonholeOptimumSearchSchemes<2>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_2_oss_parts_k_plus_1, PaperOptimumSearchSchemes<2>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_2_oss_parts_k_plus_2, PaperOptimumSearchSchemes<2>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_2_oss_parts_k_plus_3, PaperOptimumSearchSchemes<2>::VALUE_plus_three)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_2_kuc_parts_k_plus_1, KucherovOptimumSearchSchemes<2>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_2_kuc_parts_k_plus_2, KucherovOptimumSearchSchemes<2>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_2_top               , OptimalSearchSchemes<0, 2>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_2_010_ss            , VrolandOptimumSearchSchemes<2>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds    , errors_2_010_jan_uni       , (uint8_t)2, true, predictify_unidirectional)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds    , errors_2_010_jan_bi        , (uint8_t)2, true, predictify_bidirectional)->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(BM_EditDistance, errors_3_backtracking      , (uint8_t)3)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_3_pig               , PigeonholeOptimumSearchSchemes<3>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_3_oss_parts_k_plus_1, PaperOptimumSearchSchemes<3>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_3_oss_parts_k_plus_2, PaperOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_3_oss_parts_k_plus_1, PaperOptimumSearchSchemes<3>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_3_kuc_parts_k_plus_2, KucherovOptimumSearchSchemes<3>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_3_kuc_parts_k_plus_3, KucherovOptimumSearchSchemes<3>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_3_top               , OptimalSearchSchemes<0, 3>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_3_010_ss            , VrolandOptimumSearchSchemes<3>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds    , errors_3_010_jan_uni       , (uint8_t)3, true, predictify_unidirectional)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds    , errors_3_010_jan_bi        , (uint8_t)3, true, predictify_bidirectional)->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(BM_EditDistance, errors_4_backtracking      , (uint8_t)4)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_4_pig               , PigeonholeOptimumSearchSchemes<4>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_4_oss_parts_k_plus_1, PaperOptimumSearchSchemes<4>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_4_oss_parts_k_plus_2, PaperOptimumSearchSchemes<4>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_4_oss_parts_k_plus_1, PaperOptimumSearchSchemes<4>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_4_kuc_parts_k_plus_2, KucherovOptimumSearchSchemes<4>::VALUE_plus_one)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_4_kuc_parts_k_plus_3, KucherovOptimumSearchSchemes<4>::VALUE_plus_two)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_4_top               , OptimalSearchSchemes<0, 4>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_EditDistance, errors_4_010_ss            , VrolandOptimumSearchSchemes<4>::VALUE)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds    , errors_4_010_jan_uni       , (uint8_t)4, true, predictify_unidirectional)->Unit(benchmark::kMillisecond);
// BENCHMARK_CAPTURE(BM_010Seeds    , errors_4_010_jan_bi        , (uint8_t)4, true, predictify_bidirectional)->Unit(benchmark::kMillisecond);

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
