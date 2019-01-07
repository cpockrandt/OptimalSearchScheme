#include <seqan/arg_parse.h>

#include "common.h"
#include "paper_optimum_schemes.h"

using namespace seqan;

template <uint8_t errors, typename TDistanceTag>
void output_count(unsigned const read_length, TDistanceTag /**/)
{
    uint64_t trivial, pigeonhole, oss1, oss2, oss3, top, kucherov1, kucherov2, vroland;

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            trivial = countTrivialSearch<Dna>(read_length, errors, TDistanceTag());
        }

        #pragma omp section
        {
            auto pigScheme = PigeonholeOptimumSearchSchemes<errors>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(pigScheme, read_length);
            pigeonhole = countSearchScheme<Dna>(read_length, pigScheme, TDistanceTag());
            // check_OSS(pigScheme);
        }

        #pragma omp section
        {
            auto vorlandScheme = VrolandOptimumSearchSchemes<errors>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(vorlandScheme, read_length);
            vroland = countSearchScheme<Dna>(read_length, vorlandScheme, TDistanceTag());
            // check_OSS(vorlandScheme);
        }

        #pragma omp section
        {
            auto topScheme = OptimalSearchSchemes<0, errors>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(topScheme, read_length);
            top = countSearchScheme<Dna>(read_length, topScheme, TDistanceTag());
            // check_OSS(topScheme);
        }

        #pragma omp section
        {
            // K+1 parts
            auto scheme1 = PaperOptimumSearchSchemes<errors>::VALUE_plus_one;
            _optimalSearchSchemeComputeFixedBlocklength(scheme1, read_length);
            oss1 = countSearchScheme<Dna>(read_length, scheme1, TDistanceTag());
            // check_OSS(scheme1);
        }

        #pragma omp section
        {
            // K+2 parts
            auto scheme2 = PaperOptimumSearchSchemes<errors>::VALUE_plus_two;
            _optimalSearchSchemeComputeFixedBlocklength(scheme2, read_length);
            oss2 = countSearchScheme<Dna>(read_length, scheme2, TDistanceTag());
            // check_OSS(scheme2);
        }

        #pragma omp section
        {
            // K+3 parts
            auto scheme3 = PaperOptimumSearchSchemes<errors>::VALUE_plus_three;
            _optimalSearchSchemeComputeFixedBlocklength(scheme3, read_length);
            oss3 = countSearchScheme<Dna>(read_length, scheme3, TDistanceTag());
            // check_OSS(scheme3);
        }

        #pragma omp section
        {
            // K+1 parts
            if (errors > 1)
            {
                auto kscheme1 = KucherovOptimumSearchSchemes<errors>::VALUE_plus_one;
                _optimalSearchSchemeComputeFixedBlocklength(kscheme1, read_length);
                kucherov1 = countSearchScheme<Dna>(read_length, kscheme1, TDistanceTag());
                // TODO: check_OSS(kscheme1);
            }
        }

        #pragma omp section
        {
            // K+2 parts
            if (errors > 1)
            {
                auto kscheme2 = KucherovOptimumSearchSchemes<errors>::VALUE_plus_two;
                _optimalSearchSchemeComputeFixedBlocklength(kscheme2, read_length);
                kucherov2 = countSearchScheme<Dna>(read_length, kscheme2, TDistanceTag());
                // TODO: check_OSS(kscheme2);
            }
        }
    }

    std::cout << "Tri, K = " << (unsigned)errors << ", P = ***: " << trivial  << '\n';
    std::cout << "Pig, K = " << (unsigned)errors << ", P = ***: " << pigeonhole  << '\n';
    std::cout << "010, K = " << (unsigned)errors << ", P = ***: " << vroland  << '\n';
    std::cout << "OSS, K = " << (unsigned)errors << ", P = K+" << 1 << ": " << oss1 << '\n';
    std::cout << "OSS, K = " << (unsigned)errors << ", P = K+" << 2 << ": " << oss2 << '\n';
    std::cout << "OSS, K = " << (unsigned)errors << ", P = K+" << 3 << ": " << oss3 << '\n';
    std::cout << "TOP, K = " << (unsigned)errors << ", P = ***: " << top  << '\n';
    if (errors > 1)
    {
        std::cout << "Kuc, K = " << (unsigned)errors << ", P = K+" << 1 << ": " << kucherov1 << '\n';
        std::cout << "Kuc, K = " << (unsigned)errors << ", P = K+" << 2 << ": " << kucherov2 << '\n';
    }
    std::cout << "-------------------------------------------------\n";
}

int main(int argc, char *argv[])
{
    // Argument Parser
    ArgumentParser parser("Counting edges");
    addDescription(parser, "App for couting edges for trivial backtracking and Optimum Search Schemes.");

    addOption(parser, ArgParseOption("R", "read-length", "Length of the read", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "read-length", 101);

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    unsigned read_length;
    getOptionValue(read_length, parser, "read-length");

    std::cout << "*** Hamming Distance ***\n";
    output_count<1>(read_length, HammingDistance());
    output_count<2>(read_length, HammingDistance());
    output_count<3>(read_length, HammingDistance());
    output_count<4>(read_length, HammingDistance());
    // std::cout << "*** Edit Distance ***\n";
    // output_count<1>(read_length, EditDistance());
    // output_count<2>(read_length, EditDistance());
    // output_count<3>(read_length, EditDistance());
    // output_count<4>(read_length, EditDistance());

    return 0;
}
