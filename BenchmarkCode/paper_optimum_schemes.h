#include <seqan/index.h>

using namespace seqan;

template <size_t max, typename TVoidType = void>
struct BacktrackingSchemes
{
    static constexpr std::array<OptimalSearch<1>, 1> VALUE
    {{
        { {{1}}, {{0}}, {{max}}, {{0}}, 0 }
    }};
};

template <size_t max, typename TVoidType>
constexpr std::array<OptimalSearch<1>, 1> BacktrackingSchemes<max, TVoidType>::VALUE;

// Contains Optimum Search Schemes from the paper for errors K = 1, ..., 4 and parts P = K+1, ... K+3

template <size_t max, typename TVoidType = void>
struct PaperOptimumSearchSchemes;

template <typename TVoidType>
struct PaperOptimumSearchSchemes<1, TVoidType>
{
    static constexpr std::array<OptimalSearch<2>, 2> VALUE_plus_one
    {{
        { {{1, 2}}, {{0, 0}}, {{0, 1}}, {{0, 0}}, 0 },
        { {{2, 1}}, {{0, 1}}, {{0, 1}}, {{0, 0}}, 0 }
    }};
    static constexpr std::array<OptimalSearch<3>, 2> VALUE_plus_two
    {{
        { {{1, 2, 3}}, {{0, 0, 1}}, {{0, 0, 1}}, {{0, 0, 0}}, 0},
        { {{3, 2, 1}}, {{0, 0, 0}}, {{0, 1, 1}}, {{0, 0, 0}}, 0}
    }};
    static constexpr std::array<OptimalSearch<4>, 2> VALUE_plus_three
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 0, 0}}, {{0, 0, 1, 1}}, {{0, 0, 0, 0}}, 0},
        { {{4, 3, 2, 1}}, {{0, 0, 0, 1}}, {{0, 0, 1, 1}}, {{0, 0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<2>, 2> PaperOptimumSearchSchemes<1, TVoidType>::VALUE_plus_one;

template <typename TVoidType>
constexpr std::array<OptimalSearch<3>, 2> PaperOptimumSearchSchemes<1, TVoidType>::VALUE_plus_two;

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 2> PaperOptimumSearchSchemes<1, TVoidType>::VALUE_plus_three;

template <typename TVoidType>
struct PaperOptimumSearchSchemes<2, TVoidType>
{
    static constexpr std::array<OptimalSearch<3>, 3> VALUE_plus_one
    {{
        { {{1, 2, 3}}, {{0, 0, 2}}, {{0, 1, 2}}, {{0, 0, 0}}, 0},
        { {{3, 2, 1}}, {{0, 0, 0}}, {{0, 2, 2}}, {{0, 0, 0}}, 0},
        { {{2, 3, 1}}, {{0, 1, 1}}, {{0, 1, 2}}, {{0, 0, 0}}, 0}
    }};
    static constexpr std::array<OptimalSearch<4>, 3> VALUE_plus_two
    {{
        { {{2, 1, 3, 4}}, {{0, 0, 1, 1}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0},
        { {{3, 2, 1, 4}}, {{0, 0, 0, 0}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, 0},
        { {{4, 3, 2, 1}}, {{0, 0, 0, 2}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, 0}
    }};
    static constexpr std::array<OptimalSearch<5>, 3> VALUE_plus_three
    {{
        { {{2, 1, 3, 4, 5}}, {{0, 0, 0, 1, 1}}, {{0, 0, 2, 2, 2}}, {{0, 0, 0, 0, 0}}, 0},
        { {{4, 3, 2, 1, 5}}, {{0, 0, 0, 0, 0}}, {{0, 0, 1, 1, 2}}, {{0, 0, 0, 0, 0}}, 0},
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 2}}, {{0, 1, 1, 2, 2}}, {{0, 0, 0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<3>, 3> PaperOptimumSearchSchemes<2, TVoidType>::VALUE_plus_one;

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 3> PaperOptimumSearchSchemes<2, TVoidType>::VALUE_plus_two;

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> PaperOptimumSearchSchemes<2, TVoidType>::VALUE_plus_three;

template <typename TVoidType>
struct PaperOptimumSearchSchemes<3, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 3> VALUE_plus_one
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 0, 3}}, {{0, 2, 3, 3}}, {{0, 0, 0, 0}}, 0},
        { {{2, 3, 4, 1}}, {{0, 0, 0, 0}}, {{1, 2, 2, 3}}, {{0, 0, 0, 0}}, 0},
        { {{3, 4, 2, 1}}, {{0, 0, 2, 2}}, {{0, 0, 3, 3}}, {{0, 0, 0, 0}}, 0}
    }};
    static constexpr std::array<OptimalSearch<5>, 3> VALUE_plus_two
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 2, 2}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0},
        { {{4, 3, 2, 1, 5}}, {{0, 0, 0, 0, 0}}, {{1, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0},
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 3}}, {{0, 2, 2, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}
        // { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 3}}, {{0, 2, 2, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        // { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 2, 2}}, {{0, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        // { {{3, 4, 5, 2, 1}}, {{0, 0, 1, 1, 1}}, {{0, 1, 1, 2, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        // { {{4, 5, 3, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
    static constexpr std::array<OptimalSearch<6>, 3> VALUE_plus_three
    {{
        { {{1, 2, 3, 4, 5, 6}}, {{0, 0, 0, 0, 0, 3}}, {{0, 2, 2, 2, 3, 3}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{2, 3, 4, 5, 6, 1}}, {{0, 0, 0, 0, 0, 0}}, {{1, 1, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{6, 5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 2, 2}}, {{0, 0, 3, 3, 3, 3}}, {{0, 0, 0, 0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 3> PaperOptimumSearchSchemes<3, TVoidType>::VALUE_plus_one;

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> PaperOptimumSearchSchemes<3, TVoidType>::VALUE_plus_two;

template <typename TVoidType>
constexpr std::array<OptimalSearch<6>, 3> PaperOptimumSearchSchemes<3, TVoidType>::VALUE_plus_three;

template <typename TVoidType>
struct PaperOptimumSearchSchemes<4, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 3> VALUE_plus_one
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 4}}, {{0, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 0}}, {{2, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 3, 3}}, {{0, 0, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0}
    }};
    static constexpr std::array<OptimalSearch<6>, 3> VALUE_plus_two
    {{
        { {{1, 2, 3, 4, 5, 6}}, {{0, 0, 0, 0, 0, 4}}, {{0, 3, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{2, 3, 4, 5, 6, 1}}, {{0, 0, 0, 0, 0, 0}}, {{2, 2, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{6, 5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 3, 3}}, {{0, 0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}
    }};
    static constexpr std::array<OptimalSearch<7>, 3> VALUE_plus_three
    {{
        { {{1, 2, 3, 4, 5, 6, 7}}, {{0, 1, 1, 1, 1, 1, 1}}, {{3, 3, 3, 3, 3, 3, 4}}, {{0, 0, 0, 0, 0, 0, 0}}, 0},
        { {{1, 2, 3, 4, 5, 6, 7}}, {{0, 0, 0, 0, 0, 0, 0}}, {{0, 0, 4, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0, 0}}, 0},
        { {{7, 6, 5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 0, 0, 4}}, {{0, 3, 3, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 3> PaperOptimumSearchSchemes<4, TVoidType>::VALUE_plus_one;

template <typename TVoidType>
constexpr std::array<OptimalSearch<6>, 3> PaperOptimumSearchSchemes<4, TVoidType>::VALUE_plus_two;

template <typename TVoidType>
constexpr std::array<OptimalSearch<7>, 3> PaperOptimumSearchSchemes<4, TVoidType>::VALUE_plus_three;

///////////////////////////////

template <size_t max, typename TVoidType = void>
struct KucherovOptimumSearchSchemes;

template <typename TVoidType>
struct KucherovOptimumSearchSchemes<2, TVoidType>
{
    static constexpr std::array<OptimalSearch<3>, 3> VALUE_plus_one
    {{
        { {{1, 2, 3}}, {{0, 0, 0}}, {{0, 2, 2}}, {{0, 0, 0}}, 0},
        { {{3, 2, 1}}, {{0, 0, 0}}, {{0, 1, 2}}, {{0, 0, 0}}, 0},
        { {{2, 1, 3}}, {{0, 0, 1}}, {{0, 1, 2}}, {{0, 0, 0}}, 0}
    }};
    static constexpr std::array<OptimalSearch<4>, 4> VALUE_plus_two
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 0, 0}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, 0},
        { {{4, 3, 2, 1}}, {{0, 0, 0, 0}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, 0},
        { {{2, 3, 4, 1}}, {{0, 0, 0, 1}}, {{0, 0, 1, 2}}, {{0, 0, 0, 0}}, 0},
        { {{1, 2, 3, 4}}, {{0, 0, 0, 2}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<3>, 3> KucherovOptimumSearchSchemes<2, TVoidType>::VALUE_plus_one;

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 4> KucherovOptimumSearchSchemes<2, TVoidType>::VALUE_plus_two;



template <typename TVoidType>
struct KucherovOptimumSearchSchemes<3, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 4> VALUE_plus_one
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 0, 0}}, {{0, 1, 3, 3}}, {{0, 0, 0, 0}}, 0},
        { {{2, 1, 3, 4}}, {{0, 0, 1, 1}}, {{0, 1, 3, 3}}, {{0, 0, 0, 0}}, 0},
        { {{3, 4, 2, 1}}, {{0, 0, 0, 0}}, {{0, 1, 3, 3}}, {{0, 0, 0, 0}}, 0},
        { {{4, 3, 2, 1}}, {{0, 0, 1, 1}}, {{0, 1, 3, 3}}, {{0, 0, 0, 0}}, 0}
    }};
    static constexpr std::array<OptimalSearch<5>, 4> VALUE_plus_two
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 0}}, {{0, 1, 2, 3, 3}}, {{0, 0, 0, 0, 0}}, 0},
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 0}}, {{0, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0},
        { {{3, 4, 5, 2, 1}}, {{0, 0, 0, 0, 1}}, {{0, 1, 1, 3, 3}}, {{0, 0, 0, 0, 0}}, 0},
        { {{4, 5, 3, 2, 1}}, {{0, 0, 0, 1, 2}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 4> KucherovOptimumSearchSchemes<3, TVoidType>::VALUE_plus_one;

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 4> KucherovOptimumSearchSchemes<3, TVoidType>::VALUE_plus_two;



template <typename TVoidType>
struct KucherovOptimumSearchSchemes<4, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 8> VALUE_plus_one
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 0}}, {{0, 2, 2, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 1, 3, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{2, 1, 3, 4, 5}}, {{0, 0, 1, 3, 3}}, {{0, 1, 3, 3, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{1, 2, 3, 4, 5}}, {{0, 0, 1, 3, 3}}, {{0, 1, 3, 3, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{4, 3, 5, 2, 1}}, {{0, 0, 0, 1, 1}}, {{0, 1, 2, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{3, 2, 1, 4, 5}}, {{0, 0, 0, 1, 3}}, {{0, 1, 2, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{2, 1, 3, 4, 5}}, {{0, 0, 1, 2, 4}}, {{0, 1, 2, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 3, 4}}, {{0, 0, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0}
    }};

    static constexpr std::array<OptimalSearch<6>, 10> VALUE_plus_two
    {{
        { {{1, 2, 3, 4, 5, 6}}, {{0, 0, 0, 0, 0, 0,}}, {{0, 1, 2, 3, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{2, 3, 4, 5, 6, 1}}, {{0, 0, 0, 0, 0, 0,}}, {{0, 1, 2, 3, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{6, 5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 0, 1,}}, {{0, 1, 2, 2, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{4, 5, 6, 3, 2, 1}}, {{0, 0, 0, 0, 1, 2,}}, {{0, 1, 1, 3, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{3, 4, 5, 6, 2, 1}}, {{0, 0, 0, 0, 2, 3,}}, {{0, 1, 1, 2, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{5, 6, 4, 3, 2, 1}}, {{0, 0, 0, 1, 3, 3,}}, {{0, 0, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{1, 2, 3, 4, 5, 6}}, {{0, 0, 0, 3, 3, 3,}}, {{0, 0, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{1, 2, 3, 4, 5, 6}}, {{0, 0, 0, 0, 4, 4,}}, {{0, 0, 2, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{3, 4, 2, 1, 5, 6}}, {{0, 0, 0, 1, 2, 4,}}, {{0, 0, 2, 2, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0},
        { {{5, 6, 4, 3, 2, 1}}, {{0, 0, 0, 0, 4, 4,}}, {{0, 0, 1, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 8> KucherovOptimumSearchSchemes<4, TVoidType>::VALUE_plus_one;

template <typename TVoidType>
constexpr std::array<OptimalSearch<6>, 10> KucherovOptimumSearchSchemes<4, TVoidType>::VALUE_plus_two;



template <typename TVoidType>
struct KucherovOptimumSearchSchemes<1, TVoidType>
{
    static constexpr std::array<OptimalSearch<0>, 0> VALUE_plus_one
    {{
    }};
    static constexpr std::array<OptimalSearch<0>, 0> VALUE_plus_two
    {{
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<0>, 0> KucherovOptimumSearchSchemes<1, TVoidType>::VALUE_plus_one;

template <typename TVoidType>
constexpr std::array<OptimalSearch<0>, 0> KucherovOptimumSearchSchemes<1, TVoidType>::VALUE_plus_two;

///////////////////////////////

template <size_t max, typename TVoidType = void>
struct PigeonholeOptimumSearchSchemes;

template <typename TVoidType>
struct PigeonholeOptimumSearchSchemes<1, TVoidType>
{
    static constexpr std::array<OptimalSearch<2>, 2> VALUE
    {{
        { {{1, 2}}, {{0, 0}}, {{0, 1}}, {{0, 0}}, 0},
        { {{2, 1}}, {{0, 0}}, {{0, 1}}, {{0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<2>, 2> PigeonholeOptimumSearchSchemes<1, TVoidType>::VALUE;



template <typename TVoidType>
struct PigeonholeOptimumSearchSchemes<2, TVoidType>
{
    static constexpr std::array<OptimalSearch<3>, 3> VALUE
    {{
        { {{1, 2, 3}}, {{0, 0, 0}}, {{0, 2, 2}}, {{0, 0, 0}}, 0},
        { {{2, 3, 1}}, {{0, 0, 0}}, {{0, 2, 2}}, {{0, 0, 0}}, 0},
        { {{3, 2, 1}}, {{0, 0, 0}}, {{0, 2, 2}}, {{0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<3>, 3> PigeonholeOptimumSearchSchemes<2, TVoidType>::VALUE;



template <typename TVoidType>
struct PigeonholeOptimumSearchSchemes<3, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 4> VALUE
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 0, 0}}, {{0, 3, 3, 3}}, {{0, 0, 0, 0}}, 0},
        { {{2, 3, 4, 1}}, {{0, 0, 0, 0}}, {{0, 3, 3, 3}}, {{0, 0, 0, 0}}, 0},
        { {{3, 4, 2, 1}}, {{0, 0, 0, 0}}, {{0, 3, 3, 3}}, {{0, 0, 0, 0}}, 0},
        { {{4, 3, 2, 1}}, {{0, 0, 0, 0}}, {{0, 3, 3, 3}}, {{0, 0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 4> PigeonholeOptimumSearchSchemes<3, TVoidType>::VALUE;



template <typename TVoidType>
struct PigeonholeOptimumSearchSchemes<4, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 5> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 0}}, {{0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 0}}, {{0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{3, 4, 5, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{4, 5, 3, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0},
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0}}, 0}
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 5> PigeonholeOptimumSearchSchemes<4, TVoidType>::VALUE;

///////////////////////////////

template <size_t max, typename TVoidType = void>
struct VrolandOptimumSearchSchemes;

template <typename TVoidType>
struct VrolandOptimumSearchSchemes<1, TVoidType>
{
    static constexpr std::array<OptimalSearch<3>, 3> VALUE
    {{
        { {{1, 2, 3}}, {{0, 0, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, 0}, // 00*
        { {{1, 2, 3}}, {{0, 1, 1}}, {{0, 1, 1}}, {{0, 0, 0}}, 0}, // 010
        { {{2, 3, 1}}, {{0, 0, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, 0}  // *00
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<3>, 3> VrolandOptimumSearchSchemes<1, TVoidType>::VALUE;



template <typename TVoidType>
struct VrolandOptimumSearchSchemes<2, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 6> VALUE
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 0, 0}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0}, // 00**
        { {{1, 2, 3, 4}}, {{0, 1, 1, 1}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, 0}, // 010*
        { {{1, 2, 3, 4}}, {{0, 1, 2, 2}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, 0}, // 0110

        { {{2, 3, 4, 1}}, {{0, 0, 0, 0}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0}, // *00*
        { {{2, 3, 4, 1}}, {{0, 1, 1, 1}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, 0}, // *010

        { {{3, 4, 2, 1}}, {{0, 0, 0, 0}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0}  // **00
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 6> VrolandOptimumSearchSchemes<2, TVoidType>::VALUE;



template <typename TVoidType>
struct VrolandOptimumSearchSchemes<3, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 10> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 0}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // 00***
        { {{1, 2, 3, 4, 5}}, {{0, 1, 1, 1, 1}}, {{0, 1, 1, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // 010**
        { {{1, 2, 3, 4, 5}}, {{0, 1, 2, 2, 2}}, {{0, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0}, // 0110*
        { {{1, 2, 3, 4, 5}}, {{0, 1, 2, 3, 3}}, {{0, 1, 2, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // 01110

        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 0}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // *00**
        { {{2, 3, 4, 5, 1}}, {{0, 1, 1, 1, 1}}, {{0, 1, 1, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // *010*
        { {{2, 3, 4, 5, 1}}, {{0, 1, 2, 2, 2}}, {{0, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0}, // *0110

        { {{3, 4, 5, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // **00*
        { {{3, 4, 5, 2, 1}}, {{0, 1, 1, 1, 1}}, {{0, 1, 1, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // **010

        { {{4, 5, 3, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}  // ***00
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 10> VrolandOptimumSearchSchemes<3, TVoidType>::VALUE;



template <typename TVoidType>
struct VrolandOptimumSearchSchemes<4, TVoidType>
{
    static constexpr std::array<OptimalSearch<6>, 15> VALUE
    {{
        { {{1, 2, 3, 4, 5, 6}}, {{0, 0, 0, 0, 0, 0}}, {{0, 0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // 00****
        { {{1, 2, 3, 4, 5, 6}}, {{0, 1, 1, 1, 1, 1}}, {{0, 1, 1, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // 010***
        { {{1, 2, 3, 4, 5, 6}}, {{0, 1, 2, 2, 2, 2}}, {{0, 1, 2, 2, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // 0110**
        { {{1, 2, 3, 4, 5, 6}}, {{0, 1, 2, 3, 3, 3}}, {{0, 1, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // 01110*
        { {{1, 2, 3, 4, 5, 6}}, {{0, 1, 2, 3, 4, 4}}, {{0, 1, 2, 3, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // 011110

        { {{2, 3, 4, 5, 6, 1}}, {{0, 0, 0, 0, 0, 0}}, {{0, 0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // *00***
        { {{2, 3, 4, 5, 6, 1}}, {{0, 1, 1, 1, 1, 1}}, {{0, 1, 1, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // *010**
        { {{2, 3, 4, 5, 6, 1}}, {{0, 1, 2, 2, 2, 2}}, {{0, 1, 2, 2, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // *0110*
        { {{2, 3, 4, 5, 6, 1}}, {{0, 1, 2, 3, 3, 3}}, {{0, 1, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // *01110

        { {{3, 4, 5, 6, 2, 1}}, {{0, 0, 0, 0, 0, 0}}, {{0, 0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // **00**
        { {{3, 4, 5, 6, 2, 1}}, {{0, 1, 1, 1, 1, 1}}, {{0, 1, 1, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // **010*
        { {{3, 4, 5, 6, 2, 1}}, {{0, 1, 2, 2, 2, 2}}, {{0, 1, 2, 2, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // **0110

        { {{4, 5, 6, 3, 2, 1}}, {{0, 0, 0, 0, 0, 0}}, {{0, 0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // ***00*
        { {{4, 5, 6, 3, 2, 1}}, {{0, 1, 1, 1, 1, 1}}, {{0, 1, 1, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // ***010

        { {{5, 6, 4, 3, 2, 1}}, {{0, 0, 0, 0, 0, 0}}, {{0, 0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}  // ****00
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<6>, 15> VrolandOptimumSearchSchemes<4, TVoidType>::VALUE;


///////////////////////////////

template <size_t max, typename TVoidType = void>
struct VrolandMergedOptimumSearchSchemes;

template <typename TVoidType>
struct VrolandMergedOptimumSearchSchemes<1, TVoidType>
{
    static constexpr std::array<OptimalSearch<3>, 2> VALUE
    {{
        { {{1, 2, 3}}, {{0, 0, 0}}, {{0, 1, 1}}, {{0, 0, 0}}, 0}, // 00*, 010
        { {{2, 3, 1}}, {{0, 0, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, 0}  // *00
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<3>, 2> VrolandMergedOptimumSearchSchemes<1, TVoidType>::VALUE;



template <typename TVoidType>
struct VrolandMergedOptimumSearchSchemes<2, TVoidType>
{
    static constexpr std::array<OptimalSearch<4>, 3> VALUE
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 0, 0}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, 0}, // 00**,  010*, 0110
        { {{2, 3, 4, 1}}, {{0, 0, 0, 0}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, 0}, // *00*, *010
        { {{3, 4, 2, 1}}, {{0, 0, 0, 0}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0}  // **00
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<4>, 3> VrolandMergedOptimumSearchSchemes<2, TVoidType>::VALUE;



template <typename TVoidType>
struct VrolandMergedOptimumSearchSchemes<3, TVoidType>
{
    static constexpr std::array<OptimalSearch<5>, 4> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 0}}, {{0, 1, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // 00***, 010**, 0110*, 01110
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 0, 0}}, {{0, 1, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // *00**, *010*, *0110
        { {{3, 4, 5, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 1, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}, // **00*, **010
        { {{4, 5, 3, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0}  // ***00
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<5>, 4> VrolandMergedOptimumSearchSchemes<3, TVoidType>::VALUE;



template <typename TVoidType>
struct VrolandMergedOptimumSearchSchemes<4, TVoidType>
{
    static constexpr std::array<OptimalSearch<6>, 5> VALUE
    {{
        { {{1, 2, 3, 4, 5, 6}}, {{0, 0, 0, 0, 0, 0}}, {{0, 1, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // 00****, 010***, 0110**, 01110*, 011110
        { {{2, 3, 4, 5, 6, 1}}, {{0, 0, 0, 0, 0, 0}}, {{0, 1, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // *00***, *010**, *0110*, *01110
        { {{3, 4, 5, 6, 2, 1}}, {{0, 0, 0, 0, 0, 0}}, {{0, 1, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // **00**, **010*, **0110
        { {{4, 5, 6, 3, 2, 1}}, {{0, 0, 0, 0, 0, 0}}, {{0, 1, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}, // ***00*, ***010
        { {{5, 6, 4, 3, 2, 1}}, {{0, 0, 0, 0, 0, 0}}, {{0, 0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0}  // ****00
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearch<6>, 5> VrolandMergedOptimumSearchSchemes<4, TVoidType>::VALUE;
