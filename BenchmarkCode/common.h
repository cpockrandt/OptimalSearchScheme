#include <type_traits>

#include <seqan/index.h>

using namespace seqan;

namespace seqan {

template <typename TChar, typename TAlloc, typename TOwner>
struct SAValue<StringSet<String<TChar, TAlloc>, TOwner> >
{
    typedef Pair<uint16_t, uint32_t> Type;
};

template <typename TChar, typename TAlloc>
struct SAValue<String<TChar, TAlloc> >
{
    typedef uint32_t Type;
};

};

// Index type
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
typedef BidirectionalIndex<FMIndex<void, TMyFastConfig> > TIndexConfig;

template <typename TChar, typename TDistanceTag>
inline void _countTrivialSearch(unsigned const needleLength,
                                unsigned const needlePos,
                                uint8_t const errorsLeft, TDistanceTag const & /**/,
                                unsigned long long & counts)
{
    if (errorsLeft == 0 || needleLength == needlePos)
    {
        counts += needleLength - needlePos;
        return;
    }

    // Match / Mismatch
    unsigned long long countsMatch = 1, countsMismatchOrInsertion = 0;
    _countTrivialSearch<TChar>(needleLength, needlePos + 1, errorsLeft, TDistanceTag(), countsMatch);
    _countTrivialSearch<TChar>(needleLength, needlePos + 1, errorsLeft - 1, TDistanceTag(), countsMismatchOrInsertion);
    counts += countsMatch + (ValueSize<TChar>::VALUE - 1) * (countsMismatchOrInsertion + 1);

    if (std::is_same<TDistanceTag, EditDistance>::value)
    {
        // Deletion
        unsigned long long countsDeletion = 0; // TODO: do we count the deletion itself?
        _countTrivialSearch<TChar>(needleLength, needlePos, errorsLeft - 1, TDistanceTag(), countsDeletion);
        counts += ValueSize<TChar>::VALUE * countsDeletion;
        // Insertion
        counts += countsMismatchOrInsertion; // TODO: do we count the deletion itself?
    }
}

template <typename TChar, typename TDistanceTag>
inline unsigned long long countTrivialSearch(unsigned int needleLength, uint8_t const errors, TDistanceTag const & /**/)
{
    unsigned long long counts = 0;
    _countTrivialSearch<TChar>(needleLength, 0, errors, TDistanceTag(), counts);
    return counts;
}

// TODO: rename to needlePos!
template <typename TChar, typename TOptimumSearch, typename TDistanceTag>
inline void _countSearch(uint32_t const needleLength, uint32_t const needleLeftIt, uint32_t const needleRightIt,
                         uint8_t const errors, TOptimumSearch const & s,
                         bool const goToRight, uint8_t const blockIndex,
                         uint64_t & count, TDistanceTag const & /**/)
{
    // Done.
    if (needleLeftIt == 0 && needleRightIt == needleLength + 1)
        return;

    // Exact search in current block.
    uint8_t maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
    if (maxErrorsLeftInBlock == 0)
    {
        signed infixPosLeft, infixPosRight;

        if (goToRight)
        {
            infixPosLeft = needleRightIt + - 1;
            infixPosRight = needleLeftIt + s.blocklength[blockIndex] - 1;
        }
        else
        {
            infixPosLeft = needleRightIt - s.blocklength[blockIndex] - 1;
            infixPosRight = needleLeftIt - 1;
        }

        bool goToRight2 = s.pi[blockIndex + 1] > s.pi[blockIndex];
        count += infixPosRight - infixPosLeft + 1;
        if (goToRight)
            return _countSearch<TChar>(needleLength, needleLeftIt, infixPosRight + 2, errors, s, goToRight2, blockIndex + 1, count, TDistanceTag());
        else
            return _countSearch<TChar>(needleLength, infixPosLeft, needleRightIt, errors, s, goToRight2, blockIndex + 1, count, TDistanceTag());
    }
    // Approximate search in current block.
    else
    {
        unsigned charsLeft = s.blocklength[blockIndex] - (needleRightIt - needleLeftIt - 1);

        // Insertion
        if (std::is_same<TDistanceTag, EditDistance>::value)
        {
            uint8_t blockIndex2 = blockIndex;
            bool goToRight2 = goToRight;
            if (needleRightIt - needleLeftIt == s.blocklength[blockIndex])
            {
                ++blockIndex2;
                goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex];
            }

            _countSearch<TChar>(needleLength, needleLeftIt - !goToRight, needleRightIt + goToRight, errors + 1, s, goToRight2, blockIndex2, count, TDistanceTag());
        }

        for (unsigned i = 0; i < ValueSize<TChar>::VALUE; ++i)
        {
            bool delta = (i != 0);

            // Deletion
            if (std::is_same<TDistanceTag, EditDistance>::value)
                _countSearch<TChar>(needleLength, needleLeftIt, needleRightIt, errors + 1, s, goToRight, blockIndex, count, TDistanceTag());

            if (minErrorsLeftInBlock > 0 && charsLeft - 1 < minErrorsLeftInBlock - delta)
                continue;

            uint8_t blockIndex2 = blockIndex;
            bool goToRight2 = goToRight;
            if (needleRightIt - needleLeftIt == s.blocklength[blockIndex])
            {
                ++blockIndex2;
                goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex];
            }

            ++count;
            if (goToRight)
                _countSearch<TChar>(needleLength, needleLeftIt, needleRightIt + 1, errors + delta, s, goToRight2, blockIndex2, count, TDistanceTag());
            else
                _countSearch<TChar>(needleLength, needleLeftIt - 1, needleRightIt, errors + delta, s, goToRight2, blockIndex2, count, TDistanceTag());
        }
    }
}

template <typename TChar, typename TOptimumSearch, typename TDistanceTag>
inline uint64_t countSearch(unsigned int needleLength, TOptimumSearch const & s, TDistanceTag const & /**/)
{
    uint64_t count = 0; // TODO: removed initialDirection. correct?
    _countSearch<TChar>(needleLength, s.startPos, s.startPos + 1, 0, s, true /*TODO:initialDirection Rrev()*/, 0, count, TDistanceTag());
    return count;
}

template <typename TChar, typename TOptimumSearchScheme, typename TDistanceTag>
inline uint64_t countSearchScheme(unsigned const needleLength, TOptimumSearchScheme const & scheme, TDistanceTag const & /**/)
{
    uint64_t count = 0;
    for (auto const & s : scheme)
        count += countSearch<TChar>(needleLength, s, TDistanceTag());
    return count;
}

template <typename TChar, typename TOptimumSearchScheme>
inline void _optimalSearchSchemeComputeOptimalBlocklength(TOptimumSearchScheme scheme,
                                                          unsigned const maxErrors, uint32_t const needleLength)
{
    // count edges for all blocklengths and choose theoretically optimal one
    uint64_t countEdges, countEdgesOptimal = static_cast<uint64_t>(-1); // maximum value of unsigned int
    std::vector<uint32_t> blocklength, blocklengthOptimal;

    for (unsigned i = 1; i < scheme[0].pi.size(); ++i)
    {
        blocklength[i-1] = i * maxErrors;
    }
    back(blocklength) = needleLength;

    while (true)
    {
        std::vector<uint32_t> _blocklengthAbsolute = blocklength;
        for (unsigned i = _blocklengthAbsolute.size(); i > 0; --i)
        {
            _blocklengthAbsolute[i] -= _blocklengthAbsolute[i - 1];
        }

        _optimalSearchSchemeSetBlockLength(scheme, _blocklengthAbsolute);
        _optimalSearchSchemeInit(scheme);
        countEdges = countSearchScheme<TChar>(needleLength, scheme, HammingDistance()); // TODO: document that indels are always ignored!
        if (countEdges < countEdgesOptimal)
        {
            countEdgesOptimal = countEdges;
            blocklengthOptimal = _blocklengthAbsolute;
        }

        // compute next blockpos
        signed i;
        for (i = blocklength.size() - 2; i >= 0; --i)
        {
            // find rightmost element in blockpos that we can increment
            if (blocklength[i] < blocklength[i + 1] - maxErrors)
            {
                ++blocklength[i];
                // reset all elements after pos i
                for (unsigned j = i + 1; j < blocklength.size() - 1; ++j)
                {
                    blocklength[j] = blocklength[j - 1] + maxErrors;
                }
                break;
            }
        }
        if (i < 0)
            break;
    }

    _optimalSearchSchemeSetBlockLength(scheme, blocklengthOptimal);
    _optimalSearchSchemeInit(scheme);
}

typedef std::vector<bool> TBitvector;
typedef std::vector<bool> TSupport;

enum class ReturnCode {
	NOMAPPABILITY, DIRECTSEARCH, COMPMAPPABLE, ONEDIRECTION, MAPPABLE, FINISHED, UNIDIRECTIONAL, SUSPECTUNIDIRECTIONAL, FILTER, ERROR
};

template <size_t nbrBlocks, size_t N>
inline void calcConstParameters(std::array<OptimalSearch<nbrBlocks>, N> & ss)
{
    uint8_t id = 0;
    for (OptimalSearch<nbrBlocks> & s : ss){
        s.id = id;
        ++id;
//         int bsize = s.pi.size();
        uint8_t min = s.pi[0];
        uint8_t max = s.pi[0];
        // maybe < N?
        for(int i = 0; i < nbrBlocks; ++i){
            if(min > s.pi[i])
                min = s.pi[i];
            if(max < s.pi[i])
                max = s.pi[i];
            s.min[i] = min;
            s.max[i] = max;
        }
        uint8_t lastValue = s.pi[nbrBlocks - 1];
        int k = (nbrBlocks > 1) ? nbrBlocks - 2 : 0;
        while(k >= 0){
            if(s.pi[k] == lastValue - 1 || s.pi[k] == lastValue + 1)
            {
                lastValue = s.pi[k];
                --k;
            }else{
                s.startUniDir = k + 1;
                break;
            }
        }
    }
}

template <size_t nbrBlocks, size_t N>
inline void _optimalSearchSchemeComputeChronBlocklength(std::array<OptimalSearch<nbrBlocks>, N> & ss)
{
    for (OptimalSearch<nbrBlocks> & s : ss){
        s.chronBL[s.pi[0] - 1]  = s.blocklength[0];
        for(int j = 1; j < nbrBlocks; ++j)
            s.chronBL[s.pi[j] - 1] = s.blocklength[j] -  s.blocklength[j - 1];
        for(int j = 1; j < nbrBlocks; ++j)
            s.chronBL[j] += s.chronBL[j - 1];

        s.revChronBL[s.pi[nbrBlocks - 1] - 1]  = s.blocklength[nbrBlocks - 1] - s.blocklength[nbrBlocks - 2];
        for(int i = static_cast<int> (nbrBlocks) - 2; i >= 0; --i){
            s.revChronBL[s.pi[i] - 1] = s.blocklength[i] - ((i > 0) ? s.blocklength[i - 1] : 0);
        }
        for(int i = static_cast<int> (nbrBlocks) - 2; i >= 0; --i)
            s.revChronBL[i] += s.revChronBL[i + 1];
    }
    for (OptimalSearch<nbrBlocks> & s : ss){
        for (uint8_t j = 0; j < s.pi.size(); ++j)
        {
            s.blockStarts[j] = (s.pi[j] - 1 == 0) ? 0 : s.chronBL[s.pi[j] - 2];
            s.blockEnds[j] = s.chronBL[s.pi[j] - 1];
        }

        for(uint8_t j = 0; j < s.pi.size(); ++j){
            s.revblockStarts[j] = (s.pi[j] == s.pi.size()) ? 0 : s.revChronBL[s.pi[j]];
            s.revblockEnds[j] = s.revChronBL[s.pi[j] - 1];
        }
    }
}
