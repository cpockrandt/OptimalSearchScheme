#ifndef SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_
#define SEQAN_INDEX_FIND2_INDEX_APPROX_EXTENSION_H_

#include <iostream>
// #include <sdsl/bit_vectors.hpp>
// #include "common.h"
//#include "find2_index_approx_unidirectional.h"


namespace seqan{

template <typename T>
void printv(T a){
    for(int i = 0; i < a.size(); ++i){
        std::cout << static_cast<int> (a.at(i)) << ", ";
    }
    std::cout << "\n";
}

template <size_t nbrBlocks>
void print_search(OptimalSearch<nbrBlocks> const & s){
        std::cout << "Search sscheme: " << (int)s.id << "\n";
        std::cout << "Permutation: " << "\n";
        printv(s.pi);
        std::cout << "Lower bound: " << "\n";
        printv(s.l);
        std::cout << "Upper bound: " << "\n";
        printv(s.u);
        std::cout << "blockLengths: " << "\n";
        printv(s.blocklength);
        std::cout << "chronblockLengths: " << "\n";
        printv(s.chronBL);
        std::cout << "revchronblockLengths: " << "\n";
        printv(s.revChronBL);
        std::cout << "start Pos: " << "\n";
        std::cout << s.startPos << "\n";
        std::cout << "minMax: " << "\n";
        printv(s.min);
        printv(s.max);
        std::cout << "OneDirection" << "\n" << (int)s.startUniDir << "\n";
        std::cout << "\n";
}

template<typename TVector, typename TVSupport,
         typename TSAValue>
inline TVector & getTVector(std::vector<std::pair<TVector, TVSupport> > & bitvectors,
           Pair<uint8_t, Pair<TSAValue, TSAValue>> const & brange)
{
    return bitvectors[brange.i1].first;
}

template<typename TVector, typename TVSupport,
         typename TSAValue>
inline TVSupport & getTVSupport(std::vector<std::pair<TVector, TVSupport> > & bitvectors,
           Pair<uint8_t, Pair<TSAValue, TSAValue>> const & brange)
{
    return bitvectors[brange.i1].second;
}

template<typename TVector, typename TVSupport,
         typename TSAValue>
inline TVector & getTVector(std::vector<std::pair<TVector, TVSupport>* > & bitvectors,
           Pair<uint8_t, Pair<TSAValue, TSAValue>> const & brange)
{
    return bitvectors[brange.i1]->first;
}

template<typename TVector, typename TVSupport,
         typename TSAValue>
inline TVSupport & getTVSupport(std::vector<std::pair<TVector, TVSupport>* > & bitvectors,
           Pair<uint8_t, Pair<TSAValue, TSAValue>> const & brange)
{
    return bitvectors[brange.i1]->second;
}


template <typename TBitvector>
inline void getConsOnes(std::vector<TBitvector> & bitvectors,
                        Pair<uint8_t, Pair<uint32_t, uint32_t>> & inside_bit_interval,
                        uint32_t const intervalsize,
                        std::vector<std::pair<uint32_t, uint32_t>> & consOnesOutput)
{
    auto & b = getTVector(bitvectors, inside_bit_interval);//bitvectors[inside_bit_interval.i1].first;
    uint32_t k = inside_bit_interval.i2.i1;
    uint32_t startOneInterval = inside_bit_interval.i2.i1;
    while(k < inside_bit_interval.i2.i2){
        uint32_t interval = 0;
        //TODO delete second condition it should end with 1
        while(b[k + interval] == 0 && (k + interval) < inside_bit_interval.i2.i2){
            ++interval;
        }
        if(interval >= intervalsize){
            consOnesOutput.push_back(std::make_pair(startOneInterval, k));
            startOneInterval = k + interval;
        }
        k += interval;
        interval = 0;
        ++k;
    }
    consOnesOutput.push_back(std::make_pair(startOneInterval, k));
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TBitvector,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void filter_interval(TContex & ossContext,
                            TDelegate & delegate,
                            TDelegateD & delegateDirect,
                            Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                            TNeedle const & needle,
                            uint32_t needleId,
                            std::vector<TBitvector> & bitvectors,
                            uint32_t const needleLeftPos,
                            uint32_t const needleRightPos,
                            uint8_t const errors,
                            OptimalSearch<nbrBlocks> const & s,
                            uint8_t const blockIndex,
                            Pair<uint8_t, Pair<uint32_t, uint32_t>> & inside_bit_interval,
                            TDir const & ,
                            TDistanceTag const &)
{
    std::vector<std::pair<uint32_t, uint32_t>> consOnes;
    getConsOnes(bitvectors, inside_bit_interval, ossContext.normal.intervalsize, consOnes);
    uint32_t noi = countSequences(*iter.fwdIter.index);

    //TODO shorten this
    for(uint32_t i = 0; i < consOnes.size(); ++i){
        if (std::is_same<TDir, Rev>::value){
            iter.revIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.revIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter.revIter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, false, Rev(), TDistanceTag());
        }
        else
        {
            iter.fwdIter.vDesc.range.i1 = consOnes[i].first + noi;
            iter.fwdIter.vDesc.range.i2 = consOnes[i].second + noi;
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter.fwdIter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, false, Fwd(), TDistanceTag());
        }
    }
}

template <typename TContex,
          typename TDelegateD,
          typename TNeedle,
          size_t nbrBlocks>
inline void genomeSearch(TContex & ossContext,
                         TDelegateD & delegateDirect,
                         TNeedle const & needle,
                         uint32_t needleId,
                         uint8_t errors,
                         OptimalSearch<nbrBlocks> const & s,
                         uint8_t const blockIndex,
                         auto const & genome,
                         Pair<uint16_t, uint32_t> const & sa_info,
                         std::array<uint32_t, nbrBlocks> & blockStarts,
                         std::array<uint32_t, nbrBlocks> & blockEnds)
{
    for(uint32_t j = 0; j < nbrBlocks - blockIndex; ++j){
        // compare bases to needle
        for(uint32_t k = blockStarts[j]; k <  blockEnds[j]; ++k){
            if(needle[k] != genome[sa_info.i1][sa_info.i2 + k]){
                ++errors;
            }
        }
        if(errors < s.l[blockIndex + j] || errors > s.u[blockIndex + j]){
            return;
        }
    }
    delegateDirect(ossContext, sa_info, posAdd(sa_info, length(needle)), needle, needleId, errors);

}

template<typename TBitvector>
inline bool checkSinglePos(std::vector<TBitvector> & bitvectors,
                           Pair<uint8_t, Pair<uint32_t, uint32_t> > & brange,
                           uint32_t offset)
{
    if(bitvectors.empty()){
        return true;
    }
    else
    {
        auto & b = getTVector(bitvectors, brange);
        return (b[brange.i2.i1 + offset] == 1);
//         return (bitvectors[brange.i1].first[brange.i2.i1 + offset] == 1);
    }
}

inline void saPosOnFwd(Pair<uint16_t, uint32_t> & sa_info,
                       uint32_t const genomelength,
                       uint32_t const occLength)
{
    sa_info.i2 = genomelength - sa_info.i2 - occLength;
}


template <typename TContex,
          typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TBitvector,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void directSearch(TContex & ossContext,
                         TDelegateD & delegateDirect,
                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                         TNeedle const & needle,
                         uint32_t needleId,
                         std::vector<TBitvector> & bitvectors,
                         uint32_t const needleLeftPos,
                         uint32_t const needleRightPos,
                         uint8_t const errors,
                         OptimalSearch<nbrBlocks> const & s,
                         uint8_t const blockIndex,
                         Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                         TDir const & ,
                         TDistanceTag const &)
{
//     std::cout << "ITV\n";
    auto const & genome = indexText(*iter.fwdIter.index);

    if (std::is_same<TDistanceTag, EditDistance>::value){
        std::cerr << "Edit Distance not implemented\n";
        exit(0);

//         //TODO if we are only interested in the best hit call return after delegate calls
//         uint16_t needleL = length(needle);
//         uint8_t max_e = s.u[s.u.size() - 1];
//         uint8_t intIns = 0;
//         uint8_t intDel = 0;
//         //calculate net sum of internal Insertions - Deletions
//
//         if(repLength(iter) < needleRightPos - needleLeftPos - 1)
//             intIns = needleRightPos - needleLeftPos - 1 - repLength(iter);
//         else
//             intDel = repLength(iter) - (needleRightPos - needleLeftPos - 1);
//         uint8_t overlap_l = max_e;
//         uint8_t overlap_r = max_e;
//
//         uint16_t ex_infixL = needleL + overlap_l + overlap_r;
//         for(uint32_t r = 0; r < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1; ++r)
//         {
//     //         if(bitvectors[brange.i1].first[brange.i2.i1 + r] == 1){
//             if(checkSinglePos(bitvectors, brange, r)){
//                 Pair<uint16_t, uint32_t> sa_info;
//                 uint32_t chromlength;
//                 if(std::is_same<TDir, Rev>::value){
//                     sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + r];
//                     uint32_t seqOffset = getSeqOffset(sa_info);
//                     chromlength = length(genome[getSeqNo(sa_info)]);
//                     if(!(needleLeftPos + overlap_l <= seqOffset  && chromlength - 1 >= seqOffset - needleLeftPos + needleL - 1 + overlap_r))
//                         continue;
//                     setSeqOffset(sa_info, seqOffset - needleLeftPos);
// //                     sa_info.i2 = sa_info.i2 - needleLeftPos;
//                 }
//                 else
//                 {
//                     sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + r];
//                     uint32_t seqOffset = getSeqOffset(sa_info);
//                     chromlength = length(genome[getSeqNo(sa_info)]);
//                     if(!(chromlength - 1 >= seqOffset + needleRightPos - 1 + overlap_r && seqOffset + needleRightPos - 1 - overlap_l >= length(needle) + 1))
//                         continue;
// //                     sa_info.i2 = chromlength - sa_info.i2 - needleRightPos + 1;
//                     setSeqOffset(sa_info, chromlength - seqOffset - needleRightPos + 1);
//                 }
//
//                 uint32_t seqOffset = getSeqOffset(sa_info);
//                 DnaString const & ex_infix = infix(genome[getSeqNo(sa_info)], seqOffset - overlap_l,seqOffset + needleL + overlap_r);
//                 DnaString const & n_infix = infix(genome[getSeqNo(sa_info)], seqOffset, seqOffset + needleL);
//
//                 alignmentMyersBitvector(ossContext, delegateDirect, needle, needleId, n_infix, ex_infix, chromlength, sa_info, max_e, overlap_l, overlap_r, intDel, false);
//             }
//         }
    }
    else
    {
        std::array<uint32_t, nbrBlocks> blockStarts;
        std::array<uint32_t, nbrBlocks> blockEnds;
        std::copy(std::begin(s.blockStarts) + blockIndex, std::end(s.blockStarts), std::begin(blockStarts));
        std::copy(std::begin(s.blockEnds) + blockIndex, std::end(s.blockEnds), std::begin(blockEnds));

        if(std::is_same<TDir, Rev>::value){
            //modify blockstart in case we are still inside a block
            if(needleRightPos - 1 > blockStarts[0] && needleRightPos - 1 < blockEnds[0])
                blockStarts[0] = needleRightPos - 1;

            for(uint32_t i = 0; i < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1; ++i){
//                 if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                if(checkSinglePos(bitvectors, brange, i)){
                    // mappability information is in reverse index order if we use the forward index
                    Pair<uint16_t, uint32_t> sa_info = iter.fwdIter.index->sa[iter.fwdIter.vDesc.range.i1 + i];
                    uint32_t const chromlength = length(genome[sa_info.i1]);
                    //Info make sure we dont DS search something going over the chromosom edge
                    //check left chromosom boundry && check right chromosom boundry
                    if(!(needleLeftPos <= sa_info.i2 && chromlength - 1 >= sa_info.i2 - needleLeftPos + length(needle) - 1))
                        continue;

                    sa_info.i2 = sa_info.i2 - needleLeftPos;

                    //search remaining blocks
                    genomeSearch(ossContext, delegateDirect, needle, needleId, errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
                }
            }
        }
        else
        {
            //modify blockend in case we are still inside a block
            if(needleLeftPos > blockStarts[0] && needleLeftPos < blockEnds[0])
                blockEnds[0] = needleLeftPos;

            for(uint32_t i = 0; i < iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1; ++i){
//                 if(bitvectors[brange.i1].first[brange.i2.i1 + i] == 1){
                if(checkSinglePos(bitvectors, brange, i)){
                    Pair<uint16_t, uint32_t> sa_info = iter.revIter.index->sa[iter.revIter.vDesc.range.i1 + i];
                    uint32_t const chromlength = length(genome[sa_info.i1]);
                    //check left chromosom boundry && check right chromosom boundry
                    if(!(chromlength - 1 >= sa_info.i2 + needleRightPos - 1 && sa_info.i2 + needleRightPos - 1 >= length(needle) + 1))
                        continue;
                    //calculate correct starting position of the needle  on the forward index
                    sa_info.i2 = chromlength - sa_info.i2 - needleRightPos + 1;

                    //search remaining blocks
                    genomeSearch(ossContext, delegateDirect, needle, needleId , errors, s, blockIndex, genome, sa_info, blockStarts, blockEnds);
                }
            }
        }
    }
}

template <typename TText, typename TIndex, typename TIndexSpec,
          typename TDir>
inline void request_bitvector_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                       uint8_t needed_bitvector,
                                       Pair<uint8_t, Pair<uint32_t, uint32_t>> & brangeOutput,
                                       TDir const & )
{
    Pair<uint32_t, uint32_t> dirrange = (std::is_same<TDir, Rev>::value) ? range(iter.fwdIter) : range(iter.revIter);
    uint32_t nseq = countSequences(*iter.fwdIter.index);
    dirrange.i1 = dirrange.i1 - nseq;
    dirrange.i2 = dirrange.i2 - nseq;

    brangeOutput.i1 = needed_bitvector;
    brangeOutput.i2 = dirrange;
}

//TODO load bitvectors inside a struct to make accessing the correct bitvector easier
template <typename TText, typename TIndex, typename TIndexSpec,
          typename TBitvector,
          size_t nbrBlocks>
inline void get_bitvector_interval_inside(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                          std::vector<TBitvector> & bitvectors,
                                          OptimalSearch<nbrBlocks> const & s,
                                          uint8_t const blockIndex,
                                          Pair<uint8_t, Pair<uint32_t, uint32_t>> & brangeOutput,
                                          bool const goToRight2)
{
    Pair<uint32_t, uint32_t> dirrange = (goToRight2) ? range(iter.revIter) : range(iter.fwdIter);
    uint8_t needed_bitvector;
    uint8_t size = s.pi.size();
    uint8_t bitvsize = bitvectors.size();
    if (goToRight2)
        needed_bitvector = bitvsize - s.max[blockIndex - 1];
    else
        needed_bitvector = s.min[blockIndex - 1] - 1;

    uint32_t nseq = countSequences(*iter.fwdIter.index);
    dirrange.i1 = dirrange.i1 - nseq;
    dirrange.i2 = dirrange.i2 - nseq;

    brangeOutput.i1 = needed_bitvector;
    brangeOutput.i2 = dirrange;
}


template <typename TText, typename TIndex, typename TIndexSpec,
          typename TBitvector,
          size_t nbrBlocks,
          typename TDir>
inline void get_bitvector_interval(Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                       std::vector<TBitvector> & bitvectors,
                       OptimalSearch<nbrBlocks> const & s,
                       uint8_t const blockIndex,
                       Pair<uint8_t, Pair<uint32_t, uint32_t>> & brangeOutput,
                       TDir const & )
{
    uint8_t needed_bitvector;
    if (std::is_same<TDir, Rev>::value)
        needed_bitvector = s.min[blockIndex] - 1;
    else
        needed_bitvector = bitvectors.size() - s.max[blockIndex];// + 1 - 1

    request_bitvector_interval(iter, needed_bitvector, brangeOutput, TDir());
}

template<typename TContex,
         typename TText, typename TIndex, typename TIndexSpec,
         typename TBitvector,
         size_t nbrBlocks>
inline bool testUnidirectionalFilter(TContex & ossContext,
                                     Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                     std::vector<TBitvector> & bitvectors,
                                     Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                                     OptimalSearch<nbrBlocks> const & s,
                                     uint8_t const blockIndex,
                                     bool const goToRight2)
{
    // need bitinterval from inside the pattern to filter according to the mappability form
    //therefore i also need to acces the block before because of that block i got mappability of both sides
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
    get_bitvector_interval_inside(iter, bitvectors, s, blockIndex, bit_interval, goToRight2);
    auto & b2 = getTVector(bitvectors, bit_interval);//bitvectors[bit_interval.i1].first;

    //squash interval
    uint32_t startPos = bit_interval.i2.i1, endPos = bit_interval.i2.i2;
    uint32_t startPos2 = startPos;
    uint32_t endPos2 = endPos;

    while(b2[startPos] == 0 && startPos < endPos)
        ++startPos;

    while(b2[endPos - 1] == 0 && endPos > startPos)
        --endPos;

    if(startPos > endPos){
        std::cerr << "Error bit vector has only zeroes this should have been checked by checkinterval" << "\n";
        exit(0);
    }

    float ivalSize = brange.i2.i2 - brange.i2.i1;
    uint32_t count = 0;

    if(ossContext.normal.testflipdensity){
        // order of bits
        bool last = b2[startPos];
        uint32_t pos = startPos;
        while(pos < endPos){
            if(b2[pos] != last){
                ++count;
                last = !last;
            }
            ++pos;
        }
    }
    // if next condition is true then brange will be modified!!!!
    // it will contain mappability of the bitvector anchored at the other side of already searched needle

    // only interested in changes inside the supinterval (startPos - endPos)
    // allowed flips per intervalSize
    if(!ossContext.normal.testflipdensity || ivalSize * ossContext.normal.invflipdensity - 1 > static_cast<float>(count)){
        brange.i1 = bit_interval.i1;
        brange.i2.i1 = startPos;
        brange.i2.i2 = endPos;
        return true;
    }
    return false;
}

template<typename TContex,
         typename TBitvector,
         size_t nbrBlocks>
inline ReturnCode checkInterval(TContex & ossContext,
                                std::vector<TBitvector> & bitvectors,
                                Pair<uint8_t, Pair<uint32_t, uint32_t>> & brange,
                                OptimalSearch<nbrBlocks> const & s,
                                uint8_t const blockIndex)
{
    auto & b = getTVector(bitvectors, brange);//bitvectors[brange.i1].first;
    auto & rb = getTVSupport(bitvectors, brange);//bitvectors[brange.i1].second;
    rb.set_vector(&b);

    uint32_t ivalOne = rb(brange.i2.i2) - rb(brange.i2.i1);
    uint32_t ivalSize = brange.i2.i2 - brange.i2.i1;

    if(ossContext.normal.nomappability && ivalOne == 0)
        return ReturnCode::NOMAPPABILITY;

//     ivalOne < (s.pi.size() - blockIndex - 1 + ossContext.normal.directsearchblockoffset) * ossContext.normal.directsearch_th
    if(ossContext.normal.directsearch && ossContext.itvCondition(s, blockIndex, ivalOne))
        return ReturnCode::DIRECTSEARCH;

    if(ossContext.normal.compmappable && ivalOne == (brange.i2.i2 - brange.i2.i1)) //TODO maybe allow some zeroes
        return ReturnCode::COMPMAPPABLE;

    //equal or more than half zeroes
    if(ossContext.normal.suspectunidirectional && s.startUniDir <= blockIndex && static_cast<float>(ivalOne) / static_cast<float>(ivalSize) <= ossContext.normal.filter_th)
        return ReturnCode::SUSPECTUNIDIRECTIONAL;

    return ReturnCode::MAPPABLE;
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TBitvector,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline ReturnCode checkCurrentMappability(TContex & ossContext,
                                          TDelegate & delegate,
                                          TDelegateD & delegateDirect,
                                          Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                          TNeedle const & needle,
                                          uint32_t needleId,
                                          std::vector<TBitvector> & bitvectors,
                                          uint32_t const needleLeftPos,
                                          uint32_t const needleRightPos,
                                          uint8_t const errors,
                                          OptimalSearch<nbrBlocks> const & s,
                                          uint8_t const blockIndex,
                                          uint8_t const minErrorsLeftInBlock,
                                          TDir const & ,
                                          TDistanceTag const &)
{
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
    get_bitvector_interval(iter, bitvectors, s, blockIndex, bit_interval, TDir());

    ReturnCode rcode = checkInterval(ossContext, bitvectors, bit_interval, s, blockIndex);

    switch(rcode){
        case ReturnCode::NOMAPPABILITY:
            return ReturnCode::FINISHED;

        case ReturnCode::DIRECTSEARCH:
        {
            directSearch(ossContext, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        case ReturnCode::COMPMAPPABLE:
        {
            std::vector<TBitvector> empty_bitvectors;
            _optimalSearchSchemeChildren(ossContext, delegate, delegateDirect, iter, needle, needleId, empty_bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        default:
            return ReturnCode::MAPPABLE;
    }
}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TBitvector,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline ReturnCode checkMappability(TContex & ossContext,
                                   TDelegate & delegate,
                                   TDelegateD & delegateDirect,
                                   Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                   TNeedle const & needle,
                                   uint32_t needleId,
                                   std::vector<TBitvector> & bitvectors,
                                   uint32_t const current_needleLeftPos,
                                   uint32_t const current_needleRightPos,
                                   uint8_t const errors,
                                   OptimalSearch<nbrBlocks> const & s,
                                   uint8_t const blockIndex,
                                   bool const lastEdit,
                                   TDir const & ,
                                   TDistanceTag const &)
{
    Pair<uint8_t, Pair<uint32_t, uint32_t>> bit_interval;
    get_bitvector_interval(iter, bitvectors, s, blockIndex, bit_interval, TDir());

    ReturnCode rcode = checkInterval(ossContext, bitvectors, bit_interval, s, blockIndex);
    switch(rcode)
    {
        case ReturnCode::NOMAPPABILITY:
            return ReturnCode::FINISHED;

        case ReturnCode::DIRECTSEARCH:
        {
            //search directly in Genome
            directSearch(ossContext, delegateDirect, iter, needle, needleId, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        case ReturnCode::COMPMAPPABLE:
        {
            std::vector<TBitvector> empty_bitvectors;
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, empty_bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, lastEdit, TDir(), TDistanceTag());
            return ReturnCode::FINISHED;
        }

        case ReturnCode::SUSPECTUNIDIRECTIONAL:
        {
            //test unidirectional changes iter range if true
            //TODO modfy functions for TDIR
            bool goToRight2 = std::is_same<TDir, Rev>::value;
            if(testUnidirectionalFilter(ossContext, iter, bitvectors, bit_interval, s, blockIndex, goToRight2)){
                //range on iter was changed in function before
                filter_interval(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, current_needleLeftPos, current_needleRightPos, errors, s, blockIndex, bit_interval, TDir(), TDistanceTag());
                return ReturnCode::FINISHED;
            }
        }
        default:
            return ReturnCode::MAPPABLE;
    }
}

template <typename TContex,
          typename TDelegate,
          typename TDelegateDirect,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TBitvector,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeDeletion(TContex & ossContext,
                                         TDelegate & delegate,
                                         TDelegateDirect & delegateDirect,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         uint32_t needleId,
                                         std::vector<TBitvector> & bitvectors,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         bool const lastEdit,
                                         TDir const & )
{
//     std::cout << "Del\n";
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    if (minErrorsLeftInBlock == 0)
    {

        uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
        bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

        if (goToRight2)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, lastEdit, Rev(), EditDistance());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex2, lastEdit, Fwd(), EditDistance());
    }

    if (maxErrorsLeftInBlock > 0 && goDown(iter, TDir()))
    {
        do
        {
            _optimalSearchSchemeDeletion(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, true, TDir());
        } while (goRight(iter, TDir()));
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TBitvector,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeChildren(TContex & ossContext,
                                         TDelegate & delegate,
                                         TDelegateD & delegateDirect,
                                         Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                         TNeedle const & needle,
                                         uint32_t needleId,
                                         std::vector<TBitvector> & bitvectors,
                                         uint32_t const needleLeftPos,
                                         uint32_t const needleRightPos,
                                         uint8_t const errors,
                                         OptimalSearch<nbrBlocks> const & s,
                                         uint8_t const blockIndex,
                                         uint8_t const minErrorsLeftInBlock,
                                         TDir const & ,
                                         TDistanceTag const &)
{
    bool goToRight = std::is_same<TDir, Rev>::value;
    if (goDown(iter, TDir()))
    {
        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()), needle[goToRight ? needleRightPos - 1 : needleLeftPos - 1]);
            if (!std::is_same<TDistanceTag, EditDistance>::value && minErrorsLeftInBlock > 0 && charsLeft + delta < minErrorsLeftInBlock + 1u)
            {
                continue;
            }
            int32_t needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t needleRightPos2 = needleRightPos + goToRight;
            //finished Block
            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                if (std::is_same<TDistanceTag, EditDistance>::value)
                {
                    //use delta instead of false if no mismatches are allowed
                    _optimalSearchSchemeDeletion(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, false, TDir());
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, false, Rev(), TDistanceTag());
                    else
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex2, false, Fwd(), TDistanceTag());
                }
            }
            else
            {
                //if want to disable mismatches at the start and end (!delta || not_at_end) && use delta instead of false
                _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + delta, s, blockIndex, false, TDir(), TDistanceTag());
            }

            //Deletion
            if (std::is_same<TDistanceTag, EditDistance>::value)
                _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors + 1, s, blockIndex, true, TDir(), TDistanceTag());
        } while (goRight(iter, TDir()));
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TBitvector,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeExact(TContex & ossContext,
                                      TDelegate & delegate,
                                      TDelegateD & delegateDirect,
                                      Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                      TNeedle const & needle,
                                      uint32_t needleId,
                                      std::vector<TBitvector> & bitvectors,
                                      uint32_t const needleLeftPos,
                                      uint32_t const needleRightPos,
                                      uint8_t const errors,
                                      OptimalSearch<nbrBlocks> const & s,
                                      uint8_t const blockIndex,
                                      TDir const &,
                                      TDistanceTag const &)
{
    // not in last block and next Block is larger then current block
    bool goToRight2 = (blockIndex < s.pi.size() - 1) ? s.pi[blockIndex + 1] > s.pi[blockIndex] : s.pi[blockIndex] > s.pi[blockIndex - 1];
    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
    if (std::is_same<TDir, Rev>::value)
    {
        //search take rest of the block and search it reverse
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;

        if (!goDown(iter, infix(needle, infixPosLeft, infixPosRight + 1), TDir()))
            return;

        if (goToRight2)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, infixPosRight + 2, errors, s, blockIndex2, false, Fwd(), TDistanceTag());
    }
    else
    {
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        int32_t infixPosLeft = needleRightPos - s.blocklength[blockIndex] - 1;
        int32_t infixPosRight = needleLeftPos - 1;

        while (infixPosRight >= infixPosLeft)
        {
            if (!goDown(iter, needle[infixPosRight], TDir()))
                return;
            --infixPosRight;
        }

        if (goToRight2)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, infixPosLeft, needleRightPos, errors, s, blockIndex2, false, Fwd(), TDistanceTag());
    }
}

template <typename TContex,
          typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TBitvector>
inline void filteredDelegate(TContex & ossContext,
                             TDelegate & delegate,
                             Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                             TNeedle const & needle,
                             uint32_t needleId,
                             std::vector<TBitvector> & bitvectors,
                             uint8_t const errors)
{
    Pair<uint8_t, Pair<uint32_t, uint32_t>> left_bit_interval;
    request_bitvector_interval(iter, 0, left_bit_interval, Rev());

    uint32_t rangeStart = iter.fwdIter.vDesc.range.i1;
    uint32_t rangeEnd = iter.fwdIter.vDesc.range.i2;
    uint32_t lastStart = 0;
    for(uint32_t i = 0; i < rangeEnd - rangeStart; ++i)
    {
        auto & b = getTVector(bitvectors, left_bit_interval);
        if(b[left_bit_interval.i2.i1 + i] == 0)//bitvectors[left_bit_interval.i1].first[left_bit_interval.i2.i1 + i] == 0
        {
            if(i != lastStart){
                iter.fwdIter.vDesc.range.i1 = rangeStart + lastStart;
                iter.fwdIter.vDesc.range.i2 = rangeStart + i - 1;
                delegate(ossContext, iter, needle, needleId, errors, false);
            }
            lastStart = i + 1;
        }
    }
    if(lastStart < rangeEnd - rangeStart){
        iter.fwdIter.vDesc.range.i1 = rangeStart + lastStart;
        iter.fwdIter.vDesc.range.i2 = rangeStart + rangeEnd - rangeStart;
        delegate(ossContext, iter, needle, needleId, errors, false);
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          typename TBitvector,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 std::vector<TBitvector> & bitvectors,
                                 uint32_t const needleLeftPos,
                                 uint32_t const needleRightPos,
                                 uint8_t const errors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 uint8_t const blockIndex,
                                 bool const lastEdit,
                                 TDir const & ,
                                 TDistanceTag const &)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;
    bool const done = minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1;
    bool const atBlockEnd = (blockIndex > 0) ? needleRightPos - needleLeftPos - 1 == s.blocklength[blockIndex - 1] : false;        //is not true if we finished needle
    bool const checkMappa = !bitvectors.empty();

//     print_search(s);
//     std::cout << "Check if done\n";
//     std::cout << minErrorsLeftInBlock << "\tNPL: " << needleLeftPos << "\tNRP: " << needleRightPos << "\n";
    // Done. (Last step)
    if (done)
    {
//         std::cout << "Done: " << "\n";
        //last input only matters for unidirectional searches (has to be false in this case)
        if(true){
            if(checkMappa){
                filteredDelegate(ossContext, delegate, iter, needle, needleId, bitvectors, errors);
            }
            else
            {
                delegate(ossContext, iter, needle, needleId, errors, false);
            }
        }
        return;
    }


//     if(atBlockEnd && checkMappa){
//         ReturnCode rcode = checkMappability(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, lastEdit, TDir(), TDistanceTag());
//         if(rcode == ReturnCode::FINISHED)
//             return;
//     }

    // Exact search in current block.
    if (maxErrorsLeftInBlock == 0)
    {
        _optimalSearchSchemeExact(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(), TDistanceTag());
    }
    else if(!checkMappa && ossContext.itvConditionComp(iter, needleLeftPos, needleRightPos, errors, s, blockIndex))
    {
        //give emtpy bitvector and bitvector range sine we will not check mappability
        Pair<uint8_t, Pair<uint32_t, uint32_t>> dummy_bit_interval;
         directSearch(ossContext, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, dummy_bit_interval, TDir(), TDistanceTag());
    }

    // Approximate search in current block.
    else
    {

        // Insertion
        if (std::is_same<TDistanceTag, EditDistance>::value)
        {
            bool const goToRight = std::is_same<TDir, Rev>::value;
            int32_t const needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t const needleRightPos2 = needleRightPos + goToRight;

            //if we are at the end of block we need to add possible deletions because _optimalSearchScheme does not check it
            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                uint8_t const minErrorsLeftInBlock2 = (s.l[blockIndex] > (errors + 1)) ? (s.l[blockIndex] - (errors + 1)) : 0;
                if (minErrorsLeftInBlock2 == 0)
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

                    if (goToRight2)
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, true, Rev(), TDistanceTag());
                    else
                        _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex2, true, Fwd(), TDistanceTag());
                }
            }
            else
            {
                _optimalSearchScheme(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex, true, TDir(), TDistanceTag());
            }
        }

//         //checkCurrentMappability
//         uint32_t pblocklength = (blockIndex > 0) ? s.blocklength[blockIndex - 1] : 0;
//         uint32_t step = (needleRightPos - needleLeftPos - 1);
//         if(!atBlockEnd && checkMappa && ossContext.inBlockCheckMappabilityCondition(needleLeftPos, needleRightPos, s, blockIndex))
//         {
//             ReturnCode rcode = checkCurrentMappability(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
//             if(rcode == ReturnCode::FINISHED)
//                 return;
//         }
        _optimalSearchSchemeChildren(ossContext, delegate, delegateDirect, iter, needle, needleId, bitvectors, needleLeftPos, needleRightPos, errors, s, blockIndex, minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TBitvector,
          size_t nbrBlocks,
          typename TNeedle,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 std::vector<TBitvector> & bitvectors,
                                 OptimalSearch<nbrBlocks> const & s,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 TDistanceTag const &)
{
    bool initialDirection = s.pi[1] > s.pi[0];
        if(initialDirection)
            _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Rev(), TDistanceTag());
        else
            _optimalSearchScheme(ossContext, delegate, delegateDirect, it, needle, needleId, bitvectors, s.startPos, s.startPos + 1, 0, s, 0, false, Fwd(), TDistanceTag());
}

template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TBitvector,
          size_t nbrBlocks, size_t N,
          typename TNeedle,
          typename TDistanceTag>
inline void _optimalSearchScheme(TContex & ossContext,
                                 TDelegate & delegate,
                                 TDelegateD & delegateDirect,
                                 Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                 std::vector<TBitvector> & bitvectors,
                                 std::array<OptimalSearch<nbrBlocks>, N> const & ss,
                                 TNeedle const & needle,
                                 uint32_t needleId,
                                 TDistanceTag const &)
{


    for (auto & s : ss){
//         print_search(s);
        _optimalSearchScheme(ossContext, delegate, delegateDirect, it, bitvectors, s, needle, needleId, TDistanceTag());
    }
}

template </*size_t minErrors, size_t maxErrors,*/
          typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndex, typename TIndexSpec,
//           typename TText, typename TIndexSpec,
          typename TBitvector,
          size_t nbrBlocks, size_t N,
          typename TChar, typename TStringSpec,
          typename TDistanceTag>
inline void
find(TContex & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
//      Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     std::vector<TBitvector> & bitvectors,
     std::array<OptimalSearch<nbrBlocks>, N> const & ss,
     String<TChar, TStringSpec> const & needle,
     uint32_t needleId,
     TDistanceTag const & )
{
//     auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;

    auto scheme = ss;
   _optimalSearchSchemeComputeFixedBlocklength(scheme, length(needle));
   _optimalSearchSchemeComputeChronBlocklength(scheme);

//     Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<> > > it(index);

   _optimalSearchScheme(ossContext, delegate, delegateDirect, it, bitvectors, scheme, needle, needleId, TDistanceTag());
}
/*

template <size_t minErrors, size_t maxErrors,
          typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TBitvector,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void
find(TContex & ossContext,
     TDelegate & delegate,
     TDelegateD & delegateDirect,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     std::vector<TBitvector> & bitvectors, // cant be const since TVSupport.set_vector(&TVector)
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & )
{
    auto scheme = OptimalSearchSchemes<minErrors, maxErrors>::VALUE;
    calcConstParameters(scheme);

//     uint32_t len = length(needles[0]);
//     _optimalSearchSchemeComputeFixedBlocklength(scheme, len);
//     _optimalSearchSchemeComputeChronBlocklength(scheme);
    //load Bitvectors needed for scheme (Blocklength and chronblockLengths have to be calculated therefore I need to assume needle length)
//     std::vector<TBitvector * > lbitvectors;
//     linkBitvectors(ossContext, scheme, bitvectors, lbitvectors);

    uint32_t needleId = 0;
    uint32_t lastcount = 0;
    //neede to fix ++needleID to make it parrallel
    while(needleId < length(needles))
    {
        bool skip = false;
//         if(ossContext.bestXMapper){
//             uint32_t readId = getReadId(needleId, ossContext.readCount);
//             if(isMapped(ossContext.ctx, readId)){
//                 if(getMinErrors(ossContext.ctx, readId) + ossContext.strata < minErrors){
//                     skip = true;
// //                     std::cout << "Skip: " << needleId << "\n";
//                 }
//             }
//         }
        if(!skip){
            find<minErrors, maxErrors>(ossContext, delegate, delegateDirect, index, bitvectors, scheme, needles[needleId], needleId, TDistanceTag());
        }
        //TODO fix this to make parrallelization possible
        ++needleId;
    }
}
//



template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TBitvector,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 std::vector<TBitvector> & bitvectors, // cant be const since TVSupport.set_vector(&TVector)
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 2: find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 3: find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        case 4: find<0, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

//no strata (needed for one Scheme Best X)
template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    std::vector<std::pair<TBitvector, TSupport>> empty_bitvectors;
    find(minErrors, maxErrors, ossContext, delegate, delegateDirect, index, needles, empty_bitvectors, TDistanceTag());
}


template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TBitvector,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 int strata,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 std::vector<TBitvector> & bitvectors, // cant be const since TVSupport.set_vector(&TVector
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    if(strata == 99 || ossContext.oneSSBestXMapper) //TODO this is not necessary anymore
        strata = maxErrors;

    switch (maxErrors)
    {
        case 1:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
        case 2:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
        case 3:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 3 :
                {
                    find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;

        }
        case 4:
        {
            switch (strata)
            {
                case 0 :
                {
                    find<0, 0>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<1, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 1 :
                {
                    find<0, 1>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<2, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 2 :
                {
                    find<0, 2>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<3, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 3 :
                {
                    find<0, 3>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    find<4, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                case 4 :
                {
                    find<0, 4>(ossContext, delegate, delegateDirect, index, bitvectors, needles, TDistanceTag());
                    break;
                }
                default: std::cerr << "strata = " << strata << " not possible with current maxError.\n";
                exit(1);
            }
            break;
        }
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                 exit(1);
    }
}

//no bitvectors so substitute them
template <typename TContex,
          typename TDelegate, typename TDelegateD,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
                 const int maxErrors,
                 const int strata,
                 TContex & ossContext,
                 TDelegate & delegate,
                 TDelegateD & delegateDirect,
                 Index<TText, BidirectionalIndex<TIndexSpec> > & index,
                 StringSet<TNeedle, TStringSetSpec> const & needles,
                 TDistanceTag const & )
{
    std::vector<std::pair<TBitvector, TSupport>> empty_bitvectors;
    find(minErrors, maxErrors, strata, ossContext, delegate, delegateDirect, index, empty_bitvectors, needles, TDistanceTag());
}

// for find2_index_approx.h find function
template <typename TDelegate,
          typename TText, typename TIndexSpec,
          typename TNeedle, typename TStringSetSpec,
          typename TDistanceTag>
inline void find(const int minErrors,
     const int maxErrors,
     TDelegate & delegate,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     StringSet<TNeedle, TStringSetSpec> const & needles,
     TDistanceTag const & )
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(delegate, index, needles, TDistanceTag());
                break;
        case 2: find<0, 2>(delegate, index, needles, TDistanceTag());
                break;
        case 3: find<0, 3>(delegate, index, needles, TDistanceTag());
                break;
        case 4: find<0, 4>(delegate, index, needles, TDistanceTag());
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}*/

}




#endif
