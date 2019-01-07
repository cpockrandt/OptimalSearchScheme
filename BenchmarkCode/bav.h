#include <iostream>
#include <time.h>
#include <seqan/index.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>

//#define DEBUGOUTPUT

using namespace std;
using namespace seqan;

//für Verifikations-Vorhersage (für Index auf String)
template <typename TDelegate, typename TChar, typename TStringSpec, typename TIndexSpec, typename TText>
void verify(TDelegate &delegate,
            Iter< Index<String<TChar, TStringSpec>, TIndexSpec>, VSTree<TopDown<> > > &it,
            TText const &pattern,
            unsigned long offset,
            unsigned long const errorsAllowed, //TODO zu uint8_t ändern
            unsigned const errorsLeft )  //TODO zu uint8_t ändern
{
    auto const & text = indexText(container(it.fwdIter));                                                               //Text laden
    //unsigned long blockSize = length(pattern)/(errorsAllowed +2); //TODO blcoksize ggf. durchschleppen
    //offset = offset * blockSize;

    for(auto occ : getOccurrences(it))                                                                                  //Für jedes Vorkommen...
    {
#ifdef DEBUGOUTPUT
        cout << "verify 1  pattern:" << pattern << " ,occ: " << occ <<  " rep: " << representative(it.fwdIter)<< endl;
#endif
        if(SEQAN_UNLIKELY(occ < offset))                                                                                //...prüfen ob überhaupt genug Text VOR dem Occurrence für den Rest vom Pattern noch da ist
        {
#ifdef DEBUGOUTPUT
            cout << "verify X  occ too small " << occ << " < offset: " << offset << endl << endl;
#endif
            continue;
        }

        auto const patternOcc = occ - offset;

        if(SEQAN_UNLIKELY(length(text) < (patternOcc + length(pattern))))                                               //...prüfen ob überhaupt genug Text NACH der Occurence für den Rest vom Pattern noch da ist
        {
#ifdef DEBUGOUTPUT
            cout << "verify X  text too short " << length(text) << " < patternEnd " << (patternOcc + length(pattern)) << endl << endl;
#endif
            continue;
        }

        auto const & snippet = infix(text, patternOcc, min(patternOcc + length(pattern), length(text)));                //...zu vergleichende Stelle aus Text laden
        unsigned int errors = 0;

        for(int i = 0; errors <= errorsAllowed && i < length(snippet); ++i)
            errors += snippet[i] != pattern[i];                                                                         //Buchstabenweise vergleichen
        errors += length(pattern) - length(snippet);                                                                    //Hinzuaddieren von Fehlern wenn snippet zu kurz
#ifdef DEBUGOUTPUT
        cout << "verify 2: " << pattern << ", snippet: " << snippet << ", occ: " << occ << endl;
        cout << "verify 3 errors: " << errors << ", occ: " << occ << ", offset: " << offset << ", patternOcc: " << patternOcc << endl << endl;
#endif
        if(errors <= errorsAllowed)                                                                                     //Wenn Fehler nicht überschritten...
        {
            delegate(patternOcc);                                                                                       //...delegate
        }
    }
}

//für Verifikations-Vorhersage (für Index auf StringSets)
template <typename TDelegate, typename TChar, typename TStringSpec, typename TIndexOwner, typename TIndexSpec, typename TText>
void verify(TDelegate &delegate,
            Iter< Index<StringSet<String<TChar, TStringSpec>, TIndexOwner>, TIndexSpec>, VSTree<TopDown<> > > &it,
            TText const &pattern,
            unsigned long offset,
            unsigned long const errorsAllowed, //TODO zu uint8_t ändern
            unsigned const errorsLeft )  //TODO zu uint8_t ändern
{
    auto const & text = indexText(container(it.fwdIter));                                                               //Text laden
    //unsigned long blockSize = length(pattern)/(errorsAllowed +2); //TODO blcoksize ggf. durchschleppen
    //offset = offset * blockSize;

    for(auto occPair : getOccurrences(it))                                                                              //Für jedes Vorkommen...
    {
        unsigned long occ = getSeqOffset(occPair);
        auto const & occString = text[getSeqNo(occPair)];
#ifdef DEBUGOUTPUT
        cout << "verify 1  pattern:" << pattern << " ,occ: " << occ <<  " rep: " << representative(it.fwdIter)<< endl;
#endif
        if(SEQAN_UNLIKELY(occ < offset))                                                                                //...prüfen ob überhaupt genug Text VOR dem Occurrence für den Rest vom Pattern noch da ist
        {
#ifdef DEBUGOUTPUT
            cout << "verify X  occ too small " << occ << " < offset: " << offset << endl << endl;
#endif
            continue;
        }

        auto const patternOcc = occ - offset;

        if(SEQAN_UNLIKELY(length(occString) < (patternOcc + length(pattern))))                                          //...prüfen ob überhaupt genug Text NACH der Occurence für den Rest vom Pattern noch da ist
        {
#ifdef DEBUGOUTPUT
            cout << "verify X  text too short " << length(text) << " < patternEnd " << (patternOcc + length(pattern)) << endl << endl;
#endif
            continue;
        }

        auto const & snippet = infix(occString, patternOcc, min(patternOcc + length(pattern), length(occString)));      //...zu vergleichende Stelle aus Text laden
        unsigned int errors = 0;

        for(int i = 0; errors <= errorsAllowed && i < length(snippet); ++i)
            errors += snippet[i] != pattern[i];                                                                         //Buchstabenweise vergleichen
        errors += length(pattern) - length(snippet);                                                                    //Hinzuaddieren von Fehlern wenn snippet zu kurz
#ifdef DEBUGOUTPUT
        cout << "verify 2: " << pattern << ", snippet: " << snippet << ", occ: " << occ << endl;
        cout << "verify 3 errors: " << errors << ", occ: " << occ << ", offset: " << offset << ", patternOcc: " << patternOcc << endl << endl;
#endif
        if(errors <= errorsAllowed)                                                                                     //Wenn Fehler nicht überschritten...
        {
            delegate(patternOcc);                                                                                       //...delegate
        }
    }
}

//links Sucher
template <typename TDelegate, typename TPredictify, typename TIter, typename TText>
inline void trivialSearchLeft(TDelegate & delegate,
                              TIter it,
                              TText const &pattern,
                              signed const pos,
                              TPredictify &predictify,
                              unsigned const long errorsAllowed,
                              unsigned const errorsLeft,
                              bool const indels)
{

    if(predictify(it, pattern, pos, pos, errorsAllowed, errorsLeft, indels) && pos > 0) //TODO threshold entscheidung richtig einbinden //TODO pos = 0 oder nach exit block if(el == 0...
    {
        verify(delegate, it, pattern, pos+1, errorsAllowed, errorsLeft); //TODO pos+1 (früher pos) prüfen.
        return;
    }

    if (errorsLeft == 0 || pos < 0)
    {
        for(signed long i = pos; i >= 0 ;--i){                                                                          // Abkürzung wenn keine Fehler mehr erlaubt sind...
            if(!goDown(it, pattern[i],Fwd())) return;                                                                   //... und NICHT nach links zu ende gesucht werden kann abbrechen.
        }

        for(auto occ : getOccurrences(it))                                                                               //... und nach links zu ende gesucht werden kann HIT.
        {
            delegate(occ);
        }

#ifdef DEBUGOUTPUT
        cout << "TSL: " << representative(it.fwdIter) << endl;
#endif
        if(indels && errorsLeft > 0 && goDown(it, Fwd())){                                                              // Wenn noch Fehler erlaubt pattern aber schon gefunden...
            do{
                // Deletion (left of pattern)
                trivialSearchLeft(delegate, it, pattern, pos, predictify, errorsAllowed, errorsLeft - 1, indels);                                  //... restliche Fehler nach links versuchen zu erweitern.
            } while (goRight(it, Fwd()));
        }
        return;
    }

    if (indels)
    {
        // Insertion
        trivialSearchLeft(delegate, it, pattern, pos - 1, predictify, errorsAllowed, errorsLeft - 1, indels);
    }

    if (goDown(it, Fwd()))
    {
        do
        {
            unsigned delta = !ordEqual(parentEdgeLabel(it, Fwd()), pattern[pos]);
            trivialSearchLeft(delegate, it, pattern, pos -1, predictify, errorsAllowed, errorsLeft - delta, indels);

            //Deletion
            if(indels){
                trivialSearchLeft(delegate, it, pattern, pos, predictify, errorsAllowed, errorsLeft - 1, indels);
            }
        } while (goRight(it, Fwd()));
    }
}

//rechst Sucher
template <typename TDelegate, typename TPredictify, typename TIter, typename TText>
inline void trivialSearchRight(TDelegate & delegate,
                               TIter it,
                               TText const &pattern,
                               unsigned const pos,
                               TPredictify &predictify,
                               unsigned const long errorsAllowed,
                               unsigned const errorsLeft,
                               unsigned long start,
                               bool const indels)
{

    if(predictify(it, pattern, pos, pos, errorsAllowed, errorsLeft, indels))//TODO threshold entscheidung richtig einbinden
    {
        verify(delegate, it, pattern, start, errorsAllowed, errorsLeft);
        return;
    }

    if(pos >= length(pattern) || errorsLeft == 0){
        if(goDown(it, suffix(pattern, pos), Rev())){
#ifdef DEBUGOUTPUT
            cout << "TSR: " << representative(it.fwdIter) << endl;
#endif
            trivialSearchLeft(delegate, it, pattern, start-1, predictify, errorsAllowed, errorsLeft, indels);
        }

        if(indels && errorsLeft > 0 && goDown(it, Rev())){                                                              // Wenn noch Fehler erlaubt pattern aber schon gefunden...
            do{
#ifdef DEBUGOUTPUT
                cout << "TSR right extend: " << representative(it.fwdIter) << endl;
#endif
                // Deletion (right of pattern)
                trivialSearchRight(delegate, it, pattern, pos, predictify, errorsAllowed, errorsLeft - 1, start, indels);//... restliche Fehler nach rechts versuchen zu erweitern.
            } while (goRight(it, Rev()));
        }

        return;
    }

    if (indels)
    {
        // Insertion
        trivialSearchRight(delegate, it, pattern, pos + 1, predictify, errorsAllowed, errorsLeft - 1, start, indels);
    }

    if (goDown(it, Rev()))
    {
        do
        {
            unsigned delta = !ordEqual(parentEdgeLabel(it, Rev()), pattern[pos]);
            trivialSearchRight(delegate, it, pattern, pos + 1, predictify, errorsAllowed, errorsLeft - delta, start, indels);

            //Deletion
            if(indels){
                trivialSearchRight(delegate, it, pattern, pos, predictify, errorsAllowed, errorsLeft - 1, start, indels);
            }
        } while (goRight(it, Rev()));
    }
}

void printTabs(unsigned n)
{
    for (unsigned j = 0; j < n; ++j)
        cout << "-";
    cout << " ";
}

//01*1 Sucher
template <typename TDelegate, typename TPredictify, typename TIterator, typename TString>
void searchSeed(TDelegate & delegate,
                TString & pattern,
                TIterator it,
                unsigned long pos,
                unsigned long blockSize,
                unsigned long balance,
                TPredictify &predictify,
                unsigned const long errorsAllowed,
                unsigned long errorsLeft,
                bool isZeroBlock,
                unsigned long start,
                bool const indels,
                unsigned tabs = 2)
{
    unsigned long blockNumber = (pos - balance) / blockSize;                                                            //Blocknummer bestimmen für pos in kleineren Blöcken ...
    if(pos < (balance * (blockSize + 1))) blockNumber = pos / (blockSize + 1);                                          // ... korrektur für pos in Blanacierten blöcken.

    unsigned long nextBlock = ((blockNumber + 1) * blockSize);                                                          //Position des ersten Buchstaben im nächsten Block
    if(balance) nextBlock += min(blockNumber + 1, balance);                                                             //Wenn "balance" nicht Null => korrektur

    //nextBlock = min(nextBlock, length(pattern));                                                                      //TODO dürfte nicht passieren //Wenn das pattern nicht gleichmaessig aufgeteilt werden kann wird der letzte block gekürtzt
    if(nextBlock > length(pattern)) exit(66);                                                                           //TODO not needed

    if (pos >= length(pattern))
        return;

    if(errorsLeft == 0){                                                                                                // Wenn keine Fehler mehr erlaubt
        if(goDown(it, suffix(pattern, pos), Rev()))                                                                     // ...und nach rechts alles gefunden...
        {
            trivialSearchLeft(delegate, it, pattern, start-1, predictify, errorsAllowed, 0, indels);                     // ...suche nach links ohne fehler zu ende.
        }
        return;                                                                                                         // sonst abbrechen.
    }

    if (nextBlock >= length(pattern)){                                                                                  //Wenn im letzten Block
        if (goDown(it, infix(pattern, pos, length(pattern)), Rev())) {                                                  //Fehlerfrei zu Ende des Patterns suchen...
            trivialSearchLeft(delegate, it, pattern, start-1, predictify, errorsAllowed, errorsLeft, indels);            //... wenn das geht weiter mit linkssuche.
        }
        return;
    }

    if(!isZeroBlock){                                                                                                   //Wenn schon Fehler im Block gefunden worden
        if(nextBlock - pos > 0)                                                                                         //wenn noch nicht am Ende vom Block //Damit keine leeren Kanten gelaufen werden
        {
            if (!goDown(it, infix(pattern, pos, nextBlock), Rev())) {                                                   //Fehlerfrei zu Ende suchen
                return;                                                                                                 //Abbrechen wenn das nicht möglich ist.
            }
            pos = nextBlock;                                                                                            //pos auf nächsten Blockanfang setzten
        }
        searchSeed(delegate, pattern, it, pos, blockSize, balance, predictify, errorsAllowed, errorsLeft,/*isZeroBlock=*/ true, start, indels, tabs + 1);   //starte SearchSeed für nächsten als fehlerfrei markierten Block
        return;
    }

    //insertions
    if(indels){
        if(nextBlock - pos == 1){                                                                                        //... und am Ende des Blocks
            searchSeed(delegate, pattern, it, pos + 1, blockSize, balance, predictify, errorsAllowed, errorsLeft - 1, /*isZeroBlock=*/ true, start, indels, tabs + 1);//starte SearchSeed für nächsten als fehlerfrei markierten Block
        }
        else{                                                                                                           //... aber noch mittem im Block
            searchSeed(delegate, pattern, it, pos + 1, blockSize, balance, predictify, errorsAllowed, errorsLeft - 1, /*isZeroBlock=*/ false, start, indels, tabs + 1);//starte SearchSeed für nächsten Buchstaben im mit fehler markierten Block
        }
    }

    if (goDown(it, Rev())){                                                                                             //Wenn noch weitere Kanten
        do {
#ifdef DEBUGOUTPUT
            printTabs(tabs);
            //cout << representative(it.fwdIter) << " . " << representative(it.revIter) << endl;
            TString copyRepres = representative(it.revIter);
            reverse(copyRepres);
            cout << copyRepres << endl;
#endif

            bool miss = !ordEqual(parentEdgeLabel(it, Rev()), pattern[pos]);                                            //Prüft ob Kante mit Pattern übereinstimmt

            if (miss) {                                                                                                 //Wenn Miss...
                if(nextBlock - pos == 1){                                                                               //... und am Ende des Blocks
                    searchSeed(delegate, pattern, it, pos + 1, blockSize, balance, predictify, errorsAllowed, errorsLeft - 1, /*isZeroBlock=*/ true, start, indels, tabs + 1);//starte SearchSeed für nächsten als fehlerfrei markierten Block
                }
                else{                                                                                                   //... aber noch mittem im Block
                    searchSeed(delegate, pattern, it, pos + 1, blockSize, balance, predictify, errorsAllowed, errorsLeft - 1, /*isZeroBlock=*/ false, start, indels, tabs + 1);//starte SearchSeed für nächsten Buchstaben im mit fehler markierten Block
                }
            }

            else {                                                                                                      //Wenn Hit..
                //no errors in the Block
                if ((nextBlock - pos) == 1){                                                                            //... und am Ende des Blocks (zweite 0 von 01*0)
#ifdef DEBUGOUTPUT
                    cout << "SS: " << representative(it.fwdIter) << endl; // it.revIter
#endif
                    trivialSearchRight(delegate, it , pattern, pos + 1, predictify, errorsAllowed, errorsLeft, start, indels);                             //starte trivialSearchRight (erweiterungs Suche nach rechts)
                }
                else{                                                                                                   //... aber noch mittem im Block
                    searchSeed(delegate, pattern, it, pos + 1, blockSize, balance, predictify, errorsAllowed, errorsLeft, /*isZeroBlock=*/ true, start, indels, tabs + 1);//starte SearchSeed für nächsten Buchstaben im fehlerfrei markierten Block
                }
            }

            //deletion
            if(indels){
                if(nextBlock - pos == 1){                                                                               //... und am Ende des Blocks
                    searchSeed(delegate, pattern, it, pos, blockSize, balance, predictify, errorsAllowed, errorsLeft - 1, /*isZeroBlock=*/ true, start, indels, tabs + 1);//starte SearchSeed für nächsten als fehlerfrei markierten Block
                }
                else{                                                                                                   //... aber noch mittem im Block
                    searchSeed(delegate, pattern, it, pos, blockSize, balance, predictify, errorsAllowed, errorsLeft - 1, /*isZeroBlock=*/ false, start, indels, tabs + 1);//starte SearchSeed für nächsten Buchstaben im mit fehler markierten Block
                }
            }
        } while (goRight(it, Rev()));                                                                                   //Wiederholen für jede freie Kante
    }
}

//Sucher Eingang
template <typename TDelegate, typename TPredictify, typename TIterator, typename TString>
void search(TDelegate &delegate,          //übergebene Funktion
            TIterator it,                 //index in dem gesucht wird
            TPredictify &predictify,      //Grenzwertfunktion unterwelchem die Verifikations-Vorhersage genutzt wird
            unsigned long errorsAllowed,  //maximal erlaubte Fehler
            TString &pattern,             //pattern nach dem gesucht wird
            bool const indels)            //indels an/aus
{
    unsigned long pLength = length(pattern);
    unsigned long blockCount = errorsAllowed + 2;
    unsigned long blockSize =  (pLength / blockCount);
    unsigned long balance = (pLength % blockCount);                                                                     //Rest von Verteilung pLength auf blockCount
    bool isZeroBlock = true;

#ifdef DEBUGOUTPUT
    cout << ":---------------BAV START----------------:" << endl;
    cout << "pLength:    " << pLength << endl;
    cout << "blockCount: " << blockCount << endl;
    cout << "blockSize:  " << blockSize << endl;
    cout << "balance:    " << balance << endl;
#endif

    if(pLength < blockCount){                                                                                           //Wenn mehr Blöcke als Zeichen...
        trivialSearchLeft(delegate, it, pattern, pLength-1, predictify, errorsAllowed, errorsAllowed, indels);          //... vom Ende nach links Trivial suchen.
        return;
    }

    for(unsigned long block = 0; block < (errorsAllowed+1); ++block){     //nur bis k+1 da no min 0 gefunden werden muss

        goRoot(it);
        unsigned long pos = block*blockSize;
        unsigned long nextBlock = pos + blockSize;                                                                      //Position des ersten Buchstaben im nächsten Block
        if(balance) {
            pos += min(block, balance);
            nextBlock += min(block + 1, balance);
        }

        //cout << "..." << infix(pattern, pos, pos + blockSize) << endl;
        if(goDown(it, infix(pattern, pos, nextBlock), Rev()))
        {
#ifdef DEBUGOUTPUT
            cout << "- " << representative(it.fwdIter) << endl;
#endif
            searchSeed(delegate, pattern, it, nextBlock, blockSize, balance, predictify, errorsAllowed, errorsAllowed, isZeroBlock, pos, indels);

        }
    }
}