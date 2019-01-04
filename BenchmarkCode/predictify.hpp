#include <seqan/index.h>

using namespace seqan;

auto predictify	    = [] (auto const &it,
                          DnaString const &pattern,
                          signed const pos_left,
                          signed const pos_right,
                          unsigned long const errorsAllowed,
                          unsigned const errorsLeft,
                          bool const indels)
    {
	return true; 
        
        /*
        Alternativen:
        false; 						// Predictify aus
	true; 						// Immer truncaten	
	(length(pattern)-1 == pos_right); 		// Immer statt links erweiterung truncaten
        (countOccurrences(it) < 5); 			// Anzahl der aktuellen Subreads
        (errorsLeft < 5); 				// Anzahl der verbleibenden Fehler
        (length(pattern) - (pos_right - pos_left) < 5); // Anzahl verbleibender nicht geprüfter Buchstaben
        */
    };

auto predictify_true = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return true; 
    };

auto predictify_false = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return false; 
    };

auto predictify_left = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return (length(pattern)-1 == pos_right); 
    };

auto predictify_right = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return !(length(pattern)-1 == pos_right); 
    };

auto predictify_occ320 = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return (countOccurrences(it) < 320); 
    };

auto predictify_occ160 = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return (countOccurrences(it) < 160); 
    };

auto predictify_occ80 = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return (countOccurrences(it) < 80); 
    };

auto predictify_occ40 = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return (countOccurrences(it) < 40); 
    };

auto predictify_occ20 = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return (countOccurrences(it) < 20); 
    };

auto predictify_occ5 = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return (countOccurrences(it) < 5); 
    };
auto predictify_occ4 = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return (countOccurrences(it) < 4); 
    };
auto predictify_occ3 = [] (auto const &it, DnaString const &pattern, signed const pos_left, signed const pos_right,
                          unsigned long const errorsAllowed, unsigned const errorsLeft, bool const indels)
    {
	return (countOccurrences(it) < 3); 
    };

