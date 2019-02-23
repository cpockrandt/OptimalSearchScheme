#include <type_traits>

#include <seqan/index.h>

using namespace seqan;

// Index type
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
typedef BidirectionalIndex<FMIndex<void, TMyFastConfig> > TIndexConfig;

