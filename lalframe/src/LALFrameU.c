#include "config.h"

#if defined USE_FRAMEC && ! defined USE_FRAMEL
#include "LALFrameC.c"
#elif defined USE_FRAMEL && ! defined USE_FRAMEC
#include "LALFrameL.c"
#else
#error must define one of USE_FRAMEC or USE_FRAMEL (but not both)
#endif
