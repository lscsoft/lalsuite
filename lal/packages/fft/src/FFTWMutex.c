#include <lal/FFTWMutex.h>

#ifdef LAL_PTHREAD_LOCK
pthread_mutex_t lalFFTWMutex = PTHREAD_MUTEX_INITIALIZER;
#endif
