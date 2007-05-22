#include <lal/FFTWMutex.h>

#include <lal/LALRCSID.h>
NRCSID (FFTWMUTEXC,"$Id$");

#ifdef LAL_PTHREAD_LOCK
pthread_mutex_t lalFFTWMutex = PTHREAD_MUTEX_INITIALIZER;
#endif
