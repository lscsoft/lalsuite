#ifndef _FFTWMUTEX_H
#define _FFTWMUTEX_H

#include <lal/LALConfig.h>
#include <lal/LALRCSID.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( FFTWMUTEXH, "$Id$" );

#ifdef LAL_PTHREAD_LOCK
# include <pthread.h>
  extern pthread_mutex_t lalFFTWMutex;
# define LAL_FFTW_PTHREAD_MUTEX_LOCK pthread_mutex_lock( &lalFFTWMutex )
# define LAL_FFTW_PTHREAD_MUTEX_UNLOCK pthread_mutex_unlock( &lalFFTWMutex )
#else
# define LAL_FFTW_PTHREAD_MUTEX_LOCK
# define LAL_FFTW_PTHREAD_MUTEX_UNLOCK
#endif

#ifdef  __cplusplus
}
#endif

#endif /* _FFTWMUTEX_H */
