
#ifndef XLALGSL_H
#define XLALGSL_H

#include <gsl/gsl_errno.h>
#include <lal/LALConfig.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( XLALGSLH, "$Id$" );

#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
extern pthread_mutex_t lalGSLPthreadMutex;
#define XLALGSL_PTHREAD_MUTEX_LOCK pthread_mutex_lock( &lalGSLPThreadMutex )
#define XLALGSL_PTHREAD_MUTEX_UNLOCK pthread_mutex_unlock( &lalGSLPThreadMutex )
#else
#define XLALGSL_PTHREAD_MUTEX_LOCK  ((void)(0))
#define XLALGSL_PTHREAD_MUTEX_UNLOCK  ((void)(0))
#endif

#define XLAL_CALLGSL( statement ) \
        do { \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          XLALGSL_PTHREAD_MUTEX_LOCK; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off(); \
          statement; \
          gsl_set_error_handler( saveGSLErrorHandler_ ); \
          XLALGSL_PTHREAD_MUTEX_UNLOCK; \
        } while (0)


#endif /* XLALGSL_H */
