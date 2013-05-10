/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/


#ifndef XLALGSL_H
#define XLALGSL_H

#include <gsl/gsl_errno.h>
#include <lal/LALConfig.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
extern pthread_mutex_t lalGSLPthreadMutex;
#define XLALGSL_PTHREAD_MUTEX_LOCK pthread_mutex_lock( &lalGSLPthreadMutex )
#define XLALGSL_PTHREAD_MUTEX_UNLOCK pthread_mutex_unlock( &lalGSLPthreadMutex )
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

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* XLALGSL_H */
