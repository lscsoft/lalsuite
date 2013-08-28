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

/**
 * \addtogroup LALGSL_h
 * \author Creighton, J. D. E.
 *
 * \brief Provides macros for integrating the GSL error handler with the LAL status structure.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/LALGSL.h>
 * \endcode
 *
 * This header provides macros and functions for tracking and
 * reporting the runtime status of a GSL calls.  The intent is
 * simultaneously to standardize the error reporting, and to make the
 * reporting as transparent as possible to people coding individual
 * routines.
 *
 * <em>Please always use these macros when making a GSL call
 * within LAL.  This will ensure that the LAL functions always have the
 * same behaviour and will also ensure that the LAL functions are reenterant
 * and threadsafe (when LAL is configured appropriately).</em>
 *
 * \heading{GSL function calls}
 *
 * The behaviour of GSL functions depends on the error handler that has been
 * assigned.  In order that LAL functions always have the same behaviour, it
 * is necessary to use a LAL-specific GSL error handler.  This error handler
 * populates a LAL status structure with the GSL error message and code so that
 * GSL functions behave much the same way as LAL functions.  After the GSL
 * functions are called, the error handler needs to be restored to the original
 * handler so that the program calling the LAL routine has the same error handler
 * after the LAL function was called as it did before the LAL function was called.
 *
 * This module provides a simple set of macros and the default LAL GSL error
 * handler.  The macros provide a standard way to assign the LAL GSL error
 * handler before a GSL function call and to restore the original handler after
 * the call.
 *
 * Note that changing the GSL error handler is \e not a thread-safe
 * action.  Therefore it is necessary to block other threads from performing
 * GSL function calls while one thread has changed the handler.  These macros
 * ensure that such blocking is done for GSL function calls
 * <em>within other LAL routines</em> if LAL is configured with the
 * <tt>--enable-pthread-lock</tt> flag.  See below for instructions on how
 * to make other GSL function calls outside LAL thread-safe when used with LAL.
 *
 * \code
 * ATTATCHSTATUSPTR( status );
 * CALLGSL( gsl_function( x ), status );
 * CHECKSTATUSPTR( status );
 * DETATCHSTATUSPTR( status );
 * \endcode
 * Note that the LAL function must attach (and detach) a status pointer as if
 * a LAL routine were called.
 * Note also that you need to use the \c CHECKSTATUSPTR macro to check
 * the status of the call.  The equivalent to the \c TRY macro for GSL
 * functions is the \c TRYGSL macro, which is used as follows:
 * \code
 * ATTATCHSTATUSPTR( status );
 * TRYGSL( gsl_function( x ), status );
 * DETATCHSTATUSPTR( status );
 * \endcode
 *
 * If you are using GSL functions both in LAL and in the calling program, and
 * you are worried about thread-safety, the GSL function calls outside of LAL
 * need to be blocked so that they do not access the GSL error handler while
 * it has been changed to the LAL GSL error handler in a LAL function.  To do
 * this, you need to do the following:
 * \code
 * #include<lal/LALGSL.h>
 * ...
 * LALGSL_PTHREAD_MUTEX_LOCK;
 * gsl_function( x );
 * LALGSL_PTHREAD_MUTEX_UNLOCK;
 * \endcode
 * This ensures that \c gsl_function is not called while a LAL routine
 * is calling a GSL function in a different thread.  You can do this even if
 * you don't always run your code with multiple threads.  If you configure LAL
 * without the <tt>--enable-pthread-lock</tt> flag, the macros
 * \c LALGSL_PTHREAD_MUTEX_LOCK and \c LALGSL_PTHREAD_MUTEX_UNLOCK
 * do nothing.
 *
 */

#ifndef _LALGSL_H
#define _LALGSL_H

#include <lal/LALConfig.h>
#ifdef NDEBUG
#ifndef LAL_NDEBUG
#define LAL_NDEBUG
#endif
#endif

#include <stdlib.h>
#include <string.h>

#include <lal/LALMalloc.h>
#include <lal/LALDatatypes.h>
#include <lal/LALError.h>

#include <lal/XLALGSL.h>

#include <gsl/gsl_errno.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

#ifndef SWIG    /* exclude from SWIG interface */
extern LALStatus *lalGSLGlobalStatusPtr;
#endif /* SWIG */
void
LALGSLErrorHandler(const char *reason,
                   const char *file, int line, int errnum);

#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
extern pthread_mutex_t lalGSLPthreadMutex;
#define LALGSL_PTHREAD_MUTEX_LOCK pthread_mutex_lock( &lalGSLPthreadMutex )
#define LALGSL_PTHREAD_MUTEX_UNLOCK pthread_mutex_unlock( &lalGSLPthreadMutex )
#else
#define LALGSL_PTHREAD_MUTEX_LOCK  ((void)(0))
#define LALGSL_PTHREAD_MUTEX_UNLOCK  ((void)(0))
#endif

/*
 *
 * FIXME: Must disable pthread safety commands.... Problems arise if the
 * statement calls a LAL function that tries to call a GSL function using
 * CALLGSL ... this blocks!  TODO: Instead, use thread-specific globals.
 *
 */

#define CALLGSL( statement, statusptr )                                       \
  if ( (statusptr) )                                                          \
  {                                                                           \
    LALStatus *saveLALGSLGlobalStatusPtr_;                                    \
    gsl_error_handler_t *saveGSLErrorHandler_;                                \
    if ( !( (statusptr)->statusPtr ) )                                        \
      { ABORT( (statusptr), -8, "CALLGSL: null status pointer pointer" ); }   \
    /* LALGSL_PTHREAD_MUTEX_LOCK; */                                          \
    saveGSLErrorHandler_ = gsl_set_error_handler( LALGSLErrorHandler );       \
    saveLALGSLGlobalStatusPtr_ = lalGSLGlobalStatusPtr;                       \
    lalGSLGlobalStatusPtr = (statusptr)->statusPtr;                           \
    statement;                                                                \
    lalGSLGlobalStatusPtr = saveLALGSLGlobalStatusPtr_;                       \
    gsl_set_error_handler( saveGSLErrorHandler_ );                            \
    /* LALGSL_PTHREAD_MUTEX_UNLOCK; */                                        \
  }                                                                           \
  else                                                                        \
    lalAbortHook( "Abort: CALLGSL, file %s, line %d\n"                        \
                  "       Null status pointer passed to CALLGSL\n",           \
                  __FILE__, __LINE__ )


#define TRYGSL( statement, statusptr )                                        \
  if ( (statusptr) )                                                          \
  {                                                                           \
    LALStatus *saveLALGSLGlobalStatusPtr_;                                    \
    gsl_error_handler_t *saveGSLErrorHandler_;                                \
    if ( !( (statusptr)->statusPtr ) )                                        \
      { ABORT( (statusptr), -8, "CALLGSL: null status pointer pointer" ); }   \
    /* LALGSL_PTHREAD_MUTEX_LOCK;  */                                         \
    saveGSLErrorHandler_ = gsl_set_error_handler( LALGSLErrorHandler );       \
    saveLALGSLGlobalStatusPtr_ = lalGSLGlobalStatusPtr;                       \
    lalGSLGlobalStatusPtr = (statusptr)->statusPtr;                           \
    statement;                                                                \
    lalGSLGlobalStatusPtr = saveLALGSLGlobalStatusPtr_;                       \
    gsl_set_error_handler( saveGSLErrorHandler_ );                            \
    /* LALGSL_PTHREAD_MUTEX_UNLOCK; */                                        \
    if ( (statusptr)->statusPtr->statusCode )                                 \
    {                                                                         \
      SETSTATUS( statusptr, -1, "Recursive error" );                          \
      (void) LALError( statusptr, "Statement \"" #statement "\" failed:" );   \
      (void) LALTrace( statusptr, 1 );                                        \
      return;                                                                 \
    }                                                                         \
  }                                                                           \
  else                                                                        \
    lalAbortHook( "Abort: CALLGSL, file %s, line %d\n"                        \
                  "       Null status pointer passed to CALLGSL\n",           \
                  __FILE__, __LINE__ )

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LALGSL_H */
