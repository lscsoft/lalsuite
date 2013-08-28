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

#include <lal/LALStdlib.h>
#include <lal/LALGSL.h>
#include <gsl/gsl_errno.h>

LALStatus *lalGSLGlobalStatusPtr = NULL;
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
pthread_mutex_t lalGSLPthreadMutex = PTHREAD_MUTEX_INITIALIZER;
#endif


/**
 * \ingroup LALGSL_h
 * \brief LAL GSL error handler.
 *
 * \heading{Synopsis}
 * \code
 * extern LALStatus *lalGSLGlobalStatusPtr;
 * #include <lal/LALConfig.h>
 * #ifdef LAL_PTHREAD_LOCK
 * #include <pthread.h>
 * extern pthread_mutex_t lalGSLPthreadMutex;
 * #endif
 * \endcode
 *
 * \heading{Description}
 *
 * The function \c LALGSLErrorHandler() is the standard GSL error handler
 * for GSL functions called within LAL.  Its function is to take the GSL
 * error code and information and translate them into equivalent aspects of
 * the LAL status structure.  The status structure that is currently identified
 * by \c lalGSLGlobalStatusPtr is populated.  This global variable is to
 * be set to the current status pointer prior to invocation of the GSL function.
 * In addition, the GSL error handler must be set to the LAL GSL error handler
 * prior to the invocation of the GSL function.  Both of these tasks can be
 * done with the macros provided in the header \ref LALGSL_h.
 * However, doing so is not thread-safe.  Thus the macros use the mutex
 * \c lalGSLPthreadMutex to block other threads from making GSL calls
 * <em>from within LAL functions</em> while one thread has set the GSL error
 * handler and global status pointer.  This mutex must also be used to block
 * any non-LAL GSL function calls from other threads or else they may be called
 * with the LAL GSL error handler in effect.
 */
void
LALGSLErrorHandler(const char *reason,
                   const char *file, int line, int my_gsl_error)
{
    if (!lalGSLGlobalStatusPtr) {
        lalAbortHook
            ("Abort: function LALGSLErrorHandler, file %s, line %d, %s\n"
             "       Null global status pointer\n", __FILE__, __LINE__,
             "$Id$");
    }
    lalGSLGlobalStatusPtr->statusPtr = NULL;
    INITSTATUS(lalGSLGlobalStatusPtr);
    lalGSLGlobalStatusPtr->statusDescription = gsl_strerror(my_gsl_error);
    lalGSLGlobalStatusPtr->statusCode = my_gsl_error;
    lalGSLGlobalStatusPtr->file = file;
    lalGSLGlobalStatusPtr->line = line;
    LALError(lalGSLGlobalStatusPtr, reason);
    LALTrace(lalGSLGlobalStatusPtr, 1);
    return;
}
