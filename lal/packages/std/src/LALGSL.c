/**** <lalVerbatim file="LALGSLCV">
 * Author: Creighton, J. D. E.
 * $Id$
 **** </lalVerbatim> */

/* <lalLaTeX>
\subsection{Module \texttt{LALGSL.c}}

LAL GSL error handler.

\subsection*{Prototypes}
\begin{verbatim}
extern LALStatus *lalGSLGlobalStatusPtr;
#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
extern pthread_mutex_t lalGSLPthreadMutex;
#endif
\end{verbatim}
\input{LALGSLCP}
\idx[Variable]{lalGSLGlobalStatusPtr}
\idx[Variable]{lalGSLPthreadMutex}
\idx{LALGSLErrorHandler}

\subsection*{Description}

The function \verb+LALGSLErrorHandler+ is the standard GSL error handler
for GSL functions called within LAL.  Its function is to take the GSL
error code and information and translate them into equivalent aspects of
the LAL status structure.  The status structure that is currently identified
by \verb+lalGSLGlobalStatusPtr+ is populated.  This global variable is to
be set to the current status pointer prior to invocation of the GSL function.
In addition, the GSL error handler must be set to the LAL GSL error handler
prior to the invocation of the GSL function.  Both of these tasks can be
done with the macros provided in the header \texttt{LALGSL.h} \ref{s:LALGSL.h}.
However, doing so is not thread-safe.  Thus the macros use the mutex
\verb+lalGSLPthreadMutex+ to block other threads from making GSL calls
\emph{from within LAL functions} while one thread has set the GSL error
handler and global status pointer.  This mutex must also be used to block
any non-LAL GSL function calls from other threads or else they may be called
with the LAL GSL error handler in effect.

\vfill{\footnotesize\input{LALGSLCV}}
</lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALGSL.h>
#include <gsl/gsl_errno.h>

NRCSID( LALGSLC, "$Id$" );

LALStatus * lalGSLGlobalStatusPtr = NULL;
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
pthread_mutex_t lalGSLPthreadMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

/* <lalVerbatim file="LALGSLCP"> */
void
LALGSLErrorHandler(
    const char *reason,
    const char *file,
    int line,
    int my_gsl_error
    )
{ /* </lalVerbatim> */
  if ( ! lalGSLGlobalStatusPtr )
  {
    lalAbortHook( "Abort: function LALGSLErrorHandler, file %s, line %d, %s\n"
                  "       Null global status pointer\n",
                  __FILE__, __LINE__, LALGSLC ); 
  }
  lalGSLGlobalStatusPtr->statusPtr = NULL;
  INITSTATUS( lalGSLGlobalStatusPtr, "LALGSLErrorHandler", LALGSLC );
  lalGSLGlobalStatusPtr->statusDescription = gsl_strerror( my_gsl_error );
  lalGSLGlobalStatusPtr->statusCode        = my_gsl_error;
  lalGSLGlobalStatusPtr->file              = file;
  lalGSLGlobalStatusPtr->line              = line;
  LALError( lalGSLGlobalStatusPtr, reason );
  LALTrace( lalGSLGlobalStatusPtr, 1 );
  return;
}
