/****************************** <lalVerbatim file="LALGSLHV">
Author: Creighton, J. D. E.
$Id$
******************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALGSL.h}}
\label{s:LALGSL.h}

Provides macros for integrating the GSL error handler with
the LAL status structure.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALGSL.h>
\end{verbatim}

\noindent This header provides macros and functions for tracking and
reporting the runtime status of a GSL calls.  The intent is
simultaneously to standardize the error reporting, and to make the
reporting as transparent as possible to people coding individual
routines.

\emph{Please always use these macros when making a GSL call
within LAL.  This will ensure that the LAL functions always have the
same behaviour and will also ensure that the LAL functions are reenterant
and threadsafe (when LAL is configured appropriately).}

\subsubsection{GSL function calls}

The behaviour of GSL functions depends on the error handler that has been
assigned.  In order that LAL functions always have the same behaviour, it
is necessary to use a LAL-specific GSL error handler.  This error handler
populates a LAL status structure with the GSL error message and code so that
GSL functions behave much the same way as LAL functions.  After the GSL
functions are called, the error handler needs to be restored to the original
handler so that the program calling the LAL routine has the same error handler
after the LAL function was called as it did before the LAL function was called.

This module provides a simple set of macros and the default LAL GSL error
handler.  The macros provide a standard way to assign the LAL GSL error
handler before a GSL function call and to restore the original handler after
the call.

Note that changing the GSL error handler is \emph{not} a thread-safe
action.  Therefore it is necessary to block other threads from performing
GSL function calls while one thread has changed the handler.  These macros
ensure that such blocking is done for GSL function calls
\emph{within other LAL routines} if LAL is configured with the
\verb+--enable-pthread-lock+ flag.  See below for instructions on how
to make other GSL function calls outside LAL thread-safe when used with LAL.

\begin{verbatim}
ATTATCHSTATUSPTR( status );
CALLGSL( gsl_function( x ), status );
CHECKSTATUSPTR( status );
DETATCHSTATUSPTR( status );
\end{verbatim}
Note that the LAL function must attach (and detach) a status pointer as if
a LAL routine were called.
Note also that you need to use the \verb+CHECKSTATUSPTR+ macro to check
the status of the call.  The equivalent to the \verb+TRY+ macro for GSL
functions is the \verb+TRYGSL+ macro, which is used as follows:
\begin{verbatim}
ATTATCHSTATUSPTR( status );
TRYGSL( gsl_function( x ), status );
DETATCHSTATUSPTR( status );
\end{verbatim}

If you are using GSL functions both in LAL and in the calling program, and
you are worried about thread-safety, the GSL function calls outside of LAL
need to be blocked so that they do not access the GSL error handler while
it has been changed to the LAL GSL error handler in a LAL function.  To do
this, you need to do the following:
\begin{verbatim}
#include<lal/LALGSL.h>
...
LALGSL_PTHREAD_MUTEX_LOCK;
gsl_function( x );
LALGSL_PTHREAD_MUTEX_UNLOCK;
\end{verbatim}
This ensures that \verb+gsl_function+ is not called while a LAL routine
is calling a GSL function in a different thread.  You can do this even if
you don't always run your code with multiple threads.  If you configure LAL
without the \verb+--enable-pthread-lock+ flag, the macros
\verb+LALGSL_PTHREAD_MUTEX_LOCK+ and \verb+LALGSL_PTHREAD_MUTEX_UNLOCK+
do nothing.

\vfill{\footnotesize\input{LALGSLHV}}
\newpage\input{LALGSLC}
\newpage\input{LALGSLTestC}
</lalLaTeX> */

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
#include <lal/LALRCSID.h>

#include <lal/XLALGSL.h>

#include <gsl/gsl_errno.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( LALGSLH, "$Id$" );

extern LALStatus * lalGSLGlobalStatusPtr;
void
LALGSLErrorHandler(
    const char *reason,
    const char *file,
    int line,
    int errnum
    );

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

#ifdef  __cplusplus
}
#endif

#endif /* _LALGSL_H */
