/************************************ <lalVerbatim file="LALErrorHV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALError.h}}
\label{s:LALError.h}

Provides routines to report and handle errors.

\subsection*{Synopsis}
\begin{verbatim}
#include "LALError.h"
\end{verbatim}

\noindent This header covers routines that print status messages, and
that allow functions to abort.

\vfill{\footnotesize\input{LALErrorHV}}
\newpage\input{LALErrorC}

</lalLaTeX> */

#ifndef _LALERROR_H
#define _LALERROR_H

#include "LALDatatypes.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LALERRORH, "$Id$" );

/* debuglevel bit field values: */
enum
{
  LALNDEBUG  = 0,
  LALERROR   = 1,
  LALWARNING = 2,
  LALINFO    = 4,
  LALTRACE   = 8,
  LALMEMINFO = 16,
  LALMEMDBG  = ~( -1U >> 1 ) /* any non-zero value that doesn't interfere with
                                other bits: meaningless to combine this with
                                any other bit */
};

/* composite debuglevels: */
enum { LALMSGLVL1  = LALERROR };
enum { LALMSGLVL2  = LALERROR | LALWARNING };
enum { LALMSGLVL3  = LALERROR | LALWARNING | LALINFO };
enum { LALMEMTRACE = LALTRACE | LALMEMINFO };
enum { LALALLDBG   = ~0 };

int
LALPrintError( const char *fmt, ... );

void
LALAbort( const char *fmt, ... );

int
LALError( Status *status, const char *statement );

int
LALWarning( Status *status, const char *warning );

int
LALInfo( Status *status, const char *info );

int
LALTrace( Status *status, int exit );

int
LALInitStatus( Status *status, const char *function, const char *id,
               const char *file, const int line );

int
LALPrepareReturn( Status *status, const char *file, const int line );

int
LALAttatchStatusPtr( Status *status, const char *file, const int line );

int
LALDetatchStatusPtr( Status *status, const char *file, const int line );

int
LALPrepareAbort( Status *status, const INT4 code, const char *mesg,
                 const char *file, const int line );

int
LALPrepareAssertFail( Status *status, const INT4 code, const char *mesg,
                      const char *statement, const char *file,
                      const int line );

int
LALCheckStatusPtr( Status *status, const char *statement, const char *file,
                   const int line );

void
FREESTATUSPTR( Status *status );

void
REPORTSTATUS( Status *status );

#ifdef NDEBUG

#define LALError( statusptr, statement ) 0
#define LALWarning( statusptr, warning ) 0
#define LALInfo( statusptr, info )       0
#define LALTrace( statusptr, exit )      0

#else

#ifndef NOLALMACROS

#define LALError( statusptr, statement )                                    \
  ( debuglevel & LALERROR ?                                                 \
    LALPrintError( "Error[%d] %d: function %s, file %s, line %d, %s\n"       \
        "        %s %s\n", (statusptr)->level, (statusptr)->statusCode,      \
        (statusptr)->function, (statusptr)->file, (statusptr)->line,        \
        (statusptr)->Id, (statement) ? (statement) : "",                    \
        (statusptr)->statusDescription ) : 0 )

#define LALWarning( statusptr, warning )                                    \
  ( debuglevel & LALWARNING ?                                               \
    LALPrintError( "Warning[%d]: function %s, file %s, line %d, %s\n"        \
        "        %s\n", (statusptr)->level, (statusptr)->function,          \
        (statusptr)->file, (statusptr)->line, (statusptr)->Id, (warning) )  \
      : 0 )

#define LALInfo( statusptr, info )                                          \
  ( debuglevel & LALINFO ?                                                  \
    LALPrintError( "Info[%d]: function %s, file %s, line %d, %s\n"          \
        "        %s\n", (statusptr)->level, (statusptr)->function,          \
        (statusptr)->file, (statusptr)->line, (statusptr)->Id, (info) )     \
      : 0 )

#define LALTrace( statusptr, exit ) \
  ( debuglevel & LALTRACE ? \
    LALPrintError( "%s[%d]: function %s, file %s, line %d, %s\n",      \
        (exit) ? "Leave" : "Enter", (statusptr)->level, \
        (statusptr)->function, (statusptr)->file, (statusptr)->line, \
        (statusptr)->Id )     \
      : 0 )

#endif /* NOLALMACROS */

#endif /* NDEBUG */

#ifdef  __cplusplus
}
#endif

#endif /* _LALERROR_H */
