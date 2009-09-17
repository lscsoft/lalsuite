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

/************************************ <lalVerbatim file="LALErrorHV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALError.h}}
\label{s:LALError.h}

Provides routines to report and handle errors.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALError.h>
\end{verbatim}

\noindent This header covers routines that print status messages, and
that allow functions to abort.

\vfill{\footnotesize\input{LALErrorHV}}
\newpage\input{LALErrorC}

</lalLaTeX> */

#ifndef _LALERROR_H
#define _LALERROR_H

#include <lal/LALDatatypes.h>

#include <lal/XLALError.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LALERRORH, "$Id$" );

/* lalDebugLevel bit field values: */
enum
{
  LALNDEBUG  = 0,
  LALERROR   = 1,
  LALWARNING = 2,
  LALINFO    = 4,
  LALTRACE   = 8,
  LALMEMINFO = 16,
  LALNMEMDBG = 32,
  LALNMEMPAD = 64,
  LALNMEMTRK = 128,
  LALMEMDBG  = 16384 /* convenience: don't combine with other bits */
};

/* composite lalDebugLevels: */
enum { LALMSGLVL1  = LALERROR };
enum { LALMSGLVL2  = LALERROR | LALWARNING };
enum { LALMSGLVL3  = LALERROR | LALWARNING | LALINFO };
enum { LALMEMTRACE = LALTRACE | LALMEMINFO };
enum { LALALLDBG   = ~( LALNMEMDBG | LALNMEMPAD | LALNMEMTRK ) };

extern int  ( *lalRaiseHook )( int, const char *, ... );
extern void ( *lalAbortHook )( const char *, ... );

int
LALPrintError( const char *fmt, ... );

int
LALRaise( int sig, const char *fmt, ... );

void
LALAbort( const char *fmt, ... );

int
LALError( LALStatus *status, const char *statement );

int
LALWarning( LALStatus *status, const char *warning );

int
LALInfo( LALStatus *status, const char *info );

int
LALTrace( LALStatus *status, int exitflg );

int
LALInitStatus( LALStatus *status, const char *function, const char *id,
               const char *file, const int line );

int
LALPrepareReturn( LALStatus *status, const char *file, const int line );

int
LALAttatchStatusPtr( LALStatus *status, const char *file, const int line );

int
LALDetatchStatusPtr( LALStatus *status, const char *file, const int line );

int
LALPrepareAbort( LALStatus *status, const INT4 code, const char *mesg,
                 const char *file, const int line );

int
LALPrepareAssertFail( LALStatus *status, const INT4 code, const char *mesg,
                      const char *statement, const char *file,
                      const int line );

int
LALCheckStatusPtr( LALStatus *status, const char *statement, const char *file,
                   const int line );

void
FREESTATUSPTR( LALStatus *status );

void
REPORTSTATUS( LALStatus *status );

#ifdef NDEBUG

#define LALError( statusptr, statement ) 0
#define LALWarning( statusptr, warning ) 0
#define LALInfo( statusptr, info )       0
#define LALTrace( statusptr, exitflg )   0

#else

#ifndef NOLALMACROS

#define LALError( statusptr, statement )                                    \
  ( lalDebugLevel & LALERROR ?                                                 \
    LALPrintError( "Error[%d] %d: function %s, file %s, line %d, %s\n"       \
        "        %s %s\n", (statusptr)->level, (statusptr)->statusCode,      \
        (statusptr)->function, (statusptr)->file, (statusptr)->line,        \
        (statusptr)->Id, (statement),                    \
        (statusptr)->statusDescription ) : 0 )

#define LALWarning( statusptr, warning )                                    \
  ( lalDebugLevel & LALWARNING ?                                               \
    LALPrintError( "Warning[%d]: function %s, file %s, line %d, %s\n"        \
        "        %s\n", (statusptr)->level, (statusptr)->function,          \
        (statusptr)->file, (statusptr)->line, (statusptr)->Id, (warning) )  \
      : 0 )

#define LALInfo( statusptr, info )                                          \
  ( lalDebugLevel & LALINFO ?                                                  \
    LALPrintError( "Info[%d]: function %s, file %s, line %d, %s\n"          \
        "        %s\n", (statusptr)->level, (statusptr)->function,          \
        (statusptr)->file, (statusptr)->line, (statusptr)->Id, (info) )     \
      : 0 )

#define LALTrace( statusptr, exitflg ) \
  ( lalDebugLevel & LALTRACE ? \
    LALPrintError( "%s[%d]: function %s, file %s, line %d, %s\n",      \
        (exitflg) ? "Leave" : "Enter", (statusptr)->level, \
        (statusptr)->function, (statusptr)->file, (statusptr)->line, \
        (statusptr)->Id )     \
      : 0 )

#endif /* NOLALMACROS */

#endif /* NDEBUG */

#ifdef  __cplusplus
}
#endif

#endif /* _LALERROR_H */
