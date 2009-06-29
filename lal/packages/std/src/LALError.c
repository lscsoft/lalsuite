/*
*  Copyright (C) 2007 Jolien Creighton, Bernd Machenschalk
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

/************************************ <lalVerbatim file="LALErrorCV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{LALError.c}}
\label{ss:LALError.c}

Error handling routines for LAL.  These should \emph{not} be invoked in
production code, except in very specific circumstances.

\subsubsection*{Prototypes}
\input{LALErrorCP}
\idx[Variable]{lalRaiseHook}
\idx[Variable]{lalAbortHook}
\idx{LALPrintError()}
\idx{LALRaise()}
\idx{LALAbort()}
\idx{LALError()}
\idx{LALWarning()}
\idx{LALInfo()}
\idx{LALTrace()}

\subsubsection*{Description}

These functions cause LAL to print status messages and perform basic
error handling.  Their implementation is quite simple but may be
altered in the future to provide reasonable behaviour when integrated
with other systems (e.g., LDAS).  As a general rule,
\verb@LALWarning()@ and \verb@LALInfo()@ are the only routines that
programmers should use in their own modules; the other routines are
used internally by LAL.  Descriptions of the individual functions are
as follows.

\paragraph{\texttt{LALPrintError()}} prints a formatted string to some
designated output device (usually the \verb@stderr@ stream), returning
the number of characters printed, or negative if an error occurred.
The format of the argument list is the same as for the standard C
routine \verb@printf()@.  By funneling all LAL error printing through
this one routine, it is easier to adapt LAL to implementations that
have particular I/O or error-logging requirements.  Most LAL routines
should use \verb@LALError()@, \verb@LALWarning()@, and
\verb@LALInfo()@ to report their status, rather than calling
\verb@LALPrintError()@ directly.

\paragraph{\texttt{LALRaise()}} prints a formatted string to an error
logging device, as above, and then raises the requested signal.
Standard LAL routines should \emph{not} terminate execution, but should
instead return control to the calling routine, reporting errors through their
\verb@LALStatus@ structure.  Programmers should never
invoke \verb@LALRaise()@ explicitly.
A hook to a \verb@LALRaise()@-type function, \verb@lalRaiseHook@, is provided,
should the user wish to change the default behavior of \verb@LALRaise()@
(i.e., the LAL library always uses \verb@lalRaiseHook@ rather than
\verb@LALRaise@, but \verb@lalRaiseHook@ is set to \verb@LALRaise@ by
default).

\paragraph{\texttt{LALAbort()}} prints a formatted string to an error
logging device, as above, and then terminates program execution.
Usually this is done by raising a \verb@SIGABRT@ signal, but this can
change in implementations that have different requirements.  Standard
LAL routines should \emph{not} terminate execution, but should instead
return control to the calling routine, reporting errors through their
\verb@LALStatus@ structure.  The exception is when a function receives a
\verb@NULL@ status pointer, in which case it has no option but to
abort.  This is done automatically by the \verb@INITSTATUS()@ macro
(see \verb@LALStatusMacros.h@), so programmers should never need to
invoke \verb@LALAbort()@ explicitly.
A hook to a \verb@LALAbort()@-type function, \verb@lalAbortHook@, is provided,
should the user wish to change the default behavior of \verb@LALAbort()@
(i.e., the LAL library always uses \verb@lalAbortHook@ rather than
\verb@LALAbort@, but \verb@lalAbortHook@ is set to \verb@LALAbort@ by
default).

\paragraph{\texttt{LALError()}} prints the \verb@statement@
string to the error log, provided that the value of the global
\verb@lalDebugLevel@ is set to allow error messages.  It returns the
number of characters printed.  This is the standard LAL routine for
printing error messages.  However, \verb@LALError()@ is called
automatically by the status-handling macros (see
\verb@LALStatusMacros.h@) whenever a LAL function returns with
non-zero error code.  Since an error is, by definition, a condition
that would cause a routine to terminate abnormally, LAL programmers
will generally not have to call \verb@LALError()@ explicitly.

\paragraph{\texttt{LALWarning()}} prints the \verb@warning@
string to the error log, provided that the value of the global
\verb@lalDebugLevel@ is set to allow warning messages.  It returns the
number of characters printed.  A warning message is less serious than
an error message: it indicates that computation is proceeding
successfully, but with unusual or unexpected behaviour that may
invalidate the results of the computation.

\paragraph{\texttt{LALInfo()}} prints the \verb@info@
string to the error log, provided that the value of the global
\verb@lalDebugLevel@ is set to allow information messages.  It returns
the number of characters printed.  An information message indicates
that a computation is proceding normally, and simply provides
additional information about its progress.

\paragraph{\texttt{LALTrace()}} prints a message providing
information, taken from the \verb@status@ structure, about the
function currently being executed; it is used to track the progress of
execution through nested function calls.  It returns the number of
characters printed.  The message begins with the word \verb@Enter@ (if
\verb@exitflg@ = 0) or \verb@Leave@ (if \verb@exitflg@ $\neq0$), to indicate
whether the flow of execution has just entered or is about to leave
the function.  Tracking information is printed only if the value of
the global \verb@lalDebugLevel@ is set to allow it.  \verb@LALTrace()@ is
called automatically by the status macros when entering or leaving a
function (see \verb@LALStatusMacros.h@), so LAL programmers need never
invoke it explicitly.

\subsubsection*{Algorithm}

The functions \verb@LALError()@, \verb@LALWarning()@,
\verb@LALInfo()@, and \verb@LALTrace()@ print status messages
depending on the value of the global \verb@lalDebugLevel@.  Specifically,
each type of status message is associated with a particular bit in
\verb@lalDebugLevel@.  If the value of the bit is 1, that type status
message will be printed; if it is 0, that type of message will be
suppressed.  See the documentation in \verb@LALStatusMacros.h@ for
information about how to set the value of \verb@lalDebugLevel@.

These four functions are also suppressed if a module is compiled with
the \verb@NDEBUG@ flag set.  In this case, however, the function calls
are actually \emph{removed} from the object code (i.e.\ they are
replaced with the integer 0, representing their return value).  This
is used to generate streamlined production code.  Again, see the
\verb@LALStatusMacros.h@ documentation for more discussion of this
compilation flag.

\subsubsection*{Macro replacement functions}

When a LAL module is compiled with the flag \verb@NOLALMACROS@ set,
the usual status-handling macros defined in \verb@LALStatusMacros.h@
are replaced with function calls to specialized support functions that
perform the same operations.  These functions are necessarily global
in scope, and so we provide their prototype declarations below.
However, they will never be invoked explicitly in any LAL function, so
we will not bother with additional usage information.

\vspace{1ex}
\input{LALErrorCP2}
\idx{LALInitStatus()}
\idx{LALPrepareReturn()}
\idx{LALAttatchStatusPtr()}
\idx{LALDetatchStatusPtr()}
\idx{LALPrepareAbort()}
\idx{LALPrepareAssertFail()}
\idx{LALCheckStatusPtr()}
\idx{FREESTATUSPTR()}
\idx{REPORTSTATUS()}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALErrorCV}}

</lalLaTeX> */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <signal.h>
#include <lal/LALMalloc.h>
#include <lal/LALError.h>

#undef LALError
#undef LALWarning
#undef LALInfo
#undef LALTrace

#ifdef LAL_NDEBUG
/* sorry, but it looks like we have no other choice than a special treatment for Einstein@Home here.
   we can neither afford to compile LAL in DEBUG mode nor losing valuable error information when
   distributing executables to unsafe public computers.                         Bernd Machenschalk */
#ifndef EAH_BOINC
#define vfprintf( stream, fmt, ap ) 0
#endif
#endif

NRCSID( LALERRORC, "$Id$" );

extern int lalDebugLevel;

/* <lalVerbatim file="LALErrorCP"> */
int
LALPrintError( const char *fmt, ... )
{ /* </lalVerbatim> */
  int n;
  va_list ap;
  va_start( ap, fmt );
  n = vfprintf( stderr, fmt, ap );
  va_end( ap );
  return n;
}

/* <lalVerbatim file="LALErrorCP"> */
int ( *lalRaiseHook )( int, const char *, ... ) = LALRaise;
int
LALRaise( int sig, const char *fmt, ... )
{ /* </lalVerbatim> */
  va_list ap;
  va_start( ap, fmt );
  vfprintf( stderr, fmt, ap );
  va_end( ap );
  return raise( sig );
}

/* <lalVerbatim file="LALErrorCP"> */
void ( *lalAbortHook )( const char *, ... ) = LALAbort;
void
LALAbort( const char *fmt, ... )
{ /* </lalVerbatim> */
  va_list ap;
  va_start( ap, fmt );
  (void) vfprintf( stderr, fmt, ap );
  va_end( ap );
  abort();
}


/* <lalVerbatim file="LALErrorCP"> */
int
LALError( LALStatus *status, const char *statement )
{ /* </lalVerbatim> */
  int n = 0;
  if ( lalDebugLevel & LALERROR )
  {
    n = LALPrintError( "Error[%d] %d: function %s, file %s, line %d, %s\n"
        "        %s %s\n", status->level, status->statusCode,
        status->function, status->file, status->line, status->Id,
        statement ? statement : "", status->statusDescription );
  }
  return n;
}


/* <lalVerbatim file="LALErrorCP"> */
int
LALWarning( LALStatus *status, const char *warning )
{ /* </lalVerbatim> */
  int n = 0;
  if ( lalDebugLevel & LALWARNING )
  {
    n = LALPrintError( "Warning[%d]: function %s, file %s, line %d, %s\n"
        "        %s\n", status->level, status->function, status->file,
        status->line, status->Id, warning );
  }
  return n;
}


/* <lalVerbatim file="LALErrorCP"> */
int
LALInfo( LALStatus *status, const char *info )
{ /* </lalVerbatim> */
  int n = 0;
  if ( lalDebugLevel & LALINFO )
  {
    n = LALPrintError( "Info[%d]: function %s, file %s, line %d, %s\n"
        "        %s\n", status->level, status->function, status->file,
        status->line, status->Id, info );
  }
  return n;
}


/* <lalVerbatim file="LALErrorCP"> */
int
LALTrace( LALStatus *status, int exitflg )
{ /* </lalVerbatim> */
  int n = 0;
  if ( lalDebugLevel & LALTRACE )
  {
    n = LALPrintError( "%s[%d]: function %s, file %s, line %d, %s\n",
        exitflg ? "Leave" : "Enter", status->level, status->function,
        status->file, status->line, status->Id );
  }
  return n;
}


#ifdef LAL_NDEBUG

#define LALError( statusptr, statement ) 0
#define LALWarning( statusptr, warning ) 0
#define LALInfo( statusptr, info )       0
#define LALTrace( statusptr, exitflg )   0

#endif

/* <lalVerbatim file="LALErrorCP2"> */
int
LALInitStatus( LALStatus *status, const char *function, const char *id,
	       const char *file, const int line )
{ /* </lalVerbatim> */
  int exitcode = 0;
  if ( status )
  {
    INT4 level = status->level;
    exitcode = status->statusPtr ? 1 : 0;
    memset( status, 0, sizeof( LALStatus ) ); /* possible memory leak */
    status->level    = level > 0 ? level : 1;
    status->Id       = id;
    status->function = function;
    status->file     = file;
    status->line     = line;
    (void) LALTrace( status, 0 );
    if ( exitcode )
    {
      LALPrepareAbort( status, -2, "INITSTATUS: non-null status pointer",
                       file, line );
    }
    else if ( xlalErrno )
    {
      LALPrepareAbort( status, -16, "INITSTATUS: non-zero xlalErrno",
                       file, line );
      exitcode = 1;
    }
  }
  else
  {
    lalAbortHook( "Abort: function %s, file %s, line %d, %s\n"
                  "       Null status pointer passed to function\n",
                  function, file, line, id );
  }
  return exitcode;
}


/* <lalVerbatim file="LALErrorCP2"> */
int
LALPrepareReturn( LALStatus *status, const char *file, const int line )
{ /* </lalVerbatim> */
  status->file = file;
  status->line = line;
  if ( status->statusCode )
  {
    (void) LALError( status, "RETURN:" );
  }
  (void) LALTrace( status, 1 );
  if ( xlalErrno )
  {
    LALPrepareAbort( status, -32, "RETURN: untrapped XLAL error",
        file, line );
  }
  return 1;
}


/* <lalVerbatim file="LALErrorCP2"> */
int
LALAttatchStatusPtr( LALStatus *status, const char *file, const int line )
{ /* </lalVerbatim> */
  int exitcode = 0;
  if ( status->statusPtr )
  {
    LALPrepareAbort( status, -2, "ATTATCHSTATUSPTR: non-null status pointer",
                     file, line );
    exitcode = 1;
  }
  else
  {
    status->statusPtr = (LALStatus *) LALCalloc( 1, sizeof( LALStatus ) );
    if ( !status->statusPtr )
    {
      LALPrepareAbort( status, -4, "ATTATCHSTATUSPTR: memory allocation error",
                       file, line );
      exitcode = 1;
    }
    else
    {
      status->statusPtr->level = status->level + 1;
    }
  }
  return exitcode;
}


/* <lalVerbatim file="LALErrorCP2"> */
int
LALDetatchStatusPtr( LALStatus *status, const char *file, const int line )
{ /* </lalVerbatim> */
  int exitcode = 0;
  if ( status->statusPtr )
  {
    FREESTATUSPTR( status );
    status->statusCode = 0;
    status->statusDescription = NULL;
  }
  else
  {
    LALPrepareAbort( status, -8, "DETATCHSTATUSPTR: null status pointer",
                     file, line );
    exitcode = 1;
  }
  return exitcode;
}


/* <lalVerbatim file="LALErrorCP2"> */
int
LALPrepareAbort( LALStatus *status, const INT4 code, const char *mesg,
		 const char *file, const int line )
{ /* </lalVerbatim> */
  if ( status->statusPtr )
  {
    FREESTATUSPTR( status );
  }
  status->file              = file;
  status->line              = line;
  status->statusCode        = code;
  status->statusDescription = mesg;
  if ( code )
  {
    (void) LALError( status, "ABORT:" );
  }
  (void) LALTrace( status, 1 );
  return 1;
}


/* <lalVerbatim file="LALErrorCP2"> */
int
LALPrepareAssertFail( LALStatus *status, const INT4 code, const char *mesg,
		      const char *statement, const char *file,
		      const int line )
{ /* </lalVerbatim> */
  if ( status->statusPtr )
  {
    FREESTATUSPTR( status );
  }
  status->file              = file;
  status->line              = line;
  status->statusCode        = code;
  status->statusDescription = mesg;
  (void) LALError( status, statement );
  (void) LALTrace( status, 1 );
  return 1;
}


/* <lalVerbatim file="LALErrorCP2"> */
int
LALCheckStatusPtr( LALStatus *status, const char *statement, const char *file,
		   const int line )
{ /* </lalVerbatim> */
  if ( status->statusPtr->statusCode )
  {
    status->file              = file;
    status->line              = line;
    status->statusCode        = -1;
    status->statusDescription = "Recursive error";
    (void) LALError( status, statement );
    (void) LALTrace( status, 1 );
    return 1;
  }
  return 0;
}



/*
 * This function is somewhat dangerous: need to check to see
 * if status->statusPtr is initially null before calling FREESTATUSPTR
 */
/* <lalVerbatim file="LALErrorCP2"> */
void
FREESTATUSPTR( LALStatus *status )
{ /* </lalVerbatim> */
  do
  {
    LALStatus *next = status->statusPtr->statusPtr;
    LALFree( status->statusPtr );
    status->statusPtr = next;
  }
  while ( status->statusPtr );
  return;
}


/* <lalVerbatim file="LALErrorCP2"> */
void
REPORTSTATUS( LALStatus *status )
{ /* </lalVerbatim> */
  LALStatus *ptr;
  for ( ptr = status; ptr ; ptr = ptr->statusPtr )
  {
    LALPrintError( "\nLevel %i: %s\n", ptr->level, ptr->Id );
    if ( ptr->statusCode )
    {
      LALPrintError( "\tStatus code %i: %s\n", ptr->statusCode,
                     ptr->statusDescription );
    }
    else
    {
      LALPrintError( "\tStatus code 0: Nominal\n" );
    }
    LALPrintError( "\tfunction %s, file %s, line %i\n",
                   ptr->function, ptr->file, ptr->line );
  }
  return;
}

