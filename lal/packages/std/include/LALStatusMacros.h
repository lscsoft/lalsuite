/****************************** <lalVerbatim file="LALStatusMacrosHV">
Author: Creighton, J. D. E. and Creighton, T. D.
$Id$
******************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALStatusMacros.h}}
\label{s:LALStatusMacros.h}

Provides macros for handling the LAL status structure.

\subsection*{Synopsis}
\begin{verbatim}
#include "LALStatusMacros.h"
\end{verbatim}

\noindent This header defines a number of macros for interacting with
the LAL status structure \verb@Status@.  The intent is simultaneously
to standardize the error reporting, and to make the reporting as
transparent as possible to people coding individual routines.

The following summarized everything the common programmer needs to
know in order to follow LAL standard error reporting.  It can be
treated as a primer on LAL coding conventions.

\subsection*{The \texttt{Status} structure}

The \verb@Status@ structure is the standard LAL structure for
reporting progress and errors; it is described in
Sec.~\ref{ss:status-structure} of the header \verb@LALDatatypes.h@.
For completeness, we explain the fields of this structure below:
\begin{description}
\item[\texttt{INT4 statusCode}] A code indicating the exit status of a
function.  0 represents a normal exit.  Negative values are reserved
for certain standard error types.  The authors of individual functions
should assign positive values to the various ways in which their code
can fail.
\item[\texttt{const CHAR *statusDescription}] An explanatory string
corresponding to the numerical status code.
\item[\texttt{volatile const CHAR *Id}] A character string identifying
the source file and version number of the function being reported on.
\item[\texttt{const CHAR *function}] The name of the function.
\item[\texttt{const CHAR *file}] The file name of the \verb@.c@ file
containing the function code.
\item[\texttt{INT4 line}] The line number in the \verb@.c@ file of the
instruction where any error was reported.
\item[\texttt{Status *statusPtr}] A recursive pointer to another
status pointer.  This structure is used to report an error in a
subroutine of the current function.  Thus if an error occurs in a
deeply-nested routine, the status structure returned to the main
program will be the head of a linked list of status structures, one
for each nested level, with the tail structure reporting the actual
error that caused the overlying routines to fail.
\item[\texttt{INT4 level}] The nested-function level where any error
was reported.
\end{description}
In almost all circumstances the programmer will \emph{not} have to
access this structure directly, relying instead on the macros defined
in this header.  The exception is the \verb@statusCode@ field, which
the programmer may want to query directly.

The \verb@statusCode@ field is set to a nonzero value any time an
error condition arises that would lead to abnormal termination of the
current function.  Programmers can assign positive error codes to the
various types of error that may be encountered in their routines.
Additionally, the following following status codes are reserved to
report certain standard conditions:

\begin{center}
\begin{tabular}{|cp{3.5cm}p{6.5cm}|}
\hline
Code & Message & Explanation \\
\hline

\tt 0 & & Nominal execution; the function returned
successfully. \\

\tt -1 & \vspace{-1.4ex}\tt Recursive error & The function aborted due
to failure of a subroutine. \\

\tt -2 & \vspace{-1.4ex}\tt INITSTATUS: non-null status pointer & The
status structure passed to the function had a non-\verb@NULL@
\verb@statusPtr@ field, which blocks the function from calling
subroutines (it is symptomatic of something screwy going on in the
calling routine). \\

\tt -4 & \vspace{-1.4ex}\tt ATTATCHSTATUSPTR: memory allocation error
& The function was unable to allocate a \verb@statusPtr@ field to pass
down to a subroutine. \\

\tt -8 & \vspace{-1.4ex}\tt DETATCHSTATUSPTR: null status pointer &
The \verb@statusPtr@ field could not be deallocated at the end of all
subroutine calls; one of the subroutines must have lost it or set it
to \verb@NULL@. \\

\hline
\end{tabular}
\end{center}

\subsection*{LAL function calls}

All functions should have return type void.  The first argument of any
function should be a pointer to a structure of type \verb@Status@.
Thus:
\begin{verbatim} 
void MyFunction( Status *stat, ... )
\end{verbatim}
Since the function has no return code, it must report all errors or
failure through the status structure.  A function that is passed a
\verb@NULL@ pointer in place of the status pointer should terminate
the program with a \verb@SIGABRT@ signal, as this is its only way to
report the error.  However, this is the only circumstance under which
a function sould deliberately raise a signal.  In all other cases the
error should be trapped, reported in the status structure, and control
returned to the calling routine.

\subsection*{Assigning an RCS \texttt{\$Id\$} string}

Every source file should have a unique character string identifying
that version of that file.  The standard convention, for a file
\verb@MyFile.c@, is to declare a string \verb@MYFILEC@ at the top of
the module using the macro \verb@NRCSID()@ (defined in the include
file \verb@LALRCSID.h@):

\vspace{2ex}
\noindent\texttt{NRCSID( MYFILEC, \$Id\$ );}
\vspace{2ex}

\noindent where \texttt{\$Id\$} is expanded by RCS to give the full
name and version number of the source file.

\subsection*{Initializing the status structure}

The first instruction in any function, after variable declarations,
should be the macro \verb@INITSTATUS()@, which takes three arguments:
the function's status pointer, the function name (a string literal)
and the module's RCS \texttt{\$Id\$} string.
\begin{verbatim}
INITSTATUS( stat, "MyFunction", MYFILEC );
\end{verbatim}
This macro checks that a valid status pointer has been passed to the
function, and if so, initializes the other fields to indicate (by
default) nominal execution.  If \verb@stat@ is null, the macro causes
the program to terminate with a \verb@SIGABRT@ signal, as described
above.

\subsection*{Normal return from a function}

Upon completion, the function should issue the macro \verb@RETURN()@,
which takes one argument: the function's status pointer.
\begin{verbatim}
RETURN( stat );
\end{verbatim}
This takes the place of any return statements, and may also log status
reports to some suitable log file (often \verb@stderr@), depending on
implementation and the value of a global \verb@debuglevel@ parameter.
Typically \verb@RETURN()@ is used only for successful completion, with
other macros \verb@ABORT()@, \verb@ASSERT()@, \verb@CHECKSTATUSPTR()@,
and \verb@TRY()@ being used to report failure.  However, it is
possible for the programmer to assign the fields of \verb@stat@ by
hand, and then issue \verb@RETURN()@.

\subsection*{Abnormal return from a function}

The standard method to terminate a function unsuccessfully is with the
\verb@ABORT()@ macro, which takes three arguments: the status pointer,
the status code, and the status description string.  Normally the
various error codes and descriptions will be constants defined in the
function's header file \verb@MyHeader.h@:
\begin{verbatim}
ABORT( stat, MYHEADER_EMYERR, MYHEADER_MSGEMYERR );
\end{verbatim}
where the error code \verb@MYHEADER_EMYERR@ and the error message
\verb@MYHEADER_MSGEMYERR@ are defined in \verb@MyHeader.h@.  This
standard LAL naming convention for error messages prevents namespace
conflicts between different header files.  Like \verb@RETURN()@,
\verb@ABORT()@ correctly handles any status logging required by the
implementation and the \verb@debuglevel@.  Note that \verb@ABORT()@
does \emph{not} raise a \verb@SIGABRT@ signal, but instead returns
control to the calling routine.

\subsection*{Error checking within a function}

Another way to indicate an unsuccessful termination is with the macro
\verb@ASSERT()@, which takes as arguments a test statement, a status
pointer, a status code, and a status description.  The statement
\verb@ASSERT(assertion,...);@ is in all ways equivalent to the
statement \verb@if(!assertion) ABORT(...);@, except on a failure the
\verb@ASSERT()@ macro will also log the failed assertion.  In the
above example, one might have:
\begin{verbatim}
ASSERT( assertion, stat, MYHEADER_EMYERR, MYHEADER_MSGEMYERR );
\end{verbatim}
Coding purists argue that \verb@ASSERT()@ should be used only to trap
coding errors rather than runtime errors, which would be trapped using
\verb@ABORT()@.  In other words, the assertion should always test true
in the final debugged program.  At present this coding practice is not
enforced by the LAL standard.  However, programmers should be aware of
the side effects of using \verb@ASSERT()@ to exit a function in normal
runtime.

For example, it is an error to allocate dynamic memory to local
variables in a function and then fail to free it before returning.
Thus, if you have dynamically allocated memory, you cannot then use
\verb@ASSERT()@ for runtime error checking, as this does not permit
you to free the memory before returning.  Instead, you must check the
assertion, and, if it fails, free the memory and call \verb@ABORT()@.

\subsection*{Calling subroutines}

If the function is to call other LAL functions as subroutines, four
more macros are used to report possible errors arising in these
routines.  The macros are \verb@ATTATCHSTATUSPTR()@,
\verb@DETATCHSTATUSPTR()@, \verb@CHECKSTATUSPTR()@, and \verb@TRY()@.
The usage of these macros is as follows.

\begin{enumerate}

\item First, before any subroutines are called, the function must call
the macro \verb@ATTATCHSTATUSPTR()@ which takes as its argument the
status pointer of the current function:
\begin{verbatim}
ATTATCHSTATUSPTR( stat );
\end{verbatim}
This allocates \verb@stat->statusPtr@, which is the status pointer
that will be handed down into any and all subroutines.  If the pointer
has already been allocated, \verb@ATTATCHSTATUSPTR()@ will abort, as
this is symptomatic of a coding error.

In most cases \verb@ATTATCHSTATUSPTR()@ need only be called once in a
given function, immediately after \verb@INITSTATUS()@, no matter how
many subroutine calls that function makes.  The exception is if the
function deals with (or ignores) errors reported by its subroutines.
In that case, the function should detatch the status pointer using
\verb@DETATCHSTATUSPTR()@ (below), and then re-attatch it.

The macro \verb@ATTATCHSTATUSPTR()@ sets the status code to be $-1$
and the status message to be \verb@Recursive error@.  These flags are
unset when \verb@DETATCHSTATUSPTR()@ (below) is called.  This is so
that a use of \verb@RETURN()@ prior to detatching the status pointer
will yield an error.

\item When a subroutine is called, it should be handed the
\verb@statusPtr@ field of the calling function's status structure, to
report its own errors.  The calling function should test the returned
status code, and either attempt to deal with any abnormal returns, or
abort with status code $-1$.  The macro \verb@CHECKSTATUSPTR()@
simplifies the latter case.  It takes one arguments: the status
pointer of the current function (not the subroutine).
\begin{verbatim}
MySubroutine( stat->statusPtr, ... );
CHECKSTATUSPTR( stat );
\end{verbatim}
The \verb@TRY()@ macro is a somewhat more streamlined approach but
with equivalent results.  It takes two arguments.  The first is the
subroutine call, and the second is the status pointer.  Thus:
\begin{verbatim}
TRY( MySubroutine( stat->statusPtr, ... ), stat );
\end{verbatim}
The only practical difference between these two approaches is that
\verb@TRY()@ also reports the name of the failed subroutine call when
logging errors.

Similar caveats apply when using \verb@CHECKSTATUSPTR()@ and
\verb@TRY()@ as when using \verb@ASSERT()@, in that these macros can
force an immediate return with no additional housekeeping
instructions.  For instance, if you have dynamically-allocated local
memory, you should explicitly check the \verb@statusPtr->statusCode@
field to see if a subroutine failed, then free the memory and call
\verb@ABORT()@ to exit.

If the calling routine attempts to work around an error reported from
a subroutine, and the attempt fails, the routine should \emph{not} use
\verb@CHECKSTATUSPTR()@ to exit with status code $-1$.  Instead, it
should call \verb@ABORT()@ with an appropriate (positive) code and
message to indicate how the attempted workaround failed.

\item After all subroutines have been called, but before any
\verb@RETURN()@ statement, the function must call the
\verb@DETATCHSTATUSPTR()@ macro, with the status pointer of the
current function (not the subroutines) as its argument:
\begin{verbatim}
DETATCHSTATUSPTR( stat );
\end{verbatim}
This simply deallocates \verb@stat->statusPtr@ and sets it to
\verb@NULL@.  It is an error to exit the function with non-\verb@NULL@
\verb@statusPtr@, unless the exit was due to a subroutine failure.
\verb@ABORT()@ and \verb@ASSERT()@ check for this automatically; the
only place you normally need to call \verb@DETATCHSTATUSPTR()@ is
immediately before \verb@RETURN()@.  This macro also sets the status
code and the status message to nominal values.

Additionally, if a function successfully works around an error
reported by a subroutine, it should call \verb@DETATCHSTATUSPTR()@ and
\verb@ATTATCHSTATUSPTR()@ to create a fresh status pointer before
calling another subroutine.

\end{enumerate}

\subsection*{Reporting the current status}

The \verb@REPORTSTATUS()@ macro is used to issue a current status
report from the current function; the report is printed to
\verb@stderr@ or some other implementation-specific stream used for
error reporting.  \verb@REPORTSTATUS()@ takes the current status
pointer as its argument, and iteratively reports down any chain of
non-\verb@NULL@ \verb@statusPtr@ structures passed back by
subroutines.

\subsection*{Setting the top-level \texttt{Status} structure and
	\texttt{debuglevel}}

The top-level main function must set a global variable
\verb@debuglevel@, of standard C/C++ type \verb@int@ (not
\verb@INT4@).  It is invoked in this header as an \verb@extern@
variable, and controls the amount of error logging performed by the
various status macros.  The program should also declare an empty (all
fields set to zero) status structure to pass to its LAL functions.
The \verb@Status@ structure need only be declared and initialized
once, no matter how many LAL functions are called.  For example:
\begin{verbatim}
int debuglevel = 1;

int main( int argc, char **argv )
{
  static Status stat;
  MyFunction( &stat );
  REPORTSTATUS( &stat );
  return 0;
}
\end{verbatim}

A \verb@debuglevel@ of 0 means that no error logging will occur with
any of the status macros, with the exception of \verb@REPORTSTATUS()@.
A \verb@debuglevel@ of 1 means that a status message will be logged
(usually to \verb@stderr@) whenever a function terminates abnormally
(i.e.\ with non-zero \verb@statusCode@).  Any other value of
\verb@debuglevel@ means that a status message will be logged even when
functions exit nominally.  Programmers may also explicitly invoke
\verb@extern int debuglevel@ in their code to determine the verbosity
of their own error and warning messages.

Please note that all status macros with the exception of
\verb@REPORTSTATUS()@ can force a return from the calling routine.
This is a Bad Thing if the calling routine is \verb@main()@, since
\verb@main()@ must normally return \verb@int@ rather than \verb@void@.
It is therefore recommended that none of these macros other than
\verb@REPORTSTATUS()@ be used at the top level.

\subsection*{Non-confomant functions}

These standards apply only to functions that will be publicly
available in the LAL libraries.  Within a module, a programmer may
define and use subroutines that do not conform to the LAL function
standards, provided these routines are only visible within that
module.  Such functions should be declared as \verb@static@ to ensure
this.  A publicly-visible non-conformant function requires special
dispensation.

\subsection*{Notes}

Why are the status handling routines written as macros rather than
functions?  There are three good reasons.

First, many of the handling routines must be able to force an exit
from the function calling them.  This cannot be done if the routine is
in its own function, except by raising signal flags (which is a Bad
Thing according to LAL standards).

Second, it is useful for these routines to assign a status structure's
file and line fields using the \verb@__FILE__@ and \verb@__LINE__@
macros.  If the routine is its own function, then these will just give
the file and line number where the error handling routine is defined.
If the routine is a macro, then these will give the file and line
number where the macro was called, which is much more interesting.

Third, by expanding macros at compile time, the runtime performance of
the resulting code is marginally better.  Most of these macros will,
under nominal conditions, reduce to a single conditional test of an
integer value, with no additional overhead from function calling and
parameter passing.  Thus programmers can be encouraged to include
extensive error trapping in all their routines, without having to
worry about compromising performance.

\subsection*{Example: A LAL primer}

The following sections give a sample program program
\verb@LALPrimerTest.c@, along with its supporting header file
\verb@LALPrimer.h@ and function module \verb@LALPrimer.c@.  The
program itself is trivial to the point of silliness: it takes two
arguments from the command line and computes their ratio.  (Optionally
it can take a third command line argument to set the
\verb@debuglevel@.)  It is intended simply to illustrate how to use
the LAL status structure and macros in an actual, complete piece of
code.

For a more fully developed sample program, see the package
\verb@hello@.  That package also demonstrates how to document a
package using the autodocumentation utilities, which the
\verb@LALPrimer@ routines ignore.


\vfill{\footnotesize\input{LALStatusMacrosHV}}


\newpage\subsection{Sample header: \texttt{LALPrimer.h}}
\vspace{3ex}
\input{LALPrimerH}

\newpage\subsection{Sample module: \texttt{LALPrimer.c}}
\vspace{3ex}
\input{LALPrimerC}

\newpage\subsection{Sample program: \texttt{LALPrimerTest.c}}
\vspace{3ex}
\input{LALPrimerTestC}

</lalLaTeX> */


#ifndef _LALSTATUSMACROS_H
#define _LALSTATUSMACROS_H

#include "LALConfig.h"

#ifdef STDC_HEADERS
#include <string.h>
#else
#error "ERROR: non ansi standard headers"
#endif

#include "LALMalloc.h"
#include "LALDatatypes.h"
#include "LALError.h"
#include "LALRCSID.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (LALSTATUSMACROSH, "$Id$");

extern int debuglevel;

#define INITSTATUS( statusptr, funcname, id )                         \
do                                                                    \
{                                                                     \
  INT4 level;                                                         \
  if(!(statusptr))                                                    \
    {                                                                 \
      CHAR msg[1024];                                                 \
      sprintf(msg,"Abort: function %s, file %s, line %d, %s\n"        \
	      "       Null status pointer passed to function\n",      \
	      (funcname),__FILE__,__LINE__,(id));                     \
      LALAbort(msg);                                                  \
    }                                                                 \
  if((statusptr)->statusPtr) /* this may cause a memory leak */       \
    {                                                                 \
      if(!(statusptr)->level) (statusptr)->level = 1;                 \
      (statusptr)->Id        = (id);                                  \
      (statusptr)->function  = (funcname);                            \
      (statusptr)->statusPtr = NULL; /* warning: dislocated memory */ \
      ABORT(statusptr,-2,"INITSTATUS: non-null status pointer");      \
    }                                                                 \
  level = (statusptr)->level;                                         \
  memset((statusptr),0,sizeof(Status));                               \
  (statusptr)->level    = level > 0 ? level : 1 ;                     \
  (statusptr)->Id       = (id);                                       \
  (statusptr)->function = (funcname);                                 \
} while (0)

#define RETURN(statusptr)                                             \
do                                                                    \
{                                                                     \
  (statusptr)->file=__FILE__;                                         \
  (statusptr)->line=__LINE__;                                         \
  if(debuglevel==0 ||((debuglevel==1)&&((statusptr)->statusCode==0))) \
    {                                                                 \
      return;                                                         \
    }                                                                 \
  else if((statusptr)->statusCode==0)                                 \
    {                                                                 \
      LALPrintError("Nominal[%d]: ", (statusptr)->level);             \
      LALPrintError("function %s, file %s, line %d, %s\n",            \
	      (statusptr)->function,(statusptr)->file,                \
              (statusptr)->line,(statusptr)->Id);                     \
      return;                                                         \
    }                                                                 \
  else                                                                \
    {                                                                 \
      LALPrintError("Error[%d] %d: ", (statusptr)->level,             \
              (statusptr)->statusCode);                               \
      LALPrintError("function %s, file %s, line %d, %s\n",            \
              (statusptr)->function, (statusptr)->file,               \
              (statusptr)->line, (statusptr)->Id);                    \
      LALPrintError("          %s\n",(statusptr)->statusDescription); \
      return;                                                         \
    }                                                                 \
} while (0)

#define ABORT(statusptr,code,mesg)                                    \
do                                                                    \
{                                                                     \
  (statusptr)->file=__FILE__;                                         \
  (statusptr)->line=__LINE__;                                         \
  (statusptr)->statusCode=(code);                                     \
  (statusptr)->statusDescription=(mesg);                              \
  if((statusptr)->statusPtr)                                          \
    {                                                                 \
      LALFree((statusptr)->statusPtr);                                \
      (statusptr)->statusPtr=NULL;                                    \
    }                                                                 \
  if(debuglevel==0 || ((debuglevel==1)&&((code)==0)))                 \
    {                                                                 \
      return;                                                         \
    }                                                                 \
  else if((code)==0)                                                  \
    {                                                                 \
      LALPrintError("Nominal[%d]: ", (statusptr)->level);             \
      LALPrintError("function %s, file %s, line %d, %s\n",            \
	      (statusptr)->function,(statusptr)->file,                \
              (statusptr)->line,(statusptr)->Id);                     \
      return;                                                         \
    }                                                                 \
  else                                                                \
    {                                                                 \
      LALPrintError("Error[%d] %d: ", (statusptr)->level,(code));     \
      LALPrintError("function %s, file %s, line %d, %s\n",            \
              (statusptr)->function, (statusptr)->file,               \
              (statusptr)->line, (statusptr)->Id);                    \
      LALPrintError("         %s\n", (mesg));                         \
      return;                                                         \
    }                                                                 \
} while (0)

#define ASSERT(assertion,statusptr,code,mesg)                         \
do                                                                    \
{                                                                     \
  if(!(assertion))                                                    \
    {                                                                 \
      (statusptr)->file=__FILE__;                                     \
      (statusptr)->line=__LINE__;                                     \
      (statusptr)->statusCode=(code);                                 \
      (statusptr)->statusDescription=(mesg);                          \
      if((statusptr)->statusPtr)                                      \
	{                                                             \
	  LALFree((statusptr)->statusPtr);                            \
	  (statusptr)->statusPtr=NULL;                                \
	}                                                             \
      if(debuglevel==0 || ((debuglevel==1)&&((code)==0)))             \
	{                                                             \
	  return;                                                     \
	}                                                             \
      else if((code)==0)                                              \
	{                                                             \
          LALPrintError("Nominal[%d]: ", (statusptr)->level);         \
          LALPrintError("function %s, file %s, line %d, %s\n",        \
	          (statusptr)->function,(statusptr)->file,            \
                  (statusptr)->line,(statusptr)->Id);                 \
	  return;                                                     \
	}                                                             \
      else                                                            \
	{                                                             \
          LALPrintError("Error[%d] %d: ", (statusptr)->level,(code)); \
          LALPrintError("function %s, file %s, line %d, %s\n",        \
                  (statusptr)->function, (statusptr)->file,           \
                  (statusptr)->line, (statusptr)->Id);                \
          LALPrintError("         Assertion %s failed: %s\n",         \
                  #assertion, (mesg));                                \
	  return;                                                     \
	}                                                             \
    }                                                                 \
} while (0)

#define ATTATCHSTATUSPTR(statusptr)                                   \
do                                                                    \
{                                                                     \
  ASSERT(!(statusptr)->statusPtr,statusptr,-2,                        \
	 "ATTATCHSTATUSPTR: non-null status pointer");                \
  (statusptr)->statusPtr=(Status *)LALCalloc(1,sizeof(Status));       \
  ASSERT((statusptr)->statusPtr,statusptr,-4,                         \
	 "ATTATCHSTATUSPTR: memory allocation error");                \
  (statusptr)->statusPtr->level=(statusptr)->level + 1;               \
  (statusptr)->statusCode = -1;                                       \
  (statusptr)->statusDescription="Recursive error";                   \
} while (0)

#define DETATCHSTATUSPTR(statusptr)                                   \
do                                                                    \
{                                                                     \
  Status *ptr=(statusptr)->statusPtr;                                 \
  ASSERT(ptr,statusptr,-8,"DETATCHSTATUSPTR: null status pointer");   \
  while (ptr)                                                         \
  {                                                                   \
    Status *next=ptr->statusPtr;                                      \
    LALFree(ptr);                                                     \
    ptr=next;                                                         \
  }                                                                   \
  (statusptr)->statusPtr=NULL;                                        \
  (statusptr)->statusCode = 0;                                        \
  (statusptr)->statusDescription=NULL;                                \
} while (0)

#define TRY(func,statusptr)                                           \
do                                                                    \
{                                                                     \
  (func);                                                             \
  if((statusptr)->statusPtr->statusCode)                              \
    {                                                                 \
      (statusptr)->file=__FILE__;                                     \
      (statusptr)->line=__LINE__;                                     \
      (statusptr)->statusCode= -1;                                    \
      (statusptr)->statusDescription="Recursive error";               \
      if(debuglevel>0)                                                \
	{                                                             \
          LALPrintError("Error[%d] %d: ", (statusptr)->level, -1);    \
          LALPrintError("function %s, file %s, line %d, %s\n",        \
                  (statusptr)->function, (statusptr)->file,           \
                  (statusptr)->line, (statusptr)->Id);                \
          LALPrintError("          Function call %s failed\n",        \
                  #func);                                             \
	}                                                             \
      return;                                                         \
    }                                                                 \
} while (0)

#define CHECKSTATUSPTR(statusptr)                                     \
do                                                                    \
{                                                                     \
  if((statusptr)->statusPtr->statusCode)                              \
    {                                                                 \
      (statusptr)->file=__FILE__;                                     \
      (statusptr)->line=__LINE__;                                     \
      (statusptr)->statusCode= -1;                                    \
      (statusptr)->statusDescription="Recursive error";               \
      if(debuglevel>0)                                                \
	{                                                             \
          LALPrintError("Error[%d] %d: ", (statusptr)->level, -1);    \
          LALPrintError("function %s, file %s, line %d, %s\n",        \
                  (statusptr)->function, (statusptr)->file,           \
                  (statusptr)->line, (statusptr)->Id);                \
	  LALPrintError("          Function call failed\n");          \
	}                                                             \
      return;                                                         \
    }                                                                 \
} while (0)

#define REPORTSTATUS(statusptr)                                       \
do                                                                    \
{                                                                     \
  Status *ptr;                                                        \
  for(ptr=(statusptr);ptr;ptr=(ptr->statusPtr))                       \
    {                                                                 \
      LALPrintError("\nLevel %i: %s\n",ptr->level,ptr->Id);           \
      if (ptr->statusCode)                                            \
      {                                                               \
        LALPrintError("\tStatus code %i: %s\n",ptr->statusCode,       \
	              ptr->statusDescription);                        \
      }                                                               \
      else                                                            \
      {                                                               \
        LALPrintError("\tStatus code 0: Nominal\n");                  \
      }                                                               \
      LALPrintError("\tfunction %s, file %s, line %i\n",              \
                    ptr->function,ptr->file,ptr->line);               \
    }                                                                 \
} while (0)


#ifdef  __cplusplus
}
#endif

#endif /* _LALSTATUSMACROS_H */
