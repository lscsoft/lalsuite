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
#include <lal/LALStatusMacros.h>
\end{verbatim}

\noindent This header provides macros and functions for tracking and
reporting the runtime status of a program.  The intent is
simultaneously to standardize the error reporting, and to make the
reporting as transparent as possible to people coding individual
routines.

\subsection{Status-reporting objects}
\label{ss:status-reporting-objects}

LAL routines make use of two objects in reporting their current
status: the status structure \verb@LALStatus@, and the global integer
\verb@lalDebugLevel@.  These two objects are described in the following
sections.

\subsubsection{The \texttt{LALStatus} structure}
\idx[Type]{LALStatus}

LAL routines store their current execution status in a linked list of
structures of type \verb@LALStatus@, with each node in the list
representing a subroutine in the current calling sequence.  The
\verb@LALStatus@ structure is described in Sec.~\ref{ss:status-structure}
of the header \verb@LALDatatypes.h@, but for completeness, we explain
its fields below:
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
\item[\texttt{LALStatus *statusPtr}] A recursive pointer to another
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

\tt -16 & \vspace{-1.4ex}\tt INITSTATUS: non-zero xlalErrno &
The \verb@xlalError@ variable is non-zero, which suggests that an
error in an XLAL routine has occured and has not been handled. \\

\tt -16 & \vspace{-1.4ex}\tt RETURN: untrapped XLAL error code &
The \verb@xlalError@ variable is non-zero, which indicates that an
error in an XLAL routine has occured and has not been handled. \\

\hline
\end{tabular}
\end{center}

\subsubsection{The \texttt{lalDebugLevel}}
\idx[Variable]{lalDebugLevel}

The \verb@lalDebugLevel@ is a global variable, set at runtime, that
determines how much and what kind of debugging information will be
reported.  It is declared as an \verb@extern int@ in the header
\verb@LALStatusMacros.h@, and is therefore accessible in any standard
LAL module that includes this header.  Note, however, that it is
declared to be of the C type \verb@int@, which is usually but not
always a 32-bit integer (on some systems it may only be 16 bits).

The value of \verb@lalDebugLevel@ should be thought of not as a number,
but as a \emph{bit mask}, wherein each bit in the binary
representation turns on or off a specific type of status reporting.
At present, there are five types of status reporting, each associated
with a bit in \verb@lalDebugLevel@.

\paragraph{Error messages} tell the operator that a computation has
terminated abnormally, and has failed to produce an acceptable result.
Normally this is associated with assigning a non-zero
\verb@statusCode@; an error message is printed automatically whenever
a function exits with non-zero \verb@statusCode@.

\paragraph{Warning messages} tell the user that a computation is
working, but with unusual behaviour that might indicate an unreliable
or meaningless result.  Warnings do not normally result in a non-zero
\verb@statusCode@.

\paragraph{Information messages} tell the operator that the
computation is proceeding as expected, and simply provide additional
information about its progress.

\paragraph{Tracing messages} are printed automatically a subroutine
is called or returned; they simply track the current sequence of
function calls.

\paragraph{Memory information messages} are a special type of
information message; they tell the operator when and how much memory
is allocated or freed from the memory heap.

\paragraph{}The module \verb@LALError.c@ defines functions for
printing each of these types of status message.  Each type of message
is turned on by setting the corrsponding bit in \verb@lalDebugLevel@ to
1, and is suppressed by setting the bit to 0.  This header file
\verb@#define@s flags with numerical values designed to switch on the
appropriate bits.  Combinations of bits can be switched on by
combining these flags using the bitwise-\textit{or} operator,
\verb@|@.  The flags are defined as follows:

\begin{center}
\begin{tabular}{|lccl|}
\hline
Flag & Octal & Decimal & Meaning \\
\hline
\multicolumn{4}{|l|}{\it Primitive flags} \\
\tt LALNDEBUG   & 000000 &     0 & No debugging or status messages \\
\tt LALERROR    & 000001 &     1 & Turn on error messages \\
\tt LALWARNING  & 000002 &     2 & Turn on warning messages \\
\tt LALINFO     & 000004 &     4 & Turn on info messages \\
\tt LALTRACE    & 000010 &     8 & Turn on tracing messages \\
\tt LALMEMINFO  & 000020 &    16 & Turn on memory messages \\
\tt LALNMEMDBG  & 000040 &    32 & Turn off all memory debugging \\
\tt LALNMEMPAD  & 000100 &    64 & Turn off memory padding \\
\tt LALNMEMTRK  & 000200 &   128 & Turn off memory tracking \\
\tt LALMEMDBG   & 040000 & 16384 & Turn on memory debugging without messages \\
\multicolumn{4}{|l|}{\it Combination flags} \\
\tt LALMSGLVL1  & 000001 &     1 & Error messages only \\
\tt LALMSGLVL2  & 000003 &     3 & Error and warning messages \\
\tt LALMSGLVL3  & 000007 &     7 & Error, warning, and info messages \\
\tt LALMEMTRACE & 000030 &    24 & Memory and tracing messages \\
\tt LALALLDBG   & 077437 & 32543 & All messages and debugging \\
\hline
\end{tabular}
\end{center}
\idx[Constant]{LALNDEBUG}
\idx[Constant]{LALERROR}
\idx[Constant]{LALWARNING}
\idx[Constant]{LALINFO}
\idx[Constant]{LALTRACE}
\idx[Constant]{LALMEMINFO}
\idx[Constant]{LALMEMNDBG}
\idx[Constant]{LALMEMNPAD}
\idx[Constant]{LALMEMNTRK}
\idx[Constant]{LALMEMDBG}
\idx[Constant]{LALMSGLVL1}
\idx[Constant]{LALMSGLVL2}
\idx[Constant]{LALMSGLVL3}
\idx[Constant]{LALMEMTRACE}
\idx[Constant]{LALALLDBG}

The most significant bit
of \verb@lalDebugLevel@ has a special meaning in that it is not
associated with any type of status message.  However, certain pieces
of debugging or error-tracking code --- such as the memory leak
detection code in \verb@LALMalloc.c@ --- do not write status messages
and are not associated with a \verb@lalDebugLevel@ bit; instead, these
pieces of code are turned on for \emph{any} nonzero value of
\verb@lalDebugLevel@, unless the \verb@LALNMEMDBG@ bit is set.
Switching on only the most significant bit with
\verb@LALMEMDBG@ activates this code without turning on any other
error reporting.

To turn debugging code on or off at compile time (rather than
runtime), see Sec.~\ref{ss:compilation-flags}, below.

\subsection{Using the status tools}
\label{ss:using-status-tools}

The following summarizes everything the common programmer needs to
know in order to follow LAL standard error reporting.  It can be
treated as a primer on LAL coding conventions.

\subsubsection{LAL function calls}

All functions should have return type void.  The first argument of any
function should be a pointer to a structure of type \verb@LALStatus@.
Thus:
\begin{verbatim} 
void MyFunction( LALStatus *stat, ... )
\end{verbatim}
Since the function has no return code, it must report all errors or
failure through the status structure.  A function that is passed a
\verb@NULL@ pointer in place of the status pointer should terminate
the program with a \verb@SIGABRT@ signal, as this is its only way to
report the error.  However, this is one of the few circumstances under
which a function sould deliberately raise a signal.  In all other
cases the error should be trapped, reported in the status structure,
and control returned to the calling routine.

\subsubsection{Assigning an RCS \texttt{\$Id\$} string}
\idx{NRCSID()}

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

\subsubsection{Initializing the status structure}
\idx{INITSTATUS()}

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

\subsubsection{Normal return from a function}
\idx{RETURN()}

Upon completion, the function should issue the macro \verb@RETURN()@,
which takes one argument: the function's status pointer.
\begin{verbatim}
RETURN( stat );
\end{verbatim}
This takes the place of any return statements.  If
\verb@stat->statusCode@ is non-zero, the macro calls \verb@LALError()@
(see \verb@LALError.c@) to log \verb@stat->statusDescription@ and
other information, depending on implementation and the value of
\verb@lalDebugLevel@.  Typically \verb@RETURN()@ is used only for
successful completion, with other macros \verb@ABORT()@,
\verb@ASSERT()@, \verb@CHECKSTATUSPTR()@, and \verb@TRY()@ being used
to report failure.  However, it is possible for the programmer to
assign the fields of \verb@*stat@ by hand, and then issue
\verb@RETURN()@.

\subsubsection{Abnormal return from a function}
\idx{ABORT()}

The standard method to terminate a function unsuccessfully is with the
\verb@ABORT()@ macro, which takes three arguments: the status pointer,
the status code, and the status description string.  Normally the
various error codes and descriptions will be constants defined in the
function's header file \verb@MyHeader.h@:
\begin{verbatim}
ABORT( stat, MYHEADERH_EMYERR, MYHEADERH_MSGEMYERR );
\end{verbatim}
where the error code \verb@MYHEADERH_EMYERR@ and the error message
\verb@MYHEADERH_MSGEMYERR@ are defined in \verb@MyHeader.h@.  This
standard LAL naming convention for error messages prevents namespace
conflicts between different header files.  Like \verb@RETURN()@,
\verb@ABORT()@ correctly handles any status logging required by the
implementation and the \verb@lalDebugLevel@.  Note that \verb@ABORT()@
does \emph{not} raise a \verb@SIGABRT@ signal, but instead returns
control to the calling routine.

\subsubsection{Error checking within a function}
\idx{ASSERT()}

Another way to indicate an unsuccessful termination is with the macro
\verb@ASSERT()@, which takes as arguments a test statement, a status
pointer, a status code, and a status description.  The statement
\verb@ASSERT( assertion, ... );@ is in all ways equivalent to the
statement \verb@if ( !assertion ) ABORT( ... );@, except on a failure
the \verb@ASSERT()@ macro will also report the failed assertion.  In
the above example, one might have:
\begin{verbatim}
ASSERT( assertion, stat, MYHEADERH_EMYERR, MYHEADERH_MSGEMYERR );
\end{verbatim}

One subtle but important point is that the \verb@ASSERT()@ should be
used only to trap coding errors, rather than runtime errors, which
would be trapped using \verb@ABORT()@.  In other words, the assertion
should always test true in the final debugged program.  This is vital
because certain compilation flags will remove all \verb@ASSERT()@
macros at compile time, in order to speed execution of the final code.
See Sec.~\ref{ss:compilation-flags}, below.

Programmers should also be aware that using \verb@ASSERT()@ to exit a
function in normal runtime can have serious side effects.  For
example, it is an error to allocate dynamic memory to local variables
in a function and then fail to free it before returning.  Thus, if you
have dynamically allocated memory, you cannot then use \verb@ASSERT()@
for runtime error checking, as this does not permit you to free the
memory before returning.  Instead, you must explicitly check the
assertion, and, if it fails, free the memory and call \verb@ABORT()@.

\subsubsection{Calling subroutines}
\idx{ATTATCHSTATUSPTR()}
\idx{DETATCHSTATUSPTR()}
\idx{CHECKSTATUSPTR()}
\idx[Macro]{TRY()}

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
has already been allocated, \verb@ATTATCHSTATUSPTR()@ will raise a
\verb@SIGABRT@, as this is symptomatic of a coding error.

In most cases \verb@ATTATCHSTATUSPTR()@ need only be called once in a
given function, immediately after \verb@INITSTATUS()@, no matter how
many subroutine calls that function makes.  The exception is if the
function deals with (or ignores) errors reported by its subroutines.
In that case, the function should detatch the status pointer using
\verb@DETATCHSTATUSPTR()@ (below), and then re-attatch it.

The macro \verb@ATTATCHSTATUSPTR()@ sets the status code to be $-1$
and the status message to be \verb@"Recursive error"@.  These flags
are unset when \verb@DETATCHSTATUSPTR()@ (below) is called.  This is
so that a use of \verb@RETURN()@ prior to detatching the status
pointer will yield an error.

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
This simply deallocates \verb@stat->statusPtr@ (and any subsequent
structures in the list), and sets it to \verb@NULL@.  It is an error
to exit the function with non-\verb@NULL@ \verb@statusPtr@, unless the
exit was due to a subroutine failure.  \verb@ABORT()@ and
\verb@ASSERT()@ check for this automatically; the only place you
normally need to call \verb@DETATCHSTATUSPTR()@ is immediately before
\verb@RETURN()@.  This macro also sets the status code and the status
message to nominal values.

Additionally, if a function successfully works around an error
reported by a subroutine, it should call \verb@DETATCHSTATUSPTR()@ and
\verb@ATTATCHSTATUSPTR()@ to create a fresh status pointer before
calling another subroutine.

\end{enumerate}

\subsubsection{Cleaning up after subroutine failure}
\idx[Macro]{BEGINFAIL()}
\idx[Macro]{ENDFAIL()}

Although they are convenient, the \verb@TRY()@ and
\verb@CHECKSTATUSPTR()@ macros have a serious drawback in that they
may cause the calling function to return immediately.  If the calling
function had previously allocated any local memory storage, this
memory will be cast adrift, with no means of accessing or subsequently
freeing it (short of terminating the runtime process).  Such a memory
leak is a violation of the LAL function standard.

The macros \verb@BEGINFAIL()@ and \verb@ENDFAIL()@ allow a function to
test the return code of a subroutine, and, if that indicates a
failure, to execute one or more ``cleanup'' instructions before itself
returning.  Each macro takes a single argument: the current function's
status pointer.  The macros must occur in matched pairs, and use the
same syntax as a \verb@do ... while@ statement: they either span a
single instruction, or a block of instructions enclosed in braces.

For example, if a function had allocated memory to some pointer
\verb@localPointer@, any subsequent call to a subroutine
\verb@LALSubroutine()@ would take the following form:
\begin{verbatim}
LALSubroutine( stat->statusPtr, ... );
BEGINFAIL( stat )
  LALFree( localPointer );
ENDFAIL( stat );
\end{verbatim}
For another example, if a function had to create three vectors
\verb@*vector1@, \verb@*vector2@, \verb@*vector3@, the allocation
would look something like this:
\begin{verbatim}
TRY( LALSCreateVector( stat->statusPtr, &vector1, 100 ), stat );

LALSCreateVector( stat->statusPtr, &vector2, 100 );
BEGINFAIL( stat )
  TRY( LALSDestroyVector( stat->statusPtr, &vector1 ), stat );
ENDFAIL( stat );

LALSCreateVector( stat->statusPtr, &vector3, 100 );
BEGINFAIL( stat ) {
  TRY( LALSDestroyVector( stat->statusPtr, &vector1 ), stat );
  TRY( LALSDestroyVector( stat->statusPtr, &vector2 ), stat );
} ENDFAIL( stat );
\end{verbatim}
As indicated above, the cleanup instructions can include calls to
other LAL routines.  The \verb@BEGINFAIL( stat )@ macro call first
checks \verb@stat->statusPtr@ to see if a subroutine error has
occured.  If it has, the macro detaches and saves that pointer, then
attaches a new \verb@stat->statusPtr@ to be used in calls to the
cleanup routines.  After the cleanup instructions have been executed,
the \verb@ENDFAIL( stat )@ macro call reattaches the saved status
pointer and returns with a subroutine error code.  In this way, the
returned status list indicates where the original failure occurred,
rather than giving an uninformative report from the last cleanup
routine.

Of course a \emph{second} failure in one of the cleanup routines can
cause serious problems.  If the routine was called using a
\verb@TRY()@ macro, it will force an immediate return from the calling
function, with a status code and status list indicating how the cleanp
routine failed.  The original status list saved by \verb@BEGINFAIL()@
is lost.  While this loss does constitute a memory leak, the failure
of a cleanup routine in itself indicates that there are serious
problems with the memory management.

It is possible to nest \verb@BEGINFAIL()@\ldots\verb@ENDFAIL();@
blocks, but this is unlikely to serve any useful purpose.  Once
cleanup routines start to fail, it is probably beyond the scope of the
LAL function to deal with the resulting memory leaks.

\subsubsection{Issuing status messages}
\idx{LALError()}
\idx{LALWarning()}
\idx{LALInfo()}
\idx{LALTrace()}
\idx{REPORTSTATUS()}

The module \verb@LALError.c@ defines the functions \verb@LALError()@,
\verb@LALWarning()@, \verb@LALInfo()@, and \verb@LALTrace()@ to issue
various types of status message.  This is the preferred means of
printing status messages, since each type of message can be activated
or suppressed by setting \verb@lalDebugLevel@ appropriately.  In fact,
\verb@LALError()@ and \verb@LALTrace()@ are called automatically by
the status macros whenever they are required, so most LAL modules will
explicitly invoke only the \verb@LALWarning()@ and \verb@LALInfo()@
functions.

\verb@LALStatusMacros.h@ provides a macro, \verb@REPORTSTATUS()@,
which is used to report the current state of the \verb@LALStatus@ list.
It takes a status pointer as its argument:
\begin{verbatim}
REPORTSTATUS( stat );
\end{verbatim}
This macro iteratively prints the contents of \verb@stat@ and all
subsequent structures in the list to the error log.

The action of \verb@REPORTSTATUS()@ is not suppressed by any value of
\verb@lalDebugLevel@.  Therefore, as a rule, it should only be called by
test programs, not by LAL routines intended for use in production
code.

\subsubsection{Setting the initial \texttt{LALStatus} structure and
global \texttt{lalDebugLevel}}
\idx[Type]{LALStatus}
\idx[Variable]{lalDebugLevel}

As mentioned above, any module including \verb@LALStatusMacros.h@
includes the global variable \verb@lalDebugLevel@ as an
\verb@extern int@.  At least one module in the final executable
program must have a global \emph{declaration} of \verb@int lalDebugLevel@
(not \verb@extern int@), and assign \verb@lalDebugLevel@ a value.  In
most cases \verb@lalDebugLevel@ will be declared in the module containing
the \verb@main()@ function, and will be assigned a value on
declaration or from command-line arguments to \verb@main()@.
Alternatively, if the LAL functions are to be embedded in a non-LAL
program, \verb@lalDebugLevel@ can be declared and set in the topmost
module that calls LAL functions.

A \verb@LALStatus@ structure should also be declared as a local variable
in the \verb@main()@ function of a LAL program, or in the topmost
function calling LAL functions withing a non-LAL program, to pass in
its LAL function calls.  The structure must be empty (all fields set
to zero) before being passed into a function.  The \verb@LALStatus@
structure need only be declared and initialized once, no matter how
many LAL functions are called.

Thus a typical LAL program might look something like the following:

\begin{verbatim}
int lalDebugLevel = 1;

int main( int argc, char **argv )
{
  static LALStatus stat;
  MyFunction( &stat );
  REPORTSTATUS( &stat );
  return stat.statusCode;
}
\end{verbatim}

Please note that all status macros described above can force a return
from the calling routine.  This is a Bad Thing if the calling routine
is \verb@main()@, since \verb@main()@ must normally return \verb@int@
rather than \verb@void@.  It is therefore recommended that none of
these macros except \verb@REPORTSTATUS()@ be used at the top level.

\subsubsection{Non-confomant functions}

These standards apply only to functions that will be publicly
available in the LAL libraries.  Within a module, a programmer may
define and use subroutines that do not conform to the LAL function
standards, provided these routines are only visible within that
module.  Such functions should be declared as \verb@static@ to ensure
this.  A publicly-visible non-conformant function requires special
dispensation.

\subsection{Compilation flags}
\label{ss:compilation-flags}

LAL provides two flags that can be used to exclude or modify debugging
code at compile time.  Although these flags are typically
\verb@#define@d or \verb@#undef@ined globally and can affect many
modules (notably modules in the \verb@support@ package), their primary
effect is on the debugging and status-reporting tools defined in this
header.  The two flags are named \verb@NDEBUG@ and \verb@NOLALMACROS@.

\subsubsection{The \texttt{NDEBUG} flag}

Setting the \verb@NDEBUG@ (or \verb@LAL_NDEBUG@) flag turns off debugging and
error-reporting code, in order to get condensed production-line programs.  As
far as error reporting is concerned, setting the \verb@NDEBUG@ flag at compile
time is similar to setting \verb@lalDebugLevel@ equal to zero at runtime, in
that it suppresses all status messages and memory leak detection.  However,
the \verb@NDEBUG@ flag accoplishes this by telling the compiler preprocessor
to remove the relevant code from the object file, thus eliminating frequent
and unnecessary tests on \verb@lalDebugLevel@.  When debugging is turned off,
the global integer variable \verb@lalNoDebug@ is non-zero; otherwise it is
zero.

Compiling with the \verb@NDEBUG@ flag set also removes all
\verb@ASSERT()@ macros from the object code, in keeping with the
philosophy that \verb@ASSERT()@ statements should only be used to
catch coding bugs, not runtime errors.

\subsubsection{The \texttt{NOLALMACROS} flag}

Setting the \verb@NOLALMACROS@ flag replaces the status-handling
macros described above with actual functions that accomplish the same
results.  These functions are defined in the module \verb@LALError.c@.
Function calls introduce computational and memory overheads that are
absent in macros, since macro replacement occurs at compile time.
However, there are circumstances in which one might want to use
function calls rather than macro replacement.

For example, debuggers typically cannot step through the individual
instructions within a macro.  If a conflict somehow arose between a
particular piece of code and its status macros, this conflict would be
easier to catch and resolve by replacing the macros with function
calls into which the debugger could step.

\subsubsection{Using the compilation flags}

There are three ways to set these flags when compiling LAL programs or
libraries.

When compiling your own modules, the flags can be set using one or
more \verb@#define@ statements within the module or its header file:
\begin{verbatim}
#define NDEBUG
#define NOLALMACROS
\end{verbatim}
To restrict the scope of these flags, they should later be unset using
the corresponding \verb@#undef@ statements.

Alternatively, these can be set in the \verb@Makefile@ or when
compiling.  The syntax for most UNIX C compilers is something like the
following:
\begin{verbatim}
> gcc ... -DNDEBUG -DNOLALMACROS ...
\end{verbatim}

If you want to compile a large number of modules, or the entire
library, under the effects of one or more of these flags, you will not
want to go through and modify every header or \verb@Makefile@.
Instead, you may add either \verb@-DNDEBUG@ or \verb@-DNOLALMACROS@
(or both) to the environment variable \verb@CPPFLAGS@.  They will then
automatically be set for all compilations done in that environment.
The command for doing this in \verb@sh@ or \verb@bash@ shells is:
\begin{verbatim}
> CPPFLAGS="$CPPFLAGS -DNDEBUG -DNOLALMACROS"
\end{verbatim}
while in \verb@csh@ or \verb@tcsh@ shells it is:
\begin{verbatim}
> setenv CPPFLAGS "$CPPFLAGS -DNDEBUG -DNOLALMACROS"
\end{verbatim}
Note that if you plan to do further LAL code development on the same
system, you may want to keep two versions of the library around: one
with the flag(s) set and one without.

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

It should be mentioned that, for these reasons, compiling a module
with the \verb@NOLALMACROS@ flag above does not actually eliminate the
status handling macros.  Rather, the macros are modified to call
specialized functions that do most (but not all) of the processing.

\subsection*{Example: A LAL primer}

The following sections give a sample program program
\verb@LALPrimerTest.c@, along with its supporting header file
\verb@LALPrimer.h@ and function module \verb@LALPrimer.c@.  The
program itself is trivial to the point of silliness: it takes two
arguments from the command line and computes their ratio.  (Optionally
it can take a third command line argument to set the
\verb@lalDebugLevel@.)  It is intended simply to illustrate how to use
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

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (LALSTATUSMACROSH, "$Id$");

extern int lalDebugLevel;
extern const int lalNoDebug;

#define LAL_EXLAL     16384
#define LAL_MSGEXLAL  "Failure in an XLAL routine"
#define ABORTXLAL(sp) ABORT(sp,LAL_EXLAL,LAL_MSGEXLAL)

#ifndef NOLALMACROS

#define INITSTATUS( statusptr, funcname, id )                                 \
  if ( (statusptr) )                                                          \
  {                                                                           \
    INT4 level_ = (statusptr)->level ;                                        \
    INT4 statp_ = (statusptr)->statusPtr ? 1 : 0 ;                            \
    memset( (statusptr), 0, sizeof( LALStatus ) ); /* possible memory leak */ \
    (statusptr)->level    = level_ > 0 ? level_ : 1 ;                         \
    (statusptr)->Id       = (id);                                             \
    (statusptr)->function = (funcname);                                       \
    SETSTATUSFILELINE( statusptr );                                           \
    (void) LALTrace( statusptr, 0 );                                          \
    if ( statp_ )                                                              \
    {                                                                         \
      ABORT( statusptr, -2, "INITSTATUS: non-null status pointer" );          \
    }                                                                         \
    else if ( xlalErrno )                                                     \
    {                                                                         \
      ABORT( statusptr, -16, "INITSTATUS: non-zero xlalErrno" );              \
    }                                                                         \
  }                                                                           \
  else                                                                        \
    lalAbortHook( "Abort: function %s, file %s, line %d, %s\n"                \
                  "       Null status pointer passed to function\n",          \
                  (funcname), __FILE__, __LINE__, (id) )

#define RETURN( statusptr )                                                   \
  if ( 1 )                                                                    \
  {                                                                           \
    SETSTATUSFILELINE( statusptr );                                           \
    if ( (statusptr)->statusCode )                                            \
      (void) LALError( statusptr, "RETURN:" );                                \
    (void) LALTrace( statusptr, 1 );                                          \
    if ( xlalErrno )                                                          \
    {                                                                         \
      ABORT( statusptr, -32, "RETURN: untrapped XLAL error" );                \
    }                                                                         \
    return;                                                                   \
  }                                                                           \
  else (void)(0)

#define ATTATCHSTATUSPTR(statusptr)                                           \
  if ( !(statusptr)->statusPtr )                                              \
  {                                                                           \
    (statusptr)->statusPtr = (LALStatus *)LALCalloc( 1, sizeof( LALStatus ) );\
    if ( !(statusptr)->statusPtr )                                            \
    {                                                                         \
      ABORT( statusptr, -4, "ATTATCHSTATUSPTR: memory allocation error" );    \
    }                                                                         \
    (statusptr)->statusPtr->level = (statusptr)->level + 1;                   \
  }                                                                           \
  else                                                                        \
    ABORT( statusptr, -2, "ATTATCHSTATUSPTR: non-null status pointer" )

#define DETATCHSTATUSPTR( statusptr )                                         \
  if ( (statusptr)->statusPtr )                                               \
  {                                                                           \
    FREESTATUSPTR( statusptr );                                               \
    (statusptr)->statusCode = 0;                                              \
    (statusptr)->statusDescription = NULL;                                    \
  }                                                                           \
  else                                                                        \
    ABORT( statusptr, -8, "DETATCHSTATUSPTR: null status pointer" )

#define ABORT( statusptr, code, mesg )                                        \
  if ( 1 )                                                                    \
  {                                                                           \
    if ( statusptr->statusPtr ) FREESTATUSPTR( statusptr );                   \
    SETSTATUS( statusptr, code, mesg );                                       \
    if ( code )                                                               \
      (void) LALError( statusptr, "ABORT:" );                                 \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
  else (void)(0)

#ifdef LAL_NDEBUG
#define ASSERT( assertion, statusptr, code, mesg )
#else
#define ASSERT( assertion, statusptr, code, mesg )                            \
  if ( !(assertion) )                                                         \
  {                                                                           \
    if ( statusptr->statusPtr )                                               \
      FREESTATUSPTR( statusptr );                                             \
    SETSTATUS( statusptr, code, mesg );                                       \
    (void) LALError( statusptr, "Assertion \"" #assertion "\" failed:" );     \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
  else (void)(0)
#endif

#define TRY( func, statusptr )                                                \
  if ( (func), (statusptr)->statusPtr->statusCode )                           \
  {                                                                           \
    SETSTATUS( statusptr, -1, "Recursive error" );                            \
    (void) LALError( statusptr, "Function call \"" #func "\" failed:" );      \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
  else (void)(0)

#define CHECKSTATUSPTR( statusptr )                                           \
  if ( (statusptr)->statusPtr->statusCode )                                   \
  {                                                                           \
    SETSTATUS( statusptr, -1, "Recursive error" );                            \
    (void) LALError( statusptr, "CHECKSTATUSPTR:" );                          \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
  else (void)(0)

#define FREESTATUSPTR( statusptr )                                            \
  do                                                                          \
  {                                                                           \
    LALStatus *next_ = (statusptr)->statusPtr->statusPtr;                      \
    LALFree( (statusptr)->statusPtr );                                        \
    (statusptr)->statusPtr = next_;                                            \
  }                                                                           \
  while ( (statusptr)->statusPtr )

#define REPORTSTATUS( statusptr )                                             \
  do                                                                          \
  {                                                                           \
    LALStatus *ptr_;                                                          \
    for ( ptr_ = (statusptr); ptr_; ptr_ = ptr_->statusPtr )                  \
    {                                                                         \
      LALPrintError( "\nLevel %i: %s\n", ptr_->level, ptr_->Id );             \
      if ( ptr_->statusCode )                                                 \
      {                                                                       \
        LALPrintError( "\tStatus code %i: %s\n", ptr_->statusCode,            \
                       ptr_->statusDescription );                             \
      }                                                                       \
      else                                                                    \
      {                                                                       \
        LALPrintError( "\tStatus code 0: Nominal\n" );                        \
      }                                                                       \
      LALPrintError( "\tfunction %s, file %s, line %i\n",                     \
                     ptr_->function, ptr_->file, ptr_->line );                \
    }                                                                         \
  } while ( 0 )

#else /* NOLALMACROS */

#define INITSTATUS( statusptr, funcname, id ) \
  if ( LALInitStatus( statusptr, funcname, id, __FILE__, __LINE__ ) ) return

#define RETURN( statusptr ) \
  if ( LALPrepareReturn( statusptr, __FILE__, __LINE__ ), 1 ) return

#define ATTATCHSTATUSPTR( statusptr ) \
  if ( LALAttatchStatusPtr( statusptr, __FILE__, __LINE__ ) ) return

#define DETATCHSTATUSPTR( statusptr ) \
  if ( LALDetatchStatusPtr( statusptr, __FILE__, __LINE__ ) ) return

#define ABORT( statusptr, code, mesg ) \
  if ( LALPrepareAbort( statusptr, code, mesg, __FILE__, __LINE__ ), 1 ) return

#ifdef LAL_NDEBUG
#define ASSERT( assertion, statusptr, code, mesg )
#else
#define ASSERT( assertion, statusptr, code, mesg )                            \
  if ( !(assertion) )                                                         \
  {                                                                           \
    LALPrepareAssertFail( statusptr, code, mesg,                              \
                          "Assertion \"" #assertion "\" failed:",             \
                          __FILE__, __LINE__ );                               \
    return;                                                                   \
  }                                                                           \
  else (void)(0)
#endif

#define TRY( func, statusptr )                                                \
  do                                                                          \
  {                                                                           \
    (func);                                                                   \
    if ( LALCheckStatusPtr( statusptr, "Function call \"" #func "\" failed:", \
                            __FILE__, __LINE__ ) )                            \
      return;                                                                 \
  }                                                                           \
  while ( 0 )

#define CHECKSTATUSPTR( statusptr )                                           \
  if ( LALCheckStatusPtr( statusptr, "CHECKSTATUSPTR:", __FILE__, __LINE__ ) )\
    return
  
#endif /* NOLALMACROS */

/* these just have to be macros... */

#define BEGINFAIL( statusptr )                                                \
do {                                                                          \
  if ( !(statusptr) ) {                                                       \
    ABORT( statusptr, -8, "BEGINFAIL: null status pointer" );                 \
  }                                                                           \
  if ( !( (statusptr)->statusPtr ) ) {                                        \
    ABORT( statusptr, -8, "BEGINFAIL: null status pointer pointer" );         \
  }                                                                           \
  if ( (statusptr)->statusPtr->statusCode ) {                                 \
    LALStatus *statusPtrSave_ = (statusptr)->statusPtr;                       \
    (statusptr)->statusPtr = NULL;                                            \
    ATTATCHSTATUSPTR( statusptr );                                            \
    do

#define ENDFAIL( statusptr )                                                  \
    while ( 0 );                                                              \
    DETATCHSTATUSPTR( statusptr );                                            \
    (statusptr)->statusPtr = statusPtrSave_;                                  \
    SETSTATUS( statusptr, -1, "Recursive error" );                            \
    (void) LALError( statusptr, "ENDFAIL:" );                                 \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
} while ( 0 )

#define SETSTATUSFILELINE( statusptr ) \
  ( ( void ) ( (statusptr)->file = __FILE__, (statusptr)->line = __LINE__ ) )

#define SETSTATUS( statusptr, code, mesg )                                    \
  ( SETSTATUSFILELINE( statusptr ),                                           \
    (statusptr)->statusDescription = (mesg),                                  \
    (statusptr)->statusCode = (code) )


#ifdef  __cplusplus
}
#endif

#endif /* _LALSTATUSMACROS_H */
