/********************************* <lalVerbatim file="LALRCSIDHV">
Author: Unknown?  Provided by Stuart Anderson
$Id$
********************************** </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALRCSID.h}}
\label{s:LALRCSID.h}

Provides macros for assigning an RCS ID string to a file.

\subsection*{Synopsis}
\begin{verbatim}
#include "LALRCSID.h"
\end{verbatim}

\noindent This header defines a pair of macros \verb@RCSID()@ and
\verb@NRCSID()@, which are used to assign Revision Control System
(RCS) ID strings to a given source file.  Whenever a file is added or
changed in the central LAL distribution, the RCS software searches for
symbols of the form \texttt{\$Id\$} or \texttt{\$Id:} \ldots\
\texttt{\$}, and expands them into strings that give information about
the file's name, version number, when it was last modified and by
whom.  For an example of an RCS ID string, look at the \texttt{\$Id:}
\ldots\ \texttt{\$} line at the bottom of this page.

The macro \verb@RCSID()@ is called as follows:

\vspace{2ex}
\noindent\texttt{RCSID(\$Id\$);}
\vspace{2ex}

\noindent This assigns the RCS ID string to a variable
\verb@static const char *rcsid@ in a given module or header.  This
variable will be loaded onto the stack whenever code from that module
or header is used.  This can be used as a diagnostic tool by
debuggers.

The macro \verb@NRCSID()@ is called in the following manner:

\vspace{2ex}
\noindent\texttt{NRCSID(MYFILEC,\$Id\$);}
\vspace{2ex}

\noindent This assigns the RCS ID string to the variable
\verb@MYFILEC@, as above.  Standard LAL naming conventions are that
the variable name should be the file name converted to capital
letters, \emph{with} file extensions but \emph{without} periods.  Thus
the module \verb@MyFile.c@ should store its ID in the variable
\verb@MYFILEC@, while the header \verb@MyHeader.h@ should store it in
\verb@MYHEADERH@.

LAL convention dictates that all modules and header files must call
\verb@NRCSID()@ using the naming convention above.  The call must be
made after all subsidiary header files have been included (notably
this header file, or \verb@LALStdlib.h@ which includes this header
file), but before any actual functions or function prototypes.  In
modules containing LAL functions, the RCS ID string will typically be
assigned to the \verb@Status@ structure for those functions; see the
documentation for the header \verb@LALStatusMacros.h@.

This header is included automatically by all standard LAL headers,
including \verb@LALDatatypes.h@, \verb@LALStatusMacros.h@, and
\verb@LALStdlib.h@.  Thus if you have included \emph{any} of the
standard LAL headers, you will have gotten \verb@LALRCSID.h@ as well,
and don't need to \verb@#include@ it separately.  However, including
it separately is not an error, as this and all LAL headers are
required to have double-include protection.

</lalLaTeX> */

#ifndef _LALRCSID_H
#define _LALRCSID_H

#ifdef  __cplusplus
extern "C" {
#endif

#if !defined(lint)
#  ifndef __GNUC__
#    define RCSID(id)       static volatile const char *rcsid = (id)
#    define NRCSID(name,id) static volatile const char *name  = (id)
#  else
#    define RCSID(id) \
       static volatile const char __attribute__ ((unused)) *rcsid = (id)
#    define NRCSID(name,id) \
       static volatile const char __attribute__ ((unused)) *name  = (id)
#  endif /* !__GNUC__ */
#else
#  define RCSID(id)	   typedef void useless
#  define NRCSID(name,id)  typedef void useless
#endif /* !lint */

NRCSID (LALRCSIDH, "$Id$");


/* <lalLaTeX>
\vfill{\footnotesize\input{LALRCSIDHV}}
</lalLaTeX> */


#ifdef  __cplusplus
}
#endif

#endif /* _LALRCSID_H */
