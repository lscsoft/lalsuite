/*
<lalLaTeX>
\clearpage
\section{Module \texttt{Sequence.c}}
\label{s:Sequence.c}

Author: \verb|Kipp Cannon <kipp@gravity.phys.uwm.edu>|

\noindent
Revision: \verb|$Id$|

\subsection{NAME}

\texttt{XLALCreate}\textit{sequencetype}\texttt{()},
\texttt{LALCreate}\textit{sequencetype}\texttt{()},
\texttt{XLALDestroy}\textit{sequencetyp}\texttt{()},
\texttt{LALDestroy}\textit{sequencetype}\texttt{()},
\texttt{XLALCut}\textit{sequencetype}\texttt{()},
\texttt{LALCut}\textit{sequencetype}\texttt{()},
\texttt{XLALShrink}\textit{sequencetype}\texttt{()},
\texttt{LALShrink}\textit{sequencetype}\texttt{()} --- A suite of sequence
manipulation functions.

\subsection{SYNOPSIS}

\begin{verbatim}
#include <lal/Sequence.h>

\end{verbatim}
</lalLaTeX>
 */

#include <string.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALErrno.h>
#include <lal/LALError.h>
#include <lal/Sequence.h>

NRCSID(SEQUENCEC, "$Id$");

define(`DATATYPE',REAL4)
include(SequenceC.m4)

define(`DATATYPE',REAL8)
include(SequenceC.m4)

define(`DATATYPE',COMPLEX8)
include(SequenceC.m4)

define(`DATATYPE',COMPLEX16)
include(SequenceC.m4)

define(`DATATYPE',INT2)
include(SequenceC.m4)

define(`DATATYPE',UINT2)
include(SequenceC.m4)

define(`DATATYPE',INT4)
include(SequenceC.m4)

define(`DATATYPE',UINT4)
include(SequenceC.m4)

define(`DATATYPE',INT8)
include(SequenceC.m4)

define(`DATATYPE',UINT8)
include(SequenceC.m4)

/*
<lalLaTeX>
\subsection{DESCRIPTION}

This is a suite of functions for creating, destroying, and manipulating LAL
sequences.  The following methods are implemented:
\begin{itemize}
\item \texttt{XLALCreate}\textit{sequencetype}\texttt{()}
\item \texttt{XLALDestroy}\textit{sequencetype}\texttt{()}
\item \texttt{XLALCut}\textit{sequencetype}\texttt{()}
\item \texttt{XLALShrink}\textit{sequencetype}\texttt{()}
\end{itemize}
One set of these methods is available for each defined sequence type, for
example \texttt{XLALCreateREAL4Sequence()} is available for creating
sequences of \texttt{REAL4} data.  Additionally, LAL-style wrappers are
provided for each XLAL function.  For example,
\texttt{LALCreateREAL4Sequence()} is provided which is equivalent to the
XLAL version in all respects except that it adheres to the LAL calling
conventions (i.e.\ it takes a \texttt{LALStatus} pointer as its first
argument, its return type is \texttt{void}, etc.).

The behaviour of each method is probably self-evident:  the create methods
allocate memory for a new sequence and initialize its contents, the destroy
methods free all memory associated with a sequence.  The cut functions
create a new sequence by extracting a length of data from an existing
sequence of the same type.  The shrink functions perform the same operation
as the cut functions but do not create a new sequence, instead modifying
the sequence in place.
</lalLaTeX>
 */
