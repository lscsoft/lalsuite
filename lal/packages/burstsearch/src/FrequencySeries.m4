/*
<lalLaTeX>
\clearpage
\section{Module \texttt{FrequencySeries.c}}

Author: \verb|Kipp Cannon <kipp@gravity.phys.uwm.edu>|

\noindent
Revision: \verb|$Id$|

\subsection{NAME}

\texttt{XLALCreate}\textit{seriestype}\texttt{()},
\texttt{LALCreate}\textit{seriestype}\texttt{()},
\texttt{XLALDestroy}\textit{seriestyp}\texttt{()},
\texttt{LALDestroy}\textit{seriestype}\texttt{()},
\texttt{XLALCut}\textit{seriestype}\texttt{()},
\texttt{LALCut}\textit{seriestype}\texttt{()},
\texttt{XLALShrink}\textit{seriestype}\texttt{()},
\texttt{LALShrink}\textit{seriestype}\texttt{()} --- A suite of frequency
series manipulation functions.

\subsection{SYNOPSIS}

\begin{verbatim}
#include <lal/FrequencySeries.h>

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
#include <lal/FrequencySeries.h>

NRCSID(FREQUENCYSERIESC, "$Id$");

define(`DATATYPE',REAL4)
include(FrequencySeriesC.m4)

define(`DATATYPE',REAL8)
include(FrequencySeriesC.m4)

define(`DATATYPE',COMPLEX8)
include(FrequencySeriesC.m4)

define(`DATATYPE',COMPLEX16)
include(FrequencySeriesC.m4)

define(`DATATYPE',INT2)
include(FrequencySeriesC.m4)

define(`DATATYPE',UINT2)
include(FrequencySeriesC.m4)

define(`DATATYPE',INT4)
include(FrequencySeriesC.m4)

define(`DATATYPE',UINT4)
include(FrequencySeriesC.m4)

define(`DATATYPE',INT8)
include(FrequencySeriesC.m4)

define(`DATATYPE',UINT8)
include(FrequencySeriesC.m4)

/*
<lalLaTeX>
\subsection{DESCRIPTION}

This is a suite of functions for creating, destroying, and manipulating LAL
frequency series.  The following methods are implemented:
\begin{itemize}
\item \texttt{XLALCreate}\textit{seriestype}\texttt{()}
\item \texttt{XLALDestroy}\textit{seriestype}\texttt{()}
\item \texttt{XLALCut}\textit{seriestype}\texttt{()}
\item \texttt{XLALShrink}\textit{seriestype}\texttt{()}
\end{itemize}
One set of these methods is available for each defined frequency series
type, for example \texttt{XLALCreateREAL4TimeSeries()} is available for
creating frequency series of \texttt{REAL4} data.  Additionally, LAL-style
wrappers are provided for each XLAL function.  For example,
\texttt{LALCreateREAL4TimeSeries()} is provided which is equivalent to the
XLAL version in all respects except that it adheres to the LAL calling
conventions (i.e.\ it takes a \texttt{LALStatus} pointer as its first
argument, its return type is \texttt{void}, etc.).

The behaviour of each method is probably self-evident:  the create methods
allocate memory for a new frequency series and initialize its contents, the
destroy methods free all memory associated with a frequency series.  The
cut functions create a new frequency series by extracting a length of data
from an existing frequency series of the same type and updating all meta
data appropriately.  The shrink functions perform the same operation as the
cut functions but do not create a new frequency series, instead modifying
the frequency series in place.
</lalLaTeX>
 */
