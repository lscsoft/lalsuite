/*
<lalLaTeX>
\clearpage
\section{Module \texttt{TimeSeries.c}}

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
\texttt{LALShrink}\textit{seriestype}\texttt{()} --- A suite of time series
manipulation functions.

\subsection{SYNOPSIS}

\begin{verbatim}
#include <lal/TimeSeries.h>

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
#include <lal/TimeSeries.h>

NRCSID(TIMESERIESC, "$Id$");

define(`DATATYPE',REAL4)
include(TimeSeriesC.m4)

define(`DATATYPE',REAL8)
include(TimeSeriesC.m4)

define(`DATATYPE',COMPLEX8)
include(TimeSeriesC.m4)

define(`DATATYPE',COMPLEX16)
include(TimeSeriesC.m4)

define(`DATATYPE',INT2)
include(TimeSeriesC.m4)

define(`DATATYPE',UINT2)
include(TimeSeriesC.m4)

define(`DATATYPE',INT4)
include(TimeSeriesC.m4)

define(`DATATYPE',UINT4)
include(TimeSeriesC.m4)

define(`DATATYPE',INT8)
include(TimeSeriesC.m4)

define(`DATATYPE',UINT8)
include(TimeSeriesC.m4)

/*
<lalLaTeX>
\subsection{DESCRIPTION}

This is a suite of functions for creating, destroying, and manipulating LAL
time series.  The following methods are implemented:
\begin{itemize}
\item \texttt{XLALCreate}\textit{seriestype}\texttt{()}
\item \texttt{XLALDestroy}\textit{seriestype}\texttt{()}
\item \texttt{XLALCut}\textit{seriestype}\texttt{()}
\item \texttt{XLALShrink}\textit{seriestype}\texttt{()}
\end{itemize}
One set of these methods is available for each defined time series type,
for example \texttt{XLALCreateREAL4TimeSeries()} is available for creating
time series of \texttt{REAL4} data.  Additionally, LAL-style wrappers are
provided for each XLAL function.  For example,
\texttt{LALCreateREAL4TimeSeries()} is provided which is equivalent to the
XLAL version in all respects except that it adheres to the LAL calling
conventions (i.e.\ it takes a \texttt{LALStatus} pointer as its first
argument, its return type is \texttt{void}, etc.).

The behaviour of each method is probably self-evident:  the create methods
allocate memory for a new time series and initialize its contents, the
destroy methods free all memory associated with a time series.  The cut
functions create a new time series by extracting a length of data from an
existing time series of the same type and updating all meta data
appropriately.  The shrink functions perform the same operation as the cut
functions but do not create a new time series, instead modifying the time
series in place.
</lalLaTeX>
 */
