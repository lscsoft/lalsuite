dnl $Id$
/************************************ <lalVerbatim file="PrintTimeSeriesCV">
Author: J.T. Whelan
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{PrintTimeSeries.c}}
\label{ss:PrintTimeSeries.c}

Print a $\langle\mbox{datatype}\rangle$TimeSeries object into a
file.  For use in non-production and test code only.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{PrintTimeSeriesCP}
\index{\verb&LALZPrintTimeSeries()&}
\index{\verb&LALCPrintTimeSeries()&}
\index{\verb&LALDPrintTimeSeries()&}
\index{\verb&LALSPrintTimeSeries()&}
\index{\verb&LALI2PrintTimeSeries()&}
\index{\verb&LALI4PrintTimeSeries()&}
\index{\verb&LALI8PrintTimeSeries()&}
\index{\verb&LALU2PrintTimeSeries()&}
\index{\verb&LALU4PrintTimeSeries()&}
\index{\verb&LALU8PrintTimeSeries()&}
\index{\verb&LALPrintTimeSeries()&}

\subsubsection*{Description}

Each member of this family of functions prints the elements of
$\langle\mbox{datatype}\rangle$\verb+TimeSeries+ into a file.  Note:
the file name is specified using a character string.  This function is
for debugging use only: its arguments do not conform to LAL standards
so it should not be used in any real analysis codes.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

This function's arguments do not conform to the LAL spec.  For this
reason it should only be used for debugging purposes in test
functions, not in any production code.

Additionally, since printf cannot handle INT8 as integers, the
functions \verb&LALI8PrintTimeSeries()& and
\verb&LALU8PrintTimeSeries()& use a typecast to REAL8 and are thus
only valid for numbers between around $-10^{15}$ and $10^{15}$.

The output format is three or four space-separated columns: the first
column is the time of the beginning of the time series, in seconds
after the GPS reference epoch (1980 January 6), and is the same for
each entry; the second column is the time into the series (so that the
sum of the first two columns is the time of the entry in question);
for real and integer time series, the third column is the value of the
series; for complex time series, the third column is the real part and
the fourth the imaginary part of the value.

\vfill{\footnotesize\input{PrintTimeSeriesCV}}

</lalLaTeX> */


#include "LALStdlib.h"
#include <stdio.h>
#include "LALDatatypes.h"
#include "PrintFTSeries.h"

/* <lalVerbatim file="PrintTimeSeriesNRCSID"> */
NRCSID( PRINTTIMESERIESC, "$Id$" );
/* </lalVerbatim> */

define(`TYPECODE',`Z')
include(`LALPrintTimeSeries.m4')

define(`TYPECODE',`C')
include(`LALPrintTimeSeries.m4')

define(`TYPECODE',`D')
include(`LALPrintTimeSeries.m4')

define(`TYPECODE',`S')
include(`LALPrintTimeSeries.m4')

define(`TYPECODE',`I2')
include(`LALPrintTimeSeries.m4')

define(`TYPECODE',`I4')
include(`LALPrintTimeSeries.m4')

define(`TYPECODE',`I8')
include(`LALPrintTimeSeries.m4')
/* Note that LALI8PrintTimeSeries does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */

define(`TYPECODE',`U2')
include(`LALPrintTimeSeries.m4')

define(`TYPECODE',`U4')
include(`LALPrintTimeSeries.m4')

define(`TYPECODE',`U8')
include(`LALPrintTimeSeries.m4')
/* Note that LALU8PrintTimeSeries does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */

define(`TYPECODE',`')
include(`LALPrintTimeSeries.m4')
