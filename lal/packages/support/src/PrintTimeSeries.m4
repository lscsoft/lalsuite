dnl $Id$
/************************************ <lalVerbatim file="PrintTimeSeriesCV">
Author: Whelan, J. T.
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

\begin{verbatim}
LALCHARCreateVector()
LALCHARDestroyVector()
LALUnitAsString()
\end{verbatim}

\subsubsection*{Notes}

This function's arguments do not conform to the LAL spec.  For this
reason it should only be used for debugging purposes in test
functions, not in any production code.

Additionally, since printf cannot handle INT8 as integers, the
functions \verb&LALI8PrintTimeSeries()& and
\verb&LALU8PrintTimeSeries()& use a typecast to REAL8 and are thus
only valid for numbers between around $-10^{15}$ and $10^{15}$.

The first five lines of the file are a header containing:
\begin{enumerate}
\item the name of the series
\item the starting epoch 
\item the units expressed in terms of the basic SI units
\item column labels
\end{enumerate}
after which come the data, one per line.


The output format is two or three tab-separated columns: the first
column is the time corresponding to the row in question, in seconds
after the series' starting epoch; for real and integer time
series, the second column is the value of the series; for complex time
series, the second column is the real part and the third the imaginary
part of the value.

\vfill{\footnotesize\input{PrintTimeSeriesCV}}

</lalLaTeX> */


#include <lal/LALStdlib.h>
#include <stdio.h>
#include <lal/LALDatatypes.h>
#include <lal/PrintFTSeries.h>

void LALCHARCreateVector( LALStatus *, CHARVector **, UINT4 );
void LALCHARDestroyVector( LALStatus *, CHARVector ** );
void LALUnitAsString( LALStatus *status, CHARVector *output,
                      const LALUnit *input );

enum { LALUnitTextSize = sizeof("10^-32768 m^-32768/32767 kg^-32768/32767 "
				"s^-32768/32767 A^-32768/32767 " 
				"K^-32768/32767 strain^-32768/32767 "
				"count^-32768/32767") };

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
