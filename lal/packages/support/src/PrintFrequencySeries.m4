dnl $Id$
/************************************ <lalVerbatim file="PrintFrequencySeriesCV">
Author: J.T. Whelan
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{PrintFrequencySeries.c}}
\label{ss:PrintFrequencySeries.c}

Print a $\langle\mbox{datatype}\rangle$FrequencySeries object into a
file.  For use in non-production and test code only.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{PrintFrequencySeriesCP}
\index{\texttt{LALZPrintFrequencySeries()}}
\index{\texttt{LALCPrintFrequencySeries()}}
\index{\texttt{LALDPrintFrequencySeries()}}
\index{\texttt{LALSPrintFrequencySeries()}}
\index{\texttt{LALI2PrintFrequencySeries()}}
\index{\texttt{LALI4PrintFrequencySeries()}}
\index{\texttt{LALI8PrintFrequencySeries()}}
\index{\texttt{LALU2PrintFrequencySeries()}}
\index{\texttt{LALU4PrintFrequencySeries()}}
\index{\texttt{LALU8PrintFrequencySeries()}}
\index{\texttt{LALPrintFrequencySeries()}}

\subsubsection*{Description}

Each member of this family of functions prints the elements of
$\langle\mbox{datatype}\rangle$\verb+FrequencySeries+ into a file.
Note: the file name is specified using a character string.  This
function is for debugging use only: its arguments do not conform to
LAL standards so it should not be used in any real analysis codes.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}


This function's arguments do not conform to the LAL spec.  For this
reason it should only be used for debugging purposes in test
functions, not in any production code.

Additionally, since printf cannot handle INT8 as integers, the
functions \verb&LALI8PrintFrequencySeries()& and
\verb&LALU8PrintFrequencySeries()& use a typecast to REAL8 and are
thus only valid for numbers between around $-10^{15}$ and $10^{15}$.

The output format is two or four space-separated columns: the first
column is the frequency in hertz corresponfing to the row in question;
for real and integer frequency series, the second column is the value
of the series; for complex time series, the second column is the real
part and the third the imaginary part of the value.

\vfill{\footnotesize\input{PrintFrequencySeriesCV}}

</lalLaTeX> */


#include "LALStdlib.h"
#include <stdio.h>
#include "LALDatatypes.h"
#include "PrintFTSeries.h"

/* <lalVerbatim file="PrintFrequencySeriesNRCSID"> */
NRCSID( PRINTFREQUENCYSERIESC, "$Id$" );
/* </lalVerbatim> */

define(`TYPECODE',`Z')
include(`LALPrintFrequencySeries.m4')

define(`TYPECODE',`C')
include(`LALPrintFrequencySeries.m4')

define(`TYPECODE',`D')
include(`LALPrintFrequencySeries.m4')

define(`TYPECODE',`S')
include(`LALPrintFrequencySeries.m4')

define(`TYPECODE',`I2')
include(`LALPrintFrequencySeries.m4')

define(`TYPECODE',`I4')
include(`LALPrintFrequencySeries.m4')

define(`TYPECODE',`I8')
include(`LALPrintFrequencySeries.m4')
/* Note that LALI8PrintFrequencySeries does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */

define(`TYPECODE',`U2')
include(`LALPrintFrequencySeries.m4')

define(`TYPECODE',`U4')
include(`LALPrintFrequencySeries.m4')

define(`TYPECODE',`U8')
include(`LALPrintFrequencySeries.m4')
/* Note that LALU8PrintFrequencySeries does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */

define(`TYPECODE',`')
include(`LALPrintFrequencySeries.m4')
