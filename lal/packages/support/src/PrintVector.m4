dnl $Id$
/************************************ <lalVerbatim file="PrintVectorCV">
Author: Allen, B.; generalized by J.T. Whelan
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{PrintVector.c}}
\label{ss:PrintVector.c}

Print a $\langle\mbox{datatype}\rangle$Vector object into a file.  For
use in non-production and test code only.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{PrintVectorCP}
\idx{LALZPrintVector()}
\idx{LALCPrintVector()}
\idx{LALDPrintVector()}
\idx{LALSPrintVector()}
\idx{LALI2PrintVector()}
\idx{LALI4PrintVector()}
\idx{LALI8PrintVector()}
\idx{LALU2PrintVector()}
\idx{LALU4PrintVector()}
\idx{LALU8PrintVector()}
\idx{LALCHARPrintVector()}
\idx{LALPrintVector()}

\subsubsection*{Description}

Each member of this family of functions prints the elements of
$\langle\mbox{datatype}\rangle$\verb+Vector+ into a file.  Note: the
file names are $\langle\mbox{DT}\rangle$\verb+PrintVector.000+,
$\langle\mbox{DT}\rangle$\verb+PrintVector.001+, and so on.
($\langle\mbox{DT}\rangle$ is the abbreviation for the datatype,
included in the function names above.) The file numbers are
incremented with each additional call.  This function is for debugging
use only: it uses a static internal variable to keep track of the file
number so it should not be used in any real analysis codes.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
LALFopen()
LALFclose()
\end{verbatim}

\subsubsection*{Notes}

This function uses an internal static variable to keep track of file
numbers.  For this reason it should only be used for debugging
purposes in test functions, not in any production code.

Additionally, since printf cannot handle INT8 as integers, the
functions \verb&LALI8PrintVector()& and \verb&LALU8PrintVector()& use
a typecast to REAL8 and are thus only valid for numbers between around
$-10^{15}$ and $10^{15}$.

The output format is two or three space-separated columns: the first
column is the index of the element; the second is the element itself
for real and integer vectors and the real part of the element for
complex vectors; complex vectors have a third column containing the
imaginary part of the element.

\vfill{\footnotesize\input{PrintVectorCV}}

</lalLaTeX> */


#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/PrintVector.h>

/* <lalVerbatim file="VectorFactoriesNRCSID"> */
NRCSID( PRINTVECTORC, "$Id$" );
/* </lalVerbatim> */

define(`TYPECODE',`Z')
include(`LALPrintVector.m4')

define(`TYPECODE',`C')
include(`LALPrintVector.m4')

define(`TYPECODE',`D')
include(`LALPrintVector.m4')

define(`TYPECODE',`S')
include(`LALPrintVector.m4')

define(`TYPECODE',`I2')
include(`LALPrintVector.m4')

define(`TYPECODE',`I4')
include(`LALPrintVector.m4')

define(`TYPECODE',`I8')
include(`LALPrintVector.m4')
/* Note that LALI8PrintVector does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */

define(`TYPECODE',`U2')
include(`LALPrintVector.m4')

define(`TYPECODE',`U4')
include(`LALPrintVector.m4')

define(`TYPECODE',`U8')
include(`LALPrintVector.m4')
/* Note that LALU8PrintVector does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */

define(`TYPECODE',`CHAR')
include(`LALPrintVector.m4')

define(`TYPECODE',`')
include(`LALPrintVector.m4')
