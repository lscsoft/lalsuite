dnl $Id$
/************************************ <lalVerbatim file="ReadFrequencySeriesCV">
Author: Torres, C. W.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{ReadFrequencySeries.c}}
\label{ss:ReadFrequencySeries.c}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ReadFrequencySeriesCP}
\index{\texttt{LALZReadFrequencySeries()}}
\index{\texttt{LALCReadFrequencySeries()}}
\index{\texttt{LALDReadFrequencySeries()}}
\index{\texttt{LALSReadFrequencySeries()}}
\index{\texttt{LALReadFrequencySeries()}}

\subsubsection*{Description}

Each member of this family of functions reads from a file the output
of the corresponding \texttt{PrintFrequencySeries} routine.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
LALOpenDataFile()
LALParseUnitString()
LALCHARCreateVector()
LALCHARDestroyVector()
LALDCreateVector()
LALDDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

These functions perform I/O operations, which are not a part of LAL
proper. They should only be used for debugging purposes in test
functions, not in any production code.

\vfill{\footnotesize\input{ReadFrequencySeriesCV}}

</lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <string.h>
#include <lal/LALDatatypes.h>
#include <lal/ReadFTSeries.h>

void LALCHARCreateVector( LALStatus *, CHARVector **, UINT4 );
void LALCHARDestroyVector( LALStatus *, CHARVector ** );
void LALDCreateVector( LALStatus *, REAL8Vector **, UINT4 );
void LALDDestroyVector( LALStatus *, REAL8Vector ** );
void LALParseUnitString ( LALStatus *status,
	       	          LALUnit *output,
		          const CHARVector *input );
static const LALUnit lalDimensionlessUnit 
	= {  0, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };

/* <lalVerbatim file="ReadFrequencySeriesNRCSID"> */
NRCSID( READFREQUENCYSERIESC, "$Id$" );
/* </lalVerbatim> */

/* Change the first instance of the target to '\0'; returns 0 on success,
   1 on failure */
static INT2 changeCharToNull (CHAR *string, CHAR target, CHAR *offEndOfString)
{
  CHAR *charPtr;

  for ( charPtr=string; charPtr<offEndOfString; ++charPtr )
  {
    if ( *charPtr == target )
    {
      *charPtr = '\0';
      return 0;
    }
  }
  return 1;
}


define(`TYPECODE',`Z')
include(`LALReadFrequencySeries.m4')

define(`TYPECODE',`C')
include(`LALReadFrequencySeries.m4')

define(`TYPECODE',`D')
include(`LALReadFrequencySeries.m4')

define(`TYPECODE',`S')
include(`LALReadFrequencySeries.m4')

define(`TYPECODE',`')
include(`LALReadFrequencySeries.m4')
