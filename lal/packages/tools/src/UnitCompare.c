/************************************ <lalVerbatim file="UnitCompareCV">
Author: J. T. Whelan <whelan@oates.utb.edu>
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UnitCompare.c}}
\label{tools:ss:UnitCompare.c}

Compares two \texttt{LALUnit} structures, returning true if they are
equivalent, false otherwise.

\subsubsection*{Prototypes}
\input{UnitCompareCP}
\idx{LALUnitCompare()}

\subsubsection*{Description}

This function determines whether the units represented by
\texttt{input->unitOne} and \texttt{input->unitTwo} are the same (both
dimensionally and in the power-of-ten prefactor).  In this way,
programs and programmers can verify that quantities have the expected
units.
 
\subsubsection*{Algorithm}

The function first uses \texttt{LALUnitNormalize()} to bring both unit
structures into standard form, then compares the powers of ten and the
numerator and denominator of each exponent of a fundamental unit in
turn.

\subsubsection*{Uses}

\texttt{LALUnitNormalize()}

\subsubsection*{Notes}

\vfill{\footnotesize\input{UnitCompareCV}}

******************************************************* </lalLaTeX> */ 
#define TRUE 1
#define FALSE 0

#include <lal/LALStdlib.h>
#include <lal/Units.h>

NRCSID( UNITCOMPAREC, "$Id$" );

/* <lalVerbatim file="UnitCompareCP"> */
void
LALUnitCompare (LALStatus *status, BOOLEAN *output, const LALUnitPair *input)
/* </lalVerbatim> */
     /* Compare two units structures, returning false if they differ */
{
  INT2        i;
  LALUnit     unitOne, unitTwo;

  INITSTATUS( status, "LALUnitCompare", UNITCOMPAREC );
  ATTATCHSTATUSPTR (status);


  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  LALUnitNormalize(status->statusPtr, &unitOne, &(input->unitOne));
  LALUnitNormalize(status->statusPtr, &unitTwo, &(input->unitTwo));

  if (unitOne.powerOfTen != unitTwo.powerOfTen)
  {
    *output = FALSE;
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }
  for (i=0; i<LALNumUnits; ++i) {
    if (
	unitOne.unitNumerator[i] != unitTwo.unitNumerator[i]
	|| unitOne.unitDenominatorMinusOne[i] 
	   != unitTwo.unitDenominatorMinusOne[i]
	)
    {
      *output = FALSE;
      DETATCHSTATUSPTR(status);
      RETURN(status);
    } /* if */
  } /* for i */
  *output = TRUE;
  DETATCHSTATUSPTR(status);
  RETURN(status);  
}
