/************************************ <lalVerbatim file="UnitMultiplyCV">
Author: Whelan, J. T.
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UnitMultiply.c}}
\label{ss:UnitMultiply.c}

Multiplies two \texttt{LALUnit} structures.

\subsubsection*{Prototypes}
\input{UnitMultiplyCP}
\index{\texttt{LALUnitMultiply()}}

This function multiplies together the \texttt{LALUnit} structures
\texttt{input->unitOne} and \texttt{input->unitTwo}, thus allowing a
module to \textit{e.g.}, multiply two \texttt{REAL8TimeSeries} and
give the resulting \texttt{REAL8TimeSeries} the correct units.

\subsubsection*{Algorithm}

The function first adds together the overall powers of ten in the two
input unit structures, then adds each of the corresponding rational
powers in \texttt{input->unitOne} and \texttt{input->unitTwo} by na\"{\i}ve
addition of rational numbers
$$
\frac{N_1}{1+D_1} + \frac{N_2}{1+D_2} = 
\frac{N_1 (1+D_2) +  N_2(1+D_1)}{1 + (1+D_1)(1+D_2)-1}
$$
and then calls \texttt{LALUnitNormalize()} to bring the result into
standard form.

\subsubsection*{Uses}

\texttt{LALUnitNormalize()}

\subsubsection*{Notes}

\vfill{\footnotesize\input{UnitMultiplyCV}}

******************************************************* </lalLaTeX> */ 
#define TRUE 1
#define FALSE 0

#include <lal/LALStdlib.h>
#include <lal/Units.h>

NRCSID( UNITMULTIPLYC, "$Id$" );

/* <lalVerbatim file="UnitMultiplyCP"> */
void 
LALUnitMultiply (LALStatus *status, LALUnit *output, const LALUnitPair *input)
/* </lalVerbatim> */
{
  LALUnit     unReduced;
  UINT2        i;
  INT4         numer;
  UINT4        denom, denom1, denom2;

  INITSTATUS( status, "LALUnitMultiply", UNITMULTIPLYC );
  ATTATCHSTATUSPTR (status);

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  numer = input->unitOne.powerOfTen + input->unitTwo.powerOfTen;

  ASSERT(numer < 32767L && numer > -32768L, status, UNITSH_EOVERFLOW,
	 UNITSH_MSGEOVERFLOW);

  unReduced.powerOfTen = numer;
  for (i=0; i<LALNumUnits; ++i) {
    denom1 = 1 + input->unitOne.unitDenominatorMinusOne[i];
    denom2 = 1 + input->unitTwo.unitDenominatorMinusOne[i];
    denom = denom1 * denom2;

    ASSERT(denom - 1 < 65535L, status, UNITSH_EOVERFLOW,
	   UNITSH_MSGEOVERFLOW);

    /* One could use the gcd function to find the common factors of
       denom1 and denom2, but we have to reduce the fractions after
       addition anyway; consider e.g., 1/6 + 1/10 = (5+3)/30 = 4/15 */
    unReduced.unitDenominatorMinusOne[i] = denom - 1;
    numer = ((INT4) denom2) * input->unitOne.unitNumerator[i]
      + ((INT4) denom1) * input->unitTwo.unitNumerator[i];

    ASSERT(numer < 32767L && numer > -32768L, status, UNITSH_EOVERFLOW,
	   UNITSH_MSGEOVERFLOW);

    unReduced.unitNumerator[i] = numer;
  } /* for i */
  LALUnitNormalize(status->statusPtr, output, &unReduced);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
