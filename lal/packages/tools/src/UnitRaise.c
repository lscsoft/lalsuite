/************************************ <lalVerbatim file="UnitRaiseCV">
Author: Whelan, J. T.
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UnitRaise.c}}
\label{tools:ss:UnitRaise.c}

Raises an \texttt{LALUnit} structure to a specified rational power.

\subsubsection*{Prototypes}
\input{UnitRaiseCP}
\idx{LALUnitRaise()}

\subsubsection*{Description}

This function raises the \texttt{LALUnit} structure \texttt{*input} to
the rational power \texttt{*power}.  In this way, units such as
s$^{1/2}$ and m$^{-1}$ can be created using existing units.

\subsubsection*{Algorithm}

The function first multiplies the overall power of ten
\texttt{input->powerOfTen} by the rational number \texttt{*power},
checking to make sure that the resulting power is still an integer.
It then multiplies each of the rational powers in \texttt{*input} by
\texttt{*power} by na\"{\i}ve multiplication of rational numbers
$$
\left(\frac{N_1}{1+D_1}\right)\left( \frac{N_2}{1+D_2} \right)
= \frac{N_1 N_2}{1 + (1+D_1)(1+D_2)-1} 
$$ 
and then calls \texttt{LALUnitNormalize()} to bring the result into
standard form.

\subsubsection*{Uses}

\texttt{LALUnitNormalize()}

\subsubsection*{Notes}

\vfill{\footnotesize\input{UnitRaiseCV}}

******************************************************* </lalLaTeX> */ 
#define TRUE 1
#define FALSE 0

#include <lal/LALStdlib.h>
#include <lal/Units.h>

NRCSID( UNITRAISEC, "$Id$" );

/* <lalVerbatim file="UnitRaiseCP"> */
void 
LALUnitRaise (LALStatus *status, LALUnit *output, const LALUnit *input, const RAT4 *power)
/* </lalVerbatim> */
     /* Raise a Unit variable to a rational power */
{
  LALUnit     unReduced;
  UINT2       i;
  INT4        numer;
  UINT4       denom, denom1, denom2;

  INITSTATUS( status, "LALUnitRaise", UNITRAISEC );
  ATTATCHSTATUSPTR (status);


  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( power != NULL, status, UNITSH_ENULLPPARAM, UNITSH_MSGENULLPPARAM );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  denom2 = power->denominatorMinusOne + 1;

  ASSERT( input->powerOfTen % denom2 == 0, status, 
	  UNITSH_ENONINT, UNITSH_MSGENONINT);

  numer = (input->powerOfTen / (INT4) denom2) * power->numerator;

  ASSERT(numer < 32767L && numer > -32768L, status, UNITSH_EOVERFLOW,
	 UNITSH_MSGEOVERFLOW);

  unReduced.powerOfTen = numer;

  for (i=0; i<LALNumUnits; ++i) {
    denom1 = 1 + input->unitDenominatorMinusOne[i];
    denom = denom1 * denom2;

    ASSERT(denom - 1 < 65535L, status, UNITSH_EOVERFLOW,
	   UNITSH_MSGEOVERFLOW);

    unReduced.unitDenominatorMinusOne[i] = denom - 1;

    numer = input->unitNumerator[i] * power->numerator;

    ASSERT(numer < 32767L && numer > -32768L, status, UNITSH_EOVERFLOW,
	   UNITSH_MSGEOVERFLOW);

    unReduced.unitNumerator[i] = numer;
  } /* for i */

  LALUnitNormalize(status->statusPtr, output, &unReduced);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
