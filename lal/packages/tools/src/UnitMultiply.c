/************************************ <lalVerbatim file="UnitMultiplyCV">
Author: J. T. Whelan <jtwhelan@loyno.edu>
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UnitMultiply.c}}
\label{tools:ss:UnitMultiply.c}

Multiplies two \texttt{LALUnit} structures.

\subsubsection*{Prototypes}
\input{UnitMultiplyCP}
\idx{LALUnitMultiply()}

This function multiplies together the \texttt{LALUnit} structures
\texttt{*(input->unitOne)} and \texttt{*(input->unitTwo)}, thus allowing a
module to \textit{e.g.}, multiply two \texttt{REAL8TimeSeries} and
give the resulting \texttt{REAL8TimeSeries} the correct units.

\subsubsection*{Algorithm}

The function first adds together the overall powers of ten in the two
input unit structures, then adds each of the corresponding rational
powers in \texttt{*(input->unitOne)} and \texttt{*(input->unitTwo)} by na\"{\i}ve
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

LALUnit * XLALUnitMultiply( LALUnit *output, const LALUnit *unit1, const LALUnit *unit2 )
{
  static const char *func = "XLALUnitMultiply";
  LALUnit     unReduced;
  UINT2        i;
  INT4         numer;
  UINT4        denom, denom1, denom2;

  if ( ! output || ! unit1 || ! unit2 )
    XLAL_ERROR_NULL( func, XLAL_EFAULT );

  numer = unit1->powerOfTen + unit2->powerOfTen;
  if ( numer >= 32767L || numer <= -32768L )
    XLAL_ERROR_NULL( func, XLAL_ERANGE );

  unReduced.powerOfTen = numer;
  for (i=0; i<LALNumUnits; ++i) {
    denom1 = 1 + unit1->unitDenominatorMinusOne[i];
    denom2 = 1 + unit2->unitDenominatorMinusOne[i];
    denom = denom1 * denom2;

    if ( denom >= 65535L )
      XLAL_ERROR_NULL( func, XLAL_ERANGE );

    /* One could use the gcd function to find the common factors of
       denom1 and denom2, but we have to reduce the fractions after
       addition anyway; consider e.g., 1/6 + 1/10 = (5+3)/30 = 4/15 */
    unReduced.unitDenominatorMinusOne[i] = denom - 1;
    numer = ((INT4) denom2) * unit1->unitNumerator[i]
      + ((INT4) denom1) * unit2->unitNumerator[i];

    if ( numer >= 32767L || numer <= -32768L )
      XLAL_ERROR_NULL( func, XLAL_ERANGE );

    unReduced.unitNumerator[i] = numer;
  } /* for i */

  *output = unReduced;
  if ( XLALUnitNormalize( output ) == XLAL_FAILURE )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  return output;
}


LALUnit * XLALUnitDivide( LALUnit *output, const LALUnit *unit1, const LALUnit *unit2 )
{
  static const char *func = "XLALUnitDivide";
  /* invert unit2 and then multiply by unit1 ... use output as tmp var */
  if ( ! XLALUnitInvert( output, unit2 ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  if ( ! XLALUnitMultiply( output, unit1, output ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return output;
}


/* <lalVerbatim file="UnitMultiplyCP"> */
void 
LALUnitMultiply (LALStatus *status, LALUnit *output, const LALUnitPair *input)
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALUnitMultiply", UNITMULTIPLYC );

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  
  if ( ! XLALUnitMultiply( output, input->unitOne, input->unitTwo ) )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ERANGE:
        ABORT(status, UNITSH_EOVERFLOW, UNITSH_MSGEOVERFLOW);
      default:
        ABORTXLAL(status);
    }
  }
  RETURN(status);
}
