/************************************ <lalVerbatim file="UnitRaiseCV">
Author: J. T. Whelan <jtwhelan@loyno.edu>
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

LALUnit * XLALUnitRaiseRAT4( LALUnit *output, const LALUnit *input,
    const RAT4 *power )
{
  static const char *func = "XLALUnitRaiseRAT4";
  LALUnit     unReduced;
  UINT2       i;
  INT4        numer;
  UINT4       denom, denom1, denom2;
  
  if ( ! output || ! input || ! power )
    XLAL_ERROR_NULL( func, XLAL_EFAULT );

  denom2 = power->denominatorMinusOne + 1;

  if ( input->powerOfTen % denom2 )
    XLAL_ERROR_NULL( func, XLAL_EINVAL );

  numer = (input->powerOfTen / (INT4) denom2) * power->numerator;

  if ( numer >= 32767L || numer <= -32768L )
    XLAL_ERROR_NULL( func, XLAL_ERANGE );

  unReduced.powerOfTen = numer;

  for (i=0; i<LALNumUnits; ++i) {
    denom1 = 1 + input->unitDenominatorMinusOne[i];
    denom = denom1 * denom2;

    if ( denom - 1 >= 65535L )
      XLAL_ERROR_NULL( func, XLAL_ERANGE );

    unReduced.unitDenominatorMinusOne[i] = denom - 1;

    numer = input->unitNumerator[i] * power->numerator;

    if ( numer >= 32767L || numer <= -32768L )
      XLAL_ERROR_NULL( func, XLAL_ERANGE );

    unReduced.unitNumerator[i] = numer;
  } /* for i */

  *output = unReduced;
  if ( XLALUnitNormalize( output ) == XLAL_FAILURE )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  return output;
}


LALUnit * XLALUnitRaiseINT2( LALUnit *output, const LALUnit *input,
    INT2 power )
{
  static const char *func = "XLALUnitRaiseINT2";
  RAT4 pow;
  pow.numerator = power;
  pow.denominatorMinusOne = 0;
  if ( ! XLALUnitRaiseRAT4( output, input, &pow ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return output;
}


LALUnit * XLALUnitSquare( LALUnit *output, const LALUnit *input )
{
  static const char *func = "XLALUnitRaiseSquare";
  RAT4 pow;
  pow.numerator = 2;
  pow.denominatorMinusOne = 0;
  if ( ! XLALUnitRaiseRAT4( output, input, &pow ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return output;
}


LALUnit * XLALUnitSqrt( LALUnit *output, const LALUnit *input )
{
  static const char *func = "XLALUnitRaiseSqrt";
  RAT4 pow;
  pow.numerator = 1;
  pow.denominatorMinusOne = 1;
  if ( ! XLALUnitRaiseRAT4( output, input, &pow ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return output;
}


LALUnit * XLALUnitInvert( LALUnit *output, const LALUnit *input )
{
  static const char *func = "XLALUnitInvert";
  RAT4 pow;
  pow.numerator = -1;
  pow.denominatorMinusOne = 0;
  if ( ! XLALUnitRaiseRAT4( output, input, &pow ) )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return output;
}


/* <lalVerbatim file="UnitRaiseCP"> */
void 
LALUnitRaise (LALStatus *status, LALUnit *output, const LALUnit *input, const RAT4 *power)
/* </lalVerbatim> */
     /* Raise a Unit variable to a rational power */
{
  UINT4       denom2;

  INITSTATUS( status, "LALUnitRaise", UNITRAISEC );

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( power != NULL, status, UNITSH_ENULLPPARAM, UNITSH_MSGENULLPPARAM );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  denom2 = power->denominatorMinusOne + 1;

  ASSERT( input->powerOfTen % denom2 == 0, status, 
	  UNITSH_ENONINT, UNITSH_MSGENONINT);


  if ( ! XLALUnitRaiseRAT4( output, input, power ) )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT( status, UNITSH_ENONINT, UNITSH_MSGENONINT);
      case XLAL_ERANGE:
        ABORT( status, UNITSH_EOVERFLOW, UNITSH_MSGEOVERFLOW);
      default:
        ABORTXLAL( status );
    }
  }

  RETURN(status);
}
