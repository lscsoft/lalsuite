/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#define TRUE 1
#define FALSE 0

#include <lal/LALStdlib.h>
#include <lal/Units.h>

/**
 * \addtogroup UnitRaise_c
 * \author J. T. Whelan <john.whelan@ligo.org>
 *
 * \brief Raises an \c LALUnit structure to a specified rational power.
 *
 * This function raises the \c LALUnit structure <tt>*input</tt> to
 * the rational power <tt>*power</tt>.  In this way, units such as
 * \f$\mathrm{s}^{1/2}\f$ and \f$\mathrm{m}^{-1}\f$ can be created using existing units.
 *
 * ### Algorithm ###
 *
 * The function first multiplies the overall power of ten
 * <tt>input-\>powerOfTen</tt> by the rational number <tt>*power</tt>,
 * checking to make sure that the resulting power is still an integer.
 * It then multiplies each of the rational powers in <tt>*input</tt> by
 * <tt>*power</tt> by naïve multiplication of rational numbers
 * \f[
 * \left(\frac{N_1}{1+D_1}\right)\left( \frac{N_2}{1+D_2} \right)
 * = \frac{N_1 N_2}{1 + (1+D_1)(1+D_2)-1}
 * \f]
 * and then calls <tt>LALUnitNormalize()</tt> to bring the result into
 * standard form.
 *
 */
/*@{*/

/**
 * Raises a ::LALUnit structure to a rational power given by the ::RAT4 structure \c power.
 */
LALUnit * XLALUnitRaiseRAT4( LALUnit *output, const LALUnit *input,
    const RAT4 *power )
{
  LALUnit     unReduced;
  UINT2       i;
  INT4        numer;
  UINT4       denom, denom1, denom2;

  if ( ! output || ! input || ! power )
    XLAL_ERROR_NULL( XLAL_EFAULT );

  denom2 = power->denominatorMinusOne + 1;

  if ( input->powerOfTen % denom2 )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  numer = (input->powerOfTen / (INT4) denom2) * power->numerator;

  if ( numer >= 32767L || numer <= -32768L )
    XLAL_ERROR_NULL( XLAL_ERANGE );

  unReduced.powerOfTen = numer;

  for (i=0; i<LALNumUnits; ++i) {
    denom1 = 1 + input->unitDenominatorMinusOne[i];
    denom = denom1 * denom2;

    if ( denom - 1 >= 65535L )
      XLAL_ERROR_NULL( XLAL_ERANGE );

    unReduced.unitDenominatorMinusOne[i] = denom - 1;

    numer = input->unitNumerator[i] * power->numerator;

    if ( numer >= 32767L || numer <= -32768L )
      XLAL_ERROR_NULL( XLAL_ERANGE );

    unReduced.unitNumerator[i] = numer;
  } /* for i */

  *output = unReduced;
  if ( XLALUnitNormalize( output ) == XLAL_FAILURE )
    XLAL_ERROR_NULL( XLAL_EFUNC );

  return output;
}

/**
 * Raises a ::LALUnit structure to an integer power \c power.
 */
LALUnit * XLALUnitRaiseINT2( LALUnit *output, const LALUnit *input,
    INT2 power )
{
  RAT4 pow;
  pow.numerator = power;
  pow.denominatorMinusOne = 0;
  if ( ! XLALUnitRaiseRAT4( output, input, &pow ) )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return output;
}

/**
 * Produces the square of a ::LALUnit structure.
 */
LALUnit * XLALUnitSquare( LALUnit *output, const LALUnit *input )
{
  RAT4 pow;
  pow.numerator = 2;
  pow.denominatorMinusOne = 0;
  if ( ! XLALUnitRaiseRAT4( output, input, &pow ) )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return output;
}

/**
 * Produces the square-root of a ::LALUnit structure.
 */
LALUnit * XLALUnitSqrt( LALUnit *output, const LALUnit *input )
{
  RAT4 pow;
  pow.numerator = 1;
  pow.denominatorMinusOne = 1;
  if ( ! XLALUnitRaiseRAT4( output, input, &pow ) )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return output;
}

/** UNDOCUMENTED */
LALUnit * XLALUnitInvert( LALUnit *output, const LALUnit *input )
{
  RAT4 pow;
  pow.numerator = -1;
  pow.denominatorMinusOne = 0;
  if ( ! XLALUnitRaiseRAT4( output, input, &pow ) )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return output;
}

/**
 * This function raises the \c LALUnit structure <tt>*input</tt> to
 * the rational power <tt>*power</tt>.  In this way, units such as
 * \f$\mathrm{s}^{1/2}\f$ and \f$\mathrm{m}^{-1}\f$ can be created using existing units.
 *
 * \deprecated Use XLALUnitRaise() instead.
 */
void
LALUnitRaise (LALStatus *status, LALUnit *output, const LALUnit *input, const RAT4 *power)
{
  INITSTATUS(status);

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( power != NULL, status, UNITSH_ENULLPPARAM, UNITSH_MSGENULLPPARAM );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  ASSERT( input->powerOfTen % (power->denominatorMinusOne + 1) == 0, status,
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
/*@}*/
