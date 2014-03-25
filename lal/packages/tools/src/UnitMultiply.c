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
 * \addtogroup UnitMultiply_c
 * \author J. T. Whelan <john.whelan@ligo.org>
 *
 * \brief Multiplies two \c LALUnit structures.
 *
 * This function multiplies together the \c LALUnit structures
 * <tt>*(input-\>unitOne)</tt> and <tt>*(input-\>unitTwo)</tt>, thus allowing a
 * module to <em>e.g.</em>, multiply two \c REAL8TimeSeries and
 * give the resulting \c REAL8TimeSeries the correct units.
 *
 * ### Algorithm ###
 *
 * The function first adds together the overall powers of ten in the two
 * input unit structures, then adds each of the corresponding rational
 * powers in <tt>*(input-\>unitOne)</tt> and <tt>*(input-\>unitTwo)</tt> by naïve
 * addition of rational numbers
 * \f[
 * \frac{N_1}{1+D_1} + \frac{N_2}{1+D_2} =
 * \frac{N_1 (1+D_2) +  N_2(1+D_1)}{1 + (1+D_1)(1+D_2)-1}
 * \f]
 * and then calls <tt>LALUnitNormalize()</tt> to bring the result into
 * standard form.
 *
 */
/*@{*/

/**
 * This function multiplies together the \c LALUnit structures
 * <tt>*(input->unitOne)</tt> and <tt>*(input->unitTwo)</tt>, thus allowing a
 * module to eg, multiply two \c REAL8TimeSeries and
 * give the resulting \c REAL8TimeSeries the correct units.
 *
 * ### Uses ###
 *
 * <tt>XLALUnitNormalize()</tt>
 *
 */
LALUnit * XLALUnitMultiply( LALUnit *output, const LALUnit *unit1, const LALUnit *unit2 )
{
  LALUnit     unReduced;
  UINT2        i;
  INT4         numer;
  UINT4        denom, denom1, denom2;

  if ( ! output || ! unit1 || ! unit2 )
    XLAL_ERROR_NULL( XLAL_EFAULT );

  numer = unit1->powerOfTen + unit2->powerOfTen;
  if ( numer >= 32767L || numer <= -32768L )
    XLAL_ERROR_NULL( XLAL_ERANGE );

  unReduced.powerOfTen = numer;
  for (i=0; i<LALNumUnits; ++i) {
    denom1 = 1 + unit1->unitDenominatorMinusOne[i];
    denom2 = 1 + unit2->unitDenominatorMinusOne[i];
    denom = denom1 * denom2;

    if ( denom >= 65535L )
      XLAL_ERROR_NULL( XLAL_ERANGE );

    /* One could use the gcd function to find the common factors of
       denom1 and denom2, but we have to reduce the fractions after
       addition anyway; consider e.g., 1/6 + 1/10 = (5+3)/30 = 4/15 */
    unReduced.unitDenominatorMinusOne[i] = denom - 1;
    numer = ((INT4) denom2) * unit1->unitNumerator[i]
      + ((INT4) denom1) * unit2->unitNumerator[i];

    if ( numer >= 32767L || numer <= -32768L )
      XLAL_ERROR_NULL( XLAL_ERANGE );

    unReduced.unitNumerator[i] = numer;
  } /* for i */

  *output = unReduced;
  if ( XLALUnitNormalize( output ) == XLAL_FAILURE )
    XLAL_ERROR_NULL( XLAL_EFUNC );

  return output;
}

/** UNDOCUMENTED */
LALUnit * XLALUnitDivide( LALUnit *output, const LALUnit *unit1, const LALUnit *unit2 )
{
  LALUnit scratch;
  /* invert unit2 and then multiply by unit1 */
  if ( ! XLALUnitInvert( &scratch, unit2 ) )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  if ( ! XLALUnitMultiply( output, unit1, &scratch ) )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return output;
}

/**
 * DEPRECATED.
 * \deprecated Use XLALUnitMultiply() instead.
 */
void
LALUnitMultiply (LALStatus *status, LALUnit *output, const LALUnitPair *input)

{
  XLAL_PRINT_DEPRECATION_WARNING("XLALUnitMultiply");
  INITSTATUS(status);

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
/*@}*/
