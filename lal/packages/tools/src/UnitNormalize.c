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
 * \addtogroup UnitNormalize_c
 * \author J. T. Whelan <john.whelan@ligo.org>
 *
 * \brief Brings an \c LALUnit structure into standard form by reducing all of the rational exponents into LCD form.
 *
 * Since the \c LALUnit structure stores the rational powers of the
 * fundamental units as numerator and denominator, it is possible to
 * represent the same units in different ways, <em>e.g.</em>, \f$\mathrm{m}^2\f$
 * versus \f$\mathrm{m}^{4/2}\f$.  This function reduces all of those fractions to
 * convert the structure to its simplest form.
 *
 * ### Algorithm ###
 *
 * The rational powers are reduced using Euclid's algorithm \cite Geddes_1992.
 *
 * ### Notes ###
 *
 * Note that the functions <tt>LALUnitRaise()</tt>,
 * <tt>LALUnitMultiply()</tt>, and <tt>LALUnitCompare()</tt> all call
 * <tt>LALUnitNormalize()</tt> themselves, so there is usually no need to
 * call it explicitly.
 *
 */
/*@{*/

/* this is an implementation of the Euclidean Algorithm */
static UINT2
gcd(INT2 numer, UINT2 denom)
{
   UINT2 next_numer,next_denom,remainder;
   next_numer=abs(numer);
   next_denom=denom;
   while(next_denom != 0){
      remainder=next_numer%next_denom;
      next_numer=next_denom;
      next_denom=remainder;
   }
   return abs(next_numer);
}

/**
 * Returns 0 upon success or \c #XLAL_FAILURE
 * if the input pointer is \c NULL, in which case \c ::xlalErrno
 * is set to \c #XLAL_EFAULT.
 */
int XLALUnitNormalize( LALUnit *unit )
{
  UINT2 commonFactor;
  UINT2 i;

  if ( ! unit )
    XLAL_ERROR( XLAL_EFAULT );

  for (i=0; i<LALNumUnits; ++i)
  {
    commonFactor = gcd ( unit->unitNumerator[i], unit->unitDenominatorMinusOne[i] + 1 );
    unit->unitNumerator[i] /= commonFactor;
    unit->unitDenominatorMinusOne[i] = ( unit->unitDenominatorMinusOne[i] + 1 ) / commonFactor - 1;
  } /* for i */

  return 0;
}


/**
 * Reduce all the rational powers in an LALUnit structure
 * \deprecated Use XLALUnitNormalize() instead.
 */
void
LALUnitNormalize (LALStatus *status, LALUnit *output, const LALUnit *input)
{

  INITSTATUS(status);

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );

  *output = *input;
  if ( XLALUnitNormalize( output ) == XLAL_FAILURE )
  {
    ABORTXLAL( status );
  }

  RETURN( status );
}
/*@}*/
