/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon, John Whelan
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

#include <math.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/XLALError.h>

/**
 * \addtogroup UnitCompare_c
 * \author J. T. Whelan <john.whelan@ligo.org>
 *
 * \brief Function to compare two \c LALUnit structures.
 */
/*@{*/

/** Return 1 if a unit is dimensionless, 0 otherwise */
int XLALUnitIsDimensionless(const LALUnit *unit)
{
  int i;

  if(!unit)
    XLAL_ERROR(XLAL_EFAULT);

  for(i = 0; i < LALNumUnits; i++)
    if(unit->unitNumerator[i] != unit->unitDenominatorMinusOne[i])
      return 0;
  return 1;
}


/** Return the unit's prefactor */
REAL8 XLALUnitPrefactor(const LALUnit *unit)
{
  if(!unit)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  return pow(10.0, unit->powerOfTen);
}


/**
 * Return the ratio unit1 / unit2
 */
REAL8 XLALUnitRatio(const LALUnit *unit1, const LALUnit *unit2)
{
  LALUnit tmp;

  if(!unit1 || !unit2)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  XLALUnitDivide(&tmp, unit1, unit2);
  if(XLALUnitIsDimensionless(&tmp))
    return XLALUnitPrefactor(&tmp);
  XLAL_ERROR_REAL8(XLAL_EDIMS);
}



/**
 * Returns 0 if the the normal form of the two unit
 * structures are the same or > 0 if they are different.  It returns
 * \c #XLAL_FAILURE and \c ::xlalErrno is set to \c #XLAL_EFAULT if
 * one of the input pointers is \c NULL.
 *
 * Example:
 * \code
 * if(XLALUnitCompare(&unit1, &unit2)) {
 *   units_are_not_equal();
 * }
 * \endcode
 */
int XLALUnitCompare( const LALUnit *unit1, const LALUnit *unit2 )
{
  LALUnit  unitOne, unitTwo;

  if ( ! unit1 || ! unit2 )
    XLAL_ERROR( XLAL_EFAULT );

  unitOne = *unit1;
  unitTwo = *unit2;

  /* normalize the units */
  if ( XLALUnitNormalize( &unitOne ) )
    XLAL_ERROR( XLAL_EFUNC );
  if ( XLALUnitNormalize( &unitTwo ) )
    XLAL_ERROR( XLAL_EFUNC );

  /* factors of 10 disagree? */
  if ( unitOne.powerOfTen != unitTwo.powerOfTen )
    return 1;

  /* powers of dimensions disagree? use memcmp() to compare the arrays */
  if ( memcmp( unitOne.unitNumerator, unitTwo.unitNumerator, LALNumUnits * sizeof( *unitOne.unitNumerator ) ) )
    return 1;
  if ( memcmp( unitOne.unitDenominatorMinusOne, unitTwo.unitDenominatorMinusOne, LALNumUnits * sizeof( *unitOne.unitDenominatorMinusOne ) ) )
    return 1;

  /* agree in all possible ways */
  return 0;
}

/**
 * Compare two units structures, returning true if they are equivalent, false otherwise.
 * This function determines whether the units represented by
 * <tt>*(input->unitOne)</tt> and <tt>*(input->unitTwo)</tt> are the same (both
 * dimensionally and in the power-of-ten prefactor).  In this way,
 * programs and programmers can verify that quantities have the expected
 * units.
 *
 * ### Algorithm ###
 *
 * The function first uses <tt>LALUnitNormalize()</tt> to bring both unit
 * structures into standard form, then compares the powers of ten and the
 * numerator and denominator of each exponent of a fundamental unit in
 * turn.
 *
 * ### Uses ###
 *
 * <tt>LALUnitNormalize()</tt>
 *
 * \deprecated Use XLALUnitCompare() instead.
 *
 */
void
LALUnitCompare (LALStatus *status, BOOLEAN *output, const LALUnitPair *input)
{
  int code;

  INITSTATUS(status);

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );
  XLAL_PRINT_DEPRECATION_WARNING("XLALUnitCompare");

  code = XLALUnitCompare( input->unitOne, input->unitTwo );
  if ( code < 0 )
  {
    XLALClearErrno();
    ABORTXLAL( status );
  }

  *output = code ? 0 : 1;

  RETURN(status);
}
/*@}*/
