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

/************************************ <lalVerbatim file="UnitCompareCV">
Author: J. T. Whelan <john.whelan@ligo.org>
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
\texttt{*(input->unitOne)} and \texttt{*(input->unitTwo)} are the same (both
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

#include <math.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/XLALError.h>

NRCSID( UNITCOMPAREC, "$Id$" );


/* Return 1 if a unit is dimensionless, 0 otherwise */
int XLALUnitIsDimensionless(const LALUnit *unit)
{
  int i;

  if(!unit)
    XLAL_ERROR("XLALUnitIsDimensionless", XLAL_EFAULT);

  for(i = 0; i < LALNumUnits; i++)
    if(unit->unitNumerator[i] != unit->unitDenominatorMinusOne[i])
      return 0;
  return 1;
}


/* Return the unit's prefactor */
REAL8 XLALUnitPrefactor(const LALUnit *unit)
{
  if(!unit)
    XLAL_ERROR_REAL8("XLALUnitPrefactor", XLAL_EFAULT);
  return pow(10.0, unit->powerOfTen);
}


/**
 * Return the ratio unit1 / unit2
 */

REAL8 XLALUnitRatio(const LALUnit *unit1, const LALUnit *unit2)
{
  static const char func[] = "XLALUnitRatio";
  LALUnit tmp;

  if(!unit1 || !unit2)
    XLAL_ERROR_REAL8(func, XLAL_EFAULT);

  XLALUnitDivide(&tmp, unit1, unit2);
  if(XLALUnitIsDimensionless(&tmp))
    return XLALUnitPrefactor(&tmp);
  XLAL_ERROR_REAL8(func, XLAL_EDIMS);
}


/**
 * Returns 0 if units are the same (after normalization), or > 0 if they
 * are different.  Returns < 0 on error.
 *
 * Example:
 *
 * if(XLALUnitCompare(&unit1, &unit2)) {
 * 	units_are_not_equal();
 * }
 */

int XLALUnitCompare( const LALUnit *unit1, const LALUnit *unit2 )
{
  static const char func[] = "XLALUnitCompare";
  LALUnit  unitOne, unitTwo;

  if ( ! unit1 || ! unit2 )
    XLAL_ERROR( func, XLAL_EFAULT );

  unitOne = *unit1;
  unitTwo = *unit2;

  /* normalize the units */
  if ( XLALUnitNormalize( &unitOne ) )
    XLAL_ERROR( func, XLAL_EFUNC );
  if ( XLALUnitNormalize( &unitTwo ) )
    XLAL_ERROR( func, XLAL_EFUNC );

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


/* <lalVerbatim file="UnitCompareCP"> */
void
LALUnitCompare (LALStatus *status, BOOLEAN *output, const LALUnitPair *input)
/* </lalVerbatim> */
     /* Compare two units structures, returning 0 if they differ */
{
  int code;

  INITSTATUS( status, "LALUnitCompare", UNITCOMPAREC );

  ASSERT( input != NULL, status, UNITSH_ENULLPIN, UNITSH_MSGENULLPIN );

  ASSERT( output != NULL, status, UNITSH_ENULLPOUT, UNITSH_MSGENULLPOUT );
  XLALPrintDeprecationWarning( "LALUnitCompare", "XLALUnitCompare" );

  code = XLALUnitCompare( input->unitOne, input->unitTwo );
  if ( code < 0 )
  {
    XLALClearErrno();
    ABORTXLAL( status );
  }

  *output = code ? 0 : 1;

  RETURN(status);
}
