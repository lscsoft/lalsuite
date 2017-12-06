/*
*  Copyright (C) 2011 Kipp Cannon
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


#include <stdio.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

#define MAX_FRACTIONAL_ERROR 1e-15

#if defined(NDEBUG) || defined(LAL_NDEBUG) /* debugging is turned off */
int main( void )
{
	return 77; /* don't do any testing */
}
#else

static int compare_within_fraction(double a, double b, double max_fractional_difference)
{
	if(fabs(a - b) / fabs(a > b ? a : b) <= fabs(max_fractional_difference))
		return 0;
	return a > b ? +1 : -1;
}


#define REQUIRE_EQUAL(a, b) \
	if(compare_within_fraction(a, b, MAX_FRACTIONAL_ERROR)) { \
		fprintf(stderr, #a " != " #b " within %g (" #a " = %.17g, " #b " = %.17g)\n", MAX_FRACTIONAL_ERROR, a, b); \
		return 1; \
	}

int main( void )
{
	if(lalNoDebug) /* library was not compiled with debugging */
		return 77; /* don't do any testing */

	REQUIRE_EQUAL(sqrt(2), LAL_SQRT2);
	REQUIRE_EQUAL(1/sqrt(2), LAL_SQRT1_2);
	REQUIRE_EQUAL(exp(LAL_GAMMA), LAL_EXPGAMMA);
	REQUIRE_EQUAL(2 * LAL_PI, LAL_TWOPI);
	REQUIRE_EQUAL(LAL_PI / 2, LAL_PI_2);
	REQUIRE_EQUAL(LAL_PI / 4, LAL_PI_4);
	REQUIRE_EQUAL(1 / LAL_PI, LAL_1_PI);
	REQUIRE_EQUAL(2 / LAL_PI, LAL_2_PI);
	REQUIRE_EQUAL(2 / sqrt(LAL_PI), LAL_2_SQRTPI);
	REQUIRE_EQUAL(LAL_PI / 180, LAL_PI_180);
	REQUIRE_EQUAL(180 / LAL_PI, LAL_180_PI);

	REQUIRE_EQUAL(LAL_H_SI / LAL_TWOPI, LAL_HBAR_SI);

	REQUIRE_EQUAL(sqrt(LAL_HBAR_SI * LAL_C_SI / LAL_G_SI), LAL_MPL_SI);
	REQUIRE_EQUAL(sqrt(LAL_HBAR_SI * LAL_G_SI / pow(LAL_C_SI, 3)), LAL_LPL_SI);
	REQUIRE_EQUAL(sqrt(LAL_HBAR_SI * LAL_G_SI / pow(LAL_C_SI, 5)), LAL_TPL_SI);

	REQUIRE_EQUAL(LAL_MSUN_SI / LAL_MPL_SI * LAL_LPL_SI, LAL_MRSUN_SI);
	REQUIRE_EQUAL(LAL_MSUN_SI / LAL_MPL_SI * LAL_TPL_SI, LAL_MTSUN_SI);

	return 0;
}
#endif
