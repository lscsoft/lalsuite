/*
*  Copyright (C) 2011 Jolien Creighton
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

#include <lal/LALConstants.h>
#include <lal/LALDetectors.h>
#include <lal/Units.h>
#include <lal/LALSimSGWB.h>

/**
 * @addtogroup LALSimSGWBORF_c
 * @brief Routines to compute the Overlap Reduction Function for stochastic
 * background gravitational waves between two detectors.
 * @{
 */

/**
 * Computes the overlap reduction function between two detectors at a specified
 * frequency.
 *
 * Implements the formulae given in Allen & Romano (1999).
 */
double XLALSimSGWBOverlapReductionFunction(
	double f,			/**< [in] frequency (Hz) */
	const LALDetector *detector1,	/**< [in] 1st detector */
	const LALDetector *detector2	/**< [in] 2nd detector */
)
{
	double d;
	double s[3];
	double alpha, alphasq, sinalpha, cosalpha;
	double j_0, j_1, j_2;
	double rho1, rho2, rho3;
	double dd, sdds, sds1, sds2;
	double gam;
	size_t i, j;

	/* compute vector between detectors */
	s[0] = detector2->location[0] - detector1->location[0];
	s[1] = detector2->location[1] - detector1->location[1];
	s[2] = detector2->location[2] - detector1->location[2];
	d = sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);

	/* if d is zero (less than 1 meter), detectors are at the same site */
	if (d < 1.0) {
		dd = 0.0;
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
				dd += detector1->response[i][j] * detector2->response[i][j];
		gam = 2.0*dd;
		return gam;
	}

	/* change s[] into a unit vector */
	s[0] /= d;
	s[1] /= d;
	s[2] /= d;

	/* Eq. (3.32) of Allen and Romano (1999) */
	alpha = 2.0*LAL_PI*f*d/LAL_C_SI;

	alphasq = alpha*alpha;
	sinalpha = sin(alpha);
	cosalpha = cos(alpha);

	/* Eq. (3.40) of Allen and Romano (1999) */
	j_0 = sinalpha/alpha;
	j_1 = sinalpha/alphasq - cos(alpha)/alpha;
	j_2 = 3.0*sinalpha/(alpha*alphasq) - 3.0*cosalpha/alphasq - sinalpha/alpha;

	/* Eq. (3.44) of Allen and Romano (1999) */
	rho1 = 0.5*( 10.0*alphasq*j_0 - 20.0*alpha*j_1 +  10.0*j_2)/alphasq;
	rho2 = 0.5*(-20.0*alphasq*j_0 + 80.0*alpha*j_1 - 100.0*j_2)/alphasq;
	rho3 = 0.5*(  5.0*alphasq*j_0 - 50.0*alpha*j_1 + 175.0*j_2)/alphasq;

	/* Compute d1:d2, (s.d1).(d2.s), and (s.d1.s)(s.d2.s) */
	dd = sdds = sds1 = sds2 = 0.0;
	for (i = 0; i < 3; ++i) {
		double sd1 = 0.0;
		double sd2 = 0.0;
		for (j = 0; j < 3; ++j) {
			dd  += detector1->response[i][j] * detector2->response[i][j];
			sd1 += detector1->response[i][j] * s[j];
			sd2 += detector2->response[i][j] * s[j];
		}
		sds1 += sd1 * s[i];
		sds2 += sd2 * s[i];
		sdds += sd1 * sd2;
	}

	/* Eq (3.43) of Allen and Romano (1999) */
	gam = rho1*dd + rho2*sdds + rho3*sds1*sds2;

	return gam;
}

/** @} */

/*
 *
 * TEST CODE
 *
 */

#if 0
#include <stdio.h>

/*
 * Computes the overlap reduction function between the Hanford and Livingston,
 * Virgo, and TAMA for comparison with Allen & Romano.
 */
int test_orf(void)
{
	LALDetector H = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
	LALDetector L = lalCachedDetectors[LAL_LLO_4K_DETECTOR];
	LALDetector V = lalCachedDetectors[LAL_VIRGO_DETECTOR];
	LALDetector T = lalCachedDetectors[LAL_TAMA_300_DETECTOR];
	double f;
	FILE *fp;

	fp = fopen("orf.dat", "w");
	for (f = 1.0; f < 10000.0; f *= 1.01)
		fprintf(fp, "%f %f %f %f\n", f,
			XLALSimSGWBOverlapReductionFunction(f, &H, &L),
			XLALSimSGWBOverlapReductionFunction(f, &H, &V),
			XLALSimSGWBOverlapReductionFunction(f, &H, &T)
			);
	fclose(fp);

	return 0;
}

int main(void)
{
	XLALSetErrorHandler(XLALAbortErrorHandler);
	test_orf();
	LALCheckMemoryLeaks();
	return 0;
}
#endif
