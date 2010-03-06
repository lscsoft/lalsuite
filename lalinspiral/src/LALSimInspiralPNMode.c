/*
 * Copyright (C) 2008 J. Creighton
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
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/LALComplex.h>
#include <lal/LALConstants.h>

NRCSID(LALSIMINSPIRALPNMODEC, "$Id$");

/**
 * Computes h(2,2) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (79) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16 XLALSimInspiralPNMode22(
		REAL8 x,      /**< post-Newtonian parameter */
	       	REAL8 phi,    /**< orbital phase */
	       	REAL8 logx,   /**< log(x/x0) tail gauge parameter */
	       	REAL8 m1,     /**< mass of companion 1 */
	       	REAL8 m2,     /**< mass of companion 2 */
		REAL8 r,      /**< distance of source */
		int O         /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNMode22";
	REAL8 fac = -8.0*sqrt(LAL_PI/5.0)*LAL_G_SI*pow(LAL_C_SI, -2.0);
	REAL8 m = m1 + m2;
	REAL8 mu = m1*m2/m;
	REAL8 nu = mu/m;
	COMPLEX16 ans;
	REAL8 re = 0.0;
	REAL8 im = 0.0;
	switch (O) {
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
		case -1: /* use highest available pN order */
		case 6:
			re += ((27027409.0/646800.0) - (856.0/105.0)*LAL_GAMMA + (2.0/3.0)*pow(LAL_PI, 2.0) - (1712.0/105.0)*log(2.0) - (428.0/105.0)*log(x) - 18.0*pow(logx, 2.0) - ((278185.0/33264.0) - (41.0/96.0)*pow(LAL_PI, 2.0))*nu - (20261.0/2772.0)*pow(nu, 2.0) + (114635.0/99792.0)*pow(nu, 3.0))*pow(x, 3.0);
			im += ((428.0/105.0) + 12.0*logx)*LAL_PI*pow(x, 3.0);
		case 5:
			re -= ((107.0/21.0) - (34.0/21.0)*nu)*LAL_PI*pow(x, 2.5);
			im -= (24.0*nu + ((107.0/7.0) - (34.0/7.0)*nu)*logx)*pow(x, 2.5);
		case 4:
			re -= ((2173.0/1512.0) + (1069.0/216.0)*nu - (2047.0/1512.0)*pow(nu, 2.0))*pow(x, 2.0);
		case 3:
			re += 2.0*LAL_PI*pow(x, 1.5);
			im += 6.0*logx*pow(x, 1.5);
		case 2:
			re -= ((107.0/42.0) - (55.0/42.0)*nu)*x;
		case 1:
		case 0:
			re += 1.0;
	}
	ans = cmul(cpolar(1.0, -2.0*phi), crect(re, im));
	ans = cmulr(ans, (fac*nu*m/r)*x);
	return ans;
}

/**
 * Computes h(2,1) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (80) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16 XLALSimInspiralPNMode21(
		REAL8 x,      /**< post-Newtonian parameter */
	       	REAL8 phi,    /**< orbital phase */
	       	REAL8 logx,   /**< log(x/x0) tail gauge parameter */
	       	REAL8 m1,     /**< mass of companion 1 */
	       	REAL8 m2,     /**< mass of companion 2 */
		REAL8 r,      /**< distance of source */
		int O         /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNMode21";
	REAL8 fac = -8.0*sqrt(LAL_PI/5.0)*LAL_G_SI*pow(LAL_C_SI, -2.0);
	REAL8 m = m1 + m2;
	REAL8 dm = m1 - m2;
	REAL8 mu = m1*m2/m;
	REAL8 nu = mu/m;
	COMPLEX16 ans;
	REAL8 re = 0.0;
	REAL8 im = 0.0;
	switch (O) {
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
		case -1: /* use highest available pN order */
		case 5:
			re -= ((43.0/126.0) + (509.0/126.0)*nu - (79.0/168.0)*pow(nu, 2.0))*pow(x, 2.0);
		case 4:
			re += LAL_PI*pow(x, 1.5);
			im += (-(1.0/2.0) - 2.0*log(2.0) + 3.0*logx)*pow(x, 1.5);
		case 3:
			re -= ((17.0/28.0) - (5.0/7.0)*nu)*x;
		case 2:
		case 1:
			re += 1.0;
		case 0:
			re += 0.0;
	}
	ans = cmul(cpolar(1.0, -phi), crect(re, im));
	ans = cmuli(ans, (fac*nu*dm/r)*pow(x, 1.5));
	return ans;
}

/**
 * Computes h(3,3) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (82) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16 XLALSimInspiralPNMode33(
		REAL8 x,      /**< post-Newtonian parameter */
	       	REAL8 phi,    /**< orbital phase */
	       	REAL8 logx,   /**< log(x/x0) tail gauge parameter */
	       	REAL8 m1,     /**< mass of companion 1 */
	       	REAL8 m2,     /**< mass of companion 2 */
		REAL8 r,      /**< distance of source */
		int O         /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNMode33";
	REAL8 fac = 3.0*sqrt(6.0*LAL_PI/7.0)*LAL_G_SI*pow(LAL_C_SI, -2.0);
	REAL8 m = m1 + m2;
	REAL8 dm = m1 - m2;
	REAL8 mu = m1*m2/m;
	REAL8 nu = mu/m;
	COMPLEX16 ans;
	REAL8 re = 0.0;
	REAL8 im = 0.0;
	switch (O) {
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
		case -1: /* use highest available pN order */
		case 5:
			re += ((123.0/110.0) - (1838.0/165.0)*nu - (887.0/330.0)*pow(nu, 2.0))*pow(x, 2.0);
		case 4:
			re += 3.0*LAL_PI*pow(x, 1.5);
			im += (-(21.0/5.0) + 6.0*log(3.0/2.0) + 9.0*logx)*pow(x, 1.5);
		case 3:
			re -= (4.0 - 2.0*nu)*x;
		case 2:
		case 1:
			re += 1.0;
		case 0:
			re += 0.0;
	}
	ans = cmul(cpolar(1.0, -3.0*phi), crect(re, im));
	ans = cmuli(ans, (fac*nu*dm/r)*pow(x, 1.5));
	return ans;
}


/**
 * Computes h(3,2) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (83) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16 XLALSimInspiralPNMode32(
		REAL8 x,      /**< post-Newtonian parameter */
	       	REAL8 phi,    /**< orbital phase */
	       	REAL8 logx,   /**< log(x/x0) tail gauge parameter */
	       	REAL8 m1,     /**< mass of companion 1 */
	       	REAL8 m2,     /**< mass of companion 2 */
		REAL8 r,      /**< distance of source */
		int O         /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNMode32";
	REAL8 fac = -(8.0/3.0)*sqrt(LAL_PI/7.0)*LAL_G_SI*pow(LAL_C_SI, -2.0);
	REAL8 m = m1 + m2;
	REAL8 mu = m1*m2/m;
	REAL8 nu = mu/m;
	COMPLEX16 ans;
	REAL8 re = 0.0;
	REAL8 im = 0.0;
	switch (O) {
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
		case -1: /* use highest available pN order */
		case 5:
			re += 2.0*LAL_PI*(1.0 - 3.0*nu)*pow(x, 1.5);
			im += (-3.0 + (66.0/5.0)*nu + 6.0*(1.0 - 3.0*nu)*logx)*pow(x, 1.5);
		case 4:
			re -= ((193.0/90.0) - (145.0/18.0)*nu + (73.0/18.0)*pow(nu, 2.0))*x;
		case 3:
		case 2:
			re += 1.0 - 3.0*nu;
		case 1:
		case 0:
			re += 0.0;
	}
	ans = cmul(cpolar(1.0, -2.0*phi), crect(re, im));
	ans = cmulr(ans, (fac*nu*m/r)*pow(x, 2.0));
	return ans;
}

/**
 * Computes h(3,1) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (84) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16 XLALSimInspiralPNMode31(
		REAL8 x,      /**< post-Newtonian parameter */
	       	REAL8 phi,    /**< orbital phase */
	       	REAL8 logx,   /**< log(x/x0) tail gauge parameter */
	       	REAL8 m1,     /**< mass of companion 1 */
	       	REAL8 m2,     /**< mass of companion 2 */
		REAL8 r,      /**< distance of source */
		int O         /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNMode31";
	REAL8 fac = -(1.0/3.0)*sqrt(2.0*LAL_PI/35.0)*LAL_G_SI*pow(LAL_C_SI, -2.0);
	REAL8 m = m1 + m2;
	REAL8 dm = m1 - m2;
	REAL8 mu = m1*m2/m;
	REAL8 nu = mu/m;
	COMPLEX16 ans;
	REAL8 re = 0.0;
	REAL8 im = 0.0;
	switch (O) {
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
		case -1: /* use highest available pN order */
		case 5:
			re += ((607.0/198.0) - (136.0/99.0)*nu - (247.0/198.0)*pow(nu, 2.0))*pow(x, 2.0);
		case 4:
			re += LAL_PI*pow(x, 1.5);
			im += (-(7.0/5.0) - 2.0*log(2.0) + 3.0*logx)*pow(x, 1.5);
		case 3:
			re -= ((8.0/3.0) + (2.0/3.0)*nu)*x;
		case 2:
		case 1:
			re += 1.0;
		case 0:
			re += 0.0;
	}
	ans = cmul(cpolar(1.0, -phi), crect(re, im));
	ans = cmuli(ans, (fac*nu*dm/r)*pow(x, 1.5));
	return ans;
}
