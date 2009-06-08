/*
 * Copyright (C) 2008 J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria
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

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <lal/LALSimInspiral.h>
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/LALComplex.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include "check_series_macros.h"

NRCSID(LALSIMINSPIRALC, "$Id$");

/**
 * Computes the (s)Y(l,m) spin-weighted spherical harmonic.
 *
 * From somewhere ....
 *
 * See also:
 * Implements Equations (II.9)-(II.13) of
 * D. A. Brown, S. Fairhurst, B. Krishnan, R. A. Mercer, R. K. Kopparapu,
 * L. Santamaria, and J. T. Whelan,
 * "Data formats for numerical relativity waves",
 * arXiv:0709.0093v1 (2007).
 *
 * Currently only supports s=-2, l=2,3,4,5 modes.
 */
COMPLEX16 XLALSpinWeightedSphericalHarmonic(
		REAL8 theta,  /**< polar angle (rad) */
	       	REAL8 phi,    /**< azimuthal angle (rad) */
	       	int s,        /**< spin weight */
	       	int l,        /**< mode number l */
	       	int m         /**< mode number m */
		)
{
	static const char *func = "XLALSpinWeightedSphericalHarmonic";
	REAL8 fac;
	COMPLEX16 ans;

	/* sanity checks ... */
	if ( l < abs(s) ) {
		XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |s| <= l\n", func, s, l, m );
		XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
	}
	if ( l < abs(m) ) {
		XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
		XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
	}

	if ( s == -2 ) {
		if ( l == 2 ) {
			switch ( m ) {
				case -2:
					fac = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta ));
					break;
				case -1:
					fac = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 - cos( theta ));
					break;

				case 0:
					fac = sqrt( 15.0 / ( 32.0 * LAL_PI ) ) * sin( theta )*sin( theta );
					break;

				case 1:
					fac = sqrt( 5.0 / ( 16.0 * LAL_PI ) ) * sin( theta )*( 1.0 + cos( theta ));
					break;

				case 2:
					fac = sqrt( 5.0 / ( 64.0 * LAL_PI ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta ));
					break;
				default:
					XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
					XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
					break;
			} /*  switch (m) */
		}  /* l==2*/
		else if ( l == 3 ) {
			switch ( m ) {
				case -3:
					fac = sqrt(21.0/(2.0*LAL_PI))*cos(theta/2.0)*pow(sin(theta/2.0),5.0);
					break;
				case -2:
					fac = sqrt(7.0/4.0*LAL_PI)*(2.0 + 3.0*cos(theta))*pow(sin(theta/2.0),4.0);
					break;
				case -1:
					fac = sqrt(35.0/(2.0*LAL_PI))*(sin(theta) + 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
					break;
				case 0:
					fac = (sqrt(105.0/(2.0*LAL_PI))*cos(theta)*pow(sin(theta),2.0))/4.0;
					break;
				case 1:
					fac = -sqrt(35.0/(2.0*LAL_PI))*(sin(theta) - 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0;
					break;

				case 2:
					fac = sqrt(7.0/LAL_PI)*pow(cos(theta/2.0),4.0)*(-2.0 + 3.0*cos(theta))/2.0;
					break;

				case 3:
					fac = -sqrt(21.0/(2.0*LAL_PI))*pow(cos(theta/2.0),5.0)*sin(theta/2.0);
					break;

				default:
					XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
					XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
					break;
			}
		}   /* l==3 */
		else if ( l == 4 ) {
			switch ( m ) {
				case -4:
					fac = 3.0*sqrt(7.0/LAL_PI)*pow(cos(theta/2.0),2.0)*pow(sin(theta/2.0),6.0);
					break;
				case -3:
					fac = 3.0*sqrt(7.0/(2.0*LAL_PI))*cos(theta/2.0)*(1.0 + 2.0*cos(theta))*pow(sin(theta/2.0),5.0);
					break;

				case -2:
					fac = (3.0*(9.0 + 14.0*cos(theta) + 7.0*cos(2.0*theta))*pow(sin(theta/2.0),4.0))/(4.0*sqrt(LAL_PI));
					break;
				case -1:
					fac = (3.0*(3.0*sin(theta) + 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) - 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*LAL_PI));
					break;
				case 0:
					fac = (3.0*sqrt(5.0/(2.0*LAL_PI))*(5.0 + 7.0*cos(2.0*theta))*pow(sin(theta),2.0))/16.0;
					break;
				case 1:
					fac = (3.0*(3.0*sin(theta) - 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) + 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*LAL_PI));
					break;
				case 2:
					fac = (3.0*pow(cos(theta/2.0),4.0)*(9.0 - 14.0*cos(theta) + 7.0*cos(2.0*theta)))/(4.0*sqrt(LAL_PI));
					break;
				case 3:
					fac = -3.0*sqrt(7.0/(2.0*LAL_PI))*pow(cos(theta/2.0),5.0)*(-1.0 + 2.0*cos(theta))*sin(theta/2.0);
					break;
				case 4:
					fac = 3.0*sqrt(7.0/LAL_PI)*pow(cos(theta/2.0),6.0)*pow(sin(theta/2.0),2.0);
					break;
				default:
					XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
					XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
					break;
			}
		}    /* l==4 */
		else if ( l == 5 ) {
			switch ( m ) {
				case -5:
					fac = sqrt(330.0/LAL_PI)*pow(cos(theta/2.0),3.0)*pow(sin(theta/2.0),7.0);
					break;
				case -4:
					fac = sqrt(33.0/LAL_PI)*pow(cos(theta/2.0),2.0)*(2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),6.0);
					break;
				case -3:
					fac = (sqrt(33.0/(2.0*LAL_PI))*cos(theta/2.0)*(17.0 + 24.0*cos(theta) + 15.0*cos(2.0*theta))*pow(sin(theta/2.0),5.0))/4.0;
					break;
				case -2:
					fac = (sqrt(11.0/LAL_PI)*(32.0 + 57.0*cos(theta) + 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))*pow(sin(theta/2.0),4.0))/8.0;
					break;
				case -1:
					fac = (sqrt(77.0/LAL_PI)*(2.0*sin(theta) + 8.0*sin(2.0*theta) + 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) - 15.0*sin(5.0*theta)))/256.0;
					break;
				case 0:
					fac = (sqrt(1155.0/(2.0*LAL_PI))*(5.0*cos(theta) + 3.0*cos(3.0*theta))*pow(sin(theta),2.0))/32.0;
					break;
				case 1:
					fac = sqrt(77.0/LAL_PI)*(-2.0*sin(theta) + 8.0*sin(2.0*theta) - 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) + 15.0*sin(5.0*theta))/256.0;
					break;
				case 2:
					fac = sqrt(11.0/LAL_PI)*pow(cos(theta/2.0),4.0)*(-32.0 + 57.0*cos(theta) - 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))/8.0;
					break;
				case 3:
					fac = -sqrt(33.0/(2.0*LAL_PI))*pow(cos(theta/2.0),5.0)*(17.0 - 24.0*cos(theta) + 15.0*cos(2.0*theta))*sin(theta/2.0)/4.0;
					break;
				case 4:
					fac = sqrt(33.0/LAL_PI)*pow(cos(theta/2.0),6.0)*(-2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),2.0);
					break;
				case 5:
					fac = -sqrt(330.0/LAL_PI)*pow(cos(theta/2.0),7.0)*pow(sin(theta/2.0),3.0);
					break;
				default:
					XLALPrintError("XLAL Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n", func, s, l, m );
					XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
					break;
			}
		}  /* l==5 */
		else {
			XLALPrintError("XLAL Error - %s: Unsupported mode l=%d (only l in [2,5] implemented)\n", func, s);
			XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
		}
	}
       	else {
		XLALPrintError("XLAL Error - %s: Unsupported mode s=%d (only s=-2 implemented)\n", func, s);
		XLAL_ERROR_VAL(func, XLAL_EINVAL, czero);
	}
	if (m)
		ans = cmulr(cpolar(1.0, m*phi), fac);
	else
		ans = csetr(fac);
	return ans;
}


/**
 * Multiplies a mode h(l,m) by a spin-2 weighted spherical harmonic
 * to obtain hplus - i hcross, which is added to the time series.
 *
 * Implements the sum of a single term of Eq. (11) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 *
 * If sym is non-zero, symmetrically add the m and -m terms assuming
 * that h(l,-m) = (-1)^l h(l,m)*; see Eq. (78) ibid.
 */
int XLALSimAddMode(
		REAL8TimeSeries *hplus,      /**< +-polarization waveform */
	       	REAL8TimeSeries *hcross,     /**< x-polarization waveform */
	       	COMPLEX16TimeSeries *hmode,  /**< complex mode h(l,m) */
	       	REAL8 theta,                 /**< polar angle (rad) */
	       	REAL8 phi,                   /**< azimuthal angle (rad) */
	       	int l,                       /**< mode number l */
	       	int m,                       /**< mode number m */
	       	int sym                      /**< flag to add -m mode too */
		)
{
	static const char *func = "XLALSimAddMode";
	COMPLEX16 Y;
	UINT4 j;

	LAL_CHECK_VALID_SERIES(hmode, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hplus, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hcross, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hmode, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hcross, hmode, XLAL_FAILURE);

	Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
	for ( j = 0; j < hmode->data->length; ++j ) {
		COMPLEX16 hpc;
		hpc = cmul(Y, hmode->data->data[j]);
		hplus->data->data[j] += creal(hpc);
		hcross->data->data[j] += -cimag(hpc);
	}
	if ( sym ) { /* equatorial symmetry: add in -m mode */
		Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m);
		if ( l % 2 ) /* l is odd */
			Y = cneg(Y);
		for ( j = 0; j < hmode->data->length; ++j ) {
			COMPLEX16 hpc;
			hpc = cmul(Y, conj(hmode->data->data[j]));
			hplus->data->data[j] += creal(hpc);
			hcross->data->data[j] += -cimag(hpc);
		}
	}
	return 0;
}


/**
 * Computes the rate of increase of the orbital frequency for a post-Newtonian
 * inspiral.  This function returns dx/dt rather than the true angular
 * acceleration.
 *
 * Implements Equation (6) of
 * Yi Pan, Alessandra Buonanno, Yanbei Chen, and Michele Vallisneri,
 * "A physical template family for gravitational waves from precessing
 * binaries of spinning compact objects: Application to single-spin binaries"
 * arXiv:gr-qc/0310034v3 (2007).
 *
 * Note: this equation is actually dx/dt rather than (domega/dt)/(omega)^2
 * so the leading coefficient is different.  Also, this function applies
 * for non-spinning objects.
 *
 * Compare the overall coefficient, with nu=1/4, to Equation (45) of
 * Michael Boyle, Duncan A. Brown, Lawrence E. Kidder, Abdul H. Mroue,
 * Harald P. Pfeiﬀer, Mark A. Scheel, Gregory B. Cook, and Saul A. Teukolsky
 * "High-accuracy comparison of numerical relativity simulations with
 * post-Newtonian expansions"
 * arXiv:0710.0158v1 (2007).
 */
REAL8 XLALSimInspiralPNAngularAcceleration(
		REAL8 x,  /**< post-Newtonian parameter */
	       	REAL8 m1, /**< mass of companion 1 */
	       	REAL8 m2, /**< mass of companion 2 */
	       	int O     /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNAngularAcceleration";
	const REAL8 theta = (1039.0/4620.0);
	REAL8 m = m1 + m2;
	REAL8 mu = m1*m2/m;
	REAL8 nu = mu/m;
	REAL8 ans = 0;

	switch (O) {
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_REAL8(func, XLAL_EINVAL);
		case -1: /* use highest available pN order */
		case 7:
			ans += (-(4415.0/4032.0) + (358675.0/6048.0)*nu + (91495.0/1512.0)*pow(nu, 2.0))*LAL_PI*pow(x, 3.5);
		case 6:
			ans += ( ( (16447322263.0/139708800.0) - (1712.0/105.0)*LAL_GAMMA + (16.0/3.0)*pow(LAL_PI, 2.0) ) + ( -(273811877.0/1088640.0) + (451.0/48.0)*pow(LAL_PI, 2.0) - (88.0/3.0)*theta )*nu + (541.0/896.0)*pow(nu, 2.0) - (5605.0/2592.0)*pow(nu, 3.0) - (856.0/105.0)*log(16.0*x) )*pow(x, 3.0);
		case 5:
			ans += ( -(4159.0/672.0) - (15876.0/672.0)*nu )*LAL_PI*pow(x, 2.5);
		case 4:
			ans += ( (34103.0/18144.0) + (13661.0/2016.0)*nu + (59.0/18.0)*pow(nu, 2.0) )*pow(x, 2.0);
		case 3:
			ans += 4.0*LAL_PI*pow(x, 1.5);
		case 2:
			ans += ( -(743.0/336.0) - (924.0/336.0)*nu )*x;
		case 1:
		case 0:
			ans += 1.0;
	}
	ans *= ((64.0*pow(LAL_C_SI, 3.0))/(5.0*LAL_G_SI*m))*pow(x, 5.0)*nu;
	return ans;
}


/**
 * Computes the orbital angular velocity from the quantity x.
 * This is from the definition of x.
 *
 * Implements Equation (46) of
 * Michael Boyle, Duncan A. Brown, Lawrence E. Kidder, Abdul H. Mroue,
 * Harald P. Pfeiﬀer, Mark A. Scheel, Gregory B. Cook, and Saul A. Teukolsky
 * "High-accuracy comparison of numerical relativity simulations with
 * post-Newtonian expansions"
 * arXiv:0710.0158v1 (2007).
 */
REAL8 XLALSimInspiralPNAngularVelocity(
		REAL8 x,  /**< post-Newtonian parameter */
	       	REAL8 m1, /**< mass of companion 1 */
	       	REAL8 m2  /**< mass of companion 2 */
		)
{
	REAL8 m = m1 + m2;
	REAL8 ans;
	ans = (pow(LAL_C_SI, 3.0)*pow(x, 1.5))/(LAL_G_SI*m);
	return ans;
}


/**
 * Computes the orbital energy at a fixed frequency and pN order.
 *
 * Implements Equation (152) of
 * Luc Blanchet,
 * "Gravitational Radiation from Post-Newtonian Sources and Inspiralling
 * Compact Binaries",
 * http://www.livingreviews.org/lrr-2006-4/index.html
 *
 * This is the same as Equation (10) (where the spin of the objects
 * is zero) of:
 * Yi Pan, Alessandra Buonanno, Yanbei Chen, and Michele Vallisneri,
 * "A physical template family for gravitational waves from precessing
 * binaries of spinning compact objects: Application to single-spin binaries"
 * arXiv:gr-qc/0310034v3 (2007).
 * Note: this equation is actually dx/dt rather than (domega/dt)/(omega)^2
 * so the leading coefficient is different.
 */
REAL8 XLALSimInspiralPNEnergy(
		REAL8 x,  /**< post-Newtonian parameter */
	       	REAL8 m1, /**< mass of companion 1 */
	       	REAL8 m2, /**< mass of companion 2 */
	       	int O     /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNEnergy";
	const REAL8 lambda = -(1987.0/3080.0);
	REAL8 m = m1 + m2;
	REAL8 mu = m1*m2/m;
	REAL8 nu = mu/m;
	REAL8 ans = 0;
	switch (O) { /* note: fall through at each level */
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_REAL8(func, XLAL_EINVAL);
		case -1: /* use highest available pN order */
		case 7:
		case 6:
			ans += ( -(675.0/64.0) + ((209323.0/4032.0) - (205.0/96.0)*pow(LAL_PI, 2.0) - (110.0/9.0)*lambda)*nu - (155.0/96.0)*pow(nu, 2.0) - (35.0/5184.0)*pow(nu, 3.0) )*pow(x, 3.0);
		case 5:
		case 4:
			ans += ( -(27.0/8.0) + (19.0/8.0)*nu - (1.0/24.0)*pow(nu, 2.0) )*pow(x, 2.0);
		case 3:
		case 2:
			ans += ( -(3.0/4.0) - (1.0/12.0)*nu )*x;
		case 1:
		case 0:
			ans += 1;
	}
	ans *= -0.5*mu*pow(LAL_C_SI, 2.0)*x;
	return ans;
}


/* GSL integrand routine */
struct XLALSimInspiralPNEvolveOrbitParams { REAL8 m1; REAL8 m2; int O; };
static int XLALSimInspiralPNEvolveOrbitTaylorT4Integrand(double t, const double y[], double ydot[], void *params)
{
	struct XLALSimInspiralPNEvolveOrbitParams *p = params;
	ydot[0] = XLALSimInspiralPNAngularAcceleration(y[0], p->m1, p->m2, p->O);
	ydot[1] = XLALSimInspiralPNAngularVelocity(y[0], p->m1, p->m2);
	t = 0.0;
	return GSL_SUCCESS;
}


/**
 * Evolves a post-Newtonian orbit using the Taylor T4 method.
 *
 * See:
 * Michael Boyle, Duncan A. Brown, Lawrence E. Kidder, Abdul H. Mroue,
 * Harald P. Pfeiﬀer, Mark A. Scheel, Gregory B. Cook, and Saul A. Teukolsky
 * "High-accuracy comparison of numerical relativity simulations with
 * post-Newtonian expansions"
 * arXiv:0710.0158v1 (2007).
 */
int XLALSimInspiralPNEvolveOrbitTaylorT4(
		REAL8TimeSeries **x,   /**< post-Newtonian parameter [returned] */
	       	REAL8TimeSeries **phi, /**< orbital phase [returned] */
	       	LIGOTimeGPS *tc,       /**< coalescence time */
	       	REAL8 phic,            /**< coalescence phase */
	       	REAL8 deltaT,          /**< sampling interval */
		REAL8 m1,              /**< mass of companion 1 */
		REAL8 m2,              /**< mass of companion 2 */
		REAL8 fmin,            /**< start frequency */
		int O                  /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNEvolveOrbitTaylorT4";
	const UINT4 blocklen = 1024;
	const REAL8 xisco = 1.0/6.0;
	struct XLALSimInspiralPNEvolveOrbitParams params;
	REAL8 m = m1 + m2;
	REAL8 E;
	UINT4 j;
	double y[2];
	double yerr[2];
	const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
	gsl_odeiv_step *s;
	gsl_odeiv_system sys;

	/* setup ode system */
	sys.function = XLALSimInspiralPNEvolveOrbitTaylorT4Integrand;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = &params;

	/* allocate memory */
	*x = XLALCreateREAL8TimeSeries( "ORBITAL_FREQUENCY_PARAMETER", tc, 0.0, deltaT, &lalDimensionlessUnit, blocklen );
	*phi = XLALCreateREAL8TimeSeries( "ORBITAL_PHASE", tc, 0.0, deltaT, &lalDimensionlessUnit, blocklen );
	if ( !x || !phi )
		XLAL_ERROR(func, XLAL_EFUNC);

	params.m1 = m1;
	params.m2 = m2;
	params.O  = O;

	y[0] = (*x)->data->data[0] = pow(LAL_PI*LAL_G_SI*m*fmin/pow(LAL_C_SI,3.0), 2.0/3.0);
	y[1] = (*phi)->data->data[0] = 0.0;
	E = XLALSimInspiralPNEnergy(y[0], m1, m2, O);
	if (XLALIsREAL8FailNaN(E))
		XLAL_ERROR(func, XLAL_EFUNC);
	j = 0;

	s = gsl_odeiv_step_alloc(T, 2);
	while (1) {
		REAL8 dE;
		gsl_odeiv_step_apply(s, j*deltaT, deltaT, y, yerr, NULL, NULL, &sys);
		/* MECO termination condition */
		dE = -E;
		dE += E = XLALSimInspiralPNEnergy(y[0], m1, m2, O);
		if (XLALIsREAL8FailNaN(E))
			XLAL_ERROR(func, XLAL_EFUNC);
		if ( dE > 0.0 ) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at MECO\n", func);
			break;
		}
		/* ISCO termination condition for quadrupole, 1pN, 2.5pN */
		if ( (O == 0 || O == 1 || O == 2 || O == 5) && y[0] > xisco ) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at ISCO\n", func);
			break;
		}
		(*x)->data->data[j] = y[0];
		(*phi)->data->data[j] = y[1];
		++j;
		if ( j >= (*x)->data->length ) {
			if ( ! XLALResizeREAL8TimeSeries(*x, 0, (*x)->data->length + blocklen) )
				XLAL_ERROR(func, XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*phi, 0, (*phi)->data->length + blocklen) )
				XLAL_ERROR(func, XLAL_EFUNC);
		}
	}
	gsl_odeiv_step_free(s);

	/* make the correct length */
	if ( ! XLALResizeREAL8TimeSeries(*x, 0, j) )
		XLAL_ERROR(func, XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*phi, 0, j) )
		XLAL_ERROR(func, XLAL_EFUNC);

	/* adjust to correct tc and phic */
	XLALGPSAdd(&(*phi)->epoch, -1.0*j*deltaT);
	XLALGPSAdd(&(*x)->epoch, -1.0*j*deltaT);

	phic -= (*phi)->data->data[(*phi)->data->length - 1];
	for (j = 0; j < (*phi)->data->length; ++j)
		(*phi)->data->data[j] += phic;

	return (int)(*x)->data->length;
}


/**
 * Computes h(l,m) mode timeseries of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * See Eqns. (79)-(116) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(
		REAL8TimeSeries *x,   /**< post-Newtonian parameter */
	       	REAL8TimeSeries *phi, /**< orbital phase */
	       	REAL8 x0,             /**< tail-term gauge choice thing (if you don't know, just set it to zero) */
	       	REAL8 m1,             /**< mass of companion 1 */
	       	REAL8 m2,             /**< mass of companion 2 */
	       	REAL8 r,              /**< distance of source */
	       	int O,                /**< twice post-Newtonain order */
	       	int l,                /**< mode number l */
	       	int m                 /**< mode number m */
		)
{
	static const char *func = "XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries";
	COMPLEX16TimeSeries *h;
	UINT4 j;
	LAL_CHECK_VALID_SERIES(x, NULL);
	LAL_CHECK_VALID_SERIES(phi, NULL);
	LAL_CHECK_CONSISTENT_TIME_SERIES(x, phi, NULL);
	h = XLALCreateCOMPLEX16TimeSeries( "H_MODE", &x->epoch, 0.0, x->deltaT, &lalStrainUnit, x->data->length );
	if ( !h )
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	if ( l == 2 && abs(m) == 2 )
		for ( j = 0; j < h->data->length; ++j )
			h->data->data[j] = XLALSimInspiralPNMode22(x->data->data[j], phi->data->data[j], x0 > 0.0 ? log(x->data->data[j]/x0) : 0.0, m1, m2, r, O);
	else if ( l == 2 && abs(m) == 1 )
		for ( j = 0; j < h->data->length; ++j )
			h->data->data[j] = XLALSimInspiralPNMode21(x->data->data[j], phi->data->data[j], x0 > 0.0 ? log(x->data->data[j]/x0) : 0.0, m1, m2, r, O);
	else if ( l == 3 && abs(m) == 3 )
		for ( j = 0; j < h->data->length; ++j )
			h->data->data[j] = XLALSimInspiralPNMode33(x->data->data[j], phi->data->data[j], x0 > 0.0 ? log(x->data->data[j]/x0) : 0.0, m1, m2, r, O);
	else if ( l == 3 && abs(m) == 2 )
		for ( j = 0; j < h->data->length; ++j )
			h->data->data[j] = XLALSimInspiralPNMode32(x->data->data[j], phi->data->data[j], x0 > 0.0 ? log(x->data->data[j]/x0) : 0.0, m1, m2, r, O);
	else if ( l == 3 && abs(m) == 1 )
		for ( j = 0; j < h->data->length; ++j )
			h->data->data[j] = XLALSimInspiralPNMode31(x->data->data[j], phi->data->data[j], x0 > 0.0 ? log(x->data->data[j]/x0) : 0.0, m1, m2, r, O);
	else {
		XLALDestroyCOMPLEX16TimeSeries(h);
		XLALPrintError("XLAL Error - %s: Unsupported mode l=%d, m=%d\n", func, l, m );
		XLAL_ERROR_NULL(func, XLAL_EINVAL);
	}
	if ( m < 0 ) {
		REAL8 sign = l % 2 ? -1.0 : 1.0;
		for ( j = 0; j < h->data->length; ++j )
			h->data->data[j] = cmulr(conj(h->data->data[j]), sign);
	}
	return h;
}


/**
 * Given an orbit evolution phasing, construct the waveform h+ and hx.
 *
 * Implements Equation (11) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
int XLALSimInspiralPNPolarizationWaveforms(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform [returned] */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform [returned] */
	       	REAL8TimeSeries *x,       /**< post-Newtonian parameter */
	       	REAL8TimeSeries *phi,     /**< orbital phase */
	       	REAL8 x0,                 /**< tail-term gauge choice thing (if you don't know, just set it to zero) */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNPolarizationWaveforms";
	int l, m;
	LAL_CHECK_VALID_SERIES(x, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(phi, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(x, phi, XLAL_FAILURE);
	*hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &x->epoch, 0.0, x->deltaT, &lalStrainUnit, x->data->length );
	*hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &x->epoch, 0.0, x->deltaT, &lalStrainUnit, x->data->length );
	if ( ! hplus || ! hcross )
		XLAL_ERROR(func, XLAL_EFUNC);
	memset((*hplus)->data->data, 0, (*hplus)->data->length*sizeof(*(*hplus)->data->data));
	memset((*hcross)->data->data, 0, (*hcross)->data->length*sizeof(*(*hcross)->data->data));
	for ( l = 2; l <= LAL_PN_MODE_L_MAX; ++l ) {
		for ( m = 1; m <= l; ++m ) {
			COMPLEX16TimeSeries *hmode;
			hmode = XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(x, phi, x0, m1, m2, r, O, l, m);
			if ( ! hmode )
				XLAL_ERROR(func, XLAL_EFUNC);
			if ( XLALSimAddMode(*hplus, *hcross, hmode, i, 0.0, l, m, 1) < 0 )
				XLAL_ERROR(func, XLAL_EFUNC);
			XLALDestroyCOMPLEX16TimeSeries(hmode);
		}
	}
	return 0;
}

/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine allows the user to specify different pN orders
 * for phasing calcuation vs. amplitude calculations.
 */
int XLALSimInspiralPNGenerator(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 x0,                 /**< tail-term gauge choice thing (if you don't know, just set it to zero) */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 fmin,               /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int amplitudeO,           /**< twice post-Newtonian amplitude order */
	       	int phaseO                /**< twice post-Newtonian phase order */
		)
{
	static const char *func = "XLALSimInspiralPNGenerator";
	REAL8TimeSeries *x;
	REAL8TimeSeries *phi;
	int status;
	int n;
	n = XLALSimInspiralPNEvolveOrbitTaylorT4(&x, &phi, tc, phic, deltaT, m1, m2, fmin, phaseO);
	if ( n < 0 )
		XLAL_ERROR(func, XLAL_EFUNC);
	status = XLALSimInspiralPNPolarizationWaveforms(hplus, hcross, x, phi, x0, m1, m2, r, i, amplitudeO);
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(x);
	if ( status < 0 )
		XLAL_ERROR(func, XLAL_EFUNC);
	return n;
}

/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine uses the same pN order for phasing and amplitude
 * (unless the order is -1 in which case the highest available
 * order is used for both of these -- which might not be the same).
 *
 * Log terms in amplitudes are ignored.  This is a gauge choice.
 */
int XLALSimInspiralPN(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 fmin,               /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian order */
		)
{
	/* set x0=0 to ignore log terms */
	return XLALSimInspiralPNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, fmin, r, i, O, O);
}

/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 *
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 *
 * Log terms in amplitudes are ignored.  This is a gauge choice.
 */
int XLALSimInspiralPNRestricted(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 fmin,               /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian phase order */
		)
{
	/* use Newtonian order for amplitude */
	/* set x0=0 to ignore log terms */
	return XLALSimInspiralPNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, fmin, r, i, O, 0);
}


#if 0
#include <lal/PrintFTSeries.h>
int lalDebugLevel = 7;
int main(void)
{
	LIGOTimeGPS tc = { 888888888, 222222222 };
	REAL8 phic = 1.0;
	REAL8 deltaT = 1.0/16384.0;
	REAL8 m1 = 1.4*LAL_MSUN_SI;
	REAL8 m2 = 1.4*LAL_MSUN_SI;
	REAL8 r = 1e6*LAL_PC_SI;
	REAL8 i = 0.5*LAL_PI;
	REAL8 fmin = 100.0;
	int O = -1;
	REAL8TimeSeries *hplus;
	REAL8TimeSeries *hcross;
	XLALSimInspiralPN(&hplus, &hcross, &tc, phic, deltaT, m1, m2, fmin, r, i, O);
	LALDPrintTimeSeries(hplus, "hp.dat");
	LALDPrintTimeSeries(hcross, "hc.dat");
	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);
	LALCheckMemoryLeaks();
	return 0;
}
#endif
