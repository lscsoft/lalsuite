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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

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
 * Given time series for a binary's orbital dynamical variables, 
 * construct the waveform polarizations h+ and hx as a sum of 
 * -2 spin-weighted spherical harmonic modes, h_lm.
 * NB: Valid only for non-precessing systems!
 *
 * Implements Equation (11) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 *
 * FIXME: change the PN variable from x to v = \sqrt{x}
 */
int XLALSimInspiralPNPolarizationWaveformsFromModes(
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
 * Given time series for a binary's orbital dynamical variables, 
 * construct the waveform polarizations h+ and hx directly.
 * NB: Valid only for non-precessing binaries!
 *
 * Implements Equations (5.7) - (5.10) of:
 * K. G. Arun, Bala R Iyer, Moh'd S S Qusailah, "The 2.5PN gravitational wave 
 * polarisations from inspiralling compact binaries in circular orbits"
 * Class. Quant. Grav. 21 (2004) 3771-3802; Erratum-ibid. 22 (2005) 3115
 * arXiv: gr-qc/0404085
 * 
 * Note however, that we do not include the constant "memory" terms
 */
int XLALSimInspiralPNPolarizationWaveforms(
	REAL8TimeSeries **hplus,  /**< +-polarization waveform [returned] */
	REAL8TimeSeries **hcross, /**< x-polarization waveform [returned] */
	REAL8TimeSeries *V,       /**< post-Newtonian (PN) parameter */
	REAL8TimeSeries *Phi,     /**< orbital phase */
	REAL8 x0,                 /**< tail-term gauge choice (default = 0) */
	REAL8 m1,                 /**< mass of companion 1 (kg) */
	REAL8 m2,                 /**< mass of companion 2 (kg) */
	REAL8 r,                  /**< distance of source (m) */
	REAL8 i,                  /**< inclination of source (rad) */
	int ampO                  /**< twice PN order of the amplitude */
	)
{
    REAL8 M, eta, eta2, dm, dist, ampfac, phi, phiShift, v, v2, v3, v4, v5;
    REAL8 hp0, hp05, hp1, hp15, hp2, hp25;
    REAL8 hc0, hc05, hc1, hc15, hc2, hc25;
    REAL8 ci, si, ci2, ci4, ci6, si2, si4, si5;
    INT4 idx, len;

    /* Sanity check input time series */
    LAL_CHECK_VALID_SERIES(V, XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(Phi, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, XLAL_FAILURE);

    /* Allocate polarization vectors and set to 0 */
    *hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &V->epoch, 0.0, 
            V->deltaT, &lalStrainUnit, V->data->length );
    *hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &V->epoch, 0.0, 
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( ! hplus || ! hcross )
        XLAL_ERROR(__func__, XLAL_EFUNC);
    memset((*hplus)->data->data, 0, (*hplus)->data->length 
            * sizeof(*(*hplus)->data->data));
    memset((*hcross)->data->data, 0, (*hcross)->data->length 
            * sizeof(*(*hcross)->data->data));

    M = m1 + m2;
    eta = m1 * m2 / M / M; // symmetric mass ratio - '\nu' in the paper
    eta2 = eta*eta;
    dm = abs(m1 - m2) / M; // frac. mass difference - \delta m/m in the paper
    dist = r / LAL_C_SI;   // r (m) / c (m/s) --> dist in units of seconds
    /* convert mass from kg to s, so ampfac ~ M/dist is dimensionless */
    ampfac = 2. * M * LAL_G_SI * pow(LAL_C_SI, -3) * eta / dist;
    
    /** 
     * cosines and sines of inclination between 
     * line of sight (N) and binary orbital angular momentum (L_N)
     */
    ci = cos(i);  	si = sin(i);
    ci2 = ci*ci;  	ci4 = ci2*ci2;  	ci6 = ci2*ci4;
    si2 = si*si;  	si4 = si2*si2;  	si5 = si*si4;

    /* loop over time steps and compute polarizations h+ and hx */
    len = V->data->length;
    for(idx = 0; idx < len; idx++)
    {
        /* Abbreviated names in lower case for time series at this sample */
        phi = Phi->data->data[idx]; 	v = V->data->data[idx];   
        v2 = v * v; 	v3 = v * v2; 	v4 = v2*v2;	v5 = v * v4;

        /** 
         * As explained in Arun et al, a phase shift can be applied 
         * to make log terms vanish which would appear in the amplitude 
         * at 1.5PN and 2.5PN orders. This shift is given in Eq. (5.6)
         * We apply the shift only for the PN orders which need it.
         */
        if( (ampO == -1) || ampO >= 5 )
            phiShift = 2.*v3*(1. - v2*eta/2.)*log( pow(v2 / x0 , 2./3.) );
        else if( ampO >= 3 )
            phiShift = 2.*v3*log( pow(v2 / x0 , 2./3.) );
        else
            phiShift = 0.;

        phi = phi - phiShift;

        /** 
         * First set all h+/x coefficients to 0. Then use a switch to
         * set proper non-zero values up to order ampO. Note we
         * fall through the PN orders and break only after Newt. order
         */
        hp0 = hp05 = hp1 = hp15 = hp2 = hp25 = 0.;
        hc0 = hc05 = hc1 = hc15 = hc2 = hc25 = 0.;

        switch( ampO )
        {
            /* case LAL_PNORDER_THREE_POINT_FIVE: */
            case 7:
            /* case LAL_PNORDER_THREE: */
            case 6:
                XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %s\n", __func__, ampO );
                XLAL_ERROR(__func__, XLAL_EINVAL);
                break;
            case -1: // Highest known PN order - move if higher terms added!
            /* case LAL_PNORDER_TWO_POINT_FIVE: */
            case 5:
                hp25 = cos(phi)*si*dm*(1771./5120. - ci2*1667./5120. 
                        + ci4*217./9216. - ci6/9126. + eta*(681./256. 
                        + ci2*13./768. - ci4*35./768. + ci6/2304.)
                        + eta2*(-3451./9216. + ci2*673./3072. - ci4*5./9216.
                        - ci6/3072.)) + cos(2.*phi)*LAL_PI*(19./3. + 3.*ci2 
                        - ci4*2./3. + eta*(-16./3. + ci2*14./3. + 2.*ci4))
                        + cos(3.*phi)*si*dm*(3537./1024. - ci2*22977./5120. 
                        - ci4*15309./5120. + ci6*729./5120. 
                        + eta*(-23829./1280. + ci2*5529./1280. 
                        + ci4*7749./1280. - ci6*729./1280.) 
                        + eta2*(29127./5120. - ci2*27267./5120. 
                        - ci4*1647./5120. + ci6*2187./5120.)) + cos(4.*phi)
                        * (-16.*LAL_PI/3.*(1. + ci2)*si2*(1. - 3.*eta))
                        + cos(5.*phi)*si*dm*(108125./9216. + ci2*40625./9216. 
                        + ci4*83125./9216. - ci6*15625./9216. 
                        + eta*(8125./256. - ci2*40625./2304. - ci4*48125./2304.
                        + ci6*15625./2304.) + eta2*(-119375./9216. 
                        + ci2*40625./3072. + ci4*44375./9216. 
                        - ci6*15625./3072.)) + cos(7.*phi)*dm
                        * (117649./46080.*si5*(1. + ci2)*(1. - 4.*eta 
                        + 3.*eta2)) + sin(2.*phi)*(-9./5. + ci2*14./5. 
                        + ci4*7./5. + eta*(96./5. - ci2*8./5. - ci4*28./5.)) 
                        + sin(4.*phi)*si2*(1. + ci2)*(56./5. - 32.*log(2.)/3. 
                        - eta*(1193./30. - 32.*log(2.)));
                hc25 = cos(2.*phi)*ci*(2. - ci2*22./5. + eta*(-154./5. 
                        + ci2*94./5.)) + cos(4.*phi)*ci*si2*(-112./5. 
                        + 64.*log(2.)/3. + eta*(1193./15. - 64.*log(2.)))
                        + sin(phi)*si*ci*dm*(-913./7680. + ci2*1891./11520. 
                        - ci4*7./4608. + eta*(1165./384. - ci2*235./576. 
                        + ci4*7./1152.) + eta2*(-1301./4608. + ci2*301./2304.
                        - ci4*7./1536.)) + sin(2.*phi)*LAL_PI*ci*(34./3. 
                        - ci2*8./3. - eta*(20./3. - 8.*ci2)) 
                        + sin(3.*phi)*si*ci*dm*(12501./2560. - ci2*12069./1280.
                        + ci4*1701./2560. + eta*(-19581./640. + ci2*7821./320.
                        - ci4*1701./640.) + eta2*(18903./2560. 
                        - ci2*11403./1280. + ci4*5103./2560.)) 
                        + sin(4.*phi)*si2*ci*(-32.*LAL_PI/3.*(1. - 3.*eta))
                        + sin(5.*phi)*si*ci*dm*(-101875./4608. + ci2*6875./256.
                        - ci4*21875./4608. + eta*(66875./1152. 
                        - ci2*44375./576. + ci4*21875./1152.) 
                        + eta2*(-100625./4608. + ci2*83125./2304. 
                        - ci4*21875./1536.)) + sin(7.*phi)*si5*ci*dm
                        * (117649./23040.*(1. - 4.*eta + 3.*eta2));
            /* case LAL_PNORDER_TWO: */
            case 4:
                hp2 = cos(phi)*LAL_PI*si*dm*(-5./8. - ci2*5./8.) 
                        + cos(2.*phi)*(11./60. + ci2*33./10. + ci4*29./24. 
                        + ci6/24. + eta*(353./36. - 3.*ci2 - ci4*251./72. 
                        + ci6*5./24.) + eta2*(-49./12. + ci2*9./2. 
                        - ci4*7./24. - ci6*5./24.)) + cos(3.*phi)*LAL_PI*si*dm
                        * (27./8.*(1 + ci2)) + cos(4.*phi)*(118./15. 
                        - ci2*16./5. - ci4*86./15. + ci6*16./15. 
                        + eta*(-262./9. + 16.*ci2 + ci4*166./9. - ci6*16./3.)
                        + eta2*(14. - 16.*ci2 - ci4*10./3. + ci6*16./3.))
                        + cos(6.*phi)*(-81./40.*si4*(1. + ci2)
                        * (1. - 5.*eta + 5.*eta2)) + sin(phi)*si*dm
                        * (11./40. + 5.*log(2)/4. + ci2*(7./40. + log(2)/4.))
                        + sin(3.*phi)*si*dm*((-189./40. + 27./4.*log(3./2.))
                        * (1. + ci2));
                hc2 = cos(phi)*si*ci*dm*(-9./20. - 3./2.*log(2.)) 
                        + cos(3.*phi)*si*ci*dm*(189./20. - 27./2.*log(3./2.))
                        - sin(phi)*si*ci*dm*3.*LAL_PI/4. 
                        + sin(2.*phi)*ci*(17./15. + ci2*113./30. - ci4/4.
                        + eta*(143./9. - ci2*245./18. + ci4*5./4.)
                        + eta2*(-14./3. + ci2*35./6. - ci4*5./4.))
                        + sin(3.*phi)*si*ci*dm*27.*LAL_PI/4.
                        + sin(4.*phi)*ci*(44./3. - ci2*268./15. + ci4*16./5. 
                        + eta*(-476./9. + ci2*620./9. - 16.*ci4)
                        + eta2*(68./3. - ci2*116./3. + 16.*ci4))
                        + sin(6.*phi)*ci*(-81./20.*si4
                        * (1. - 5.*eta + 5.*eta2));
            /* case LAL_PNORDER_ONE_POINT_FIVE: */
            case 3:
                hp15 = cos(phi)*si*dm*(19./64. + ci2*5./16. - ci4/192. 
                        + eta*(-49./96. + ci2/8. + ci4/96.))
                        + cos(2.*phi)*(-2.*LAL_PI*(1. + ci2))
                        + cos(3.*phi)*si*dm*(-657./128. - ci2*45./16. 
                        + ci4*81./128. + eta*(225./64. - ci2*9./8. 
                        - ci4*81./64.)) + cos(5.*phi)*si*dm*(625./384.*si2
                        * (1. + ci2)*(1. - 2.*eta));
                hc15 = sin(phi)*si*ci*dm*(21./32. - ci2*5./96. 
                        + eta*(-23./48. + ci2*5./48.)) 
                        - 4.*LAL_PI*ci*sin(2.*phi) + sin(3.*phi) * si * ci * dm
                        * (-603./64. + ci2*135./64. 
                        + eta*(171./32. - ci2*135./32.)) 
                        + sin(5.*phi)*si*ci*dm*(625./192.*si2*(1. - 2.*eta));
            /* case LAL_PNORDER_ONE: */
            case 2:
                hp1 = cos(2.*phi)*(19./6. + 3./2.*ci2 - ci4/3. 
                        + eta*(-19./6. + ci2*11./6. + ci4)) 
                        - cos(4.*phi) * (4./3.*si2*(1. + ci2)*(1. - 3.*eta));
                hc1 = sin(2.*phi)*ci*(17./3. - ci2*4./3. 
                        + eta*(-13./3. + 4.*ci2)) 
                        + sin(4.*phi) * ci * si2*(-8./3.*(1. - 3.*eta));
            /*case LAL_PNORDER_HALF:*/
            case 1:
                hp05 = - si*dm*(cos(phi)*(5./8. + ci2/8.) 
                        - cos(3.*phi)*(9./8. + 9.*ci2/8.));
                hc05 = si*ci*dm*(-sin(phi)*3./4. + sin(3.*phi)*9./4.);
            case 0:
                hp0 = -(1. + ci2)*cos(2.*phi);
                hc0 = -2.*ci*sin(2.*phi);
                break;
            /*case LAL_PNORDER_NEWTONIAN:*/
            default:
                XLALPrintError("XLAL Error - %s: Invalid amp. PN order %s\n",
                        __func__, ampO );
                XLAL_ERROR(__func__, XLAL_EINVAL);
                break;
        } /* End switch on ampO */

        /* Fill the output polarization arrays */
        (*hplus)->data->data[idx] = ampfac * v2 * ( hp0 + v * ( hp05 
                + v * ( hp1 + v * ( hp15 + v * ( hp2 + v * ( hp25 ) ) ) ) ) );
        (*hcross)->data->data[idx] = ampfac * v2 * ( hc0 + v * ( hc05 
                + v * ( hc1 + v * ( hc15 + v * ( hc2 + v * ( hc25 ) ) ) ) ) );

    } /* end loop over time series samples idx */
    return XLAL_SUCCESS;
}
