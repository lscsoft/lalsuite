/*
 * Copyright (C) 2008 J. Creighton, Evan Ochsner, Chris Pankow, 2018 Riccardo Sturani
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
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define EPS 1.e-6

/**
 * @addtogroup LALSimInspiralPNMode_c
 * @brief Routines for mode decompositions of post-Newtonian waveforms
 *
 * @{
 */

/**
 * Computes h(l,m) mode timeseries of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * See Eqns. (79)-(116) of:
 * Lawrence E. Kidder, \"Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit\", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(
		REAL8TimeSeries *v,   /**< post-Newtonian parameter */
	       	REAL8TimeSeries *phi, /**< orbital phase */
	       	REAL8 v0,             /**< tail gauge parameter (default = 1) */
	       	REAL8 m1,             /**< mass of companion 1 (kg) */
	       	REAL8 m2,             /**< mass of companion 2 (kg) */
	       	REAL8 r,              /**< distance of source (m) */
	       	int O,                /**< twice post-Newtonain order */
	       	int l,                /**< mode number l */
	       	int m                 /**< mode number m */
		)
{
	LAL_CHECK_VALID_SERIES(v, NULL);
	LAL_CHECK_VALID_SERIES(phi, NULL);
	LAL_CHECK_CONSISTENT_TIME_SERIES(v, phi, NULL);
	COMPLEX16TimeSeries *hlm;
	UINT4 j;
	if ( l == 2 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode22(v, phi, v0, m1, m2, r, O);
	else if ( l == 2 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode21(v, phi, v0, m1, m2, r, O);
	else if ( l == 2 && m == 0 )
		hlm = XLALSimInspiralPNMode20(v, phi, v0, m1, m2, r, O);
	else if ( l == 3 && abs(m) == 3 )
		hlm = XLALSimInspiralPNMode33(v, phi, v0, m1, m2, r, O);
	else if ( l == 3 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode32(v, phi, v0, m1, m2, r, O);
	else if ( l == 3 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode31(v, phi, v0, m1, m2, r, O);
	else if ( l == 3 && m == 0 )
		hlm = XLALSimInspiralPNMode30(v, phi, v0, m1, m2, r, O);
	else if ( l == 4 && abs(m) == 4 )
		hlm = XLALSimInspiralPNMode44(v, phi, v0, m1, m2, r, O);
	else if ( l == 4 && abs(m) == 3 )
		hlm = XLALSimInspiralPNMode43(v, phi, v0, m1, m2, r, O);
	else if ( l == 4 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode42(v, phi, v0, m1, m2, r, O);
	else if ( l == 4 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode41(v, phi, v0, m1, m2, r, O);
	else if ( l == 4 && m == 0 )
		hlm = XLALSimInspiralPNMode40(v, phi, v0, m1, m2, r, O);
	else if ( l == 5 && abs(m) == 5 )
		hlm = XLALSimInspiralPNMode55(v, phi, v0, m1, m2, r, O);
	else if ( l == 5 && abs(m) == 4 )
		hlm = XLALSimInspiralPNMode54(v, phi, v0, m1, m2, r, O);
	else if ( l == 5 && abs(m) == 3 )
		hlm = XLALSimInspiralPNMode53(v, phi, v0, m1, m2, r, O);
	else if ( l == 5 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode52(v, phi, v0, m1, m2, r, O);
	else if ( l == 5 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode51(v, phi, v0, m1, m2, r, O);
	else if ( l == 5 && m == 0 )
		hlm = XLALSimInspiralPNMode50(v, phi, v0, m1, m2, r, O);
	else if ( l == 6 && abs(m) == 6 )
		hlm = XLALSimInspiralPNMode66(v, phi, v0, m1, m2, r, O);
	else if ( l == 6 && abs(m) == 5 )
		hlm = XLALSimInspiralPNMode65(v, phi, v0, m1, m2, r, O);
	else if ( l == 6 && abs(m) == 4 )
		hlm = XLALSimInspiralPNMode64(v, phi, v0, m1, m2, r, O);
	else if ( l == 6 && abs(m) == 3 )
		hlm = XLALSimInspiralPNMode63(v, phi, v0, m1, m2, r, O);
	else if ( l == 6 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode62(v, phi, v0, m1, m2, r, O);
	else if ( l == 6 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode61(v, phi, v0, m1, m2, r, O);
	else if ( l == 6 && m == 0 )
		hlm = XLALSimInspiralPNMode60(v, phi, v0, m1, m2, r, O);
	else {
		XLALPrintError("XLAL Error - %s: Unsupported mode l=%d, m=%d\n", __func__, l, m );
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}
	if ( !hlm )
		XLAL_ERROR_NULL(XLAL_EFUNC);
	if ( m < 0 ) {
		REAL8 sign = l % 2 ? -1.0 : 1.0;
		for ( j = 0; j < hlm->data->length; ++j )
			hlm->data->data[j] = sign * conj(hlm->data->data[j]);
	}
	return hlm;
}

/**
 * Computes h(l,m) mode timeseries of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * See Eqns. (79)-(116) of:
 * Lawrence E. Kidder, \"Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit\", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALCreateSimInspiralPNModeCOMPLEX16TimeSeriesLALConvention(
       REAL8TimeSeries *v,   /**< post-Newtonian parameter */
       REAL8TimeSeries *phi, /**< orbital phase */
       REAL8 m1,             /**< mass of companion 1 (kg) */
       REAL8 m2,             /**< mass of companion 2 (kg) */
       REAL8 r,              /**< distance of source (m) */
       int O,                /**< twice post-Newtonain order */
       int l,                /**< mode number l */
       int m                 /**< mode number m */
										 )
{
	LAL_CHECK_VALID_SERIES(v, NULL);
	LAL_CHECK_VALID_SERIES(phi, NULL);
	LAL_CHECK_CONSISTENT_TIME_SERIES(v, phi, NULL);
	COMPLEX16TimeSeries *hlm;
	UINT4 j;
	if ( l == 2 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode22(v, phi, 1., m1, m2, r, O);
	else if ( l == 2 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode21(v, phi, 1., m1, m2, r, O);
	else if ( l == 2 && m == 0 )
		hlm = XLALSimInspiralPNMode20(v, phi, 1., m1, m2, r, O);
	else if ( l == 3 && abs(m) == 3 )
		hlm = XLALSimInspiralPNMode33(v, phi, 1., m1, m2, r, O);
	else if ( l == 3 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode32(v, phi, 1., m1, m2, r, O);
	else if ( l == 3 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode31(v, phi, 1., m1, m2, r, O);
	else if ( l == 3 && m == 0 )
		hlm = XLALSimInspiralPNMode30(v, phi, 1., m1, m2, r, O);
	else if ( l == 4 && abs(m) == 4 )
		hlm = XLALSimInspiralPNMode44(v, phi, 1., m1, m2, r, O);
	else if ( l == 4 && abs(m) == 3 )
		hlm = XLALSimInspiralPNMode43(v, phi, 1., m1, m2, r, O);
	else if ( l == 4 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode42(v, phi, 1., m1, m2, r, O);
	else if ( l == 4 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode41(v, phi, 1., m1, m2, r, O);
	else if ( l == 4 && m == 0 )
		hlm = XLALSimInspiralPNMode40(v, phi, 1., m1, m2, r, O);
	else if ( l == 5 && abs(m) == 5 )
		hlm = XLALSimInspiralPNMode55(v, phi, 1., m1, m2, r, O);
	else if ( l == 5 && abs(m) == 4 )
		hlm = XLALSimInspiralPNMode54(v, phi, 1., m1, m2, r, O);
	else if ( l == 5 && abs(m) == 3 )
		hlm = XLALSimInspiralPNMode53(v, phi, 1., m1, m2, r, O);
	else if ( l == 5 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode52(v, phi, 1., m1, m2, r, O);
	else if ( l == 5 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode51(v, phi, 1., m1, m2, r, O);
	else if ( l == 5 && m == 0 )
		hlm = XLALSimInspiralPNMode50(v, phi, 1., m1, m2, r, O);
	else if ( l == 6 && abs(m) == 6 )
		hlm = XLALSimInspiralPNMode66(v, phi, 1., m1, m2, r, O);
	else if ( l == 6 && abs(m) == 5 )
		hlm = XLALSimInspiralPNMode65(v, phi, 1., m1, m2, r, O);
	else if ( l == 6 && abs(m) == 4 )
		hlm = XLALSimInspiralPNMode64(v, phi, 1., m1, m2, r, O);
	else if ( l == 6 && abs(m) == 3 )
		hlm = XLALSimInspiralPNMode63(v, phi, 1., m1, m2, r, O);
	else if ( l == 6 && abs(m) == 2 )
		hlm = XLALSimInspiralPNMode62(v, phi, 1., m1, m2, r, O);
	else if ( l == 6 && abs(m) == 1 )
		hlm = XLALSimInspiralPNMode61(v, phi, 1., m1, m2, r, O);
	else if ( l == 6 && m == 0 )
		hlm = XLALSimInspiralPNMode60(v, phi, 1., m1, m2, r, O);
	else {
		XLALPrintError("XLAL Error - %s: Unsupported mode l=%d, m=%d\n", __func__, l, m );
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}
	if ( !hlm )
		XLAL_ERROR_NULL(XLAL_EFUNC);
	REAL8 sign=-1;
	if ( m < 0 ) {
		sign*= l % 2 ? -1.0 : 1.0;
		for ( j = 0; j < hlm->data->length; ++j )
		  hlm->data->data[j] = sign * conj(hlm->data->data[j]);
	}
	else
	        for ( j = 0; j < hlm->data->length; ++j )
		  hlm->data->data[j] = -hlm->data->data[j];

	return hlm;
}

/**
 * Computes h(2,2) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (79) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode22(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 v0,             /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_22 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = -8.0*sqrt(LAL_PI/5.0)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 nu3 = nu*nu2;
    REAL8 pi2 = LAL_PI*LAL_PI;
    REAL8 phi, v, v2, logv, logv_v0, logv0 = log(v0);
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re0 = 1., re2 = 0., re3 = 0., im3log = 0., re4 = 0., re5 = 0.;
    REAL8 im5 = 0., im5log = 0., re6 = 0., im6 = 0., re6log = 0.;
    REAL8 re6logsq = 0., im6log = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
            re6 = (27027409.0/646800.0) - (856.0/105.0)*LAL_GAMMA 
                    + (2.0/3.0)*pi2 - (1712.0/105.0)*log(2.0)
                    - ((278185.0/33264.0) - (41.0/96.0)*pi2)*nu
                    - (20261.0/2772.0)*nu2 + (114635.0/99792.0)*nu3;
            re6log = - (856.0/105.0);
            re6logsq = - 72.0;
            im6 = (428.0/105.0)*LAL_PI;
            im6log = 24.0*LAL_PI;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 5:
            re5 = - ((107.0/21.0) - (34.0/21.0)*nu)*LAL_PI;
            im5 = - 24.0*nu;
            im5log = - ((107.0/7.0) - (34.0/7.0)*nu)*2.0;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
            re4 = - ((2173.0/1512.0) + (1069.0/216.0)*nu 
                    - (2047.0/1512.0)*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
            re3 = 2.0*LAL_PI;
            im3log = 12.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 2:
            re2 = - ((107.0/42.0) - (55.0/42.0)*nu);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        logv = log(v);
        logv_v0 = logv - logv0;
        phi = Phi->data->data[j];
        re = re0 + v2*(re2 + v*(re3 + v*(re4 + v*(re5 + v*(re6
                + re6log*logv + re6logsq*logv_v0*logv_v0)))));
        im = v*v2*(im3log*logv_v0 + v2*(im5 + im5log*logv_v0
                + v*(im6 + im6log*logv_v0)));
        ans = cpolar(1.0, -2.0*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*m/r)*v2);
    }
    return hlm;
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
COMPLEX16TimeSeries *XLALSimInspiralPNMode21(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 v0,             /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_21 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = -(8.0/3.0)*sqrt(LAL_PI/5.0)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re1 = 0., re3 = 0., re4 = 0., im4 = 0., im4log = 0., re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = -((43.0/126.0)+(509.0/126.0)*nu-(79.0/168.0)*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
            re4 = LAL_PI;
            im4 = -(1.0/2.0) - 2.0*log(2.0);
            im4log = 6.0;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
            re3 = - ((17.0/28.0) - (5.0/7.0)*nu);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
            re1 = 1.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re1 + v2 * (re3 + v * (re4 + v * re5));
        im = v*v2 * (im4 + im4log * log(v/v0));
        ans = cpolar(1.0, -phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v2*v);
    }
    return hlm;
}

/**
 * Computes h(2,0) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (81) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode20(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int UNUSED O          /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_20 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = (2./7.)*sqrt(10.*LAL_PI/3.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 v, v2;
    COMPLEX16 ans;
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        ans = 1.0;
        hlm->data->data[j] = ans * ((fac*mu/r)*v2);
    }
    return hlm;
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
COMPLEX16TimeSeries *XLALSimInspiralPNMode33(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 v0,             /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_33 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = 3.0*sqrt(6.0*LAL_PI/7.0)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re1 = 0., re3 = 0., re4 = 0., im4 = 0., im4log = 0., re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = ((123.0/110.0) - (1838.0/165.0)*nu
                    - (887.0/330.0)*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
            re4 = 3.0*LAL_PI;
            im4 = -(21.0/5.0) + 6.0*log(3.0/2.0);
            im4log = 18.0;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
            re3 = - (4.0 - 2.0*nu);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
            re1 = 1.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re1 + v2 * (re3 + v * (re4 + v * re5));
        im = v*v2 * (im4 + im4log * log(v/v0));
        ans = cpolar(1.0, -3.0*phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v2*v);
    }
    return hlm;
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
COMPLEX16TimeSeries *XLALSimInspiralPNMode32(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 v0,             /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_32 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = -(8.0/3.0)*sqrt(LAL_PI/7.0)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 nuterm = (1. - 3.*nu);
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re2 = 0., re4 = 0., re5 = 0., im5 = 0., im5log = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = 2.0*LAL_PI*nuterm;
            im5 = -3.0 + (66.0/5.0)*nu;
            im5log = 12.0*nuterm;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
            re4 = - ((193.0/90.0) - (145.0/18.0)*nu 
                    + (73.0/18.0)*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
        case 2:
            re2 = nuterm;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re2 + v2 * (re4 + v * re5);
        im = v*v2 * (im5 + im5log * log(v/v0));
        ans = cpolar(1.0, -2.0*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*m/r)*v2*v2);
    }
    return hlm;
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
COMPLEX16TimeSeries *XLALSimInspiralPNMode31(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 v0,             /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_31 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = -(1.0/3.0)*sqrt(2.0*LAL_PI/35.0)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re1 = 0., re3 = 0., re4 = 0., im4 = 0., im4log = 0., re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = ((607.0/198.0) - (136.0/99.0)*nu 
                    - (247.0/198.0)*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
            re4 = LAL_PI;
            im4 = -(7.0/5.0) - 2.0*log(2.0);
            im4log = 6.0;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
            re3 = - ((8.0/3.0) + (2.0/3.0)*nu);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
            re1 = 1.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re1 + v2 * (re3 + v * (re4 + v * re5));
        im = v*v2 * (im4 + im4log * log(v/v0));
        ans = cpolar(1.0, -phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v2*v);
    }
    return hlm;
}

/**
 * Computes h(3,0) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (85) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode30(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_30 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = (16./5.)*sqrt(6.*LAL_PI/35.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 v, v7;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re = 1.;
        case 4:
        case 3:
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v7 = v*v*v*v*v*v*v;
        ans = crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*m*nu2/r)*v7);
    }
    return hlm;
}

/**
 * Computes h(4,4) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (86) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode44(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 v0,             /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_44 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = (64./9.)*sqrt(LAL_PI/7.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 nuterm = 1. - 3.*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re2 = 0., re4 = 0., re5 = 0., im5 = 0., im5log = 0., re6 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
            re6 = 1068671./200200. - (1088119./28600.)*nu
                    + (146879./2340.)*nu2 - (226097./17160.)*nu2*nu;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 5:
            re5 = 4.*LAL_PI*nuterm;
            im5 = -42./5. + (1193./40.) *nu + 8.*nuterm*log(2.);
            im5log = 24.*nuterm;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
            re4 = 593./110. - (1273./66.)*nu + (175./22.)*nu2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
        case 2:
            re2 = nuterm;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re2 + v2 * (re4 + v * (re5 + v * re6));
        im = v*v2 * (im5 + im5log * log(v/v0));
        ans = cpolar(1.0, -4.*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*m/r)*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(4,3) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (87) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode43(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_43 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = (9./5.)*sqrt(2.*LAL_PI/7.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re3 = 0., re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = 39./11. - (1267./132.)*nu + (131./33.)*nu2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re3 + v2 * re5;
        ans = cpolar(1.0, -3.*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*dm/r)*v*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(4,2) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (88) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode42(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 v0,             /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_42 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = - (8./63.)*sqrt(LAL_PI)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 nuterm = 1. - 3.*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re2 = 0., re4 = 0., re5 = 0., im5 = 0., im5log = 0., re6 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
            re6 = 1038039./200200. - (606751./28600.)*nu
                    + (400453./25740.)*nu2 + (25783./17160.)*nu*nu2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 5:
            re5 = 2.*LAL_PI*nuterm;
            im5 = -21./5. + (84./5.) *nu + 8.*nuterm*log(2.);
            im5log = 12.*nuterm;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
            re4 = 437./110. - (805./66.)*nu + (19./22.)*nu2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
        case 2:
            re2 = nuterm;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re2 + v2 * (re4 + v * (re5 + v * re6));
        im = v*v2 * (im5 + im5log * log(v/v0));
        ans = cpolar(1.0, -2.*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*m/r)*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(4,1) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (89) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode41(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_41 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = - (1./105.)*sqrt(2.*LAL_PI)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re3 = 0., re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = - (101./33. - (337./44.)*nu + (83./33.)*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re3 + v2 * re5;
        ans = cpolar(1.0, -phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(4,0) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (90) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode40(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int UNUSED O          /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_40 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = (1./63.)*sqrt(LAL_PI/10.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 v, v2;
    COMPLEX16 ans;
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        ans = 1.0;
        hlm->data->data[j] = ans * ((fac*mu/r)*v2);
    }
    return hlm;
}

/**
 * Computes h(5,5) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (91) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode55(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_55 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = - (125./12.)*sqrt(5.*LAL_PI/66.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re3 = 0., re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = - (263./39. - (688./39.)*nu + (256./39.)*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re3 + v2 * re5;
        ans = cpolar(1.0, -5.*phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(5,4) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (92) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode54(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_54 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = (256./45.)*sqrt(LAL_PI/33.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re4 = 0., re6 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
            re6 = - (4451./910. - (3619./130.)*nu + (521./13.)*nu2
                    - (339./26.)*nu*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re4 + v2 * re6;
        ans = cpolar(1.0, -4.*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*m/r)*v2*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(5,3) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (93) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode53(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_53 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = - (9./20.)*sqrt(2.*LAL_PI/22.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re3 = 0., re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = - (69./13. - (464./39.)*nu + (88./39.)*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re3 + v2 * re5;
        ans = cpolar(1.0, -3.*phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(5,2) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (94) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode52(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_52 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = - (16./135.)*sqrt(LAL_PI/11.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re4 = 0., re6 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
            re6 = - (3911./910. - (3079./130.)*nu + (413./13.)*nu2
                    - (231./26.)*nu*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re4 + v2 * re6;
        ans = cpolar(1.0, -2.*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*m/r)*v2*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(5,1) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (95) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode51(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_51 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = - (1./180.)*sqrt(LAL_PI/77.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re3 = 0., re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = - (179./39. - (352./39.)*nu + (4./39.)*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re3 + v2 * re5;
        ans = cpolar(1.0, -1.*phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(5,0) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * THIS MODE IS ZERO TO THE ORDER CONSIDERED IN:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode50(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 UNUSED m1,      /**< mass of companion 1 (kg) */
        REAL8 UNUSED m2,      /**< mass of companion 2 (kg) */
        REAL8 UNUSED r,       /**< distance of source (m) */
        int UNUSED O          /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_50 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        hlm->data->data[j] = 0.0;
    }
    return hlm;
}

/**
 * Computes h(6,6) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (96) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode66(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_66 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = - (432./5.)*sqrt(LAL_PI/715.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re4 = 0., re6 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
            re6 = - (113./14. - (91./2.)*nu + (64.)*nu2
                    - (39./2.)*nu*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re4 + v2 * re6;
        ans = cpolar(1.0, -6.*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*m/r)*v2*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(6,5) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (97) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode65(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_65 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = -(625./63.)*sqrt(5.*LAL_PI/429.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = 1. - 4.*nu + 3.*nu2;
        case 4:
        case 3:
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re5;
        ans = cpolar(1.0, -5.*phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v*v2*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(6,4) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (98) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode64(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_64 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac =(1024./495.)*sqrt(2.*LAL_PI/195.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re4 = 0., re6 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
            re6 = - (113./14. - (91./2.)*nu + (64.)*nu2
                    - (39./2.)*nu*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re4 + v2 * re6;
        ans = cpolar(1.0, -4.*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*m/r)*v2*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(6,3) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (99) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode63(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_63 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = (81./385.)*sqrt(LAL_PI/13.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = 1. - 4.*nu + 3.*nu2;
        case 4:
        case 3:
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re5;
        ans = cpolar(1.0, -3.*phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v*v2*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(6,2) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (100) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode62(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_62 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = - (16./1485.)*sqrt(LAL_PI/13.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re4 = 0., re6 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
            re6 = - (81./14. - (59./2.)*nu + 32.*nu2
                    - (7./2.)*nu*nu2);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case 3:
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re4 + v2 * re6;
        ans = cpolar(1.0, -2.*phi) * crect(re, im);
        hlm->data->data[j] = ans * ((fac*nu*m/r)*v2*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(6,1) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * Implements Equation (101) of:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode61(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 m1,             /**< mass of companion 1 (kg) */
        REAL8 m2,             /**< mass of companion 2 (kg) */
        REAL8 r,              /**< distance of source (m) */
        int O                 /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_61 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    REAL8 fac = - (1./2079.)*sqrt(2.*LAL_PI/65.)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
    REAL8 m = m1 + m2;
    REAL8 dm = m1 - m2;
    REAL8 mu = m1*m2/m;
    REAL8 nu = mu/m;
    REAL8 nu2 = nu*nu;
    REAL8 phi, v, v2;
    COMPLEX16 ans;
    REAL8 re = 0.0;
    REAL8 im = 0.0;
    REAL8 re5 = 0.;
    /* Set coefficients for requested PN (amplitude) order */
    switch (O) {
        default: /* unsupported pN order */
            XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, O/2, O%2?".5":"" );
            XLAL_ERROR_NULL(XLAL_EINVAL);
        case -1: /* use highest available pN order */
        case 6:
        case 5:
            re5 = 1. - 4.*nu + 3.*nu2;
        case 4:
        case 3:
        case 2:
        case 1:
        case 0:
            break;
    }
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        v = V->data->data[j];
        v2 = v*v;
        phi = Phi->data->data[j];
        re = re5;
        ans = cpolar(1.0, -phi) * crect(re, im);
        hlm->data->data[j] = ans * I * ((fac*nu*dm/r)*v*v2*v2*v2);
    }
    return hlm;
}

/**
 * Computes h(6,0) mode of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform.
 *
 * THIS MODE IS ZERO TO THE ORDER CONSIDERED IN:
 * Lawrence E. Kidder, "Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
COMPLEX16TimeSeries *XLALSimInspiralPNMode60(
        REAL8TimeSeries *V,   /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi, /**< orbital phase */
        REAL8 UNUSED v0,      /**< tail gauge parameter (default=1) */
        REAL8 UNUSED m1,      /**< mass of companion 1 (kg) */
        REAL8 UNUSED m2,      /**< mass of companion 2 (kg) */
        REAL8 UNUSED r,       /**< distance of source (m) */
        int UNUSED O          /**< twice post-Newtonian order */
        )
{
    LAL_CHECK_VALID_SERIES(V, NULL);
    LAL_CHECK_VALID_SERIES(Phi, NULL);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, NULL);
    COMPLEX16TimeSeries *hlm;
    hlm = XLALCreateCOMPLEX16TimeSeries( "H_60 MODE", &V->epoch, 0.0,
            V->deltaT, &lalStrainUnit, V->data->length );
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);
    UINT4 j;
    /* Loop over time samples, compute hlm(t) */
    for (j=0; j < V->data->length; j++) {
        hlm->data->data[j] = 1.0;
    }
    return hlm;
}

/**
 * Utility type and function to compute trigonometric functions from the cosine of an angle
 */

typedef struct tagLALSimInspiralInclAngle {
  REAL8 ci;
  REAL8 si;
  REAL8 ciSq;
  REAL8 siSq;
  REAL8 c2i;
  REAL8 s2i;
  REAL8 c3i;
  REAL8 s3i;
  REAL8 ciBy2;
  REAL8 siBy2;
  REAL8 ciBy2Sq;
  REAL8 siBy2Sq;
  REAL8 ciBy2Qu;
  REAL8 siBy2Qu;
  REAL8 ciBy2Sx;
  REAL8 siBy2Sx;
  REAL8 ciBy2Et;
  REAL8 siBy2Et;
} LALSimInspiralInclAngle;

static
LALSimInspiralInclAngle *XLALSimInspiralComputeInclAngle(
        REAL8 ciota /** INPUT */
        )
{
  LALSimInspiralInclAngle *angle=(LALSimInspiralInclAngle *) LALMalloc(sizeof(LALSimInspiralInclAngle));

  angle->ci=ciota;
  angle->ciSq=ciota*ciota;
  angle->siSq=1.-ciota*ciota;
  angle->si=sqrt(angle->siSq);
  angle->c2i=angle->ciSq-angle->siSq;
  angle->s2i=2.*ciota*angle->si;
  angle->c3i=angle->c2i*ciota-angle->s2i*angle->si;
  angle->s3i=ciota*angle->s2i+angle->si*angle->c2i;
  angle->ciBy2Sq=(1.+ciota)/2.;
  angle->siBy2Sq=(1.-ciota)/2.;
  angle->ciBy2=sqrt(angle->ciBy2Sq);
  angle->siBy2=sqrt(angle->siBy2Sq);
  angle->ciBy2Qu=angle->ciBy2Sq*angle->ciBy2Sq;
  angle->siBy2Qu=angle->siBy2Sq*angle->siBy2Sq;
  angle->ciBy2Sx=angle->ciBy2Qu*angle->ciBy2Sq;
  angle->siBy2Sx=angle->siBy2Qu*angle->siBy2Sq;
  angle->ciBy2Et=angle->ciBy2Qu*angle->ciBy2Qu;
  angle->siBy2Et=angle->siBy2Qu*angle->siBy2Qu;
  return angle;

} /* End of XLALSimInspiralComputeInclAngle*/

/**
 * Computes the 5 l=2 modes of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform for spinning-precessing waveforms.
 *
 * For reference see eq. B1,2 of
 * Arun et al., "Higher-order spin effects in the amplitude and phase of gravitational waveforms
 * emitted by inspiraling compact binaries: Ready-to-use gravitational waveforms",
 * Physical Review D 79, 104023 (2009), Erratum PRD 84, 049901 (2011), arXiv:0810.5336v3 [gr-qc].
 * All variable are inputs.
 *
 * !!BEWARE!! Spin components are given wrt LNh. LNh components define angles
 * incl and alpha (polarization angle) entering the mode analytic form.
 *
 */

INT4 XLALSimInspiralSpinPNMode2m(SphHarmTimeSeries **h2ms, /**< OUTPUT*/
				 REAL8TimeSeries *V,    /**< post-Newtonian parameter*/
				 REAL8TimeSeries *Phi,  /**< orbital phase */
				 REAL8TimeSeries *LNhx, /**< angular momentum unit vector x component */
				 REAL8TimeSeries *LNhy, /**< angular momentum unit vector y component */
				 REAL8TimeSeries *LNhz, /**< angular momentum unit vector z components */
				 REAL8TimeSeries *e1x,  /**< x-axis tetrad x component*/
				 REAL8TimeSeries *e1y,  /**< x-axis tetrad y component*/
				 REAL8TimeSeries *e1z,  /**< x-axis tetrad z component*/
				 REAL8TimeSeries *S1x,  /**< spin1-x component */
				 REAL8TimeSeries *S1y,  /**< spin1-y component */
				 REAL8TimeSeries *S1z,  /**< spin1-z component */
				 REAL8TimeSeries *S2x,  /**< spin2-x component */
				 REAL8TimeSeries *S2y,  /**< spin2-y component */
				 REAL8TimeSeries *S2z,  /**< spin2-z component */
				 REAL8 m1,              /**< mass of companion 1 (kg) */
				 REAL8 m2,              /**< mass of companion 2 (kg) */
				 REAL8 distance,        /**< distance of source (m) */
				 int ampO               /**< twice post-Newtonian amp-order */
				 )
{

  LAL_CHECK_VALID_SERIES(V, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(Phi, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(LNhx, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(LNhy, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(LNhz, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(e1x, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(e1y, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(e1z, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhx, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhy, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhz, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, e1x, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, e1y, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, e1z, XLAL_FAILURE);

  UINT4 idx;
  REAL8 eta=m1/(m1+m2)*m2/(m1+m2);
  REAL8 dm=(m1-m2)/(m1+m2);

  REAL8 amp1   =0.;
  REAL8 amp2   =0.;
  REAL8 amp2S  =0.;
  REAL8 amp3   =0.;
  REAL8 amp3pi =0.;
  REAL8 amp3S  =0.;

  switch (ampO) {
  case -1: /* use highest available pN order */
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
    __attribute__ ((fallthrough));
#endif
  case 3:
    amp3S  = 1.;
    amp3   = (-1.7+2.*eta)/2.8;
    amp3pi = 2.*LAL_PI;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
    __attribute__ ((fallthrough));
#endif
  case 2:
    amp2  = -10.7/4.2 +5.5/4.2*eta;
    amp2S = 1.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
    __attribute__ ((fallthrough));
#endif
  case 1:
    amp1 = 1.;
  case 0:
    break;
  default: /* unsupported pN order */
    XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, ampO/2, ampO%2?".5":"" );
    XLAL_ERROR(XLAL_EINVAL);
  }

  REAL8TimeSeries *S1xtmp=XLALCreateREAL8TimeSeries("S1x",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S1ytmp=XLALCreateREAL8TimeSeries("S1y",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S1ztmp=XLALCreateREAL8TimeSeries("S1z",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S2xtmp=XLALCreateREAL8TimeSeries("S2x",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S2ytmp=XLALCreateREAL8TimeSeries("S2y",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S2ztmp=XLALCreateREAL8TimeSeries("S2z",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  for (idx=0;idx<V->data->length;idx++) {
    S1xtmp->data->data[idx]=0.;
    S1ytmp->data->data[idx]=0.;
    S1ztmp->data->data[idx]=0.;
    S2xtmp->data->data[idx]=0.;
    S2ytmp->data->data[idx]=0.;
    S2ztmp->data->data[idx]=0.;
  }

  REAL8 checkL,checke;
  for (idx=0;idx<LNhx->data->length;idx++){
    checkL=LNhx->data->data[idx]*LNhx->data->data[idx]+LNhy->data->data[idx]*LNhy->data->data[idx]+LNhz->data->data[idx]*LNhz->data->data[idx];
    if ((checkL<1.-EPS)||(checkL>1.+EPS)) {
      XLALPrintError("XLAL Error - %s: LNhat unit vector has not unit norm %12.4e at idx %d\n", __func__,checkL,idx);
      XLAL_ERROR(XLAL_EINVAL);
    }
    checke=e1x->data->data[idx]*e1x->data->data[idx]+e1y->data->data[idx]*e1y->data->data[idx]+e1z->data->data[idx]*e1z->data->data[idx];
    if ((checke<1.-EPS)||(checke>1.+EPS)) {
      XLALPrintError("XLAL Error - %s: e1 unit vector has not unit norm %12.4e at idx %d\n", __func__,checke,idx);
      XLAL_ERROR(XLAL_EINVAL);
    }
  }

  if ( S1x ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1x, XLAL_FAILURE);
    memcpy(S1xtmp->data->data, S1x->data->data, S1xtmp->data->length*sizeof(REAL8) );
  }

  if ( S1y ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1y, XLAL_FAILURE);
    memcpy(S1ytmp->data->data, S1y->data->data, S1ytmp->data->length*sizeof(REAL8) );
  }

  if ( S1z ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1z, XLAL_FAILURE);
    memcpy(S1ztmp->data->data, S1z->data->data, S1ztmp->data->length*sizeof(REAL8) );
  }

  if ( S2x ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2x, XLAL_FAILURE);
    memcpy(S2xtmp->data->data, S2x->data->data, S2xtmp->data->length*sizeof(REAL8) );
  }

  if ( S2y ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2y, XLAL_FAILURE);
    memcpy(S2ytmp->data->data, S2y->data->data, S2ytmp->data->length*sizeof(REAL8) );
  }

  if ( S2z ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2z, XLAL_FAILURE);
    memcpy(S2ztmp->data->data, S2z->data->data, S2ztmp->data->length*sizeof(REAL8) );
  }

  REAL8 Psi,v,v2,v3;
  REAL8 c2Pp3a,c2Pm3a,c2Pp2a,c2Pm2a,cPp2a,cPm2a,c2Ppa,c2Pma,cPpa,cPma,c2P,cP;
  REAL8 s2Pp3a,s2Pm3a,s2Pp2a,s2Pm2a,sPp2a,sPm2a,s2Ppa,s2Pma,sPpa,sPma,s2P,sP;
  REAL8 ca,sa,c2a,s2a,c3a,s3a;
  REAL8 Sax,Ssx,Say,Ssy;
  REAL8 Saz,Ssz;
  //REAL8 e2x,nx,lx;
  REAL8 e2y,e2z,ny,nz,ly,lz;
  LALSimInspiralInclAngle *an=NULL;
  COMPLEX16TimeSeries *h22 =XLALCreateCOMPLEX16TimeSeries( "h22",  &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h2m2=XLALCreateCOMPLEX16TimeSeries( "h2-2", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h21 =XLALCreateCOMPLEX16TimeSeries( "h21",  &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h2m1=XLALCreateCOMPLEX16TimeSeries( "h2-1", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h20 =XLALCreateCOMPLEX16TimeSeries( "h20",  &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  REAL8 amp22=-8.*eta*(m1+m2)*LAL_G_SI/LAL_C_SI/LAL_C_SI / distance * sqrt(LAL_PI / 5.);
  REAL8 re023p,re13,re2S,re3S,im023p,im13,im2S,im3S;

  REAL8 const sqrt1p5=sqrt(1.5);

  for (idx=0; idx<V->data->length; idx++) {
    v=V->data->data[idx];
    v2=v*v;
    v3=v2*v;
    Psi=Phi->data->data[idx];
    an=XLALSimInspiralComputeInclAngle(LNhz->data->data[idx]);
    //e2x=LNhy->data->data[idx]*e1z->data->data[idx]-LNhz->data->data[idx]*e1y->data->data[idx];
    e2y=LNhz->data->data[idx]*e1x->data->data[idx]-LNhx->data->data[idx]*e1z->data->data[idx];
    e2z=LNhx->data->data[idx]*e1y->data->data[idx]-LNhy->data->data[idx]*e1x->data->data[idx];
    //nx=e1x->data->data[idx]*cos(Psi)+e2x*sin(Psi);
    ny=e1y->data->data[idx]*cos(Psi)+e2y*sin(Psi);
    nz=e1z->data->data[idx]*cos(Psi)+e2z*sin(Psi);
    //lx=e2x*cos(Psi)-e1x->data->data[idx]*sin(Psi);
    ly=e2y*cos(Psi)-e1y->data->data[idx]*sin(Psi);
    lz=e2z*cos(Psi)-e1z->data->data[idx]*sin(Psi);
    if (an->si>EPS) {
      ca=LNhx->data->data[idx]/an->si;
      sa=LNhy->data->data[idx]/an->si;
      sP=nz/an->si;
      cP=lz/an->si;
    }
    else {
      // If L//N then alpha is ill defined, only alpha + phi is defined
      ca=1.;
      sa=0.;
      cP=ca*ny-sa*ly;
      sP=-sa*ny-ca*ly;
    }
    c2a=2.*ca*ca-1.;
    s2a=2.*ca*sa;
    c3a=c2a*ca-s2a*sa;
    s3a=c2a*sa+s2a*ca;
    c2P=cP*cP-sP*sP;
    s2P=2.*cP*sP;
    c2Pp2a=c2P*c2a-s2P*s2a;
    c2Pm2a=c2P*c2a+s2P*s2a;
    s2Pp2a=s2P*c2a+c2P*s2a;
    s2Pm2a=s2P*c2a-c2P*s2a;
    c2Pp3a=c2P*c3a-s2P*s3a;
    c2Pm3a=c2P*c3a+s2P*s3a;
    s2Pp3a=c2P*s3a+s2P*c3a;
    s2Pm3a=-c2P*s3a+s2P*c3a;
    c2Ppa=c2P*ca-s2P*sa;
    c2Pma=c2P*ca+s2P*sa;
    s2Ppa=s2P*ca+c2P*sa;
    s2Pma=s2P*ca-c2P*sa;
    cPpa=cP*ca-sP*sa;
    cPma=cP*ca+sP*sa;
    sPpa=sP*ca+cP*sa;
    sPma=sP*ca-cP*sa;
    cPp2a=cP*c2a-sP*s2a;
    cPm2a=cP*c2a+sP*s2a;
    sPp2a=sP*c2a+cP*s2a;
    sPm2a=sP*c2a-cP*s2a;

    Sax=( S1xtmp->data->data[idx]-S2xtmp->data->data[idx])/2.;
    Ssx=( S1xtmp->data->data[idx]+S2xtmp->data->data[idx])/2.;
    Say=( S1ytmp->data->data[idx]-S2ytmp->data->data[idx])/2.;
    Ssy=( S1ytmp->data->data[idx]+S2ytmp->data->data[idx])/2.;
    Saz=( S1ztmp->data->data[idx]-S2ztmp->data->data[idx])/2.;
    Ssz=( S1ztmp->data->data[idx]+S2ztmp->data->data[idx])/2.;

    //l=2, m=2,-2
    re023p = (1. + v2*amp2 + v3*amp3pi ) * ( c2Pp2a*an->ciBy2Qu + c2Pm2a*an->siBy2Qu );
    re13 = v * dm/3.*amp1*an->si*(1.+v2*amp3)*( cPp2a*an->ciBy2Sq + cPm2a*an->siBy2Sq);
    re2S = v2*amp2S * 0.5 *( (-cPpa*an->ciBy2Sq - cPma*an->siBy2Sq )*(Sax+dm*Ssx) + (-sPma*an->siBy2Sq + sPpa*an->ciBy2Sq)*(Say+dm*Ssy) );
    re3S = ( ca*( 1.-an->siSq*c2a) - 5./3.*( an->ciBy2Qu*c2Pp3a+an->siBy2Qu*c2Pm3a) + 1./6.*( (1.-5.*an->ci)*an->ciBy2Sq*c2Ppa+(1.+5.*an->ci)*an->siBy2Sq*c2Pma) )*(Ssx+dm*Sax);
    re3S+= eta*(-0.5*ca*(1.-an->siSq*c2a) -1./6.*( an->ciBy2Qu*c2Pp3a+an->siBy2Qu*c2Pm3a) + 1./12.*(( 9.-an->ci)*an->ciBy2Sq*c2Ppa+(9.+an->ci)*an->siBy2Sq*c2Pma) )*Ssx;
    re3S+= (-sa*(1.+an->siSq*c2a) + 5./3.*(-an->ciBy2Qu*s2Pp3a+an->siBy2Qu*s2Pm3a) + 1./6.*( (-1.+5.*an->ci)*an->ciBy2Sq*s2Ppa+(1.+5.*an->ci)*an->siBy2Sq*s2Pma) )*(Ssy+dm*Say);
    re3S+= eta*( 0.5*sa*(1.+c2a*an->siSq) + 1./6.*(-an->ciBy2Qu*s2Pp3a+an->siBy2Qu*s2Pm3a) + 1./12.*((-9.+an->ci)*an->ciBy2Sq*s2Ppa+(9.+an->ci)*an->siBy2Sq*s2Pma) )*Ssy;
    re3S*= an->si;
    re3S+= (-c2a*an->ci*an->siSq + 2.*((1.-5./3.*an->ci)*an->ciBy2Qu*c2Pp2a-(1.+5./3.*an->ci)*an->siBy2Qu*c2Pm2a))*(Ssz+dm*Saz);
    re3S+= eta*(0.5*c2a*an->siSq*an->ci + 1./3.*((5.-an->ci)*an->ciBy2Qu*c2Pp2a - (5.+an->ci)*an->siBy2Qu*c2Pm2a) )*Ssz;
    re3S*=amp3S*v3;

    im023p = (1. + v2*amp2 + v3*amp3pi ) * (-s2Pp2a*an->ciBy2Qu + s2Pm2a*an->siBy2Qu );
    im13 = v * dm/3.*amp1*an->si*(1.+v2*amp3)*(-sPp2a*an->ciBy2Sq + sPm2a*an->siBy2Sq);
    im2S = v2*amp2S * 0.5 *( ( sPpa*an->ciBy2Sq - sPma*an->siBy2Sq)*(Sax+dm*Ssx) + ( cPpa*an->ciBy2Sq + cPma*an->siBy2Sq )*(Say+dm*Ssy) );
    im3S = ( -sa*(1.-2.*an->siSq*ca*ca) + 5./3.*( an->ciBy2Qu*s2Pp3a-an->siBy2Qu*s2Pm3a) + 1./6.*(-(1.-5.*an->ci)*an->ciBy2Sq*s2Ppa+(1.+5.*an->ci)*an->siBy2Sq*s2Pma) )*(Ssx+dm*Sax);
    im3S+= eta*(0.5*sa*(1.-2.*an->siSq*ca*ca)+1./6.*( an->ciBy2Qu*s2Pp3a-an->siBy2Qu*s2Pm3a) + 1./12.*((-9.+an->ci)*an->ciBy2Sq*s2Ppa+(9.+an->ci)*an->siBy2Sq*s2Pma) )*Ssx;
    im3S+= ( -ca*(1.-2.*an->siSq*sa*sa) - 5./3.*(an->ciBy2Qu*c2Pp3a+an->siBy2Qu*c2Pm3a) - 1./6.*((1.-5.*an->ci)*an->ciBy2Sq*c2Ppa+(1.+5.*an->ci)*an->siBy2Sq*c2Pma) )*(Ssy+dm*Say);
    im3S+= eta*( 0.5*ca*(1.-2.*an->siSq*sa*sa) - 1./6.*(an->ciBy2Qu*c2Pp3a + an->siBy2Qu*c2Pm3a) + 1./12.*((-9.+an->ci)*an->ciBy2Sq*c2Ppa-(9.+an->ci)*an->siBy2Sq*c2Pma) )*Ssy;
    im3S*= an->si;
    im3S+= ( s2a*an->ci*an->siSq - 2.*((1.-5./3.*an->ci)*an->ciBy2Qu*s2Pp2a + (1.+5./3.*an->ci)*an->siBy2Qu*s2Pm2a) ) * (Ssz+dm*Saz);
    im3S+= eta*(-0.5*s2a*an->ci*an->siSq - 1./3.*((5.-an->ci)*an->ciBy2Qu*s2Pp2a + (5.+an->ci)*an->siBy2Qu*s2Pm2a ) )*Ssz;
    im3S*=amp3S*v3;

    h22->data->data[idx] =amp22*v2*(re023p+re13+re2S+re3S+I*( im023p+im13+im2S+im3S));
    h2m2->data->data[idx]=amp22*v2*(re023p-re13-re2S+re3S+I*(-im023p+im13+im2S-im3S));

    //l=2, m=1,-1
    re023p = an->si * ( 1. + v2*amp2 + v3*amp3pi ) * ( c2Ppa * an->ciBy2Sq - c2Pma * an->siBy2Sq );
    re13 = v * dm/3.*amp1*(1.+v2*amp3)*(-cPpa*(an->ci+an->c2i)/2. -cPma*an->siBy2Sq*(1.+2.*an->ci) );
    re2S = v2*amp2S * 0.5 * ( an->si*sP*(Say+dm*Ssy) + (cPma*an->siBy2Sq+cPpa*an->ciBy2Sq)*(Saz+dm*Ssz) );
    re3S = an->si*(ca*an->c2i + (1.-10./3.*an->ci)*an->ciBy2Sq*c2Ppa + (1.+10./3.*an->ci)*an->siBy2Sq*c2Pma)*(Ssz+dm*Saz);
    re3S+= an->si*eta*( (5.-2.*an->ci)/6.*an->ciBy2Sq*c2Ppa+an->siBy2Sq*(5.+2.*an->ci)/6.*c2Pma - 0.5*ca*an->c2i )*Ssz;
    re3S+= 1./3.*((-7.+10.*an->ci)*an->ciBy2Qu*s2Pp2a-(7.+10.*an->ci)*an->siBy2Qu*s2Pm2a + .5*an->siSq*s2P + 3.*an->ci*an->siSq*s2a) * (dm*Say+Ssy);
    re3S+= eta*((0.5+an->ci/3.)*an->ciBy2Qu*s2Pp2a+(0.5-an->ci/3.)*an->siBy2Qu*s2Pm2a -1.3/1.2*an->siSq*s2P - 0.5*an->ci*an->siSq*s2a) * Ssy;
    re3S+= 1./3.*((-7.+10.*an->ci)*an->ciBy2Qu*c2Pp2a+(7.+10.*an->ci)*an->siBy2Qu*c2Pm2a - 5.*an->siSq*an->ci*c2P + 3.*an->ci*(an->siSq*c2a-an->ciSq)) * (dm*Sax+Ssx);
    re3S+= eta*((0.5+an->ci/3.)*an->ciBy2Qu*c2Pp2a-(0.5-an->ci/3.)*an->siBy2Qu*c2Pm2a -1./6.*an->ci*an->siSq*c2P - 0.5*an->ci*(an->siSq*c2a-an->ciSq)) * Ssx;
    re3S*=amp3S*v3;

    im023p = an->si * ( 1. + v2*amp2 + v3*amp3pi ) * (-s2Ppa * an->ciBy2Sq - s2Pma * an->siBy2Sq );
    im13 = v * dm/3.*amp1*(1.+v2*amp3)*( sPpa*(an->ci+an->c2i)/2. -sPma*an->siBy2Sq*(1.+2.*an->ci) );
    im2S = v2*amp2S * 0.5 * ( an->si*sP*(Sax+dm*Ssx) + (sPma*an->siBy2Sq-sPpa*an->ciBy2Sq)*(Saz+dm*Ssz) );
    im3S = an->si*(-sa*an->c2i + (-1.+10./3.*an->ci)*an->ciBy2Sq*s2Ppa + (1.+10./3.*an->ci)*an->siBy2Sq*s2Pma)*(Ssz+dm*Saz);
    im3S+= an->si*eta*( (-5.+2.*an->ci)/6.*an->ciBy2Sq*s2Ppa+an->siBy2Sq*(5.+2.*an->ci)/6.*s2Pma+0.5*sa*an->c2i)*Ssz;
    im3S+= 1./3.*((-7.+10.*an->ci)*an->ciBy2Qu*c2Pp2a+(7.+10.*an->ci)*an->siBy2Qu*c2Pm2a + 5.*an->siSq*an->ci*c2P + 3.*an->ci*(an->siSq*c2a+an->ciSq)) * (dm*Say+Ssy);
    im3S+= eta*((0.5+an->ci/3.)*an->ciBy2Qu*c2Pp2a+(-0.5+an->ci/3.)*an->siBy2Qu*c2Pm2a + 1./6.*an->ci*an->siSq*c2P - 0.5*an->ci*(an->siSq*c2a+an->ciSq)) * Ssy;
    im3S+= 1./3.*((7.-10.*an->ci)*an->ciBy2Qu*s2Pp2a+(7.+10.*an->ci)*an->siBy2Qu*s2Pm2a + .5*an->siSq*s2P - 3.*an->ci*an->siSq*s2a) * (dm*Sax+Ssx);
    im3S+= eta*(-(0.5+an->ci/3.)*an->ciBy2Qu*s2Pp2a-(0.5-an->ci/3.)*an->siBy2Qu*s2Pm2a -1.3/1.2*an->siSq*s2P + 0.5*an->ci*an->siSq*s2a) * Ssx;
    im3S*=amp3S*v3;

    h21->data->data[idx] =amp22*v2*( re023p+re13+re2S+re3S+I*( im023p+im13+im2S+im3S));
    h2m1->data->data[idx]=amp22*v2*(-re023p+re13+re2S-re3S+I*( im023p-im13-im2S+im3S));

    //l=2, m=0
    //Real components
    re023p = ( 1. + v2*amp2 + v3*amp3pi) * ( an->siSq*c2P );
    re13 = 0.;
    re2S = 0.;
    re3S = an->s2i/3.*( (3.-5.*c2P)*(Ssz+dm*Saz) - eta/2.*(3.+c2P)*Ssz );
    re3S+=(-2.*ca*an->ciSq + 2./3.*(5.*an->ci-2.)*an->ciBy2Sq*c2Ppa - 2./3.*(5.*an->ci+2.)*an->siBy2Sq*c2Pma)*(Ssx+dm*Sax);
    re3S+=eta*( ca*an->ciSq + an->ciBy2Sq/3.*(an->ci+4.)*c2Ppa + an->siBy2Sq/3.*(4.-an->ci)*c2Pma)*Ssx;
    re3S+=(-2.*sa*an->ciSq + 2./3.*(5.*an->ci-2.)*an->ciBy2Sq*s2Ppa + 2./3.*(5.*an->ci+2.)*an->siBy2Sq*s2Pma)*(Ssy+dm*Say);
    re3S+=eta*(sa*an->ciSq + 1./3.*(4.+an->ci)*an->ciBy2Sq*s2Ppa - an->siBy2Sq*(4.-an->ci)/3.*s2Pma)*Ssy;
    re3S*=amp3S*v3*an->si;
    //Imaginary components
    im023p = 0.;
    im13 = v*amp1*dm/3.*(1.+v2*amp3)*an->s2i*sP;
    im2S = v2*amp2S/3.* (-2.*sP*an->si*(Saz+dm*Ssz) + (-sPpa*an->ciBy2Sq+sPma*an->siBy2Sq)*(Sax+dm*Ssx) + (cPpa*an->ciBy2Sq+cPma*an->siBy2Sq)*(Say+dm*Ssy) );
    im3S = 0.;
    //h20
    h20->data->data[idx]=amp22*v2*sqrt1p5*(re023p+re13+re2S+re3S+I*(im023p+im13+im2S+im3S));

    LALFree(an); an=NULL;

  }

  XLALDestroyREAL8TimeSeries(S1xtmp);
  XLALDestroyREAL8TimeSeries(S1ytmp);
  XLALDestroyREAL8TimeSeries(S1ztmp);
  XLALDestroyREAL8TimeSeries(S2xtmp);
  XLALDestroyREAL8TimeSeries(S2ytmp);
  XLALDestroyREAL8TimeSeries(S2ztmp);

  *h2ms=XLALSphHarmTimeSeriesAddMode(*h2ms, h22, 2, 2);
  XLALDestroyCOMPLEX16TimeSeries(h22);
  *h2ms=XLALSphHarmTimeSeriesAddMode(*h2ms, h21, 2, 1);
  XLALDestroyCOMPLEX16TimeSeries(h21);
  *h2ms=XLALSphHarmTimeSeriesAddMode(*h2ms, h20, 2, 0);
  XLALDestroyCOMPLEX16TimeSeries(h20);
  *h2ms=XLALSphHarmTimeSeriesAddMode(*h2ms, h2m1, 2, -1);
  XLALDestroyCOMPLEX16TimeSeries(h2m1);
  *h2ms=XLALSphHarmTimeSeriesAddMode(*h2ms, h2m2, 2, -2);
  XLALDestroyCOMPLEX16TimeSeries(h2m2);

  return XLAL_SUCCESS;

} /* End of XLALSimSpinInspiralPNMode2m*/

/**
 * Computes all l=3 modes of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform for spinning-precessing waveforms.
 *
 * For reference see eq. B3 of
 * Arun et al., "Higher-order spin effects in the amplitude and phase of gravitational waveforms
 * emitted by inspiraling compact binaries: Ready-to-use gravitational waveforms",
 * Physical Review D 79, 104023 (2009), Erratum PRD 84, 049901 (2011), arXiv:0810.5336v3 [gr-qc].
 * All variable are inputs.
 *
 * !!BEWARE!! Spin components are given wrt LNh. LNh components define angles
 * incl and alpha (polarization angle) entering the mode analytic form.
 *
 */

INT4 XLALSimInspiralSpinPNMode3m(SphHarmTimeSeries **h3ms, /**< OUTPUT*/
				 REAL8TimeSeries *V,    /**< post-Newtonian parameter */
				 REAL8TimeSeries *Phi,  /**< orbital phase */
				 REAL8TimeSeries *LNhx, /**< angular momentum unit vector x components */
				 REAL8TimeSeries *LNhy, /**< angular momentum unit vector y components */
				 REAL8TimeSeries *LNhz, /**< angular momentum unit vector z components */
				 REAL8TimeSeries *e1x,  /**< x-axis tetrad x component*/
				 REAL8TimeSeries *e1y,  /**< x-axis tetrad y component*/
				 REAL8TimeSeries *e1z,  /**< x-axis tetrad z component*/
				 REAL8TimeSeries *S1x,  /**< spin1-x component */
				 REAL8TimeSeries *S1y,  /**< spin1-y component */
				 REAL8TimeSeries *S1z,  /**< spin1-z component */
				 REAL8TimeSeries *S2x,  /**< spin2-x component */
				 REAL8TimeSeries *S2y,  /**< spin2-y component */
				 REAL8TimeSeries *S2z,  /**< spin2-z component */
				 REAL8 m1,              /**< mass of companion 1 (kg) */
				 REAL8 m2,              /**< mass of companion 2 (kg) */
				 REAL8 distance,        /**< distance of source (m) */
				 int ampO               /**< twice post-Newtonian amplitude order */)
{
  LAL_CHECK_VALID_SERIES(V, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(Phi, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(LNhx, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(LNhy, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(LNhz, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(e1x, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(e1y, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(e1z, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhx, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhy, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhz, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, e1x, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, e1y, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, e1z, XLAL_FAILURE);

  UINT4 idx;
  REAL8 eta=m1/(m1+m2)*m2/(m1+m2);
  REAL8 dm=(m1-m2)/(m1+m2);

  REAL8 amp3 =0.;
  REAL8 amp3S=0.;
  REAL8 amp2 =0.;
  REAL8 amp1 =0.;

  switch (ampO) {
  case -1: /* use highest available pN order */
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
    __attribute__ ((fallthrough));
#endif
  case 3:
    amp3 = 1.;
    amp3S= 1.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
    __attribute__ ((fallthrough));
#endif
  case 2:
    amp2 = (1.-3.*eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
    __attribute__ ((fallthrough));
#endif
  case 1:
    amp1 = 1.;
  case 0:
    break;
  default: /* unsupported pN order */
    XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, ampO/2, ampO%2?".5":"" );
    XLAL_ERROR(XLAL_EINVAL);
  }

  REAL8 checkL,checke;
  for (idx=0;idx<LNhx->data->length;idx++){
    checkL=LNhx->data->data[idx]*LNhx->data->data[idx]+LNhy->data->data[idx]*LNhy->data->data[idx]+LNhz->data->data[idx]*LNhz->data->data[idx];
    if ((checkL<1.-EPS)||(checkL>1.+EPS)) {
      XLALPrintError("XLAL Error - %s: LNhat unit vector has not unit norm %12.4e\n", __func__,checkL);
      XLAL_ERROR(XLAL_EINVAL);
    }
    checke=e1x->data->data[idx]*e1x->data->data[idx]+e1y->data->data[idx]*e1y->data->data[idx]+e1z->data->data[idx]*e1z->data->data[idx];
    if ((checke<1.-EPS)||(checke>1.+EPS)) {
      XLALPrintError("XLAL Error - %s: e1 unit vector has not unit norm %12.4e at idx %d\n", __func__,checke,idx);
      XLAL_ERROR(XLAL_EINVAL);
    }
  }

  REAL8TimeSeries *S1xtmp=XLALCreateREAL8TimeSeries("S1x",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S1ytmp=XLALCreateREAL8TimeSeries("S1y",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S1ztmp=XLALCreateREAL8TimeSeries("S1z",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S2xtmp=XLALCreateREAL8TimeSeries("S2x",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S2ytmp=XLALCreateREAL8TimeSeries("S2y",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S2ztmp=XLALCreateREAL8TimeSeries("S2z",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  for (idx=0;idx<V->data->length;idx++) {
    S1xtmp->data->data[idx]=0.;
    S1ytmp->data->data[idx]=0.;
    S1ztmp->data->data[idx]=0.;
    S2xtmp->data->data[idx]=0.;
    S2ytmp->data->data[idx]=0.;
    S2ztmp->data->data[idx]=0.;
  }

  if ( S1x ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1x, XLAL_FAILURE);
    memcpy(S1xtmp->data->data, S1x->data->data, S1xtmp->data->length*sizeof(REAL8) );
  }

  if ( S1y ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1y, XLAL_FAILURE);
    memcpy(S1ytmp->data->data, S1y->data->data, S1ytmp->data->length*sizeof(REAL8) );
  }

  if ( S1z ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1z, XLAL_FAILURE);
    memcpy(S1ztmp->data->data, S1z->data->data, S1ztmp->data->length*sizeof(REAL8) );
  }

  if ( S2x ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2x, XLAL_FAILURE);
    memcpy(S2xtmp->data->data, S2x->data->data, S2xtmp->data->length*sizeof(REAL8) );
  }

  if ( S2y ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2y, XLAL_FAILURE);
    memcpy(S2ytmp->data->data, S2y->data->data, S2ytmp->data->length*sizeof(REAL8) );
  }

  if ( S2z ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2z, XLAL_FAILURE);
    memcpy(S2ztmp->data->data, S2z->data->data, S2ztmp->data->length*sizeof(REAL8) );
  }

  REAL8 Psi,v,v2,v3;
  REAL8 c3Pp3a,c3Pm3a,c3Pp2a,c3Pm2a,c3Ppa,c3Pma,c2Pp3a,c2Pm3a,c2Pp2a,c2Pm2a,c2Ppa,c2Pma,cPp3a,cPm3a,cPp2a,cPm2a,cPpa,cPma,cP,c2P,c3P;
  REAL8 s3Pp3a,s3Pm3a,s3Pp2a,s3Pm2a,s3Ppa,s3Pma,s2Pp3a,s2Pm3a,s2Pp2a,s2Pm2a,s2Ppa,sPp2a,sPm2a,s2Pma,sPpa,sPma,sPp3a,sPm3a,sP,s2P,s3P;
  REAL8 ca,sa,c2a,s2a,c3a,s3a;
  //REAL8 Sax, Say, Saz;
  REAL8 Ssx,Ssy,Ssz;
  //REAL8 e2x,nx,lx;
  REAL8 e2y,e2z,ny,nz,ly,lz;
  LALSimInspiralInclAngle *an=NULL;
  COMPLEX16TimeSeries *h33=XLALCreateCOMPLEX16TimeSeries( "h33", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h3m3=XLALCreateCOMPLEX16TimeSeries( "h3-3", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h32=XLALCreateCOMPLEX16TimeSeries( "h32", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h3m2=XLALCreateCOMPLEX16TimeSeries( "h3-2", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h31=XLALCreateCOMPLEX16TimeSeries( "h31", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h3m1=XLALCreateCOMPLEX16TimeSeries( "h3-1", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h30=XLALCreateCOMPLEX16TimeSeries( "h30", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);

  REAL8 amp33= eta*(m1+m2)*LAL_G_SI/pow(LAL_C_SI,2) / distance * sqrt(2.*LAL_PI / 21.);
  REAL8 re3,re4,re5,re5S,im3,im4,im5,im5S;

  for (idx=0; idx<V->data->length; idx++) {
    v=V->data->data[idx];
    v2=v*v;
    v3=v2*v;
    Psi=Phi->data->data[idx];
    an=XLALSimInspiralComputeInclAngle(LNhz->data->data[idx]);
    //e2x=LNhy->data->data[idx]*e1z->data->data[idx]-LNhz->data->data[idx]*e1y->data->data[idx];
    e2y=LNhz->data->data[idx]*e1x->data->data[idx]-LNhx->data->data[idx]*e1z->data->data[idx];
    e2z=LNhx->data->data[idx]*e1y->data->data[idx]-LNhy->data->data[idx]*e1x->data->data[idx];
    //nx=e1x->data->data[idx]*cos(Psi)+e2x*sin(Psi);
    ny=e1y->data->data[idx]*cos(Psi)+e2y*sin(Psi);
    nz=e1z->data->data[idx]*cos(Psi)+e2z*sin(Psi);
    //lx=e2x*cos(Psi)-e1x->data->data[idx]*sin(Psi);
    ly=e2y*cos(Psi)-e1y->data->data[idx]*sin(Psi);
    lz=e2z*cos(Psi)-e1z->data->data[idx]*sin(Psi);
    if (an->si>EPS) {
      ca=LNhx->data->data[idx]/an->si;
      sa=LNhy->data->data[idx]/an->si;
      sP=nz/an->si;
      cP=lz/an->si;
    }
    else {
      // If L//N then alpha is ill defined, only alpha + phi is defined
      ca=1.;
      sa=0.;
      cP=ca*ny-sa*ly;
      sP=-sa*ny-ca*ly;
    }
    c2a=2.*ca*ca-1.;
    s2a=2.*ca*sa;
    c3a=c2a*ca-s2a*sa;
    s3a=c2a*sa+s2a*ca;
    c2P=cP*cP-sP*sP;
    s2P=2.*cP*sP;
    c3P=c2P*cP-s2P*sP;
    s3P=c2P*sP+s2P*cP;
    c3Pp3a=c3P*c3a-s3P*s3a;
    c3Pm3a=c3P*c3a+s3P*s3a;
    s3Pp3a=s3P*c3a+c3P*s3a;
    s3Pm3a=s3P*c3a-c3P*s3a;
    c3Pp2a=c3P*c2a-s3P*s2a;
    c3Pm2a=c3P*c2a+s3P*s2a;
    s3Pp2a=s3P*c2a+c3P*s2a;
    s3Pm2a=s3P*c2a-c3P*s2a;
    c3Ppa=c3P*ca-s3P*sa;
    c3Pma=c3P*ca+s3P*sa;
    s3Ppa=s3P*ca+c3P*sa;
    s3Pma=s3P*ca-c3P*sa;
    c2Pp3a=c2P*c3a-s2P*s3a;
    c2Pm3a=c2P*c3a+s2P*s3a;
    s2Pp3a=s2P*c3a+c2P*s3a;
    s2Pm3a=s2P*c3a-c2P*s3a;
    c2Pp2a=c2P*c2a-s2P*s2a;
    c2Pm2a=c2P*c2a+s2P*s2a;
    s2Pp2a=s2P*c2a+c2P*s2a;
    s2Pm2a=s2P*c2a-c2P*s2a;
    c2Ppa=c2P*ca-s2P*sa;
    c2Pma=c2P*ca+s2P*sa;
    s2Ppa=s2P*ca+c2P*sa;
    s2Pma=s2P*ca-c2P*sa;
    cPp3a=cP*c3a-sP*s3a;
    cPm3a=cP*c3a+sP*s3a;
    sPp3a=sP*c3a+cP*s3a;
    sPm3a=sP*c3a-cP*s3a;
    cPp2a=cP*c2a-sP*s2a;
    cPm2a=cP*c2a+sP*s2a;
    sPp2a=sP*c2a+cP*s2a;
    sPm2a=sP*c2a-cP*s2a;
    cPpa=cP*ca-sP*sa;
    cPma=cP*ca+sP*sa;
    sPpa=sP*ca+cP*sa;
    sPma=sP*ca-cP*sa;

    //Sax=S1xtmp->data->data[idx]-S2xtmp->data->data[idx];
    Ssx=S1xtmp->data->data[idx]+S2xtmp->data->data[idx];
    //Say=S1ytmp->data->data[idx]-S2ytmp->data->data[idx];
    Ssy=S1ytmp->data->data[idx]+S2ytmp->data->data[idx];
    //Saz=S1ztmp->data->data[idx]-S2ztmp->data->data[idx];
    Ssz=S1ztmp->data->data[idx]+S2ztmp->data->data[idx];

    //l=3, m=3,-3
    re3 = amp1*v*dm* (-9.*c3Pm3a*an->siBy2Sx - cPm3a*an->siBy2Qu*an->ciBy2Sq + cPp3a*an->siBy2Sq*an->ciBy2Qu + 9.*c3Pp3a*an->ciBy2Sx);
    re4 = -4.*amp2*v2*an->si*(c2Pm3a*an->siBy2Qu - c2Pp3a*an->ciBy2Qu);
    re5 = amp3*v3*dm* (36.*(1.-0.5*eta)*(-an->ciBy2Sx*c3Pp3a + an->siBy2Sx*c3Pm3a) + 2./3.*(1.+0.25*eta)*an->siSq*(-an->ciBy2Sq*cPp3a+an->siBy2Sq*cPm3a) );
    re5S= amp3S*v3*eta*16.*( (-an->siBy2Qu*c2Pm2a + an->ciBy2Qu*c2Pp2a)*Ssx + (-an->siBy2Qu*s2Pm2a - an->ciBy2Qu*s2Pp2a)*Ssy );
    im3 = amp1*v*dm* (-9.*s3Pm3a*an->siBy2Sx - sPm3a*an->siBy2Qu*an->ciBy2Sq - sPp3a*an->siBy2Sq*an->ciBy2Qu - 9.*s3Pp3a*an->ciBy2Sx);
    im4 = -4.*amp2*v2*an->si*(s2Pm3a*an->siBy2Qu + s2Pp3a*an->ciBy2Qu);
    im5 = amp3*v3*dm* (36.*(1.-0.5*eta)*( an->ciBy2Sx*s3Pp3a + an->siBy2Sx*s3Pm3a) + 2./3.*(1.+0.25*eta)*an->siSq*( an->ciBy2Sq*sPp3a+an->siBy2Sq*sPm3a) );
    im5S= amp3S*v3*eta*16.*( (-an->siBy2Qu*s2Pm2a - an->ciBy2Qu*s2Pp2a)*Ssx + ( an->siBy2Qu*c2Pm2a - an->ciBy2Qu*c2Pp2a)*Ssy );
    h33->data->data[idx] =amp33*v2*( re3+re4+re5+re5S+I*(im3+im4+im5+im5S));
    h3m3->data->data[idx]=amp33*v2*(-re3+re4-re5+re5S+I*(im3-im4+im5-im5S));

    //l=3, m=2,-2
    re3 = amp1*v*dm*an->si/4.*( 18.*c3Pm2a*an->siBy2Qu + cPm2a/3.*(1.+3.*an->ci)*an->siBy2Sq + cPp2a/3.*(1.-3.*an->ci)*an->ciBy2Sq + 18.*c3Pp2a*an->ciBy2Qu );
    re4 = amp2*v2*4.*( an->ciBy2Qu*c2Pp2a*(2./3.-an->ci) + an->siBy2Qu*c2Pm2a*(2./3.+an->ci) );
    re5 = amp3*v3*dm*an->si*(-18.*(1.-0.5*eta)*(an->ciBy2Qu*c3Pp2a + an->siBy2Qu*c3Pm2a) + (1.+0.25*eta)*2./9.*((-1.+3.*an->ci)*an->ciBy2Sq*cPp2a - (1.+3.*an->ci)*an->siBy2Sq*cPm2a ) );
    re5S= amp3S*v3*eta*16./3.*( an->si*(( an->ciBy2Sq*c2Ppa + an->siBy2Sq*c2Pma)*Ssx + (-an->ciBy2Sq*s2Ppa + an->siBy2Sq*s2Pma)*Ssy) + (-an->ciBy2Qu*c2Pp2a + an->siBy2Qu*c2Pm2a)*Ssz );
    im3 = amp1*v*dm*an->si/4.*( 18.*s3Pm2a*an->siBy2Qu + sPm2a/3.*(1.+3.*an->ci)*an->siBy2Sq - sPp2a/3.*(1.-3.*an->ci)*an->ciBy2Sq - 18.*s3Pp2a*an->ciBy2Qu );
    im4 = amp2*v2*4.*(-an->ciBy2Qu*s2Pp2a*(2./3.-an->ci) + an->siBy2Qu*s2Pm2a*(2./3.+an->ci) );
    im5 = amp3*v3*dm*an->si*( 18.*(1.-0.5*eta)*(an->ciBy2Qu*s3Pp2a - an->siBy2Qu*s3Pm2a) + (1.+0.25*eta)*2./9.*(( 1.-3.*an->ci)*an->ciBy2Sq*sPp2a - (1.+3.*an->ci)*an->siBy2Sq*sPm2a ) );
    im5S= amp3S*v3*eta*16./3.*( an->si*((-an->ciBy2Sq*s2Ppa + an->siBy2Sq*s2Pma)*Ssx + (-an->ciBy2Sq*c2Ppa - an->siBy2Sq*c2Pma)*Ssy) + ( an->ciBy2Qu*s2Pp2a + an->siBy2Qu*s2Pm2a)*Ssz );
    h32->data->data[idx] =amp33*sqrt(6.)*v2*(re3+re4+re5+re5S+I*( im3+im4+im5+im5S));
    h3m2->data->data[idx]=amp33*sqrt(6.)*v2*(re3-re4+re5-re5S+I*(-im3+im4-im5+im5S));

    //l=3, m=1,-1
    re3 = amp1*v*dm/4.*( 45.*c3Ppa*an->siSq*an->ciBy2Sq - 180.*c3Pma*an->siBy2Qu*an->ciBy2Sq - cPma/2.*an->siBy2Sq*(13./3.+20./3.*an->ci+5.*an->c2i) + cPpa/2.*an->ciBy2Sq*(13./3.-20./3.*an->ci+5.*an->c2i) );
    re4 = amp2*v2*10./3.*an->si*(an->ciBy2Sq*c2Ppa*(1.-3.*an->ci) - an->siBy2Sq*c2Pma*(1.+3.*an->ci));
    re5 = amp3*v3*dm*( -180.*(1.-0.5*eta)*an->ciBy2Sq*an->siBy2Sq*( an->ciBy2Sq*c3Ppa - an->siBy2Sq*c3Pma) + 8./9.*(1.+0.25*eta)*(-an->ciBy2Sq*(13./8.-2.5*an->ci+15./8.*an->c2i)*cPpa + an->siBy2Sq*(13./8.+2.5*an->ci+15./8.*an->c2i)*cPma ) );
    re5S= amp3S*v3*eta*16./3.*( 4.*an->si*(-an->ciBy2Sq*c2Ppa - an->siBy2Sq*c2Pma)*Ssz + (-an->ciBy2Qu*c2Pp2a + an->siBy2Qu*c2Pm2a)*Ssx + (-an->ciBy2Qu*s2Pp2a - an->siBy2Qu*s2Pm2a - 3.*an->siSq*s2P)*Ssy );
    im3 = amp1*v*dm/4.*(-45.*s3Ppa*an->siSq*an->ciBy2Sq - 180.*s3Pma*an->siBy2Qu*an->ciBy2Sq - sPma/2.*an->siBy2Sq*(13./3.+20./3.*an->ci+5.*an->c2i) - sPpa/2.*an->ciBy2Sq*(13./3.-20./3.*an->ci+5.*an->c2i) );
    im4 = amp2*v2*10./3.*an->si*(-an->ciBy2Sq*s2Ppa*(1.-3.*an->ci) - an->siBy2Sq*s2Pma*(1.+3.*an->ci));
    im5 = amp3*v3*dm*( -180.*(1.-0.5*eta)*an->ciBy2Sq*an->siBy2Sq*(-an->ciBy2Sq*s3Ppa - an->siBy2Sq*s3Pma) + 8./9.*(1.+0.25*eta)*( an->ciBy2Sq*(13./8.-2.5*an->ci+15./8.*an->c2i)*sPpa + an->siBy2Sq*(13./8.+2.5*an->ci+15./8.*an->c2i)*sPma ) );
    im5S= amp3S*v3*eta*16./3.*( 4*an->si*( an->ciBy2Sq*s2Ppa - an->siBy2Sq*s2Pma)*Ssz + ( an->ciBy2Qu*s2Pp2a + an->siBy2Qu*s2Pm2a -3.*an->siSq*s2P)*Ssx + (-an->ciBy2Qu*c2Pp2a + an->siBy2Qu*c2Pm2a)*Ssy );
    h31->data->data[idx] =amp33*sqrt(.6)*v2*( re3+re4+re5+re5S+I*(im3+im4+im5+im5S));
    h3m1->data->data[idx]=amp33*sqrt(.6)*v2*(-re3+re4-re5+re5S+I*(im3-im4+im5-im5S));

    //l=3, m=0
    re3 = amp1*v*dm*an->si*( cP/4.*(3.+5.*an->c2i) +22.5*c3P*an->siSq);
    re4 = 0.;
    re5 = amp3*v3*dm*an->si*(-90.*(1.-0.5*eta)*an->siSq*c3P - 2.*(1.+0.25*eta)*(1.+5./3.*an->c2i)*cP);
    re5S= 0.;
    im3 = 0.;
    im4 = amp2*v2*20.*an->ci*an->siSq*s2P;
    im5 = 0.;
    im5S= amp3S*v3*eta*an->si*16.*( 2.*(an->ciBy2Sq*s2Ppa - an->siBy2Sq*s2Pma)*Ssx - 2.*( an->ciBy2Sq*c2Ppa + an->siBy2Sq*c2Pma)*Ssy + 3.*(an->si*s2P)*Ssz );
    h30->data->data[idx]=amp33/sqrt(5.)*v2*(re3+re4+re5+re5S+I*(im3+im4+im5+im5S));

    LALFree(an);
    an=NULL;

  }

  XLALDestroyREAL8TimeSeries(S1xtmp);
  XLALDestroyREAL8TimeSeries(S1ytmp);
  XLALDestroyREAL8TimeSeries(S1ztmp);
  XLALDestroyREAL8TimeSeries(S2xtmp);
  XLALDestroyREAL8TimeSeries(S2ytmp);
  XLALDestroyREAL8TimeSeries(S2ztmp);

  *h3ms=XLALSphHarmTimeSeriesAddMode(*h3ms, h33, 3, 3);
  XLALDestroyCOMPLEX16TimeSeries(h33);
  *h3ms=XLALSphHarmTimeSeriesAddMode(*h3ms, h32, 3, 2);
  XLALDestroyCOMPLEX16TimeSeries(h32);
  *h3ms=XLALSphHarmTimeSeriesAddMode(*h3ms, h31, 3, 1);
  XLALDestroyCOMPLEX16TimeSeries(h31);
  *h3ms=XLALSphHarmTimeSeriesAddMode(*h3ms, h30, 3, 0);
  XLALDestroyCOMPLEX16TimeSeries(h30);
  *h3ms=XLALSphHarmTimeSeriesAddMode(*h3ms, h3m1, 3, -1);
  XLALDestroyCOMPLEX16TimeSeries(h3m1);
  *h3ms=XLALSphHarmTimeSeriesAddMode(*h3ms, h3m2, 3, -2);
  XLALDestroyCOMPLEX16TimeSeries(h3m2);
  *h3ms=XLALSphHarmTimeSeriesAddMode(*h3ms, h3m3, 3, -3);
  XLALDestroyCOMPLEX16TimeSeries(h3m3);

  return XLAL_SUCCESS;

}  /* End of XLALSimSpinInspiralPNMode3m*/

/**
 * Computes all l=3 modes of spherical harmonic decomposition of
 * the post-Newtonian inspiral waveform for spinning-precessing waveforms.
 *
 * For reference see eq. B3 of
 * Arun et al., "Higher-order spin effects in the amplitude and phase of gravitational waveforms
 * emitted by inspiraling compact binaries: Ready-to-use gravitational waveforms",
 * Physical Review D 79, 104023 (2009), Erratum PRD 84, 049901 (2011), arXiv:0810.5336v3 [gr-qc].
 * All variable are inputs.
 *
 * !!BEWARE!! Spin components are given wrt LNh. LNh components define angles
 * incl and alpha (polarization angle) entering the mode analytic form.
 *
 */
INT4 XLALSimInspiralSpinPNMode4m(SphHarmTimeSeries **h4ms, /**< OUTPUT */
				 REAL8TimeSeries *V,    /**< post-Newtonian parameter */
				 REAL8TimeSeries *Phi,  /**< orbital phase */
				 REAL8TimeSeries *LNhx, /**< angular momentum unit vector x components */
				 REAL8TimeSeries *LNhy, /**< angular momentum unit vector y components */
				 REAL8TimeSeries *LNhz, /**< angular momentum unit vector z components */
				 REAL8TimeSeries *e1x,  /**< x-axis tetrad x component*/
				 REAL8TimeSeries *e1y,  /**< x-axis tetrad y component*/
				 REAL8TimeSeries *e1z,  /**< x-axis tetrad z component*/
				 UNUSED REAL8TimeSeries *S1x,  /**< spin1-x component */
				 UNUSED REAL8TimeSeries *S1y,  /**< spin1-y component */
				 UNUSED REAL8TimeSeries *S1z,  /**< spin1-z component */
				 UNUSED REAL8TimeSeries *S2x,  /**< spin2-x component */
				 UNUSED REAL8TimeSeries *S2y,  /**< spin2-y component */
				 UNUSED REAL8TimeSeries *S2z,  /**< spin2-z component */
				 REAL8 m1,              /**< mass of companion 1 (kg) */
				 REAL8 m2,              /**< mass of companion 2 (kg) */
				 REAL8 distance,        /**< distance of source (m) */
				 int ampO               /**< twice post-Newtonian amplitude order */)
{
  LAL_CHECK_VALID_SERIES(V, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(Phi, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(LNhx, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(LNhy, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(LNhz, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(e1x, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(e1y, XLAL_FAILURE);
  LAL_CHECK_VALID_SERIES(e1z, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhx, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhy, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhz, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, e1x, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, e1y, XLAL_FAILURE);
  LAL_CHECK_CONSISTENT_TIME_SERIES(V, e1z, XLAL_FAILURE);

  UINT4 idx;
  REAL8 eta=m1/(m1+m2)*m2/(m1+m2);
  REAL8 dm=(m1-m2)/(m1+m2);

  REAL8 amp3 =0.;
  REAL8 amp2 =0.;

  switch (ampO) {
  case -1: /* use highest available pN order */
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
    __attribute__ ((fallthrough));
#endif
  case 3:
    amp3 = 1.-2.*eta;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
    __attribute__ ((fallthrough));
#endif
  case 2:
    amp2 = 1.-3.*eta;
  case 1:
  case 0:
    break;
  default: /* unsupported pN order */
    XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__, ampO/2, ampO%2?".5":"" );
    XLAL_ERROR(XLAL_EINVAL);
  }

  REAL8 checkL,checke;
  for (idx=0;idx<LNhx->data->length;idx++){
    checkL=LNhx->data->data[idx]*LNhx->data->data[idx]+LNhy->data->data[idx]*LNhy->data->data[idx]+LNhz->data->data[idx]*LNhz->data->data[idx];
    if ((checkL<1.-EPS)||(checkL>1.+EPS)) {
      XLALPrintError("XLAL Error - %s: LNhat unit vector has not unit norm %12.4e\n", __func__,checkL);
      XLAL_ERROR(XLAL_EINVAL);
    }
    checke=e1x->data->data[idx]*e1x->data->data[idx]+e1y->data->data[idx]*e1y->data->data[idx]+e1z->data->data[idx]*e1z->data->data[idx];
    if ((checke<1.-EPS)||(checke>1.+EPS)) {
      XLALPrintError("XLAL Error - %s: e1 unit vector has not unit norm %12.4e at idx %d\n", __func__,checke,idx);
      XLAL_ERROR(XLAL_EINVAL);
    }
  }

  /*REAL8TimeSeries *S1xtmp=XLALCreateREAL8TimeSeries("S1x",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S1ytmp=XLALCreateREAL8TimeSeries("S1y",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S1ztmp=XLALCreateREAL8TimeSeries("S1z",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S2xtmp=XLALCreateREAL8TimeSeries("S2x",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S2ytmp=XLALCreateREAL8TimeSeries("S2y",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  REAL8TimeSeries *S2ztmp=XLALCreateREAL8TimeSeries("S2z",&V->epoch,0.,V->deltaT,&lalDimensionlessUnit,V->data->length);
  for (idx=0;idx<V->data->length;idx++) {
    S1xtmp->data->data[idx]=0.;
    S1ytmp->data->data[idx]=0.;
    S1ztmp->data->data[idx]=0.;
    S2xtmp->data->data[idx]=0.;
    S2ytmp->data->data[idx]=0.;
    S2ztmp->data->data[idx]=0.;
  }

  if ( S1x ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1x, XLAL_FAILURE);
    memcpy(S1xtmp->data->data, S1x->data->data, S1xtmp->data->length*sizeof(REAL8) );
  }

  if ( S1y ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1y, XLAL_FAILURE);
    memcpy(S1ytmp->data->data, S1y->data->data, S1ytmp->data->length*sizeof(REAL8) );
  }

  if ( S1z ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1z, XLAL_FAILURE);
    memcpy(S1ztmp->data->data, S1z->data->data, S1ztmp->data->length*sizeof(REAL8) );
  }

  if ( S2x ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2x, XLAL_FAILURE);
    memcpy(S2xtmp->data->data, S2x->data->data, S2xtmp->data->length*sizeof(REAL8) );
  }

  if ( S2y ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2y, XLAL_FAILURE);
    memcpy(S2ytmp->data->data, S2y->data->data, S2ytmp->data->length*sizeof(REAL8) );
  }

  if ( S2z ) {
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2z, XLAL_FAILURE);
    memcpy(S2ztmp->data->data, S2z->data->data, S2ztmp->data->length*sizeof(REAL8) );
  }*/

  REAL8 Psi,v,v2,v3;
  REAL8 cP,sP,c2P,s2P,c3P,s3P,c4P,s4P;
  REAL8 ca,sa,c2a,s2a,c3a,s3a,c4a,s4a;
  //REAL8 e2x,nx,lx;
  REAL8 ny,nz,ly,lz,e2y,e2z;
  REAL8 c4Pm4a,c3Pm4a,c2Pm4a,cPm4a,cPp4a,c2Pp4a,c3Pp4a,c4Pp4a,c4Pm3a,c3Pm3a,c2Pm3a,cPm3a,cPp3a,c2Pp3a,c3Pp3a,c4Pp3a,c4Pm2a,c3Pm2a,c2Pm2a,cPm2a,cPp2a,c2Pp2a,c3Pp2a,c4Pp2a,c4Pma,c3Pma,c2Pma,cPma,cPpa,c2Ppa,c3Ppa,c4Ppa;
  REAL8 s4Pm4a,s3Pm4a,s2Pm4a,sPm4a,sPp4a,s2Pp4a,s3Pp4a,s4Pp4a,s4Pm3a,s3Pm3a,s2Pm3a,sPm3a,sPp3a,s2Pp3a,s3Pp3a,s4Pp3a,s4Pm2a,s3Pm2a,s2Pm2a,sPm2a,sPp2a,s2Pp2a,s3Pp2a,s4Pp2a,s4Pma,s3Pma,s2Pma,sPma,sPpa,s2Ppa,s3Ppa,s4Ppa;
  //REAL8 Ssx, Ssy, Ssz;
  //REAL8 Sax, Say, Saz;
  LALSimInspiralInclAngle *an=NULL;
  COMPLEX16TimeSeries *h44 =XLALCreateCOMPLEX16TimeSeries( "h44", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h4m4=XLALCreateCOMPLEX16TimeSeries( "h4-4", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h43 =XLALCreateCOMPLEX16TimeSeries( "h43", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h4m3=XLALCreateCOMPLEX16TimeSeries( "h4-3", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h42 =XLALCreateCOMPLEX16TimeSeries( "h42", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h4m2=XLALCreateCOMPLEX16TimeSeries( "h4-2", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h41 =XLALCreateCOMPLEX16TimeSeries( "h41", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h4m1=XLALCreateCOMPLEX16TimeSeries( "h4-1", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);
  COMPLEX16TimeSeries *h40 =XLALCreateCOMPLEX16TimeSeries( "h40", &V->epoch, 0., V->deltaT, &lalStrainUnit, V->data->length);

  REAL8 amp44= 2.*eta*(m1+m2)*LAL_G_SI/pow(LAL_C_SI,2) / distance * sqrt(LAL_PI / 7.);
  REAL8 re4,re5,im4,im5;

  for (idx=0; idx<V->data->length; idx++) {
    v=V->data->data[idx];
    v2=v*v;
    v3=v2*v;
    Psi=Phi->data->data[idx];
    an=XLALSimInspiralComputeInclAngle(LNhz->data->data[idx]);
    //e2x=LNhy->data->data[idx]*e1z->data->data[idx]-LNhz->data->data[idx]*e1y->data->data[idx];
    e2y=LNhz->data->data[idx]*e1x->data->data[idx]-LNhx->data->data[idx]*e1z->data->data[idx];
    e2z=LNhx->data->data[idx]*e1y->data->data[idx]-LNhy->data->data[idx]*e1x->data->data[idx];
    //nx=e1x->data->data[idx]*cos(Psi)+e2x*sin(Psi);
    ny=e1y->data->data[idx]*cos(Psi)+e2y*sin(Psi);
    nz=e1z->data->data[idx]*cos(Psi)+e2z*sin(Psi);
    //lx=e2x*cos(Psi)-e1x->data->data[idx]*sin(Psi);
    ly=e2y*cos(Psi)-e1y->data->data[idx]*sin(Psi);
    lz=e2z*cos(Psi)-e1z->data->data[idx]*sin(Psi);
    if (an->si>EPS) {
      ca=LNhx->data->data[idx]/an->si;
      sa=LNhy->data->data[idx]/an->si;
      sP=nz/an->si;
      cP=lz/an->si;
    }
    else {
      ca=1.;
      sa=0.;
      cP=ca*ny-sa*ly;
      sP=-sa*ny-ca*ly;
    }
    c2a=2.*ca*ca-1.;
    s2a=2.*ca*sa;
    c3a=c2a*ca-s2a*sa;
    s3a=c2a*sa+s2a*ca;
    c4a=c2a*c2a-s2a*s2a;
    s4a=c2a*s2a+s2a*c2a;
    c2P=cP*cP-sP*sP;
    s2P=2.*cP*sP;
    c3P=c2P*cP-s2P*sP;
    s3P=c2P*sP+s2P*cP;
    c4P=c2P*c2P-s2P*s2P;
    s4P=c2P*s2P+s2P*c2P;
    c4Pp4a=c4P*c4a-s4P*s4a;
    c4Pm4a=c4P*c4a+s4P*s4a;
    s4Pp4a=s4P*c4a+c4P*s4a;
    s4Pm4a=s4P*c4a-c4P*s4a;
    c4Pp3a=c4P*c3a-s4P*s3a;
    c4Pm3a=c4P*c3a+s4P*s3a;
    s4Pp3a=s4P*c3a+c4P*s3a;
    s4Pm3a=s4P*c3a-c4P*s3a;
    c4Pp2a=c4P*c2a-s4P*s2a;
    c4Pm2a=c4P*c2a+s4P*s2a;
    s4Pp2a=s4P*c2a+c4P*s2a;
    s4Pm2a=s4P*c2a-c4P*s2a;
    c4Ppa=c4P*ca-s4P*sa;
    c4Pma=c4P*ca+s4P*sa;
    s4Ppa=s4P*ca+c4P*sa;
    s4Pma=s4P*ca-c4P*sa;
    c3Pp4a=c3P*c4a-s3P*s4a;
    c3Pm4a=c3P*c4a+s3P*s4a;
    s3Pp4a=s3P*c4a+c3P*s4a;
    s3Pm4a=s3P*c4a-c3P*s4a;
    c3Pp3a=c3P*c3a-s3P*s3a;
    c3Pm3a=c3P*c3a+s3P*s3a;
    s3Pp3a=s3P*c3a+c3P*s3a;
    s3Pm3a=s3P*c3a-c3P*s3a;
    c3Pp2a=c3P*c2a-s3P*s2a;
    c3Pm2a=c3P*c2a+s3P*s2a;
    s3Pp2a=s3P*c2a+c3P*s2a;
    s3Pm2a=s3P*c2a-c3P*s2a;
    c3Ppa=c3P*ca-s3P*sa;
    c3Pma=c3P*ca+s3P*sa;
    s3Ppa=s3P*ca+c3P*sa;
    s3Pma=s3P*ca-c3P*sa;
    c2Pp4a=c2P*c4a-s2P*s4a;
    c2Pm4a=c2P*c4a+s2P*s4a;
    s2Pp4a=s2P*c4a+c2P*s4a;
    s2Pm4a=s2P*c4a-c2P*s4a;
    c2Pp3a=c2P*c3a-s2P*s3a;
    c2Pm3a=c2P*c3a+s2P*s3a;
    s2Pp3a=s2P*c3a+c2P*s3a;
    s2Pm3a=s2P*c3a-c2P*s3a;
    c2Pp2a=c2P*c2a-s2P*s2a;
    c2Pm2a=c2P*c2a+s2P*s2a;
    s2Pp2a=s2P*c2a+c2P*s2a;
    s2Pm2a=s2P*c2a-c2P*s2a;
    c2Ppa=c2P*ca-s2P*sa;
    c2Pma=c2P*ca+s2P*sa;
    s2Ppa=s2P*ca+c2P*sa;
    s2Pma=s2P*ca-c2P*sa;
    cPp4a=cP*c4a-sP*s4a;
    cPm4a=cP*c4a+sP*s4a;
    sPp4a=sP*c4a+cP*s4a;
    sPm4a=sP*c4a-cP*s4a;
    cPp3a=cP*c3a-sP*s3a;
    cPm3a=cP*c3a+sP*s3a;
    sPp3a=sP*c3a+cP*s3a;
    sPm3a=sP*c3a-cP*s3a;
    cPp2a=cP*c2a-sP*s2a;
    cPm2a=cP*c2a+sP*s2a;
    sPp2a=sP*c2a+cP*s2a;
    sPm2a=sP*c2a-cP*s2a;
    cPpa=cP*ca-sP*sa;
    cPma=cP*ca+sP*sa;
    sPpa=sP*ca+cP*sa;
    sPma=sP*ca-cP*sa;

    //Sax=S1xtmp->data->data[idx]-S2xtmp->data->data[idx];
    //Ssx=S1xtmp->data->data[idx]+S2xtmp->data->data[idx];
    //Say=S1ytmp->data->data[idx]-S2ytmp->data->data[idx];
    //Ssy=S1ytmp->data->data[idx]+S2ytmp->data->data[idx];
    //Saz=S1ztmp->data->data[idx]-S2ztmp->data->data[idx];
    //Ssz=S1ztmp->data->data[idx]+S2ztmp->data->data[idx];

    //l=4, m=4,-4
    re4 = -8./9.*amp2*v2*(4.*c4Pm4a*an->siBy2Et + c2Pm4a*an->siBy2Sx*an->ciBy2Sq + c2Pp4a*an->siBy2Sq*an->ciBy2Sx + 4.*c4Pp4a*an->ciBy2Et);
    re5 = amp3*v3*dm*an->si*(-9./5.*(c3Pp4a*an->ciBy2Sx + c3Pm4a*an->siBy2Sx) -1./60.*an->siSq*( cPp4a*an->ciBy2Sq + cPm4a*an->siBy2Sq) );
    im4 = -8./9.*amp2*v2*(4.*s4Pm4a*an->siBy2Et + s2Pm4a*an->siBy2Sx*an->ciBy2Sq - s2Pp4a*an->siBy2Sq*an->ciBy2Sx - 4.*s4Pp4a*an->ciBy2Et);
    im5 = amp3*v3*dm*an->si*( 9./5.*(s3Pp4a*an->ciBy2Sx - s3Pm4a*an->siBy2Sx) +1./60.*an->siSq*( sPp4a*an->ciBy2Sq - sPm4a*an->siBy2Sq) );
    h44->data->data[idx] =amp44*v2*( re4+re5+I*( im4+im5));
    h4m4->data->data[idx]=amp44*v2*( re4-re5+I*(-im4+im5));

    //l=4, m=3,-3
    re4 = 8./9.*amp2*v2*an->si*( 4.*(c4Pm3a*an->siBy2Sx - c4Pp3a*an->ciBy2Sx) + 0.5*(an->siBy2Qu*(0.5+an->ci)*c2Pm3a - an->ciBy2Qu*(0.5-an->ci)*c2Pp3a ) );
    re5 = amp3*v3*dm* (9./5.*( an->ciBy2Sx*(-1.5+2.*an->ci)*c3Pp3a + an->siBy2Sx*(1.5+2.*an->ci)*c3Pm3a ) + 1./80.*an->siSq*( cPp3a*(an->ci+1./3.+2./3.*an->c2i) + cPm3a*(an->ci-1./3.-2./3.*an->c2i) ) );
    im4 = 8./9.*amp2*v2*an->si*( 4.*(s4Pm3a*an->siBy2Sx + s4Pp3a*an->ciBy2Sx) + 0.5*(an->siBy2Qu*(0.5+an->ci)*s2Pm3a + an->ciBy2Qu*(0.5-an->ci)*s2Pp3a ) );
    im5 = amp3*v3*dm* (9./5.*(-an->ciBy2Sx*(-1.5+2.*an->ci)*s3Pp3a + an->siBy2Sx*(1.5+2.*an->ci)*s3Pm3a ) + 1./80.*an->siSq*(-sPp3a*(an->ci+1./3.+2./3.*an->c2i) + sPm3a*(an->ci-1./3.-2./3.*an->c2i) ) );
    h43->data->data[idx] =amp44*sqrt(2.)*v2*( re4+re5+I*( im4+im5));
    h4m3->data->data[idx]=amp44*sqrt(2.)*v2*(-re4+re5+I*( im4-im5));

    //l=4, m=2,-2
    re4 = -1./63.*amp2*v2*( 112.*an->siSq*(c4Pm2a*an->siBy2Qu + c4Pp2a*an->ciBy2Qu ) + an->siBy2Qu*c2Pm2a*2.*(7.*an->c2i+14.*an->ci+9.) + an->ciBy2Qu*c2Pp2a*2.*(7.*an->c2i-14.*an->ci+9.) );
    re5 = amp3*v3*dm*an->si*( 9./10.*( c3Pp2a*an->ciBy2Qu*(-1.+2.*an->ci) - c3Pm2a*an->siBy2Qu*(1.+2.*an->ci) ) - 1./70.*( cPp2a*an->ciBy2Sq*(1.-7./6.*an->ci+7./6.*an->c2i) + cPm2a*an->siBy2Sq*(1.+7./6.*an->ci+7./6.*an->c2i) ) );
    im4 = -1./63.*amp2*v2*( 112.*an->siSq*(s4Pm2a*an->siBy2Qu - s4Pp2a*an->ciBy2Qu ) + an->siBy2Qu*s2Pm2a*2.*(7.*an->c2i+14.*an->ci+9.) - an->ciBy2Qu*s2Pp2a*2.*(7.*an->c2i-14.*an->ci+9.) );
    im5 = amp3*v3*dm*an->si*( 9./10.*( s3Pp2a*an->ciBy2Qu*( 1.-2.*an->ci) - s3Pm2a*an->siBy2Qu*(1.+2.*an->ci) ) + 1./70.*( sPp2a*an->ciBy2Sq*(1.-7./6.*an->ci+7./6.*an->c2i) - sPm2a*an->siBy2Sq*(1.+7./6.*an->ci+7./6.*an->c2i) ) );
    h42->data->data[idx] =amp44*sqrt(7.)*v2*( re4+re5+I*( im4+im5));
    h4m2->data->data[idx]=amp44*sqrt(7.)*v2*( re4-re5+I*(-im4+im5));

    //l=4, m=1,-1
    re4 = amp2*v2/9.*an->si*( 56.*an->siSq*( c4Pma*an->siBy2Sq - c4Ppa*an->ciBy2Sq) + 2.*( c2Pma*an->siBy2Sq*(3.+3.5*an->ci+3.5*an->c2i) - c2Ppa*an->ciBy2Sq*(3.-3.5*an->ci+3.5*an->c2i) ) );
    re5 = amp3*v3*dm*( 63./40.*an->siSq*( c3Ppa*an->ciBy2Sq*(-1.+4.*an->ci) + c3Pma*an->siBy2Sq*(1.+4.*an->ci) ) + 1./30.*( cPpa*an->ciBy2Sq*(-15./8.+15./4.*an->ci-21./8.*an->c2i+7./4.*an->c3i) + cPma*an->siBy2Sq*( 15./8.+15./4.*an->ci+21./8.*an->c2i+7./4.*an->c3i) ) );
    im4 = amp2*v2/9.*an->si*( 56.*an->siSq*( s4Pma*an->siBy2Sq + s4Ppa*an->ciBy2Sq) + 2.*( s2Pma*an->siBy2Sq*(3.+3.5*an->ci+3.5*an->c2i) + s2Ppa*an->ciBy2Sq*(3.-3.5*an->ci+3.5*an->c2i) ) );
    im5 = amp3*v3*dm*( 63./40.*an->siSq*( s3Ppa*an->ciBy2Sq*( 1.-4.*an->ci) + s3Pma*an->siBy2Sq*(1.+4.*an->ci) ) + 1./30.*( sPpa*an->ciBy2Sq*( 15./8.-15./4.*an->ci+21./8.*an->c2i-7./4.*an->c3i) + sPma*an->siBy2Sq*( 15./8.+15./4.*an->ci+21./8.*an->c2i+7./4.*an->c3i) ) );
    h41->data->data[idx] =amp44*sqrt(2./7.)*v2*( re4+re5+I*(im4+im5));
    h4m1->data->data[idx]=amp44*sqrt(2./7.)*v2*(-re4+re5+I*(im4-im5));

    //l=4, m=0
    re4 = -amp2*v2/18.*an->siSq*( 56.*an->siSq*c4P + (5.+7.*an->c2i)*c2P );
    re5 = 0.;
    im4 = 0.;
    im5 = -amp3*v3*dm*an->s2i/40.*( 63.*an->siSq*s3P + (1.+7.*an->c2i)/6.*sP );
    h40->data->data[idx]=amp44*sqrt(10./7.)*v2*(re4+re5+I*(im4+im5));

    LALFree(an);
    an=NULL;

  }

  /*XLALDestroyREAL8TimeSeries(S1xtmp);
  XLALDestroyREAL8TimeSeries(S1ytmp);
  XLALDestroyREAL8TimeSeries(S1ztmp);
  XLALDestroyREAL8TimeSeries(S2xtmp);
  XLALDestroyREAL8TimeSeries(S2ytmp);
  XLALDestroyREAL8TimeSeries(S2ztmp);*/

  *h4ms=XLALSphHarmTimeSeriesAddMode(*h4ms, h44, 4, 4);
  XLALDestroyCOMPLEX16TimeSeries(h44);
  *h4ms=XLALSphHarmTimeSeriesAddMode(*h4ms, h43, 4, 3);
  XLALDestroyCOMPLEX16TimeSeries(h43);
  *h4ms=XLALSphHarmTimeSeriesAddMode(*h4ms, h42, 4, 2);
  XLALDestroyCOMPLEX16TimeSeries(h42);
  *h4ms=XLALSphHarmTimeSeriesAddMode(*h4ms, h41, 4, 1);
  XLALDestroyCOMPLEX16TimeSeries(h41);
  *h4ms=XLALSphHarmTimeSeriesAddMode(*h4ms, h40, 4, 0);
  XLALDestroyCOMPLEX16TimeSeries(h40);
  *h4ms=XLALSphHarmTimeSeriesAddMode(*h4ms, h4m1, 4, -1);
  XLALDestroyCOMPLEX16TimeSeries(h4m1);
  *h4ms=XLALSphHarmTimeSeriesAddMode(*h4ms, h4m2, 4, -2);
  XLALDestroyCOMPLEX16TimeSeries(h4m2);
  *h4ms=XLALSphHarmTimeSeriesAddMode(*h4ms, h4m3, 4, -3);
  XLALDestroyCOMPLEX16TimeSeries(h4m3);
  *h4ms=XLALSphHarmTimeSeriesAddMode(*h4ms, h4m4, 4, -4);
  XLALDestroyCOMPLEX16TimeSeries(h4m4);

  return XLAL_SUCCESS;

}  /* End of XLALSimSpinInspiralPNMode4m*/

/** @} */
