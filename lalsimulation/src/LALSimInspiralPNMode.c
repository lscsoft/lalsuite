/*
 * Copyright (C) 2008 J. Creighton, Evan Ochsner, Chris Pankow
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
        case 5:
            re5 = - ((107.0/21.0) - (34.0/21.0)*nu)*LAL_PI;
            im5 = - 24.0*nu;
            im5log = - ((107.0/7.0) - (34.0/7.0)*nu)*2.0;
        case 4:
            re4 = - ((2173.0/1512.0) + (1069.0/216.0)*nu 
                    - (2047.0/1512.0)*nu2);
        case 3:
            re3 = 2.0*LAL_PI;
            im3log = 12.;
        case 2:
            re2 = - ((107.0/42.0) - (55.0/42.0)*nu);
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
    REAL8 fac = -8.0*sqrt(LAL_PI/5.0)*LAL_G_SI/LAL_C_SI/LAL_C_SI;
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
        case 4:
            re4 = LAL_PI;
            im4 = -(1.0/2.0) - 2.0*log(2.0);
            im4log = 6.0;
        case 3:
            re3 = - ((17.0/28.0) - (5.0/7.0)*nu);
        case 2:
        case 1:
            re1 = 1.;
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
        ans = crect(1., 0.);
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
        case 4:
            re4 = 3.0*LAL_PI;
            im4 = -(21.0/5.0) + 6.0*log(3.0/2.0);
            im4log = 18.0;
        case 3:
            re3 = - (4.0 - 2.0*nu);
        case 2:
        case 1:
            re1 = 1.;
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
        case 4:
            re4 = - ((193.0/90.0) - (145.0/18.0)*nu 
                    + (73.0/18.0)*nu2);
        case 3:
        case 2:
            re2 = nuterm;
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
        case 4:
            re4 = LAL_PI;
            im4 = -(7.0/5.0) - 2.0*log(2.0);
            im4log = 6.0;
        case 3:
            re3 = - ((8.0/3.0) + (2.0/3.0)*nu);
        case 2:
        case 1:
            re1 = 1.;
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
        case 5:
            re5 = 4.*LAL_PI*nuterm;
            im5 = -42./5. + (1193./40.) *nu + 8.*nuterm*log(2.);
            im5log = 24.*nuterm;
        case 4:
            re4 = 593./110. - (1273./66.)*nu + (175./22.)*nu2;
        case 3:
        case 2:
            re2 = nuterm;
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
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
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
        case 5:
            re5 = 2.*LAL_PI*nuterm;
            im5 = -21./5. + (84./5.) *nu + 8.*nuterm*log(2.);
            im5log = 12.*nuterm;
        case 4:
            re4 = 437./110. - (805./66.)*nu + (19./22.)*nu2;
        case 3:
        case 2:
            re2 = nuterm;
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
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
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
        ans = crect(1., 0.);
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
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
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
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
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
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
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
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
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
        case 4:
        case 3:
            re3 = 1. - 2.*nu;
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
        hlm->data->data[j] = crect(0., 0.);
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
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
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
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
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
        case 5:
        case 4:
            re4 = 1. - 5.*nu + 5.*nu2;
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
        hlm->data->data[j] = crect(1., 0.);
    }
    return hlm;
}

