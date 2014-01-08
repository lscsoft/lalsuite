/*
 * Copyright (C) 2008 J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria, D. Keppel, Evan Ochsner, C. Pankow
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

#include <complex.h>
#include <math.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimBlackHoleRingdown.h>

#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * (Twice) the highest known PN order of amplitude correction for
 * non-precessing binaries.
 */
#define MAX_NONPRECESSING_AMP_PN_ORDER 6

/**
 * (Twice) the highest known PN order of amplitude correction for
 * precessing binaries.
 */
#define MAX_PRECESSING_AMP_PN_ORDER 3

/* Macro functions to rotate the components of a vector about an axis */
#define ROTATEZ(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) - vy*sin(angle);\
	tmp2 = vx*sin(angle) + vy*cos(angle);\
	vx = tmp1;\
	vy = tmp2

#define ROTATEY(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) + vz*sin(angle);\
	tmp2 = - vx*sin(angle) + vz*cos(angle);\
	vx = tmp1;\
	vz = tmp2


/* Internal utility function to check all spin components are zero
   returns 1 if all spins zero, otherwise returns 0 */
static int checkSpinsZero(REAL8 s1x, REAL8 s1y, REAL8 s1z,
        REAL8 s2x, REAL8 s2y, REAL8 s2z);

static int checkSpinsZero(REAL8 s1x, REAL8 s1y, REAL8 s1z,
        REAL8 s2x, REAL8 s2y, REAL8 s2z)
{
    if( s1x != 0. || s1y != 0. || s1z != 0.
            || s2x != 0. || s2y != 0. || s2z != 0. )
        return 0;
    else
        return 1;
}

/* Internal utility function to check transverse spins are zero
   returns 1 if x and y components of spins are zero, otherwise returns 0 */
static int checkTransverseSpinsZero(REAL8 s1x, REAL8 s1y, REAL8 s2x, REAL8 s2y);

static int checkTransverseSpinsZero(REAL8 s1x, REAL8 s1y, REAL8 s2x, REAL8 s2y)
{
    if( s1x != 0. || s1y != 0. || s2x != 0. || s2y != 0. )
        return 0;
    else
        return 1;
}

/* Internal utility function to check tidal parameters are zero
   returns 1 if both tidal parameters zero, otherwise returns 0 */
static int checkTidesZero(REAL8 lambda1, REAL8 lambda2);

static int checkTidesZero(REAL8 lambda1, REAL8 lambda2)
{
    if( lambda1 != 0. || lambda2 != 0. )
        return 0;
    else
        return 1;
}

/**
 * Macro procedure for aborting if non-default LALSimInspiralWaveformFlags
 * struct was provided, but that approximant does not use the struct
 * and only has a single default use case.
 *
 * The ChooseWaveform functions will fail in such a case, so the user does not
 * think they are including features that are unavailable.
 *
 * All of the macros below will destroy the LALSimInspiralWaveformFlags struct,
 * print a specific warning and raise a general XLAL_ERROR for invalid argument.
 */
#define ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/**
 * Same as above macro, but returns a null pointer rather than XLAL_FAILURE int
 */
#define ABORT_NONDEFAULT_WAVEFORM_FLAGS_NULL(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.\n", __func__);\
	XLAL_ERROR_NULL(XLAL_EINVAL);\
	} while (0)

/**
 * Macro procedure for aborting if non-zero spins
 * given to a non-spinning approximant
 */
#define ABORT_NONZERO_SPINS(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-zero spins were given, but this is a non-spinning approximant.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/**
 * Macro procedure for aborting if non-zero transverse spin
 * components given to a non-precessing approximant
 */
#define ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-zero transverse spins were given, but this is a non-precessing approximant.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/**
 * Macro procedure for aborting if non-zero tidal parameters
 * given to an approximant with no tidal corrections
 */
#define ABORT_NONZERO_TIDES(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/**
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralSpinOrder is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_SPIN_ORDER(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/**
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralTidalOrder is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_TIDAL_ORDER(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralTidalOrder provided, but this approximant does not use that flag.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/**
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralInteraction is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_INTERACTION(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralInteraction provided, but this approximant does not use that flag.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/**
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralFrameAxis is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_FRAME_AXIS(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/**
 * Same as above macro, but returns a null pointer rather than XLAL_FAILURE int
 */
#define ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.\n", __func__);\
	XLAL_ERROR_NULL(XLAL_EINVAL);\
	} while (0)

/**
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralModesChoice is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_MODES_CHOICE(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/**
 * Same as above macro, but returns a null pointer rather than XLAL_FAILURE int
 */
#define ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.\n", __func__);\
	XLAL_ERROR_NULL(XLAL_EINVAL);\
	} while (0)

/**
 * Function that gives the default ending frequencies of the given approximant.
 *
 */
double XLALSimInspiralGetFinalFreq(
    REAL8 m1,                               /**< mass of companion 1 (kg) */
    REAL8 m2,                               /**< mass of companion 2 (kg) */
    const REAL8 S1x,                              /**< x-component of the dimensionless spin of object 1 */
    const REAL8 S1y,                              /**< y-component of the dimensionless spin of object 1 */
    const REAL8 S1z,                              /**< z-component of the dimensionless spin of object 1 */
    const REAL8 S2x,                              /**< x-component of the dimensionless spin of object 2 */
    const REAL8 S2y,                              /**< y-component of the dimensionless spin of object 2 */
    const REAL8 S2z,                              /**< z-component of the dimensionless spin of object 2 */
    Approximant approximant                 /**< post-Newtonian approximant to use for waveform production */
    )
{
    double  freq;   /* The return value */

    /* internal variables */
    /* we will use m1, m2 in solar masses */
    m1 /= LAL_MSUN_SI;
    m2 /= LAL_MSUN_SI;

    /* needed for Phenom */
    double chi;

    /* Needed for EOBNR */
    REAL8 spin1[3];
    REAL8 spin2[3];
    int     modeL;
    int     modeM;
    COMPLEX16Vector modefreqVec;
    COMPLEX16      modeFreq;

    switch (approximant)
    {
        /* non-spinning inspiral-only models */
        // CHECKME: do they really all use Schwarzschild ISCO? */
        case TaylorEt:
        case TaylorT1:
        case TaylorT2:
        case TaylorT3:
        case TaylorT4:
        case TaylorF2:
            /* Check that spins are zero */
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
            {
                XLALPrintError("Non-zero spins were given, but this is a non-spinning approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            }
        case TaylorF2RedSpin:
        case TaylorF2RedSpinTidal:
            /* Schwarzschild ISCO */
	        freq = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_MSUN_SI*LAL_G_SI);
            break;

        /* IMR models */
        /* EOBNR models all call the same code, just with different inputs */
        case EOBNRv2HM:
        case EOBNRv2:
        case SEOBNRv1:
            // FIXME: Probably shouldn't hard code the modes.
            if ( approximant == EOBNRv2HM )
            {
                modeL = 5;
                modeM = 5;
            }
            else
            {
                modeL = 2;
                modeM = 2;
            }
            if ( approximant == EOBNRv2 || approximant == EOBNRv2HM )
            {
                /* Check that spins are zero */
                if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                    XLALPrintError("Non-zero spins were given, but this is a non-spinning approximant.\n");
                    XLAL_ERROR(XLAL_EINVAL);
                spin1[0] = 0.; spin1[1] = 0.; spin1[2] = 0.;
                spin2[0] = 0.; spin2[1] = 0.; spin2[2] = 0.;
            }
            else
            {
                spin1[0] = S1x; spin1[1] = S1y; spin1[2] = S1z;
                spin2[0] = S2x; spin2[1] = S2y; spin2[2] = S2z;
            }

            modefreqVec.length = 1;
            modefreqVec.data   = &modeFreq;
            if ( XLALSimIMREOBGenerateQNMFreqV2( &modefreqVec, m1, m2, spin1, spin2, modeL, modeM, 1, approximant) != XLAL_SUCCESS )
            {
                XLAL_ERROR( XLAL_EFUNC );
            }

            freq = creal(modeFreq) / (2 * LAL_PI);
            break;

        case IMRPhenomA:
            /* Check that spins are zero */
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                XLALPrintError("Non-zero spins were given, but this is a non-spinning approximant.\n");
                XLAL_ERROR(XLAL_EINVAL);
            freq = XLALSimIMRPhenomAGetFinalFreq(m1, m2);
            break;

        case IMRPhenomB:
            chi = XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z);
            freq = XLALSimIMRPhenomBGetFinalFreq(m1, m2, chi);
            break;

        case IMRPhenomC:
            chi = XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z);
            freq = XLALSimIMRPhenomCGetFinalFreq(m1, m2, chi);
            break;


        // FIXME: Following I don't know how to calculate */
        /* Spinning inspiral-only time domain */
        case SpinTaylorT2:
        case SpinTaylorT4:
        case PhenSpinTaylor:
        /* Spinning with ringdown attachment */
        case PhenSpinTaylorRD:
        /* Spinning inspiral-only frequency domain */
        case SpinTaylorF2:
            XLALPrintError("I don't know how to calculate final freq. for this approximant, sorry!\n");
            XLAL_ERROR(XLAL_EINVAL);
            break;


        default:
            XLALPrintError("Unsupported approximant\n");
            XLAL_ERROR(XLAL_EINVAL);
    }

    return freq;
}
    

/**
 * Compute the polarizations from all the -2 spin-weighted spherical harmonic
 * modes stored in 'hlms'. Be sure that 'hlms' is the head of the linked list!
 *
 * The computation done is:
 * \f$hp(t) - i hc(t) = \sum_l \sum_m h_lm(t) -2Y_lm(iota,psi)\f$
 *
 * iota and psi are the inclination and polarization angle of the observer
 * relative to the source of GWs.
 */
int XLALSimInspiralPolarizationsFromSphHarmTimeSeries(
    REAL8TimeSeries **hp, /**< Plus polarization time series [returned] */
    REAL8TimeSeries **hc, /**< Cross polarization time series [returned] */
    SphHarmTimeSeries *hlms, /**< Head of linked list of waveform modes */
    REAL8 iota, /**< inclination of viewer to source frame (rad) */
    REAL8 psi /**< polarization angle (rad) */
    )
{
    int ret;
    SphHarmTimeSeries *ts = hlms;
    size_t length = ts->mode->data->length;
    // Destroy hp, hc TimeSeries if they already exist
    if( (*hp) ) XLALDestroyREAL8TimeSeries( *hp );
    if( (*hc) ) XLALDestroyREAL8TimeSeries( *hc );
    *hp = XLALCreateREAL8TimeSeries("hplus", &(ts->mode->epoch), ts->mode->f0,
                ts->mode->deltaT, &lalStrainUnit, length);
    *hc = XLALCreateREAL8TimeSeries("hplus", &(ts->mode->epoch), ts->mode->f0,
                ts->mode->deltaT, &lalStrainUnit, length);
    memset( (*hp)->data->data, 0, (*hp)->data->length*sizeof(REAL8) );
    memset( (*hc)->data->data, 0, (*hc)->data->length*sizeof(REAL8) );
    while (ts) { // Add the contribution from the current mode to hp, hx...
        // This function adds hlm(t) * Y_lm(incl,psi) to (h+ - i hx)(t)
        ret = XLALSimAddMode(*hp, *hc, ts->mode, iota, psi, ts->l, ts->m, 0);
        if( ret != XLAL_SUCCESS ) XLAL_ERROR(XLAL_EFUNC);
        ts = ts->next;
    }

    return XLAL_SUCCESS;
}

/**
 * Multiplies a mode h(l,m) by a spin-2 weighted spherical harmonic
 * to obtain hplus - i hcross, which is added to the time series.
 *
 * Implements the sum of a single term of Eq. (11) of:
 * Lawrence E. Kidder, \"Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit\", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 *
 * If sym is non-zero, symmetrically add the m and -m terms assuming
 * that \f$h(l,-m) = (-1)^l h(l,m)*\f$; see Eq. (78) ibid.
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
		hpc = Y * hmode->data->data[j];
		hplus->data->data[j] += creal(hpc);
		hcross->data->data[j] += -cimag(hpc);
	}
	if ( sym ) { /* equatorial symmetry: add in -m mode */
		Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m);
		if ( l % 2 ) /* l is odd */
			Y = -Y;
		for ( j = 0; j < hmode->data->length; ++j ) {
			COMPLEX16 hpc;
			hpc = Y * conj(hmode->data->data[j]);
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
 * Given time series for a binary's orbital dynamical variables,
 * construct the waveform polarizations h+ and hx as a sum of
 * -2 spin-weighted spherical harmonic modes, h_lm.
 * NB: Valid only for non-precessing systems!
 *
 * Implements Equation (11) of:
 * Lawrence E. Kidder, \"Using Full Information When Computing Modes of
 * Post-Newtonian Waveforms From Inspiralling Compact Binaries in Circular
 * Orbit\", Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
 */
int XLALSimInspiralPNPolarizationWaveformsFromModes(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform [returned] */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform [returned] */
	       	REAL8TimeSeries *v,       /**< post-Newtonian parameter */
	       	REAL8TimeSeries *phi,     /**< orbital phase */
	       	REAL8 v0,                 /**< tail-term gauge choice (default = 1) */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian order */
		)
{
	int l, m;
	LAL_CHECK_VALID_SERIES(v, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(phi, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(v, phi, XLAL_FAILURE);
	*hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &v->epoch, 0.0, v->deltaT, &lalStrainUnit, v->data->length );
	*hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &v->epoch, 0.0, v->deltaT, &lalStrainUnit, v->data->length );
	if ( ! hplus || ! hcross )
		XLAL_ERROR(XLAL_EFUNC);
	memset((*hplus)->data->data, 0, (*hplus)->data->length*sizeof(*(*hplus)->data->data));
	memset((*hcross)->data->data, 0, (*hcross)->data->length*sizeof(*(*hcross)->data->data));
	for ( l = 2; l <= LAL_PN_MODE_L_MAX; ++l ) {
		for ( m = 1; m <= l; ++m ) {
			COMPLEX16TimeSeries *hmode;
			hmode = XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(v, phi, v0, m1, m2, r, O, l, m);
			if ( ! hmode )
				XLAL_ERROR(XLAL_EFUNC);
			if ( XLALSimAddMode(*hplus, *hcross, hmode, i, 0.0, l, m, 1) < 0 )
				XLAL_ERROR(XLAL_EFUNC);
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
 * Implements Equations (8.8) - (8.10) of:
 * Luc Blanchet, Guillaume Faye, Bala R. Iyer and Siddhartha Sinha,
 * \"The third post-Newtonian gravitational wave polarisations
 * and associated spherical harmonic modes for inspiralling compact binaries
 * in quasi-circular orbits\", Class. Quant. Grav. 25 165003 (2008);
 * arXiv:0802.1249.
 * NB: Be sure to check arXiv:0802.1249v3, which corrects a typo!
 *
 * Note however, that we do not include the constant \"memory\" terms
 */
int XLALSimInspiralPNPolarizationWaveforms(
	REAL8TimeSeries **hplus,  /**< +-polarization waveform [returned] */
	REAL8TimeSeries **hcross, /**< x-polarization waveform [returned] */
	REAL8TimeSeries *V,       /**< post-Newtonian (PN) parameter */
	REAL8TimeSeries *Phi,     /**< orbital phase */
	REAL8 v0,                 /**< tail-term gauge choice (default = 1) */
	REAL8 m1,                 /**< mass of companion 1 (kg) */
	REAL8 m2,                 /**< mass of companion 2 (kg) */
	REAL8 r,                  /**< distance of source (m) */
	REAL8 i,                  /**< inclination of source (rad) */
	int ampO                  /**< twice PN order of the amplitude */
	)
{
    REAL8 M, eta, eta2, eta3, dm, dist, ampfac, phi, phiShift, v, v2, v3;
    REAL8 hp0, hp05, hp1, hp15, hp2, hp25, hp3;
    REAL8 hc0, hc05, hc1, hc15, hc2, hc25, hc3;
    REAL8 ci, si, ci2, ci4, ci6, ci8, si2, si3, si4, si5, si6;
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
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hplus)->data->data, 0, (*hplus)->data->length 
            * sizeof(*(*hplus)->data->data));
    memset((*hcross)->data->data, 0, (*hcross)->data->length 
            * sizeof(*(*hcross)->data->data));

    M = m1 + m2;
    eta = m1 * m2 / M / M; // symmetric mass ratio - '\nu' in the paper
    eta2 = eta*eta;	eta3 = eta2*eta;
    dm = (m1 - m2) / M; // frac. mass difference - \delta m/m in the paper
    dist = r / LAL_C_SI;   // r (m) / c (m/s) --> dist in units of seconds
    /* convert mass from kg to s, so ampfac ~ M/dist is dimensionless */
    ampfac = 2. * M * LAL_G_SI * pow(LAL_C_SI, -3) * eta / dist;
    
    /* 
     * cosines and sines of inclination between 
     * line of sight (N) and binary orbital angular momentum (L_N)
     */
    ci = cos(i);  	si = sin(i);
    ci2 = ci*ci;  ci4 = ci2*ci2;  ci6 = ci2*ci4;  ci8 = ci6*ci2;
    si2 = si*si;  si3 = si2*si;  si4 = si2*si2;  si5 = si*si4; si6 = si4*si2;

    /* loop over time steps and compute polarizations h+ and hx */
    len = V->data->length;
    for(idx = 0; idx < len; idx++)
    {
        /* Abbreviated names in lower case for time series at this sample */
        phi = Phi->data->data[idx]; 	v = V->data->data[idx];   
        v2 = v * v; 	v3 = v * v2;

        /* 
         * As explained in Blanchet et al, a phase shift can be applied 
         * to make log terms vanish which would appear in the amplitude 
         * at 1.5PN and 2.5PN orders. This shift is given in Eq. (8.8)
         * We apply the shift only for the PN orders which need it.
         */
        if( (ampO == -1) || ampO >= 5 )
            phiShift = 3.*v3*(1. - v2*eta/2.)*log( v2 / v0 / v0  );
        else if( ampO >= 3 )
            phiShift = 3.*v3*log( v2 / v0 / v0 );
        else
            phiShift = 0.;

        phi = phi - phiShift;

        /*
         * First set all h+/x coefficients to 0. Then use a switch to
         * set proper non-zero values up to order ampO. Note we
         * fall through the PN orders and break only after Newt. order
         */
        hp0 = hp05 = hp1 = hp15 = hp2 = hp25 = hp3 = 0.;
        hc0 = hc05 = hc1 = hc15 = hc2 = hc25 = hc3 = 0.;

        switch( ampO )
        {
            /* case LAL_PNORDER_THREE_POINT_FIVE: */
            case 7:
                XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %d\n", __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            case -1: // Highest known PN order - move if higher terms added!
            /* case LAL_PNORDER_THREE: */
            case 6:
		/* The reference had a typo in the 3PN terms and needed an errata.
		 * These should match arXiv:0802.1249v3, which has the fix. */
                hp3 = LAL_PI*dm*si*cos(phi)*(19./64. + ci2*5./16. - ci4/192.
                        + eta*(-19./96. + ci2*3./16. + ci4/96.)) + cos(2.*phi)
                        * (-465497./11025. + (LAL_GAMMA*856./105. 
                        - 2.*LAL_PI*LAL_PI/3. + log(16.*v2)*428./105.)
                        * (1. + ci2) - ci2*3561541./88200. - ci4*943./720.
                        + ci6*169./720. - ci8/360. + eta*(2209./360.
                        - LAL_PI*LAL_PI*41./96.*(1. + ci2) + ci2*2039./180.
                        + ci4*3311./720. - ci6*853./720. + ci8*7./360.)
                        + eta2*(12871./540. - ci2*1583./60. - ci4*145./108.
                        + ci6*56./45. - ci8*7./180.) + eta3*(-3277./810.
                        + ci2*19661./3240. - ci4*281./144. - ci6*73./720.
                        + ci8*7./360.)) + LAL_PI*dm*si*cos(3.*phi)*(-1971./128.
                        - ci2*135./16. + ci4*243./128. + eta*(567./64.
                        - ci2*81./16. - ci4*243./64.)) + si2*cos(4.*phi)
                        * (-2189./210. + ci2*1123./210. + ci4*56./9. 
                        - ci6*16./45. + eta*(6271./90. - ci2*1969./90.
                        - ci4*1432./45. + ci6*112./45.) + eta2*(-3007./27.
                        + ci2*3493./135. + ci4*1568./45. - ci6*224./45.)
                        + eta3*(161./6. - ci2*1921./90. - ci4*184./45.
                        + ci6*112./45.)) + dm*cos(5.*phi)*(LAL_PI*3125./384.
                        * si3*(1. + ci2)*(1. - 2.*eta)) + si4*cos(6.*phi)
                        * (1377./80. + ci2*891./80. - ci4*729./280. 
                        + eta*(-7857./80. - ci2*891./16. + ci4*729./40.)
                        + eta2*(567./4. + ci2*567./10. - ci4*729./20.)
                        + eta3*(-729./16. - ci2*243./80. + ci4*729./40.)) 
                        + cos(8.*phi)*(-1024./315.*si6*(1. + ci2)*(1. - 7.*eta 
                        + 14.*eta2 - 7.*eta3)) + dm*si*sin(phi)*(-2159./40320.
                        - log(2.)*19./32. + (-95./224. - log(2.)*5./8.)*ci2
                        + (181./13440. + log(2.)/96.)*ci4 + eta*(1369./160.
                        + log(2.)*19./48. + (-41./48. - log(2.)*3./8.)*ci2
                        + (-313./480. - log(2.)/48.)*ci4)) + sin(2.*phi)
                        * (-428.*LAL_PI/105.*(1. + ci2)) + dm*si*sin(3.*phi)
                        * (205119./8960. - log(3./2.)*1971./64. 
                        + (1917./224. - log(3./2.)*135./8.)*ci2
                        + (-43983./8960. + log(3./2.)*243./64.)*ci4 + eta
                        * (-54869./960. + log(3./2.)*567./32. 
                        + (-923./80. - log(3./2.)*81./8.)*ci2 
                        + (41851./2880. - log(3./2.)*243./32.)*ci4)) 
                        + dm*si3*(1. + ci2)*sin(5.*phi)*(-113125./5376. 
                        + log(5./2.)*3125./192. 
                        + eta*(17639./320. - log(5./2.)*3125./96.));
                hc3 = dm*si*ci*cos(phi)*(11617./20160. + log(2.)*21./16.
                        + (-251./2240. - log(2.)*5./48.)*ci2 
                        + eta*(-2419./240. - log(2.)*5./24.
                        + (727./240. + log(2.)*5./24.)*ci2)) + ci*cos(2.*phi)
                        * (LAL_PI*856./105.) + dm*si*ci*cos(3.*phi)
                        * (-36801./896. + log(3./2.)*1809./32.
                        + (65097./4480. - log(3./2.)*405./32.)*ci2 
                        + eta*(28445./288. - log(3./2.)*405./16. 
                        + (-7137./160. + log(3./2.)*405./16.)*ci2)) 
                        + dm*si3*ci*cos(5.*phi)*(113125./2688. 
                        - log(5./2.)*3125./96. + eta*(-17639./160. 
                        + log(5./2.)*3125./48.)) + LAL_PI*dm*si*ci*sin(phi)
                        * (21./32. - ci2*5./96. + eta*(-5./48. + ci2*5./48.))
                        + ci*sin(2.*phi)*(-3620761./44100. 
                        + LAL_GAMMA*1712./105. - 4.*LAL_PI*LAL_PI/3.
                        + log(16.*v2)*856./105. - ci2*3413./1260. 
                        + ci4*2909./2520. - ci6/45. + eta*(743./90. 
                        - 41.*LAL_PI*LAL_PI/48. + ci2*3391./180. 
                        - ci4*2287./360. + ci6*7./45.) + eta2*(7919./270.
                        - ci2*5426./135. + ci4*382./45. - ci6*14./45.) 
                        + eta3*(-6457./1620. + ci2*1109./180. - ci4*281./120.
                        + ci6*7./45.)) + LAL_PI*dm*si*ci*sin(3.*phi)
                        * (-1809./64. + ci2*405./64. + eta*(405./32. 
                        - ci2*405./32.)) + si2*ci*sin(4.*phi)*(-1781./105.
                        + ci2*1208./63. - ci4*64./45. + eta*(5207./45. 
                        - ci2*536./5. + ci4*448./45.) + eta2*(-24838./135.
                        + ci2*2224./15. - ci4*896./45.) + eta3*(1703./45.
                        - ci2*1976./45. + ci4*448./45.)) + dm*sin(5.*phi)
                        * (3125.*LAL_PI/192.*si3*ci*(1. - 2.*eta)) 
                        + si4*ci*sin(6.*phi)*(9153./280. - ci2*243./35. 
                        + eta*(-7371./40. + ci2*243./5.) + eta2*(1296./5. 
                        - ci2*486./5.) + eta3*(-3159./40. + ci2*243./5.))
                        + sin(8.*phi)*(-2048./315.*si6*ci*(1. - 7.*eta 
                        + 14.*eta2 - 7.*eta3));
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
                        + cos(5.*phi)*si*dm*(-108125./9216. + ci2*40625./9216. 
                        + ci4*83125./9216. - ci6*15625./9216. 
                        + eta*(8125./256. - ci2*40625./2304. - ci4*48125./2304.
                        + ci6*15625./2304.) + eta2*(-119375./9216. 
                        + ci2*40625./3072. + ci4*44375./9216. 
                        - ci6*15625./3072.)) + cos(7.*phi)*dm
                        * (117649./46080.*si5*(1. + ci2)*(1. - 4.*eta 
                        + 3.*eta2)) + sin(2.*phi)*(-9./5. + ci2*14./5. 
                        + ci4*7./5. + eta*(32. + ci2*56./5. - ci4*28./5.)) 
                        + sin(4.*phi)*si2*(1. + ci2)*(56./5. - 32.*log(2.)/3. 
                        + eta*(-1193./30. + 32.*log(2.)));
                /* below would have a constant memory term of si2*ci*eta*6./5. */
                hc25 = cos(2.*phi)*ci*(2. - ci2*22./5. + eta*(-282./5. 
                        + ci2*94./5.)) + cos(4.*phi)*ci*si2*(-112./5. 
                        + 64.*log(2.)/3. + eta*(1193./15. - 64.*log(2.)))
                        + sin(phi)*si*ci*dm*(-913./7680. + ci2*1891./11520. 
                        - ci4*7./4608. + eta*(1165./384. - ci2*235./576. 
                        + ci4*7./1152.) + eta2*(-1301./4608. + ci2*301./2304.
                        - ci4*7./1536.)) + sin(2.*phi)*LAL_PI*ci*(34./3. 
                        - ci2*8./3. + eta*(-20./3. + 8.*ci2)) 
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
                hp2 = cos(phi)*LAL_PI*si*dm*(-5./8. - ci2/8.) 
                        + cos(2.*phi)*(11./60. + ci2*33./10. + ci4*29./24. 
                        - ci6/24. + eta*(353./36. - 3.*ci2 - ci4*251./72. 
                        + ci6*5./24.) + eta2*(-49./12. + ci2*9./2. 
                        - ci4*7./24. - ci6*5./24.)) + cos(3.*phi)*LAL_PI*si*dm
                        * (27./8.*(1 + ci2)) + cos(4.*phi)*si2*2./15.*(59. 
                        + ci2*35. - ci4*8. 
                        - eta*5./3.*(131. + 59.*ci2 + 24.*ci4)
                        + eta2*5.*(21. - 3.*ci2 - 8.*ci4))
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
                        + sin(4.*phi)*ci*si2*4./15.*(55. - 12.*ci2 
                        - eta*5./3.*(119. - 36.*ci2)
                        + eta2*5.*(17. - 12.*ci2))
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
                        - 4.*LAL_PI*ci*sin(2.*phi) + sin(3.*phi)*si*ci*dm
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
                        + sin(4.*phi)*ci*si2*(-8./3.*(1. - 3.*eta));
            /*case LAL_PNORDER_HALF:*/
            case 1:
                hp05 = - si*dm*(cos(phi)*(5./8. + ci2/8.) 
                        - cos(3.*phi)*(9./8. + 9.*ci2/8.));
                hc05 = si*ci*dm*(-sin(phi)*3./4. + sin(3.*phi)*9./4.);
            case 0:
                /* below would have a constant memory term of -si2/96.*(17. + ci2) */
                hp0 = -(1. + ci2)*cos(2.*phi);
                hc0 = -2.*ci*sin(2.*phi);
                break;
            /*case LAL_PNORDER_NEWTONIAN:*/
            default:
                XLALPrintError("XLAL Error - %s: Invalid amp. PN order %s\n",
                        __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        } /* End switch on ampO */

        /* Fill the output polarization arrays */
        (*hplus)->data->data[idx] = ampfac * v2 * ( hp0 + v * ( hp05 
                + v * ( hp1 + v * ( hp15 + v * ( hp2 + v * ( hp25 + v * hp3 
                ) ) ) ) ) );
        (*hcross)->data->data[idx] = ampfac * v2 * ( hc0 + v * ( hc05 
                + v * ( hc1 + v * ( hc15 + v * ( hc2 + v * ( hc25 + v * hc3 
                ) ) ) ) ) );

    } /* end loop over time series samples idx */
    return XLAL_SUCCESS;
}

/**
 * Computes polarizations h+ and hx for a spinning, precessing binary
 * when provided time series of all the dynamical quantities.
 * Amplitude can be chosen between 1.5PN and Newtonian orders (inclusive).
 *
 * Based on K.G. Arun, Alesssandra Buonanno, Guillaume Faye and Evan Ochsner
 * \"Higher-order spin effects in the amplitude and phase of gravitational
 * waveforms emitted by inspiraling compact binaries: Ready-to-use
 * gravitational waveforms\", Phys Rev. D 79, 104023 (2009), arXiv:0810.5336
 *
 * HOWEVER, the formulae have been adapted to use the output of the so-called
 * \"Frameless\" convention for evolving precessing binary dynamics,
 * which is not susceptible to hitting coordinate singularities.
 *
 * FIXME: Clean up and commit Mathematica NB Showing correctness. Cite here.
 *
 * NOTE: The vectors MUST be given in the so-called radiation frame where
 * Z is the direction of propagation, X is the principal '+' axis and Y = Z x X
 * For different convention (Z is the direction of initial total angular
 * momentum, useful for GRB and comparison to NR, see XLALSimSpinInspiralGenerator())
 */
int XLALSimInspiralPrecessingPolarizationWaveforms(
	REAL8TimeSeries **hplus,  /**< +-polarization waveform [returned] */
	REAL8TimeSeries **hcross, /**< x-polarization waveform [returned] */
	REAL8TimeSeries *V,       /**< post-Newtonian parameter */
	REAL8TimeSeries *Phi,     /**< orbital phase */
	REAL8TimeSeries *S1x,	  /**< Spin1 vector x component */
	REAL8TimeSeries *S1y,	  /**< Spin1 vector y component */
	REAL8TimeSeries *S1z,	  /**< Spin1 vector z component */
	REAL8TimeSeries *S2x,	  /**< Spin2 vector x component */
	REAL8TimeSeries *S2y,	  /**< Spin2 vector y component */
	REAL8TimeSeries *S2z,	  /**< Spin2 vector z component */
	REAL8TimeSeries *LNhatx,  /**< unit orbital ang. mom. x comp. */
	REAL8TimeSeries *LNhaty,  /**< unit orbital ang. mom. y comp. */
	REAL8TimeSeries *LNhatz,  /**< unit orbital ang. mom. z comp. */
	REAL8TimeSeries *E1x,	  /**< orbital plane basis vector x comp. */
	REAL8TimeSeries *E1y,	  /**< orbital plane basis vector y comp. */
	REAL8TimeSeries *E1z,	  /**< orbital plane basis vector z comp. */
	REAL8 m1,                 /**< mass of companion 1 (kg) */
	REAL8 m2,                 /**< mass of companion 2 (kg) */
	REAL8 r,                  /**< distance of source (m) */
	REAL8 v0,                 /**< tail-term gauge choice (default = 1) */
	INT4 ampO	 	  /**< twice amp. post-Newtonian order */
	)
{
    REAL8 s1x, s1y, s1z, s2x, s2y, s2z, lnhx, lnhy, lnhz;
    REAL8 e1x, e1y, e1z, e2x, e2y, e2z, nx, ny, nz, lx, ly, lz;
    REAL8 nx2, ny2, nz2, nz3, lx2, ly2, lz2, lz3;
    REAL8 hplus0, hcross0, hplus05, hcross05, hplus1, hcross1;
    REAL8 hplus15, hcross15, hplusSpin1, hcrossSpin1;
    REAL8 hplusSpin15, hcrossSpin15, hplusTail15, hcrossTail15; 
    REAL8 M, eta, dm, phi, v, v2, dist, ampfac, logfac = 0.;
    INT4 idx, len;

    /* Macros to check time series vectors */
    LAL_CHECK_VALID_SERIES(V, 			XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(Phi, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1x, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1y, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1z, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2x, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2y, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2z, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatx, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhaty, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatz, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1x, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1y, 		XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1z, 		XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1x, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1y, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1z, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2x, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2y, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2z, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatx, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhaty, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatz, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1x, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1y, 	XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1z, 	XLAL_FAILURE);

    /* Allocate polarization vectors and set to 0 */
    *hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &V->epoch, 
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &V->epoch, 
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    if ( ! hplus || ! hcross )
        XLAL_ERROR(XLAL_EFUNC);
    memset((*hplus)->data->data, 0, 
            (*hplus)->data->length*sizeof(*(*hplus)->data->data));
    memset((*hcross)->data->data, 0, 
            (*hcross)->data->length*sizeof(*(*hcross)->data->data));

    M = m1 + m2;
    eta = m1 * m2 / M / M; // symmetric mass ratio - '\nu' in the paper
    dm = (m1 - m2) / M;    // frac. mass difference - \delta m/m in the paper
    dist = r / LAL_C_SI;   // r (m) / c (m/s) --> dist in units of seconds
    /* convert mass from kg to s, so ampfac ~ M/dist is dimensionless */
    ampfac = 2. * M * LAL_G_SI * pow(LAL_C_SI, -3) * eta / dist;
    
    /* loop over time steps and compute polarizations h+ and hx */
    len = V->data->length;
    for(idx = 0; idx < len; idx++)
    {
        /* Abbreviated names in lower case for time series at this sample */
        phi  = Phi->data->data[idx]; 	v = V->data->data[idx];     v2 = v * v;
        lnhx = LNhatx->data->data[idx]; e1x = E1x->data->data[idx];
        lnhy = LNhaty->data->data[idx];	e1y = E1y->data->data[idx];
        lnhz = LNhatz->data->data[idx];	e1z = E1z->data->data[idx];
        s1x  = S1x->data->data[idx];	s2x = S2x->data->data[idx];
        s1y  = S1y->data->data[idx];	s2y = S2y->data->data[idx];
        s1z  = S1z->data->data[idx];	s2z = S2z->data->data[idx];

        /* E2 = LNhat x E1 */
        e2x = lnhy*e1z - lnhz*e1y;
        e2y = lnhz*e1x - lnhx*e1z;
        e2z = lnhx*e1y - lnhy*e1x;

        /* Unit orbital separation vector */
        nx = e1x*cos(phi) + e2x*sin(phi);
        ny = e1y*cos(phi) + e2y*sin(phi);
        nz = e1z*cos(phi) + e2z*sin(phi);

        /* Unit inst. orbital velocity vector */
        lx = e2x*cos(phi) - e1x*sin(phi);
        ly = e2y*cos(phi) - e1y*sin(phi);
        lz = e2z*cos(phi) - e1z*sin(phi);

        /* Powers of vector components */
        nx2 = nx*nx;	ny2 = ny*ny;	nz2 = nz*nz;	nz3 = nz*nz2;
        lx2 = lx*lx;	ly2 = ly*ly;	lz2 = lz*lz;	lz3 = lz*lz2;

        /* 
         * First set all h+/x coefficients to 0. Then use a switch to
         * set proper non-zero values up to order ampO. Note we
         * fall through the PN orders and break only after Newt. order
         */
        hplus0 = hplus05 = hplus1 = hplus15 = hplusTail15 = 0.;
        hcross0 = hcross05 = hcross1 = hcross15 = hcrossTail15 = 0.;
        hplusSpin1 = hplusSpin15 = hcrossSpin1 = hcrossSpin15 = 0.;

        switch( ampO )
        {
            /*
             * case LAL_PNORDER_THREE_POINT_FIVE:
             * case LAL_PNORDER_THREE:
             * case LAL_PNORDER_TWO_POINT_FIVE:
             * case LAL_PNORDER_TWO:
             */
            case 7:
            case 6:
            case 5:
            case 4:
                XLALPrintError("XLAL Error - %s: Amp. corrections not known "
                        "to PN order %d, highest is %d\n", __func__, ampO,
                        MAX_PRECESSING_AMP_PN_ORDER );
                XLAL_ERROR(XLAL_EINVAL);
                break;
            case -1: /* Use highest known PN order - move if new orders added */
            /*case LAL_PNORDER_ONE_POINT_FIVE:*/
            case 3:
                /* 1.5PN non-spinning amp. corrections */
                hplus15 = (dm*(2*lx*nx*nz*(-95 + 90*lz2 - 65*nz2 
                        - 2*eta*(-9 + 90*lz2 - 65*nz2)) - 2*ly*ny*nz
                        * (-95 + 90*lz2 - 65*nz2 - 2*eta*(-9 + 90*lz2 - 65*nz2))
                        + 6*lx2*lz*(13 - 4*lz2 + 29*nz2 + eta*(-2 + 8*lz2 
                        - 58*nz2)) - 6*ly2*lz*(13 - 4*lz2 + 29*nz2 + eta
                        * (-2 + 8*lz2 - 58*nz2)) - lz*(nx2 - ny2)*(83 - 6*lz2 
                        + 111*nz2 + 6*eta*(-1 + 2*lz2 - 37*nz2))))/24.;
                hcross15 = (dm*(lz*(6*(19 - 4*eta)*lx*ly + (-101 + 12*eta)
                        * nx*ny) + (-149 + 36*eta) * (ly*nx + lx*ny)*nz 
                        + 6*(-3 + eta) * (2*lx*ly*lz - lz*nx*ny - 3*ly*nx*nz 
                        - 3*lx*ny*nz) + (1 - 2*eta) * (6*lz3*(-4*lx*ly + nx*ny) 
                        + 90*lz2*(ly*nx + lx*ny)*nz + 3*lz*(58*lx*ly 
                        - 37*nx*ny)*nz2 - 65*(ly*nx + lx*ny)*nz3)))/12.;
                /* 1.5PN spinning amp. corrections */
                hplusSpin15 = (6*lz*ny*s1x + 6*dm*lz*ny*s1x - 3*eta*lz*ny*s1x 
                        + 2*ly2*lnhy*s1y + 2*dm*ly2*lnhy*s1y 
                        + 2*eta*ly2*lnhy*s1y + 6*lz*nx*s1y + 6*dm*lz*nx*s1y 
                        - 3*eta*lz*nx*s1y + 8*lnhy*nx2*s1y + 8*dm*lnhy*nx2*s1y 
                        - eta*lnhy*nx2*s1y - 8*lnhy*ny2*s1y - 8*dm*lnhy*ny2*s1y
                        + eta*lnhy*ny2*s1y + 2*ly2*lnhz*s1z + 2*dm*ly2*lnhz*s1z
                        + 2*eta*ly2*lnhz*s1z - 6*ly*nx*s1z - 6*dm*ly*nx*s1z 
                        - 9*eta*ly*nx*s1z + 8*lnhz*nx2*s1z + 8*dm*lnhz*nx2*s1z 
                        - eta*lnhz*nx2*s1z - 8*lnhz*ny2*s1z - 8*dm*lnhz*ny2*s1z
                        + eta*lnhz*ny2*s1z + 6*lz*ny*s2x - 6*dm*lz*ny*s2x 
                        - 3*eta*lz*ny*s2x + lnhx*(2*ly2*((1 + dm + eta)*s1x 
                        + (1 - dm + eta)*s2x) + (nx2 - ny2)*((8 + 8*dm - eta)
                        * s1x - (-8 + 8*dm + eta)*s2x)) + 2*ly2*lnhy*s2y 
                        - 2*dm*ly2*lnhy*s2y + 2*eta*ly2*lnhy*s2y + 6*lz*nx*s2y 
                        - 6*dm*lz*nx*s2y - 3*eta*lz*nx*s2y + 8*lnhy*nx2*s2y 
                        - 8*dm*lnhy*nx2*s2y - eta*lnhy*nx2*s2y - 8*lnhy*ny2*s2y 
                        + 8*dm*lnhy*ny2*s2y + eta*lnhy*ny2*s2y + 2*ly2*lnhz*s2z 
                        - 2*dm*ly2*lnhz*s2z + 2*eta*ly2*lnhz*s2z - 6*ly*nx*s2z 
                        + 6*dm*ly*nx*s2z - 9*eta*ly*nx*s2z + 8*lnhz*nx2*s2z 
                        - 8*dm*lnhz*nx2*s2z - eta*lnhz*nx2*s2z - 8*lnhz*ny2*s2z 
                        + 8*dm*lnhz*ny2*s2z + eta*lnhz*ny2*s2z - 3*lx*ny 
                        * ((2 + 2*dm + 3*eta)*s1z + (2 - 2*dm + 3*eta)*s2z)
                        - 2*lx2*(lnhx*((1 + dm + eta)*s1x + (1 - dm + eta)*s2x) 
                        + lnhy*((1 + dm + eta)*s1y + (1 - dm + eta)*s2y) 
                        + lnhz*((1 + dm + eta)*s1z + (1 - dm + eta)*s2z)))/3.;
                hcrossSpin15 = (-3*lz*(nx*((2 + 2*dm - eta)*s1x 
                        - (-2 + 2*dm + eta)*s2x) + ny*((-2 - 2*dm + eta)*s1y 
                        + (-2 + 2*dm + eta)*s2y)) + ny*(-6*ly*s1z - 6*dm*ly*s1z 
                        - 9*eta*ly*s1z + 16*lnhz*nx*s1z + 16*dm*lnhz*nx*s1z 
                        - 2*eta*lnhz*nx*s1z + 2*lnhx*nx*((8 + 8*dm - eta)*s1x 
                        - (-8 + 8*dm + eta)*s2x) + 2*lnhy*nx*((8 + 8*dm - eta)
                        * s1y - (-8 + 8*dm + eta)*s2y) - 6*ly*s2z + 6*dm*ly*s2z 
                        - 9*eta*ly*s2z + 16*lnhz*nx*s2z - 16*dm*lnhz*nx*s2z 
                        - 2*eta*lnhz*nx*s2z) - lx*(4*lnhx*ly*((1 + dm + eta)*s1x
                        + (1 - dm + eta)*s2x) - 3*nx*((2 + 2*dm + 3*eta)*s1z 
                        + (2 - 2*dm + 3*eta)*s2z) + 4*ly*(lnhy*((1 + dm + eta)
                        * s1y + (1 - dm + eta)*s2y) + lnhz*((1 + dm + eta)*s1z 
                        + (1 - dm + eta)*s2z))))/3.;
                /* 1.5PN tail amp. corrections */
                logfac = log(v/v0);
                hplusTail15 = 2*((lx2 - ly2 - nx2 + ny2)*LAL_PI 
                        + 12*(lx*nx - ly*ny)*logfac);
                hcrossTail15 = 4*((lx*ly - nx*ny)*LAL_PI 
                        + 6*(ly*nx + lx*ny)*logfac);

            /*case LAL_PNORDER_ONE:*/
            case 2:
                /* 1PN non-spinning amp. corrections */
                hplus1 = (-13*lx2 + 13*ly2 + 6*lx2*lz2 - 6*ly2*lz2 
                        + 13*(nx2 - ny2) - 2*lz2*(nx2 - ny2) - 32*lx*lz*nx*nz 
                        + 32*ly*lz*ny*nz - 14*lx2*nz2 + 14*ly2*nz2 
                        + 10*(nx2 - ny2)*nz2)/6. + (eta*(lx2 - 18*lx2*lz2 
                        + 96*lx*lz*nx*nz - 96*ly*lz*ny*nz + 42*lx2*nz2 
                        + ly2*(-1 + 18*lz2 - 42*nz2) + (nx2 - ny2) 
                        * (-1 + 6*lz2 - 30*nz2)))/6.;
                hcross1 = (eta*(lx*ly - nx*ny - 6*(lz2*(3*lx*ly - nx*ny) 
                        - 8*lz*(ly*nx + lx*ny)*nz + (-7*lx*ly 
                        + 5*nx*ny)*nz2)))/3. + (-13*(lx*ly - nx*ny) 
                        + 2*(lz2*(3*lx*ly - nx*ny) - 8*lz*(ly*nx + lx*ny)*nz 
                        + (-7*lx*ly + 5*nx*ny)*nz2))/3.;
                /* 1PN spinning amp. corrections */
                hplusSpin1 = (-(ny*((1 + dm)*s1x + (-1 + dm)*s2x)) 
                        - nx*((1 + dm)*s1y + (-1 + dm)*s2y))/2.;
                hcrossSpin1 = (nx*((1 + dm)*s1x + (-1 + dm)*s2x) 
                        - ny*((1 + dm)*s1y + (-1 + dm)*s2y))/2.;

            /*case LAL_PNORDER_HALF:*/
            case 1:
                /* 0.5PN non-spinning amp. corrections */
                hplus05 = (dm*(-2*lx2*lz + 2*ly2*lz + lz*(nx2 - ny2) 
                        + 6*lx*nx*nz - 6*ly*ny*nz))/2.;
                hcross05 = dm*(-2*lx*ly*lz + lz*nx*ny 
					+ 3*ly*nx*nz + 3*lx*ny*nz);

            /*case LAL_PNORDER_NEWTONIAN:*/
            case 0:
                /* Newtonian order polarizations */
                hplus0 = lx2 - ly2 - nx2 + ny2;
                hcross0 = 2*lx*ly - 2*nx*ny;
                break;
            default: 
                XLALPrintError("XLAL Error - %s: Invalid amp. PN order %s\n", 
                        __func__, ampO );
                XLAL_ERROR(XLAL_EINVAL);
                break;
        } /* End switch on ampO */

        /* Fill the output polarization arrays */
        (*hplus)->data->data[idx] = ampfac * v2 * ( hplus0 
                + v * ( hplus05 + v * ( hplus1 + hplusSpin1 
                + v * ( hplus15 + hplusSpin15 + hplusTail15 ) ) ) );
        (*hcross)->data->data[idx] = ampfac * v2 * ( hcross0 
                + v * ( hcross05 + v * ( hcross1 + hcrossSpin1 
                + v * ( hcross15 + hcrossSpin15 + hcrossTail15 ) ) ) );
    } /* end of loop over time samples, idx */
    return XLAL_SUCCESS;
}

/**
 * Compute the physical template family "Q" vectors for a spinning, precessing
 * binary when provided time series of all the dynamical quantities.
 * These vectors always supplied to dominant order.
 *
 * Based on Pan, Buonanno, Chan and Vallisneri PRD69 104017, (see also theses
 * of Diego Fazi and Ian Harry)
 *
 * NOTE: The vectors MUST be given in the so-called radiation frame where
 * Z is the direction of propagation, X is the principal '+' axis and Y = Z x X
 */


int XLALSimInspiralPrecessingPTFQWaveforms(
        REAL8TimeSeries **Q1,     /**< PTF-Q1 waveform [returned] */
        REAL8TimeSeries **Q2,     /**< PTF-Q2 waveform [returned] */
        REAL8TimeSeries **Q3,     /**< PTF-Q2 waveform [returned] */
        REAL8TimeSeries **Q4,     /**< PTF-Q2 waveform [returned] */
        REAL8TimeSeries **Q5,     /**< PTF-Q2 waveform [returned] */
        REAL8TimeSeries *V,       /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi,     /**< orbital phase */
        REAL8TimeSeries *S1x,     /**< Spin1 vector x component */
        REAL8TimeSeries *S1y,     /**< Spin1 vector y component */
        REAL8TimeSeries *S1z,     /**< Spin1 vector z component */
        REAL8TimeSeries *S2x,     /**< Spin2 vector x component */
        REAL8TimeSeries *S2y,     /**< Spin2 vector y component */
        REAL8TimeSeries *S2z,     /**< Spin2 vector z component */
        REAL8TimeSeries *LNhatx,  /**< unit orbital ang. mom. x comp. */
        REAL8TimeSeries *LNhaty,  /**< unit orbital ang. mom. y comp. */
        REAL8TimeSeries *LNhatz,  /**< unit orbital ang. mom. z comp. */
        REAL8TimeSeries *E1x,     /**< orbital plane basis vector x comp. */
        REAL8TimeSeries *E1y,     /**< orbital plane basis vector y comp. */
        REAL8TimeSeries *E1z,     /**< orbital plane basis vector z comp. */
        REAL8 m1,                 /**< mass of companion 1 (kg) */
        REAL8 m2,                 /**< mass of companion 2 (kg) */
        REAL8 r                  /**< distance of source (m) */
        )
{
    REAL8 lnhx, lnhy, lnhz;
    REAL8 e1x, e1y, e1z, e2x, e2y, e2z, nx, ny, nz, lx, ly, lz;
    REAL8 nx2, ny2, nz2, lx2, ly2, lz2;
    REAL8 q1tmp, q2tmp, q3tmp, q4tmp, q5tmp;
    REAL8 M, eta, phi, v, v2, dist, ampfac;
    INT4 idx, len;
    REAL8 sqrt_three = pow(3,0.5);

    /* Macros to check time series vectors */
    LAL_CHECK_VALID_SERIES(V,                   XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(Phi,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1x,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1y,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1z,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2x,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2y,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2z,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatx,              XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhaty,              XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatz,              XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1x,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1y,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1z,                 XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1x,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1y,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1z,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2x,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2y,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2z,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatx, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhaty, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatz, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1x,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1y,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1z,    XLAL_FAILURE);


    /* Allocate polarization vectors and set to 0 */
    *Q1 = XLALCreateREAL8TimeSeries( "PTF_Q_1", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *Q2 = XLALCreateREAL8TimeSeries( "PTF_Q_2", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *Q3 = XLALCreateREAL8TimeSeries( "PTF_Q_3", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *Q4 = XLALCreateREAL8TimeSeries( "PTF_Q_4", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *Q5 = XLALCreateREAL8TimeSeries( "PTF_Q_5", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );

    if ( ! Q1 || ! Q2 || !Q3 || !Q4 || !Q5 )
        XLAL_ERROR(XLAL_EFUNC);
    memset((*Q1)->data->data, 0,
            (*Q1)->data->length*sizeof(*(*Q1)->data->data));
    memset((*Q2)->data->data, 0,
            (*Q2)->data->length*sizeof(*(*Q2)->data->data));
    memset((*Q3)->data->data, 0,
            (*Q3)->data->length*sizeof(*(*Q3)->data->data));
    memset((*Q4)->data->data, 0,
            (*Q4)->data->length*sizeof(*(*Q4)->data->data));
    memset((*Q5)->data->data, 0,
            (*Q5)->data->length*sizeof(*(*Q5)->data->data));

    M = m1 + m2;
    eta = m1 * m2 / M / M; // symmetric mass ratio - '\nu' in the paper
    dist = r / LAL_C_SI;   // r (m) / c (m/s) --> dist in units of seconds
    /* convert mass from kg to s, so ampfac ~ M/dist is dimensionless */
    ampfac = 2. * M * LAL_G_SI * pow(LAL_C_SI, -3) * eta / dist;

    /* loop over time steps and compute the Qi */
    len = V->data->length;
    for(idx = 0; idx < len; idx++)
    {
        /* Abbreviated names in lower case for time series at this sample */
        phi  = Phi->data->data[idx];    v = V->data->data[idx];     v2 = v * v;
        lnhx = LNhatx->data->data[idx]; e1x = E1x->data->data[idx];
        lnhy = LNhaty->data->data[idx]; e1y = E1y->data->data[idx];
        lnhz = LNhatz->data->data[idx]; e1z = E1z->data->data[idx];

        /* E2 = LNhat x E1 */
        e2x = lnhy*e1z - lnhz*e1y;
        e2y = lnhz*e1x - lnhx*e1z;
        e2z = lnhx*e1y - lnhy*e1x;

        /* Unit orbital separation vector */
        nx = e1x*cos(phi) + e2x*sin(phi);
        ny = e1y*cos(phi) + e2y*sin(phi);
        nz = e1z*cos(phi) + e2z*sin(phi);

        /* Unit inst. orbital velocity vector */
        lx = e2x*cos(phi) - e1x*sin(phi);
        ly = e2y*cos(phi) - e1y*sin(phi);
        lz = e2z*cos(phi) - e1z*sin(phi);

        /* Powers of vector components */
        nx2 = nx*nx;    ny2 = ny*ny;    nz2 = nz*nz;
        lx2 = lx*lx;    ly2 = ly*ly;    lz2 = lz*lz;

        /* 
         * NOTE: For PTF waveforms, we must use only the dominant amplitude
         * The Q values are computed from equations 13,14,17, 46 + 47 in PBCV or
         * more simply from equations (3.10) in Diego Fazi's thesis.
         * Note that Q5 is simplified from that shown in Fazi's thsis
         * by using traceless condition. As demonstrated in (6.29)
         * in Ian Harry's thesis.
         */

        q1tmp = lx2 - ly2 - nx2 + ny2;
        q2tmp = 2*lx*ly - 2*nx*ny;
        q3tmp = 2*lx*lz - 2*nx*nz;
        q4tmp = 2*ly*lz - 2*ny*nz;
        q5tmp = sqrt_three * (nz2 - lz2);

        /* Fill the output vectors */
        (*Q1)->data->data[idx] = ampfac * v2 * q1tmp;
        (*Q2)->data->data[idx] = ampfac * v2 * q2tmp;
        (*Q3)->data->data[idx] = ampfac * v2 * q3tmp;
        (*Q4)->data->data[idx] = ampfac * v2 * q4tmp;
        (*Q5)->data->data[idx] = ampfac * v2 * q5tmp;
    }
    return XLAL_SUCCESS;
}

/**
 * Function to specify the desired orientation of a precessing binary in terms
 * of several angles and then compute the vector components in the so-called
 * \"radiation frame\" (with the z-axis along the direction of propagation) as
 * needed to specify binary configuration for ChooseTDWaveform.
 *
 * Input:
 * thetaJN is the inclination between total angular momentum (J) and the
 * direction of propagation (N)
 * theta1 and theta2 are the inclinations of S1 and S2
 * measured from the Newtonian orbital angular momentum (L_N)
 * phi12 is the difference in azimuthal angles of S1 and S2.
 * chi1, chi2 are the dimensionless spin magnitudes ( \f$0 \le chi1,2 \le 1\f$)
 * phiJL is the azimuthal angle of L_N on its cone about J.
 * m1, m2, f_ref are the component masses and reference GW frequency,
 * they are needed to compute the magnitude of L_N, and thus J.
 *
 * Output:
 * incl - inclination angle of L_N relative to N
 * x, y, z components of E1 (unit vector in the initial orbital plane)
 * x, y, z components S1 and S2 (unit spin vectors times their
 * dimensionless spin magnitudes - i.e. they have unit magnitude for
 * extremal BHs and smaller magnitude for slower spins).
 *
 * NOTE: Here the \"total\" angular momentum is computed as
 * J = L_N + S1 + S2
 * where L_N is the Newtonian orbital angular momentum. In fact, there are
 * PN corrections to L which contribute to J that are NOT ACCOUNTED FOR
 * in this function. This is done so the function does not need to know about
 * the PN order of the system and to avoid subtleties with spin-orbit
 * contributions to L. Also, it is believed that the difference in Jhat
 * with or without these PN corrections to L is quite small.
 *
 * NOTE: fRef = 0 is not a valid choice. If you will pass fRef=0 into
 * ChooseWaveform, then here pass in f_min, the starting GW frequency
 *
 * The various rotations in this transformation are described in more detail
 * in a Mathematica notebook available here:
 * https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/Waveforms/TransformPrecessingInitialConditions
 */
int XLALSimInspiralTransformPrecessingInitialConditions(
		REAL8 *incl,	/**< Inclination angle of L_N (returned) */
		REAL8 *S1x,	/**< S1 x component (returned) */
		REAL8 *S1y,	/**< S1 y component (returned) */
		REAL8 *S1z,	/**< S1 z component (returned) */
		REAL8 *S2x,	/**< S2 x component (returned) */
		REAL8 *S2y,	/**< S2 y component (returned) */
		REAL8 *S2z,	/**< S2 z component (returned) */
		REAL8 thetaJN, 	/**< zenith angle between J and N (rad) */
		REAL8 phiJL,  	/**< azimuthal angle of L_N on its cone about J (rad) */
		REAL8 theta1,  	/**< zenith angle between S1 and LNhat (rad) */
		REAL8 theta2,  	/**< zenith angle between S2 and LNhat (rad) */
		REAL8 phi12,  	/**< difference in azimuthal angle btwn S1, S2 (rad) */
		REAL8 chi1,	/**< dimensionless spin of body 1 */
		REAL8 chi2,	/**< dimensionless spin of body 2 */
		REAL8 m1,	/**< mass of body 1 (kg) */
		REAL8 m2,	/**< mass of body 2 (kg) */
		REAL8 fRef	/**< reference GW frequency (Hz) */
		)
{
	/* Check that fRef is sane */
	if( fRef == 0. )
	{
		XLALPrintError("XLAL Error - %s: fRef=0 is invalid. Please pass in the starting GW frequency instead.\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
	}

	REAL8 omega0, M, eta, theta0, phi0, Jnorm, tmp1, tmp2;
	REAL8 Jhatx, Jhaty, Jhatz, LNhx, LNhy, LNhz, Jx, Jy, Jz, LNmag;
	REAL8 s1hatx, s1haty, s1hatz, s2hatx, s2haty, s2hatz;
	REAL8 s1x, s1y, s1z, s2x, s2y, s2z;

	/* Starting frame: LNhat is along the z-axis and the unit
	 * spin vectors are defined from the angles relative to LNhat.
	 * Note that we put s1hat in the x-z plane, and phi12
	 * sets the azimuthal angle of s2hat measured from the x-axis.
	 */
	LNhx = 0.;
	LNhy = 0.;
	LNhz = 1.;
	s1hatx = sin(theta1);
	s1haty = 0.;
	s1hatz = cos(theta1);
	s2hatx = sin(theta2) * cos(phi12);
	s2haty = sin(theta2) * sin(phi12);
	s2hatz = cos(theta2);

	/* Define several internal variables needed for magnitudes */
	omega0 = LAL_PI * fRef; // orbital angular frequency at reference point
	m1 *= LAL_G_SI / LAL_C_SI / LAL_C_SI / LAL_C_SI;
	m2 *= LAL_G_SI / LAL_C_SI / LAL_C_SI / LAL_C_SI;
	M = m1 + m2;
	eta = m1 * m2 / M / M;
	
	/* Define S1, S2, J with proper magnitudes */
	LNmag = pow(M, 5./3.) * eta * pow(omega0, -1./3.);
	s1x = m1 * m1 * chi1 * s1hatx;
	s1y = m1 * m1 * chi1 * s1haty;
	s1z = m1 * m1 * chi1 * s1hatz;
	s2x = m2 * m2 * chi2 * s2hatx;
	s2y = m2 * m2 * chi2 * s2haty;
	s2z = m2 * m2 * chi2 * s2hatz;
	Jx = s1x + s2x;
	Jy = s1y + s2y;
	Jz = LNmag * LNhz + s1z + s2z;

	/* Normalize J to Jhat, find it's angles in starting frame */
	Jnorm = sqrt( Jx*Jx + Jy*Jy + Jz*Jz);
	Jhatx = Jx / Jnorm;
	Jhaty = Jy / Jnorm;
	Jhatz = Jz / Jnorm;
	theta0 = acos(Jhatz);
	phi0 = atan2(Jhaty, Jhatx);

	/* Rotation 1: Rotate about z-axis by -phi0 to put Jhat in x-z plane */
	ROTATEZ(-phi0, LNhx, LNhy, LNhz);
	ROTATEZ(-phi0, s1hatx, s1haty, s1hatz);
	ROTATEZ(-phi0, s2hatx, s2haty, s2hatz);
	ROTATEZ(-phi0, Jhatx, Jhaty, Jhatz);

	/* Rotation 2: Rotate about new y-axis by -theta0
	 * to put Jhat along z-axis
	 */
	ROTATEY(-theta0, LNhx, LNhy, LNhz);
	ROTATEY(-theta0, s1hatx, s1haty, s1hatz);
	ROTATEY(-theta0, s2hatx, s2haty, s2hatz);
	ROTATEY(-theta0, Jhatx, Jhaty, Jhatz);

	/* Rotation 3: Rotate about new z-axis by phiJL to put L at desired
	 * azimuth about J. Note that is currently in x-z plane towards -x
	 * (i.e. azimuth=pi). Hence we rotate about z by phiJL - LAL_PI
	 */
	ROTATEZ(phiJL - LAL_PI, LNhx, LNhy, LNhz);
	ROTATEZ(phiJL - LAL_PI, s1hatx, s1haty, s1hatz);
	ROTATEZ(phiJL - LAL_PI, s2hatx, s2haty, s2hatz);
	ROTATEZ(phiJL - LAL_PI, Jhatx, Jhaty, Jhatz);

	/* Rotation 4: Let N be in x-z plane, inclined from J by thetaJN.
	 * We don't need to explicitly construct it, but rotating the system
	 * about the y-axis by - thetaJN will bring N onto the z-axis.
	 */
	ROTATEY(-thetaJN, LNhx, LNhy, LNhz);
	ROTATEY(-thetaJN, s1hatx, s1haty, s1hatz);
	ROTATEY(-thetaJN, s2hatx, s2haty, s2hatz);
	ROTATEY(-thetaJN, Jhatx, Jhaty, Jhatz);

	/* Rotation 5: We already have N along z. To complete the
	 * transformation into the radiation frame we rotate s.t.
	 * LNhat is in the x-z plane. Thus, if LNhat has azimuth phi0,
	 * we rotate about the z-axis by -phi0.
	 */
	phi0 = atan2(LNhy, LNhx);
	ROTATEZ(-phi0, LNhx, LNhy, LNhz);
	ROTATEZ(-phi0, s1hatx, s1haty, s1hatz);
	ROTATEZ(-phi0, s2hatx, s2haty, s2hatz);
	ROTATEZ(-phi0, Jhatx, Jhaty, Jhatz);

	/* We have completed the transformation. Now find the inclination
	 * of LN relative to N (the final z-axis). */
	*incl = acos(LNhz);
	
	/* Multiply spin unit vectors by chi magnitude (but NOT m_i^2) */
	s1hatx *= chi1;
	s1haty *= chi1;
	s1hatz *= chi1;
	s2hatx *= chi2;
	s2haty *= chi2;
	s2hatz *= chi2;

	/* Set pointers to rotated spin vectors */
	*S1x = s1hatx;
	*S1y = s1haty;
	*S1z = s1hatz;
	*S2x = s2hatx;
	*S2y = s2haty;
	*S2z = s2hatz;

	return XLAL_SUCCESS;
}

/**
 * DEPRECATED: USE XLALSimInspiralChooseTDWaveform() INSTEAD
 *
 * Chooses between different approximants when requesting a waveform to be generated
 * For spinning waveforms, all known spin effects up to given PN order are included
 *
 * The parameters passed must be in SI units.
 */
int XLALSimInspiralChooseWaveform(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 i,                                    /**< inclination of source (rad) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    XLAL_PRINT_DEPRECATION_WARNING("XLALSimInspiralChooseTDWaveform");

    return XLALSimInspiralChooseTDWaveform(hplus, hcross, phiRef, deltaT, m1,
        m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i, lambda1, lambda2,
        waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
}

/**
 * Chooses between different approximants when requesting a waveform to be generated
 * For spinning waveforms, all known spin effects up to given PN order are included
 * Returns the waveform in the time domain.
 *
 * The parameters passed must be in SI units.
 */
int XLALSimInspiralChooseTDWaveform(
    REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 S1x,                                  /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                                  /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                                  /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                                  /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                                  /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                                  /**< z-component of the dimensionless spin of object 2 */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 i,                                    /**< inclination of source (rad) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    REAL8 LNhatx, LNhaty, LNhatz, E1x, E1y, E1z;
    int ret;
    /* N.B. the quadrupole of a spinning compact body labeled by A is 
     * Q_A = - quadparam_A chi_A^2 m_A^3 (see gr-qc/9709032)
     * where quadparam = 1 for BH ~= 4-8 for NS.
     * This affects the quadrupole-monopole interaction.
     * For now, hardcode quadparam1,2 = 1.
     * Will later add ability to set via LALSimInspiralTestGRParam
     */
    REAL8 v0 = 1., quadparam1 = 1., quadparam2 = 1.;

    /* General sanity checks that will abort */

    if( nonGRparams && approximant != PhenSpinTaylor && approximant != PhenSpinTaylorRD )
    {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* General sanity check the input parameters - only give warnings! */
    if( deltaT > 1. )
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaT = %e requested.\nPerhaps sample rate and time step size were swapped?\n", __func__, deltaT);
    if( deltaT < 1./16385. )
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaT = %e requested.\nCheck for errors, this could create very large time series.\n", __func__, deltaT);
    if( m1 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m1, m1/LAL_MSUN_SI);
    if( m2 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m2, m2/LAL_MSUN_SI);
    if( m1 + m2 > 1000. * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested.\nSignal not likely to be in band of ground-based detectors.\n", __func__, m1+m2, (m1+m2)/LAL_MSUN_SI);
    if( S1x*S1x + S1y*S1y + S1z*S1z > 1.000001 )
        XLALPrintWarning("XLAL Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, S1x, S1y, S1z);
    if( S2x*S2x + S2y*S2y + S2z*S2z > 1.000001 )
        XLALPrintWarning("XLAL Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested.\nAre you sure you want to violate the Kerr bound?\n", __func__, S2x, S2y, S2z);
    if( f_min < 1. )
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested.\nCheck for errors, this could create a very long waveform.\n", __func__, f_min);
    if( f_min > 40.000001 )
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested.\nCheck for errors, the signal will start in band.\n", __func__, f_min);

    switch (approximant)
    {
        /* non-spinning inspiral-only models */
        case TaylorEt:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorEtPNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, r, i, amplitudeO, phaseO);
            break;

        case TaylorT1:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                ABORT_NONDEFAULT_SPIN_ORDER(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, i, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO, phaseO);
            break;

        case TaylorT2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                ABORT_NONDEFAULT_SPIN_ORDER(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT2PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, i, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO, phaseO);
            break;

        case TaylorT3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                ABORT_NONDEFAULT_SPIN_ORDER(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, i, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO, phaseO);
            break;

        case TaylorT4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !XLALSimInspiralSpinOrderIsDefault( XLALSimInspiralGetSpinOrder(waveFlags)) )
                ABORT_NONDEFAULT_SPIN_ORDER(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorT4PNGenerator(hplus, hcross, phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, i, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO, phaseO);
            break;

        /* non-spinning inspiral-merger-ringdown models */
        case IMRPhenomA:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // NB: f_max = 0 will generate up to the ringdown cut-off frequency
            ret = XLALSimIMRPhenomAGenerateTD(hplus, hcross, phiRef, deltaT,
                    m1, m2, f_min, 0., r, i);
            break;

        case EOBNRv2HM:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // FIXME: need to create a function to take in different modes or produce an error if all modes not given
            ret = XLALSimIMREOBNRv2AllModes(hplus, hcross, phiRef, deltaT,
                    m1, m2, f_min, r, i);
            break;

        case EOBNRv2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            ret = XLALSimIMREOBNRv2DominantMode(hplus, hcross, phiRef, deltaT,
                    m1, m2, f_min, r, i);
            break;

        /* spinning inspiral-only models */
        case SpinTaylorT2:
            /* Waveform-specific sanity checks */
            /* Sanity check unused fields of waveFlags */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            E1x = cos(i);
            E1y = 0.;
            E1z = - sin(i);
            /* Maximum PN amplitude order for precessing waveforms is 
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ? 
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT2(hplus, hcross, phiRef, v0, deltaT,
                    m1, m2, f_min, f_ref, r, S1x, S1y, S1z, S2x, S2y, S2z,
                    LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2,
                    quadparam1, quadparam2,
                    XLALSimInspiralGetSpinOrder(waveFlags),
                    XLALSimInspiralGetTidalOrder(waveFlags),
                    phaseO, amplitudeO);
            break;

        // need to make a consistent choice for SpinTaylorT4 and PSpinInspiralRD waveform inputs
        // proposal: TotalJ frame of PSpinInspiralRD
        // inclination denotes the angle between the view direction 
        // and J (J is constant during the evolution, J//z, both N and initial 
        // L are in the x-z plane) and the spin coordinates are given wrt 
        // initial ** L **.
        case SpinTaylorT4:
            /* Waveform-specific sanity checks */
            /* Sanity check unused fields of waveFlags */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            E1x = cos(i);
            E1y = 0.;
            E1z = - sin(i);
            /* Maximum PN amplitude order for precessing waveforms is 
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ? 
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER;
            /* Call the waveform driver routine */
            ret = XLALSimInspiralSpinTaylorT4(hplus, hcross, phiRef, v0, deltaT,
                    m1, m2, f_min, f_ref, r, S1x, S1y, S1z, S2x, S2y, S2z,
                    LNhatx, LNhaty, LNhatz, E1x, E1y, E1z, lambda1, lambda2,
                    quadparam1, quadparam2,
                    XLALSimInspiralGetSpinOrder(waveFlags),
                    XLALSimInspiralGetTidalOrder(waveFlags),
                    phaseO, amplitudeO);
            break;

        /* spin aligned inspiral-merger-ringdown models */
        case IMRPhenomB:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            // NB: f_max = 0 will generate up to the ringdown cut-off frequency
            ret = XLALSimIMRPhenomBGenerateTD(hplus, hcross, phiRef, deltaT,
                    m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, 0., r, i);
            break;

        case PhenSpinTaylor:
            /* Waveform-specific sanity checks */
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimSpinInspiralGenerator(hplus, hcross, phiRef,
					       deltaT, m1, m2, f_min, f_ref, r, i, S1x, S1y, S1z, S2x, S2y, S2z,
					       phaseO, amplitudeO, waveFlags, nonGRparams);
            break;

        case PhenSpinTaylorRD:
            /* Waveform-specific sanity checks */
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at the start.\n", __func__);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenSpinInspiralRDGenerator(hplus, hcross, phiRef,
							deltaT, m1, m2, f_min, f_ref, r, i, S1x, S1y, S1z, S2x, S2y, S2z,
							phaseO, amplitudeO,  waveFlags, nonGRparams);
            break;

        case SEOBNRv1:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            if( f_ref != 0.)
                XLALPrintWarning("XLAL Warning - %s: This approximant does use f_ref. The reference phase will be defined at coalescence.\n", __func__);
            /* Call the waveform driver routine */
            ret = XLALSimIMRSpinAlignedEOBWaveform(hplus, hcross, phiRef, 
                    deltaT, m1, m2, f_min, r, i, S1z, S2z);
            break;

        default:
            XLALPrintError("TD version of approximant not implemented in lalsimulation\n");
            XLAL_ERROR(XLAL_EINVAL);
    }

    if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);

    return ret;
}

/**
 * Chooses between different approximants when requesting a waveform to be generated
 * For spinning waveforms, all known spin effects up to given PN order are included
 * Returns the waveform in the frequency domain.
 */
int XLALSimInspiralChooseFDWaveform(
    COMPLEX16FrequencySeries **hptilde,     /**< FD plus polarization */
    COMPLEX16FrequencySeries **hctilde,     /**< FD cross polarization */
    REAL8 phiRef,                           /**< reference orbital phase (rad) */
    REAL8 deltaF,                           /**< sampling interval (Hz) */
    REAL8 m1,                               /**< mass of companion 1 (kg) */
    REAL8 m2,                               /**< mass of companion 2 (kg) */
    REAL8 S1x,                              /**< x-component of the dimensionless spin of object 1 */
    REAL8 S1y,                              /**< y-component of the dimensionless spin of object 1 */
    REAL8 S1z,                              /**< z-component of the dimensionless spin of object 1 */
    REAL8 S2x,                              /**< x-component of the dimensionless spin of object 2 */
    REAL8 S2y,                              /**< y-component of the dimensionless spin of object 2 */
    REAL8 S2z,                              /**< z-component of the dimensionless spin of object 2 */
    REAL8 f_min,                            /**< starting GW frequency (Hz) */
    REAL8 f_max,                            /**< ending GW frequency (Hz) */
    REAL8 r,                                /**< distance of source (m) */
    REAL8 i,                                /**< inclination of source (rad) */
    REAL8 lambda1,                          /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                          /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags, /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, /**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                         /**< twice post-Newtonian amplitude order */
    int phaseO,                             /**< twice post-Newtonian order */
    Approximant approximant                 /**< post-Newtonian approximant to use for waveform production */
    )
{
    REAL8 LNhatx, LNhaty, LNhatz;
    int ret;
    unsigned int j;
    REAL8 pfac, cfac;

    /* General sanity checks that will abort */
    /*
     * If non-GR approximants are added, change the below to
     * if( nonGRparams && approximant != nonGR1 && approximant != nonGR2 )
     */
    if( nonGRparams )
    {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* General sanity check the input parameters - only give warnings! */
    if( deltaF > 1. )
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaF = %e requested...This corresponds to a very short TD signal (with padding). Consider a smaller value.\n", __func__, deltaF);
    if( deltaF < 1./4096. )
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaF = %e requested...This corresponds to a very long TD signal. Consider a larger value.\n", __func__, deltaF);
    if( m1 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested...Perhaps you have a unit conversion error?\n", __func__, m1, m1/LAL_MSUN_SI);
    if( m2 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested...Perhaps you have a unit conversion error?\n", __func__, m2, m2/LAL_MSUN_SI);
    if( m1 + m2 > 1000. * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested...Signal not likely to be in band of ground-based detectors.\n", __func__, m1+m2, (m1+m2)/LAL_MSUN_SI);
    if( S1x*S1x + S1y*S1y + S1z*S1z > 1.000001 )
        XLALPrintWarning("XLAL Warning - %s: S1 = (%e,%e,%e) with norm > 1 requested...Are you sure you want to violate the Kerr bound?\n", __func__, S1x, S1y, S1z);
    if( S2x*S2x + S2y*S2y + S2z*S2z > 1.000001 )
        XLALPrintWarning("XLAL Warning - %s: S2 = (%e,%e,%e) with norm > 1 requested...Are you sure you want to violate the Kerr bound?\n", __func__, S2x, S2y, S2z);
    if( f_min < 1. )
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested...Check for errors, this could create a very long waveform.\n", __func__, f_min);
    if( f_min > 40.000001 )
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested...Check for errors, the signal will start in band.\n", __func__, f_min);

    switch (approximant)
    {
        /* inspiral-only models */
        case TaylorF2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2(hptilde, phiRef, deltaF, m1, m2,
                    S1z, S2z, f_min, f_max, r, lambda1, lambda2,
                    XLALSimInspiralGetSpinOrder(waveFlags),
                    XLALSimInspiralGetTidalOrder(waveFlags),
                    phaseO, amplitudeO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* The above returns h(f) for optimal orientation (i=0, Fp=1, Fc=0)
             * To get generic polarizations we multiply by incl. dependence
             * and note hc(f) \propto I * hp(f)
             */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            cfac = cos(i);
            pfac = 0.5 * (1. + cfac*cfac);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        /* non-spinning inspiral-merger-ringdown models */
        case IMRPhenomA:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkSpinsZero(S1x, S1y, S1z, S2x, S2y, S2z) )
                ABORT_NONZERO_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomAGenerateFD(hptilde, phiRef, deltaF, m1, m2,
                    f_min, f_max, r);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* The above returns h(f) for optimal orientation (i=0, Fp=1, Fc=0)
             * To get generic polarizations we multiply by incl. dependence
             * and note hc(f) \propto I * hp(f)
             */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            cfac = cos(i);
            pfac = 0.5 * (1. + cfac*cfac);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        /* spinning inspiral-only models */
        case SpinTaylorF2:
            /* Waveform-specific sanity checks */
            /* Sanity check unused fields of waveFlags */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE(waveFlags);
            if( S2z != 0. ) // This is a single-spin model
                ABORT_NONZERO_SPINS(waveFlags);
            LNhatx = sin(i);
            LNhaty = 0.;
            LNhatz = cos(i);
            /* Maximum PN amplitude order for precessing waveforms is 
             * MAX_PRECESSING_AMP_PN_ORDER */
            amplitudeO = 0; /* amplitudeO <= MAX_PRECESSING_AMP_PN_ORDER ? 
                    amplitudeO : MAX_PRECESSING_AMP_PN_ORDER */;
            /* Call the waveform driver routine */
            // FIXME: Note the HACK to use lambda1 as polarization angle psi!!
            ret = XLALSimInspiralSpinTaylorF2(hptilde, 0., phiRef, deltaF,
                    m1, m2, f_min, r, S1x, S1y, S1z,
                    LNhatx, LNhaty, LNhatz, phaseO, amplitudeO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            ret = XLALSimInspiralSpinTaylorF2(hctilde, LAL_PI/4.,phiRef, deltaF,
                    m1, m2, f_min, r, S1x, S1y, S1z,
                    LNhatx, LNhaty, LNhatz, phaseO, amplitudeO);
            break;

        /* FIXME: Comment out this case, as I don't have its source code */
        //case TaylorR2F4:
        //    /* Waveform-specific sanity checks */
        //    if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
        //        ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
        //    if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
        //        ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
        //    /* Call the waveform driver routine */
        //    ret = XLALSimInspiralTaylorR2F4(hptilde, phiRef, deltaF, m1, m2,
        //            S1z, S2z, f_min, r, phaseO, amplitudeO);
        //    break;

        case TaylorF2RedSpin:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2ReducedSpin(hptilde, phiRef, deltaF,
                    m1, m2, XLALSimInspiralTaylorF2ReducedSpinComputeChi(m1, m2, S1z, S2z),
                    f_min, f_max, r, phaseO, amplitudeO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* The above returns h(f) for optimal orientation (i=0, Fp=1, Fc=0)
             * To get generic polarizations we multiply by incl. dependence
             * and note hc(f) \propto I * hp(f)
             */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            cfac = cos(i);
            pfac = 0.5 * (1. + cfac*cfac);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        case TaylorF2RedSpinTidal:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimInspiralTaylorF2ReducedSpinTidal(hptilde,phiRef,deltaF,
                    m1, m2, XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    lambda1, lambda2, f_min, f_max, r, phaseO, amplitudeO);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* The above returns h(f) for optimal orientation (i=0, Fp=1, Fc=0)
             * To get generic polarizations we multiply by incl. dependence
             * and note hc(f) \propto I * hp(f)
             */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            cfac = cos(i);
            pfac = 0.5 * (1. + cfac*cfac);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        /* spinning inspiral-merger-ringdown models */
        case IMRPhenomB:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomBGenerateFD(hptilde, phiRef, deltaF, m1, m2,
                    XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, f_max, r);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* The above returns h(f) for optimal orientation (i=0, Fp=1, Fc=0)
             * To get generic polarizations we multiply by incl. dependence
             * and note hc(f) \propto I * hp(f)
             */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            cfac = cos(i);
            pfac = 0.5 * (1. + cfac*cfac);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        case IMRPhenomC:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralWaveformFlagsIsDefault(waveFlags) )
                ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags);
            if( !checkTransverseSpinsZero(S1x, S1y, S2x, S2y) )
                ABORT_NONZERO_TRANSVERSE_SPINS(waveFlags);
            if( !checkTidesZero(lambda1, lambda2) )
                ABORT_NONZERO_TIDES(waveFlags);
            /* Call the waveform driver routine */
            ret = XLALSimIMRPhenomCGenerateFD(hptilde, phiRef, deltaF, m1, m2,
                    XLALSimIMRPhenomBComputeChi(m1, m2, S1z, S2z),
                    f_min, f_max, r);
            if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);
            /* The above returns h(f) for optimal orientation (i=0, Fp=1, Fc=0)
             * To get generic polarizations we multiply by incl. dependence
             * and note hc(f) \propto I * hp(f)
             */
            *hctilde = XLALCreateCOMPLEX16FrequencySeries("FD hcross",
                    &((*hptilde)->epoch), (*hptilde)->f0, (*hptilde)->deltaF,
                    &((*hptilde)->sampleUnits), (*hptilde)->data->length);
            cfac = cos(i);
            pfac = 0.5 * (1. + cfac*cfac);
            for(j = 0; j < (*hptilde)->data->length; j++) {
                (*hctilde)->data->data[j] = I*cfac * (*hptilde)->data->data[j];
                (*hptilde)->data->data[j] *= pfac;
            }
            break;

        default:
            XLALPrintError("FD version of approximant not implemented in lalsimulation\n");
            XLAL_ERROR(XLAL_EINVAL);
    }

    if (ret == XLAL_FAILURE) XLAL_ERROR(XLAL_EFUNC);

    return ret;
}


/**
 * Interface to compute a set of -2 spin-weighted spherical harmonic modes
 * for a binary inspiral of any available amplitude and phase PN order.
 * The phasing is computed with any of the TaylorT1, T2, T3, T4 methods.
 *
 * It can also return the (2,2), (2,1), (3,3), (4,4), (5,5) modes of the EOBNRv2
 * model. Note that EOBNRv2 will ignore ampO, phaseO, lmax and f_ref arguments.
 */
SphHarmTimeSeries *XLALSimInspiralChooseTDModes(
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    int lmax,                                   /**< generate all modes with l <= lmax */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    REAL8 v0 = 1.;
    SphHarmTimeSeries *hlm = NULL;

    /* General sanity checks that will abort */
    /*
     * If non-GR approximants are added, change the below to
     * if( nonGRparams && approximant != nonGR1 && approximant != nonGR2 )
     */
    if( nonGRparams )
    {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* General sanity check the input parameters - only give warnings! */
    if( deltaT > 1. )
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaT = %e requested.\nPerhaps sample rate and time step size were swapped?\n", __func__, deltaT);
    if( deltaT < 1./16385. )
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaT = %e requested.\nCheck for errors, this could create very large time series.\n", __func__, deltaT);
    if( m1 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m1, m1/LAL_MSUN_SI);
    if( m2 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m2, m2/LAL_MSUN_SI);
    if( m1 + m2 > 1000. * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested.\nSignal not likely to be in band of ground-based detectors.\n", __func__, m1+m2, (m1+m2)/LAL_MSUN_SI);
    if( f_min < 1. )
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested.\nCheck for errors, this could create a very long waveform.\n", __func__, f_min);
    if( f_min > 40.000001 )
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested.\nCheck for errors, the signal will start in band.\n", __func__, f_min);

    switch (approximant)
    {
        case TaylorT1:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT1PNModes(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, lmax);
            break;
        case TaylorT2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT2PNModes(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, lmax);
            break;
        case TaylorT3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT3PNModes(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, lmax);
            break;
        case TaylorT4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT4PNModes(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, lmax);
            break;
        case EOBNRv2:
        case EOBNRv2HM:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimIMREOBNRv2Modes(phiRef, deltaT, m1, m2, f_min, r);
            // EOB driver only outputs modes with m>0, add m<0 modes by symmetry
            size_t l, j;
            int m;
            for( l=2; l <= XLALSphHarmTimeSeriesGetMaxL( hlm ); l++ ) {
                for( m=-l; m<0; m++){
                    COMPLEX16TimeSeries* inmode = XLALSphHarmTimeSeriesGetMode(
                            hlm, l, -m );
                    if( !inmode ) continue;
                    COMPLEX16TimeSeries* tmpmode = XLALCutCOMPLEX16TimeSeries(
                            inmode, 0, inmode->data->length );
                    for( j=0; j < tmpmode->data->length; j++ ){
                        tmpmode->data->data[j] = cpow(-1, l)
                            * conj( tmpmode->data->data[j] );
                    }
                    hlm = XLALSphHarmTimeSeriesAddMode( hlm, tmpmode, l, m );
                    XLALDestroyCOMPLEX16TimeSeries( tmpmode );
                }
            }
            break;

        default:
            XLALPrintError("Cannot generate modes for this approximant\n");
            XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);

    return hlm;
}

/**
 * Interface to compute a single -2 spin-weighted spherical harmonic mode
 * for a binary inspiral of any available amplitude and phase PN order.
 * The phasing is computed with any of the TaylorT1, T2, T3, T4 methods.
 */
COMPLEX16TimeSeries *XLALSimInspiralChooseTDMode(
    REAL8 phiRef,                               /**< reference orbital phase (rad) */
    REAL8 deltaT,                               /**< sampling interval (s) */
    REAL8 m1,                                   /**< mass of companion 1 (kg) */
    REAL8 m2,                                   /**< mass of companion 2 (kg) */
    REAL8 f_min,                                /**< starting GW frequency (Hz) */
    REAL8 f_ref,                                /**< reference GW frequency (Hz) */
    REAL8 r,                                    /**< distance of source (m) */
    REAL8 lambda1,                              /**< (tidal deformability of mass 1) / m1^5 (dimensionless) */
    REAL8 lambda2,                              /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
    LALSimInspiralWaveformFlags *waveFlags,     /**< Set of flags to control special behavior of some waveform families. Pass in NULL (or None in python) for default flags */
    LALSimInspiralTestGRParam *nonGRparams, 	/**< Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */
    int amplitudeO,                             /**< twice post-Newtonian amplitude order */
    int phaseO,                                 /**< twice post-Newtonian order */
    int l,                                      /**< l index of mode */
    int m,                                      /**< m index of mode */
    Approximant approximant                     /**< post-Newtonian approximant to use for waveform production */
    )
{
    REAL8 v0 = 1.;
    COMPLEX16TimeSeries *hlm;
    SphHarmTimeSeries *ts;

    /* General sanity checks that will abort */
    /*
     * If non-GR approximants are added, change the below to
     * if( nonGRparams && approximant != nonGR1 && approximant != nonGR2 )
     */
    if( nonGRparams )
    {
        XLALPrintError("XLAL Error - %s: Passed in non-NULL pointer to LALSimInspiralTestGRParam for an approximant that does not use LALSimInspiralTestGRParam\n", __func__);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }

    /* General sanity check the input parameters - only give warnings! */
    if( deltaT > 1. )
        XLALPrintWarning("XLAL Warning - %s: Large value of deltaT = %e requested.\nPerhaps sample rate and time step size were swapped?\n", __func__, deltaT);
    if( deltaT < 1./16385. )
        XLALPrintWarning("XLAL Warning - %s: Small value of deltaT = %e requested.\nCheck for errors, this could create very large time series.\n", __func__, deltaT);
    if( m1 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m1 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m1, m1/LAL_MSUN_SI);
    if( m2 < 0.09 * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Small value of m2 = %e (kg) = %e (Msun) requested.\nPerhaps you have a unit conversion error?\n", __func__, m2, m2/LAL_MSUN_SI);
    if( m1 + m2 > 1000. * LAL_MSUN_SI )
        XLALPrintWarning("XLAL Warning - %s: Large value of total mass m1+m2 = %e (kg) = %e (Msun) requested.\nSignal not likely to be in band of ground-based detectors.\n", __func__, m1+m2, (m1+m2)/LAL_MSUN_SI);
    if( f_min < 1. )
        XLALPrintWarning("XLAL Warning - %s: Small value of fmin = %e requested.\nCheck for errors, this could create a very long waveform.\n", __func__, f_min);
    if( f_min > 40.000001 )
        XLALPrintWarning("XLAL Warning - %s: Large value of fmin = %e requested.\nCheck for errors, the signal will start in band.\n", __func__, f_min);

    switch (approximant)
    {
        case TaylorT1:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT1PNMode(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case TaylorT2:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT2PNMode(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case TaylorT3:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault( XLALSimInspiralGetFrameAxis(waveFlags)) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault( XLALSimInspiralGetModesChoice(waveFlags)) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT3PNMode(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case TaylorT4:
            /* Waveform-specific sanity checks */
            if( !XLALSimInspiralFrameAxisIsDefault(
                    XLALSimInspiralGetFrameAxis(waveFlags) ) )
                ABORT_NONDEFAULT_FRAME_AXIS_NULL(waveFlags);
            if( !XLALSimInspiralModesChoiceIsDefault(
                    XLALSimInspiralGetModesChoice(waveFlags) ) )
                ABORT_NONDEFAULT_MODES_CHOICE_NULL(waveFlags);
            /* Call the waveform driver routine */
            hlm = XLALSimInspiralTaylorT4PNMode(phiRef, v0,
                    deltaT, m1, m2, f_min, f_ref, r, lambda1, lambda2,
                    XLALSimInspiralGetTidalOrder(waveFlags), amplitudeO,
                    phaseO, l, m);
            break;
        case EOBNRv2:
        case EOBNRv2HM:
            ts = XLALSimIMREOBNRv2Modes(phiRef, deltaT, m1, m2, f_min, r);
            hlm = XLALSphHarmTimeSeriesGetMode(ts, l, m);
            break;

        default:
            XLALPrintError("Cannot generate modes for this approximant\n");
            XLAL_ERROR_NULL(XLAL_EINVAL);
    }
    if ( !hlm )
        XLAL_ERROR_NULL(XLAL_EFUNC);

    return hlm;
}


/**
 * Checks whether the given approximant is implemented in lalsimulation's XLALSimInspiralChooseTDWaveform().
 *
 * returns 1 if the approximant is implemented, 0 otherwise.
 */
int XLALSimInspiralImplementedTDApproximants(
    Approximant approximant /**< post-Newtonian approximant for use in waveform production */
    )
{
    switch (approximant)
    {
        case TaylorEt:
        case TaylorT1:
        case TaylorT2:
        case TaylorT3:
        case TaylorT4:
        case EOBNRv2:
        case IMRPhenomA:
        case EOBNRv2HM:
        case SpinTaylorT2:
        case SpinTaylorT4:
        case IMRPhenomB:
        case PhenSpinTaylor:
        case PhenSpinTaylorRD:
        case SEOBNRv1:
            return 1;

        default:
            return 0;
    }
}

/**
 * Checks whether the given approximant is implemented in lalsimulation's XLALSimInspiralChooseFDWaveform().
 *
 * returns 1 if the approximant is implemented, 0 otherwise.
 */
int XLALSimInspiralImplementedFDApproximants(
    Approximant approximant /**< post-Newtonian approximant for use in waveform production */
    )
{
    switch (approximant)
    {
        case IMRPhenomA:
        case IMRPhenomB:
        case IMRPhenomC:
        //case TaylorR2F4:
        case TaylorF2:
        case SpinTaylorF2:
        case TaylorF2RedSpin:
        case TaylorF2RedSpinTidal:
            return 1;

        default:
            return 0;
    }
}

/**
 * XLAL function to determine approximant from a string.  The string need not
 * match exactly, only contain a member of the Approximant enum.
 */
int XLALGetApproximantFromString(const CHAR *inString)
{
#ifndef LAL_NDEBUG
  if ( !inString )
    XLAL_ERROR( XLAL_EFAULT );
#endif

  if ( strstr(inString, "TaylorF2RedSpinTidal" ) )
  {
    return TaylorF2RedSpinTidal;
  }
  else if ( strstr(inString, "TaylorF2RedSpin" ) )
  {
    return TaylorF2RedSpin;
  }
  else if ( strstr(inString, "SpinTaylorF2" ) )
  {
    return SpinTaylorF2;
  }
  else if ( strstr(inString, "TaylorF2" ) )
  {
    return TaylorF2;
  }
  else if ( strstr(inString, "TaylorR2F4" ) )
  {
    return TaylorR2F4;
  }
  else if ( strstr(inString, "PhenSpinTaylorRD" ) )
  {
    return PhenSpinTaylorRD;
  }
  else if ( strstr(inString, "PhenSpinTaylor" ) )
  {
    return PhenSpinTaylor;
  }
  else if ( strstr(inString, "SpinTaylorT2" ) )
  {
    return SpinTaylorT2;
  }
  else if ( strstr(inString, "SpinTaylorT4" ) )
  {
    return SpinTaylorT4;
  }
  else if ( strstr(inString, "SpinTaylorFrameless" ) )
  {
    return SpinTaylorFrameless;
  }
  else if ( strstr(inString, "SpinTaylorT3" ) )
  {
    return SpinTaylorT3;
  }
  else if ( strstr(inString, "SpinTaylor" ) )
  {
    return SpinTaylor;
  }
  else if ( strstr(inString, "SpinQuadTaylor" ) )
  {
    return SpinQuadTaylor;
  }
  else if ( strstr(inString, "TaylorT1" ) )
  {
    return TaylorT1;
  }
  else if ( strstr(inString, "TaylorT2" ) )
  {
    return TaylorT2;
  }
  else if ( strstr(inString, "TaylorT3" ) )
  {
    return TaylorT3;
  }
  else if ( strstr(inString, "TaylorT4" ) )
  {
    return TaylorT4;
  }
  else if ( strstr(inString, "IMRPhenomA" ) )
  {
    return IMRPhenomA;
  }
  else if ( strstr(inString, "IMRPhenomB" ) )
  {
    return IMRPhenomB;
  }
  else if ( strstr(inString, "IMRPhenomC" ) )
  {
    return IMRPhenomC;
  }
  else if ( strstr(inString, "IMRPhenomFA" ) )
  {
    return IMRPhenomFA;
  }
  else if ( strstr(inString, "IMRPhenomFB" ) )
  {
    return IMRPhenomFB;
  }
  else if ( strstr(inString, "SEOBNRv1" ) )
  {
    return SEOBNRv1;
  }
  else if ( strstr(inString, "EOBNRv2HM" ) )
  {
    return EOBNRv2HM;
  }
  else if ( strstr(inString, "EOBNRv2" ) )
  {
    return EOBNRv2;
  }
  else if ( strstr(inString, "EOBNR" ) )
  {
    return EOBNR;
  }
  else if ( strstr(inString, "EOB" ) )
  {
    return EOB;
  }
  else if ( strstr(inString, "AmpCorPPN" ) )
  {
    return AmpCorPPN;
  }
  else if ( strstr(inString, "GeneratePPN" ) )
  {
    return GeneratePPN;
  }
  else if ( strstr(inString, "NumRelNinja2" ) )
  {
    return NumRelNinja2;
  }
  else if ( strstr(inString, "NumRel" ) )
  {
    return NumRel;
  }
  else if ( strstr(inString, "Ninja2" ) )
  {
    return NumRelNinja2;
  }
  else if ( strstr(inString, "FindChirpSP" ) )
  {
    return FindChirpSP;
  }
  else if ( strstr(inString, "FindChirpPTF" ) )
  {
    return FindChirpPTF;
  }
  else if ( strstr(inString, "TaylorEt" ) )
  {
    return TaylorEt;
  }
  else if ( strstr(inString, "TaylorN" ) )
  {
    return TaylorN;
  }
  else if ( strstr(inString, "TaylorF1" ) )
  {
    return TaylorF1;
  }
  else if ( strstr(inString, "PadeT1" ) )
  {
    return PadeT1;
  }
  else if ( strstr(inString, "PadeF1" ) )
  {
    return PadeF1;
  }
  else if ( strstr(inString, "BCVSpin" ) )
  {
    return BCVSpin;
  }
  else if ( strstr(inString, "BCVC" ) )
  {
    return BCVC;
  }
  else if ( strstr(inString, "BCV" ) )
  {
    return BCV;
  }
  else if ( strstr(inString, "FrameFile" ) )
  {
    return FrameFile;
  }
  else if ( strstr(inString, "Eccentricity" ) )
  {
    return Eccentricity;
  }
  else
  {
    XLALPrintError( "Cannot parse approximant from string: %s \n", inString );
    XLAL_ERROR( XLAL_EINVAL );
  }
}

/**
 * XLAL function to determine string from approximant enum.
 * This function needs to be updated when new approximants are added.
 */
char* XLALGetStringFromApproximant(Approximant approximant)
{
  switch (approximant)
  {
    case TaylorF2RedSpinTidal:
      return strdup("TaylorF2RedSpinTidal");
    case TaylorF2RedSpin:
      return strdup("TaylorF2RedSpin");
    case TaylorF2:
      return strdup("TaylorF2");
    case PhenSpinTaylor:
      return strdup("PhenSpinTaylor");
    case TaylorR2F4:
      return strdup("TaylorR2F4");
    case PhenSpinTaylorRD:
      return strdup("PhenSpinTaylorRD");
    case SpinTaylorF2:
      return strdup("SpinTaylorF2");
    case SpinTaylorT2:
      return strdup("SpinTaylorT2");
    case SpinTaylorT4:
      return strdup("SpinTaylorT4");
    case SpinTaylorFrameless:
      return strdup("SpinTaylorFrameless");
    case SpinTaylorT3:
      return strdup("SpinTaylorT3");
    case SpinTaylor:
      return strdup("SpinTaylor");
    case SpinQuadTaylor:
      return strdup("SpinQuadTaylor");
    case TaylorT1:
      return strdup("TaylorT1");
    case TaylorT2:
      return strdup("TaylorT2");
    case TaylorT3:
      return strdup("TaylorT3");
    case TaylorT4:
      return strdup("TaylorT4");
    case IMRPhenomA:
      return strdup("IMRPhenomA");
    case IMRPhenomB:
      return strdup("IMRPhenomB");
    case IMRPhenomC:
      return strdup("IMRPhenomC");
    case IMRPhenomFA:
      return strdup("IMRPhenomFA");
    case IMRPhenomFB:
      return strdup("IMRPhenomFB");
    case SEOBNRv1:
      return strdup("SEOBNRv1");
    case EOBNRv2HM:
      return strdup("EOBNRv2HM");
    case EOBNRv2:
      return strdup("EOBNRv2");
    case EOBNR:
      return strdup("EOBNR");
    case EOB:
      return strdup("EOB");
    case AmpCorPPN:
      return strdup("AmpCorPPN");
    case GeneratePPN:
      return strdup("GeneratePPN");
    case NumRelNinja2:
      return strdup("NumRelNinja2");
    case NumRel:
      return strdup("NumRel");
    case FindChirpSP:
      return strdup("FindChirpSP");
    case FindChirpPTF:  
      return strdup("FindChirpPTF");
    case TaylorEt:
      return strdup("TaylorET");
    case TaylorN:  
      return strdup("TaylorN");
    case TaylorF1:
      return strdup("TaylorF1");
    case PadeT1:
      return strdup("PadeT1");
    case PadeF1:
      return strdup("PadeF1");
    case BCVSpin:
      return strdup("BCVSpin");
    case BCVC:
      return strdup("BCVC");
    case BCV:
      return strdup("BCV");
    case FrameFile:
      return strdup("FrameFile");
    case Eccentricity:
      return strdup("Eccentricity");
    default:
        XLALPrintError("Not a valid approximant\n");
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
}

/**
 * XLAL function to determine PN order from a string.  The string need not
 * match exactly, only contain a member of the LALPNOrder enum.
 */
int XLALGetOrderFromString(const CHAR *inString)
{

#ifndef LAL_NDEBUG
  if ( !inString )
    XLAL_ERROR( XLAL_EFAULT );
#endif

  if ( strstr(inString, "newtonian") )
  {
    return LAL_PNORDER_NEWTONIAN;
  }
  else if ( strstr(inString, "oneHalfPN") )
  {
    return LAL_PNORDER_HALF;
  }
  else if ( strstr(inString, "onePN") )
  {
    return LAL_PNORDER_ONE;
  }
  else if ( strstr(inString, "onePointFivePN") )
  {
    return LAL_PNORDER_ONE_POINT_FIVE;
  }
  else if ( strstr(inString, "twoPN") )
  {
    return LAL_PNORDER_TWO;
  }
  else if ( strstr(inString, "twoPointFivePN") )
  {
    return LAL_PNORDER_TWO_POINT_FIVE;
  }
  else if (strstr(inString, "threePN") )
  {
    return LAL_PNORDER_THREE;
  }
  else if ( strstr(inString, "threePointFivePN") )
  {
    return LAL_PNORDER_THREE_POINT_FIVE;
  }
  else if ( strstr(inString, "pseudoFourPN") )
  {
    return LAL_PNORDER_PSEUDO_FOUR;
  }
  else
  {
    XLALPrintError( "Cannot parse order from string: %s\n", inString );
    XLAL_ERROR( XLAL_EINVAL );
  }
}

/**
 * XLAL function to determine tapering flag from a string.  The string must
 * match exactly with a member of the LALSimInspiralApplyTaper enum.
 */
int XLALGetTaperFromString(const CHAR *inString)
{
  if ( ! strcmp( "TAPER_NONE", inString ) )
  {
    return LAL_SIM_INSPIRAL_TAPER_NONE;
  }
  else if ( ! strcmp( "TAPER_START", inString ) )
  {
    return LAL_SIM_INSPIRAL_TAPER_START;
  }
  else if ( ! strcmp( "TAPER_END", inString ) )
  {
    return LAL_SIM_INSPIRAL_TAPER_END;
  }
  else if ( ! strcmp( "TAPER_STARTEND", inString ) )
  {
    return LAL_SIM_INSPIRAL_TAPER_STARTEND;
  }
  else
  {
    XLALPrintError( "Invalid injection tapering option specified: %s\n", inString );
    XLAL_ERROR( XLAL_EINVAL );
  }
}

/**
 * XLAL function to determine LALSimInspiralInteraction from a string.
 *
 * TODO: return the bit sum if the string is a concatenation of several
 * interaction terms. Also make names match cases of enum.
 */
int XLALGetInteractionFromString(const CHAR *inString) 
{
  if (strstr(inString, "NO")) {
    return LAL_SIM_INSPIRAL_INTERACTION_NONE;
  } else if (strstr(inString, "SO15")) {
    return LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN;
  } else if (strstr(inString,"SS")) {
    return LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN;
  } else if (strstr(inString,"SELF")) {
    return LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN;
  } else if (strstr(inString, "QM")) {
    return LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN;
  } else if (strstr(inString, "SO25")) {
    return LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN;
  } else if (strstr(inString, "SO")) {
    return LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_3PN;
  } else if (strstr(inString, "ALL_SPIN")) {
    return LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN;
  } else if (strstr(inString, "TIDAL5PN")) {
    return LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN;
  } else if (strstr(inString, "TIDAL")) {
    return LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN;
  } else if (strstr(inString, "ALL")){
    return LAL_SIM_INSPIRAL_INTERACTION_ALL;
  } else {
    XLALPrintError( "Cannot parse LALSimInspiralInteraction from string: %s\n Please add 'ALL' to the above string for including all spin interactions\n", inString );
    XLAL_ERROR( XLAL_EINVAL );
  }
}

/**
 * XLAL function to determine axis choice flag from a string.
 * The string need not match exactly, only contain a member
 * of the LALSimInspiralFrameAxis enum.  Will return default case
 * 'View' (line of sight) if the string contains no match.
 */
int XLALGetFrameAxisFromString(const CHAR *inString) 
{
  if (strstr(inString, "TotalJ"))
    return LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J;
  else if (strstr(inString, "OrbitalL"))
    return LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L;
  else
    return LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW;
}

/**
 * XLAL function to determine modes choice from a string.
 */
int XLALGetHigherModesFromString(const CHAR *inString)
{
  if (strstr(inString, "L2"))
    return LAL_SIM_INSPIRAL_MODES_CHOICE_RESTRICTED;
  else if (strstr(inString, "L3"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_3L;
  else if (strstr(inString, "L4"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_3L;
  else if (strstr(inString, "L23"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_2AND3L;
  else if (strstr(inString, "L24"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_2AND4L;
  else if (strstr(inString, "L34"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_3AND4L;
  else if (strstr(inString, "L234"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_2AND3AND4L;
  else if (strstr(inString, "L5"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_5L;
  else if (strstr(inString, "L25"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_2AND5L;
  else if (strstr(inString, "L35"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_3AND5L;
  else if (strstr(inString, "L45"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_4AND5L;
  else if (strstr(inString, "L235"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_2AND3AND5L;
  else if (strstr(inString, "L245"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_2AND4AND5L;
  else if (strstr(inString, "L345"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_3AND4AND5L;
  else if (strstr(inString, "ALL"))
    return  LAL_SIM_INSPIRAL_MODES_CHOICE_ALL;
  else {
    XLALPrintError(" Error: invalid value %s for mode choice\n",inString);
    return 0;
  }
}

int XLALSimInspiralGetSpinSupportFromApproximant(Approximant approx){

  SpinSupport spin_support=LAL_SIM_INSPIRAL_NUMSPINSUPPORT;
  switch (approx)
  {
    case SpinTaylor:
    case SpinTaylorFrameless:
    case SpinTaylorT4:
    case SpinTaylorT2:
    case PhenSpinTaylor:
    case PhenSpinTaylorRD:
    case SpinTaylorT3:
      spin_support=LAL_SIM_INSPIRAL_PRECESSINGSPIN;
      break;
    case SpinTaylorF2:
    case FindChirpPTF:
      spin_support=LAL_SIM_INSPIRAL_SINGLESPIN;
      break;
    case TaylorF2:
    case TaylorF2RedSpin:
    case TaylorF2RedSpinTidal:
    case IMRPhenomB:
    case IMRPhenomC:
    case SEOBNRv1:
    case TaylorR2F4:
    case IMRPhenomFB:
    case FindChirpSP:
      spin_support=LAL_SIM_INSPIRAL_ALIGNEDSPIN;
      break;
    case TaylorEt:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case IMRPhenomA:
    case EOBNRv2HM:
    case EOBNRv2:
    case EOBNR:
    case EOB:
    case IMRPhenomFA:
    case GeneratePPN:
      spin_support=LAL_SIM_INSPIRAL_SPINLESS;
      break;
    default:
      XLALPrintError("Approximant not supported by lalsimuation TD/FD routines \n");
      XLAL_ERROR(XLAL_EINVAL);
    }

    return spin_support;

}
