/*
 * Copyright (C) 2008 J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria, D. Keppel
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
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include "LALSimInspiraldEnergyFlux.c"
#include "LALSimInspiralPNCoefficients.c"
#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* v at isco */
#define LALSIMINSPIRAL_T1_VISCO 1.L/sqrt(6.L)
/* use error codes above 1024 to avoid conflicts with GSL */
#define LALSIMINSPIRAL_T1_TEST_ISCO 1025
/* Number of variables used for these waveforms */
#define LALSIMINSPIRAL_NUM_T1_VARIABLES 2
/* absolute and relative tolerance for adaptive Runge-Kutta ODE integrator */
#define LALSIMINSPIRAL_T1_ABSOLUTE_TOLERANCE 1.e-12
#define LALSIMINSPIRAL_T1_RELATIVE_TOLERANCE 1.e-12

/**
 * This structure contains the intrinsic parameters and post-newtonian 
 * co-efficients for the denergy/dv and flux expansions. 
 * These are computed by XLALSimInspiralTaylorT1Setup routine.
 */

typedef struct
{
	/* Angular velocity coefficient */
	REAL8 av;

	/* Taylor expansion coefficents in domega/dt */
	expnCoeffsdEnergyFlux akdEF;

	/* symmetric mass ratio and total mass */
	REAL8 nu,m,mchirp,chi1,chi2,lambda1,lambda2;
} expnCoeffsTaylorT1;

typedef REAL8 (SimInspiralTaylorT1Energy)(
	REAL8 v, /**< post-Newtonian parameter */
	expnCoeffsdEnergyFlux *ak
);

typedef REAL8 (SimInspiralTaylorT1dEnergy)(
	REAL8 v, /**< post-Newtonian parameter */
	expnCoeffsdEnergyFlux *ak
);

typedef REAL8 (SimInspiralTaylorT1Flux)(
	REAL8 v, /**< post-Newtonian parameter */
	expnCoeffsdEnergyFlux *ak
);

/**
 * This strucuture contains pointers to the functions for calculating
 * the post-newtonian terms at the desired order. They can be set by
 * XLALSimInspiralTaylorT1Setup by passing an appropriate PN order. 
 */

typedef struct
tagexpnFuncTaylorT1
{
	SimInspiralTaylorT1Energy *energy;
	SimInspiralTaylorT1dEnergy *dEnergy;
	SimInspiralTaylorT1Flux *flux;
} expnFuncTaylorT1;

typedef struct
{
	REAL8 (*dEdv)(REAL8 v, expnCoeffsdEnergyFlux *ak);
	REAL8 (*flux)(REAL8 v, expnCoeffsdEnergyFlux *ak);
	expnCoeffsTaylorT1 ak;
}XLALSimInspiralTaylorT1PNEvolveOrbitParams;

/** 
 * This function is used in the call to the integrator.
 */
static int 
XLALSimInspiralTaylorT1PNEvolveOrbitIntegrand(double UNUSED t, const double y[], double ydot[], void *params)
{
	XLALSimInspiralTaylorT1PNEvolveOrbitParams* p = (XLALSimInspiralTaylorT1PNEvolveOrbitParams*)params;
	ydot[0] = -p->ak.av*p->flux(y[0],&p->ak.akdEF)/p->dEdv(y[0],&p->ak.akdEF);
	ydot[1] = y[0]*y[0]*y[0]*p->ak.av;
	return GSL_SUCCESS;
}


/**
 * This function is used in the call to the integrator to determine the stopping condition.
 */
static int
XLALSimInspiralTaylorT1StoppingTest(double UNUSED t, const double y[], double UNUSED ydot[], void UNUSED *params)
{
	if (y[0] >= LALSIMINSPIRAL_T1_VISCO) /* frequency above ISCO */
		return LALSIMINSPIRAL_T1_TEST_ISCO;
	else /* Step successful, continue integrating */
		return GSL_SUCCESS;
}


/**
 * Set up the expnCoeffsTaylorT1 and expnFuncTaylorT1 structures for
 * generating a TaylorT1 waveform and select the post-newtonian
 * functions corresponding to the desired order. 
 *
 * Inputs given in SI units.
 */
static int 
XLALSimInspiralTaylorT1Setup(
    expnCoeffsTaylorT1 *ak,			/**< coefficients for TaylorT1 evolution [modified] */
    expnFuncTaylorT1 *f,			/**< functions for TaylorT1 evolution [modified] */
    REAL8 m1,					/**< mass of companion 1 (kg) */
    REAL8 m2,					/**< mass of companion 2 (kg) */
    REAL8 lambda1,				/**< (tidal deformability of body 1)/(mass of body 1)^5 */
    REAL8 lambda2,				/**< (tidal deformability of body 2)/(mass of body 2)^5 */
    LALSimInspiralInteraction interactionFlags,	/**< flag to control spin and tidal effects */
    int O					/**< twice post-Newtonian order */
)
{
    ak->m = m1 + m2;
    REAL8 mu = m1 * m2 / ak->m;
    ak->nu = mu/ak->m;
    ak->chi1 = m1/ak->m;
    ak->chi2 = m2/ak->m;
    ak->mchirp = ak->m * pow(ak->nu, 0.6);
    /* convert mchirp from kg to s */
    ak->mchirp *= LAL_G_SI / pow(LAL_C_SI, 3.0);

    /* Angular velocity co-efficient */
    ak->av = pow(LAL_C_SI, 3.0)/(LAL_G_SI*ak->m);

    /* Taylor co-efficients for E(v). */
    ak->akdEF.ETaN = XLALSimInspiralPNEnergy_0PNCoeff(ak->nu);
    ak->akdEF.ETa1 = XLALSimInspiralPNEnergy_2PNCoeff(ak->nu);
    ak->akdEF.ETa2 = XLALSimInspiralPNEnergy_4PNCoeff(ak->nu);
    ak->akdEF.ETa3 = XLALSimInspiralPNEnergy_6PNCoeff(ak->nu);

    /* Taylor co-efficients for dE(v)/dv. */
    ak->akdEF.dETaN = 2.0 * ak->akdEF.ETaN;
    ak->akdEF.dETa1 = 2.0 * ak->akdEF.ETa1;
    ak->akdEF.dETa2 = 3.0 * ak->akdEF.ETa2;
    ak->akdEF.dETa3 = 4.0 * ak->akdEF.ETa3;

    /* Taylor co-efficients for flux. */
    ak->akdEF.FTaN = XLALSimInspiralTaylorT1Flux_0PNCoeff(ak->nu);
    ak->akdEF.FTa2 = XLALSimInspiralTaylorT1Flux_2PNCoeff(ak->nu);
    ak->akdEF.FTa3 = XLALSimInspiralTaylorT1Flux_3PNCoeff(ak->nu);
    ak->akdEF.FTa4 = XLALSimInspiralTaylorT1Flux_4PNCoeff(ak->nu);
    ak->akdEF.FTa5 = XLALSimInspiralTaylorT1Flux_5PNCoeff(ak->nu);
    ak->akdEF.FTa6 = XLALSimInspiralTaylorT1Flux_6PNCoeff(ak->nu);
    ak->akdEF.FTl6 = XLALSimInspiralTaylorT1Flux_6PNLogCoeff(ak->nu);
    ak->akdEF.FTa7 = XLALSimInspiralTaylorT1Flux_7PNCoeff(ak->nu);

    /* Tidal co-efficients for E(v), dE/dv, and flux */
    ak->akdEF.ETa5  = 0.;
    ak->akdEF.ETa6  = 0.;
    ak->akdEF.FTa10 = 0.;
    ak->akdEF.FTa12 = 0.;
    if( interactionFlags >= LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN)
    {
        ak->akdEF.ETa5  = XLALSimInspiralPNEnergy_10PNTidalCoeff(ak->chi2,ak->chi1,lambda1)
                        + XLALSimInspiralPNEnergy_10PNTidalCoeff(ak->chi1,ak->chi2,lambda2);
        ak->akdEF.FTa10 = XLALSimInspiralTaylorT1Flux_10PNTidalCoeff(ak->chi1,lambda1)
                        + XLALSimInspiralTaylorT1Flux_10PNTidalCoeff(ak->chi2,lambda2);
    }
    if( interactionFlags >= LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN )
    {
        ak->akdEF.ETa6  = XLALSimInspiralPNEnergy_12PNTidalCoeff(ak->chi2,ak->chi1,lambda1)
                        + XLALSimInspiralPNEnergy_12PNTidalCoeff(ak->chi1,ak->chi2,lambda2);
        ak->akdEF.FTa12 = XLALSimInspiralTaylorT1Flux_12PNTidalCoeff(ak->chi1,lambda1)
                        + XLALSimInspiralTaylorT1Flux_12PNTidalCoeff(ak->chi2,lambda2);
    }
    ak->akdEF.dETa5 = 6.0 * ak->akdEF.ETa5;
    ak->akdEF.dETa6 = 7.0 * ak->akdEF.ETa6;

    switch (O)
    {
        case 0:
            f->energy = &XLALSimInspiralEt0;
            f->dEnergy = &XLALSimInspiraldEt0;
            f->flux = &XLALSimInspiralFt0;
            break;
        case 1:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for PN order %d\n", __func__,O);
            XLAL_ERROR(XLAL_EINVAL);
            break;
        case 2:
            f->energy = &XLALSimInspiralEt2;
            f->dEnergy = &XLALSimInspiraldEt2;
            f->flux = &XLALSimInspiralFt2;
            break;
        case 3:
            f->energy = &XLALSimInspiralEt2;
            f->dEnergy = &XLALSimInspiraldEt2;
            f->flux = &XLALSimInspiralFt3;
            break;
        case 4:
            f->energy = &XLALSimInspiralEt4;
            f->dEnergy = &XLALSimInspiraldEt4;
            f->flux = &XLALSimInspiralFt4;
            break;
        case 5:
            f->energy = &XLALSimInspiralEt4;
            f->dEnergy = &XLALSimInspiraldEt4;
            f->flux = &XLALSimInspiralFt5;
            break;
        case 6:
            f->energy = &XLALSimInspiralEt6;
            f->dEnergy = &XLALSimInspiraldEt6;
            f->flux = &XLALSimInspiralFt6;
            break;
        case 7:
        case -1:
            f->energy = &XLALSimInspiralEt6;
            f->dEnergy = &XLALSimInspiraldEt6;
            f->flux = &XLALSimInspiralFt7;
            break;
        case 8:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for PN order %d\n", __func__,O);
            XLAL_ERROR(XLAL_EINVAL);
            break;
        default:
            XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
    }

  return 0;
}


/**
 * Evolves a post-Newtonian orbit using the Taylor T1 method.
 *
 * See Section IIIA: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */
int XLALSimInspiralTaylorT1PNEvolveOrbit(
		REAL8TimeSeries **V,                        /**< post-Newtonian parameter [returned] */
		REAL8TimeSeries **phi,                      /**< orbital phase [returned] */
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int O                                       /**< twice post-Newtonian order */
		)
{
	double lengths, VRef = 0.;
	int len, intreturn, idx, idxRef = 0;
	XLALSimInspiralTaylorT1PNEvolveOrbitParams params;
	ark4GSLIntegrator *integrator = NULL;
	expnFuncTaylorT1 expnfunc;
	expnCoeffsTaylorT1 ak;

	if(XLALSimInspiralTaylorT1Setup(&ak,&expnfunc,m1,m2,lambda1,lambda2,interactionFlags,O))
		XLAL_ERROR(XLAL_EFUNC);

	params.flux=expnfunc.flux;
	params.dEdv=expnfunc.dEnergy;
	params.ak=ak;

	LIGOTimeGPS tc = LIGOTIMEGPSZERO;
	double yinit[LALSIMINSPIRAL_NUM_T1_VARIABLES];
	REAL8Array *yout;

	/* length estimation (Newtonian) */
	/* since integration is adaptive, we could use a better estimate */
	lengths = (5.0/256.0) * pow(LAL_PI, -8.0/3.0) * pow(f_min * ak.mchirp, -5.0/3.0) / f_min;

	yinit[0] = cbrt(LAL_PI * LAL_G_SI * ak.m * f_min) / LAL_C_SI;
	yinit[1] = 0.;

	/* initialize the integrator */
	integrator = XLALAdaptiveRungeKutta4Init(LALSIMINSPIRAL_NUM_T1_VARIABLES,
		XLALSimInspiralTaylorT1PNEvolveOrbitIntegrand,
		XLALSimInspiralTaylorT1StoppingTest,
		LALSIMINSPIRAL_T1_ABSOLUTE_TOLERANCE, LALSIMINSPIRAL_T1_RELATIVE_TOLERANCE);
	if( !integrator )
	{
		XLALPrintError("XLAL Error - %s: Cannot allocate integrator\n", __func__);
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* stop the integration only when the test is true */
	integrator->stopontestonly = 1;

	/* run the integration */
	len = XLALAdaptiveRungeKutta4Hermite(integrator, (void *) &params, yinit, 0.0, lengths, deltaT, &yout);

	intreturn = integrator->returncode;
	XLALAdaptiveRungeKutta4Free(integrator);

	if (!len) 
	{
		XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intreturn);
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* Adjust tStart so last sample is at time=0 */
	XLALGPSAdd(&tc, -1.0*(len-1)*deltaT);

	/* allocate memory for output vectors */
	*V = XLALCreateREAL8TimeSeries( "ORBITAL_VELOCITY_PARAMETER", &tc, 0., deltaT, &lalDimensionlessUnit, len);
	*phi = XLALCreateREAL8TimeSeries( "ORBITAL_PHASE", &tc, 0., deltaT, &lalDimensionlessUnit, len);

	if ( !V || !phi )
	{
		XLALDestroyREAL8Array(yout);
		XLAL_ERROR(XLAL_EFUNC);
	}

	/* Do a constant phase shift to get desired value of phiRef */
	/* For fRef==0, phiRef is phase of last sample */
	if( fRef == 0. )
		phiRef -= yout->data[3*len-1];
	/* For fRef==fmin, phiRef is phase of first sample */
	else if( fRef == f_min )
		phiRef -= yout->data[2*len];
	/* phiRef is phase when f==fRef */
	else
	{
		VRef = pow(LAL_PI * LAL_G_SI*(m1+m2) * fRef, 1./3.) / LAL_C_SI;
		idx = 0;
		do {
			idxRef = idx;
			idx++;
		} while (yout->data[len+idx] <= VRef);
		phiRef -= yout->data[2*len+idxRef];
	}

	/* Copy time series of dynamical variables */
	/* from yout array returned by integrator to output time series */
	/* Note the first 'len' members of yout are the time steps */
	for( idx = 0; idx < len; idx++ )
	{	
		(*V)->data->data[idx]   = yout->data[len+idx];
		(*phi)->data->data[idx]	= yout->data[2*len+idx] + phiRef;
	}

	XLALDestroyREAL8Array(yout);

	return (int)(*V)->data->length;
}


/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine allows the user to specify different pN orders
 * for phasing calcuation vs. amplitude calculations.
 */
int XLALSimInspiralTaylorT1PNGenerator(
		REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
		REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 v0,                                   /**< tail-term gauge choice (default = 1) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 r,                                    /**< distance of source (m) */
		REAL8 i,                                    /**< inclination of source (rad) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int amplitudeO,                             /**< twice post-Newtonian amplitude order */
		int phaseO                                  /**< twice post-Newtonian phase order */
		)
{
	/* The Schwarzschild ISCO frequency - for sanity checking fRef */
	REAL8 fISCO = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_G_SI);

	/* Sanity check fRef value */
	if( fRef < 0. )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be >= 0\n", 
				__func__, fRef);
		XLAL_ERROR(XLAL_EINVAL);
	}
	if( fRef != 0. && fRef < f_min )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be > fStart = %f\n", 
				__func__, fRef, f_min);
		XLAL_ERROR(XLAL_EINVAL);
	}
	if( fRef >= fISCO )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be < Schwar. ISCO=%f\n",
				__func__, fRef, fISCO);
		XLAL_ERROR(XLAL_EINVAL);
	}

	REAL8TimeSeries *V;
	REAL8TimeSeries *phi;
	int status;
	int n;
	n = XLALSimInspiralTaylorT1PNEvolveOrbit(&V, &phi, phiRef, deltaT,
			m1, m2, f_min, fRef, lambda1, lambda2, interactionFlags, phaseO);
	if ( n < 0 )
		XLAL_ERROR(XLAL_EFUNC);
	status = XLALSimInspiralPNPolarizationWaveforms(hplus, hcross, V, phi,
			v0, m1, m2, r, i, amplitudeO);
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(V);
	if ( status < 0 )
		XLAL_ERROR(XLAL_EFUNC);
	return n;
}

SphHarmTimeSeries *XLALSimInspiralTaylorT1PNModes(
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 v0,                                   /**< tail-term gauge choice (default = 1) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 r,                                    /**< distance of source (m) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(individual mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(individual mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int amplitudeO,                             /**< twice post-Newtonian amplitude order */
		int phaseO,                                 /**< twice post-Newtonian phase order */
		int l                                       /**< generate all modes with l <= lmax */
		)
{
	SphHarmTimeSeries *hlm = NULL;
	/* The Schwarzschild ISCO frequency - for sanity checking fRef */
	REAL8 fISCO = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_G_SI);

	/* Sanity check fRef value */
	if( fRef < 0. )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be >= 0\n", 
				__func__, fRef);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}
	if( fRef != 0. && fRef < f_min )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be > fStart = %f\n", 
				__func__, fRef, f_min);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}
	if( fRef >= fISCO )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be < Schwar. ISCO=%f\n",
				__func__, fRef, fISCO);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}

	REAL8TimeSeries *V;
	REAL8TimeSeries *phi;
	int n;
	n = XLALSimInspiralTaylorT1PNEvolveOrbit(&V, &phi, phiRef, deltaT,
			m1, m2, f_min, fRef, lambda1, lambda2,
			interactionFlags, phaseO);
	if ( n < 0 )
		XLAL_ERROR_NULL(XLAL_EFUNC);
	int mi, li;
	COMPLEX16TimeSeries *hxx;
	for(li=0; li<=l; li++){
		for(mi=-l; mi<=l; mi++){
			hxx = XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(V, phi,
				v0, m1, m2, r, amplitudeO, li, mi);
			if ( !hxx ){
				XLAL_ERROR_NULL(XLAL_EFUNC);
			}
	 		XLALSphHarmTimeSeriesAddMode(hlm, hxx, li, mi);
			XLALDestroyCOMPLEX16TimeSeries(hxx);
		}
	}
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(V);
	return hlm;
}

/**
 * Driver routine to compute the -2 spin-weighted spherical harmonic mode
 * using TaylorT1 phasing.
 */
COMPLEX16TimeSeries *XLALSimInspiralTaylorT1PNMode(
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 v0,                                   /**< tail-term gauge choice (default = 1) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 r,                                    /**< distance of source (m) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(individual mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(individual mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int amplitudeO,                             /**< twice post-Newtonian amplitude order */
		int phaseO,                                 /**< twice post-Newtonian phase order */
		int l,                                      /**< l index of mode */
		int m                                       /**< m index of mode */
		)
{
	COMPLEX16TimeSeries *hlm;
	/* The Schwarzschild ISCO frequency - for sanity checking fRef */
	REAL8 fISCO = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_G_SI);

	/* Sanity check fRef value */
	if( fRef < 0. )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be >= 0\n", 
				__func__, fRef);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}
	if( fRef != 0. && fRef < f_min )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be > fStart = %f\n", 
				__func__, fRef, f_min);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}
	if( fRef >= fISCO )
	{
		XLALPrintError("XLAL Error - %s: fRef = %f must be < Schwar. ISCO=%f\n",
				__func__, fRef, fISCO);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}

	REAL8TimeSeries *V;
	REAL8TimeSeries *phi;
	int n;
	n = XLALSimInspiralTaylorT1PNEvolveOrbit(&V, &phi, phiRef, deltaT,
			m1, m2, f_min, fRef, lambda1, lambda2,
			interactionFlags, phaseO);
	if ( n < 0 )
		XLAL_ERROR_NULL(XLAL_EFUNC);
	hlm = XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(V, phi,
			v0, m1, m2, r, amplitudeO, l, m);
	if ( !hlm )
		XLAL_ERROR_NULL(XLAL_EFUNC);
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(V);
	return hlm;
}


/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine uses the same pN order for phasing and amplitude
 * (unless the order is -1 in which case the highest available
 * order is used for both of these -- which might not be the same).
 *
 * Constant log term in amplitude set to 1.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT1PN(
		REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
		REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 r,                                    /**< distance of source (m) */
		REAL8 i,                                    /**< inclination of source (rad) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int O                                       /**< twice post-Newtonian order */
		)
{
	/* set v0 to default value 1 */
	return XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, phiRef, 1.0,
			deltaT, m1, m2, f_min, fRef, r, i, lambda1, lambda2,
			interactionFlags, O, O);
}


/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 *
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 *
 * Constant log term in amplitude set to 1.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT1PNRestricted(
		REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
		REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
		REAL8 phiRef,                               /**< reference orbital phase (rad) */
		REAL8 deltaT,                               /**< sampling interval (s) */
		REAL8 m1,                                   /**< mass of companion 1 (kg) */
		REAL8 m2,                                   /**< mass of companion 2 (kg) */
		REAL8 f_min,                                /**< starting GW frequency (Hz) */
		REAL8 fRef,                                 /**< reference GW frequency (Hz) */
		REAL8 r,                                    /**< distance of source (m)*/
		REAL8 i,                                    /**< inclination of source (rad) */
		REAL8 lambda1,                              /**< (tidal deformability of body 1)/(mass of body 1)^5 */
		REAL8 lambda2,                              /**< (tidal deformability of body 2)/(mass of body 2)^5 */
		LALSimInspiralInteraction interactionFlags, /**< flag to control spin and tidal effects */
		int O                                       /**< twice post-Newtonian phase order */
		)
{
	/* use Newtonian order for amplitude */
	/* set v0 to default value 1 */
	return XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, phiRef, 1.0,
			deltaT, m1, m2, f_min, fRef, r, i, lambda1, lambda2,
			interactionFlags, 0, O);
}


#if 0
#include <lal/PrintFTSeries.h>
#include <lal/PrintFTSeries.h>
extern int lalDebugLevel;
int main(void)
{
	LIGOTimeGPS tc = { 888888888, 222222222 };
	REAL8 phic = 1.0;
	REAL8 deltaT = 1.0/16384.0;
	REAL8 m1 = 1.4*LAL_MSUN_SI;
	REAL8 m2 = 1.4*LAL_MSUN_SI;
	REAL8 r = 1e6*LAL_PC_SI;
	REAL8 i = 0.5*LAL_PI;
	REAL8 f_min = 100.0;
	REAL8 fRef = 0.;
	int O = -1;
	REAL8TimeSeries *hplus;
	REAL8TimeSeries *hcross;
	lalDebugLevel = 7;
	XLALSimInspiralTaylorT1PN(&hplus, &hcross, &tc, phic, deltaT, m1, m2, f_min, fRef, r, i, lambda1, lambda2, interactionFlags, O);
	LALDPrintTimeSeries(hplus, "hp.dat");
	LALDPrintTimeSeries(hcross, "hc.dat");
	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);
	LALCheckMemoryLeaks();
	return 0;
}
#endif
