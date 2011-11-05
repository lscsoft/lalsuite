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
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/LALComplex.h>
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

NRCSID(LALSIMINSPIRALTAYLORT1C, "$Id$");

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
	REAL8 nu,m;
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
 * This function is used in the call to the GSL integrator.
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
 * Set up the expnCoeffsTaylorT1 and expnFuncTaylorT1 structures for
 * generating a TaylorT1 waveform and select the post-newtonian
 * functions corresponding to the desired order. 
 *
 * Inputs given in SI units.
 */
static int 
XLALSimInspiralTaylorT1Setup(
    expnCoeffsTaylorT1 *ak,		/**< coefficients for TaylorT1 evolution [modified] */
    expnFuncTaylorT1 *f,		/**< functions for TaylorT1 evolution [modified] */
    REAL8 m1,				/**< mass of companion 1 */
    REAL8 m2,				/**< mass of companion 2 */
    int O				/**< twice post-Newtonian order */
)
{
    ak->m = m1 + m2;
    REAL8 mu = m1 * m2 / ak->m;
    ak->nu = mu/ak->m;

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
		REAL8TimeSeries **v,   /**< post-Newtonian parameter [returned] */
		REAL8TimeSeries **phi, /**< orbital phase [returned] */
		LIGOTimeGPS *tc,       /**< coalescence time */
		REAL8 phic,            /**< coalescence phase */
		REAL8 deltaT,          /**< sampling interval */
		REAL8 m1,              /**< mass of companion 1 */
		REAL8 m2,              /**< mass of companion 2 */
		REAL8 f_min,           /**< start frequency */
		int O                  /**< twice post-Newtonian order */
		)
{
	const UINT4 blocklen = 1024;
	const REAL8 visco = 1./sqrt(6.);
	XLALSimInspiralTaylorT1PNEvolveOrbitParams params;
	expnFuncTaylorT1 expnfunc;
	expnCoeffsTaylorT1 ak;

	if(XLALSimInspiralTaylorT1Setup(&ak,&expnfunc,m1,m2,O))
		XLAL_ERROR(XLAL_EFUNC);

	params.flux=expnfunc.flux;
	params.dEdv=expnfunc.dEnergy;
	params.ak=ak;

	REAL8 E;
	UINT4 j;
	double y[2];
	double yerr[2];
	const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
	gsl_odeiv_step *s;
	gsl_odeiv_system sys;

	/* setup ode system */
	sys.function = XLALSimInspiralTaylorT1PNEvolveOrbitIntegrand;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = &params;

	/* allocate memory */
	*v = XLALCreateREAL8TimeSeries( "ORBITAL_VELOCITY_PARAMETER", tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	*phi = XLALCreateREAL8TimeSeries( "ORBITAL_PHASE", tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	if ( !v || !phi )
		XLAL_ERROR(XLAL_EFUNC);

	y[0] = (*v)->data->data[0] = cbrt(LAL_PI*LAL_G_SI*ak.m*f_min)/LAL_C_SI;
	y[1] = (*phi)->data->data[0] = 0.;
	E = expnfunc.energy(y[0],&(ak.akdEF));
	if (XLALIsREAL8FailNaN(E))
		XLAL_ERROR(XLAL_EFUNC);
	j = 0;

	s = gsl_odeiv_step_alloc(T, 2);
	while (1) {
		REAL8 dE;
		++j;
		gsl_odeiv_step_apply(s, j*deltaT, deltaT, y, yerr, NULL, NULL, &sys);
		/* MECO termination condition */
		dE = -E;
		dE += E = expnfunc.energy(y[0],&(ak.akdEF));
		if (XLALIsREAL8FailNaN(E))
			XLAL_ERROR(XLAL_EFUNC);
		if ( dE > 0.0 ) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at MECO\n", __func__);
			break;
		}
		/* ISCO termination condition for quadrupole, 1pN, 2.5pN */
		if ( (O == 0 || O == 1 || O == 2 || O == 5 || O == 7) && y[0] > visco ) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at ISCO\n", __func__);
			break;
		}
		if ( j >= (*v)->data->length ) {
			if ( ! XLALResizeREAL8TimeSeries(*v, 0, (*v)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*phi, 0, (*phi)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
		}
		(*v)->data->data[j] = y[0];
		(*phi)->data->data[j] = y[1];
	}
	gsl_odeiv_step_free(s);

	/* make the correct length */
	if ( ! XLALResizeREAL8TimeSeries(*v, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*phi, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);

	/* adjust to correct tc and phic */
	XLALGPSAdd(&(*v)->epoch, -1.0*j*deltaT);
	XLALGPSAdd(&(*phi)->epoch, -1.0*j*deltaT);

	/* phi here is the orbital phase = 1/2 * GW phase.
	 * End GW phase specified on command line.
	 * Adjust phase so phi = phi_end/2 at the end */

	phic /= 2.;
	phic -= (*phi)->data->data[j-1];
	for (j = 0; j < (*phi)->data->length; ++j)
		(*phi)->data->data[j] += phic;

	return (int)(*v)->data->length;
}


/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine allows the user to specify different pN orders
 * for phasing calcuation vs. amplitude calculations.
 */
int XLALSimInspiralTaylorT1PNGenerator(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 x0,                 /**< tail-term gauge choice thing (if you don't know, just set it to zero) */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int amplitudeO,           /**< twice post-Newtonian amplitude order */
	       	int phaseO                /**< twice post-Newtonian phase order */
		)
{
	REAL8TimeSeries *v;
	REAL8TimeSeries *phi;
	int status;
	int n;
	n = XLALSimInspiralTaylorT1PNEvolveOrbit(&v, &phi, tc, phic, deltaT, m1, m2, f_min, phaseO);
	if ( n < 0 )
		XLAL_ERROR(XLAL_EFUNC);
	status = XLALSimInspiralPNPolarizationWaveforms(hplus, hcross, v, phi, x0, m1, m2, r, i, amplitudeO);
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(v);
	if ( status < 0 )
		XLAL_ERROR(XLAL_EFUNC);
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
int XLALSimInspiralTaylorT1PN(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian order */
		)
{
	/* set x0=0 to ignore log terms */
	return XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, f_min, r, i, O, O);
}


/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 *
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 *
 * Log terms in amplitudes are ignored.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT1PNRestricted(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian phase order */
		)
{
	/* use Newtonian order for amplitude */
	/* set x0=0 to ignore log terms */
	return XLALSimInspiralTaylorT1PNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, f_min, r, i, 0, O);
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
	int O = -1;
	REAL8TimeSeries *hplus;
	REAL8TimeSeries *hcross;
	lalDebugLevel = 7;
	XLALSimInspiralTaylorT1PN(&hplus, &hcross, &tc, phic, deltaT, m1, m2, f_min, r, i, O);
	LALDPrintTimeSeries(hplus, "hp.dat");
	LALDPrintTimeSeries(hcross, "hc.dat");
	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);
	LALCheckMemoryLeaks();
	return 0;
}
#endif
