/*
 * Copyright (C) 2011 D. Keppel
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
#include <lal/FindRoot.h>
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

/**
 * This structure contains the intrinsic parameters and post-newtonian 
 * co-efficients for the denergy/dv and flux expansions. 
 * These are computed by XLALSimInspiralTaylorT1Setup routine.
 */

typedef struct
{
	/* coefficients for the TaylorEt phase evolution */
	REAL8 phN, ph1, ph2, ph3;

	/* coefficients for the TaylorEt zeta evolution */
	REAL8 zN, z2, z3, z4, z5, z6, z6l, z7;

	/* coefficients for the TaylorEt v(zeta) function */
	REAL8 v1, v2, v3;
} expnCoeffsTaylorEt;

typedef REAL8 (SimInspiralTaylorEtdPhase)(
	REAL8 zeta, /**< post-Newtonian parameter */
	expnCoeffsTaylorEt *ak
);

typedef REAL8 (SimInspiralTaylorEtdZeta)(
	REAL8 zeta, /**< post-Newtonian parameter */
	expnCoeffsTaylorEt *ak
);

typedef REAL8 (SimInspiralTaylorEtVOfZeta)(
	REAL8 zeta, /**< post-Newtonian parameter */
	expnCoeffsTaylorEt *ak
);

/**
 * This strucuture contains pointers to the functions for calculating
 * the post-newtonian terms at the desired order. They can be set by
 * XLALSimInspiralTaylorT1Setup by passing an appropriate PN order. 
 */

typedef struct
{
	SimInspiralTaylorEtdPhase *dphase;
	SimInspiralTaylorEtdZeta *dzeta;
	SimInspiralTaylorEtVOfZeta *vOfZeta;
} expnFuncTaylorEt;


/**
 *Computes v(zeta) for a post-Newtonian inspiral.
 *
 * Implements the square root of Equation 3.11 of: Alessandra Buonanno, Bala R
 * Iyer, Evan Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of
 * post-Newtonian templates for compact binary inspiral signals in
 * gravitational-wave detectors", Phys. Rev. D 80, 084043 (2009),
 * arXiv:0907.0700v1
 */

static REAL8
XLALSimInspiralTaylorEtVOfZeta_0PN(
		REAL8 zeta,			/**< post-Newtonian parameter */
		expnCoeffsTaylorEt UNUSED *ak	/**< coefficients for TaylorEt evolution */
		)
{
	return sqrt(zeta);
}

static REAL8
XLALSimInspiralTaylorEtVOfZeta_2PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	return sqrt(zeta * (1.0
		+ ak->v1 * zeta));
}

static REAL8
XLALSimInspiralTaylorEtVOfZeta_4PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta2 = zeta*zeta;
	return sqrt(zeta * (1.0
		+ ak->v1 * zeta
		+ ak->v2 * zeta2));
}

static REAL8
XLALSimInspiralTaylorEtVOfZeta_6PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta2 = zeta*zeta;
	REAL8 zeta3 = zeta2*zeta;
	return sqrt(zeta * (1.0
		+ ak->v1 * zeta
		+ ak->v2 * zeta2
		+ ak->v3 * zeta3));
}


/**
 * Computes the rate of increase of the phase for a post-Newtonian
 * inspiral.
 *
 * Implements Equation 3.13a of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8
XLALSimInspiralTaylorEtPhasing_0PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta3_2 = pow(zeta, 3.0/2.0);
	return ak->phN * zeta3_2;
}

static REAL8
XLALSimInspiralTaylorEtPhasing_2PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta3_2 = pow(zeta, 3.0/2.0);
	return ak->phN * zeta3_2 * (1.0
		+ ak->ph1 * zeta);
}

static REAL8
XLALSimInspiralTaylorEtPhasing_4PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta3_2 = pow(zeta, 3.0/2.0);
	REAL8 zeta2 = zeta*zeta;
	return ak->phN * zeta3_2 * (1.0
		+ ak->ph1 * zeta
		+ ak->ph2 * zeta2);
}

static REAL8
XLALSimInspiralTaylorEtPhasing_6PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta3_2 = pow(zeta, 3.0/2.0);
	REAL8 zeta2 = zeta*zeta;
	REAL8 zeta3 = zeta*zeta2;
	return ak->phN * zeta3_2 * (1.0
		+ ak->ph1 * zeta
		+ ak->ph2 * zeta2
		+ ak->ph3 * zeta3);
}


/**
 * Computes the rate of increase of zeta for a post-Newtonian
 * inspiral.
 *
 * Implements Equation 3.13b of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */

static REAL8
XLALSimInspiralTaylorEtZeta_0PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta5 = pow(zeta, 5.0);
	return ak->zN * zeta5;
}

static REAL8
XLALSimInspiralTaylorEtZeta_2PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta5 = pow(zeta, 5.0);
	return ak->zN * zeta5 * (1.0
		+ ak->z2 * zeta);
}

static REAL8
XLALSimInspiralTaylorEtZeta_3PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta1_2 = sqrt(zeta);
	REAL8 zeta5 = pow(zeta, 5.0);
	return ak->zN * zeta5 * (1.0
		+ ak->z2 * zeta
		+ ak->z3 * zeta * zeta1_2);
}

static REAL8
XLALSimInspiralTaylorEtZeta_4PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta2 = zeta*zeta;
	REAL8 zeta1_2 = sqrt(zeta);
	REAL8 zeta5 = pow(zeta, 5.0);
	return ak->zN * zeta5 * (1.0
		+ ak->z2 * zeta
		+ ak->z3 * zeta * zeta1_2
		+ ak->z4 * zeta2);
}

static REAL8
XLALSimInspiralTaylorEtZeta_5PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta2 = zeta*zeta;
	REAL8 zeta1_2 = sqrt(zeta);
	REAL8 zeta5 = pow(zeta, 5.0);
	return ak->zN * zeta5 * (1.0
		+ ak->z2 * zeta
		+ ak->z3 * zeta * zeta1_2
		+ ak->z4 * zeta2
		+ ak->z5 * zeta2 * zeta1_2);
}

static REAL8
XLALSimInspiralTaylorEtZeta_6PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta2 = zeta*zeta;
	REAL8 zeta3 = zeta2*zeta;
	REAL8 zeta1_2 = sqrt(zeta);
	REAL8 zeta5 = zeta2*zeta3;
	return ak->zN * zeta5 * (1.0
		+ ak->z2 * zeta
		+ ak->z3 * zeta * zeta1_2
		+ ak->z4 * zeta2
		+ ak->z5 * zeta2 * zeta1_2
		+ (ak->z6 + ak->z6l * log(16.0 * zeta)) * zeta3);
}

static REAL8
XLALSimInspiralTaylorEtZeta_7PN(
		REAL8 zeta,		/**< post-Newtonian parameter */
		expnCoeffsTaylorEt *ak	/**< coefficients for TaylorEt evolution */
		)
{
	REAL8 zeta2 = zeta*zeta;
	REAL8 zeta3 = zeta2*zeta;
	REAL8 zeta1_2 = sqrt(zeta);
	REAL8 zeta5 = zeta2*zeta3;
	return ak->zN * zeta5 * (1.0
		+ ak->z2 * zeta
		+ ak->z3 * zeta * zeta1_2
		+ ak->z4 * zeta2
		+ ak->z5 * zeta2 * zeta1_2
		+ (ak->z6 + ak->z6l * log(16.0 * zeta)) * zeta3
		+ ak->z7 * zeta3 * zeta1_2);
}


typedef struct
{
	REAL8 (*dphase)(REAL8 zeta, expnCoeffsTaylorEt *ak);
	REAL8 (*dzeta)(REAL8 zeta, expnCoeffsTaylorEt *ak);
	expnCoeffsTaylorEt ak;
}XLALSimInspiralTaylorEtPNEvolveOrbitParams;

/**
 * This function is used in the zeta intialization.
 */
typedef struct
{
	REAL8 f_target;
	XLALSimInspiralTaylorEtPNEvolveOrbitParams *eoparams;
}XLALSimInspiralTaylorEtPhasingWrapperParams;

static REAL8
XLALSimInspiralTaylorEtPhasingWrapper(
		REAL8 zeta,
		void *params
		)
{
	XLALSimInspiralTaylorEtPhasingWrapperParams *wrapperparams = (XLALSimInspiralTaylorEtPhasingWrapperParams*) params;
	return wrapperparams->f_target - wrapperparams->eoparams->dphase(zeta, &(wrapperparams->eoparams->ak));
}

/** 
 * This function is used in the call to the GSL integrator.
 */
static int 
XLALSimInspiralTaylorEtPNEvolveOrbitIntegrand(double UNUSED t, const double y[], double ydot[], void *params)
{
	XLALSimInspiralTaylorEtPNEvolveOrbitParams* p = (XLALSimInspiralTaylorEtPNEvolveOrbitParams*)params;
	ydot[0] = p->dzeta(y[0],&p->ak);
	ydot[1] = p->dphase(y[0],&p->ak);
	return GSL_SUCCESS;
}


/**
 * Set up the expnCoeffsTaylorEt and expnFuncTaylorEt structures for
 * generating a TaylorEt waveform and select the post-newtonian
 * functions corresponding to the desired order. 
 *
 * Inputs given in SI units.
 */
static int 
XLALSimInspiralTaylorEtSetup(
		expnCoeffsTaylorEt *ak,	/**< coefficients for TaylorEt evolution [modified] */
		expnFuncTaylorEt *f,	/**< functions for TaylorEt evolution [modified] */
		REAL8 m1,		/**< mass of companion 1 */
		REAL8 m2,		/**< mass of companion 2 */
		int O			/**< twice post-Newtonian order */
		)
{
	REAL8 m = m1 + m2;
	REAL8 nu = m1 * m2 / (m * m);
	m *= LAL_G_SI * pow(LAL_C_SI, -3.0);

	/* Taylor co-efficients for dPhase/dt. */
	ak->phN = XLALSimInspiralTaylorEtPhasing_0PNCoeff(m);
	ak->ph1 = XLALSimInspiralTaylorEtPhasing_2PNCoeff(nu);
	ak->ph2 = XLALSimInspiralTaylorEtPhasing_4PNCoeff(nu);
	ak->ph3 = XLALSimInspiralTaylorEtPhasing_6PNCoeff(nu);

	/* Taylor co-efficients for dZeta/dt. */
	ak->zN = XLALSimInspiralTaylorEtZeta_0PNCoeff(m, nu);
	ak->z2 = XLALSimInspiralTaylorEtZeta_2PNCoeff(nu);
	ak->z3 = XLALSimInspiralTaylorEtZeta_3PNCoeff(nu);
	ak->z4 = XLALSimInspiralTaylorEtZeta_4PNCoeff(nu);
	ak->z5 = XLALSimInspiralTaylorEtZeta_5PNCoeff(nu);
	ak->z6 = XLALSimInspiralTaylorEtZeta_6PNCoeff(nu);
	ak->z6l = XLALSimInspiralTaylorEtZeta_6PNLogCoeff(nu);
	ak->z7 = XLALSimInspiralTaylorEtZeta_7PNCoeff(nu);
    
	/* Taylor co-efficients for v(zeta). */
	ak->v1 = XLALSimInspiralTaylorEtVOfZeta_2PNCoeff(nu);
	ak->v2 = XLALSimInspiralTaylorEtVOfZeta_4PNCoeff(nu);
	ak->v3 = XLALSimInspiralTaylorEtVOfZeta_6PNCoeff(nu);

	switch (O)
	{
		case 0:
			f->dphase = &XLALSimInspiralTaylorEtPhasing_0PN;
			f->dzeta = &XLALSimInspiralTaylorEtZeta_0PN;
			f->vOfZeta = &XLALSimInspiralTaylorEtVOfZeta_0PN;
			break;
		case 1:
			XLALPrintError("XLAL Error - %s: PN approximant not supported for PN order %d\n", __func__,O);
			XLAL_ERROR(XLAL_EINVAL);
			break;
		case 2:
			f->dphase = &XLALSimInspiralTaylorEtPhasing_2PN;
			f->dzeta = &XLALSimInspiralTaylorEtZeta_2PN;
			f->vOfZeta = &XLALSimInspiralTaylorEtVOfZeta_2PN;
			break;
		case 3:
			f->dphase = &XLALSimInspiralTaylorEtPhasing_2PN;
			f->dzeta = &XLALSimInspiralTaylorEtZeta_3PN;
			f->vOfZeta = &XLALSimInspiralTaylorEtVOfZeta_2PN;
			break;
		case 4:
			f->dphase = &XLALSimInspiralTaylorEtPhasing_4PN;
			f->dzeta = &XLALSimInspiralTaylorEtZeta_4PN;
			f->vOfZeta = &XLALSimInspiralTaylorEtVOfZeta_4PN;
			break;
		case 5:
			f->dphase = &XLALSimInspiralTaylorEtPhasing_4PN;
			f->dzeta = &XLALSimInspiralTaylorEtZeta_5PN;
			f->vOfZeta = &XLALSimInspiralTaylorEtVOfZeta_4PN;
			break;
		case 6:
			f->dphase = &XLALSimInspiralTaylorEtPhasing_6PN;
			f->dzeta = &XLALSimInspiralTaylorEtZeta_6PN;
			f->vOfZeta = &XLALSimInspiralTaylorEtVOfZeta_6PN;
			break;
		case 7:
		case -1:
			f->dphase = &XLALSimInspiralTaylorEtPhasing_6PN;
			f->dzeta = &XLALSimInspiralTaylorEtZeta_7PN;
			f->vOfZeta = &XLALSimInspiralTaylorEtVOfZeta_6PN;
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
 * Evolves a post-Newtonian orbit using the Taylor Et method.
 *
 * See Section IIIE: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 */
int XLALSimInspiralTaylorEtPNEvolveOrbit(
		REAL8TimeSeries **V,   /**< post-Newtonian parameter [returned] */
		REAL8TimeSeries **phi, /**< orbital phase [returned] */
		REAL8 phic,            /**< orbital phase at end */
		REAL8 deltaT,          /**< sampling interval */
		REAL8 m1,              /**< mass of companion 1 */
		REAL8 m2,              /**< mass of companion 2 */
		REAL8 f_min,           /**< start frequency */
		int O                  /**< twice post-Newtonian order */
		)
{
	const UINT4 blocklen = 1024;
	const REAL8 visco = 1./sqrt(6.);
	XLALSimInspiralTaylorEtPNEvolveOrbitParams params;
	expnFuncTaylorEt expnfunc;
	expnCoeffsTaylorEt ak;

	if(XLALSimInspiralTaylorEtSetup(&ak,&expnfunc,m1,m2,O))
		XLAL_ERROR(XLAL_EFUNC);

	params.dphase = expnfunc.dphase;
	params.dzeta = expnfunc.dzeta;
	params.ak = ak;

	UINT4 j;
	LIGOTimeGPS tc = LIGOTIMEGPSZERO;
	double y[2];
	double yerr[2];
	const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
	gsl_odeiv_step *s;
	gsl_odeiv_system sys;

	/* setup ode system */
	sys.function = XLALSimInspiralTaylorEtPNEvolveOrbitIntegrand;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = &params;

	/* allocate memory */
	*V = XLALCreateREAL8TimeSeries( "ORBITAL_VELOCITY_PARAMETER", &tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	*phi = XLALCreateREAL8TimeSeries( "ORBITAL_PHASE", &tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	if ( !V || !phi )
		XLAL_ERROR(XLAL_EFUNC);

	XLALSimInspiralTaylorEtPhasingWrapperParams wrapperparams;
	wrapperparams.f_target = LAL_PI * f_min;
	wrapperparams.eoparams = &params;
	REAL8 v_min = cbrt(f_min * LAL_PI * LAL_G_SI * (m1 + m2))/LAL_C_SI;
	REAL8 xmax = 10.0 * v_min*v_min;
	REAL8 xacc = 1.0e-8;
	REAL8 xmin = 0.1 * v_min*v_min;

	y[0] = XLALDBisectionFindRoot(XLALSimInspiralTaylorEtPhasingWrapper, xmin, xmax, xacc, (void*) (&wrapperparams));
	y[1] = (*phi)->data->data[0] = 0.;
	(*V)->data->data[0] = expnfunc.vOfZeta(y[0], &ak);
	j = 0;

	s = gsl_odeiv_step_alloc(T, 2);
	while (1) {
		REAL8 tmpv;
		++j;
		gsl_odeiv_step_apply(s, j*deltaT, deltaT, y, yerr, NULL, NULL, &sys);
		tmpv = expnfunc.vOfZeta(y[0], &ak);
		if ( tmpv > visco ) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at ISCO\n", __func__);
			break;
		}
		if ( j >= (*V)->data->length ) {
			if ( ! XLALResizeREAL8TimeSeries(*V, 0, (*V)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*phi, 0, (*phi)->data->length + blocklen) )
				XLAL_ERROR(XLAL_EFUNC);
		}
		(*V)->data->data[j] = tmpv;
		(*phi)->data->data[j] = y[1];
	}
	gsl_odeiv_step_free(s);

	/* make the correct length */
	if ( ! XLALResizeREAL8TimeSeries(*V, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*phi, 0, j) )
		XLAL_ERROR(XLAL_EFUNC);

	/* adjust to correct time */
	XLALGPSAdd(&(*V)->epoch, -1.0*j*deltaT);
	XLALGPSAdd(&(*phi)->epoch, -1.0*j*deltaT);

	/* Shift phase so phi = phic at the end */
	phic -= (*phi)->data->data[j-1];
	for (j = 0; j < (*phi)->data->length; ++j)
		(*phi)->data->data[j] += phic;

	return (int)(*V)->data->length;
}


/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine allows the user to specify different pN orders
 * for phasing calcuation vs. amplitude calculations.
 */
int XLALSimInspiralTaylorEtPNGenerator(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	REAL8 phic,               /**< orbital phase at end */
	       	REAL8 v0,                 /**< tail-term gauge choice (default = 1) */
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
	REAL8TimeSeries *V;
	REAL8TimeSeries *phi;
	int status;
	int n;
	n = XLALSimInspiralTaylorEtPNEvolveOrbit(&V, &phi, phic, deltaT, m1, m2, f_min, phaseO);
	if ( n < 0 )
		XLAL_ERROR(XLAL_EFUNC);
	status = XLALSimInspiralPNPolarizationWaveforms(hplus, hcross, V, phi, v0, m1, m2, r, i, amplitudeO);
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(V);
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
 * Constant log term in amplitude set to 1.  This is a gauge choice.
 */
int XLALSimInspiralTaylorEtPN(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	REAL8 phic,               /**< orbital phase at end */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian order */
		)
{
	/* set v0 to default value 1 */
	return XLALSimInspiralTaylorEtPNGenerator(hplus, hcross, phic, 1.0, deltaT, m1, m2, f_min, r, i, O, O);
}


/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 *
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 *
 * Constant log term in amplitude set to 1.  This is a gauge choice.
 */
int XLALSimInspiralTaylorEtPNRestricted(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	REAL8 phic,               /**< orbital phase at end */
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
	/* set v0 to default value 1 */
	return XLALSimInspiralTaylorEtPNGenerator(hplus, hcross, phic, 1.0, deltaT, m1, m2, f_min, r, i, 0, O);
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
	XLALSimInspiralTaylorEtPN(&hplus, &hcross, &tc, phic, deltaT, m1, m2, f_min, r, i, O);
	LALDPrintTimeSeries(hplus, "hp.dat");
	LALDPrintTimeSeries(hcross, "hc.dat");
	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);
	LALCheckMemoryLeaks();
	return 0;
}
#endif
