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

#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

NRCSID(LALSIMINSPIRALTAYLORT4C, "$Id$");

/**
 * This structure contains the intrinsic parameters and post-newtonian 
 * co-efficients for the energy and angular acceleration expansions. 
 * These are computed by XLALSimInspiralTaylorT4Setup routine.
 */

typedef struct
tagexpnCoeffsTaylorT4 {
   
   /* */
   REAL8 lambda,theta;
   /* Taylor expansion coefficents in domega/dt*/
   REAL8 aatN,aat2,aat3,aat4,aat5,aat6a,aat6b,aat7;
   /* Taylor expansion coefficents for energy*/
   REAL8 etN,et2,et4,et6;

    /*Angular velocity co-efficient*/
    REAL8 av;

   /* symmetric mass ratio, total mass, component masses*/
   REAL8 nu,m,m1,m2,mu;

   
}expnCoeffsTaylorT4;

typedef REAL8 (SimInspiralEnergy4)(
   REAL8 x, /**< post-Newtonian parameter */
   expnCoeffsTaylorT4 *ak
);

typedef REAL8 (SimInspiralAngularAcceleration4)(
   REAL8 x, /**< post-Newtonian parameter */
   expnCoeffsTaylorT4 *ak
);

/**
 * This strucuture contains pointers to the functions for calculating
 * the post-newtonian terms at the desired order. They can be set by
 * XLALSimInspiralTaylorT4Setup by passing an appropriate PN order. 
 */

typedef struct
tagexpnFuncTaylorT4
{
   SimInspiralEnergy4 *energy4;
   SimInspiralAngularAcceleration4 *angacc4;
} expnFuncTaylorT4;

/**
 * Computes the orbital energy at a fixed frequency and pN order.
 *
 * Implements Equation (152) of
 * Luc Blanchet,
 * "Gravitational Radiation from Post-Newtonian Sources and Inspiralling
 * Compact Binaries",
 * <a href="http://www.livingreviews.org/lrr-2006-4">lrr-2006-4</a>.
 *
 * This is the same as Equation (10) (where the spin of the objects
 * is zero) of:
 * Yi Pan, Alessandra Buonanno, Yanbei Chen, and Michele Vallisneri,
 * "A physical template family for gravitational waves from precessing
 * binaries of spinning compact objects: Application to single-spin binaries"
 * <a href="http://arxiv.org/abs/gr-qc/0310034">arXiv:gr-qc/0310034v3 </a>.
 * Note: this equation is actually dx/dt rather than (domega/dt)/(omega)^2
 * so the leading coefficient is different.
 */

static REAL8 
XLALSimInspiralEnergy4_4PN(
        REAL8 x,  /**< post-Newtonian parameter */
        expnCoeffsTaylorT4 *ak /**< PN co-efficients and intrinsic parameters */
)
{
    REAL8 ans = 0;

    ans += ak->etN
    //1PN
    +ak->et2*x
    //2PN
    +ak->et4*x*x;

    ans*=x;

    return ans;
}

static REAL8 
XLALSimInspiralEnergy4_6PN(
        REAL8 x,  /**< post-Newtonian parameter */
        expnCoeffsTaylorT4 *ak /**< PN co-efficients and intrinsic parameters */
)
{

    REAL8 ans = 0;

    ans += ak->etN
    //1PN
    +ak->et2*x
    //2PN
    +ak->et4*x*x
    //3PN
    +ak->et6*x*x*x;
    
    ans*=x;
    
    return ans;
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
 * <a href="http://arxiv.org/abs/gr-qc/0310034">arXiv:gr-qc/0310034v3 </a>.
 *
 * Note: this equation is actually dx/dt rather than (domega/dt)/(omega)^2
 * so the leading coefficient is different.  Also, this function applies
 * for non-spinning objects.
 *
 * The overall co-efficient with nu=0.25 is the same as to Equation (45) of
 * Michael Boyle, Duncan A. Brown, Lawrence E. Kidder, Abdul H. Mroue,
 * Harald P. Pfeiﬀer, Mark A. Scheel, Gregory B. Cook, and Saul A. Teukolsky
 * "High-accuracy comparison of numerical relativity simulations with
 * post-Newtonian expansions"
 * <a href="http://arxiv.org/abs/0710.0158v2">arXiv:0710.0158v2</a>.
 */
 
static REAL8 
XLALSimInspiralAngularAcceleration4_4PN(
  REAL8       x, /**< post-Newtonian parameter */
  expnCoeffsTaylorT4 *ak /**< PN co-efficients and intrinsic parameters */
)
{
    REAL8 ans;
    REAL8 y,x2;
    
    y=sqrt(x);
    x2=x*x;
    
    ans=ak->aatN
    /*1PN*/
    + ak->aat2*x
    /*1.5PN*/
    + ak->aat3*x*y
    /*2PN*/
    + ak->aat4*x2;

    ans*=x2*x2*x;

    return ans;
}


static REAL8 
XLALSimInspiralAngularAcceleration4_5PN(
  REAL8       x, /**< post-Newtonian parameter */
  expnCoeffsTaylorT4 *ak /**< PN co-efficients and intrinsic parameters */
)
{
    REAL8 ans;
    REAL8 y,x2;
    
    y=sqrt(x);
    x2=x*x;
    
    ans= ak->aatN
    /*1PN*/
    + ak->aat2*x
    /*1.5PN*/
    + ak->aat3*x*y
    /*2PN*/
    + ak->aat4*x2
    /*2.5PN*/
    + ak->aat5*x2*y;

    ans*=x2*x2*x;

    return ans;
}

static REAL8 
XLALSimInspiralAngularAcceleration4_6PN(
  REAL8       x, /**< post-Newtonian parameter */
  expnCoeffsTaylorT4 *ak /**< PN co-efficients and intrinsic parameters */
)
{
    REAL8 ans;
    REAL8 y,x2,x3;
    
    y=sqrt(x);
    x2=x*x;
    x3=x*x2;
    
    ans= ak->aatN
    /*1PN*/
    + ak->aat2*x
    /*1.5PN*/
    + ak->aat3*x*y
    /*2PN*/
    + ak->aat4*x2
    /*2.5PN*/
    + ak->aat5*x2*y
    /*3PN*/
    + (ak->aat6a - ak->aat6b*log(16.*x))*x3;

    ans*=x3*x2;

    return ans;
}

static REAL8 
XLALSimInspiralAngularAcceleration4_7PN(
  REAL8       x, /**< post-Newtonian parameter */
  expnCoeffsTaylorT4 *ak /**< PN co-efficients and intrinsic parameters */
)
{
    REAL8 ans;
    REAL8 y,x2,x3;
    
    y=sqrt(x);
    x2=x*x;
    x3=x*x2;
    
    ans= ak->aatN
    /*1PN*/
    + ak->aat2*x
    /*1.5PN*/
    + ak->aat3*x*y
    /*2PN*/
    + ak->aat4*x2
    /*2.5PN*/
    + ak->aat5*x2*y
    /*3PN*/
    + (ak->aat6a - ak->aat6b*log(16.*x))*x3
    /*3.5PN*/
    + ak->aat7*x3*y;

    ans*=x3*x2;

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
 * <a href="http://arxiv.org/abs/0710.0158v2">arXiv:0710.0158v2</a>.
 */
static REAL8 
XLALSimInspiralTaylorT4PNAngularVelocity(
    REAL8 x,  /**< post-Newtonian parameter */
    expnCoeffsTaylorT4 *ak /**< PN co-efficients and intrinsic parameters */
)
{
    return ak->av*pow(x, 1.5);
}

typedef struct
{
    REAL8 (*func)(REAL8 x, expnCoeffsTaylorT4 *ak);
    expnCoeffsTaylorT4 ak;
}XLALSimInspiralTaylorT4PNEvolveOrbitParams;

/** 
 * This function is used in the call to the GSL integrator.
 */
static int 
XLALSimInspiralTaylorT4PNEvolveOrbitIntegrand(double UNUSED t, const double y[], double ydot[], void *params)
{
	XLALSimInspiralTaylorT4PNEvolveOrbitParams* p = (XLALSimInspiralTaylorT4PNEvolveOrbitParams*)params;
	ydot[0] = p->func(y[0],&p->ak);
	ydot[1] = XLALSimInspiralTaylorT4PNAngularVelocity(y[0],&p->ak);
	t = 0.0;
	return GSL_SUCCESS;
}


/**
 * Set up the expnCoeffsTaylorT4 and expnFuncTaylorT4 structures for
 * generating a TaylorT4 waveform and select the post-newtonian
 * functions corresponding to the desired order. 
 *
 * Inputs given in SI units.
 */
static int 
XLALSimInspiralTaylorT4Setup(
    expnCoeffsTaylorT4 *ak,	/**< coefficients for TaylorT4 evolution [modified] */
    expnFuncTaylorT4 *f,	/**< functions for TaylorT4 evolution [modified] */
    REAL8 m1,		/**< mass of companion 1 */
    REAL8 m2,		/**< mass of companion 2 */
    int O			/**< twice post-Newtonian order */
)
{
    ak->m1 = m1;
    ak->m2 = m2;
    ak->m = ak->m1 + ak->m2;
    ak->mu = m1 * m2 / ak->m;
    ak->nu = ak->mu/ak->m;
    ak->lambda = -(1987.0/3080.0);
    ak->theta = (1039.0/4620.0);
    
    /* PN co-efficients for energy*/
    //N
    ak->etN=-0.5*ak->mu*pow(LAL_C_SI, 2.0);
    //1PN
    ak->et2=ak->etN*( -(3.0/4.0) - (1.0/12.0)*ak->nu );
    //2PN
    ak->et4=ak->etN*( -(27.0/8.0) + (19.0/8.0)*ak->nu - (1.0/24.0)*pow(ak->nu, 2.0) );
    //3PN
    ak->et6=ak->etN*( -(675.0/64.0) + ((209323.0/4032.0) - (205.0/96.0)*pow(LAL_PI, 2.0) - (110.0/9.0)*ak->lambda)*ak->nu - (155.0/96.0)*pow(ak->nu, 2.0) - (35.0/5184.0)*pow(ak->nu, 3.0) );
    
    /* PN co-efficients for angular acceleration*/
    //N
    ak->aatN = ((64.0*pow(LAL_C_SI, 3.0))/(5.0*LAL_G_SI*ak->m))*ak->nu;
    //1PN
    ak->aat2 = ak->aatN*( -(743.0/336.0) - (924.0/336.0)*ak->nu );
    //1.5PN
    ak->aat3 = ak->aatN*4.0*LAL_PI;
    //2PN
    ak->aat4 = ak->aatN*( (34103.0/18144.0) + (13661.0/2016.0)*ak->nu + (59.0/18.0)*pow(ak->nu, 2.0) );
    //2.5PN
    ak->aat5 = ak->aatN*( -(4159.0/672.0) - (15876.0/672.0)*ak->nu )*LAL_PI;
    //3PN
    ak->aat6a = ak->aatN*( ( (16447322263.0/139708800.0) - (1712.0/105.0)*LAL_GAMMA + (16.0/3.0)*pow(LAL_PI, 2.0) ) + ( -(273811877.0/1088640.0) + (451.0/48.0)*pow(LAL_PI, 2.0) - (88.0/3.0)*ak->theta )*ak->nu + (541.0/896.0)*ak->nu*ak->nu - (5605.0/2592.0)*ak->nu*ak->nu*ak->nu );
    ak->aat6b = ak->aatN*(856.0/105.0);
    //3.5PN
    ak->aat7 = ak->aatN*(-(4415.0/4032.0) + (358675.0/6048.0)*ak->nu + (91495.0/1512.0)*pow(ak->nu, 2.0))*LAL_PI;
    
    /*Angular velocity co-efficient*/
    ak->av = pow(LAL_C_SI, 3.0)/(LAL_G_SI*ak->m);
    
    switch (O)
    {
        case 0:
        case 1:
        case 2:
        case 3:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for PN order %d\n", __func__,O);
            XLAL_ERROR(__func__, XLAL_EINVAL);
            break;
        case 4:
            f->energy4 = &XLALSimInspiralEnergy4_4PN;
            f->angacc4 = &XLALSimInspiralAngularAcceleration4_4PN;
            break;
        case 5:
            f->energy4 = &XLALSimInspiralEnergy4_4PN;
            f->angacc4 = &XLALSimInspiralAngularAcceleration4_5PN;
            break;
        case 6:
            f->energy4 = &XLALSimInspiralEnergy4_6PN;
            f->angacc4 = &XLALSimInspiralAngularAcceleration4_6PN;
            break;
        case 7:
            f->energy4 = &XLALSimInspiralEnergy4_6PN;
            f->angacc4 = &XLALSimInspiralAngularAcceleration4_7PN;
            break;
        case 8:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for PN order %d\n", __func__,O);
            XLAL_ERROR(__func__, XLAL_EINVAL);
            break;
        default:
            XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
            XLAL_ERROR(__func__, XLAL_EINVAL);
    }
  
  return 0;
}


/**
 * Evolves a post-Newtonian orbit using the Taylor T4 method.
 *
 * See:
 * Michael Boyle, Duncan A. Brown, Lawrence E. Kidder, Abdul H. Mroue,
 * Harald P. Pfeiﬀer, Mark A. Scheel, Gregory B. Cook, and Saul A. Teukolsky
 * "High-accuracy comparison of numerical relativity simulations with
 * post-Newtonian expansions"
 * <a href="http://arxiv.org/abs/0710.0158v2">arXiv:0710.0158v2</a>.
 */
int XLALSimInspiralTaylorT4PNEvolveOrbit(
		REAL8TimeSeries **x,   /**< post-Newtonian parameter [returned] */
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
	static const char *func = "XLALSimInspiralTaylorT4PNEvolveOrbit";
	const UINT4 blocklen = 1024;
	const REAL8 xisco = 1./6.;
	XLALSimInspiralTaylorT4PNEvolveOrbitParams params;
    expnFuncTaylorT4 expnfunc;
    expnCoeffsTaylorT4 ak;
    
    if(XLALSimInspiralTaylorT4Setup(&ak,&expnfunc,m1,m2,O))
        XLAL_ERROR(__func__, XLAL_EFUNC);
    
    params.func=expnfunc.angacc4;
    params.ak=ak;
    
	REAL8 E;
	UINT4 j;
	double y[2];
	double yerr[2];
	const gsl_odeiv_step_type *T = gsl_odeiv_step_rk4;
	gsl_odeiv_step *s;
	gsl_odeiv_system sys;

	/* setup ode system */
	sys.function = XLALSimInspiralTaylorT4PNEvolveOrbitIntegrand;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = &params;

	/* allocate memory */
	*x = XLALCreateREAL8TimeSeries( "ORBITAL_FREQUENCY_PARAMETER", tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	*phi = XLALCreateREAL8TimeSeries( "ORBITAL_PHASE", tc, 0., deltaT, &lalDimensionlessUnit, blocklen );
	if ( !x || !phi )
		XLAL_ERROR(func, XLAL_EFUNC);

	y[0] = (*x)->data->data[0] = pow(LAL_PI*LAL_G_SI*ak.m*f_min/pow(LAL_C_SI,3.), 2./3.);
	y[1] = (*phi)->data->data[0] = 0.;
	E = expnfunc.energy4(y[0],&ak);
	if (XLALIsREAL8FailNaN(E))
		XLAL_ERROR(func, XLAL_EFUNC);
	j = 0;

	s = gsl_odeiv_step_alloc(T, 2);
	while (1) {
		REAL8 dE;
		++j;
		gsl_odeiv_step_apply(s, j*deltaT, deltaT, y, yerr, NULL, NULL, &sys);
		/* MECO termination condition */
		dE = -E;
		dE += E = expnfunc.energy4(y[0],&ak);
		if (XLALIsREAL8FailNaN(E))
			XLAL_ERROR(func, XLAL_EFUNC);
		if ( dE > 0.0 ) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at MECO\n", func);
			break;
		}
		/* ISCO termination condition for quadrupole, 1pN, 2.5pN */
		if ( (O == 0 || O == 1 || O == 2 || O == 5 || O == 7) && y[0] > xisco ) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at ISCO\n", func);
			break;
		}
		if ( j >= (*x)->data->length ) {
			if ( ! XLALResizeREAL8TimeSeries(*x, 0, (*x)->data->length + blocklen) )
				XLAL_ERROR(func, XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*phi, 0, (*phi)->data->length + blocklen) )
				XLAL_ERROR(func, XLAL_EFUNC);
		}
		(*x)->data->data[j] = y[0];
		(*phi)->data->data[j] = y[1];
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
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine allows the user to specify different pN orders
 * for phasing calcuation vs. amplitude calculations.
 */
int XLALSimInspiralTaylorT4PNGenerator(
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
	static const char *func = "XLALSimInspiralPNGenerator";
	REAL8TimeSeries *x;
	REAL8TimeSeries *phi;
	int status;
	int n;
	n = XLALSimInspiralTaylorT4PNEvolveOrbit(&x, &phi, tc, phic, deltaT, m1, m2, f_min, phaseO);
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
int XLALSimInspiralTaylorT4PN(
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
	return XLALSimInspiralTaylorT4PNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, f_min, r, i, O, O);
}


/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 *
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 *
 * Log terms in amplitudes are ignored.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT4PNRestricted(
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
	return XLALSimInspiralTaylorT4PNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, f_min, r, i, 0, O);
}


#if 0
#include <lal/PrintFTSeries.h>
extern int lalDebugLevel;
int main(void)
{
	LIGOTimeGPS tc = LIGOTIMEGPSZERO;
	REAL8 phic = 0.0;
	REAL8 deltaT = 1.0/4096.;
	REAL8 m1 = 5.*LAL_MSUN_SI;
	REAL8 m2 = 5.*LAL_MSUN_SI;
	REAL8 r = 1e6*LAL_PC_SI;
	REAL8 i = 0.5*LAL_PI;
	REAL8 f_min = 40.0;
	int O = 7;
	REAL8TimeSeries *hplus;
	REAL8TimeSeries *hcross;
	lalDebugLevel = 7;
	XLALSimInspiralTaylorT4PNGenerator(&hplus, &hcross, &tc, phic,0., deltaT, m1, m2, f_min, r, i,0, O);
	LALDPrintTimeSeries(hplus, "hp.dat");
	LALDPrintTimeSeries(hcross, "hc.dat");
	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);
	LALCheckMemoryLeaks();
	return 0;
}
#endif
