/*
 * Copyright (C) 2012 P. Ajith 
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
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimInspiral.h>
#include <lal/VectorOps.h>
#include "check_series_macros.h"
#define LAL_PISQR 9.869604401089358

#define UNUSED(expr) do { (void)(expr); } while (0)

/* use error codes above 1024 to avoid conflicts with GSL */
#define LALSIMINSPIRAL_ST5_TEST_ENERGY 			1025
#define LALSIMINSPIRAL_ST5_TEST_VDOT 			1026
#define LALSIMINSPIRAL_ST5_TEST_COORDINATE 		1027
#define LALSIMINSPIRAL_ST5_TEST_VNAN			1028
#define LALSIMINSPIRAL_ST5_TEST_FREQBOUND 		1029
#define LALSIMINSPIRAL_ST5_DERIVATIVE_OMEGANONPOS 	1030

/* (2x) Highest available PN order - UPDATE IF NEW ORDERS ADDED!!*/
#define LAL_MAX_PN_ORDER 8
/* Number of variables used for precessing waveforms */
#define LAL_NUM_ST5_VARIABLES 10
/* absolute and relative tolerance for adaptive Runge-Kutta ODE integrator */
#define LAL_ST5_ABSOLUTE_TOLERANCE 1.e-12
#define LAL_ST5_RELATIVE_TOLERANCE 1.e-12

/**
 * Structure containing the non-dynamical coefficients needed
 * to evolve a spinning, precessing binary and produce a waveform.
 * This struct is passed to the 2 static functions below
 */

typedef struct SpinTaylorT5StructParams {
    double mass1;
    double mass2;
    double totalMass;
    double eta;
    double delta;
    int phaseO;
    double massRatio;
    double m1Sqr;
    double m2Sqr;
    double mSqr;
    double etaSqr;
    double deltaEta;
	double dEbF0nonSpin;
	double dEbF1nonSpin;
	double dEbF2nonSpin;
	double dEbF3nonSpin;
	double dEbF4nonSpin;
	double dEbF5nonSpin;
	double dEbF6nonSpin;
	double dEbF6LogNonSpin;
	double dEbF7nonSpin;
	double phiOrb0nonSpin;
    double phiOrb1nonSpin;
    double phiOrb2nonSpin;
	double phiOrb3nonSpin;
	double phiOrb4nonSpin;
	double phiOrb5nonSpin;
	double phiOrb6nonSpin;
	double phiOrb6LognonSpin;
    double phiOrb7nonSpin;
    double v0;
    double vMax;
	double phiRef; 
} SpinTaylorT5Params;

/* declaration of static functions */ 
static int XLALSimInspiralSpinTaylorT5Derivatives(
		double t, 
		const double y[],
		double dy[], 
		void *mParams
		);

static REAL8 XLALdEnergyByFluxSpinPrec(
		double v, 
		double *chi1, 
		double *chi2, 
        double *LNh, 
		SpinTaylorT5Params *params
		);

static int polarizationsInRadiationFrame(
		REAL8TimeSeries **hplus,  	/**< +-polarization waveform [returned] */
		REAL8TimeSeries **hcross, 	/**< x-polarization waveform [returned] */
		REAL8TimeSeries *V,       	/**< post-Newtonian parameter */
		REAL8TimeSeries *Phi,     	/**< orbital phase */
		REAL8TimeSeries *LNhxVec,  	/**< unit orbital ang. mom. x comp. */
		REAL8TimeSeries *LNhyVec,  	/**< unit orbital ang. mom. y comp. */
		REAL8TimeSeries *LNhzVec,  	/**< unit orbital ang. mom. z comp. */
		REAL8 m1, 					/**< mass1 (KG) */
		REAL8 m2, 					/**< mass2 (KG) */
		REAL8 r, 					/**< distance (meter) */
		REAL8 Theta					/**< angle between the initial total ang momentum and line of sight. */
	);

static int computeOrbitalPhase(
        SpinTaylorT5Params *params, 
		REAL8TimeSeries *V, 
		REAL8TimeSeries *LNhxVec, 
		REAL8TimeSeries *LNhyVec, 
		REAL8TimeSeries *LNhzVec, 
		REAL8TimeSeries *S1xVec, 
		REAL8TimeSeries *S1yVec, 
		REAL8TimeSeries *S1zVec, 
		REAL8TimeSeries *S2xVec, 
		REAL8TimeSeries *S2yVec, 
		REAL8TimeSeries *S2zVec, 
		REAL8TimeSeries *orbPhase);

static int spinTaylorT5Init(
		REAL8 m1, 
		REAL8 m2, 
		REAL8 fStart, 
		REAL8 deltaT, 
		REAL8 phiRef, 
		UINT4 phaseO, 
		SpinTaylorT5Params *mParams
	);

static int XLALSimInspiralSpinTaylorT5StoppingTest(
		double t, 
		const double y[],
		double dy[], 
		void *mParams
	);

static double dotProduct(double a[3], double b[3]);
static void crossProduct(double a[3], double b[3], double c[3]);
static void rotateVector(double a[3], double b[3], double phi, double theta, double psi);


/**
 * Generate time-domain generic spinnin PN waveforms in the SpinTaylorT5 approximaton.
 */
int XLALSimInspiralSpinTaylorT5 (
		REAL8TimeSeries **hplus,        /**< +-polarization waveform */
		REAL8TimeSeries **hcross,       /**< x-polarization waveform */
		REAL8 phiRef,                   /**< orbital phase at reference pt. */
		REAL8 deltaT,                   /**< sampling interval (s) */
		REAL8 m1,                       /**< mass of companion 1 (kg) */
		REAL8 m2,                       /**< mass of companion 2 (kg) */
		REAL8 fStart,                   /**< start GW frequency (Hz) */
		REAL8 r,                        /**< distance of source (m) */
		REAL8 s1x,                      /**< initial value of S1x */
		REAL8 s1y,                      /**< initial value of S1y */
		REAL8 s1z,                      /**< initial value of S1z */
		REAL8 s2x,                      /**< initial value of S2x */
		REAL8 s2y,                      /**< initial value of S2y */
		REAL8 s2z,                      /**< initial value of S2z */
		REAL8 incAngle, 				/**< inclination angle with J_ini */
		int phaseO,                     /**< twice PN phase order */
		int amplitudeO                  /**< twice PN amplitude order */
	) {


	UINT4 i, intStatus, lenReturn;
	LIGOTimeGPS tStart = LIGOTIMEGPSZERO;
	REAL8TimeSeries *orbPhase=NULL, *V=NULL, *LNhxVec=NULL, *LNhyVec=NULL, *LNhzVec=NULL; 
	REAL8TimeSeries *S1xVec=NULL, *S1yVec=NULL, *S1zVec=NULL, *S2xVec=NULL, *S2yVec=NULL, *S2zVec=NULL;
	ark4GSLIntegrator *integrator = NULL;     /* GSL integrator object */
    SpinTaylorT5Params *mParams;
    SpinTaylorT5Params SpinTaylorParameters;
    REAL8 LNh[3], S1[3], S2[3], J[3], Jh[3], Lmag, Jmag, v; 
    
    mParams = &SpinTaylorParameters;		/* structure containing various parameters and coefficients */
	REAL8 yinit[LAL_NUM_ST5_VARIABLES];     /* initial values of parameters */
    REAL8Array *yout=NULL;	 				/* time series of variables returned from integrator */

	/* initialize the mParams containing the the PN expansion parameters */
	if (spinTaylorT5Init(m1, m2, fStart, deltaT, phiRef, phaseO, mParams) != XLAL_SUCCESS) {
		XLALPrintError("XLAL Error - %s: Unable to initialize the mParams structure.\n", __func__); 
        return XLAL_FAILURE; 
	}
		
    /* Estimate length of waveform using Newtonian chirp time formula. */
	UINT4 dataLength = pow(2, ceil(log2(5*mParams->totalMass/(256.*pow(mParams->v0,8.)*mParams->eta)/deltaT)));

	/* Adjust tStart so last sample is at time=0 */
    XLALGPSAdd(&tStart, -1.0*(dataLength-1)*deltaT);

	/* allocate memory for vectors storing the PN parameter, orbital phase
	and the unit vector along the Newtonian orbital angular momentum  */
    V = XLALCreateREAL8TimeSeries( "PN Parameter", &tStart, 0.0, deltaT, &lalStrainUnit, dataLength);
    orbPhase = XLALCreateREAL8TimeSeries( "Orbital Phase", &tStart, 0.0, deltaT, &lalStrainUnit, dataLength);
    LNhxVec = XLALCreateREAL8TimeSeries( "Unit vec along Newt. ang. mom (X comp)", &tStart, 0., deltaT, &lalStrainUnit, dataLength);
    LNhyVec = XLALCreateREAL8TimeSeries( "Unit vec along Newt. ang. mom (Y comp)", &tStart, 0., deltaT, &lalStrainUnit, dataLength);
    LNhzVec = XLALCreateREAL8TimeSeries( "Unit vec along Newt. ang. mom (Z comp)", &tStart, 0., deltaT, &lalStrainUnit, dataLength);
    S1xVec = XLALCreateREAL8TimeSeries( "Spin ang mom of body 1 (X comp)", &tStart, 0., deltaT, &lalStrainUnit, dataLength);
    S1yVec = XLALCreateREAL8TimeSeries( "Spin ang mom of body 1 (Y comp)", &tStart, 0., deltaT, &lalStrainUnit, dataLength);
    S1zVec = XLALCreateREAL8TimeSeries( "Spin ang mom of body 1 (Z comp)", &tStart, 0., deltaT, &lalStrainUnit, dataLength);
    S2xVec = XLALCreateREAL8TimeSeries( "Spin ang mom of body 2 (X comp)", &tStart, 0., deltaT, &lalStrainUnit, dataLength);
    S2yVec = XLALCreateREAL8TimeSeries( "Spin ang mom of body 2 (Y comp)", &tStart, 0., deltaT, &lalStrainUnit, dataLength);
    S2zVec = XLALCreateREAL8TimeSeries( "Spin ang mom of body 2 (Z comp)", &tStart, 0., deltaT, &lalStrainUnit, dataLength);

	memset(V->data->data, 0, dataLength*sizeof(REAL8));
	memset(orbPhase->data->data, 0, dataLength*sizeof(REAL8));
	memset(LNhxVec->data->data, 0, dataLength*sizeof(REAL8));
	memset(LNhyVec->data->data, 0, dataLength*sizeof(REAL8));
	memset(LNhzVec->data->data, 0, dataLength*sizeof(REAL8));
	memset(S1xVec->data->data, 0, dataLength*sizeof(REAL8));
	memset(S1yVec->data->data, 0, dataLength*sizeof(REAL8));
	memset(S1zVec->data->data, 0, dataLength*sizeof(REAL8));
	memset(S2xVec->data->data, 0, dataLength*sizeof(REAL8));
	memset(S2yVec->data->data, 0, dataLength*sizeof(REAL8));
	memset(S2zVec->data->data, 0, dataLength*sizeof(REAL8));

    /* binary parameters in a coordinate system whose z axis is along the initial 
    orbital angular momentum  */
    LNh[0] = 0.;                		/* unit vector along Newt. ang. momtm  	*/
    LNh[1] = 0.;                		/* unit vector along Newt. ang. momtm  	*/
    LNh[2] = 1.;                		/* unit vector along Newt. ang. momtm 	*/ 
    S1[0]  = s1x*mParams->m1Sqr;     	/* spin vector of the first object		*/
    S1[1]  = s1y*mParams->m1Sqr;     	/* spin vector of the first object		*/
    S1[2]  = s1z*mParams->m1Sqr;     	/* spin vector of the first object		*/
    S2[0]  = s2x*mParams->m2Sqr;     	/* spin vector of the second object 	*/
    S2[1]  = s2y*mParams->m2Sqr;     	/* spin vector of the second object 	*/
    S2[2]  = s2z*mParams->m2Sqr;     	/* spin vector of the second object 	*/

    /* compute the initial total angular momentum  */
    Lmag = (mParams->eta*mParams->mSqr/mParams->v0)*(1. + (1.5+mParams->eta/6.)*mParams->v0*mParams->v0);     // magnitde of orb ang momentum (1PN)
    //Lmag = (mParams->eta*mParams->mSqr/mParams->v0);											/* magnitude of orb ang momentum (Newtonian) */
    for (i=0; i<3; i++) J[i] = Lmag*LNh[i] + S1[i] + S2[i];           // total angular momentum vec
    Jmag = sqrt(J[0]*J[0] + J[1]*J[1] + J[2]*J[2]);                   // magnitude of tot ang momentum 
    for (i=0; i<3; i++) Jh[i] = J[i]/Jmag;                            // unit vector along J 
    
	/* transform to the coordinate system defined by the initial total angular momentum */
	REAL8 JIniTheta = -acos(Jh[2]);
    REAL8 JIniPhi = -atan2(Jh[1], Jh[0]);

	/* if both objects are non-spinning, the coordinate system defined by the total angular
	momentum (= orbital angular momentum) will cause the dynamical variables to hit 
	coordinate singularities. Thus, for the case of non-spinning binaries, we avoid this
	situlation by choosing a coordinate system which is misaligned by 45 degs with the 
	total angular momentum direction. We also shift the inclination angle by this amount 
	so that, at the end, the waveforms are not affected by this choice of coordinate 
	system */
	if (JIniTheta == 0) {
		JIniTheta = LAL_PI/4.;
		incAngle = incAngle+JIniTheta; 
	}

	rotateVector(LNh, LNh, JIniPhi, JIniTheta, 0.);
    rotateVector(S1, S1, JIniPhi, JIniTheta, 0.);
    rotateVector(S2, S2, JIniPhi, JIniTheta, 0.);
    rotateVector(Jh, Jh, JIniPhi, JIniTheta, 0.);

    // initialise the evolution variables. 
    yinit[0]  = v = mParams->v0;    // PN expansion parameter
    yinit[1]  = LNh[0];    // unit vector along Newt. ang. momtm 
    yinit[2]  = LNh[1];    // unit vector along Newt. ang. momtm 
    yinit[3]  = LNh[2];    // unit vector along Newt. ang. momtm 
    yinit[4]  = S1[0];     // spin vector of the first object
    yinit[5]  = S1[1];     // spin vector of the first object
    yinit[6]  = S1[2];     // spin vector of the first object
    yinit[7]  = S2[0];     // spin vector of the second object 
    yinit[8]  = S2[1];     // spin vector of the second object 
    yinit[9]  = S2[2];     // spin vector of the second object 
	
    /* initialize the integrator */
    integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST5_VARIABLES,
            XLALSimInspiralSpinTaylorT5Derivatives,
            XLALSimInspiralSpinTaylorT5StoppingTest,
            LAL_ST5_ABSOLUTE_TOLERANCE, LAL_ST5_RELATIVE_TOLERANCE);

    if( !integrator ) {
        XLALPrintError("XLAL Error - %s: Cannot allocate integrator\n", __func__);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* stop the integration only when the test is true */
    integrator->stopontestonly = 0;

    /* run the ntegration; note: time is measured in \hat{t} = t / M */
    lenReturn = XLALAdaptiveRungeKutta4Hermite(integrator, (void *) mParams, yinit,
            0.0, dataLength*deltaT, deltaT, &yout);

    intStatus = integrator->returncode;
    XLALAdaptiveRungeKutta4Free(integrator);

    if (!lenReturn) {
        XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intStatus);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* Print warning about abnormal termination */
    if (intStatus != 0 && intStatus != LALSIMINSPIRAL_ST5_TEST_ENERGY
		&& intStatus != LALSIMINSPIRAL_ST5_TEST_FREQBOUND) {
			XLALPrintWarning("XLAL Warning - %s: integration terminated with code %d.\n \
			Waveform parameters were m1 = %e, m2 = %e, S1 = [%e,%e,%e], S2 = [%e,%e,%e], LNh = [%e,%e,%e]\n", 
				__func__, intStatus, m1/LAL_MSUN_SI, m2/LAL_MSUN_SI, S1[0], S1[1], S1[2], S2[0], 
				S2[1], S2[2], LNh[0], LNh[1], LNh[2]);
    }

	/* Copy dynamical variables from yout array to output time series.
	note that yout[0:dataLength=1] is time */
	for (i = 0; i<lenReturn; i++) {
        V->data->data[i] 	   = yout->data[lenReturn+ i]; 
        LNhxVec->data->data[i] = yout->data[2*lenReturn+ i];  
        LNhyVec->data->data[i] = yout->data[3*lenReturn+ i];
        LNhzVec->data->data[i] = yout->data[4*lenReturn+ i];
		S1xVec->data->data[i]  = yout->data[5*lenReturn+ i];  
        S1yVec->data->data[i]  = yout->data[6*lenReturn+ i];
        S1zVec->data->data[i]  = yout->data[7*lenReturn+ i];
		S2xVec->data->data[i]  = yout->data[8*lenReturn+ i];  
        S2yVec->data->data[i]  = yout->data[9*lenReturn+ i];
        S2zVec->data->data[i]  = yout->data[10*lenReturn+ i];
	}

	XLALDestroyREAL8Array(yout);

	/* compute the orbital phase */
	if (computeOrbitalPhase(mParams, V, LNhxVec, LNhyVec, LNhzVec, S1xVec, S1yVec, S1zVec, 
			S2xVec, S2yVec, S2zVec, orbPhase) != XLAL_SUCCESS) {
		XLALPrintError("XLAL Error - %s: Unable to compute the orbital phase .\n", __func__); 
        return XLAL_FAILURE; 
	}

	/* project the polarizations into the radiation frame */
	switch (amplitudeO) {
		case 0: 
			if (polarizationsInRadiationFrame(hplus, hcross, V, orbPhase, LNhxVec, 
					LNhyVec, LNhzVec, m1, m2, r, incAngle) != XLAL_SUCCESS) {
				XLALPrintError("XLAL Error - %s: Unable to project the waveforms into radiation frame.\n", __func__); 
				return XLAL_FAILURE; 
			}
		break; 
		default: 
			XLALPrintError("XLAL Error - %s: Polarizations are currently computed only in amplitudeO = 0\n", __func__);
        	return XLAL_FAILURE; 
	}

	/* clear memory */
	XLALDestroyREAL8TimeSeries(V);
	XLALDestroyREAL8TimeSeries(orbPhase);
	XLALDestroyREAL8TimeSeries(LNhxVec);
	XLALDestroyREAL8TimeSeries(LNhyVec);
	XLALDestroyREAL8TimeSeries(LNhzVec);
	XLALDestroyREAL8TimeSeries(S1xVec);
	XLALDestroyREAL8TimeSeries(S1yVec);
	XLALDestroyREAL8TimeSeries(S1zVec);
	XLALDestroyREAL8TimeSeries(S2xVec);
	XLALDestroyREAL8TimeSeries(S2yVec);
	XLALDestroyREAL8TimeSeries(S2zVec);

	return XLAL_SUCCESS;
}

/**
 * Evolution of dynamical variables in the SpinTaylorT5
 */
static int XLALSimInspiralSpinTaylorT5Derivatives(
		double t, 
		const double y[],
		double dydt[], 
		void *mParams
	) {

    REAL8 LNh[3], S1[3], S2[3], chi1[3], chi2[3], Omega1[3], Omega2[3];
    REAL8 dS1byDt[3], dS2byDt[3], dLNhByDt[3];
    REAL8 m, eta, q, delta, v, dEbyF, om0, v2, S1dotLNh, S2dotLNh, Lmag;
    UINT4 phaseO, i, onePFpnFlag = 1, onePNFlag = 1;
    SpinTaylorT5Params *params = (SpinTaylorT5Params*) mParams;
	UNUSED(t);

    /* binary parameters */
    m       = params->totalMass;    /* in seconds */
    eta     = params->eta;
    delta   = params->delta;
    phaseO  = params->phaseO;
    q       = params->massRatio;

    /* evolution variables. */
    v      = y[0];   /* PN expansion parameter: (m \omega_orb)^(1/3) */
    LNh[0] = y[1];   /* unit vector in the direction of Newt. ang. momentum (LNh_x) */
    LNh[1] = y[2];   /* unit vector in the direction of Newt. ang. momentum (LNh_y) */
    LNh[2] = y[3];   /* unit vector in the direction of Newt. ang. momentum (LNh_z) */
    S1[0]  = y[4];   /* spin vector of the first object */
    S1[1]  = y[5];   /* spin vector of the first object */
    S1[2]  = y[6];   /* spin vector of the first object */
    S2[0]  = y[7];   /* spin vector of the second object  */
    S2[1]  = y[8];   /* spin vector of the second object  */
    S2[2]  = y[9];   /* spin vector of the second object  */
    
	/* comptute the dimensionless spins */
	chi1[0] = S1[0]/params->m1Sqr; 
	chi1[1] = S1[1]/params->m1Sqr; 
	chi1[2] = S1[2]/params->m1Sqr; 
	chi2[0] = S2[0]/params->m2Sqr; 
	chi2[1] = S2[1]/params->m2Sqr; 
	chi2[2] = S2[2]/params->m2Sqr; 

    /* compute energy and flux function */
    dEbyF = XLALdEnergyByFluxSpinPrec(v, chi1, chi2, LNh, params);
        
    /* spin evolution equations - BBF Eq.(7.5 - 7.8) and Kesden et al 2010 Eq (2.2) */
    v2 = v*v;
	om0 = v2*v2*v/m;
    S1dotLNh = dotProduct(S1, LNh);
    S2dotLNh = dotProduct(S2, LNh);
    
    /* if low PN orders are used, set the higher order PN terms to zero */
    if (phaseO < 6) onePFpnFlag = 0;   /* 1.5 PN correction to leading order spin (3PN)  */
    if (phaseO < 5) onePNFlag = 0;     /* next to leading order in spin (2.5PN) */
    if (phaseO < 3) om0 = 0;           /* leading order spin (1.5PN) */

    for (i=0; i<3; i++) {
        Omega1[i] = om0*((0.75 + eta/2. - 0.75*delta) * LNh[i] 
             + onePNFlag * v *(-3.*(S2dotLNh + S1dotLNh*q) * LNh[i] + S2[i])/(2.*params->mSqr) 
             + onePFpnFlag * v2 *(0.5625 + 1.25*eta - params->etaSqr/24. - 0.5625*delta 
                 + 0.625*params->deltaEta)*LNh[i]);
         
        Omega2[i] = om0*((0.75 + eta/2. + 0.75*delta) * LNh[i] 
             + onePNFlag * v *(-3.*(S1dotLNh + S2dotLNh/q) * LNh[i] + S1[i])/(2.*params->mSqr) 
             + onePFpnFlag * v2 *(0.5625 + 1.25*eta - params->etaSqr/24. + 0.5625*delta 
                 - 0.625*params->deltaEta)*LNh[i]);
    }

    crossProduct(Omega1, S1, dS1byDt);
    crossProduct(Omega2, S2, dS2byDt);
    
    /* angular momentum evolution eqn (BCV2, Eq (9) */
    Lmag = (eta*params->mSqr/v)*(1. + (1.5+eta/6.)*v2);
	dLNhByDt[0] = -(dS1byDt[0] + dS2byDt[0])/Lmag; 
	dLNhByDt[1] = -(dS1byDt[1] + dS2byDt[1])/Lmag; 
	dLNhByDt[2] = -(dS1byDt[2] + dS2byDt[2])/Lmag; 

    /* evolve the variables  */
    dydt[0]  = -1./(dEbyF*m);
    dydt[1]  = dLNhByDt[0];
    dydt[2]  = dLNhByDt[1];
    dydt[3]  = dLNhByDt[2];
    dydt[4]  = dS1byDt[0];
    dydt[5]  = dS1byDt[1];
    dydt[6]  = dS1byDt[2];
    dydt[7]  = dS2byDt[0];
    dydt[8]  = dS2byDt[1];
    dydt[9]  = dS2byDt[2];
	
    if (dEbyF > 0.0) 								/* energy test fails! (dE/F > 0) */
        return LALSIMINSPIRAL_ST5_TEST_ENERGY;
	else 
		return GSL_SUCCESS;

}


/**
 * Compute the re-expanded (dEnergy/dv)/Flux for generic spinning binaries
 */
static REAL8 XLALdEnergyByFluxSpinPrec(
		double v, 
		double *chi1, 
		double *chi2, 
        double *LNh, 
		SpinTaylorT5Params *params
	) {

    REAL8 chi_s[3], chi_a[3];
    REAL8 chisDotLNh, chiaDotLNh, chisDotChia, chisSqr, chiaSqr; 
	REAL8 dEbF0 = 0., dEbF1 = 0., dEbF2 = 0., dEbF3 = 0., dEbF4 = 0., dEbF5 = 0., dEbF6 = 0.; 
	REAL8 dEbF6L = 0., dEbF7 = 0.;
	REAL8 v2, v3, v4, v5, v6, v7, v9; 

	REAL8 eta = params->eta; 
	REAL8 delta = params->delta; 
	UINT4 phaseO = params->phaseO;

    /* symmetric and anti-symmetric combinations of spin vectors */
	chi_s[0] = 0.5*(chi1[0]+chi2[0]);
	chi_s[1] = 0.5*(chi1[1]+chi2[1]);
	chi_s[2] = 0.5*(chi1[2]+chi2[2]);
	chi_a[0] = 0.5*(chi1[0]-chi2[0]);
	chi_a[1] = 0.5*(chi1[1]-chi2[1]);
	chi_a[2] = 0.5*(chi1[2]-chi2[2]);
    
    /* dot products  */
    chisDotLNh  = dotProduct(chi_s, LNh);
    chiaDotLNh  = dotProduct(chi_a, LNh);
    chisDotChia = dotProduct(chi_s, chi_a);
    chisSqr     = dotProduct(chi_s, chi_s);
    chiaSqr     = dotProduct(chi_a, chi_a);

	switch (phaseO) {
		case -1: 	/* highest available PN order */
		case 8: 	/* pseudo 4PN */
		case 7: 	/* 3.5 PN */
			dEbF7 = params->dEbF7nonSpin; 
		case 6: 	/* 3 PN */
			dEbF6 = params->dEbF6nonSpin; 
			dEbF6L = params->dEbF6LogNonSpin; 
		case 5: 	/* 2.5 PN */
			dEbF5 = -((7729*LAL_PI)/672. - (13*LAL_PI*eta)/8. + delta*(-72.96676587301587 - (13*eta)/4.)*chiaDotLNh 
				+ delta*(-0.46875 + (15*eta)/32.)*chiaDotLNh*chiaDotLNh*chiaDotLNh 
				+ delta*(-0.28125 + (9*eta)/32.)*chiaDotLNh*chiaSqr 
				+ (-72.96676587301587 + (2453*eta)/36. + (17*eta*eta)/2.)*chisDotLNh + (-1.40625 
				+ (135*eta)/32.)*chiaDotLNh*chiaDotLNh*chisDotLNh + (-0.28125 + (27*eta)/32.)*chiaSqr*chisDotLNh
				+ delta*(-1.40625 + (45*eta)/32.)*chiaDotLNh*chisDotLNh*chisDotLNh + (-0.46875 
				+ (45*eta)/32.)*chisDotLNh*chisDotLNh*chisDotLNh + (-0.5625 + (27*eta)/16.)*chiaDotLNh*chisDotChia 
				+ delta*(-0.5625 + (9*eta)/16.)*chisDotLNh*chisDotChia + delta*(-0.28125
				+ (9*eta)/32.)*chiaDotLNh*chisSqr + (-0.28125 + (27*eta)/32.)*chisDotLNh*chisSqr);

		case 4: 	/* 2 PN */
			dEbF4 = 3.010315295099521 + (233*chisDotChia*delta)/48. - (719*chiaDotLNh*chisDotLNh*delta)/48. 
				+ chiaSqr*(2.4270833333333335 - 10*eta) + chisDotLNh*chisDotLNh*(-7.489583333333333 - eta/24.) 
				+ chisSqr*(2.4270833333333335 + (7*eta)/24.) + (5429*eta)/1008. + (617*eta*eta)/144. 
				+ chiaDotLNh*chiaDotLNh*(-7.489583333333333 + 30*eta); 
		case 3: 	/* 1.5 PN */
			dEbF3 = (113*chiaDotLNh*delta)/12. + chisDotLNh*(9.416666666666666 - (19*eta)/3.) - 4*LAL_PI; 
		case 2: 	/* 1 PN */
			dEbF2 = params->dEbF2nonSpin; 
		case 1: 	/* 0.5 PN */
			dEbF1 = params->dEbF1nonSpin;
		case 0: 	/* Newtonian */
			dEbF0 = params->dEbF0nonSpin; 
			break; 
		default: 
            XLALPrintError("XLAL Error - %s: Invalid phase. PN order %s\n", __func__, phaseO);
            XLAL_ERROR(XLAL_EINVAL);
		break;
    }

	v2 = v*v; v3 = v2*v; v4 = v2*v2; v5 = v4*v; v6 = v3*v3; v7 = v6*v; v9 = v7*v2; 

	/* compute the re-expanded dEnerby/Flux function */
	return (dEbF0/v9)*(1. + dEbF1*v + dEbF2*v2 + dEbF3*v3 + dEbF4*v4 + dEbF5*v5 
			+ (dEbF6 + dEbF6L*log(4.*v))*v6 + dEbF7*v7);

}

/**
 * Compute the GW polarizations in the radiation frame
 */
static int polarizationsInRadiationFrame(
		REAL8TimeSeries **hplus,  	/**< +-polarization waveform [returned] */
		REAL8TimeSeries **hcross, 	/**< x-polarization waveform [returned] */
		REAL8TimeSeries *V,       	/**< post-Newtonian parameter */
		REAL8TimeSeries *Phi,     	/**< orbital phase */
		REAL8TimeSeries *LNhxVec,  	/**< unit orbital ang. mom. x comp. */
		REAL8TimeSeries *LNhyVec,  	/**< unit orbital ang. mom. y comp. */
		REAL8TimeSeries *LNhzVec,  	/**< unit orbital ang. mom. z comp. */
		REAL8 m1, 					/**< mass1 (KG) */
		REAL8 m2, 					/**< mass2 (KG) */
		REAL8 r, 					/**< distance (meter) */
		REAL8 Theta					/**< angle between the initial total ang momentum and line of sight. */
	) {

    REAL8 deltaPhi, twoPhiS, alphaDot;
    REAL8 LNx, LNy, LNz, v2, LNx_p2, LNy_p2, LNz_p2, LN_xz;
    UINT4 i; 
    REAL8Vector *alphaVec=NULL, *iVec=NULL;
	REAL8 cos_TwoPhiS, sin_TwoPhiS;

    UINT4 dataLength = V->data->length; 
	REAL8 deltaT = V->deltaT; 
	REAL8 m = (m1+m2)*LAL_MTSUN_SI/LAL_MSUN_SI;
	REAL8 eta = m1*m2/(m1+m2)/(m1+m2);
	REAL8 amp = 4.*eta*m/(r/LAL_C_SI); 

    /* allocate memory for the arrays */
    alphaVec = XLALCreateREAL8Vector(dataLength);
    iVec = XLALCreateREAL8Vector(dataLength);
	memset(alphaVec->data, 0, alphaVec->length * sizeof( REAL8 ));
	memset(iVec->data, 0, iVec->length * sizeof( REAL8 ));

    /* angles describing the time-evolution of the orientation of the orbital*/
    /* angular momentum in the source frame */
    for (i=0; i<dataLength; i++) {
        if (V->data->data[i]) {
            alphaVec->data[i] = atan2(LNhyVec->data->data[i], LNhxVec->data->data[i]);  // called alpha in BCV2
            iVec->data[i] = acos(LNhzVec->data->data[i]);                         // called i in BCV2
        }
    }
    
    /* unwrap the angles */
    XLALREAL8VectorUnwrapAngle(alphaVec, alphaVec);
    XLALREAL8VectorUnwrapAngle(iVec, iVec);

    deltaPhi = 0.;

    /* Allocate polarization vectors and set to 0 */
    *hplus = XLALCreateREAL8TimeSeries( "H_PLUS", &V->epoch, 
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *hcross = XLALCreateREAL8TimeSeries( "H_CROSS", &V->epoch, 
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    if ( ! hplus || ! hcross )	XLAL_ERROR(XLAL_EFUNC);
    memset((*hplus)->data->data, 0, (*hplus)->data->length*sizeof(*(*hplus)->data->data));
    memset((*hcross)->data->data, 0, (*hcross)->data->length*sizeof(*(*hcross)->data->data));

	REAL8 cos_Theta = cos(Theta);
	REAL8 cos_TwoTheta = cos(2.*Theta);
	REAL8 sin_Theta = sin(Theta);
	REAL8 sin_TwoTheta = sin(2.*Theta);

    for (i=0; i<dataLength; i++) {
    
        if (V->data->data[i]) {

          //REAL8 alpha = alphaVec->data[i]; // Warning: set-but-not-used!
            LNx = LNhxVec->data->data[i];
            LNy = LNhyVec->data->data[i];
            LNz = LNhzVec->data->data[i];
            v2 = V->data->data[i]*V->data->data[i];

			LNx_p2 = LNx*LNx; 
			LNy_p2 = LNy*LNy; 
			LNz_p2 = LNz*LNz; 
			LN_xz = LNx*LNz; 
        
            /* compute d alpha/dt*/
            if (i == dataLength-1) alphaDot = (alphaVec->data[i] - alphaVec->data[i-1])/deltaT;
            else alphaDot = (alphaVec->data[i+1] - alphaVec->data[i])/deltaT;

            /* \int \cos(i) \alpha_dot deltaT  -- integrating Eq.(18) of BCV2 */
            deltaPhi += LNz*alphaDot*deltaT;   

            /* the phase measured in the detector */
            twoPhiS = 2.*(Phi->data->data[i] - deltaPhi); 

			cos_TwoPhiS = cos(twoPhiS);
			sin_TwoPhiS = sin(twoPhiS);
    
			/* compute the polarizations. For derivation, see the Mathematica notebook Precession_ProjectiontoRadiationFrame.nb*/
			/* (based on Sec IIC of BCV2)*/
            (*hplus)->data->data[i] = -(0.5*amp*v2*(cos_TwoPhiS*(-LNx_p2 + LNy_p2*LNz_p2 + ((LNy 
					+ LN_xz)*cos_Theta + sin_Theta 
					- LNz_p2*sin_Theta)*((LNy - LN_xz)*cos_Theta + (-1. + LNz_p2)*sin_Theta)) 
					+ LNy*sin_TwoPhiS*(LN_xz*(3. + cos_TwoTheta) - (-1. + LNz_p2)*sin_TwoTheta)))/(-1. + LNz_p2);
            
            (*hcross)->data->data[i] = -(amp*v2*(cos_Theta*(-(LNx*LNy*(1 + LNz_p2)*cos_TwoPhiS) 
					+ (-LNx_p2 + LNy_p2)*LNz*sin_TwoPhiS) 
					+ (-1. + LNz_p2)*(LNy*LNz*cos_TwoPhiS + LNx*sin_TwoPhiS)*sin_Theta))/(-1. + LNz_p2);

        } 
			
    }

    XLALDestroyREAL8Vector(alphaVec);
    XLALDestroyREAL8Vector(iVec);

	return XLAL_SUCCESS; 
}

/**
 * Compute the orbital phase as an explicit function of v (TaylorT2 approximant)
 */
static int computeOrbitalPhase(
        SpinTaylorT5Params *params, 
		REAL8TimeSeries *V, 
		REAL8TimeSeries *LNhxVec, 
		REAL8TimeSeries *LNhyVec, 
		REAL8TimeSeries *LNhzVec, 
		REAL8TimeSeries *S1xVec, 
		REAL8TimeSeries *S1yVec, 
		REAL8TimeSeries *S1zVec, 
		REAL8TimeSeries *S2xVec, 
		REAL8TimeSeries *S2yVec, 
		REAL8TimeSeries *S2zVec, 
		REAL8TimeSeries *orbPhase
	) {

    REAL8 chi_s[3], chi_a[3], LNh[3];
    REAL8 chisDotLNh, chiaDotLNh, chisDotChia, chisSqr, chiaSqr; 
	REAL8 v, v2, v3, v4, v5, v6, v7;
	REAL8 phi0 = 0., phi1 = 0., phi2 = 0., phi3 = 0., phi4 = 0., phi5 = 0.;
	REAL8 phi6 = 0., phi6L = 0., phi7 = 0.;
	UINT4 i; 

	REAL8 eta = params->eta; 
	REAL8 delta = params->delta; 
	REAL8 phiOrb0 = params->phiRef/2. - orbPhase->data->data[0]; 
	UINT4 phaseO = params->phaseO;


	for (i = 0; i<V->data->length; i++) {

		if (V->data->data[i]) {

		   /* symmetric and anti-symmetric combinations of spin vectors */
			chi_s[0] = 0.5*(S1xVec->data->data[i]/params->m1Sqr + S2xVec->data->data[i]/params->m2Sqr); 
			chi_s[1] = 0.5*(S1yVec->data->data[i]/params->m1Sqr + S2yVec->data->data[i]/params->m2Sqr); 
			chi_s[2] = 0.5*(S1zVec->data->data[i]/params->m1Sqr + S2zVec->data->data[i]/params->m2Sqr); 
			chi_a[0] = 0.5*(S1xVec->data->data[i]/params->m1Sqr - S2xVec->data->data[i]/params->m2Sqr); 
			chi_a[1] = 0.5*(S1yVec->data->data[i]/params->m1Sqr - S2yVec->data->data[i]/params->m2Sqr); 
			chi_a[2] = 0.5*(S1zVec->data->data[i]/params->m1Sqr - S2zVec->data->data[i]/params->m2Sqr); 

			LNh[0] = LNhxVec->data->data[i];
			LNh[1] = LNhyVec->data->data[i];
			LNh[2] = LNhzVec->data->data[i];

			/* dot products */
			chisDotLNh  = dotProduct(chi_s, LNh);
			chiaDotLNh  = dotProduct(chi_a, LNh);
			chisDotChia = dotProduct(chi_s, chi_a);
			chisSqr     = dotProduct(chi_s, chi_s);
			chiaSqr     = dotProduct(chi_a, chi_a);
			
			v = V->data->data[i];

			switch (phaseO) {
				case -1: 	/* highest available PN order */
				case 8: 	/* pseudo 4PN */
				case 7: 	/* 3.5 PN */
					phi7 = params->phiOrb7nonSpin; 
				case 6: 	/* 3 PN */
					phi6 = params->phiOrb6nonSpin; 
					phi6L = params->phiOrb6LognonSpin; 
				case 5: 	/* 2.5 PN */
					phi5 = (chiaDotLNh*delta*(-363.58382936507934 - (35*eta)/2.) + chisDotLNh*(-363.58382936507934 + 
						(6065*eta)/18. + (85*eta*eta)/2.) + (38645*LAL_PI)/672. - (65*eta*LAL_PI)/8.)*log(v);
				case 4: 	/* 2 PN */
					phi4 = 15.051576475497606 + (1165*chisDotChia*delta)/48. - (3595*chiaDotLNh*chisDotLNh*delta)/48. + 
						chiaSqr*(12.135416666666666 - 50*eta) + 
						chisDotLNh*chisDotLNh*(-37.447916666666664 - (5*eta)/24.) + (27145*eta)/1008. + (3085*eta*eta)/144. + 
						chisSqr*(12.135416666666666 + (35*eta)/24.) + chiaDotLNh*chiaDotLNh*(-37.447916666666664 + 150*eta); 
				case 3: 	/* 1.5 PN */
					phi3 = (565*chiaDotLNh*delta)/24. + chisDotLNh*(23.541666666666668 - (95*eta)/6.) - 10*LAL_PI; 
				case 2: 	/* 1 PN */
					phi2 = params->phiOrb2nonSpin; 
				case 1: 	/* 0.5 PN */
					phi1 = params->phiOrb1nonSpin;
				case 0: 	/* Newtonian */
					phi0 = params->phiOrb0nonSpin; 
					break; 
				default: 
					XLALPrintError("XLAL Error - %s: Invalid phase. PN order %s\n", __func__, phaseO);
					XLAL_ERROR(XLAL_EINVAL);
				break;
			}

			v2 = v*v; v3 = v2*v; v4 = v2*v2; v5 = v4*v; v6 = v3*v3; v7 = v6*v;

			/* compute the re-expanded dEnerby/Flux function */
			orbPhase->data->data[i] = (phi0/v5)*(1. + phi1*v + phi2*v2 + phi3*v3 + phi4*v4 + phi5*v5 
					+ (phi6 + phi6L*log(4.*v))*v6 + phi7*v7)  + phiOrb0;

		}
	}
    
	return XLAL_SUCCESS; 
}

/**
 * Initialize the non-dynamical variables required for the SpinTaylorT5 evolution
 */
static int spinTaylorT5Init(
		REAL8 m1, 
		REAL8 m2, 
		REAL8 fStart, 
		REAL8 deltaT, 
		REAL8 phiRef, 
		UINT4 phaseO, 
		SpinTaylorT5Params *mParams
	) {

	/* various parameters describing the binary */
    mParams->mass1 = m1*LAL_MTSUN_SI/LAL_MSUN_SI;								/* m1 in seconds */ 
    mParams->mass2 = m2*LAL_MTSUN_SI/LAL_MSUN_SI;								/* m2 in seconds */
    mParams->totalMass = mParams->mass1+mParams->mass2;							/* m1+m2 in seconds */
    mParams->eta = mParams->mass1*mParams->mass2/(mParams->totalMass*mParams->totalMass);
    mParams->delta = (mParams->mass1-mParams->mass2)/mParams->totalMass;		/* (m1-m2)/(m1+m2) */
    mParams->phaseO = phaseO; 													/* twice the PN phase order */
    mParams->massRatio = mParams->mass2/mParams->mass1;							/* q = m2/m1, where m2 > m1 */
	mParams->m1Sqr = mParams->mass1*mParams->mass1; 							/* m1^2 in seconds^2 */
	mParams->m2Sqr = mParams->mass2*mParams->mass2; 							/* m2^2 in seconds^2 */
    mParams->mSqr = mParams->totalMass*mParams->totalMass;						/* m^2 in seconds^2 */
    mParams->etaSqr = mParams->eta*mParams->eta; 								/* eta^2 */
    mParams->deltaEta = mParams->delta*mParams->eta; 							/* delta * eta */
	mParams->v0 = cbrt(LAL_PI*mParams->totalMass*fStart);						/* starting value for the PN parameter */
    mParams->vMax = cbrt(0.45*LAL_PI*mParams->totalMass/deltaT);				/* set an emperical maximum on the PN parameter (0.9 f_Nyquist) */
	mParams->phiRef = phiRef; 

	/* coefficients of the reexpanded dEnergy/flux function */
	mParams->dEbF0nonSpin = -5./(32.*mParams->eta);
	mParams->dEbF1nonSpin = 0.;
	mParams->dEbF2nonSpin = 2.2113095238095237 + (11*mParams->eta)/4.; 
	mParams->dEbF3nonSpin = 0.;					/* currently spinning and non-spinning terms are not separated FIXME*/
	mParams->dEbF4nonSpin = 0.;					/* currently spinning and non-spinning terms are not separated FIXME*/
	mParams->dEbF5nonSpin = 0.;					/* currently spinning and non-spinning terms are not separated FIXME*/
	mParams->dEbF6nonSpin = -115.2253249962622 - (15211*mParams->etaSqr)/6912. 
			+ (25565*mParams->etaSqr*mParams->eta)/5184. 
			+ (32*LAL_PISQR)/3. + mParams->eta*(258.1491854023631 
			- (451*LAL_PISQR)/48.) + (1712*LAL_GAMMA)/105.; 
	mParams->dEbF6LogNonSpin = -1712./105.; 	/* coefficient of log(4v) at 3PN */
	mParams->dEbF7nonSpin = (-15419335*LAL_PI)/1.016064e6 
			- (75703*mParams->eta*LAL_PI)/6048. + (14809*mParams->etaSqr*LAL_PI)/3024.; 

	mParams->phiOrb0nonSpin = -1./(32.*mParams->eta);
    mParams->phiOrb1nonSpin = 0.;
    mParams->phiOrb2nonSpin = 3.685515873015873 + (55*mParams->eta)/12.; 
	mParams->phiOrb3nonSpin = 0.;
	mParams->phiOrb4nonSpin = 0.;
	mParams->phiOrb5nonSpin = 0.;
	mParams->phiOrb6nonSpin = 657.6504345051205 - (15737765635*mParams->eta)/1.2192768e7 + (76055*mParams->etaSqr)/6912. - 
                (127825*mParams->etaSqr*mParams->eta)/5184. - (160*LAL_PISQR)/3. + 
                (2255*mParams->eta*LAL_PISQR)/48. - (1712*LAL_GAMMA)/21.;
	mParams->phiOrb6LognonSpin = -1712/21.;
    mParams->phiOrb7nonSpin = (77096675*LAL_PI)/2.032128e6 + (378515*mParams->eta*LAL_PI)/12096. - (74045*mParams->etaSqr*LAL_PI)/6048.; 

	return XLAL_SUCCESS; 
}

/**
 * Internal function called by the integration routine.
 * Stops the integration if
 * 1) The energy decreases with increasing orbital frequency
 * 3) The orbital frequency becomes infinite
 */
static int XLALSimInspiralSpinTaylorT5StoppingTest(
		double t, 
		const double y[],
		double dy[], 
		void *mParams
	) {

    SpinTaylorT5Params *params = (SpinTaylorT5Params*) mParams;
    UNUSED(t);
	REAL8 LNh[3], S1[3], S2[3], chi1[3], chi2[3], v, dEbyF;

    /* evolution variables. */
    v      = y[0];   /* PN expansion parameter: (m \omega_orb)^(1/3) */
    LNh[0] = y[1];   /* unit vector in the direction of Newt. ang. momentum (LNh_x) */
    LNh[1] = y[2];   /* unit vector in the direction of Newt. ang. momentum (LNh_y) */
    LNh[2] = y[3];   /* unit vector in the direction of Newt. ang. momentum (LNh_z) */
    S1[0]  = y[4];   /* spin vector of the first object */
    S1[1]  = y[5];   /* spin vector of the first object */
    S1[2]  = y[6];   /* spin vector of the first object */
    S2[0]  = y[7];   /* spin vector of the second object  */
    S2[1]  = y[8];   /* spin vector of the second object  */
    S2[2]  = y[9];   /* spin vector of the second object  */
    
	/* comptute the dimensionless spins */
	chi1[0] = S1[0]/params->m1Sqr; 
	chi1[1] = S1[1]/params->m1Sqr; 
	chi1[2] = S1[2]/params->m1Sqr; 
	chi2[0] = S2[0]/params->m2Sqr; 
	chi2[1] = S2[1]/params->m2Sqr; 
	chi2[2] = S2[2]/params->m2Sqr; 

    /* compute energy and flux function */
    dEbyF = XLALdEnergyByFluxSpinPrec(v, chi1, chi2, LNh, params);

    if (dEbyF > -1.e-10) 							/* energy test fails! (dE/F > 0) */
        return LALSIMINSPIRAL_ST5_TEST_ENERGY;
    else if (dy[0] < 0.0) 							/* dv/dt < 0! */
        return LALSIMINSPIRAL_ST5_TEST_VDOT;		
    else if (isnan(v) || isinf(v))		 			/* v is nan! */
        return LALSIMINSPIRAL_ST5_TEST_VNAN;
    else if ((v < params->v0) || (v > params->vMax) || (v < 0.) || (v >= 1.))
        return LALSIMINSPIRAL_ST5_TEST_FREQBOUND;
    else 											/* Step successful, continue integrating */
        return GSL_SUCCESS;
}

static double dotProduct(double a[3], double b[3]) {

	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void crossProduct(double a[3], double b[3], double c[3]) {

    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static void rotateVector(double a[3], double b[3], double phi, double theta, double psi) {
/* rotate vector "a" by three Euler angles phi, theta, psi around axes 
 * z, y and x, respectively. b = R_x(psi) R_y(theta) R_z(phi) a */

    double x = a[0]; 
    double y = a[1];
    double z = a[2];

    b[0] = x*cos(phi)*cos(theta) - y*cos(theta)*sin(phi) 
            + z*sin(theta);
    
    b[1] = -(z*cos(theta)*sin(psi)) + x*(cos(psi)*sin(phi) 
            + cos(phi)*sin(psi)*sin(theta)) 
            + y*(cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta));
    
    b[2] = z*cos(psi)*cos(theta) + x*(sin(phi)*sin(psi) 
            - cos(phi)*cos(psi)*sin(theta)) 
            + y*(cos(phi)*sin(psi) + cos(psi)*sin(phi)*sin(theta));

}


