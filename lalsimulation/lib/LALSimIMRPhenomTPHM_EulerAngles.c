/*
 * Copyright (C) 2020 Hector Estelles
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


/**
 * \author Hector Estelles
 */

#include <math.h>

/* LAL Header Files */
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>

/* GSL Header Files */
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>

#include <lal/XLALGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

#include <lal/LALSimInspiralPrecess.h>

#include "LALSimIMRPhenomTHM_internals.h"

#define MIN(a,b) (((a)<(b))?(a):(b))

// FIXME: reference equations to PhenomTP paper and future TPHM paper.

/* Routines and wrappers to compute precessing Euler angles for 'twisting-up' IMRPhenomTHM into IMRPhenomTPHM */

/* ************************************************************* */
/* ************** STRUCTS and AUXILIARY FUNCTIONS ************** */
/* ************************************************************* */

/* Collection of functions and structs needed for the Euler angle routines */

// Tolerance for integrator in numerical angles routine
#define LAL_ST4_ABSOLUTE_TOLERANCE 1.e-12
#define LAL_ST4_RELATIVE_TOLERANCE 1.e-12

/* MECO time fit: Needed for option EulerRDVersion=2. */
static double IMRPhenomT_MECOTime(double eta, double S, double dchi, double delta);
static double IMRPhenomT_MECOTime(double eta, double S, double dchi, double delta){

	double nospin, eqspin, uneqspin;

   	nospin = (-236.05916243036953 - 831.5431113584278*eta)/(1 + 20.613030474069898*eta);

   	eqspin = (48.728774732041 + 8.161848266333303*eta + 146.22536542867712*eta*eta)*S + (-13.633462771084439 + 97.744032521123*eta - 333.74293640880865*eta*eta)*S*S + 
   (-76.67450724471107 + 1030.4869060625915*eta - 3160.810490374918*eta*eta)*S*S*S + (-20.941621754455277 - 26.348664719217858*eta + 778.3469028250508*eta*eta)*S*S*S*S + 
   (86.59936006891306 - 1115.1092303071634*eta + 3348.2374777841*eta*eta)*S*S*S*S*S + (13.413110005766034 + 46.573063129481895*eta - 631.7603046833175*eta*eta)*S*S*S*S*S*S;

   	uneqspin = -1385.6109494038105*dchi*delta*(1 - 2.465600020183968*eta)*eta*eta + 10.837064426546098*dchi*dchi*eta*eta*eta + 355.1229694251773*dchi*delta*(1 - 4.1057183984004695*eta)*eta*eta*S;

   	return (nospin + eqspin + uneqspin);
}

/* Struct needed for computing gamma from the minimal rotation condition. */
typedef struct taggammaIntegration
{
   gsl_spline *alpha_spline;
   gsl_spline *cosbeta_spline;
   gsl_interp_accel *alpha_acc;
   gsl_interp_accel *cosbeta_acc;
}
gammaIntegration;

/* Function to return RHS of minimal rotation condition equation, for computing gamma angle. */
static double f_alphadotcosi( double x, void * inparams )
{
  gammaIntegration* params = (gammaIntegration*) inparams;
  REAL8 alphadot = gsl_spline_eval_deriv( params->alpha_spline, x, params->alpha_acc );
  REAL8 cosbeta = gsl_spline_eval( params->cosbeta_spline, x, params->cosbeta_acc );

  return -1. * alphadot * cosbeta;
}


/* Functions to perform a ZYZ rotation of a given vector by some given angles. */
void IMRPhenomT_rotate_z(REAL8 cosangle, REAL8 sinangle, REAL8 *vx, REAL8 *vy, REAL8 *vz);
void IMRPhenomT_rotate_y(REAL8 cosangle, REAL8 sinangle, REAL8 *vx, REAL8 *vy, REAL8 *vz);

void IMRPhenomT_rotate_z(REAL8 cosangle, REAL8 sinangle, REAL8 *vx, REAL8 *vy, REAL8 *vz)
{
  const REAL8 tmpx = *vx;
  const REAL8 tmpy = *vy;
  const REAL8 tmpz = *vz;

  REAL8 tmp1 = tmpx*cosangle - tmpy*sinangle;
  REAL8 tmp2 = tmpx*sinangle + tmpy*cosangle;

  *vx = tmp1;
  *vy = tmp2;
  *vz = tmpz;
}

void IMRPhenomT_rotate_y(REAL8 cosangle, REAL8 sinangle, REAL8 *vx, REAL8 *vy, REAL8 *vz)
{
  const REAL8 tmpx = *vx;
  const REAL8 tmpy = *vy;
  const REAL8 tmpz = *vz;

  REAL8 tmp1 = + tmpx*cosangle + tmpz*sinangle;
  REAL8 tmp2 = - tmpx*sinangle + tmpz*cosangle;

  *vx = tmp1;
  *vy = tmpy;
  *vz = tmp2;
}

/* Function to unwrap a time series that contains an angle, to obtain a continuous time series. */
void unwrap_array(double *in, double *out, int len);
void unwrap_array(double *in, double *out, int len) {
    out[0] = in[0];
    for (int i = 1; i < len; i++) {
        double d = in[i] - in[i-1];
        d = d > M_PI ? d - 2 * M_PI : (d < -M_PI ? d + 2 * M_PI : d);
        out[i] = out[i-1] + d;
    }
}

/* Struct needed for numerical integration of the spin precessing equations */
typedef struct tagPhenomTPHMEvolution
{
   gsl_spline *v_spline;
   gsl_interp_accel *v_acc;
   XLALSimInspiralSpinTaylorTxCoeffs* Tparams;
   REAL8 ToffSign;
}PhenomTPHMEvolution;

/*******************************************************/
/*************** ANALYTICAL EULER ANGLES ***************/
/*******************************************************/


/** This function provides a wrapper to the analytical PN angle descriptions computed by
the IMRPhenomXP(HM) waveform models and routines to incorporate a simple approximation in the ringdown
for the precessing angles (see https://arxiv.org/pdf/1209.3712.pdf, https://arxiv.org/abs/1806.10734, https://arxiv.org/abs/2004.08302). 
It provides also additional ways of treating the plunge-merger angles and possibility to compute third
Euler angle gamma directly from the minimal rotation condition. **/

/* Summary of options:
- EulerRDVersion: 0 (no RD angles attached, MSA/NNLO computed for the whole waveform)
                  1 (MSA/NNLO computed up to the peak time of the 22 mode, then RD angles are attached.)
                  2 (MSA/NNLO computed up to MECO time, linear continuation performed up to the peak time, then RD angles are attached.)

- GammaVersion:   0 (Gamma is computed from the MSA/NNLO PN expressions.)
                  1 (Gamma computed numerically from minimal rotation condition.)
*/  

int PNAnalyticalInspiralEulerAngles(
  REAL8TimeSeries **alphaTS,          /* Alpha angle time series [out] */ 
  REAL8TimeSeries **cosbetaTS,        /* cos(Beta) angle time series [out] */
  REAL8TimeSeries **gammaTS,          /* Gamma angle time series [out] */
  REAL8Sequence *xorb,                /* x(t)=v(t)^2 time series [in] */
  IMRPhenomTWaveformStruct *pWF,      /* IMRPhenomT* waveform struct */
  IMRPhenomTPhase22Struct *pPhase,    /* IMRPhenomT phase struct */
  IMRPhenomXWaveformStruct *pWFX,     /* IMRPhenomX waveform struct */
  IMRPhenomXPrecessionStruct *pPrec,  /* IMRPhenomX precessing struct */
  INT4 EulerRDVersion,                /* Version of ringdown angles attachment */
  INT4 GammaVersion);                 /* Version of gamma evaluation */

int PNAnalyticalInspiralEulerAngles(
  REAL8TimeSeries **alphaTS,
  REAL8TimeSeries **cosbetaTS,
  REAL8TimeSeries **gammaTS,
  REAL8Sequence *xorb,
  IMRPhenomTWaveformStruct *pWF,
  IMRPhenomTPhase22Struct *pPhase,
  IMRPhenomXWaveformStruct *pWFX,
  IMRPhenomXPrecessionStruct *pPrec,
  INT4 EulerRDVersion,
  INT4 GammaVersion)
{

  /***** Set up analytical PN evaluation length for the different versions ****/
  
  size_t length = 0; // Initialize length
  REAL8 tMECO = IMRPhenomT_MECOTime(pWF->eta, pWF->Shat, pWF->dchi, pWF->delta); //MECO time needed for EulerRDVersion=2

  /* If minimum time is greater that tMECO, EulerRDVersion=2 is not allowed. Defaulting to EulerRDVersion=1. */
  if(pPhase->tmin>=tMECO - 10*pPhase->dtM && EulerRDVersion == 2)
    {
      EulerRDVersion = 1;
      XLAL_PRINT_WARNING("Waveform is too short for EulerRDVersion=2. Defaulting to EulerRDVersion=1.\n");
    }
  if(pPhase->tmin>= - 10*pPhase->dtM && GammaVersion == 1)
    {
      GammaVersion = 0;
      XLAL_PRINT_WARNING("Waveform is too short for GammaVersion=1. Defaulting to GammaVersion=0.\n");
    }

  /* Check if spins are almost aligned. In this case, disable ringdown angles addition */
  REAL8 AStol = 1E-4;
  if((sqrt(pPrec->chi1x * pPrec->chi1x + pPrec->chi1y * pPrec->chi1y) < AStol) &&
      (sqrt(pPrec->chi2x * pPrec->chi2x + pPrec->chi2y * pPrec->chi2y) < AStol))
  {
    EulerRDVersion = 0;
    XLAL_PRINT_WARNING("System is almost aligned spin. Defaulting to EulerRDVersion=0.\n");
  }

  /* Different versions of the ringdown Euler angles attachment.
  Depending on the version, inspiral angles need to be computed only up to a certain length. */
  if(EulerRDVersion == 0)
  {
   	length = pPhase->wflength; // Whole waveform length
  }
  else if(EulerRDVersion == 1)
  {
    if(pPhase->tmin>=0) // If waveform contains only ringdown, set PN length to 0.
    {
      length = 0;
    }
    else{
      length = floor((0.0 - pPhase->tmin)/pWF->dtM); // If not, set PN length from tmin to t=0 (peak time)
    }
  }
  else if(EulerRDVersion == 2)
  {
    REAL8 t0 = tMECO;
    length = floor((t0 - pPhase->tmin)/pWF->dtM); // set PN length from tmin to tMECO
  }

  /* Offset for alpha angle. Essentially is the initial azimuthal angle of Lhat in the J frame */
  REAL8 alphaOff = atan2(pPrec->S1y + pPrec->S2y, pPrec->S1x + pPrec->S2x) - LAL_PI;

  REAL8 tt;

  /* Analytical Euler angles from PhenomX implementation */
  switch(pPrec->IMRPhenomXPrecVersion) 
  {
    /* ~~~~~ Use NNLO PN Euler Angles - Appendix G of arXiv:2004.06503 and https://dcc.ligo.org/LIGO-T1500602 ~~~~~ */
    case 101:
    case 102:
    case 103:
    case 104:
    {
      REAL8 s, s2, cos_beta, omega, omega_cbrt2, omega_cbrt, logomega, v, L;

      /* Compute values at reference frequency */

      omega_cbrt2 = pow(LAL_PI*pWF->MfRef,2./3);
      omega_cbrt = sqrt(omega_cbrt2);
      omega = omega_cbrt2*omega_cbrt;
      logomega = log(omega);

      REAL8 alphaRef         = IMRPhenomX_PN_Euler_alpha_NNLO(pPrec,omega,omega_cbrt2,omega_cbrt,logomega) ;
      REAL8 gammaRef       = -IMRPhenomX_PN_Euler_epsilon_NNLO(pPrec,omega,omega_cbrt2,omega_cbrt,logomega);

      /* Loop for evaluating NNLO angles in the desired points */
      for(UINT4 jdx = 0; jdx < length; jdx++)
      {
        omega_cbrt2 = xorb->data[jdx];
        omega_cbrt = sqrt(omega_cbrt2);
        v = omega_cbrt;
        omega = omega_cbrt2*omega_cbrt;
        logomega = log(omega);

        (*alphaTS)->data->data[jdx]         = IMRPhenomX_PN_Euler_alpha_NNLO(pPrec,omega,omega_cbrt2,omega_cbrt,logomega) - alphaRef + alphaOff ;
        (*gammaTS)->data->data[jdx]       = -IMRPhenomX_PN_Euler_epsilon_NNLO(pPrec,omega,omega_cbrt2,omega_cbrt,logomega) - gammaRef - alphaOff ;

        L = XLALSimIMRPhenomXLPNAnsatz(v, pWFX->eta/v, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);

        s        = pPrec->Sperp / (L + pPrec->SL);
        s2       = s*s;
        cos_beta = copysign(1.0, L + pPrec->SL) / sqrt(1.0 + s2);

        (*cosbetaTS)->data->data[jdx] = cos_beta;
      }

     break;
    }
    case 220:
    case 221:
    case 222:
    case 223:
    case 224:
    {
      vector vangles = {0.,0.,0.};
      REAL8 v, cos_beta;

      /* ~~~~~ Euler Angles from Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967 ~~~~~ */
      for(UINT4 jdx = 0; jdx < length; jdx++)
      {
        v = sqrt(xorb->data[jdx]);

        vangles  = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(v,pWFX,pPrec);

        (*alphaTS)->data->data[jdx]    = vangles.x + alphaOff;
        (*gammaTS)->data->data[jdx]  = -vangles.y - alphaOff;
        cos_beta = vangles.z;
        (*cosbetaTS)->data->data[jdx] = cos_beta;
      }

     break;
    }
    default:
    {
      XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPrecessionVersion not recognized. Recommended default is 223.\n");
      break;
    }
  }

  /**** COMPUTE GAMMA FROM MINIMAL ROTATION CONDITION (if specified) ****/
  /* alpha and cosbeta are interpolated to compute the RHS of the minimal rotation condition.
  Then gamma is integrated using Boole's rule. */
  if(GammaVersion == 1 && length>0)
  {
    //printf("length: %zu\n",length);
  	// Generate time array
  	REAL8Sequence *timesPN = NULL;
  	timesPN = XLALCreateREAL8Sequence(length);
  	for(UINT8 i=0; i < length; i++) 
  	{
    	timesPN->data[i] = i*pWF->dtM; 
  	}

  	// Interpolate alpha
  	gsl_interp_accel *accel_alpha;
  	gsl_spline *spline_alpha;

  	accel_alpha = gsl_interp_accel_alloc();
  	spline_alpha = gsl_spline_alloc(gsl_interp_cspline, length);
 	  gsl_spline_init(spline_alpha, timesPN->data, (*alphaTS)->data->data, length);

  	// Interpolate cosbeta
  	gsl_interp_accel *accel_cosb;
  	gsl_spline *spline_cosb;

  	accel_cosb = gsl_interp_accel_alloc();
  	spline_cosb = gsl_spline_alloc(gsl_interp_cspline, length);
  	gsl_spline_init(spline_cosb, timesPN->data, (*cosbetaTS)->data->data, length);

  	// Setup helper struct for gamma integration
  	gammaIntegration gammaStruct;
  	gammaStruct.alpha_spline = spline_alpha;
  	gammaStruct.alpha_acc    = accel_alpha;
  	gammaStruct.cosbeta_spline  = spline_cosb;
  	gammaStruct.cosbeta_acc     = accel_cosb;

  	// Integrate gamma using Boole's rule. (https://en.wikipedia.org/wiki/Boole%27s_rule)
    // Gamma is obtained through the minimal rotation condition: \dot{gamma} = -\dot{alpha(t)}*cos[beta(t)]
    for( UINT4 i = 1; i < timesPN->length; i++ ){
  	double t1 = timesPN->data[i-1];
  	double t2 = timesPN->data[i];
  	(*gammaTS)->data->data[i] = (*gammaTS)->data->data[i-1] + (1.0/90.0)*(t2 - t1)*(7.0*f_alphadotcosi(t1,&gammaStruct)
                  + 32.0*f_alphadotcosi((t1+3.0*t2)/4.0,&gammaStruct)
                  + 12.0*f_alphadotcosi((t1+t2)/2.0,&gammaStruct)
                  + 32.0*f_alphadotcosi((3.0*t1+t2)/4.0,&gammaStruct)
                  + 7.0*f_alphadotcosi(t2,&gammaStruct));/* Boole's Rule */

      }

    // Free gsl spline objects  
    gsl_spline_free(spline_alpha);
    gsl_spline_free(spline_cosb);
    gsl_interp_accel_free(accel_alpha);
    gsl_interp_accel_free(accel_cosb);

    XLALDestroyREAL8Sequence( timesPN);

  }

  /**** RINGDOWN APPROXIMATION FOR EULER ANGLES ****/
  /* (see https://arxiv.org/pdf/1209.3712.pdf, https://arxiv.org/abs/1806.10734, https://arxiv.org/abs/2004.08302)
  Analytical ringdown approximation for the Euler angles is discussed in the above papers.
  Here, that approximation is attached to the Euler angles time series at t=0, i.e the peak time of the non-precessing 22 mode.
  For EulerRDVersion=2, the angles inthe strong field regime between the MECO time (end of validity of PN framework) and the peak time
  are modelled also by a simple analytical linear continuation, to try to mitigative possible PN pathologies in this region. */

  if(EulerRDVersion == 1)
  {
    /* Offsets for the ringdown angles, in order to impose continuity with inspiral angles */
    
    // Offsets initialization
    REAL8 alphaRD0 = 0.0;
    REAL8 cosbetaRD0 = 0.0;
    REAL8 gammaRD0 = 0.0;

    /* If waveform too short to contain prepeak cycles, impose offsets */
    if(pPhase->tmin>=-pWF->dtM)
    {
      alphaRD0 = alphaOff;
      cosbetaRD0 = cos(pPrec->thetaJ_Sf);
      gammaRD0 = -alphaOff;
      length = 0;
    }
    /* If waveform long enough, set offsets for continuity with preRD region */
    else{
      alphaRD0 = (*alphaTS)->data->data[length-1];
      cosbetaRD0 = (*cosbetaTS)->data->data[length-1];
      gammaRD0 = (*gammaTS)->data->data[length-1];
    }

   /* Compute Euler angles in the RD approximation */
   for(UINT4 jdx = length; jdx < pPhase->wflength; jdx++)
    {
      	tt = (jdx - length)*pWF->dtM;
        (*alphaTS)->data->data[jdx] = alphaRD0 + pPhase->EulerRDslope*tt;
        (*cosbetaTS)->data->data[jdx] = cosbetaRD0;
        (*gammaTS)->data->data[jdx] = gammaRD0 - (*alphaTS)->data->data[jdx]*(*cosbetaTS)->data->data[jdx] + alphaRD0*cosbetaRD0;
    }

  }

  /*** Alternative way of treating Euler angles during plunge-merger ***/
  /* MSA/NNLO angles computed up tp the MECO time, and then linear continuation to the peak time. Then Euler RD angles are attached. */
  else if(EulerRDVersion == 2)
  {
    
    REAL8 t0 = tMECO;
  	size_t lengthInt = floor((t0 - pPhase->tmin)/pWF->dtM); // Length from tmin to tMECO
  	length = floor((0.0 - pPhase->tmin)/pWF->dtM); // Length from tmin to t=0

    /* beta values for computing beta derivative */
    REAL8 beta1 = acos((*cosbetaTS)->data->data[lengthInt-1]); 
    REAL8 beta2 = acos((*cosbetaTS)->data->data[lengthInt-2]);
    REAL8 beta3 = acos((*cosbetaTS)->data->data[lengthInt-3]);

    // Compute alpha and cosbeta derivatives at tMECO using three point backwards discretization.
    REAL8 alphader = (3.0*(*alphaTS)->data->data[lengthInt-1] - 4.0*(*alphaTS)->data->data[lengthInt-2] + (*alphaTS)->data->data[lengthInt-3])/(2*pWF->dtM);
    REAL8 betader = (3.0*beta1 - 4.0*beta2 + beta3)/(2*pWF->dtM);

    //printf("alphader: %.16f, betader: %.16f\n",alphader, betader);

    // Angle values at tMECO.
    REAL8 alphaInt0 = (*alphaTS)->data->data[lengthInt-1];
    REAL8 betaInt0 = beta1;
    REAL8 gammaInt0 = (*gammaTS)->data->data[lengthInt-1];


    // Compute alpha and cosbeta as a linear continuation from tMECO to peak time. Gamma analytically computed from minimal rotation condition.
    int kk = 0;
    REAL8 betalin = 0.0;

    if(fabs(betader)>1E-15)
    {
      for(UINT4 jdx = lengthInt; jdx < length; jdx++)
      {
          (*alphaTS)->data->data[jdx]    = alphaInt0 + alphader*kk*pWF->dtM;
          betalin = betaInt0 + betader*kk*pWF->dtM;
          if(betalin > LAL_PI){betalin=LAL_PI; XLAL_PRINT_INFO("Warning: beta>180º reached, correcting to 180º.");}
          if(betalin < 0.0){betalin=0.0;       XLAL_PRINT_INFO("Warning: beta<0º reached, correcting to 0º.");}
          (*cosbetaTS)->data->data[jdx]    = cos(betalin);
          (*gammaTS)->data->data[jdx]  = gammaInt0  - (alphader/betader)*sin(betaInt0 + betader*kk*pWF->dtM) + (alphader/betader)*sin(betaInt0);
          kk++;
      }
    }
    else{
      for(UINT4 jdx = lengthInt; jdx < length; jdx++)
      {
          tt = (jdx - lengthInt)*pWF->dtM;
          (*alphaTS)->data->data[jdx]    = alphaInt0 + alphader*kk*pWF->dtM;
          betalin = betaInt0 + betader*kk*pWF->dtM;
          if(betalin > LAL_PI){betalin=LAL_PI; XLAL_PRINT_INFO("Warning: beta>pi reached, correcting to pi.");}
          if(betalin < 0.0){betalin=0.0;       XLAL_PRINT_INFO("Warning: beta<0 reached, correcting to 0.");}
          (*cosbetaTS)->data->data[jdx]    = cos(betalin);
          (*gammaTS)->data->data[jdx]  = gammaInt0  - alphader*tt*cos(betaInt0);
          kk++;
      }
    }

    // At peak time, attach Euler angles computed from the RD approximation.
    for(UINT4 jdx = length; jdx < pPhase->wflength; jdx++)
    {
      	tt = (jdx - length)*pWF->dtM;
        (*alphaTS)->data->data[jdx] = (*alphaTS)->data->data[length-1] + pPhase->EulerRDslope*tt;
        (*cosbetaTS)->data->data[jdx] = (*cosbetaTS)->data->data[length-1];
        (*gammaTS)->data->data[jdx] = (*gammaTS)->data->data[length-1] - (*alphaTS)->data->data[jdx]*(*cosbetaTS)->data->data[jdx] + (*alphaTS)->data->data[length-1]*(*cosbetaTS)->data->data[length-1];
    }

  }

  return XLAL_SUCCESS;

}

/*******************************************************/
/*************** NUMERICAL EULER ANGLES ****************/
/*******************************************************/

/* This section provides the needed functions to perform the evolution of the spin precessing equations, that are specified in XLALSimInspiralSpinDerivatives and 
related SimInspiral code employing the PN expansion parameter v(t) computed from the orbital frequency predicted by IMRPhenomT, instead of evolving v from 
a TaylorT* approximant. */ 

/* Function needed by RungeKutta integrator, we just provide a trivial GSL_SUCCESS since TaylorT* stopping conditions do not apply here. */
int StoppingTest(double UNUSED t, const double UNUSED values[], double UNUSED dvalues[], void UNUSED *mparams);

int StoppingTest(double UNUSED t, const double UNUSED values[], double UNUSED dvalues[], void UNUSED *mparams)
{
  return GSL_SUCCESS;
}

/* Function to provide RHS of the spin precessing equations to the integrator, analogous to XLALSimInspiralSpinTaylorT{1/4/5}Derivatives function.
Main difference is that here v(t) is read from an interpolated spline, instead of evolved from a vdot expression. */
INT4 XLALSimIMRPhenomTPHMSpinDerivatives(
  REAL8 UNUSED t,
  const REAL8 values[],
  REAL8 dvalues[],
  void *mparams
  );

INT4 XLALSimIMRPhenomTPHMSpinDerivatives(
  REAL8 UNUSED t,
  const REAL8 values[],
  REAL8 dvalues[],
  void *mparams
  )
{
    INT4 status;

    /* coordinates and derivatives */
    REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z, E1x, E1y, E1z;
    REAL8 dLNhx, dLNhy, dLNhz, dS1x, dS1y, dS1z, dS2x, dS2y, dS2z, dE1x, dE1y, dE1z;

    /* auxiliary variables */
    REAL8 v;
    REAL8 LNhdotS1, LNhdotS2;

    /* Structs */
    PhenomTPHMEvolution *params = mparams;
    XLALSimInspiralSpinTaylorTxCoeffs *Tparams = params->Tparams;

    /* copy variables */
    LNhx = values[0] ; LNhy     = values[1] ; LNhz  = values[2] ;
    S1x  = values[3] ; S1y      = values[4] ; S1z   = values[5] ;
    S2x  = values[6] ; S2y      = values[7] ; S2z   = values[8];
    E1x  = values[9]; E1y     = values[10]; E1z   = values[11];

    /* Variables for relating integration time t with time employed in the interpolation of v(t) */
    REAL8 toff = fabs(params->ToffSign);
    REAL8 sign = copysign(1.0, params->ToffSign);

    /* Evaluate v(t) from an interpolation of IMRPhenomT v(t) */
    v = gsl_spline_eval( params->v_spline, toff + sign*t, params->v_acc );

    /* Compute derivatives of quantities that are evolved */
    LNhdotS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
    LNhdotS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);

    status = XLALSimInspiralSpinDerivatives(&dLNhx,&dLNhy,&dLNhz,&dE1x,&dE1y,&dE1z,&dS1x,&dS1y,&dS1z,&dS2x,&dS2y,&dS2z,v,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,LNhdotS1,LNhdotS2,Tparams);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function XLALSimInspiralSpinDerivatives has failed.");

    dvalues[0]    = sign*dLNhx; dvalues[1]     = sign*dLNhy ; dvalues[2]    = sign*dLNhz;
    dvalues[3]    = sign*dS1x ; dvalues[4]     = sign*dS1y  ; dvalues[5]    = sign*dS1z ;
    dvalues[6]    = sign*dS2x ; dvalues[7]     = sign*dS2y  ; dvalues[8]   = sign*dS2z ;
    dvalues[9]   = sign*dE1x ; dvalues[10]    = sign*dE1y  ; dvalues[11]   = sign*dE1z ;

    return GSL_SUCCESS;
}


/* Function to evolve LNhat and individual spins using IMRPhenomT phase evolution.
 Similar to XLALSimInspiralSpinTaylorPNEvolveOrbit. */

int IMRPhenomTPHM_EvolveOrbit(
  REAL8TimeSeries **V,            /**< post-Newtonian parameter [returned]*/
  REAL8TimeSeries **S1x,          /**< Spin1 vector x component [returned]*/
  REAL8TimeSeries **S1y,          /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **S1z,          /**< "    "    "  z component [returned]*/
  REAL8TimeSeries **S2x,          /**< Spin2 vector x component [returned]*/
  REAL8TimeSeries **S2y,          /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **S2z,          /**< "    "    "  z component [returned]*/
  REAL8TimeSeries **LNhatx,       /**< unit orbital ang. mom. x [returned]*/
  REAL8TimeSeries **LNhaty,       /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **LNhatz,       /**< "    "    "  z component [returned]*/
  REAL8TimeSeries **E1x,          /**< orb. plane basis vector x[returned]*/
  REAL8TimeSeries **E1y,          /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **E1z,          /**< "    "    "  z component [returned]*/
  REAL8Sequence *xorb,            /**< squared velocity from IMRPhenomT phase evolution*/
  REAL8 m1_SI,                    /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                    /**< Mass of companion 2 (kg) */
  REAL8 s1x,                      /**< x component of primary spin at reference frequency*/
  REAL8 s1y,                      /**< y component of primary spin at reference frequency*/
  REAL8 s1z,                      /**< z component of primary spin at reference frequency*/
  REAL8 s2x,                      /**< x component of secondary spin at reference frequency*/
  REAL8 s2y,                      /**< y component of secondary spin at reference frequency*/
  REAL8 s2z,                      /**< z component of secondary spin at reference frequency*/
  IMRPhenomTWaveformStruct *pWF,  /**< PhenomTHM waveform struct*/
  IMRPhenomTPhase22Struct *pPhase /**< PhenomT phase struct*/
  );

/**
 * @addtogroup LALSimIMRPhenomT_c
 * @{
 *
 * @name Routines for spin evolution
 * @{
 *
 * @author Héctor Estellés
 *
 * @brief C Code for the spin evolution of the IMRPhenomTP(HM) models.
 */
/**
 * Function to return the time series for the evolution of the individual spins, the Newtonian angular momentum direction and the orbital plane basis vector E.
 * It computes the frequency of the non-precessing 22 mode between minimum frequency and the peak time of the 22,
 * and approximates the PN parameter v(t) from it. Thenn it calls the internal function IMRPhenomTPHM_EvolveOrbit which evolves the spin precessing PN evolution equations
 * using the precomputed v(t) as a driver for the evolution. 
 */

int XLALSimIMRPhenomTPHM_EvolveOrbit( //FIXME: maybe move to Euler angles? Maria will take a look
  REAL8TimeSeries **V,            /**< post-Newtonian parameter [returned]*/
  REAL8TimeSeries **S1x,          /**< Spin1 vector x component [returned]*/
  REAL8TimeSeries **S1y,          /**< Spin1 vector y component [returned]*/
  REAL8TimeSeries **S1z,          /**< Spin1 vector z component [returned]*/
  REAL8TimeSeries **S2x,          /**< Spin2 vector x component [returned]*/
  REAL8TimeSeries **S2y,          /**< Spin2 vector y component [returned]*/
  REAL8TimeSeries **S2z,          /**< Spin2 vector z component [returned]*/
  REAL8TimeSeries **LNhatx,       /**< unit orbital ang. mom. x [returned]*/
  REAL8TimeSeries **LNhaty,       /**< unit orbital ang. mom. y [returned]*/
  REAL8TimeSeries **LNhatz,       /**< unit orbital ang. mom. z [returned]*/
  REAL8TimeSeries **E1x,          /**< orb. plane basis vector x [returned]*/
  REAL8TimeSeries **E1y,          /**< orb. plane basis vector y [returned]*/
  REAL8TimeSeries **E1z,          /**< orb. plane basis vector z [returned]*/
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams       /**< LAL dictionary containing accessory parameters */
  )
{

  int status;

  /* Sanity checks */
  if(fRef  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");  }
  if(deltaT   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaT must be positive.\n");                         }
  if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                             }
  if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                             }
  if(fmin    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "f_min must be positive.\n");                          }
  if(fRef > 0.0 && fRef < fmin){ XLAL_ERROR(XLAL_EDOM, "fRef must be >= f_min or =0 to use f_min.\n"); }

  /* Swap components if m2>m1 */  
  if(m1_SI < m2_SI)
  {
    REAL8 auxswap;
    
    auxswap = m2_SI;
    m2_SI = m1_SI;
    m1_SI = auxswap;

    auxswap = chi2x;
    chi2x = chi1x;
    chi1x = auxswap;

    auxswap = chi2y;
    chi2y = chi1y;
    chi1y = auxswap;

    auxswap = chi2z;
    chi2z = chi1z;
    chi1z = auxswap;
  }
  
  REAL8 mass_ratio = m1_SI / m2_SI;
  REAL8 chi1L = chi1z;
  REAL8 chi2L = chi2z;

  /*
    Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
      - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
      - For 200 > mass ratio > 20 and spins <= 0.9: print a warning message that model can be pathological.
      - For mass ratios > 200: throw a hard error that model is not valid.
      - For spins > 0.99: throw a warning that we are extrapolating the model to extremal
  */

  if(mass_ratio > 20.0 && chi1L < 0.9 && m2_SI/LAL_MSUN_SI >= 0.5  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
  if(mass_ratio > 20.0 && (chi1L >= 0.9 || m2_SI/LAL_MSUN_SI < 0.5) ) { XLAL_PRINT_INFO("Warning: Model can be pathological at these parameters"); }
  if(mass_ratio > 200. && fabs(mass_ratio - 200) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 200."); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }


  IMRPhenomTWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomTWaveformStruct));
  status = IMRPhenomTSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, 1.0, deltaT, fmin, fRef, phiRef, lalParams);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetWaveformVariables has failed.");

  IMRPhenomTPhase22Struct *pPhase;
  pPhase = XLALMalloc(sizeof(IMRPhenomTPhase22Struct));
  status   = IMRPhenomTSetPhase22Coefficients(pPhase, pWF);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetPhase22Coefficients has failed.");

  /* Length of the required sequence and length of the different inspiral regions of the 22 phase. */
  size_t length = pPhase->wflength;
  size_t length_insp_early = pPhase->wflength_insp_early;
  size_t length_insp_late = pPhase->wflength_insp_late;
  size_t length_insp = length_insp_early + length_insp_late;

  /*Initialize time */
  REAL8 t, thetabar, theta, w22;
  REAL8 factheta = pow(5.0,1./8);

  /* Initialize REAL8 sequences for storing the phase and the frequency of the 22 */
  REAL8Sequence *xorb = NULL;
  xorb = XLALCreateREAL8Sequence(length);

  /* Compute x=(0.5*omega_22)^(2/3). 
  Depending on the inspiral region, theta is defined with a fitted t0 parameter or with t0=0. */

  if(pWF->inspVersion!=0) // If reconstruction is non-default, this will compute frequency, phase and PN expansion parameter x from TaylorT3 with fitted t0 for the early inspiral region-
  {
    for(UINT4 jdx = 0; jdx < length_insp_early; jdx++)
    {
      t = pPhase->tmin + jdx*pWF->dtM;

      thetabar = pow(pWF->eta*(pPhase->tt0-t),-1./8);
    
      theta = factheta*thetabar;

      w22 = IMRPhenomTomega22(t, theta, pWF, pPhase);
      xorb->data[jdx] = pow(0.5*w22,2./3);
    }

    for(UINT4 jdx = length_insp_early; jdx < length_insp; jdx++) // For times later than the early-late inspiral boundary, it computes phase, frequency and x with the extended TaylorT3 (with tt0=0)
    {
      t = pPhase->tmin + jdx*pWF->dtM;

      thetabar = pow(-pWF->eta*t,-1./8);
    
      theta = factheta*thetabar;

      w22 = IMRPhenomTomega22(t, theta, pWF, pPhase);
      xorb->data[jdx] = pow(0.5*w22,2./3);
    }

  }

  else // If default reconstruction, only one inspiral region with extended TaylorT3 (with tt0=0) is employed up to the inspiral-merger boundary time
  {
      for(UINT4 jdx = 0; jdx < length_insp; jdx++)
      {
        t = pPhase->tmin + jdx*pWF->dtM;

        thetabar = pow(-pWF->eta*t,-1./8);
    
        theta = factheta*thetabar;

        w22 = IMRPhenomTomega22(t, theta, pWF, pPhase);
        xorb->data[jdx] = pow(0.5*w22,2./3);
      }

  }

  /* During the 22 phase merger region, phase and x are also computed because the higher modes end the inspiral at a fixed time,
  which can occur later than the end of the 22 inspiral. */
  for(UINT4 jdx = length_insp; jdx < length; jdx++)
  {
    t = pPhase->tmin + jdx*pWF->dtM;

    w22 = IMRPhenomTomega22(t, 0.0, pWF, pPhase);
    xorb->data[jdx] = pow(0.5*w22,2./3);
  }

  /* Compute evolution of LNhat and individual spins */
  status = IMRPhenomTPHM_EvolveOrbit(V,S1x,S1y,S1z,S2x,S2y,S2z,LNhatx,LNhaty,LNhatz,E1x,E1y,E1z,xorb,m1_SI,m2_SI,chi1x,chi1y,chi1z,chi2x,chi2y,chi2z,pWF,pPhase);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTPHM_EvolveOrbit has failed.");

  LALFree(pWF);
  LALFree(pPhase);
  XLALDestroyREAL8Sequence(xorb);

  return status;
}
/** @} */
/** @} */

/* Evolve the spin precessing PN evolution equations using the precomputed v(t) as a driver for the evolution. */
int IMRPhenomTPHM_EvolveOrbit(
  REAL8TimeSeries **V,            /**< post-Newtonian parameter [returned]*/
  REAL8TimeSeries **S1x,          /**< Spin1 vector x component [returned]*/
  REAL8TimeSeries **S1y,          /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **S1z,          /**< "    "    "  z component [returned]*/
  REAL8TimeSeries **S2x,          /**< Spin2 vector x component [returned]*/
  REAL8TimeSeries **S2y,          /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **S2z,          /**< "    "    "  z component [returned]*/
  REAL8TimeSeries **LNhatx,       /**< unit orbital ang. mom. x [returned]*/
  REAL8TimeSeries **LNhaty,       /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **LNhatz,       /**< "    "    "  z component [returned]*/
  REAL8TimeSeries **E1x,          /**< orb. plane basis vector x[returned]*/
  REAL8TimeSeries **E1y,          /**< "    "    "  y component [returned]*/
  REAL8TimeSeries **E1z,          /**< "    "    "  z component [returned]*/
  REAL8Sequence *xorb,            /**< squared velocity from IMRPhenomT phase evolution*/
  REAL8 m1_SI,                    /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                    /**< Mass of companion 2 (kg) */
  REAL8 s1x,                      /**< x component of primary spin at reference frequency*/
  REAL8 s1y,                      /**< y component of primary spin at reference frequency*/
  REAL8 s1z,                      /**< z component of primary spin at reference frequency*/
  REAL8 s2x,                      /**< x component of secondary spin at reference frequency*/
  REAL8 s2y,                      /**< y component of secondary spin at reference frequency*/
  REAL8 s2z,                      /**< z component of secondary spin at reference frequency*/
  IMRPhenomTWaveformStruct *pWF,  /**< PhenomTHM waveform struct*/
  IMRPhenomTPhase22Struct *pPhase /**< PhenomT phase struct*/
  )
{
    /* Initialize objects and variables needed for the RK integrator */
    
    LALAdaptiveRungeKuttaIntegrator *integrator = NULL;     /* GSL integrator object */
    INT4 intreturn;         /* Variable for storing status returned by the integrator */
    REAL8 yinit[12];        /* initial values of parameters */
    REAL8Array *yout;       /* time series of variables returned from integrator */
    
    /* intermediate variables */
    int i, len, len2;
    REAL8 norm1, norm2, m1sec, m2sec, Msec;

    /* Standard tests for the input/ouput time series */
    if ( !V || !S1x || !S1y || !S1z || !S2x || !S2y || !S2z
            || !LNhatx || !LNhaty || !LNhatz || !E1x || !E1y || !E1z )
    {
        XLALPrintError("XLAL Error - %s: NULL(s) in output parameters\n",
                       __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* Set integration length and interpolation length */
    //size_t lengthAll = floor((0.0 - pPhase->tmin)/pWF->dtM); // Waveform length up to the peak time: Integration length
    size_t lengthAll2 = xorb->length; // Full waveform length

    /* Interpolation of "orbital" v */
    REAL8Sequence *timesPN = NULL;
    REAL8Sequence *vorb = NULL;
    timesPN = XLALCreateREAL8Sequence(lengthAll2);
    vorb = XLALCreateREAL8Sequence(lengthAll2);

    for(UINT4 idx=0; idx<lengthAll2; idx++)
    {
      timesPN->data[idx] = idx*pWF->dtM;
      vorb->data[idx] = sqrt(xorb->data[idx]);
    }

    // Initialization of gsl_spline objects
    gsl_interp_accel *accel_v;
    gsl_spline *spline_v;
    accel_v = gsl_interp_accel_alloc();
    spline_v = gsl_spline_alloc(gsl_interp_cspline, vorb->length);

    // Interpolation of v(t)
    gsl_spline_init(spline_v, timesPN->data, vorb->data, vorb->length);

    /* Fill needed structs.
    PhenomTPHMEvolution is a struct for the SpinDerivatives function that contains the SpinTaylor struct
    (spin derivatives are universal for all SpinTaylor approximants), the interpolated v(t) and other needed parameters. */

    PhenomTPHMEvolution params;

    XLALSimInspiralSpinTaylorTxCoeffs *Tparams=NULL; // PN coefficients struct

    /* We select fixed values for the different PN orders: No tidal interactions (it is a BBH model), highest PN order for 
    spins but perpendicular corrections to Lhat disabled, since it was not clear the effect on this on the derived Euler angles from Lhat. */

    REAL8 lambda1 = 0.0; REAL8 lambda2 = 0.0; REAL8 quadparam1 = 0.0; REAL8 quadparam2 = 0.0;
    LALSimInspiralSpinOrder spinO = -1;
    LALSimInspiralTidalOrder tideO = 0;
    INT4 phaseO = -1;
    INT4 lscorr = 0;

    XLALSimInspiralSpinTaylorT4Setup(&Tparams, m1_SI, m2_SI, pWF->fRef, 0.0, lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, lscorr);

    params.v_spline = spline_v;
    params.v_acc = accel_v;
    params.Tparams = Tparams;
    params.ToffSign = 0.0;

    /* Initial quantities for evolution */

    // Needed for getting dimensionful spins 
    m1sec = m1_SI / LAL_MSUN_SI * LAL_MTSUN_SI;
    m2sec = m2_SI / LAL_MSUN_SI * LAL_MTSUN_SI;
    Msec = m1sec + m2sec;

    /* Put initial values into a single array for the integrator */

    /* LNh(x,y,z) */
    yinit[0] = 0.0;
    yinit[1] = 0.0;
    yinit[2] = 1.0;
    /* S1(x,y,z) */
    norm1 = m1sec * m1sec / Msec / Msec;
    yinit[3] = norm1 * s1x;
    yinit[4] = norm1 * s1y;
    yinit[5] = norm1 * s1z;
    /* S2(x,y,z) */
    norm2 = m2sec * m2sec / Msec / Msec;
    yinit[6] = norm2 * s2x;
    yinit[7] = norm2 * s2y;
    yinit[8]= norm2 * s2z;
    /* E1(x,y,z) */
    yinit[9] = 1.0;
    yinit[10] = 0.0;
    yinit[11] = 0.0;

    /* initialize the integrator */
    integrator = XLALAdaptiveRungeKutta4Init(12,
                XLALSimIMRPhenomTPHMSpinDerivatives,
                StoppingTest,
                LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    
    if( !integrator )
    {
        XLALPrintError("XLAL Error - %s: Cannot allocate integrator\n", 
                __func__);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* stop the integration only when ending time is achieved. */
    integrator->stopontestonly = 0;

    /********* Integration from tref to t=0 ******/
    /* Formally integration will perform forwards with a positive time array,
    then ending time is select as -tmin and initial time is selected as -tmin + tref */

    REAL8 tend = -pPhase->tmin;
    REAL8 tint1 = fabs(tend + pPhase->tRef);
    size_t length1 = floor(tint1/pWF->dtM);

    /* run the integration; note: time is measured in \hat{t} = t / M.
     Integration is indeed perform until tend - dtM to avoid an extra final point. */
    len = XLALAdaptiveRungeKutta4Hermite(integrator, &params, yinit,
            tint1, tend - pWF->dtM, pWF->dtM, &yout);

    // Check if integration failed.
    if (!yout)
    {
        XLALPrintError("XLAL Error - %s: integration failed (yout == NULL)\n",
                       __func__);
        XLAL_ERROR(XLAL_EFUNC);
    }

    intreturn = integrator->returncode;
    if (!len) 
    {
        XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intreturn);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* epoch of evolved time series same as waveform epoch */
    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO;
    XLALGPSAdd(&ligotimegps_zero, pPhase->tminSec);

    //size_t lengthAll = floor((tend - tint1)/pWF->dtM) + length1;
    size_t lengthAll = len + length1;

    /* allocate memory for output vectors */
    *V = XLALCreateREAL8TimeSeries( "PN_EXPANSION_PARAMETER", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll);  
    *S1x = XLALCreateREAL8TimeSeries( "SPIN1_X_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *S1y = XLALCreateREAL8TimeSeries( "SPIN1_Y_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *S1z = XLALCreateREAL8TimeSeries( "SPIN1_Z_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *S2x = XLALCreateREAL8TimeSeries( "SPIN2_X_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *S2y = XLALCreateREAL8TimeSeries( "SPIN2_Y_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *S2z = XLALCreateREAL8TimeSeries( "SPIN2_Z_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *LNhatx = XLALCreateREAL8TimeSeries( "LNHAT_X_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *LNhaty = XLALCreateREAL8TimeSeries( "LNHAT_Y_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *LNhatz = XLALCreateREAL8TimeSeries( "LNHAT_Z_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *E1x = XLALCreateREAL8TimeSeries( "E1_BASIS_X_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *E1y = XLALCreateREAL8TimeSeries( "E1_BASIS_Y_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    *E1z = XLALCreateREAL8TimeSeries( "E1_BASIS_Z_COMPONENT", &ligotimegps_zero, 0., 
            pWF->deltaT, &lalDimensionlessUnit, lengthAll); 
    if ( !*V || !*S1x || !*S1y || !*S1z || !*S2x || !*S2y || !*S2z
             || !*LNhatx || !*LNhaty || !*LNhatz || !*E1x || !*E1y || !*E1z )
    {
        XLALDestroyREAL8Array(yout);
        LALFree(Tparams);
        XLALAdaptiveRungeKuttaFree(integrator);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* Copy dynamical variables from yout array to output time series.
     * Note the first 'len' members of yout are the time steps. 
     */
    for( i = 0; i < len; i++ )
    { 
        int j = i + length1;
        (*V)->data->data[j]     = vorb->data[j];
        (*LNhatx)->data->data[j]  = yout->data[len+i];
        (*LNhaty)->data->data[j]  = yout->data[2*len+i];
        (*LNhatz)->data->data[j]  = yout->data[3*len+i];
        (*S1x)->data->data[j]     = yout->data[4*len+i]/norm1;
        (*S1y)->data->data[j]     = yout->data[5*len+i]/norm1;
        (*S1z)->data->data[j]     = yout->data[6*len+i]/norm1;
        (*S2x)->data->data[j]     = yout->data[7*len+i]/norm2;
        (*S2y)->data->data[j]     = yout->data[8*len+i]/norm2;
        (*S2z)->data->data[j]     = yout->data[9*len+i]/norm2;
        (*E1x)->data->data[j]     = yout->data[10*len+i];
        (*E1y)->data->data[j]     = yout->data[11*len+i];
        (*E1z)->data->data[j]     = yout->data[12*len+i];
    }

    // Free PN Taylor struct and integrator. If tmin!=tref and a second integration step is needed, will be defined again, to avoid possible overwritting.
    LALFree(Tparams);
    XLALAdaptiveRungeKuttaFree(integrator);
    XLALDestroyREAL8Array(yout);

    /********* Integration from tref to tmin ******/
    /* If tref!=tmin, quantities have to be evolved backwards from tref to tmin.
    Formally this is done by a forward integration, but physically it corresponds to a backward integration. */
    
    if(fabs(pPhase->tRef-pPhase->tmin)>=pWF->dtM)
    {

      /* Initialize objects and variables needed for the RK integrator */

      LALAdaptiveRungeKuttaIntegrator *integrator2 = NULL;
      REAL8Array *yout2;   // time series of variables returned from integrator 
      REAL8 yinit2[12];

      /* Put initial values into a single array for the integrator */

      // LNh(x,y,z) 
      yinit2[0] = 0.0;
      yinit2[1] = 0.0;
      yinit2[2] = 1.0;
      // S1(x,y,z) 
      norm1 = m1sec * m1sec / Msec / Msec;
      yinit2[3] = norm1 * s1x;
      yinit2[4] = norm1 * s1y;
      yinit2[5] = norm1 * s1z;
      // S2(x,y,z) 
      norm2 = m2sec * m2sec / Msec / Msec;
      yinit2[6] = norm2 * s2x;
      yinit2[7] = norm2 * s2y;
      yinit2[8]= norm2 * s2z;
      // E1(x,y,z) 
      yinit2[9] = 1.0;
      yinit2[10] = 0.0;
      yinit2[11] = 0.0;


      /* Set up PN coefficients struct */
      XLALSimInspiralSpinTaylorT4Setup(&Tparams, m1_SI, m2_SI, pWF->fmin, 0.0, lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, lscorr);
      params.Tparams = Tparams;
      // Passed to actually performed backward integration. For t' in the integration, v(t) will be evaluated at t = tint1 - t'
      params.ToffSign = -tint1; 

      /* Set up integrator */
      integrator2 = XLALAdaptiveRungeKutta4Init(12,
                XLALSimIMRPhenomTPHMSpinDerivatives,
                StoppingTest,
                LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);

      /* Stop integration only when ending time is achieved */
      integrator2->stopontestonly = 0;

      /* Perform integration. Formally is done from dtM to tint1 = -tmin + tref, but physically is done from tref to tmin. */
      len2 = XLALAdaptiveRungeKutta4Hermite(integrator2, &params, yinit2,
             pWF->dtM, tint1, pWF->dtM, &yout2);

      /* Check if integration failed */
      if (!yout2)
      {
          XLALPrintError("XLAL Error - %s: integration failed (yout == NULL)\n",
                       __func__);
          XLAL_ERROR(XLAL_EFUNC);
      }
      intreturn = integrator->returncode;
      if (!len2) 
      {
          XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intreturn);
          XLAL_ERROR(XLAL_EFUNC);
      }

      /* Copy dynamical variables from yout array to output time series.
        * Note the first 'len' members of yout are the time steps. */
      for( i = 0; i < len2; i++ )
      { 
          int j = len2 - 1 - i; // Time series are filled backwards from the index corresponding to tref to fill all the elements.
          (*V)->data->data[j]     = vorb->data[j];
          (*LNhatx)->data->data[j]  = yout2->data[len2+i];
          (*LNhaty)->data->data[j]  = yout2->data[2*len2+i];
          (*LNhatz)->data->data[j]  = yout2->data[3*len2+i];
          (*S1x)->data->data[j]     = yout2->data[4*len2+i]/norm1;
          (*S1y)->data->data[j]     = yout2->data[5*len2+i]/norm1;
          (*S1z)->data->data[j]     = yout2->data[6*len2+i]/norm1;
          (*S2x)->data->data[j]     = yout2->data[7*len2+i]/norm2;
          (*S2y)->data->data[j]     = yout2->data[8*len2+i]/norm2;
          (*S2z)->data->data[j]     = yout2->data[9*len2+i]/norm2;
          (*E1x)->data->data[j]     = yout2->data[10*len2+i];
          (*E1y)->data->data[j]     = yout2->data[11*len2+i];
          (*E1z)->data->data[j]     = yout2->data[12*len2+i];
      }

      // Free integrator, PN coef struct and ouput array.
      XLALAdaptiveRungeKuttaFree(integrator2);
      XLALDestroyREAL8Array(yout2);
      LALFree(Tparams);
    }
    
    // Free gsl spline objects
    gsl_spline_free(spline_v);
    gsl_interp_accel_free(accel_v);

    XLALDestroyREAL8Sequence( timesPN);
    XLALDestroyREAL8Sequence( vorb);

    return XLAL_SUCCESS;
}


/* Function to compute Euler angles from the evolution of LNhat:
    By definition:
    - alpha = arctan(LNhat_y / LNhat_x)
    - cosbeta = LNhat*Jhat = LNhat_z (in J-frame)
    - gamma = - \int \dot{alpha}*cosbeta (minimal rotation condition)
 */

int IMRPhenomTPHM_NumericalEulerAngles( 
  REAL8TimeSeries **alphaInt,               /* Alpha angle time series [out] */ 
  REAL8TimeSeries **cosbetaInt,             /* cos(Beta) angle time series [out] */           
  REAL8TimeSeries **gammaInt,               /* Gamma angle time series [out] */
  REAL8 *af_evolved,                        /* Final spin predicted by the evolved spin values at the peak time */
  REAL8Sequence *xorb,                      /* x(t)=v(t)^2 time series [in] */
  REAL8 deltaT,                             /* Sampling interval (s) */
  REAL8 m1_SI,                              /* Mass of companion 1 (kg) */
  REAL8 m2_SI,                              /* Mass of companion 2 (kg) */
  REAL8 s1x,                                /* x component of primary spin */
  REAL8 s1y,                                /* y component of primary spin */
  REAL8 s1z,                                /* z component of primary spin */
  REAL8 s2x,                                /* x component of secondary spin */
  REAL8 s2y,                                /* y component of secondary spin */
  REAL8 s2z,                                /* z component of secondary spin */
  REAL8 epochT,                             /* initial time in mass units */
  IMRPhenomTWaveformStruct *pWF,            /* PhenomTHM waveform struct */
  IMRPhenomTPhase22Struct *pPhase,          /* PhenomTHM phase struct */
  UNUSED IMRPhenomXPrecessionStruct *pPrec  /* PhenomX precessing struct */
  );

int IMRPhenomTPHM_NumericalEulerAngles(
  REAL8TimeSeries **alphaTS,                /* Alpha angle time series [out] */ 
  REAL8TimeSeries **cosbetaTS,              /* cos(Beta) angle time series [out] */
  REAL8TimeSeries **gammaTS,                /* Gamma angle time series [out] */
  REAL8 *af_evolved,                        /* Final spin predicted by the evolved spin values at the peak time */
  REAL8Sequence *xorb,                      /* x(t)=v(t)^2 time series [in] */
  REAL8 deltaT,                             /* Sampling interval (s) */
  REAL8 m1_SI,                              /* Mass of companion 1 (kg) */
  REAL8 m2_SI,                              /* Mass of companion 2 (kg) */
  REAL8 s1x,                                /* x component of primary spin */
  REAL8 s1y,                                /* y component of primary spin */
  REAL8 s1z,                                /* z component of primary spin */
  REAL8 s2x,                                /* x component of secondary spin */
  REAL8 s2y,                                /* y component of secondary spin */
  REAL8 s2z,                                /* z component of secondary spin */
  REAL8 epochT,                             /* initial time in mass units */
  IMRPhenomTWaveformStruct *pWF,            /* PhenomTHM waveform struct */
  IMRPhenomTPhase22Struct *pPhase,          /* PhenomTHM phase struct */
  UNUSED IMRPhenomXPrecessionStruct *pPrec  /* PhenomX precessing struct */
  )
{
  int status = XLAL_SUCCESS;

  /* Initialize needed time series for calling IMRPhenomTPHM_EvolveOrbit */
  REAL8TimeSeries *V = NULL;
  REAL8TimeSeries *S1x = NULL;
  REAL8TimeSeries *S1y = NULL;
  REAL8TimeSeries *S1z = NULL;
  REAL8TimeSeries *S2x = NULL;
  REAL8TimeSeries *S2y = NULL;
  REAL8TimeSeries *S2z = NULL;
  REAL8TimeSeries *LNhatx = NULL;
  REAL8TimeSeries *LNhaty = NULL;
  REAL8TimeSeries *LNhatz = NULL;
  REAL8TimeSeries *E1x = NULL;
  REAL8TimeSeries *E1y = NULL;
  REAL8TimeSeries *E1z = NULL;

  /* Evolve LNhat and spins */
  status = IMRPhenomTPHM_EvolveOrbit(&V, &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z,\
    xorb, m1_SI, m2_SI, s1x, s1y, s1z, s2x, s2y, s2z, pWF, pPhase);

  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTPHM_EvolveOrbit has failed.");

  size_t lenPN = LNhatx->data->length;

  // Compute final spin as predicted by spins at tpeak

  // Needed for getting dimensionful spins 
    REAL8 m1sec = m1_SI / LAL_MSUN_SI * LAL_MTSUN_SI;
    REAL8 m2sec = m2_SI / LAL_MSUN_SI * LAL_MTSUN_SI;
    REAL8 Msec = m1sec + m2sec;
    REAL8 norm1 = m1sec * m1sec / Msec / Msec;
    REAL8 norm2 = m2sec * m2sec / Msec / Msec;

    REAL8 lnxpeak, lnypeak, lnzpeak;
    REAL8 s1xpeak,s1ypeak,s1zpeak,s2xpeak,s2ypeak,s2zpeak;

    // Values of spins and Lhat at the peak time

    s1xpeak = norm1*(S1x)->data->data[lenPN-1];
    s1ypeak = norm1*(S1y)->data->data[lenPN-1];
    s1zpeak = norm1*(S1z)->data->data[lenPN-1];

    s2xpeak = norm2*(S2x)->data->data[lenPN-1];
    s2ypeak = norm2*(S2y)->data->data[lenPN-1];
    s2zpeak = norm2*(S2z)->data->data[lenPN-1];

    lnxpeak = (LNhatx)->data->data[lenPN-1];
    lnypeak = (LNhaty)->data->data[lenPN-1];
    lnzpeak = (LNhatz)->data->data[lenPN-1];

    // Projection of the spins in the Lhat direction and perpendicular direction, following the right hand rule.

    REAL8 s1Lpeak = s1xpeak*lnxpeak + s1ypeak*lnypeak + s1zpeak*lnzpeak;
    REAL8 s2Lpeak = s2xpeak*lnxpeak + s2ypeak*lnypeak + s2zpeak*lnzpeak;

    REAL8 s1xparallel = s1Lpeak*lnxpeak;
    REAL8 s1yparallel = s1Lpeak*lnypeak;
    REAL8 s1zparallel = s1Lpeak*lnzpeak;

    REAL8 s1xperp = s1xpeak - s1xparallel;
    REAL8 s1yperp = s1ypeak - s1yparallel;
    REAL8 s1zperp = s1zpeak - s1zparallel;

    REAL8 s2xparallel = s2Lpeak*lnxpeak;
    REAL8 s2yparallel = s2Lpeak*lnypeak;
    REAL8 s2zparallel = s2Lpeak*lnzpeak;

    REAL8 s2xperp = s2xpeak - s2xparallel;
    REAL8 s2yperp = s2ypeak - s2yparallel;
    REAL8 s2zperp = s2zpeak - s2zparallel;

    // Final spìn from common geometrical formula employing the evolved spins at merger.

    REAL8 Sperp = sqrt( pow(s1xperp + s2xperp, 2) +  pow(s1yperp + s2yperp, 2) + pow(s1zperp + s2zperp, 2) );

    REAL8 af_nonprec = XLALSimIMRPhenomXFinalSpin2017(pWF->eta, s1Lpeak/norm1, s2Lpeak/norm2);

    REAL8 Sf = copysign(1.0,af_nonprec)*sqrt(Sperp*Sperp + pow(af_nonprec,2)); 

    if(Sf>1.0){Sf = 1.0;};
    if(Sf<-1.0){Sf = -1.0;};

    (*af_evolved) = Sf;

  /* Initialize needed sequences for storing times and angles */

  REAL8Sequence *alpha = NULL;
  REAL8Sequence *alphaaux = NULL;
  REAL8Sequence *cosbeta = NULL;
  REAL8Sequence *gamma = NULL;
  REAL8Sequence *times = NULL;

  // Offset for alpha angle, corresponding to the angle of the perpendicular to J projection of L in the plane perpendicular to J at tref.
  REAL8 alphaOff = atan2(pPrec->S1y + pPrec->S2y, pPrec->S1x + pPrec->S2x) - LAL_PI;

  // Waveform lenght and index of tref
  size_t lengthAll = pPhase->wflength;
  REAL8 tend = -pPhase->tmin;
  REAL8 tint1 = fabs(tend + pPhase->tRef);
  size_t length1 = floor(tint1/pWF->dtM); // This corresponds to the index of tref/fref.

  /* Setup sequences for angles and time */
  alpha = XLALCreateREAL8Sequence(lengthAll);
  alphaaux = XLALCreateREAL8Sequence(lengthAll);
  cosbeta = XLALCreateREAL8Sequence(lengthAll);
  gamma = XLALCreateREAL8Sequence(lengthAll);
  times = XLALCreateREAL8Sequence(lengthAll);


  /**** ROTATION OF LNHAT FROM L0 frame to J frame ****/
  /* Rotation described in Appendix C eq. C13-14 of IMRPhenomXPHM paper: https://arxiv.org/pdf/2004.06503.pdf
  Applied to the whole time evolution of LNhat */

  REAL8 cosPhiJ = cos(-pPrec->phiJ_Sf);
  REAL8 sinPhiJ = sin(-pPrec->phiJ_Sf);

  REAL8 cosThetaJ = cos(-pPrec->thetaJ_Sf);
  REAL8 sinThetaJ = sin(-pPrec->thetaJ_Sf);

  REAL8 cosKappa = cos(-pPrec->kappa);
  REAL8 sinKappa = sin(-pPrec->kappa);

  for(UINT8 i=0; i < lenPN; i++) 
  {
    
    IMRPhenomT_rotate_z(cosPhiJ,    sinPhiJ,    &(LNhatx->data->data[i]), &(LNhaty->data->data[i]), &(LNhatz->data->data[i]));
    IMRPhenomT_rotate_y(cosThetaJ,  sinThetaJ,  &(LNhatx->data->data[i]), &(LNhaty->data->data[i]), &(LNhatz->data->data[i]));
    IMRPhenomT_rotate_z(cosKappa,   sinKappa,   &(LNhatx->data->data[i]), &(LNhaty->data->data[i]), &(LNhatz->data->data[i]));

    /* Compute Euler angles in the J frame */
    alphaaux->data[i] = atan2(LNhaty->data->data[i], LNhatx->data->data[i]);
    cosbeta->data[i] = LNhatz->data->data[i];
    times->data[i] = i*deltaT + epochT; 
  }

  // Substract alpha value at tref. Added to alphaOff for comodity.
  alphaOff -= atan2(LNhaty->data->data[length1], LNhatx->data->data[length1]);
  for(UINT8 i=0; i < lenPN; i++)
  {
    alphaaux->data[i] += alphaOff;
  }

  // Unwrap alpha array to have a continuous time series.
  unwrap_array(alphaaux->data, alpha->data, lenPN);

  REAL8 EulerRDslope = GetEulerSlope(Sf, pWF->Mfinal);

  // Compute RD angles approximation. Same as in PNAnalyticalInspiralEulerAngles.
  for(UINT8 jdx = lenPN-1; jdx < lengthAll; jdx++)
   {
        times->data[jdx] = jdx*deltaT + epochT;
        alpha->data[jdx] = alpha->data[lenPN-2] + EulerRDslope*times->data[jdx];
        cosbeta->data[jdx] = cosbeta->data[lenPN-2];
   }

  /**** COMPUTE GAMMA FROM MINIMAL ROTATION CONDITION ****/
  /* alpha and cosbeta are interpolated to compute the RHS of the minimal rotation condition.
  Then gamma is integrated using Boole's rule. */

  // Interpolate alpha
  gsl_interp_accel *accel_alpha;
  gsl_spline *spline_alpha;

  accel_alpha = gsl_interp_accel_alloc();
  spline_alpha = gsl_spline_alloc(gsl_interp_cspline, lengthAll);

  gsl_spline_init(spline_alpha, times->data, alpha->data, lengthAll);

  // Interpolate cosbeta
  gsl_interp_accel *accel_cosb;
  gsl_spline *spline_cosb;

  accel_cosb = gsl_interp_accel_alloc();
  spline_cosb = gsl_spline_alloc(gsl_interp_cspline, lengthAll);

  gsl_spline_init(spline_cosb, times->data, cosbeta->data, lengthAll);

  // Setup helper struct for gamma integration
  gammaIntegration gammaStruct;

  gammaStruct.alpha_spline = spline_alpha;
  gammaStruct.alpha_acc    = accel_alpha;
  gammaStruct.cosbeta_spline  = spline_cosb;
  gammaStruct.cosbeta_acc     = accel_cosb;

  // Offset for gamma (this step is also employed to fill the output time series for alpha and cosbeta from the sequences)
  gamma->data[0] = -alpha->data[0];
  (*gammaTS)->data->data[0] = gamma->data[0];
  (*alphaTS)->data->data[0] = alpha->data[0];
  (*cosbetaTS)->data->data[0] = cosbeta->data[0];

  // Integrate gamma
  for( UINT4 i = 1; i < lengthAll; i++ ){
  double t1 = times->data[i-1];
  double t2 = times->data[i];
  gamma->data[i] = gamma->data[i-1] + (1.0/90.0)*(t2 - t1)*(7.0*f_alphadotcosi(t1,&gammaStruct)
                  + 32.0*f_alphadotcosi((t1+3.0*t2)/4.0,&gammaStruct)
                  + 12.0*f_alphadotcosi((t1+t2)/2.0,&gammaStruct)
                  + 32.0*f_alphadotcosi((3.0*t1+t2)/4.0,&gammaStruct)
                  + 7.0*f_alphadotcosi(t2,&gammaStruct));/* Boole's Rule */

  // Fill output time series
  (*gammaTS)->data->data[i] = gamma->data[i];
  (*alphaTS)->data->data[i] = alpha->data[i];
  (*cosbetaTS)->data->data[i] = cosbeta->data[i];
  }

  // Substract gamma value at tref and apply corresponding offset.
  for( UINT4 i = 1; i < lengthAll; i++ ){
  (*gammaTS)->data->data[i] += - gamma->data[length1] - alpha->data[length1];
  }


  // Destroy time series and sequences
  XLALDestroyREAL8TimeSeries( V );
  XLALDestroyREAL8TimeSeries( S1x );
  XLALDestroyREAL8TimeSeries( S1y );
  XLALDestroyREAL8TimeSeries( S1z );
  XLALDestroyREAL8TimeSeries( S2x );
  XLALDestroyREAL8TimeSeries( S2y );
  XLALDestroyREAL8TimeSeries( S2z );
  XLALDestroyREAL8TimeSeries( LNhatx );
  XLALDestroyREAL8TimeSeries( LNhaty );
  XLALDestroyREAL8TimeSeries( LNhatz );
  XLALDestroyREAL8TimeSeries( E1x );
  XLALDestroyREAL8TimeSeries( E1y );
  XLALDestroyREAL8TimeSeries( E1z );

  XLALDestroyREAL8Sequence( alpha);
  XLALDestroyREAL8Sequence( alphaaux);
  XLALDestroyREAL8Sequence( cosbeta);
  XLALDestroyREAL8Sequence( gamma);
  XLALDestroyREAL8Sequence( times);

  // Free gsl spline objects
  gsl_spline_free(spline_alpha);
  gsl_spline_free(spline_cosb);
  gsl_interp_accel_free(accel_alpha);
  gsl_interp_accel_free(accel_cosb);

return status;
}