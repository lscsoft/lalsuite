/*
*  Copyright (C) 2011 Craig Robinson, Enrico Barausse, Yi Pan, Prayush Kumar
*  (minor changes), Andrea Taracchini
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
#include <complex.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimSphHarmMode.h>
#include <LALSimInspiralWaveformFlags.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"
#include "LALSimInspiralPrecess.h"

/* Include all the static function files we need */
#include "LALSimIMREOBHybridRingdown.c"
#include "LALSimIMREOBFactorizedWaveform.c"
#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMRSpinEOBInitialConditions.c"
#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinAlignedEOBHcapDerivative.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBFactorizedWaveform.c"
#include "LALSimIMRSpinEOBFactorizedFlux.c"
/* OPTIMIZED */
#include "LALSimIMRSpinEOBHamiltonianOptimized.c"
#include "LALSimIMRSpinEOBComputeAmpPhasefromEOMSoln.c"
#include "LALSimIMRSpinAlignedEOBGSLOptimizedInterpolation.c"
#include "LALSimIMRSpinAlignedEOBHcapDerivativeOptimized.c"
/* END OPTIMIZED */


#define debugOutput 0
#define SEOBNRv4HM 41


//static int debugPK = 0;

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


/**
 * ModeArray is a structure which allows to select the modes to include
 * in the waveform.
 * This function will create a structure with the default modes for every model
 */
static INT4 XLALSetup_EOB__std_mode_array_structure(
    LALValue *ModeArray, UINT4 SpinAlignedEOBversion)
{

  /* setup ModeArray */

      if (SpinAlignedEOBversion == SEOBNRv4HM) {
        /* Adding all the modes of SEOBNRv4HM
        * i.e. [(2,2),(2,1),(3,3),(4,4),(5,5)]
        the relative -m modes are added automatically*/
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 5, 5);
      }
      else{
        /*All the other spin aligned model so far only have the 22 mode
        */
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
      }


    return XLAL_SUCCESS;
}

/**
 * ModeArray is a structure which allows to select the modes to include
 * in the waveform.
 * This function check if the selected modes are available for a given model
 */
static INT4 XLALCheck_EOB_mode_array_structure(
    LALValue *ModeArray, UINT4 SpinAlignedEOBversion)
{
  INT4 flagTrue = 0;
  UINT4 modeL;
  UINT4 modeM;
  UINT4 nModes;
  const UINT4 lmModes[5][2] = {{2, 2}, {3, 3}, {2, 1}, {4, 4}, {5, 5}};
  if (SpinAlignedEOBversion == SEOBNRv4HM) {
    /*If one select SEOBNRv4HM all the modes above are selected to check
    */
    nModes = 5;
  }
  else{
    /*If not only the 22 mode is selected to check
    */
    nModes = 1;
  }
  /*Loop over all the possible modes
  *we only check +m modes, when one select (l,m) mode is actually
  *selecting (l,|m|) mode
  */
  for (UINT4 ELL = 2; ELL <= LAL_SIM_L_MAX_MODE_ARRAY; ELL++) {
    for (UINT4 EMM = 0; EMM <= ELL; EMM++) {
      if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ELL, EMM) == 1) {
            for (UINT4 k = 0; k < nModes; k++) {
              modeL  = lmModes[k][0];
              modeM = lmModes[k][1];
              if ((modeL == ELL)&&(modeM == EMM)) {
                flagTrue=1;
              }
            }
            /*For each active mode check if is available for the selected model
            */
            if (flagTrue != 1) {
              XLALPrintError ("Mode (%d,%d) is not available by the model %d.\n", ELL,
                 EMM, SpinAlignedEOBversion);
              return XLAL_FAILURE;
            }
            flagTrue = 0;
          }
      }
    }

    return XLAL_SUCCESS;
}




/**
 * NR fit to the geometric GW frequency M_{total}omega_{22} of a BNS merger,
 * defined by the time when the (2,2) amplitude peaks
 * See Eq.(2) in https://arxiv.org/pdf/1504.01764.pdf with coefficients
 * given by the 3rd row of Table II therein. Compared to NR for 0 <= kappa2T <= 500
 */
static REAL8 XLALSimNSNSMergerFreq(
                                       TidalEOBParams *tidal1, /**< Tidal parameters of body 1 */
                                       TidalEOBParams *tidal2  /**< Tidal parameters of body 2 */
)
{
    REAL8 X1 = tidal1->mByM;
    REAL8 X2 = tidal2->mByM;
    REAL8 lambda1 = tidal1->lambda2Tidal; // Dimensionless quadrupolar tidal deformability normalized to M^5
    REAL8 lambda2 = tidal2->lambda2Tidal; // Dimensionless quadrupolar tidal deformability normalized to M^5

    REAL8 X1v, X2v, lambda1v, lambda2v;
    if ( X1 >= X2 ) {
        X1v = X1;
        X2v = X2;
        lambda1v = lambda1;
        lambda2v = lambda2;
    }
    else {
        X1v = X2;
        X2v = X1;
        lambda1v = lambda2;
        lambda2v = lambda1;
    }
    REAL8 kappa2T = 3*(X2v/X1v*lambda1v + X1v/X2v*lambda2v);
    if ( kappa2T < 0. ) {
        XLAL_ERROR (XLAL_EFUNC);
    }
    else {
        if ( kappa2T > 500. ) {
            kappa2T = 500.;
        }
       return 0.3596*(1. + 0.024384*kappa2T - 0.000017167*kappa2T*kappa2T)/(1. + 0.068865*kappa2T);
    }
}

static int UNUSED
XLALEOBSpinStopCondition (double UNUSED t,/**<< UNUSED */
			  const double values[],  /**<< Dynamical variables */
			  double dvalues[],   /**<<1st time-derivatives of dynamical variables */
			  void *funcParams  /**<< physical parameters*/
  )
{

  SpinEOBParams *params = (SpinEOBParams *) funcParams;
  double omega_x, omega_y, omega_z, omega;
  double r2;

  omega_x = values[1] * dvalues[2] - values[2] * dvalues[1];
  omega_y = values[2] * dvalues[0] - values[0] * dvalues[2];
  omega_z = values[0] * dvalues[1] - values[1] * dvalues[0];

  r2 = values[0] * values[0] + values[1] * values[1] + values[2] * values[2];
  omega =
    sqrt (omega_x * omega_x + omega_y * omega_y + omega_z * omega_z) / r2;

  /* Terminate when omega reaches peak, and separation is < 6M */
  //if ( omega < params->eobParams->omega )
  if (r2 < 36. && omega < params->eobParams->omega)
    {
      return 1;
    }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}


/**
 * Stopping condition for the regular resolution EOB orbital evolution
 * -- stop when reaching max orbital frequency in strong field.
 * At each test,
 * if omega starts to decrease, return 1 to stop evolution;
 * if not, update omega with current value and return GSL_SUCCESS to continue evolution.
 */
static int
XLALEOBSpinAlignedStopCondition (double UNUSED t, /**< UNUSED */
				 const double values[],
						  /**< dynamical variable values */
				 double dvalues[],/**< dynamical variable time derivative values */
				 void *funcParams /**< physical parameters */
  )
{

  REAL8 omega, r;
  SpinEOBParams *params = (SpinEOBParams *) funcParams;

  r = values[0];
  omega = dvalues[1];

//  printf("function 1: r = %.16e, omega = %.16e, pr = %.16e, dpr = %.16e, t = %.16e \n",values[0],dvalues[1],values[2],dvalues[2],t);
  //if ( omega < params->eobParams->omega )
  if (r < 6 && (omega < params->eobParams->omega || dvalues[2] >= 0))
    {
      params->eobParams->omegaPeaked = 0;
      return 1;
    }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

/*
 * Stopping condition for the high resolution EOB orbital evolution
 * -- stop when reaching a minimum radius 0.3M out of the EOB horizon (Eqs. 9b, 37)
 * or when getting nan in any of the four ODE equations
 * At each test,
 * if conditions met, return 1 to stop evolution;
 * if not, return GSL_SUCCESS to continue evolution.
 */
static int
XLALSpinAlignedHiSRStopCondition (double UNUSED t, /**< UNUSED */
				  const double UNUSED values[],/**< dynamical variable values */
				  double dvalues[],	/**< dynamical variable time derivative values */
				  void UNUSED * funcParams     /**< physical parameters */
  )
{
  if (dvalues[0] >= 0. || dvalues[2] >= 0. || isnan (dvalues[3]) || isnan (dvalues[2])
      || isnan (dvalues[1]) || isnan (dvalues[0]))
    {
      return 1;
    }
  return GSL_SUCCESS;
}

/**
 * This function defines the stopping criteria for the tidal EOB model
 * Note that here
 * values[0] = r
 * values[1] = phi
 * values[2] = pr
 * values[3] = pphi
 * dvalues[0] = dr/dt
 * dvalues[1] = dphi/dt
 * dvalues[2] = dpr/dt
 * dvalues[3] = dpphi/dt = omega
 */
static int
XLALSpinAlignedNSNSStopCondition (double UNUSED t, /**< UNUSED */
				    const double UNUSED values[],
							       /**< dynamical variable values */
				    double dvalues[],	/**< dynamical variable time derivative values */
				    void UNUSED * funcParams   /**< physical parameters */
  )
{
  REAL8 omega, r;
  UINT4 counter;
  SpinEOBParams *params = (SpinEOBParams *) funcParams;
  REAL8 rMerger = pow ( params->eobParams->omegaMerger/2., -2./3. );
//  printf("rMerger %.16e\n", rMerger);
  r = values[0];
  omega = dvalues[1];
  counter = params->eobParams->omegaPeaked;
    REAL8 eta = params->eobParams->eta;
  if (debugOutput)   printf("function NSNS: r = %.16e, omega = %.16e, pr = %.16e, dpr = %.16e, count = %.16u\n",values[0],dvalues[1],values[2],dvalues[2],counter);
//  printf("%.16e %.16e %.16e %.16e\n",values[0],dvalues[1],values[2],dvalues[2]);
  REAL8 rCheck = 1.5*rMerger;
  if (r < rCheck && omega < params->eobParams->omega)
    {
      if (debugOutput) printf("Peak detection %.16e %.16e\n", omega, params->eobParams->omega);
      params->eobParams->omegaPeaked = counter + 1;
    }
  if ( omega >= params->eobParams->omegaMerger/2. ) {
      if (debugOutput) printf("Stop at Tim's freq at r=%.16e\n", r);
      return 1;
  }
   if ( r < rCheck && values[2] >= 0 ) {
        if (debugOutput) printf("Stop at pr >= 0 at r=%.16e\n", r);
        return 1;
   }
    if ( r < rCheck && dvalues[0] >= 0 ) {
        if (debugOutput) printf("Stop at dr/dt >= 0 at r=%.16e\n", r);
        return 1;
    }
   if ( r < rCheck && dvalues[2] >= 0 ) {
       params->eobParams->omegaPeaked = counter + 1;
       if (debugOutput) printf("Stop at dpr/dt >= 0 at r=%.16e\n", r);
       // return 1;
  }
  if (r < rCheck && params->eobParams->omegaPeaked == 3 )
    {
     params->eobParams->omegaPeaked = 0;
     if (debugOutput) printf("params->eobParams->omegaPeaked == 3 at r=%.16e\n", r);
     return 1;
  }
    if ( isnan (dvalues[3]) || isnan (dvalues[2]) || isnan (dvalues[1]) || isnan (dvalues[0]) ) {
        if (debugOutput) printf("Stop at nan's at r=%.16e\n", r);
        return 1;
    }
    if ( fabs(r/params->eobParams->rad - 1.) < 1.e-3*64./5.*eta/(r*r*r*r)*0.02 && fabs(r/params->eobParams->rad - 1.) > 0.) {
        if (debugOutput) printf("Radius is stalling at r=%.16e and rad=%.16e\n", r, params->eobParams->rad);
        return 1;
    }
  params->eobParams->omega = omega;
  params->eobParams->rad = r;
    if( LAL_PI/params->deltaT <= 2.*omega ) {
        params->eobParams->NyquistStop = 1;
        if (debugOutput) printf("Stop at Nyquist at r=%.16e\n", r);
        XLAL_PRINT_WARNING ("Waveform will be generated only up to half the sampling frequency, thus discarding any physical higher-frequency contect above that!\n");
        return 1;
    }
  return GSL_SUCCESS;
}

static int
XLALSpinAlignedHiSRStopConditionV4 (double UNUSED t, /**< UNUSED */
                                    const double UNUSED values[],
                                    /**< dynamical variable values */
                                    double dvalues[],	/**< dynamical variable time derivative values */
                                    void UNUSED * funcParams   /**< physical parameters */
)
{
    REAL8 omega, r;
    UINT4 counter;
    SpinEOBParams *params = (SpinEOBParams *) funcParams;
    r = values[0];
    omega = dvalues[1];
    counter = params->eobParams->omegaPeaked;
    //printf("function 2: r = %.16e, omega = %.16e, pr = %.16e, dpr = %.16e, count = %.16u \n",values[0],dvalues[1],values[2],dvalues[2],counter);
    if (r < 6. && omega < params->eobParams->omega)
    {
        //        printf("Peak detection %.16e %.16e\n", omega, params->eobParams->omega);
        params->eobParams->omegaPeaked = counter + 1;
    }
    if (dvalues[2] >= 0. || params->eobParams->omegaPeaked == 5
        || isnan (dvalues[3]) || isnan (dvalues[2]) || isnan (dvalues[1])
        || isnan (dvalues[0]))
    {
        //        if ( dvalues[2] >= 0 ) printf("dvalues[2] >= 0\n");
        //        if ( params->eobParams->omegaPeaked == 5 ) printf("params->eobParams->omegaPeaked == 5\n");
        //        if ( isnan( dvalues[3] ) || isnan (dvalues[2]) || isnan (dvalues[1]) || isnan (dvalues[0]) ) printf("%.16e %.16e %.16e %.16e\n", dvalues[0], dvalues[1], dvalues[2], dvalues[3]);
        return 1;
    }
    params->eobParams->omega = omega;
    return GSL_SUCCESS;
}

/**
 * @addtogroup LALSimIMRSpinAlignedEOB_c
 *
 * @author Craig Robinson, Yi Pan
 *
 * @brief Functions for producing SEOBNRv1 and v2 waveforms.
 *
 * Functions for producing SEOBNRv1 waveforms for
 * spinning binaries, as described in
 * Taracchini et al. ( PRD 86, 024011 (2012), arXiv 1202.0790 ).
 * All equation numbers in this file refer to equations of this paper,
 * unless otherwise specified.
 *
 * @review SEOBNRv1 has been reviewd by Riccardo Sturani, B. Sathyaprakash and Prayush Kumar.
 * The review concluded in fall 2012.
 *
 * Functions for producing SEOBNRv2 waveforms for
 * spinning binaries, as described in
 * Taracchini et al. ( arXiv 1311.2544 ).
 *
 * @review SEOBNRv2 has been reviewed by Riccardo Sturani, Prayush Kumar and Stas Babak.
 * The review concluded with git hash 5bc6bb861de2eb72ca403b9e0f529d83080490fe (August 2014).
 *
 * @{
 */

/**
 * This function returns the frequency at which the peak amplitude occurs
 * in SEOBNRv(x)
 *
 */
double
XLALSimIMRSpinAlignedEOBPeakFrequency (REAL8 m1SI,
				/**< mass of companion 1 (kg) */
				       REAL8 m2SI,
				/**< mass of companion 2 (kg) */
				       const REAL8 spin1z,
				/**< z-component of the dimensionless spin of object 1 */
				       const REAL8 spin2z,
				/**< z-component of the dimensionless spin of object 2 */
				       UINT4 SpinAlignedEOBversion
				/**< 1 for SEOBNRv1, 2 for SEOBNRv2, 4 for SEOBNRv4 */
  )
{

  /* The return variable, frequency in Hz */
  double retFreq;

  REAL8 nrOmega;

  /* The mode to use; currently, only l=2, m=2 is supported */
  INT4 ll = 2;
  INT4 mm = 2;

  /* Mass parameters (we'll work in units of solar mass) */
  REAL8 m1 = m1SI / LAL_MSUN_SI;
  REAL8 m2 = m2SI / LAL_MSUN_SI;
  REAL8 Mtotal = m1 + m2;
  REAL8 eta = m1 * m2 / (Mtotal * Mtotal);

  /* We need spin vectors for the SigmaKerr function */
  REAL8Vector *sigmaKerr = XLALCreateREAL8Vector (3);
  int ii;
  REAL8 aa;

  REAL8Vector s1Vec, s2Vec;
  s1Vec.length = s2Vec.length = 3;
  REAL8 spin1[3] = { 0., 0., spin1z };
  REAL8 spin2[3] = { 0., 0., spin2z };
  s1Vec.data = spin1;
  s2Vec.data = spin2;
  /* the SigmaKerr function uses spins, not chis */
  for (ii = 0; ii < 3; ii++)
    {
      s1Vec.data[ii] *= m1 * m1;
      s2Vec.data[ii] *= m2 * m2;
    }

  /* Calculate the normalized spin of the deformed-Kerr background */
  if (XLALSimIMRSpinEOBCalculateSigmaKerr (sigmaKerr, m1, m2, &s1Vec, &s2Vec)
      == XLAL_FAILURE)
    {
      XLALDestroyREAL8Vector (sigmaKerr);
      XLAL_ERROR (XLAL_EFUNC);
    }

  aa = sigmaKerr->data[2];

  /* Now get the frequency at the peak amplitude */
  switch (SpinAlignedEOBversion)
    {
    case 1:
      nrOmega = XLALSimIMREOBGetNRSpinPeakOmega (ll, mm, eta, aa);
      break;
    case 2:
      nrOmega = XLALSimIMREOBGetNRSpinPeakOmegav2 (ll, mm, eta, aa);
      break;
    case 4:
      nrOmega = XLALSimIMREOBGetNRSpinPeakOmegaV4 (ll, mm, eta, aa/(1. - 2.*eta));
      break;
    default:
      XLALPrintError
	("XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1, v2 and v4 are available.\n",
	 __func__);
      XLAL_ERROR (XLAL_EINVAL);
      break;
    }

  retFreq = nrOmega / (2 * LAL_PI * Mtotal * LAL_MTSUN_SI);

  /* Free memory */
  XLALDestroyREAL8Vector (sigmaKerr);

  return retFreq;
}



int
XLALSimIMRSpinAlignedEOBWaveform (REAL8TimeSeries ** hplus,	     /**<< OUTPUT, +-polarization waveform */
				  REAL8TimeSeries ** hcross,	     /**<< OUTPUT, x-polarization waveform */
				  const REAL8 phiC,		     /**<< coalescence orbital phase (rad) */
				  REAL8 deltaT,			     /**<< sampling time step */
				  const REAL8 m1SI,		     /**<< mass-1 in SI unit */
				  const REAL8 m2SI,		     /**<< mass-2 in SI unit */
				  const REAL8 fMin,		     /**<< starting frequency of the 22 mode (Hz) */
				  const REAL8 r,		     /**<< distance in SI unit */
				  const REAL8 inc,		     /**<< inclination angle */
				  const REAL8 spin1z,		     /**<< z-component of spin-1, dimensionless */
				  const REAL8 spin2z,		      /**<< z-component of spin-2, dimensionless */
				  UINT4 SpinAlignedEOBversion,		      /**<< 1 for SEOBNRv1, 2 for SEOBNRv2, 4 for SEOBNRv4, 201 for SEOBNRv2T, 401 for SEOBNRv4T, 41 for SEOBNRv4HM */
                  LALDict *LALParams /**<< Dictionary of additional wf parameters, including tidal and nonGR */
  )
{
  int ret;

  REAL8 lambda2Tidal1 = 0;
  REAL8 omega02Tidal1 = 0;
  REAL8 lambda3Tidal1 = 0;
  REAL8 omega03Tidal1 = 0;
  REAL8 lambda2Tidal2 = 0;
  REAL8 omega02Tidal2 = 0;
  REAL8 lambda3Tidal2 = 0;
  REAL8 omega03Tidal2 = 0;
  REAL8 quadparam1 = 0;
  REAL8 quadparam2 = 0;

  lambda2Tidal1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(LALParams);
  lambda2Tidal2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(LALParams);
  if ( (SpinAlignedEOBversion == 201 || SpinAlignedEOBversion == 401) && lambda2Tidal1 != 0. ) {
      omega02Tidal1 = XLALSimInspiralWaveformParamsLookupTidalQuadrupolarFMode1(LALParams);
      if ( omega02Tidal1 == 0. ) {
          XLALPrintError ("XLAL Error - %s: Tidal parameters are not set correctly! Always provide non-zero f-mode frequency when lambda2 is non-zero!\n", __func__);
          XLAL_ERROR (XLAL_EDOM);
      }
      lambda3Tidal1 = XLALSimInspiralWaveformParamsLookupTidalOctupolarLambda1(LALParams);
      omega03Tidal1 = XLALSimInspiralWaveformParamsLookupTidalOctupolarFMode1(LALParams);
      if ( lambda3Tidal1 != 0. && omega03Tidal1 == 0. ) {
          XLALPrintError ("XLAL Error - %s: Tidal parameters are not set correctly! Always provide non-zero octupolar f-mode frequency when lambda3 is non-zero!\n", __func__);
          XLAL_ERROR (XLAL_EDOM);
      }
  }
  if ( (SpinAlignedEOBversion == 201 || SpinAlignedEOBversion == 401) && lambda2Tidal2 != 0. ) {
      omega02Tidal2 = XLALSimInspiralWaveformParamsLookupTidalQuadrupolarFMode2(LALParams);
      if ( omega02Tidal2 == 0. ) {
          XLALPrintError ("XLAL Error - %s: Tidal parameters are not set correctly! Always provide non-zero  f-mode frequency when lambda2 is non-zero!\n", __func__);
          XLAL_ERROR (XLAL_EDOM);
      }
      lambda3Tidal2 = XLALSimInspiralWaveformParamsLookupTidalOctupolarLambda2(LALParams);
      omega03Tidal2 = XLALSimInspiralWaveformParamsLookupTidalOctupolarFMode2(LALParams);
      if ( lambda3Tidal2 != 0. && omega03Tidal2 == 0. ) {
          XLALPrintError ("XLAL Error - %s: Tidal parameters are not set correctly! Always provide non-zero octupolar f-mode frequency when lambda3 is non-zero!\n", __func__);
          XLAL_ERROR (XLAL_EDOM);
      }
    }
  quadparam1 = 1. + XLALSimInspiralWaveformParamsLookupdQuadMon1(LALParams);
  quadparam2 = 1. + XLALSimInspiralWaveformParamsLookupdQuadMon2(LALParams);

  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(LALParams);
  /*ModeArray includes the modes chosen by the user
  */
  if (ModeArray == NULL) {
    /*If the user doesn't choose the modes, use the standard modes
    */
    ModeArray = XLALSimInspiralCreateModeArray();
    XLALSetup_EOB__std_mode_array_structure(ModeArray, SpinAlignedEOBversion);
  }
  /*Check that the modes chosen are available for the given model
  */
  if(XLALCheck_EOB_mode_array_structure(ModeArray, SpinAlignedEOBversion) == XLAL_FAILURE){
    XLALPrintError ("Not available mode chosen.\n");
    XLAL_ERROR(XLAL_EFUNC);
  }



  REAL8Vector *nqcCoeffsInput = XLALCreateREAL8Vector(10);
  INT4 nqcFlag = 0;


  if ( SpinAlignedEOBversion == 401 ) {
      nqcFlag = 1;
      REAL8 m1BH, m2BH;
      m1BH = m1SI / (m1SI + m2SI) * 50. * LAL_MSUN_SI;
      m2BH = m2SI / (m1SI + m2SI) * 50. * LAL_MSUN_SI;
#if debugOutput
      printf("First run SEOBNRv4 to compute NQCs\n");
#endif
      ret = XLALSimIMRSpinAlignedEOBWaveformAll (hplus, hcross, phiC, 1./32768, m1BH, m2BH, 2*pow(10.,-1.5)/(2.*LAL_PI)/((m1BH + m2BH)*LAL_MTSUN_SI/LAL_MSUN_SI), r, inc, spin1z, spin2z, 400,
					 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, nqcCoeffsInput, nqcFlag, ModeArray);
      if (ret == XLAL_FAILURE){
        if ( nqcCoeffsInput ) XLALDestroyREAL8Vector( nqcCoeffsInput );
        if(ModeArray) XLALDestroyValue(ModeArray);
        XLAL_ERROR(XLAL_EFUNC);
      }
      nqcFlag = 2;
  }
#if debugOutput
    printf("Generate EOB wf\n");
#endif
    {
      //REAL8Vector *nqcCoeffsInput = XLALCreateREAL8Vector(10);
      //INT4 nqcFlag = 0;
      ret = XLALSimIMRSpinAlignedEOBWaveformAll (hplus, hcross,
                                                 phiC, deltaT, m1SI, m2SI, fMin, r, inc, spin1z, spin2z, SpinAlignedEOBversion,
                                                 lambda2Tidal1, lambda2Tidal2,
                                                 omega02Tidal1, omega02Tidal2,
                                                 lambda3Tidal1, lambda3Tidal2,
                                                 omega03Tidal1, omega03Tidal2,
                                                 quadparam1, quadparam2,
                                                 nqcCoeffsInput, nqcFlag, ModeArray);
     if (ret == XLAL_FAILURE){
       if ( nqcCoeffsInput ) XLALDestroyREAL8Vector( nqcCoeffsInput );
       if (ModeArray) XLALDestroyValue(ModeArray);
       XLAL_ERROR(XLAL_EFUNC);
     }
      if(ModeArray) XLALDestroyValue(ModeArray); 
          
      if ( nqcCoeffsInput )
        XLALDestroyREAL8Vector( nqcCoeffsInput );
    }
  return ret;
}

/**
 * This function generates spin-aligned SEOBNRv1,2,2opt,4,4opt,2T,4T,4HM complex modes hlm.
 * Currently, only the h22 harmonic is available for all the models with the exception of SEOBNRv4HM
 * which contains also the modes hlm =  ((2,1),(3,3),(4,4),(5,5)) besides the (2,2). For this model
 * all available harmonics are generated, one cannot choose which harmonic to generate.
 * STEP 0) Prepare parameters, including pre-computed coefficients
 * for EOB Hamiltonian, flux and waveform
 * STEP 1) Solve for initial conditions
 * STEP 2) Evolve EOB trajectory until reaching the peak of orbital frequency
 * STEP 3) Step back in time by tStepBack and volve EOB trajectory again
 * using high sampling rate, stop at 0.3M out of the "EOB horizon".
 * STEP 4) Locate the peak of orbital frequency for NQC and QNM calculations
 * STEP 5) Calculate NQC correction using hi-sampling data
 * STEP 6) Calculate QNM excitation coefficients using hi-sampling data
 * STEP 7) Generate full inspiral mode using desired sampling frequency
 * STEP 8) Generate full IMR modes -- attaching ringdown to inspiral
 * STEP 9) Generate full IMR modes
 * Note that sanity checks on merger for SEOBNRv4 have revealed that for
 * eta<=0.15 and chi1>0.95 about 0.04% of the waveforms display either
 * very shallow double amplitude peak or slightly negave time-derivative of
 * the GW freq at merger.
 * Note that SEOBNRv2T and SEOBNRv4T can display similar features. The model was
 * validated on the range q=[1,3], Sz=[-0.5,0.5], Lambda2=[0,5000]. Waveforms
 * will not fail on Sz=[-0.7,0.7], but with possibly stronger unwanted features.
 * The initial conditions solver can also fail for low starting frequencies,
 * with a failure rate of ~0.3% at fmin=10Hz for M=3Msol.
 */
int
XLALSimIMRSpinAlignedEOBModes (SphHarmTimeSeries ** hlmmode,
				     /**<< OUTPUT, mode hlm */
             //SM
             REAL8Vector ** dynamics_out, /**<< OUTPUT, low-sampling dynamics */
             REAL8Vector ** dynamicsHi_out, /**<< OUTPUT, high-sampling dynamics */
             //SM
				     REAL8 deltaT,
				     /**<< sampling time step */
				     const REAL8 m1SI,
				     /**<< mass-1 in SI unit */
				     const REAL8 m2SI,
				     /**<< mass-2 in SI unit */
				     const REAL8 fMin,
				     /**<< starting frequency of the 22 mode (Hz) */
				     const REAL8 r,
				     /**<< distance in SI unit */
				     const REAL8 spin1z,
				     /**<< z-component of spin-1, dimensionless */
				     const REAL8 spin2z,
				      /**<< z-component of spin-2, dimensionless */
                     UINT4 SpinAlignedEOBversion,
                     /**<< 1 for SEOBNRv1, 2 for SEOBNRv2, 4 for SEOBNRv4, 201 for SEOBNRv2T, 401 for SEOBNRv4T, 41 for SEOBNRv4HM */
				     const REAL8 lambda2Tidal1,
                     /**<< dimensionless adiabatic quadrupole tidal deformability for body 1 (2/3 k2/C^5) */
				     const REAL8 lambda2Tidal2,
                     /**<< dimensionless adiabatic quadrupole tidal deformability for body 2 (2/3 k2/C^5) */
				     const REAL8 omega02Tidal1,
                     /**<< quadrupole f-mode angular freq for body 1 m_1*omega_{02,1}*/
				     const REAL8 omega02Tidal2,
                      /**<< quadrupole f-mode angular freq for body 2 m_2*omega_{02,2}*/
				     const REAL8 lambda3Tidal1,
                     /**<< dimensionless adiabatic octupole tidal deformability for body 1 (2/15 k3/C^7) */
				     const REAL8 lambda3Tidal2,
                     /**<< dimensionless adiabatic octupole tidal deformability for body 2 (2/15 k3/C^7) */
				     const REAL8 omega03Tidal1,
                     /**<< octupole f-mode angular freq for body 1 m_1*omega_{03,1}*/
				     const REAL8 omega03Tidal2,
                     /**<< octupole f-mode angular freq for body 2 m_2*omega_{03,2}*/
             const REAL8 quadparam1,
                     /**<< parameter kappa_1 of the spin-induced quadrupole for body 1, quadrupole is Q_A = -kappa_A m_A^3 chi_A^2 */
				     const REAL8 quadparam2,
                     /**<< parameter kappa_2 of the spin-induced quadrupole for body 2, quadrupole is Q_A = -kappa_A m_A^3 chi_A^2 */
                     REAL8Vector *nqcCoeffsInput,
                     /**<< Input NQC coeffs */
                     const INT4 nqcFlag
                     /**<< Flag to tell the code to use the NQC coeffs input thorugh nqcCoeffsInput */
  )
{
  REAL8 STEP_SIZE = STEP_SIZE_CALCOMEGA;
  INT4 use_tidal = 0;
  if ( (lambda3Tidal1 != 0. && lambda2Tidal1 == 0.) || (lambda3Tidal2 != 0. && lambda2Tidal2 == 0.) ) {
      XLALPrintError ("XLAL Error - %s: Tidal parameters are not set correctly! You must have a non-zero lambda2 if you provide a non-zero lambda3!\n",
                      __func__);
      XLAL_ERROR (XLAL_EDOM);
  }
  if ( SpinAlignedEOBversion==201 || SpinAlignedEOBversion==401 )
    {
        if ( (lambda2Tidal1 != 0. && omega02Tidal1 == 0.)
            || (lambda2Tidal2 != 0. && omega02Tidal2 == 0.) ) {
            XLALPrintError ("XLAL Error - %s: Tidal parameters are not set correctly! Always provide non-zero compactness and f-mode frequency when lambda2 is non-zero!\n",
                 __func__);
            XLAL_ERROR (XLAL_EDOM);
        }
        if ( (lambda3Tidal1 != 0. && omega03Tidal1 == 0.)
            || (lambda3Tidal2 != 0. && omega03Tidal2 == 0.) ) {
            XLALPrintError ("XLAL Error - %s: Tidal parameters are not set correctly! Always provide non-zero compactness and f-mode frequency when lambda3 is non-zero!\n",
                            __func__);
            XLAL_ERROR (XLAL_EDOM);
        }
        if (fmax(m1SI/m2SI, m2SI/m1SI)>3.) {
            XLALPrintWarning("XLAL Warning - %s: Generating waveform with mass ratio outside the recommended range q=[1,3].\n", __func__);
        }
        if (fabs(spin1z)>0.5 || fabs(spin2z)>0.5) {
            XLALPrintWarning("XLAL Warning - %s: Generating waveform with at least one aligned spin component outside the recommended range chi=[-0.5,0.5].\n", __func__);
        }
        if (lambda2Tidal1>5000. || lambda2Tidal2>5000.) {
            XLALPrintWarning("XLAL Warning - %s: Generating waveform with at least one quadrupole tidal deformability outside the recommended range Lambda2=[0,5000].\n", __func__);
        }
        if ((m1SI+m2SI)/LAL_MSUN_SI*LAL_MTSUN_SI*fMin<1.477e-4) {
            XLALPrintWarning("XLAL Warning - %s: Generating waveform with a low starting frequency. Initial conditions solver can fail when pushing the model to lower frequencies, with a rate of failure ~0.3%% at Mf~1.5e-4 (10Hz for M=3Msol).\n", __func__);
        }
      use_tidal = 1;
      if ( SpinAlignedEOBversion==201 )
          SpinAlignedEOBversion=2;
      if ( SpinAlignedEOBversion==401 )
          SpinAlignedEOBversion=4;
    }
  INT4 use_optimized_v2_or_v4 = 0;
  /* If we want SEOBNRv2_opt, then reset SpinAlignedEOBversion=2 and set use_optimized_v2_or_v4=1 */
  if (SpinAlignedEOBversion == 200)
    {
      SpinAlignedEOBversion = 2;
      use_optimized_v2_or_v4 = 1;
    }
  /* If we want SEOBNRv4_opt, then reset SpinAlignedEOBversion=4 and set use_optimized_v4=1 */
  if (SpinAlignedEOBversion == 400)
    {
      SpinAlignedEOBversion = 4;
      use_optimized_v2_or_v4 = 1;
    }

   INT4 use_hm = 0;
    /* The list of available modes */
//    const UINT4 lmModes[5][2] = {{2, 2}, {2, 1}, {3, 3}, {4, 4}, {5, 5}};
    const UINT4 lmModes[5][2] = {{2, 2}, {3, 3}, {2, 1}, {4, 4},{5, 5}};
    REAL8Vector *hLMAllHi = NULL;
    REAL8Vector *hLMAll = NULL;
    UINT4 nModes = 1;
    /* If we want SEOBNRv4HM, then reset SpinAlignedEOBversion=4 and set use_hm=1 */
    if (SpinAlignedEOBversion == 41)
    {
        SpinAlignedEOBversion = 4;
        use_hm = 1;
        nModes = 5;
    }

  /* If the EOB version flag is neither 1, 2, nor 4, exit */
  if (SpinAlignedEOBversion != 1 && SpinAlignedEOBversion != 2
      && SpinAlignedEOBversion != 4)
    {
      XLALPrintError
	("XLAL Error - %s: SEOBNR version flag incorrectly set to %u\n",
	 __func__, SpinAlignedEOBversion);
      XLAL_ERROR (XLAL_EERR);
    }

  Approximant SpinAlignedEOBapproximant;
  switch (SpinAlignedEOBversion)
    {
    case 1:
      SpinAlignedEOBapproximant = SEOBNRv1;
      break;
    case 2:
      SpinAlignedEOBapproximant = SEOBNRv2;
      break;
    case 4:
      SpinAlignedEOBapproximant = SEOBNRv4;
      break;
    default:
      XLALPrintError
	("XLAL Error - %s: SEOBNR version flag incorrectly set to %u\n",
	 __func__, SpinAlignedEOBversion);
      XLAL_ERROR (XLAL_EERR);
      break;
    }

  /*
   * Check spins
   */
  if (spin1z < -1.0 || spin2z < -1.0)
    {
      XLALPrintError ("XLAL Error - %s: Component spin less than -1!\n",
		      __func__);
      XLAL_ERROR (XLAL_EINVAL);
    }
  if (spin1z > 1.0 || spin2z > 1.0)
    {
      XLALPrintError ("XLAL Error - %s: Component spin bigger than 1!\n",
          __func__);
      XLAL_ERROR (XLAL_EINVAL);
    }
  /* If either spin > 0.6, model not available, exit */
  if (SpinAlignedEOBversion == 1 && (spin1z > 0.6 || spin2z > 0.6))
    {
      XLALPrintError
	("XLAL Error - %s: Component spin larger than 0.6!\nSEOBNRv1 is only available for spins in the range -1 < a/M < 0.6.\n",
	 __func__);
      XLAL_ERROR (XLAL_EINVAL);
    }
  /* For v2 the upper bound is 0.99 */
  if ((SpinAlignedEOBversion == 2) && (spin1z > 0.99 || spin2z > 0.99))
    {
      XLALPrintError
	("XLAL Error - %s: Component spin larger than 0.99!\nSEOBNRv2 and SEOBNRv2_opt are only available for spins in the range -1 < a/M < 0.99.\n",
	 __func__);
      XLAL_ERROR (XLAL_EINVAL);
    }
  //R.C: region where SEOBNRv4HM can be generated, see the test on the review page
  //here  https://git.ligo.org/waveforms/reviews/SEOBNRv4HM/wikis/visual-inspection#check-if-there-are-cases-for-which-the-maximum-amplitude-of-the-22-mode-is-smaller-than-the-one-of-any-of-the-higher-order-modes 
  if(use_hm){
    if(m1SI>m2SI){
      if((m1SI/m2SI > 57.) && (spin1z> 0.96)){
        XLALPrintError
	      ("XLAL Error - %s: Component spin of the more massive BH larger than 0.96 for mass ratio bigger than 57.!\nSEOBNRv4HM \
        for mass ratio bigger than 57 is only available when the spin on the more massive BH in the range -1 < a/M < 0.96.\n",__func__);
        XLAL_ERROR (XLAL_EINVAL);
      }
    }
    if(m2SI>m1SI){
      if((m2SI/m1SI > 57.) && (spin2z> 0.96)){
        XLALPrintError
	      ("XLAL Error - %s: Component spin of the more massive BH larger than 0.96 for mass ratio bigger than 57.!\nSEOBNRv4HM \
        for mass ratio bigger than 57 is only available when the spin on the more massive BH in the range -1 < a/M < 0.96.\n",__func__);
        XLAL_ERROR (XLAL_EINVAL);
      }
    }  
  }  

  INT4 i;

  REAL8Vector *values = NULL;

  /* EOB spin vectors used in the Hamiltonian */
  REAL8Vector *sigmaStar = NULL;
  REAL8Vector *sigmaKerr = NULL;
  REAL8 a, tplspin;
  REAL8 chiS, chiA;

  /* Wrapper spin vectors used to calculate sigmas */
  REAL8Vector s1Vec, s1VecOverMtMt;
  REAL8Vector s2Vec, s2VecOverMtMt;
  REAL8 spin1[3] = { 0, 0, spin1z };
  REAL8 spin2[3] = { 0, 0, spin2z };
  REAL8 s1Data[3], s2Data[3], s1DataNorm[3], s2DataNorm[3];

  /* Parameters of the system */
  REAL8 m1, m2, mTotal, eta, mTScaled;
  REAL8 amp0;
  LIGOTimeGPS tc = LIGOTIMEGPSZERO;

  /* Dynamics of the system */
  REAL8Vector rVec, phiVec, prVec, pPhiVec;
  /* OPTIMIZED */
  REAL8Vector ampVec, phaseVec;
  ampVec.data = NULL;
  phaseVec.data = NULL;
  /* END OPTIMIZED */

  REAL8 omega, v;
  REAL8Vector *hamVHi;
  REAL8Vector *hamV = NULL;
  REAL8Vector *omegaVec = NULL;

  /* SEOBNRv4HM modes */
  INT4 modeL, modeM;

  /* Cartesian vectors needed to calculate Hamiltonian */
  REAL8Vector cartPosVec, cartMomVec;
  REAL8 cartPosData[3], cartMomData[3];

  /* Signal mode */
  COMPLEX16 hLM, hT;
  REAL8Vector *sigReVec = NULL, *sigImVec = NULL;

  /* Non-quasicircular correction */
  EOBNonQCCoeffs nqcCoeffs;
  COMPLEX16 hNQC;
  REAL8Vector *ampNQC = NULL, *phaseNQC = NULL;

  /* Ringdown freq used to check the sample rate */
  COMPLEX16Vector modefreqVec;
  COMPLEX16 modeFreq;

  /* We will have to switch to a high sample rate for ringdown attachment */
  REAL8 deltaTHigh;
  UINT4 resampFac;
  UINT4 resampPwr;
  REAL8 resampEstimate;

  /* How far will we have to step back to attach the ringdown? */
  REAL8 tStepBack;
  INT4 nStepBack;

  /* Dynamics and details of the high sample rate part used to attach the ringdown */
  UINT4 hiSRndx;
  REAL8Vector timeHi, rHi, phiHi, prHi, pPhiHi;
  REAL8Vector *sigReHi = NULL, *sigImHi = NULL;
  REAL8Vector *omegaHi = NULL;

  /* Indices of peak frequency and final point */
  /* Needed to attach ringdown at the appropriate point */
  UINT4 peakIdx = 0, finalIdx = 0;

  /* Variables for the integrator */
  LALAdaptiveRungeKuttaIntegrator *integrator = NULL;
  REAL8Array *dynamics = NULL;
  REAL8Array *dynamicsHi = NULL;

  REAL8Array *dynamicstmp = NULL;	// DAVIDS: MOVING THIS OUTSIDE IF-BLOCK FOR NOW
  REAL8Array *dynamicsHitmp = NULL;	// DAVIDS: MOVING THE VARIABLE DECLARATION OUTSIDE THE IF-BLOCK FOR NOW
  INT4 retLen_fromOptStep2 = 0;	// DAVIDS: ONLY SETTING TO ZERO TO SUPRESS COMPILER WARNING
  INT4 retLen_fromOptStep3 = 0;	// DAVIDS: ^ DITTO

  INT4 retLen;
  REAL8 UNUSED tMax;

  /* Accuracies of adaptive Runge-Kutta integrator */
  const REAL8 EPS_ABS = 1.0e-10;
  const REAL8 EPS_REL = 1.0e-9;	// Davids: changed exponent from -9 to -10

  /*
   * STEP 0) Prepare parameters, including pre-computed coefficients
   *         for EOB Hamiltonian, flux and waveform
   */

  /* Parameter structures containing important parameters for the model */
  SpinEOBParams seobParams;
  SpinEOBHCoeffs seobCoeffs;
  EOBParams eobParams;
  FacWaveformCoeffs hCoeffs;
  NewtonMultipolePrefixes prefixes;
  TidalEOBParams tidal1, tidal2;

  /* fStart is the start frequency of the 22 mode generation */
  REAL8 fStart;

  /* Initialize parameters */
  m1 = m1SI / LAL_MSUN_SI;
  m2 = m2SI / LAL_MSUN_SI;
  mTotal = m1 + m2;
  mTScaled = mTotal * LAL_MTSUN_SI;
  eta = m1 * m2 / (mTotal * mTotal);

  /* For v2 the upper bound is mass ratio 100 */
  if ((SpinAlignedEOBversion == 2 || SpinAlignedEOBversion == 4)
      && eta < 100. / 101. / 101.)
    {
      XLALPrintError
	("XLAL Error - %s: Mass ratio larger than 100!\nSEOBNRv2, SEOBNRv2_opt, SEOBNRv4, and SEOBNRv4_opt are only available for mass ratios up to 100.\n",
	 __func__);
      XLAL_ERROR (XLAL_EINVAL);
    }

  amp0 = mTotal * LAL_MRSUN_SI / r;

  #if debugOutput
  printf("Amp = %.16e\n", amp0);
  #endif

  fStart = fMin;

   /* If fMin is too high, then for safety in the initial conditions we integrate nonetheless from r=10M */
  if (pow (10., -3. / 2.)/(LAL_PI * mTScaled) < fMin) {
    fStart = pow (10., -3. / 2.)/(LAL_PI * mTScaled);
    //FP    XLAL_PRINT_WARNING ("Waveform will be generated from %f and stored from %f.",fStart,fMin);
  }

  /*
  if (pow (LAL_PI * fMin * mTScaled, -2. / 3.) < 10.0)
    {
      XLAL_PRINT_WARNING
	("Waveform generation may fail due to high starting frequency. The starting frequency corresponds to a small initial radius of %.2fM. We recommend a lower starting frequency that corresponds to an estimated starting radius > 10M.",
	 pow (LAL_PI * fMin * mTScaled, -2.0 / 3.0));
    }
  */

  /* TODO: Insert potentially necessary checks on the arguments */

  /* Calculate the time we will need to step back for ringdown */
  tStepBack = 100. * mTScaled;
  if (SpinAlignedEOBversion == 4)
    {
      tStepBack = 150. * mTScaled;
    }
  nStepBack = ceil (tStepBack / deltaT);


  /* Calculate the resample factor for attaching the ringdown */
  /* We want it to be a power of 2 */
  /* If deltaT > Mtot/50, reduce deltaT by the smallest power of two for which deltaT < Mtot/50 */
  resampEstimate = 50. * deltaT / mTScaled;
  resampFac = 1;
  //resampFac = 1 << (UINT4)ceil(log2(resampEstimate));

  if (resampEstimate > 1.)
    {
      resampPwr = (UINT4) ceil (log2 (resampEstimate));
      while (resampPwr--)
	{
	  resampFac *= 2u;
	}
    }


  /* Allocate the values vector to contain the initial conditions */
  /* Since we have aligned spins, we can use the 4-d vector as in the non-spin case */
  UINT4 num_elements_in_values_vector = 4;
  if (use_optimized_v2_or_v4)
    /* In v2opt/v4opt, we add an additional two elements to this vector, to store amplitude & phase.
       After ODE solves EOMs (which involves 4 coupled ODEs & thus 4 solution vectors),
       amp & phase are computed explicitly from sparse EOM solution, then interpolated in v2opt/v4opt. */
    num_elements_in_values_vector = 6;
  if (!(values = XLALCreateREAL8Vector (num_elements_in_values_vector)))
    {
      XLAL_ERROR (XLAL_ENOMEM);
    }
  memset (values->data, 0, values->length * sizeof (REAL8));

  /* Set up structures and calculate necessary PN parameters */
  /* Unlike the general case, we only need to calculate these once */
  memset (&seobParams, 0, sizeof (seobParams));
  memset (&seobCoeffs, 0, sizeof (seobCoeffs));
  memset (&eobParams, 0, sizeof (eobParams));
  memset (&hCoeffs, 0, sizeof (hCoeffs));
  memset (&prefixes, 0, sizeof (prefixes));

  /* Before calculating everything else, check sample freq is high enough */
  modefreqVec.length = 1;
  modefreqVec.data = &modeFreq;

  UINT4 mode_highest_freqL = 2;
  UINT4 mode_highest_freqM = 2;

  if (use_hm) {
    //RC: if we are using SEOBNRv4HM, the check for the Nyquist frequency
    //should be done for the 55 mode, the frequency of the RD scales with l
    mode_highest_freqL = 5;
    mode_highest_freqM = 5;
  }

  if (XLALSimIMREOBGenerateQNMFreqV2
      (&modefreqVec, m1, m2, spin1, spin2, mode_highest_freqL, mode_highest_freqM, 1,
       SpinAlignedEOBapproximant) == XLAL_FAILURE)
    {
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_EFUNC);
    }

    /* Check if fMin exceeds 95% the ringdown frequency; if so, don't generate a wf */
    REAL8 fRD = 0.95*creal (modeFreq)/(2.*LAL_PI);
//    UNUSED REAL8 fMerger = GetNRSpinPeakOmegaV4 (2, 2, eta, 0.5*(spin1z + spin2z) + 0.5*(spin1z - spin2z)*(m1 - m2)/(m1 + m2)/(1. - 2.*eta))/(2.*LAL_PI*mTScaled);
//    printf("fMin %.16e\n", fMin);
//    printf("fStart %.16e\n", fStart);
//    printf("fMerger %.16e\n", fMerger);
//    printf("fRD %.16e\n", fRD);
    if ( fMin > fRD ) {
        XLALPrintError
        ("XLAL Error - Starting frequency is above ringdown frequency!\n");
        XLALDestroyREAL8Vector (values);
        XLAL_ERROR (XLAL_EINVAL);
    }

  /* If Nyquist freq < 220 QNM freq for SEOBNRv1/2/4 and freq < 550 QNM for SEOBNRv4HM, exit */
  if (use_tidal == 0 && deltaT > LAL_PI / creal (modeFreq))
    {
      XLALPrintError
	("XLAL Error - %s: Ringdown frequency > Nyquist frequency!\nAt present this situation is not supported.\n",
	 __func__);
      XLALDestroyREAL8Vector (values);
     XLAL_ERROR (XLAL_EDOM);
    }

  if (!(sigmaStar = XLALCreateREAL8Vector (3)))
    {
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_ENOMEM);
    }

  if (!(sigmaKerr = XLALCreateREAL8Vector (3)))
    {
      XLALDestroyREAL8Vector (sigmaStar);
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_ENOMEM);
    }

// Rescale tidal polarizabilites by powers of mNS/M
// Rescale f-mode freqs by M/mNS
// No rescaling needed for spin-induced quadrupole parameters
  tidal1.mByM = m1SI / (m1SI + m2SI);
  tidal1.lambda2Tidal = lambda2Tidal1 * pow(tidal1.mByM,5);
  tidal1.omega02Tidal = omega02Tidal1 / tidal1.mByM;
  tidal1.lambda3Tidal = lambda3Tidal1 * pow(tidal1.mByM,7);
  tidal1.omega03Tidal = omega03Tidal1 / tidal1.mByM;
  tidal1.quadparam = quadparam1;

  tidal2.mByM = m2SI / (m1SI + m2SI);
  tidal2.lambda2Tidal = lambda2Tidal2 * pow(tidal2.mByM,5);
  tidal2.omega02Tidal = omega02Tidal2 / tidal2.mByM;
  tidal2.lambda3Tidal = lambda3Tidal2 * pow(tidal2.mByM,7);
  tidal2.omega03Tidal = omega03Tidal2 / tidal2.mByM;
  tidal2.quadparam = quadparam2;

  seobCoeffs.tidal1 = &tidal1;
  seobCoeffs.tidal2 = &tidal2;

  hCoeffs.tidal1 = &tidal1;
  hCoeffs.tidal2 = &tidal2;

  seobParams.deltaT = deltaT /( (m1 + m2) * LAL_MTSUN_SI );
  seobParams.alignedSpins = 1;
  seobParams.tortoise = 1;
  seobParams.sigmaStar = sigmaStar;
  seobParams.sigmaKerr = sigmaKerr;
  seobParams.seobCoeffs = &seobCoeffs;
  seobParams.eobParams = &eobParams;
  seobParams.nqcCoeffs = &nqcCoeffs;
  seobParams.use_hm = use_hm;
  eobParams.hCoeffs = &hCoeffs;
  eobParams.prefixes = &prefixes;

  eobParams.m1 = m1;
  eobParams.m2 = m2;
  eobParams.eta = eta;

  s1Vec.length = s2Vec.length = 3;
  s1VecOverMtMt.length = s2VecOverMtMt.length = 3;
  s1Vec.data = s1Data;
  s2Vec.data = s2Data;
  s1VecOverMtMt.data = s1DataNorm;
  s2VecOverMtMt.data = s2DataNorm;

    if ( use_tidal == 1 ) {
        REAL8 omegaMerger = XLALSimNSNSMergerFreq( &tidal1, &tidal2 );
        REAL8 rMerger = pow ( omegaMerger/2., -2./3. );
        if ( pow( fStart*LAL_PI*mTScaled, -2./3. ) <= 2.*rMerger ) {
            XLALPrintError
            ("XLAL Error - %s: fmin is too high for a tidal run, it should be at most %.16e Hz\n", __func__, pow (2.*rMerger, -3. / 2.)/(LAL_PI * mTScaled));
            XLAL_ERROR (XLAL_EINVAL);
        }
    }

  /* copy the spins into the appropriate vectors, and scale them by the mass */
  memcpy (s1Data, spin1, sizeof (s1Data));
  memcpy (s2Data, spin2, sizeof (s2Data));
  memcpy (s1DataNorm, spin1, sizeof (s1DataNorm));
  memcpy (s2DataNorm, spin2, sizeof (s2DataNorm));

  /* Calculate chiS and chiA */

  chiS = 0.5 * (spin1[2] + spin2[2]);
  chiA = 0.5 * (spin1[2] - spin2[2]);

  for (i = 0; i < 3; i++)
    {
      s1Data[i] *= m1 * m1;
      s2Data[i] *= m2 * m2;
    }
  for (i = 0; i < 3; i++)
    {
      s1DataNorm[i] = s1Data[i] / mTotal / mTotal;
      s2DataNorm[i] = s2Data[i] / mTotal / mTotal;
    }
  seobParams.s1Vec = &s1VecOverMtMt;
  seobParams.s2Vec = &s2VecOverMtMt;

  cartPosVec.length = cartMomVec.length = 3;
  cartPosVec.data = cartPosData;
  cartMomVec.data = cartMomData;
  memset (cartPosData, 0, sizeof (cartPosData));
  memset (cartMomData, 0, sizeof (cartMomData));

  /* Populate the initial structures */
  if (XLALSimIMRSpinEOBCalculateSigmaStar (sigmaStar, m1, m2, &s1Vec, &s2Vec)
      == XLAL_FAILURE)
    {
      XLALDestroyREAL8Vector (sigmaKerr);
      XLALDestroyREAL8Vector (sigmaStar);
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_EFUNC);
    }

  if (XLALSimIMRSpinEOBCalculateSigmaKerr (sigmaKerr, m1, m2, &s1Vec, &s2Vec)
      == XLAL_FAILURE)
    {
      XLALDestroyREAL8Vector (sigmaKerr);
      XLALDestroyREAL8Vector (sigmaStar);
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_EFUNC);
    }

  /* Calculate the value of a */
  /* XXX I am assuming that, since spins are aligned, it is okay to just use the z component XXX */
  /* TODO: Check this is actually the way it works in LAL */
  switch (SpinAlignedEOBversion)
    {
    case 1:
      tplspin = 0.0;
      break;
    case 2:
      tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;
      break;
    case 4:
      tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;
      break;
    default:
      XLALPrintError
	("XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1, v2, and v4 are available.\n",
	 __func__);
      if(sigmaKerr){
        XLALDestroyREAL8Vector (sigmaKerr);
      }
      if(sigmaStar){
        XLALDestroyREAL8Vector (sigmaStar);
      }
      if(values){
        XLALDestroyREAL8Vector (values);
      }
      XLAL_ERROR (XLAL_EINVAL);
      break;
    }
  /*for ( i = 0; i < 3; i++ )
     {
     a += sigmaKerr->data[i]*sigmaKerr->data[i];
     }
     a = sqrt( a ); */
  seobParams.a = a = sigmaKerr->data[2];
  seobParams.chi1 = spin1[2];
  seobParams.chi2 = spin2[2];

  /* Now compute the spinning H coefficients and store them in seobCoeffs */
    if (XLALSimIMRCalculateSpinEOBHCoeffs
      (&seobCoeffs, eta, a, SpinAlignedEOBversion) == XLAL_FAILURE)
    {
      XLALDestroyREAL8Vector (sigmaKerr);
      XLALDestroyREAL8Vector (sigmaStar);
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_EFUNC);
    }
  //RC: SEOBNRv4HM uses the same dynamics of SEOBNRv4, we momentarily put use_hm = 0 such that it uses the SEOBNRv4 coefficients to compute the flux LALSimIMRSpinEOBFactorizedFlux
  if(use_hm == 1){
    seobParams.use_hm = 0;
  }
  if (XLALSimIMREOBCalcSpinFacWaveformCoefficients
      (&hCoeffs, &seobParams, m1, m2, eta, tplspin, chiS, chiA,
       SpinAlignedEOBversion) == XLAL_FAILURE)
    {
      XLALDestroyREAL8Vector (sigmaKerr);
      XLALDestroyREAL8Vector (sigmaStar);
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_EFUNC);
    }
  if(use_hm == 1){
    seobParams.use_hm = 1;
  }

  if (XLALSimIMREOBComputeNewtonMultipolePrefixes
      (&prefixes, eobParams.m1, eobParams.m2) == XLAL_FAILURE)
    {
      XLALDestroyREAL8Vector (sigmaKerr);
      XLALDestroyREAL8Vector (sigmaStar);
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_EFUNC);
    }

  switch (SpinAlignedEOBversion)
    {
    case 1:
      if (XLALSimIMRGetEOBCalibratedSpinNQC (&nqcCoeffs, 2, 2, eta, a) ==
	  XLAL_FAILURE)
	{
    if(sigmaKerr){
      XLALDestroyREAL8Vector (sigmaKerr);
    }
    if(sigmaStar){
      XLALDestroyREAL8Vector (sigmaStar);
    }
    if(values){
      XLALDestroyREAL8Vector (values);
    }
	  XLAL_ERROR (XLAL_EFUNC);
	}
      break;
    case 2:
      if (XLALSimIMRGetEOBCalibratedSpinNQC3D
	  (&nqcCoeffs, 2, 2, m1, m2, a, chiA) == XLAL_FAILURE)
	{
    if(sigmaKerr){
      XLALDestroyREAL8Vector (sigmaKerr);
    }
    if(sigmaStar){
      XLALDestroyREAL8Vector (sigmaStar);
    }
    if(values){
      XLALDestroyREAL8Vector (values);
    }
	  XLAL_ERROR (XLAL_EFUNC);
	}
      break;
    case 4:
      nqcCoeffs.a1 = 0.;
      nqcCoeffs.a2 = 0.;
      nqcCoeffs.a3 = 0.;
      nqcCoeffs.a3S = 0.;
      nqcCoeffs.a4 = 0.;
      nqcCoeffs.a5 = 0.;
      nqcCoeffs.b1 = 0.;
      nqcCoeffs.b2 = 0.;
      nqcCoeffs.b3 = 0.;
      nqcCoeffs.b4 = 0.;
      break;
    default:
      XLALPrintError
	("XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1, v2, and v4 are available.\n",
	 __func__);
      if(sigmaKerr){
        XLALDestroyREAL8Vector (sigmaKerr);
      }
      if(sigmaStar){
        XLALDestroyREAL8Vector (sigmaStar);
      }
      if(values){
        XLALDestroyREAL8Vector (values);
      }
      XLAL_ERROR (XLAL_EINVAL);
      break;
    }

  /*
   * STEP 1) Solve for initial conditions
   */

  /* Set the initial conditions. For now we use the generic case */
  /* Can be simplified if spin-aligned initial conditions solver available. The cost of generic code is negligible though. */
  REAL8Vector *tmpValues = XLALCreateREAL8Vector (14);
  if (!tmpValues)
    {
      XLALDestroyREAL8Vector (sigmaKerr);
      XLALDestroyREAL8Vector (sigmaStar);
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_ENOMEM);
    }

  memset (tmpValues->data, 0, tmpValues->length * sizeof (REAL8));

  /* We set inc zero here to make it easier to go from Cartesian to spherical coords */
  /* No problem setting inc to zero in solving spin-aligned initial conditions. */
  /* inc is not zero in generating the final h+ and hx */
  if (XLALSimIMRSpinEOBInitialConditions
      (tmpValues, m1, m2, fStart, 0, s1Data, s2Data, &seobParams,
       use_optimized_v2_or_v4) == XLAL_FAILURE)
    {
      XLALDestroyREAL8Vector (tmpValues);
      XLALDestroyREAL8Vector (sigmaKerr);
      XLALDestroyREAL8Vector (sigmaStar);
      XLALDestroyREAL8Vector (values);
      XLAL_ERROR (XLAL_EFUNC);
    }

//  printf( "ICs = %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", tmpValues->data[0], tmpValues->data[1], tmpValues->data[2],
//     tmpValues->data[3], tmpValues->data[4], tmpValues->data[5], tmpValues->data[6], tmpValues->data[7], tmpValues->data[8],
//     tmpValues->data[9], tmpValues->data[10], tmpValues->data[11] );

  /* Taken from Andrea's code */
/*  memset( tmpValues->data, 0, tmpValues->length*sizeof(tmpValues->data[0]));*/
#if 0
  tmpValues->data[0] = 20;
  tmpValues->data[3] = -0.0001008684225106323/eta;
  tmpValues->data[4] = 1.16263606988612/eta / tmpValues->data[0];	// q=8 chi1=0.5
#endif

  /* Now convert to Spherical */
  /* The initial conditions code returns Cartesian components of four vectors x, p, S1 and S2,
   * in the special case that the binary starts on the x-axis and the two spins are aligned
   * with the orbital angular momentum along the z-axis.
   * Therefore, in spherical coordinates the initial conditions are
   * r = x; phi = 0.; pr = px; pphi = r * py.
   */
  values->data[0] = tmpValues->data[0];
  values->data[1] = 0.;
  values->data[2] = tmpValues->data[3];
  values->data[3] = tmpValues->data[0] * tmpValues->data[4];

  eobParams.rad = values->data[0];
  eobParams.omegaPeaked = 0;
  eobParams.omegaMerger = XLALSimNSNSMergerFreq( &tidal1, &tidal2 );
  eobParams.NyquistStop = 0;

  //fprintf( stderr, "Spherical initial conditions: %e %e %e %e\n", values->data[0], values->data[1], values->data[2], values->data[3] );

  /*
   * STEP 2) Evolve EOB trajectory until reaching the peak of orbital frequency
   */

  /* Now we have the initial conditions, we can initialize the adaptive integrator */
#if debugOutput
    printf("Begin integration\n");
#endif
    if ( use_tidal == 1 ) {
        if (!
            (integrator =
             XLALAdaptiveRungeKutta4Init (4, XLALSpinAlignedHcapDerivative,
                                          XLALSpinAlignedNSNSStopCondition,
                                          EPS_ABS, EPS_REL)))
        {
            XLALDestroyREAL8Vector (values);
            XLAL_ERROR (XLAL_EFUNC);
        }
    }
    else {
        if (use_optimized_v2_or_v4)
        {
            if (!
                (integrator =
                 XLALAdaptiveRungeKutta4InitEighthOrderInstead (4,
							  XLALSpinAlignedHcapDerivativeOptimized,
							  XLALEOBSpinAlignedStopCondition,
							  EPS_ABS, EPS_REL)))
            {
                XLALDestroyREAL8Vector (values);
                XLAL_ERROR (XLAL_EFUNC);
            }
        }
        else
        {
            if (!
                (integrator =
                 XLALAdaptiveRungeKutta4Init (4, XLALSpinAlignedHcapDerivative,
					XLALEOBSpinAlignedStopCondition,
					EPS_ABS, EPS_REL)))
            {
                XLALDestroyREAL8Vector (values);
                XLAL_ERROR (XLAL_EFUNC);
            }
        }
    }

  integrator->stopontestonly = 1;
  integrator->retries = 1;

  if (use_optimized_v2_or_v4)
    {
      /* BEGIN OPTIMIZED */
      retLen_fromOptStep2 =
	XLALAdaptiveRungeKutta4NoInterpolate (integrator, &seobParams,
					      values->data, 0.,
					      20. / mTScaled,
					      deltaT / mTScaled,
					      &dynamicstmp,2);/* Last parameter added when funcions were combined in LALAdaptiveRungeKuttaIntegrator.c*/
      if (retLen_fromOptStep2 == XLAL_FAILURE || !dynamicstmp)
        {
          XLAL_ERROR (XLAL_EFUNC);
        }
      retLen =
	SEOBNRv2OptimizedInterpolatorNoAmpPhase (dynamicstmp, 0.,
						 deltaT / mTScaled,
						 retLen_fromOptStep2,
						 &dynamics);
      /* END OPTIMIZED */
    }
  else
    {
      retLen =
	XLALAdaptiveRungeKutta4 (integrator, &seobParams, values->data, 0.,
				 20. / mTScaled, deltaT / mTScaled,
				 &dynamics);
    }
  if (retLen == XLAL_FAILURE || dynamics == NULL)
    {
      XLAL_ERROR (XLAL_EFUNC);
    }

  /* Set up pointers to the dynamics */
  // REAL8Vector tVec;
  // tVec.data = dynamics->data;
  // tVec.length = rVec.length = phiVec.length = prVec.length = pPhiVec.length = retLen;
  rVec.length = phiVec.length = prVec.length = pPhiVec.length = retLen;
  rVec.data = dynamics->data + retLen;
  phiVec.data = dynamics->data + 2 * retLen;
  prVec.data = dynamics->data + 3 * retLen;
  pPhiVec.data = dynamics->data + 4 * retLen;

  //printf( "We think we hit the peak at time %e\n", dynamics->data[retLen-1] );

  /* TODO : Insert high sampling rate / ringdown here */


  if (tStepBack > retLen * deltaT)
    {
      tStepBack = 0.5 * retLen * deltaT;	//YPnote: if 100M of step back > actual time of evolution, step back 50% of the later
      nStepBack = ceil (tStepBack / deltaT);
    }

  //SM
  // retLen is going to be reused for the length of the high-sampling dynamics
  INT4 retLen_out = retLen;
  //SM

  /*
   * STEP 3) Step back in time by tStepBack and volve EOB trajectory again
   *         using high sampling rate, stop at 0.3M out of the "EOB horizon".
   */

  /* Set up the high sample rate integration */
  hiSRndx = retLen - nStepBack;
  deltaTHigh = deltaT / (REAL8) resampFac;

#if debugOutput
  printf (
	   "Stepping back %d points - we expect %d points at high SR\n",
	   nStepBack, nStepBack * resampFac);
  printf (
	   "Commencing high SR integration... from %.16e %.16e %.16e %.16e %.16e\n",
	   (dynamics->data)[hiSRndx], rVec.data[hiSRndx],
	   phiVec.data[hiSRndx], prVec.data[hiSRndx], pPhiVec.data[hiSRndx]);
#endif

  values->data[0] = rVec.data[hiSRndx];
  values->data[1] = phiVec.data[hiSRndx];
  values->data[2] = prVec.data[hiSRndx];
  values->data[3] = pPhiVec.data[hiSRndx];
  eobParams.rad = values->data[0];
  eobParams.omegaPeaked = 0;
  eobParams.NyquistStop = 0;



  /* For HiSR evolution, we stop at a radius 0.3M from the deformed Kerr singularity,
   * or when any derivative of Hamiltonian becomes nan */
  integrator->stop = XLALSpinAlignedHiSRStopCondition;
  if (SpinAlignedEOBversion == 4)
    {
      integrator->stop = XLALSpinAlignedHiSRStopConditionV4;
    }
  if ( use_tidal == 1 ) {
      integrator->stop = XLALSpinAlignedNSNSStopCondition;
    }

  if (use_optimized_v2_or_v4)
    {
      /* BEGIN OPTIMIZED: */
      retLen_fromOptStep3 =
	XLALAdaptiveRungeKutta4NoInterpolate (integrator, &seobParams,
					      values->data, 0.,
					      20. / mTScaled,
					      deltaTHigh / mTScaled,
					      &dynamicsHitmp,2);/* Last parameter added when funcions were combined in LALAdaptiveRungeKuttaIntegrator.c*/
      if (retLen_fromOptStep3 == XLAL_FAILURE || !dynamicsHitmp)
        {
          XLAL_ERROR (XLAL_EFUNC);
        }
      retLen =
	SEOBNRv2OptimizedInterpolatorNoAmpPhase (dynamicsHitmp, 0.,
						 deltaTHigh / mTScaled,
						 retLen_fromOptStep3,
						 &dynamicsHi);
      /* END OPTIMIZED */
    }
  else
    {
      retLen =
	XLALAdaptiveRungeKutta4 (integrator, &seobParams, values->data, 0.,
				 20. / mTScaled, deltaTHigh / mTScaled,
				 &dynamicsHi);
    }
  if (retLen == XLAL_FAILURE || dynamicsHi == NULL)
    {
      if(tmpValues){
        XLALDestroyREAL8Vector (tmpValues);
      }
      if(sigmaKerr){
        XLALDestroyREAL8Vector (sigmaKerr);
      }
      if(sigmaStar){
        XLALDestroyREAL8Vector (sigmaStar);
      }
      if(values){
        XLALDestroyREAL8Vector (values);
      }
      XLAL_ERROR (XLAL_EFUNC);
    }

  //SM
  // retLen now means the length of the high-sampling dynamics
  // We also keep track of the starting time of the high-sampling dynamics
  INT4 retLenHi_out = retLen;
  REAL8 tstartHi = hiSRndx * deltaT / mTScaled;
  //SM

//  fprintf( stderr, "We got %d points at high SR\n", retLen );

  /* Set up pointers to the dynamics */
  rHi.length = phiHi.length = prHi.length = pPhiHi.length = timeHi.length =
    retLen;
  timeHi.data = dynamicsHi->data;
  rHi.data = dynamicsHi->data + retLen;
  phiHi.data = dynamicsHi->data + 2 * retLen;
  prHi.data = dynamicsHi->data + 3 * retLen;
  pPhiHi.data = dynamicsHi->data + 4 * retLen;


    /* Allocate the high sample rate vectors */
  sigReHi =
    XLALCreateREAL8Vector (retLen +
			   (UINT4) ceil (20 /
					 (cimag (modeFreq) * deltaTHigh)));
  sigImHi =
    XLALCreateREAL8Vector (retLen +
			   (UINT4) ceil (20 /
					 (cimag (modeFreq) * deltaTHigh)));
  omegaHi =
    XLALCreateREAL8Vector (retLen +
			   (UINT4) ceil (20 /
					 (cimag (modeFreq) * deltaTHigh)));
  ampNQC = XLALCreateREAL8Vector (retLen);
  phaseNQC = XLALCreateREAL8Vector (retLen);
  hamVHi = XLALCreateREAL8Vector (retLen);
  //RC: we save the Hamiltonian in a vector such that we can compute it only once
  //and use it for calculate all the modes

  if (!sigReHi || !sigImHi || !omegaHi || !ampNQC || !phaseNQC|| !hamVHi)
    {
      if(tmpValues){
        XLALDestroyREAL8Vector (tmpValues);
      }
      if(sigmaKerr){
        XLALDestroyREAL8Vector (sigmaKerr);
      }
      if(sigmaStar){
        XLALDestroyREAL8Vector (sigmaStar);
      }
      if(values){
        XLALDestroyREAL8Vector (values);
      }
      XLAL_ERROR (XLAL_ENOMEM);
    }

  memset (sigReHi->data, 0, sigReHi->length * sizeof (sigReHi->data[0]));
  memset (sigImHi->data, 0, sigImHi->length * sizeof (sigImHi->data[0]));

  /* Populate the high SR waveform */
  REAL8 omegaOld = 0.0;
  INT4 phaseCounter = 0;
  for (i = 0; i < retLen; i++)
    {
      values->data[0] = rHi.data[i];
      values->data[1] = phiHi.data[i];
      values->data[2] = prHi.data[i];
      values->data[3] = pPhiHi.data[i];

      if (use_optimized_v2_or_v4)
	{
	  /* OPTIMIZED: */
	  omega =
	    XLALSimIMRSpinAlignedEOBCalcOmegaOptimized (values->data,
							&seobParams);
	  /* END OPTIMIZED: */
	}
      else
	{
	  omega =
	    XLALSimIMRSpinAlignedEOBCalcOmega (values->data, &seobParams, STEP_SIZE);
	}

      if (omega < 1.0e-15)
	omega = 1.0e-9;		//YPnote: make sure omega>0 during very-late evolution when numerical errors are huge.
      omegaHi->data[i] = omega;	//YPnote: omega<0 is extremely rare and had only happenned after relevant time interval.
      v = cbrt (omega);

      /* Calculate the value of the Hamiltonian */
      cartPosVec.data[0] = values->data[0];
      cartMomVec.data[0] = values->data[2];
      cartMomVec.data[1] = values->data[3] / values->data[0];

      if (use_optimized_v2_or_v4)
	{
	  /* OPTIMIZED: */
	  hamVHi->data[i] =
	    XLALSimIMRSpinEOBHamiltonianOptimized (eta, &cartPosVec,
						   &cartMomVec,
						   &s1VecOverMtMt,
						   &s2VecOverMtMt, sigmaKerr,
						   sigmaStar,
						   seobParams.tortoise,
						   &seobCoeffs);
	  /* END OPTIMIZED: */
	}
      else
	{
	  hamVHi->data[i] =
	    XLALSimIMRSpinEOBHamiltonian (eta, &cartPosVec, &cartMomVec,
					  &s1VecOverMtMt, &s2VecOverMtMt,
					  sigmaKerr, sigmaStar,
					  seobParams.tortoise, &seobCoeffs);
	}

      if (omega <= omegaOld && !peakIdx)
	{
//      printf( "Have we got the peak? omegaOld = %.16e, omega = %.16e\n", omegaOld, omega );
	  peakIdx = i;
	}
      omegaOld = omega;
    }

  finalIdx = retLen - 1;
  if (!peakIdx)
    peakIdx = finalIdx;

//  printf( "We now think the peak is at %d\n", peakIdx );
#if debugOutput
    FILE *out = fopen ("saDynamicsHi.dat", "w");
    for (i = 0; i < retLen; i++)
    {
        fprintf (out, "%.16e %.16e %.16e %.16e %.16e %.16e\n", timeHi.data[i],
                 rHi.data[i], phiHi.data[i], prHi.data[i], pPhiHi.data[i], omegaHi->data[i]);
    }
    fclose (out);
#endif

  /*
   * STEP 4) Locate the peak of orbital frequency for NQC and QNM calculations
   */

  /* Stuff to find the actual peak time */
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;
  REAL8 omegaDeriv1;		//, omegaDeriv2;
  REAL8 time1, time2;
  REAL8 timePeak, timewavePeak = 0., omegaDerivMid;
  REAL8 sigAmpSqHi = 0., oldsigAmpSqHi = 0.;
  INT4 peakCount = 0;

  spline = gsl_spline_alloc (gsl_interp_cspline, retLen);
  acc = gsl_interp_accel_alloc ();

  time1 = dynamicsHi->data[peakIdx];

  gsl_spline_init (spline, dynamicsHi->data, omegaHi->data, retLen);
  omegaDeriv1 = gsl_spline_eval_deriv (spline, time1, acc);
  if (omegaDeriv1 > 0.)
    {
      time2 = dynamicsHi->data[peakIdx + 1];
      //omegaDeriv2 = gsl_spline_eval_deriv( spline, time2, acc );
    }
  else
    {
      //omegaDeriv2 = omegaDeriv1;
      time2 = time1;
      time1 = dynamicsHi->data[peakIdx - 1];
      peakIdx--;
      omegaDeriv1 = gsl_spline_eval_deriv (spline, time1, acc);
    }

  do
    {
      timePeak = (time1 + time2) / 2.;
      omegaDerivMid = gsl_spline_eval_deriv (spline, timePeak, acc);

      if (omegaDerivMid * omegaDeriv1 < 0.0)
	{
	  //omegaDeriv2 = omegaDerivMid;
	  time2 = timePeak;
	}
      else
	{
	  omegaDeriv1 = omegaDerivMid;
	  time1 = timePeak;
	}
    }
  while (time2 - time1 > 1.0e-5);

//    if (use_tidal == 1)
//        timePeak = dynamicsHi->data[retLen-1];


gsl_spline_free( spline );
gsl_interp_accel_free( acc );


  //XLALPrintInfo( "Estimation of the peak is now at time %.16e\n", timePeak );


  /*
   * STEP 5) Calculate NQC correction using hi-sampling data
   */

  /* Calculate nonspin and amplitude NQC coefficients from fits and interpolation table */
  /*switch ( SpinAlignedEOBversion )
     {
     case 1:
     if ( XLALSimIMRGetEOBCalibratedSpinNQC( &nqcCoeffs, 2, 2, eta, a ) == XLAL_FAILURE )
     {
     XLAL_ERROR( XLAL_EFUNC );
     }
     break;
     case 2:
     if ( XLALSimIMRGetEOBCalibratedSpinNQC3D( &nqcCoeffs, 2, 2, eta, a, chiA ) == XLAL_FAILURE )
     {
     XLAL_ERROR( XLAL_EFUNC );
     }
     break;
     default:
     XLALPrintError( "XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n", __func__);
     XLAL_ERROR( XLAL_EINVAL );
     break;
     } */

    /* Calculate the time of amplitude peak. Despite the name, this is in fact the shift in peak time from peak orb freq time */
    switch (SpinAlignedEOBversion)
    {
        case 1:
            timewavePeak = XLALSimIMREOBGetNRSpinPeakDeltaT (2, 2, eta, a);
            break;
        case 2:
            timewavePeak = XLALSimIMREOBGetNRSpinPeakDeltaTv2 (2, 2, m1, m2, spin1z, spin2z);	// David debug: we need to be using v2 for SpinAlignedEOBversion 2, right?
            break;
        case 4:
            timewavePeak =
            XLALSimIMREOBGetNRSpinPeakDeltaTv4 (2, 2, m1, m2, spin1z, spin2z);
            break;
        default:
            XLALPrintError
            ("XLAL Error - %s: Unknown SEOBNR version!\nAt present only v1 and v2 are available.\n",
             __func__);
             if(tmpValues){
               XLALDestroyREAL8Vector (tmpValues);
             }
             if(sigmaKerr){
               XLALDestroyREAL8Vector (sigmaKerr);
             }
             if(sigmaStar){
               XLALDestroyREAL8Vector (sigmaStar);
             }
             if(values){
               XLALDestroyREAL8Vector (values);
             }
             if(sigReHi){
               XLALDestroyREAL8Vector(sigReHi);
             }
             if(sigImHi){
               XLALDestroyREAL8Vector(sigImHi);
             }
             if(omegaHi){
               XLALDestroyREAL8Vector(omegaHi);
             }
             if(ampNQC){
               XLALDestroyREAL8Vector(ampNQC);
             }
             if(phaseNQC){
               XLALDestroyREAL8Vector(phaseNQC);
             }
             if(hamVHi){
               XLALDestroyREAL8Vector(hamVHi);
             }
            XLAL_ERROR (XLAL_EINVAL);
            break;
    }
    if (use_tidal == 1) {
        timewavePeak = 0.;
    }    /*                      */

    /* Evaluating the modes */

    /*The for over the modes should start here */
    // REAL8 nqcCoeffsMatrix[nModes][10];    //RC: Andrea coded the 2d array in this way. It works in C99 and indeed I don't get any errors when compiling, but I think this should be deprecated, it's better the gsl matrix
    gsl_matrix * nqcCoeffsMatrix = gsl_matrix_alloc (nModes, 10);


  if(use_hm == 1)
    {
    //RC for the SEOBNRv4HM model we need to compute again the coefficients for the modes because they are different from those using for computing the flux LALSimIMRSpinEOBFactorizedFlux.
    if (XLALSimIMREOBCalcSpinFacWaveformCoefficients
        (&hCoeffs, &seobParams, m1, m2, eta, tplspin, chiS, chiA,
         SpinAlignedEOBversion) == XLAL_FAILURE)
      {
        if(tmpValues){
          XLALDestroyREAL8Vector (tmpValues);
        }
        if(sigmaKerr){
          XLALDestroyREAL8Vector (sigmaKerr);
        }
        if(sigmaStar){
          XLALDestroyREAL8Vector (sigmaStar);
        }
        if(values){
          XLALDestroyREAL8Vector (values);
        }
        if(sigReHi){
          XLALDestroyREAL8Vector(sigReHi);
        }
        if(sigImHi){
          XLALDestroyREAL8Vector(sigImHi);
        }
        if(omegaHi){
          XLALDestroyREAL8Vector(omegaHi);
        }
        if(ampNQC){
          XLALDestroyREAL8Vector(ampNQC);
        }
        if(phaseNQC){
          XLALDestroyREAL8Vector(phaseNQC);
        }
        if(hamVHi){
          XLALDestroyREAL8Vector(hamVHi);
        }
        XLAL_ERROR (XLAL_EFUNC);
      }
      if((fabs(m1-m2) == 0) && (chiA == 0)){
        //RC: in the case of equal mass equal spins (where odd m modes are 0), the calibration parameter is 0
      }
      else{
      //RC: in addition to the PN coefficient, here we also need to evaluate a spinning pseudo-PN coefficient for the 21 and the 55 mode. These coefficients do not break
      //the odd mode's symmetry and for this reason they are 0 if chiS == chiA == 0. This of course hold also for non-spinning configurations.
        XLALSimIMREOBCalcCalibCoefficientHigherModes(&hCoeffs, &seobParams,2,1,&phiHi,&rHi,&prHi,
              omegaHi,hamVHi,&pPhiHi,timePeak -timewavePeak,m1,m2,chiS,chiA,deltaTHigh / mTScaled);
        XLALSimIMREOBCalcCalibCoefficientHigherModes(&hCoeffs, &seobParams,5,5,&phiHi,&rHi,&prHi,
                    omegaHi,hamVHi,&pPhiHi,timePeak -timewavePeak- 10,m1,m2,chiS,chiA,deltaTHigh / mTScaled);
        //RC: The - 10 here is because for the 55 mode the attachment is done at timePeak -timewavePeak- 10. timewavePeak is computed here outside
        //the for loop for the modes, for this reason I'm putting the -10 by hand, if it was inside the loop I could have used XLALSimIMREOBGetNRSpinPeakDeltaTv4 (5, 5, m1, m2, spin1z, spin2z)
      }
    }




    hLMAllHi = XLALCreateREAL8Vector((UINT4)2*sigReHi->length*nModes);
    memset(hLMAllHi->data, 0, hLMAllHi->length*sizeof (REAL8));
for ( UINT4 k = 0; k<nModes; k++) {
    modeL  = lmModes[k][0];
    modeM = lmModes[k][1];
#if debugOutput
    char filename2[sizeof "saModesXXHi.dat"];
    sprintf(filename2,"saModes%01d%01dHi.dat",modeL,modeM);
    out = fopen (filename2, "w");
#endif
    for(i=0; i<retLen; i++){
        values->data[0] = rHi.data[i];
        values->data[1] = phiHi.data[i];
        values->data[2] = prHi.data[i];
        values->data[3] = pPhiHi.data[i];
        v = cbrt (omegaHi->data[i]);
        if (XLALSimIMRSpinEOBGetSpinFactorizedWaveform
            (&hLM, values, v, hamVHi->data[i], modeL, modeM, &seobParams,
             use_optimized_v2_or_v4) == XLAL_FAILURE)
        {
            if(tmpValues){
              XLALDestroyREAL8Vector (tmpValues);
            }
            if(sigmaKerr){
              XLALDestroyREAL8Vector (sigmaKerr);
            }
            if(sigmaStar){
              XLALDestroyREAL8Vector (sigmaStar);
            }
            if(values){
              XLALDestroyREAL8Vector (values);
            }
            if(sigReHi){
              XLALDestroyREAL8Vector(sigReHi);
            }
            if(sigImHi){
              XLALDestroyREAL8Vector(sigImHi);
            }
            if(omegaHi){
              XLALDestroyREAL8Vector(omegaHi);
            }
            if(ampNQC){
              XLALDestroyREAL8Vector(ampNQC);
            }
            if(phaseNQC){
              XLALDestroyREAL8Vector(phaseNQC);
            }
            if(hamVHi){
              XLALDestroyREAL8Vector(hamVHi);
            }
            if(hLMAllHi){
              XLALDestroyREAL8Vector(hLMAllHi);
            }
            XLAL_ERROR (XLAL_EFUNC);
        }
#if debugOutput
        fprintf (out, "%.16e %.16e %.16e\n", timeHi.data[i],
                 creal (hLM) , cimag (hLM ) );
#endif
        ampNQC->data[i] = cabs (hLM);
        sigReHi->data[i] = (REAL8) (amp0 * creal (hLM));
        sigImHi->data[i] = (REAL8) (amp0 * cimag (hLM));
        phaseNQC->data[i] = carg (hLM) + phaseCounter * LAL_TWOPI;

        if (i && phaseNQC->data[i] > phaseNQC->data[i - 1])
        {
            phaseCounter--;
            phaseNQC->data[i] -= LAL_TWOPI;
        }

    }
#if debugOutput
    fclose (out);
#endif




  if (SpinAlignedEOBversion == 4)
    {
        if ( use_tidal == 1 && nqcFlag == 2 ) {
            nqcCoeffs.a1 = nqcCoeffsInput->data[0];
            nqcCoeffs.a2 = nqcCoeffsInput->data[1];
            nqcCoeffs.a3 = nqcCoeffsInput->data[2];
            nqcCoeffs.a3S = nqcCoeffsInput->data[3];
            nqcCoeffs.a4 = nqcCoeffsInput->data[4];
            nqcCoeffs.a5 = nqcCoeffsInput->data[5];
            nqcCoeffs.b1 = nqcCoeffsInput->data[6];
            nqcCoeffs.b2 = nqcCoeffsInput->data[7];
            nqcCoeffs.b3 = nqcCoeffsInput->data[8];
            nqcCoeffs.b4 = nqcCoeffsInput->data[9];
        }
        else {
          if(eta == 0.25 && chiA == 0 && ( (modeL == 2 && modeM == 1) || (modeL == 3 && modeM == 3)|| (modeL == 5 && modeM == 5) ) ) { //RC:Since the mode is 0 by symmetry, we don't need to compute the NQCs
            nqcCoeffs.a1 = 0.;
            nqcCoeffs.a2 = 0.;
            nqcCoeffs.a3 = 0.;
            nqcCoeffs.a3S = 0.;
            nqcCoeffs.a4 = 0.;
            nqcCoeffs.a5 = 0.;
            nqcCoeffs.b1 = 0.;
            nqcCoeffs.b2 = 0.;
            nqcCoeffs.b3 = 0.;
            nqcCoeffs.b4 = 0.;
          }
          else{
            if (XLALSimIMRSpinEOBCalculateNQCCoefficientsV4
                (ampNQC, phaseNQC, &rHi, &prHi, omegaHi, modeL, modeM, timePeak,
                 deltaTHigh / mTScaled, m1, m2, a, chiA, chiS, &nqcCoeffs,
                 SpinAlignedEOBversion) == XLAL_FAILURE)
            {
              if(tmpValues){
                XLALDestroyREAL8Vector (tmpValues);
              }
              if(sigmaKerr){
                XLALDestroyREAL8Vector (sigmaKerr);
              }
              if(sigmaStar){
                XLALDestroyREAL8Vector (sigmaStar);
              }
              if(values){
                XLALDestroyREAL8Vector (values);
              }
              if(sigReHi){
                XLALDestroyREAL8Vector(sigReHi);
              }
              if(sigImHi){
                XLALDestroyREAL8Vector(sigImHi);
              }
              if(omegaHi){
                XLALDestroyREAL8Vector(omegaHi);
              }
              if(ampNQC){
                XLALDestroyREAL8Vector(ampNQC);
              }
              if(phaseNQC){
                XLALDestroyREAL8Vector(phaseNQC);
              }
              if(hamVHi){
                XLALDestroyREAL8Vector(hamVHi);
              }
              if(hLMAllHi){
                XLALDestroyREAL8Vector(hLMAllHi);
              }
              XLAL_ERROR (XLAL_EFUNC);
            }
          }
        }
    }
    if ( SpinAlignedEOBversion == 4 && nqcFlag == 1 ) {
        nqcCoeffsInput->data[0] = nqcCoeffs.a1;
        nqcCoeffsInput->data[1] = nqcCoeffs.a2;
        nqcCoeffsInput->data[2] = nqcCoeffs.a3;
        nqcCoeffsInput->data[3] = nqcCoeffs.a3S;
        nqcCoeffsInput->data[4] = nqcCoeffs.a4;
        nqcCoeffsInput->data[5] = nqcCoeffs.a5;
        nqcCoeffsInput->data[6] = nqcCoeffs.b1;
        nqcCoeffsInput->data[7] = nqcCoeffs.b2;
        nqcCoeffsInput->data[8] = nqcCoeffs.b3;
        nqcCoeffsInput->data[9] = nqcCoeffs.b4;


        // FINISHED COMPUTING NQC. NOW MUST FREE ALLOCATED MEMORY!

        XLALDestroyREAL8Vector (tmpValues);
        XLALDestroyREAL8Vector (sigmaKerr);
        XLALDestroyREAL8Vector (sigmaStar);
        XLALDestroyREAL8Vector (values);
        XLALDestroyREAL8Vector (ampNQC);
        XLALDestroyREAL8Vector (phaseNQC);
        XLALDestroyREAL8Vector (sigReVec);
        XLALDestroyREAL8Vector (sigImVec);
        XLALAdaptiveRungeKuttaFree (integrator);
        XLALDestroyREAL8Array (dynamics);
        XLALDestroyREAL8Array (dynamicsHi);

        if (dynamicstmp)
          {
            XLALDestroyREAL8Array (dynamicstmp);
          }
        if (dynamicsHitmp)
          {
            XLALDestroyREAL8Array (dynamicsHitmp);
          }

        XLALDestroyREAL8Vector (sigReHi);
        XLALDestroyREAL8Vector (sigImHi);
        XLALDestroyREAL8Vector (omegaHi);

#if debugOutput
        printf
        ("Tidal point-mass NQC should not be 0 here: %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
         nqcCoeffs.a1, nqcCoeffs.a2, nqcCoeffs.a3, nqcCoeffs.a3S, nqcCoeffs.a4,
         nqcCoeffs.a5, nqcCoeffs.b1, nqcCoeffs.b2, nqcCoeffs.b3, nqcCoeffs.b4);
#endif
        return XLAL_SUCCESS;
    }
    /* Here we store the NQC coefficients for the different modes in some matrices */
    gsl_matrix_set(nqcCoeffsMatrix,k,0, nqcCoeffs.a1);
    gsl_matrix_set(nqcCoeffsMatrix,k,1, nqcCoeffs.a2);
    gsl_matrix_set(nqcCoeffsMatrix,k,2, nqcCoeffs.a3);
    gsl_matrix_set(nqcCoeffsMatrix,k,3, nqcCoeffs.a3S);
    gsl_matrix_set(nqcCoeffsMatrix,k,4, nqcCoeffs.a4);
    gsl_matrix_set(nqcCoeffsMatrix,k,5, nqcCoeffs.a5);
    gsl_matrix_set(nqcCoeffsMatrix,k,6, nqcCoeffs.b1);
    gsl_matrix_set(nqcCoeffsMatrix,k,7, nqcCoeffs.b2);
    gsl_matrix_set(nqcCoeffsMatrix,k,8, nqcCoeffs.b3);
    gsl_matrix_set(nqcCoeffsMatrix,k,9, nqcCoeffs.b4);
#if debugOutput
  printf
    ("(%d,%d)-mode NQC should not be 0 here: %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
     modeL, modeM, nqcCoeffs.a1, nqcCoeffs.a2, nqcCoeffs.a3, nqcCoeffs.a3S, nqcCoeffs.a4,
     nqcCoeffs.a5, nqcCoeffs.b1, nqcCoeffs.b2, nqcCoeffs.b3, nqcCoeffs.b4);
     FILE *NQC = NULL;
     sprintf(filename2,"NQC%01d%01dHi.dat",modeL,modeM);
     NQC = fopen (filename2, "w");
     fprintf(NQC, "%.16e \n%.16e \n%.16e \n%.16e \n%.16e\n", nqcCoeffs.a1, nqcCoeffs.a2, nqcCoeffs.a3,nqcCoeffs.b1, nqcCoeffs.b2);
     fclose(NQC);
#endif


/* Apply to the high sampled part */
#if debugOutput
    char filename[sizeof "saModesXXHiNQC.dat"];
    sprintf(filename,"saModes%01d%01dHiNQC.dat",modeL,modeM);
    out = fopen (filename, "w");
#endif
  for (i = 0; i < retLen; i++)
    {
      values->data[0] = rHi.data[i];
      values->data[1] = phiHi.data[i];
      values->data[2] = prHi.data[i];
      values->data[3] = pPhiHi.data[i];

      if (XLALSimIMREOBNonQCCorrection
	  (&hNQC, values, omegaHi->data[i], &nqcCoeffs) == XLAL_FAILURE)
	{
    if(tmpValues){
      XLALDestroyREAL8Vector (tmpValues);
    }
    if(sigmaKerr){
      XLALDestroyREAL8Vector (sigmaKerr);
    }
    if(sigmaStar){
      XLALDestroyREAL8Vector (sigmaStar);
    }
    if(values){
      XLALDestroyREAL8Vector (values);
    }
    if(sigReHi){
      XLALDestroyREAL8Vector(sigReHi);
    }
    if(sigImHi){
      XLALDestroyREAL8Vector(sigImHi);
    }
    if(omegaHi){
      XLALDestroyREAL8Vector(omegaHi);
    }
    if(ampNQC){
      XLALDestroyREAL8Vector(ampNQC);
    }
    if(phaseNQC){
      XLALDestroyREAL8Vector(phaseNQC);
    }
    if(hamVHi){
      XLALDestroyREAL8Vector(hamVHi);
    }
    if(hLMAllHi){
      XLALDestroyREAL8Vector(hLMAllHi);
    }
	  XLAL_ERROR (XLAL_EFUNC);
	}
      hLM = sigReHi->data[i];
      hLM += I * sigImHi->data[i];
      hLM *= hNQC;
        if ( (lambda2Tidal1 != 0. && omega02Tidal1 != 0.) || (lambda2Tidal2 != 0. && omega02Tidal2 != 0.) ) {
            if (XLALSimIMRSpinEOBWaveformTidal
                (&hT, values, cbrt(omegaHi->data[i]), 2, 2, &seobParams) )
            {
                XLAL_ERROR (XLAL_EFUNC);
            }
            hLM += amp0*hT;
        }

      sigReHi->data[i] = (REAL8) creal (hLM);
      sigImHi->data[i] = (REAL8) cimag (hLM);



        hLMAllHi->data[2*k*sigReHi->length + i] = sigReHi->data[i];
        hLMAllHi->data[(1+2*k)*sigReHi->length + i] = sigImHi->data[i];
#if debugOutput
        fprintf (out, "%.16e %.16e %.16e %.16e %.16e\n", timeHi.data[i],
                 creal (hLM) / amp0, cimag (hLM) / amp0,
                 hLMAllHi->data[2*k*sigReHi->length + i]/ amp0 ,hLMAllHi->data[(1+2*k)*sigReHi->length + i]/ amp0
                 );
#endif

      if (SpinAlignedEOBversion==1 || SpinAlignedEOBversion==2) {
          sigAmpSqHi = creal (hLM) * creal (hLM) + cimag (hLM) * cimag (hLM);
          if (sigAmpSqHi < oldsigAmpSqHi && peakCount == 0 && (i - 1) * deltaTHigh / mTScaled < timePeak-timewavePeak)
          {
              timewavePeak = (i - 1) * deltaTHigh / mTScaled;
              peakCount += 1;
          }
          oldsigAmpSqHi = sigAmpSqHi;
      }
    }
#if debugOutput
  fclose (out);
  printf ("NQCs entering hNQC: %f, %f, %f, %f, %f, %f\n", nqcCoeffs.a1,
	  nqcCoeffs.a2, nqcCoeffs.a3, nqcCoeffs.a3S, nqcCoeffs.a4,
	  nqcCoeffs.a5);
  printf ("NQCs entering hNQC: %f, %f, %f, %f\n", nqcCoeffs.b1, nqcCoeffs.b2,
	  nqcCoeffs.b3, nqcCoeffs.b4);
  printf ("Stas, again: timePeak = %.16f, timewavePeak = %.16f \n", timePeak,
	  timewavePeak);
#endif

}




/* REMOVE THIS */


    if (timewavePeak < 1.0e-16 || peakCount == 0)
    {
      //printf("YP::warning: could not locate mode peak, use calibrated time shift of amplitude peak instead.\n");
      /* NOTE: instead of looking for the actual peak, use the calibrated value,    */
      /*       ignoring the error in using interpolated NQC instead of iterated NQC */
      timewavePeak = timePeak - timewavePeak;
    }

  /*
   * STEP 6) Calculate QNM excitation coefficients using hi-sampling data
   */

  /*out = fopen( "saInspWaveHi.dat", "w" );
     for ( i = 0; i < retLen; i++ )
     {
     fprintf( out, "%.16e %.16e %.16e\n", timeHi.data[i], sigReHi->data[i], sigImHi->data[i] );
     }
     fclose( out ); */

  /* Attach the ringdown at the time of amplitude peak */
  REAL8 combSize = 7.5;		/* Eq. 34 */
  REAL8 chi =
    (spin1[2] + spin2[2]) / 2. +
    ((spin1[2] - spin2[2]) / 2.) * ((m1 - m2) / (m1 + m2)) / (1. - 2. * eta);

  /* Modify the combsize for SEOBNRv2 */
  /* If chi1=chi2=0, comb = 11. if chi < 0.8, comb = 12. if chi >= 0.8, comb =
   * 13.5 */
   //RC: the combsize is not used for the RD attachment for v4, I don't understand why the code enters inside this loop even for v4...
  if (SpinAlignedEOBversion == 2 || SpinAlignedEOBversion == 4)
    {
      combSize = (spin1[2] == 0.
		  && spin2[2] == 0.) ? 11. : ((eta > 10. / 121.
					       && chi >= 0.8) ? 8.5 : 12.);
      if ((eta > 30. / 31. / 31. && eta <= 10. / 121. && chi >= 0.8)
	  || (eta <= 30. / 31. / 31. && chi >= 0.8 && chi < 0.9))
	combSize = 13.5;
    }

  REAL8 timeshiftPeak;
  timeshiftPeak = timePeak - timewavePeak;
  if (SpinAlignedEOBversion == 2 || SpinAlignedEOBversion == 4)
    {
      timeshiftPeak =
	(timePeak - timewavePeak) > 0. ? (timePeak - timewavePeak) : 0.;
    }

  /*printf("YP::timePeak and timewavePeak: %.16e and %.16e\n",timePeak,timewavePeak);
     printf("YP::timeshiftPeak and combSize: %.16e and %.16e\n",timeshiftPeak,combSize);
     printf("PK::chi and SpinAlignedEOBversion: %.16e and %u\n\n", chi,SpinAlignedEOBversion); */

  REAL8Vector *rdMatchPoint = XLALCreateREAL8Vector (4);
  if (!rdMatchPoint)
    {
      if(tmpValues){
        XLALDestroyREAL8Vector (tmpValues);
      }
      if(sigmaKerr){
        XLALDestroyREAL8Vector (sigmaKerr);
      }
      if(sigmaStar){
        XLALDestroyREAL8Vector (sigmaStar);
      }
      if(values){
        XLALDestroyREAL8Vector (values);
      }
      if(sigReHi){
        XLALDestroyREAL8Vector(sigReHi);
      }
      if(sigImHi){
        XLALDestroyREAL8Vector(sigImHi);
      }
      if(omegaHi){
        XLALDestroyREAL8Vector(omegaHi);
      }
      if(ampNQC){
        XLALDestroyREAL8Vector(ampNQC);
      }
      if(phaseNQC){
        XLALDestroyREAL8Vector(phaseNQC);
      }
      if(hamVHi){
        XLALDestroyREAL8Vector(hamVHi);
      }
      if(hLMAllHi){
        XLALDestroyREAL8Vector(hLMAllHi);
      }
      XLAL_ERROR (XLAL_ENOMEM);
    }

  if (combSize > timePeak - timeshiftPeak)
    {
      XLALPrintError ("The comb size looks to be too big!!!\n");
    }

    REAL8Vector *OmVec = NULL;
    if (use_tidal == 1) {
        timeshiftPeak = 0.;
        UINT4 indAmax = 0;
        REAL8 Anew, Aval = sqrt( sigReHi->data[0]*sigReHi->data[0] +sigImHi->data[0]*sigImHi->data[0] );
        for (i =  1; i < (INT4) timeHi.length; i++) {
            Anew = sqrt( sigReHi->data[i]*sigReHi->data[i] +sigImHi->data[i]*sigImHi->data[i]);
            if ( Aval > Anew ) {
                indAmax = i-1;
                break;
            } else {
                indAmax = i;
                Aval = Anew;
            }
        }
#if debugOutput
        printf("indAmax %d\n",indAmax);
#endif

        if ( (INT4)indAmax== (INT4) timeHi.length-1 ) {
            INT4 peakDet = 0;
            REAL8 dAnew, dAval;
            REAL8 Avalp = sqrt( sigReHi->data[1]*sigReHi->data[1] +sigImHi->data[1]*sigImHi->data[1] );
            REAL8 Avalm = sqrt( sigReHi->data[0]*sigReHi->data[0] +sigImHi->data[0]*sigImHi->data[0] );
            dAval= Avalp - Avalm;
            INT4 iSkip = (INT4) 1./(deltaTHigh / mTScaled);
            for (i=(INT4) timeHi.length/2; i<(INT4) timeHi.length; i=i+iSkip) {
                Avalm = sqrt( sigReHi->data[i-1]*sigReHi->data[i-1] +sigImHi->data[i-1]*sigImHi->data[i-1]);
                Avalp = sqrt( sigReHi->data[i]*sigReHi->data[i] +sigImHi->data[i]*sigImHi->data[i]);
                dAnew = Avalp - Avalm;
#if debugOutput
                printf("%.16e %.16e %.16e\n", i*deltaTHigh / mTScaled, dAnew, dAval);
#endif
                if ( peakDet==0) {
                    if ( dAnew<dAval) peakDet++;
                    dAval = dAnew;
                }
                else {
                    if ( dAnew>dAval ) {
                        indAmax = i-1;
                        break;
                    }
                    else {
                        dAval = dAnew;
                    }
                }
            }
        }
#if debugOutput
        printf("indAmax %d\n",indAmax);
#endif
        UINT4 indOmax = 0;
        REAL8 dt = timeHi.data[1] - timeHi.data[0];
        REAL8 re = sigReHi->data[0], im = sigImHi->data[0];
        REAL8 red = ( sigReHi->data[1] - sigReHi->data[0] ) / dt, imd = ( sigImHi->data[1] - sigImHi->data[0] ) / dt;
        REAL8 Onew, Oval = - (imd*re - red*im) / (re*re + im*im);

        OmVec = XLALCreateREAL8Vector (sigReHi->length);
        memset (OmVec->data, 0, OmVec->length * sizeof (REAL8));

        for (i = 1; i < (INT4) timeHi.length-1; i++) {
            re = sigReHi->data[i];
            im = sigImHi->data[i];
            red = ( sigReHi->data[i+1] - sigReHi->data[i] ) / dt;
            imd = ( sigImHi->data[i+1] - sigImHi->data[i] ) / dt;
            Onew = - (imd*re - red*im) / (re*re + im*im);
            OmVec->data[i] = Onew;
            if ( Onew >= Oval ) {
                indOmax = i;
                Oval = Onew;
            }
        }

        if ( indAmax>timeHi.length/2 && timeHi.data[indAmax] <= timeHi.data[indOmax] ) {
            timePeak = timeHi.data[indAmax] ;
        }
        else {
            timePeak = timeHi.data[indOmax];
        }
        if ( eobParams.NyquistStop ==1 ) timePeak = timeHi.data[timeHi.length - 1];
    }

    /* Having located the peak of orbital frequency, we set time and phase of coalescence */
    XLALGPSAdd (&tc, -mTScaled * (dynamics->data[hiSRndx] + timePeak));



  rdMatchPoint->data[0] =
    combSize <
    timePeak - timeshiftPeak ? timePeak - timeshiftPeak - combSize : 0;
  //rdMatchPoint->data[0] = timePeak + 2.0*deltaTHigh;
  rdMatchPoint->data[1] = timePeak - timeshiftPeak;
  rdMatchPoint->data[2] = dynamicsHi->data[finalIdx];
#if debugOutput
  printf ("YP::comb range: %f, %f\n", rdMatchPoint->data[0],
	  rdMatchPoint->data[1]);
#endif
  rdMatchPoint->data[0] -=
    fmod (rdMatchPoint->data[0], deltaTHigh / mTScaled);
  rdMatchPoint->data[1] -=
    fmod (rdMatchPoint->data[1], deltaTHigh / mTScaled);
#if debugOutput
    printf("tattach = %.16f\n", rdMatchPoint->data[1]);
#endif
UINT4 indAmpMax = 0; //RC this variable is only used for SEOBNRv4HM. It stores the index of the attaching point of the 22 mode

for ( UINT4 k = 0; k<nModes; k++) {

    modeL  = lmModes[k][0];
    modeM = lmModes[k][1];

    //RC: Here I modify the attachment point for the 55 mode, which is tpeak22 -10M and is passed with the variable rdMatchPoint->data[3]
    if((modeL == 5)&&(modeM==5)){
      rdMatchPoint->data[3] = timePeak - timeshiftPeak-10;
      rdMatchPoint->data[3] -= fmod (rdMatchPoint->data[3], deltaTHigh / mTScaled);
    }

    for ( i=0; i<retLen; i++ ) {
        sigReHi->data[i] = hLMAllHi->data[i + 2*k*sigReHi->length];
        sigImHi->data[i] = hLMAllHi->data[i + (2*k+1)*sigReHi->length];
    }
    for ( i=retLen; i<(INT4)sigReHi->length; i++ ) {
        sigReHi->data[i] = 0.;
        sigImHi->data[i] = 0.;
    }

    if ( use_tidal == 1 ) {
        INT4 kount;
        REAL8 dtGeom = deltaTHigh / mTScaled;
        REAL8Vector *ampl, *phtmp;
        ampl = XLALCreateREAL8Vector (sigReHi->length);
        memset (ampl->data, 0, ampl->length * sizeof (REAL8));
        phtmp = XLALCreateREAL8Vector (sigReHi->length);
        memset (phtmp->data, 0,phtmp->length * sizeof (REAL8));
        INT4 iEnd= (INT4)rdMatchPoint->data[1]/dtGeom;
        UINT4 iM = iEnd - 1000;
        REAL8 omega0 = OmVec->data[iEnd];
        REAL8 omegaM = OmVec->data[iM];
        REAL8 omegaMdot = (OmVec->data[iM]-OmVec->data[iM-10])/(10*dtGeom);
        REAL8 tau = 0.5*LAL_PI/omega0;

        for (kount=0; kount<iEnd; kount++){
            ampl->data[kount] = sqrt(sigReHi->data[kount]/amp0*sigReHi->data[kount]/amp0 + sigImHi->data[kount]/amp0*sigImHi->data[kount]/amp0)/(1. + exp((kount*dtGeom-rdMatchPoint->data[1]-15)/tau));
            phtmp->data[kount] = carg(sigReHi->data[kount] + I * sigImHi->data[kount]);
        }
        REAL8 amplEnd = sqrt(sigReHi->data[iEnd]/amp0*sigReHi->data[iEnd]/amp0 + sigImHi->data[iEnd]/amp0*sigImHi->data[iEnd]/amp0);
        REAL8 amplEndM1 = sqrt(sigReHi->data[iEnd-1]/amp0*sigReHi->data[iEnd-1]/amp0 + sigImHi->data[iEnd-1]/amp0*sigImHi->data[iEnd-1]/amp0);
        REAL8 ampldotEnd = (amplEnd-amplEndM1)/dtGeom;
        if ( ampldotEnd < 0. ) ampldotEnd = 0.;
        for (kount=iEnd;kount<(INT4)(sigReHi->length);kount++){
            ampl->data[kount] = (amplEnd + ampldotEnd*(kount - iEnd)*dtGeom)/(1. + exp((kount*dtGeom-rdMatchPoint->data[1]-15)/tau));
        }
        for (kount=(INT4)iM;kount<(INT4)(sigReHi->length);kount++){
            phtmp->data[kount] = phtmp->data[iM] - ( omega0*(kount-(INT4)iM)*dtGeom + (omega0-omegaM)*(omega0-omegaM)/omegaMdot*(exp(-(kount-(INT4)iM)*dtGeom/(omega0-omegaM)*omegaMdot) - 1.) );
        }
#if debugOutput
        FILE *testout = fopen ("test.dat", "w");
        for (kount=0;kount<(INT4)(sigReHi->length);kount++){
            fprintf(testout, "%.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", kount*dtGeom,ampl->data[kount],phtmp->data[kount],sigReHi->data[kount],sigImHi->data[kount],sqrt(sigReHi->data[kount]/amp0*sigReHi->data[kount]/amp0 + sigImHi->data[kount]/amp0*sigImHi->data[kount]/amp0),carg(sigReHi->data[kount] + I * sigImHi->data[kount]));
        }
        fclose(testout);
#endif
        for (kount=0;kount<(INT4)(sigReHi->length);kount++){
            sigReHi->data[kount] = amp0*ampl->data[kount]*cos(phtmp->data[kount]);
            sigImHi->data[kount] = amp0*ampl->data[kount]*sin(phtmp->data[kount]);
        }
        XLALDestroyREAL8Vector(ampl);
        XLALDestroyREAL8Vector(phtmp);
                            /*
        if (XLALSimIMREOBTaper (sigReHi, sigImHi, 2, 2,
                                            deltaTHigh, m1, m2, spin1[0],
                                            spin1[1], spin1[2], spin2[0],
                                            spin2[1], spin2[2], &timeHi,
                                            rdMatchPoint,
                                            SpinAlignedEOBapproximant) ==
            XLAL_FAILURE)
        {
            XLAL_ERROR (XLAL_EFUNC);
        }
                             */
    }
    else {
        if (SpinAlignedEOBversion == 1 || SpinAlignedEOBversion == 2)
        {
            if (XLALSimIMREOBHybridAttachRingdown (sigReHi, sigImHi, 2, 2,
					     deltaTHigh, m1, m2, spin1[0],
					     spin1[1], spin1[2], spin2[0],
					     spin2[1], spin2[2], &timeHi,
					     rdMatchPoint,
					     SpinAlignedEOBapproximant) ==
                XLAL_FAILURE)
            {
                XLAL_ERROR (XLAL_EFUNC);
            }
        }
        else if (SpinAlignedEOBversion == 4)
        {

            if (XLALSimIMREOBAttachFitRingdown (sigReHi, sigImHi, modeL, modeM,
					  deltaTHigh, m1, m2, spin1[0],
					  spin1[1], spin1[2], spin2[0],
					  spin2[1], spin2[2], &timeHi,
					  rdMatchPoint,
					  SpinAlignedEOBapproximant, &indAmpMax) ==
                XLAL_FAILURE)
            {
              if(tmpValues){
                XLALDestroyREAL8Vector (tmpValues);
              }
              if(sigmaKerr){
                XLALDestroyREAL8Vector (sigmaKerr);
              }
              if(sigmaStar){
                XLALDestroyREAL8Vector (sigmaStar);
              }
              if(values){
                XLALDestroyREAL8Vector (values);
              }
              if(sigReHi){
                XLALDestroyREAL8Vector(sigReHi);
              }
              if(sigImHi){
                XLALDestroyREAL8Vector(sigImHi);
              }
              if(omegaHi){
                XLALDestroyREAL8Vector(omegaHi);
              }
              if(ampNQC){
                XLALDestroyREAL8Vector(ampNQC);
              }
              if(phaseNQC){
                XLALDestroyREAL8Vector(phaseNQC);
              }
              if(hamVHi){
                XLALDestroyREAL8Vector(hamVHi);
              }
              if(hLMAllHi){
                XLALDestroyREAL8Vector(hLMAllHi);
              }
              XLAL_ERROR (XLAL_EFUNC);
            }

        }

    }
    for ( i=0; i<(INT4)sigReHi->length; i++) {
        hLMAllHi->data[2*k*sigReHi->length + i] = sigReHi->data[i];
        hLMAllHi->data[(1+2*k)*sigReHi->length + i] = sigImHi->data[i];
    }
#if debugOutput
    char filename[sizeof "saModesXXHiIMR.dat"];
    sprintf(filename,"saModes%01d%01dHiIMR.dat",modeL,modeM);
    out = fopen (filename, "w");
    for ( i=0; i<(INT4)sigReHi->length; i++) {
        fprintf (out, "%.16e %.16e %.16e\n", i*deltaTHigh / mTScaled,hLMAllHi->data[2*k*sigReHi->length + i]/amp0,hLMAllHi->data[(1+2*k)*sigReHi->length + i]/amp0);
    }
    fclose(out);
#endif


}


  /*
   * STEP 7) Generate full inspiral waveform using desired sampling frequency
   */

  if (use_optimized_v2_or_v4)
    {
      // maybe dynamicstmp and dynamicsHitmp should be called "intermediateDynamics(Hi)" now since they aren't so temporary anymore?
      GenerateAmpPhaseFromEOMSoln (retLen_fromOptStep2, dynamicstmp->data,
				   &seobParams);
      /*
       * We used dynamics and dynamicsHi to store solution to equations of motion.
       *   The solution was needed to find, e.g., the time at peak freq (STEP 4).
       *   At this point, the solution to the EOMs is no longer needed, so we
       *   now repurpose dynamics and dynamicsHi to store amplitude and phase
       *   information only.
`      * This is the most efficient solution, as it frees up unused memory.
       */
      XLALDestroyREAL8Array (dynamics);
      XLALDestroyREAL8Array (dynamicsHi);
      dynamics = NULL;
      dynamicsHi = NULL;
      retLen =
	SEOBNRv2OptimizedInterpolatorOnlyAmpPhase (dynamicstmp, 0.,
						   deltaT / mTScaled,
						   retLen_fromOptStep2,
						   &dynamics);

      ampVec.length = phaseVec.length = retLen;
      ampVec.data = dynamics->data + 5 * retLen;
      phaseVec.data = dynamics->data + 6 * retLen;

      GenerateAmpPhaseFromEOMSoln (retLen_fromOptStep3, dynamicsHitmp->data,
				   &seobParams);
      retLen =
	SEOBNRv2OptimizedInterpolatorOnlyAmpPhase (dynamicsHitmp, 0.,
						   deltaTHigh / mTScaled,
						   retLen_fromOptStep3,
						   &dynamicsHi);
    }




  /* Now create vectors at the correct sample rate, and compile the complete waveform */
  sigReVec =
    XLALCreateREAL8Vector (rVec.length + ceil (sigReHi->length / resampFac));
  sigImVec = XLALCreateREAL8Vector (sigReVec->length);

  memset (sigReVec->data, 0, sigReVec->length * sizeof (REAL8));
  memset (sigImVec->data, 0, sigImVec->length * sizeof (REAL8));

  hLMAll = XLALCreateREAL8Vector((UINT4)2*sigReVec->length*nModes);
  memset(hLMAll->data, 0, hLMAll->length*sizeof (REAL8));


  /* Generate full inspiral waveform using desired sampling frequency */
  if (use_optimized_v2_or_v4)
    {
      for (i = 0; i < (INT4) rVec.length; i++)
	{

	  hLM = ampVec.data[i] * cexp (I * phaseVec.data[i]);

	  sigReVec->data[i] = amp0 * creal (hLM);
	  sigImVec->data[i] = amp0 * cimag (hLM);
    hLMAll->data[i] = sigReVec->data[i];
    hLMAll->data[sigReVec->length + i] = sigImVec->data[i];
	}
    }
  else
    {
#if debugOutput
        out = fopen ("saDynamics.dat", "w");
#endif
        hamV = XLALCreateREAL8Vector(rVec.length);
        memset(hamV->data, 0., hamV->length*sizeof(REAL8));
        omegaVec = XLALCreateREAL8Vector(rVec.length);
        memset(omegaVec->data, 0., omegaVec->length*sizeof(REAL8));
        if (!omegaVec|| !hamV)
        {
          if(tmpValues){
            XLALDestroyREAL8Vector (tmpValues);
          }
          if(sigmaKerr){
            XLALDestroyREAL8Vector (sigmaKerr);
          }
          if(sigmaStar){
            XLALDestroyREAL8Vector (sigmaStar);
          }
          if(values){
            XLALDestroyREAL8Vector (values);
          }
          if(sigReHi){
            XLALDestroyREAL8Vector(sigReHi);
          }
          if(sigImHi){
            XLALDestroyREAL8Vector(sigImHi);
          }
          if(omegaHi){
            XLALDestroyREAL8Vector(omegaHi);
          }
          if(ampNQC){
            XLALDestroyREAL8Vector(ampNQC);
          }
          if(phaseNQC){
            XLALDestroyREAL8Vector(phaseNQC);
          }
          if(hamVHi){
            XLALDestroyREAL8Vector(hamVHi);
          }
          if(hLMAllHi){
            XLALDestroyREAL8Vector(hLMAllHi);
          }
          if(sigReVec){
            XLALDestroyREAL8Vector(sigReVec);
          }
          if(sigImVec){
            XLALDestroyREAL8Vector(sigImVec);
          }
          if(hLMAll){
            XLALDestroyREAL8Vector(hLMAll);
          }
          XLAL_ERROR (XLAL_ENOMEM);
        }
      for (i = 0; i < (INT4) rVec.length; i++)
	{
	  values->data[0] = rVec.data[i];
	  values->data[1] = phiVec.data[i];
	  values->data[2] = prVec.data[i];
	  values->data[3] = pPhiVec.data[i];

	  /* Do not need to add an if(use_optimized_v2_or_v4), since this is strictly unoptimized code (see if(use_optimized_v2_or_v4) above) */
	  omegaVec->data[i] =
        XLALSimIMRSpinAlignedEOBCalcOmega (values->data, &seobParams, STEP_SIZE);
#if debugOutput
        fprintf (out, "%.16e %.16e %.16e %.16e %.16e %.16e\n", dynamics->data[i],
                 rVec.data[i], phiVec.data[i], prVec.data[i], pPhiVec.data[i],omegaVec->data[i]);
#endif
	  /* Calculate the value of the Hamiltonian */
	  cartPosVec.data[0] = values->data[0];
	  cartMomVec.data[0] = values->data[2];
	  cartMomVec.data[1] = values->data[3] / values->data[0];

	  hamV->data[i] =
	    XLALSimIMRSpinEOBHamiltonian (eta, &cartPosVec, &cartMomVec,
					  &s1VecOverMtMt, &s2VecOverMtMt,
					  sigmaKerr, sigmaStar,
					  seobParams.tortoise, &seobCoeffs);
    }



    for ( UINT4 k = 0; k<nModes; k++) {
        modeL  = lmModes[k][0];
        modeM = lmModes[k][1];

        nqcCoeffs.a1 = gsl_matrix_get(nqcCoeffsMatrix,k,0);
        nqcCoeffs.a2 = gsl_matrix_get(nqcCoeffsMatrix,k,1);
        nqcCoeffs.a3 = gsl_matrix_get(nqcCoeffsMatrix,k,2);
        nqcCoeffs.a3S = gsl_matrix_get(nqcCoeffsMatrix,k,3);
        nqcCoeffs.a4 = gsl_matrix_get(nqcCoeffsMatrix,k,4);
        nqcCoeffs.a5 = gsl_matrix_get(nqcCoeffsMatrix,k,5);
        nqcCoeffs.b1 = gsl_matrix_get(nqcCoeffsMatrix,k,6);
        nqcCoeffs.b2 = gsl_matrix_get(nqcCoeffsMatrix,k,7);
        nqcCoeffs.b3 = gsl_matrix_get(nqcCoeffsMatrix,k,8);
        nqcCoeffs.b4 = gsl_matrix_get(nqcCoeffsMatrix,k,9);
#if debugOutput
        printf
        ("(%d,%d)-mode NQC should not be 0 here: %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
         modeL, modeM, nqcCoeffs.a1, nqcCoeffs.a2, nqcCoeffs.a3, nqcCoeffs.a3S, nqcCoeffs.a4,
         nqcCoeffs.a5, nqcCoeffs.b1, nqcCoeffs.b2, nqcCoeffs.b3, nqcCoeffs.b4);
#endif
        for (i = 0; i < (INT4) rVec.length; i++)
        {
            values->data[0] = rVec.data[i];
            values->data[1] = phiVec.data[i];
            values->data[2] = prVec.data[i];
            values->data[3] = pPhiVec.data[i];
            if (XLALSimIMRSpinEOBGetSpinFactorizedWaveform
                (&hLM, values,  cbrt (omegaVec->data[i]), hamV->data[i], modeL, modeM, &seobParams, 0 /*use_optimized_v2_or_v4 */ ) == XLAL_FAILURE)
            {
              if(tmpValues){
                XLALDestroyREAL8Vector (tmpValues);
              }
              if(sigmaKerr){
                XLALDestroyREAL8Vector (sigmaKerr);
              }
              if(sigmaStar){
                XLALDestroyREAL8Vector (sigmaStar);
              }
              if(values){
                XLALDestroyREAL8Vector (values);
              }
              if(sigReHi){
                XLALDestroyREAL8Vector(sigReHi);
              }
              if(sigImHi){
                XLALDestroyREAL8Vector(sigImHi);
              }
              if(omegaHi){
                XLALDestroyREAL8Vector(omegaHi);
              }
              if(ampNQC){
                XLALDestroyREAL8Vector(ampNQC);
              }
              if(phaseNQC){
                XLALDestroyREAL8Vector(phaseNQC);
              }
              if(hamVHi){
                XLALDestroyREAL8Vector(hamVHi);
              }
              if(hLMAllHi){
                XLALDestroyREAL8Vector(hLMAllHi);
              }
              if(sigReVec){
                XLALDestroyREAL8Vector(sigReVec);
              }
              if(sigImVec){
                XLALDestroyREAL8Vector(sigImVec);
              }
              if(hLMAll){
                XLALDestroyREAL8Vector(hLMAll);
              }
              XLAL_ERROR (XLAL_EFUNC);
            }
            hT = 0.;
            if ( (lambda2Tidal1 != 0. && omega02Tidal1 != 0.) || (lambda2Tidal2 != 0. && omega02Tidal2 != 0.) ) {
            if (XLALSimIMRSpinEOBWaveformTidal
                (&hT, values, cbrt (omegaVec->data[i]), 2, 2, &seobParams)
                == XLAL_FAILURE)
                {
                    XLAL_ERROR (XLAL_EFUNC);
                }
            }

            if (XLALSimIMREOBNonQCCorrection (&hNQC, values, omegaVec->data[i], &nqcCoeffs)
                == XLAL_FAILURE)
            {
              if(tmpValues){
                XLALDestroyREAL8Vector (tmpValues);
              }
              if(sigmaKerr){
                XLALDestroyREAL8Vector (sigmaKerr);
              }
              if(sigmaStar){
                XLALDestroyREAL8Vector (sigmaStar);
              }
              if(values){
                XLALDestroyREAL8Vector (values);
              }
              if(sigReHi){
                XLALDestroyREAL8Vector(sigReHi);
              }
              if(sigImHi){
                XLALDestroyREAL8Vector(sigImHi);
              }
              if(omegaHi){
                XLALDestroyREAL8Vector(omegaHi);
              }
              if(ampNQC){
                XLALDestroyREAL8Vector(ampNQC);
              }
              if(phaseNQC){
                XLALDestroyREAL8Vector(phaseNQC);
              }
              if(hamVHi){
                XLALDestroyREAL8Vector(hamVHi);
              }
              if(hLMAllHi){
                XLALDestroyREAL8Vector(hLMAllHi);
              }
              if(sigReVec){
                XLALDestroyREAL8Vector(sigReVec);
              }
              if(sigImVec){
                XLALDestroyREAL8Vector(sigImVec);
              }
              if(hLMAll){
                XLALDestroyREAL8Vector(hLMAll);
              }
              XLAL_ERROR (XLAL_EFUNC);
            }

            hLM *= hNQC;
            hLM += hT;

            if (use_tidal==1) {
                REAL8 dtGeom = deltaTHigh / mTScaled;
                INT4 iEnd= (INT4)rdMatchPoint->data[1]/dtGeom;
                REAL8 omega0 = OmVec->data[iEnd];
                REAL8 tau = 0.5*LAL_PI/omega0;
                REAL8 dtGeomLow = deltaT / mTScaled;
                sigReVec->data[i] = amp0 * creal (hLM)/(1.  + exp(( i*dtGeomLow - (rdMatchPoint->data[1]+15 + (dynamics->data)[hiSRndx]) )/tau));
                sigImVec->data[i] = amp0 * cimag (hLM)/(1. + exp(( i*dtGeomLow - (rdMatchPoint->data[1] +15 + (dynamics->data)[hiSRndx]))/tau));
            }
            else {
                sigReVec->data[i] = amp0 * creal (hLM);
                sigImVec->data[i] = amp0 * cimag (hLM);
            }
            hLMAll->data[2*k*sigReVec->length + i] = sigReVec->data[i];
            hLMAll->data[(1+2*k)*sigReVec->length + i] = sigImVec->data[i];
        }

#if outputDebug
         fclose (out);
        fclose(out2);
#endif
  }
}
    if ( OmVec )
        XLALDestroyREAL8Vector(OmVec);
    if ( omegaVec )
         XLALDestroyREAL8Vector(omegaVec);

  /*
   * STEP 8) Generate full IMR modes -- attaching ringdown to inspiral
   */

  /* Attach the ringdown part to the inspiral */
for ( UINT4 k = 0; k<nModes; k++) {
  for (i = 0; i < (INT4) (sigReHi->length / resampFac); i++)
    {
      hLMAll->data[2*k*sigReVec->length + i + hiSRndx] = hLMAllHi->data[2*k*sigReHi->length + i*resampFac];
      hLMAll->data[(2*k+1)*sigReVec->length + i + hiSRndx] = hLMAllHi->data[(2*k+1)*sigReHi->length + i*resampFac];
    }
}



    /* Cut wf if fMin requested by user was high */
    INT4 kMin = 0;
    if ( fStart != fMin ) {
        REAL8 finst;
        gsl_spline *splineRe = NULL;
        gsl_interp_accel *accRe = NULL;
        gsl_spline *splineIm = NULL;
        gsl_interp_accel *accIm = NULL;
        //RC: since the cut of the waveform is done on the 22 mode, we assign here the 22 mode to sigReVec and sigImVec
        for ( i=0; i<(INT4)sigReVec->length; i++) {
          sigReVec->data[i] = hLMAll->data[i];
          sigImVec->data[i] = hLMAll->data[sigReVec->length + i];
        }
        REAL8Vector *tmpRe = XLALCreateREAL8Vector(sigReVec->length), *tmpIm = XLALCreateREAL8Vector(sigReVec->length);
        for ( i=0; i < (INT4) sigReVec->length; i++) {
            tmpRe->data[i] = sigReVec->data[i] / amp0;
            tmpIm->data[i] = sigImVec->data[i] / amp0;
        }
        splineRe = gsl_spline_alloc (gsl_interp_cspline, sigReVec->length);
        splineIm = gsl_spline_alloc (gsl_interp_cspline, sigImVec->length);
        accRe = gsl_interp_accel_alloc ();
        accIm = gsl_interp_accel_alloc ();
        REAL8 dRe, dIm;
        REAL8Vector *timeList;
        timeList = XLALCreateREAL8Vector (sigReVec->length);
        for ( i=0; i < (INT4) sigReVec->length; i++) {
            timeList->data[i] = i*deltaT/mTScaled;
        }
        gsl_spline_init (splineRe, timeList->data, tmpRe->data, tmpRe->length);
        gsl_spline_init (splineIm, timeList->data, tmpIm->data, tmpIm->length);
        REAL8 norm;
        for ( i=1; i < (INT4) tmpRe->length - 1; i++) {
            norm = tmpRe->data[i]*tmpRe->data[i] + tmpIm->data[i]*tmpIm->data[i];
            if ( norm > 0. ) {
                dRe = gsl_spline_eval_deriv (splineRe, timeList->data[i], accRe);
                dIm = gsl_spline_eval_deriv (splineIm, timeList->data[i], accIm);
                finst = (dRe*tmpIm->data[i] - dIm*tmpRe->data[i])/norm;
//                printf("%.16e %.16e\n", timeList->data[i], finst);
                finst = finst/(2.*LAL_PI*mTScaled);
//                printf("%.16e %.16e %.16e\n", timeList->data[i], finst, fMin);
                if ( finst > fMin ) {
                    kMin = i;
                    break;
                }
            }
            else {
                continue;
            }
        }
        gsl_spline_free( splineRe );
        gsl_interp_accel_free( accRe );
        gsl_spline_free( splineIm );
        gsl_interp_accel_free( accIm );
        XLALDestroyREAL8Vector( tmpRe );
        XLALDestroyREAL8Vector( tmpIm );
        XLALDestroyREAL8Vector(timeList);
    }

#if debugOutput
    for ( UINT4 k = 0; k<nModes; k++) {
        modeL  = lmModes[k][0];
        modeM = lmModes[k][1];
        char filename[sizeof "saModesXXIMR.dat"];
        sprintf(filename,"saModes%01d%01dIMR.dat",modeL,modeM);
        out = fopen (filename, "w");
        for ( i=0; i<(INT4)sigReVec->length; i++) {
            fprintf (out, "%.16e %.16e %.16e\n", i*deltaT / mTScaled,hLMAll->data[2*k*sigReVec->length + i]/amp0,hLMAll->data[(1+2*k)*sigReVec->length + i]/amp0);
        }
        fclose(out);
    }
#endif

XLALGPSAdd (&tc, deltaT * (REAL8) kMin);

/*
 * STEP 9) Generate full IMR hp and hx waveforms
 */
//RC: this function stops now here and return the array with the modes

/*
 * STEP 9) Return real and imaginary part of the modes
 */
  SphHarmTimeSeries *hlms = NULL;
  char mode_name[32];

for ( UINT4 k = 0; k<nModes; k++) {
  modeL  = lmModes[k][0];
  modeM = lmModes[k][1];
  snprintf(mode_name, sizeof(mode_name), "(%d, %d) mode", modeL, modeM);
  COMPLEX16TimeSeries *tmp_mode = XLALCreateCOMPLEX16TimeSeries(mode_name, &tc, 0.0,
    deltaT, &lalStrainUnit, sigReVec->length - kMin);
  for (UINT4 t = kMin; t< (UINT4) sigReVec->length; t++)
      {
            tmp_mode->data->data[t - kMin]  = hLMAll->data[2*k*sigReVec->length + t];
            tmp_mode->data->data[t - kMin]  += I * hLMAll->data[(2*k+1)*sigReVec->length + t];
        }
        hlms = XLALSphHarmTimeSeriesAddMode(hlms, tmp_mode, modeL, modeM);
        XLALDestroyCOMPLEX16TimeSeries(tmp_mode);
      }
      /* Point the output pointers to the relevant time series and return */
      (*hlmmode) = hlms;

      /* Free memory */
      XLALDestroyREAL8Vector (tmpValues);
      XLALDestroyREAL8Vector (sigmaKerr);
      XLALDestroyREAL8Vector (sigmaStar);
      XLALDestroyREAL8Vector (values);
      XLALDestroyREAL8Vector (rdMatchPoint);
      XLALDestroyREAL8Vector (ampNQC);
      XLALDestroyREAL8Vector (phaseNQC);
      XLALDestroyREAL8Vector (sigReVec);
      XLALDestroyREAL8Vector (sigImVec);
      XLALAdaptiveRungeKuttaFree (integrator);
      //SM
      //XLALDestroyREAL8Array (dynamics);
      //XLALDestroyREAL8Array (dynamicsHi);
      //SM
      gsl_matrix_free (nqcCoeffsMatrix);

      if (dynamicstmp)
        {
          XLALDestroyREAL8Array (dynamicstmp);	// DAVIDS: We are done with these now
        }
      if (dynamicsHitmp)
        {
          XLALDestroyREAL8Array (dynamicsHitmp);	// DAVIDS: Done with these now
        }

      XLALDestroyREAL8Vector (sigReHi);
      XLALDestroyREAL8Vector (sigImHi);
      XLALDestroyREAL8Vector (omegaHi);
      if ( hLMAllHi )
          XLALDestroyREAL8Vector (hLMAllHi);
      if ( hLMAll )
          XLALDestroyREAL8Vector (hLMAll);
      if ( hamV )
          XLALDestroyREAL8Vector (hamV);
      if ( hamVHi )
          XLALDestroyREAL8Vector (hamVHi);

      //SM
      // Copy dynamics to output in the form of a REAL8Vector (required for SWIG wrapping, REAL8Array does not work)
      *dynamics_out = XLALCreateREAL8Vector(5 * retLen_out);
      *dynamicsHi_out = XLALCreateREAL8Vector(5 * retLenHi_out);
      for (i = 0; i < 5*retLen_out; i++) (*dynamics_out)->data[i] = dynamics->data[i];
      // We have to add the starting time of the high-sampling dynamics, as the output of the integrator starts with 0
      for (i = 0; i < retLenHi_out; i++) (*dynamicsHi_out)->data[i] = tstartHi + dynamicsHi->data[i];
      for (i = retLenHi_out; i < 5*retLenHi_out; i++) (*dynamicsHi_out)->data[i] = dynamicsHi->data[i];
      XLALDestroyREAL8Array (dynamics);
      XLALDestroyREAL8Array (dynamicsHi);
      //SM

      return XLAL_SUCCESS;
    }

/**
 * This function takes the modes from the function XLALSimIMRSpinAlignedEOBModes and combine them into h+ and hx
 */

int
XLALSimIMRSpinAlignedEOBWaveformAll (REAL8TimeSeries ** hplus,
				     /**<< OUTPUT, real part of the modes */
				     REAL8TimeSeries ** hcross,
				     /**<< OUTPUT, complex part of the modes */
				     const REAL8 phiC,
				     /**<< coalescence orbital phase (rad) */
				     REAL8 deltaT,
				     /**<< sampling time step */
				     const REAL8 m1SI,
				     /**<< mass-1 in SI unit */
				     const REAL8 m2SI,
				     /**<< mass-2 in SI unit */
				     const REAL8 fMin,
				     /**<< starting frequency of the 22 mode (Hz) */
				     const REAL8 r,
				     /**<< distance in SI unit */
				     const REAL8 inc,
				     /**<< inclination angle */
				     const REAL8 spin1z,
				     /**<< z-component of spin-1, dimensionless */
				     const REAL8 spin2z,
				      /**<< z-component of spin-2, dimensionless */
                     UINT4 SpinAlignedEOBversion,
                     /**<< 1 for SEOBNRv1, 2 for SEOBNRv2, 4 for SEOBNRv4, 201 for SEOBNRv2T, 401 for SEOBNRv4T, 41 for SEOBNRv4HM */
				     const REAL8 lambda2Tidal1,
                     /**<< dimensionless adiabatic quadrupole tidal deformability for body 1 (2/3 k2/C^5) */
				     const REAL8 lambda2Tidal2,
                     /**<< dimensionless adiabatic quadrupole tidal deformability for body 2 (2/3 k2/C^5) */
				     const REAL8 omega02Tidal1,
                     /**<< quadrupole f-mode angular freq for body 1 m_1*omega_{02,1}*/
				     const REAL8 omega02Tidal2,
                      /**<< quadrupole f-mode angular freq for body 2 m_2*omega_{02,2}*/
				     const REAL8 lambda3Tidal1,
                     /**<< dimensionless adiabatic octupole tidal deformability for body 1 (2/15 k3/C^7) */
				     const REAL8 lambda3Tidal2,
                     /**<< dimensionless adiabatic octupole tidal deformability for body 2 (2/15 k3/C^7) */
				     const REAL8 omega03Tidal1,
                     /**<< octupole f-mode angular freq for body 1 m_1*omega_{03,1}*/
				     const REAL8 omega03Tidal2,
                     /**<< octupole f-mode angular freq for body 2 m_2*omega_{03,2}*/
             const REAL8 quadparam1,
                     /**<< parameter kappa_1 of the spin-induced quadrupole for body 1, quadrupole is Q_A = -kappa_A m_A^3 chi_A^2 */
				     const REAL8 quadparam2,
                     /**<< parameter kappa_2 of the spin-induced quadrupole for body 2, quadrupole is Q_A = -kappa_A m_A^3 chi_A^2 */
                     REAL8Vector *nqcCoeffsInput,
                     /**<< Input NQC coeffs */
                     const INT4 nqcFlag,
                     /**<< Flag to tell the code to use the NQC coeffs input thorugh nqcCoeffsInput */
            LALValue *ModeArray
            /**<< Structure containing the modes to use in the waveform */
  )
  {

    REAL8 coa_phase = phiC;

    SphHarmTimeSeries *hlms = NULL;
    //SM
    REAL8Vector *dynamics = NULL;
    REAL8Vector *dynamicsHi = NULL;
    //SM

    //RC: XLALSimIMRSpinAlignedEOBModes computes the modes and put them into hlm

    if(XLALSimIMRSpinAlignedEOBModes (&hlms,
                                   //SM
                                   &dynamics, &dynamicsHi,
                                   //SM
                                   deltaT, m1SI, m2SI, fMin, r, spin1z, spin2z, SpinAlignedEOBversion,
                                               lambda2Tidal1, lambda2Tidal2,
                                               omega02Tidal1, omega02Tidal2,
                                               lambda3Tidal1, lambda3Tidal2,
                                               omega03Tidal1, omega03Tidal2,
                                               quadparam1, quadparam2,
                                               nqcCoeffsInput, nqcFlag) == XLAL_FAILURE){
                                                 if(dynamics) XLALDestroyREAL8Vector(dynamics);
                                                 if(dynamicsHi) XLALDestroyREAL8Vector(dynamicsHi);
                                                 XLAL_ERROR (XLAL_EFUNC);
                                               };

    //RC: For SEOBNRv4T we also need to exit from this function when  nqcFlag == 1 because when this flag is 1, it is only computing the NQCs and not the wf
    if (nqcFlag == 1){
      if(hlms)
        XLALDestroySphHarmTimeSeries(hlms);
      return XLAL_SUCCESS;
    }

    //RC: Here we read lenght and epoch of the modes. They are all the same by definition.
    INT4 len = hlms->mode->data->length;
    LIGOTimeGPS tc = hlms->mode->epoch;

    //RC: defining and initializing to 0 the hp and hc vectors
    REAL8TimeSeries *hPlusTS =
      XLALCreateREAL8TimeSeries ("H_PLUS", &tc, 0.0, deltaT, &lalStrainUnit,
               len);
    REAL8TimeSeries *hCrossTS =
      XLALCreateREAL8TimeSeries ("H_CROSS", &tc, 0.0, deltaT, &lalStrainUnit,
              len);
    memset( hPlusTS->data->data, 0, hPlusTS->data->length * sizeof(REAL8) );
    memset( hCrossTS->data->data, 0, hCrossTS->data->length * sizeof(REAL8) );

    //RC: adding all the modes in the SphHarmTimeSeries hlms to hp and hc
    SphHarmTimeSeries *hlms_temp = hlms;
    while ( hlms_temp ) {
      if (XLALSimInspiralModeArrayIsModeActive(ModeArray, hlms_temp->l, hlms_temp->m) == 1){
        /*Here we check if the mode generated is in the ModeArray structure
        */
          //R.C. the angles in the spin-weighted spherical harmonics are set accordingly to the document https://dcc.ligo.org/DocDB/0152/T1800226/003/LAL_GW_Frames.pdf
          //for the conventions in master
          XLALSimAddMode(hPlusTS, hCrossTS, hlms_temp->mode, inc, LAL_PI/2. - coa_phase, hlms_temp->l, hlms_temp->m, 1);
          /*When the function XLALSimAddMode is called with the last argument == 1
          *is using both +m and -m modes
          */
          //printf("Mode (%d,%d) active\n", hlms_temp->l,hlms_temp->m);
        }
        hlms_temp = hlms_temp->next;
    }
    


    /* Point the output pointers to the relevant time series and return */
    (*hplus) = hPlusTS;
    (*hcross) = hCrossTS;

    if(hlms)
      XLALDestroySphHarmTimeSeries(hlms);

    //SM
    if(dynamics) XLALDestroyREAL8Vector(dynamics);
    if(dynamicsHi) XLALDestroyREAL8Vector(dynamicsHi);
    //SM

    return XLAL_SUCCESS;
  }

/** @} */
