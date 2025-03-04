
/*
 * Copyright (C) 2018/2019 Geraint Pratten
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

 #ifdef __cplusplus
 extern "C" {
 #endif

/**
 * \author Geraint Pratten
 *
 */

#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_sf_ellint.h>

#include <lal/SphericalHarmonics.h>
#include "LALSimIMR.h"

#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomX_precession.h"
#include "LALSimIMRPhenomX_PNR.h"
#include "LALSimIMRPhenomXHM_qnm.h"

#include <lal/LALAdaptiveRungeKuttaIntegrator.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


#ifndef _OPENMP
#define omp ignore
#endif

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif

/* ~~~~~~~~~~ Functions to Perform Frame Transformations and Populate Structs ~~~~~~~~~~ */
 /**
     Function to populate the IMRPhenomXPrecessionStruct:
       - Calculates frame transformation
       - PN Euler angles
       - Frame transformations
       - Orbital angular momenta etc
 */
int IMRPhenomXGetAndSetPrecessionVariables(
  IMRPhenomXWaveformStruct *pWF,
  IMRPhenomXPrecessionStruct *pPrec,
  REAL8 m1_SI,
  REAL8 m2_SI,
  REAL8 chi1x,
  REAL8 chi1y,
  REAL8 chi1z,
  REAL8 chi2x,
  REAL8 chi2y,
  REAL8 chi2z,
  LALDict *lalParams,
  INT4 debug_flag
)
{
  /*
      Here we assume m1 > m2, q > 1, dm = m1 - m2 = delta = sqrt(1-4eta) > 0
  */
  pWF->LALparams = lalParams;

  /* Pre-cache useful powers here */
  pPrec->sqrt2   = 1.4142135623730951;
  pPrec->sqrt5   = 2.23606797749978981;
  pPrec->sqrt6   = 2.44948974278317788;
  pPrec->sqrt7   = 2.64575131106459072;
  pPrec->sqrt10  = 3.16227766016838;
  pPrec->sqrt14  = 3.74165738677394133;
  pPrec->sqrt15  = 3.87298334620741702;
  pPrec->sqrt70  = 8.36660026534075563;
  pPrec->sqrt30  = 5.477225575051661;
  pPrec->sqrt2p5 = 1.58113883008419;

  pPrec->debug_prec = debug_flag;

  // Get IMRPhenomX precession version from LAL dictionary
  pPrec->IMRPhenomXPrecVersion = XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams);
  if (pPrec->IMRPhenomXPrecVersion == 300) pPrec->IMRPhenomXPrecVersion = 223;

  REAL8 chi_in_plane = sqrt(chi1x*chi1x+chi1y*chi1y+chi2x*chi2x+chi2y*chi2y);
  if(chi_in_plane<1e-7 && (pPrec->IMRPhenomXPrecVersion==320||pPrec->IMRPhenomXPrecVersion==321||pPrec->IMRPhenomXPrecVersion==310||pPrec->IMRPhenomXPrecVersion==311))
  {
  pPrec->IMRPhenomXPrecVersion=102;
  }

  if( (pPrec->IMRPhenomXPrecVersion==320||pPrec->IMRPhenomXPrecVersion==321||pPrec->IMRPhenomXPrecVersion==310||pPrec->IMRPhenomXPrecVersion==311))
  {

  int status=XLAL_SUCCESS;
  pPrec->PNarrays = XLALMalloc(sizeof(PhenomXPInspiralArrays));

  // check mode array to estimate frequency range over which splines will need to be evaluated
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);
  if (ModeArray != NULL)
  {
   IMRPhenomX_GetandSetModes(ModeArray,pPrec);
   XLALDestroyValue(ModeArray);
  }

  // buffer for GSL interpolation to succeed
  REAL8 buffer=(pWF->deltaF>0.)? 3.*pWF->deltaF: 0.5;
  REAL8 flow=(pWF->fMin-buffer)*2./pPrec->M_MAX;
  XLAL_CHECK(flow>0.,XLAL_EDOM,"Error in %s: starting frequency for SpinTaylor angles must be positive!",__func__);
  status=IMRPhenomX_InspiralAngles_SpinTaylor(pPrec->PNarrays,chi1x,chi1y,chi1z,chi2x,chi2y,chi2z,flow,pPrec->IMRPhenomXPrecVersion,pWF,lalParams);

  // if PN numerical integration fails, default to MSA+fallback to NNLO
  if(status==XLAL_FAILURE) {
                            LALFree(pPrec->PNarrays);
                            XLAL_PRINT_WARNING("Warning: due to a failure in the SpinTaylor routines, the model will default to MSA angles.");
                            pPrec->IMRPhenomXPrecVersion=223;
                            }
 // end of SpinTaylor code

            }

  // Get expansion order for MSA system of equations. Default is taken to be 5.
  pPrec->ExpansionOrder        = XLALSimInspiralWaveformParamsLookupPhenomXPExpansionOrder(lalParams);

  // Get toggle for PNR angles
  INT4 PNRUseTunedAngles = XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(lalParams);
  pPrec->IMRPhenomXPNRUseTunedAngles = PNRUseTunedAngles;

  // Get PNR angle interpolation tolerance
  pPrec->IMRPhenomXPNRInterpTolerance = XLALSimInspiralWaveformParamsLookupPhenomXPNRInterpTolerance(lalParams);

  // Get toggle for symmetric waveform
  INT4 AntisymmetricWaveform = XLALSimInspiralWaveformParamsLookupPhenomXAntisymmetricWaveform(lalParams);
  pPrec->IMRPhenomXAntisymmetricWaveform = AntisymmetricWaveform;

  // Set toggle for polarization calculation: +1 for symmetric waveform (default), -1 for antisymmetric waveform; refer to XXXX.YYYYY for details
  pPrec->PolarizationSymmetry = 1.0;

  //
  int pflag = pPrec->IMRPhenomXPrecVersion;
  if(pflag != 101 && pflag != 102 && pflag != 103 && pflag != 104 && pflag != 220 && pflag != 221 && pflag != 222 && pflag != 223 && pflag != 224 && pflag!=310 && pflag!=311 && pflag!=320 && pflag!=321)
  {
    XLAL_ERROR(XLAL_EINVAL, "Error in IMRPhenomXGetAndSetPrecessionVariables: Invalid precession flag. Allowed versions are 101, 102, 103, 104, 220, 221, 222, 223, 224, 310, 311, 320 or 321.\n");
  }

  switch( pflag )
    {
        case 101: // NNLO single spin PNEuler angles + 2PN non-spinning L
        case 102: // NNLO single spin PNEuler angles + 3PN spinning L
        case 103: // NNLO single spin PNEuler angles + 4PN spinning L
        case 104: // NNLO single spin PNEuler angles + 4PN spinning L + LOS terms in L
    {
      break;
    }
    case 220: // MSA using expressions as detailed in arXiv:1703.03967. Defaults to NNLO v102 if MSA fails.
    case 221: // MSA using expressions as detailed in arXiv:1703.03967. Terminal failure if MSA fails.
    case 222: // MSA using expressions as implemented in LALSimInspiralFDPrecAngles. Terminal failure if MSA fails.
    case 223: // MSA using expressions as implemented in LALSimInspiralFDPrecAngles. Defaults to NNLO v102 if MSA fails.
    case 224: // MSA using expressions as detailed in arXiv:1703.03967, with \zeta_0 and \phi_{z,0} as in LALSimInspiralFDPrecAngles.  Defaults to NNLO v102 if MSA fails.
    {
       /*
          Double-spin model using angles from Chatziioannou et al, PRD, 95, 104004, (2017), arXiv:1703.03967
          Uses 3PN L
       */
       #if DEBUG == 1
        printf("Initializing MSA system...\n");
       #endif

       if(pPrec->ExpansionOrder < -1 || pPrec->ExpansionOrder > 5)
       {
         XLAL_ERROR(XLAL_EINVAL, "Error in IMRPhenomXGetAndSetPrecessionVariables: Invalid expansion order for MSA corrections. Default is 5, allowed values are [-1,0,1,2,3,4,5].\n");
       }
       break;

    }

    case 310: // Numerical integration of SpinTaylor equations, constant angles in MRD
    case 311: // Numerical integration of SpinTaylor equations, constant angles in MRD, BBH precession
    case 320: // Numerical integration of SpinTaylor equations, analytical continuation in MRD
    case 321: // Numerical integration of SpinTaylor equations, analytical continuation in MRD, BBH precession
        {
           break;
        }


        default:
        {
            XLAL_ERROR(XLAL_EINVAL, "Error in IMRPhenomXGetAndSetPrecessionVariables: IMRPhenomXPrecessionVersion not recognized.\n");
      break;
        }
    }

  // get first digit of precessing version: this tags the method employed to compute the Euler angles
  // 1: NNLO; 2: MSA; 3: SpinTaylor (numerical)
  int precversionTag=(pPrec->IMRPhenomXPrecVersion-(pPrec->IMRPhenomXPrecVersion%100))/100;
  pPrec->precessing_tag=precversionTag;

    /* Define a number of convenient local parameters */
  const REAL8 m1        = m1_SI / pWF->Mtot_SI;   /* Normalized mass of larger companion:   m1_SI / Mtot_SI */
  const REAL8 m2        = m2_SI / pWF->Mtot_SI;   /* Normalized mass of smaller companion:  m2_SI / Mtot_SI */
  const REAL8 M         = (m1 + m2);              /* Total mass in solar units */

  // Useful powers of mass
  const REAL8 m1_2      = m1 * m1;
  const REAL8 m1_3      = m1 * m1_2;
  const REAL8 m1_4      = m1 * m1_3;
  const REAL8 m1_5      = m1 * m1_4;
  const REAL8 m1_6      = m1 * m1_5;
  const REAL8 m1_7      = m1 * m1_6;
  const REAL8 m1_8      = m1 * m1_7;

  const REAL8 m2_2      = m2 * m2;

  pWF->M = M;
  pWF->m1_2 = m1_2;
  pWF->m2_2 = m2_2;

  const REAL8 q         = m1/m2; // q = m1 / m2 > 1.0

  // Powers of eta
  const REAL8 eta       = pWF->eta;
  const REAL8 eta2      = eta*eta;
  const REAL8 eta3      = eta*eta2;
  const REAL8 eta4      = eta*eta3;
  const REAL8 eta5      = eta*eta4;
  const REAL8 eta6      = eta*eta5;

  // \delta in terms of q > 1
  const REAL8 delta     = pWF->delta;
  const REAL8 delta2    = delta*delta;
  const REAL8 delta3    = delta*delta2;

  // Cache these powers, as we use them regularly
  pPrec->eta            = eta;
  pPrec->eta2           = eta2;
  pPrec->eta3           = eta3;
  pPrec->eta4           = eta4;

  pPrec->inveta         = 1.0 / eta;
  pPrec->inveta2        = 1.0 / eta2;
  pPrec->inveta3        = 1.0 / eta3;
  pPrec->inveta4        = 1.0 / eta4;
  pPrec->sqrt_inveta    = 1.0 / sqrt(eta);

  const REAL8 chi_eff   = pWF->chiEff;

  pPrec->twopiGM        = LAL_TWOPI * LAL_G_SI * (m1_SI + m2_SI) / LAL_C_SI / LAL_C_SI / LAL_C_SI;
  pPrec->piGM           = LAL_PI * (m1_SI + m2_SI) * (LAL_G_SI / LAL_C_SI) / (LAL_C_SI * LAL_C_SI);

  /* Set spin variables in pPrec struct */
  pPrec->chi1x          = chi1x;
  pPrec->chi1y          = chi1y;
  pPrec->chi1z          = chi1z;
  pPrec->chi1_norm      = sqrt(chi1x*chi1x + chi1y*chi1y + chi1z*chi1z);

  pPrec->chi2x          = chi2x;
  pPrec->chi2y          = chi2y;
  pPrec->chi2z          = chi2z;
  pPrec->chi2_norm      = sqrt(chi2x*chi2x + chi2y*chi2y + chi2z*chi2z);

  /* Check that spins obey Kerr bound */
  if((!PNRUseTunedAngles)||(pWF->PNR_SINGLE_SPIN != 1)){ /*Allow the single-spin mapping for PNR to break the Kerr limit*/
    XLAL_CHECK(fabs(pPrec->chi1_norm) <= 1.0, XLAL_EDOM, "Error in IMRPhenomXSetPrecessionVariables: |S1/m1^2| must be <= 1.\n");
    XLAL_CHECK(fabs(pPrec->chi2_norm) <= 1.0, XLAL_EDOM, "Error in IMRPhenomXSetPrecessionVariables: |S2/m2^2| must be <= 1.\n");
  }

  /* Calculate dimensionful spins */
  pPrec->S1x        = chi1x * m1_2;
  pPrec->S1y        = chi1y * m1_2;
  pPrec->S1z        = chi1z * m1_2;
  pPrec->S1_norm    = fabs(pPrec->chi1_norm) * m1_2;

  pPrec->S2x        = chi2x * m2_2;
  pPrec->S2y        = chi2y * m2_2;
  pPrec->S2z        = chi2z * m2_2;
  pPrec->S2_norm    = fabs(pPrec->chi2_norm) * m2_2;

  // Useful powers
  pPrec->S1_norm_2  = pPrec->S1_norm * pPrec->S1_norm;
  pPrec->S2_norm_2  = pPrec->S2_norm * pPrec->S2_norm;

  pPrec->chi1_perp  = sqrt(chi1x*chi1x + chi1y*chi1y);
  pPrec->chi2_perp  = sqrt(chi2x*chi2x + chi2y*chi2y);

  /* Get spin projections */
  pPrec->S1_perp    = (m1_2) * sqrt(chi1x*chi1x + chi1y*chi1y);
  pPrec->S2_perp    = (m2_2) * sqrt(chi2x*chi2x + chi2y*chi2y);

  /* Norm of in-plane vector sum: Norm[ S1perp + S2perp ] */
  pPrec->STot_perp     = sqrt( (pPrec->S1x+pPrec->S2x)*(pPrec->S1x+pPrec->S2x) + (pPrec->S1y+pPrec->S2y)*(pPrec->S1y+pPrec->S2y) );

  /* This is called chiTot_perp to distinguish from Sperp used in contrusction of chi_p. For normalization, see Sec. IV D of arXiv:2004.06503 */
  pPrec->chiTot_perp   = pPrec->STot_perp * (M*M) / m1_2;
  /* Store pPrec->chiTot_perp to pWF so that it can be used in XCP modifications (PNRUseTunedCoprec) */
  pWF->chiTot_perp = pPrec->chiTot_perp;

  /* disable tuned PNR angles, tuned coprec and mode asymmetries in low in-plane spin limit */
  if((chi_in_plane < 1e-7)&&(pPrec->IMRPhenomXPNRUseTunedAngles == 1)&&(pWF->PNR_SINGLE_SPIN != 1)){
    XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedAngles(lalParams, 0);
    PNRUseTunedAngles = 0;
    pPrec->IMRPhenomXPNRUseTunedAngles = 0;
    pPrec->IMRPhenomXAntisymmetricWaveform = 0;
    AntisymmetricWaveform = 0;
    XLALSimInspiralWaveformParamsInsertPhenomXAntisymmetricWaveform(lalParams, 0);
    XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedCoprec(lalParams, 0);
  }

  /*
    Calculate the effective precessing spin parameter (Schmidt et al, PRD 91, 024043, 2015):
      - m1 > m2, so body 1 is the larger black hole
  */
  pPrec->A1             = 2.0 + (3.0 * m2) / (2.0 * m1);
  pPrec->A2             = 2.0 + (3.0 * m1) / (2.0 * m2);
  pPrec->ASp1           = pPrec->A1 * pPrec->S1_perp;
  pPrec->ASp2           = pPrec->A2 * pPrec->S2_perp;

  /* S_p = max(A1 S1_perp, A2 S2_perp) */
  const REAL8 num       = (pPrec->ASp2 > pPrec->ASp1) ? pPrec->ASp2 : pPrec->ASp1;
  const REAL8 den       = (m2 > m1) ? pPrec->A2*(m2_2) : pPrec->A1*(m1_2);

  /* chi_p = max(A1 * Sp1 , A2 * Sp2) / (A_i * m_i^2) where i is the index of the larger BH */
  const REAL8 chip      = num / den;
  const REAL8 chi1L     = chi1z;
  const REAL8 chi2L     = chi2z;


  pPrec->chi_p          = chip;
  // (PNRUseTunedCoprec)
  pWF->chi_p             = pPrec->chi_p;
  pPrec->phi0_aligned   = pWF->phi0;

  /* Effective (dimensionful) aligned spin */
  pPrec->SL             = chi1L*m1_2 + chi2L*m2_2;

  /* Effective (dimensionful) in-plane spin */
  pPrec->Sperp          = chip * m1_2;                  /* m1 > m2 */

  pPrec->MSA_ERROR      = 0;

  pPrec->pWF22AS = NULL;

  /* Calculate parameter for two-spin to single-spin map used in PNR and XCP */
  /* Initialize PNR variables */
  pPrec->chi_singleSpin = 0.0;
  pPrec->costheta_singleSpin = 0.0;
  pPrec->costheta_final_singleSpin = 0.0;
  pPrec->chi_singleSpin_antisymmetric = 0.0;
  pPrec->theta_antisymmetric = 0.0;
  pPrec->PNR_HM_Mflow = 0.0;
  pPrec->PNR_HM_Mfhigh = 0.0;

  pPrec->PNR_q_window_lower = 0.0;
  pPrec->PNR_q_window_upper = 0.0;
  pPrec->PNR_chi_window_lower = 0.0;
  pPrec->PNR_chi_window_upper = 0.0;
  // pPrec->PNRInspiralScaling = 0;

  UINT4 status = IMRPhenomX_PNR_GetAndSetPNRVariables(pWF, pPrec);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_PNR_GetAndSetPNRVariables failed in IMRPhenomXGetAndSetPrecessionVariables.\n");

  pPrec->alphaPNR = 0.0;
  pPrec->betaPNR = 0.0;
  pPrec->gammaPNR = 0.0;

  /*...#...#...#...#...#...#...#...#...#...#...#...#...#...#.../
  /      Get and/or store CoPrec params into pWF and pPrec     /
  /...#...#...#...#...#...#...#...#...#...#...#...#...#...#...*/

  status = IMRPhenomX_PNR_GetAndSetCoPrecParams(pWF,pPrec,lalParams);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC,
  "Error: IMRPhenomX_PNR_GetAndSetCoPrecParams failed \
  in IMRPhenomXGetAndSetPrecessionVariables.\n");

  /*..#...#...#...#...#...#...#...#...#...#...#...#...#...#...*/



  //
  if( pflag == 220 || pflag == 221 || pflag == 222 || pflag == 223 || pflag == 224 )
    {
      #if DEBUG == 1
        printf("Evaluating MSA system.\n");
        printf("Expansion Order : %d\n",pPrec->ExpansionOrder);
      #endif

      IMRPhenomX_Initialize_MSA_System(pWF,pPrec,pPrec->ExpansionOrder);

      if(pPrec->MSA_ERROR == 1)
      {
        // In version 220, 223 and 224 if the MSA system fails to initialize we default to the NNLO PN angles using the 3PN aligned-spin orbital angular momentum
        if(pflag == 220 || pflag == 223 || pflag == 224)
        {
          XLAL_PRINT_WARNING("Warning: Initialization of MSA system failed. Defaulting to NNLO angles using 3PN aligned-spin approximation.");
          pPrec->IMRPhenomXPrecVersion = 102;
          pflag  = pPrec->IMRPhenomXPrecVersion;
        }
        else // Otherwise, if the MSA system fails to initialize we trigger a terminal error
        {
          XLAL_ERROR(XLAL_EDOM,"Error: IMRPhenomX_Initialize_MSA_System failed to initialize. Terminating.\n");
        }
      }
    }

  #if DEBUG == 1
    printf("In IMRPhenomXSetPrecessionVariables... \n\n");
    printf("chi_p   : %e\n",pPrec->chi_p);
    printf("phic    : %e\n",pPrec->phi0_aligned);
    printf("SL      : %e\n",pPrec->SL);
    printf("Sperp   : %e\n\n",pPrec->Sperp);
  #endif

  /*...#...#...#...#...#...#...#...#...#...#...#...#...#...#.../
  /      Compute and set final spin and RD frequency           /
  /...#...#...#...#...#...#...#...#...#...#...#...#...#...#...*/
  IMRPhenomX_SetPrecessingRemnantParams(pWF,pPrec,lalParams);
  /*..#...#...#...#...#...#...#...#...#...#...#...#...#...#...*/

  /* Useful powers of \chi_p */
  const REAL8 chip2    = chip * chip;

  /* Useful powers of spins aligned with L */
  const REAL8 chi1L2   = chi1L * chi1L;
  const REAL8 chi2L2   = chi2L * chi2L;

  const REAL8 log16    = 2.772588722239781;

  /*  Cache the orbital angular momentum coefficients for future use.

      References:
        - Kidder, PRD, 52, 821-847, (1995), arXiv:gr-qc/9506022
        - Blanchet, LRR, 17, 2, (2014), arXiv:1310.1528
        - Bohe et al, 1212.5520v2
        - Marsat, CQG, 32, 085008, (2015), arXiv:1411.4118
  */
  switch( pflag )
  {
    /* 2PN non-spinning orbital angular momentum (as per IMRPhenomPv2) */
    case 101:
    {
      pPrec->L0   = 1.0;
      pPrec->L1   = 0.0;
      pPrec->L2   = ((3.0/2.0) + (eta/6.0));
      pPrec->L3   = 0.0;
      pPrec->L4   = (81.0 + (-57.0 + eta)*eta)/24.;
      pPrec->L5   = 0.0;
      pPrec->L6   = 0.0;
      pPrec->L7   = 0.0;
      pPrec->L8   = 0.0;
      pPrec->L8L  = 0.0;
      break;
    }
    /* 3PN orbital angular momentum */
    case 102:
    case 220:
    case 221:
    case 224:
    case 310:
    case 311:
    case 320:
    case 321:
    {
      pPrec->L0   = 1.0;
      pPrec->L1   = 0.0;
      pPrec->L2   = 3.0/2. + eta/6.0;
      pPrec->L3   = (5*(chi1L*(-2 - 2*delta + eta) + chi2L*(-2 + 2*delta + eta)))/6.;
      pPrec->L4   = (81 + (-57 + eta)*eta)/24.;
      pPrec->L5   = (-7*(chi1L*(72 + delta*(72 - 31*eta) + eta*(-121 + 2*eta)) + chi2L*(72 + eta*(-121 + 2*eta) + delta*(-72 + 31*eta))))/144.;
      pPrec->L6   = (10935 + eta*(-62001 + eta*(1674 + 7*eta) + 2214*powers_of_lalpi.two))/1296.;
      pPrec->L7   = 0.0;
      pPrec->L8   = 0.0;

      // This is the log(x) term
      pPrec->L8L  = 0.0;
      break;

    }
    /* 3PN orbital angular momentum using non-conserved spin norms as per LALSimInspiralFDPrecAngles.c  */
    case 222:
    case 223:
    {
      pPrec->L0   = 1.0;
      pPrec->L1   = 0.0;
      pPrec->L2   = 3.0/2. + eta/6.0;
      pPrec->L3   = (-7*(chi1L + chi2L + chi1L*delta - chi2L*delta) + 5*(chi1L + chi2L)*eta)/6.;
      pPrec->L4   = (81 + (-57 + eta)*eta)/24.;
      pPrec->L5   = (-1650*(chi1L + chi2L + chi1L*delta - chi2L*delta) + 1336*(chi1L + chi2L)*eta + 511*(chi1L - chi2L)*delta*eta + 28*(chi1L + chi2L)*eta2)/600.;
      pPrec->L6   = (10935 + eta*(-62001 + 1674*eta + 7*eta2 + 2214*powers_of_lalpi.two))/1296.;
      pPrec->L7   = 0.0;
      pPrec->L8   = 0.0;

      // This is the log(x) term
      pPrec->L8L  = 0.0;
      break;
    }
    /* 4PN orbital angular momentum */
    case 103:
    {
      pPrec->L0   = 1.0;
      pPrec->L1   = 0.0;
      pPrec->L2   = 3.0/2. + eta/6.0;
      pPrec->L3   = (5*(chi1L*(-2 - 2*delta + eta) + chi2L*(-2 + 2*delta + eta)))/6.;
      pPrec->L4   = (81 + (-57 + eta)*eta)/24.;
      pPrec->L5   = (-7*(chi1L*(72 + delta*(72 - 31*eta) + eta*(-121 + 2*eta)) + chi2L*(72 + eta*(-121 + 2*eta) + delta*(-72 + 31*eta))))/144.;
      pPrec->L6   = (10935 + eta*(-62001 + eta*(1674 + 7*eta) + 2214*powers_of_lalpi.two))/1296.;
      pPrec->L7   = (chi2L*(-324 + eta*(1119 - 2*eta*(172 + eta)) + delta*(324 + eta*(-633 + 14*eta)))
                          - chi1L*(324 + eta*(-1119 + 2*eta*(172 + eta)) + delta*(324 + eta*(-633 + 14*eta))))/32.;
      pPrec->L8   = 2835/128. - (eta*(-10677852 + 100*eta*(-640863 + eta*(774 + 11*eta))
                      + 26542080*LAL_GAMMA + 675*(3873 + 3608*eta)*powers_of_lalpi.two))/622080. - (64*eta*log16)/3.;

      pPrec->L8L  = -(64.0/3.0) * eta;
      break;
    }
    /*
        4PN orbital angular momentum + leading order in spin at all PN orders terms.
          - Marsat, CQG, 32, 085008, (2015), arXiv:1411.4118
          - Siemonsen et al, PRD, 97, 064010, (2018), arXiv:1606.08832
    */
    case 104:
    {
      pPrec->L0   = 1.0;
      pPrec->L1   = 0.0;
      pPrec->L2   = 3.0/2. + eta/6.0;
      pPrec->L3   = (5*(chi1L*(-2 - 2*delta + eta) + chi2L*(-2 + 2*delta + eta)))/6.;
      pPrec->L4   = (81 + (-57 + eta)*eta)/24.;
      pPrec->L5   = (-7*(chi1L*(72 + delta*(72 - 31*eta) + eta*(-121 + 2*eta)) + chi2L*(72 + eta*(-121 + 2*eta) + delta*(-72 + 31*eta))))/144.;
      pPrec->L6   = (10935 + eta*(-62001 + eta*(1674 + 7*eta) + 2214*powers_of_lalpi.two))/1296.;
      pPrec->L7   = (chi2L*(-324 + eta*(1119 - 2*eta*(172 + eta)) + delta*(324 + eta*(-633 + 14*eta)))
                          - chi1L*(324 + eta*(-1119 + 2*eta*(172 + eta)) + delta*(324 + eta*(-633 + 14*eta))))/32.;
      pPrec->L8   = 2835/128. - (eta*(-10677852 + 100*eta*(-640863 + eta*(774 + 11*eta))
                      + 26542080*LAL_GAMMA + 675*(3873 + 3608*eta)*powers_of_lalpi.two))/622080. - (64*eta*log16)/3.;

      // This is the log(x) term at 4PN, x^4/2 * log(x)
      pPrec->L8L  = -(64.0/3.0) * eta;

      // Leading order in spin at all PN orders, note that the 1.5PN terms are already included. Here we have additional 2PN and 3.5PN corrections.
      pPrec->L4  += (chi1L2*(1 + delta - 2*eta) + 4*chi1L*chi2L*eta - chi2L2*(-1 + delta + 2*eta))/2.;
      pPrec->L7  +=  (3*(chi1L + chi2L)*eta*(chi1L2*(1 + delta - 2*eta) + 4*chi1L*chi2L*eta - chi2L2*(-1 + delta + 2*eta)))/4.;

      break;
    }

    default:
    {
      XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPrecVersion not recognized. Requires version 101, 102, 103, 104, 220, 221, 222, 223, 224, 310, 311, 320 or 321.\n");
      break;
    }
  }

  /* Reference orbital angular momentum */
  pPrec->LRef = M * M * XLALSimIMRPhenomXLPNAnsatz(pWF->v_ref, pWF->eta / pWF->v_ref, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);

  /*
    In the following code block we construct the convetions that relate the source frame and the LAL frame.

    A detailed discussion of the conventions can be found in Appendix C and D of arXiv:2004.06503 and https://dcc.ligo.org/LIGO-T1500602
  */

  /* Get source frame (*_Sf) J = L + S1 + S2. This is an instantaneous frame in which L is aligned with z */
  pPrec->J0x_Sf = (m1_2)*chi1x + (m2_2)*chi2x;
  pPrec->J0y_Sf = (m1_2)*chi1y + (m2_2)*chi2y;
  pPrec->J0z_Sf = (m1_2)*chi1z + (m2_2)*chi2z + pPrec->LRef;

  pPrec->J0     = sqrt(pPrec->J0x_Sf*pPrec->J0x_Sf + pPrec->J0y_Sf*pPrec->J0y_Sf + pPrec->J0z_Sf*pPrec->J0z_Sf);

  /* Get angle between J0 and LN (z-direction) */
  if(pPrec->J0 < 1e-10)
  {
    XLAL_PRINT_WARNING("Warning: |J0| < 1e-10. Setting thetaJ = 0.\n");
    pPrec->thetaJ_Sf = 0.0;
  }
  else
  {
    pPrec->thetaJ_Sf = acos(pPrec->J0z_Sf / pPrec->J0);
  }

  const double phiRef = pWF->phiRef_In;

  INT4 convention     = XLALSimInspiralWaveformParamsLookupPhenomXPConvention(lalParams);

  if ( !(convention == 0 || convention == 1 || convention == 5 || convention == 6 || convention == 7))
  {
    XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPConvention not recognized. Requires version 0, 1, 5, 6 or 7.\n");
  }

  #if DEBUG == 1
    printf("\n*** Convention = %i\n", convention);
  #endif

  /* Get azimuthal angle of J0 in the source frame */
  if(fabs(pPrec->J0x_Sf) < MAX_TOL_ATAN && fabs(pPrec->J0y_Sf) < MAX_TOL_ATAN)
  {
      #if DEBUG == 1
        printf("\nAligned spin limit!\n");
      #endif

      /* Impose the aligned spin limit */
      switch(convention)
      {
        case 0:
        case 5:
        {
          pPrec->phiJ_Sf = LAL_PI/2.0 - phiRef;
          break;
        }
        case 1:
        case 6:
        case 7:
        {
          pPrec->phiJ_Sf = 0;
          break;
        }

      }
  }
  else
  {
      pPrec->phiJ_Sf = atan2(pPrec->J0y_Sf, pPrec->J0x_Sf); /* azimuthal angle of J0 in the source frame */
  }
  pPrec->phi0_aligned = - pPrec->phiJ_Sf;

  switch(convention)
  {
    case 0:
    {
      pWF->phi0 = pPrec->phi0_aligned;
      break;
    }
    case 1:
    {
      pWF->phi0 = 0;
      break;
    }
    case 5:
    case 6:
    case 7:
    {
      break;
    }
  }

  /*
      Here we follow the same prescription as in IMRPhenomPv2:

      Now rotate from SF to J frame to compute alpha0, the azimuthal angle of LN, as well as
      thetaJ, the angle between J and N.

      The J frame is defined by imposing that J points in the z-direction and the line of sight N is in the xz-plane
      (with positive projection along x).

      The components of any vector in the (new) J-frame can be obtained by rotation from the (old) source frame (SF).
      This is done by multiplying by: RZ[-kappa].RY[-thetaJ].RZ[-phiJ]

      Note that kappa is determined by rotating N with RY[-thetaJ].RZ[-phiJ], which brings J to the z-axis, and
      taking the opposite of the azimuthal angle of the rotated N.
  */

  /* Determine kappa via rotations, as above */
  pPrec->Nx_Sf = sin(pWF->inclination)*cos((LAL_PI / 2.0) - phiRef);
  pPrec->Ny_Sf = sin(pWF->inclination)*sin((LAL_PI / 2.0) - phiRef);
  pPrec->Nz_Sf = cos(pWF->inclination);

  REAL8 tmp_x = pPrec->Nx_Sf;
  REAL8 tmp_y = pPrec->Ny_Sf;
  REAL8 tmp_z = pPrec->Nz_Sf;

  IMRPhenomX_rotate_z(-pPrec->phiJ_Sf,   &tmp_x, &tmp_y, &tmp_z);
  IMRPhenomX_rotate_y(-pPrec->thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);

  /* Note difference in overall - sign w.r.t PhenomPv2 code */
  pPrec->kappa = XLALSimIMRPhenomXatan2tol(tmp_y,tmp_x, MAX_TOL_ATAN);

  /* Now determine alpha0 by rotating LN. In the source frame, LN = {0,0,1} */
  tmp_x = 0.0;
  tmp_y = 0.0;
  tmp_z = 1.0;
  IMRPhenomX_rotate_z(-pPrec->phiJ_Sf,   &tmp_x, &tmp_y, &tmp_z);
  IMRPhenomX_rotate_y(-pPrec->thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);
  IMRPhenomX_rotate_z(-pPrec->kappa,     &tmp_x, &tmp_y, &tmp_z);

  if (fabs(tmp_x) < MAX_TOL_ATAN && fabs(tmp_y) < MAX_TOL_ATAN)
  {
      /* This is the aligned spin case */
      #if DEBUG == 1
        printf("\nAligned-spin case.\n");
      #endif

      switch(convention)
      {
        case 0:
        case 5:
        {
          pPrec->alpha0 = LAL_PI;
          break;
        }
        case 1:
        case 6:
        case 7:
        {
          pPrec->alpha0 = LAL_PI - pPrec->kappa;
          break;
        }
      }
  }
  else
  {
      switch(convention)
      {
        case 0:
        case 5:
        {
          pPrec->alpha0 = atan2(tmp_y,tmp_x);
          break;
        }
        case 1:
        case 6:
        case 7:
        {
          pPrec->alpha0 = LAL_PI - pPrec->kappa;
          break;
        }
      }
  }


  switch(convention)
  {
    case 0:
    case 5:
    {
        /* Now determine thetaJN by rotating N */
        tmp_x = pPrec->Nx_Sf;
        tmp_y = pPrec->Ny_Sf;
        tmp_z = pPrec->Nz_Sf;
        IMRPhenomX_rotate_z(-pPrec->phiJ_Sf,   &tmp_x, &tmp_y, &tmp_z);
        IMRPhenomX_rotate_y(-pPrec->thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);
        IMRPhenomX_rotate_z(-pPrec->kappa,     &tmp_x, &tmp_y, &tmp_z);

        /* We don't need the y-component but we will store it anyway */
        pPrec->Nx_Jf = tmp_x;
        pPrec->Ny_Jf = tmp_y;
        pPrec->Nz_Jf = tmp_z;

        /* This is a unit vector, so no normalization */
        pPrec->thetaJN = acos(pPrec->Nz_Jf);
        break;
    }
    case 1:
    case 6:
    case 7:
    {
        REAL8 J0dotN     = (pPrec->J0x_Sf * pPrec->Nx_Sf) + (pPrec->J0y_Sf * pPrec->Ny_Sf) + (pPrec->J0z_Sf * pPrec->Nz_Sf);
        pPrec->thetaJN   = acos( J0dotN / pPrec->J0 );
        pPrec->Nz_Jf     = cos(pPrec->thetaJN);
        pPrec->Nx_Jf     = sin(pPrec->thetaJN);
        break;
    }
  }


  /*
      Define the polarizations used. This follows the conventions adopted for IMRPhenomPv2.

      The IMRPhenomP polarizations are defined following the conventions in Arun et al (arXiv:0810.5336),
      i.e. projecting the metric onto the P, Q, N triad defining where: P = (N x J) / |N x J|.

      However, the triad X,Y,N used in LAL (the "waveframe") follows the definition in the
      NR Injection Infrastructure (Schmidt et al, arXiv:1703.01076).

      The triads differ from each other by a rotation around N by an angle \zeta. We therefore need to rotate
      the polarizations by an angle 2 \zeta.
  */
  pPrec->Xx_Sf = -cos(pWF->inclination) * sin(phiRef);
  pPrec->Xy_Sf = -cos(pWF->inclination) * cos(phiRef);
  pPrec->Xz_Sf = +sin(pWF->inclination);

  tmp_x = pPrec->Xx_Sf;
  tmp_y = pPrec->Xy_Sf;
  tmp_z = pPrec->Xz_Sf;

  IMRPhenomX_rotate_z(-pPrec->phiJ_Sf,   &tmp_x, &tmp_y, &tmp_z);
  IMRPhenomX_rotate_y(-pPrec->thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);
  IMRPhenomX_rotate_z(-pPrec->kappa,     &tmp_x, &tmp_y, &tmp_z);


  /*
      The components tmp_i are now the components of X in the J frame.

      We now need the polar angle of this vector in the P, Q basis of Arun et al:

          P = (N x J) / |NxJ|

      Note, that we put N in the (pos x)z half plane of the J frame
  */

  switch(convention)
  {
    case 0:
    case 5:
    {
      /* Get polar angle of X vector in J frame in the P,Q basis of Arun et al */
      pPrec->PArunx_Jf = +0.0;
      pPrec->PAruny_Jf = -1.0;
      pPrec->PArunz_Jf = +0.0;

      /* Q = (N x P) by construction */
      pPrec->QArunx_Jf =  pPrec->Nz_Jf;
      pPrec->QAruny_Jf =  0.0;
      pPrec->QArunz_Jf = -pPrec->Nx_Jf;
      break;
    }
    case 1:
    case 6:
    case 7:
    {
      /* Get polar angle of X vector in J frame in the P,Q basis of Arun et al */
      pPrec->PArunx_Jf = pPrec->Nz_Jf;
      pPrec->PAruny_Jf = 0;
      pPrec->PArunz_Jf = -pPrec->Nx_Jf;

      /* Q = (N x P) by construction */
      pPrec->QArunx_Jf =  0;
      pPrec->QAruny_Jf =  1;
      pPrec->QArunz_Jf =  0;
      break;
    }
  }

  // (X . P)
  pPrec->XdotPArun = (tmp_x * pPrec->PArunx_Jf) + (tmp_y * pPrec->PAruny_Jf) + (tmp_z * pPrec->PArunz_Jf);

  // (X . Q)
  pPrec->XdotQArun = (tmp_x * pPrec->QArunx_Jf) + (tmp_y * pPrec->QAruny_Jf) + (tmp_z * pPrec->QArunz_Jf);

  /* Now get the angle zeta */
  pPrec->zeta_polarization = atan2(pPrec->XdotQArun, pPrec->XdotPArun);

  /* ********** PN Euler Angle Coefficients ********** */
  /*
      This uses the single spin PN Euler angles as per IMRPhenomPv2
  */

  /* ********** PN Euler Angle Coefficients ********** */
  switch( pflag )
  {
    case 101:
    case 102:
    case 103:
    case 104:
    {
      /*
          This uses the single spin PN Euler angles as per IMRPhenomPv2
      */

      /* Post-Newtonian Euler Angles: alpha */
      REAL8 chiL       = (1.0 + q) * (chi_eff / q);
      REAL8 chiL2      = chiL * chiL;

      pPrec->alpha1    = -35/192. + (5*delta)/(64.*m1);

      pPrec->alpha2    = ((15*chiL*delta*m1)/128. - (35*chiL*m1_2)/128.)/eta;

      pPrec->alpha3    = -5515/3072. + eta*(-515/384. - (15*delta2)/(256.*m1_2)
                          + (175*delta)/(256.*m1)) + (4555*delta)/(7168.*m1)
                          + ((15*chip2*delta*m1_3)/128. - (35*chip2*m1_4)/128.)/eta2;

      /* This is the term proportional to log(w) */
      pPrec->alpha4L   = ((5*chiL*delta2)/16. - (5*chiL*delta*m1)/3. + (2545*chiL*m1_2)/1152.
                          + ((-2035*chiL*delta*m1)/21504.
                          + (2995*chiL*m1_2)/9216.)/eta + ((5*chiL*chip2*delta*m1_5)/128.
                          - (35*chiL*chip2*m1_6)/384.)/eta3
                          - (35*LAL_PI)/48. + (5*delta*LAL_PI)/(16.*m1));

      pPrec->alpha5    = (5*(-190512*delta3*eta6 + 2268*delta2*eta3*m1*(eta2*(323 + 784*eta)
                          + 336*(25*chiL2 + chip2)*m1_4) + 7*m1_3*(8024297*eta4 + 857412*eta5
                          + 3080448*eta6 + 143640*chip2*eta2*m1_4
                          - 127008*chip2*(-4*chiL2 + chip2)*m1_8
                          + 6048*eta3*((2632*chiL2 + 115*chip2)*m1_4 - 672*chiL*m1_2*LAL_PI))
                          + 3*delta*m1_2*(-5579177*eta4 + 80136*eta5 - 3845520*eta6
                          + 146664*chip2*eta2*m1_4 + 127008*chip2*(-4*chiL2 + chip2)*m1_8
                          - 42336*eta3*((726*chiL2 + 29*chip2)*m1_4
                          - 96*chiL*m1_2*LAL_PI))))/(6.5028096e7*eta4*m1_3);

      /* Post-Newtonian Euler Angles: epsilon */
      pPrec->epsilon1  = -35/192. + (5*delta)/(64.*m1);

      pPrec->epsilon2  = ((15*chiL*delta*m1)/128. - (35*chiL*m1_2)/128.)/eta;

      pPrec->epsilon3  = -5515/3072. + eta*(-515/384. - (15*delta2)/(256.*m1_2)
                          + (175*delta)/(256.*m1)) + (4555*delta)/(7168.*m1);

      /* This term is proportional to log(w) */
      pPrec->epsilon4L = (5*chiL*delta2)/16. - (5*chiL*delta*m1)/3. + (2545*chiL*m1_2)/1152.
                          + ((-2035*chiL*delta*m1)/21504. + (2995*chiL*m1_2)/9216.)/eta - (35*LAL_PI)/48.
                          + (5*delta*LAL_PI)/(16.*m1);

      pPrec->epsilon5  = (5*(-190512*delta3*eta3 + 2268*delta2*m1*(eta2*(323 + 784*eta)
                        + 8400*chiL2*m1_4)
                        - 3*delta*m1_2*(eta*(5579177 + 504*eta*(-159 + 7630*eta))
                        + 254016*chiL*m1_2*(121*chiL*m1_2 - 16*LAL_PI))
                        + 7*m1_3*(eta*(8024297 + 36*eta*(23817 + 85568*eta))
                        + 338688*chiL*m1_2*(47*chiL*m1_2 - 12*LAL_PI))))/(6.5028096e7*eta*m1_3);

      break;
    }
    case 220:
    case 221:
    case 222:
    case 223:
    case 224:
    case 310:
    case 311:
    case 320:
    case 321:
    {
      pPrec->alpha1    = 0;
      pPrec->alpha2    = 0;
      pPrec->alpha3    = 0;
      pPrec->alpha4L   = 0;
      pPrec->alpha5    = 0;
      pPrec->epsilon1  = 0;
      pPrec->epsilon2  = 0;
      pPrec->epsilon3  = 0;
      pPrec->epsilon4L = 0;
      pPrec->epsilon5  = 0;
      break;
    }
    default:
    {
      XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPrecVersion not recognized. Requires version 101, 102, 103, 104, 220, 221, 222, 223, 224, 310, 311, 320 or 321.\n");
      break;
    }
  }

  REAL8 alpha_offset = 0, epsilon_offset = 0;

  #if DEBUG == 1
      printf("thetaJN             : %e\n",   pPrec->thetaJN);
      printf("phiJ_Sf             : %e\n", pPrec->phiJ_Sf);
      printf("alpha0              : %e\n", pPrec->alpha0);
      printf("pi-kappa            : %e\n", LAL_PI-pPrec->kappa);
      printf("kappa               : %e\n", pPrec->kappa);
      printf("pi/2 - phiRef       : %e\n", LAL_PI_2 - phiRef);
      printf("zeta_polarization   : %.16e\n", pPrec->zeta_polarization);
      printf("zeta_polarization   : %.16e\n", acos(pPrec->XdotPArun));
      printf("zeta_polarization   : %.16e\n", asin(pPrec->XdotQArun));
      printf("zeta_polarization   : %.16e\n\n", LAL_PI_2 - acos(pPrec->XdotQArun));
      printf("alpha1              : %e\n",  pPrec->alpha1);
      printf("alpha2              : %e\n",  pPrec->alpha2);
      printf("alpha3              : %e\n",  pPrec->alpha3);
      printf("alpha4L             : %e\n",  pPrec->alpha4L);
      printf("alpha5              : %e\n\n",  pPrec->alpha5);
  #endif


  switch(convention)
  {
    case 0:
      pPrec->epsilon0 = 0;
      break;
    case 1:
    case 6:
      pPrec->epsilon0 = pPrec->phiJ_Sf - LAL_PI;
      break;
    case 5:
    case 7:
      pPrec->epsilon0 = 0;
      break;
  }

  if(convention == 5 || convention == 7)
  {
    pPrec->alpha_offset = -pPrec->alpha0;
    pPrec->epsilon_offset = 0;
    pPrec->alpha_offset_1 = -pPrec->alpha0;
    pPrec->epsilon_offset_1 = 0;
    pPrec->alpha_offset_3 = -pPrec->alpha0;
    pPrec->epsilon_offset_3 = 0;
    pPrec->alpha_offset_4 = -pPrec->alpha0;
    pPrec->epsilon_offset_4 = 0;
  }
  else
  {
    /* Get initial Get \alpha and \epsilon offsets at \omega = pi * M * f_{Ref} */
    Get_alphaepsilon_atfref(&alpha_offset, &epsilon_offset, 2, pPrec, pWF);
    pPrec->alpha_offset       = alpha_offset;
    pPrec->epsilon_offset     = epsilon_offset;
    pPrec->alpha_offset_1     = alpha_offset;
    pPrec->epsilon_offset_1   = epsilon_offset;
    pPrec->alpha_offset_3     = alpha_offset;
    pPrec->epsilon_offset_3   = epsilon_offset;
    pPrec->alpha_offset_4     = alpha_offset;
    pPrec->epsilon_offset_4   = epsilon_offset;
  }

  pPrec->cexp_i_alpha   = 0.;
  pPrec->cexp_i_epsilon = 0.;
  pPrec->cexp_i_betah   = 0.;

  /* Activate multibanding for Euler angles it threshold !=0. Only for PhenomXPHM. */
  if(XLALSimInspiralWaveformParamsLookupPhenomXPHMThresholdMband(lalParams)==0)
  {
    /* User switched off multibanding */
    pPrec->MBandPrecVersion = 0;
  }
  else
  {
    /* User requested multibanding */
    pPrec->MBandPrecVersion = 1;

    /* Switch off multiband for very high mass as in IMRPhenomXHM. */
    if(pWF->Mtot > 500)
    {
      XLAL_PRINT_WARNING("Very high mass, only merger in frequency band, multibanding not efficient, switching off for non-precessing modes and Euler angles.");
      pPrec->MBandPrecVersion = 0;
      XLALSimInspiralWaveformParamsInsertPhenomXHMThresholdMband(lalParams, 0.);
    }

    if(pPrec->IMRPhenomXPrecVersion < 200)
    {
      /* The NNLO angles can have a worse, even pathological, behaviour for high mass ratio and double spin cases.
       The waveform will look noisy, we switch off the multibanding for mass ratio above 8 to avoid worsen even more the waveform. */
      if(pWF->q > 8)
      {
        XLAL_PRINT_WARNING("Very high mass ratio, NNLO angles may become pathological, switching off multibanding for angles.\n");
        XLALSimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(lalParams, 0.);
        pPrec->MBandPrecVersion = 0;
      }
    }
    /* The MSA angles give quite 'noisy' waveforms in this corner of parameter space so we switch off multibanding to avoid worsen the waveform. */
    else if ( pWF->q > 50 && pWF->Mtot > 100 )
    {
      XLALSimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(lalParams, 0.);
      pPrec->MBandPrecVersion = 0;
    }

  }


  const REAL8 ytheta  = pPrec->thetaJN;
  const REAL8 yphi    = 0.0;
  pPrec->Y2m2         = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, -2);
  pPrec->Y2m1         = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2, -1);
  pPrec->Y20          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2,  0);
  pPrec->Y21          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2,  1);
  pPrec->Y22          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 2,  2);
  pPrec->Y3m3         = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, -3);
  pPrec->Y3m2         = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, -2);
  pPrec->Y3m1         = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3, -1);
  pPrec->Y30          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3,  0);
  pPrec->Y31          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3,  1);
  pPrec->Y32          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3,  2);
  pPrec->Y33          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 3,  3);
  pPrec->Y4m4         = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, -4);
  pPrec->Y4m3         = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, -3);
  pPrec->Y4m2         = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, -2);
  pPrec->Y4m1         = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4, -1);
  pPrec->Y40          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4,  0);
  pPrec->Y41          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4,  1);
  pPrec->Y42          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4,  2);
  pPrec->Y43          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4,  3);
  pPrec->Y44          = XLALSpinWeightedSphericalHarmonic(ytheta, yphi, -2, 4,  4);

  /*
      Check whether maximum opening angle becomes larger than \pi/2 or \pi/4.

      If (L + S_L) < 0, then Wigner-d Coefficients will not track the angle between J and L, meaning
      that the model may become pathological as one moves away from the aligned-spin limit.

      If this does not happen, then max_beta will be the actual maximum opening angle.

      This function uses a 2PN non-spinning approximation to the orbital angular momentum L, as
      the roots can be analytically derived.

      Returns XLAL_PRINT_WARNING if model is in a pathological regime.
  */
  IMRPhenomXPCheckMaxOpeningAngle(pWF,pPrec);

  return XLAL_SUCCESS;
}
/* Function to set remnant quantities related to final
spin within IMRPhenomXGetAndSetPrecessionVariables */
INT4 IMRPhenomX_SetPrecessingRemnantParams(
  IMRPhenomXWaveformStruct *pWF,
  IMRPhenomXPrecessionStruct *pPrec,
  LALDict *lalParams
){

  /*

  The strategy for setting the precessing remnant spin in PhenomX
  is multifaceted becuase there are not only the various options
  implemented for the original PhenomXP, but there are also the
  options needed for PNR's construction, and its turning off of
  its CoPrecessing model outside of the PNR calibration region.
  In that latter case, the final spin transitions from the
  non-precessing final spin, where PNR's CoPrecessing model is
  tuned, to the precessing final spin, as is needed for the EZH
  effective ringdown result.

  A layer on top of this, is the need for the l=m=2 fundamental
  QNM frequency to be computed, in some way, amid all of these
  scenarios.

  This function exists to draw a conceptual circle around all of
  this mess.

  Most comments below are cogent, but some have been left for
  historical record?

  */

  /* Define shorthand variables for PNR CoPrec options */
  INT4 status = 0;

  // Toggle for PNR coprecessing tuning
  INT4 PNRUseInputCoprecDeviations = pPrec->IMRPhenomXPNRUseInputCoprecDeviations;

  // Toggle for enforced use of non-precessing spin as is required during tuning of PNR's coprecessing model
  INT4 PNRUseTunedCoprec = pPrec->IMRPhenomXPNRUseTunedCoprec;

  // High-level toggle for whether to apply deviations
  INT4 APPLY_PNR_DEVIATIONS = pWF->APPLY_PNR_DEVIATIONS;

  /*
  HISTOICAL COMMENT
      Update the final spin in pWF to account for precessing spin effects:

      XLALSimIMRPhenomXPrecessingFinalSpin2017(eta,chi1L,chi2L,chi_perp);

      Note that chi_perp gets weighted by an appropriate mass factor

        q_factor      = m1/M (when m1 > m2)
        REAL8 Sperp   = chip * q_factor * q_factor;
        REAL8 af      = copysign(1.0, af_parallel) * sqrt(Sperp*Sperp + af_parallel*af_parallel);
  */
  REAL8 M = pWF->M;
  REAL8 af_parallel = XLALSimIMRPhenomXFinalSpin2017(pWF->eta,pPrec->chi1z,pPrec->chi2z);
  double Lfinal     = M*M*af_parallel - pWF->m1_2*pPrec->chi1z - pWF->m2_2*pPrec->chi2z;

  int fsflag = XLALSimInspiralWaveformParamsLookupPhenomXPFinalSpinMod(lalParams);
  if (fsflag == 4 && pPrec->precessing_tag!=3) fsflag = 3;

  /* For PhenomPNR, we wil use the PhenomPv2 final spin function's result, modified such that its sign is given by sign( cos(betaRD) ). See the related fsflag case below. */
  if (PNRUseTunedCoprec) fsflag = 5;

  /* When tuning the coprecessing model, we wish to enforce use of the non-precessing final spin. See the related fsflag case below. */
  if (PNRUseInputCoprecDeviations) fsflag = 6;

  /* Generate and store ringdown value of precession angle beta. This is to be used for e.g. setting the sign of the final spin, and calculating the effective ringdown frequency */
  if( PNRUseTunedCoprec ){
    pWF->betaRD = IMRPhenomX_PNR_GenerateRingdownPNRBeta( pWF, pPrec);
  }

  //
  double chi1L = pPrec->chi1z;
  double chi2L = pPrec->chi2z;

  //
  int pflag = pPrec->IMRPhenomXPrecVersion;

  //
  switch(fsflag)
  {
    case 0:
      pWF->afinal_prec    = XLALSimIMRPhenomXPrecessingFinalSpin2017(pWF->eta,chi1L,chi2L,pPrec->chi_p);
      break;
    case 1:
      pWF->afinal_prec    = XLALSimIMRPhenomXPrecessingFinalSpin2017(pWF->eta,chi1L,chi2L,pPrec->chi1x);
      break;
    case 2:
    case 4:
      pWF->afinal_prec    = XLALSimIMRPhenomXPrecessingFinalSpin2017(pWF->eta,chi1L,chi2L,pPrec->chiTot_perp);
      break;
    case 3:
      if( pflag == 220 || pflag == 221 || pflag == 222 || pflag == 223 || pflag == 224 )
      {
        if(pPrec->MSA_ERROR == 1 )
        {
          XLAL_PRINT_WARNING("Initialization of MSA system failed. Defaulting to final spin version 0.\n");
          pWF->afinal_prec    = XLALSimIMRPhenomXPrecessingFinalSpin2017(pWF->eta,chi1L,chi2L,pPrec->chi_p);
        }
        else
        {

          INT2 sign = 1;
          if (XLALSimInspiralWaveformParamsLookupPhenomXPTransPrecessionMethod(lalParams) == 1 ){
            sign = copysign(1, af_parallel);
          }
          pWF->afinal_prec    = sign * sqrt( pPrec->SAv2 + Lfinal*Lfinal + 2.0*Lfinal*(pPrec->S1L_pav + pPrec->S2L_pav) ) / (M*M);

        }
      }
      else
      {
        XLAL_PRINT_WARNING("Error: XLALSimInspiralWaveformParamsLookupPhenomXPFinalSpinMod version 3 requires PrecVersion 220, 221, 222, 223 or 224. Defaulting to version 0.\n");
        pWF->afinal_prec    = XLALSimIMRPhenomXPrecessingFinalSpin2017(pWF->eta,chi1L,chi2L,pPrec->chi_p);
      }
      break;

    case 5:
      {
      /*-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~*
      Implement Pv2 final spin but with sign derived from EZH's model for ringdown beta.
      *-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~*/

      // Use these as input into method for effective RD frequency
      // * get the value of ringdown beta
      INT2 sign = 1;
      sign = copysign(1, cos(pWF->betaRD) );

      /* Calculate Pv2 final spin without alteration. Below we alter it. NOTE that XLALSimIMRPhenomXPrecessingFinalSpin2017 appears to be an explicit copy of the actual Pv2 final spin function, FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH. The original PhenomX have not referenced this code duplication.  */
      double afinal_prec_Pv2 = XLALSimIMRPhenomXPrecessingFinalSpin2017(pWF->eta,chi1L,chi2L,pPrec->chi_p);

      // The eqiuvalent PhenomPv2 code reference would look like the commented line below.
      // double afinal_prec_Pv2 = FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH( m1, m2, chi1L, chi2L, chi_p );

      // Define the PNR final spin to be the Pv2 final spin magnitude, with direction given by sign of cos betaRD
      pWF->afinal_prec = sign * fabs(afinal_prec_Pv2);

      // // Experimental version of final spin
      // pWF->afinal_prec = cos(pWF->betaRD) * fabs(afinal_prec_Pv2);
      }
      break;

    case 6:
      {
        /*-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~*
         During PNR tuning, we wish to evaluate the coprecessing model with the same final spin that would be used in PhenomXHM. We implement this here.
         *-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~*/
        pWF->afinal_prec = pWF->afinal_nonprec;
      }
      break;
    default:
    {
      XLAL_ERROR(XLAL_EDOM,"Error: XLALSimInspiralWaveformParamsLookupPhenomXPFinalSpinMod version not recognized. Requires PhenomXPFinalSpinMod of 0, 1, 2, 3, or 5.\n");
    }
  }

  /* (PNRUseTunedCoprec) When not generating PNR make NO destinction between afinal and afinal_prec */
  if( !PNRUseTunedCoprec ){
    pWF->afinal = pWF->afinal_prec;
  } else {
    /*  ELSE, use the non-precessing final spin defined in IMRPhenomXSetWaveformVariables. XCP uses the non-precessing parameters as a base upon which to add precessing deviations. The line below is added only for clarity: pWF->afinal is already equal to pWF->afinal_nonprec as is set in IMRPhenomXSetWaveformVariables */

    /* pWF->afinal = pWF->afinal_nonprec */
    // ABOVE but commented out, we see what the final spin assignment WOULD BE if NO WINDOWING of coprec tuning were used

    pWF->afinal = (pWF->pnr_window)*pWF->afinal_nonprec + (1.0-pWF->pnr_window)*pWF->afinal_prec;
    // Above: NOTE that as PNR is turned off outside of its calibration window, we want to turn on use of the precessing final spin as defined in the code section above

    #if DEBUG == 1
      printf("*** PNR Co-precessing model in use (l=m=2) ***\n\n");
      printf("PNR window   : %e\n", pWF->pnr_window);
      printf("pWF->afinal      : %e\n",pWF->afinal);
      printf("pWF->afinal_prec  : %e\n",pWF->afinal_prec);
      printf("pWF->afinal_nonprec  : %e\n",pWF->afinal_nonprec);
    #endif
  }


  if( fabs(pWF->afinal) > 1.0 )
  {
        XLAL_PRINT_WARNING("Warning: Final spin magnitude %g > 1. Setting final spin magnitude = 1.", pWF->afinal);
        pWF->afinal = copysign(1.0, pWF->afinal);
  }

  /* Update ringdown and damping frequency: no precession; to be used for PNR tuned deviations */
  pWF->fRING     = evaluate_QNMfit_fring22(pWF->afinal) / (pWF->Mfinal);
  pWF->fDAMP     = evaluate_QNMfit_fdamp22(pWF->afinal) / (pWF->Mfinal);
  //pWF->fISCO     = XLALSimIMRPhenomXfISCO(pWF->afinal);

  #if DEBUG == 1
    printf("afinal (prec)  : %e\n",pWF->afinal);
    printf("fring  (prec)  : %e\n",pWF->fRING);
    printf("fdamp  (prec)  : %e\n\n",pWF->fDAMP);
  #endif

  // Copy IMRPhenomXReturnCoPrec to pWF
  pWF->IMRPhenomXReturnCoPrec = pPrec->IMRPhenomXReturnCoPrec;
  #if DEBUG == 1
    printf("pPrec->IMRPhenomXReturnCoPrec : %i\n",pPrec->IMRPhenomXReturnCoPrec);
    printf("pWF->IMRPhenomXReturnCoPrec   : %i\n",pWF->IMRPhenomXReturnCoPrec);
  #endif

  //
  if( APPLY_PNR_DEVIATIONS )
  {
    /* Add an overall deviation to the high-level ringdown frequency (PNRUseTunedCoprec) */
    pWF->fRING = pWF->fRING - (pWF->PNR_DEV_PARAMETER * pWF->NU5);
    pWF->fDAMP = pWF->fDAMP + (pWF->PNR_DEV_PARAMETER * pWF->NU6);
  }

  // we want to define the quantities below if PNR is used (PNRUseTunedCoprec). In particular,  pWF->fRINGEffShiftDividedByEmm is used by the HMs
  // Define identifiers for perturbation theory frequencies
  if( PNRUseTunedCoprec && (pWF->PNR_SINGLE_SPIN != 1))
  {
    const REAL8 fRING22_prec = evaluate_QNMfit_fring22(pWF->afinal_prec) / (pWF->Mfinal);
    const REAL8 fRING21_prec = evaluate_QNMfit_fring21(pWF->afinal_prec) / (pWF->Mfinal);
    pWF->fRING22_prec = fRING22_prec;

    // * Calculate and store single quantity needed to determine effective ringdown frequencies for all QNMs
    pWF->fRINGEffShiftDividedByEmm = (1.0-fabs(cos(pWF->betaRD))) * ( fRING22_prec  -  fRING21_prec );


    // As we turn off PNR tuning, we want to turn on use of the effective ringdown frequency
    // NOTE that when pWF->pnr_window=0, this should reduce to the def of pWF->fRINGCP below

    const INT4 emm = 2;
    // NOTE that we use pWF->fRING and not fRING22_prec below because of how pWF->afinal is defined using (1-pWF->pnr_window)
    pWF->fRING = pWF->fRING  -  (1.0-pWF->pnr_window) * emm * pWF->fRINGEffShiftDividedByEmm;
    // pWF->fRING = (pWF->pnr_window)*pWF->fRING  -  (1-pWF->pnr_window) * emm * pWF->fRINGEffShiftDividedByEmm;

    #if DEBUG == 1
      printf("pflag                         : %i\n",pflag);
      printf("pWF->betaRD                   : %e\n",pWF->betaRD);
      printf("pWF->fRINGEffShiftDividedByEm : %e\n", pWF->fRINGEffShiftDividedByEmm);
      printf("fring22 (prec)                : %e\n",fRING22_prec);
      printf("fRING                         : %e\n",pWF->fRING);
    #endif

  }

  //
  return status;

};

/** Get alpha and epsilon offset depending of the mprime (second index of the non-precessing mode) */
void Get_alphaepsilon_atfref(REAL8 *alpha_offset, REAL8 *epsilon_offset, UINT4 mprime, IMRPhenomXPrecessionStruct *pPrec, IMRPhenomXWaveformStruct *pWF)
{
  /* Compute the offsets due to the choice of integration constant in alpha and epsilon PN formula */
  double omega_ref       = pWF->piM * pWF->fRef * 2./mprime;

  int pflag = pPrec->IMRPhenomXPrecVersion;

  /* Explicitly enumerate MSA flags */
  if(pflag == 220 || pflag == 221 || pflag == 222 || pflag == 223 || pflag == 224)
  {
    const double v        = cbrt ( omega_ref );
    const vector vangles  = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(v,pWF,pPrec);

    *alpha_offset    = vangles.x - pPrec->alpha0;
    *epsilon_offset  = vangles.y - pPrec->epsilon0;
  }

  else
  {
    double logomega_ref    = log(omega_ref);
    double omega_ref_cbrt  = cbrt(omega_ref);
    double omega_ref_cbrt2 = omega_ref_cbrt * omega_ref_cbrt;

    *alpha_offset = (
        pPrec->alpha1  / omega_ref
      + pPrec->alpha2  / omega_ref_cbrt2
      + pPrec->alpha3  / omega_ref_cbrt
      + pPrec->alpha4L * logomega_ref
      + pPrec->alpha5  * omega_ref_cbrt - pPrec->alpha0
    );

    *epsilon_offset =  (
          pPrec->epsilon1  / omega_ref
        + pPrec->epsilon2  / omega_ref_cbrt2
        + pPrec->epsilon3  / omega_ref_cbrt
        + pPrec->epsilon4L * logomega_ref
        + pPrec->epsilon5  * omega_ref_cbrt - pPrec->epsilon0
    );
  }

  return;
};


/**
  This is a convenient wrapper function for PN orbital angular momentum.
*/
REAL8 XLALSimIMRPhenomXLPNAnsatz(
  REAL8 v,        /**< Input velocity  */
  REAL8 LNorm,    /**< Orbital angular momentum normalization */
  REAL8 L0,       /**< Newtonian orbital angular momentum (i.e. LN = 1.0*LNorm) */
  REAL8 L1,       /**< 0.5PN Orbital angular momentum */
  REAL8 L2,       /**< 1.0PN Orbital angular momentum */
  REAL8 L3,       /**< 1.5PN Orbital angular momentum */
  REAL8 L4,       /**< 2.0PN Orbital angular momentum */
  REAL8 L5,       /**< 2.5PN Orbital angular momentum */
  REAL8 L6,       /**< 3.0PN Orbital angular momentum */
  REAL8 L7,       /**< 3.5PN Orbital angular momentum */
  REAL8 L8,       /**< 4.0PN Orbital angular momentum */
  REAL8 L8L       /**< 4.0PN logarithmic orbital angular momentum term */
)
{
  const REAL8 x   = v*v;
  const REAL8 x2  = x*x;
  const REAL8 x3  = x*x2;
  const REAL8 x4  = x*x3;
  const REAL8 sqx = sqrt(x);

  /*
      Here LN is the Newtonian pre-factor: LN = \eta / \sqrt{x} :

      L = L_N \sum_a L_a x^{a/2}
        = L_N [ L0 + L1 x^{1/2} + L2 x^{2/2} + L3 x^{3/2} + ... ]

  */
  return LNorm * (L0 + L1*sqx + L2*x + L3*(x*sqx) + L4*x2 + L5*(x2*sqx) + L6*x3 + L7*(x3*sqx) + L8*x4 + L8L*x4*log(x));
}


/**
    2PN non-spinning orbital angular momentum as a function of x = v^2 = (Pi M f)^{2/3}

    - Bohe et al, 1212.5520v2, Eq 4.7
*/
REAL8 XLALSimIMRPhenomXL2PNNS(const REAL8 v, const REAL8 eta)
{
  const REAL8 eta2  = eta*eta;
  const REAL8 x     = v*v;
  const REAL8 x2    = x*x;
  const REAL8 sqx   = v;

  return (eta / sqx) * ( 1.0 + x * (3/2. + eta/6.) + x2 * (27/8. - (19*eta)/8. + eta2/24.) );
}

/**
    3PN orbital angular momentum as a function of x = v^2 = (Pi M f)^{2/3}

    Includes linear in spin corrections up to 3.5PN

    - Bohe et al, 1212.5520v2, Eq 4.7
*/
REAL8 XLALSimIMRPhenomXL3PNAS(const REAL8 v, const REAL8 eta, const REAL8 chi1L, const REAL8 chi2L, const REAL8 delta)
{
  const REAL8 eta2   = eta*eta;

  const REAL8 x      = v*v;
  const REAL8 x2     = x*x;
  const REAL8 x3     = x2*x;
  const REAL8 sqx    = v;

  return (eta/sqx) * (
    1.0   // Newtonian
    + (3/2. + eta/6.) * x // 1PN
    + (chi1L*(-5/3. - (5*delta)/3. + (5*eta)/6.) + chi2L*(-5/3. + (5*delta)/3. + (5*eta)/6.)) * x * sqx // 1.5PN
    + (27/8. - (19*eta)/8. + eta2/24.) * x2 // 2PN
    + ((-7*(chi1L*(72 + delta*(72 - 31*eta) + eta*(-121 + 2*eta)) + chi2L*(72 + eta*(-121 + 2*eta) + delta*(-72 + 31*eta))))/144.) * x2 * sqx // 2.5PN
    + ((10935 + eta*(-62001 + eta*(1674 + 7*eta) + 2214*LAL_TWOPI))/1296.) * x3 // 3PN
    );
}

/**
    4PN orbital angular momentum as a function of x = v^2 = (Pi M f)^{2/3}

    - Bohe et al, 1212.5520v2, Eq 4.7
    - Marsat, CQG, 32, 085008, (2015), arXiv:1411.4118
    - Siemonsen et al, PRD, 97, 064010, (2018), arXiv:1606.08832
*/
REAL8 XLALSimIMRPhenomXL4PNAS(const REAL8 v, const REAL8 eta, const REAL8 chi1L, const REAL8 chi2L, const REAL8 delta)
{
  const REAL8 x      = v*v;
  const REAL8 x2     = x*x;
  const REAL8 x3     = x2*x;
  const REAL8 x4     = x3*x;
  const REAL8 sqx    = sqrt(x);
  const REAL8 logx   = log(x);

  return (eta/sqx) * (
            1.0 //Newtonian
            + ((9 + eta)/6.) * x // 1PN
            + ((5*chi1L*(-2 - 2*delta + eta) + 5*chi2L*(-2 + 2*delta + eta))/6.) * x * sqx // 1.5PN
            + ((81 + (-57 + eta)*eta)/24.) * x2 // 2PN
            + ((-7*(chi1L*(72 + delta*(72 - 31*eta) + eta*(-121 + 2*eta)) + chi2L*(72 + eta*(-121 + 2*eta) + delta*(-72 + 31*eta))))/144.) * x2 * sqx //2.5PN
            + ((10935 + eta*(-62001 + eta*(1674 + 7*eta) + 2214*LAL_TWOPI))/1296.) * x3 //3PN
            + ((chi2L*(-324 + eta*(1119 - 2*eta*(172 + eta)) + delta*(324 + eta*(-633 + 14*eta)))
                      - chi1L*(324 + eta*(-1119 + 2*eta*(172 + eta)) + delta*(324 + eta*(-633 + 14*eta))))/32.) * x3 * sqx //3.5PN
            + (2835/128. - (eta*(-10677852 + 100*eta*(-640863 + eta*(774 + 11*eta)) + 26542080*LAL_GAMMA
                        + 675*(3873 + 3608*eta)*LAL_TWOPI))/622080. - (64*eta*log(16.))/3.) * x4 //4PN
            + (-64*eta/3.) * x4 * logx // 4PN log[x] term
    );
}

/**
    4PN orbital angular momentum as a function of x = v^2 = (Pi M f)^{2/3}

    - Bohe et al, 1212.5520v2, Eq 4.7
    - Marsat, CQG, 32, 085008, (2015), arXiv:1411.4118
    - Siemonsen et al, PRD, 97, 064010, (2018), arXiv:1606.08832
*/
REAL8 XLALSimIMRPhenomXL4PNLOSIAS(const REAL8 v, const REAL8 eta, const REAL8 chi1L, const REAL8 chi2L, const REAL8 delta)
{
  const REAL8 chi1L2 = chi1L * chi1L;
  const REAL8 chi2L2 = chi2L * chi2L;

  const REAL8 x      = v*v;
  const REAL8 x2     = x*x;
  const REAL8 x3     = x2*x;
  const REAL8 x4     = x3*x;
  const REAL8 sqx    = sqrt(x);
  const REAL8 logx   = log(x);

  return (eta/sqx) * (
            1.0 //Newtonian
            + ((9 + eta)/6.) * x // 1PN
            + ((5*chi1L*(-2 - 2*delta + eta) + 5*chi2L*(-2 + 2*delta + eta))/6.) * x * sqx // 1.5PN
            + ((81 + (-57 + eta)*eta)/24.) * x2 // 2PN
            + ((-7*(chi1L*(72 + delta*(72 - 31*eta) + eta*(-121 + 2*eta)) + chi2L*(72 + eta*(-121 + 2*eta) + delta*(-72 + 31*eta))))/144.) * x2 * sqx //2.5PN
            + ((10935 + eta*(-62001 + eta*(1674 + 7*eta) + 2214*LAL_TWOPI))/1296.) * x3 //3PN
            + ((chi2L*(-324 + eta*(1119 - 2*eta*(172 + eta)) + delta*(324 + eta*(-633 + 14*eta)))
                      - chi1L*(324 + eta*(-1119 + 2*eta*(172 + eta)) + delta*(324 + eta*(-633 + 14*eta))))/32.) * x3 * sqx //3.5PN
            + (2835/128. - (eta*(-10677852 + 100*eta*(-640863 + eta*(774 + 11*eta)) + 26542080*LAL_GAMMA
                        + 675*(3873 + 3608*eta)*LAL_TWOPI))/622080. - (64*eta*log(16.))/3.) * x4 //4PN
            + (-64*eta/3.) * x4 * logx // 4PN log[x] term
            + ((chi1L2*(1 + delta - 2*eta) + 4*chi1L*chi2L*eta - chi2L2*(-1 + delta + 2*eta))/2.) * x2 // 2PN LOS term, quadratic in spin
            + ((3*(chi1L + chi2L)*eta*(chi1L2*(1 + delta - 2*eta) + 4*chi1L*chi2L*eta - chi2L2*(-1 + delta + 2*eta)))/4.) * x3 // 3.5PN LOS term, cubic in spin
  );
}

/**
    External wrapper function to next-to-next-to-leading (NNLO) in spin-orbit
    expression for the PN Euler angle alpha. This expression is derived by PN
    re-expanding and averaging over the orientation of the spin in the orbital plane.

    - P Schmidt, 2014, http://orca.cf.ac.uk/64062/
    - A Bohé et al, https://dcc.ligo.org/LIGO-T1500602
    - Discussion in arXiv:2004.06503 Sec IV A
*/
/* Wrapper to NNLO PN alpha angle */
REAL8 XLALSimIMRPhenomXPNEuleralphaNNLO(
  const REAL8 f,          /**< Geometric frequency                                                          */
  const REAL8 eta,        /**< Symmetric mass rato                                                          */
  const REAL8 chi1L,      /**< Dimensionless aligned spin of larger BH                                      */
  const REAL8 chi2L,      /**< Dimensionless aligned spin of smaller BH                                     */
  const REAL8 chip,       /**< Effective precession parameter: Schmidt, Ohme, Hannam, PRD, 91,024043 (2015) */
  const REAL8 alpha0      /**< Euler angle at reference Frequency, defines a constant offset                */
)
{
  const REAL8 omega       = LAL_PI * f;
  const REAL8 logomega    = log(omega);
  const REAL8 omega_cbrt  = cbrt(omega);
  const REAL8 omega_cbrt2 = omega_cbrt*omega_cbrt;

  /* Note that we assume: m1 > m2 and dm = m1 - m2 > 0 */
  const REAL8 delta   = sqrt(1.0 - 4.0*eta);
  const REAL8 delta2  = delta*delta;
  const REAL8 delta3  = delta*delta2;

  const REAL8 m1      = 0.5*(1.0 + delta);
  const REAL8 m1_2    = m1*m1;
  const REAL8 m1_3    = m1*m1_2;
  const REAL8 m1_4    = m1*m1_3;
  const REAL8 m1_5    = m1*m1_4;
  const REAL8 m1_6    = m1_3*m1_3;
  const REAL8 m1_8    = m1_4*m1_4;

  const REAL8 m2      = 0.5*(1.0 - delta);
  const REAL8 q       = m1/m2;

  const REAL8 eta2    = eta*eta;
  const REAL8 eta3    = eta*eta2;
  const REAL8 eta4    = eta*eta3;
  const REAL8 eta5    = eta*eta4;
  const REAL8 eta6    = eta*eta5;

  const REAL8 chi_eff = m1*chi1L + m2*chi2L;
  const REAL8 chiL    = (1.0 + q) * chi_eff / q;
  const REAL8 chiL2   = chiL*chiL;
  const REAL8 chip2   = chip*chip;

  const REAL8 alpha1  = (-35/192. - (5*delta)/(64.*m1));
  const REAL8 alpha2  = (((-15*chiL*delta*m1)/128. - (35*chiL*m1_2)/128.)/eta);
  const REAL8 alpha3  = (-5515/3072. + eta*(-515/384. - (15*delta2)/(256.*m1_2) - (175*delta)/(256.*m1)) - (4555*delta)/(7168.*m1)
                        + ((-15*chip2*delta*m1_3)/128. - (35*chip2*m1_4)/128.)/eta2);
  const REAL8 alpha4  = (((5*chiL*delta2)/16. + (5*chiL*delta*m1)/3. + (2545*chiL*m1_2)/1152. + ((2035*chiL*delta*m1)/21504.
                        + (2995*chiL*m1_2)/9216.)/eta + ((-5*chiL*chip2*delta*m1_5)/128.
                        - (35*chiL*chip2*m1_6)/384.)/eta3 - (35*LAL_PI)/48. - (5*delta*LAL_PI)/(16.*m1)));
  const REAL8 alpha5  = ((5*(190512*delta3*eta6 + 2268*delta2*eta3*m1*(eta2*(323 + 784*eta) + 336*(25*chiL2 + chip2)*m1_4)
                      + 7*m1_3*(8024297*eta4 + 857412*eta5 + 3080448*eta6 + 143640*chip2*eta2*m1_4 - 127008*chip2*(-4*chiL2 + chip2)*m1_8
                      + 6048*eta3*((2632*chiL2 + 115*chip2)*m1_4 - 672*chiL*m1_2*LAL_PI)) + 3*delta*m1_2*(5579177*eta4
                      - 80136*eta5 + 3845520*eta6 - 146664*chip2*eta2*m1_4 - 127008*chip2*(-4*chiL2 + chip2)*m1_8
                      + 42336*eta3*((726*chiL2 + 29*chip2)*m1_4 - 96*chiL*m1_2*LAL_PI))))/(6.5028096e7*eta4*m1_3));

  return (alpha1/omega + alpha2/omega_cbrt2 + alpha3/omega_cbrt + alpha4*logomega + alpha5*omega_cbrt + alpha0);
}


/** External wrapper to NNLO PN epsilon angle. See documentation above XLALSimIMRPhenomXPNEuleralphaNNLO. */
REAL8 XLALSimIMRPhenomXPNEulerepsilonNNLO(
  REAL8 f,          /**< Geometric frequency                                                          */
  REAL8 eta,        /**< Symmetric mass rato                                                          */
  REAL8 chi1L,      /**< Dimensionless aligned spin of larger BH                                      */
  REAL8 chi2L,      /**< Dimensionless aligned spin of smaller BH                                     */
  REAL8 chip,       /**< Effective precession parameter: Schmidt, Ohme, Hannam, PRD, 91,024043 (2015) */
  REAL8 epsilon0    /**< Euler angle at reference Frequency, defines a constant offset                */
)
{
  const REAL8 omega       = LAL_PI * f;
  const REAL8 logomega    = log(omega);
  const REAL8 omega_cbrt  = cbrt(omega);
  const REAL8 omega_cbrt2 = omega_cbrt*omega_cbrt;

  /* Note that we assume: m1 > m2 and dm = m1 - m2 > 0 */
  const REAL8 delta   = sqrt(1.0 - 4.0*eta);
  const REAL8 delta2  = delta*delta;
  const REAL8 delta3  = delta*delta2;

  const REAL8 m1      = 0.5*(1.0 + delta);
  const REAL8 m1_2    = m1*m1;
  const REAL8 m1_3    = m1*m1_2;
  const REAL8 m1_4    = m1*m1_3;
  const REAL8 m1_5    = m1*m1_4;
  const REAL8 m1_6    = m1_3*m1_3;
  const REAL8 m1_8    = m1_4*m1_4;

  const REAL8 m2      = 0.5*(1.0 - delta);
  const REAL8 q       = m1/m2;

  const REAL8 eta2    = eta*eta;
  const REAL8 eta3    = eta*eta2;
  const REAL8 eta4    = eta*eta3;
  const REAL8 eta5    = eta*eta4;
  const REAL8 eta6    = eta*eta5;

  const REAL8 chi_eff = m1*chi1L + m2*chi2L;
  const REAL8 chiL    = (1.0 + q) * chi_eff / q;
  const REAL8 chiL2   = chiL*chiL;
  const REAL8 chip2   = chip*chip;

  const REAL8 epsilon1   = (-35/192. - (5*delta)/(64.*m1));
  const REAL8 epsilon2   = (((-15*chiL*delta*m1)/128. - (35*chiL*m1_2)/128.)/eta);
  const REAL8 epsilon3   = (-5515/3072. + eta*(-515/384. - (15*delta2)/(256.*m1_2) - (175*delta)/(256.*m1)) - (4555*delta)/(7168.*m1)
                              + ((-15*chip2*delta*m1_3)/128. - (35*chip2*m1_4)/128.)/eta2);
  const REAL8 epsilon4L  = (((5*chiL*delta2)/16. + (5*chiL*delta*m1)/3. + (2545*chiL*m1_2)/1152. + ((2035*chiL*delta*m1)/21504.
                              + (2995*chiL*m1_2)/9216.)/eta + ((-5*chiL*chip2*delta*m1_5)/128.
                              - (35*chiL*chip2*m1_6)/384.)/eta3 - (35*LAL_PI)/48. - (5*delta*LAL_PI)/(16.*m1)));
  const REAL8 epsilon5   = ((5*(190512*delta3*eta6 + 2268*delta2*eta3*m1*(eta2*(323 + 784*eta) + 336*(25*chiL2 + chip2)*m1_4)
                              + 7*m1_3*(8024297*eta4 + 857412*eta5 + 3080448*eta6 + 143640*chip2*eta2*m1_4 - 127008*chip2*(-4*chiL2 + chip2)*m1_8
                              + 6048*eta3*((2632*chiL2 + 115*chip2)*m1_4 - 672*chiL*m1_2*LAL_PI)) + 3*delta*m1_2*(5579177*eta4
                              - 80136*eta5 + 3845520*eta6 - 146664*chip2*eta2*m1_4 - 127008*chip2*(-4*chiL2 + chip2)*m1_8
                              + 42336*eta3*((726*chiL2 + 29*chip2)*m1_4 - 96*chiL*m1_2*LAL_PI))))/(6.5028096e7*eta4*m1_3));

  return (epsilon1/omega + epsilon2/omega_cbrt2 + epsilon3/omega_cbrt + epsilon4L*logomega + epsilon5*omega_cbrt + epsilon0);

}

/** Internal function to calculate alpha using pre-cached NNLO PN expressions */
double IMRPhenomX_PN_Euler_alpha_NNLO(
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomX Precession Struct */
  const double omega,                       /**< Orbital frequency */
  const double omega_cbrt2,                 /**< Orbital frequency */
  const double omega_cbrt,                  /**< Cubic root of orbital frequency */
  const double logomega                     /**< Natural logarithm of orbital frequency */
)
{
  double alpha = 0.0;

  /* Note that alpha_offset already includes the term alpha0 */
  alpha = (
      pPrec->alpha1  / omega
    + pPrec->alpha2  / omega_cbrt2
    + pPrec->alpha3  / omega_cbrt
    + pPrec->alpha4L * logomega
    + pPrec->alpha5  * omega_cbrt
    - pPrec->alpha_offset
  );

  return alpha;
}

/** Internal function to calculate epsilon using pre-cached NNLO PN expressions */
double IMRPhenomX_PN_Euler_epsilon_NNLO(
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomX Precession Struct */
  const double omega,                       /**< Orbital frequency */
  const double omega_cbrt2,                 /**< Orbital frequency */
  const double omega_cbrt,                  /**< Cubic root of orbital frequency */
  const double logomega                     /**< Natural logarithm of orbital frequency */
)
{
  double epsilon = 0.0;

  /* Note that epsilon_offset already includes the term epsilon0 */
  epsilon =  (
        pPrec->epsilon1  / omega
      + pPrec->epsilon2  / omega_cbrt2
      + pPrec->epsilon3  / omega_cbrt
      + pPrec->epsilon4L * logomega
      + pPrec->epsilon5  * omega_cbrt
      - pPrec->epsilon_offset
  );

  return epsilon;
}

/** Core twisting up routine, see Section III A of arXiv:2004.06503 */
int IMRPhenomXPTwistUp22(
  const REAL8 Mf,                           /**< Frequency (Hz) */
  const COMPLEX16 hAS,                      /**< Underlying aligned-spin IMRPhenomXAS strain */
  IMRPhenomXWaveformStruct *pWF,            /**< IMRPhenomX Waveform Struct    */
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomXP Precession Struct */
  COMPLEX16 *hp,                            /**< [out] h_+ polarization \f$\tilde h_+\f$ */
  COMPLEX16 *hc                             /**< [out] h_x polarization \f$\tilde h_x\f$ */
)
{
  XLAL_CHECK(hp  != NULL, XLAL_EFAULT);
  XLAL_CHECK(hc  != NULL, XLAL_EFAULT);

  /* Euler angles */
  double alpha       = 0.0;
  double epsilon     = 0.0;

  double cBetah      = 0.0;
  double sBetah      = 0.0;

  const double omega       = LAL_PI * Mf;
  const double logomega    = log(omega);
  const double omega_cbrt  = cbrt(omega);
  const double omega_cbrt2 = omega_cbrt * omega_cbrt;

  const double v           = omega_cbrt;

  double s, s2, cos_beta;

  if(pPrec->IMRPhenomXPNRUseTunedAngles)
  {
    alpha       = pPrec->alphaPNR - pPrec->alpha_offset;
    epsilon     = -1.0 * pPrec->gammaPNR - pPrec->epsilon_offset;
    cos_beta    = cos(pPrec->betaPNR);
  }
  else
  {
    switch(pPrec->IMRPhenomXPrecVersion)
    {
      /* ~~~~~ Use NNLO PN Euler Angles - Appendix G of arXiv:2004.06503 and https://dcc.ligo.org/LIGO-T1500602 ~~~~~ */
      case 101:
      case 102:
      case 103:
      case 104:
      {
      alpha         = IMRPhenomX_PN_Euler_alpha_NNLO(pPrec,omega,omega_cbrt2,omega_cbrt,logomega);
      epsilon       = IMRPhenomX_PN_Euler_epsilon_NNLO(pPrec,omega,omega_cbrt2,omega_cbrt,logomega);

      const REAL8 L = XLALSimIMRPhenomXLPNAnsatz(v, pWF->eta/v, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);

      /*
          We ignore the sign of L + SL below:
            s := Sp / (L + SL)
      */
      s        = pPrec->Sperp / (L + pPrec->SL);
      s2       = s*s;
      cos_beta = copysign(1.0, L + pPrec->SL) / sqrt(1.0 + s2);

      break;
      }
      case 220:
      case 221:
      case 222:
      case 223:
      case 224:
      {
      vector vangles = {0.,0.,0.};

      /* ~~~~~ Euler Angles from Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967 ~~~~~ */
      vangles  = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(v,pWF,pPrec);

      alpha    = vangles.x - pPrec->alpha_offset;
      epsilon  = vangles.y - pPrec->epsilon_offset;
      cos_beta = vangles.z;

      break;
      }
      default:
      {
        XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPrecessionVersion not recognized. Recommended default is 223.\n");
        break;
      }
    }
  }


  INT4 status = 0;
  status = IMRPhenomXWignerdCoefficients_cosbeta(&cBetah, &sBetah, cos_beta);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Call to IMRPhenomXWignerdCoefficients_cosbeta failed.");

  /* Useful powers of the Wigner coefficients */
  const REAL8 cBetah2 = cBetah * cBetah;
  const REAL8 cBetah3 = cBetah * cBetah2;
  const REAL8 cBetah4 = cBetah * cBetah3;
  const REAL8 sBetah2 = sBetah * sBetah;
  const REAL8 sBetah3 = sBetah * sBetah2;
  const REAL8 sBetah4 = sBetah * sBetah3;


  /*
      Compute the Wigner d coefficients, see Appendix A of arXiv:2004.06503
        d22  = Table[WignerD[{2, mp, 2}, 0, -\[Beta], 0], {mp, -2, 2}]
        d2m2 = Table[WignerD[{2, mp, -2}, 0, -\[Beta], 0], {mp, -2, 2}]
  */
  // d22  = {d^2_{-2,2} , d^2_{-1,2}, d^2_{0,2}, d^2_{1,2}, d^2_{2,2} }
  const COMPLEX16 d22[5]   = {sBetah4, 2.0*cBetah*sBetah3, pPrec->sqrt6*sBetah2*cBetah2, 2.0*cBetah3*sBetah, cBetah4};

  // d2m2 = {d^2_{-2,-2} , d^2_{-1,-2}, d^2_{0,-2}, d^2_{1,-2}, d^2_{2,-2} }
  const COMPLEX16 d2m2[5]  = {d22[4], -d22[3], d22[2], -d22[1], d22[0]}; /* Exploit symmetry d^2_{-2,m} = (-1)^m d^2_{2,-m}*/

  const COMPLEX16 Y2mA[5]  = {pPrec->Y2m2, pPrec->Y2m1, pPrec->Y20, pPrec->Y21, pPrec->Y22};

  /* Precompute powers of e^{i m alpha} */
  COMPLEX16 cexp_i_alpha         = cexp(+I*alpha);
  COMPLEX16 cexp_2i_alpha        = cexp_i_alpha  * cexp_i_alpha;
  COMPLEX16 cexp_mi_alpha        = 1.0 / cexp_i_alpha;
  COMPLEX16 cexp_m2i_alpha       = cexp_mi_alpha * cexp_mi_alpha;
  COMPLEX16 cexp_im_alpha_l2[5]  = {cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha};

  COMPLEX16 hp_sum               = 0;
  COMPLEX16 hc_sum               = 0;

  REAL8 polarizationSymmetry = pPrec->PolarizationSymmetry;

  /* Loop over m' modes and perform the actual twisting up */
  for(int m=-2; m<=2; m++)
  {
    /* Transfer functions, see Eq. 3.5 and 3.6 of arXiv:2004.06503 */
    COMPLEX16 A2m2emm    = cexp_im_alpha_l2[-m+2] * d2m2[m+2]  * Y2mA[m+2];        /*  = cexp(I*m*alpha) * d22[m+2]   * Y2mA[m+2] */
    COMPLEX16 A22emmstar = cexp_im_alpha_l2[m+2]  * d22[m+2]   * conj(Y2mA[m+2]);  /*  = cexp(-I*m*alpha) * d2m2[m+2] * conj(Y2mA[m+2])  */
    hp_sum +=    A2m2emm + polarizationSymmetry * A22emmstar;
    hc_sum += I*(A2m2emm - polarizationSymmetry * A22emmstar);
  }

  /* Note that \gamma = - \epsilon */
  COMPLEX16 eps_phase_hP = cexp(-2.0*I*epsilon) * hAS / 2.0;

  /* Return h_+ and h_x */
  *hp = eps_phase_hP * hp_sum;
  *hc = eps_phase_hP * hc_sum;

  /* When debugging, save angles to output file. */
  #if DEBUG == 1
      FILE *fileangle;
      char fileSpec[40];

      sprintf(fileSpec, "angles_XP.dat");

      fileangle = fopen(fileSpec,"a");

      // COMPLEX16 cexp_i_epsilon = cexp(I*epsilon);
      // COMPLEX16 cexp_i_betah = cBetah + I*sBetah;
      //fprintf(fileangle, "%.16e  %.16e  %.16e  %.16e  %.16e  %.16e  %.16e\n",  XLALSimIMRPhenomXUtilsMftoHz(Mf, pWF->Mtot), creal(cexp_i_alpha), cimag(cexp_i_alpha), creal(cexp_i_epsilon), cimag(cexp_i_epsilon), creal(cexp_i_betah), cimag(cexp_i_betah));
      fprintf(fileangle, "%.16e  %.16e  %.16e  %.16e\n", XLALSimIMRPhenomXUtilsMftoHz(Mf, pWF->Mtot), alpha, epsilon, cos_beta);
      fclose(fileangle);
  #endif

  return XLAL_SUCCESS;
}

int IMRPhenomXWignerdCoefficients_cosbeta(
  REAL8 *cos_beta_half, /**< [out] cos(beta/2) */
  REAL8 *sin_beta_half, /**< [out] sin(beta/2) */
  const REAL8 cos_beta  /**< cos(beta) */
)
{
  /* Note that the results here are indeed always non-negative */
  *cos_beta_half = + sqrt( fabs(1.0 + cos_beta) / 2.0 );  /* cos(beta/2) */
  *sin_beta_half = + sqrt( fabs(1.0 - cos_beta) / 2.0 );  /* sin(beta/2) */

  return XLAL_SUCCESS;
}

int IMRPhenomXWignerdCoefficients(
  REAL8 *cos_beta_half, /**< [out] cos(beta/2) */
  REAL8 *sin_beta_half, /**< [out] sin(beta/2) */
  const REAL8 v,        /**< Cubic root of (Pi * Frequency (geometric)) */
  IMRPhenomXWaveformStruct *pWF,     /**< IMRPhenomX waveform struct   */
  IMRPhenomXPrecessionStruct *pPrec  /**< IMRPhenomX precession struct */
)
{
  /* Orbital angular momentum */
  const REAL8 L = XLALSimIMRPhenomXLPNAnsatz(v, pWF->eta/v, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);

  /*
      We ignore the sign of L + SL below:
        s := Sp / (L + SL)
  */
  const REAL8 s        = pPrec->Sperp / (L + pPrec->SL);
  const REAL8 s2       = s*s;
  const REAL8 cos_beta = copysign(1.0, L + pPrec->SL) / sqrt(1.0 + s2);

  *cos_beta_half = + sqrt( fabs(1.0 + cos_beta) / 2.0 );  /* cos(beta/2) */
  *sin_beta_half = + sqrt( fabs(1.0 - cos_beta) / 2.0 );  /* sin(beta/2) */

  return XLAL_SUCCESS;
}

/**
    Helper function to check if maximum opening angle > pi/2 or pi/4 and issues a warning. See discussion in https://dcc.ligo.org/LIGO-T1500602
*/
int IMRPhenomXPCheckMaxOpeningAngle(
  IMRPhenomXWaveformStruct *pWF,      /**< IMRPhenomX Waveform Struct */
  IMRPhenomXPrecessionStruct *pPrec   /**< IMRPhenomXP Precession Struct */
)
{
    const REAL8 eta           = pWF->eta;

    /* For now, use the 2PN non-spinning maximum opening angle */
    const REAL8 v_at_max_beta = sqrt(2.0 / 3.0) * sqrt( (-9.0 - eta + sqrt(1539.0 - 1008.0*eta + 19.0*eta*eta)) / (81 - 57*eta + eta*eta) );

    REAL8 cBetah = 0.0;
    REAL8 sBetah = 0.0;

    INT4 status;
    status = IMRPhenomXWignerdCoefficients(&cBetah, &sBetah, v_at_max_beta, pWF, pPrec);
    XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Call to IMRPhenomXWignerdCoefficients failed.");

    const REAL8 L_min    = XLALSimIMRPhenomXL2PNNS(v_at_max_beta,eta);
    const REAL8 max_beta = 2.0 * acos(cBetah);

    /*
      If L + SL becomes < 0, WignerdCoefficients does not track the angle between J and L.
      The model may become pathological as one moves away from the aligned spin limit.

      If this does not happen, then max_beta is the actual maximum opening angle as predicted by the model.
    */
    if ((L_min + pPrec->SL) < 0. && pPrec->chi_p > 0.)
    {
      XLAL_PRINT_WARNING("The maximum opening angle exceeds Pi/2.\nThe model may be pathological in this regime.");
    }
    else if (max_beta > LAL_PI_4)
    {
      XLAL_PRINT_WARNING("The maximum opening angle %g is larger than Pi/4.\nThe model has not been tested against NR in this regime.", max_beta);
    }

    return XLAL_SUCCESS;
}




/* ~~~~~~~~~~ Routines to Calculate PN Angles following the MSA approach of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967 ~~~~~~~~~~*/
/** Returns the 3PN accurate orbital angular momentum as implemented in LALSimInspiralFDPrecAngles_internals.c */
double IMRPhenomX_L_norm_3PN_of_v(const double v, const double v2, const double L_norm, IMRPhenomXPrecessionStruct *pPrec)
{
    return L_norm*(1. + v2*(pPrec->constants_L[0] + v*pPrec->constants_L[1] + v2*(pPrec->constants_L[2] + v*pPrec->constants_L[3] + v2*(pPrec->constants_L[4]))));
}


/** Wrapper to generate \f$\phi_z\f$, \f$\zeta\f$ and \f$\cos \theta_L\f$ at a given frequency */
vector IMRPhenomX_Return_phi_zeta_costhetaL_MSA(
  const double v,                   /**< Velocity                       */
  IMRPhenomXWaveformStruct *pWF,    /**< IMRPhenomX waveform struct     */
  IMRPhenomXPrecessionStruct *pPrec /**< IMRPhenomX precession struct   */
)
{
  vector vout = {0.,0.,0.};

  /* Change code here to determine PN order passed for L */
  const double L_norm    = pWF->eta / v;
  const double J_norm    = IMRPhenomX_JNorm_MSA(L_norm,pPrec);

  double L_norm3PN       = 0.0;

  /* Orbital angular momentum at 3PN, coefficients are pre-cached when initializing precession struct */
  if(pPrec->IMRPhenomXPrecVersion == 222 || pPrec->IMRPhenomXPrecVersion == 223)
  {
    L_norm3PN = IMRPhenomX_L_norm_3PN_of_v(v, v*v, L_norm, pPrec);
  }
  else
  {
    L_norm3PN = XLALSimIMRPhenomXLPNAnsatz(v, L_norm, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);
  }
  const double J_norm3PN = IMRPhenomX_JNorm_MSA(L_norm3PN,pPrec);

  /*
      Get roots to S^2 equation :
          vroots.x = A1 = S_{3}^2
          vroots.y = A2 = S_{-}^2
          vroots.z = A3 = S_{+}^2
  */
  const vector vRoots    = IMRPhenomX_Return_Roots_MSA(L_norm,J_norm,pPrec);

  pPrec->S32  = vRoots.x;
  pPrec->Smi2 = vRoots.y;
  pPrec->Spl2 = vRoots.z;

  pPrec->Spl2mSmi2   = pPrec->Spl2 - pPrec->Smi2;
  pPrec->Spl2pSmi2   = pPrec->Spl2 + pPrec->Smi2;
  pPrec->Spl         = sqrt(pPrec->Spl2);
  pPrec->Smi         = sqrt(pPrec->Smi2);

  const double SNorm = IMRPhenomX_Return_SNorm_MSA(v,pPrec);
  pPrec->S_norm      = SNorm;
  pPrec->S_norm_2    = SNorm * SNorm;

  vector vMSA = {0.,0.,0.};
  if(fabs(pPrec->Smi2 - pPrec->Spl2) > 1.e-5)
  {
    /* Get phiz_0_MSA and zeta_0_MSA */
    vMSA = IMRPhenomX_Return_MSA_Corrections_MSA(v, L_norm, J_norm, pPrec);
  }

  const double phiz_MSA     = vMSA.x;
  const double zeta_MSA     = vMSA.y;

  const double phiz         = IMRPhenomX_Return_phiz_MSA(v,J_norm,pPrec);
  const double zeta         = IMRPhenomX_Return_zeta_MSA(v,pPrec);
  double cos_theta_L        = IMRPhenomX_costhetaLJ(L_norm3PN,J_norm3PN,SNorm);

  vout.x = phiz + phiz_MSA;
  vout.y = zeta + zeta_MSA;
  vout.z = cos_theta_L;

  return vout;
}

/** This function initializes all the core variables required for the MSA system. This will be called first. */
int IMRPhenomX_Initialize_MSA_System(IMRPhenomXWaveformStruct *pWF, IMRPhenomXPrecessionStruct *pPrec, int ExpansionOrder)
{
    /*
        Sanity check on the precession version
    */
    int pflag = pPrec->IMRPhenomXPrecVersion;
    if(pflag != 220 && pflag != 221 && pflag != 222 && pflag != 223 && pflag != 224)
    {
        XLAL_ERROR(XLAL_EINVAL,"Error: MSA system requires IMRPhenomXPrecVersion 220, 221, 222, 223 or 224.\n");
    }

    /*
      First initialize the system of variables needed for Chatziioannou et al, PRD, 88, 063011, (2013), arXiv:1307.4418:
        - Racine et al, PRD, 80, 044010, (2009), arXiv:0812.4413
        - Favata, PRD, 80, 024002, (2009), arXiv:0812.0069
        - Blanchet et al, PRD, 84, 064041, (2011), arXiv:1104.5659
        - Bohe et al, CQG, 30, 135009, (2013), arXiv:1303.7412
    */
    const double eta          = pPrec->eta;
    const double eta2         = pPrec->eta2;
    const double eta3         = pPrec->eta3;
    const double eta4         = pPrec->eta4;

    const double m1 = pWF->m1;
    const double m2 = pWF->m2;

    /* PN Coefficients for d \omega / d t as per LALSimInspiralFDPrecAngles_internals.c */
    const double domegadt_constants_NS[17] = {96./5.,-1486./35.,-264./5.,384.*LAL_PI/5.,34103./945.,13661./105.,944./15.,LAL_PI*(-4159./35.),LAL_PI*(-2268./5.),(16447322263./7276500. + LAL_PI*LAL_PI*512./5. - LAL_LN2*109568./175. -LAL_GAMMA*54784./175.),(-56198689./11340. + LAL_PI*LAL_PI*902./5.),1623./140.,-1121./27.,-54784./525.,-LAL_PI*883./42.,LAL_PI*71735./63.,LAL_PI*73196./63.};
    const double domegadt_constants_SO[18] = {-904./5.,-120.,-62638./105.,4636./5.,-6472./35.,3372./5.,-LAL_PI*720.,-LAL_PI*2416./5.,-208520./63.,796069./105.,-100019./45.,-1195759./945.,514046./105.,-8709./5.,-LAL_PI*307708./105.,LAL_PI*44011./7.,-LAL_PI*7992./7.,LAL_PI*151449./35.};
    const double domegadt_constants_SS[4]  = {-494./5.,-1442./5.,-233./5.,-719./5.};

    const double L_csts_nonspin[9]         = {3./2.,1./6.,27./8.,-19./8.,1./24.,135./16.,-6889/144.+ 41./24.*LAL_PI*LAL_PI,31./24.,7./1296.};
    const double L_csts_spinorbit[6]       = {-14./6.,-3./2.,-11./2.,133./72.,-33./8.,7./4.};


    /*
        Note that Chatziioannou et al use q = m2/m1, where m1 > m2 and therefore q < 1
        IMRPhenomX assumes m1 > m2 and q > 1. For the internal MSA code, flip q and
        dump this to pPrec->qq, where qq explicitly dentoes that this is 0 < q < 1.
    */
    const double q            = m2 / m1; // m2 / m1, q < 1, m1 > m2
    const double invq         = 1.0 / q; // m2 / m1, q < 1, m1 > m2
    pPrec->qq                 = q;
    pPrec->invqq              = invq;

    const double mu           = (m1 * m2) / (m1 + m2);

    #if DEBUG == 1
      printf("m1                = %.6f\n\n",pWF->m1);
      printf("m2                = %.6f\n\n",pWF->m2);
      printf("q (<1)            = %.6f\n\n",pPrec->qq);
    #endif

    /* \delta and powers of \delta in terms of q < 1, should just be m1 - m2 */
    pPrec->delta_qq   = (1.0 - pPrec->qq) / (1.0 + pPrec->qq);
    pPrec->delta2_qq  = pPrec->delta_qq * pPrec->delta_qq;
    pPrec->delta3_qq  = pPrec->delta_qq * pPrec->delta2_qq;
    pPrec->delta4_qq  = pPrec->delta_qq * pPrec->delta3_qq;

    /* Initialize empty vectors */
    vector S1v  = {0.,0.,0.};
    vector S2v  = {0.,0.,0.};

    /* Define source frame such that \hat{L} = {0,0,1} with L_z pointing along \hat{z}. */
    vector Lhat = {0.,0.,1.};

    /* Set LHat variables - these are fixed. */
    pPrec->Lhat_cos_theta = 1.0; /* Cosine of Polar angle of orbital angular momentum */
    pPrec->Lhat_phi       = 0.0; /* Azimuthal angle of orbital angular momentum       */
    pPrec->Lhat_theta     = 0.0; /* Polar angle of orbital angular momentum           */

    /* Dimensionful spin vectors, note eta = m1 * m2 and q = m2/m1  */
    S1v.x = pPrec->chi1x * eta/q; /* eta / q = m1^2 */
    S1v.y = pPrec->chi1y * eta/q;
    S1v.z = pPrec->chi1z * eta/q;

    S2v.x = pPrec->chi2x * eta*q; /* eta * q = m2^2 */
    S2v.y = pPrec->chi2y * eta*q;
    S2v.z = pPrec->chi2z * eta*q;

    REAL8 S1_0_norm = IMRPhenomX_vector_L2_norm(S1v);
    REAL8 S2_0_norm = IMRPhenomX_vector_L2_norm(S2v);

    /* Initial dimensionful spin vectors at reference frequency */
    /* S1 = {S1x,S1y,S1z} */
    pPrec->S1_0.x = S1v.x;
    pPrec->S1_0.y = S1v.y;
    pPrec->S1_0.z = S1v.z;

    /* S2 = {S2x,S2y,S2z} */
    pPrec->S2_0.x = S2v.x;
    pPrec->S2_0.y = S2v.y;
    pPrec->S2_0.z = S2v.z;

    /* Reference velocity v and v^2 */
    pPrec->v_0      = cbrt( pPrec->piGM * pWF->fRef );
    pPrec->v_0_2    = pPrec->v_0 * pPrec->v_0;

    /* Reference orbital angular momenta */
    vector L_0       = {0.,0.,0.};
    L_0              = IMRPhenomX_vector_scalar(Lhat,pPrec->eta / pPrec->v_0);
    pPrec->L_0       = L_0;

    #if DEBUG == 1
      printf("v_0                = %.6f\n\n",pPrec->v_0);

      printf("chi1x              = %.6f\n",pPrec->chi1x);
      printf("chi1y              = %.6f\n",pPrec->chi1y);
      printf("chi1z              = %.6f\n\n",pPrec->chi1z);

      printf("chi2x              = %.6f\n",pPrec->chi2x);
      printf("chi2y              = %.6f\n",pPrec->chi2y);
      printf("chi2z              = %.6f\n\n",pPrec->chi2z);

      printf("S1_0.x             = %.6f\n",pPrec->S1_0.x);
      printf("S1_0.y             = %.6f\n",pPrec->S1_0.y);
      printf("S1_0.z             = %.6f\n",pPrec->S1_0.z);
      printf("S1_0               = %.6f\n\n",S1_0_norm);

      printf("S2_0.x             = %.6f\n",pPrec->S2_0.x);
      printf("S2_0.y             = %.6f\n",pPrec->S2_0.y);
      printf("S2_0.z             = %.6f\n",pPrec->S2_0.z);
      printf("S2_0               = %.6f\n\n",S2_0_norm);
    #endif

    /* Inner products used in MSA system */
    double dotS1L, dotS2L, dotS1Ln, dotS2Ln, dotS1S2;

    dotS1L  = IMRPhenomX_vector_dot_product(S1v,Lhat);
    dotS2L  = IMRPhenomX_vector_dot_product(S2v,Lhat);
    dotS1S2 = IMRPhenomX_vector_dot_product(S1v,S2v);
    dotS1Ln = dotS1L / S1_0_norm;
    dotS2Ln = dotS2L / S2_0_norm;

    /* Add dot products to struct */
    pPrec->dotS1L  = dotS1L;
    pPrec->dotS2L  = dotS2L;
    pPrec->dotS1S2 = dotS1S2;
    pPrec->dotS1Ln = dotS1Ln;
    pPrec->dotS2Ln = dotS2Ln;

    #if DEBUG == 1
      printf("Lhat_0.x               = %.6f\n",Lhat.x);
      printf("Lhat_0.y               = %.6f\n",Lhat.y);
      printf("Lhat_0.z               = %.6f\n\n",Lhat.z);

      printf("dotS1L                 = %.6f\n",pPrec->dotS1L);
      printf("dotS2L                 = %.6f\n",pPrec->dotS2L);
      printf("dotS1Ln                = %.6f\n",pPrec->dotS1Ln);
      printf("dotS2Ln                = %.6f\n",pPrec->dotS2Ln);
      printf("dotS1S2                = %.6f\n\n",pPrec->dotS1S2);
    #endif

    /* Coeffcients for PN orbital angular momentum at 3PN, as per LALSimInspiralFDPrecAngles_internals.c */
    pPrec->constants_L[0] = (L_csts_nonspin[0] + eta*L_csts_nonspin[1]);
    pPrec->constants_L[1] = IMRPhenomX_Get_PN_beta(L_csts_spinorbit[0], L_csts_spinorbit[1], pPrec);
    pPrec->constants_L[2] = (L_csts_nonspin[2] + eta*L_csts_nonspin[3] + eta*eta*L_csts_nonspin[4]);
    pPrec->constants_L[3] = IMRPhenomX_Get_PN_beta((L_csts_spinorbit[2]+L_csts_spinorbit[3]*eta), (L_csts_spinorbit[4]+L_csts_spinorbit[5]*eta), pPrec);
    pPrec->constants_L[4] = (L_csts_nonspin[5]+L_csts_nonspin[6]*eta +L_csts_nonspin[7]*eta*eta+L_csts_nonspin[8]*eta*eta*eta);

    /* Effective total spin */
    const double Seff  = (1.0 + q) * pPrec->dotS1L + (1 + (1.0/q))*pPrec->dotS2L;
    const double Seff2 =  Seff * Seff;

    pPrec->Seff    = Seff;
    pPrec->Seff2   = Seff2;

    #if DEBUG == 1
      printf("Seff             = %.6f\n\n",pPrec->Seff);
    #endif

    /* Initial total spin, S = S1 + S2 */
    vector S0 = {0.,0.,0.};
    S0 = IMRPhenomX_vector_sum(S1v,S2v);

    /* Cache total spin in the precession struct */
    pPrec->S_0 = S0;

    #if DEBUG == 1
      printf("S_0_x             = %.6f\n",pPrec->S_0.x);
      printf("S_0_y             = %.6f\n",pPrec->S_0.y);
      printf("S_0_z             = %.6f\n\n",pPrec->S_0.z);
    #endif

    /* Initial total angular momentum, J = L + S1 + S2 */
    pPrec->J_0 = IMRPhenomX_vector_sum(pPrec->L_0,pPrec->S_0);

    #if DEBUG == 1
      printf("J_0_x             = %.6f\n",pPrec->J_0.x);
      printf("J_0_y             = %.6f\n",pPrec->J_0.y);
      printf("J_0_z             = %.6f\n\n",pPrec->J_0.z);
    #endif

    /* Norm of total initial spin */
    pPrec->S_0_norm   = IMRPhenomX_vector_L2_norm(S0);
    pPrec->S_0_norm_2 = pPrec->S_0_norm * pPrec->S_0_norm;

    /* Norm of orbital and total angular momenta */
    pPrec->L_0_norm   = IMRPhenomX_vector_L2_norm(pPrec->L_0);
    pPrec->J_0_norm   = IMRPhenomX_vector_L2_norm(pPrec->J_0);

    const double L0norm = pPrec->L_0_norm;
    const double J0norm = pPrec->J_0_norm;

    #if DEBUG == 1
      printf("L_0_norm             = %.6f\n",pPrec->L_0_norm);
      printf("J_0_norm             = %.6f\n\n",pPrec->J_0_norm);
    #endif

    /* Useful powers */
    pPrec->S_0_norm_2 = pPrec->S_0_norm * pPrec->S_0_norm;
    pPrec->J_0_norm_2 = pPrec->J_0_norm * pPrec->J_0_norm;
    pPrec->L_0_norm_2 = pPrec->L_0_norm * pPrec->L_0_norm;

    /* Vector for obtaining B, C, D coefficients */
    UNUSED vector vBCD;
    vBCD = IMRPhenomX_Return_Spin_Evolution_Coefficients_MSA(pPrec->L_0_norm,pPrec->J_0_norm, pPrec);

    #if DEBUG == 1
      printf("B             = %.6f\n",vBCD.x);
      printf("C             = %.6f\n",vBCD.y);
      printf("D             = %.6f\n\n",vBCD.z);
    #endif

    /*
        Get roots to S^2 equation : S^2_+, S^2_-, S^2_3
            vroots.x = A1 = S_{3}^2
            vroots.y = A2 = S_{-}^2
            vroots.z = A3 = S_{+}^2
    */
    vector vRoots = {0.,0.,0.};

    vRoots = IMRPhenomX_Return_Roots_MSA(pPrec->L_0_norm,pPrec->J_0_norm,pPrec);

    // Set roots
    pPrec->Spl2 = vRoots.z;
    pPrec->Smi2 = vRoots.y;
    pPrec->S32  = vRoots.x;

    // S^2_+ + S^2_-
    pPrec->Spl2pSmi2 = pPrec->Spl2 + pPrec->Smi2;

    // S^2_+ - S^2_-
    pPrec->Spl2mSmi2 = pPrec->Spl2 - pPrec->Smi2;

    // S_+ and S_-
    pPrec->Spl       = sqrt(pPrec->Spl2);
    pPrec->Smi       = sqrt(pPrec->Smi2);

    /* Eq. 45 of PRD 95, 104004, (2017), arXiv:1703.03967, set from initial conditions */
    pPrec->SAv2      = 0.5 * ( pPrec->Spl2pSmi2 );
    pPrec->SAv       = sqrt(pPrec->SAv2);
    pPrec->invSAv2   = 1.0 / pPrec->SAv2;
    pPrec->invSAv    = 1.0 / pPrec->SAv;

    #if DEBUG == 1
      printf("From vRoots... \n");
      printf("Spl2             = %.6f\n",pPrec->Spl2);
      printf("Smi2             = %.6f\n",pPrec->Smi2);
      printf("S32              = %.6f\n",pPrec->S32);
      printf("SAv2             = %.6f\n",pPrec->SAv2);
      printf("SAv              = %.6f\n\n",pPrec->SAv);
    #endif

    /* c_1 is determined by Eq. 41 of PRD, 95, 104004, (2017), arXiv:1703.03967 */
    const double c_1          = 0.5 * (J0norm*J0norm - L0norm*L0norm - pPrec->SAv2) / pPrec->L_0_norm * eta;
    const double c1_2         = c_1 * c_1;

    /* Useful powers and combinations of c_1 */
    pPrec->c1          = c_1;
    pPrec->c12         = c_1 * c_1;
    pPrec->c1_over_eta = c_1 / eta;

    /* Average spin couplings over one precession cycle: A9 - A14 of arXiv:1703.03967  */
    const double omqsq = (1.0 - q) * (1.0 - q) + 1e-16;
    const double omq2  = (1.0 - q * q) + 1e-16;

    /* Precession averaged spin couplings, Eq. A9 - A14 of arXiv:1703.03967, note that we only use the initial values  */
    pPrec->S1L_pav     = (c_1 * (1.0 + q) - q * eta * Seff) / (eta * omq2);
    pPrec->S2L_pav     = - q * (c_1 * (1.0 + q) - eta * Seff) / (eta * omq2);
    pPrec->S1S2_pav    = 0.5 * pPrec->SAv2 - 0.5*(pPrec->S1_norm_2 + pPrec->S2_norm_2);
    pPrec->S1Lsq_pav   = (pPrec->S1L_pav*pPrec->S1L_pav) + ((pPrec->Spl2mSmi2)*(pPrec->Spl2mSmi2) * pPrec->v_0_2) / (32.0 * eta2 * omqsq);
    pPrec->S2Lsq_pav   = (pPrec->S2L_pav*pPrec->S2L_pav) + (q*q*(pPrec->Spl2mSmi2)*(pPrec->Spl2mSmi2) * pPrec->v_0_2) / (32.0 * eta2 * omqsq);
    pPrec->S1LS2L_pav  = pPrec->S1L_pav*pPrec->S2L_pav - q * (pPrec->Spl2mSmi2)*(pPrec->Spl2mSmi2)*pPrec->v_0_2 / (32.0 * eta2 * omqsq);

    /* Spin couplings in arXiv:1703.03967 */
    pPrec->beta3       = ( (113./12.) + (25./4.)*(m2/m1) )*pPrec->S1L_pav + ( (113./12.) + (25./4.)*(m1/m2) )*pPrec->S2L_pav;

    pPrec->beta5       = ( ( (31319./1008.) - (1159./24.)*eta) + (m2/m1)*((809./84) - (281./8.)*eta) )*pPrec->S1L_pav
                          + ( ( (31319./1008.) - (1159./24.)*eta) + (m1/m2)*((809./84) - (281./8.)*eta) )*pPrec->S2L_pav;

    pPrec->beta6       = LAL_PI * ( ( (75./2.) + (151./6.)*(m2/m1))*pPrec->S1L_pav + ( (75./2.) + (151./6.)*(m1/m2))*pPrec->S2L_pav );

    pPrec->beta7       =   (
                              ( (130325./756) - (796069./2016)*eta + (100019./864.)*eta2 ) + (m2/m1)*( (1195759./18144) - (257023./1008.)*eta + (2903/32.)*eta2 ) * pPrec->S1L_pav
                           +  ( (130325./756) - (796069./2016)*eta + (100019./864.)*eta2 ) + (m1/m2)*( (1195759./18144) - (257023./1008.)*eta + (2903/32.)*eta2 ) * pPrec->S2L_pav
                          );

    pPrec->sigma4      =    (1.0/mu) * ( (247./48.)*pPrec->S1S2_pav - (721./48.)*pPrec->S1L_pav*pPrec->S2L_pav )
                        + (1.0/(m1*m1)) * ( (233./96.)*pPrec->S1_norm_2 - (719./96.)*pPrec->S1Lsq_pav )
                        + (1.0/(m2*m2)) * ( (233./96.)*pPrec->S2_norm_2 - (719./96.)*pPrec->S2Lsq_pav );

    /* Compute PN coefficients using precession-averaged spin couplings */
    pPrec->a0          = 96.0 * eta / 5.0;

    /* These are all normalized by a factor of a0 */
    pPrec->a2          = -(743./336.) - (11.0/4.)*eta;
    pPrec->a3          = 4.0 * LAL_PI - pPrec->beta3;
    pPrec->a4          = (34103./18144.) + (13661./2016.)*eta + (59./18.)*eta2 - pPrec->sigma4;
    pPrec->a5          = -(4159./672.)*LAL_PI - (189./8.)*LAL_PI*eta - pPrec->beta5;
    pPrec->a6          = (16447322263./139708800.) + (16./3.)*LAL_PI*LAL_PI - (856./105)*log(16.) - (1712./105.)*LAL_GAMMA - pPrec->beta6
                          + eta*( (451./48)*LAL_PI*LAL_PI - (56198689./217728.) ) + eta2*(541./896.) - eta3*(5605./2592.);
    pPrec->a7          = -(4415./4032.)*LAL_PI + (358675./6048.)*LAL_PI*eta + (91495./1512.)*LAL_PI*eta2 - pPrec->beta7;

    // Coefficients are weighted by an additional factor of a_0
    pPrec->a2         *= pPrec->a0;
    pPrec->a3         *= pPrec->a0;
    pPrec->a4         *= pPrec->a0;
    pPrec->a5         *= pPrec->a0;
    pPrec->a6         *= pPrec->a0;
    pPrec->a7         *= pPrec->a0;

    #if DEBUG == 1
        printf("a0     = %.6f\n",pPrec->a0);
        printf("a2     = %.6f\n",pPrec->a2);
        printf("a3     = %.6f\n",pPrec->a3);
        printf("a4     = %.6f\n",pPrec->a4);
        printf("a5     = %.6f\n\n",pPrec->a5);
    #endif

    /* For versions 222 and 223, we compute PN coefficients using initial spin couplings, as per LALSimInspiralFDPrecAngles_internals.c  */
    if(pflag == 222 || pflag == 223)
    {
      pPrec->a0 = eta*domegadt_constants_NS[0];
      pPrec->a2 = eta*(domegadt_constants_NS[1] + eta*(domegadt_constants_NS[2]));
      pPrec->a3 = eta*(domegadt_constants_NS[3] + IMRPhenomX_Get_PN_beta(domegadt_constants_SO[0], domegadt_constants_SO[1], pPrec));
      pPrec->a4 = eta*(domegadt_constants_NS[4] + eta*(domegadt_constants_NS[5] + eta*(domegadt_constants_NS[6])) + IMRPhenomX_Get_PN_sigma(domegadt_constants_SS[0], domegadt_constants_SS[1], pPrec) + IMRPhenomX_Get_PN_tau(domegadt_constants_SS[2], domegadt_constants_SS[3], pPrec));
      pPrec->a5 = eta*(domegadt_constants_NS[7] + eta*(domegadt_constants_NS[8]) + IMRPhenomX_Get_PN_beta((domegadt_constants_SO[2] + eta*(domegadt_constants_SO[3])), (domegadt_constants_SO[4] + eta*(domegadt_constants_SO[5])), pPrec));
    }

    /* Debugging */
    #if DEBUG == 1
        printf("Using list of coefficients... \n");
        printf("a0     = %.6f\n",pPrec->a0);
        printf("a2     = %.6f\n",pPrec->a2);
        printf("a3     = %.6f\n",pPrec->a3);
        printf("a4     = %.6f\n",pPrec->a4);
        printf("a5     = %.6f\n\n",pPrec->a5);
    #endif

    /* Useful powers of a_0 */
    pPrec->a0_2 = pPrec->a0*pPrec->a0;
    pPrec->a0_3 = pPrec->a0_2*pPrec->a0;
    pPrec->a2_2 = pPrec->a2*pPrec->a2;

    /*
      Calculate g coefficients as in Appendix A of Chatziioannou et al, PRD, 95, 104004, (2017), arXiv:1703.03967.
      These constants are used in TaylorT2 where domega/dt is expressed as an inverse polynomial
    */
    pPrec->g0   = 1. / pPrec->a0;

    // Eq. A2 (1703.03967)
    pPrec->g2   = -(pPrec->a2 / pPrec->a0_2);

    // Eq. A3 (1703.03967)
    pPrec->g3   = -(pPrec->a3/pPrec->a0_2);

    // Eq.A4 (1703.03967)
    pPrec->g4   = -(pPrec->a4*pPrec->a0 - pPrec->a2_2) / pPrec->a0_3;

    // Eq. A5 (1703.03967)
    pPrec->g5   = -(pPrec->a5*pPrec->a0 - 2.0*pPrec->a3*pPrec->a2) / pPrec->a0_3;

    #if DEBUG == 1
      printf("g0     = %.6f\n",pPrec->g0);
      printf("g2     = %.6f\n",pPrec->g2);
      printf("g3     = %.6f\n",pPrec->g3);
      printf("g4     = %.6f\n",pPrec->g4);
      printf("g5     = %.6f\n\n",pPrec->g5);
    #endif

    // Useful powers of delta
    const double delta  = pPrec->delta_qq;
    const double delta2 = delta * delta;
    const double delta3 = delta * delta2;
    const double delta4 = delta * delta3;

    // These are the phase coefficients of Eq. 51 of PRD, 95, 104004, (2017), arXiv:1703.03967
    pPrec->psi0   = 0.0;
    pPrec->psi1   = 0.0;
    pPrec->psi2   = 0.0;

    /* \psi_1 is defined in Eq. C1 of Appendix C in PRD, 95, 104004, (2017), arXiv:1703.03967  */
    pPrec->psi1   = 3.0 * (2.0 * eta2 * Seff - c_1) / (eta * delta2);

    double c_1_over_nu   = pPrec->c1_over_eta;
    double c_1_over_nu_2 = c_1_over_nu * c_1_over_nu;
    double one_p_q_sq    = (1.+q) * (1.+q);
    double Seff_2        = Seff * Seff;
    double q_2           = q * q;
    double one_m_q_sq    = (1.-q)*(1.-q);
    double one_m_q2_2    = (1. - q_2) * (1. - q_2);
    double one_m_q_4     = one_m_q_sq * one_m_q_sq;

    double term1, term2, term3, term4, term5, term6, term7, term8;

    /*  This implements the Delta term as in LALSimInspiralFDPrecAngles.c
        c.f. https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiralFDPrecAngles_internals.c#L145
    */
    if(pflag == 222 || pflag == 223)
    {
      const double Del1 = 4. * c_1_over_nu_2 * one_p_q_sq;
      const double Del2 = 8. * c_1_over_nu * q * (1. + q) * Seff;
      const double Del3 = 4. * (one_m_q2_2 * pPrec->S1_norm_2 - q_2 * Seff_2);
      const double Del4 = 4. * c_1_over_nu_2 * q_2 * one_p_q_sq;
      const double Del5 = 8. * c_1_over_nu * q_2 * (1. + q) * Seff;
      const double Del6 = 4. * (one_m_q2_2 * pPrec->S2_norm_2 - q_2 * Seff_2);
      pPrec->Delta      = sqrt( fabs( (Del1 - Del2 - Del3) * (Del4 - Del5 - Del6) ));
    }
    else
    {
      /* Coefficients of \Delta as defined in Eq. C3 of Appendix C in PRD, 95, 104004, (2017), arXiv:1703.03967. */
      term1  = c1_2 * eta / (q * delta4);
      term2  = -2.0 * c_1 * eta3 * (1.0 + q) * Seff / (q * delta4);
      term3  = -eta2 * (delta2 * pPrec->S1_norm_2 - eta2 * Seff2) / delta4;
      /*
          Is this 1) (c1_2 * q * eta / delta4) or 2) c1_2*eta2/delta4?

            - In paper.pdf, the expression 1) is used.

          Using eta^2 leads to higher frequency oscillations, use q * eta
      */
      term4  = c1_2 * eta * q / delta4;
      term5  = -2.0*c_1*eta3*(1.0 + q)*Seff / delta4;
      term6  = -eta2 * (delta2*pPrec->S2_norm_2 - eta2*Seff2) / delta4;

      /* \Delta as in Eq. C3 of Appendix C in PRD, 95, 104004, (2017) */
      pPrec->Delta  = sqrt( fabs( (term1 + term2 + term3) * (term4 + term5 + term6) ) );
    }

    /*  This implements the Delta term as in LALSimInspiralFDPrecAngles.c
        c.f. https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiralFDPrecAngles_internals.c#L160
    */
    if(pflag == 222 || pflag == 223)
    {
      const double u1 = 3. * pPrec->g2 / pPrec->g0;
      const double u2 = 0.75 * one_p_q_sq / one_m_q_4;
      const double u3 = -20. * c_1_over_nu_2 * q_2 * one_p_q_sq;
      const double u4 = 2. * one_m_q2_2 * (q * (2. + q) * pPrec->S1_norm_2 + (1. + 2. * q) * pPrec->S2_norm_2 - 2. * q * pPrec->SAv2);
      const double u5 = 2. * q_2 * (7. + 6. * q + 7. * q_2) * 2. * c_1_over_nu * Seff;
      const double u6 = 2. * q_2 * (3. + 4. * q + 3. * q_2) * Seff_2;
      const double u7 = q * pPrec->Delta;

      /* Eq. C2 (1703.03967) */
      pPrec->psi2 = u1 + u2*(u3 + u4 + u5 - u6 + u7);
    }
    else
    {
      /* \psi_2 is defined in Eq. C2 of Appendix C in PRD, 95, 104004, (2017). Here we implement system of equations as in paper.pdf */
      term1         = 3.0 * pPrec->g2 / pPrec->g0;

      /* q^2 or no q^2 in term2? Consensus on retaining q^2 term: https://git.ligo.org/waveforms/reviews/phenompv3hm/issues/7 */
      term2         = 3.0 * q * q / (2.0 * eta3);
      term3         = 2.0 * pPrec->Delta;
      term4         = -2.0*eta2*pPrec->SAv2 / delta2;
      term5         = -10.*eta*c1_2 / delta4;
      term6         = 2.0 * eta2 * (7.0 + 6.0*q + 7.0*q*q) * c_1 * Seff / (omqsq * delta2);
      term7         = -eta3 * (3.0 + 4.0*q + 3.0*q*q) * Seff2 / (omqsq * delta2);
      term8         = eta * (q * (2.0+q)*pPrec->S1_norm_2 + (1.0 + 2.0*q)*pPrec->S2_norm_2) / ( omqsq );

      /* \psi_2, C2 of Appendix C of PRD, 95, 104004, (2017)  */
      pPrec->psi2  = term1 + term2 * (term3 + term4 + term5 + term6 + term7 + term8);
    }

    #if DEBUG == 1
      printf("psi1     = %.6f\n",pPrec->psi1);
      printf("psi2     = %.6f\n\n",pPrec->psi2);
    #endif

    /* Eq. D1 of PRD, 95, 104004, (2017), arXiv:1703.03967  */
    const double Rm       = pPrec->Spl2 - pPrec->Smi2;
    const double Rm_2     = Rm * Rm;

    /* Eq. D2 and D3 Appendix D of PRD, 95, 104004, (2017), arXiv:1703.03967   */
    const double cp       = pPrec->Spl2 * eta2 - c1_2;
    double cm             = pPrec->Smi2 * eta2 - c1_2;

    /*
      Check if cm goes negative, this is likely pathological. If so, set MSA_ERROR to 1, so that waveform generator can handle
      the error approriately
    */
    // if(cm < 0.0)
    // {
    //   pPrec->MSA_ERROR = 1;
    //   XLAL_PRINT_ERROR("Error, coefficient cm = %.16f, which is negative and likely to be pathological. Triggering MSA failure.\n",cm);
    // }

    /* fabs is here to help enforce positive definite cpcm */
    const double cpcm      = fabs( cp * cm );
    const double sqrt_cpcm = sqrt(cpcm);

    /* Eq. D4 in PRD, 95, 104004, (2017), arXiv:1703.03967 ; Note difference to published version.  */
    const double a1dD  = 0.5 + 0.75/eta;

    /* Eq. D5 in PRD, 95, 104004, (2017), arXiv:1703.03967  */
    const double a2dD  = -0.75*Seff/eta;

    /* Eq. E3 in PRD, 95, 104004, (2017), arXiv:1703.03967 ; Note that this is Rm * D2   */
    const double D2RmSq = (cp - sqrt_cpcm) / eta2 ;

    /* Eq. E4 in PRD, 95, 104004, (2017), arXiv:1703.03967 ; Note that this is Rm^2 * D4  */
    const double D4RmSq = -0.5*Rm*sqrt_cpcm/eta2 - cp/eta4*(sqrt_cpcm - cp);

    const double S0m  = pPrec->S1_norm_2 - pPrec->S2_norm_2;

    /* Difference of spin norms squared, as used in Eq. D6 of PRD, 95, 104004, (2017), arXiv:1703.03967  */
    double aw   = (-3.*(1. + q)/q*(2.*(1. + q)*eta2*Seff*c_1 - (1. + q)*c1_2 + (1. - q)*eta2*S0m));
    double cw   = 3./32./eta*Rm_2;
    double dw   = 4.0*cp - 4.0*D2RmSq*eta2;
    double hw   = -2.0*(2.0*D2RmSq - Rm)*c_1;
    double fw   = Rm*D2RmSq - D4RmSq - 0.25*Rm_2;

    const double adD  = aw / dw;
    const double hdD  = hw / dw;
    const double cdD  = cw / dw;
    const double fdD  = fw / dw;

    const double gw   = 3./16./eta2/eta*Rm_2*(c_1 - eta2*Seff);
    const double gdD  = gw / dw;

    /* Useful powers of the coefficients */
    const double hdD_2        = hdD   * hdD;
    const double adDfdD       = adD * fdD;
    const double adDfdDhdD    = adDfdD * hdD;
    const double adDhdD_2     = adD * hdD_2;

    #if DEBUG == 1
        printf("\na1dD      = %.6f\n",a1dD);
        printf("a2dD      = %.6f\n",a2dD);
        printf("adD       = %.6f\n",adD);
        printf("cdD       = %.6f\n",cdD);
        printf("hdD       = %.6f\n",hdD);
        printf("fdD       = %.6f\n",fdD);
        printf("Rm        = %.6f\n",Rm);
        printf("Delta     = %.6f\n",pPrec->Delta);
        printf("sqrt_cpcm = %.6f\n",sqrt_cpcm);
        printf("c1        = %.6f\n",pPrec->c1);
        printf("gdD       = %.6f\n\n",gdD);
    #endif

    // Eq. D10 in PRD, 95, 104004, (2017), arXiv:1703.03967
    pPrec->Omegaz0 = a1dD + adD;

    // Eq. D11 in PRD, 95, 104004, (2017), arXiv:1703.03967
    pPrec->Omegaz1 = a2dD - adD*Seff - adD*hdD;

    // Eq. D12 in PRD, 95, 104004, (2017), arXiv:1703.03967
    pPrec->Omegaz2 = adD*hdD*Seff + cdD - adD*fdD + adD*hdD_2;

    // Eq. D13 in PRD, 95, 104004, (2017), arXiv:1703.03967
    pPrec->Omegaz3 = (adDfdD - cdD - adDhdD_2)*(Seff + hdD) + adDfdDhdD;

    // Eq. D14 in PRD, 95, 104004, (2017), arXiv:1703.03967
    pPrec->Omegaz4 = (cdD + adDhdD_2 - 2.0*adDfdD)*(hdD*Seff + hdD_2 - fdD) - adD*fdD*fdD;

    // Eq. D15 in PRD, 95, 104004, (2017), arXiv:1703.03967
    pPrec->Omegaz5 = (cdD - adDfdD + adDhdD_2) * fdD * (Seff + 2.0*hdD) - (cdD + adDhdD_2 - 2.0*adDfdD) * hdD_2 * (Seff + hdD) - adDfdD*fdD*hdD;

    #if DEBUG == 1
        printf("Omegaz0     = %.6f\n",pPrec->Omegaz0);
        printf("Omegaz1     = %.6f\n",pPrec->Omegaz1);
        printf("Omegaz2     = %.6f\n",pPrec->Omegaz2);
        printf("Omegaz3     = %.6f\n",pPrec->Omegaz3);
        printf("Omegaz4     = %.6f\n",pPrec->Omegaz4);
        printf("Omegaz5     = %.6f\n\n",pPrec->Omegaz5);
    #endif

    /*
        If Omegaz5 > 1000, this is larger than we expect and the system may be pathological.
        - Set MSA_ERROR = 1 to trigger an error
    */
    if(fabs(pPrec->Omegaz5) > 1000.0)
    {
      pPrec->MSA_ERROR = 1;
      XLAL_PRINT_WARNING("Warning, |Omegaz5| = %.16f, which is larger than expected and may be pathological. Triggering MSA failure.\n",pPrec->Omegaz5);
    }

    const double g0 = pPrec->g0;

    /* Coefficients of Eq. 65, as defined in Equations D16 - D21 of PRD, 95, 104004, (2017), arXiv:1703.03967 */
    pPrec->Omegaz0_coeff = 3.0 * g0 * pPrec->Omegaz0;
    pPrec->Omegaz1_coeff = 3.0 * g0 * pPrec->Omegaz1;
    pPrec->Omegaz2_coeff = 3.0 * (g0 * pPrec->Omegaz2 + pPrec->g2*pPrec->Omegaz0);
    pPrec->Omegaz3_coeff = 3.0 * (g0 * pPrec->Omegaz3 + pPrec->g2*pPrec->Omegaz1 + pPrec->g3*pPrec->Omegaz0);
    pPrec->Omegaz4_coeff = 3.0 * (g0 * pPrec->Omegaz4 + pPrec->g2*pPrec->Omegaz2 + pPrec->g3*pPrec->Omegaz1 + pPrec->g4*pPrec->Omegaz0);
    pPrec->Omegaz5_coeff = 3.0 * (g0 * pPrec->Omegaz5 + pPrec->g2*pPrec->Omegaz3 + pPrec->g3*pPrec->Omegaz2 + pPrec->g4*pPrec->Omegaz1 + pPrec->g5*pPrec->Omegaz0);

    /* Coefficients of zeta: in Appendix E of PRD, 95, 104004, (2017), arXiv:1703.03967  */
    const double c1oveta2 = c_1 / eta2;
    pPrec->Omegazeta0 = pPrec->Omegaz0;
    pPrec->Omegazeta1 = pPrec->Omegaz1 + pPrec->Omegaz0 * c1oveta2;
    pPrec->Omegazeta2 = pPrec->Omegaz2 + pPrec->Omegaz1 * c1oveta2;
    pPrec->Omegazeta3 = pPrec->Omegaz3 + pPrec->Omegaz2 * c1oveta2 + gdD;
    pPrec->Omegazeta4 = pPrec->Omegaz4 + pPrec->Omegaz3 * c1oveta2 - gdD*Seff - gdD*hdD;
    pPrec->Omegazeta5 = pPrec->Omegaz5 + pPrec->Omegaz4 * c1oveta2 + gdD*hdD*Seff + gdD*(hdD_2 - fdD);

    #if DEBUG == 1
        printf("Omegazeta0     = %.6f\n",pPrec->Omegazeta0);
        printf("Omegazeta1     = %.6f\n",pPrec->Omegazeta1);
        printf("Omegazeta2     = %.6f\n",pPrec->Omegazeta2);
        printf("Omegazeta3     = %.6f\n",pPrec->Omegazeta3);
        printf("Omegazeta4     = %.6f\n",pPrec->Omegazeta4);
        printf("Omegazeta5     = %.6f\n\n",pPrec->Omegazeta5);
    #endif

    pPrec->Omegazeta0_coeff = -pPrec->g0 * pPrec->Omegazeta0;
    pPrec->Omegazeta1_coeff = -1.5 * pPrec->g0 * pPrec->Omegazeta1;
    pPrec->Omegazeta2_coeff = -3.0*(pPrec->g0 * pPrec->Omegazeta2 + pPrec->g2*pPrec->Omegazeta0);
    pPrec->Omegazeta3_coeff = 3.0*(pPrec->g0 * pPrec->Omegazeta3 + pPrec->g2*pPrec->Omegazeta1 + pPrec->g3*pPrec->Omegazeta0);
    pPrec->Omegazeta4_coeff = 3.0*(pPrec->g0 * pPrec->Omegazeta4 + pPrec->g2*pPrec->Omegazeta2 + pPrec->g3*pPrec->Omegazeta1 + pPrec->g4*pPrec->Omegazeta0);
    pPrec->Omegazeta5_coeff = 1.5*(pPrec->g0*pPrec->Omegazeta5 + pPrec->g2*pPrec->Omegazeta3 + pPrec->g3*pPrec->Omegazeta2 + pPrec->g4*pPrec->Omegazeta1 + pPrec->g5*pPrec->Omegazeta0);

    /* Expansion order of corrections to retain */
    switch(ExpansionOrder)
    {
      /* Generate all orders */
      case -1:
      {
        break;
      }
      case 1:
      {
        pPrec->Omegaz1_coeff    = 0.0;
        pPrec->Omegazeta1_coeff = 0.0;
        #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                    __attribute__ ((fallthrough));
        #endif

      }
      case 2:
      {
        pPrec->Omegaz2_coeff    = 0.0;
        pPrec->Omegazeta2_coeff = 0.0;
        #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                    __attribute__ ((fallthrough));
        #endif

      }
      case 3:
      {
        pPrec->Omegaz3_coeff    = 0.0;
        pPrec->Omegazeta3_coeff = 0.0;
        #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                    __attribute__ ((fallthrough));
        #endif

      }
      case 4:
      {
        pPrec->Omegaz4_coeff    = 0.0;
        pPrec->Omegazeta4_coeff = 0.0;
        #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                    __attribute__ ((fallthrough));
        #endif

      }
      case 5:
      {
        pPrec->Omegaz5_coeff    = 0.0;
        pPrec->Omegazeta5_coeff = 0.0;
        break;
      }
      default:
      {
        XLAL_ERROR(XLAL_EDOM, "Expansion order for MSA corrections = %i not recognized. Default is 5. Allowed values are: [-1,1,2,3,4,5].",ExpansionOrder);
      }
    }

    #if DEBUG == 1
      printf("Omegaz0_coeff     = %.6f\n",pPrec->Omegaz0_coeff);
      printf("Omegaz1_coeff     = %.6f\n",pPrec->Omegaz1_coeff);
      printf("Omegaz2_coeff     = %.6f\n",pPrec->Omegaz2_coeff);
      printf("Omegaz3_coeff     = %.6f\n",pPrec->Omegaz3_coeff);
      printf("Omegaz4_coeff     = %.6f\n",pPrec->Omegaz4_coeff);
      printf("Omegaz5_coeff     = %.6f\n\n",pPrec->Omegaz5_coeff);

      printf("Omegazeta0_coeff     = %.6f\n",pPrec->Omegazeta0_coeff);
      printf("Omegazeta1_coeff     = %.6f\n",pPrec->Omegazeta1_coeff);
      printf("Omegazeta2_coeff     = %.6f\n",pPrec->Omegazeta2_coeff);
      printf("Omegazeta3_coeff     = %.6f\n",pPrec->Omegazeta3_coeff);
      printf("Omegazeta4_coeff     = %.6f\n",pPrec->Omegazeta4_coeff);
      printf("Omegazeta5_coeff     = %.6f\n\n",pPrec->Omegazeta5_coeff);
    #endif

    /* Get psi0 term */
    double psi_of_v0      = 0.0;
    double mm             = 0.0;
    double tmpB           = 0.0;
    double volume_element = 0.0;
    double vol_sign       = 0.0;

    #if DEBUG == 1
      printf("psi1     = %.6f\n",pPrec->psi1);
      printf("psi2     = %.6f\n\n",pPrec->psi2);
      printf("S_0_norm = %.6f\n\n",pPrec->S_0_norm);
    #endif

    /* Tolerance chosen to be consistent with implementation in LALSimInspiralFDPrecAngles */
    if( fabs(pPrec->Smi2 - pPrec->Spl2) < 1.0e-5)
    {
      pPrec->psi0  = 0.0;
    }
    else
    {
      mm      = sqrt( (pPrec->Smi2 - pPrec->Spl2) / (pPrec->S32 - pPrec->Spl2) );
      tmpB    = (pPrec->S_0_norm*pPrec->S_0_norm - pPrec->Spl2) / (pPrec->Smi2 - pPrec->Spl2);

      volume_element  = IMRPhenomX_vector_dot_product( IMRPhenomX_vector_cross_product(L_0,S1v), S2v);
      vol_sign        = (volume_element > 0) - (volume_element < 0);

      psi_of_v0       = IMRPhenomX_psiofv(pPrec->v_0, pPrec->v_0_2, 0.0, pPrec->psi1, pPrec->psi2, pPrec);

      if( tmpB < 0. || tmpB > 1. )
      {
        if(tmpB > 1.0 && (tmpB - 1.) < 0.00001)
        {
          pPrec->psi0 = gsl_sf_ellint_F(asin(vol_sign*sqrt(1.)) , mm, GSL_PREC_DOUBLE ) - psi_of_v0;
        }
        if(tmpB < 0.0 && tmpB > -0.00001)
        {
          pPrec->psi0 = gsl_sf_ellint_F(asin(vol_sign*sqrt(0.)), mm, GSL_PREC_DOUBLE ) - psi_of_v0;
        }
      }
      else
      {
        pPrec->psi0   = gsl_sf_ellint_F(asin( vol_sign * sqrt(tmpB) ), mm, GSL_PREC_DOUBLE) - psi_of_v0;
      }
    }

    #if DEBUG == 1
      printf("psi0_of_v0  = %.6f\n",psi_of_v0);
      printf("tmpB        = %.6f\n",tmpB);
      printf("psi0        = %.6f\n\n",pPrec->psi0);
    #endif

    vector vMSA = {0.,0.,0.};

    double phiz_0     = 0.0;
    UNUSED double phiz_0_MSA = 0.0;

    double zeta_0     = 0.0;
    UNUSED double zeta_0_MSA = 0.0;

    /* Tolerance chosen to be consistent with implementation in LALSimInspiralFDPrecAngles */
    if( fabs(pPrec->Spl2 - pPrec->Smi2) > 1.e-5 )
    {
      vMSA = IMRPhenomX_Return_MSA_Corrections_MSA(pPrec->v_0,pPrec->L_0_norm,pPrec->J_0_norm,pPrec);

      phiz_0_MSA = vMSA.x;
      zeta_0_MSA = vMSA.y;
    }

    // Initial \phi_z
    pPrec->phiz_0 = 0.0;
    phiz_0        = IMRPhenomX_Return_phiz_MSA(pPrec->v_0,pPrec->J_0_norm,pPrec);

    // Initial \zeta
    pPrec->zeta_0 = 0.0;
    zeta_0        = IMRPhenomX_Return_zeta_MSA(pPrec->v_0,pPrec);

    pPrec->phiz_0    = - phiz_0 - vMSA.x;
    pPrec->zeta_0    = - zeta_0 - vMSA.y;

    #if DEBUG == 1
      printf("v_0            = %.6f\n",pPrec->v_0);
      printf("c1             = %.6f\n\n",pPrec->c1);

      printf("eta            = %.6f\n",pPrec->eta);
      printf("eta2           = %.6f\n",pPrec->eta2);
      printf("eta3           = %.6f\n",pPrec->eta3);
      printf("eta4           = %.6f\n",pPrec->eta4);
      printf("ieta           = %.6f\n",pPrec->inveta);
      printf("ieta2          = %.6f\n",pPrec->inveta2);
      printf("ieta3          = %.6f\n\n",pPrec->inveta3);

      printf("SAv            = %.6f\n",pPrec->SAv);
      printf("SAv2           = %.6f\n\n",pPrec->SAv2);
      printf("invSAv2        = %.6f\n\n",pPrec->invSAv2);
      printf("invSAv         = %.6f\n\n",pPrec->invSAv);

      printf("J_0_norm       = %.6f\n",pPrec->J_0_norm);
      printf("L_0_norm       = %.6f\n\n",pPrec->L_0_norm);

      printf("phiz_0         = %.6f\n",pPrec->phiz_0);
      printf("zeta_0         = %.6f\n",pPrec->zeta_0);
      printf("phiz_0_MSA     = %.6f\n",phiz_0_MSA);
      printf("zeta_0_MSA     = %.6f\n\n",zeta_0_MSA);
    #endif

    return XLAL_SUCCESS;
}

double IMRPhenomX_psiofv(const double v, const double v2, const double psi0, const double psi1, const double psi2, const IMRPhenomXPrecessionStruct *pPrec)
{
  /* Equation 51 in arXiv:1703.03967 */
  return ( psi0 - 0.75*pPrec->g0 * pPrec->delta_qq * (1.0 + psi1*v + psi2*v2) / (v2*v) );
}

/**
    Here we solve for the roots of Eq 21 in Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967:
      - Roots for (d S^2)/(d t^2) = -A^2 (S^2 -S+^2)(S^2 - S-^2)(S^2 - S3^2)
      - Returns Spl2 (S+^2), Smi2 (S-^2) and S3^2

      - Note: agrees with independent implementation in Mathematica
*/
vector IMRPhenomX_Return_Roots_MSA(double LNorm, double JNorm, const IMRPhenomXPrecessionStruct *pPrec)
{
  vector vout;

  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  tmp1 = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  tmp4 = 0.0;
  tmp5 = 0.0;
  tmp6 = 0.0;

  vector vBCD;
  vBCD = IMRPhenomX_Return_Spin_Evolution_Coefficients_MSA(LNorm, JNorm, pPrec);

  /* Update struct. Note, that this agreed with independent implementation in Mathematica. */
  const double B  = vBCD.x;
  const double C  = vBCD.y;
  const double D  = vBCD.z;

  const double S1Norm2 = pPrec->S1_norm_2;
  const double S2Norm2 = pPrec->S2_norm_2;

  const double S0Norm2 = pPrec->S_0_norm_2;

  const double B2 = B  * B;
  const double B3 = B2 * B;
  const double BC = B  * C;

  const double p  = C - B2 / 3;
  const double qc = (2.0/27.0)*B3 - BC/3.0 + D;

  const double sqrtarg = sqrt(-p/3.0);
  double acosarg = 1.5 * qc/p/sqrtarg;

  // Make sure that acosarg is appropriately bounded
  if(acosarg < -1)
  {
    acosarg = -1;
  }
  if(acosarg > 1)
  {
    acosarg = +1;
  }
  const double theta     = acos(acosarg) / 3.0;
  const double cos_theta = cos(theta);

  const double dotS1Ln = pPrec->dotS1Ln;
  const double dotS2Ln = pPrec->dotS2Ln;

  double S32, Spl2, Smi2;

  // tmp1 = S32
  // tmp2 = Smi2
  // tmp3 = Spl2
  if(theta != theta || sqrtarg!=sqrtarg || dotS1Ln == 1 || dotS2Ln == 1 || dotS1Ln == -1 || dotS2Ln == -1 || S1Norm2 == 0
  || S2Norm2 == 0)
  {
    S32  = 0.0;

    Smi2 = S0Norm2;

    /*
      Add a numerical perturbation to prevent azimuthal precession angle
      from diverging.
    */

    // Smi2 = S02^2 + epsilon perturbation
    Spl2 = Smi2 + 1e-9;
  }
  else
  {
    /* E.g. see discussion on elliptic functions in arXiv:0711.4064 */
    tmp1  = 2.0*sqrtarg*cos(theta - 2.0*LAL_TWOPI/3.0) - B/3.0;
    tmp2  = 2.0*sqrtarg*cos(theta - LAL_TWOPI/3.0) - B/3.0;
    tmp3  = 2.0*sqrtarg*cos_theta - B/3.0;

    tmp4 = fmax(fmax(tmp1,tmp2),tmp3);
    tmp5 = fmin(fmin(tmp1,tmp2),tmp3);

    // As tmp1 and tmp3 are set by findind the min and max, find the remaining root
    if( (tmp4 - tmp3) > 0.0 && (tmp5 - tmp3) < 0.0)
    {
      tmp6 = tmp3;
    }
    else if( (tmp4 - tmp1) > 0.0 && (tmp5 - tmp1) < 0.0)
    {
      tmp6 = tmp1;
    }
    else
    {
      tmp6 = tmp2;
    }

    /*
        When Spl2 ~ 0 to numerical roundoff then Smi2 can sometimes be ~ negative causing NaN's.
        This occurs in a very limited portion of the parameter space where spins are ~ 0 to numerical roundoff.
        We can circumvent by enforcing +ve definite behaviour when tmp4 ~ 0. Note that S32 can often be negative, this is fine.
    */
    tmp4 = fabs(tmp4);
    tmp6 = fabs(tmp6);

    // Return the roots
    Spl2 = tmp4;
    S32  = tmp5;
    Smi2 = tmp6;
  }

  vout.x = S32;
  vout.y = Smi2;
  vout.z = Spl2;

  return vout;
}

/**
    Get norm of J using Eq 41 of Chatziioannou et al, PRD 95, 104004, (2017)
*/
double IMRPhenomX_JNorm_MSA(const double LNorm, IMRPhenomXPrecessionStruct *pPrec)
{
  const double JNorm2 = (LNorm*LNorm + (2.0 * LNorm * pPrec->c1_over_eta) + pPrec->SAv2);
  return sqrt(JNorm2);
}

/**
    Get norm of S, see PRD 95, 104004, (2017)
*/
double IMRPhenomX_Return_SNorm_MSA(const double v, IMRPhenomXPrecessionStruct *pPrec)
{
  /*
    sn, cn are Jacobi elliptic functions
    psi is the phase and m a parameter entering the Jacobi elliptic functions
  */
  double sn, cn, dn, m, psi;
  double v2 = v*v;

  /*
    If spin norms ~ cancel then we do not need to evaluate the Jacobi elliptic function.
    Check tolerance?
  */
  if( fabs(pPrec->Smi2 - pPrec->Spl2) < 1.0e-5 )
  {
    sn = 0.0;
  }
  else
  {
    /* Equation 25 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967 */
    m   = (pPrec->Smi2 - pPrec->Spl2) / (pPrec->S32 - pPrec->Spl2);

    psi = IMRPhenomX_psiofv(v, v2, pPrec->psi0, pPrec->psi1, pPrec->psi2, pPrec);

    /* Evaluate the Jacobi ellptic functions */
    gsl_sf_elljac_e(psi, m, &sn, &cn, &dn);
  }

  /* Equation 23 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967 */
  const double SNorm2 = pPrec->Spl2 + (pPrec->Smi2 - pPrec->Spl2)*sn*sn;

  return sqrt(SNorm2);
}


/**
    Get coefficients for Eq 21 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
*/
vector IMRPhenomX_Return_Spin_Evolution_Coefficients_MSA(const double LNorm, const double JNorm, const IMRPhenomXPrecessionStruct *pPrec)
{
  vector vout;

  // Total angular momenta: J = L + S1 + S2
  const double JNorm2  = JNorm * JNorm;

  // Orbital angular momenta
  const double LNorm2  = LNorm * LNorm;

  // Dimensionfull spin angular momenta
  const double S1Norm2 = pPrec->S1_norm_2;
  const double S2Norm2 = pPrec->S2_norm_2;

  const double q       = pPrec->qq;
  const double eta     = pPrec->eta;

  const double J2mL2   = (JNorm2 - LNorm2);
  const double J2mL2Sq = J2mL2 * J2mL2;

  const double delta   = pPrec->delta_qq;
  const double deltaSq = delta*delta;


  /*
      Note:
        S_{eff} \equiv \xi = (1 + q)(S1.L) + (1 + 1/q)(S2.L)
  */
  const double Seff    = pPrec->Seff;

  // Note that we do not evaluate Eq. B1 here as it is v dependent whereas B, C and D are not

  // Set Eq. B2, B_coeff
  vout.x = (LNorm2 + S1Norm2)*q + 2.0*LNorm*Seff - 2.0*JNorm2 -
                  S1Norm2 - S2Norm2 + (LNorm2 + S2Norm2)/q;

  // Set Eq. B3, C_coeff
  vout.y = J2mL2Sq - 2.0*LNorm*Seff*J2mL2 - 2.0*((1.0 - q)/q)*LNorm2*(S1Norm2 - q*S2Norm2) +
                      4.0*eta*LNorm2*Seff*Seff - 2.0*delta*(S1Norm2 - S2Norm2)*Seff*LNorm +
                      2.0*((1.0 - q)/q)*(q*S1Norm2 - S2Norm2)*JNorm2;

  // Set Eq. B4, D_coeff
  vout.z = ((1.0 - q)/q)*(S2Norm2 - q*S1Norm2)*J2mL2Sq
                      + deltaSq*(S1Norm2 - S2Norm2)*(S1Norm2 - S2Norm2)*LNorm2/eta
                      + 2.0*delta*LNorm*Seff*(S1Norm2 - S2Norm2)*J2mL2;




  return vout;
}



/**
    Get c constants from Appendix B (B6, B7, B8) of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
*/
vector IMRPhenomX_Return_Constants_c_MSA(const double v, const double JNorm, const IMRPhenomXPrecessionStruct *pPrec)
{
  const double v2 = v*v;
  const double v3 = v*v2;
  const double v4 = v2*v2;
  const double v6 = v3*v3;

  const double JNorm2 = JNorm * JNorm;

  vector vout = {0.,0.,0.};

  const double Seff = pPrec->Seff;

  if(pPrec->IMRPhenomXPrecVersion != 220)
  {
    // Equation B6 of Chatziioannou et al, PRD 95, 104004, (2017)
    vout.x = JNorm * ( 0.75*(1.0 - Seff*v) * v2 * (pPrec->eta3 + 4.0*pPrec->eta3*Seff*v
                  - 2.0*pPrec->eta*(JNorm2 - pPrec->Spl2 + 2.0*(pPrec->S1_norm_2 - pPrec->S2_norm_2)*pPrec->delta_qq)*v2
                  - 4.0*pPrec->eta*Seff*(JNorm2 - pPrec->Spl2)*v3 + (JNorm2 - pPrec->Spl2)*(JNorm2 - pPrec->Spl2)*v4*pPrec->inveta) );

    // Equation B7 of Chatziioannou et al, PRD 95, 104004, (2017)
    vout.y = JNorm * ( -1.5 * pPrec->eta * (pPrec->Spl2 - pPrec->Smi2)*(1.0 + 2.0*Seff*v - (JNorm2 - pPrec->Spl2)*v2*pPrec->inveta2) * (1.0 - Seff*v)*v4 );

    // Equation B8 of Chatziioannou et al, PRD 95, 104004, (2017)
    vout.z = JNorm * ( 0.75 * pPrec->inveta * (pPrec->Spl2 - pPrec->Smi2)*(pPrec->Spl2 - pPrec->Smi2)*(1.0 - Seff * v)*v6 );
  }
  else
  {
    /*  This is as implemented in LALSimInspiralFDPrecAngles, should be equivalent to above code.
        c.f. https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiralFDPrecAngles_internals.c#L578
    */
    double v_2    = v * v;
    double v_3    = v * v_2;
    double v_4    = v_2 * v_2;
    double v_6    = v_2 * v_4;
    double J_norm = JNorm;
    double delta  = pPrec->delta_qq;
    double eta    = pPrec->eta;
    double eta_2  = eta * eta;

    vout.x = -0.75*((JNorm2-pPrec->Spl2)*(JNorm2-pPrec->Spl2)*v_4/(pPrec->eta) - 4.*(pPrec->eta)*(pPrec->Seff)*(JNorm2-pPrec->Spl2)*v_3-2.*(JNorm2-pPrec->Spl2+2*((pPrec->S1_norm_2)-(pPrec->S2_norm_2))*(delta))*(pPrec->eta)*v_2+(4.*(pPrec->Seff)*v+1)*(pPrec->eta)*(eta_2)) *J_norm*v_2*((pPrec->Seff)*v-1.);
    vout.y = 1.5*(pPrec->Smi2-pPrec->Spl2)*J_norm*((JNorm2-pPrec->Spl2)/(pPrec->eta)*v_2-2.*(pPrec->eta)*(pPrec->Seff)*v-(pPrec->eta))*((pPrec->Seff)*v-1.)*v_4;
    vout.z = -0.75*J_norm*((pPrec->Seff)*v-1.)*(pPrec->Spl2-pPrec->Smi2)*(pPrec->Spl2-pPrec->Smi2)*v_6/(pPrec->eta);
  }

  return vout;

}

/**
    Get d constants from Appendix B (B9, B10, B11) of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
*/
vector IMRPhenomX_Return_Constants_d_MSA(const double LNorm, const double JNorm, const IMRPhenomXPrecessionStruct *pPrec)
{

  const double LNorm2 = LNorm * LNorm;
  const double JNorm2 = JNorm * JNorm;

  vector vout = {0.,0.,0.};

  vout.x = -( JNorm2 - (LNorm + pPrec->Spl)*(LNorm + pPrec->Spl))
                  * ( (JNorm2 - (LNorm - pPrec->Spl)*(LNorm - pPrec->Spl)) );
  vout.y = -2.0*(pPrec->Spl2 - pPrec->Smi2)*(JNorm2 + LNorm2 - pPrec->Spl2);
  vout.z = -(pPrec->Spl2 - pPrec->Smi2)*(pPrec->Spl2 - pPrec->Smi2);

  return vout;
}



/**
    Calculate (L dot J)
*/
double IMRPhenomX_costhetaLJ(const double L_norm, const double J_norm, const double S_norm)
{

  double costhetaLJ = 0.5*(J_norm*J_norm + L_norm*L_norm - S_norm*S_norm)/(L_norm * J_norm);

  if (costhetaLJ >  1.0) costhetaLJ = +1.0;
  if (costhetaLJ < -1.0) costhetaLJ = -1.0;

  return costhetaLJ;
}

/**
    Get \f$\psi\f$ using Eq 51 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967:
      - Here \f$\psi\f$ is the phase of S as in Eq 23
      - Note that the coefficients are defined in Appendix C (C1 and C2)
*/
double IMRPhenomX_Return_Psi_MSA(double v, double v2, const IMRPhenomXPrecessionStruct *pPrec)
{
  return ( -0.75 * pPrec->g0 * pPrec->delta_qq * (1.0 + pPrec->psi1*v + pPrec->psi2*v2) / (v2*v) );
}

/**
    Get \f$\dot{\psi}\f$ using Eq 24 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967:
*/
double IMRPhenomX_Return_Psi_dot_MSA(const double v, const IMRPhenomXPrecessionStruct *pPrec)
{
  /*
      \frac{d \psi}{d t} = \frac{A}{2} \sqrt{S+^2 - S3^2}
  */

  const double v2 = v*v;

  const double A_coeff = -1.5 * (v2 * v2 * v2) * (1.0 - v*pPrec->Seff) * pPrec->sqrt_inveta;
  const double psi_dot = 0.5 * A_coeff * sqrt(pPrec->Spl2 - pPrec->S32);

  return psi_dot;
}

/**
    Get \f$\phi_z\f$ using Eq 66 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967:
      - The coefficients are given in Appendix D (D15 - D26)
*/
double IMRPhenomX_Return_phiz_MSA(const double v, const double JNorm, const IMRPhenomXPrecessionStruct *pPrec)
{
  const double invv    = 1.0/v;
  const double invv2   = invv * invv;

  const double LNewt   = (pPrec->eta/v);

  const double c1      = pPrec->c1;
  const double c12     = c1 * c1;

  const double SAv2    = pPrec->SAv2;
  const double SAv     = pPrec->SAv;
  const double invSAv  = pPrec->invSAv;
  const double invSAv2 = pPrec->invSAv2;

  // These are log functions defined in Eq. D27 and D28 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  const double log1    = log( fabs(c1 + JNorm*pPrec->eta   + pPrec->eta*LNewt)   );
  const double log2    = log( fabs(c1 + JNorm*SAv*v + SAv2*v) );

  // Eq. D22 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  const double phiz_0_coeff = (JNorm * pPrec->inveta4) * (0.5*c12 - c1*pPrec->eta2*invv/6.0 - SAv2*pPrec->eta2/3.0 - pPrec->eta4*invv2/3.0)
                      - (c1 * 0.5 * pPrec->inveta)*(c12 * pPrec->inveta4 - SAv2 * pPrec->inveta2)*log1;

  // Eq. D23 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967;
  // Note the factor of c12 in the second term
  const double phiz_1_coeff = - 0.5 * JNorm * pPrec->inveta2 * (c1 + pPrec->eta * LNewt) + 0.5*pPrec->inveta3*(c12 - pPrec->eta2*SAv2)*log1;

  // Eq. D24 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  const double phiz_2_coeff = -JNorm + SAv*log2 - c1*log1*pPrec->inveta;

  // Eq. D25 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  const double phiz_3_coeff = JNorm*v - pPrec->eta*log1 + c1*log2*pPrec->invSAv;

  // Eq. D26 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  const double phiz_4_coeff = (0.5*JNorm*invSAv2*v)*(c1 + v*SAv2) - (0.5*invSAv2*invSAv)*(c12 - pPrec->eta2*SAv2)*log2;

  // Eq. D27 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  const double phiz_5_coeff = -JNorm*v*( (0.5*c12*invSAv2*invSAv2) - (c1*v*invSAv2/6.0) - v*v/3.0 - pPrec->eta2*invSAv2/3.0)
                    + (0.5*c1*invSAv2*invSAv2*invSAv)*(c12 - pPrec->eta2*SAv2)*log2;

  /*
      Eq. 66 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967

      \phi_{z,-1} = \sum^5_{n=0} <\Omega_z>^(n) \phi_z^(n) + \phi_{z,-1}^0

      Note that the <\Omega_z>^(n) are given by pPrec->Omegazn_coeff's as in Eqs. D15-D20
  */
  double phiz_out   = (   phiz_0_coeff*pPrec->Omegaz0_coeff
                              + phiz_1_coeff*pPrec->Omegaz1_coeff
                              + phiz_2_coeff*pPrec->Omegaz2_coeff
                              + phiz_3_coeff*pPrec->Omegaz3_coeff
                              + phiz_4_coeff*pPrec->Omegaz4_coeff
                              + phiz_5_coeff*pPrec->Omegaz5_coeff
                              + pPrec->phiz_0
                            );

  if (phiz_out != phiz_out) phiz_out = 0;

  return (phiz_out);
}

/**
    Get \f$\zeta\f$ using Eq F5 in Appendix F of Chatziioannou et al, PRD 95, 104004, (2017):
*/
double IMRPhenomX_Return_zeta_MSA(const double v, const IMRPhenomXPrecessionStruct *pPrec)
{
  /*
      Eq. F5 of Chatziioannou et al, PRD 95, 104004, (2017)

      \zeta_{z,-1} = \eta v^{-3} \sum^5_{n=0} <\Omega_{\zeta}>^(n) v^(n) + \zeta_{-1}^0

      Note that the <\Omega_{\eta}>^(n) are given by pPrec->Omegazetan_coeff's as in Eqs. F6-F11
  */

  const double invv        = 1.0/v;
  const double invv2       = invv*invv;
  const double invv3       = invv*invv2;
  const double v2          = v*v;
  const double logv        = log(v);

  // \zeta_{z,-1}
  // Note factor of log(v) as per LALSimInspiralFDPrecAngles_internals.c, https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiralFDPrecAngles_internals.c#L718
  double zeta_out = pPrec->eta * (
        pPrec->Omegazeta0_coeff*invv3
      + pPrec->Omegazeta1_coeff*invv2
      + pPrec->Omegazeta2_coeff*invv
      + pPrec->Omegazeta3_coeff*logv
      + pPrec->Omegazeta4_coeff*v
      + pPrec->Omegazeta5_coeff*v2
    )
      + pPrec->zeta_0;

  if (zeta_out != zeta_out) zeta_out = 0;

  // \zeta_{z,-1}
  return zeta_out;
}

/*
    Get MSA corrections to \zeta and \phi_z using Eq. F19 in Appendix F and Eq. 67 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967 respectively:
      -
*/
vector IMRPhenomX_Return_MSA_Corrections_MSA(double v, double LNorm, double JNorm, const IMRPhenomXPrecessionStruct *pPrec)
{
  vector c_vec = {0.,0.,0.};
  vector d_vec = {0.,0.,0.};

  int pflag    = pPrec->IMRPhenomXPrecVersion;

  double v2    = v * v;

  vector vMSA  = {0.,0.,0.};

  // Sets c0, c2 and c4 in pPrec as per Eq. B6-B8 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  c_vec = IMRPhenomX_Return_Constants_c_MSA(v, JNorm, pPrec);

  // Sets d0, d2 and d4 in pPrec as per Eq. B9-B11 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  d_vec = IMRPhenomX_Return_Constants_d_MSA(LNorm, JNorm, pPrec);

  const double c0 = c_vec.x;
  const double c2 = c_vec.y;
  const double c4 = c_vec.z;

  const double d0 = d_vec.x;
  const double d2 = d_vec.y;
  const double d4 = d_vec.z;

  // Pre-cache a bunch of useful variables
  const double two_d0    = 2.0 * d0;

  // Eq. B20 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  const double sd        = sqrt( fabs(d2*d2 - 4.0*d0*d4) );

  // Eq. F20 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  const double A_theta_L = 0.5 * ( (JNorm/LNorm) + (LNorm/JNorm) - (pPrec->Spl2 / (JNorm * LNorm)) );

  // Eq. F21 of Chatziioannou et al, PRD 95, 104004, (2017), arXiv:1703.03967
  const double B_theta_L = 0.5 * pPrec->Spl2mSmi2 / (JNorm * LNorm);

  // Coefficients for B16
  const double nc_num    = 2.0*(d0 + d2 + d4);
  const double nc_denom  = two_d0 + d2 + sd;

  // Equations B16 and B17 respectively
  const double nc        = nc_num   / nc_denom;
  const double nd        = nc_denom / two_d0;

  const double sqrt_nc   = sqrt(fabs(nc));
  const double sqrt_nd   = sqrt(fabs(nd));

  // Get phase and phase evolution of S
  const double psi     = IMRPhenomX_Return_Psi_MSA(v,v2,pPrec) + pPrec->psi0;
  const double psi_dot = IMRPhenomX_Return_Psi_dot_MSA(v,pPrec);

  // Trigonometric calls are expensive, pre-cache them
  // Note: arctan(tan(x)) = 0 if and only if x \in (−pi/2,pi/2).
  const double tan_psi     = tan(psi);
  const double atan_psi    = atan(tan_psi);

  double C1, C2, C2num, C2den;

  // Eq. B18
  C1 = -0.5 * (c0/d0 - 2.0*(c0+c2+c4)/nc_num);

  // Eq. B19
  C2num = c0*( -2.0*d0*d4 + d2*d2 + d2*d4 ) - c2*d0*( d2 + 2.0*d4 ) + c4*d0*( two_d0 + d2 );
  C2den = 2.0 * d0 * sd * (d0 + d2 + d4);
  C2    = C2num / C2den;

  // These are defined in Appendix B, B14 and B15 respectively
  double Cphi, Dphi;
  Cphi = (C1 + C2);
  Dphi = (C1 - C2);

  double phiz_0_MSA_Cphi_term, phiz_0_MSA_Dphi_term;

  /* Calculate C_phi term in Eq. 67 */
  if(nc == 1.0)
  {
    // Limit implied by Eq. 67
    phiz_0_MSA_Cphi_term = 0.0;
  }
  else
  {
    if(pflag == 222 || pflag == 223)
    {
      // As implemented in LALSimInspiralFDPrecAngles.c, c.f. https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiralFDPrecAngles_internals.c#L772
      phiz_0_MSA_Cphi_term = fabs( (c4 * d0 * ((2*d0+d2) + sd) - c2 * d0 * ((d2+2.*d4) - sd) - c0 * ((2*d0*d4) - (d2+d4) * (d2 - sd))) / (C2den)) * (sqrt_nc / (nc - 1.) * (atan_psi - atan(sqrt_nc * tan_psi))) / psi_dot;
    }
    else
    {
      // First term in Eq. 67
      phiz_0_MSA_Cphi_term = ( (Cphi / psi_dot) * sqrt_nc / (nc - 1.0) ) * atan( ( (1.0 - sqrt_nc) * tan_psi ) / ( 1.0 + (sqrt_nc * tan_psi * tan_psi) ) );
    }
  }

  /* Calculate D_phi term in Eq. 67 */
  if(nd == 1.0)
  {
    // Limit implied by Eq. 67
    phiz_0_MSA_Dphi_term = 0.0;
  }
  else
  {
    if(pflag == 222 || pflag == 223)
    {
      // As implemented in LALSimInspiralFDPrecAngles.c, c.f. https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiralFDPrecAngles_internals.c#L779
      phiz_0_MSA_Dphi_term = fabs( (-c4 * d0 * ((2*d0+d2) - sd) + c2 * d0 * ((d2+2.*d4) + sd) - c0 * (-(2*d0*d4) + (d2+d4) * (d2 + sd)))) / (C2den) * (sqrt_nd / (nd - 1.) * (atan_psi - atan(sqrt_nd * tan_psi))) / psi_dot;
    }
    else
    {
      // Second term in Eq. 67
      phiz_0_MSA_Dphi_term = ( (Dphi / psi_dot) * sqrt_nd / (nd - 1.0) ) * atan( ( (1.0 - sqrt_nd) * tan_psi ) / ( 1.0 + (sqrt_nd * tan_psi * tan_psi) ) );
    }
  }

  // Eq. 67
  vMSA.x = (  phiz_0_MSA_Cphi_term + phiz_0_MSA_Dphi_term );

  /*
    The first MSA correction to \zeta as given in Eq. F19
  */
  if(pflag == 222 || pflag == 223 || pflag == 224)
  {
    /*
        As implemented in LALSimInspiralFDPrecAngles.c, c.f. https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiralFDPrecAngles_internals.c#L786
        Note that Cphi and Dphi are *not* used but phiz_0_MSA_Cphi_term and phiz_0_MSA_Dphi_term are
    */
    vMSA.y   =  A_theta_L*vMSA.x +  2.*B_theta_L*d0*(phiz_0_MSA_Cphi_term/(sd-d2) - phiz_0_MSA_Dphi_term/(sd+d2));
  }
  else
  {
    /* Eq. F19 as in arXiv:1703.03967 */
    vMSA.y = ( ( A_theta_L * (Cphi + Dphi) ) + (2.0 * d0 * B_theta_L) * ( ( Cphi / (sd - d2) ) - ( Dphi / (sd + d2) ) ) ) / psi_dot;
  }

  // Return 0 if the angles are NAN
  if (vMSA.x != vMSA.x) vMSA.x = 0;
  if (vMSA.y != vMSA.y) vMSA.y = 0;

  // Obsolete component that we initialize to zero just in case
  vMSA.z = 0.0;

  return vMSA;
}


/**
 * Internal function to computes the PN spin-orbit couplings. As in LALSimInspiralFDPrecAngles.c
 * cf https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiralFDPrecAngles_internals.c#L798
 */
double IMRPhenomX_Get_PN_beta(const double a, const double b, const IMRPhenomXPrecessionStruct *pPrec)
{
    return (pPrec->dotS1L* (a + b*pPrec->qq) + pPrec->dotS2L*(a + b/pPrec->qq));
}

/**
 * Internal function to compute PN spin-spin couplings. As in LALSimInspiralFDPrecAngles.c
 * cf https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiralFDPrecAngles_internals.c#L806
 */
double IMRPhenomX_Get_PN_sigma(const double a, const double b, const IMRPhenomXPrecessionStruct *pPrec)
{
    return pPrec->inveta * (a*pPrec->dotS1S2 - b*pPrec->dotS1L*pPrec->dotS2L);
}

/**
 * Internal function to computes PN spin-spin couplings. As in LALSimInspiralFDPrecAngles.c
 */
double IMRPhenomX_Get_PN_tau(const double a, const double b, const IMRPhenomXPrecessionStruct *pPrec)
{
    return (pPrec->qq * ( (pPrec->S1_norm_2 * a) - b*pPrec->dotS1L*pPrec->dotS1L) + (a*pPrec->S2_norm_2 - b*pPrec->dotS2L*pPrec->dotS2L) / pPrec->qq) / pPrec->eta;
}



/* ~~~~~~~~~~ Vector Utility Functions ~~~~~~~~~~ */
double IMRPhenomX_vector_dot_product(const vector v1, const vector v2)
{
  return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
}


vector IMRPhenomX_vector_cross_product(const vector v1, const vector v2)
{
  vector v3;
  v3.x = v1.y*v2.z - v1.z*v2.y;
  v3.y = v1.z*v2.x - v1.x*v2.z;
  v3.z = v1.x*v2.y - v1.y*v2.x;

  return v3;
}

double IMRPhenomX_vector_L2_norm(const vector v1)
{
  const double dot_product = (v1.x*v1.x) + (v1.y*v1.y) + (v1.z*v1.z);
  return sqrt(dot_product);
}

vector IMRPhenomX_vector_scalar(const vector v1, const double a)
{
    vector v2;
    v2.x = a * v1.x;
    v2.y = a * v1.y;
    v2.z = a * v1.z;

    return v2;
}

vector IMRPhenomX_vector_sum(const vector v1, const vector v2)
{
  vector v3;
  v3.x = v1.x + v2.x;
  v3.y = v1.y + v2.y;
  v3.z = v1.z + v2.z;
  return v3;
}

vector IMRPhenomX_vector_diff(const vector v1, const vector v2)
{
  vector v3;
  v3.x = v1.x - v2.x;
  v3.y = v1.y - v2.y;
  v3.z = v1.z - v2.z;
  return v3;
}

vector IMRPhenomX_vector_PolarToCartesian(const sphpolvector v1)
{
    vector v2;
    const double rsinth = v1.r * sin(v1.theta);

    v2.x = rsinth * cos(v1.phi);
    v2.y = rsinth * sin(v1.phi);
    v2.z = v1.r   * cos(v1.theta);

    return v2;
}

sphpolvector IMRPhenomX_vector_CartesianToPolar(const vector v1)
{
    sphpolvector v2;

    v2.r     = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
    v2.theta = acos(v1.z / v2.r);
    v2.phi   = XLALSimIMRPhenomXatan2tol(v1.y,v1.x,1e-15);

    return v2;
}

/** Function to rotate vector about z axis by given angle */
vector IMRPhenomX_vector_rotate_z(const REAL8 angle, const vector v1)
{
  vector v2;

  v2.x = v1.x * cos(angle) - v1.y * sin(angle);
  v2.y = v1.x * sin(angle) + v1.y * cos(angle);
  v2.z = v1.z;

  return v2;
}

/** Function to rotate vector about y axis by given angle */
vector IMRPhenomX_vector_rotate_y(const REAL8 angle, const vector v1)
{
  vector v2;

  v2.x = + v1.x * cos(angle) + v1.z * sin(angle);
  v2.y =   v1.y;
  v2.z = - v1.x * sin(angle) + v1.z * cos(angle);

  return v2;
}


/** Function to rotate vector about z axis by given angle */
void IMRPhenomX_rotate_z(const REAL8 angle, REAL8 *vx, REAL8 *vy, REAL8 *vz)
{
  const REAL8 tmpx = *vx;
  const REAL8 tmpy = *vy;
  const REAL8 tmpz = *vz;

  REAL8  cosa=cos(angle), sina=sin(angle);

  REAL8 tmp1 = tmpx*cosa - tmpy*sina;
  REAL8 tmp2 = tmpx*sina + tmpy*cosa;

  *vx = tmp1;
  *vy = tmp2;
  *vz = tmpz;
}

/** Function to rotate vector about y axis by given angle */
void IMRPhenomX_rotate_y(REAL8 angle, REAL8 *vx, REAL8 *vy, REAL8 *vz)
{
  const REAL8 tmpx = *vx;
  const REAL8 tmpy = *vy;
  const REAL8 tmpz = *vz;

  REAL8  cosa=cos(angle), sina=sin(angle);

  REAL8 tmp1 = + tmpx*cosa + tmpz*sina;
  REAL8 tmp2 = - tmpx*sina + tmpz*cosa;

  *vx = tmp1;
  *vy = tmpy;
  *vz = tmp2;
}

REAL8 IMRPhenomX_Cartesian_to_SphericalPolar_theta(const double x, const double y, const UNUSED double z)
{
  REAL8 theta;
  const REAL8 norm = sqrt(x*x + y*y + z*z);

  /* Note that we add a numerical offset of 1e-16 to prevent accidental division by zero */
  theta = acos(z / (norm + 1e-16));

  return theta;
}

REAL8 IMRPhenomX_Cartesian_to_SphericalPolar_phi(const double x, const double y, const UNUSED double z)
{
  REAL8 phi;

  phi = XLALSimIMRPhenomXatan2tol(y,x,1e-16);

  return phi;
}

vector IMRPhenomX_vector_PolarToCartesian_components(const REAL8 mag, const REAL8 theta, const REAL8 phi)
{
  vector v1;

  const double rsintheta = mag * sin(theta);

  v1.x = rsintheta * cos(phi);
  v1.y = rsintheta * sin(phi);
  v1.z = mag * cos(theta);

  return v1;
}


/** used in numerical evaluation of Euler angles */
static REAL8TimeSeries *appendTS(REAL8TimeSeries *start, REAL8TimeSeries *end) {
    UINT4 origlen = start->data->length;
    start = XLALResizeREAL8TimeSeries(start, 0,
            start->data->length + end->data->length - 1);

    memcpy(start->data->data + origlen -2, end->data->data,
            (end->data->length)*sizeof(REAL8));

    XLALGPSAdd(&(start->epoch), -end->deltaT*(end->data->length - 1));
    XLALDestroyREAL8TimeSeries(end);


    return start;
}


/** Analytical continuation for alpha angle in MRD */
int alphaMRD_coeff(gsl_spline spline_alpha, gsl_interp_accel accel_alpha, double fmaxPN, IMRPhenomXWaveformStruct *pWF, PhenomXPalphaMRD *alpha_params){

    int success = GSL_SUCCESS;

    double ftrans = fmaxPN;
    double f1 = 0.97 *ftrans ;
    double f2 = 0.99 *ftrans ;
    double f1sq = f1*f1, f2sq = f2*f2;
    double f1cube = f1sq*f1;

    double alpha1, alpha2, dalpha1;

    success = gsl_spline_eval_e(&spline_alpha, f1, &accel_alpha,&alpha1);
    alpha1 = -alpha1;
    if(success != GSL_SUCCESS)
        {
        XLALPrintError("XLAL Error - %s: Alpha could not be interpolated at f=%.5f\n",__func__,XLALSimIMRPhenomXUtilsMftoHz(f1,pWF->Mtot));

        }

    success = success + gsl_spline_eval_deriv_e (&spline_alpha, f1, &accel_alpha,&dalpha1);
    dalpha1 = -dalpha1;
    if(success != GSL_SUCCESS)
        {
        XLALPrintError("XLAL Error - %s: dalpha/df could not be interpolated at f=%.5f\n",__func__,XLALSimIMRPhenomXUtilsMftoHz(f1,pWF->Mtot));

        }


    success = success + gsl_spline_eval_e (&spline_alpha, f2, &accel_alpha,&alpha2);
    alpha2 = -alpha2;
    if(success != GSL_SUCCESS)
        {
        XLALPrintError("XLAL Error - %s: Alpha could not be interpolated at f=%.5f\n",__func__,XLALSimIMRPhenomXUtilsMftoHz(f2,pWF->Mtot));

        }

    double aC=0., bC=0., cC=0.;

    if(success == GSL_SUCCESS){

     aC = (f1cube*(f1 - f2)*(f1 + f2)*dalpha1 + 2*(pow(f1,4) - 2*f1sq*f2sq)*alpha1 + 2*pow(f2,4)*alpha2)/(2.*pow(f1sq - f2sq,2));

     bC = (pow(f1,4)*f2sq*(f1*(f1 - f2)*(f1 + f2)*dalpha1 + 2*f2sq*(-alpha1 + alpha2)))/(2.*pow(f1sq - f2sq,2));

     cC = (f1sq*(f1*(-pow(f1,4) + pow(f2,4))*dalpha1 + 4*pow(f2,4)*(alpha1 - alpha2)))/(2.*pow(f1sq - f2sq,2));

    }
    alpha_params->aRD = aC;
    alpha_params->bRD = bC;
    alpha_params->cRD = cC;

    return success;


          }

double alphaMRD(double Mf, PhenomXPalphaMRD *alpha_params){

    double Mf2 = Mf*Mf;
    double Mf4 = Mf2*Mf2;
    double minalpha=alpha_params->aRD + alpha_params->bRD/Mf4 + alpha_params->cRD/Mf2;

    return (-minalpha);

          }

double dalphaMRD(double Mf, PhenomXPalphaMRD *alpha_params){

    double Mf2 = Mf*Mf;
    double Mf3 = Mf2*Mf;
    double Mf5 = Mf2*Mf3;

    return (4.*alpha_params->bRD/Mf5+2.*alpha_params->cRD/Mf3);

          }

/** Function to determine coefficients of analytical continuation of beta through MRD */
int betaMRD_coeff(gsl_spline spline_cosb, gsl_interp_accel accel_cosb, double fmaxPN, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPrecessionStruct *pPrec){

    int success = GSL_SUCCESS;
    double fdamp = pWF->fDAMP;

    QNMFits *qnms = (QNMFits *) XLALMalloc(sizeof(QNMFits));
    IMRPhenomXHM_Initialize_QNMs(qnms);
    pPrec->beta_params->dfdamp = qnms->fdamp_lm[0](pWF->afinal)/pWF->Mfinal-fdamp;
    LALFree(qnms);

    double kappa = 2.*LAL_PI*pPrec->beta_params->dfdamp;

    double f1 = 0.97 *fmaxPN ;
    double f2 = 0.98 *fmaxPN ;
    double f1sq = f1*f1, f2sq = f2*f2 ;
    double ef1 = exp(kappa*f1), ef2 = exp(kappa*f2);

    double cosbeta1;
    success = gsl_spline_eval_e(&spline_cosb, f1, &accel_cosb, &cosbeta1);

    double dcosbeta2;
    success=gsl_spline_eval_deriv_e (&spline_cosb, f2, &accel_cosb,&dcosbeta2);

    double cosbeta2;
    success = gsl_spline_eval_e(&spline_cosb, f2, &accel_cosb,&cosbeta2);

    double cosbetamax;
    success = gsl_spline_eval_e(&spline_cosb, pPrec->fmax_inspiral, &accel_cosb, &cosbetamax);

    double aC,bC,cC,dC;

    if(fabs(cosbeta1)>1 || fabs(cosbeta2)>1 || fabs(cosbetamax)>1||success!=GSL_SUCCESS)
     {

       aC=0.;
       bC=0.;
       cC=0.;
       dC=0. ;
       pPrec->beta_params->flat_RD=true;
       success = GSL_SUCCESS;

     }

    else{

    double beta1 =  acos(cosbeta1);
    double beta2 =  acos(cosbeta2);
    double sqrtarg = 1.-cosbeta2*cosbeta2;
    double dbeta2 = - dcosbeta2/sqrt((sqrtarg <= 0. ? 1. : sqrtarg));


    double off = (cosbetamax < 0. ? LAL_PI : 0.);

    aC= (-(ef1*pow(f1,4)*(off - beta1)) + ef2*pow(f2,3)*(f2*(-f1 + f2)*dbeta2 + (-(f2*(3 + f2*kappa)) + f1*(4 + f2*kappa))*(off - beta2)))/pow(f1 - f2,2);

    bC =(2*ef1*pow(f1,4)*f2*(off - beta1) + ef2*pow(f2,3)*((f1 - f2)*f2*(f1 + f2)*dbeta2 - (-(f2sq*(2 + f2*kappa)) + f1sq*(4 + f2*kappa))*(off - beta2)))/pow(f1 - f2,2);

    cC =(-(ef1*pow(f1,4)*f2sq*(off - beta1)) + ef2*f1*pow(f2,4)*(f2*(-f1 + f2)*dbeta2 + (-(f2*(2 + f2*kappa)) + f1*(3 + f2*kappa))*(off - beta2)))/pow(f1 - f2,2);

    dC = off;


    pPrec->beta_params->flat_RD=false;

     }

    pPrec->beta_params->aRD = aC;
    pPrec->beta_params->bRD = bC;
    pPrec->beta_params->cRD = cC;
    pPrec->beta_params->dRD = dC;

    pPrec->beta_params->cosbeta_sign =  copysign(1.0, cosbetamax);

    return success;

          }



double betaMRD(double Mf, UNUSED IMRPhenomXWaveformStruct *pWF,PhenomXPbetaMRD *beta_params){

    double beta=0.;

    if(beta_params->flat_RD)
        beta=acos(beta_params->cosbeta_sign);
    else
    {
    double kappa = 2.*LAL_PI*beta_params->dfdamp;
    double Mf2 = Mf*Mf;
    double Mf3 = Mf2* Mf;
        beta = exp(-Mf*kappa)/Mf*(beta_params->aRD/Mf+beta_params->bRD/(Mf2)+beta_params->cRD/Mf3)+beta_params->dRD;

    }

    return(beta);

          }


/* Integrate minimal rotation condition to compute gamma once alpha and beta are known
   Uses Boole's rule
*/
int gamma_from_alpha_cosbeta(double *gamma, double Mf, double deltaMf,IMRPhenomXWaveformStruct *pWF,IMRPhenomXPrecessionStruct *pPrec){

    REAL8 Mf_high = Mf+deltaMf;
    REAL8 alphadoti, Mf_aux, gammai;
    REAL8 step = (Mf_high-Mf)*0.25;

    double alphadotcosbeta[5], cosbeta_aux[5];

    int status_alpha=GSL_SUCCESS, status_beta=GSL_SUCCESS;

    if(Mf<=pPrec->ftrans_MRD)
    {
        for(UINT4 jdx=0; jdx<5; jdx++){

            Mf_aux=Mf+jdx*step;
            status_beta=gsl_spline_eval_e(pPrec->cosbeta_spline, Mf_aux, pPrec->cosbeta_acc,&cosbeta_aux[jdx]);
            status_alpha=gsl_spline_eval_deriv_e(pPrec->alpha_spline, Mf_aux, pPrec->alpha_acc,&alphadoti);
            XLAL_CHECK(status_alpha == GSL_SUCCESS && status_beta==GSL_SUCCESS, XLAL_EFUNC, "Error in %s: could not evaluate splines for alpha and/or gamma angles.\n",__func__);
            alphadotcosbeta[jdx]=cosbeta_aux[jdx]*alphadoti;
                                      }

    }

    else
    {
        for(UINT4 jdx=0; jdx<5; jdx++){

            Mf_aux=Mf+jdx*step;
            cosbeta_aux[jdx] = cos(betaMRD(Mf_aux,pWF,pPrec->beta_params));
            alphadoti = dalphaMRD(Mf_aux,pPrec->alpha_params);
            alphadotcosbeta[jdx]=cosbeta_aux[jdx]*alphadoti;
                                  }
    }

    gammai= -2.*step/45.*(7.*alphadotcosbeta[0]+32.*alphadotcosbeta[1]+12.*alphadotcosbeta[2]+32.*alphadotcosbeta[3]+7.*alphadotcosbeta[4]);
    *gamma = gammai;

    return(status_beta+status_alpha);


}

/** This function evaluates the SpinTaylor Euler angles on a frequency grid passed by the user. Used in LALSimIMRPhenomX.c. */
int IMRPhenomXPSpinTaylorAnglesIMR(
       REAL8Sequence **alphaFS,              /**< [out] Alpha angle frequency series [out] */
       REAL8Sequence **cosbetaFS,            /**< [out]  cos(Beta) angle frequency series [out] */
       REAL8Sequence **gammaFS,              /**< [out] Gamma angle frequency series [out] */
       REAL8Sequence *freqsIN,               /**< [in]  Frequency grid on which Euler angles  will be evaluated [in]*/
       IMRPhenomXWaveformStruct *pWF,        /**< [in] Waveform structure [in]*/
       IMRPhenomXPrecessionStruct *pPrec,    /**< [in] Precession structure [in]*/
       LALDict *LALparams        /**< LAL Dictionary struct */
       )
     {

         int status = XLAL_SUCCESS;
         REAL8 fRef= pWF->fRef;
         REAL8 fmin = freqsIN->data[0];
         REAL8 fmax = freqsIN->data[freqsIN->length-1];
         // some useful quantities in geom units
         REAL8 Mfmin=XLALSimIMRPhenomXUtilsHztoMf(fmin,pWF->Mtot);

         /* Sanity checks */
         if(fmin    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "fmin must be positive.\n");                          }
         if(fmax    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "fmax must be positive.\n");                          }
         if(fmax    <= fmin) { XLAL_ERROR(XLAL_EDOM, "fmax must be larger than fmin.\n");                          }
         if(fRef < fmin){ XLAL_ERROR(XLAL_EDOM, "fRef must be >= fmin.\n"); }

         size_t output_length = freqsIN->length;

         //Evaluate splines for alpha and cosbeta
         status=IMRPhenomX_InterpolateAlphaBeta_SpinTaylor(pWF,pPrec,LALparams);
         XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Error in %s: IMRPhenomX_InterpolateAlphaBeta_SpinTaylor failed.\n",__func__);

         //Evaluate splines for gamma
         status = IMRPhenomX_InterpolateGamma_SpinTaylor(fmin,fmax,pWF,pPrec);
         XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Error in %s: IMRPhenomX_InterpolateGamma_SpinTaylor failed.\n",__func__);

         REAL8 alphamin=0., cosbetamin=0.;

         status = gsl_spline_eval_e(pPrec->alpha_spline, Mfmin, pPrec->alpha_acc,&alphamin);
         if( status != GSL_SUCCESS)
         {
         XLALPrintError("Error in %s: could not evaluate alpha(f_min). Got alpha=%.4f \n",__func__,alphamin);
         XLAL_ERROR(XLAL_EFUNC);}

         status = gsl_spline_eval_e(pPrec->cosbeta_spline, Mfmin, pPrec->cosbeta_acc,&cosbetamin);
         if ( status != GSL_SUCCESS)
         {XLALPrintError("Error in %s: could not evaluate cosbeta(f_min). Got cosbeta=%.4f \n",__func__,cosbetamin);
         XLAL_ERROR(XLAL_EFUNC);}

         // check if we can evaluate alpha(fref) using PN or MRD analytical continuation

         if(pWF->MfRef<=pPrec->ftrans_MRD)
             status=gsl_spline_eval_e(pPrec->alpha_spline, pWF->MfRef, pPrec->alpha_acc, &pPrec->alpha_ref);
         else if(pPrec->IMRPhenomXPrecVersion==320||pPrec->IMRPhenomXPrecVersion==321)
             pPrec->alpha_ref = alphaMRD(pWF->MfRef,pPrec->alpha_params);
         else
             status=gsl_spline_eval_e(pPrec->alpha_spline, pPrec->ftrans_MRD, pPrec->alpha_acc, &pPrec->alpha_ref);

         XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Error in %s: could not evaluate alpha(f_ref).\n",__func__);

         *alphaFS = XLALCreateREAL8Sequence(output_length);
         *cosbetaFS = XLALCreateREAL8Sequence(output_length);
         *gammaFS = XLALCreateREAL8Sequence(output_length);

         REAL8 alphaOff = pPrec->alpha0;
         pPrec->alpha_offset=-pPrec->alpha_ref+alphaOff;

         (*gammaFS)->data[0] = 0.;
         (*alphaFS)->data[0] = alphamin-pPrec->alpha_ref+alphaOff;
         (*cosbetaFS)->data[0] = cosbetamin;


         REAL8 alphai=0., cosbetai=0., gammai=0.;
         REAL8 Mf;


         for( UINT4 i = 1; i < output_length; i++ ){

             Mf=XLALSimIMRPhenomXUtilsHztoMf(freqsIN->data[i],pWF->Mtot);
             REAL8 deltaMF =XLALSimIMRPhenomXUtilsHztoMf(freqsIN->data[i]-freqsIN->data[i-1],pWF->Mtot);

             if(Mf<pPrec->ftrans_MRD)

                     {

                         alphai = gsl_spline_eval(pPrec->alpha_spline, Mf, pPrec->alpha_acc);
                         cosbetai = gsl_spline_eval(pPrec->cosbeta_spline, Mf, pPrec->cosbeta_acc);
                         gammai = gsl_spline_eval(pPrec->gamma_spline, Mf, pPrec->gamma_acc);

                         (*alphaFS)->data[i] = alphai+pPrec->alpha_offset;
                         (*cosbetaFS)->data[i] = cosbetai;
                         (*gammaFS)-> data[i] = gammai;

                     }


             else {

                 if(pPrec->IMRPhenomXPrecVersion==320 || pPrec->IMRPhenomXPrecVersion==321 ){

                     (*alphaFS)->data[i]=alphaMRD(Mf,pPrec->alpha_params)+pPrec->alpha_offset;
                     REAL8 beta_MRD=betaMRD(Mf,pWF,pPrec->beta_params);
                     (*cosbetaFS)->data[i]=cos(beta_MRD);
                     REAL8 deltagamma = 0.;
                     status = gamma_from_alpha_cosbeta(&deltagamma,Mf, deltaMF,pWF,pPrec);
                     if(status==GSL_SUCCESS)
                     (*gammaFS) -> data[i] = (*gammaFS)-> data[i-1]+deltagamma;
                     else (*gammaFS) -> data[i] = (*gammaFS)-> data[i-1];
                     }

                 // just continue angles into RD with constant values

                 else{
                     (*alphaFS)->data[i]=(*alphaFS)->data[i-1];
                     (*cosbetaFS)->data[i]=(*cosbetaFS)-> data[i-1];
                     (*gammaFS) -> data[i]=(*gammaFS)-> data[i-1];
                     }
             }

         }

         return status;

     }

/**  This function builds and stores splines for  \f$\alpha\f$  and \f$\cos\beta\f$ in the frequency range covered by PN, and computes a spline for \f$\gamma\f$ between fmin and fmax   */
int IMRPhenomX_SpinTaylorAnglesSplinesAll(
       REAL8 fmin,                           /**< [in]  Minimum frequency of the gamma spline [in]*/
       REAL8 fmax,                           /**< [in]  Maximum frequency of the gamma spline [in]*/
       IMRPhenomXWaveformStruct *pWF,        /**< [in] Waveform structure [in]*/
       IMRPhenomXPrecessionStruct *pPrec,    /**< [in] Precession structure [in]*/
       LALDict *LALparams        /**< LAL Dictionary struct */
       )
     {

         int status = XLAL_SUCCESS;
         REAL8 fRef= pWF->fRef;

         /* Sanity checks */
         if(fmin    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "fmin must be positive.\n");                          }
         if(fmax    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "fmax must be positive.\n");                          }
         if(fmax    <= fmin) { XLAL_ERROR(XLAL_EDOM, "fmax must be larger than fmin.\n");                          }
         if(fRef < fmin){ XLAL_ERROR(XLAL_EDOM, "fRef must be >= fmin.\n"); }

         //Evaluate splines for alpha and cosbeta
         status=IMRPhenomX_InterpolateAlphaBeta_SpinTaylor(pWF,pPrec,LALparams);
         XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Error in %s: IMRPhenomX_InterpolateAlphaBeta_SpinTaylor failed.\n",__func__);

         //Evaluate splines for gamma
         status = IMRPhenomX_InterpolateGamma_SpinTaylor(fmin,fmax,pWF,pPrec);
         XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Error in %s: IMRPhenomX_InterpolateGamma_SpinTaylor failed.\n",__func__);

         // check if we can evaluate alpha(fref) using PN or an analytical continuation

         if(pWF->MfRef<=pPrec->ftrans_MRD)
             status=gsl_spline_eval_e(pPrec->alpha_spline, pWF->MfRef, pPrec->alpha_acc, &pPrec->alpha_ref);
         else if(pPrec->IMRPhenomXPrecVersion==320||pPrec->IMRPhenomXPrecVersion==321)
             pPrec->alpha_ref = alphaMRD(pWF->MfRef,pPrec->alpha_params);
         else
             status=gsl_spline_eval_e(pPrec->alpha_spline, pPrec->ftrans_MRD, pPrec->alpha_acc, &pPrec->alpha_ref);

         XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Error in %s: could not evaluate alpha(f_ref).\n",__func__);

         return status;

     }

/** This function computes gamma from the minimal rotation condition and stores a spline for it */
int IMRPhenomX_InterpolateGamma_SpinTaylor(
       REAL8 fmin,               /**< starting frequency (Hz) */
       REAL8 fmax,               /**< maximum frequency (Hz) */
       IMRPhenomXWaveformStruct *pWF,        /**< [in] Waveform structure [in]*/
       IMRPhenomXPrecessionStruct *pPrec    /**< [in] Precession structure [in]*/
       )
     {

         int status = XLAL_SUCCESS;
         REAL8 fRef= pWF->fRef;

         // some useful quantities in geom units
         REAL8 Mfmin=XLALSimIMRPhenomXUtilsHztoMf(fmin,pWF->Mtot);

         /* Sanity checks */
         if(fmin    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "fmin must be positive.\n");                          }
         if(fmax    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "fmax must be positive.\n");                          }
         if(fmax    <= fmin) { XLAL_ERROR(XLAL_EDOM, "fmax must be larger than fmin.\n");                          }
         if(fRef < fmin){ XLAL_ERROR(XLAL_EDOM, "fRef must be >= fmin.\n"); }

         REAL8 seglen=XLALSimInspiralChirpTimeBound(pWF->fRef, pWF->m1_SI, pWF->m2_SI, pWF->chi1L,pWF->chi2L);
         REAL8 deltaFv1= 1./MAX(4.,pow(2, ceil(log(seglen)/log(2))));
         REAL8 deltaF = MIN(deltaFv1,0.1);
         REAL8 deltaMF = XLALSimIMRPhenomXUtilsHztoMf(deltaF,pWF->Mtot);
         // length of frequency series expected by the user
         size_t iStart = (size_t) (fmin / deltaF);
         size_t iStop  = (size_t) (fmax / deltaF) + 1;
         size_t output_length = iStop-iStart;
          if(output_length<4)
         {
         XLALPrintError("Error in %s: no. of points is insufficient for spline interpolation of gamma",__func__);
         XLAL_ERROR(XLAL_EFUNC);
         }


         REAL8Sequence *frequencies = NULL;
         frequencies=XLALCreateREAL8Sequence(output_length);
         frequencies->data[0]=Mfmin;

         REAL8Sequence *gamma_array = NULL;
         gamma_array=XLALCreateREAL8Sequence(output_length);
         gamma_array->data[0]=0.;
         REAL8 Mf;
         int gamma_status=GSL_SUCCESS;

         for( UINT4 i = 1; i < output_length; i++ )
         {

             Mf=Mfmin+i*deltaMF;
             frequencies->data[i]=Mf;
             REAL8 deltagamma=0.;

             if(Mf<pPrec->ftrans_MRD){
                         gamma_status=gamma_from_alpha_cosbeta(&deltagamma,Mf, deltaMF,pWF,pPrec);
                         gamma_array->data[i]=gamma_array->data[i-1]+deltagamma;
                     }
             else {
                 if(pPrec->IMRPhenomXPrecVersion==320 || pPrec->IMRPhenomXPrecVersion==321 ){
                     gamma_status = gamma_from_alpha_cosbeta(&deltagamma, Mf, deltaMF,pWF,pPrec);
                     gamma_array->data[i]=gamma_array->data[i-1]+deltagamma;
                     }
                 // just continue angles into RD with constant values
                 else{
                     gamma_array->data[i]=gamma_array->data[i-1];
                     }
                  }

                  if(gamma_status!=GSL_SUCCESS) status=gamma_status;

         }

         if(status==GSL_SUCCESS)
         {
            pPrec->gamma_acc = gsl_interp_accel_alloc();
            pPrec->gamma_spline = gsl_spline_alloc(gsl_interp_cspline, output_length);
            gsl_spline_init(pPrec->gamma_spline, frequencies->data, gamma_array -> data, output_length);
            pPrec->gamma_ref=gsl_spline_eval(pPrec->gamma_spline, pWF->MfRef, pPrec->gamma_acc);
         }

         else{

          gsl_spline_free(pPrec->alpha_spline);
          gsl_interp_accel_free(pPrec->alpha_acc);

          gsl_spline_free(pPrec->cosbeta_spline);
          gsl_interp_accel_free(pPrec->cosbeta_acc);


         }


         XLALDestroyREAL8Sequence(frequencies);
         XLALDestroyREAL8Sequence(gamma_array);

         return status;

     }



/** This function computes cubic splines of the alpha and beta inspiral Euler angles, which are then stored into a IMRPhenomXPrecessionStruct structure.
 - If the user passed PhenomXPFinalSpinMod=4, the function corrects the estimate for the final precessing spin based on the result of the PN integration.
 - For versions 32*, the function also computes the parameters needed to obtain a smooth MRD continuation of alpha/beta.
 - The memory allocated for the PN arrays is freed at the end of this function.
*/
int IMRPhenomX_InterpolateAlphaBeta_SpinTaylor(
            IMRPhenomXWaveformStruct *pWF,        /**< [in] Waveform structure [in]*/
            IMRPhenomXPrecessionStruct *pPrec,    /**< [in] Precession structure [in]*/
            LALDict *LALparams        /**< LAL Dictionary struct */
            )
{

              int status = XLAL_SUCCESS;

              size_t lenPN = pPrec->PNarrays->V_PN->data->length;

              REAL8Sequence *fgw =NULL ;
              /* Setup sequences for angles*/
              REAL8Sequence *alpha = NULL;
              REAL8Sequence *alphaaux = NULL;
              REAL8Sequence *cosbeta = NULL;


              fgw=XLALCreateREAL8Sequence(lenPN);
              alpha=XLALCreateREAL8Sequence(lenPN);
              alphaaux=XLALCreateREAL8Sequence(lenPN);
              cosbeta=XLALCreateREAL8Sequence(lenPN);


              REAL8 fgw_Mf, fgw_Hz, Mfmax_PN=0.;
              // i_max is used to discard possibly unphysical points in the calculation of the final spin
              UINT8 i_max=0;
              REAL8 LNhatx_temp,LNhaty_temp,LNhatz_temp;

              for(UINT8 i=0; i < lenPN; i++){

                  LNhatx_temp = (pPrec->PNarrays->LNhatx_PN->data->data[i]);
                  LNhaty_temp = (pPrec->PNarrays->LNhaty_PN->data->data[i]);
                  LNhatz_temp = (pPrec->PNarrays->LNhatz_PN->data->data[i]);

                  IMRPhenomX_rotate_z(-pPrec->phiJ_Sf,  &LNhatx_temp, &LNhaty_temp, &LNhatz_temp);
                  IMRPhenomX_rotate_y(-pPrec->thetaJ_Sf, &LNhatx_temp, &LNhaty_temp, &LNhatz_temp);
                  IMRPhenomX_rotate_z(-pPrec->kappa,  &LNhatx_temp, &LNhaty_temp, &LNhatz_temp);

                  fgw_Hz= pow(pPrec->PNarrays->V_PN->data->data[i],3.)/pPrec->piGM;
                  fgw_Mf= XLALSimIMRPhenomXUtilsHztoMf(fgw_Hz,pWF->Mtot);

                  if(fgw_Hz>0.){

                  /* Compute Euler angles in the J frame */
                  alphaaux->data[i] = atan2(LNhaty_temp, LNhatx_temp);
                  cosbeta->data[i] = LNhatz_temp;
                  fgw->data[i] = fgw_Mf;

                  Mfmax_PN = fgw_Mf;
                  i_max = i;
                  }

                  else
                      break;


              }

            REAL8 fmax_inspiral = Mfmax_PN-pWF->deltaMF;
            if(fmax_inspiral > pWF->fRING-pWF->fDAMP) fmax_inspiral = 1.020 * pWF->fMECO;

            pPrec->ftrans_MRD = 0.98*fmax_inspiral;
            pPrec->fmax_inspiral= fmax_inspiral;

            // Interpolate alpha
            XLALSimIMRPhenomXUnwrapArray(alphaaux->data, alpha->data, lenPN);
            pPrec->alpha_acc = gsl_interp_accel_alloc();
            pPrec->alpha_spline = gsl_spline_alloc(gsl_interp_cspline, lenPN);


            status = gsl_spline_init(pPrec->alpha_spline, fgw->data, alpha->data, lenPN);

            if (status != GSL_SUCCESS)
            {
                 XLALPrintError("Error in %s: error in computing gsl spline for alpha.\n",__func__);
            }

            // Interpolate cosbeta
            pPrec->cosbeta_acc = gsl_interp_accel_alloc();
            pPrec->cosbeta_spline = gsl_spline_alloc(gsl_interp_cspline, lenPN);
            status =gsl_spline_init(pPrec->cosbeta_spline, fgw->data, cosbeta->data, lenPN);

            if (status != GSL_SUCCESS)
            {
                 XLALPrintError("Error in %s: error in computing gsl spline for cos(beta).\n",__func__);
            }

            REAL8 cosbetamax;

            status = gsl_spline_eval_e(pPrec->cosbeta_spline, fmax_inspiral, pPrec->cosbeta_acc,&cosbetamax);
            if(status != GSL_SUCCESS)
            {
                XLALPrintError("Error in %s: error in computing cosbeta.\n",__func__);
            }

            // estimate final spin using spins at the end of the PN integration

            if(XLALSimInspiralWaveformParamsLookupPhenomXPFinalSpinMod(LALparams)==4){

            REAL8 m1 = pWF->m1_SI / pWF->Mtot_SI;
            REAL8 m2 = pWF->m2_SI / pWF->Mtot_SI;

            vector Lnf  = {pPrec->PNarrays->LNhatx_PN->data->data[i_max],pPrec->PNarrays->LNhaty_PN->data->data[i_max],pPrec->PNarrays->LNhatz_PN->data->data[i_max]};
            REAL8 Lnorm = sqrt(IMRPhenomX_vector_dot_product(Lnf,Lnf));
            vector S1f  = {pPrec->PNarrays->S1x_PN->data->data[i_max],pPrec->PNarrays->S1y_PN->data->data[i_max],pPrec->PNarrays->S1z_PN->data->data[i_max]};
            vector S2f  = {pPrec->PNarrays->S2x_PN->data->data[i_max],pPrec->PNarrays->S2y_PN->data->data[i_max],pPrec->PNarrays->S2z_PN->data->data[i_max]};


            REAL8 dotS1L = IMRPhenomX_vector_dot_product(S1f,Lnf)/Lnorm;
            REAL8 dotS2L  = IMRPhenomX_vector_dot_product(S2f,Lnf)/Lnorm;
            vector S1_perp = IMRPhenomX_vector_diff(S1f,IMRPhenomX_vector_scalar(Lnf, dotS1L));
            S1_perp = IMRPhenomX_vector_scalar(S1_perp,m1*m1);
            vector S2_perp = IMRPhenomX_vector_diff(S2f,IMRPhenomX_vector_scalar(Lnf, dotS2L));
            S2_perp = IMRPhenomX_vector_scalar(S2_perp,m2*m2);
            vector Stot_perp = IMRPhenomX_vector_sum(S1_perp,S2_perp);
            REAL8 S_perp_norm = sqrt(IMRPhenomX_vector_dot_product(Stot_perp,Stot_perp));
            REAL8 chi_perp_norm = S_perp_norm *pow(m1 + m2,2)/pow(m1,2);

            pWF->afinal= copysign(1.0, cosbetamax)* XLALSimIMRPhenomXPrecessingFinalSpin2017(pWF->eta,dotS1L,dotS2L,chi_perp_norm);

            pWF->fRING     = evaluate_QNMfit_fring22(pWF->afinal) / (pWF->Mfinal);
            pWF->fDAMP     = evaluate_QNMfit_fdamp22(pWF->afinal) / (pWF->Mfinal);
            }

            // initialize parameters for RD continuation
            pPrec->alpha_params    = XLALMalloc(sizeof(PhenomXPalphaMRD));
            pPrec->beta_params    = XLALMalloc(sizeof(PhenomXPbetaMRD));

            if(pPrec->IMRPhenomXPrecVersion==320 || pPrec->IMRPhenomXPrecVersion==321){

            status = alphaMRD_coeff(*pPrec->alpha_spline, *pPrec->alpha_acc, pPrec->fmax_inspiral, pWF, pPrec->alpha_params);
            if(status!=XLAL_SUCCESS) XLALPrintError("XLAL Error in %s: error in computing parameters for MRD continuation of Euler angles.\n",__func__);


            status = betaMRD_coeff(*pPrec->cosbeta_spline, *pPrec->cosbeta_acc, pPrec->fmax_inspiral, pWF, pPrec);
             if(status!=XLAL_SUCCESS) XLALPrintError("XLAL Error in %s: error in computing parameters for MRD continuation of Euler angles.\n",__func__);

            }


            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->V_PN);
            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->S1x_PN);
            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->S1y_PN);
            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->S1z_PN);
            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->S2x_PN);
            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->S2y_PN);
            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->S2z_PN);
            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->LNhatx_PN);
            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->LNhaty_PN);
            XLALDestroyREAL8TimeSeries(pPrec->PNarrays->LNhatz_PN);
            XLALFree(pPrec->PNarrays);



            XLALDestroyREAL8Sequence(fgw);
            XLALDestroyREAL8Sequence(alphaaux);
            XLALDestroyREAL8Sequence(cosbeta);
            XLALDestroyREAL8Sequence(alpha);

            if(status != GSL_SUCCESS){

            gsl_spline_free(pPrec->alpha_spline);
            gsl_spline_free(pPrec->cosbeta_spline);
            gsl_interp_accel_free(pPrec->alpha_acc);
            gsl_interp_accel_free(pPrec->cosbeta_acc);


            }



            return status;

          }




/** Wrapper of  XLALSimInspiralSpinTaylorPNEvolveOrbit : if integration is successful, stores arrays containing PN solution in  a PhenomXPInspiralArrays struct  */
int IMRPhenomX_InspiralAngles_SpinTaylor(
            PhenomXPInspiralArrays *arrays, /**< [out] Struct containing solutions returned by PNEvolveOrbit   */
            REAL8 chi1x,   /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) at fRef */
            REAL8 chi1y,   /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) at fRef */
            REAL8 chi1z,   /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) at fRef */
            REAL8 chi2x,   /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) at fRef */
            REAL8 chi2y,   /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) at fRef */
            REAL8 chi2z,   /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) at fRef */
            REAL8 fmin,    /**< minimum frequency (Hz) */
            int PrecVersion, /**< precessing version (int) */
            IMRPhenomXWaveformStruct *pWF,        /**< Waveform structure [in]*/
            LALDict *LALparams        /**< LAL Dictionary struct */
            )
{

              int status = XLAL_SUCCESS;
              REAL8 fRef=pWF->fRef;

              /* Sanity checks */
              if(fmin    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "fmin must be positive.\n");                          }
              if(fRef < fmin){ XLAL_ERROR(XLAL_EDOM, "fRef must be >= fmin.\n"); }


              REAL8 m1_SI=pWF->m1_SI,m2_SI=pWF->m2_SI;
              REAL8 s1x=chi1x, s1y=chi1y, s1z=chi1z;
              REAL8 s2x=chi2x, s2y=chi2y, s2z=chi2z;

              UNUSED REAL8 piGM = LAL_PI * (pWF->m1_SI + pWF->m2_SI) * (LAL_G_SI / LAL_C_SI) / (LAL_C_SI * LAL_C_SI);

              REAL8 quadparam1=pWF->quadparam1;
              REAL8 quadparam2=pWF->quadparam2;
              REAL8 lambda1=pWF->lambda1;
              REAL8 lambda2=pWF->lambda2;


              // Tidal parameters are preloaded in aligned spin waveform structure. If using BH values in the twisting up, overwrite these. Note that the aligned spin model (i.e. the co-precessing modes) will still be using the NS tidal parameters
              if (PrecVersion==311 || PrecVersion==321)
              {
                quadparam1 = 1.;
                quadparam2 = 1.;
                lambda1 = 0.;
                lambda2 = 0.;
              }


              int phaseO = XLALSimInspiralWaveformParamsLookupPNPhaseOrder(LALparams);
              int spinO = XLALSimInspiralWaveformParamsLookupPNSpinOrder(LALparams);
              int tideO = XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALparams);
              int lscorr = XLALSimInspiralWaveformParamsLookupLscorr(LALparams);
              if(lscorr!=0) XLAL_PRINT_WARNING("IMRPhenomXP with SpinTaylor angles was only reviewed for lscorr=0.\n");


              // enforce default orders if user does not specify any: this is added to ensure backward compatibility should default PN orders change in the future
              if(phaseO == -1)  XLALSimInspiralWaveformParamsInsertPNPhaseOrder(LALparams,7);
              if(spinO == -1)  XLALSimInspiralWaveformParamsInsertPNSpinOrder(LALparams,6);
              if(tideO == -1)  XLALSimInspiralWaveformParamsInsertPNTidalOrder(LALparams,12);

               /* Initialize needed time series */
              REAL8TimeSeries *V = NULL;
              REAL8TimeSeries *Phi = NULL;
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



              REAL8 lnhatx,lnhaty,lnhatz, e1y,e1z,e1x;
              lnhatx = lnhaty = e1y = e1z = 0;
              lnhatz = e1x = 1.;

              // if the user does not specify any SpinTaylor approximant, default to SpinTaylorT4
              const char * approx_name=XLALSimInspiralWaveformParamsLookupPhenomXPSpinTaylorVersion(LALparams);
              if(approx_name==NULL)
                  approx_name="SpinTaylorT4";
              Approximant approx = XLALSimInspiralGetApproximantFromString(approx_name);

              REAL8 fS,fE;

              REAL8 fMECO_Hz=XLALSimIMRPhenomXUtilsMftoHz(pWF->fMECO,pWF->Mtot);

              // for versions 32* to work, we need to start integration in the inspiral
               if(fmin>fMECO_Hz &&(PrecVersion==320 || PrecVersion==321))
              {
                  fmin=fMECO_Hz;
              }

              REAL8 fCut = XLALSimIMRPhenomXUtilsMftoHz(pWF->fRING+8.*pWF->fDAMP,pWF->Mtot);
              int n;

              REAL8 coarse_fac = (float)XLALSimInspiralWaveformParamsLookupPhenomXPSpinTaylorCoarseFactor(LALparams);
              if(coarse_fac  < 1) { XLAL_ERROR(XLAL_EDOM, "Coarse factor must be >= 1!\n");}

              REAL8 deltaT_coarse = 0.5*coarse_fac/(fCut);

              fS=fmin;
              fE=fCut;


              if( fRef < LAL_REAL4_EPS  || fabs(fRef - fmin) < LAL_REAL4_EPS )

              {
                  fRef = fmin;

                  n=XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi, &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z, deltaT_coarse, m1_SI, m2_SI,fS,fE,s1x,s1y,s1z,s2x,s2y,s2z,lnhatx,lnhaty,lnhatz,e1x,e1y,e1z,lambda1,lambda2,quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approx);

                  if( n < 0 )
                          XLAL_ERROR(XLAL_EFUNC);
              }


              else if((fRef - fmin) > LAL_REAL4_EPS )
              {
                  /* Integrate backward to fStart */
                  fS = fRef;
                  fE = fmin-0.5;
                  n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi,
                          &S1x, &S1y, &S1z, &S2x, &S2y, &S2z,
                          &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z,
                          deltaT_coarse, m1_SI, m2_SI, fS, fE, s1x, s1y, s1z, s2x, s2y,
                          s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                  quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approx);


                  if( n < 0 )  XLAL_ERROR(XLAL_EFUNC);


                  if(V->data->length>1){

                  REAL8TimeSeries *V2=NULL, *Phi2=NULL, *S1x2=NULL, *S1y2=NULL, *S1z2=NULL, *S2x2=NULL, *S2y2=NULL, *S2z2=NULL;
                  REAL8TimeSeries *LNhatx2=NULL, *LNhaty2=NULL, *LNhatz2=NULL, *E1x2=NULL, *E1y2=NULL, *E1z2=NULL;


                  /* Integrate forward to end of waveform */
                  fS = fRef;
                  fE = fCut;
                  n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V2, &Phi2,
                          &S1x2, &S1y2, &S1z2, &S2x2, &S2y2, &S2z2,
                          &LNhatx2, &LNhaty2, &LNhatz2, &E1x2, &E1y2, &E1z2,
                          deltaT_coarse, m1_SI, m2_SI, fS, fE, s1x, s1y, s1z, s2x, s2y,
                          s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                  quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approx);
                  if( n < 0 )  XLAL_ERROR(XLAL_EFUNC);



                  // Stitch 2nd set of vectors onto 1st set.
                  V = appendTS(V, V2);
                  Phi = appendTS(Phi, Phi2);
                  S1x = appendTS(S1x, S1x2);
                  S1y = appendTS(S1y, S1y2);
                  S1z = appendTS(S1z, S1z2);
                  S2x = appendTS(S2x, S2x2);
                  S2y = appendTS(S2y, S2y2);
                  S2z = appendTS(S2z, S2z2);
                  LNhatx = appendTS(LNhatx, LNhatx2);
                  LNhaty = appendTS(LNhaty, LNhaty2);
                  LNhatz = appendTS(LNhatz, LNhatz2);
                  E1x = appendTS(E1x, E1x2);
                  E1y = appendTS(E1y, E1y2);
                  E1z = appendTS(E1z, E1z2);
                  }

                  else{

                XLALDestroyREAL8TimeSeries( V );
                XLALDestroyREAL8TimeSeries( Phi );
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
                // the failure will be caught by GetAndSetPrecessionVariables. If the integration fails, the waveform generation will default to prec. version 223
                return(XLAL_FAILURE);

                       }

                }

              else
              {
                  XLAL_ERROR(XLAL_EFUNC);
              }

              size_t copyLength;

              if(coarse_fac>1)
              {
              // at high frequencies, perform a fine integration
              // stop some bins before the last one

              int lenLow=V->data->length;
              int nbuffer=MIN(9,lenLow-1);

              if(lenLow-1-nbuffer<0) nbuffer=lenLow-1;

              copyLength=lenLow-1-nbuffer;

              REAL8 vtrans = V->data->data[lenLow-1-nbuffer];
              REAL8 ftrans = pow(vtrans,3.)/piGM;

              REAL8 LNhatx_trans=LNhatx->data->data[lenLow-1-nbuffer];
              REAL8 LNhaty_trans=LNhaty->data->data[lenLow-1-nbuffer];
              REAL8 LNhatz_trans=LNhatz->data->data[lenLow-1-nbuffer];

              REAL8 E1x_trans, E1y_trans, E1z_trans;
              E1x_trans = E1x->data->data[lenLow-1-nbuffer];
              E1y_trans = E1y->data->data[lenLow-1-nbuffer];
              E1z_trans = E1z->data->data[lenLow-1-nbuffer];

              REAL8 S1x_trans, S1y_trans, S1z_trans, S2x_trans, S2y_trans, S2z_trans;
              S1x_trans = S1x->data->data[lenLow-1-nbuffer];
              S1y_trans = S1y->data->data[lenLow-1-nbuffer];
              S1z_trans = S1z->data->data[lenLow-1-nbuffer];

              S2x_trans = S2x->data->data[lenLow-1-nbuffer];
              S2y_trans = S2y->data->data[lenLow-1-nbuffer];
              S2z_trans = S2z->data->data[lenLow-1-nbuffer];

              fS=ftrans;
              fE=fCut;
              REAL8 deltaT = 0.5/(fCut);

              XLAL_CHECK(fS > 0., XLAL_EFUNC, "Error: Transition frequency in PN integration is not positive.\n");

              REAL8TimeSeries *Phi_PN=NULL, *E1x_PN=NULL, *E1y_PN=NULL, *E1z_PN=NULL;

              n=XLALSimInspiralSpinTaylorPNEvolveOrbit(&arrays->V_PN, &Phi_PN, &arrays->S1x_PN, &arrays->S1y_PN, &arrays->S1z_PN, &arrays->S2x_PN, &arrays->S2y_PN, &arrays->S2z_PN, &arrays->LNhatx_PN, &arrays->LNhaty_PN, &arrays->LNhatz_PN, &E1x_PN, &E1y_PN, &E1z_PN, deltaT, m1_SI, m2_SI,fS,fE,S1x_trans,S1y_trans,S1z_trans,S2x_trans,S2y_trans,S2z_trans,LNhatx_trans,LNhaty_trans,LNhatz_trans,E1x_trans, E1y_trans, E1z_trans,lambda1,lambda2,quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approx);

              // free high frequency arrays
              XLALDestroyREAL8TimeSeries( Phi_PN );
              XLALDestroyREAL8TimeSeries( E1x_PN );
              XLALDestroyREAL8TimeSeries( E1y_PN );
              XLALDestroyREAL8TimeSeries( E1z_PN );

              if( n < 0 )
                          XLAL_ERROR(XLAL_EFUNC);

              size_t lenPN=lenLow-nbuffer-1+arrays->V_PN->data->length;

              if(lenPN < 4) {
                XLALPrintError("Error in %s: no. of points is insufficient for spline interpolation",__func__);
                XLAL_ERROR(XLAL_EFUNC);
				}

              // resize series to full length (coarse+fine)
              XLALResizeREAL8TimeSeries(arrays->V_PN,-(lenLow-nbuffer-1),lenPN);
              XLALResizeREAL8TimeSeries(arrays->LNhatx_PN,-(lenLow-nbuffer-1),lenPN);
              XLALResizeREAL8TimeSeries(arrays->LNhaty_PN,-(lenLow-nbuffer-1),lenPN);
              XLALResizeREAL8TimeSeries(arrays->LNhatz_PN,-(lenLow-nbuffer-1),lenPN);
              XLALResizeREAL8TimeSeries(arrays->S1x_PN,-(lenLow-nbuffer-1),lenPN);
              XLALResizeREAL8TimeSeries(arrays->S1y_PN,-(lenLow-nbuffer-1),lenPN);
              XLALResizeREAL8TimeSeries(arrays->S1z_PN,-(lenLow-nbuffer-1),lenPN);
              XLALResizeREAL8TimeSeries(arrays->S2x_PN,-(lenLow-nbuffer-1),lenPN);
              XLALResizeREAL8TimeSeries(arrays->S2y_PN,-(lenLow-nbuffer-1),lenPN);
              XLALResizeREAL8TimeSeries(arrays->S2z_PN,-(lenLow-nbuffer-1),lenPN);

              }

              else{

                copyLength=V->data->length-1;
                if(copyLength < 4) {
                XLALPrintError("Error in %s: no. of points is insufficient for spline interpolation",__func__);
                XLAL_ERROR(XLAL_EFUNC);
				}

                LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO;

                /* allocate memory for output vectors */
                arrays->V_PN = XLALCreateREAL8TimeSeries( "PN_EXPANSION_PARAMETER", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);
                arrays->S1x_PN = XLALCreateREAL8TimeSeries( "SPIN1_X_COMPONENT", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);
                arrays->S1y_PN = XLALCreateREAL8TimeSeries( "SPIN1_Y_COMPONENT", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);
                arrays->S1z_PN = XLALCreateREAL8TimeSeries( "SPIN1_Z_COMPONENT", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);
                arrays->S2x_PN = XLALCreateREAL8TimeSeries( "SPIN2_X_COMPONENT", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);
                arrays->S2y_PN = XLALCreateREAL8TimeSeries( "SPIN2_Y_COMPONENT", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);
                arrays->S2z_PN = XLALCreateREAL8TimeSeries( "SPIN2_Z_COMPONENT", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);
                arrays->LNhatx_PN = XLALCreateREAL8TimeSeries( "LNHAT_X_COMPONENT", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);
                arrays->LNhaty_PN = XLALCreateREAL8TimeSeries( "LNHAT_Y_COMPONENT", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);
                arrays->LNhatz_PN = XLALCreateREAL8TimeSeries( "LNHAT_Z_COMPONENT", &ligotimegps_zero, 0.,
                        deltaT_coarse, &lalDimensionlessUnit, copyLength);



              }


              // copy coarse-grid data into fine-grid arrays
              memcpy(arrays->V_PN->data->data,V->data->data,(copyLength)*sizeof(REAL8));
              memcpy(arrays->LNhatx_PN->data->data,LNhatx->data->data,(copyLength)*sizeof(REAL8));
              memcpy(arrays->LNhaty_PN->data->data,LNhaty->data->data,(copyLength)*sizeof(REAL8));
              memcpy(arrays->LNhatz_PN->data->data,LNhatz->data->data,(copyLength)*sizeof(REAL8));
              memcpy(arrays->S1x_PN->data->data,S1x->data->data,(copyLength)*sizeof(REAL8));
              memcpy(arrays->S1y_PN->data->data,S1y->data->data,(copyLength)*sizeof(REAL8));
              memcpy(arrays->S1z_PN->data->data,S1z->data->data,(copyLength)*sizeof(REAL8));
              memcpy(arrays->S2x_PN->data->data,S2x->data->data,(copyLength)*sizeof(REAL8));
              memcpy(arrays->S2y_PN->data->data,S2y->data->data,(copyLength)*sizeof(REAL8));
              memcpy(arrays->S2z_PN->data->data,S2z->data->data,(copyLength)*sizeof(REAL8));


              XLALDestroyREAL8TimeSeries( V );
              XLALDestroyREAL8TimeSeries( Phi );
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

              // check that the first frequency node returned is indeed below the fmin requested, to avoid interpolation errors. If not return an error which will trigger the fallback to MSA
              REAL8 fminPN=pow(arrays->V_PN->data->data[0],3.)/piGM;
              if(fminPN<0.||fminPN>fmin) return(XLAL_FAILURE);

              return status;

          }



/** Wrapper of IMRPhenomX_SpinTaylorAnglesSplinesAll: fmin and fmax are determined by the function based on the mode content and binary's parameters . Used in LALSimIMRPhenomXPHM.c */
int IMRPhenomX_Initialize_Euler_Angles(
    IMRPhenomXWaveformStruct *pWF,
    IMRPhenomXPrecessionStruct *pPrec,
    LALDict *lalParams
)

{
      int status = XLAL_SUCCESS;
      REAL8 thresholdPMB  = XLALSimInspiralWaveformParamsLookupPhenomXPHMThresholdMband(lalParams);

      // start below fMin to avoid interpolation artefacts
      REAL8 buffer = (pWF->deltaF>0.) ? 3.*pWF->deltaF : 0.5;
      REAL8 fminAngles = (pWF->fMin-buffer)*2./pPrec->M_MAX;
      // check we still pass a meaningful fmin
      XLAL_CHECK(fminAngles > 0., XLAL_EFUNC, "Error - %s: fMin is too low and numerical angles could not be computed.\n",__func__);

      // If MB is on, we take advantage of the fact that we can compute angles on an array

      if(thresholdPMB>0.)
        pPrec->Mfmax_angles = pWF->fRING+4.*pWF->fDAMP;
      else
        pPrec->Mfmax_angles = (MAX(pWF->MfMax,pWF->fRING+4.*pWF->fDAMP)+XLALSimIMRPhenomXUtilsHztoMf(buffer,pWF->Mtot))*2./pPrec->M_MIN;
      REAL8 fmaxAngles = XLALSimIMRPhenomXUtilsMftoHz(pPrec->Mfmax_angles,pWF->Mtot);

      // we add a few bins to fmax to make sure we do not run into interpolation errors
      status = IMRPhenomX_SpinTaylorAnglesSplinesAll(fminAngles,fmaxAngles,pWF,pPrec,lalParams);
      XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "%s: IMRPhenomX_SpinTaylorAnglesSplinesAll failed.",__func__);

      status = gsl_spline_eval_e(pPrec->alpha_spline, pPrec->ftrans_MRD, pPrec->alpha_acc,&pPrec->alpha_ftrans);
      XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "%s: could not compute alpha et the end of inspiral.",__func__);

      status = gsl_spline_eval_e(pPrec->cosbeta_spline, pPrec->ftrans_MRD, pPrec->cosbeta_acc,&pPrec->cosbeta_ftrans);
      XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "%s: could not compute cosbeta et the end of inspiral.",__func__);

      status = gsl_spline_eval_e(pPrec->gamma_spline, pPrec->ftrans_MRD, pPrec->gamma_acc,&pPrec->gamma_ftrans);
      XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "%s: could not compute gamma et the end of inspiral.",__func__);

      return status;

}





/** XLAL function that  evaluates the SpinTaylor Euler angles on a frequency grid passed by the user. Used in LALSimIMRPhenomX.c. */
int XLALSimIMRPhenomXPSpinTaylorAngles(
            REAL8Sequence **alphaFS,              /**< [out] Alpha angle frequency series [out] */
            REAL8Sequence **cosbetaFS,            /**< [out]  cos(Beta) angle frequency series [out] */
            REAL8Sequence **gammaFS,              /**< [out] Gamma angle frequency series [out] */
            REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
            REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
            REAL8 s1x,                /**< x component of primary spin*/
            REAL8 s1y,                /**< y component of primary spin*/
            REAL8 s1z,                /**< z component of primary spin */
            REAL8 s2x,                /**< x component of secondary spin*/
            REAL8 s2y,                /**< y component of secondary spin*/
            REAL8 s2z,                /**< z component of secondary spin */
            REAL8 fmin,               /**< starting GW frequency (Hz) */
            REAL8 fmax,               /**< maximum GW frequency (Hz) */
            REAL8 deltaF,             /**< starting GW frequency (Hz) */
            REAL8 fRef,               /**< reference GW frequency (Hz) */
            REAL8 phiRef,               /**< reference orbital phase (rad) */
            LALDict *LALparams        /**< LAL Dictionary struct */
            )
          {

              int status = XLAL_SUCCESS;
              int pversion = XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(LALparams);

              /* Sanity checks */
              if(fRef  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef must be positive or set to 0 to ignore.\n");  }
              if(deltaF   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaF must be positive.\n");                         }
              if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1_SI must be positive.\n");                             }
              if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2_SI must be positive.\n");                             }
              if(fmin    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "fmin must be positive.\n");                          }
              if(fmax    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "fmax must be positive.\n");                          }
              if(fRef > 0.0 && fRef < fmin){ XLAL_ERROR(XLAL_EDOM, "fRef must be >= fmin or =0 to use fmin.\n"); }

              // length of frequency series expected by the user
              size_t iStart = (size_t) (fmin / deltaF);
              size_t iStop  = (size_t) (fmax / deltaF) + 1;
              size_t output_length = iStop-iStart;

              status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI,&m2_SI,&s1x,&s1y,&s1z,&s2x,&s2y,&s2z);
              XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

              /* Initialize IMR PhenomX Waveform struct and check that it initialized correctly */
               IMRPhenomXWaveformStruct *pWF;
               pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
               status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, s1z, s2z, deltaF, fRef, phiRef, fmin, fmax, 1e6*LAL_PC_SI, 0., LALparams, 0);
               XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

              IMRPhenomXPrecessionStruct *pPrec;
              pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
              pPrec->M_MIN = 2, pPrec->M_MAX = 2;
              status = IMRPhenomXGetAndSetPrecessionVariables(pWF,pPrec,m1_SI,m2_SI,s1x,s1y,s1z,s2x,s2y,s2z,LALparams,0);
                   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetAndSetPrecessionVariables failed.\n");

              // if version has been changed because of fallback, raise error
              XLAL_CHECK(pPrec->IMRPhenomXPrecVersion == pversion, XLAL_EFUNC, "Error: %s failed.\n",__func__);

              //Evaluate splines for alpha and cosbeta
              status=IMRPhenomX_InterpolateAlphaBeta_SpinTaylor(pWF,pPrec,LALparams);
              XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Error in %s: IMRPhenomX_InterpolateAlphaBeta_SpinTaylor failed.\n",__func__);

              REAL8 alphamin=0., cosbetamin=0.;
              int success=XLAL_SUCCESS;
              REAL8 Mfmin=XLALSimIMRPhenomXUtilsHztoMf(fmin,pWF->Mtot);

              success = gsl_spline_eval_e(pPrec->alpha_spline, Mfmin, pPrec->alpha_acc,&alphamin);
              success = success + gsl_spline_eval_e(pPrec->cosbeta_spline, Mfmin, pPrec->cosbeta_acc,&cosbetamin);
              XLAL_CHECK(XLAL_SUCCESS == success, XLAL_EFUNC, "Error: %s: could not evaluate angles at fMin.\n",__func__);


              // check if we can evaluate alpha(fref) using PN or an analytical continuation

              if(pWF->MfRef<pPrec->ftrans_MRD)
                  success=gsl_spline_eval_e(pPrec->alpha_spline, pWF->MfRef, pPrec->alpha_acc, &pPrec->alpha_ref);
              else if(pPrec->IMRPhenomXPrecVersion==320||pPrec->IMRPhenomXPrecVersion==321)
                  pPrec->alpha_ref = alphaMRD(pWF->MfRef,pPrec->alpha_params);
              else
                  success=gsl_spline_eval_e(pPrec->alpha_spline, pPrec->ftrans_MRD, pPrec->alpha_acc, &pPrec->alpha_ref);

              XLAL_CHECK(XLAL_SUCCESS == success, XLAL_EFUNC, "Error: %s: could not evaluate angles at fRef.\n",__func__);



              *alphaFS = XLALCreateREAL8Sequence(output_length);
              *cosbetaFS = XLALCreateREAL8Sequence(output_length);
              *gammaFS = XLALCreateREAL8Sequence(output_length);


              REAL8 alphaOff = pPrec->alpha0;
              pPrec->alpha_offset=-pPrec->alpha_ref+alphaOff;

              // determine offset for gamma et the end
              (*gammaFS)->data[0] = 0.;
              (*alphaFS)->data[0] = alphamin-pPrec->alpha_ref+alphaOff;
              (*cosbetaFS)->data[0] = cosbetamin;


              REAL8 alphai=0., cosbetai=0., gammai=0.;
              REAL8 Mf;
              REAL8Sequence *frequencies = NULL;
              frequencies=XLALCreateREAL8Sequence(output_length);
              frequencies->data[0]=Mfmin;

              for( UINT4 i = 1; i < output_length; i++ ){

                  Mf=Mfmin+i*pWF->deltaMF;
                  frequencies->data[i]=Mf;
                  REAL8 deltagamma=0.;

                  if(Mf<pPrec->ftrans_MRD)


                          {

                              success = gsl_spline_eval_e(pPrec->alpha_spline, Mf, pPrec->alpha_acc,&alphai);
                              success = success + gsl_spline_eval_e(pPrec->cosbeta_spline, Mf, pPrec->cosbeta_acc,&cosbetai);
                              success = success + gamma_from_alpha_cosbeta(&deltagamma,Mf, pWF->deltaMF,pWF,pPrec);
                               if(success != XLAL_SUCCESS)
                                    XLALPrintError("%s: Interpolation of SpinTaylor angles failed at f=%.5f)\n",__func__,XLALSimIMRPhenomXUtilsMftoHz(Mf,pWF->Mtot));
                              gammai=(*gammaFS)->data[i-1]+deltagamma;

                              (*alphaFS)->data[i] = alphai+pPrec->alpha_offset;
                              (*cosbetaFS)->data[i] = cosbetai;
                              (*gammaFS)-> data[i] = gammai;


                          }


                  else {

                      if(pPrec->IMRPhenomXPrecVersion==320 || pPrec->IMRPhenomXPrecVersion==321){

                          (*alphaFS)->data[i]=alphaMRD(Mf,pPrec->alpha_params)+pPrec->alpha_offset;
                          REAL8 beta_MRD=betaMRD(Mf,pWF,pPrec->beta_params);
                          (*cosbetaFS)->data[i]=cos(beta_MRD);
                          status=gamma_from_alpha_cosbeta(&deltagamma,Mf, pWF->deltaMF,pWF,pPrec);
                           if(success != XLAL_SUCCESS)
                                    XLALPrintError("%s: could not integrate minimal rotation condition at f=%.5f)\n",__func__,XLALSimIMRPhenomXUtilsMftoHz(Mf,pWF->Mtot));
                          gammai=(*gammaFS)->data[i-1]+deltagamma;
                          (*gammaFS)-> data[i] =gammai;

                          }

                      // just continue angles into RD with constant values

                      else{
                          (*alphaFS)->data[i]=(*alphaFS)->data[i-1];
                          (*cosbetaFS)->data[i]=(*cosbetaFS)-> data[i-1];
                          (*gammaFS) -> data[i]=(*gammaFS)-> data[i-1];
                          }
                  }


              }


              pPrec->gamma_acc = gsl_interp_accel_alloc();
              pPrec->gamma_spline = gsl_spline_alloc(gsl_interp_cspline, output_length);
              gsl_spline_init(pPrec->gamma_spline, frequencies->data, (*gammaFS) -> data, output_length);
              pPrec->gamma_ref=gsl_spline_eval(pPrec->gamma_spline, pWF->MfRef, pPrec->gamma_acc);

              for( UINT4 i = 0; i < output_length; i++ ) (*gammaFS) -> data[i]=(*gammaFS)-> data[i]-pPrec->gamma_ref-pPrec->epsilon0;

              LALFree(pPrec->alpha_params);
              LALFree(pPrec->beta_params);

              gsl_spline_free(pPrec->alpha_spline);
              gsl_interp_accel_free(pPrec->alpha_acc);

              gsl_spline_free(pPrec->cosbeta_spline);
              gsl_interp_accel_free(pPrec->cosbeta_acc);

              gsl_spline_free(pPrec->gamma_spline);
              gsl_interp_accel_free(pPrec->gamma_acc);


              LALFree(pWF);
              LALFree(pPrec);
              XLALDestroyREAL8Sequence(frequencies);


              return status;

          }



void IMRPhenomX_GetandSetModes(LALValue *ModeArray,IMRPhenomXPrecessionStruct *pPrec){

    INT2Sequence *modeseq;
    modeseq=XLALSimInspiralModeArrayReadModes(ModeArray);

    int nmodes=modeseq->length/2;
    float M_MAX=1., M_MIN=4.;
    for(int jj=0; jj<nmodes; jj++)
    {
      if(modeseq->data[2*jj+1]>M_MAX) M_MAX=(float)(modeseq->data[2*jj+1]);
      if(abs(modeseq->data[2*jj+1])<M_MIN) M_MIN=(float)abs(modeseq->data[2*jj+1]);
      }

    XLALDestroyINT2Sequence(modeseq);

    pPrec->M_MIN = M_MIN; pPrec->M_MAX=M_MAX;
    return;

}

/** Core twisting up routine for SpinTaylor angles */
int IMRPhenomXPTwistUp22_NumericalAngles(
  const COMPLEX16 hAS,    /**< Underlying aligned-spin IMRPhenomXAS strain */
  REAL8 alpha, /**< cosbeta Euler angle series */
  REAL8 cos_beta, /**< cosbeta Euler angle series */
  REAL8 gamma, /**< gamma Euler angle series */
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomXP Precession Struct */
  COMPLEX16 *hp,                            /**< [out] h_+ polarization \f$\tilde h_+\f$ */
  COMPLEX16 *hc                             /**< [out] h_x polarization \f$\tilde h_x\f$ */

)
{
  XLAL_CHECK(hp  != NULL, XLAL_EFAULT);
  XLAL_CHECK(hc  != NULL, XLAL_EFAULT);


  double cBetah      = 0.0;
  double sBetah      = 0.0;

  double epsilon=-(gamma-pPrec->gamma_ref)+pPrec->epsilon0;

  if(pPrec->IMRPhenomXPrecVersion!=310 && pPrec->IMRPhenomXPrecVersion!=311 && pPrec->IMRPhenomXPrecVersion!=320 && pPrec->IMRPhenomXPrecVersion!=321) XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXPTwistUp22_NumericalAngles, incorrect precessing version passed (must be 310, 311, 320,or 321)!\n");

  INT4 status = 0;
  status = IMRPhenomXWignerdCoefficients_cosbeta(&cBetah, &sBetah, cos_beta);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "Call to IMRPhenomXWignerdCoefficients_cosbeta failed.");

  /* Useful powers of the Wigner coefficients */
  const REAL8 cBetah2 = cBetah * cBetah;
  const REAL8 cBetah3 = cBetah * cBetah2;
  const REAL8 cBetah4 = cBetah * cBetah3;
  const REAL8 sBetah2 = sBetah * sBetah;
  const REAL8 sBetah3 = sBetah * sBetah2;
  const REAL8 sBetah4 = sBetah * sBetah3;

  /*
      Compute the Wigner d coefficients, see Appendix A of arXiv:2004.06503
        d22  = Table[WignerD[{2, mp, 2}, 0, -\[Beta], 0], {mp, -2, 2}]
        d2m2 = Table[WignerD[{2, mp, -2}, 0, -\[Beta], 0], {mp, -2, 2}]
  */
  // d22  = {d^2_{-2,2} , d^2_{-1,2}, d^2_{0,2}, d^2_{1,2}, d^2_{2,2} }
  const COMPLEX16 d22[5]   = {sBetah4, 2.0*cBetah*sBetah3, pPrec->sqrt6*sBetah2*cBetah2, 2.0*cBetah3*sBetah, cBetah4};

  // d2m2 = {d^2_{-2,-2} , d^2_{-1,-2}, d^2_{0,-2}, d^2_{1,-2}, d^2_{2,-2} }
  const COMPLEX16 d2m2[5]  = {d22[4], -d22[3], d22[2], -d22[1], d22[0]}; /* Exploit symmetry d^2_{-2,m} = (-1)^m d^2_{2,-m}*/

  const COMPLEX16 Y2mA[5]  = {pPrec->Y2m2, pPrec->Y2m1, pPrec->Y20, pPrec->Y21, pPrec->Y22};

  /* Precompute powers of e^{i m alpha} */
  COMPLEX16 cexp_i_alpha         = cexp(+I*alpha);
  COMPLEX16 cexp_2i_alpha        = cexp_i_alpha  * cexp_i_alpha;
  COMPLEX16 cexp_mi_alpha        = 1.0 / cexp_i_alpha;
  COMPLEX16 cexp_m2i_alpha       = cexp_mi_alpha * cexp_mi_alpha;
  COMPLEX16 cexp_im_alpha_l2[5]  = {cexp_m2i_alpha, cexp_mi_alpha, 1.0, cexp_i_alpha, cexp_2i_alpha};

  COMPLEX16 hp_sum               = 0;
  COMPLEX16 hc_sum               = 0;

  /* Loop over m' modes and perform the actual twisting up */
  for(int m=-2; m<=2; m++)
  {
    /* Transfer functions, see Eqs. 3.5-3.7 of arXiv:2004.06503 */
    COMPLEX16 A2m2emm    = cexp_im_alpha_l2[-m+2] * d2m2[m+2]  * Y2mA[m+2];        /*  = cexp(I*m*alpha) * d22[m+2]   * Y2mA[m+2] */
    COMPLEX16 A22emmstar = cexp_im_alpha_l2[m+2]  * d22[m+2]   * conj(Y2mA[m+2]);  /*  = cexp(-I*m*alpha) * d2m2[m+2] * conj(Y2mA[m+2])  */
    hp_sum +=    A2m2emm + A22emmstar;
    hc_sum += I*(A2m2emm - A22emmstar);
  }

  /* Note that \gamma = - \epsilon */
  COMPLEX16 eps_phase_hP = cexp(-2.0*I*epsilon) * hAS / 2.0;

  /* Return h_+ and h_x */
  *hp = eps_phase_hP * hp_sum;
  *hc = eps_phase_hP * hc_sum;

  return XLAL_SUCCESS;
}


#ifdef __cplusplus
}
#endif
