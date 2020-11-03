/*
* Copyright (C) 2019 Marta Colleoni, Cecilio García Quirós
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
//
//  Created by Marta on 13/02/2019.
//

#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>
#include <lal/Sequence.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/FrequencySeries.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/XLALError.h>

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _OPENMP
#define omp ignore
#endif

#define L_MAX 4

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1 //print debugging info
#endif

#ifndef PHENOMXHMMBAND
#define MBAND 0
#else
#define MBAND PHENOMXHMMBAND //use multibanding
#endif

#include "LALSimIMRPhenomXHM_internals.h"
#include "LALSimIMRPhenomXHM_internals.c"

#include "LALSimIMRPhenomXHM_structs.h"
#include "LALSimIMRPhenomXHM_qnm.h"
#include "LALSimIMRPhenomXHM_multiband.c"

//#include "LALSimIMRPhenomXHM.h" (non-precessing review version )

#include "LALSimIMRPhenomXPHM.c"

/* Note: This is declared in LALSimIMRPhenomX_internals.c and avoids namespace clash */
IMRPhenomX_UsefulPowers powers_of_lalpiHM;


//This is a wrapper function for adding higher modes to the ModeArray
static LALDict *IMRPhenomXHM_setup_mode_array(LALDict *lalParams);

/*
* Helper function to multiple hlm with Ylm.
* Adapted from LALSimIMREOBNRv2HMROMUtilities.c
*/
static int IMRPhenomXHMFDAddMode(
  COMPLEX16FrequencySeries *hptilde,  /**<[out] hp series*/
  COMPLEX16FrequencySeries *hctilde,  /**<[out] hc series */
  COMPLEX16FrequencySeries *hlmtilde, /**< hlm mode to add */
  REAL8 theta,                        /**< Inclination [rad] */
  REAL8 phi,                          /**< Azimuthal angle [rad]*/
  INT4 l,                             /**< l index of the lm mode */
  INT4 m,                             /**< m index of the lm mode */
  INT4 sym                            /**< Equatorial symmetry */
);


/* Return hptilde and hctilde from a sum of modes */
static int IMRPhenomXHM_MultiMode(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency domain h+ GW strain */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency domain hx GW strain */
  REAL8 m1_SI,                        /**< primary mass [kg] */
  REAL8 m2_SI,                        /**< secondary mass [kg] */
  REAL8 chi1z,                        /**< aligned spin of primary */
  REAL8 chi2z,                        /**< aligned spin of secondary */
  REAL8 f_min,                        /**< Starting GW frequency (Hz) */
  REAL8 f_max,                        /**< End frequency; 0 defaults to Mf = 0.3 */
  REAL8 deltaF,                       /**< Sampling frequency (Hz) */
  REAL8 distance,                     /**< distance of source (m) */
  REAL8 inclination,                  /**< inclination of source (rad) */
  REAL8 phiRef,                       /**< reference orbital phase (rad) */
  REAL8 fRef_In,                      /**< Reference frequency */
  LALDict *lalParams                  /**< LALDict struct */
);

/* Return hptilde and hctilde from a sum of modes */
static int IMRPhenomXHM_MultiMode2(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency domain h+ GW strain */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency domain hx GW strain */
  REAL8Sequence *freqs_In,            /**< min and max frequency [Hz] */
  IMRPhenomXWaveformStruct *pWF,      /**< waveform parameters */
  LALDict *lalParams                  /**< LALDict struct */
);


/* Definitions of functions */


/* Wrapper function for adding higher modes to the ModeArray */
static LALDict *IMRPhenomXHM_setup_mode_array(LALDict *lalParams)
{
  /* setup ModeArray */
  INT4 lalParams_In = 0;
  if (lalParams == NULL)
  {
    lalParams_In = 1;
    lalParams = XLALCreateDict();
  }
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams);

  /* If the mode array is empty, populate using a default choice of modes */
  if (ModeArray == NULL)
  {
    /* Default behaviour */
    XLAL_PRINT_INFO("Using default modes for IMRPhenomXHM.\n");
    ModeArray = XLALSimInspiralCreateModeArray();

    /* Only define +m modes as we get -m modes for free */
    /* IMRPhenomXHM has the following calibrated modes. 22 mode taken from IMRPhenomXAS */
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -1);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -3);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -4);

    XLALSimInspiralWaveformParamsInsertModeArray(lalParams, ModeArray);
  }
  else {XLAL_PRINT_INFO("Using custom modes for PhenomXHM.\n"); }

  XLALDestroyValue(ModeArray);
  if(lalParams_In == 1)
  {
    XLALDestroyDict(lalParams);
  }

  return lalParams;
}


/* Compute the frequency array and initialize htildelm to the corresponding length. */
int SetupWFArrays(
  REAL8Sequence **freqs,                 /**< [out] frequency grid [Hz */
  COMPLEX16FrequencySeries **htildelm,  /**< [out] Frequency domain hlm GW strain */
  REAL8Sequence *freqs_In,              /**< fmin, fmax [Hz] */
  IMRPhenomXWaveformStruct *pWF,        /**< Waveform structure with parameters */
  LIGOTimeGPS ligotimegps_zero          /**< = {0,0} */
)
{

  /* Inherit minimum and maximum frequencies to generate wavefom from input frequency grid */
  double f_min = freqs_In->data[0];
  double f_max = freqs_In->data[freqs_In->length - 1];

  /* Size of array */
  size_t npts = 0;

  /* Index shift between freqs and the frequency series */
  UINT4 offset = 0;

  /* If deltaF is non-zero then we need to generate a uniformly sampled frequency grid of spacing deltaF. Start at f = 0. */
  if(pWF->deltaF > 0)
  {
    /* Return the closest power of 2 */
    npts  = (size_t) (f_max / pWF->deltaF) + 1;

    /* Debug information */
    #if DEBUG == 1
    printf("npts     = %zu\n",npts);
    printf("fMin     = %.6f\n",f_min);
    printf("fMax     = %.6f\n",f_max);
    printf("dF       = %.6f\n",pWF->deltaF);
    printf("f_max / deltaF = %.6f\n", f_max / pWF->deltaF);
    #endif

    /* Coalescence time is fixed to t=0, shift by overall length in time. Model is calibrated such that it peaks approx 500M before the end of the waveform, add this time to the epoch. */
    XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / pWF->deltaF), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.",pWF->deltaF);
    /* Initialize the htilde frequency series */
    *htildelm = XLALCreateCOMPLEX16FrequencySeries("htildelm: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,npts);
    /* Check that frequency series generated okay */
    XLAL_CHECK(*htildelm,XLAL_ENOMEM,"Failed to allocate COMPLEX16FrequencySeries of length %zu for f_max = %f, deltaF = %g.\n",npts,f_max,pWF->deltaF);

    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t) (f_min / pWF->deltaF);
    size_t iStop  = (size_t) (f_max / pWF->deltaF) + 1;

    XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
    "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);

    #if DEBUG == 1
	  printf("f_min asked returned = %.16e %.16e \n",f_min, iStart*pWF->deltaF);
	  printf("f_max asked returned = %.16e %.16e \n",f_max, iStop*pWF->deltaF);
    #endif

    /* Allocate memory for frequency array and terminate if this fails */
    (*freqs) = XLALCreateREAL8Sequence(iStop - iStart);
    if (!(*freqs))
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
    /* Populate frequency array */
    for (UINT4 i = iStart; i < iStop; i++)
    {
      (*freqs)->data[i-iStart] = i * pWF->deltaF;
    }
    offset = iStart;
  }
  else
  {
    /* freqs is a frequency grid with non-uniform spacing, so we start at the lowest given frequency */
    npts      = freqs_In->length;
    *htildelm = XLALCreateCOMPLEX16FrequencySeries("htildelm: FD waveform, 22 mode", &ligotimegps_zero, f_min, pWF->deltaF, &lalStrainUnit, npts);
    XLAL_CHECK (*htildelm, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu from sequence.", npts);
    offset = 0;
    (*freqs)  = XLALCreateREAL8Sequence(freqs_In->length);

    /* Allocate memory for frequency array and terminate if this fails */
    if (!(*freqs))
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }

    /* Populate frequency array */
    for (UINT4 i = 0; i < freqs_In->length; i++)
    {
      (*freqs)->data[i] = freqs_In->data[i];
    }
  }//end loop of freqs

  memset((*htildelm)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htildelm)->sampleUnits), &((*htildelm)->sampleUnits), &lalSecondUnit);

  return offset;
}


/**
 * @addtogroup LALSimIMRPhenomX_c
 * @{
 *
 * @name Routines for IMRPhenomXHM
 * @{
 */



/********************************/
/*                              */
/*        SINGLE MODE           */
/*                              */
/********************************/

/* Functions to compute the strain of just one mode htildelm */

/** Returns the Fourier domain strain of just one mode: h_lm. Supports positive and negative m.
Notice than even when the specified frequencies are positives, the m>0 only lives in the negative frequencies regime.
With m>0 the mode h_lm is zero for positive frequencies and for the negative frequencies is equal to (-1)^l h*_l-m(-f).
In the contrary, h_l-m is zero for negative frequencies and only lives for positive frequencies.
This is a wrapper function that uses XLALSimIMRPhenomXASGenerateFD for the 22 mode and IMRPhenomXHMGenerateFDOneMode for the higher modes.
 */
 int XLALSimIMRPhenomXHMGenerateFDOneMode(
   COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
   REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
   REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
   REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
   REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
   UINT4 ell,                           /**< l index of the mode */
   INT4 emm,                            /**< m index of the mode */
   REAL8 distance,                      /**< Luminosity distance (m) */
   REAL8 f_min,                         /**< Starting GW frequency (Hz) */
   REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
   REAL8 deltaF,                        /**< Sampling frequency (Hz) */
   REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
   REAL8 fRef_In,                       /**< Reference frequency (Hz) */
   LALDict *lalParams                   /**< Extra params */
 )
 {

   /* If the 22 is required, call to PhenomX. */
   if(ell == 2 && abs(emm) == 2){
     XLALSimIMRPhenomXASGenerateFD(
       htildelm,
       m1_SI,
       m2_SI,
       chi1L,
       chi2L,
       distance,
       f_min,
       f_max,
       deltaF,
       phiRef,
       fRef_In,
       lalParams
     );
     if(emm>0){
       #if DEBUG == 1
       printf("\nTransforming to positive m by doing (-1)^l*Conjugate, frequencies must be negatives.\n");
       #endif
       for(UINT4 idx=0; idx<(*htildelm)->data->length; idx++){
         (*htildelm)->data->data[idx] = conj((*htildelm)->data->data[idx]);
       }
     }
     size_t iStart = (size_t) (f_min / deltaF);
     return iStart;
   }
   /***** Higher modes ****/
   else{

     UINT4 status;

     #if DEBUG == 1
     printf("\n**********************************************************************\n");
     printf("\n*                IMRPhenomXHMGenereateFDOneMode        %i%i            *\n", ell, abs(emm));
     printf("\n**********************************************************************\n");
     printf("fRef_In : %e\n",fRef_In);
     printf("m1_SI   : %e\n",m1_SI);
     printf("m2_SI   : %e\n",m2_SI);
     printf("chi1L   : %e\n",chi1L);
     printf("chi2L   : %e\n\n",chi2L);
     printf("Performing sanity checks...\n");
     #endif

     /* Sanity checks */
     if(*htildelm)       { XLAL_CHECK(NULL != htildelm, XLAL_EFAULT);                                   }
     if(fRef_In  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");  }
     if(deltaF   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaF must be positive.\n");                         }
     if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                             }
     if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                             }
     if(f_min    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "f_min must be positive.\n");                          }
     if(f_max    <  0.0) { XLAL_ERROR(XLAL_EDOM, "f_max must be non-negative.\n");                      }
     if(distance <  0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");    }

     /*
     	Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
     		- For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
     		- For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
     		- For mass ratios > 1000: throw a hard error that model is not valid.
     		- For spins > 0.99: throw a warning that we are extrapolating the model to extremal

     */
     REAL8 mass_ratio;
     if(m1_SI > m2_SI)
     {
   	  mass_ratio = m1_SI / m2_SI;
     }
     else
     {
   	  mass_ratio = m2_SI / m1_SI;
     }
     if(mass_ratio > 20.0  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
     if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); } // The 1e-12 is to avoid rounding errors
     if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

     /* Use an auxiliar laldict to not overwrite the input argument */
     LALDict *lalParams_aux;
     /* setup mode array */
     if (lalParams == NULL)
     {
         lalParams_aux = XLALCreateDict();
     }
     else{
         lalParams_aux = XLALDictDuplicate(lalParams);
     }
     lalParams_aux = IMRPhenomXHM_setup_mode_array(lalParams_aux);
     LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);

     /* first check if (l,m) mode is 'activated' in the ModeArray */
     /* if activated then generate the mode, else skip this mode. */
     if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm) != 1 )
     { /* skip mode */
       XLALPrintError("XLAL Error - %i%i mode is not included\n", ell, emm);
       XLAL_ERROR(XLAL_EDOM);
     } /* else: generate mode */


     /* If no reference frequency is given, we will set it to the starting gravitational wave frequency */
     REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;

     #if DEBUG == 1
     printf("\n\n **** Initializing waveform struct... **** \n\n");
     #endif

     /* Initialize the useful powers of LAL_PI */
     status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
     status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
     XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");


     /* Initialize IMRPhenomX Waveform struct and check that it generated successfully */
     IMRPhenomXWaveformStruct *pWF;
     pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
     status = IMRPhenomXSetWaveformVariables(pWF,m1_SI, m2_SI, chi1L, chi2L, deltaF, fRef, phiRef, f_min, f_max, distance, 0.0, lalParams_aux, DEBUG);
     XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

     #if DEBUG == 1
     printf("\n f_max_prime = %.16f, fMax = %.16f \n",pWF->f_max_prime, pWF->fMax);
     #endif


     /* Return the closest power of 2 */
     size_t npts = NextPow2(pWF->f_max_prime / deltaF) + 1;
     /* Frequencies will be set using only the lower and upper bounds that we passed */
     size_t iStart = (size_t) (f_min / deltaF);
     size_t iStop  = (size_t) (pWF->f_max_prime / deltaF);
     XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
     "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);

     /* Create a REAL8 frequency series only passing the boundaries (fMin, fMax).  */
     REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
     freqs->data[0] = pWF->fMin;
     freqs->data[1] = pWF->f_max_prime;


     #if DEBUG == 1
     printf("\n\nfstart, fend = %.16f %.16f\n\n", freqs->data[0], freqs->data[1]);
     printf("\n\n **** Calling IMRPhenomXHMGenerateFDOneMode... **** \n\n");
     printf("\n f_max_prime = %.16f, fMax = %.16f \n", pWF->f_max_prime, pWF->fMax);
     #endif

     /*** Call core single mode waveform generator ***/
     status = IMRPhenomXHMGenerateFDOneMode(htildelm, freqs, pWF, ell, abs(emm), lalParams_aux);
     XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHMGenerateFDOneMode failed to generate IMRPhenomXHM waveform.");

     #if DEBUG == 1
     printf("\n\n **** Call to IMRPhenomXHMGenerateFD complete. **** \n\n");
     printf("\n f_max_prime = %.16f, fMax = %.16f \n",pWF->f_max_prime, pWF->fMax);
     #endif


     if(emm>0){
       #if DEBUG == 1
       printf("\nTransforming to positive m by doing (-1)^l*Conjugate, frequencies must be negatives.\n");
       #endif
       INT4 minus1l = 1;
       if(ell%2!=0){
         #if DEBUG == 1
         printf("\nl odd\n");
         #endif
         minus1l = -1;
       }
       for(UINT4 idx=0; idx<(*htildelm)->data->length; idx++){
         (*htildelm)->data->data[idx] = minus1l*conj((*htildelm)->data->data[idx]);
       }
     }

     /* Resize htildelm if needed */
     REAL8 lastfreq;
     if (pWF->f_max_prime < pWF->fMax)
     {
       /* The user has requested a higher f_max than Mf = fCut.
       Resize the frequency series to fill with zeros beyond the cutoff frequency. */
       lastfreq = pWF->fMax;
     }
     else
     {
       lastfreq = pWF->f_max_prime;
     }
     size_t n = (*htildelm)->data->length;
     // We want to have the length be a power of 2 + 1
     size_t n_full = NextPow2(lastfreq / pWF->deltaF) + 1;

     /* Resize the COMPLEX16 frequency series */
     *htildelm = XLALResizeCOMPLEX16FrequencySeries(*htildelm, 0, n_full);
     XLAL_CHECK (*htildelm, XLAL_ENOMEM, "Failed to resize waveform COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->f_max_prime, n_full, pWF->fMax );


     /* Free allocated memory */
     LALFree(pWF);
     XLALDestroyREAL8Sequence(freqs);
     XLALDestroyValue(ModeArray);
     XLALDestroyDict(lalParams_aux);

     #if DEBUG == 1
     printf("\n Leaving XLALSimIMRPhenomXHMGenerateFDOneMode \n");
     #endif

     return iStart;
   }//Higher modes
 }
/** @} **
* @} **/

/* Core function to generate the waveform of one mode: htildelm. */
int IMRPhenomXHMGenerateFDOneMode(
  COMPLEX16FrequencySeries **htildelm,  /**< [out] hlm for one mode **/
  REAL8Sequence *freqs_In,              /**< fmin, fmax [Hz] **/
  IMRPhenomXWaveformStruct *pWF,        /**< structure of the 22 mode **/
  UINT4 ell,                            /**< first index of the mode **/
  UINT4 emm,                            /**< second index of the mode **/
  LALDict *lalParams                    /**< extra params **/
)
{

  #if DEBUG == 1
  printf("\n\n ***** IMRPhenomXHMGenerateFDOneMode **** \n\n");
  #endif

  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  UINT4 initial_status = XLAL_SUCCESS;

  REAL8Sequence *freqs;
  UINT4 offset = SetupWFArrays(&freqs, htildelm, freqs_In, pWF, ligotimegps_zero);

  #if DEBUG == 1
  printf("\nIMRPhenomXHMGenerateFDOneMode: Length@freqs, offset %i %u",freqs->length,offset);
  printf("\nIMRPhenomXHMGenerateFDOneMode: fstart, fend = %.16f %.16f\n\n", freqs->data[0], freqs->data[freqs->length-1]);
  #endif

  /* Check if LAL dictionary exists. If not, create a LAL dictionary. */
  INT4 lalParams_In = 0;
  if(lalParams == NULL)
  {
    lalParams_In = 1;
    lalParams = XLALCreateDict();
  }

  // allocate qnm struct
  QNMFits *qnms = (QNMFits *) XLALMalloc(sizeof(QNMFits));
  IMRPhenomXHM_Initialize_QNMs(qnms);

  // Populate pWFHM
  IMRPhenomXHMWaveformStruct *pWFHM = (IMRPhenomXHMWaveformStruct *) XLALMalloc(sizeof(IMRPhenomXHMWaveformStruct));
  IMRPhenomXHM_SetHMWaveformVariables(ell, emm, pWFHM, pWF, qnms, lalParams);
  LALFree(qnms);

  /* Fill only the modes that are not zero. Odd modes for equal mass, equal spins are zero. */
  /* Since hlmtilde is initialized with zeros, for zero modes we just go to the return line. */
  if(pWFHM->Ampzero==0){

    /* Allocate coefficients of 22 mode */
    IMRPhenomXAmpCoefficients *pAmp22=(IMRPhenomXAmpCoefficients *) XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
    IMRPhenomXPhaseCoefficients *pPhase22=(IMRPhenomXPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
    IMRPhenomXGetPhaseCoefficients(pWF, pPhase22);

    /* Allocate and initialize the PhenomXHM lm amplitude and phae coefficients struct */
    IMRPhenomXHMAmpCoefficients *pAmp = (IMRPhenomXHMAmpCoefficients*) XLALMalloc(sizeof(IMRPhenomXHMAmpCoefficients));
    IMRPhenomXHMPhaseCoefficients *pPhase = (IMRPhenomXHMPhaseCoefficients*) XLALMalloc(sizeof(IMRPhenomXHMPhaseCoefficients));

    /* Allocate and initialize the PhenomXHM lm phase and amp coefficients struct */
    IMRPhenomXHM_FillAmpFitsArray(pAmp);
    IMRPhenomXHM_FillPhaseFitsArray(pPhase);

    /* Get coefficients for Amplitude and phase */
    if (pWFHM->MixingOn == 1) {
      // For mode with mixing we need the spheroidal coeffs of the 32 phase and the 22 amplitude coeffs.
      GetSpheroidalCoefficients(pPhase, pPhase22, pWFHM, pWF);
      IMRPhenomXGetAmplitudeCoefficients(pWF, pAmp22);
    }
    IMRPhenomXHM_GetAmplitudeCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF);
    IMRPhenomXHM_GetPhaseCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF,lalParams);


    IMRPhenomX_UsefulPowers powers_of_Mf;
    REAL8 Msec = pWF->M_sec;    // Variable to transform Hz to Mf
    REAL8 Amp0 = pWFHM->Amp0;   // Transform amplitude from NR to physical units
    REAL8 amp, phi;

    /* Multiply by (-1)^l to get the true h_l-m(f) */
    if(ell%2 != 0){
      Amp0 = -Amp0;
    }

    #if DEBUG == 1
    //Initialize file to print h_l-m with freqs, amplitude and phase in "Physical" units
    FILE *file;
    char fileSpec[40];
    sprintf(fileSpec, "simulation%i_SM.dat", pWFHM->modeTag);
    printf("\nOutput file: %s\r\n",fileSpec);
    file = fopen(fileSpec,"w");
    fprintf(file,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i Mtot = %.16f distance = %.16f\n", pWF->q, pWF->chi1L, pWF->chi2L, 22, pWF->Mtot, pWF->distance/LAL_PC_SI/pow(10.,6));
    fprintf(file,"# Frequency (Hz)    Amplitude    Phase\n");
    #endif


    /* Loop over frequencies to generate waveform */
    /* Modes with mixing */
    if(pWFHM->MixingOn==1){

      REAL8 Mf;
      for (UINT4 idx = 0; idx < freqs->length; idx++)
      {
        Mf    = Msec * freqs->data[idx];
        initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
        if(initial_status != XLAL_SUCCESS)
        {
          XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
        }
        else
        {
          amp = IMRPhenomXHM_Amplitude_ModeMixing(Mf, &powers_of_Mf, pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF);
          phi = IMRPhenomXHM_Phase_ModeMixing(Mf, &powers_of_Mf, pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF);
          /* Reconstruct waveform: h_l-m(f) = A(f) * Exp[I phi(f)] */
          ((*htildelm)->data->data)[idx+offset] = Amp0 * amp * cexp(I * phi);

          #if DEBUG == 1
          fprintf(file, "%.16f  %.16e %.16f\n",  freqs->data[idx], Amp0*amp, phi);
          #endif
        }
      }
    }  /* Modes without mixing */
    else{
      for (UINT4 idx = 0; idx < freqs->length; idx++)
      {
        REAL8 Mf    = Msec * freqs->data[idx];
        initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
        if(initial_status != XLAL_SUCCESS)
        {
          XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
        }
        else
        {
          amp = IMRPhenomXHM_Amplitude_noModeMixing(Mf, &powers_of_Mf, pAmp, pWFHM);
          phi = IMRPhenomXHM_Phase_noModeMixing(Mf, &powers_of_Mf, pPhase, pWFHM, pWF);
          /* Reconstruct waveform: h_l-m(f) = A(f) * Exp[I phi(f)] */
          ((*htildelm)->data->data)[idx+offset] = Amp0 * amp * cexp(I * phi);
          #if DEBUG == 1
          fprintf(file, "%.16f  %.16e %.16f\n",  freqs->data[idx], Amp0*amp, phi);
          #endif
        }
      }
    }
    #if DEBUG == 1
    fclose(file);
    #endif


    // Free allocated memory
    LALFree(pAmp);
    LALFree(pPhase);
    LALFree(pAmp22);
    LALFree(pPhase22);
  } // End of non-vanishing modes

  // Free allocated memory
  LALFree(pWFHM);
  XLALDestroyREAL8Sequence(freqs);
  if(lalParams_In == 1)
  {
    XLALDestroyDict(lalParams);
  }

  #if DEBUG == 1
  printf("\nIMRPhenomXHMGenerateFDOneMode: Leaving... \n");
  #endif


  return initial_status;

}


/** @addtogroup LALSimIMRPhenomX_c
* @{
* @name Routines for IMRPhenomXHM
* @{ **/

/** Get the model evaluated in a prebuilt frequency array.
    It does not use Multibanding.
 */
int XLALSimIMRPhenomXHMFrequencySequenceOneMode(
  COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
  REAL8Sequence *freqs,                /**< frequency array to evaluate model (positives) */
  REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
  UINT4 ell,                           /**< l index of the mode */
  INT4 emm,                            /**< m index of the mode */
  REAL8 distance,                      /**< Luminosity distance (m) */
  REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
  REAL8 fRef_In,                       /**< Reference frequency (Hz) */
  LALDict *lalParams
)
{

  if(ell == 2 && abs(emm) == 2){
    INT4 status = XLALSimIMRPhenomXASFrequencySequence(
      htildelm,
      freqs,
      m1_SI,
      m2_SI,
      chi1L,
      chi2L,
      distance,
      phiRef,
      fRef_In,
      lalParams
    );
    XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "XLALSimIMRPhenomXHMFrequencySequenceOneMode failed to generate IMRPhenomXHM waveform.");
    /* Transfor to positive m mode if needed. Do (-1)^l*Conjugate[htildelm]. The frequencies must be negative. */
    if(emm>0){
      for(UINT4 idx=0; idx<(*htildelm)->data->length; idx++){
        (*htildelm)->data->data[idx] = conj((*htildelm)->data->data[idx]);
      }
    }
    return status;
  }
  else{ //Higher modes

    /* Variable to check correct calls to functions. */
    INT4 status = 0;

    /* Sanity checks */
    if(*htildelm)       { XLAL_CHECK(NULL != htildelm, XLAL_EFAULT);                                   }
    if(fRef_In  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");  }
    if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                             }
    if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                             }
    if(distance <  0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");    }

    /*
      Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
        - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
        - For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
        - For mass ratios > 1000: throw a hard error that model is not valid.
        - For spins > 0.99: throw a warning that we are extrapolating the model to extremal

    */
    REAL8 mass_ratio;
    if(m1_SI > m2_SI)
    {
      mass_ratio = m1_SI / m2_SI;
    }
    else
    {
      mass_ratio = m2_SI / m1_SI;
    }
    if(mass_ratio > 20.0  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
    if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); } // The 1e-12 is to avoid rounding errors
    if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

    /* Use an auxiliar laldict to not overwrite the input argument */
    LALDict *lalParams_aux;
    /* setup mode array */
    if (lalParams == NULL)
    {
        lalParams_aux = XLALCreateDict();
    }
    else{
        lalParams_aux = XLALDictDuplicate(lalParams);
    }
    lalParams_aux = IMRPhenomXHM_setup_mode_array(lalParams_aux);
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);

    /* first check if (l,m) mode is 'activated' in the ModeArray */
    /* if activated then generate the mode, else skip this mode. */
    if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm) != 1 )
    { /* skip mode */
      XLALPrintError("XLAL Error - %i%i mode is not included\n", ell, emm);
      XLAL_ERROR(XLAL_EDOM);
    } /* else: generate mode */


    /* If fRef is not provided (i.e. set to 0), then take fRef to be the starting GW Frequency. */
    REAL8 fRef = (fRef_In == 0.0) ? freqs->data[0] : fRef_In;


    /* Initialize the useful powers of LAL_PI */
    status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Failed to initialize useful powers of LAL_PI.");
    status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Failed to initialize useful powers of LAL_PI.");

    /* Get minimum and maximum frequencies. */
    REAL8 f_min_In  = freqs->data[0];
    REAL8 f_max_In  = freqs->data[freqs->length - 1];

    /*
      Initialize IMRPhenomX waveform struct and perform sanity check.
      Passing deltaF = 0 tells us that freqs contains a frequency grid with non-uniform spacing.
      The function waveform start at lowest given frequency.
    */
    IMRPhenomXWaveformStruct *pWF;
    pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
    status = IMRPhenomXSetWaveformVariables(pWF,m1_SI, m2_SI, chi1L, chi2L, 0.0, fRef, phiRef, f_min_In, f_max_In, distance, 0.0, lalParams_aux, PHENOMXDEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error:  failed.\n");

    /* Now call the IMRPhenomXHM one mode waveform generator */
    status = IMRPhenomXHMGenerateFDOneMode(
      htildelm,
      freqs,
      pWF,
      ell,
      abs(emm),
      lalParams_aux
    );
    XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "XLALSimIMRPhenomXHMFrequencySequenceOneMode failed to generate IMRPhenomXHM waveform.");

    /* Transfor to positive m mode if needed. Do (-1)^l*Conjugate[htildelm]. The frequencies must be negative. */
    if(emm>0){
      INT4 minus1l = 1;
      if(ell%2!=0){
        minus1l = -1;
      }
      for(UINT4 idx=0; idx<(*htildelm)->data->length; idx++){
        (*htildelm)->data->data[idx] = minus1l*conj((*htildelm)->data->data[idx]);
      }
    }

    /* Free memory */
    LALFree(pWF);
    XLALDestroyValue(ModeArray);
    XLALDestroyDict(lalParams_aux);


    return XLAL_SUCCESS;

  }//Higher modes
}

/*********************************************/
/*                                           */
/*          MULTIMODE WAVEFORM               */
/*                                           */
/*********************************************/

/** Returns the hptilde and hctilde of the multimode waveform for positive frequencies.

XLALSimIMRPhenomXHM calls the function for a single mode that can be XLALSimIMRPhenomXHMGenerateFDOneMode or XLALSimIMRPhenomXHMMultiBandOneMode,
depending on if the Multibanding is active or not.

By default XLALSimIMRPhenomXHM is only used when the Multibanding is activated, since each mode has a different coarse frequency array and we can not recycle the array.

This is just a wrapper of the function that actually carry out the calculations: IMRPhenomXHM_MultiMode2.

*/

/* Return hptilde, hctilde */
int XLALSimIMRPhenomXHM(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx */
  REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
  REAL8 chi1L,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2L,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 f_min,                        /**< Starting GW frequency (Hz) */
  REAL8 f_max,                        /**< End frequency; 0 defaults to Mf = 0.3 */
  REAL8 deltaF,                       /**< Sampling frequency (Hz) */
  REAL8 distance,                     /**< distance of source (m) */
  REAL8 inclination,                  /**< inclination of source (rad) */
  REAL8 phiRef,                       /**< reference orbital phase (rad) */
  REAL8 fRef_In,                      /**< Reference frequency */
  LALDict *lalParams                  /**<linked list containing the extra parameters */
)
{
  /* define and init return code for this function */
  int retcode;

  /* Sanity checks on input parameters: check pointers, etc.
  More checks are done inside XLALSimIMRPhenomXHMGenerateFDOneMode/XLALSimIMRPhenomXHMMultiBandOneMode */
  XLAL_CHECK(NULL != hptilde, XLAL_EFAULT);
  XLAL_CHECK(NULL != hctilde, XLAL_EFAULT);
  XLAL_CHECK(*hptilde == NULL, XLAL_EFAULT);
  XLAL_CHECK(*hctilde == NULL, XLAL_EFAULT);
  XLAL_CHECK(distance > 0, XLAL_EDOM, "distance must be positive.\n");
  /*
  	Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
  		- For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
  		- For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
  		- For mass ratios > 1000: throw a hard error that model is not valid.
  		- For spins > 0.99: throw a warning that we are extrapolating the model to extremal

  */
  REAL8 mass_ratio;
  if(m1_SI > m2_SI)
  {
	  mass_ratio = m1_SI / m2_SI;
  }
  else
  {
	  mass_ratio = m2_SI / m1_SI;
  }
  if(mass_ratio > 20.0  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
  if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

  /* Check that the modes chosen are available for the model */
  XLAL_CHECK(check_input_mode_array(lalParams) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");


  /* Evaluate the model */
  retcode = IMRPhenomXHM_MultiMode(
    hptilde,
    hctilde,
    m1_SI,
    m2_SI,
    chi1L,
    chi2L,
    f_min,
    f_max,
    deltaF,
    distance,
    inclination,
    phiRef,
    fRef_In,
    lalParams
  );
  XLAL_CHECK(retcode == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHM_MultiMode failed in XLALSimIMRPhenomXHM.");

  #if DEBUG == 1
  printf("\n******Leaving XLALSimIMRPhenomXHM*****\n");
  #endif

  return retcode;
}

/** Returns the hptilde and hctilde of the multimode waveform for positive frequencies.

XLALSimIMRPhenomXHM2 builds each mode explicitly in the loop over modes, recycling some common quantities between modes like
the frequency array, the powers of frequencies, structures, etc so it has less overhead than XLALSimIMRPhenomXHM.

By default XLALSimIMRPhenomXHM2 is only used for the version without Multibanding.

 This is just a wrapper of the function that actually carry out the calculations: IMRPhenomXHM_MultiMode2.
 */
int XLALSimIMRPhenomXHM2(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency domain h+ GW strain */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency domain hx GW strain */
  REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
  REAL8 f_min,                         /**< Starting GW frequency (Hz) */
  REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
  REAL8 deltaF,                        /**< Sampling frequency (Hz) */
  REAL8 distance,                      /**< Luminosity distance (m) */
  REAL8 inclination,                   /**< Inclination of the source */
  REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
  REAL8 fRef_In,                       /**< Reference frequency (Hz) */
  LALDict *lalParams                   /**< LAL Dictionary */
)
{
  UINT4 debug = DEBUG;

  UINT4 status;

  #if DEBUG == 1
  printf("\n*** IMRPhenomXHM2 ***\n");
  printf("fRef_In : %e\n",fRef_In);
  printf("m1_SI   : %.16e\n",m1_SI);
  printf("m2_SI   : %.16e\n",m2_SI);
  printf("chi1L   : %.16e\n",chi1L);
  printf("chi2L   : %.16e\n",chi2L);
  printf("f_min   : %.16e  \n", f_min);
  printf("Performing sanity checks...\n");
  #endif

  /* Perform initial sanity checks */
  if(*hptilde)       { XLAL_CHECK(NULL != hptilde, XLAL_EFAULT);                                   }
  if(*hctilde)       { XLAL_CHECK(NULL != hctilde, XLAL_EFAULT);                                   }
  if(fRef_In  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");  }
  if(deltaF   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaF must be positive.\n");                         }
  if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                             }
  if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                             }
  if(f_min    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "f_min must be positive.\n");                          }
  if(f_max    <  0.0) { XLAL_ERROR(XLAL_EDOM, "f_max must be non-negative.\n");                      }
  if(distance <  0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");    }
  /*
  	Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
  		- For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
  		- For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
  		- For mass ratios > 1000: throw a hard error that model is not valid.
  		- For spins > 0.99: throw a warning that we are extrapolating the model to extremal

  */
  REAL8 mass_ratio;
  if(m1_SI > m2_SI)
  {
	  mass_ratio = m1_SI / m2_SI;
  }
  else
  {
	  mass_ratio = m2_SI / m1_SI;
  }
  if(mass_ratio > 20.0  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
  if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

  /* Check that the modes chosen are available for the model */
  XLAL_CHECK(check_input_mode_array(lalParams) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");

  /* If no reference frequency is given, set it to the starting gravitational wave frequency */
  REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;


  #if DEBUG == 1
  printf("\nInitializing waveform struct...\n");
  #endif

  /* Initialize the useful powers of LAL_PI */
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

  /* Initialize IMR PhenomX Waveform struct and check that it initialized correctly */
  IMRPhenomXWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, deltaF, fRef, phiRef, f_min, f_max, distance, inclination, lalParams, debug);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error:  failed.\n");


  /* Check that the frequency array will be consistent: fmin < fmax_prime */
  /* Return the closest power of 2 */
  size_t npts = NextPow2(pWF->f_max_prime / deltaF) + 1;
  /* Frequencies will be set using only the lower and upper bounds that we passed */
  size_t iStart = (size_t) (f_min / deltaF);
  size_t iStop  = (size_t) (pWF->f_max_prime / deltaF);
  XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
  "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);

  /*  Create a REAL8 frequency series.
  Use fLow, fHigh, deltaF to compute frequency sequence. Only pass the boundaries (fMin, f_max_prime).   */
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = pWF->fMin;
  freqs->data[1] = pWF->f_max_prime;


  /* We now call the core IMRPhenomXHM_Multimode2 waveform generator. */
  status = IMRPhenomXHM_MultiMode2(hptilde, hctilde, freqs, pWF, lalParams);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHM_MultiMode2 failed to generate IMRPhenomXHM waveform.");

  #if DEBUG == 1
  printf("\n\n **** Call to IMRPhenomXHM_MultiMode2 complete. **** \n\n");
  #endif


  /* Resize hptilde, hctilde */
  REAL8 lastfreq;
  if (pWF->f_max_prime < pWF->fMax)
  {
    /* The user has requested a higher f_max than Mf = fCut.
    Resize the frequency series to fill with zeros beyond the cutoff frequency. */
    lastfreq = pWF->fMax;
  }
  else{  // We have to look for a power of 2 anyway.
    lastfreq = pWF->f_max_prime;
  }
  // We want to have the length be a power of 2 + 1
  size_t n_full = NextPow2(lastfreq / deltaF) + 1;
  size_t n = (*hptilde)->data->length;

  /* Resize the COMPLEX16 frequency series */
  *hptilde = XLALResizeCOMPLEX16FrequencySeries(*hptilde, 0, n_full);
  XLAL_CHECK (*hptilde, XLAL_ENOMEM, "Failed to resize h_+ COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );

  /* Resize the COMPLEX16 frequency series */
  *hctilde = XLALResizeCOMPLEX16FrequencySeries(*hctilde, 0, n_full);
  XLAL_CHECK (*hctilde, XLAL_ENOMEM, "Failed to resize h_x COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );


  /* Free memory */
  LALFree(pWF);
  XLALDestroyREAL8Sequence(freqs);

  return XLAL_SUCCESS;
}


/**
 * Returns hptilde and hctilde as a complex frequency series with entries exactly at the frequencies specified in
 * the REAL8Sequence *freqs (which can be unequally spaced). No zeros are added. Assumes positive frequencies.
 */

 int XLALSimIMRPhenomXHMFrequencySequence(
    COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+  */
    COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx  */
    REAL8Sequence *freqs,               /**< Input Frequency series [Hz]         */
    REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
    REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
    REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 distance,                     /**< Distance of source (m) */
    REAL8 inclination,                  /**< inclination of source (rad) */
    REAL8 phiRef,                       /**< Orbital phase (rad) at reference frequency */
    REAL8 fRef_In,                         /**< Reference frequency (Hz) */
    LALDict *lalParams                  /**< LAL Dictionary struct */
 )
 {
   /* Variable to check correct calls to functions. */
   INT4 status;

   /* Perform initial sanity checks */
   XLAL_CHECK(NULL != hptilde, XLAL_EFAULT, "Error: hptilde already defined.                        \n");
   XLAL_CHECK(NULL != hctilde, XLAL_EFAULT, "Error: hctilde already defined.                        \n");
   XLAL_CHECK(NULL != freqs,   XLAL_EFAULT, "Error: Input freq array must be defined.               \n");
   XLAL_CHECK(fRef_In  >= 0, XLAL_EFUNC,    "Error: fRef must be positive and greater than 0.       \n");
   XLAL_CHECK(m1_SI    >  0, XLAL_EFUNC,    "Error: m1 must be positive and greater than 0.         \n");
   XLAL_CHECK(m2_SI    >  0, XLAL_EFUNC,    "Error: m2 must be positive and greater than 0.         \n");
   XLAL_CHECK(distance >  0, XLAL_EFUNC,    "Error: Distance must be positive and greater than 0.   \n");

   /*
   Perform a basic sanity check on the region of the parameter space in which model is evaluated.
   Behaviour is as follows, consistent with the choices for IMRPhenomXAS/IMRPhenomXHM
     - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
     - For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
     - For mass ratios > 1000: throw a hard error that model is not valid.
     - For spins > 0.99: throw a warning that we are extrapolating the model to extremal

   */
   REAL8 mass_ratio;
   if(m1_SI > m2_SI)
   {
     mass_ratio = m1_SI / m2_SI;
   }
   else
   {
     mass_ratio = m2_SI / m1_SI;
   }
   if(mass_ratio > 20.0  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
   if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); } // The 1e-12 is to avoid rounding errors
   if(fabs(chi1z) > 0.99 || fabs(chi2z) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

   /* If no reference frequency is given, set it to the starting gravitational wave frequency */
   REAL8 fRef = (fRef_In == 0.0) ? freqs->data[0] : fRef_In; //It is giving valgrind error, but it is not needed. f_ref = f_min in WaveformCache.c and SimInspiral.c.

   /* Get minimum and maximum frequencies. */
   REAL8 f_min_In  = freqs->data[0];
   REAL8 f_max_In  = freqs->data[freqs->length - 1];

   /*
     Passing deltaF = 0 implies that freqs is a frequency grid with non-uniform spacing.
     The function waveform then start at lowest given frequency.
     The Multibanding has to be switched off since it is built only for equally spaced frequency grid.
   */
   /* Use an auxiliar laldict to not overwrite the input argument */
   LALDict *lalParams_aux;
   /* setup mode array */
   if (lalParams == NULL)
   {
       lalParams_aux = XLALCreateDict();
   }
   else{
       lalParams_aux = XLALDictDuplicate(lalParams);
   }
   if(XLALSimInspiralWaveformParamsLookupPhenomXHMThresholdMband(lalParams_aux)!=0)
   {
     XLAL_PRINT_WARNING("Warning: Function is aimed for non-uniform frequency grid, switching off Multibanding.");
     XLALSimInspiralWaveformParamsInsertPhenomXHMThresholdMband(lalParams_aux, 0);
   }

   /* Initialize IMRPhenomX waveform struct and perform sanity check. */
   IMRPhenomXWaveformStruct *pWF;
   pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
   status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, 0.0, fRef, phiRef, f_min_In, f_max_In, distance, inclination, lalParams_aux, DEBUG);
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

   /* Call the core IMRPhenomXHM waveform generator without multibanding. */
   status = IMRPhenomXHM_MultiMode2(hptilde, hctilde, freqs, pWF, lalParams_aux);
   XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXPHM_hplushcross failed to generate IMRPhenomXPHM waveform.");

   /* Free memory */
   LALFree(pWF);
   XLALDestroyDict(lalParams_aux);

   return status;
 }

/** @}
* @} **/


/* Core function of XLALSimIMRPhenomXHM, returns hptilde, hctilde corresponding to a sum of modes.
The default modes are 22, 21, 33, 32 and 44. It returns also the contribution of the corresponding negatives modes. */
static int IMRPhenomXHM_MultiMode(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency domain h+ GW strain */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency domain hx GW strain */
  REAL8 m1_SI,                        /**< primary mass [kg] */
  REAL8 m2_SI,                        /**< secondary mass [kg] */
  REAL8 chi1z,                        /**< aligned spin of primary */
  REAL8 chi2z,                        /**< aligned spin of secondary */
  REAL8 f_min,                        /**< Starting GW frequency (Hz) */
  REAL8 f_max,                        /**< End frequency; 0 defaults to Mf = 0.3 */
  REAL8 deltaF,                       /**< Sampling frequency (Hz) */
  REAL8 distance,                     /**< distance of source (m) */
  REAL8 inclination,                  /**< inclination of source (rad) */
  REAL8 phiRef,                       /**< reference orbital phase (rad) */
  REAL8 fRef_In,                      /**< Reference frequency */
  LALDict *lalParams                  /**< LALDict struct */
)
{
  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  INT4 sym; /* sym will decide whether to add the -m mode (when equatorial symmetry is present) */

  /* Use an auxiliar laldict to not overwrite the input argument */
  LALDict *lalParams_aux;
  /* setup mode array */
  if (lalParams == NULL)
  {
      lalParams_aux = XLALCreateDict();
  }
  else{
      lalParams_aux = XLALDictDuplicate(lalParams);
  }
  lalParams_aux = IMRPhenomXHM_setup_mode_array(lalParams_aux);
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);

  /* At this point ModeArray should contain the list of modes
  and therefore if NULL then something is wrong and abort. */
  if (ModeArray == NULL)
  {
    XLAL_ERROR(XLAL_EDOM, "ModeArray is NULL when it shouldn't be. Aborting.\n");
  }

  INT4 status = 0;

  INT4 init = 0; //Variable to control initialization of hptilde and hctilde

  /* When calling only the 32 mode, we need to call also the 22 for the mixing, but not sum it for hp, hc */
  INT4 add22 = 1;


  COMPLEX16FrequencySeries *htilde22 = NULL;
  INT4 posMode, negMode;
  /* Take input/default value for the threshold of the Multibanding. If = 0 then do not use Multibanding. */
  REAL8 resTest  = XLALSimInspiralWaveformParamsLookupPhenomXHMThresholdMband(lalParams_aux);

  /***** Loop over modes ******/
  for (UINT4 ell = 2; ell <= L_MAX; ell++)
  {
    for (UINT4 emm = 1; emm <= ell; emm++)
    { /* Loop over only positive m is intentional.
      The single mode function returns the negative mode h_l-m, and the positive is added automatically in IMRPhenomXHMFDAddMode. */
      /* First check if (l,m) mode is 'activated' in the ModeArray */
      /* if activated then generate the mode, else skip this mode. */
      posMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm);
      negMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, -emm);
      if ( posMode != 1 && negMode != 1)
      { /* skip mode */
        continue;
      } /* else: generate mode */
      #if DEBUG == 1
      printf("\n Mode %i%i\n",ell, emm);
      #endif

      // Variable to store the strain of only one (negative) mode: h_l-m
      COMPLEX16FrequencySeries *htildelm = NULL;


      if (resTest == 0){  // No multibanding
        XLALSimIMRPhenomXHMGenerateFDOneMode(&htildelm, m1_SI, m2_SI, chi1z, chi2z, ell, -emm, distance, f_min, f_max, deltaF, phiRef, fRef_In, lalParams_aux);
      }
      else{               // With multibanding
        if(ell==3 && emm==2){  // mode with mixing
          XLALSimIMRPhenomXHMMultiBandOneModeMixing(&htildelm, htilde22, m1_SI, m2_SI, chi1z, chi2z, ell, -emm, distance, f_min, f_max, deltaF, phiRef, fRef_In, lalParams_aux);
        }
        else{                  // modes without mixing
          XLALSimIMRPhenomXHMMultiBandOneMode(&htildelm, m1_SI, m2_SI, chi1z, chi2z, ell, -emm, distance, f_min, f_max, deltaF, phiRef, fRef_In, lalParams_aux);
        }
        // If the 22 mode is active we will recycle for the mixing of the 32, we save it in another variable: htilde22.
        if(ell==2 && emm==2 && htildelm){
          htilde22 = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &(ligotimegps_zero), 0.0, deltaF, &lalStrainUnit, htildelm->data->length);
          for(UINT4 idx = 0; idx < htildelm->data->length; idx++){
            htilde22->data->data[idx] = htildelm->data->data[idx];
          }
        }
      }

      /**** For debugging ****/
      #if DEBUG == 1
      //Save mode l-m in file
      FILE *file;
      char fileSpec[40];
      sprintf(fileSpec, "simulation%i%i_MM.dat", ell,emm);
      printf("\nOutput file: %s\r\n",fileSpec);
      file = fopen(fileSpec,"w");
      REAL8 m2_SI_a, m1_SI_a, chi1z_a, chi2z_a;
      if(m1_SI < m2_SI){
        m2_SI_a = m1_SI;
        m1_SI_a = m2_SI;
        chi1z_a = chi2z;
        chi2z_a = chi1z;
      }
      else{
        m1_SI_a = m1_SI;
        m2_SI_a = m2_SI;
        chi1z_a = chi1z;
        chi2z_a = chi2z;
      }
      fprintf(file,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i Mtot = %.16f distance = %.16f\n", m1_SI_a/m2_SI_a, chi1z_a, chi2z_a, ell, emm, (m1_SI+m2_SI)/ LAL_MSUN_SI, distance/LAL_PC_SI/pow(10.,6));
      fprintf(file,"# Frequency (Hz)    Real     Imaginary\n");
      COMPLEX16 data;
      for(UINT4 idx = 0; idx < ((htildelm)->data->length); idx++)
      {
        data = ((htildelm)->data->data)[idx];
        fprintf(file, "%.16f  %.16e %.16e\n",  idx*deltaF, creal(data), cimag(data));
      }
      fclose(file);
      #endif
      /**** End debugging ****/


      if (!(htildelm)){ XLAL_ERROR(XLAL_EFUNC); return XLAL_FAILURE;}

      /* We test for hypothetical m=0 modes */
      if (emm == 0)
      {
        sym = 0;
      }
      else
      {
        sym = 1;
      }
      /* Allocate hptilde and hctilde if they were not yet. */
      /* Taken from IMRPhenomHM */
      if( init == 0){
        /* Coalescence time is fixed to t=0, shift by overall length in time. Model is calibrated such that it peaks approx 500M before the end of the waveform, add this time to the epoch. */
        XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / deltaF), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.", deltaF);
        size_t n = (htildelm)->data->length;
        *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &(ligotimegps_zero), 0.0, deltaF, &lalStrainUnit, n);
        if (!(hptilde)){  XLAL_ERROR(XLAL_EFUNC);}
        memset((*hptilde)->data->data, 0, n * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
        *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform",  &(ligotimegps_zero), 0.0, deltaF, &lalStrainUnit, n);
        if (!(hctilde)){ XLAL_ERROR(XLAL_EFUNC);}
        memset((*hctilde)->data->data, 0, n * sizeof(COMPLEX16));
        XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);
        init = 1;
      }
      /* Skip 22 mode if it was only required for the mixing of the 32 */
      if(ell == 2 && emm == 2 && add22 == 0){
        status = XLAL_SUCCESS;
      }
      /* Add the hl-m mode to hptilde and hctilde. */
      else{     // According to LAL documentation the azimuthal angle phi = pi/2 - phiRef. This is not taken into account in PhenomD and that's why there is an dephasing of Pi/2 between PhenomD and SEOBNRv4
          if(posMode==1 && negMode!=1){
            status = IMRPhenomXHMFDAddMode(*hptilde, *hctilde, htildelm, inclination, LAL_PI_2 , ell, emm, 0);  // add only the positive mode
          }
          else if(posMode!=1 && negMode==1){
            status = IMRPhenomXHMFDAddMode(*hptilde, *hctilde, htildelm, inclination, LAL_PI_2 , ell, -emm, 0);  // add only the negative mode
          }
          else{
            status = IMRPhenomXHMFDAddMode(*hptilde, *hctilde, htildelm, inclination, LAL_PI_2 , ell, emm, sym);    // add both positive and negative modes
          }
    }
    XLALDestroyCOMPLEX16FrequencySeries(htildelm);
  }//Loop over emm
}//Loop over ell
XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHM_Multimode failed to generate IMRPhenomXHM waveform.");

/* Free memory */
XLALDestroyCOMPLEX16FrequencySeries(htilde22);
XLALDestroyValue(ModeArray);
XLALDestroyDict(lalParams_aux);

#if DEBUG == 1
printf("\n******Leaving IMRPhenomXHM_MultiMode*****\n");
#endif

return XLAL_SUCCESS;
}



/* Core function of XLALSimIMRPhenomXHM2, returns hptilde, hctilde corresponding to a sum of modes.
The default modes are 22, 21, 33, 32 and 44 with the respective negatives ones. */
static int IMRPhenomXHM_MultiMode2(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency domain h+ GW strain */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency domain hx GW strain */
  REAL8Sequence *freqs_In,            /**< min and max frequency [Hz] */
  IMRPhenomXWaveformStruct *pWF,      /**< waveform parameters */
  LALDict *lalParams                  /**< LALDict struct */
)
{
  #if DEBUG == 1
  printf("\n\n**** IMRPhenomXHM_MultiMode2 **** \n\n");
  #endif

  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  /* Initialize useful powers of LAL_PI */
  int status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

  /* Build the frequency array and initialize htildelm to the length of freqs. */
  REAL8Sequence *freqs;
  COMPLEX16FrequencySeries *htildelm;
  // offset is the number of frequency points between 0 and f_min.
  UINT4 offset = SetupWFArrays(&freqs, &htildelm, freqs_In, pWF, ligotimegps_zero);

  UINT4 len = freqs->length;

  #if DEBUG == 1
  printf("\n\nfstart, fend, lenfreqs lenhtildelm = %.16f %.16f %i %i\n\n", freqs->data[0], freqs->data[len-1], len, htildelm->data->length);
  #endif

  // Allocate qnm struct, it contains ringdown and damping frequencies.
  QNMFits *qnms = (QNMFits *) XLALMalloc(sizeof(QNMFits));
  IMRPhenomXHM_Initialize_QNMs(qnms);

  UINT4 initial_status = XLAL_SUCCESS;

  // Variable to transform Hz in adimensional frequencies Mf:  Mf = Msec * Hz.
  REAL8 Msec  = pWF->M_sec;

  /* Transform the frequency array to adimensional frequencies to evaluate the model.
  Compute and array of the structure IMRPhenomX_UsefulPowers, with the useful powers of each frequency. */
  REAL8 *Mf = (REAL8 *)XLALMalloc(len * sizeof(REAL8));
  IMRPhenomX_UsefulPowers *powers_of_Mf = (IMRPhenomX_UsefulPowers *)XLALMalloc(len * sizeof(IMRPhenomX_UsefulPowers));

  for (UINT4 idx = 0; idx < len; idx++){
    Mf[idx] = Msec * freqs->data[idx];
    initial_status     = IMRPhenomX_Initialize_Powers(&(powers_of_Mf[idx]), Mf[idx]);
    if(initial_status != XLAL_SUCCESS)
    {
      status = initial_status;
      XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
    }
  }

  INT4 sym = 1; /* sym will decide whether to add the m mode (when equatorial symmetry is present) */

  /* Use an auxiliar laldict to not overwrite the input argument */
  LALDict *lalParams_aux;
  /* setup mode array */
  if (lalParams == NULL)
  {
      lalParams_aux = XLALCreateDict();
  }
  else{
      lalParams_aux = XLALDictDuplicate(lalParams);
  }
  lalParams_aux = IMRPhenomXHM_setup_mode_array(lalParams_aux);
  LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);

  /* At this point ModeArray should contain the list of modes
  and therefore if NULL then something is wrong and abort. */
  if (ModeArray == NULL)
  {
    XLAL_ERROR(XLAL_EDOM, "ModeArray is NULL when it shouldn't be. Aborting.\n");
  }

  /* Allocate hptilde and hctilde */
  /* Coalescence time is fixed to t=0, shift by overall length in time. */
  //XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / pWF->deltaF), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.", pWF->deltaF);
  size_t n = (htildelm)->data->length;
  XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / pWF->deltaF), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.", pWF->deltaF);
  *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &(ligotimegps_zero), 0.0, pWF->deltaF, &lalStrainUnit, n);
  if (!(hptilde)){   XLAL_ERROR(XLAL_EFUNC);}
  memset((*hptilde)->data->data, 0, n * sizeof(COMPLEX16));  // what is this for??
  XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit); // what does it do?

  *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform",  &(ligotimegps_zero), 0.0, pWF->deltaF, &lalStrainUnit, n);
  if (!(hctilde)){ XLAL_ERROR(XLAL_EFUNC);}
  memset((*hctilde)->data->data, 0, n * sizeof(COMPLEX16));
  XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

  #if DEBUG == 1
  printf("\nlength htildelm = %zu\n", n);
  #endif

  /******* Compute and add 22 Mode if active **********/
  COMPLEX16FrequencySeries *htilde22 = NULL;
  INT4 pos22, neg22;
  pos22 = XLALSimInspiralModeArrayIsModeActive(ModeArray, 2, 2);
  neg22 = XLALSimInspiralModeArrayIsModeActive(ModeArray, 2, -2);
  if (pos22 == 1 || neg22 == 1){
    #if DEBUG == 1
    printf("\n\nCalling IMRPhenomXASFrequencySequence...\n\n");
    #endif
    COMPLEX16FrequencySeries *htilde22tmp = NULL;
    XLALSimIMRPhenomXASFrequencySequence(&htilde22tmp, freqs, pWF->m1_SI, pWF->m2_SI, pWF->chi1L, pWF->chi2L, pWF->distance, pWF->phiRef_In, pWF->fRef, lalParams_aux);
    //If htilde22tmp is shorter than hptilde, need to resize
    htilde22 = XLALCreateCOMPLEX16FrequencySeries("htilde22: FD waveform", &(ligotimegps_zero), 0.0, pWF->deltaF, &lalStrainUnit, n);
    for(UINT4 idx = 0; idx < offset; idx++)
    {
      (htilde22->data->data)[idx] = 0.;
    }
    for(UINT4 idx = 0; idx < ((htilde22tmp)->data->length); idx++)
    {
      ((htilde22)->data->data)[idx+offset] = ((htilde22tmp)->data->data)[idx];
    }
    for(UINT4 idx = ((htilde22tmp)->data->length); idx < ((htildelm)->data->length) - offset; idx++)
    {
      (htilde22->data->data)[idx+offset] = 0.;
    }
    if(pos22==1 && neg22!=1){
      status = IMRPhenomXHMFDAddMode(*hptilde, *hctilde, htilde22, pWF->inclination, LAL_PI_2, 2, 2, 0);  // add only the positive mode
    }
    else if(pos22!=1 && neg22==1){
      status = IMRPhenomXHMFDAddMode(*hptilde, *hctilde, htilde22, pWF->inclination, LAL_PI_2, 2, -2, 0);  // add only the negative mode
    }
    else{
      status = IMRPhenomXHMFDAddMode(*hptilde, *hctilde, htilde22, pWF->inclination, LAL_PI_2, 2, 2, sym); // add both positive and negative modes
    }

    #if DEBUG == 1
    //Save 22 mode in file
    printf("phifRef22=%f\n",pWF->phifRef);
    FILE *file22;
    char fileSpec22[40];

    sprintf(fileSpec22, "simulation%i_MM2.dat", 22);
    printf("\nOutput file: %s\r\n",fileSpec22);
    file22 = fopen(fileSpec22,"w");

    fprintf(file22,"# q = %.16f m1 = %.16f m2 = %.16f chi1 = %.16f chi2 = %.16f lm = %i Mtot = %.16f distance = %.16f\n", pWF->q, pWF->m1, pWF->m2, pWF->chi1L, pWF->chi2L, 22, pWF->Mtot, pWF->distance/LAL_PC_SI/pow(10.,6));
    fprintf(file22,"# Frequency (Hz)    Real    Imaginary\n");

    COMPLEX16 data22;
    for(UINT4 idx = 0; idx < ((htilde22)->data->length); idx++)
    {
      data22 = ((htilde22)->data->data)[idx];
      fprintf(file22, "%.16f  %.16e %.16e\n",  idx*pWF->deltaF, creal(data22), cimag(data22));
    }
    fclose(file22);
    printf("\n lenhtilde22 =   %i\n",  htilde22->data->length);
    #endif

    XLALDestroyCOMPLEX16FrequencySeries(htilde22tmp);
  } // End of 22 mode


  #if DEBUG == 1
  printf("\n\nlenfreqs lenhtildelm =  %i %i\n\n",  len, htildelm->data->length);
  #endif

  // Initialize Amplitude and phase coefficients of the 22. pAmp22 will be filled only for the 32. pPhase22 is used for all the modes, that is why is computed outside the loop.
  IMRPhenomXAmpCoefficients *pAmp22=(IMRPhenomXAmpCoefficients *) XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
  IMRPhenomXPhaseCoefficients *pPhase22=(IMRPhenomXPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
  IMRPhenomXGetPhaseCoefficients(pWF, pPhase22);


  /****** Loop over higher modes ******/
  REAL8 amp, phi;  // Amplitude and phase of the hlm mode
  INT4 minus1l, posMode, negMode;
  REAL8 Amp0;
  for (UINT4 ell = 2; ell <= L_MAX; ell++)
  {
    /* Multiply by (-1)^l to get the true h_l-m(f) */
    if(ell%2 != 0){
      minus1l = -1;
    }
    else{
      minus1l = 1;
    }
    Amp0 = minus1l * pWF->ampNorm * pWF->amp0;

    for (INT4 emm = 1; emm < (INT4)ell + 1; emm++){
      /* Loop over only positive m is intentional. The single mode function returns the negative mode h_l-m, and the positive is added automatically in IMRPhenomXHMFDAddMode. */
      /* First check if (l,m) mode is 'activated' in the ModeArray */
      /* if activated then generate the mode, else skip this mode. */
      posMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm);
      negMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, -emm);
      if ((posMode != 1 && negMode != 1) || (ell == 2 && abs(emm) == 2))
      { /* skip mode */
        continue;
      } /* else: generate mode */

      /* We test for hypothetical m=0 modes */
      if (emm == 0)
      {
        sym = 0;
      }
      else
      {
        sym = 1;
      }

      /* Now build the corresponding hlm mode */

      // Populate pWFHM with useful parameters of each mode
      IMRPhenomXHMWaveformStruct *pWFHM = (IMRPhenomXHMWaveformStruct *) XLALMalloc(sizeof(IMRPhenomXHMWaveformStruct));
      IMRPhenomXHM_SetHMWaveformVariables(ell, emm, pWFHM, pWF, qnms, lalParams_aux);


      // If mode is even and black holes are equal we skip this part and return an array of zeros.
      if(pWFHM->Ampzero==0){

        #if DEBUG == 1
        printf("\n\nInitializing amplitude struct of %i...\n\n",pWFHM->modeTag);
        #endif

        /* Allocate and initialize the PhenomXHM lm amplitude and phase coefficients struct */
        IMRPhenomXHMAmpCoefficients *pAmp;
        pAmp = XLALMalloc(sizeof(IMRPhenomXHMAmpCoefficients));
        IMRPhenomXHMPhaseCoefficients *pPhase;
        pPhase = XLALMalloc(sizeof(IMRPhenomXHMPhaseCoefficients));
        IMRPhenomXHM_FillAmpFitsArray(pAmp);
        IMRPhenomXHM_FillPhaseFitsArray(pPhase);

        /* Get coefficients for Amplitude and Phase */
        if (pWFHM->MixingOn == 1)  {
          // For the 32 we need the speroidal coefficients for the phase and the amplitude coefficients of the 22.
          GetSpheroidalCoefficients(pPhase, pPhase22, pWFHM, pWF);
          IMRPhenomXGetAmplitudeCoefficients(pWF, pAmp22);
        }
        IMRPhenomXHM_GetAmplitudeCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF);
        IMRPhenomXHM_GetPhaseCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF,lalParams_aux);

        #if DEBUG == 1
        printf("\n\n **** Amplitude and Phase struct initialized. **** \n\n");
        #endif

        /* Loop over frequencies to generate waveform */
        /* With mode mixng */
        if(pWFHM->MixingOn==1){
          // If the 22 mode has been already computed we use it for the mixing of the 32.
          if(pos22 == 1 || neg22 == 1){
            for (UINT4 idx = 0; idx < len; idx++)
            {
              COMPLEX16 wf22 = htilde22->data->data[idx + offset]; //This will be rescaled inside SpheroidalToSphericalRecycle for the rotation
              amp = IMRPhenomXHM_Amplitude_ModeMixingRecycle(Mf[idx], &powers_of_Mf[idx], wf22, pAmp, pPhase, pWFHM);
              phi = IMRPhenomXHM_Phase_ModeMixingRecycle(Mf[idx], &powers_of_Mf[idx], wf22, pAmp, pPhase, pWFHM);
              /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
              ((htildelm)->data->data)[idx+offset] = Amp0 * amp * cexp(I * phi);
            }
          }
          // If the 22 has not been computed, its ringdown part is computed internally using pAmp22 and pPhase22.
          else{
            for (UINT4 idx = 0; idx < len; idx++)
            {
              amp = IMRPhenomXHM_Amplitude_ModeMixing(Mf[idx], &powers_of_Mf[idx], pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF);
              phi = IMRPhenomXHM_Phase_ModeMixing(Mf[idx], &powers_of_Mf[idx], pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF);
              /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
              ((htildelm)->data->data)[idx+offset] = Amp0 * amp * cexp(I * phi);
            }
          }
        }     /* No mode mixing */
        else{
          for (UINT4 idx = 0; idx < len; idx++)
          {
            amp = IMRPhenomXHM_Amplitude_noModeMixing(Mf[idx], &powers_of_Mf[idx], pAmp, pWFHM);
            phi = IMRPhenomXHM_Phase_noModeMixing(Mf[idx], &powers_of_Mf[idx], pPhase, pWFHM, pWF);
            /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
            ((htildelm)->data->data)[idx+offset] = Amp0 * amp * cexp(I * phi);
          }
        }/* End of loop over frequencies */
        /* In GetPhaseCoefficients we call IMRPhenomX_Phase_22_ConnectionCoefficients so the coefficients below are initialized,
         however every time a new mode is called we need to have these coefficients to be initially zero.*/
        pPhase22->C1Int = 0;
        pPhase22->C2Int = 0;
        pPhase22->C1MRD = 0;
        pPhase22->C2MRD = 0;

        #if DEBUG == 1
        ParametersToFile(pWF, pWFHM, pAmp, pPhase);
        #endif
        /* Free memory */
        LALFree(pAmp);
        LALFree(pPhase);
      }
      // Return array of zeros if the mode is zero
      else{
        for (UINT4 idx = 0; idx < len; idx++)
        {
          ((htildelm)->data->data)[idx+offset] = 0.;
        }
      }// end Amp zero


      #if DEBUG == 1
      // Save the hlm mode into a file
      FILE *file;
      char fileSpec[40];

      sprintf(fileSpec, "simulation%i_MM2.dat", pWFHM->modeTag);
      printf("\nOutput file: %s\r\n",fileSpec);
      file = fopen(fileSpec,"w");

      fprintf(file,"# q = %.16f m1 = %.16f m2 = %.16f chi1 = %.16f chi2 = %.16f lm = %i Mtot = %.16f distance = %.16f\n", pWF->q, pWF->m1, pWF->m2, pWF->chi1L, pWF->chi2L, pWFHM->modeTag, pWF->Mtot, pWF->distance/LAL_PC_SI/pow(10.,6));
      fprintf(file,"# Frequency (Hz)    Amplitude    Phase\n");

      COMPLEX16 data0;
      for(UINT4 idx = 0; idx < htildelm->data->length; idx++)
      {
        data0 = ((htildelm)->data->data)[idx];
        fprintf(file, "%.16f  %.16e %.16e\n",  idx*pWF->deltaF, creal(data0), cimag(data0));
      }
      fclose(file);
      #endif

      // Add the recent mode to hptilde and hctilde. beta is the azimuthal angle = pi/2 - phiRef, it is computed in IMRPhenomX_internals.c
      if(posMode==1 && negMode!=1){
        status = IMRPhenomXHMFDAddMode(*hptilde, *hctilde, htildelm, pWF->inclination, LAL_PI_2, ell, emm, 0);  // add only the positive mode
      }
      else if(posMode!=1 && negMode==1){
        status = IMRPhenomXHMFDAddMode(*hptilde, *hctilde, htildelm, pWF->inclination, LAL_PI_2, ell, -emm, 0);  // add only the negative mode
      }
      else{
        status = IMRPhenomXHMFDAddMode(*hptilde, *hctilde, htildelm, pWF->inclination, LAL_PI_2, ell, emm, sym); // add both positive and negative modes
      }
      LALFree(pWFHM);
    }
  }// End loop of higher modes


  /* Free allocated memory */
  XLALDestroyCOMPLEX16FrequencySeries(htilde22);
  XLALDestroyCOMPLEX16FrequencySeries(htildelm);
  XLALDestroyREAL8Sequence(freqs);
  XLALDestroyValue(ModeArray);
  LALFree(powers_of_Mf);
  LALFree(pAmp22);
  LALFree(pPhase22);
  LALFree(qnms);
  LALFree(Mf);
  XLALDestroyDict(lalParams_aux);


  return status;
}


/* Function to sum one mode (htildelm) to hp/c tilde */
static int IMRPhenomXHMFDAddMode(
  COMPLEX16FrequencySeries *hptilde,  /**<[out] hp series*/
  COMPLEX16FrequencySeries *hctilde,  /**<[out] hc series */
  COMPLEX16FrequencySeries *htildelm, /**< hlm mode to add */
  REAL8 theta,                        /**< Inclination [rad] */
  REAL8 phi,                          /**< Azimuthal angle [rad]*/
  INT4 l,                             /**< l index of the lm mode */
  INT4 m,                             /**< m index of the lm mode */
  INT4 sym                            /**< Equatorial symmetry */
){
  COMPLEX16 Ystar, Ym;
  UINT4 j;
  COMPLEX16 hlm; /* helper variable that contains a single point of hlmtilde */

  INT4 minus1l; /* (-1)^l */
  if (l % 2 !=0)
  {
    minus1l = -1;
  }
  else
  {
    minus1l = 1;
  }

  Ym = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m);

  if (sym)
  {
	/* Equatorial symmetry: add in -m and m mode */
    Ystar = conj(XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m));
    COMPLEX16 factorp = 0.5 * (Ym + minus1l * Ystar);
    COMPLEX16 factorc = I * 0.5 * ( Ym - minus1l * Ystar);

	#if DEBUG == 1
    printf("\nfactorp = %.16e", cabs(factorp));
    printf("\nfactorc = %.16e",cabs(factorc));
    printf("\nhtildelm length = %i\n", htildelm->data->length);
    printf("\nhp/hc tilde length = %i %i\n", hptilde->data->length,hptilde->data->length);
	#endif

    for (j = 0; j < htildelm->data->length; ++j) {
      hlm = htildelm->data->data[j];
      hptilde->data->data[j] += (factorp * hlm);
      hctilde->data->data[j] += (factorc * hlm);
    }
  }
  else  // only positives or negative modes
  {
     COMPLEX16 Ylm, factorp, factorc;
     if (m >0){  // only m>0
       Ylm = conj(XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m));
       factorp = 0.5*minus1l*Ylm;
       factorc = -I*0.5*minus1l*Ylm;
     }
     else{   // only m<0
       Ylm = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
       factorp = 0.5*Ylm;
       factorc = I*0.5*Ylm;
     }

    COMPLEX16* datap = hptilde->data->data;
    COMPLEX16* datac = hctilde->data->data;

    for ( j = 0; j < htildelm->data->length; ++j ) {
      hlm = (htildelm->data->data[j]);
      datap[j] += factorp*hlm;
      datac[j] += factorc*hlm;
    }
  }
  return XLAL_SUCCESS;
}

/** @addtogroup LALSimIMRPhenomX_c
* @{
*
* @name Routines for IMRPhenomXHM
* @{
*
* @author Cecilio García Quirós, Marta Colleoni, Sascha Husa
*
* @brief C code for IMRPhenomXHM phenomenological waveform model.
*
* This is an aligned-spin frequency domain model wit subdomimant modes: 2|2|, 2|1|, 3|3|, 3|2|, 4|4| .
* See C.García-Quirós et al for details. Any studies that use this waveform model should include
* a reference to this paper.
*
* @note The model was calibrated to mass-ratios from 1 to 1000.
* The calibration points will be given in forthcoming papers.
*
* @attention The model is usable outside this parameter range,
* and in tests to date gives sensible physical results,
* but conclusive statements on the physical fidelity of
* the model for these parameters await comparisons against further
* numerical-relativity simulations. For more information, see the review wiki
* under https://git.ligo.org/waveforms/reviews/imrphenomx/wikis/home.
* DCC link to the paper and supplementary material: https://dcc.ligo.org/P2000011-v2
*
**/

/** Returns amplitude of one single mode. It does not support the 22 mode. */
int XLALSimIMRPhenomXHMAmplitude(
    REAL8FrequencySeries **amplitude, /**< [out] FD amp */
    REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
    REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
    REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
    REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
    UINT4 ell,                           /**< l index of the mode */
    INT4 emm,                            /**< m index of the mode */
    REAL8 distance,                      /**< Luminosity distance (m) */
    REAL8 f_min,                         /**< Starting GW frequency (Hz) */
    REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
    REAL8 deltaF,                        /**< Sampling frequency (Hz) */
    REAL8 phiRef,                        /**< Orbital amp at fRef (rad) */
    REAL8 fRef_In,                       /**< Reference frequency (Hz) */
    LALDict *lalParams                   /**< Extra params */
  )
  {
    UINT4 status;

    #if DEBUG == 1
    printf("\nI am in %i %i mode\n",ell,emm);
    printf("fRef_In : %e\n",fRef_In);
    printf("m1_SI   : %e\n",m1_SI);
    printf("m2_SI   : %e\n",m2_SI);
    printf("chi1L   : %e\n",chi1L);
    printf("chi2L   : %e\n\n",chi2L);
    printf("Performing sanity checks...\n");
    #endif

    /* Sanity checks */
    if(*amplitude)       { XLAL_CHECK(NULL != amplitude, XLAL_EFAULT);                                   }
    if(fRef_In  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");  }
    if(deltaF   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaF must be positive.\n");                         }
    if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                             }
    if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                             }
    if(f_min    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "f_min must be positive.\n");                          }
    if(f_max    <  0.0) { XLAL_ERROR(XLAL_EDOM, "f_max must be non-negative.\n");                      }
    if(distance <  0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");    }
    /*
      Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
        - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
        - For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
        - For mass ratios > 1000: throw a hard error that model is not valid.
        - For spins > 0.99: throw a warning that we are extrapolating the model to extremal

    */
    REAL8 mass_ratio;
    if(m1_SI > m2_SI)
    {
      mass_ratio = m1_SI / m2_SI;
    }
    else
    {
      mass_ratio = m2_SI / m1_SI;
    }
    if(mass_ratio > 20.0  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
    if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); } // The 1e-12 is to avoid rounding errors
    if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }


      /* Set LIGOTimeGPS */
      LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

      /* Use an auxiliar laldict to not overwrite the input argument */
      LALDict *lalParams_aux;
      /* setup mode array */
      if (lalParams == NULL)
      {
          lalParams_aux = XLALCreateDict();
      }
      else{
          lalParams_aux = XLALDictDuplicate(lalParams);
      }
      lalParams_aux = IMRPhenomXHM_setup_mode_array(lalParams_aux);
      LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);

      /* first check if (l,m) mode is 'activated' in the ModeArray */
      /* if activated then generate the mode, else skip this mode. */
      if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm) != 1)
      { /* skip mode */
        XLALPrintError("XLAL Error - %i%i mode is not included\n", ell, emm);
        XLAL_ERROR(XLAL_EDOM);
      } /* else: generate mode */
      XLALDestroyValue(ModeArray);


      /* If no reference frequency is given, we will set it to the starting gravitational wave frequency */
      REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;

      /* Initialize the useful powers of LAL_PI */
      status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
      status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
      XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

      /* Initialize IMRPhenomX Waveform struct and check that it generated successfully */
      IMRPhenomXWaveformStruct *pWF;
      pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
      status = IMRPhenomXSetWaveformVariables(pWF,m1_SI, m2_SI, chi1L, chi2L, deltaF, fRef, phiRef, f_min, f_max, distance, 0.0, lalParams_aux, DEBUG);
      XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");


      /* Check that the frequency array will be consistent: fmin < fmax_prime */
      /* Return the closest power of 2 */
      size_t npts = NextPow2(pWF->f_max_prime / deltaF) + 1;
      /* Frequencies will be set using only the lower and upper bounds that we passed */
      size_t iStart = (size_t) (f_min / deltaF);
      size_t iStop  = (size_t) (pWF->f_max_prime / deltaF);
      XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
      "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);


      /*
      Create a REAL8 frequency series.
      Use fLow, fHigh, deltaF to compute frequency sequence. Only pass the boundaries (fMin, fMax).
      */
      REAL8Sequence *freqs_In = XLALCreateREAL8Sequence(2);
      freqs_In->data[0] = pWF->fMin;
      freqs_In->data[1] = pWF->f_max_prime;

      REAL8Sequence *freqs;
      UINT4 offset = SetupWFArraysReal(&freqs, amplitude, freqs_In, pWF, ligotimegps_zero);

      // allocate qnm struct
      QNMFits *qnms = (QNMFits *) XLALMalloc(sizeof(QNMFits));
      IMRPhenomXHM_Initialize_QNMs(qnms);
      // Populate pWFHM
      IMRPhenomXHMWaveformStruct *pWFHM = (IMRPhenomXHMWaveformStruct *) XLALMalloc(sizeof(IMRPhenomXHMWaveformStruct));
      IMRPhenomXHM_SetHMWaveformVariables(ell, abs(emm), pWFHM, pWF, qnms, lalParams_aux);
      LALFree(qnms);

      //populate coefficients of 22 mode, to rotate to spherical
      IMRPhenomXAmpCoefficients *pAmp22=(IMRPhenomXAmpCoefficients *) XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
      IMRPhenomXPhaseCoefficients *pPhase22=(IMRPhenomXPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));

      /* Allocate and initialize the PhenomXHM lm amplitude coefficients struct */
      IMRPhenomXHMAmpCoefficients *pAmp = (IMRPhenomXHMAmpCoefficients*) XLALMalloc(sizeof(IMRPhenomXHMAmpCoefficients));
      IMRPhenomXHMPhaseCoefficients *pPhase = (IMRPhenomXHMPhaseCoefficients*) XLALMalloc(sizeof(IMRPhenomXHMPhaseCoefficients));

      /* Allocate and initialize the PhenomXHM lm phase and amp coefficients struct */
      IMRPhenomXHM_FillAmpFitsArray(pAmp);


      /* Get coefficients for Amplitude and phase */
      if (pWFHM->MixingOn == 1)  {
        IMRPhenomXHM_FillPhaseFitsArray(pPhase);
        IMRPhenomXGetPhaseCoefficients(pWF, pPhase22);
        GetSpheroidalCoefficients(pPhase, pPhase22, pWFHM, pWF);
        IMRPhenomXGetAmplitudeCoefficients(pWF, pAmp22);
      }
      IMRPhenomXHM_GetAmplitudeCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF);

      REAL8 Amp0 = pWFHM->Amp0, amp;
      IMRPhenomX_UsefulPowers powers_of_Mf;
      /* Loop over frequencies to generate waveform */
      /* Modes with mixing */
      if(pWFHM->MixingOn==1){

        REAL8 Mf;
        for (UINT4 idx = 0; idx < freqs->length; idx++)
        {
          Mf    = pWF->M_sec * freqs->data[idx];
          INT4 initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf, Mf);
          if(initial_status != XLAL_SUCCESS)
          {
            status = initial_status;
            XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
          }
          else
          {
            amp = IMRPhenomXHM_Amplitude_ModeMixing(Mf, &powers_of_Mf, pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF);
            ((*amplitude)->data->data)[idx+offset] = Amp0 * amp;
          }
        }
      }  /* Modes without mixing */
      else{
        for (UINT4 idx = 0; idx < freqs->length; idx++)
        {
          REAL8 Mf    = pWF->M_sec * freqs->data[idx];
          INT4 initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
          if(initial_status != XLAL_SUCCESS)
          {
            status = initial_status;
            XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
          }
          else
          {
            amp = IMRPhenomXHM_Amplitude_noModeMixing(Mf, &powers_of_Mf, pAmp, pWFHM);
            ((*amplitude)->data->data)[idx+offset] = Amp0 * amp;
          }
        }
      }

      /* Resize amplitude if needed */
      if (pWF->f_max_prime < pWF->fMax)
      {
        /*
        The user has requested a higher f_max than Mf = fCut.
        Resize the frequency series to fill with zeros beyond the cutoff frequency.
        */
        size_t n = (*amplitude)->data->length;

        // We want to have the length be a power of 2 + 1
        size_t n_full = NextPow2(pWF->fMax / pWF->deltaF) + 1;

        /* Resize the COMPLEX16 frequency series */
        *amplitude = XLALResizeREAL8FrequencySeries(*amplitude, 0, n_full);
        XLAL_CHECK (*amplitude, XLAL_ENOMEM, "Failed to resize waveform REAL8FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );
      }

      LALFree(pWFHM);
      LALFree(pWF);
      LALFree(pAmp22);
      LALFree(pAmp);
      LALFree(pPhase22);
      LALFree(pPhase);
      XLALDestroyDict(lalParams_aux);

      return XLAL_SUCCESS;
  }

/** Returns phase of one single mode. It does not support the 22 mode. */
  int XLALSimIMRPhenomXHMPhase(
      REAL8FrequencySeries **phase,        /**< [out] FD amp */
      REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
      REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
      REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
      REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
      UINT4 ell,                           /**< l index of the mode */
      INT4 emm,                            /**< m index of the mode */
      REAL8 distance,                      /**< Luminosity distance (m) */
      REAL8 f_min,                         /**< Starting GW frequency (Hz) */
      REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
      REAL8 deltaF,                        /**< Sampling frequency (Hz) */
      REAL8 phiRef,                        /**< Orbital amp at fRef (rad) */
      REAL8 fRef_In,                       /**< Reference frequency (Hz) */
      LALDict *lalParams                   /**< Extra params */
    )
    {
      UINT4 status;

      #if DEBUG == 1
      printf("\nI am in %i %i mode\n",ell,emm);
      printf("fRef_In : %e\n",fRef_In);
      printf("m1_SI   : %e\n",m1_SI);
      printf("m2_SI   : %e\n",m2_SI);
      printf("chi1L   : %e\n",chi1L);
      printf("chi2L   : %e\n\n",chi2L);
      printf("Performing sanity checks...\n");
      #endif

      /* Sanity checks */
      if(*phase)       { XLAL_CHECK(NULL != phase, XLAL_EFAULT);                                   }
      if(fRef_In  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");  }
      if(deltaF   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaF must be positive.\n");                         }
      if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                             }
      if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                             }
      if(f_min    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "f_min must be positive.\n");                          }
      if(f_max    <  0.0) { XLAL_ERROR(XLAL_EDOM, "f_max must be non-negative.\n");                      }
      if(distance <  0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");    }
      /*
        Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
          - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
          - For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
          - For mass ratios > 1000: throw a hard error that model is not valid.
          - For spins > 0.99: throw a warning that we are extrapolating the model to extremal

      */
      REAL8 mass_ratio;
      if(m1_SI > m2_SI)
      {
        mass_ratio = m1_SI / m2_SI;
      }
      else
      {
        mass_ratio = m2_SI / m1_SI;
      }
      if(mass_ratio > 20.0  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
      if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); } // The 1e-12 is to avoid rounding errors
      if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }


        /* Set LIGOTimeGPS */
        LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

        /* Use an auxiliar laldict to not overwrite the input argument */
        LALDict *lalParams_aux;
        /* setup mode array */
        if (lalParams == NULL)
        {
            lalParams_aux = XLALCreateDict();
        }
        else{
            lalParams_aux = XLALDictDuplicate(lalParams);
        }
        lalParams_aux = IMRPhenomXHM_setup_mode_array(lalParams_aux);
        LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);

        /* first check if (l,m) mode is 'activated' in the ModeArray */
        /* if activated then generate the mode, else skip this mode. */
        if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm) != 1)
        { /* skip mode */
          XLALPrintError("XLAL Error - %i%i mode is not included\n", ell, emm);
          XLAL_ERROR(XLAL_EDOM);
        } /* else: generate mode */
        XLALDestroyValue(ModeArray);


        /* If no reference frequency is given, we will set it to the starting gravitational wave frequency */
        REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;

        /* Initialize the useful powers of LAL_PI */
        status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
        status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

        /* Initialize IMRPhenomX Waveform struct and check that it generated successfully */
        IMRPhenomXWaveformStruct *pWF;
        pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
        status = IMRPhenomXSetWaveformVariables(pWF,m1_SI, m2_SI, chi1L, chi2L, deltaF, fRef, phiRef, f_min, f_max, distance, 0.0, lalParams_aux, DEBUG);
        XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");


        /* Check that the frequency array will be consistent: fmin < fmax_prime */
        /* Return the closest power of 2 */
        size_t npts = NextPow2(pWF->f_max_prime / deltaF) + 1;
        /* Frequencies will be set using only the lower and upper bounds that we passed */
        size_t iStart = (size_t) (f_min / deltaF);
        size_t iStop  = (size_t) (pWF->f_max_prime / deltaF);
        XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
        "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);


        /*
        Create a REAL8 frequency series.
        Use fLow, fHigh, deltaF to compute frequency sequence. Only pass the boundaries (fMin, fMax).
        */
        REAL8Sequence *freqs_In = XLALCreateREAL8Sequence(2);
        freqs_In->data[0] = pWF->fMin;
        freqs_In->data[1] = pWF->f_max_prime;

        REAL8Sequence *freqs;
        UINT4 offset = SetupWFArraysReal(&freqs, phase, freqs_In, pWF, ligotimegps_zero);

        // allocate qnm struct
        QNMFits *qnms = (QNMFits *) XLALMalloc(sizeof(QNMFits));
        IMRPhenomXHM_Initialize_QNMs(qnms);
        // Populate pWFHM
        IMRPhenomXHMWaveformStruct *pWFHM = (IMRPhenomXHMWaveformStruct *) XLALMalloc(sizeof(IMRPhenomXHMWaveformStruct));
        IMRPhenomXHM_SetHMWaveformVariables(ell, abs(emm), pWFHM, pWF, qnms, lalParams_aux);
        LALFree(qnms);

        //populate coefficients of 22 mode, to rotate to spherical
        IMRPhenomXAmpCoefficients *pAmp22=(IMRPhenomXAmpCoefficients *) XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
        IMRPhenomXPhaseCoefficients *pPhase22=(IMRPhenomXPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
        IMRPhenomXGetPhaseCoefficients(pWF, pPhase22);

        /* Allocate and initialize the PhenomXHM lm amplitude coefficients struct */
        IMRPhenomXHMAmpCoefficients *pAmp = (IMRPhenomXHMAmpCoefficients*) XLALMalloc(sizeof(IMRPhenomXHMAmpCoefficients));
        IMRPhenomXHMPhaseCoefficients *pPhase = (IMRPhenomXHMPhaseCoefficients*) XLALMalloc(sizeof(IMRPhenomXHMPhaseCoefficients));

        /* Allocate and initialize the PhenomXHM lm phase and amp coefficients struct */
        IMRPhenomXHM_FillPhaseFitsArray(pPhase);


        /* Get coefficients for Amplitude and phase */
        if (pWFHM->MixingOn == 1)  {
          GetSpheroidalCoefficients(pPhase, pPhase22, pWFHM, pWF);
          IMRPhenomXGetAmplitudeCoefficients(pWF, pAmp22);
          IMRPhenomXHM_FillAmpFitsArray(pAmp);
          IMRPhenomXHM_GetAmplitudeCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF);
        }
        IMRPhenomXHM_GetPhaseCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF,lalParams_aux);

        IMRPhenomX_UsefulPowers powers_of_Mf;
        REAL8 addpi = 0;   // Add pi to the phase if (-1)^l is negative

        /* Multiply by (-1)^l to get the true negative mode */
        if(ell%2 != 0){
          addpi = LAL_PI;
        }
        /* Loop over frequencies to generate waveform */
        /* Modes with mixing */
        if(pWFHM->MixingOn==1){

          REAL8 Mf;
          for (UINT4 idx = 0; idx < freqs->length; idx++)
          {
            Mf    = pWF->M_sec * freqs->data[idx];
            INT4 initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
            if(initial_status != XLAL_SUCCESS)
            {
              status = initial_status;
              XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
            }
            else
            {
              REAL8 phi = IMRPhenomXHM_Phase_ModeMixing(Mf, &powers_of_Mf, pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF);
              ((*phase)->data->data)[idx+offset] = phi + addpi;
            }
          }
        }  /* Modes without mixing */
        else{
          for (UINT4 idx = 0; idx < freqs->length; idx++)
          {
            REAL8 Mf    = pWF->M_sec * freqs->data[idx];
            INT4 initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
            if(initial_status != XLAL_SUCCESS)
            {
              status = initial_status;
              XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d", initial_status);
            }
            else
            {
              REAL8 phi = IMRPhenomXHM_Phase_noModeMixing(Mf, &powers_of_Mf, pPhase, pWFHM, pWF);
              ((*phase)->data->data)[idx+offset] = phi + addpi;
            }
          }
        }

        /* If the mode is positive, the wf is (-1)^l*conjugate(wf). For the phase this is translated in adding or not pi and changing the sign. */
        if(emm > 0){
            for (UINT4 idx = 0; idx < freqs->length; idx++)
            {
                ((*phase)->data->data)[idx+offset] = addpi - ((*phase)->data->data)[idx+offset];
            }
        }

        /* Resize htildelm if needed */
        if (pWF->f_max_prime < pWF->fMax)
        {
          /*
          The user has requested a higher f_max than Mf = fCut.
          Resize the frequency series to fill with zeros beyond the cutoff frequency.
          */
          size_t n = (*phase)->data->length;

          // We want to have the length be a power of 2 + 1
          size_t n_full = NextPow2(pWF->fMax / pWF->deltaF) + 1;

          /* Resize the COMPLEX16 frequency series */
          *phase = XLALResizeREAL8FrequencySeries(*phase, 0, n_full);
          XLAL_CHECK (*phase, XLAL_ENOMEM, "Failed to resize waveform REAL8FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );
        }

        LALFree(pWFHM);
        LALFree(pWF);
        LALFree(pAmp22);
        LALFree(pAmp);
        LALFree(pPhase22);
        LALFree(pPhase);
        XLALDestroyDict(lalParams_aux);

        return XLAL_SUCCESS;
    }
/** @}
* @} **/
