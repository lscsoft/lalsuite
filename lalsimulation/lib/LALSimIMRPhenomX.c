/*
 *  Copyright (C) 2018 Geraint Pratten
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

/* Standard LAL */
#include <lal/Sequence.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/XLALError.h>

/* LAL datatypes and constants */
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

/* Time series, frequency series and spherical harmonics */
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/SphericalHarmonics.h>
#include <lal/FrequencySeries.h>

/* LALSimulation */
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>

/* Standard C */
#include <math.h>
#include <complex.h>
#include <stdbool.h>

/* GSL */
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#define PHENOMXDEBUG 0
#else
#define DEBUG 1 //print debugging info
#define PHENOMXDEBUG 1
#endif

/* Link IMRPhenomX routines */
#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomX_ringdown.h"
#include "LALSimIMRPhenomX_intermediate.h"
#include "LALSimIMRPhenomX_inspiral.h"
#include "LALSimIMRPhenomX_internals.c"
//#include "LALSimIMRPhenomX_tidal.c"
//#include "LALSimIMRPhenomX_precession.c"

/* Note: This is declared in LALSimIMRPhenomX_internals.c and avoids namespace clashes */
IMRPhenomX_UsefulPowers powers_of_lalpi;

#ifndef _OPENMP
#define omp ignore
#endif

/* ******** ALIGNED SPIN IMR PHENOMENOLOGICAL WAVEFORM: IMRPhenomXAS ********* */

/* EXTERNAL ROUTINES */

/**
 * @addtogroup LALSimIMRPhenomX_c
 * @brief Routines to produce IMRPhenomX-family of phenomenological
 * inspiral-merger-ringdown waveforms.
 *
 * These are frequency-domain models for compact binaries at comparable and extreme mass ratios,
 * tuned to numerical-relativity simulations.
 *  * IMRPhenomXAS model for 22 mode non-precessing binaries. https://arxiv.org/abs/2001.11412 ; DCC link: https://dcc.ligo.org/LIGO-P2000018
 *  * IMRPhenomXHM model with subdominant modes non-precessing binaries. https://arxiv.org/abs/2001.10914 ; DCC link: https://dcc.ligo.org/P1900393-v1
 *  * Multibanding for IMRPhenomXHM. https://arxiv.org/abs/2001.10897 ; DCC link: https://dcc.ligo.org/LIGO-P1900391
 *
 * @review IMRPhenomXAS & IMRPhenomXHM reviewed by maria.haney, patricia-schmidt, roberto.cotesta, anuradha.samajdar, jonathan.thompson, nv.krishnendu
 * Review wiki: https://git.ligo.org/waveforms/reviews/imrphenomx/wikis/home
 *
 *
 */

 /**
  * @addtogroup LALSimIMRPhenomX_c
  * @{
  *
  * @name Routines for IMRPhenomXAS
  * @{
  *
  * @author Geraint Pratten
  *
  * @brief C code for IMRPhenomXAS phenomenological waveform model.
  *
  * This is an aligned-spin frequency domain model for the 22 mode.
  * See G.Pratten et al for details. Any studies that use this waveform model should include
  * a reference to both of this paper.
  *
  * @note The model was calibrated to mass-ratios from 1 to 1000.
  * The calibration points will be given in forthcoming papers.
  *
  * @attention The model is usable outside this parameter range,
  * and in tests to date gives sensible physical results,
  * but conclusive statements on the physical fidelity of
  * the model for these parameters await comparisons against further
  * numerical-relativity simulations. For more information, see the review wiki
  * under https://git.ligo.org/waveforms/reviews/imrphenomx/wikis/home
  *
  *
  * Waveform flags:
  * 	InsPhaseVersion: Determines the inspiral phase model.
  *			- 104 : Canonical TaylorF2 at 3.5PN including cubic-in-spin and quadratic-in-spin corrections. Uses 4 pseudo-PN coefficients. (RECOMMENDED).
  *			- 105 : Canonical TaylorF2 at 3.5PN including cubic-in-spin and quadratic-in-spin corrections. Uses 5 pseudo-PN coefficients.
  *			- 114 : Extended TaylorF2. Includes cubic-in-spin, quadratic-in-spin corrections, 4PN and 4.5PN orbital corrections. Uses 4 pseudo-PN coefficients.
  *			- 115 : Extended TaylorF2. Includes cubic-in-spin, quadratic-in-spin corrections, 4PN and 4.5PN orbital corrections. Uses 5 pseudo-PN coefficients.
  *
  *     IntPhaseVersion: Determines the intermediate phase model.
  *			-	104 : 4th order polynomial ansatz.
  *			-   105 : 5th order polynomial ansatz. (RECOMMENDED).
  *
  *		RDPhaseVersion: Determines the merger-ringdown phase model.
  *			-	105 : Deformed Lorentzian using 5 coefficients. (RECOMMENDED).
  *
  *		InsAmpVersion : Determines inspiral amplitude model.
  *			-	103 : Canonical PN re-expanded TaylorF2 amplitude with pseudo-PN corrections. (RECOMMENDED).
  *
  *		IntAmpVersion : Determines intermediate amplitude model.
  *			-	104 : Based on a 4th order polynomial ansatz. Less accurate but stable extrapolation. (RECOMMENDED).
  *			-	105 : Based on 5th order polynomial ansatz. More accurate in calibration domain, more unstable extrapolation.
  *
  *		RDAmpVersion : Determines the merger-ringdown amplitude model.
  *			-	103 : Deformed Lorentzian with 3 free coefficients. Uses 1 calibrated collocation point and 2 calibrated phenomenological coefficients. (RECOMMENDED).
  */


/**
 *  Driver routine to calculate an IMRPhenomX aligned-spin,
 *  inspiral-merger-ringdown phenomenological waveform model
 *  in the frequency domain.
 *
 *  arXiv:2001.11412, https://arxiv.org/abs/2001.11412
 *
 *  All input parameters should be in SI units. Angles should be in radians.
 *
 *  XLALSimIMRPhenomXASGenerateFD() returns the strain of the 2-2 mode as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 */
int XLALSimIMRPhenomXASGenerateFD(
  COMPLEX16FrequencySeries **htilde22, /**< [out] FD waveform */
  REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,                      /**< Luminosity distance (m) */
  REAL8 f_min,                         /**< Starting GW frequency (Hz) */
  REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = 0.3 */
  REAL8 deltaF,                        /**< Sampling frequency (Hz) */
  REAL8 phi0,                          /**< Orbital phase at fRef (rad) */
  REAL8 fRef_In,                       /**< Reference frequency (Hz) */
  LALDict *lalParams                   /**< LAL Dictionary */
)
{
  UINT4 status;

  /* Set debug status here */
  UINT4 debug = PHENOMXDEBUG;

  if(debug)
  {
    printf("fRef_In : %e\n",fRef_In);
    printf("m1_SI   : %e\n",m1_SI);
    printf("m2_SI   : %e\n",m2_SI);
    printf("chi1L   : %e\n",chi1L);
    printf("chi2L   : %e\n\n",chi2L);
    printf("Performing sanity checks...\n");
  }

  /* Perform initial sanity checks */
  if(*htilde22)       { XLAL_CHECK(NULL != htilde22, XLAL_EFAULT);                                   }
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

  /* If no reference frequency is given, set it to the starting gravitational wave frequency */
  REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;


  if(debug)
  {
    printf("\n\n **** Initializing waveform struct... **** \n\n");
  }


  /* Initialize the useful powers of LAL_PI */
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

  /* Initialize IMR PhenomX Waveform struct and check that it initialized correctly */
  IMRPhenomXWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, deltaF, fRef, phi0, f_min, f_max, distance, 0.0, lalParams, debug);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

  /*
      Create a REAL8 frequency series.
      Use fLow, fHigh, deltaF to compute frequency sequence. Only pass the boundaries (fMin, fMax).
  */
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = pWF->fMin;
  freqs->data[1] = pWF->f_max_prime;


  if(debug)
  {
    printf("\n\n **** Calling IMRPhenomXASGenerateFD... **** \n\n");
  }

  /* We now call the core IMRPhenomXAS waveform generator */
  status = IMRPhenomXASGenerateFD(htilde22, freqs, pWF, lalParams);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXASFDCore failed to generate IMRPhenomX waveform.");

  if(debug)
  {
    printf("\n\n **** Call to IMRPhenomXASGenerateFD complete. **** \n\n");
  }

  /*
      We now resize htilde22 if our waveform was generated to a cut-off frequency below
      the desired maximum frequency. Simply fill the remaining frequencies with zeros.
  */
  if (pWF->f_max_prime < pWF->fMax)
  {
    /*
        As the user has requested an f_max > Mf = fCut,
        we resize the frequency series to fill with zeros beyond the cutoff frequency.
    */
    size_t n = (*htilde22)->data->length;

    /* Enforce length to be a power of 2 + 1 */
    size_t n_full = NextPow2(pWF->fMax / pWF->deltaF) + 1;

    /* Resize the COMPLEX16 frequency series */
    *htilde22 = XLALResizeCOMPLEX16FrequencySeries(*htilde22, 0, n_full);
    XLAL_CHECK (*htilde22, XLAL_ENOMEM, "Failed to resize waveform COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );
  }
  LALFree(pWF);
  XLALDestroyREAL8Sequence(freqs);
  return XLAL_SUCCESS;
}


/**
 * Compute waveform in LAL format at specified frequencies for the IMRPhenomX model.
 *
 *  All input parameters should be in SI units. Angles should be in radians.
 *
 * XLALSimIMRPhenomXASFrequencySequence() returns the strain of the 2-2 mode as a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added. Assumes positive frequencies.
 */
 int XLALSimIMRPhenomXASFrequencySequence(
     COMPLEX16FrequencySeries **htilde22, /**< [out] FD waveform */
     const REAL8Sequence *freqs,          /**< [out] Frequency series [Hz] */
     REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
     REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
     REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
     REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
     REAL8 distance,                      /**< Luminosity distance (m) */
     REAL8 phi0,                          /**< Phase at reference frequency */
     REAL8 fRef_In,                       /**< Reference frequency (Hz) */
     LALDict *lalParams                   /**< LAL Dictionary */
 )
 {
   INT4 return_code = 0;

   /* Sanity checks */
   if(*htilde22)       { XLAL_CHECK(NULL != htilde22, XLAL_EFAULT);                                   }
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

   // Check on the mass-ratio with a 1e-12 tolerance to avoid rounding errors
   if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000."); }
   if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

   // If fRef is not provided, then set fRef to be the starting GW Frequency
   REAL8 fRef = (fRef_In == 0.0) ? freqs->data[0] : fRef_In;

   UINT4 status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
   XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

   /*
      This routine automatically performs sanity checks on the masses, spins, etc.
   */
   REAL8 f_min_In  = freqs->data[0];
   REAL8 f_max_In  = freqs->data[freqs->length - 1];

   /*
      Passing deltaF = 0 implies that freqs is a frequency grid with non-uniform spacing.
      The function waveform then start at lowest given frequency.
   */

   /* Initialize IMRPhenomX waveform struct and perform sanity check. */
   IMRPhenomXWaveformStruct *pWF;
   pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
   return_code = IMRPhenomXSetWaveformVariables(pWF,m1_SI, m2_SI, chi1L, chi2L, 0.0, fRef, phi0, f_min_In, f_max_In, distance, 0.0, lalParams, 0);
   XLAL_CHECK(XLAL_SUCCESS == return_code, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

   /* Now call the core IMRPhenomX waveform generator */
   return_code = IMRPhenomXASGenerateFD(
     htilde22,
     freqs,
     pWF,
     lalParams
   );
   XLAL_CHECK(return_code == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXASFDCore failed to generate IMRPhenomX waveform.");
   LALFree(pWF);

   return XLAL_SUCCESS;
 }

 /** @} */
 /** @} */


 /* *********************************************************************************
  *
  * The following private function generates an IMRPhenomX frequency-domain waveform
  *   - Only aligned-spin
  *   - Only the 22-mode
  *   - Physical parameters are passed via the waveform struct
  * *********************************************************************************
  */
int IMRPhenomXASGenerateFD(
  COMPLEX16FrequencySeries **htilde22, /**< [out] FD waveform           */
  const REAL8Sequence *freqs_In,       /**< Input frequency grid        */
  IMRPhenomXWaveformStruct *pWF,       /**< IMRPhenomX Waveform Struct  */
  LALDict *lalParams                   /**< LAL Dictionary Structure    */
)
{
  /* Inherits debug flag from waveform struct */
  UINT4 debug = PHENOMXDEBUG;

  if(debug)
  {
    printf("\n **** Now in IMRPhenomXASGenerateFD... **** \n");
  }

  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  /* Initialize useful powers of LAL_PI */
  int status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");

  /* Inherit minimum and maximum frequencies to generate wavefom from input frequency grid */
  double f_min = freqs_In->data[0];
  double f_max = freqs_In->data[freqs_In->length - 1];

  /* Size of array */
  size_t npts     = 0;

  /* Index shift between freqs and the frequency series */
  UINT4 offset    = 0;

  /* Initialize frequency sequence */
  REAL8Sequence *freqs = NULL;

  /* If deltaF is non-zero then we need to generate a uniformly sampled frequency grid of spacing deltaF. Start at f = 0. */
  if(pWF->deltaF > 0)
  {
    /* Return the closest power of 2 */
    npts = NextPow2(f_max / pWF->deltaF) + 1;

    /* Debug information */
    if(debug)
    {
      printf("npts     = %zu\n",npts);
      printf("fMin     = %.4f\n",f_min);
      printf("fMax     = %.4f\n",f_max);
      printf("dF       = %.4f\n",pWF->deltaF);
    }

    XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / pWF->deltaF ), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.", pWF->deltaF);

    /* Initialize the htilde frequency series */
    *htilde22 = XLALCreateCOMPLEX16FrequencySeries("htilde22: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,npts);

    /* Check that frequency series generated okay */
    XLAL_CHECK(*htilde22,XLAL_ENOMEM,"Failed to allocate COMPLEX16FrequencySeries of length %zu for f_max = %f, deltaF = %g.\n",npts,f_max,pWF->deltaF);

    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t) (f_min / pWF->deltaF);
    size_t iStop  = (size_t) (f_max / pWF->deltaF) + 1;

    XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
          "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);

    /* Allocate memory for frequency array and terminate if this fails */
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!freqs)
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }

    /* Populate frequency array */
    for (UINT4 i = iStart; i < iStop; i++)
    {
      freqs->data[i-iStart] = i * pWF->deltaF;
    }
    offset = iStart;
  }
  else
  {
    /* freqs is a frequency grid with non-uniform spacing, so we start at the lowest given frequency */
    npts      = freqs_In->length;
    *htilde22 = XLALCreateCOMPLEX16FrequencySeries("htilde22: FD waveform, 22 mode", &ligotimegps_zero, f_min, pWF->deltaF, &lalStrainUnit, npts);

    XLAL_CHECK (*htilde22, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu from sequence.", npts);

    offset = 0;
    freqs  = XLALCreateREAL8Sequence(freqs_In->length);

    /* Allocate memory for frequency array and terminate if this fails */
    if (!freqs)
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }

    /* Populate frequency array */
    for (UINT4 i = 0; i < freqs_In->length; i++)
    {
      freqs->data[i] = freqs_In->data[i];
    }
  }

  memset((*htilde22)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htilde22)->sampleUnits), &((*htilde22)->sampleUnits), &lalSecondUnit);

  /* Check if LAL dictionary exists. If not, create a LAL dictionary. */
  INT4 lalParams_In = 0;
  if(lalParams == NULL)
  {
    lalParams_In = 1;
    lalParams = XLALCreateDict();
  }

  if(debug)
  {
    printf("\n\n **** Initializing amplitude struct... **** \n\n");
  }

  /* Allocate and initialize the PhenomX 22 amplitude coefficients struct */
  IMRPhenomXAmpCoefficients *pAmp22;
  pAmp22 = XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
  status = IMRPhenomXGetAmplitudeCoefficients(pWF,pAmp22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetAmplitudeCoefficients failed.\n");

  if(debug)
  {
    printf("\n\n **** Amplitude struct initialized. **** \n\n");
    printf("\n\n **** Initializing phase struct... **** \n\n");
  }

  /* Allocate and initialize the PhenomX 22 phase coefficients struct */
  IMRPhenomXPhaseCoefficients *pPhase22;
  pPhase22 = XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
  status   = IMRPhenomXGetPhaseCoefficients(pWF,pPhase22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetPhaseCoefficients failed.\n");

  if(debug)
  {
    printf("\n\n **** Phase struct initialized. **** \n\n");
  }

  /*
      Apply time shifts so peak amplitude is near t ~ 0.
  */

  /* Initialize a struct containing useful powers of Mf at fRef */
  IMRPhenomX_UsefulPowers powers_of_MfRef;
  status = IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF->MfRef);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "IMRPhenomX_Initialize_Powers failed for MfRef.\n");

  /* Linear time and phase shifts so that model peaks near t ~ 0 */
  REAL8 lina = 0;

  /* Get phase connection coefficients */
  IMRPhenomX_Phase_22_ConnectionCoefficients(pWF,pPhase22);
  double linb=IMRPhenomX_TimeShift_22(pPhase22, pWF);

  /* 1/eta is used to re-scale phase */
  REAL8 inveta    = (1.0 / pWF->eta);

  /* Calculate phase at reference frequency: phifRef = 2.0*phi0 + LAL_PI_4 + PhenomXPhase(fRef) */
  pWF->phifRef = -(inveta * IMRPhenomX_Phase_22(pWF->MfRef, &powers_of_MfRef, pPhase22, pWF) + linb*pWF->MfRef + lina) + 2.0*pWF->phi0 + LAL_PI_4;

  /*
      Here we declare explicit REAL8 variables for main loop in order to avoid numerous
      pointer calls.
  */
  //REAL8 MfRef     = pWF->MfRef;
  REAL8 Msec      = pWF->M_sec;

  REAL8 C1IM      = pPhase22->C1Int;
  REAL8 C2IM      = pPhase22->C2Int;
  REAL8 C1RD      = pPhase22->C1MRD;
  REAL8 C2RD      = pPhase22->C2MRD;

  REAL8 fPhaseIN  = pPhase22->fPhaseMatchIN;
  REAL8 fPhaseIM  = pPhase22->fPhaseMatchIM;
  REAL8 fAmpIN    = pAmp22->fAmpMatchIN;
  REAL8 fAmpIM    = pAmp22->fAmpRDMin;

  if(debug)
  {
    printf("\n\n **** Phase struct initialized. **** \n\n");
    printf("C1IM     = %.4f\n",C1IM);
    printf("C2IM     = %.4f\n",C2IM);
    printf("C1RD     = %.4f\n",C1RD);
    printf("C2RD     = %.4f\n",C2RD);
    printf("fIN      = %.4f\n",fPhaseIN);
    printf("fIM      = %.4f\n",fPhaseIM);
  }

  REAL8 Amp0      = pWF->amp0 * pWF->ampNorm;

  /* initial_status used to track  */
  UINT4 initial_status = XLAL_SUCCESS;

  /* Now loop over main driver to generate waveform:  h(f) = A(f) * Exp[I phi(f)] */
  #pragma omp parallel for
  for (UINT4 idx = 0; idx < freqs->length; idx++)
  {
    double Mf    = Msec * freqs->data[idx];   // Mf is declared locally inside the loop
    UINT4 jdx    = idx  + offset;             // jdx is declared locally inside the loop

    /* Initialize a struct containing useful powers of Mf */
    IMRPhenomX_UsefulPowers powers_of_Mf;
    initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
    if(initial_status != XLAL_SUCCESS)
    {
      status = initial_status;
      XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
    }
    else
    {
      /* Generate amplitude and phase at MfRef */
      REAL8 amp = 0.0;
      REAL8 phi = 0.0;

      /* The functions in this routine are inlined to help performance. */
      /* Construct phase */
      if(Mf < fPhaseIN)
      {
        phi = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pPhase22);
      }
      else if(Mf > fPhaseIM)
      {
        phi = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1RD + (C2RD * Mf);
      }
      else
      {
        phi = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1IM + (C2IM * Mf);
      }

	  /* Scale phase by 1/eta */
	  phi  *= inveta;
      phi  += linb*Mf + lina + pWF->phifRef;

	  /* Construct amplitude */
	  if(Mf < fAmpIN)
	  {
		  amp = IMRPhenomX_Inspiral_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
	  }
	  else if(Mf > fAmpIM)
	  {
		  amp = IMRPhenomX_Ringdown_Amp_22_Ansatz(Mf, pWF, pAmp22);
	  }
	  else
	  {
        amp = IMRPhenomX_Intermediate_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
      }

	  /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
      ((*htilde22)->data->data)[jdx] = Amp0 * powers_of_Mf.m_seven_sixths * amp * cexp(I * phi);
    }
  }

  // Free allocated memory
  LALFree(pAmp22);
  LALFree(pPhase22);
  XLALDestroyREAL8Sequence(freqs);
  if(lalParams_In == 1)
  {
    XLALDestroyDict(lalParams);
  }

  return status;
}


/* Useful wrapper to check for uniform frequency grids taken from LALSimIMRPhenomHM.c */
int IMRPhenomXCheckForUniformFrequencies(
  REAL8Sequence *frequencies,
  REAL8 df
)
{
  INT4 IsUniform = 0;

  /* If the frequency series has length of 2 and a df > 0 then it is uniformly sampled */
  if( (frequencies->length == 2) && (df > 0.) )
  {
    IsUniform = 1;
  }

  return IsUniform;
};
