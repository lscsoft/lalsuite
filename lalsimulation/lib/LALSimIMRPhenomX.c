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

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#define PHENOMXDEBUG 0
#define PHENOMXPDEBUG 0
#else
#define DEBUG 1 //print debugging info
#define PHENOMXDEBUG 1
#define PHENOMXPDEBUG 1
#endif

/* Link IMRPhenomX routines */
#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomX_ringdown.h"
#include "LALSimIMRPhenomX_intermediate.h"
#include "LALSimIMRPhenomX_inspiral.h"
#include "LALSimIMRPhenomX_internals.c"
//#include "LALSimIMRPhenomX_tidal.c"
#include "LALSimIMRPhenomX_precession.c"

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
 *  * IMRPhenomXAS model for the 22 mode of non-precessing binaries. https://arxiv.org/abs/2001.11412 ; DCC link: https://dcc.ligo.org/P2000018
 *  * IMRPhenomXHM model with subdominant modes non-precessing binaries. https://arxiv.org/abs/2001.10914 ; DCC link: https://dcc.ligo.org/P1900393
 *  * Multibanding for IMRPhenomXHM. https://arxiv.org/abs/2001.10897 ; DCC link: https://dcc.ligo.org/P1900391
 *  * IMRPhenomXP model for the 22 mode of non-precessing binaries. https://arxiv.org/abs/2004.06503 ; DCC link: https://dcc.ligo.org/P2000039
 *  * IMRPhenomXPHM model with subdominant modes for precessing binaries. https://arxiv.org/abs/2004.06503 ; DCC link: https://dcc.ligo.org/P2000039
 *
 * @review IMRPhenomXAS & IMRPhenomXHM reviewed by Maria Haney, Patricia Schmidt,
 * Roberto Cotesta, Anuradha Samajdar, Jonathan Thompson, N.V. Krishnendu.
 * IMRPhenomXP & IMRPhenomXPHM reviewed by Maria Haney, Jonathan Thompson,
 * Marta Colleoni, David Keitel.
 * Combined review wiki:
 * https://git.ligo.org/waveforms/reviews/imrphenomx/-/wikis/home
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
  *
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
  REAL8 lastfreq;
  if (pWF->f_max_prime < pWF->fMax)
  {
    /*
        As the user has requested an f_max > Mf = fCut,
        we resize the frequency series to fill with zeros beyond the cutoff frequency.
    */
    lastfreq = pWF->fMax;
  }
  else{  // We have to look for a power of 2 anyway.
    lastfreq = pWF->f_max_prime;
  }
  /* Enforce length to be a power of 2 + 1 */
  size_t n_full = NextPow2(lastfreq / pWF->deltaF) + 1;
  size_t n = (*htilde22)->data->length;

  /* Resize the COMPLEX16 frequency series */
  *htilde22 = XLALResizeCOMPLEX16FrequencySeries(*htilde22, 0, n_full);
  XLAL_CHECK (*htilde22, XLAL_ENOMEM, "Failed to resize waveform COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );


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
    npts = (size_t) (f_max/pWF->deltaF) + 1;

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




/* ******** PRECESSING IMR PHENOMENOLOGICAL WAVEFORM: IMRPhenomXP ********* */


/*
    Decleration for internal function to generate precessing-spin, 22 only IMRPhenomXP waveform.
    Declared here and not in header due to IMRPhenomXPrecessionStruct.
*/
int IMRPhenomXPGenerateFD(
  COMPLEX16FrequencySeries **hptilde,       /**< [out] FD waveform           */
  COMPLEX16FrequencySeries **hctilde,       /**< [out] FD waveform           */
  const REAL8Sequence *freqs_In,            /**< Input frequency grid        */
  IMRPhenomXWaveformStruct *pWF,            /**< IMRPhenomX Waveform Struct  */
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomXP Waveform Struct */
  LALDict *lalParams                        /**< LAL Dictionary Structure    */
);

/**
 * @addtogroup LALSimIMRPhenomX_c
 * @{
 *
 * @name Routines for IMRPhenomXP
 * @{
 *
 * @author Geraint Pratten, Cecilio García Quirós
 *
 * @brief C code for IMRPhenomXP phenomenological waveform model.
 *
 * This is a precessing frequency domain model.
 * See Pratten, García-Quirós, Colleoni et al arXiv:XXXX.YYYY for details.
*  Studies using this model are kindly asked
 * to cite Pratten et al arXiv:2001.11412, García-Quirós et al arXiv:2001.10914
 * and Pratten, García-Quirós, Colleoni et al arXiv:XXXX.YYYY.
 *
 * @note The underlying aligned-spin model was calibrated for
 * mass-ratios 1 to 1000.
 *
 * @attention The model can be called outside this parameter space
 * and in all tests to date gives sensible physical results
 * but conclusive statements on the physical fidelity of
 * the model requires further detailed comparisons against
 * numerical-relativity simulations. For more information, see the review wiki
 * under https://git.ligo.org/waveforms/reviews/imrphenomx/wikis/home
 *
 *  IMRPhenomXP/HM is based on a modular framework. User can specify flags to control how Euler angles are calculated, the final spin
 *  parameterization and the conventions used in reconstructing the waveform in the LAL frame. A detailed discussion can be
 *  found in arXiv:XXXX.YYYY. The various flags are detailed below.
 *
 *  Precession flags:
 *   PhenomXPrecVersion:
 *     - 101 : NNLO PN Euler angles and a 2PN non-spinning approximation to L
 *     - 102 : NNLO PN Euler angles and a 3PN spinning approximation to L
 *     - 103 : NNLO PN Euler angles and a 4PN spinning approximation to L
 *     - 104 : NNLO PN Euler angles and a 4PN spinning approximation to L augmeneted with leading PN order at all order in spin terms. See N. Siemonsen et al, PRD, 97, 124046, (2018), arXiv:1712.08603
 *     - 220 : MSA Euler angles and a 3PN spinning approximation to L, see K. Chatziioannou et al, PRD, 95, 104004, (2017), arXiv:1703.03967. Defaults to NNLO version 102 if MSA fails to initialize.
 *     - 221 : MSA Euler angles and a 3PN spinning approximation to L.
 *     - 222 : MSA Euler angles as implemented in LALSimInspiralFDPrecAngles.
 *     - 223 : MSA Euler angles as implemented in LALSimInspiralFDPrecAngles. Defaults to NNLO version 102 if MSA fails to initialize. [Default].
 *     - 224 : As version 220 but using the \f$\phi_{z,0}\f$ and \f$\zeta_{z,0}\f$ prescription from 223.
 *
 *   PhenomXPExpansionOrder:
 *     - -1, 0, 1, 2, 3, 4, 5. Controls the expansion order of the leading-order MSA terms for both \f$\zeta\f$ and \f$\phi_z\f$. [Default is 5].
 *
 *   PhenomXPFinalSpinMod:
 *     - 0 : Modify final spin based on \f$\chi_p\f$. [Recommended default for NNLO angles].
 *     - 1 : Modify final spin using \f$\chi_{1x}\f$. This is pathological. Do not use. Implemented to compare to PhenomPv3 pre-bug fix.
 *     - 2 : Modify final spin using norm of total in-plane spin vector.
 *     - 3 : Modify final spin using precession-averaged couplings from MSA analysis. Only works with MSA Euler angles (versions 220, 221, 222, 223 and 224). If MSA fails to initialize
 *           or called with NNLO angles, default to version 0. [Default]
 *
 *   PhenomXPConvention (App. C and Table II in App. F of arXiv:XXXX.YYYY):
 *     - 0 : Conventions defined as following https://dcc.ligo.org/LIGO-T1500602
 *     - 1 : Convention defined following App. C, see Table II of arXiv:XXXX.YYYY for specific details. [Default]
 *     - 5 : Conventions as used in PhenomPv3/HM
 *     - 6 : Conventions defined following App. C, see Table II of arXiv:XXXX.YYYY for specific details.
 *     - 7 : Conventions defined following App. C, see Table II of arXiv:XXXX.YYYY for specific details.
 */


/*
 *  Prototype wrapper function:

 *  Driver routine to calculate an IMRPhenomX precessing,
 *  inspiral-merger-ringdown phenomenological waveform model
 *  in the frequency domain.
 *
 *  All input parameters should be in SI units. Angles should be in radians.
 *
 *  Returns the plus and cross polarizations as a complex frequency series with
 *  equal spacing deltaF and contains zeros from zero frequency to the starting
 *  frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 *
 */
int XLALSimIMRPhenomXPGenerateFD(
  COMPLEX16FrequencySeries **hptilde, /**< [out] Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde, /**< [out] Frequency-domain waveform hx */
  REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
  REAL8 chi1x,                        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  const REAL8 distance,               /**< Distance of source (m) */
  const REAL8 inclination,            /**< inclination of source (rad) */
  const REAL8 phiRef,                 /**< Orbital phase (rad) at reference frequency */
  REAL8 f_min,                        /**< Starting GW frequency (Hz) */
  REAL8 f_max,                        /**< Ending GW frequency (Hz); Defaults to Mf = 0.3 if no f_max is specified. */
  const REAL8 deltaF,                 /**< Sampling frequency (Hz). To use non-uniform frequency grid, set deltaF <= 0. */
  REAL8 fRef_In,                      /**< Reference frequency (Hz) */
  LALDict *lalParams                  /**< LAL Dictionary struct */
)
{
  UINT4 status;
  UINT4 debug = PHENOMXPDEBUG;

  status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI,&m2_SI,&chi1x,&chi1y,&chi1z,&chi2x,&chi2y,&chi2z);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

  #if PHENOMXPDEBUG == 1
    printf("fRef_In : %e\n",fRef_In);
    printf("m1_SI   : %e\n",m1_SI);
    printf("m2_SI   : %e\n",m2_SI);
    printf("chi1z   : %e\n",chi1z);
    printf("chi2z   : %e\n",chi2z);
    printf("phiRef  : %e\n",phiRef);
    printf("Prec V. : %d\n\n",XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams));
    printf("Performing sanity checks...\n");
  #endif

  /* Perform initial sanity checks */
  if(*hptilde)        { XLAL_CHECK(NULL != hptilde, XLAL_EFAULT);                                                     }
  if(*hctilde)        { XLAL_CHECK(NULL != hctilde, XLAL_EFAULT);                                                     }
  if(fRef_In  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");                   }
  if(deltaF   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaF must be positive.\n");                                          }
  if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                                              }
  if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                                              }
  if(f_min    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "f_min must be positive.\n");                                           }
  if(f_max    <  0.0) { XLAL_ERROR(XLAL_EDOM, "f_max must be non-negative.\n");                                       }
  if(distance <= 0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");                     }

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
  if(mass_ratio > 20.0  ) { XLAL_PRINT_WARNING("Warning: Extrapolating outside of Numerical Relativity calibration domain. NNLO angles may become pathological at large mass ratios.\n"); }
  if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000.\n"); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1z) > 0.99 || fabs(chi2z) > 0.99) { XLAL_PRINT_WARNING("Warning: Extrapolating to extremal spins, aligned spin model is not trusted.\n"); }

  /* If no reference frequency is given, set it to the starting gravitational wave frequency */
  const REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;
  
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

  /* Spins aligned with the orbital angular momenta */
  const REAL8 chi1L = chi1z;
  const REAL8 chi2L = chi2z;

  #if PHENOMXPDEBUG == 1
      printf("\n\n **** Initializing waveform struct... **** \n\n");
  #endif

  /* Initialize useful powers of LAL_PI - this is used in the code called by IMRPhenomXPGenerateFD */
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.\n");

  /* Initialize IMR PhenomX Waveform struct and check that it initialized correctly */
  IMRPhenomXWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, deltaF, fRef, phiRef, f_min, f_max, distance, inclination, lalParams_aux, debug);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");


  /*
      Create a REAL8 frequency series.
      Use fLow, fHigh, deltaF to compute frequency sequence. Only pass the boundaries (fMin, fMax).
  */
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = pWF->fMin;
  freqs->data[1] = pWF->f_max_prime;


  #if PHENOMXPDEBUG == 1
      printf("\n\n **** Initializing precession struct... **** \n\n");
  #endif


  /* Initialize IMR PhenomX Precession struct and check that it generated successfully */
  IMRPhenomXPrecessionStruct *pPrec;
  pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));

  status = IMRPhenomXGetAndSetPrecessionVariables(
           pWF,
           pPrec,
           m1_SI,
           m2_SI,
           chi1x,
           chi1y,
           chi1z,
           chi2x,
           chi2y,
           chi2z,
           lalParams_aux,
           PHENOMXPDEBUG
         );
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

  #if PHENOMXPDEBUG == 1
      printf("\n\n **** Calling IMRPhenomXPGenerateFD... **** \n\n");
  #endif

  /* We now call the core IMRPhenomXP waveform generator */
  status = IMRPhenomXPGenerateFD(hptilde, hctilde, freqs, pWF, pPrec, lalParams_aux);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXPGenerateFD failed to generate IMRPhenomX waveform.\n");

  #if PHENOMXPDEBUG == 1
      printf("\n\n **** Call to IMRPhenomXPGenerateFD complete. **** \n\n");
  #endif

  /* Resize hptilde & hctilde */
  REAL8 lastfreq;
  if (pWF->f_max_prime < pWF->fMax)
  {
    /*
        The user has requested a higher f_max than Mf = fCut.
        Resize the frequency series to fill with zeros beyond the cutoff frequency.
    */
    lastfreq = pWF->fMax;
  }
  else
  {
     lastfreq = pWF->f_max_prime;
  }
  /* Enforce length to be a power of 2 + 1 */
  size_t n_full = NextPow2(lastfreq / pWF->deltaF) + 1;
  size_t n = (*hptilde)->data->length;

  /* Resize the COMPLEX16 frequency series */
  *hptilde = XLALResizeCOMPLEX16FrequencySeries(*hptilde, 0, n_full);
  XLAL_CHECK (*hptilde, XLAL_ENOMEM, "Failed to resize h_+ COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).\n", n, pWF->fCut, n_full, pWF->fMax );

  /* Resize the COMPLEX16 frequency series */
  *hctilde = XLALResizeCOMPLEX16FrequencySeries(*hctilde, 0, n_full);
  XLAL_CHECK (*hctilde, XLAL_ENOMEM, "Failed to resize h_x COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).\n", n, pWF->fCut, n_full, pWF->fMax );


  /* Destroy structs and arrays */
  LALFree(pWF);
  LALFree(pPrec);
  XLALDestroyREAL8Sequence(freqs);
  XLALDestroyDict(lalParams_aux);

  return XLAL_SUCCESS;
}


/**
 * Compute waveform in LAL format at specified frequencies for the IMRPhenomXP model.
 *
 * XLALSimIMRPhenomXPGenerateFD() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRPhenomXPFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 */
 int XLALSimIMRPhenomXPFrequencySequence(
  COMPLEX16FrequencySeries **hptilde,   /**< [out] Frequency-domain waveform h+  */
  COMPLEX16FrequencySeries **hctilde,   /**< [out] Frequency-domain waveform hx  */
  const REAL8Sequence *freqs,           /**< input Frequency series [Hz]         */
  REAL8 m1_SI,                        /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                        /**< mass of companion 2 (kg) */
  REAL8 chi1x,                        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  const REAL8 distance,               /**< Distance of source (m) */
  const REAL8 inclination,            /**< inclination of source (rad) */
  const REAL8 phiRef,                 /**< Orbital phase (rad) at reference frequency */
  REAL8 fRef_In,                      /**< Reference frequency (Hz) */
  LALDict *lalParams                  /**< LAL Dictionary struct */
)
 {
   UINT4 status = 0;

   status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI,&m2_SI,&chi1x,&chi1y,&chi1z,&chi2x,&chi2y,&chi2z);
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

   XLAL_CHECK(freqs != NULL, XLAL_EFAULT, "Error: XLALSimIMRPhenomXPFrequencySequence *freqs is null.\n");

   /* Perform initial sanity checks */
   if(*hptilde)        { XLAL_CHECK(NULL != hptilde, XLAL_EFAULT);                                                     }
   if(*hctilde)        { XLAL_CHECK(NULL != hctilde, XLAL_EFAULT);                                                     }
   if(fRef_In  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");                   }
   if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                                              }
   if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                                              }
   if(distance <  0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");                     }

   /*
     Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
       - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
       - For 1000 > mass ratio > 20 and spins <= 0.99: print a warning message that we are extrapolating outside of *NR* calibration domain.
       - For mass ratios > 1000: throw a hard error that model is not valid.
       - For spins > 0.99: throw a warning that we are extrapolating the model to extremal
       - At high mass ratios > 20, the NNLO Euler angles may develop known pathologies.
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
   if(mass_ratio > 20.0  ) { XLAL_PRINT_WARNING("Warning: Extrapolating outside of Numerical Relativity calibration domain. NNLO angles may become pathological at large mass ratios.\n"); }
   if(mass_ratio > 1000. && fabs(mass_ratio - 1000) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 1000.\n"); } // The 1e-12 is to avoid rounding errors
   if(sqrt(chi1x*chi1x + chi1y*chi1y + chi1z*chi1z) > 0.99 || sqrt(chi2x*chi2x + chi2y*chi2y + chi2z*chi2z) > 0.99) { XLAL_PRINT_WARNING("Warning: Extrapolating to extremal spins, model is not trusted.\n"); }

   REAL8 chi1L, chi2L;
   chi1L = chi1z;
   chi2L = chi2z;

   /* If fRef is not provided, then set fRef to be the starting GW Frequency */
   const REAL8 fRef = (fRef_In == 0.0) ? freqs->data[0] : fRef_In;

   const REAL8 f_min_In  = freqs->data[0];
   const REAL8 f_max_In  = freqs->data[freqs->length - 1];
   
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

   /*
      Passing deltaF = 0 implies that freqs is a frequency grid with non-uniform spacing.
      The function waveform then start at lowest given frequency.
   */
   status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
   XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.\n");

   /* Initialize IMRPhenomX waveform struct and perform sanity check. */
   IMRPhenomXWaveformStruct *pWF;
   pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
   status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, 0.0, fRef, phiRef, f_min_In, f_max_In, distance, inclination, lalParams_aux, 0);
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

   /* Initialize IMR PhenomX Precession struct and check that it generated successfully */
   IMRPhenomXPrecessionStruct *pPrec;
   pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));

   status = IMRPhenomXGetAndSetPrecessionVariables(
              pWF,
              pPrec,
              m1_SI,
              m2_SI,
              chi1x,
              chi1y,
              chi1z,
              chi2x,
              chi2y,
              chi2z,
              lalParams_aux,
              PHENOMXDEBUG
            );
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetAndSetPrecessionVariables failed.\n");

   /* Now call the core IMRPhenomXP waveform generator */
   status = IMRPhenomXPGenerateFD(hptilde, hctilde, freqs, pWF, pPrec, lalParams_aux);
   XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXASFDCore failed to generate IMRPhenomX waveform.\n");

   LALFree(pPrec);
   LALFree(pWF);
   XLALDestroyDict(lalParams);

   return XLAL_SUCCESS;
 }

 /*
  *  Adapted from IMRPhenomPv2. Avoids clashes with IMRPhenomD/P namespace.
  *
  *  Function to map LAL parameters
  *    - masses,
  *    - 6 spin components
  *    - phiRef at f_ref
  *    - inclination at f_ref
  *
  *  Assumed to be in the source frame (L frame):
  *    - LN vector points in the z direction, i.e. lnhat = (0,0,1)
  *    - Separation vector n is in the x direction
  *    - Spherical angles of the line of sight N are (incl,Pi/2-phiRef))
  *
  *  The IMRPhenomXP intrinsic parameters inherit those from IMRPhenomP:
  *  [chi1_l, chi2_l, chip, thetaJN, alpha0 and phi_aligned]
  *
  *  All input masses and frequencies should be in SI units.
  *
  */
  int XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame(
       REAL8 *chi1L,                     /**< [out] Dimensionless aligned spin on companion 1 */
       REAL8 *chi2L,                     /**< [out] Dimensionless aligned spin on companion 2 */
       REAL8 *chi_p,                     /**< [out] Effective precession parameter: Schmidt, Ohme, Hannam, PRD, 91,024043 (2015) */
       REAL8 *thetaJN,                   /**< [out] Angle between J0 and line of sight (z-direction)            */
       REAL8 *alpha0,                    /**< [out] Initial value of alpha angle (azimuthal precession angle)   */
       REAL8 *phi_aligned,               /**< [out] Initial phase to feed the underlying aligned-spin model     */
       REAL8 *zeta_polarization,         /**< [out] Angle to rotate the polarizations                           */
       const REAL8 m1_SI,                /**< Mass of companion 1 (kg)    */
       const REAL8 m2_SI,                /**< Mass of companion 2 (kg)    */
       const REAL8 f_ref,                /**< Reference GW frequency (Hz) */
       const REAL8 phiRef,               /**< Reference phase (Hz)        */
       const REAL8 incl,                 /**< Inclination : angle between LN and the line of sight */
       const REAL8 chi1x,                /**< Initial value of chi1x: dimensionless spin of BH 1 in L frame    */
       const REAL8 chi1y,                /**< Initial value of chi1y: dimensionless spin of BH 1 in L frame    */
       const REAL8 chi1z,                /**< Initial value of chi1z: dimensionless spin of BH 1 in L frame    */
       const REAL8 chi2x,                /**< Initial value of chi2x: dimensionless spin of BH 2 in L frame    */
       const REAL8 chi2y,                /**< Initial value of chi2y: dimensionless spin of BH 2 in L frame    */
       const REAL8 chi2z,                /**< Initial value of chi2z: dimensionless spin of BH 2 in L frame    */
       LALDict *lalParams                /**< LAL Dictionary */
   )
   {
     /* Perform the usual sanity checks... */
     XLAL_CHECK(chi1L       != NULL, XLAL_EFAULT);
     XLAL_CHECK(chi2L       != NULL, XLAL_EFAULT);
     XLAL_CHECK(chi_p       != NULL, XLAL_EFAULT);
     XLAL_CHECK(thetaJN     != NULL, XLAL_EFAULT);
     XLAL_CHECK(alpha0      != NULL, XLAL_EFAULT);
     XLAL_CHECK(phi_aligned != NULL, XLAL_EFAULT);

     XLAL_CHECK(f_ref > 0, XLAL_EDOM, "Error in XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame: Reference frequency must be positive.\n");
     XLAL_CHECK(m1_SI > 0, XLAL_EDOM, "Error in XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame: m1 must be positive.\n");
     XLAL_CHECK(m2_SI > 0, XLAL_EDOM, "Error in XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame: m2 must be positive.\n");
     XLAL_CHECK(fabs(chi1x*chi1x + chi1y*chi1y + chi1z*chi1z) <= 1.0, XLAL_EDOM, "Error in XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame: |S1/m1^2| must be <= 1.\n");
     XLAL_CHECK(fabs(chi2x*chi2x + chi2y*chi2y + chi2z*chi2z) <= 1.0, XLAL_EDOM, "Error in XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame: |S2/m2^2| must be <= 1.\n");

     /* Get Precession version */
     const int IMRPhenomXPrecVersion     = XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams);

     const REAL8 m1    = m1_SI / LAL_MSUN_SI;    /* Component mass 1 in solar masses */
     const REAL8 m2    = m2_SI / LAL_MSUN_SI;    /* Component mass 2 in solar masses */
     const REAL8 M     = m1+m2;                  /* Total mass in solar masses */
     const REAL8 m1_2  = m1*m1;
     const REAL8 m2_2  = m2*m2;
     const REAL8 eta   = m1 * m2 / (M*M);        /* Symmetric mass-ratio  */
     const REAL8 delta = sqrt(1.0 - 4.0*eta);    /* PN symmetry parameter */

     /*
         From the components in the source frame, we can determine the IMRPhenomXP intrinsic parameters:
         chi1L, chi2L, chi_p and phi_aligned.

         We also compute the spherical angles of J, which are used to transform to the J frame.
     */

     /* Aligned spins */
     *chi1L = chi1z; /* Dimensionless aligned spin on BH 1 */
     *chi2L = chi2z; /* Dimensionless aligned spin on BH 2 */

     /* Magnitude of the spin projections in the orbital plane */
     const REAL8 S1_perp = m1_2*sqrt(chi1x*chi1x + chi1y*chi1y);
     const REAL8 S2_perp = m2_2*sqrt(chi2x*chi2x + chi2y*chi2y);

     /*
         Calculate the single effective precession parameter:
           - Schmidt, Ohme, Hannam, PRD, 91,024043 (2015), arXiv:1408.1810
     */
     const REAL8 A1    = 2 + (3*m2) / (2*m1);
     const REAL8 A2    = 2 + (3*m1) / (2*m2);
     const REAL8 ASp1  = A1*S1_perp;
     const REAL8 ASp2  = A2*S2_perp;
     const REAL8 num   = (ASp2 > ASp1) ? ASp2 : ASp1;
     const REAL8 den   = (m2 > m1) ? A2*m2_2 : A1*m1_2;

     /*  chi_p = max(A1 Sp1, A2 Sp2) / (A_i m_i^2) for i index of larger BH. See Eq. 3.3 and 3.4 of arXiv:1408.1810 */
     *chi_p             = num / den;

     #if DEBUG == 1
     printf("\nCalling: XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame..\n");
     printf("m1    : %e\n",m1);
     printf("m2    : %e\n",m2);
     printf("eta   : %e\n",eta);
     printf("chi_p : %e\n",num/den);
     #endif

     /* Define phase, velocity, etc at reference frequency */
     const REAL8 m_sec = M * LAL_MTSUN_SI;   /* Total mass in seconds */
     const REAL8 piM   = LAL_PI * m_sec;
     const REAL8 v_ref = cbrt(piM * f_ref);

     /* Compute L - we implement different alternative approximations */
     REAL8 L0 = 0.0;

     switch(IMRPhenomXPrecVersion)
     {
       case 101:
       {
         /* Orbital angular momentum L to 2PN, non-spinning. This is as implemented in IMRPhenomPv2 */
         L0 = M*M * XLALSimIMRPhenomXL2PNNS(v_ref, eta);
         break;
       }
       case 102:
       {
         /* Orbital angular momentum L to 3PN with spinning terms. */
         L0 = M*M * XLALSimIMRPhenomXL3PNAS(v_ref, eta, *chi1L, *chi2L, delta);
         break;
       }
       case 103:
       {
         /* 4PN orbital angular momentum with spin terms */
         L0 = M*M * XLALSimIMRPhenomXL4PNAS(v_ref, eta, *chi1L, *chi2L, delta);
         break;
       }
       case 104:
       {
         /* 4PN+ orbital angular momentum */
         L0 = M*M * XLALSimIMRPhenomXL4PNLOSIAS(v_ref, eta, *chi1L, *chi2L, delta);
         break;
       }
       case 220:
       case 221:
       case 222:
       case 223:
       case 224:
       {
         /* 3PN orbital angular momentum  */
         L0 = M*M * XLALSimIMRPhenomXL3PNAS(v_ref, eta, *chi1L, *chi2L, delta);
         break;
       }
       default:
       {
         XLAL_ERROR( XLAL_EINVAL, "Error in XLALSimIMRPhenomXPCalculateModelParametersFromSourceFrame: IMRPhenomXPrecVersion not valid.\n" );
         break;
       }
     }

     /*
       In the following code block we construct the convetions that relate the source frame and the LAL frame.

       A detailed discussion of the conventions can be found in Appendix C and D of arXiv:XXXX.YYYY and https://dcc.ligo.org/LIGO-T1500602
     */
     INT4 convention     = XLALSimInspiralWaveformParamsLookupPhenomXPConvention(lalParams);
     if ( !(convention == 0 || convention == 1 || convention == 5 || convention == 6 || convention == 7) )
     {
       XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPConvention not recognized.\n");
     }

     /* Below, _Sf indicates source frame (L frame) components. Elsewhere we also use _Jf for J frame components */
     const REAL8 J0x_Sf = m1_2*chi1x + m2_2*chi2x;
     const REAL8 J0y_Sf = m1_2*chi1y + m2_2*chi2y;
     const REAL8 J0z_Sf = L0 + (m1_2*chi1z) + (m2_2*chi2z);
     const REAL8 J0     = sqrt(J0x_Sf*J0x_Sf + J0y_Sf*J0y_Sf + J0z_Sf*J0z_Sf);

     /* Compute thetaJ_Sf, the angle between J0 and LN (z-direction) */
     REAL8 thetaJ_Sf;
     if (J0 < 1e-10)
     {
       XLAL_PRINT_WARNING("Warning: |J0| < 1e-10. Setting thetaJ_Sf = 0.\n");
       thetaJ_Sf = 0;
     }
     else
     {
       thetaJ_Sf = acos(J0z_Sf / J0);
     }

     double phiRefIn = phiRef;

     /* The azimuthal angle of J0 in the source frame */
     REAL8 phiJ_Sf;
     if (fabs(J0x_Sf) < MAX_TOL_ATAN && fabs(J0y_Sf) < MAX_TOL_ATAN)
     {
       // Then we are ~ aligned spin
       switch(convention)
       {
         case 0:
         case 5:
         {
           phiJ_Sf = LAL_PI/2. - phiRefIn;
           break;
         }
         case 1:
         case 6:
         case 7:
         {
           phiJ_Sf = 0;
           break;
         }
         default:
         {
           XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPConvention not recognized.\n");
           break;
         }
       }
     }
     else
     {
       phiJ_Sf = atan2(J0y_Sf, J0x_Sf);
     }

     switch(convention)
     {
       case 0:
       {
         *phi_aligned = - phiJ_Sf;
         break;
       }
       case 1:
       {
         *phi_aligned = 0.0;
         break;
       }
       case 5:
       case 6:
       case 7:
       {
         *phi_aligned = phiRef;
         break;
       }
       default:
       {
         XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPConvention not recognized.\n");
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
         This is done by multiplying by: RZ[kappa].RY[-thetaJ].RZ[-phiJ]

         Note that kappa is determined by rotating N with RY[-thetaJ].RZ[-phiJ], which brings J to the z-axis, and
         taking the opposite of the azimuthal angle of the rotated N.
     */

     /*
         First determine kappa in the source frame.

         The components of N are given in Eq (35c) of T1500606-v6:
     */
     REAL8 Nx_Sf = sin(incl)*cos(LAL_PI/2. - phiRef);
     REAL8 Ny_Sf = sin(incl)*sin(LAL_PI/2. - phiRef);
     REAL8 Nz_Sf = cos(incl);

     REAL8 tmp_x = Nx_Sf;
     REAL8 tmp_y = Ny_Sf;
     REAL8 tmp_z = Nz_Sf;

     IMRPhenomX_rotate_z(-phiJ_Sf,   &tmp_x, &tmp_y, &tmp_z);
     IMRPhenomX_rotate_y(-thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);

     REAL8 kappa;
     kappa = XLALSimIMRPhenomXatan2tol(tmp_y, tmp_x, MAX_TOL_ATAN);

     /* Now determine alpha0 by rotating LN. In the source frame, LN = {0,0,1} */
     tmp_x = 0.0;
     tmp_y = 0.0;
     tmp_z = 1.0;

     IMRPhenomX_rotate_z(-phiJ_Sf,   &tmp_x, &tmp_y, &tmp_z);
     IMRPhenomX_rotate_y(-thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);
     IMRPhenomX_rotate_z(-kappa,     &tmp_x, &tmp_y, &tmp_z);

     if (fabs(tmp_x) < MAX_TOL_ATAN && fabs(tmp_y) < MAX_TOL_ATAN)
     {
       // This is the aligned spin case
       switch(convention)
       {
         case 0:
         case 5:
         {
           *alpha0 = LAL_PI;
           break;
         }
         case 1:
         case 6:
         case 7:
         {
           *alpha0 = LAL_PI - kappa;
           break;
         }
         default:
         {
           XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPConvention not recognized.\n");
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
           *alpha0 = atan2(tmp_y,tmp_x);
           break;
         }
         case 1:
         case 6:
         case 7:
         {
           *alpha0 = LAL_PI - kappa;
           break;
         }
         default:
         {
           XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPConvention not recognized.\n");
           break;
         }
       }
     }

     REAL8 Nx_Jf = 0.0;
     REAL8 Nz_Jf = 0.0;
     switch(convention)
     {
       case 0:
       case 5:
       {
         // Finally we determine thetaJ, by rotating N
         tmp_x = Nx_Sf;
         tmp_y = Ny_Sf;
         tmp_z = Nz_Sf;
         IMRPhenomX_rotate_z(-phiJ_Sf,     &tmp_x, &tmp_y, &tmp_z);
         IMRPhenomX_rotate_y(-thetaJ_Sf,   &tmp_x, &tmp_y, &tmp_z);
         IMRPhenomX_rotate_z(kappa,        &tmp_x, &tmp_y, &tmp_z);

         // Do not need the y-component
         Nx_Jf = tmp_x;
         Nz_Jf = tmp_z;
         *thetaJN    = acos(Nz_Jf);
         break;
       }
       case 1:
       case 6:
       case 7:
       {
          REAL8 J0dotN = (J0x_Sf * Nx_Sf) + (J0y_Sf * Ny_Sf) + (J0z_Sf * Nz_Sf);
          *thetaJN     = acos( J0dotN / J0 );
          Nz_Jf        = cos(*thetaJN);
          Nx_Jf        = sin(*thetaJN);
          break;
       }
       default:
       {
         XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPConvention not recognized.\n");
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
     REAL8 Xx_Sf = -cos(incl)*sin(phiRef);
     REAL8 Xy_Sf = -cos(incl)*cos(phiRef);
     REAL8 Xz_Sf = sin(incl);

     tmp_x = Xx_Sf;
     tmp_y = Xy_Sf;
     tmp_z = Xz_Sf;
     IMRPhenomX_rotate_z(-phiJ_Sf,   &tmp_x, &tmp_y, &tmp_z);
     IMRPhenomX_rotate_y(-thetaJ_Sf, &tmp_x, &tmp_y, &tmp_z);
     IMRPhenomX_rotate_z(kappa,      &tmp_x, &tmp_y, &tmp_z);

     REAL8 PArunx_Jf, PAruny_Jf, PArunz_Jf;
     REAL8 QArunx_Jf, QAruny_Jf, QArunz_Jf;

     switch(convention)
     {
       case 0:
       case 5:
       {
         /*
             The components tmp_i are now the components of X in the J frame.

             We now need the polar angle of this vector in the P, Q basis of Arun et al:

                 P = (N x J) / |NxJ|

             Note, that we put N in the (pos x)z half plane of the J frame
         */

         /* Get the polar angle of X vector in the J frame using the P, Q basis of Arun et al */
         PArunx_Jf = 0.0;
         PAruny_Jf = -1.0;
         PArunz_Jf = 0.0;

         // Q = N x P
         QArunx_Jf = Nz_Jf;
         QAruny_Jf = 0.0;
         QArunz_Jf = -Nx_Jf;
         break;
       }
       case 1:
       case 6:
       case 7:
       {
         PArunx_Jf = +Nz_Jf;
         PAruny_Jf = 0.0;
         PArunz_Jf = -Nx_Jf;

         QArunx_Jf = 0.0;
         QAruny_Jf = 1.0;
         QArunz_Jf = 0.0;
         break;
       }
       default:
       {
         XLAL_ERROR(XLAL_EINVAL,"Error: IMRPhenomXPConvention not recognized.\n");
         break;
       }
     }

     REAL8 XdotPArun     = tmp_x*PArunx_Jf + tmp_y*PAruny_Jf + tmp_z*PArunz_Jf;
     REAL8 XdotQArun     = tmp_x*QArunx_Jf + tmp_y*QAruny_Jf + tmp_z*QArunz_Jf;

     /* Return polarization angle */
     *zeta_polarization  = atan2(XdotQArun , XdotPArun);

     return XLAL_SUCCESS;
 }


 /*
  *  Prototype wrapper function:

  *  Driver routine to calculate the MSA Euler angles in the frequency domain.
  *
  *  All input parameters should be in SI units.
  *
  *  Returns \f$\phi_z\f$, \f$\zeta\f$ and \f$\cos \theta_L \f$
  */
 int XLALSimIMRPhenomXPMSAAngles(
  REAL8Sequence *phiz_of_f,        /**< [out] The azimuthal angle of L around J */
  REAL8Sequence *zeta_of_f,        /**< [out] The third Euler angle describing L with respect to J. Fixed by minmal rotation condition. */
  REAL8Sequence *costhetaL_of_f,   /**< [out]  Cosine of polar angle between L and J */
  const REAL8Sequence *freqs,      /**< [out] Frequency series [Hz]         */
  REAL8 m1_SI,               /**< mass of companion 1 (kg) */
  REAL8 m2_SI,               /**< mass of companion 2 (kg) */
  REAL8 chi1x,               /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,               /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,               /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,               /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,               /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,               /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 fRef_In,             /**< Reference frequency (Hz) */
  LALDict *lalParams               /**< LAL Dictionary struct */
)
{
   /*
      Passing deltaF = 0 implies that freqs is a frequency grid with non-uniform spacing.
      The function waveform then start at lowest given frequency.
   */

   UINT4 status = 0;

   status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI,&m2_SI,&chi1x,&chi1y,&chi1z,&chi2x,&chi2y,&chi2z);
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

   /* Perform initial sanity checks */
   if(fRef_In  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");                   }
   if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                                              }
   if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                                              }

   REAL8 chi1L, chi2L;
   chi1L = chi1z;
   chi2L = chi2z;

   /* If fRef is not provided, then set fRef to be the starting GW Frequency */
   const REAL8 fRef = (fRef_In == 0.0) ? freqs->data[0] : fRef_In;

   const REAL8 f_min_In  = freqs->data[0];
   const REAL8 f_max_In  = freqs->data[freqs->length - 1];
   
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

   /* Initialize IMRPhenomX waveform struct and perform sanity check. */
   IMRPhenomXWaveformStruct *pWF;
   pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
   status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, 0.0, fRef, 0.0, f_min_In, f_max_In, 1.0, 0.0, lalParams_aux, 0);
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");


   /* Initialize IMR PhenomX Precession struct and check that it generated successfully */
   IMRPhenomXPrecessionStruct *pPrec;
   pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));

   const int pflag     = XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams_aux);

   if(pflag != 220 && pflag != 221 && pflag != 222 && pflag != 223 && pflag != 224)
   {
     XLAL_ERROR(XLAL_EDOM, "Error: MSA system currently only supported for IMRPhenomXPrecVersion 220, 221, 222, 223 or 224.\n");
   }

   status = IMRPhenomXGetAndSetPrecessionVariables(
              pWF,
              pPrec,
              m1_SI,
              m2_SI,
              chi1x,
              chi1y,
              chi1z,
              chi2x,
              chi2y,
              chi2z,
              lalParams_aux,
              PHENOMXDEBUG
            );
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetAndSetPrecessionVariables failed.\n");

   REAL8 v        = 0.0;
   vector vangles = {0.,0.,0.};

   for(UINT4 i = 0; i < (*freqs).length; i++)
   {
     // Input list of *orbital* frequencies not *gravitational-wave* frequencies*
     // v     = cbrt( ((*freqs).data[i]) * pPrec->twopiGM );

     // Input list of *gravitational-wave* frequencies not *orbital* frequencies*
     v       = cbrt( ((*freqs).data[i]) * pPrec->piGM );
     vangles = IMRPhenomX_Return_phi_zeta_costhetaL_MSA(v,pWF,pPrec);

     (*phiz_of_f).data[i]      = vangles.x - pPrec->alpha_offset;
     (*zeta_of_f).data[i]      = vangles.y - pPrec->epsilon_offset;
     (*costhetaL_of_f).data[i] = vangles.z;
   }

   LALFree(pPrec);
   LALFree(pWF);
   XLALDestroyDict(lalParams);

   return XLAL_SUCCESS;
 }

 /*
  *  Prototype wrapper function:

  *  Driver routine to calculate the NNLO PN Euler angles in the frequency domain.
  *
  *  All input parameters should be in SI units.
  *
  *  Returns \f$\alpha\f$, \f$\cos \beta\f$ and \f$\gamma\f$
  */
 int XLALSimIMRPhenomXPPNAngles(
  REAL8Sequence *alpha_of_f,                /**< [out] Azimuthal angle of L w.r.t J */
  REAL8Sequence *gamma_of_f,                /**< [out] Third Euler angle describing L w.r.t J, fixed by minimal rotation condition */
  REAL8Sequence *cosbeta_of_f,              /**< [out] Cosine of polar angle between L and J */
  const REAL8Sequence *freqs,               /**< [out] Frequency series [Hz]         */
  REAL8 m1_SI,                              /**< mass of companion 1 (kg) */
  REAL8 m2_SI,                              /**< mass of companion 2 (kg) */
  REAL8 chi1x,                              /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1y,                              /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi1z,                              /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2x,                              /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2y,                              /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 chi2z,                              /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
  REAL8 fRef_In,                            /**< Reference frequency (Hz) */
  LALDict *lalParams                        /**< LAL Dictionary struct */
)
{
   /*
      Passing deltaF = 0 implies that freqs is a frequency grid with non-uniform spacing.
      The function waveform then start at lowest given frequency.
   */

   UINT4 status = 0;

   status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI,&m2_SI,&chi1x,&chi1y,&chi1z,&chi2x,&chi2y,&chi2z);
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: XLALIMRPhenomXPCheckMassesAndSpins failed.\n");

   #if PHENOMXPDEBUG == 1
    printf("\nm1  : %.6f\n",m1_SI);
    printf("m2    : %.6f\n",m2_SI);
    printf("chi1x : %.6f\n",chi1x);
    printf("chi1y : %.6f\n",chi1y);
    printf("chi1z : %.6f\n",chi1z);
    printf("chi2x : %.6f\n",chi2x);
    printf("chi2y : %.6f\n",chi2y);
    printf("chi2z : %.6f\n\n",chi2z);
   #endif

   /* Perform initial sanity checks */
   if(fRef_In  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");                   }
   if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                                              }
   if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                                              }

   REAL8 chi1L, chi2L;
   chi1L = chi1z;
   chi2L = chi2z;

   /* If fRef is not provided, then set fRef to be the starting GW Frequency */
   const REAL8 fRef = (fRef_In == 0.0) ? freqs->data[0] : fRef_In;

   const REAL8 f_min_In  = freqs->data[0];
   const REAL8 f_max_In  = freqs->data[freqs->length - 1];
   
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

   /* Initialize IMRPhenomX waveform struct and perform sanity check. */
   IMRPhenomXWaveformStruct *pWF;
   pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
   status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, 0.0, fRef, 0.0, f_min_In, f_max_In, 1.0, 0.0, lalParams_aux, 0);
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

   /* Initialize IMR PhenomX Precession struct and check that it generated successfully */
   IMRPhenomXPrecessionStruct *pPrec;
   pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));

   status = IMRPhenomXGetAndSetPrecessionVariables(
              pWF,
              pPrec,
              m1_SI,
              m2_SI,
              chi1x,
              chi1y,
              chi1z,
              chi2x,
              chi2y,
              chi2z,
              lalParams_aux,
              PHENOMXDEBUG
            );
   XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetAndSetPrecessionVariables failed.\n");

   /* Spin variable used to calculate \cos \beta */
   REAL8 s, s2;

   /* Orbital velocity and orbital frequency */
   REAL8 v, omega, logomega, omega_cbrt, omega_cbrt2, f;

   /* PN Orbital angular momenta */
   REAL8 L = 0.0;

   for(UINT4 i = 0; i < (*freqs).length; i++)
   {
     f           = ((*freqs).data[i]);

     /* Orbital frequency and velocity */
     omega       = f * pPrec->piGM;
     logomega    = log(omega);
     omega_cbrt  = cbrt(omega);
     omega_cbrt2 = omega_cbrt * omega_cbrt;
     v           = omega_cbrt;

     L = XLALSimIMRPhenomXLPNAnsatz(v, pWF->eta/v, pPrec->L0, pPrec->L1, pPrec->L2, pPrec->L3, pPrec->L4, pPrec->L5, pPrec->L6, pPrec->L7, pPrec->L8, pPrec->L8L);

     (*alpha_of_f).data[i]      =  IMRPhenomX_PN_Euler_alpha_NNLO(pPrec,omega,omega_cbrt2,omega_cbrt,logomega);

     /* \gamma = - \epsilon */
     (*gamma_of_f).data[i]      = -IMRPhenomX_PN_Euler_epsilon_NNLO(pPrec,omega,omega_cbrt2,omega_cbrt,logomega);

     s        = pPrec->Sperp / (L + pPrec->SL);
     s2       = s*s;

     (*cosbeta_of_f).data[i]    = 1.0 / sqrt(1.0 + s2);
   }

   LALFree(pPrec);
   LALFree(pWF);
   XLALDestroyDict(lalParams_aux);

   return XLAL_SUCCESS;
 }

 /** @} */
 /** @} */

 /* *********************************************************************************
  *
  * The following private function generates an IMRPhenomXP frequency-domain waveform
  *   - Precessing spins
  *   - Only the 22 mode in the co-precessing frame
  *   - Physical parameters are passed via the waveform struct
  * *********************************************************************************
  */
int IMRPhenomXPGenerateFD(
  COMPLEX16FrequencySeries **hptilde,       /**< [out] FD waveform           */
  COMPLEX16FrequencySeries **hctilde,       /**< [out] FD waveform           */
  const REAL8Sequence *freqs_In,            /**< Input frequency grid        */
  IMRPhenomXWaveformStruct *pWF,            /**< IMRPhenomX Waveform Struct  */
  IMRPhenomXPrecessionStruct *pPrec,        /**< IMRPhenomXP Waveform Struct */
  LALDict *lalParams                        /**< LAL Dictionary Structure    */
)
{
  #if PHENOMXPDEBUG == 1
    printf("\n **** Now in IMRPhenomXPGenerateFD... **** \n");
  #endif

  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  /* Initialize useful powers of LAL_PI */
  int status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.\n");

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
    npts = (size_t) (f_max / pWF->deltaF) + 1;

    /* Debug information */
    #if PHENOMXPDEBUG == 1
      printf("npts     = %zu\n",npts);
      printf("fMin     = %.4f\n",f_min);
      printf("fMax     = %.4f\n",f_max);
      printf("dF       = %.4f\n",pWF->deltaF);
    #endif

    /* Coalescence time is fixed to t=0, shift by overall length in time.  */
    XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / pWF->deltaF), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.\n", pWF->deltaF);

    /* Initialize the htilde frequency series */
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,npts);
    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,npts);

    /* Check that frequency series generated okay */
    XLAL_CHECK(*hptilde,XLAL_ENOMEM,"Failed to allocate COMPLEX16FrequencySeries h_+ of length %zu for f_max = %f, deltaF = %g.\n",npts,f_max,pWF->deltaF);
    XLAL_CHECK(*hctilde,XLAL_ENOMEM,"Failed to allocate COMPLEX16FrequencySeries h_x of length %zu for f_max = %f, deltaF = %g.\n",npts,f_max,pWF->deltaF);

    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t) (f_min / pWF->deltaF);
    size_t iStop  = (size_t) (f_max / pWF->deltaF) + 1;

    XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
          "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.\n", iStart, iStop, npts);

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
    *hptilde  = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &ligotimegps_zero, f_min, pWF->deltaF, &lalStrainUnit, npts);
    *hctilde  = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &ligotimegps_zero, f_min, pWF->deltaF, &lalStrainUnit, npts);

    XLAL_CHECK (*hptilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries h_+ of length %zu from sequence.\n", npts);
    XLAL_CHECK (*hctilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries h_x of length %zu from sequence.\n", npts);

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

  memset((*hptilde)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*hptilde)->sampleUnits), &((*hptilde)->sampleUnits), &lalSecondUnit);

  memset((*hctilde)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*hctilde)->sampleUnits), &((*hctilde)->sampleUnits), &lalSecondUnit);

  /* Check if LAL dictionary exists. If not, create a LAL dictionary. */
  INT4 lalParams_In = 0;
  if(lalParams == NULL)
  {
    lalParams_In = 1;
    lalParams    = XLALCreateDict();
  }

  #if PHENOMXPDEBUG == 1
    printf("\n\n **** Initializing amplitude struct... **** \n\n");
  #endif

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

  /* Allocate and initialize the PhenomX 22 amplitude coefficients struct */
  IMRPhenomXAmpCoefficients *pAmp22;
  pAmp22 = XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
  status = IMRPhenomXGetAmplitudeCoefficients(pWF,pAmp22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetAmplitudeCoefficients failed.\n");

  #if PHENOMXPDEBUG == 1
    printf("\n\n **** Amplitude struct initialized. **** \n\n");
    printf("\n\n **** Initializing phase struct... **** \n\n");
  #endif

  /* Allocate and initialize the PhenomX 22 phase coefficients struct */
  IMRPhenomXPhaseCoefficients *pPhase22;
  pPhase22 = XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
  status   = IMRPhenomXGetPhaseCoefficients(pWF,pPhase22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetPhaseCoefficients failed.\n");

  #if PHENOMXPDEBUG == 1
    printf("\n\n **** Phase struct initialized. **** \n\n");
  #endif

  /* Initialize a struct containing useful powers of Mf at fRef */
  IMRPhenomX_UsefulPowers powers_of_MfRef;
  status = IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF->MfRef);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "IMRPhenomX_Initialize_Powers failed for MfRef.\n");

  /*
      Time shift so peak amplitude is near t ~ 0.
  */

  /* Linear time and phase shifts so that model peaks near t ~ 0                        */
  /* lina provides an additional phase shift, which we currently don't use, but we keep */
  /* the variable in case we change our mind about phase alignement in the future       */
  REAL8 lina = 0;

  /* Get phase connection coefficients = phase offsets between calibration regions to achieve C^1 phase */

  IMRPhenomX_Phase_22_ConnectionCoefficients(pWF,pPhase22);

  /* Apply time shift, see IMRPhenomX_TimeShift_22 for details of current implementation */
  double linb = IMRPhenomX_TimeShift_22(pPhase22, pWF);

  /* Inverse of the symmetric mass ratio */
  REAL8 inveta    = (1.0 / pWF->eta);

  /* Construct reference phase, see discussion in Appendix A of arXiv:2001.11412 */
  pWF->phifRef    = -(inveta * IMRPhenomX_Phase_22(pWF->MfRef, &powers_of_MfRef, pPhase22, pWF) + linb*pWF->MfRef + lina) + 2.0*pWF->phi0 + LAL_PI_4;

  /*
      Here we declare explicit REAL8 variables for main loop in order to avoid numerous
      pointer calls.
  */
  REAL8 C1IM              = pPhase22->C1Int;
  REAL8 C2IM              = pPhase22->C2Int;
  REAL8 C1RD              = pPhase22->C1MRD;
  REAL8 C2RD              = pPhase22->C2MRD;

  REAL8 fPhaseIN          = pPhase22->fPhaseMatchIN;
  REAL8 fPhaseIM          = pPhase22->fPhaseMatchIM;
  REAL8 fAmpIN            = pAmp22->fAmpMatchIN;
  REAL8 fAmpIM            = pAmp22->fAmpRDMin;

  #if PHENOMXPDEBUG == 1
    printf("\n\n **** Phase struct initialized. **** \n\n");
    printf("C1IM     = %.4f\n",C1IM);
    printf("C2IM     = %.4f\n",C2IM);
    printf("C1RD     = %.4f\n",C1RD);
    printf("C2RD     = %.4f\n",C2RD);
    printf("fIN      = %.4f\n",fPhaseIN);
    printf("fIM      = %.4f\n",fPhaseIM);
    printf("thetaJN  = %.4f\n",pPrec->thetaJN);
  #endif

  /* amplitude definition as in XAS */
  REAL8 Amp0          = pWF->amp0 * pWF->ampNorm;

  /* initial_status used to track  */
  UINT4 initial_status = XLAL_SUCCESS;

  /* Approximate frequency of the peak */
  REAL8 fmax         = pAmp22->fAmpRDMin;

  /* Initialize a struct containing useful powers of fmax */
  IMRPhenomX_UsefulPowers powers_of_fmax;
  initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_fmax,fmax);
  if(initial_status != XLAL_SUCCESS)
  {
    status = initial_status;
    XLALPrintError("IMRPhenomX_Initialize_Powers failed for fmax, initial_status=%d",initial_status);
  }

  #if DEBUG == 1
  // Save anles into a file
  FILE *fileangle;
  char fileSpec[40];
  sprintf(fileSpec, "angles_XP.dat");

  printf("\nOutput angle file: %s\r\n", fileSpec);
  fileangle = fopen(fileSpec,"w");

  fprintf(fileangle,"# q = %.16e m1 = %.16e m2 = %.16e chi1 = %.16e chi2 = %.16e  Mtot = %.16e distance = %.16e\n", pWF->q, pWF->m1, pWF->m2, pWF->chi1L, pWF->chi2L,  pWF->Mtot, pWF->distance/LAL_PC_SI/1e6);
  fprintf(fileangle,"#fHz   cexp_i_alpha(re im)   cexp_i_epsilon(re im)    cexp_i_betah(re im)\n");

  fclose(fileangle);
  #endif

  /* Now loop over frequencies to generate waveform:  h(f) = A(f) * Exp[I phi(f)] */
  #pragma omp parallel for
  for (UINT4 idx = 0; idx < freqs->length; idx++)
  {
    double Mf    = pWF->M_sec * freqs->data[idx];
    UINT4 jdx    = idx  + offset;

    COMPLEX16 hcoprec     = 0.0;  /* Co-precessing waveform */
    COMPLEX16 hplus       = 0.0;  /* h_+ */
    COMPLEX16 hcross      = 0.0;  /* h_x */

    /* Initialize a struct containing useful powers of Mf */
    IMRPhenomX_UsefulPowers powers_of_Mf;
    initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
    if(initial_status != XLAL_SUCCESS)
    {
      status = initial_status;
      XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d\n",initial_status);
    }
    else
    {
      /* Generate amplitude and phase at Mf */

      // initialize amplitude and phase
      REAL8 amp = 0.0;
      REAL8 phi = 0.0;

      /* Here we explicitly call the functions which treat the non-overlapping */
      /* inspiral, intermediate and ringdown frequency regions for the non-precessing waveform. */

      /* Get phase */
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
      /* Scale phase by 1/eta and apply phase and time shifts */
      phi  *= inveta;
      phi  += linb*Mf + lina + pWF->phifRef;

      /* Get amplitude */
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

      /* Waveform in co-precessing frame: h(f) = A(f) * Exp[I phi(f)] */
      hcoprec = Amp0 * powers_of_Mf.m_seven_sixths * amp * cexp(I * phi);

      /* Transform modes from co-precessing frame to inertial frame */
      IMRPhenomXPTwistUp22(Mf,hcoprec,pWF,pPrec,&hplus,&hcross);

      /* Populate h_+ and h_x */
      ((*hptilde)->data->data)[jdx] = hplus;
      ((*hctilde)->data->data)[jdx] = hcross;
    }
  }

  /*
      Loop over h+ and hx and rotate waveform by 2 \zeta.

      Note that we do this here rather than in LALSimInpsiral.c as
      \zeta_polarization is part of the pPrec struct.

  */
  COMPLEX16 PhPpolp, PhPpolc;
  REAL8 cosPolFac, sinPolFac;

  /* If \zeta is non-zero, we need to rotate h+ and hx */
  if(fabs(pPrec->zeta_polarization) > 0.0)
  {
    cosPolFac = cos(2.0 * pPrec->zeta_polarization);
    sinPolFac = sin(2.0 * pPrec->zeta_polarization);

    for (UINT4 i = 0; i < (*hptilde)->data->length; i++)
    {
        PhPpolp = (*hptilde)->data->data[i];
        PhPpolc = (*hctilde)->data->data[i];

        (*hptilde)->data->data[i] = (cosPolFac * PhPpolp) + (sinPolFac * PhPpolc);
        (*hctilde)->data->data[i] = (cosPolFac * PhPpolc) - (sinPolFac * PhPpolp);
    }
  }

  /* Free allocated memory */
  LALFree(pAmp22);
  LALFree(pPhase22);
  XLALDestroyREAL8Sequence(freqs);

  if(lalParams_In == 1)
  {
    XLALDestroyDict(lalParams);
  }

  return status;
}
