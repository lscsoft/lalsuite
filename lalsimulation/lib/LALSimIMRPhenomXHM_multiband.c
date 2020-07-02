/*
 * Copyright (C) 2019 Cecilio García Quirós, Sascha Husa
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
 /* This code applies the multibanding technique described in arXiv:2001.10897 to the model IMRPhenomXHM described in arXiv:2001.10914. */

#include <complex.h>

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_spline.h>

#include "LALSimIMRPhenomXHM_multiband.h"


/*****************************************************/
/*                                                   */
/*                MULTIBANDING GRIDS                 */
/*                                                   */
/*****************************************************/

/* Build IMRPhenomXMultiBandingGridStruct object, equally-spaced grid */
IMRPhenomXMultiBandingGridStruct XLALSimIMRPhenomXGridComp(
  REAL8 fSTART,  /**< Starting frequency of the uniform bin **/
  REAL8 fEND,    /**< Ending frequency of the uniform bin **/
  REAL8 mydf     /**< Frequency spacing of the bin **/
) {

  IMRPhenomXMultiBandingGridStruct myGrid;

  REAL8 Deltaf, maxfreq;
  INT4 nCoarseIntervals;

  Deltaf = fEND - fSTART;
  nCoarseIntervals = (int)ceil(Deltaf/mydf);
  maxfreq = fSTART + mydf * nCoarseIntervals;

  /* Model Version Parameters */
  myGrid.debug = DEBUG;
  myGrid.nIntervals = nCoarseIntervals;
  myGrid.xStart = fSTART;
  myGrid.xEndRequested = fEND;
  myGrid.xEndFrom_xStart_dx = maxfreq;
  myGrid.xMax = maxfreq;
  myGrid.deltax = mydf;
  myGrid.Length = nCoarseIntervals + 1;

  #if DEBUG == 1
  printf("\nGridComp: fSTART = %.16f",fSTART);
  printf("\nGridComp: xMax = %.16f",maxfreq);
  printf("\nGridComp: Length = %i",nCoarseIntervals+1);
  printf("\nGridComp: x0 = %e : x1 = %e : length = %d : mydf = %e\n", myGrid.xStart, myGrid.xMax, myGrid.Length, mydf);
  #endif

  return myGrid;
}

/* Build non equally-spaced coarse frequency grid */
INT4 XLALSimIMRPhenomXMultibandingGrid(
  REAL8 fstartIn,                             /**< Minimun frequency in NR unit s**/
  REAL8 fend,                                 /**< End of inspiral frequency bins **/
  REAL8 MfLorentzianEnd,                      /**< Determines the last frequency bin **/
  REAL8 Mfmax,                                /**< Maximun frequency in NR units **/
  REAL8 evaldMf,                              /**< Spacing of the uniform frequency grid (NR units) **/
  REAL8 dfpower,                              /**< decaying frequency power to estimate frequency spacing **/
  REAL8 dfcoefficient,                        /**< multiplying factor to the estimate of the frequency spacing **/
  IMRPhenomXMultiBandingGridStruct *allGrids, /**<[out] list of non-uniform frequency bins**/
  REAL8 dfmerger,                             /**<[out] Spacing merger bin**/
  REAL8 dfringdown                            /**<[out] Spacing ringdown bin**/
){
  // the df power law: df = dfcoefficient f^dfpower

  INT4 index, intdfRatio;
  // Variables to control the subgrids to be computed
  INT4 preComputeFirstGrid, nMergerGrid, nRingdownGrid, nDerefineInspiralGrids;

  /* Object to store information about one subgrid */
  IMRPhenomXMultiBandingGridStruct coarseGrid;
  coarseGrid.xMax = fstartIn;

  REAL8 df0, FrequencyFactor=1, nextfSTART, origLogFreqFact;
  REAL8 fSTART, mydf = evaldMf, fEND, fEndGrid0,fStartInspDerefinement, fEndInsp;

  /* Number of fine freq points between two coarse freq points */
  /* The numerator of this quantity corresponds to eqs. 2.8, 2.9 in arXiv:2001.10897. */
  REAL8 const dfRatio = dfcoefficient * pow(fstartIn, dfpower)/evaldMf;

  REAL8 df0original = dfRatio*evaldMf;
  coarseGrid.deltax = df0original;

  #if DEBUG == 1
  printf("\ndfcoefficient = %.16e", dfcoefficient);
  printf("\nfstartIn = %.16e", fstartIn);
  printf("\ndfpower = %.16e", dfpower);
  printf("\nevaldMf = %.16e", evaldMf);
  printf("\ndfRatio = %.16e\n", dfRatio);
  #endif

  if (dfRatio < 1.0) {
    /* User asks for a df that is coarser than one predicted by the multibanding criteria, so we take the users's df  */
    #if DEBUG == 1
    printf("\n****Adjusting frequency factors!****\n");
    #endif
    preComputeFirstGrid = 1;
    intdfRatio = 1;
    df0 = evaldMf;
    fEndGrid0 = pow(evaldMf/dfcoefficient,1./dfpower);
    fStartInspDerefinement = fEndGrid0 + 2 * df0;
    #if DEBUG == 1
    printf("\nevaldMf = %.6f\n", fEndGrid0);
    printf("\ndfcoefficient = %.6f\n", dfcoefficient);
    printf("\ndfpower = %.6f\n", dfpower);
    printf("\nfEndGrid0 = %.6f\n", fEndGrid0);
    #endif
  }
  else {
    #if DEBUG == 1
    printf("proceed without preComputeFirstGrid!");
    #endif
    preComputeFirstGrid = 0;
    intdfRatio = (int)floor(dfRatio);
    fStartInspDerefinement = fstartIn;
  }

  df0 = evaldMf * intdfRatio;

  #if DEBUG == 1
  printf("\nintdfRatio = %d\n", intdfRatio);
  printf("\nfStartInspDerefinement  = %e\n", fStartInspDerefinement );
  printf("\ndf0, df0original, evaldMf  = %.6e %.6e %.6e\n", df0, df0original, evaldMf );
  #endif


  if (fStartInspDerefinement >= fend) {
    fEndInsp = fStartInspDerefinement;
    nDerefineInspiralGrids = 0;
  }
  else
  {
    // Compute the number of inspiral subgrids needed and the ending frequency
    FrequencyFactor = pow(2., (1./dfpower));
    origLogFreqFact = logbase(FrequencyFactor, fend/fStartInspDerefinement);   // This is eq 2.40 in arXiv:2001.10897

    nDerefineInspiralGrids=(int)(ceil(origLogFreqFact));
    fEndInsp = fStartInspDerefinement * pow(FrequencyFactor,nDerefineInspiralGrids);  // estimate, could change due to boundary effects? maybe need one grid less

    #if DEBUG == 1
    printf("FrequencyFactor, fend, fStartInspDerefinement, origLogFreqFact: %e : %e : %e :%e\n", FrequencyFactor, fend, fStartInspDerefinement, origLogFreqFact);
    printf("df0/evaldMf = %e\n", df0/evaldMf);
    printf("Factor in frequency between adjacent inspiral grids = %e\n",  FrequencyFactor);
    printf("Number of subgrids required = %d :  unrounded =  : %e\n",  nDerefineInspiralGrids, origLogFreqFact);
    #endif
  }

  #if DEBUG == 1
  printf("\nMfMECO = %.16e\n fEndInsp = %.16e\n MfLorentzianEnd = %.16e\n Mfmax = %.16e\n", fend, fEndInsp, MfLorentzianEnd, Mfmax);
  #endif

  /* Adjust transition frequencies for special cases. */
  if (fEndInsp + evaldMf >= MfLorentzianEnd) {
    nMergerGrid   = 0;
    if (fEndInsp + evaldMf >= Mfmax) {
      nRingdownGrid = 0;
    } else {
      nRingdownGrid = 1;
    }
  } else {
    nMergerGrid   = 1;
    nRingdownGrid = 1;
    if(MfLorentzianEnd > Mfmax){
      nRingdownGrid = 0;
    }
  }


  #if DEBUG == 1
  printf("nMergerGrid = %d\n", nMergerGrid);
  printf("nRingdownGrid = %d\n", nRingdownGrid);
  printf("fStartInspDerefinement = %e\n",  fStartInspDerefinement);
  printf("fEndInsp = %e\n",  fEndInsp);
  #endif


  // Precompute the first grid if needed
  if (preComputeFirstGrid > 0) {

    mydf = evaldMf;
    coarseGrid = XLALSimIMRPhenomXGridComp(fstartIn,fEndGrid0,mydf);

    allGrids[0] = coarseGrid;
    allGrids[0].intdfRatio = 1;
    #if DEBUG == 1
    printf("\nAdding preComputeFirstGrid %i\n",preComputeFirstGrid);
    printf("xStart: %.6f\n", allGrids[0].xStart);
    printf("xEnd: %.6f\n", allGrids[0].xEndRequested);
    printf("Length: %i\n", allGrids[0].Length);
    printf("deltax: %.6e\n", allGrids[0].deltax);
    printf("xMax: %.6f\n", allGrids[0].xMax);
    #endif

    fStartInspDerefinement = coarseGrid.xMax;
    df0  = 2*coarseGrid.deltax;
    df0original = 2*df0original;
  }


  // Loop over inspiral derefinement grids
  if (nDerefineInspiralGrids > 0) {
    index = 0;
    nextfSTART = fStartInspDerefinement;
    #if DEBUG == 1
    printf("nDerefineInspiralGrids before loop = %d\n", nDerefineInspiralGrids);
    #endif

    while (index < nDerefineInspiralGrids) {
      if(df0original < evaldMf){
        #if DEBUG == 1
          printf("\nAdjusting freq factors!!\n");
        #endif
        mydf = evaldMf;
        intdfRatio = 1;
      }
      else{
        intdfRatio = (int)floor(df0original/evaldMf);
        mydf = evaldMf*intdfRatio;
      }
      if(index + preComputeFirstGrid == 0){
         fSTART = nextfSTART;
      }
      else{
        fSTART = nextfSTART + mydf;
      }
      fEND = fSTART * FrequencyFactor;

      #if DEBUG == 1
      printf("\n(index, fSTART, fEND) = (%d, %e, %e, %e)\n", index + preComputeFirstGrid, fSTART, fEND, mydf);
      #endif

      coarseGrid = XLALSimIMRPhenomXGridComp(fSTART, fEND, mydf);

      #if DEBUG == 1
      printf("xStart: %.16e\n", coarseGrid.xStart);
      printf("xEnd: %.16e\n", coarseGrid.xEndRequested);
      printf("Length: %i\n", coarseGrid.Length);
      printf("deltax: %.16e\n", coarseGrid.deltax);
      printf("mydf: %.16e\n", mydf);
      printf("xMax: %.16e\n", coarseGrid.xMax);
      printf("intdfRatio = %i\n", intdfRatio);
      #endif

      df0original = 2*df0original;
      //nextmydf   = 2 * mydf;
      nextfSTART = coarseGrid.xMax;

      allGrids[index + preComputeFirstGrid] = coarseGrid;
      allGrids[index+preComputeFirstGrid].intdfRatio = intdfRatio;

      index = index + 1;
    }
    fEndInsp = coarseGrid.xMax;
  }
  else{
    #if DEBUG == 1
    printf("\nSkipping Inspiral Loop %i\n", nDerefineInspiralGrids);
    #endif
    if (preComputeFirstGrid > 0){
      fEndInsp = coarseGrid.xMax;
    }
  }

  #if DEBUG == 1
  printf("\nfStartInspDerefinement after loop = %e\n",  fStartInspDerefinement);
  printf("fEndInsp after loop = %e\n",  fEndInsp);
  printf("nDerefineInspiralGrids = %i\n", nDerefineInspiralGrids);
  #endif

  // Add merger grid
  if (nMergerGrid > 0) {
    df0original = dfmerger; //check if the Delta_f given by the Lorentzian is smaller than at the beginning of the merger bin.
    if(2*coarseGrid.deltax < dfmerger){
      df0original = 2*coarseGrid.deltax;
    }
    if(df0original < evaldMf){
      #if DEBUG == 1
        printf("\nAdjusting freq factors!!\n");
      #endif
      mydf = evaldMf;
      intdfRatio = 1;
    }
    else{
      intdfRatio = (int)floor(df0original/evaldMf);
      mydf = evaldMf*intdfRatio;
    }
    fSTART = fEndInsp + mydf;

    if(fEndInsp == fstartIn){
      fSTART = fEndInsp;
    }

    INT4 mergerIndex = 0;

    if(fSTART > MfLorentzianEnd){
      nMergerGrid = 0;
      #if DEBUG == 1
      printf("\nNOT adding merger grid\n");
      #endif
    }
    else{
      coarseGrid = XLALSimIMRPhenomXGridComp(fSTART, MfLorentzianEnd, mydf);
      mergerIndex = preComputeFirstGrid + nDerefineInspiralGrids;

      df0original = 2*df0original;

      allGrids[mergerIndex] = coarseGrid;
      allGrids[mergerIndex].intdfRatio = intdfRatio;

      #if DEBUG == 1
      printf("\nadding merger grid\n");
      printf("fSTART = %.6f\n", fSTART);
      printf("fEND = %.6f\n", coarseGrid.xMax);
      printf("MfLorentzianEnd = %.6f\n", MfLorentzianEnd);
      printf("mydf = %.16e\n", mydf);
      printf("mergerIndex = %i\n", mergerIndex);
      printf("intdfRatio = %i\n", intdfRatio);
      printf("# fine points float %.16f\n", allGrids[mergerIndex].deltax/evaldMf);
      #endif
    }


  }

  // Add RD grid
  if (nRingdownGrid > 0) {
      df0original = dfringdown;
      if(df0original < evaldMf){
        #if DEBUG == 1
          printf("\nAdjusting freq factors!!\n");
        #endif
        mydf = evaldMf;
        intdfRatio = 1;
      }
      else{
        intdfRatio = (int)floor(df0original/evaldMf);
        mydf = evaldMf*intdfRatio;
      }
      fSTART = coarseGrid.xMax + mydf;
      if(coarseGrid.xMax == fstartIn){
        fSTART = fEndInsp;
      }

      INT4 RDindex = 0;

      if(fSTART > Mfmax){
        nRingdownGrid = 0;
        #if DEBUG == 1
        printf("\nNOT adding RD grid\n");
        #endif
      }
      else{
        coarseGrid = XLALSimIMRPhenomXGridComp(fSTART, Mfmax, mydf);

        RDindex = preComputeFirstGrid + nDerefineInspiralGrids + nMergerGrid;

        df0original = 2*df0original;

        allGrids[RDindex] = coarseGrid;
        allGrids[RDindex].intdfRatio = intdfRatio;

        #if DEBUG == 1
        printf("\nadding RD grid\n");
        printf("Mfmax = %e\n", Mfmax);
        printf("fSTART = %.6f\n", fSTART);
        printf("fEND = %.6f\n", coarseGrid.xMax);
        printf("mydf = %.16e\n", mydf);
        printf("RDIndex = %i\n", RDindex);
        printf("intdfRatio = %i\n", intdfRatio);
        printf("# fine points float %.16f\n", allGrids[RDindex].deltax/evaldMf);
        #endif
      }
  }

  INT4 nGridsUsed = preComputeFirstGrid + nDerefineInspiralGrids + nMergerGrid+nRingdownGrid;
  #if DEBUG == 1
  printf("final grid length = %d\n", nGridsUsed);
  printf("intdfRatio = %d\n", intdfRatio);
  #endif

  return nGridsUsed;
}

/**
 * @addtogroup LALSimIMRPhenomX_c
 * @{
 * @name Routines for IMRPhenomXHM Multibanding
 * @{
 *
 */

/*****************************************************/
/*                                                   */
/*            MULTIBANDING WAVEFORMS                 */
/*                                                   */
/*****************************************************/

/** Return htildelm, the waveform of one mode without mode-mixing.
 e^(I phi) is interpolated with linear order and an iterative procedure.
   Amplitude uses standard 1st (by default) or 3rd order interpolation with gsl. */
int XLALSimIMRPhenomXHMMultiBandOneMode(
  COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
  REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
  UINT4 ell,                           /**< l index of the mode */
  INT4 emmIn,                            /**< m index of the mode */
  REAL8 distance,                      /**< Luminosity distance (m) */
  REAL8 f_min,                         /**< Starting GW frequency (Hz) */
  REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = \ref f_CUT */
  REAL8 deltaF,                        /**< Sampling frequency (Hz) */
  REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
  REAL8 fRef_In,                       /**< Reference frequency (Hz) */
  LALDict *lalParams                   /**< Extra params */
){
  UINT4 emm = abs(emmIn);
  #if DEBUG == 1
  printf("\nMode %i %i \n",ell,emm);
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


  #if DEBUG == 1
  printf("\n**********************************************************************\n");
  printf("\n*                IMRPhenomXHMMultiBandOneMode        %i%i            *\n", ell, emm);
  printf("\n**********************************************************************\n");
  printf("\nm1, m2, chi1, chi2 %.16f %.16f %.16f %.16f\n", m1_SI/LAL_MSUN_SI, m2_SI/LAL_MSUN_SI, chi1L, chi2L);
  #endif


  /* When mode does not vanishes */
  int debug = DEBUG;

  // Define two powers of pi to avoid clashes between PhenomX and PhenomXHM files.
  int status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PIHM.");

  /* Initialize IMRPhenomX Waveform struct and check that it initialized correctly */
  IMRPhenomXWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, deltaF, fRef_In, phiRef, f_min, f_max, distance, 0.0, lalParams, debug);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");


  status = IMRPhenomXHMMultiBandOneMode(htildelm, pWF, ell, emm, lalParams);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHMMultiBandOneMode failed to generate IMRPhenomXHM waveform.");

  INT4 offset = (size_t) (pWF->fMin / deltaF);

  if(emmIn>0){
    /* (-1)^l */
   INT4 minus1l = 1;
   if (ell % 2 !=0){
     minus1l = -1;
   }
   #if DEBUG == 1
    printf("\nTransforming to positive m by doing (-1)^l*Conjugate, frequencies must be negatives.\n");
    #endif
    for(UINT4 idx=offset; idx<(*htildelm)->data->length; idx++){
      (*htildelm)->data->data[idx] = minus1l*conj((*htildelm)->data->data[idx]);
    }
  }

  REAL8 lastfreq;
  /* Resize htildelm if needed */
  if (pWF->f_max_prime < pWF->fMax)
  { /* The user has requested a higher f_max than Mf = fCut.
    Resize the frequency series to fill with zeros beyond the cutoff frequency. */
    lastfreq = pWF->fMax;
  }
  else{  // We have to look for a power of 2 anyway. Without MBAND this step is already satisfied
    lastfreq = pWF->f_max_prime;
  }
  // We want to have the length be a power of 2 + 1
  size_t n_full = NextPow2(lastfreq / deltaF) + 1;
  size_t n = (*htildelm)->data->length;

  /* Resize the COMPLEX16 frequency series */
  *htildelm = XLALResizeCOMPLEX16FrequencySeries(*htildelm, 0, n_full);
  XLAL_CHECK (*htildelm, XLAL_ENOMEM, "Failed to resize waveform COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );
  XLALUnitMultiply(&((*htildelm)->sampleUnits), &((*htildelm)->sampleUnits), &lalSecondUnit);
  
  LALFree(pWF);

  return offset;

}
/** @}
 @} **/


int IMRPhenomXHMMultiBandOneMode(
    COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform **/
    IMRPhenomXWaveformStruct *pWF,       /**< Waveform structure 22 mode **/
    UINT4 ell,                           /**< First index (l,m) mode **/
    UINT4 emm,                           /**< Second incex (l,m) mode **/
    LALDict *lalParams                   /**< LAL dictionary **/
)
{
  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  REAL8 deltaF = pWF->deltaF;
  pWF->deltaF = 0; // Needed for SetupWFArraysReal function, if not introduces an extra offset in amplitude and phase.

  /* Check that the frequency array will be consistent: fmin < fmax_prime */
  /* Return the closest power of 2 */
  size_t npts = NextPow2(pWF->f_max_prime / deltaF) + 1;
  /* Frequencies will be set using only the lower and upper bounds that we passed */
  size_t iStart = (size_t) (pWF->fMin / deltaF);
  size_t iStop  = (size_t) (pWF->f_max_prime / deltaF) + 1;
  XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
  "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);

  #if DEBUG == 1
  printf("\n***********************\n");
  printf("pWF->fMin, deltaF, iStart = %.16e %.16e %zu", pWF->fMin, deltaF, iStart);
  printf("\n***********************\n");
  #endif


  size_t offset = iStart;


  /* If it is odd mode with equal black holes then return array of zeros 0 */
  if(pWF->q == 1 && pWF->chi1L == pWF->chi2L && emm%2!=0){  // Mode zero
    *htildelm = XLALCreateCOMPLEX16FrequencySeries("htildelm: FD waveform", &(ligotimegps_zero), 0.0, deltaF, &lalStrainUnit, iStop);
    for(unsigned int idx = 0; idx < iStop; idx++){
      ((*htildelm)->data->data)[idx] = 0.;
    }

    return XLAL_SUCCESS;
  }

  /** Mode non-zero **/


  /* Final grid spacing, adimensional (NR) units */
  REAL8 evaldMf = XLALSimIMRPhenomXUtilsHztoMf(deltaF, pWF->Mtot);

  /* Threshold for the Multibanding. It is a measure of how restrictive it is. Smaller value implies stronger restriction, more points where evaluate the model. */
  REAL8 resTest  = XLALSimInspiralWaveformParamsLookupPhenomXHMThresholdMband(lalParams);
  UINT4 ampinterpolorder = XLALSimInspiralWaveformParamsLookupPhenomXHMAmpInterpolMB(lalParams);
  #if DEBUG == 1
  printf("\n***** MBAND = %i, resTest = %.16f, ampIntorder = %i\n", MBAND, resTest, ampinterpolorder);
  #endif
  /* Variable for the Multibanding criteria */
  REAL8 dfpower = 11./6.;
  REAL8 dfcoefficient = 8. * sqrt(3./5.) * LAL_PI * powers_of_lalpi.m_one_sixth * sqrt(2.)*cbrt(2) /(cbrt(emm)*emm) * sqrt(resTest * pWF->eta);
  /* Variables for the coarse frequency grid */
  REAL8 Mfmin = XLALSimIMRPhenomXUtilsHztoMf(iStart*deltaF, pWF->Mtot);
  REAL8 Mfmax = XLALSimIMRPhenomXUtilsHztoMf(pWF->f_max_prime, pWF->Mtot);
  REAL8 MfMECO, MfLorentzianEnd;
  REAL8 dfmerger = 0., dfringdown = 0.;//, relerror = 0.001;


  IMRPhenomXHMWaveformStruct *pWFHM = (IMRPhenomXHMWaveformStruct *) XLALMalloc(sizeof(IMRPhenomXHMWaveformStruct));
  //Initialize pWFHM, pAmp(22), pPhase(22) and pass to the functions. No l, m needed
  //populate coefficients of 22 mode, to rotate to spherical
  IMRPhenomXAmpCoefficients   *pAmp22   = (IMRPhenomXAmpCoefficients *) XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
  IMRPhenomXPhaseCoefficients *pPhase22 = (IMRPhenomXPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
  IMRPhenomXGetPhaseCoefficients(pWF, pPhase22);

  /* Allocate and initialize the PhenomXHM lm amplitude coefficients struct */
  IMRPhenomXHMAmpCoefficients *pAmp = (IMRPhenomXHMAmpCoefficients*)XLALMalloc(sizeof(IMRPhenomXHMAmpCoefficients));
  IMRPhenomXHMPhaseCoefficients *pPhase = (IMRPhenomXHMPhaseCoefficients*)XLALMalloc(sizeof(IMRPhenomXHMPhaseCoefficients));

  if(ell == 2 && emm ==2){
    MfMECO = pWF->fMECO;
    #if DEBUG == 1
    printf("\nfRING = %e\n",pWF->fRING);
    printf("fDAMP = %e\n",pWF->fDAMP);
    printf("alphaL22 = %.16e", pPhase22->cLovfda/pWF->eta);
    #endif
    MfLorentzianEnd = pWF->fRING + 2*pWF->fDAMP;
    IMRPhenomXGetAmplitudeCoefficients(pWF, pAmp22);
    dfmerger = deltaF_mergerBin(pWF->fDAMP, pPhase22->cLovfda/pWF->eta, resTest);
    dfringdown = deltaF_ringdownBin(pWF->fDAMP, pPhase22->cLovfda/pWF->eta, pAmp22->gamma2/(pAmp22->gamma3*pWF->fDAMP),resTest);
  }
  else{
    // allocate qnm struct
    QNMFits *qnms = (QNMFits *) XLALMalloc(sizeof(QNMFits));
    IMRPhenomXHM_Initialize_QNMs(qnms);
    // Populate pWFHM
    IMRPhenomXHM_SetHMWaveformVariables(ell, emm, pWFHM, pWF,qnms, lalParams);
    LALFree(qnms);

    /* Allocate and initialize the PhenomXHM lm phase and amp coefficients struct */
    IMRPhenomXHM_FillAmpFitsArray(pAmp);
    IMRPhenomXHM_FillPhaseFitsArray(pPhase);

    /* Get coefficients for Amplitude and phase */
    IMRPhenomXHM_GetAmplitudeCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF);
    IMRPhenomXHM_GetPhaseCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF,lalParams);

    MfMECO = pWFHM->fMECOlm;
    MfLorentzianEnd = pWFHM->fRING + 2*pWFHM->fDAMP;
    #if DEBUG == 1
    printf("\nfRING = %e\n",pWFHM->fRING);
    printf("fDAMP = %e\n",pWFHM->fDAMP);
    printf("alphaL = %.16e", pPhase->alphaL);
    #endif
    dfmerger = deltaF_mergerBin(pWFHM->fDAMP, pPhase->alphaL, resTest);
    dfringdown = deltaF_ringdownBin(pWFHM->fDAMP, pPhase->alphaL, pAmp->lambda/(pAmp->sigma*pWFHM->fDAMP), resTest);
  }

  /* Allocate memory for the list of grids. The number of grids must be less than lengthallGrids. */
  UINT4 lengthallGrids = 20;
  IMRPhenomXMultiBandingGridStruct *allGrids = (IMRPhenomXMultiBandingGridStruct*)XLALMalloc(lengthallGrids * sizeof(IMRPhenomXMultiBandingGridStruct));

  if (allGrids == NULL)
  {
    #if DEBUG == 1
    printf("Malloc of allGrids failed!\n");
    #endif
    return -1;
  }

  #if DEBUG == 1
  printf("\nMfmin = %.6f\n", Mfmin);
  printf("MfMECO = %.6f\n", MfMECO);
  printf("MfLorentzianEnd = %.6f\n", MfLorentzianEnd);
  printf("Mfmax = %.6f\n", Mfmax);
  printf("evaldMf = %.6e\n", evaldMf);
  printf("dfpower = %.6e\n", dfpower);
  printf("dfcoefficient = %.6e\n", dfcoefficient);
  printf("dfmerger = %.6e\n", dfmerger);
  printf("dfringdown = %.6e\n", dfringdown);
  #endif

  /* Compute the coarse frequency array. It is stored in a list of grids. */
  UINT4 nGridsUsed = XLALSimIMRPhenomXMultibandingGrid(Mfmin, MfMECO, MfLorentzianEnd, Mfmax, evaldMf, dfpower, dfcoefficient, allGrids, dfmerger, dfringdown);

  #if DEBUG == 1
  printf("allGrids[1].Length = %i\n", allGrids[0].Length);
  #endif

  /* Number of fine frequencies per coarse interval in every coarse grid */
  INT4 mydfRatio[lengthallGrids];
  /* Actual number of subgrids to be used. We allocated more than needed. */
  UINT4 actualnumberofGrids = 0;
  /* Length of coarse frequency array */
  UINT4 lenCoarseArray = 0;

  /* Transform the coarse frequency array to 1D array. */
  // Take only the subgrids needed
  for(UINT4 kk = 0; kk < nGridsUsed; kk++){
    lenCoarseArray = lenCoarseArray + allGrids[kk].Length;
    actualnumberofGrids++;

    mydfRatio[kk] = allGrids[kk].intdfRatio;

    #if DEBUG == 1
    printf("\nkk = %i\n",kk);
    printf("xStart: %.6e\n", allGrids[kk].xStart);
    printf("xEnd: %.6e\n", allGrids[kk].xEndRequested);
    printf("Length: %i\n", allGrids[kk].Length);
    printf("deltax, Hz: %.6e %.6e\n", allGrids[kk].deltax, XLALSimIMRPhenomXUtilsMftoHz(allGrids[kk].deltax, pWF->Mtot));
    printf("evaldMf, Hz: %.6e %.6e\n", evaldMf, XLALSimIMRPhenomXUtilsMftoHz(evaldMf, pWF->Mtot));
    printf("xMax: %.16e\n", allGrids[kk].xMax);
    printf("Mfmax: %.16e\n", Mfmax);
    printf("# fine points %i\n", mydfRatio[kk]);
    printf("# fine points float %.16f\n", allGrids[kk].deltax/evaldMf);
    #endif

    if(allGrids[kk].xMax + evaldMf >= Mfmax){
      break;
    }
  }

  // Add extra points to the coarse grid if the last freq is lower than Mfmax
  while(allGrids[actualnumberofGrids-1].xMax < Mfmax){
    allGrids[actualnumberofGrids-1].xMax =   allGrids[actualnumberofGrids-1].xMax +  allGrids[actualnumberofGrids-1].deltax;
    allGrids[actualnumberofGrids-1].Length =  allGrids[actualnumberofGrids-1].Length + 1;
    lenCoarseArray++;
  }

  #if DEBUG == 1
  if(ell==2 && emm==2){
    printf("\nfDAMP = %.16e\n", XLALSimIMRPhenomXUtilsMftoHz(pWF->fDAMP, pWF->Mtot));
  }else{
    printf("\nfDAMP = %.16e\n", XLALSimIMRPhenomXUtilsMftoHz(pWFHM->fDAMP, pWF->Mtot));
  }
  printf("actualnumberofGrids = %i\n", actualnumberofGrids);
  printf("lenCoarseArray = %i\n", lenCoarseArray);
  printf("Last grid.xMax = %.16f", allGrids[actualnumberofGrids-1].xMax);
  #endif

  // Transform coarse frequency array to 1D vector
  REAL8 *IntLawpoints = (REAL8*)XLALMalloc(lenCoarseArray * sizeof(REAL8));
  UINT4 lenIntLawpoints = 0;

  for(UINT4 kk = 0; kk < actualnumberofGrids; kk++){
    for(INT4 ll = 0; ll < allGrids[kk].Length; ll++){
      IntLawpoints[lenIntLawpoints] = (allGrids[kk].xStart + allGrids[kk].deltax*ll);
      lenIntLawpoints++;
    }
  }
  /* End of coarse frequency array. */

  #if DEBUG == 1
  printf("\n******** Coarse frequencies array done ********* \n");
  printf("\nlenIntLawpoints, coarse[0], coarse[-1], Mfmax, M_sec = %i %.16e %.16e %.16e %.16e\n",lenIntLawpoints, IntLawpoints[0], IntLawpoints[lenIntLawpoints-1], Mfmax, pWF->M_sec);
  #endif

  /* IntLawpoints stores the adimensional frequencies (Mf) */
  /* coarseFreqs will store the same frequency array but in Hz */
  /* Allocate memory for frequency array and terminate if this fails */
  REAL8Sequence *coarseFreqs;
  coarseFreqs = XLALCreateREAL8Sequence(lenCoarseArray);
  if (!coarseFreqs) {XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed."); }

  /* Populate frequency array */
  #if DEBUG == 1
  printf("\n***** Coarse freqs *****\n");
  #endif
  REAL8 divis = 1./pWF->M_sec;
  for (UINT4 ii = 0; ii < lenCoarseArray; ii++)
  {
    coarseFreqs->data[ii] = IntLawpoints[ii]*divis;
  }
  #if DEBUG == 1
  printf("\nFirst/Last Coarse freqs *****%.16e %.16e\n", coarseFreqs->data[0], coarseFreqs->data[coarseFreqs->length-1]);
  #endif

  /* Next we will compute amplitude and phase of one mode in the coarse frequency array. */

  REAL8FrequencySeries *amplitude, *phase;

  /** Compute 22 using PhenomX functions **/
  // This is copied from IMRPhenomXASGenerateFD.
  if(ell == 2 && emm ==2){
    #if DEBUG == 1
    printf("\n** Computing Amplitude and Phase of PhenomX %i **********\n", lenCoarseArray);
    #endif

    amplitude = XLALCreateREAL8FrequencySeries("amplitude22: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,coarseFreqs->length);
    phase     = XLALCreateREAL8FrequencySeries("phase22: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,coarseFreqs->length);

    IMRPhenomXGetAmplitudeCoefficients(pWF, pAmp22);

    /* Initialize a struct containing useful powers of Mf at fRef */
    IMRPhenomX_UsefulPowers powers_of_MfRef;
    int status = IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF->MfRef);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "IMRPhenomX_Initialize_Powers failed for MfRef.\n");


    /* Linear time and phase shifts so that model peaks near t ~ 0 */
    REAL8 lina = 0;

    IMRPhenomX_Phase_22_ConnectionCoefficients(pWF,pPhase22);
    double linb=IMRPhenomX_TimeShift_22(pPhase22, pWF);

    // Calculate IMRPhenomX phase at reference frequency
    REAL8 phiref22 = -1./pWF->eta*IMRPhenomX_Phase_22(pWF->MfRef, &powers_of_MfRef, pPhase22, pWF) - linb*pWF->MfRef - lina + 2.0*pWF->phi0 + LAL_PI_4;


    for(UINT4 kk = 0; kk < (coarseFreqs)->length; kk++){
      REAL8 Mff = coarseFreqs->data[kk]*pWF->M_sec;
      IMRPhenomX_UsefulPowers powers_of_f;
      IMRPhenomX_Initialize_Powers(&powers_of_f,Mff);
      amplitude->data->data[kk] = IMRPhenomX_Amplitude_22(Mff, &powers_of_f, pAmp22, pWF) * pWF->amp0;
      phase->data->data[kk] = 1./pWF->eta*IMRPhenomX_Phase_22(Mff, &powers_of_f, pPhase22, pWF) + linb*Mff + lina + phiref22;
    }
  }
  /** Higher modes **/
  else{

    #if DEBUG == 1
    printf("\n******* Computing Coarse Amplitude And Phase ****************\n");
    #endif
    /* Compute coarse amplitude and phase */
    IMRPhenomXHM_Amplitude(&amplitude, coarseFreqs, pWF, pAmp22, pPhase22, pWFHM, pAmp, pPhase);
    IMRPhenomXHM_Phase(&phase, coarseFreqs, pWF, pAmp22, pPhase22, pWFHM, pAmp, pPhase);
  }

  #if DEBUG == 1
  printf("\n******* Computed Coarse Amp and Phase **************** %i\n", coarseFreqs->length);
  #endif

  /* Transform the REAL8FrequencySeries to vector since gsl uses vectors for the interpolation.
     gsl is only used for the amplitude but to keep the code more symmetric we transform also the phase. */
  REAL8 *ILamplm = (REAL8*)XLALMalloc(lenCoarseArray * sizeof(REAL8));
  REAL8 *ILphaselm = (REAL8*)XLALMalloc(lenCoarseArray * sizeof(REAL8));

  for (UINT4 ii = 0; ii < lenCoarseArray; ii++)
  {
    ILamplm[ii] = ((amplitude)->data->data)[ii];
    ILphaselm[ii] = ((phase)->data->data)[ii];
  }

  #if DEBUG == 1
  //Save coarse amplitude and phase in file
  FILE *file0;
  char fileSpec0[40];
  sprintf(fileSpec0, "coarseamplitude%i%i.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec0);
  file0 = fopen(fileSpec0,"w");
  fprintf(file0,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->q, pWF->chi1L, pWF->chi2L, ell, emm);
  FILE *file3;
  char fileSpec3[40];
  sprintf(fileSpec3, "coarsephase%i%i.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec3);
  file3 = fopen(fileSpec3,"w");
  fprintf(file3,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->q, pWF->chi1L, pWF->chi2L, ell, emm);
  for(UINT4 idx = 0; idx < lenCoarseArray; idx++)
  {
    fprintf(file0, "%.16f  %.16e \n",  coarseFreqs->data[idx], ILamplm[idx]);
    fprintf(file3, "%.16f  %.16e \n",  coarseFreqs->data[idx], ILphaselm[idx]);
  }
  fclose(file0);
  fclose(file3);
  #endif

  /* Free allocated memory */
  XLALDestroyREAL8FrequencySeries(amplitude);
  XLALDestroyREAL8FrequencySeries(phase);
  XLALDestroyREAL8Sequence(coarseFreqs);

  /***** Linear Interpolation of e^(I*phi) with the iterative procedure *****/

  /* Estimation of the length of the fineGrid */
  INT4 lenWF = 0;

  for(UINT4 kk = 0; kk < actualnumberofGrids; kk++){
    lenWF = lenWF + (allGrids[kk].Length -1) * mydfRatio[kk] + 2*mydfRatio[kk];
    #if DEBUG == 1
    printf("\nmydfRatio[%i] = %i %i %i", kk, mydfRatio[kk],allGrids[kk].Length, (allGrids[kk].Length -1) * mydfRatio[kk] + 2*mydfRatio[kk]);
    printf("\nlenWF = %i\n\n", lenWF);
    #endif
  }

  // Variable to control the number of coarse points we have used
  int pointsPrecessedSoFar = 0;

  /* Allocate memory for the complex exponential and the fine equally-spaced frequency array */
  COMPLEX16 *expphi = (COMPLEX16*)XLALMalloc(lenWF*sizeof(COMPLEX16));
  if (expphi == NULL){
    return -1;
  }

  REAL8 *finefreqs = (REAL8*)XLALMalloc(lenWF*sizeof(REAL8));
  if (finefreqs == NULL){
    return -1;
  }

  UINT4 count = 0; // Variable to track the point being filled in the fine frequency grid

  /* Loop over allgrids */
  bool stop = false;
  for (UINT4 i = 0; i<actualnumberofGrids && !stop; i++){
    #if DEBUG == 1
    printf("\ni = %i\n", i);
    #endif

    UINT4 lcoarseGrid = allGrids[i].Length;
    if(lcoarseGrid == 0){
      break;
    }

    /* Linear Interpolation and iteration, here I get fineGridResult */
    if(i==actualnumberofGrids-1){
      lcoarseGrid--;
    }

    REAL8 Omega, phi0, Mfhere = 0, Mfnext = 0;
    COMPLEX16 h0, Q;

    /* Loop over the coarse points of a subgrid */
    for(UINT4 j = 0; j < lcoarseGrid && !stop ; j++){
      Mfhere = IntLawpoints[pointsPrecessedSoFar + j];
      Mfnext = IntLawpoints[pointsPrecessedSoFar + j + 1];

      INT4 ratio;
      // If we are in the last coarse point of a sublist we use the number of fine points of the next subgrid.
      if(j==lcoarseGrid-1 && i < actualnumberofGrids-1){
        ratio = mydfRatio[i+1];
      }
      else{
        ratio = mydfRatio[i];
      }
      if(Mfnext + evaldMf >= Mfmax){
        double dratio = (Mfmax - Mfhere)/evaldMf + 1;
        ratio = (int) dratio;
        int roundratio = round((Mfmax - Mfhere)/evaldMf) + 1;
        if(fabs(dratio-roundratio) < 0.0001) ratio = roundratio; //To get the correct rounded integer
        /* Break the loop if you overpass the maximun frequency */
        stop = true;
        #if DEBUG == 1
        printf("\nMfmax, Mfhere, evaldMf, ratio, lastfreqHz= %.16f %.16f %.16f %i %.16e\n", Mfmax, Mfhere, evaldMf, ratio, XLALSimIMRPhenomXUtilsMftoHz(Mfhere+(ratio-1)*evaldMf,pWF->Mtot) );
        #endif
      }

      /********************************************/
      /*     Inner Loop: linear interpolation     */
      /********************************************/
      UINT4 jjdx = j + pointsPrecessedSoFar;
      if(jjdx < lenCoarseArray){
        Omega = (ILphaselm[jjdx+ 1] - ILphaselm[jjdx])/(IntLawpoints[jjdx + 1] - IntLawpoints[jjdx]);
      }
      else{
        Omega = (ILphaselm[jjdx] - ILphaselm[jjdx -1])/(IntLawpoints[jjdx] - IntLawpoints[jjdx -1]);
      }
      phi0 = ILphaselm[jjdx];

      h0   = cexp(I*phi0);
      Q    = cexp(I*evaldMf*Omega);

      finefreqs[count] = Mfhere;
      expphi[count]    = h0;
      count++;

      /* This loop carry out the eq. 2.32 in arXiv:2001.10897 */
      for(int kk = 1; kk < ratio; kk++){       // Compute finefreqs and fine expphi
        finefreqs[count] = Mfhere + evaldMf*kk;
        expphi[count]    = Q*expphi[count-1];
        count++;
      }

    }
    pointsPrecessedSoFar = pointsPrecessedSoFar + lcoarseGrid;
  }// End loop over ILgrids. count should be aprox = to lenWF

  #if DEBUG == 1
  printf("\ncount = %i\n", count);
  #endif

  #if DEBUG == 1
  printf("\n******* Interpolate Amplitude **********\n");
  #endif

  /*** Interpolation and evaluation in fineFreqs of the amplitude ***/
  REAL8 *fineAmp = (REAL8*)XLALMalloc(count * sizeof(REAL8));
  interpolateAmplitude(fineAmp, IntLawpoints, ILamplm, finefreqs, lenCoarseArray, count, ampinterpolorder);

  /**** Build the waveform ****/
  // Due to round of erros, the last freq may be greater Mfmax. The difference should be less than the frequency step.
  // Remove that extra frequency point
  while(finefreqs[count-1] > Mfmax){
    count--;
  }
  #if DEBUG == 1
  printf("\n******* Building Waveform 2**********\n");
  printf("\nfinefreqs[0]Hz = %.16f\n", finefreqs[0]/pWF->M_sec);
  printf("\nfinefreqs[%i]Hz = %.16f\n", count-1, finefreqs[count-1]/pWF->M_sec);
  printf("Mfmax in Hz = %.16f\n", Mfmax/pWF->M_sec);
  printf("count, offset = %i %zu\n", count-1, offset);
  printf("Mfmaxtheory = %.16f\n", (count-1+offset)*deltaF);
  #endif

  /* Intialize FrequencySeries for htildelm */
  while(count+offset > iStop){ count--; }
  size_t n = iStop;
  XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / deltaF), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.", deltaF);
  *htildelm = XLALCreateCOMPLEX16FrequencySeries("htildelm: FD waveform", &(ligotimegps_zero), 0.0, deltaF, &lalStrainUnit, n);

  for(int idx = 0; idx < (int) offset; idx++){
    ((*htildelm)->data->data)[idx] = 0.;
  }

   /* (-1)^l */
  INT4 minus1l;
  if (ell % 2 !=0){
    minus1l = -1;
  }
  else{
    minus1l = +1;
  }

  for(UINT4 idx = 0; idx < count; idx++){
    /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
    ((*htildelm)->data->data)[idx + offset] = minus1l * fineAmp[idx] * expphi[idx];
  }

  /* Sometimes rounding error can make that count+offset < iStop, meaning that we have computed less points than those we are gonna output,
    here we make sure that all the elements of htildelm are initialized and we put the extra point(s) to 0. */
  for(UINT4 idx = count-1; idx < iStop-offset; idx++){
    /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
    ((*htildelm)->data->data)[idx + offset] = 0.;
  }

  /* Free allocated memory */
  LALFree(pAmp);
  LALFree(pAmp22);
  LALFree(pPhase);
  LALFree(pPhase22);


  #if DEBUG == 1
  //Save hlm in file
  FILE *file;
  char fileSpec[40];
  sprintf(fileSpec, "simulation%i%i_Multiband.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec);
  file = fopen(fileSpec,"w");
  fprintf(file,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->q, pWF->chi1L, pWF->chi2L, ell, emm);
  fprintf(file,"# Frequency (Hz)    Real    Imaginary\n");
  COMPLEX16 data;
  for(UINT4 idx = 0; idx < count; idx++)
  {
    data = expphi[idx];
    fprintf(file, "%.16f  %.16e %.16e\n",  (idx+offset)*deltaF, creal(data), cimag(data));
  }
  fclose(file);
  #endif

  #if DEBUG == 1
  //Save fine amplitude and freq array in file
  FILE *file2;
  char fileSpec2[40];
  sprintf(fileSpec2, "amplitude%i%i_fine.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec2);
  file2 = fopen(fileSpec2,"w");
  fprintf(file2,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->q, pWF->chi1L, pWF->chi2L, ell, emm);
  printf("\ncount, count + offset, len htildelm = %i %i %i\n", count, (int)(count + offset), (*htildelm)->data->length);
  for(UINT4 idx = 0; idx < count; idx++)
  {
    fprintf(file2, "%.16f  %.16f %.16e\n",  finefreqs[idx], (idx+offset)*evaldMf, fineAmp[idx]);
  }
  fclose(file2);
  #endif

  /* Free allocated memory */
  LALFree(expphi);
  LALFree(finefreqs);
  LALFree(allGrids);
  LALFree(fineAmp);
  LALFree(pWFHM);
  LALFree(ILamplm);
  LALFree(ILphaselm);
  LALFree(IntLawpoints);

  return XLAL_SUCCESS;
}


/** @addtogroup LALSimIMRPhenomX_c
* @{
* @name Routines for IMRPhenomXHM Multibanding
* @{
* @author Cecilio García Quirós, Sascha Husa
*
* @brief C code for applying Multibanding to IMRPhenomXHM_Multimode
*
* This is a technique to make the evaluation of waveform models faster by evaluating the model
* in a coarser non-uniform grid and interpolate this to the final fine uniform grid.
* We apply this technique to the fourier domain model IMRPhenomXHM as describe in this paper: https://arxiv.org/abs/2001.10897
*
* Multibanding flags:
*   ThresholdMband: Determines the strength of the Multibanding algorithm.
*   The lower this value is, the slower is the evaluation but more accurate is the final waveform compared to the one without multibanding.
*         - 0.001: DEFAULT value
*         - 0: switch off the multibanding
*
*   AmpInterpol: Determines the gsl interpolation order for the amplitude.
*         - 1: linear interpolation (DEFAULT)
*         - 3: cubic interpolation
*
*/

/** Returns htildelm the waveform of one mode that present mode-mixing.
The multibanding is applied to the spherical part (inspiral and intermediate) and to the spheroidal part (ringdown).
Both are interpolated and evaluated in the fine frequency grid and the ringdown part is rotated back to spherical */
int XLALSimIMRPhenomXHMMultiBandOneModeMixing(
  COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
  COMPLEX16FrequencySeries *htilde22,  /**< Precomputed FD waveform of dominant mode */
  REAL8 m1_SI,                         /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                         /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                         /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                         /**< Dimensionless aligned spin of companion 2 */
  UINT4 ell,                           /**< l index of the mode */
  INT4 emmIn,                          /**< m index of the mode */
  REAL8 distance,                      /**< Luminosity distance (m) */
  REAL8 f_min,                         /**< Starting GW frequency (Hz) */
  REAL8 f_max,                         /**< End frequency; 0 defaults to Mf = \ref f_CUT */
  REAL8 deltaF,                        /**< Sampling frequency (Hz) */
  REAL8 phiRef,                        /**< Orbital phase at fRef (rad) */
  REAL8 fRef_In,                       /**< Reference frequency (Hz) */
  LALDict *lalParams                   /**< Extra params */
){
  UINT4 emm = abs(emmIn);
  #if DEBUG == 1
  printf("\nMode %i %i \n",ell,emm);
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

  #if DEBUG == 1
  printf("\n**********************************************************************\n");
  printf("\n*                IMRPhenomXHMMultiBandOneModeMixing  %i%i            *\n", ell, emm);
  printf("\n**********************************************************************\n");
  printf("\nm1, m2, chi1, chi2, f_min, f_max, deltaF  %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", m1_SI/LAL_MSUN_SI, m2_SI/LAL_MSUN_SI, chi1L, chi2L, f_min, f_max, deltaF);
  if(htilde22 == NULL){
    printf("*** 22 mode not computed before ***\n\n");
  }
  #endif


  int debug = DEBUG;

  // Define two powers of pi to avoid clashes between PhenomX and PhenomXHM files.
  int status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.");
  status = IMRPhenomX_Initialize_Powers(&powers_of_lalpiHM, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PIHM.");

  /* Initialize IMRPhenomX Waveform struct and check that it initialized correctly */
  IMRPhenomXWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, deltaF, fRef_In, phiRef, f_min, f_max, distance, 0.0, lalParams, debug);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");


  //int offset = IMRPhenomXHMMultiBandOneModeMixing(htildelm, htilde22, pWF, ell, emm, lalParams);
  status = IMRPhenomXHMMultiBandOneModeMixing(htildelm, htilde22, pWF, ell, emm, lalParams);
  XLAL_CHECK(status == XLAL_SUCCESS, XLAL_EFUNC, "IMRPhenomXHMMultiBandOneModeMixing failed to generate IMRPhenomXHM waveform.");

  INT4 offset = (size_t) (pWF->fMin / deltaF);

  if(emmIn>0){
    /* (-1)^l */
    INT4 minus1l = 1;
    if (ell % 2 != 0){
        minus1l = -1;
    }
    #if DEBUG == 1
    printf("\nTransforming to positive m by doing (-1)^l*Conjugate, frequencies must be negatives.\n");
    #endif
    for(UINT4 idx=offset; idx<(*htildelm)->data->length; idx++){
      (*htildelm)->data->data[idx] = minus1l*conj((*htildelm)->data->data[idx]);
    }
  }

  REAL8 lastfreq;
  /* Resize htildelm if needed */
  if (pWF->f_max_prime < pWF->fMax)
  { /* The user has requested a higher f_max than Mf = fCut.
    Resize the frequency series to fill with zeros beyond the cutoff frequency. */
    lastfreq = pWF->fMax;
  }
  else{  // We have to look for a power of 2 anyway. Without MBAND this step is already satisfied
    lastfreq = pWF->f_max_prime;
  }
  // We want to have the length be a power of 2 + 1
  size_t n_full = NextPow2(lastfreq / deltaF) + 1;
  size_t n = (*htildelm)->data->length;

  /* Resize the COMPLEX16 frequency series */
  *htildelm = XLALResizeCOMPLEX16FrequencySeries(*htildelm, 0, n_full);
  XLAL_CHECK (*htildelm, XLAL_ENOMEM, "Failed to resize waveform COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, pWF->fCut, n_full, pWF->fMax );
  XLALUnitMultiply(&((*htildelm)->sampleUnits), &((*htildelm)->sampleUnits), &lalSecondUnit);
  
  LALFree(pWF);

  return offset;

}

/** @}
* @}
*/


int IMRPhenomXHMMultiBandOneModeMixing(
    COMPLEX16FrequencySeries **htildelm, /**< [out] FD waveform */
    COMPLEX16FrequencySeries *htilde22,  /**< Recycle the 22 mode if previously computed **/
    IMRPhenomXWaveformStruct *pWF,       /**< Structure of 22 mode **/
    UINT4 ell,                           /**< First index (l,m) mode **/
    UINT4 emm,                           /**< Second incex (l,m) mode **/
    LALDict *lalParams                   /**< LAL dictionary **/
)
{
  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  REAL8 deltaF = pWF->deltaF;
  pWF->deltaF = 0; // Needed for SetupWFArraysReal function, if not introduces an extra offset in amplitude and phase.

  /* Check that the frequency array will be consistent: fmin < fmax_prime */
  /* Return the closest power of 2 */
  size_t npts = NextPow2(pWF->f_max_prime / deltaF) + 1;
  /* Frequencies will be set using only the lower and upper bounds that we passed */
  size_t iStart = (size_t) (pWF->fMin / deltaF);
  size_t iStop  = (size_t) (pWF->f_max_prime / deltaF) + 1;
  XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
  "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);

  // Allocate qnm struct
  QNMFits *qnms = (QNMFits *) XLALMalloc(sizeof(QNMFits));
  IMRPhenomXHM_Initialize_QNMs(qnms);

  // Populate pWFHM with useful parameters of each mode
  IMRPhenomXHMWaveformStruct *pWFHM = (IMRPhenomXHMWaveformStruct *) XLALMalloc(sizeof(IMRPhenomXHMWaveformStruct));
  IMRPhenomXHM_SetHMWaveformVariables(ell, emm, pWFHM, pWF,qnms, lalParams);
  LALFree(qnms);

  /* Final grid spacing, adimensional (NR) units */
  REAL8 evaldMf = XLALSimIMRPhenomXUtilsHztoMf(deltaF, pWF->Mtot);

  /* Threshold for the Multibanding. It is a measure of how restrictive it is. Smaller value implies stronger restriction, more points where evaluate the model. */
  REAL8 resTest  = XLALSimInspiralWaveformParamsLookupPhenomXHMThresholdMband(lalParams);
  UINT4 ampinterpolorder = XLALSimInspiralWaveformParamsLookupPhenomXHMAmpInterpolMB(lalParams);
  #if DEBUG == 1
  printf("\n***** MBAND = %i, resTest = %.16f, ampIntorder = %i\n", MBAND, resTest, ampinterpolorder);
  #endif
  /* Variable for the Multibanding criteria */
  REAL8 dfpower = 11./6.;
  REAL8 dfcoefficient = 8. * sqrt(3./5.) * LAL_PI * powers_of_lalpi.m_one_sixth * sqrt(2.)*cbrt(2) /(cbrt(emm)*emm) * sqrt(resTest * pWF->eta);
  /* Variables for the coarse frequency grid */
  REAL8 Mfmin = XLALSimIMRPhenomXUtilsHztoMf(iStart*deltaF, pWF->Mtot);
  REAL8 Mfmax = XLALSimIMRPhenomXUtilsHztoMf(pWF->f_max_prime, pWF->Mtot);
  REAL8 MfMECO = pWFHM->fMECOlm;                          // Separate the inspiral grid from the merger grid
  REAL8 MfLorentzianEnd = pWFHM->fRING + 2*pWFHM->fDAMP;  // Separate the merger grid from the ringdown grid
  REAL8 dfmerger = 0., dfringdown = 0.;

  //Initialize pWFHM, pAmp(22), pPhase(22) and pass to the functions. No l, m needed
  //populate coefficients of 22 mode, to rotate to spherical
  IMRPhenomXAmpCoefficients   *pAmp22   = (IMRPhenomXAmpCoefficients *) XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
  IMRPhenomXPhaseCoefficients *pPhase22 = (IMRPhenomXPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
  IMRPhenomXGetPhaseCoefficients(pWF, pPhase22);
  //IMRPhenomX_Phase_22_ConnectionCoefficients(pWF,pPhase22);//ceci where should this go? discontinuity

  /* Allocate and initialize the PhenomXHM lm amplitude coefficients struct */
  IMRPhenomXHMAmpCoefficients *pAmp = (IMRPhenomXHMAmpCoefficients*)XLALMalloc(sizeof(IMRPhenomXHMAmpCoefficients));
  IMRPhenomXHMPhaseCoefficients *pPhase = (IMRPhenomXHMPhaseCoefficients*)XLALMalloc(sizeof(IMRPhenomXHMPhaseCoefficients));

  /* Initialize the PhenomXHM lm phase and amp coefficients struct */
  IMRPhenomXHM_FillAmpFitsArray(pAmp);
  IMRPhenomXHM_FillPhaseFitsArray(pPhase);

  /* Get coefficients for Amplitude and phase */
  GetSpheroidalCoefficients(pPhase, pPhase22, pWFHM, pWF);
  IMRPhenomXGetAmplitudeCoefficients(pWF, pAmp22);

  IMRPhenomXHM_GetAmplitudeCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF);
  IMRPhenomXHM_GetPhaseCoefficients(pAmp, pPhase, pAmp22, pPhase22, pWFHM, pWF,lalParams);

  dfmerger = deltaF_mergerBin(pWFHM->fDAMP, pPhase->alphaL_S, resTest);
  dfringdown = deltaF_ringdownBin(pWFHM->fDAMP, pPhase->alphaL_S, pAmp->lambda/(pAmp->sigma*pWFHM->fDAMP), resTest);

  #if DEBUG == 1
  printf("f_min = %.6f, Mfmin = %.6f\n", pWF->fMin, Mfmin);
  printf("f_max = %.6f, f_max_prime = %.6f, Mfmax = %.6f\n", pWF->fMax, pWF->f_max_prime, Mfmax);
  printf("\nfRING = %e\n",pWFHM->fRING);
  printf("fDAMP = %e\n",pWFHM->fDAMP);
  printf("\nMfmin = %.6f\n", Mfmin);
  printf("MfMECO = %.6f\n", MfMECO);
  printf("MfLorentzianEnd = %.6f\n", MfLorentzianEnd);
  printf("Mfmax = %.6f\n", Mfmax);
  printf("evaldMf = %.6e\n", evaldMf);
  printf("dfpower = %.6e\n", dfpower);
  printf("dfcoefficient = %.6e\n", dfcoefficient);
  printf("dfmerger = %.6e\n", dfmerger);
  printf("dfringdown = %.6e\n", dfringdown);
  printf("alphaL_S = %.16e\n", pPhase->alphaL_S);
  #endif

  /* Allocate memory for the list of grids. The number of grids must be less than lengthallGrids. */
  UINT4 lengthallGrids = 20;
  IMRPhenomXMultiBandingGridStruct *allGrids = (IMRPhenomXMultiBandingGridStruct*)XLALMalloc(lengthallGrids * sizeof(IMRPhenomXMultiBandingGridStruct));

  if (allGrids == NULL)
  {
    #if DEBUG == 1
    printf("Malloc of allGrids failed!\n");
    #endif
    return -1;
  }

  /* Compute the coarse frequency array. It is stored in a list of grids. */
  UINT4 nGridsUsed = XLALSimIMRPhenomXMultibandingGrid(Mfmin, MfMECO, MfLorentzianEnd, Mfmax, evaldMf, dfpower, dfcoefficient, allGrids, dfmerger, dfringdown);

  #if DEBUG == 1
  printf("allGrids[0].Length = %i\n", allGrids[0].Length);
  #endif

  /* Number of fine frequencies per coarse interval in every coarse grid */
  INT4 mydfRatio[lengthallGrids];
  /* Actual number of subgrids to be used. We allocated more than needed. */
  UINT4 actualnumberofGrids = 0;
  /* Length of coarse frequency array */
  UINT4 lenCoarseArray = 0;

  /* Transform the coarse frequency array to 1D array. */
  // Take only the subgrids needed
  for(UINT4 kk = 0; kk < nGridsUsed; kk++){
    lenCoarseArray = lenCoarseArray + allGrids[kk].Length;
    actualnumberofGrids++;

    mydfRatio[kk] = allGrids[kk].intdfRatio;

    #if DEBUG == 1
    printf("\nkk = %i\n",kk);
    printf("xStart: %.16e\n", allGrids[kk].xStart);
    printf("xEnd: %.16e\n", allGrids[kk].xEndRequested);
    printf("Length: %i\n", allGrids[kk].Length);
    printf("deltax: %.16e\n", allGrids[kk].deltax);
    printf("evaldMf: %.16e\n", evaldMf);
    printf("xMax: %.16e\n", allGrids[kk].xMax);
    printf("# fine points %i\n", mydfRatio[kk]);
    printf("Last grid.xMax = %.16f\n", allGrids[actualnumberofGrids-1].xMax);
    #endif

    if(allGrids[kk].xMax + evaldMf >= Mfmax){
      break;
    }
  }

  // Add extra points to the coarse grid if the last freq is lower than Mfmax
  while(allGrids[actualnumberofGrids-1].xMax < Mfmax){
    allGrids[actualnumberofGrids-1].xMax =   allGrids[actualnumberofGrids-1].xMax +  allGrids[actualnumberofGrids-1].deltax;
    allGrids[actualnumberofGrids-1].Length =  allGrids[actualnumberofGrids-1].Length + 1;
    lenCoarseArray++;
  }

  #if DEBUG == 1
  printf("\nfDAMP = %.16e\n", XLALSimIMRPhenomXUtilsMftoHz(pWFHM->fDAMP, pWF->Mtot));
  printf("\nactualnumberofGrids = %i\n", actualnumberofGrids);
  printf("lenCoarseArray = %i\n", lenCoarseArray);
  printf("Last grid.xMax = %.16f", allGrids[actualnumberofGrids-1].xMax);
  #endif

  // Transform coarse frequency array to 1D vector
  REAL8 *IntLawpoints = (REAL8*)XLALMalloc(lenCoarseArray * sizeof(REAL8));
  UINT4 lenIntLawpoints = 0;

  for(UINT4 kk = 0; kk < actualnumberofGrids; kk++){
    for(INT4 ll = 0; ll < allGrids[kk].Length; ll++){
      IntLawpoints[lenIntLawpoints] = (allGrids[kk].xStart + allGrids[kk].deltax*ll);
      lenIntLawpoints++;
    }
  }



  /* End of coarse frequency array. */

  #if DEBUG == 1
  printf("\n******** Coarse frequencies array done ********* \n");
  printf("\nlenIntLawpoints, coarse[0], coarse[-1], Mfmax, M_sec = %i %.16e %.16e %.16e %.16e\n",lenIntLawpoints, IntLawpoints[0], IntLawpoints[lenIntLawpoints-1], Mfmax, pWF->M_sec);
  #endif

  /* IntLawpoints stores the adimensional frequencies (Mf) */
  /* coarseFreqs will store the same frequency array but in Hz */
  /* Allocate memory for frequency array and terminate if this fails */
  REAL8Sequence *coarseFreqs;
  coarseFreqs = XLALCreateREAL8Sequence(lenCoarseArray);
  if (!coarseFreqs) { XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed."); }

  /* Populate frequency array */
  #if DEBUG == 1
  printf("\n***** Coarse freqs *****\n");
  #endif
  REAL8 divis = 1./pWF->M_sec;
  for (UINT4 ii = 0; ii < lenCoarseArray; ii++)
  {
    coarseFreqs->data[ii] = IntLawpoints[ii]*divis;
  }
  #if DEBUG == 1
  printf("\nFirst/Last Coarse freqs *****%.16f %.16f\n", coarseFreqs->data[0], coarseFreqs->data[coarseFreqs->length-1]);
  #endif

  /* Next we will compute amplitude and phase of one mode in the coarse frequency array. For that we need to split the frequency array in the spherical and spheroidal part. */
  /* We will compute 2 waveforms: one in the spherical part and another in the spheroidal. */

  #if DEBUG == 1
  printf("\n******* Splitting spherical/spheroidal coarseFreqs %i****************\n",lenCoarseArray);
  #endif


  REAL8 *IntLawpointsS = (REAL8*)XLALMalloc(lenCoarseArray*sizeof(REAL8));
  REAL8 *IntLawpointsSS = (REAL8*)XLALMalloc(lenCoarseArray*sizeof(REAL8));

  double MfRDcutMin, MfRDcutMax;
  if(pPhase->fPhaseMatchIM < pAmp->fAmpMatchIM){
    MfRDcutMin = pPhase->fPhaseMatchIM;
    MfRDcutMax = pAmp->fAmpMatchIM;
  }else{
    MfRDcutMin = pAmp->fAmpMatchIM;
    MfRDcutMax = pPhase->fPhaseMatchIM;
  }

  #if DEBUG == 1
  printf("\nMfRDcutMin = %.16f, MfRDcutMax = %.16f, lastcoarseArray = %.16f\n", MfRDcutMin, MfRDcutMax, IntLawpoints[lenCoarseArray-1]);
  #endif

  /* Compute the coarse frequencies in the spherical and in the spheroidal part. */

  unsigned int lencoarseS = 0, lencoarseSS = 0, enter = 0;
  for(UINT4 idx = 0; idx < lenCoarseArray; idx++){
    double ff =  IntLawpoints[idx];
    if(ff > MfRDcutMin){
      if(lencoarseSS < 1 && lencoarseS>0){ // This numbers tells how many points we add to spherical after MfRDcutMin and how many to spheroidal before MfRDcutMin
        IntLawpointsSS[lencoarseSS] = IntLawpointsS[lencoarseS-1];
        lencoarseSS++;

        if(idx == lenCoarseArray-1){
          IntLawpointsS[lencoarseS] = ff;
          lencoarseS++;
        }
      }
      IntLawpointsSS[lencoarseSS] = ff;
      lencoarseSS++;

      if(idx < lenCoarseArray-1 && IntLawpoints[idx+1] < MfRDcutMax){
        IntLawpointsS[lencoarseS] = ff;
        lencoarseS++;
      }
      if(idx < lenCoarseArray-1 && IntLawpoints[idx+1] > MfRDcutMax && enter<2){
        IntLawpointsS[lencoarseS] = ff;
        lencoarseS++;
        enter++;
      }
    }
    else{
      IntLawpointsS[lencoarseS] = ff;
      lencoarseS++;
    }
  }

  #if DEBUG == 1
  printf("\n******* Coarse Freqs Spherical/Spheroidal ****************\n");
  printf("%i ", lencoarseS);
  printf("%i ", lencoarseSS);
  #endif

  /* Transform spherical and spheroidal frequencies to Hz */
  REAL8Sequence *coarseFreqsS = XLALCreateREAL8Sequence(lencoarseS);
  REAL8Sequence *coarseFreqsSS = XLALCreateREAL8Sequence(lencoarseSS);
  for (UINT4 ii = 0; ii < lencoarseS; ii++)
  {
    coarseFreqsS->data[ii] = IntLawpointsS[ii]/pWF->M_sec;
  }
  for (UINT4 ii = 0; ii < lencoarseSS; ii++)
  {
    coarseFreqsSS->data[ii] = IntLawpointsSS[ii]/pWF->M_sec;
  }


  #if DEBUG == 1
  printf("\n******* Computing Coarse Phase and Amp ****************\n");
  printf("%i ", coarseFreqsS->length);
  printf("%i ", coarseFreqsSS->length);
  #endif

  /* Compute coarse amplitude and phase in the spherical and spheroidal part */
  /* Declare amplitude and phase variables */
  REAL8FrequencySeries *amplitude, *phase;
  REAL8FrequencySeries *phaseSS, *amplitudeSS;


  if(lencoarseS > ampinterpolorder){
    IMRPhenomXHM_Phase(&phase, coarseFreqsS, pWF, pAmp22, pPhase22, pWFHM, pAmp, pPhase);
    IMRPhenomXHM_Amplitude(&amplitude, coarseFreqsS, pWF, pAmp22, pPhase22, pWFHM, pAmp, pPhase);
  }
  if(lencoarseSS > ampinterpolorder){
    IMRPhenomXHM_PhaseMixing(&phaseSS, coarseFreqsSS, pWF, pWFHM, pPhase);
    IMRPhenomXHM_AmplitudeMixing(&amplitudeSS, coarseFreqsSS, pWF, pWFHM, pAmp, pPhase);
    #if DEBUG == 1
    printf("\nLength@phaseSS = %i\n",phaseSS->data->length);
    #endif
  }



  #if DEBUG == 1
  printf("\n******* Computed Coarse Amp and Phase **************** %i\n", coarseFreqs->length);
  #endif

  /* Transform the REAL8FrequencySeries to vector since gsl uses vectors for the interpolation.
     gsl is only used for the amplitude but to keep the code more symmetric we transform also the phase. */
  REAL8 *ILamplm = (REAL8*)XLALMalloc(lencoarseS * sizeof(REAL8));
  REAL8 *ILphaselm = (REAL8*)XLALMalloc(lencoarseS * sizeof(REAL8));
  REAL8 *ILamplmSS = (REAL8*)XLALMalloc(lencoarseSS * sizeof(REAL8));
  REAL8 *ILphaselmSS = (REAL8*)XLALMalloc(lencoarseSS * sizeof(REAL8));

  if(lencoarseS > ampinterpolorder){
    for (UINT4 ii = 0; ii < lencoarseS; ii++)
    {
      ILamplm[ii] = ((amplitude)->data->data)[ii];
      ILphaselm[ii] = ((phase)->data->data)[ii];
    }
  }
  if(lencoarseSS > ampinterpolorder){
    for (UINT4 ii = 0; ii < lencoarseSS; ii++)
    {
      ILamplmSS[ii] = ((amplitudeSS)->data->data)[ii];
      ILphaselmSS[ii] = ((phaseSS)->data->data)[ii];
    }
  }

  #if DEBUG == 1
    if(lencoarseS > ampinterpolorder) printf("\nLast spherical phases %.16f %.16f", ILphaselm[lencoarseS-2], ILphaselm[lencoarseS-1]);
    if(lencoarseSS > ampinterpolorder) printf("\nLast spherical freqs %.16f %.16f", IntLawpoints[lencoarseS-2], IntLawpoints[lencoarseS-1]);
  #endif

  #if DEBUG == 1
  //Save coarse amplitude and phase (spherical and spheroidal) in file
  FILE *file0;
  char fileSpec0[40];
  sprintf(fileSpec0, "coarseamplitude%i%i.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec0);
  file0 = fopen(fileSpec0,"w");
  fprintf(file0,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->q, pWF->chi1L, pWF->chi2L, ell, emm);
  FILE *file3;
  char fileSpec3[40];
  sprintf(fileSpec3, "coarsephase%i%i.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec3);
  file3 = fopen(fileSpec3,"w");
  fprintf(file3,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->q, pWF->chi1L, pWF->chi2L, ell, emm);
  for(UINT4 idx = 0; idx < (UINT4)lencoarseS && lencoarseS>ampinterpolorder; idx++)
  {
    fprintf(file0, "%.16f  %.16e \n",  coarseFreqsS->data[idx], ILamplm[idx]);
    fprintf(file3, "%.16f  %.16e \n",  coarseFreqsS->data[idx]*pWF->M_sec, ILphaselm[idx]);
  }
  fprintf(file3, "\n\n");
  if(lencoarseSS > ampinterpolorder){
    for(UINT4 idx = 0; idx < (UINT4)lencoarseSS && lencoarseSS>ampinterpolorder; idx++)
    {
      fprintf(file0, "%.16f  %.16e \n",  coarseFreqsSS->data[idx], ILamplmSS[idx]);
      fprintf(file3, "%.16f  %.16e \n",  coarseFreqsSS->data[idx]*pWF->M_sec, ILphaselmSS[idx]);
    }
  }
  fclose(file0);
  fclose(file3);
  #endif

  /* Free allocated memory */
  if(lencoarseS > ampinterpolorder){
    XLALDestroyREAL8FrequencySeries(amplitude);
    XLALDestroyREAL8FrequencySeries(phase);
  }
  if(lencoarseSS > ampinterpolorder){
    XLALDestroyREAL8FrequencySeries(phaseSS);
    XLALDestroyREAL8FrequencySeries(amplitudeSS);
  }
  XLALDestroyREAL8Sequence(coarseFreqs);
  XLALDestroyREAL8Sequence(coarseFreqsS);
  XLALDestroyREAL8Sequence(coarseFreqsSS);


  /***** Linear Interpolation of e^(I*phi) with the iterative procedure *****/

  /* Estimation of the length of the fineGrid */
  INT4 lenWF = 0;

  for(UINT4 kk = 0; kk < actualnumberofGrids; kk++){
    lenWF = lenWF + (allGrids[kk].Length -1) * mydfRatio[kk] + 2*mydfRatio[kk];
    #if DEBUG == 1
    printf("\nmydfRatio[%i] = %i %i %i", kk, mydfRatio[kk],allGrids[kk].Length, (allGrids[kk].Length -1) * mydfRatio[kk] + 2*mydfRatio[kk]);
    printf("\nlenWF = %i", lenWF);
    #endif
  }

  // Variable to control the number of coarse points we have used
  int pointsPrecessedSoFar = 0;

  /* Allocate memory for the complex exponential and the fine equally-spaced frequency array */
  COMPLEX16 *expphi = (COMPLEX16*)XLALMalloc(lenWF*sizeof(COMPLEX16));
  if (expphi == NULL) {
    return -1;
  }

  REAL8 *finefreqs = (REAL8*)XLALMalloc(lenWF*sizeof(REAL8));
  if (finefreqs == NULL){
    return -1;
  }

  UINT4 count = 0;        // Variable to track the point being filled in the fine frequency grid
  UINT4 coarsecount = 0;  // Variable to track the point used in the coarse frequency grid
  UINT4 RDcutMin = 0, RDcutMax = 0;  // Variables to know where to separate the fine spherical and spheroidal part. RDcutMin is for the phase and RDcutMax for the amplitude.
  COMPLEX16 Q32 = 0;
  INT4 lenRD = round((Mfmax - MfRDcutMin)/evaldMf) + 3; //+3 to be safe and no run out of memory because of rounding erros
  if(lenRD <= 0) {
    lenRD = 1;
  }

  /* After the rotation from spheroidal to spherical we will have to add a linear part to the phase of the 32. This part is also computed using Multibanding. */
  COMPLEX16 *linear32 = (COMPLEX16*)XLALMalloc( lenRD * sizeof(COMPLEX16));

  /* Loop over allgrids */
  bool stop = false;
  for (UINT4 i = 0; i<actualnumberofGrids && !stop; i++){
    #if DEBUG == 1
    printf("\ni = %i\n", i);
    #endif


    /* Compute mydfRatio */
    UINT4 lcoarseGrid = allGrids[i].Length;
    if(lcoarseGrid == 0){
      break;
    }

    /* Linear Interpolation and iteration */
    if(i==actualnumberofGrids-1){
      lcoarseGrid--;
    }

    REAL8 Omega, phi0, Mfhere = 0, Mfnext=0;
    COMPLEX16 h0, Q;

    /* Loop over the coarse points of a subgrid */
    for(UINT4 j = 0; j < lcoarseGrid && !stop; j++){
      Mfhere = IntLawpoints[pointsPrecessedSoFar + j];
      Mfnext = IntLawpoints[pointsPrecessedSoFar + j + 1] ;

      INT4 ratio;
      // If we are in the last coarse point of a sublist we use the number of fine points of the next subgrid.
      if(j==lcoarseGrid-1 && i < actualnumberofGrids-1){
        ratio = mydfRatio[i+1];
      }
      else{
        ratio = mydfRatio[i];
      }
      if(Mfnext + evaldMf >= Mfmax){
        double dratio = (Mfmax - Mfhere)/evaldMf + 1;
        ratio = (int) dratio;
        int roundratio = round((Mfmax - Mfhere)/evaldMf) + 1;
        if(fabs(dratio-roundratio) < 0.0001) ratio = roundratio;  // To get the correct rounded integer
        /* Break the loop if you overpass the maximun frequency */
        stop = true;
        #if DEBUG == 1
        printf("\nMfmax, Mfhere, evaldMf, ratio= %.16f %.16f %.16f %i\n", Mfmax, Mfhere, evaldMf, ratio);
        #endif
      }

      /********************************************/
      /*     Inner Loop: linear interpolation     */
      /********************************************/
      /****  Spheroidal part ****/
      if(Mfhere > MfRDcutMin && lencoarseSS > ampinterpolorder) {
        if(RDcutMin == 0 ){
          MfRDcutMin = Mfhere;
          RDcutMin = count;
          linear32[0] = cexp( I * (pPhase->C1RD*Mfhere+pPhase->CRD + pPhase->deltaphiLM));
          Q32 = cexp( I * evaldMf * (pPhase->C1RD));
          #if DEBUG == 1
              printf("\n*** Starting spheroidal part ****\n");
              printf("Mfhere, MfRDcutMin, RDcutMin = %.16e %.16e %i\n", Mfhere, MfRDcutMin, RDcutMin);
          #endif
          if(lencoarseS <= ampinterpolorder){
             coarsecount++; //When the coarse array does not have spherical part.
          }
        }

        INT4 jdx = pointsPrecessedSoFar + j - coarsecount + 1;

        phi0 = ILphaselmSS[jdx];

        if(j + pointsPrecessedSoFar < lenCoarseArray){
          Omega = (ILphaselmSS[ jdx + 1] - ILphaselmSS[jdx])/(IntLawpoints[pointsPrecessedSoFar + j + 1] - IntLawpoints[pointsPrecessedSoFar + j]);
        }
        else{
          Omega = (ILphaselmSS[jdx] - ILphaselmSS[jdx - 1])/(IntLawpoints[pointsPrecessedSoFar + j] - IntLawpoints[pointsPrecessedSoFar + j -1]);
        }

        if(Mfhere > pAmp->fAmpMatchIM && RDcutMax < 1){
          RDcutMax = count;
        }
      }
      /****  Spherical part ****/
      else{
        UINT4 jjdx = j + pointsPrecessedSoFar;
        if(jjdx < lenCoarseArray){
          Omega = (ILphaselm[jjdx+ 1] - ILphaselm[jjdx])/(IntLawpoints[jjdx + 1] - IntLawpoints[jjdx]);
        }
        else{
          Omega = (ILphaselm[jjdx] - ILphaselm[jjdx -1])/(IntLawpoints[jjdx] - IntLawpoints[jjdx -1]);

        }
        phi0 = ILphaselm[jjdx];
        coarsecount++;
      }

      h0   = cexp(I*phi0);
      Q    = cexp(I*evaldMf*Omega);

      if(RDcutMin == 0 && lencoarseS>ampinterpolorder){ // Spherical part
        finefreqs[count] = Mfhere;
        expphi[count]    = h0;
        count++;
        /* This loop carry out the eq. 2.32 in arXiv:2001.10897 */
        for(int kk = 1; kk < ratio; kk++){       // Compute finefreqs and fine expphi
          finefreqs[count] = Mfhere + evaldMf*kk;
          expphi[count]    = Q*expphi[count-1];
          count++;
        }
      }
      else{   // Spheroidal part
        finefreqs[count] = Mfhere;
        expphi[count]    = h0;
        linear32[count - RDcutMin +1] = Q32*linear32[count - RDcutMin ];
        count++;
        /* This loop carry out the eq. 2.32 in arXiv:2001.10897 but for spheroidals */
        for(int kk = 1; kk < ratio; kk++){        // Compute finefreqs and fine expphi
          finefreqs[count] = Mfhere + evaldMf*kk;
          expphi[count]    = Q*expphi[count-1];
          linear32[count - RDcutMin +1] = Q32*linear32[count - RDcutMin ];
          count++;
        }
      }
    }
    pointsPrecessedSoFar = pointsPrecessedSoFar + lcoarseGrid;
  }// End loop over ILgrids. count should be aprox = to lenWF

  if(RDcutMin == 0){
    RDcutMin = count;
  }
  if(lencoarseS <= ampinterpolorder){
    RDcutMin = 0;
  }
  if(RDcutMax == 0){
    RDcutMax = count;
  }


  #if DEBUG == 1
  printf("\nTheory and practice (should not be equal) %i %i \n", lencoarseS, coarsecount );
  printf("\n******* Interpolate Amplitude **********\n");
  #endif

  /*** Interpolation and evaluation in fineFreqs of the amplitude ***/
  REAL8 *fineAmp = (REAL8*)XLALMalloc(count * sizeof(REAL8));
  REAL8 *fineAmpSS = (REAL8*)XLALMalloc(count * sizeof(REAL8));

  interpolateAmplitudeMixing(fineAmp, fineAmpSS, IntLawpointsS, IntLawpointsSS, ILamplm, ILamplmSS, finefreqs, lencoarseS, lencoarseSS, count, RDcutMin, RDcutMax, ampinterpolorder);


  /**** Build the waveform ****/
  size_t offset = iStart;

  #if DEBUG ==1
  printf("\nfinefreqs[0], evaldMf, quotient, offset = %.16e %.16e %.16e %zu\n", finefreqs[0], evaldMf, finefreqs[0]/evaldMf, offset);
  #endif

  // Due to round of erros, the last freq may be greater Mfmax. The difference should be less than the frequency step.
  // Remove that extra frequency point
  while(finefreqs[count-1] > Mfmax){
    count--;
    if(RDcutMin>count) RDcutMin = count;
    if(RDcutMax>count) RDcutMax = count;
  }

  /* Intialize FrequencySeries for htildelm */
  while(count + offset > iStop){ count--;}
  size_t n = iStop;
  XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / deltaF ), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.", deltaF);
  //XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / deltaF + 500*pWF->Mtot*LAL_MTSUN_SI), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df + 500M with df = %g.", deltaF);
  *htildelm = XLALCreateCOMPLEX16FrequencySeries("htildelm: FD waveform", &(ligotimegps_zero), 0.0, deltaF, &lalStrainUnit, n);
  for(int idx = 0; idx < (int) offset; idx++){
    ((*htildelm)->data->data)[idx] = 0.;
  }

  /* Sometimes rounding error can make that count+offset < iStop, meaning that we have computed less points than those we are gonna output,
    here we make sure that all the elements of htildelm are initialized and we put the extra point(s) to 0. */
  for(UINT4 idx = count-1; idx < iStop-offset; idx++){
    /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
    ((*htildelm)->data->data)[idx + offset] = 0.;
  }

  #if DEBUG == 1
  printf("\n**** pPhase->fPhaseMatchIM, pPhase->fAmpMatchIM, MfRDcutMin = %.16f %.16f %.16f \n",pPhase->fPhaseMatchIM, pAmp->fAmpMatchIM, MfRDcutMin);
  printf("\ncount RDcutMin RDcutMax: %i %i %i\n", count, RDcutMin, RDcutMax);
  printf("\nMfRDcutMin MfRDcutMax: %.6f %.6f\n", MfRDcutMin, MfRDcutMax);
  //Save lm in file
  FILE *file5;
  char fileSpec5[40];
  sprintf(fileSpec5, "sphericalfine%i%i.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec5);
  file5 = fopen(fileSpec5,"w");
  fprintf(file5,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->q, pWF->chi1L, pWF->chi2L, ell, emm);
  fprintf(file5,"# Frequency (Hz)   Fine Spherical Phase \n");
  //Save lm in file
  FILE *file6;
  char fileSpec6[40];
  sprintf(fileSpec6, "spheroidalfine%i%i.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec6);
  file6 = fopen(fileSpec6,"w");
  fprintf(file6,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->q, pWF->chi1L, pWF->chi2L, ell, emm);
  fprintf(file6,"# Frequency (Hz)   Fine Spheroidal Phase \n");
  #endif

  /* (-1)^l */
  INT4 minus1l;
  if (ell % 2 != 0)
  minus1l = -1;
  else
  minus1l = +1;

  /******** Spherical part **********/
  #if DEBUG == 1
  printf("\nRDcutMin = %i", RDcutMin);
  #endif
  for(UINT4 idx = 0; idx < RDcutMin; idx++){
    COMPLEX16 data5 =  minus1l * fineAmp[idx] * expphi[idx];
    ((*htildelm)->data->data)[idx + offset] = data5;
    #if DEBUG == 1
    fprintf(file5, "%.16f  %.16e %.16e\n",(idx + offset)*deltaF, creal(data5), cimag(data5));
    #endif
  }
  #if DEBUG ==1
  fclose(file5);
  #endif


  /********** Rotate Spheroidal part ******************/
  // Compute htilde22 for the ringdown part if the 22 mode was not computed  before.
  COMPLEX16FrequencySeries *htilde22tmp = NULL;

  if(htilde22 == NULL && count>RDcutMin){
    REAL8Sequence *freqs = XLALCreateREAL8Sequence(count - RDcutMin);
    #if DEBUG == 1
    printf("\nBuilding RD 22 mode\n");
    printf("\nlen@freqs, RDcutMin, count = %i %i %i\n", freqs->length, RDcutMin, count);
    #endif
    for(UINT4 idx = RDcutMin; idx < count; idx++){
      freqs->data[idx-RDcutMin] = finefreqs[idx]/pWF->M_sec;
    }
    XLALSimIMRPhenomXASFrequencySequence(&htilde22tmp, freqs, pWF->m1_SI, pWF->m2_SI, pWF->chi1L, pWF->chi2L, pWF->distance, pWF->phiRef_In, pWF->fRef, lalParams);
    XLALDestroyREAL8Sequence(freqs);
  }
  else{
    htilde22tmp = XLALCreateCOMPLEX16FrequencySeries("htilde22: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,count-RDcutMin);
    for(UINT4 idx = RDcutMin; idx < count; idx++){
      htilde22tmp->data->data[idx - RDcutMin] = htilde22->data->data[idx + offset];
    }
  }



  COMPLEX16 shift32 = 0;

  #if DEBUG == 1
  printf("\nLength@htilde22, count, RDcutMin, RDcutMax, offset %i %i %i %i %zu\n",htilde22tmp->data->length, (int)(count+offset), RDcutMin, RDcutMax, offset);
  #endif
  for(UINT4 idx = RDcutMin; idx < count; idx++){
    double Mf = finefreqs[idx];
    IMRPhenomX_UsefulPowers powers_of_f;
    IMRPhenomX_Initialize_Powers(&powers_of_f, Mf);
    //22 waveform complete (without rescaling)
    //COMPLEX16 wf22 = htilde22->data->data[idx + offset];
    COMPLEX16 wf22 = htilde22tmp->data->data[idx - RDcutMin];
    //32 waveform in spheroidal
    REAL8 amplm = fineAmpSS[idx-RDcutMin] * (powers_of_f.m_seven_sixths*pWFHM->Amp0);
    COMPLEX16 expphilm = expphi[idx];
    //Rotation to spherical with Berti's coefficients
    shift32 = linear32[idx-RDcutMin];
    COMPLEX16 sphericalWF_32 = (conj(pWFHM->mixingCoeffs[2])*wf22 + conj(pWFHM->mixingCoeffs[3])*amplm*expphilm)*shift32;
    COMPLEX16 data;

    // Use spheroidal phase but spherical amplitude
    if(Mf < pAmp->fAmpMatchIM){
      data = sphericalWF_32/cabs(sphericalWF_32) * fineAmp[idx];
      ((*htildelm)->data->data)[idx + offset] = data * minus1l;
    }
    // Use spheroidal amplitude and phase
    else{
      data = sphericalWF_32;
      ((*htildelm)->data->data)[idx + offset] = data * minus1l;
    }
    #if DEBUG == 1
    data = ((*htildelm)->data->data)[idx + offset];
    fprintf(file6, "%.16f  %.16e %.16e\n",(idx + offset)*deltaF, creal(data), cimag(data));
    #endif
  }
  #if DEBUG == 1
  fclose(file6);
  #endif

  /* Free allocated memory */
  XLALDestroyCOMPLEX16FrequencySeries(htilde22tmp);
  LALFree(pAmp);
  LALFree(pAmp22);
  LALFree(pPhase);
  LALFree(pPhase22);
  LALFree(linear32);


  #if DEBUG == 1
  printf("\n******* Building Waveform 2**********\n");
  printf("\nfinefreqs[0]Hz = %.16f\n", finefreqs[0]/pWF->M_sec);
  printf("\nfinefreqs[%i]Hz = %.16f\n", count-1, finefreqs[count-1]/pWF->M_sec);
  printf("Mfmax in Hz = %.16f\n", Mfmax/pWF->M_sec);
  printf("count-1, offset = %i %zu\n", count-1, offset);
  printf("Mfmaxtheory = %.16f\n", (count-1+offset)*deltaF);
  printf("\nlength@htildelm, count+offset = %i %zu", (*htildelm)->data->length, count+offset);
  printf("\nf_max_prime = %.16e", pWF->f_max_prime);
  printf("\nlastfreq true = %.16e", (*htildelm)->data->length * deltaF);
  printf("\nlastfreq true = %.16e %.16e", creal((*htildelm)->data->data[count-1]), cimag((*htildelm)->data->data[count-1]));
  #endif


  #if DEBUG == 1
  //Save hlm mode in file
  FILE *file;
  char fileSpec[40];
  sprintf(fileSpec, "simulation%i%i_Multiband.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec);
  file = fopen(fileSpec,"w");
  fprintf(file,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->m1_SI/pWF->m2_SI, pWF->chi1L, pWF->chi2L, ell, emm);
  fprintf(file,"# Frequency (Hz)    Real    Imaginary\n");
  COMPLEX16 data;
  for(UINT4 idx = 0; idx < ((*htildelm)->data->length); idx++)
  {
    data = ((*htildelm)->data->data)[idx];
    fprintf(file, "%.16f  %.16e %.16e\n",  idx*deltaF, creal(data), cimag(data));
  }
  fclose(file);
  #endif

  #if DEBUG == 1
  //Save fine amplitude and freq array in file
  FILE *file2;
  char fileSpec2[40];
  sprintf(fileSpec2, "amplitude%i%i_fine.dat", ell,emm);
  printf("\nOutput file: %s\r\n",fileSpec2);
  file2 = fopen(fileSpec2,"w");
  fprintf(file2,"# q = %.16f chi1 = %.16f chi2 = %.16f lm = %i%i\n", pWF->m1_SI/pWF->m2_SI, pWF->chi1L, pWF->chi2L, ell, emm);
  for(UINT4 idx = 0; idx < RDcutMax && lencoarseS>ampinterpolorder; idx++)
  {
	  fprintf(file2, "%.16f  %.16e\n",  finefreqs[idx], fineAmp[idx]);
  }
  fclose(file2);
  #endif

  /* Free allocated memory */
  LALFree(expphi);
  LALFree(finefreqs);
  LALFree(allGrids);
  LALFree(fineAmp);
  LALFree(fineAmpSS);
  LALFree(pWFHM);
  LALFree(ILamplm);
  LALFree(ILamplmSS);
  LALFree(ILphaselm);
  LALFree(ILphaselmSS);
  LALFree(IntLawpoints);
  LALFree(IntLawpointsS);
  LALFree(IntLawpointsSS);

  return XLAL_SUCCESS;
}


/**************************************/
/*      INTERPOLATING FUNCTIONS       */
/**************************************/

/* Interpolate the 2D array [coarseFreqs, coarseAmp] and evaluate in finefreqs */
/* Interpolation uses third order */
static int interpolateAmplitude(
  double *fineAmp,      /**<[out] amplitude in the fine uniform grid **/
  double coarsefreqs[], /**< non-uniform frequency array**/
  double coarseAmp[],   /**< amplitude in the non-uniform frequency array **/
  double finefreqs[],   /**< uniform fine frequency grid**/
  int lengthCoarse,     /**< length of non-uniform freq array **/
  int lengthFine,       /**< length of uniform fine freq array **/
  int ampinterpolorder     /**< order of the gsl interpolation **/
){

  #if DEBUG == 1
  printf("\n****Building interpolants*****\n");
  printf("Number of points to interpolate = %i\r\n", lengthCoarse);
  #endif

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline;
  switch(ampinterpolorder){
    case 1:{
      spline = gsl_spline_alloc(gsl_interp_linear, lengthCoarse);
      break;
    }
    case 3:{
      spline = gsl_spline_alloc(gsl_interp_cspline, lengthCoarse);
      break;
    }
    default:{spline = gsl_spline_alloc(gsl_interp_cspline, lengthCoarse);}
  }
  gsl_spline_init(spline, coarsefreqs, coarseAmp, lengthCoarse);

  #if DEBUG == 1
  printf("\n****Loop for fine freqs*****\n");
  #endif

  for(INT4 kk = 0; kk < lengthFine; kk++){
    if(finefreqs[kk] < coarsefreqs[0] || finefreqs[kk] > coarsefreqs[lengthCoarse-1]){
      #if DEBUG == 1
      printf("\nOut of coarse range: coarse[0], coarse[-1] fine[%i] %.16e %.16e %.16e\n", kk, coarsefreqs[0], coarsefreqs[lengthCoarse-1], finefreqs[kk]);
      #endif
      fineAmp[kk] = fineAmp[kk-1];
    }
    else{
      fineAmp[kk] = gsl_spline_eval(spline, finefreqs[kk], acc);
    }
  }
  #if DEBUG == 1
  printf("\n****Free memory*****\n");
  #endif
  /* Free memory */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return 0;

}

/* Do 2 interpolations, one in the spherical regime and other in the spheroidal. */
static int interpolateAmplitudeMixing(
  double *fineAmp,            /**<[out] spherical amplitude in the fine uniform grid **/
  double *fineAmpSS,          /**<[out] spheroidal amplitude in the fine uniform grid **/
  double coarsefreqs[],       /**< non-uniform frequency array**/
  double coarsefreqsSS[],     /**< non-uniform frequency array spheroidal **/
  double coarseAmp[],         /**< amplitude in the non-uniform frequency array **/
  double coarseAmpSS[],       /**< spheroidal amplitude in the non-uniform frequency array **/
  double finefreqs[],         /**< uniform fine frequency grid**/
  int lengthCoarse,           /**< length of non-uniform freq array **/
  int lengthCoarseSS,         /**< length of non-uniform freq array Sphroidal **/
  int lengthFine,             /**< length of uniform fine freq array **/
  int sphericalfinecount,     /**< length of spherical fine grid **/
  int sphericalfinecountMax,  /**< length of spherical fine grid **/
  int ampinterpolorder           /**< order of interpolation **/
)
{
  int spheroidalfinecount = lengthFine - sphericalfinecount;

  REAL8 *sphericalfineFreqs = (REAL8*)XLALMalloc(sphericalfinecountMax * sizeof(REAL8));
  REAL8 *spheroidalfineFreqs = (REAL8*)XLALMalloc(spheroidalfinecount * sizeof(REAL8));

  for(int i = 0; i < sphericalfinecountMax; i++){
    sphericalfineFreqs[i] = finefreqs[i];
  }
  for(int i = 0; i < spheroidalfinecount; i++){
    spheroidalfineFreqs[i] = finefreqs[i + sphericalfinecount];
  }

  if(lengthCoarse > ampinterpolorder) interpolateAmplitude(fineAmp, coarsefreqs, coarseAmp, sphericalfineFreqs, lengthCoarse, sphericalfinecountMax, ampinterpolorder);
  if(lengthCoarseSS > ampinterpolorder) interpolateAmplitude(fineAmpSS, coarsefreqsSS, coarseAmpSS, spheroidalfineFreqs,  lengthCoarseSS, spheroidalfinecount, ampinterpolorder);

  #if DEBUG == 1
  //Save lm in file
  FILE *file2;
  char fileSpec2[40];
  sprintf(fileSpec2, "amplitude%i%i_fineSS.dat", 3,2);
  printf("\nOutput file: %s\r\n",fileSpec2);
  file2 = fopen(fileSpec2,"w");
  REAL8 data2;
  for(long int idx = 0; idx < spheroidalfinecount; idx++)
  {
    data2 = fineAmpSS[idx]; // uninitialised value
    fprintf(file2, "%.16f  %.16e\n",  finefreqs[ sphericalfinecount + idx], data2);
  }
  fclose(file2);
  #endif

  LALFree(sphericalfineFreqs);
  LALFree(spheroidalfineFreqs);

  return 0;

}


/****************************************/
/*                                      */
/*         AUXILIARY FUNCTIONS          */
/*                                      */
/****************************************/

/* Set up frequency array and frequency series for amplitude or phase */
// We pass it the coarse frequency array so it just initialize amplitude or phase.
static int SetupWFArraysReal(
  REAL8Sequence **freqs,           /**<[out] Frequency array to evaluate model **/
  REAL8FrequencySeries **amphase,  /**<[out] Initialize amplitude or phase with the length of freqs **/
  REAL8Sequence *freqs_In,         /**< Input frequency array or fmin, fmax **/
  IMRPhenomXWaveformStruct *pWF,   /**< Structure of the 22 mode **/
  LIGOTimeGPS ligotimegps_zero     /**< Needed to initialize amphase **/
){

  /* Inherit minimum and maximum frequencies to generate wavefom from input frequency grid */
  double f_min = freqs_In->data[0];
  double f_max = freqs_In->data[freqs_In->length - 1];

  /* Size of array */
  size_t npts     = 0;

  /* Index shift between freqs and the frequency series */
  UNUSED UINT4 offset    = 0;

  #if DEBUG == 1
  printf("f_min, f_max = %.6f %.6f \n",f_min,f_max);
  #endif

  /* If deltaF is non-zero then we need to generate a uniformly sampled frequency grid of spacing deltaF. Start at f = 0. */
  if(pWF->deltaF > 0)
  {
    #if DEBUG == 1
    printf("\n******* deltaF > 0 ************\n");
    #endif
    /* Return the closest power of 2 */
    npts = NextPow2(f_max / pWF->deltaF) + 1;

    /* Debug information */
    if(pWF->debug)
    {
      #if DEBUG == 1
      printf("npts     = %zu\n",npts);
      printf("fMin     = %.4f\n",f_min);
      printf("fMax     = %.4f\n",f_max);
      printf("dF       = %.4f\n",pWF->deltaF);
      #endif
    }

    /* Coalescence time is fixed to t=0, shift by overall length in time.  */
    XLAL_CHECK(XLALGPSAdd(&ligotimegps_zero, -1. / pWF->deltaF), XLAL_EFUNC, "Failed to shift the coalescence time to t=0. Tried to apply a shift of -1/df with df = %g.",pWF->deltaF);
    #if DEBUG == 1
    printf("f_min, f_max = %.6f %.6f \n",f_min,f_max);
    #endif
    /* Initialize the htilde frequency series */
    *amphase = XLALCreateREAL8FrequencySeries("amphase: FD waveform",&ligotimegps_zero,0.0,pWF->deltaF,&lalStrainUnit,npts);
    /* Check that frequency series generated okay */
    XLAL_CHECK(*amphase,XLAL_ENOMEM,"Failed to allocate REAL8FrequencySeries of length %zu for f_max = %f, deltaF = %g.\n",npts,f_max,pWF->deltaF);

    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t) (f_min / pWF->deltaF);
    size_t iStop  = (size_t) (f_max / pWF->deltaF) + 1;

    XLAL_CHECK ( (iStop <= npts) && (iStart <= iStop), XLAL_EDOM,
    "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);
    #if DEBUG == 1
    printf("f_min, f_max = %.6f %.6f \n",f_min,f_max);
    #endif
    /* Allocate memory for frequency array and terminate if this fails */
    (*freqs) = XLALCreateREAL8Sequence(iStop - iStart);
    #if DEBUG == 1
    printf("f_min, f_max = %.6f %.6f \n",f_min,f_max);
    #endif
    if (!(*freqs))
    {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
    #if DEBUG == 1
    printf("f_min, f_max = %.6f %.6f \n",f_min,f_max);
    #endif
    /* Populate frequency array */
    for (UINT4 i = iStart; i < iStop; i++)
    {
      (*freqs)->data[i-iStart] = i * pWF->deltaF;
    }
    offset = iStart;
  }
  else
  {
    #if DEBUG == 1
    printf("\n******* deltaF = 0 ************\n");
    #endif
    /* freqs is a frequency grid with non-uniform spacing, so we start at the lowest given frequency */
    npts      = freqs_In->length;
    *amphase = XLALCreateREAL8FrequencySeries("amphase: FD waveform, 22 mode", &ligotimegps_zero, f_min, pWF->deltaF, &lalStrainUnit, npts);

    XLAL_CHECK (*amphase, XLAL_ENOMEM, "Failed to allocated waveform REAL8FrequencySeries of length %zu from sequence.", npts);

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
  }//end freqs
  memset((*amphase)->data->data, 0, npts * sizeof(REAL8));
  XLALUnitMultiply(&((*amphase)->sampleUnits), &((*amphase)->sampleUnits), &lalSecondUnit);

  return offset;
}


/** Functions to compute coarse amplitude and phase **/

/* Spherical amplitude evaluated in an input frequency array */
static int IMRPhenomXHM_Amplitude(
  REAL8FrequencySeries **amplm,           /**<[out] amplitude of hlm mode **/
  REAL8Sequence *freqs_In,                /**< Frequency array to evaluate model or fmin, fmax  **/
  IMRPhenomXWaveformStruct *pWF,          /**< Structure of the 22 mode **/
  IMRPhenomXAmpCoefficients *pAmp22,      /**< Amplitude coefficients 22 */
  IMRPhenomXPhaseCoefficients *pPhase22,  /**< Phase coefficients 22 */
  IMRPhenomXHMWaveformStruct *pWFHM,      /**< waveform parameters lm mode */
  IMRPhenomXHMAmpCoefficients *pAmp,      /**< Amplitude coefficients lm */
  IMRPhenomXHMPhaseCoefficients *pPhase  /**< Phase coefficients 22 */
)
{

  #if DEBUG == 1
  printf("\n **** IMRPhenomXHM_Amplitude **** \n");
  printf("\nf_min, f_max = %.16e %.16e\n", freqs_In->data[0], freqs_In->data[freqs_In->length-1]);
  #endif

  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  int status = 0;
  REAL8Sequence *freqs;
  UINT4 offset = SetupWFArraysReal(&freqs, amplm, freqs_In, pWF, ligotimegps_zero);

  #if DEBUG == 1
  printf("\n***Length@freqs, offset %i %i",freqs_In->length,offset);
  printf("\n\nfstart, fend = %.16f %.16f\n\n", freqs_In->data[0], freqs_In->data[freqs_In->length-1]);
  printf("\n***Length@freqs, offset %i %i",freqs->length,offset);
  printf("\n\nfstart, fend = %.16f %.16f\n\n", freqs->data[0], freqs->data[freqs->length-1]);
  #endif

  if(pWFHM->Ampzero==0){
    IMRPhenomX_UsefulPowers powers_of_Mf;
    UINT4 initial_status = XLAL_SUCCESS;
    REAL8 Msec = pWF->M_sec;
    REAL8 amp;

    /* Loop over frequencies to generate waveform */
    if(pWFHM->MixingOn==1){
      for (UINT4 idx = 0; idx < freqs->length; idx++)
      {
        REAL8 Mf = Msec * freqs->data[idx];
        initial_status = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
        if(initial_status != XLAL_SUCCESS)
        {
          status = initial_status;
          XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
        }
        else
        {
          amp = IMRPhenomXHM_Amplitude_ModeMixing(Mf, &powers_of_Mf, pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF);
          /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
          ((*amplm)->data->data)[idx+offset] = pWFHM->Amp0 * amp;
        }
      }
    }
    else{
      for (UINT4 idx = 0; idx < freqs->length; idx++)
      {
        REAL8 Mf = Msec * freqs->data[idx];
        initial_status = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
        if(initial_status != XLAL_SUCCESS)
        {
          status = initial_status;
          XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
        }
        else
        {
          amp = IMRPhenomXHM_Amplitude_noModeMixing(Mf, &powers_of_Mf, pAmp, pWFHM);
          /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
          ((*amplm)->data->data)[idx+offset] = pWFHM->Amp0 * amp;
        }
      }
    }
  }
  /* Free allocated memory */
  XLALDestroyREAL8Sequence(freqs);

  return status;
}


/* Ringdown amplitude ansatz evaluated in an input frequency array */
static int IMRPhenomXHM_AmplitudeMixing(
  REAL8FrequencySeries **amplm,           /**<[out] amplitude of hlm mode **/
  REAL8Sequence *freqs_In,                /**< Frequency array to evaluate model or fmin, fmax  **/
  IMRPhenomXWaveformStruct *pWF,          /**< Structure of the 22 mode **/
  IMRPhenomXHMWaveformStruct *pWFHM,      /**< waveform parameters lm mode */
  IMRPhenomXHMAmpCoefficients *pAmp,      /**< Amplitude coefficients lm */
  UNUSED IMRPhenomXHMPhaseCoefficients *pPhase  /**< Phase coefficients 22 */
)
{
  #if DEBUG == 1
  printf("\n **** IMRPhenomXHM_Amplitude **** \n");
  #endif

  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  int status = 0;
  REAL8Sequence *freqs;
  UINT4 offset = SetupWFArraysReal(&freqs, amplm, freqs_In, pWF, ligotimegps_zero);

  #if DEBUG == 1
  printf("\n***Length@freqs, offset %i %i",freqs_In->length,offset);
  printf("\n\nfstart, fend = %.16f %.16f\n\n", freqs_In->data[0], freqs_In->data[freqs_In->length-1]);
  printf("\n***Length@freqs, offset %i %i",freqs->length,offset);
  printf("\n\nfstart, fend = %.16f %.16f\n\n", freqs->data[0], freqs->data[freqs->length-1]);
  #endif

  /* Loop over frequencies to generate waveform */
  if(pWFHM->Ampzero==0){
    IMRPhenomX_UsefulPowers powers_of_Mf;
    UINT4 initial_status = XLAL_SUCCESS;
    REAL8 Msec = pWF->M_sec;
    REAL8 amp;

    for (UINT4 idx = 0; idx < freqs->length; idx++)
    {
      REAL8 Mf = Msec * freqs->data[idx];
      initial_status = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
      if(initial_status != XLAL_SUCCESS)
      {
        status = initial_status;
        XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
      }
      else
      {
        amp = IMRPhenomXHM_RD_Amp_Ansatz(powers_of_Mf.itself, pWFHM, pAmp);
        /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
        ((*amplm)->data->data)[idx+offset] =  amp;
      }
    }
  }
  /* Free allocated memory */
  XLALDestroyREAL8Sequence(freqs);

  return status;

}

/* Spherical phase evaluated in an input frequency array */
static int IMRPhenomXHM_Phase(
  REAL8FrequencySeries **phaselm,         /**<[out] phase of hlm mode **/
  REAL8Sequence *freqs_In,                /**< Frequency array to evaluate model or fmin, fmax  **/
  IMRPhenomXWaveformStruct *pWF,          /**< Structure of the 22 mode **/
  IMRPhenomXAmpCoefficients *pAmp22,      /**< Amplitude coefficients 22 */
  IMRPhenomXPhaseCoefficients *pPhase22,  /**< Phase coefficients 22 */
  IMRPhenomXHMWaveformStruct *pWFHM,      /**< waveform parameters lm mode */
  IMRPhenomXHMAmpCoefficients *pAmp,      /**< Amplitude coefficients lm */
  IMRPhenomXHMPhaseCoefficients *pPhase  /**< Phase coefficients 22 */
)
{
  #if DEBUG == 1
  printf("\n **** IMRPhenomXHM_Phase **** \n");
  #endif

  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  int status = 0;
  REAL8Sequence *freqs;
  UINT4 offset = SetupWFArraysReal(&freqs, phaselm, freqs_In, pWF, ligotimegps_zero);

  #if DEBUG == 1
  printf("\n***Length@freqs, offset %i %i",freqs->length,offset);
  printf("\n\nfstart, fend = %.16f %.16f\n\n", freqs->data[0], freqs->data[freqs->length-1]);
  #endif

  if(pWFHM->Ampzero==0){
    IMRPhenomX_UsefulPowers powers_of_Mf;
    UINT4 initial_status = XLAL_SUCCESS;
    REAL8 Msec = pWF->M_sec;
    REAL8 phi;

    /* Loop over frequencies to generate waveform */
    if(pWFHM->MixingOn==1){
      //double twopi_m1=1./(2.*LAL_PI);
      UINT4 enter = 1;
      REAL8 Mf = Msec * freqs->data[0];
      initial_status = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
      if(initial_status != XLAL_SUCCESS)
      {
        status = initial_status;
        XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
      }
      else
      {
        //REAL8 phiprevious = 0;
        for (UINT4 idx = 0; idx < freqs->length; idx++)
        {
          Mf = Msec * freqs->data[idx];
          initial_status = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
          if(initial_status != XLAL_SUCCESS)
          {
            status = initial_status;
            XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
          }
          else
          {
            phi = IMRPhenomXHM_Phase_ModeMixing(Mf, &powers_of_Mf, pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF);
            /* Only the first coarse point in the RD needs the unwrapping.
            If we remove the enter condition we would get a nice and smooth coarse phase up to pAmp->fAmpMatchIM. */
            if(Mf > pPhase->fPhaseMatchIM && idx>0 && enter == 1){
              #if DEBUG == 1
                  printf("\nExtrapolating intermediate phase out of its range for one point\n");
              #endif
              phi = IMRPhenomXHM_Inter_Phase_AnsatzInt(Mf, &powers_of_Mf, pWFHM, pPhase);
              phi = phi + pPhase->deltaphiLM;
              //phi = phi+2.*LAL_PI*round((phiprevious-phi)*twopi_m1);
              enter = 0;
            }
            //phiprevious = phi;
            /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
            ((*phaselm)->data->data)[idx+offset] = phi;
          }
        }
      }
    }
    else{
      /* Loop over frequencies to generate waveform */
      for (UINT4 idx = 0; idx < freqs->length; idx++)
      {
        REAL8 Mf    = Msec * freqs->data[idx];
        initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
        if(initial_status != XLAL_SUCCESS)
        {
          status = initial_status;
          XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
        }
        else
        {
          phi = IMRPhenomXHM_Phase_noModeMixing(Mf, &powers_of_Mf, pPhase, pWFHM, pWF);
          /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
          ((*phaselm)->data->data)[idx+offset] = phi;
        }
      }
    }
  }
  /* Free allocated memory */
  XLALDestroyREAL8Sequence(freqs);

  return status;
}


/* Ringdown phase ansatz evaluated in an input frequency array */
static int IMRPhenomXHM_PhaseMixing(
  REAL8FrequencySeries **phaselm,         /**<[out] phase of hlm mode **/
  REAL8Sequence *freqs_In,                /**< Frequency array to evaluate model or fmin, fmax  **/
  IMRPhenomXWaveformStruct *pWF,          /**< Structure of the 22 mode **/
  IMRPhenomXHMWaveformStruct *pWFHM,      /**< waveform parameters lm mode */
  IMRPhenomXHMPhaseCoefficients *pPhase  /**< Phase coefficients 22 */
)
{
  #if DEBUG == 1
  printf("\n **** IMRPhenomXHM_PhaseMixing **** \n");
  #endif

  /* Set LIGOTimeGPS */
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}

  int status = 0;
  REAL8Sequence *freqs;
  UINT4 offset = SetupWFArraysReal(&freqs, phaselm, freqs_In, pWF, ligotimegps_zero);

  #if DEBUG == 1
  printf("\n\nfstart, fend = %.16f %.16f\n\n", freqs_In->data[0], freqs_In->data[freqs_In->length-1]);
  printf("\n***Length@freqs, offset %i %i",freqs->length,offset);
  printf("\n\nfstart, fend = %.16f %.16f\n\n", freqs->data[0], freqs->data[freqs->length-1]);
  #endif

  if(pWFHM->Ampzero==0){
    IMRPhenomX_UsefulPowers powers_of_Mf;
    UINT4 initial_status = XLAL_SUCCESS;
    REAL8 Msec = pWF->M_sec;
    REAL8 phi;

    /* Loop over frequencies to generate waveform */
    for (UINT4 idx = 0; idx < freqs->length; idx++)
    {
      REAL8 Mf    = Msec * freqs->data[idx];
      initial_status     = IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);
      if(initial_status != XLAL_SUCCESS)
      {
        status = initial_status;
        XLALPrintError("IMRPhenomX_Initialize_Powers failed for Mf, initial_status=%d",initial_status);
      }
      else
      {
        phi = IMRPhenomXHM_RD_Phase_AnsatzInt(Mf, &powers_of_Mf, pWFHM, pPhase);
        /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
        ((*phaselm)->data->data)[idx+offset] = phi;
      }
    }
  }
  /* Free allocated memory */
  XLALDestroyREAL8Sequence(freqs);

  return status;
}

/* Log function */
static double logbase(double base, double x) {
  return (log(x) / log(base));
}

/* Right hand of eq. 2.27 in arXiv:2001.10897. */
static double deltaF_mergerBin(REAL8 fdamp, REAL8 alpha4, REAL8 abserror)
{
      double aux = sqrt(sqrt(3.)*3);
      return 4. * fdamp * sqrt(abserror/fabs(alpha4)) / aux;
}

/* Correspond to eqs. 2.28 and 2.31 in arXiv:2001.10897 */
static double deltaF_ringdownBin(REAL8 fdamp, REAL8 alpha4, REAL8 LAMBDA, REAL8 abserror){
  double dfphase = 5*fdamp*sqrt(abserror*0.5/fabs(alpha4));
  double dfamp   = sqrt(2*abserror)/fabs(LAMBDA);
  if (dfphase <= dfamp){
    return dfphase;
  }
  else{
    return dfamp;
  }
}

// static double deltaF_ringdownBinAmp(REAL8 fdamp, REAL8 lambda, REAL8 sigma, REAL8 relerror){
//    return sqrt(sqrt(24.*relerror))*sigma*fdamp/lambda;
// }
