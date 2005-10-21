/*********************************** <lalVerbatim file="StackSlideIsolatedCV">
Authors: Mendell, Gregory
$Id$
************************************ </lalVerbatim> */

/* REVISIONS: */
/* 04/12/05 gam; Move from using StackSlideOld to using StackSlide function. */
/* 04/12/05 gam; add RunStackSlideIsolatedMonteCarloSimulation to StackSlideIsolated.c. Previous REVISIONS: */
/* 05/11/04 gam; Add RunStackSlideMonteCarloSimulation to software inject signals into the SFTs for Monte Carlo simulations */
/* 05/21/04 gam; Normalize the SFTs from LALSignalToSFTs as per the normalization in makefakedata_v2.c and lalapps/src/pulsar/make_sfts.c. */
/* 05/21/04 gam; Save the input SFTs and Parameter space; run Monte Carlo simulation on the parameter space with random offsets */
/* 05/26/04 gam; Continue work on Monto Carlo code. */
/* 05/26/04 gam; Change finishedSUMs to finishedSUMs; add startSUMs; defaults are TRUE; use to control I/O during Monte Carlo */
/* 05/26/04 gam; Add whichMCSUM = which Monte Carlo SUM; default is -1. */
/* 05/28/04 gam; Use LALUniformDeviate from LAL utilities package (include <lal/Random.h>) to generate random mismatch during Monte Carlo. */
/* 06/01/04 gam; Make sure CHECKSTATUSPTR called after every LAL function call */
/* 06/01/04 gam; pPulsarSignalParams->startTimeGPS and duration should be based on gpsStartTime and duration, not actualStartTime */
/* 06/01/04 gam; For now fix DeltaRA, DeltaDec,and DeltaFDeriv1 to default values when finding mismatch during Monte Carlo */
/* 07/14/04 gam; add option to use function LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs; see LAL inject packge GeneratePulsarSignal.c and .h */
/* 07/19/04 gam; if (params->testFlag & 4) > 0 use LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs for Monte Carlo simulations */
/* 08/02/04 gam; if (params->testFlag & 4) > 0 ComputeSky uses reference time GPSin: params->timeStamps[0].gpsSeconds, params->timeStamps[0].gpsNanoSeconds */
/* 08/02/04 gam; Set pSFTandSignalParams->resTrig = 0 to avoid serious bug in LALFastGeneratePulsarSFTs when using lookup tables (LUTs) for trig calls. */
/* 12/06/04 gam; get params->sampleRate, = effective sample rate, from the SFTs; calculate params->deltaT after reading SFTs. */
/* 12/06/04 gam; add params->gpsEpochStartTimeNan; get gpsEpochStartTime, gpsEpochStartTimeNan, and gpsStartTime from command line; */ 
/* 12/06/04 gam; if (params->testFlag & 8) > 0 use fixed values for psi and cosIota during Monte Carlo simulations */
/* 05/06/05 gam; If params->debugOptionFlag & 128 > 0 and isolated case, just creates a SUM from the STKs without sliding */
/* 05/19/05 gam; Add INT4 *sumBinMask; params->sumBinMask == 0 if bin should be excluded from search or Monte Carlo due to cleaning */
/* 05/24/05 gam; make maxPower and totalEventCount part of params; change finishSUMs to finishPeriodicTable; end xml in FinalizeSearch */
/* 05/24/05 gam; if (params->testFlag & 16 > 0) use results from prior jobs in the pipeline and report on current MC results */
/* 05/24/05 gam; if (params->testFlag & 32 > 0) iterate MC up to 10 times to converge on desired confidence */
/* 05/24/05 gam; if (params->debugOptionFlag & 32 > 0) print Monte Carlo Simulation results to stdout */
/* 05/24/05 gam; add StackSlideMonteCarloResultsTable */
/* 07/15/2005 gam; need to save params->sumBinMask[iFreq] as savSumBinMask[iFreq] and inject where savSumBinMask[iFreq] !=0, and reset params->sumBinMask[0] == 1 */
/* 07/15/2005 gam; Change RunStackSlideIsolatedMonteCarloSimulation to inject at a */
/*                 random point in the parameters space a random number of times.  */
/* 07/29/05 gam; if (params->testFlag & 64) > 0 set searchSurroundingPts == 1 and  */
/*               search surrounding parameters space pts; else search nearest only */
/* 08/31/05 gam; In StackSlideComputeSky set ssbT0 to gpsStartTime, which is gpsEpochStartTime in this code; this now gives the epoch that defines T_0 at the SSB! */
/*               This affects GPSin and tRef as well. */
/* 09/01/05 gam; If params->numMCRescalings > 0 use this and params->rescaleMCFractionSFTs to */
/*               rescale SFTs to run numMCRescalings Monte Carlo simulations in parallel.     */
/* 09/05/05 gam; To inject isotropically on the sky, inject uniform distribution for sin(DEC). */
/* 09/06/05 gam; Change params->maxMCfracErr to params->maxMCErr, the absolute error in confidence for convergence. */
/* 09/08/05 gam; Add Dterms to params used with LALFastGeneratePulsarSFTs. */
/* 09/12/05 gam; if (params->testFlag & 128) > 0 make BLKs and STKs narrower band based on extra bins */
/* 09/16/06 gam; In CountOrAssignSkyPosData and MC code, for each DEC adjust deltaRA and numRA to evenly  */
/*               space grid points so that deltaRA used is <= input deltaRA from the command line divided */
/*               by cos(DEC). This fixes a problem where the last grid point could wrap around the sphere */
/*               and be closer to the first grid point that the spacing between other grid points.        */
/* 10/20/06 gam; if ( (params->weightFlag & 8) > 0 ) renorm the BLK (SFT) data up front with the inverse mean absolute value of the BLK data = params->blkRescaleFactor */

/*********************************************/
/*                                           */
/* START SECTION: define preprocessor flags  */
/*                                           */
/*********************************************/
#define INCLUDE_PRINT_PEAKS_TABLE_CODE
#define INCLUDE_PRINT_REPORTMCRESULTS_CODE
/* #define DEBUG_STACKSLIDE_ISOLATED */
/* #define DEBUG_SUM_TEMPLATEPARAMS */
/* #define PRINT_STACKSLIDE_BINOFFSETS */
/* #define PRINT_STACKSLIDE_BINMISMATCH */
/* #define PRINT_STACKSLIDE_SPECIAL_DEBUGGING_INFO */
/* #define DEBUG_MONTECARLOTIMEDOMAIN_DATA */
/* #define PRINT_MONTECARLOTIMEDOMAIN_DATA */
/* #define DEBUG_MONTECARLOSFT_DATA */
/* #define PRINT_MONTECARLOSFT_DATA */
/* #define PRINTCOMPARISON_INPUTVSMONTECARLOSFT_DATA */
/* #define DEBUG_RANDOMTRIALPARAMETERS */
/* #define DEBUG_LALFASTGENERATEPULSARSFTS */
/* #define DEBUG_SETFIXED_RANDVAL */
/* #define PRINT_ONEMONTECARLO_OUTPUTSFT */
/* #define PRINT_MAXPOWERANDBINEACHSFT */
#define PRINT_SLIDINGTURNEDOFF_WARNING
/*********************************************/
/*                                           */
/* END SECTION: define preprocessor flags    */
/*                                           */
/*********************************************/
/*********************************************/
/*                                           */
/* START SECTION: include header files       */
/* (Note that most include statements        */
/*  are in the header files below.)          */
/*                                           */
/*********************************************/
#include "StackSlideIsolated.h"
/*********************************************/
/*                                           */
/* END SECTION: include header files         */
/*                                           */
/*********************************************/

NRCSID( STACKSLIDEISOLATEDC,  "$Id$");

/*********************************************************************************/
/*              START FUNCTION: StackSlideIsolated                               */
/*********************************************************************************/
void StackSlideIsolated (
    LALStatus                        *status,
    SnglStackSlidePeriodicTable      *loudestPeaksArray,
    LALFindStackSlidePeakOutputs     *pLALFindStackSlidePeakOutputs,
    LALFindStackSlidePeakParams      *pLALFindStackSlidePeakParams,
    LALUpdateLoudestStackSlideParams *pLALUpdateLoudestStackSlideParams,
    LALDetector                      *cachedDetector,
    StackSlideParams                 *stksldParams,
    StackSlideSearchParams           *params
)
{
 INT4 iSky, iFreqDeriv, k;
 INT4 kSUM = -1;
  
 SnglStackSlidePeriodicTable *peakSav;
  
 /* 04/12/05 gam; add parameters for StackSlideComputeSky and StackSlide */
 StackSlideSkyParams *csParams;
 TdotsAndDeltaTs *pTdotsAndDeltaTs;
 EarthState earth;
 EmissionTime emit;
  
 #ifdef DEBUG_STACKSLIDE_ISOLATED
   fprintf(stdout, "\nStart StackSlideIsolated\n\n");
   fflush(stdout);
 #endif    

 INITSTATUS( status, "StackSlideIsolated", STACKSLIDEISOLATEDC );
 ATTATCHSTATUSPTR(status);
  
 /* 04/12/05 gam; allocate memory and initialize parameters for StackSlideComputeSky and StackSlide */
 csParams=(StackSlideSkyParams *)LALMalloc(sizeof(StackSlideSkyParams));
 csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
 /* 12/05/05 gam; IMPORTANT, must set these to the epoch that gives T0 at SSB; this differs when running Monte Carlo */
 if ( (params->testFlag & 4) > 0 ) {
      stksldParams->gpsStartTimeSec = (UINT4)params->timeStamps[0].gpsSeconds;
      stksldParams->gpsStartTimeNan = (UINT4)params->timeStamps[0].gpsNanoSeconds;
      csParams->gpsStartTimeSec     = (UINT4)params->timeStamps[0].gpsSeconds;
      csParams->gpsStartTimeNan     = (UINT4)params->timeStamps[0].gpsNanoSeconds;
 } else {
      stksldParams->gpsStartTimeSec = params->gpsEpochStartTimeSec;
      stksldParams->gpsStartTimeNan = params->gpsEpochStartTimeNan;
      csParams->gpsStartTimeSec     = params->gpsEpochStartTimeSec;
      csParams->gpsStartTimeNan     = params->gpsEpochStartTimeNan;
 }
 csParams->tGPS=params->timeStamps;
 csParams->spinDwnOrder=params->numSpinDown;
 csParams->mObsSFT=params->numSTKs;
 csParams->tSFT=params->tSTK;
 csParams->edat=params->edat;
 csParams->emit=&emit;
 csParams->earth=&earth;
 csParams->baryinput = &params->baryinput;
 pTdotsAndDeltaTs=(TdotsAndDeltaTs *)LALMalloc(sizeof(TdotsAndDeltaTs));
 pTdotsAndDeltaTs->vecTDots = (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSTKs);
 if (params->numSpinDown>0) {
    pTdotsAndDeltaTs->vecDeltaTs = (REAL8 **)LALMalloc(sizeof(REAL8 *)*params->numSTKs);          
    for(k=0;k<params->numSTKs;k++) {
        pTdotsAndDeltaTs->vecDeltaTs[k]  = (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSpinDown);          
    }
 }

 /* Loop over sky positions */
 for(iSky=0;iSky<params->numSkyPosTotal;iSky++) {

   /* 10/28/04 gam; weight STKs; the weights depend on F_+ and F_x when params->weightSTKsIncludingBeamPattern is True */
   if (params->weightSTKsIncludingBeamPattern) {       
        /* 10/28/04 gam; get squared detector response F_+^2 or F_x^2 for one sky position, one polarization angle, for midpoints of a timeStamps */
        GetDetResponseTStampMidPts(status->statusPtr, params->detResponseTStampMidPts, params->timeStamps, params->numSTKs, params->tSTK,
             cachedDetector, params->skyPosData[iSky], params->orientationAngle, COORDINATESYSTEM_EQUATORIAL, params->plusOrCross);
         
        /* 10/28/04 gam; apply powerFlux style weights including detector beam pattern response */
        WeightSTKsIncludingBeamPattern(params->STKData,
              params->savSTKData,
              params->inverseSquareMedians,
              params->sumInverseSquareMedians,
              params->detResponseTStampMidPts,
              params->numSTKs, params->nBinsPerSTK, params->tSTK);
   }
   
   /* setup to call StackSlideComputeSky */

   csParams->skyPos[0]=params->skyPosData[iSky][0];
   csParams->skyPos[1]=params->skyPosData[iSky][1];
   csParams->baryinput->alpha=params->skyPosData[iSky][0];
   csParams->baryinput->delta=params->skyPosData[iSky][1];

   /* get Tdot and DeltaT's for this sky position */
   StackSlideComputeSky(status->statusPtr, pTdotsAndDeltaTs, csParams);
   CHECKSTATUSPTR (status);
      
   /* Loop over frequency derivatives */
   for(iFreqDeriv=0;iFreqDeriv<params->numFreqDerivIncludingNoSpinDown;iFreqDeriv++) {
    
    if (params->numSpinDown>0) {
       stksldParams->freqDerivData = params->freqDerivData[iFreqDeriv];
    }
    
    kSUM++; /* 04/12/05 gam; increment sum number */
    
    /* StackSlideOld(status->statusPtr,params->SUMData,params->STKData,stksldParams);*/ /*Feb 05/14 vir*/
    /* 04/12/05 gam; Move from using StackSlideOld to using StackSlide function. */
    /* 05/06/05 gam; If params->debugOptionFlag & 128 > 0 just creates a SUM from the STKs without sliding */
    if ((params->debugOptionFlag & 128) > 0 ) {
      #ifdef PRINT_SLIDINGTURNEDOFF_WARNING
        fprintf(stdout,"\n\n !!! WARNING: SLIDING TURNED OFF !!! \n\n");
        fflush(stdout);
      #endif
      SumStacks(status->statusPtr,params->SUMData,params->STKData,stksldParams);
    } else {
      StackSlide(status->statusPtr,params->SUMData,params->STKData,pTdotsAndDeltaTs,stksldParams);
    }
    CHECKSTATUSPTR (status);
    
    /* 02/25/05 gam; utility for printing one SUM to a file */
    if (params->outputSUMFlag > 0) {
      printOneStackSlideSUM(params->SUMData[0],params->outputSUMFlag,params->outputFile,kSUM,params->f0SUM,params->dfSUM,params->nBinsPerSUM,params->numSUMsTotal);
    }

    #ifdef DEBUG_SUM_TEMPLATEPARAMS
        fprintf(stdout,"\nTemplate parameters for SUM #%i:\n",kSUM);
        fprintf(stdout,"RA = %18.10f\n",params->skyPosData[iSky][0]);
        fprintf(stdout,"DEC = %18.10f\n",params->skyPosData[iSky][1]);
        for(k=0;k<params->numSpinDown;k++)
        {
           fprintf(stdout,"FREQDERIV%i = %18.10f\n",k+1,params->freqDerivData[iFreqDeriv][k]);
        } /* END for(k=0;k<params->numSpinDown;k++) */          
    #endif
    
    /* 01/21/04 gam; Call LALFindStackSlidePeaks only if thresholdFlag > 0 */
    if (params->thresholdFlag > 0) {
         /* 02/04/04 gam; Initialize list */ /* 02/11/04 gam; initial params->peaks here: */
         params->peaks = (SnglStackSlidePeriodicTable *)LALMalloc(sizeof(SnglStackSlidePeriodicTable));      
         params->peaks->next = NULL;
         peakSav = params->peaks; /* Save pointer to first peak */
         if (params->whichMCSUM > 0) {
           pLALFindStackSlidePeakParams->iSUM = params->whichMCSUM; /* 05/26/04 gam; this is the number of the Monte Carlo SUM. */
         } else {
           pLALFindStackSlidePeakParams->iSUM = kSUM; /* 02/04/04 gam added this and next 2 lines */
         }
         pLALFindStackSlidePeakParams->iSky = iSky;
         pLALFindStackSlidePeakParams->iDeriv = iFreqDeriv;
         pLALFindStackSlidePeakOutputs->pntrToPeaksPntr = &(params->peaks);  /* 03/02/04 gam; added this and changed next line */
         LALFindStackSlidePeaks(pLALFindStackSlidePeakOutputs, params->SUMData[0], pLALFindStackSlidePeakParams);
                 
         params->peaks = peakSav->next;  /* Go back to the beginning of the list */
         LALFree(peakSav);               /* No longer need this place holder */

         while (params->peaks) {
           peakSav = params->peaks;
           params->totalEventCount++; /* 05/25/05 gam */
           if (pLALFindStackSlidePeakParams->updateMeanStdDev) {
              /* 03/02/04 gam; update pw_mean_thissum, pw_stddev_thissum, pwr_snr using pwMeanWithoutPeaks and pwStdDevWithoutPeaks */              
              if ( (*(pLALFindStackSlidePeakOutputs->binsWithoutPeaks)) > 0) {
                /* Note that LALFindStackSlidePeaks ensures that pwStdDevWithoutPeaks is never 0. */
                params->peaks->pw_mean_thissum = *(pLALFindStackSlidePeakOutputs->pwMeanWithoutPeaks);
                params->peaks->pw_stddev_thissum = *(pLALFindStackSlidePeakOutputs->pwStdDevWithoutPeaks);
                params->peaks->pwr_snr = params->peaks->power / (*(pLALFindStackSlidePeakOutputs->pwStdDevWithoutPeaks));
              }
           }
           if (params->peaks->power > params->maxPower) {
              params->maxPower = params->peaks->power; /* 05/25/05 gam */ /* 02/09/04 gam; keep track for search summary */
           }
           
           /* 02/17/04 gam; keep an array of loudest peaks */
           if (params->outputLoudestFromPeaks) {
               LALUpdateLoudestFromPeaks(loudestPeaksArray, params->peaks, params->nBinsPerOutputEvent);
           } else {
              #ifdef INCLUDE_PRINT_PEAKS_TABLE_CODE
                if ((params->debugOptionFlag & 2) > 0 ) {
                  fprintf(stdout,"%11i %11i %10.6f %10.6f %14.6e %11.6f %10.3f %8.3f %12i \n",params->totalEventCount,params->peaks->sum_no,params->peaks->sky_ra,params->peaks->sky_dec,params->peaks->fderiv_1,params->peaks->frequency,params->peaks->power,params->peaks->pwr_snr,params->peaks->width_bins);
                  fflush(stdout);
                }
              #endif
              if (params->outputEventFlag > 0) {
                 /*  end previous row with comma */
                 if ( params->xmlStream->first ) {
                    params->xmlStream->first = 0;
                 } else { 
                    fprintf( params->xmlStream->fp, ",\n" );
                 }
                 /* print out the row */
                 /* 02/09/04 gam; remove fderiv_2-5; remove pw_max_thissum and freq_max_thissum; change false_alarm_prob_upperlimit to false_alarm_prob */
                 /* 02/12/04 gam; Add freq_index to SnglStackSlidePeriodicTable */
                 fprintf( params->xmlStream->fp, SNGL_LOCAL_STACKSLIDEPERIODIC_ROW,
		   params->peaks->sum_no,
		   params->peaks->event_no_thissum,
		   params->peaks->overlap_event,
		   params->peaks->frequency,
		   params->peaks->freq_index,
		   params->peaks->power,
		   params->peaks->width_bins,
		   params->peaks->num_subpeaks,
		   params->peaks->pwr_snr,
		   params->peaks->false_alarm_prob,
		   params->peaks->goodness_of_fit,
		   params->peaks->sky_ra,
		   params->peaks->sky_dec,
		   params->peaks->fderiv_1,
		   params->peaks->pw_mean_thissum,
		   params->peaks->pw_stddev_thissum
                );
              } /* end if (params->outputEventFlag > 0 */
	   } /* end if (params->outputLoudestFromPeaks) */
           params->peaks = params->peaks->next;
           LALFree(peakSav);
         } /* end while (params->peaks) */
    /* 02/11/04 gam; next is end if (params->thresholdFlag > 0) */       
    } else if (params->outputLoudestFromSUMs) {
         /* 02/17/04 gam */
         if (params->whichMCSUM > 0) {
           pLALUpdateLoudestStackSlideParams->iSUM = params->whichMCSUM; /* 05/26/04 gam; this is the number of the Monte Carlo SUM. */
         } else {
           pLALUpdateLoudestStackSlideParams->iSUM = kSUM;
         }
         pLALUpdateLoudestStackSlideParams->iSky = iSky;
         pLALUpdateLoudestStackSlideParams->iDeriv = iFreqDeriv;
         LALUpdateLoudestFromSUMs(loudestPeaksArray, params->SUMData[0], pLALUpdateLoudestStackSlideParams);
    } /* end else if (params->outputLoudestFromSUMs) */
  } /* 04/12/05 gam; END for(iFreqDeriv=0;iFreqDeriv<params->numFreqDerivIncludingNoSpinDown;iFreqDeriv++) */
 } /* 04/12/05 gam; END for(iSky=0;iSky<params->numSkyPosTotal;iSky++) */

 /* 02/17/04 gam; output the loudest events */
 if ( (params->outputLoudestFromPeaks || params->outputLoudestFromSUMs) && (params->outputEventFlag > 0) ) {
    for (k=0;k<params->keepThisNumber; k++) {
       #ifdef INCLUDE_PRINT_PEAKS_TABLE_CODE
         /* 04/15/04 gam */
         if ((params->debugOptionFlag & 2) > 0 ) {
          fprintf(stdout,"%11i %11i %10.6f %10.6f %14.6e %11.6f %10.3f %8.3f %12i \n",k+1,loudestPeaksArray[k].sum_no,loudestPeaksArray[k].sky_ra,loudestPeaksArray[k].sky_dec,loudestPeaksArray[k].fderiv_1,loudestPeaksArray[k].frequency,loudestPeaksArray[k].power,loudestPeaksArray[k].pwr_snr,loudestPeaksArray[k].width_bins);
          fflush(stdout);
         }
       #endif    
       if (loudestPeaksArray[k].power > params->maxPower) {
              params->maxPower = loudestPeaksArray[k].power; /* 05/25/05 gam */ /* 02/20/04 gam; keep track for search summary */
       }       
       /*  end previous row with comma */
       if ( params->xmlStream->first ) {
          params->xmlStream->first = 0;
       } else { 
          fprintf( params->xmlStream->fp, ",\n" );
       }
       fprintf( params->xmlStream->fp, SNGL_LOCAL_STACKSLIDEPERIODIC_ROW,
        loudestPeaksArray[k].sum_no,
        loudestPeaksArray[k].event_no_thissum,
        loudestPeaksArray[k].overlap_event,
        loudestPeaksArray[k].frequency,
        loudestPeaksArray[k].freq_index,
        loudestPeaksArray[k].power,
        loudestPeaksArray[k].width_bins,
        loudestPeaksArray[k].num_subpeaks,
        loudestPeaksArray[k].pwr_snr,
        loudestPeaksArray[k].false_alarm_prob,
        loudestPeaksArray[k].goodness_of_fit,
        loudestPeaksArray[k].sky_ra,
        loudestPeaksArray[k].sky_dec,
        loudestPeaksArray[k].fderiv_1,
        loudestPeaksArray[k].pw_mean_thissum,
        loudestPeaksArray[k].pw_stddev_thissum
       );
    } /* end for (k = 0; k < arraySize; k++) */
 } /* end if ( (params->outputLoudestFromPeaks || params->outputLoudestFromSUMs) && (params->outputEventFlag > 0) ) */

  /* 04/12/05 gam; allocate memory for parameters for StackSlideComputeSky and StackSlide */
 LALFree(csParams->skyPos);
 LALFree(csParams);
 if (params->numSpinDown>0) {
     for(k=0;k<params->numSTKs;k++) {
        LALFree(pTdotsAndDeltaTs->vecDeltaTs[k]);
     }
     LALFree(pTdotsAndDeltaTs->vecDeltaTs);
 }
 LALFree(pTdotsAndDeltaTs->vecTDots);
 LALFree(pTdotsAndDeltaTs);
  
 CHECKSTATUSPTR (status);
 DETATCHSTATUSPTR (status);

} /* END StackSlideIsolated */
/*********************************************************************************/
/*              END FUNCTION: StackSlideIsolated                                 */
/*********************************************************************************/

/*************************************************************/
/*                                                           */
/* START FUNCTION: RunStackSlideIsolatedMonteCarloSimulation */
/* Injects fake signals and runs Monte Carlo simulation.     */
/*                                                           */
/*************************************************************/
void RunStackSlideIsolatedMonteCarloSimulation(LALStatus *status, StackSlideSearchParams *params, INT4 nSamples)
{
  INT4    i = 0; /* all purpose index */
  INT4    j = 0; /* all purpose index */
  INT4    k = 0; /* all purpose index */
  PulsarSignalParams *pPulsarSignalParams = NULL;
  REAL4TimeSeries *signal = NULL;
  LIGOTimeGPS GPSin;            /* reference-time for pulsar parameters at the detector; will convert to SSB! */
  LALDetector cachedDetector;
  REAL8 cosIota;                    /* cosine of inclination angle iota of the source */
  REAL8 h_0 = params->threshold4;   /* Source amplitude */
  SFTParams *pSFTParams = NULL;
  LIGOTimeGPSVector timestamps;
  SFTVector *outputSFTs = NULL;
  /* 07/14/04 gam; next two are needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  SkyConstAndZeroPsiAMResponse *pSkyConstAndZeroPsiAMResponse;
  SFTandSignalParams *pSFTandSignalParams;
  REAL4 renorm;               /* Need to renormalize SFTs to account for different sample rates */  
  FFT **savBLKData = NULL;    /* 05/21/04 gam; save the input noise SFTs */
  REAL8 **savSkyPosData;      /* 05/26/04 gam; save Sky Position Data */
  REAL8 **savFreqDerivData;   /* 05/26/04 gam; save Frequency Derivative Data */
  INT4 numSkyPosTotal = 0;    /* 05/26/04 gam; total number of Sky positions to cover */
  INT4 numFreqDerivTotal = 0; /* 05/26/04 gam; total number of Frequency evolution models to cover */
  INT4 numParamSpacePts = 0;  /* 05/26/04 gam; total number of points in the parameter space to cover */
  INT4 numSUMsTotal = 0;      /* 05/26/04 gam; total Number of Sums output = numSUMsPerParamPt*numParamSpacePts */
  INT4 numSUMsTotalm1 = -1;                 /* 05/26/04 gam */
  INT4 numFreqDerivIncludingNoSpinDown = 0; /* 05/26/04 gam */
  INT4 kSUM = 0;                            /* 05/26/04 gam */
  INT4 iFreq = 0;                           /* 05/26/04 gam */  
  REAL8 f0SUM = 0.0;                        /* 05/26/04 gam */
  REAL8 bandSUM = 0.0;                      /* 05/26/04 gam */
  INT4 nBinsPerSUM = 0;                     /* 05/26/04 gam */
  INT4 nBinsPerSUMm1 = -1;                  /* 05/26/04 gam */  
  REAL8 real8NBinsPerSUM;                   /* 07/15/05 gam */
  INT4 seed=0; /* 05/28/04 gam; next 14 are for computing random mismatch */
  REAL4 randval;
  REAL8 real8RandVal;
  RandomParams *randPar=NULL;
  FILE *fpRandom;
  INT4 rndCount;
  REAL8 cosTmpDEC;
  REAL8 tmpDeltaRA;
  INT4 tmpNumRA; /* 09/16/05 gam */
  REAL8 tmpRA;
  REAL8 tmpDec;
  INT4 iRA;
  INT4 iDec;
  INT4 iDeriv;  
  REAL8 DeltaRA  = params->stksldSkyPatchData->deltaRA;
  REAL8 DeltaDec = params->stksldSkyPatchData->deltaDec;
  REAL8 DeltaFDeriv1 = params->deltaFDeriv1;
  REAL8 DeltaFDeriv2 = params->deltaFDeriv2;
  REAL8 DeltaFDeriv3 = params->deltaFDeriv3;
  REAL8 DeltaFDeriv4 = params->deltaFDeriv4;
  REAL8 DeltaFDeriv5 = params->deltaFDeriv5;
  REAL8 startRA  = params->stksldSkyPatchData->startRA;
  REAL8 stopRA   = params->stksldSkyPatchData->stopRA;  /* 09/16/05 gam */
  REAL8 startDec = params->stksldSkyPatchData->startDec;
  REAL8 startFDeriv1 = params->startFDeriv1;
  REAL8 startFDeriv2 = params->startFDeriv2;
  REAL8 startFDeriv3 = params->startFDeriv3;
  REAL8 startFDeriv4 = params->startFDeriv4;
  REAL8 startFDeriv5 = params->startFDeriv5;
  REAL8 rangeRA  = params->stksldSkyPatchData->stopRA - params->stksldSkyPatchData->startRA;
  /* REAL8 rangeDec = params->stksldSkyPatchData->stopDec - params->stksldSkyPatchData->startDec; */
  REAL8 rangeDec = sin(params->stksldSkyPatchData->stopDec) - sin(params->stksldSkyPatchData->startDec); /* 09/05/05 gam */
  REAL8 sinStartDec = sin(startDec);                                                                     /* 09/05/05 gam */
  REAL8 sinStartDecmhalfDeltaDec = sin(startDec-(0.5*DeltaDec));                                         /* 09/05/05 gam */
  REAL8 rangeDecStartDecpmhDeltaDec = sin(startDec+(0.5*DeltaDec)) - sin(startDec-(0.5*DeltaDec));       /* 09/05/05 gam */
  REAL8 rangeFDeriv1 = params->stopFDeriv1 - params->startFDeriv1;
  REAL8 rangeFDeriv2 = params->stopFDeriv2 - params->startFDeriv2;
  REAL8 rangeFDeriv3 = params->stopFDeriv3 - params->startFDeriv3;
  REAL8 rangeFDeriv4 = params->stopFDeriv4 - params->startFDeriv4;
  REAL8 rangeFDeriv5 = params->stopFDeriv5 - params->startFDeriv5;
  INT4 firstNonExcludedBin = -1; /* 05/19/05 gam */
  INT4 lastNonExcludedBin  = -1; /* 05/19/05 gam */
  /* 05/24/05 gam; use results from prior jobs in the pipeline and report on current MC results */
  BOOLEAN reportMCResults = 0;
  REAL4 priorLoudestEvent = 0.0;
  REAL8 priorStartFreq = 0.0;
  REAL8 priorBand = 0.0;
  REAL8 priorConfidence = 0.0;
  REAL8 priorUL = 0.0;
  REAL8 priorUncertainty = 0.0;
  /* REAL8 nAboveLE = 0.0;
  REAL8 nTrials = 0.0; */ /* 09/01/05 gam */
  REAL8 *nAboveLE = NULL;
  REAL8 *nTrials = NULL;
  REAL8 *h_0Rescaled = NULL;  
  INT4 iRescale = 0;
  INT4 iRescaleBest = 0;
  REAL8 *arrayULs = NULL;
  REAL8 *arrayConfs = NULL;
  INT4 *arrayConverged = NULL;
  REAL8 *arrayIteratedULs = NULL;   /* 09/01/05 gam */
  REAL8 *arrayIteratedConfs = NULL; /* 09/01/05 gam */ 
  INT4 iMC = 0;             /* which entire MC simulation is currently running */
  INT4 countMC = 0;         /* count of number of completed MC simulations, including rescalings */
  INT4 maxMC = 1;           /* default number of MC iterations to perform; set equal to params->maxMCinterations when interating MC */
  BOOLEAN breakOutOfLoop;   /* 09/01/05 gam; if interated MC converges or "runs off the tracks" break out of loop */
  INT4 onePlusNumMCRescalings = 1 + params->numMCRescalings; /* 09/01/05 gam; one plus number of times to rescale fake SFTs */
  INT4 twoPlusNumMCRescalings = 2 + params->numMCRescalings; /* 09/01/05 gam; one plus number of times to rescale fake SFTs */
  REAL8 maxMCErr = params->maxMCErr;  /* maximum absolute error in confidence to accept.  */
  REAL8 thisMCErr = 1;           /* keep track of last absolute error in confidence. */
  REAL8 lastMCErr = 1;           /* keep track of last absolute error in confidence. */
  INT4  *savSumBinMask;          /* 07/15/2005 gam */
  INT4  *arrayAvailableFreqBins; /* 07/15/2005 gam */
  INT4  nAvailableBins;          /* 07/15/2005 gam */
  REAL8 real8NAvailableBins;     /* 07/15/2005 gam */
  BOOLEAN searchSurroundingPts;  /* 07/29/05 gam; if (params->testFlag & 64) > 0 search surrounding parameters space pts */
  BOOLEAN narrowBLKandSTKband;   /* 09/12/05 gam; if (params->testFlag & 128) > 0 make BLKs and STKs narrower band based on extra bins */  
  INT4 jStart = 0;
  INT4 jSavedBLK = 0;
  INT4 extraBins = 0;
  REAL8 deltaBLKf0 = 0.0;
  REAL8 f0BLK, f0STK, bandBLK, bandSTK;
  INT4 nBinsPerBLK, nBinsPerSTK;
  
  INITSTATUS( status, "RunStackSlideIsolatedMonteCarloSimulation", STACKSLIDEISOLATEDC );
  ATTATCHSTATUSPTR(status);

  /* 07/29/05 gam; check that parameters are set to search SUM for loudest event; not search above threshold needed */
  if ( !( ((params->outputEventFlag & 2) > 0) && (params->thresholdFlag <= 0) ) ) {
      ABORT( status, STACKSLIDEISOLATEDH_EMCEVENTTHRESHOLDFLAGS, STACKSLIDEISOLATEDH_MSGEMCEVENTTHRESHOLDFLAGS);
  }

  /* 09/01/05 gam */
  if (params->numMCRescalings < 0) {
      ABORT( status, STACKSLIDEISOLATEDH_ENUMMCRESCALINGS, STACKSLIDEISOLATEDH_MSGENUMMCRESCALINGS);
  }
  h_0Rescaled = (REAL8 *)LALMalloc(onePlusNumMCRescalings*sizeof(REAL8));

  if ( (params->testFlag & 32) > 0 ) {
      /* 05/24/05 gam; iterate MC to converge on desired confidence */
      maxMC = params->maxMCinterations;
  }
  if (maxMC < 1) {
      ABORT( status, STACKSLIDEISOLATEDH_EMAXMC, STACKSLIDEISOLATEDH_MSGEMAXMC);
  }

  if ( (params->testFlag & 16) > 0 ) {
    /* 05/24/05 gam; use results from prior jobs in the pipeline and report on current MC results */
    reportMCResults = 1;

    if (maxMCErr <= 0.0) {
       ABORT( status, STACKSLIDEISOLATEDH_EMAXMCERR, STACKSLIDEISOLATEDH_MSGEMAXMCERR);
    }

    /* 09/01/05 gam; allocate memory dependent on the number of rescalings */
    nAboveLE = (REAL8 *)LALMalloc(onePlusNumMCRescalings*sizeof(REAL8));
    nTrials  = (REAL8 *)LALMalloc(onePlusNumMCRescalings*sizeof(REAL8));

    /* get results for previous jobs in the pipeline */
    getStackSlidePriorResults(status->statusPtr,
                               &priorLoudestEvent,
                               &priorStartFreq,
                               &priorBand,
                               &priorConfidence,
                               &priorUL,
                               &priorUncertainty,
                               params->priorResultsFile);
    CHECKSTATUSPTR (status);
    /* 10/20/05 gam; use params->blkRescaleFactor to rescale priorUL for internal use. */
    params->threshold4 = params->blkRescaleFactor*priorUL; /* Change injection amplitude to estimate from prior jobs in the pipeline */
    h_0 = params->threshold4;     /* Note that h_0 == params->threshold4 */
    #ifdef INCLUDE_PRINT_REPORTMCRESULTS_CODE
      if ((params->debugOptionFlag & 32) > 0 ) {
        fprintf(stdout,"Loudest Event, Start Frequency, Search Band, Confidence, Estimated Upper Limit, Uncertainty = \n");
        fprintf(stdout,"%20.6f %20.10f %20.10f %8.2f %20.10e %20.10e\n",priorLoudestEvent,priorStartFreq,priorBand,priorConfidence,priorUL,priorUncertainty);
        fflush(stdout);
      }
    #endif
    if (params->numMCRescalings > 0) {
      arrayULs   = (REAL8 *)LALMalloc((twoPlusNumMCRescalings*maxMC)*sizeof(REAL8));
      arrayConfs = (REAL8 *)LALMalloc((twoPlusNumMCRescalings*maxMC)*sizeof(REAL8));
      arrayConverged = (INT4 *)LALMalloc((twoPlusNumMCRescalings*maxMC)*sizeof(INT4));
    } else {
      arrayULs   = (REAL8 *)LALMalloc(maxMC*sizeof(REAL8));
      arrayConfs = (REAL8 *)LALMalloc(maxMC*sizeof(REAL8));
      arrayConverged = (INT4 *)LALMalloc(maxMC*sizeof(INT4));
    }
    arrayIteratedULs   = (REAL8 *)LALMalloc(maxMC*sizeof(REAL8));
    arrayIteratedConfs = (REAL8 *)LALMalloc(maxMC*sizeof(REAL8));
  } /* END if ( (params->testFlag & 16) > 0 ) */
  if ( (params->numMCRescalings > 0) && (reportMCResults != 1)) {
      ABORT( status, STACKSLIDEISOLATEDH_E2NUMMCRESCALINGS, STACKSLIDEISOLATEDH_MSGE2NUMMCRESCALINGS);
  }

  /* 07/29/05 gam; if testFlag if (params->testFlag & 64) > 0 search all parameters spaces points surrounding the injection; else search nearest only */
  if ( (params->testFlag & 64) > 0 ) {
       searchSurroundingPts = 1;  /* Search only the surrounding paremeters space point, in the Euclidean sense. */
  } else {
       searchSurroundingPts = 0;  /* Search only the nearest paremeters space point, in the Euclidean sense. */
  }

  /* 05/21/04 gam; save the input noise SFTs */
  savBLKData=(FFT **)LALMalloc(params->numBLKs*sizeof(FFT *));
  for (i=0;i<params->numBLKs;i++)
  {
        savBLKData[i]=(FFT *)LALMalloc(sizeof(FFT));
        savBLKData[i]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
        savBLKData[i]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
        savBLKData[i]->fft->data->data=(COMPLEX8 *)LALMalloc(params->nBinsPerBLK*sizeof(COMPLEX8));

        savBLKData[i]->fft->epoch = params->BLKData[i]->fft->epoch;
        savBLKData[i]->fft->f0 = params->BLKData[i]->fft->f0;
        savBLKData[i]->fft->deltaF = params->BLKData[i]->fft->deltaF;
        savBLKData[i]->fft->data->length = params->BLKData[i]->fft->data->length;
        for(j=0;j<params->nBinsPerBLK;j++) {
           savBLKData[i]->fft->data->data[j].re = params->BLKData[i]->fft->data->data[j].re;
           savBLKData[i]->fft->data->data[j].im = params->BLKData[i]->fft->data->data[j].im;
        }
  }  

  /* 09/12/05 gam; save the original values and check if using use narrow band BLKs and STKS */
  nBinsPerBLK = params->nBinsPerBLK; 
  nBinsPerSTK = params->nBinsPerSTK;
  f0BLK = params->f0BLK;
  f0STK = params->f0STK;
  bandBLK = params->bandBLK;
  bandSTK = params->bandSTK;
  if ( (params->testFlag & 128) > 0 ) {    
     extraBins = params->nBinsPerBLK - params->nBinsPerSUM;  /* The difference should be twice the maximum bandwidth of a signal during a run */
     extraBins = extraBins/2 + 1;  /* Want half the extra number of bins plus 1 for safety. */
     deltaBLKf0 = ((REAL8)(extraBins))/params->tEffBLK; /* Shift start frequency of BLKs and STKs from injected frequency by this much. */
     params->nBinsPerBLK = 2*extraBins;
     if ( (extraBins < 1) || (params->nBinsPerBLK > nBinsPerBLK) ) {
        params->nBinsPerBLK = nBinsPerBLK; /* reset this and just continue; should add ABORT */
        narrowBLKandSTKband = 0;
     } else {     
       narrowBLKandSTKband = 1;
       params->nBinsPerSTK = params->nBinsPerBLK;
       params->bandBLK = ((REAL8)(params->nBinsPerBLK))/params->tEffBLK;
       params->bandSTK = params->bandBLK;
       params->iMinBLK = 0;
       params->iMaxBLK = params->nBinsPerBLK - 1;

       /* reallocate BLK data */
       for (i=0;i<params->numBLKs;i++) {
          LALFree(params->BLKData[i]->fft->data->data);
          LALFree(params->BLKData[i]->fft->data);
          LALFree(params->BLKData[i]->fft);
          LALFree(params->BLKData[i]);
       }
       LALFree(params->BLKData);
       params->BLKData=(FFT **)LALMalloc(params->numBLKs*sizeof(FFT *));
       for (i=0;i<params->numBLKs;i++) {
           params->BLKData[i]=(FFT *)LALMalloc(sizeof(FFT));
           params->BLKData[i]->fft=(COMPLEX8FrequencySeries *)LALMalloc(sizeof(COMPLEX8FrequencySeries));
           params->BLKData[i]->fft->data=(COMPLEX8Vector *)LALMalloc(sizeof(COMPLEX8Vector));
           params->BLKData[i]->fft->data->data=(COMPLEX8 *)LALMalloc(params->nBinsPerBLK*sizeof(COMPLEX8));
           params->BLKData[i]->fft->epoch = savBLKData[i]->fft->epoch;
           params->BLKData[i]->fft->f0 = savBLKData[i]->fft->f0; /* This will get updated with each injection */
           params->BLKData[i]->fft->deltaF = savBLKData[i]->fft->deltaF;
           params->BLKData[i]->fft->data->length = params->nBinsPerBLK;
       }     
     }     
  } else {
     narrowBLKandSTKband = 0;
     jStart = 0;
  }
  
  /* 05/21/04 gam; save the input skyPosData */
  savSkyPosData=(REAL8 **)LALMalloc(params->numSkyPosTotal*sizeof(REAL8 *));
  for(i=0;i<params->numSkyPosTotal;i++)
  {
        savSkyPosData[i] = (REAL8 *)LALMalloc(2*sizeof(REAL8));
        savSkyPosData[i][0] = params->skyPosData[i][0];
        savSkyPosData[i][1] = params->skyPosData[i][1];
  }
  /* reallocate memory for the skyPosData structure */
  for(i=0;i<params->numSkyPosTotal;i++)
  {
      LALFree(params->skyPosData[i]);
  }
  LALFree(params->skyPosData);
  numSkyPosTotal = params->numSkyPosTotal;  /* sav original number of sky positions */
  /* 07/29/05 gam */
  if (searchSurroundingPts) {
    params->numSkyPosTotal = 4;               /* Set up to search 4 surrounding sky position */
  } else {
    params->numSkyPosTotal = 1;               /* Set up to run on nearest sky position */
  }
  params->skyPosData=(REAL8 **)LALMalloc(params->numSkyPosTotal*sizeof(REAL8 *));
  /* params->skyPosData[0] = (REAL8 *)LALMalloc(2*sizeof(REAL8)); */ /* 07/29/05 gam */
  for(i=0;i<params->numSkyPosTotal;i++) {
        params->skyPosData[i] = (REAL8 *)LALMalloc(2*sizeof(REAL8));
  }

  /* 05/21/04 gam; save the input freqDerivData and related parameters: */
  numParamSpacePts = params->numParamSpacePts;
  numSUMsTotal = params->numSUMsTotal;
  numFreqDerivTotal = params->numFreqDerivTotal;
  numFreqDerivIncludingNoSpinDown = params->numFreqDerivIncludingNoSpinDown;
  if (params->numSpinDown > 0) {
    savFreqDerivData=(REAL8 **)LALMalloc(params->numFreqDerivTotal*sizeof(REAL8 *));
    for(i=0;i<params->numFreqDerivTotal;i++)
    {
        savFreqDerivData[i] = (REAL8 *)LALMalloc(params->numSpinDown*sizeof(REAL8));
        for(j=0;j<params->numSpinDown;j++)
        {
          savFreqDerivData[i][j] = params->freqDerivData[i][j];
        }
    }
    /* reallocate memory for the freqDerivData structure */
    for(i=0;i<params->numFreqDerivTotal;i++)
    {
        LALFree(params->freqDerivData[i]);
    }
    LALFree(params->freqDerivData);
    /* 07/29/05 gam */
    if (searchSurroundingPts) {
       params->numFreqDerivTotal = (INT4)pow(2.0,((REAL8)params->numSpinDown)); /* Set up to search surrounding frequency deriviatives */
    } else {
       params->numFreqDerivTotal = 1;                                  /* Set up to search nearest frequency deriviatives */
    }
    params->freqDerivData=(REAL8 **)LALMalloc(params->numFreqDerivTotal*sizeof(REAL8 *));
    /* params->freqDerivData[0] = (REAL8 *)LALMalloc(params->numSpinDown*sizeof(REAL8)); */ /* 07/29/05 gam */
    for(i=0;i<params->numFreqDerivTotal;i++) {
      params->freqDerivData[i] = (REAL8 *)LALMalloc(params->numSpinDown*sizeof(REAL8));
    }    
  } /* END if (params->numSpinDown > 0) */
  /* 07/29/05 gam */
  if (params->numFreqDerivTotal != 0) {
    params->numFreqDerivIncludingNoSpinDown = params->numFreqDerivTotal;
  } else {
    params->numFreqDerivIncludingNoSpinDown = 1;  /* Even if numSpinDown = 0 still need to count case of zero spindown. */
  }
  /* params->numParamSpacePts = 1;
  params->numSUMsTotal = 1;  */ /* 07/29/05 gam */
  params->numParamSpacePts = params->numSkyPosTotal*params->numFreqDerivIncludingNoSpinDown;
  params->numSUMsTotal = params->numSUMsPerParamSpacePt*params->numParamSpacePts;
   
  /* Allocate memory for PulsarSignalParams and initialize */
  pPulsarSignalParams = (PulsarSignalParams *)LALMalloc(sizeof(PulsarSignalParams));
  pPulsarSignalParams->pulsar.position.system = COORDINATESYSTEM_EQUATORIAL;
  pPulsarSignalParams->pulsar.spindown = NULL;
  if (params->numSpinDown > 0) {
    LALDCreateVector(status->statusPtr, &(pPulsarSignalParams->pulsar.spindown),((UINT4)params->numSpinDown));
    CHECKSTATUSPTR (status);
  }
  pPulsarSignalParams->orbit = NULL; 
  pPulsarSignalParams->transfer = NULL;  
  /* Set up pulsar site */
  if (strstr(params->IFO, "LHO")) {
       cachedDetector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  } else if (strstr(params->IFO, "LLO")) {
       cachedDetector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  } else if (strstr(params->IFO, "GEO")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  } else if (strstr(params->IFO, "VIRGO")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  } else if (strstr(params->IFO, "TAMA")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  } else {
      /* "Invalid or null IFO" */
      ABORT( status, STACKSLIDEISOLATEDH_EIFO, STACKSLIDEISOLATEDH_MSGEIFO);
  }    
  pPulsarSignalParams->site = &cachedDetector;     
  pPulsarSignalParams->ephemerides = params->edat;
  /* pPulsarSignalParams->startTimeGPS = params->actualStartTime; */ /* 06/01/04 gam; uses gpsStartTime requested; not actualStartTime */
  pPulsarSignalParams->startTimeGPS.gpsSeconds = (INT4)params->gpsStartTimeSec;
  pPulsarSignalParams->startTimeGPS.gpsNanoSeconds = (INT4)params->gpsStartTimeNan;
  pPulsarSignalParams->duration = (UINT4)params->duration;
  pPulsarSignalParams->samplingRate = (REAL8)ceil(2.0*params->bandBLK); /* Make sampleRate an integer so that T*samplingRate = integer for integer T */
  pPulsarSignalParams->fHeterodyne = params->f0BLK;
  /* Find the time at the SSB that corresponds to the arrive time at the detector of first data requested. */
  GPSin.gpsSeconds = (INT4)params->gpsEpochStartTimeSec;     /* 12/06/04 gam; GPS epoch Sec that gives reference time in SSB */
  GPSin.gpsNanoSeconds = (INT4)params->gpsEpochStartTimeNan; /* 12/06/04 gam; GPS epoch Nan that gives reference time in SSB */
  
  /* Allocate memory for SFTParams and initialize */
  pSFTParams = (SFTParams *)LALMalloc(sizeof(SFTParams));
  pSFTParams->Tsft = params->tBLK;
  timestamps.length = params->numBLKs;
  timestamps.data = params->timeStamps; 
  pSFTParams->timestamps = &timestamps;
  pSFTParams->noiseSFTs = NULL;  
  
  /* 05/26/04 gam; initialize variables that keep track of which SUM we are working on */
  params->startSUMs = 1;   /* 05/26/04 gam; use to control I/O during Monte Carlo. Default is TRUE. */
  params->finishPeriodicTable = 0; /* 05/24/05 gam */ /* 05/26/04 gam; use to control I/O during Monte Carlo. Default is TRUE. */
  params->whichMCSUM = -1; /* 05/26/04 gam; which SUM the Monte Carlo Simulation is running on. Default is -1 */
  numSUMsTotalm1 = numSUMsTotal - 1;  /* Index of last SUM */

  /* 05/26/04 gam; initialize variables that keep track of which frequency we are working on */  
  f0SUM = params->f0SUM;
  bandSUM = params->bandSUM;
  nBinsPerSUM = params->nBinsPerSUM;
  real8NBinsPerSUM = ((REAL8)nBinsPerSUM); /* 07/15/05 gam */
  nBinsPerSUMm1 = nBinsPerSUM - 1;
  params->keepThisNumber = 1;      /* During MC only keep 1 event; this will be the injected event */

  /* 07/29/05 gam */
  if (searchSurroundingPts) {
      params->bandSUM = 2.0*params->dfSUM;
      params->nBinsPerSUM = 2;              /* Search surrounding frequencies */
  } else {
      params->bandSUM = params->dfSUM;      /* Search nearest frequencies */
      params->nBinsPerSUM = 1;
  }

  /* 05/28/04 gam; Initial seed and randPar to use LALUniformDeviate to generate random mismatch during Monte Carlo. */
  fpRandom = fopen("/dev/urandom","r");
  rndCount = fread(&seed, sizeof(INT4),1, fpRandom);
  fclose(fpRandom);
  /* seed = 1234; */ /* Test value */
  LALCreateRandomParams(status->statusPtr, &randPar, seed);
  CHECKSTATUSPTR (status);
  
  /* 07/14/04 gam; allocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  if ( (params->testFlag & 4) > 0 ) {
     #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
          fprintf(stdout,"allocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs \n");
          fflush(stdout);
     #endif
     pSkyConstAndZeroPsiAMResponse = (SkyConstAndZeroPsiAMResponse *)LALMalloc(sizeof(SkyConstAndZeroPsiAMResponse));
     pSkyConstAndZeroPsiAMResponse->skyConst = (REAL8 *)LALMalloc((2*params->numSpinDown*(params->numBLKs+1)+2*params->numBLKs+3)*sizeof(REAL8));
     pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi = (REAL4 *)LALMalloc(params->numBLKs*sizeof(REAL4));
     pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi = (REAL4 *)LALMalloc(params->numBLKs*sizeof(REAL4));
     pSFTandSignalParams = (SFTandSignalParams *)LALMalloc(sizeof(SFTandSignalParams));
     /* create lookup table (LUT) values for doing trig */
     /* pSFTandSignalParams->resTrig = 64; */ /* length sinVal and cosVal; resolution of trig functions = 2pi/resTrig */
     pSFTandSignalParams->resTrig = 0; /* 08/02/04 gam; avoid serious bug when using LUTs for trig calls */
     pSFTandSignalParams->trigArg = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
     pSFTandSignalParams->sinVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
     pSFTandSignalParams->cosVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
     for (k=0; k<=pSFTandSignalParams->resTrig; k++) {
       pSFTandSignalParams->trigArg[k]= ((REAL8)LAL_TWOPI) * ((REAL8)k) / ((REAL8)pSFTandSignalParams->resTrig);
       pSFTandSignalParams->sinVal[k]=sin( pSFTandSignalParams->trigArg[k] );
       pSFTandSignalParams->cosVal[k]=cos( pSFTandSignalParams->trigArg[k] );
     }
     pSFTandSignalParams->pSigParams = pPulsarSignalParams;
     pSFTandSignalParams->pSFTParams = pSFTParams;
     pSFTandSignalParams->nSamples = nSamples;
     pSFTandSignalParams->Dterms = params->Dterms; /* 09/08/05 gam */
     pPulsarSignalParams->samplingRate = 2.0*params->bandBLK; /* can set samplingRate to exactly 2.0*bandBLK when using LALFastGeneratePulsarSFTs */
     GPSin.gpsSeconds = params->timeStamps[0].gpsSeconds;         /* 08/02/04 gam */
     GPSin.gpsNanoSeconds = params->timeStamps[0].gpsNanoSeconds; /* 08/02/04 gam */
     TRY (LALCreateSFTVector (status->statusPtr, &outputSFTs, params->numBLKs, params->nBinsPerBLK), status); /* 09/08/05 gam */
  } else { 
     pSkyConstAndZeroPsiAMResponse = NULL;
     pSFTandSignalParams = NULL;
  }

  /* 05/19/05 gam; find first and last nonexcluded frequency bins */
  savSumBinMask=(INT4 *)LALMalloc(nBinsPerSUM*sizeof(INT4)); /* 07/15/2005 gam */
  for(iFreq=0;iFreq<nBinsPerSUM;iFreq++) {
     savSumBinMask[iFreq] = params->sumBinMask[iFreq]; /* 07/15/2005 gam; save this value */
     /* 05/19/05 gam; Only inject into bins that are not exclude */
     if (params->sumBinMask[iFreq] != 0) {
        if (firstNonExcludedBin == -1) {
           firstNonExcludedBin = iFreq;
        }
        lastNonExcludedBin = iFreq;
     }
  }
  LALFree(params->sumBinMask); /* 07/15/2005 gam; reallocate and save just 1 */
  /* 07/29/05 gam */
  if (searchSurroundingPts) {
    params->sumBinMask=(INT4 *)LALMalloc(2*sizeof(INT4));
    params->sumBinMask[0] = 1;   /* Default value; will check when injected frequency is set up below */
    params->sumBinMask[1] = 1;   /* Default value; will check when injected frequency is set up below */
  } else {
    params->sumBinMask=(INT4 *)LALMalloc(sizeof(INT4));
    params->sumBinMask[0] = 1;   /* Will only inject next bins where savSumBinMask != 0. */
  }

  /* 07/15/2005 gam; initialize array with available frequency; i.e., bins not excluded due to instrument lines */ 
  arrayAvailableFreqBins=(INT4 *)LALMalloc(nBinsPerSUM*sizeof(INT4));
  i = 0;
  for(iFreq=0;iFreq<nBinsPerSUM;iFreq++) {
     arrayAvailableFreqBins[iFreq] = iFreq; /* for safety initialize with default value; should never access this array past nAvailableBins - 1 */
     if (savSumBinMask[iFreq] != 0) {
        arrayAvailableFreqBins[i] = iFreq;
        i++;
     }
  }
  nAvailableBins = i;
  if (nAvailableBins <= 0) {
     /* 07/29/05 gam; should ABORT if cleaning or avoiding all freq bins */
     ABORT( status, STACKSLIDEISOLATEDH_ENOFREQBINS, STACKSLIDEISOLATEDH_MSGENOFREQBINS);
  }
  real8NAvailableBins = ((REAL8)nAvailableBins);

 /***********************************************************/
 /*                                                         */
 /* START SECTION: LOOP OVER ENTIRE MONTE CARLO SIMULATIONS */
 /*                                                         */
 /***********************************************************/
 for(iMC=0;iMC<maxMC;iMC++) {
   
  /* This is the big loop that reruns the entire MC maxMC times or until MC has converged */
  
  /* initialize at beginning of each MC */
  /* nAboveLE = 0.0;
  nTrials = 0.0; */ /* 09/01/05 gam */
  if (reportMCResults) {  
    for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) {
      nAboveLE[iRescale] = 0.0;
      nTrials[iRescale] = 0.0;
    }
  }
  
  /* 07/15/05 gam; change Monte Carlo Simulation to inject random signals from the parameter space */
  /**********************************************/
  /*                                            */
  /* START SECTION: MONTE CARLO INJECTIONS LOOP */
  /*                                            */
  /**********************************************/
  for(kSUM=0;kSUM<params->numMCInjections;kSUM++) {

    params->whichMCSUM = kSUM; /* kSUM is which injection we are working on */  

    /* 07/15/05 gam; Find random sky position: RA and Dec */
    if (rangeRA != 0.0) {
       StackSlideGetUniformDeviate(status->statusPtr, &tmpRA, startRA, rangeRA, randPar); CHECKSTATUSPTR (status);
    } else if (DeltaRA != 0.0) {
       StackSlideGetUniformDeviate(status->statusPtr, &tmpRA, startRA-(0.5*DeltaRA), DeltaRA, randPar); CHECKSTATUSPTR (status);
    } else {
       tmpRA = startRA;  /* Use fixed RA */
    }
    if (rangeDec != 0.0) {
       /* StackSlideGetUniformDeviate(status->statusPtr, &tmpDec, startDec, rangeDec, randPar); CHECKSTATUSPTR (status); */ /* 09/05/05 gam */
       StackSlideGetUniformDeviate(status->statusPtr, &tmpDec, sinStartDec, rangeDec, randPar); CHECKSTATUSPTR (status);
       tmpDec = asin(tmpDec); /* 09/05/05 gam; found uniformly distributed sin(DEC); need DEC in radians */
    } else if (DeltaDec != 0.0) {
       /* StackSlideGetUniformDeviate(status->statusPtr, &tmpDec, startDec-(0.5*DeltaDec), DeltaDec, randPar); CHECKSTATUSPTR (status); */ /* 09/05/05 gam */
       StackSlideGetUniformDeviate(status->statusPtr, &tmpDec, sinStartDecmhalfDeltaDec, rangeDecStartDecpmhDeltaDec, randPar); CHECKSTATUSPTR (status);
       tmpDec = asin(tmpDec); /* 09/05/05 gam; found uniformly distributed sin(DEC); need DEC in radians */
    } else {
       tmpDec = startDec;  /* Use fixed Dec */
    }
    
    /* Set up params->skyPosData with position(s) nearest injected sky position */
    if (DeltaDec != 0.0) {
       iDec = floor( (tmpDec - startDec)/DeltaDec + 0.5 );
    } else {
       iDec = 0;
    }
    params->skyPosData[0][1] = startDec + ((REAL8)iDec)*DeltaDec;
    cosTmpDEC = cos(params->skyPosData[0][1]);
    if (cosTmpDEC != 0.0) {
          tmpDeltaRA = DeltaRA/cosTmpDEC;
          /* 09/16/05 gam */
          if ( (tmpDeltaRA != 0.0) && (stopRA > startRA) ) {
            tmpNumRA = ceil((stopRA - startRA)/tmpDeltaRA);
            tmpDeltaRA = (stopRA - startRA)/((REAL8)tmpNumRA);
          } else {
            tmpDeltaRA = 0.0;
          }
    } else {
          tmpDeltaRA = 0.0; /* We are at a celestial pole */
    }
    if (tmpDeltaRA != 0.0) {
       iRA = floor( (tmpRA - startRA)/tmpDeltaRA + 0.5 );
    } else {
       iRA = 0;
    }
    params->skyPosData[0][0] = startRA + ((REAL8)iRA)*tmpDeltaRA;

    /* 07/29/05 gam; find 3 other surrounding sky positions */
    if (searchSurroundingPts) {
       /* !!! WARNING: Signs are chooses such that DeltaDec and DeltaRA are always positive !!! */
       /* Go to next nearest RA; keep DEC the same */
       if (params->skyPosData[0][0] < tmpRA) {
          params->skyPosData[1][0] = params->skyPosData[0][0] + tmpDeltaRA;
       } else {
          params->skyPosData[1][0] = params->skyPosData[0][0] - tmpDeltaRA;
       }
       params->skyPosData[1][1] = params->skyPosData[0][1];
       
       /* Go to next nearest DEC; keep RA the same */
       if (params->skyPosData[0][1] < tmpDec) {
          params->skyPosData[2][1] = params->skyPosData[0][1] + DeltaDec;
       } else {
          params->skyPosData[2][1] = params->skyPosData[0][1] - DeltaDec;
       }
       params->skyPosData[2][0] = params->skyPosData[0][0];
       
       /* Find tmpDeltaRA at next nearest DEC: */
       cosTmpDEC = cos(params->skyPosData[2][1]);
       if (cosTmpDEC != 0.0) {
          tmpDeltaRA = DeltaRA/cosTmpDEC;
          /* 09/16/05 gam */
          if ( (tmpDeltaRA != 0.0) && (stopRA > startRA) ) {
            tmpNumRA = ceil((stopRA - startRA)/tmpDeltaRA);
            tmpDeltaRA = (stopRA - startRA)/((REAL8)tmpNumRA);
          } else {
            tmpDeltaRA = 0.0;
          }
       } else {
          tmpDeltaRA = 0.0; /* We are at a celestial pole */
       }
       
       /* Change both DEC and RA to find last corner of square that surrounds injected sky position */
       params->skyPosData[3][1] = params->skyPosData[2][1]; /* same DEC as skyPosData[2] */
       if (params->skyPosData[2][0] < tmpRA) {
          params->skyPosData[3][0] = params->skyPosData[2][0] + tmpDeltaRA;
       } else {
          params->skyPosData[3][0] = params->skyPosData[2][0] - tmpDeltaRA;
       }
    } /* END if (searchSurroundingPts) */
    
    /* Check if sky frame is rotated and act accordingly */
    if ( (params->parameterSpaceFlag & 1) >  0) {
      /* rotate skyPosData into coordinates with Earth's average acceleration at the pole. */
      RotateSkyCoordinates(status->statusPtr, &tmpRA, &tmpDec, tmpRA, tmpDec, params->aveEarthAccRA, params->aveEarthAccDEC, 0.0);
      CHECKSTATUSPTR (status);
      RotateSkyPosData(status->statusPtr, params->skyPosData, params->numSkyPosTotal, params->aveEarthAccRA, params->aveEarthAccDEC, 0.0);
      CHECKSTATUSPTR (status);
    } else if ( (params->parameterSpaceFlag & 2) >  0) {
      /* rotate skyPosData into galactic plane */
      RotateSkyCoordinates(status->statusPtr, &tmpRA, &tmpDec, tmpRA, tmpDec, 4.65, ((REAL8)LAL_PI_2)-0.4, 0.0);
      CHECKSTATUSPTR (status);
      RotateSkyPosData(status->statusPtr, params->skyPosData, params->numSkyPosTotal, 4.65,((REAL8)LAL_PI_2)-0.4, 0.0);
      CHECKSTATUSPTR (status);
    }
    pPulsarSignalParams->pulsar.position.longitude = tmpRA;
    pPulsarSignalParams->pulsar.position.latitude = tmpDec;

    #ifdef DEBUG_RANDOMTRIALPARAMETERS
        fprintf(stdout,"\nkSUM = %i, search  RA[0] = %23.10e, inject RA = %23.10e \n",kSUM,params->skyPosData[0][0],pPulsarSignalParams->pulsar.position.longitude);
        fprintf(stdout,"kSUM = %i, search DEC[0] = %23.10e, inject DEC = %23.10e \n",kSUM,params->skyPosData[0][1],pPulsarSignalParams->pulsar.position.latitude);
        if (searchSurroundingPts) {
          fprintf(stdout,"search  RA[1] = %23.10e \n",params->skyPosData[1][0]);
          fprintf(stdout,"search DEC[1] = %23.10e \n",params->skyPosData[1][1]);
          fprintf(stdout,"search  RA[2] = %23.10e \n",params->skyPosData[2][0]);
          fprintf(stdout,"search DEC[2] = %23.10e \n",params->skyPosData[2][1]);
          fprintf(stdout,"search  RA[3] = %23.10e \n",params->skyPosData[3][0]);
          fprintf(stdout,"search DEC[3] = %23.10e \n",params->skyPosData[3][1]);
        }
        fflush(stdout);      
    #endif

    if (params->numSpinDown > 0) {
      for(k=0;k<params->numSpinDown;k++) {
        if (k == 0) {
          if (rangeFDeriv1 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv1, rangeFDeriv1, randPar); CHECKSTATUSPTR (status);
          } else if (DeltaFDeriv1 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv1-(0.5*DeltaFDeriv1), DeltaFDeriv1, randPar); CHECKSTATUSPTR (status);
          } else {
            real8RandVal = startFDeriv1;  /* Use fixed value */
          }
          pPulsarSignalParams->pulsar.spindown->data[k] = real8RandVal;
          /* Set up params->freqDerivData with deriv(s) nearest injected deriv */
          if (DeltaFDeriv1 != 0.0) {
             iDeriv = floor( (real8RandVal - startFDeriv1)/DeltaFDeriv1 + 0.5 );
          } else {
             iDeriv = 0;
          }
          params->freqDerivData[0][k] = startFDeriv1 + ((REAL8)iDeriv)*DeltaFDeriv1;
        } else if (k == 1) {
          if (rangeFDeriv2 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv2, rangeFDeriv2, randPar); CHECKSTATUSPTR (status);
          } else if (DeltaFDeriv2 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv2-(0.5*DeltaFDeriv2), DeltaFDeriv2, randPar); CHECKSTATUSPTR (status);
          } else {
            real8RandVal = startFDeriv2;  /* Use fixed value */
          }
          pPulsarSignalParams->pulsar.spindown->data[k] = real8RandVal;
          /* Set up params->freqDerivData with deriv(s) nearest injected deriv */
          if (DeltaFDeriv2 != 0.0) {
             iDeriv = floor( (real8RandVal - startFDeriv2)/DeltaFDeriv2 + 0.5 );
          } else {
             iDeriv = 0;
          }
          params->freqDerivData[0][k] = startFDeriv2 + ((REAL8)iDeriv)*DeltaFDeriv2;
        } else if (k == 2) {
          if (rangeFDeriv3 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv3, rangeFDeriv3, randPar); CHECKSTATUSPTR (status);
          } else if (DeltaFDeriv3 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv3-(0.5*DeltaFDeriv3), DeltaFDeriv3, randPar); CHECKSTATUSPTR (status);
          } else {
            real8RandVal = startFDeriv3;  /* Use fixed value */
          }
          pPulsarSignalParams->pulsar.spindown->data[k] = real8RandVal;
          /* Set up params->freqDerivData with deriv(s) nearest injected deriv */
          if (DeltaFDeriv3 != 0.0) {
             iDeriv = floor( (real8RandVal - startFDeriv3)/DeltaFDeriv3 + 0.5 );
          } else {
             iDeriv = 0;
          }
          params->freqDerivData[0][k] = startFDeriv3 + ((REAL8)iDeriv)*DeltaFDeriv3;
        } else if (k == 3) {
          if (rangeFDeriv4 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv4, rangeFDeriv4, randPar); CHECKSTATUSPTR (status);
          } else if (DeltaFDeriv4 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv4-(0.5*DeltaFDeriv4), DeltaFDeriv4, randPar); CHECKSTATUSPTR (status);
          } else {
            real8RandVal = startFDeriv4;  /* Use fixed value */
          }
          pPulsarSignalParams->pulsar.spindown->data[k] = real8RandVal;
          /* Set up params->freqDerivData with deriv(s) nearest injected deriv */
          if (DeltaFDeriv4 != 0.0) {
             iDeriv = floor( (real8RandVal - startFDeriv4)/DeltaFDeriv4 + 0.5 );
          } else {
             iDeriv = 0;
          }
          params->freqDerivData[0][k] = startFDeriv4 + ((REAL8)iDeriv)*DeltaFDeriv4;
        } else if (k == 4) {
          if (rangeFDeriv5 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv5, rangeFDeriv5, randPar); CHECKSTATUSPTR (status);
          } else if (DeltaFDeriv5 != 0.0) {
            StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, startFDeriv5-(0.5*DeltaFDeriv5), DeltaFDeriv5, randPar); CHECKSTATUSPTR (status);
          } else {
            real8RandVal = startFDeriv5;  /* Use fixed value */
          }
          pPulsarSignalParams->pulsar.spindown->data[k] = real8RandVal;
          /* Set up params->freqDerivData with deriv(s) nearest injected deriv */
          if (DeltaFDeriv5 != 0.0) {
             iDeriv = floor( (real8RandVal - startFDeriv5)/DeltaFDeriv5 + 0.5 );
          } else {
             iDeriv = 0;
          }
          params->freqDerivData[0][k] = startFDeriv5 + ((REAL8)iDeriv)*DeltaFDeriv5;
        } /* END if (k == 0) ELSE ... */

        #ifdef DEBUG_RANDOMTRIALPARAMETERS
          fprintf(stdout,"kSUM = %i, search fDeriv[%i] = %23.10e, inject fDeriv[%i] = %23.10e \n",kSUM,k,params->freqDerivData[0][k],k,pPulsarSignalParams->pulsar.spindown->data[k]);
          fflush(stdout);
        #endif
      } /* END for(k=0;k<params->numSpinDown;k++) */

      /* 07/29/05 gam; find other surrounding frequency derivatives */
      if (searchSurroundingPts) {
         /* Note i starts at 1; add in surrounding points */
         for(i=1;i<params->numFreqDerivTotal;i++) {
            for(k=0;k<params->numSpinDown;k++) {
                /* !!! WARNING: Signs are chooses such that deltaFDeriv is negative for odd derivatives; positive for even derivative !!! */
                if (k == 0) {
                   if ( (i % 2) > 0 ) {
                     if (params->freqDerivData[0][k] < pPulsarSignalParams->pulsar.spindown->data[k]) {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] - params->deltaFDeriv1;
                     } else {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] + params->deltaFDeriv1;
                     }
                   } else {
                     params->freqDerivData[i][k] = params->freqDerivData[0][k];
                   }
                } else if (k == 1) {
                   if ( floor((i % 4)/2) > 0 ) {
                     if (params->freqDerivData[0][k] < pPulsarSignalParams->pulsar.spindown->data[k]) {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] + params->deltaFDeriv2;
                     } else {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] - params->deltaFDeriv2;
                     }
                   } else {
                     params->freqDerivData[i][k] = params->freqDerivData[0][k];
                   }
                } else if (k == 2) {
                   if ( floor((i % 8)/4) > 0 ) {
                     if (params->freqDerivData[0][k] < pPulsarSignalParams->pulsar.spindown->data[k]) {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] - params->deltaFDeriv3;
                     } else {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] + params->deltaFDeriv3;
                     }
                   } else {
                     params->freqDerivData[i][k] = params->freqDerivData[0][k];
                   }
                } else if (k == 3) {
                   if ( floor((i % 16)/8) > 0 ) {
                     if (params->freqDerivData[0][k] < pPulsarSignalParams->pulsar.spindown->data[k]) {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] + params->deltaFDeriv4;
                     } else {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] - params->deltaFDeriv4;
                     }
                   } else {
                     params->freqDerivData[i][k] = params->freqDerivData[0][k];
                   }
                } else if (k == 4) {
                   if ( floor((i % 32)/16) > 0 ) {
                     if (params->freqDerivData[0][k] < pPulsarSignalParams->pulsar.spindown->data[k]) {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] - params->deltaFDeriv5;
                     } else {
                       params->freqDerivData[i][k] = params->freqDerivData[0][k] + params->deltaFDeriv5;
                     }
                   } else {
                     params->freqDerivData[i][k] = params->freqDerivData[0][k];
                   }
                } /* END if (k == 0) ELSE ... */
                #ifdef DEBUG_RANDOMTRIALPARAMETERS
                   fprintf(stdout,"next fDeriv[%i][%i] = %23.10e\n",i,k,params->freqDerivData[i][k]);
                   fflush(stdout);
                #endif
            } /* for(k=0;k<params->numSpinDown;k++) */
         } /* END for(i=1;i<params->numFreqDerivTotal;i++) */
      } /* END if (searchSurroundingPts) */

    } /* END if (params->numSpinDown > 0) */

    /* update the reference time for the injected sky position */
    /* LALConvertGPS2SSB(status->statusPtr,&(pPulsarSignalParams->pulsar.tRef), GPSin, pPulsarSignalParams);
    CHECKSTATUSPTR (status);  */ /* 08/31/05 gam */
    pPulsarSignalParams->pulsar.tRef.gpsSeconds = GPSin.gpsSeconds;
    pPulsarSignalParams->pulsar.tRef.gpsNanoSeconds = GPSin.gpsNanoSeconds;
    
    /* fill in SkyConstAndZeroPsiAMResponse for use with LALFastGeneratePulsarSFTs */
    if ( (params->testFlag & 4) > 0 ) {
        #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
            fprintf(stdout,"about to call LALComputeSkyAndZeroPsiAMResponse \n");
            fflush(stdout);
        #endif
        LALComputeSkyAndZeroPsiAMResponse (status->statusPtr, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
        CHECKSTATUSPTR (status);
    }

    /* FIND RANDOM FREQUENCY THAT IS NOT EXCLUDED */
    StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, 0.0, real8NAvailableBins, randPar); CHECKSTATUSPTR (status);
    i = floor( real8RandVal );    /* nearest bin arrayAvailableFreqBins */
    if (i >= nAvailableBins) i = nAvailableBins - 1; if (i < 0) i = 0; /* make sure i in range */
    iFreq = arrayAvailableFreqBins[i];  /* randomly chosen available frequency bin */
    /* if (savSumBinMask[iFreq] == 0) probably should ABORT ! */
    params->f0SUM = f0SUM + iFreq*params->dfSUM;
    StackSlideGetUniformDeviate(status->statusPtr, &real8RandVal, params->f0SUM-(0.5*params->dfSUM), params->dfSUM, randPar); CHECKSTATUSPTR (status);
    pPulsarSignalParams->pulsar.f0 = real8RandVal; /* injected frequency including random mismatch */
    
    /* 07/29/05 gam; search is over 2 freq bins; find the other bin */
    if (searchSurroundingPts) {    
        if (params->f0SUM < pPulsarSignalParams->pulsar.f0) {
           /* Check if iFreq + 1 is a valid frequency bin */
           if ( (iFreq+1) < nBinsPerSUM ) {
              /* Set up to search the closest bin and the one AFTER it */
              params->bandSUM = 2.0*params->dfSUM;
              params->nBinsPerSUM = 2;
              params->sumBinMask[0] = 1; /* By contruction above it is an available bin */
              params->sumBinMask[1] = savSumBinMask[iFreq+1];
           } else {
              /* Only the closest bin is in the search parameter space */
              params->bandSUM = params->dfSUM;
              params->nBinsPerSUM = 1;
              params->sumBinMask[0] = 1; /* By contruction above it is an available bin */
              params->sumBinMask[1] = 0; 
           }
        } else {
           /* Check if iFreq - 1 is a valid frequency bin */
           if ( (iFreq-1) >= 0 ) {
              /* Set up to search the closest bin and the one BEFORE it */
              params->bandSUM = 2.0*params->dfSUM;
              params->nBinsPerSUM = 2;
              params->f0SUM = params->f0SUM - params->dfSUM; /* start search one bin before injected freq */
              params->sumBinMask[0] = savSumBinMask[iFreq-1];
              params->sumBinMask[1] = 1; /* By contruction above it is an available bin */
           } else {
              /* Only the closest bin is in the search parameter space */
              params->bandSUM = params->dfSUM;
              params->nBinsPerSUM = 1;
              params->sumBinMask[0] = 1; /* By contruction above it is an available bin */
              params->sumBinMask[1] = 0;
           }
        } /* END if (params->f0SUM < pPulsarSignalParams->pulsar.f0) */
    } /* END if (searchSurroundingPts) */

    /* 09/12/05 gam */
    if (narrowBLKandSTKband) {
       params->f0BLK = params->f0SUM - deltaBLKf0;
       params->f0STK = params->f0BLK;
       pPulsarSignalParams->fHeterodyne = params->f0BLK;
       /* what is the index to params->f0BLK in the savBLKs */
       jStart = floor((params->f0BLK - f0BLK)*params->tEffBLK + 0.5);  /* else jStart is initialized to zero above */
    }
    #ifdef DEBUG_RANDOMTRIALPARAMETERS
                   fprintf(stdout,"narrowBLKandSTKband = %i, extraBins = %i, params->f0BLK = %23.10e, deltaBLKf0 = %23.10e, params->f0SUM = = %23.10e, pPulsarSignalParams->pulsar.f0 = %23.10e\n",narrowBLKandSTKband,extraBins,params->f0BLK,deltaBLKf0,params->f0SUM,pPulsarSignalParams->pulsar.f0);
                   fprintf(stdout,"jStart = %i,\n",jStart);
                   fflush(stdout);
    #endif

    /* 12/06/04 gam */
    if ( (params->testFlag & 8) > 0 ) {
         pPulsarSignalParams->pulsar.psi = params->orientationAngle;
         cosIota =params->cosInclinationAngle;
    } else {
         /* get random value for psi */
         LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
         #ifdef DEBUG_SETFIXED_RANDVAL
            randval = params->threshold5; /* Temporarily use threshold5 for this; need to add testParameters to commandline. */
         #endif
         pPulsarSignalParams->pulsar.psi = (randval - 0.5) * ((REAL4)LAL_PI_2);

         /* get random value for cosIota */
         LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
         #ifdef DEBUG_SETFIXED_RANDVAL
            randval = params->threshold5; /* Temporarily use threshold5 for this; need to add testParameters to commandline. */
         #endif
         cosIota = 2.0*((REAL8)randval) - 1.0;
    } /* END if ( (params->testFlag & 8) > 0 ) */
    /* get A_+ and A_x from h_0 and cosIota */
    pPulsarSignalParams->pulsar.aPlus = (REAL4)(0.5*h_0*(1.0 + cosIota*cosIota));
    pPulsarSignalParams->pulsar.aCross = (REAL4)(h_0*cosIota);

    /* Last: get random value for phi0 */
    LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status); /* 05/28/04 gam */
    pPulsarSignalParams->pulsar.phi0 = ((REAL8)randval) * ((REAL8)LAL_TWOPI);

    #ifdef DEBUG_RANDOMTRIALPARAMETERS
        fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",iFreq,h_0);
        fprintf(stdout,"savSumBinMask[%i] = %i \n",iFreq,savSumBinMask[iFreq]);
        fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",iFreq,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
        fprintf(stdout,"iFreq = %i, inject psi = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.psi);
        fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.phi0);
        fprintf(stdout,"iFreq = %i, search start f0 = %23.10e, closest f0 = %23.10e, inject f0 = %23.10e \n",iFreq,params->f0SUM,f0SUM + iFreq*params->dfSUM,pPulsarSignalParams->pulsar.f0);
        fflush(stdout);
    #endif
    #ifdef DEBUG_MONTECARLOTIMEDOMAIN_DATA
        fprintf(stdout,"iFreq = %i, params->f0SUM = %23.10e, pPulsarSignalParams->pulsar.f0 = %23.10e \n",iFreq,params->f0SUM,pPulsarSignalParams->pulsar.f0);
        fprintf(stdout,"nBinsPerSUM = %i, params->nBinsPerSUM = %i\n",nBinsPerSUM,params->nBinsPerSUM);
        fflush(stdout);	
    #endif

    /* 07/14/04 gam; check if using LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
    if ( (params->testFlag & 4) > 0 ) {
        #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
          fprintf(stdout,"About to call LALFastGeneratePulsarSFTs \n");
          fflush(stdout);
        #endif

        /* 09/08/05 gam; LALFastGeneratePulsarSFTs only fills in data in a band 2*Dterms wide, for safety reinitialize. */
        for(i=0;i<params->numBLKs;i++) {
             for(j=0;j<params->nBinsPerBLK;j++) {
                 outputSFTs->data[i].data->data[j].re = 0.0;
                 outputSFTs->data[i].data->data[j].im = 0.0;
               }
        }

        /* 07/14/04 gam; use SkyConstAndZeroPsiAMResponse from LALComputeSkyAndZeroPsiAMResponse and SFTandSignalParams to generate SFTs fast. */
        LALFastGeneratePulsarSFTs (status->statusPtr, &outputSFTs, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
        CHECKSTATUSPTR (status);

        #ifdef PRINT_ONEMONTECARLO_OUTPUTSFT
          if (params->whichMCSUM == 0) {
            REAL4  fPlus;
            REAL4  fCross;
            i=0; /* index of which outputSFT to output */
            fPlus = pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi[i]*cos(pPulsarSignalParams->pulsar.psi) + pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi[i]*sin(pPulsarSignalParams->pulsar.psi);
            fCross = pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi[i]*cos(pPulsarSignalParams->pulsar.psi) - pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi[i]*sin(pPulsarSignalParams->pulsar.psi);
            fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",iFreq,h_0);
            fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",iFreq,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
            fprintf(stdout,"iFreq = %i, inject psi = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.psi);
            fprintf(stdout,"iFreq = %i, fPlus, fCross = %23.10e,  %23.10e \n",iFreq,fPlus,fCross);
            fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.phi0);
            fprintf(stdout,"iFreq = %i, search f0 = %23.10e, inject f0 = %23.10e \n",iFreq,f0SUM + iFreq*params->dfSUM,pPulsarSignalParams->pulsar.f0);
            fprintf(stdout,"outputSFTs->data[%i].data->length = %i \n",i,outputSFTs->data[i].data->length);
            fflush(stdout);
            for(j=0;j<params->nBinsPerBLK;j++) {
               /* fprintf(stdout,"%i %g %g \n",j,outputSFTs->data[i].data->data[j].re,outputSFTs->data[i].data->data[j].im); */
               fprintf(stdout,"%i %g \n",j,outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im);
               fflush(stdout);
            }
          }
        #endif

        #ifdef PRINT_MAXPOWERANDBINEACHSFT
         {
           INT4 jMaxPwr;
           REAL4 maxPwr;
           REAL4 pwr;
           for(i=0;i<params->numBLKs;i++) {
             maxPwr = 0.0;
             jMaxPwr = -1;
             for(j=0;j<params->nBinsPerBLK;j++) {
               pwr = outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im;
               if (pwr > maxPwr) {
                  maxPwr = pwr;
                  jMaxPwr = j;
               }
             }
             fprintf(stdout,"Max power for SFT %i is in bin %i = %g \n",i,jMaxPwr,maxPwr);
             fflush(stdout);
           }
         }
        #endif

        /* Add outputSFTs with injected signal to input noise SFTs; no renorm should be needed */
        for(i=0;i<params->numBLKs;i++) {
          params->BLKData[i]->fft->f0 = params->f0BLK; /* update with each injection */
          for(j=0;j<params->nBinsPerBLK;j++) {
             jSavedBLK = j + jStart; /* 09/12/05 gam */
             #ifdef PRINTCOMPARISON_INPUTVSMONTECARLOSFT_DATA
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].re, outputSFTs->data[%i].data->data[%i].re = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[jSavedBLK].re,outputSFTs->data[i].data->data[j].re);
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].im, outputSFTs->data[%i].data->data[%i].im = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[jSavedBLK].im,outputSFTs->data[i].data->data[j].im);
               fflush(stdout);
             #endif
             params->BLKData[i]->fft->data->data[j].re = savBLKData[i]->fft->data->data[jSavedBLK].re + outputSFTs->data[i].data->data[j].re;
             params->BLKData[i]->fft->data->data[j].im = savBLKData[i]->fft->data->data[jSavedBLK].im + outputSFTs->data[i].data->data[j].im;

          }
        }
    /* next is else for if ( (params->testFlag & 4) > 0 ) else... */
    } else {
        #ifdef PRINT_ONEMONTECARLO_OUTPUTSFT
          if (params->whichMCSUM == 0) {
            /* To compare with LALFastGeneratePulsarSFTs, which uses ComputeSky, use first time stamp as reference time */
            GPSin.gpsSeconds = params->timeStamps[0].gpsSeconds;
            GPSin.gpsNanoSeconds = params->timeStamps[0].gpsNanoSeconds;
            /* LALConvertGPS2SSB(status->statusPtr,&(pPulsarSignalParams->pulsar.tRef), GPSin, pPulsarSignalParams); CHECKSTATUSPTR (status); */
            /* 08/31/05 gam */
            pPulsarSignalParams->pulsar.tRef.gpsSeconds = GPSin.gpsSeconds;
            pPulsarSignalParams->pulsar.tRef.gpsNanoSeconds = GPSin.gpsNanoSeconds;
          }
        #endif

        signal = NULL; /* Call LALGeneratePulsarSignal to generate signal */
        LALGeneratePulsarSignal(status->statusPtr, &signal, pPulsarSignalParams);
        CHECKSTATUSPTR (status);

        #ifdef DEBUG_MONTECARLOTIMEDOMAIN_DATA
          fprintf(stdout,"signal->deltaT = %23.10e \n",signal->deltaT);
          fprintf(stdout,"signal->epoch.gpsSeconds = %i \n",signal->epoch.gpsSeconds);
          fprintf(stdout,"signal->epoch.gpsNanoSeconds = %i \n",signal->epoch.gpsNanoSeconds);
          fprintf(stdout,"pPulsarSignalParams->duration*pPulsarSignalParams->samplingRate, signal->data->length = %g %i\n",pPulsarSignalParams->duration*pPulsarSignalParams->samplingRate,signal->data->length);
          fflush(stdout);
        #endif
        #ifdef PRINT_MONTECARLOTIMEDOMAIN_DATA
          for(i=0;i<signal->data->length;i++)  {
            fprintf(stdout,"signal->data->data[%i] = %g\n",i,signal->data->data[i]);
            fflush(stdout);
          }
        #endif

        outputSFTs = NULL; /* Call LALSignalToSFTs to generate outputSFTs with injection data */
        LALSignalToSFTs(status->statusPtr, &outputSFTs, signal, pSFTParams);
        CHECKSTATUSPTR (status);

        #ifdef PRINT_ONEMONTECARLO_OUTPUTSFT
           if (params->whichMCSUM == 0) {
              REAL4 realHetPhaseCorr;
              REAL4 imagHetPhaseCorr;
              REAL8 deltaTHeterodyne;
              i=0; /* index of which outputSFT to output */
              fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",iFreq,h_0);
              fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",iFreq,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
              fprintf(stdout,"iFreq = %i, inject psi = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.psi);
              fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.phi0);
              fprintf(stdout,"iFreq = %i, search f0 = %23.10e, inject f0 = %23.10e \n",iFreq,f0SUM + iFreq*params->dfSUM,pPulsarSignalParams->pulsar.f0);
              fprintf(stdout,"outputSFTs->data[%i].data->length = %i \n",i,outputSFTs->data[i].data->length);
              /* deltaTHeterodyne = the time from the refererence time to the start of an SFT */
              GPSin.gpsSeconds = params->timeStamps[i].gpsSeconds;
              GPSin.gpsNanoSeconds = params->timeStamps[i].gpsNanoSeconds;
              TRY (LALDeltaFloatGPS (status->statusPtr, &deltaTHeterodyne, &GPSin, &(pPulsarSignalParams->pulsar.tRef)), status);
              /* deltaTHeterodyne = 0.0; */ /* UNCOMMENT TO TURN OFF CORRECTION OF PHASE */
              fprintf(stdout,"t_h = %23.10e\n",deltaTHeterodyne);
              fflush(stdout);
              realHetPhaseCorr = cos(LAL_TWOPI*pPulsarSignalParams->fHeterodyne*deltaTHeterodyne);
              imagHetPhaseCorr = sin(LAL_TWOPI*pPulsarSignalParams->fHeterodyne*deltaTHeterodyne);
              renorm = ((REAL4)nSamples)/((REAL4)(outputSFTs->data[i].data->length - 1));
              for(j=0;j<params->nBinsPerBLK;j++) {
                /* fprintf(stdout,"%i %g %g \n",j,renorm*outputSFTs->data[i].data->data[j].re*realHetPhaseCorr - renorm*outputSFTs->data[i].data->data[j].im*imagHetPhaseCorr,renorm*outputSFTs->data[i].data->data[j].re*imagHetPhaseCorr + renorm*outputSFTs->data[i].data->data[j].im*realHetPhaseCorr); */
                fprintf(stdout,"%i %g \n",j,renorm*renorm*outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + renorm*renorm*outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im);
                fflush(stdout);
              }
           }
        #endif

        #ifdef PRINT_MAXPOWERANDBINEACHSFT
         {
           INT4 jMaxPwr;
           REAL4 maxPwr;
           REAL4 pwr;
           for(i=0;i<params->numBLKs;i++) {
             maxPwr = 0.0;
             jMaxPwr = -1;
             renorm = ((REAL4)nSamples)/((REAL4)(outputSFTs->data[i].data->length - 1));
             for(j=0;j<params->nBinsPerBLK;j++) {
               pwr = outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im;
               pwr = renorm*renorm*pwr;
               if (pwr > maxPwr) {
                  maxPwr = pwr;
                  jMaxPwr = j;
               }
             }
             fprintf(stdout,"Max power for SFT %i is in bin %i = %g \n",i,jMaxPwr,maxPwr);
             fflush(stdout);
           }
         }
        #endif
        
        /* Add outputSFTs with injected signal to input noise SFTs; renorm is needed. */
        /* 05/21/04 gam; Normalize the SFTs from LALSignalToSFTs as per the normalization in makefakedata_v2.c and lalapps/src/pulsar/make_sfts.c. */
        renorm = ((REAL4)nSamples)/((REAL4)(outputSFTs->data[0].data->length - 1)); /* 05/21/04 gam; should be the same for all outputSFTs; note minus 1 is needed */
        #ifdef DEBUG_MONTECARLOSFT_DATA  
          fprintf(stdout,"Mutiplying outputSFTs->data[%i] with renorm = %g \n",i,renorm);
          fflush(stdout);
        #endif
        for(i=0;i<params->numBLKs;i++) {
          params->BLKData[i]->fft->f0 = params->f0BLK; /* update with each injection */
          for(j=0;j<params->nBinsPerBLK;j++) {
             jSavedBLK = j + jStart; /* 09/12/05 gam */
             #ifdef PRINTCOMPARISON_INPUTVSMONTECARLOSFT_DATA
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].re, outputSFTs->data[%i].data->data[%i].re = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[jSavedBLK].re,renorm*outputSFTs->data[i].data->data[j].re);
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].im, outputSFTs->data[%i].data->data[%i].im = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[jSavedBLK].im,renorm*outputSFTs->data[i].data->data[j].im);
               fflush(stdout);
             #endif
             params->BLKData[i]->fft->data->data[j].re = savBLKData[i]->fft->data->data[jSavedBLK].re + renorm*outputSFTs->data[i].data->data[j].re;
             params->BLKData[i]->fft->data->data[j].im = savBLKData[i]->fft->data->data[jSavedBLK].im + renorm*outputSFTs->data[i].data->data[j].im;
          }
        }
    } /* END if ( (params->testFlag & 4) > 0 ) else */

    #ifdef DEBUG_MONTECARLOSFT_DATA  
        fprintf(stdout,"params->numBLKs, outputSFTs->length = %i, %i \n",params->numBLKs,outputSFTs->length);
        fflush(stdout);
        for(i=0;i<outputSFTs->length;i++) {
          fprintf(stdout,"params->f0BLK, outputSFTs->data[%i].f0 = %g, %g\n",i,params->f0BLK,outputSFTs->data[i].f0);
          fprintf(stdout,"params->tBLK, 1.0/outputSFTs->data[%i].deltaF = %g, %g\n",i,params->tBLK,1.0/outputSFTs->data[i].deltaF);
          fprintf(stdout,"params->timeStamps[%i].gpsSeconds, outputSFTs->data[%i].epoch.gpsSeconds = %i, %i\n",i,i,params->timeStamps[i].gpsSeconds,outputSFTs->data[i].epoch.gpsSeconds);
          fprintf(stdout,"params->timeStamps[%i].gpsNanoSeconds, outputSFTs->data[%i].epoch.gpsNanoSeconds = %i, %i\n",i,i,params->timeStamps[i].gpsNanoSeconds,outputSFTs->data[i].epoch.gpsNanoSeconds);
          fprintf(stdout,"params->nBinsPerBLK, outputSFTs->data[%i].data->length = %i, %i\n",i,params->nBinsPerBLK,outputSFTs->data[i].data->length);
          fflush(stdout);
        }
    #endif
    #ifdef PRINT_MONTECARLOSFT_DATA
        for(i=0;i<outputSFTs->length;i++) {
          for(j=0;j<outputSFTs->data[i].data->length;j++) {
            fprintf(stdout,"outputSFTs->data[%i].data->data[%i].re, outputSFTs->data[%i].data->data[%i].im = %g, %g\n",i,j,i,j,outputSFTs->data[i].data->data[j].re,outputSFTs->data[i].data->data[j].im);
            fflush(stdout);
          }
        }
    #endif

    /* 09/01/05 gam; this loop will run at least once; if params->numMCRescalings > 0 rescale the SFTs and rerun the injection params->numMCRescalings times */
    for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) {

      /* if (kSUM > 0 ) */ /* 09/01/05 gam */
      if ( (kSUM > 0)  || (iRescale > 0) ) {
       params->startSUMs = 0;  /* Indicate that SUMs already started */
      }
      
      if ( (iRescale % 2) == 1 ) {
        h_0Rescaled[iRescale] = h_0 + ((REAL8)((iRescale + 1)/2))*params->rescaleMCFraction*h_0;
      } else {
        h_0Rescaled[iRescale] = h_0 - ((REAL8)(iRescale/2))*params->rescaleMCFraction*h_0;
      }      
      params->threshold4 = h_0Rescaled[iRescale];
      
      if (iRescale > 0) {
         /* Rescale the SFTs using the ratio h_0Rescaled[iRescale]/h_0 */
         REAL4 rescaleFactor = ((REAL4)(h_0Rescaled[iRescale]/h_0));
         if ( !( (params->testFlag & 4) > 0 ) ) {
           /* Not using LALFastGeneratePulsarSFTs! Need to include renorm factor in RescaleFactor */
           renorm = ((REAL4)nSamples)/((REAL4)(outputSFTs->data[0].data->length - 1)); /* 05/21/04 gam; should be the same for all outputSFTs; note minus 1 is needed */
           rescaleFactor = rescaleFactor*renorm;
         } /* END if ( (params->testFlag & 4) > 0 ) else ... */
         /* Add rescaled outputSFTs with injected signal to input noise SFTs */
         for(i=0;i<params->numBLKs;i++) {
           for(j=0;j<params->nBinsPerBLK;j++) {
                jSavedBLK = j + jStart; /* 09/12/05 gam */
                params->BLKData[i]->fft->data->data[j].re = savBLKData[i]->fft->data->data[jSavedBLK].re + rescaleFactor*outputSFTs->data[i].data->data[j].re;
                params->BLKData[i]->fft->data->data[j].im = savBLKData[i]->fft->data->data[jSavedBLK].im + rescaleFactor*outputSFTs->data[i].data->data[j].im;
           }
         }
      } /* END if (iRescale > 0) */

      /**************************************************************/
      /* Call StackSlideApplySearch to analyze this MC injection!!! */
      /**************************************************************/
      params->maxPower = 0.0; /* 05/24/05 gam; need to reinitialize */
      StackSlideApplySearch(status->statusPtr,params);
      CHECKSTATUSPTR (status);

      /* 05/24/05 gam */
      if (reportMCResults) {
         nTrials[iRescale] += 1.0;  /* count total number of MC trials run */
         /* priorLoudestEvent is the Loudest event from prior jobs in the pipeline */
         if (params->maxPower >= priorLoudestEvent) {
            nAboveLE[iRescale] += 1.0;  /* count number of MC trials that result in power >= to the loudest event */
         }
      }

    } /* END for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) */ /* 09/01/05 gam */

    /* 09/01/06 gam; moved after StackSlideApplySearch and outside iRescale loop */
    /* 07/14/04 gam; check if using LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
    if ( !((params->testFlag & 4) > 0) ) {
        LALDestroySFTVector(status->statusPtr, &outputSFTs);
        CHECKSTATUSPTR (status); /* 06/01/04 gam */
      
        LALFree(signal->data->data);
        LALFree(signal->data);
        LALFree(signal);
    }

  } /* END for(kSUM=0;kSUM<numSUMsTotal;kSUM++) */
  /**********************************************/
  /*                                            */
  /* END SECTION: MONTE CARLO INJECTIONS LOOP   */
  /*                                            */
  /**********************************************/

  /* 05/24/05 gam */  
  if (reportMCResults) {

     /* 09/01/05 gam */
     breakOutOfLoop = 0; /* 09/01/05 gam; if we converge, then exit interated MC loop */     
     iRescaleBest = 0;   /* keep track of which rescaling came closest to desired confidence. */             
     lastMCErr = 1;      /* keep track of last absolute error in confidence. */
     for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) {
       /* Compute confidence; note that h_0 == params->threshold4 == injection amplitude */
       arrayULs[countMC] = h_0Rescaled[iRescale];
       if (nTrials[iRescale] > 0.0) {
         arrayConfs[countMC] = nAboveLE[iRescale]/nTrials[iRescale];
         thisMCErr = fabs(priorConfidence - arrayConfs[countMC]);
         /* Check, are we closer? */
         if ( thisMCErr < lastMCErr ) {
            iRescaleBest = iRescale;
            lastMCErr = thisMCErr;
         }
         /* Check, are we done? */
         /* if ( fabs(priorConfidence - arrayConfs[countMC])/priorConfidence <= maxMCfracErr ) */ /* 09/06/05 gam */
         if ( thisMCErr <= maxMCErr ) {
            arrayConverged[countMC] = 1;  /* converged on desired confidence */
            breakOutOfLoop = 1; /* We are within maxMCErr; break out of the iMC loop */
         } else {
            arrayConverged[countMC] = 0; /* not converged on desired confidence */
         }       
       } else {
         arrayConfs[countMC] = -1;
         arrayConverged[countMC] = 0; /* not converged on desired confidence */
         breakOutOfLoop = 1; /* MC not producing any trials; should not happen but stop if it does */
       }
       countMC++; /* Initialzed to 0 above; count how many MC simulations have been completed */
     } /* END for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) */

     if (breakOutOfLoop) {
         break; /* MC converged or not producing any trials so break out of iMC loop. */
     }

     if (params->numMCRescalings > 0) {
       StackSlideULLeastSquaresLinFit(&(arrayULs[countMC]),arrayULs,arrayConfs,priorConfidence,(countMC - onePlusNumMCRescalings), onePlusNumMCRescalings);
       arrayConfs[countMC] = -1; /* We hope this is priorConfidence, but we do not know yet what the confidence is for this interpolated value. */
       arrayConverged[countMC] = -1; /* We do not know yet whether this interpolated UL converges on the desired confidence. */
       countMC++; /* Initialzed to 0 above; count how many MC simulations have been completed */
     } 

     /* 09/01/05 gam; when iterating the entire MC use best values from each iteration */
     arrayIteratedULs[iMC] = h_0Rescaled[iRescaleBest];
     arrayIteratedConfs[iMC] = nAboveLE[iRescaleBest]/nTrials[iRescaleBest];
          
     /* if (maxMC > 1) then interate entire MC to converge on desired confidence */
     if (maxMC > 1) {
        /* Compute new guess for h_0 */
        if (params->numMCRescalings > 0) {
            /* 09/01/05 gam; if just finished first iteration; next try is interpolated value from rescaled MC */
            if (arrayULs[countMC-1] > 0.0) {
                 h_0 = arrayULs[countMC-1];  /* Note that the minus 1 is needed because we did countMC++ above */
            } else {
                 h_0 = arrayIteratedULs[iMC]; /* The least squares fit failed for some reason; go with the best value from above */
            }
        } else {
          if ( iMC > 0 ) {
            /* linearly interpolate to next guess */
            /* if ( (arrayIteratedConfs[iMC] - arrayIteratedConfs[iMC-1]) != 0.0 ) {
              h_0 = arrayIteratedULs[iMC] + (priorConfidence - arrayIteratedConfs[iMC])*(arrayIteratedULs[iMC]
                   - arrayIteratedULs[iMC-1])/(arrayIteratedConfs[iMC] - arrayIteratedConfs[iMC-1]);
            } */  /* 09/01/05 gam; find the least square linear fit to the previous iteration to get next guess */
            StackSlideULLeastSquaresLinFit(&h_0,arrayIteratedULs,arrayIteratedConfs,priorConfidence,0,iMC+1);
            if (h_0 <= 0.0) {
               h_0 = h_0Rescaled[0]; /* The least squares fit failed for some reason; will force new guess just below. */
            }
          } else {
            h_0 = h_0Rescaled[0]; /* Will force new guess just below. */
          }
        }
        if ( h_0Rescaled[0] == h_0 ) {
           /* Make sure guess for new iteration differs from previous guess. Note that */
           /* if nTrials[0] == 0 then we would have broken out of the iMC loop above. */
           if ( (priorConfidence - (nAboveLE[0]/nTrials[0])) > 0.0 ) {
              h_0 = h_0Rescaled[0] + priorUncertainty; /* need to make h_0 larger */
           } else {
              h_0 = h_0Rescaled[0] - priorUncertainty; /* need to make h_0 smaller */
           }
        }
        if (h_0 <= 0.0) {
           h_0 = priorUncertainty; /* prevent making h_0 too small */
        }
        params->threshold4 = h_0; /* set params->threshold4 == h_0 to new guess */
     }
  } /* END if (reportMCResults) */

 } /* END for(iMC=0;iMC<maxMC;iMC++) */
 /***********************************************************/
 /*                                                         */
 /* END SECTION: LOOP OVER ENTIRE MONTE CARLO SIMULATIONS   */
 /*                                                         */
 /***********************************************************/

  /* 05/24/05 gam; always finishPeriodicTable table here */
  if (params->outputEventFlag > 0) {
      fprintf( params->xmlStream->fp, STACKSLIDE_XML_TABLE_FOOTER );
      params->xmlStream->table = no_table;
  }
  
  /* 05/24/05 gam */  
  if (reportMCResults) {
    if (params->outputEventFlag > 0) {
       /* write to the StackSlide Monte Carlo Results Table. */
       StackSlideMonteCarloResultsTable *stksldMonteCarloResultsTable;
       stksldMonteCarloResultsTable = (StackSlideMonteCarloResultsTable *) LALCalloc( 1, sizeof(StackSlideMonteCarloResultsTable) );
            
       /* print the table header */
       fprintf( params->xmlStream->fp, LIGOLW_XML_LOCAL_SEARCHRESULTS_STACKSLIDEMONTECARLO );

       /* print out the rows */
       for(iMC=0;iMC<countMC;iMC++) {
          stksldMonteCarloResultsTable->loudest_event = priorLoudestEvent; /* from prior jobs in the pipeline */
          stksldMonteCarloResultsTable->start_freq    = f0SUM;             /* saved value MC actually ran on  */
          stksldMonteCarloResultsTable->band          = bandSUM;           /* saved value MC actually ran on  */
          stksldMonteCarloResultsTable->upper_limit   = arrayULs[iMC]/params->blkRescaleFactor; /* 10/20/05 gam; need to remove blkRescaleFactor which has a default value of 1.0 */
          stksldMonteCarloResultsTable->confidence    = arrayConfs[iMC];
          stksldMonteCarloResultsTable->converged     = arrayConverged[iMC];
          stksldMonteCarloResultsTable->next          = NULL;
          if (iMC > 0) {
                 fprintf( params->xmlStream->fp, ",\n" ); /* end previous row with a comma */
          }
          fprintf( params->xmlStream->fp, LOCAL_SEARCHRESULTS_STACKSLIDEMONTECARLO_ROW,
            stksldMonteCarloResultsTable->loudest_event,
            stksldMonteCarloResultsTable->start_freq,
            stksldMonteCarloResultsTable->band,
            stksldMonteCarloResultsTable->upper_limit,
            stksldMonteCarloResultsTable->confidence,
            stksldMonteCarloResultsTable->converged
          );
       } /* END for(iMC=0;iMC<countMC;iMC++) */

       /* End the table: */
       fprintf( params->xmlStream->fp, STACKSLIDE_XML_TABLE_FOOTER );
       params->xmlStream->table = no_table;
       LALFree(stksldMonteCarloResultsTable);
    } /* END if (params->outputEventFlag > 0) */
    #ifdef INCLUDE_PRINT_REPORTMCRESULTS_CODE
      for(iMC=0;iMC<countMC;iMC++) {
        if ((params->debugOptionFlag & 32) > 0 ) {
              fprintf(stdout,"MC Trial #%i UL = %20.10e, Conf = %20.10f, Converged = %i\n",iMC,arrayULs[iMC]/params->blkRescaleFactor,arrayConfs[iMC],arrayConverged[iMC]);
              fflush(stdout);
        }
      }
    #endif    
    LALFree(arrayULs);
    LALFree(arrayConfs);
    LALFree(arrayConverged);
    LALFree(arrayIteratedULs);
    LALFree(arrayIteratedConfs);    
    LALFree(nAboveLE);   /* 09/01/05 gam */
    LALFree(nTrials);    /* 09/01/05 gam */
  } /* END if (reportMCResults) */
  LALFree(h_0Rescaled); /* 09/01/05 gam */
      
  LALDestroyRandomParams(status->statusPtr, &randPar); /* 05/28/04 gam */
  CHECKSTATUSPTR (status);
  
  LALFree(pSFTParams);
      
  if (params->numSpinDown > 0) {
    LALDDestroyVector(status->statusPtr, &(pPulsarSignalParams->pulsar.spindown));
    CHECKSTATUSPTR (status);
  }
  LALFree(pPulsarSignalParams);

  /* 07/14/04 gam; deallocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  if ( (params->testFlag & 4) > 0 ) {
    #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
      fprintf(stdout,"deallocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs \n");
      fflush(stdout);
    #endif
    LALDestroySFTVector(status->statusPtr, &outputSFTs);
    CHECKSTATUSPTR (status);
    LALFree(pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi);
    LALFree(pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi);
    LALFree(pSkyConstAndZeroPsiAMResponse->skyConst);
    LALFree(pSkyConstAndZeroPsiAMResponse);
    LALFree(pSFTandSignalParams->trigArg);
    LALFree(pSFTandSignalParams->sinVal);
    LALFree(pSFTandSignalParams->cosVal);
    LALFree(pSFTandSignalParams);
  }
          
  /* 05/21/04 gam; free memory of saved data */
  for (i=0;i<params->numBLKs;i++) {
      LALFree(savBLKData[i]->fft->data->data);
      LALFree(savBLKData[i]->fft->data);
      LALFree(savBLKData[i]->fft);
      LALFree(savBLKData[i]);
  }
  LALFree(savBLKData);
  
  /* deallocate memory for the savSkyPosData structure */
  for(i=0;i<numSkyPosTotal;i++)
  {
      LALFree(savSkyPosData[i]);
  }
  LALFree(savSkyPosData);
  
  if (params->numSpinDown > 0) {
    /* deallocate memory for the savFreqDerivData structure */
    for(i=0;i<numFreqDerivTotal;i++)
    {
        LALFree(savFreqDerivData[i]);
    }
    LALFree(savFreqDerivData);
  }
  
  LALFree(savSumBinMask);          /* 07/15/05 gam */
  LALFree(arrayAvailableFreqBins); /* 07/15/05 gam */
  
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);
}
/*************************************************************/
/*                                                           */
/* END FUNCTION: RunStackSlideIsolatedMonteCarloSimulation   */
/* Injects fake signals and runs Monte Carlo simulation.     */
/*                                                           */
/*************************************************************/

/* 05/24/05 gam; Function that reads results from previous jobs in the pipeline */
void getStackSlidePriorResults(LALStatus *status,
                               REAL4 *priorLoudestEvent,
                               REAL8 *priorStartFreq,
                               REAL8 *priorBand,
                               REAL8 *priorConfidence,
                               REAL8 *priorUL,
                               REAL8 *priorUncertainty,
                               CHAR  *priorResultsFile)
{
  INT4 i;
  INT4 numWordsToIgnore = 43;
  CHAR buffer[256];
  FILE *fpPrior = NULL;
  
  INITSTATUS( status, "getStackSlidePriorResults", STACKSLIDEISOLATEDC );
  ATTATCHSTATUSPTR (status);
     
  fpPrior = fopen( priorResultsFile, "r");
  if (fpPrior == NULL) {
     ABORT( status, STACKSLIDEISOLATEDH_EBADRESULTSFILE, STACKSLIDEISOLATEDH_MSGEBADRESULTSFILE);
  }

  /* skip first numWordsToIgnore words */
  for(i=0;i<numWordsToIgnore;i++){
     fscanf(fpPrior,"%s",buffer);
  }
  /* data was written using matlab format string '%20.6f %20.10f %20.10f %8.2f %20.10e %20.10e\n' */
  fscanf(fpPrior,"%f%lf%lf%lf%le%le\n",priorLoudestEvent,priorStartFreq,priorBand,priorConfidence,priorUL,priorUncertainty);
  fclose(fpPrior);

  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);
}

/* 07/15/2005 gam; Change RunStackSlideIsolatedMonteCarloSimulation */
void StackSlideGetUniformDeviate(LALStatus *status, REAL8 *returnVal, REAL8 startVal, REAL8 range, RandomParams *randPar)
{  
  REAL4 randval;
  
  INITSTATUS( status, "StackSlideGetUniformDeviate", STACKSLIDEISOLATEDC );
  ATTATCHSTATUSPTR (status);

  LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status);
  *returnVal = startVal + ((REAL8)randval)*range; /* return value in the interval (startVal, startVal+range) */
  
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);
}

/* 09/01/05 gam */
void StackSlideULLeastSquaresLinFit(REAL8 *interpolatedUL, const REAL8 *arrayULs, const REAL8 *arrayConfs, REAL8 desiredConf, INT4 startIndex, INT4 numVals)
{
     /* Use standard least square linear regression fit of data :         
          A*slope + B*intercept = E,
          C*slope + D*intercept = F,
        where A,B,C,D come from minimizing sum of square deviations,
          sum_i (slope*arrayULs[i] + intercept - arrayConfs[i])^2 
     */

     INT4 i;
     REAL8 A,B,C,D,E,F,determ,slope,intercept,startIndexPlusNumVals;
     startIndexPlusNumVals = startIndex+numVals;
     if (numVals < 2) {
        *interpolatedUL = -1;
        return;
     }

     A = 0.0;
     B = 0.0;
     C = 0.0;
     D = numVals;
     E = 0;
     F = 0;     
     for(i=startIndex;i<startIndexPlusNumVals;i++) {
        A += arrayULs[i]*arrayULs[i];
        B += arrayULs[i];
        E += arrayULs[i]*arrayConfs[i];
        F += arrayConfs[i];
     }
     C = B;
     determ = A*D-B*C;
     if (determ !=0) {
       slope = (D*E-B*F)/determ;
       intercept = (A*F-C*E)/determ;
       if (slope !=0) {
          *interpolatedUL = (desiredConf - intercept)/slope;
          return;
       } else {
          *interpolatedUL = -1;
          return;
       }
     } else {
       *interpolatedUL = -1;
       return;
     }
} /* END INT4 StackSlideULLeastSquaresLinFit */
