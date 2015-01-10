/*
*  Copyright (C) 2007 Gregory Mendell
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

/**
 * \author Gregory Mendell
 */

/* REVISIONS: */

/*********************************************/
/*                                           */
/* START SECTION: define preprocessor flags  */
/*                                           */
/*********************************************/
#define INCLUDE_PRINT_PEAKS_TABLE_CODE
#define INCLUDE_PRINT_REPORTMCRESULTS_CODE
/* #define DEBUG_STACKSLIDE_FSTAT_ISOLATED */
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
#include "StackSlideFStatIsolated.h"
/*********************************************/
/*                                           */
/* END SECTION: include header files         */
/*                                           */
/*********************************************/

/*********************************************************************************/
/*              START FUNCTION: StackSlideFStatIsolated                               */
/*********************************************************************************/
void StackSlideFStatIsolated (
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
 INT4 iSky, iFreqDeriv, i,j,k;
 INT4 kSUM = -1;
  
 SnglStackSlidePeriodicTable *peakSav;
  
 /* add parameters for StackSlideComputeSky and StackSlide */
 StackSlideSkyParams *csParams;
 TdotsAndDeltaTs *pTdotsAndDeltaTs;
 EarthState earth;
 EmissionTime emit;
 
 /* Copied from DriveHoughFStat.c: */ 
 INT4 *mCohSft, nSFTs;                   /* number of SFTs in each stack and total number of SFTs */
 INT4 nStacks;                           /* number of stacks */
 INT4 sftlength;                         /* number of bins in each sft */ 
 SFTVector *inputSFTs=NULL;              /* vector of SFTtypes and SFTtype is COMPLEX8FrequencySeries */
 SFTVectorSequence stackSFTs;            /* sequence of sft vectors -- one for each stack */ 
 LIGOTimeGPSVector midTstack, startTsft;  
 REAL8  *tStack, tStackAvg;              /* duration of each stack; average duration */ 
 REAL8 deltaFstack;                      /* frequency resolution of Fstat calculation */ 
 REAL8FrequencySeriesVector FstatVect;   /* Fstatistic vectors for each stack */
 FstatStackParams FstatPar;
 REAL8 refTime;
 LIGOTimeGPS refTimeGPS;
    
 #ifdef DEBUG_STACKSLIDE_ISOLATED
   fprintf(stdout, "\nStart StackSlideFStatIsolated\n\n");
   fflush(stdout);
 #endif    

 INITSTATUS(status);
 ATTATCHSTATUSPTR(status);
 
 /* allocate memory for inputSFTs */
 TRY (LALCreateSFTVector (status->statusPtr, &inputSFTs, params->numBLKs, params->nBinsPerBLK), status); 
 
 /* Copy SFT data into SFTVector *inputSFTs struct */
 for(i=0;i<params->numBLKs;i++) {
      inputSFTs->data[i]->epoch = BLKData[i]->fft->epoch;
      inputSFTs->data[i]->f0 = BLKData[i]->fft->f0;
      inputSFTs->data[i]->deltaF = BLKData[i]->fft->deltaF;
      inputSFTs->data[i]->data->length = BLKData[i]->fft->data->length;
    for(j=0;j<params->nBinsPerBLK;j++) {
             inputSFTs->data[i].data->data[j].re = params->BLKData[i]->fft->data->data[j].re;
             inputSFTs->data[i].data->data[j].im = params->BLKData[i]->fft->data->data[j].im;
    }
 }

 /* normalize sfts */
 LAL_CALL( LALNormalizeSFTVect (&status, inputSFTs, (UINT4)params->nBinsPerNRM, 0), &status);

 /* Set up the reference time at the SSB */
 /* IMPORTANT, must set these to the epoch that gives T0 at SSB; this differs when running Monte Carlo */
 if ( (params->testFlag & 4) > 0 ) {
      stksldParams->gpsStartTimeSec = (UINT4)params->timeStamps[0].gpsSeconds;
      stksldParams->gpsStartTimeNan = (UINT4)params->timeStamps[0].gpsNanoSeconds;
 } else {
      stksldParams->gpsStartTimeSec = params->gpsEpochStartTimeSec;
      stksldParams->gpsStartTimeNan = params->gpsEpochStartTimeNan;
 }
 refTimeGPS.gpsSeconds = (INT4)stksldParams->gpsStartTimeSec;
 refTimeGPS.gpsNanoSeconds = (INT4)stksldParams->gpsStartTimeNan;
 LALGPStoFloat(status->statusPtr,&refTime,&refTimeGPS); CHECKSTATUSPTR(status);  

 /* set up the stacks */
 /* Note that params->numSTKsPerSUM is the number of stacks to make if there are not gaps in the data. */ 
 /* if sfts are to be split evenly among stacks */
 LAL_CALL( SetUpStacks1( &status, &stackSFTs, inputSFTs, params->numSTKsPerSUM), &status);
 /* if time is to be split evenly between stacks */
 /* LAL_CALL( SetUpStacks2( &status, &stackSFTs, inputSFTs, &startTsft, params->numSTKsPerSUM), &status); */

 /* This is the actual number of STKs */
 nStacks = stackSFTs.length; 
 params->numSTKs  = stackSFTs.length;  /* also update parameter that hold actual number of STKs */
 
 /* set up vector of stack durations */
 tStack = NULL;
 tStack = (REAL8 *)LALMalloc( nStacks * sizeof(REAL8));

 /* set up vector of number of sfts in each stack */
 mCohSft = NULL;
 mCohSft = (INT4 *)LALMalloc( nStacks * sizeof(INT4));

 /* set up vector containing mid times of stacks */    
 midTstack.length = nStacks;
 midTstack.data = (LIGOTimeGPS *)LALMalloc( nStacks * sizeof(LIGOTimeGPS));

 for (k=0; k<nStacks; k++) {
    LIGOTimeGPS tempT1, tempT2;
    INT4 tempInt;
    REAL8 tempF1, tempF2, tempMid;

    /* number of sfts in stack */
    tempInt = stackSFTs.data[k].length;
    mCohSft[k] = tempInt;

    /* duration of each stack */
    tempT1 = stackSFTs.data[k].data[0].epoch;
    tempT2 = stackSFTs.data[k].data[tempInt-1].epoch;
    LAL_CALL ( LALGPStoFloat ( &status, &tempF1, &tempT1), &status);
    LAL_CALL ( LALGPStoFloat ( &status, &tempF2, &tempT2), &status);
    tStack[k] = tempF2 + timeBase - tempF1;
    
    /* mid timestamp of each stack */
    tempMid = tempF1 + 0.5 * tStack[k];
    LAL_CALL ( LALFloatToGPS ( & status, midTstack.data + k, &tempMid), &status);
 }

 /* use stacks info to calculate search frequencies */
 /* Fstat is calculated at the frequency resolutions of the stacks
     here the stacks may be of different durations so we take the average
     This is valid if the stack durations are not very different which 
     will, hopefully, be true */
 tStackAvg = 0.0;
 for (k=0; k<nStacks; k++)
    tStackAvg += tStack[k];
 tStackAvg /= nStacks;
 deltaFstack = 1.0/tStackAvg; 
 /* update effective time baseline of STKs and the number of frequency bins per STK*/ 
 params->tEffSTK = tStackAvg;
 params->nBinsPerSTK = floor(params->bandSTK*params->tEffSTK + 0.5);

 /* set up parameters for Fstat calculation */
 FstatPar.nStacks = nStacks;
 FstatPar.tsStack = &midTstack;
 FstatPar.tStackAvg = tStackAvg;
 FstatPar.fBand = params->bandSTK;
 FstatPar.fStart = params->f0STK;
 FstatPar.nfSizeCylinder = 0;        /* Used by hough search */
 FstatPar.mCohSft = mCohSft;
 FstatPar.refTime = refTime; 
 FstatPar.SSBprecision = SSBPREC_RELATIVISTIC; 
 FstatPar.Dterms = params->Dterms;
 FstatPar.detector = *cachedDetector;
 FstatPar.edat = params->edat;
 FstatPar.ts = &startTsft;
 FstatPar.fdot = NULL; 
 if (params->numSpinDown>0) {
    LAL_CALL ( LALDCreateVector( &status, &(FstatPar.fdot), params->numSpinDown), &status);
 }

 /* set up memory for Fstat vectors */ 
 FstatVect.length = nStacks;
 FstatVect.data = NULL;
 FstatVect.data = (REAL8FrequencySeries *)LALMalloc(nStacks * sizeof(REAL8FrequencySeries));
 for (k=0; k<nStacks; k++) { 
      FstatVect.data[k].epoch = stackSFTs->data[k]->epoch;  /* WARNING: IS THIS SET CORRECTLY !!!! */
      FstatVect.data[k].deltaF = deltaFstack;
      FstatVect.data[k].f0 = params->f0STK;
      FstatVect.data[k].data = (REAL8Sequence *)LALMalloc(sizeof(REAL8Sequence));
      FstatVect.data[k].data->length = params->nBinsPerSTK;
      FstatVect.data[k].data->data = (REAL8 *)LALMalloc( params->nBinsPerSTK * sizeof(REAL8));
    }
 } /* end of Fstat memory allocation */

 /* Loop over sky positions */
 for(iSky=0;iSky<params->numSkyPosTotal;iSky++) {

   FstatPar.alpha = params->skyPosData[iSky][0];
   FstatPar.delta = params->skyPosData[iSky][1];

   /* Loop over frequency derivatives */
   for(iFreqDeriv=0;iFreqDeriv<params->numFreqDerivIncludingNoSpinDown;iFreqDeriv++) {
    
    if (params->numSpinDown>0) {
        for(k=0;k<params->numSpinDown;k++)
        {
           FstatPar.fdot->data[k] = params->freqDerivData[iFreqDeriv][k]);
        }
    }
    
    kSUM++; /* increment sum number */
    
    LAL_CALL(ComputeFstatStack( &status, &FstatVect, &stackSFTs, &FstatPar), &status);    
    
    /* If params->debugOptionFlag & 128 > 0 just creates a SUM from the STKs without sliding */
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

} /* END StackSlideFStatIsolated */
/*********************************************************************************/
/*              END FUNCTION: StackSlideFStatIsolated                                 */
/*********************************************************************************/
