/*********************************** <lalVerbatim file="StackSlideIsolatedCV">
Authors: Mendell, Gregory
$Id$
************************************ </lalVerbatim> */

/* REVISIONS: */

/*********************************************/
/*                                           */
/* START SECTION: define preprocessor flags  */
/*                                           */
/*********************************************/
#define INCLUDE_PRINT_PEAKS_TABLE_CODE
/* #define DEBUG_STACKSLIDE_ISOLATED */
/* #define DEBUG_SUM_TEMPLATEPARAMS */
/* #define PRINT_STACKSLIDE_BINOFFSETS */
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
    REAL4                            *maxPower,
    INT4                             *totalEventCount,
    SnglStackSlidePeriodicTable      *loudestPeaksArray,
    LALFindStackSlidePeakOutputs     *pLALFindStackSlidePeakOutputs,
    LALFindStackSlidePeakParams      *pLALFindStackSlidePeakParams,
    LALUpdateLoudestStackSlideParams *pLALUpdateLoudestStackSlideParams,
    LALDetector                      *cachedDetector,
    StackSlideParams                 *stksldParams,
    StackSlideSearchParams           *params
)
{

  INT4 i, j, k, kSUM;
  SnglStackSlidePeriodicTable *peakSav;
  INT4 isav = -1;      /* 10/28/04 gam; saved index value */ 
  
  #ifdef DEBUG_STACKSLIDE_ISOLATED
    fprintf(stdout, "\nStart StackSlideIsolated\n\n");
    fflush(stdout);
  #endif    

  INITSTATUS( status, "StackSlideIsolated", STACKSLIDEISOLATEDC );
  ATTATCHSTATUSPTR(status);
  
  for(kSUM=0;kSUM<params->numSUMsTotal;kSUM++) {
    
    /* 02/11/04 gam; set up pointer to current sky position and spindown parameters */
    i = kSUM/params->numFreqDerivIncludingNoSpinDown;   /* 01/28/04 gam; index to params->skyPosData for this SUM; */
    j = kSUM % params->numFreqDerivIncludingNoSpinDown; /* 01/28/04 gam; index to params->freqDerivData for this SUM; */    
    stksldParams->skyPosData = &(params->skyPosData[i]);
    if (params->numFreqDerivTotal != 0) {
         stksldParams->freqDerivData = &(params->freqDerivData[j]);
    }    

    /* 10/28/04 gam; weight STKs; the weights depend on F_+ and F_x when params->weightSTKsIncludingBeamPattern is True */
    if (params->weightSTKsIncludingBeamPattern) {
        if (i != isav) {
       
         /* 10/28/04 gam; get squared detector response F_+^2 or F_x^2 for one sky position, one polarization angle, for midpoints of a timeStamps */
         GetDetResponseTStampMidPts(status->statusPtr, params->detResponseTStampMidPts, params->timeStamps, params->numSTKs, params->tSTK,
              cachedDetector, params->skyPosData[i], params->orientationAngle, COORDINATESYSTEM_EQUATORIAL, params->plusOrCross);
         
         /* 10/28/04 gam; apply powerFlux style weights including detector beam pattern response */
         WeightSTKsIncludingBeamPattern(params->STKData,
               params->savSTKData,
               params->inverseSquareMedians,
               params->sumInverseSquareMedians,
               params->detResponseTStampMidPts,
               params->numSTKs, params->nBinsPerSTK, params->tSTK);            
        }
        isav = i;
    }
    
     /* 11/08/03 gam; Add in first version of StackSlide written by Mike Landry */     
    StackSlideOld(status->statusPtr,params->SUMData,params->STKData,stksldParams); /*Feb 05/14 vir*/
    CHECKSTATUSPTR (status);
    
    /* 02/25/05 gam; utility for printing one SUM to a file */
    if (params->outputSUMFlag > 0) {
      printOneStackSlideSUM(params->SUMData[0],params->outputSUMFlag,params->outputFile,kSUM,params->f0SUM,params->dfSUM,params->nBinsPerSUM,params->numSUMsTotal);
    }

    #ifdef DEBUG_SUM_TEMPLATEPARAMS
        fprintf(stdout,"\nTemplate parameters for SUM #%i:\n",kSUM);	
        fprintf(stdout,"RA = %18.10f\n",params->skyPosData[i][0]);
        fprintf(stdout,"DEC = %18.10f\n",params->skyPosData[i][1]);
        for(k=0;k<params->numSpinDown;k++)
        {
           fprintf(stdout,"FREQDERIV%i = %18.10f\n",k+1,params->freqDerivData[j][k]);
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
         pLALFindStackSlidePeakParams->iSky = i;
         pLALFindStackSlidePeakParams->iDeriv = j;
         pLALFindStackSlidePeakOutputs->pntrToPeaksPntr = &(params->peaks);  /* 03/02/04 gam; added this and changed next line */
         LALFindStackSlidePeaks(pLALFindStackSlidePeakOutputs, params->SUMData[0], pLALFindStackSlidePeakParams);
                 
         params->peaks = peakSav->next;  /* Go back to the beginning of the list */
         LALFree(peakSav);               /* No longer need this place holder */

         while (params->peaks) {
           peakSav = params->peaks;
           (*totalEventCount)++;
           if (pLALFindStackSlidePeakParams->updateMeanStdDev) {
              /* 03/02/04 gam; update pw_mean_thissum, pw_stddev_thissum, pwr_snr using pwMeanWithoutPeaks and pwStdDevWithoutPeaks */              
              if ( (*(pLALFindStackSlidePeakOutputs->binsWithoutPeaks)) > 0) {
                /* Note that LALFindStackSlidePeaks ensures that pwStdDevWithoutPeaks is never 0. */
                params->peaks->pw_mean_thissum = *(pLALFindStackSlidePeakOutputs->pwMeanWithoutPeaks);
                params->peaks->pw_stddev_thissum = *(pLALFindStackSlidePeakOutputs->pwStdDevWithoutPeaks);
                params->peaks->pwr_snr = params->peaks->power / (*(pLALFindStackSlidePeakOutputs->pwStdDevWithoutPeaks));
              }
           }
           if (params->peaks->power > *maxPower) {
              *maxPower = params->peaks->power; /* 02/09/04 gam; keep track for search summary */
           }
           
           /* 02/17/04 gam; keep an array of loudest peaks */
           if (params->outputLoudestFromPeaks) {
               LALUpdateLoudestFromPeaks(loudestPeaksArray, params->peaks, params->nBinsPerOutputEvent);
           } else {
              #ifdef INCLUDE_PRINT_PEAKS_TABLE_CODE
                if ((params->debugOptionFlag & 2) > 0 ) {
                  fprintf(stdout,"%11i %11i %10.6f %10.6f %14.6e %11.6f %10.3f %8.3f %12i \n",*totalEventCount,params->peaks->sum_no,params->peaks->sky_ra,params->peaks->sky_dec,params->peaks->fderiv_1,params->peaks->frequency,params->peaks->power,params->peaks->pwr_snr,params->peaks->width_bins);
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
           pLALUpdateLoudestStackSlideParams->iSUM = kSUM;  /* 02/17/04 gam; see top of for loop for kSUM, i, and j */
         }
         pLALUpdateLoudestStackSlideParams->iSky = i;
         pLALUpdateLoudestStackSlideParams->iDeriv = j;              
         LALUpdateLoudestFromSUMs(loudestPeaksArray, params->SUMData[0], pLALUpdateLoudestStackSlideParams);
    } /* end else if (params->outputLoudestFromSUMs) */
      
  } /* 02/11/04 gam; new end for(kSUM=0;kSUM<params->numSUMsTotal;kSUM++) */
  
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
       if (loudestPeaksArray[k].power > *maxPower) {
              *maxPower = loudestPeaksArray[k].power; /* 02/20/04 gam; keep track for search summary */
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
  
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);

} /* END StackSlideIsolated */
/*********************************************************************************/
/*              END FUNCTION: StackSlideIsolated                                 */
/*********************************************************************************/

/* <lalVerbatim file="StackSlideCP"> */
void StackSlideOld(	LALStatus *status, 
			REAL4FrequencySeries **SUMData, 
			REAL4FrequencySeries **STKData, 
			StackSlideParams *params)
{  /* </lalVerbatim> */                        
    
    /* Declare and variables needed for making SUMs */
    INT4 iSky = 0;       /* index gives which Sky Position */
    INT4 iFreqDeriv = 0; /* index give which Freq Derivs  */
    INT4 kSUM = 0;       /* index gives which SUM */
    INT4 iSUM = 0;       /* index gives which SUM bin */
    
    INT4 numLoopOnSpindown; /* 1 if numSpindown is 0, else = params->numFreqDerivTotal */
    
    INT4 i,k,m;
    INT4 binoffset = 0;

    StackSlideSkyParams *csParams;  /* Container for StackSlideComputeSky */
    
    TdotsAndDeltaTs *pTdotsAndDeltaTs;
    EarthState earth;
    EmissionTime emit;

    REAL8 f_t;

    REAL4 invNumSTKs = 1.0/((REAL4)params->numSTKs);  /* 12/03/04 gam */
    
    INT4 iMinSTK = floor((params->f0SUM-params->f0STK)*params->tSTK + 0.5); /* Index of mimimum frequency to include when making SUMs from STKs */
    INT4 iMaxSTK = iMinSTK + params->nBinsPerSUM - 1;                       /* Index of maximum frequency to include when making SUMs from STKs */
    
    REAL8 refFreq = params->f0SUM + ((REAL8)(params->nBinsPerSUM/2))*params->dfSUM; /* 12/06/04 gam */
    
    INITSTATUS (status, "StackSlide", STACKSLIDEISOLATEDC);
    ATTATCHSTATUSPTR(status);

    /* Allocate space and set quantities for call to StackSlideComputeSky() */
    csParams=(StackSlideSkyParams *)LALMalloc(sizeof(StackSlideSkyParams));
    csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
    pTdotsAndDeltaTs=(TdotsAndDeltaTs *)LALMalloc(sizeof(TdotsAndDeltaTs));
    pTdotsAndDeltaTs->vecTDots  = (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSTKs);

    if (params->numSpinDown>0) {
       pTdotsAndDeltaTs->vecDeltaTs  = (REAL8 **)LALMalloc(sizeof(REAL8 *)*params->numSTKs);          
       for(i=0;i<params->numSTKs;i++) {
           pTdotsAndDeltaTs->vecDeltaTs[i]  = (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSpinDown);          
       }
    }

    /* if no_spindowns > 0, only pTdotsAndDeltaTs->vecDeltaTs= (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSTKs); 
       i.e. if no spindowns then you don't need to allocate the second of the two lines above */
    /* 02/02/04 gam; moved next 5 lines from inside loops: */
    if (params->numFreqDerivTotal==0) {
        numLoopOnSpindown = 1;  /* Even if no spindown, need to execute iFreqDeriv loop once to handle case of zero spindown */
    } else {
        numLoopOnSpindown = params->numFreqDerivTotal;
    }

    for(iSky=0;iSky<params->numSkyPosTotal;iSky++) {
        /* call ComputeSky for every sky position */

        csParams->skyPos[0]=params->skyPosData[iSky][0];
        csParams->skyPos[1]=params->skyPosData[iSky][1];

        params->baryinput->alpha=params->skyPosData[iSky][0];
        params->baryinput->delta=params->skyPosData[iSky][1];

        csParams->tGPS=params->timeStamps;
        csParams->gpsStartTimeSec = params->gpsStartTimeSec; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
        csParams->gpsStartTimeNan = params->gpsStartTimeNan; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
        csParams->spinDwnOrder=params->numSpinDown;
        csParams->mObsSFT=params->numSTKs;
        csParams->tSFT=params->tSTK;
        csParams->edat=params->edat;
        csParams->emit=&emit;
        csParams->earth=&earth;
        csParams->baryinput=params->baryinput;

        for(i=0;i<params->numSTKs;i++) {
          pTdotsAndDeltaTs->vecTDots[i] = 0.0; /* Initialize */
          for(m=0;m<params->numSpinDown;m++) {
              pTdotsAndDeltaTs->vecDeltaTs[i][m] = 0.0; /* Initialize */
          }
        }

        /* get Tdot and DeltaT's for this sky position */
        StackSlideComputeSky(status->statusPtr, pTdotsAndDeltaTs, csParams);
        CHECKSTATUSPTR (status);

        /* for all spindowns, loop */
        for(iFreqDeriv=0;iFreqDeriv<numLoopOnSpindown;iFreqDeriv++) {
            /* call computesky for all stacks at once */
            /* compute offset  for each stack */

			
		kSUM = iSky*numLoopOnSpindown + iFreqDeriv; /* 01/28/04 gam; need to find index to which to SUM */

            /* loop over STKs */ 
            for(k=0;k<params->numSTKs;k++) {
               /* compute frequency */
               /* f_t = params->f0STK * params->tSTK; */
               /* f_t = params->f0STK; */ /* 12/06/04 gam */
               f_t = refFreq;
               for (m=0; m<params->numSpinDown; m++) {
                   /* f_t += ( params->freqDerivData[iFreqDeriv][m] * pTdotsAndDeltaTs[2*k*(params->numSpinDown)+2*(INT4)m+1]); */
                   f_t += params->freqDerivData[iFreqDeriv][m] * pTdotsAndDeltaTs->vecDeltaTs[k][m];
               }
               f_t = f_t * pTdotsAndDeltaTs->vecTDots[k];
               /* binoffset = floor(( (f_t - params->f0STK) * params->tSTK) + 0.5 ); */ /* 12/06/04 gam */
               binoffset = floor(( (f_t - refFreq) * params->tSTK) + 0.5 );

               #ifdef PRINT_STACKSLIDE_BINOFFSETS
                  fprintf(stdout, "In StackSlide for SFT #%i binoffset = %i \n",k,binoffset);
                  fflush(stdout);   
               #endif

               /* Add the power from this STK into the SUM with the appropriate binoffset */
               for(i=iMinSTK;i<=iMaxSTK; i++) {
                iSUM = i - iMinSTK;
                if (k==0) {
                   /* Starting a new SUM: initialize */
                   SUMData[kSUM]->data->data[iSUM] = STKData[k]->data->data[i+binoffset];
                   SUMData[kSUM]->epoch.gpsSeconds = params->gpsStartTimeSec;
                   SUMData[kSUM]->epoch.gpsNanoSeconds = params->gpsStartTimeNan;
                   SUMData[kSUM]->f0=params->f0SUM;
                   SUMData[kSUM]->deltaF=params->dfSUM;
                   SUMData[kSUM]->data->length=params->nBinsPerSUM;
                } else {
                   SUMData[kSUM]->data->data[iSUM] += STKData[k]->data->data[i+binoffset];
                }
              } /* END for(i=iMinSTK;i<=iMaxSTK; i++) */
            } /* END for(k=0;k<params->numSTKs;k++) */

            /* 12/03/04 gam; added params->divideSUMsByNumSTKs */
            if (params->divideSUMsByNumSTKs) {
               /* Normalize the SUMs with params->numSTKs*/
               for(i=0;i<params->nBinsPerSUM; i++) {
                  /* SUMData[kSUM]->data->data[i] =  SUMData[kSUM]->data->data[i]/((REAL4)params->numSTKs); */
                  SUMData[kSUM]->data->data[i] =  SUMData[kSUM]->data->data[i]*invNumSTKs;  /* 12/03/04 gam; multiply by 1.0/numSTKs */
               }
            }
            
        }  /* for iFreqDeriv = 0 to numFreqDerivTotal */      
    } /* for iSky = `0 to numSkyPosTotal */

    /* Deallocate space for ComputeSky parameters */
    LALFree(csParams->skyPos);
    LALFree(csParams);

    /* Deallocate memory */

    if (params->numSpinDown>0) {
        for(i=0;i<params->numSTKs;i++) {
           LALFree(pTdotsAndDeltaTs->vecDeltaTs[i]);          
        }
        LALFree(pTdotsAndDeltaTs->vecDeltaTs);
    }
    LALFree(pTdotsAndDeltaTs->vecTDots);
    LALFree(pTdotsAndDeltaTs);
  
    CHECKSTATUSPTR (status);
    DETATCHSTATUSPTR (status);
} /* END StackSlideOld() */
