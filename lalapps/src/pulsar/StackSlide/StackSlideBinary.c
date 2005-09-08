/*********************************** <lalVerbatim file="StackSlideBinaryCV">
Author: Virginia Re
$Id$
************************************ </lalVerbatim> */

/* Revisions: */
/* 04/12/05 gam; Add StackSlideSearchParams *params to StackSlideBinary. */
/* 04/12/05 gam; Remove from StackSlideParams *stksldParams, those already in StackSlideSearchParams *params */
/* 04/12/05 gam; Clean up indentation; remove obsolete code */
/*********************************************/
/*                                           */
/* START SECTION: define preprocessor flags  */
/*                                           */
/*********************************************/
 #define DEBUG_BINARY_CODE 
 #define DEBUG_LALFASTGENERATEPULSARSFTS
 #define PRINT_MAXPOWERANDBINEACHSFT
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

#include "StackSlideBinary.h"
/*********************************************/
/*                                           */
/* END SECTION: include header files         */
/*                                           */
/*********************************************/

NRCSID( STACKSLIDEBINARYC,  "$Id$");

/* 04/12/05 gam; add StackSlideSearchParams *params */
void StackSlideBinary(  LALStatus *status,
                        StackSlideParams *stksldParams,
                        StackSlideSearchParams *params
                      )
{
    INT4 i,k,iSMA;
    INT4 kSUM = -1;         /* 04/12/05 gam; counter that keeps track of which SUM is being computed */
    static INT4 iSky = 0;   /* index gives which Sky Position */
    INT8 iSkyCoh=0;
    INT4 iFreqDeriv;
   
    /* 04/12/05 gam; moved from DriveStackSlide.c to StackSlideBinary.c */
    FILE *binaryfp; /* 05/02/17 vir: pointer to output file for sums*/
    FILE *binaryLE;
    FILE *fpSavedEvents; /*05/09/06 vir: file containing events above threshold for each search*/

    char filename[]="myoutbinary.txt";
    char filename2[]="outLE.txt";

    StackSlideSkyParams *csParams; /* Container for ComputeSky */
 
    TdotsAndDeltaTs *pTdotsAndDeltaTs;

    EarthState earth;
    EmissionTime emit;
   
    REAL4 randval;
    RandomParams *randPar=NULL;
    INT4 seed=0;
    
    /* 04/12/05 gam; These all have moved from StackSlideParams struct to here */
    REAL8 *ParamsSMA;
    /* REAL8 *ParamsTperi; */ /* currently unused */
    REAL8 SemiMajorAxis;  /* orbital radius of binary (in sec) */
    REAL8 LoudestEvent;
    REAL8 peakFreq;
    
    printf("Start function StackSlideBinary\n");

    INITSTATUS (status, "StackSlideBinary", STACKSLIDEBINARYC); 
    ATTATCHSTATUSPTR(status);
  
    ParamsSMA=(REAL8 *)LALMalloc(params->numSUMsTotal*sizeof(REAL8)); /*stores SMA for each sum in parameter space*/
  
    /* Allocate space and set quantities for call to ComputeSky() */
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
    pTdotsAndDeltaTs->vecTDots  = (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSTKs);

    if (params->numSpinDown>0) {
        pTdotsAndDeltaTs->vecDeltaTs  = (REAL8 **)LALMalloc(sizeof(REAL8 *)*params->numSTKs);          
        for(i=0;i<params->numSTKs;i++) {
           pTdotsAndDeltaTs->vecDeltaTs[i]  = (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSpinDown);          
        }
    }
    
    /* loop over sky postions */
    for(iSky=0;iSky<params->numSkyPosTotal;iSky++) {
        binaryLE=fopen(filename2, "w");

        params->skyPosData[iSky][0]=params->alphaSX1;
        params->skyPosData[iSky][1]=params->deltaSX1;

        csParams->skyPos[0]=params->skyPosData[iSky][0];/*remember:this must be the same as baryinput.alpha*/
        csParams->skyPos[1]=params->skyPosData[iSky][1];/*remember:this must be the same as baryinput.delta*/
        csParams->baryinput->alpha=params->skyPosData[iSky][0];
        csParams->baryinput->delta=params->skyPosData[iSky][1];
        
        /*ADD HERE BINARY PARAMETERS to be input by means of greg's command line arg.*/
          csParams->OrbitalPeriod=66960; /*no need to go into the commandline*/
          csParams->OrbitalEccentricity=params->OrbitalEccentricity; /*remember: params must be the comm line variable*/
          csParams->ArgPeriapse=params->ArgPeriapse;
          csParams->TperiapseSSB.gpsSeconds=params->TperiapseSSBSec;
          csParams->TperiapseSSB.gpsNanoSeconds=params->TperiapseSSBNanoSec;
        /*these are passed through the command line, MAKE SURE ABOUT CALLING params structure*/
    
        	fpSavedEvents=fopen("SavedEvents.txt","w+");


	  
        /* loop over orbit semimajor axis */
        for(iSMA=0; iSMA < params->nMaxSMA; iSMA++){ /*05/02/18 vir: for each point generate a new random*/
         
         /* loop over Tperi not implemented */
         /* for (iT=0; iT < params->nMaxTperi; iT++)*/
            LALCreateRandomParams(status->statusPtr, &randPar, seed); /*05/02/17 vir: create random params*/
            LALUniformDeviate(status->statusPtr, &randval, randPar);
            /*CHECKSTATUSPTR (status); */

  /*!*/          /*SemiMajorAxis=params->SMAcentral+((REAL8)randval-0.5)*(params->deltaSMA);*/
	    SemiMajorAxis=params->SMAcentral + iSMA*(params->deltaSMA); /*05/09/05 vir */
            
            #ifdef DEBUG_BINARY_CODE  
	    fprintf(stdout,"SMA is %f\n",SemiMajorAxis);
            fflush(stdout);
            #endif
            
	    csParams->SemiMajorAxis=SemiMajorAxis;
            /*params->TperiapseSSBsec=params->Tpericentral+((REAL8)randval-0.5)*(params->deltaTperi); */
            /*05/02/18 vir: csParams->TperiapseSSBsec=params->TperiapseSSBsec; this is now assigned in CL*/

            /* Call STACKSLIDECOMPUTESKYBINARY() */
            StackSlideComputeSkyBinary(status->statusPtr, pTdotsAndDeltaTs, iSkyCoh, csParams);
            CHECKSTATUSPTR (status);

            /* Call STACKSLIDE() for every point in spindown and SMA parameter space*/

            /* for all spindowns, loop */
            for(iFreqDeriv=0;iFreqDeriv<params->numFreqDerivIncludingNoSpinDown;iFreqDeriv++) {

               if (params->numSpinDown>0) {
                   stksldParams->freqDerivData = params->freqDerivData[iFreqDeriv];
               }

               kSUM++; /* increment which SUM we are computing */
               StackSlide(status->statusPtr,params->SUMData,params->STKData,pTdotsAndDeltaTs,stksldParams);
               CHECKSTATUSPTR (status);
               
               ParamsSMA[kSUM]=SemiMajorAxis; /*05/02/18 vir:*/ /* 04/12/05 gam; moved out of StackSlide function */
       
               if((params->nMaxSMA ==1)&&(params->nMaxTperi ==1)) {
                  binaryfp=fopen(filename, "w+");
                  for (k=0; k< params->nBinsPerSUM ; k++) {
                      fprintf(binaryfp,"%f %g\n", params->f0SUM + (REAL8)k*params->dfSUM, params->SUMData[0]->data->data[k]);
                  }
                  fclose(binaryfp);
                  /* 05/02/17 vir: write sum on the file only for no mismatch*/
               } else {
                  /*Look for loudest event in SUMData*/
                  /*05/02/17 vir: for mismatched params, output the loudest event in binaryLE and all events above threshold                       in fpSavedEvents*/
                  #ifdef DEBUG_BINARY_CODE
                     fprintf(stdout,"the SemiMajorAxis is %f\n", SemiMajorAxis);
                     fflush(stdout);
                  #endif
		     
	            FindBinaryLoudest(&LoudestEvent, &peakFreq, params->SUMData, stksldParams, SemiMajorAxis, fpSavedEvents);       		  
                    fprintf(binaryLE,"%f %f %f\n", peakFreq, LoudestEvent, ParamsSMA[kSUM]);
               
	       }/*end of else deltaSMA > 0*/
            
	    } /* for iFreqDeriv = 0 to numFreqDerivTotal */

            LALDestroyRandomParams(status->statusPtr, &randPar );

         /*}end of for iT*/
        }/*end of for iSMA 05/02/17 vir*/
        fclose(binaryLE);    
       fclose(fpSavedEvents);
   
}/*end of for iSky*/
           
    /* Deallocate space for ComputeSky parameters */
    LALFree(csParams->skyPos);
    LALFree(csParams);
        
   /* Deallocate memory */
   /* pTdotsAndDeltaTs->vecDeltaTs is allocated only if numSpinDown > 0. */      
   if (params->numSpinDown>0) {
     for(i=0;i<params->numSTKs;i++) {
       LALFree(pTdotsAndDeltaTs->vecDeltaTs[i]);          
     }
     LALFree(pTdotsAndDeltaTs->vecDeltaTs);
   }
   LALFree(pTdotsAndDeltaTs->vecTDots);
   LALFree(pTdotsAndDeltaTs);
        
   LALFree(ParamsSMA);
   
   CHECKSTATUSPTR (status);
   DETATCHSTATUSPTR (status);

}/*end of void StackSlideBinary()*/

/*Feb 14/05 vir */
/*Start Function FindBinaryLoudest*/
void FindBinaryLoudest(REAL8 *LoudestEvent, REAL8 *peakFreq, REAL4FrequencySeries **SUMData, StackSlideParams *stksldParams, REAL8 SemiMajorAxis, FILE *fpSavedEvents)
{
	
        REAL8 max=0.0;
 
	/*05/09/05 vir: save peaks above threshold*/
/*!*/	REAL8 threshold = 2.8;
/*!*/   REAL8 peak = 0.0;
	
        INT4 iMinSTK=0; 
        INT4 iMaxSTK=stksldParams->nBinsPerSUM;

        INT4 i;
        INT4 indexFreq=0;

        for (i=iMinSTK; i<iMaxSTK; i++)
        {
                if(SUMData[0]->data->data[i] > max)
	       /*	if(SUMData[0]->data->data[i] > threshold )*/
                { 
		  max = SUMData[0]->data->data[i];
                  indexFreq=i;
                }
                *peakFreq=stksldParams->f0SUM + indexFreq*stksldParams->dfSUM;
        }
        *LoudestEvent=max;
       
 for (i=iMinSTK; i<iMaxSTK; i++)
        {
                if(SUMData[0]->data->data[i] > threshold)
                { 
		  peak = SUMData[0]->data->data[i];
                  indexFreq=i;
                
/*                *peakFreq=stksldParams->f0SUM + indexFreq*stksldParams->dfSUM;*/
		  
             fprintf(fpSavedEvents,"%f %f %f \n",stksldParams->f0SUM + indexFreq*stksldParams->dfSUM , peak, SemiMajorAxis);/*05/09/06 vir*/     
	
        #ifdef DEBUG_BINARY_CODE
         fprintf(stdout, " The peak above threshold is %f at freq  %f\n", peak, stksldParams->f0SUM + indexFreq*stksldParams->dfSUM );
         fflush(stdout);
        #endif
		}
	} /*end of for =1 iMinSTk to iMaxSTK */
	
 /*fclose(fpSavedEvents);*/

	}
/**************************************************************************/
/*                 End  Function FindBinaryLoudest                        */
/**************************************************************************/


/**********************************************************************************/
/*         Start Function RunStackSlideBinaryMonteCarloSimulation                 */
/**********************************************************************************/
void RunStackSlideBinaryMonteCarloSimulation(LALStatus *status, StackSlideSearchParams *params, INT4 nSamples)
{
   /*define variables...*/
  
  INT4    i = 0; /* all purpose index */
  INT4    j = 0; /* all purpose index */
  INT4    k = 0; /* all purpose index */

  INT4 iFreq = 0;                           /* 05/26/04 gam */  

  REAL8 cosIota;                    /* cosine of inclination angle iota of the source */
  REAL8 h_0 = params->threshold4;   /* Source amplitude: not sure if threshold 4*/

  LIGOTimeGPS GPSin;            /* reference-time for pulsar parameters at the detector; will convert to SSB! */
  LIGOTimeGPS TperiLIGO;        /*time of the periapse passage in LIGOtime format*/
  LALDetector cachedDetector;

  PulsarSignalParams *pPulsarSignalParams = NULL; /*pointer to the pulsar parameters*/

        REAL4TimeSeries *signal = NULL;
        SFTParams *pSFTParams = NULL;
        LIGOTimeGPSVector timestamps;
	SFTVector *outputSFTs = NULL;  /*pointer to SFTs generated by LALFastgeneratePulsarSignal*/
        FFT **savBLKData = NULL; 
        INT4  *savSumBinMask;          /* 07/15/2005 gam */

	/* 07/14/04 gam; next two are needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  SkyConstAndZeroPsiAMResponse *pSkyConstAndZeroPsiAMResponse;
  SFTandSignalParams *pSFTandSignalParams;


  REAL8 f0SUM = 0.0;                        /* 05/26/04 gam */
  REAL8 bandSUM = 0.0;                      /* 05/26/04 gam */
  INT4 nBinsPerSUM = 0;                     /* 05/26/04 gam */
  INT4 nBinsPerSUMm1 = -1;                  /* 05/26/04 gam */  
  REAL8 real8NBinsPerSUM;                   /* 07/15/05 gam */
  INT4 firstNonExcludedBin = -1;            /* 05/19/05 gam */
  INT4 lastNonExcludedBin  = -1;            /* 05/19/05 gam */

  INT4  *arrayAvailableFreqBins; /* 07/15/2005 gam */
  INT4  nAvailableBins;          /* 07/15/2005 gam */
  REAL8 real8NAvailableBins;     /* 07/15/2005 gam */

  

   INITSTATUS( status, "RunStackSlideBinaryMonteCarloSimulation", STACKSLIDEBINARYC );
   ATTATCHSTATUSPTR(status);

printf("!!!!!!!!!start function MC!!!!!!!!!!!!!!!\n");

/*allocate memory for initial BLKs of data 05/09/07 vir*/
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


/* Allocate memory for PulsarSignalParams and initialize */
  pPulsarSignalParams = (PulsarSignalParams *)LALMalloc(sizeof(PulsarSignalParams));
  pPulsarSignalParams->pulsar.position.system = COORDINATESYSTEM_EQUATORIAL;
  pPulsarSignalParams->pulsar.spindown = NULL;

/*Allocate memory for orbital parameters*/
  pPulsarSignalParams->orbit=(BinaryOrbitParams *)LALMalloc(sizeof(BinaryOrbitParams));

  pPulsarSignalParams->transfer = NULL;  

  
  /*set detector site*/
  
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
      ABORT( status, STACKSLIDEBINARYH_EIFO, STACKSLIDEBINARYH_MSGEIFO);
  }    
  pPulsarSignalParams->site = &cachedDetector;     
  pPulsarSignalParams->ephemerides = params->edat;

  pPulsarSignalParams->startTimeGPS.gpsSeconds = (INT4)params->gpsStartTimeSec;
  pPulsarSignalParams->startTimeGPS.gpsNanoSeconds = (INT4)params->gpsStartTimeNan;
  pPulsarSignalParams->duration = (UINT4)params->duration;
  pPulsarSignalParams->samplingRate = (REAL8)ceil(2.0*params->bandBLK); /* Make sampleRate an integer so that T*samplingRate = integer for integer T */
  pPulsarSignalParams->fHeterodyne = params->f0BLK;  
  /* Find the time at the SSB that corresponds to the arrive time at the detector of first data requested. */
  GPSin.gpsSeconds = (INT4)params->gpsEpochStartTimeSec;     /* 12/06/04 gam; GPS epoch Sec that gives reference time in SSB */
  GPSin.gpsNanoSeconds = (INT4)params->gpsEpochStartTimeNan; /* 12/06/04 gam; GPS epoch Nan that gives reference time in SSB */
TperiLIGO.gpsSeconds=731163327;
TperiLIGO.gpsNanoSeconds = 0;
  
  /* Allocate memory for SFTParams and initialize */
  pSFTParams = (SFTParams *)LALMalloc(sizeof(SFTParams));
  pSFTParams->Tsft = params->tBLK;
  timestamps.length = params->numBLKs;
  timestamps.data = params->timeStamps; 
  pSFTParams->timestamps = &timestamps;
  pSFTParams->noiseSFTs = NULL;  

  /*NOT SURE if needed*/
  /* 05/26/04 gam; initialize variables that keep track of which frequency we are working on */  
  f0SUM = params->f0SUM;
  bandSUM = params->bandSUM;
  nBinsPerSUM = params->nBinsPerSUM;
  real8NBinsPerSUM = ((REAL8)nBinsPerSUM); /* 07/15/05 gam */
  nBinsPerSUMm1 = nBinsPerSUM - 1;
  params->keepThisNumber = 1;      /* During MC only keep 1 event; this will be the injected event */


/*generate random numeber seed here ?*/

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
     pPulsarSignalParams->samplingRate = 2.0*params->bandBLK; /* can set samplingRate to exactly 2.0*bandBLK when using LALFastGeneratePulsarSFTs */
     GPSin.gpsSeconds = params->timeStamps[0].gpsSeconds;         /* 08/02/04 gam */
     GPSin.gpsNanoSeconds = params->timeStamps[0].gpsNanoSeconds; /* 08/02/04 gam */
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

  /*......Search surrounding points.........*/
  
/*Store frequency bins where a signal is injected into an array */

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
     ABORT( status, STACKSLIDEBINARYH_ENOFREQBINS, STACKSLIDEBINARYH_MSGENOFREQBINS);
  }
  real8NAvailableBins = ((REAL8)nAvailableBins);
  
#ifdef STACKSLIDEBINARY_DEBUG
printf("%f\n",real8NAvailableBins);
#endif


 /***********************************************************/
 /*                                                         */
 /* START SECTION: LOOP OVER ENTIRE MONTE CARLO SIMULATIONS */
 /*                                                         */
 /***********************************************************/

/*start loop over points in par space kSUM */

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

/*choose a frequency for injection */

    pPulsarSignalParams->pulsar.f0 = 480.0; /*will become random*/

    /*nuisance parameters*/
         pPulsarSignalParams->pulsar.psi = params->orientationAngle;
         cosIota =params->cosInclinationAngle;        /*will become random if testFlag=8*/

    /* h_0 is fixed equal to params->threshold4 above; get A_+ and A_x from h_0 and random cosIota */
    pPulsarSignalParams->pulsar.aPlus = (REAL4)(0.5*h_0*(1.0 + cosIota*cosIota));
    pPulsarSignalParams->pulsar.aCross = (REAL4)(h_0*cosIota);

 
 
    /*phi0, SMA Tperi ecc ArgPeri*/
    pPulsarSignalParams->pulsar.phi0 = 0.0;
    pPulsarSignalParams->pulsar.position.longitude = params->alphaSX1;
    pPulsarSignalParams->pulsar.position.latitude = params->deltaSX1;

   
           if (pPulsarSignalParams->orbit != NULL)
                  {
                  pPulsarSignalParams->orbit->omega = 1.0;
                  pPulsarSignalParams->orbit->orbitEpoch = TperiLIGO;
                  pPulsarSignalParams->orbit->oneMinusEcc = 1.0;
                  pPulsarSignalParams->orbit->angularSpeed = 0.00001;
                  pPulsarSignalParams->orbit->rPeriNorm = 1.44;
#ifdef DEBUG_LALFASTGENERATEPULSARSFTS
		  fprintf(stdout,"omega %f, Tperi %i, ecc %f, angspeed %f \n", pPulsarSignalParams->orbit->omega,pPulsarSignalParams->orbit->orbitEpoch.gpsSeconds , pPulsarSignalParams->orbit->oneMinusEcc,  pPulsarSignalParams->orbit->angularSpeed = 0.00001 );
                  }
#endif

	   /* 07/14/04 gam; check if using LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
    if ( (params->testFlag & 4) > 0 ) {
        #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
          fprintf(stdout,"About to call LALFastGeneratePulsarSFTs \n");
          fflush(stdout);
        #endif
        /* 07/14/04 gam; use SkyConstAndZeroPsiAMResponse from LALComputeSkyAndZeroPsiAMResponse and SFTandSignalParams to generate SFTs fast. */
        LALFastGeneratePulsarSFTs(status->statusPtr, &outputSFTs, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
        CHECKSTATUSPTR (status);
    }
          /*fprintf(stdout,"%i %g \n",j,outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im);
          fflush(stdout);*/

/*should check memory allocation for outputSFTs */

         
        #ifdef PRINT_MAXPOWERANDBINEACHSFT
        /*{
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
         }*/
        #endif
 
    
        /* Add outputSFTs with injected signal to input noise SFTs; no renorm should be needed */
/*        for(i=0;i<params->numBLKs;i++) {
          for(j=0;j<params->nBinsPerBLK;j++) {
             #ifdef PRINTCOMPARISON_INPUTVSMONTECARLOSFT_DATA
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].re, outputSFTs->data[%i].data->data[%i].re = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[j].re,outputSFTs->data[i].data->data[j].re);
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].im, outputSFTs->data[%i].data->data[%i].im = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[j].im,outputSFTs->data[i].data->data[j].im);
               fflush(stdout);
             #endif
             params->BLKData[i]->fft->data->data[j].re = savBLKData[i]->fft->data->data[j].re + outputSFTs->data[i].data->data[j].re;
             params->BLKData[i]->fft->data->data[j].im = savBLKData[i]->fft->data->data[j].im + outputSFTs->data[i].data->data[j].im;
          }
        }*/
    /* next is else for if ( (params->testFlag & 4) > 0 ) else... */

/* 
  for iSMA = iSMAmin to iSMAmax { 
             for iT = iTmin to iTmax {  
             
	     SMA = rand1;
	     Tperi=rand2;
	          LALFastGeneratePulsarSFTs() 
                                      }
                                 }      
*/
CHECKSTATUSPTR(status);
DETATCHSTATUSPTR(status);
}

/**********************************************************************************/
/*         End Function RunStackSlideBinaryMonteCarloSimulation                   */
/**********************************************************************************/




