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
 #define PRINT_RES
 #define PRINT_INJECTION_VALUES
/* #define DEBUG_LALFASTGENERATEPULSARSFTS*/
 #define DEBUG_LALGENERATEPULSARSFTS
/*#define DEBUG_CONFIDENCE*/
/* #define PRINT_MAXPOWERANDBINEACHSFT*/
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
    char filename2[]="outLE.txt"; /* this must be passed from CLA */

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
    REAL4 SNR;
    
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
          csParams->OrbitalPeriod=68023; /*66960 no need to go into the commandline*/
          csParams->OrbitalEccentricity=params->OrbitalEccentricity; /*remember: params must be the comm line variable*/
          csParams->ArgPeriapse=params->ArgPeriapse;
          csParams->TperiapseSSB.gpsSeconds=params->TperiapseSSBSec;
          csParams->TperiapseSSB.gpsNanoSeconds=params->TperiapseSSBNanoSec;
        /*these are passed through the command line, MAKE SURE ABOUT CALLING params structure*/
    
        	fpSavedEvents=fopen("SavedEvents.txt","w+");


	  
        /* loop over orbit semimajor axis */
        for(iSMA=0; iSMA < params->nMaxSMA; iSMA++){ /*05/02/18 vir: for each point generate a new random*/
         
         /* loop over Tperi not implemented for S4 search */
         /* for (iT=0; iT < params->nMaxTperi; iT++)*/
            LALCreateRandomParams(status->statusPtr, &randPar, seed); /*05/02/17 vir: create random params*/
            LALUniformDeviate(status->statusPtr, &randval, randPar);
            /*CHECKSTATUSPTR (status); */

  /*!*/          /*SemiMajorAxis=params->SMAcentral+((REAL8)randval-0.5)*(params->deltaSMA);*/
	    SemiMajorAxis=params->SMAcentral + (iSMA - 2.0 )*(params->deltaSMA); /*05/09/05 vir */
            
                 /*this is so because in a multiple filter search we need at most 4/5 filters to*/
	         /*cover >1 sigma parameter space, these filters will range between 1.24 to 1.54/1.64 at deltaSMA=0.1*/
	    

            #ifdef DEBUG_BINARY_CODE  
	    fprintf(stdout,"SMA is %f\n",SemiMajorAxis);
            /*fflush(stdout);*/
            #endif
            
	      
	    csParams->SemiMajorAxis=SemiMajorAxis;
            /*params->TperiapseSSBsec=params->Tpericentral+((REAL8)randval-0.5)*(params->deltaTperi); */
            /*05/02/18 vir: csParams->TperiapseSSBsec=params->TperiapseSSBsec; this is now assigned in CL*/

	    
	    
            /* Call STACKSLIDECOMPUTESKYBINARY() */
            StackSlideComputeSkyBinary(status->statusPtr, pTdotsAndDeltaTs, iSkyCoh, csParams);
            CHECKSTATUSPTR (status);

	    

            /* Call STACKSLIDE() for every point in spindown and SMA parameter space*/
	    
         /*printf("numFreqDerivIncludingNoSpinDown %i\n", params->numFreqDerivIncludingNoSpinDown);*/
            
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
                 
		  /*Look for loudest also in this case ....*/

		  FindBinaryLoudest(&LoudestEvent, &peakFreq, params->SUMData, stksldParams, SemiMajorAxis, fpSavedEvents);                    fprintf(binaryLE,"%f %f %f\n", peakFreq, LoudestEvent, ParamsSMA[kSUM]);
		  
        
        	 /*05/11/17 add calculation of standard deviation only if MC active*/
                   if( (params->testFlag & 2) > 0 ){  
	           ValidateMCResults(status->statusPtr, params->SUMData[0], params, &SNR);
                   CHECKSTATUSPTR (status);
                   		   
		   }
		 
		  /* 05/02/17 vir: write sum on the file only for no mismatch*/
               } else {
                  /*Look for loudest event in SUMData*/
                  /*05/02/17 vir: for mismatched params, output the loudest event in binaryLE and all events above threshold                       in fpSavedEvents*/
                  #ifdef DEBUG_BINARY_CODE
                     fprintf(stdout,"the SemiMajorAxis is %f\n", SemiMajorAxis);
                     fflush(stdout);
                  #endif
                    /* this finds the peaks above threshold*/ 
	            FindBinaryLoudest(&LoudestEvent, &peakFreq, params->SUMData, stksldParams, SemiMajorAxis, fpSavedEvents);       		  /*the file outLE created by this function contains the loudest for each filter and corr. freq. */
                    fprintf(binaryLE,"%f %f %f\n", peakFreq, LoudestEvent, ParamsSMA[kSUM]);

		   
	       }/*end of else nMaxSMA == 1 */
          
	       /*Allocate here space for xml table*/
	       
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
                   fprintf(stdout,"end function StackSlideBinary\n");
		   fflush(stdout);

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
       
	printf("LE %f\n", max);
	
 for (i=iMinSTK; i<iMaxSTK; i++)
        {
                if(SUMData[0]->data->data[i] > threshold)
                { 
		  peak = SUMData[0]->data->data[i];
                  indexFreq=i;
                
/*                *peakFreq=stksldParams->f0SUM + indexFreq*stksldParams->dfSUM;*/
		  
             fprintf(fpSavedEvents,"%f %f %f \n",stksldParams->f0SUM + indexFreq*stksldParams->dfSUM , peak, SemiMajorAxis);/*05/09/06 vir*/     
	
        #ifdef DEBUG_BINARY_CODE_1
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
  
  
  INT4    i = 0; 
  INT4    j = 0; 
  INT4    k = 0; 
  
  INT4 kSUM = 0;    
  
  FILE *outfp;
  
  INT4 iFreq = 0;                           

  REAL8 cosIota;                    
  REAL8 h_0 = 0;  /* params->threshold4;    Source amplitude: not sure if threshold 4*/
  
  REAL4 renorm;                 /* Need to renormalize SFTs to account for different sample rates */  
  LIGOTimeGPS GPSin;            /* reference-time for pulsar parameters at the detector; will convert to SSB! */
  LIGOTimeGPS TperiLIGO;        /*time of the periapse passage in LIGOtime format*/
  LALDetector cachedDetector;

  PulsarSignalParams *pPulsarSignalParams = NULL; /*pointer to the pulsar parameters*/

        REAL4TimeSeries *signal = NULL;
        SFTParams *pSFTParams = NULL;
        LIGOTimeGPSVector timestamps;
	SFTVector *outputSFTs = NULL;  /*pointer to SFTs generated by LALFastgeneratePulsarSignal*/
        FFT **savBLKData = NULL; 
        INT4  *savSumBinMask;          

   /* 07/14/04 gam; next two are needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  SkyConstAndZeroPsiAMResponse *pSkyConstAndZeroPsiAMResponse;
  SFTandSignalParams *pSFTandSignalParams;

  BOOLEAN BinsCanBeExcluded = 0;
  BOOLEAN searchSurroundingPts = 0;
  BOOLEAN searchSurroundingOrbitalPts = 0;
  /*use results from a previous job in the pipeline*/
  BOOLEAN reportMCResults = 0;
 /* REAL4 priorLoudestEvent = 0.0;
  REAL8 priorStartFreq = 0.0;
  REAL8 priorBand = 0.0;*/
  REAL8 tmp_Confidence = 0.0;
  REAL8 tmp_conf_err = 0.0;
  REAL8 SearchLoudestEvent = 0.0;
  
  /*REAL8 priorConfidence = 0.0;
  REAL8 priorUL = 0.0;
  REAL8 priorUncertainty = 0.0;*/

  REAL8 *arrayULs = NULL;
  REAL8 *arrayConfs = NULL;
  INT4 *arrayConverged = NULL;

  
  BOOLEAN MCInjectRandomSMA = 1;                    /* 05/09/12 vir : will be moved into command line */
  BOOLEAN MCInjectRandomTperi = 1;                  /* 05/09/12 vir : will be moved into command line */
  BOOLEAN MCInjectRandomSignalPars = 1;  
  
  REAL8 SemiMajorAxis;
  REAL8 HalfAmplSMAParSpace = 0.18;         /* 05/09/12 vir :in sec */
  REAL8 HalfAmplTperiParSpace = 300;        /* 05/09/12 vir :in sec */
  REAL8 ClosestSMATemplate;
  UINT4 ClosestTperiTemplate;

  REAL8 f0SUM = 0.0;                        
  REAL8 bandSUM = 0.0;                      
  INT4 nBinsPerSUM = 0;                     
  INT4 nBinsPerSUMm1 = -1;                    
  REAL8 real8NBinsPerSUM;                   
  INT4 firstNonExcludedBin = -1;            
  INT4 lastNonExcludedBin  = -1;            
  INT4 InjFreqBin = 0;
  
  REAL4 randval;
  REAL8 real8RandVal;
  RandomParams *randPar=NULL;

  INT4  *arrayAvailableFreqBins; 
  INT4  nAvailableBins;          
  REAL8 real8NAvailableBins;     

  REAL8 *nAboveLE = NULL;
  REAL8 *nTrials = NULL;

  
  
 INT4 onePlusNumMCRescalings = 1 + params->numMCRescalings; /* 09/01/05 gam; one plus number of times to rescale fake SFTs */
 INT4 twoPlusNumMCRescalings = 2 + params->numMCRescalings; 

 
  FILE *fph0;
/*  CHAR h0file[256];*/
  INT4 Nh0;
  REAL8 *h0data;

   REAL8 *h_0Rescaled = NULL; 
   
   INT4 iRescale = 0;
    INT4 jStart = 0;
    INT4 jSavedBLK = 0;

  BOOLEAN breakOutOfLoop;
  REAL8 maxMCErr = params->maxMCErr;
  INT4 seed=0;
  REAL8 thisMCErr = 1;           /* keep track of last absolute error in confidence. */
  REAL8 lastMCErr = 1;
  INT4 iMC = 0;             /* which entire MC simulation is currently running */
  INT4 countMC = 0;         /* count of number of completed MC simulations, including rescalings */
  INT4 maxMC = 1;  /* default number of MC iterations to perform; set equal to params->maxMCinterations when interating MC */
  REAL8 minh_0 = 0; 
 

        FILE *fpMCInjections; /*modify this into a name containing which injection we are working on*/
        FILE *fpReadFromMCInjectionFile;
        REAL8 StartFreq;
        REAL4 LEvent;
	REAL8 OrbRadius;
        REAL8 Conf = 0.95;        

	/*store final results*/
        FILE *fpMonteCarloResults = NULL;
	CHAR MonteCarloResults[256];
	
  
   INITSTATUS( status, "RunStackSlideBinaryMonteCarloSimulation", STACKSLIDEBINARYC );
   ATTATCHSTATUSPTR(status);

printf("!!!!!!!!!start function MC!!!!!!!!!!!!!!!\n");


if (maxMC < 1) {
      ABORT( status, STACKSLIDEBINARYH_EMAXMC, STACKSLIDEBINARYH_MSGEMAXMC);
  }

if ( (params->deltaSMA > 0) && (params->nMaxSMA == 1) ){
     ABORT (status, STACKSLIDEBINARYH_EDELTASMAnSMA, STACKSLIDEBINARYH_MSGDELTASMAnSMA);
  }


if ( (params->testFlag & 16) > 0 ) {
    
	reportMCResults = 1;

    /* 09/01/05 gam; allocate memory dependent on the number of rescalings */
    nAboveLE = (REAL8 *)LALMalloc(onePlusNumMCRescalings*sizeof(REAL8)); /*REMOVE!!!!!!!!!!!!!!!!!!!!!!!!*/
    nTrials  = (REAL8 *)LALMalloc(onePlusNumMCRescalings*sizeof(REAL8));


}

    /*open MonteCarloResults file*/
     strcpy(MonteCarloResults, params->priorResultsFile);
     fpMonteCarloResults=fopen(MonteCarloResults, "w+");
     
    /*open h0File if you don't want injection from command line*/
             if( params->threshold4 == 0 ){
       Readh0File(status->statusPtr, fph0, &Nh0, &h0data);
       CHECKSTATUSPTR (status);
	     }else{
	     h_0=params->threshold4;
	     }

if ( (params->testFlag & 32) > 0 ) {
      /* 05/24/05 gam; iterate MC to converge on desired confidence */
      maxMC = params->maxMCinterations;
  }



   if( (params->threshold4) == 0){
    maxMC = Nh0;     /*05/11/27 I want as many iterations as the number of h_0 values*/
    printf("MaxMC %i\n", maxMC);
   }

if (params->numMCRescalings > 0) {
      arrayULs   = (REAL8 *)LALMalloc((twoPlusNumMCRescalings*maxMC)*sizeof(REAL8));
      arrayConfs = (REAL8 *)LALMalloc((twoPlusNumMCRescalings*maxMC)*sizeof(REAL8));
      arrayConverged = (INT4 *)LALMalloc((twoPlusNumMCRescalings*maxMC)*sizeof(INT4));
    } else {

      arrayULs   = (REAL8 *)LALMalloc(maxMC*sizeof(REAL8));
      arrayConfs = (REAL8 *)LALMalloc(maxMC*sizeof(REAL8));
      arrayConverged = (INT4 *)LALMalloc(maxMC*sizeof(INT4));
    }

/*  if ( (params->numMCRescalings > 0) && (reportMCResults != 1)) {
      ABORT( status, STACKSLIDEBINARYH_E2NUMMCRESCALINGS, STACKSLIDEBINARYH_MSGE2NUMMCRESCALINGS);
  }*/


h_0Rescaled = (REAL8 *)LALMalloc(onePlusNumMCRescalings*sizeof(REAL8));



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

  TperiLIGO.gpsSeconds = params->TperiapseSSBSec;
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
/*     printf("params->sumBinMask[%d], %d\n", iFreq, params->sumBinMask[iFreq]);*/
     if (params->sumBinMask[iFreq] != 0) {
        if (firstNonExcludedBin == -1) {
           firstNonExcludedBin = iFreq;
        }
        lastNonExcludedBin = iFreq;
     }
  }
  printf("nBinsPerSUM %i firstNonExcludedBin %i lastNonExcludedBin %i\n",nBinsPerSUM,firstNonExcludedBin, lastNonExcludedBin);
  LALFree(params->sumBinMask); /* 07/15/2005 gam; reallocate and save just 1 */


  
  /*if(Search surrounding points)*/
  /* 07/29/05 gam */
  if (searchSurroundingPts) {
    params->sumBinMask=(INT4 *)LALMalloc(2*sizeof(INT4));
    params->sumBinMask[0] = 1;   /* Default value; will check when injected frequency is set up below */
    params->sumBinMask[1] = 1;   /* Default value; will check when injected frequency is set up below */
  } else {
    params->sumBinMask=(INT4 *)LALMalloc(sizeof(INT4));
    params->sumBinMask[0] = 1;   /* Will only inject next bins where savSumBinMask != 0. */
  }  /* end if(Search surrounding points)*/
 

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
fprintf(stdout, "real8NAvailableBins %f\n",real8NAvailableBins);
fflush(stdout);
#endif


 /***********************************************************/
 /*                                                         */
 /* START SECTION: LOOP OVER ENTIRE MONTE CARLO SIMULATIONS */
 /*                                                         */
 /***********************************************************/

printf("!!!!!!!!!!! MaxMC %i\n\n", maxMC);


for(iMC=0;iMC<maxMC;iMC++) {

	/* This is the big loop that reruns the entire MC maxMC times or until MC has converged */


/*05/10/18 vir: commented out*/
    if (reportMCResults) {  
        for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) {
        nAboveLE[iRescale] = 0.0;
        nTrials[iRescale] = 0.0;
       }
     }


    /*05/10/03 vir: open file that will store the MCInjectionresults for every set of injections*/  
    fpMCInjections=fopen("MCInjectionResults.txt", "w+");

   
   if (iMC > Nh0) break; 
   
   if ( params->threshold4 == 0 )
   {
	  
	       h_0 = h0data[iMC];
	    printf("xxxxxxxxxxxxx h_0 = %g\n\n", h0data[iMC]);

	    
      }else{
     h_0=params->threshold4 ;
      }
 
   /*  printf("h_0 %g\n", h_0);*/
      
  /**********************************************/
  /*                                            */
  /* START SECTION: MONTE CARLO INJECTIONS LOOP */
  /*                                            */
  /**********************************************/
  
    for(kSUM=0;kSUM<params->numMCInjections;kSUM++) {

	  
    params->whichMCSUM = kSUM; /* kSUM is which injection we are working on */  

    if (MCInjectRandomSMA) { /*then generate a random SMA for injection*/
      LALCreateRandomParams(status->statusPtr, &randPar, seed); 
      LALUniformDeviate(status->statusPtr, &randval, randPar); 
      SemiMajorAxis=params->SMAcentral+((REAL8)randval-0.5)*(HalfAmplSMAParSpace);
    }else{
       /* SemiMajorAxis=params->SMAcentral;  */
       SemiMajorAxis = 1.44;
    }
 

    if (MCInjectRandomTperi ) {
	  LALCreateRandomParams(status->statusPtr, &randPar, seed); 
          LALUniformDeviate(status->statusPtr, &randval, randPar); 
          TperiLIGO.gpsSeconds=params->TperiapseSSBSec+((REAL8)randval-0.5)*(HalfAmplTperiParSpace);
    }else{
       /* TperiLIGO.gpsSeconds=params->TperiapseSSBSec;*/
	  TperiLIGO.gpsSeconds = 731163327;
    }
	/*..............................*/
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
    


    
    if(BinsCanBeExcluded){
 /* FIND RANDOM FREQUENCY THAT IS NOT EXCLUDED */
     LALCreateRandomParams(status->statusPtr, &randPar, seed); 
     LALUniformDeviate(status->statusPtr, &randval, randPar); 
     real8RandVal=randval*real8NAvailableBins;
     i = floor( real8RandVal );    /* nearest bin arrayAvailableFreqBins */
     if (i >= nAvailableBins) i = nAvailableBins - 1; if (i < 0) i = 0; /* make sure i in range */
     iFreq = arrayAvailableFreqBins[i];  /* randomly chosen available frequency bin */
     fprintf(stdout,"iFreq %i real8RandVal %f\n", iFreq, real8RandVal);  
    
     params->f0SUM = f0SUM + iFreq*params->dfSUM; /* search f0*/
    
    }else{
    
	    /*randomly a bin in the SUM where to inject signal*/
	LALCreateRandomParams(status->statusPtr, &randPar, seed); 
        LALUniformDeviate(status->statusPtr, &randval, randPar); 
        InjFreqBin = floor(randval*(params->nBinsPerSUM));
	if(InjFreqBin >= params->nBinsPerSUM) InjFreqBin = params->nBinsPerSUM - 1;
       
	params->f0SUM=f0SUM ; /*05/10/19 vir: the start of the sum remains the same*/
    }
         LALCreateRandomParams(status->statusPtr, &randPar, seed); 
         LALUniformDeviate(status->statusPtr, &randval, randPar); 
         real8RandVal=params->f0SUM - (0.5*params->dfSUM) + randval*params->dfSUM + InjFreqBin*params->dfSUM;
	 pPulsarSignalParams->pulsar.f0 = real8RandVal; 

     #ifdef PRINT_INJECTION_VALUES
	 fprintf(stdout, "f0SUM %f, pPulsarSignalParams->pulsar.f0 %f InjFreqBin %d\n", f0SUM, pPulsarSignalParams->pulsar.f0, InjFreqBin);
     #endif
	 
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



	 
                               /*nuisance parameters*/

	 if ( (params->testFlag & 8) > 0 || (MCInjectRandomSignalPars == 0) ) {
         pPulsarSignalParams->pulsar.psi = params->orientationAngle;
         cosIota =params->cosInclinationAngle;
    } else {
          /*rand psi*/
	  LALCreateRandomParams(status->statusPtr, &randPar, seed); 
          LALUniformDeviate(status->statusPtr, &randval, randPar); 
          pPulsarSignalParams->pulsar.psi = (randval - 0.5) * ((REAL4)LAL_PI_2);
          /*rand cosiota*/
	  LALCreateRandomParams(status->statusPtr, &randPar, seed); 
          LALUniformDeviate(status->statusPtr, &randval, randPar); 
          cosIota = 2.0*((REAL8)randval) - 1.0;
         
    }
    
	 pPulsarSignalParams->pulsar.aPlus = (REAL4)(0.5*h_0*(1.0 + cosIota*cosIota));
         pPulsarSignalParams->pulsar.aCross = (REAL4)(h_0*cosIota);

 
 
    /*orbital parameters*/
    if(MCInjectRandomSignalPars){
    LALUniformDeviate(status->statusPtr, &randval, randPar); 
    pPulsarSignalParams->pulsar.phi0 = ((REAL8)randval) * ((REAL8)LAL_TWOPI);
    }else{
      pPulsarSignalParams->pulsar.phi0 = 0.0;
    }
    
    pPulsarSignalParams->pulsar.position.longitude = params->alphaSX1;
    pPulsarSignalParams->pulsar.position.latitude = params->deltaSX1;

   
           if (pPulsarSignalParams->orbit != NULL)
                  {
                  pPulsarSignalParams->orbit->omega = 0.0;
                  pPulsarSignalParams->orbit->orbitEpoch = TperiLIGO;
                  pPulsarSignalParams->orbit->oneMinusEcc = 1.0;
                  pPulsarSignalParams->orbit->angularSpeed = 9.2368E-5; /*06/02/27 changed*/ 
                  pPulsarSignalParams->orbit->rPeriNorm = SemiMajorAxis;
    #ifdef PRINT_INJECTION_VALUES
		  fprintf(stdout,"inj pars are omega %f, Tperi %i, ecc %f, angspeed %f, rPeriNorm %f, phi0 %f\n", pPulsarSignalParams->orbit->omega,pPulsarSignalParams->orbit->orbitEpoch.gpsSeconds , pPulsarSignalParams->orbit->oneMinusEcc,  pPulsarSignalParams->orbit->angularSpeed , pPulsarSignalParams->orbit->rPeriNorm,  pPulsarSignalParams->pulsar.phi0 );
    #endif
                  }


	   /* 07/14/04 gam; check if using LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
    if ( (params->testFlag & 4) > 0 ) {  
        #ifdef DEBUG_LALFASTGENERATEPULSARSFTS
          fprintf(stdout,"About to call LALFastGeneratePulsarSFTs \n");
          fflush(stdout);
        #endif
        /* 07/14/04 gam; use SkyConstAndZeroPsiAMResponse from LALComputeSkyAndZeroPsiAMResponse 
	 * and SFTandSignalParams to generate SFTs fast. */
        LALFastGeneratePulsarSFTs(status->statusPtr, &outputSFTs, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
        CHECKSTATUSPTR (status);
	
    
 
         
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
/*        for(i=0;i<params->numBLKs;i++) {
          for(j=0;j<params->nBinsPerBLK;j++) {
             #ifdef PRINTCOMPARISON_INPUTVSMONTECARLOSFT_DATA
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].re, 
	       outputSFTs->data[%i].data->data[%i].re = %g, %g\n",i,j,i,j,
	       savBLKData[i]->fft->data->data[j].re,outputSFTs->data[i].data->data[j].re);
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].im,
	       outputSFTs->data[%i].data->data[%i].im = %g, %g\n",i,j,i,j,
	       savBLKData[i]->fft->data->data[j].im,outputSFTs->data[i].data->data[j].im);
               fflush(stdout);
             #endif
             params->BLKData[i]->fft->data->data[j].re = savBLKData[i]->fft->data->data[j].re + outputSFTs->data[i].data->data[j].re;
             params->BLKData[i]->fft->data->data[j].im = savBLKData[i]->fft->data->data[j].im + outputSFTs->data[i].data->data[j].im;
          }*/
        }else{  /* next is else for if ( (params->testFlag & 4) > 0 ) else... */

       #ifdef PRINT_INJECTION_VALUES

	printf("use LALGeneratePulsarSignal and LALSignalToSFTs\n");
        printf("GPSin.Seconds %i\n", GPSin.gpsSeconds);
	printf(" pPulsarSignalParams->pulsar.tRef.gpsSeconds %i\n",  pPulsarSignalParams->pulsar.tRef.gpsSeconds);
       #endif
	
        signal = NULL; /* Call LALGeneratePulsarSignal to generate signal */
        LALGeneratePulsarSignal(status->statusPtr, &signal, pPulsarSignalParams);
        CHECKSTATUSPTR (status);
       
	#ifdef PRINT_INJECTION_VALUES
	fprintf(stdout,"signal->deltaT = %23.10e \n",signal->deltaT);
        fprintf(stdout,"signal->epoch.gpsSeconds = %i \n",signal->epoch.gpsSeconds);
        fprintf(stdout,"signal->epoch.gpsNanoSeconds = %i \n",signal->epoch.gpsNanoSeconds);
        fprintf(stdout,"pPulsarSignalParams->duration*pPulsarSignalParams->samplingRate (# of points) %g, signal->data->length  %i\n",pPulsarSignalParams->duration*pPulsarSignalParams->samplingRate,signal->data->length);
        fprintf(stdout,"pPulsarSignalParams->samplingRate = %g \n",pPulsarSignalParams->samplingRate);
        fprintf(stdout,"pSFTandSignalParams->nSamples = %i \n",nSamples);
        #endif    
        
	outputSFTs = NULL; /* Call LALSignalToSFTs to generate outputSFTs with injection data */
        LALSignalToSFTs(status->statusPtr, &outputSFTs, signal, pSFTParams);
        CHECKSTATUSPTR (status);

	/*(whichMCsum =0 )*/

              #ifdef PRINT_INJECTION_VALUES
              fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",InjFreqBin,h_0);
              fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",
			      InjFreqBin,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
              fprintf(stdout,"iFreq = %i, inject psi = %23.10e \n",InjFreqBin,pPulsarSignalParams->pulsar.psi);
              fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",InjFreqBin,pPulsarSignalParams->pulsar.phi0);
              fprintf(stdout,"iFreq = %i, search f0 = %23.10e, inject f0 = %23.10e \n",
			      InjFreqBin,f0SUM + InjFreqBin*params->dfSUM,pPulsarSignalParams->pulsar.f0);
             fprintf(stdout,"outputSFTs->data[0].data->length = %i \n",outputSFTs->data[0].data->length);
               #endif

	     /* CHECK the index i*/
             
	     i=0;
	     renorm = ((REAL4)nSamples)/((REAL4)(outputSFTs->data[i].data->length - 1));
             printf("renorm %g\n", renorm);
	     

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
        /*renorm = ((REAL4)nSamples)/((REAL4)(outputSFTs->data[0].data->length - 1));*/ /* 05/21/04 gam; should be the same for all outputSFTs; note minus 1 is needed */
	     
          for(i=0;i<params->numBLKs;i++) {
          for(j=0;j<params->nBinsPerBLK;j++) {
             #ifdef PRINTCOMPARISON_INPUTVSMONTECARLOSFT_DATA
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].re, outputSFTs->data[%i].data->data[%i].re = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[j].re,renorm*outputSFTs->data[i].data->data[j].re);
               fprintf(stdout,"savBLKData[%i]->fft->data->data[%i].im, outputSFTs->data[%i].data->data[%i].im = %g, %g\n",i,j,i,j,savBLKData[i]->fft->data->data[j].im,renorm*outputSFTs->data[i].data->data[j].im);
               fflush(stdout);
             #endif
             params->BLKData[i]->fft->data->data[j].re = savBLKData[i]->fft->data->data[j].re + renorm*outputSFTs->data[i].data->data[j].re;
             params->BLKData[i]->fft->data->data[j].im = savBLKData[i]->fft->data->data[j].im + renorm*outputSFTs->data[i].data->data[j].im;
          }
        }
	  
fprintf(stdout,"nBinsPerBLK %i, deltaF %f, f0 %g\n",params->nBinsPerBLK, params->BLKData[0]->fft->deltaF, params->BLKData[0]->fft->f0 );
	  
	  /* if only 1 MC sum {*/
  
    #ifdef PRINT_RES   
   outfp=fopen("OutputBLKplusInj.txt", "w+");

    for (j=0;j<params->nBinsPerBLK;j++){
/*fprintf(outfp,"%g %g\n",params->BLKData[0]->fft->f0+j*params->BLKData[0]->fft->deltaF,
 * params->BLKData[0]->fft->data->data[j].re * params->BLKData[0]->fft->data->data[j].re + params->BLKData[0]->fft->data->data[j].im * params->BLKData[0]->fft->data->data[j].im);*/
/*fprintf(outfp,"%g %g\n",savBLKData[0]->fft->f0+j*savBLKData[0]->fft->deltaF,
 * savBLKData[0]->fft->data->data[j].re * savBLKData[0]->fft->data->data[j].re + savBLKData[0]->fft->data->data[j].im * savBLKData[0]->fft->data->data[j].im);*/

	    fprintf(outfp,"%g %g\n",savBLKData[0]->fft->f0+j*savBLKData[0]->fft->deltaF, 
		outputSFTs->data[0].data->data[j].re * outputSFTs->data[0].data->data[j].re + outputSFTs->data[0].data->data[j].im * outputSFTs->data[0].data->data[j].im);		  
    }

    fclose(outfp);
 
    #endif

    /*}*/
      
    
    if (pPulsarSignalParams->orbit != NULL) {
                 #ifdef DEBUG_LALGENERATEPULSARSFTS
		  fprintf(stdout,"inj pars are omega %f, Tperi %i, ecc %f, angspeed %f, rPeriNorm %f, phi0 %f\n", pPulsarSignalParams->orbit->omega,pPulsarSignalParams->orbit->orbitEpoch.gpsSeconds , pPulsarSignalParams->orbit->oneMinusEcc,  pPulsarSignalParams->orbit->angularSpeed , pPulsarSignalParams->orbit->rPeriNorm,  pPulsarSignalParams->pulsar.phi0 );
                #endif
                  }


    
	} /* END if ( (params->testFlag & 4) > 0 ) else */

    
/* ~~~~ Search will be done using nearest template or the nearest neighbours? ~~~ */
    

if(searchSurroundingOrbitalPts){
fprintf(stdout,"option still not implemented\n");
} else{
fprintf(stdout,"will search only for the nearest neighbour in orbital par space\n");

/*05/09/13 vir: open the orbital parameter file and estimate which is the closest template*/

FILE *fpOrbParamsTemplates;

CHAR filename[256];
strcpy(filename, params->parameterSpaceFile);
fpOrbParamsTemplates = fopen(filename,"r");

if (fpOrbParamsTemplates == NULL){
fprintf(stdout,"Could not find file %s\n", filename);
fflush(stdout);
}

float col1 = 0.0;
UINT4 col2 = 0 ;

INT4 numParams = 0 ;


  while( fscanf(fpOrbParamsTemplates, "%g %i", &col1, &col2)!= EOF){
   numParams = numParams +1;
    }
   printf("numPars %i\n", numParams);
   fclose(fpOrbParamsTemplates);
 

   REAL8 *TemplateSMA = NULL;
   UINT4 *TemplateTperi = NULL;

   TemplateSMA = (REAL8 *)LALMalloc(numParams*sizeof(REAL8));
   TemplateTperi = (UINT4 *)LALMalloc(numParams*sizeof(UINT4));

	      fpOrbParamsTemplates = fopen("OrbParamsTemplates.txt","r");
                  INT4 k=0;
      while( (fscanf(fpOrbParamsTemplates, "%g %i", &col1, &col2)) != EOF){
          
	        TemplateSMA[k]=col1;
		TemplateTperi[k]=col2;    
			k++;
              
              }
        fclose(fpOrbParamsTemplates);	

	REAL8 minDiff = 1.0E5 ;
	REAL8 Diff = 0.0;
    
	for(k=0; k<numParams; k++){
  
		Diff=fabs(TemplateSMA[k] - SemiMajorAxis);


		if(Diff < minDiff){	
		ClosestSMATemplate=TemplateSMA[k];
		minDiff=Diff;
			}
			
           }
	INT4 minDiffT = 10000;
	INT4 DiffT = 0;
        for(k=0; k<numParams; k++){
         
		 DiffT=abs(TemplateTperi[k] - TperiLIGO.gpsSeconds);
		
		if( DiffT < minDiffT ){                            
        	      ClosestTperiTemplate=TemplateTperi[k];
		      minDiffT=DiffT;
					}
			
           }

	fprintf(stdout,"ClosestSMATemplate %f InjectionSMA  %f ClosestTperiTemplate %i InjectionTperi %i\n",ClosestSMATemplate, SemiMajorAxis, ClosestTperiTemplate, TperiLIGO.gpsSeconds);

/*Store That value for the search*/

params->SMAcentral = ClosestSMATemplate;
params->TperiapseSSBSec = ClosestTperiTemplate;
/*remember: when running MonteCarlo, deltaSMA must be set = 0*/



LALFree(TemplateSMA);
LALFree(TemplateTperi);

	
	}/*end of if(nearestneighbours) else...*/

/*...end if (SearchSurroundingPoints)*/
    

/*!*/
/*05/10/18 vir: commented out*/
/**/
    for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) {

      
      if ( (kSUM > 0)  || (iRescale > 0) ) {
       params->startSUMs = 0;  
      }
      
      if ( (iRescale % 2) == 1 ) {
        h_0Rescaled[iRescale] = h_0 + ((REAL8)((iRescale + 1)/2))*params->rescaleMCFraction*h_0;
      } else {
        h_0Rescaled[iRescale] = h_0 - ((REAL8)(iRescale/2))*params->rescaleMCFraction*h_0;
      }      
      
      if(params->numMCRescalings > 0){
      params->threshold4 = h_0Rescaled[iRescale];
      }
         /*printf("h_0Rescaled %g\n", h_0Rescaled[iRescale]);*/
      
              
      if (iRescale > 0) {
         
         REAL4 rescaleFactor = ((REAL4)(h_0Rescaled[iRescale]/h_0));
         if ( !( (params->testFlag & 4) > 0 ) ) {
           
		 renorm = ((REAL4)nSamples)/((REAL4)(outputSFTs->data[0].data->length - 1)); 
		 rescaleFactor = rescaleFactor*renorm;
         } 
	 
	 for(i=0;i<params->numBLKs;i++) {
           for(j=0;j<params->nBinsPerBLK;j++) {
                jSavedBLK = j + jStart; 
		params->BLKData[i]->fft->data->data[j].re = savBLKData[i]->fft->data->data[jSavedBLK].re + rescaleFactor*outputSFTs->data[i].data->data[j].re;
                params->BLKData[i]->fft->data->data[j].im = savBLKData[i]->fft->data->data[jSavedBLK].im + rescaleFactor*outputSFTs->data[i].data->data[j].im;
           }
         }
      } /* END if (iRescale > 0) */

 /*}end of for iRescale < onePlusnumMCRescalings */

/*!*/


    /**************************************************************/
    /* Call StackSlideApplySearch to analyze this MC injection!!! */
    /**************************************************************/


      params->maxPower = 0.0; /* 05/24/05 gam; need to reinitialize */
      StackSlideApplySearch(status->statusPtr,params);
      CHECKSTATUSPTR (status);
  
     
      
                /*05/10/03 may go into a function*/  
         fpReadFromMCInjectionFile=fopen("outLE.txt", "r");
	 fscanf(fpReadFromMCInjectionFile,"%lf%f%lf\n", &StartFreq, &LEvent, &OrbRadius );
	 fprintf(fpMCInjections, "%lf %f %lf %lf %le\n", StartFreq, LEvent, pPulsarSignalParams->pulsar.f0, SemiMajorAxis, h_0);
    
	 /*06/01/09 vir: this should go up, before this loop starts */
         if (params->numMCRescalings){  
	 getStackSlideBinarySearchResults(status->statusPtr, params, &SearchLoudestEvent);
         CHECKSTATUSPTR (status);
         }

    	     if ( (reportMCResults) && (params->numMCRescalings) > 0) {
             nTrials[iRescale] += 1.0;  
	       if (LEvent >= SearchLoudestEvent) {
                 nAboveLE[iRescale] += 1.0;           
	        }
             	/* fprintf(fpMonteCarloResults,"%g %f %f %f %f\n", h_0Rescaled[iRescale], params->f0SUM, SearchLoudestEvent, LEvent, nAboveLE[iRescale]/nTrials[iRescale] );*/
    
	     }
    fprintf(stdout,"Test3 \n\n");

/*else if(reportMCResults)
                  {
      		ComputeConfidence(status->statusPtr, SearchLoudestEvent, &tmp_Confidence, &tmp_conf_err);
                CHECKSTATUSPTR (status);
             fprintf(fpMonteCarloResults,"%g %f %f %f %f\n", h0data[countMC], params->f0SUM, SearchLoudestEvent, tmp_Confidence, tmp_conf_err);

		     }*/

          if (params->numMCRescalings){  
          printf("confidence at rescaling %g is %f\n", h_0Rescaled[iRescale], nAboveLE[iRescale]/nTrials[iRescale]);
	  }
	    }   /* END for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) */ /* 09/01/05 gam */

      
	    if ( !((params->testFlag & 4) > 0) ) {
        LALDestroySFTVector(status->statusPtr, &outputSFTs);
        CHECKSTATUSPTR (status); /* 06/01/04 gam */
      
        LALFree(signal->data->data);
        LALFree(signal->data);
        LALFree(signal);
    }

      printf("kSUM %i \n\n", kSUM);

  } /*end of for kSUM = 0 to numInjectionsTotal*/

/*05/10/03 vir*/
  fclose(fpMCInjections);
  fclose(fpReadFromMCInjectionFile); 


  /**********************************************/
  /*                                            */
  /* END SECTION: MONTE CARLO INJECTIONS LOOP   */
  /*                                            */
  /**********************************************/

      /*StackSlideApplySearch will generate 2 files : one containg the peaks above threshold for this injection and one 
       containing the LE for this injection*/
      /* 1) open previous file containg the loudest event from the search */
      /* 2) Compare each LE in MC with previousLE and find confidence (it will be one for each injection) */
      /* 3) set a desired confidence Conf=0.95 */
      /* 4) compare each confidence with the desired one and choose the closest points to spot to which signal Conf corresponds*/
      
      if ((reportMCResults) && (params->numMCRescalings) == 0){ /*the first one is true for testFlag & 16 > 0 !!*/
      
	      breakOutOfLoop = 0;
       
     /* getStackSlideBinaryPriorResults(status->statusPtr,
                               &priorLoudestEvent,
                               &priorStartFreq,
                               &priorBand,
                               params->priorResultsFile);
      CHECKSTATUSPTR (status);
      */
                    
      getStackSlideBinarySearchResults(status->statusPtr, params, &SearchLoudestEvent);
      CHECKSTATUSPTR (status);

                /*fprintf(stdout,"search loudest event %lf\n", SearchLoudestEvent);
                  fflush(stdout);*/

      		 	  
      ComputeConfidence(status->statusPtr, SearchLoudestEvent, &tmp_Confidence, &tmp_conf_err);
      CHECKSTATUSPTR (status);

                   /*copy final results into a file*/
      
                  fprintf(fpMonteCarloResults,"%g %f %f %f %f\n", h0data[iMC], params->f0SUM, SearchLoudestEvent, tmp_Confidence, tmp_conf_err);
		  
/*! 06/10/01 !*/          }

/*! 06/10/01! */ if(reportMCResults){
		  
		  
	          for(iRescale = 0; iRescale < onePlusNumMCRescalings; iRescale++){
		  arrayULs[countMC] = h_0Rescaled[iRescale];

		  if (nTrials[iRescale] > 0.0) {
                      arrayConfs[countMC] = nAboveLE[iRescale]/nTrials[iRescale];
                      thisMCErr = fabs(Conf - arrayConfs[countMC]);
		  }
		  
	printf("*** arrayULs[countMC] %g, h_0Rescaled[iRescale] %g arrayConfs[countMC] %f\n\n ",  arrayULs[countMC], h_0Rescaled[iRescale], arrayConfs[countMC] );

	if (nTrials[iRescale] > 0.0) {

	fprintf(fpMonteCarloResults,"%g %f %f %f %f\n", h_0Rescaled[iRescale], params->f0SUM, SearchLoudestEvent, LEvent, nAboveLE[iRescale]/nTrials[iRescale] );
	
	}
	
	/*arrayConfs[countMC] = tmp_Confidence;  this contains the confidences for all h_0 
                  thisMCErr = fabs(Conf - arrayConfs[countMC]);
                  
		    
		   fprintf(stdout, "countMC %i thisMCErr %f tmp_Confidence %f \n", countMC, thisMCErr, tmp_Confidence);
		   fflush(stdout);*/
		
		   /*check: if converged ok, else run another set of MCInjections*/
   
                if ( thisMCErr <= maxMCErr ) {
                      arrayConverged[countMC] = 1; 
		   breakOutOfLoop = 1; 

		  	 /*if ( thisMCErr < lastMCErr ) {
                           lastMCErr = thisMCErr;
                           minh_0=h0data[countMC]; 
         	                 arrayULs[countMC]=minh_0;

			   printf("last MCErr %f minh_0 %g\n", lastMCErr , minh_0);
			 }*/
		}else {
                    arrayConverged[countMC] = 0; 
		    arrayULs[countMC] = 0;
       
	       }       
	 
			 
         /* Check, are we done? if none of the injections has produced a good confidence break loop*/
                  
				
         /*if ( thisMCErr <= maxMCErr ) {
            arrayConverged[countMC] = 1;  
            breakOutOfLoop = 1; 
                    
	    
                   arrayULs[countMC]=h_0;
		   printf("h0 %g\n", arrayULs[countMC]);
		   printf("Converged in confidence at first set of iteration...now rescale\n");
         } else {
            arrayConverged[countMC] = 0; 
            arrayULs[countMC] = 0;
	    printf("not converged\n");
           
	 } */      
        
      countMC++; /* Initialized to 0 above; count how many MC simulations have been completed */
      }/* END for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) */

      if (breakOutOfLoop) {
         break; /* MC converged or not producing any trials so break out of iMC loop. */
     }

      if (params->numMCRescalings > 0) {

      /*ComputeUpperLimit(status->statusPtr);
      CHECKSTATUSPTR (status);*/
      }
     /*recombine the sums and look for triggering event...*/
   

      
      }/*end of if (reportMCResults)*/


  }/* end of for iMC =0 to maxMC*/

fclose(fpMonteCarloResults);

/*if( minh_0 == 0 ){
printf("not converged: confidences are too low. Must run another set of injections\n\n");
}else{
	h_0 = minh_0;
        params->threshold4 = minh_0;

	}*/

 /***********************************************************/
 /*                                                         */
 /* END SECTION: LOOP OVER ENTIRE MONTE CARLO SIMULATIONS   */
 /*                                                         */
 /***********************************************************/

CHECKSTATUSPTR(status);
DETATCHSTATUSPTR(status);
}

/**********************************************************************************/
/*         End Function RunStackSlideBinaryMonteCarloSimulation                   */
/**********************************************************************************/




/***********************************************************************************/
/*             Start function getStackSlidePriorResults                            */
/***********************************************************************************/


void getStackSlideBinaryPriorResults(LALStatus *status,
                               REAL4 *priorLoudestEvent,
                               REAL8 *priorStartFreq,
                               REAL8 *priorBand,
                               CHAR  *priorResultsFile)
{
  /*INT4 i;
  INT4 numWordsToIgnore = 43;
  CHAR buffer[256];*/
	
  
  FILE *fpPrior = NULL;
  
  INITSTATUS( status, "getStackSlideBinaryPriorResults", STACKSLIDEBINARYC );
  ATTATCHSTATUSPTR (status);
     
  printf("get prior results\n");
  
  
  fpPrior = fopen( priorResultsFile, "r");
  fscanf(fpPrior, "%lf%f%lf\n",priorStartFreq, priorLoudestEvent,  priorBand);
    
  fclose(fpPrior);
  
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);

  
}


/***********************************************************************************/
/*           Start function Compute Confidence                                     */
/***********************************************************************************/


void ComputeConfidence(LALStatus *status, REAL4 priorLoudestEvent, REAL8 *Confidence, REAL8 *conf_err){
  INITSTATUS( status, "ComputeConfidence", STACKSLIDEBINARYC );
  ATTATCHSTATUSPTR (status);


  REAL8 freq;
  REAL4 LE;
  REAL8 OrbRadius;
  REAL8 a;
  REAL8 b;
  
  INT4 count = 0;
  INT4 tot = 0;
 
 
  
  /*there should be an option to open many files containes in a result directory*/
  FILE *fpMCInjectionResults;
  
  fprintf(stdout, "Compute confidence\n");
  fflush(stdout);

 

  fpMCInjectionResults = fopen("MCInjectionResults.txt", "r");
  
  if (fpMCInjectionResults == NULL){
  fprintf(stdout, "cannot open file MCInjectionResults.txt \n\n" );
  fflush(stdout);
  }
	  
  /*LE has to be the loudest of all */
    while( fscanf(fpMCInjectionResults, "%lf%f%lf%lf%lf\n", &freq, &LE, &OrbRadius, &a, &b) != EOF){
  
	 	  if (LE >= priorLoudestEvent) count++;
           	  tot++;

                #ifdef DEBUG_CONFIDENCE
		  printf("MCLE %f tot %i count %i SearchLE %f\n", LE, tot, count, priorLoudestEvent);
                #endif

  }

  fclose(fpMCInjectionResults);
  
  *Confidence=(REAL8)((REAL8)count/(REAL8)tot);
  *conf_err=(REAL8)(1.0/sqrt((REAL8)tot));

  fprintf(stdout,"end function ComputeConfidence \n\n");
  fflush(stdout);

  
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);
}

/***********************************************************************************/
/*             End function Compute Confidence                                     */
/***********************************************************************************/


/***********************************************************************************/
/*             Start function getStackSlideBinarySearchResults                     */
/***********************************************************************************/


void getStackSlideBinarySearchResults(LALStatus *status, StackSlideSearchParams *params, REAL8 *SearchLoudestEvent){
 
	INITSTATUS( status, "getStackSlideBinarySearchResults", STACKSLIDEBINARYC );
        ATTATCHSTATUSPTR (status);
        
       #define TEST_MODE
	
       /* the purpose of this function is to open a previous file in the pipeline containg the results af the search on */
       /* that particular frequency band and for each of the filters in the grid*/
	
        FILE *fpSearch;
	REAL8 a;
        REAL8 b;
	REAL8 c;
	REAL8 MinLoudest = 0.0 ;
	
	fprintf(stdout,"start function getStackSlideBinarySearchResults\n ");
        fflush(stdout);

     /*	printf("params->NormalizationParam %f\n", params->normalizationParameter);*/
        		
        char filename[512];
     /*sprintf(filename,"PriorResultsFile_%.3f-%.2f-%i.txt", params->f0SUM, params->SMAcentral, params->TperiapseSSBSec);*/
        
	sprintf(filename,"PriorResultsFile_%.3f.txt", params->f0SUM);
        
        #ifdef DEBUG_GETSTACKSLIDEBINARYSEARCHRESULTS
	fprintf(stdout,"filename %s\n", filename);
	fflush(stdout);
        #endif
        
	fpSearch=fopen(filename,"r");
	if ( fpSearch == NULL ) 
	{
		fprintf(stdout,"The search file %s doesn't exist!!\n\n", filename);
	        fflush(stdout);
        #ifdef TEST_MODE
        fprintf(stdout, "WARNING: we are opening an archive file for compiling purposes only!!\n\n");
	fflush(stdout);
	
		fpSearch=fopen("PriorResultsFile_599.000-1.26-731163327.txt","r");
        
        #endif	
	}
	
	
	
	/*fscanf(fpSearch, "%lf %lf %lf", &a, SearchLoudestEvent, &c) ;*/
	 
        /* determine the loudest of all here and call THAT SearchLoudestEvent....do that when the condor script is set*/
	  
	while ( (fscanf(fpSearch, "%lf %lf %lf", &a, &b, &c)) != EOF ){
	 
		if (b >= MinLoudest){
		MinLoudest = b;
		}
	
	}
  
	*SearchLoudestEvent = MinLoudest;
	printf("SearchLoudestEvent %f \n\n", *SearchLoudestEvent);

	fprintf(stdout,"end function getStackSlideBinarySearchResults\n ");
        fflush(stdout);


	fclose(fpSearch);
	 
	

	CHECKSTATUSPTR (status);
        DETATCHSTATUSPTR (status);

}

void ValidateMCResults(LALStatus *status, const REAL4FrequencySeries *oneSUM, StackSlideSearchParams *params, REAL4 *SNR){
	INITSTATUS( status, "ValidateMCResults", STACKSLIDEBINARYC );
        ATTATCHSTATUSPTR (status);

/*        #define PRINT_VALIDATEMC */
 
	fprintf(stdout,"start function ValidateMCResults\n ");
        fflush(stdout);

	
	INT4 i=0;
        INT4 iPwrMax = 0;
	REAL4 pwrStdDev = 0.0;
	REAL4 pwrMax = 0.0;
        REAL4 pwrSNR = 0.0;
	
	REAL4 pwrSum= 0.0;    /* sum of power */
        REAL4 pwrSum2 = 0.0;  /* sum of power squared = sum of power*/
        REAL4 pwrMean = 0.0;
      
        

	for(i=0;i<oneSUM->data->length;i++) {

       	pwrSum += oneSUM->data->data[i];
        pwrSum2 += oneSUM->data->data[i]*oneSUM->data->data[i];
	}

        for(i=0; i<oneSUM->data->length; i++){
        	if (oneSUM->data->data[i] > pwrMax) {
                     pwrMax = oneSUM->data->data[i];
                     iPwrMax = i;
                  }
	} 
        
	pwrSNR = pwrMax*sqrt(params->numBLKs);

	printf("pwrMax %f iPwrMax %d pwrSNR %f \n", pwrMax, iPwrMax, pwrSNR);

	pwrMean = pwrSum/(oneSUM->data->length); /* Mean pwr */
        pwrStdDev = pwrSum2/(oneSUM->data->length) - pwrMean*pwrMean;   

	if ( pwrStdDev > 0.0 ) {
             pwrStdDev = sqrt( pwrStdDev );
           } else {
            pwrStdDev = -1.0; /* greg: Should not happen except due to round off error in test data */
           }

        
	pwrMax=pwrMax/pwrStdDev;
        
	*SNR = pwrMax;
    
        #ifdef PRINT_VALIDATEMC
	fprintf(stdout,"length %d pwrStdDev %f\n", oneSUM->data->length, pwrStdDev);
        fflush(stdout);
        #endif
	
	CHECKSTATUSPTR (status);
        DETATCHSTATUSPTR (status);

	
}

void ComputeUpperLimit(LALStatus *status, const REAL8 *arrayULs, const REAL8 *arrayConfs, REAL8 desiredConf){
	INITSTATUS( status, "ComputeUL", STACKSLIDEBINARYC );
        ATTATCHSTATUSPTR (status);
        
        #define DEBUG_COMPUTEUPPERLIMIT
	
        #ifdef DEBUG_COMPUTEUPPERLIMIT
        printf("Start function ComputeUL\n");
        #endif

 
	
        #ifdef DEBUG_COMPUTEUPPERLIMIT
	printf("End function ComputeUL\n");
        #endif
        CHECKSTATUSPTR (status);
        DETATCHSTATUSPTR (status);


}


void Readh0File(LALStatus *status, FILE *fph0, INT4 *N, REAL8 **h0){

        INITSTATUS( status, "Readh0File", STACKSLIDEBINARYC );
        ATTATCHSTATUSPTR (status);
	
      	printf("Readh0File\n");
        INT4 i=0;
	REAL8 dummy = 0.0;
        char filename[]="h0File.txt";
	fph0=fopen(filename,"r");

	 if (fph0==NULL) {
         printf("ERROR : cannot open file %s\n",filename);
         exit(1);
        }
      
      
	 while(fscanf( fph0,"%lf\n",&dummy)!=EOF) {
		 i++;
		 
	 }
        	 fclose(fph0);
               
         (*N)=i;

        (*h0)=(REAL8 *)LALMalloc((i)*sizeof(REAL8));

         i=0;
         fph0=fopen(filename,"r");
         while (fscanf(fph0,"%lf\n",&dummy)!=EOF) {
         (*h0)[i]=dummy;
         i++;
	     }

        fclose(fph0);


	 
        CHECKSTATUSPTR (status);
        DETATCHSTATUSPTR (status);

}
