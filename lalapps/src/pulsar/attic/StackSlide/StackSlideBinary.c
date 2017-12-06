/*
*  Copyright (C) 2007 Gregory Mendell, Virginia Re
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
 * \author Virginia Re
 */

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
 /*#define PRINT_RES*/
 #define PRINT_INJECTION_VALUES_MC
/* #define DEBUG_LALFASTGENERATEPULSARSFTS*/
 #define DEBUG_LALGENERATEPULSARSFTS
/*#define DEBUG_CONFIDENCE*/
/* #define PRINT_MAXPOWERANDBINEACHSFT*/
#define PRINT_ESTIMATEDSTACKSLIDEPOWER_WARNING

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

    /*char filename[]="myoutbinary.txt";*/
    char filename[512];
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

    /*REAL4 SNR;*/
    
    fprintf(stdout,"Start function StackSlideBinary\n");
    fflush(stdout);
    INITSTATUS(status); 
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
        
	    if (params->numMCInjections == 0) { /*06/04/07 vir*/	
	       binaryLE=fopen(filename2, "w");
               fpSavedEvents=fopen("SavedEvents.txt","w+");
	       }
	
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
            fflush(stdout);
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

              /*06/04/24 vir*/
	       if ((params-> debugOptionFlag & 64) >0 ){
	       #ifdef PRINT_ESTIMATEDSTACKSLIDEPOWER_WARNING
               fprintf(stdout,"\n\n !!! WARNING: FINDING ESTIMATED STACKSLIDE POWER !!! \n\n");
               fflush(stdout);
               #endif    
		   ValidateMCResults(status->statusPtr,params->SUMData,params->STKData,pTdotsAndDeltaTs,params,iSky,stksldParams);    
	       }else{
	       
		       StackSlide(status->statusPtr,params->SUMData,params->STKData,pTdotsAndDeltaTs,stksldParams);
	       }
	       
	      CHECKSTATUSPTR (status);
                
	       
               ParamsSMA[kSUM]=SemiMajorAxis; /*05/02/18 vir:*/ /* 04/12/05 gam; moved out of StackSlide function */
                
                 
               if((params->nMaxSMA ==1)&&(params->nMaxTperi ==1)) {

        if (params->numMCInjections < 2) { /*06/04/03 vir*/		    
  
		                strcpy(filename, "myoutbinary.txt");

				binaryfp=fopen(filename, "w");
		  
                  for (k=0; k< params->nBinsPerSUM ;k++) {
			  
  	            fprintf(binaryfp,"%f %g\n", params->f0SUM + (REAL8)k*params->dfSUM, params->SUMData[0]->data->data[k]);
  /*....maybe when in MC mode you don't need to open so many files because FindbinaryLoudest looks for the max in the data*/
		  }
		  /*fflush(binaryfp);*/
                  fclose(binaryfp);
  }  /*06/04/03 vir*/	  /*end if (params->numMCInjections == 0)*/                

	
		
	/*Look for loudest also in this case ....*/

		  FindBinaryLoudest(&LoudestEvent, &peakFreq, params->SUMData, stksldParams, SemiMajorAxis, fpSavedEvents);                   /*06/03/14 vir*/
		 params->maxPower=LoudestEvent; 
		 /**/
		  		   
                  if (params->numMCInjections == 0) { /*06/04/03 vir*/	
		  fprintf(binaryLE,"%f %f %f\n", peakFreq, LoudestEvent, ParamsSMA[kSUM]);
                   }/*06/03/04 vir*/
		      
        	 /*05/11/17 add calculation of standard deviation only if MC active*/
                   /*if( (params->testFlag & 2) > 0 ){  
	           ValidateMCResults(status->statusPtr, params->SUMData[0], params, &SNR);
                   CHECKSTATUSPTR (status);
                   		   
		   } */

		   
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

		    /*06/03/27 vir: added this option to out myoutbinary file for each filter */
		    if (params->numMCInjections == 0 ){
		    
			    sprintf(filename, "myoutbinary_f%f-R%f.txt", params->f0SUM, SemiMajorAxis);

		            binaryfp=fopen(filename, "w");
		  
                  for (k=0; k< params->nBinsPerSUM ;k++) {
			  
  	            fprintf(binaryfp,"%f %g\n", params->f0SUM + (REAL8)k*params->dfSUM, params->SUMData[0]->data->data[k]);
  		
		        }
                     fclose(binaryfp);
		    }
		   
	       }/*end of else nMaxSMA == 1 */
          
	       /*Allocate here space for xml table*/
	       
	    } /* for iFreqDeriv = 0 to numFreqDerivTotal */

            LALDestroyRandomParams(status->statusPtr, &randPar );
   
	   
         /*}end of for iT*/
        }/*end of for iSMA 05/02/17 vir*/
     
      fprintf(stdout,"params->numMCInjections %d\n",params->numMCInjections );
	           fflush(stdout);

if (params->numMCInjections == 0) { /*06/04/07 vir*/	
       fclose(binaryLE);    
       fclose(fpSavedEvents);
}
             
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
       
	fprintf(stdout,"Loudest Event from sum (or MC) %f\n", max);
	fflush(stdout);
	
 for (i=iMinSTK; i<iMaxSTK; i++)
        {
                if(SUMData[0]->data->data[i] > threshold)
                { 
		  peak = SUMData[0]->data->data[i];
                  indexFreq=i;
                
/*                *peakFreq=stksldParams->f0SUM + indexFreq*stksldParams->dfSUM;*/
		  
             /*fprintf(fpSavedEvents,"%f %f %f \n",stksldParams->f0SUM + indexFreq*stksldParams->dfSUM , peak, SemiMajorAxis);  05/09/06 vir*/     
	
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

  INT4 Vcount ;
  INT4 Vtot ;
  
 /* FILE *outfp;*/
  
  INT4 iFreq = 0;                           

  REAL8 cosIota;                    
  REAL8 h_0 = 0.0;  /* params->threshold4;    Source amplitude: not sure if threshold 4*/
  
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

  BOOLEAN BinsCanBeExcluded = 1;
  BOOLEAN searchSurroundingPts = 0;
  BOOLEAN searchSurroundingOrbitalPts = 0;
  /*use results from a previous job in the pipeline*/
  BOOLEAN reportMCResults = 0;
 /* REAL4 priorLoudestEvent = 0.0;
  REAL8 priorStartFreq = 0.0;
  REAL8 priorBand = 0.0;*/
  /*REAL8 tmp_Confidence = 0.0;*/
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
  REAL8 AmplSMAParSpace = 0.36;         /* 05/09/12 vir :in sec ; 06/03/15 vir : modified, shold be larger*/
  REAL8 AmplTperiParSpace = 600;        /* 05/09/12 vir :in sec */
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
  FILE *fpRandom;
  INT4 rndCount;

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
  
  INT4 iMC = 0;             /* which entire MC simulation is currently running */
  INT4 countMC = 0;         /* count of number of completed MC simulations, including rescalings */
  INT4 maxMC = 1;  /* default number of MC iterations to perform; set equal to params->maxMCinterations when interating MC */
  
 
        FILE *fpOrbParamsTemplates;
        /*FILE *fpMCInjections; modify this into a name containing which injection we are working on*/
        /*FILE *fpReadFromMCInjectionFile;*/
        /*REAL8 StartFreq;*/
        /*REAL4 LEvent;*/
	/*REAL8 OrbRadius;*/
        REAL8 Conf = 0.95;        

	/*store final results*/
        FILE *fpMonteCarloResults = NULL; /*used only when thresholdFlag==0 and numRescalings>0*/
	CHAR MonteCarloResults[256];
	
  
   INITSTATUS(status);
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
    fprintf(stdout,"MaxMC %i\n", maxMC);
    fflush(stdout);
   }

if(reportMCResults){/*06/06/12 vir*/

if (params->numMCRescalings > 0) {
      arrayULs   = (REAL8 *)LALMalloc((twoPlusNumMCRescalings*maxMC)*sizeof(REAL8));
      arrayConfs = (REAL8 *)LALMalloc((twoPlusNumMCRescalings*maxMC)*sizeof(REAL8));
      arrayConverged = (INT4 *)LALMalloc((twoPlusNumMCRescalings*maxMC)*sizeof(INT4));
    } else {

      arrayULs   = (REAL8 *)LALMalloc(maxMC*sizeof(REAL8));
      arrayConfs = (REAL8 *)LALMalloc(maxMC*sizeof(REAL8));
      arrayConverged = (INT4 *)LALMalloc(maxMC*sizeof(INT4));
    }
}/*06/06/12 vir*/

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
  pPulsarSignalParams->dtDelayBy2 = 0;
  pPulsarSignalParams->dtPosBy2 = 0;
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


/* 06/06/13 vir */
  if(BinsCanBeExcluded){
   if (searchSurroundingPts) {
      params->bandSUM = 2.0*params->dfSUM;
      params->nBinsPerSUM = 2;              /* Search surrounding frequencies */
   } else {
      params->bandSUM = params->dfSUM;      /* Search nearest frequencies */
      params->nBinsPerSUM = 1;
   }
  }
  

   /*generate random numeber seed here */
  fpRandom = fopen("/dev/urandom","r");
  rndCount = fread(&seed, sizeof(INT4),1, fpRandom);
  fclose(fpRandom);

   LALCreateRandomParams(status->statusPtr, &randPar, seed); /*06/03/20 vir*/

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

/*06/03/07 vir: open file containing orbital parameters and store them into array */

/*FILE *fpOrbParamsTemplates;*/

CHAR filename[256];
strcpy(filename, params->parameterSpaceFile);
fpOrbParamsTemplates = fopen(filename,"r");

if (fpOrbParamsTemplates == NULL){
fprintf(stdout,"Could not find file %s\n", filename);
fflush(stdout);
}

REAL8 col1 = 0; /*06/03/24 vir: float to REAL8*/
UINT4 col2 = 0 ;

INT4 numParams = 0 ;


  while( fscanf(fpOrbParamsTemplates, "%lf %i", &col1, &col2)!= EOF){
   numParams = numParams +1;
    }
   fprintf(stdout,"numPars %i\n", numParams);
   fflush(stdout);
   fclose(fpOrbParamsTemplates);
 

   REAL8 *TemplateSMA = NULL;
   UINT4 *TemplateTperi = NULL;

   TemplateSMA = (REAL8 *)LALMalloc(numParams*sizeof(REAL8));
   TemplateTperi = (UINT4 *)LALMalloc(numParams*sizeof(UINT4));

	      fpOrbParamsTemplates = fopen("OrbParamsTemplates.txt","r");
                  INT4 z=0;
      while( (fscanf(fpOrbParamsTemplates, "%lf %i", &col1, &col2)) != EOF){
          
	        TemplateSMA[z]=col1;
		TemplateTperi[z]=col2;    
			z++;
              
              }
        fclose(fpOrbParamsTemplates);	

/*end open orbital parameter file */

/*put here getStackSlideBinarySearchResults*/
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
/*    fpMCInjections=fopen("MCInjectionResults.txt", "w+"); 06/4/24 removed because not used*/

   
   if (iMC > Nh0) break; 
   
   if ( params->threshold4 == 0 )
   {
	    h_0 = h0data[iMC];
	    fprintf(stdout," h_0 = %g\n\n", h0data[iMC]);
	    fflush(stdout);
      }else{
     h_0=params->threshold4 ;
      }
 
   /*  printf("h_0 %g\n", h_0);*/
      
  /**********************************************/
  /*                                            */
  /* START SECTION: MONTE CARLO INJECTIONS LOOP */
  /*                                            */
  /**********************************************/
     Vcount=0; /*06/03/15 vir: need to initialize here */
     Vtot=0;
   
    for(kSUM=0;kSUM<params->numMCInjections;kSUM++) {
     params->SMAcentral=1.44; /*06/03/20 vir: need to reinitialize here*/
	  
    params->whichMCSUM = kSUM; /* kSUM is which injection we are working on */  

    if (MCInjectRandomSMA) { /*then generate a random SMA for injection*/
      /*LALCreateRandomParams(status->statusPtr, &randPar, seed); */ /*06/03/20 vir*/
      LALUniformDeviate(status->statusPtr, &randval, randPar); 
      SemiMajorAxis=params->SMAcentral+((REAL8)randval-0.5)*(AmplSMAParSpace);
#ifdef PRINT_INJECTION_VALUES
      fprintf(stdout,"params->SMAcentral %lf SemimajorAxis for Injection %f randval %f\n\n",params->SMAcentral, SemiMajorAxis, (REAL8)randval);
      fflush(stdout);
#endif
      
    }else{
       /* SemiMajorAxis=params->SMAcentral;  */
       SemiMajorAxis = 1.44;
    }
 

    if (MCInjectRandomTperi ) {
	  /*LALCreateRandomParams(status->statusPtr, &randPar, seed); */ /*06/03/20 vir*/
          LALUniformDeviate(status->statusPtr, &randval, randPar); 
          TperiLIGO.gpsSeconds=params->TperiapseSSBSec+((REAL8)randval-0.5)*(AmplTperiParSpace);
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
     /*LALCreateRandomParams(status->statusPtr, &randPar, seed); */ /*06/03/20 vir*/
     LALUniformDeviate(status->statusPtr, &randval, randPar); 
     real8RandVal=randval*real8NAvailableBins;
     i = floor( real8RandVal );    /* nearest bin arrayAvailableFreqBins */
     if (i >= nAvailableBins) i = nAvailableBins - 1; if (i < 0) i = 0; /* make sure i in range */
     iFreq = arrayAvailableFreqBins[i];  /* randomly chosen available frequency bin */
          params->f0SUM = f0SUM + iFreq*params->dfSUM; /* search f0*/
    fprintf(stdout,"iFreq %i real8RandVal %f, params->f0SUM %f\n", iFreq, real8RandVal, params->f0SUM);  
     fflush(stdout);

    }else{
    
	    /*randomly a bin in the SUM where to inject signal*/
	/*LALCreateRandomParams(status->statusPtr, &randPar, seed); */ /*06/03/20 vir*/

        LALUniformDeviate(status->statusPtr, &randval, randPar); 
        InjFreqBin = floor(randval*(params->nBinsPerSUM));
	if(InjFreqBin >= params->nBinsPerSUM) InjFreqBin = params->nBinsPerSUM - 1;
       
	params->f0SUM=f0SUM ; /*05/10/19 vir: the start of the sum remains the same*/
    }
         /*LALCreateRandomParams(status->statusPtr, &randPar, seed); */ /*06/03/20 vir*/
         LALUniformDeviate(status->statusPtr, &randval, randPar); 
         real8RandVal=params->f0SUM - (0.5*params->dfSUM) + randval*params->dfSUM + InjFreqBin*params->dfSUM;
	 pPulsarSignalParams->pulsar.f0 = real8RandVal; 

     #ifdef PRINT_INJECTION_VALUES_MC
	 fprintf(stdout, "f0SUM %f, pPulsarSignalParams->pulsar.f0 %f InjFreqBin %d, params->bandSUM %f\n", f0SUM, pPulsarSignalParams->pulsar.f0, InjFreqBin, params->bandSUM);
	 fflush(stdout);
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
	 /* LALCreateRandomParams(status->statusPtr, &randPar, seed); */
          LALUniformDeviate(status->statusPtr, &randval, randPar); 
          pPulsarSignalParams->pulsar.psi = 4*(randval - 0.5) * ((REAL4)LAL_PI_2); /*06/03/20 vir: -pi < psi< pi*/
          /*rand cosiota*/
	 /* LALCreateRandomParams(status->statusPtr, &randPar, seed);*/ 
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
fflush(stdout);    
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

	fprintf(stdout,"use LALGeneratePulsarSignal and LALSignalToSFTs\n");
        fprintf(stdout,"GPSin.Seconds %i\n", GPSin.gpsSeconds);
	fprintf(stdout,"pPulsarSignalParams->pulsar.tRef.gpsSeconds %i\n",  pPulsarSignalParams->pulsar.tRef.gpsSeconds);
        fflush(stdout);
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
        fflush(stdout);
        #endif    
        
	outputSFTs = NULL; /* Call LALSignalToSFTs to generate outputSFTs with injection data */
        LALSignalToSFTs(status->statusPtr, &outputSFTs, signal, pSFTParams);
        CHECKSTATUSPTR (status);

       fprintf(stdout,"Test4\n\n");
       fflush(stdout);  


	/*(whichMCsum =0 )*/

              #ifdef PRINT_INJECTION_VALUES_MC
              fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",InjFreqBin,h_0);
              fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",
			      InjFreqBin,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
              fprintf(stdout,"iFreq = %i, inject psi = %23.10e, LAL_PI_2 %g, LAL_TWOPI %g\n",InjFreqBin,pPulsarSignalParams->pulsar.psi, LAL_PI_2, LAL_TWOPI);
              fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",InjFreqBin,pPulsarSignalParams->pulsar.phi0);
              fprintf(stdout,"iFreq = %i, search f0 = %23.10e, inject f0 = %23.10e \n",
			      InjFreqBin,f0SUM + InjFreqBin*params->dfSUM,pPulsarSignalParams->pulsar.f0);
             fprintf(stdout,"outputSFTs->data[0].data->length = %i \n",outputSFTs->data[0].data->length);
             fflush(stdout);       
             #endif

	     /* CHECK the index i*/
             
	     i=0;
	     renorm = ((REAL4)nSamples)/((REAL4)(outputSFTs->data[i].data->length - 1));
             fprintf(stdout,"renorm %g\n", renorm);
	     fflush(stdout);

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
fflush(stdout);

	  /* if only 1 MC sum {*/
  
    #ifdef PRINT_RES   
   FILE *outfp;
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
fflush(stdout);   
#endif
                  }


    
	} /* END if ( (params->testFlag & 4) > 0 ) else */

    
/* ~~~~ Search will be done using nearest template or the nearest neighbours? ~~~ */
    

if(searchSurroundingOrbitalPts){
fprintf(stdout,"option still not implemented\n");
fflush(stdout);
} else{
fprintf(stdout,"will search only for the nearest neighbour in orbital par space\n");
fflush(stdout);
/*05/09/13 vir: open the orbital parameter file and estimate which is the closest template*/
/*moved up*/
/*end moved up*/

	REAL8 minDiff = 1.0E5 ;
	REAL8 Diff = 0.0;
        INT4 k=0;
	for(k=0; k<numParams; k++){
  
		Diff=fabs(TemplateSMA[k] - SemiMajorAxis);
                 #ifdef PRINT_INJECTION_VALUES
		    fprintf(stdout,"templateSMA[%i] %f \n\n",k, TemplateSMA[k] );
                    fflush(stdout);
                 #endif
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

  #ifdef PRINT_INJECTION_VALUES
	fprintf(stdout,"ClosestSMATemplate %f InjectionSMA  %f ClosestTperiTemplate %i InjectionTperi %i\n",ClosestSMATemplate, SemiMajorAxis, ClosestTperiTemplate, TperiLIGO.gpsSeconds);
         fflush(stdout);
  #endif

	 /*Store That value for the search*/

params->SMAcentral = ClosestSMATemplate;
params->TperiapseSSBSec = ClosestTperiTemplate;
/*remember: when running MonteCarlo, deltaSMA must be set = 0*/


	
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
         fprintf(stdout,"h_0Rescaled %g\n", h_0Rescaled[iRescale]);
         fflush(stdout);
              
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

      fprintf(stdout,"Test5\n\n");
      fflush(stdout);
 /*}end of for iRescale < onePlusnumMCRescalings */

/*!*/


    /**************************************************************/
    /* Call StackSlideApplySearch to analyze this MC injection!!! */
    /**************************************************************/
   
      params->maxPower = 0.0; /* 05/24/05 gam; need to reinitialize */
      StackSlideApplySearch(status->statusPtr,params);
      CHECKSTATUSPTR (status);
  
     
          /*05/10/03 may go into a function*/  
        /* fpReadFromMCInjectionFile=fopen("outLE.txt", "r");
	 fscanf(fpReadFromMCInjectionFile,"%lf%f%lf\n", &StartFreq, &LEvent, &OrbRadius );
	 fprintf(fpMCInjections, "%lf %f %lf %lf %le\n", StartFreq, LEvent, pPulsarSignalParams->pulsar.f0, SemiMajorAxis, h_0);   06/03/07 vir: commented out*/
    
	fprintf(stdout,"params->maxPower %lf\n\n", params->maxPower);
        fflush(stdout);
	 /*06/01/09 vir: this should go up, before this loop starts */
         if (params->numMCRescalings){  
	 getStackSlideBinarySearchResults(status->statusPtr, params, &SearchLoudestEvent);
         CHECKSTATUSPTR (status);
         }

    	     if ( (reportMCResults) && (params->numMCRescalings) > 0) {
             nTrials[iRescale] += 1.0;  
	       if (params->maxPower >= SearchLoudestEvent) {
                 nAboveLE[iRescale] += 1.0;           
	        }
             	/* fprintf(fpMonteCarloResults,"%g %f %f %f %f\n", h_0Rescaled[iRescale], params->f0SUM, SearchLoudestEvent, LEvent, nAboveLE[iRescale]/nTrials[iRescale] );*/
    
	     }
  

    /*06/03/14 vir: moved here from below, but must go still up. Also, take into account that if numrescalings >0....*/
      getStackSlideBinarySearchResults(status->statusPtr, params, &SearchLoudestEvent); 
      CHECKSTATUSPTR (status);
      
      /*ComputeConfidence(status->statusPtr, SearchLoudestEvent, params->maxPower, &tmp_Confidence, &tmp_conf_err);
      CHECKSTATUSPTR (status);*/
    
    /**/


      /*06/03/15 vir: recalculate conf*/
    if ((reportMCResults) && (params->numMCRescalings) == 0){ 
      if (params->maxPower > SearchLoudestEvent){
      Vcount++;
      }
      Vtot++;
    }  
            
      /**/
    
    /*else if(reportMCResults)
                  {


      		ComputeConfidence(status->statusPtr, SearchLoudestEvent, &tmp_Confidence, &tmp_conf_err);
                CHECKSTATUSPTR (status);
             fprintf(fpMonteCarloResults,"%g %f %f %f %f\n", h0data[countMC], params->f0SUM, SearchLoudestEvent, tmp_Confidence, tmp_conf_err);

		     }*/

          if (params->numMCRescalings){  
          fprintf(stdout,"confidence at rescaling %g is %f nTrials %g\n", h_0Rescaled[iRescale], nAboveLE[iRescale]/nTrials[iRescale], nTrials[iRescale] );
	  fflush(stdout);
	  }
	    }   /* END for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) */ /* 09/01/05 gam */

      
	    if ( !((params->testFlag & 4) > 0) ) {
        LALDestroySFTVector(status->statusPtr, &outputSFTs);
        CHECKSTATUSPTR (status); /* 06/01/04 gam */
      
        LALFree(signal->data->data);
        LALFree(signal->data);
        LALFree(signal);
    }

      fprintf(stdout,"kSUM %i \n\n", kSUM);
      fflush(stdout);

  } /*end of for kSUM = 0 to numInjectionsTotal*/
 
if ((reportMCResults) && (params->numMCRescalings) == 0){ 
      /*tmp_Confidence = (REAL8)((REAL8)Vcount/(REAL8)Vtot);*/
	arrayConfs[iMC]=(REAL8)((REAL8)Vcount/(REAL8)Vtot);
      tmp_conf_err=(REAL8)(1.0/sqrt((REAL8)Vtot));
      fprintf(stdout,"Vtot %i Vcount %i iMC %i, New tmp_Confidence %lf Conf_err %f\n", Vtot, Vcount, iMC, arrayConfs[iMC], tmp_conf_err);
 }


/*05/10/03 vir*/
 /* fclose(fpMCInjections);
  fclose(fpReadFromMCInjectionFile); 06/03/07 vir: commented out as above*/


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
                    
 
                 if( params->threshold4 == 0){
               
	      fprintf(stdout,"h0data[iMC] %le \n\n", h0data[iMC]);
              fflush(stdout);

	      /*copy final results of MC Injections into a file: this might be called MCInjectionResults.txt from CLA*/

                  fprintf(fpMonteCarloResults,"%g %f %f %f %f\n", h0data[iMC], params->f0SUM, SearchLoudestEvent, arrayConfs[iMC], tmp_conf_err);
      }	  

      /*! 06/10/01 !*/          }

/*! 06/10/01! */ if(reportMCResults){
		  
		  
	          for(iRescale = 0; iRescale < onePlusNumMCRescalings; iRescale++){
		  arrayULs[countMC] = h_0Rescaled[iRescale];

		  if (nTrials[iRescale] > 0.0) {
                      arrayConfs[countMC] = nAboveLE[iRescale]/nTrials[iRescale];
                      thisMCErr = fabs(Conf - arrayConfs[countMC]);
		  
printf("*** arrayULs[countMC] %g, h_0Rescaled[iRescale] %g arrayConfs[countMC] %g\n\n ",  arrayULs[countMC], h_0Rescaled[iRescale], arrayConfs[countMC]  );
		  
/*copy final results of MC Injections into a file: this might be called MCInjectionResults_Rescalings.txt from CLA*/
	fprintf(fpMonteCarloResults,"%g %f %f %g %f\n", h_0Rescaled[iRescale], params->f0SUM, SearchLoudestEvent, params->maxPower, arrayConfs[countMC] );
	
	}
	
		
		   /*check: if converged ok, else run another set of MCInjections*/
   
                if ( thisMCErr <= maxMCErr ) {
                      arrayConverged[countMC] = 1; 
		   breakOutOfLoop = 1;
		   fprintf(stdout,"Converged in confidence at first set of iteration...now rescale\n");
		   fflush(stdout);
		   if (params->threshold4 == 0){
                    arrayULs[countMC]=h0data[countMC];
		   }
		    else if(params->numMCRescalings){
		    arrayULs[countMC]=h_0Rescaled[iRescale];
		    }
		   
		  	 }else {
                    arrayConverged[countMC] = 0; 
		    arrayULs[countMC] = 0;
                    fprintf(stdout,"not converged\n");
		    fflush(stdout);

	       }       
	 
			 
                
      countMC++; /* Initialized to 0 above; count how many MC simulations have been completed */
      }/* END for(iRescale=0;iRescale<onePlusNumMCRescalings;iRescale++) */

      if (breakOutOfLoop) {
         break; /* MC converged or not producing any trials so break out of iMC loop. */
     }

      if (params->numMCRescalings > 0) {

      /*ComputeUpperLimit(status->statusPtr);
      CHECKSTATUSPTR (status);*/
      }
     
      
      }/*end of if (reportMCResults)*/


  }/* end of for iMC =0 to maxMC*/

if(reportMCResults){
    LALFree(arrayULs);
    LALFree(arrayConfs);
    LALFree(arrayConverged);
    LALFree(nAboveLE);   
    LALFree(nTrials); 

}/*06/06/12 vir*/

LALFree(h_0Rescaled); 

  LALDestroyRandomParams(status->statusPtr, &randPar); /* 06/06/12 vir */
  CHECKSTATUSPTR (status);

LALFree(pSFTParams);

LALFree(pPulsarSignalParams->orbit);
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

/* 06/06/12 free memory of saved data */
  for (i=0;i<params->numBLKs;i++) {
      LALFree(savBLKData[i]->fft->data->data);
      LALFree(savBLKData[i]->fft->data);
      LALFree(savBLKData[i]->fft);
      LALFree(savBLKData[i]);
  }
  LALFree(savBLKData);

  LALFree(savSumBinMask);          /* 06/06/12 vir */
  LALFree(arrayAvailableFreqBins); 


  
LALFree(TemplateSMA);
LALFree(TemplateTperi);

fclose(fpMonteCarloResults);


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


/*void getStackSlideBinaryPriorResults(LALStatus *status,
                               REAL4 *priorLoudestEvent,
                               REAL8 *priorStartFreq,
                               REAL8 *priorBand,
                               CHAR  *priorResultsFile)
{*/
  /*INT4 i;
  INT4 numWordsToIgnore = 43;
  CHAR buffer[256];*/
	
 /* 
  FILE *fpPrior = NULL;
  
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
     
  printf("get prior results\n");
  
  
  fpPrior = fopen( priorResultsFile, "r");
  fscanf(fpPrior, "%lf%f%lf\n",priorStartFreq, priorLoudestEvent,  priorBand);
    
  fclose(fpPrior);
  
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);

  
}*/


/***********************************************************************************/
/*           Start function Compute Confidence                                     */
/***********************************************************************************/


void ComputeConfidence(LALStatus *status, REAL4 priorLoudestEvent, REAL4 maxPower, REAL8 *Confidence, REAL8 *conf_err){
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);


  /*REAL8 freq;
  REAL4 LE;
  REAL8 OrbRadius;
  REAL8 a;
  REAL8 b;*/
  
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
    /*while( fscanf(fpMCInjectionResults, "%lf%f%lf%lf%lf\n", &freq, &LE, &OrbRadius, &a, &b) != EOF){*/
  
	 	  if (maxPower >= priorLoudestEvent) count++;
           	  tot++;

              /*  #ifdef DEBUG_CONFIDENCE*/
		  fprintf(stdout,"MCLE %f tot %i count %i SearchLE %f\n", maxPower, tot, count, priorLoudestEvent);
                fflush(stdout);
                   /*#endif*/

 /* }*/

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
 
	INITSTATUS(status);
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

void ValidateMCResults(LALStatus *status,
			REAL4FrequencySeries **SUMData,
			REAL4FrequencySeries **STKData,
			TdotsAndDeltaTs *pTdotsAndDeltaTs,
			StackSlideSearchParams *searchParams,
			INT4 iSky,
			StackSlideParams *params){
	

	 
  INT4 kSUM = 0; /* index gives which SUM */ /* HARD CODED TO 0; This function only returns one SUM! */
  INT4 iSUM = 0; /* index gives which SUM bin */  
  REAL8 f_t;

  REAL4 invNumSTKs = 1.0/((REAL4)params->numSTKs); /* 12/03/04 gam */ /* 02/21/05 gam */

  INT4 i,k,m;  
  INT4 binoffset = 0;  
  REAL8 deltaKappa;
  REAL8 piDeltaKappa;
  
  INT4 iMinSTK = floor((params->f0SUM-params->f0STK)*params->tSTK + 0.5); /* Index of mimimum frequency to include when making SUMs from STKs */
  INT4 iMaxSTK = iMinSTK + params->nBinsPerSUM - 1;                       /* Index of maximum frequency to include when making SUMs from STKs */

  REAL8 refFreq = params->f0SUM + ((REAL8)(params->nBinsPerSUM/2))*params->dfSUM; /* 12/06/04 gam */ /* 02/21/05 gam */

  REAL4Vector *fPlus2;  /* square of F_+ */
  REAL4Vector *fCross2; /* square of F_+ */
  REAL4 nrmFactor = searchParams->threshold3;
  
  REAL4 Aplus = searchParams->threshold4/searchParams->blkRescaleFactor;
  REAL4 Across = searchParams->threshold5;
  REAL4 invSqrtSn;  
  REAL4 d2;
  REAL8 smallFactor = 1.e-3;

  LALDetector cachedDetector ;

	INITSTATUS(status);
        ATTATCHSTATUSPTR (status);

        fprintf(stdout,"start function ValidateMCResults\n ");
        fflush(stdout);
	

  
/*set detector site*/
  
  if (strstr(searchParams->IFO, "LHO")) {
       cachedDetector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  } else if (strstr(searchParams->IFO, "LLO")) {
       cachedDetector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  } else if (strstr(searchParams->IFO, "GEO")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  } else if (strstr(searchParams->IFO, "VIRGO")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  } else if (strstr(searchParams->IFO, "TAMA")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  } else {
      
      ABORT( status, STACKSLIDEBINARYH_EIFO, STACKSLIDEBINARYH_MSGEIFO);
  }    


/* Compute F_+ and F_x at the midpoint time of each SFT */
  fPlus2=(REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  fPlus2->data=(REAL4 *)LALMalloc(searchParams->numSTKs*sizeof(REAL4));
  fCross2=(REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  fCross2->data=(REAL4 *)LALMalloc(searchParams->numSTKs*sizeof(REAL4));  
  
GetDetResponseTStampMidPts(status->statusPtr, fPlus2, searchParams->timeStamps, searchParams->numSTKs, searchParams->tSTK,
           &cachedDetector, searchParams->skyPosData[iSky], searchParams->orientationAngle, COORDINATESYSTEM_EQUATORIAL, 1);
CHECKSTATUSPTR (status);
GetDetResponseTStampMidPts(status->statusPtr, fCross2, searchParams->timeStamps, searchParams->numSTKs, searchParams->tSTK,
           &cachedDetector, searchParams->skyPosData[iSky], searchParams->orientationAngle, COORDINATESYSTEM_EQUATORIAL, 0);
  CHECKSTATUSPTR (status);	


printf("skyPosData[%d][0] %f skyPosData[%d][1] %f\n",iSky, searchParams->skyPosData[iSky][0], iSky, searchParams->skyPosData[iSky][1]);

	 
for(k=0;k<params->numSTKs;k++) {

        /* compute the frequency */
        f_t = refFreq;

        for (m=0; m<params->numSpinDown; m++) {
            f_t += params->freqDerivData[m] * pTdotsAndDeltaTs->vecDeltaTs[k][m];
        }
        f_t *= pTdotsAndDeltaTs->vecTDots[k];
        binoffset = floor(( (f_t - refFreq) * params->tSTK) + 0.5 );

  /* deltaKappa = kappa - k = mismatch in bins between actual frequency and freqency used. */
        deltaKappa = fabs( ((f_t - refFreq) * params->tSTK) - (REAL8)binoffset );
        piDeltaKappa = ((REAL8)LAL_PI)*deltaKappa;

      for (i=iMinSTK;i<=iMaxSTK; i++) {
            iSUM = i - iMinSTK;

 /* properly normalized inverse amplitude spectral density in 1/strain/Hz^{1/2} */
            invSqrtSn = sqrt(searchParams->inverseMedians[k]->data[i+binoffset])/sqrt(nrmFactor);

	    /* optimal signal to noise ratio, d2 */
            d2 = ( (Aplus*searchParams->blkRescaleFactor*invSqrtSn)*(Aplus*searchParams->blkRescaleFactor*invSqrtSn)*fPlus2->data[k] +
                 (Across*searchParams->blkRescaleFactor*invSqrtSn)*(Across*searchParams->blkRescaleFactor*invSqrtSn)*fCross2->data[k] )*((REAL4)searchParams->tBLK);

	    /* include frequency mismatch factor if not basically equal to 1 */ 
	    
	    if (piDeltaKappa > smallFactor) {
                d2 *= ((REAL4)(sin(piDeltaKappa)*sin(piDeltaKappa)/(piDeltaKappa*piDeltaKappa)));
            }
            
      /* Print out info about the beam pattern, noise amplitude spectral density, mismatch (deltaKappa), SNR, and binoffset (doppler modulation) */
            if ((searchParams->debugOptionFlag & 2) > 0 ) {
               if (i == (iMinSTK + iMaxSTK)/2) {
                  if (k == 0) {
                     fprintf(stdout,"\n   gpsSeconds     A_+        |F_+|      A_x       |F_x|     sqrt(Sn)  deltaKappa   d2   binoffset\n\n");
                     fflush(stdout);
                  }
                  fprintf(stdout,"%12i %12.4e %8.4f %12.4e %8.4f %12.4e %8.4f %8.4f %5i\n",searchParams->timeStamps[k].gpsSeconds,Aplus,sqrt(fPlus2->data[k]),Across,sqrt(fCross2->data[k]),1.0/(searchParams->blkRescaleFactor*invSqrtSn),deltaKappa,d2,binoffset);
                  fflush(stdout);
               }
            }

            if (k==0) {
               /* Starting a new SUM: initialize */
               SUMData[kSUM]->data->data[iSUM] = d2;
               SUMData[kSUM]->epoch.gpsSeconds = params->gpsStartTimeSec;
               SUMData[kSUM]->epoch.gpsNanoSeconds = params->gpsStartTimeNan;
               SUMData[kSUM]->f0=params->f0SUM;
               SUMData[kSUM]->deltaF=params->dfSUM;
               SUMData[kSUM]->data->length=params->nBinsPerSUM;
            } else {
               SUMData[kSUM]->data->data[iSUM] += d2;
            }
        }/*end of for iMinSTK*/

	    
}/*end of for k =1 to numSTKs*/
/* Find average d2 and then estimate of StackSlide Power */  
  for(i=0;i<params->nBinsPerSUM; i++) {
        SUMData[kSUM]->data->data[i] = SUMData[kSUM]->data->data[i]*invNumSTKs;
        SUMData[kSUM]->data->data[i] = 1.0 + 0.5*SUMData[kSUM]->data->data[i];
  }

  LALFree(fPlus2->data);
  LALFree(fPlus2);
  LALFree(fCross2->data);
  LALFree(fCross2);


printf("nstk %i sqrt(nrmFactor) %f invSqrtSn %g d2 %g\n", searchParams->numSTKs, sqrt(nrmFactor), invSqrtSn, d2 );
printf("tstk %f, psi %f SUMData[kSUM]->data->data[%i] %f \n",  params->tSTK, searchParams->orientationAngle,i, SUMData[kSUM]->data->data[i-1] );
	


/*               #define PRINT_VALIDATEMC 
 
	
	
	INT4 i=0;
        INT4 iPwrMax = 0;
	REAL4 pwrStdDev = 0.0;
	REAL4 pwrMax = 0.0;
        REAL4 pwrSNR = 0.0;
	
	REAL4 pwrSum= 0.0;    
        REAL4 pwrSum2 = 0.0;  
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

	pwrMean = pwrSum/(oneSUM->data->length); 
        pwrStdDev = pwrSum2/(oneSUM->data->length) - pwrMean*pwrMean;   

	if ( pwrStdDev > 0.0 ) {
             pwrStdDev = sqrt( pwrStdDev );
           } else {
            pwrStdDev = -1.0;            }

        
	pwrMax=pwrMax/pwrStdDev;
        
	*SNR = pwrMax;
    
        #ifdef PRINT_VALIDATEMC
	fprintf(stdout,"length %d pwrStdDev %f\n", oneSUM->data->length, pwrStdDev);
        fflush(stdout);
        #endif
*/	
	CHECKSTATUSPTR (status);
        DETATCHSTATUSPTR (status);

	
}

void ComputeUpperLimit(LALStatus *status, const REAL8 *arrayULs, const REAL8 *arrayConfs, REAL8 desiredConf){
	INITSTATUS(status);
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

        INITSTATUS(status);
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
