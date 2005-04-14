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
/* #define DEBUG_BINARY_CODE */
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
        
        /* loop over orbit semimajor axis */
        for(iSMA=0; iSMA < params->nMaxSMA; iSMA++){ /*05/02/18 vir: for each point generate a new random*/
         
         /* loop over Tperi not implemented */
         /* for (iT=0; iT < params->nMaxTperi; iT++)*/
            LALCreateRandomParams(status->statusPtr, &randPar, seed); /*05/02/17 vir: create random params*/
            LALUniformDeviate(status->statusPtr, &randval, randPar);
            /*CHECKSTATUSPTR (status); */

            SemiMajorAxis=params->SMAcentral+((REAL8)randval-0.5)*(params->deltaSMA);
            printf("SMA is %f\n",SemiMajorAxis);
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
                      fprintf(binaryfp,"%f %f\n", params->f0SUM + (REAL8)k*params->dfSUM, params->SUMData[0]->data->data[k]);
                  }
                  fclose(binaryfp);
                  /* 05/02/17 vir: write sum on the file only for no mismatch*/
               } else {
                  /*Look for loudest event in SUMData*/
                  /*05/02/17 vir: for mismatched params, output only the loudest event */
                  #ifdef DEBUG_BINARY_CODE
                     fprintf(stdout,"the SemiMajorAxis is %f\n", SemiMajorAxis);
                     fflush(stdout);
                  #endif
                  FindBinaryLoudest(&LoudestEvent, &peakFreq, params->SUMData, stksldParams);
                  fprintf(binaryLE,"%f %f %f\n", peakFreq, LoudestEvent, ParamsSMA[kSUM]);
               }/*end of else deltaSMA > 0*/
            } /* for iFreqDeriv = 0 to numFreqDerivTotal */

            LALDestroyRandomParams(status->statusPtr, &randPar );

         /*}end of for iT*/
        }/*end of for iSMA 05/02/17 vir*/
        fclose(binaryLE);    
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
void FindBinaryLoudest(REAL8 *LoudestEvent, REAL8 *peakFreq, REAL4FrequencySeries **SUMData, StackSlideParams *stksldParams)
{
        REAL8 max=0;

        INT4 iMinSTK=0; 
        INT4 iMaxSTK=stksldParams->nBinsPerSUM;

        INT4 i;
        INT4 indexFreq=0;

        for (i=iMinSTK; i<iMaxSTK; i++)
        {
                if(SUMData[0]->data->data[i] > max)
                { max = SUMData[0]->data->data[i];
                  indexFreq=i;
                }
                *peakFreq=stksldParams->f0SUM + indexFreq*stksldParams->dfSUM;
        }
        *LoudestEvent=max;
       
        #ifdef DEBUG_BINARY_CODE
         fprintf(stdout, "Loudest binary peak is %f and corr freq is %f SMA %f\n", max, *peakFreq, SemiMajorAxis);
         fflush(stdout);
        #endif
}
