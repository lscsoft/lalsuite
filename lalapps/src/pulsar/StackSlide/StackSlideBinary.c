/*********************************** <lalVerbatim file="StackSlideCV">
Author:  Landry, M., and Mendell, G.
$Id$
************************************ </lalVerbatim> */


#include "StackSlideBinary.h"
NRCSID( STACKSLIDEC,  "$Id$");


/*call the function!!*/
void StackSlideBinary(      LALStatus *status,
                        StackSlideParams *stksldParams,
			 REAL4FrequencySeries **STKData,
			REAL4FrequencySeries **SUMData
			)
			
{
	
    INT4 i,m,iSMA;
    static INT4 iSky = 0;       /* index gives which Sky Position */    
    INT8 iSkyCoh=0;
    /*INT4 iFreqDeriv;*/
    INT4 numLoopOnSpindown; 



printf("Start function StackSlideBinary\n");    
	StackSlideSkyParams *csParams;           /* Container for ComputeSky */
 
      
	
	TdotsAndDeltaTs *pTdotsAndDeltaTs;

	EarthState earth;
	EmissionTime emit;

	INITSTATUS (status, "StackSlideBinary", STACKSLIDEC); 
        ATTATCHSTATUSPTR(status);

	
	/* Allocate space and set quantities for call to ComputeSky() */
	csParams=(StackSlideSkyParams *)LALMalloc(sizeof(StackSlideSkyParams));
	csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));

   

  	pTdotsAndDeltaTs=(TdotsAndDeltaTs *)LALMalloc(sizeof(TdotsAndDeltaTs));
  	pTdotsAndDeltaTs->vecTDots  = (REAL8 *)LALMalloc(sizeof(REAL8)*stksldParams->numSTKs);
  	
	if (stksldParams->numSpinDown>0) {
	  pTdotsAndDeltaTs->vecDeltaTs  = (REAL8 **)LALMalloc(sizeof(REAL8 *)*stksldParams->numSTKs);          
	  for(i=0;i<stksldParams->numSTKs;i++)
              {
       		pTdotsAndDeltaTs->vecDeltaTs[i]  = (REAL8 *)LALMalloc(sizeof(REAL8)*stksldParams->numSpinDown);          
              }
	    }

  REAL4 randval;
   RandomParams *randPar=NULL;
   INT4 seed=0;

 for(iSky=0;iSky<stksldParams->numSkyPosTotal;iSky++)
  	{
	 
	csParams->skyPos[0]=stksldParams->skyPosData[iSky][0];/*remember:this must be the same as baryinput.alpha*/
	csParams->skyPos[1]=stksldParams->skyPosData[iSky][1];/*remember:this must be the same as baryinput.delta*/
	
	stksldParams->baryinput->alpha=stksldParams->skyPosData[iSky][0];
	stksldParams->baryinput->delta=stksldParams->skyPosData[iSky][1]; /*these are specifyed in SetupBaryInput function*/
	
	csParams->tGPS=stksldParams->timeStamps;
        csParams->gpsStartTimeSec = stksldParams->gpsStartTimeSec; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
        csParams->gpsStartTimeNan = stksldParams->gpsStartTimeNan; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
	csParams->spinDwnOrder=stksldParams->numSpinDown;
	csParams->mObsSFT=stksldParams->numSTKs;
	csParams->tSFT=stksldParams->tSTK;
        csParams->edat=stksldParams->edat;
        csParams->emit=&emit;
        csParams->earth=&earth;
        csParams->baryinput=stksldParams->baryinput;

              	
	/* what the hell does this number mean? */
        for(i=0;i<stksldParams->numSTKs;i++)
         {
      	  pTdotsAndDeltaTs->vecTDots[i] = 0.0; /* Initialize for each sky position*/
          for(m=0;m<stksldParams->numSpinDown;m++)
            {
      	     pTdotsAndDeltaTs->vecDeltaTs[i][m] = 0.0; /* Initialize for each sky position */

            }
        }




/*ADD HERE BINARY PARAMETERS to be input by means of greg's command line arg.*/

	csParams->OrbitalPeriod=66960; /*no need to go into the commandline*/
	csParams->OrbitalEccentricity=stksldParams->OrbitalEccentricity; /*remember: params must be the comm line variable*/
        csParams->ArgPeriapse=stksldParams->ArgPeriapse;
        csParams->TperiapseSSB.gpsSeconds=stksldParams->TperiapseSSBSec;
        csParams->TperiapseSSB.gpsNanoSeconds=stksldParams->TperiapseSSBNanoSec;
/*these are passed through the command line, MAKE SURE ABOUT CALLING params structure*/

	 for(iSMA=0; iSMA < stksldParams->nMaxSMA; iSMA++){ /*05/02/18 vir: for each point generate a new random*/
	/* for (iT=0; iT < stksldParams->nMaxTperi; iT++)*/
LALCreateRandomParams(status->statusPtr, &randPar, seed); /*05/02/17 vir: create random params*/
              LALUniformDeviate(status->statusPtr, &randval, randPar);
	      /*CHECKSTATUSPTR (status); */
	      
	      stksldParams->SemiMajorAxis=stksldParams->SMAcentral+((REAL8)randval-0.5)*(stksldParams->deltaSMA);                         
	      csParams->SemiMajorAxis=stksldParams->SemiMajorAxis;
	      /*stksldParams->TperiapseSSBsec=params->Tpericentral+((REAL8)randval-0.5)*(stksldParams->deltaTperi); */
              
	      /*05/02/18 vir: csParams->TperiapseSSBsec=stksldParams->TperiapseSSBsec; this is now assigned in CL*/	 
	
	      /* Call STACKSLIDECOMPUTESKYBINARY() */
	StackSlideComputeSkyBinary(status->statusPtr, pTdotsAndDeltaTs, iSkyCoh, csParams);


   	/*CHECKSTATUSPTR (status);*/

	         /* Call STACKSLIDE() for every point in spindown and SMA parameter space*/

	if (stksldParams->numFreqDerivTotal==0) {
	   numLoopOnSpindown = 1;  /* Even if no spindown, need to execute iFreqDeriv loop once to handle case of zero spindown */
	} else {
	   numLoopOnSpindown = stksldParams->numFreqDerivTotal;
	}

	/* for all spindowns, loop */
 	for(stksldParams->iFreqDeriv=0;stksldParams->iFreqDeriv<numLoopOnSpindown;stksldParams->iFreqDeriv++)/*put this outside*/
  	{
 		


/*for (iFreqDeriv = 0; iFreqDeriv<stksldParams->numSpinDown; iFreqDeriv++)
{*/
	
	StackSlide(status->statusPtr,SUMData,STKData,pTdotsAndDeltaTs,stksldParams); /*temporary definad below but afterwords called in from lal using #include <lal/LALStackSlide.h>*/
/*}*/
	} /* for iFreqDeriv = 0 to numFreqDerivTotal */


LALDestroyRandomParams(status->statusPtr, &randPar );

           }/*end of for iSMA 05/02/17 vir*/
	 /*}end of for iT*/
	}/*end of for iSky*/
           /* Deallocate space for ComputeSky parameters */
	LALFree(csParams->skyPos);
	LALFree(csParams);
      /* Deallocate memory */
        /* pTdotsAndDeltaTs->vecDeltaTs is allocated only if numSpinDown > 0. */      
	if (stksldParams->numSpinDown>0) {

	  for(i=0;i<stksldParams->numSTKs;i++) 
	  { 

     		LALFree(pTdotsAndDeltaTs->vecDeltaTs[i]);          
          }
  	 
	  LALFree(pTdotsAndDeltaTs->vecDeltaTs);
	                           }
  	LALFree(pTdotsAndDeltaTs->vecTDots);
	LALFree(pTdotsAndDeltaTs);
  
        

}/*end of void StackSlideBinary()*/

