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
	
    INT4 i,m;
    static INT4 iSky = 0;       /* index gives which Sky Position */    
    INT8 iSkyCoh=0;
    /*INT4 iFreqDeriv;*/
    INT4 numLoopOnSpindown; 



printf("Start function StackSlideBinary\n");    
	StackSlideSkyParams *csParams;           /* Container for ComputeSky */
	/* INT4 iSkyCoh;                           For ComputeSky */
	
 
       /* REAL4FrequencySeries **STKData;  */
       /* REAL4FrequencySeries **SUMData;*/

	
	TdotsAndDeltaTs *pTdotsAndDeltaTs;

	EarthState earth;
	EmissionTime emit;

	INITSTATUS (status, "StackSlideBinary", STACKSLIDEC); 
        ATTATCHSTATUSPTR(status);

	
	/* Allocate space and set quantities for call to ComputeSky() */
	csParams=(StackSlideSkyParams *)LALMalloc(sizeof(StackSlideSkyParams));
	csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));

        /*SUMData=(REAL4FrequencySeries **)LALMalloc(sizeof(REAL4FrequencySeries *));
        SUMData[0]=(REAL4FrequencySeries *)LALMalloc(sizeof(REAL4FrequencySeries));
        SUMData[0]->data=(REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
        SUMData[0]->data->data=(REAL4 *)LALMalloc(stksldParams->nBinsPerSUM*sizeof(REAL4));*/


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
	

/*        for(iSMA=0; iSMA < iSMAmax; iSMA++){*/
              LALUniformDeviate(status->statusPtr, &randval, randPar);
	      /*CHECKSTATUSPTR (status); */
	      
	      stksldParams->SemiMajorAxis=stksldParams->SMAcentral+((REAL8)randval-0.5)*(stksldParams->deltaSMA);                      csParams->SemiMajorAxis=stksldParams->SemiMajorAxis;

		 
		 /* Call STACKSLIDECOMPUTESKYBINARY() */
	StackSlideComputeSkyBinary(status->statusPtr, pTdotsAndDeltaTs, iSkyCoh, csParams);


   	/*CHECKSTATUSPTR (status);*/

	         /* Call STACKSLIDE() for every point in spindown parameter space*/

	if (stksldParams->numFreqDerivTotal==0) {
	   numLoopOnSpindown = 1;  /* Even if no spindown, need to execute iFreqDeriv loop once to handle case of zero spindown */
	} else {
	   numLoopOnSpindown = stksldParams->numFreqDerivTotal;
	}

	/* for all spindowns, loop */
 	for(stksldParams->iFreqDeriv=0;stksldParams->iFreqDeriv<numLoopOnSpindown;stksldParams->iFreqDeriv++)/*put this outside*/
  	{
 		
printf("numLoop on SpinDown is %d\n",numLoopOnSpindown);


/*for (iFreqDeriv = 0; iFreqDeriv<stksldParams->numSpinDown; iFreqDeriv++)
{*/
	
	StackSlide(status->statusPtr,SUMData,STKData,pTdotsAndDeltaTs,stksldParams); /*temporary definad below but afterwords called in from lal using #include <lal/LALStackSlide.h>*/
/*}*/
	
	} /* for iFreqDeriv = 0 to numFreqDerivTotal */

	/*}end of for iT*/
	/*}end of for iSMA*/
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

/*********************************************************************************/
/*********************************************************************************/
/*              DECLARATION OF FUNCTIONS                                         */
/*********************************************************************************/
/*********************************************************************************/
                                          


                 /* %%%%%%%%%%%%%%%%%% */


/*********************************************************************************/
/*              START function: StackSlide                                       */
/*********************************************************************************/

	
       void StackSlide(	LALStatus *status, 
			REAL4FrequencySeries **SUMData, 
			REAL4FrequencySeries **STKData,
                        TdotsAndDeltaTs *pTdotsAndDeltaTs,
			StackSlideParams *params)
{
       
    /*INT4 iFreqDeriv */; /* index give which Freq Derivs  */
    INT4 kSUM ;       /* index gives which SUM */
    INT4 iSUM ;       /* index gives which SUM bin */
    REAL8 f_t;
    

    INT4 numLoopOnSpindown; 
    
    INT4 i,k,m, iSky;
    INT4 binoffset = 0;  
 
    printf("start function StackSlide\n");
        	
       /*INITSTATUS (status, "StackSlide", STACKSLIDEC); 
       ATTATCHSTATUSPTR(status);*/

  
	INT4 iMinSTK = floor((params->f0SUM-params->f0STK)*params->tSTK + 0.5); /* Index of mimimum frequency                                                                       to include when making SUMs from STKs */
	
        INT4 iMaxSTK = iMinSTK + params->nBinsPerSUM - 1;                       /* Index of maximum frequency                                                                       to include when making SUMs from STKs */
  printf("iFreqDeriv is %d\n",params->iFreqDeriv);



printf("iMaxSTK is %d\n", iMaxSTK);
 	 	
		kSUM = iSky*numLoopOnSpindown + params->iFreqDeriv; /* 01/28/04 gam; need to find index to which to SUM */
 printf("iMinSTK is %d\n",iMinSTK);
      for(k=0;k<params->numSTKs;k++) {
      	    printf("numSTKs is %d\n", params->numSTKs);
      /* compute frequency */
	/* f_t = params->f0STK * params->tSTK; */
	f_t = params->f0STK;
	printf("start ft is %f\n",f_t);
	for (m=0; m<params->numSpinDown; m++) 
       	   {
		   printf("freqDerivData[%d][%d] %f\n",params->iFreqDeriv,m, params->freqDerivData[params->iFreqDeriv][m]);
		   printf("numspindown  is %d m is %d k is %d \n",params->numSpinDown, m, k );
		f_t += params->freqDerivData[params->iFreqDeriv][m] * pTdotsAndDeltaTs->vecDeltaTs[k][m];
		/*params->freqDerivData[m] * pTdotsAndDeltaTs->vecDeltaTs[k][m]; */
	
       	   }	   
	   f_t = f_t * pTdotsAndDeltaTs->vecTDots[k];
	   	
	   printf("corrected ft is %f\n",f_t);

	   binoffset = floor(( (f_t - params->f0STK) * params->tSTK) + 0.5 );
 printf("binoffset is %d\n", binoffset);
           #ifdef DEBUG_STACKSLIDE_FNC
              fprintf(stdout, "In StackSlide for SFT #%i binoffset = %i \n",k,binoffset);
              fflush(stdout);   
           #endif
	printf("start Normalization\n");


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

            }/*end of for iMinSTK*/
      } /* END for(k=0;k<params->numSTKs;k++) */
                               /* Normalize the SUMs with params->numSTKs*/
      
      for(i=0;i<params->nBinsPerSUM; i++) {
               SUMData[kSUM]->data->data[i] =  SUMData[kSUM]->data->data[i]/((REAL4)params->numSTKs);
      printf("Normalized sum is %f\n",SUMData[kSUM]->data->data[i]);

                          	
      }	
      
       printf("end function StackSlide\n");     
	
  /*CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);*/
 }   /*end of void StackSlide*/





/*********************************************************************************/
/*              END function: StackSlide                                         */
/*********************************************************************************/



                         /* %%%%%%%%%%%%%%%%%% */



/*********************************************************************************/
/*              START function: StackSlideComputeSkyBinary                       */
/*********************************************************************************/



static void Ft(LALStatus *status, REAL8 *tr, REAL8 t, void *tr0);


/*Redefinition of all the parameters needed for orbital motion for easier use*/

static REAL8 a;      /* semi major axis -maybe a random*/
static REAL8 Period;    /* Period */
static REAL8 ecc;    /* eccentricity */
static REAL8 parg;   /* argument of periapse -maybe a random*/
static LIGOTimeGPS Tperi;  /* Time of periapse passage as measured in SSB frame */
static REAL8 p,q;    /* coefficients of phase model */
static REAL8 E;      /* eccentric anomaly */



void StackSlideComputeSkyBinary( LALStatus		*status, 
				 TdotsAndDeltaTs 	*pTdotsAndDeltaTs, 
				 INT8 		iSkyCoh,
				 /*StackSlideBinarySkyParams 	*params*/
				 StackSlideSkyParams *params
			        )

{  
printf("start function ComputeSkyBinary\n");

	INT4	m, n, nP;
	REAL8	dTbary;
	LALTimeInterval dTBaryInterval;
	LALTimeInterval HalfSFT;
	REAL8 HalfSFTfloat;
	REAL8   dTbarySP;
	REAL8   dTperi;
	REAL8   dTcoord;
	REAL8   Tdotbin;
	REAL8   basedTperi;
	DFindRootIn input;
	REAL8 tr0;
	REAL8 acc;

	
 INITSTATUS (status, "StackSlideComputeSkyBinary", STACKSLIDEC);
 ATTATCHSTATUSPTR(status);
 

/*for (i=0;i< SSparams->numSTKs; i++)
        	{
			tGPS[i].gpsSeconds=SSparams->gpsStartTimeSec +(UINT4)(i*(SSparams->tSTK));
	        	tGPS[i].gpsNanoSeconds=SSparams->gpsStartTimeNan;
	                       
       	        }*/
 /*is this SSparams or params?*/
  
 
/*#ifdef 0 */
/* Check for non-negativity of sky positions in SkyCoh[] */
 ASSERT(iSkyCoh>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);

/* Check to make sure sky positions are loaded */
 ASSERT(params->skyPos!=NULL, status, COMPUTESKYBINARYH_ENULL, COMPUTESKYBINARYH_MSGENULL);
 ASSERT(params->skyPos!=NULL, status, COMPUTESKYBINARYH_ENULL, COMPUTESKYBINARYH_MSGENULL);
 /* Check to make sure parameters are loaded and reasonable */
 ASSERT(params->spinDwnOrder>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT(params->mObsSFT>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT(params->tSFT>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
	 
/* Check to make sure orbital parameters are loaded and reasonable */
ASSERT(params->SemiMajorAxis>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT(params->OrbitalPeriod>0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT(params->OrbitalEccentricity>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT((params->ArgPeriapse>=0)&&(params->ArgPeriapse<=LAL_TWOPI), status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
 ASSERT(params->TperiapseSSB.gpsSeconds>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA); 
 ASSERT((params->TperiapseSSB.gpsNanoSeconds>=0)&&(params->TperiapseSSB.gpsNanoSeconds<1e9), status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
/*#endif*/ 

  /* Here we redefine the orbital variables for ease of use */
printf("leap %d\n",params->edat->leap);

 a=params->SemiMajorAxis;  /* This is the projected semi-major axis of the orbit normalised by the speed of light */
 Period=params->OrbitalPeriod;  /* This is the period of the orbit in seconds */
 ecc=params->OrbitalEccentricity;  /* This is the eccentricity of the orbit */
 parg=params->ArgPeriapse;  /* This is the argument of periapse defining the angular location of the source at periapsis */
                            /* measured relative to the ascending node */
 Tperi.gpsSeconds=params->TperiapseSSB.gpsSeconds;  /* This is the GPS time as measured in the SSB of the observed */
 Tperi.gpsNanoSeconds=params->TperiapseSSB.gpsNanoSeconds;  /* periapse passage of the source */



 /* Convert half the SFT length to a LALTimeInterval for later use */
 HalfSFTfloat=params->tSFT/2.0;
 fprintf(stdout,"HalfSFTfl %f\n",HalfSFTfloat);
 LALFloatToInterval(status->statusPtr,&HalfSFT,&HalfSFTfloat);

/* fprintf(stdout,"HalfSFT %f\n",HalfSFT.seconds);*/


 /* Here we check that the GPS timestamps are greater than zero */
  for(n=0;n<params->mObsSFT;n++)
   { 

      ASSERT(params->tGPS[n].gpsSeconds>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);
/*ASSERT(tGPS[n].gpsSeconds>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);*/


   }
 /* Check to make sure pointer to output is not NULL */
         /*ASSERT(skyConst!=NULL, status, COMPUTESKYBINARYH_ENNUL, COMPUTESKYBINARYH_MSGENNUL);*/
 
	 ASSERT(pTdotsAndDeltaTs!=NULL, status, COMPUTESKYBINARYH_ENNUL, COMPUTESKYBINARYH_MSGENNUL);

 
 /* prepare params input sky position structure */
/* params->baryinput->alpha=params->skyPos[iSkyCoh];
 params->baryinput->delta=params->skyPos[iSkyCoh+1];*/
 /* params->baryinput->alpha=SSparams->skyPosData[iSky][iSkyCoh];
  params->baryinput->delta=SSparams->skyPosData[iSky][iSkyCoh+1]; maybe put it back*/
 



  /* calculate phase model coefficients p and q which are defined in the ComputeSkyBinary header LAL documentation  */
 
 p=((LAL_TWOPI/(Period*1.0))*a*sqrt(1-(ecc*ecc))*cos(parg))-ecc;
 q=(LAL_TWOPI/(Period*1.0))*a*sin(parg);

 /* Calculate the required accuracy for the root finding procedure in the main loop */
 acc=LAL_TWOPI*(REAL8)ACC/Period;   /* ACC is defined in ComputeSkyBinary.h and represents the required */
                                    /* timing precision in seconds (roughly)*/
 

 /* begin loop over SFT's */
 for (n=0; n<params->mObsSFT; n++) 
   {

     /* Calculate the detector time at the mid point of current SFT ( T(i)+(tsft/2) ) using LAL functions */
 /*!*/ /* LALIncrementGPS(status->statusPtr,&(params->baryinput->tgps),&params->tGPS[n],&HalfSFT);*/
/*	   LALIncrementGPS(status->statusPtr,&(params->baryinput->tgps),&tGPS[n],&HalfSFT);*/
	  

	   
 LALIncrementGPS(status->statusPtr, &(params->baryinput->tgps),&params->tGPS[n],&HalfSFT);
 fprintf(stdout,"tgps %d\n",params->baryinput->tgps.gpsSeconds);


fprintf(stdout,"leap %d\n",params->edat->leap);

     /* Convert this mid point detector time into barycentric time (SSB) */

               /* GetSSBTime(&(baryinput.tgps), &TmidSSB) ;*/

 
       LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat);    
       LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);   

       CHECKSTATUSPTR(status);printf("status %d\n",status->statusCode);
	     
   /*   LALBarycenterEarth(status->statusPtr, csParams->earth, &(csParams->baryinput->tgps), csParams->edat);    
     LALBarycenter(status->statusPtr, csParams->emit, csParams->baryinput, csParams->earth);*/       




     
     /* Calculate the time difference since the observed periapse passage in barycentric time (SSB). */ 
     /* This time difference, when converted to REAL8, should lose no precision unless we are dealing */
     /* with periods >~ 1 Year */
     
       LALDeltaGPS(status->statusPtr,&dTBaryInterval,&(params->emit->te),&Tperi); 
               
       LALIntervalToFloat(status->statusPtr,&dTbary,&dTBaryInterval);

      /*     LALDeltaGPS(status->statusPtr,&dTBaryInterval,&(csParams->emit->te),&Tperi); 
     LALIntervalToFloat(status->statusPtr,&dTbary,&dTBaryInterval);*/

    
     /* Calculate the time since the last periapse passage ( < Single Period (SP) ) */
     dTbarySP=Period*((dTbary/(1.0*Period))-(REAL8)floor(dTbary/(1.0*Period)));
    printf("dTbarySP is %f\n",dTbarySP);
     /* Calculate number of full orbits completed since the input observed periapse passage */
     nP=(INT4)floor(dTbary/(1.0*Period));
      printf("nP is %i\n",nP);

     /* begin root finding procedure */
     tr0 = dTbarySP;        /* we wish to find the value of the eccentric anomaly E corresponding to the time */
                            /* since the last periapse passage */
     input.function = Ft;   /* This is the name of the function we must solve to find E */
     input.xmin = 0.0;      /* We know that E will be found between 0 and 2PI */
     input.xmax = LAL_TWOPI;
     input.xacc = acc;      /* The accuracy of the root finding procedure */

                                                    
     /* expand domain until a root is bracketed */
     LALDBracketRoot(status->statusPtr,&input,&tr0); 
  
     /* bisect domain to find eccentric anomoly E corresponding to the current midpoint timestamp */
     LALDBisectionFindRoot(status->statusPtr,&E,&input,&tr0); 
 
     /* Now we calculate the time interval since the input periapse passage as measured at the source */ 
     dTperi=(Period/LAL_TWOPI)*(E-(ecc*sin(E)))+((REAL8)nP*Period); 
     printf("dTperi is %f\n",dTperi);
   
     /* The following quantity is the derivative of the time coordinate measured at the source with */
     /* respect to the time coordinate measured in the SSB : dt_(source)/dt_(SSB) */
     dTcoord=(1.0-(ecc*cos(E)))/(1.0+(p*cos(E))-(q*sin(E)));  
     printf("dTcoord is %f\n",dTcoord);

   
     
     /* The following quantity is the derivitive of the time coordinate measured in the SSB with */
     /* respect to the time coordinate measured at the chosen detector.  It was calculated via the */
     /* last call to LALBarycenter : dt_(SSB)/dt_(detector)  */
      Tdotbin = params->emit->tDot*dTcoord; /*=dt_source/dt_detector*/ 
    printf("Tdotbin is %f\n",Tdotbin);


      /*    Tdotbin = csParams->emit->tDot*dTcoord; */
     
     /* Loop over all spin down orders plus 0th order (f0) */
     /* In this loop we calculate the SkyConstants defined in the documentation as A_{s,alpha} and B_{s,alpha} */

      
      /*!*/ pTdotsAndDeltaTs->vecTDots[n]= Tdotbin;
 printf("pTdotsAndDeltaTs->vecTDots is %f\n",pTdotsAndDeltaTs->vecTDots[n]);

 
     
        for (m=0; m<params->spinDwnOrder; m++)     
       {
	 /* raise the quantity dTperi to the power m */
	 basedTperi = pow(dTperi, (REAL8)m+1);
	 printf("basedTperi %f\n",basedTperi);
/*!*/	 /*the 2 lines below must be changed */
	 printf("spindownorder %d \n",params->spinDwnOrder);
	 /* Calculate A coefficients */
	/* skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m]=1.0/((REAL8)m+1.0)*basedTperi*dTperi-0.5*params->tSFT*basedTperi*Tdotbin;*/
	 /* Calculate B coefficients */
/*	 skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m+1]= params->tSFT*basedTperi*Tdotbin;*/
	 
 /*expressing the time difference in the SSB */
/*!*/        /*pTdotsAndDeltaTs->vecDeltaTs[m][n]= basedTperi*pTdotsAndDeltaTs->vecTDots[n]; */
	 pTdotsAndDeltaTs->vecDeltaTs[n][m]= basedTperi; /*in the Source*/ 
          printf("vecDeltaTs %f\n", pTdotsAndDeltaTs->vecDeltaTs[n][m]);

       }    
   printf("end function StackSlideComputeSkyBinary\n");    
} 
/*LALFree(params->edat->ephemE);
LALFree(params->edat->ephemS);
LALFree(params->edat);*/

 /* Normal Exit */
 DETATCHSTATUSPTR(status);
 RETURN(status);
}

/**************************************************************************************************/

static void Ft(LALStatus *status, REAL8 *tr, REAL8 lE, void *tr0)
{
  INITSTATUS(status, "Ft", "Function Ft()");
  ASSERT(tr0,status, 1, "Null pointer");

  /* this is the function relating the observed time since periapse in the SSB to the true eccentric anomoly E */

  *tr = *(REAL8 *)tr0*(-1.0) + (Period/LAL_TWOPI)*(lE+(p*sin(lE))+q*(cos(lE)-1.0));
  RETURN(status);
}



/*********************************************************************************/
/*              END function: StackSlideComputeSkyBinary                         */
/*********************************************************************************/
	

