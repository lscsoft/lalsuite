/*
*  Copyright (C) 2007 Virginia Re
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

/*Test code for stackslide binary. It uses fake hard-coded Fourier transforms */
#include <StackSlideBinaryTest.h>

NRCSID(STACKSLIDEBINARYTESTC , "$Id$");

INT4 lalDebugLevel=7;

              /******************/

    REAL8 y;  /*temporary variables to read into eph. files*/
    INT4 x,z;
    INT2 leap;
    INT4 iSFT;           /*index to which fakeSFT*/
    INT4 i,l;
    REAL8 try;
    extern char *optarg;
          
    INT4 iSky = 0;
    INT8 iSkyCoh=0;
      
    LALDetector cachedDetector;

    StackSlideBinarySkyParams *params;     /* Container for ComputeSky/ComputeSkyBinary Parameters*/
    TdotsAndDeltaTs *pTdotsAndDeltaTs;
    StackSlideParams *SSparams;

    BarycenterInput baryinput;
    EarthState earth;
    EmissionTime emit;
    EphemerisData *edat=NULL;
    
    FILE *fpE = NULL;
    FILE *fpS = NULL;


    REAL8 **SUMData;
    REAL8 **fakeSFT;   /*2-d array containing fake SFT data*/
    static LIGOTimeGPS	*tGPS;	
    LIGOTimeGPS TstartDet;
    LIGOTimeGPS TmidSSB;

    static LALStatus status;  /*status structure*/

    
/***************************************************************************/
/*                                  					   */
/*                  START Function : main				   */
/* 									   */
/***************************************************************************/

int main()
{
status.statusCode = 0;
status.statusPtr = NULL;
status.statusDescription = NULL;



/***************************************************************************************/	
          /* Allocate space and set quantities for call to ComputeSky() */
/***************************************************************************************/
    
	params=(StackSlideBinarySkyParams *)LALMalloc(sizeof(StackSlideBinarySkyParams));
	params->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
        pTdotsAndDeltaTs=(TdotsAndDeltaTs *)LALMalloc(sizeof(TdotsAndDeltaTs));
	SSparams=(StackSlideParams *)LALMalloc(sizeof(StackSlideParams));
        	
	pTdotsAndDeltaTs->vecTDots  = (REAL8 *)LALMalloc(6*sizeof(REAL8 *));

       /* params->tGPS=(LIGOTimeGPS *)LALMalloc(6*sizeof(LIGOTimeGPS));*/
	tGPS=(LIGOTimeGPS *)LALMalloc(6*sizeof(LIGOTimeGPS));
             
                   /*allocate memory for fake SFTs arrays*/
    
   SUMData=(REAL8 **)LALMalloc(sizeof(REAL8 *)*6);
    
    for (iSFT=0; iSFT<6; iSFT++){
	   SUMData[iSFT]=(REAL8 *)LALMalloc(7*sizeof(REAL8 ));
                                }
   
  
   fakeSFT=(REAL8 **)LALMalloc(sizeof(REAL8 *)*6);
 
    for (iSFT=0; iSFT<6; iSFT++){
   
   fakeSFT[iSFT]=(REAL8 *)LALMalloc(7*sizeof(REAL8 ));
       if (iSFT==0){
    fakeSFT[iSFT][0]=1;
    fakeSFT[iSFT][1]=2;
    fakeSFT[iSFT][2]=4;
    fakeSFT[iSFT][3]=10;
    fakeSFT[iSFT][4]=3;
    fakeSFT[iSFT][5]=0;
    fakeSFT[iSFT][6]=0;

    }
    else if (iSFT==1){
    fakeSFT[iSFT][0]=0;
    fakeSFT[iSFT][1]=0;
    fakeSFT[iSFT][2]=2;
    fakeSFT[iSFT][3]=0;
    fakeSFT[iSFT][4]=0;
    fakeSFT[iSFT][5]=0;
    fakeSFT[iSFT][6]=0;

    }

    else if (iSFT==2){
    fakeSFT[iSFT][0]=0;
    fakeSFT[iSFT][1]=1;
    fakeSFT[iSFT][2]=0;
    fakeSFT[iSFT][3]=0;
    fakeSFT[iSFT][4]=0;
    fakeSFT[iSFT][5]=0;
    fakeSFT[iSFT][6]=0;

    }
    else if (iSFT==3){
    fakeSFT[iSFT][0]=0;
    fakeSFT[iSFT][1]=2;
    fakeSFT[iSFT][2]=0;
    fakeSFT[iSFT][3]=0;
    fakeSFT[iSFT][4]=0;
    fakeSFT[iSFT][5]=0;
    fakeSFT[iSFT][6]=0;
    }
    else if (iSFT==4){
    fakeSFT[iSFT][0]=0;
    fakeSFT[iSFT][1]=0;
    fakeSFT[iSFT][2]=0;
    fakeSFT[iSFT][3]=0;
    fakeSFT[iSFT][4]=1;
    fakeSFT[iSFT][5]=0;
    fakeSFT[iSFT][6]=0;
    }
    else if (iSFT==5){
    fakeSFT[iSFT][0]=0;
    fakeSFT[iSFT][1]=0;
    fakeSFT[iSFT][2]=0;
    fakeSFT[iSFT][3]=1;
    fakeSFT[iSFT][4]=0;
    fakeSFT[iSFT][5]=0;
    fakeSFT[iSFT][6]=0;

    }
    }/*end of for iSFT=1 to 7*/

       
	
	/*assign values to some key parameters*/

	SSparams->numSTKs=6;
	SSparams->numSpinDown=0;
        SSparams->numSkyPosTotal=1;
        SSparams->gpsStartTimeSec=721163027;/*714153733+2000*60;*/ /*very first gps start time at the detector*/
	SSparams->gpsStartTimeNan=0;
        SSparams->tSTK=60.0;
	SSparams->numFreqDerivTotal=1;
        SSparams->f0STK=6.0400000000000000e+02;
	params->mObsSFT=6;
	TstartDet.gpsSeconds=714153733;
	TstartDet.gpsNanoSeconds=0;
	 
	/* BINARY PARAMETERS */
	
        params->SemiMajorAxis=1.44;
        params->OrbitalPeriod=66960;
	params->OrbitalEccentricity=0;
        params->ArgPeriapse=0;	
	params->TperiapseSSB.gpsSeconds=714153733;
	params->TperiapseSSB.gpsNanoSeconds=0;
            
	
	  if (SSparams->numSpinDown>0) 
	{
 	  pTdotsAndDeltaTs->vecDeltaTs  = (REAL8 **)LALMalloc(sizeof(REAL8 *)*SSparams->numSTKs);          
	  for(i=0;i<SSparams->numSTKs;i++)
              {
  		pTdotsAndDeltaTs->vecDeltaTs[i]  = (REAL8 *)LALMalloc(sizeof(REAL8)*SSparams->numSpinDown);          
              }
        }

 	/* what the hell does this number mean? */
      for(i=0;i<SSparams->numSTKs;i++)
         {
      	  pTdotsAndDeltaTs->vecTDots[i] = 0.0; /* Initialize */
          for(l=0;l<SSparams->numSpinDown;l++)
            {
      	     pTdotsAndDeltaTs->vecDeltaTs[i][l] = 0.0; /* Initialize */
            }
        }



/***************************************************************************************/
/**          End of memory allocation and parameters assignment                       **/	
/***************************************************************************************/


         params->gpsStartTimeSec = SSparams->gpsStartTimeSec; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
         params->gpsStartTimeNan = SSparams->gpsStartTimeNan; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
         params->spinDwnOrder=SSparams->numSpinDown;
         params->mObsSFT=SSparams->numSTKs;
         params->tSFT=SSparams->tSTK;
         params->edat=SSparams->edat;
         params->gpsTperiapseSSBSec=params->TperiapseSSB.gpsSeconds;

    
        	 
	/* LALInitBarycenter(&status,SSparams->edat);*/
         
             /*CALL SetupBaryInput to read in the ephemeris files*/
		 
          if (SetupBaryInput()) printf("ok\n");

	  
       
          /* call STACKSLIDECOMPUTESKYBINARY() */

       	 StackSlideComputeSkyBinary(&status, pTdotsAndDeltaTs, iSkyCoh, params);

       
       
	/* call STACKSLIDEBINARY() the name will change into stackslide*/


        StackSlideBinary(&status,SUMData, fakeSFT, pTdotsAndDeltaTs , SSparams );

    return status.statusCode;



         /*Free allocated memory*/

         if (FreeMem())return 3;


}/*END of main*/



/***************************************************************************/
/*                                  					   */
/*                  END Function : main				   */
/* 									   */
/***************************************************************************/



	 
/*********************************************************************************/
/*********************************************************************************/
/*              DECLARATION OF FUNCTIONS                                         */
/*********************************************************************************/
/*********************************************************************************/
                                          


                                  /* %%%%%%%%%%%%%%%%%% */

/*********************************************************************************/
/*              START function: StackSlideBinary                                 */
/*********************************************************************************/

 void StackSlideBinary(	LALStatus *status, 
			REAL8 **SUMData,
			REAL8 **fakeSFT,
			TdotsAndDeltaTs *pTdotsAndDeltaTs,
			StackSlideParams *SSparams)
{
                       printf("start function StackSlidebinary\n");
		       
	REAL8 f_t;
	 
	
	SSparams->f0SUM = 604.02;
  	SSparams->nBinsPerSUM = 4;
  	SSparams->numSTKs = 6;
        SSparams->dfSUM = 1/60;
  REAL8 refFreq = SSparams->f0SUM + ((REAL8)(SSparams->nBinsPerSUM/2))*SSparams->dfSUM; /* 04/04/05 gam */

	INT4 iMinSTK = floor((SSparams->f0SUM-SSparams->f0STK)*SSparams->tSTK + 0.5); /* Index of mimimum frequency                                                                       to include when making SUMs from STKs */
	
        INT4 iMaxSTK = iMinSTK + SSparams->nBinsPerSUM - 1;                    /* Index of maximum frequency                                                                       to include when making SUMs from STKs */
#ifdef DEBUG_STACKSLIDE_FNC
printf("iMinSTK is %d and iMAxSTK is %d\n",iMinSTK,iMaxSTK);
#endif

    INT4 numLoopOnSpindown; /* 1 if numSpindown is 0, else = params->numFreqDerivTotal */
    INT4 i,k,m, iFreqDeriv, iSUM, kSUM, kmax;
    INT4 binoffset = 0;

     
   
    REAL8 FDeriv1 = 0;
    REAL8 FDeriv2 = 0;
    REAL8 FDeriv3 = 0;
  /*  REAL8 FDeriv4 = 0;*/
    
    
	  	
     INITSTATUS (status, "StackSlideBinary", STACKSLIDEBINARYTESTC); 
      ATTATCHSTATUSPTR(status);
	
        	

        /* 02/02/04 gam; moved next 5 lines from inside loops: */
  	if (SSparams->numFreqDerivTotal==0) {
	   numLoopOnSpindown = 1;  /* Even if no spindown, need to execute iFreqDeriv loop once to handle case of zero spindown */
	} else {
  	   numLoopOnSpindown = SSparams->numFreqDerivTotal;
	}
#ifdef DEBUG_STACKSLIDE_FNC
	printf("numloopOnSpinDown %d and numSpinDown is %d\n",numLoopOnSpindown, SSparams->numSpinDown);
#endif


if (SSparams->numSpinDown > 0) {
  /* SSparams->freqDerivData=(REAL8 **)LALMalloc(SSparams->numFreqDerivTotal*sizeof(REAL8 *));
   for(i=0;i<SSparams->numFreqDerivTotal;i++)*/
	SSparams->freqDerivData=(REAL8 **)LALMalloc(SSparams->numSTKs*sizeof(REAL8 *));
   for(i=0;i<SSparams->numSTKs;i++)
   {
        SSparams->freqDerivData[i] = (REAL8 *)LALMalloc(SSparams->numSpinDown*sizeof(REAL8));
   }
 }

if(SSparams->numSpinDown > 0){
	for (k=0; k<SSparams->numSTKs; k++){
        	for (m=0 ; m<SSparams->numSpinDown ; m++ )
	      {
		if (m==0){
		     SSparams->freqDerivData[k][m]=FDeriv1;  
		     } else if(m==1){
		     SSparams->freqDerivData[k][m]=FDeriv2;  
		     } else if(m==2){
		     SSparams->freqDerivData[k][m]=FDeriv3;  
		     }else {
			     SSparams->freqDerivData[k][m]=0;
		     }

	      }

}
}
                         /*START STACKSLIDING*/	
	
	/* for all spindowns, loop */
 	for(iFreqDeriv=0;iFreqDeriv<numLoopOnSpindown;iFreqDeriv++)
  	{
 	/* call computesky for all stacks at once */
 	/* compute offset  for each stack */
	
	kSUM = iSky*numLoopOnSpindown + iFreqDeriv; /* 01/28/04 gam; need to find index to which to SUM */
	
 #ifdef DEBUG_STACKSLIDE_FNC
	fprintf(stdout,"kSUM %d \n",kSUM);
 #endif

        for(k=0;k<SSparams->numSTKs;k++) {
      	    
      /* compute frequency */
	
	/* f_t = SSparams->f0STK; 04/04/05 vir */
          f_t=refFreq;
	
	
	for (m=0; m<SSparams->numSpinDown; m++) 
       	   {
		/* f_t += ( params->freqDerivData[iFreqDeriv][m] * pTdotsAndDeltaTs[2*k*(params->numSpinDown)+2*(INT4)m+1]); */            
  #ifdef DEBUG_STACKSLIDE_FNC 
   printf("freqDerivData[%d][%d] is %f\n", iFreqDeriv,m , SSparams->freqDerivData[iFreqDeriv][m] );
  #endif
/*!*/		f_t += SSparams->freqDerivData[iFreqDeriv][m] * pTdotsAndDeltaTs->vecDeltaTs[k][m];
       	   }	   
	      f_t = f_t * pTdotsAndDeltaTs->vecTDots[k];
	      
      fprintf(stdout,"ft after correcting %f \n",f_t);
   
	   /*binoffset = floor(( (f_t - SSparams->f0STK) * SSparams->tSTK) + 0.5 );04/04/05 vir*/
           binoffset = floor(( (f_t - refFreq) * SSparams->tSTK) + 0.5 ); /*04/04/05 vir*/

	   fprintf(stdout,"binoffset  %d \n",binoffset);
           
           #ifdef DEBUG_STACKSLIDE_FNC
              fprintf(stdout, "In StackSlide for SFT #%i binoffset = %i \n",k,binoffset);
              fflush(stdout);   
           #endif
 	
	    for(i=iMinSTK;i<=iMaxSTK; i++) {
                iSUM = i - iMinSTK;
	        if (k==0) {
		   SUMData[k][iSUM]=fakeSFT[k][i+binoffset];
		   
  #ifdef DEBUG_STACKSLIDE_FNC	  
    printf("at the step k= %d and iSUM=%d the partial sum SUMData[%d][%d] is %f\n",k ,iSUM,k,iSUM,SUMData[k][iSUM]);
  #endif
	/* Starting a new SUM: initialize */
	          /* SUMData[kSUM]->data->data[iSUM] = STKData[k]->data->data[i+binoffset];
                   SUMData[kSUM]->epoch.gpsSeconds = params->gpsStartTimeSec;
		   SUMData[kSUM]->epoch.gpsNanoSeconds = params->gpsStartTimeNan;
		   SUMData[kSUM]->f0=params->f0SUM;
		   SUMData[kSUM]->deltaF=params->dfSUM;
		   SUMData[kSUM]->data->length=params->nBinsPerSUM;*/
	        } else {
		      SUMData[k][iSUM]=SUMData[k-1][iSUM]+fakeSFT[k][i+binoffset];
#ifdef DEBUG_STACKSLIDE_FNC
 printf("at the step k= %d and iSUM=%d the partial sum SUMData[%d][%d] is %f\n",k ,iSUM,k,iSUM,SUMData[k][iSUM]);
#endif
		               }
		
            }/* end of for i=iMinSTK to iMaxSTK*/
      } /* END for(k=0;k<params->numSTKs;k++) */
 
	kmax=SSparams->numSTKs-1;	    
      
	/* Normalize the SUMs with params->numSTKs*/
              for(i=0;i<SSparams->nBinsPerSUM; i++) {
        
		      /* SUMData[kSUM]->data->data[i] =  SUMData[kSUM]->data->data[i]/((REAL4)params->numSTKs);*/
      	SUMData[kmax][i] =  SUMData[kmax][i]/((REAL4)SSparams->numSTKs) ;
	
   
     printf("Final Normalized  sum is SUMData[%d][%d] is %f\n",kmax, i,SUMData[kmax][i]);
   
	      }	
 

      } /* for iFreqDeriv = 0 to numFreqDerivTotal */

	/* Deallocate memory */
        /*  pTdotsAndDeltaTs->vecDeltaTs is allocated only if numSpinDown > 0. */      
	

	if (SSparams->numSpinDown > 0){
		for(i=1; i< SSparams->numFreqDerivTotal; i++)
		{
				LALFree(SSparams->freqDerivData[i]);
		}
        LALFree(SSparams->freqDerivData);

	}


  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);

    
printf("end of function StackSlidebinary\n");	
 }/*end of StackSlide function*/


/*********************************************************************************/
/*              END function: StackSlideBinary                                   */
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
				 StackSlideBinarySkyParams 	*params
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

	
 INITSTATUS (status, "StackSlideComputeSkyBinary", STACKSLIDEBINARYTESTC);
 ATTATCHSTATUSPTR(status);
 
/*&&&*/
for (i=0;i< SSparams->numSTKs; i++)
        	{
			tGPS[i].gpsSeconds=SSparams->gpsStartTimeSec +(UINT4)(i*(SSparams->tSTK));
	        	tGPS[i].gpsNanoSeconds=0;
	                       
       	        }
 
 
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

 a=params->SemiMajorAxis;  /* This is the projected semi-major axis of the orbit normalised by the speed of light */
 Period=params->OrbitalPeriod;  /* This is the period of the orbit in seconds */
 ecc=params->OrbitalEccentricity;  /* This is the eccentricity of the orbit */
 parg=params->ArgPeriapse;  /* This is the argument of periapse defining the angular location of the source at periapsis */
                            /* measured relative to the ascending node */
 Tperi.gpsSeconds=params->TperiapseSSB.gpsSeconds;  /* This is the GPS time as measured in the SSB of the observed */
 Tperi.gpsNanoSeconds=params->TperiapseSSB.gpsNanoSeconds;  /* periapse passage of the source */



 /* Convert half the SFT length to a LALTimeInterval for later use */
 HalfSFTfloat=params->tSFT/2.0;
/* fprintf(stdout,"HalfSFTfl %f\n",HalfSFTfloat);*/
 LALFloatToInterval(status->statusPtr,&HalfSFT,&HalfSFTfloat);

/* fprintf(stdout,"HalfSFT %f\n",HalfSFT.seconds);*/


 /* Here we check that the GPS timestamps are greater than zero */
  for(n=0;n<params->mObsSFT;n++)
   { 

      /*ASSERT(params->tGPS[n].gpsSeconds>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);*/
ASSERT(tGPS[n].gpsSeconds>=0, status, COMPUTESKYBINARYH_ENEGA, COMPUTESKYBINARYH_MSGENEGA);


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
 

 /* begin loop over SFT's */ /*NOTE: in case the stack duration is > SFT duration, this must be params->numSTKs*/
 for (n=0; n<params->mObsSFT; n++) 
   {

     /* Calculate the detector time at the mid point of current SFT ( T(i)+(tsft/2) ) using LAL functions */
 /*!*/ /* LALIncrementGPS(status->statusPtr,&(params->baryinput->tgps),&params->tGPS[n],&HalfSFT);*/
/*	   LALIncrementGPS(status->statusPtr,&(params->baryinput->tgps),&tGPS[n],&HalfSFT);*/
	   
	   
 LALIncrementGPS(status->statusPtr,&(baryinput.tgps),&tGPS[n],&HalfSFT);

 /*fprintf(stdout,"half tgps in detector time %d\n",baryinput.tgps.gpsSeconds);*/


     /* Convert this mid point detector time into barycentric time (SSB) */

                GetSSBTime(&(baryinput.tgps), &TmidSSB) ;

 
   /*LALBarycenterEarth(status->statusPtr, &earth, &(params->baryinput->tgps), params->edat);    */
           CHECKSTATUSPTR(status);printf("status %d\n",status->statusCode);
	     
  /* LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);    */
  /*   LALBarycenterEarth(status->statusPtr, csParams->earth, &(csParams->baryinput->tgps), csParams->edat);    
     LALBarycenter(status->statusPtr, csParams->emit, csParams->baryinput, csParams->earth);*/       




     
     /* Calculate the time difference since the observed periapse passage in barycentric time (SSB). */ 
     /* This time difference, when converted to REAL8, should lose no precision unless we are dealing */
     /* with periods >~ 1 Year */
     LALDeltaGPS(status->statusPtr,&dTBaryInterval,&(emit.te),&Tperi); 
        
  
     LALIntervalToFloat(status->statusPtr,&dTbary,&dTBaryInterval);
/*     LALDeltaGPS(status->statusPtr,&dTBaryInterval,&(csParams->emit->te),&Tperi); 
     LALIntervalToFloat(status->statusPtr,&dTbary,&dTBaryInterval);*/

#ifdef DEBUG_STACKSLIDE_FNC
     printf("dTbary is %f\n",dTbary);
#endif
     /* Calculate the time since the last periapse passage ( < Single Period (SP) ) */
     dTbarySP=Period*((dTbary/(1.0*Period))-(REAL8)floor(dTbary/(1.0*Period)));

#ifdef DEBUG_STACKSLIDE_FNC
     printf("dTbarySP is %f\n",dTbarySP);
#endif
     
     /* Calculate number of full orbits completed since the input observed periapse passage */
     nP=(INT4)floor(dTbary/(1.0*Period));
 
#ifdef DEBUG_STACKSLIDE_FNC     
     printf("nP is %i\n",nP);
#endif
     
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
 
  #ifdef DEBUG_STACKSLIDE_FNC
   printf("dTperi is %f\n",dTperi);
  #endif
  
   /* The following quantity is the derivative of the time coordinate measured at the source with */
     /* respect to the time coordinate measured in the SSB : dt_(source)/dt_(SSB) */
     dTcoord=(1.0-(ecc*cos(E)))/(1.0+(p*cos(E))-(q*sin(E)));  
   
  #ifdef DEBUG_STACKSLIDE_FNC
   printf("dTcoord is %f\n",dTcoord);
  #endif
   
     
     /* The following quantity is the derivitive of the time coordinate measured in the SSB with */
     /* respect to the time coordinate measured at the chosen detector.  It was calculated via the */
     /* last call to LALBarycenter : dt_(SSB)/dt_(detector)  */
      Tdotbin = emit.tDot*dTcoord; /*=dt_source/dt_detector*/ 

 #ifdef DEBUG_STACKSLIDE_FNC     
      printf("n is %d\n",n);
      printf("Tdotbin is %f\n",Tdotbin);
      printf("Tdot is %f\n",emit.tDot);
#endif

      /*    Tdotbin = csParams->emit->tDot*dTcoord; */
     
     /* Loop over all spin down orders plus 0th order (f0) */
     /* In this loop we calculate the SkyConstants defined in the documentation as A_{s,alpha} and B_{s,alpha} */

      
      /*!*/ pTdotsAndDeltaTs->vecTDots[n]= Tdotbin;
 #ifdef DEBUG_STACKSLIDE_FNC
      printf("pTdotsAndDeltaTs->vecTDots is %f\n",pTdotsAndDeltaTs->vecTDots[n]);
#endif
 
     
        for (m=0; m<params->spinDwnOrder; m++)     
       {
	 /* raise the quantity dTperi to the power m */
	 basedTperi = pow(dTperi, (REAL8)m+1);
	 
#ifdef DEBUG_STACKSLIDE_FNC
	 printf("basedTperi %f\n",basedTperi);
#endif
	
	 
	 /* Calculate A coefficients */
	/* skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m]=1.0/((REAL8)m+1.0)*basedTperi*dTperi-0.5*params->tSFT*basedTperi*Tdotbin;*/
	 /* Calculate B coefficients */
/*	 skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m+1]= params->tSFT*basedTperi*Tdotbin;*/
	 
 /*expressing the time difference in the SSB */
/*!*/        /*pTdotsAndDeltaTs->vecDeltaTs[m][n]= basedTperi*pTdotsAndDeltaTs->vecTDots[n]; */
	 pTdotsAndDeltaTs->vecDeltaTs[n][m]= basedTperi;   /* *dTcoord; in the source*/ 

  #ifdef DEBUG_STACKSLIDE_FNC
    printf("vecDeltaTs is %f\n",pTdotsAndDeltaTs->vecDeltaTs[n][m]);
  #endif
       }    
       
} 
LALFree(edat->ephemE);
LALFree(edat->ephemS);
LALFree(edat);

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

/*********************************************************************************/
/*              START function: FreeMem                                          */
/*********************************************************************************/
	
int FreeMem(void)
{
	if (SSparams->numSpinDown>0) {
	  for(i=0;i<SSparams->numSTKs;i++) {
     		LALFree(pTdotsAndDeltaTs->vecDeltaTs[i]);          
                                            }
  	  LALFree(pTdotsAndDeltaTs->vecDeltaTs);
	                              }
  	LALFree(pTdotsAndDeltaTs->vecTDots);
        LALFree(pTdotsAndDeltaTs);
	
	

	for (i=0; i<7 ; i++)
               {
	          LALFree(fakeSFT[i]);
               }
        LALFree(fakeSFT);
		
        LALFree(SSparams);
    
	LALFree(params);
        
	
return 0;

}

/*********************************************************************************/
/*              END function: FreeMem                                          */
/*********************************************************************************/

/*********************************************************************************/
/*              START function: SetupBaryInput                                   */
/*********************************************************************************/

int SetupBaryInput()
{
FILE *fpE;
FILE *fpS;

              /*allocate memory for the ephemeris*/
         edat=(EphemerisData *)LALMalloc(sizeof(EphemerisData));

         (*edat).ephiles.earthEphemeris="/home/re/data/LAL_tmp/earth00-04.dat";
         (*edat).ephiles.sunEphemeris="/home/re/data/LAL_tmp/sun00-04.dat";
         	 
	 /*Check that the ephemeris files exist*/
	 fpE=fopen((*edat).ephiles.earthEphemeris,"r");
         fpS=fopen((*edat).ephiles.sunEphemeris,"r");

	 if ((fpE==NULL) || (fpS==NULL))
	   {	 
	      fprintf(stdout,"impossible to find ephemeris files ");
           } else 
	 {
		 fprintf(stdout,"ephemeris file %s %s\n found!\n", (*edat).ephiles.sunEphemeris, (*edat).ephiles.earthEphemeris);
         }

	 fscanf(fpS,"%d \t %le \t %d\n",&x, &y, &z);
	 
	 /*fprintf(stdout,"x y z are %d \t %12.12le \t %d\n", x,y,z);*/

         /*Read in ephemeris files*/
        LALInitBarycenter(&status, edat);

#ifdef DEBUG_STACKSLIDE_FNC
	fprintf(stdout,"first entry %f\n",(*edat).ephemE[0].gps);
#endif
	
	printf("end of function SetupBaryInput\n");

 return 3;
  /* convert the observation start time at the detector into start time at the SSB */

	
}
/*********************************************************************************/
/*              END function: SetupBaryInput                                   */
/*********************************************************************************/

	
/*********************************************************************************/
/*              START function: GetSSBTime                                   */
/*********************************************************************************/

int GetSSBTime(LIGOTimeGPS *tdet, LIGOTimeGPS *tssb)
{
static LALStatus status;  
printf("start function GetSSBTime\n");
         baryinput.tgps.gpsSeconds=tdet->gpsSeconds;
	 baryinput.tgps.gpsNanoSeconds=tdet->gpsNanoSeconds;


	/*Get coords for detector LHO*/
       
         cachedDetector=lalCachedDetectors[LALDetectorIndexLHODIFF];
         baryinput.site.location[0]=cachedDetector.location[0]/LAL_C_SI;
         baryinput.site.location[1]=cachedDetector.location[1]/LAL_C_SI;
         baryinput.site.location[2]=cachedDetector.location[2]/LAL_C_SI;
         baryinput.alpha=0;
	 baryinput.delta=0;
	 baryinput.dInv=0.e0;
	 
	 LALBarycenterEarth(&status, &earth, tdet, edat);

	 REPORTSTATUS(&status);
LALBarycenter(&status, &emit, &baryinput, &earth);

tssb->gpsSeconds=emit.te.gpsSeconds;
tssb->gpsNanoSeconds=emit.te.gpsNanoSeconds;
printf("time at the SSB is : %d\n",tssb->gpsSeconds);
printf("end function GetSSBTime");
return 4;
}
	
/*********************************************************************************/
/*              END function: GetSSBTime                                   */
/*********************************************************************************/






