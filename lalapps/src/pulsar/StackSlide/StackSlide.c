/*********************************** <lalVerbatim file="StackSlideCV">
Author:  Landry, M., and Mendell, G.
$Id$
************************************ </lalVerbatim> */

/* REVISIONS: */
/* 01/28/04 gam; need to find index to which to SUM */
/* 02/02/04 gam; Further changes to handle case numSpinDown = 0 and numFreqDerivTotal = 0. */
/* 02/02/04 gam; Add preprocessor flag section. When flags defined will print debugging info; added/modified some printf statements for debugging. */
/* 04/26/04 gam; Change LALStackSlide to StackSlide and LALSTACKSLIDE to STACKSLIDE for initial entry to LALapps. */
/* 06/05/04 gam; Add gpsStartTimeSec and gpsStartTimeNan to StackSlideSkyParams; set these to epoch that gives T0 at SSB. */
/* 12/03/04 gam; Clean up indentation; remove extraneous or obsolete comments. */
/* 12/03/04 gam; Add parameter: BOOLEAN divideSUMsByNumSTKs. */
/* 12/06/04 gam; use refFreq = STK frequency in the middle of the SUM band as the reference to calculate how many bins to slide, rather than f0STK. */

/*********************************************/
/*                                           */
/* START SECTION: define preprocessor flags  */
/*                                           */
/*********************************************/
/* #define DEBUG_STACKSLIDE_FNC */
/* #define DEBUG_STACKSLIDECOMPUTESKY_FNC */
/*********************************************/
/*                                           */
/* END SECTION: define preprocessor flags    */
/*                                           */
/*********************************************/

/* <lalLaTeX> 
\subsection{Module \texttt{StackSlide.c}}\label{ss:StackSlide.c}
Computes frequency model, slide stacks accordingly, and sum them.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{StackSlideCP}
\index{\texttt{StackSlideComputeSky()}}

\subsubsection*{Description}
Slide 'em.

\subsubsection*{Algorithm}
Sum them.

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StackSlideCV}}

</lalLaTeX> */

#include "StackSlide.h"/*This includes isolated/binary info*/


NRCSID( STACKSLIDEC,  "$Id$");

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
    
    INITSTATUS (status, "StackSlide", STACKSLIDEC); 
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

               #ifdef DEBUG_STACKSLIDE_FNC
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
        	
       INITSTATUS (status, "StackSlide", STACKSLIDEC); 
       ATTATCHSTATUSPTR(status);

  
	INT4 iMinSTK = floor((params->f0SUM-params->f0STK)*params->tSTK + 0.5); /* Index of mimimum frequency                                                                       to include when making SUMs from STKs */
	
        INT4 iMaxSTK = iMinSTK + params->nBinsPerSUM - 1;                       /* Index of maximum frequency                                                                       to include when making SUMs from STKs */
  printf("iFreqDeriv is %d\n",params->iFreqDeriv);

  printf("iMaxSTK is %d\n", iMaxSTK);
 	 	

  kSUM = iSky*numLoopOnSpindown + params->iFreqDeriv; /* 01/28/04 gam; need to find index to which to SUM */
 /*05/02/18 vir kSUM=iSky*(numLoopOnSpinDown*params->nMaxSMA*params->nMaxTperi)+...*/
  /*... + params->iFreqDeriv*(params->nMaxSMA*params->nMaxTperi)+params->iSMA*(params->nMaxTperi) + params->iT;*/
  
  printf("iMinSTK is %d\n",iMinSTK);
  
  for(k=0;k<params->numSTKs;k++) {
  
	  printf("numSTKs is %d\n", params->numSTKs);
                               
	                               /* compute frequency */
	
	f_t = params->f0STK; /*this should be refFreq;*/
	printf("start ft is %f\n",f_t);
	for (m=0; m<params->numSpinDown; m++) 
       	   {
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
  /*    printf("Normalized sum is %f\n",SUMData[kSUM]->data->data[i]);*/
                           	          }	
      /*05/02/18 vir:*/ params->ParamsSMA[kSUM]=params->SemiMajorAxis;
      
      /*05/02/18 vir: params->ParamsTperi[kSUM]=params->TperiapseSSB*/
       printf("end function StackSlide\n");     
	
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);
 }   /*end of void StackSlide*/


/*********************************************************************************/
/*              END function: StackSlide                                         */
/*********************************************************************************/



/*********************************************************************************/
/*              START function: StackSlideComputeSky                             */
/*********************************************************************************/



/* Internal routines copied from ComputeSky.c */ 
static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps);
static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f);
 
/* <lalVerbatim file="StackSlideCP"> */
void StackSlideComputeSky (LALStatus            *status,
			   TdotsAndDeltaTs      *pTdotsAndDeltaTs,
			   StackSlideSkyParams  *params)
{  /* </lalVerbatim> */
  
  INT4 m, n;
  REAL8 t;
  REAL8 basedTbary;

  REAL8 dTbary;
  REAL8 tBary;
  REAL8 tB0;

  INITSTATUS (status, "StackSlideComputeSky", STACKSLIDEC);
  #ifdef DEBUG_STACKSLIDECOMPUTESKY_FNC
     fprintf(stdout, "Start StackSlideComputeSky\n"); 
     fflush(stdout);
  #endif  
  ATTATCHSTATUSPTR(status);
 
  /* Check for non-negativity of sky positions in SkyCoh[] */
  /* removed this, so is there another check we should do on sky position validity? */
  /* ASSERT(iSkyCoh>=0, status, STACKSLIDECOMPUTESKYH_ENEGA, STACKSLIDECOMPUTESKYH_MSGENEGA); */
 
  /* Check to make sure sky positions are loaded */
  ASSERT(params->skyPos!=NULL, status, STACKSLIDECOMPUTESKYH_ENULL, STACKSLIDECOMPUTESKYH_MSGENULL);
  ASSERT(params->skyPos!=NULL, status, STACKSLIDECOMPUTESKYH_ENULL, STACKSLIDECOMPUTESKYH_MSGENULL);
 
  /* Check to make sure parameters are loaded and reasonable */
  ASSERT(params->spinDwnOrder>=0, status, STACKSLIDECOMPUTESKYH_ENEGA, STACKSLIDECOMPUTESKYH_MSGENEGA);
  ASSERT(params->mObsSFT>=0, status, STACKSLIDECOMPUTESKYH_ENEGA, STACKSLIDECOMPUTESKYH_MSGENEGA);
  ASSERT(params->tSFT>=0, status, STACKSLIDECOMPUTESKYH_ENEGA, STACKSLIDECOMPUTESKYH_MSGENEGA);
 
  for(n=0;n<params->mObsSFT;n++) {
     ASSERT(params->tGPS[n].gpsSeconds>=0, status, STACKSLIDECOMPUTESKYH_ENEGA, STACKSLIDECOMPUTESKYH_MSGENEGA);
  }
 
  /* Check to make sure pointer to output is not NULL */
  ASSERT(pTdotsAndDeltaTs!=NULL, status, STACKSLIDECOMPUTESKYH_ENNUL, STACKSLIDECOMPUTESKYH_MSGENNUL);
 
  params->baryinput->tgps.gpsSeconds = params->gpsStartTimeSec; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
  params->baryinput->tgps.gpsNanoSeconds = params->gpsStartTimeNan; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
 
  LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat);
  
  LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);
 
  TimeToFloat(&tB0, &(params->emit->te));
  for (n=0; n<params->mObsSFT; n++) {
     t=(REAL8)(params->tGPS[n].gpsSeconds)+(REAL8)(params->tGPS[n].gpsNanoSeconds)*1.0E-9+0.5*params->tSFT; 
     
     FloatToTime(&(params->baryinput->tgps), &t);

     LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat);
     LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);
     
     TimeToFloat(&tBary, &(params->emit->te));

     dTbary = tBary-tB0;		
     
     pTdotsAndDeltaTs->vecTDots[n]= params->emit->tDot;
     
     /* for (m=0; m<params->spinDwnOrder+1; m++) */
      for (m=0; m<params->spinDwnOrder; m++) {
         basedTbary = pow(dTbary, (REAL8)m+1);
         pTdotsAndDeltaTs->vecDeltaTs[n][m]= basedTbary;
      }
  }
  
  /* Normal Exit */
  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* END StackSlideComputeSky() */

static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps)
{
  INT4 x, y;

  x=tgps->gpsSeconds;
  y=tgps->gpsNanoSeconds;
  *f=(REAL8)x+(REAL8)y*1.e-9;
}

static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f)
{
  REAL8 temp0, temp2, temp3;
  REAL8 temp1, temp4;
  
  temp0 = floor(*f);     /* this is tgps.S */
  temp1 = (*f) * 1.e10;
  temp2 = fmod(temp1, 1.e10);
  temp3 = fmod(temp1, 1.e2); 
  temp4 = (temp2-temp3) * 0.1;

  tgps->gpsSeconds = (INT4)temp0;
  tgps->gpsNanoSeconds = (INT4)temp4;
}
/*********************************************************************************/
/*              END function: StackSlideComputeSkY                               */
/*********************************************************************************/



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
	


