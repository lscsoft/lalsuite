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
Computes the phase model coefficients necessary for a successful demodulation.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{StackSlideCP}
\index{\texttt{ComputeSky()}}

\subsubsection*{Description}
Slide 'em.

\subsubsection*{Algorithm}
Shift 'em.

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StackSlideCV}}

</lalLaTeX> */

#include "StackSlide.h"

/* NRCSID( STACKSLIDEC, "$Id$"); */
NRCSID( STACKSLIDEC,  "$Id$");
/* <lalVerbatim file="ComputeSkyCP"> */

void StackSlide(	LALStatus *status, 
			REAL4FrequencySeries **SUMData, 
			REAL4FrequencySeries **STKData, 
			StackSlideParams *params)
{
    /* Declare and variables needed for making SUMs */
    /* UINT4 iSky = 0;  */ /* index gives which Sky Position */
    /* UINT4 iFreqDeriv = 0; */  /* index give which Freq Derivs  */
    /* UINT4 kSUM = 0;   */     /* index gives which SUM */
    /* UINT4 iSUM = 0;   */     /* index gives which SUM bin */
    
    INT4 iSky = 0;       /* index gives which Sky Position */
    INT4 iFreqDeriv = 0; /* index give which Freq Derivs  */
    INT4 kSUM = 0;       /* index gives which SUM */
    INT4 iSUM = 0;       /* index gives which SUM bin */
    
    INT4 numLoopOnSpindown; /* 1 if numSpindown is 0, else = params->numFreqDerivTotal */
    
    INT4 i,k,m;
    INT4 binoffset = 0;
	
	StackSlideSkyParams *csParams;           /* Container for ComputeSky */
	/* INT4 iSkyCoh;                  For ComputeSky */
	
	TdotsAndDeltaTs *pTdotsAndDeltaTs;

	EarthState earth;
	EmissionTime emit;

	REAL8 f_t;
	
    /* UINT4 iMinSTK = floor((params->f0SUM-params->f0STK)*params->tSTK + 0.5); */ /* Index of mimimum frequency to include when making SUMs from STKs */
    /* UINT4 iMaxSTK = iMinSTK + params->nBinsPerSUM - 1;                       */ /* Index of maximum frequency to include when making SUMs from STKs */
    
    INT4 iMinSTK = floor((params->f0SUM-params->f0STK)*params->tSTK + 0.5); /* Index of mimimum frequency to include when making SUMs from STKs */
    INT4 iMaxSTK = iMinSTK + params->nBinsPerSUM - 1;                       /* Index of maximum frequency to include when making SUMs from STKs */

INITSTATUS (status, "StackSlide", STACKSLIDEC); 
ATTATCHSTATUSPTR(status);
	
  	/* Allocate space and set quantities for call to ComputeSky() */
	csParams=(StackSlideSkyParams *)LALMalloc(sizeof(StackSlideSkyParams));
	csParams->skyPos=(REAL8 *)LALMalloc(2*sizeof(REAL8));
  	pTdotsAndDeltaTs=(TdotsAndDeltaTs *)LALMalloc(sizeof(TdotsAndDeltaTs));
  	pTdotsAndDeltaTs->vecTDots  = (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSTKs);
  	
	if (params->numSpinDown>0) {
	  pTdotsAndDeltaTs->vecDeltaTs  = (REAL8 **)LALMalloc(sizeof(REAL8 *)*params->numSTKs);          
	  for(i=0;i<params->numSTKs;i++)
              {
       		pTdotsAndDeltaTs->vecDeltaTs[i]  = (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSpinDown);          
              }
	    }

 	/* if no_spindowns > 0, only pTdotsAndDeltaTs->vecDeltaTs= (REAL8 *)LALMalloc(sizeof(REAL8)*params->numSTKs); 
	i.e. if no spindowns then you don't need to allocate the second of the two lines above
	*/
        /* 02/02/04 gam; moved next 5 lines from inside loops: */
	if (params->numFreqDerivTotal==0) {
	   numLoopOnSpindown = 1;  /* Even if no spindown, need to execute iFreqDeriv loop once to handle case of zero spindown */
	} else {
	   numLoopOnSpindown = params->numFreqDerivTotal;
	}

  for(iSky=0;iSky<params->numSkyPosTotal;iSky++)
  	{
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

	 /* Index for sky patch : remove as not used */
	/* iSkyCoh = params->patchNumber; */

	/* what the hell does this number mean? */
        for(i=0;i<params->numSTKs;i++)
         {
      	  pTdotsAndDeltaTs->vecTDots[i] = 0.0; /* Initialize */
          for(m=0;m<params->numSpinDown;m++)
            {
      	     pTdotsAndDeltaTs->vecDeltaTs[i][m] = 0.0; /* Initialize */
            }
        }

	/* Call STACKSLIDECOMPUTESKY() */
	/* StackSlideComputeSky(status->statusPtr, pTdotsAndDeltaTs, iSkyCoh, csParams); */
	
	StackSlideComputeSky(status->statusPtr, pTdotsAndDeltaTs, csParams);
	/* ComputeSky(&status, pTdotsAndDeltaTs, iSkyCoh, csParams); */
    	CHECKSTATUSPTR (status);

	
	/* as we must execute the for loop below, even if there are no frequency derivatives, 
	   set numLoopOnSpindown, which governs how many times the for loop is executed */
	/* if (params->numFreqDerivTotal==0) {
	   numLoopOnSpindown = 1;
	   } else {
	   numLoopOnSpindown = params->numFreqDerivTotal;
 	   } */ /* 02/02/04 gam; moved above loops. */
	
	/* for all spindowns, loop */
 	for(iFreqDeriv=0;iFreqDeriv<numLoopOnSpindown;iFreqDeriv++)
  	{
 	/* call computesky for all stacks at once */
 	/* compute offset  for each stack */
	
	kSUM = iSky*numLoopOnSpindown + iFreqDeriv; /* 01/28/04 gam; need to find index to which to SUM */
 
      for(k=0;k<params->numSTKs;k++) {
      	    
      /* compute frequency */
	/* f_t = params->f0STK * params->tSTK; */
	f_t = params->f0STK;
	for (m=0; m<params->numSpinDown; m++) 
       	   {
		/* f_t += ( params->freqDerivData[iFreqDeriv][m] * pTdotsAndDeltaTs[2*k*(params->numSpinDown)+2*(INT4)m+1]); */
		f_t += params->freqDerivData[iFreqDeriv][m] * pTdotsAndDeltaTs->vecDeltaTs[k][m];
       	   }	   
	f_t = f_t * pTdotsAndDeltaTs->vecTDots[k];
	
	   
	   binoffset = floor(( (f_t - params->f0STK) * params->tSTK) + 0.5 );
 
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
            }
      } /* END for(k=0;k<params->numSTKs;k++) */
      /* Normalize the SUMs with params->numSTKs*/
      for(i=0;i<params->nBinsPerSUM; i++) {
               SUMData[kSUM]->data->data[i] =  SUMData[kSUM]->data->data[i]/((REAL4)params->numSTKs);
      	}	
      } /* for iFreqDeriv = 0 to numFreqDerivTotal */
      } /* for iSky = `0 to numSkyPosTotal */
      	
	/* Deallocate space for ComputeSky parameters */
	LALFree(csParams->skyPos);
	LALFree(csParams);

      /* Deallocate memory */
        /* 02/02/04 gam; pTdotsAndDeltaTs->vecDeltaTs is allocated only if numSpinDown > 0. */      
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
 }
 
/* Internal routines copied from ComputeSky.c */ 
static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps);
static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f);

/* <lalVerbatim file="ComputeSkyCP"> */
void StackSlideComputeSky	(LALStatus		*status, 
				 TdotsAndDeltaTs 	*pTdotsAndDeltaTs, 
				 /* INT8 		iSkyCoh, */
				 StackSlideSkyParams 	*params)
{  /* </lalVerbatim> */
  
  INT4	m, n;
	REAL8	t;
	REAL8	basedTbary;
	
	REAL8	dTbary;
	REAL8	tBary;
	REAL8	tB0;
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
 
 for(n=0;n<params->mObsSFT;n++)
   {
     ASSERT(params->tGPS[n].gpsSeconds>=0, status, STACKSLIDECOMPUTESKYH_ENEGA, 	STACKSLIDECOMPUTESKYH_MSGENEGA);
   }
 
 /* Check to make sure pointer to output is not NULL */
 ASSERT(pTdotsAndDeltaTs!=NULL, status, STACKSLIDECOMPUTESKYH_ENNUL, STACKSLIDECOMPUTESKYH_MSGENNUL);
 
 /*params->baryinput->tgps.gpsSeconds=params->tGPS[0].gpsSeconds;
 params->baryinput->tgps.gpsNanoSeconds=params->tGPS[0].gpsNanoSeconds; */ /* 06/05/04 gam */
 params->baryinput->tgps.gpsSeconds = params->gpsStartTimeSec; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
 params->baryinput->tgps.gpsNanoSeconds = params->gpsStartTimeNan; /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
 
 /* params->baryinput->alpha=params->skyPos[iSkyCoh];
 params->baryinput->delta=params->skyPos[iSkyCoh+1]; */
 
 LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat);
  
 LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);
 
 TimeToFloat(&tB0, &(params->emit->te));
 for (n=0; n<params->mObsSFT; n++) 
   {
     t=(REAL8)(params->tGPS[n].gpsSeconds)+(REAL8)(params->tGPS[n].gpsNanoSeconds)*1.0E-9+0.5*params->tSFT; 
     
     FloatToTime(&(params->baryinput->tgps), &t);
	
     LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat);
     LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);
     
     TimeToFloat(&tBary, &(params->emit->te));

     dTbary = tBary-tB0;		
     
     pTdotsAndDeltaTs->vecTDots[n]= params->emit->tDot;
     
     /* for (m=0; m<params->spinDwnOrder+1; m++) */
      for (m=0; m<params->spinDwnOrder; m++) 
       {
	 basedTbary = pow(dTbary, (REAL8)m+1);
	 /* pTdotsAndDeltaTs[2*n*(params->spinDwnOrder+1)+2*(INT4)m]=1.0/((REAL8)m+1.0)*basedTbary*dTbary-0.5*params->tSFT*params->emit->tDot*basedTbary;
	    pTdotsAndDeltaTs[2*n*(params->spinDwnOrder+1)+2*(INT4)m+1]= params->tSFT*params->emit->tDot*basedTbary; */
	 
	 pTdotsAndDeltaTs->vecDeltaTs[n][m]= basedTbary;
        }
   } 
 /* Normal Exit */
 DETATCHSTATUSPTR(status);
 RETURN(status);
}

/* More internal routines: these from ComputeSky.c */
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
