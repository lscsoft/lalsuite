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
 * \author Michael Landry, Gregory Mendell, Virginia Re
 */

/* REVISIONS: */
/* 01/28/04 gam; need to find index to which to SUM */
/* 02/02/04 gam; Further changes to handle case numSpinDown = 0 and numFreqDerivTotal = 0. */
/* 02/02/04 gam; Add preprocessor flag section. When flags defined will print debugging info; added/modified some printf statements for debugging. */
/* 04/26/04 gam; Change LALStackSlide to StackSlide and LALSTACKSLIDE to STACKSLIDE for initial entry to LALapps. */
/* 06/05/04 gam; Add gpsStartTimeSec and gpsStartTimeNan to StackSlideSkyParams; set these to epoch that gives T0 at SSB. */
/* 12/03/04 gam; Clean up indentation; remove extraneous or obsolete comments. */
/* 12/03/04 gam; Add parameter: BOOLEAN divideSUMsByNumSTKs. */
/* 12/06/04 gam; use refFreq = STK frequency in the middle of the SUM band as the reference to calculate how many bins to slide, rather than f0STK. */
/* 02/21/05 gam; propagate divideSUMsByNumSTKs and refFreq from StackSlideOld to new StackSlide function; */
/* 02/22/05 gam; Fix bug, kSUM not set correctly in StackSlide function; clean up indentation everywhere */
/* 04/12/05 gam; move variables not used by StackSlide function into StackSlideIsolated or StackSlideBinary */
/* 04/12/05 gam; Simplify StackSlideParams struct; change REAL8 **freqDerivData to REAL8 *freqDerivData; */
/* 05/06/05 gam; Add function SumStacks with just creates a SUM from the STKs without sliding */
/* 07/27/05 gam; Replace TimeToFloat and FloatToTime in StackSlideComputeSky with LAL functions used in StackSlideComputeSkyBinary */
/* 07/27/05 gam; In StackSlideComputeSky should ASSERT params->mObsSFT>0, not params->mObsSFT>=0. */
/* 07/27/05 gam; Add missing check of status pointer, CHECKSTATUSPTR(status), to StackSlideComputeSky */
/* 08/31/05 gam; In StackSlideComputeSky set ssbT0 to gpsStartTime; this now gives the epoch that defines T_0 at the SSB! */
  
/*********************************************/
/*                                           */
/* START SECTION: define preprocessor flags  */
/*                                           */
/*********************************************/
/* #define DEBUG_STACKSLIDE_FNC */
/* #define PRINT_STACKSLIDE_BINOFFSETS */
/* #define PRINT_STACKSLIDE_BINMISMATCH */
/* #define PRINT_STACKSLIDE_SPECIAL_DEBUGGING_INFO */
/* #define DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC*/
/* #define DEBUG_STACKSLIDECOMPUTESKY_FNC */
/*********************************************/
/*                                           */
/* END SECTION: define preprocessor flags    */
/*                                           */
/*********************************************/

/**
 *
 * ### Module \ref StackSlide.c ###
 *
 * Computes frequency model, slide stacks accordingly, and sum them.
 *
 * ### Prototypes ###
 *
 * <tt>StackSlideComputeSky()</tt>
 *
 * ### Description ###
 *
 * Slide 'em.
 *
 * ### Algorithm ###
 *
 * Sum them.
 *
 * ### Uses ###
 *
 */

#include "StackSlide.h"/*This includes isolated/binary info*/

/*********************************************************************************/
/*              START function: StackSlide                                       */
/*********************************************************************************/
void StackSlide(	LALStatus *status,
			REAL4FrequencySeries **SUMData,
			REAL4FrequencySeries **STKData,
                        TdotsAndDeltaTs *pTdotsAndDeltaTs,
			StackSlideParams *params)
{
       
  INT4 kSUM = 0; /* index gives which SUM */ /* 02/22/05 gam; HARD CODED TO 0; This function only returns one SUM! */
  INT4 iSUM = 0; /* index gives which SUM bin */  
  REAL8 f_t;

  REAL4 invNumSTKs = 1.0/((REAL4)params->numSTKs); /* 12/03/04 gam */ /* 02/21/05 gam */

  /* INT4 i,k,m, iSky; */ /* 02/22/05 gam; iSky is obsolete */
  INT4 i,k,m;  
  INT4 binoffset = 0;  
  
  #ifdef DEBUG_STACKSLIDE_FNC
    fprintf(stdout,"start function StackSlide\n");
    fflush(stdout);
  #endif   

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  INT4 iMinSTK = floor((params->f0SUM-params->f0STK)*params->tSTK + 0.5); /* Index of mimimum frequency to include when making SUMs from STKs */
  INT4 iMaxSTK = iMinSTK + params->nBinsPerSUM - 1;                       /* Index of maximum frequency to include when making SUMs from STKs */

  REAL8 refFreq = params->f0SUM + ((REAL8)(params->nBinsPerSUM/2))*params->dfSUM; /* 12/06/04 gam */ /* 02/21/05 gam */
  
  #ifdef DEBUG_STACKSLIDE_FNC
    fprintf(stdout,"invNumSTKs %g\n",invNumSTKs);
    fprintf(stdout,"params->divideSUMsByNumSTKs = %i \n",params->divideSUMsByNumSTKs);  
    fprintf(stdout,"refFreq %g\n",refFreq);
    fprintf(stdout,"iMaxSTK is %d\n",iMaxSTK);
    fprintf(stdout,"iMinSTK is %d\n",iMinSTK);
    fprintf(stdout,"numSTKs is %d\n", params->numSTKs);    
    fflush(stdout);
  #endif

  for(k=0;k<params->numSTKs;k++) {

        /* compute frequency */
        f_t = refFreq; /* 12/06/04 gam */ /* 02/21/05 gam */
	
        for (m=0; m<params->numSpinDown; m++) {
            f_t += params->freqDerivData[m] * pTdotsAndDeltaTs->vecDeltaTs[k][m];
        }
        f_t *= pTdotsAndDeltaTs->vecTDots[k]; /* 04/12/05 gam; *= same as f_t = f_t *... but faster */
        binoffset = floor(( (f_t - refFreq) * params->tSTK) + 0.5 ); /* 12/06/04 gam */ /* 02/21/05 gam */
        
        #ifdef PRINT_STACKSLIDE_BINOFFSETS
               fprintf(stdout, "In StackSlide for SFT #%i  binoffset = %i \n",k,binoffset);
               fflush(stdout);
        #endif

        #ifdef PRINT_STACKSLIDE_BINMISMATCH
	       /*choose one of the following*/
	       
 /*! 24/10*/        /*  fprintf(stdout, "%i binmismatch = %g \n",STKData[k]->epoch.gpsSeconds,fabs( ((f_t - refFreq) * params->tSTK) - (REAL8)binoffset ));*/
	 
	       /*      fprintf(stdout, "%i %g \n",STKData[k]->epoch.gpsSeconds,fabs( ((f_t - refFreq) * params->tSTK) - (REAL8)binoffset ));*/
   
	       fprintf(stdout, "%i %g %g %g %i\n",STKData[k]->epoch.gpsSeconds,fabs( f_t - refFreq),  pTdotsAndDeltaTs->vecTDots[k], refFreq , binoffset);

	       
           fflush(stdout);
        #endif
 
        #ifdef PRINT_STACKSLIDE_SPECIAL_DEBUGGING_INFO
            /* A comparison is made between the STK bin with maximum power
               and the STK power when demodulating for the SUM midpoint
               frequency, i.e., the bin for that freq + binoffset.
               For an injected signal into the middle of the SUM band
               these should agree if StackSlide is working properly. */
            {
               INT4 iPwrMax = -1;
               REAL4 pwrMax = 0.0;
               for(i=0;i<STKData[k]->data->length; i++) {
                 if (pwrMax < STKData[k]->data->data[i]) {
                   iPwrMax = i;
                   pwrMax = STKData[k]->data->data[i];
                 }
               }

               if (iPwrMax != iMinSTK + params->nBinsPerSUM/2+binoffset) {
                 fprintf(stdout, "Missed max pwr: SFT #%i, iPwrMax = %i, N/2+binoffset = %i, pwrMax = %g, STK(N/2+binoffset) = %g\n",k,iPwrMax,iMinSTK + params->nBinsPerSUM/2+binoffset,pwrMax,STKData[k]->data->data[iMinSTK + params->nBinsPerSUM/2+binoffset]);
                 fprintf(stdout, "for SFT #%i STK(N/2+binoffset-1) = %g \n",k,STKData[k]->data->data[iMinSTK + params->nBinsPerSUM/2+binoffset-1]);
                 fprintf(stdout, "for SFT #%i STK(N/2+binoffset)   = %g \n",k,STKData[k]->data->data[iMinSTK + params->nBinsPerSUM/2+binoffset]);
                 fprintf(stdout, "for SFT #%i STK(N/2+binoffset+1) = %g \n",k,STKData[k]->data->data[iMinSTK + params->nBinsPerSUM/2+binoffset+1]);
                 fprintf(stdout, "Compare for SFT #%i binoffset vs REAL8 binoffset: %i, %g \n",k,binoffset,((f_t - refFreq) * params->tSTK));
               } else {
                 fprintf(stdout, "Found max pwr: SFT #%i, iPwrMax = %i, N/2+binoffset = %i, pwrMax = %g, STK(N/2+binoffset) = %g\n",k,iPwrMax,iMinSTK + params->nBinsPerSUM/2+binoffset,pwrMax,STKData[k]->data->data[iMinSTK + params->nBinsPerSUM/2+binoffset]);
               }
               fflush(stdout);
            }
        #endif

        for (i=iMinSTK;i<=iMaxSTK; i++) {
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
      
  /* 12/03/04 gam; added params->divideSUMsByNumSTKs */ /* 02/21/05 gam */
  
  


if (params->divideSUMsByNumSTKs) {
     /* Normalize the SUMs with params->numSTKs*/
     for(i=0;i<params->nBinsPerSUM; i++) {
        /* SUMData[kSUM]->data->data[i] =  SUMData[kSUM]->data->data[i]/((REAL4)params->numSTKs); */
        SUMData[kSUM]->data->data[i] = SUMData[kSUM]->data->data[i]*invNumSTKs;  /* 12/03/04 gam; multiply by 1.0/numSTKs */
     }
  }

  #ifdef DEBUG_STACKSLIDE_FNC
    fprintf(stdout,"end function StackSlide\n");
    fflush(stdout);
  #endif
  
  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);
 }   /*end of void StackSlide*/
/*********************************************************************************/
/*              END function: StackSlide                                         */
/*********************************************************************************/

/* 05/06/05 gam; Add function SumStacks with just creates a SUM from the STKs without sliding */
/*********************************************************************************/
/*              START function: SumStacks                                        */
/*********************************************************************************/
void SumStacks( 	LALStatus *status,
			REAL4FrequencySeries **SUMData,
			REAL4FrequencySeries **STKData,
			StackSlideParams *params)
{  
  INT4 kSUM = 0; /* index gives which SUM     */ /* HARD CODED TO 0; This function only returns one SUM! */
  INT4 iSUM = 0; /* index gives which SUM bin */  

  REAL4 invNumSTKs = 1.0/((REAL4)params->numSTKs);

  INT4 i,k;

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  INT4 iMinSTK = floor((params->f0SUM-params->f0STK)*params->tSTK + 0.5); /* Index of mimimum frequency to include when making SUMs from STKs */
  INT4 iMaxSTK = iMinSTK + params->nBinsPerSUM - 1;                       /* Index of maximum frequency to include when making SUMs from STKs */

  for (k=0;k<params->numSTKs;k++) {

        for (i=iMinSTK;i<=iMaxSTK; i++) {
            iSUM = i - iMinSTK;
            if (k==0) {
               /* Starting a new SUM: initialize */
               SUMData[kSUM]->data->data[iSUM] = STKData[k]->data->data[i];
               SUMData[kSUM]->epoch.gpsSeconds = params->gpsStartTimeSec;
               SUMData[kSUM]->epoch.gpsNanoSeconds = params->gpsStartTimeNan;
               SUMData[kSUM]->f0=params->f0SUM;
               SUMData[kSUM]->deltaF=params->dfSUM;
               SUMData[kSUM]->data->length=params->nBinsPerSUM;
            } else {
               SUMData[kSUM]->data->data[iSUM] += STKData[k]->data->data[i];
            }
        }/*end of for iMinSTK*/

  } /* END for(k=0;k<params->numSTKs;k++) */
  /* Normalize the SUMs with params->numSTKs*/

  if (params->divideSUMsByNumSTKs) {
     /* Normalize the SUMs with params->numSTKs*/
     for(i=0;i<params->nBinsPerSUM; i++) {
        SUMData[kSUM]->data->data[i] = SUMData[kSUM]->data->data[i]*invNumSTKs;
     }
  }

  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);
 }   /*end of void SumStacks*/
/*********************************************************************************/
/*              END function: SumStacks                                          */
/*********************************************************************************/

/*********************************************************************************/
/*              START function: StackSlideComputeSky                             */
/*********************************************************************************/

/* Internal routines copied from ComputeSky.c */ 
/* static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps);
static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f); */ /* 07/27/05 gam */
 

void StackSlideComputeSky (LALStatus            *status,
			   TdotsAndDeltaTs      *pTdotsAndDeltaTs,
			   StackSlideSkyParams  *params)
{  
  
  INT4 m, n;
  /* REAL8 t; */                  /* 07/27/05 gam; this and next variable are no longer used */
  /* REAL8 tBary; */
  REAL8 dTbary;                   /* T - T0 */
  REAL8 basedTbary;               /* temporary container for (T - T0)^m */
  /* REAL8 tB0; */                /* 07/27/05 gam; use next four variable to compute T - T0 in SSB using LAL functions */
  LIGOTimeGPS ssbT0;              /* T0 in the Solar System BaryCenter */
  LALTimeInterval dTBaryInterval; /* T - T0 in the Solar System BaryCenter */
  LALTimeInterval HalfSFT;        /* timebase of SFT/2 as GPS interval  */
  REAL8 HalfSFTfloat;             /* timebase of SFT/2 as REAL8 number */

  INITSTATUS(status);
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
  /* ASSERT(params->mObsSFT>=0, status, STACKSLIDECOMPUTESKYH_ENEGA, STACKSLIDECOMPUTESKYH_MSGENEGA); */ /* 07/27/05 gam */
  ASSERT(params->mObsSFT>0, status, STACKSLIDECOMPUTESKYH_ENEGA, STACKSLIDECOMPUTESKYH_MSGENEGA);
  ASSERT(params->tSFT>=0, status, STACKSLIDECOMPUTESKYH_ENEGA, STACKSLIDECOMPUTESKYH_MSGENEGA);
 
  for(n=0;n<params->mObsSFT;n++) {
     ASSERT(params->tGPS[n].gpsSeconds>=0, status, STACKSLIDECOMPUTESKYH_ENEGA, STACKSLIDECOMPUTESKYH_MSGENEGA);
  }
 
  /* Check to make sure pointer to output is not NULL */
  ASSERT(pTdotsAndDeltaTs!=NULL, status, STACKSLIDECOMPUTESKYH_ENNUL, STACKSLIDECOMPUTESKYH_MSGENNUL);
 
  /* 07/27/05 gam; note that the gpsStartTime refers to the start of the epoch that defines the template spindown parameters. */
  /* Next lines of code find ssbT0 == T0, the SSB time that defines this epoch */  
  /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */  
  /* params->baryinput->tgps.gpsSeconds = params->gpsStartTimeSec;
  params->baryinput->tgps.gpsNanoSeconds = params->gpsStartTimeNan;
  LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat); CHECKSTATUSPTR(status);
  LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth); CHECKSTATUSPTR(status); */
  /* TimeToFloat(&tB0, &(params->emit->te)); */ /* 07/27/05 gam; replace tB0 with ssbT0 */
  /* ssbT0.gpsSeconds = params->emit->te.gpsSeconds;
  ssbT0.gpsNanoSeconds = params->emit->te.gpsNanoSeconds; */
  /* 08/31/05 gam; set ssbT0 to gpsStartTime; this now gives the epoch that defines T0 at the SSB! */
  ssbT0.gpsSeconds = ((INT4)params->gpsStartTimeSec);
  ssbT0.gpsNanoSeconds = ((INT4)params->gpsStartTimeNan);

  /* 07/27/05 gam; Find GPS interval the represent 1/2 the SFT time baseline */
  HalfSFTfloat=params->tSFT/2.0;
  LALFloatToInterval(status->statusPtr,&HalfSFT,&HalfSFTfloat); CHECKSTATUSPTR(status);

  /* Find T - T0, dT/dt, powers of T - T0 for midpoint time of each SFT; save in pTdotsAndDeltaTs struct */
  for (n=0; n<params->mObsSFT; n++) {     
     /* 07/27/05 gam; set params->baryinput->tgps to the midpoint time of the current SFT */
     /*t=(REAL8)(params->tGPS[n].gpsSeconds)+(REAL8)(params->tGPS[n].gpsNanoSeconds)*1.0E-9+0.5*params->tSFT;
     FloatToTime(&(params->baryinput->tgps), &t);*/ /* Instead use LALIncrementGPS: */
     LALIncrementGPS(status->statusPtr, &(params->baryinput->tgps),&(params->tGPS[n]),&HalfSFT);
          
     /* 07/27/05 gam; find T - T0 in SSB and tDot = dT/dt for midpoint time of current SFT */
     LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat); CHECKSTATUSPTR(status);
     LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth); CHECKSTATUSPTR(status);
     /* TimeToFloat(&tBary, &(params->emit->te));  
     dTbary = tBary-tB0; */ /* 07/27/05 gam; replace tB0 with ssbT0 */
     /* Subtract two LIGOTimeGPS to get LALTimeInterval; convert to REAL8 */
     LALDeltaGPS(status->statusPtr,&dTBaryInterval,&(params->emit->te),&ssbT0); CHECKSTATUSPTR(status);
     LALIntervalToFloat(status->statusPtr,&dTbary,&dTBaryInterval); CHECKSTATUSPTR(status);
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

/* static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps)
{
  INT4 x, y;

  x=tgps->gpsSeconds;
  y=tgps->gpsNanoSeconds;
  *f=(REAL8)x+(REAL8)y*1.e-9;
} */ /* 07/27/05 gam; obsolete */

/* static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f)
{
  REAL8 temp0, temp2, temp3;
  REAL8 temp1, temp4;
  
  temp0 = floor(*f); */    /* this is tgps.S */
/*  temp1 = (*f) * 1.e10;
  temp2 = fmod(temp1, 1.e10);
  temp3 = fmod(temp1, 1.e2); 
  temp4 = (temp2-temp3) * 0.1;

  tgps->gpsSeconds = (INT4)temp0;
  tgps->gpsNanoSeconds = (INT4)temp4;
} */ /* 07/27/05 gam; obsolete */

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

void StackSlideComputeSkyBinary( LALStatus 			*status, 
				 TdotsAndDeltaTs 		*pTdotsAndDeltaTs, 
				 INT8 				iSkyCoh,
				 /*StackSlideBinarySkyParams 	*params*/
				 StackSlideSkyParams 		*params
			        )

{
  INT4 m, n, nP;
  REAL8 dTbary;
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

  #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
    fprintf(stdout,"start function ComputeSkyBinary\n");
    fflush(stdout);   
  #endif

  INITSTATUS(status);
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
  #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
    fprintf(stdout,"leap %d\n",params->edat->leap); 
    fflush(stdout);   
  #endif

  a=params->SemiMajorAxis;  /* This is the projected semi-major axis of the orbit normalised by the speed of light */
  Period=params->OrbitalPeriod;  /* This is the period of the orbit in seconds */
  ecc=params->OrbitalEccentricity;  /* This is the eccentricity of the orbit */
  parg=params->ArgPeriapse;  /* This is the argument of periapse defining the angular location of the source at periapsis */
                            /* measured relative to the ascending node */
  Tperi.gpsSeconds=params->TperiapseSSB.gpsSeconds;  /* This is the GPS time as measured in the SSB of the observed */
  Tperi.gpsNanoSeconds=params->TperiapseSSB.gpsNanoSeconds;  /* periapse passage of the source */

  /* Convert half the SFT length to a LALTimeInterval for later use */
  HalfSFTfloat=params->tSFT/2.0;
  #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
    fprintf(stdout,"HalfSFTfl %f\n",HalfSFTfloat);
    fflush(stdout);   
  #endif 

  LALFloatToInterval(status->statusPtr,&HalfSFT,&HalfSFTfloat);

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
     /*  LALIncrementGPS(status->statusPtr,&(params->baryinput->tgps),&tGPS[n],&HalfSFT);*/

     LALIncrementGPS(status->statusPtr, &(params->baryinput->tgps),&params->tGPS[n],&HalfSFT);
 
     #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
       fprintf(stdout,"tgps %d\n",params->baryinput->tgps.gpsSeconds);
       fprintf(stdout,"leap %d\n",params->edat->leap);
       fflush(stdout);   
     #endif

     /* Convert this mid point detector time into barycentric time (SSB) */

     /* GetSSBTime(&(baryinput.tgps), &TmidSSB) ;*/

     LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat);    
     LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);   

     CHECKSTATUSPTR(status);
     #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
       fprintf(stdout,"status %d\n",status->statusCode); 
       fflush(stdout);   
     #endif

     /* LALBarycenterEarth(status->statusPtr, csParams->earth, &(csParams->baryinput->tgps), csParams->edat);    
        LALBarycenter(status->statusPtr, csParams->emit, csParams->baryinput, csParams->earth);*/       

     /* Calculate the time difference since the observed periapse passage in barycentric time (SSB). */ 
     /* This time difference, when converted to REAL8, should lose no precision unless we are dealing */
     /* with periods >~ 1 Year */
     
     LALDeltaGPS(status->statusPtr,&dTBaryInterval,&(params->emit->te),&Tperi); 
               
     LALIntervalToFloat(status->statusPtr,&dTbary,&dTBaryInterval);

     /* LALDeltaGPS(status->statusPtr,&dTBaryInterval,&(csParams->emit->te),&Tperi); 
        LALIntervalToFloat(status->statusPtr,&dTbary,&dTBaryInterval);*/

     /* Calculate the time since the last periapse passage ( < Single Period (SP) ) */
     dTbarySP=Period*((dTbary/(1.0*Period))-(REAL8)floor(dTbary/(1.0*Period)));
     
     /* Calculate number of full orbits completed since the input observed periapse passage */
     nP=(INT4)floor(dTbary/(1.0*Period));
     
     #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
       fprintf(stdout,"dTbary is %f\n",dTbary);
       fprintf(stdout,"dTbarySP is %f\n",dTbarySP);
       fprintf(stdout,"nP is %i\n",nP);
       fflush(stdout);   
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

     #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
       fprintf(stdout,"dTperi is %f\n",dTperi);     
       fflush(stdout);   
     #endif

     /* The following quantity is the derivative of the time coordinate measured at the source with */
     /* respect to the time coordinate measured in the SSB : dt_(source)/dt_(SSB) */
     dTcoord=(1.0-(ecc*cos(E)))/(1.0+(p*cos(E))-(q*sin(E)));  

     /* The following quantity is the derivitive of the time coordinate measured in the SSB with */
     /* respect to the time coordinate measured at the chosen detector.  It was calculated via the */
     /* last call to LALBarycenter : dt_(SSB)/dt_(detector)  */
     Tdotbin = params->emit->tDot*dTcoord; /*=dt_source/dt_detector puts you back into the source*/ 

     /*fprintf(stdout,"%g %g\n", params->emit->tDot, dTcoord);*/

     #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
        fprintf(stdout,"tDot is %f\n",params->emit->tDot); 
       fprintf(stdout,"dTcoord is %f\n",dTcoord);     
       fprintf(stdout,"Tdotbin is %f\n",Tdotbin);
       fflush(stdout);   
     #endif

     /*    Tdotbin = csParams->emit->tDot*dTcoord; */
     
     /* Loop over all spin down orders plus 0th order (f0) */
     /* In this loop we calculate the SkyConstants defined in the documentation as A_{s,alpha} and B_{s,alpha} */

      
     /*!*/ pTdotsAndDeltaTs->vecTDots[n]= Tdotbin;
     
     #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
       fprintf(stdout,"pTdotsAndDeltaTs->vecTDots is %f\n",pTdotsAndDeltaTs->vecTDots[n]);
       fflush(stdout);   
     #endif

     for (m=0; m<params->spinDwnOrder; m++)
     {
       /* raise the quantity dTperi to the power m */
       basedTperi = pow(dTperi, (REAL8)m+1);

       
       #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
         fprintf(stdout,"basedTperi %f\n",basedTperi);
         fflush(stdout);
       #endif         

       /*!*/ /*the 2 lines below must be changed */

       /* Calculate A coefficients */
       /* skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m]=1.0/((REAL8)m+1.0)*basedTperi*dTperi-0.5*params->tSFT*basedTperi*Tdotbin;*/
       /* Calculate B coefficients */
       /* skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m+1]= params->tSFT*basedTperi*Tdotbin;*/

       /*expressing the time difference in the SSB */
       /*!*/ /*pTdotsAndDeltaTs->vecDeltaTs[m][n]= basedTperi*pTdotsAndDeltaTs->vecTDots[n]; */
       pTdotsAndDeltaTs->vecDeltaTs[n][m]= basedTperi; /*in the Source*/ 

       #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
         fprintf(stdout,"vecDeltaTs %f\n", pTdotsAndDeltaTs->vecDeltaTs[n][m]);
         fflush(stdout);
       #endif
    
     } /* END for (m=0; m<params->spinDwnOrder; m++) */
    
     #ifdef DEBUG_STACKSLIDECOMPUTESKYBINARY_FNC
       fprintf(stdout,"end function StackSlideComputeSkyBinary\n");
       fflush(stdout);
     #endif
  } /* for (n=0; n<params->mObsSFT; n++) */
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
  INITSTATUS(status);
  ASSERT(tr0,status, 1, "Null pointer");

  /* this is the function relating the observed time since periapse in the SSB to the true eccentric anomoly E */

  *tr = *(REAL8 *)tr0*(-1.0) + (Period/LAL_TWOPI)*(lE+(p*sin(lE))+q*(cos(lE)-1.0));
  RETURN(status);
}

/*********************************************************************************/
/*              END function: StackSlideComputeSkyBinary                         */
/*********************************************************************************/
