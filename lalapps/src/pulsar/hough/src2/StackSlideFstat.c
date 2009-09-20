/*
 *  Copyright (C) 2005 Gregory Mendell
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
 * \author Gregory Mendell
 * \file StackSlideFstat.c
 * \brief Module with functions that StackSlide a vector of Fstat values or any REAL8FrequencySeriesVector.
*/

/* define preprocessor flags:  */
/* #define PRINT_STACKSLIDE_BINOFFSETS */
/* Uncomment the next flag when using a threshold to only accept local maximum that stand high enough above neighboring bins */
/* #define INCLUDE_EXTRA_STACKSLIDE_THREHOLDCHECK */
/* THRESHOLDFRACSS is the fraction of power that must be in the local maxima when INCLUDE_EXTRA_STACKSLIDE_THREHOLDCHECK is defined */
#define THRESHOLDFRACSS 0.667

/* include files: */
#include "./StackSlideFstat.h"

RCSID( "$Id$");

#define TRUE (1==1)
#define FALSE (1==0)

#define SSMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define SSMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))

extern int lalDebugLevel;

#define BLOCKSIZE_REALLOC 50

static int smallerStackSlide(const void *a,const void *b) {
  SemiCohCandidate a1, b1;
  a1 = *((const SemiCohCandidate *)a);
  b1 = *((const SemiCohCandidate *)b);
  
  if( a1.significance < b1.significance )
    return(1);
  else if( a1.significance > b1.significance)
    return(-1);
  else
    return(0);
}

/** \brief Function StackSlides a vector of Fstat frequency series or any REAL8FrequencySeriesVector.
    \param out SemiCohCandidateList is a list of candidates
    \param vecF is a vector of Fstat frequency series or any REAL8FrequencySeriesVector.
    \param params is a pointer to SemiCoherentParams
    \out SemiCohCandidateList is a list of candidates
*/
void StackSlideVecF(LALStatus *status,
                    SemiCohCandidateList  *out,        /* output candidates */
                    REAL4FrequencySeriesVector *vecF,  /* vector with Fstat values or any REAL8FrequencySeriesVector */
                    SemiCoherentParams *params)        /* input parameters  */
{
  REAL8FrequencySeries stackslideSum;  /* The output of StackSliding the vecF values */
  
  REAL8 *pstackslideData;  /* temporary pointer */
  REAL4 *pFVecData;        /* temporary pointer */

  REAL8 pixelFactor;
  REAL8 alphaStart, alphaEnd, dAlpha, thisAlpha;
  REAL8 deltaStart, deltaEnd, dDelta, thisDelta;
  REAL8 fdotStart, fdotEnd, dfdot, thisFdot;
  UINT4 ialpha,nalpha,idelta,ndelta,ifdot,nfdot;
  UINT2 numSpindown;

  /* INT4 fBinIni, fBinFin, nSearchBins, nSearchBinsm1; */ 
  INT4 fBinIni, fBinFin, nSearchBins, nStackBinsm1; /* 12/14/06 gm; now account for extra bins in stack  */
  INT4 j, k, nStacks, offset, offsetj;
  INT4 extraBinsFstat, halfExtraBinsFstat; /* 12/14/06 gm; new parameter and half this parameter */
  UINT4 uk;
  REAL8 f0, deltaF, tEffSTK, fmid, alpha, delta;
  /* REAL8 patchSizeX, patchSizeY, f1jump; */  /* 12/14/06 gm; f1jump no longer needed */
  REAL8 patchSizeX, patchSizeY;
  REAL8VectorSequence *vel, *pos;
  REAL8 fdot, refTime;
  LIGOTimeGPS refTimeGPS;
  LIGOTimeGPSVector   *tsMid;
  REAL8Vector *timeDiffV=NULL;
  REAL8Vector *weightsV=NULL;
  REAL8 thisWeight;
  REAL8 threshold;
  
  toplist_t *stackslideToplist;

  PulsarDopplerParams outputPoint, outputPointUnc, inputPoint;

  /* Add error checking here: */
  ASSERT ( vecF != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );

  INITSTATUS( status, "StackSlideVecF", rcsid );
  ATTATCHSTATUSPTR (status);

  /* create toplist of candidates */
  if (params->useToplist) {
    create_toplist(&stackslideToplist, out->length, sizeof(SemiCohCandidate), smallerStackSlide);
  }

  /* copy some params to local variables */
  pixelFactor = params->pixelFactor;   
  nStacks = vecF->length;
  extraBinsFstat = (INT4)params->extraBinsFstat;      /* 12/14/06 gm; new parameter */
  halfExtraBinsFstat = extraBinsFstat/2;              /* 12/14/06 gm; half the new parameter */
  nSearchBins = vecF->data->data->length;
  nSearchBins -= extraBinsFstat;                      /* 12/14/06 gm; subtract extra bins */
  /* nSearchBinsm1 = nSearchBins - 1; */
  nStackBinsm1  = vecF->data->data->length - 1;       /* 12/14/06 gm; need to know this for the entire stack */
  deltaF = vecF->data[0].deltaF;
  tEffSTK = 1.0/deltaF; /*  Effective time baseline a stack */
  /* fBinIni = floor( f0*tEffSTK + 0.5);
  fBinFin = fBinIni + nSearchBins - 1; */
  fBinIni = floor( vecF->data[0].f0*tEffSTK + 0.5);   /* 12/14/06 gm; defined on the full stack band */
  fBinFin = fBinIni + vecF->data->data->length - 1;   /* 12/14/06 gm; defined on the full stack band */
  /* fmid = vecF->data[0].f0 + ((REAL8)(nSearchBins/2))*deltaF;*/
  /* f0 = vecF->data[0].f0; */  
  fmid = ((REAL8)((fBinIni+fBinFin)/2))*deltaF;       /* 12/14/06 gm; defined on the full stack band */
  f0 = ((REAL8)(fBinIni+halfExtraBinsFstat))*deltaF;  /* 12/14/06 gm; start frequency of stackslide output */
  weightsV=params->weightsV; /* Needs to be NULL or normalized stack weights */
  threshold = params->threshold;
    
  numSpindown = 1;       /* Current search is over 1 spindown value, df/dt */
  nfdot = params->nfdot; /* Number of df/dt values to search over */
  alpha = params->alpha;
  delta = params->delta;
  vel = params->vel;
  pos = params->pos;
  fdot = params->fdot;
  tsMid = params->tsMid;
  refTimeGPS = params->refTime;  
  refTime = XLALGPSGetREAL8(&refTimeGPS);

  /* allocate memory for StackSlide of Fvec values */
  stackslideSum.epoch = vecF->data[0].epoch; /* just use the epoch of the first stack */
  stackslideSum.deltaF = deltaF;
  stackslideSum.f0 = f0;
  stackslideSum.data = (REAL8Sequence *)LALCalloc( 1, sizeof(REAL8Sequence));
  stackslideSum.data->length = nSearchBins;
  stackslideSum.data->data = (REAL8 *)LALCalloc( 1, nSearchBins * sizeof(REAL8));
  pstackslideData = stackslideSum.data->data;
  
  /* set patch size */
  /* this is supposed to be the "educated guess" 
     delta theta = 1.0 / (Tcoh * f0 * Vepi )
     where Tcoh is coherent time baseline, 
     f0 is frequency and Vepi is rotational velocity 
     of detector */
  patchSizeX = params->patchSizeX;
  patchSizeY = params->patchSizeY;
  /* if patchsize is negative, then set default value */
  if ( patchSizeX < 0 )
    patchSizeX = 0.5 / ( fBinFin * VEPI ); 

  if ( patchSizeY < 0 )
    patchSizeY = 0.5 / ( fBinFin * VEPI ); 
  
  LogPrintf(LOG_DEBUG,"StackSlide patchsize is %f rad x %f rad\n", patchSizeX, patchSizeY);
    
  /* set up sky grid */    
  alphaStart = alpha - patchSizeX/2.0;
  alphaEnd   = alpha + patchSizeX/2.0;
  deltaStart = delta + patchSizeY/2.0;
  deltaEnd   = delta + patchSizeY/2.0; 
  dDelta = 1.0/(VTOT *pixelFactor * fBinFin);
  ndelta = floor( (deltaEnd - deltaStart)/dDelta + 0.5);
  if (ndelta < 1) ndelta = 1;  /* do at least one value of delta below */
  /* Will compute dAlpha and nalpha for each delta in the loops below */

  /* calculate time differences from start of observation time */
  TRY( LALDCreateVector( status->statusPtr, &timeDiffV, nStacks), status);

  for (k=0; k<nStacks; k++) {
    REAL8 tMidStack;
    tMidStack = XLALGPSGetREAL8(tsMid->data + k);
    timeDiffV->data[k] = tMidStack - refTime;
  }

  /* if there are residual spindowns */
  /* f1jump = 1.0 / timeDiffV->data[nStacks - 1]; */ /* resolution in residual fdot */  
  /* dfdot     = deltaF * f1jump; */
  dfdot     = params->dfdot; /* 12/14/06 gm; dfdot now a user parameter; default is df1dot/nStacks1; not nfdot set above. */
  fdotStart = fdot - dfdot*(REAL8)(nfdot/2);
  fdotEnd   = fdot + dfdot*(REAL8)(nfdot/2);
  if (nfdot < 1) nfdot = 1; /* do at least one value of fdot below */

  /* The input parameter space point */
  inputPoint.refTime = refTimeGPS;
  inputPoint.orbit = NULL;
  INIT_MEM ( inputPoint.fkdot );
  inputPoint.fkdot[0] = fmid;
  inputPoint.fkdot[1] = fdot;
  inputPoint.Alpha = alpha;
  inputPoint.Delta = delta;

  /* Values for output parameter space point that do not change */  
  outputPoint.refTime = refTimeGPS;
  outputPoint.orbit = NULL;
  INIT_MEM ( outputPoint.fkdot );

  /* uncertainties in the output parameter space point */
  outputPointUnc.refTime = refTimeGPS;
  outputPointUnc.orbit = NULL;
  INIT_MEM ( outputPointUnc.fkdot );
  outputPointUnc.fkdot[0] = deltaF;
  outputPointUnc.fkdot[1] = dfdot;
  outputPointUnc.Delta = dDelta;

  /* Let the loops over sky position and spindown values begin! */

  /* loop over delta */
  for (idelta = 0; idelta<ndelta; idelta++) {

      thisDelta = deltaStart + ((REAL8)idelta)*dDelta;

      if ( (thisDelta < LAL_PI_2) && (thisDelta > -1.0*LAL_PI_2) ) {
         /* Find the spacing in alpha for thisDelta and adjust this to */
         /* fit evenly spaced points between alphaEnd and alphaStart */
         dAlpha = dDelta/cos(thisDelta); 
         nalpha = ceil( (alphaEnd - alphaStart)/dAlpha );
         if (nalpha < 1) nalpha = 1; /* do at least one value of alpha */
         dAlpha = (alphaEnd - alphaStart)/((REAL8)nalpha);
      } else {
         dAlpha = 0.0;
         nalpha = 1;
      }
      outputPointUnc.Alpha = dAlpha;

      /* loop over alpha */
      for (ialpha = 0; ialpha<nalpha; ialpha++) {
          thisAlpha = alphaStart + ((REAL8)ialpha)*dAlpha;

          /* loop over fdot */
          for (ifdot = 0; ifdot<nfdot; ifdot++) {
              thisFdot = fdotStart + ((REAL8)ifdot)*dfdot;

              /* for thisAlpha, thisDelta, and thisFdot, compute the StackSlide Sum (average) of the Fvecs */
              outputPoint.fkdot[0] = fmid; /* This will get changed by solving the Master Equation */
              outputPoint.fkdot[1] = thisFdot;
              outputPoint.Alpha = thisAlpha;
              outputPoint.Delta = thisDelta;

              /* loop over each stack  */
              for (k=0; k<nStacks; k++) {

                  /* COMPUTE f(t) using the master equation and find bin offset for  */
                  /* the frequency in the middle of the band, fmid.                  */
                  /* ASSUMES same offset for entire band, assumed to be very narrow. */
                  TRY ( LALappsFindFreqFromMasterEquation(status->statusPtr,&outputPoint,&inputPoint,(vel->data + 3*k),timeDiffV->data[k],numSpindown), status);

                  offset = floor( (outputPoint.fkdot[0] - fmid)*tEffSTK + 0.5 );
                  
                  #ifdef PRINT_STACKSLIDE_BINOFFSETS
                     LogPrintf(LOG_DETAIL,"offset = %i for stack %i.\n",offset,k);
                  #endif

                  pFVecData = vecF->data[k].data->data;
                  /* loop over frequency bins */
                  if (weightsV == NULL) {
                   /* WITHOUT WEIGHTS */
                   for(j=0; j<nSearchBins; j++) {
                     /* offsetj = j + offset; */
                     offsetj = j + halfExtraBinsFstat + offset;    /* 12/14/06 gm; account for extra bins */
                     if (offsetj < 0) j = 0;                       /* 12/14/06 gm; do not go off the end of the stack */
                     /* if (offsetj > nSearchBinsm1) j = nSearchBinsm1; */
                     if (offsetj > nStackBinsm1) j = nStackBinsm1; /* 12/14/06 gm; do not go off the end of the stack */

                     if (k == 0) {
                        pstackslideData[j] = pFVecData[offsetj];
                     } else {
                        pstackslideData[j] += pFVecData[offsetj];
                     } /* END if (k == 0) */
                   } /* END for(j=0; j<nSearchBins; j++) */
                 } else {
                   /* WITH WEIGHTS */
                   thisWeight = weightsV->data[k];
                   for(j=0; j<nSearchBins; j++) {
                     /* offsetj = j + offset; */
                     offsetj = j + halfExtraBinsFstat + offset;    /* 12/14/06 gm; account for extra bins */
                     if (offsetj < 0) j = 0;                       /* 12/14/06 gm; do not go off the end of the stack */
                     /* if (offsetj > nSearchBinsm1) j = nSearchBinsm1; */
                     if (offsetj > nStackBinsm1) j = nStackBinsm1; /* 12/14/06 gm; do not go off the end of the stack */

                     if (k == 0) {
                        pstackslideData[j] = thisWeight*pFVecData[offsetj];
                     } else {
                        pstackslideData[j] += thisWeight*pFVecData[offsetj];
                     } /* END if (k == 0) */
                   } /* END for(j=0; j<nSearchBins; j++) */
                 } /* END if (weightsV == NULL) */

              } /* END for (k=0; k<nStacks; k++) */

             /* TO DO: get candidates */
             /* 12/18/06 gm; Even if using a toplist, get candidate list based on threshold if thresold > 0.0 */
             if ( params->useToplist && threshold <= 0.0 ) {
               TRY(GetStackSlideCandidates_toplist( status->statusPtr, stackslideToplist, &stackslideSum, &outputPoint, &outputPointUnc), status);
             } else {
               TRY(GetStackSlideCandidates_threshold( status->statusPtr, out, &stackslideSum, &outputPoint, &outputPointUnc, threshold), status);
             }

          } /* END for (ifdot = 0; ifdot<nfdot; ifdot++) */

      } /* END for (ialpha = 0; ialpha<nalpha; ialpha++) */

  } /* END for (idelta = 0; idelta<ndelta; idelta++) */
  
  /* free remaining memory */
  TRY( LALDDestroyVector( status->statusPtr, &timeDiffV), status);
  
  LALFree(stackslideSum.data->data);
  LALFree(stackslideSum.data);

  /* copy toplist candidates to output structure if necessary */
  if ( params->useToplist && threshold <= 0.0 ) {
    for ( uk=0; uk<stackslideToplist->elems; uk++) {
      out->list[uk] = *((SemiCohCandidate *)(toplist_elem(stackslideToplist, uk)));
    }
    out->nCandidates = stackslideToplist->elems;
    free_toplist(&stackslideToplist);
  } else if ( params->useToplist ) {
    /* 12/18/06 gm; Even if using a toplist, get candidate list based on threshold if thresold > 0.0 */
    /* First put candidates into the top list: */
    for ( uk=0; uk<out->nCandidates; uk++) {
        insert_into_toplist(stackslideToplist, out->list + uk);
    }
    /* Now put back into the output list: */
    for ( uk=0; uk<stackslideToplist->elems; uk++) {
      out->list[uk] = *((SemiCohCandidate *)(toplist_elem(stackslideToplist, uk)));
    }
    out->nCandidates = stackslideToplist->elems;
    free_toplist(&stackslideToplist);
  }

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* END StackSlideVecF */




/** \brief Function StackSlides a vector of Fstat frequency series or any REAL8FrequencySeriesVector.
           This is similar to StackSlideVecF but adapted to calculate the hough number count and to be as 
	   similar to Hough as possible but without using the hough look-up-tables.
    \param out SemiCohCandidateList is a list of candidates
    \param vecF is a vector of Fstat frequency series or any REAL4FrequencySeriesVector.
    \param params is a pointer to SemiCoherentParams
    \out SemiCohCandidateList is a list of candidates
*/
void StackSlideVecF_HoughMode(LALStatus *status,
			      SemiCohCandidateList  *out,        /**< output candidates */
			      REAL4FrequencySeriesVector *vecF,  /**< vector with Fstat values or any REAL8FrequencySeriesVector */
			      SemiCoherentParams *params)        /**< input parameters  */
{

  UINT2  xSide, ySide, maxNBins, maxNBorders;
  INT8  fBinIni, fBinFin, fBin;
  INT4  iHmap, nfdot;
  UINT4 k, nStacks ;
  REAL8 deltaF, dfdot, alpha, delta;
  REAL8 patchSizeX, patchSizeY;
  REAL8VectorSequence *vel, *pos;
  REAL8 fdot, refTime;
  LIGOTimeGPS refTimeGPS;
  LIGOTimeGPSVector   *tsMid;
  REAL8Vector *timeDiffV=NULL;
  UINT8Vector hist; /* histogram vector */ 
  UINT8Vector histTotal; /* total histogram vector */
  HoughStats stats; /* statistics struct */
  CHAR *fileStats = NULL;
  FILE *fpStats = NULL;

  /* a minimal number of hough structs needed */
  HOUGHMapTotal ht;
  HOUGHDemodPar   parDem;  /* demodulation parameters */
  HOUGHSizePar    parSize; 
  HOUGHResolutionPar parRes;   /* patch grid information */
  HOUGHPatchGrid  patch;   /* Patch description */ 
  INT4 nfSize;

  toplist_t *houghToplist;
  UINT8FrequencyIndexVector freqInd; /* for trajectory in time-freq plane */

  INITSTATUS( status, "StackSlideVecF_HoughMode", rcsid );
  ATTATCHSTATUSPTR (status);


  /* check input is not null */
  if ( out == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  
  if ( out->length == 0 ) {
    ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  }  
  if ( out->list == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  }  
  if ( vecF == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  
  if ( vecF->length == 0 ) {
    ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  }  
  if ( vecF->data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  
  if ( params == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  



  /* copy some parameters from peakgram vector */
  deltaF = vecF->data->deltaF;
  nStacks = vecF->length;
  fBinIni =  (UINT4)(vecF->data[0].f0/deltaF + 0.5);
  fBinFin = fBinIni +  vecF->data->data->length - 1;

  
  /* copy some params to local variables */
  nfdot = params->nfdot;
  dfdot = params->dfdot;
  alpha = params->alpha;
  delta = params->delta;
  vel = params->vel;
  pos = params->pos;
  fdot = params->fdot;
  tsMid = params->tsMid;
  refTimeGPS = params->refTime;
  refTime = XLALGPSGetREAL8(&refTimeGPS);

  /* set patch size */
  /* this is supposed to be the "educated guess"
     delta theta = 1.0 / (Tcoh * f0 * Vepi )
     where Tcoh is coherent time baseline,
     f0 is frequency and Vepi is rotational velocity
     of detector */
  patchSizeX = params->patchSizeX;
  patchSizeY = params->patchSizeY;

  /* calculate time differences from start of observation time for each stack*/
  TRY( LALDCreateVector( status->statusPtr, &timeDiffV, nStacks), status);
  
  for (k=0; k<nStacks; k++) {
    REAL8 tMidStack;
    tMidStack = XLALGPSGetREAL8(tsMid->data + k);
    timeDiffV->data[k] = tMidStack - refTime;
  }

  /* residual spindown trajectory */
  freqInd.deltaF = deltaF;
  freqInd.length = nStacks;
  freqInd.data = NULL;
  freqInd.data =  ( UINT8 *)LALCalloc(1,nStacks*sizeof(UINT8));
  if ( freqInd.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  


  /* resolution in space of residual spindowns */
  ht.dFdot.length = 1;
  ht.dFdot.data = NULL;
  ht.dFdot.data = (REAL8 *)LALCalloc( 1, ht.dFdot.length * sizeof(REAL8));
  if ( ht.dFdot.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  

  /* the residual spindowns */
  ht.spinRes.length = 1;
  ht.spinRes.data = NULL;
  ht.spinRes.data = (REAL8 *)LALCalloc( 1, ht.spinRes.length*sizeof(REAL8));
  if ( ht.spinRes.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  

  /* the residual spindowns */
  ht.spinDem.length = 1;
  ht.spinDem.data = NULL;
  ht.spinDem.data = (REAL8 *)LALCalloc( 1, ht.spinRes.length*sizeof(REAL8));
  if ( ht.spinDem.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  

  /* the demodulation params */
  parDem.deltaF = deltaF;
  parDem.skyPatch.alpha = alpha;
  parDem.skyPatch.delta = delta;
  parDem.spin.length = 1;
  parDem.spin.data = NULL;
  parDem.spin.data = (REAL8 *)LALCalloc(1, sizeof(REAL8));
  if ( parDem.spin.data == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  
  parDem.spin.data[0] = fdot;

  /* the skygrid resolution params */
  parRes.deltaF = deltaF;
  parRes.patchSkySizeX  = patchSizeX;
  parRes.patchSkySizeY  = patchSizeY;
  parRes.pixelFactor = params->pixelFactor;
  parRes.pixErr = PIXERR;
  parRes.linErr = LINERR;
  parRes.vTotC = VTOT;


  {
    REAL8 maxTimeDiff, startTimeDiff, endTimeDiff;

    startTimeDiff = fabs(timeDiffV->data[0]);
    endTimeDiff = fabs(timeDiffV->data[timeDiffV->length - 1]);
    maxTimeDiff = SSMAX( startTimeDiff, endTimeDiff);

    /* set number of freq. bins for which LUTs will be calculated */
    /* this sets the range of residual spindowns values */
    /* phmdVS.nfSize  = 2*nfdotBy2 + 1; */
    nfSize  = 2 * floor((nfdot-1) * (REAL4)(dfdot * maxTimeDiff / deltaF) + 0.5f) + 1; 
  }

  /* adjust fBinIni and fBinFin to take maxNBins into account */
  /* and make sure that we have fstat values for sufficient number of bins */
  parRes.f0Bin =  fBinIni;      

  fBinIni += params->extraBinsFstat;
  fBinFin -= params->extraBinsFstat;
  /* this is not very clean -- the Fstat calculation has to know how many extra bins are needed */

  LogPrintf(LOG_DETAIL, "Freq. range analyzed by Hough = [%fHz - %fHz] (%d bins)\n", 
	    fBinIni*deltaF, fBinFin*deltaF, fBinFin - fBinIni + 1);
  ASSERT ( fBinIni < fBinFin, status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );

  /* initialise number of candidates -- this means that any previous candidates 
     stored in the list will be lost for all practical purposes*/
  out->nCandidates = 0; 
  
  /* create toplist of candidates */
  if (params->useToplist) {
    create_toplist(&houghToplist, out->length, sizeof(SemiCohCandidate), smallerStackSlide);
  }
  else { 
    /* if no toplist then use number of hough maps */
    INT4 numHmaps = (fBinFin - fBinIni + 1) * nfSize;
    if (out->length != numHmaps) {
      out->length = numHmaps;
      out->list = (SemiCohCandidate *)LALRealloc( out->list, out->length * sizeof(SemiCohCandidate));
      if ( out->list == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }  
    }
  }

  /*------------------ start main Hough calculation ---------------------*/

  /* initialization */  
  fBin= fBinIni; /* initial search bin */
  iHmap = 0; /* hough map index */

  while( fBin <= fBinFin ){
    INT8 fBinSearch, fBinSearchMax;
    UINT4 i,j; 
    REAL8UnitPolarCoor sourceLocation;
    	
    parRes.f0Bin =  fBin;      
    TRY( LALHOUGHComputeSizePar( status->statusPtr, &parSize, &parRes ),  status );
    xSide = parSize.xSide;
    ySide = parSize.ySide;

    maxNBins = parSize.maxNBins;
    maxNBorders = parSize.maxNBorders;
	
    /*------------------ create patch grid at fBin ----------------------*/
    patch.xSide = xSide;
    patch.ySide = ySide;
    patch.xCoor = NULL;
    patch.yCoor = NULL;
    patch.xCoor = (REAL8 *)LALCalloc(1,xSide*sizeof(REAL8));
    if ( patch.xCoor == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }  

    patch.yCoor = (REAL8 *)LALCalloc(1,ySide*sizeof(REAL8));
    if ( patch.yCoor == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }  
    TRY( LALHOUGHFillPatchGrid( status->statusPtr, &patch, &parSize ), status );
  
    /*-------------- initializing the Total Hough map space ------------*/   
    ht.xSide = xSide;
    ht.ySide = ySide;
    ht.skyPatch.alpha = alpha;
    ht.skyPatch.delta = delta;
    ht.mObsCoh = nStacks;
    ht.deltaF = deltaF;
    ht.spinDem.data[0] = fdot;
    ht.patchSizeX = patchSizeX;
    ht.patchSizeY = patchSizeY;
    ht.dFdot.data[0] = dfdot;
    ht.map   = NULL;
    ht.map   = (HoughTT *)LALCalloc(1,xSide*ySide*sizeof(HoughTT));
    if ( ht.map == NULL ) {
      ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
    }  
    
    TRY( LALHOUGHInitializeHT( status->statusPtr, &ht, &patch), status); /*not needed */
 
    /*  Search frequency interval possible using the same LUTs */
    fBinSearch = fBin;
    fBinSearchMax = fBin + parSize.nFreqValid - 1;   
    
    /* loop over frequencies */    
    while ( (fBinSearch <= fBinFin) && (fBinSearch < fBinSearchMax) )  { 

      /* finally we can construct the hough maps and select candidates */
      {
	INT4   n, nfdotBy2;

	nfdotBy2 = nfdot/2;
	ht.f0Bin = fBinSearch;

	/*loop over all values of residual spindown */
	/* check limits of loop */
	for( n = -nfdotBy2; n <= nfdotBy2 ; n++ ){ 

	  ht.spinRes.data[0] =  n*dfdot; 
	  
	  for (j=0; j < (UINT4)nStacks; j++) {
	    freqInd.data[j] = fBinSearch + floor( (REAL4)(timeDiffV->data[j]*n*dfdot/deltaF) + 0.5f);
	  }






	  /* get candidates */
	  if ( params->useToplist ) 
	    {
	      TRY(GetHoughCandidates_toplist( status->statusPtr, houghToplist, &ht, &patch, &parDem), status);
	    }
	  else 
	    {
	      TRY(GetHoughCandidates_threshold( status->statusPtr, out, &ht, &patch, &parDem, params->threshold), status);
	    }
	  
	  /* increment hough map index */ 	  
	  ++iHmap;
	  
	} /* end loop over spindown trajectories */

      } /* end of block for calculating total hough maps */
      

      /*------ shift the search freq. & PHMD structure 1 freq.bin -------*/

      ++fBinSearch;


    }   /* closing while loop over fBinSearch */


  }  


  /* copy toplist candidates to output structure if necessary */
  if ( params->useToplist ) {
    for ( k=0; k<houghToplist->elems; k++) {
      out->list[k] = *((SemiCohCandidate *)(toplist_elem(houghToplist, k)));
    }
    out->nCandidates = houghToplist->elems;
    free_toplist(&houghToplist);
  }

  LALFree(ht.map);

  /* free remaining memory */
  LALFree(ht.spinRes.data);
  LALFree(ht.spinDem.data);
  LALFree(ht.dFdot.data);

  TRY( LALDDestroyVector( status->statusPtr, &timeDiffV), status);
  LALFree(freqInd.data);

  LALFree(patch.xCoor);
  LALFree(patch.yCoor);

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* END StackSlideVecF */









/* Calculate f(t) using the master equation given by Eq. 6.18 in gr-qc/0407001 */
/* Returns f(t) in outputPoint.fkdot[0] */
void LALappsFindFreqFromMasterEquation(LALStatus *status, 
                                       PulsarDopplerParams *outputPoint,  /* outputs f(t) for output sky position and spindown values                       */
                                       PulsarDopplerParams *inputPoint,   /* input demodulation f0, sky position, and spindown values                       */
                                       REAL8 *vel,                        /* vx = vel[0], vy = vel[1], vz = vel[2] = ave detector velocity                  */
                                       REAL8 deltaT,                      /* time since the reference time                                                  */
                                       UINT2 numSpindown)                 /* Number of spindown values == high deriv. of include == 1 if just df/dt, etc... */
{
                  UINT2 k;
                  REAL8 f0, F0, F0zeta, alpha, delta, cosAlpha, cosDelta, sinAlpha, sinDelta;
                  REAL8 nx, ny, nz, ndx, ndy, ndz;
                  REAL8 vx, vy, vz;
                  REAL8 kFact, deltaTPowk;
                  PulsarSpins inputfkdot; /* input demodulation spindown values */
                  PulsarSpins deltafkdot; /* residual spindown values */

                  INITSTATUS( status, "LALappsFindFreqFromMasterEquation", rcsid );
                  ATTATCHSTATUSPTR (status);
  
                  ASSERT ( outputPoint != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );
                  ASSERT ( inputPoint != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );
                  ASSERT ( vel != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );

                  /* the x, y, and z components of the input demodulation sky position: */
                  alpha = inputPoint->Alpha;
                  delta = inputPoint->Delta;
                  cosAlpha = cos(alpha);
                  cosDelta = cos(delta);
                  sinAlpha = sin(alpha);
                  sinDelta = sin(delta);
                  ndx = cosDelta*cosAlpha;
                  ndy = cosDelta*sinDelta;
                  ndz = sinDelta;

                  /* the x, y, and z components of the output sky position: */
                  alpha = outputPoint->Alpha;
                  delta = outputPoint->Delta;
                  cosAlpha = cos(alpha);
                  cosDelta = cos(delta);
                  sinAlpha = sin(alpha);
                  sinDelta = sin(delta);
                  nx = cosDelta*cosAlpha;
                  ny = cosDelta*sinDelta;
                  nz = sinDelta;
                  
                  f0 = inputPoint->fkdot[0];  /* input f0 */

                  /* input and residual spindown values: */
                  deltafkdot[0] = 0;  /* unused */
                  inputfkdot[0] = f0; /* unused */
                  for (k=1; k<=numSpindown; k++) {
                      inputfkdot[k] = inputPoint->fkdot[k];
                      deltafkdot[k] = outputPoint->fkdot[k] - inputPoint->fkdot[k];
                  }

                  /* the x, y, and z components of velocity of the detector(s) for this stack */
                  vx = vel[0];
                  vy = vel[1];
                  vz = vel[2];

                  /* Compute F0 */
                  F0 = f0;
                  kFact = 1.0;
                  deltaTPowk = 1.0;
                  for (k=1; k<=numSpindown; k++) {
                      kFact *= ((REAL8)k);
                      deltaTPowk *= deltaT;
                      F0 += deltafkdot[k]*deltaTPowk/kFact;
                  }
                  
                  /* Compute F0 plus spindown; call this F0zeta.  See the master equation. */                  
                  F0zeta = F0;
                  kFact = 1.0;
                  deltaTPowk = 1.0;
                  for (k=1; k<=numSpindown; k++) {
                      kFact *= ((REAL8)k);
                      deltaTPowk *= deltaT;
                      F0zeta += inputfkdot[k]*deltaTPowk/kFact;
                  }

                  /* Compute the output frequency. */
                  /* NOTE that the small correction fkdot[k]*(deltaT^(k-1)/(k-1)!)*(r - r0)/c is ignored */
                  outputPoint->fkdot[0] = F0 + F0zeta*( vx*(nx-ndx) + vy*(ny-ndy) + vz*(nz-ndz) );

                  DETATCHSTATUSPTR (status);
                  RETURN(status);
} /* END LALappsFindFreqFromMasterEquation */

/* Get StackSlide candidates using a fixed threshold */
void GetStackSlideCandidates_threshold(LALStatus *status,
                                       SemiCohCandidateList *out,            /* output list of candidates */
                                       REAL8FrequencySeries *stackslideSum,  /* input stackslide sum of F stat values */
                                       PulsarDopplerParams *outputPoint,     /* parameter space point for which to output candidate */
                                       PulsarDopplerParams *outputPointUnc,  /* uncertainties in parameter space point for which to output candidate */
                                       REAL8 threshold)                      /* threshold on significance */
{
  REAL8 deltaF, f0, freq;
  INT4 j, jminus1, jplus1, nSearchBins, nSearchBinsm1, numCandidates;
  SemiCohCandidate thisCandidate;
  BOOLEAN isLocalMax = TRUE;
  REAL8 thisSig, thisSigMinus1, thisSigPlus1;
  REAL8 *pstackslideData;  /* temporary pointer */
  
  INITSTATUS( status, "GetStackSlideCandidates_threshold", rcsid );
  ATTATCHSTATUSPTR (status);

  ASSERT ( out != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );
  ASSERT ( out->length > 0, status, STACKSLIDEFSTAT_EVAL, STACKSLIDEFSTAT_MSGEVAL );
  ASSERT ( out->list != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );
  ASSERT ( stackslideSum != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );
  ASSERT ( stackslideSum->data->data != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );
  ASSERT ( outputPoint != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );

  pstackslideData = stackslideSum->data->data;

  f0 = stackslideSum->f0;
  deltaF = stackslideSum->deltaF;
  thisCandidate.dFreq = deltaF; 
  thisCandidate.fdot = outputPoint->fkdot[1];
  thisCandidate.dFdot = outputPointUnc->fkdot[1];
  thisCandidate.alpha = outputPoint->Alpha;
  thisCandidate.delta = outputPoint->Delta;
  thisCandidate.dAlpha = outputPointUnc->Alpha;
  thisCandidate.dDelta = outputPointUnc->Delta;

  numCandidates = out->nCandidates;  
  nSearchBins = stackslideSum->data->length;
  nSearchBinsm1 = nSearchBins - 1;

  /* Search frequencies for candidates above threshold that are local maxima */
  for(j=0; j<nSearchBins; j++) {
     freq = f0 + ((REAL8)j)*deltaF;
          
     /* realloc list if necessary */
     if (numCandidates >= out->length) {
         out->length += BLOCKSIZE_REALLOC;
         out->list = (SemiCohCandidate *)LALRealloc( out->list, out->length * sizeof(SemiCohCandidate));
         LogPrintf(LOG_DEBUG, "Need to realloc StackSlide candidate list to %d entries\n", out->length);
     } /* need a safeguard to ensure that the reallocs don't happen too often */

     thisSig = pstackslideData[j]; /* Should we do more than this to find significance? */

     /* candidates are above threshold and local maxima */
     if (thisSig > threshold) {
        /* 12/14/06 gm; improve local maxima checking */
        jminus1 = j - 1;
        jplus1 = j + 1;
        if (jminus1 < 0) {
           if (jplus1 > nSearchBinsm1) jplus1 = nSearchBinsm1; /* just to be safe */
           thisSigPlus1 = pstackslideData[jplus1];
           isLocalMax = (thisSig > thisSigPlus1);
           #ifdef INCLUDE_EXTRA_STACKSLIDE_THREHOLDCHECK
              isLocalMax = ( isLocalMax && ( thisSig > 2.0*THRESHOLDFRACSS*thisSigPlus1 ) );
           #endif
        } else if (jplus1 > nSearchBinsm1) {
           if (jminus1 < 0) jminus1 = 0; /* just to be safe */
           thisSigMinus1 = pstackslideData[jminus1];
           isLocalMax = (thisSig > thisSigMinus1);
           #ifdef INCLUDE_EXTRA_STACKSLIDE_THREHOLDCHECK
              isLocalMax = ( isLocalMax && ( thisSig > 2.0*THRESHOLDFRACSS*thisSigMinus1 ) );
           #endif
        } else {
           thisSigMinus1 = pstackslideData[jminus1];
           thisSigPlus1 = pstackslideData[jplus1];
           isLocalMax = (thisSig > thisSigMinus1) && (thisSig > thisSigPlus1);
           #ifdef INCLUDE_EXTRA_STACKSLIDE_THREHOLDCHECK
              isLocalMax = ( isLocalMax && ( thisSig > THRESHOLDFRACSS*(thisSigMinus1 + thisSigPlus1) ) );
           #endif
        }
        if ( (numCandidates < out->length) && isLocalMax ) {
           thisCandidate.significance = thisSig;
           thisCandidate.freq =  freq;
           out->list[numCandidates] = thisCandidate;
           numCandidates++;
           out->nCandidates = numCandidates;
        }
     }
  } /* END for(j=0; j<nSearchBins; j++) */

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* END GetStackSlideCandidates_threshold */

/* Get StackSlide candidates as a toplist */
void GetStackSlideCandidates_toplist(LALStatus *status,
                                     toplist_t *list,
                                     REAL8FrequencySeries *stackslideSum,  /* input stackslide sum of F stat values */
                                     PulsarDopplerParams *outputPoint,     /* parameter space point for which to output candidate */
                                     PulsarDopplerParams *outputPointUnc)  /* uncertainties in parameter space point for which to output candidate */
{
  REAL8 deltaF, f0, freq;
  INT4 j, nSearchBins;
  SemiCohCandidate thisCandidate;
  REAL8 *pstackslideData;  /* temporary pointer */

  INITSTATUS( status, "GetStackSlideCandidates_toplist", rcsid );
  ATTATCHSTATUSPTR (status);

  ASSERT ( stackslideSum != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );
  ASSERT ( stackslideSum->data->data != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );
  ASSERT ( outputPoint != NULL, status, STACKSLIDEFSTAT_ENULL, STACKSLIDEFSTAT_MSGENULL );

  pstackslideData = stackslideSum->data->data;

  f0 = stackslideSum->f0;
  deltaF = stackslideSum->deltaF;
  thisCandidate.dFreq = deltaF; 
  thisCandidate.fdot = outputPoint->fkdot[1];
  thisCandidate.dFdot = outputPointUnc->fkdot[1];
  thisCandidate.alpha = outputPoint->Alpha;
  thisCandidate.delta = outputPoint->Delta;
  thisCandidate.dAlpha = outputPointUnc->Alpha;
  thisCandidate.dDelta = outputPointUnc->Delta;

  nSearchBins = stackslideSum->data->length;

  /* Search frequencies for candidates above threshold that are local maxima */
  for(j=0; j<nSearchBins; j++) {

     freq = f0 + ((REAL8)j)*deltaF;
     thisCandidate.freq =  freq;
     thisCandidate.significance = pstackslideData[j]; /* Should we do more than this to find significance? */
     
     insert_into_toplist(list, &thisCandidate);

  } /* END for(j=0; j<nSearchBins; j++) */

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* END GetStackSlideCandidates_toplist */
