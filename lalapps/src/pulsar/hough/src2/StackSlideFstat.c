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

/* include files: */
#include "./StackSlideFstat.h"

RCSID( "$Id$");

#define TRUE (1==1)
#define FALSE (1==0)

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
                    REAL8FrequencySeriesVector *vecF,  /* vector with Fstat values or any REAL8FrequencySeriesVector */
                    SemiCoherentParams *params)        /* input parameters  */
{
  REAL8FrequencySeries stackslideSum;  /* The output of StackSliding the vecF values */
  
  REAL8 *pstackslideData;  /* temporary pointer */
  REAL8 *pFVecData;        /* temporary pointer */

  REAL8 pixelFactor;
  REAL8 alphaStart, alphaEnd, dAlpha, thisAlpha;
  REAL8 deltaStart, deltaEnd, dDelta, thisDelta;
  REAL8 fdotStart, fdotEnd, dfdot, thisFdot;
  UINT4 ialpha,nalpha,idelta,ndelta,ifdot,nfdot;
  UINT2 numSpindown;

  INT4 fBinIni, fBinFin, nSearchBins, nSearchBinsm1;
  INT4 j, k, nStacks, offset, offsetj;
  UINT4 uk;
  REAL8 f0, deltaF, tEffSTK, fmid, alpha, delta;
  REAL8 patchSizeX, patchSizeY, f1jump;
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
  nSearchBins = vecF->data->data->length;
  nSearchBinsm1 = nSearchBins - 1;
  f0 = vecF->data[0].f0;
  deltaF = vecF->data[0].deltaF;
  tEffSTK = 1.0/deltaF; /*  Effective time baseline a stack */
  fBinIni = floor( f0*tEffSTK + 0.5);
  fBinFin = fBinIni + nSearchBins - 1;
  fmid = vecF->data[0].f0 + ((REAL8)(nSearchBins/2))*deltaF;
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
  TRY ( LALGPStoFloat( status->statusPtr, &refTime, &refTimeGPS), status);

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
    TRY ( LALGPStoFloat ( status->statusPtr, &tMidStack, tsMid->data + k), status);
    timeDiffV->data[k] = tMidStack - refTime;
  }

  /* if there are residual spindowns */
  f1jump = 1.0 / timeDiffV->data[nStacks - 1]; /* resolution in residual fdot */  
  dfdot     = deltaF * f1jump;
  fdotStart = fdot - dfdot*(REAL8)(nfdot/2);
  fdotEnd   = fdot + dfdot*(REAL8)(nfdot/2);
  if (nfdot < 1) nfdot = 1; /* do at least one value of fdot below */

  /* The input parameter space point */
  inputPoint.refTime = refTimeGPS;
  inputPoint.orbit = NULL;
  inputPoint.fkdot[0] = fmid;
  inputPoint.fkdot[1] = fdot;
  inputPoint.fkdot[2] = 0.0;
  inputPoint.fkdot[3] = 0.0;
  inputPoint.Alpha = alpha;
  inputPoint.Delta = delta;

  /* Values for output parameter space point that do not change */  
  outputPoint.refTime = refTimeGPS;
  outputPoint.orbit = NULL;
  outputPoint.fkdot[2] = 0.0;
  outputPoint.fkdot[3] = 0.0;

  /* uncertainties in the output parameter space point */
  outputPointUnc.refTime = refTimeGPS;
  outputPointUnc.orbit = NULL;
  outputPointUnc.fkdot[0] = deltaF;
  outputPointUnc.fkdot[1] = dfdot;
  outputPointUnc.fkdot[2] = 0.0;
  outputPointUnc.fkdot[3] = 0.0;
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
         dAlpha = (alphaEnd - alphaStart)/((REAL8)nalpha);
         if (nalpha < 1) nalpha = 1; /* do at least one value of alpha */
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
                     LogPrintf(LOG_DETAIL,"offset = %i for stack %i.",offset,k);
                  #endif

                  pFVecData = vecF->data[k].data->data;
                  /* loop over frequency bins */
                  if (weightsV == NULL) {
                   /* WITHOUT WEIGHTS */
                   for(j=0; j<nSearchBins; j++) {
                     offsetj = j + offset;
                     if (offsetj < 0) j = 0;  /* TO DO: NEED EXTRA BINS IN STACKS TO ALLOW FOR SLIDING */
                     if (offsetj > nSearchBinsm1) j = nSearchBinsm1;

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
                     offsetj = j + offset;
                     if (offsetj < 0) j = 0;  /* TO DO: NEED EXTRA BINS IN STACKS TO ALLOW FOR SLIDING */
                     if (offsetj > nSearchBinsm1) j = nSearchBinsm1;

                     if (k == 0) {
                        pstackslideData[j] = thisWeight*pFVecData[offsetj];
                     } else {
                        pstackslideData[j] += thisWeight*pFVecData[offsetj];
                     } /* END if (k == 0) */
                   } /* END for(j=0; j<nSearchBins; j++) */
                 } /* END if (weightsV == NULL) */

              } /* END for (k=0; k<nStacks; k++) */

             /* TO DO: get candidates */
             if ( params->useToplist ) {
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
  if ( params->useToplist ) {
    for ( uk=0; uk<stackslideToplist->elems; uk++) {
      out->list[uk] = *((SemiCohCandidate *)(toplist_elem(stackslideToplist, uk)));
    }
    out->nCandidates = stackslideToplist->elems;
    free_toplist(&stackslideToplist);
  }

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
  REAL8 thisSig;
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
        jminus1 = j - 1; if (j < 0) j = 0;
        jplus1 = j + 1;  if (j > nSearchBinsm1) j = nSearchBinsm1;
        isLocalMax = (thisSig > pstackslideData[jminus1]) && (thisSig > pstackslideData[jplus1]);       
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
