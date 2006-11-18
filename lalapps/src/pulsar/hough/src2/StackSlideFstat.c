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

#include "./StackSlideFstat.h"

RCSID( "$Id$");

#define TRUE (1==1)
#define FALSE (1==0)

extern int lalDebugLevel;

/** \brief Function StackSlides a vector of Fstat frequency series or any REAL8FrequencySeriesVector.
    \param out SemiCohCandidateList is a list of candidates
    \param vecF is a vector of Fstat frequency series or any REAL8FrequencySeriesVector.
    \param params is a pointer to SemiCoherentParams
    \out SemiCohCandidateList is a list of candidates
*/
void LALappsStackSlideVecF(LALStatus *status,
			  SemiCohCandidateList  *out,        /* output candidates */
			  REAL8FrequencySeriesVector *vecF,  /* vector with Fstat values or any REAL8FrequencySeriesVector */
			  SemiCoherentParams *params)        /* input parameters  */
{

  /* parameters for the semicoherent stage -- hough or stackslide */
  typedef struct tagSemiCoherentParams {
    LIGOTimeGPSVector *tsMid;  /**< timestamps of mid points of stacks */
    LIGOTimeGPS refTime;       /**< reference time for f, fdot definition */
    REAL8VectorSequence *vel;  /**< detector velocity for each stack */
    REAL8VectorSequence *pos;  /**< detector position for each stack */
    REAL8 alpha;               /**< right ascension of demodulation point */
    REAL8 delta;               /**< declination of demodulation point*/
    REAL8 pixelFactor;         /**< Resolution of semicoherent sky-grid */
    REAL8 patchSizeX;          /**< Size of semicoherent sky-patch */
    REAL8 patchSizeY;          /**< Size of semicoherent sky-patch */
    REAL8 fdot;                /**< spindown value of demodulation point */
    UINT4 nfdot;               /**< number of fdot values to search over */ 
    CHAR *outBaseName;         /**< file for writing output -- if chosen */
    BOOLEAN useToplist;        /**< Use a toplist for producing candidates? */
    REAL8  threshold;          /**< Threshold for candidate selection */
  } SemiCoherentParams;

  /** one candidate */
  typedef struct tagSemiCohCandidate {
    REAL8 freq;        /**< frequency */
    REAL8 alpha;       /**< right ascension */
    REAL8 delta;       /**< declination */
    REAL8 fdot;        /**< spindown */
    REAL8 dFreq;       /**< frequency error */
    REAL8 dAlpha;      /**< alpha error */
    REAL8 dDelta ;     /**< delta error */
    REAL8 dFdot;       /**< fdot error */
    REAL8 significance;/**< significance */
  } SemiCohCandidate;  

  /** structure for storing candidates */
  typedef struct tagSemiCohCandidateList {
    LIGOTimeGPS refTime;       /**< reference time for candidates */
    INT4 length;               /**< maximum allowed length of vectors */
    INT4 nCandidates;          /**< number of candidates -- must be less than length */
    SemiCohCandidate *list;    /**> list of candidates */
  } SemiCohCandidateList;

  REAL8FrequencySeries stackslideSum;  /* The output of StackSliding the vecF values */
  
  REAL8 *pstackslideData;  /* temporary pointer */
  REAL8 *pFVecData;  /* temporary pointer */
  
  REAL8 pixelFactor;
  REAL8 alphaStart, alphaEnd, dAlpha;
  REAL8 deltaStart, deltaEnd, dDelta;
  REAL8 fdotStart, fdotEnd, dfdot;
  UINT4 ialpha,nalpha,idelta,ndelta,ifdot,nfdot;

  INT4 fBinIni, fBinFin, fBin, nSearchBins;
  INT4 j, k, nStacks, offset;
  REAL8 f0,deltaF, alpha, delta;
  REAL8 patchSizeX, patchSizeY, f1jump;
  REAL8VectorSequence *vel, *pos;
  REAL8 fdot, refTime;
  LIGOTimeGPS refTimeGPS;
  LIGOTimeGPSVector   *tsMid;
  REAL8Vector *timeDiffV=NULL;

  /* toplist_t *stackslideToplist; */ /* TO BE ADDED IN */

  INITSTATUS( status, "LALappsStackSlideVecF", rcsid );
  ATTATCHSTATUSPTR (status);

  /* create toplist of candidates */
  /* if (params->useToplist) {
    create_toplist(&stackslideToplist, out->length, sizeof(SemiCohCandidate), smallerHough);
  } */

  /* copy some params to local variables */
  pixelFactor = params->pixelFactor;   
  nStacks = vecF->length;  
  nSearchBins = vecF->data->data->length;  
  f0 = vecF->data[0].f0;
  deltaF = vecF->data[0].deltaF;  
  fBinIni = floor( f0/deltaF + 0.5);
  fBinFin = fBinIni + nSearchBins - 1;
  nfdot = params->nfdot;
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
  /* if patchsize is not negative, then set default value */
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
  if (nfdot < 1); nfdot = 1; /* do at least one value of fdto below */
  
  offset = 0;
  
  /* Let the loops over sky position and spindown values begin! */
  
  /* loop over delta */
  for (idelta = 0; idelta<ndelta; idelta++) {
      
      thisDelta = deltaStart + ((REAL8)idelta)*dDelta;
      
      if ( (thisDelta < LAL_PI_2) && (thisDelta > -1.0*LAL_PI_2) ) {
         dAlpha = dDelta/cos(thisDelta);
         nalpha = floor( (alphaEnd - alphaStart)/dAlpha + 0.5);
         if (nalpha < 1) nalpha = 1; /* do at least one value of alpha */
      } else {
         dAlpha = 0.0;
         nalpha = 1;
      }
      
      /* loop over alpha */
      for (ialpha = 0; ialpha<nalpha; ialpha++) {
          thisAlpha = alphaStart + ((REAL8)ialpha)*dAlpha;

          /* loop over fdot */
          for (ifdot = 0; ifdot<nfdot; ifdot++) {
              thisFdot = fdotStart + ((REAL8)ifdot)*dfdot;

              /* for thisAlpha, thisDelta, and thisFdot, compute the StackSlide Sum (average) of the Fvecs */
              
              /* loop over each stack  */
              for (k=0; k<nStacks; k++) {
                  pV = FstatVect->data[k].data->data;
                  
                  /* TO DO: COMPUTE f(t) using master equation and find bin offset */
                  /* ASSUMES same offset for entire band, assumed to be very narrow */

                  /* loop over frequency bins */
                  for(j=0; j<nSearchBins; j++) {
                     if (k == 0) {
                        pstackslideData[j] += pV[j + offset];
                     } else {
                        pstackslideData[j] += pV[j + offset];
                     } /* END if (k == 0) */
                  } /* END for(j=0; j<nSearchBins; j++) */

              } /* END for (k=0; k<nStacks; k++) */
              
              /* TO DO: find candidates */

          } /* END for (ifdot = 0; ifdot<nfdot; ifdot++) */

      } /* END for (ialpha = 0; ialpha<nalpha; ialpha++) */

  } /* END for (idelta = 0; idelta<ndelta; idelta++) */
  
  /* free remaining memory */
  TRY( LALDDestroyVector( status->statusPtr, &timeDiffV), status);
  
  LALFree(stackslideSum.data->data);
  LALFree(stackslideSum.data);

  /* copy toplist candidates to output structure if necessary */
  /* if ( params->useToplist ) {
    for ( k=0; k<stackslideToplist->elems; k++) {
      out->list[k] = *((SemiCohCandidate *)(toplist_elem(stackslideToplist, k)));
    }
    out->nCandidates = stackslideToplist->elems;
    free_toplist(&stackslideToplist);
  } */

  DETATCHSTATUSPTR (status);
  RETURN(status);

}
