/*  
 *  Copyright (C) 2005-2008 Badri Krishnan, Alicia Sintes, Bernd Machenschalk
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
 * 
 */

/*
   This is a local copy of ComputeFstatHoughMap() of HierarchicalSearch.c
   See that file for details

   $Id$
*/

#include "HierarchicalSearch.h"
#include "LocalOptimizationFlags.h"

#if __x86_64__
#include "EinsteinAtHome/hough_x64.ci"
#warning including "EinsteinAtHome/hough_x64.ci"
#elif __SSE2__
#include "EinsteinAtHome/hough_sse2.ci"
#warning including "EinsteinAtHome/hough_sse2.ci"
#elif __i386__
#include "EinsteinAtHome/hough_x87.ci"
#warning including "EinsteinAtHome/hough_x87.ci"
#else
#warning ADDPHMD2HD_WLR_LOOP unset
#endif

RCSID( "$Id$");

#define HSMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define HSMIN(x,y) ( (x) < (y) ? (x) : (y) )
#define INIT_MEM(x) memset(&(x), 0, sizeof((x)))


#ifdef OUTPUT_TIMING
extern time_t clock0;
extern UINT4 nSFTs;
extern UINT4 nStacks;
extern UINT4 nSkyRefine;
#endif


/* comparison function for the toplist */
static int smallerHough(const void *a,const void *b) {
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


/* we point nearly all LocalHOUGH functions back to LAL */

#define LocalHOUGHComputeSizePar     LALHOUGHComputeSizePar
#define LocalHOUGHFillPatchGrid      LALHOUGHFillPatchGrid
#define LocalHOUGHParamPLUT          LALHOUGHParamPLUT
#define LocalHOUGHConstructPLUT      LALHOUGHConstructPLUT
#define LocalHOUGHConstructSpacePHMD LALHOUGHConstructSpacePHMD
#define LocalHOUGHWeighSpacePHMD     LALHOUGHWeighSpacePHMD
#define LocalHOUGHInitializeHT       LALHOUGHInitializeHT
#define LocalHOUGHupdateSpacePHMDup  LALHOUGHupdateSpacePHMDup
#define LocalHOUGHWeighSpacePHMD     LALHOUGHWeighSpacePHMD


/* possibly optimized local copies of LALHOUGH functions */

/* stupid MSC doesn't understand inline function (at least not in VS2003) */
#ifdef _MSC_VER
#define INLINE static
#else
#define INLINE inline
#endif

#ifdef __GNUC__
#define ALWAYS_INLINE __attribute__((always_inline))
#endif

static void
LocalHOUGHConstructHMT_W  (LALStatus                  *status, 
			   HOUGHMapTotal              *ht     , /**< The output hough map */
			   UINT8FrequencyIndexVector  *freqInd, /**< time-frequency trajectory */ 
			   PHMDVectorSequence         *phmdVS); /**< set of partial hough map derivatives */

static void
LocalHOUGHAddPHMD2HD_W    (LALStatus      *status, /**< the status pointer */
			   HOUGHMapDeriv  *hd,     /**< the Hough map derivative */
			   HOUGHphmd      *phmd);  /**< info from a partial map */ 

/* this is the only function that's actually changed for optimization */
INLINE void
LocalHOUGHAddPHMD2HD_Wlr  (LALStatus*    status,
			   HoughDT*      map,
			   HOUGHBorder** pBorderP,
			   INT4         length,
			   HoughDT       weight,
			   INT4         xSide, 
			   INT4         ySide) ALWAYS_INLINE;


/* this function is identical to LALComputeFstatHoughMap() in HierarchicalSearch.c,
   except for calling the LocalHOUGH* functions instead of the LALHOUGH* ones
*/
void
LocalComputeFstatHoughMap (LALStatus            *status,
			   SemiCohCandidateList *out,    /* output candidates */
			   HOUGHPeakGramVector  *pgV,    /* peakgram vector */
			   SemiCoherentParams   *params)
{
  /* hough structures */
  HOUGHMapTotal ht;
  HOUGHptfLUTVector   lutV; /* the Look Up Table vector*/
  PHMDVectorSequence  phmdVS;  /* the partial Hough map derivatives */
  UINT8FrequencyIndexVector freqInd; /* for trajectory in time-freq plane */
  HOUGHResolutionPar parRes;   /* patch grid information */
  HOUGHPatchGrid  patch;   /* Patch description */ 
  HOUGHParamPLUT  parLut;  /* parameters needed to build lut  */
  HOUGHDemodPar   parDem;  /* demodulation parameters */
  HOUGHSizePar    parSize; 

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

  toplist_t *houghToplist;

  INITSTATUS( status, "LocalComputeFstatHoughMap", rcsid );
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
  if ( pgV == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  
  if ( pgV->length == 0 ) {
    ABORT ( status, HIERARCHICALSEARCH_EVAL, HIERARCHICALSEARCH_MSGEVAL );
  }  
  if ( pgV->pg == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  
  if ( params == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  



  /* copy some parameters from peakgram vector */
  deltaF = pgV->pg->deltaF;
  nStacks = pgV->length;
  fBinIni = pgV->pg[0].fBinIni;
  fBinFin = pgV->pg[0].fBinFin;

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
  TRY ( LALGPStoFloat( status->statusPtr, &refTime, &refTimeGPS), status);

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

    TRY ( LALGPStoFloat ( status->statusPtr, &tMidStack, tsMid->data + k), status);
    timeDiffV->data[k] = tMidStack - refTime;
  }



  /*--------------- first memory allocation --------------*/
  /* look up table vector */
  lutV.length = nStacks;
  lutV.lut = NULL;
  lutV.lut = (HOUGHptfLUT *)LALCalloc(1,nStacks*sizeof(HOUGHptfLUT));
  if ( lutV.lut == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
  }  


  /* partial hough map derivative vector */
  phmdVS.length  = nStacks;

  {
    REAL8 maxTimeDiff, startTimeDiff, endTimeDiff;

    startTimeDiff = fabs(timeDiffV->data[0]);
    endTimeDiff = fabs(timeDiffV->data[timeDiffV->length - 1]);
    maxTimeDiff = HSMAX( startTimeDiff, endTimeDiff);

    /* set number of freq. bins for which LUTs will be calculated */
    /* this sets the range of residual spindowns values */
    /* phmdVS.nfSize  = 2*nfdotBy2 + 1; */
    phmdVS.nfSize  = 2 * floor((nfdot-1) * (REAL4)(dfdot * maxTimeDiff / deltaF) + 0.5f) + 1; 
  }

  phmdVS.deltaF  = deltaF;
  phmdVS.phmd = NULL;
  phmdVS.phmd=(HOUGHphmd *)LALCalloc( 1,phmdVS.length * phmdVS.nfSize *sizeof(HOUGHphmd));
  if ( phmdVS.phmd == NULL ) {
    ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
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
    create_toplist(&houghToplist, out->length, sizeof(SemiCohCandidate), smallerHough);
  }
  else { 
    /* if no toplist then use number of hough maps */
    INT4 numHmaps = (fBinFin - fBinIni + 1)*phmdVS.nfSize;
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
    	
    parRes.f0Bin =  fBin;      
    TRY( LocalHOUGHComputeSizePar( status->statusPtr, &parSize, &parRes ),  status );
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
    TRY( LocalHOUGHFillPatchGrid( status->statusPtr, &patch, &parSize ), status );
    
    /*------------- other memory allocation and settings----------------- */
    for(j=0; j<lutV.length; ++j){
      lutV.lut[j].maxNBins = maxNBins;
      lutV.lut[j].maxNBorders = maxNBorders;
      lutV.lut[j].border = (HOUGHBorder *)LALCalloc(1,maxNBorders*sizeof(HOUGHBorder));
      if ( lutV.lut[j].border == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }  

      lutV.lut[j].bin =	(HOUGHBin2Border *)LALCalloc(1,maxNBins*sizeof(HOUGHBin2Border));
      if ( lutV.lut[j].bin == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }  

      for (i=0; i<maxNBorders; ++i){
	lutV.lut[j].border[i].ySide = ySide;
	lutV.lut[j].border[i].xPixel = (COORType *)LALCalloc(1,ySide*sizeof(COORType));
	if ( lutV.lut[j].border[i].xPixel == NULL ) {
	  ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
	}  
      }
    }

    for(j = 0; j < phmdVS.length * phmdVS.nfSize; ++j){
      phmdVS.phmd[j].maxNBorders = maxNBorders;
      phmdVS.phmd[j].leftBorderP = (HOUGHBorder **)LALCalloc(1,maxNBorders*sizeof(HOUGHBorder *));
      if ( phmdVS.phmd[j].leftBorderP == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }  

      phmdVS.phmd[j].rightBorderP = (HOUGHBorder **)LALCalloc(1,maxNBorders*sizeof(HOUGHBorder *));
      if ( phmdVS.phmd[j].rightBorderP == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }  

      phmdVS.phmd[j].ySide = ySide;
      phmdVS.phmd[j].firstColumn = NULL;
      phmdVS.phmd[j].firstColumn = (UCHAR *)LALCalloc(1,ySide*sizeof(UCHAR));
      if ( phmdVS.phmd[j].firstColumn == NULL ) {
	ABORT ( status, HIERARCHICALSEARCH_ENULL, HIERARCHICALSEARCH_MSGENULL );
      }  
    }
    
    /*------------------- create all the LUTs at fBin ---------------------*/  
    for (j=0; j < (UINT4)nStacks; j++){  /* create all the LUTs */
      parDem.veloC.x = vel->data[3*j];
      parDem.veloC.y = vel->data[3*j + 1];
      parDem.veloC.z = vel->data[3*j + 2];      
      parDem.positC.x = pos->data[3*j];
      parDem.positC.y = pos->data[3*j + 1];
      parDem.positC.z = pos->data[3*j + 2];
      parDem.timeDiff = timeDiffV->data[j];

      /* calculate parameters needed for buiding the LUT */
      TRY( LocalHOUGHParamPLUT( status->statusPtr, &parLut, &parSize, &parDem), status);

      /* build the LUT */
      TRY( LocalHOUGHConstructPLUT( status->statusPtr, &(lutV.lut[j]), &patch, &parLut ), status);

    }
    

    
    /*--------- build the set of  PHMD centered around fBin -------------*/     
    phmdVS.fBinMin = fBin - phmdVS.nfSize/2;
    TRY( LocalHOUGHConstructSpacePHMD(status->statusPtr, &phmdVS, pgV, &lutV), status );
    TRY( LocalHOUGHWeighSpacePHMD(status->statusPtr, &phmdVS, params->weightsV), status);
    
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

    TRY( LocalHOUGHInitializeHT( status->statusPtr, &ht, &patch), status); /*not needed */
    
    /*  Search frequency interval possible using the same LUTs */
    fBinSearch = fBin;
    fBinSearchMax = fBin + parSize.nFreqValid - 1;
     
    /* Study all possible frequencies with one set of LUT */    
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

	  TRY( LocalHOUGHConstructHMT_W(status->statusPtr, &ht, &freqInd, &phmdVS),status );

	  /* get candidates */
	  if ( params->useToplist ) {
	    TRY(GetHoughCandidates_toplist( status->statusPtr, houghToplist, &ht, &patch, &parDem), status);
	  }
	  else {
	    TRY(GetHoughCandidates_threshold( status->statusPtr, out, &ht, &patch, &parDem, params->threshold), status);
	  }
	  
	  /* increment hough map index */ 	  
	  ++iHmap;
	  
	} /* end loop over spindown trajectories */

      } /* end of block for calculating total hough maps */
      

      /*------ shift the search freq. & PHMD structure 1 freq.bin -------*/
      ++fBinSearch;
      TRY( LocalHOUGHupdateSpacePHMDup(status->statusPtr, &phmdVS, pgV, &lutV), status );
      TRY( LocalHOUGHWeighSpacePHMD(status->statusPtr, &phmdVS, params->weightsV), status);      

    }   /* closing while loop over fBinSearch */

#ifdef OUTPUT_TIMING
    /* printf ("xside x yside = %d x %d = %d\n", parSize.xSide, parSize.ySide, parSize.xSide * parSize.ySide ); */
    nSkyRefine = parSize.xSide * parSize.ySide;
#endif
    
    fBin = fBinSearch;
    
    /*--------------  Free partial memory -----------------*/
    LALFree(patch.xCoor);
    LALFree(patch.yCoor);
    LALFree(ht.map);

    for (j=0; j<lutV.length ; ++j){
      for (i=0; i<maxNBorders; ++i){
	LALFree( lutV.lut[j].border[i].xPixel);
      }
      LALFree( lutV.lut[j].border);
      LALFree( lutV.lut[j].bin);
    }
    for(j=0; j<phmdVS.length * phmdVS.nfSize; ++j){
      LALFree( phmdVS.phmd[j].leftBorderP);
      LALFree( phmdVS.phmd[j].rightBorderP);
      LALFree( phmdVS.phmd[j].firstColumn);
    }
    
  } /* closing first while */

  
  /* free remaining memory */
  LALFree(ht.spinRes.data);
  LALFree(ht.spinDem.data);
  LALFree(ht.dFdot.data);
  LALFree(lutV.lut);
  LALFree(phmdVS.phmd);
  LALFree(freqInd.data);
  LALFree(parDem.spin.data);

  TRY( LALDDestroyVector( status->statusPtr, &timeDiffV), status);

  /* copy toplist candidates to output structure if necessary */
  if ( params->useToplist ) {
    for ( k=0; k<houghToplist->elems; k++) {
      out->list[k] = *((SemiCohCandidate *)(toplist_elem(houghToplist, k)));
    }
    out->nCandidates = houghToplist->elems;
    free_toplist(&houghToplist);
  }

  DETATCHSTATUSPTR (status);
  RETURN(status);

}



/* this function is identical to LALHOUGHConstructHMT_W in DriveHough.c,
   except for that it calls LocalHOUGHAddPHMD2HD_W instead of LALHOUGHAddPHMD2HD_W
*/
static void
LocalHOUGHConstructHMT_W (LALStatus                  *status, 
			  HOUGHMapTotal              *ht, /**< The output hough map */
			  UINT8FrequencyIndexVector  *freqInd, /**< time-frequency trajectory */ 
			  PHMDVectorSequence         *phmdVS) /**< set of partial hough map derivatives */
{

  UINT4    k,j;
  UINT4    breakLine;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  UINT8    fBinMin;   /* present minimum frequency bin */ 
  INT8     fBin;      /* present frequency bin */
  UINT2    xSide,ySide;
 
  HOUGHMapDeriv hd; /* the Hough map derivative */

  /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHConstructHMT_W", rcsid);
  ATTATCHSTATUSPTR (status); 

  /*   Make sure the arguments are not NULL: */ 
  ASSERT (phmdVS,  status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (ht,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (freqInd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd,  status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (freqInd->data, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (freqInd->length == phmdVS->length, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (freqInd->deltaF == phmdVS->deltaF, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements  */
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  /* -------------------------------------------   */
  
   /* Make sure the ht map contains some pixels */
  ASSERT (ht->xSide, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (ht->ySide, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);

  length = phmdVS->length;
  nfSize = phmdVS->nfSize; 
  
  fBinMin = phmdVS->fBinMin; /* initial frequency value  od the cilinder*/
  
  breakLine = phmdVS->breakLine;

  /* number of physical pixels */
  xSide = ht->xSide;
  ySide = ht->ySide;
  
  /* Make sure initial breakLine is in [0,nfSize)  */
  ASSERT ( breakLine < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
  
  /* -------------------------------------------   */
  
  /* Initializing  hd map and memory allocation */
  hd.xSide = xSide;
  hd.ySide = ySide;
  hd.map = (HoughDT *)LALMalloc(ySide*(xSide+1)*sizeof(HoughDT));
  if (hd. map == NULL) {
    ABORT( status, LALHOUGHH_EMEM, LALHOUGHH_MSGEMEM); 
  }

  /* -------------------------------------------   */
 
  TRY( LALHOUGHInitializeHD(status->statusPtr, &hd), status);
  for ( k=0; k<length; ++k ){ 
    /* read the frequency index and make sure is in the proper interval*/
    fBin =freqInd->data[k] -fBinMin;

    ASSERT ( fBin < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
    ASSERT ( fBin >= 0,     status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
 
    /* find index */
    j = (fBin + breakLine) % nfSize;

    /* Add the corresponding PHMD to HD */
    TRY( LocalHOUGHAddPHMD2HD_W(status->statusPtr,
				&hd, &(phmdVS->phmd[j*length+k]) ), status);
  }

  TRY( LALHOUGHIntegrHD2HT(status->statusPtr, ht, &hd), status);
  
  /* Free memory and exit */
  LALFree(hd.map);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}



/* this function is derived from LALHOUGHAddPHMD2HD_W in HoughMap.c.
   The two originally almost identical loops were put into the inline function
   LocalHOUGHAddPHMD2HD_Wlr() that is augmented with prefetching code to
   speed up memory access. The parameterization is such that the originally
   two loops become actually identical in LocalHOUGHAddPHMD2HD_Wlr
   (leftBorder <-> rightBorder, lengthLeft <-> lengthRight, weight <-> -weight)
*/
static void
LocalHOUGHAddPHMD2HD_W (LALStatus      *status, /**< the status pointer */
			HOUGHMapDeriv  *hd,     /**< the Hough map derivative */
			HOUGHphmd      *phmd)   /**< info from a partial map */ 
{

  INT2     k;
  UINT2    xSide,ySide;
  HoughDT  weight;

  INITSTATUS (status, "LALHOUGHAddPHMD2HD_W", rcsid);
  ATTATCHSTATUSPTR (status); 

  /*   Make sure the arguments are not NULL: */ 
  ASSERT (hd,   status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  ASSERT (phmd, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);

  /* Make sure the map contains some pixels */
  ASSERT (hd->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (hd->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);

  /* aliases */
  weight = phmd->weight;
  xSide = hd->xSide;
  ySide = hd->ySide;
  
  /* first column correction */
  for ( k=0; k< ySide; ++k ){
    hd->map[k*(xSide+1) + 0] += phmd->firstColumn[k] * weight;
  }

  /* left borders =>  increase according to weight */
  TRY ( LocalHOUGHAddPHMD2HD_Wlr (status,
				  hd->map,
				  phmd->leftBorderP,
				  phmd->lengthLeft,
				  weight,
				  xSide,
				  ySide), status );
  
  /* right borders => decrease according to weight */
  TRY ( LocalHOUGHAddPHMD2HD_Wlr (status,
				  hd->map,
				  phmd->rightBorderP,
				  phmd->lengthRight,
				  - weight,
				  xSide,
				  ySide), status );
  
  /* cleanup */
  DETATCHSTATUSPTR (status);
  
  /* normal exit */
  RETURN (status);
}



/* prefetch compiler directives */
#if EAH_HOUGH_PREFETCH > EAH_HOUGH_PREFETCH_NONE
#if defined(__INTEL_COMPILER) ||  defined(_MSC_VER)
// not tested yet with icc or MS Visual C 
#include "xmmintrin.h"
#define PREFETCH(a) _mm_prefetch((char *)(void *)(a),_MM_HINT_T0)
#elif defined(__GNUC__)
#define PREFETCH(a) __builtin_prefetch(a)
#else
#define PREFETCH(a) a
#endif
#else
#define PREFETCH(a) a
#endif

#if ( (EAH_HOUGH_PREFETCH == EAH_HOUGH_PREFETCH_X87) || (EAH_HOUGH_PREFETCH == EAH_HOUGH_PREFETCH_AMD64) )

INLINE void __attribute__((always_inline))
LocalHOUGHAddPHMD2HD_Wlr  (LALStatus*    status,
			   HoughDT*      map,
			   HOUGHBorder** pBorderP,
			   INT4         length,
			   HoughDT       weight,
			   INT4         xSide, 
			   INT4         ySide) {

  INT4  xSideP1 = xSide +1; /* avoid 16 bit types */
  INT4  yLower, yUpper;
  INT4  k;
  COORType     *xPixel;
  HOUGHBorder* borderP;
   
  for (k=0; k< length; ++k){

    borderP = pBorderP[k];
    xPixel =  &( (*borderP).xPixel[0] );

    yLower = (*borderP).yLower;
    yUpper = (*borderP).yUpper;

    if(k < length-1) {
	INT4 ylkp1 = pBorderP[k+1]->yLower;
	PREFETCH(&(pBorderP[k+1]->xPixel[ylkp1]));
    } 	

    if(k < length-2) {
	PREFETCH(pBorderP[k+2]);
    } 	
   
    if (yLower < 0) {
      fprintf(stderr,"WARNING: Fixing yLower (%d -> 0) [HoughMap.c %d]\n",
	      yLower, __LINE__);
      yLower = 0;
    }
    if (yUpper >= ySide) {
      fprintf(stderr,"WARNING: Fixing yUpper (%d -> %d) [HoughMap.c %d]\n",
	      yUpper, ySide-1, __LINE__);
      yUpper = ySide - 1;
    }

    ADDPHMD2HD_WLR_LOOP(xPixel,yLower,yUpper,xSideP1,map,weight)

  };
    
};


#elif EAH_HOUGH_PREFETCH == EAH_HOUGH_PREFETCH_DIRECT

#ifndef EAH_HOUGH_BATCHSIZE_LOG2
#define EAH_HOUGH_BATCHSIZE_LOG2 2
#endif
#define EAH_HOUGH_BATCHSIZE (1 << EAH_HOUGH_BATCHSIZE_LOG2)

#define CHECK_INDEX(IDX,OFFSET) \
      if ((IDX < 0) || ( IDX >= ySide*(xSide+1))|| xPixel[OFFSET] < 0 || xPixel[OFFSET] >= xSideP1) { \
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d xSide:%d\n", \
		__FILE__,__LINE__,IDX,ySide*(xSide+1),OFFSET,xPixel[OFFSET],xSide ); \
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE); \
      } 

INLINE void
LocalHOUGHAddPHMD2HD_Wlr (LALStatus*    status,
			  HoughDT*      map,
			  HOUGHBorder** pBorderP,
			  INT4         length,
			  HoughDT      weight,
			  INT4         xSide, 
			  INT4         ySide)
{

  INT4        k,j;
  INT4        yLower, yUpper;
  register    HoughDT    tempM0,tempM1,tempM2,tempM3;
  INT4        sidx,sidx0,sidx1,sidx2,sidx3,sidxBase, sidxBase_n; /* pre-calcuted array index for sanity check */
  INT4        c_c,c_n,offs;
  COORType    *xPixel;
  HOUGHBorder *borderP;
  INT4        xSideP1 = xSide +1; /* avoid 16 bit types */
  INT4        xSideP1_2=xSideP1+xSideP1;
  INT4        xSideP1_3=xSideP1_2+xSideP1;

  HoughDT     *pf_addr[8]; 

  for (k=0; k< length; ++k) {

    /*  Make sure the arguments are not NULL: (Commented for performance) */ 
    /*  ASSERT (phmd->leftBorderP[k], status, HOUGHMAPH_ENULL,
	HOUGHMAPH_MSGENULL); */

    borderP = pBorderP[k];
    xPixel =  &( (*borderP).xPixel[0] );

    yLower = borderP->yLower;
    yUpper = borderP->yUpper;

    sidxBase=yLower*xSideP1;
    sidxBase_n = sidxBase+(xSideP1 << EAH_HOUGH_BATCHSIZE_LOG2);


    if(k < length-1) {
	INT4 ylkp1 = pBorderP[k+1]->yLower;
	PREFETCH(&(pBorderP[k+1]->xPixel[ylkp1]));
    } 	

    if(k < length-2) {
	PREFETCH(pBorderP[k+2]);
    } 	

   
    if (yLower < 0) {
      fprintf(stderr,"WARNING: Fixing yLower (%d -> 0) [HoughMap.c %d]\n",
	      yLower, __LINE__);
      yLower = 0;
    }
    if (yUpper >= ySide) {
      fprintf(stderr,"WARNING: Fixing yUpper (%d -> %d) [HoughMap.c %d]\n",
	      yUpper, ySide-1, __LINE__);
      yUpper = ySide - 1;
    }


    /* fill first prefetch address array entries */
    c_c =0;
    c_n =EAH_HOUGH_BATCHSIZE;

    offs = yUpper - yLower+1;
    if (offs > EAH_HOUGH_BATCHSIZE) {
	offs = EAH_HOUGH_BATCHSIZE; 
    }	
	
    	
    for(j=yLower; j < yLower+offs; j++) {
        sidx0=xPixel[j]+ j*xSideP1;
        PREFETCH(pf_addr[c_c++] = map + sidx0);
#ifndef LAL_NDEBUG
        CHECK_INDEX(sidx0,j);
#endif
    }		
		
    c_c=0;
    for(j=yLower; j<=yUpper-(2*EAH_HOUGH_BATCHSIZE-1);j+=EAH_HOUGH_BATCHSIZE){

      sidx0 = xPixel[j+EAH_HOUGH_BATCHSIZE]+sidxBase_n;; 
      sidx1 = xPixel[j+EAH_HOUGH_BATCHSIZE+1]+sidxBase_n+xSideP1;
      sidx2 = xPixel[j+EAH_HOUGH_BATCHSIZE+2]+sidxBase_n+xSideP1_2;
      sidx3 = xPixel[j+EAH_HOUGH_BATCHSIZE+3]+sidxBase_n+xSideP1_3;;
	
      PREFETCH(xPixel +(j+(EAH_HOUGH_BATCHSIZE+EAH_HOUGH_BATCHSIZE)));

      PREFETCH(pf_addr[c_n] = map + sidx0);
      PREFETCH(pf_addr[c_n+1] = map + sidx1);
      PREFETCH(pf_addr[c_n+2] = map + sidx2);
      PREFETCH(pf_addr[c_n+3] = map + sidx3);

#ifndef LAL_NDEBUG 
      CHECK_INDEX(sidx0,j+EAH_HOUGH_BATCHSIZE);
      CHECK_INDEX(sidx1,j+EAH_HOUGH_BATCHSIZE+1);
      CHECK_INDEX(sidx2,j+EAH_HOUGH_BATCHSIZE+2);
      CHECK_INDEX(sidx3,j+EAH_HOUGH_BATCHSIZE+3);
#endif

      tempM0 = *(pf_addr[c_c]) +weight;
      tempM1 = *(pf_addr[c_c+1]) +weight;
      tempM2 = *(pf_addr[c_c+2]) +weight;
      tempM3 = *(pf_addr[c_c+3]) +weight;

      sidxBase = sidxBase_n;
      sidxBase_n+=xSideP1 << EAH_HOUGH_BATCHSIZE_LOG2;

      (*(pf_addr[c_c]))=tempM0;
      (*(pf_addr[c_c+1]))=tempM1;
      (*(pf_addr[c_c+2]))=tempM2;
      (*(pf_addr[c_c+3]))=tempM3;

      c_c ^= EAH_HOUGH_BATCHSIZE;
      c_n ^= EAH_HOUGH_BATCHSIZE;
    }

    sidxBase=j*xSideP1;
    for(; j<=yUpper;++j){
      sidx = sidxBase + xPixel[j];
#ifndef LAL_NDEBUG
      CHECK_INDEX(sidx0,j);
#endif
      map[sidx] += weight;
      sidxBase+=xSideP1;
    }

  }
}


#else /* EAH_HOUGH_PREFETCH */

/* original generic version w/o prefetching */

INLINE void
LocalHOUGHAddPHMD2HD_Wlr (LALStatus*    status,
			  HoughDT*      map,
			  HOUGHBorder** pBorderP,
			  INT4         length,
			  HoughDT      weight,
			  INT4         xSide, 
			  INT4         ySide)
{
  INT4        k,j;
  INT4        yLower, yUpper;
  COORType    *xPixel;
  HOUGHBorder *borderP;
  INT4        sidx;

  for (k=0; k< length; ++k) {

    /* local aliases */
    borderP = pBorderP[k];
    yLower = (*borderP).yLower;
    yUpper = (*borderP).yUpper;
    xPixel =  &((*borderP).xPixel[0]);
   
    /* check boundary conditions */
    if (yLower < 0) {
      fprintf(stderr,"WARNING: Fixing yLower (%d -> 0) [HoughMap.c %d]\n",
	      yLower, __LINE__);
      yLower = 0;
    }
    if (yUpper >= ySide) {
      fprintf(stderr,"WARNING: Fixing yUpper (%d -> %d) [HoughMap.c %d]\n",
	      yUpper, ySide-1, __LINE__);
      yUpper = ySide - 1;
    }

    /* increase / decrease according to weight */
    for(j = yLower; j <= yUpper; j++){
      sidx = j*(xSide+1) + xPixel[j];
      if ((sidx < 0) || (sidx >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j,xPixel[j] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
      map[sidx] += weight;
    }
  }
}

#endif /* EAH_HOUGH_PREFETCH */
