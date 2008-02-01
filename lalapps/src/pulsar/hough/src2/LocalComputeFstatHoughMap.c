/*
   This is a local copy of ComputeFstatHoughMap() of HierarchicalSearch.c
   See that file for details (short documentation, copyright etc.

   $Id$
*/

#include"HierarchicalSearch.h"

void LocalComputeFstatHoughMap(LALStatus *status,
			       SemiCohCandidateList  *out,   /* output candidates */
			       HOUGHPeakGramVector *pgV, /* peakgram vector */
			       SemiCoherentParams *params)
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
  UINT8Vector hist; /* histogram vector */ 
  UINT8Vector histTotal; /* total histogram vector */
  HoughStats stats; /* statistics struct */
  CHAR *fileStats = NULL;
  FILE *fpStats = NULL;

  toplist_t *houghToplist;

  INITSTATUS( status, "ComputeFstatHoughMap", rcsid );
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
    REAL8UnitPolarCoor sourceLocation;
    	
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
