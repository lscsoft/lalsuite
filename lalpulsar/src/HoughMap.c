/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes, Bernd Machenschalk
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

#include <lal/HoughMap.h>

/*
 * The functions that make up the guts of this module
 */

/**
 * This function initializes the Hough map derivative
 * space  HOUGHMapDeriv *hd to zero. Note that the length of the map
 * hd->map should be hd->ySide * (hd->xSide + 1).
 */
void LALHOUGHInitializeHD (LALStatus      *status,
			  HOUGHMapDeriv   *hd) /* the Hough map derivative */
{

  INT4     k, maxk;
  HoughDT  *pointer;

   /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (hd,    status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  /* Make sure the map contains some pixels */
  ASSERT (hd->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (hd->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);

  /* -------------------------------------------   */

  /* initializing the Hough map derivative space */
  pointer = &( hd->map[0]);
  maxk    = hd->ySide*(hd->xSide+1);

  for ( k=0; k< maxk; ++k ){
    *pointer = 0;
    ++pointer;
  }

  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}

/**
 * This function initializes the total Hough map
 * HOUGHMapTotal *ht  to zero and checks consistency between
 * the number of physical pixels in the
 * map  and  those given by the grid information structure
 * HOUGHPatchGrid *patch.
 */
void LALHOUGHInitializeHT (LALStatus      *status,
			   HOUGHMapTotal   *ht,     /* the total Hough map */
			   HOUGHPatchGrid  *patch) /* patch information */
{

  INT4     k,maxk;
  HoughTT  *pointer;

   /* --------------------------------------------- */
  INITSTATUS(status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (ht, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  ASSERT (patch, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  /* Make sure the map contains some pixels */
  ASSERT (ht->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (ht->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (patch->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (patch->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  /*   Check consistency between patch and ht size  (size mismatch) */
  ASSERT (ht->xSide == patch->xSide, status, HOUGHMAPH_ESZMM, HOUGHMAPH_MSGESZMM);
  ASSERT (ht->ySide == patch->ySide, status, HOUGHMAPH_ESZMM, HOUGHMAPH_MSGESZMM);

  /* -------------------------------------------   */

  /* number of physical pixels */
  maxk = ht->ySide * ht->xSide;

  /* initializing the Hough map space */
   pointer = &(ht->map[0]);
   for ( k=0; k< maxk; ++k ){
     *pointer = 0;
     ++pointer;
   }

  (void)patch;

  /* normal exit */
  RETURN (status);
}


/**
 * Given an initial Hough map derivative HOUGHMapDeriv *hd and a representation
 * of a phmd HOUGHphmd *phmd, the function  LALHOUGHAddPHMD2HD() accumulates
 * the partial Hough map derivative *phmd to *hd by adding +1 or
 * -1 to the pixels corresponding to the left or right borders respectively.
 * It takes into account corrections due to border effects as well.
 */
void LALHOUGHAddPHMD2HD (LALStatus      *status, 	/**< the status pointer */
			 HOUGHMapDeriv  *hd,  		/**< the Hough map derivative */
			 HOUGHphmd      *phmd) 		/**< info from a partial map */
{

  INT2     k,j;
  INT2     yLower, yUpper;
  UINT2    lengthLeft,lengthRight, xSide,ySide;
  COORType     *xPixel;
  HOUGHBorder  *borderP;

   /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (hd,   status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  ASSERT (phmd, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  /* -------------------------------------------   */
  /* Make sure the map contains some pixels */
  ASSERT (hd->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (hd->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);

  xSide = hd->xSide;
  ySide = hd->ySide;

  /* first column correction */
  for ( k=0; k< ySide; ++k ){
    hd->map[k*(xSide+1) + 0] += phmd->firstColumn[k];
  }

  lengthLeft = phmd->lengthLeft;
  lengthRight= phmd->lengthRight;

  /* left borders =>  +1 increase */
  for (k=0; k< lengthLeft; ++k){

    /*  Make sure the arguments are not NULL: (Commented for performance) */
    /*  ASSERT (phmd->leftBorderP[k], status, HOUGHMAPH_ENULL,
	HOUGHMAPH_MSGENULL); */

    borderP = phmd->leftBorderP[k];

    yLower = (*borderP).yLower;
    yUpper = (*borderP).yUpper;
    xPixel =  &( (*borderP).xPixel[0] );


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

    for(j=yLower; j<=yUpper;++j){
      hd->map[j *(xSide+1) + xPixel[j] ] += 1;
    }
  }

  /* right borders =>  -1 decrease */
  for (k=0; k< lengthRight; ++k){

    /*  Make sure the arguments are not NULL: (Commented for performance) */
    /*  ASSERT (phmd->rightBorderP[k], status, HOUGHMAPH_ENULL,
	HOUGHMAPH_MSGENULL); */

    borderP = phmd->rightBorderP[k];

    yLower = (*borderP).yLower;
    yUpper = (*borderP).yUpper;
    xPixel =  &( (*borderP).xPixel[0] );


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

    for(j=yLower; j<=yUpper;++j){
      hd->map[j*(xSide+1) + xPixel[j] ] -= 1;
    }
  }


  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}


/**
 * Adds a hough map derivative into a total hough map derivative taking into
 * account the weight of the partial hough map
 */
void LALHOUGHAddPHMD2HD_W (LALStatus      *status, 	/**< the status pointer */
			   HOUGHMapDeriv  *hd,  	/**< the Hough map derivative */
			   HOUGHphmd      *phmd) 	/**< info from a partial map */
{

  INT2     k,j;
  INT2     yLower, yUpper;
  UINT2    lengthLeft,lengthRight, xSide,ySide;
  COORType     *xPixel;
  HOUGHBorder  *borderP;
  HoughDT    weight;
  INT4       sidx; /* pre-calcuted array index for sanity check */

   /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (hd,   status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  ASSERT (phmd, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  /* -------------------------------------------   */
  /* Make sure the map contains some pixels */
  ASSERT (hd->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (hd->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);

  weight = phmd->weight;

  xSide = hd->xSide;
  ySide = hd->ySide;

  /* first column correction */
  for ( k=0; k< ySide; ++k ){
    hd->map[k*(xSide+1) + 0] += phmd->firstColumn[k] * weight;
  }

  lengthLeft = phmd->lengthLeft;
  lengthRight= phmd->lengthRight;

  /* left borders =>  increase according to weight*/
  for (k=0; k< lengthLeft; ++k){

    /*  Make sure the arguments are not NULL: (Commented for performance) */
    /*  ASSERT (phmd->leftBorderP[k], status, HOUGHMAPH_ENULL,
	HOUGHMAPH_MSGENULL); */

    borderP = phmd->leftBorderP[k];

    yLower = (*borderP).yLower;
    yUpper = (*borderP).yUpper;
    xPixel =  &( (*borderP).xPixel[0] );

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

    for(j=yLower; j<=yUpper;++j){
      sidx = j *(xSide+1) + xPixel[j];
      if ((sidx < 0) || (sidx >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j,xPixel[j] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
      hd->map[sidx] += weight;
    }
  }

  /* right borders => decrease according to weight*/
  for (k=0; k< lengthRight; ++k){

    /*  Make sure the arguments are not NULL: (Commented for performance) */
    /*  ASSERT (phmd->rightBorderP[k], status, HOUGHMAPH_ENULL,
	HOUGHMAPH_MSGENULL); */

    borderP = phmd->rightBorderP[k];

    yLower = (*borderP).yLower;
    yUpper = (*borderP).yUpper;
    xPixel =  &( (*borderP).xPixel[0] );

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

    for(j=yLower; j<=yUpper;++j){
      sidx = j*(xSide+1) + xPixel[j];
      if ((sidx < 0) || (sidx >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j,xPixel[j] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
      hd->map[sidx] -= weight;
    }
  }


  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}

/**
 * This function constructs a total Hough map
 * HOUGHMapTotal *ht from its derivative HOUGHMapDeriv *hd by
 * integrating each row (x-direction).
 */
void LALHOUGHIntegrHD2HT (LALStatus       *status,
			  HOUGHMapTotal   *ht,     	/* the total Hough map */
			  HOUGHMapDeriv   *hd) 		/* the Hough map derivative */
{

  INT2    i,j;
  UINT2   xSide,ySide;
  HoughTT accumulator;

   /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (hd, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  ASSERT (ht, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  /* Make sure the map contains some pixels */
  ASSERT (hd->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (hd->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (ht->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (ht->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  /*   Check consistency between hd and ht size  (size mismatch) */
  ASSERT (ht->xSide == hd->xSide, status, HOUGHMAPH_ESZMM, HOUGHMAPH_MSGESZMM);
  ASSERT (ht->ySide == hd->ySide, status, HOUGHMAPH_ESZMM, HOUGHMAPH_MSGESZMM);
  /* -------------------------------------------   */

  /* number of physical pixels */
  xSide = ht->xSide;
  ySide = ht->ySide;

  /* To construct the Hough map from the derivative,
     the latter must be integrated row-wise (x direction) */

  /* Loop on the rows */
  for (j=0; j< ySide; ++j){
    accumulator = 0;
    for ( i=0; i<xSide; ++i){
      ht->map[j*xSide +i] = ( accumulator += hd->map[j*(xSide+1) +i]);
    }
  }



  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}

/**  Find source sky location given stereographic coordinates indexes */
void LALStereo2SkyLocation (LALStatus  *status,
         REAL8UnitPolarCoor *sourceLocation, /* output*/
	 UINT2              xPos,
	 UINT2              yPos,
	 HOUGHPatchGrid    *patch,
	 HOUGHDemodPar     *parDem)
{

  REAL8Cart2Coor        sourceProjected;
  REAL8UnitPolarCoor    sourceRotated;
  REAL8UnitPolarCoor    skyPatchCenter;
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (sourceLocation, status, HOUGHMAPH_ENULL,HOUGHMAPH_MSGENULL);
  ASSERT (patch , status, HOUGHMAPH_ENULL,HOUGHMAPH_MSGENULL);
  ASSERT (parDem, status, HOUGHMAPH_ENULL,HOUGHMAPH_MSGENULL);

  sourceProjected.x = patch->xCoor[xPos];
  sourceProjected.y = patch->yCoor[yPos];

  skyPatchCenter.alpha = parDem->skyPatch.alpha;
  skyPatchCenter.delta = parDem->skyPatch.delta;

  /* invert the stereographic projection for a point on the projected plane */
  TRY( LALStereoInvProjectCart( status->statusPtr,
				&sourceRotated, &sourceProjected ), status );

  /* undo roation in case the patch is not centered at the south pole */
  TRY( LALInvRotatePolarU( status->statusPtr,
       sourceLocation, &sourceRotated, &skyPatchCenter ), status );

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

