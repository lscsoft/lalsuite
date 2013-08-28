/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes
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

/*-----------------------------------------------------------------------
 *
 * File Name: PatchGrid.c
 *
 * Authors: Sintes, A.M., Krishnan, B.
 *
 *
 * History:   Created by Sintes May 14, 2001
 *            Modified by Badri Krishnan Feb 2003
 *            Modified by Sintes May 2003
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * PatchGrid.c
 *
 * SYNOPSIS
 *
 * DESCRIPTION
 *
 * DIAGNOSTICS
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Alicia Sintes, Badri Krishnan
 * \file
 * \ingroup LUT_h
 * \brief Function for tiling  the sky-patch (on the projected plane).
 *
 * ### Description ###
 *
 * This  is a \c provisionalfinal now ? routine for tiling a  sky-pacth
 * on the projected plane.
 *
 * Patch size specified by user
 *
 * == doc needs to be updated ==
 *
 * The reason to call it  \c provisional  is because
 * the size of the patch depends on the grid used in the
 * demodulation stage. Neighbour sky-patches should not be separated
 * nor overlapping too much.
 * Here for setting the patch size, we  consider only \f$v_{epicycle}\f$,
 * the frequency \c f0 and  \c deltaF so that the ` longitudinal'  size
 * of the patch is given by <tt>side == deltaF/f0 * c/v_epi</tt>.
 * By taking \c f0 to be the maximun frequency considered in that step,
 * the patch-size is valid for a whole frequency range.\   \
 *
 * Given input parameters,  the function LALHOUGHPatchGrid() provides
 * patch information.
 *
 * The input <tt>*in1</tt> is a structure of type  \c HOUGHResolutionPar containing
 * some resolution parameters such as:
 * <tt>in1->f0</tt> a frequency, <tt>in1->deltaF</tt> the frequency resolution, and
 * <tt>in1->minWidthRatio</tt>  the ratio between the minimum  annulus width
 * for this search and the minimun  annulus width for  1 year integration time.
 * This value should be in the interval  [1.0, 25.0].
 *
 */


#include <lal/LUT.h>


/*
 * The functions that make up the guts of this module
 */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALHOUGHComputeSizePar (LALStatus  *status, /* demodulated case */
                   HOUGHSizePar        *out,
                   HOUGHResolutionPar  *in1
                   )
{ 


  INT8    f0Bin;        /* corresponding freq. bin  */
  REAL8   deltaF;    /* frequency resolution  df=1/TCOH*/
  REAL8   pixelFactor; /* number of pixel that fit in the thinnest annulus*/
  REAL8   pixErr;   /* for validity of LUT as PIXERR */
  REAL8   linErr;   /* as LINERR circle ->line */
  REAL8   vTotC;    /* estimate value of v-total/C as VTOT */

  REAL8   deltaX; /* pixel size in the projected plane */
  REAL8   deltaY;
  UINT2   xSide;    /* number of pixels in the x direction (projected plane)*/
  UINT2   ySide;    /* number of pixels in the y direction */
  UINT2   maxSide;    /* number of pixels in the y direction */

  UINT2   maxNBins;    /* maximum number of bins affecting the patch. For
                               memory allocation */
  REAL8   patchSizeX;  /* size of sky patch projected */
  REAL8   patchSizeY;  /* size of sky patch prijected */
  REAL8   patchSkySizeX;     /* Size of sky patch in radians */
  REAL8   patchSkySizeY;

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in1 , status, LUTH_ENULL, LUTH_MSGENULL);

  patchSkySizeX = in1->patchSkySizeX;
  patchSkySizeY = in1->patchSkySizeY;

 /* Make sure the user chose a sky patch smaller that pi of the sky */
  ASSERT (patchSkySizeX > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (patchSkySizeX <= LAL_PI,status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (patchSkySizeY > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (patchSkySizeY <= LAL_PI,status, LUTH_EVAL, LUTH_MSGEVAL);

  pixelFactor = in1->pixelFactor;
  pixErr = in1->pixErr;
  linErr = in1->linErr;

  /* Make sure the parameters make sense */
  ASSERT (pixelFactor > 0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (pixErr < 1.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (linErr < 1.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (pixErr > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (linErr > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);

  f0Bin  = in1->f0Bin;
  deltaF = in1->deltaF;
  vTotC  = in1->vTotC;

  ASSERT (f0Bin  > 0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (deltaF > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (vTotC  > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (vTotC  < 0.001, status, LUTH_EVAL, LUTH_MSGEVAL);

  /*  ************ THIS IS FOR THE DEMODULATED CASE ONLY ***********   */
  out->deltaF = deltaF;
  out->f0Bin = f0Bin;

  deltaX = deltaY = 1.0/(vTotC *pixelFactor * f0Bin);
  out->deltaX = deltaX;
  out->deltaY = deltaY;

  patchSizeX = 4.0 * tan(0.25*patchSkySizeX);
  patchSizeY = 4.0 * tan(0.25*patchSkySizeY);
  xSide = out->xSide = ceil( patchSizeX/deltaX );
  ySide = out->ySide = ceil( patchSizeY/deltaY );
  maxSide =  MAX( xSide , ySide );

  /* the max number of bins that can fit in the diagonal and a bit more
      1.41->1.5  */
  maxNBins = out->maxNBins = ceil( (1.5* maxSide)/pixelFactor);

  /* maximum number of borders affecting the patch +1. Each Bin has up to 4
  borders, but they are common (they share 2 with the next bin). Depending on
  the orientation could be equal to maxNBins. In other cases twice  maxNBins.
  Here we are very conservative*/

  out->maxNBorders = 1 + 2*maxNBins;

  /* LUT validity */
  /* From equation: nFreqValid = pixErr/(pixelFactor* vTotC * 0.5 * maxSide * deltaX);
      or the same is */
  out->nFreqValid = (pixErr * 2 * f0Bin)/(pixelFactor * maxNBins);

  /* max. angle (rad.) from the pole to consider a circle as a line in the projected plane */
  /* epsilon ~ 2/radius circle; radius ~h*h/(2b), epsilon = 4b/(h*h) */

  out->epsilon=8.0* linErr/(maxSide*deltaX*maxSide) ;

 /* -------------------------------------------   */


  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALHOUGHComputeNDSizePar (LALStatus  *status, /* non-demod. case */
                   HOUGHSizePar        *out,
                   HOUGHResolutionPar  *in1
                   )
{ 

  INT8    f0Bin;        /* corresponding freq. bin  */
  REAL8   deltaF;    /* frequency resolution  df=1/TCOH*/
  REAL8   pixelFactor; /* number of pixel that fit in the thinnest annulus*/
  REAL8   pixErr;   /* for validity of LUT as PIXERR */
  REAL8   linErr;   /* as LINERR circle ->line */
  REAL8   vTotC;    /* estimate value of v-total/C as VTOT */

  REAL8   deltaX; /* pixel size in the projected plane */
  REAL8   deltaY;
  UINT2   xSide;    /* number of pixels in the x direction (projected plane)*/
  UINT2   ySide;    /* number of pixels in the y direction */
  UINT2   maxDopplerBin;
  UINT2   maxSide;    /* number of pixels in the y direction */

  UINT2   maxNBins1,maxNBins2;
  UINT2   maxNBins;    /* maximum number of bins affecting the patch. For
                               memory allocation */
  REAL8   patchSizeX;  /* size of sky patch projected */
  REAL8   patchSizeY;  /* size of sky patch prijected */
  REAL8   patchSkySizeX;     /* Size of sky patch in radians */
  REAL8   patchSkySizeY;

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in1 , status, LUTH_ENULL, LUTH_MSGENULL);

  patchSkySizeX = in1->patchSkySizeX;
  patchSkySizeY = in1->patchSkySizeY;

 /* Make sure the user chose a sky patch smaller that pi of the sky */
  ASSERT (patchSkySizeX > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (patchSkySizeX <= LAL_PI,status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (patchSkySizeY > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (patchSkySizeY <= LAL_PI,status, LUTH_EVAL, LUTH_MSGEVAL);

  pixelFactor = in1->pixelFactor;
  pixErr = in1->pixErr;
  linErr = in1->linErr;

  /* Make sure the parameters make sense */
  ASSERT (pixErr < 1.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (linErr < 1.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (pixErr > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (linErr > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);

  f0Bin  = in1->f0Bin;
  deltaF = in1->deltaF;
  vTotC  = in1->vTotC;

  ASSERT (f0Bin  > 0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (deltaF > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (vTotC  > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (vTotC  < 0.001, status, LUTH_EVAL, LUTH_MSGEVAL);

  /*  ************ THIS IS FOR THE *NON* DEMODULATED CASE ONLY ***********   */
  out->deltaF = deltaF;
  out->f0Bin = f0Bin;

  deltaX = deltaY = 1.0/(vTotC *pixelFactor * f0Bin);
  out->deltaX = deltaX;
  out->deltaY = deltaY;

  patchSizeX = 4.0 * tan(0.25*patchSkySizeX);
  patchSizeY = 4.0 * tan(0.25*patchSkySizeY);
  xSide = out->xSide = ceil( patchSizeX/deltaX );
  ySide = out->ySide = ceil( patchSizeY/deltaY );
  maxSide =  MAX( xSide , ySide );

  /* the max number of bins that can fit in the full sky  */
  maxDopplerBin = floor( f0Bin * vTotC +0.5);
  maxNBins1 = 1+ 2* maxDopplerBin;

  /* estimate of the max number of bins that can fit in the patch
  (over-estimated) as the max number of bins that can fit in the diagonal
  and a bit more 1.41->1.5  */
  maxNBins2 = ceil( (1.5* maxSide)/pixelFactor);

  maxNBins = out->maxNBins = MIN( maxNBins1 , maxNBins2 );

  /* maximum number of borders affecting the patch +1. Each Bin has up to 4
  borders, but they are common (they share 2 with the next bin). Depending on
  the orientation could be equal to maxNBins. In other cases twice  maxNBins.
  Here we are very conservative*/

  out->maxNBorders = 1 + 2*maxNBins;

  /* LUT validity */
  if(maxDopplerBin < 2){
    /* the answer is infinity more or less */
    out->nFreqValid = f0Bin;
  }
  else {
    REAL8 num, den;
    num =  maxDopplerBin *maxDopplerBin;
    den = (maxDopplerBin -1)*(maxDopplerBin -1);
    out->nFreqValid =pixErr/(pixelFactor*vTotC)*sqrt(num/den -1.0);
  }

  /* max. angle (rad.) from the pole to consider a circle as a line in the projected plane */
  /* epsilon ~ 2/radius circle; radius ~h*h/(2b), epsilon = 4b/(h*h) */

  out->epsilon=8.0* linErr/(maxSide*deltaX*maxSide) ;

 /* -------------------------------------------   */


  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}





/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALHOUGHFillPatchGrid (LALStatus      *status,
                   HOUGHPatchGrid      *out,  /* */
                   HOUGHSizePar        *in1)
{ 

  REAL8   deltaX;
  REAL8   deltaY;
  UINT2   xSide;    /* number of pixels in the x direction (projected plane)*/
  UINT2   ySide;    /* number of pixels in the y direction */

  INT4    i;
  REAL8   xMin, xMax, x1;
  REAL8   yMin, yMax, myy1;
  REAL8   *xCoord;
  REAL8   *yCoord;
  /* --------------------------------------------- */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in1, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (out->xCoor, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (out->yCoor, status, LUTH_ENULL, LUTH_MSGENULL);

  xCoord = out->xCoor;
  yCoord = out->yCoor;

  xSide = out->xSide;
  ySide = out->ySide;

  /* Make sure there are physical pixels in that patch */
  ASSERT (xSide, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (ySide, status, LUTH_EVAL, LUTH_MSGEVAL);

  deltaX = out->deltaX = in1->deltaX;
  deltaY = out->deltaY = in1->deltaY;

  /* Make sure pixel sizes are positive definite */
  ASSERT (deltaX > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (deltaY > 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);

  out->f0 = in1->f0Bin * in1->deltaF;
  out->deltaF = in1->deltaF;

  /* -------------------------------------------   */
  /* Calculation of the patch limits with respect to the centers of the
     last pixels.  Note xSide and ySide are integers */
  xMax = out->xMax = deltaX*0.5*(xSide-1);
  yMax = out->yMax = deltaY*0.5*(ySide-1);

  xMin = out->xMin = -xMax;
  yMin = out->yMin = -yMax;

  /* -------------------------------------------   */

  /* Coordiantes of the pixel centers, in the projected plane */
  x1=xMin;
  for (i=0;i<xSide;++i){
    xCoord[i] = x1;
    x1+= deltaX;
  }

  myy1=yMin;
  for (i=0;i<ySide;++i){
    yCoord[i] = myy1;
    myy1+= deltaY;
  }

 /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}









