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
 * File Name: ParamPLUT.c
 *
 * Authors: Sintes, A.M., Krishnan, B.,
 *
 *
 * History:   Created by Sintes May 15, 2001
 *            Modified by Badri Krishnan Feb 2003
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * ParamPLUT.c
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
 * \author Sintes, A. M., Krishnan, B.
 * \file
 * \ingroup LUT_h
 * \brief Function that calculates the parameters needed for generating the look-up-table.
 *
 * \heading{Description}
 * This routine calculates the parameters needed for generating the look-up-table.
 * It is valid for all cases in which the Hough transform
 * master equation is of the form:
 * \f$f(t)\f$-\c f0 = \f$\vec\xi \cdot (\hat n-\hat N)\f$, or
 * equivalently,
 * \f$\cos(\phi)\f$ = (\f$f(t)-\f$\c f0 + \f$\vec\xi \cdot\hat N\f$)/\f$\vert\vec\xi\vert\f$.
 * \f$\vec\xi\f$, hereafter \c xi, is calculated according to the demodulation procedure used in a
 * first stage.\\
 *
 * \heading{Uses}
 * \code
 * LALRotatePolarU()
 * \endcode
 *
 */

#include <lal/LUT.h>

void LALHOUGHCalcParamPLUT (LALStatus    *status,
                   HOUGHParamPLUT    *out, /* parameters needed build LUT*/
                   HOUGHSizePar      *size,
                   HOUGHDemodPar     *par)  /* demodulation parameters */
{

  /* --------------------------------------------- */

  REAL8   f0;  /* frequency corresponding to f0Bin */
  INT8    f0Bin;
  REAL8   deltaF;  /*  df=1/TCOH  */
  REAL8   delta;
  REAL8   vFactor, xFactor;
  REAL8   xiX, xiY, xiZ;
  REAL8   modXi,invModXi;
  REAL8UnitPolarCoor   xiInit;
  UINT4   spinOrder, i;
  REAL8   *spinF;
  REAL8   timeDiff;    /*  T(t)-T(t0) */
  REAL8   timeDiffProd;
  /* --------------------------------------------- */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (par, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (size, status, LUTH_ENULL, LUTH_MSGENULL);

  /*   Make sure f0Bin  is not zero: */
  f0Bin = size->f0Bin;
  ASSERT (f0Bin, status, LUTH_EFREQ, LUTH_MSGEFREQ);

  out->f0Bin = f0Bin;
  deltaF = out->deltaF = size->deltaF;

  f0 = f0Bin * deltaF;

  out->epsilon = size->epsilon;
  out->nFreqValid = size->nFreqValid;
  /* -------------------------------------------   */


  /* -------------------------------------------   */
  /* *********** xi calculation *****************  */

  vFactor = f0;
  xFactor = 0.0;

  spinOrder = par->spin.length;

  if(spinOrder){
    ASSERT (par->spin.data , status, LUTH_ENULL, LUTH_MSGENULL);
    timeDiff = par->timeDiff;
    timeDiffProd = 1.0;
    spinF = par->spin.data;

    for (i=0; i<spinOrder; ++i ){
      xFactor += spinF[i] * timeDiffProd * (i+1.0);
      timeDiffProd *= timeDiff;
      vFactor += spinF[i] * timeDiffProd;
    }
  }

  xiX = vFactor * (par->veloC.x) + xFactor * (par->positC.x);
  xiY = vFactor * (par->veloC.y) + xFactor * (par->positC.y);
  xiZ = vFactor * (par->veloC.z) + xFactor * (par->positC.z);

  /* -------------------------------------------   */
  /* ***** convert xi into Polar coordinates ***** */

  modXi = sqrt(xiX*xiX + xiY*xiY + xiZ*xiZ);
  /* for testing we used:  modXi = F0* 1.06e-4; */
  invModXi = 1./modXi;

  xiInit.delta = asin( xiZ*invModXi);
  /* the arc sine is in the interval [-pi/2,pi/2] */

  if( xiX || xiY ){
    xiInit.alpha = atan2(xiY, xiX);
  }else{
    xiInit.alpha = 0.0;
  }

  /*   if( (xiX == 0.0 ) && (xiY == 0.0 ) ){ */
  /*     xiInit.alpha = 0.0; */
  /*   }else{  xiInit.alpha = atan2(xiY, xiX);  } */

  /* -------------------------------------------   */
  /* **** Rotate Patch, so that its center becomes */
  /* **** the south pole {x,y,z} = {0,0,-1}  ***** */
  /*  Calculate xi in the new coordinate system.   */

  TRY(LALRotatePolarU(status->statusPtr,
		      &(*out).xi ,&xiInit, &(*par).skyPatch), status);

  /* -------------------------------------------   */
  delta = out->xi.delta;

  out->cosDelta = deltaF*invModXi;
  out->cosPhiMax0 = deltaF*0.5*invModXi -sin(delta);
  out->cosPhiMin0 = (out->cosPhiMax0) -(out->cosDelta);
  out->offset = 0;

  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}




