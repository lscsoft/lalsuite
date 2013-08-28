/*
*  Copyright (C) 2007 Jolien Creighton, Alicia Sintes Olives
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
 * File Name: Stereographic.c
 *
 * Authors: Sintes, A.M.,
 *
 *
 * History:   Created by Sintes May 15, 2001
 *            Modified...
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * Stereographic.c
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
 * \author Sintes, A. M.
 * \file
 * \ingroup LUT_h
 * \brief Routines to perform rotations on the celestial sphere and stereographic projection.
 *
 * \heading{Description}
 *
 * The function LALRotatePolarU() rotates the celestial sphere so that
 * a given point, in the rotated coordinates, corresponds to (\f$\alpha = 0\f$, \f$\delta = -\pi/2\f$).
 * The inputs are:
 * <tt>*par</tt> the reference point (e.g., the center of the sky-patch) of type
 * \c REAL8UnitPolarCoor and
 * <tt>*in</tt> the point on the celestial sphere we want to rotate. The output is
 * <tt>*out</tt> of type \c REAL8UnitPolarCoor containing the coordinates of the
 * point in the rotated reference frame.
 *
 * The function LALInvRotatePolarU() does the inverse rotation. Given the
 * reference point <tt>*par</tt> (e.g., the center of the sky-patch) of type
 * \c REAL8UnitPolarCoor  and a point <tt>*in</tt> in the rotated reference
 * frame, the output <tt>*out</tt> are the coordinates of the point is the
 * same reference system as <tt>*par</tt>. All inputs and output being
 * of type \c REAL8UnitPolarCoor.
 *
 * Given a point on the celestial sphere <tt>*in</tt> of type
 * \c REAL8UnitPolarCoor, the function LALStereoProjectPolar()
 * returns <tt>*out</tt>,
 * of type \c REAL8Polar2Coor, the stereographic projection of that point
 * in polar coordinates, with the particularity  that <tt>out->radius</tt>  can be positive
 * or negative. <tt>in->delta</tt>=\f$\pi/2\f$ is an invalid argument  and an error will
 * output.
 *
 * Given a point on the celestial sphere <tt>*in</tt> of type
 * \c REAL8UnitPolarCoor, the function LALStereoProjectCart()
 * returns <tt>*out</tt>, of type \c REAL8Cart2Coor, the stereographic projection of that point
 * in Cartesian coordinates. <tt>in->delta</tt>=\f$\pi/2\f$ is an invalid argument  and an error will
 * output.
 *
 * Given a point on the projected plane <tt>*in</tt> , the functions
 * LALStereoInvProjectPolar()  and LALStereoInvProjectCart()
 * provide the corresponding point on the sphere <tt>*out</tt> (corresponding to the inverse
 * stereographic  projection) of type \c REAL8UnitPolarCoor.
 *
 */

#include <lal/LUT.h>

/*
 * The functions that make up the guts of this module
 */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALRotatePolarU(LALStatus            *status,
		     REAL8UnitPolarCoor   *out,
		     REAL8UnitPolarCoor   *in,
		     REAL8UnitPolarCoor   *par)
{

  REAL8 Xx, Xy, Xz;
  REAL8 Yx, Yy, Yz;
  REAL8 Zx, Zy, Zz;
  REAL8 alphaN, deltaN;
  REAL8 alphaIn, deltaIn;

  REAL8 VIx, VIy, VIz;
  REAL8 Vx, Vy, Vz;

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (par, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in,  status, LUTH_ENULL, LUTH_MSGENULL);
  /* -------------------------------------------   */

  alphaN = par->alpha;
  deltaN = par->delta;

  alphaIn = in->alpha;
  deltaIn = in->delta;

  /* -------------------------------------------   */
  /* Initial vector */
  VIx = cos(deltaIn) * cos(alphaIn) ;
  VIy = cos(deltaIn) * sin(alphaIn);
  VIz = sin(deltaIn);

  /* -------------------------------------------   */
  /* Rotation matrix */
  Xx =  -sin(deltaN) * cos(alphaN) ;
  Xy =  -sin(deltaN) * sin(alphaN);
  Xz =   cos(deltaN);
  Yx =  -sin(alphaN);
  Yy =   cos(alphaN);
  Yz =    0.0;
  Zx =  -cos(deltaN) * cos(alphaN);
  Zy =  -cos(deltaN) * sin(alphaN);
  Zz =  -sin(deltaN);

  /* -------------------------------------------   */
  /* Final vector */

  Vx = Xx*VIx + Xy*VIy + Xz*VIz ;
  Vy = Yx*VIx + Yy*VIy + Yz*VIz ;
  Vz = Zx*VIx + Zy*VIy + Zz*VIz ;

  /* -------------------------------------------   */

  /* output */

  if( Vx || Vy ){
    out->alpha = atan2(Vy,Vx);
  } else {
    out->alpha =0.0;
  }
  out->delta = asin(Vz);

  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}




/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALInvRotatePolarU(LALStatus            *status,
		     REAL8UnitPolarCoor   *out,
		     REAL8UnitPolarCoor   *in,
		     REAL8UnitPolarCoor   *par)
{

  REAL8 Xx, Xy, Xz;
  REAL8 Yx, Yy, Yz;
  REAL8 Zx, Zy, Zz;
  REAL8 alphaN, deltaN;
  REAL8 alphaIn, deltaIn;

  REAL8 VIx, VIy, VIz;
  REAL8 Vx, Vy, Vz;

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (par, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in,  status, LUTH_ENULL, LUTH_MSGENULL);
  /* -------------------------------------------   */

  alphaN = par->alpha;
  deltaN = par->delta;

  alphaIn = in->alpha;
  deltaIn = in->delta;

  /* -------------------------------------------   */
  /* Initial vector */
  VIx = cos(deltaIn) * cos(alphaIn) ;
  VIy = cos(deltaIn) * sin(alphaIn);
  VIz = sin(deltaIn);

  /* -------------------------------------------   */
  /* Inverse Rotation matrix */
  Xx =  - cos(alphaN) * sin(deltaN) ;
  Xy =  - sin(alphaN);
  Xz =  - cos(alphaN) * cos(deltaN);
  Yx =  - sin(alphaN) * sin(deltaN);
  Yy =    cos(alphaN);
  Yz =  - sin(alphaN) * cos(deltaN) ;
  Zx =     cos(deltaN);
  Zy =     0.0;
  Zz =    -sin(deltaN);

  /* -------------------------------------------   */
  /* Final vector */

  Vx = Xx*VIx + Xy*VIy + Xz*VIz ;
  Vy = Yx*VIx + Yy*VIy + Yz*VIz ;
  Vz = Zx*VIx + Zy*VIy + Zz*VIz ;

  /* -------------------------------------------   */

  /* output */

  if( Vx || Vy ){
    out->alpha = atan2(Vy,Vx);
  } else {
    out->alpha =0.0;
  }
  out->delta = asin(Vz);

  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}

/* -------------------------------------------   */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALStereoProjectPolar(LALStatus           *status,
			   REAL8Polar2Coor     *out,
			   REAL8UnitPolarCoor  *in)
{

  REAL8   mygamma;
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in,  status, LUTH_ENULL, LUTH_MSGENULL);
  /*   Make sure delta is not pi/2: */
  ASSERT (in->delta - LAL_PI_2 , status, LUTH_EVAL, LUTH_MSGEVAL);
  /* -------------------------------------------   */

  out->alpha = in->alpha;
  mygamma = LAL_PI_4 + 0.5*(in->delta);
  out->radius = 2.0 * tan(mygamma); /* positive or negative ! */
  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}




/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALStereoProjectCart(LALStatus           *status,
			  REAL8Cart2Coor      *out,
			  REAL8UnitPolarCoor  *in )
{

  REAL8   mygamma;
  REAL8   alpha, radius;
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in,  status, LUTH_ENULL, LUTH_MSGENULL);
  /*   Make sure delta is not  pi/2: */
  ASSERT (in->delta - LAL_PI_2 , status, LUTH_EVAL, LUTH_MSGEVAL);
  /* -------------------------------------------   */

  alpha = in->alpha;
  mygamma = LAL_PI_4 + 0.5*(in->delta);
  radius = 2.0 * tan(mygamma); /* positive or negative */

  out->x = radius * cos(alpha);
  out->y = radius * sin(alpha);
  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALStereoInvProjectPolar(LALStatus        *status,
			   REAL8UnitPolarCoor  *out,
			   REAL8Polar2Coor     *in)
{


  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in,  status, LUTH_ENULL, LUTH_MSGENULL);
  /* -------------------------------------------   */

  out->alpha = in->alpha;
  out->delta = 2.0* atan(0.5*(in->radius) ) - LAL_PI_2;

  /*  Note: since I have not ask for a positive radius input, */
  /*  delta in principle is not cofined to (-pi/2, pi/2)      */
  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALStereoInvProjectCart(LALStatus           *status,
			     REAL8UnitPolarCoor  *out,
			     REAL8Cart2Coor      *in)
{

  REAL8 x,y,radius;
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in,  status, LUTH_ENULL, LUTH_MSGENULL);
  /* -------------------------------------------   */

  x = in->x;
  y = in->y;
  if ( x || y ){
    radius = sqrt(x*x + y*y);

    out->alpha = atan2(y,x);
    out->delta = 2.0* atan(0.5*(radius) ) - LAL_PI_2;
    /*  Note: since I have not ask for a positive radius input, */
    /*  delta in principle is not cofined to (-pi/2, pi/2)      */

  } else { /* x=y =0 */
    out->alpha = 0.0;
    out->delta = - LAL_PI_2;
  }

  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}



