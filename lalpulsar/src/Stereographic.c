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
 * Revision: $Id$
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

/************************************ <lalVerbatim file="StereographicCV">
Author: Sintes, A. M.
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Module \texttt{Stereographic.c}}
\label{ss:Stereographic.c}
Routines to perform rotations on the celestial sphere and stereographic
projection.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{StereographicD}
\index{\verb&LALRotatePolarU()&}
\index{\verb&LALInvRotatePolarU()&}
\index{\verb&LALStereoProjectPolar()&}
\index{\verb&LALStereoProjectCart()&}
\index{\verb&LALStereoInvProjectPolar()&}
\index{\verb&LALStereoInvProjectCart()&}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

${}$

 \indent The function \verb&LALRotatePolarU()& rotates the celestial sphere so that
 a given point,
 in the rotated coordinates, corresponds to ($\alpha = 0$, $\delta = -\pi/2$).
The inputs are:
\verb@*par@ the reference point (e.g., the center of the sky-patch) of type
 \verb@REAL8UnitPolarCoor@ and
\verb@*in@ the point on the celestial sphere we want to rotate. The output is
\verb@*out@ of type \verb@REAL8UnitPolarCoor@ containing the coordinates of the
point in the rotated reference frame.\\

 The function \verb&LALInvRotatePolarU()& does the inverse rotation. Given the
reference point \verb@*par@ (e.g., the center of the sky-patch) of type
 \verb@REAL8UnitPolarCoor@  and a point \verb@*in@ in the rotated reference
 frame, the output \verb@*out@ are the coordinates of the point is the
 same reference system as \verb@*par@. All inputs and output being
  of type \verb@REAL8UnitPolarCoor@. \\

Given a point on the celestial sphere \verb@*in@ of type
\verb@REAL8UnitPolarCoor@, the function \verb&LALStereoProjectPolar()&
returns \verb@*out@,
of type \verb@REAL8Polar2Coor@, the stereographic projection of that point
in polar coordinates, with the particularity  that \verb@out->radius@  can be positive
or negative. \verb@in->delta@=$\pi/2$ is an invalid argument  and an error will
output. \\

Given a point on the celestial sphere \verb@*in@ of type
\verb@REAL8UnitPolarCoor@, the function \verb@LALStereoProjectCart()@
returns \verb@*out@, of type \verb@REAL8Cart2Coor@, the stereographic projection of that point
in Cartesian coordinates. \verb@in->delta@=$\pi/2$ is an invalid argument  and an error will
output. \\

Given a point on the projected plane \verb@*in@ , the functions
\verb&LALStereoInvProjectPolar()&  and \verb&LALStereoInvProjectCart()&
provide the corresponding point on the sphere \verb@*out@ (corresponding to the inverse
stereographic  projection) of type \verb@REAL8UnitPolarCoor@.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
%%\begin{verbatim}
%%LALZDestroyVector()
%%\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{StereographicCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
</lalLaTeX> */



#include <lal/LUT.h>

NRCSID (STEREOGRAPHICC, "$Id$");


/*
 * The functions that make up the guts of this module
 */



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*  ***************************** <lalVerbatim file="StereographicD"> */
void LALRotatePolarU(LALStatus            *status,
		     REAL8UnitPolarCoor   *out,
		     REAL8UnitPolarCoor   *in,
		     REAL8UnitPolarCoor   *par)
{ /*  ************************************************ </lalVerbatim> */

  REAL8 Xx, Xy, Xz;
  REAL8 Yx, Yy, Yz;
  REAL8 Zx, Zy, Zz;
  REAL8 alphaN, deltaN;
  REAL8 alphaIn, deltaIn;

  REAL8 VIx, VIy, VIz;
  REAL8 Vx, Vy, Vz;

  /* --------------------------------------------- */
  INITSTATUS (status, " LALRotatePolarU", STEREOGRAPHICC);
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
/*  ***************************** <lalVerbatim file="StereographicD"> */
void LALInvRotatePolarU(LALStatus            *status,
		     REAL8UnitPolarCoor   *out,
		     REAL8UnitPolarCoor   *in,
		     REAL8UnitPolarCoor   *par)
{ /*  ************************************************ </lalVerbatim> */

  REAL8 Xx, Xy, Xz;
  REAL8 Yx, Yy, Yz;
  REAL8 Zx, Zy, Zz;
  REAL8 alphaN, deltaN;
  REAL8 alphaIn, deltaIn;

  REAL8 VIx, VIy, VIz;
  REAL8 Vx, Vy, Vz;

  /* --------------------------------------------- */
  INITSTATUS (status, " LALInvRotatePolarU", STEREOGRAPHICC);
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
/*  ***************************** <lalVerbatim file="StereographicD"> */
void LALStereoProjectPolar(LALStatus           *status,
			   REAL8Polar2Coor     *out,
			   REAL8UnitPolarCoor  *in)
{ /*  ************************************************ </lalVerbatim> */

  REAL8   mygamma;
  /* --------------------------------------------- */
  INITSTATUS (status, "LALStereoProjectPolar", STEREOGRAPHICC);
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
/*  ***************************** <lalVerbatim file="StereographicD"> */
void LALStereoProjectCart(LALStatus           *status,
			  REAL8Cart2Coor      *out,
			  REAL8UnitPolarCoor  *in )
{ /*  ************************************************ </lalVerbatim> */

  REAL8   mygamma;
  REAL8   alpha, radius;
  /* --------------------------------------------- */
  INITSTATUS (status, "LALStereoProjectCart", STEREOGRAPHICC);
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
/*  ***************************** <lalVerbatim file="StereographicD"> */
void LALStereoInvProjectPolar(LALStatus        *status,
			   REAL8UnitPolarCoor  *out,
			   REAL8Polar2Coor     *in)
{ /*  ************************************************ </lalVerbatim> */


  /* --------------------------------------------- */
  INITSTATUS (status, "LALStereoInvProjectPolar", STEREOGRAPHICC);
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
/*  ***************************** <lalVerbatim file="StereographicD"> */
void LALStereoInvProjectCart(LALStatus           *status,
			     REAL8UnitPolarCoor  *out,
			     REAL8Cart2Coor      *in)
{ /*  ************************************************ </lalVerbatim> */

  REAL8 x,y,radius;
  /* --------------------------------------------- */
  INITSTATUS (status, "LALStereoInvProjectCart", STEREOGRAPHICC);
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



