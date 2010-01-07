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

/**
 * \author Alicia Sintes, Badri Krishnan
 * \ingroup pulsarHough
 * \file ParamPLUT.c
 * \brief Routines for creating Look-Up-Tables for demodulated data
 */

/*-----------------------------------------------------------------------
 *
 * File Name: ParamPLUT.c
 *
 * Authors: Sintes, A.M., Krishnan, B.,
 *
 * Revision: $Id$
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

/************************************ <lalVerbatim file="ParamPLUTCV">
Author: Sintes, A. M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Module \texttt{ParamPLUT.c}}
\label{ss:ParamPLUT.c}
Function that calculates the parameters needed for generating the look-up-table.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ParamPLUTD}
\index{\verb&LALHOUGHParamPLUT()&}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}
This routine calculates the parameters needed for generating the look-up-table.
It is valid for all cases in which the Hough transform
 master equation is of the form:
$f(t)$-\verb@f0@ = $\vec\xi \cdot (\hat n-\hat N)$, or
equivalently,
$\cos(\phi)$ = ($f(t)-$\verb@f0@ + $\vec\xi \cdot\hat N$)/$\vert\vec\xi\vert$.
$\vec\xi$, hereafter \verb@xi@, is calculated according to the demodulation procedure used in a
 first stage.\\



The inputs are:
\begin{description}
\item[\texttt{INT8   f0Bin}]: The frequency bin to construct the {\sc lut}.
\item[\texttt{HOUGHDemodPar  *par}]:  The demodulation parameters:
\begin{description}
\item[\texttt{par->deltaF }]: Frequency resolution: \texttt{df=1/TCOH}.
\item[\texttt{par->skyPatch }]:  $N_{center}$ (alpha, delta):
position of the center of the patch.
\item[\texttt{par->patchSizeX }]: Size of sky patch along x-axis measured in radians.
\item[\texttt{par->patchSizeY }]: Size of sky patch along y-axis measured in radians.
\item[\texttt{par->veloC }]:  $v(t)/c$ (x,y,z): relative detector
velocity.
\item[\texttt{par->positC }]: $(\vec x(t)-\vec x(\hat t_0))/c$ (x,y,z). Position
of the detector.
\item[\texttt{par->timeDiff }]:  $T_{\hat N}(t)-T_{\hat N}(\hat t_0)$: Time difference.
\item[\texttt{par->spin }]:  \texttt{length}: Maximum order of
spin-down parameter.
 \texttt{*data}: Pointer to spin-down parameter set $F_k$.
\end{description}
\end{description}

The output \verb@*out@ of type \verb@HOUGHParamPLUT@ contains
all the parameters needed to build the look-up-table for constructing
the partial Hough maps. Those are:
\begin{description}
\item[\texttt{out->f0Bin }]: Frequency  bin for which it has been constructed.
\item[\texttt{out->deltaF }]: Frequency resolution: \texttt{df=1/TCOH}.
\item[\texttt{out->xi }]: Center of the circle on the celestial
sphere, xi(alpha,delta) in the rotated coordinates.
\item[\texttt{out->cosDelta }]: $\Delta \cos(\phi)$ corresponding to
one annulus: \verb@deltaF/|xi|@.
\item[\texttt{out->cosPhiMax0 }]: $\max(\cos(\phi))$ of the
\texttt{f0Bin}  : \verb@(xi*N +deltaF/2)/|xi|@.
\item[\texttt{out->cosPhiMin0 }]:  $\min(\cos(\phi))$ of the
\texttt{f0Bin}  : \verb@cosPhiMax0-cosDelta@.
\item[\texttt{out->epsilon }]: Maximum angle (distance in radians) from the pole
to consider  a circle as a line in the projected plane:
\verb@8.* LINERR * f0Bin* VEPI * VEPI / VTOT@. For explanations see Sintes'
notes.
\end{description}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALRotatePolarU()
\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{ParamPLUTCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
</lalLaTeX> */


#include <lal/LUT.h>

NRCSID (PARAMPLUTC, "$Id$");

/* <lalVerbatim file="ParamPLUTD"> */
void LALHOUGHParamPLUT (LALStatus    *status,
                   HOUGHParamPLUT    *out, /* parameters needed build LUT*/
                   HOUGHSizePar      *size,
                   HOUGHDemodPar     *par)  /* demodulation parameters */
{ /* </lalVerbatim> */

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

  INITSTATUS (status, "LALHOUGHParamPLUT", PARAMPLUTC);
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




