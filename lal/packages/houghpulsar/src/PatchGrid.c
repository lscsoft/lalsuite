/*-----------------------------------------------------------------------
 *
 * File Name: PatchGrid.c
 *
 * Authors: Sintes, A.M., Krishnan, B.
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes May 14, 2001
 *            Modified by Badri Krishnan Feb 2003
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

/************************************ <lalVerbatim file="PatchGridCV">
Author: Sintes, A. M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Module \texttt{PatchGrid.c}}
\label{ss:PatchGrid.c}
Function for tiling  the sky-patch (on the projected plane).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{PatchGridD}
\index{\verb&LALHOUGHPatchGrid()&}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

This  is a {\sc provisional} routine for tiling a  sky-pacth
on the projected plane.
The reason to call it  {\sc provisional}  is because
the size of the patch depends on the grid used in the 
demodulation stage. Neighbour sky-patches should not be separated
nor overlapping too much.
Here for setting the patch size, we  consider only $v_{epicycle}$, 
the frequency \verb@f0@ and  \verb@deltaF@ so that the \lq longitudinal'  size
of the patch is given by \verb@side == deltaF/f0 * c/v_epi@.
 By taking \verb@f0@ to be the maximun frequency considered in that step, 
 the patch-size is valid for a whole frequency range.\\
  
  
 Given input parameters,  the function \verb&LALHOUGHPatchGrid()& provides  
  patch information.
  
 The input \verb@*in1@ is a structure of type  \verb@HOUGHResolutionPar@ 
 containing
  some resolution parameters such as:
 \texttt{in1->f0} a frequency, \texttt{in1->deltaF} the frequency resolution, and
 \texttt{in1->minWidthRatio}  the ratio between the minimum  annulus width
for this search and the minimun  annulus width for  1 year integration time.
This value should be in the interval  [1.0, 25.0].
 
 The output structure \verb@*out@  of type \verb@HOUGHPatchGrid@ 
 stores patch grid information. The fields are:
\begin{description}
\item[\texttt{out->f0 }]  The frequency to construct grid 
\item[\texttt{out->deltaF }]  The frequency resolution: \texttt{df=1/TCOH}
\item[\texttt{out->minWidthRatio }]  Same as \texttt{in1->minWidthRatio}
\item[\texttt{out->deltaX }] Space discretization in the  x-direction
\verb@=deltaF * minWidthRatio/(f0 * VTOT * PIXELFACTORX)@.
\item[\texttt{out->xMin }]  Minimum  x value allowed, given as the  center of the first pixel.
\item[\texttt{out->xMax }] Maximun x value allowed, given as the  center of the last pixel.
\item[\texttt{out->xSide }] Real number of pixels in the x direction (in the
projected plane). It should be smaller or equal to \texttt{xSideMax} 
(\verb@xSide = VTOT * PIXELFACTORX / (VEPI * minWidthRatio)@ ).
\item[\texttt{out->xSideMax }]  Length of \texttt{xCoor}.
\item[\texttt{out->xCoor }] Coordinates of the pixel centers in the
x-direction.
\item[\texttt{out->deltaY }] Space resolution in the  y-direction
\verb@=deltaF * minWidthRatio/(f0 * VTOT * PIXELFACTORY)@.
\item[\texttt{out->yMin }]  Minimum  y value allowed, given as the  center of the first pixel.
\item[\texttt{out->yMax }]  Maximun y value allowed, given as the  center of the last pixel.
\item[\texttt{out->ySide }] Real number of pixels in the y-direction. It should
be smaller or  equal to \texttt{ySideMax} 
 (\verb@ySide = VTOT * PIXELFACTORY / (VEPI * minWidthRatio)@ ).
\item[\texttt{out->ySideMax }]  Length of \texttt{yCoor}.
\item[\texttt{out->yCoor }] Coordinates of the pixel centers in the
y-direction.
\end{description}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
%%\begin{verbatim}
%%LALZDestroyVector()
%%\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{PatchGridCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
</lalLaTeX> */


#include <lal/LUT.h>

NRCSID (PATCHGRIDC, "$Id$");

/* <lalVerbatim file="PatchGridD"> */
void LALHOUGHPatchGrid (LALStatus      *status,
                   HOUGHPatchGrid      *out,  /* */
                   HOUGHResolutionPar  *in1)  /* information */
{ /* </lalVerbatim> */

  /* vvvvvvvvvvvvvv */
  /* to be modified */
  /* ^^^^^^^^^^^^^^ */

  /* --------------------------------------------- */

  INT4    i;
  REAL8   f0; /* frequency to construct grid */
  REAL8   deltaF;  /* df=1/TCOH */
  REAL8   minWidthRatio;
  /*(min annuli width in this search)/(min annuli width in 1 year) 
    [1.0, 25.0]*/
  REAL8   deltaX;
  REAL8   xMin, xMax, x1;
  REAL8   patchSizeX, patchSizeY;  /* Size of sky patch in radians */
  UINT2   xSide;
  REAL8   *xCoord;

  REAL8   deltaY;
  REAL8   yMin, yMax, y1;
  UINT2   ySide;
  REAL8   *yCoord;
  /* --------------------------------------------- */

  INITSTATUS (status, "LALHOUGHPatchGrid", PATCHGRIDC);
  ATTATCHSTATUSPTR (status); 

  /*   Make sure the arguments are not NULL: */ 
  ASSERT (out, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (in1, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (out->xCoor, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (out->yCoor, status, LUTH_ENULL, LUTH_MSGENULL);
  /* Make sure the patch contains some pixels */
  ASSERT (out->xSideMax, status, LUTH_ESIZE, LUTH_MSGESIZE);
  ASSERT (out->ySideMax, status, LUTH_ESIZE, LUTH_MSGESIZE);
  
  ASSERT (in1->minWidthRatio >= 1.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (in1->minWidthRatio <= 25.0,status, LUTH_EVAL, LUTH_MSGEVAL);

  /* Make sure the user chose a sky patch smaller that 1/6 of the sky */ 
  ASSERT (in1->patchSizeX >= 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (in1->patchSizeX <= 1.447,status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (in1->patchSizeY >= 0.0, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (in1->patchSizeY <= 1.447,status, LUTH_EVAL, LUTH_MSGEVAL);
  /* -------------------------------------------   */
  /* The size of the patch depends on the grid used for the 
     demodulation part. Neigbour sky-patches should not be separated,
     nor overlapping too much.
     The patch size here considers only v_epicicle, f0 & deltaF:
                  side == deltaF/f0 * c/v_epi
     But the patch-size should be valid for a certain frequency range. 
     Note:  The above equation  will be used ONLY if the user did 
            not specify the sky-patch size      */
 /*  ------------------------------------------ */

  f0 = out->f0 = in1->f0;
  deltaF = out->deltaF = in1->deltaF;
  patchSizeX = out->patchSizeX = in1->patchSizeX;
  patchSizeY = out->patchSizeY = in1->patchSizeY;
  minWidthRatio = out->minWidthRatio = in1->minWidthRatio;

  deltaX =  out->deltaX =  deltaF * minWidthRatio/ (f0 * VTOT*PIXELFACTORX );
  deltaY =  out->deltaY =  deltaF * minWidthRatio/ (f0 * VTOT*PIXELFACTORY );


  if ( patchSizeX*patchSizeY ) {
      xSide = out->xSide = floor( patchSizeX/deltaX );
      ySide = out->ySide = floor( patchSizeY/deltaY );
  } else {
      xSide =  out->xSide =  VTOT * PIXELFACTORX / ( VEPI * minWidthRatio);
      ySide =  out->ySide =  VTOT * PIXELFACTORY / ( VEPI * minWidthRatio);
      patchSizeX = out->patchSizeX = xSide*deltaX;
      patchSizeY = out->patchSizeY = ySide*deltaY;
  }

  /* -------------------------------------------   */
  /* check for size mismatch */
  ASSERT (xSide <= out->xSideMax, status, LUTH_ESZMM, LUTH_MSGESZMM);
  ASSERT (ySide <= out->ySideMax, status, LUTH_ESZMM, LUTH_MSGESZMM);
  /* Make sure there are physical pixels in that patch */
  ASSERT (xSide, status, LUTH_EVAL, LUTH_MSGEVAL);
  ASSERT (ySide, status, LUTH_EVAL, LUTH_MSGEVAL);
  
  
  /* -------------------------------------------   */
  /* Calculation of the patch limits with respect to the centers of the 
     last pixels.  Note xSide and ySide are integers */
  xMax = out->xMax = deltaX*0.5*(xSide-1);
  yMax = out->yMax = deltaY*0.5*(ySide-1);
  
  xMin = out->xMin = -xMax;
  yMin = out->yMin = -yMax;  

  /* -------------------------------------------   */  
  /* Coordiantes of the pixel centers, in the projected plane */ 
 
  /*  xCoord = &(*out).xCoor[0];  or &(out->xCoor[0]) or out->xCoor */
  xCoord = out->xCoor;
  yCoord = out->yCoor;

  x1=xMin;
  for (i=0;i<xSide;++i){
    xCoord[i] = x1;
    x1+= deltaX;
  }
  
  y1=yMin;
  for (i=0;i<ySide;++i){
    yCoord[i] = y1;
    y1+= deltaY;
  }

 /* -------------------------------------------   */
  
  DETATCHSTATUSPTR (status);
  
  /* normal exit */
  RETURN (status);
}
