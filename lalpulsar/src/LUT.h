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

/** \file
 * \addtogroup moduleHoughPulsar
 * \date $Date$
 * \brief Provides structures and function prototypes required for the construction of look up tables
   that are the core for building the Hough maps.
 *  \author Sintes, A.M, Papa, M.A. and Krishnan, B.
 * File Name: LUT.h
 *
 *
 *
 * $Id$
 *
 * History:   Created by Sintes May 11, 2001
 *            Modified by Badri Krishnan Feb 2003
 *

\par Description

 Our goal is the construction of Hough maps. In order to produce them
 efficiently, the present implemetation makes use of {\sc lut}s.
 Here we provide the necessary routines for their construction and use.\\

 In principle, the subroutines provided are valid for
  {\sc any} Hough master equation of the form:
   \f$ \nu-F_0  =\vec\xi (t) \cdot (\hat n -\hat N )\, ,\f$
 where \f$\nu\f$ is the measured frequency of the signal at time \f$t\f$,
  \f$F_0\f$ intrinsic frequency of the signal at that time, \f$\hat n\f$ location of the
  souce in the sky, \f$\hat N\f$ the center of the sky patch used in the
  demodulation procedure,
   and
  \f$\vec\xi (t)\f$ any vector.

  The  form of this vector \f$\vec\xi (t)\f$
  depends on the demodulation procedure
   used in the previous  step.  In our case this corresponds to
\f$\vec\xi (t) = \left( F_0+
\sum_k F_k \left[ \Delta T  \right]^k\right) \frac{\vec v(t)}{c}
 + \left( \sum_k k F_k \left[ \Delta  T \right]^{k-1}\right)
\frac {\vec x(t)- \vec x(\hat t_0)}{c}\, ,\f$
and
   \f$F_0 \equiv  f_0 + \sum_k \Delta f_k
 \left[ \Delta T \right]^k \, , \f$
 where
 \f$\vec v(t)\f$ is the velocity of the detector, \f$\vec x(t)\f$ is the detector
 position,
 \f$ T_{\hat N}(t)\f$ is the  time  at
 the solar system barycenter (for a given sky
  location \f$\hat N\f$),
   \f$\Delta T \equiv T_{\hat N}(t)-T_{\hat N}(\hat t_0)\f$,
  \f$\Delta f_k = f_k-F_k\f$ the  residual spin-down parameter,
\f$F_k\f$ the  spin-down parameter used in the demodulation, and \f$f_0\f$, \f$f_k\f$ the
intrinsic  frequency and spin-down parameters of the source at time \f$\hat t_0\f$.\\

 Looking at the generic Hough master equation, one realizes that
 for a fixed  time, a given value of \f$F_0\f$, and a measured frequency \f$\nu\f$
  (from a selected peak), the source could be located anywhere on a circle (whose
 center points in the same direction of \f$\vec\xi (t)\f$ and is characterized by
 \f$\phi\f$, the angle between \f$\hat n\f$ and \f$\vec\xi\f$).  Since the Hough transform is performed on a
  set of spectra with discrete frequencies, a peak on the spectrum appearing at
  \f$\nu\f$  could correspond to any source with a demodulated
   frequency in a certain interval. As a consequence, the location of the sources
   compatible with \f$F_0\f$ and \f$\nu\f$ is not a circle but an annulus with a certain
   width.\\

   Our purpose is to map these annuli on a discrete space. An estimation of the
   average thickness of the annuli  tells us that the vast majority of
    annuli will be very thin, and therefore our algorithm should not be
    optimized for drawing
   thick annuli but for thin ones. Also, the mapping implementation should be one with a uniform
    probability distribution in order to avoid discretization errors.
    In order to remove border effects, we use a biunivocal mapping, which
    requires
    that a pixel in a partial Hough map can belong only to one annulus,
    just touched by one peak of the spectrum. The criteria for the biunivocal mapping
    is that if and only if the center of the pixel is inside the annulus, then the pixel
   will be enhanced.\\

  In order to simplify  (reduce the computational cost of) this task we
  construct look up tables ({\sc lut}) where the  borders of these annuli are
  marked for any possible \f$\nu
   -F_0\f$ value. Since we work on a discrete space these {\sc lut} are
   valid for many \f$F_0\f$ values.\\

\par Uses
\code
LALHOUGHComputeSizePar ()
LALHOUGHComputeNDSizePar ()
LALHOUGHFillPatchGrid ()
LALHOUGHParamPLUT ()
LALNDHOUGHParamPLUT ()
LALRotatePolarU()
LALInvRotatePolarU()
LALStereoProjectPolar()
LALStereoProjectCart()
LALStereoInvProjectPolar()
LALStereoInvProjectCart()
LALHOUGHConstructPLUT()
\endcode

 */


/************************************ <lalVerbatim file="LUTHV">
Author: Sintes, A.M.,  Papa, M.A.,  Krishnan, B.
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Header \texttt{LUT.h}}
\label{s:LUT.h}

Provides structures and functions required for the construction of look up tables
{\sc (lut)} that are the core for building the Hough maps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Synopsis}

\begin{verbatim}
#include <lal/LUT.h>
\end{verbatim}



 Our goal is the construction of Hough maps. In order to produce them
 efficiently, the present implemetation makes use of {\sc lut}s.
 Here we provide the necessary routines for their construction and use.\\

 In principle, the subroutines provided are valid for
  {\sc any} Hough master equation of the form:
   $$ \nu-F_0  =\vec\xi (t) \cdot (\hat n -\hat N )\, ,$$
 where $\nu$ is the measured frequency of the signal at time $t$,
  $F_0$ intrinsic frequency of the signal at that time, $\hat n$ location of the
  souce in the sky, $\hat N$ the center of the sky patch used in the
  demodulation procedure,
   and
  $\vec\xi (t)$ any vector.

  The  form of this vector $\vec\xi (t)$
  depends on the demodulation procedure
   used in the previous  step.  In our case this corresponds to
$$\vec\xi (t) = \left( F_0+
\sum_k F_k \left[ \Delta T  \right]^k\right) \frac{\vec v(t)}{c}
 + \left( \sum_k k F_k \left[ \Delta  T \right]^{k-1}\right)
\frac {\vec x(t)- \vec x(\hat t_0)}{c}\, ,$$
and
   $$F_0 \equiv  f_0 + \sum_k \Delta f_k
 \left[ \Delta T \right]^k \, , $$
 where
 $\vec v(t)$ is the velocity of the detector, $\vec x(t)$ is the detector
 position,
 $ T_{\hat N}(t)$ is the  time  at
 the solar system barycenter (for a given sky
  location $\hat N$),
   $\Delta T \equiv T_{\hat N}(t)-T_{\hat N}(\hat t_0)$,
  $\Delta f_k = f_k-F_k$ the  residual spin-down parameter,
$F_k$ the  spin-down parameter used in the demodulation, and $f_0$, $f_k$ the
intrinsic  frequency and spin-down parameters of the source at time $\hat t_0$.\\

 Looking at the generic Hough master equation, one realizes that
 for a fixed  time, a given value of $F_0$, and a measured frequency $\nu$
  (from a selected peak), the source could be located anywhere on a circle (whose
 center points in the same direction of $\vec\xi (t)$ and is characterized by
 $\phi$, the angle between $\hat n$ and $\vec\xi$).  Since the Hough transform is performed on a
  set of spectra with discrete frequencies, a peak on the spectrum appearing at
  $\nu$  could correspond to any source with a demodulated
   frequency in a certain interval. As a consequence, the location of the sources
   compatible with $F_0$ and $\nu$ is not a circle but an annulus with a certain
   width.\\

   Our purpose is to map these annuli on a discrete space. An estimation of the
   average thickness of the annuli  tells us that the vast majority of
    annuli will be very thin, and therefore our algorithm should not be
    optimized for drawing
   thick annuli but for thin ones. Also, the mapping implementation should be one with a uniform
    probability distribution in order to avoid discretization errors.
    In order to remove border effects, we use a biunivocal mapping, which
    requires
    that a pixel in a partial Hough map can belong only to one annulus,
    just touched by one peak of the spectrum. The criteria for the biunivocal mapping
    is that if and only if the center of the pixel is inside the annulus, then the pixel
   will be enhanced.\\

  In order to simplify  (reduce the computational cost of) this task we
  construct look up tables ({\sc lut}) where the  borders of these annuli are
  marked for any possible $\nu
   -F_0$ value. Since we work on a discrete space these {\sc lut} are
   valid for many $F_0$ values.\\

%\begin{figure}
%\centerline{\vbox{
%\psfig{figure=stereo.ps,height=5.0cm,width=6.0cm}
%}}
%\vspace*{10pt}
%\caption[]{Stereographic projection.
%}
%\vspace*{10pt}
%\end{figure}
%

\begin{figure}
\noindent\includegraphics{LUTstereo}
%\noindent\includegraphics{width=6cm,angle=0}{LUTstereo}
\caption[]{Stereographic projection. [Credit to: D.Hilbert, S.Cohn-Vossen, P.
Nemenyi, {\it \lq\lq Geometry and Imagination"}, Chelsea Publishing Company, New
York 1952.]}
\end{figure}


%\noindent\includegraphics{LUTstereo}

 At this point we have already chosen a sky tiling to produce the Hough map  efficiently.
 It consists of changing coordinates so that the center of the   patch is located at
  $(0,-\pi/2)$ in ($\alpha-\delta$) (or in any other coordinate system), then we make
  use
  of the stereographic projection and we take horizontal and vertical lines on
  the projected plane at constant space separation.
  This projection is advantageous because it avoids distortions, i.e. the
pixel size is almost constant independently of the sky location, and makes the
algorithm simpler. The stereographic projection has the property to map circles
on the celestial sphere to circles on the projected plane.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Error conditions}
\vspace{0.1in}
\input{LUTHErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Constant declarations}
The initial settings at compile time are the following:\\


\texttt{LINERR     0.001}: Maximum \lq\lq error'',
in terms  of the width of the thinnest annulus, which allows one to
represent/approximate a circle by  a line. \\

\texttt{PIXERR     0.1}: Maximum \lq\lq error'', in terms of
the width of the thinnest annulus,
 which allows one to consider two border equivalents.
It is only relevant for determining the {\sc lut} frequency range validity.\\


\texttt{PIXELFACTORX  2}: Width of the thinnest annulus in terms of pixels
 (x-direction).\\

\texttt{PIXELFACTORY  2}: Width of the thinnest annulus in terms of pixels
(y-direction).\\

\texttt{VEPI 1.0e-06}: Earth epicycle velocity divided by the light velocity
$c$. TO BE CHANGED DEPENDING ON THE DETECTOR. \\

\texttt{VTOT 1.0e-04}: Earth total velocity divided by $c$. TO BE CHANGED
DEPENDING ON THE DETECTOR. \\

\texttt{SIDEX   (100* PIXELFACTORX)}, \\

\texttt{SIDEY   (100* PIXELFACTORY)}:  The maximun number of pixels in a
given patch is  \texttt{SIDEX}$\times$\texttt{SIDEY}. These numbers should be
integers.
They  come from assuming that: the \lq longitudinal' patch size is $dF*10^6/f_0$,
where $10^6 = c/v_{epicycle}$, the thinnest possible annulus is $dF*10^4/f_0$,
where $10^4 = c/v_{tot}$. Therefore, the ratio (longitudinal patch size)/(width
 thinnest annuli) = 100
 (100 TO BE CHANGED DEPENDING ON THE DETECTOR).\\


{\verb@MAX_N_BINS   150@}: Maximun number of frequency bins that can affect a
patch.\\


{\verb@MAX_N_BORDERS   208@}: Maximun number of borders in a patch.\\

Note several of these constants will be removed in a next release. They are
used only for tiling the sky-patch or as default parameters in the test codes.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Structures and type definitions}

\begin{verbatim}
typedef UCHAR COORType
\end{verbatim}
To be changed to {\texttt{INT2 COORType}} if the number of pixels in the
x-direction exceeds 255.

\begin{verbatim}
struct HOUGHBorder
\end{verbatim}
\index{\texttt{HOUGHBorder}}

\noindent This structure stores the border of a circle clipped on the projected plane.
 The fields are:

\begin{description}
\item[\texttt{INT4  yUpper}]  Upper y pixel index affected by this border.
\item[\texttt{INT4  yLower}]  Lower y pixel index affected by this border.
 \texttt{yUpper}$<$\texttt{yLower} or \texttt{yUpper}$<0$ are possible.
\item[\texttt{INT4  yCenter}] y pixel value of the center of the circle.
\item[\texttt{UINT2 ySide}] Length of xPixel.
\item[\texttt{COORType  *xPixel}]  x pixel index of the border to be
marked.
\end{description}



\begin{verbatim}
struct HOUGHBin2Border
\end{verbatim}
\index{\texttt{HOUGHBin2Border}}

\noindent This structure stores the border indexes corresponding to one frequency
bin plus the corrections to be added to the first column of the patch.
The fields are:

\begin{description}
\item[\texttt{INT2   leftB1}]  Border index to be used ({\it start-border \lq $+1$'}).
\item[\texttt{INT2   rightB1}] Border index to be used ({\it stop-border \lq $-1$'}).
\item[\texttt{INT2   leftB2}]  Border index  to be used ({\it start-border \lq $+1$'}).
\item[\texttt{INT2   rightB2}]  Border index  to be used ({\it stop-border \lq $-1$'}).
\item[\texttt{INT2   piece1max}]:  Interval limits of the (first piece) correction
 to the first column.
\item[\texttt{INT2   piece1min}]:  If \texttt{piece1min > piece1max} no
corrections should be added.
\item[\texttt{INT2   piece2max}]:  Interval limits of the (second piece)  correction
 to the first column.
\item[\texttt{INT2   piece2min}]:  If \texttt{piece2min > piece2max} no
corrections should be added.
\end{description}

\begin{verbatim}
struct HOUGHptfLUT
\end{verbatim}
\index{\texttt{HOUGHptfLUT}}

\noindent This structure stores the patch-time-frequency {\it look up table}.
The fields are:

\begin{description}
\item[\texttt{INT2    timeIndex }]  Time index of the {\sc lut}.
\item[\texttt{INT8    f0Bin }]  Frequency bin for which it has been
constructed.
\item[\texttt{REAL8   deltaF }]  Frequency resolution \texttt{df=1/TCOH}, where \texttt{1/TCOH}
is the coherent integration time used in teh demodulation procedure.
\item[\texttt{INT4    nFreqValid }] Number of frequencies where the {\sc lut} is
valid.
\item[\texttt{INT4    iniBin }]  First bin affecting the patch with respect to
\verb@f0@.
\item[\texttt{INT4    nBin }]  Exact number of bins affecting the patch.
\item[\texttt{INT4    offset }]  Frequency bin corresponding to center of patch measured with respect to \texttt{f0Bin} (zero in modulated case)
\item[\texttt{UINT2   maxNBins}] Maximum number of bins affecting the patch (for
                               memory allocation purposes).
\item[\texttt{UINT2   maxNBorders}]  Maximum number of borders affecting the
				patch (for memory allocation purposes).
\item[\texttt{HOUGHBorder      *border} ]  The annulus borders.
\item[\texttt{HOUGHBin2Border  *bin} ]  Bin to border
correspondence.
\end{description}

\begin{verbatim}
struct HOUGHPatchGrid
\end{verbatim}
\index{\texttt{HOUGHPatchGrid}}

\noindent This structure stores patch-frequency {\it grid} information. The fields are:

\begin{description}
\item[\texttt{REAL8   f0 }]  Frequency to construct grid.
\item[\texttt{REAL8   deltaF }]  Frequency resolution: \texttt{df=1/TCOH}.
\item[\texttt{REAL8   patchSizeX }] Patch size in radians along x-axis.
\item[\texttt{REAL8   patchSizeY }] Patch size in radians along y-axis.
\item[\texttt{REAL8   minWidthRatio }]  Ratio between the minimum  annulus width
for this search, and the minimun  annulus width for  1 year integration time.
This value should be in the interval  [1.0, 25.0].
\item[\texttt{REAL8   deltaX }] Longitudinal space resolution, x-direction.
\item[\texttt{REAL8   xMin }]  Patch limit, as the coordinate of the center of the first pixel.
\item[\texttt{REAL8   xMax }]  Patch limit, as the coordinate of the center  of the last pixel.
\item[\texttt{UINT2   xSide }] Real number of pixels in the x direction (in the
projected plane). It should be less than or equal to \texttt{xSideMax}.
\item[\texttt{UINT2   xSideMax }]  Maximun number of pixels in the x direction
                        (for memory allocation), i.e. length of \verb@xCoor@.
\item[\texttt{REAL8   *xCoor }] Coordinates of the pixel centers.
\item[\texttt{REAL8   deltaY }] Longitudinal space resolution, y-direction.
\item[\texttt{REAL8   yMin }]  Patch limit, as center of the first pixel.
\item[\texttt{REAL8   yMax }]  Patch limit, as center of the last pixel.
\item[\texttt{UINT2   ySide }] Real number of pixels in the y-direction (in the
projected plane). It should be less than or equal to \texttt{ySideMax}.
\item[\texttt{UINT2   ySideMax }]  Maximun number of pixels in the y direction
                        (for memory allocation), i.e. length of \verb@yCoor@.
\item[\texttt{REAL8   *yCoor }] Coordinates of the pixel centers.
\end{description}


\begin{verbatim}
struct HOUGHResolutionPar
\end{verbatim}
\index{\texttt{HOUGHResolutionPar}}

 \noindent This structure holds the parameters needed for gridding the patch
The fields are:

\begin{description}
\item[\texttt{REAL8   f0 }]  Frequency at which construct the grid.
\item[\texttt{REAL8   deltaF }]  Frequency resolution: \texttt{df=1/TCOH}.
\item[\texttt{REAL8   patchSizeX }] Patch size in radians along x-axis.
\item[\texttt{REAL8   patchSizeY }] Patch size in radians along y-axis.
\item[\texttt{REAL8   minWidthRatio }] Ratio between the minimum  annulus width
for this search and the minimun  annulus width for one year integration time.
This value should be in the interval  [1.0, 25.0].
\end{description}

\begin{verbatim}
struct REAL8Cart3Coor
\end{verbatim}
\index{\texttt{REAL8Cart3Coor }}

\noindent Three dimensional Cartessian coordinates.
The fields are:

\begin{description}
\item[\texttt{REAL8  x }]
\item[\texttt{REAL8  y }]
\item[\texttt{REAL8  z }]
\end{description}

\begin{verbatim}
struct REAL8Cart2Coor
\end{verbatim}
\index{\texttt{REAL8Cart2Coor}}

\noindent Two dimensional Cartessian coordinates.
The fields are:

\begin{description}
\item[\texttt{REAL8  x }]
\item[\texttt{REAL8  y }]
\end{description}


\begin{verbatim}
struct REAL8Polar2Coor
\end{verbatim}
\index{\texttt{REAL8Polar2Coor}}

\noindent Two dimensional polar coordinates.
The fields are:

\begin{description}
\item[\texttt{REAL8  alpha }]
\item[\texttt{REAL8  radius }]
\end{description}


\begin{verbatim}
struct REAL8UnitPolarCoor
\end{verbatim}
\index{\texttt{REAL8UnitPolarCoor}}

\noindent Polar coordinates of a unitary vector on the sphere.
The fields are:

\begin{description}
\item[\texttt{REAL8  alpha }] Any value
\item[\texttt{REAL8  delta }]  In the interval [$-\pi/2, \,  \pi/2$]
\end{description}

\begin{verbatim}
struct HOUGHDemodPar
\end{verbatim}
\index{\texttt{HOUGHDemodPar}}

\noindent Demodulation parameters needed for the Hough transform. All
coordinates are assumed to be with respect to the same reference system. The fields are:

\begin{description}
\item[\texttt{REAL8               deltaF }]: Frequency resolution: \texttt{df=1/TCOH}.
\item[\texttt{REAL8UnitPolarCoor  skyPatch }]:  $N_{center}$ (alpha, delta):
position of the center of the patch.
\item[\texttt{REAL8   patchSizeX }] Patch size in radians along x-axis.
\item[\texttt{REAL8   patchSizeY }] Patch size in radians along y-axis.
\item[\texttt{REAL8Cart3Coor      veloC }]:  $v(t)/c$ (x,y,z): Relative detector
velocity
\item[\texttt{REAL8Cart3Coor      positC }]: $(x(t)-x(t0))/c$ (x,y,z): Position
of the detector.
\item[\texttt{REAL8               timeDiff }]:  $T_{\hat N} (t)-T_{\hat N} (\hat t_0)$:
 Time difference.
\item[\texttt{REAL8Vector         spin }]: Spin down information. It includes
the fields: \texttt{length}: maximum order of
spin-down parameter, and
 \texttt{*data}: pointer to spin-down parameter set $F_k$.
\end{description}

\begin{verbatim}
struct HOUGHParamPLUT
\end{verbatim}
\index{\texttt{HOUGHParamPLUT}}

\noindent Parameters needed to construct the partial look up table. The fields are:

\begin{description}
\item[\texttt{INT8             f0Bin }] Frequency  bin for which it has been constructed
\item[\texttt{REAL8            deltaF }] Frequency resolution: \texttt{df=1/TCOH}.
\item[\texttt{REAL8UnitPolarCoor   xi }] Center of the circle on the celestial
sphere, $\xi$(alpha,delta) in the rotated coordinates.
\item[\texttt{REAL8            cosDelta }] $\Delta \cos(\phi)$ corresponding to
one annulus.
\item[\texttt{INT4             offset }] Frequency bin corresponding to center of patch; measured w.r.t. \texttt{f0Bin}.
\item[\texttt{INT4             nFreqValid }] Number of frequency bins for which the LUT is valid.
\item[\texttt{REAL8            cosPhiMax0 }] $\max(\cos(\phi))$ of the
\texttt{f0Bin}.
\item[\texttt{REAL8            cosPhiMin0 }]  $\min(\cos(\phi))$ of the
\texttt{f0Bin}.
\item[\texttt{REAL8            epsilon }] maximum angle (distance in radians) from the pole
to consider  a circle as a line in the projected plane.
\end{description}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{LUTHV}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage\input{StereographicC}
\newpage\input{PatchGridC}
\newpage\input{ParamPLUTC}
\newpage\input{NDParamPLUTC}
\newpage\input{ConstructPLUTC}

%%%%%%%%%%Test program. %%
\newpage\input{TestConstructPLUTC}
\newpage\input{TestNDConstructPLUTC}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

</lalLaTeX> */



/*
 * 4.  Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _LUT_H
#define _LUT_H

/*
 * 5. Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
# include <stdlib.h>
# include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

/*
 *  #include "LALRCSID.h"
 *  not needed, it is included in "LALConstants.h"
 */



/*
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * 6. Assignment of Id string using NRCSID()
 */

NRCSID (LUTH, "$Id$");

/*
 * 7. Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/* <lalErrTable file="LUTHErrorTable"> */

#define LUTH_ENULL 1
#define LUTH_ESIZE 2
#define LUTH_ESZMM 4
#define LUTH_EINT  6
#define LUTH_ESAME 8
#define LUTH_EFREQ 10
#define LUTH_EVAL 12

#define LUTH_MSGENULL "Null pointer"
#define LUTH_MSGESIZE "Invalid input size"
#define LUTH_MSGESZMM "Size mismatch"
#define LUTH_MSGEINT  "Invalid interval"
#define LUTH_MSGESAME "Input/Output data vectors are the same"
#define LUTH_MSGEFREQ "Invalid frequency"
#define LUTH_MSGEVAL  "Invalid value"

/* </lalErrTable>  */


/* ******************************************************
 * 8. Macros. But, note that macros are deprecated.
 *    They could be moved to the modules where are needed
 */


#define MAX(A, B)  (((A) < (B)) ? (B) : (A))
#define MIN(A, B)  (((A) < (B)) ? (A) : (B))
#define cot(A)  (1./tan(A))



/* *******************************************************
 * 9. Constant Declarations. (discouraged)
 */

/** Maximum ``error'' (as a fraction of the width of the thinnest annulus)
    which allows to represent a circle by  a line. */
#define LINERR     0.001

/** Maximum ``error'' (as a fraction of the width of the thinnest annulus)
    which allows to consider two border equivalents!
    It is relevant for determining the LUT frequency range validity */
#define PIXERR     0.5


/** Width of the thinnest annulus in terms of pixels */
/* #define PIXELFACTOR  2 */

/** Earth v_epicycle/c  TO BE CHANGED DEPENDING ON DETECTOR */
/* #define VEPI 1.0e-06 */
/* vepi = 2*pi*R_earth /(day *c) */
#define VEPI LAL_TWOPI * LAL_REARTH_SI / ( LAL_DAYSID_SI * LAL_C_SI )

/** Total detector velocity/c TO BE CHANGED DEPENDING ON DETECTOR */
/* #define VTOT 1.06e-04 */
/* vtot = 2*pi* 1AU / (year * c) */
#define VTOT LAL_TWOPI * LAL_AU_SI / ( LAL_YRSID_SI * LAL_C_SI )


/* **************************************************************
 * 10. Structure, enum, union, etc., typdefs.
 */

typedef INT2 COORType;
 /* typedef INT4 COORType; */
  /* typedef  UCHAR COORType; */

  /** This structure stores the border of a circle clipped on the projected plane. */
typedef struct tagHOUGHBorder{
  INT4  yUpper;    /**< upper y pixel affected by this border */
  INT4  yLower;    /**< lower y pixel affected by this border and
		      yUpper<yLower or yUpper<0 are possible */
  INT4  yCenter;   /**< y pixel value of the center of the circle */
  UINT2     ySide; /**< length of xPixel */
  COORType *xPixel; /**< x pixel index to be marked */
} HOUGHBorder;


  /** This structure stores the border indexes corresponding to one frequency
      bin plus the corrections to be added to the first column of the patch. */
typedef struct tagHOUGHBin2Border{
  INT2   leftB1;     /**< index  of the border[xxx] to be used */
  INT2   rightB1;
  INT2   leftB2;
  INT2   rightB2;
  INT2   piece1max;  /**< interval limits of the first column to be add */
  INT2   piece1min;
  INT2   piece2max;
  INT2   piece2min;
} HOUGHBin2Border;

  /** Patch-Time-Frequency Look Up table*/
typedef struct tagHOUGHptfLUT{
  INT2    timeIndex;  /**< time index of the LUT */
  INT8    f0Bin;      /**< freq. bin for which it has been constructed */
  REAL8   deltaF;     /**< df=1/TCOH */
  INT8    nFreqValid; /**< number of frequencies where the LUT is valid */
  INT4    iniBin;     /**< first bin affecting the patch with respect to f0 */
  INT4    nBin;       /**< number of bins affecting the patch */
  INT4    offset;      /**< freq. bin (wrt f0Bin) containing center of patch  */
  UINT2   maxNBins;    /**< maximum number of bins affecting the patch. For
                               memory allocation */
  UINT2   maxNBorders; /**< maximum number of borders affecting the patch. For
                               memory allocation */
  HOUGHBorder      *border; /**< the annulus borders */
  HOUGHBin2Border  *bin;    /**< Bin to Border correspondence */
} HOUGHptfLUT;


  /** Patch-Frequency Grid information */
typedef struct tagHOUGHPatchGrid{
  REAL8   f0;         /**< frequency to construct grid */
  REAL8   deltaF;     /**< df=1/TCOH */
  REAL8   deltaX;
  REAL8   xMin;     /**< patch limits, as centers of the last pixels */
  REAL8   xMax;
  UINT2   xSide;    /**< number of pixels in the x direction (projected plane)*/
  REAL8   *xCoor;   /**< coordinates of the pixel centers */
  REAL8   deltaY;
  REAL8   yMin;     /**< patch limits,as centers of the last pixels */
  REAL8   yMax;
  UINT2   ySide;    /**< number of pixels in the y direction */
  REAL8   *yCoor;   /**< coordinates of the pixel centers  */
} HOUGHPatchGrid;

  /** parameters needed for gridding the patch */
typedef struct tagHOUGHResolutionPar{
  INT8    f0Bin; /**< frequency bin */
  REAL8   deltaF;        /**< df=1/TCOH */
  REAL8   patchSkySizeX;     /**< Size of sky patch in radians */
  REAL8   patchSkySizeY;
  REAL8   pixelFactor; /**< number of pixel that fit in the thinnest annulus*/
  REAL8   pixErr;   /**< for validity of LUT as PIXERR */
  REAL8   linErr;   /**< as LINERR circle ->line */
  REAL8   vTotC;    /**< estimate value of v-total/C as VTOT */
} HOUGHResolutionPar;

  /** required for constructing patch */
typedef struct tagHOUGHSizePar{
  INT8    f0Bin; /**< corresponding freq. bin  */
  REAL8   deltaF;        /**< df=1/TCOH */
  REAL8   deltaX; /**< pixel size in the projected plane */
  REAL8   deltaY;
  UINT2   xSide;    /**< number of pixels in the x direction (projected plane)*/
  UINT2   ySide;    /**< number of pixels in the y direction */
  UINT2   maxNBins;    /**< maximum number of bins affecting the patch. For
                               memory allocation */
  UINT2   maxNBorders; /**< maximum number of borders affecting the patch. For
                               memory allocation */
  INT8    nFreqValid; /**< number of frequencies where the LUT is valid */
  REAL8   epsilon; /**< max. angle (rad.) from the pole to consider
			       a circle as a line in the projected plane */
} HOUGHSizePar;

  /** cartesian coordinates on a plane */
typedef struct tagREAL8Cart3Coor{
  REAL8  x;
  REAL8  y;
  REAL8  z;
} REAL8Cart3Coor;

  /** cartesian coordinates on a plane */
typedef struct tagREAL8Cart2Coor{
  REAL8  x;
  REAL8  y;
} REAL8Cart2Coor;

  /** polar coordinates on a plane */
typedef struct tagREAL8Polar2Coor{
  REAL8  alpha;
  REAL8  radius;
} REAL8Polar2Coor;

  /** right ascension and declination */
typedef struct tagREAL8UnitPolarCoor{
  REAL8  alpha;  /**< any value */
  REAL8  delta;  /**< \f$ -\pi/2, \pi/2 \f$ */
} REAL8UnitPolarCoor;

  /** Parameters required for constructing LUT */
typedef struct tagHOUGHParamPLUT{
  INT8             f0Bin;   /**< freq. bin for which it has been constructed */
  REAL8            deltaF;  /**< df=1/TCOH */
  REAL8UnitPolarCoor   xi;  /**< \f$ \xi(\alpha,\delta)\f$ in rotated coordinates */
  REAL8            cosDelta;    /**< \f$ \Delta \cos(\phi)\f$ for one annulus */
  INT4             offset;
  INT8             nFreqValid;
  REAL8            cosPhiMax0;
  REAL8            cosPhiMin0;
  REAL8            epsilon; /**< max. angle (rad.) from the pole to consider
			       a circle as a line in the projected plane */
} HOUGHParamPLUT;

  /** Demodulation parameters needed for the Hough transform (all in same coordinate system */
typedef struct tagHOUGHDemodPar{
  /* all coordinates with respect the same reference system */
  REAL8               deltaF;   /**<   df=1/TCOH */
  REAL8UnitPolarCoor  skyPatch; /**<   N_center \f$(alpha, delta)\f$ */
  REAL8Cart3Coor      veloC;    /**<   \f$ v(t)/c {x,y,z} \f$ */
  REAL8Cart3Coor      positC;   /**<  \f$ (x(t)-x(t_0))/c {x,y,z} \f$*/
  REAL8               timeDiff; /**<  \f$ T(t)-T(t_0) \f$ */
  REAL8Vector         spin; /**< spindown parameters */
} HOUGHDemodPar;



/*
 * 11. Extern Global variables. (discouraged)
 */



/*
 * 12. Functions Declarations (i.e., prototypes).
 */
void LALHOUGHComputeSizePar (LALStatus  *status, /* demod case */
                   HOUGHSizePar        *out,
                   HOUGHResolutionPar  *in
		   );

void LALHOUGHComputeNDSizePar (LALStatus  *status, /* non -demod case*/
                   HOUGHSizePar          *out,
                   HOUGHResolutionPar    *in
		   );


void LALHOUGHFillPatchGrid (LALStatus   *status,
		   HOUGHPatchGrid      *out,
                   HOUGHSizePar        *par
		       );

void LALHOUGHParamPLUT (LALStatus   *status, /* Demod. case */
                   HOUGHParamPLUT   *out,  /* parameters needed build LUT*/
		   HOUGHSizePar     *sizePar,
                   HOUGHDemodPar    *par  /* demodulation parameters */
			);

void LALNDHOUGHParamPLUT (LALStatus   *status, /* non-demod. case */
                   HOUGHParamPLUT   *out,  /* parameters needed build LUT*/
		   HOUGHSizePar     *sizePar,
                   HOUGHDemodPar    *par  /* demodulation parameters */
			);

void LALRotatePolarU(LALStatus            *status,
		     REAL8UnitPolarCoor   *out,
		     REAL8UnitPolarCoor   *in,
		     REAL8UnitPolarCoor   *par
		     );

void LALInvRotatePolarU(LALStatus            *status,
			REAL8UnitPolarCoor   *out,
			REAL8UnitPolarCoor   *in,
			REAL8UnitPolarCoor   *par
			);

void LALStereoProjectPolar(LALStatus           *status,
			   REAL8Polar2Coor     *out,
			   REAL8UnitPolarCoor  *in
			   );

void LALStereoProjectCart(LALStatus           *status,
			  REAL8Cart2Coor      *out,
			  REAL8UnitPolarCoor  *in
			  );

void LALStereoInvProjectPolar(LALStatus           *status,
			      REAL8UnitPolarCoor  *out,
			      REAL8Polar2Coor     *in
			      );

void LALStereoInvProjectCart(LALStatus           *status,
			     REAL8UnitPolarCoor  *out,
			     REAL8Cart2Coor      *in
			     );

void LALHOUGHConstructPLUT(LALStatus       *status,
			   HOUGHptfLUT     *lut,
			   HOUGHPatchGrid  *patch,
			   HOUGHParamPLUT  *par
			   );


#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection _LUT_H */













