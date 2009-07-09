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

/** \file  ConstructPLUT.c
 * \ingroup pulsarHough
 * \author Sintes, A and Krishnan, B
 * \brief Core routines for constructing the Hough Look-Up-Tables
 * \date $Date$
 * Revision: $Id$
 *
 * History:   Created by Sintes June 7, 2001
 *            Modified by Badri Krishnan Feb 2003
 *-----------------------------------------------------------------------
 *
\par Description

This module is the core of the Hough transform. The LAL function
LALHOUGHConstructPLUT()
constructs the look up tables that will be used to build the so-called
partial-Hough maps. Each look up table is valid for a given sky-patch, time, and
frequency range around a certain  \verb@f0@ value. The look up table contains
all the necessary information regarding the borders of the annuli clipped on
the \lq projected' two dimensional sky-patch.

The inputs are:  HOUGHPatchGrid   containing the grid patch
information. This is independent of the sky location of the
patch, And  HOUGHParamPLUT  with all the other parameters needed.

The output is the look up table  HOUGHptfLUT


\par Uses
\code
PLUTInitialize(HOUGHptfLUT  *)
FillPLUT(HOUGHParamPLUT *, HOUGHptfLUT *, HOUGHPatchGrid *)
CheckLeftCircle(REAL8, REAL8, REAL8, INT4 *, INT4 *, INT4 *,
                HOUGHPatchGrid *)
CheckRightCircle(REAL8, REAL8, REAL8, INT4 *, INT4 *, INT4 *,
                 HOUGHPatchGrid *)
DrawRightCircle(REAL8, REAL8, REAL8, INT4, INT4, COORType *,
                HOUGHPatchGrid *)
DrawLeftCircle(REAL8, REAL8, REAL8, INT4, INT4, COORType *,
               HOUGHPatchGrid *)
CheckLineCase(REAL8, REAL8, REAL8, REAL8 *, INT4 *)
FindExactLine(REAL8, REAL8, REAL8 *, REAL8 *)
FindLine(REAL8, REAL8, REAL8, REAL8 *, REAL8 *)
CheckLineIntersection(REAL8, REAL8, REAL8, INT4 *, INT4 *, INT4 *,
                      HOUGHPatchGrid *)
DrawLine(REAL8, REAL8, REAL8,INT4, INT4, COORType *, HOUGHPatchGrid *)
Fill1Column(INT4, UINT4*, HOUGHptfLUT *, HOUGHPatchGrid *)
FillCaseN1(INT4, INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *)
FillCaseN2(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *)
FillCaseN3(INT4, INT4, INT4, INT4 *, HOUGHptfLUT *, HOUGHPatchGrid *)
FillCaseN4(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *)
FillCaseN5(INT4, INT4, INT4, HOUGHptfLUT *)
FillCaseN6(INT4, INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *)
FillCaseN7(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *)
FillCaseN8(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *)
Fill1ColumnAnor(INT4, HOUGHptfLUT *, HOUGHPatchGrid *)
FillCaseA1(INT4, INT4, INT4, HOUGHptfLUT *)
FillCaseA2(INT4, INT4, INT4, HOUGHptfLUT *)
FillCaseA3(INT4, INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *)
InitialCircleCase(UINT4 *,REAL8, REAL8, REAL8, REAL8 *, INT4 *, INT4 *,
                  HOUGHptfLUT *, HOUGHPatchGrid *)
SecondCircleCase(INT4, UINT4*, REAL8, REAL8, REAL8, INT4, REAL8 *,
                 INT4*,INT4 *,INT4*, HOUGHptfLUT *, HOUGHPatchGrid *)
FollowCircleCase(INT4,UINT4 *,REAL8,REAL8,REAL8,REAL8,REAL8,INT4 *,
                 INT4 *,INT4 *, HOUGHptfLUT *, HOUGHPatchGrid *)
InitialLineCase(UINT4 *, REAL8, REAL8, REAL8, INT4 *, HOUGHptfLUT *,
                HOUGHPatchGrid *)
SecondLineCase(INT4, UINT4 *, REAL8, REAL8, REAL8, INT4 *, HOUGHptfLUT *,
               HOUGHPatchGrid *)
FollowLineCase(INT4, UINT4 *, REAL8, REAL8, REAL8, REAL8, INT4, INT4 *,
                           HOUGHptfLUT *, HOUGHPatchGrid *)
\endcode


*/


/*
 * 1.  An author and Id block
 */

/************************************ <lalVerbatim file="ConstructPLUTCV">
Author: Sintes, A. M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */

/*
 * 2. Commented block with the documetation of this module
 */


/* <lalLaTeX>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Module \texttt{ConstructPLUT.c}}
\label{ss:ConstructPLUT.c}
Construction of the look up table for generating partial Hough maps assuming the
use of the stereographic projection.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ConstructPLUTD}
\index{\verb&LALHOUGHConstructPLUT()&}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Macros (used only internally)}

\begin{verbatim}
#define MAX(A, B)  (((A) < (B)) ? (B) : (A))
#define MIN(A, B)  (((A) < (B)) ? (A) : (B))
#define cot(A)  (1./tan(A))
#define rint(x) floor((x)+0.5)
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Static function declarations}
\begin{verbatim}
static void PLUTInitialize(HOUGHptfLUT  *);
static void FillPLUT(HOUGHParamPLUT *, HOUGHptfLUT *, HOUGHPatchGrid *);
static void CheckLeftCircle(REAL8, REAL8, REAL8, INT4 *, INT4 *, INT4 *,
                            HOUGHPatchGrid *);
static void CheckRightCircle(REAL8, REAL8, REAL8, INT4 *, INT4 *, INT4 *,
                             HOUGHPatchGrid *);
static void DrawRightCircle(REAL8, REAL8, REAL8, INT4, INT4, COORType *,
                            HOUGHPatchGrid *);
static void DrawLeftCircle(REAL8, REAL8, REAL8, INT4, INT4, COORType *,
                           HOUGHPatchGrid *);
static void CheckLineCase(REAL8, REAL8, REAL8, REAL8 *, INT4 *);
static void FindExactLine(REAL8, REAL8, REAL8 *, REAL8 *);
static void FindLine(REAL8, REAL8, REAL8, REAL8 *, REAL8 *);
static void CheckLineIntersection(REAL8, REAL8, REAL8, INT4 *, INT4 *, INT4 *,
                                  HOUGHPatchGrid *);
static void DrawLine(REAL8, REAL8, REAL8,INT4, INT4, COORType *, HOUGHPatchGrid *);
static void Fill1Column(INT4, UINT4*, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN1(INT4, INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN2(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN3(INT4, INT4, INT4, INT4 *, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN4(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN5(INT4, INT4, INT4, HOUGHptfLUT *);
static void FillCaseN6(INT4, INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN7(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN8(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void Fill1ColumnAnor(INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseA1(INT4, INT4, INT4, HOUGHptfLUT *);
static void FillCaseA2(INT4, INT4, INT4, HOUGHptfLUT *);
static void FillCaseA3(INT4, INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void InitialCircleCase(UINT4 *,REAL8, REAL8, REAL8, REAL8 *, INT4 *, INT4 *,
                              HOUGHptfLUT *, HOUGHPatchGrid *);
static void SecondCircleCase(INT4, UINT4*, REAL8, REAL8, REAL8, INT4, REAL8 *,
                             INT4*,INT4 *,INT4*, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FollowCircleCase(INT4,UINT4 *,REAL8,REAL8,REAL8,REAL8,REAL8,INT4 *,
                             INT4 *,INT4 *, HOUGHptfLUT *, HOUGHPatchGrid *);
static void InitialLineCase(UINT4 *, REAL8, REAL8, REAL8, INT4 *, HOUGHptfLUT *,
                            HOUGHPatchGrid *);
static void SecondLineCase(INT4, UINT4 *, REAL8, REAL8, REAL8, INT4 *, HOUGHptfLUT *,
                           HOUGHPatchGrid *);
static void FollowLineCase(INT4, UINT4 *, REAL8, REAL8, REAL8, REAL8, INT4, INT4 *,
                           HOUGHptfLUT *, HOUGHPatchGrid *);
\end{verbatim}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

This module is the core of the Hough transform. The LAL function \verb&LALHOUGHConstructPLUT()&
constructs the look up tables that will be used to build the so-called
partial-Hough maps. Each look up table is valid for a given sky-patch, time, and
frequency range around a certain  \verb@f0@ value. The look up table contains
all the necessary information regarding the borders of the annuli clipped on
the \lq projected' two dimensional sky-patch.

The inputs are:  \verb@HOUGHPatchGrid  *patch@ containing the grid patch
information. This is independent of the sky location of the
patch. And  \verb@HOUGHParamPLUT  *par@ with all the other parameters needed.

The output is:  \verb@HOUGHptfLUT   *lut@. The fields are:

\begin{description}
\item[\texttt{lut-> timeIndex }]  Time index of the {\sc lut}.
\item[\texttt{lut->f0Bin }]  Frequency bin for which it has been
constructed.
\item[\texttt{lut->deltaF }]  Frequency resolution: \texttt{df=1/TCOH}
\item[\texttt{lut->nFreqValid }]  \verb@ = PIXERR * f0Bin * VEPI / VTOT@:
Number of frequencies where the {\sc lut} is
valid.
\item[\texttt{lut->iniBin }]  First bin affecting the patch with respect to
\verb@f0@.
\item[\texttt{lut->nBin }]  Exact number of bins affecting the patch.
\item[\texttt{lut->maxNBins }] Maximum number of bins affecting the patch (for
memory allocation purposes), i.e. length of \texttt{lut->bin}.
\item[\texttt{lut->maxNorders }] Maximum number of borders affecting the patch (for
memory allocation purposes), i.e. length of \texttt{lut->border}.
\item[\texttt{lut->border} ]  The annulus borders.
\item[\texttt{lut->bin} ]  Bin to border correspondence.
\end{description}



More detailed documentation can be found in the source code itself.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
%\begin{verbatim}
%LALZDestroyVector()
%\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{ConstructPLUTCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
</lalLaTeX> */


/*
 * 3. Includes. These should appear in the following order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */


#include <lal/LUT.h>


/*
 * 4. Assignment of Id string using NRCSID()
 */

NRCSID (CONSTRUCTPLUTC, "$Id$");



/*
 * 5.a) Constants, structures (used only internally in this module)
 */

/*
 * 5.b) Type declarations (used only internally)
 */

/*
 * 5.c) Macros (used only internally)
 */

#define MAX(A, B)  (((A) < (B)) ? (B) : (A))
#define MIN(A, B)  (((A) < (B)) ? (A) : (B))
#define cot(A)  (1./tan(A))
#define rint(x) floor((x)+0.5)


/*
 * 5.d) Extern global variable declarations (strongly discouraged)
 */


/*
 * 5.e) static Global variables
 */


/*
 * 5.f) Static function declarations
 */

static void PLUTInitialize(HOUGHptfLUT  *);
static void FillPLUT(HOUGHParamPLUT *, HOUGHptfLUT *, HOUGHPatchGrid *);

static void CheckLeftCircle(REAL8,REAL8,REAL8,INT4 *,INT4 *,INT4 *,
                          HOUGHPatchGrid *);
static void CheckRightCircle(REAL8,REAL8,REAL8,INT4 *,INT4 *,INT4 *,
                           HOUGHPatchGrid *);
static void DrawRightCircle(REAL8,REAL8,REAL8,INT4,INT4, COORType *,
                          HOUGHPatchGrid *);
static void DrawLeftCircle(REAL8,REAL8,REAL8,INT4,INT4,COORType *,
                         HOUGHPatchGrid *);
static void CheckLineCase(REAL8, REAL8, REAL8,REAL8 *, INT4 *);
/* static void FindExactLine(REAL8,REAL8,REAL8 *,REAL8 *); */
static void FindLine(REAL8,REAL8,REAL8,REAL8 *,REAL8 *);
static void CheckLineIntersection(REAL8,REAL8,REAL8,INT4 *,INT4 *,INT4 *,
                                HOUGHPatchGrid *);
static void DrawLine(REAL8, REAL8, REAL8,INT4, INT4, COORType *, HOUGHPatchGrid *);
static void Fill1Column(INT4, UINT4*, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN1(INT4, INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN2(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN3(INT4, INT4, INT4, INT4 *, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN4(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN5(INT4, INT4, INT4, HOUGHptfLUT *);
static void FillCaseN6(INT4, INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN7(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseN8(INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void Fill1ColumnAnor(INT4, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FillCaseA1(INT4, INT4, INT4, HOUGHptfLUT *);
static void FillCaseA2(INT4, INT4, INT4, HOUGHptfLUT *);
static void FillCaseA3(INT4, INT4, INT4, HOUGHptfLUT *, HOUGHPatchGrid *);

static void InitialCircleCase(UINT4 *,REAL8, REAL8, REAL8, REAL8 *, INT4 *, INT4 *,
                            HOUGHptfLUT *, HOUGHPatchGrid *);
static void SecondCircleCase(INT4, UINT4*,REAL8,REAL8, REAL8,INT4,REAL8*,
                           INT4*,INT4 *,INT4*, HOUGHptfLUT *, HOUGHPatchGrid *);
static void FollowCircleCase(INT4,UINT4 *,REAL8,REAL8,REAL8,REAL8,REAL8,INT4 *,
                           INT4 *,INT4 *, HOUGHptfLUT *, HOUGHPatchGrid *);
static void InitialLineCase(UINT4 *, REAL8, REAL8, REAL8, INT4 *, HOUGHptfLUT *,
                          HOUGHPatchGrid *);
static void SecondLineCase(INT4, UINT4 *, REAL8, REAL8, REAL8, INT4 *, HOUGHptfLUT *,
                         HOUGHPatchGrid *);
static void FollowLineCase(INT4, UINT4 *,REAL8, REAL8, REAL8,REAL8,INT4,INT4 *,
                         HOUGHptfLUT *, HOUGHPatchGrid *);


/*
 * 5.g)  The functions that make up the guts of this module
 */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*  ***************************** <lalVerbatim file="ConstructPLUTD"> */
void LALHOUGHConstructPLUT(LALStatus       *status,
			   HOUGHptfLUT     *lut,
			   HOUGHPatchGrid  *patch,
			   HOUGHParamPLUT  *par)
{ /*  ************************************************ </lalVerbatim> */

  INT8    f0Bin;

  /* --------------------------------------------- */
  INITSTATUS (status, " LALHOUGHConstructPLUT", CONSTRUCTPLUTC);
  /*  ATTATCHSTATUSPTR (status); */

  /*   Make sure the arguments are not NULL: */
  /* use aborts instead of asserts */
  if (lut == NULL) {
    ABORT( status, LUTH_ENULL, LUTH_MSGENULL);
  }
  /* ASSERT (lut, status, LUTH_ENULL, LUTH_MSGENULL); */

  if (patch == NULL) {
    ABORT( status, LUTH_ENULL, LUTH_MSGENULL);
  }
  /* ASSERT (patch, status, LUTH_ENULL, LUTH_MSGENULL); */

  if (par == NULL) {
    ABORT( status, LUTH_ENULL, LUTH_MSGENULL);
  }
  /* ASSERT (par,  status, LUTH_ENULL, LUTH_MSGENULL); */

  if (  fabs((REAL4)par->deltaF - (REAL4)patch->deltaF) > 1.0e-6) {
    ABORT( status, LUTH_EVAL, LUTH_MSGEVAL);
  }
  /* ASSERT (par->deltaF == patch->deltaF,  status, LUTH_EVAL, LUTH_MSGEVAL); */

  /* -------------------------------------------   */

  f0Bin =  par->f0Bin;

  lut->deltaF = par->deltaF;
  lut->f0Bin  = f0Bin;

  lut->nFreqValid = par->nFreqValid;
  /* lut->nFreqValid = PIXERR * f0Bin *VEPI/VTOT; */

  /* -------------------------------------------   */

  PLUTInitialize(lut);
  FillPLUT(par, lut, patch);


  /* Make sure number of bins makes sense with the dimensions  */
  if (lut->nBin <= 0) {
    ABORT( status, LUTH_ESIZE, LUTH_MSGESIZE);
  }
  /* ASSERT (lut->nBin >0, status, LUTH_ESIZE, LUTH_MSGESIZE); */

  if (lut->nBin > lut->maxNBins) {
    ABORT( status, LUTH_ESIZE, LUTH_MSGESIZE);
  }
  /* ASSERT (lut->nBin <= lut->maxNBins, status, LUTH_ESIZE, LUTH_MSGESIZE); */


  /* -------------------------------------------   */
  /* DETATCHSTATUSPTR (status); */

  /* normal exit */
  RETURN (status);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/*        >>>>>>>>     FUNCTIONS  ONLY  USED INTERNALLY     <<<<<<    */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */




/***************************************************************/
/*     Subroutine to initialize the partial LUT                */
/*     Authors:  Sintes, A.M                                   */
/*               June 7, 2001                                  */
/***************************************************************/

static void  PLUTInitialize(HOUGHptfLUT  *lut){
  UINT4  i;
  UINT4 maxNBins, maxNBorders;

  maxNBins = lut->maxNBins;
  maxNBorders = lut->maxNBorders;

  for(i=0;i<maxNBins;++i){

    /* check lut->bin is not null */
    if ( lut->bin + i == NULL ) {
      fprintf(stderr," pointer lut->bin+%d is null [ConstructPLUT.c %d]\n", i, __LINE__);
    }

    lut->bin[i].leftB1  = 0;
    lut->bin[i].rightB1 = 0;
    lut->bin[i].leftB2  = 0;
    lut->bin[i].rightB2 = 0;
    lut->bin[i].piece1max = -1;
    lut->bin[i].piece1min = 0;
    lut->bin[i].piece2max = -1;
    lut->bin[i].piece2min = 0;
  }

  for(i=0;i<maxNBorders;++i){

    /* check lut->border is not null */
    if ( lut->border + i == NULL ) {
      fprintf(stderr," pointer lut->border+%d is null [ConstructPLUT.c %d]\n", i, __LINE__);
    }

    lut->border[i].yUpper = -1;
    lut->border[i].yLower = 0;
  }
  return;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/****************************************************************/
/*       subroutine to fill the  LUT      			*/
/*                                                        	*/
/*     Authors:  Sintes, A.M             			*/
/*		 March 16, 2001		     	                */
/****************************************************************/

/* ******************* Some explanations *********************** */
/* The program should go like this: */

/* First call subroutine XXXXXXXXXX.c */
/*    to calculate xi(t) for a given f0 (and demodulation parameters). */
/*    Then rotate that vector to our new coordinates, if it was not  */
/*    constructed in that way. */

/*    We need here |xi|, alpha, delta[-pi/2,pi/2],  on the sphere. */
/*       or its cartesian coordiantes  */

/*       The master equation to solve is: */
/*           f(t)-f0 = xi*n - xi*N */
/*        or */
/*           cos( phi)= [f(t)-f0 +xi*N]/|xi| */

/*        where:       xi*N= -|xi|* sin(delta) */

/*        Given a bin in the peakgram with */
/*                  Deltaf= f(t)(center bin)- f0= k df  */
/*                         (df => frequency resolution =1/TCOH) */
/*        Then */
/*           [cos(phi)]max = min[ 1, k (df/|xi|) + (xi*N +df/2)/|xi| ] */
/*        and */
/*           [cos(phi)]min = max[-1, k (df/|xi|) + (xi*N -df/2)/|xi| ] */
/*       Note  */
/*         if  [cos(phi)]max = 1,  do  not increase k !!! */
/*         if  [cos(phi)]min =-1,  do  not decrease k !!! */

/*           [cos(phi)]max(k) = min[ 1, df/|xi|+[cos(phi)]max(k-1) ] */
/*           [cos(phi)]min(k) = max[ 1, df/|xi|+[cos(phi)]min(k-1) ]  */
/* ********************************************************************* */
/*     for a bin,peak or border: */

/*     if cos(phi)=+-1  => nothing to be drawn! */
/*     otherwise */
/*         obtain phi in the interval [0, LAL_PI] */
/*         if delta+phi= pi/2, or delta-phi=pi/2 (carefull, with +2*pi) */
/*                   draw line */
/*         otherwise draw a circle */
/* ******************************************************************* */
/* a pixel j , corresponds to x(j)= patch->deltaX*(0.5+j- patch->xSide/2.) */
/* the nearest center of a pixel is j = round[ x/patch->deltaX +patch->xSide/2.-0.5]  */
/* ******************************************************************* */
/* The way to convert k into the bin index is the following:    */
/*    for k>=0, binindex= k                                     */
/*    for k<0,  binindex= nBin+iniBin-1-k                       */
/*   f0 corresponds to k=0                                      */
/* This will be used when reading the peakgram!                 */
/* ******************************************************************* */

static void  FillPLUT(HOUGHParamPLUT  *par, HOUGHptfLUT  *lut,
                    HOUGHPatchGrid  *patch){

  /********************************************************/
  /*    variables that need to be calculated before       */
  /********************************************************/

  REAL8 cosDelta;   /* = df/|xi|                                             */
  REAL8 cosPhiMax0; /* max val of cosPhi of freq bin containing patch center */
  REAL8 cosPhiMin0; /* mix val of cosPhi of freq bin containing patch center */
  REAL8 alpha;      /* = xi.alpha in the rotated coordinates                 */
  REAL8 delta;      /* = xi.delta                                            */
  REAL8 epsilon;    /* = 8 * LINERR/(f0Bin * VTOT * patchSize^2)             */


  /********************************************************/
  UINT4 maxNBins;
  UINT4 maxNBorders;

  UINT4 lastBorder =0;  /* counter of the last build border */
  UINT4 currentBin =0;  /* counter of the bin studied       */

  INT4 ifailPlus  = 1; /* =1 (ok, continue to next bin), =0 (stop) */
  INT4 ifailMinus = 1; /* =1 (ok, continue to previous bin), =0 (stop) */

  INT4 directionPlus=-1;  /* = +1,or -1 */
  INT4 directionPlusZero=-1;  /* = +1,or -1 */
  INT4 directionMinus; /* = +1,or -1 */

  INT4 pathology=1;   /* =1 (normal), =0 (anormal case) */
  INT4 lineCase;    /* =1 line Case, =0 circle case */
  INT4 nBinPos;

  /********************************************************/

  REAL8 rCritic;
  REAL8 lambda;
  REAL8 rcOldPlus;
  REAL8 rcOldMinus;
  REAL8 cosPhiMax, cosPhiMin, phi;
  REAL8 ang1,ang2;
  REAL8 eps;

  /********************************************************/
  maxNBins = lut->maxNBins;
  maxNBorders = lut->maxNBorders;

  alpha = par->xi.alpha;
  delta = par->xi.delta;
  cosDelta = par->cosDelta;
  cosPhiMax0 = par->cosPhiMax0;
  cosPhiMin0 = par->cosPhiMin0;
  epsilon = par->epsilon;

  /********************************************************/

  /* Copy value of offset */
  lut->offset = par->offset;

  /********************************************************/

  lambda = 2* delta -LAL_PI*0.5;
  rCritic = 2* cos(lambda) /(1 - sin(lambda) );

  /* Initializing variables to irrelevant values,
     since these values should never be used */
  rcOldPlus = rCritic;
  rcOldMinus = rCritic;
  /********************************************************/
  /*  starting with the (central) bin corresponding to:   */
  /*            Delta_f(t) =0 (k=0), border cosPhiMax     */
  /********************************************************/

  /* cosPhiMax = MIN(1, cosPhiMax0); */
  cosPhiMax = cosPhiMax0;

  if(cosPhiMax > 0.99999999){ /* avoid points or small circles */
    ifailPlus = 0;            /* do not go to the next bin */
    directionPlus = -1;
  } else{

    phi  = acos(cosPhiMax);   /* in the interval (0,PI) */
    ang1 = delta + phi;
    ang2 = delta - phi;

    /* check for lines, or numerical lines! */
    CheckLineCase(epsilon, ang1, ang2, &eps, &lineCase);

    if( lineCase ){
 	/* line case */
      InitialLineCase(&lastBorder,alpha,delta,eps, &ifailPlus,lut, patch);
      directionPlus = -1;
    } else{
      /* circle case */
      InitialCircleCase(&lastBorder,alpha, ang1, ang2,
			&rcOldPlus, &directionPlus, &ifailPlus,lut, patch);
    }
  }

  directionPlusZero= directionPlus;

  /********************************************************/
  /* moving to the other bins                             */
  /********************************************************/

  /********************************************************/
  /* Plus direction: increasing cosPhiMax  */
  /********************************************************/



  while (ifailPlus){

    ++currentBin;
    pathology = 1;

    /* Some checks */
    /* #ifdef CHECKHOUGHINDEX */
    if (currentBin > maxNBins || lastBorder>= maxNBorders ){
      fprintf(stderr,"currentBin=%d not in range 1 to maxNBins=%d\n"
	      "or lastborder=%d >= maxNBorders=%d [ConstructPLUT.c %d]\n",
	      currentBin,maxNBins,lastBorder,maxNBorders, __LINE__);
    }
    /* #endif */

    cosPhiMax = cosPhiMax + cosDelta;
    /* or cosPhiMax = MIN(1,cosPhiMax + cosDelta ); */

    if( cosPhiMax > 0.99999999){ /* check appropiate value */
      ifailPlus = 0;
    } else {

      phi  = acos(cosPhiMax); /* it is in the interval (0,pi) */
      ang1 = delta + phi;
      ang2 = delta - phi;

      /* check for lines, or numerical lines! */
      CheckLineCase(epsilon, ang1, ang2, &eps, &lineCase);

      if( lineCase ){
 	/* line case */
	FollowLineCase(currentBin, &lastBorder,alpha,delta,eps,
		       rcOldPlus, directionPlus, &ifailPlus,lut, patch);
      } else{
	/* circle case */
	FollowCircleCase(currentBin, &lastBorder,alpha, ang1, ang2,rCritic,
			 rcOldPlus, &pathology, &directionPlus,
			 &ifailPlus, lut, patch);
      }
    }

    if(pathology){
      Fill1Column(currentBin, &lastBorder, lut, patch);
    }else{
      Fill1ColumnAnor(currentBin, lut, patch);
    }
  }

  nBinPos = currentBin;


  /********************************************************/
  /*  starting with the (central) bin corresponding to:   */
  /*            Delta_f(t) =0 (k=0), border cosPhiMin     */
  /********************************************************/

  /* cosPhiMin = MAX(-1, cosPhiMin0); */
  cosPhiMin = cosPhiMin0;

  if(cosPhiMin < -0.99999999){ /* avoid points or small circles */
    ifailMinus = 0;            /* do not go to the next bin */
  } else{

    phi  = acos(cosPhiMin);   /* in the interval (0,PI) */
    ang1 = delta + phi;
    ang2 = delta - phi;

    /* check for lines, or numerical lines! */
    CheckLineCase(epsilon, ang1, ang2, &eps, &lineCase);

    if( lineCase ){
      /* line case */
      SecondLineCase(currentBin,&lastBorder,alpha,delta,eps, &ifailMinus,
                     lut, patch);
      directionMinus = -1;
      pathology = 0;
    } else{
      /* circle case */
      pathology = 1; /* provisionally */
      SecondCircleCase(currentBin, &lastBorder,alpha, ang1, ang2,
		       directionPlusZero, &rcOldMinus,
		       &pathology, &directionMinus,
		       &ifailMinus, lut, patch);
    }
  }

  /********************************************************/
  /*  the way to identify initial pathologies:            */
  /*    initial or second case being a line !             */
  /*    or two circles, both with direction=-1            */
  /********************************************************/

  if(pathology){
    Fill1Column(0, &lastBorder,lut, patch);
  }else{
    Fill1ColumnAnor(0, lut, patch);
  }

  /********************************************************/
  /* moving to the other bins                             */
  /********************************************************/

  /********************************************************/
  /* Minus direction: decreasing cosPhiMin                */
  /********************************************************/

  while(ifailMinus){

    ++currentBin;
    pathology = 1;

    /* Some checks */
    /* #ifdef CHECKHOUGHINDEX */
    if (currentBin > maxNBins || lastBorder>= maxNBorders ){
      fprintf(stderr,"currentBin=%d not in range 1 to maxNBins=%d\n"
	      "or lastborder=%d >= maxNBorders=%d [ConstructPLUT.c %d]\n",
	      currentBin,maxNBins,lastBorder,maxNBorders, __LINE__);
    }
    /* #endif */


    cosPhiMin = cosPhiMin - cosDelta;

    if( cosPhiMin < -0.99999999){ /* check appropiate value */
      ifailMinus = 0;
    } else {

      phi  = acos(cosPhiMin); /* it is in the interval (0,pi) */
      ang1 = delta + phi;
      ang2 = delta - phi;

      /* check for lines, or numerical lines! */
      CheckLineCase(epsilon, ang1, ang2, &eps, &lineCase);

      if( lineCase ){
    	/* line case */
	FollowLineCase(currentBin, &lastBorder,alpha,delta,eps,
		       rcOldMinus, directionMinus, &ifailMinus,lut, patch);
      } else{
	/* circle case */
	FollowCircleCase(currentBin, &lastBorder,alpha, ang1, ang2,rCritic,
			 rcOldMinus, &pathology, &directionMinus,
			 &ifailMinus, lut, patch);
      }
    }

    if(pathology){
      Fill1Column(currentBin, &lastBorder, lut, patch);
    }else{
      Fill1ColumnAnor(currentBin, lut, patch);
    }
  }

  /********************************************************/
  /* set iniBin,nBin  etc */
  /********************************************************/

  lut->nBin = currentBin + 1;
  lut->iniBin = nBinPos - currentBin;

  return;
}

/* end of the subroutine*/
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/****************************************************************/
/*								*/
/*     from initialLineCase.c          March 16, 2001          	*/
/*                                                         	*/
/*       program to find and identify bins and borders      	*/
/*	 in the initial line case                             	*/
/*                                                        	*/
/*     Authors:  Sintes, A.M             			*/
/*								*/
/****************************************************************/


/****************************************************************/

static void InitialLineCase(UINT4  *lastBorderP, REAL8 alpha, REAL8  delta,
		     REAL8 eps, INT4 *ifailP, HOUGHptfLUT  *lut,
		     HOUGHPatchGrid  *patch){


  INT4 lastBorder;


  REAL8  xA,yA;
  INT4 yymin,yymax;
  REAL8 xRel, slope;
  INT4 noIn;   /* if no intersection occurs noIn=0 */

  /* some paranoid error checking */
  if (patch == NULL ) {
    fprintf(stderr, "pointer patch is null [ConstructPLUT.c %d]\n", __LINE__);
  }

  if (lut == NULL ) {
    fprintf(stderr, "pointer lut is null [ConstructPLUT.c %d]\n",__LINE__);
  }

  if (ifailP == NULL ) {
    fprintf(stderr, "pointer ifailP is null [ConstructPLUT.c %d]\n",__LINE__);
  }

  if (lastBorderP == NULL ) {
    fprintf(stderr, "pointer lastBorderP is null [ConstructPLUT.c %d]\n",__LINE__);
  }
  lastBorder = *lastBorderP;



 /************************************************/
  FindLine(alpha, delta, eps, &xA, &yA);
  CheckLineIntersection(alpha, xA, yA, &yymin, &yymax, &noIn, patch);

  if( noIn ==0 ){
    *ifailP = 0;    /* =1 (ok), =0 (stop) */
    return;
  }
  ++lastBorder;

  if (yymin < 0) {
    fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
	    yymin, __LINE__);
    yymin = 0;
  }
  if (yymax >= patch->ySide) {
    fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
	    yymax, patch->ySide-1, __LINE__);
    yymax = patch->ySide - 1;
  }

  lut->border[lastBorder].yUpper = yymax;
  lut->border[lastBorder].yLower = yymin;

  DrawLine(alpha, xA, yA, yymin, yymax,
	   &lut->border[lastBorder].xPixel[0] , patch);

  /************************************************/

  if( (alpha == LAL_PI*0.5) || (alpha == LAL_PI*1.5) ||
      (alpha == -LAL_PI*0.5) ){       /* horizontal line */

    if( yA < 0 ){ /* convention */
      lut->bin[0].rightB1 = lastBorder;
      lut->bin[1].leftB1  = lastBorder;
      lut->border[lastBorder].yCenter = 0; /*or smaller*/
    } else {
      lut->bin[0].leftB2  = lastBorder;
      lut->bin[1].rightB2 = lastBorder;
      lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
    }
  }
  else { /* non- horizontal line */
    xRel = xA + tan(alpha)*yA;

    if( (alpha == 0) || (alpha == LAL_PI) || (alpha == -LAL_PI) ){
      /* vertical line */
      slope = +1; /* arbitrary number, just to avoid overflow */
    } else {
      slope = - cot(alpha);
    }

    if( xRel < 0 ){
      lut->bin[0].leftB2 = lastBorder;
      lut->bin[1].rightB2= lastBorder;

      if ( slope < 0 ){
	lut->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
      }

    } else {
      lut->bin[0].rightB1 = lastBorder;
      lut->bin[1].leftB1  = lastBorder;

      if ( slope > 0 ){
	lut->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
      }
    }
  } /* end non-horizontal line */


  /************************************************/

  *lastBorderP =  lastBorder;

 return;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/****************************************************************/
/*								*/
/*     from secondLineCase.c          March 16, 2001          	*/
/*                                                         	*/
/*       program to find and identify bins and borders      	*/
/*	 in the second line case                             	*/
/*                                                        	*/
/*     Authors:  Sintes, A.M             			*/
/*								*/
/****************************************************************/


/****************************************************************/

static void SecondLineCase(INT4 currentBin, UINT4  *lastBorderP,
		    REAL8 alpha, REAL8  delta,
		     REAL8 eps, INT4 *ifailP, HOUGHptfLUT *lut,
		     HOUGHPatchGrid  *patch){

  /* we are for sure in a pathological case. Border names are
     changed accordingly */

  INT4 lastBorder;


  REAL8  xA,yA;
  INT4 yymin,yymax;
  REAL8 xRel, slope;
  INT4 noIn;   /* if no intersection occurs noIn=0 */


  /* some paranoid error checking */
  if (patch == NULL ) {
    fprintf(stderr, "pointer patch is null [ConstructPLUT.c %d]\n",__LINE__);
  }

  if (lut == NULL ) {
    fprintf(stderr, "pointer lut is null [ConstructPLUT.c %d]\n",__LINE__);
  }

  if (ifailP == NULL ) {
    fprintf(stderr, "pointer ifailP is null [ConstructPLUT.c %d]\n",__LINE__);
  }

  if (lastBorderP == NULL ) {
    fprintf(stderr, "pointer lastBorderP is null [ConstructPLUT.c %d]\n",__LINE__);
  }
  lastBorder = *lastBorderP;

 /************************************************/
  FindLine(alpha, delta, eps, &xA, &yA);
  CheckLineIntersection(alpha, xA, yA, &yymin, &yymax, &noIn, patch);

  if( noIn ==0 ){
    *ifailP = 0;  /* =1 (ok), =0 (stop) */
    return;
  }
  ++lastBorder;

  if (yymin < 0) {
    fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
	    yymin, __LINE__);
    yymin = 0;
  }
  if (yymax >= patch->ySide) {
    fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
	    yymax, patch->ySide-1, __LINE__);
    yymax = patch->ySide - 1;
  }

  lut->border[lastBorder].yUpper = yymax;
  lut->border[lastBorder].yLower = yymin;

  DrawLine(alpha, xA, yA, yymin, yymax,
	   &lut->border[lastBorder].xPixel[0], patch );

  /************************************************/
  /* all are pathological cases. The code is similar to the
     InitialLineCase, except for the modifications:
     rightB1 -> rightB2, leftB2 -> leftB1, where marked */

  if( (alpha == LAL_PI*0.5) || (alpha == LAL_PI*1.5) ||
      (alpha == -LAL_PI*0.5) ){       /* horizontal line */

    if( yA < 0 ){ /* convention */
      lut->bin[0].rightB2= lastBorder; /* modified */
      lut->bin[currentBin+1].leftB1 = lastBorder;
      lut->border[lastBorder].yCenter = 0; /*or smaller*/
    } else {
      lut->bin[0].leftB1 = lastBorder; /* modified */
      lut->bin[currentBin+1].rightB2= lastBorder;
      lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
    }
  }
  else { /* non- horizontal line */
    xRel = xA + tan(alpha)*yA;

    if( (alpha == 0) || (alpha == LAL_PI) || (alpha == -LAL_PI) ){
      /* vertical line */
      slope = +1; /* arbitrary number, just to avoid overflow */
    } else {
      slope = - cot(alpha);
    }

    if( xRel < 0 ){
      lut->bin[0].leftB1   = lastBorder; /* modified */
      lut->bin[currentBin+1].rightB2= lastBorder;

      if ( slope < 0 ){
	lut->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
      }

    } else {
      lut->bin[0].rightB2  = lastBorder; /* modified */
      lut->bin[currentBin+1].leftB1 = lastBorder;

      if ( slope > 0 ){
	lut->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
      }
    }
  } /* end non-horizontal line */


  /************************************************/

  *lastBorderP =  lastBorder;

  /************************************************/
 return;
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/****************************************************************/
/*								*/
/*     from followLineCase.c          March 16, 2001          	*/
/*                                                         	*/
/*       program to find and identify bins and borders      	*/
/*	 in the line case                                	*/
/*                                                        	*/
/*     Authors:  Sintes, A.M             			*/
/*								*/
/****************************************************************/


/****************************************************************/

static void FollowLineCase(INT4 currentBin, UINT4  *lastBorderP,
		    REAL8 alpha, REAL8  delta, REAL8 eps,REAL8 rcOld,
		    INT4 direction, INT4 *ifailP, HOUGHptfLUT *lut,
		    HOUGHPatchGrid  *patch){

  INT4 lastBorder;


  REAL8  xA,yA;
  INT4 yymin,yymax;
  INT4 noIn;   /* if no intersection occurs noIn=0 */


  /* some paranoid error checking */
  if (patch == NULL ) {
    fprintf(stderr, "pointer patch is null [ConstructPLUT.c %d]\n",__LINE__);
  }

  if (lut == NULL ) {
    fprintf(stderr, "pointer lut is null [ConstructPLUT.c %d]\n",__LINE__);
  }

  if (ifailP == NULL ) {
    fprintf(stderr, "pointer ifailP is null [ConstructPLUT.c %d]\n",__LINE__);
  }

  if (lastBorderP == NULL ) {
    fprintf(stderr, "pointer lastBorderP is null [ConstructPLUT.c %d]\n",__LINE__);
  }
  lastBorder = *lastBorderP;

  /************************************************/
  FindLine(alpha, delta, eps, &xA, &yA);
  CheckLineIntersection(alpha, xA, yA, &yymin, &yymax, &noIn, patch);

  if( noIn ==0 ){
    *ifailP = 0; /* =1 (ok), =0 (stop) */
    return;
  }
  ++lastBorder;

  if (yymin < 0) {
    fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
	    yymin, __LINE__);
    yymin = 0;
  }
  if (yymax >= patch->ySide) {
    fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
	    yymax, patch->ySide-1, __LINE__);
    yymax = patch->ySide - 1;
  }

  lut->border[lastBorder].yUpper = yymax;
  lut->border[lastBorder].yLower = yymin;

  DrawLine(alpha, xA, yA, yymin, yymax,
	   &lut->border[lastBorder].xPixel[0] , patch);

  /************************************************/
  if( direction == 1 ){
    /* This means that the initial case was a circle
       that deformed into a line when increasing radius,
       an old value of rcOld does exist */
    REAL8 xcOld, ycOld, xRel, slope;

    xcOld = rcOld *cos(alpha);
    ycOld = rcOld *sin(alpha);

    if( (alpha == LAL_PI*0.5) || (alpha == LAL_PI*1.5) ||
	(alpha == -LAL_PI*0.5) ){       /* horizontal line */

      lut->bin[currentBin].leftB1   = lastBorder;
      lut->bin[currentBin+1].rightB1= lastBorder;

      /* alternatively, one can also set two extra borders */
      /*  lut->bin[currentBin].rightB2  = lastBorder; */
      /*  lut->bin[currentBin+1].leftB2 = lastBorder; */

      if( ycOld < yA ){
	lut->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
      }

    }
    else { /* non- horizontal line */
      xRel = xA + tan(alpha)*(yA -ycOld);

      if( (alpha == 0) || (alpha == LAL_PI) || (alpha == -LAL_PI) ){
	/* vertical line */
	slope = +1; /* arbitrary number, just to avoid overflow */
      } else {
	slope = - cot(alpha);
      }

      if( xRel < xcOld ){
	lut->bin[currentBin].leftB1   = lastBorder;
	lut->bin[currentBin+1].rightB1= lastBorder;

	if ( slope > 0 ){
	  lut->border[lastBorder].yCenter = 0; /*or smaller*/
	} else {
	  lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
	}

      } else {
	lut->bin[currentBin].rightB2  = lastBorder;
	lut->bin[currentBin+1].leftB2 = lastBorder;

	if ( slope < 0 ){
	  lut->border[lastBorder].yCenter = 0; /*or smaller*/
	} else {
	  lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
	}
      }
    } /* end non-horizontal line */
  } /* end of +1 direction*/
  else {
    /************************************************/
    /*  This means ( direction == -1)               */
    /************************************************/
    /* This means the first case was a line.
       Latter on, it can be transformed into circles,
       when radius decreases. No value of rcOld set */

    REAL8 xRel, slope;

    if( (alpha == LAL_PI*0.5) || (alpha == LAL_PI*1.5) ||
	(alpha == -LAL_PI*0.5) ){       /* horizontal line */

      if( yA < 0 ){ /* convention */
	lut->bin[currentBin].rightB1  = lastBorder;
	lut->bin[currentBin+1].leftB1 = lastBorder;
	lut->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	lut->bin[currentBin].leftB2   = lastBorder;
	lut->bin[currentBin+1].rightB2= lastBorder;
	lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
      }
    }
    else { /* non- horizontal line */
      xRel = xA + tan(alpha)*yA;

      if( (alpha == 0) || (alpha == LAL_PI) || (alpha == -LAL_PI) ){
	/* vertical line */
	slope = +1; /* arbitrary number, just to avoid overflow */
      } else {
	slope = - cot(alpha);
      }

      if( xRel < 0 ){
	lut->bin[currentBin].leftB2   = lastBorder;
	lut->bin[currentBin+1].rightB2= lastBorder;

	if ( slope < 0 ){
	  lut->border[lastBorder].yCenter = 0; /*or smaller*/
	} else {
	  lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
	}

      } else {
	lut->bin[currentBin].rightB1  = lastBorder;
	lut->bin[currentBin+1].leftB1 = lastBorder;

	if ( slope > 0 ){
	  lut->border[lastBorder].yCenter = 0; /*or smaller*/
	} else {
	  lut->border[lastBorder].yCenter = patch->ySide -1;/*or bigger */
	}
      }
    } /* end non-horizontal line */

  }/* end of -1 direction*/


  /************************************************/

  *lastBorderP =  lastBorder;

  return;
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/****************************************************************/
/*         from lines.c         March 19, 2001                	*/
/*                                                         	*/
/*   Subroutines check, find and clip lines on the patch       	*/
/*		                                        	*/
/*                                                        	*/
/*     Authors:   Sintes, A.M           			*/
/*								*/
/****************************************************************/


/* First version: to be improved thinking about
   substituting "ceil" and "floor" for the corresponding values
   assuming positive arguments. */

/****************************************************************/
/* Subroutine to check if we are in the line case, numerical
   lines inluded */
/****************************************************************/

static void CheckLineCase(REAL8 epsilon, REAL8 ang1, REAL8 ang2,
		    REAL8 *epsP, INT4 *lineCaseP){

  /* lineCaseP =1 line Case, =0 circle case */

  /* exact lines occur if  */
  /*   ((ang1==LAL_PI/2.) || (ang2==LAL_PI/2.) ||(ang2==-3.*LAL_PI/2.)) */
  /* other cases correspond to numerical lines */

  REAL8 e1,e21,e22;


  /* some paranoid error checking */
  if (epsP == NULL ) {
    fprintf(stderr, "pointer epsP is null [ConstructPLUT.c %d]\n",__LINE__);
  }

  if (lineCaseP == NULL ) {
    fprintf(stderr, "pointer lineCaseP is null [ConstructPLUT.c %d]\n",__LINE__);
  }



  e1  = ang1 - LAL_PI*0.5;
  e21 = ang2 - LAL_PI*0.5;
  e22 = e21 + 2*LAL_PI;

  if(fabs(e1) < epsilon ){
    *lineCaseP = 1;
    *epsP = e1;
    return;
  }

  if( fabs(e21) < epsilon ){
    *lineCaseP = 1;
    *epsP = e21;
    return;
  }

  if( fabs(e22) < epsilon ){
    *lineCaseP = 1;
    *epsP = e22;
    return;
  }


  *lineCaseP = 0;
  return;
}


/****************************************************************/
/* Subroutine to find the parameters of the exact line!!!       */
/*      The line equation:                                      */
/*                y = cotg(alpha) [xA - x ] +yA                 */
/****************************************************************/

#if 0 /* NOT USED */
static void FindExactLine(REAL8 alpha, REAL8 delta,
		   REAL8 *xA, REAL8 *yA){

  REAL8 lambda, rA;

  lambda =  2.*delta -LAL_PI*0.5;

  if (fabs(1.-sin(lambda) < 1.0e-6) {
    fprintf(stderr, "warning: possible division by 0 [ConstructPLUT.c %d]\n",__LINE__);
  }
  rA = 2.* cos(lambda)/(1.-sin(lambda) );
  *xA = rA* cos(alpha);
  *yA = rA* sin(alpha);

return;
}
#endif

/***************************************************************/
/* Subroutine to find the parameters of 'approximated' lines   */
/*      The line equation:                                     */
/*                 y = cotg(alpha) [xA - x ] +yA               */
/***************************************************************/

static void FindLine(REAL8 alpha, REAL8 delta, REAL8 eps,
	      REAL8 *xA, REAL8 *yA){

  REAL8 lambda, rA;

  lambda =  2.*delta -LAL_PI*0.5-eps;
  rA = 2.* cos(lambda)/(1.-sin(lambda) );
  *xA = rA* cos(alpha);
  *yA = rA* sin(alpha);

return;
}


/*********************************************************/
/* Subroutine to check if the line intersects the patch. */
/*     Output:                                           */
/*          information if intersects the patch:         */
/*                 noIn =0 no intersection               */
/* 		   noIn =1 intersection                  */
/*          and the y-pixel range of the intersection    */
/*********************************************************/

static void CheckLineIntersection(REAL8 alpha, REAL8 xA, REAL8 yA,
		      INT4 *yyminP, INT4 *yymaxP, INT4 *noInP,
		      HOUGHPatchGrid  *patch){

  INT4 yymin,yymax,noIn;
  volatile REAL4 kkk;

  yymin = 0;
  yymax = 0;
  noIn = 0;

  if ((alpha == 0) || (alpha == LAL_PI) || (alpha == -LAL_PI) ){
    /* vertical line */
    if ((patch->xMin <= xA) && (patch->xMax >= xA)){
      noIn  = 1;
      yymin = 0;
      yymax = patch->ySide-1;
    }
  }
  else{
    if ( (alpha == LAL_PI*0.5) || (alpha == LAL_PI*1.5) ||
	 (alpha == -LAL_PI*0.5) ){
      /* horizontal line */
      if ((patch->yMin <= yA) && (patch->yMax >= yA)){
	noIn  = 1;
	/* yymin =  ceil((REAL4)(yA/patch->deltaY-0.5) + (REAL4)(patch->ySide*0.5)); */
	kkk =  yA/patch->deltaY-0.5;
	kkk += patch->ySide*0.5;
	yymin =  ceil(kkk);

	/* yymax = floor((REAL4)(yA/patch->deltaY-0.5) +(REAL4)(patch->ySide*0.5)); */
	kkk =  yA/patch->deltaY-0.5;
	kkk += patch->ySide*0.5;
	yymax = floor(kkk);
	/* Note  yymax < yymin,   to identify an horizontal line!
	   If the area to mark is above (arriba), mark yymin and higher.
	   If the area to mark is below (abajo), mark yymax and lower.*/
      }
    }
    else{
	/* generic case */
      REAL8 myy1,y2,slope;
      slope = cot(alpha);
      myy1  = slope*(xA-patch->xMin)+yA;
      y2  = slope*(xA-patch->xMax)+yA;

      if (  ( (myy1 >= patch->yMin) || (y2 >= patch->yMin) )
	    && ( (myy1 <= patch->yMax) || (y2 <= patch->yMax) ) ){
	REAL8 yupper,ylower;
	noIn   = 1;
	yupper = MAX(myy1,y2);
	ylower = MIN(myy1,y2); /* or  ylower=myy1+y2-yupper  */
	yupper = MIN(yupper,patch->yMax);
	ylower = MAX(ylower,patch->yMin);
	/* yymin  = ceil((REAL4)(ylower/patch->deltaY -0.5)+(REAL4)(patch->ySide*0.5)); */
	kkk =  ylower/patch->deltaY -0.5;
	kkk += patch->ySide*0.5;
	yymin  = ceil(kkk);
	/* yymax  =floor((REAL4)(yupper/patch->deltaY-0.5)+(REAL4)(patch->ySide*0.5)); */
	kkk =  yupper/patch->deltaY-0.5;
	kkk += patch->ySide*0.5;
	yymax = floor(kkk);
	/* in case of a pseudo-horizontal line (yymax < yymin) */
	/* we use the same convention as in the horizontal case */
      }
    }
  }

  *yyminP = yymin;
  *yymaxP = yymax;
  *noInP  = noIn;

  return;
}


/*****************************************************************/
/*   Subroutine to draw non-horizontal lines!                    */
/*****************************************************************/
static void DrawLine(REAL8 alpha, REAL8 xA, REAL8 yA,
	      INT4 yymin, INT4 yymax, COORType  *column,
	      HOUGHPatchGrid  *patch){

  INT4 jj;
  volatile REAL4 kkk;

  column[yymin] = patch->xSide;
  column[yymax] = patch->xSide;

  /* the if-else, is not really needed, just to avoid repeating
     the ceil when not necessary! */

 if ((alpha == 0) || (alpha == LAL_PI)){
   /* vertical line */
   INT4  xpixel;

   /* xpixel = ceil( (REAL4)(xA/patch->deltaX-0.5) +(REAL4)(patch->xSide*0.5)); */
   kkk =  xA/patch->deltaX-0.5;
   kkk += patch->xSide*0.5;
   xpixel = ceil(kkk);

   if (xpixel < 0) {
     xpixel = 0;
   }
   if (xpixel > patch->xSide) {
     xpixel = patch->xSide;
   }

   for(jj=yymin;jj<=yymax;++jj){
     column[jj] = xpixel;
   }
 }
 else{
   /* remaining cases */
   REAL8 tanalpha;
   REAL8 xofy;
   INT4  xpixel;

   tanalpha = tan(alpha);

   for(jj=yymin;jj<=yymax;++jj){
     /* will not be executed in the horizontal case */
     xofy = xA + tanalpha*( yA- patch->yCoor[jj] );
     kkk =  xofy/patch->deltaX-0.5;
     kkk += patch->xSide*0.5;
     xpixel = ceil(kkk);

     if (xpixel < 0) {
       xpixel = 0;
     }
     if (xpixel > patch->xSide) {
       xpixel = patch->xSide;
     }

     column[jj] = xpixel;
   }
 }

 return;
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/****************************************************************/
/*								*/
/*     from initialCircleCase.c          March 16, 2001       	*/
/*                                                         	*/
/*       program to find and identify bins and borders      	*/
/*	 in the circle case                                	*/
/*                                                        	*/
/*     Authors:  Sintes, A.M             			*/
/*								*/
/****************************************************************/



/****************************************************************/

static void InitialCircleCase(UINT4  *lastBorderP, REAL8 alpha,
		      REAL8 ang1, REAL8  ang2,
		      REAL8 *rcOldP, INT4 *directionP,
		      INT4 *ifailP, HOUGHptfLUT *lut,
		      HOUGHPatchGrid  *patch){

  INT4 lastBorder;
  INT4 direction;    /* +1, or -1 */


  REAL8 rho1,rho2,radius;
  REAL8 xc,yc,rc; /* coordinates of the center of the circle */
  INT4 pieces;


  lastBorder = *lastBorderP;

  rho1 = cos(ang1)/(1. -sin(ang1));
  rho2 = cos(ang2)/(1. -sin(ang2));
  rc   = rho1 + rho2;  /* positive or negative */
  xc   = rc*cos(alpha);
  yc   = rc*sin(alpha);
  radius = fabs( rho1 - rho2 ); /* abs could be avoided */

 /************************************************/
  /* check for exclusion of intersection of the circle with the patch*/
  if( ( yc+radius < patch->yMin ) || ( yc-radius > patch->yMax )  ||
      ( xc+radius < patch->xMin ) || ( xc-radius > patch->xMax )  ||
      ( sqrt(patch->xMax*patch->xMax + patch->yMax*patch->yMax) + radius < fabs(rc) ) ){
    /* no intersection */
    *ifailP = 0;  /* =1 (ok), =0 (stop) */
    return;
  }

  /************************************************/
  /* possible intersection case:                  */
  /************************************************/

  /* direction = +1 (+,-) (left, right) */
  /* direction = -1 (-,+) (right, left) */
  /************************************************/
  direction = ( (fabs(rc) < radius)  ? (1) : (-1) );
  /************************************************/


  if( direction == -1){
    pieces = 0;

    /* check left circle */
    if( xc >= patch->xMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
		floor( yc/patch->deltaY + patch->ySide*0.5);
	           /*  rint( yc/patch->deltaY + patch->ySide*0.5 -0.5); */
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0] , patch);
	lut->bin[0].rightB1= lastBorder;
	lut->bin[1].leftB1 = lastBorder;
      }
    }

    /* check right circle */
    if( xc <= patch->xMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	              floor( yc/patch->deltaY + patch->ySide*0.5);
		      /*rint( yc/patch->deltaY + patch->ySide*0.5 -0.5); */
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0], patch );
	lut->bin[0].leftB2  = lastBorder;
	lut->bin[1].rightB2 = lastBorder;
      }
    }

    /* real no intersection case */
    if( pieces ==0 ){
      *ifailP = 0;
      return;
    }

  } else {

    /************************************************/
    /*  This means ( direction == +1)               */
    /************************************************/
    /* no  pathologies   */

    pieces = 0;
    /* check left circle */
    if( xc >= patch->xMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	            floor( yc/patch->deltaY + patch->ySide*0.5);
		    /* rint( yc/patch->deltaY + patch->ySide*0.5 -0.5); */
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0] , patch);

	lut->bin[0].leftB1 = lastBorder;
	lut->bin[1].rightB1= lastBorder;
      }
    }

    /* check right circle */
    if( xc <= patch->xMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	            floor( yc/patch->deltaY + patch->ySide*0.5);
		    /*  rint( yc/patch->deltaY + patch->ySide*0.5 -0.5); */
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0], patch );

	lut->bin[0].rightB2= lastBorder;
	lut->bin[1].leftB2 = lastBorder;
      }
    }

    /* real no intersection case */
    if( pieces ==0){
      *ifailP = 0;
      return;
    }

  }

 /************************************************/
  *lastBorderP =  lastBorder;
  *directionP = direction;
  *rcOldP = rc;


  return;
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/****************************************************************/
/*								*/
/*     from secondCircleCase.c          March 16, 2001         	*/
/*                                                         	*/
/*       program to find and identify bins and borders      	*/
/*	 in the circle case                                	*/
/*                                                        	*/
/*     Authors:  Sintes, A.M             			*/
/*								*/
/****************************************************************/


/****************************************************************/

static void SecondCircleCase(INT4 currentBin, UINT4  *lastBorderP,
		      REAL8 alpha, REAL8 ang1, REAL8  ang2,
		      INT4 directionPlus, REAL8 *rcOldP,
		      INT4 *pathologyP,INT4 *directionP,
		      INT4 *ifailP, HOUGHptfLUT *lut,
		      HOUGHPatchGrid  *patch){

  INT4 lastBorder;
  INT4 pathology;    /* =1 (normal), =0 (anormal) */
  INT4 direction;    /* +1, or -1 */


  REAL8 rho1,rho2,radius;
  REAL8 xc,yc,rc; /* coordinates of the center of the circle */
  INT4 pieces;

  pathology = *pathologyP;
  lastBorder = *lastBorderP;

  rho1 = cos(ang1)/(1. -sin(ang1));
  rho2 = cos(ang2)/(1. -sin(ang2));
  rc   = rho1 + rho2;  /* positive or negative */
  xc   = rc*cos(alpha);
  yc   = rc*sin(alpha);
  radius = fabs( rho1 - rho2 ); /* abs could be avoided */

 /************************************************/
  /* check for exclusion of intersection of the circle with the patch*/
  if( ( yc+radius < patch->yMin ) || ( yc-radius > patch->yMax )  ||
      ( xc+radius < patch->xMin ) || ( xc-radius > patch->xMax )  ||
      ( sqrt(patch->xMax*patch->xMax + patch->yMax*patch->yMax) + radius < fabs(rc) ) ){
    /* no intersection */
    *ifailP = 0;  /* =1 (ok), =0 (stop) */
    return;
  }

 /************************************************/
  /* possible intersection case:                  */
  /************************************************/

  /* direction = +1 (+,-) (left, right) */
  /* direction = -1 (-,+) (right, left) */
  /************************************************/
  direction = ( (fabs(rc) < radius)  ? (1) : (-1) );
  /************************************************/


  if( direction == -1){ /* possible pathologies */

    if( directionPlus == -1){ pathology = 0; }

    pieces = 0;

    /* check left circle */
    if( xc >= patch->xMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	           floor( yc/patch->deltaY + patch->ySide*0.5);
		   /*   rint( yc/patch->deltaY + patch->ySide*0.5 -0.5); */
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0] , patch);
	if ( pathology ){
	  lut->bin[0].rightB1 = lastBorder;
	} else{
	  lut->bin[0].rightB2 = lastBorder; /* modified */
	}
	lut->bin[currentBin+1].leftB1 = lastBorder;
      }
    }

    /* check right circle */
    if( xc <= patch->xMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	             floor( yc/patch->deltaY + patch->ySide*0.5);
		     /* rint( yc/patch->deltaY + patch->ySide*0.5 -0.5);*/
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0], patch );
	if ( pathology ){
	  lut->bin[0].leftB2 = lastBorder;
	} else{
	  lut->bin[0].leftB1 = lastBorder; /*modified */
	}
	lut->bin[currentBin+1].rightB2 = lastBorder;
      }
    }

    /* real no intersection case */
    if( pieces ==0 ){
      *ifailP = 0;
      return;
    }

  } else {

    /************************************************/
    /*  This means ( direction == +1)               */
    /************************************************/
    /* no  pathologies   */

    pieces = 0;
    /* check left circle */
    if( xc >= patch->xMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	           floor( yc/patch->deltaY + patch->ySide*0.5);
		   /*   rint( yc/patch->deltaY + patch->ySide*0.5 -0.5); */
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0] , patch);

	lut->bin[0].leftB1   = lastBorder;
	lut->bin[currentBin+1].rightB1= lastBorder;
      }
    }

    /* check right circle */
    if( xc <= patch->xMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	           floor( yc/patch->deltaY + patch->ySide*0.5 );
		   /*  rint( yc/patch->deltaY + patch->ySide*0.5 -0.5);*/
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0], patch );

	lut->bin[0].rightB2  = lastBorder;
	lut->bin[currentBin+1].leftB2 = lastBorder;
      }
    }

    /* real no intersection case */
    if( pieces ==0){
      *ifailP = 0;
      return;
    }

  }


  /************************************************/

  *lastBorderP =  lastBorder;
  *pathologyP = pathology;
  *directionP = direction;
  *rcOldP = rc;


  return;
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/****************************************************************/
/*								*/
/*     from followCircleCase.c          March 16, 2001         	*/
/*                                                         	*/
/*       program to find and identify bins and borders      	*/
/*	 in the circle case                                	*/
/*                                                        	*/
/*     Authors:  Sintes, A.M             			*/
/*								*/
/****************************************************************/


/****************************************************************/

static void FollowCircleCase(INT4 currentBin, UINT4  *lastBorderP, REAL8 alpha,
		      REAL8 ang1, REAL8  ang2, REAL8 rCritic,
		      REAL8 rcOld, INT4 *pathologyP,
		      INT4 *directionP, INT4 *ifailP, HOUGHptfLUT *lut,
		      HOUGHPatchGrid  *patch){

  INT4 lastBorder;
  INT4 pathology;    /* =1 (normal), =0 (anormal) */
  INT4 direction;    /* +1, or -1 */


  REAL8 rho1,rho2,radius;
  REAL8 xc,yc,rc; /* coordinates of the center of the circle */
  INT4 pieces;


  lastBorder = *lastBorderP;
  pathology = *pathologyP;
  direction = *directionP;

  rho1 = cos(ang1)/(1. -sin(ang1));
  rho2 = cos(ang2)/(1. -sin(ang2));
  rc   = rho1 + rho2;  /* positive or negative */
  xc   = rc*cos(alpha);
  yc   = rc*sin(alpha);
  radius = fabs( rho1 - rho2 ); /* abs could be avoided */

  /************************************************/
  /* check for exclusion of intersection of the circle with the patch*/
  if( ( yc+radius < patch->yMin ) || ( yc-radius > patch->yMax )  ||
      ( xc+radius < patch->xMin ) || ( xc-radius > patch->xMax )  ||
      ( sqrt(patch->xMax*patch->xMax + patch->yMax*patch->yMax) + radius < fabs(rc) ) ){
    /* no intersection */
    *ifailP = 0;  /* =1 (ok), =0 (stop) */
    return;
  }

  /************************************************/
  /* possible intersection case:                  */
  /************************************************/

  if( direction == -1){
    /* no pathological cases can happen */
    pieces = 0;

    /* check left circle */
    if( xc >= patch->xMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	           floor( yc/patch->deltaY + patch->ySide*0.5);
		   /*   rint( yc/patch->deltaY + patch->ySide*0.5 -0.5); */
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0] , patch);
	lut->bin[currentBin].rightB1  = lastBorder;
	lut->bin[currentBin+1].leftB1 = lastBorder;
      }
    }

    /* check right circle */
    if( xc <= patch->xMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	            floor( yc/patch->deltaY + patch->ySide*0.5 );
		    /*  rint( yc/patch->deltaY + patch->ySide*0.5 -0.5); */
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0], patch );
	lut->bin[currentBin].leftB2    = lastBorder;
	lut->bin[currentBin+1].rightB2 = lastBorder;
      }
    }

    /* real no intersection case */
    if( pieces ==0 ){
      *ifailP = 0;
      return;
    }

  } else {

    /************************************************/
    /*  This means ( direction == 1)                */
    /************************************************/
    /* pathologies can happen, i.e. circles can convert into lines
       and the two circles describing the borders of an annulus
       do not need to be concentric but:  ( ) ( ) */

    /* check for pathologies */
    if( (rcOld -rCritic)*(rc -rCritic) < 0 ){
      pathology = 0;
      direction = -1; /* for the next step */
    }

    pieces = 0;
    /* check left circle */
    if( xc >= patch->xMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	           floor( yc/patch->deltaY + patch->ySide*0.5 );
		   /*   rint( yc/patch->deltaY + patch->ySide*0.5 -0.5); */
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0] , patch);
	if( pathology ){
	  lut->bin[currentBin].leftB1   = lastBorder;
	  lut->bin[currentBin+1].rightB1= lastBorder;
	} else {
	  lut->bin[currentBin].rightB2  = lastBorder;
	  lut->bin[currentBin+1].leftB1 = lastBorder;
	}
      }
    }

    /* check right circle */
    if( xc <= patch->xMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn, patch);
      if(noIn){
	++pieces;
	++lastBorder;

	if (yymin < 0) {
	  fprintf(stderr,"WARNING: Fixing yymin (%d -> 0) [ConstructPLUT.c %d]\n",
		  yymin, __LINE__);
	  yymin = 0;
	}
	if (yymax >= patch->ySide) {
	  fprintf(stderr,"WARNING: Fixing yymax (%d -> %d) [ConstructPLUT.c %d]\n",
		  yymax, patch->ySide-1, __LINE__);
	  yymax = patch->ySide - 1;
	}

	lut->border[lastBorder].yUpper = yymax;
	lut->border[lastBorder].yLower = yymin;
	lut->border[lastBorder].yCenter =
	             floor( yc/patch->deltaY + patch->ySide*0.5 );
		     /*  rint( yc/patch->deltaY + patch->ySide*0.5 -0.5);*/
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &lut->border[lastBorder].xPixel[0] , patch);
	if( pathology ){
	  lut->bin[currentBin].rightB2  = lastBorder;
	  lut->bin[currentBin+1].leftB2 = lastBorder;
	} else {
	  lut->bin[currentBin].leftB1   = lastBorder;
	  lut->bin[currentBin+1].rightB2= lastBorder;
	}
      }
    }

    /* real no intersection case */
    if( pieces ==0){
      *ifailP = 0;
      return;
    }

  }

  /************************************************/

  *lastBorderP =  lastBorder;
  *pathologyP = pathology;
  *directionP = direction;

  return;
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/****************************************************************/
/*          from checkCircles.c        March 30, 2001          	*/
/*                                                         	*/
/*       subroutines to check if a circle intersects         	*/
/*	 the patch						*/
/*		                                        	*/
/*                                                        	*/
/*     Authors:  Sintes, A.M            			*/
/*								*/
/****************************************************************/


/*********************************************************/
/*  Subroutine to check if the left part of a circle     */
/*  intersects the patch (assuming xc>=patch->xMin)          */
/*  Output:                                              */
/*     information if intersects the patch:              */
/*               noIn =0 no intersection                  */
/* 		 noIn =1 intersection                    */
/*     and the y-pixel range of the intersection         */
/*********************************************************/

static void CheckLeftCircle(REAL8 xc, REAL8 yc, REAL8 radius,
		     INT4 *yyminP, INT4 *yymaxP, INT4 *noInP,
		     HOUGHPatchGrid  *patch){

  INT4 yymin,yymax;
  INT4 noIn, noIn1, noIn2;
  REAL8 ylower;
  REAL8 yupper,kkpos;
  volatile REAL4 kkk;


  /*********************************************************/
  noIn = 0;
  noIn1= 0; /* upper quadrant */
  noIn2= 0; /* lower quadrant */
  *noInP = noIn;

  if(xc < patch->xMin) return;  /* non optimized */

  /*********************************************************/
  /* Reduction of the interval to look at */
  /*********************************************************/

  ylower = MAX( patch->yMin, yc-radius);
  yupper = MIN( patch->yMax, yc+radius);

  /* convert to the values of the 'near' y-pixel */
  /* yymax  = floor((REAL4) (yupper/patch->deltaY -0.5) + (REAL4) (0.5*patch->ySide)); */
  kkk =  yupper/patch->deltaY -0.5;
  kkk += 0.5*patch->ySide;
  yymax  = floor(kkk);

  kkk = ylower/patch->deltaY - 0.5;
  kkk += 0.5*patch->ySide;
  yymin  = ceil(kkk);
  /* yymin  = ceil( (ylower/patch->deltaY-0.5) +patch->ySide/2.); */

  /* NEVER try here getting the pixel center, problems looking like
     horizontal lines */
/*
 *   ylower = patch->yCoor[yymin];
 *   yupper = patch->yCoor[yymax];
 *
 */


  /*********************************************************/
  /* looking at the upper-left quadrant */
  /*********************************************************/

  if( yc < patch->yMax ){
    REAL8 x1,myy1;

    /* setting the yymax value */
    /* in case of no over-flow problems !!! */

    kkpos= (radius-(yupper-yc) )*(radius+(yupper-yc) );
    kkpos= fabs(kkpos);

     x1 = xc - sqrt( kkpos);
     /* x1 = xc - sqrt( radius*radius -(yupper-yc)*(yupper-yc)); */

    if( (x1 >= patch->xMin) && (x1 <= patch->xMax) ){
      noIn1= 1; /* does intersect and yymax is the value to return! */
    }
    else{
      if( x1 > patch->xMax ){
        kkpos = (radius-(patch->xMax-xc))*(radius+(patch->xMax-xc));
	kkpos= fabs(kkpos);
	myy1 = yc + sqrt(kkpos);

	/* myy1 = yc + sqrt( radius*radius -(patch->xMax-xc)*(patch->xMax-xc));*/
	if( myy1>=patch->yMin ){
	  noIn1 = 1;
	  /* yymax = floor((REAL4) (myy1/patch->deltaY-0.5) + (REAL4) (0.5*patch->ySide)); */
	  kkk =  myy1/patch->deltaY-0.5;
	  kkk += 0.5*patch->ySide;
	  yymax = floor(kkk);
	}
      }
    }
  }
  /*********************************************************/
  /* looking at the lower-left quadrant */
  /*********************************************************/

  if( yc > patch->yMin ){
    REAL8 x1,myy1;

    /*setting the yymin value */
    kkpos= (radius-(yc-ylower) )*(radius+(yc-ylower) );
    kkpos= fabs(kkpos);
    x1 = xc - sqrt( kkpos);
    /* x1 = xc - sqrt( radius*radius -(yc-ylower)*(yc-ylower)); */

    if( (x1 >= patch->xMin) && (x1 <= patch->xMax) ){
      noIn2 = 1; /* does intersect and yymin is the value to return! */
    }
    else{
      if( x1 > patch->xMax ){
        kkpos = (radius-(patch->xMax-xc))*(radius+(patch->xMax-xc));
	kkpos= fabs(kkpos);
	myy1 = yc - sqrt(kkpos);

	/* myy1 = yc - sqrt( radius*radius -(patch->xMax-xc)*(patch->xMax-xc));*/

	if( myy1<=patch->yMax ){
	  noIn2 = 1;
	  /* yymin = ceil((REAL4) (myy1/patch->deltaY-0.5)+(REAL4)(0.5*patch->ySide)); */
	  kkk =  myy1/patch->deltaY-0.5;
	  kkk += 0.5*patch->ySide;
	  yymin = ceil(kkk);
	}
      }
    }
  }

  /*********************************************************/
  /* correcting yymax, yymin values */
  /*********************************************************/

  if ( noIn1 && (!noIn2) ){
    REAL8 x1,myy1;

    noIn = 1;
    kkpos= (radius-(yc-ylower) )*(radius+(yc-ylower) );
    kkpos= fabs(kkpos);
    x1 = xc - sqrt( kkpos);
    /* x1 = xc - sqrt( radius*radius -(yc-ylower)*(yc-ylower));*/

    if( x1 < patch->xMin ){
      /* the value yymin, needs to be set correctly */
      kkpos = (radius-(patch->xMin-xc))*(radius+(patch->xMin-xc));
      kkpos= fabs(kkpos);
      myy1 = yc + sqrt(kkpos);
      /* myy1 = yc + sqrt( radius*radius -(patch->xMin-xc)*(patch->xMin-xc));*/
      /* yymin = ceil( (REAL4) (myy1/patch->deltaY-0.5)+ (REAL4)(0.5*patch->ySide) ); */
      kkk =  myy1/patch->deltaY-0.5;
      kkk += 0.5*patch->ySide;
      yymin = ceil(kkk);
    }
  }

  /*********************************************************/

  if ( (! noIn1) && noIn2 ){
    /* the value yymax, needs to be set correctly */
    REAL8 x1,myy1;

    noIn = 1;
    kkpos= (radius-(yupper-yc) )*(radius+(yupper-yc) );
    kkpos= fabs(kkpos);
    x1 = xc - sqrt( kkpos);

    /* x1 = xc - sqrt( radius*radius -(yupper-yc)*(yupper-yc));*/
    if( x1 < patch->xMin ){
      /* the value yymin, needs to be set correctly */
      kkpos = (radius-(patch->xMin-xc))*(radius+(patch->xMin-xc));
      kkpos= fabs(kkpos);
      myy1 = yc - sqrt(kkpos);

      /*myy1 = yc - sqrt( radius*radius -(patch->xMin-xc)*(patch->xMin-xc));*/
      kkk =  myy1/patch->deltaY-0.5;
      kkk += 0.5*patch->ySide;
      yymax = floor(kkk);
    }
  }

  /* we could also set noIn = noIn1+noIn2; */
  /*********************************************************/
  if ( noIn1 && noIn2 ){
    /* the values yymin, yymax are correct */
    noIn = 1;
  }
  /*********************************************************/

  *yyminP = yymin;
  *yymaxP = yymax;
  *noInP  = noIn;

return;
}


/*********************************************************/
/*  Subroutine to check if the right part of a circle    */
/*  intersects the patch.                                */
/*  Output:                                              */
/*     information if intersects the patch:              */
/*               noIn =0 nointersection                  */
/* 		 noIn =1 intersection                    */
/*     and the y-pixel range of the intersection         */
/*********************************************************/


static void CheckRightCircle(REAL8 xc, REAL8 yc, REAL8 radius,
		     INT4 *yyminP, INT4 *yymaxP, INT4 *noInP,
		     HOUGHPatchGrid  *patch){

  INT4 yymin,yymax;
  INT4 noIn, noIn1, noIn2;
  REAL8 ylower;
  volatile REAL4 kkk;
  REAL8 yupper, kkpos;


 /*********************************************************/
  noIn = 0;
  noIn1= 0; /* upper quadrant */
  noIn2= 0; /* lower quadrant */
  *noInP = noIn;
  if( xc > patch->xMax ) return;  /* non optimized */

  /*********************************************************/
  /* Reduction of the interval to look at */
  /*********************************************************/
  ylower = MAX( patch->yMin, yc-radius);
  yupper = MIN( patch->yMax, yc+radius);

  /* convert to the value of the 'near' y-pixel */
  /* yymax = floor((REAL4) (yupper/patch->deltaY -0.5) + (REAL4) (0.5*patch->ySide)); */
  kkk =  yupper/patch->deltaY -0.5;
  kkk += 0.5*patch->ySide;
  yymax = floor(kkk);

  /* yymin  = ceil( (ylower/patch->deltaY-0.5) +patch->ySide/2.); */
  kkk = (ylower/patch->deltaY - 0.5);
  kkk += 0.5*patch->ySide;
  yymin  = ceil(kkk);


  /*********************************************************/
  /* looking at the upper-right quadrant */
  /*********************************************************/

  if( yc < patch->yMax ){
    REAL8 x1,myy1;

    /* in case of no over-flow problems !!! */

    kkpos= (radius-(yupper-yc) )*(radius+(yupper-yc) );
    kkpos= fabs(kkpos);
    x1 = xc + sqrt( kkpos);

    /* x1 = xc + sqrt( radius*radius -(yupper-yc)*(yupper-yc));*/

    if( (x1 >= patch->xMin) && (x1 <= patch->xMax) ){
      noIn1 = 1; /* does intersect and yymax is the value to return! */
    }
    else{
      if( x1 < patch->xMin ){
        kkpos = (radius-(patch->xMin-xc))*(radius+(patch->xMin-xc));
        kkpos= fabs(kkpos);
        myy1 = yc + sqrt(kkpos);

	/*myy1 = yc + sqrt(radius*radius -(patch->xMin-xc)*(patch->xMin-xc));*/
	if( myy1>=patch->yMin ){
	  noIn1  = 1;
	  /* yymax = floor(myy1/patch->deltaY+patch->ySide/2.-0.5);*/
	  /* yymax = floor((REAL4) (myy1/patch->deltaY-0.5) + (REAL4) (0.5*patch->ySide)); */
	  kkk =  myy1/patch->deltaY-0.5;
	  kkk += 0.5*patch->ySide;
	  yymax = floor(kkk);
	}
      }
    }
  }

 /*********************************************************/
  /* looking at the lower-right quadrant */
 /*********************************************************/

  if( yc > patch->yMin ){
    REAL8 x1,myy1;

    kkpos= (radius-(yc-ylower) )*(radius+(yc-ylower) );
    kkpos= fabs(kkpos);
    x1 = xc + sqrt( kkpos);

    /* x1 = xc + sqrt( radius*radius -(yc-ylower)*(yc-ylower));*/
    if( (x1 >= patch->xMin) && (x1 <= patch->xMax) ){
      noIn2 = 1; /* does intersect and yymin is the value to return! */
    }
    else{
      if( x1 < patch->xMin ){
        kkpos = (radius-(patch->xMin-xc))*(radius+(patch->xMin-xc));
        kkpos= fabs(kkpos);
        myy1 = yc - sqrt(kkpos);

	/* myy1 = yc - sqrt( radius*radius -(patch->xMin-xc)*(patch->xMin-xc));*/
	if( myy1<=patch->yMax ){
	  noIn2  = 1;
	  /* yymin = ceil(myy1/patch->deltaY+patch->ySide/2.-0.5);*/
	  /* yymin = ceil((REAL4) (myy1/patch->deltaY-0.5)+(REAL4)(0.5*patch->ySide)); */
	  kkk =  myy1/patch->deltaY-0.5;
	  kkk += 0.5*patch->ySide;
	  yymin = ceil(kkk);
	}
      }
    }
  }

  /*********************************************************/
  /* correcting yymax, yymin values */
  /*********************************************************/

  if ( noIn1 && (!noIn2) ){
    REAL8 x1,myy1;

    noIn = 1;
    kkpos= (radius-(yc-ylower) )*(radius+(yc-ylower) );
    kkpos= fabs(kkpos);
    x1 = xc + sqrt( kkpos);

    /* x1 = xc + sqrt( radius*radius -(yc-ylower)*(yc-ylower));*/
    if( x1 > patch->xMax ){
      /* the value yymin, needs to be set correctly */
      kkpos = (radius-(patch->xMax-xc))*(radius+(patch->xMax-xc));
      kkpos= fabs(kkpos);
      myy1 = yc + sqrt(kkpos);

      /* myy1 = yc + sqrt( radius*radius -(patch->xMax-xc)*(patch->xMax-xc));*/
      /* yymin = ceil( myy1/patch->deltaY+patch->ySide/2.-0.5 ); */
      /* yymin = ceil( (REAL4) (myy1/patch->deltaY-0.5)+ (REAL4)(0.5*patch->ySide) ); */
      kkk =  myy1/patch->deltaY-0.5;
      kkk += 0.5*patch->ySide;
      yymin = ceil(kkk);
    }
  }

  /*********************************************************/

  if ( (! noIn1) && noIn2 ){
    /* the value yymax, needs to be set correctly */
    REAL8 x1,myy1;

    noIn = 1;
    kkpos= (radius-(yupper-yc) )*(radius+(yupper-yc) );
    kkpos= fabs(kkpos);
    x1 = xc + sqrt( kkpos);
    /* x1 = xc + sqrt( radius*radius -(yupper-yc)*(yupper-yc)); */
    if( x1 > patch->xMax ){
      /* the value yymin, needs to be set correctly */
      kkpos = (radius-(patch->xMax-xc))*(radius+(patch->xMax-xc));
      kkpos= fabs(kkpos);
      myy1 = yc - sqrt(kkpos);

      /* myy1 = yc - sqrt( radius*radius -(patch->xMax-xc)*(patch->xMax-xc));*/
      /* yymax = floor( myy1/patch->deltaY+patch->ySide/2.-0.5 ); */
      /* yymax = floor( (REAL4) (myy1/patch->deltaY-0.5)+ (REAL4)(0.5*patch->ySide)); */
      kkk =  myy1/patch->deltaY-0.5;
      kkk += 0.5*patch->ySide;
      yymax = floor(kkk);

    }
  }

  /*********************************************************/

  if ( noIn1 && noIn2 ){
    /* the values yymin, yymax are correct */
    noIn = 1;
  }

  /*********************************************************/

  *yyminP = yymin;
  *yymaxP = yymax;
  *noInP  = noIn;

return;
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/****************************************************************/
/*      from drawCircles.c         Marc 16, 2001          	*/
/*                                                         	*/
/*       subroutines to clip circles on the patch         	*/
/*		                                        	*/
/*                                                        	*/
/*     Authors:  Sintes, A.M     		        	*/
/*								*/
/****************************************************************/


/****************************************************************/
/*   using sqrt!!!,                                             */
/*   to be substitute by the Bresenham algorithm                */
/****************************************************************/

static void DrawLeftCircle(REAL8 xc, REAL8 yc, REAL8 radius,
		   INT4 yymin, INT4 yymax, COORType *column,
		   HOUGHPatchGrid  *patch){
  INT4 jj;

  for(jj=yymin; jj<=yymax; ++jj){
    REAL8   realx, kkpos;
    INT4    myindex;
    volatile REAL4   kkk;

    kkpos = (radius-(yc-patch->yCoor[jj]) )*(radius+ (yc-patch->yCoor[jj]));
    kkpos= fabs(kkpos);
    realx = xc - sqrt(kkpos);

    /*realx = xc - sqrt( radius*radius
    	- (yc-patch->yCoor[jj])*(yc-patch->yCoor[jj]) ); */
    /* myindex = ceil((REAL4) (realx/patch->deltaX-0.5)+ (REAL4)(0.5*patch->xSide) ); */
    kkk =  realx/patch->deltaX-0.5;
    kkk += 0.5*patch->xSide;
    myindex = ceil(kkk);
    if( myindex<0 ) myindex=0;
    if (myindex > patch->xSide) myindex = patch->xSide;

    column[jj] = myindex;
  }

 return;
}


/****************************************************************/
/*   using sqrt!!!,                                             */
/*   to be substitute by the Bresenham algorithm                */
/****************************************************************/

static void DrawRightCircle(REAL8 xc, REAL8 yc, REAL8 radius,
		    INT4 yymin, INT4 yymax, COORType *column,
		    HOUGHPatchGrid  *patch){
  INT4 jj;

  for(jj=yymin; jj<=yymax; ++jj){
    REAL8  realx,kkpos;
    INT4   myindex;
    volatile REAL4   kkk;

    kkpos = (radius-(yc-patch->yCoor[jj]) )*(radius+ (yc-patch->yCoor[jj]));
    kkpos= fabs(kkpos);
    realx = xc + sqrt(kkpos);

    /*realx = xc + sqrt( radius*radius
    	- (yc-patch->yCoor[jj])*( yc-patch->yCoor[jj]) ); */
    /* myindex = ceil((REAL4)(realx/patch->deltaX-0.5)+(REAL4)(0.5*patch->xSide)); */
    kkk =  realx/patch->deltaX-0.5;
    kkk += 0.5*patch->xSide;
    myindex = ceil(kkk);
    if( myindex<0 ) myindex=0;
    if( myindex > patch->xSide-1 ) myindex=patch->xSide;
    column[jj] = myindex;
  }

 return;
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/****************************************************************/
/*     from  fill1Column.c            March 21, 2001     	*/
/*                                                         	*/
/*       Subroutines to correct border effects           	*/
/*								*/
/*                                                        	*/
/*     Authors:   Sintes, A.M                                   */
/*								*/
/****************************************************************/


/****************************************************************/
/* Normal case: for circles and lines. NOT the pathological one */
/****************************************************************/

static void Fill1Column(INT4 currentBin, UINT4 *lastBorderP,
                      HOUGHptfLUT *lut, HOUGHPatchGrid  *patch){

  INT4 lb1,rb1,lb2,rb2; /* The border index. If zero means that */
        /* it does not intersect the patch, or nothing to clip */
  INT4 lastBorder;

  lastBorder = *lastBorderP;

  lb1 = lut->bin[currentBin].leftB1;
  rb1 = lut->bin[currentBin].rightB1;
  lb2 = lut->bin[currentBin].leftB2;
  rb2 = lut->bin[currentBin].rightB2;

  /************************************************************/
  /* we want to set the values                                */
  /*              lut->bin[currentBin].piece*          */
  /* and if needed  set:                                      */
  /*      ++lastBorder;                                       */
  /*      lut->bin[currentBin].rightB1 = lastBorder;   */
  /*      lut->border[lastBorder].yUpper               */
  /*      lut->border[lastBorder].yLower               */
  /*      lut->border[lastBorder].xPixel[***] = 0;     */
  /************************************************************/


  /* call the proper case to correct the border effect */
  /* could be improved using nested if's */

  if( lb1 && rb1 ){
    FillCaseN1(lb1,rb1,currentBin, lut, patch);
    return;
  }

  if(lb1 && !(rb1) && !(lb2) ){
    FillCaseN2(lb1, currentBin, lut, patch);
    return;
  }

  if(lb1 && !(rb1) && lb2 ){
    FillCaseN3(lb1,lb2, currentBin, &lastBorder, lut, patch);
    *lastBorderP = lastBorder;
    return;
  }

  if(!(lb1) && rb1 && !(rb2) ){
    FillCaseN4(rb1, currentBin, lut, patch);
    return;
  }

  if(!(lb1) && rb1 && rb2){
    FillCaseN5(rb1,rb2, currentBin, lut);
    return;
  }

  if(!(lb1) && !(rb1) && lb2 && rb2){
    FillCaseN6(lb2,rb2, currentBin, lut, patch);
    return;
  }

  if(!(lb1) && !(rb1) && lb2 && !(rb2)){
    FillCaseN7(lb2, currentBin, lut, patch);
    return;
  }

  if(!(lb1) && !(rb1) && !(lb2) && rb2){
    FillCaseN8(rb2, currentBin, lut, patch);
    return;
  }

  if(!(currentBin) &&!(lb1) && !(rb1) && !(lb2) && !(rb2)){
    lut->bin[0].piece1max = patch->ySide-1;
    return;
  }


  return;
}


/***************************************************************/
/* case: ( lb1 && rb1 )  */
/***************************************************************/
static void FillCaseN1(INT4 lb1, INT4 rb1, INT4 currentBin,
                      HOUGHptfLUT *lut, HOUGHPatchGrid  *patch){
  INT4 lb1UpY,lb1LoY;
  INT4 rb1UpY,rb1LoY;
  INT4 yCl,yCr;

  lb1UpY = lut->border[lb1].yUpper;
  lb1LoY = lut->border[lb1].yLower;

  if( (lb1UpY == patch->ySide-1) && (lb1LoY == 0) ) return;

  rb1UpY = lut->border[rb1].yUpper;
  rb1LoY = lut->border[rb1].yLower;

  yCl = lut->border[lb1].yCenter;
  yCr = lut->border[rb1].yCenter;

  if( lb1LoY > yCl ){
    lut->bin[currentBin].piece1max = lb1LoY-1;
    if( rb1LoY > yCr ){
      lut->bin[currentBin].piece1min = rb1LoY;
    }
    /*     else{  */   /* already initialized */
    /*       lut->bin[currentBin].piece1min = 0; } */
    return;
  }

  if( lb1UpY < yCl ){
    lut->bin[currentBin].piece1min = lb1UpY+1;
    if( rb1UpY < yCr ){
      lut->bin[currentBin].piece1max = rb1UpY;
    }
    else{
      lut->bin[currentBin].piece1max = patch->ySide-1;
    }
   }

  return;
}

/***************************************************************/
/* case: (lb1 && !(rb1) && !(lb2) ) */
/***************************************************************/
static void FillCaseN2(INT4 lb1, INT4 currentBin,
                     HOUGHptfLUT *lut, HOUGHPatchGrid  *patch){
  INT4 lb1UpY,lb1LoY;
  INT4 yCl;

  lb1UpY = lut->border[lb1].yUpper;
  lb1LoY = lut->border[lb1].yLower;
  yCl = lut->border[lb1].yCenter;

  if( lb1LoY > yCl ){
      lut->bin[currentBin].piece1max = lb1LoY-1;
      /*  lut->bin[currentBin].piece1min = 0; */
      return;
  }

  if( lb1UpY < yCl ){
      lut->bin[currentBin].piece1max = patch->ySide-1;
      lut->bin[currentBin].piece1min = lb1UpY+1;
  }
  return;
}

/***************************************************************/
/* case: (lb1 && !(rb1) && lb2 ) */
/***************************************************************/
static void FillCaseN3(INT4 lb1, INT4 lb2, INT4 currentBin, INT4 *lastBorderP,
                     HOUGHptfLUT *lut, HOUGHPatchGrid  *patch){
  INT4 lastBorder;
  INT4 lb1UpY,lb1LoY;
  INT4 lb2UpY,lb2LoY;
  INT4 yCl;
  INT4 k;

  lastBorder = *lastBorderP;
  lb1UpY = lut->border[lb1].yUpper;
  lb1LoY = lut->border[lb1].yLower;
  yCl = lut->border[lb1].yCenter;
  lb2UpY = lut->border[lb2].yUpper;
  lb2LoY = lut->border[lb2].yLower;

  if( lb1LoY > yCl ){
    lut->bin[currentBin].piece1max = lb1LoY-1;
    lut->bin[currentBin].piece1min = lb2UpY+1;
    lut->bin[currentBin].piece2max = lb2LoY-1;
    /*  lut->bin[currentBin].piece2min = 0; */
    return;
  }

  if( lb1UpY < yCl ){
    lut->bin[currentBin].piece1max = patch->ySide-1;
    lut->bin[currentBin].piece1min = lb2UpY+1;
    lut->bin[currentBin].piece2max = lb2LoY-1;
    lut->bin[currentBin].piece2min = lb1UpY+1;
    return;
  }

  /* adding a ficticious border rb1 to correct the lb1 effects */
  ++lastBorder;
  lut->bin[currentBin].rightB1 = lastBorder;
  lut->border[lastBorder].yUpper = lb2UpY;
  lut->border[lastBorder].yLower = lb2LoY;
  for( k=lb2LoY; k<=lb2UpY;++k){
    lut->border[lastBorder].xPixel[k] = 0;
  }
  *lastBorderP = lastBorder;

  return;
}

/***************************************************************/
/* case: (!(lb1) && rb1 && !(rb2) ) */
/***************************************************************/
static void FillCaseN4(INT4 rb1, INT4 currentBin,
                     HOUGHptfLUT *lut, HOUGHPatchGrid  *patch){
  INT4 rb1UpY,rb1LoY;
  INT4 yCr;

  rb1UpY = lut->border[rb1].yUpper;
  rb1LoY = lut->border[rb1].yLower;
  yCr = lut->border[rb1].yCenter;

  if( rb1UpY < yCr ){
    lut->bin[currentBin].piece1max = rb1UpY;
    /*  lut->bin[currentBin].piece1min = 0; */
    return;
  }

  if( rb1LoY > yCr ){
    lut->bin[currentBin].piece1max = patch->ySide-1;
    lut->bin[currentBin].piece1min = rb1LoY;
    return;
  }

  lut->bin[currentBin].piece1max = patch->ySide-1;
  /*  lut->bin[currentBin].piece1min = 0; */

  return;
}

/***************************************************************/
/* case: (!(lb1) && rb1 && rb2 ) */
/***************************************************************/
static void  FillCaseN5(INT4 rb1, INT4 rb2, INT4 currentBin,
                     HOUGHptfLUT *lut){
  INT4 rb1UpY,rb1LoY;
  INT4 rb2UpY,rb2LoY;
  INT4 yCr;

  rb1UpY = lut->border[rb1].yUpper;
  rb1LoY = lut->border[rb1].yLower;
  yCr = lut->border[rb1].yCenter;
  rb2UpY = lut->border[rb2].yUpper;
  rb2LoY = lut->border[rb2].yLower;

  if( rb1UpY < yCr ){
    lut->bin[currentBin].piece1max = rb1UpY;
    lut->bin[currentBin].piece1min = rb2LoY;
    return;
  }

  if( rb1LoY > yCr ){
    lut->bin[currentBin].piece1max = rb2UpY;
    lut->bin[currentBin].piece1min = rb1LoY;
    return;
  }

  lut->bin[currentBin].piece1max = rb2UpY;
  lut->bin[currentBin].piece1min = rb2LoY;

  return;
}


/***************************************************************/
/* case: (!(lb1) && !(rb1) && lb2 && rb2) */
/***************************************************************/
static void  FillCaseN6(INT4 lb2, INT4 rb2, INT4 currentBin,
                     HOUGHptfLUT *lut, HOUGHPatchGrid  *patch){
  INT4 lb2UpY,lb2LoY;
  INT4 rb2UpY,rb2LoY;
  INT4 yCl,yCr;

  lb2UpY = lut->border[lb2].yUpper;
  lb2LoY = lut->border[lb2].yLower;

  if( (lb2UpY == patch->ySide-1) && (lb2LoY == 0) ) return;

  rb2UpY = lut->border[rb2].yUpper;
  rb2LoY = lut->border[rb2].yLower;
  yCl = lut->border[lb2].yCenter;
  yCr = lut->border[rb2].yCenter;

  if( lb2UpY >= yCl ){
    lut->bin[currentBin].piece1min = lb2UpY+1;
    if( rb2UpY >= yCr ){
      lut->bin[currentBin].piece1max = rb2UpY;
    }
    else{
      lut->bin[currentBin].piece1max = patch->ySide-1;
    }
  }

  if( lb2LoY <= yCl ){
    lut->bin[currentBin].piece2max = lb2LoY-1;
    if( rb2LoY <= yCr){
      lut->bin[currentBin].piece2min = rb2LoY;
    }
    /*   else{ */
    /*    lut->bin[currentBin].piece2min = 0; } */
  }

  return;
}

/***************************************************************/
/* case:  (!(lb1) && !(rb1) && lb2 && !(rb2)) */
/***************************************************************/
static void FillCaseN7(INT4 lb2, INT4 currentBin,
            HOUGHptfLUT *lut, HOUGHPatchGrid  *patch){
  INT4 lb2UpY,lb2LoY;
  INT4 yCl;

  lb2UpY = lut->border[lb2].yUpper;
  lb2LoY = lut->border[lb2].yLower;
  yCl = lut->border[lb2].yCenter;

  if( lb2UpY >= yCl ){
    lut->bin[currentBin].piece1max = patch->ySide-1;
    lut->bin[currentBin].piece1min = lb2UpY+1;
  }

  if( lb2LoY <= yCl ){
    lut->bin[currentBin].piece2max = lb2LoY-1;
    /* lut->bin[currentBin].piece2min = 0; */
  }

  return;
}

/***************************************************************/
/* case:  (!(lb1) && !(rb1) && !(lb2) && rb2) */
/***************************************************************/
static void  FillCaseN8(INT4 rb2, INT4 currentBin,
                      HOUGHptfLUT *lut, HOUGHPatchGrid  *patch){
  INT4 rb2UpY,rb2LoY;
  INT4 yCr;

  rb2UpY = lut->border[rb2].yUpper;
  rb2LoY = lut->border[rb2].yLower;
  yCr = lut->border[rb2].yCenter;

  lut->bin[currentBin].piece1max = patch->ySide-1;
  /* lut->bin[currentBin].piece1min = 0; */

  if( rb2UpY >= yCr){
    lut->bin[currentBin].piece1max = rb2UpY;
  }

  if( rb2LoY <= yCr){
    lut->bin[currentBin].piece1min = rb2LoY;
  }

  return;
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/****************************************************************/
/*      from fill1ColumnAnormal.c          March 23, 2001    	*/
/*                                                         	*/
/*       Subroutines to correct border effects           	*/
/*								*/
/*                                                        	*/
/*     Authors:   Sintes, A.M                                   */
/*								*/
/****************************************************************/


/****************************************************************/
/* Anormal case: for circles and lines. The pathological case */
/****************************************************************/

static void Fill1ColumnAnor(INT4 currentBin, HOUGHptfLUT *lut,
                          HOUGHPatchGrid  *patch){

  /* It will be similar to the non-pathological case.
     Note that here  the circles have very large radius.
     The y.coor of the center is very large or set to the border.
     Logic not defined yet! */


  INT4 lb1,rb1,lb2,rb2; /* The border index. If zero means that */
        /* it does not intersect the patch, or nothing to clip */

   /* The pathological case implies the name convention:
      lb1 is equivalent to lb2 ), rb2 equivalent to rb1 (,
      just to avoid border having the same name */


  lb1 = lut->bin[currentBin].leftB1;
  rb1 = lut->bin[currentBin].rightB1;
  lb2 = lut->bin[currentBin].leftB2;
  rb2 = lut->bin[currentBin].rightB2;

  /************************************************************/
  /* we want to set the values                                */
  /*              lut->bin[currentBin].piece*          */
  /************************************************************/



  /* call the proper case to correct the border effect */

 if( rb1 && !(rb2) && !(lb1) ){
    FillCaseN4(rb1,currentBin, lut, patch);
    return;
  }

 if( !(rb1)&& lb2 && !(rb2) && !(lb1) ){
    FillCaseN7(lb2,currentBin, lut, patch);
    return;
  }

 if( !(rb1)&& !(lb2) && rb2 ){
    FillCaseN4(rb2,currentBin, lut, patch);
    return;
  }

 if( !(rb1)&& !(lb2) && !(rb2) && lb1 ){
    FillCaseN7(lb1,currentBin, lut, patch);
    return;
  }
 /************************************************************/

 if( rb1 && rb2 ){
    FillCaseA1(rb1,rb2,currentBin, lut);
    return;
  }

 if( !(rb1)&& lb2 && !(rb2) && lb1 ){
    FillCaseA2(lb1,lb2,currentBin, lut);
    return;
  }

  if( !(rb2) && lb1 && rb1 ){
    FillCaseA3(lb1,rb1,currentBin, lut, patch);
    return;
  }

  if( !(rb1)&& lb2 && rb2 ){
    FillCaseA3(lb2,rb2,currentBin, lut, patch);
    return;
  }

  if(!(currentBin) && !(lb1) && !(rb1) && !(lb2) && !(rb2)){
    lut->bin[0].piece1max = patch->ySide-1;
    return;
  }

  return;
}

/***************************************************************/
/* case: ( rb1 && rb2 )  */
/***************************************************************/
static void FillCaseA1(INT4 rb1, INT4 rb2, INT4 currentBin,
                     HOUGHptfLUT *lut){
  INT4 rb1UpY,rb1LoY;
  INT4 rb2UpY,rb2LoY;
  INT4 yCr1;

  rb1UpY = lut->border[rb1].yUpper;
  rb1LoY = lut->border[rb1].yLower;
  yCr1 = lut->border[rb1].yCenter;

  rb2UpY = lut->border[rb2].yUpper;
  rb2LoY = lut->border[rb2].yLower;

  /* since the circles are big, number of case are simplifyed */
  /* more cases will appear in a general case with small radius */

  if( yCr1 < rb1UpY ){
    lut->bin[currentBin].piece1max = rb2UpY;
    lut->bin[currentBin].piece1min = rb1LoY;
  } else {
    lut->bin[currentBin].piece1max = rb1UpY;
    lut->bin[currentBin].piece1min = rb2LoY;
  }

  return;
}

/***************************************************************/
/* case: ( !(rb1)&& lb2 && !(rb2) && lb1 )  */
/***************************************************************/
static void FillCaseA2(INT4 lb1, INT4 lb2, INT4 currentBin,
                     HOUGHptfLUT *lut){
  INT4 lb1UpY,lb1LoY;
  INT4 lb2UpY,lb2LoY;
  INT4 yCl1;

  lb1UpY = lut->border[lb1].yUpper;
  lb1LoY = lut->border[lb1].yLower;
  yCl1 = lut->border[lb1].yCenter;

  lb2UpY = lut->border[lb2].yUpper;
  lb2LoY = lut->border[lb2].yLower;

  /* since the circles are big, number of case are simplifyed */
  /* more cases will appear in the general case */

  if( yCl1 < lb1UpY ){
    lut->bin[currentBin].piece1max = lb2LoY-1;
    lut->bin[currentBin].piece1min = lb1UpY+1;
  } else {
    lut->bin[currentBin].piece1max = lb1LoY-1;
    lut->bin[currentBin].piece1min = lb2UpY+1;
  }

  return;
}

/***************************************************************/
/* case: ( !(rb2) && lb1 && rb1 ) or  ( !(rb1)&& lb2 && rb2 ) */
/***************************************************************/
static void FillCaseA3(INT4 lb1, INT4 rb1, INT4 currentBin,
                     HOUGHptfLUT *lut, HOUGHPatchGrid  *patch){

  /* here we should code all possible cases with no exceptions */
  INT4 rb1UpY,rb1LoY;
  INT4 lb1UpY,lb1LoY;
  INT4 yCr1, yCl1;

  rb1UpY = lut->border[rb1].yUpper;
  rb1LoY = lut->border[rb1].yLower;
  yCr1 = lut->border[rb1].yCenter;

  lb1UpY = lut->border[lb1].yUpper;
  lb1LoY = lut->border[lb1].yLower;
  yCl1 = lut->border[lb1].yCenter;


  if( (lb1UpY == patch->ySide-1) && (lb1LoY == 0)) return;

  /* provisional settings */
  if( lb1UpY > yCl1 ){
    lut->bin[currentBin].piece1max = patch->ySide-1;
    lut->bin[currentBin].piece1min = lb1UpY+1;
  }

  if( lb1LoY < yCl1 ){
    lut->bin[currentBin].piece2max = lb1LoY-1;
    /* lut->bin[currentBin].piece2min = 0; */
  }

  /* corrections */
  if ( rb1LoY > yCr1 ){
    lut->bin[currentBin].piece2min = rb1LoY;
    return;
  }

  if ( rb1UpY < yCr1 ) {
    lut->bin[currentBin].piece1max = rb1UpY;
  }

  return;

}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */




#undef MIN
#undef MAX
