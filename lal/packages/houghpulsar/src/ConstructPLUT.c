/*-----------------------------------------------------------------------
 *
 * File Name: ConstructPLUT.c
 *
 * Authors: Sintes, A.M., 
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes June 7, 2001
 *
 *-----------------------------------------------------------------------
 */

/*
 * 1.  An author and Id block
 */

/************************************ <lalVerbatim file="ConstructPLUTCV">
Author: Sintes, A. M. 
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
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Static global variables}

The patch-time-frequency Look Up Table\\
\verb@static HOUGHptfLUT   *LookUpTable;@\\

 \noindent The patch limits with respect to the centers of the last pixels\\
\verb@static REAL8   LUTxMin;@\\
\verb@static REAL8   LUTxMax;@\\
\verb@static REAL8   LUTyMin;@\\
\verb@static REAL8   LUTyMax;@\\

 \noindent Coordiantes of the pixel centers, in the projected plane \\
\verb@static REAL8    *Xcoor;@\\
\verb@static REAL8    *Ycoor;@\\
\verb@static REAL8    DIFFX, DIFFY, INVDIFFX, INVDIFFY;@\\

\noindent \verb@static UINT2   LUTxSide, LUTySide;@


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Static function declarations}
\begin{verbatim}
static void PLUTInitialize(UINT2);
static void FillPLUT(HOUGHParamPLUT *,UINT2, UINT2);
static void CheckLeftCircle(REAL8, REAL8, REAL8, INT4 *, INT4 *, INT4 *);
static void CheckRightCircle(REAL8, REAL8, REAL8, INT4 *, INT4 *, INT4 *);
static void DrawRightCircle(REAL8, REAL8, REAL8, INT4, INT4, COORType *);
static void DrawLeftCircle(REAL8, REAL8, REAL8, INT4, INT4, COORType *);
static void CheckLineCase(REAL8, REAL8, REAL8, REAL8 *, INT4 *);
static void FindExactLine(REAL8, REAL8, REAL8 *, REAL8 *);
static void FindLine(REAL8, REAL8, REAL8, REAL8 *, REAL8 *);
static void CheckLineIntersection(REAL8, REAL8, REAL8, INT4 *, INT4 *, INT4 *);
static void DrawLine(REAL8, REAL8, REAL8, INT4, INT4, COORType *);
static void Fill1Column(INT4, INT4*);
static void FillCaseN1(INT4, INT4, INT4);
static void FillCaseN2(INT4, INT4);
static void FillCaseN3(INT4, INT4, INT4, INT4 *);
static void FillCaseN4(INT4, INT4);
static void FillCaseN5(INT4, INT4, INT4);
static void FillCaseN6(INT4, INT4, INT4);
static void FillCaseN7(INT4, INT4);
static void FillCaseN8(INT4, INT4);
static void Fill1ColumnAnor(INT4);
static void FillCaseA1(INT4, INT4, INT4);
static void FillCaseA2(INT4, INT4, INT4);
static void FillCaseA3(INT4, INT4, INT4);
static void InitialCircleCase(INT4 *,REAL8, REAL8, REAL8, REAL8, REAL8 *,INT4 *,INT4 *);
static void SecondCircleCase(INT4, INT4*, REAL8, REAL8, REAL8, REAL8, INT4,
                             REAL8*, INT4*, INT4 *, INT4*);
static void FollowCircleCase(INT4, INT4 *, REAL8, REAL8, REAL8, REAL8, REAL8,
                             INT4 *,INT4 *,INT4 *);
static void InitialLineCase(INT4 *, REAL8, REAL8, REAL8, INT4 *);
static void SecondLineCase(INT4, INT4 *, REAL8, REAL8, REAL8, INT4 *);
static void FollowLineCase(INT4, INT4 *, REAL8, REAL8, REAL8, REAL8, INT4, INT4 *);
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

#include <lal/LALConstants.h>
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


/* the patch-time-frequency Look Up Table */
static HOUGHptfLUT   *LookUpTable;

/* the patch limits with respect to the centers of the last pixels */
static REAL8   LUTxMin;
static REAL8   LUTxMax;
static REAL8   LUTyMin;
static REAL8   LUTyMax;

/* coordiantes of the pixel centers, in the projected plane */
static REAL8    *Xcoor;
static REAL8    *Ycoor;

static REAL8    DIFFX, DIFFY, INVDIFFX, INVDIFFY;

/* to be changed */
static UINT2   LUTxSide, LUTySide;


/*
 * 5.f) Static function declarations
 */

static void PLUTInitialize(UINT2);
static void FillPLUT(HOUGHParamPLUT *, UINT2, UINT2);

static void CheckLeftCircle(REAL8,REAL8,REAL8,INT4 *,INT4 *,INT4 *);
static void CheckRightCircle(REAL8,REAL8,REAL8,INT4 *,INT4 *,INT4 *);
static void DrawRightCircle(REAL8,REAL8,REAL8,INT4,INT4, COORType *);
static void DrawLeftCircle(REAL8,REAL8,REAL8,INT4,INT4,COORType *);
static void CheckLineCase(REAL8, REAL8, REAL8,REAL8 *, INT4 *);
static void FindExactLine(REAL8,REAL8,REAL8 *,REAL8 *);
static void FindLine(REAL8,REAL8,REAL8,REAL8 *,REAL8 *);
static void CheckLineIntersection(REAL8,REAL8,REAL8,INT4 *,INT4 *,INT4 *);
static void DrawLine(REAL8, REAL8, REAL8,INT4, INT4, COORType *);
static void Fill1Column(INT4, INT4*);
static void FillCaseN1(INT4, INT4, INT4);
static void FillCaseN2(INT4, INT4);
static void FillCaseN3(INT4, INT4, INT4, INT4 *);
static void FillCaseN4(INT4, INT4);
static void FillCaseN5(INT4, INT4, INT4);
static void FillCaseN6(INT4, INT4, INT4);
static void FillCaseN7(INT4, INT4);
static void FillCaseN8(INT4, INT4);
static void Fill1ColumnAnor(INT4);
static void FillCaseA1(INT4, INT4, INT4);
static void FillCaseA2(INT4, INT4, INT4);
static void FillCaseA3(INT4, INT4, INT4);

static void InitialCircleCase(INT4 *,REAL8, REAL8, REAL8, REAL8, REAL8 *, INT4 *, INT4 *);
static void SecondCircleCase(INT4, INT4*,REAL8,REAL8, REAL8,REAL8,INT4,REAL8*,INT4*,INT4 *,INT4*);
static void FollowCircleCase(INT4,INT4 *,REAL8,REAL8,REAL8,REAL8,REAL8,INT4 *,INT4 *,INT4 *);
static void InitialLineCase(INT4 *, REAL8, REAL8, REAL8, INT4 *);
static void SecondLineCase(INT4, INT4 *, REAL8, REAL8, REAL8, INT4 *);
static void FollowLineCase(INT4, INT4 *,REAL8, REAL8, REAL8,REAL8,INT4,INT4 *);


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
  INT2    i;
  /* --------------------------------------------- */
  INITSTATUS (status, " LALHOUGHConstructPLUT", CONSTRUCTPLUTC);
  /*  ATTATCHSTATUSPTR (status); */

  /*   Make sure the arguments are not NULL: */ 
  ASSERT (lut, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (patch, status, LUTH_ENULL, LUTH_MSGENULL);
  ASSERT (par,  status, LUTH_ENULL, LUTH_MSGENULL);
  
  ASSERT (par->deltaF == patch->deltaF,  status, LUTH_EVAL, LUTH_MSGEVAL);
  /* -------------------------------------------   */

  f0Bin =  par->f0Bin;

 /* -------------------------------------------   */
 /* setting Global variables & fields  */

   LookUpTable = lut;
   
 /*LookUpTable->timeIndex = ; */
  LookUpTable->deltaF = par->deltaF;
  LookUpTable->f0Bin  = f0Bin;
  LookUpTable->nFreqValid = PIXERR * f0Bin *VEPI/VTOT;
  
  LUTxMin = patch->xMin;
  LUTxMax = patch->xMax;
  LUTyMin = patch->yMin;
  LUTyMax = patch->yMax;

  Xcoor = patch->xCoor; 
  Ycoor = patch->yCoor;

  DIFFX = patch->deltaX;
  DIFFY = patch->deltaY;
  INVDIFFX = 1./ DIFFX;
  INVDIFFY = 1./ DIFFY;

  LUTxSide = patch->xSide;
  LUTySide = patch->ySide;

  /* -------------------------------------------   */
  
  PLUTInitialize(lut->maxNBins);  
  FillPLUT(par, lut->maxNBins, lut->maxNBorders);
  
  /* *lut= LookUpTable; old!!!  */
  
  /* Make sure number of bins makes sense with the dimensions  */
  ASSERT (lut->nBin >0, status, LUTH_ESIZE, LUTH_MSGESIZE);
  ASSERT (lut->nBin <= lut->maxNBins, status, LUTH_ESIZE, LUTH_MSGESIZE);
  
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

static void  PLUTInitialize(UINT2 maxNBins){
  UINT2  i;
  for(i=0;i<maxNBins;++i){
    LookUpTable->bin[i].leftB1  = 0;
    LookUpTable->bin[i].rightB1 = 0;
    LookUpTable->bin[i].leftB2  = 0;
    LookUpTable->bin[i].rightB2 = 0;
    LookUpTable->bin[i].piece1max = -1; 
    LookUpTable->bin[i].piece1min = 0;
    LookUpTable->bin[i].piece2max = -1;
    LookUpTable->bin[i].piece2min = 0;
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
/* a pixel j , corresponds to x(j)= DIFFX*(0.5+j- LUTxSide/2.) */
/* the nearest center of a pixel is j = round[ x/DIFFX +LUTxSide/2.-0.5]  */
/* ******************************************************************* */
/* The way to convert k into the bin index is the following:    */
/*    for k>=0, binindex= k                                     */
/*    for k<0,  binindex= nBin+iniBin-1-k                       */
/*   f0 corresponds to k=0                                      */
/* This will be used when reading the peakgram!                 */ 
/* ******************************************************************* */

static void  FillPLUT(HOUGHParamPLUT  *par, UINT2 maxNBins, UINT2 maxNBorders){

  /********************************************************/
  /*    variables that need to be calculated before       */
  /********************************************************/

  REAL8 cosDelta;   /* = df/|xi|                         */
  REAL8 cosPhiMax0; /* = (xi*N +df/2)/|xi|               */
  REAL8 cosPhiMin0; /* = cosPhiMax0-deltaCos             */
  REAL8 alpha;      /* = xi.alpha in the rotated coordinates */
  REAL8 delta;      /* = xi.delta                        */
  REAL8 epsilon;    /* = LINERR *8.0d-8 * f0/df          */

  /********************************************************/

  INT4 lastBorder =0;  /* counter of the last build border */ 
  INT4 currentBin =0;  /* counter of the bin studied       */

  INT4 ifailPlus  = 1; /* =1 (ok, continue to next bin), =0 (stop) */
  INT4 ifailMinus = 1; /* =1 (ok, continue to previous bin), =0 (stop) */

  INT4 directionPlus;  /* = +1,or -1 */
  INT4 directionMinus; /* = +1,or -1 */

  INT4 pathology;   /* =1 (normal), =0 (anormal case) */
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

  alpha = par->xi.alpha;
  delta = par->xi.delta;
  cosDelta = par->cosDelta;
  cosPhiMax0 = par->cosPhiMax0;
  cosPhiMin0 = par->cosPhiMin0; 
  epsilon = par->epsilon;
  
  /********************************************************/

  lambda = 2* delta -LAL_PI*0.5;
  rCritic = 2* cos(lambda) /(1 - sin(lambda) );

  /********************************************************/
  /*  starting with the (central) bin corresponding to:   */
  /*            Delta_f(t) =0 (k=0), border cosPhiMax     */
  /********************************************************/ 

  /* cosPhiMax = MIN(1, cosPhiMax0); */
  cosPhiMax = cosPhiMax0;
  
  if(cosPhiMax > 0.99999999){ /* avoid points or small circles */
    ifailPlus = 0;            /* do not go to the next bin */
  } else{
    
    phi  = acos(cosPhiMax);   /* in the interval (0,PI) */
    ang1 = delta + phi;
    ang2 = delta - phi;
    
    /* check for lines, or numerical lines! */
    CheckLineCase(epsilon, ang1, ang2, &eps, &lineCase);

    if( lineCase ){
 	/* line case */
      InitialLineCase(&lastBorder,alpha,delta,eps, &ifailPlus);
      directionPlus = -1;
    } else{
      /* circle case */
      InitialCircleCase(&lastBorder,alpha, ang1, ang2, rCritic, 
			&rcOldPlus, &directionPlus, &ifailPlus);
    }    
  }
  

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
#ifdef CHECKHOUGHINDEX
    if (currentBin > maxNBins || lastBorder>= maxNBorders ){ 
      fprintf(stderr,"currentBin=%d not in range 1 to maxNBins=%d\n"
	      "or lastborder=%d >= maxNBorders=%d\n",
	      currentBin,maxNBins,lastBorder,maxNBorders);
      abort();
    }
#endif

    cosPhiMax = cosPhiMax + cosDelta;
    /* or cosPhiMax = MIN(1,cosPhiMax + cosDelta ); */

    if( cosPhiMax > 0.999999){ /* check appropiate value */
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
		       rcOldPlus, directionPlus, &ifailPlus);
      } else{
	/* circle case */
	FollowCircleCase(currentBin, &lastBorder,alpha, ang1, ang2,rCritic, 
			 rcOldPlus, &pathology, &directionPlus, &ifailPlus);
      }    
    }
    
    if(pathology){
      Fill1Column(currentBin, &lastBorder);
    }else{
      Fill1ColumnAnor(currentBin); /* not yet coded */
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
      SecondLineCase(currentBin,&lastBorder,alpha,delta,eps, &ifailMinus);
      directionMinus = -1;
      pathology = 0;  
    } else{
      /* circle case */
      pathology = 1; /* provisionally */
      SecondCircleCase(currentBin, &lastBorder,alpha, ang1, ang2, 
		       rCritic, directionPlus, &rcOldMinus, 
		       &pathology, &directionMinus, &ifailMinus);
    }    
  }

  /********************************************************/
  /*  the way to identify initial pathologies:            */
  /*    initial or second case being a line !             */
  /*    or two circles, both with direction=-1            */
  /********************************************************/

  if(pathology){
    Fill1Column(0, &lastBorder);
  }else{
    Fill1ColumnAnor(0); /* not yet coded */
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
#ifdef CHECKHOUGHINDEX
    if (currentBin > maxNBins || lastBorder>= maxNBorders ){ 
      fprintf(stderr,"currentBin=%d not in range 1 to maxNBins=%d\n"
	      "or lastborder=%d >= maxNBorders=%d\n",
	      currentBin,maxNBins,lastBorder,maxNBorders);
      abort();
    }
#endif
    

    cosPhiMin = cosPhiMin - cosDelta;

    if( cosPhiMin < -0.999999){ /* check appropiate value */
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
		       rcOldMinus, directionMinus, &ifailMinus);
      } else{
	/* circle case */
	FollowCircleCase(currentBin, &lastBorder,alpha, ang1, ang2,rCritic, 
			 rcOldMinus, &pathology, &directionMinus, &ifailMinus);
      }    
    }
    
    if(pathology){
      Fill1Column(currentBin, &lastBorder);
    }else{
      Fill1ColumnAnor(currentBin); /* not yet coded */
    }
  }
  
  /********************************************************/
  /* set iniBin,nBin  etc */
  /********************************************************/

  LookUpTable->nBin = currentBin+1;
  LookUpTable->iniBin = nBinPos-currentBin;

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

static void InitialLineCase(INT4  *lastBorderP, REAL8 alpha, REAL8  delta, 
		     REAL8 eps, INT4 *ifailP){


  INT4 lastBorder;
  INT4 ifail;        /* =1 (ok), =0 (stop) */

  REAL8  xA,yA;
  INT4 yymin,yymax;
  REAL8 xRel, slope;
  INT4 noIn;   /* if no intersection occurs noIn=0 */

  lastBorder = *lastBorderP;

 /************************************************/
  FindLine(alpha, delta, eps, &xA, &yA);
  CheckLineIntersection(alpha, xA, yA, &yymin, &yymax, &noIn);

  if( noIn ==0 ){
    *ifailP = 0;
    return;
  }
  ++lastBorder;

  LookUpTable->border[lastBorder].yUpper = yymax;
  LookUpTable->border[lastBorder].yLower = yymin;

  DrawLine(alpha, xA, yA, yymin, yymax,
	   &LookUpTable->border[lastBorder].xPixel[0] );

  /************************************************/

  if( (alpha == LAL_PI*0.5) || (alpha == LAL_PI*1.5) ||
      (alpha == -LAL_PI*0.5) ){       /* horizontal line */

    if( yA < 0 ){ /* convention */
      LookUpTable->bin[0].rightB1 = lastBorder;
      LookUpTable->bin[1].leftB1  = lastBorder;
      LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
    } else {
      LookUpTable->bin[0].leftB2  = lastBorder;
      LookUpTable->bin[1].rightB2 = lastBorder;
      LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
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
      LookUpTable->bin[0].leftB2 = lastBorder;
      LookUpTable->bin[1].rightB2= lastBorder;

      if ( slope < 0 ){
	LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
      }

    } else {
      LookUpTable->bin[0].rightB1 = lastBorder;
      LookUpTable->bin[1].leftB1  = lastBorder;
      
      if ( slope > 0 ){
	LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
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

static void SecondLineCase(INT4 currentBin, INT4  *lastBorderP, 
		    REAL8 alpha, REAL8  delta, 
		     REAL8 eps, INT4 *ifailP){

  /* we are for sure in a pathological case. Border names are
     changed accordingly */

  INT4 lastBorder;
  INT4 ifail;        /* =1 (ok), =0 (stop) */

  REAL8  xA,yA;
  INT4 yymin,yymax;
  REAL8 xRel, slope;
  INT4 noIn;   /* if no intersection occurs noIn=0 */

  lastBorder = *lastBorderP;

 /************************************************/
  FindLine(alpha, delta, eps, &xA, &yA);
  CheckLineIntersection(alpha, xA, yA, &yymin, &yymax, &noIn);

  if( noIn ==0 ){
    *ifailP = 0;
    return;
  }
  ++lastBorder;

  LookUpTable->border[lastBorder].yUpper = yymax;
  LookUpTable->border[lastBorder].yLower = yymin;

  DrawLine(alpha, xA, yA, yymin, yymax,
	   &LookUpTable->border[lastBorder].xPixel[0] );

  /************************************************/
  /* all are pathological cases. The code is similar to the 
     InitialLineCase, except for the modifications:
     rightB1 -> rightB2, leftB2 -> leftB1, where marked */

  if( (alpha == LAL_PI*0.5) || (alpha == LAL_PI*1.5) ||
      (alpha == -LAL_PI*0.5) ){       /* horizontal line */

    if( yA < 0 ){ /* convention */
      LookUpTable->bin[0].rightB2= lastBorder; /* modified */
      LookUpTable->bin[currentBin+1].leftB1 = lastBorder;
      LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
    } else {
      LookUpTable->bin[0].leftB1 = lastBorder; /* modified */
      LookUpTable->bin[currentBin+1].rightB2= lastBorder;
      LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
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
      LookUpTable->bin[0].leftB1   = lastBorder; /* modified */
      LookUpTable->bin[currentBin+1].rightB2= lastBorder;

      if ( slope < 0 ){
	LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
      }

    } else {
      LookUpTable->bin[0].rightB2  = lastBorder; /* modified */
      LookUpTable->bin[currentBin+1].leftB1 = lastBorder;
      
      if ( slope > 0 ){
	LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
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

static void FollowLineCase(INT4 currentBin, INT4  *lastBorderP,
		    REAL8 alpha, REAL8  delta, REAL8 eps,REAL8 rcOld, 
		    INT4 direction, INT4 *ifailP){

  INT4 lastBorder;
  INT4 ifail;        /* =1 (ok), =0 (stop) */

  REAL8  xA,yA;
  INT4 yymin,yymax;
  INT4 noIn;   /* if no intersection occurs noIn=0 */

  lastBorder = *lastBorderP;

  /************************************************/
  FindLine(alpha, delta, eps, &xA, &yA);
  CheckLineIntersection(alpha, xA, yA, &yymin, &yymax, &noIn);

  if( noIn ==0 ){
    *ifailP = 0;
    return;
  }
  ++lastBorder;

  LookUpTable->border[lastBorder].yUpper = yymax;
  LookUpTable->border[lastBorder].yLower = yymin;

  DrawLine(alpha, xA, yA, yymin, yymax,
	   &LookUpTable->border[lastBorder].xPixel[0] );

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

      LookUpTable->bin[currentBin].leftB1   = lastBorder;
      LookUpTable->bin[currentBin+1].rightB1= lastBorder;

      /* alternatively, one can also set two extra borders */
      /*  LookUpTable->bin[currentBin].rightB2  = lastBorder; */
      /*  LookUpTable->bin[currentBin+1].leftB2 = lastBorder; */

      if( ycOld < yA ){
	LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
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
	LookUpTable->bin[currentBin].leftB1   = lastBorder;
	LookUpTable->bin[currentBin+1].rightB1= lastBorder;

	if ( slope > 0 ){
	  LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
	} else {
	  LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
	}

      } else {
	LookUpTable->bin[currentBin].rightB2  = lastBorder;
	LookUpTable->bin[currentBin+1].leftB2 = lastBorder;

	if ( slope < 0 ){
	  LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
	} else {
	  LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
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
	LookUpTable->bin[currentBin].rightB1  = lastBorder;
	LookUpTable->bin[currentBin+1].leftB1 = lastBorder;
	LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
      } else {
	LookUpTable->bin[currentBin].leftB2   = lastBorder;
	LookUpTable->bin[currentBin+1].rightB2= lastBorder;
	LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
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
	LookUpTable->bin[currentBin].leftB2   = lastBorder;
	LookUpTable->bin[currentBin+1].rightB2= lastBorder;

	if ( slope < 0 ){
	  LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
	} else {
	  LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
	}

      } else {
	LookUpTable->bin[currentBin].rightB1  = lastBorder;
	LookUpTable->bin[currentBin+1].leftB1 = lastBorder;

	if ( slope > 0 ){
	  LookUpTable->border[lastBorder].yCenter = 0; /*or smaller*/
	} else {
	  LookUpTable->border[lastBorder].yCenter = LUTySide -1;/*or bigger */
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

static void FindExactLine(REAL8 alpha, REAL8 delta, 
		   REAL8 *xA, REAL8 *yA){

  REAL8 lambda, rA;
  
  lambda =  2.*delta -LAL_PI*0.5;
  rA = 2.* cos(lambda)/(1.-sin(lambda) ); 
  *xA = rA* cos(alpha);
  *yA = rA* sin(alpha);

return;
}

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
		      INT4 *yyminP, INT4 *yymaxP, INT4 *noInP){
 
  INT4 yymin,yymax,noIn;

  noIn = 0;

  if ((alpha == 0) || (alpha == LAL_PI) || (alpha == -LAL_PI) ){ 
    /* vertical line */
    if ((LUTxMin <= xA) && (LUTxMax >= xA)){
      noIn  = 1;
      yymin = 0;
      yymax = LUTySide-1;
    }
  }
  else{
    if ( (alpha == LAL_PI*0.5) || (alpha == LAL_PI*1.5) ||
	 (alpha == -LAL_PI*0.5) ){
      /* horizontal line */
      if ((LUTyMin <= yA) && (LUTyMax >= yA)){
	noIn  = 1;
	yymin =  ceil(yA/DIFFY+LUTySide*0.5-0.5);
	yymax = floor(yA/DIFFY+LUTySide*0.5-0.5);
	/* Note  yymax < yymin,   to identify an horizontal line!
	   If the area to mark is above (arriba), mark yymin and higher. 
	   If the area to mark is below (abajo), mark yymax and lower.*/
      }
    }
    else{ 
	/* generic case */
      REAL8 y1,y2,slope;
      slope = cot(alpha);
      y1  = slope*(xA-LUTxMin)+yA;
      y2  = slope*(xA-LUTxMax)+yA;
      
      if (  ( (y1 >= LUTyMin) || (y2 >= LUTyMin) ) 
	    && ( (y1 <= LUTyMax) || (y2 <= LUTyMax) ) ){ 
	REAL8 yupper,ylower;
	noIn   = 1;
	yupper = MAX(y1,y2);
	ylower = MIN(y1,y2); /* or  ylower=y1+y2-yupper  */
	yupper = MIN(yupper,LUTyMax);
	ylower = MAX(ylower,LUTyMin);
	yymin  = ceil(ylower/DIFFY+LUTySide*0.5-0.5);
	yymax  =floor(yupper/DIFFY+LUTySide*0.5-0.5);
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
	      INT4 yymin, INT4 yymax, COORType  *column){

  INT4 jj;
  
  column[yymin] = LUTxSide;
  column[yymax] = LUTxSide;

  /* the if-else, is not really needed, just to avoid repeating
     the ceil when not necessary! */
     
 if ((alpha == 0) || (alpha == LAL_PI)){ 
   /* vertical line */
   INT4  xpixel;

   xpixel = ceil(xA/DIFFX+LUTxSide*0.5-0.5);
   for(jj=yymin;jj<=yymax;++jj){
     column[jj] = xpixel;
   }
 }
 else{
   /* remaining cases */
   REAL8 tanalpha;
   REAL8 xofy;

   tanalpha = tan(alpha);

   for(jj=yymin;jj<=yymax;++jj){  
     /* will not be executed in the horizontal case */    
     xofy = xA + tanalpha*( yA- Ycoor[jj] );
     column[jj] = ceil(xofy/DIFFX+LUTxSide*0.5-0.5);
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

static void InitialCircleCase(INT4  *lastBorderP, REAL8 alpha,
		      REAL8 ang1, REAL8  ang2, REAL8 rCritic, 
		      REAL8 *rcOldP, INT4 *directionP, INT4 *ifailP){

  INT4 lastBorder;
  INT4 direction;    /* +1, or -1 */
  INT4 ifail;        /* =1 (ok), =0 (stop) */

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
  if( ( yc+radius < LUTyMin ) || ( yc-radius > LUTyMax )  ||
      ( xc+radius < LUTxMin ) || ( xc-radius > LUTxMax )  ||
      ( sqrt(LUTxMax*LUTxMax + LUTyMax*LUTyMax) + radius < fabs(rc) ) ){
    /* no intersection */
    *ifailP = 0;
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
    if( xc >= LUTxMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	LookUpTable->bin[0].rightB1= lastBorder;
	LookUpTable->bin[1].leftB1 = lastBorder;
      }
    }
    
    /* check right circle */
    if( xc <= LUTxMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	LookUpTable->bin[0].leftB2  = lastBorder;
	LookUpTable->bin[1].rightB2 = lastBorder;
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
    if( xc >= LUTxMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
       
	LookUpTable->bin[0].leftB1 = lastBorder;
	LookUpTable->bin[1].rightB1= lastBorder;
      }
    }
    
    /* check right circle */
    if( xc <= LUTxMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	
	LookUpTable->bin[0].rightB2= lastBorder;
	LookUpTable->bin[1].leftB2 = lastBorder;
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

static void SecondCircleCase(INT4 currentBin, INT4  *lastBorderP, 
		      REAL8 alpha, REAL8 ang1, REAL8  ang2, 
		      REAL8 rCritic, INT4 directionPlus, REAL8 *rcOldP, 
		      INT4 *pathologyP,INT4 *directionP, INT4 *ifailP){

  INT4 lastBorder;
  INT4 pathology;    /* =1 (normal), =0 (anormal) */
  INT4 direction;    /* +1, or -1 */
  INT4 ifail;        /* =1 (ok), =0 (stop) */

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
  if( ( yc+radius < LUTyMin ) || ( yc-radius > LUTyMax )  ||
      ( xc+radius < LUTxMin ) || ( xc-radius > LUTxMax )  ||
      ( sqrt(LUTxMax*LUTxMax + LUTyMax*LUTyMax) + radius < fabs(rc) ) ){
    /* no intersection */
    *ifailP = 0;
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
    if( xc >= LUTxMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	if ( pathology ){
	  LookUpTable->bin[0].rightB1 = lastBorder;
	} else{
	  LookUpTable->bin[0].rightB2 = lastBorder; /* modified */
	}
	LookUpTable->bin[currentBin+1].leftB1 = lastBorder;
      }
    }
    
    /* check right circle */
    if( xc <= LUTxMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	if ( pathology ){
	  LookUpTable->bin[0].leftB2 = lastBorder;
	} else{
	  LookUpTable->bin[0].leftB1 = lastBorder; /*modified */
	}
	LookUpTable->bin[currentBin+1].rightB2 = lastBorder;
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
    if( xc >= LUTxMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
       
	LookUpTable->bin[0].leftB1   = lastBorder;
	LookUpTable->bin[currentBin+1].rightB1= lastBorder;
      }
    }
    
    /* check right circle */
    if( xc <= LUTxMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	
	LookUpTable->bin[0].rightB2  = lastBorder;
	LookUpTable->bin[currentBin+1].leftB2 = lastBorder;
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

static void FollowCircleCase(INT4 currentBin, INT4  *lastBorderP, REAL8 alpha,
		      REAL8 ang1, REAL8  ang2, REAL8 rCritic, 
		      REAL8 rcOld, INT4 *pathologyP, 
		      INT4 *directionP, INT4 *ifailP){

  INT4 lastBorder;
  INT4 pathology;    /* =1 (normal), =0 (anormal) */
  INT4 direction;    /* +1, or -1 */
  INT4 ifail;        /* =1 (ok), =0 (stop) */

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
  if( ( yc+radius < LUTyMin ) || ( yc-radius > LUTyMax )  ||
      ( xc+radius < LUTxMin ) || ( xc-radius > LUTxMax )  ||
      ( sqrt(LUTxMax*LUTxMax + LUTyMax*LUTyMax) + radius < fabs(rc) ) ){
    /* no intersection */
    *ifailP = 0;
    return;
  }

  /************************************************/
  /* possible intersection case:                  */
  /************************************************/

  if( direction == -1){ 
    /* no pathological cases can happen */
    pieces = 0;

    /* check left circle */
    if( xc >= LUTxMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	LookUpTable->bin[currentBin].rightB1  = lastBorder;
	LookUpTable->bin[currentBin+1].leftB1 = lastBorder;
      }
    }
    
    /* check right circle */
    if( xc <= LUTxMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	LookUpTable->bin[currentBin].leftB2    = lastBorder;
	LookUpTable->bin[currentBin+1].rightB2 = lastBorder;
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
    if( xc >= LUTxMin ){ /* draw ( ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckLeftCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawLeftCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	if( pathology ){
	  LookUpTable->bin[currentBin].leftB1   = lastBorder;
	  LookUpTable->bin[currentBin+1].rightB1= lastBorder;
	} else {
	  LookUpTable->bin[currentBin].rightB2  = lastBorder;
	  LookUpTable->bin[currentBin+1].leftB1 = lastBorder;
	}
      }
    }
    
    /* check right circle */
    if( xc <= LUTxMax){ /* draw ) ? */
      INT4 yymax;
      INT4 yymin;
      INT4 noIn; /* if no intersection occurs noIn=0 */

      CheckRightCircle(xc,yc,radius, &yymin, &yymax, &noIn);
      if(noIn){
	++pieces;
	++lastBorder;

	LookUpTable->border[lastBorder].yUpper = yymax;
	LookUpTable->border[lastBorder].yLower = yymin;
	LookUpTable->border[lastBorder].yCenter = 
	             rint( yc/DIFFY + LUTySide*0.5 -0.5);
	DrawRightCircle(xc,yc,radius,yymin,yymax,
		       &LookUpTable->border[lastBorder].xPixel[0] );
	if( pathology ){
	  LookUpTable->bin[currentBin].rightB2  = lastBorder;
	  LookUpTable->bin[currentBin+1].leftB2 = lastBorder;
	} else {
	  LookUpTable->bin[currentBin].leftB1   = lastBorder;
	  LookUpTable->bin[currentBin+1].rightB2= lastBorder;
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
/*  intersects the patch (assuming xc>=LUTxMin)          */
/*  Output:                                              */
/*     information if intersects the patch:              */
/*               noIn =0 no intersection                  */
/* 		 noIn =1 intersection                    */
/*     and the y-pixel range of the intersection         */
/*********************************************************/

static void CheckLeftCircle(REAL8 xc, REAL8 yc, REAL8 radius,
		     INT4 *yyminP, INT4 *yymaxP, INT4 *noInP){
  
  INT4 yymin,yymax;
  INT4 noIn, noIn1, noIn2;
  REAL8 ylower;
  REAL8 yupper;
 
  if(xc < LUTxMin) return;  /* non optimized */
  
  /*********************************************************/
  noIn = 0;
  noIn1= 0; /* upper quadrant */
  noIn2= 0; /* lower quadrant */

  /*********************************************************/
  /* Reduction of the interval to look at */
  /*********************************************************/

  ylower = MAX( LUTyMin, yc-radius);
  yupper = MIN( LUTyMax, yc+radius);

 /* convert  to the values of the 'near' y-pixel */  
  yymax  = floor(yupper/DIFFY+LUTySide/2.-0.5);
  /*  yupper = Ycoor[yymax];  /*  yupper = DIFFY*(0.5 + yymax - LUTySide/2.) */
  /* this gives problems when circles are almost horizontal lines */

  yymin  = ceil(ylower/DIFFY+LUTySide/2.-0.5);
  /*  ylower = Ycoor[yymin];   /*ylower = DIFFY*(0.5 + yymin - LUTySide/2.) */

  /*********************************************************/
  /* looking at the upper-left quadrant */
  /*********************************************************/

  if( yc < LUTyMax ){
    REAL8 x1,y1;
    
    /* setting the yymax value */
    /* in case of no over-flow problems !!! */
    x1 = xc - sqrt( radius*radius -(yupper-yc)*(yupper-yc));

    if( (x1 >= LUTxMin) && (x1 <= LUTxMax) ){
      noIn1= 1; /* does intersect and yymax is the value to return! */
    }
    else{
      if( x1 > LUTxMax ){
	y1 = yc + sqrt( radius*radius -(LUTxMax-xc)*(LUTxMax-xc));
	if( y1>=LUTyMin ){
	  noIn1 = 1;
	  yymax = floor(y1/DIFFY+LUTySide/2.-0.5);
	}
      }
    }
  }
  /*********************************************************/
  /* looking at the lower-left quadrant */
  /*********************************************************/

  if( yc > LUTyMin ){
    REAL8 x1,y1;

    /*setting the yymin value */
    x1 = xc - sqrt( radius*radius -(yc-ylower)*(yc-ylower));
    if( (x1 >= LUTxMin) && (x1 <= LUTxMax) ){
      noIn2 = 1; /* does intersect and yymin is the value to return! */
    }
    else{
      if( x1 > LUTxMax ){
	y1 = yc - sqrt( radius*radius -(LUTxMax-xc)*(LUTxMax-xc));
	if( y1<=LUTyMax ){
	  noIn2 = 1;
	  yymin = ceil(y1/DIFFY+LUTySide/2.-0.5);
	}
      }
    }
  }

  /*********************************************************/
  /* correcting yymax, yymin values */
  /*********************************************************/

  if ( noIn1 && (!noIn2) ){ 
    REAL8 x1,y1;
    
    noIn = 1;
    x1 = xc - sqrt( radius*radius -(yc-ylower)*(yc-ylower));
    if( x1 < LUTxMin ){    
      /* the value yymin, needs to be set correctly */
      y1 = yc + sqrt( radius*radius -(LUTxMin-xc)*(LUTxMin-xc));
      yymin = ceil( y1/DIFFY+LUTySide/2.-0.5 );
    }
  }

  /*********************************************************/

  if ( (! noIn1) && noIn2 ){ 
    /* the value yymax, needs to be set correctly */
    REAL8 x1,y1;
 
   noIn = 1;
    x1 = xc - sqrt( radius*radius -(yupper-yc)*(yupper-yc));
    if( x1 < LUTxMin ){    
      /* the value yymin, needs to be set correctly */
      y1 = yc - sqrt( radius*radius -(LUTxMin-xc)*(LUTxMin-xc));
      yymax = floor( y1/DIFFY+LUTySide/2.-0.5 );
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
		     INT4 *yyminP, INT4 *yymaxP, INT4 *noInP){

  INT4 yymin,yymax;
  INT4 noIn, noIn1, noIn2;
  REAL8 ylower;
  REAL8 yupper;
 
  if( xc > LUTxMax ) return;  /* non optimized */

 /*********************************************************/
  noIn = 0;
  noIn1= 0; /* upper quadrant */
  noIn2= 0; /* lower quadrant */

  /*********************************************************/
  /* Reduction of the interval to look at */
  /*********************************************************/
  ylower = MAX( LUTyMin, yc-radius);
  yupper = MIN( LUTyMax, yc+radius);

  /* convert to the value of the 'near' y-pixel */ 
  yymax  = floor(yupper/DIFFY+LUTySide/2.-0.5);
  /* yupper = Ycoor[yymax];  /* yupper = DIFFY*(0.5 + yymax - LUTySide/2.) */
  
  yymin  = ceil(ylower/DIFFY+LUTySide/2.-0.5);
  /* ylower = Ycoor[yymin];  /* ylower = DIFFY*(0.5 + yymin - LUTySide/2.) */
  
  /*********************************************************/
  /* looking at the upper-right quadrant */
  /*********************************************************/
  
  if( yc < LUTyMax ){
    REAL8 x1,y1;

    /* in case of no over-flow problems !!! */
    x1 = xc + sqrt( radius*radius -(yupper-yc)*(yupper-yc));

    if( (x1 >= LUTxMin) && (x1 <= LUTxMax) ){
      noIn1 = 1; /* does intersect and yymax is the value to return! */
    }
    else{
      if( x1 < LUTxMin ){
	y1 = yc + sqrt(radius*radius -(LUTxMin-xc)*(LUTxMin-xc));
	if( y1>=LUTyMin ){
	  noIn1  = 1;
	  yymax = floor(y1/DIFFY+LUTySide/2.-0.5);
	}
      }
    }
  }

 /*********************************************************/
    /* looking at the lower-right quadrant */
 /*********************************************************/

  if( yc > LUTyMin ){
    REAL8 x1,y1;

    x1 = xc + sqrt( radius*radius -(yc-ylower)*(yc-ylower));
    if( (x1 >= LUTxMin) && (x1 <= LUTxMax) ){
      noIn2 = 1; /* does intersect and yymin is the value to return! */
    }
    else{
      if( x1 < LUTxMin ){
	y1 = yc - sqrt( radius*radius -(LUTxMin-xc)*(LUTxMin-xc));
	if( y1<=LUTyMax ){
	  noIn2  = 1;
	  yymin = ceil(y1/DIFFY+LUTySide/2.-0.5);
	}
      }
    }
  }

  /*********************************************************/
  /* correcting yymax, yymin values */
  /*********************************************************/
  
  if ( noIn1 && (!noIn2) ){ 
    REAL8 x1,y1;
    
    noIn = 1;
    x1 = xc + sqrt( radius*radius -(yc-ylower)*(yc-ylower));
    if( x1 > LUTxMax ){    
      /* the value yymin, needs to be set correctly */
      y1 = yc + sqrt( radius*radius -(LUTxMax-xc)*(LUTxMax-xc));
      yymin = ceil( y1/DIFFY+LUTySide/2.-0.5 );
    }
  }

  /*********************************************************/

  if ( (! noIn1) && noIn2 ){ 
    /* the value yymax, needs to be set correctly */
    REAL8 x1,y1;
 
   noIn = 1;
    x1 = xc + sqrt( radius*radius -(yupper-yc)*(yupper-yc));
    if( x1 > LUTxMax ){    
      /* the value yymin, needs to be set correctly */
      y1 = yc - sqrt( radius*radius -(LUTxMax-xc)*(LUTxMax-xc));
      yymax = floor( y1/DIFFY+LUTySide/2.-0.5 );
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
		   INT4 yymin, INT4 yymax, COORType *column){
  INT4 jj;

  for(jj=yymin; jj<=yymax; ++jj){
    REAL8   realx;
    INT4    index;
    realx = xc - sqrt( radius*radius - (yc-Ycoor[jj])*( yc-Ycoor[jj]) );
    index = ceil(realx*INVDIFFX+LUTxSide/2. -0.5);
    if( index<0 ) index=0;
    column[jj] = index;
  }
 
 return;
}


/****************************************************************/
/*   using sqrt!!!,                                             */
/*   to be substitute by the Bresenham algorithm                */
/****************************************************************/

static void DrawRightCircle(REAL8 xc, REAL8 yc, REAL8 radius,
		    INT4 yymin, INT4 yymax, COORType *column){
  INT4 jj;

  for(jj=yymin; jj<=yymax; ++jj){
    REAL8  realx;
    INT4   index;
    realx = xc + sqrt( radius*radius - (yc-Ycoor[jj])*( yc-Ycoor[jj]) );
    index = ceil(realx*INVDIFFX+LUTxSide/2. -0.5);
    if( index > LUTxSide-1 ) index=LUTxSide;
    column[jj] = index;
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

static void Fill1Column(INT4 currentBin, INT4 *lastBorderP){

  INT4 lb1,rb1,lb2,rb2; /* The border index. If zero means that */
        /* it does not intersect the patch, or nothing to clip */
  INT4 lastBorder;

  lastBorder = *lastBorderP;

  lb1 = LookUpTable->bin[currentBin].leftB1;
  rb1 = LookUpTable->bin[currentBin].rightB1;
  lb2 = LookUpTable->bin[currentBin].leftB2;
  rb2 = LookUpTable->bin[currentBin].rightB2;

  /************************************************************/
  /* we want to set the values                                */
  /*              LookUpTable->bin[currentBin].piece*          */
  /* and if needed  set:                                      */
  /*      ++lastBorder;                                       */
  /*      LookUpTable->bin[currentBin].rightB1 = lastBorder;   */
  /*      LookUpTable->border[lastBorder].yUpper               */
  /*      LookUpTable->border[lastBorder].yLower               */
  /*      LookUpTable->border[lastBorder].xPixel[***] = 0;     */
  /************************************************************/


  /* call the proper case to correct the border effect */
  /* could be improved using nested if's */

  if( lb1 && rb1 ){
    FillCaseN1(lb1,rb1,currentBin);
    return;
  }

  if(lb1 && !(rb1) && !(lb2) ){
    FillCaseN2(lb1, currentBin);
    return;
  }

  if(lb1 && !(rb1) && lb2 ){
    FillCaseN3(lb1,lb2, currentBin, &lastBorder);
    *lastBorderP = lastBorder;
    return;
  }

  if(!(lb1) && rb1 && !(rb2) ){
    FillCaseN4(rb1, currentBin);
    return;
  }

  if(!(lb1) && rb1 && rb2){
    FillCaseN5(rb1,rb2, currentBin);
    return;
  }

  if(!(lb1) && !(rb1) && lb2 && rb2){
    FillCaseN6(lb2,rb2, currentBin);
    return;
  }

  if(!(lb1) && !(rb1) && lb2 && !(rb2)){
    FillCaseN7(lb2, currentBin);
    return;
  }

  if(!(lb1) && !(rb1) && !(lb2) && rb2){
    FillCaseN8(rb2, currentBin);
    return;
  }

  if(!(currentBin) &&!(lb1) && !(rb1) && !(lb2) && !(rb2)){
    LookUpTable->bin[0].piece1max = LUTySide-1;
    return;
  }


  return;
}


/***************************************************************/
/* case: ( lb1 && rb1 )  */
/***************************************************************/
static void FillCaseN1(INT4 lb1, INT4 rb1, INT4 currentBin){
  INT4 lb1UpY,lb1LoY;
  INT4 rb1UpY,rb1LoY;
  INT4 yCl,yCr;

  lb1UpY = LookUpTable->border[lb1].yUpper;
  lb1LoY = LookUpTable->border[lb1].yLower;

  if( (lb1UpY == LUTySide-1) && (lb1LoY == 0) ) return;

  rb1UpY = LookUpTable->border[rb1].yUpper;
  rb1LoY = LookUpTable->border[rb1].yLower;
 
  yCl = LookUpTable->border[lb1].yCenter;
  yCr = LookUpTable->border[rb1].yCenter;

  if( lb1LoY > yCl ){
    LookUpTable->bin[currentBin].piece1max = lb1LoY-1;
    if( rb1LoY > yCr ){
      LookUpTable->bin[currentBin].piece1min = rb1LoY;
    }
    /*     else{  */   /* already initialized */
    /*       LookUpTable->bin[currentBin].piece1min = 0; } */
    return;
  }

  if( lb1UpY < yCl ){
    LookUpTable->bin[currentBin].piece1min = lb1UpY+1;
    if( rb1UpY < yCr ){
      LookUpTable->bin[currentBin].piece1max = rb1UpY;
    }
    else{ 
      LookUpTable->bin[currentBin].piece1max = LUTySide-1;
    }
   }

  return;
}

/***************************************************************/
/* case: (lb1 && !(rb1) && !(lb2) ) */
/***************************************************************/
static void FillCaseN2(INT4 lb1, INT4 currentBin){
  INT4 lb1UpY,lb1LoY;
  INT4 yCl;

  lb1UpY = LookUpTable->border[lb1].yUpper;
  lb1LoY = LookUpTable->border[lb1].yLower;
  yCl = LookUpTable->border[lb1].yCenter;

  if( lb1LoY > yCl ){
      LookUpTable->bin[currentBin].piece1max = lb1LoY-1;
      /*  LookUpTable->bin[currentBin].piece1min = 0; */
      return;
  }

  if( lb1UpY < yCl ){
      LookUpTable->bin[currentBin].piece1max = LUTySide-1;
      LookUpTable->bin[currentBin].piece1min = lb1UpY+1;
  }
  return;
}

/***************************************************************/
/* case: (lb1 && !(rb1) && lb2 ) */
/***************************************************************/
static void FillCaseN3(INT4 lb1, INT4 lb2, INT4 currentBin, INT4 *lastBorderP){
  INT4 lastBorder;
  INT4 lb1UpY,lb1LoY;
  INT4 lb2UpY,lb2LoY;
  INT4 yCl;
  INT4 k;

  lastBorder = *lastBorderP;
  lb1UpY = LookUpTable->border[lb1].yUpper;
  lb1LoY = LookUpTable->border[lb1].yLower;
  yCl = LookUpTable->border[lb1].yCenter;
  lb2UpY = LookUpTable->border[lb2].yUpper;
  lb2LoY = LookUpTable->border[lb2].yLower;

  if( lb1LoY > yCl ){
    LookUpTable->bin[currentBin].piece1max = lb1LoY-1;
    LookUpTable->bin[currentBin].piece1min = lb2UpY+1;
    LookUpTable->bin[currentBin].piece2max = lb2LoY-1;
    /*  LookUpTable->bin[currentBin].piece2min = 0; */
    return;
  } 

  if( lb1UpY < yCl ){
    LookUpTable->bin[currentBin].piece1max = LUTySide-1;
    LookUpTable->bin[currentBin].piece1min = lb2UpY+1;
    LookUpTable->bin[currentBin].piece2max = lb2LoY-1;
    LookUpTable->bin[currentBin].piece2min = lb1UpY+1; 
    return;
  }
  
  /* adding a ficticious border rb1 to correct the lb1 effects */
  ++lastBorder;                                      
  LookUpTable->bin[currentBin].rightB1 = lastBorder;  
  LookUpTable->border[lastBorder].yUpper = lb2UpY;              
  LookUpTable->border[lastBorder].yLower = lb2LoY;  
  for( k=lb2LoY; k<=lb2UpY;++k){          
    LookUpTable->border[lastBorder].xPixel[k] = 0;    
  }
  *lastBorderP = lastBorder;

  return;
}

/***************************************************************/
/* case: (!(lb1) && rb1 && !(rb2) ) */
/***************************************************************/
static void FillCaseN4(INT4 rb1, INT4 currentBin){
  INT4 rb1UpY,rb1LoY;
  INT4 yCr;
  
  rb1UpY = LookUpTable->border[rb1].yUpper;
  rb1LoY = LookUpTable->border[rb1].yLower;
  yCr = LookUpTable->border[rb1].yCenter;
  
  if( rb1UpY < yCr ){
    LookUpTable->bin[currentBin].piece1max = rb1UpY;
    /*  LookUpTable->bin[currentBin].piece1min = 0; */
    return;
  }
  
  if( rb1LoY > yCr ){
    LookUpTable->bin[currentBin].piece1max = LUTySide-1;
    LookUpTable->bin[currentBin].piece1min = rb1LoY;
    return;
  } 

  LookUpTable->bin[currentBin].piece1max = LUTySide-1;
  /*  LookUpTable->bin[currentBin].piece1min = 0; */

  return;
}

/***************************************************************/
/* case: (!(lb1) && rb1 && rb2 ) */
/***************************************************************/
static void  FillCaseN5(INT4 rb1, INT4 rb2, INT4 currentBin){
  INT4 rb1UpY,rb1LoY;
  INT4 rb2UpY,rb2LoY;
  INT4 yCr;
  
  rb1UpY = LookUpTable->border[rb1].yUpper;
  rb1LoY = LookUpTable->border[rb1].yLower;
  yCr = LookUpTable->border[rb1].yCenter;
  rb2UpY = LookUpTable->border[rb2].yUpper;
  rb2LoY = LookUpTable->border[rb2].yLower;

  if( rb1UpY < yCr ){
    LookUpTable->bin[currentBin].piece1max = rb1UpY;
    LookUpTable->bin[currentBin].piece1min = rb2LoY; 
    return;
  }
  
  if( rb1LoY > yCr ){
    LookUpTable->bin[currentBin].piece1max = rb2UpY;
    LookUpTable->bin[currentBin].piece1min = rb1LoY;
    return;
  } 

  LookUpTable->bin[currentBin].piece1max = rb2UpY;
  LookUpTable->bin[currentBin].piece1min = rb2LoY;

  return;
}


/***************************************************************/
/* case: (!(lb1) && !(rb1) && lb2 && rb2) */
/***************************************************************/
static void  FillCaseN6(INT4 lb2, INT4 rb2, INT4 currentBin){
  INT4 lb2UpY,lb2LoY;
  INT4 rb2UpY,rb2LoY;
  INT4 yCl,yCr;

  lb2UpY = LookUpTable->border[lb2].yUpper;
  lb2LoY = LookUpTable->border[lb2].yLower;

  if( (lb2UpY == LUTySide-1) && (lb2LoY == 0) ) return;

  rb2UpY = LookUpTable->border[rb2].yUpper;
  rb2LoY = LookUpTable->border[rb2].yLower;
  yCl = LookUpTable->border[lb2].yCenter;
  yCr = LookUpTable->border[rb2].yCenter;

  if( lb2UpY >= yCl ){
    LookUpTable->bin[currentBin].piece1min = lb2UpY+1;
    if( rb2UpY >= yCr ){
      LookUpTable->bin[currentBin].piece1max = rb2UpY;
    }
    else{
      LookUpTable->bin[currentBin].piece1max = LUTySide-1;
    }
  }

  if( lb2LoY <= yCl ){
    LookUpTable->bin[currentBin].piece2max = lb2LoY-1;
    if( rb2LoY <= yCr){
      LookUpTable->bin[currentBin].piece2min = rb2LoY;
    }
    /*   else{ */
    /*    LookUpTable->bin[currentBin].piece2min = 0; } */
  }
  
  return;
}

/***************************************************************/
/* case:  (!(lb1) && !(rb1) && lb2 && !(rb2)) */
/***************************************************************/
static void FillCaseN7(INT4 lb2, INT4 currentBin){
  INT4 lb2UpY,lb2LoY;
  INT4 yCl;

  lb2UpY = LookUpTable->border[lb2].yUpper;
  lb2LoY = LookUpTable->border[lb2].yLower;
  yCl = LookUpTable->border[lb2].yCenter;

  if( lb2UpY >= yCl ){
    LookUpTable->bin[currentBin].piece1max = LUTySide-1;
    LookUpTable->bin[currentBin].piece1min = lb2UpY+1;
  }

  if( lb2LoY <= yCl ){
    LookUpTable->bin[currentBin].piece2max = lb2LoY-1;
    /* LookUpTable->bin[currentBin].piece2min = 0; */
  }

  return;
}

/***************************************************************/
/* case:  (!(lb1) && !(rb1) && !(lb2) && rb2) */
/***************************************************************/
static void  FillCaseN8(INT4 rb2, INT4 currentBin){
  INT4 rb2UpY,rb2LoY;
  INT4 yCr;

  rb2UpY = LookUpTable->border[rb2].yUpper;
  rb2LoY = LookUpTable->border[rb2].yLower;
  yCr = LookUpTable->border[rb2].yCenter;

  LookUpTable->bin[currentBin].piece1max = LUTySide-1;
  /* LookUpTable->bin[currentBin].piece1min = 0; */

  if( rb2UpY >= yCr){
    LookUpTable->bin[currentBin].piece1max = rb2UpY;
  }

  if( rb2LoY <= yCr){
    LookUpTable->bin[currentBin].piece1min = rb2LoY;
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

static void Fill1ColumnAnor(INT4 currentBin){

  /* It will be similar to the non-pathological case.
     Note that here  the circles have very large radius.
     The y.coor of the center is very large or set to the border. 
     Logic not defined yet! */


  INT4 lb1,rb1,lb2,rb2; /* The border index. If zero means that */
        /* it does not intersect the patch, or nothing to clip */

   /* The pathological case implies the name convention:
      lb1 is equivalent to lb2 ), rb2 equivalent to rb1 (,
      just to avoid border having the same name */


  lb1 = LookUpTable->bin[currentBin].leftB1;
  rb1 = LookUpTable->bin[currentBin].rightB1;
  lb2 = LookUpTable->bin[currentBin].leftB2;
  rb2 = LookUpTable->bin[currentBin].rightB2;

  /************************************************************/
  /* we want to set the values                                */
  /*              LookUpTable->bin[currentBin].piece*          */
  /************************************************************/



  /* call the proper case to correct the border effect */
 
 if( rb1 && !(rb2) && !(lb1) ){
    FillCaseN4(rb1,currentBin);
    return;
  }

 if( !(rb1)&& lb2 && !(rb2) && !(lb1) ){
    FillCaseN7(lb2,currentBin);
    return;
  }
 
 if( !(rb1)&& !(lb2) && rb2 ){
    FillCaseN4(rb2,currentBin);
    return;
  }

 if( !(rb1)&& !(lb2) && !(rb2) && lb1 ){
    FillCaseN7(lb1,currentBin);
    return;
  }
 /************************************************************/
   
 if( rb1 && rb2 ){
    FillCaseA1(rb1,rb2,currentBin);
    return;
  }

 if( !(rb1)&& lb2 && !(rb2) && lb1 ){
    FillCaseA2(lb1,lb2,currentBin);
    return;
  }
 
  if( !(rb2) && lb1 && rb1 ){
    FillCaseA3(lb1,rb1,currentBin);
    return;
  }

  if( !(rb1)&& lb2 && rb2 ){
    FillCaseA3(lb2,rb2,currentBin);
    return;
  }
  
  if(!(currentBin) && !(lb1) && !(rb1) && !(lb2) && !(rb2)){
    LookUpTable->bin[0].piece1max = LUTySide-1;
    return;
  }
  
  return;
}

/***************************************************************/
/* case: ( rb1 && rb2 )  */
/***************************************************************/
static void FillCaseA1(INT4 rb1, INT4 rb2, INT4 currentBin){
  INT4 rb1UpY,rb1LoY;
  INT4 rb2UpY,rb2LoY;
  INT4 yCr1;

  rb1UpY = LookUpTable->border[rb1].yUpper;
  rb1LoY = LookUpTable->border[rb1].yLower;
  yCr1 = LookUpTable->border[rb1].yCenter;

  rb2UpY = LookUpTable->border[rb2].yUpper;
  rb2LoY = LookUpTable->border[rb2].yLower;

  /* since the circles are big, number of case are simplifyed */
  /* more cases will appear in a general case with small radius */

  if( yCr1 < rb1UpY ){
    LookUpTable->bin[currentBin].piece1max = rb2UpY;
    LookUpTable->bin[currentBin].piece1min = rb1LoY;
  } else {
    LookUpTable->bin[currentBin].piece1max = rb1UpY;
    LookUpTable->bin[currentBin].piece1min = rb2LoY;
  }
  
  return;
}

/***************************************************************/
/* case: ( !(rb1)&& lb2 && !(rb2) && lb1 )  */
/***************************************************************/
static void FillCaseA2(INT4 lb1, INT4 lb2, INT4 currentBin){
  INT4 lb1UpY,lb1LoY;
  INT4 lb2UpY,lb2LoY;
  INT4 yCl1;

  lb1UpY = LookUpTable->border[lb1].yUpper;
  lb1LoY = LookUpTable->border[lb1].yLower;
  yCl1 = LookUpTable->border[lb1].yCenter;

  lb2UpY = LookUpTable->border[lb2].yUpper;
  lb2LoY = LookUpTable->border[lb2].yLower;

  /* since the circles are big, number of case are simplifyed */
  /* more cases will appear in the general case */

  if( yCl1 < lb1UpY ){
    LookUpTable->bin[currentBin].piece1max = lb2LoY-1;
    LookUpTable->bin[currentBin].piece1min = lb1UpY+1;
  } else {
    LookUpTable->bin[currentBin].piece1max = lb1LoY-1;
    LookUpTable->bin[currentBin].piece1min = lb2UpY+1;
  }

  return;
}

/***************************************************************/
/* case: ( !(rb2) && lb1 && rb1 ) or  ( !(rb1)&& lb2 && rb2 ) */
/***************************************************************/
static void FillCaseA3(INT4 lb1, INT4 rb1, INT4 currentBin){

  /* here we should code all possible cases with no exceptions */
  INT4 rb1UpY,rb1LoY;
  INT4 lb1UpY,lb1LoY;
  INT4 yCr1, yCl1;

  rb1UpY = LookUpTable->border[rb1].yUpper;
  rb1LoY = LookUpTable->border[rb1].yLower;
  yCr1 = LookUpTable->border[rb1].yCenter;

  lb1UpY = LookUpTable->border[lb1].yUpper;
  lb1LoY = LookUpTable->border[lb1].yLower;
  yCl1 = LookUpTable->border[lb1].yCenter;


  if( (lb1UpY == LUTySide-1) && (lb1LoY == 0)) return;

  /* provisional settings */
  if( lb1UpY > yCl1 ){
    LookUpTable->bin[currentBin].piece1max = LUTySide-1;
    LookUpTable->bin[currentBin].piece1min = lb1UpY+1;
  }

  if( lb1LoY < yCl1 ){
    LookUpTable->bin[currentBin].piece2max = lb1LoY-1;
    /* LookUpTable->bin[currentBin].piece2min = 0; */
  }

  /* corrections */
  if ( rb1LoY > yCr1 ){
    LookUpTable->bin[currentBin].piece2min = rb1LoY;
    return;
  }

  if ( rb1UpY < yCr1 ) {
    LookUpTable->bin[currentBin].piece1max = rb1UpY;
  }

  return;

}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */




#undef MIN
#undef MAX
