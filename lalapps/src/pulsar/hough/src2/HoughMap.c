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

/**
 *
 * \file HoughMap.c
 * \brief Subroutines for initialization and construction of Hough-map derivatives and total Hough-maps.
 * \author Sintes, A.M., 
 * \ingroup pulsarHough
 * Revision: $Id$
 *
 * History:   Created by Sintes June 22, 2001
 *            Modified    "    August 6, 2001
 *
 
\par Description 

The function  LALHOUGHInitializeHD() initializes the Hough map derivative
space  HOUGHMapDeriv *hd to zero. Note that the length of the map
hd->map should be hd->ySide * (hd->xSide + 1).

The function  LALHOUGHInitializeHT() initializes the total Hough map 
HOUGHMapTotal *ht  to zero and checks consistency between 
the number of physical pixels in the
map  and  those given by the grid information structure 
HOUGHPatchGrid *patch.

Given an initial Hough map derivative HOUGHMapDeriv *hd and a representation
of a phmd HOUGHphmd *phmd, the function  LALHOUGHAddPHMD2HD() accumulates
the partial Hough map derivative *phmd to *hd by adding +1 or
-1 to the pixels corresponding to the left or right borders respectively. 
It takes into account corrections due to border effects as well.
 
The function  LALHOUGHIntegrHD2HT() constructs a total Hough map 
HOUGHMapTotal *ht from its derivative HOUGHMapDeriv *hd by
integrating each row (x-direction).


 */

/************************************ <lalVerbatim file="HoughMapCV">
Author: Sintes, A. M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Module \texttt{HoughMap.c}}
\label{ss:HoughMap.c}
Subroutines for 
initialization and construction of Hough-map derivatives and total Hough-maps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{HoughMapD}
\index{\verb&LALHOUGHInitializeHD()&}
\index{\verb&LALHOUGHInitializeHT()&}
\index{\verb&LALHOUGHAddPHMD2HD()&}
\index{\verb&LALHOUGHIntegrHD2HT()&}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

${}$


The function \verb&LALHOUGHInitializeHD()& initializes the Hough map derivative
space \verb&HOUGHMapDeriv *hd& to zero. Note that the length of the map
\verb@hd->map@ should be \verb@hd->ySide * (hd->xSide + 1)@.\\

The function \verb&LALHOUGHInitializeHT()& initializes the total Hough map 
\verb&HOUGHMapTotal *ht&  to zero and checks consistency between 
the number of physical pixels in the
map  and  those given by the grid information structure 
\verb&HOUGHPatchGrid *patch&.\\

 Given an initial Hough map derivative \verb@HOUGHMapDeriv *hd@ and a representation
  of a {\sc phmd}
 \verb@HOUGHphmd *phmd@, the function \verb&LALHOUGHAddPHMD2HD()& accumulates
 the partial Hough map derivative \verb@*phmd@ to \verb@*hd@ by adding $+1$ or
 $-1$ to
 the pixels corresponding to the {\it left} or {\it right} borders respectively.
 It takes into account corrections due to {\it border} effects as well.\\
 
The function \verb&LALHOUGHIntegrHD2HT()& constructs a total Hough map 
\verb&HOUGHMapTotal *ht& from its derivative \verb@HOUGHMapDeriv *hd@ by
integrating each row (x-direction).
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
%%\begin{verbatim}
%%LALZDestroyVector()
%%\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{HoughMapCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
</lalLaTeX> */

#define PREFETCH_OPT 1
#define EAH_OPTIMIZATION 1
#define EAH_OPTIMIZATION_HM 1


#include <lal/HoughMap.h>

NRCSID (HOUGHMAPC, "$Id$");

/* specials for Apples assembler */ 
#ifdef __APPLE__ 
#define AD_FLOAT ".single " 
#define AD_ASCII ".ascii " 
#define AD_ALIGN16 ".align 4" 
#define AD_ALIGN32 ".align 5" 
#define AD_ALIGN64 ".align 6" 
#else /* x86 gas */ 
#define AD_FLOAT ".float " 
#define AD_ASCII ".string " 
#define AD_ALIGN16 ".align 16" 
#define AD_ALIGN32 ".align 32" 
#define AD_ALIGN64 ".align 64" 
#endif 

#if defined(PREFETCH_OPT)

#if defined(__INTEL_COMPILER) ||  defined(_MSC_VER)
// not tested yet with icc or MS Visual C 
#include <mmintrin.h>
#define PREFETCH(a) _mm_prefetch(a,_MM_HINT_T0)
#elif defined(__GNUC__)
#define PREFETCH(a) __builtin_prefetch(a)
#else
#define PREFETCH(a) a
#endif
#else 
#define PREFETCH(a) a
#endif

/*
 * The functions that make up the guts of this module
 */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="HoughMapD"> */
void LALHOUGHInitializeHD (LALStatus      *status,
			  HOUGHMapDeriv   *hd) /* the Hough map derivative */
{ /*   *********************************************  </lalVerbatim> */

  INT4     k, maxk;
  HoughDT  *pointer;

   /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHInitializeHD", HOUGHMAPC);
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

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="HoughMapD"> */
void LALHOUGHInitializeHT (LALStatus      *status,
			   HOUGHMapTotal   *ht,     /* the total Hough map */
			   HOUGHPatchGrid  *patch) /* patch information */
{ /*   *********************************************  </lalVerbatim> */

  INT4     k,maxk;
  HoughTT  *pointer;

   /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHInitializeHT", HOUGHMAPC);

  /*   Make sure the arguments are not NULL: */ 
  ASSERT (ht, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  ASSERT (patch, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  /* Make sure the map contains some pixels */
  ASSERT (ht->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (ht->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (patch->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (patch->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  /*   Check consistency between patch and ht size  (size mismatch) */
  ASSERT (ht->xSide==patch->xSide, status, HOUGHMAPH_ESZMM, HOUGHMAPH_MSGESZMM);
  ASSERT (ht->ySide==patch->ySide, status, HOUGHMAPH_ESZMM, HOUGHMAPH_MSGESZMM);

  /* -------------------------------------------   */

  /* number of physical pixels */
  maxk = ht->ySide * ht->xSide;
 
  /* initializing the Hough map space */
   pointer = &(ht->map[0]); 
   for ( k=0; k< maxk; ++k ){ 
     *pointer = 0; 
     ++pointer;  
   } 

  /* normal exit */
  RETURN (status);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/** Adds a hough map derivative into a total hough map derivative WITHOUT 
    taking into account the weight of the partial hough map */
/* *******************************  <lalVerbatim file="HoughMapD"> */
void LALHOUGHAddPHMD2HD (LALStatus      *status, /**< the status pointer */
			 HOUGHMapDeriv  *hd,  /**< the Hough map derivative */
			 HOUGHphmd      *phmd) /**< info from a partial map */ 
{ /*   *********************************************  </lalVerbatim> */

  INT2     k,j;
  INT2     yLower, yUpper;
  UINT2    lengthLeft,lengthRight, xSide,ySide;
  COORType     *xPixel;
  HOUGHBorder  *borderP;

   /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHAddPHMD2HD", HOUGHMAPC);
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

  lengthLeft = phmd->lengthLeft;
  lengthRight= phmd->lengthRight;
  

  if(lengthLeft > 0) {
      PREFETCH(&(phmd->leftBorderP[0]->xPixel[0]));
  }	

  /* first column correction */
  for ( k=0; k< ySide; ++k ){
    hd->map[k*(xSide+1) + 0] += phmd->firstColumn[k];
  }

  /* left borders =>  +1 increase */
  for (k=0; k< lengthLeft; ++k){

    /*  Make sure the arguments are not NULL: (Commented for performance) */ 
    /*  ASSERT (phmd->leftBorderP[k], status, HOUGHMAPH_ENULL,
	HOUGHMAPH_MSGENULL); */

    borderP = phmd->leftBorderP[k];
    xPixel =  &( (*borderP).xPixel[0] );

    yLower = (*borderP).yLower;
    yUpper = (*borderP).yUpper;
    
    if(k < lengthLeft-1) {
	PREFETCH(&(phmd->leftBorderP[k+1]->xPixel[0]));
    } 	

   
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


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/** Adds a hough map derivative into a total hough map derivative taking into
    account the weight of the partial hough map */
/* *******************************  <lalVerbatim file="HoughMapD"> */
void LALHOUGHAddPHMD2HD_W (LALStatus      *status, /**< the status pointer */
			   HOUGHMapDeriv  *hd,  /**< the Hough map derivative */
			   HOUGHphmd      *phmd) /**< info from a partial map */ 
{ /*   *********************************************  </lalVerbatim> */

  INT4     k,j;
  INT4     yLower, yUpper;
  UINT4    lengthLeft,lengthRight, xSide,ySide,xSideP1,xSideP1_2,xSideP1_3;
  COORType     *xPixel;
  HOUGHBorder  *borderP;
  HoughDT    weight;
  register HoughDT    tempM0,tempM1,tempM2,tempM3;
  INT4       sidx,sidx0,sidx1,sidx2,sidx3,sidxBase, sidxBase_n; /* pre-calcuted array index for sanity check */
  INT4	     c_c, c_n ,offs;

  HoughDT  *map_pointer;
  HoughDT  *pf_addr[8]; 



   /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHAddPHMD2HD_W", HOUGHMAPC);
  ATTATCHSTATUSPTR (status); 

  /*   Make sure the arguments are not NULL: */ 
  ASSERT (hd,   status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);
  ASSERT (phmd, status, HOUGHMAPH_ENULL, HOUGHMAPH_MSGENULL);

  PREFETCH(phmd->leftBorderP);
  PREFETCH(phmd->rightBorderP);


  /* -------------------------------------------   */
  /* Make sure the map contains some pixels */
  ASSERT (hd->xSide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
  ASSERT (hd->ySide, status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);

  weight = phmd->weight;

  xSide = hd->xSide;
  ySide = hd->ySide;
  xSideP1=xSide+1;
  xSideP1_2=xSideP1+xSideP1;
  xSideP1_3=xSideP1_2+xSideP1;
  map_pointer = &( hd->map[0]);

  lengthLeft = phmd->lengthLeft;
  lengthRight= phmd->lengthRight;
  

  if(lengthLeft > 0) {
      borderP = phmd->leftBorderP[0];
#ifdef PREFETCH_OPT
      PREFETCH(&(borderP->xPixel[borderP->yLower]));
#endif
  }	

  if(lengthRight > 0) {
      borderP = phmd->rightBorderP[0];
#ifdef PREFETCH_OPT
      PREFETCH(&(borderP->xPixel[borderP->yLower]));
#endif
  }	
  
  /* first column correction */
  for ( k=0; k< ySide; ++k ){
    map_pointer[k*(xSide+1) + 0] += phmd->firstColumn[k] * weight;
  }


  /* left borders =>  increase according to weight*/
  for (k=0; k< lengthLeft; ++k){

    /*  Make sure the arguments are not NULL: (Commented for performance) */ 
    /*  ASSERT (phmd->leftBorderP[k], status, HOUGHMAPH_ENULL,
	HOUGHMAPH_MSGENULL); */

    borderP = phmd->leftBorderP[k];
    xPixel =  &( (*borderP).xPixel[0] );

    yLower = (*borderP).yLower;
    yUpper = (*borderP).yUpper;

#ifdef PREFETCH_OPT
    if(k < lengthLeft-1) {
	INT4 ylkp1 = phmd->leftBorderP[k+1]->yLower;
	PREFETCH(&(phmd->leftBorderP[k+1]->xPixel[ylkp1]));
    } 	

    if(k < lengthLeft-2) {
	PREFETCH(phmd->leftBorderP[k+2]);
    } 	
#endif

   
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


#if defined(EAH_OPTIMIZATION) && (EAH_OPTIMIZATION_HM == 1)



/* don't clobber ebx , used for PIC on Mac OS */

__asm __volatile (
	"push %%ebx				\n\t"
	"mov %[xPixel], %%eax  			\n\t"
	"mov %[yLower], %%ebx  			\n\t"
	"lea (%%eax,%%ebx,0x2), %%esi  		\n\t"
	"mov %[xSideP1], %%edx   		\n\t"

	"mov %[yUpper] , %%edi  			\n\t"
	"lea -0x2(%%eax,%%edi,0x2),%%eax  	\n\t"
	
	"mov %[map] , %%edi  			\n\t"
	"mov %%ebx,%%ecx  			\n\t"
	"imul %%edx, %%ecx  			\n\t"	
	"lea (%%edi, %%ecx, 0x8), %%edi  	\n\t"
	"fldl %[w]  				\n\t"

	"cmp  %%eax,%%esi  			\n\t"
	"jmp  2f 				\n\t"

	AD_ALIGN32 "\n"
	"1:  					\n\t"
	
	"movzwl (%%esi),%%ebx			\n\t"
	"movzwl 2(%%esi),%%ecx			\n\t"
		
	"lea (%%edi, %%ebx, 0x8) , %%ebx  	\n\t"
	"fldl (%%ebx)  				\n\t"	
	"lea (%%edi,%%edx,0x8) , %%edi  	\n\t"
	"lea (%%edi,%%ecx,0x8) , %%ecx   	\n\t"
	"fldl (%%ecx)  				\n\t"

	"fxch %%st(1)			\n\t"
	"fadd %%st(2),%%st  			\n\t"
	"fstpl (%%ebx)  			\n\t"
	"fadd %%st(1),%%st	  		\n\t"	
	"fstpl (%%ecx)  			\n\t"
	"lea (%%edi,%%edx,0x8), %%edi  	\n\t"	

	"lea 4(%%esi) , %%esi  		\n\t"
	"cmp  %%eax,%%esi  		\n\t"
	"2:	  				\n\t"

	"jbe 1b	  				\n\t"
	"add $0x2,%%eax				\n\t"
	"cmp %%eax,%%esi			\n\t"
	"jne 3f  				\n\t"


	"movzwl (%%esi) , %%ebx  		\n\t"
	"lea (%%edi, %%ebx, 0x8) , %%ebx  	\n\t"
	"fldl (%%ebx)  				\n\t"
	"fadd %%st(1),%%st  			\n\t"
	"fstpl (%%ebx)  			\n\t"
	
	"3:  					\n\t"
	
	"fstp %%st  				\n\t"
	"pop %%ebx				\n\t"
	
: 
: [xPixel] "m" (xPixel) , [yLower] "m" (yLower) , [yUpper] "m" (yUpper), [xSideP1] "m" (xSideP1) , [map] "m" (map_pointer) , [w] "m" (weight)
: "eax", "ecx", "edx", "esi","edi","cc", "st(7)","st(6)", "st(5)"

);


#elif defined(EAH_OPTIMIZATION_HM) && (EAH_OPTIMIZATION_HM == 2) 
#if  defined(EAH_HM_BATCHSIZE) && (EAH_HM_BATCHSIZE == 4)
#define BATCHSIZE 4
#define BATCHSIZE_LOG2 2
#elif defined(EAH_HM_BATCHSIZE) && (EAH_HM_BATCHSIZE == 2)
#define BATCHSIZE 2
#define BATCHSIZE_LOG2 1
#else 
#define BATCHSIZE 1
#define BATCHSIZE_LOG2 0
#endif

    sidxBase=yLower*xSideP1;
    sidxBase_n = sidxBase+(xSideP1 << BATCHSIZE_LOG2);	
    /* fill first cache entries */	

    
    c_c =0;
    c_n =BATCHSIZE;

    offs = yUpper - yLower+1;
    if (offs > BATCHSIZE) {
	offs = BATCHSIZE; 
    }	
	
    	
    for(j=yLower; j < yLower+offs; j++) {
        PREFETCH(pf_addr[c_c++] = map_pointer + xPixel[j] + j*xSideP1);			
    }		
		
    c_c=0;
    for(j=yLower; j<=yUpper-(2*BATCHSIZE-1);j+=BATCHSIZE){

      	
      sidx0 = xPixel[j+BATCHSIZE]+sidxBase_n;; 
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      sidx1 = xPixel[j+(BATCHSIZE+1)]+sidxBase_n+xSideP1;
#endif
#if (BATCHSIZE == 4)
      sidx2 = xPixel[j+(BATCHSIZE+2)]+sidxBase_n+xSideP1_2;
      sidx3 = xPixel[j+(BATCHSIZE+3)]+sidxBase_n+xSideP1_3;;
#endif
	
      PREFETCH(xPixel +(j+(BATCHSIZE+BATCHSIZE)));

      PREFETCH(pf_addr[c_n] = map_pointer+sidx0);
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      PREFETCH(pf_addr[c_n+1] = map_pointer+sidx1);
#endif
#if (BATCHSIZE == 4)
      PREFETCH(pf_addr[c_n+2] = map_pointer+sidx2);
      PREFETCH(pf_addr[c_n+3] = map_pointer+sidx3);
#endif

#ifndef LAL_NDEBUG 
      if ((sidx0 < 0) || (sidx0 >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j,xPixel[j] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      if ((sidx1 < 0) || (sidx1 >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j+1,xPixel[j+1] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
#endif
#if (BATCHSIZE == 4)
      if ((sidx2 < 0) || (sidx2 >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx2,ySide*(xSide+1),j+2,xPixel[j+2] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
      if ((sidx3 < 0) || (sidx3 >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx3,ySide*(xSide+1),j+3,xPixel[j+3] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
#endif
#endif 

      tempM0 = *(pf_addr[c_c]) +weight;
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      tempM1 = *(pf_addr[c_c+1]) +weight;
#endif
#if (BATCHSIZE == 4)
      tempM2 = *(pf_addr[c_c+2]) +weight;
      tempM3 = *(pf_addr[c_c+3]) +weight;
#endif
      sidxBase = sidxBase_n;
      sidxBase_n+=xSideP1 << BATCHSIZE_LOG2;
      
      (*(pf_addr[c_c]))=tempM0;
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      (*(pf_addr[c_c+1]))=tempM1;
#endif
#if (BATCHSIZE == 4)
      (*(pf_addr[c_c+2]))=tempM2;
      (*(pf_addr[c_c+3]))=tempM3;
#endif 

      c_c ^= BATCHSIZE;
      c_n ^= BATCHSIZE;
    }


    for(; j<=yUpper;++j){
      sidx = j*xSideP1 + xPixel[j];
#ifndef LAL_NDEBUG
      if ((sidx < 0) || (sidx >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j,xPixel[j] );

	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
#endif
      map_pointer[sidx] += weight;
    }


#else
    for(j=yLower; j<=yUpper;++j){
      sidx = j *(xSide+1) + xPixel[j];
      if ((sidx < 0) || (sidx >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j,xPixel[j] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
      map_pointer[sidx] += weight;
    }

#endif

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

#ifdef PREFETCH_OPT
    if(k < lengthRight-1) {
	INT4 ylkp1 = phmd->rightBorderP[k+1]->yLower;
	PREFETCH(&(phmd->rightBorderP[k+1]->xPixel[ylkp1]));
    } 	

    if(k < lengthRight-2) {
	PREFETCH(phmd->rightBorderP[k+2]);
    } 	
#endif
   
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


#if defined(EAH_OPTIMIZATION_HM) && (EAH_OPTIMIZATION_HM == 1)

__asm __volatile (
	"push %%ebx				\n\t"
	"mov %[xPixel], %%eax  			\n\t"
	"mov %[yLower], %%ebx  			\n\t"
	"mov %[xSideP1], %%edx   		\n\t"
	"lea (%%eax,%%ebx,0x2), %%esi  		\n\t"

	"mov %[yUpper] , %%edi  		\n\t"
	"lea -0x2(%%eax,%%edi,0x2),%%eax  	\n\t"
	
	"mov %[map] , %%edi  			\n\t"
	"mov %%ebx,%%ecx  			\n\t"
	"imul %%edx, %%ecx  			\n\t"	
	"lea (%%edi, %%ecx, 0x8), %%edi  	\n\t"
	"fldl %[w]  				\n\t"

	"cmp  %%eax,%%esi  			\n\t"
	"jmp  2f 				\n\t"

	AD_ALIGN32 "\n"
	"1:  					\n\t"

	"movzwl (%%esi),%%ebx			\n\t"
	"movzwl 2(%%esi),%%ecx			\n\t"

	"lea (%%edi, %%ebx, 0x8) , %%ebx  	\n\t"
	"fldl (%%ebx)  				\n\t"	
	"lea (%%edi,%%edx,0x8) , %%edi  	\n\t"
	"lea (%%edi,%%ecx,0x8) , %%ecx   	\n\t"
	"fldl (%%ecx)  				\n\t"


	"fxch %%st(1)			\n\t"
	"fsub %%st(2),%%st  			\n\t"
	"fstpl (%%ebx)  			\n\t"
	"fsub %%st(1),%%st	  		\n\t"	
	"fstpl (%%ecx)  			\n\t"
	"lea (%%edi,%%edx,0x8), %%edi  	\n\t"	
	"lea 4(%%esi) , %%esi  		\n\t"
	"cmp  %%eax,%%esi  		\n\t"
	"2:	  				\n\t"

	"jbe 1b	  				\n\t"
	"add $0x2,%%eax				\n\t"
	"cmp %%eax,%%esi			\n\t"
	"jne 3f  				\n\t"


	"movzwl (%%esi) , %%ebx  		\n\t"
	"lea (%%edi, %%ebx, 0x8) , %%ebx  	\n\t"
	"fldl (%%ebx)  				\n\t"
	"fsub %%st(1),%%st  			\n\t"
	"fstpl (%%ebx)  			\n\t"
	
	"3:  					\n\t"
	
	"fstp %%st  				\n\t"
	"pop %%ebx				\n\t"
	
: 
: [xPixel] "m" (xPixel) , [yLower] "m" (yLower) , [yUpper] "m" (yUpper), [xSideP1] "m" (xSideP1) , [map] "m" (map_pointer) , [w] "m" (weight)
: "eax", "ecx", "edx", "esi", "edi","cc", "st(7)","st(6)", "st(5)"

);


#elif defined(EAH_OPTIMIZATION_HM) && (EAH_OPTIMIZATION_HM == 2)
#if  defined(EAH_HM_BATCHSIZE) && (EAH_HM_BATCHSIZE == 4)
#define BATCHSIZE 4
#define BATCHSIZE_LOG2 2
#elif defined(EAH_HM_BATCHSIZE) && (EAH_HM_BATCHSIZE == 2)
#define BATCHSIZE 2
#define BATCHSIZE_LOG2 1
#else 
#define BATCHSIZE 1
#define BATCHSIZE_LOG2 0
#endif

    sidxBase=yLower*xSideP1;
    sidxBase_n = sidxBase+(xSideP1 << BATCHSIZE_LOG2);	
    /* fill first cache entries */	

    
    c_c =0;
    c_n =BATCHSIZE;

    offs = yUpper - yLower+1;
    if (offs > BATCHSIZE) {
	offs = BATCHSIZE; 
    }	
	
    	
    for(j=yLower; j < yLower+offs; j++) {
        PREFETCH(pf_addr[c_c++] = map_pointer + xPixel[j] + j*xSideP1);			
    }		
		
    c_c=0;
    for(j=yLower; j<=yUpper-(BATCHSIZE*2-1);j+=BATCHSIZE){

      	
      sidx0 = xPixel[j+BATCHSIZE]+sidxBase_n;; 
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      sidx1 = xPixel[j+(BATCHSIZE+1)]+sidxBase_n+xSideP1;
#endif
#if (BATCHSIZE == 4) 
      sidx2 = xPixel[j+(BATCHSIZE+2)]+sidxBase_n+xSideP1_2;
      sidx3 = xPixel[j+(BATCHSIZE+3)]+sidxBase_n+xSideP1_3;;
#endif
	
      PREFETCH(xPixel +(j+(BATCHSIZE+BATCHSIZE)));

      PREFETCH(pf_addr[c_n] = map_pointer+sidx0);
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      PREFETCH(pf_addr[c_n+1] = map_pointer+sidx1);
#endif
#if (BATCHSIZE == 4) 
      PREFETCH(pf_addr[c_n+2] = map_pointer+sidx2);
      PREFETCH(pf_addr[c_n+3] = map_pointer+sidx3);
#endif

#ifndef LAL_NDEBUG 
      if ((sidx0 < 0) || (sidx0 >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j,xPixel[j] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      if ((sidx1 < 0) || (sidx1 >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j+1,xPixel[j+1] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
#endif
#if (BATCHSIZE == 4) 
      if ((sidx2 < 0) || (sidx2 >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx2,ySide*(xSide+1),j+2,xPixel[j+2] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
      if ((sidx3 < 0) || (sidx3 >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx3,ySide*(xSide+1),j+3,xPixel[j+3] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
#endif

#endif 

      tempM0 = *(pf_addr[c_c]) -weight;
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      tempM1 = *(pf_addr[c_c+1]) -weight;
#endif
#if (BATCHSIZE == 4) 
      tempM2 = *(pf_addr[c_c+2]) -weight;
      tempM3 = *(pf_addr[c_c+3]) -weight;
#endif

      sidxBase = sidxBase_n;
      sidxBase_n+=xSideP1 << BATCHSIZE_LOG2;
      
      (*(pf_addr[c_c]))=tempM0;
#if (BATCHSIZE == 4) || (BATCHSIZE == 2)
      (*(pf_addr[c_c+1]))=tempM1;
#endif
#if (BATCHSIZE == 4) 
      (*(pf_addr[c_c+2]))=tempM2;
      (*(pf_addr[c_c+3]))=tempM3;
#endif
      c_c ^= BATCHSIZE;
      c_n ^= BATCHSIZE;
    }


    for(; j<=yUpper;++j){
      sidx = j*xSideP1 + xPixel[j];
#ifndef LAL_NDEBUG
      if ((sidx < 0) || (sidx >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j,xPixel[j] );

	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
#endif
      map_pointer[sidx] -= weight;
    }

#else

    for(j=yLower; j<=yUpper;++j){
      sidx = j*(xSide+1) + xPixel[j];
      if ((sidx < 0) || (sidx >= ySide*(xSide+1))) {
	fprintf(stderr,"\nERROR: %s %d: map index out of bounds: %d [0..%d] j:%d xp[j]:%d\n",
		__FILE__,__LINE__,sidx,ySide*(xSide+1),j,xPixel[j] );
	ABORT(status, HOUGHMAPH_ESIZE, HOUGHMAPH_MSGESIZE);
      }
      map_pointer[sidx] -= weight;
    }
#endif

  }


  /* -------------------------------------------   */
  
  DETATCHSTATUSPTR (status);
  
  /* normal exit */
  RETURN (status);
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="HoughMapD"> */
void LALHOUGHIntegrHD2HT (LALStatus       *status,
			  HOUGHMapTotal   *ht,     /* the total Hough map */
			  HOUGHMapDeriv   *hd) /* the Hough map derivative */
{ /*   *********************************************  </lalVerbatim> */

  INT2    i,j;
  UINT2   xSide,ySide;
  HoughTT accumulator;

   /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHIntegrHD2HT", HOUGHMAPC);
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
  ASSERT (ht->xSide==hd->xSide, status, HOUGHMAPH_ESZMM, HOUGHMAPH_MSGESZMM);
  ASSERT (ht->ySide==hd->ySide, status, HOUGHMAPH_ESZMM, HOUGHMAPH_MSGESZMM);
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

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="HoughMapD"> */
/**  Find source sky location given stereographic coordinates indexes */
void LALStereo2SkyLocation (LALStatus  *status,
         REAL8UnitPolarCoor *sourceLocation, /* output*/
	 UINT2              xPos,
	 UINT2              yPos,
	 HOUGHPatchGrid    *patch,
	 HOUGHDemodPar     *parDem)
{  /*   *********************************************  </lalVerbatim> */
    
  REAL8Cart2Coor        sourceProjected;
  REAL8UnitPolarCoor    sourceRotated;
  REAL8UnitPolarCoor    skyPatchCenter;
  /* --------------------------------------------- */
  INITSTATUS (status, "Stereo2SkyLocation", HOUGHMAPH);
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
