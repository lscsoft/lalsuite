/**************************** <lalVerbatim file="BilinearTransformCV">
Author: Creighton, T. D.
$Id$
***************************** </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{BilinearTransform.c}}
\label{ss:BilinearTransform.c}

Transforms the complex frequency coordinate of a ZPG filter.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{BilinearTransformCP}
\index{\verb&WToZCOMPLEX8ZPGFilter()&}
\index{\verb&WToZCOMPLEX16ZPGFilter()&}

\subsubsection*{Description}

These functions perform an in-place bilinear transformation on an
object \verb@*filter@ of type \verb@<datatype>ZPGFilter@, transforming
from $w$ to $z=(1+iw)/(1-iw)$.  Care is taken to ensure that zeros and
poles at $w=\infty$ are correctly transformed to $z=-1$, and zeros and
poles at $w=-i$ are correctly transformed to $z=\infty$.  In addition
to simply relocating the zeros and poles, residual factors are also
incorporated into the gain of the filter (i.e.\ the leading
coefficient of the rational function).

\subsubsection*{Algorithm}

The vectors \verb@filter->zeros@ and \verb@filter->poles@ only record
those zeros and poles that have finite value.  If one includes the
point $\infty$ on the complex plane, then a rational function always
has the same number of zeros and poles: a number \verb@num@ that is
the larger of \verb@z->zeros->length@ or \verb@z->poles->length@.  If
one or the other vector has a smaller length, then after the
transformation that vector will receive additional elements, with a
complex value of $z=-1$, to bring its length up to \verb@num@.
However, each vector will then \emph{lose} those elements that
previously had values $w=-i$, (which are sent to $z=\infty$,) thus
possibly decreasing the length of the vector.  These routines handle
this by simply allocating a new vector for the transformed data, and
freeing the old vector after the transformation.

When transforming a zero $w_k$ on the complex plane, one makes use of
the identity:
$$
(w - w_k) = -(w_k + i)\times\frac{z-z_k}{z+1} \; ,
$$
and similarly, when transforming a pole at $w_k$,
$$
(w - w_k)^{-1} = -(w_k + i)^{-1}\times\frac{z+1}{z-z_k} \; ,
$$
where $z=(1+iw)/(1-iw)$ and $z_k=(1+iw_k)/(1-iw_k)$.  If there are an
equal number of poles and zeros being transformed, then the factors of
$z+1$ will cancel; otherwise, the remaining factors correspond to the
zeros or poles at $z=-1$ brought in from $w=\infty$.  The factor
$(z-z_k)$ represents the transformation of the zero or pole at $w_k$.
The important factor to note, though, is the factor $-(w_k+i)^{\pm1}$.
This factor represents the change in the gain \verb@filter->gain@.
When $w_k=-i$, the transformation is slightly different:
$$
(w + i) = \frac{2i}{z+1} \; ;
$$
thus the gain correction factor is $2i$ (rather than 0) in this case.

The algorithms in this module computes and stores all the gain
correction factors before applying them to the gain.  The correction
factors are sorted in order of absolute magnitude, and are multiplied
together in small- and large-magnitude pairs.  In this way one reduces
the risk of overrunning the floating-point dynamical range during
intermediate calculations.

As a similar precaution, the routines in this module use the algorithm
discussed in the \verb@VectorOps@ package whenever they perform
complex division, to avoid intermediate results that mey be the
product of two large numbers.  When transforming $z=(1+iw)/(1-iw)$,
these routines also test for special cases (such as $w$ purely
imaginary) that have qualitatively significant results ($z$ purely
real), so that one doesn't end up with, for instance, an imaginary
part of $10^{-12}$ instead of 0.

\subsubsection*{Uses}
\begin{verbatim}
I4CreateVector()
SCreateVector()
DCreateVector()
CCreateVector()
ZCreateVector()
I4DestroyVector()
SDestroyVector()
DDestroyVector()
CDestroyVector()
ZDestroyVector()
SHeapIndex()
DHeapIndex()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{BilinearTransformCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include "AVFactories.h"
#include "VectorOps.h"
#include "Sort.h"
#include <math.h>
#include "ZPGFilter.h"

NRCSID(BILINEARTRANSFORMC,"$Id$");


/* <lalVerbatim file="BilinearTransformCP"> */
void WToZCOMPLEX8ZPGFilter(Status            *stat,
			   COMPLEX8ZPGFilter *filter)
{ /* </lalVerbatim> */
  INT4 i;        /* A counter. */
  INT4 j;        /* Another counter. */
  INT4 num;      /* The total number of zeros or poles. */
  INT4 numZeros; /* The number of finite zeros. */
  INT4 numPoles; /* The number of finite poles. */
  COMPLEX8 *a;   /* A zero or pole in the w plane. */
  COMPLEX8 *b;   /* A zero or pole in the z plane. */
  COMPLEX8 *g;   /* A gain correction factor. */
  COMPLEX8Vector *z=NULL;   /* Vector of zeros or poles in z plane. */
  COMPLEX8Vector *gain=NULL; /* Vector of gain correction factors. */
  COMPLEX8Vector null;       /* A vector of zero length. */
  REAL4Vector *absGain=NULL; /* Magnitudes of gain corrections. */
  INT4Vector *index=NULL;    /* Index array for sorting absGain. */

  INITSTATUS(stat,"WToZCOMPLEX8ZPGFilter",BILINEARTRANSFORMC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure the filter pointer is non-null. */
  ASSERT(filter,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);

  /* If the filter->zeros or filter->poles pointers is null, this
     means that there are no zeros or no poles.  For simplicity, we
     set the pointer to point to a vector of zero length. */
  null.length=0;
  null.data=NULL;
  if(!filter->zeros)
    filter->zeros=&null;
  if(!filter->poles)
    filter->poles=&null;

  /* Check that the vector lengths are non-negative, and, if positive,
     that the vector data pointer is non-null. */
  numZeros=filter->zeros->length;
  ASSERT(numZeros>=0,stat,ZPGFILTER_EBAD,ZPGFILTER_MSGEBAD);
  if(numZeros>0)
    ASSERT(filter->zeros->data,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);
  numPoles=filter->poles->length;
  ASSERT(numPoles>=0,stat,ZPGFILTER_EBAD,ZPGFILTER_MSGEBAD);
  if(numPoles>0)
    ASSERT(filter->poles->data,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);

  /* Compute the total number of zeros and poles in the w-plane,
     including those at w=infinity. */
  num = (numZeros>numPoles) ? numZeros : numPoles;
  numZeros=numPoles=num;

  /* If there are neither zeros nor poles, then there is nothing to
     transform. */
  if(num==0){
    filter->zeros=NULL;
    filter->poles=NULL;
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
  }

  /* Compute the revised number of zeros and poles in the z-plane,
     excluding those at z=infinity (w=-i). */
  for(i=0,a=filter->zeros->data;i<filter->zeros->length;i++,a++)
    if((a->re==0.0)&&(a->im==-1.0))
      numZeros--;
  for(i=0,a=filter->poles->data;i<filter->poles->length;i++,a++)
    if((a->re==0.0)&&(a->im==-1.0))
      numPoles--;

  /* Create the vector of gain correction factors. */
  TRY(CCreateVector(stat->statusPtr,&gain,filter->zeros->length+
		    filter->poles->length),stat);
  g=gain->data;

  /* Create the new vector of zeros. */
  TRY(CCreateVector(stat->statusPtr,&z,numZeros),stat);
  b=z->data;
  /* Transform existing zeros from w to z, except for those at w=-i,
     which are mapped to z=infinity.  At the same time, compute the
     gain correction factors. */
  for(i=0,j=0,a=filter->zeros->data;i<filter->zeros->length;
      i++,a++,g++){
    REAL4 ar=a->re;
    REAL4 ai=a->im;
    if(ar==0.0){
      if(ai==-1.0){
	/* w=-i is mapped to z=infinity. */
	g->re=0.0;
	g->im=2.0;
      }else{
	/* w=i*y is mapped to z=(1-y)/(1+y). */
	b->re=(1.0-ai)/(1.0+ai);
	b->im=0.0;
	g->re=0.0;
	g->im=-(1.0+ai);
	b++;
	j++;
      }
    }else if(fabs(1.0+ai)>fabs(ar)){
      REAL8 ratio = -ar/(1.0+ai);
      REAL8 denom = 1.0+ai - ratio*ar;

      b->re = (1.0-ai + ratio*ar)/denom;
      b->im = (ar - ratio*(1.0-ai))/denom;
      g->re = -ar;
      g->im = -(1.0+ai);
      b++;
      j++;
    }else{
      REAL8 ratio = -(1.0+ai)/ar;
      REAL8 denom = -ar + ratio*(1.0+ai);

      b->re = ((1.0-ai)*ratio + ar)/denom;
      b->im = (ar*ratio - 1.0+ai)/denom;
      g->re = -ar;
      g->im = -(1.0+ai);
      b++;
      j++;
    }
  }
  /* Transform any remaining zeros at w=infinity to z=-1. */
  for(;j<numZeros;b++,j++){
    b->re = -1.0;
    b->im = 0.0;
  }
  /* Replace the old filter zeros with the new ones. */
  if(filter->zeros->length>0)
    TRY(CDestroyVector(stat->statusPtr,&(filter->zeros)),stat);
  filter->zeros=z;
  z=NULL;

  /* Create the new vector of poles. */
  TRY(CCreateVector(stat->statusPtr,&z,numPoles),stat);
  b=z->data;
  /* Transform existing poles from w to z, except for those at w=-i,
     which are mapped to z=infinity.  At the same time, compute the
     gain correction factors. */
  for(i=0,j=0,a=filter->poles->data;i<filter->poles->length;
      i++,a++,g++){
    REAL4 ar=a->re;
    REAL4 ai=a->im;
    if(ar==0.0){
      if(ai==-1.0){
	/* w=-i is mapped to z=infinity. */
	g->re=0.0;
	g->im=-0.5;
      }else{
	/* w=i*y is mapped to z=(1-y)/(1+y). */
	b->re=(1.0-ai)/(1.0+ai);
	b->im=0.0;
	g->re=0.0;
	g->im=1.0/(1.0+ai);
	b++;
	j++;
      }
    }else if(fabs(1.0+ai)>fabs(ar)){
      REAL8 ratio = -ar/(1.0+ai);
      REAL8 denom = 1.0+ai - ratio*ar;

      b->re = (1.0-ai + ratio*ar)/denom;
      b->im = (ar - ratio*(1.0-ai))/denom;
      g->re = ratio/denom;
      g->im = 1.0/denom;
      b++;
      j++;
    }else{
      REAL8 ratio = -(1.0+ai)/ar;
      REAL8 denom = -ar + ratio*(1.0+ai);

      b->re = ((1.0-ai)*ratio + ar)/denom;
      b->im = (ar*ratio - 1.0+ai)/denom;
      g->re = ratio/denom;
      g->im = 1.0/denom;
      b++;
      j++;
    }
  }
  /* Transform any remaining poles at w=infinity to z=-1. */
  for(;j<numPoles;b++,j++){
    b->re = -1.0;
    b->im = 0.0;
  }
  /* Replace the old filter poles with the new ones. */
  if(filter->poles->length>0)
    TRY(CDestroyVector(stat->statusPtr,&(filter->poles)),stat);
  filter->poles=z;
  z=NULL;

  /* To avoid numerical overflow when applying the gain correction
     factors, we should multiply alternately by large and small
     factors.  First, create a vector of the factors' magnitudes. */
  TRY(SCreateVector(stat->statusPtr,&absGain,gain->length),stat);
  TRY(CVectorAbs(stat->statusPtr,absGain,gain),stat);

  /* Now create an index vector that indexes the magnitudes from small
     to large, and free the magnitude vector. */
  TRY(I4CreateVector(stat->statusPtr,&index,gain->length),stat);
  TRY(SHeapIndex(stat->statusPtr,index,absGain),stat);
  TRY(SDestroyVector(stat->statusPtr,&absGain),stat);

  /* Now multiply the gain alternately by small and large correction
     factors. */
  for(i=0,j=gain->length-1;i<j;i++,j--){
    /* Multiply the small and largest factors together. */
    REAL4 ar=gain->data[index->data[i]].re;
    REAL4 ai=gain->data[index->data[i]].im;
    REAL4 br=gain->data[index->data[j]].re;
    REAL4 bi=gain->data[index->data[j]].im;
    REAL4 cr=ar*br-ai*bi;
    REAL4 ci=ar*bi+ai*br;

    /* Multiply the gain by the combined factor. */
    br=filter->gain.re;
    bi=filter->gain.im;
    filter->gain.re=br*cr-bi*ci;
    filter->gain.im=br*ci+bi*cr;
  }
  if(i==j){
    /* Multiply by the remaining odd factor. */
    REAL4 cr=gain->data[index->data[i]].re;
    REAL4 ci=gain->data[index->data[i]].im;
    REAL4 br=filter->gain.re;
    REAL4 bi=filter->gain.im;

    filter->gain.re=br*cr-bi*ci;
    filter->gain.im=br*ci+bi*cr;
  }

  /* Free remaining temporary vectors, and exit. */
  TRY(CDestroyVector(stat->statusPtr,&gain),stat);
  TRY(I4DestroyVector(stat->statusPtr,&index),stat);
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="BilinearTransformCP"> */
void WToZCOMPLEX16ZPGFilter(Status             *stat,
			    COMPLEX16ZPGFilter *filter)
{ /* </lalVerbatim> */
  INT4 i;        /* A counter. */
  INT4 j;        /* Another counter. */
  INT4 num;      /* The total number of zeros or poles. */
  INT4 numZeros; /* The number of finite zeros. */
  INT4 numPoles; /* The number of finite poles. */
  COMPLEX16 *a;  /* A zero or pole in the w plane. */
  COMPLEX16 *b;  /* A zero or pole in the z plane. */
  COMPLEX16 *g;  /* A gain correction factor. */
  COMPLEX16Vector *z=NULL;  /* Vector of zeros or poles in z plane. */
  COMPLEX16Vector *gain=NULL; /* Vector of gain correction factors. */
  COMPLEX16Vector null;       /* A vector of zero length. */
  REAL8Vector *absGain=NULL;  /* Magnitudes of gain corrections. */
  INT4Vector *index=NULL;     /* Index array for sorting absGain. */

  INITSTATUS(stat,"WToZCOMPLEX16ZPGFilter",BILINEARTRANSFORMC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure the filter pointer is non-null. */
  ASSERT(filter,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);

  /* If the filter->zeros or filter->poles pointers is null, this
     means that there are no zeros or no poles.  For simplicity, we
     set the pointer to point to a vector of zero length. */
  null.length=0;
  null.data=NULL;
  if(!filter->zeros)
    filter->zeros=&null;
  if(!filter->poles)
    filter->poles=&null;

  /* Check that the vector lengths are non-negative, and, if positive,
     that the vector data pointer is non-null. */
  numZeros=filter->zeros->length;
  ASSERT(numZeros>=0,stat,ZPGFILTER_EBAD,ZPGFILTER_MSGEBAD);
  if(numZeros>0)
    ASSERT(filter->zeros->data,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);
  numPoles=filter->poles->length;
  ASSERT(numPoles>=0,stat,ZPGFILTER_EBAD,ZPGFILTER_MSGEBAD);
  if(numPoles>0)
    ASSERT(filter->poles->data,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);

  /* Compute the total number of zeros and poles in the w-plane,
     including those at w=infinity. */
  num = (numZeros>numPoles) ? numZeros : numPoles;
  numZeros=numPoles=num;

  /* If there are neither zeros nor poles, then there is nothing to
     transform. */
  if(num==0){
    filter->zeros=NULL;
    filter->poles=NULL;
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
  }

  /* Compute the revised number of zeros and poles in the z-plane,
     excluding those at z=infinity (w=-i). */
  for(i=0,a=filter->zeros->data;i<filter->zeros->length;i++,a++)
    if((a->re==0.0)&&(a->im==-1.0))
      numZeros--;
  for(i=0,a=filter->poles->data;i<filter->poles->length;i++,a++)
    if((a->re==0.0)&&(a->im==-1.0))
      numPoles--;

  /* Create the vector of gain correction factors. */
  TRY(ZCreateVector(stat->statusPtr,&gain,filter->zeros->length+
		    filter->poles->length),stat);
  g=gain->data;

  /* Create the new vector of zeros. */
  TRY(ZCreateVector(stat->statusPtr,&z,numZeros),stat);
  b=z->data;
  /* Transform existing zeros from w to z, except for those at w=-i,
     which are mapped to z=infinity.  At the same time, compute the
     gain correction factors. */
  for(i=0,j=0,a=filter->zeros->data;i<filter->zeros->length;
      i++,a++,g++){
    REAL8 ar=a->re;
    REAL8 ai=a->im;
    if(ar==0.0){
      if(ai==-1.0){
	/* w=-i is mapped to z=infinity. */
	g->re=0.0;
	g->im=2.0;
      }else{
	/* w=i*y is mapped to z=(1-y)/(1+y). */
	b->re=(1.0-ai)/(1.0+ai);
	b->im=0.0;
	g->re=0.0;
	g->im=-(1.0+ai);
	b++;
	j++;
      }
    }else if(fabs(1.0+ai)>fabs(ar)){
      REAL8 ratio = -ar/(1.0+ai);
      REAL8 denom = 1.0+ai - ratio*ar;

      b->re = (1.0-ai + ratio*ar)/denom;
      b->im = (ar - ratio*(1.0-ai))/denom;
      g->re = -ar;
      g->im = -(1.0+ai);
      b++;
      j++;
    }else{
      REAL8 ratio = -(1.0+ai)/ar;
      REAL8 denom = -ar + ratio*(1.0+ai);

      b->re = ((1.0-ai)*ratio + ar)/denom;
      b->im = (ar*ratio - 1.0+ai)/denom;
      g->re = -ar;
      g->im = -(1.0+ai);
      b++;
      j++;
    }
  }
  /* Transform any remaining zeros at w=infinity to z=-1. */
  for(;j<numZeros;b++,j++){
    b->re = -1.0;
    b->im = 0.0;
  }
  /* Replace the old filter zeros with the new ones. */
  if(filter->zeros->length>0)
    TRY(ZDestroyVector(stat->statusPtr,&(filter->zeros)),stat);
  filter->zeros=z;
  z=NULL;

  /* Create the new vector of poles. */
  TRY(ZCreateVector(stat->statusPtr,&z,numPoles),stat);
  b=z->data;
  /* Transform existing poles from w to z, except for those at w=-i,
     which are mapped to z=infinity.  At the same time, compute the
     gain correction factors. */
  for(i=0,j=0,a=filter->poles->data;i<filter->poles->length;
      i++,a++,g++){
    REAL8 ar=a->re;
    REAL8 ai=a->im;
    if(ar==0.0){
      if(ai==-1.0){
	/* w=-i is mapped to z=infinity. */
	g->re=0.0;
	g->im=-0.5;
      }else{
	/* w=i*y is mapped to z=(1-y)/(1+y). */
	b->re=(1.0-ai)/(1.0+ai);
	b->im=0.0;
	g->re=0.0;
	g->im=1.0/(1.0+ai);
	b++;
	j++;
      }
    }else if(fabs(1.0+ai)>fabs(ar)){
      REAL8 ratio = -ar/(1.0+ai);
      REAL8 denom = 1.0+ai - ratio*ar;

      b->re = (1.0-ai + ratio*ar)/denom;
      b->im = (ar - ratio*(1.0-ai))/denom;
      g->re = ratio/denom;
      g->im = 1.0/denom;
      b++;
      j++;
    }else{
      REAL8 ratio = -(1.0+ai)/ar;
      REAL8 denom = -ar + ratio*(1.0+ai);

      b->re = ((1.0-ai)*ratio + ar)/denom;
      b->im = (ar*ratio - 1.0+ai)/denom;
      g->re = ratio/denom;
      g->im = 1.0/denom;
      b++;
      j++;
    }
  }
  /* Transform any remaining poles at w=infinity to z=-1. */
  for(;j<numPoles;b++,j++){
    b->re = -1.0;
    b->im = 0.0;
  }
  /* Replace the old filter poles with the new ones. */
  if(filter->poles->length>0)
    TRY(ZDestroyVector(stat->statusPtr,&(filter->poles)),stat);
  filter->poles=z;
  z=NULL;

  /* To avoid numerical overflow, we should multiply alternately by
     large and small factors.  First, create a vector of the factors'
     magnitudes. */
  TRY(DCreateVector(stat->statusPtr,&absGain,gain->length),stat);
  TRY(ZVectorAbs(stat->statusPtr,absGain,gain),stat);

  /* Now create an index vector that indexes the magnitudes from small
     to large, and free the magnitude vector. */
  TRY(I4CreateVector(stat->statusPtr,&index,gain->length),stat);
  TRY(DHeapIndex(stat->statusPtr,index,absGain),stat);
  TRY(DDestroyVector(stat->statusPtr,&absGain),stat);

  /* Now multiply the gain alternately by small and large correction
     factors. */
  for(i=0,j=gain->length-1;i<j;i++,j--){
    /* Multiply by the small factor. */
    REAL8 ar=gain->data[index->data[i]].re;
    REAL8 ai=gain->data[index->data[i]].im;
    REAL8 gr=filter->gain.re;
    REAL8 gi=filter->gain.im;

    filter->gain.re=gr*ar-gi*gi;
    filter->gain.im=gr*ai+gi*ar;

    /* Multiply by the large factor. */
    ar=gain->data[index->data[j]].re;
    ai=gain->data[index->data[j]].im;
    gr=filter->gain.re;
    gi=filter->gain.im;

    filter->gain.re=gr*ar-gi*gi;
    filter->gain.im=gr*ai+gi*ar;
  }
  if(i==j){
    /* Multiply by the remaining odd factor. */
    REAL8 ar=gain->data[index->data[i]].re;
    REAL8 ai=gain->data[index->data[i]].im;
    REAL8 gr=filter->gain.re;
    REAL8 gi=filter->gain.im;

    filter->gain.re=gr*ar-gi*gi;
    filter->gain.im=gr*ai+gi*ar;
  }

  /* Free remaining temporary vectors, and exit. */
  TRY(ZDestroyVector(stat->statusPtr,&gain),stat);
  TRY(I4DestroyVector(stat->statusPtr,&index),stat);
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
