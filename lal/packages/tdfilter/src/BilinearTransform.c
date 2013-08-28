/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <lal/Sort.h>
#include <math.h>
#include <lal/ZPGFilter.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * \addtogroup BilinearTransform_c
 * \author Creighton, T. D.
 *
 * \brief Transforms the complex frequency coordinate of a ZPG filter.
 *
 * \heading{Description}
 *
 * These functions perform an in-place bilinear transformation on an
 * object <tt>*filter</tt> of type <tt>\<datatype\>ZPGFilter</tt>, transforming
 * from \f$w\f$ to \f$z=(1+iw)/(1-iw)\f$.  Care is taken to ensure that zeros and
 * poles at \f$w=\infty\f$ are correctly transformed to \f$z=-1\f$, and zeros and
 * poles at \f$w=-i\f$ are correctly transformed to \f$z=\infty\f$.  In addition
 * to simply relocating the zeros and poles, residual factors are also
 * incorporated into the gain of the filter (i.e.\ the leading
 * coefficient of the rational function).
 *
 * \heading{Algorithm}
 *
 * The vectors <tt>filter->zeros</tt> and <tt>filter->poles</tt> only record
 * those zeros and poles that have finite value.  If one includes the
 * point \f$\infty\f$ on the complex plane, then a rational function always
 * has the same number of zeros and poles: a number \c num that is
 * the larger of <tt>z->zeros->length</tt> or <tt>z->poles->length</tt>.  If
 * one or the other vector has a smaller length, then after the
 * transformation that vector will receive additional elements, with a
 * complex value of \f$z=-1\f$, to bring its length up to \c num.
 * However, each vector will then \e lose those elements that
 * previously had values \f$w=-i\f$, (which are sent to \f$z=\infty\f$,) thus
 * possibly decreasing the length of the vector.  These routines handle
 * this by simply allocating a new vector for the transformed data, and
 * freeing the old vector after the transformation.
 *
 * When transforming a zero \f$w_k\f$ on the complex plane, one makes use of
 * the identity:
 * \f[
 * (w - w_k) = -(w_k + i)\times\frac{z-z_k}{z+1} \; ,
 * \f]
 * and similarly, when transforming a pole at \f$w_k\f$,
 * \f[
 * (w - w_k)^{-1} = -(w_k + i)^{-1}\times\frac{z+1}{z-z_k} \; ,
 * \f]
 * where \f$z=(1+iw)/(1-iw)\f$ and \f$z_k=(1+iw_k)/(1-iw_k)\f$.  If there are an
 * equal number of poles and zeros being transformed, then the factors of
 * \f$z+1\f$ will cancel; otherwise, the remaining factors correspond to the
 * zeros or poles at \f$z=-1\f$ brought in from \f$w=\infty\f$.  The factor
 * \f$(z-z_k)\f$ represents the new position of the transformed zero or pole.
 * The important factor to note, though, is the factor \f$-(w_k+i)^{\pm1}\f$.
 * This factor represents the change in the gain <tt>filter->gain</tt>.
 * When \f$w_k=-i\f$, the transformation is slightly different:
 * \f[
 * (w + i) = \frac{2i}{z+1} \; ;
 * \f]
 * thus the gain correction factor is \f$2i\f$ (rather than 0) in this case.
 *
 * The algorithm in this module computes and stores all the gain
 * correction factors before applying them to the gain.  The correction
 * factors are sorted in order of absolute magnitude, and are multiplied
 * together in small- and large-magnitude pairs.  In this way one reduces
 * the risk of overrunning the floating-point dynamical range during
 * intermediate calculations.
 *
 * As a similar precaution, the routines in this module use the algorithm
 * discussed in the \c VectorOps package whenever they perform
 * complex division, to avoid intermediate results that may be the
 * product of two large numbers.  When transforming \f$z=(1+iw)/(1-iw)\f$,
 * these routines also test for special cases (such as \f$w\f$ purely
 * imaginary) that have qualitatively significant results (\f$z\f$ purely
 * real), so that one doesn't end up with, for instance, an imaginary
 * part of \f$10^{-12}\f$ instead of 0.
 */
/*@{*/

/*
 * WARNING: NOT A PROPER COMPARE FUNCTION
 * this returns -1 if abs(a)<abs(b), otherwise it returns 1 (don't check equal)
 * also: don't worry about possible overflow
 */
static int CompareCOMPLEX8Abs( void UNUSED *p, const void *a, const void *b )
{
  REAL8 ar = crealf(*(const COMPLEX8 *)a);
  REAL8 ai = cimagf(*(const COMPLEX8 *)a);
  REAL8 br = crealf(*(const COMPLEX8 *)b);
  REAL8 bi = cimagf(*(const COMPLEX8 *)b);
  if ( (ar*ar+ai*ai) < (br*br+bi*bi) )
    return -1;
  return 1;
}

/*
 * WARNING: NOT A PROPER COMPARE FUNCTION
 * this returns -1 if abs(a)<abs(b), otherwise it returns 1 (don't check equal)
 * also: don't worry about possible overflow
 */
static int CompareCOMPLEX16Abs( void UNUSED *p, const void *a, const void *b )
{
  REAL8 ar = creal(*(const COMPLEX16 *)a);
  REAL8 ai = cimag(*(const COMPLEX16 *)a);
  REAL8 br = creal(*(const COMPLEX16 *)b);
  REAL8 bi = cimag(*(const COMPLEX16 *)b);
  if ( (ar*ar+ai*ai) < (br*br+bi*bi) )
    return -1;
  return 1;
}

/** \see See \ref BilinearTransform_c for documentation */
int XLALWToZCOMPLEX8ZPGFilter( COMPLEX8ZPGFilter *filter )
{
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
  INT4Vector *idx=NULL;    /* Index array for sorting absGain. */

  /* Make sure the filter pointer is non-null. */
  if ( ! filter )
    XLAL_ERROR( XLAL_EFAULT );

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
  if (numZeros<0)
    XLAL_ERROR(XLAL_EINVAL);
  if(numZeros>0)
    if (!filter->zeros->data)
      XLAL_ERROR(XLAL_EFAULT);
  numPoles=filter->poles->length;
  if (numPoles<0)
    XLAL_ERROR(XLAL_EINVAL);
  if(numPoles>0)
    if (!filter->poles->data)
      XLAL_ERROR(XLAL_EFAULT);

  /* Compute the total number of zeros and poles in the w-plane,
     including those at w=infinity. */
  num = (numZeros>numPoles) ? numZeros : numPoles;
  numZeros=numPoles=num;

  /* If there are neither zeros nor poles, then there is nothing to
     transform.  (The <0 case should never occur if the ASSERT()
     macros have done their job, but is included for extra safety.) */
  if(num<=0){
    filter->zeros=NULL;
    filter->poles=NULL;
    return 0;
  }

  /* Compute the revised number of zeros and poles in the z-plane,
     excluding those at z=infinity (w=-i). */
  for(i=0,a=filter->zeros->data;i<(INT4)filter->zeros->length;i++,a++)
    if((crealf(*a)==0.0)&&(cimagf(*a)==-1.0))
      numZeros--;
  for(i=0,a=filter->poles->data;i<(INT4)filter->poles->length;i++,a++)
    if((crealf(*a)==0.0)&&(cimagf(*a)==-1.0))
      numPoles--;

  /* Create the vector of gain correction factors. */
  /* Create the new vector of zeros. */
  gain=XLALCreateCOMPLEX8Vector(filter->zeros->length+filter->poles->length);
  z=XLALCreateCOMPLEX8Vector(numZeros);
  if (!gain||!z)
  {
    XLALDestroyCOMPLEX8Vector(gain);
    XLALDestroyCOMPLEX8Vector(z);
    XLAL_ERROR(XLAL_EFUNC);
  }
  g=gain->data;
  b=z->data;

  /* Transform existing zeros from w to z, except for those at w=-i,
     which are mapped to z=infinity.  At the same time, compute the
     gain correction factors. */
  for(i=0,j=0,a=filter->zeros->data;i<(INT4)filter->zeros->length;
      i++,a++,g++){
    REAL4 ar=crealf(*a);
    REAL4 ai=cimagf(*a);
    if(ar==0.0){
      if(ai==-1.0){
	/* w=-i is mapped to z=infinity. */
	*g=2.0*I;
      }else{
	/* w=i*y is mapped to z=(1-y)/(1+y). */
	*b=(1.0-ai)/(1.0+ai);
	*g=-(1.0+ai)*I;
	b++;
	j++;
      }
    }else if(fabs(1.0+ai)>fabs(ar)){
      REAL8 ratio = -ar/(1.0+ai);
      REAL8 denom = 1.0+ai - ratio*ar;

      *b = (1.0-ai + ratio*ar)/denom;
      *b += I*(ar - ratio*(1.0-ai))/denom;
      *g = -ar;
      *g += -(1.0+ai)*I;
      b++;
      j++;
    }else{
      REAL8 ratio = -(1.0+ai)/ar;
      REAL8 denom = -ar + ratio*(1.0+ai);

      *b = ((1.0-ai)*ratio + ar)/denom;
      *b += I*(ar*ratio - 1.0+ai)/denom;
      *g = -ar;
      *g = -(1.0+ai)*I;
      b++;
      j++;
    }
  }
  /* Transform any remaining zeros at w=infinity to z=-1. */
  for(;j<numZeros;b++,j++){
    *b = -1.0;
  }
  /* Replace the old filter zeros with the new ones. */
  if(filter->zeros->length>0)
    XLALDestroyCOMPLEX8Vector(filter->zeros);
  filter->zeros=z;
  z=NULL;

  /* Create the new vector of poles. */
  z=XLALCreateCOMPLEX8Vector(numPoles);
  if (!gain||!z)
  {
    XLALDestroyCOMPLEX8Vector(gain);
    XLAL_ERROR(XLAL_EFUNC);
  }
  b=z->data;
  /* Transform existing poles from w to z, except for those at w=-i,
     which are mapped to z=infinity.  At the same time, compute the
     gain correction factors. */
  for(i=0,j=0,a=filter->poles->data;i<(INT4)filter->poles->length;
      i++,a++,g++){
    REAL4 ar=crealf(*a);
    REAL4 ai=cimagf(*a);
    if(ar==0.0){
      if(ai==-1.0){
	/* w=-i is mapped to z=infinity. */
	*g=-0.5*I;
      }else{
	/* w=i*y is mapped to z=(1-y)/(1+y). */
	*b=(1.0-ai)/(1.0+ai);
	*g=I*1.0/(1.0+ai);
	b++;
	j++;
      }
    }else if(fabs(1.0+ai)>fabs(ar)){
      REAL8 ratio = -ar/(1.0+ai);
      REAL8 denom = 1.0+ai - ratio*ar;

      *b = (1.0-ai + ratio*ar)/denom;
      *b += I*(ar - ratio*(1.0-ai))/denom;
      *g = ratio/denom;
      *g += I*1.0/denom;
      b++;
      j++;
    }else{
      REAL8 ratio = -(1.0+ai)/ar;
      REAL8 denom = -ar + ratio*(1.0+ai);

      *b = ((1.0-ai)*ratio + ar)/denom;
      *b += I*(ar*ratio - 1.0+ai)/denom;
      *g = ratio/denom;
      *g += I*1.0/denom;
      b++;
      j++;
    }
  }
  /* Transform any remaining poles at w=infinity to z=-1. */
  for(;j<numPoles;b++,j++){
    *b = -1.0;
  }
  /* Replace the old filter poles with the new ones. */
  if(filter->poles->length>0)
    XLALDestroyCOMPLEX8Vector(filter->poles);
  filter->poles=z;
  z=NULL;

  /* To avoid numerical overflow when applying the gain correction
     factors, we should multiply alternately by large and small
     factors.  Create an idx vector that indexes the magnitudes
     from small to large. */
  idx=XLALCreateINT4Vector(gain->length);
  if(!idx||XLALHeapIndex(idx->data,gain->data,gain->length,sizeof(*gain->data),NULL,CompareCOMPLEX8Abs)<0)
  {
    XLALDestroyCOMPLEX8Vector(gain);
    XLALDestroyINT4Vector(idx);
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* Now multiply the gain alternately by small and large correction
     factors. */
  for(i=0,j=gain->length-1;i<j;i++,j--){
    /* Multiply the small and largest factors together. */
    /* Multiply the gain by the combined factor. */
    filter->gain *= gain->data[idx->data[i]] * gain->data[idx->data[j]];
  }
  if(i==j){
    /* Multiply by the remaining odd factor. */
    filter->gain *= gain->data[idx->data[i]];
  }

  /* Free remaining temporary vectors, and exit. */
  XLALDestroyCOMPLEX8Vector(gain);
  XLALDestroyINT4Vector(idx);
  return 0;
}

/** \see See \ref BilinearTransform_c for documentation */
int XLALWToZCOMPLEX16ZPGFilter( COMPLEX16ZPGFilter *filter )
{
  INT4 i;        /* A counter. */
  INT4 j;        /* Another counter. */
  INT4 num;      /* The total number of zeros or poles. */
  INT4 numZeros; /* The number of finite zeros. */
  INT4 numPoles; /* The number of finite poles. */
  COMPLEX16 *a;   /* A zero or pole in the w plane. */
  COMPLEX16 *b;   /* A zero or pole in the z plane. */
  COMPLEX16 *g;   /* A gain correction factor. */
  COMPLEX16Vector *z=NULL;   /* Vector of zeros or poles in z plane. */
  COMPLEX16Vector *gain=NULL; /* Vector of gain correction factors. */
  COMPLEX16Vector null;       /* A vector of zero length. */
  INT4Vector *idx=NULL;    /* Index array for sorting absGain. */

  /* Make sure the filter pointer is non-null. */
  if ( ! filter )
    XLAL_ERROR( XLAL_EFAULT );

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
  if (numZeros<0)
    XLAL_ERROR(XLAL_EINVAL);
  if(numZeros>0)
    if (!filter->zeros->data)
      XLAL_ERROR(XLAL_EFAULT);
  numPoles=filter->poles->length;
  if (numPoles<0)
    XLAL_ERROR(XLAL_EINVAL);
  if(numPoles>0)
    if (!filter->poles->data)
      XLAL_ERROR(XLAL_EFAULT);

  /* Compute the total number of zeros and poles in the w-plane,
     including those at w=infinity. */
  num = (numZeros>numPoles) ? numZeros : numPoles;
  numZeros=numPoles=num;

  /* If there are neither zeros nor poles, then there is nothing to
     transform.  (The <0 case should never occur if the ASSERT()
     macros have done their job, but is included for extra safety.) */
  if(num<=0){
    filter->zeros=NULL;
    filter->poles=NULL;
    return 0;
  }

  /* Compute the revised number of zeros and poles in the z-plane,
     excluding those at z=infinity (w=-i). */
  for(i=0,a=filter->zeros->data;i<(INT4)filter->zeros->length;i++,a++)
    if((creal(*a)==0.0)&&(cimag(*a)==-1.0))
      numZeros--;
  for(i=0,a=filter->poles->data;i<(INT4)filter->poles->length;i++,a++)
    if((creal(*a)==0.0)&&(cimag(*a)==-1.0))
      numPoles--;

  /* Create the vector of gain correction factors. */
  /* Create the new vector of zeros. */
  gain=XLALCreateCOMPLEX16Vector(filter->zeros->length+filter->poles->length);
  z=XLALCreateCOMPLEX16Vector(numZeros);
  if (!gain||!z)
  {
    XLALDestroyCOMPLEX16Vector(gain);
    XLALDestroyCOMPLEX16Vector(z);
    XLAL_ERROR(XLAL_EFUNC);
  }
  g=gain->data;
  b=z->data;

  /* Transform existing zeros from w to z, except for those at w=-i,
     which are mapped to z=infinity.  At the same time, compute the
     gain correction factors. */
  for(i=0,j=0,a=filter->zeros->data;i<(INT4)filter->zeros->length;
      i++,a++,g++){
    REAL8 ar=creal(*a);
    REAL8 ai=cimag(*a);
    if(ar==0.0){
      if(ai==-1.0){
	/* w=-i is mapped to z=infinity. */
	*g=2.0*I;
      }else{
	/* w=i*y is mapped to z=(1-y)/(1+y). */
	*b=(1.0-ai)/(1.0+ai);
	*g=-(1.0+ai)*I;
	b++;
	j++;
      }
    }else if(fabs(1.0+ai)>fabs(ar)){
      REAL8 ratio = -ar/(1.0+ai);
      REAL8 denom = 1.0+ai - ratio*ar;

      *b = (1.0-ai + ratio*ar)/denom;
      *b += I*(ar - ratio*(1.0-ai))/denom;
      *g = -ar;
      *g += -(1.0+ai)*I;
      b++;
      j++;
    }else{
      REAL8 ratio = -(1.0+ai)/ar;
      REAL8 denom = -ar + ratio*(1.0+ai);

      *b = ((1.0-ai)*ratio + ar)/denom;
      *b += I*(ar*ratio - 1.0+ai)/denom;
      *g = -ar;
      *g += -(1.0+ai)*I;
      b++;
      j++;
    }
  }
  /* Transform any remaining zeros at w=infinity to z=-1. */
  for(;j<numZeros;b++,j++){
    *b = -1.0;
  }
  /* Replace the old filter zeros with the new ones. */
  if(filter->zeros->length>0)
    XLALDestroyCOMPLEX16Vector(filter->zeros);
  filter->zeros=z;
  z=NULL;

  /* Create the new vector of poles. */
  z=XLALCreateCOMPLEX16Vector(numPoles);
  if (!gain||!z)
  {
    XLALDestroyCOMPLEX16Vector(gain);
    XLAL_ERROR(XLAL_EFUNC);
  }
  b=z->data;
  /* Transform existing poles from w to z, except for those at w=-i,
     which are mapped to z=infinity.  At the same time, compute the
     gain correction factors. */
  for(i=0,j=0,a=filter->poles->data;i<(INT4)filter->poles->length;
      i++,a++,g++){
    REAL8 ar=creal(*a);
    REAL8 ai=cimag(*a);
    if(ar==0.0){
      if(ai==-1.0){
	/* w=-i is mapped to z=infinity. */
	*g=-0.5*I;
      }else{
	/* w=i*y is mapped to z=(1-y)/(1+y). */
	*b=(1.0-ai)/(1.0+ai);
	*g=I*1.0/(1.0+ai);
	b++;
	j++;
      }
    }else if(fabs(1.0+ai)>fabs(ar)){
      REAL8 ratio = -ar/(1.0+ai);
      REAL8 denom = 1.0+ai - ratio*ar;

      *b = (1.0-ai + ratio*ar)/denom;
      *b += I*(ar - ratio*(1.0-ai))/denom;
      *g = ratio/denom;
      *g += I*1.0/denom;
      b++;
      j++;
    }else{
      REAL8 ratio = -(1.0+ai)/ar;
      REAL8 denom = -ar + ratio*(1.0+ai);

      *b = ((1.0-ai)*ratio + ar)/denom;
      *b += I*(ar*ratio - 1.0+ai)/denom;
      *g = ratio/denom;
      *g += I*1.0/denom;
      b++;
      j++;
    }
  }
  /* Transform any remaining poles at w=infinity to z=-1. */
  for(;j<numPoles;b++,j++){
    *b = -1.0;
  }
  /* Replace the old filter poles with the new ones. */
  if(filter->poles->length>0)
    XLALDestroyCOMPLEX16Vector(filter->poles);
  filter->poles=z;
  z=NULL;

  /* To avoid numerical overflow when applying the gain correction
     factors, we should multiply alternately by large and small
     factors.  Create an idx vector that indexes the magnitudes
     from small to large. */
  idx=XLALCreateINT4Vector(gain->length);
  if(!idx||XLALHeapIndex(idx->data,gain->data,gain->length,sizeof(*gain->data),NULL,CompareCOMPLEX16Abs)<0)
  {
    XLALDestroyCOMPLEX16Vector(gain);
    XLALDestroyINT4Vector(idx);
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* Now multiply the gain alternately by small and large correction
     factors. */
  for(i=0,j=gain->length-1;i<j;i++,j--){
    /* Multiply the small and largest factors together. */
    /* Multiply the gain by the combined factor. */
    filter->gain *= gain->data[idx->data[i]] * gain->data[idx->data[j]];
  }
  if(i==j){
    /* Multiply by the remaining odd factor. */
    filter->gain *= gain->data[idx->data[i]];
  }

  /* Free remaining temporary vectors, and exit. */
  XLALDestroyCOMPLEX16Vector(gain);
  XLALDestroyINT4Vector(idx);
  return 0;
}


/**
 * Deprecated.
 * \deprecated Use XLALWToZCOMPLEX8ZPGFilter() instead
 */
void
LALWToZCOMPLEX8ZPGFilter( LALStatus         *stat,
			  COMPLEX8ZPGFilter *filter )
{
  INITSTATUS(stat);

  if(XLALWToZCOMPLEX8ZPGFilter(filter)<0)
  {
    int code=xlalErrno;
    XLALClearErrno();
    switch(code){
      case XLAL_EFAULT:
        ABORT(stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
      case XLAL_EINVAL:
        ABORT(stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);
      default:
        ABORTXLAL(stat);
    }
  }

  RETURN(stat);
}


/**
 * Deprecated.
 * \deprecated Use XLALWToZCOMPLEX16ZPGFilter() instead
 */
void
LALWToZCOMPLEX16ZPGFilter( LALStatus          *stat,
			   COMPLEX16ZPGFilter *filter )
{
  INITSTATUS(stat);

  if(XLALWToZCOMPLEX16ZPGFilter(filter)<0)
  {
    int code=xlalErrno;
    XLALClearErrno();
    switch(code){
      case XLAL_EFAULT:
        ABORT(stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
      case XLAL_EINVAL:
        ABORT(stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);
      default:
        ABORTXLAL(stat);
    }
  }

  RETURN(stat);
}
/*@}*/
