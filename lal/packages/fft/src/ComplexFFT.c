/*
*  Copyright (C) 2007 Jolien Creighton
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

#include <config.h>

#include <fftw3.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FFTWMutex.h>

/**
 * \addtogroup ComplexFFT_h
 *
 * \heading{Description}
 *
 * This package provides a (X)LAL-style interface with the FFTW fast Fourier
 * transform package [\ref fj_1998].
 *
 * The routines LALCreateForwardComplexFFTPlan() and
 * LALCreateReverseComplexFFTPlan() create plans for computing the
 * forward and reverse FFTs of a given size.  The optimum plan is either
 * estimated (reasonably fast) if the measure flag is zero, or measured (can be
 * time-consuming, but gives better performance) if the measure flag is
 * non-zero.  The routine LALDestroyComplexFFTPlan() destroys either
 * of these flavours of plans.
 *
 * The routine LALCOMPLEX8VectorFFT() performs either the forward or
 * reverse FFT depending on the plan.  The discrete Fourier transform \f$H_k\f$,
 * \f$k=0\ldots n-1\f$ of a vector \f$h_j\f$, \f$j=0\ldots n-1\f$, of length \f$n\f$ is defined
 * by
 * \f[
 *   H_k = \sum_{j=0}^{n-1} h_j e^{-2\pi ijk/n}
 * \f]
 * and, similarly, the \e inverse Fourier transform is defined by
 * \f[
 *   h_j = \frac{1}{n}\sum_{k=0}^{n-1} H_k e^{2\pi ijk/n}.
 * \f]
 * However, the present implementation of the \e reverse FFT omits the
 * factor of \f$1/n\f$.  The input and output vectors must be distinct.
 *
 * \heading{Operating Instructions}
 *
 * \code
 * const UINT4 n = 17;
 * static LALStatus status;
 * ComplexFFTPlan *pfwd = NULL;
 * ComplexFFTPlan *prev = NULL;
 * COMPLEX8Vector *avec = NULL;
 * COMPLEX8Vector *bvec = NULL;
 * COMPLEX8Vector *cvec = NULL;
 *
 * LALCreateForwardComplexFFTPlan( &status, &pfwd, n, 0 );
 * LALCreateReverseComplexFFTPlan( &status, &prev, n, 0 );
 * LALCCreateVector( &status, &avec, n );
 * LALCCreateVector( &status, &bvec, n );
 * LALCCreateVector( &status, &cvec, n );
 *
 * <assign data>
 *
 * LALCOMPLEX8VectorFFT( &status, bvec, avec, pfwd );
 * LALCOMPLEX8VectorFFT( &status, cvec, bvec, prev );
 *
 * LALDestroyComplexFFTPlan( &status, &pfwd );
 * LALDestroyComplexFFTPlan( &status, &prev );
 * LALCDestroyVector( &status, &avec );
 * LALCDestroyVector( &status, &bvec );
 * LALCDestroyVector( &status, &cvec );
 * \endcode
 *
 * \heading{Algorithm}
 *
 * The FFTW [\ref fj_1998] is used.
 *
 * \heading{Uses}
 *
 * \heading{Notes}
 *
 * <ol>
 * <li> The sign convention used here is the opposite of the definition in
 * <em>Numerical Recipes</em> [\ref ptvf1992], but agrees with the one used
 * by FFTW [\ref fj_1998] and the other LIGO software components.
 * </li><li> The result of the inverse FFT must be multiplied by \f$1/n\f$ to recover
 * the original vector.  This is different from the \c datacondAPI where
 * the factor is applied by default.
 * </li><li> The size \f$n\f$ of the transform can be any positive integer; the
 * performance is \f$O(n\log n)\f$.  However, better performance is obtained if \f$n\f$
 * is the product of powers of 2, 3, 5, 7, and zero or one power of either 11
 * or 13.  Transforms when \f$n\f$ is a power of 2 are especially fast.  See
 * Ref. [\ref fj_1998].
 * </li><li> LALMalloc() is used by all the fftw routines.
 * </li><li> The input and output vectors for LALCOMPLEX8VectorFFT() must
 * be distinct.
 * </li></ol>
 *
 *
 *
*/
/*@{*/

/** Plan to perform an FFT of COMPLEX8 data
 */
struct
tagCOMPLEX8FFTPlan
{
  INT4       sign; /**< sign in transform exponential, -1 for forward, +1 for reverse */
  UINT4      size; /**< length of the complex data vector for this plan */
  fftwf_plan plan; /**< the FFTW plan */
};

/** Plan to perform an FFT of COMPLEX16 data
 */
struct
tagCOMPLEX16FFTPlan
{
  INT4       sign; /**< sign in transform exponential, -1 for forward, +1 for reverse */
  UINT4      size; /**< length of the complex data vector for this plan */
  fftw_plan  plan; /**< the FFTW plan */
};


/*
 *
 * XLAL COMPLEX8 functions
 *
 */


COMPLEX8FFTPlan * XLALCreateCOMPLEX8FFTPlan( UINT4 size, int fwdflg, int measurelvl )
{
  COMPLEX8FFTPlan *plan;
  COMPLEX8 *tmp1;
  COMPLEX8 *tmp2;
  int flags = FFTW_UNALIGNED;

  if ( ! size )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  /* based on measurement level, set fftw3 flags to perform
   * requested degree of measurement */
  switch ( measurelvl )
  {
    case 0: /* estimate */
      flags |= FFTW_ESTIMATE;
      break;
    default: /* exhaustive measurement */
      flags |= FFTW_EXHAUSTIVE;
      /* fall-through */
    case 2: /* lengthy measurement */
      flags |= FFTW_PATIENT;
      /* fall-through */
    case 1: /* measure the best plan */
      flags |= FFTW_MEASURE;
      break;
  }

  /* allocate memory for the plan and the temporary arrays */
  plan = XLALMalloc( sizeof( *plan ) );
  tmp1 = XLALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = XLALMalloc( size * sizeof( *tmp2 ) );
  if ( ! plan || ! tmp1 || ! tmp2 )
  {
    XLALFree( plan );
    XLALFree( tmp1 );
    XLALFree( tmp2 );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  /* create the plan */
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  plan->plan = fftwf_plan_dft_1d( size,
      (fftwf_complex *)tmp1, (fftwf_complex *)tmp2,
      fwdflg ? FFTW_FORWARD : FFTW_BACKWARD, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* free temporary arrays */
  XLALFree( tmp2 );
  XLALFree( tmp1 );

  /* check to see success of plan creation */
  if ( ! plan->plan )
  {
    XLALFree( plan );
    XLAL_ERROR_NULL( XLAL_EFAILED );
  }

  /* now set remaining plan fields */
  plan->size = size;
  plan->sign = ( fwdflg ? -1 : 1 );

  return plan;
}


COMPLEX8FFTPlan * XLALCreateForwardCOMPLEX8FFTPlan( UINT4 size, int measurelvl )
{
  COMPLEX8FFTPlan *plan;
  plan = XLALCreateCOMPLEX8FFTPlan( size, 1, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return plan;
}


COMPLEX8FFTPlan * XLALCreateReverseCOMPLEX8FFTPlan( UINT4 size, int measurelvl )
{
  COMPLEX8FFTPlan *plan;
  plan = XLALCreateCOMPLEX8FFTPlan( size, 0, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return plan;
}


void XLALDestroyCOMPLEX8FFTPlan( COMPLEX8FFTPlan *plan )
{
  if ( plan )
  {
    if ( plan->plan )
    {
      LAL_FFTW_PTHREAD_MUTEX_LOCK;
      fftwf_destroy_plan( plan->plan );
      LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
    }
    memset( plan, 0, sizeof( *plan ) );
    XLALFree( plan );
  }
  return;
}


int XLALCOMPLEX8VectorFFT( COMPLEX8Vector *output, COMPLEX8Vector *input,
    const COMPLEX8FFTPlan *plan )
{
  if ( ! output || ! input || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! output->data || ! input->data || output->data == input->data )
    XLAL_ERROR( XLAL_EINVAL ); /* note: must be out-of-place */
  if ( output->length != plan->size || input->length != plan->size )
    XLAL_ERROR( XLAL_EBADLEN );

  /* do the fft */
  fftwf_execute_dft(
      plan->plan,
      (fftwf_complex *)input->data,
      (fftwf_complex *)output->data
      );
  return 0;
}



/*
 *
 * XLAL COMPLEX16 functions
 *
 */


COMPLEX16FFTPlan * XLALCreateCOMPLEX16FFTPlan( UINT4 size, int fwdflg, int measurelvl )
{
  COMPLEX16FFTPlan *plan;
  COMPLEX16 *tmp1;
  COMPLEX16 *tmp2;
  int flags = FFTW_UNALIGNED;

  if ( ! size )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  /* based on measurement level, set fftw3 flags to perform
   * requested degree of measurement */
  switch ( measurelvl )
  {
    case 0: /* estimate */
      flags |= FFTW_ESTIMATE;
      break;
    default: /* exhaustive measurement */
      flags |= FFTW_EXHAUSTIVE;
      /* fall-through */
    case 2: /* lengthy measurement */
      flags |= FFTW_PATIENT;
      /* fall-through */
    case 1: /* measure the best plan */
      flags |= FFTW_MEASURE;
      break;
  }

  /* allocate memory for the plan and the temporary arrays */
  plan = XLALMalloc( sizeof( *plan ) );
  tmp1 = XLALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = XLALMalloc( size * sizeof( *tmp2 ) );
  if ( ! plan || ! tmp1 || ! tmp2 )
  {
    XLALFree( plan );
    XLALFree( tmp1 );
    XLALFree( tmp2 );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  /* create the plan */
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  plan->plan = fftw_plan_dft_1d( size,
      (fftw_complex *)tmp1, (fftw_complex *)tmp2,
      fwdflg ? FFTW_FORWARD : FFTW_BACKWARD, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* free temporary arrays */
  XLALFree( tmp2 );
  XLALFree( tmp1 );

  /* check to see success of plan creation */
  if ( ! plan->plan )
  {
    XLALFree( plan );
    XLAL_ERROR_NULL( XLAL_EFAILED );
  }

  /* now set remaining plan fields */
  plan->size = size;
  plan->sign = ( fwdflg ? -1 : 1 );

  return plan;
}


COMPLEX16FFTPlan * XLALCreateForwardCOMPLEX16FFTPlan( UINT4 size, int measurelvl )
{
  COMPLEX16FFTPlan *plan;
  plan = XLALCreateCOMPLEX16FFTPlan( size, 1, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return plan;
}


COMPLEX16FFTPlan * XLALCreateReverseCOMPLEX16FFTPlan( UINT4 size, int measurelvl )
{
  COMPLEX16FFTPlan *plan;
  plan = XLALCreateCOMPLEX16FFTPlan( size, 0, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return plan;
}


void XLALDestroyCOMPLEX16FFTPlan( COMPLEX16FFTPlan *plan )
{
  if ( plan )
  {
    if ( plan->plan )
    {
      LAL_FFTW_PTHREAD_MUTEX_LOCK;
      fftw_destroy_plan( plan->plan );
      LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
    }
    memset( plan, 0, sizeof( *plan ) );
    XLALFree( plan );
  }
  return;
}


int XLALCOMPLEX16VectorFFT( COMPLEX16Vector *output, COMPLEX16Vector *input,
    const COMPLEX16FFTPlan *plan )
{
  if ( ! output || ! input || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! output->data || ! input->data || output->data == input->data )
    XLAL_ERROR( XLAL_EINVAL ); /* note: must be out-of-place */
  if ( output->length != plan->size || input->length != plan->size )
    XLAL_ERROR( XLAL_EBADLEN );

  /* do the fft */
  fftw_execute_dft(
      plan->plan,
      (fftw_complex *)input->data,
      (fftw_complex *)output->data
      );
  return 0;
}



/*
 *
 * LAL COMPLEX8 functions
 *
 */



void
LALCreateForwardCOMPLEX8FFTPlan(
    LALStatus       *status,
    COMPLEX8FFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCreateForwardCOMPLEX8FFTPlan", "XLALCreateForwardCOMPLEX8FFTPlan");

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  *plan = XLALCreateCOMPLEX8FFTPlan( size, 1, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}



void
LALCreateReverseCOMPLEX8FFTPlan(
    LALStatus       *status,
    COMPLEX8FFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCreateReverseCOMPLEX8FFTPlan", "XLALCreateReverseCOMPLEX8FFTPlan");

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  *plan = XLALCreateCOMPLEX8FFTPlan( size, 0, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}



void
LALDestroyCOMPLEX8FFTPlan (
    LALStatus       *status,
    COMPLEX8FFTPlan **plan
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALDestroyCOMPLEX8FFTPlan", "XLALDestroyCOMPLEX8FFTPlan");
  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( *plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  XLALDestroyCOMPLEX8FFTPlan( *plan );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT( status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}



void
LALCOMPLEX8VectorFFT (
    LALStatus      *status,
    COMPLEX8Vector *output,
    COMPLEX8Vector *input,
    COMPLEX8FFTPlan *plan
    )
{
  int code;
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCOMPLEX8VectorFFT", "XLALCOMPLEX8VectorFFT");

  ASSERT( output, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( input, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  ASSERT( output->data, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( input->data, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* make sure that it is not the same data! */
  ASSERT( output->data != input->data, status,
          COMPLEXFFTH_ESAME, COMPLEXFFTH_MSGESAME );

  /* make sure that the lengths agree */
  ASSERT( plan->size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
  ASSERT( output->length == plan->size, status,
          COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );
  ASSERT( input->length == plan->size, status,
          COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );

  code = XLALCOMPLEX8VectorFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        if ( ! plan->size ) /* plan size was invalid */
        {
          ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
        }
        else if ( output->data == input->data ) /* same data pointers */
        {
          ABORT( status, COMPLEXFFTH_ESAME, COMPLEXFFTH_MSGESAME );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/*
 *
 * LAL COMPLEX16 functions
 *
 */



void
LALCreateForwardCOMPLEX16FFTPlan(
    LALStatus       *status,
    COMPLEX16FFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCreateForwardCOMPLEX16FFTPlan", "XLALCreateForwardCOMPLEX16FFTPlan");

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  *plan = XLALCreateCOMPLEX16FFTPlan( size, 1, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}



void
LALCreateReverseCOMPLEX16FFTPlan(
    LALStatus       *status,
    COMPLEX16FFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCreateReverseCOMPLEX16FFTPlan", "XLALCreateReverseCOMPLEX16FFTPlan");

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  *plan = XLALCreateCOMPLEX16FFTPlan( size, 0, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}



void
LALDestroyCOMPLEX16FFTPlan (
    LALStatus       *status,
    COMPLEX16FFTPlan **plan
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALDestroyCOMPLEX16FFTPlan", "XLALDestroyCOMPLEX16FFTPlan");
  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( *plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  XLALDestroyCOMPLEX16FFTPlan( *plan );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT( status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}



void
LALCOMPLEX16VectorFFT (
    LALStatus      *status,
    COMPLEX16Vector *output,
    COMPLEX16Vector *input,
    COMPLEX16FFTPlan *plan
    )
{
  int code;
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCOMPLEX16VectorFFT", "XLALCOMPLEX16VectorFFT");

  ASSERT( output, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( input, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  ASSERT( output->data, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( input->data, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* make sure that it is not the same data! */
  ASSERT( output->data != input->data, status,
          COMPLEXFFTH_ESAME, COMPLEXFFTH_MSGESAME );

  /* make sure that the lengths agree */
  ASSERT( plan->size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
  ASSERT( output->length == plan->size, status,
          COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );
  ASSERT( input->length == plan->size, status,
          COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );

  code = XLALCOMPLEX16VectorFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        if ( ! plan->size ) /* plan size was invalid */
        {
          ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
        }
        else if ( output->data == input->data ) /* same data pointers */
        {
          ABORT( status, COMPLEXFFTH_ESAME, COMPLEXFFTH_MSGESAME );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}
/*@}*/
