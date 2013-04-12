/*
*  Copyright (C) 2010 Karsten Wiesner, Jolien Creighton, Kipp Cannon, Shin Kee Chung
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

#include <string.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/RealFFT.h>
#include <lal/XLALError.h>
#include <lal/FFTWMutex.h>
#include <fftw3.h>

#include "CudaFunctions.h"
#include "CudaFFT.h"

struct
tagREAL4FFTPlan
{
  INT4       sign;
  UINT4      size;
  cufftHandle plan;
  REAL4	    *d_real;
  COMPLEX8  *d_complex;
};

struct
tagREAL8FFTPlan
{
  INT4       sign;
  UINT4      size;
  fftw_plan  plan;
};


/*
 *
 * REAL4 XLAL Functions
 *
 */


REAL4FFTPlan * XLALCreateREAL4FFTPlan( UINT4 size, int fwdflg, int measurelvl )
{
  UINT4 createSize;
  REAL4FFTPlan *plan;

  if ( ! size )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  /* "use" measurelvl */
  measurelvl = 0;

  /* allocate memory for the plan and the temporary arrays */
  plan = XLALMalloc( sizeof( *plan ) );
  if ( ! plan )
  {
    XLALFree( plan );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  /*
   * Use a different size for plan creation to avoid current CUDA bug
   * in performing FFT to array with size 1
   */

  int retval;
	
  if( size == 1 ) createSize = 2;
  else	createSize = size;
  /* LAL_FFTW_PTHREAD_MUTEX_LOCK; */
  if ( fwdflg ) /* forward */
    retval= cufftPlan1d( &plan->plan, createSize, CUFFT_R2C, 1 );
  else /* reverse */
    retval= cufftPlan1d( &plan->plan, createSize, CUFFT_C2R, 1 );
  /* LAL_FFTW_PTHREAD_MUTEX_UNLOCK; */
  /* check to see success of plan creation */
	
  /* "Plan=0" Bugfix by Wiesner, K.: plan->plan is an integer handle not a pointer and 0 is a valid handle
      So checking against 0 and occasionaly destroy the plan is a bug.

      if ( ! plan->plan )
      {
         XLALFree( plan );
         XLAL_ERROR_NULL( XLAL_EFAILED );
      }
  */
	
  /* Allocate memory in the GPU */
  plan->d_real = XLALCudaMallocReal(size);
  plan->d_complex = XLALCudaMallocComplex(size/2 + 1);

  /* now set remaining plan fields */
  plan->size = size;
  plan->sign = ( fwdflg ? -1 : 1 );

  return plan;
}


REAL4FFTPlan * XLALCreateForwardREAL4FFTPlan( UINT4 size, int measurelvl )
{
  REAL4FFTPlan *plan;
  plan = XLALCreateREAL4FFTPlan( size, 1, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return plan;
}


REAL4FFTPlan * XLALCreateReverseREAL4FFTPlan( UINT4 size, int measurelvl )
{
  REAL4FFTPlan *plan;
  plan = XLALCreateREAL4FFTPlan( size, 0, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return plan;
}


void XLALDestroyREAL4FFTPlan( REAL4FFTPlan *plan )
{
  if ( ! plan )
    XLAL_ERROR_VOID( XLAL_EFAULT );

  /* Plan=0 Bugfix
    if ( ! plan->plan )
      XLAL_ERROR_VOID( XLAL_EINVAL );
  */

  /* LAL_FFTW_PTHREAD_MUTEX_LOCK; */
  /* Free the Cuda specific variables */
  XLALCudaFree( plan->d_real );
  cufftDestroy( plan->plan );
  XLALCudaFree( plan->d_complex );
  /* LAL_FFTW_PTHREAD_MUTEX_UNLOCK; */
  memset( plan, 0, sizeof( *plan ) );
  XLALFree( plan );
  return;
}

int XLALREAL4ForwardFFT( COMPLEX8Vector *output, const REAL4Vector *input,
    const REAL4FFTPlan *plan )
{
  if ( ! output || ! input || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  /* Plan=0 Bugfix
    if ( ! plan->plan || ! plan->size || plan->sign != -1 )
  */
  if ( ! plan->size || plan->sign != -1 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! output->data || ! input->data )
    XLAL_ERROR( XLAL_EINVAL );
  if ( input->length != plan->size || output->length != plan->size/2 + 1 )
    XLAL_ERROR( XLAL_EBADLEN );

  /* do the fft */
  /*
   * Special check for vector with size 1 to avoid CUDA bug
   */
  if( plan->size == 1 )
  {
    output->data[0].realf_FIXME = input->data[0];
    output->data[0].imagf_FIXME = 0.0;
  }
  else
    cudafft_execute_r2c( plan->plan,
    (cufftComplex *)(output->data), (cufftReal *)(input->data),
    (cufftComplex *)plan->d_complex, (cufftReal *)plan->d_real, plan->size );

  /* Nyquist frequency */
  if( plan->size%2 == 0 )
  {
    output->data[plan->size/2].imagf_FIXME = 0.0;
  }

  return 0;
}


int XLALREAL4ReverseFFT( REAL4Vector *output, const COMPLEX8Vector *input,
    const REAL4FFTPlan *plan )
{
  if ( ! output || ! input || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  /* Plan=0 Bugfix 
     if ( ! plan->plan || ! plan->size || plan->sign != 1 )
   */
  if ( ! plan->size || plan->sign != 1 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! output->data || ! input->data )
    XLAL_ERROR( XLAL_EINVAL );
  if ( output->length != plan->size || input->length != plan->size/2 + 1 )
    XLAL_ERROR( XLAL_EBADLEN );
  if ( cimagf(input->data[0]) != 0.0 )
    XLAL_ERROR( XLAL_EDOM );  /* imaginary part of DC must be zero */
  if ( ! plan->size % 2 && cimagf(input->data[plan->size/2]) != 0.0 )
    XLAL_ERROR( XLAL_EDOM );  /* imaginary part of Nyquist must be zero */

  /* perform the fft */
  /*
   * Special check to handle array of size 1
   */
  if( plan->size == 1 )
  {
    output->data[0] = crealf(input->data[0]);
  }
  else
    cudafft_execute_c2r( plan->plan,
    (cufftReal *)(output->data), (cufftComplex *)(input->data),
    (cufftReal *)plan->d_real, (cufftComplex *)plan->d_complex, plan->size );

  return 0;
}


int XLALREAL4VectorFFT( REAL4Vector *output, const REAL4Vector *input,
    const REAL4FFTPlan *plan )
{
  COMPLEX8 *tmp;
  UINT4 k;

  if ( ! output || ! input || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  /* Plan=0 Bugfix
     if ( ! plan->plan || ! plan->size )
  */
  if ( ! plan->size )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! output->data || ! input->data || output->data == input->data )
    XLAL_ERROR( XLAL_EINVAL ); /* note: must be out-of-place */
  if ( output->length != plan->size || input->length != plan->size )
    XLAL_ERROR( XLAL_EBADLEN );

  if( plan->size == 1 )
  {
    output->data[0] = input->data[0];
  }
  else
  {
    tmp = XLALMalloc( plan->size * sizeof(*tmp) );
    /* do the fft */
    if( plan->sign == -1 )
    {
      cudafft_execute_r2c( plan->plan,
	  (cufftComplex *)tmp, (cufftReal *)(input->data),
	  (cufftComplex *)plan->d_complex, (cufftReal *)plan->d_real, plan->size );
      output->data[0] = crealf(tmp[0]);

      for( k = 1; k < (plan->size + 1)/2; k++ )
      {
	output->data[k] = crealf(tmp[k]);
	output->data[plan->size - k] = cimagf(tmp[k]);
      }

      if( plan->size % 2 == 0 )
      {
	output->data[plan->size/2] = crealf(tmp[plan->size/2]);
      }
    }
    else
    {
      tmp[0].realf_FIXME = input->data[0];
      tmp[0].imagf_FIXME = 0.0;

      for( k = 1; k < (plan->size + 1)/2; k++ )
      {
	tmp[k].realf_FIXME = input->data[k];
	tmp[k].imagf_FIXME = input->data[plan->size - k];
      }

      if( plan->size%2 == 0 )
      {
	tmp[plan->size/2].realf_FIXME = input->data[plan->size/2];
	tmp[plan->size/2].imagf_FIXME = 0.0;
      }

      cudafft_execute_c2r( plan->plan,
	  (cufftReal *)(output->data), (cufftComplex *)tmp,
	  (cufftReal *)plan->d_real, (cufftComplex *)plan->d_complex, plan->size );
    }
    XLALFree(tmp);
  }

  return 0;
}


int XLALREAL4PowerSpectrum( REAL4Vector *spec, const REAL4Vector *data,
    const REAL4FFTPlan *plan )
{
  COMPLEX8 *tmp;
  UINT4 k;

  if ( ! spec || ! data || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  /* Plan=0 Bugfix  
     if ( ! plan->plan || ! plan->size )
  */
  if (! plan->size )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! spec->data || ! data->data )
    XLAL_ERROR( XLAL_EINVAL );
  if ( data->length != plan->size || spec->length != plan->size/2 + 1 )
    XLAL_ERROR( XLAL_EBADLEN );

  /* allocate temporary storage space */
  tmp = XLALMalloc( (plan->size/2 + 1) * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( XLAL_ENOMEM );

  /* Check for size 1 to avoid the CUDA bug */
  if( plan->size == 1 )
    tmp[0].realf_FIXME = data->data[0];
  /* transform the data */
  else
    cudafft_execute_r2c( plan->plan,
    (cufftComplex *)(tmp), (cufftReal *)(data->data),
    (cufftComplex *)plan->d_complex, (cufftReal *)plan->d_real, plan->size );

  /* now reconstruct the spectrum from the temporary storage */

  /* dc component */
  spec->data[0] = crealf(tmp[0]) * crealf(tmp[0]);

  /* other components */
  for (k = 1; k < (plan->size + 1)/2; ++k) /* k < size/2 rounded up */
  {
    REAL4 re = crealf(tmp[k]);
    REAL4 im = cimagf(tmp[k]);
    spec->data[k]  = re * re + im * im;
    spec->data[k] *= 2.0; /* accounts for negative frequency part */
  }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* size is even */
    spec->data[k] = crealf(tmp[k]) * crealf(tmp[k]);

  /* clenup and exit */
  XLALFree( tmp );
  return 0;
}



/*
 *
 * REAL8 XLAL Functions
 *
 */



REAL8FFTPlan * XLALCreateREAL8FFTPlan( UINT4 size, int fwdflg, int measurelvl )
{
  REAL8FFTPlan *plan;
  REAL8 *tmp1;
  REAL8 *tmp2;
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

  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  if ( fwdflg ) /* forward */
    plan->plan = fftw_plan_r2r_1d( size, tmp1, tmp2, FFTW_R2HC, flags );
  else /* reverse */
    plan->plan = fftw_plan_r2r_1d( size, tmp1, tmp2, FFTW_HC2R, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* free temporary arrays */
  XLALFree( tmp2 );
  XLALFree( tmp1 );

  /* check to see success of plan creation */
  /* Plan=0 Bugfix
  if ( ! plan->plan )
  {
    XLALFree( plan );
    XLAL_ERROR_NULL( XLAL_EFAILED );
  }
  */

  /* now set remaining plan fields */
  plan->size = size;
  plan->sign = ( fwdflg ? -1 : 1 );

  return plan;
}


REAL8FFTPlan * XLALCreateForwardREAL8FFTPlan( UINT4 size, int measurelvl )
{
  REAL8FFTPlan *plan;
  plan = XLALCreateREAL8FFTPlan( size, 1, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return plan;
}


REAL8FFTPlan * XLALCreateReverseREAL8FFTPlan( UINT4 size, int measurelvl )
{
  REAL8FFTPlan *plan;
  plan = XLALCreateREAL8FFTPlan( size, 0, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return plan;
}


void XLALDestroyREAL8FFTPlan( REAL8FFTPlan *plan )
{
  if ( ! plan )
    XLAL_ERROR_VOID( XLAL_EFAULT );
  /* Plan=0 Bugfix
     if ( ! plan->plan )
     XLAL_ERROR_VOID( XLAL_EINVAL );
  */
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  fftw_destroy_plan( plan->plan );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
  memset( plan, 0, sizeof( *plan ) );
  XLALFree( plan );

  return;
}

int XLALREAL8ForwardFFT( COMPLEX16Vector *output, REAL8Vector *input,
    const REAL8FFTPlan *plan )
{
  REAL8 *tmp;
  UINT4 k;

  if ( ! output || ! input || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  /* Plan=0 Bugfix
    if ( ! plan->plan || ! plan->size || plan->sign != -1 )
  */
  if ( ! plan->size || plan->sign != -1 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! output->data || ! input->data )
    XLAL_ERROR( XLAL_EINVAL );
  if ( input->length != plan->size || output->length != plan->size/2 + 1 )
    XLAL_ERROR( XLAL_EBADLEN );

  /* create temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( XLAL_ENOMEM );

  /* do the fft */
  fftw_execute_r2r( plan->plan, input->data, tmp );

  /* now unpack the results into the output vector */

  /* dc component */
  output->data[0].re = tmp[0];
  output->data[0].im = 0.0;

  /* other components */
  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size/2 rounded up */
  {
    output->data[k].re = tmp[k];
    output->data[k].im = tmp[plan->size - k];
  }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* n is even */
  {
    output->data[plan->size/2].re = tmp[plan->size/2];
    output->data[plan->size/2].im = 0.0;
  }

  XLALFree( tmp );
  return 0;
}


int XLALREAL8ReverseFFT( REAL8Vector *output, COMPLEX16Vector *input,
    const REAL8FFTPlan *plan )
{
  REAL8 *tmp;
  UINT4 k;

  if ( ! output || ! input || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  /* Plan=0 Bugfix:
     if ( ! plan->plan || ! plan->size || plan->sign != 1 )
  */
  if ( ! plan->size || plan->sign != 1 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! output->data || ! input->data )
    XLAL_ERROR( XLAL_EINVAL );
  if ( output->length != plan->size || input->length != plan->size/2 + 1 )
    XLAL_ERROR( XLAL_EBADLEN );
  if ( input->data[0].im != 0.0 )
    XLAL_ERROR( XLAL_EDOM );  /* imaginary part of DC must be zero */
  if ( ! plan->size % 2 && input->data[plan->size/2].im != 0.0 )
    XLAL_ERROR( XLAL_EDOM );  /* imaginary part of Nyquist must be zero */

  /* create temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( XLAL_ENOMEM );

  /* unpack input into temporary array */

  /* dc component */
  tmp[0] = input->data[0].re;

  /* other components */
  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size / 2 rounded up */
  {
    tmp[k]              = input->data[k].re;
    tmp[plan->size - k] = input->data[k].im;
  }

  /* Nyquist component */
  if ( plan->size%2 == 0 ) /* n is even */
    tmp[plan->size/2] = input->data[plan->size/2].re;

  /* perform the fft */
  fftw_execute_r2r( plan->plan, tmp, output->data );

  /* cleanup and exit */
  XLALFree( tmp );
  return 0;
}


int XLALREAL8VectorFFT( REAL8Vector *output, REAL8Vector *input,
    const REAL8FFTPlan *plan )
{
  if ( ! output || ! input || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  /* Plan=0 Bugfix
     if ( ! plan->plan || ! plan->size )
  */
  if (! plan->size )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! output->data || ! input->data || output->data == input->data )
    XLAL_ERROR( XLAL_EINVAL ); /* note: must be out-of-place */
  if ( output->length != plan->size || input->length != plan->size )
    XLAL_ERROR( XLAL_EBADLEN );

  /* do the fft */
  fftw_execute_r2r( plan->plan, input->data, output->data );
  return 0;
}


int XLALREAL8PowerSpectrum( REAL8Vector *spec, REAL8Vector *data,
    const REAL8FFTPlan *plan )
{
  REAL8 *tmp;
  UINT4 k;

  if ( ! spec || ! data || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  /* Plan=0 Bugfix  
     if ( ! plan->plan || ! plan->size )
  */
  if ( ! plan->size )
    XLAL_ERROR( XLAL_EINVAL );
  if ( ! spec->data || ! data->data )
    XLAL_ERROR( XLAL_EINVAL );
  if ( data->length != plan->size || spec->length != plan->size/2 + 1 )
    XLAL_ERROR( XLAL_EBADLEN );

  /* allocate temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( XLAL_ENOMEM );

  /* transform the data */
  fftw_execute_r2r( plan->plan, data->data, tmp );

  /* now reconstruct the spectrum from the temporary storage */

  /* dc component */
  spec->data[0] = tmp[0] * tmp[0];

  /* other components */
  for (k = 1; k < (plan->size + 1)/2; ++k) /* k < size/2 rounded up */
  {
    REAL8 re = tmp[k];
    REAL8 im = tmp[plan->size - k];
    spec->data[k]  = re * re + im * im;
    spec->data[k] *= 2.0; /* accounts for negative frequency part */
  }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* size is even */
    spec->data[plan->size/2] = tmp[plan->size/2] * tmp[plan->size/2];

  /* clenup and exit */
  XLALFree( tmp );
  return 0;
}




/*
 *
 * LAL Functions
 *
 */





void
LALCreateForwardREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCreateForwardREAL4FFTPlan", "XLALCreateForwardREAL4FFTPlan");

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  *plan = XLALCreateREAL4FFTPlan( size, 1, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}



void
LALCreateReverseREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCreateReverseREAL4FFTPlan", "XLALCreateReverseREAL4FFTPlan");

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  *plan = XLALCreateREAL4FFTPlan( size, 0, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}



void
LALDestroyREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALDestroyREAL4FFTPlan", "XLALDestroyREAL4FFTPlan");
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( *plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  XLALDestroyREAL4FFTPlan( *plan );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }
  *plan = NULL;
  RETURN( status );
}



void
LALForwardREAL4FFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    REAL4Vector    *input,
    REAL4FFTPlan    *plan
    )
{
  int code;
  UINT4 n;
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALForwardREAL4FFT", "XLALForwardREAL4FFT");

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  /* Plan=0 Bugfix
     ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  */

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( input->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( output->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  ASSERT( plan->sign == -1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  code = XLALREAL4ForwardFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( plan->sign != -1 ) /* plan sign was wrong */
        {
          ABORT( status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}



void
LALReverseREAL4FFT(
    LALStatus      *status,
    REAL4Vector    *output,
    COMPLEX8Vector *input,
    REAL4FFTPlan    *plan
    )
{
  int code;
  UINT4 n;
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALReverseREAL4FFT", "XLALReverseREAL4FFT");

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  /* Plan=0 Bugfix
     ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  */

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( output->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( cimagf(input->data[0]) == 0, status, REALFFTH_EDATA, REALFFTH_MSGEDATA );
  ASSERT( n % 2 || cimagf(input->data[n / 2]) == 0, status,
      REALFFTH_EDATA, REALFFTH_MSGEDATA );

  ASSERT( plan->sign == 1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  code = XLALREAL4ReverseFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( plan->sign != 1 ) /* plan sign was wrong */
        {
          ABORT( status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      case XLAL_EDOM: /* either DC or Nyquist imaginary part was non zero */
        ABORT( status, REALFFTH_EDATA, REALFFTH_MSGEDATA );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}



void
LALREAL4PowerSpectrum (
    LALStatus   *status,
    REAL4Vector *spec,
    REAL4Vector *data,
    REAL4FFTPlan *plan
    )
{
  int code;
  UINT4 n;

  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALREAL4PowerSpectrum", "XLALREAL4PowerSpectrum");

  ASSERT( spec, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( spec->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( data->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  /* Plan=0 Bugfix
    ASSERT( plan->plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  */

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( data->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( spec->length == n/2 + 1, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  code = XLALREAL4PowerSpectrum( spec, data, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}



void
LALREAL4VectorFFT(
    LALStatus   *status,
    REAL4Vector *output,
    REAL4Vector *input,
    REAL4FFTPlan *plan
    )
{
  int code;
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALREAL4VectorFFT", "XLALREAL4VectorFFT");

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  /* Plan=0 Bugfix
     ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  */

  /* make sure that it is not the same data! */
  ASSERT( output->data != input->data, status,
      REALFFTH_ESAME, REALFFTH_MSGESAME );

  ASSERT( plan->size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( output->length == plan->size, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->length == plan->size, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  code = XLALREAL4VectorFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        if ( ! plan->size ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( output->data == input->data ) /* same data pointers */
        {
          ABORT( status, REALFFTH_ESAME, REALFFTH_MSGESAME );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/*
 *
 * LAL REAL8 Functions
 *
 */



void
LALCreateForwardREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCreateForwardREAL8FFTPlan", "XLALCreateForwardREAL8FFTPlan");

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  *plan = XLALCreateREAL8FFTPlan( size, 1, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}



void
LALCreateReverseREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALCreateReverseREAL8FFTPlan", "XLALCreateReverseREAL8FFTPlan");

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  *plan = XLALCreateREAL8FFTPlan( size, 0, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}



void
LALDestroyREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan
    )
{
    INITSTATUS(status);
  XLALPrintDeprecationWarning("LALDestroyREAL8FFTPlan", "XLALDestroyREAL8FFTPlan");
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( *plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  XLALDestroyREAL8FFTPlan( *plan );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }
  *plan = NULL;
  RETURN( status );
}



void
LALForwardREAL8FFT(
    LALStatus      *status,
    COMPLEX16Vector *output,
    REAL8Vector    *input,
    REAL8FFTPlan    *plan
    )
{
  int code;
  UINT4 n;
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALForwardREAL8FFT", "XLALForwardREAL8FFT");

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  /* Plan=0 Bugfix
     ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  */

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( input->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( output->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  ASSERT( plan->sign == -1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  code = XLALREAL8ForwardFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( plan->sign != -1 ) /* plan sign was wrong */
        {
          ABORT( status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}



void
LALReverseREAL8FFT(
    LALStatus      *status,
    REAL8Vector    *output,
    COMPLEX16Vector *input,
    REAL8FFTPlan    *plan
    )
{
  int code;
  UINT4 n;
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALReverseREAL8FFT", "XLALReverseREAL8FFT");

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  /* Plan=0 Bugfix
     ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  */

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( output->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->data[0].im == 0, status, REALFFTH_EDATA, REALFFTH_MSGEDATA );
  ASSERT( n % 2 || input->data[n / 2].im == 0, status,
      REALFFTH_EDATA, REALFFTH_MSGEDATA );

  ASSERT( plan->sign == 1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  code = XLALREAL8ReverseFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( plan->sign != 1 ) /* plan sign was wrong */
        {
          ABORT( status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      case XLAL_EDOM: /* either DC or Nyquist imaginary part was non zero */
        ABORT( status, REALFFTH_EDATA, REALFFTH_MSGEDATA );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}



void
LALREAL8PowerSpectrum (
    LALStatus   *status,
    REAL8Vector *spec,
    REAL8Vector *data,
    REAL8FFTPlan *plan
    )
{
  int code;
  UINT4 n;

  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALREAL8PowerSpectrum", "XLALREAL8PowerSpectrum");

  ASSERT( spec, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( spec->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( data->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  /* Plan=0 Bugfix
     ASSERT( plan->plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  */

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( data->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( spec->length == n/2 + 1, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  code = XLALREAL8PowerSpectrum( spec, data, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}



void
LALREAL8VectorFFT(
    LALStatus   *status,
    REAL8Vector *output,
    REAL8Vector *input,
    REAL8FFTPlan *plan
    )
{
  int code;
  INITSTATUS(status);
  XLALPrintDeprecationWarning("LALREAL8VectorFFT", "XLALREAL8VectorFFT");

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  /* Plan=0 Bugfix
     ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  */

  /* make sure that it is not the same data! */
  ASSERT( output->data != input->data, status,
      REALFFTH_ESAME, REALFFTH_MSGESAME );

  ASSERT( plan->size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( output->length == plan->size, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->length == plan->size, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  code = XLALREAL8VectorFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        if ( ! plan->size ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( output->data == input->data ) /* same data pointers */
        {
          ABORT( status, REALFFTH_ESAME, REALFFTH_MSGESAME );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}
