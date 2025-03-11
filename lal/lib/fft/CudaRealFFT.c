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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#include <string.h>

#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/RealFFT.h>
#include <lal/XLALError.h>
#include <lal/FFTWMutex.h>
#include <fftw3.h>

#include "CudaFunctions.h"
#include "CudaFFT.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

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


REAL4FFTPlan * XLALCreateREAL4FFTPlan( UINT4 size, int fwdflg, UNUSED int measurelvl )
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

  if( size == 1 ) createSize = 2;
  else	createSize = size;
  /* LAL_FFTW_WISDOM_LOCK; */
  if ( fwdflg ) /* forward */
    cufftPlan1d( &plan->plan, createSize, CUFFT_R2C, 1 );
  else /* reverse */
    cufftPlan1d( &plan->plan, createSize, CUFFT_C2R, 1 );
  /* LAL_FFTW_WISDOM_UNLOCK; */
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

  /* LAL_FFTW_WISDOM_LOCK; */
  /* Free the Cuda specific variables */
  XLALCudaFree( plan->d_real );
  cufftDestroy( plan->plan );
  XLALCudaFree( plan->d_complex );
  /* LAL_FFTW_WISDOM_UNLOCK; */
  memset( plan, 0, sizeof( *plan ) );
  XLALFree( plan );
  return;
}

int XLALREAL4ForwardFFT( COMPLEX8Vector *output, const REAL4Vector *input, const REAL4FFTPlan *plan )
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
    output->data[0] = crectf( input->data[0], 0.0 );
  }
  else
    cudafft_execute_r2c( plan->plan,
    (cufftComplex *)(output->data), (cufftReal *)(input->data),
    (cufftComplex *)plan->d_complex, (cufftReal *)plan->d_real, plan->size );

  /* Nyquist frequency */
  if( plan->size%2 == 0 )
  {
    output->data[plan->size/2] = crectf( crealf(output->data[plan->size/2]), 0.0 );
  }

  return 0;
}


int XLALREAL4ReverseFFT( REAL4Vector *output, const COMPLEX8Vector *input, const REAL4FFTPlan *plan )
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


int XLALREAL4VectorFFT( REAL4Vector * _LAL_RESTRICT_ output, const REAL4Vector * _LAL_RESTRICT_ input, const REAL4FFTPlan *plan )
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
      tmp[0] = crectf( input->data[0], 0.0 );

      for( k = 1; k < (plan->size + 1)/2; k++ )
      {
	tmp[k] = crectf( input->data[k], input->data[plan->size - k] );
      }

      if( plan->size%2 == 0 )
      {
	tmp[plan->size/2] = crectf( input->data[plan->size/2], 0.0 );
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
    tmp[0] = crectf( data->data[0], cimagf(tmp[0]) );
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

  LAL_FFTW_WISDOM_LOCK;
  if ( fwdflg ) /* forward */
    plan->plan = fftw_plan_r2r_1d( size, tmp1, tmp2, FFTW_R2HC, flags );
  else /* reverse */
    plan->plan = fftw_plan_r2r_1d( size, tmp1, tmp2, FFTW_HC2R, flags );
  LAL_FFTW_WISDOM_UNLOCK;

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
  LAL_FFTW_WISDOM_LOCK;
  fftw_destroy_plan( plan->plan );
  LAL_FFTW_WISDOM_UNLOCK;
  memset( plan, 0, sizeof( *plan ) );
  XLALFree( plan );

  return;
}

int XLALREAL8ForwardFFT( COMPLEX16Vector *output, const REAL8Vector *input, const REAL8FFTPlan *plan )
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
  output->data[0] = crect( tmp[0], 0.0 );

  /* other components */
  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size/2 rounded up */
  {
    output->data[k] = crect( tmp[k], tmp[plan->size - k] );
  }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* n is even */
  {
    output->data[plan->size/2] = crect( tmp[plan->size/2], 0.0 );
  }

  XLALFree( tmp );
  return 0;
}


int XLALREAL8ReverseFFT( REAL8Vector *output, const COMPLEX16Vector *input, const REAL8FFTPlan *plan )
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
  if ( cimag(input->data[0]) != 0.0 )
    XLAL_ERROR( XLAL_EDOM );  /* imaginary part of DC must be zero */
  if ( ! plan->size % 2 && cimag(input->data[plan->size/2]) != 0.0 )
    XLAL_ERROR( XLAL_EDOM );  /* imaginary part of Nyquist must be zero */

  /* create temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( XLAL_ENOMEM );

  /* unpack input into temporary array */

  /* dc component */
  tmp[0] = creal(input->data[0]);

  /* other components */
  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size / 2 rounded up */
  {
    tmp[k]              = creal(input->data[k]);
    tmp[plan->size - k] = cimag(input->data[k]);
  }

  /* Nyquist component */
  if ( plan->size%2 == 0 ) /* n is even */
    tmp[plan->size/2] = creal(input->data[plan->size/2]);

  /* perform the fft */
  fftw_execute_r2r( plan->plan, tmp, output->data );

  /* cleanup and exit */
  XLALFree( tmp );
  return 0;
}


int XLALREAL8VectorFFT( REAL8Vector * _LAL_RESTRICT_ output, const REAL8Vector * _LAL_RESTRICT_ input, const REAL8FFTPlan *plan )
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


int XLALREAL8PowerSpectrum( REAL8Vector *spec, const REAL8Vector *data,
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
