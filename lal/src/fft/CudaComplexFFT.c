/*
*  Copyright (C) 2010 Karsten Wiesner, Jolien Creighton
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

#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/XLALError.h>
#include <lal/ComplexFFT.h>
#include <lal/CudaPlan.h>
#include <lal/FFTWMutex.h>
#include <fftw3.h>

#include "CudaFunctions.h"
#include "CudaFFT.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*
 *
 * XLAL COMPLEX8 functions
 *
 */


COMPLEX8FFTPlan * XLALCreateCOMPLEX8FFTPlan( UINT4 size, int fwdflg, UNUSED int measurelvl )
{
  COMPLEX8FFTPlan *plan;
  UINT4 createSize;

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

  /* create the plan */
  /* LAL_FFTW_WISDOM_LOCK; */
  /*
   * Change the size to avoid CUDA bug with FFT of size 1
   */
  if( size == 1 ) createSize = 2;
  else createSize = size;
  cufftPlan1d( &plan->plan, createSize, CUFFT_C2C, 1 );
  /* LAL_FFTW_WISDOM_UNLOCK; */

  /* "Plan=0" Bugfix by Wiesner, K.: plan->plan is an integer handle not a pointer and 0 is a valid handle
      So checking against 0 and occasionaly destroy the plan is a bug.
    if ( ! plan->plan )
   {
    XLALFree( plan );
    XLAL_ERROR_NULL( XLAL_EFAILED );
  }
  */

  plan->d_input = XLALCudaMallocComplex(size);
  plan->d_output = XLALCudaMallocComplex(size);
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
  if ( ! plan )
    XLAL_ERROR_VOID( XLAL_EFAULT );
  /* Plan=0 Bugfix
   if ( ! plan->plan )
      XLAL_ERROR_VOID( XLAL_EINVAL );
  */
  //LAL_FFTW_WISDOM_LOCK;
  cufftDestroy( plan->plan );
  XLALCudaFree(plan->d_input);
  XLALCudaFree(plan->d_output);
  //LAL_FFTW_WISDOM_UNLOCK;
  memset( plan, 0, sizeof( *plan ) );
  XLALFree( plan );
  return;
}


int XLALCOMPLEX8VectorFFT( COMPLEX8Vector * _LAL_RESTRICT_ output, const COMPLEX8Vector * _LAL_RESTRICT_ input, const COMPLEX8FFTPlan *plan )
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
  if( plan->size == 1 )
  {
    output->data[0] = input->data[0];
  }
  else
    cudafft_execute_c2c( plan->plan,
	(cufftComplex *)(output->data), (cufftComplex *)(input->data),
	(cufftComplex *)plan->d_output, (cufftComplex *)plan->d_input,
	plan->sign, plan->size );
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
  LAL_FFTW_WISDOM_LOCK;
  plan->plan = fftw_plan_dft_1d( size,
      (fftw_complex *)tmp1, (fftw_complex *)tmp2,
      fwdflg ? FFTW_FORWARD : FFTW_BACKWARD, flags );
  LAL_FFTW_WISDOM_UNLOCK;

  /* free temporary arrays */
  XLALFree( tmp2 );
  XLALFree( tmp1 );

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
  if ( ! plan )
    XLAL_ERROR_VOID( XLAL_EFAULT );
  /*Plan=0 Bugfix
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


int XLALCOMPLEX16VectorFFT( COMPLEX16Vector * _LAL_RESTRICT_ output, const COMPLEX16Vector * _LAL_RESTRICT_ input, const COMPLEX16FFTPlan *plan )
{
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

  /* do the fft */
  fftw_execute_dft(
      plan->plan,
      (fftw_complex *)input->data,
      (fftw_complex *)output->data
      );
  return 0;
}
