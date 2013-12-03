/*
*  Copyright (C) 2007 Duncan Brown
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

/*-----------------------------------------------------------------------
 *
 * File Name: IntelComplexFFT.c
 *
 * Author: Brown D. A.
 *
 *
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <mkl_dfti.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>

#ifdef LAL_QTHREAD
extern int dummy_have_qthread;
#endif

struct
tagComplexFFTPlan
{
  INT4       sign;
  UINT4      size;
  DFTI_DESCRIPTOR *plan;
};


#define CHECKINTELFFTSTATUS( fstat )                            \
  if ( (fstat) != DFTI_NO_ERROR )                               \
  {                                                             \
    char *errmsg = DftiErrorMessage( fftStat );              \
    LALPrintError( errmsg );                                    \
    ABORT( status, COMPLEXFFTH_EINTL, COMPLEXFFTH_MSGEINTL );   \
  }                                                             \
  else (void)(0)


void
LALCreateForwardComplexFFTPlan(
    LALStatus       *status,
    ComplexFFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{
  INT8 fftStat;

#ifdef LAL_QTHREAD
  dummy_have_qthread = 0;
#endif

  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALCreateForwardComplexFFTPlan");

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  /* allocate memory */
  *plan = LALMalloc( sizeof( **plan ) );
  if ( ! *plan )
  {
    ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
  }

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;

  /* create complex intel fft descriptor */
  fftStat = DftiCreateDescriptor( &((*plan)->plan),
      DFTI_SINGLE, DFTI_COMPLEX, 1, size);
  CHECKINTELFFTSTATUS( fftStat );

  /* configure intel fft descriptor */
  fftStat = DftiSetValue( (*plan)->plan, DFTI_PLACEMENT,
      DFTI_NOT_INPLACE );
  CHECKINTELFFTSTATUS( fftStat );

  /* commit the intel fft descriptor */
  fftStat = DftiCommitDescriptor( (*plan)->plan );
  CHECKINTELFFTSTATUS( fftStat );

  RETURN( status );
}


void
LALCreateReverseComplexFFTPlan(
    LALStatus       *status,
    ComplexFFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{
  INT8 fftStat;

#ifdef LAL_QTHREAD
    dummy_have_qthread = 0;
#endif

  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALCreateReverseComplexFFTPlan");

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  /* allocate memory */
  *plan = LALMalloc( sizeof( **plan ) );
  if ( ! *plan )
  {
    ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
  }

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;

  /* create complex intel fft descriptor */
  fftStat = DftiCreateDescriptor( &((*plan)->plan),
      DFTI_SINGLE, DFTI_COMPLEX, 1, size);
  CHECKINTELFFTSTATUS( fftStat );

  /* configure intel fft descriptor */
  fftStat = DftiSetValue( (*plan)->plan, DFTI_PLACEMENT,
      DFTI_NOT_INPLACE );
  CHECKINTELFFTSTATUS( fftStat );

  /* commit the intel fft descriptor */
  fftStat = DftiCommitDescriptor( (*plan)->plan );
  CHECKINTELFFTSTATUS( fftStat );

  RETURN( status );
}


void
LALDestroyComplexFFTPlan (
    LALStatus       *status,
    ComplexFFTPlan **plan
    )
{
  INT8 fftStat;

#ifdef LAL_QTHREAD
    dummy_have_qthread = 0;
#endif

  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALDestroyComplexFFTPlan");

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( *plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* free intel fft descriptor and set to NULL pointer */
  fftStat = DftiFreeDescriptor( &((*plan)->plan) );
  CHECKINTELFFTSTATUS( fftStat );

  LALFree( *plan );
  *plan = NULL;

  RETURN( status );
}

void
LALCOMPLEX8VectorFFT (
    LALStatus      *status,
    COMPLEX8Vector *output,
    COMPLEX8Vector *input,
    ComplexFFTPlan *plan
    )
{
  INT8 fftStat;

#ifdef LAL_QTHREAD
    dummy_have_qthread = 0;
#endif

  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALCOMPLEX8VectorFFT");

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

  /* complex intel fft */
  if ( plan->sign == -1 )
  {
    fftStat = DftiComputeForward( plan->plan, input->data, output->data );
    CHECKINTELFFTSTATUS( fftStat );
  }
  else if ( plan->sign = 1 )
  {
    fftStat = DftiComputeBackward( plan->plan, input->data, output->data );
    CHECKINTELFFTSTATUS( fftStat );
  }
  else
  {
    ABORT( status, COMPLEXFFTH_ESIGN, COMPLEXFFTH_MSGESIGN );
  }

  RETURN( status );
}

/* double precision routines */

#undef CHECKINTELFFTSTATUS
