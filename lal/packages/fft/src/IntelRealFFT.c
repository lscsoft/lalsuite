/*----------------------------------------------------------------------- 
 * 
 * File Name: IntelRealFFT.c
 *
 * Author: Brown D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <mkl_dfti.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>

NRCSID( REALFFTC, "$Id$" );

extern int dummy_have_qthread;

struct
tagRealFFTPlan
{
  INT4       sign;
  UINT4      size;
  DFTI_DESCRIPTOR *plan;
  REAL4     *tmp;
};


#define CHECKINTELFFTSTATUS( fstat )                    \
  if ( (fstat) != DFTI_NO_ERROR )                       \
  {                                                     \
    char *errmsg = DftiErrorMessage( fftStat );         \
    LALPrintError( errmsg );                            \
    ABORT( status, REALFFTH_EINTL, REALFFTH_MSGEINTL ); \
  }                                                     \
  else (void)(0)


/* NOTE: this function declaration, and others in this file, are overridden by
 * a macros in RealFFT.h.  If you expect to use these functions, you'll first
 * need to sort out the header file so that calls to these functions are not
 * redirected by the macros.  The recommended course of action would be to
 * rename them "Real" --> "REAL4", following the names of the equivalent
 * functions in RealFFT.h */
void
LALCreateForwardRealFFTPlan(
    LALStatus    *status,
    RealFFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{
  INT8  fftStat;

#ifdef LAL_QTHREAD
  dummy_have_qthread = 0;
#endif

  INITSTATUS( status, "LALCreateForwardRealFFTPlan", REALFFTC );

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  /* allocate memory */
  *plan = LALMalloc( sizeof( **plan ) );
  if ( ! *plan )
  {
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;

  /* make the intel fft descriptor */
  fftStat = DftiCreateDescriptor( &((*plan)->plan),
      DFTI_SINGLE, DFTI_REAL, 1, size);
  CHECKINTELFFTSTATUS( fftStat );

  /* configure intel fft descriptor */
  fftStat = DftiSetValue( (*plan)->plan, DFTI_PLACEMENT, 
      DFTI_NOT_INPLACE );
  CHECKINTELFFTSTATUS( fftStat );
  fftStat = DftiSetValue( (*plan)->plan, DFTI_PACKED_FORMAT, 
      DFTI_PACK_FORMAT );
  CHECKINTELFFTSTATUS( fftStat );

  /* commit the intel fft descriptor */
  fftStat = DftiCommitDescriptor( (*plan)->plan );
  CHECKINTELFFTSTATUS( fftStat );

  /* create workspace to do the fft into */
  (*plan)->tmp = LALMalloc( size * sizeof( REAL4 ) );
  if ( ! (*plan)->tmp )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }

  RETURN( status );
}


void
LALCreateReverseRealFFTPlan(
    LALStatus    *status,
    RealFFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{
  INT8  fftStat;

#ifdef LAL_QTHREAD
  dummy_have_qthread = 0;
#endif

  INITSTATUS( status, "LALCreateReverseRealFFTPlan", REALFFTC );

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  /* allocate memory */
  *plan = LALMalloc( sizeof( **plan ) );
  if ( ! *plan )
  {
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;

  /* make the intel fft descriptor */
  fftStat = DftiCreateDescriptor( &((*plan)->plan),
      DFTI_SINGLE, DFTI_REAL, 1, size);
  CHECKINTELFFTSTATUS( fftStat );

  /* configure intel fft descriptor */
  fftStat = DftiSetValue( (*plan)->plan, DFTI_PLACEMENT, 
      DFTI_NOT_INPLACE );
  CHECKINTELFFTSTATUS( fftStat );
  fftStat = DftiSetValue( (*plan)->plan, DFTI_PACKED_FORMAT, 
      DFTI_PACK_FORMAT );
  CHECKINTELFFTSTATUS( fftStat );

  /* commit the intel fft descriptor */
  fftStat = DftiCommitDescriptor( (*plan)->plan );
  CHECKINTELFFTSTATUS( fftStat );

  /* create workspace to do the fft into */
  (*plan)->tmp = LALMalloc( size * sizeof( REAL4 ) );
  if ( ! (*plan)->tmp )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }

  RETURN( status );
}


void
LALDestroyRealFFTPlan(
    LALStatus    *status,
    RealFFTPlan **plan
    )
{
  INT8  fftStat;

#ifdef LAL_QTHREAD
  dummy_have_qthread = 0;
#endif

  INITSTATUS( status, "LALDestroyRealFFTPlan", REALFFTC );

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( *plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  /* destroy intel fft descriptor */
  fftStat = DftiFreeDescriptor( &((*plan)->plan) );
  CHECKINTELFFTSTATUS( fftStat );

  LALFree( (*plan)->tmp );

  LALFree( *plan );
  *plan = NULL;

  RETURN( status );
}


void
LALForwardRealFFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    REAL4Vector    *input,
    RealFFTPlan    *plan
    )
{
  INT8 fftStat;
  UINT4 n;
  UINT4 k;

#ifdef LAL_QTHREAD
  dummy_have_qthread = 0;
#endif

  INITSTATUS( status, "LALForwardRealFFT", REALFFTC );
  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->tmp, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( input->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( output->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  ASSERT( plan->sign == -1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  /* execute intel fft */
  fftStat = DftiComputeForward( plan->plan, input->data, plan->tmp );
  CHECKINTELFFTSTATUS( fftStat );

  /* dc component */
  output->data[0].re = plan->tmp[0];
  output->data[0].im = 0;

  /* other components */
  for ( k = 1; k < ( n + 1 ) / 2; ++k ) /* k < n/2 rounded up */
  {
    output->data[k].re = plan->tmp[2 * k - 1];
    output->data[k].im = plan->tmp[2 * k];
  }

  /* Nyquist frequency */
  if ( n % 2 == 0 ) /* n is even */
  {
    output->data[n / 2].re = plan->tmp[n - 1];
    output->data[n / 2].im = 0;
  }

  RETURN( status );
}


void
LALReverseRealFFT(
    LALStatus      *status,
    REAL4Vector    *output,
    COMPLEX8Vector *input,
    RealFFTPlan    *plan
    )
{
  INT8 fftStat;
  UINT4 n;
  UINT4 k;

#ifdef LAL_QTHREAD
  dummy_have_qthread = 0;
#endif

  INITSTATUS( status, "LALReverseRealFFT", REALFFTC );

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->tmp, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( input->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( output->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  ASSERT( plan->sign == 1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  /* make sure that Nyquist is purely real if number of points is even */
  ASSERT( n % 2 || input->data[n / 2].im == 0, status,
      REALFFTH_EDATA, REALFFTH_MSGEDATA );

  /* dc component */
  plan->tmp[0] = input->data[0].re;

  /* other components */
  for ( k = 1; k < ( n + 1 ) / 2; ++k ) /* k < n / 2 rounded up */
  {
    plan->tmp[2 * k - 1] = input->data[k].re;
    plan->tmp[2 * k]     = input->data[k].im;
  }

  /* Nyquist component */
  if ( n % 2 == 0 ) /* n is even */
  {
    plan->tmp[n - 1] = input->data[n / 2].re;
  }

  /* execute intel fft */
  fftStat = DftiComputeBackward( plan->plan, plan->tmp, output->data );
  CHECKINTELFFTSTATUS( fftStat );

  RETURN( status );
}


void
LALRealPowerSpectrum (
    LALStatus   *status,
    REAL4Vector *spec,
    REAL4Vector *data,
    RealFFTPlan *plan
    )
{
  INT8  fftStat;
  UINT4 n;
  UINT4 k;

#ifdef LAL_QTHREAD
  dummy_have_qthread = 0;
#endif

  INITSTATUS( status, "LALRealPowerSpectrum", REALFFTC );
  
  ASSERT( spec, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( spec->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( data->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( plan->plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( plan->tmp, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( data->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( spec->length == n/2 + 1, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  /* execute intel fft */
  fftStat = DftiComputeForward( plan->plan, data->data, plan->tmp );
  CHECKINTELFFTSTATUS( fftStat );

  /* dc component */
  spec->data[0] = plan->tmp[0] * plan->tmp[0];

  /* other components */
  for (k = 1; k < (n + 1)/2; ++k) /* k < n/2 rounded up */
  {
    REAL4 re = plan->tmp[2 * k - 1];
    REAL4 im = plan->tmp[2 * k];
    spec->data[k] = re * re + im * im;
  }

  /* Nyquist frequency */
  if (n % 2 == 0) /* n is even */
    spec->data[n / 2] = plan->tmp[n / 2] * plan->tmp[n / 2];

  RETURN( status );
}


void
LALREAL4VectorFFT(
    LALStatus   *status,
    REAL4Vector *output,
    REAL4Vector *input,
    RealFFTPlan *plan
    )
{
  INT8 fftStat;
  UINT4 n;
  UINT4 k;

#ifdef LAL_QTHREAD
  dummy_have_qthread = 0;
#endif

  INITSTATUS( status, "LALREAL4VectorFFT", REALFFTC );

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->tmp, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  /* make sure that it is not the same data! */
  ASSERT( output->data != input->data, status,
      REALFFTH_ESAME, REALFFTH_MSGESAME );

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( output->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  /* complex intel fft */
  if ( plan->sign == -1 )
  {
    fftStat = DftiComputeForward( plan->plan, input->data, plan->tmp );
    CHECKINTELFFTSTATUS( fftStat );

    output->data[0] = plan->tmp[0];
    for ( k = 1 ; k < (n + 1) / 2 ; ++k ) /* k < n/2 rounded up */
    {
      output->data[k]     = plan->tmp[2 * k - 1];
      output->data[n - k] = plan->tmp[2 * k];
    }
    if ( n % 2 == 0) /* n is even */
    {
      output->data[n / 2] = plan->tmp[n - 1];
    }
  }
  else if ( plan->sign = 1 )
  {
    plan->tmp[0] = input->data[0];
    for ( k = 1 ; k < (n + 1) / 2 ; ++k ) /* k < n/2 rounded up */
    {
      plan->tmp[2 * k - 1] = input->data[k];
      plan->tmp[2 * k]     = input->data[n - k];
    }
    if ( n % 2 == 0) /* n is even */
    {
      plan->tmp[n - 1] = input->data[n / 2];
    }

    fftStat = DftiComputeBackward( plan->plan, plan->tmp, output->data );
    CHECKINTELFFTSTATUS( fftStat );
  }
  else
  {
    ABORT( status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );
  }

  RETURN( status );
}

/* double precision routines */

#undef CHECKINTELFFTSTATUS
