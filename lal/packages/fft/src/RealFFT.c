/*----------------------------------------------------------------------- 
 * 
 * File Name: RealFFT.c
 * 
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * RealFFT
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Creates plans for forward and inverse real FFTs.
 * Performs foward and inverse real FFTs on vectors and sequences of vectors.
 * 
 * DIAGNOSTICS 
 *
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>

/*
#ifdef HAVE_SRFFTW_H
#include <srfftw.h>
#elif HAVE_RFFTW_H
#include <rfftw.h>
#else
#error "don't have either srfftw.h or rfftw.h"
#endif
*/

/* don't actually include sfftw.h or fftw.h because these are broken */
#ifndef HAVE_SRFFTW_H
#ifndef HAVE_RFFTW_H
#error "don't have either srfftw.h or rfftw.h"
#endif
#endif

#include <pthread.h>

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>

/* here are the nearly equivalent fftw prototypes */
#ifndef FFTW_H
typedef REAL4 fftw_real;
typedef COMPLEX8 fftw_complex;
typedef void *fftw_plan;
typedef void *( *fftw_malloc_type_function )( size_t );
typedef void  ( *fftw_free_type_function )( void * );
fftw_plan fftw_create_plan( int, int, int );
void fftw_destroy_plan( fftw_plan );
void fftw_one( fftw_plan, fftw_complex *, fftw_complex * );
extern fftw_malloc_type_function fftw_malloc_hook;
extern fftw_free_type_function fftw_free_hook;
#define FFTW_ESTIMATE       (0)
#define FFTW_MEASURE        (1)
#define FFTW_OUT_OF_PLACE   (0)
#define FFTW_IN_PLACE       (8)
#define FFTW_USE_WISDOM    (16)
#define FFTW_THREADSAFE   (128)
#define FFTW_FORWARD       (-1)
#define FFTW_BACKWARD       (1)
#endif

#ifndef RFFTW_H
typedef fftw_plan rfftw_plan;
#define FFTW_REAL_TO_COMPLEX FFTW_FORWARD
#define FFTW_COMPLEX_TO_REAL FFTW_BACKWARD
rfftw_plan rfftw_create_plan( int, int, int );
void rfftw_destroy_plan( rfftw_plan );
void rfftw( rfftw_plan, int, fftw_real *, int, int, fftw_real *, int, int );
void rfftw_one( rfftw_plan, fftw_real *, fftw_real * );
#endif

NRCSID (REALFFTC, "$Id$");

static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;

/* tell FFTW to use LALMalloc and LALFree */
#define FFTWHOOKS \
  do { fftw_malloc_hook = LALMalloc; fftw_free_hook = LALFree; } while(0)

void
LALEstimateFwdRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{
  INITSTATUS (stat, "LALEstimateFwdRealFFTPlan", REALFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the plan has not been previously defined */
  ASSERT (*plan == NULL, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that the requested size is valid */
  ASSERT (size > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc (sizeof(RealFFTPlan));
  ASSERT (*plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;
  pthread_mutex_lock( &mut );  
  (*plan)->plan = (void *)
    rfftw_create_plan (
        size,
        FFTW_REAL_TO_COMPLEX,
        FFTW_THREADSAFE | FFTW_ESTIMATE /* | FFTW_USE_WISDOM */
        );
  pthread_mutex_unlock( &mut );  

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, REALFFT_ENULL, REALFFT_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

void
LALEstimateInvRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{
  INITSTATUS (stat, "LALEstimateInvRealFFTPlan", REALFFTC);
 
  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the plan has not been previously defined */
  ASSERT (*plan == NULL, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that the requested size is valid */
  ASSERT (size > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc (sizeof(RealFFTPlan));
  ASSERT (*plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;
  pthread_mutex_lock( &mut );
  (*plan)->plan = (void *)
    rfftw_create_plan (
        size,
        FFTW_COMPLEX_TO_REAL,
        FFTW_THREADSAFE | FFTW_ESTIMATE /* | FFTW_USE_WISDOM */
        );
  pthread_mutex_unlock( &mut );

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, REALFFT_ENULL, REALFFT_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

void
LALMeasureFwdRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{
  INITSTATUS (stat, "LALMeasureFwdRealFFTPlan", REALFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the plan has not been previously defined */
  ASSERT (*plan == NULL, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that the requested size is valid */
  ASSERT (size > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc (sizeof(RealFFTPlan));
  ASSERT (*plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;
  pthread_mutex_lock( &mut );  
  (*plan)->plan = (void *)
    rfftw_create_plan (
        size,
        FFTW_REAL_TO_COMPLEX,
        FFTW_THREADSAFE | FFTW_MEASURE /* | FFTW_USE_WISDOM */
        );
  pthread_mutex_unlock( &mut );  

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, REALFFT_ENULL, REALFFT_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

void
LALMeasureInvRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{
  INITSTATUS (stat, "LALMeasureInvRealFFTPlan", REALFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the plan has not been previously defined */
  ASSERT (*plan == NULL, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that the requested size is valid */
  ASSERT (size > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc (sizeof(RealFFTPlan));
  ASSERT (*plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;
  pthread_mutex_lock( &mut );  
  (*plan)->plan = (void *)
    rfftw_create_plan (
        size,
        FFTW_COMPLEX_TO_REAL,
        FFTW_THREADSAFE | FFTW_MEASURE /* | FFTW_USE_WISDOM */
        );
  pthread_mutex_unlock( &mut );  

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, REALFFT_ENULL, REALFFT_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

void
LALDestroyRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan
    )
{
  INITSTATUS (stat, "LALDestroyRealFFTPlan", REALFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* check that the plan is not NULL */
  ASSERT (*plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* destroy plan and set to NULL pointer */
  pthread_mutex_lock( &mut );  
  rfftw_destroy_plan ((rfftw_plan)(*plan)->plan);
  pthread_mutex_unlock( &mut );  
  LALFree (*plan);
  *plan = NULL;

  /* normal exit */
  RETURN (stat);
}


void
LALREAL4VectorFFT (
    LALStatus      *stat,
    REAL4Vector *vout,
    REAL4Vector *vinp,
    RealFFTPlan *plan
    )
{
  INITSTATUS (stat, "LALREAL4VectorFFT", REALFFTC);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the data exists */
  ASSERT (vout->data, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp->data, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that it is not the same data! */
  ASSERT (vout->data != vinp->data, stat, REALFFT_ESAME, REALFFT_MSGESAME);

  /* make sure that the lengths agree */
  ASSERT (plan->size > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);
  ASSERT (vout->length == plan->size, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vinp->length == plan->size, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);

  rfftw_one (
      (rfftw_plan) plan->plan,
      (fftw_real *)vinp->data,
      (fftw_real *)vout->data
      );

  /* normal exit */
  RETURN (stat);
}


void
LALREAL4VectorSequenceFFT (
    LALStatus              *stat,
    REAL4VectorSequence *vout,
    REAL4VectorSequence *vinp,
    RealFFTPlan         *plan
    )
{
  INITSTATUS (stat, "LALREAL4VectorSequenceFFT", REALFFTC);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the data exists */
  ASSERT (vout->data, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp->data, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that it is not the same data! */
  ASSERT (vout->data != vinp->data, stat, REALFFT_ESAME, REALFFT_MSGESAME);

  /* make sure that the lengths agree */
  ASSERT (plan->size > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);
  ASSERT (vout->vectorLength == plan->size, stat,
          REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vinp->vectorLength == plan->size, stat,
          REALFFT_ESZMM, REALFFT_MSGESZMM);

  ASSERT (vout->length > 0, stat, REALFFT_ESLEN, REALFFT_MSGESLEN);
  ASSERT (vinp->length == vout->length, stat, REALFFT_ESLEN, REALFFT_MSGESLEN);
  rfftw (
      (rfftw_plan) plan->plan, vout->length,
      (fftw_real *)vinp->data, 1, plan->size,
      (fftw_real *)vout->data, 1, plan->size
      );

  /* normal exit */
  RETURN (stat);
}


void
LALFwdRealFFT (
    LALStatus         *stat,
    COMPLEX8Vector *vout,
    REAL4Vector    *vinp,
    RealFFTPlan    *plan
    )
{
  REAL4Vector *vtmp = NULL;
  UINT4        n;
  UINT4        k;

  INITSTATUS (stat, "LALFwdRealFFT", REALFFTC);
  ATTATCHSTATUSPTR (stat);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that sizes agree */
  n = plan->size;
  ASSERT (n > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);
  ASSERT (vinp->length == n, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vout->length == n/2 + 1, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == 1, stat, REALFFT_ESIGN, REALFFT_MSGESIGN);

  /* create temporary vector and check that it was created */
  TRY (LALCreateVector (stat->statusPtr, &vtmp, n), stat);

  /* perform the FFT */
  TRY (LALREAL4VectorFFT (stat->statusPtr, vtmp, vinp, plan), stat);

  /* DC component */
  vout->data[0].re = vtmp->data[0];
  vout->data[0].im = 0;

  /* other components */
  for (k = 1; k < (n + 1)/2; ++k) /* k < n/2 rounded up */
  {
    vout->data[k].re =  vtmp->data[k];
    vout->data[k].im = -vtmp->data[n - k]; /* correct sign */
  }

  /* Nyquist frequency */
  if (n % 2 == 0) /* n is even */
  {
    vout->data[n/2].re = vtmp->data[n/2];
    vout->data[n/2].im = 0;
  }

  /* destroy temporary vector and check that it was destroyed */
  TRY (LALDestroyVector (stat->statusPtr, &vtmp), stat);

  DETATCHSTATUSPTR (stat);
  /* normal exit */
  RETURN (stat);
}


void
LALInvRealFFT (
    LALStatus         *stat,
    REAL4Vector    *vout,
    COMPLEX8Vector *vinp,
    RealFFTPlan    *plan
    )
{
  REAL4Vector *vtmp = NULL;
  UINT4        n;
  UINT4        k;

  INITSTATUS (stat, "LALInvRealFFT", REALFFTC);
  ATTATCHSTATUSPTR (stat);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that sizes agree */
  n = plan->size;
  ASSERT (n > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);
  ASSERT (vinp->length == n/2 + 1, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vout->length == n, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == -1, stat, REALFFT_ESIGN, REALFFT_MSGESIGN);

  /* create temporary vector and check that it was created */
  TRY (LALCreateVector (stat->statusPtr, &vtmp, n), stat);

  /* DC component */
  vtmp->data[0] = vinp->data[0].re;

  /* other components */
  for ( k = 1; k < (n + 1)/2; ++k ) /* k < n/2 rounded up */
  {
    vtmp->data[k]     =  vinp->data[k].re;
    vtmp->data[n - k] = -vinp->data[k].im;   /* correct sign */
  }

  /* Nyquist component */
  if ( n % 2 == 0 ) /* n is even */
  {
    vtmp->data[n/2] = vinp->data[n/2].re;

    /* make sure the Nyquist component is purely real */
    ASSERT (vinp->data[n/2].im == 0, stat, REALFFT_EDATA, REALFFT_MSGEDATA);
  }

  /* perform the FFT */
  TRY (LALREAL4VectorFFT (stat->statusPtr, vout, vtmp, plan), stat);

  /* destroy temporary vector and check that it was destroyed */
  TRY (LALDestroyVector (stat->statusPtr, &vtmp), stat);

  DETATCHSTATUSPTR (stat);
  /* normal exit */
  RETURN (stat);
}


void
LALRealPowerSpectrum (
    LALStatus      *stat,
    REAL4Vector *vout,
    REAL4Vector *vinp,
    RealFFTPlan *plan
    )
{
  REAL4Vector *vtmp = NULL;
  UINT4        n;
  UINT4        k;

  INITSTATUS (stat, "LALRealPowerSpectrum", REALFFTC);
  ATTATCHSTATUSPTR (stat);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that sizes agree */
  n = plan->size;
  ASSERT (n > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);
  ASSERT (vinp->length == n, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vout->length == n/2 + 1, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == 1, stat, REALFFT_ESIGN, REALFFT_MSGESIGN);

  /* create temporary vector and check that it was created */
  TRY (LALCreateVector (stat->statusPtr, &vtmp, n), stat);

  /* perform the FFT */
  TRY (LALREAL4VectorFFT (stat->statusPtr, vtmp, vinp, plan), stat);

  /* DC component */
  vout->data[0] = (vtmp->data[0])*(vtmp->data[0]);

  /* other components */
  for (k = 1; k < (n + 1)/2; ++k) /* k < n/2 rounded up */
  {
    REAL4 re = vtmp->data[k];
    REAL4 im = vtmp->data[n - k];
    vout->data[k] = re*re + im*im;
  }

  /* Nyquist frequency */
  if (n % 2 == 0) /* n is even */
    vout->data[n/2] = (vtmp->data[n/2])*(vtmp->data[n/2]);

  /* destroy temporary vector and check that it was destroyed */
  TRY (LALDestroyVector (stat->statusPtr, &vtmp), stat);

  DETATCHSTATUSPTR (stat);
  /* normal exit */
  RETURN (stat);
}


void
LALFwdRealSequenceFFT (
    LALStatus                 *stat,
    COMPLEX8VectorSequence *vout,
    REAL4VectorSequence    *vinp,
    RealFFTPlan            *plan
    )
{
  REAL4VectorSequence    *vtmp = NULL;
  CreateVectorSequenceIn  sqin; 

  UINT4 m;
  UINT4 n;
  UINT4 j;
  UINT4 k;

  INITSTATUS (stat, "LALFwdRealSequenceFFT", REALFFTC);
  ATTATCHSTATUSPTR (stat);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that sizes agree */
  m = sqin.length = vinp->length;
  n = sqin.vectorLength = plan->size;
  ASSERT (m > 0, stat, REALFFT_ESLEN, REALFFT_MSGESLEN);
  ASSERT (n > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);
  ASSERT (vinp->vectorLength == n, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vout->vectorLength == n/2 + 1, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vout->length == m, stat, REALFFT_ESLEN, REALFFT_MSGESLEN);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == 1, stat, REALFFT_ESIGN, REALFFT_MSGESIGN);

  /* create temporary vector sequence and check that it was created */
  TRY (LALCreateVectorSequence (stat->statusPtr, &vtmp, &sqin), stat);

  /* perform the FFT */
  TRY (LALREAL4VectorSequenceFFT (stat->statusPtr, vtmp, vinp, plan), stat);

  /* loop over transforms */
  for (j = 0; j < m; ++j)
  {
    COMPLEX8 *z = vout->data + j*(n/2 + 1);
    REAL4    *x = vtmp->data + j*n;

    /* DC component */
    z[0].re = x[0];
    z[0].im = 0;

    /* other components */
    for (k = 1; k < (n + 1)/2; ++k) /* k < n/2 rounded up */
    {
      z[k].re =  x[k];
      z[k].im = -x[n - k];              /* correct sign */
    }

    /* Nyquist frequency */
    if (n % 2 == 0) /* n is even */
    {
      z[n/2].re = x[n/2];
      z[n/2].im = 0;
    }
  }

  /* destroy temporary vector and check that it was destroyed */
  TRY (LALDestroyVectorSequence (stat->statusPtr, &vtmp), stat);

  DETATCHSTATUSPTR (stat);
  /* normal exit */
  RETURN (stat);
}


void
LALInvRealSequenceFFT (
    LALStatus                 *stat,
    REAL4VectorSequence    *vout,
    COMPLEX8VectorSequence *vinp,
    RealFFTPlan            *plan
    )
{
  REAL4VectorSequence *vtmp = NULL;
  CreateVectorSequenceIn  sqin; 

  UINT4 m;
  UINT4 n;
  UINT4 j;
  UINT4 k;

  INITSTATUS (stat, "LALInvRealSequenceFFT", REALFFTC);
  ATTATCHSTATUSPTR (stat);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that sizes agree */
  m = sqin.length = vinp->length;
  n = sqin.vectorLength = plan->size;
  ASSERT (m > 0, stat, REALFFT_ESLEN, REALFFT_MSGESLEN);
  ASSERT (n > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);
  ASSERT (vinp->vectorLength == n/2 + 1, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vout->vectorLength == n, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vout->length == m, stat, REALFFT_ESLEN, REALFFT_MSGESLEN);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == -1, stat, REALFFT_ESIGN, REALFFT_MSGESIGN);

  /* create temporary vector sequence and check that it was created */
  TRY (LALCreateVectorSequence (stat->statusPtr, &vtmp, &sqin), stat);

  /* loop over transforms */
  for (j = 0; j < m; ++j)
  {
    COMPLEX8 *z = vinp->data + j*(n/2 + 1);
    REAL4    *x = vtmp->data + j*n;

    /* DC component */
    x[0] = z[0].re;

    /* make sure the DC component is purely real */
    ASSERT (z[0].im == 0, stat, REALFFT_EDATA, REALFFT_MSGEDATA);

    /* other components */
    for (k = 1; k < (n + 1)/2; ++k) /* k < n/2 rounded up */
    {
      x[k]     =  z[k].re;
      x[n - k] = -z[k].im;              /* correct sign */
    }

    /* Nyquist component */
    if (n % 2 == 0) /* n is even */
    {
      x[n/2] = z[n/2].re;
      /* make sure Nyquist component is purely real */
      ASSERT (z[n/2].im == 0, stat, REALFFT_EDATA, REALFFT_MSGEDATA);
    }
  }

  /* perform the FFT */
  TRY (LALREAL4VectorSequenceFFT (stat->statusPtr, vout, vtmp, plan), stat);

  /* destroy temporary vector sequence and check that it was destroyed */
  TRY (LALDestroyVectorSequence (stat->statusPtr, &vtmp), stat);

  DETATCHSTATUSPTR (stat);
  /* normal exit */
  RETURN (stat);
}


void
LALRealSequencePowerSpectrum (
    LALStatus              *stat,
    REAL4VectorSequence *vout,
    REAL4VectorSequence *vinp,
    RealFFTPlan         *plan
    )
{
  REAL4VectorSequence    *vtmp = NULL;
  CreateVectorSequenceIn  sqin; 

  UINT4 m;
  UINT4 n;
  UINT4 j;
  UINT4 k;

  INITSTATUS (stat, "LALRealSequencePowerSpectrum", REALFFTC);
  ATTATCHSTATUSPTR (stat);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (vinp, stat, REALFFT_ENULL, REALFFT_MSGENULL);
  ASSERT (plan, stat, REALFFT_ENULL, REALFFT_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFT_ENNUL, REALFFT_MSGENNUL);

  /* make sure that sizes agree */
  m = sqin.length = vinp->length;
  n = sqin.vectorLength = plan->size;
  ASSERT (m > 0, stat, REALFFT_ESLEN, REALFFT_MSGESLEN);
  ASSERT (n > 0, stat, REALFFT_ESIZE, REALFFT_MSGESIZE);
  ASSERT (vinp->vectorLength == n, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vout->vectorLength == n/2 + 1, stat, REALFFT_ESZMM, REALFFT_MSGESZMM);
  ASSERT (vout->length == m, stat, REALFFT_ESLEN, REALFFT_MSGESLEN);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == 1, stat, REALFFT_ESIGN, REALFFT_MSGESIGN);

  /* create temporary vector sequence and check that it was created */
  TRY (LALCreateVectorSequence (stat->statusPtr, &vtmp, &sqin), stat);

  /* perform the FFT */
  TRY (LALREAL4VectorSequenceFFT (stat->statusPtr, vtmp, vinp, plan), stat);

  /* loop over transforms */
  for (j = 0; j < m; ++j)
  {
    REAL4 *s = vout->data + j*(n/2 + 1);
    REAL4 *x = vtmp->data + j*n;

    /* DC component */
    s[0] = x[0]*x[0];

    /* other components */
    for (k = 1; k < (n + 1)/2; ++k) /* k < n/2 rounded up */
    {
      REAL4 re = x[k];
      REAL4 im = x[n - k];
      s[k] = re*re + im*im;
    }

    /* Nyquist frequency */
    if (n % 2 == 0) /* n is even */
      s[n/2] = x[n/2]*x[n/2];
  }

  /* destroy temporary vector and check that it was destroyed */
  TRY (LALDestroyVectorSequence (stat->statusPtr, &vtmp), stat);

  DETATCHSTATUSPTR (stat);
  /* normal exit */
  RETURN (stat);
}

#undef FFTWHOOKS
