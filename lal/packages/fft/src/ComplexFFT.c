/*----------------------------------------------------------------------- 
 * 
 * File Name: ComplexFFT.c
 * 
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * ComplexFFT
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Creates plans for forward and inverse complex FFTs.
 * Performs foward and inverse complex FFTs on vectors.
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
#ifdef HAVE_SFFTW_H
#include <sfftw.h>
#elif HAVE_FFTW_H
#include <fftw.h>
#else
#error "don't have either sfftw.h or fftw.h"
#endif
*/

/* don't actually include sfftw.h or fftw.h because these are broken */
#ifndef HAVE_SFFTW_H
#ifndef HAVE_FFTW_H
#error "don't have either sfftw.h or fftw.h"
#endif
#endif

#include <pthread.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>

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

NRCSID( COMPLEXFFTC, "$Id$" );

static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;

/* tell FFTW to use LALMalloc and LALFree */
#define FFTWHOOKS \
  do { fftw_malloc_hook = LALMalloc; fftw_free_hook = LALFree; } while(0)


void
LALEstimateFwdComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{
  INITSTATUS( stat, "LALEstimateFwdComplexFFTPlan", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* make sure that the plan has not been previously defined */
  ASSERT( *plan == NULL, stat, COMPLEXFFT_ENNUL, COMPLEXFFT_MSGENNUL );

  /* make sure that the requested size is valid */
  ASSERT( size > 0, stat, COMPLEXFFT_ESIZE, COMPLEXFFT_MSGESIZE );

  /* allocate memory */
  *plan = (ComplexFFTPlan *) LALMalloc( sizeof( ComplexFFTPlan ) );
  ASSERT( *plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;
  pthread_mutex_lock( &mut );  
  (*plan)->plan = (void *)
    fftw_create_plan( size, 1, FFTW_THREADSAFE | FFTW_ESTIMATE );
  pthread_mutex_unlock( &mut );  

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );
  }

  /* normal exit */
  RETURN( stat );
}


void
LALEstimateInvComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{
  INITSTATUS( stat, "LALEstimateInvComplexFFTPlan", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* make sure that the plan has not been previously defined */
  ASSERT( *plan == NULL, stat, COMPLEXFFT_ENNUL, COMPLEXFFT_MSGENNUL );

  /* make sure that the requested size is valid */
  ASSERT( size > 0, stat, COMPLEXFFT_ESIZE, COMPLEXFFT_MSGESIZE );

  /* allocate memory */
  *plan = (ComplexFFTPlan *) LALMalloc( sizeof( ComplexFFTPlan ) );
  ASSERT( *plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;
  pthread_mutex_lock( &mut );  
  (*plan)->plan = (void *)
    fftw_create_plan( size, -1, FFTW_THREADSAFE | FFTW_ESTIMATE );
  pthread_mutex_unlock( &mut );  

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );
  }

  /* normal exit */
  RETURN( stat );
}


void
LALMeasureFwdComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{
  INITSTATUS( stat, "LALMeasureFwdComplexFFTPlan", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* make sure that the plan has not been previously defined */
  ASSERT( *plan == NULL, stat, COMPLEXFFT_ENNUL, COMPLEXFFT_MSGENNUL );

  /* make sure that the requested size is valid */
  ASSERT( size > 0, stat, COMPLEXFFT_ESIZE, COMPLEXFFT_MSGESIZE );

  /* allocate memory */
  *plan = (ComplexFFTPlan *) LALMalloc( sizeof( ComplexFFTPlan ) );
  ASSERT( *plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;
  pthread_mutex_lock( &mut );  
  (*plan)->plan = (void *)
    fftw_create_plan( size, 1, FFTW_THREADSAFE | FFTW_MEASURE );
  pthread_mutex_unlock( &mut );  

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );
  }

  /* normal exit */
  RETURN( stat );
}


void
LALMeasureInvComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{
  INITSTATUS( stat, "LALMeasureInvComplexFFTPlan", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* make sure that the plan has not been previously defined */
  ASSERT( *plan == NULL, stat, COMPLEXFFT_ENNUL, COMPLEXFFT_MSGENNUL );

  /* make sure that the requested size is valid */
  ASSERT( size > 0, stat, COMPLEXFFT_ESIZE, COMPLEXFFT_MSGESIZE );

  /* allocate memory */
  *plan = (ComplexFFTPlan *) LALMalloc( sizeof( ComplexFFTPlan ) );
  ASSERT( *plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;
  pthread_mutex_lock( &mut );  
  (*plan)->plan = (void *)
    fftw_create_plan( size, -1, FFTW_THREADSAFE | FFTW_MEASURE );
  pthread_mutex_unlock( &mut );

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );
  }

  /* normal exit */
  RETURN( stat );
}


void
LALDestroyComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan
    )
{
  INITSTATUS (stat, "LALDestroyComplexFFTPlan", COMPLEXFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* check that the plan is not NULL */
  ASSERT( *plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* destroy plan and set to NULL pointer */
  pthread_mutex_lock( &mut );  
  fftw_destroy_plan( (fftw_plan)(*plan)->plan );
  pthread_mutex_unlock( &mut );  
  LALFree( *plan );
  *plan = NULL;

  /* normal exit */
  RETURN( stat );
}


void
LALCOMPLEX8VectorFFT (
    LALStatus         *stat,
    COMPLEX8Vector *vout,
    COMPLEX8Vector *vinp,
    ComplexFFTPlan *plan
    )
{
  INITSTATUS( stat, "LALCOMPLEX8VectorFFT", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT( vout, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );
  ASSERT( vinp, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );
  ASSERT( plan, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* make sure that the data exists */
  ASSERT( vout->data, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );
  ASSERT( vinp->data, stat, COMPLEXFFT_ENULL, COMPLEXFFT_MSGENULL );

  /* make sure that it is not the same data! */
  ASSERT( vout->data != vinp->data, stat,
          COMPLEXFFT_ESAME, COMPLEXFFT_MSGESAME );

  /* make sure that the lengths agree */
  ASSERT( plan->size > 0, stat, COMPLEXFFT_ESIZE, COMPLEXFFT_MSGESIZE );
  ASSERT( vout->length == plan->size, stat,
          COMPLEXFFT_ESZMM, COMPLEXFFT_MSGESZMM );
  ASSERT( vinp->length == plan->size, stat,
          COMPLEXFFT_ESZMM, COMPLEXFFT_MSGESZMM );

  fftw_one(
      (fftw_plan) plan->plan,
      (fftw_complex *) vinp->data,
      (fftw_complex *) vout->data
      );

  /* normal exit */
  RETURN( stat );
}

