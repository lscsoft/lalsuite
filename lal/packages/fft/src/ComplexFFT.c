/**** <lalVerbatim file="ComplexFFTCV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Module \texttt{ComplexFFT.c}}
 * \label{ss:ComplexFFT.c}
 * 
 * Functions for performing complex FFTs.
 * 
 * \subsubsection*{Prototypes}
 * \vspace{0.1in}
 * \input{ComplexFFTCP}
 * \idx{LALCreateForwardComplexFFTPlan()}
 * \idx{LALCreateReverseComplexFFTPlan()}
 * \idx{LALDestroyComplexFFTPlan()}
 * \idx{LALCOMPLEX8VectorFFT()}
 * 
 * \subsubsection*{Description}
 * 
 * This package provides a LAL-style interface with the FFTW fast Fourier
 * transform package~\cite{fj:1998}.
 * 
 * The routines \texttt{LALCreateForwardComplexFFTPlan()} and
 * \texttt{LALCreateReverseComplexFFTPlan()} create plans for computing the
 * forward and reverse FFTs of a given size.  The optimum plan is either
 * estimated (reasonably fast) if the measure flag is zero, or measured (can be
 * time-consuming, but gives better performance) if the measure flag is
 * non-zero.  The routine \texttt{LALDestroyComplexFFTPlan()} destroys either
 * of these flavours of plans.
 * 
 * The routine \texttt{LALCOMPLEX8VectorFFT()} performs either the forward or
 * reverse FFT depending on the plan.  The discrete Fourier transform $H_k$,
 * $k=0\ldots n-1$ of a vector $h_j$, $j=0\ldots n-1$, of length $n$ is defined
 * by
 * \[
 *   H_k = \sum_{j=0}^{n-1} h_j e^{-2\pi ijk/n}
 * \]
 * and, similarly, the \emph{inverse} Fourier transform is defined by
 * \[
 *   h_j = \frac{1}{n}\sum_{k=0}^{n-1} H_k e^{2\pi ijk/n}.
 * \]
 * However, the present implementation of the \emph{reverse} FFT omits the
 * factor of $1/n$.  The input and output vectors must be distinct.
 * 
 * \subsubsection*{Operating Instructions}
 * 
 * \begin{verbatim}
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
 * \end{verbatim}
 * 
 * \subsubsection*{Algorithm}
 * 
 * The FFTW~\cite{fj:1998} is used.
 * 
 * \subsubsection*{Uses}
 * 
 * \subsubsection*{Notes}
 * 
 * \begin{enumerate}
 * \item The sign convention used here is the opposite of the definition in
 * \textit{Numerical Recipes}~\cite{ptvf:1992}, but agrees with the one used
 * by FFTW~\cite{fj:1998} and the other LIGO software components.
 * \item The result of the inverse FFT must be multiplied by $1/n$ to recover
 * the original vector.  This is different from the \texttt{datacondAPI} where
 * the factor is applied by default.
 * \item The size $n$ of the transform can be any positive integer; the
 * performance is $O(n\log n)$.  However, better performance is obtained if $n$
 * is the product of powers of 2, 3, 5, 7, and zero or one power of either 11
 * or 13.  Transforms when $n$ is a power of 2 are especially fast.  See
 * Ref.~\cite{fj:1998}.
 * \item LALMalloc() is used by all the fftw routines.
 * \item The input and output vectors for \texttt+LALCOMPLEX8VectorFFT()+ must
 * be distinct.
 * \end{enumerate}
 * 
 * \vfill{\footnotesize\input{ComplexFFTCV}}
 * 
 **** </lalLaTeX> */

#include <config.h>

#ifdef HAVE_SFFTW_H
#include <sfftw.h>
#elif HAVE_FFTW_H
#include <fftw.h>
#else
#error "don't have either sfftw.h or fftw.h"
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FFTWMutex.h>

NRCSID( COMPLEXFFTC, "$Id$" );

/* tell FFTW to use LALMalloc and LALFree */
#define FFTWHOOKS \
  do { fftw_malloc_hook = LALMallocShort; fftw_free_hook = LALFree; } while(0)

struct
tagComplexFFTPlan
{
  INT4      sign;
  UINT4     size;
  fftw_plan plan;
};

/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCreateForwardComplexFFTPlan(
    LALStatus       *status,
    ComplexFFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{ /* </lalVerbatim> */
  int flags = FFTW_THREADSAFE | ( measure ? FFTW_MEASURE : FFTW_ESTIMATE );
  INITSTATUS( status, "LALCreateForwardComplexFFTPlan", COMPLEXFFTC );
  FFTWHOOKS;

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

  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = fftw_create_plan( size, FFTW_FORWARD, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
  }

  RETURN( status );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCreateReverseComplexFFTPlan(
    LALStatus       *status,
    ComplexFFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{ /* </lalVerbatim> */
  int flags = FFTW_THREADSAFE | ( measure ? FFTW_MEASURE : FFTW_ESTIMATE );
  INITSTATUS( status, "LALCreateReverseComplexFFTPlan", COMPLEXFFTC );
  FFTWHOOKS;

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

  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = fftw_create_plan( size, FFTW_BACKWARD, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
  }

  RETURN( status );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALDestroyComplexFFTPlan (
    LALStatus       *status,
    ComplexFFTPlan **plan
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALDestroyComplexFFTPlan", COMPLEXFFTC );
  FFTWHOOKS;

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( *plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* destroy plan and set to NULL pointer */
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  fftw_destroy_plan( (*plan)->plan );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
  LALFree( *plan );
  *plan = NULL;

  RETURN( status );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCOMPLEX8VectorFFT (
    LALStatus      *status,
    COMPLEX8Vector *output,
    COMPLEX8Vector *input,
    ComplexFFTPlan *plan
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCOMPLEX8VectorFFT", COMPLEXFFTC );
  FFTWHOOKS;

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

  fftw_one(
      plan->plan,
      (fftw_complex *)input->data,
      (fftw_complex *)output->data
      );

  RETURN( status );
}

#ifdef KEEP_OLD_COMPLEX_FFT

void
LALEstimateFwdComplexFFTPlan (
    LALStatus       *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{
  INITSTATUS( stat, "LALEstimateFwdComplexFFTPlan", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* make sure that the plan has not been previously defined */
  ASSERT( *plan == NULL, stat, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );

  /* make sure that the requested size is valid */
  ASSERT( size > 0, stat, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  /* allocate memory */
  *plan = (ComplexFFTPlan *) LALMalloc( sizeof( ComplexFFTPlan ) );
  ASSERT( *plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = (void *)
    fftw_create_plan( size, 1, FFTW_THREADSAFE | FFTW_ESTIMATE );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN( stat );
}


void
LALEstimateInvComplexFFTPlan (
    LALStatus       *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{
  INITSTATUS( stat, "LALEstimateInvComplexFFTPlan", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* make sure that the plan has not been previously defined */
  ASSERT( *plan == NULL, stat, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );

  /* make sure that the requested size is valid */
  ASSERT( size > 0, stat, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  /* allocate memory */
  *plan = (ComplexFFTPlan *) LALMalloc( sizeof( ComplexFFTPlan ) );
  ASSERT( *plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = (void *)
    fftw_create_plan( size, -1, FFTW_THREADSAFE | FFTW_ESTIMATE );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN( stat );
}


void
LALMeasureFwdComplexFFTPlan (
    LALStatus       *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{
  INITSTATUS( stat, "LALMeasureFwdComplexFFTPlan", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* make sure that the plan has not been previously defined */
  ASSERT( *plan == NULL, stat, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );

  /* make sure that the requested size is valid */
  ASSERT( size > 0, stat, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  /* allocate memory */
  *plan = (ComplexFFTPlan *) LALMalloc( sizeof( ComplexFFTPlan ) );
  ASSERT( *plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = (void *)
    fftw_create_plan( size, 1, FFTW_THREADSAFE | FFTW_MEASURE );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN( stat );
}


void
LALMeasureInvComplexFFTPlan (
    LALStatus       *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{
  INITSTATUS( stat, "LALMeasureInvComplexFFTPlan", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* make sure that the plan has not been previously defined */
  ASSERT( *plan == NULL, stat, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );

  /* make sure that the requested size is valid */
  ASSERT( size > 0, stat, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  /* allocate memory */
  *plan = (ComplexFFTPlan *) LALMalloc( sizeof( ComplexFFTPlan ) );
  ASSERT( *plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = (void *)
    fftw_create_plan( size, -1, FFTW_THREADSAFE | FFTW_MEASURE );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN( stat );
}

#endif /* KEEP_OLD_COMPLEX_FFT */
