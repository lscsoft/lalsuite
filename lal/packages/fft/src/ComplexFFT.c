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

#if defined HAVE_LIBFFTW3F && defined HAVE_FFTW3_H
#define USE_FFTW3
#include <fftw3.h>
#elif defined HAVE_SFFTW_H
#include <sfftw.h>
#elif defined HAVE_FFTW_H
#include <fftw.h>
#ifndef FFTW_ENABLE_FLOAT
#error "included fftw.h is not for single-precision"
#endif
#else
#error "don't have sfftw.h, fftw.h, or fftw3.h"
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FFTWMutex.h>

NRCSID( COMPLEXFFTC, "$Id$" );

/* tell FFTW to use LALMalloc and LALFree */
#ifdef USE_FFTW3
#define FFTWHOOKS ((void)0)
#else
#define FFTWHOOKS \
  do { fftw_malloc_hook = LALMallocShort; fftw_free_hook = LALFree; } while(0)
#endif

struct
tagComplexFFTPlan
{
  INT4       sign;
  UINT4      size;
#ifdef USE_FFTW3
  fftwf_plan plan;
#else
  fftw_plan  plan;
#endif
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
#ifdef USE_FFTW3
  COMPLEX8 *tmp1;
  COMPLEX8 *tmp2;
  int flags = FFTW_UNALIGNED;
  if ( measure == 0 ) flags |= FFTW_ESTIMATE;
  else if ( measure == 1 ) flags |= FFTW_MEASURE;
  else if ( measure == 2 ) flags |= FFTW_PATIENT;
  else flags |= FFTW_EXHAUSTIVE;
#else
  int flags = FFTW_THREADSAFE | ( measure ? FFTW_MEASURE : FFTW_ESTIMATE );
#endif
  INITSTATUS( status, "LALCreateForwardComplexFFTPlan", COMPLEXFFTC );
  FFTWHOOKS;

#ifndef USE_FFTW3
  ASSERT( fftw_sizeof_fftw_real() == 4, status,
      COMPLEXFFTH_ESNGL, COMPLEXFFTH_MSGESNGL );
#endif
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

#ifdef USE_FFTW3
  tmp1 = LALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = LALMalloc( size * sizeof( *tmp2 ) );
  if ( !tmp1 || !tmp2 )
  {
    if ( tmp1 ) LALFree( tmp1 );
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
  }
#endif
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
#ifdef USE_FFTW3
  (*plan)->plan = fftwf_plan_dft_1d( size,
      (fftwf_complex *)tmp1, (fftwf_complex *)tmp2, FFTW_FORWARD, flags );
#else
  (*plan)->plan = fftw_create_plan( size, FFTW_FORWARD, flags );
#endif
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
#ifdef USE_FFTW3
  LALFree( tmp2 );
  LALFree( tmp1 );
#endif
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
#ifdef USE_FFTW3
  COMPLEX8 *tmp1;
  COMPLEX8 *tmp2;
  int flags = FFTW_UNALIGNED;
  if ( measure == 0 ) flags |= FFTW_ESTIMATE;
  else if ( measure == 1 ) flags |= FFTW_MEASURE;
  else if ( measure == 2 ) flags |= FFTW_PATIENT;
  else flags |= FFTW_EXHAUSTIVE;
#else
  int flags = FFTW_THREADSAFE | ( measure ? FFTW_MEASURE : FFTW_ESTIMATE );
#endif
  INITSTATUS( status, "LALCreateReverseComplexFFTPlan", COMPLEXFFTC );
  FFTWHOOKS;

#ifndef USE_FFTW3
  ASSERT( fftw_sizeof_fftw_real() == 4, status,
      COMPLEXFFTH_ESNGL, COMPLEXFFTH_MSGESNGL );
#endif
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

#ifdef USE_FFTW3
  tmp1 = LALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = LALMalloc( size * sizeof( *tmp2 ) );
  if ( !tmp1 || !tmp2 )
  {
    if ( tmp1 ) LALFree( tmp1 );
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
  }
#endif
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
#ifdef USE_FFTW3
  (*plan)->plan = fftwf_plan_dft_1d( size,
      (fftwf_complex *)tmp1, (fftwf_complex *)tmp2, FFTW_BACKWARD, flags );
#else
  (*plan)->plan = fftw_create_plan( size, FFTW_BACKWARD, flags );
#endif
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
#ifdef USE_FFTW3
  LALFree( tmp2 );
  LALFree( tmp1 );
#endif
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
#ifdef USE_FFTW3
  fftwf_destroy_plan( (*plan)->plan );
#else
  fftw_destroy_plan( (*plan)->plan );
#endif
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

#ifdef USE_FFTW3
  fftwf_execute_dft(
      plan->plan,
      (fftwf_complex *)input->data,
      (fftwf_complex *)output->data
      );
#else
  fftw_one(
      plan->plan,
      (fftw_complex *)input->data,
      (fftw_complex *)output->data
      );
#endif

  RETURN( status );
}

/* double precision routines if they are available */
#if defined HAVE_LIBFFTW3 && defined HAVE_FFTW3_H
#endif
