/**** <lalVerbatim file="RealFFTCV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Module \texttt{RealFFT.c}}
 * \label{ss:RealFFT.c}
 * 
 * Functions for performing real FFTs.
 * 
 * \subsubsection*{Prototypes}
 * \vspace{0.1in}
 * \input{RealFFTCP}
 * \idx{LALCreateForwardRealFFTPlan()}
 * \idx{LALCreateReverseRealFFTPlan()}
 * \idx{LALDestroyRealFFTPlan()}
 * \idx{LALForwardRealFFT()}
 * \idx{LALReverseRealFFT()}
 * \idx{LALRealPowerSpectrum()}
 * \idx{LALREAL4VectorFFT()}
 * 
 * \subsubsection*{Description}
 * 
 * This package provides a LAL-style interface with the FFTW fast Fourier
 * transform package~\cite{fj:1998}.
 * 
 * The routines \texttt{LALCreateForwardRealFFTPlan()} and
 * \texttt{LALCreateReverseRealFFTPlan()} create plans for computing the
 * forward (real-to-complex) and reverse (complex-to-real) FFTs of a specified
 * size.  The optimum plan is either estimated (reasonably fast) if the measure
 * flag is zero, or measured (can be time-consuming, but gives better
 * performance) if the measure flag is non-zero.  The routine
 * \texttt{LALDestroyRealFFTPlan()} destroys any of these flavours of plans.
 * 
 * The routines \texttt{LALForwardRealFFT()} and \texttt{LALReverseRealFFT()}
 * perform the forward (real-to-complex) and reverse (complex-to-real) FFTs
 * using the plans.  The discrete Fourier transform $H_k$,
 * $k=0\ldots\lfloor{n/2}\rfloor$ ($n/2$ rounded down), of a vector $h_j$,
 * $j=0\ldots n-1$, of length $n$ is defined by
 * \[
 *   H_k = \sum_{j=0}^{n-1} h_j e^{-2\pi ijk/n}
 * \]
 * and, similarly, the \emph{inverse} Fourier transform is defined by
 * \[
 *   h_j = \frac{1}{n} \sum_{k=0}^{n-1} H_k e^{2\pi ijk/n}
 * \]
 * where $H_k$ for $\lfloor{n/2}\rfloor<k<n$ can be obtained from the relation
 * $H_k=H_{n-k}^\ast$.  The present implementation of the \emph{reverse} FFT
 * omits the factor of $1/n$.
 * 
 * The routines in this package require that the vector $h_j$, $j=0\ldots n-1$
 * be real; consequently, $H_k=H_{n-k}^\ast$ ($0\le k\le\lfloor n/2\rfloor$),
 * i.e., the negative frequency Fourier components are the complex conjugate of
 * the positive frequency Fourier components when the data is real.  Therefore,
 * one need compute and store only the first $\lfloor n/2\rfloor+1$ components
 * of $H_k$; only the values of $H_k$ for $k=0\ldots \lfloor n/2\rfloor$ are
 * returned (integer division is rounded down, e.g., $\lfloor 7/2\rfloor=3$).
 * 
 * The routine \texttt{LALRealPowerSpectrum()} computes the power spectrum
 * $P_k=|H_k|^2$, $k=0\ldots \lfloor n/2\rfloor$, of the data $h_j$,
 * $j=0\ldots n-1$.
 * 
 * The routine \texttt{LALREAL4VectorFFT()} is essentially a direct calls to
 * FFTW routines without any re-packing of the data.  This routine should not
 * be used unless the user understands the packing used in FFTW.
 * 
 * \subsubsection*{Operating Instructions}
 * 
 * \begin{verbatim}
 * const UINT4 n = 32;
 * static LALStatus status; 
 * RealFFTPlan            *pfwd = NULL;
 * RealFFTPlan            *prev = NULL;
 * REAL4Vector            *hvec = NULL;
 * COMPLEX8Vector         *Hvec = NULL;
 * REAL4Vector            *Pvec = NULL;
 * 
 * LALCreateForwardRealFFTPlan( &status, &pfwd, n );
 * LALCreateReverseRealFFTPlan( &status, &prev, n );
 * LALSCreateVector( &status, &hvec, n );
 * LALCCreateVector( &status, &Hvec, n/2 + 1 );
 * LALSCreateVector( &status, &Pvec, n/2 + 1 );
 * 
 * <assign data>
 *
 * LALRealPowerSpectrum( &status, Pvec, hvec, pfwd );
 * LALForwardRealFFT( &status, Hvec, hvec, pfwd );
 * LALReverseRealFFT( &status, hvec, Hvec, pinv );
 * 
 * LALDestroyRealFFTPlan( &status, &pfwd );
 * LALDestroyRealFFTPlan( &status, &pinv );
 * LALSDestroyVector( &status, &hvec );
 * LALCDestroyVector( &status, &Hvec );
 * LALSDestroyVector( &status, &Pvec );
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
 * \item The sign convention used here is the opposite of
 * \textit{Numerical Recipes}~\cite{ptvf:1992}, but agrees with the one used
 * by FFTW~\cite{fj:1998} and the other LIGO software components.
 * \item The result of the reverse FFT must be multiplied by $1/n$ to recover
 * the original vector.  This is unlike the \textit{Numerical
 * Recipes}~\cite{ptvf:1992} convension where the factor is $2/n$ for real
 * FFTs.  This is different from the \texttt{datacondAPI} where the
 * normalization constant is applied by default.
 * \item The size $n$ of the transform can be any positive integer; the
 * performance is $O(n\log n)$.  However, better performance is obtained if $n$
 * is the product of powers of 2, 3, 5, 7, and zero or one power of either 11
 * or 13.  Transforms when $n$ is a power of 2 are especially fast.  See
 * Ref.~\cite{fj:1998}.
 * \item All of these routines leave the input array undamaged.  (Except for
 * \verb+LALREAL4VectorFFT+.)
 * \item LALMalloc() is used by all the fftw routines.
 * \end{enumerate}
 * 
 * \vfill{\footnotesize\input{RealFFTCV}}
 * 
 **** </lalLaTeX> */


#include <config.h>

#if defined HAVE_LIBFFTW3F && defined HAVE_FFTW3_H
#include <fftw3.h>
#define USE_FFTW3
#elif defined HAVE_SRFFTW_H
#include <srfftw.h>
#elif defined HAVE_RFFTW_H
#include <rfftw.h>
#ifndef FFTW_ENABLE_FLOAT
#error "included rfftw.h is not for single-precision"
#endif
#else
#error "don't have srfftw.h, rfftw.h, or fftw3.h"
#endif

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/FFTWMutex.h>

/* tell FFTW to use LALMalloc and LALFree */
#ifdef USE_FFTW3
#define FFTWHOOKS ((void)0)
#else
#define FFTWHOOKS \
  do { fftw_malloc_hook = LALMallocShort; fftw_free_hook = LALFree; } while(0)
#endif

NRCSID( REALFFTC, "$Id$" );


struct
tagRealFFTPlan
{
  INT4       sign;
  UINT4      size;
#ifdef USE_FFTW3
  fftwf_plan  plan;
#else
  rfftw_plan plan;
#endif
};


/* <lalVerbatim file="RealFFTCP"> */  
void
LALCreateForwardRealFFTPlan(
    LALStatus    *status,
    RealFFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{ /* </lalVerbatim> */
#ifdef USE_FFTW3
  REAL4 *tmp1;
  REAL4 *tmp2;
  int flags = FFTW_UNALIGNED;
  if ( measure == 0 ) flags |= FFTW_ESTIMATE;
  else if ( measure == 1 ) flags |= FFTW_MEASURE;
  else if ( measure == 2 ) flags |= FFTW_PATIENT;
  else flags |= FFTW_EXHAUSTIVE;
#else
  int flags = FFTW_THREADSAFE | ( measure ? FFTW_MEASURE : FFTW_ESTIMATE );
#endif
  INITSTATUS( status, "LALCreateForwardRealFFTPlan", REALFFTC );
  FFTWHOOKS;

#ifndef USE_FFTW3
  ASSERT( fftw_sizeof_fftw_real() == 4, status,
      REALFFTH_ESNGL, REALFFTH_MSGESNGL );  
#endif
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

#ifdef USE_FFTW3
  tmp1 = LALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = LALMalloc( size * sizeof( *tmp2 ));
  if ( !tmp1 || !tmp2 )
  {
    if ( tmp1 ) LALFree( tmp1 );
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }
#endif
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
#ifdef USE_FFTW3
  (*plan)->plan = fftwf_plan_r2r_1d( size, tmp1, tmp2, FFTW_R2HC, flags );
#else
  (*plan)->plan = rfftw_create_plan( size, FFTW_REAL_TO_COMPLEX, flags );
#endif
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
#ifdef USE_FFTW3
  LALFree( tmp2 );
  LALFree( tmp1 );
#endif
  if ( ! (*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */  
void
LALCreateReverseRealFFTPlan(
    LALStatus    *status,
    RealFFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{ /* </lalVerbatim> */
#ifdef USE_FFTW3
  REAL4 *tmp1;
  REAL4 *tmp2;
  int flags = FFTW_UNALIGNED;
  if ( measure == 0 ) flags |= FFTW_ESTIMATE;
  else if ( measure == 1 ) flags |= FFTW_MEASURE;
  else if ( measure == 2 ) flags |= FFTW_PATIENT;
  else flags |= FFTW_EXHAUSTIVE;
#else
  int flags = FFTW_THREADSAFE | ( measure ? FFTW_MEASURE : FFTW_ESTIMATE );
#endif
  INITSTATUS( status, "LALCreateReverseRealFFTPlan", REALFFTC );
  FFTWHOOKS;

#ifndef USE_FFTW3
  ASSERT( fftw_sizeof_fftw_real() == 4, status,
      REALFFTH_ESNGL, REALFFTH_MSGESNGL );  
#endif
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

#ifdef USE_FFTW3
  tmp1 = LALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = LALMalloc( size * sizeof( *tmp2 ) );
  if ( !tmp1 || !tmp2 )
  {
    if ( tmp1 ) LALFree( tmp1 );
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }
#endif
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
#ifdef USE_FFTW3
  (*plan)->plan = fftwf_plan_r2r_1d( size, tmp1, tmp2, FFTW_HC2R, flags );
#else
  (*plan)->plan = rfftw_create_plan( size, FFTW_COMPLEX_TO_REAL, flags );
#endif
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
#ifdef USE_FFTW3
  LALFree( tmp2 );
  LALFree( tmp1 );
#endif
  if ( ! (*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
  }

  RETURN( status );
}
  

/* <lalVerbatim file="RealFFTCP"> */  
void
LALDestroyRealFFTPlan(
    LALStatus    *status,
    RealFFTPlan **plan
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALDestroyRealFFTPlan", REALFFTC );
  FFTWHOOKS;
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( *plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  /* destroy plan and set to NULL pointer */
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
#ifdef USE_FFTW3
  fftwf_destroy_plan( (*plan)->plan );
#else
  rfftw_destroy_plan( (*plan)->plan );
#endif
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
  LALFree( *plan );
  *plan = NULL;

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */  
void
LALForwardRealFFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    REAL4Vector    *input,
    RealFFTPlan    *plan
    )
{ /* </lalVerbatim> */
  REAL4 *tmp;
  UINT4 n;
  UINT4 k;

  INITSTATUS( status, "LALForwardRealFFT", REALFFTC );

  FFTWHOOKS;

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( input->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( output->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  ASSERT( plan->sign == -1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  tmp = LALMalloc( n * sizeof( *tmp ) );
  if ( ! tmp )
  {
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }

#ifdef USE_FFTW3
  fftwf_execute_r2r( plan->plan, input->data, tmp );
#else
  rfftw_one( plan->plan, (fftw_real *) input->data, tmp );
#endif

  /* dc component */
  output->data[0].re = tmp[0];
  output->data[0].im = 0;

  /* other components */
  for ( k = 1; k < ( n + 1 ) / 2; ++k ) /* k < n/2 rounded up */
  {
    output->data[k].re = tmp[k];
    output->data[k].im = tmp[n - k];
  }

  /* Nyquist frequency */
  if ( n % 2 == 0 ) /* n is even */
  {
    output->data[n / 2].re = tmp[n / 2];
    output->data[n / 2].im = 0;
  }

  LALFree( tmp );
  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */  
void
LALReverseRealFFT(
    LALStatus      *status,
    REAL4Vector    *output,
    COMPLEX8Vector *input,
    RealFFTPlan    *plan
    )
{ /* </lalVerbatim> */
  REAL4 *tmp;
  UINT4 n;
  UINT4 k;

  INITSTATUS( status, "LALReverseRealFFT", REALFFTC );
  FFTWHOOKS;

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( input->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( output->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  ASSERT( plan->sign == 1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  /* make sure that Nyquist is purely real if number of points is even */
  ASSERT( n % 2 || input->data[n / 2].im == 0, status,
      REALFFTH_EDATA, REALFFTH_MSGEDATA );

  tmp = LALMalloc( n * sizeof( *tmp ) );
  if ( ! tmp )
  {
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }

  /* dc component */
  tmp[0] = input->data[0].re;

  /* other components */
  for ( k = 1; k < ( n + 1 ) / 2; ++k ) /* k < n / 2 rounded up */
  {
    tmp[k]     = input->data[k].re;
    tmp[n - k] = input->data[k].im;
  }

  /* Nyquist component */
  if ( n % 2 == 0 ) /* n is even */
  {
    tmp[n / 2] = input->data[n / 2].re;
  }

#ifdef USE_FFTW3
  fftwf_execute_r2r( plan->plan, tmp, output->data );
#else
  rfftw_one( plan->plan, tmp, (fftw_real *) output->data );
#endif

  LALFree( tmp );
  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */  
void
LALRealPowerSpectrum (
    LALStatus   *status,
    REAL4Vector *spec,
    REAL4Vector *data,
    RealFFTPlan *plan
    )
{ /* </lalVerbatim> */
  REAL4 *tmp;
  UINT4 n;
  UINT4 k;

  INITSTATUS( status, "LALRealPowerSpectrum", REALFFTC );

  FFTWHOOKS;

  ASSERT( spec, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( spec->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( data->data, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( plan->plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( data->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( spec->length == n/2 + 1, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  tmp = LALMalloc( n * sizeof( *tmp ) );
  if ( ! tmp )
  {
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }

#ifdef USE_FFTW3
  fftwf_execute_r2r( plan->plan, data->data, tmp );
#else
  rfftw_one( plan->plan, (fftw_real *) data->data, tmp );
#endif

  /* dc component */
  spec->data[0] = tmp[0] * tmp[0];

  /* other components */
  for (k = 1; k < (n + 1)/2; ++k) /* k < n/2 rounded up */
  {
    REAL4 re = tmp[k];
    REAL4 im = tmp[n - k];
    spec->data[k] = re * re + im * im;
  }

  /* Nyquist frequency */
  if (n % 2 == 0) /* n is even */
    spec->data[n / 2] = tmp[n / 2] * tmp[n / 2];

  LALFree( tmp );
  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */  
void
LALREAL4VectorFFT(
    LALStatus   *status,
    REAL4Vector *output,
    REAL4Vector *input,
    RealFFTPlan *plan
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALREAL4VectorFFT", REALFFTC );

  FFTWHOOKS;

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  /* make sure that it is not the same data! */
  ASSERT( output->data != input->data, status,
      REALFFTH_ESAME, REALFFTH_MSGESAME );

  ASSERT( plan->size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( output->length == plan->size, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->length == plan->size, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );

#ifdef USE_FFTW3
  fftwf_execute_r2r( plan->plan, input->data, output->data );
#else
  rfftw_one( plan->plan, (fftw_real *)input->data, (fftw_real *)output->data );
#endif

  RETURN( status );
}

#undef FFTWHOOKS

/* double precision routines if they are available */
#if defined HAVE_LIBFFTW3 && defined HAVE_FFTW3_H
#endif
