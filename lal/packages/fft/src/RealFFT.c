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

#ifdef HAVE_SRFFTW_H
#include <srfftw.h>
#elif HAVE_RFFTW_H
#include <rfftw.h>
#else
#error "don't have either srfftw.h or rfftw.h"
#endif

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/FFTWMutex.h>

/* tell FFTW to use LALMalloc and LALFree */
#define FFTWHOOKS \
  do { fftw_malloc_hook = LALMalloc; fftw_free_hook = LALFree; } while(0)

NRCSID( REALFFTC, "$Id$" );


struct
tagRealFFTPlan
{
  INT4       sign;
  UINT4      size;
  rfftw_plan plan;
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
  int flags = FFTW_THREADSAFE | ( measure ? FFTW_MEASURE : FFTW_ESTIMATE );
  INITSTATUS( status, "LALCreateForwardRealFFTPlan", REALFFTC );
  FFTWHOOKS;

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc( sizeof( **plan ) );
  if ( ! *plan )
  {
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;

  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = rfftw_create_plan( size, FFTW_REAL_TO_COMPLEX, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
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
  int flags = FFTW_THREADSAFE | ( measure ? FFTW_MEASURE : FFTW_ESTIMATE );
  INITSTATUS( status, "LALCreateReverseRealFFTPlan", REALFFTC );
  FFTWHOOKS;

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc( sizeof( **plan ) );
  if ( ! *plan )
  {
    ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
  }

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;

  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = rfftw_create_plan( size, FFTW_COMPLEX_TO_REAL, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
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
  rfftw_destroy_plan( (*plan)->plan );
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
  fftw_real *tmp;
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

  rfftw_one( plan->plan, (fftw_real *) input->data, tmp );

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
  fftw_real *tmp;
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

  rfftw_one( plan->plan, tmp, (fftw_real *) output->data );

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
  fftw_real *tmp;
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

  rfftw_one( plan->plan, (fftw_real *) data->data, tmp );

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

  rfftw_one( plan->plan, (fftw_real *)input->data, (fftw_real *)output->data );

  RETURN( status );
}


#ifdef KEEP_OLD_REAL_FFT

void
LALREAL4VectorSequenceFFT(
    LALStatus           *stat,
    REAL4VectorSequence *vout,
    REAL4VectorSequence *vinp,
    RealFFTPlan         *plan
    )
{
  INITSTATUS( stat, "LALREAL4VectorSequenceFFT", REALFFTC );

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT( vout, stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( vinp, stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL );

  /* make sure that the data exists */
  ASSERT( vout->data, stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( vinp->data, stat, REALFFTH_ENULL, REALFFTH_MSGENULL );

  /* make sure that it is not the same data! */
  ASSERT( vout->data != vinp->data, stat, REALFFTH_ESAME, REALFFTH_MSGESAME );

  /* make sure that the lengths agree */
  ASSERT( plan->size > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( vout->vectorLength == plan->size, stat,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( vinp->vectorLength == plan->size, stat,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );

  ASSERT( vout->length > 0, stat, REALFFTH_ESLEN, REALFFTH_MSGESLEN );
  ASSERT( vinp->length == vout->length, stat,
      REALFFTH_ESLEN, REALFFTH_MSGESLEN );
  rfftw(
      (rfftw_plan) plan->plan, vout->length,
      (fftw_real *)vinp->data, 1, plan->size,
      (fftw_real *)vout->data, 1, plan->size
      );

  /* normal exit */
  RETURN( stat );
}


void
LALEstimateFwdRealFFTPlan (
    LALStatus    *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{
  INITSTATUS (stat, "LALEstimateFwdRealFFTPlan", REALFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the plan has not been previously defined */
  ASSERT (*plan == NULL, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

  /* make sure that the requested size is valid */
  ASSERT (size > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc (sizeof(RealFFTPlan));
  ASSERT (*plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = (void *)
    rfftw_create_plan (
        size,
        FFTW_REAL_TO_COMPLEX,
        FFTW_THREADSAFE | FFTW_ESTIMATE /* | FFTW_USE_WISDOM */
        );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

void
LALEstimateInvRealFFTPlan (
    LALStatus    *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{
  INITSTATUS (stat, "LALEstimateInvRealFFTPlan", REALFFTC);
 
  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the plan has not been previously defined */
  ASSERT (*plan == NULL, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

  /* make sure that the requested size is valid */
  ASSERT (size > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc (sizeof(RealFFTPlan));
  ASSERT (*plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = (void *)
    rfftw_create_plan (
        size,
        FFTW_COMPLEX_TO_REAL,
        FFTW_THREADSAFE | FFTW_ESTIMATE /* | FFTW_USE_WISDOM */
        );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

void
LALMeasureFwdRealFFTPlan (
    LALStatus    *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{
  INITSTATUS (stat, "LALMeasureFwdRealFFTPlan", REALFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the plan has not been previously defined */
  ASSERT (*plan == NULL, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

  /* make sure that the requested size is valid */
  ASSERT (size > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc (sizeof(RealFFTPlan));
  ASSERT (*plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = 1;
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = (void *)
    rfftw_create_plan (
        size,
        FFTW_REAL_TO_COMPLEX,
        FFTW_THREADSAFE | FFTW_MEASURE /* | FFTW_USE_WISDOM */
        );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

void
LALMeasureInvRealFFTPlan (
    LALStatus    *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{
  INITSTATUS (stat, "LALMeasureInvRealFFTPlan", REALFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the plan has not been previously defined */
  ASSERT (*plan == NULL, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

  /* make sure that the requested size is valid */
  ASSERT (size > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);

  /* allocate memory */
  *plan = (RealFFTPlan *) LALMalloc (sizeof(RealFFTPlan));
  ASSERT (*plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* assign plan fields */
  (*plan)->size = size;
  (*plan)->sign = -1;
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  (*plan)->plan = (void *)
    rfftw_create_plan (
        size,
        FFTW_COMPLEX_TO_REAL,
        FFTW_THREADSAFE | FFTW_MEASURE /* | FFTW_USE_WISDOM */
        );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* check that the plan is not NULL */
  if ( !(*plan)->plan )
  {
    LALFree( *plan );
    *plan = NULL;
    ABORT( stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}



void
LALFwdRealFFT (
    LALStatus      *stat,
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
  ASSERT (vout, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (vinp, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

  /* make sure that sizes agree */
  n = plan->size;
  ASSERT (n > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
  ASSERT (vinp->length == n, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vout->length == n/2 + 1, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == 1, stat, REALFFTH_ESIGN, REALFFTH_MSGESIGN);

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
    LALStatus      *stat,
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
  ASSERT (vout, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (vinp, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

  /* make sure that sizes agree */
  n = plan->size;
  ASSERT (n > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
  ASSERT (vinp->length == n/2 + 1, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vout->length == n, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == -1, stat, REALFFTH_ESIGN, REALFFTH_MSGESIGN);

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
    ASSERT (vinp->data[n/2].im == 0, stat, REALFFTH_EDATA, REALFFTH_MSGEDATA);
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
LALFwdRealSequenceFFT (
    LALStatus              *stat,
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
  ASSERT (vout, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (vinp, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

  /* make sure that sizes agree */
  m = sqin.length = vinp->length;
  n = sqin.vectorLength = plan->size;
  ASSERT (m > 0, stat, REALFFTH_ESLEN, REALFFTH_MSGESLEN);
  ASSERT (n > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
  ASSERT (vinp->vectorLength == n, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vout->vectorLength == n/2 + 1, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vout->length == m, stat, REALFFTH_ESLEN, REALFFTH_MSGESLEN);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == 1, stat, REALFFTH_ESIGN, REALFFTH_MSGESIGN);

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
    LALStatus              *stat,
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
  ASSERT (vout, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (vinp, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

  /* make sure that sizes agree */
  m = sqin.length = vinp->length;
  n = sqin.vectorLength = plan->size;
  ASSERT (m > 0, stat, REALFFTH_ESLEN, REALFFTH_MSGESLEN);
  ASSERT (n > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
  ASSERT (vinp->vectorLength == n/2 + 1, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vout->vectorLength == n, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vout->length == m, stat, REALFFTH_ESLEN, REALFFTH_MSGESLEN);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == -1, stat, REALFFTH_ESIGN, REALFFTH_MSGESIGN);

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
    ASSERT (z[0].im == 0, stat, REALFFTH_EDATA, REALFFTH_MSGEDATA);

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
      ASSERT (z[n/2].im == 0, stat, REALFFTH_EDATA, REALFFTH_MSGEDATA);
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
    LALStatus           *stat,
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
  ASSERT (vout, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (vinp, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the data pointers are not NULL */
  ASSERT (vout->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (vinp->data, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);
  ASSERT (plan->plan, stat, REALFFTH_ENNUL, REALFFTH_MSGENNUL);

  /* make sure that sizes agree */
  m = sqin.length = vinp->length;
  n = sqin.vectorLength = plan->size;
  ASSERT (m > 0, stat, REALFFTH_ESLEN, REALFFTH_MSGESLEN);
  ASSERT (n > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
  ASSERT (vinp->vectorLength == n, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vout->vectorLength == n/2 + 1, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vout->length == m, stat, REALFFTH_ESLEN, REALFFTH_MSGESLEN);

  /* make sure that the correct plan is being used */
  ASSERT (plan->sign == 1, stat, REALFFTH_ESIGN, REALFFTH_MSGESIGN);

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
#endif /* KEEP_OLD_REAL_FFT */

#undef FFTWHOOKS
