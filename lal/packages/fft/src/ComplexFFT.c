#if 0  /* autodoc block */

<lalVerbatim file="ComplexFFTCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{ComplexFFT.c}}
\label{ss:ComplexFFT.c}

Functions for performing complex FFTs.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ComplexFFTCP}
\index{\texttt{LALEstimateFwdComplexFFTPlan()}}
\index{\texttt{LALEstimateInvComplexFFTPlan()}}
\index{\texttt{LALMeasureFwdComplexFFTPlan()}}
\index{\texttt{LALMeasureInvComplexFFTPlan()}}
\index{\texttt{LALDestroyComplexFFTPlan()}}
\index{\texttt{LALCOMPLEX8VectorFFT()}}

\subsubsection*{Description}

This package provides a LAL-style interface with the FFTW fast Fourier
transform package~\cite{fj:1998}.

The routines \texttt{LALEstimateFwdComplexFFTPlan()},
\texttt{LALEstimateInvComplexFFTPlan()}, \texttt{LALMeasureFwdComplexFFTPlan()}, and
\texttt{LALMeasureInvComplexFFTPlan()} create plans for computing the forward and
inverse FFTs.  The optimum plan is either estimated (reasonably fast) or
measured (can be time-consuming, but gives better performance).  The routine
\texttt{LALDestroyComplexFFTPlan()} destroys any of these flavours of plans.

The routine \texttt{LALCOMPLEX8VectorFFT()} performs either the forward or
inverse FFT depending on the plan.  The discrete Fourier transform $H_k$,
$k=0\ldots n-1$ of a vector $h_j$, $j=0\ldots n-1$, of length $n$ is defined
by
\[
  H_k = \sum_{j=0}^{n-1} h_j e^{2\pi ijk/n}
\]
and, similarly, the inverse Fourier transform is defined by
\[
  h_j = \frac{1}{n}\sum_{k=0}^{n-1} H_k e^{-2\pi ijk/n}.
\]
However, the present implementation of the inverse FFT omits the factor of
$1/n$.


\subsubsection*{Operating Instructions}

\begin{verbatim}
const UINT4 n = 17;   /* example length of sequence of vectors */

static LALStatus status; 

ComplexFFTPlan *pfwd = NULL;
ComplexFFTPlan *pinv = NULL;
COMPLEX8Vector *avec = NULL;
COMPLEX8Vector *bvec = NULL;
COMPLEX8Vector *cvec = NULL;

/* create data vector and sequence */
LALEstimateFwdComplexFFTPlan( &status, &pfwd, n );
LALEstimateInvComplexFFTPlan( &status, &pinv, n );
LALCCreateVector( &status, &avec, n );
LALCCreateVector( &status, &bvec, n );
LALCCreateVector( &status, &cvec, n );

/* assign data ... */

/* perform FFTs */
LALCOMPLEX8VectorFFT( &status, bvec, avec, pfwd );
LALCOMPLEX8VectorFFT( &status, cvec, bvec, pinv );

/* destroy plans, vectors, and sequences */
LALDestroyComplexFFTPlan( &status, &pfwd );
LALDestroyComplexFFTPlan( &status, &pinv );
LALCDestroyVector( &status, &avec );
LALCDestroyVector( &status, &bvec );
LALCDestroyVector( &status, &cvec );
\end{verbatim}

\subsubsection*{Algorithm}

The FFTW~\cite{fj:1998} is used.

\subsubsection*{Uses}

\subsubsection*{Notes}

\begin{enumerate}
\item The sign convention used here agrees with the definition in
\textit{Numerical Recipes}~\cite{ptvf:1992}, but is opposite from the one used
by FFTW~\cite{fj:1998}.  It is also the opposite of that used by the
\texttt{datacondAPI}.  To convert, use the relation:
\verb+H_LAL[k] = H_datacond[n-k]+ (for $0<\verb+k+\le\verb+n/2+$).
\item The result of the inverse FFT must be multiplied by $1/n$ to recover the
original vector.  This is different from the convension used in the
\texttt{datacondAPI} where the factor is applied by default.
\item The size $n$ of the transform can be any positive integer; the
performance is $O(n\log n)$.  However, better performance is obtained if $n$
is the product of powers of 2, 3, 5, 7, and zero or one power of either 11 or
13.  Transforms when $n$ is a power of 2 are especially fast.  See
Ref.~\cite{fj:1998}.
\item LALMalloc() is used by all the fftw routines.
\end{enumerate}

\vfill{\footnotesize\input{ComplexFFTCV}}

</lalLaTeX>

#endif /* autodoc block */


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
  do { fftw_malloc_hook = LALMalloc; fftw_free_hook = LALFree; } while(0)


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALEstimateFwdComplexFFTPlan (
    LALStatus       *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{ /* </lalVerbatim> */
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


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALEstimateInvComplexFFTPlan (
    LALStatus       *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{ /* </lalVerbatim> */
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


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALMeasureFwdComplexFFTPlan (
    LALStatus       *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{ /* </lalVerbatim> */
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


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALMeasureInvComplexFFTPlan (
    LALStatus       *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    )
{ /* </lalVerbatim> */
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


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALDestroyComplexFFTPlan (
    LALStatus       *stat,
    ComplexFFTPlan **plan
    )
{ /* </lalVerbatim> */
  INITSTATUS (stat, "LALDestroyComplexFFTPlan", COMPLEXFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT( plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* check that the plan is not NULL */
  ASSERT( *plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* destroy plan and set to NULL pointer */
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  fftw_destroy_plan( (fftw_plan)(*plan)->plan );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
  LALFree( *plan );
  *plan = NULL;

  /* normal exit */
  RETURN( stat );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCOMPLEX8VectorFFT (
    LALStatus      *stat,
    COMPLEX8Vector *vout,
    COMPLEX8Vector *vinp,
    ComplexFFTPlan *plan
    )
{ /* </lalVerbatim> */
  INITSTATUS( stat, "LALCOMPLEX8VectorFFT", COMPLEXFFTC );

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT( vout, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( vinp, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( plan, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* make sure that the data exists */
  ASSERT( vout->data, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( vinp->data, stat, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );

  /* make sure that it is not the same data! */
  ASSERT( vout->data != vinp->data, stat,
          COMPLEXFFTH_ESAME, COMPLEXFFTH_MSGESAME );

  /* make sure that the lengths agree */
  ASSERT( plan->size > 0, stat, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
  ASSERT( vout->length == plan->size, stat,
          COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );
  ASSERT( vinp->length == plan->size, stat,
          COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );

  fftw_one(
      (fftw_plan) plan->plan,
      (fftw_complex *) vinp->data,
      (fftw_complex *) vout->data
      );

  /* normal exit */
  RETURN( stat );
}

