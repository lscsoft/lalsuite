#if 0  /* autodoc block */

<lalVerbatim file="RealFFTCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{RealFFT.c}}
\label{ss:RealFFT.c}

Functions for performing real FFTs.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{RealFFTCP}
\index{\texttt{LALEstimateFwdRealFFTPlan()}}
\index{\texttt{LALEstimateInvRealFFTPlan()}}
\index{\texttt{LALMeasureFwdRealFFTPlan()}}
\index{\texttt{LALMeasureInvRealFFTPlan()}}
\index{\texttt{LALDestroyRealFFTPlan()}}
\index{\texttt{LALREAL4VectorFFT()}}
\index{\texttt{LALREAL4VectorSequenceFFT()}}
\index{\texttt{LALFwdRealFFT()}}
\index{\texttt{LALInvRealFFT()}}
\index{\texttt{LALRealPowerSpectrum()}}
\index{\texttt{LALFwdRealSequenceFFT()}}
\index{\texttt{LALInvRealSequenceFFT()}}
\index{\texttt{LALRealSequencePowerSpectrum()}}

\subsubsection*{Description}

This package provides a LAL-style interface with the FFTW fast Fourier
transform package~\cite{fj:1998}.

The routines \texttt{LALEstimateFwdRealFFTPlan()},
\texttt{LALEstimateInvRealFFTPlan()}, \texttt{LALMeasureFwdRealFFTPlan()}, and
\texttt{LALMeasureInvRealFFTPlan()} create plans for computing the forward
(real-to-complex) and inverse (complex-to-real) FFTs.  The optimum plan is
either estimated (reasonably fast) or measured (can be time-consuming, but
gives better performance).  The routine \texttt{LALDestroyRealFFTPlan()}
destroys any of these flavours of plans.

The routines \texttt{LALREAL4VectorFFT()} and
\texttt{LALREAL4VectorSequenceFFT()} are essentially direct calls to FFTW
routines without any re-packing of the data.  These routines should not be used
unless the user understands the packing used in FFTW.

The routines \texttt{LALFwdRealFFT()} and \texttt{LALInvRealFFT()} perform the
forward (real-to-complex) and inverse (complex-to-real) FFTs using the plans.
The discrete Fourier transform $H_k$, $k=0\ldots \lfloor n/2\rfloor$ ($n/2$
rounded down), of a vector $h_j$, $j=0\ldots n-1$, of length $n$ is defined by
\[
  H_k = \sum_{j=0}^{n-1} h_j e^{2\pi ijk/n}
\]
and, similarly, the inverse Fourier transform is defined by
\[
  h_j = \frac{1}{n} \sum_{k=0}^{n-1} H_k e^{-2\pi ijk/n}
\]
where $H_k$ for $\lfloor n/2\rfloor<k<n$ can be obtained from the relation
$H_k=H_{n-k}^\ast$.  The present implementation of the inverse FFT omits the
factor of $1/n$.

The routines in this package require that the vector $h_j$, $j=0\ldots n-1$ be
real; consequently, $H_k=H_{n-k}^\ast$ ($0\le k\le\lfloor n/2\rfloor$), i.e.,
the negative frequency Fourier components are the complex conjugate of the
positive frequency Fourier components when the data is real.  Therefore, one
need compute and store only the first $\lfloor n/2\rfloor+1$ components of
$H_k$; only the values of $H_k$ for $k=0\ldots \lfloor n/2\rfloor$ are
returned (integer division is rounded down, e.g., $\lfloor 7/2\rfloor=3$).

The routine \texttt{LALRealPowerSpectrum()} computes the power spectrum
$P_k=|H_k|^2$, $k=0\ldots \lfloor n/2\rfloor$, of the data $h_j$,
$j=0\ldots n-1$.

The routines \texttt{LALFwdRealSequenceFFT()} and
\texttt{LALInvRealSequenceFFT()} operate on sequences of vectors
$\{h^{(0)}_{j_0},\ldots,h^{(m-1)}_{j_{m-1}}\}$, $j_0,\ldots,j_{m-1}=0\ldots
n-1$, to produce the Fourier components
$\{H^{(0)}_{k_0},\ldots,H^{(m-1)}_{k_{m-1}}\}$, $k_0,\ldots,k_{m-1}=0\ldots
\lfloor n/2\rfloor$, and vice versa respectively.  Here, $m$ is the length of the sequence
of vectors and $n$ and $\lfloor n/2\rfloor+1$ are the lengths of the vectors themselves.
Similarly, the routine \texttt{LALRealSequencePowerSpectrum()} computes the power
spectra $\{P^{(0)}_{k_0}=|H^(0)_{k_0}|^2,\ldots,P^{(m-1)}_{k_{m-1}}=
|H^{(m-1)}_{k_{m-1}}|^2\}$, $k_0,\ldots,k_{m-1}=0\ldots n/2$, of the sequence
of vectors $\{h^{(0)}_{j_0},\ldots,h^{(m-1)}_{j_{m-1}}\}$,
$j_0,\ldots,j_{m-1}=0\ldots n-1$.

\subsubsection*{Operating Instructions}

\begin{verbatim}
const INT4 m = 4;   /* example length of sequence of vectors */
const INT4 n = 32;  /* example vector length */

static LALStatus status; 

RealFFTPlan            *pfwd = NULL;
RealFFTPlan            *pinv = NULL;
REAL4Vector            *hvec = NULL;
COMPLEX8Vector         *Hvec = NULL;
REAL4Vector            *Pvec = NULL;
REAL4VectorSequence    *hseq = NULL;
COMPLEX8VectorSequence *Hseq = NULL;
REAL4VectorSequence    *Pseq = NULL;
CreateVectorSequenceIn  seqinp;
CreateVectorSequenceIn  seqout;

/* create data vector and sequence */
seqinp.length       = m;
seqinp.vectorLength = n;
seqout.length       = m;
seqout.vectorLength = n/2 + 1;

LALEstimateFwdRealFFTPlan( &status, &pfwd, n );
LALEstimateInvRealFFTPlan( &status, &pinv, n );
LALSCreateVector( &status, &hvec, n );
LALCCreateVector( &status, &Hvec, n/2 + 1 );
LALSCreateVector( &status, &Pvec, n/2 + 1 );
LALSCreateVectorSequence( &status, &hseq, &seqinp );
LALCCreateVectorSequence( &status, &Hseq, &seqout );
LALSCreateVectorSequence( &status, &Pseq, &seqout );

/* assign data ... */

/* perform FFTs */
LALFwdRealFFT( &status, Hvec, hvec, pfwd );
LALInvRealFFT( &status, hvec, Hvec, pinv );
LALRealPowerSpectrum( &status, Pvec, hvec, pfwd );
LALFwdRealSequenceFFT( &status, Hseq, hseq, pfwd );
LALInvRealSequenceFFT( &status, hseq, Hseq, pinv );
LALRealSequencePowerSpectrum( &status, Pseq, hseq, pfwd );

/* destroy plans, vectors, and sequences */
LALDestroyRealFFTPlan( &status, &pfwd );
LALDestroyRealFFTPlan( &status, &pinv );
LALSDestroyVector( &status, &hvec );
LALCDestroyVector( &status, &Hvec );
LALSDestroyVector( &status, &Pvec );
LALSDestroyVectorSequence( &status, &hseq );
LALCDestroyVectorSequence( &status, &Hseq );
LALSDestroyVectorSequence( &status, &Pseq );
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
\verb+H_LAL[k].im = -H_datacondAPI[k].im+.
\item The result of the inverse FFT must be multiplied by $1/n$ to recover the
original vector.  This is unlike the \textit{Numerical
Recipes}~\cite{ptvf:1992} convension where the factor is $2/n$ for real FFTs.
This is also different from the \texttt{datacondAPI} where the normalization
constant is applied by default.
\item The size $n$ of the transform can be any positive integer; the
performance is $O(n\log n)$.  However, better performance is obtained if $n$
is the product of powers of 2, 3, 5, 7, and zero or one power of either 11 or
13.  Transforms when $n$ is a power of 2 are especially fast.  See
Ref.~\cite{fj:1998}.
\item All of these routines leave the input array undamaged.  (Except for
\verb+LALREAL4VectorFFT+ and \verb+LALREAL4VectorSequenceFFT+.)
\item LALMalloc() is used by all the fftw routines.
\end{enumerate}

\vfill{\footnotesize\input{RealFFTCV}}

</lalLaTeX>

#endif /* autodoc block */


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

#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#else
#define pthread_mutex_lock( pmut )
#define pthread_mutex_unlock( pmut )
#endif

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

/* tell FFTW to use LALMalloc and LALFree */
#define FFTWHOOKS \
  do { fftw_malloc_hook = LALMalloc; fftw_free_hook = LALFree; } while(0)

/* <lalVerbatim file="RealFFTCP"> */  
void
LALEstimateFwdRealFFTPlan (
    LALStatus    *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{ /* </lalVerbatim> */
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
    ABORT( stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

/* <lalVerbatim file="RealFFTCP"> */  
void
LALEstimateInvRealFFTPlan (
    LALStatus    *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{ /* </lalVerbatim> */
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
    ABORT( stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

/* <lalVerbatim file="RealFFTCP"> */  
void
LALMeasureFwdRealFFTPlan (
    LALStatus    *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{ /* </lalVerbatim> */
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
    ABORT( stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

/* <lalVerbatim file="RealFFTCP"> */  
void
LALMeasureInvRealFFTPlan (
    LALStatus    *stat,
    RealFFTPlan **plan,
    UINT4         size
    )
{ /* </lalVerbatim> */
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
    ABORT( stat, REALFFTH_ENULL, REALFFTH_MSGENULL );
  }

  /* normal exit */
  RETURN (stat);
}

/* <lalVerbatim file="RealFFTCP"> */  
void
LALDestroyRealFFTPlan (
    LALStatus    *stat,
    RealFFTPlan **plan
    )
{ /* </lalVerbatim> */
  INITSTATUS (stat, "LALDestroyRealFFTPlan", REALFFTC);

  FFTWHOOKS;

  /* make sure that the argument is not NULL */
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* check that the plan is not NULL */
  ASSERT (*plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* destroy plan and set to NULL pointer */
  pthread_mutex_lock( &mut );  
  rfftw_destroy_plan ((rfftw_plan)(*plan)->plan);
  pthread_mutex_unlock( &mut );  
  LALFree (*plan);
  *plan = NULL;

  /* normal exit */
  RETURN (stat);
}


/* <lalVerbatim file="RealFFTCP"> */  
void
LALREAL4VectorFFT (
    LALStatus   *stat,
    REAL4Vector *vout,
    REAL4Vector *vinp,
    RealFFTPlan *plan
    )
{ /* </lalVerbatim> */
  INITSTATUS (stat, "LALREAL4VectorFFT", REALFFTC);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (vinp, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the data exists */
  ASSERT (vout->data, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (vinp->data, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that it is not the same data! */
  ASSERT (vout->data != vinp->data, stat, REALFFTH_ESAME, REALFFTH_MSGESAME);

  /* make sure that the lengths agree */
  ASSERT (plan->size > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
  ASSERT (vout->length == plan->size, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vinp->length == plan->size, stat, REALFFTH_ESZMM, REALFFTH_MSGESZMM);

  rfftw_one (
      (rfftw_plan) plan->plan,
      (fftw_real *)vinp->data,
      (fftw_real *)vout->data
      );

  /* normal exit */
  RETURN (stat);
}


/* <lalVerbatim file="RealFFTCP"> */  
void
LALREAL4VectorSequenceFFT (
    LALStatus           *stat,
    REAL4VectorSequence *vout,
    REAL4VectorSequence *vinp,
    RealFFTPlan         *plan
    )
{ /* </lalVerbatim> */
  INITSTATUS (stat, "LALREAL4VectorSequenceFFT", REALFFTC);

  FFTWHOOKS;

  /* make sure that the arguments are not NULL */
  ASSERT (vout, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (vinp, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (plan, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that the data exists */
  ASSERT (vout->data, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);
  ASSERT (vinp->data, stat, REALFFTH_ENULL, REALFFTH_MSGENULL);

  /* make sure that it is not the same data! */
  ASSERT (vout->data != vinp->data, stat, REALFFTH_ESAME, REALFFTH_MSGESAME);

  /* make sure that the lengths agree */
  ASSERT (plan->size > 0, stat, REALFFTH_ESIZE, REALFFTH_MSGESIZE);
  ASSERT (vout->vectorLength == plan->size, stat,
          REALFFTH_ESZMM, REALFFTH_MSGESZMM);
  ASSERT (vinp->vectorLength == plan->size, stat,
          REALFFTH_ESZMM, REALFFTH_MSGESZMM);

  ASSERT (vout->length > 0, stat, REALFFTH_ESLEN, REALFFTH_MSGESLEN);
  ASSERT (vinp->length == vout->length, stat, REALFFTH_ESLEN, REALFFTH_MSGESLEN);
  rfftw (
      (rfftw_plan) plan->plan, vout->length,
      (fftw_real *)vinp->data, 1, plan->size,
      (fftw_real *)vout->data, 1, plan->size
      );

  /* normal exit */
  RETURN (stat);
}


/* <lalVerbatim file="RealFFTCP"> */  
void
LALFwdRealFFT (
    LALStatus      *stat,
    COMPLEX8Vector *vout,
    REAL4Vector    *vinp,
    RealFFTPlan    *plan
    )
{ /* </lalVerbatim> */
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


/* <lalVerbatim file="RealFFTCP"> */  
void
LALInvRealFFT (
    LALStatus      *stat,
    REAL4Vector    *vout,
    COMPLEX8Vector *vinp,
    RealFFTPlan    *plan
    )
{ /* </lalVerbatim> */
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


/* <lalVerbatim file="RealFFTCP"> */  
void
LALRealPowerSpectrum (
    LALStatus   *stat,
    REAL4Vector *vout,
    REAL4Vector *vinp,
    RealFFTPlan *plan
    )
{ /* </lalVerbatim> */
  REAL4Vector *vtmp = NULL;
  UINT4        n;
  UINT4        k;

  INITSTATUS (stat, "LALRealPowerSpectrum", REALFFTC);
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


/* <lalVerbatim file="RealFFTCP"> */  
void
LALFwdRealSequenceFFT (
    LALStatus              *stat,
    COMPLEX8VectorSequence *vout,
    REAL4VectorSequence    *vinp,
    RealFFTPlan            *plan
    )
{ /* </lalVerbatim> */
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


/* <lalVerbatim file="RealFFTCP"> */  
void
LALInvRealSequenceFFT (
    LALStatus              *stat,
    REAL4VectorSequence    *vout,
    COMPLEX8VectorSequence *vinp,
    RealFFTPlan            *plan
    )
{ /* </lalVerbatim> */
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


/* <lalVerbatim file="RealFFTCP"> */  
void
LALRealSequencePowerSpectrum (
    LALStatus           *stat,
    REAL4VectorSequence *vout,
    REAL4VectorSequence *vinp,
    RealFFTPlan         *plan
    )
{ /* </lalVerbatim> */
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

#undef FFTWHOOKS
