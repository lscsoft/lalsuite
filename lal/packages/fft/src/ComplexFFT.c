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

#ifdef LAL_FFTW3_ENABLED
#include <fftw3.h>
#else
#error "FFTW3 is **REQUIRED**"
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FFTWMutex.h>

NRCSID( COMPLEXFFTC, "$Id$" );


struct
tagCOMPLEX8FFTPlan
{
  INT4       sign;
  UINT4      size;
  fftwf_plan plan;
};

struct
tagCOMPLEX16FFTPlan
{
  INT4       sign;
  UINT4      size;
  fftw_plan  plan;
};


/*
 *
 * XLAL COMPLEX8 functions
 *
 */


COMPLEX8FFTPlan * XLALCreateCOMPLEX8FFTPlan( UINT4 size, int fwdflg, int measurelvl )
{
  static const char *func = "XLALCreateCOMPLEX8FFTPlan";
  COMPLEX8FFTPlan *plan;
  COMPLEX8 *tmp1;
  COMPLEX8 *tmp2;
  int flags = FFTW_UNALIGNED;

  if ( ! size )
    XLAL_ERROR_NULL( func, XLAL_EBADLEN );

  /* based on measurement level, set fftw3 flags to perform 
   * requested degree of measurement */
  switch ( measurelvl )
  {
    case 0: /* estimate */
      flags |= FFTW_ESTIMATE;
      break;
    default: /* exhaustive measurement */
      flags |= FFTW_EXHAUSTIVE;
      /* fall-through */
    case 2: /* lengthy measurement */
      flags |= FFTW_PATIENT;
      /* fall-through */
    case 1: /* measure the best plan */
      flags |= FFTW_MEASURE;
      break;
  }

  /* allocate memory for the plan and the temporary arrays */
  plan = LALMalloc( sizeof( *plan ) );
  tmp1 = LALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = LALMalloc( size * sizeof( *tmp2 ) );
  if ( ! plan || ! tmp1 || ! tmp2 )
  {
    if ( plan ) LALFree( plan );
    if ( tmp1 ) LALFree( tmp1 );
    if ( tmp2 ) LALFree( tmp2 );
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  /* create the plan */
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  plan->plan = fftwf_plan_dft_1d( size,
      (fftwf_complex *)tmp1, (fftwf_complex *)tmp2,
      fwdflg ? FFTW_FORWARD : FFTW_BACKWARD, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* free temporary arrays */
  LALFree( tmp2 );
  LALFree( tmp1 );

  /* check to see success of plan creation */
  if ( ! plan->plan )
  {
    LALFree( plan );
    XLAL_ERROR_NULL( func, XLAL_EFAILED );
  }

  /* now set remaining plan fields */
  plan->size = size;
  plan->sign = ( fwdflg ? -1 : 1 );

  return plan;
}


COMPLEX8FFTPlan * XLALCreateForwardCOMPLEX8FFTPlan( UINT4 size, int measurelvl )
{
  static const char *func = "XLALCreateForwardCOMPLEX8FFTPlan";  
  COMPLEX8FFTPlan *plan;
  plan = XLALCreateCOMPLEX8FFTPlan( size, 1, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return plan;
}


COMPLEX8FFTPlan * XLALCreateReverseCOMPLEX8FFTPlan( UINT4 size, int measurelvl )
{
  static const char *func = "XLALCreateReverseCOMPLEX8FFTPlan";  
  COMPLEX8FFTPlan *plan;
  plan = XLALCreateCOMPLEX8FFTPlan( size, 0, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return plan;
}


void XLALDestroyCOMPLEX8FFTPlan( COMPLEX8FFTPlan *plan )
{
  static const char *func = "XLALDestroyCOMPLEX8FFTPlan";
  if ( ! plan )
    XLAL_ERROR_VOID( func, XLAL_EFAULT );
  if ( ! plan->plan )
    XLAL_ERROR_VOID( func, XLAL_EINVAL );
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  fftwf_destroy_plan( plan->plan );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
  memset( plan, 0, sizeof( *plan ) );
  LALFree( plan );
  return;
}


int XLALCOMPLEX8VectorFFT( COMPLEX8Vector *output, COMPLEX8Vector *input,
    COMPLEX8FFTPlan *plan )
{
  static const char *func = "XLALCOMPLEX8VectorFFT";
  if ( ! output || ! input || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! output->data || ! input->data || output->data == input->data )
    XLAL_ERROR( func, XLAL_EINVAL ); /* note: must be out-of-place */
  if ( output->length != plan->size || input->length != plan->size )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* do the fft */
  fftwf_execute_dft(
      plan->plan,
      (fftwf_complex *)input->data,
      (fftwf_complex *)output->data
      );
  return 0;
}



/*
 *
 * XLAL COMPLEX16 functions
 *
 */


COMPLEX16FFTPlan * XLALCreateCOMPLEX16FFTPlan( UINT4 size, int fwdflg, int measurelvl )
{
  static const char *func = "XLALCreateCOMPLEX16FFTPlan";
  COMPLEX16FFTPlan *plan;
  COMPLEX16 *tmp1;
  COMPLEX16 *tmp2;
  int flags = FFTW_UNALIGNED;

  if ( ! size )
    XLAL_ERROR_NULL( func, XLAL_EBADLEN );

  /* based on measurement level, set fftw3 flags to perform 
   * requested degree of measurement */
  switch ( measurelvl )
  {
    case 0: /* estimate */
      flags |= FFTW_ESTIMATE;
      break;
    default: /* exhaustive measurement */
      flags |= FFTW_EXHAUSTIVE;
      /* fall-through */
    case 2: /* lengthy measurement */
      flags |= FFTW_PATIENT;
      /* fall-through */
    case 1: /* measure the best plan */
      flags |= FFTW_MEASURE;
      break;
  }

  /* allocate memory for the plan and the temporary arrays */
  plan = LALMalloc( sizeof( *plan ) );
  tmp1 = LALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = LALMalloc( size * sizeof( *tmp2 ) );
  if ( ! plan || ! tmp1 || ! tmp2 )
  {
    if ( plan ) LALFree( plan );
    if ( tmp1 ) LALFree( tmp1 );
    if ( tmp2 ) LALFree( tmp2 );
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  /* create the plan */
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  plan->plan = fftw_plan_dft_1d( size,
      (fftw_complex *)tmp1, (fftw_complex *)tmp2,
      fwdflg ? FFTW_FORWARD : FFTW_BACKWARD, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* free temporary arrays */
  LALFree( tmp2 );
  LALFree( tmp1 );

  /* check to see success of plan creation */
  if ( ! plan->plan )
  {
    LALFree( plan );
    XLAL_ERROR_NULL( func, XLAL_EFAILED );
  }

  /* now set remaining plan fields */
  plan->size = size;
  plan->sign = ( fwdflg ? -1 : 1 );

  return plan;
}


COMPLEX16FFTPlan * XLALCreateForwardCOMPLEX16FFTPlan( UINT4 size, int measurelvl )
{
  static const char *func = "XLALCreateForwardCOMPLEX16FFTPlan";  
  COMPLEX16FFTPlan *plan;
  plan = XLALCreateCOMPLEX16FFTPlan( size, 1, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return plan;
}


COMPLEX16FFTPlan * XLALCreateReverseCOMPLEX16FFTPlan( UINT4 size, int measurelvl )
{
  static const char *func = "XLALCreateReverseCOMPLEX16FFTPlan";  
  COMPLEX16FFTPlan *plan;
  plan = XLALCreateCOMPLEX16FFTPlan( size, 0, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return plan;
}


void XLALDestroyCOMPLEX16FFTPlan( COMPLEX16FFTPlan *plan )
{
  static const char *func = "XLALDestroyCOMPLEX16FFTPlan";
  if ( ! plan )
    XLAL_ERROR_VOID( func, XLAL_EFAULT );
  if ( ! plan->plan )
    XLAL_ERROR_VOID( func, XLAL_EINVAL );
  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  fftw_destroy_plan( plan->plan );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
  memset( plan, 0, sizeof( *plan ) );
  LALFree( plan );
  return;
}


int XLALCOMPLEX16VectorFFT( COMPLEX16Vector *output, COMPLEX16Vector *input,
    COMPLEX16FFTPlan *plan )
{
  static const char *func = "XLALCOMPLEX16VectorFFT";
  if ( ! output || ! input || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! output->data || ! input->data || output->data == input->data )
    XLAL_ERROR( func, XLAL_EINVAL ); /* note: must be out-of-place */
  if ( output->length != plan->size || input->length != plan->size )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* do the fft */
  fftw_execute_dft(
      plan->plan,
      (fftw_complex *)input->data,
      (fftw_complex *)output->data
      );
  return 0;
}



/*
 *
 * LAL COMPLEX8 functions
 *
 */


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCreateForwardCOMPLEX8FFTPlan(
    LALStatus       *status,
    COMPLEX8FFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCreateForwardCOMPLEX8FFTPlan", COMPLEXFFTC );

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  *plan = XLALCreateCOMPLEX8FFTPlan( size, 1, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCreateReverseCOMPLEX8FFTPlan(
    LALStatus       *status,
    COMPLEX8FFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCreateReverseCOMPLEX8FFTPlan", COMPLEXFFTC );

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  *plan = XLALCreateCOMPLEX8FFTPlan( size, 0, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALDestroyCOMPLEX8FFTPlan (
    LALStatus       *status,
    COMPLEX8FFTPlan **plan
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALDestroyCOMPLEX8FFTPlan", COMPLEXFFTC );
  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( *plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  XLALDestroyCOMPLEX8FFTPlan( *plan );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT( status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCOMPLEX8VectorFFT (
    LALStatus      *status,
    COMPLEX8Vector *output,
    COMPLEX8Vector *input,
    COMPLEX8FFTPlan *plan
    )
{ /* </lalVerbatim> */
  int code;
  INITSTATUS( status, "LALCOMPLEX8VectorFFT", COMPLEXFFTC );

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

  code = XLALCOMPLEX8VectorFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        if ( ! plan->size ) /* plan size was invalid */
        {
          ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
        }
        else if ( output->data == input->data ) /* same data pointers */
        {
          ABORT( status, COMPLEXFFTH_ESAME, COMPLEXFFTH_MSGESAME );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/*
 *
 * LAL COMPLEX16 functions
 *
 */


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCreateForwardCOMPLEX16FFTPlan(
    LALStatus       *status,
    COMPLEX16FFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCreateForwardCOMPLEX16FFTPlan", COMPLEXFFTC );

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  *plan = XLALCreateCOMPLEX16FFTPlan( size, 1, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCreateReverseCOMPLEX16FFTPlan(
    LALStatus       *status,
    COMPLEX16FFTPlan **plan,
    UINT4            size,
    INT4             measure
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCreateReverseCOMPLEX16FFTPlan", COMPLEXFFTC );

  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( ! *plan, status, COMPLEXFFTH_ENNUL, COMPLEXFFTH_MSGENNUL );
  ASSERT( size > 0, status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );

  *plan = XLALCreateCOMPLEX16FFTPlan( size, 0, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, COMPLEXFFTH_EALOC, COMPLEXFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, COMPLEXFFTH_EFFTW, COMPLEXFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALDestroyCOMPLEX16FFTPlan (
    LALStatus       *status,
    COMPLEX16FFTPlan **plan
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALDestroyCOMPLEX16FFTPlan", COMPLEXFFTC );
  ASSERT( plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  ASSERT( *plan, status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
  XLALDestroyCOMPLEX16FFTPlan( *plan );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT( status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }
  RETURN( status );
}


/* <lalVerbatim file="ComplexFFTCP"> */  
void
LALCOMPLEX16VectorFFT (
    LALStatus      *status,
    COMPLEX16Vector *output,
    COMPLEX16Vector *input,
    COMPLEX16FFTPlan *plan
    )
{ /* </lalVerbatim> */
  int code;
  INITSTATUS( status, "LALCOMPLEX16VectorFFT", COMPLEXFFTC );

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

  code = XLALCOMPLEX16VectorFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        if ( ! plan->size ) /* plan size was invalid */
        {
          ABORT( status, COMPLEXFFTH_ESIZE, COMPLEXFFTH_MSGESIZE );
        }
        else if ( output->data == input->data ) /* same data pointers */
        {
          ABORT( status, COMPLEXFFTH_ESAME, COMPLEXFFTH_MSGESAME );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, COMPLEXFFTH_ENULL, COMPLEXFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, COMPLEXFFTH_ESZMM, COMPLEXFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}

