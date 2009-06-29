/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
 * \idx{LALCreateForwardREAL4FFTPlan()}
 * \idx{LALCreateReverseREAL4FFTPlan()}
 * \idx{LALDestroyREAL4FFTPlan()}
 * \idx{LALForwardREAL4FFT()}
 * \idx{LALReverseREAL4FFT()}
 * \idx{LALREAL4PowerSpectrum()}
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
 * $P_k=2|H_k|^2$, $k=1\ldots \lfloor (n-1)/2\rfloor$,
 * $P_0=|H_0|^2$, and $P_{n/2}=|H_{n/2}|^2$ if $n$ is even, of the data $h_j$,
 * $j=0\ldots n-1$.  The factor of two except at DC and Nyquist accounts for
 * the power in negative frequencies.
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

#include <fftw3.h>

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/FFTWMutex.h>


NRCSID( REALFFTC, "$Id$" );


struct
tagREAL4FFTPlan
{
  INT4       sign;
  UINT4      size;
  fftwf_plan plan;
};

struct
tagREAL8FFTPlan
{
  INT4       sign;
  UINT4      size;
  fftw_plan  plan;
};


/*
 *
 * REAL4 XLAL Functions
 *
 */


REAL4FFTPlan * XLALCreateREAL4FFTPlan( UINT4 size, int fwdflg, int measurelvl )
{
  static const char *func = "XLALCreateREAL4FFTPlan";
  REAL4FFTPlan *plan;
  REAL4 *tmp1;
  REAL4 *tmp2;
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
  plan = XLALMalloc( sizeof( *plan ) );
  tmp1 = XLALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = XLALMalloc( size * sizeof( *tmp2 ) );
  if ( ! plan || ! tmp1 || ! tmp2 )
  {
    XLALFree( plan );
    XLALFree( tmp1 );
    XLALFree( tmp2 );
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  if ( fwdflg ) /* forward */
    plan->plan = fftwf_plan_r2r_1d( size, tmp1, tmp2, FFTW_R2HC, flags );
  else /* reverse */
    plan->plan = fftwf_plan_r2r_1d( size, tmp1, tmp2, FFTW_HC2R, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* free temporary arrays */
  XLALFree( tmp2 );
  XLALFree( tmp1 );

  /* check to see success of plan creation */
  if ( ! plan->plan )
  {
    XLALFree( plan );
    XLAL_ERROR_NULL( func, XLAL_EFAILED );
  }

  /* now set remaining plan fields */
  plan->size = size;
  plan->sign = ( fwdflg ? -1 : 1 );

  return plan;
}


REAL4FFTPlan * XLALCreateForwardREAL4FFTPlan( UINT4 size, int measurelvl )
{
  static const char *func = "XLALCreateForwardREAL4FFTPlan";
  REAL4FFTPlan *plan;
  plan = XLALCreateREAL4FFTPlan( size, 1, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return plan;
}


REAL4FFTPlan * XLALCreateReverseREAL4FFTPlan( UINT4 size, int measurelvl )
{
  static const char *func = "XLALCreateReverseREAL4FFTPlan";
  REAL4FFTPlan *plan;
  plan = XLALCreateREAL4FFTPlan( size, 0, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return plan;
}


void XLALDestroyREAL4FFTPlan( REAL4FFTPlan *plan )
{
  if ( plan )
  {
    if ( plan->plan )
    {
      LAL_FFTW_PTHREAD_MUTEX_LOCK;
      fftwf_destroy_plan( plan->plan );
      LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
    }
    memset( plan, 0, sizeof( *plan ) );
    XLALFree( plan );
  }
  return;
}

int XLALREAL4ForwardFFT( COMPLEX8Vector *output, const REAL4Vector *input,
    const REAL4FFTPlan *plan )
{
  static const char *func = "XLALREAL4ForwardFFT";
  REAL4 *tmp;
  UINT4 k;

  if ( ! output || ! input || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size || plan->sign != -1 )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! output->data || ! input->data )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( input->length != plan->size || output->length != plan->size/2 + 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* create temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( func, XLAL_ENOMEM );

  /* do the fft */
  fftwf_execute_r2r( plan->plan, input->data, tmp );

  /* now unpack the results into the output vector */

  /* dc component */
  output->data[0].re = tmp[0];
  output->data[0].im = 0.0;

  /* other components */
  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size/2 rounded up */
  {
    output->data[k].re = tmp[k];
    output->data[k].im = tmp[plan->size - k];
  }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* n is even */
  {
    output->data[plan->size/2].re = tmp[plan->size/2];
    output->data[plan->size/2].im = 0.0;
  }

  XLALFree( tmp );
  return 0;
}


int XLALREAL4ReverseFFT( REAL4Vector *output, const COMPLEX8Vector *input,
    const REAL4FFTPlan *plan )
{
  static const char *func = "XLALREAL4ReverseFFT";
  REAL4 *tmp;
  UINT4 k;

  if ( ! output || ! input || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size || plan->sign != 1 )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! output->data || ! input->data )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( output->length != plan->size || input->length != plan->size/2 + 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );
  if ( input->data[0].im != 0.0 )
    XLAL_ERROR( func, XLAL_EDOM );  /* imaginary part of DC must be zero */
  if ( ! plan->size % 2 && input->data[plan->size/2].im != 0.0 )
    XLAL_ERROR( func, XLAL_EDOM );  /* imaginary part of Nyquist must be zero */

  /* create temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( func, XLAL_ENOMEM );

  /* unpack input into temporary array */

  /* dc component */
  tmp[0] = input->data[0].re;

  /* other components */
  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size / 2 rounded up */
  {
    tmp[k]              = input->data[k].re;
    tmp[plan->size - k] = input->data[k].im;
  }

  /* Nyquist component */
  if ( plan->size%2 == 0 ) /* n is even */
    tmp[plan->size/2] = input->data[plan->size/2].re;

  /* perform the fft */
  fftwf_execute_r2r( plan->plan, tmp, output->data );

  /* cleanup and exit */
  XLALFree( tmp );
  return 0;
}


int XLALREAL4VectorFFT( REAL4Vector *output, const REAL4Vector *input,
    const REAL4FFTPlan *plan )
{
  static const char *func = "XLALREAL4VectorFFT";
  if ( ! output || ! input || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! output->data || ! input->data || output->data == input->data )
    XLAL_ERROR( func, XLAL_EINVAL ); /* note: must be out-of-place */
  if ( output->length != plan->size || input->length != plan->size )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* do the fft */
  fftwf_execute_r2r( plan->plan, input->data, output->data );
  return 0;
}


int XLALREAL4PowerSpectrum( REAL4Vector *spec, const REAL4Vector *data,
    const REAL4FFTPlan *plan )
{
  static const char *func = "XLALREAL4PowerSpectrum";
  REAL4 *tmp;
  UINT4 k;

  if ( ! spec || ! data || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! spec->data || ! data->data )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( data->length != plan->size || spec->length != plan->size/2 + 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* allocate temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( func, XLAL_ENOMEM );

  /* transform the data */
  fftwf_execute_r2r( plan->plan, data->data, tmp );

  /* now reconstruct the spectrum from the temporary storage */

  /* dc component */
  spec->data[0] = tmp[0] * tmp[0];

  /* other components */
  for (k = 1; k < (plan->size + 1)/2; ++k) /* k < size/2 rounded up */
  {
    REAL4 re = tmp[k];
    REAL4 im = tmp[plan->size - k];
    spec->data[k]  = re * re + im * im;
    spec->data[k] *= 2.0; /* accounts for negative frequency part */
  }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* size is even */
    spec->data[k] = tmp[k] * tmp[k];

  /* clenup and exit */
  XLALFree( tmp );
  return 0;
}



/*
 *
 * REAL8 XLAL Functions
 *
 */



REAL8FFTPlan * XLALCreateREAL8FFTPlan( UINT4 size, int fwdflg, int measurelvl )
{
  static const char *func = "XLALCreateREAL8FFTPlan";
  REAL8FFTPlan *plan;
  REAL8 *tmp1;
  REAL8 *tmp2;
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
  plan = XLALMalloc( sizeof( *plan ) );
  tmp1 = XLALMalloc( size * sizeof( *tmp1 ) );
  tmp2 = XLALMalloc( size * sizeof( *tmp2 ) );
  if ( ! plan || ! tmp1 || ! tmp2 )
  {
    XLALFree( plan );
    XLALFree( tmp1 );
    XLALFree( tmp2 );
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  }

  LAL_FFTW_PTHREAD_MUTEX_LOCK;
  if ( fwdflg ) /* forward */
    plan->plan = fftw_plan_r2r_1d( size, tmp1, tmp2, FFTW_R2HC, flags );
  else /* reverse */
    plan->plan = fftw_plan_r2r_1d( size, tmp1, tmp2, FFTW_HC2R, flags );
  LAL_FFTW_PTHREAD_MUTEX_UNLOCK;

  /* free temporary arrays */
  XLALFree( tmp2 );
  XLALFree( tmp1 );

  /* check to see success of plan creation */
  if ( ! plan->plan )
  {
    XLALFree( plan );
    XLAL_ERROR_NULL( func, XLAL_EFAILED );
  }

  /* now set remaining plan fields */
  plan->size = size;
  plan->sign = ( fwdflg ? -1 : 1 );

  return plan;
}


REAL8FFTPlan * XLALCreateForwardREAL8FFTPlan( UINT4 size, int measurelvl )
{
  static const char *func = "XLALCreateForwardREAL8FFTPlan";
  REAL8FFTPlan *plan;
  plan = XLALCreateREAL8FFTPlan( size, 1, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return plan;
}


REAL8FFTPlan * XLALCreateReverseREAL8FFTPlan( UINT4 size, int measurelvl )
{
  static const char *func = "XLALCreateReverseREAL8FFTPlan";
  REAL8FFTPlan *plan;
  plan = XLALCreateREAL8FFTPlan( size, 0, measurelvl );
  if ( ! plan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  return plan;
}


void XLALDestroyREAL8FFTPlan( REAL8FFTPlan *plan )
{
  if ( plan )
  {
    if ( plan->plan )
    {
      LAL_FFTW_PTHREAD_MUTEX_LOCK;
      fftw_destroy_plan( plan->plan );
      LAL_FFTW_PTHREAD_MUTEX_UNLOCK;
    }
    memset( plan, 0, sizeof( *plan ) );
    XLALFree( plan );
  }
  return;
}

int XLALREAL8ForwardFFT( COMPLEX16Vector *output, REAL8Vector *input,
    const REAL8FFTPlan *plan )
{
  static const char *func = "XLALREAL8ForwardFFT";
  REAL8 *tmp;
  UINT4 k;

  if ( ! output || ! input || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size || plan->sign != -1 )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! output->data || ! input->data )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( input->length != plan->size || output->length != plan->size/2 + 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* create temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( func, XLAL_ENOMEM );

  /* do the fft */
  fftw_execute_r2r( plan->plan, input->data, tmp );

  /* now unpack the results into the output vector */

  /* dc component */
  output->data[0].re = tmp[0];
  output->data[0].im = 0.0;

  /* other components */
  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size/2 rounded up */
  {
    output->data[k].re = tmp[k];
    output->data[k].im = tmp[plan->size - k];
  }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* n is even */
  {
    output->data[plan->size/2].re = tmp[plan->size/2];
    output->data[plan->size/2].im = 0.0;
  }

  XLALFree( tmp );
  return 0;
}


int XLALREAL8ReverseFFT( REAL8Vector *output, COMPLEX16Vector *input,
    const REAL8FFTPlan *plan )
{
  static const char *func = "XLALREAL8ReverseFFT";
  REAL8 *tmp;
  UINT4 k;

  if ( ! output || ! input || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size || plan->sign != 1 )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! output->data || ! input->data )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( output->length != plan->size || input->length != plan->size/2 + 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );
  if ( input->data[0].im != 0.0 )
    XLAL_ERROR( func, XLAL_EDOM );  /* imaginary part of DC must be zero */
  if ( ! plan->size % 2 && input->data[plan->size/2].im != 0.0 )
    XLAL_ERROR( func, XLAL_EDOM );  /* imaginary part of Nyquist must be zero */

  /* create temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( func, XLAL_ENOMEM );

  /* unpack input into temporary array */

  /* dc component */
  tmp[0] = input->data[0].re;

  /* other components */
  for ( k = 1; k < (plan->size + 1)/2; ++k ) /* k < size / 2 rounded up */
  {
    tmp[k]              = input->data[k].re;
    tmp[plan->size - k] = input->data[k].im;
  }

  /* Nyquist component */
  if ( plan->size%2 == 0 ) /* n is even */
    tmp[plan->size/2] = input->data[plan->size/2].re;

  /* perform the fft */
  fftw_execute_r2r( plan->plan, tmp, output->data );

  /* cleanup and exit */
  XLALFree( tmp );
  return 0;
}


int XLALREAL8VectorFFT( REAL8Vector *output, REAL8Vector *input,
    const REAL8FFTPlan *plan )
{
  static const char *func="XLALREAL8VectorFFT";
  if ( ! output || ! input || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! output->data || ! input->data || output->data == input->data )
    XLAL_ERROR( func, XLAL_EINVAL ); /* note: must be out-of-place */
  if ( output->length != plan->size || input->length != plan->size )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* do the fft */
  fftw_execute_r2r( plan->plan, input->data, output->data );
  return 0;
}


int XLALREAL8PowerSpectrum( REAL8Vector *spec, REAL8Vector *data,
    const REAL8FFTPlan *plan )
{
  static const char *func = "XLALREAL8PowerSpectrum";
  REAL8 *tmp;
  UINT4 k;

  if ( ! spec || ! data || ! plan )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! plan->plan || ! plan->size )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! spec->data || ! data->data )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( data->length != plan->size || spec->length != plan->size/2 + 1 )
    XLAL_ERROR( func, XLAL_EBADLEN );

  /* allocate temporary storage space */
  tmp = XLALMalloc( plan->size * sizeof( *tmp ) );
  if ( ! tmp )
    XLAL_ERROR( func, XLAL_ENOMEM );

  /* transform the data */
  fftw_execute_r2r( plan->plan, data->data, tmp );

  /* now reconstruct the spectrum from the temporary storage */

  /* dc component */
  spec->data[0] = tmp[0] * tmp[0];

  /* other components */
  for (k = 1; k < (plan->size + 1)/2; ++k) /* k < size/2 rounded up */
  {
    REAL8 re = tmp[k];
    REAL8 im = tmp[plan->size - k];
    spec->data[k]  = re * re + im * im;
    spec->data[k] *= 2.0; /* accounts for negative frequency part */
  }

  /* Nyquist frequency */
  if ( plan->size%2 == 0 ) /* size is even */
    spec->data[plan->size/2] = tmp[plan->size/2] * tmp[plan->size/2];

  /* clenup and exit */
  XLALFree( tmp );
  return 0;
}




/*
 *
 * LAL Functions
 *
 */




/* <lalVerbatim file="RealFFTCP"> */
void
LALCreateForwardREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCreateForwardREAL4FFTPlan", REALFFTC );
  XLALPrintDeprecationWarning("LALCreateForwardREAL4FFTPlan", "XLALCreateForwardREAL4FFTPlan");

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  *plan = XLALCreateREAL4FFTPlan( size, 1, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALCreateReverseREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCreateReverseREAL4FFTPlan", REALFFTC );
  XLALPrintDeprecationWarning("LALCreateReverseREAL4FFTPlan", "XLALCreateReverseREAL4FFTPlan");

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  *plan = XLALCreateREAL4FFTPlan( size, 0, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALDestroyREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALDestroyREAL4FFTPlan", REALFFTC );
  XLALPrintDeprecationWarning("LALDestroyREAL4FFTPlan", "XLALDestroyREAL4FFTPlan");
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( *plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  XLALDestroyREAL4FFTPlan( *plan );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }
  *plan = NULL;
  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALForwardREAL4FFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    REAL4Vector    *input,
    REAL4FFTPlan    *plan
    )
{ /* </lalVerbatim> */
  int code;
  UINT4 n;
  INITSTATUS( status, "LALForwardREAL4FFT", REALFFTC );
  XLALPrintDeprecationWarning("LALForwardREAL4FFT", "XLALForwardREAL4FFT");

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

  code = XLALREAL4ForwardFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( plan->sign != -1 ) /* plan sign was wrong */
        {
          ABORT( status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALReverseREAL4FFT(
    LALStatus      *status,
    REAL4Vector    *output,
    COMPLEX8Vector *input,
    REAL4FFTPlan    *plan
    )
{ /* </lalVerbatim> */
  int code;
  UINT4 n;
  INITSTATUS( status, "LALReverseREAL4FFT", REALFFTC );
  XLALPrintDeprecationWarning("LALReverseREAL4FFT", "XLALReverseREAL4FFT");

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( output->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->data[0].im == 0, status, REALFFTH_EDATA, REALFFTH_MSGEDATA );
  ASSERT( n % 2 || input->data[n / 2].im == 0, status,
      REALFFTH_EDATA, REALFFTH_MSGEDATA );

  ASSERT( plan->sign == 1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  code = XLALREAL4ReverseFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( plan->sign != 1 ) /* plan sign was wrong */
        {
          ABORT( status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      case XLAL_EDOM: /* either DC or Nyquist imaginary part was non zero */
        ABORT( status, REALFFTH_EDATA, REALFFTH_MSGEDATA );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALREAL4PowerSpectrum (
    LALStatus   *status,
    REAL4Vector *spec,
    REAL4Vector *data,
    REAL4FFTPlan *plan
    )
{ /* </lalVerbatim> */
  int code;
  UINT4 n;

  INITSTATUS( status, "LALREAL4PowerSpectrum", REALFFTC );
  XLALPrintDeprecationWarning("LALREAL4PowerSpectrum", "XLALREAL4PowerSpectrum");

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

  code = XLALREAL4PowerSpectrum( spec, data, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALREAL4VectorFFT(
    LALStatus   *status,
    REAL4Vector *output,
    REAL4Vector *input,
    REAL4FFTPlan *plan
    )
{ /* </lalVerbatim> */
  int code;
  INITSTATUS( status, "LALREAL4VectorFFT", REALFFTC );
  XLALPrintDeprecationWarning("LALREAL4VectorFFT", "XLALREAL4VectorFFT");

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

  code = XLALREAL4VectorFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        if ( ! plan->size ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( output->data == input->data ) /* same data pointers */
        {
          ABORT( status, REALFFTH_ESAME, REALFFTH_MSGESAME );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/*
 *
 * LAL REAL8 Functions
 *
 */


/* <lalVerbatim file="RealFFTCP"> */
void
LALCreateForwardREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCreateForwardREAL8FFTPlan", REALFFTC );
  XLALPrintDeprecationWarning("LALCreateForwardREAL8FFTPlan", "XLALCreateForwardREAL8FFTPlan");

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  *plan = XLALCreateREAL8FFTPlan( size, 1, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALCreateReverseREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan,
    UINT4         size,
    INT4          measure
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCreateReverseREAL8FFTPlan", REALFFTC );
  XLALPrintDeprecationWarning("LALCreateReverseREAL8FFTPlan", "XLALCreateReverseREAL8FFTPlan");

  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( ! *plan, status, REALFFTH_ENNUL, REALFFTH_MSGENNUL );
  ASSERT( size > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );

  *plan = XLALCreateREAL8FFTPlan( size, 0, measure );
  if ( ! *plan )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EBADLEN:
        ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EFAILED:
        ABORT( status, REALFFTH_EFFTW, REALFFTH_MSGEFFTW );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALDestroyREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALDestroyREAL8FFTPlan", REALFFTC );
  XLALPrintDeprecationWarning("LALDestroyREAL8FFTPlan", "XLALDestroyREAL8FFTPlan");
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( *plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  XLALDestroyREAL8FFTPlan( *plan );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
      default:
        ABORTXLAL( status );
    }
  }
  *plan = NULL;
  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALForwardREAL8FFT(
    LALStatus      *status,
    COMPLEX16Vector *output,
    REAL8Vector    *input,
    REAL8FFTPlan    *plan
    )
{ /* </lalVerbatim> */
  int code;
  UINT4 n;
  INITSTATUS( status, "LALForwardREAL8FFT", REALFFTC );
  XLALPrintDeprecationWarning("LALForwardREAL8FFT", "XLALForwardREAL8FFT");

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

  code = XLALREAL8ForwardFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( plan->sign != -1 ) /* plan sign was wrong */
        {
          ABORT( status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALReverseREAL8FFT(
    LALStatus      *status,
    REAL8Vector    *output,
    COMPLEX16Vector *input,
    REAL8FFTPlan    *plan
    )
{ /* </lalVerbatim> */
  int code;
  UINT4 n;
  INITSTATUS( status, "LALReverseREAL8FFT", REALFFTC );
  XLALPrintDeprecationWarning("LALReverseREAL8FFT", "XLALReverseREAL8FFT");

  ASSERT( output, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  ASSERT( output->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( input->data, status, REALFFTH_ENULL, REALFFTH_MSGENULL );
  ASSERT( plan->plan, status, REALFFTH_ENULL, REALFFTH_MSGENULL );

  n = plan->size;
  ASSERT( n > 0, status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
  ASSERT( output->length == n, status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->length == n / 2 + 1, status,
      REALFFTH_ESZMM, REALFFTH_MSGESZMM );
  ASSERT( input->data[0].im == 0, status, REALFFTH_EDATA, REALFFTH_MSGEDATA );
  ASSERT( n % 2 || input->data[n / 2].im == 0, status,
      REALFFTH_EDATA, REALFFTH_MSGEDATA );

  ASSERT( plan->sign == 1, status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );

  code = XLALREAL8ReverseFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( plan->sign != 1 ) /* plan sign was wrong */
        {
          ABORT( status, REALFFTH_ESIGN, REALFFTH_MSGESIGN );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      case XLAL_EDOM: /* either DC or Nyquist imaginary part was non zero */
        ABORT( status, REALFFTH_EDATA, REALFFTH_MSGEDATA );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALREAL8PowerSpectrum (
    LALStatus   *status,
    REAL8Vector *spec,
    REAL8Vector *data,
    REAL8FFTPlan *plan
    )
{ /* </lalVerbatim> */
  int code;
  UINT4 n;

  INITSTATUS( status, "LALREAL8PowerSpectrum", REALFFTC );
  XLALPrintDeprecationWarning("LALREAL8PowerSpectrum", "XLALREAL8PowerSpectrum");

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

  code = XLALREAL8PowerSpectrum( spec, data, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_ENOMEM:
        ABORT( status, REALFFTH_EALOC, REALFFTH_MSGEALOC );
      case XLAL_EINVAL:
        if ( ! n ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="RealFFTCP"> */
void
LALREAL8VectorFFT(
    LALStatus   *status,
    REAL8Vector *output,
    REAL8Vector *input,
    REAL8FFTPlan *plan
    )
{ /* </lalVerbatim> */
  int code;
  INITSTATUS( status, "LALREAL8VectorFFT", REALFFTC );
  XLALPrintDeprecationWarning("LALREAL8VectorFFT", "XLALREAL8VectorFFT");

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

  code = XLALREAL8VectorFFT( output, input, plan );
  if ( code )
  {
    code = xlalErrno;
    XLALClearErrno();
    switch ( code )
    {
      case XLAL_EINVAL:
        if ( ! plan->size ) /* plan size was invalid */
        {
          ABORT( status, REALFFTH_ESIZE, REALFFTH_MSGESIZE );
        }
        else if ( output->data == input->data ) /* same data pointers */
        {
          ABORT( status, REALFFTH_ESAME, REALFFTH_MSGESAME );
        }
        else /* one of the data pointers was NULL */
        {
          ABORT( status, REALFFTH_ENULL, REALFFTH_MSGENULL );
        }
      case XLAL_EBADLEN: /* size mismatch */
        ABORT( status, REALFFTH_ESZMM, REALFFTH_MSGESZMM );
      default:
        ABORTXLAL( status );
    }
  }

  RETURN( status );
}

