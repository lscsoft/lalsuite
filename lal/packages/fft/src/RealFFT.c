/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon, Josh Willis
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

#include <config.h>

#include <complex.h>
#include <fftw3.h>

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/FFTWMutex.h>
#include <lal/LALConfig.h> /* Needed to know whether aligning memory */

/**
 * \addtogroup RealFFT_h
 *
 * \section sec_RealFFT_LAL LAL-style functions [DEPRECATED]
 *
 * This package also provides a (deprecated!) LAL-style interface with the FFTW fast Fourier
 * transform package \cite fj_1998.
 *
 * The routines LALCreateForwardRealFFTPlan() and
 * LALCreateReverseRealFFTPlan() create plans for computing the
 * forward (real-to-complex) and reverse (complex-to-real) FFTs of a specified
 * size.  The optimum plan is either estimated (reasonably fast) if the measure
 * flag is zero, or measured (can be time-consuming, but gives better
 * performance) if the measure flag is non-zero.  The routine
 * LALDestroyRealFFTPlan() destroys any of these flavours of plans.
 *
 * The routines LALForwardRealFFT() and LALReverseRealFFT()
 * perform the forward (real-to-complex) and reverse (complex-to-real) FFTs
 * using the plans.  The discrete Fourier transform \f$H_k\f$,
 * \f$k=0\ldots\lfloor{n/2}\rfloor\f$ (\f$n/2\f$ rounded down), of a vector \f$h_j\f$,
 * \f$j=0\ldots n-1\f$, of length \f$n\f$ is defined by
 * \f[
 * H_k = \sum_{j=0}^{n-1} h_j e^{-2\pi ijk/n}
 * \f]
 * and, similarly, the \e inverse Fourier transform is defined by
 * \f[
 * h_j = \frac{1}{n} \sum_{k=0}^{n-1} H_k e^{2\pi ijk/n}
 * \f]
 * where \f$H_k\f$ for \f$\lfloor{n/2}\rfloor<k<n\f$ can be obtained from the relation
 * \f$H_k=H_{n-k}^\ast\f$.  The present implementation of the \e reverse FFT
 * omits the factor of \f$1/n\f$.
 *
 * The routines in this package require that the vector \f$h_j\f$, \f$j=0\ldots n-1\f$
 * be real; consequently, \f$H_k=H_{n-k}^\ast\f$ (\f$0\le k\le\lfloor n/2\rfloor\f$),
 * i.e., the negative frequency Fourier components are the complex conjugate of
 * the positive frequency Fourier components when the data is real.  Therefore,
 * one need compute and store only the first \f$\lfloor n/2\rfloor+1\f$ components
 * of \f$H_k\f$; only the values of \f$H_k\f$ for \f$k=0\ldots \lfloor n/2\rfloor\f$ are
 * returned (integer division is rounded down, e.g., \f$\lfloor 7/2\rfloor=3\f$).
 *
 * The routine LALRealPowerSpectrum() computes the power spectrum
 * \f$P_k=2|H_k|^2\f$, \f$k=1\ldots \lfloor (n-1)/2\rfloor\f$,
 * \f$P_0=|H_0|^2\f$, and \f$P_{n/2}=|H_{n/2}|^2\f$ if \f$n\f$ is even, of the data \f$h_j\f$,
 * \f$j=0\ldots n-1\f$.  The factor of two except at DC and Nyquist accounts for
 * the power in negative frequencies.
 *
 * The routine LALREAL4VectorFFT() is essentially a direct calls to
 * FFTW routines without any re-packing of the data.  This routine should not
 * be used unless the user understands the packing used in FFTW.
 *
 * \subsection ss_RealFFT_OP Operating Instructions
 *
 * \code
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
 * \endcode
 *
 * ### Algorithm ###
 *
 * The FFTW \cite fj_1998 is used.
 *
 * ### Uses ###
 *
 *
 * ### Notes ###
 *
 * <ol>
 * <li> The sign convention used here is the opposite of
 * <em>Numerical Recipes</em> \cite ptvf1992, but agrees with the one used
 * by FFTW \cite fj_1998 and the other LIGO software components.
 * </li>
 * <li> The result of the reverse FFT must be multiplied by \f$1/n\f$ to recover
 * the original vector.  This is unlike the <em>Numerical
 * Recipes</em> \cite ptvf1992 convension where the factor is \f$2/n\f$ for real
 * FFTs.  This is different from the \c datacondAPI where the
 * normalization constant is applied by default.
 * </li>
 * <li> The size \f$n\f$ of the transform can be any positive integer; the
 * performance is \f$O(n\log n)\f$.  However, better performance is obtained if \f$n\f$
 * is the product of powers of 2, 3, 5, 7, and zero or one power of either 11
 * or 13.  Transforms when \f$n\f$ is a power of 2 are especially fast.  See
 * \cite fj_1998.
 * </li>
 * <li> All of these routines leave the input array undamaged.  (Except for LALREAL4VectorFFT().)
 * </li>
 * <li> LALMalloc() is used by all the fftw routines.
 * </li>
 * </ol>
 *
 */
/*@{*/

/**
 * \brief Plan to perform FFT of REAL4 data.
 */
struct
tagREAL4FFTPlan
{
  INT4       sign; /**< sign in transform exponential, -1 for forward, +1 for reverse */
  UINT4      size; /**< length of the real data vector for this plan */
  fftwf_plan plan; /**< the FFTW plan */
};

/**
 * \brief Plan to perform FFT of REAL8 data.
 */
struct
tagREAL8FFTPlan
{
  INT4       sign; /**< sign in transform exponential, -1 for forward, +1 for reverse */
  UINT4      size; /**< length of the real data vector for this plan */
  fftw_plan  plan; /**< the FFTW plan */
};



/* single- and double-precision routines */

#define SINGLE_PRECISION
#include "RealFFT_source.c"
#undef SINGLE_PRECISION
#include "RealFFT_source.c"
