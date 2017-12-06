/*
*  Copyright (C) 2007 Jolien Creighton, Josh Willis
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

#include <fftw3.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FFTWMutex.h>
#include <lal/LALConfig.h> /* Needed to know whether aligning memory */

/**
 * \addtogroup ComplexFFT_h
 *
 * ### Description ###
 *
 * This package provides a (X)LAL-style interface with the FFTW fast Fourier
 * transform package \cite fj_1998.
 *
 * The routines LALCreateForwardComplexFFTPlan() and
 * LALCreateReverseComplexFFTPlan() create plans for computing the
 * forward and reverse FFTs of a given size.  The optimum plan is either
 * estimated (reasonably fast) if the measure flag is zero, or measured (can be
 * time-consuming, but gives better performance) if the measure flag is
 * non-zero.  The routine LALDestroyComplexFFTPlan() destroys either
 * of these flavours of plans.
 *
 * The routine LALCOMPLEX8VectorFFT() performs either the forward or
 * reverse FFT depending on the plan.  The discrete Fourier transform \f$H_k\f$,
 * \f$k=0\ldots n-1\f$ of a vector \f$h_j\f$, \f$j=0\ldots n-1\f$, of length \f$n\f$ is defined
 * by
 * \f[
 * H_k = \sum_{j=0}^{n-1} h_j e^{-2\pi ijk/n}
 * \f]
 * and, similarly, the \e inverse Fourier transform is defined by
 * \f[
 * h_j = \frac{1}{n}\sum_{k=0}^{n-1} H_k e^{2\pi ijk/n}.
 * \f]
 * However, the present implementation of the \e reverse FFT omits the
 * factor of \f$1/n\f$.  The input and output vectors must be distinct.
 *
 * ### Operating Instructions ###
 *
 * \code
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
 * <li> The sign convention used here is the opposite of the definition in
 * <em>Numerical Recipes</em> \cite ptvf1992, but agrees with the one used
 * by FFTW \cite fj_1998 and the other LIGO software components.
 * </li><li> The result of the inverse FFT must be multiplied by \f$1/n\f$ to recover
 * the original vector.  This is different from the \c datacondAPI where
 * the factor is applied by default.
 * </li><li> The size \f$n\f$ of the transform can be any positive integer; the
 * performance is \f$O(n\log n)\f$.  However, better performance is obtained if \f$n\f$
 * is the product of powers of 2, 3, 5, 7, and zero or one power of either 11
 * or 13.  Transforms when \f$n\f$ is a power of 2 are especially fast.  See
 * Ref. \cite fj_1998.
 * </li><li> LALMalloc() is used by all the fftw routines.
 * </li><li> The input and output vectors for LALCOMPLEX8VectorFFT() must
 * be distinct.
 * </li></ol>
 *
 */
/*@{*/

/**
 * Plan to perform an FFT of COMPLEX8 data
 */
struct
tagCOMPLEX8FFTPlan
{
  INT4       sign; /**< sign in transform exponential, -1 for forward, +1 for reverse */
  UINT4      size; /**< length of the complex data vector for this plan */
  fftwf_plan plan; /**< the FFTW plan */
};

/**
 * Plan to perform an FFT of COMPLEX16 data
 */
struct
tagCOMPLEX16FFTPlan
{
  INT4       sign; /**< sign in transform exponential, -1 for forward, +1 for reverse */
  UINT4      size; /**< length of the complex data vector for this plan */
  fftw_plan  plan; /**< the FFTW plan */
};

/* single- and double-precision routines */

#define SINGLE_PRECISION
#include "ComplexFFT_source.c"
#undef SINGLE_PRECISION
#include "ComplexFFT_source.c"
