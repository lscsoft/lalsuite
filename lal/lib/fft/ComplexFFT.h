/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef _COMPLEXFFT_H
#define _COMPLEXFFT_H

#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \defgroup ComplexFFT_h Header ComplexFFT_h
 * \ingroup lal_fft
 *
 * \brief Performs complex-to-complex FFTs.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/ComplexFFT.h>
 * \endcode
 *
 * Perform complex-to-complex fast Fourier transforms of vectors using the
 * package FFTW \cite fj_1998 .
 *
 */
/** @{ */

/** \name Error Codes */
/**@{*/
#define COMPLEXFFTH_ENULL 1	/**< Null pointer */
#define COMPLEXFFTH_ENNUL 2	/**< Non-null pointer */
#define COMPLEXFFTH_ESIZE 4	/**< Invalid input size */
#define COMPLEXFFTH_ESZMM 8	/**< Size mismatch */
#define COMPLEXFFTH_ESLEN 16	/**< Invalid/mismatched sequence lengths */
#define COMPLEXFFTH_ESAME 32	/**< Input/Output data vectors are the same */
#define COMPLEXFFTH_EALOC 64	/**< Memory allocation failed */
#define COMPLEXFFTH_EFFTW 128	/**< Error in FFTW */
#define COMPLEXFFTH_ESNGL 256	/**< FFTW library is not single-precision */
#define COMPLEXFFTH_EINTL 512	/**< Error in Intel FFT library */
#define COMPLEXFFTH_ESIGN 1024	/**< Unknown sign of transform in plan */
/**@}*/

/** \cond DONT_DOXYGEN */
#define COMPLEXFFTH_MSGENULL "Null pointer"
#define COMPLEXFFTH_MSGENNUL "Non-null pointer"
#define COMPLEXFFTH_MSGESIZE "Invalid input size"
#define COMPLEXFFTH_MSGESZMM "Size mismatch"
#define COMPLEXFFTH_MSGESLEN "Invalid/mismatched sequence lengths"
#define COMPLEXFFTH_MSGESAME "Input/Output data vectors are the same"
#define COMPLEXFFTH_MSGEALOC "Memory allocation failed"
#define COMPLEXFFTH_MSGEFFTW "Error in FFTW"
#define COMPLEXFFTH_MSGESNGL "FFTW library is not single-precision"
#define COMPLEXFFTH_MSGEINTL "Error in Intel FFT library"
#define COMPLEXFFTH_MSGESIGN "Unknown sign of transform in plan"
/** \endcond */

/** Plan to perform FFT of COMPLEX8 data */
typedef struct tagCOMPLEX8FFTPlan COMPLEX8FFTPlan;
/** Plan to perform FFT of COMPLEX16 data */
typedef struct tagCOMPLEX16FFTPlan COMPLEX16FFTPlan;
#define tagComplexFFTPlan tagCOMPLEX8FFTPlan
#define ComplexFFTPlan COMPLEX8FFTPlan

#ifdef SWIG /* SWIG interface directives */
SWIGLAL(VIEWIN_ARRAYS(COMPLEX8Vector, output));
SWIGLAL(VIEWIN_ARRAYS(COMPLEX16Vector, output));
#endif /* SWIG */

/*
 *
 * XLAL COMPLEX8 functions
 *
 */

/**
 * Returns a new COMPLEX8FFTPlan
 * A COMPLEX8FFTPlan is required to perform a FFT that involves complex data.
 * A different plan is required for each size of the complex data vectors
 * and for each direction of transform (forward or reverse).
 * A forward transform performs
 * \f[Z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,z[j]\f]
 * where N, the size of the transform, is the length of the vectors z and Z.
 * A reverse transform performs
 * \f[w[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,W[k]\f]
 * where N, the size of the transform, is the length of the vectors w and W.
 *
 * @note
 * The reverse transform of the forward transform of some data is
 * equal to N times the original data (we therefore call it a "reverse"
 * transform rather than an "inverse" transform).
 *
 * @param[in] size The number of points in the complex data.
 * @param[in] fwdflg Set non-zero for a forward FFT plan;
 * otherwise create a reverse plan
 * @param[in] measurelvl Measurement level for plan creation:
 * - 0: no measurement, just estimate the plan;
 * - 1: measure the best plan;
 * - 2: perform a lengthy measurement of the best plan;
 * - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c COMPLEX8FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateCOMPLEX8Plan() function shall fail if:
 * - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 * - [\c XLAL_ENOMEM] Insufficient storage space is available.
 * - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 * .
 */
COMPLEX8FFTPlan * XLALCreateCOMPLEX8FFTPlan( UINT4 size, int fwdflg, int measurelvl );

/**
 * Returns a new COMPLEX8FFTPlan for a forward transform
 *
 * A COMPLEX8FFTPlan is required to perform a FFT that involves complex data.
 * A different plan is required for each size of the complex data vectors.
 * A forward transform performs
 * \f[Z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,z[j]\f]
 * where N, the size of the transform, is the length of the vector z and Z.
 *
 * @param[in] size The number of points in the complex data.
 * @param[in] measurelvl Measurement level for plan creation:
 * - 0: no measurement, just estimate the plan;
 * - 1: measure the best plan;
 * - 2: perform a lengthy measurement of the best plan;
 * - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c COMPLEX8FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateForwardCOMPLEX8Plan() function shall fail if:
 * - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 * - [\c XLAL_ENOMEM] Insufficient storage space is available.
 * - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 * .
 */
COMPLEX8FFTPlan * XLALCreateForwardCOMPLEX8FFTPlan( UINT4 size, int measurelvl );

/**
 * Returns a new COMPLEX8FFTPlan for a reverse transform
 *
 * A COMPLEX8FFTPlan is required to perform a FFT that involves complex data.
 * A reverse transform performs
 * \f[w[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,W[k]\f]
 * where N, the size of the transform, is the length of the vectors w and W.
 *
 * @note
 * The reverse transform of the forward transform of some data is
 * equal to N times the original data (we therefore call it a "reverse"
 * transform rather than an "inverse" transform).
 *
 * @param[in] size The number of points in the complex data.
 * @param[in] measurelvl Measurement level for plan creation:
 * - 0: no measurement, just estimate the plan;
 * - 1: measure the best plan;
 * - 2: perform a lengthy measurement of the best plan;
 * - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c COMPLEX8FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateReverseCOMPLEX8Plan() function shall fail if:
 * - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 * - [\c XLAL_ENOMEM] Insufficient storage space is available.
 * - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 * .
 */
COMPLEX8FFTPlan * XLALCreateReverseCOMPLEX8FFTPlan( UINT4 size, int measurelvl );

/**
 * Destroys a COMPLEX8FFTPlan
 * @param[in] plan A pointer to the COMPLEX8FFTPlan to be destroyed.
 */
void XLALDestroyCOMPLEX8FFTPlan( COMPLEX8FFTPlan *plan );

/**
 * Perform a COMPLEX8Vector to COMPLEX8Vector FFT
 *
 * This routine computes
 * \f[Z[k] = \sum_{j=0}^{N-1} e^{\mp2\pi ijk/N}\,z[j],\f]
 * and where the minus sign is used if a forward plan is provided as the argument
 * and the plus sign is used if a reverse plan is provided as the argument;
 * here N is the length of the input and output vectors z and Z.
 *
 * @param[out] output The complex output data vector Z of length N
 * @param[in] input The input complex data vector z of length N
 * @param[in] plan The FFT plan to use for the transform
 * @note
 * The input and output vectors must be distinct.
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALCOMPLEX8VectorFFT() function shall fail if:
 * - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 * - [\c XLAL_EINVAL] A argument is invalid or the input and output data
 * vectors are the same.
 * - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 * incompatible.
 * - [\c XLAL_ENOMEM] Insufficient storage space is available.
 * .
 */
int XLALCOMPLEX8VectorFFT( COMPLEX8Vector * _LAL_RESTRICT_ output, const COMPLEX8Vector * _LAL_RESTRICT_ input, const COMPLEX8FFTPlan *plan );

/*
 *
 * XLAL COMPLEX16 functions
 *
 */

/**
 * Returns a new COMPLEX16FFTPlan
 *
 * A COMPLEX16FFTPlan is required to perform a FFT that involves complex data.
 * A different plan is required for each size of the complex data vectors
 * and for each direction of transform (forward or reverse).
 * A forward transform performs
 * \f[Z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,z[j]\f]
 * where N, the size of the transform, is the length of the vectors z and Z.
 * A reverse transform performs
 * \f[w[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,W[k]\f]
 * where N, the size of the transform, is the length of the vectors w and W.
 *
 * @note
 * The reverse transform of the forward transform of some data is
 * equal to N times the original data (we therefore call it a "reverse"
 * transform rather than an "inverse" transform).
 *
 * @param[in] size The number of points in the complex data.
 * @param[in] fwdflg Set non-zero for a forward FFT plan;
 * otherwise create a reverse plan
 * @param[in] measurelvl Measurement level for plan creation:
 * - 0: no measurement, just estimate the plan;
 * - 1: measure the best plan;
 * - 2: perform a lengthy measurement of the best plan;
 * - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c COMPLEX16FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateCOMPLEX16Plan() function shall fail if:
 * - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 * - [\c XLAL_ENOMEM] Insufficient storage space is available.
 * - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 * .
 */
COMPLEX16FFTPlan * XLALCreateCOMPLEX16FFTPlan( UINT4 size, int fwdflg, int measurelvl );

/**
 * Returns a new COMPLEX16FFTPlan for a forward transform
 *
 * A COMPLEX16FFTPlan is required to perform a FFT that involves complex data.
 * A different plan is required for each size of the complex data vectors.
 * A forward transform performs
 * \f[Z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,z[j]\f]
 * where N, the size of the transform, is the length of the vector z and Z.
 *
 * @param[in] size The number of points in the complex data.
 * @param[in] measurelvl Measurement level for plan creation:
 * - 0: no measurement, just estimate the plan;
 * - 1: measure the best plan;
 * - 2: perform a lengthy measurement of the best plan;
 * - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c COMPLEX16FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateForwardCOMPLEX16Plan() function shall fail if:
 * - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 * - [\c XLAL_ENOMEM] Insufficient storage space is available.
 * - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 * .
 */
COMPLEX16FFTPlan * XLALCreateForwardCOMPLEX16FFTPlan( UINT4 size, int measurelvl );

/**
 * Returns a new COMPLEX16FFTPlan for a reverse transform
 *
 * A COMPLEX16FFTPlan is required to perform a FFT that involves complex data.
 * A reverse transform performs
 * \f[w[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,W[k]\f]
 * where N, the size of the transform, is the length of the vectors w and W.
 *
 * @note
 * The reverse transform of the forward transform of some data is
 * equal to N times the original data (we therefore call it a "reverse"
 * transform rather than an "inverse" transform).
 *
 * @param[in] size The number of points in the complex data.
 * @param[in] measurelvl Measurement level for plan creation:
 * - 0: no measurement, just estimate the plan;
 * - 1: measure the best plan;
 * - 2: perform a lengthy measurement of the best plan;
 * - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c COMPLEX16FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateReverseCOMPLEX16Plan() function shall fail if:
 * - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 * - [\c XLAL_ENOMEM] Insufficient storage space is available.
 * - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 * .
 */
COMPLEX16FFTPlan * XLALCreateReverseCOMPLEX16FFTPlan( UINT4 size, int measurelvl );

/**
 * Destroys a COMPLEX16FFTPlan
 * @param[in] plan A pointer to the COMPLEX16FFTPlan to be destroyed.
 */
void XLALDestroyCOMPLEX16FFTPlan( COMPLEX16FFTPlan *plan );

/**
 * Perform a COMPLEX16Vector to COMPLEX16Vector FFT
 *
 * This routine computes
 * \f[Z[k] = \sum_{j=0}^{N-1} e^{\mp2\pi ijk/N}\,z[j],\f]
 * and where the minus sign is used if a forward plan is provided as the argument
 * and the plus sign is used if a reverse plan is provided as the argument;
 * here N is the length of the input and output vectors z and Z.
 *
 * @param[out] output The complex output data vector Z of length N
 * @param[in] input The input complex data vector z of length N
 * @param[in] plan The FFT plan to use for the transform
 * @note
 * The input and output vectors must be distinct.
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALCOMPLEX16VectorFFT() function shall fail if:
 * - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 * - [\c XLAL_EINVAL] A argument is invalid or the input and output data
 * vectors are the same.
 * - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 * incompatible.
 * - [\c XLAL_ENOMEM] Insufficient storage space is available.
 * .
 */
int XLALCOMPLEX16VectorFFT( COMPLEX16Vector * _LAL_RESTRICT_ output, const COMPLEX16Vector * _LAL_RESTRICT_ input, const COMPLEX16FFTPlan *plan );

/** @} */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _COMPLEXFFT_H */
