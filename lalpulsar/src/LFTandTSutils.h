/*
 *  Copyright (C) 2009 Reinhard Prix, Chris Messenger, Pinkesh Patel
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

/**
 * \defgroup LFTandTSutils_h Header LFTandTSutils.h
 * \ingroup pkg_SFTIO
 * \author Reinhard Prix, Chris Messenger
 * \date 2009
 * \brief Utility functions for working with Long Fourier Transforms and Time Series.
 */
/*@{*/

#ifndef _LFTANDTSUTILS_H  /* Double-include protection. */
#define _LFTANDTSUTILS_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/*---------- exported INCLUDES ----------*/
#include <gsl/gsl_spline.h>

#include <lal/SFTutils.h>
#include <lal/LALDatatypes.h>
#include <lal/LALComputeAM.h>
#include <lal/SSBtimes.h>

/*---------- exported DEFINES ----------*/

/*---------- exported MACROS ----------*/
#define NhalfPosDC(N) ( (UINT4)( N / 2 ) + 1 )
#define NhalfNeg(N) ( (N) - NhalfPosDC(N) )		/* making sure NhalfPosDC(N) + NhalfNeg(N) = N */

/*---------- exported types ----------*/

/** Multi-IFO container for COMPLEX8 resampled timeseries */
typedef struct tagMultiCOMPLEX8TimeSeries
{
  UINT4 length;                         /**< number of IFOs */
  COMPLEX8TimeSeries **data;        	/**< array of COMPLEX8TimeSeries (pointers) */
} MultiCOMPLEX8TimeSeries;


/** Struct holding the results of comparing two floating-point vectors (real-valued or complex),
 * using various different comparison metrics
 */
typedef struct tagVectorComparison
{
  REAL4 relErr_L1;		///< relative error between vectors using L1 norm \f$r_1(x,y) \equiv \frac{|x - y|_1}{0.5(|x|_1 + |y|_1)}\f$
  REAL4 relErr_L2;		///< relative error between vectors using L2 norm \f$r_2(x,y) \equiv \frac{|x - y|_2}{0.5(|x|_2 + |y|_2)}\f$
  REAL4 angleV;			///< angle between the two vectors, according to \f$\cos\theta = \frac{\Re(x\cdot y^*)}{|x|_2 \, |y|_2}\f$
  REAL4 relErr_atMaxAbsx;	///< single-sample relative error *at* maximum |sample-value| of first vector 'x'
  REAL4 relErr_atMaxAbsy;	///< single-sample relative error *at* maximum |sample-value| of second vector 'x'
} VectorComparison;

/*---------- exported Global variables ----------*/

/*---------- exported prototypes [API] ----------*/

COMPLEX8TimeSeries *XLALSFTVectorToCOMPLEX8TimeSeries ( const SFTVector *sfts, const LIGOTimeGPS *start_in, const LIGOTimeGPS *end_in );
MultiCOMPLEX8TimeSeries *XLALMultiSFTVectorToCOMPLEX8TimeSeries ( const MultiSFTVector *multisfts );
SFTtype *XLALSFTVectorToLFT ( SFTVector *sfts, REAL8 upsampling );

int XLALReorderFFTWtoSFT (COMPLEX8Vector *X);
int XLALReorderSFTtoFFTW (COMPLEX8Vector *X);
int XLALTimeShiftSFT ( SFTtype *sft, REAL8 shift );

int XLALGSLInterpolateREAL8Vector ( REAL8Vector **yi, REAL8Vector *xi, gsl_spline *spline );

int XLALGSLInitInterpolateREAL8Vector ( gsl_spline **spline, REAL8Vector *x, REAL8Vector *y );

int XLALFFTShiftCOMPLEX8Vector ( COMPLEX8Vector **x );

int XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **x, const REAL8 shift );

int XLALFrequencyShiftCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **x, const REAL8 shift );

int XLALSpinDownCorrectionMultiFaFb ( MultiCOMPLEX8TimeSeries **Fa, MultiCOMPLEX8TimeSeries **Fb, const PulsarDopplerParams *doppler );

void XLALDestroyMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes );
COMPLEX8TimeSeries *XLALDuplicateCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *times );
MultiCOMPLEX8TimeSeries *XLALDuplicateMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes );

int XLALSincInterpolateCOMPLEX8TimeSeries ( COMPLEX8Vector **y_out, const REAL8Vector *t_out, const COMPLEX8TimeSeries *ts_in, UINT4 Dterms );
int XLALDirichletInterpolateCOMPLEX8FrequencySeries ( COMPLEX8Vector *y_out, const REAL8Vector *f_out, const COMPLEX8FrequencySeries *fs_in, UINT4 Dterms );
SFTtype *XLALDirichletInterpolateSFT ( const SFTtype *sft_in, REAL8 f0Out, REAL8 dfOut, UINT4 numBinsOut, UINT4 Dterms );

int XLALCompareCOMPLEX8Vectors ( VectorComparison *result, const COMPLEX8Vector *x, const COMPLEX8Vector *y, const VectorComparison *tol );
int XLALCompareREAL4Vectors    ( VectorComparison *result, const REAL4Vector *x, const REAL4Vector *y, const VectorComparison *tol );
int XLALCheckVectorComparisonTolerances ( const VectorComparison *result, const VectorComparison *tol );

/*@}*/

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */
