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

/*---------- exported types ----------*/

/** Multi-IFO container for COMPLEX8 resampled timeseries */
typedef struct tagMultiCOMPLEX8TimeSeries {
  UINT4 length;                         /**< number of IFOs */
  COMPLEX8TimeSeries **data;        /**< array of COMPLEX8TimeSeries (pointers) */
} MultiCOMPLEX8TimeSeries;

/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */

/*---------- exported prototypes [API] ----------*/

COMPLEX8TimeSeries *XLALSFTVectorToCOMPLEX8TimeSeries ( SFTVector *sfts,      /**< input SFT vector, gets modified! */
							const LIGOTimeGPS *start_in,    /**< input start time */
							const LIGOTimeGPS *end_in       /**< input end time */
							);

SFTtype *XLALSFTVectorToLFT ( const SFTVector *sfts, REAL8 upsampling );


int XLALReorderFFTWtoSFT (COMPLEX8Vector *X);
int XLALReorderSFTtoFFTW (COMPLEX8Vector *X);
int XLALTimeShiftSFT ( SFTtype *sft, REAL8 shift );
int XLALMultiplySFTbyCOMPLEX8 ( SFTtype *sft, COMPLEX8 factor );

MultiCOMPLEX8TimeSeries *XLALMultiSFTVectorToCOMPLEX8TimeSeries (
                         MultiSFTVector *multisfts  /**< [in] multi SFT vector, gets modified! */
                         );


int XLALAntennaWeightCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **Faoft,                         /**< [out] the timeseries weighted by a(t) */
                                          COMPLEX8TimeSeries **Fboft,                         /**< [out] the timeseries weighted by b(t) */
                                          const COMPLEX8TimeSeries *timeseries,         /**< [in] the input timeseries */
                                          const AMCoeffs *AMcoef,                       /**< [in] the AM coefficients */
                                          const SFTVector *sfts                         /**< [in] the SFT data */
                                          );

int XLALAntennaWeightMultiCOMPLEX8TimeSeries(MultiCOMPLEX8TimeSeries **Faoft,                        /**< [out] the timeseries weighted by a(t) */
                                             MultiCOMPLEX8TimeSeries **Fboft,                         /**< [out] the timeseries weighted by b(t) */
                                             const MultiCOMPLEX8TimeSeries *multiTimeseries,         /**< [in] the input multi-detector timeseries */
                                             const MultiAMCoeffs *multiAMcoef,                        /**< [in] the multi-detector AM coefficients */
                                             const MultiSFTVector *multisfts                        /**< [in] the multi-detector SFT data */
                                             );

int XLALBarycentricResampleMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **Faoft_RS,                         /**< [out] the resampled timeseries Fa(t_SSB) */
                                                     MultiCOMPLEX8TimeSeries **Fboft_RS,                         /**< [out] the resampled timeseries Fb(t_SSB) */
                                                     const MultiCOMPLEX8TimeSeries *Faoft,                       /**< [in] the detector frame timeseries Fa(t) */
                                                     const MultiCOMPLEX8TimeSeries *Fboft,                       /**< [in] the detector frame timeseries Fb(t) */
                                                     const MultiSSBtimes *multiSSB,                              /**< [in] the multi-detector SSB times data */
                                                     const MultiSFTVector *multiSFTs,                            /**< [in] the multi-detector SFT data */
                                                     const REAL8 deltaF                                          /**< [in] the user defined output frequency resolution */
                                                     );

int XLALBarycentricResampleCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **Faoft_RS,                         /**< [out] the resampled timeseries Fa(t_SSB) */
                                                COMPLEX8TimeSeries **Fboft_RS,                         /**< [out] the resampled timeseries Fb(t_SSB) */
                                                const COMPLEX8TimeSeries *Faoft,                       /**< [in] the input detector frame timeseries Fa(t) */
                                                const COMPLEX8TimeSeries *Fboft,                       /**< [in] the input detector frame timeseries Fb(t) */
                                                const SSBtimes *SSB,                                   /**< [in] the SSB times at the midpoints of each SFT */
                                                const SFTVector *SFTs                                  /**< [in] the SFT data */
                                                );

int XLALGSLInterpolateREAL8Vector ( REAL8Vector **yi,            /* output interpolated timeseries */
                                    REAL8Vector *xi,              /* input interpolation points */
                                    gsl_spline *spline         /* [in] pre-computed spline data */
                                    );

int XLALGSLInitInterpolateREAL8Vector( gsl_spline **spline,
                                       REAL8Vector *x,
                                       REAL8Vector *y
                                       );

int XLALFFTShiftCOMPLEX8Vector(COMPLEX8Vector **x);

int XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **x,	/**< [in/out] timeseries to time-shift */
                                                const REAL8 shift );            /**< [in] freq-shift in Hz */

int XLALFrequencyShiftCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **x,	/**< [in/out] timeseries to time-shift */
                                           const REAL8 shift );         /**< [in] freq-shift in Hz */

int XLALSpinDownCorrectionMultiFaFb ( MultiCOMPLEX8TimeSeries **Fa,	/**< [in/out] timeseries to time-shift */
                                      MultiCOMPLEX8TimeSeries **Fb,	/**< [in/out] timeseries to time-shift */
                                      const PulsarDopplerParams *doppler		/**< parameter-space point to correct for */
                                      );

void XLALDestroyMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes );
COMPLEX8TimeSeries *XLALDuplicateCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *times );
MultiCOMPLEX8TimeSeries *XLALDuplicateMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes );

/*@}*/

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */
