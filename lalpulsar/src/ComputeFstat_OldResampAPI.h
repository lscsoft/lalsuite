/*
 * Copyright (C) 2009 Chris Messenger, Pinkesh Patel
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
 * \author Chris Messenger
 * \date 2005
 * \ingroup pulsarCoherent
 * \file
 * \brief Header-file defining the API for the F-statistic functions.
 *
 * This code is (partly) a descendant of an earlier implementation found in
 * LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
 * ComputSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
 * LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
 *
 */

#ifndef _COMPUTEFSTAT_OLDRESAMPAPI_H  /* Double-include protection. */
#define _COMPUTEFSTAT_OLDRESAMPAPI_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- exported INCLUDES ----------*/
#include <lal/LALComputeAM.h>
#include <lal/ComputeFstat.h>
#include <lal/PulsarDataTypes.h>
#include <lal/DetectorStates.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/*---------- exported DEFINES ----------*/

/*----- Error-codes -----*/
#define COMPUTEFSTATRSC_ENULL 		1
#define COMPUTEFSTATRSC_ENONULL 	2
#define COMPUTEFSTATRSC_EINPUT   	3
#define COMPUTEFSTATRSC_EMEM   		4
#define COMPUTEFSTATRSC_EXLAL		5
#define COMPUTEFSTATRSC_EIEEE		6

#define COMPUTEFSTATRSC_MSGENULL 	"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATRSC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPUTEFSTATRSC_MSGEINPUT   	"Invalid input"
#define COMPUTEFSTATRSC_MSGEMEM   	"Out of memory. Bad."
#define COMPUTEFSTATRSC_MSGEXLAL	"XLAL function call failed"
#define COMPUTEFSTATRSC_MSGEIEEE	"Floating point failure"

/*---------- exported types ----------*/

/** Multi-IFO container for COMPLEX8 resampled timeseries */
typedef struct tagMultiCOMPLEX8TimeSeries {
  UINT4 length;                         /**< number of IFOs */
  COMPLEX8TimeSeries **data;        /**< array of COMPLEX8TimeSeries (pointers) */
} MultiCOMPLEX8TimeSeries;

/**
 * Struct holding buffered ComputeFStat()-internal quantities to avoid unnecessarily
 * recomputing things that depend ONLY on the skyposition and detector-state series (but not on the spins).
 * For the first call of ComputeFStatFreqBand_RS() the pointer-entries should all be NULL.
 */
struct tagComputeFBuffer_RS {
  MultiDetectorStateSeries *multiDetStates;             /**< buffer for each detStates (store pointer) and skypos */
  REAL8 Alpha, Delta;				              /**< skyposition of candidate */
  LIGOTimeGPS segstart;                                       /**< the start time of the first SFT of the first detector (used to check if the segment has changed) */
  MultiSSBtimes *multiSSB;
  MultiSSBtimes *multiBinary;
  MultiAMCoeffs *multiAMcoef;
  MultiCmplxAMCoeffs *multiCmplxAMcoef;
  MultiCOMPLEX8TimeSeries *multiTimeseries;                   /**< the buffered unweighted multi-detector timeseries */
  MultiCOMPLEX8TimeSeries *multiFa_resampled;                 /**< the buffered multi-detector resampled timeseries weighted by a(t) */
  MultiCOMPLEX8TimeSeries *multiFb_resampled;                 /**< the buffered multi-detector resampled timeseries weighted by b(t) */
};

/**
 * Struct holding a vector of buffered ComputeFStat()-internal quantities to avoid unnecessarily
 * recomputing things that depend ONLY on the skyposition and detector-state series (but not on the spins).
 */
typedef struct tagComputeFBufferVector_RS {
  ComputeFBuffer_RS **data;                                    /**< pointer to a series of ComputeFBuffer_RS structures */
  UINT4 length;                                               /**< the length of the vector */
} ComputeFBufferVector_RS;

/*---------- exported prototypes [API] ----------*/

void ComputeFStatFreqBand_RS ( LALStatus *status,
			       REAL4FrequencySeries *fstatVector,
			       const PulsarDopplerParams *doppler,
			       MultiSFTVector *multiSFTs,                   /* modified */
			       const MultiNoiseWeights *multiWeights,
			       ComputeFParams *params
			       );

void ResampleMultiSFTs ( LALStatus *status,
			 MultiCOMPLEX8TimeSeries **multitimeseries,
			 REAL8 deltaF,
			 const MultiAMCoeffs *multiAMcoef,
			 const MultiSSBtimes *multiSSB,
			 const MultiSFTVector *multiSFTs
			 );

MultiCOMPLEX8TimeSeries *XLALMultiSFTVectorToCOMPLEX8TimeSeries (
                         MultiSFTVector *multisfts  /**< [in] multi SFT vector, gets modified! */
			 );


int XLALEarliestMultiSFTsample ( LIGOTimeGPS *out,
				 MultiSFTVector *multisfts
				 );

int XLALLatestMultiSFTsample ( LIGOTimeGPS *out,
			       MultiSFTVector *multisfts
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

int XLALEarliestMultiSSBtime ( LIGOTimeGPS *out,              /**< output GPS time */
			       const MultiSSBtimes *multiSSB,        /**< input multi SSB SFT-midpoint timestamps */
			       const REAL8 Tsft                     /**< the length of an SFT */
			       );

int XLALLatestMultiSSBtime ( LIGOTimeGPS *out,                   /**< output latest GPS time */
			     const MultiSSBtimes *multiSSB,      /**< input multi SSB SFT-midpoint timestamps */
			     const REAL8 Tsft                    /**< the length of an SFT */
			     );

int XLALGSLInterpolateREAL8Vector ( REAL8Vector **yi,            /**< output interpolated timeseries */
				    REAL8Vector *xi,              /**< input interpolation points */
				    gsl_spline *spline		/**< [in] pre-computed spline data */
				    );

int XLALGSLInitInterpolateREAL8Vector( gsl_spline **spline,
				       REAL8Vector *x,
				       REAL8Vector *y
				       );

int XLALFFTShiftCOMPLEX8Vector(COMPLEX8Vector **x);

int XLALFrequencyShiftMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries **x,	/**< [in/out] timeseries to time-shift */
						const REAL8 shift );	        /**< [in] freq-shift in Hz */

int XLALFrequencyShiftCOMPLEX8TimeSeries ( COMPLEX8TimeSeries **x,	/**< [in/out] timeseries to time-shift */
					   const REAL8 shift );	        /**< [in] freq-shift in Hz */

int XLALSpinDownCorrectionMultiFaFb ( MultiCOMPLEX8TimeSeries **Fa,	/**< [in/out] timeseries to time-shift */
				      MultiCOMPLEX8TimeSeries **Fb,	/**< [in/out] timeseries to time-shift */
				      const PulsarDopplerParams *doppler		/**< parameter-space point to correct for */
				      );

void XLALDestroyMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes );
COMPLEX8TimeSeries *XLALDuplicateCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *times );
MultiCOMPLEX8TimeSeries *XLALDuplicateMultiCOMPLEX8TimeSeries ( MultiCOMPLEX8TimeSeries *multiTimes );

void XLALEmptyComputeFBuffer_RS ( ComputeFBuffer_RS *cfb);

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
