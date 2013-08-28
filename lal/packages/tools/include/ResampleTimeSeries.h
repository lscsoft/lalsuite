/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Patrick Brady
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

#include <lal/LALDatatypes.h>
#include <lal/BandPassTimeSeries.h>

#ifndef _RESAMPLETIMESERIES_H
#define _RESAMPLETIMESERIES_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup ResampleTimeSeries_h
 * \author Brown, D. A.
 *
 * \brief Provides routines to resample a time series.
 *
 * At present only integer downsampling of ::REAL4TimeSeries by a power of two is supported.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/ResampleTimeSeries.h>
 * \endcode
 *
 * This header covers routines that resample time series by applying
 * a low pass filter and decimating the resulting time series. Further
 * documentation is given in the individual routines' modules.
 *
 */
/*@{*/

/**\name Error Codes */
/*@{*/
#define RESAMPLETIMESERIESH_ENULL 1	/**< Null pointer */
#define RESAMPLETIMESERIESH_ENNUL 2	/**< Non-null pointer */
#define RESAMPLETIMESERIESH_EZERO 3	/**< Length of input time series is zero */
#define RESAMPLETIMESERIESH_ERATE 4	/**< Sample rate is zero */
#define RESAMPLETIMESERIESH_EUPSM 5	/**< Cannot upsample */
#define RESAMPLETIMESERIESH_EHIGH 6	/**< Input sample rate is greater than 32kHz */
#define RESAMPLETIMESERIESH_ELOG2 7	/**< Only power-of-two resampling is avaliable */
#define RESAMPLETIMESERIESH_EFILT 8	/**< Unknown filter type */
#define RESAMPLETIMESERIESH_EINVD 9	/**< Invalid or non-integer resample factor */
#define RESAMPLETIMESERIESH_ELDAS 10	/**< Input resample factor with LDAS FIR */
/*@}*/

/** \cond DONT_DOXYGEN */
#define RESAMPLETIMESERIESH_MSGENULL "Null pointer"
#define RESAMPLETIMESERIESH_MSGENNUL "Non-null pointer"
#define RESAMPLETIMESERIESH_MSGEZERO "Length of input time series is zero"
#define RESAMPLETIMESERIESH_MSGERATE "Sample rate is zero"
#define RESAMPLETIMESERIESH_MSGEUPSM "Cannot upsample"
#define RESAMPLETIMESERIESH_MSGEHIGH "Input sample rate is greater than 32kHz"
#define RESAMPLETIMESERIESH_MSGELOG2 "Only power-of-two resampling is avaliable"
#define RESAMPLETIMESERIESH_MSGEFILT "Unknown filter type"
#define RESAMPLETIMESERIESH_MSGEINVD "Invalid or non-integer resample factor"
#define RESAMPLETIMESERIESH_MSGELDAS "Input resample factor with LDAS FIR"
/** \endcond */

/**
 * This enum type contains the different low pass filters available to
 * prevent power above the new Nyquist frequency entering the resampled
 * time series due to aliasing.
 */
typedef enum
{
  defaultButterworth,	/**< An IIR butterwoth filter of order 20 with attenuation 0.1 at the new Nyquist frequency.
                         * See the package tdfilters for documentation of butterworth filters in LAL.
                         */
  LDASfirLP		/**< For downsampling by a factor of 2, 4 or 8 an
                         * implementation of the FIR filter used by the LDAS datacondAPI resample().
                         * This is provided for testing the result of standalone codes and codes
                         * running under LDAS. The LDAS filter provided here has filter order parameter
                         * 10, so the order of the filter is \f$2 \times 10 \times q\f$ where \f$q\f$ is the
                         * resampling ratio.
                         */
}
ResampleTSFilter;


/**
 * This union is provided so that the code can store the parameters of the
 * filter in a place accessible by the user for user designed low pass filters.
 * This is not presently implemented and this structure may be ignored.
 */
typedef union
tagResampleTSFilterParams
{
  PassBandParamStruc    butterworth;	/**< A structure of type \c PassBandParamStruc used to store the parameters
                                         * of the butterworth filter used to perform low pass filtering
                                         */
  REAL8IIRFilter        iirfilter;	/**< A structure of type \c REAL8IIRFilter used to store the parameters of
                                         * the IIR or FIR filter used to perform low pass filtering
                                         */
}
ResampleTSFilterParams;

/**
 * This structure controls the behaviour of the resampling function.
 */
typedef struct
tagResampleTSParams
{
  REAL8                   deltaT;	/**< The sample interval desired in the down sampled time series */
  ResampleTSFilter        filterType;	/**< The type of filter with which to perform the low pass filtering */
  ResampleTSFilterParams  filterParams;	/**< Filter parameters for the low pass filter (Presently ignored) */
}
ResampleTSParams;

/*@}*/

/* ---------- Function prototypes ---------- */

int XLALResampleREAL4TimeSeries( REAL4TimeSeries *series, REAL8 dt );
int XLALResampleREAL8TimeSeries( REAL8TimeSeries *series, REAL8 dt );

void
LALResampleREAL4TimeSeries(
    LALStatus          *status,
    REAL4TimeSeries    *ts,
    ResampleTSParams   *params
    );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _RESAMPLETIMESERIES_H */

