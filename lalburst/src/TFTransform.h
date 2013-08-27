/*
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#ifndef _TFTRANSFORM_H
#define _TFTRANSFORM_H


#include <gsl/gsl_matrix.h>


#include <lal/LALDatatypes.h>
#include <lal/Window.h>
#include <lal/RealFFT.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Sequence.h>


#ifdef  __cplusplus		/* C++ protection. */
extern "C" {
#endif

/** \defgroup TFTransform_h Header TFTransform.h
    \ingroup pkg_burstsearch

Provides routines to to compute time-frequency planes from either
time-domain or frequency-domain data, for use in the excess
power search technique.

\heading{Synopsis}
\code
#include "TFTransform.h"
\endcode

This package provides a suite for functions for computing time-frequency
representations of data using stacked Fourier transforms.

The first few functions simply combine functionality from the packages
\c fft and \c Window in a convenient way.  They are designed to
streamline the task of setting up structures to prepare for taking many
discrete Fourier transforms (DFTs), including windows and plans for FFTW.

A general description of the time-frequency (TF) transform provided by
TFTransform is as follows.  Suppose one starts with some data \f$h_j\f$, \f$0 \le j < n\f$
in the time domain, with sampling time \f$\Delta t\f$, so that the data point
\f$h_j\f$ corresponds to a time \f$t_j = t_\textrm{start} + j \Delta t\f$.  Taking the
standard DFT yields complex data
\anchor standarddft
\f{equation}{
{\tilde h}_\gamma = \sum_{j=0}^{n-1} \, e^{-2 \pi i j \gamma / n} \, h_j
\tag{standarddft}
\f}
in the Fourier domain, for \f$0 \le \gamma \le [n/2]+1\f$.  Here the data point
\f${\tilde h}_\gamma\f$ corresponds to a frequency \f$f_\gamma = \gamma \Delta f\f$,
where \f$\Delta f= 1/(n \Delta t)\f$ is the frequency resolution.


Now suppose that we can factorize the number \f$n\f$ of data points as
\f{equation}{
n = 2 N_T N_F.
\f}
Then, by a time-frequency plane we shall mean a set of \f$N_T N_F\f$ complex
numbers \f$H_{I\Gamma}\f$ with \f$0 \le I < N_T\f$ and \f$0 \le \Gamma < N_F\f$, obtained
by an invertible linear transformation from the original data, such  that the
data point \f$H_{I\Gamma}\f$ corresponds approximately to a time \f$t_I = t_\textrm{
start} + I {\overline {\Delta t}}\f$ and to a frequency \f$f_\Gamma = \Gamma
{\overline {\Delta f}}\f$.  Here \f$N_F\f$ is the number of frequency bins in the TF
plane, and \f$N_T\f$ is the number of time bins.  The time resolution \f${\overline
{\Delta t}}\f$ and frequency resolution \f${\overline {\Delta f}}\f$ are related by
\f${\overline {\Delta t}} \ {\overline {\Delta f}} =1\f$, and are given by
\f${\overline {\Delta t}} = 2 N_F \Delta t\f$ and \f${\overline {\Delta f}} = N_T
\Delta f\f$.  Note that there are many other time-frequency representations
of data that are not of this type; see [\ref ab1999].


There are many possible choices of linear transformations from the data \f$h_j\f$
to data \f$H_{J\Gamma}\f$ satisfying the above properties.  Here we have
implemented two simple choices.  The first choice consists of dividing the
time-domain data \f$h_j\f$ into \f$N_T\f$ equal-sized chunks, each of length \f$n/N_T\f$,
and then taking the forward DFT of each chunk.  Then, \f$H_{J\Gamma}\f$ is just
the \f$\Gamma\f$th element of the \f$J\f$th chunk.  In terms of formulae this
corresponds to
\anchor verticalTFP \f{equation}{
H_{J\Sigma} = \sum_{k=0}^{2 N_F-1} \, \exp \left[ 2 \pi i k \Sigma / (2
N_F) \right] \, h_{2 N_F J + k},
\tag{verticalTFP}
\f}
for \f$0 \le J < N_T\f$ and \f$0 \le \Sigma < N_F\f$.  We call this first type
of TF plane a vertical TF plane, since it corresponds to a series of
vertical lines if the time axis is horizontal and the frequency axis
vertical.

The second type of TF plane is obtained by first taking a DFT of all the time
domain data to obtain frequency domain data, then dividing the frequency
domain data into \f$N_F\f$ equal-sized chunks, then taking the inverse DFT of each
chunk.  We call the resulting TF plane a horizontal TF plane. In terms of
formulae the TF plane elements are \anchor horizontalTFP
\f{equation}{
H_{J\Sigma} =
\sum_{\gamma=0}^{N_T-1} \, \exp \left[ -2 \pi i J \gamma / N_T \right] \,
{\tilde h}_{N_T \Sigma + \gamma}, \tag{horizontalTFP}
\f}
for \f$0 \le J < N_T\f$ and \f$0 \le \Sigma < N_F\f$, where \f${\tilde h}_\gamma\f$ is given by
Eq. \eqref{standarddft}.

*/
/*@{*/


/* ---------- Data types ---------- */

/**
 * A time-frequency plane
 */
typedef struct tagREAL8TimeFrequencyPlaneTiles {
  unsigned max_length;
  unsigned min_channels;
  unsigned max_channels;
  unsigned tiling_start;
  unsigned tiling_end;
  unsigned inv_fractional_stride;
  double dof_per_pixel;
} REAL8TimeFrequencyPlaneTiles;
typedef struct tagREAL8TimeFrequencyPlane {
  CHAR name[LALNameLength];	/**< name of data from which this was computed */
  LIGOTimeGPS epoch;		/**< epoch of data from which this was computed */
  REAL8 deltaT;			/**< time resolution of the plane */
  REAL8 fseries_deltaF;		/**< input frequency series' resolution */
  REAL8 deltaF;			/**< TF plane's frequency resolution (channel spacing) */
  REAL8 flow;			/**< low frequency boundary of TF plane */
  gsl_matrix *channel_data;   	/**< channel data.  each channel is placed into its own column.
                                 * channel_data[i * channels + j] corresponds to time
                                 *
                                 * epoch + i * deltaT
                                 *
                                 * and the frequency band
                                 *
                                 * [flow + j * deltaF, flow + (j + 1) * deltaF)
                                 */

  REAL8Sequence *channel_buffer;/**< re-usable holding area for the data for a single channel */
  REAL8Sequence *unwhitened_channel_buffer;	/**< UNDOCUMENTED */
  REAL8TimeFrequencyPlaneTiles tiles;	  	/**< time-frequency plane's tiling information */

	REAL8Window *window;	/**< time-domain window applied to input time series for tapering edges to 0 */
	INT4 window_shift;	/**< by how many samples a window's start should be shifted from the start of the window preceding it */
	REAL8Sequence *two_point_spectral_correlation;		/**< two-point spectral correlation of the whitened frequency series, computed from the time-domain window function */
} REAL8TimeFrequencyPlane;

/** UNDOCUMENTED */
typedef struct tagExcessPowerFilter {
  COMPLEX16FrequencySeries *fseries;
  REAL8 unwhitened_rms;			/**< root mean square of the unwhitened time series corresponding to this filter */
} ExcessPowerFilter;
typedef struct tagLALExcessPowerFilterBank {
  int n_filters;
  ExcessPowerFilter *basis_filters;

  REAL8Sequence *twice_channel_overlap;	  /**< twice the inner product of filters for neighbouring channels;
                                           * twice_channel_overlap[0] is twice the inner product of the
                                           * filters for channels 0 and 1, and so on (for n channels, there
                                           * are n - 1 channel_overlaps) */
  REAL8Sequence *unwhitened_cross;	/**< the mean square cross terms for wide channels (indices same as for twice_channel_overlap) */
} LALExcessPowerFilterBank;


/* ---------- Function prototypes ---------- */

double XLALExcessPowerFilterInnerProduct(
	const COMPLEX16FrequencySeries *filter1,
	const COMPLEX16FrequencySeries *filter2,
	const REAL8Sequence *correlation,
	const REAL8FrequencySeries *psd
);


COMPLEX16FrequencySeries *XLALCreateExcessPowerFilter(
	REAL8 channel_flow,
	REAL8 channel_width,
	const REAL8FrequencySeries *psd,
	const REAL8Sequence *correlation
);


LALExcessPowerFilterBank *XLALCreateExcessPowerFilterBank(
	double filter_deltaF,
	double flow,
	double channel_bandwidth,
	int n_channels,
	const REAL8FrequencySeries *psd,
	const REAL8Sequence *two_point_spectral_correlation
);


void XLALDestroyExcessPowerFilterBank(
	LALExcessPowerFilterBank *bank
);


REAL8Sequence *XLALREAL8WindowTwoPointSpectralCorrelation(
	const REAL8Window *window,
	const REAL8FFTPlan *plan
);


REAL8TimeFrequencyPlane *XLALCreateTFPlane(
	UINT4 tseries_length,
	REAL8 tseries_deltaT,
	REAL8 flow,
	REAL8 bandwidth,
	REAL8 tiling_fractional_stride,
	REAL8 tiling_max_bandwidth,
	REAL8 tiling_max_duration,
	const REAL8FFTPlan *plan
);


void XLALDestroyTFPlane(
	REAL8TimeFrequencyPlane *plane
);


int XLALFreqSeriesToTFPlane(
	REAL8TimeFrequencyPlane *tfplane,
	const LALExcessPowerFilterBank *filter_bank,
	const COMPLEX16FrequencySeries *fseries,
	const REAL8FFTPlan *reverseplan
);


SnglBurst *XLALComputeExcessPower(
	const REAL8TimeFrequencyPlane *plane,
	const LALExcessPowerFilterBank *filter_bank,
	SnglBurst *head,
	double confidence_threshold
);


INT4 XLALOverlappedSegmentsCommensurate(
	INT4 target_length,
	INT4 segment_length,
	INT4 segment_shift
);


INT4 XLALEPGetTimingParameters(
	INT4 window_length,
	INT4 max_tile_length,
	REAL8 fractional_tile_stride,
	INT4 *psd_length,
	INT4 *psd_shift,
	INT4 *window_shift,
	INT4 *window_pad,
	INT4 *tiling_length
);

/*@}*/

#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
