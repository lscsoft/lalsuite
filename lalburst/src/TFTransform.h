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


#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>


#ifdef  __cplusplus		/* C++ protection. */
extern "C" {
#endif

/**
 * \defgroup TFTransform_h Header TFTransform.h
 * \ingroup lalburst_burstsearch
 *
 * Provides routines to to compute time-frequency planes from either
 * time-domain or frequency-domain data, for use in the excess
 * power search technique.
 *
 * ### Synopsis ###
 *
 * \code
 * #include "TFTransform.h"
 * \endcode
 *
 * This package provides a suite for functions for computing time-frequency
 * representations of data using stacked Fourier transforms.
 *
 * The first few functions simply combine functionality from the packages
 * \c fft and \c Window in a convenient way.  They are designed to
 * streamline the task of setting up structures to prepare for taking many
 * discrete Fourier transforms (DFTs), including windows and plans for FFTW.
 *
 * A general description of the time-frequency (TF) transform provided by
 * TFTransform is as follows.  Suppose one starts with some data \f$h_j\f$, \f$0 \le j < n\f$
 * in the time domain, with sampling time \f$\Delta t\f$, so that the data point
 * \f$h_j\f$ corresponds to a time \f$t_j = t_\textrm{start} + j \Delta t\f$.  Taking the
 * standard DFT yields complex data
 * \f{equation}{
 * \label{standarddft}
 * {\tilde h}_\gamma = \sum_{j=0}^{n-1} \, e^{-2 \pi i j \gamma / n} \, h_j
 * \f}
 * in the Fourier domain, for \f$0 \le \gamma \le [n/2]+1\f$.  Here the data point
 * \f${\tilde h}_\gamma\f$ corresponds to a frequency \f$f_\gamma = \gamma \Delta f\f$,
 * where \f$\Delta f= 1/(n \Delta t)\f$ is the frequency resolution.
 *
 * Now suppose that we can factorize the number \f$n\f$ of data points as
 * \f{equation}{
 * n = 2 N_T N_F.
 * \f}
 * Then, by a time-frequency plane we shall mean a set of \f$N_T N_F\f$ complex
 * numbers \f$H_{I\Gamma}\f$ with \f$0 \le I < N_T\f$ and \f$0 \le \Gamma < N_F\f$, obtained
 * by an invertible linear transformation from the original data, such  that the
 * data point \f$H_{I\Gamma}\f$ corresponds approximately to a time \f$t_I = t_\textrm{
 * start} + I {\overline {\Delta t}}\f$ and to a frequency \f$f_\Gamma = \Gamma
 * {\overline {\Delta f}}\f$.  Here \f$N_F\f$ is the number of frequency bins in the TF
 * plane, and \f$N_T\f$ is the number of time bins.  The time resolution \f${\overline
 * {\Delta t}}\f$ and frequency resolution \f${\overline {\Delta f}}\f$ are related by
 * \f${\overline {\Delta t}} \ {\overline {\Delta f}} =1\f$, and are given by
 * \f${\overline {\Delta t}} = 2 N_F \Delta t\f$ and \f${\overline {\Delta f}} = N_T
 * \Delta f\f$.  Note that there are many other time-frequency representations
 * of data that are not of this type; see \cite ab1999 .
 *
 * There are many possible choices of linear transformations from the data \f$h_j\f$
 * to data \f$H_{J\Gamma}\f$ satisfying the above properties.  Here we have
 * implemented two simple choices.  The first choice consists of dividing the
 * time-domain data \f$h_j\f$ into \f$N_T\f$ equal-sized chunks, each of length \f$n/N_T\f$,
 * and then taking the forward DFT of each chunk.  Then, \f$H_{J\Gamma}\f$ is just
 * the \f$\Gamma\f$th element of the \f$J\f$th chunk.  In terms of formulae this
 * corresponds to
 * \f{equation}{
 * \label{verticalTFP}
 * H_{J\Sigma} = \sum_{k=0}^{2 N_F-1} \, \exp \left[ 2 \pi i k \Sigma / (2
 * N_F) \right] \, h_{2 N_F J + k},
 * \f}
 * for \f$0 \le J < N_T\f$ and \f$0 \le \Sigma < N_F\f$.  We call this first type
 * of TF plane a vertical TF plane, since it corresponds to a series of
 * vertical lines if the time axis is horizontal and the frequency axis
 * vertical.
 *
 * The second type of TF plane is obtained by first taking a DFT of all the time
 * domain data to obtain frequency domain data, then dividing the frequency
 * domain data into \f$N_F\f$ equal-sized chunks, then taking the inverse DFT of each
 * chunk.  We call the resulting TF plane a horizontal TF plane. In terms of
 * formulae the TF plane elements are
 * \f{equation}{
 * \label{horizontalTFP}
 * H_{J\Sigma} =
 * \sum_{\gamma=0}^{N_T-1} \, \exp \left[ -2 \pi i J \gamma / N_T \right] \,
 * {\tilde h}_{N_T \Sigma + \gamma},
 * \f}
 * for \f$0 \le J < N_T\f$ and \f$0 \le \Sigma < N_F\f$, where \f${\tilde h}_\gamma\f$ is given by
 * \eqref{standarddft}.
 *
 */
/*@{*/


/* ---------- Function prototypes ---------- */

double XLALExcessPowerFilterInnerProduct(
	const COMPLEX16FrequencySeries *filter1,
	const COMPLEX16FrequencySeries *filter2,
	const REAL8Sequence *correlation,
	const REAL8FrequencySeries *psd
);


COMPLEX16FrequencySeries *XLALCreateExcessPowerFilter(
	double channel_flow,
	double channel_width,
	const REAL8FrequencySeries *psd,
	const REAL8Sequence *correlation
);


int XLALOverlappedSegmentsCommensurate(
	int target_length,
	int segment_length,
	int segment_shift
);


#ifdef SWIG /* SWIG interface directives */
SWIGLAL(INOUT_SCALARS(int*, psd_length));
#endif
int XLALEPGetTimingParameters(
	int window_length,
	int max_tile_length,
	double fractional_tile_stride,
	int *psd_length,
	int *psd_shift,
	int *window_shift,
	int *window_pad,
	int *tiling_length
);

/*@}*/

#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
