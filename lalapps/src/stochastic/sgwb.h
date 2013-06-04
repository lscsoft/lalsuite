/*
 * sgwb.h - SGWB Standalone Analysis Pipeline
 *        - Main Stochastic Search Function Prototypes
 *
 * Copyright (C) 2002-2006,2010 Adam Mercer
 * Copyright (C) 2003-2004 Tania Regimbau
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 *
 */

#ifndef SGWB_H
#define SGWB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <getopt.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/FrameCalibration.h>
#include <lal/LALFrStream.h>
#include <lal/LALStdio.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/PrintFTSeries.h>

#include <lalapps.h>

/* generate a data window */
REAL4Window *data_window(REAL8 delta_t,
    INT4 length,
    INT4 hann_duration);

/* return the spectrum */
REAL4FrequencySeries *omega_gw(LALStatus *status,
    REAL4 exponent,
    REAL8 f_ref,
    REAL4 omega_ref,
    UINT4 length,
    REAL8 f0,
    REAL8 delta_f);

/* return the overlap reduction function */
REAL4FrequencySeries *overlap_reduction_function(LALStatus *status,
    UINT4 length,
    REAL8 f0,
    REAL8 delta_f,
    INT4 site_one,
    INT4 site_two);

/* calculate and return the inverse noise */
REAL4FrequencySeries *inverse_noise(LALStatus *status,
    REAL4FrequencySeries *psd,
    COMPLEX8FrequencySeries *response);

/* calculate and return the optimal filter */
REAL4FrequencySeries *optimal_filter(LALStatus *status,
    REAL4FrequencySeries *overlap,
    REAL4FrequencySeries *omega,
    REAL4FrequencySeries *psd_one,
    REAL4FrequencySeries *psd_two,
    REAL4Window *window,
    REAL8 *sigma,
    REAL8 f_ref,
    INT4 segment_duration);

/* estimate the PSD */
REAL4FrequencySeries *estimate_psd(REAL4TimeSeries *series,
    REAL8 f0,
    INT4 shrink_length,
    INT4 psd_window_duration);

/* return a unity response function, for use with calibrated data */
COMPLEX8FrequencySeries *unity_response(LIGOTimeGPS epoch,
    REAL8 f0,
    REAL8 delta_f,
    LALUnit units,
    INT4 length);

/* generate response function for LIGO data */
COMPLEX8FrequencySeries *ligo_response(LALStatus *status,
    CHAR *ifo,
    CHAR *cache_file,
    LIGOTimeGPS epoch,
    REAL8 f0,
    REAL8 delta_f,
    LALUnit units,
    INT4 length,
    INT4 offset);

/* wrapper to unity_response() and ligo_response() for generating the
 * appropriate response for the given detector */
COMPLEX8FrequencySeries *generate_response(LALStatus *status,
    CHAR *ifo,
    CHAR *cache_file,
    LIGOTimeGPS epoch,
    REAL8 f0,
    REAL8 delta_f,
    LALUnit units,
    INT4 length,
    INT4 offset);

/* return the frequency mask */
REAL4FrequencySeries *frequency_mask(REAL8 f0,
    REAL8 delta_f,
    INT4 length,
    INT4 bins);

/* zero pad and fft */
COMPLEX8FrequencySeries *zero_pad_and_fft(LALStatus *status,
    REAL4TimeSeries *series,
    REAL8 delta_f,
    INT4 length,
    REAL4Window *window);

/* generate the cross correlation spectra */
COMPLEX8FrequencySeries *cc_spectrum(LALStatus *status,
    COMPLEX8FrequencySeries *zero_pad_one,
    COMPLEX8FrequencySeries *zero_pad_two,
    COMPLEX8FrequencySeries *response_one,
    COMPLEX8FrequencySeries *response_two,
    REAL4FrequencySeries *opt_filter);

/* generate cross correlation statistic from cross correlation spectra */
REAL8 cc_statistic(COMPLEX8FrequencySeries *cc_spectra);

/* main stochastic search */
StochasticTable *stochastic_search(LALStatus *status,
    REAL4TimeSeries *series_one,
    REAL4TimeSeries *series_two,
    REAL4FrequencySeries *overlap,
    REAL4FrequencySeries *omega,
    REAL4Window *window,
    INT4 num_segments,
    INT4 filter_length,
    INT4 segment_duration,
    INT4 segs_in_interval,
    INT4 segment_length,
    INT4 psd_window_duration,
    CHAR *ifo_one,
    CHAR *ifo_two,
    CHAR *channel_one,
    CHAR *channel_two,
    CHAR *calibration_cache_one,
    CHAR *calibration_cache_two,
    INT4 calibration_offset,
    REAL8 f_min,
    REAL8 f_max,
    REAL8 f_ref);

#endif /* SGWB_H */

/*
 * vim: et
 */
