/*
 * Copyright (C) 2015 Jolien Creighton
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

#ifndef _LALSIMSPHHARMMODE_H
#define _LALSIMSPHHARMMODE_H

#include <lal/LALDatatypes.h>
#include <lal/LALSimSphHarmSeries.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * @defgroup LALSimSphHarmMode_h Header LALSimSphHarmMode.h
 * @ingroup lalsimulation_general
 * @brief Routines to construct waveforms from spherical harmonic mode
 * decompositions.
 * @details
 * The gravitational wave polarizations depend on the location of the
 * observer relative to the source frame, defined by the spherical-polar
 * coordinates @p theta and @p phi (equivalently, inclination and azimuthal
 * phase).  Waveforms are sometimes decomposed into spin -2 weighted
 * spherical harmonics, indexed with mode quantum numbers @p l and @p m,
 * with the resulting waveform modes stored as COMPLEX16TimeSeries.
 * These routines reconstruct a waveform from these mode decompositions
 * for a given inclination and azimuthal phase.
 */

int XLALSimAddMode(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, COMPLEX16TimeSeries *hmode, REAL8 theta, REAL8 phi, int l, int m, int sym);
int XLALSimAddModeFromModes(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, SphHarmTimeSeries *hmode, REAL8 theta, REAL8 phi);
int XLALSimNewTimeSeriesFromModes(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, SphHarmTimeSeries *hmode, REAL8 theta, REAL8 phi);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMSPHHARMMODE_H */
