/*
*  Copyright (C) 2007 Jolien Creighton
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

#ifndef SEGMENT_H
#define SEGMENT_H

/*
 *
 * Routine to create a single overwhitened data segment from a time series:
 * this performs the FFT to get the data into the frequency domain and then
 * multiplies it by the inverse power spectrum to overwhiten it.
 *
 */

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>

int compute_data_segment(
    COMPLEX8FrequencySeries  *segment,
    UINT4                     segmentNumber,
    REAL4TimeSeries          *series,
    REAL4FrequencySeries     *invspec,
    COMPLEX8FrequencySeries  *response,
    REAL8                     segmentDuration,
    REAL8                     strideDuration,
    REAL4FFTPlan             *fwdPlan
    );

#endif /* SEGMENT_H */
