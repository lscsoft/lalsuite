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

/* vim: set noet ts=4 sw=4: */
/** \file
 * \ingroup std
 * \author Creighton, J. D. E.
 * \brief Routines for exporting time series data as sound files.
 *
 * Supported formats are AU and WAVE.  AIFF is not supported at this time.
 *
 * This code is based on code written by Paul Bourke.  See:
 * http://local.wasp.uwa.edu.au/~pbourke/dataformats/au/
 * http://local.wasp.uwa.edu.au/~pbourke/dataformats/wave/
 *
 */

#ifndef AUDIO_H
#define AUDIO_H

#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/** Records a time series as a .wav audio file */
int XLALAudioWAVRecordREAL4TimeSeries( FILE *fp, REAL4TimeSeries *series );

/** Records a time series as a .wav audio file */
int XLALAudioWAVRecordREAL8TimeSeries( FILE *fp, REAL8TimeSeries *series );

/** Records a time series as a .au audio file */
int XLALAudioAURecordREAL4TimeSeries( FILE *fp, REAL4TimeSeries *series );

/** Records a time series as a .au audio file */
int XLALAudioAURecordREAL8TimeSeries( FILE *fp, REAL8TimeSeries *series );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* AUDIO_H */
