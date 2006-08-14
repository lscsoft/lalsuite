/* vim: set noet ts=4 sw=4: */
/** \file
 * \ingroup std
 * \author Creighton, J. D. E.
 * \date $Date$
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

NRCSID( AUDIOH, "$Id$" );

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

/** Records a time series as a .wav audio file */
int XLALAudioWAVRecordREAL4TimeSeries( FILE *fp, REAL4TimeSeries *series );

/** Records a time series as a .au audio file */
int XLALAudioAURecordREAL4TimeSeries( FILE *fp, REAL4TimeSeries *series );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* AUDIO_H */
