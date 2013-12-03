/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

#ifndef _STREAMOUTPUT_H
#define _STREAMOUTPUT_H

#include <lal/LALStdlib.h>
#include <lal/Grid.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \addtogroup StreamOutput_h
 * \author Creighton, T. D.
 *
 * \brief Provides routines to write data from LAL data structures to an open stream.
 *
 * ### Synopsis ###
 *
 * \code
 * #include "StreamOutput.h"
 * \endcode
 *
 * This header provides prototypes for routines that write the
 * contents of LAL time/frequency series structures or grid structures to
 * a file (or other I/O) stream, in a standard format.  The routines do
 * not provide a system-level interface to create files and open or close
 * file streams; they simply assume that they have been passed an open,
 * writeable stream.  Nonetheless, because they involve I/O stream
 * manipulation, these routines are placed in the \c lalsupport
 * library rather than in \c lal proper.
 *
 */
/*@{*/

/** \name Error Codes */ /*@{*/
#define STREAMOUTPUTH_ENUL 1    /**< Unexpected null pointer in arguments */
#define STREAMOUTPUTH_EPRN 2    /**< Print statement failed */
/*@}*/
/*@}*/

#define STREAMOUTPUTH_MSGENUL "Unexpected null pointer in arguments"
#define STREAMOUTPUTH_MSGEPRN "Print statement failed"


/* Function prototypes. */
void
LALI2WriteTSeries( LALStatus  *status, FILE *stream, INT2TimeSeries *series );
void
LALI4WriteTSeries( LALStatus  *status, FILE *stream, INT4TimeSeries *series );
void
LALI8WriteTSeries( LALStatus  *status, FILE *stream, INT8TimeSeries *series );
void
LALU2WriteTSeries( LALStatus  *status, FILE *stream, UINT2TimeSeries *series );
void
LALU4WriteTSeries( LALStatus  *status, FILE *stream, UINT4TimeSeries *series );
void
LALU8WriteTSeries( LALStatus  *status, FILE *stream, UINT8TimeSeries *series );
void
LALSWriteTSeries( LALStatus  *status, FILE *stream, REAL4TimeSeries *series );
void
LALDWriteTSeries( LALStatus  *status, FILE *stream, REAL8TimeSeries *series );
void
LALCWriteTSeries( LALStatus  *status, FILE *stream, COMPLEX8TimeSeries *series );
void
LALZWriteTSeries( LALStatus  *status, FILE *stream, COMPLEX16TimeSeries *series );

void
LALI2WriteTVectorSeries( LALStatus  *status, FILE *stream, INT2TimeVectorSeries *series );
void
LALI4WriteTVectorSeries( LALStatus  *status, FILE *stream, INT4TimeVectorSeries *series );
void
LALI8WriteTVectorSeries( LALStatus  *status, FILE *stream, INT8TimeVectorSeries *series );
void
LALU2WriteTVectorSeries( LALStatus  *status, FILE *stream, UINT2TimeVectorSeries *series );
void
LALU4WriteTVectorSeries( LALStatus  *status, FILE *stream, UINT4TimeVectorSeries *series );
void
LALU8WriteTVectorSeries( LALStatus  *status, FILE *stream, UINT8TimeVectorSeries *series );
void
LALSWriteTVectorSeries( LALStatus  *status, FILE *stream, REAL4TimeVectorSeries *series );
void
LALDWriteTVectorSeries( LALStatus  *status, FILE *stream, REAL8TimeVectorSeries *series );
void
LALCWriteTVectorSeries( LALStatus  *status, FILE *stream, COMPLEX8TimeVectorSeries *series );
void
LALZWriteTVectorSeries( LALStatus  *status, FILE *stream, COMPLEX16TimeVectorSeries *series );

void
LALI2WriteTArraySeries( LALStatus  *status, FILE *stream, INT2TimeArraySeries *series );
void
LALI4WriteTArraySeries( LALStatus  *status, FILE *stream, INT4TimeArraySeries *series );
void
LALI8WriteTArraySeries( LALStatus  *status, FILE *stream, INT8TimeArraySeries *series );
void
LALU2WriteTArraySeries( LALStatus  *status, FILE *stream, UINT2TimeArraySeries *series );
void
LALU4WriteTArraySeries( LALStatus  *status, FILE *stream, UINT4TimeArraySeries *series );
void
LALU8WriteTArraySeries( LALStatus  *status, FILE *stream, UINT8TimeArraySeries *series );
void
LALSWriteTArraySeries( LALStatus  *status, FILE *stream, REAL4TimeArraySeries *series );
void
LALDWriteTArraySeries( LALStatus  *status, FILE *stream, REAL8TimeArraySeries *series );
void
LALCWriteTArraySeries( LALStatus  *status, FILE *stream, COMPLEX8TimeArraySeries *series );
void
LALZWriteTArraySeries( LALStatus  *status, FILE *stream, COMPLEX16TimeArraySeries *series );

void
LALI2WriteFSeries( LALStatus  *status, FILE *stream, INT2FrequencySeries *series );
void
LALI4WriteFSeries( LALStatus  *status, FILE *stream, INT4FrequencySeries *series );
void
LALI8WriteFSeries( LALStatus  *status, FILE *stream, INT8FrequencySeries *series );
void
LALU2WriteFSeries( LALStatus  *status, FILE *stream, UINT2FrequencySeries *series );
void
LALU4WriteFSeries( LALStatus  *status, FILE *stream, UINT4FrequencySeries *series );
void
LALU8WriteFSeries( LALStatus  *status, FILE *stream, UINT8FrequencySeries *series );
void
LALSWriteFSeries( LALStatus  *status, FILE *stream, REAL4FrequencySeries *series );
void
LALDWriteFSeries( LALStatus  *status, FILE *stream, REAL8FrequencySeries *series );
void
LALCWriteFSeries( LALStatus  *status, FILE *stream, COMPLEX8FrequencySeries *series );
void
LALZWriteFSeries( LALStatus  *status, FILE *stream, COMPLEX16FrequencySeries *series );




void
LALI2WriteGrid( LALStatus  *status, FILE *stream, INT2Grid *grid );
void
LALI4WriteGrid( LALStatus  *status, FILE *stream, INT4Grid *grid );
void
LALI8WriteGrid( LALStatus  *status, FILE *stream, INT8Grid *grid );
void
LALU2WriteGrid( LALStatus  *status, FILE *stream, UINT2Grid *grid );
void
LALU4WriteGrid( LALStatus  *status, FILE *stream, UINT4Grid *grid );
void
LALU8WriteGrid( LALStatus  *status, FILE *stream, UINT8Grid *grid );
void
LALSWriteGrid( LALStatus  *status, FILE *stream, REAL4Grid *grid );
void
LALDWriteGrid( LALStatus  *status, FILE *stream, REAL8Grid *grid );
void
LALCWriteGrid( LALStatus  *status, FILE *stream, COMPLEX8Grid *grid );
void
LALZWriteGrid( LALStatus  *status, FILE *stream, COMPLEX16Grid *grid );

#ifdef __cplusplus
}
#endif

#endif /* _STREAMOUTPUT_H */
