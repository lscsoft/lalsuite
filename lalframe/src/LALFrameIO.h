/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Robert Adam Mercer, Xavier Siemens
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

#ifndef _LALFRAMEIO_H
#define _LALFRAMEIO_H

#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/LALFrameU.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif
struct tagLALFrFile;

/* alias some LALFrameU types */
typedef LALFrameUFrameH LALFrameH;
typedef struct tagLALFrFile LALFrFile;

int XLALFrFileClose(LALFrFile * frfile);
LALFrFile *XLALFrFileOpenURL(const char *url);
size_t XLALFrFileQueryNFrame(const LALFrFile * frfile);
LIGOTimeGPS *XLALFrFileQueryGTime(LIGOTimeGPS * start,
    const LALFrFile * frfile, size_t pos);
double XLALFrFileQueryDt(const LALFrFile * frfile, size_t pos);
LALTYPECODE XLALFrFileQueryChanType(const LALFrFile * frfile,
    const char *chname, size_t pos);
size_t XLALFrFileQueryChanVectorLength(const LALFrFile * frfile,
    const char *chname, size_t pos);
int XLALFrFileCksumValid(LALFrFile * frfile);

INT2TimeSeries *XLALFrFileReadINT2TimeSeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);
INT4TimeSeries *XLALFrFileReadINT4TimeSeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);
INT8TimeSeries *XLALFrFileReadINT8TimeSeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);
UINT2TimeSeries *XLALFrFileReadUINT2TimeSeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);
UINT4TimeSeries *XLALFrFileReadUINT4TimeSeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);
UINT8TimeSeries *XLALFrFileReadUINT8TimeSeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);
REAL4TimeSeries *XLALFrFileReadREAL4TimeSeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);
REAL8TimeSeries *XLALFrFileReadREAL8TimeSeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);
COMPLEX8TimeSeries *XLALFrFileReadCOMPLEX8TimeSeriesMetadata(LALFrFile *
    stream, const char *chname, size_t pos);
COMPLEX16TimeSeries *XLALFrFileReadCOMPLEX16TimeSeriesMetadata(LALFrFile *
    stream, const char *chname, size_t pos);
REAL4FrequencySeries *XLALFrFileReadREAL4FrequencySeriesMetadata(LALFrFile *
    stream, const char *chname, size_t pos);
REAL8FrequencySeries *XLALFrFileReadREAL8FrequencySeriesMetadata(LALFrFile *
    stream, const char *chname, size_t pos);
COMPLEX8FrequencySeries
    *XLALFrFileReadCOMPLEX8FrequencySeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);
COMPLEX16FrequencySeries
    *XLALFrFileReadCOMPLEX16FrequencySeriesMetadata(LALFrFile * stream,
    const char *chname, size_t pos);

/* routines to read channels from a frame at position pos in a frame file */
INT2TimeSeries *XLALFrFileReadINT2TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
INT4TimeSeries *XLALFrFileReadINT4TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
INT8TimeSeries *XLALFrFileReadINT8TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
UINT2TimeSeries *XLALFrFileReadUINT2TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
UINT4TimeSeries *XLALFrFileReadUINT4TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
UINT8TimeSeries *XLALFrFileReadUINT8TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
REAL4TimeSeries *XLALFrFileReadREAL4TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
REAL8TimeSeries *XLALFrFileReadREAL8TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
COMPLEX8TimeSeries *XLALFrFileReadCOMPLEX8TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
COMPLEX16TimeSeries *XLALFrFileReadCOMPLEX16TimeSeries(LALFrFile * stream,
    const char *chname, size_t pos);
REAL4FrequencySeries *XLALFrFileReadREAL4FrequencySeries(LALFrFile * stream,
    const char *chname, size_t pos);
REAL8FrequencySeries *XLALFrFileReadREAL8FrequencySeries(LALFrFile * stream,
    const char *chname, size_t pos);
COMPLEX8FrequencySeries *XLALFrFileReadCOMPLEX8FrequencySeries(LALFrFile *
    stream, const char *chname, size_t pos);
COMPLEX16FrequencySeries *XLALFrFileReadCOMPLEX16FrequencySeries(LALFrFile *
    stream, const char *chname, size_t pos);

int XLALFrameAddFrHistory(LALFrameH * frame, const char *name,
    const char *comment);
int XLALFrameAddFrDetector(LALFrameH * frame, const LALFrDetector * detector);
void XLALFrameFree(LALFrameH * frame);
LALFrameH *XLALFrameNew(const LIGOTimeGPS * epoch, double duration,
    const char *project, int run, int frnum, int detectorFlags);

int XLALFrameAddINT2TimeSeriesAdcData(LALFrameH * frame,
    const INT2TimeSeries * series);
int XLALFrameAddINT4TimeSeriesAdcData(LALFrameH * frame,
    const INT4TimeSeries * series);
int XLALFrameAddREAL4TimeSeriesAdcData(LALFrameH * frame,
    const REAL4TimeSeries * series);
int XLALFrameAddREAL8TimeSeriesAdcData(LALFrameH * frame,
    const REAL8TimeSeries * series);

int XLALFrameAddINT2TimeSeriesSimData(LALFrameH * frame,
    const INT2TimeSeries * series);
int XLALFrameAddINT4TimeSeriesSimData(LALFrameH * frame,
    const INT4TimeSeries * series);
int XLALFrameAddREAL4TimeSeriesSimData(LALFrameH * frame,
    const REAL4TimeSeries * series);
int XLALFrameAddREAL8TimeSeriesSimData(LALFrameH * frame,
    const REAL8TimeSeries * series);

int XLALFrameAddINT2TimeSeriesProcData(LALFrameH * frame,
    const INT2TimeSeries * series);
int XLALFrameAddINT4TimeSeriesProcData(LALFrameH * frame,
    const INT4TimeSeries * series);
int XLALFrameAddINT8TimeSeriesProcData(LALFrameH * frame,
    const INT8TimeSeries * series);
int XLALFrameAddUINT2TimeSeriesProcData(LALFrameH * frame,
    const UINT2TimeSeries * series);
int XLALFrameAddUINT4TimeSeriesProcData(LALFrameH * frame,
    const UINT4TimeSeries * series);
int XLALFrameAddUINT8TimeSeriesProcData(LALFrameH * frame,
    const UINT8TimeSeries * series);
int XLALFrameAddREAL4TimeSeriesProcData(LALFrameH * frame,
    const REAL4TimeSeries * series);
int XLALFrameAddREAL8TimeSeriesProcData(LALFrameH * frame,
    const REAL8TimeSeries * series);
int XLALFrameAddCOMPLEX8TimeSeriesProcData(LALFrameH * frame,
    const COMPLEX8TimeSeries * series);
int XLALFrameAddCOMPLEX16TimeSeriesProcData(LALFrameH * frame,
    const COMPLEX16TimeSeries * series);

int XLALFrameAddREAL4FrequencySeriesProcData(LALFrameH * frame,
    const REAL4FrequencySeries * series, int subtype);
int XLALFrameAddREAL8FrequencySeriesProcData(LALFrameH * frame,
    const REAL8FrequencySeries * series, int subtype);
int XLALFrameAddCOMPLEX8FrequencySeriesProcData(LALFrameH * frame,
    const COMPLEX8FrequencySeries * series, int subtype);
int XLALFrameAddCOMPLEX16FrequencySeriesProcData(LALFrameH * frame,
    const COMPLEX16FrequencySeries * series, int subtype);

/* frame writing function */

int XLALFrameWrite(LALFrameH * frame, const char *fname);

/* direct output functions */

int XLALFrWriteINT2TimeSeries(const INT2TimeSeries * series, int fnum);
int XLALFrWriteINT4TimeSeries(const INT4TimeSeries * series, int fnum);
int XLALFrWriteINT8TimeSeries(const INT8TimeSeries * series, int fnum);
int XLALFrWriteREAL4TimeSeries(const REAL4TimeSeries * series, int fnum);
int XLALFrWriteREAL8TimeSeries(const REAL8TimeSeries * series, int fnum);
int XLALFrWriteCOMPLEX8TimeSeries(const COMPLEX8TimeSeries * series,
    int fnum);
int XLALFrWriteCOMPLEX16TimeSeries(const COMPLEX16TimeSeries * series,
    int fnum);

int XLALFrWriteREAL4FrequencySeries(const REAL4FrequencySeries * series,
    int fnum, int subtype);
int XLALFrWriteREAL8FrequencySeries(const REAL8FrequencySeries * series,
    int fnum, int subtype);
int XLALFrWriteCOMPLEX8FrequencySeries(const COMPLEX8FrequencySeries *
    series, int fnum, int subtype);
int XLALFrWriteCOMPLEX16FrequencySeries(const COMPLEX16FrequencySeries *
    series, int fnum, int subtype);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
