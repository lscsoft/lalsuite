/*
*  Copyright (C) 2014 Jolien Creighton
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

#ifndef _LALFRAMEUFRAMEC_H
#define _LALFRAMEUFRAMEC_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

void XLALFrameUFrFileClose_FrameC_(LALFrameUFrFile * stream);
LALFrameUFrFile *XLALFrameUFrFileOpen_FrameC_(const char *filename, const char *mode);
int XLALFrameUFileCksumValid_FrameC_(LALFrameUFrFile * stream);
void XLALFrameUFrTOCFree_FrameC_(LALFrameUFrTOC * toc);
LALFrameUFrTOC *XLALFrameUFrTOCRead_FrameC_(LALFrameUFrFile * stream);
size_t XLALFrameUFrTOCQueryNFrame_FrameC_(const LALFrameUFrTOC * toc);
double XLALFrameUFrTOCQueryGTimeModf_FrameC_(double *iptr, const LALFrameUFrTOC * toc, size_t pos);
double XLALFrameUFrTOCQueryDt_FrameC_(const LALFrameUFrTOC * toc, size_t pos);
size_t XLALFrameUFrTOCQueryAdcN_FrameC_(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQueryAdcName_FrameC_(const LALFrameUFrTOC * toc, size_t adc);
size_t XLALFrameUFrTOCQuerySimN_FrameC_(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQuerySimName_FrameC_(const LALFrameUFrTOC * toc, size_t sim);
size_t XLALFrameUFrTOCQueryProcN_FrameC_(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQueryProcName_FrameC_(const LALFrameUFrTOC * toc, size_t proc);
size_t XLALFrameUFrTOCQueryDetectorN_FrameC_(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQueryDetectorName_FrameC_(const LALFrameUFrTOC * toc, size_t det);
void XLALFrameUFrameHFree_FrameC_(LALFrameUFrameH * frame);
LALFrameUFrameH *XLALFrameUFrameHAlloc_FrameC_(const char *name, double start, double dt, int frnum);
LALFrameUFrameH *XLALFrameUFrameHRead_FrameC_(LALFrameUFrFile * stream, int pos);
int XLALFrameUFrameHWrite_FrameC_(LALFrameUFrFile * stream, LALFrameUFrameH * frame);
int XLALFrameUFrameHFrChanAdd_FrameC_(LALFrameUFrameH * frame, LALFrameUFrChan * channel);
int XLALFrameUFrameHFrDetectorAdd_FrameC_(LALFrameUFrameH * frame, LALFrameUFrDetector * detector);
int XLALFrameUFrameHFrHistoryAdd_FrameC_(LALFrameUFrameH * frame, LALFrameUFrHistory * history);
const char *XLALFrameUFrameHQueryName_FrameC_(const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryRun_FrameC_(const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryFrame_FrameC_(const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryDataQuality_FrameC_(const LALFrameUFrameH * frame);
double XLALFrameUFrameHQueryGTimeModf_FrameC_(double *iptr, const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryULeapS_FrameC_(const LALFrameUFrameH * frame);
double XLALFrameUFrameHQueryDt_FrameC_(const LALFrameUFrameH * frame);
int XLALFrameUFrameHSetRun_FrameC_(LALFrameUFrameH * frame, int run);
void XLALFrameUFrChanFree_FrameC_(LALFrameUFrChan * channel);
LALFrameUFrChan *XLALFrameUFrChanRead_FrameC_(LALFrameUFrFile * stream, const char *name, size_t pos);
LALFrameUFrChan *XLALFrameUFrAdcChanAlloc_FrameC_(const char *name, int dtype, size_t ndata);
LALFrameUFrChan *XLALFrameUFrSimChanAlloc_FrameC_(const char *name, int dtype, size_t ndata);
LALFrameUFrChan *XLALFrameUFrProcChanAlloc_FrameC_(const char *name, int type, int subtype, int dtype, size_t ndata);
const char *XLALFrameUFrChanQueryName_FrameC_(const LALFrameUFrChan * channel);
double XLALFrameUFrChanQueryTimeOffset_FrameC_(const LALFrameUFrChan * channel);
int XLALFrameUFrChanSetSampleRate_FrameC_(LALFrameUFrChan * channel, double sampleRate);
int XLALFrameUFrChanSetTimeOffset_FrameC_(LALFrameUFrChan * channel, double timeOffset);
int XLALFrameUFrChanVectorAlloc_FrameC_(LALFrameUFrChan * channel, int dtype, size_t ndata);
int XLALFrameUFrChanVectorCompress_FrameC_(LALFrameUFrChan * channel, int compressLevel);
int XLALFrameUFrChanVectorExpand_FrameC_(LALFrameUFrChan * channel);
const char *XLALFrameUFrChanVectorQueryName_FrameC_(const LALFrameUFrChan * channel);
int XLALFrameUFrChanVectorQueryCompress_FrameC_(const LALFrameUFrChan * channel);
int XLALFrameUFrChanVectorQueryType_FrameC_(const LALFrameUFrChan * channel);
void *XLALFrameUFrChanVectorQueryData_FrameC_(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNBytes_FrameC_(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNData_FrameC_(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNDim_FrameC_(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNx_FrameC_(const LALFrameUFrChan * channel, size_t dim);
double XLALFrameUFrChanVectorQueryDx_FrameC_(const LALFrameUFrChan * channel, size_t dim);
double XLALFrameUFrChanVectorQueryStartX_FrameC_(const LALFrameUFrChan * channel, size_t dim);
const char *XLALFrameUFrChanVectorQueryUnitX_FrameC_(const LALFrameUFrChan * channel, size_t dim);
const char *XLALFrameUFrChanVectorQueryUnitY_FrameC_(const LALFrameUFrChan * channel);
int XLALFrameUFrChanVectorSetName_FrameC_(LALFrameUFrChan * channel, const char *name);
int XLALFrameUFrChanVectorSetDx_FrameC_(LALFrameUFrChan * channel, double dx);
int XLALFrameUFrChanVectorSetStartX_FrameC_(LALFrameUFrChan * channel, double x0);
int XLALFrameUFrChanVectorSetUnitX_FrameC_(LALFrameUFrChan * channel, const char *unit);
int XLALFrameUFrChanVectorSetUnitY_FrameC_(LALFrameUFrChan * channel, const char *unit);
void XLALFrameUFrDetectorFree_FrameC_(LALFrameUFrDetector * detector);
LALFrameUFrDetector *XLALFrameUFrDetectorRead_FrameC_(LALFrameUFrFile * stream, const char *name);
LALFrameUFrDetector *XLALFrameUFrDetectorAlloc_FrameC_(const char *name, const char *prefix, double latitude, double longitude, double elevation, double azimuthX, double azimuthY, double altitudeX, double altitudeY, double midpointX, double midpointY, int localTime);
const char *XLALFrameUFrDetectorQueryName_FrameC_(const LALFrameUFrDetector * detector);
const char *XLALFrameUFrDetectorQueryPrefix_FrameC_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryLongitude_FrameC_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryLatitude_FrameC_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryElevation_FrameC_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmXAzimuth_FrameC_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmYAzimuth_FrameC_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmXAltitude_FrameC_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmYAltitude_FrameC_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmXMidpoint_FrameC_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmYMidpoint_FrameC_(const LALFrameUFrDetector * detector);
int XLALFrameUFrDetectorQueryLocalTime_FrameC_(const LALFrameUFrDetector * detector);
void XLALFrameUFrHistoryFree_FrameC_(LALFrameUFrHistory * history);
LALFrameUFrHistory *XLALFrameUFrHistoryAlloc_FrameC_(const char *name, double gpssec, const char *comment);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif /* _LALFRAMEUFRAMEC_H */
