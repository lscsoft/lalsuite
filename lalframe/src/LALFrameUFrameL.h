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

#ifndef _LALFRAMEUFRAMEL_H
#define _LALFRAMEUFRAMEL_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

void XLALFrameUFrFileClose_FrameL_(LALFrameUFrFile * stream);
LALFrameUFrFile *XLALFrameUFrFileOpen_FrameL_(const char *filename, const char *mode);
int XLALFrameUFileCksumValid_FrameL_(LALFrameUFrFile * stream);
void XLALFrameUFrTOCFree_FrameL_(LALFrameUFrTOC * toc);
LALFrameUFrTOC *XLALFrameUFrTOCRead_FrameL_(LALFrameUFrFile * stream);
size_t XLALFrameUFrTOCQueryNFrame_FrameL_(const LALFrameUFrTOC * toc);
double XLALFrameUFrTOCQueryGTimeModf_FrameL_(double *iptr, const LALFrameUFrTOC * toc, size_t pos);
double XLALFrameUFrTOCQueryDt_FrameL_(const LALFrameUFrTOC * toc, size_t pos);
size_t XLALFrameUFrTOCQueryAdcN_FrameL_(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQueryAdcName_FrameL_(const LALFrameUFrTOC * toc, size_t adc);
size_t XLALFrameUFrTOCQuerySimN_FrameL_(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQuerySimName_FrameL_(const LALFrameUFrTOC * toc, size_t sim);
size_t XLALFrameUFrTOCQueryProcN_FrameL_(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQueryProcName_FrameL_(const LALFrameUFrTOC * toc, size_t proc);
size_t XLALFrameUFrTOCQueryDetectorN_FrameL_(const LALFrameUFrTOC * toc);
const char *XLALFrameUFrTOCQueryDetectorName_FrameL_(const LALFrameUFrTOC * toc, size_t det);
void XLALFrameUFrameHFree_FrameL_(LALFrameUFrameH * frame);
LALFrameUFrameH *XLALFrameUFrameHAlloc_FrameL_(const char *name, double start, double dt, int frnum);
LALFrameUFrameH *XLALFrameUFrameHRead_FrameL_(LALFrameUFrFile * stream, int pos);
int XLALFrameUFrameHWrite_FrameL_(LALFrameUFrFile * stream, LALFrameUFrameH * frame);
int XLALFrameUFrameHFrChanAdd_FrameL_(LALFrameUFrameH * frame, LALFrameUFrChan * channel);
int XLALFrameUFrameHFrDetectorAdd_FrameL_(LALFrameUFrameH * frame, LALFrameUFrDetector * detector);
int XLALFrameUFrameHFrHistoryAdd_FrameL_(LALFrameUFrameH * frame, LALFrameUFrHistory * history);
const char *XLALFrameUFrameHQueryName_FrameL_(const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryRun_FrameL_(const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryFrame_FrameL_(const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryDataQuality_FrameL_(const LALFrameUFrameH * frame);
double XLALFrameUFrameHQueryGTimeModf_FrameL_(double *iptr, const LALFrameUFrameH * frame);
int XLALFrameUFrameHQueryULeapS_FrameL_(const LALFrameUFrameH * frame);
double XLALFrameUFrameHQueryDt_FrameL_(const LALFrameUFrameH * frame);
int XLALFrameUFrameHSetRun_FrameL_(LALFrameUFrameH * frame, int run);
void XLALFrameUFrChanFree_FrameL_(LALFrameUFrChan * channel);
LALFrameUFrChan *XLALFrameUFrChanRead_FrameL_(LALFrameUFrFile * stream, const char *name, size_t pos);
LALFrameUFrChan *XLALFrameUFrAdcChanAlloc_FrameL_(const char *name, int dtype, size_t ndata);
LALFrameUFrChan *XLALFrameUFrSimChanAlloc_FrameL_(const char *name, int dtype, size_t ndata);
LALFrameUFrChan *XLALFrameUFrProcChanAlloc_FrameL_(const char *name, int type, int subtype, int dtype, size_t ndata);
const char *XLALFrameUFrChanQueryName_FrameL_(const LALFrameUFrChan * channel);
double XLALFrameUFrChanQueryTimeOffset_FrameL_(const LALFrameUFrChan * channel);
int XLALFrameUFrChanSetSampleRate_FrameL_(LALFrameUFrChan * channel, double sampleRate);
int XLALFrameUFrChanSetTimeOffset_FrameL_(LALFrameUFrChan * channel, double timeOffset);
int XLALFrameUFrChanVectorAlloc_FrameL_(LALFrameUFrChan * channel, int dtype, size_t ndata);
int XLALFrameUFrChanVectorCompress_FrameL_(LALFrameUFrChan * channel, int compressLevel);
int XLALFrameUFrChanVectorExpand_FrameL_(LALFrameUFrChan * channel);
const char *XLALFrameUFrChanVectorQueryName_FrameL_(const LALFrameUFrChan * channel);
int XLALFrameUFrChanVectorQueryCompress_FrameL_(const LALFrameUFrChan * channel);
int XLALFrameUFrChanVectorQueryType_FrameL_(const LALFrameUFrChan * channel);
void *XLALFrameUFrChanVectorQueryData_FrameL_(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNBytes_FrameL_(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNData_FrameL_(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNDim_FrameL_(const LALFrameUFrChan * channel);
size_t XLALFrameUFrChanVectorQueryNx_FrameL_(const LALFrameUFrChan * channel, size_t dim);
double XLALFrameUFrChanVectorQueryDx_FrameL_(const LALFrameUFrChan * channel, size_t dim);
double XLALFrameUFrChanVectorQueryStartX_FrameL_(const LALFrameUFrChan * channel, size_t dim);
const char *XLALFrameUFrChanVectorQueryUnitX_FrameL_(const LALFrameUFrChan * channel, size_t dim);
const char *XLALFrameUFrChanVectorQueryUnitY_FrameL_(const LALFrameUFrChan * channel);
int XLALFrameUFrChanVectorSetName_FrameL_(LALFrameUFrChan * channel, const char *name);
int XLALFrameUFrChanVectorSetDx_FrameL_(LALFrameUFrChan * channel, double dx);
int XLALFrameUFrChanVectorSetStartX_FrameL_(LALFrameUFrChan * channel, double x0);
int XLALFrameUFrChanVectorSetUnitX_FrameL_(LALFrameUFrChan * channel, const char *unit);
int XLALFrameUFrChanVectorSetUnitY_FrameL_(LALFrameUFrChan * channel, const char *unit);
void XLALFrameUFrDetectorFree_FrameL_(LALFrameUFrDetector * detector);
LALFrameUFrDetector *XLALFrameUFrDetectorRead_FrameL_(LALFrameUFrFile * stream, const char *name);
LALFrameUFrDetector *XLALFrameUFrDetectorAlloc_FrameL_(const char *name, const char *prefix, double latitude, double longitude, double elevation, double azimuthX, double azimuthY, double altitudeX, double altitudeY, double midpointX, double midpointY, int localTime);
const char *XLALFrameUFrDetectorQueryName_FrameL_(const LALFrameUFrDetector * detector);
const char *XLALFrameUFrDetectorQueryPrefix_FrameL_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryLongitude_FrameL_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryLatitude_FrameL_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryElevation_FrameL_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmXAzimuth_FrameL_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmYAzimuth_FrameL_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmXAltitude_FrameL_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmYAltitude_FrameL_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmXMidpoint_FrameL_(const LALFrameUFrDetector * detector);
double XLALFrameUFrDetectorQueryArmYMidpoint_FrameL_(const LALFrameUFrDetector * detector);
int XLALFrameUFrDetectorQueryLocalTime_FrameL_(const LALFrameUFrDetector * detector);
void XLALFrameUFrHistoryFree_FrameL_(LALFrameUFrHistory * history);
LALFrameUFrHistory *XLALFrameUFrHistoryAlloc_FrameL_(const char *name, double gpssec, const char *comment);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif /* _LALFRAMEUFRAMEL_H */
