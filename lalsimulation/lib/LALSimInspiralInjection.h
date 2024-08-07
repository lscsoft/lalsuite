/*
*  Copyright (C) 2023 Jolien Creighton
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef _LALSIMINSPIRALINJECTION_H
#define _LALSIMINSPIRALINJECTION_H

#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/LALDict.h>
#include <lal/LALDictSequence.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

LALDictSequence * XLALSimInspiralInjectionSequenceFromH5File(const char *fname);
int XLALSimInspiralInjectionSequenceToH5File(const LALDictSequence *injseq, const char *fname);

LIGOTimeGPS * XLALSimInspiralInjectionEndTime(LIGOTimeGPS *epoch, LALDict *injparams);
LIGOTimeGPS * XLALSimInspiralInjectionStartTime(LIGOTimeGPS *epoch, LALDict *injparams);

int XLALSimInspiralInjectionSequenceIsEndTimeOrdered(LALDictSequence *injseq);
int XLALSimInspiralInjectionSequenceIsStartTimeOrdered(LALDictSequence *injseq);
int XLALSimInspiralInjectionSequenceOrderByEndTime(LALDictSequence *injseq);
int XLALSimInspiralInjectionSequenceOrderByStartTime(LALDictSequence *injseq);
LALDictSequence * XLALSimInspiralInjectionSequenceInInterval(const LALDictSequence *injseq, const LIGOTimeGPS *start, const LIGOTimeGPS *end);

int XLALSimInspiralInjectionTDWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LALDict *injparams, REAL8 deltaT);

REAL8TimeSeries * XLALSimInspiralInjectionStrain(LALDict *injparams, REAL8 deltaT, const LALDetector *detector);

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LALSIMINSPIRALINJECTION_H */
