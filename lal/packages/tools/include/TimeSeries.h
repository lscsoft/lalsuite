/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#ifndef _TIMESERIES_H
#define _TIMESERIES_H


#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/* COMPLEX8 prototypes */

void XLALDestroyCOMPLEX8TimeSeries (
	COMPLEX8TimeSeries *series
);

COMPLEX8TimeSeries *XLALCreateCOMPLEX8TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

COMPLEX8TimeSeries *XLALCutCOMPLEX8TimeSeries (
	const COMPLEX8TimeSeries *series,
	size_t first,
	size_t length
);

COMPLEX8TimeSeries *XLALResizeCOMPLEX8TimeSeries (
	COMPLEX8TimeSeries *series,
	int first,
	size_t length
);

COMPLEX8TimeSeries *XLALShrinkCOMPLEX8TimeSeries (
	COMPLEX8TimeSeries *series,
	size_t first,
	size_t length
);

COMPLEX8TimeSeries *XLALAddCOMPLEX8TimeSeries (
	COMPLEX8TimeSeries *arg1,
	const COMPLEX8TimeSeries *arg2
);

/* COMPLEX16 prototypes */

void XLALDestroyCOMPLEX16TimeSeries (
	COMPLEX16TimeSeries *series
);

COMPLEX16TimeSeries *XLALCreateCOMPLEX16TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

COMPLEX16TimeSeries *XLALCutCOMPLEX16TimeSeries (
	const COMPLEX16TimeSeries *series,
	size_t first,
	size_t length
);

COMPLEX16TimeSeries *XLALResizeCOMPLEX16TimeSeries (
	COMPLEX16TimeSeries *series,
	int first,
	size_t length
);

COMPLEX16TimeSeries *XLALShrinkCOMPLEX16TimeSeries (
	COMPLEX16TimeSeries *series,
	size_t first,
	size_t length
);

COMPLEX16TimeSeries *XLALAddCOMPLEX16TimeSeries (
	COMPLEX16TimeSeries *arg1,
	const COMPLEX16TimeSeries *arg2
);

/* REAL4 prototypes */

void XLALDestroyREAL4TimeSeries (
	REAL4TimeSeries *series
);

REAL4TimeSeries *XLALCreateREAL4TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

REAL4TimeSeries *XLALCutREAL4TimeSeries (
	const REAL4TimeSeries *series,
	size_t first,
	size_t length
);

REAL4TimeSeries *XLALResizeREAL4TimeSeries (
	REAL4TimeSeries *series,
	int first,
	size_t length
);

REAL4TimeSeries *XLALShrinkREAL4TimeSeries (
	REAL4TimeSeries *series,
	size_t first,
	size_t length
);

REAL4TimeSeries *XLALAddREAL4TimeSeries (
	REAL4TimeSeries *arg1,
	const REAL4TimeSeries *arg2
);

/* REAL8 prototypes */

void XLALDestroyREAL8TimeSeries (
	REAL8TimeSeries *series
);

REAL8TimeSeries *XLALCreateREAL8TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

REAL8TimeSeries *XLALCutREAL8TimeSeries (
	const REAL8TimeSeries *series,
	size_t first,
	size_t length
);

REAL8TimeSeries *XLALResizeREAL8TimeSeries (
	REAL8TimeSeries *series,
	int first,
	size_t length
);

REAL8TimeSeries *XLALShrinkREAL8TimeSeries (
	REAL8TimeSeries *series,
	size_t first,
	size_t length
);

REAL8TimeSeries *XLALAddREAL8TimeSeries (
	REAL8TimeSeries *arg1,
	const REAL8TimeSeries *arg2
);

/* INT2 prototypes */

void XLALDestroyINT2TimeSeries (
	INT2TimeSeries *series
);

INT2TimeSeries *XLALCreateINT2TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

INT2TimeSeries *XLALCutINT2TimeSeries (
	const INT2TimeSeries *series,
	size_t first,
	size_t length
);

INT2TimeSeries *XLALResizeINT2TimeSeries (
	INT2TimeSeries *series,
	int first,
	size_t length
);

INT2TimeSeries *XLALShrinkINT2TimeSeries (
	INT2TimeSeries *series,
	size_t first,
	size_t length
);

INT2TimeSeries *XLALAddINT2TimeSeries (
	INT2TimeSeries *arg1,
	const INT2TimeSeries *arg2
);

/* INT4 prototypes */

void XLALDestroyINT4TimeSeries (
	INT4TimeSeries *series
);

INT4TimeSeries *XLALCreateINT4TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

INT4TimeSeries *XLALCutINT4TimeSeries (
	const INT4TimeSeries *series,
	size_t first,
	size_t length
);

INT4TimeSeries *XLALResizeINT4TimeSeries (
	INT4TimeSeries *series,
	int first,
	size_t length
);

INT4TimeSeries *XLALShrinkINT4TimeSeries (
	INT4TimeSeries *series,
	size_t first,
	size_t length
);

INT4TimeSeries *XLALAddINT4TimeSeries (
	INT4TimeSeries *arg1,
	const INT4TimeSeries *arg2
);

/* INT8 prototypes */

void XLALDestroyINT8TimeSeries (
	INT8TimeSeries *series
);

INT8TimeSeries *XLALCreateINT8TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

INT8TimeSeries *XLALCutINT8TimeSeries (
	const INT8TimeSeries *series,
	size_t first,
	size_t length
);

INT8TimeSeries *XLALResizeINT8TimeSeries (
	INT8TimeSeries *series,
	int first,
	size_t length
);

INT8TimeSeries *XLALShrinkINT8TimeSeries (
	INT8TimeSeries *series,
	size_t first,
	size_t length
);

INT8TimeSeries *XLALAddINT8TimeSeries (
	INT8TimeSeries *arg1,
	const INT8TimeSeries *arg2
);

/* UINT2 prototypes */

void XLALDestroyUINT2TimeSeries (
	UINT2TimeSeries *series
);

UINT2TimeSeries *XLALCreateUINT2TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

UINT2TimeSeries *XLALCutUINT2TimeSeries (
	const UINT2TimeSeries *series,
	size_t first,
	size_t length
);

UINT2TimeSeries *XLALResizeUINT2TimeSeries (
	UINT2TimeSeries *series,
	int first,
	size_t length
);

UINT2TimeSeries *XLALShrinkUINT2TimeSeries (
	UINT2TimeSeries *series,
	size_t first,
	size_t length
);

UINT2TimeSeries *XLALAddUINT2TimeSeries (
	UINT2TimeSeries *arg1,
	const UINT2TimeSeries *arg2
);

/* UINT4 prototypes */

void XLALDestroyUINT4TimeSeries (
	UINT4TimeSeries *series
);

UINT4TimeSeries *XLALCreateUINT4TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

UINT4TimeSeries *XLALCutUINT4TimeSeries (
	const UINT4TimeSeries *series,
	size_t first,
	size_t length
);

UINT4TimeSeries *XLALResizeUINT4TimeSeries (
	UINT4TimeSeries *series,
	int first,
	size_t length
);

UINT4TimeSeries *XLALShrinkUINT4TimeSeries (
	UINT4TimeSeries *series,
	size_t first,
	size_t length
);

UINT4TimeSeries *XLALAddUINT4TimeSeries (
	UINT4TimeSeries *arg1,
	const UINT4TimeSeries *arg2
);

/* UINT8 prototypes */

void XLALDestroyUINT8TimeSeries (
	UINT8TimeSeries *series
);

UINT8TimeSeries *XLALCreateUINT8TimeSeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaT,
	const LALUnit *sampleUnits,
	size_t length
);

UINT8TimeSeries *XLALCutUINT8TimeSeries (
	const UINT8TimeSeries *series,
	size_t first,
	size_t length
);

UINT8TimeSeries *XLALResizeUINT8TimeSeries (
	UINT8TimeSeries *series,
	int first,
	size_t length
);

UINT8TimeSeries *XLALShrinkUINT8TimeSeries (
	UINT8TimeSeries *series,
	size_t first,
	size_t length
);

UINT8TimeSeries *XLALAddUINT8TimeSeries (
	UINT8TimeSeries *arg1,
	const UINT8TimeSeries *arg2
);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif  /* _TIMESERIES_H */
