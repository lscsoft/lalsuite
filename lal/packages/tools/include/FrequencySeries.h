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


#ifndef _FREQUENCYSERIES_H
#define _FREQUENCYSERIES_H


#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/* COMPLEX8 prototypes */

void XLALDestroyCOMPLEX8FrequencySeries (
	COMPLEX8FrequencySeries *series
);

COMPLEX8FrequencySeries *XLALCreateCOMPLEX8FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);

COMPLEX8FrequencySeries *XLALCutCOMPLEX8FrequencySeries (
	const COMPLEX8FrequencySeries *series,
	size_t first,
	size_t length
);

COMPLEX8FrequencySeries *XLALResizeCOMPLEX8FrequencySeries (
	COMPLEX8FrequencySeries *series,
	int first,
	size_t length
);

COMPLEX8FrequencySeries *XLALShrinkCOMPLEX8FrequencySeries (
	COMPLEX8FrequencySeries *series,
	size_t first,
	size_t length
);

COMPLEX8FrequencySeries *XLALAddCOMPLEX8FrequencySeries (
	COMPLEX8FrequencySeries *arg1,
	const COMPLEX8FrequencySeries *arg2
);

COMPLEX8FrequencySeries *XLALConjugateCOMPLEX8FrequencySeries (
        COMPLEX8FrequencySeries *series
);

COMPLEX8FrequencySeries *XLALMultiplyCOMPLEX8FrequencySeries (
	COMPLEX8FrequencySeries *arg1,
	const COMPLEX8FrequencySeries *arg2
);

/* COMPLEX 16 prototypes */

void XLALDestroyCOMPLEX16FrequencySeries (
	COMPLEX16FrequencySeries *series
);

COMPLEX16FrequencySeries *XLALCreateCOMPLEX16FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);

COMPLEX16FrequencySeries *XLALCutCOMPLEX16FrequencySeries (
	const COMPLEX16FrequencySeries *series,
	size_t first,
	size_t length
);

COMPLEX16FrequencySeries *XLALResizeCOMPLEX16FrequencySeries (
	COMPLEX16FrequencySeries *series,
	int first,
	size_t length
);

COMPLEX16FrequencySeries *XLALShrinkCOMPLEX16FrequencySeries (
	COMPLEX16FrequencySeries *series,
	size_t first,
	size_t length
);

COMPLEX16FrequencySeries *XLALAddCOMPLEX16FrequencySeries (
	COMPLEX16FrequencySeries *arg1,
	const COMPLEX16FrequencySeries *arg2
);

COMPLEX16FrequencySeries *XLALConjugateCOMPLEX16FrequencySeries (
        COMPLEX16FrequencySeries *series
);

COMPLEX16FrequencySeries *XLALMultiplyCOMPLEX16FrequencySeries (
	COMPLEX16FrequencySeries *arg1,
	const COMPLEX16FrequencySeries *arg2
);

/* REAL4 prototypes */

void XLALDestroyREAL4FrequencySeries (
	REAL4FrequencySeries *series
);

REAL4FrequencySeries *XLALCreateREAL4FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);

REAL4FrequencySeries *XLALCutREAL4FrequencySeries (
	const REAL4FrequencySeries *series,
	size_t first,
	size_t length
);

REAL4FrequencySeries *XLALResizeREAL4FrequencySeries (
	REAL4FrequencySeries *series,
	int first,
	size_t length
);

REAL4FrequencySeries *XLALShrinkREAL4FrequencySeries (
	REAL4FrequencySeries *series,
	size_t first,
	size_t length
);

REAL4FrequencySeries *XLALAddREAL4FrequencySeries (
	REAL4FrequencySeries *arg1,
	const REAL4FrequencySeries *arg2
);

REAL4FrequencySeries *XLALMultiplyREAL4FrequencySeries (
	REAL4FrequencySeries *arg1,
	const REAL4FrequencySeries *arg2
);

/* REAL8 prototypes */

void XLALDestroyREAL8FrequencySeries (
	REAL8FrequencySeries *series
);

REAL8FrequencySeries *XLALCreateREAL8FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);

REAL8FrequencySeries *XLALCutREAL8FrequencySeries (
	const REAL8FrequencySeries *series,
	size_t first,
	size_t length
);

REAL8FrequencySeries *XLALResizeREAL8FrequencySeries (
	REAL8FrequencySeries *series,
	int first,
	size_t length
);

REAL8FrequencySeries *XLALShrinkREAL8FrequencySeries (
	REAL8FrequencySeries *series,
	size_t first,
	size_t length
);

REAL8FrequencySeries *XLALAddREAL8FrequencySeries (
	REAL8FrequencySeries *arg1,
	const REAL8FrequencySeries *arg2
);

REAL8FrequencySeries *XLALMultiplyREAL8FrequencySeries (
	REAL8FrequencySeries *arg1,
	const REAL8FrequencySeries *arg2
);

/* INT2 prototypes */

void XLALDestroyINT2FrequencySeries (
	INT2FrequencySeries *series
);


INT2FrequencySeries *XLALCreateINT2FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);


INT2FrequencySeries *XLALCutINT2FrequencySeries (
	const INT2FrequencySeries *series,
	size_t first,
	size_t length
);

INT2FrequencySeries *XLALResizeINT2FrequencySeries (
	INT2FrequencySeries *series,
	int first,
	size_t length
);

INT2FrequencySeries *XLALShrinkINT2FrequencySeries (
	INT2FrequencySeries *series,
	size_t first,
	size_t length
);

INT2FrequencySeries *XLALAddINT2FrequencySeries (
	INT2FrequencySeries *arg1,
	const INT2FrequencySeries *arg2
);

INT2FrequencySeries *XLALMultiplyINT2FrequencySeries (
	INT2FrequencySeries *arg1,
	const INT2FrequencySeries *arg2
);

/* INT4 prototypes */

void XLALDestroyINT4FrequencySeries (
	INT4FrequencySeries *series
);

INT4FrequencySeries *XLALCreateINT4FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);

INT4FrequencySeries *XLALCutINT4FrequencySeries (
	const INT4FrequencySeries *series,
	size_t first,
	size_t length
);

INT4FrequencySeries *XLALResizeINT4FrequencySeries (
	INT4FrequencySeries *series,
	int first,
	size_t length
);

INT4FrequencySeries *XLALShrinkINT4FrequencySeries (
	INT4FrequencySeries *series,
	size_t first,
	size_t length
);

INT4FrequencySeries *XLALAddINT4FrequencySeries (
	INT4FrequencySeries *arg1,
	const INT4FrequencySeries *arg2
);

INT4FrequencySeries *XLALMultiplyINT4FrequencySeries (
	INT4FrequencySeries *arg1,
	const INT4FrequencySeries *arg2
);

/* INT8 prototypes */

void XLALDestroyINT8FrequencySeries (
	INT8FrequencySeries *series
);

INT8FrequencySeries *XLALCreateINT8FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);

INT8FrequencySeries *XLALCutINT8FrequencySeries (
	const INT8FrequencySeries *series,
	size_t first,
	size_t length
);


INT8FrequencySeries *XLALResizeINT8FrequencySeries (
	INT8FrequencySeries *series,
	int first,
	size_t length
);

INT8FrequencySeries *XLALShrinkINT8FrequencySeries (
	INT8FrequencySeries *series,
	size_t first,
	size_t length
);

INT8FrequencySeries *XLALAddINT8FrequencySeries (
	INT8FrequencySeries *arg1,
	const INT8FrequencySeries *arg2
);

INT8FrequencySeries *XLALMultiplyINT8FrequencySeries (
	INT8FrequencySeries *arg1,
	const INT8FrequencySeries *arg2
);

/* UINT2 prototypes */

void XLALDestroyUINT2FrequencySeries (
	UINT2FrequencySeries *series
);

UINT2FrequencySeries *XLALCreateUINT2FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);

UINT2FrequencySeries *XLALCutUINT2FrequencySeries (
	const UINT2FrequencySeries *series,
	size_t first,
	size_t length
);

UINT2FrequencySeries *XLALResizeUINT2FrequencySeries (
	UINT2FrequencySeries *series,
	int first,
	size_t length
);

UINT2FrequencySeries *XLALShrinkUINT2FrequencySeries (
	UINT2FrequencySeries *series,
	size_t first,
	size_t length
);

UINT2FrequencySeries *XLALAddUINT2FrequencySeries (
	UINT2FrequencySeries *arg1,
	const UINT2FrequencySeries *arg2
);

UINT2FrequencySeries *XLALMultiplyUINT2FrequencySeries (
	UINT2FrequencySeries *arg1,
	const UINT2FrequencySeries *arg2
);

/* UINT4 prototypes */

void XLALDestroyUINT4FrequencySeries (
	UINT4FrequencySeries *series
);

UINT4FrequencySeries *XLALCreateUINT4FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);

UINT4FrequencySeries *XLALCutUINT4FrequencySeries (
	const UINT4FrequencySeries *series,
	size_t first,
	size_t length
);

UINT4FrequencySeries *XLALResizeUINT4FrequencySeries (
	UINT4FrequencySeries *series,
	int first,
	size_t length
);

UINT4FrequencySeries *XLALShrinkUINT4FrequencySeries (
	UINT4FrequencySeries *series,
	size_t first,
	size_t length
);

UINT4FrequencySeries *XLALAddUINT4FrequencySeries (
	UINT4FrequencySeries *arg1,
	const UINT4FrequencySeries *arg2
);

UINT4FrequencySeries *XLALMultiplyUINT4FrequencySeries (
	UINT4FrequencySeries *arg1,
	const UINT4FrequencySeries *arg2
);


/* UINT8 prototypes */

void XLALDestroyUINT8FrequencySeries (
	UINT8FrequencySeries *series
);

UINT8FrequencySeries *XLALCreateUINT8FrequencySeries (
	const CHAR *name,
	const LIGOTimeGPS *epoch,
	REAL8 f0,
	REAL8 deltaF,
	const LALUnit *sampleUnits,
	size_t length
);

UINT8FrequencySeries *XLALCutUINT8FrequencySeries (
	const UINT8FrequencySeries *series,
	size_t first,
	size_t length
);

UINT8FrequencySeries *XLALResizeUINT8FrequencySeries (
	UINT8FrequencySeries *series,
	int first,
	size_t length
);

UINT8FrequencySeries *XLALShrinkUINT8FrequencySeries (
	UINT8FrequencySeries *series,
	size_t first,
	size_t length
);

UINT8FrequencySeries *XLALAddUINT8FrequencySeries (
	UINT8FrequencySeries *arg1,
	const UINT8FrequencySeries *arg2
);

UINT8FrequencySeries *XLALMultiplyUINT8FrequencySeries (
	UINT8FrequencySeries *arg1,
	const UINT8FrequencySeries *arg2
);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif  /* _FREQUENCYSERIES_H */
