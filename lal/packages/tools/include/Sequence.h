/*
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


#ifndef _SEQUENCE_H
#define _SEQUENCE_H


#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/* COMPLEX8 prototypes */

void XLALDestroyCOMPLEX8Sequence (
	COMPLEX8Sequence *sequence
);

COMPLEX8Sequence *XLALCreateCOMPLEX8Sequence (
	size_t length
);

COMPLEX8Sequence *XLALCutCOMPLEX8Sequence (
	COMPLEX8Sequence *sequence,
	size_t first,
	size_t length
);

COMPLEX8Sequence *XLALCopyCOMPLEX8Sequence (
	COMPLEX8Sequence *sequence
);

void XLALShiftCOMPLEX8Sequence (
	COMPLEX8Sequence *sequence,
	int count
);

COMPLEX8Sequence *XLALResizeCOMPLEX8Sequence (
	COMPLEX8Sequence *sequence,
	int first,
	size_t length
);

COMPLEX8Sequence *XLALShrinkCOMPLEX8Sequence (
	COMPLEX8Sequence *sequence,
	size_t first,
	size_t length
);

COMPLEX8 XLALCOMPLEX8Sum (
	const COMPLEX8 *data,
	size_t first,
	size_t count
);

REAL4 XLALCOMPLEX8SumSquares (
	const COMPLEX8 *data,
	size_t first,
	size_t count
);

COMPLEX8 XLALCOMPLEX8SequenceSum (
	const COMPLEX8Sequence *sequence,
	size_t first,
	size_t count
);

REAL4 XLALCOMPLEX8SequenceSumSquares (
	const COMPLEX8Sequence *sequence,
	size_t first,
	size_t count
);

COMPLEX8Sequence *XLALConjugateCOMPLEX8Sequence (
        COMPLEX8Sequence *series
);

/* COMPLEX16 prototypes */

void XLALDestroyCOMPLEX16Sequence (
	COMPLEX16Sequence *sequence
);

COMPLEX16Sequence *XLALCreateCOMPLEX16Sequence (
	size_t length
);

COMPLEX16Sequence *XLALCutCOMPLEX16Sequence (
	COMPLEX16Sequence *sequence,
	size_t first,
	size_t length
);

COMPLEX16Sequence *XLALCopyCOMPLEX16Sequence (
	COMPLEX16Sequence *sequence
);

void XLALShiftCOMPLEX16Sequence (
	COMPLEX16Sequence *sequence,
	int count
);

COMPLEX16Sequence *XLALResizeCOMPLEX16Sequence (
	COMPLEX16Sequence *sequence,
	int first,
	size_t length
);

COMPLEX16Sequence *XLALShrinkCOMPLEX16Sequence (
	COMPLEX16Sequence *sequence,
	size_t first,
	size_t length
);

COMPLEX16 XLALCOMPLEX16Sum (
	const COMPLEX16 *data,
	size_t first,
	size_t count
);

REAL8 XLALCOMPLEX16SumSquares (
	const COMPLEX16 *data,
	size_t first,
	size_t count
);

COMPLEX16 XLALCOMPLEX16SequenceSum (
	const COMPLEX16Sequence *sequence,
	size_t first,
	size_t count
);

REAL8 XLALCOMPLEX16SequenceSumSquares (
	const COMPLEX16Sequence *sequence,
	size_t first,
	size_t count
);

COMPLEX16Sequence *XLALConjugateCOMPLEX16Sequence (
        COMPLEX16Sequence *series
);

/* REAL4 prototypes */

void XLALDestroyREAL4Sequence (
	REAL4Sequence *sequence
);

REAL4Sequence *XLALCreateREAL4Sequence (
	size_t length
);

REAL4Sequence *XLALCutREAL4Sequence (
	REAL4Sequence *sequence,
	size_t first,
	size_t length
);

REAL4Sequence *XLALCopyREAL4Sequence (
	REAL4Sequence *sequence
);

void XLALShiftREAL4Sequence (
	REAL4Sequence *sequence,
	int count
);

REAL4Sequence *XLALResizeREAL4Sequence (
	REAL4Sequence *sequence,
	int first,
	size_t length
);

REAL4Sequence *XLALShrinkREAL4Sequence (
	REAL4Sequence *sequence,
	size_t first,
	size_t length
);

REAL4 XLALREAL4Sum (
	const REAL4 *data,
	size_t first,
	size_t count
);

REAL4 XLALREAL4SumSquares (
	const REAL4 *data,
	size_t first,
	size_t count
);

REAL4 XLALREAL4SequenceSum (
	const REAL4Sequence *sequence,
	size_t first,
	size_t count
);

REAL4 XLALREAL4SequenceSumSquares (
	const REAL4Sequence *sequence,
	size_t first,
	size_t count
);

/* REAL8 prototypes */

void XLALDestroyREAL8Sequence (
	REAL8Sequence *sequence
);

REAL8Sequence *XLALCreateREAL8Sequence (
	size_t length
);

REAL8Sequence *XLALCutREAL8Sequence (
	REAL8Sequence *sequence,
	size_t first,
	size_t length
);

REAL8Sequence *XLALCopyREAL8Sequence (
	REAL8Sequence *sequence
);

void XLALShiftREAL8Sequence (
	REAL8Sequence *sequence,
	int count
);

REAL8Sequence *XLALResizeREAL8Sequence (
	REAL8Sequence *sequence,
	int first,
	size_t length
);

REAL8Sequence *XLALShrinkREAL8Sequence (
	REAL8Sequence *sequence,
	size_t first,
	size_t length
);

REAL8 XLALREAL8Sum (
	const REAL8 *data,
	size_t first,
	size_t count
);

REAL8 XLALREAL8SumSquares (
	const REAL8 *data,
	size_t first,
	size_t count
);

REAL8 XLALREAL8SequenceSum (
	const REAL8Sequence *sequence,
	size_t first,
	size_t count
);

REAL8 XLALREAL8SequenceSumSquares (
	const REAL8Sequence *sequence,
	size_t first,
	size_t count
);

/* INT2 prototypes */

void XLALDestroyINT2Sequence (
	INT2Sequence *sequence
);

INT2Sequence *XLALCreateINT2Sequence (
	size_t length
);

INT2Sequence *XLALCutINT2Sequence (
	INT2Sequence *sequence,
	size_t first,
	size_t length
);

INT2Sequence *XLALCopyINT2Sequence (
	INT2Sequence *sequence
);

void XLALShiftINT2Sequence (
	INT2Sequence *sequence,
	int count
);

INT2Sequence *XLALResizeINT2Sequence (
	INT2Sequence *sequence,
	int first,
	size_t length
);

INT2Sequence *XLALShrinkINT2Sequence (
	INT2Sequence *sequence,
	size_t first,
	size_t length
);

INT2 XLALINT2Sum (
	const INT2 *data,
	size_t first,
	size_t count
);

UINT2 XLALINT2SumSquares (
	const INT2 *data,
	size_t first,
	size_t count
);

INT2 XLALINT2SequenceSum (
	const INT2Sequence *sequence,
	size_t first,
	size_t count
);

UINT2 XLALINT2SequenceSumSquares (
	const INT2Sequence *sequence,
	size_t first,
	size_t count
);

/* INT4 prototypes */

void XLALDestroyINT4Sequence (
	INT4Sequence *sequence
);

INT4Sequence *XLALCreateINT4Sequence (
	size_t length
);

INT4Sequence *XLALCutINT4Sequence (
	INT4Sequence *sequence,
	size_t first,
	size_t length
);

INT4Sequence *XLALCopyINT4Sequence (
	INT4Sequence *sequence
);

void XLALShiftINT4Sequence (
	INT4Sequence *sequence,
	int count
);

INT4Sequence *XLALResizeINT4Sequence (
	INT4Sequence *sequence,
	int first,
	size_t length
);

INT4Sequence *XLALShrinkINT4Sequence (
	INT4Sequence *sequence,
	size_t first,
	size_t length
);

INT4 XLALINT4Sum (
	const INT4 *data,
	size_t first,
	size_t count
);

UINT4 XLALINT4SumSquares (
	const INT4 *data,
	size_t first,
	size_t count
);

INT4 XLALINT4SequenceSum (
	const INT4Sequence *sequence,
	size_t first,
	size_t count
);

UINT4 XLALINT4SequenceSumSquares (
	const INT4Sequence *sequence,
	size_t first,
	size_t count
);

/* INT8 prototypes */

void XLALDestroyINT8Sequence (
	INT8Sequence *sequence
);

INT8Sequence *XLALCreateINT8Sequence (
	size_t length
);

INT8Sequence *XLALCutINT8Sequence (
	INT8Sequence *sequence,
	size_t first,
	size_t length
);

INT8Sequence *XLALCopyINT8Sequence (
	INT8Sequence *sequence
);

void XLALShiftINT8Sequence (
	INT8Sequence *sequence,
	int count
);

INT8Sequence *XLALResizeINT8Sequence (
	INT8Sequence *sequence,
	int first,
	size_t length
);

INT8Sequence *XLALShrinkINT8Sequence (
	INT8Sequence *sequence,
	size_t first,
	size_t length
);

INT8 XLALINT8Sum (
	const INT8 *data,
	size_t first,
	size_t count
);

UINT8 XLALINT8SumSquares (
	const INT8 *data,
	size_t first,
	size_t count
);

INT8 XLALINT8SequenceSum (
	const INT8Sequence *sequence,
	size_t first,
	size_t count
);

UINT8 XLALINT8SequenceSumSquares (
	const INT8Sequence *sequence,
	size_t first,
	size_t count
);

/* UINT2 prototypes */

void XLALDestroyUINT2Sequence (
	UINT2Sequence *sequence
);

UINT2Sequence *XLALCreateUINT2Sequence (
	size_t length
);

UINT2Sequence *XLALCutUINT2Sequence (
	UINT2Sequence *sequence,
	size_t first,
	size_t length
);

UINT2Sequence *XLALCopyUINT2Sequence (
	UINT2Sequence *sequence
);

void XLALShiftUINT2Sequence (
	UINT2Sequence *sequence,
	int count
);

UINT2Sequence *XLALResizeUINT2Sequence (
	UINT2Sequence *sequence,
	int first,
	size_t length
);

UINT2Sequence *XLALShrinkUINT2Sequence (
	UINT2Sequence *sequence,
	size_t first,
	size_t length
);

UINT2 XLALUINT2Sum (
	const UINT2 *data,
	size_t first,
	size_t count
);

UINT2 XLALUINT2SumSquares (
	const UINT2 *data,
	size_t first,
	size_t count
);

UINT2 XLALUINT2SequenceSum (
	const UINT2Sequence *sequence,
	size_t first,
	size_t count
);

UINT2 XLALUINT2SequenceSumSquares (
	const UINT2Sequence *sequence,
	size_t first,
	size_t count
);

/* UINT4 prototypes */

void XLALDestroyUINT4Sequence (
	UINT4Sequence *sequence
);

UINT4Sequence *XLALCreateUINT4Sequence (
	size_t length
);

UINT4Sequence *XLALCutUINT4Sequence (
	UINT4Sequence *sequence,
	size_t first,
	size_t length
);

UINT4Sequence *XLALCopyUINT4Sequence (
	UINT4Sequence *sequence
);

void XLALShiftUINT4Sequence (
	UINT4Sequence *sequence,
	int count
);

UINT4Sequence *XLALResizeUINT4Sequence (
	UINT4Sequence *sequence,
	int first,
	size_t length
);

UINT4Sequence *XLALShrinkUINT4Sequence (
	UINT4Sequence *sequence,
	size_t first,
	size_t length
);

UINT4 XLALUINT4Sum (
	const UINT4 *data,
	size_t first,
	size_t count
);

UINT4 XLALUINT4SumSquares (
	const UINT4 *data,
	size_t first,
	size_t count
);

UINT4 XLALUINT4SequenceSum (
	const UINT4Sequence *sequence,
	size_t first,
	size_t count
);

UINT4 XLALUINT4SequenceSumSquares (
	const UINT4Sequence *sequence,
	size_t first,
	size_t count
);

/* UINT8 prototypes */

void XLALDestroyUINT8Sequence (
	UINT8Sequence *sequence
);

UINT8Sequence *XLALCreateUINT8Sequence (
	size_t length
);

UINT8Sequence *XLALCutUINT8Sequence (
	UINT8Sequence *sequence,
	size_t first,
	size_t length
);

UINT8Sequence *XLALCopyUINT8Sequence (
	UINT8Sequence *sequence
);

void XLALShiftUINT8Sequence (
	UINT8Sequence *sequence,
	int count
);

UINT8Sequence *XLALResizeUINT8Sequence (
	UINT8Sequence *sequence,
	int first,
	size_t length
);

UINT8Sequence *XLALShrinkUINT8Sequence (
	UINT8Sequence *sequence,
	size_t first,
	size_t length
);

UINT8 XLALUINT8Sum (
	const UINT8 *data,
	size_t first,
	size_t count
);

UINT8 XLALUINT8SumSquares (
	const UINT8 *data,
	size_t first,
	size_t count
);

UINT8 XLALUINT8SequenceSum (
	const UINT8Sequence *sequence,
	size_t first,
	size_t count
);

UINT8 XLALUINT8SequenceSumSquares (
	const UINT8Sequence *sequence,
	size_t first,
	size_t count
);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif  /* _SEQUENCE_H */
