/*
*  Copyright (C) 2016 Jolien Creighton
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

#ifndef _LAL_VALUE_H
#define _LAL_VALUE_H

#include <stddef.h>
#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

struct tagLALValue;
typedef struct tagLALValue LALValue;

LALValue *XLALValueAlloc(size_t size);
LALValue *XLALValueRealloc(LALValue *value, size_t size);
LALValue *XLALValueDuplicate(const LALValue *value);
LALValue *XLALValueCopy(LALValue *copy, const LALValue *orig);
LALValue *XLALValueSet(LALValue *value, const void *data, size_t size, LALTYPECODE type);

void XLALDestroyValue(LALValue *value);

LALValue *XLALCreateValue(const void * data, size_t size, LALTYPECODE type);
LALValue *XLALCreateStringValue(const char *value);
LALValue *XLALCreateCHARValue(CHAR value);
LALValue *XLALCreateINT2Value(INT2 value);
LALValue *XLALCreateINT4Value(INT4 value);
LALValue *XLALCreateINT8Value(INT8 value);
LALValue *XLALCreateUCHARValue(UCHAR value);
LALValue *XLALCreateUINT2Value(UINT2 value);
LALValue *XLALCreateUINT4Value(UINT4 value);
LALValue *XLALCreateUINT8Value(UINT8 value);
LALValue *XLALCreateREAL4Value(REAL4 value);
LALValue *XLALCreateREAL8Value(REAL8 value);
LALValue *XLALCreateCOMPLEX8Value(COMPLEX8 value);
LALValue *XLALCreateCOMPLEX16Value(COMPLEX16 value);

LALTYPECODE XLALValueGetType(const LALValue *value);
size_t XLALValueGetSize(const LALValue *value);
/* warning: shallow pointer */
const void * XLALValueGetDataPtr(const LALValue *value);
void * XLALValueGetData(void *data, size_t size, LALTYPECODE type, const LALValue *value);
int XLALValueEqual(const LALValue *value1, const LALValue *value2);

/* warning: shallow pointer */
const char * XLALValueGetString(const LALValue *value);
CHAR XLALValueGetCHAR(const LALValue *value);
INT2 XLALValueGetINT2(const LALValue *value);
INT4 XLALValueGetINT4(const LALValue *value);
INT8 XLALValueGetINT8(const LALValue *value);
UCHAR XLALValueGetUCHAR(const LALValue *value);
UINT2 XLALValueGetUINT2(const LALValue *value);
UINT4 XLALValueGetUINT4(const LALValue *value);
UINT8 XLALValueGetUINT8(const LALValue *value);
REAL4 XLALValueGetREAL4(const LALValue *value);
REAL8 XLALValueGetREAL8(const LALValue *value);
COMPLEX8 XLALValueGetCOMPLEX8(const LALValue *value);
COMPLEX16 XLALValueGetCOMPLEX16(const LALValue *value);

REAL8 XLALValueGetAsREAL8(const LALValue *value);

void XLALValuePrint(const LALValue *value, int fd);

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LAL_VALUE_H */
