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

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALValue.h>
#include "LALValue_private.h"

/* warning: garbage until set */
LALValue * XLALValueAlloc(size_t size)
{
	LALValue *v = XLALMalloc(sizeof(*v) + size);
	if (!v)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	v->size = size;
	return v;
}

/* warning: garbage until set */
LALValue * XLALValueRealloc(LALValue *value, size_t size)
{
	if (value == NULL)
		return XLALValueAlloc(size);
	value = XLALRealloc(value->data, sizeof(*value) + size);
	if (!value)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	value->size = size;
	return value;
}

LALValue * XLALValueDuplicate(const LALValue *value)
{
	size_t size = sizeof(LALValue) + value->size;
	LALValue *copy = XLALMalloc(size);
	if (!copy)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	return memcpy(copy, value, size);
}

LALValue * XLALValueCopy(LALValue *copy, const LALValue *orig)
{
	LALTYPECODE type = XLALValueGetType(orig);
	size_t size = XLALValueGetSize(orig);
	const void * data = XLALValueGetDataPtr(orig);
	return XLALValueSet(copy, data, size, type);
}

LALValue * XLALValueSet(LALValue *value, const void *data, size_t size, LALTYPECODE type)
{
	XLAL_CHECK_NULL(value != NULL, XLAL_EFAULT);

	/* sanity check type-size relation */
	switch (type) {
	case LAL_CHAR_TYPE_CODE:
		/* variable length strings allowed */
		break;
	case LAL_I2_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(INT2), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_I4_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(INT4), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_I8_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(INT8), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_UCHAR_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(UCHAR), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_U2_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(UINT2), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_U4_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(UINT4), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_U8_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(UINT8), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_S_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(REAL4), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_D_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(REAL8), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_C_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(COMPLEX8), XLAL_ETYPE, "Wrong size for type");
		break;
	case LAL_Z_TYPE_CODE:
		XLAL_CHECK_NULL(size == sizeof(COMPLEX16), XLAL_ETYPE, "Wrong size for type");
		break;
	default:
                XLAL_ERROR_NULL(XLAL_ETYPE, "Unsupported LALTYPECODE value 0%o", (unsigned int)type);
	}

	/* make sure sizes are compatible */
	XLAL_CHECK_NULL(size == value->size, XLAL_ESIZE, "Value has incompatible size");

	value->type = type;
	memcpy(value->data, data, size);
	return value;
}

void XLALDestroyValue(LALValue *value)
{
	XLALFree(value);
	return;
}

LALValue *XLALCreateValue(const void * data, size_t size, LALTYPECODE type)
{
	LALValue *v = XLALValueAlloc(size);
	if (v == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	v = XLALValueSet(v, data, size, type);
	if (v == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return v;
}

LALValue *XLALCreateStringValue(const char *value)
{
	size_t size = strlen(value) + 1;
	LALValue *v = XLALCreateValue(value, size, LAL_CHAR_TYPE_CODE);
	if (v == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return v;
}

#define DEFINE_CREATE_FUNC(TYPE, TCODE) \
	LALValue *XLALCreate ## TYPE ## Value(TYPE value) \
	{ \
		LALValue *v = XLALCreateValue(&value, sizeof(value), TCODE); \
		if (v == NULL) \
			XLAL_ERROR_NULL(XLAL_EFUNC); \
		return v; \
	}

DEFINE_CREATE_FUNC(CHAR, LAL_CHAR_TYPE_CODE)
DEFINE_CREATE_FUNC(INT2, LAL_I2_TYPE_CODE)
DEFINE_CREATE_FUNC(INT4, LAL_I4_TYPE_CODE)
DEFINE_CREATE_FUNC(INT8, LAL_I8_TYPE_CODE)
DEFINE_CREATE_FUNC(UCHAR, LAL_UCHAR_TYPE_CODE)
DEFINE_CREATE_FUNC(UINT2, LAL_U2_TYPE_CODE)
DEFINE_CREATE_FUNC(UINT4, LAL_U4_TYPE_CODE)
DEFINE_CREATE_FUNC(UINT8, LAL_U8_TYPE_CODE)
DEFINE_CREATE_FUNC(REAL4, LAL_S_TYPE_CODE)
DEFINE_CREATE_FUNC(REAL8, LAL_D_TYPE_CODE)
DEFINE_CREATE_FUNC(COMPLEX8, LAL_C_TYPE_CODE)
DEFINE_CREATE_FUNC(COMPLEX16, LAL_Z_TYPE_CODE)

#undef DEFINE_CREATE_FUNC

LALTYPECODE XLALValueGetType(const LALValue *value)
{
	return value->type;
}

size_t XLALValueGetSize(const LALValue *value)
{
	return value->size;
}

/* warning: shallow pointer */
const void * XLALValueGetDataPtr(const LALValue *value)
{
	return value->data;
}

void * XLALValueGetData(void *data, size_t size, LALTYPECODE type, const LALValue *value)
{
	if (value->size != size || value->type != type)
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect value type");
	return memcpy(data, value->data, size);
}

int XLALValueEqual(const LALValue *value1, const LALValue *value2)
{
	if (value1->size == value2->size && value1->type == value2->type)
		return memcmp(value1->data, value2->data, value1->size) == 0;
	return 0;
}

/* warning: shallow pointer */
const char * XLALValueGetString(const LALValue *value)
{
	/* sanity check the type */
	if (value->type != LAL_CHAR_TYPE_CODE)
		XLAL_ERROR_NULL(XLAL_ETYPE, "Value is not a string");
	/* make sure this is a nul-terminated string */
	if (value->size == 0 || ((const char *)(value->data))[value->size - 1] != '\0')
		XLAL_ERROR_NULL(XLAL_ETYPE, "Value is not a string");
	return (const char *)(value->data);
}

#define DEFINE_GET_FUNC(TYPE, TCODE, FAILVAL) \
	TYPE XLALValueGet ## TYPE (const LALValue *value) \
	{ \
		XLAL_CHECK_VAL(FAILVAL, value->type == TCODE, XLAL_ETYPE); \
		return *(const TYPE *)(value->data); \
	}

DEFINE_GET_FUNC(CHAR, LAL_CHAR_TYPE_CODE, XLAL_FAILURE)
DEFINE_GET_FUNC(INT2, LAL_I2_TYPE_CODE, XLAL_FAILURE)
DEFINE_GET_FUNC(INT4, LAL_I4_TYPE_CODE, XLAL_FAILURE)
DEFINE_GET_FUNC(INT8, LAL_I8_TYPE_CODE, XLAL_FAILURE)
DEFINE_GET_FUNC(UCHAR, LAL_UCHAR_TYPE_CODE, XLAL_FAILURE)
DEFINE_GET_FUNC(UINT2, LAL_U2_TYPE_CODE, XLAL_FAILURE)
DEFINE_GET_FUNC(UINT4, LAL_U4_TYPE_CODE, XLAL_FAILURE)
DEFINE_GET_FUNC(UINT8, LAL_U8_TYPE_CODE, XLAL_FAILURE)
DEFINE_GET_FUNC(REAL4, LAL_S_TYPE_CODE, XLAL_REAL4_FAIL_NAN)
DEFINE_GET_FUNC(REAL8, LAL_D_TYPE_CODE, XLAL_REAL8_FAIL_NAN)
DEFINE_GET_FUNC(COMPLEX8, LAL_C_TYPE_CODE, XLAL_REAL4_FAIL_NAN)
DEFINE_GET_FUNC(COMPLEX16, LAL_Z_TYPE_CODE, XLAL_REAL8_FAIL_NAN)

#undef DEFINE_GET_FUNC

REAL8 XLALValueGetAsREAL8(const LALValue *value)
{
	const INT8 max_as_double = LAL_INT8_C(9007199254740992);
	const UINT8 umax_as_double = LAL_INT8_C(9007199254740992);
	INT8 i;
	UINT8 u;
	REAL8 result;
	switch (value->type) {
	case LAL_CHAR_TYPE_CODE:
		if (value->size == 1)
			result = *(const CHAR *)value;
		else
			XLAL_ERROR_REAL8(XLAL_ETYPE, "Cannot convert string to float");
		break;
	case LAL_I2_TYPE_CODE:
		result = XLALValueGetINT2(value);
		break;
	case LAL_I4_TYPE_CODE:
		result = XLALValueGetINT4(value);
		break;
	case LAL_I8_TYPE_CODE:
		result = i = XLALValueGetINT8(value);
		if (i > max_as_double || -i > max_as_double)
			XLAL_PRINT_WARNING("Loss of precision");
		break;
	case LAL_UCHAR_TYPE_CODE:
		result = XLALValueGetUCHAR(value);
		break;
	case LAL_U2_TYPE_CODE:
		result = XLALValueGetUINT2(value);
		break;
	case LAL_U4_TYPE_CODE:
		result = XLALValueGetUINT4(value);
		break;
	case LAL_U8_TYPE_CODE:
		result = u = XLALValueGetUINT8(value);
		if (u > umax_as_double)
			XLAL_PRINT_WARNING("Loss of precision");
		break;
	case LAL_S_TYPE_CODE:
		result = XLALValueGetREAL4(value);
		break;
	case LAL_D_TYPE_CODE:
		result = XLALValueGetREAL8(value);
		break;
	case LAL_C_TYPE_CODE:
	case LAL_Z_TYPE_CODE:
		XLAL_ERROR_REAL8(XLAL_ETYPE, "Cannot convert complex to float");
	default:
                XLAL_ERROR_REAL8(XLAL_ETYPE, "Unsupported LALTYPECODE value 0%o", (unsigned int)value->type);
	}
	return result;
}

void XLALValuePrint(const LALValue *value, int fd)
{
	int fd2 = dup(fd);
	FILE *fp = fdopen(fd2, "w");
	COMPLEX8 c;
	COMPLEX16 z;
	switch (value->type) {
	case LAL_CHAR_TYPE_CODE:
		if (value->size == 1)
			fprintf(fp, "'%c'", *(const CHAR *)(value->data));
		else
			fprintf(fp, "\"%s\"", (const CHAR *)(value->data));
		break;
	case LAL_I2_TYPE_CODE:
		if (value->size != sizeof(INT2)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%" LAL_INT2_PRId, XLALValueGetINT2(value));
		break;
	case LAL_I4_TYPE_CODE:
		if (value->size != sizeof(INT4)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%" LAL_INT4_PRId, XLALValueGetINT4(value));
		break;
	case LAL_I8_TYPE_CODE:
		if (value->size != sizeof(INT8)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%" LAL_INT8_PRId, XLALValueGetINT8(value));
		break;
	case LAL_UCHAR_TYPE_CODE:
		if (value->size != sizeof(UCHAR)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "0x%x", XLALValueGetUCHAR(value));
		break;
	case LAL_U2_TYPE_CODE:
		if (value->size != sizeof(UINT2)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%" LAL_INT2_PRIu, XLALValueGetUINT2(value));
		break;
	case LAL_U4_TYPE_CODE:
		if (value->size != sizeof(UINT4)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%" LAL_INT4_PRIu, XLALValueGetUINT4(value));
		break;
	case LAL_U8_TYPE_CODE:
		if (value->size != sizeof(UINT8)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%" LAL_INT8_PRIu, XLALValueGetUINT8(value));
		break;
	case LAL_S_TYPE_CODE:
		if (value->size != sizeof(REAL4)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%g", XLALValueGetREAL4(value));
		break;
	case LAL_D_TYPE_CODE:
		if (value->size != sizeof(REAL8)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%g", XLALValueGetREAL8(value));
		break;
	case LAL_C_TYPE_CODE:
		c = XLALValueGetCOMPLEX8(value);
		if (value->size != sizeof(COMPLEX8)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%g + %gj", crealf(c), cimagf(c));
		break;
	case LAL_Z_TYPE_CODE:
		z = XLALValueGetCOMPLEX16(value);
		if (value->size != sizeof(COMPLEX16)) {
			fclose(fp);
			XLAL_ERROR_VOID(XLAL_ESIZE, "Value has incorrect size for type");
		}
		fprintf(fp, "%g + %gj", creal(z), cimag(z));
		break;
	default:
		XLAL_ERROR_VOID(XLAL_ETYPE, "Unsupported LALTYPECODE value 0%o", (unsigned int)value->type);
	}
	fclose(fp);
	return;
}
