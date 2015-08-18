/*
*  Copyright (C) 2015 Jolien Creighton
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

#ifndef _H5FILEIO_H
#define _H5FILEIO_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include <stddef.h>
#include <lal/LALDatatypes.h>

/**
 * @defgroup H5FileIO_h Header H5FileIO.h
 * @ingroup lal_support
 *
 * @brief Provides standard LAL support for reading and
 * writing HDF5 files.
 *
 * ### EXAMPLE ###
 * 
 * @code
 * #include <lal/H5FileIO.h>
 * ...
 * LALH5File *file;
 * LALINT2Array *orig;
 * LALINT2Array *copy;
 * ...
 * // write an array to a file as dataset "testset"
 * // assume orig has been created and its data set
 * file = XLALH5FileOpen("testfile.h5", "w");
 * XLALH5FileWriteINT2Array(file, "testset", orig);
 * XLALH5FileClose(file);
 * ...
 * // read the array back from the file
 * file = XLALH5FileOpen("testfile.h5", "r");
 * copy = XLALH5FileReadINT2Array(file, "testset");
 * XLALH5FileClose(file);
 * @endcode
 *
 * @sa https://www.hdfgroup.org/HDF5/
 *
 * @{
 * @defgroup H5FileIOLowLevel_c Low-Level Routines
 * @defgroup H5FileIOMidLevel_c Mid-Level Routines
 * @defgroup H5FileIOHighLevel_c High-Level Routines
 * @}
 */

struct tagLALH5File;
/**
 * @brief Incomplete type for a HDF5 file.
 * @details
 * The #LALH5File is a structure that is associated with a HDF5 file or group.
 *
 * Allocate #LALH5File structures using XLALH5FileOpen() or
 * XLALH5FileGroupOpen().
 *
 * Deallocate #LALH5File structures using XLALH5FileClose().
 */
typedef struct tagLALH5File LALH5File;

struct tagLALH5Dataset;
/** 
 * @brief Incomplete type for a HDF5 dataset.
 * @details
 * The #LALH5Dataset is a structure that is associated with a HDF5 dataset.
 * Deallocate #LALH5Dataset structures using XLALH5DatasetFree().
 */
typedef struct tagLALH5Dataset LALH5Dataset;

void XLALH5FileClose(LALH5File *file);
LALH5File * XLALH5FileOpen(const char *path, const char *mode);
LALH5File * XLALH5GroupOpen(LALH5File *file, const char *name);

/* LOW-LEVEL ROUTINES */

void XLALH5DatasetFree(LALH5Dataset *dset);

LALH5Dataset * XLALH5DatasetAlloc(LALH5File *file, const char *name, LALTYPECODE dtype, UINT4Vector *dimLength);
LALH5Dataset * XLALH5DatasetAlloc1D(LALH5File *file, const char *name, LALTYPECODE dtype, size_t length);
int XLALH5DatasetWrite(LALH5Dataset *dset, void *data);

int XLALH5FileAddScalarAttribute(LALH5File *file, const char *key, const void *value, LALTYPECODE dtype);
int XLALH5FileAddStringAttribute(LALH5File *file, const char *key, const char *value);
int XLALH5FileAddLIGOTimeGPSAttribute(LALH5File *file, const char *key, const LIGOTimeGPS *value);

LALTYPECODE XLALH5FileQueryScalarAttributeType(LALH5File *file, const char *key);
int XLALH5FileQueryScalarAttributeValue(void *value, LALH5File *file, const char *key);
int XLALH5FileQueryStringAttributeValue(char *value, size_t size, LALH5File *file, const char *key);
LIGOTimeGPS * XLALH5FileQueryLIGOTimeGPSAttributeValue(LIGOTimeGPS *value, LALH5File *file, const char *key);

LALH5Dataset * XLALH5DatasetRead(LALH5File *file, const char *name);
size_t XLALH5DatasetQueryNPoints(LALH5Dataset *dset);
size_t XLALH5DatasetQueryNBytes(LALH5Dataset *dset);
LALTYPECODE XLALH5DatasetQueryType(LALH5Dataset *dset);
int XLALH5DatasetQueryNDim(LALH5Dataset *dset);
UINT4Vector * XLALH5DatasetQueryDims(LALH5Dataset *dset);
int XLALH5DatasetQueryData(void *data, LALH5Dataset *dset);

int XLALH5DatasetAddScalarAttribute(LALH5Dataset *dset, const char *key, const void *value, LALTYPECODE dtype);
int XLALH5DatasetAddStringAttribute(LALH5Dataset *dset, const char *key, const char *value);
int XLALH5DatasetAddLIGOTimeGPSAttribute(LALH5Dataset *dset, const char *key, const LIGOTimeGPS *value);

LALTYPECODE XLALH5DatasetQueryScalarAttributeType(LALH5Dataset *dset, const char *key);
int XLALH5DatasetQueryScalarAttributeValue(void *value, LALH5Dataset *dset, const char *key);
int XLALH5DatasetQueryStringAttributeValue(char *value, size_t size, LALH5Dataset *dset, const char *key);
LIGOTimeGPS * XLALH5DatasetQueryLIGOTimeGPSAttributeValue(LIGOTimeGPS *value, LALH5Dataset *dset, const char *key);

/* MID-LEVEL ROUTINES */

int XLALH5DatasetAddCHARAttribute(LALH5Dataset *dset, const char *key, CHAR value);
int XLALH5DatasetAddINT2Attribute(LALH5Dataset *dset, const char *key, INT2 value);
int XLALH5DatasetAddINT4Attribute(LALH5Dataset *dset, const char *key, INT4 value);
int XLALH5DatasetAddINT8Attribute(LALH5Dataset *dset, const char *key, INT8 value);
int XLALH5DatasetAddUCHARAttribute(LALH5Dataset *dset, const char *key, UCHAR value);
int XLALH5DatasetAddUINT2Attribute(LALH5Dataset *dset, const char *key, UINT2 value);
int XLALH5DatasetAddUINT4Attribute(LALH5Dataset *dset, const char *key, UINT4 value);
int XLALH5DatasetAddUINT8Attribute(LALH5Dataset *dset, const char *key, UINT8 value);
int XLALH5DatasetAddREAL4Attribute(LALH5Dataset *dset, const char *key, REAL4 value);
int XLALH5DatasetAddREAL8Attribute(LALH5Dataset *dset, const char *key, REAL8 value);
int XLALH5DatasetAddCOMPLEX8Attribute(LALH5Dataset *dset, const char *key, COMPLEX8 value);
int XLALH5DatasetAddCOMPLEX16Attribute(LALH5Dataset *dset, const char *key, COMPLEX16 value);

CHAR XLALH5DatasetQueryCHARAttributeValue(LALH5Dataset *dset, const char *key);
INT2 XLALH5DatasetQueryINT2AttributeValue(LALH5Dataset *dset, const char *key);
INT4 XLALH5DatasetQueryINT4AttributeValue(LALH5Dataset *dset, const char *key);
INT8 XLALH5DatasetQueryINT8AttributeValue(LALH5Dataset *dset, const char *key);
UCHAR XLALH5DatasetQueryUCHARAttributeValue(LALH5Dataset *dset, const char *key);
UINT2 XLALH5DatasetQueryUINT2AttributeValue(LALH5Dataset *dset, const char *key);
UINT4 XLALH5DatasetQueryUINT4AttributeValue(LALH5Dataset *dset, const char *key);
UINT8 XLALH5DatasetQueryUINT8AttributeValue(LALH5Dataset *dset, const char *key);
REAL4 XLALH5DatasetQueryREAL4AttributeValue(LALH5Dataset *dset, const char *key);
REAL8 XLALH5DatasetQueryREAL8AttributeValue(LALH5Dataset *dset, const char *key);
COMPLEX8 XLALH5DatasetQueryCOMPLEX8AttributeValue(LALH5Dataset *dset, const char *key);
COMPLEX16 XLALH5DatasetQueryCOMPLEX16AttributeValue(LALH5Dataset *dset, const char *key);

LALH5Dataset *XLALH5DatasetAllocCHARVector(LALH5File *file, const char *name, CHARVector *vector);
LALH5Dataset *XLALH5DatasetAllocINT2Vector(LALH5File *file, const char *name, INT2Vector *vector);
LALH5Dataset *XLALH5DatasetAllocINT4Vector(LALH5File *file, const char *name, INT4Vector *vector);
LALH5Dataset *XLALH5DatasetAllocINT8Vector(LALH5File *file, const char *name, INT8Vector *vector);
LALH5Dataset *XLALH5DatasetAllocUINT2Vector(LALH5File *file, const char *name, UINT2Vector *vector);
LALH5Dataset *XLALH5DatasetAllocUINT4Vector(LALH5File *file, const char *name, UINT4Vector *vector);
LALH5Dataset *XLALH5DatasetAllocUINT8Vector(LALH5File *file, const char *name, UINT8Vector *vector);
LALH5Dataset *XLALH5DatasetAllocREAL4Vector(LALH5File *file, const char *name, REAL4Vector *vector);
LALH5Dataset *XLALH5DatasetAllocREAL8Vector(LALH5File *file, const char *name, REAL8Vector *vector);
LALH5Dataset *XLALH5DatasetAllocCOMPLEX8Vector(LALH5File *file, const char *name, COMPLEX8Vector *vector);
LALH5Dataset *XLALH5DatasetAllocCOMPLEX16Vector(LALH5File *file, const char *name, COMPLEX16Vector *vector);

LALH5Dataset *XLALH5DatasetAllocINT2Array(LALH5File *file, const char *name, INT2Array *array);
LALH5Dataset *XLALH5DatasetAllocINT4Array(LALH5File *file, const char *name, INT4Array *array);
LALH5Dataset *XLALH5DatasetAllocINT8Array(LALH5File *file, const char *name, INT8Array *array);
LALH5Dataset *XLALH5DatasetAllocUINT2Array(LALH5File *file, const char *name, UINT2Array *array);
LALH5Dataset *XLALH5DatasetAllocUINT4Array(LALH5File *file, const char *name, UINT4Array *array);
LALH5Dataset *XLALH5DatasetAllocUINT8Array(LALH5File *file, const char *name, UINT8Array *array);
LALH5Dataset *XLALH5DatasetAllocREAL4Array(LALH5File *file, const char *name, REAL4Array *array);
LALH5Dataset *XLALH5DatasetAllocREAL8Array(LALH5File *file, const char *name, REAL8Array *array);
LALH5Dataset *XLALH5DatasetAllocCOMPLEX8Array(LALH5File *file, const char *name, COMPLEX8Array *array);
LALH5Dataset *XLALH5DatasetAllocCOMPLEX16Array(LALH5File *file, const char *name, COMPLEX16Array *array);

CHARVector *XLALH5DatasetReadCHARVector(LALH5Dataset *dset);
INT2Vector *XLALH5DatasetReadINT2Vector(LALH5Dataset *dset);
INT4Vector *XLALH5DatasetReadINT4Vector(LALH5Dataset *dset);
INT8Vector *XLALH5DatasetReadINT8Vector(LALH5Dataset *dset);
UINT2Vector *XLALH5DatasetReadUINT2Vector(LALH5Dataset *dset);
UINT4Vector *XLALH5DatasetReadUINT4Vector(LALH5Dataset *dset);
UINT8Vector *XLALH5DatasetReadUINT8Vector(LALH5Dataset *dset);
REAL4Vector *XLALH5DatasetReadREAL4Vector(LALH5Dataset *dset);
REAL8Vector *XLALH5DatasetReadREAL8Vector(LALH5Dataset *dset);
COMPLEX8Vector *XLALH5DatasetReadCOMPLEX8Vector(LALH5Dataset *dset);
COMPLEX16Vector *XLALH5DatasetReadCOMPLEX16Vector(LALH5Dataset *dset);

INT2Array *XLALH5DatasetReadINT2Array(LALH5Dataset *dset);
INT4Array *XLALH5DatasetReadINT4Array(LALH5Dataset *dset);
INT8Array *XLALH5DatasetReadINT8Array(LALH5Dataset *dset);
UINT2Array *XLALH5DatasetReadUINT2Array(LALH5Dataset *dset);
UINT4Array *XLALH5DatasetReadUINT4Array(LALH5Dataset *dset);
UINT8Array *XLALH5DatasetReadUINT8Array(LALH5Dataset *dset);
REAL4Array *XLALH5DatasetReadREAL4Array(LALH5Dataset *dset);
REAL8Array *XLALH5DatasetReadREAL8Array(LALH5Dataset *dset);
COMPLEX8Array *XLALH5DatasetReadCOMPLEX8Array(LALH5Dataset *dset);
COMPLEX16Array *XLALH5DatasetReadCOMPLEX16Array(LALH5Dataset *dset);

/* HIGH-LEVEL ROUTINES */

int XLALH5FileWriteCHARVector(LALH5File *file, const char *name, CHARVector *vector);
int XLALH5FileWriteINT2Vector(LALH5File *file, const char *name, INT2Vector *vector);
int XLALH5FileWriteINT4Vector(LALH5File *file, const char *name, INT4Vector *vector);
int XLALH5FileWriteINT8Vector(LALH5File *file, const char *name, INT8Vector *vector);
int XLALH5FileWriteUINT2Vector(LALH5File *file, const char *name, UINT2Vector *vector);
int XLALH5FileWriteUINT4Vector(LALH5File *file, const char *name, UINT4Vector *vector);
int XLALH5FileWriteUINT8Vector(LALH5File *file, const char *name, UINT8Vector *vector);
int XLALH5FileWriteREAL4Vector(LALH5File *file, const char *name, REAL4Vector *vector);
int XLALH5FileWriteREAL8Vector(LALH5File *file, const char *name, REAL8Vector *vector);
int XLALH5FileWriteCOMPLEX8Vector(LALH5File *file, const char *name, COMPLEX8Vector *vector);
int XLALH5FileWriteCOMPLEX16Vector(LALH5File *file, const char *name, COMPLEX16Vector *vector);

int XLALH5FileWriteINT2Array(LALH5File *file, const char *name, INT2Array *array);
int XLALH5FileWriteINT4Array(LALH5File *file, const char *name, INT4Array *array);
int XLALH5FileWriteINT8Array(LALH5File *file, const char *name, INT8Array *array);
int XLALH5FileWriteUINT2Array(LALH5File *file, const char *name, UINT2Array *array);
int XLALH5FileWriteUINT4Array(LALH5File *file, const char *name, UINT4Array *array);
int XLALH5FileWriteUINT8Array(LALH5File *file, const char *name, UINT8Array *array);
int XLALH5FileWriteREAL4Array(LALH5File *file, const char *name, REAL4Array *array);
int XLALH5FileWriteREAL8Array(LALH5File *file, const char *name, REAL8Array *array);
int XLALH5FileWriteCOMPLEX8Array(LALH5File *file, const char *name, COMPLEX8Array *array);
int XLALH5FileWriteCOMPLEX16Array(LALH5File *file, const char *name, COMPLEX16Array *array);

int XLALH5FileWriteINT2TimeSeries(LALH5File *file, const char *name, INT2TimeSeries *series);
int XLALH5FileWriteINT4TimeSeries(LALH5File *file, const char *name, INT4TimeSeries *series);
int XLALH5FileWriteINT8TimeSeries(LALH5File *file, const char *name, INT8TimeSeries *series);
int XLALH5FileWriteUINT2TimeSeries(LALH5File *file, const char *name, UINT2TimeSeries *series);
int XLALH5FileWriteUINT4TimeSeries(LALH5File *file, const char *name, UINT4TimeSeries *series);
int XLALH5FileWriteUINT8TimeSeries(LALH5File *file, const char *name, UINT8TimeSeries *series);
int XLALH5FileWriteREAL4TimeSeries(LALH5File *file, const char *name, REAL4TimeSeries *series);
int XLALH5FileWriteREAL8TimeSeries(LALH5File *file, const char *name, REAL8TimeSeries *series);
int XLALH5FileWriteCOMPLEX8TimeSeries(LALH5File *file, const char *name, COMPLEX8TimeSeries *series);
int XLALH5FileWriteCOMPLEX16TimeSeries(LALH5File *file, const char *name, COMPLEX16TimeSeries *series);

int XLALH5FileWriteREAL4FrequencySeries(LALH5File *file, const char *name, REAL4FrequencySeries *series);
int XLALH5FileWriteREAL8FrequencySeries(LALH5File *file, const char *name, REAL8FrequencySeries *series);
int XLALH5FileWriteCOMPLEX8FrequencySeries(LALH5File *file, const char *name, COMPLEX8FrequencySeries *series);
int XLALH5FileWriteCOMPLEX16FrequencySeries(LALH5File *file, const char *name, COMPLEX16FrequencySeries *series);

CHARVector *XLALH5FileReadCHARVector(LALH5File *file, const char *name);
INT2Vector *XLALH5FileReadINT2Vector(LALH5File *file, const char *name);
INT4Vector *XLALH5FileReadINT4Vector(LALH5File *file, const char *name);
INT8Vector *XLALH5FileReadINT8Vector(LALH5File *file, const char *name);
UINT2Vector *XLALH5FileReadUINT2Vector(LALH5File *file, const char *name);
UINT4Vector *XLALH5FileReadUINT4Vector(LALH5File *file, const char *name);
UINT8Vector *XLALH5FileReadUINT8Vector(LALH5File *file, const char *name);
REAL4Vector *XLALH5FileReadREAL4Vector(LALH5File *file, const char *name);
REAL8Vector *XLALH5FileReadREAL8Vector(LALH5File *file, const char *name);
COMPLEX8Vector *XLALH5FileReadCOMPLEX8Vector(LALH5File *file, const char *name);
COMPLEX16Vector *XLALH5FileReadCOMPLEX16Vector(LALH5File *file, const char *name);

INT2Array *XLALH5FileReadINT2Array(LALH5File *file, const char *name);
INT4Array *XLALH5FileReadINT4Array(LALH5File *file, const char *name);
INT8Array *XLALH5FileReadINT8Array(LALH5File *file, const char *name);
UINT2Array *XLALH5FileReadUINT2Array(LALH5File *file, const char *name);
UINT4Array *XLALH5FileReadUINT4Array(LALH5File *file, const char *name);
UINT8Array *XLALH5FileReadUINT8Array(LALH5File *file, const char *name);
REAL4Array *XLALH5FileReadREAL4Array(LALH5File *file, const char *name);
REAL8Array *XLALH5FileReadREAL8Array(LALH5File *file, const char *name);
COMPLEX8Array *XLALH5FileReadCOMPLEX8Array(LALH5File *file, const char *name);
COMPLEX16Array *XLALH5FileReadCOMPLEX16Array(LALH5File *file, const char *name);

INT2TimeSeries *XLALH5FileReadINT2TimeSeries(LALH5File *file, const char *name);
INT4TimeSeries *XLALH5FileReadINT4TimeSeries(LALH5File *file, const char *name);
INT8TimeSeries *XLALH5FileReadINT8TimeSeries(LALH5File *file, const char *name);
UINT2TimeSeries *XLALH5FileReadUINT2TimeSeries(LALH5File *file, const char *name);
UINT4TimeSeries *XLALH5FileReadUINT4TimeSeries(LALH5File *file, const char *name);
UINT8TimeSeries *XLALH5FileReadUINT8TimeSeries(LALH5File *file, const char *name);
REAL4TimeSeries *XLALH5FileReadREAL4TimeSeries(LALH5File *file, const char *name);
REAL8TimeSeries *XLALH5FileReadREAL8TimeSeries(LALH5File *file, const char *name);
COMPLEX8TimeSeries *XLALH5FileReadCOMPLEX8TimeSeries(LALH5File *file, const char *name);
COMPLEX16TimeSeries *XLALH5FileReadCOMPLEX16TimeSeries(LALH5File *file, const char *name);

REAL4FrequencySeries *XLALH5FileReadREAL4FrequencySeries(LALH5File *file, const char *name);
REAL8FrequencySeries *XLALH5FileReadREAL8FrequencySeries(LALH5File *file, const char *name);
COMPLEX8FrequencySeries *XLALH5FileReadCOMPLEX8FrequencySeries(LALH5File *file, const char *name);
COMPLEX16FrequencySeries *XLALH5FileReadCOMPLEX16FrequencySeries(LALH5File *file, const char *name);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
