#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/Units.h>
#include <lal/H5FileIO.h>

#define TYPECODE CHAR
#define TYPE CHAR
#include "H5FileIOVectorHL_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U2
#define TYPE UINT2
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U4
#define TYPE UINT4
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U8
#define TYPE UINT8
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#include "H5FileIOFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE D
#define TYPE REAL8
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#include "H5FileIOFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE C 
#define TYPE COMPLEX8
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#include "H5FileIOFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE Z 
#define TYPE COMPLEX16
#include "H5FileIOVectorHL_source.c"
#include "H5FileIOArrayHL_source.c"
#include "H5FileIOTimeSeries_source.c"
#include "H5FileIOFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE

/**
 * @addtogroup H5FileIOHighLevel_c
 * @brief High-level routines for reading/writing HDF5 files.
 * @details
 * These routines are high-level routines for accessing HDF files.
 * @{
 */

/**
 * @name Routines to Write Vectors to HDF5 files
 * @{
 */

/**
 * @fn int XLALH5FileWriteCHARVector(LALH5File *file, const char *name, CHARVector *vector)
 * @brief Writes a vector to a #LALH5File
 * @details
 * Creates a new HDF5 dataset with name @p name within a HDF5 file
 * associated with the #LALH5File @p file structure and allocates a
 * The data in the dataset is set to the data in the vector @p vector.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for writing.
 *
 * @param file Pointer to a #LALH5File structure in which to create the dataset.
 * @param name Pointer to a string with the name of the dataset to create.
 * @param vector Pointer to a vector structure containing the data.
 * @retval 0 Success.
 * @retval -1 Failure.
 */

/**
 * @fn int XLALH5FileWriteINT2Vector(LALH5File *file, const char *name, INT2Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/**
 * @fn int XLALH5FileWriteINT4Vector(LALH5File *file, const char *name, INT4Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/**
 * @fn int XLALH5FileWriteINT8Vector(LALH5File *file, const char *name, INT8Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/**
 * @fn int XLALH5FileWriteUINT2Vector(LALH5File *file, const char *name, UINT2Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/**
 * @fn int XLALH5FileWriteUINT4Vector(LALH5File *file, const char *name, UINT4Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/**
 * @fn int XLALH5FileWriteUINT8Vector(LALH5File *file, const char *name, UINT8Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/**
 * @fn int XLALH5FileWriteREAL4Vector(LALH5File *file, const char *name, REAL4Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/**
 * @fn int XLALH5FileWriteREAL8Vector(LALH5File *file, const char *name, REAL8Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/**
 * @fn int XLALH5FileWriteCOMPLEX8Vector(LALH5File *file, const char *name, COMPLEX8Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/**
 * @fn int XLALH5FileWriteCOMPLEX16Vector(LALH5File *file, const char *name, COMPLEX16Vector *vector)
 * @copydoc XLALH5FileWriteCHARVector()
 */

/** @} */

/**
 * @name Routines to Write Arrays to HDF5 Files
 * @{
 */

/**
 * @fn int XLALH5FileWriteINT2Array(LALH5File *file, const char *name, INT2Array *array)
 * @brief Writes an array to a #LALH5File
 * @details
 * Creates a new HDF5 dataset with name @p name within a HDF5 file
 * associated with the #LALH5File @p file structure and allocates a
 * The data in the dataset is set to the data in the array @p array.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for writing.
 *
 * @param file Pointer to a #LALH5File structure in which to create the dataset.
 * @param name Pointer to a string with the name of the dataset to create.
 * @param array Pointer to an array structure containing the data.
 * @retval 0 Success.
 * @retval -1 Failure.
 */

/**
 * @fn int XLALH5FileWriteINT4Array(LALH5File *file, const char *name, INT4Array *array)
 * @copydoc XLALH5FileWriteINT2Array()
 */

/**
 * @fn int XLALH5FileWriteINT8Array(LALH5File *file, const char *name, INT8Array *array)
 * @copydoc XLALH5FileWriteINT2Array()
 */

/**
 * @fn int XLALH5FileWriteUINT2Array(LALH5File *file, const char *name, UINT2Array *array)
 * @copydoc XLALH5FileWriteINT2Array()
 */

/**
 * @fn int XLALH5FileWriteUINT4Array(LALH5File *file, const char *name, UINT4Array *array)
 * @copydoc XLALH5FileWriteINT2Array()
 */

/**
 * @fn int XLALH5FileWriteUINT8Array(LALH5File *file, const char *name, UINT8Array *array)
 * @copydoc XLALH5FileWriteINT2Array()
 */

/**
 * @fn int XLALH5FileWriteREAL4Array(LALH5File *file, const char *name, REAL4Array *array)
 * @copydoc XLALH5FileWriteINT2Array()
 */

/**
 * @fn int XLALH5FileWriteREAL8Array(LALH5File *file, const char *name, REAL8Array *array)
 * @copydoc XLALH5FileWriteINT2Array()
 */

/**
 * @fn int XLALH5FileWriteCOMPLEX8Array(LALH5File *file, const char *name, COMPLEX8Array *array)
 * @copydoc XLALH5FileWriteINT2Array()
 */

/**
 * @fn int XLALH5FileWriteCOMPLEX16Array(LALH5File *file, const char *name, COMPLEX16Array *array)
 * @copydoc XLALH5FileWriteINT2Array()
 */

/** @} */

/**
 * @name Routines to Write Time Series to HDF5 Files
 * @{
 */

/**
 * @fn int XLALH5FileWriteINT2TimeSeries(LALH5File *file, const char *name, INT2TimeSeries *series)
 * @brief Writes a time series to a #LALH5File
 * @details
 * Creates a new HDF5 dataset with name @p name within a HDF5 file
 * associated with the #LALH5File @p file structure and allocates a
 * The data in the dataset is set to the data in the time series @p series.
 * Other metadata (name, epoch, deltaT, f0, and sampleUnits) are set as
 * attributes.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for writing.
 *
 * @param file Pointer to a #LALH5File structure in which to create the dataset.
 * @param name Pointer to a string with the name of the dataset to create.
 * @param series Pointer to time series structure containing the data.
 * @retval 0 Success.
 * @retval -1 Failure.
 */

/**
 * @fn int XLALH5FileWriteINT4TimeSeries(LALH5File *file, const char *name, INT4TimeSeries *series)
 * @copydoc XLALH5FileWriteINT2TimeSeries()
 */

/**
 * @fn int XLALH5FileWriteINT8TimeSeries(LALH5File *file, const char *name, INT8TimeSeries *series)
 * @copydoc XLALH5FileWriteINT2TimeSeries()
 */

/**
 * @fn int XLALH5FileWriteUINT2TimeSeries(LALH5File *file, const char *name, UINT2TimeSeries *series)
 * @copydoc XLALH5FileWriteINT2TimeSeries()
 */

/**
 * @fn int XLALH5FileWriteUINT4TimeSeries(LALH5File *file, const char *name, UINT4TimeSeries *series)
 * @copydoc XLALH5FileWriteINT2TimeSeries()
 */

/**
 * @fn int XLALH5FileWriteUINT8TimeSeries(LALH5File *file, const char *name, UINT8TimeSeries *series)
 * @copydoc XLALH5FileWriteINT2TimeSeries()
 */

/**
 * @fn int XLALH5FileWriteREAL4TimeSeries(LALH5File *file, const char *name, REAL4TimeSeries *series)
 * @copydoc XLALH5FileWriteINT2TimeSeries()
 */

/**
 * @fn int XLALH5FileWriteREAL8TimeSeries(LALH5File *file, const char *name, REAL8TimeSeries *series)
 * @copydoc XLALH5FileWriteINT2TimeSeries()
 */

/**
 * @fn int XLALH5FileWriteCOMPLEX8TimeSeries(LALH5File *file, const char *name, COMPLEX8TimeSeries *series)
 * @copydoc XLALH5FileWriteINT2TimeSeries()
 */

/**
 * @fn int XLALH5FileWriteCOMPLEX16TimeSeries(LALH5File *file, const char *name, COMPLEX16TimeSeries *series)
 * @copydoc XLALH5FileWriteINT2TimeSeries()
 */

/** @} */

/**
 * @name Routines to Write Frequency Series to HDF5 Files
 * @{
 */

/**
 * @fn int XLALH5FileWriteREAL4FrequencySeries(LALH5File *file, const char *name, REAL4FrequencySeries *series)
 * @brief Writes a frequency series to a #LALH5File
 * @details
 * Creates a new HDF5 dataset with name @p name within a HDF5 file
 * associated with the #LALH5File @p file structure and allocates a
 * The data in the dataset is set to the data in the frequency series @p series.
 * Other metadata (name, epoch, deltaF, f0, and sampleUnits) are set as
 * attributes.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for writing.
 *
 * @param file Pointer to a #LALH5File structure in which to create the dataset.
 * @param name Pointer to a string with the name of the dataset to create.
 * @param series Pointer to frequency series structure containing the data.
 * @retval 0 Success.
 * @retval -1 Failure.
 */

/**
 * @fn int XLALH5FileWriteREAL8FrequencySeries(LALH5File *file, const char *name, REAL8FrequencySeries *series)
 * @copydoc XLALH5FileWriteREAL4FrequencySeries()
 */

/**
 * @fn int XLALH5FileWriteCOMPLEX8FrequencySeries(LALH5File *file, const char *name, COMPLEX8FrequencySeries *series)
 * @copydoc XLALH5FileWriteREAL4FrequencySeries()
 */

/**
 * @fn int XLALH5FileWriteCOMPLEX16FrequencySeries(LALH5File *file, const char *name, COMPLEX16FrequencySeries *series)
 * @copydoc XLALH5FileWriteREAL4FrequencySeries()
 */

/** @} */

/**
 * @name Routines to Read Vectors from HDF5 Files
 * @{
 */

/**
 * @fn CHARVector *XLALH5FileReadCHARVector(LALH5File *file, const char *name)
 * @brief Reads a vector from a #LALH5File
 * @details
 * Reads a vector of data from a dataset named @p name in an HDF5 file
 * associated with the #LALH5File @p file.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for reading.
 * @param file Pointer to a #LALH5File to be read.
 * @param name Pointer to a string with the name of the dataset to read.
 * @returns Pointer to a vector containing the data in the dataset.
 * @retval NULL Failure.
 */

/**
 * @fn INT2Vector *XLALH5FileReadINT2Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/**
 * @fn INT4Vector *XLALH5FileReadINT4Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/**
 * @fn INT8Vector *XLALH5FileReadINT8Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/**
 * @fn UINT2Vector *XLALH5FileReadUINT2Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/**
 * @fn UINT4Vector *XLALH5FileReadUINT4Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/**
 * @fn UINT8Vector *XLALH5FileReadUINT8Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/**
 * @fn REAL4Vector *XLALH5FileReadREAL4Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/**
 * @fn REAL8Vector *XLALH5FileReadREAL8Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/**
 * @fn COMPLEX8Vector *XLALH5FileReadCOMPLEX8Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/**
 * @fn COMPLEX16Vector *XLALH5FileReadCOMPLEX16Vector(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadCHARVector()
 */

/** @} */

/**
 * @name Routines to Read Arrays from HDF5 Files
 * @{
 */

/**
 * @fn INT2Array *XLALH5FileReadINT2Array(LALH5File *file, const char *name)
 * @brief Reads an array from a #LALH5File
 * @details
 * Reads an array of data from a dataset named @p name in an HDF5 file
 * associated with the #LALH5File @p file.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for reading.
 * @param file Pointer to a #LALH5File to be read.
 * @param name Pointer to a string with the name of the dataset to read.
 * @returns Pointer to an array containing the data in the dataset.
 * @retval NULL Failure.
 */

/**
 * @fn INT4Array *XLALH5FileReadINT4Array(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2Array()
 */

/**
 * @fn INT8Array *XLALH5FileReadINT8Array(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2Array()
 */

/**
 * @fn UINT2Array *XLALH5FileReadUINT2Array(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2Array()
 */

/**
 * @fn UINT4Array *XLALH5FileReadUINT4Array(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2Array()
 */

/**
 * @fn UINT8Array *XLALH5FileReadUINT8Array(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2Array()
 */

/**
 * @fn REAL4Array *XLALH5FileReadREAL4Array(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2Array()
 */

/**
 * @fn REAL8Array *XLALH5FileReadREAL8Array(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2Array()
 */

/**
 * @fn COMPLEX8Array *XLALH5FileReadCOMPLEX8Array(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2Array()
 */

/**
 * @fn COMPLEX16Array *XLALH5FileReadCOMPLEX16Array(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2Array()
 */

/** @} */

/**
 * @name Routines to Read Time Series from HDF5 Files
 * @{
 */

/**
 * @fn INT2TimeSeries *XLALH5FileReadINT2TimeSeries(LALH5File *file, const char *name)
 * @brief Reads a time series from a #LALH5File
 * @details
 * Reads a time series of data from a dataset named @p name in an HDF5 file
 * associated with the #LALH5File @p file.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for reading.
 * @param file Pointer to a #LALH5File to be read.
 * @param name Pointer to a string with the name of the dataset to read.
 * @returns Pointer to a time series containing the data in the dataset.
 * @retval NULL Failure.
 */

/**
 * @fn INT4TimeSeries *XLALH5FileReadINT4TimeSeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2TimeSeries()
 */

/**
 * @fn INT8TimeSeries *XLALH5FileReadINT8TimeSeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2TimeSeries()
 */

/**
 * @fn UINT2TimeSeries *XLALH5FileReadUINT2TimeSeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2TimeSeries()
 */

/**
 * @fn UINT4TimeSeries *XLALH5FileReadUINT4TimeSeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2TimeSeries()
 */

/**
 * @fn UINT8TimeSeries *XLALH5FileReadUINT8TimeSeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2TimeSeries()
 */

/**
 * @fn REAL4TimeSeries *XLALH5FileReadREAL4TimeSeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2TimeSeries()
 */

/**
 * @fn REAL8TimeSeries *XLALH5FileReadREAL8TimeSeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2TimeSeries()
 */

/**
 * @fn COMPLEX8TimeSeries *XLALH5FileReadCOMPLEX8TimeSeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2TimeSeries()
 */

/**
 * @fn COMPLEX16TimeSeries *XLALH5FileReadCOMPLEX16TimeSeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadINT2TimeSeries()
 */

/** @} */

/**
 * @name Routines to Read Frequency Series from HDF5 Files
 * @{
 */

/**
 * @fn REAL4FrequencySeries *XLALH5FileReadREAL4FrequencySeries(LALH5File *file, const char *name)
 * @brief Reads a frequency series from a #LALH5File
 * @details
 * Reads a frequency series of data from a dataset named @p name in an HDF5 file
 * associated with the #LALH5File @p file.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for reading.
 * @param file Pointer to a #LALH5File to be read.
 * @param name Pointer to a string with the name of the dataset to read.
 * @returns Pointer to a frequency series containing the data in the dataset.
 * @retval NULL Failure.
 */

/**
 * @fn REAL8FrequencySeries *XLALH5FileReadREAL8FrequencySeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadREAL4FrequencySeries()
 */

/**
 * @fn COMPLEX8FrequencySeries *XLALH5FileReadCOMPLEX8FrequencySeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadREAL4FrequencySeries()
 */

/**
 * @fn COMPLEX16FrequencySeries *XLALH5FileReadCOMPLEX16FrequencySeries(LALH5File *file, const char *name)
 * @copydoc XLALH5FileReadREAL4FrequencySeries()
 */

/** @} */

/** @} */
