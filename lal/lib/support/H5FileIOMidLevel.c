#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/StringVector.h>
#include <lal/H5FileIO.h>

#define FAILVAL XLAL_FAILURE

#define TYPECODE CHAR
#define TYPE CHAR
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE UCHAR
#define TYPE UCHAR
#include "H5FileIOScalar_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U2
#define TYPE UINT2
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U4
#define TYPE UINT4
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U8
#define TYPE UINT8
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#undef FAILVAL

#define FAILVAL XLAL_REAL4_FAIL_NAN

#define TYPECODE S
#define TYPE REAL4
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#undef FAILVAL

#define FAILVAL XLAL_REAL4_FAIL_NAN

#define TYPECODE D
#define TYPE REAL8
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#undef FAILVAL

#define FAILVAL XLAL_REAL4_FAIL_NAN

#define TYPECODE C 
#define TYPE COMPLEX8
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#undef FAILVAL

#define FAILVAL XLAL_REAL8_FAIL_NAN

#define TYPECODE Z 
#define TYPE COMPLEX16
#include "H5FileIOScalar_source.c"
#include "H5FileIOVector_source.c"
#include "H5FileIOArray_source.c"
#undef TYPECODE
#undef TYPE

#undef FAILVAL

/**
 * @addtogroup H5FileIOMidLevel_c
 * @brief Mid-level routines for reading/writing HDF5 files.
 * @details
 * These routines are mid-level routines for accessing HDF files.
 * They make use of the #LALH5Dataset structures.  For most
 * purposes, it is more convenient to use the high-level routines
 * that directly access #LALH5File structures.
 * @{
 */

/**
 * @name Routines to Add Attributes
 * @{
 */

/**
 * @fn int XLALH5DatasetAddCHARAttribute(LALH5Dataset *dset, const char *key, CHAR value)
 * @brief Adds a scalar attribute to a #LALH5Dataset
 * @details
 * This routine adds a scalar-valued attribute with name @p key and value @p
 * value to a HDF5 dataset associated with the #LALH5Dataset @p dset.
 * @param dset Pointer to a #LALH5Dataset to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Value of the scalar attribute to be added.
 * @retval 0 Success.
 * @retval -1 Failure.
 */

/**
 * @fn int XLALH5DatasetAddINT2Attribute(LALH5Dataset *dset, const char *key, INT2 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddINT4Attribute(LALH5Dataset *dset, const char *key, INT4 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddINT8Attribute(LALH5Dataset *dset, const char *key, INT8 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddUCHARAttribute(LALH5Dataset *dset, const char *key, UCHAR value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddUINT2Attribute(LALH5Dataset *dset, const char *key, UINT2 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddUINT4Attribute(LALH5Dataset *dset, const char *key, UINT4 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddUINT8Attribute(LALH5Dataset *dset, const char *key, UINT8 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddREAL4Attribute(LALH5Dataset *dset, const char *key, REAL4 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddREAL8Attribute(LALH5Dataset *dset, const char *key, REAL8 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddCOMPLEX8Attribute(LALH5Dataset *dset, const char *key, COMPLEX8 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/**
 * @fn int XLALH5DatasetAddCOMPLEX16Attribute(LALH5Dataset *dset, const char *key, COMPLEX16 value)
 * @copydoc XLALH5DatasetAddCHARAttribute()
 */

/** @} */

/**
 * @name Routines to Query Attributes
 * @{
 */

/**
 * @fn XLALH5DatasetQueryCHARAttributeValue(LALH5Dataset *dset, const char *key)
 * @brief Gets the value of a scalar attribute in a #LALH5Dataset
 * @details
 * This routine queries the datatype of a scalar attribute with
 * name @p key in a HDF5 dataset associated with the #LALH5Dataset @p dset.
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns Value of the scalar attribute.
 */

/**
 * @fn XLALH5DatasetQueryINT2AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryINT4AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryINT8AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryUCHARAttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryUINT2AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryUINT4AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryUINT8AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryREAL4AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryREAL8AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryCOMPLEX8AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/**
 * @fn XLALH5DatasetQueryCOMPLEX16AttributeValue(LALH5Dataset *dset, const char *key)
 * @copydoc XLALH5DatasetQueryCHARAttributeValue()
 */

/** @} */

/**
 * @name Routines to Allocate Vector Datasets
 * @{
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocCHARVector(LALH5File *file, const char *name, CHARVector *vector)
 * @brief Allocates a #LALH5Dataset
 * @details
 * Creates a new HDF5 dataset with name @p name within a HDF5 file
 * associated with the #LALH5File @p file structure and allocates a
 * #LALH5Dataset structure associated with the dataset.  The data
 * in the dataset is set to the data in the vector @p vector.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for writing.
 *
 * @param file Pointer to a #LALH5File structure in which to create the dataset.
 * @param name Pointer to a string with the name of the dataset to create.
 * @param vector Pointer to a vector structure containing the data.
 * @returns A pointer to a #LALH5Dataset structure associated with the
 * specified dataset within a HDF5 file.
 * @retval NULL An error occurred creating the dataset.
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocINT2Vector(LALH5File *file, const char *name, INT2Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocINT4Vector(LALH5File *file, const char *name, INT4Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocINT8Vector(LALH5File *file, const char *name, INT8Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocUINT2Vector(LALH5File *file, const char *name, UINT2Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocUINT4Vector(LALH5File *file, const char *name, UINT4Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocUINT8Vector(LALH5File *file, const char *name, UINT8Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocREAL4Vector(LALH5File *file, const char *name, REAL4Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocREAL8Vector(LALH5File *file, const char *name, REAL8Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocCOMPLEX8Vector(LALH5File *file, const char *name, COMPLEX8Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocCOMPLEX16Vector(LALH5File *file, const char *name, COMPLEX16Vector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocStringVector(LALH5File *file, const char *name, LALStringVector *vector)
 * @copydoc XLALH5DatasetAllocCHARVector()
 */
LALH5Dataset *XLALH5DatasetAllocStringVector(LALH5File *file, const char *name, LALStringVector *vector)
{
	LALH5Dataset *dataset;
	if (!file || !name || !vector)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (!vector->length || !vector->data)
		XLAL_ERROR_NULL(XLAL_EINVAL);
		dataset = XLALH5DatasetAllocStringData(file, name, vector->length);
	if (!dataset)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	if (XLALH5DatasetWrite(dataset, vector->data) < 0) {
		XLALH5DatasetFree(dataset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
        return dataset;
}

/** @} */

/**
 * @name Routines to Allocate Array Datasets
 * @{
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocINT2Array(LALH5File *file, const char *name, INT2Array *array)
 * @brief Allocates a #LALH5Dataset
 * @details
 * Creates a new HDF5 dataset with name @p name within a HDF5 file
 * associated with the #LALH5File @p file structure and allocates a
 * #LALH5Dataset structure associated with the dataset.  The data
 * in the dataset is set to the data in the array @p array.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for writing.
 *
 * @param file Pointer to a #LALH5File structure in which to create the dataset.
 * @param name Pointer to a string with the name of the dataset to create.
 * @param array Pointer to an array structure containing the data.
 * @returns A pointer to a #LALH5Dataset structure associated with the
 * specified dataset within a HDF5 file.
 * @retval NULL An error occurred creating the dataset.
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocINT4Array(LALH5File *file, const char *name, INT4Array *array)
 * @copydoc XLALH5DatasetAllocINT2Array()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocINT8Array(LALH5File *file, const char *name, INT8Array *array)
 * @copydoc XLALH5DatasetAllocINT2Array()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocUINT2Array(LALH5File *file, const char *name, UINT2Array *array)
 * @copydoc XLALH5DatasetAllocINT2Array()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocUINT4Array(LALH5File *file, const char *name, UINT4Array *array)
 * @copydoc XLALH5DatasetAllocINT2Array()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocUINT8Array(LALH5File *file, const char *name, UINT8Array *array)
 * @copydoc XLALH5DatasetAllocINT2Array()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocREAL4Array(LALH5File *file, const char *name, REAL4Array *array)
 * @copydoc XLALH5DatasetAllocINT2Array()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocREAL8Array(LALH5File *file, const char *name, REAL8Array *array)
 * @copydoc XLALH5DatasetAllocINT2Array()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocCOMPLEX8Array(LALH5File *file, const char *name, COMPLEX8Array *array)
 * @copydoc XLALH5DatasetAllocINT2Array()
 */

/**
 * @fn LALH5Dataset *XLALH5DatasetAllocCOMPLEX16Array(LALH5File *file, const char *name, COMPLEX16Array *array)
 * @copydoc XLALH5DatasetAllocINT2Array()
 */

/** @} */

/**
 * @name Routines to Read Vector Datasets
 * @{
 */

/**
 * @fn CHARVector *XLALH5DatasetReadCHARVector(LALH5Dataset *dset)
 * @brief Reads a #LALH5Dataset
 * @param dset Pointer to a #LALH5Dataset to be read.
 * @returns Pointer to a vector containing the data in the dataset.
 * @retval NULL Failure.
 */

/**
 * @fn INT2Vector *XLALH5DatasetReadINT2Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/**
 * @fn INT4Vector *XLALH5DatasetReadINT4Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/**
 * @fn INT8Vector *XLALH5DatasetReadINT8Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/**
 * @fn UINT2Vector *XLALH5DatasetReadUINT2Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/**
 * @fn UINT4Vector *XLALH5DatasetReadUINT4Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/**
 * @fn UINT8Vector *XLALH5DatasetReadUINT8Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/**
 * @fn REAL4Vector *XLALH5DatasetReadREAL4Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/**
 * @fn REAL8Vector *XLALH5DatasetReadREAL8Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/**
 * @fn COMPLEX8Vector *XLALH5DatasetReadCOMPLEX8Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/**
 * @fn COMPLEX16Vector *XLALH5DatasetReadCOMPLEX16Vector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */

/* helper routine to read a sequence of fixed-length character vectors */
static CHARVectorSequence *XLALH5DatasetReadCHARVectorSequence(LALH5Dataset *dset)
{
    CHARVectorSequence *sequence;
    size_t npoints;
    size_t nbytes;
    size_t size;
    int ndim;

    XLAL_CHECK_NULL(dset, XLAL_EFAULT);

	if (!XLALH5DatasetCheckFixedLengthStringData(dset))
		XLAL_ERROR_NULL(XLAL_ETYPE);

    ndim = XLALH5DatasetQueryNDim(dset);
    XLAL_CHECK_NULL(ndim == 1, XLAL_EDIMS);
   
    npoints = XLALH5DatasetQueryNPoints(dset);
    XLAL_CHECK_NULL(npoints != (size_t)(-1), XLAL_EFUNC);

    nbytes = XLALH5DatasetQueryNBytes(dset);
    XLAL_CHECK_NULL(npoints != (size_t)(-1), XLAL_EFUNC);

    size = nbytes / npoints;
    XLAL_CHECK_NULL(nbytes == size * npoints, XLAL_EBADLEN, "Number of bytes not divisible by number of points");

    sequence = XLALCreateCHARVectorSequence(npoints, size);
    XLAL_CHECK_NULL(sequence, XLAL_EFUNC);

    if (XLALH5DatasetQueryData(sequence->data, dset) == -1) {
         XLALDestroyCHARVectorSequence(sequence);
         XLAL_ERROR_NULL(XLAL_EFUNC);
    }

    return sequence;
}


#if 0
    UINT4 length; /**< The number \a l of vectors. */
    UINT4 vectorLength; /**< The length \a n of each vector. */
    CHAR *data; /**< Pointer to the data array.  Element \a i of vector \a j is \c data[ \a jn + \a i \c ]. */
#endif

/**
 * @fn LALStringVector *XLALH5DatasetReadStringVector(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadCHARVector()
 */
LALStringVector *XLALH5DatasetReadStringVector(LALH5Dataset *dset)
{
	LALStringVector *vector;
	size_t npoints;
	int ndim;

	if (!dset)
		XLAL_ERROR_NULL(XLAL_EFAULT);

	if (!XLALH5DatasetCheckStringData(dset)) {
		/* it might be fixed-length string data... */
		CHARVectorSequence *sequence;
		sequence = XLALH5DatasetReadCHARVectorSequence(dset);
		if (!sequence)
			XLAL_ERROR_NULL(XLAL_EFUNC);
		/* translate CHARVectorSequence to StringVector */
		vector = NULL;
		for (size_t j = 0; j < sequence->length; ++j) {
			LALStringVector *new = XLALAppendString2Vector(vector, sequence->data + j * sequence->vectorLength);
			if (!new) {
         			XLALDestroyCHARVectorSequence(sequence);
				XLALDestroyStringVector(vector);
				XLAL_ERROR_NULL(XLAL_EFUNC);
			}
			vector = new;
		}
         	XLALDestroyCHARVectorSequence(sequence);
		return vector;
	}
 
	ndim = XLALH5DatasetQueryNDim(dset);
	if (ndim != 1)
		XLAL_ERROR_NULL(XLAL_EDIMS);
 
	npoints = XLALH5DatasetQueryNPoints(dset);
	if (npoints == (size_t)(-1))
		XLAL_ERROR_NULL(XLAL_EFUNC);

	vector = XLALCreateEmptyStringVector(npoints);
	if (!vector)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	if (XLALH5DatasetQueryData(vector->data, dset) == -1) {
		XLALDestroyStringVector(vector);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	return vector;
}

/** @} */

/**
 * @name Routines to Read Array Datasets
 * @{
 */

/**
 * @fn INT2Array *XLALH5DatasetReadINT2Array(LALH5Dataset *dset)
 * @brief Reads a #LALH5Dataset
 * @param dset Pointer to a #LALH5Dataset to be read.
 * @returns Pointer to an array containing the data in the dataset.
 * @retval NULL Failure.
 */

/**
 * @fn INT4Array *XLALH5DatasetReadINT4Array(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadINT2Array()
 */

/**
 * @fn INT8Array *XLALH5DatasetReadINT8Array(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadINT2Array()
 */

/**
 * @fn UINT2Array *XLALH5DatasetReadUINT2Array(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadINT2Array()
 */

/**
 * @fn UINT4Array *XLALH5DatasetReadUINT4Array(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadINT2Array()
 */

/**
 * @fn UINT8Array *XLALH5DatasetReadUINT8Array(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadINT2Array()
 */

/**
 * @fn REAL4Array *XLALH5DatasetReadREAL4Array(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadINT2Array()
 */

/**
 * @fn REAL8Array *XLALH5DatasetReadREAL8Array(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadINT2Array()
 */

/**
 * @fn COMPLEX8Array *XLALH5DatasetReadCOMPLEX8Array(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadINT2Array()
 */

/**
 * @fn COMPLEX16Array *XLALH5DatasetReadCOMPLEX16Array(LALH5Dataset *dset)
 * @copydoc XLALH5DatasetReadINT2Array()
 */

/** @} */

/** @} */
