#include <config.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include <stdio.h>
#include <string.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/AVFactories.h>
#include <lal/H5FileIO.h>

/* INTERNAL */


#ifndef HAVE_HDF5
#pragma GCC diagnostic ignored "-Wunused-parameter"
#else

#define LAL_H5_FILE_MODE_READ  H5F_ACC_RDONLY
#define LAL_H5_FILE_MODE_WRITE H5F_ACC_TRUNC

struct tagLALH5File {
	char fname[FILENAME_MAX];
	hid_t file_id;
	unsigned int mode;
	int is_a_group;
};

struct tagLALH5Dataset {
	hid_t dataset_id;
	hid_t space_id;
	hid_t dtype_id; /* note: this is the in-memory type */
};

/* creates HDF5 float complex data type; use H5Tclose() to free */
typedef struct { float re; float im; } internal_float_complex_type;
static hid_t XLALH5TypeNativeFloatComplex(void)
{
	hid_t dtype_id;
	dtype_id = H5Tcreate(H5T_COMPOUND, sizeof(internal_float_complex_type));
	H5Tinsert(dtype_id, "r", HOFFSET(internal_float_complex_type, re), H5T_NATIVE_FLOAT);
	H5Tinsert(dtype_id, "i", HOFFSET(internal_float_complex_type, im), H5T_NATIVE_FLOAT);
	return dtype_id;
}

/* creates HDF5 double complex data type; use H5Tclose() to free */
typedef struct { double re; double im; } internal_double_complex_type;
static hid_t XLALH5TypeNativeDoubleComplex(void)
{
	hid_t dtype_id;
	dtype_id = H5Tcreate(H5T_COMPOUND, sizeof(internal_double_complex_type));
	H5Tinsert(dtype_id, "r", HOFFSET(internal_double_complex_type, re), H5T_NATIVE_DOUBLE);
	H5Tinsert(dtype_id, "i", HOFFSET(internal_double_complex_type, im), H5T_NATIVE_DOUBLE);
	return dtype_id;
}

/* creates HDF5 LIGOTimeGPS type; use H5Tclose() to free */
static hid_t XLALH5TypeNativeLIGOTimeGPS(void)
{
	hid_t dtype_id;
	hid_t int32_dtype_id;
	int32_dtype_id = H5Tcopy(H5T_NATIVE_INT);
	H5Tset_size(int32_dtype_id, 4); // 4 bytes = 32 bits
	dtype_id = H5Tcreate(H5T_COMPOUND, sizeof(LIGOTimeGPS));
	H5Tinsert(dtype_id, "gpsSeconds", HOFFSET(LIGOTimeGPS, gpsSeconds), int32_dtype_id);
	H5Tinsert(dtype_id, "gpsNanoSeconds", HOFFSET(LIGOTimeGPS, gpsNanoSeconds), int32_dtype_id);
	H5Tclose(int32_dtype_id);
	return dtype_id;
}

/* converts LALTYPECODE to HDF5 type; use H5Tclose() to free */
static hid_t XLALH5TypeFromLALType(LALTYPECODE dtype)
{
	hid_t dtype_id;
	size_t size = 1U << (dtype & LAL_TYPE_SIZE_MASK);

	/* allow only supported data types */
	switch (dtype) {
		case LAL_CHAR_TYPE_CODE:
		case LAL_I2_TYPE_CODE:
		case LAL_I4_TYPE_CODE:
		case LAL_I8_TYPE_CODE:
		case LAL_UCHAR_TYPE_CODE:
		case LAL_U2_TYPE_CODE:
		case LAL_U4_TYPE_CODE:
		case LAL_U8_TYPE_CODE:
		case LAL_S_TYPE_CODE:
		case LAL_D_TYPE_CODE:
		case LAL_C_TYPE_CODE:
		case LAL_Z_TYPE_CODE:
			/* these are all supported */
			break;
		default:
			/* anything else is not */
			XLAL_ERROR(XLAL_ETYPE, "Unsupported LALTYPECODE value 0%o", (unsigned int)dtype);
	}
	
	if ((dtype & LAL_CMPLX_TYPE_FLAG)) {
		/* complex numbers are treated as compound types */
		switch (size) {
		case 8:
			dtype_id = XLALH5TypeNativeFloatComplex();
			break;
		case 16:
			dtype_id = XLALH5TypeNativeDoubleComplex();
			break;
		default:
			/* it should be impossible to get here */
			XLAL_ERROR(XLAL_ETYPE, "Not reached");
		}
	} else if ((dtype & LAL_FLTPT_TYPE_FLAG)) {
		/* floating point number */
		switch (size) {
		case 4:
			dtype_id = H5Tcopy(H5T_NATIVE_FLOAT);
			break;
		case 8:
			dtype_id = H5Tcopy(H5T_NATIVE_DOUBLE);
			break;
		default:
			/* it should be impossible to get here */
			XLAL_ERROR(XLAL_ETYPE, "Not reached");
		}
	} else if ((dtype & LAL_UNSGN_TYPE_FLAG)) {
		/* unsigned integer */
		switch (size) {
		case 1:
			dtype_id = H5Tcopy(H5T_NATIVE_UCHAR);
			break;
		case 2:
			dtype_id = H5Tcopy(H5T_NATIVE_UINT16);
			break;
		case 4:
			dtype_id = H5Tcopy(H5T_NATIVE_UINT32);
			break;
		case 8:
			dtype_id = H5Tcopy(H5T_NATIVE_UINT64);
			break;
		default:
			/* it should be impossible to get here */
			XLAL_ERROR(XLAL_ETYPE, "Not reached");
		}
	} else {
		switch (size) {
		case 1:
			dtype_id = H5Tcopy(H5T_NATIVE_CHAR);
			break;
		case 2:
			dtype_id = H5Tcopy(H5T_NATIVE_INT16);
			break;
		case 4:
			dtype_id = H5Tcopy(H5T_NATIVE_INT32);
			break;
		case 8:
			dtype_id = H5Tcopy(H5T_NATIVE_INT64);
			break;
		default:
			/* it should be impossible to get here */
			XLAL_ERROR(XLAL_ETYPE, "Not reached");
		}
	}

	return dtype_id;
}

/* converts HDF5 type to LALTYPECODE */
static LALTYPECODE XLALTypeFromH5Type(hid_t dtype_id)
{
	LALTYPECODE dtype = 0;
	switch (H5Tget_class(dtype_id)) {
	case H5T_INTEGER:
		switch (H5Tget_sign(dtype_id)) {
		case H5T_SGN_NONE:
			dtype |= LAL_UNSGN_TYPE_FLAG;
			break;
		default:
			/* do nothing */
			break;
		}
		switch (H5Tget_size(dtype_id)) {
		case 1:
			dtype |= LAL_1_BYTE_TYPE_SIZE;
			break;
		case 2:
			dtype |= LAL_2_BYTE_TYPE_SIZE;
			break;
		case 4:
			dtype |= LAL_4_BYTE_TYPE_SIZE;
			break;
		case 8:
			dtype |= LAL_8_BYTE_TYPE_SIZE;
			break;
		default:
			XLAL_ERROR(XLAL_ETYPE, "Unsupported data type\n");
			break;
		}
		break;
	case H5T_FLOAT:
		dtype |= LAL_FLTPT_TYPE_FLAG;
		switch (H5Tget_size(dtype_id)) {
		case 4:
			dtype |= LAL_4_BYTE_TYPE_SIZE;
			break;
		case 8:
			dtype |= LAL_8_BYTE_TYPE_SIZE;
			break;
		default:
			XLAL_ERROR(XLAL_ETYPE, "Unsupported data type\n");
			break;
		}
		break;
	case H5T_COMPOUND: {
		/* note: complex numbers are the only compound type supported */
		char *s;

		/* sanity check the dtype_id */

		/* must have 2 members ... */
		if (H5Tget_nmembers(dtype_id) != 2)
			XLAL_ERROR(XLAL_ETYPE, "Unsupported data type\n");

		/* both members must be floating-point type ... */
		if (H5Tget_member_class(dtype_id, 0) != H5T_FLOAT)
			XLAL_ERROR(XLAL_ETYPE, "Unsupported data type\n");
		if (H5Tget_member_class(dtype_id, 1) != H5T_FLOAT)
			XLAL_ERROR(XLAL_ETYPE, "Unsupported data type\n");

		/* first member name must be something like "real" ... */
		s = H5Tget_member_name(dtype_id, 0);
		if (XLALStringNCaseCompare(s, "real", strlen(s)) != 0) {
			free(s);
			XLAL_ERROR(XLAL_ETYPE, "Unsupported data type\n");
		}
		free(s);

		/* second member name must be something like "imaginary" ... */
		s = H5Tget_member_name(dtype_id, 1);
		if (XLALStringNCaseCompare(s, "imaginary", strlen(s)) != 0) {
			free(s);
			XLAL_ERROR(XLAL_ETYPE, "Unsupported data type\n");
		}
		free(s);

		dtype |= LAL_CMPLX_TYPE_FLAG | LAL_FLTPT_TYPE_FLAG;
		switch (H5Tget_size(dtype_id)) {
		case 8:
			dtype |= LAL_8_BYTE_TYPE_SIZE;
			break;
		case 16:
			dtype |= LAL_16_BYTE_TYPE_SIZE;
			break;
		default:
			XLAL_ERROR(XLAL_ETYPE, "Unsupported data type\n");
			break;
		}
		break;
	}
	default:
		XLAL_ERROR(XLAL_ETYPE, "Unsupported data type\n");
		break;
	}
	return dtype;
}

/* creates a HDF5 file for writing */
static LALH5File * XLALH5FileCreate(const char *path)
{
	char tmpfname[FILENAME_MAX];
	LALH5File *file;
	if (snprintf(tmpfname, sizeof(tmpfname), "%s.tmp", path) < 0)
		XLAL_ERROR_NULL(XLAL_EINVAL);
	file = LALCalloc(1, sizeof(*file));
	if (!file)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	XLALStringCopy(file->fname, path, sizeof(file->fname));
	file->file_id = H5Fcreate(tmpfname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file->file_id < 0) {
		LALFree(file);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not create HDF5 file `%s'", path);
	}
	file->mode = LAL_H5_FILE_MODE_WRITE;
	return file;
}

/* opens a HDF5 file for reading */
static LALH5File * XLALH5FileOpenRead(const char *path)
{
	LALH5File *file;
	file = LALCalloc(1, sizeof(*file));
	if (!file)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	file->file_id = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file->file_id < 0) {
		LALFree(file);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not open HDF5 file `%s'", path);
	}
	file->mode = LAL_H5_FILE_MODE_READ;
	return file;
}


/* adds a scalar attribute */
static int XLALH5GenericAddScalarAttribute(hid_t hid, const char *key, const void *value, LALTYPECODE dtype)
{
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;

	dtype_id = XLALH5TypeFromLALType(dtype);
	if (dtype_id < 0)
		XLAL_ERROR(XLAL_EFUNC);

	space_id = H5Screate(H5S_SCALAR);
	if (space_id < 0) {
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create dataspace for attribute `%s'", key);
	}

	attr_id = H5Acreate(hid, key, dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(space_id);
	if (attr_id < 0) {
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create attribute `%s'", key);
	}

	if (H5Awrite(attr_id, dtype_id, value) < 0) {
		H5Aclose(attr_id);
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not write attribute `%s'", key);
	}
	
	H5Aclose(attr_id);
	H5Tclose(dtype_id);
	return 0;
}

/* adds a string attribute */
static int XLALH5GenericAddStringAttribute(hid_t hid, const char *key, const char *value)
{
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;

	dtype_id = H5Tcopy(H5T_C_S1);
	if (dtype_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not create datatype for attribure `%s'", key);
	if (H5Tset_size(dtype_id, H5T_VARIABLE) < 0)
	{
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create datatype for attribure `%s'", key);
	}

	space_id = H5Screate(H5S_SCALAR);
	if (space_id < 0) {
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create dataspace for attribute `%s'", key);
	}

	attr_id = H5Acreate(hid, key, dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(space_id);
	if (attr_id < 0) {
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create attribute `%s'", key);
	}

	if (H5Awrite(attr_id, dtype_id, &value) < 0) {
		H5Aclose(attr_id);
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not write attribute `%s'", key);
	}
	
	H5Aclose(attr_id);
	H5Tclose(dtype_id);
	return 0;
}

/* adds a LIGOTimeGPS attribute */
static int XLALH5GenericAddLIGOTimeGPSAttribute(hid_t hid, const char *key, const LIGOTimeGPS *value)
{
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;

	dtype_id = XLALH5TypeNativeLIGOTimeGPS();
	if (dtype_id < 0)
		XLAL_ERROR(XLAL_EFUNC);

	space_id = H5Screate(H5S_SCALAR);
	if (space_id < 0) {
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create dataspace for attribute `%s'", key);
	}

	attr_id = H5Acreate(hid, key, dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(space_id);
	if (attr_id < 0) {
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create attribute `%s'", key);
	}

	if (H5Awrite(attr_id, dtype_id, value) < 0) {
		H5Aclose(attr_id);
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not write attribute `%s'", key);
	}
	
	H5Aclose(attr_id);
	H5Tclose(dtype_id);
	return 0;
}

/* gets the datatype of an attribute */
static LALTYPECODE XLALH5GenericQueryScalarAttributeType(hid_t hid, const char *key)
{
	LALTYPECODE dtype;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;
	hid_t memtype_id;

	attr_id = H5Aopen(hid, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read dataset attribute `%s'", key);

	/* sanity check: make sure this is a scalar attribute */
	space_id = H5Aget_space(attr_id);
	if (space_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read dataspace of dataset attribute `%s'", key);
	}
	if (H5Sget_simple_extent_ndims(space_id) != 0 || H5Sget_simple_extent_npoints(space_id) != 1) {
		H5Sclose(space_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Attribute `%s' is not a scalar attribute", key);
	}

	dtype_id = H5Aget_type(attr_id);
	H5Aclose(attr_id);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read type of dataset attribute `%s'", key);

	memtype_id = H5Tget_native_type(dtype_id, H5T_DIR_ASCEND);
	H5Tclose(dtype_id);
	if (memtype_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not get native type of dataset attribute `%s'", key);

	dtype = XLALTypeFromH5Type(memtype_id);
	H5Tclose(memtype_id);
	if ((int)dtype < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return dtype;
}

/* gets the value of a scalar attribute */
static int XLALH5GenericQueryScalarAttributeValue(void *value, hid_t hid, const char *key)
{
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;
	hid_t memtype_id;

	attr_id = H5Aopen(hid, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read dataset attribute `%s'", key);

	/* sanity check: make sure this is a scalar attribute */
	space_id = H5Aget_space(attr_id);
	if (space_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read dataspace of dataset attribute `%s'", key);
	}
	if (H5Sget_simple_extent_ndims(space_id) != 0 || H5Sget_simple_extent_npoints(space_id) != 1) {
		H5Sclose(space_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Attribute `%s' is not a scalar attribute", key);
	}
	H5Sclose(space_id);

	dtype_id = H5Aget_type(attr_id);
	if (dtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type of dataset attribute `%s'", key);
	}

	memtype_id = H5Tget_native_type(dtype_id, H5T_DIR_ASCEND);
	H5Tclose(dtype_id);
	if (memtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not get native type of dataset attribute `%s'", key);
	}

	if (H5Aread(attr_id, memtype_id, value) < 0) {
		H5Tclose(memtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read data from dataset attribute `%s'", key);
	}

	H5Tclose(memtype_id);
	H5Aclose(attr_id);
	return 0;
}

/* gets the value of a string attribute */
static int XLALH5GenericQueryStringAttributeValue(char *value, size_t size, hid_t hid, const char *key)
{
	char *str;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;
	hid_t memtype_id;
	int n;

	attr_id = H5Aopen(hid, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read dataset attribute `%s'", key);

	/* sanity check: make sure this is just one variable-length string */
	space_id = H5Aget_space(attr_id);
	if (space_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read dataspace of dataset attribute `%s'", key);
	}
	if (H5Sget_simple_extent_ndims(space_id) != 0) {
		H5Sclose(space_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Attribute `%s' is not a string", key);
	}

	dtype_id = H5Aget_type(attr_id);
	if (dtype_id < 0) {
		H5Sclose(space_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type of dataset attribute `%s'", key);
	}
	H5Tclose(dtype_id);

	memtype_id = H5Tcopy(H5T_C_S1);
	if (memtype_id < 0 || H5Tset_size(memtype_id, H5T_VARIABLE)) {
		H5Sclose(space_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO);
	}

	if (H5Aread(attr_id, memtype_id, &str) < 0) {
		H5Tclose(memtype_id);
		H5Sclose(space_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read data from dataset attribute `%s'", key);
	}

	n = snprintf(value, value == NULL ? 0 : size, "%s", str);

	H5Dvlen_reclaim(memtype_id, space_id, H5P_DEFAULT, &str);

	H5Tclose(memtype_id);
	H5Sclose(space_id);
	H5Aclose(attr_id);
	return n;
}

/* gets the value of a LIGOTimeGPS attribute */
static LIGOTimeGPS * XLALH5GenericQueryLIGOTimeGPSAttributeValue(LIGOTimeGPS *value, hid_t hid, const char *key)
{
	char *s;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;
	hid_t memtype_id;

	attr_id = H5Aopen(hid, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read dataset attribute `%s'", key);

	/* sanity check: make sure this is a scalar attribute */
	space_id = H5Aget_space(attr_id);
	if (space_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read dataspace of dataset attribute `%s'", key);
	}
	if (H5Sget_simple_extent_ndims(space_id) != 0 || H5Sget_simple_extent_npoints(space_id) != 1) {
		H5Sclose(space_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_EIO, "Attribute `%s' is not a scalar attribute", key);
	}
	H5Sclose(space_id);

	/* sanity check the dtype_id */

	dtype_id = H5Aget_type(attr_id);
	if (dtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read type of dataset attribute `%s'", key);
	}

	/* must be compound data type ... */
	if (H5Tget_class(dtype_id) != H5T_COMPOUND) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for dataset attribute `%s'", key);
	}

	/* must have 2 members ... */
	if (H5Tget_nmembers(dtype_id) != 2) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for dataset attribute `%s'", key);
	}

	/* both members must be integer type ... */
	if (H5Tget_member_class(dtype_id, 0) != H5T_INTEGER) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for dataset attribute `%s'", key);
	}
	if (H5Tget_member_class(dtype_id, 1) != H5T_INTEGER) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for dataset attribute `%s'", key);
	}

	/* FIXME: should also check member sizes? */

	/* first member name must be "gpsSeconds" ... */
	s = H5Tget_member_name(dtype_id, 0);
	if (strcmp(s, "gpsSeconds") != 0) {
		free(s);
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for dataset attribute `%s'", key);
	}
	free(s);

	/* second member name must be "gpsNanoSeconds" ... */
	s = H5Tget_member_name(dtype_id, 1);
	if (strcmp(s, "gpsNanoSeconds") != 0) {
		free(s);
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for dataset attribute `%s'", key);
	}
	free(s);

	H5Tclose(dtype_id);

	memtype_id = XLALH5TypeNativeLIGOTimeGPS();
	if (memtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not get native type of dataset attribute `%s'", key);
	}

	if (H5Aread(attr_id, memtype_id, value) < 0) {
		H5Tclose(memtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read data from dataset attribute `%s'", key);
	}

	H5Tclose(memtype_id);
	H5Aclose(attr_id);
	return value;
}

#endif /* HAVE_HDF5 */

/* EXPORTED ROUTINES */

/**
 * @addtogroup H5FileIOLowLevel_c
 * @brief Low-level routines for reading/writing HDF5 files.
 * @details
 * These routines are basic routines for accessing HDF5 files.
 * The mid-level and high-level routines, which are based on
 * these routines, are typically more convenient routines to use.
 * @{
 */

/**
 * @name File Opening/Closing Routines
 * @{
 */

/**
 * @brief Closes a #LALH5File
 * @details
 * This routine closes a #LALH5File and deallocates resources
 * associated with it.  If the file was opened for writing, this
 * routine also renames the temporary file as the actual file.
 *
 * @param file A pointer to a #LALH5File structure to close.
 */
void XLALH5FileClose(LALH5File *file)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_VOID(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (file) {
		if (file->is_a_group)
			H5Gclose(file->file_id);
		else {
			if (file->mode == LAL_H5_FILE_MODE_WRITE) {
				char tmpfname[FILENAME_MAX];
				size_t namelen;
				namelen = H5Fget_name(file->file_id, NULL, 0);
				if (sizeof(tmpfname) <= namelen) {
					H5Fclose(file->file_id);
					LALFree(file);
					XLAL_ERROR_VOID(XLAL_EIO, "Failed to move temporary file");
				}
				H5Fget_name(file->file_id, tmpfname, sizeof(tmpfname));
				if (rename(tmpfname, file->fname) < 0) {
					H5Fclose(file->file_id);
					LALFree(file);
					XLAL_ERROR_VOID(XLAL_EIO, "Failed to move temporary file");
				}
			}
			H5Fclose(file->file_id);
		}
		LALFree(file);
	}
	return;
#endif
}

/**
 * @brief Opens a #LALH5File
 * @details
 * Opens a HDF5 file with pathname @p path and creates a #LALH5File structure
 * associated with it.
 *
 * The @p mode parameter points to a string that determines whether the
 * file is being opened for reading and writing.  Allowed strings are:
 *
 * <dl>
 * <dt>r</dt><dd>Open file for reading.</dd>
 * <dt>w</dt><dd>Truncate to zero length or create file for writing.</dd>
 * </dl>
 *
 * If a file is opened for writing then data is initially written to a
 * temporary file, and this file is renamed once the #LALH5File structure
 * is closed with XLALH5FileClose().
 *
 * @param path Pointer to a string containing the path of the file to open.
 * @param mode Mode to open the file, either "r" or "w".
 * @returns A pointer to a #LALH5File structure associated with the
 * specified HDF5 file.
 * @retval NULL An error occurred opening the file.
 */
LALH5File * XLALH5FileOpen(const char *path, const char *mode)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (path == NULL || mode == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (strcmp(mode, "r") == 0)
		return XLALH5FileOpenRead(path);
	else if (strcmp(mode, "w") == 0)
		return XLALH5FileCreate(path);
	XLAL_ERROR_NULL(XLAL_EINVAL, "Invalid mode \"%s\": must be either \"r\" or \"w\"", mode);
#endif
}

/**
 * @brief Opens a group in a #LALH5File
 * @details
 * Opens a HDF5 group with name @p name contained in the HDF5 file
 * associated with the #LALH5File @p file.  If the HDF5 file is
 * being read, the specified group must exist in that file.  If
 * the HDF5 file is being written, the specified group is created
 * within the file.
 *
 * @param file Pointer to a #LALH5File structure in which to open the group.
 * @param name Pointer to a string with the name of the group to open.
 * @returns A pointer to a #LALH5File structure associated with the
 * specified group within a HDF5 file.
 * @retval NULL An error occurred opening the group.
 */
LALH5File * XLALH5GroupOpen(LALH5File *file, const char *name)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALH5File *group;
	if (file == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	group = LALCalloc(1, sizeof(*group));
	if (!group)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	group->is_a_group = 1;
	group->mode = file->mode;
	if (!name) /* this is the same as the file */
		group->file_id = file->file_id;
	else if (group->mode == LAL_H5_FILE_MODE_READ)
		group->file_id = H5Gopen(file->file_id, name, H5P_DEFAULT);
	else if (group->mode == LAL_H5_FILE_MODE_WRITE) {
		hid_t gcpl; /* property list to allow intermediate groups to be created */
		gcpl = H5Pcreate(H5P_LINK_CREATE);
		if (gcpl < 0 || H5Pset_create_intermediate_group(gcpl, 1) < 0) {
			LALFree(group);
			XLAL_ERROR_NULL(XLAL_EIO);
		}
    		group->file_id = H5Gcreate(file->file_id, name, gcpl, H5P_DEFAULT, H5P_DEFAULT);
		H5Pclose(gcpl);
	} else 
		XLAL_ERROR_NULL(XLAL_EINVAL, "Corrupted file structure");
	if (group->file_id < 0) {
		LALFree(group);
		XLAL_ERROR_NULL(XLAL_EIO, "Failed to open group `%s'", name ? name : "(null)");
	}
	return group;
#endif
}

/** @} */

/**
 * @name File Attribute Routines
 * @{
 */

/**
 * @brief Adds a scalar attribute to a #LALH5File
 * @details
 * This routine adds a scalar-valued attribute with name @p key
 * and value given by the memory addressed by @p value to a
 * HDF5 file associated with the #LALH5File @p file.
 * The data type of the scalar value is given by the #LALTYPECODE
 * @p dtype.
 * @param file Pointer to a #LALH5File to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to the value of the scalar attribute to be added.
 * @param dtype #LALTYPECODE value specifying the data type of the attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileAddScalarAttribute(LALH5File *file, const char *key, const void *value, LALTYPECODE dtype)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (file == NULL || key == NULL || value == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (XLALH5GenericAddScalarAttribute(file->file_id, key, value, dtype) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief Adds a string attribute to a #LALH5File
 * @details
 * This routine adds a NUL-terminated variable-length string @p value
 * attribute with name @p key to a HDF5 file associated with the
 * #LALH5File @p file.
 * @param file Pointer to a #LALH5File to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a string with the value of the new attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileAddStringAttribute(LALH5File *file, const char *key, const char *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (file == NULL || key == NULL || value == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (XLALH5GenericAddStringAttribute(file->file_id, key, value) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief Adds a LIGOTimeGPS attribute to a #LALH5File
 * @details
 * This routine adds a LIGOTimeGPS @p value attribute with name @p key to a
 * HDF5 file associated with the #LALH5File @p file.
 * @param file Pointer to a #LALH5File to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a LIGOTimeGPS structure with the value of the new
 * attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileAddLIGOTimeGPSAttribute(LALH5File *file, const char *key, const LIGOTimeGPS *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (file == NULL || key == NULL || value == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (XLALH5GenericAddLIGOTimeGPSAttribute(file->file_id, key, value) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief Gets the datatype of an attribute in a #LALH5File
 * @details
 * This routine queries the datatype of a scalar attribute with
 * name @p key in a HDF5 file associated with the #LALH5File @p file.
 * @param file Pointer to a #LALH5File to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns #LALTYPECODE value of the datatype of the scalar attribute.
 * @retval -1 Failure.
 */
LALTYPECODE XLALH5FileQueryScalarAttributeType(LALH5File *file, const char *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALTYPECODE dtype;

	if (file == NULL || key == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	dtype = XLALH5GenericQueryScalarAttributeType(file->file_id, key);
	if ((int)(dtype) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return dtype;
#endif
}

/**
 * @brief Gets the value of a scalar attribute in a #LALH5File
 * @details
 * This routine queries the value of a scalar attribute with
 * name @p key in a HDF5 file associated with the #LALH5File @p file.
 * The value is stored in memory pointed to by the pointer @p value.
 *
 * @attention
 * This routine does not allocate memory for @p value.  The calling
 * routine must ensure that the memory addressed by the pointer @p value
 * is sufficient to hold the value in the attribute.
 *
 * @param value Pointer to memory in which the value will be stored.
 * @param file Pointer to a #LALH5File to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileQueryScalarAttributeValue(void *value, LALH5File *file, const char *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (value == NULL || file == NULL || key == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (XLALH5GenericQueryScalarAttributeValue(value, file->file_id, key) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief Gets the value of a string attribute in a #LALH5File
 * @details
 * This routine queries the value of a string attribute with
 * name @p key in a HDF5 file associated with the #LALH5File @p file.
 * The result is written into the buffer pointed to by @p value, the size
 * of which is @p size bytes.  If @p value is NULL, no data is copied but
 * the routine returns the length of the string.  Therefore, this routine
 * can be called once to determine the amount of memory required, the
 * memory can be allocated, and then it can be called a second time to
 * read the string.  If the parameter @p size is less than or equal to
 * the string length then only $p size-1 bytes of the string are copied
 * to the buffer @p value.
 * @note The return value is the length of the string, not including the
 * terminating NUL character; thus the buffer @p value should be allocated
 * to be one byte larger.
 * @param value Pointer to a buffer into which the string will be written.
 * @param size Size in bytes of the buffer into which the string will be
 * written.
 * @param file Pointer to a #LALH5File to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns The number of bytes that would be written to @p value had @p size
 * been sufficiently large excluding the terminating NUL byte.
 * @retval NULL Failure.
 */
int XLALH5FileQueryStringAttributeValue(char *value, size_t size, LALH5File *file, const char *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	int n;

	if (file == NULL || key == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	n = XLALH5GenericQueryStringAttributeValue(value, size, file->file_id, key);
	if (n < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return n;
#endif
}

/**
 * @brief Gets the value of a LIGOTimeGPS attribute in a #LALH5File
 * @details
 * This routine queries the value of a LIGOTimeGPS attribute with
 * name @p key in a HDF5 file associated with the #LALH5File @p file.
 * The value is stored in memory pointed to by the pointer @p value.
 * @param value Pointer to a LIGOTimeGPS structure in which the attribute
 * value will be stored.
 * @param file Pointer to a #LALH5File to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns Pointer to the LIGOTimeGPS structure passed to this routine.
 * @retval NULL Failure.
 */
LIGOTimeGPS * XLALH5FileQueryLIGOTimeGPSAttributeValue(LIGOTimeGPS *value, LALH5File *file, const char *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (value == NULL || file == NULL || key == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (XLALH5GenericQueryLIGOTimeGPSAttributeValue(value, file->file_id, key) == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return value;
#endif
}

/** @} */

/**
 * @name Dataset Routines
 * @{
 */

/**
 * @brief Frees a #LALH5Dataset
 * @details
 * Closes a HDF5 dataset associated with the #LALH5Dataset @p dset
 * and deallocates memory of the #LALH5Dataset structure.
 * @param dset Pointer to a #LALH5Dataset structure to close.
 */
void XLALH5DatasetFree(LALH5Dataset *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_VOID(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (dset) {
		H5Tclose(dset->dtype_id);
		H5Sclose(dset->space_id);
		H5Dclose(dset->dataset_id);
		LALFree(dset);
	}
	return;
#endif
}

/**
 * @brief Allocates a multi-dimensional #LALH5Dataset
 * @details
 * Creates a new HDF5 dataset with name @p name within a HDF5 file
 * associated with the #LALH5File @p file structure and allocates a
 * #LALH5Dataset structure associated with the dataset.  The type
 * of data to be stored in the dataset is given by the #LALTYPECODE
 * @p dtype and the rank array dimensions of the dataset is given by
 * the UINT4Vector @p dimLength.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for writing.
 *
 * @param file Pointer to a #LALH5File structure in which to create the dataset.
 * @param name Pointer to a string with the name of the dataset to create.
 * @param dtype #LALTYPECODE value specifying the data type.
 * @param dimLength Pointer to a UINT4Vector specifying the dataspace
 * dimensions.
 * @returns A pointer to a #LALH5Dataset structure associated with the
 * specified dataset within a HDF5 file.
 * @retval NULL An error occurred creating the dataset.
 */
LALH5Dataset * XLALH5DatasetAlloc(LALH5File *file, const char *name, LALTYPECODE dtype, UINT4Vector *dimLength)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALH5Dataset *dset;
	hsize_t *dims;
	UINT4 dim;

	if (name == NULL || file == NULL || dimLength == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (file->mode != LAL_H5_FILE_MODE_WRITE)
		XLAL_ERROR_NULL(XLAL_EINVAL, "Attempting to write to a read-only HDF5 file");

	dset = LALCalloc(1, sizeof(*dset));
	if (!dset)
		XLAL_ERROR_NULL(XLAL_ENOMEM);

	/* create datatype */
	dset->dtype_id = XLALH5TypeFromLALType(dtype);
	if (dset->dtype_id < 0) {
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* copy dimensions to HDF5 type */
	dims = LALCalloc(dimLength->length, sizeof(*dims));
	if (!dims) {
		H5Tclose(dset->dtype_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	}
	for (dim = 0; dim < dimLength->length; ++dim)
		dims[dim] = dimLength->data[dim];

	/* create dataspace */
	dset->space_id = H5Screate_simple(dimLength->length, dims, NULL);
	LALFree(dims);
	if (dset->space_id < 0) {
		H5Tclose(dset->dtype_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not create dataspace for dataset `%s'", name);
	}

	/* create dataset */
	dset->dataset_id = H5Dcreate(file->file_id, name, dset->dtype_id, dset->space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset->dataset_id < 0) {
		H5Tclose(dset->dtype_id);
		H5Sclose(dset->space_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not create dataset `%s'", name);
	}
	
	return dset;
#endif
}

/**
 * @brief Allocates a 1-dimensional #LALH5Dataset
 * @details
 * Creates a new HDF5 dataset with name @p name within a HDF5 file
 * associated with the #LALH5File @p file structure and allocates a
 * #LALH5Dataset structure associated with the dataset.  The type
 * of data to be stored in the dataset is given by the #LALTYPECODE
 * @p dtype and the number of points in the dataset is given by
 * the @p length parameter.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for writing.
 *
 * @param file Pointer to a #LALH5File structure in which to create the dataset.
 * @param name Pointer to a string with the name of the dataset to create.
 * @param dtype #LALTYPECODE value specifying the data type.
 * @param length The number of points of data in the dataset.
 * @returns A pointer to a #LALH5Dataset structure associated with the
 * specified dataset within a HDF5 file.
 * @retval NULL An error occurred creating the dataset.
 */
LALH5Dataset * XLALH5DatasetAlloc1D(LALH5File *file, const char *name, LALTYPECODE dtype, size_t length)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALH5Dataset *dset;
	hsize_t npoints = length;

	if (name == NULL || file == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (file->mode != LAL_H5_FILE_MODE_WRITE)
		XLAL_ERROR_NULL(XLAL_EINVAL, "Attempting to write to a read-only HDF5 file");

	dset = LALCalloc(1, sizeof(*dset));
	if (!dset)
		XLAL_ERROR_NULL(XLAL_ENOMEM);

	/* create datatype */
	dset->dtype_id = XLALH5TypeFromLALType(dtype);
	if (dset->dtype_id < 0) {
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* create dataspace */
	dset->space_id = H5Screate_simple(1, &npoints, NULL);
	if (dset->space_id < 0) {
		H5Tclose(dset->dtype_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not create dataspace for dataset `%s'", name);
	}

	/* create dataset */
	dset->dataset_id = H5Dcreate(file->file_id, name, dset->dtype_id, dset->space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset->dataset_id < 0) {
		H5Tclose(dset->dtype_id);
		H5Sclose(dset->space_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not create dataset `%s'", name);
	}
	
	return dset;
#endif
}

/**
 * @brief Writes data to a #LALH5Dataset
 * @details
 * Writes the data contained in @p data to a HDF5 dataset associated
 * with the #LALH5Dataset @p dset structure.
 * @param dset Pointer to a #LALH5Dataset structure to which to write the data.
 * @param data Pointer to the data buffer to be written.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetWrite(LALH5Dataset *dset, void *data)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (dset == NULL || data == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (H5Dwrite(dset->dataset_id, dset->dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
		XLAL_ERROR(XLAL_EIO, "Could not write data to dataset");
	return 0;
#endif
}

/**
 * @brief Reads a #LALH5Dataset
 * @details
 * Opens an existing HDF5 dataset with name @p name within a HDF5 file
 * associated with the #LALH5File @p file structure and allocates a
 * #LALH5Dataset structure associated with the dataset.
 *
 * The #LALH5File @p file passed to this routine must be a file
 * opened for reading.
 *
 * @param file Pointer to a #LALH5File structure containing the dataset
 * to be opened.
 * @param name Pointer to a string with the name of the dataset to open.
 * @returns A pointer to a #LALH5Dataset structure associated with the
 * specified dataset within a HDF5 file.
 * @retval NULL An error occurred creating the dataset.
 */
LALH5Dataset * XLALH5DatasetRead(LALH5File *file, const char *name)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t dtype_id;
	LALH5Dataset *dset;
	if (name == NULL || file == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (file->mode != LAL_H5_FILE_MODE_READ)
		XLAL_ERROR_NULL(XLAL_EINVAL, "Attempting to read a write-only HDF5 file");
	dset = LALCalloc(1, sizeof(*dset));
	if (!dset)
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	dset->dataset_id = H5Dopen(file->file_id, name, H5P_DEFAULT);
	if (dset->dataset_id < 0) {
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read dataset `%s'", name);
	}
	dset->space_id = H5Dget_space(dset->dataset_id);
	if (dset->space_id < 0) {
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read dataspace of dataset `%s'", name);
	}
	dtype_id = H5Dget_type(dset->dataset_id);
	if (dtype_id < 0) {
		H5Sclose(dset->space_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read datatype of dataset `%s'", name);
	}
	/* convert type to native type */
	dset->dtype_id = H5Tget_native_type(dtype_id, H5T_DIR_ASCEND);
	H5Tclose(dtype_id);
	if (dset->dtype_id < 0) {
		H5Sclose(dset->space_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not get native type for dataset `%s'", name);
	}
	return dset;
#endif
}

/**
 * @brief Gets the number of points in a #LALH5Dataset
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @returns The number of points in the HDF5 dataset associated
 * with the specified #LALH5Dataset.
 * @retval (size_t)(-1) Failure.
 */
size_t XLALH5DatasetQueryNPoints(LALH5Dataset *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hssize_t npoints;
	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	npoints = H5Sget_simple_extent_npoints(dset->space_id);
	if (npoints < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read number of points in dataspace");
	return npoints;
#endif
}

/**
 * @brief Gets the number of bytes in a #LALH5Dataset
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @returns The number of bytes in the HDF5 dataset associated
 * with the specified #LALH5Dataset.  This is the number of
 * bytes required to hold the data in that dataset.
 * @retval (size_t)(-1) Failure.
 */
size_t XLALH5DatasetQueryNBytes(LALH5Dataset *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	size_t size;
	size_t npoints;
	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	size = H5Tget_size(dset->dtype_id);
	if (size == 0)
		XLAL_ERROR(XLAL_EIO, "Could not read size of datatype");
	npoints = XLALH5DatasetQueryNPoints(dset);
	if (npoints == (size_t)(-1))
		XLAL_ERROR(XLAL_EFUNC);
	return size * npoints;
#endif
}

/**
 * @brief Gets the number of type of data in a #LALH5Dataset
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @returns The number of #LALTYPECODE of the datatype in the
 * HDF5 dataset associated with the specified #LALH5Dataset.
 * @retval -1 Failure.
 */
LALTYPECODE XLALH5DatasetQueryType(LALH5Dataset *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALTYPECODE dtype;
	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	dtype = XLALTypeFromH5Type(dset->dtype_id);
	if ((int)dtype < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return dtype;
#endif
}

/**
 * @brief Gets the number of dimensions of the dataspace in a #LALH5Dataset
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @returns The number of dimensions in the dataspace in the
 * HDF5 dataset associated with the specified #LALH5Dataset.
 * @retval -1 Failure.
 */
int XLALH5DatasetQueryNDim(LALH5Dataset *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	int rank;
	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	rank = H5Sget_simple_extent_ndims(dset->space_id);
	if (rank < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read rank of dataset");
	return rank;
#endif
}

/**
 * @brief Gets the dimensions of the dataspace in a #LALH5Dataset
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @returns A pointer to a newly-allocated UINT4Vector containing
 * the dimensions of the dataspace in the HDF5 dataset associated
 * with the * specified #LALH5Dataset.
 * @retval NULL Failure.
 */
UINT4Vector * XLALH5DatasetQueryDims(LALH5Dataset *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	UINT4Vector *dimLength;
	hsize_t *dims;
	int rank;
	int dim;

	if (dset == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);

	rank = XLALH5DatasetQueryNDim(dset);
	if (rank < 0)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	dims = LALCalloc(rank, sizeof(*dims));
	if (H5Sget_simple_extent_dims(dset->space_id, dims, NULL) < 0) {
		LALFree(dims);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read dimensions of dataspace");
	}

	dimLength = XLALCreateUINT4Vector(rank);
	if (dimLength == NULL) {
		LALFree(dims);
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	}

	for (dim = 0; dim < rank; ++dim)
		dimLength->data[dim] = dims[dim];
	LALFree(dims);

	return dimLength;
#endif
}

/**
 * @brief Gets the data contained in a #LALH5Dataset
 * @details
 * This routine reads data from a HDF5 dataset associated with
 * the #LALH5Dataset @p dset and stores the data in the buffer
 * @p data.  This buffer should be sufficiently large to hold
 * the entire contents of the dataset; this size can be determined
 * with the routine XLALH5DatasetQueryNBytes().
 * @param data Pointer to a memory in which to store the data.
 * @param dset Pointer to a #LALH5Dataset from which to extract the data.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetQueryData(void *data, LALH5Dataset *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (data == NULL || dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (H5Dread(dset->dataset_id, dset->dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
		LALFree(data);
		XLAL_ERROR(XLAL_EIO, "Could not read data from dataset");
	}
	return 0;
#endif
}

/** @} */

/**
 * @name Dataset Attribute Routines
 * @{
 */

/**
 * @brief Adds a scalar attribute to a #LALH5Dataset
 * @details
 * This routine adds a scalar-valued attribute with name @p key
 * and value given by the memory addressed by @p value to a
 * HDF5 dataset associated with the #LALH5Dataset @p dset.
 * The data type of the scalar value is given by the #LALTYPECODE
 * @p dtype.
 * @param dset Pointer to a #LALH5Dataset to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to the value of the scalar attribute to be added.
 * @param dtype #LALTYPECODE value specifying the data type of the attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetAddScalarAttribute(LALH5Dataset *dset, const char *key, const void *value, LALTYPECODE dtype)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (dset == NULL || key == NULL || value == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (XLALH5GenericAddScalarAttribute(dset->dataset_id, key, value, dtype) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief Adds a string attribute to a #LALH5Dataset
 * @details
 * This routine adds a NUL-terminated variable-length string @p value
 * attribute with name @p key to a HDF5 dataset associated with the
 * #LALH5Dataset @p dset.
 * @param dset Pointer to a #LALH5Dataset to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a string with the value of the new attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetAddStringAttribute(LALH5Dataset *dset, const char *key, const char *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (dset == NULL || key == NULL || value == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (XLALH5GenericAddStringAttribute(dset->dataset_id, key, value) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief Adds a LIGOTimeGPS attribute to a #LALH5Dataset
 * @details
 * This routine adds a LIGOTimeGPS @p value attribute with name @p key to a
 * HDF5 dataset associated with the #LALH5Dataset @p dset.
 * @param dset Pointer to a #LALH5Dataset to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a LIGOTimeGPS structure with the value of the new
 * attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetAddLIGOTimeGPSAttribute(LALH5Dataset *dset, const char *key, const LIGOTimeGPS *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (dset == NULL || key == NULL || value == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (XLALH5GenericAddLIGOTimeGPSAttribute(dset->dataset_id, key, value) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief Gets the datatype of an attribute in a #LALH5Dataset
 * @details
 * This routine queries the datatype of a scalar attribute with
 * name @p key in a HDF5 dataset associated with the #LALH5Dataset @p dset.
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns #LALTYPECODE value of the datatype of the scalar attribute.
 * @retval -1 Failure.
 */
LALTYPECODE XLALH5DatasetQueryScalarAttributeType(LALH5Dataset *dset, const char *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALTYPECODE dtype;

	if (dset == NULL || key == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	dtype = XLALH5GenericQueryScalarAttributeType(dset->dataset_id, key);
	if ((int)(dtype) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return dtype;
#endif
}

/**
 * @brief Gets the value of a scalar attribute in a #LALH5Dataset
 * @details
 * This routine queries the value of a scalar attribute with
 * name @p key in a HDF5 dataset associated with the #LALH5Dataset @p dset.
 * The value is stored in memory pointed to by the pointer @p value.
 *
 * @attention
 * This routine does not allocate memory for @p value.  The calling
 * routine must ensure that the memory addressed by the pointer @p value
 * is sufficient to hold the value in the attribute.
 *
 * @param value Pointer to memory in which the value will be stored.
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetQueryScalarAttributeValue(void *value, LALH5Dataset *dset, const char *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (value == NULL || dset == NULL || key == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (XLALH5GenericQueryScalarAttributeValue(value, dset->dataset_id, key) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief Gets the value of a string attribute in a #LALH5Dataset
 * @details
 * This routine queries the value of a string attribute with
 * name @p key in a HDF5 dataset associated with the #LALH5Dataset @p dset.
 * The result is written into the buffer pointed to by @p value, the size
 * of which is @p size bytes.  If @p value is NULL, no data is copied but
 * the routine returns the length of the string.  Therefore, this routine
 * can be called once to determine the amount of memory required, the
 * memory can be allocated, and then it can be called a second time to
 * read the string.  If the parameter @p size is less than or equal to
 * the string length then only $p size-1 bytes of the string are copied
 * to the buffer @p value.
 * @note The return value is the length of the string, not including the
 * terminating NUL character; thus the buffer @p value should be allocated
 * to be one byte larger.
 * @param value Pointer to a buffer into which the string will be written.
 * @param size Size in bytes of the buffer into which the string will be
 * written.
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns The number of bytes that would be written to @p value had @p size
 * been sufficiently large excluding the terminating NUL byte.
 * @retval NULL Failure.
 */
int XLALH5DatasetQueryStringAttributeValue(char *value, size_t size, LALH5Dataset *dset, const char *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	int n;

	if (dset == NULL || key == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	n = XLALH5GenericQueryStringAttributeValue(value, size, dset->dataset_id, key);
	if (n < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return n;
#endif
}

/**
 * @brief Gets the value of a LIGOTimeGPS attribute in a #LALH5Dataset
 * @details
 * This routine queries the value of a LIGOTimeGPS attribute with
 * name @p key in a HDF5 dataset associated with the #LALH5Dataset @p dset.
 * The value is stored in memory pointed to by the pointer @p value.
 * @param value Pointer to a LIGOTimeGPS structure in which the attribute
 * value will be stored.
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns Pointer to the LIGOTimeGPS structure passed to this routine.
 * @retval NULL Failure.
 */
LIGOTimeGPS * XLALH5DatasetQueryLIGOTimeGPSAttributeValue(LIGOTimeGPS *value, LALH5Dataset *dset, const char *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (value == NULL || dset == NULL || key == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (XLALH5GenericQueryLIGOTimeGPSAttributeValue(value, dset->dataset_id, key) == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return value;
#endif
}

/** @} */

/** @} */
