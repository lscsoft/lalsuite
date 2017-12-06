#include <config.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/AVFactories.h>
#include <lal/H5FileIO.h>

/* INTERNAL */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#ifndef HAVE_HDF5
#pragma GCC diagnostic ignored "-Wunused-parameter"
#else

#define LAL_H5_FILE_MODE_READ  H5F_ACC_RDONLY
#define LAL_H5_FILE_MODE_WRITE H5F_ACC_TRUNC

struct tagLALH5Object {
	hid_t object_id; /* this object's id must be first */
};

struct tagLALH5File {
	hid_t file_id; /* this object's id must be first */
	unsigned int mode;
	int is_a_group;
	char fname[FILENAME_MAX];
};

struct tagLALH5Dataset {
	hid_t dataset_id; /* this object's id must be first */
	hid_t parent_id;
	hid_t space_id;
	hid_t dtype_id; /* note: this is the in-memory type */
	char name[]; /* flexible array member must be last */
};

/* creates HDF5 enum data type; use H5Tclose() to free */
static hid_t XLALH5TypeEnum(const char *names[], const int values[], size_t length)
{
	hid_t dtype_id;
	size_t i;
	dtype_id = H5Tenum_create(H5T_NATIVE_INT);
	if (dtype_id < 0)
		XLAL_ERROR(XLAL_EIO);
	for (i = 0; i < length; ++i) {
		herr_t status = H5Tenum_insert(dtype_id, names[i], values + i);
		if (status < 0) {
			H5Tclose(dtype_id);
			XLAL_ERROR(XLAL_EIO);
		}
	}
	return dtype_id;
}

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

#if 0
static hid_t XLALGetObjectIdentifier(const void *ptr)
{
	/*
	 * use type punning to get the identifier
	 * note that the identifier is always the
	 * first member of the LALH5 types
	 */
	union { const void *ptr; const hid_t *hid; } id = {ptr};
	if (ptr == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	return *id.hid;
}
#endif

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
 * @name File/Group Routines
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
void XLALH5FileClose(LALH5File UNUSED *file)
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
                H5Fflush(file->file_id , H5F_SCOPE_GLOBAL);
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
LALH5File * XLALH5FileOpen(const char UNUSED *path, const char UNUSED *mode)
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
LALH5File * XLALH5GroupOpen(LALH5File UNUSED *file, const char UNUSED *name)
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
		group->file_id = H5Gopen2(file->file_id, name, H5P_DEFAULT);
	else if (group->mode == LAL_H5_FILE_MODE_WRITE) {
		hid_t gcpl; /* property list to allow intermediate groups to be created */
		gcpl = H5Pcreate(H5P_LINK_CREATE);
		if (gcpl < 0 || H5Pset_create_intermediate_group(gcpl, 1) < 0) {
			LALFree(group);
			XLAL_ERROR_NULL(XLAL_EIO);
		}
    		group->file_id = H5Gcreate2(file->file_id, name, gcpl, H5P_DEFAULT, H5P_DEFAULT);
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

/**
 * @brief Checks for existence of a group in a #LALH5File
 * @details
 * Checks if group with name @p name exists in the HDF5 file associated with the
 * #LALH5File file @p file. If the group exists the return code is 1.
 * If the group does not exist a return code value of 0 is used.
 *
 * @attention Read failure results in a returned value of 0.
 *
 * @param file Pointer to a #LALH5File structure to check for group.
 * @param name Pointer to a string with the name of the group to check.
 * @retval 0 Group does not exist or failure.
 * @retval 1 Group exists.
 */
int XLALH5FileCheckGroupExists(const LALH5File UNUSED *file, const char UNUSED *name)
{
#ifndef HAVE_HDF5
        XLAL_ERROR_VAL(0, XLAL_EFAILED, "HDF5 support not implemented");
#else
	H5G_info_t group_info;
	hsize_t i;

	if (file == NULL || name == NULL)
		XLAL_ERROR_VAL(0, XLAL_EFAULT);

	if (H5Gget_info(file->file_id, &group_info) < 0)
		XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read group info");

	for (i = 0; i < group_info.nlinks; ++i) {
		H5O_info_t obj_info;
		if (H5Oget_info_by_idx(file->file_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &obj_info, H5P_DEFAULT) < 0)
			XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read object info");
		if (obj_info.type == H5O_TYPE_GROUP) {
			char *base;
			hid_t obj_id;
			int n;
			obj_id = H5Oopen_by_addr(file->file_id, obj_info.addr);
			if (obj_id < 0)
				XLAL_ERROR_VAL(0, XLAL_EIO, "Could not open object");
			n = H5Iget_name(obj_id, NULL, 0);
			if (n < 0) {
				H5Oclose(obj_id);
				XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read object");
			}
			char group_name[n + 1];
			n = H5Iget_name(obj_id, group_name, n + 1);
			if (n < 0) {
				H5Oclose(obj_id);
				XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read object");
			}
			/* get basename */
			if ((base = strrchr(group_name, '/')))
				++base;
			else
				base = group_name;

			if (strcmp(name, base) == 0)
				return 1; /* found it */
		}
	}

	/* no matching group found */
	return 0;
#endif
}

/**
 * @brief Checks for existence of a dataset in a #LALH5File
 * @details
 * Checks if dataset with name @p name exists in the HDF5 file associated with the
 * #LALH5File file @p file. If the dataset exists the return code is 1.
 * If the dataset does not exist a return code value of 0 is used.
 *
 * @attention Read failure results in a returned value of 0.
 *
 * @param file Pointer to a #LALH5File structure to check for dataset.
 * @param name Pointer to a string with the name of the dataset to check.
 * @retval 0 Dataset does not exist or failure.
 * @retval 1 Dataset exists.
 */
int XLALH5FileCheckDatasetExists(const LALH5File UNUSED *file, const char UNUSED *name)
{
#ifndef HAVE_HDF5
        XLAL_ERROR_VAL(0, XLAL_EFAILED, "HDF5 support not implemented");
#else
	H5G_info_t group_info;
	hsize_t i;

	if (file == NULL || name == NULL)
		XLAL_ERROR_VAL(0, XLAL_EFAULT);

	if (H5Gget_info(file->file_id, &group_info) < 0)
		XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read group info");

	for (i = 0; i < group_info.nlinks; ++i) {
		H5O_info_t obj_info;
		if (H5Oget_info_by_idx(file->file_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &obj_info, H5P_DEFAULT) < 0)
			XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read object info");
		if (obj_info.type == H5O_TYPE_DATASET) {
			char *base;
			hid_t obj_id;
			int n;
			obj_id = H5Oopen_by_addr(file->file_id, obj_info.addr);
			if (obj_id < 0)
				XLAL_ERROR_VAL(0, XLAL_EIO, "Could not open object");
			n = H5Iget_name(obj_id, NULL, 0);
			if (n < 0) {
				H5Oclose(obj_id);
				XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read object");
			}
			char dset_name[n + 1];
			n = H5Iget_name(obj_id, dset_name, n + 1);
			H5Oclose(obj_id);
			if (n < 0)
				XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read object");
			/* get basename */
			if ((base = strrchr(dset_name, '/')))
				++base;
			else
				base = dset_name;

			if (strcmp(name, base) == 0)
				return 1; /* found it */
		}
	}

	/* no matching dataset found */
	return 0;
#endif
}

/**
 * @brief DEPRECATED: Checks for existence of a Group in a LALH5File object #LALH5File
 * @details
 * Checks if group with name @p exists in the HDF5 files associated with the
 * #LALH5File file. If the group exists the return code is 1
 * if the group does not exist a return code value of 0 is used.
 *
 * @deprecated
 * Use XLALH5FileCheckGroupExists instead.
 *
 * @param file Pointer to a #LALH5File structure to check for group in
 * @param name Pointer to a string with the name of the group to check.
 * @returns int, 1 if group exists. 0 if not.
 */
int XLALH5CheckGroupExists(LALH5File UNUSED *file, const char UNUSED *name)
{
#ifndef HAVE_HDF5
        XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	XLAL_PRINT_DEPRECATION_WARNING("XLALH5FileCheckGroupExists");
        if (H5Gget_objinfo(file->file_id, name, 0, NULL))
        {
          return 0;
        }
        else
        {
          return 1;
        }
#endif
}

/**
 * @brief Gets the number of groups contained in a #LALH5File
 * @details
 * This routines returns the number of groups contained in a
 * an #LALH5File @p file which can be either a file or a group.
 * This routine does not recursively count subgroups of the
 * groups found.
 *
 * @param file Pointer to a #LALH5File file or group to be queried.
 * @returns The number of groups contained in the file or group.
 * @retval -1 Failure.
 */
size_t XLALH5FileQueryNGroups(const LALH5File UNUSED *file)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	H5G_info_t group_info;
	size_t num = 0;
	hsize_t i;

	if (file == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5Gget_info(file->file_id, &group_info) < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read group info");

	for (i = 0; i < group_info.nlinks; ++i) {
		H5O_info_t obj_info;
		if (H5Oget_info_by_idx(file->file_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &obj_info, H5P_DEFAULT) < 0)
			XLAL_ERROR(XLAL_EIO, "Could not read object info");
		if (obj_info.type == H5O_TYPE_GROUP)
			++num;
	}

	return num;
#endif
}

/**
 * @brief Gets the name of a group contained in a #LALH5File
 * @details
 * This routines gets the name of a group contained in a #LALH5File
 * @p file which can be either a file or a group.
 * The index @p pos identifies which group's name is returned.
 * The result is written into the buffer pointed to by @p name, the size
 * of which is @p size bytes.  If @p name is NULL, no data is copied but
 * the routine returns the length of the string.  Therefore, this routine
 * can be called once to determine the amount of memory required, the
 * memory can be allocated, and then it can be called a second time to
 * read the string.  If the parameter @p size is less than or equal to
 * the string length then only $p size-1 bytes of the string are copied
 * to the buffer @p name.
 * @note The return value is the length of the string, not including the
 * terminating NUL character; thus the buffer @p name should be allocated
 * to be one byte larger.
 * @param name Pointer to a buffer into which the string will be written.
 * @param size Size in bytes of the buffer into which the string will be
 * written.
 * @param file Pointer to a #LALH5File file or group to be queried.
 * @param pos The index identifying which group contained in the file.
 * @retval  0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileQueryGroupName(char UNUSED *name, size_t UNUSED size, const LALH5File UNUSED *file, int UNUSED pos)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	H5G_info_t group_info;
	size_t num = 0;
	hsize_t i;

	if (file == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5Gget_info(file->file_id, &group_info) < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read group info");

	for (i = 0; i < group_info.nlinks; ++i) {
		H5O_info_t obj_info;
		if (H5Oget_info_by_idx(file->file_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &obj_info, H5P_DEFAULT) < 0)
			XLAL_ERROR(XLAL_EIO, "Could not read object info");
		if (obj_info.type == H5O_TYPE_GROUP) {
			if (num == (size_t)pos) { /* found it */
				hid_t obj_id;
				int n;
				obj_id = H5Oopen_by_addr(file->file_id, obj_info.addr);
				if (obj_id < 0)
					XLAL_ERROR(XLAL_EIO, "Could not open object");
				n = H5Iget_name(obj_id, name, size);
				H5Oclose(obj_id);
				if (n < 0)
					XLAL_ERROR(XLAL_EIO, "Could not read object name");
				return n;
			} else
				++num;
		}
	}

	/* failed to find position */
	XLAL_ERROR(XLAL_EINVAL, "No group associated with given position");
#endif
}

/**
 * @brief Gets the number of datasets contained in a #LALH5File
 * @details
 * This routines returns the number of datasets contained in a
 * an #LALH5File @p file which can be either a file or a group.
 *
 * @param file Pointer to a #LALH5File file or group to be queried.
 * @returns The number of datasets contained in the file or group.
 * @retval -1 Failure.
 */
size_t XLALH5FileQueryNDatasets(const LALH5File UNUSED *file)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	H5G_info_t group_info;
	size_t num = 0;
	hsize_t i;

	if (file == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5Gget_info(file->file_id, &group_info) < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read group info");

	for (i = 0; i < group_info.nlinks; ++i) {
		H5O_info_t obj_info;
		if (H5Oget_info_by_idx(file->file_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &obj_info, H5P_DEFAULT) < 0)
			XLAL_ERROR(XLAL_EIO, "Could not read object info");
		if (obj_info.type == H5O_TYPE_DATASET)
			++num;
	}

	return num;
#endif
}

/**
 * @brief Gets the name of a dataset contained in a #LALH5File
 * @details
 * This routines gets the name of a dataset contained in a #LALH5File
 * @p file which can be either a file or a group.
 * The index @p pos identifies which dataset's name is returned.
 * The result is written into the buffer pointed to by @p name, the size
 * of which is @p size bytes.  If @p name is NULL, no data is copied but
 * the routine returns the length of the string.  Therefore, this routine
 * can be called once to determine the amount of memory required, the
 * memory can be allocated, and then it can be called a second time to
 * read the string.  If the parameter @p size is less than or equal to
 * the string length then only $p size-1 bytes of the string are copied
 * to the buffer @p name.
 * @note The return value is the length of the string, not including the
 * terminating NUL character; thus the buffer @p name should be allocated
 * to be one byte larger.
 * @param name Pointer to a buffer into which the string will be written.
 * @param size Size in bytes of the buffer into which the string will be
 * written.
 * @param file Pointer to a #LALH5File file or group to be queried.
 * @param pos The index identifying which dataset contained in the file.
 * @retval  0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileQueryDatasetName(char UNUSED *name, size_t UNUSED size, const LALH5File UNUSED *file, int UNUSED pos)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	H5G_info_t group_info;
	size_t num = 0;
	hsize_t i;

	if (file == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5Gget_info(file->file_id, &group_info) < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read group info");

	for (i = 0; i < group_info.nlinks; ++i) {
		H5O_info_t obj_info;
		if (H5Oget_info_by_idx(file->file_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &obj_info, H5P_DEFAULT) < 0)
			XLAL_ERROR(XLAL_EIO, "Could not read object info");
		if (obj_info.type == H5O_TYPE_DATASET) {
			if (num == (size_t)pos) { /* found it */
				hid_t obj_id;
				int n;
				obj_id = H5Oopen_by_addr(file->file_id, obj_info.addr);
				if (obj_id < 0)
					XLAL_ERROR(XLAL_EIO, "Could not open object");
				n = H5Iget_name(obj_id, name, size);
				H5Oclose(obj_id);
				if (n < 0)
					XLAL_ERROR(XLAL_EIO, "Could not read object name");
				return n;
			} else
				++num;
		}
	}

	/* failed to find position */
	XLAL_ERROR(XLAL_EINVAL, "No dataset associated with given position");
#endif
}

/**
 * @brief DEPRECATED: Gets dataset names from a #LALH5File
 * @details
 * This routine returns the names of all datasets in a #LALH5File.
 *
 * @deprecated
 * Use XLALH5AttributeQueryN() and XLALH5AttributeQueryName() instead.
 *
 * @param names Pointer a list of strings to be returned to the user. Memory
 * should be freed by the caller.
 * @param N Pointer to a UINT4 where the number of datasets will be recorded
 * @param file #LALH5File from which to read datasets
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileGetDatasetNames(LALH5File UNUSED *file, char UNUSED *** names, UINT4 UNUSED *N)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
    hsize_t ng;
    int ns;
    int otype;
    int i;

	XLAL_PRINT_DEPRECATION_WARNING("XLALH5FileQueryNDatasets() and XLALH5FileQueryDatasetName");

	if (file == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (names == NULL)
		XLAL_ERROR(XLAL_EFAULT);

    H5Gget_num_objs(file->file_id, &ng);
    *N = ng;

    /*
     *  Filter objects that don't meet our requirements
     */
    for (i = 0; i < (int)ng; i++) {
        otype = H5Gget_objtype_by_idx(file->file_id, (size_t)i);
        if (otype != H5G_DATASET) {
            (*N)--;
        }
    }
    char ** namelist = (char**) XLALMalloc((*N) * sizeof(*names));

    for (i = 0; i < (int)ng; i++) {
        otype = H5Gget_objtype_by_idx(file->file_id, (size_t)i);
        if (otype != H5G_DATASET) {
            continue;
        }

        ns = H5Gget_objname_by_idx(file->file_id, (size_t)i, NULL, 0) + 1;
        namelist[i] = (char*) XLALMalloc(ns * sizeof(namelist[i]));
        H5Gget_objname_by_idx(file->file_id, (size_t)i, namelist[i], ns);
    }

  *names = namelist;
  return(XLAL_SUCCESS);
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
void XLALH5DatasetFree(LALH5Dataset UNUSED *dset)
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
LALH5Dataset * XLALH5DatasetAlloc(LALH5File UNUSED *file, const char UNUSED *name, LALTYPECODE UNUSED dtype, UINT4Vector UNUSED *dimLength)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALH5Dataset *dset;
	hsize_t *dims;
	UINT4 dim;
	size_t namelen;

	if (name == NULL || file == NULL || dimLength == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (file->mode != LAL_H5_FILE_MODE_WRITE)
		XLAL_ERROR_NULL(XLAL_EINVAL, "Attempting to write to a read-only HDF5 file");

	namelen = strlen(name);
	dset = LALCalloc(1, sizeof(*dset) + namelen + 1);  /* use flexible array member to record name */
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
	dset->dataset_id = H5Dcreate2(file->file_id, name, dset->dtype_id, dset->space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset->dataset_id < 0) {
		H5Tclose(dset->dtype_id);
		H5Sclose(dset->space_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not create dataset `%s'", name);
	}

	/* record name of dataset and parent id */
	snprintf(dset->name, namelen + 1, "%s", name);
	dset->parent_id = file->file_id;

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
LALH5Dataset * XLALH5DatasetAlloc1D(LALH5File UNUSED *file, const char UNUSED *name, LALTYPECODE UNUSED dtype, size_t UNUSED length)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALH5Dataset *dset;
	hsize_t npoints = length;
	size_t namelen;

	if (name == NULL || file == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (file->mode != LAL_H5_FILE_MODE_WRITE)
		XLAL_ERROR_NULL(XLAL_EINVAL, "Attempting to write to a read-only HDF5 file");

	namelen = strlen(name);
	dset = LALCalloc(1, sizeof(*dset) + namelen + 1);  /* use flexible array member to record name */
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
	dset->dataset_id = H5Dcreate2(file->file_id, name, dset->dtype_id, dset->space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset->dataset_id < 0) {
		H5Tclose(dset->dtype_id);
		H5Sclose(dset->space_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not create dataset `%s'", name);
	}
	
	/* record name of dataset and parent id */
	snprintf(dset->name, namelen + 1, "%s", name);
	dset->parent_id = file->file_id;

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
int XLALH5DatasetWrite(LALH5Dataset UNUSED *dset, void UNUSED *data)
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
LALH5Dataset * XLALH5DatasetRead(LALH5File UNUSED *file, const char UNUSED *name)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t dtype_id;
	LALH5Dataset *dset;
	size_t namelen;
	if (name == NULL || file == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (file->mode != LAL_H5_FILE_MODE_READ)
		XLAL_ERROR_NULL(XLAL_EINVAL, "Attempting to read a write-only HDF5 file");

	namelen = strlen(name);
	dset = LALCalloc(1, sizeof(*dset) + namelen + 1);  /* use flexible array member to record name */
	if (!dset)
		XLAL_ERROR_NULL(XLAL_ENOMEM);

	dset->dataset_id = H5Dopen2(file->file_id, name, H5P_DEFAULT);
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

	/* record name of dataset and parent id */
	snprintf(dset->name, namelen + 1, "%s", name);
	dset->parent_id = file->file_id;
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
size_t XLALH5DatasetQueryNPoints(LALH5Dataset UNUSED *dset)
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
size_t XLALH5DatasetQueryNBytes(LALH5Dataset UNUSED *dset)
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
 * @brief Gets the type of data in a #LALH5Dataset
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @returns The #LALTYPECODE of the datatype in the
 * HDF5 dataset associated with the specified #LALH5Dataset.
 * @retval -1 Failure.
 */
LALTYPECODE XLALH5DatasetQueryType(LALH5Dataset UNUSED *dset)
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
int XLALH5DatasetQueryNDim(LALH5Dataset UNUSED *dset)
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
UINT4Vector * XLALH5DatasetQueryDims(LALH5Dataset UNUSED *dset)
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
int XLALH5DatasetQueryData(void UNUSED *data, LALH5Dataset UNUSED *dset)
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
 * @name Attribute Routines
 * @anchor attribute_routines
 * @details
 * These routines allow for reading or writing attributes to #LALH5File
 * or #LALH5Dataset objects.  To use these, the pointer to the object
 * should be cast to a #LALH5Generic type as in the following example.
 * @code
 * LALH5File *file = XLALH5FileOpen("example.h5", "r");
 * LALH5Dataset *dset = XLALH5DatasetRead(file, "example_dataset");
 * size_t num_file_attrs;
 * size_t num_dset_attrs;
 * num_file_attrs = XLALH5AttributeQueryN((LALH5Generic)file);
 * num_dset_attrs = XLALH5AttributeQueryN((LALH5Generic)dset);
 * @endcode
 * @{
 */


/**
 * @brief Checks for existence of an attribute associated with a #LALH5File
 * or #LALH5Dataset
 * @details
 * Checks if there is an attribute with name @p name associated with
 * an HDF5 object @p object that is either a #LALH5File or
 * #LALH5Dataset.
 * If the attribute exists the return code is 1;
 * if the dataset does not exist a return code value of 0 is used.
 *
 * @attention Read failure results in a returned value of 0.
 *
 * @param object Pointer to a #LALH5File or #LALH5Dataset to check for
 * attribute.
 * @param name Pointer to a string with the name of the attribute to check.
 * @retval 0 Attribute does not exist or failure.
 * @retval 1 Attribute exists.
 */
size_t XLALH5AttributeCheckExists(const LALH5Generic UNUSED object, const char UNUSED *name)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	H5O_info_t obj_info;
	hsize_t i;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR_VAL(0, XLAL_EINVAL);

	if (H5Oget_info(obj_id, &obj_info) < 0)
		XLAL_ERROR_VAL(0, XLAL_EIO, "Could not get HDF5 object info");

	for (i = 0; i < obj_info.num_attrs; ++i) {
		hid_t attr_id;
		int n;

		attr_id = H5Aopen_idx(obj_id, i);
		if (attr_id < 0)
			XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read attribute");

		n = H5Aget_name(attr_id, 0, NULL);
		if (n < 0) {
			H5Aclose(attr_id);
			XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read attribute name");
		}
		char attr_name[n + 1];
		n = H5Aget_name(attr_id, n + 1, attr_name);
		H5Aclose(attr_id);
		if (n < 0)
			XLAL_ERROR_VAL(0, XLAL_EIO, "Could not read attribute name");
		if (strcmp(name, attr_name) == 0)
			return 1; /* found it */
	}

	/* no matching attribute found */
	return 0;
#endif
}

/**
 * @brief Gets the number of attributes associated with a #LALH5File
 * or #LALH5Dataset
 * @details
 * This routines returns the number of attributes associated with
 * an HDF5 object @p object that is either a #LALH5File or
 * #LALH5Dataset.
 *
 * @param object Pointer to a #LALH5File or #LALH5Dataset that will be
 * queried.
 * @returns The number of attributes associated with the object.
 * @retval -1 Failure.
 */
size_t XLALH5AttributeQueryN(const LALH5Generic UNUSED object)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	H5O_info_t obj_info;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	if (H5Oget_info(obj_id, &obj_info) < 0)
		XLAL_ERROR(XLAL_EIO, "Could not get HDF5 object info");

	return obj_info.num_attrs;
#endif
}

/**
 * @brief Gets the name of an attribute associated with a #LALH5File
 * or #LALH5Dataset
 * @details
 * This routines gets the name of an attributes associated with
 * an HDF5 object @p object that is either a #LALH5File or
 * #LALH5Dataset.
 * The index @p pos identifies which attribute's name is returned.
 * The result is written into the buffer pointed to by @p name, the size
 * of which is @p size bytes.  If @p name is NULL, no data is copied but
 * the routine returns the length of the string.  Therefore, this routine
 * can be called once to determine the amount of memory required, the
 * memory can be allocated, and then it can be called a second time to
 * read the string.  If the parameter @p size is less than or equal to
 * the string length then only $p size-1 bytes of the string are copied
 * to the buffer @p name.
 * @note The return value is the length of the string, not including the
 * terminating NUL character; thus the buffer @p name should be allocated
 * to be one byte larger.
 * @param name Pointer to a buffer into which the string will be written.
 * @param size Size in bytes of the buffer into which the string will be
 * written.
 * @param object Pointer to a #LALH5File or #LALH5Dataset that will be
 * queried.
 * @param pos The index identifying which attribute associated with the object.
 * @retval  0 Success.
 * @retval -1 Failure.
 */
int XLALH5AttributeQueryName(char UNUSED *name, size_t UNUSED size, const LALH5Generic UNUSED object, int UNUSED pos)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	hid_t attr_id;
	H5O_info_t obj_info;
	int n;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	if (H5Oget_info(obj_id, &obj_info) < 0)
		XLAL_ERROR(XLAL_EIO, "Could not get HDF5 object info");

	if ((hsize_t)pos >= obj_info.num_attrs)
		XLAL_ERROR(XLAL_EINVAL, "No attribute associated with given position");

	attr_id = H5Aopen_idx(obj_id, pos);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute");

	n = H5Aget_name(attr_id, size, name);
	H5Aclose(attr_id);
	if (n < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute name");

	return n;
#endif
}

/**
 * @brief Adds a scalar attribute to a #LALH5File or #LALH5Dataset
 * @details
 * This routine adds a scalar-valued attribute with name @p key
 * and value given by the memory addressed by @p value to a
 * HDF5 object associated with the #LALH5File or #LALH5Dataset
 * @p object.
 * The data type of the scalar value is given by the #LALTYPECODE
 * @p dtype.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to which the attribute
 * will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to the value of the scalar attribute to be added.
 * @param dtype #LALTYPECODE value specifying the data type of the attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5AttributeAddScalar(LALH5Generic UNUSED object, const char UNUSED *key, const void UNUSED *value, LALTYPECODE UNUSED dtype)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	dtype_id = XLALH5TypeFromLALType(dtype);
	if (dtype_id < 0)
		XLAL_ERROR(XLAL_EFUNC);

	space_id = H5Screate(H5S_SCALAR);
	if (space_id < 0) {
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create dataspace for attribute `%s'", key);
	}

	attr_id = H5Acreate2(obj_id, key, dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
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
#endif
}

/**
 * @brief Adds a string attribute to a #LALH5File or #LALH5Dataset
 * @details
 * This routine adds a NUL-terminated variable-length string @p value
 * attribute with name @p key to a HDF5 object associated with the
 * #LALH5File or #LALH5Dataset @p object.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to which the
 * attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a string with the value of the new attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5AttributeAddString(LALH5Generic UNUSED object, const char UNUSED *key, const char UNUSED *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

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

	attr_id = H5Acreate2(obj_id, key, dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
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
#endif
}

/**
 * @brief Adds a LIGOTimeGPS attribute to a #LALH5File or #LALH5Dataset
 * @details
 * This routine adds a LIGOTimeGPS @p value attribute with name @p key to a
 * HDF5 object associated with the #LALH5File or #LALH5Dataset @p object.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to which the
 * attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a LIGOTimeGPS structure with the value of the new
 * attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5AttributeAddLIGOTimeGPS(LALH5Generic UNUSED object, const char UNUSED *key, const LIGOTimeGPS UNUSED *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	dtype_id = XLALH5TypeNativeLIGOTimeGPS();
	if (dtype_id < 0)
		XLAL_ERROR(XLAL_EFUNC);

	space_id = H5Screate(H5S_SCALAR);
	if (space_id < 0) {
		H5Tclose(dtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create dataspace for attribute `%s'", key);
	}

	attr_id = H5Acreate2(obj_id, key, dtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
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
#endif
}

/**
 * @brief Adds a 1d enum array attribute to a #LALH5File or #LALH5Dataset
 * @details
 * This routine adds the 1d enum array @p value of length @p length attribute
 * with name @p key to a HDF5 object associated with the #LALH5File or
 * #LALH5Dataset @p object.
 * The names and values of the @p nenum enumeration constants are provided
 * in the arrays @p enumnames and @p enumvals respectively.
 * 
 * @param object Pointer to a #LALH5File or #LALH5Dataset to which the
 * attribute will be added.
 * @param enumnames Pointer to an array of names of the enum constants.
 * @param enumvals Pointer to an array of values of the enum constants.
 * @param nenum Number of enum constants.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to an array of enum values.
 * @param length Length of the array of enum values.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5AttributeAddEnumArray1D(LALH5Generic UNUSED object, const char UNUSED *enumnames[], const int UNUSED enumvals[], size_t UNUSED nenum, const char UNUSED *key, const int UNUSED value[], size_t UNUSED length)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hsize_t dims[] = {length};
	hid_t obj_id;
	hid_t attr_id;
	hid_t space_id;
	hid_t edtype_id;
	hid_t adtype_id;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	edtype_id = XLALH5TypeEnum(enumnames, enumvals, nenum);
	if (edtype_id < 0)
		XLAL_ERROR(XLAL_EFUNC, "Could not create enum type for attribute `%s'", key);

	adtype_id = H5Tarray_create2(edtype_id, 1, dims);
	H5Tclose(edtype_id);
	if (adtype_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not create array type for attribute `%s'", key);

	space_id = H5Screate(H5S_SCALAR);
	if (space_id < 0) {
		H5Tclose(adtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create dataspace for attribute `%s'", key);
	}

	attr_id = H5Acreate2(obj_id, key, adtype_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(space_id);
	if (attr_id < 0) {
		H5Tclose(adtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not create attribute `%s'", key);
	}

	if (H5Awrite(attr_id, adtype_id, value) < 0) {
		H5Aclose(attr_id);
		H5Tclose(adtype_id);
		XLAL_ERROR(XLAL_EIO, "Could not write attribute `%s'", key);
	}
	
	H5Aclose(attr_id);
	H5Tclose(adtype_id);
	return 0;
#endif
}

/**
 * @brief Gets the datatype of an attribute in a #LALH5File or #LALH5Dataset
 * @details
 * This routine queries the datatype of a scalar attribute with
 * name @p key in a HDF5 object associated with the #LALH5File or #LALH5Dataset
 * @p object.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns #LALTYPECODE value of the datatype of the scalar attribute.
 * @retval -1 Failure.
 */
LALTYPECODE XLALH5AttributeQueryScalarType(const LALH5Generic UNUSED object, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALTYPECODE dtype;
	hid_t obj_id;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;
	hid_t memtype_id;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	attr_id = H5Aopen(obj_id, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute `%s'", key);

	/* sanity check: make sure this is a scalar attribute */
	space_id = H5Aget_space(attr_id);
	if (space_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read dataspace of attribute `%s'", key);
	}
	if (H5Sget_simple_extent_ndims(space_id) != 0 || H5Sget_simple_extent_npoints(space_id) != 1) {
		H5Sclose(space_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Attribute `%s' is not a scalar attribute", key);
	}

	dtype_id = H5Aget_type(attr_id);
	H5Aclose(attr_id);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read type of attribute `%s'", key);

	memtype_id = H5Tget_native_type(dtype_id, H5T_DIR_ASCEND);
	H5Tclose(dtype_id);
	if (memtype_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not get native type of attribute `%s'", key);

	dtype = XLALTypeFromH5Type(memtype_id);
	H5Tclose(memtype_id);
	if ((int)dtype < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return dtype;
#endif
}

/**
 * @brief Gets the value of a scalar attribute in a #LALH5File or #LALH5Dataset
 * @details
 * This routine queries the value of a scalar attribute with
 * name @p key in a HDF5 object associated with the #LALH5File or #LALH5Dataset
 * @p object.
 * The value is stored in memory pointed to by the pointer @p value.
 *
 * @attention
 * This routine does not allocate memory for @p value.  The calling
 * routine must ensure that the memory addressed by the pointer @p value
 * is sufficient to hold the value in the attribute.
 *
 * @param value Pointer to memory in which the value will be stored.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5AttributeQueryScalarValue(void UNUSED *value, const LALH5Generic UNUSED object, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;
	hid_t memtype_id;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	attr_id = H5Aopen(obj_id, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute `%s'", key);

	/* sanity check: make sure this is a scalar attribute */
	space_id = H5Aget_space(attr_id);
	if (space_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read dataspace of attribute `%s'", key);
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
		XLAL_ERROR(XLAL_EIO, "Could not read type of attribute `%s'", key);
	}

	memtype_id = H5Tget_native_type(dtype_id, H5T_DIR_ASCEND);
	H5Tclose(dtype_id);
	if (memtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not get native type of attribute `%s'", key);
	}

	if (H5Aread(attr_id, memtype_id, value) < 0) {
		H5Tclose(memtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read data from attribute `%s'", key);
	}

	H5Tclose(memtype_id);
	H5Aclose(attr_id);
	return 0;
#endif
}

/**
 * @brief Gets the value of a string attribute in a #LALH5File or #LALH5Dataset
 * @details
 * This routine queries the value of a string attribute with
 * name @p key in a HDF5 object associated with the #LALH5File or #LALH5Dataset
 * @p object.
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
 * @param object Pointer to a #LALH5File or #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns The number of bytes that would be written to @p value had @p size
 * been sufficiently large excluding the terminating NUL byte.
 * @retval -1 Failure.
 */
int XLALH5AttributeQueryStringValue(char UNUSED *value, size_t UNUSED size, const LALH5Generic UNUSED object, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	char *str;
	hid_t obj_id;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;
	hid_t memtype_id;
	int n;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	attr_id = H5Aopen(obj_id, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute `%s'", key);

	/* sanity check: make sure this is just one variable-length string */
	space_id = H5Aget_space(attr_id);
	if (space_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read dataspace of attribute `%s'", key);
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
		XLAL_ERROR(XLAL_EIO, "Could not read type of attribute `%s'", key);
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
		XLAL_ERROR(XLAL_EIO, "Could not read data from attribute `%s'", key);
	}

	n = snprintf(value, value == NULL ? 0 : size, "%s", str);

	H5Dvlen_reclaim(memtype_id, space_id, H5P_DEFAULT, &str);

	H5Tclose(memtype_id);
	H5Sclose(space_id);
	H5Aclose(attr_id);
	return n;
#endif
}

/**
 * @brief Gets the value of a LIGOTimeGPS attribute in a #LALH5File or #LALH5Dataset
 * @details
 * This routine queries the value of a LIGOTimeGPS attribute with
 * name @p key in a HDF5 object associated with the #LALH5File or #LALH5Dataset
 * @p object.
 * The value is stored in memory pointed to by the pointer @p value.
 * @param value Pointer to a LIGOTimeGPS structure in which the attribute
 * value will be stored.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns Pointer to the LIGOTimeGPS structure passed to this routine.
 * @retval NULL Failure.
 */
LIGOTimeGPS * XLALH5AttributeQueryLIGOTimeGPSValue(LIGOTimeGPS UNUSED *value, const LALH5Generic UNUSED object, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	char *s;
	hid_t obj_id;
	hid_t attr_id;
	hid_t space_id;
	hid_t dtype_id;
	hid_t memtype_id;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR_NULL(XLAL_EINVAL);

	attr_id = H5Aopen(obj_id, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read attribute `%s'", key);

	/* sanity check: make sure this is a scalar attribute */
	space_id = H5Aget_space(attr_id);
	if (space_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read dataspace of attribute `%s'", key);
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
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read type of attribute `%s'", key);
	}

	/* must be compound data type ... */
	if (H5Tget_class(dtype_id) != H5T_COMPOUND) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}

	/* must have 2 members ... */
	if (H5Tget_nmembers(dtype_id) != 2) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}

	/* both members must be integer type ... */
	if (H5Tget_member_class(dtype_id, 0) != H5T_INTEGER) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}
	if (H5Tget_member_class(dtype_id, 1) != H5T_INTEGER) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}

	/* FIXME: should also check member sizes? */

	/* first member name must be "gpsSeconds" ... */
	s = H5Tget_member_name(dtype_id, 0);
	if (strcmp(s, "gpsSeconds") != 0) {
		free(s);
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}
	free(s);

	/* second member name must be "gpsNanoSeconds" ... */
	s = H5Tget_member_name(dtype_id, 1);
	if (strcmp(s, "gpsNanoSeconds") != 0) {
		free(s);
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}
	free(s);

	H5Tclose(dtype_id);

	memtype_id = XLALH5TypeNativeLIGOTimeGPS();
	if (memtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not get native type of attribute `%s'", key);
	}

	if (H5Aread(attr_id, memtype_id, value) < 0) {
		H5Tclose(memtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read data from attribute `%s'", key);
	}

	H5Tclose(memtype_id);
	H5Aclose(attr_id);
	return value;
#endif
}

/**
 * @brief Gets the length of a 1D enum array attribute in a #LALH5File or
 * #LALH5Dataset
 * @details
 * This routine queries the length of a 1D enum array with name @p key in a
 * HDF5 object associated with the #LALH5File or #LALH5Dataset @p object.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns Length of the 1D enum array.
 * @retval -1 Failure.
 */
size_t XLALH5AttributeQueryEnumArray1DLength(const LALH5Generic UNUSED object, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hsize_t dims[] = {-1};
	hid_t obj_id;
	hid_t attr_id;
	hid_t dtype_id;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	attr_id = H5Aopen(obj_id, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute `%s'", key);

	dtype_id = H5Aget_type(attr_id);
	if (dtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type of attribute `%s'", key);
	}

	/* must be array data type ... */
	if (H5Tget_class(dtype_id) != H5T_ARRAY) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}

	/* must be a 1d array data type ... */
	if (H5Tget_array_ndims(dtype_id) != 1) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}

	if (H5Tget_array_dims2(dtype_id, dims) < 0) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_ETYPE, "Could not get array length of attribute `%s'", key);
	}

	H5Tclose(dtype_id);
	H5Aclose(attr_id);
	return dims[0];
#endif
}

/**
 * @brief Gets the values in a 1D enum array attribute in a #LALH5File or
 * #LALH5Dataset
 * @details
 * This routine reads the values of a 1D enum array with name @p key in a
 * HDF5 object associated with the #LALH5File or #LALH5Dataset @p object
 * into the array @p value.
 *
 * @attention
 * This routine does not allocate memory for @p value.  The calling
 * routine must ensure that the memory addressed by the pointer @p value
 * is sufficient to hold the value in the attribute.
 * The routine XLALH5AttributeQueryEnumArray1DLength() will yield the
 * length that this array must have.
 * 
 * @param value Pointer to an array where then enum values will be stored.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5AttributeQueryEnumArray1DValue(int UNUSED value[], const LALH5Generic UNUSED object, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	hid_t attr_id;
	hid_t dtype_id;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EINVAL);

	attr_id = H5Aopen(obj_id, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute `%s'", key);

	dtype_id = H5Aget_type(attr_id);
	if (dtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type of attribute `%s'", key);
	}

	/* must be array data type ... */
	if (H5Tget_class(dtype_id) != H5T_ARRAY) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}

	/* must be a 1d array data type ... */
	if (H5Tget_array_ndims(dtype_id) != 1) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_ETYPE, "Incorrect data type for attribute `%s'", key);
	}

	if (H5Aread(attr_id, dtype_id, value) < 0) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_ETYPE, "Could not read data from attribute `%s'", key);
	}

	H5Tclose(dtype_id);
	H5Aclose(attr_id);
	return 0;
#endif
}

/**
 * @brief Gets the number of constants in an enum type associated with an
 * attribute in a #LALH5File or #LALH5Dataset
 * @details
 * This routine queries the number constants in an enum type associated
 * with the attribute with name @p key in a
 * HDF5 object associated with the #LALH5File or #LALH5Dataset @p object.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns Number of constants in enum type.
 * @retval -1 Failure.
 */
size_t XLALH5AttributeQueryNEnum(const LALH5Generic UNUSED object, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	hid_t attr_id;
	hid_t dtype_id;
	H5T_class_t dtype_class;
	int nenum;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EFUNC);

	attr_id = H5Aopen(obj_id, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute `%s'", key);

	dtype_id = H5Aget_type(attr_id);
	if (dtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type of attribute `%s'", key);
	}

	dtype_class = H5Tget_class(dtype_id);
	if (dtype_class == H5T_ARRAY) {
		/* get base data type */
		hid_t super_id;
		super_id = H5Tget_super(dtype_id);
		H5Tclose(dtype_id);
		dtype_id = super_id;
		dtype_class = H5Tget_class(dtype_id);
	}
	if (dtype_class == H5T_NO_CLASS) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type class of attribute `%s'", key);
	}
	if (dtype_class != H5T_ENUM) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Attribute `%s' is not an enum type", key);
	}

	nenum = H5Tget_nmembers(dtype_id);
	if (nenum < 0) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read number of enum constants of attribute `%s'", key);
	}

	H5Tclose(dtype_id);
	H5Aclose(attr_id);
	return nenum;
#endif
}

/**
 * @brief Gets the name of a constants in an enum type associated with an
 * attribute in a #LALH5File or #LALH5Dataset
 * @details
 * This routine queries the name of a constant in an enum type associated
 * with the attribute with name @p key in a
 * HDF5 object associated with the #LALH5File or #LALH5Dataset @p object.
 * The index @p pos identifies which constant's name is returned.
 * The result is written into the buffer pointed to by @p name, the size
 * of which is @p size bytes.  If @p name is NULL, no data is copied but
 * the routine returns the length of the string.  Therefore, this routine
 * can be called once to determine the amount of memory required, the
 * memory can be allocated, and then it can be called a second time to
 * read the string.  If the parameter @p size is less than or equal to
 * the string length then only $p size-1 bytes of the string are copied
 * to the buffer @p name.
 * @note The return value is the length of the string, not including the
 * terminating NUL character; thus the buffer @p name should be allocated
 * to be one byte larger.
 * @param name Pointer to a buffer into which the string will be written.
 * @param size Size in bytes of the buffer into which the string will be
 * written.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @param pos The index identifying which enum constant.
 * @returns The number of bytes that would be written to @p name had @p size
 * been sufficiently large excluding the terminating NUL byte.
 * @retval -1 Failure.
 */
int XLALH5AttributeQueryEnumName(char UNUSED *name, size_t UNUSED size, const LALH5Generic UNUSED object, const char UNUSED *key, int UNUSED pos)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	hid_t attr_id;
	hid_t dtype_id;
	H5T_class_t dtype_class;
	int nenum;
	int n;
	char *s;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EFUNC);

	attr_id = H5Aopen(obj_id, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute `%s'", key);

	dtype_id = H5Aget_type(attr_id);
	if (dtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type of attribute `%s'", key);
	}

	dtype_class = H5Tget_class(dtype_id);
	if (dtype_class == H5T_ARRAY) {
		/* get base data type */
		hid_t super_id;
		super_id = H5Tget_super(dtype_id);
		H5Tclose(dtype_id);
		dtype_id = super_id;
		dtype_class = H5Tget_class(dtype_id);
	}
	if (dtype_class == H5T_NO_CLASS) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type class of attribute `%s'", key);
	}
	if (dtype_class != H5T_ENUM) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Attribute `%s' is not an enum type", key);
	}

	nenum = H5Tget_nmembers(dtype_id);
	if (nenum < 0) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read number of enum constants of attribute `%s'", key);
	}
	if (pos >= nenum) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EINVAL, "Requested enum constant of attribute `%s' does not exist", key);
	}

	s = H5Tget_member_name(dtype_id, pos);
	if (s == NULL) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read enum constant name of attribute `%s'", key);
	}

	n = snprintf(name, name == NULL ? 0 : size, "%s", s);

	free(s);
	H5Tclose(dtype_id);
	H5Aclose(attr_id);
	return n;
#endif
}

/**
 * @brief Gets the value of a constants in an enum type associated with an
 * attribute in a #LALH5File or #LALH5Dataset
 * @details
 * This routine queries the value of a constant in an enum type associated
 * with the attribute with name @p key in a
 * HDF5 object associated with the #LALH5File or #LALH5Dataset @p object.
 * The index @p pos identifies which constant's name is returned.
 * @param object Pointer to a #LALH5File or #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @param pos The index identifying which enum constant.
 * @returns The value of the enum constant.
 * @retval -1 Failure.  Note that -1 might also be a valid enum constant!
 */
int XLALH5AttributeQueryEnumValue(const LALH5Generic UNUSED object, const char UNUSED *key, int UNUSED pos)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_VAL(INT_MAX, XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t obj_id;
	hid_t attr_id;
	hid_t dtype_id;
	H5T_class_t dtype_class;
	int nenum;
	int val;

	obj_id = object.generic->object_id;
	if (obj_id < 0)
		XLAL_ERROR(XLAL_EFUNC);

	attr_id = H5Aopen(obj_id, key, H5P_DEFAULT);
	if (attr_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read attribute `%s'", key);

	dtype_id = H5Aget_type(attr_id);
	if (dtype_id < 0) {
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type of attribute `%s'", key);
	}

	dtype_class = H5Tget_class(dtype_id);
	if (dtype_class == H5T_ARRAY) {
		/* get base data type */
		hid_t super_id;
		super_id = H5Tget_super(dtype_id);
		H5Tclose(dtype_id);
		dtype_id = super_id;
		dtype_class = H5Tget_class(dtype_id);
	}
	if (dtype_class == H5T_NO_CLASS) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read type class of attribute `%s'", key);
	}
	if (dtype_class != H5T_ENUM) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Attribute `%s' is not an enum type", key);
	}

	nenum = H5Tget_nmembers(dtype_id);
	if (nenum < 0) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read number of enum constants of attribute `%s'", key);
	}
	if (pos >= nenum) {
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EINVAL, "Requested enum constant of attribute `%s' does not exist", key);
	}

	if (H5Tget_member_value(dtype_id, pos, &val) < 0)
	{
		H5Tclose(dtype_id);
		H5Aclose(attr_id);
		XLAL_ERROR(XLAL_EIO, "Could not read enum constant value of attribute `%s'", key);
	}

	H5Tclose(dtype_id);
	H5Aclose(attr_id);
	return val;
#endif
}

/** @} */

/**
 * @name File Attribute Routines
 *
 * @deprecated
 * Use the general @ref attribute_routines "Attribute Routines"
 * instead.
 *
 * @{
 */

/**
 * @brief DEPRECATED: Adds a scalar attribute to a #LALH5File
 * @details
 * This routine adds a scalar-valued attribute with name @p key
 * and value given by the memory addressed by @p value to a
 * HDF5 file associated with the #LALH5File @p file.
 * The data type of the scalar value is given by the #LALTYPECODE
 * @p dtype.
 *
 * @deprecated
 * Use XLALH5AttributeAddScalar() instead.
 * 
 * @param file Pointer to a #LALH5File to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to the value of the scalar attribute to be added.
 * @param dtype #LALTYPECODE value specifying the data type of the attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileAddScalarAttribute(LALH5File UNUSED *file, const char UNUSED *key, const void UNUSED *value, LALTYPECODE UNUSED dtype)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
        XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeAddScalar");
	if (XLALH5AttributeAddScalar((LALH5Generic)file, key, value, dtype) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief DEPRECATED: Gets attribute names from a #LALH5File
 * @details
 * This routine returns the names of all attributes from a HDF5 Dataset
 *
 * @deprecated
 * Use XLALH5AttributeQueryN() and XLALH5AttributeQueryName() instead.
 *
 * @param names Pointer a list of strings to be returned to the user. Memory
 * should be freed by the caller.
 * @param file Pointer to a #LALH5File from which the attributes will be added.
 * @param N Pointer to a UINT4 where the number of datasets will be recorded
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileGetAttributeNames(LALH5File UNUSED *file, char UNUSED *** names, UINT4 UNUSED *N)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
    int na;
    int ns;
    int i;
    hid_t aid;

	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeQueryN() and XLALH5AttributeQueryName");
	if (file == NULL)
		XLAL_ERROR(XLAL_EFAULT);
	if (names == NULL)
		XLAL_ERROR(XLAL_EFAULT);

    na = H5Aget_num_attrs(file->file_id);
    char ** namelist = (char**) XLALMalloc(na * sizeof(*names));

    for (i = 0; i <na; i++) {
        aid = H5Aopen_idx(file->file_id, (unsigned int)i);
        ns = H5Aget_name(aid, 0, NULL) + 1;
        namelist[i] = (char*) XLALMalloc(ns * sizeof(namelist[i]));
        H5Aget_name(aid, ns, namelist[i]);
        H5Aclose(aid);
    }

  *N=na;
  *names = namelist;
  return(XLAL_SUCCESS);
#endif
}


/**
 * @brief DEPRECATED: Adds a string attribute to a #LALH5File
 * @details
 * This routine adds a NUL-terminated variable-length string @p value
 * attribute with name @p key to a HDF5 file associated with the
 * #LALH5File @p file.
 *
 * @deprecated
 * Use XLALH5AttributeAddString() instead.
 *
 * @param file Pointer to a #LALH5File to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a string with the value of the new attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileAddStringAttribute(LALH5File UNUSED *file, const char UNUSED *key, const char UNUSED *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeAddString");
	if (XLALH5AttributeAddString((LALH5Generic)file, key, value) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief DEPRECATED: Adds a LIGOTimeGPS attribute to a #LALH5File
 * @details
 * This routine adds a LIGOTimeGPS @p value attribute with name @p key to a
 * HDF5 file associated with the #LALH5File @p file.
 *
 * @deprecated
 * Use XLALH5AttributeAddLIGOTimeGPS() instead.
 *
 * @param file Pointer to a #LALH5File to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a LIGOTimeGPS structure with the value of the new
 * attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileAddLIGOTimeGPSAttribute(LALH5File UNUSED *file, const char UNUSED *key, const LIGOTimeGPS UNUSED *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeAddLIGOTimeGPS");
	if (XLALH5AttributeAddLIGOTimeGPS((LALH5Generic)file, key, value) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief DEPRECATED: Gets the datatype of an attribute in a #LALH5File
 * @details
 * This routine queries the datatype of a scalar attribute with
 * name @p key in a HDF5 file associated with the #LALH5File @p file.
 *
 * @deprecated
 * Use XLALH5AttributeQueryScalarType() instead.
 *
 * @param file Pointer to a #LALH5File to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns #LALTYPECODE value of the datatype of the scalar attribute.
 * @retval -1 Failure.
 */
LALTYPECODE XLALH5FileQueryScalarAttributeType(LALH5File UNUSED *file, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALTYPECODE dtype;

	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeQueryScalarType");
	dtype = XLALH5AttributeQueryScalarType((LALH5Generic)file, key);
	if ((int)(dtype) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return dtype;
#endif
}

/**
 * @brief DEPRECATED: Gets the value of a scalar attribute in a #LALH5File
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
 * @deprecated
 * Use XLALH5AttributeQueryScalarValue() instead.
 *
 * @param value Pointer to memory in which the value will be stored.
 * @param file Pointer to a #LALH5File to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5FileQueryScalarAttributeValue(void UNUSED *value, LALH5File UNUSED *file, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeQueryScalarValue");
	if (XLALH5AttributeQueryScalarValue(value, (LALH5Generic)file, key) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief DEPRECATED: Gets the value of a string attribute in a #LALH5File
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
 *
 * @deprecated
 * Use XLALH5AttributeQueryStringValue() instead.
 *
 * @param value Pointer to a buffer into which the string will be written.
 * @param size Size in bytes of the buffer into which the string will be
 * written.
 * @param file Pointer to a #LALH5File to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns The number of bytes that would be written to @p value had @p size
 * been sufficiently large excluding the terminating NUL byte.
 * @retval -1 Failure.
 */
int XLALH5FileQueryStringAttributeValue(char UNUSED *value, size_t UNUSED size, LALH5File UNUSED *file, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	int n;

	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeQueryStringValue");
	n = XLALH5AttributeQueryStringValue(value, size, (LALH5Generic)file, key);
	if (n < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return n;
#endif
}

/**
 * @brief DEPRECATED: Gets the value of a LIGOTimeGPS attribute in a #LALH5File
 * @details
 * This routine queries the value of a LIGOTimeGPS attribute with
 * name @p key in a HDF5 file associated with the #LALH5File @p file.
 * The value is stored in memory pointed to by the pointer @p value.
 *
 * @deprecated
 * Use XLALH5AttributeQueryLIGOTimeGPSValue() instead.
 *
 * @param value Pointer to a LIGOTimeGPS structure in which the attribute
 * value will be stored.
 * @param file Pointer to a #LALH5File to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns Pointer to the LIGOTimeGPS structure passed to this routine.
 * @retval NULL Failure.
 */
LIGOTimeGPS * XLALH5FileQueryLIGOTimeGPSAttributeValue(LIGOTimeGPS UNUSED *value, LALH5File UNUSED *file, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeQueryLIGOTimeGPSValue");
	if (XLALH5AttributeQueryLIGOTimeGPSValue(value, (LALH5Generic)file, key) == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return value;
#endif
}

/** @} */

/**
 * @name Dataset Attribute Routines
 *
 * @deprecated
 * Use the general @ref attribute_routines "Attribute Routines"
 * instead.
 *
 * @{
 */

/**
 * @brief DEPRECATED: Adds a scalar attribute to a #LALH5Dataset
 * @details
 * This routine adds a scalar-valued attribute with name @p key
 * and value given by the memory addressed by @p value to a
 * HDF5 dataset associated with the #LALH5Dataset @p dset.
 * The data type of the scalar value is given by the #LALTYPECODE
 * @p dtype.
 *
 * @deprecated
 * Use XLALH5AttributeAddScalar() instead.
 * 
 * @param dset Pointer to a #LALH5Dataset to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to the value of the scalar attribute to be added.
 * @param dtype #LALTYPECODE value specifying the data type of the attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetAddScalarAttribute(LALH5Dataset UNUSED *dset, const char UNUSED *key, const void UNUSED *value, LALTYPECODE UNUSED dtype)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
        XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeAddScalar");
	if (XLALH5AttributeAddScalar((LALH5Generic)dset, key, value, dtype) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief DEPRECATED: Adds a string attribute to a #LALH5Dataset
 * @details
 * This routine adds a NUL-terminated variable-length string @p value
 * attribute with name @p key to a HDF5 dataset associated with the
 * #LALH5Dataset @p dset.
 *
 * @deprecated
 * Use XLALH5AttributeAddString() instead.
 *
 * @param dset Pointer to a #LALH5Dataset to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a string with the value of the new attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetAddStringAttribute(LALH5Dataset UNUSED *dset, const char UNUSED *key, const char UNUSED *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeAddString");
	if (XLALH5AttributeAddString((LALH5Generic)dset, key, value) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief DEPRECATED: Adds a LIGOTimeGPS attribute to a #LALH5Dataset
 * @details
 * This routine adds a LIGOTimeGPS @p value attribute with name @p key to a
 * HDF5 dataset associated with the #LALH5Dataset @p dset.
 *
 * @deprecated
 * Use XLALH5AttributeAddLIGOTimeGPS() instead.
 *
 * @param dset Pointer to a #LALH5Dataset to which the attribute will be added.
 * @param key Pointer to a string with the name of the new attribute.
 * @param value Pointer to a LIGOTimeGPS structure with the value of the new
 * attribute.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetAddLIGOTimeGPSAttribute(LALH5Dataset UNUSED *dset, const char UNUSED *key, const LIGOTimeGPS UNUSED *value)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeAddLIGOTimeGPS");
	if (XLALH5AttributeAddLIGOTimeGPS((LALH5Generic)dset, key, value) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief DEPRECATED: Gets the datatype of an attribute in a #LALH5Dataset
 * @details
 * This routine queries the datatype of a scalar attribute with
 * name @p key in a HDF5 dataset associated with the #LALH5Dataset @p dset.
 *
 * @deprecated
 * Use XLALH5AttributeQueryScalarType() instead.
 *
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns #LALTYPECODE value of the datatype of the scalar attribute.
 * @retval -1 Failure.
 */
LALTYPECODE XLALH5DatasetQueryScalarAttributeType(LALH5Dataset UNUSED *dset, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	LALTYPECODE dtype;

	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeQueryScalarType");
	dtype = XLALH5AttributeQueryScalarType((LALH5Generic)dset, key);
	if ((int)(dtype) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return dtype;
#endif
}

/**
 * @brief DEPRECATED: Gets the value of a scalar attribute in a #LALH5Dataset
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
 * @deprecated
 * Use XLALH5AttributeQueryScalarValue() instead.
 *
 * @param value Pointer to memory in which the value will be stored.
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5DatasetQueryScalarAttributeValue(void UNUSED *value, LALH5Dataset UNUSED *dset, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeQueryScalarValue");
	if (XLALH5AttributeQueryScalarValue(value, (LALH5Generic)dset, key) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return 0;
#endif
}

/**
 * @brief DEPRECATED: Gets the value of a string attribute in a #LALH5Dataset
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
 *
 * @deprecated
 * Use XLALH5AttributeQueryStringValue() instead.
 *
 * @param value Pointer to a buffer into which the string will be written.
 * @param size Size in bytes of the buffer into which the string will be
 * written.
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns The number of bytes that would be written to @p value had @p size
 * been sufficiently large excluding the terminating NUL byte.
 * @retval NULL Failure.
 */
int XLALH5DatasetQueryStringAttributeValue(char UNUSED *value, size_t UNUSED size, LALH5Dataset UNUSED *dset, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	int n;

	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeQueryStringValue");
	n = XLALH5AttributeQueryStringValue(value, size, (LALH5Generic)dset, key);
	if (n < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return n;
#endif
}

/**
 * @brief DEPRECATED: Gets the value of a LIGOTimeGPS attribute in a #LALH5Dataset
 * @details
 * This routine queries the value of a LIGOTimeGPS attribute with
 * name @p key in a HDF5 dataset associated with the #LALH5Dataset @p dset.
 * The value is stored in memory pointed to by the pointer @p value.
 *
 * @deprecated
 * Use XLALH5AttributeQueryLIGOTimeGPSValue() instead.
 *
 * @param value Pointer to a LIGOTimeGPS structure in which the attribute
 * value will be stored.
 * @param dset Pointer to a #LALH5Dataset to be queried.
 * @param key Pointer to a string with the name of the attribute to query.
 * @returns Pointer to the LIGOTimeGPS structure passed to this routine.
 * @retval NULL Failure.
 */
LIGOTimeGPS * XLALH5DatasetQueryLIGOTimeGPSAttributeValue(LIGOTimeGPS UNUSED *value, LALH5Dataset UNUSED *dset, const char UNUSED *key)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	XLAL_PRINT_DEPRECATION_WARNING("XLALH5AttributeQueryLIGOTimeGPSValue");
	if (XLALH5AttributeQueryLIGOTimeGPSValue(value, (LALH5Generic)dset, key) == NULL)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return value;
#endif
}

/** @} */

/** @} */

/**
 * @addtogroup H5Table_group
 * @brief Routines for reading/writing tables in HDF5 files.
 * @details
 * These routines are basic routines for reading and writing
 * tables in HDF5 files.
 * @{
 */

/**
 * @brief Allocates a #LALH5Dataset dataset to hold a table.
 * @details
 * This routine creates a dataset with name @p name within an HDF5
 * file associated with the #LALH5File @p file structure and allocates
 * a #LALH5Dataset structure associated with the dataset.  The type
 * of data to be stored are table data comprised of rows of size
 * @p rowsz; each row having @p ncols columns with names @p cols, types @p
 * types, and offsets @p offsets.  

 * with the #LALH5Dataset @p dset and stores the data in the buffer @p data.
 * This buffer should be sufficiently large to hold the requested rows
 * of the dataset, which is comprised of a number of rows each of size
 * @p rowsz.  Each row is comprised of a number of columns
 * which are named in the comma-separated list in string @p cols.
 * The offsets of the columns and the sizes of the columns in each row
 * of data is provided by the arrays @p offsets and @p colsz respectively.
 *
 * The following example shows how to create a dataset named @a particles
 * in which particle data can be stored.
 * @code
 * #include <stddef.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/LALH5FileIO.h>
 *
 * struct body {REAL8 x; REAL8 y; REAL8 z; REAL8 vx; REAL8 vy; REAL8 vz;};
 * size_t ncols = 6;
 * size_t types[6] = {LAL_D_TYPE_CODE, LAL_D_TYPE_CODE, LAL_D_TYPE_CODE, LAL_D_TYPE_CODE, LAL_D_TYPE_CODE, LAL_D_TYPE_CODE};
 * size_t offsets[6] = {offsetof(struct body, x), offsetof(struct body, y), offsetof(struct body, z), offsetof(struct body, vx), offsetof(struct body, vy), offsetof(struct body, vz)};
 * size_t rowsz = sizeof(*data);
 * LALH5File *file = XLALH5FileOpen("example.h5", "w");
 * LALH5Dataset *dset = XLALH5TableAlloc(file, "particles", ncols, cols, types, offsets, rowsz);
 * @endcode
 *
 * @param file Pointer to a #LALH5File in which to create the dataset.
 * @param name Pointer to a string with the name of the dataset to create (also
 * the table name).
 * @param ncols Number of columns in each row.
 * @param cols Pointer to an array of strings giving the column names.
 * @param types Pointer to an array of #LALTYPECODE values specifying the data
 * type of each column.
 * @param offsets Pointer to an array of offsets for each column.
 * @param rowsz Size of each row of data.
 * @returns A pointer to a #LALH5Dataset structure associated with the specified
 * dataset within a HDF5 file.
 * @retval NULL An error occurred creating the dataset.
 */
LALH5Dataset * XLALH5TableAlloc(LALH5File UNUSED *file, const char UNUSED *name, size_t UNUSED ncols, const char UNUSED **cols, const LALTYPECODE UNUSED *types, const size_t UNUSED *offsets, size_t UNUSED rowsz)
{
#ifndef HAVE_HDF5
	XLAL_ERROR_NULL(XLAL_EFAILED, "HDF5 support not implemented");
#else
	const size_t chunk_size = 32;
	hid_t dtype_id[ncols];
	hid_t tdtype_id;
	size_t col;
	size_t namelen;
	herr_t status;
	LALH5Dataset *dset;

	if (file == NULL || cols == NULL || types == NULL || offsets == NULL)
		XLAL_ERROR_NULL(XLAL_EFAULT);

	if (file->mode != LAL_H5_FILE_MODE_WRITE)
		XLAL_ERROR_NULL(XLAL_EINVAL, "Attempting to write to a read-only HDF5 file");

	/* map the LAL types to HDF5 types */
	for (col = 0; col < ncols; ++col) {
		dtype_id[col] = XLALH5TypeFromLALType(types[col]);
		if (dtype_id[col] < 0)
			XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* make empty table */
	/* note: table title and dataset name are the same */
	status = H5TBmake_table(name, file->file_id, name, ncols, 0, rowsz, cols, offsets, dtype_id, chunk_size, NULL, 0, NULL);
	for (col = 0; col < ncols; ++col)
		H5Tclose(dtype_id[col]);
	
	if (status < 0)
		XLAL_ERROR_NULL(XLAL_EIO);

	/* now read dataset from file to return */

	namelen = strlen(name);
	dset = LALCalloc(1, sizeof(*dset) + namelen + 1);  /* use flexible array member to record name */
	if (!dset)
		XLAL_ERROR_NULL(XLAL_ENOMEM);

	dset->dataset_id = H5Dopen2(file->file_id, name, H5P_DEFAULT);
	if (dset->dataset_id < 0) {
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read dataset `%s'", name);
	}

	dset->space_id = H5Dget_space(dset->dataset_id);
	if (dset->space_id < 0) {
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read dataspace of dataset `%s'", name);
	}

	tdtype_id = H5Dget_type(dset->dataset_id);
	if (tdtype_id < 0) {
		H5Sclose(dset->space_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not read datatype of dataset `%s'", name);
	}

	/* convert type to native type */
	dset->dtype_id = H5Tget_native_type(tdtype_id, H5T_DIR_ASCEND);
	H5Tclose(tdtype_id);
	if (dset->dtype_id < 0) {
		H5Sclose(dset->space_id);
		LALFree(dset);
		XLAL_ERROR_NULL(XLAL_EIO, "Could not get native type for dataset `%s'", name);
	}

	/* record name of dataset and parent id */
	snprintf(dset->name, namelen + 1, "%s", name);
	dset->parent_id = file->file_id;
	return dset;
#endif
}

/**
 * @brief Appends rows of data to a #LALH5Dataset dataset holding a table.
 * @details
 * This routine appends @p nrows rows of data each having size @p rowsz
 * to a HDF5 table dataset associated with the #LALH5Dataset structure
 * @p dset.  The data is contained in @p data and the offsets and sizes
 * of each column are given by the arrays @p offsets and @p colsz.
 *
 * The following example shows how to append rows to a dataset named
 * @a particles in which particle data is stored.
 * @code
 * #include <stddef.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/LALH5FileIO.h>
 *
 * struct body {REAL8 x; REAL8 y; REAL8 z; REAL8 vx; REAL8 vy; REAL8 vz;} data[] = {
 * 	{ .x = +0.83, .y = -0.13, .z = -0.24, .vx = +0.60, .vy = -0.79, .vz = -0.98 },
 * 	{ .x = -0.59, .y = -0.12, .z = -0.60, .vx = +0.38, .vy = +0.61, .vz = -0.75 },
 * 	{ .x = -0.48, .y = -0.78, .z = +0.83, .vx = -0.58, .vy = -0.09, .vz = +0.44 },
 * 	{ .x = +0.99, .y = -0.13, .z = -0.48, .vx = +0.44, .vy = +0.06, .vz = -0.58 },
 * 	{ .x = -0.26, .y = +0.92, .z = -0.98, .vx = -0.75, .vy = +0.45, .vz = +0.30 },
 * 	{ .x = +0.88, .y = -0.50, .z = +0.53, .vx = -0.82, .vy = +0.50, .vz = -0.36 },
 * 	{ .x = -0.20, .y = +0.55, .z = +0.85, .vx = -0.87, .vy = +0.41, .vz = +0.38 },
 * 	{ .x = +0.48, .y = -0.09, .z = +0.29, .vx = -0.76, .vy = +0.27, .vz = -0.38 },
 * 	{ .x = +0.11, .y = -0.13, .z = -0.02, .vx = +0.69, .vy = +0.64, .vz = +0.06 },
 * 	{ .x = -0.85, .y = +0.08, .z = -0.25, .vx = -0.20, .vy = -0.69, .vz = +0.77 }
 * };
 * size_t nrows = sizeof(data)/sizeof(*data);
 * size_t ncols = 6;
 * size_t types[6] = {LAL_D_TYPE_CODE, LAL_D_TYPE_CODE, LAL_D_TYPE_CODE, LAL_D_TYPE_CODE, LAL_D_TYPE_CODE, LAL_D_TYPE_CODE};
 * size_t offsets[6] = {offsetof(struct body, x), offsetof(struct body, y), offsetof(struct body, z), offsetof(struct body, vx), offsetof(struct body, vy), offsetof(struct body, vz)};
 * size_t colsz[6] = {sizeof(data->x), sizeof(data->y), sizeof(data->z), sizeof(data->vx), sizeof(data->vy), sizeof(data->vz)};
 * size_t rowsz = sizeof(*data);
 * LALH5File *file = XLALH5FileOpen("example.h5", "w");
 * LALH5Dataset *dset = XLALH5TableAlloc(file, "particles", ncols, cols, types, offsets, rowsz);
 * XLALH5TableAppend(dset, offsets, colsz, nrows, rowsz, data);
 * @endcode
 *
 * @param dset Pointer to a #LALH5Dataset containing the table into which the
 * rows will be appended.
 * @param offsets Pointer to an array of offsets for each column.
 * @param colsz Pointer to an array of sizes of each column.
 * @param nrows Number of rows of data that will be appended.
 * @param rowsz Size of each row of data.
 * @param data Pointer to a memory in which contains the data to be appended.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5TableAppend(LALH5Dataset UNUSED *dset, const size_t UNUSED *offsets, const size_t UNUSED *colsz, size_t UNUSED nrows, size_t UNUSED rowsz, const void UNUSED *data)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (dset == NULL || offsets == NULL || colsz == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5TBappend_records(dset->parent_id, dset->name, nrows, rowsz, offsets, colsz, data) < 0)
		XLAL_ERROR(XLAL_EIO);

	return 0;
#endif
}

/**
 * @brief Reads table data from a #LALH5Dataset
 * @details
 * This routine reads the table data from a HDF5 dataset associated with the
 * #LALH5Dataset @p dset and stores the data in the buffer @p data.
 * This buffer should be sufficiently large to hold the entire contents
 * of the dataset, which is comprised of a number of rows each of size
 * @p rowsz.  The number of rows can be determined with the routine
 * XLALH5TableQueryNRows().  Each row is comprised of a number of columns,
 * which can be determined with the routine XLALH5TableQueryNColumns().
 * The offsets of the columns and the sizes of the columns in each row
 * of data is provided by the arrays @p offsets and @p colsz respectively.
 *
 * The following example shows how to read the data from a table
 * named @a particles.
 * @code
 * #include <stddef.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/LALH5FileIO.h>
 *
 * struct body {REAL8 x; REAL8 y; REAL8 z; REAL8 vx; REAL8 vy; REAL8 vz;} *data;
 * size_t offsets[] = {offsetof(struct body, x), offsetof(struct body, y), offsetof(struct body, z), offsetof(struct body, vx), offsetof(struct body, vy), offsetof(struct body, vz)};
 * size_t colsz[] = {sizeof(data->x), sizeof(data->y), sizeof(data->z), sizeof(data->vx), sizeof(data->vy), sizeof(data->vz)};
 * size_t rowsz = sizeof(*data);
 * LALH5File *file = XLALH5FileOpen("example.h5", "r");
 * LALH5Dataset *dset = XLALH5DatasetRead(file, "particles");
 * size_t nrows = XLALH5TableQueryNRows(dset);
 * data = LALMalloc(nrows * rowsz);
 * XLALH5TableRead(data, dset, offsets, colsz, rowsz);
 * @endcode
 *
 * @param data Pointer to a memory in which to store the data.
 * @param dset Pointer to a #LALH5Dataset from which to extract the data.
 * @param offsets Pointer to an array of offsets for each column.
 * @param colsz Pointer to an array of sizes of each column.
 * @param rowsz Size of each row of data.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5TableRead(void UNUSED *data, const LALH5Dataset UNUSED *dset, const size_t UNUSED *offsets, const size_t UNUSED *colsz, size_t UNUSED rowsz)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (data == NULL || dset == NULL || offsets == NULL || colsz == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5TBread_table(dset->parent_id, dset->name, rowsz, offsets, colsz, data) < 0)
		XLAL_ERROR(XLAL_EIO);

	return 0;
#endif
}

/**
 * @brief Reads certain rows of table data from a #LALH5Dataset
 * @details
 * This routine reads certain rows of table data from a HDF5 dataset associated
 * with the #LALH5Dataset @p dset and stores the data in the buffer @p data.
 * This buffer should be sufficiently large to hold the requested @p nrows rows
 * of the dataset, starting with row @p row0 (the first row is row 0), which is
 * comprised of a number of rows each of size @p rowsz.  Each row is comprised
 * of a number of columns, which can be determined with the routine
 * XLALH5TableQueryNColumns().  The offsets of the columns and the sizes of the
 * columns in each row of data is provided by the arrays @p offsets and
 * @p colsz respectively.
 *
 * The following example shows how to read rows 2 to 5 inclusive
 * (where the first row is row 0) of data from a table named
 * @a particles.
 * @code
 * #include <stddef.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/LALH5FileIO.h>
 *
 * struct body {REAL8 x; REAL8 y; REAL8 z; REAL8 vx; REAL8 vy; REAL8 vz;} data[4];
 * size_t offsets[] = {offsetof(struct body, x), offsetof(struct body, y), offsetof(struct body, z), offsetof(struct body, vx), offsetof(struct body, vy), offsetof(struct body, vz)};
 * size_t colsz[] = {sizeof(data->x), sizeof(data->y), sizeof(data->z), sizeof(data->vx), sizeof(data->vy), sizeof(data->vz)};
 * size_t rowsz = sizeof(*data);
 * size_t row0 = 2;
 * size_t nrows = sizeof(data)/sizeof(*data);
 * LALH5File *file = XLALH5FileOpen("example.h5", "r");
 * LALH5Dataset *dset = XLALH5DatasetRead(file, "particles");
 * XLALH5TableReadRows(data, dset, offsets, colsz, row0, nrows, rowsz);
 * @endcode
 *
 * @param data Pointer to a memory in which to store the data.
 * @param dset Pointer to a #LALH5Dataset from which to extract the data.
 * @param offsets Pointer to an array of offsets for each column.
 * @param colsz Pointer to an array of sizes of each column.
 * @param row0 The first row to read.
 * @param nrows The Number of rows to read.
 * @param rowsz Size of each row of data.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5TableReadRows(void UNUSED *data, const LALH5Dataset UNUSED *dset, const size_t UNUSED *offsets, const size_t UNUSED *colsz, size_t UNUSED row0, size_t UNUSED nrows, size_t UNUSED rowsz)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (data == NULL || dset == NULL || offsets == NULL || colsz == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5TBread_records(dset->parent_id, dset->name, row0, nrows, rowsz, offsets, colsz, data) < 0)
		XLAL_ERROR(XLAL_EIO);

	return 0;
#endif
}

/**
 * @brief Reads certain columns of table data from a #LALH5Dataset
 * @details
 * This routine reads certain columns of table data from a HDF5 dataset associated
 * with the #LALH5Dataset @p dset and stores the data in the buffer @p data.
 * This buffer should be sufficiently large to hold the requested rows
 * of the dataset, which is comprised of a number of rows each of size
 * @p rowsz.  Each row is comprised of a number of columns
 * which are named in the comma-separated list in string @p cols.
 * The offsets of the columns and the sizes of the columns in each row
 * of data is provided by the arrays @p offsets and @p colsz respectively.
 *
 * The following example shows how to read rows 2 to 5 inclusive
 * (where the first row is row 0) of certain columns of data from a table named
 * @a particles.
 * @code
 * #include <stddef.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/LALH5FileIO.h>
 *
 * struct position {REAL8 x; REAL8 y; REAL8 z;} data[4];
 * const char *cols = "x,y,z";
 * size_t offsets[] = {offsetof(struct position, x), offsetof(struct position, y), offsetof(struct position, z)};
 * size_t colsz[] = {sizeof(data->x), sizeof(data->y), sizeof(data->z)};
 * size_t rowsz = sizeof(*data);
 * size_t row0 = 2;
 * size_t nrows = sizeof(data)/sizeof(*data);
 * LALH5File *file = XLALH5FileOpen("example.h5", "r");
 * LALH5Dataset *dset = XLALH5DatasetRead(file, "particles");
 * XLALH5TableReadColumns(data, dset, cols, offsets, colsz, row0, nrows, rowsz);
 * @endcode
 *
 * @param data Pointer to a memory in which to store the data.
 * @param dset Pointer to a #LALH5Dataset from which to extract the data.
 * @param cols Pointer to an string listing the column names separated by commas.
 * @param offsets Pointer to an array of offsets for each column.
 * @param colsz Pointer to an array of sizes of each column.
 * @param row0 The first row to read.
 * @param nrows The Number of rows to read.
 * @param rowsz Size of each row of data.
 * @retval 0 Success.
 * @retval -1 Failure.
 */
int XLALH5TableReadColumns(void UNUSED *data, const LALH5Dataset UNUSED *dset, const char UNUSED *cols, const size_t UNUSED *offsets, const size_t UNUSED *colsz, size_t UNUSED row0, size_t UNUSED nrows, size_t UNUSED rowsz)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	if (data == NULL || dset == NULL || offsets == NULL || colsz == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5TBread_fields_name(dset->parent_id, dset->name, cols, row0, nrows, rowsz, offsets, colsz, data) < 0)
		XLAL_ERROR(XLAL_EIO);

	return 0;
#endif
}

/**
 * @brief Gets the number of rows in a #LALH5Dataset containing table data.
 * @param dset Pointer to a #LALH5Dataset containing table data to be queried.
 * @returns The number of rows of table data in the HDF5 dataset associated
 * with the specified #LALH5Dataset.
 * @retval (size_t)(-1) Failure.
 */
size_t XLALH5TableQueryNRows(const LALH5Dataset UNUSED *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hsize_t nrows;

	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5TBget_table_info(dset->parent_id, dset->name, NULL, &nrows) < 0)
		XLAL_ERROR(XLAL_EIO);

	return nrows;
#endif
}

/**
 * @brief Gets the number of columns in a #LALH5Dataset containing table data.
 * @param dset Pointer to a #LALH5Dataset containing table data to be queried.
 * @returns The number of columns of table data in the HDF5 dataset associated
 * with the specified #LALH5Dataset.
 * @retval (size_t)(-1) Failure.
 */
size_t XLALH5TableQueryNColumns(const LALH5Dataset UNUSED *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hsize_t ncols;

	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5TBget_table_info(dset->parent_id, dset->name, &ncols, NULL) < 0)
		XLAL_ERROR(XLAL_EIO);

	return ncols;
#endif
}

/**
 * @brief Gets the number of bytes in a row in a #LALH5Dataset containing table
 * data.
 * @param dset Pointer to a #LALH5Dataset containing table data to be queried.
 * @returns The number of bytes of a row of data in the HDF5 dataset associated
 * with the specified #LALH5Dataset.
 * @retval (size_t)(-1) Failure.
 */
size_t XLALH5TableQueryRowSize(const LALH5Dataset UNUSED *dset)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	size_t rowsz;

	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	if (H5TBget_field_info(dset->parent_id, dset->name, NULL, NULL, NULL, &rowsz) < 0)
		XLAL_ERROR(XLAL_EIO);

	return rowsz;
#endif
}

/**
 * @brief Gets the name of a column in a #LALH5Dataset containing table
 * data.
 * @details
 * This routines gets the name of a column of data in a #LALH5Dataset
 * @p dset that contains table data.
 * The index @p pos identifies which column's name is returned.
 * The result is written into the buffer pointed to by @p name, the size
 * of which is @p size bytes.  If @p name is NULL, no data is copied but
 * the routine returns the length of the string.  Therefore, this routine
 * can be called once to determine the amount of memory required, the
 * memory can be allocated, and then it can be called a second time to
 * read the string.  If the parameter @p size is less than or equal to
 * the string length then only $p size-1 bytes of the string are copied
 * to the buffer @p name.
 * @note The return value is the length of the string, not including the
 * terminating NUL character; thus the buffer @p name should be allocated
 * to be one byte larger.
 * @param name Pointer to a buffer into which the string will be written.
 * @param size Size in bytes of the buffer into which the string will be
 * written.
 * @param dset Pointer to a #LALH5Dataset containing table data to be
 * queried.
 * @param pos The index identifying which column of the table.
 * @retval  0 Success.
 * @retval -1 Failure.
 */
int XLALH5TableQueryColumnName(char UNUSED *name, size_t UNUSED size, const LALH5Dataset UNUSED *dset, int UNUSED pos)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	int ncols;
	char *col;
	int n;

	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	ncols = H5Tget_nmembers(dset->dtype_id);
	if (ncols < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read type members");

	if (ncols <= pos)
		XLAL_ERROR(XLAL_EINVAL, "Requested column does not exist");

	col = H5Tget_member_name(dset->dtype_id, pos);
	if (col == NULL)
		XLAL_ERROR(XLAL_EIO, "Could not read column name");

	n = snprintf(name, name == NULL ? 0 : size, "%s", col);
	free(col);
	return n;
#endif
}

/**
 * @brief Gets the size in bytes of a column in a #LALH5Dataset containing table
 * data.
 * @param dset Pointer to a #LALH5Dataset containing table data to be
 * queried.
 * @param pos The index identifying which column of the table.
 * @returns The size in bytes of the specified column entry in the HDF5 dataset
 * associated with the specified #LALH5Dataset.
 * @retval -1 Failure.
 */
size_t XLALH5TableQueryColumnSize(const LALH5Dataset UNUSED *dset, int UNUSED pos)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t member_type_id;
	hid_t native_member_type_id;
	int ncols;
	size_t size;

	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	ncols = H5Tget_nmembers(dset->dtype_id);
	if (ncols < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read type members");

	if (ncols <= pos)
		XLAL_ERROR(XLAL_EINVAL, "Requested column does not exist");

	member_type_id = H5Tget_member_type(dset->dtype_id, pos);
	if (member_type_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read column type");

	native_member_type_id = H5Tget_native_type(member_type_id, H5T_DIR_DEFAULT);
	if (native_member_type_id < 0) {
		H5Tclose(member_type_id);
		XLAL_ERROR(XLAL_EIO, "Could not read column native type");
	}

	size = H5Tget_size(native_member_type_id);
	if (size == 0) {
		H5Tclose(native_member_type_id);
		H5Tclose(member_type_id);
		XLAL_ERROR(XLAL_EIO, "Could not read column size");
	}

	H5Tclose(native_member_type_id);
	H5Tclose(member_type_id);

	return size;
#endif
}

/**
 * @brief Gets the type of data stored in a column in a #LALH5Dataset
 * containing table data.
 * @param dset Pointer to a #LALH5Dataset containing table data to be
 * queried.
 * @param pos The index identifying which column of the table.
 * @returns The #LALTYPECODE of the datatype of the specified column entry in
 * the HDF5 dataset associated with the specified #LALH5Dataset.
 * @retval -1 Failure.
 */
LALTYPECODE XLALH5TableQueryColumnType(const LALH5Dataset UNUSED *dset, int UNUSED pos)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	hid_t member_type_id;
	int ncols;
	LALTYPECODE type;

	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	ncols = H5Tget_nmembers(dset->dtype_id);
	if (ncols < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read type members");

	if (ncols <= pos)
		XLAL_ERROR(XLAL_EINVAL, "Requested column does not exist");

	member_type_id = H5Tget_member_type(dset->dtype_id, pos);
	if (member_type_id < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read column type");

	type = XLALTypeFromH5Type(member_type_id);
	H5Tclose(member_type_id);

	if ((int)type < 0)
		XLAL_ERROR(XLAL_EFUNC);
	return type;
#endif
}

/**
 * @brief Gets the offset of the data in a column in a #LALH5Dataset
 * containing table data.
 * @param dset Pointer to a #LALH5Dataset containing table data to be
 * queried.
 * @param pos The index identifying which column of the table.
 * @returns The offset of the specified column data in a row of data
 * in the HDF5 dataset associated with the specified #LALH5Dataset.
 * @retval (size_t)(-1) Failure.
 */
size_t XLALH5TableQueryColumnOffset(const LALH5Dataset UNUSED *dset, int UNUSED pos)
{
#ifndef HAVE_HDF5
	XLAL_ERROR(XLAL_EFAILED, "HDF5 support not implemented");
#else
	int ncols;
	size_t offset;

	if (dset == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	ncols = H5Tget_nmembers(dset->dtype_id);
	if (ncols < 0)
		XLAL_ERROR(XLAL_EIO, "Could not read type members");

	if (ncols <= pos)
		XLAL_ERROR(XLAL_EINVAL, "Requested column does not exist");

	/* H5Tget_member_offset fails iff H5Tget_member_class does */
	if (H5Tget_member_class(dset->dtype_id, pos) < 0)
		XLAL_ERROR(XLAL_EINVAL, "Could not read offset");

	offset = H5Tget_member_offset(dset->dtype_id, pos);
	return offset;
#endif
}

/** @} */
