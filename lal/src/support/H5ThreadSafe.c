#ifdef HAVE_HDF5

#if defined LAL_PTHREAD_LOCK && ! defined H5_HAVE_THREADSAFE

/* LAL is threadsafe but the HDF5 library is not, so make it threadsafe */

#include <pthread.h>
static pthread_mutex_t lalHDF5Mutex = PTHREAD_MUTEX_INITIALIZER;
#define LAL_HDF5_MUTEX_LOCK pthread_mutex_lock(&lalHDF5Mutex);
#define LAL_HDF5_MUTEX_UNLOCK pthread_mutex_unlock(&lalHDF5Mutex);

/* wrap HDF5 library routines with the mutex */

static inline herr_t threadsafe_H5Aclose(hid_t attr_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Aclose(attr_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Acreate2(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id, hid_t aapl_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Acreate2(loc_id, attr_name, type_id, space_id, acpl_id, aapl_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline ssize_t threadsafe_H5Aget_name(hid_t attr_id, size_t buf_size, char *buf)
{
	LAL_HDF5_MUTEX_LOCK
	ssize_t retval = H5Aget_name(attr_id, buf_size, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline int threadsafe_H5Aget_num_attrs(hid_t loc_id)
{
	LAL_HDF5_MUTEX_LOCK
	int retval = H5Aget_num_attrs(loc_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Aget_space(hid_t attr_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Aget_space(attr_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Aget_type(hid_t attr_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Aget_type(attr_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Aopen(hid_t obj_id, const char *attr_name, hid_t aapl_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Aopen(obj_id, attr_name, aapl_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Aopen_idx(hid_t loc_id, unsigned idx)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Aopen_idx(loc_id, idx);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Aread(hid_t attr_id, hid_t type_id, void *buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Aread(attr_id, type_id, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Awrite(hid_t attr_id, hid_t type_id, const void *buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Awrite(attr_id, type_id, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Dclose(hid_t dset_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Dclose(dset_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Dcreate2(hid_t loc_id, const char *name, hid_t type_id, hid_t space_id, hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Dcreate2(loc_id, name, type_id, space_id, lcpl_id, dcpl_id, dapl_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Dget_space(hid_t dset_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Dget_space(dset_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Dget_type(hid_t dset_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Dget_type(dset_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Dopen2(hid_t file_id, const char *name, hid_t dapl_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Dopen2(file_id, name, dapl_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Dread(hid_t dset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t plist_id, void *buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Dread(dset_id, mem_type_id, mem_space_id, file_space_id, plist_id, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Dvlen_reclaim(hid_t type_id, hid_t space_id, hid_t plist_id, void *buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Dvlen_reclaim(type_id, space_id, plist_id, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Dwrite(hid_t dset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t plist_id, const void *buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Dwrite(dset_id, mem_type_id, mem_space_id, file_space_id, plist_id, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Fclose(hid_t file_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Fclose(file_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Fcreate(const char *filename, unsigned flags, hid_t create_plist, hid_t access_plist)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Fcreate(filename, flags, create_plist, access_plist);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Fflush(hid_t object_id, H5F_scope_t scope)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Fflush(object_id, scope);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline ssize_t threadsafe_H5Fget_name(hid_t obj_id, char *name, size_t size)
{
	LAL_HDF5_MUTEX_LOCK
	ssize_t retval = H5Fget_name(obj_id, name, size);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Fopen(const char *filename, unsigned flags, hid_t access_plist)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Fopen(filename, flags, access_plist);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Gclose(hid_t group_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Gclose(group_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Gcreate2(hid_t loc_id, const char *name, hid_t lcpl_id, hid_t gcpl_id, hid_t gapl_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Gcreate2(loc_id, name, lcpl_id, gcpl_id, gapl_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Gget_info(hid_t loc_id, H5G_info_t *ginfo)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Gget_info(loc_id, ginfo);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Gget_info_by_name(hid_t loc_id, const char *name, H5G_info_t *ginfo, hid_t lapl_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Gget_info_by_name(loc_id, name, ginfo, lapl_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Gget_num_objs(hid_t loc_id, hsize_t *num_objs)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Gget_num_objs(loc_id, num_objs);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Gget_objinfo(hid_t loc_id, const char *name, hbool_t follow_link, H5G_stat_t *statbuf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Gget_objinfo(loc_id, name, follow_link, statbuf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline ssize_t threadsafe_H5Gget_objname_by_idx(hid_t loc_id, hsize_t idx, char* name, size_t size)
{
	LAL_HDF5_MUTEX_LOCK
	ssize_t retval = H5Gget_objname_by_idx(loc_id, idx, name, size);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline H5G_obj_t threadsafe_H5Gget_objtype_by_idx(hid_t loc_id, hsize_t idx)
{
	LAL_HDF5_MUTEX_LOCK
	H5G_obj_t retval = H5Gget_objtype_by_idx(loc_id, idx);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Gopen2(hid_t loc_id, const char *name, hid_t gapl_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Gopen2(loc_id, name, gapl_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline ssize_t threadsafe_H5Iget_name(hid_t id, char *name , size_t size)
{
	LAL_HDF5_MUTEX_LOCK
	ssize_t retval = H5Iget_name(id, name , size);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Oclose(hid_t object_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Oclose(object_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Oget_info(hid_t loc_id, H5O_info_t *oinfo)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Oget_info(loc_id, oinfo);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Oget_info_by_idx(hid_t loc_id, const char *group_name, H5_index_t idx_type, H5_iter_order_t order, hsize_t n, H5O_info_t *oinfo, hid_t lapl_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Oget_info_by_idx(loc_id, group_name, idx_type, order, n, oinfo, lapl_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Oopen_by_addr(hid_t loc_id, haddr_t addr)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Oopen_by_addr(loc_id, addr);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Pclose(hid_t plist_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Pclose(plist_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Pcreate(hid_t cls_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Pcreate(cls_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Pset_create_intermediate_group(hid_t plist_id, unsigned crt_intmd)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Pset_create_intermediate_group(plist_id, crt_intmd);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Sclose(hid_t space_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Sclose(space_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Screate(H5S_class_t type)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Screate(type);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Screate_simple(int rank, const hsize_t dims[], const hsize_t maxdims[])
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Screate_simple(rank, dims, maxdims);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline int threadsafe_H5Sget_simple_extent_dims(hid_t space_id, hsize_t dims[], hsize_t maxdims[])
{
	LAL_HDF5_MUTEX_LOCK
	int retval = H5Sget_simple_extent_dims(space_id, dims, maxdims);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline int threadsafe_H5Sget_simple_extent_ndims(hid_t space_id)
{
	LAL_HDF5_MUTEX_LOCK
	int retval = H5Sget_simple_extent_ndims(space_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hssize_t threadsafe_H5Sget_simple_extent_npoints(hid_t space_id)
{
	LAL_HDF5_MUTEX_LOCK
	hssize_t retval = H5Sget_simple_extent_npoints(space_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5TBappend_records(hid_t loc_id, const char *dset_name, hsize_t nrecords, size_t type_size, const size_t *field_offset, const size_t *dst_sizes, const void *buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5TBappend_records(loc_id, dset_name, nrecords, type_size, field_offset, dst_sizes, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5TBget_field_info(hid_t loc_id, const char *dset_name, char *field_names[], size_t *field_sizes, size_t *field_offsets, size_t *type_size)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5TBget_field_info(loc_id, dset_name, field_names, field_sizes, field_offsets, type_size);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5TBget_table_info(hid_t loc_id, const char *dset_name, hsize_t *nfields, hsize_t *nrecords)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5TBget_table_info(loc_id, dset_name, nfields, nrecords); 
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5TBmake_table(const char *table_title, hid_t loc_id, const char *dset_name, hsize_t nfields, hsize_t nrecords, size_t type_size, const char *field_names[], const size_t *field_offset, const hid_t *field_types, hsize_t chunk_size, void *fill_data, int compress, const void *buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5TBmake_table(table_title, loc_id, dset_name, nfields, nrecords, type_size, field_names, field_offset, field_types, chunk_size, fill_data, compress, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5TBread_fields_name(hid_t loc_id, const char *dset_name, const char *field_names, hsize_t start, hsize_t nrecords, size_t type_size, const size_t *field_offset, const size_t *dst_sizes, void *buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5TBread_fields_name(loc_id, dset_name, field_names, start, nrecords, type_size, field_offset, dst_sizes, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5TBread_records(hid_t loc_id, const char *dset_name, hsize_t start, hsize_t nrecords, size_t type_size, const size_t *dst_offset, const size_t *dst_sizes, void *buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5TBread_records(loc_id, dset_name, start, nrecords, type_size, dst_offset, dst_sizes, buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5TBread_table(hid_t loc_id, const char *dset_name, size_t dst_size, const size_t *dst_offset, const size_t *dst_sizes, void *dst_buf)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5TBread_table(loc_id, dset_name, dst_size, dst_offset, dst_sizes, dst_buf);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Tarray_create2(hid_t base_id, unsigned ndims, const hsize_t dim[])
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Tarray_create2(base_id, ndims, dim);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Tclose(hid_t type_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Tclose(type_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Tcopy(hid_t type_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Tcopy(type_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Tcreate(H5T_class_t type, size_t size)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Tcreate(type, size);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Tenum_create(hid_t base_id)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Tenum_create(base_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Tenum_insert(hid_t type, const char *name, const void *value)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Tenum_insert(type, name, value);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline int threadsafe_H5Tget_array_dims2(hid_t type_id, hsize_t dims[])
{
	LAL_HDF5_MUTEX_LOCK
	int retval = H5Tget_array_dims2(type_id, dims);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline int threadsafe_H5Tget_array_ndims(hid_t type_id)
{
	LAL_HDF5_MUTEX_LOCK
	int retval = H5Tget_array_ndims(type_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline H5T_class_t threadsafe_H5Tget_class(hid_t type_id)
{
	LAL_HDF5_MUTEX_LOCK
	H5T_class_t retval = H5Tget_class(type_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline H5T_class_t threadsafe_H5Tget_member_class(hid_t type_id, unsigned membno)
{
	LAL_HDF5_MUTEX_LOCK
	H5T_class_t retval = H5Tget_member_class(type_id, membno);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline char *threadsafe_H5Tget_member_name(hid_t type_id, unsigned membno)
{
	LAL_HDF5_MUTEX_LOCK
	char *retval = H5Tget_member_name(type_id, membno);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline size_t threadsafe_H5Tget_member_offset(hid_t type_id, unsigned membno)
{
	LAL_HDF5_MUTEX_LOCK
	size_t retval = H5Tget_member_offset(type_id, membno);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Tget_member_type(hid_t type_id, unsigned membno)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Tget_member_type(type_id, membno);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Tget_member_value(hid_t type_id, unsigned membno, void *value)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Tget_member_value(type_id, membno, value);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Tget_native_type(hid_t type_id, H5T_direction_t direction)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Tget_native_type(type_id, direction);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline int threadsafe_H5Tget_nmembers(hid_t type_id)
{
	LAL_HDF5_MUTEX_LOCK
	int retval = H5Tget_nmembers(type_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline H5T_sign_t threadsafe_H5Tget_sign(hid_t type_id)
{
	LAL_HDF5_MUTEX_LOCK
	H5T_sign_t retval = H5Tget_sign(type_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline size_t threadsafe_H5Tget_size(hid_t type_id)
{
	LAL_HDF5_MUTEX_LOCK
	size_t retval = H5Tget_size(type_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline hid_t threadsafe_H5Tget_super(hid_t type)
{
	LAL_HDF5_MUTEX_LOCK
	hid_t retval = H5Tget_super(type);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Tinsert(hid_t parent_id, const char *name, size_t offset, hid_t member_id)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Tinsert(parent_id, name, offset, member_id);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5Tset_size(hid_t type_id, size_t size)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5Tset_size(type_id, size);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5check_version(unsigned majnum, unsigned minnum, unsigned relnum)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5check_version(majnum, minnum, relnum);
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

static inline herr_t threadsafe_H5open(void)
{
	LAL_HDF5_MUTEX_LOCK
	herr_t retval = H5open();
	LAL_HDF5_MUTEX_UNLOCK
	return retval;
}

#undef LAL_HDF5_MUTEX_LOCK
#undef LAL_HDF5_MUTEX_UNLOCK

#else /* either the HDF5 library is threadsafe, or LAL does not require */

/* use regular HDF library routines */

#define threadsafe_H5Aclose H5Aclose
#define threadsafe_H5Acreate2 H5Acreate2
#define threadsafe_H5Aget_name H5Aget_name
#define threadsafe_H5Aget_num_attrs H5Aget_num_attrs
#define threadsafe_H5Aget_space H5Aget_space
#define threadsafe_H5Aget_type H5Aget_type
#define threadsafe_H5Aopen H5Aopen
#define threadsafe_H5Aopen_idx H5Aopen_idx
#define threadsafe_H5Aread H5Aread
#define threadsafe_H5Awrite H5Awrite
#define threadsafe_H5Dclose H5Dclose
#define threadsafe_H5Dcreate2 H5Dcreate2
#define threadsafe_H5Dget_space H5Dget_space
#define threadsafe_H5Dget_type H5Dget_type
#define threadsafe_H5Dopen2 H5Dopen2
#define threadsafe_H5Dread H5Dread
#define threadsafe_H5Dvlen_reclaim H5Dvlen_reclaim
#define threadsafe_H5Dwrite H5Dwrite
#define threadsafe_H5Fclose H5Fclose
#define threadsafe_H5Fcreate H5Fcreate
#define threadsafe_H5Fflush H5Fflush
#define threadsafe_H5Fget_name H5Fget_name
#define threadsafe_H5Fopen H5Fopen
#define threadsafe_H5Gclose H5Gclose
#define threadsafe_H5Gcreate2 H5Gcreate2
#define threadsafe_H5Gget_info H5Gget_info
#define threadsafe_H5Gget_info_by_name H5Gget_info_by_name
#define threadsafe_H5Gget_num_objs H5Gget_num_objs
#define threadsafe_H5Gget_objinfo H5Gget_objinfo
#define threadsafe_H5Gget_objname_by_idx H5Gget_objname_by_idx
#define threadsafe_H5Gget_objtype_by_idx H5Gget_objtype_by_idx
#define threadsafe_H5Gopen2 H5Gopen2
#define threadsafe_H5Iget_name H5Iget_name
#define threadsafe_H5Oclose H5Oclose
#define threadsafe_H5Oget_info H5Oget_info
#define threadsafe_H5Oget_info_by_idx H5Oget_info_by_idx
#define threadsafe_H5Oopen_by_addr H5Oopen_by_addr
#define threadsafe_H5Pclose H5Pclose
#define threadsafe_H5Pcreate H5Pcreate
#define threadsafe_H5Pset_create_intermediate_group H5Pset_create_intermediate_group
#define threadsafe_H5Sclose H5Sclose
#define threadsafe_H5Screate H5Screate
#define threadsafe_H5Screate_simple H5Screate_simple
#define threadsafe_H5Sget_simple_extent_dims H5Sget_simple_extent_dims
#define threadsafe_H5Sget_simple_extent_ndims H5Sget_simple_extent_ndims
#define threadsafe_H5Sget_simple_extent_npoints H5Sget_simple_extent_npoints
#define threadsafe_H5TBappend_records H5TBappend_records
#define threadsafe_H5TBget_field_info H5TBget_field_info
#define threadsafe_H5TBget_table_info H5TBget_table_info
#define threadsafe_H5TBmake_table H5TBmake_table
#define threadsafe_H5TBread_fields_name H5TBread_fields_name
#define threadsafe_H5TBread_records H5TBread_records
#define threadsafe_H5TBread_table H5TBread_table
#define threadsafe_H5Tarray_create2 H5Tarray_create2
#define threadsafe_H5Tclose H5Tclose
#define threadsafe_H5Tcopy H5Tcopy
#define threadsafe_H5Tcreate H5Tcreate
#define threadsafe_H5Tenum_create H5Tenum_create
#define threadsafe_H5Tenum_insert H5Tenum_insert
#define threadsafe_H5Tget_array_dims2 H5Tget_array_dims2
#define threadsafe_H5Tget_array_ndims H5Tget_array_ndims
#define threadsafe_H5Tget_class H5Tget_class
#define threadsafe_H5Tget_member_class H5Tget_member_class
#define threadsafe_H5Tget_member_name H5Tget_member_name
#define threadsafe_H5Tget_member_offset H5Tget_member_offset
#define threadsafe_H5Tget_member_type H5Tget_member_type
#define threadsafe_H5Tget_member_value H5Tget_member_value
#define threadsafe_H5Tget_native_type H5Tget_native_type
#define threadsafe_H5Tget_nmembers H5Tget_nmembers
#define threadsafe_H5Tget_sign H5Tget_sign
#define threadsafe_H5Tget_size H5Tget_size
#define threadsafe_H5Tget_super H5Tget_super
#define threadsafe_H5Tinsert H5Tinsert
#define threadsafe_H5Tset_size H5Tset_size
#define threadsafe_H5check_version H5check_version
#define threadsafe_H5open H5open

#endif

#endif /* HAVE_HDF5 */
