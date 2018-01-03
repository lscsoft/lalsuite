#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)

#define ATYPE CONCAT2(TYPE,Array)

#define WRITEFUNC CONCAT2(XLALH5FileWrite,ATYPE)
#define READFUNC CONCAT2(XLALH5FileRead,ATYPE)
#define DSETALLOCFUNC CONCAT2(XLALH5DatasetAlloc,ATYPE)
#define DSETREADFUNC CONCAT2(XLALH5DatasetRead,ATYPE)

int WRITEFUNC(LALH5File *file, const char *name, ATYPE *array)
{
	LALH5Dataset *dataset;
	if (!file || !name || !array)
		XLAL_ERROR(XLAL_EFAULT);
	if (!array->data || !array->dimLength || !array->dimLength->length || !array->dimLength->data)
		XLAL_ERROR(XLAL_EINVAL);
	dataset = DSETALLOCFUNC(file, name, array);
	if (!dataset)
		XLAL_ERROR(XLAL_EFUNC);
	XLALH5DatasetFree(dataset);
	return 0;
}

ATYPE *READFUNC(LALH5File *file, const char *dset)
{
	ATYPE *array;
	LALH5Dataset *dataset;
	if (!file || !dset)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	dataset = XLALH5DatasetRead(file, dset);
	if (!dataset)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	array = DSETREADFUNC(dataset);
	XLALH5DatasetFree(dataset);
	if (!array)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return array;
}

#undef WRITEFUNC
#undef READFUNC
#undef DSETALLOCFUNC
#undef DSETREADFUNC

#undef CONCAT2x
#undef CONCAT2

#undef ATYPE
