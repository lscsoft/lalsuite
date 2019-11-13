#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)

#define VTYPE CONCAT2(TYPE,Vector)

#define WRITEFUNC CONCAT2(XLALH5FileWrite,VTYPE)
#define READFUNC CONCAT2(XLALH5FileRead,VTYPE)
#define DSETALLOCFUNC CONCAT2(XLALH5DatasetAlloc,VTYPE)
#define DSETREADFUNC CONCAT2(XLALH5DatasetRead,VTYPE)

int WRITEFUNC(LALH5File *file, const char *name, VTYPE *vector)
{
	LALH5Dataset *dataset;
	if (!file || !name || !vector)
		XLAL_ERROR(XLAL_EFAULT);
	if (!vector->length || !vector->data)
		XLAL_ERROR(XLAL_EINVAL);
	dataset = DSETALLOCFUNC(file, name, vector);
	if (!dataset)
		XLAL_ERROR(XLAL_EFUNC);
	XLALH5DatasetFree(dataset);
	return 0;
}

VTYPE *READFUNC(LALH5File *file, const char *dset)
{
	VTYPE *vector;
	LALH5Dataset *dataset;
	if (!file || !dset)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	dataset = XLALH5DatasetRead(file, dset);
	if (!dataset)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	vector = DSETREADFUNC(dataset);
	XLALH5DatasetFree(dataset);
	if (!vector)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	return vector;
}

#undef WRITEFUNC
#undef READFUNC
#undef DSETALLOCFUNC
#undef DSETREADFUNC

#undef CONCAT2x
#undef CONCAT2

#undef VTYPE
