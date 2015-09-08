#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)

#define VTYPE CONCAT2(TYPE,Vector)
#define TCODE CONCAT3(LAL_,TYPECODE,_TYPE_CODE)

#define ALLOCFUNC CONCAT2(XLALH5DatasetAlloc,VTYPE)
#define READFUNC CONCAT2(XLALH5DatasetRead,VTYPE)

#define CREATEFUNC CONCAT2(XLALCreate,VTYPE)
#define DESTROYFUNC CONCAT2(XLALDestroy,VTYPE)

LALH5Dataset *ALLOCFUNC(LALH5File *file, const char *name, VTYPE *vector)
{
	LALH5Dataset *dataset;
	if (!file || !name || !vector)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (!vector->length || !vector->data)
		XLAL_ERROR_NULL(XLAL_EINVAL);
	dataset = XLALH5DatasetAlloc1D(file, name, TCODE, vector->length);
	if (!dataset)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	if (XLALH5DatasetWrite(dataset, vector->data) < 0) {
		XLALH5DatasetFree(dataset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	return dataset;
}

VTYPE *READFUNC(LALH5Dataset *dset)
{
	VTYPE *vector;
	LALTYPECODE type;
	size_t npoints;
	int ndim;

	/* error checking */

	if (!dset)
		XLAL_ERROR_NULL(XLAL_EFAULT);

	ndim = XLALH5DatasetQueryNDim(dset);
	if (ndim != 1)
		XLAL_ERROR_NULL(XLAL_EDIMS);

	type = XLALH5DatasetQueryType(dset);
	if (type != TCODE)
		XLAL_ERROR_NULL(XLAL_ETYPE);

	npoints = XLALH5DatasetQueryNPoints(dset);
	if (npoints == (size_t)(-1))
		XLAL_ERROR_NULL(XLAL_EFUNC);

	vector = CREATEFUNC(npoints);
	if (!vector)
		XLAL_ERROR_NULL(XLAL_ENOMEM);

	if (XLALH5DatasetQueryData(vector->data, dset) == -1) {
		DESTROYFUNC(vector);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	return vector;
}

#undef CONCAT2x
#undef CONCAT2
#undef CONCAT3x
#undef CONCAT3

#undef VTYPE
#undef TCODE

#undef ALLOCFUNC
#undef READFUNC

#undef CREATEFUNC
#undef DESTROYFUNC
