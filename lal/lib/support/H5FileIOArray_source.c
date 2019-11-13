#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)

#define ATYPE CONCAT2(TYPE,Array)
#define TCODE CONCAT3(LAL_,TYPECODE,_TYPE_CODE)

#define ALLOCFUNC CONCAT2(XLALH5DatasetAlloc,ATYPE)
#define READFUNC CONCAT2(XLALH5DatasetRead,ATYPE)

#define CREATEFUNC CONCAT2(XLALCreate,ATYPE)
#define DESTROYFUNC CONCAT2(XLALDestroy,ATYPE)

LALH5Dataset *ALLOCFUNC(LALH5File *file, const char *name, ATYPE *array)
{
	LALH5Dataset *dataset;
	if (!file || !name || !array)
		XLAL_ERROR_NULL(XLAL_EFAULT);
	if (!array->data || !array->dimLength || !array->dimLength->length || !array->dimLength->data)
		XLAL_ERROR_NULL(XLAL_EINVAL);
	dataset = XLALH5DatasetAlloc(file, name, TCODE, array->dimLength);
	if (!dataset)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	if (XLALH5DatasetWrite(dataset, array->data) < 0) {
		XLALH5DatasetFree(dataset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	return dataset;
}

ATYPE *READFUNC(LALH5Dataset *dset)
{
	ATYPE *array;
	LALTYPECODE type;
	UINT4Vector *dimLength;

	if (!dset)
		XLAL_ERROR_NULL(XLAL_EFAULT);

	type = XLALH5DatasetQueryType(dset);
	if (type != TCODE)
		XLAL_ERROR_NULL(XLAL_ETYPE);

	dimLength = XLALH5DatasetQueryDims(dset);
	if (!dimLength)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	array = CREATEFUNC(dimLength);
	XLALDestroyUINT4Vector(dimLength);
	if (!array)
		XLAL_ERROR_NULL(XLAL_ENOMEM);

	if (XLALH5DatasetQueryData(array->data, dset) == -1) {
		DESTROYFUNC(array);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	return array;
}

#undef CONCAT2x
#undef CONCAT2
#undef CONCAT3x
#undef CONCAT3

#undef ATYPE
#undef TCODE

#undef ALLOCFUNC
#undef READFUNC

#undef CREATEFUNC
#undef DESTROYFUNC
