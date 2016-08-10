#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)

#define VTYPE CONCAT2(TYPE,Vector)
#define STYPE CONCAT2(TYPE,FrequencySeries)

#define FILEWRITEFUNC CONCAT2(XLALH5FileWrite,STYPE)
#define FILEREADFUNC CONCAT2(XLALH5FileRead,STYPE)

#define DSETALLOCFUNC CONCAT2(XLALH5DatasetAlloc,VTYPE)
#define DSETREADFUNC CONCAT2(XLALH5DatasetRead,VTYPE)

int FILEWRITEFUNC(LALH5File *file, const char *name, STYPE *series)
{
	char sampleUnits[LALUnitTextSize];
	LALH5Dataset *dset;
	if (!file || !name || !series)
		XLAL_ERROR(XLAL_EFAULT);
	if (!series->data)
		XLAL_ERROR(XLAL_EINVAL);
	dset = DSETALLOCFUNC(file, name, series->data);
	if (!dset)
		XLAL_ERROR(XLAL_EFUNC);
	if (XLALH5AttributeAddString((LALH5Generic)dset, "name", series->name) < 0) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC, "Could not set name attribute");
	}
	if (XLALH5AttributeAddLIGOTimeGPS((LALH5Generic)dset, "epoch", &series->epoch) < 0) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC, "Could not set epoch attribute");
	}
	if (XLALH5DatasetAddREAL8Attribute(dset, "deltaF", series->deltaF) < 0) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC, "Could not set deltaF attribute");
	}
	if (XLALH5DatasetAddREAL8Attribute(dset, "f0", series->f0) < 0) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC, "Could not set f0 attribute");
	}
	if (XLALUnitAsString(sampleUnits, sizeof(sampleUnits), &series->sampleUnits) == NULL)
	{
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC);
	}
	if (XLALH5AttributeAddString((LALH5Generic)dset, "sampleUnits", sampleUnits) < 0) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC, "Could not set sampleUnits attribute");
	}
	XLALH5DatasetFree(dset);
	return 0;
}

STYPE *FILEREADFUNC(LALH5File *file, const char *name)
{
	char sampleUnits[LALUnitTextSize];
	STYPE *series;
	LALH5Dataset *dset;
	int n;

	if (!file || !name)
		XLAL_ERROR_NULL(XLAL_EFAULT);

	dset = XLALH5DatasetRead(file, name);
	if (!dset)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	series = XLALMalloc(sizeof(*series));
	if (!series) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR_NULL(XLAL_ENOMEM);
	}

	/* read metadata */

	n = XLALH5AttributeQueryStringValue(series->name, sizeof(series->name), (LALH5Generic)dset, "name");
	if (n < 0) {
		LALFree(series);
		XLALH5DatasetFree(dset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	if ((size_t)n >= sizeof(series->name))
		XLAL_PRINT_WARNING("Name of frequency series was truncated");

	n = XLALH5AttributeQueryStringValue(sampleUnits, sizeof(sampleUnits), (LALH5Generic)dset, "sampleUnits");
	if (n < 0) {
		LALFree(series);
		XLALH5DatasetFree(dset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	/* note: treat failure to parse sample unit string as a warning */
	if ((size_t)n >= sizeof(sampleUnits) || XLALParseUnitString(&series->sampleUnits, sampleUnits) == NULL) {
		XLAL_PRINT_WARNING("Could not parse unit string `%s'", sampleUnits);
		series->sampleUnits = lalDimensionlessUnit;
	}

	if (XLALH5AttributeQueryLIGOTimeGPSValue(&series->epoch, (LALH5Generic)dset, "epoch") == NULL) {
		LALFree(series);
		XLALH5DatasetFree(dset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	series->deltaF = XLALH5DatasetQueryREAL8AttributeValue(dset, "deltaF");
	if (XLAL_IS_REAL8_FAIL_NAN(series->deltaF)) {
		LALFree(series);
		XLALH5DatasetFree(dset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	series->f0 = XLALH5DatasetQueryREAL8AttributeValue(dset, "f0");
	if (XLAL_IS_REAL8_FAIL_NAN(series->f0)) {
		LALFree(series);
		XLALH5DatasetFree(dset);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	series->data = DSETREADFUNC(dset);
	XLALH5DatasetFree(dset);
	if (!series->data) {
		LALFree(series);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	return series;
}

#undef CONCAT2x
#undef CONCAT2

#undef VTYPE
#undef STYPE

#undef FILEWRITEFUNC
#undef FILEREADFUNC

#undef DSETALLOCFUNC
#undef DSETREADFUNC
