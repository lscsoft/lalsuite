#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)

#define TCODE CONCAT3(LAL_,TYPECODE,_TYPE_CODE)
#define ADDFUNC CONCAT3(XLALH5DatasetAdd,TYPE,Attribute)
#define QUERYFUNC CONCAT3(XLALH5DatasetQuery,TYPE,AttributeValue)


int ADDFUNC(LALH5Dataset *dset, const char *key, TYPE value)
{
	if (!dset || !key)
		XLAL_ERROR(XLAL_EFAULT);
	return XLALH5AttributeAddScalar((LALH5Generic)dset, key, &value, TCODE);
}

TYPE QUERYFUNC(LALH5Dataset *dset, const char *key)
{
	LALTYPECODE type;
	TYPE value;
	if (!dset || !key)
		XLAL_ERROR_VAL(FAILVAL, XLAL_EFAULT);
	type = XLALH5AttributeQueryScalarType((LALH5Generic)dset, key);
	if (type != TCODE)
		XLAL_ERROR_VAL(FAILVAL, XLAL_ETYPE);
	if (XLALH5AttributeQueryScalarValue(&value, (LALH5Generic)dset, key) < 0)
		XLAL_ERROR_VAL(FAILVAL, XLAL_EFUNC);
	return value;
}

#undef CONCAT2x
#undef CONCAT2
#undef CONCAT3x
#undef CONCAT3

#undef TCODE

#undef ADDFUNC
#undef QUERYFUNC
