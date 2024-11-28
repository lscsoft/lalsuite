#include <lal/LALConfig.h>

#ifndef LAL_HDF5_ENABLED
int main(void) { return 77; /* don't do any testing */ }
#else

#if defined(GENERATE_HDF5_TEST_FILE)
#define H5_USE_110_API
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/StringVector.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/H5FileIO.h>

#define FNAME "test.h5"
#define GROUP "path/to"
#define DSET  "testdset"

#define NDIM 3
#define DIM0 2
#define DIM1 4
#define DIM2 3
#define NPTS (DIM0 * DIM1 * DIM2)

#define EPOCH { 123456789, 987654321 }
static LIGOTimeGPS epoch = EPOCH;
static LIGOTimeGPS now;

static UINT2 count;

static int generate_int_data(void)
{
	return rand() % SHRT_MAX;
}

static float generate_float_data(void)
{
	return rand() / (RAND_MAX + 1.0);
}

static float complex generate_complex_data(void)
{
	return rand() / (RAND_MAX + 1.0) + I * rand() / (RAND_MAX + 1.0);
}

static char * generate_string_data(void)
{
	const char *s;
	switch (rand() % 18) {
	case 0: s = "lorem"; break;
	case 1: s = "ipsum"; break;
	case 2: s = "dolor"; break;
	case 3: s = "sit"; break;
	case 4: s = "amet"; break;
	case 5: s = "consectetur"; break;
	case 6: s = "adipiscing"; break;
	case 7: s = "elit"; break;
	case 8: s = "sed"; break;
	case 9: s = "eiusmod"; break;
	case 10: s = "tempor"; break;
	case 11: s = "incididunt"; break;
	case 12: s = "ut"; break;
	case 13: s = "labore"; break;
	case 14: s = "et"; break;
	case 15: s = "dolore"; break;
	case 16: s = "magna"; break;
	case 17: s = "aliqua"; break;
	default:
		fprintf(stderr, "Cannot get here!");
		exit(1);
	}
	return XLALStringDuplicate(s);
}

#define DEFINE_WRITE_FUNCTION(type) \
	static int write_ ## type (type *x) { \
		LALH5File *file; \
		LALH5File *group; \
		XLALGPSTimeNow(&now); \
		++count; \
		file = XLALH5FileOpen(FNAME, "w"); \
		XLALH5FileAddLIGOTimeGPSAttribute(file, "creation_time_gps", &now); \
		XLALH5FileAddScalarAttribute(file, "test_count", &count, LAL_U2_TYPE_CODE); \
		group = XLALH5GroupOpen(file, GROUP); \
		XLALH5FileAddStringAttribute(group, "test_data_type", #type); \
		XLALH5FileWrite ## type (group, DSET, x); \
		XLALH5FileClose(group); \
		XLALH5FileClose(file); \
		return 0; \
	}

#define DEFINE_READ_FUNCTION(type) \
	static type * read_ ## type (void) { \
		type *x; \
		LALH5File *file; \
		LALH5File *group; \
		LIGOTimeGPS creation_time; \
		UINT2 cnt; \
		char t[sizeof(#type)]; \
		file = XLALH5FileOpen(FNAME, "r"); \
		XLALH5FileQueryLIGOTimeGPSAttributeValue(&creation_time, file, "creation_time_gps"); \
		if (XLALGPSCmp(&creation_time, &now)) { \
			fprintf(stderr, " FAIL\n"); \
			exit(1); /* fail */ \
		} \
		XLALH5FileQueryScalarAttributeValue(&cnt, file, "test_count"); \
		if (cnt != count) { \
			fprintf(stderr, " FAIL\n"); \
			exit(1); /* fail */ \
		} \
		group = XLALH5GroupOpen(file, GROUP); \
		XLALH5FileQueryStringAttributeValue(t, sizeof(t), group, "test_data_type"); \
		XLALH5FileClose(group); \
		if (strcmp(t, #type)) { \
			fprintf(stderr, " FAIL\n"); \
			exit(1); /* fail */ \
		} \
		x = XLALH5FileRead ## type (file, GROUP "/" DSET); \
		XLALH5FileClose(file); \
		return x; \
	}

#define DEFINE_TEST_FUNCTION(type) \
	static void test_ ## type (void) { \
		type *orig; \
		type *copy; \
		fprintf(stderr, "Testing Read/Write of %s...", #type); \
		orig = create_ ## type (); \
		write_ ## type (orig); \
		copy = read_ ## type (); \
		if (compare_ ## type (orig, copy)) { \
			fprintf(stderr, " FAIL\n"); \
			exit(1); /* fail */ \
		} \
		XLALDestroy ## type (copy); \
		XLALDestroy ## type (orig); \
		fprintf(stderr, " PASS\n"); \
	}

/* VECTOR ROUTINES */

#define DEFINE_CREATE_VECTOR_FUNCTION(type) \
	static type * create_ ## type (void) { \
		type *v; \
		size_t i; \
		v = XLALCreate ## type (NPTS); \
		for (i = 0; i < NPTS; ++i) \
			v->data[i] = GENERATE_DATA(); \
		return v; \
	}

#define DEFINE_COMPARE_VECTOR_FUNCTION(type) \
	static int compare_ ## type (type *v1, type *v2) { \
		if (v1->length != v2->length) \
			return 1; \
		return memcmp(v1->data, v2->data, v1->length * sizeof(*v1->data)); \
	}

static LALStringVector * create_StringVector(void)
{
	LALStringVector *v;
	size_t i;
	v = XLALCreateEmptyStringVector(NPTS);
	for (i = 0; i < NPTS; ++i)
		v->data[i] = generate_string_data();
	return v;
}

static int compare_StringVector(LALStringVector *v1, LALStringVector *v2)
{
	if (v1->length != v2->length)
		return 1;
	for (size_t i = 0; i < v1->length; ++i)
		if (strcmp(v1->data[i], v2->data[i]) != 0)
			return 1;
	return 0;
}

/* ARRAY ROUTINES */

#define DEFINE_CREATE_ARRAY_FUNCTION(type) \
	static type * create_ ## type (void) { \
		type *a; \
		size_t i; \
		a = XLALCreate ## type ## L(NDIM, DIM0, DIM1, DIM2); \
		for (i = 0; i < NPTS; ++i) \
			a->data[i] = GENERATE_DATA(); \
		return a; \
	}

#define DEFINE_COMPARE_ARRAY_FUNCTION(type) \
	static int compare_ ## type (type *a1, type *a2) { \
		size_t sz1 = 1, sz2 = 1, i; \
		if (a1->dimLength->length != a2->dimLength->length) \
			return 1; \
		for (i = 0; i < a1->dimLength->length; ++i) { \
			if (a1->dimLength->data[i] != a2->dimLength->data[i]) \
				return 1; \
			sz1 *= a1->dimLength->data[i]; \
			sz2 *= a2->dimLength->data[i]; \
		} \
		if (sz1 != sz2) \
			return 1; \
		return memcmp(a1->data, a2->data, sz1 * sizeof(*a1->data)); \
	}

/* TIME/FREQUENCY SERIES ROUTINES */

#define DEFINE_CREATE_TIME_FREQUENCY_SERIES_FUNCTION(type) \
	static type * create_ ## type (void) { \
		type *s; \
		size_t i; \
		s = XLALCreate ## type ("test_" #type, &epoch, 0.0, 0.1, &lalStrainUnit, NPTS); \
		for (i = 0; i < NPTS; ++i) \
			s->data->data[i] = GENERATE_DATA(); \
		return s; \
	}

#define DEFINE_COMPARE_TIME_SERIES_FUNCTION(type) \
	static int compare_ ## type (type *s1, type *s2) { \
		int c; \
		if (s1->data->length != s2->data->length) \
			return 1; \
		if ((c = strcmp(s1->name, s2->name))) \
			return c; \
		if ((c = (s1->deltaT > s2->deltaT) - (s1->deltaT < s2->deltaT))) \
			return c; \
		if ((c = (s1->f0 > s2->f0) - (s1->f0 < s2->f0))) \
			return c; \
		if (XLALUnitCompare(&s1->sampleUnits, &s2->sampleUnits)) \
			return 1; \
		return memcmp(s1->data->data, s2->data->data, s1->data->length * sizeof(*s1->data->data)); \
	}

#define DEFINE_COMPARE_FREQUENCY_SERIES_FUNCTION(type) \
	static int compare_ ## type (type *s1, type *s2) { \
		int c; \
		if (s1->data->length != s2->data->length) \
			return 1; \
		if ((c = strcmp(s1->name, s2->name))) \
			return c; \
		if ((c = (s1->deltaF > s2->deltaF) - (s1->deltaF < s2->deltaF))) \
			return c; \
		if ((c = (s1->f0 > s2->f0) - (s1->f0 < s2->f0))) \
			return c; \
		if (XLALUnitCompare(&s1->sampleUnits, &s2->sampleUnits)) \
			return 1; \
		return memcmp(s1->data->data, s2->data->data, s1->data->length * sizeof(*s1->data->data)); \
	}

#define DEFINE_VECTOR_FUNCTIONS(type) \
	DEFINE_WRITE_FUNCTION(type) \
	DEFINE_READ_FUNCTION(type) \
	DEFINE_COMPARE_VECTOR_FUNCTION(type) \
	DEFINE_CREATE_VECTOR_FUNCTION(type) \
	DEFINE_TEST_FUNCTION(type)

#define DEFINE_ARRAY_FUNCTIONS(type) \
	DEFINE_WRITE_FUNCTION(type) \
	DEFINE_READ_FUNCTION(type) \
	DEFINE_COMPARE_ARRAY_FUNCTION(type) \
	DEFINE_CREATE_ARRAY_FUNCTION(type) \
	DEFINE_TEST_FUNCTION(type)

#define DEFINE_TIME_SERIES_FUNCTIONS(type) \
	DEFINE_WRITE_FUNCTION(type) \
	DEFINE_READ_FUNCTION(type) \
	DEFINE_COMPARE_TIME_SERIES_FUNCTION(type) \
	DEFINE_CREATE_TIME_FREQUENCY_SERIES_FUNCTION(type) \
	DEFINE_TEST_FUNCTION(type)

#define DEFINE_FREQUENCY_SERIES_FUNCTIONS(type) \
	DEFINE_WRITE_FUNCTION(type) \
	DEFINE_READ_FUNCTION(type) \
	DEFINE_COMPARE_FREQUENCY_SERIES_FUNCTION(type) \
	DEFINE_CREATE_TIME_FREQUENCY_SERIES_FUNCTION(type) \
	DEFINE_TEST_FUNCTION(type)

#define GENERATE_DATA generate_int_data
DEFINE_VECTOR_FUNCTIONS(CHARVector)
DEFINE_VECTOR_FUNCTIONS(INT2Vector)
DEFINE_VECTOR_FUNCTIONS(INT4Vector)
DEFINE_VECTOR_FUNCTIONS(INT8Vector)
DEFINE_VECTOR_FUNCTIONS(UINT2Vector)
DEFINE_VECTOR_FUNCTIONS(UINT4Vector)
DEFINE_VECTOR_FUNCTIONS(UINT8Vector)
#undef GENERATE_DATA
#define GENERATE_DATA generate_float_data
DEFINE_VECTOR_FUNCTIONS(REAL4Vector)
DEFINE_VECTOR_FUNCTIONS(REAL8Vector)
#undef GENERATE_DATA
#define GENERATE_DATA generate_complex_data
DEFINE_VECTOR_FUNCTIONS(COMPLEX8Vector)
DEFINE_VECTOR_FUNCTIONS(COMPLEX16Vector)
#undef GENERATE_DATA

#define StringVector LALStringVector
DEFINE_WRITE_FUNCTION(StringVector)
DEFINE_READ_FUNCTION(StringVector)
DEFINE_TEST_FUNCTION(StringVector)

#define GENERATE_DATA generate_int_data
DEFINE_ARRAY_FUNCTIONS(INT2Array)
DEFINE_ARRAY_FUNCTIONS(INT4Array)
DEFINE_ARRAY_FUNCTIONS(INT8Array)
DEFINE_ARRAY_FUNCTIONS(UINT2Array)
DEFINE_ARRAY_FUNCTIONS(UINT4Array)
DEFINE_ARRAY_FUNCTIONS(UINT8Array)
#undef GENERATE_DATA
#define GENERATE_DATA generate_float_data
DEFINE_ARRAY_FUNCTIONS(REAL4Array)
DEFINE_ARRAY_FUNCTIONS(REAL8Array)
#undef GENERATE_DATA
#define GENERATE_DATA generate_complex_data
DEFINE_ARRAY_FUNCTIONS(COMPLEX8Array)
DEFINE_ARRAY_FUNCTIONS(COMPLEX16Array)
#undef GENERATE_DATA

#define GENERATE_DATA generate_int_data
DEFINE_TIME_SERIES_FUNCTIONS(INT2TimeSeries)
DEFINE_TIME_SERIES_FUNCTIONS(INT4TimeSeries)
DEFINE_TIME_SERIES_FUNCTIONS(INT8TimeSeries)
DEFINE_TIME_SERIES_FUNCTIONS(UINT2TimeSeries)
DEFINE_TIME_SERIES_FUNCTIONS(UINT4TimeSeries)
DEFINE_TIME_SERIES_FUNCTIONS(UINT8TimeSeries)
#undef GENERATE_DATA
#define GENERATE_DATA generate_float_data
DEFINE_TIME_SERIES_FUNCTIONS(REAL4TimeSeries)
DEFINE_TIME_SERIES_FUNCTIONS(REAL8TimeSeries)
#undef GENERATE_DATA
#define GENERATE_DATA generate_complex_data
DEFINE_TIME_SERIES_FUNCTIONS(COMPLEX8TimeSeries)
DEFINE_TIME_SERIES_FUNCTIONS(COMPLEX16TimeSeries)
#undef GENERATE_DATA

#define GENERATE_DATA generate_float_data
DEFINE_FREQUENCY_SERIES_FUNCTIONS(REAL4FrequencySeries)
DEFINE_FREQUENCY_SERIES_FUNCTIONS(REAL8FrequencySeries)
#undef GENERATE_DATA
#define GENERATE_DATA generate_complex_data
DEFINE_FREQUENCY_SERIES_FUNCTIONS(COMPLEX8FrequencySeries)
DEFINE_FREQUENCY_SERIES_FUNCTIONS(COMPLEX16FrequencySeries)
#undef GENERATE_DATA

#if defined(GENERATE_HDF5_TEST_FILE)
/* this small snippet is left as an example on how to regenerate the files
   used for testing HDF5 attributes reading/writing
*/
static void create_hdf5(void) {
	char const *buf = "La soupe est pr\xc3\xaate.";

	hid_t file_create = H5Fcreate(TEST_DATA_DIR "hdf5_attr_utf8.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t fspace = H5Screate(H5S_SCALAR);
	hid_t memtype_id = H5Tcopy(H5T_C_S1);
	H5Tset_cset(memtype_id, H5T_CSET_UTF8);
	H5Tset_size(memtype_id, H5T_VARIABLE);

	hid_t attr = H5Acreate2(file_create, "CANONICAL_FILE_BASENAME", memtype_id, fspace, H5P_DEFAULT, H5P_DEFAULT);

	if(!(H5Awrite(attr, memtype_id, &buf) >= 0)){
		fprintf(stderr, "Attribute write fail\n");
		exit(1);
	}
	if(!(H5Sclose(fspace) >= 0)){
		fprintf(stderr, "Space close fail\n");
		exit(1);
	}
	if(!(H5Tclose(memtype_id) >= 0)){
		fprintf(stderr, "Memtype close fail\n");
		exit(1);

	}
	if(!(H5Aclose(attr) >= 0)){
		fprintf(stderr, "Attribute close fail\n");
		exit(1);
	}

	H5Fflush(file_create, H5F_SCOPE_GLOBAL);
	H5Fclose(file_create);
}
#endif

static void check_string_reading_from_hdf5_attr(
	char const* filename,
	char const* attribute_name,
	char const* buf,
	char const* attr_type) {

	fprintf(stderr, "Testing HDF5 attribute reading %s...\n", attr_type);
	LALH5File *file_read = XLALH5FileOpen(filename, "r");
	LALH5Generic gfile = {.file = file_read};

	int len = XLALH5AttributeQueryStringValue(NULL, 0, gfile, attribute_name);
	if(!(len > 0)) {
		fprintf(stderr, "Attribute length read incorrect\n");
		exit(1);
	}
	char *read_string_back = read_string_back = XLALCalloc(1, len + 1 );
	XLALH5AttributeQueryStringValue(read_string_back, len+1, gfile, attribute_name);

	if(strcmp(read_string_back, buf) != 0){
		fprintf(stderr, "Read string incorrect '%s'\n", read_string_back);
		exit(1);
	}

	XLALFree(read_string_back);
	XLALH5FileClose(file_read);
	fprintf(stderr, "Testing HDF5 attribute reading %s... ok!\n", attr_type);
}

int main(void)
{
	XLALSetErrorHandler(XLALAbortErrorHandler);

	test_CHARVector();
	test_INT2Vector();
	test_INT4Vector();
 	test_INT8Vector();
	test_UINT2Vector();
	test_UINT4Vector();
	test_UINT8Vector();
	test_REAL4Vector();
	test_REAL8Vector();
	test_COMPLEX8Vector();
	test_COMPLEX16Vector();
	test_StringVector();

	test_INT2Array();
	test_INT4Array();
	test_INT8Array();
	test_UINT2Array();
	test_UINT4Array();
	test_UINT8Array();
	test_REAL4Array();
	test_REAL8Array();
	test_COMPLEX8Array();
	test_COMPLEX16Array();

	test_INT2TimeSeries();
	test_INT4TimeSeries();
	test_INT8TimeSeries();
	test_UINT2TimeSeries();
	test_UINT4TimeSeries();
	test_UINT8TimeSeries();
	test_REAL4TimeSeries();
	test_REAL8TimeSeries();
	test_COMPLEX8TimeSeries();
	test_COMPLEX16TimeSeries();

	test_REAL4FrequencySeries();
	test_REAL8FrequencySeries();
	test_COMPLEX8FrequencySeries();
	test_COMPLEX16FrequencySeries();

#if defined(GENERATE_HDF5_TEST_FILE)
	/* define the macro GENERATE_HDF5_TEST_FILE to generate each file */
	create_hdf5();
#endif

	check_string_reading_from_hdf5_attr(
		TEST_DATA_DIR "hdf5_attr_ascii.h5",
		"CANONICAL_FILE_BASENAME",
		"test long string ascii",
		"ascii"
	);

	check_string_reading_from_hdf5_attr(
		TEST_DATA_DIR "hdf5_attr_utf8.h5",
		"CANONICAL_FILE_BASENAME",
		"La soupe est pr\xc3\xaate.",
		"utf8"
	);

	LALCheckMemoryLeaks();
	return 0;
}

#endif /* ! HAVE_HDF5 */
