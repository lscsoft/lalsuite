#include <math.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/Date.h>
#include <lal/Units.h>

#define LAL_CHECK_VALID_SERIES(s,val) \
	do { \
		if ( !(s) ) XLAL_ERROR_VAL( val, XLAL_EFAULT ); \
		if ( !(s)->data || !(s)->data->data || !(s)->data->length ) XLAL_ERROR_VAL( val, XLAL_EINVAL ); \
	} while (0)

#define LAL_CHECK_CONSISTENT_TIME_SERIES(s1,s2,val) \
	do { \
		if ( XLALGPSCmp( &(s1)->epoch, &(s2)->epoch ) != 0 ) XLAL_ERROR_VAL( val, XLAL_ETIME ); \
		if ( fabs( (s1)->deltaT - (s2)->deltaT ) > LAL_REAL8_EPS ) XLAL_ERROR_VAL( val, XLAL_ETIME ); \
		if ( fabs( (s1)->f0 - (s2)->f0 ) > LAL_REAL8_EPS ) XLAL_ERROR_VAL( val, XLAL_EFREQ ); \
		if ( XLALUnitCompare( &(s1)->sampleUnits, &(s2)->sampleUnits ) ) XLAL_ERROR_VAL( val, XLAL_EUNIT ); \
		if ( (s1)->data->length != (s1)->data->length ) XLAL_ERROR_VAL(val, XLAL_EBADLEN ); \
	} while (0)

#define LAL_CHECK_COMPATIBLE_TIME_SERIES(s1,s2,val) \
	do { \
		if ( fabs( (s1)->deltaT - (s2)->deltaT ) > LAL_REAL8_EPS ) XLAL_ERROR_VAL( val, XLAL_ETIME ); \
		if ( fabs( (s1)->f0 - (s2)->f0 ) > LAL_REAL8_EPS ) XLAL_ERROR_VAL( val, XLAL_EFREQ ); \
		if ( XLALUnitCompare( &(s1)->sampleUnits, &(s2)->sampleUnits ) ) XLAL_ERROR_VAL( val, XLAL_EUNIT ); \
	} while (0)
