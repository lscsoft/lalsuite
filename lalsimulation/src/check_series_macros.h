#include <math.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/Date.h>
#include <lal/Units.h>

#define LAL_CHECK_VALID_SERIES(s,val) \
	do { \
		if ( !(s) ) XLAL_ERROR_VAL( __func__, XLAL_EFAULT, val ); \
		if ( !(s)->data || !(s)->data->data || !(s)->data->length ) XLAL_ERROR_VAL( __func__, XLAL_EINVAL, val ); \
	} while (0)

#define LAL_CHECK_CONSISTENT_TIME_SERIES(s1,s2,val) \
	do { \
		if ( XLALGPSCmp( &(s1)->epoch, &(s2)->epoch ) != 0 ) XLAL_ERROR_VAL( __func__, XLAL_ETIME, val ); \
		if ( fabs( (s1)->deltaT - (s2)->deltaT ) > LAL_REAL8_EPS ) XLAL_ERROR_VAL( __func__, XLAL_ETIME, val ); \
		if ( fabs( (s1)->f0 - (s2)->f0 ) > LAL_REAL8_EPS ) XLAL_ERROR_VAL( __func__, XLAL_EFREQ, val ); \
		if ( XLALUnitCompare( &(s1)->sampleUnits, &(s2)->sampleUnits ) ) XLAL_ERROR_VAL( __func__, XLAL_EUNIT, val ); \
		if ( (s1)->data->length != (s1)->data->length ) XLAL_ERROR_VAL( __func__, XLAL_EBADLEN, val ); \
	} while (0)
