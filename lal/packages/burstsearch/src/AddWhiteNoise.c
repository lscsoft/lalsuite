/******** <lalVerbatim file="AddWhiteNoiseCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID (ADDWHITENOISEC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/LALErrno.h>
#include <lal/Random.h>
#include <lal/Sequence.h>
#include <lal/XLALError.h>


/*
 *  Add white noise to complex vector
 */

/******** <lalVerbatim file="AddWhiteNoiseCP"> ********/
int XLALAddWhiteNoise(
	COMPLEX8Vector *v,
	REAL8 amplitude
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALAddWhiteNoise";
	RandomParams *params;
	REAL4Vector *noise_r, *noise_i;
	INT4 i;

	/* no-op on NULL input vector */
	if(!v)
		return(0);

	/* Seed random number generator */
	params = XLALCreateRandomParams(0);

	/* Create temporary sequences */
	noise_r = XLALCreateREAL4Sequence(v->length);
	noise_i = XLALCreateREAL4Sequence(v->length);

	/* Check for malloc failures */
	if(!params || !noise_r || !noise_i) {
		XLALDestroyRandomParams(params);
		XLALDestroyREAL4Sequence(noise_r);
		XLALDestroyREAL4Sequence(noise_i);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* Fill temporary sequences with Gaussian deviates */
	XLALNormalDeviates(noise_r, params);
	XLALNormalDeviates(noise_i, params);

	/* Add noise to input */
	for(i = 0; i < (INT4)v->length; i++) {
		v->data[i].re += amplitude * noise_r->data[i];
		v->data[i].im += amplitude * noise_i->data[i];
	}

	/* Clean up */
	XLALDestroyRandomParams(params);
	XLALDestroyREAL4Sequence(noise_r);
	XLALDestroyREAL4Sequence(noise_i);

	/* success */
	return(0);
}


/******** <lalVerbatim file="AddWhiteNoiseCP"> ********/
void
LALAddWhiteNoise(
	LALStatus *status,
	COMPLEX8Vector *v,
	REAL8 noiseLevel
)
/******** </lalVerbatim> ********/
{
	INITSTATUS (status, "LALAddWhiteNoise", ADDWHITENOISEC);
	ATTATCHSTATUSPTR (status);

	/* make sure that arguments are not NULL */
	ASSERT(v, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(v->data, status, LAL_NULL_ERR, LAL_NULL_MSG);

	/* make sure length of series is nonzero */
	ASSERT(v->length > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

	/* Wrap XLAL call in an ASSERT() */
	if(XLALAddWhiteNoise(v, noiseLevel)) {
		XLALClearErrno();
		ABORT(status, LAL_FAIL_ERR, LAL_FAIL_MSG);
	}

	/* normal exit */
	DETATCHSTATUSPTR (status);
	RETURN (status);
}
