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
	COMPLEX8Sequence *v,
	REAL8 amplitude
)
/******** </lalVerbatim> ********/
{
	const char func[] = "XLALAddWhiteNoise";
	RandomParams *params;
	REAL4Sequence *noise_r, *noise_i;
	size_t i;

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
	for(i = 0; i < v->length; i++) {
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

