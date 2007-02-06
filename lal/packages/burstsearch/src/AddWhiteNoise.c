/******** <lalVerbatim file="AddWhiteNoiseCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID (ADDWHITENOISEC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/Random.h>
#include <lal/Sequence.h>
#include <lal/XLALError.h>


/*
 *  Add Gaussian white noise to various vectors
 */

/******** <lalVerbatim file="AddWhiteNoiseCP"> ********/
REAL4Sequence *XLALREAL4AddWhiteNoise(
	REAL4Sequence *sequence,
	REAL4 rms,
	RandomParams *params
)
/******** </lalVerbatim> ********/
{
	const char func[] = "XLALREAL4AddWhiteNoise";
	REAL4Sequence *noise;
	unsigned i;

	/* create temporary storage */
	noise = XLALCreateREAL4Sequence(sequence->length);
	if(!noise)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/* fill temporary sequence with Gaussian deviates and add to input */
	XLALNormalDeviates(noise, params);
	for(i = 0; i < sequence->length; i++)
		sequence->data[i] += rms * noise->data[i];

	/* clean up */
	XLALDestroyREAL4Sequence(noise);

	/* success */
	return sequence;
}


/******** <lalVerbatim file="AddWhiteNoiseCP"> ********/
COMPLEX8Sequence *XLALCOMPLEX8AddWhiteNoise(
	COMPLEX8Sequence *sequence,
	REAL8 rms,
	RandomParams *params
)
/******** </lalVerbatim> ********/
{
	const char func[] = "XLALCOMPLEX8AddWhiteNoise";
	REAL4Sequence *noise;
	unsigned i;

	/* create temporary storage */
	noise = XLALCreateREAL4Sequence(sequence->length);
	if(!noise)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/* fill temporary sequence with Gaussian deviates and add to input */
	XLALNormalDeviates(noise, params);
	for(i = 0; i < sequence->length; i++)
		sequence->data[i].re += rms * noise->data[i];

	/* repeat for imaginary component */
	XLALNormalDeviates(noise, params);
	for(i = 0; i < sequence->length; i++)
		sequence->data[i].im += rms * noise->data[i];

	/* clean up */
	XLALDestroyREAL4Sequence(noise);

	/* success */
	return sequence;
}
