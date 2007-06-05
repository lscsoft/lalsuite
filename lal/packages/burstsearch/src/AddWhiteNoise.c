/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#include <lal/LALRCSID.h>


NRCSID (ADDWHITENOISEC, "$Id$");


#include <lal/EPSearch.h>	/* for our own prototypes */
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
