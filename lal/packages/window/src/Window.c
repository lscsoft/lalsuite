/*
 * Copyright (C) 2007 Bruce Allen, Duncan Brown, Jolien Creighton, Kipp
 * Cannon, Patrick Brady, Teviet Creighton
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_bessel.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/Sequence.h>
#include <lal/Window.h>
#include <lal/XLALError.h>

/*
 * ============================================================================
 *
 *                             Private utilities
 *
 * ============================================================================
 */


/**
 * Constructs a REAL4Window from a REAL8Window by quantizing the
 * double-precision data to single-precision.  The REAL8Window is freed
 * unconditionally.  Intended to be used as a wrapper, to convert any
 * function that constructs a REAL8Window into a function to construct a
 * REAL4Window.
 */
static REAL4Window *XLALREAL4Window_from_REAL8Window(REAL8Window *orig)
{
	REAL4Window *new;
	REAL4Sequence *data;
	UINT4 i;

	if(!orig)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new = XLALMalloc(sizeof(*new));
	data = XLALCreateREAL4Sequence(orig->data->length);

	if(!new || !data) {
		XLALDestroyREAL8Window(orig);
		XLALFree(new);
		XLALDestroyREAL4Sequence(data);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	for(i = 0; i < data->length; i++)
		data->data[i] = orig->data->data[i];
	new->data = data;
	new->sumofsquares = orig->sumofsquares;
	new->sum = orig->sum;

	XLALDestroyREAL8Window(orig);

	return new;
}


/**
 * Maps the length of a window and the offset within the window to the "y"
 * co-ordinate of the LAL documentation.
 *
 * Input:
 * length > 0,
 * 0 <= i < length
 *
 * Output:
 * length < 2 --> return 0.0
 * i == 0 --> return -1.0
 * i == (length - 1) / 2 --> return 0.0
 * i == length - 1 --> return +1.0
 *
 * e.g., length = 5 (odd), then i == 2 --> return 0.0
 * if length = 6 (even), then i == 2.5 --> return 0.0
 *
 * (in the latter case, obviously i can't be a non-integer, but that's the
 * value it would have to be for this function to return 0.0)
 */
static double Y(int length, int i)
{
	length -= 1;
	return length > 0 ? (2 * i - length) / (double) length : 0;
}


/**
 * Computes the sum of squares, and sum, of the samples in a window
 * function.
 *
 * Two techniques are employed to achieve accurate results.  Firstly, the
 * loop iterates from the edges to the centre.  Generally, window functions
 * have smaller values at the edges and larger values in the middle, so
 * adding them in this order avoids adding small numbers to big numbers.
 * The loops also implement Kahan's compensated summation algorithm in
 * which a second variable is used to accumulate round-off errors and fold
 * them into later iterations.
 */
static REAL8 sum_squares(REAL8 *start, int length)
{

	REAL8 sum = 0.0;
	REAL8 e = 0.0;
	REAL8 *end = start + length - 1;

	/* Note:  don't try to simplify this loop, the statements must be
	 * done like this to induce the C compiler to perform the
	 * arithmetic in the correct order. */
	for(; start < end; start++, end--) {
		REAL8 temp = sum;
		/* what we want to add */
		REAL8 x = *start * *start + *end * *end + e;
		/* add */
		sum += x;
		/* negative of what was actually added */
		e = temp - sum;
		/* what didn't get added, add next time */
		e += x;
	}

	if(start == end)
		/* length is odd, get middle sample */
		sum += *start * *start + e;

	return sum;
}


static REAL8 sum_samples(REAL8 *start, int length)
{

	REAL8 sum = 0.0;
	REAL8 e = 0.0;
	REAL8 *end = start + length - 1;

	/* Note:  don't try to simplify this loop, the statements must be
	 * done like this to induce the C compiler to perform the
	 * arithmetic in the correct order. */
	for(; start < end; start++, end--) {
		REAL8 temp = sum;
		/* what we want to add */
		REAL8 x = *start + *end + e;
		/* add */
		sum += x;
		/* negative of what was actually added */
		e = temp - sum;
		/* what didn't get added, add next time */
		e += x;
	}

	if(start == end)
		/* length is odd, get middle sample */
		sum += *start + e;

	return sum;
}


/*
 * ============================================================================
 *
 *                              Public Utilities
 *
 * ============================================================================
 */


/**
 * Constructs a new REAL8Window from a REAL8Sequence.  The window "owns"
 * the sequence, when the window is destroyed the sequence will be
 * destroyed with it.  If this function fails, the sequence is destroyed.
 * The return value is the address of the newly allocated REAL8Window or
 * NULL on failure.
 */
REAL8Window *XLALCreateREAL8WindowFromSequence(REAL8Sequence *sequence)
{
	REAL8Window *new;

	new = XLALMalloc(sizeof(*new));
	if(!new) {
		XLALDestroyREAL8Sequence(sequence);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	new->data = sequence;
	new->sumofsquares = sum_squares(new->data->data, new->data->length);
	new->sum = sum_samples(new->data->data, new->data->length);

	return new;
}


/**
 * Single-precision version of XLALCreateREAL8WindowFromSequence().
 */
REAL4Window *XLALCreateREAL4WindowFromSequence(REAL4Sequence *sequence)
{
	/* a double-precision copy of the data is used from which to
	 * compute the window's metadata.  this provides for more accurate
	 * results, but mostly means I don't have to write single-precision
	 * versions of the summing loops */
	REAL8Sequence *workspace;
	REAL4Window *new;
	UINT4 i;

	workspace = XLALCreateREAL8Sequence(sequence->length);
	new = XLALMalloc(sizeof(*new));
	if(!workspace || !new) {
		XLALDestroyREAL4Sequence(sequence);
		XLALDestroyREAL8Sequence(workspace);
		XLALFree(new);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	for(i = 0; i < workspace->length; i++)
		workspace->data[i] = sequence->data[i];

	new->data = sequence;
	new->sumofsquares = sum_squares(workspace->data, workspace->length);
	new->sum = sum_samples(workspace->data, workspace->length);

	XLALDestroyREAL8Sequence(workspace);

	return new;
}


/**
 * Multiply a REAL8Sequence in-place by a REAL8Window with a normalization
 * that preserves the variance of a zero-mean stationary Gaussian random
 * process.  If the window's length is N samples and its sum-of-squares is
 * S, then the input sequence is multiplied by the window * \f$\sqrt{N / S}\f$.
 * Returns the address of the REAL8Sequence or NULL on failure.
 */
REAL8Sequence *XLALUnitaryWindowREAL8Sequence(REAL8Sequence *sequence, const REAL8Window *window)
{
	unsigned i;
	double norm = sqrt(window->data->length / window->sumofsquares);

	if(window->sumofsquares <= 0)
		XLAL_ERROR_NULL(XLAL_EDOM);
	if(sequence->length != window->data->length)
		XLAL_ERROR_NULL(XLAL_EBADLEN);

	for(i = 0; i < window->data->length; i++)
		sequence->data[i] *= window->data->data[i] * norm;

	return sequence;
}


/**
 * Double-precision complex version of XLALUnitaryWindowREAL8Sequence().
 */
COMPLEX16Sequence *XLALUnitaryWindowCOMPLEX16Sequence(COMPLEX16Sequence *sequence, const REAL8Window *window)
{
	unsigned i;
	double norm = sqrt(window->data->length / window->sumofsquares);

	if(window->sumofsquares <= 0)
		XLAL_ERROR_NULL(XLAL_EDOM);
	if(sequence->length != window->data->length)
		XLAL_ERROR_NULL(XLAL_EBADLEN);

	for(i = 0; i < window->data->length; i++)
		sequence->data[i] *= window->data->data[i] * norm;

	return sequence;
}


/**
 * Single-precision version of XLALUnitaryWindowREAL8Sequence().
 */
REAL4Sequence *XLALUnitaryWindowREAL4Sequence(REAL4Sequence *sequence, const REAL4Window *window)
{
	unsigned i;
	float norm = sqrt(window->data->length / window->sumofsquares);

	if(window->sumofsquares <= 0)
		XLAL_ERROR_NULL(XLAL_EDOM);
	if(sequence->length != window->data->length)
		XLAL_ERROR_NULL(XLAL_EBADLEN);

	for(i = 0; i < window->data->length; i++)
		sequence->data[i] *= window->data->data[i] * norm;

	return sequence;
}


/**
 * Single-precision complex version of XLALUnitaryWindowREAL8Sequence().
 */
COMPLEX8Sequence *XLALUnitaryWindowCOMPLEX8Sequence(COMPLEX8Sequence *sequence, const REAL4Window *window)
{
	unsigned i;
	double norm = sqrt(window->data->length / window->sumofsquares);

	if(window->sumofsquares <= 0)
		XLAL_ERROR_NULL(XLAL_EDOM);
	if(sequence->length != window->data->length)
		XLAL_ERROR_NULL(XLAL_EBADLEN);

	for(i = 0; i < window->data->length; i++)
		sequence->data[i] *= window->data->data[i] * norm;

	return sequence;
}


/*
 * ============================================================================
 *
 *                                REAL8Window
 *
 * ============================================================================
 */


REAL8Window *XLALCreateRectangularREAL8Window(UINT4 length)

{
	REAL8Sequence *sequence;
	UINT4 i;

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* flat, box-car, top-hat, rectangle, whatever */
	for(i = 0; i < length; i++)
		sequence->data[i] = 1;

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreateHannREAL8Window(UINT4 length)

{
	REAL8Sequence *sequence;
	UINT4 i;

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* cos^2, zero at both end points, 1 in the middle */
	for(i = 0; i < (length + 1) / 2; i++)
		sequence->data[i] = sequence->data[length - 1 - i] = pow(cos(LAL_PI_2 * Y(length, i)), 2);

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreateWelchREAL8Window(UINT4 length)

{
	REAL8Sequence *sequence;
	UINT4 i;

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* downward-opening parabola, zero at both end points, 1 in the
	 * middle */
	for(i = 0; i < (length + 1) / 2; i++)
		sequence->data[i] = sequence->data[length - 1 - i] = 1 - pow(Y(length, i), 2.0);

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreateBartlettREAL8Window(UINT4 length)

{
	REAL8Sequence *sequence;
	UINT4 i;

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* downward-opening triangle, zero at both end points (non-zero end
	 * points is a different window called the "triangle" window), 1 in
	 * the middle */
	for(i = 0; i < (length + 1) / 2; i++)
		sequence->data[i] = sequence->data[length - 1 - i] = 1 + Y(length, i);

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreateParzenREAL8Window(UINT4 length)

{
	REAL8Sequence *sequence;
	UINT4 i;

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* ?? Copied from LAL Software Description */
	for(i = 0; i < (length + 1) / 4; i++)
		sequence->data[i] = sequence->data[length - 1 - i] = 2 * pow(1 + Y(length, i), 3);
	for(; i < (length + 1) / 2; i++) {
		double y = Y(length, i);
		sequence->data[i] = sequence->data[length - 1 - i] = 1 - 6 * y * y * (1 + y);
	}

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreatePapoulisREAL8Window(UINT4 length)

{
	REAL8Sequence *sequence;
	UINT4 i;

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* ?? Copied from LAL Software Description */
	for(i = 0; i < (length + 1) / 2; i++) {
		double y = Y(length, i);
		sequence->data[i] = sequence->data[length - 1 - i] = (1 + y) * cos(LAL_PI * y) - sin(LAL_PI * y) / LAL_PI;
	}

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreateHammingREAL8Window(UINT4 length)

{
	REAL8Sequence *sequence;
	UINT4 i;

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* cos^2, like Hann window, but with a bias of 0.08 */
	for(i = 0; i < (length + 1) / 2; i++)
		sequence->data[i] = sequence->data[length - 1 - i] = 0.08 + 0.92 * pow(cos(LAL_PI_2 * Y(length, i)), 2);

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreateKaiserREAL8Window(UINT4 length, REAL8 beta)

{
	REAL8Sequence *sequence;
	REAL8 I0beta=0;
	UINT4 i;

	if(beta < 0)
		XLAL_ERROR_NULL(XLAL_ERANGE);

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* pre-compute I0(beta) */
	if(beta < 705)
		I0beta = gsl_sf_bessel_I0(beta);

	/* I0(beta sqrt(1 - y^2)) / I0(beta)
	 *
	 * note that in many places the window is defined with pi
	 * multiplying beta in the numerator and denominator, but not here.
	 *
	 * The asymptotic forms for large beta are derived from the
	 * asymptotic form of I0(x) which is
	 *
	 * I0(x) --> exp(x) / sqrt(2 pi x)
	 *
	 * Although beta may be large, beta sqrt(1 - y^2) can be small
	 * (when y ~= +/- 1), so there is a need for two asymptotic forms:
	 * one for large beta alone, and one for large beta sqrt(1 - y^2).
	 *
	 * When beta alone is large,
	 *
	 * w(y) = I0(beta sqrt(1 - y^2)) sqrt(2 pi beta) / exp(beta)
	 *
	 * and when beta sqrt(1 - y^2) is large,
	 *
	 * w(y) = exp(-beta * (1 - sqrt(1 - y^2))) / sqrt(1 - y^2)
	 *
	 * As a function of beta, the asymptotic approximation and the
	 * "exact" form are found to disagree by about 20% in the y = +/- 1
	 * bins near the edge of the "exact" form's domain of validity.  To
	 * smooth this out, a linear transition to the asymptotic form
	 * occurs between beta = 695 and beta = 705. */

	for(i = 0; i < (length + 1) / 2; i++) {
		double y = Y(length, i);
		double x = sqrt(1 - y * y);
		double w1=0, w2=0;

		if(beta < 705)
			w1 = gsl_sf_bessel_I0(beta * x) / I0beta;
		if(beta >= 695) {
			/* FIXME:  should an interpolation be done across
			 * the transition from small beta x to large beta
			 * x? */
			/* Note:  the inf * 0 when beta = inf and y = +/- 1
			 * needs to be hard-coded */
			if(beta * x < 700)
				w2 = y == 0 ? 1 : gsl_sf_bessel_I0(beta * x) * sqrt(LAL_2_PI * beta) / exp(beta);
			else
				/* Note:  when beta = inf, the inf * 0 in
				 * the y = 0 sample must be hard-coded,
				 * which we do by simply testing for y = 0.
				 * And when beta = inf and x = 0 (y = +/-
				 * 1), the conditional goes onto this
				 * branch, and results in a 0/0 which we
				 * have to hard-code */
				w2 = y == 0 ? 1 : x == 0 ? 0 : exp(-beta * (1 - x)) / sqrt(x);
		}

		if(beta < 695)
			sequence->data[i] = sequence->data[length - 1 - i] = w1;
		else if(beta < 705) {
			double r = (beta - 695) / (705 - 695);
			sequence->data[i] = sequence->data[length - 1 - i] = (1 - r) * w1 + r * w2;
		} else
			sequence->data[i] = sequence->data[length - 1 - i] = w2;
	}

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreateCreightonREAL8Window(UINT4 length, REAL8 beta)

{
	REAL8Sequence *sequence;
	UINT4 i;

	if(beta < 0)
		XLAL_ERROR_NULL(XLAL_ERANGE);

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* ?? Copied from LAL Software Description */
	for(i = 0; i < (length + 1) / 2; i++) {
		double y = Y(length, i);
		/* NOTE:  divide-by-zero in y^2 / (1 - y^2) when i = 0
		 * seems to work out OK.  It's well-defined algebraically,
		 * but I'm surprised the FPU doesn't complain.  The 0/0
		 * when beta = i = 0 has to be hard-coded, as does the inf
		 * * 0 when beta = inf and y = 0 (which is done by just
		 * checking for y = 0). The fabs() is there because Macs,
		 * with optimizations turned on, incorrectly state that
		 * 1-y^2 = -0 when y = 1, which converts the argument of
		 * exp() from -inf to +inf, and causes the window to
		 * evaluate to +inf instead of 0 at the end points. See
		 * also the -0 at the end points of the Welch window on
		 * Macs. */
		sequence->data[i] = sequence->data[length - 1 - i] = (beta == 0 && y == -1) || y == 0 ? 1 : exp(-beta * y * y / fabs(1 - y * y));
	}

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreateTukeyREAL8Window(UINT4 length, REAL8 beta)

{
	REAL8Sequence *sequence;
	UINT4 transition_length = beta * length + 0.5;
	UINT4 i;

	if(beta < 0 || beta > 1)
		XLAL_ERROR_NULL(XLAL_ERANGE);

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* 1.0 and flat in the middle, cos^2 transition at each end, zero
	 * at end points, 0.0 <= beta <= 1.0 sets what fraction of the
	 * window is transition (0 --> rectangle window, 1 --> Hann window)
	 * */
	for(i = 0; i < (transition_length + 1) / 2; i++)
		sequence->data[i] = sequence->data[length - 1 - i] = pow(cos(LAL_PI_2 * Y(transition_length, i)), 2);
	for(; i < (length + 1) / 2; i++)
		sequence->data[i] = sequence->data[length - 1 - i] = 1;

	return XLALCreateREAL8WindowFromSequence(sequence);
}



REAL8Window *XLALCreateGaussREAL8Window(UINT4 length, REAL8 beta)

{
	REAL8Sequence *sequence;
	UINT4 i;

	if(beta < 0)
		XLAL_ERROR_NULL(XLAL_ERANGE);

	sequence = XLALCreateREAL8Sequence(length);
	if(!sequence)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/* pre-compute -1/2 beta^2 */
	beta = -0.5 * beta * beta;

	/* exp(-1/2 beta^2 y^2) */
	for(i = 0; i < (length + 1) / 2; i++) {
		double y = Y(length, i);
		/* Note:  we have to hard-code the 0 * inf when y = 0 and
		 * beta = inf, which we do by simply checking for y = 0 */
		sequence->data[i] = sequence->data[length - 1 - i] = y == 0 ? 1 : exp(y * y * beta);
	}

	return XLALCreateREAL8WindowFromSequence(sequence);
}



void XLALDestroyREAL8Window(REAL8Window * window)

{
	if(window)
		XLALDestroyREAL8Sequence(window->data);
	XLALFree(window);
}


/*
 * ============================================================================
 *
 *                                REAL4Window
 *
 * ============================================================================
 */

REAL4Window *XLALCreateRectangularREAL4Window(UINT4 length)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateRectangularREAL8Window(length));
}



REAL4Window *XLALCreateHannREAL4Window(UINT4 length)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateHannREAL8Window(length));
}



REAL4Window *XLALCreateWelchREAL4Window(UINT4 length)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateWelchREAL8Window(length));
}



REAL4Window *XLALCreateBartlettREAL4Window(UINT4 length)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateBartlettREAL8Window(length));
}



REAL4Window *XLALCreateParzenREAL4Window(UINT4 length)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateParzenREAL8Window(length));
}



REAL4Window *XLALCreatePapoulisREAL4Window(UINT4 length)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreatePapoulisREAL8Window(length));
}



REAL4Window *XLALCreateHammingREAL4Window(UINT4 length)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateHammingREAL8Window(length));
}



REAL4Window *XLALCreateKaiserREAL4Window(UINT4 length, REAL4 beta)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateKaiserREAL8Window(length, beta));
}



REAL4Window *XLALCreateCreightonREAL4Window(UINT4 length, REAL4 beta)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateCreightonREAL8Window(length, beta));
}



REAL4Window *XLALCreateTukeyREAL4Window(UINT4 length, REAL4 beta)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateTukeyREAL8Window(length, beta));
}



REAL4Window *XLALCreateGaussREAL4Window(UINT4 length, REAL4 beta)

{
	return XLALREAL4Window_from_REAL8Window(XLALCreateGaussREAL8Window(length, beta));
}



void XLALDestroyREAL4Window(REAL4Window * window)

{
	if(window)
		XLALDestroyREAL4Sequence(window->data);
	XLALFree(window);
}


// ---------- some generic window-handling functions to simplify using the above window-functions ----------

typedef enum tagLALWindowType
  {
    LAL_WINDOWTYPE_RECTANGULAR = 0,
    LAL_WINDOWTYPE_HANN,
    LAL_WINDOWTYPE_WELCH,
    LAL_WINDOWTYPE_BARTLETT,
    LAL_WINDOWTYPE_PARZEN,
    LAL_WINDOWTYPE_PAPOULIS,
    LAL_WINDOWTYPE_HAMMING,
    LAL_WINDOWTYPE_KAISER,
    LAL_WINDOWTYPE_CREIGHTON,
    LAL_WINDOWTYPE_TUKEY,
    LAL_WINDOWTYPE_GAUSS,
    LAL_WINDOWTYPE_LAST
  } LALWindowType;

const struct {
  const char *const name;	/**< window name */
  const BOOLEAN hasBeta;	/**< does this window need a 'beta' parameter? */
} AllowedWindows[LAL_WINDOWTYPE_LAST] = {

  [LAL_WINDOWTYPE_RECTANGULAR] 	= { "rectangular",  	0 },
  [LAL_WINDOWTYPE_HANN]		= { "hann",		0 },
  [LAL_WINDOWTYPE_WELCH]	= { "welch",		0 },
  [LAL_WINDOWTYPE_BARTLETT]	= { "bartlett",		0 },
  [LAL_WINDOWTYPE_PARZEN]	= { "parzen",		0 },
  [LAL_WINDOWTYPE_PAPOULIS]	= { "papoulis",		0 },
  [LAL_WINDOWTYPE_HAMMING]	= { "hamming",		1 },
  [LAL_WINDOWTYPE_KAISER]	= { "kaiser",		1 },
  [LAL_WINDOWTYPE_CREIGHTON]	= { "creighton",	1 },
  [LAL_WINDOWTYPE_TUKEY]	= { "tukey",		1 },
  [LAL_WINDOWTYPE_GAUSS]	= { "gauss",		1 },
};

/**
 * Parse window-name string (case-insensitive) into an internal
 * window-type index (>=0, returned), and also check if the user-input 'beta'
 * is valid for given window. Window-types that don't take a beta
 * input parameter need to have beta==0.
 *
 * Returns XLAL_FAILURE=-1 on error
 */
static int
XLALParseWindowNameAndCheckBeta ( const char *windowName,	//< [in] window-name to parse
                                  REAL8 beta			//< [in] beta user-input, checked for validity
                                  )
{
  XLAL_CHECK ( windowName != NULL, XLAL_EINVAL );

  // convert input window-name into lower-case first
  char windowNameLC [ strlen(windowName) + 1 ];
  strcpy ( windowNameLC, windowName );
  XLALStringToLowerCase ( windowNameLC );

  for ( UINT4 i = 0; i < LAL_WINDOWTYPE_LAST; i ++ )
    {
      if ( strcmp ( windowNameLC, AllowedWindows[i].name ) == 0 )
        {
          XLAL_CHECK ( AllowedWindows[i].hasBeta || (beta == 0 ), XLAL_EINVAL, "Invalid non-zero input beta=%g for window '%s'\n", beta, windowName );
          return i;
        }
    } // for i < LAL_WINDOWTYPE_LAST

  // we only come here if no window-name matched
  XLALPrintError ("Invalid Window-name '%s', allowed are (case-insensitive):\n[%s", windowName, AllowedWindows[0].name );
  for ( UINT4 j = 1; j < LAL_WINDOWTYPE_LAST; j++ ) {
    XLALPrintError (", %s", AllowedWindows[j].name );
  }
  XLALPrintError ("]\n");

  XLAL_ERROR ( XLAL_EINVAL );

} // XLALParseWindowNameAndCheckBeta()

/**
 * Generic window-function wrapper, allowing to select a window by its name.
 * windowBeta must be set to '0' for windows without parameter.
 */
REAL8Window *
XLALCreateNamedREAL8Window ( const char *windowName, REAL8 beta, UINT4 length )
{
  XLAL_CHECK_NULL ( length > 0, XLAL_EINVAL );

  int wintype;
  XLAL_CHECK_NULL ( (wintype = XLALParseWindowNameAndCheckBeta ( windowName, beta )) >= 0, XLAL_EFUNC );

  REAL8Window *win = NULL;
  switch ( wintype )
    {
    case LAL_WINDOWTYPE_RECTANGULAR:
      win = XLALCreateRectangularREAL8Window ( length );
      break;
    case LAL_WINDOWTYPE_HANN:
      win = XLALCreateHannREAL8Window ( length );
      break;
    case LAL_WINDOWTYPE_WELCH:
      win = XLALCreateWelchREAL8Window ( length );
      break;
    case LAL_WINDOWTYPE_BARTLETT:
      win = XLALCreateBartlettREAL8Window ( length );
      break;
    case LAL_WINDOWTYPE_PARZEN:
      win = XLALCreateParzenREAL8Window ( length );
      break;
    case LAL_WINDOWTYPE_PAPOULIS:
      win = XLALCreatePapoulisREAL8Window ( length );
      break;
    case LAL_WINDOWTYPE_HAMMING:
      win = XLALCreateHammingREAL8Window ( length );
      break;
    case LAL_WINDOWTYPE_KAISER:
      win = XLALCreateKaiserREAL8Window ( length, beta );
      break;
    case LAL_WINDOWTYPE_CREIGHTON:
      win = XLALCreateCreightonREAL8Window ( length, beta );
      break;
    case LAL_WINDOWTYPE_TUKEY:
      win = XLALCreateTukeyREAL8Window ( length, beta );
      break;
    case LAL_WINDOWTYPE_GAUSS:
      win = XLALCreateGaussREAL8Window ( length, beta );
      break;
    default:
      XLAL_ERROR_NULL ( XLAL_EERR, "Internal ERROR: Invalid window-type '%d', must be within [0, %d]\n", wintype, LAL_WINDOWTYPE_LAST - 1 );
      break;
    } // switch(wintype)

  XLAL_CHECK_NULL (win != NULL, XLAL_EFUNC );

  return win;

} /* XLALCreateNamedREAL8Window() */


REAL4Window *
XLALCreateNamedREAL4Window ( const char *windowName, REAL8 beta, UINT4 length )
{
  return XLALREAL4Window_from_REAL8Window ( XLALCreateNamedREAL8Window ( windowName, beta, length ) );
}
