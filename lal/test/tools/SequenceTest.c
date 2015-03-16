/*
*  Copyright (C) 2007 Bernd Machenschalk, Kipp Cannon
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <stdio.h>
#include <stdlib.h>

#include <lal/LALDatatypes.h>
#include <lal/Sequence.h>
#include <lal/LALMalloc.h>

static REAL4Sequence *random_sequence(size_t length)
{
	REAL4Sequence *s = XLALCreateREAL4Sequence(length);

	while(length--)
		s->data[length] = rand() / (float) RAND_MAX;
	return(s);
}


static INT4Sequence *sequential_sequence(size_t length)
{
	INT4Sequence *s = XLALCreateINT4Sequence(length);

	while(length--)
		s->data[length] = length;
	return(s);
}


static int cmp_real4(REAL4 *a, REAL4 *b, size_t n)
{
	REAL4 d;

	while(n--) {
		d = *a - *b;
		if(d < 0.0)
			return(-1);
		if(d > 0.0)
			return(+1);
	}
	return(0);
}


static int cmp(REAL4Sequence *a, REAL4Sequence *b)
{
	int result;

	result = cmp_real4(a->data, b->data, (b->length < a->length) ? b->length : a->length);
	if(result)
		return(result);
	return(a->length - b->length);
}


int main(void)
{
	REAL4Sequence *x, *y;
	INT4Sequence *a;
	int i;

	/*
	 * Destroy
	 */

	/* NULL pointer */
	XLALDestroyREAL4Sequence(NULL);

	/* Incompletely/Incorrectly initialized structure */
	x = XLALCalloc(1, sizeof(*x));
	x->length = 1000;
	x->data = NULL;
	XLALDestroyREAL4Sequence(x);


	/*
	 * Create
	 */

	/* try segfaulting on array access */
	x = XLALCreateREAL4Sequence(1);
	x->data[0] = 1.0;
	if((x->length != 1) || (x->data[0] != 1.0)) {
		fprintf(stderr, "Create test 1 failed\n");
		exit(1);
	}
	XLALDestroyREAL4Sequence(x);

	/*
	 * Cut
	 */

	/* byte-by-byte compare extraction of a piece */
	x = random_sequence(1024);
	y = XLALCutREAL4Sequence(x, 256, 512);
	if(cmp_real4(x->data + 256, y->data, 512)) {
		fprintf(stderr, "Cut test 1 failed\n");
		exit(1);
	}
	XLALDestroyREAL4Sequence(y);

	/* byte-by-byte compare extraction of entire sequence */
	y = XLALCutREAL4Sequence(x, 0, 1024);
	if(cmp_real4(x->data, y->data, 1024)) {
		fprintf(stderr, "Cut test 2 failed\n");
		exit(1);
	}
	XLALDestroyREAL4Sequence(x);
	XLALDestroyREAL4Sequence(y);

	/*
	 * Copy
	 */

	/* byte-by-byte compare copy */
	x = random_sequence(1024);
	y = XLALCopyREAL4Sequence(x);
	if(cmp(x, y)) {
		fprintf(stderr, "Copy test 1 failed\n");
		exit(1);
	}
	XLALDestroyREAL4Sequence(x);
	XLALDestroyREAL4Sequence(y);

	/*
	 * Shift
	 */

	/* test a positive shift */
	x = random_sequence(1024);
	y = XLALCopyREAL4Sequence(x);
	XLALShiftREAL4Sequence(y, 512);
	for(i = 0; i < 512; i++)
		if(y->data[i] != 0.0) {
			fprintf(stderr, "Shift test 1a failed\n");
			exit(1);
		}
	if(cmp_real4(x->data, y->data + 512, 512)) {
		fprintf(stderr, "Shift test 1b failed\n");
		exit(1);
	}

	/* test a subsequent negative shift */
	XLALShiftREAL4Sequence(y, -768);
	for(i = 256; i < 1024; i++)
		if(y->data[i] != 0.0) {
			fprintf(stderr, "Shift test 2a failed\n");
			exit(1);
		}
	if(cmp_real4(x->data + 256, y->data, 256)) {
		fprintf(stderr, "Shift test 2b failed\n");
		exit(1);
	}
	XLALDestroyREAL4Sequence(x);
	XLALDestroyREAL4Sequence(y);

	/*
	 * Resize
	 */

	/* test resize to subset */
	a = sequential_sequence(1024);
	XLALResizeINT4Sequence(a, 256, 512);
	for(i = 0; i < (int) a->length; i++)
		if(a->data[i] != i + 256) {
			fprintf(stderr, "Resize test 1a failed\n");
			exit(1);
		}
	if(a->length != 512) {
		fprintf(stderr, "Resize test 1b failed\n");
		exit(1);
	}
	XLALDestroyINT4Sequence(a);

	/* test resize to superset */
	a = sequential_sequence(16);
	for(i = 0; i < (int) a->length; i++)
		fprintf(stdout, "%d: %d\n", i, a->data[i]);
	fprintf(stdout, "\n");
	XLALResizeINT4Sequence(a, -8, 32);
	for(i = 0; i < (int) a->length; i++)
		fprintf(stdout, "%d: %d\n", i, a->data[i]);
	if(a->length != 32) {
		fprintf(stderr, "Resize test 2a failed\n");
		exit(1);
	}
	for(i = 0; i < 8; i++)
		if(a->data[i] != 0) {
			fprintf(stderr, "Resize test 2b failed\n");
			exit(1);
		}
	for(; i < 24; i++)
		if(a->data[i] != i - 8) {
			fprintf(stderr, "Resize test 2c failed\n");
			exit(1);
		}
	for(; i < 32; i++)
		if(a->data[i] != 0) {
			fprintf(stderr, "Resize test 2d failed\n");
			exit(1);
		}
	XLALDestroyINT4Sequence(a);

	/*
	 * Sum
	 */

	a = sequential_sequence(1024);
	if(XLALINT4SequenceSum(a, 0, a->length) != (1023 + 0) * 1024 / 2) {
		fprintf(stderr, "Sum test 1 failed\n");
		exit(1);
	}
	XLALDestroyINT4Sequence(a);

	/*
	 * Sum squares
	 */

	a = sequential_sequence(1024);
	if(XLALINT4SequenceSumSquares(a, 0, a->length) != 1023 * (1023 + 1) * (2 * 1023 + 1) / 6) {
		fprintf(stderr, "Sum squares test 1 failed\n");
		exit(1);
	}
	XLALDestroyINT4Sequence(a);

	/*
	 * Success
	 */

	exit(0);
}
