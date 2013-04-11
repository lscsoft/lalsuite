/*
 * Copyright (C) 2007 Bruce Allen, Duncan Brown, Jolien Creighton, Kipp
 * Cannon, Teviet Creighton
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
 * 02111-1307 USA
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lal/LALDatatypes.h>
#include <lal/Window.h>
#include <lal/XLALError.h>
#include <lal/LALMalloc.h>

#define NWINDOWS 11

int lalDebugLevel = 0;

const char *names[] = {
	"rectangular",
	"Hann",
	"Welch",
	"Bartlett",
	"Parzen",
	"Papoulis",
	"Hamming",
	"Kaiser",
	"Creighton",
	"Tukey",
	"Gauss"
};

static int create_single_windows(REAL4Window **windows, int length, double kaiser_beta, double creighton_beta, double tukey_beta, double gauss_beta)
{

  for ( UINT4 i = 0; i < NWINDOWS; i ++ )
    {
      REAL8 beta;

      if ( !strcmp ( names[i], "Kaiser" ) ) {
        beta = kaiser_beta;
      } else if ( !strcmp ( names[i], "Creighton" ) ) {
        beta = creighton_beta;
      } else if ( !strcmp ( names[i], "Tukey" ) ) {
        beta = tukey_beta;
      } else if ( !strcmp ( names[i], "Gauss" ) ) {
        beta = gauss_beta;
      } else {
        beta = 0;
      }

      XLAL_CHECK ( (windows[i] = XLALCreateNamedREAL4Window ( names[i], beta, length )) != NULL, XLAL_EFUNC );

    } // for i < NWINDOWS

  return XLAL_SUCCESS;

} // create_single_windows()


static void free_single_windows(REAL4Window **windows)
{
	int i;
	for(i = 0; i < NWINDOWS; i++)
		XLALDestroyREAL4Window(windows[i]);
}


static int create_double_windows(REAL8Window **windows, int length, double kaiser_beta, double creighton_beta, double tukey_beta, double gauss_beta)
{

  for ( UINT4 i = 0; i < NWINDOWS; i ++ )
    {
      REAL8 beta;

      if ( !strcmp ( names[i], "Kaiser" ) ) {
        beta = kaiser_beta;
      } else if ( !strcmp ( names[i], "Creighton" ) ) {
        beta = creighton_beta;
      } else if ( !strcmp ( names[i], "Tukey" ) ) {
        beta = tukey_beta;
      } else if ( !strcmp ( names[i], "Gauss" ) ) {
        beta = gauss_beta;
      } else {
        beta = 0;
      }

      XLAL_CHECK ( (windows[i] = XLALCreateNamedREAL8Window ( names[i], beta, length )) != NULL, XLAL_EFUNC );

    } // for i < NWINDOWS

  return XLAL_SUCCESS;

} // create_double_windows()


static void free_double_windows(REAL8Window **windows)
{
	int i;
	for(i = 0; i < NWINDOWS; i++)
		XLALDestroyREAL8Window(windows[i]);
}


static double fractional_difference(double a, double b)
{
	if(a != 0)
		/* plan A */
		return fabs((a - b) / a);
	if(b != 0)
		/* plan B */
		return fabs((a - b) / b);
	/* both are 0 */
	return 0;
}


/*
 * Sum-of-squares test.
 */


static int _test_sum_of_squares(const double *correct, int length, double kaiser_beta, double creighton_beta, double tukey_beta, double gauss_beta)
{
	const double max_error = 1e-12;
	REAL4Window *windows1[NWINDOWS];
	REAL8Window *windows2[NWINDOWS];
	int i;

	XLAL_CHECK ( create_single_windows(windows1, length, kaiser_beta, creighton_beta, tukey_beta, gauss_beta) == XLAL_SUCCESS, XLAL_EFUNC );

	XLAL_CHECK ( create_double_windows(windows2, length, kaiser_beta, creighton_beta, tukey_beta, gauss_beta) == XLAL_SUCCESS, XLAL_EFUNC );

	for(i = 0; i < NWINDOWS; i++) {
		if(fractional_difference(windows1[i]->sumofsquares, correct[i]) > max_error) {
                  XLAL_ERROR (XLAL_EFAILED, "error: single-precision %d-sample %s window fails sum-of-squares test:  expected %.17g, got %.17g\n",
                              length, names[i], correct[i], windows1[i]->sumofsquares);
		}
		if(fractional_difference(windows2[i]->sumofsquares, correct[i]) > max_error) {
                  XLAL_ERROR (XLAL_EFAILED, "error: double-precision %d-sample %s window fails sum-of-squares test:  expected %.17g, got %.17g\n", length, names[i], correct[i], windows1[i]->sumofsquares);
		}
	}

	free_single_windows(windows1);
	free_double_windows(windows2);

	return XLAL_SUCCESS;

} // _test_sum_of_squares()


static int test_sum_of_squares(void)
{
	double correct_1024[] = {
		1024.0,			/* rectangle */
		383.625,		/* Hann */
		545.6,			/* Welch */
		340.9996741609645,	/* Bartlett */
		275.84464285585898,	/* Parzen */
		300.06446358192244,	/* Papoulis */
		406.5466,		/* Hamming */
		375.17819205246843,	/* Kaiser */
		392.64506106773848,	/* Creighton */
		703.625,		/* Tukey */
		451.20289927038817	/* Gauss */
	};
	double correct_1025[] = {
		1025.0,			/* rectangle */
		384,			/* Hann */
		546.0 + 2.0 / 15.0,	/* Welch */
		341.333984375,		/* Bartlett */
		276.1142857152779,	/* Parzen */
		300.35778172967611,	/* Papoulis */
		406.944,		/* Hamming */
		375.544934942032,	/* Kaiser */
		393.028878331734330,	/* Creighton */
		704,			/* Tukey */
		451.64394001239367	/* Gauss */
	};

	XLAL_CHECK ( _test_sum_of_squares(correct_1024, 1024, 6, 2, 0.5, 2) == XLAL_SUCCESS, XLAL_EFUNC );

	XLAL_CHECK ( _test_sum_of_squares(correct_1025, 1025, 6, 2, 0.5, 2) == XLAL_SUCCESS, XLAL_EFUNC );

	return XLAL_SUCCESS;

} // test_sum_of_squares()


/*
 * Test end- and mid-points.
 */


static int _test_end_and_midpoints(int length, double kaiser_beta, double creighton_beta, double tukey_beta, double gauss_beta)
{
	const double max_error = 1e-16;
	double correct_end[] = {
		1,			/* rectangle */
		0,			/* Hann */
		0,			/* Welch */
		0,			/* Bartlett */
		0,			/* Parzen */
		0,			/* Papoulis */
		0.08,			/* Hamming */
		0.014873337104763204,	/* Kaiser (to be adjusted below) */
		0,			/* Creighton (to be adjusted below) */
		0,			/* Tukey (to be adjusted below) */
		0			/* Gauss (to be adjusted below) */
	};
	double correct_mid[] = {
		1,	/* rectangle */
		1,	/* Hann */
		1,	/* Welch */
		1,	/* Bartlett */
		1,	/* Parzen */
		1,	/* Papoulis */
		1,	/* Hamming */
		1,	/* Kaiser */
		1,	/* Creighton */
		1,	/* Tukey */
		1	/* Gauss */
	};
	REAL4Window *windows1[NWINDOWS];
	REAL8Window *windows2[NWINDOWS];
	int i;

	/* set end value of Kaiser window */
	correct_end[7] = kaiser_beta == 0 ? 1 : kaiser_beta == HUGE_VAL ? 0 : correct_end[7];

	/* set end value of Creighton window */
	correct_end[8] = creighton_beta == 0 ? 1 : correct_end[8];

	/* set end value of Tukey window */
	correct_end[9] = tukey_beta == 0 ? 1 : correct_end[9];

	/* set end value of Gauss window */
	correct_end[10] = exp(-0.5 * gauss_beta * gauss_beta);

	XLAL_CHECK ( create_single_windows(windows1, length, kaiser_beta, creighton_beta, tukey_beta, gauss_beta) == XLAL_SUCCESS, XLAL_EFUNC );
	XLAL_CHECK ( create_double_windows(windows2, length, kaiser_beta, creighton_beta, tukey_beta, gauss_beta) == XLAL_SUCCESS, XLAL_EFUNC );

	for(i = 0; i < NWINDOWS; i++) {
		/* if length < 2, then there are no end samples */
		if(length >= 2) {
			if(fabs(windows1[i]->data->data[0] - (float) correct_end[i]) > max_error) {
                          XLAL_ERROR ( XLAL_EFAILED, "error: single-precision %d-sample %s window fails end-point test:  expected %.17g, got %.17g\n", length, names[i], correct_end[i], windows1[i]->data->data[0]);
			}
			if(fabs(windows2[i]->data->data[0] - correct_end[i]) > max_error) {
                          XLAL_ERROR ( XLAL_EFAILED, "error: double-precision %d-sample %s window fails end-point test:  expected %.17g, got %.17g\n", length, names[i], correct_end[i], windows2[i]->data->data[0]);
			}
			if(fabs(windows1[i]->data->data[length - 1] - (float) correct_end[i]) > max_error) {
                          XLAL_ERROR ( XLAL_EFAILED, "error: single-precision %d-sample %s window fails end-point test:  expected %.17g, got %.17g\n", length, names[i], correct_end[i], windows1[i]->data->data[length - 1]);
			}
			if(fabs(windows2[i]->data->data[length - 1] - correct_end[i]) > max_error) {
                          XLAL_ERROR ( XLAL_EFAILED, "error: double-precision %d-sample %s window fails end-point test:  expected %.17g, got %.17g\n", length, names[i], correct_end[i], windows2[i]->data->data[length - 1]);
			}
		}
		/* even-lengthed windows have no middle sample */
		if(length & 1) {
			if(windows1[i]->data->data[length / 2] != (float) correct_mid[i]) {
                          XLAL_ERROR ( XLAL_EFAILED, "error: single-precision %d-sample %s window fails mid-point test:  expected %.17g, got %.17g\n", length, names[i], correct_mid[i], windows1[i]->data->data[length / 2]);
			}
			if(windows2[i]->data->data[length / 2] != correct_mid[i]) {
                          XLAL_ERROR ( XLAL_EFAILED,  "error: double-precision %d-sample %s window fails mid-point test:  expected %.17g, got %.17g\n", length, names[i], correct_mid[i], windows1[i]->data->data[length / 2]);
			}
		}
	}

	free_single_windows(windows1);
	free_double_windows(windows2);

	return XLAL_SUCCESS;
}


static int test_end_and_midpoints(void)
{
	int fail = 0;

	if(_test_end_and_midpoints(1025, HUGE_VAL, HUGE_VAL, 1, HUGE_VAL))
		fail = 1;
	if(_test_end_and_midpoints(1025, 6, 2, 0.5, 2))
		fail = 1;
	if(_test_end_and_midpoints(1025, 0, 0, 0, 0))
		fail = 1;
	if(_test_end_and_midpoints(1024, HUGE_VAL, HUGE_VAL, 1, HUGE_VAL))
		fail = 1;
	if(_test_end_and_midpoints(1024, 6, 2, 0.5, 2))
		fail = 1;
	if(_test_end_and_midpoints(1024, 0, 0, 0, 0))
		fail = 1;
	if(_test_end_and_midpoints(3, HUGE_VAL, HUGE_VAL, 1, HUGE_VAL))
		fail = 1;
	if(_test_end_and_midpoints(3, 6, 2, 0.5, 2))
		fail = 1;
	if(_test_end_and_midpoints(3, 0, 0, 0, 0))
		fail = 1;
	if(_test_end_and_midpoints(1, HUGE_VAL, HUGE_VAL, 1, HUGE_VAL))
		fail = 1;
	if(_test_end_and_midpoints(1, 6, 2, 0.5, 2))
		fail = 1;
	if(_test_end_and_midpoints(1, 0, 0, 0, 0))
		fail = 1;

	return fail;
}


/*
 * Input parameter safety
 */


static int test_parameter_safety(void)
{
	REAL4Window *window1;
	REAL8Window *window2;
	int fail = 0;

	window1 = XLALCreateKaiserREAL4Window(10, -1);
	if(window1) {
		fprintf(stderr, "error: single-precision Kaiser window accepted out-of-range parameter\n");
		XLALDestroyREAL4Window(window1);
		fail = 1;
	}

	window2 = XLALCreateKaiserREAL8Window(10, -1);
	if(window1) {
		fprintf(stderr, "error: double-precision Kaiser window accepted out-of-range parameter\n");
		XLALDestroyREAL8Window(window2);
		fail = 1;
	}

	window1 = XLALCreateCreightonREAL4Window(10, -1);
	if(window1) {
		fprintf(stderr, "error: single-precision Creighton window accepted out-of-range parameter\n");
		XLALDestroyREAL4Window(window1);
		fail = 1;
	}

	window2 = XLALCreateCreightonREAL8Window(10, -1);
	if(window1) {
		fprintf(stderr, "error: double-precision Creighton window accepted out-of-range parameter\n");
		XLALDestroyREAL8Window(window2);
		fail = 1;
	}

	window1 = XLALCreateTukeyREAL4Window(10, -1);
	if(window1) {
		fprintf(stderr, "error: single-precision Tukey window accepted out-of-range parameter\n");
		XLALDestroyREAL4Window(window1);
		fail = 1;
	}

	window1 = XLALCreateTukeyREAL4Window(10, 2);
	if(window1) {
		fprintf(stderr, "error: single-precision Tukey window accepted out-of-range parameter\n");
		XLALDestroyREAL4Window(window1);
		fail = 1;
	}

	window2 = XLALCreateTukeyREAL8Window(10, -1);
	if(window1) {
		fprintf(stderr, "error: double-precision Tukey window accepted out-of-range parameter\n");
		XLALDestroyREAL8Window(window2);
		fail = 1;
	}

	window2 = XLALCreateTukeyREAL8Window(10, 2);
	if(window1) {
		fprintf(stderr, "error: double-precision Tukey window accepted out-of-range parameter\n");
		XLALDestroyREAL8Window(window2);
		fail = 1;
	}

	window1 = XLALCreateGaussREAL4Window(10, -1);
	if(window1) {
		fprintf(stderr, "error: single-precision Gauss window accepted out-of-range parameter\n");
		XLALDestroyREAL4Window(window1);
		fail = 1;
	}

	window2 = XLALCreateGaussREAL8Window(10, -1);
	if(window1) {
		fprintf(stderr, "error: double-precision Gauss window accepted out-of-range parameter\n");
		XLALDestroyREAL8Window(window2);
		fail = 1;
	}

	return fail;
}


/*
 * Display sample windows.
 */


static void _display(int n, double kaiser_beta, double creighton_beta, double tukey_beta, double gauss_beta)
{
	REAL8Window *rectangle = XLALCreateRectangularREAL8Window(n);
	REAL8Window *hann = XLALCreateHannREAL8Window(n);
	REAL8Window *welch = XLALCreateWelchREAL8Window(n);
	REAL8Window *bartlett = XLALCreateBartlettREAL8Window(n);
	REAL8Window *parzen = XLALCreateParzenREAL8Window(n);
	REAL8Window *papoulis = XLALCreatePapoulisREAL8Window(n);
	REAL8Window *hamming = XLALCreateHammingREAL8Window(n);
	REAL8Window *kaiser = XLALCreateKaiserREAL8Window(n, kaiser_beta);
	REAL8Window *creighton = XLALCreateCreightonREAL8Window(n, creighton_beta);
	REAL8Window *tukey = XLALCreateTukeyREAL8Window(n, tukey_beta);
	REAL8Window *gauss = XLALCreateGaussREAL8Window(n, gauss_beta);
	int i;

	printf("n = %d\n", n);
	printf("kaiser beta = %g\n", kaiser_beta);
	printf("creighton beta = %g\n", creighton_beta);
	printf("tukey beta = %g\n", tukey_beta);
	printf("gaussian beta = %g\n", gauss_beta);

	printf("  rect     hann     welch  bartlett  parzen  papoulis  hamming  kaiser   creight   tukey    gauss\n");
	for(i = 0; i < n; i++) {
		printf("%8.6f", rectangle->data->data[i]);
		printf(" %8.6f", hann->data->data[i]);
		printf(" %8.6f", welch->data->data[i]);
		printf(" %8.6f", bartlett->data->data[i]);
		printf(" %8.6f", parzen->data->data[i]);
		printf(" %8.6f", papoulis->data->data[i]);
		printf(" %8.6f", hamming->data->data[i]);
		printf(" %8.6f", kaiser->data->data[i]);
		printf(" %8.6f", creighton->data->data[i]);
		printf(" %8.6f", tukey->data->data[i]);
		printf(" %8.6f", gauss->data->data[i]);
		printf("\n");
	}
	printf("\n");

	XLALDestroyREAL8Window(rectangle);
	XLALDestroyREAL8Window(hann);
	XLALDestroyREAL8Window(welch);
	XLALDestroyREAL8Window(bartlett);
	XLALDestroyREAL8Window(parzen);
	XLALDestroyREAL8Window(papoulis);
	XLALDestroyREAL8Window(hamming);
	XLALDestroyREAL8Window(kaiser);
	XLALDestroyREAL8Window(creighton);
	XLALDestroyREAL8Window(tukey);
	XLALDestroyREAL8Window(gauss);
}


static void display(void)
{
	_display(14, 0, 0, 0, 0);
	_display(15, 0, 0, 0, 0);
	_display(14, 6, 2, 0.5, 2);
	_display(15, 6, 2, 0.5, 2);
	_display(14, HUGE_VAL, HUGE_VAL, 1, HUGE_VAL);
	_display(15, HUGE_VAL, HUGE_VAL, 1, HUGE_VAL);
	_display(5, 6, 2, 0.5, 2);
	_display(4, 6, 2, 0.5, 2);
	_display(3, 6, 2, 0.5, 2);
	_display(2, 6, 2, 0.5, 2);
	_display(1, 0, 0, 0, 0);
	_display(1, 6, 2, 0.5, 2);
	_display(1, HUGE_VAL, HUGE_VAL, 1.0, HUGE_VAL);
}


/*
 * Entry point.
 */


int main(void)
{
	int fail = 0;
	char *hosttype = getenv("hosttype");

	/* DEC Alpha's FPU is not capable of computing some of these window
	 * functions correctly.  This is safe because in the cases where it
	 * fails it raises SIGFPE and the program crashes.  I can't be
	 * bothered to code up the work-arounds needed to get the windows
	 * working on that platform. */

	if(!strcmp(hosttype ? hosttype : "", "alpha")) {
		fprintf(stderr, "Window functions not computable on DEC Alpha, tests skipped!  Set environment variable \"hosttype\" to something other than \"alpha\" to force tests\n");
		exit(77);
	}

	/* Numerical tests:  assume that if the end-points, mid-points, and
	 * sum-of-squares are all as expected, then the window functions
	 * are correct */

	if(test_end_and_midpoints())
		fail = 1;
	if(test_sum_of_squares())
		fail = 1;

	/* Test parameter safety */

	if(test_parameter_safety())
		fail = 1;

	/* Verbosity */

	display();

        LALCheckMemoryLeaks();

	return fail;
}
