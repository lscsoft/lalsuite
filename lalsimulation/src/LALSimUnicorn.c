#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_rng.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALSimBurst.h>

/* this is the xpm image of the unicorn! */
#define char const char /* to silence compiler warnings ... */
#include "unicorn.xpm"
#undef char

/* RGB to luminance */
static double RGB2Y(double R, double G, double B)
{
	double Y;
	Y = 0.2126*R + 0.7152*G + 0.0722*B;
	return Y;
}

/* convert radix 256 string to unsigned long */
static unsigned long radix256toul(const char *s, size_t n)
{
	unsigned long ul = 0;
	memcpy(&ul, s, n < sizeof(ul) ? n : sizeof(ul));
	return ul;
}

/* converts a RGB color string of the form #xxxxxx to a luminance value */
static double strtolum(const char *s)
{
	int R, G, B, max;
	int width;
	char fmt[64];
	if (s[0] != '#') { /* not a rgb color! */
		return -1.0;
	}
	width = strlen(s + 1);
	width /= 3;
	snprintf(fmt, sizeof(fmt), "#%%%dx%%%dx%%%dx", width, width, width);
	sscanf(s, fmt, &R, &G, &B);
	for (max = 16; width > 1; --width)
		max *= 16;
	return RGB2Y((double)R/(double)max, (double)G/(double)max, (double)B/(double)max);
}

/* routine to parse a xmp file */
static REAL8Array *xpmtoarr(const char *xpm[])
{
	REAL8Array *arr; /* the image stored as an array of luminance values */
	double *Y; /* the luminance values of the tagged colors */
	size_t nY;
	int cols, rows, colors, chars_per_pixel;
	int color, row;
	int c;

	/* parse first line of xpm */
	sscanf(xpm[0], "%d %d %d %d", &cols, &rows, &colors, &chars_per_pixel);
	if (cols < 1 || rows < 1 || colors < 1 || chars_per_pixel < 1)
		return NULL;

	/* allocate memory for lookup table of color luminance values */
	for (c = 1, nY = 256; c < chars_per_pixel; ++c)
		nY *= 256;
	Y = LALMalloc(nY * sizeof(*Y));
	if (!Y)
		return NULL;

	/* parse lines of xpm to get color luminance values */
	for (color = 0; color < colors; ++color) {
		const char *s = xpm[1 + color];
		size_t i;
		i = radix256toul(s, chars_per_pixel);
		s = strrchr(s, '#');
		Y[i] = strtolum(s);
		if (Y[i] < 0.0 || Y[i] > 1.0) {
			/* an error occurred parsing this string */
			LALFree(Y);
			return NULL;
		}
	}

	/* scan the image and store luminance values in the array */
	arr = XLALCreateREAL8ArrayL(2, rows, cols);
	if (!arr) {
		LALFree(Y);
		return NULL;
	}
	for (row = 0; row < rows; ++row) {
		const char *s = xpm[1 + colors + row];
		int col;
		for (col = 0; col < cols; ++col) {
			size_t offset = row * cols + col;
			size_t i;
			i = radix256toul(s, chars_per_pixel);
			arr->data[offset] = Y[i];
			s += chars_per_pixel;
		}
	}

	LALFree(Y);
	return arr;
}

#if 0 /* for debugging purposes ... */
static void arr2asc(REAL8Array *arr)
{
	size_t row, rows, col, cols;
	rows = arr->dimLength->data[0];
	cols = arr->dimLength->data[1];
	for (row = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col) {
			int offset = row * cols + col;
			char c;
			if (arr->data[offset] < 0.3)
				c = '.';
			else if (arr->data[offset] < 0.6)
				c = 'o';
			else
				c = 'O';
			printf("%c", c);
		}
		printf("\n");
	}
	return;
}
#endif

/**
 * Generates a time-frequency unicorn signal.
 *
 * \warning Unicorns are rare and are hard to catch.
 */
int XLALSimUnicorn(
	REAL8TimeSeries **hplus,	/**< plus-polarization waveform [returned] */
	REAL8TimeSeries **hcross, 	/**< cross-polarization waveform [returned] */
	double f_min,			/**< minimum frequency of the unicorn (Hz) */
	double f_max,			/**< maximum frequency of the unicorn (Hz) */
	double V,			/**< time-frequency volume of image pixels */
	double hrss,			/**< root-sum-squared value of the waveform (s) */
	double deltaT,			/**< sampleing interval (s) */
	gsl_rng *rng			/**< random number generator */
)
{
	REAL8Array *arr;
	size_t rows;
	double dt, df, fstart;

	/* check sanity of input parameters */
	XLAL_CHECK(V >= LAL_2_PI, XLAL_EINVAL, "Time-frequency volume must be greater than 2/pi");
	XLAL_CHECK(f_max > f_min, XLAL_EINVAL, "Maximum frequency must be greater than minimum frequency");
	XLAL_CHECK(f_max <= 0.5 / deltaT , XLAL_EINVAL, "Maximum frequency must be less than Nyquist frequency");

	/* create image array from xpm */
	arr = xpmtoarr(unicorn);
	if (!arr)
		XLAL_ERROR(XLAL_EDATA, "Could not parse .xpm image file");

	/* generate waveform from image */
	rows = arr->dimLength->data[0];
	df = (f_max - f_min) / rows;
	dt = V / df;
	fstart = f_min + 0.5 * df;
	if (XLALSimBurstImg(hplus, hcross, arr, dt, df, fstart, hrss, deltaT, rng) < 0)
		XLAL_ERROR(XLAL_EFUNC);
	XLALDestroyREAL8Array(arr);

	return 0;
}
