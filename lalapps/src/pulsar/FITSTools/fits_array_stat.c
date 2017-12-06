//
// From https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html:
//
// FITS Tools: Handy FITS Utilities that illustrate how to use CFITSIO
// -------------------------------------------------------------------
//
// These are working programs written in ANSI C that illustrate how one can
// easily read, write, and modify FITS files using the CFITSIO library. Most of
// these programs are very short, containing only a few 10s of lines of
// executable code or less, yet they perform quite useful operations on FITS
// files. Copy the programs to your local machine, then compile, and link them
// with the CFITSIO library. A short description of how to use each program can
// be displayed by executing the program without any command line arguments.
//
// You may freely modify, reuse, and redistribute these programs as you wish. It
// is often easier to use one of these programs as a template when writing a new
// program, rather than coding the new program completely from scratch.
//

/**
 * \file
 * \ingroup lalapps_pulsar_FITSTools
 */

#include <config.h>
#include <string.h>
#include <stdio.h>

#if defined(HAVE_LIBCFITSIO)
#include <fitsio.h>
#else
#error CFITSIO library is not available
#endif

int main(int argc, char *argv[])
{
  fitsfile *fptr = 0;  /* FITS file pointer */
  int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
  int hdutype = 0, naxis = 0, ii = 0;
  long naxes[2], totpix = 0, fpixel[2];
  double *pix, sum = 0., meanval = 0., minval = 1.E33, maxval = -1.E33;

  if (argc != 2) {
    fprintf(stderr, "Usage: %s array \n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Compute statistics of pixels in the input array\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  imarith array.fits                    - the whole array\n");
    fprintf(stderr, "  imarith 'array.fits[200:210,300:310]' - array section\n");
    fprintf(stderr, "  imarith 'table.fits+1[bin (X,Y) = 4]' - array constructed\n");
    fprintf(stderr, "     from X and Y columns of a table, with 4-pixel bin size\n");
    return (0);
  }

  if (!fits_open_image(&fptr, argv[1], READONLY, &status)) {
    if (fits_get_hdu_type(fptr, &hdutype, &status) || hdutype != IMAGE_HDU) {
      fprintf(stderr, "Error: this program only works on arrays, not tables\n");
      return (1);
    }

    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, 2, naxes, &status);

    if (status || naxis != 2) {
      fprintf(stderr, "Error: NAXIS = %d.  Only 2-D arrays are supported.\n", naxis);
      return (1);
    }

    pix = (double *) malloc(naxes[0] * sizeof(double)); /* memory for 1 row */

    if (pix == NULL) {
      fprintf(stderr, "Memory allocation error\n");
      return (1);
    }

    totpix = naxes[0] * naxes[1];
    fpixel[0] = 1;  /* read starting with first pixel in each row */

    /* process array one row at a time; increment row # in each loop */
    for (fpixel[1] = 1; fpixel[1] <= naxes[1]; fpixel[1]++) {
      /* give starting pixel coordinate and number of pixels to read */
      if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0],0, pix,0, &status)) {
        break;  /* jump out of loop on error */
      }

      for (ii = 0; ii < naxes[0]; ii++) {
        sum += pix[ii];                      /* accumlate sum */
        if (pix[ii] < minval) {
          minval = pix[ii];  /* find min and  */
        }
        if (pix[ii] > maxval) {
          maxval = pix[ii];  /* max values    */
        }
      }
    }

    free(pix);
    fits_close_file(fptr, &status);
  }

  if (status)  {
    fits_report_error(stderr, status); /* print any error message */
  } else {
    if (totpix > 0) {
      meanval = sum / totpix;
    }

    printf("Statistics of %ld x %ld  array\n",
           naxes[0], naxes[1]);
    printf("  sum of pixels = %g\n", sum);
    printf("  mean value    = %g\n", meanval);
    printf("  minimum value = %g\n", minval);
    printf("  maximum value = %g\n", maxval);
  }

  return (status);
}
