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

#if !defined(PAGER) || !defined(HAVE_POPEN) || !defined(HAVE_PCLOSE)
#define popen(...) stdout
#define pclose(...)
#endif

int main(int argc, char *argv[])
{
  fitsfile *fptr = 0;   /* FITS file pointer, defined in fitsio.h */
  int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
  int bitpix = 0, naxis = 0, ii = 0, d = 0;
  long naxes[9] = {1,1,1,1,1,1,1,1,1}, fpixel[9] = {1,1,1,1,1,1,1,1,1};
  double *pixels = 0;
  char format[20], hdformat[20];

  if (argc != 2) {
    fprintf(stderr, "Usage:  %s filename[ext][section filter] \n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "List the the pixel values in a FITS array \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example: \n");
    fprintf(stderr, "  %s array.fits                    - list the whole array\n", argv[0]);
    fprintf(stderr, "  %s array.fits[100:110,400:410]   - list a section\n", argv[0]);
    fprintf(stderr, "  %s table.fits[2][bin (x,y) = 32] - list the pixels in\n", argv[0]);
    fprintf(stderr, "         an array constructed from a 2D histogram of X and Y\n");
    fprintf(stderr, "         columns in a table with a binning factor = 32\n");
    return (0);
  }

  FILE *fout = popen(PAGER, "w");
  if (fout == NULL) {
    fprintf(stderr, "Could not execute '%s'\n", PAGER);
    return (1);
  }

  if (!fits_open_file(&fptr, argv[1], READONLY, &status)) {
    if (!fits_get_img_param(fptr, 9, &bitpix, &naxis, naxes, &status)) {
      if (naxis > 9 || naxis == 0) {
        fprintf(stderr, "Error: only 1- to 9-dimensional arrays are supported\n");
      } else {
        /* get memory for 1 row */
        pixels = (double *) malloc(naxes[0] * sizeof(double));

        if (pixels == NULL) {
          fprintf(stderr, "Memory allocation error\n");
          return (1);
        }

        if (bitpix > 0) {  /* set the default output format string */
          strcpy(hdformat, "   %7d");
          strcpy(format,   "   %7.0g");
        } else {
          strcpy(hdformat, "   %15d");
          strcpy(format,   "   %15.5g");
        }

        if (naxis > 2) {  /* label higher dimensions */
          fprintf(fout, "#");
          for (d = naxis - 1; d > 0; d--) {
            fprintf(fout, "%1iD ", d+1);
          }
          fprintf(fout, "\n");
        }

        /* loop over all the rows in the array */
        for (fpixel[8] = 1; fpixel[8] <= naxes[8]; fpixel[8]++) {
          for (fpixel[7] = 1; fpixel[7] <= naxes[7]; fpixel[7]++) {
            for (fpixel[6] = 1; fpixel[6] <= naxes[6]; fpixel[6]++) {
              for (fpixel[5] = 1; fpixel[5] <= naxes[5]; fpixel[5]++) {
                for (fpixel[4] = 1; fpixel[4] <= naxes[4]; fpixel[4]++) {
                  for (fpixel[3] = 1; fpixel[3] <= naxes[3]; fpixel[3]++) {
                    for (fpixel[2] = 1; fpixel[2] <= naxes[2]; fpixel[2]++) {
                      for (fpixel[1] = 1; fpixel[1] <= naxes[1]; fpixel[1]++) {
                        if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL,
                                          pixels, NULL, &status)) { /* read row of pixels */
                          break;  /* jump out of loop on error */
                        }

                        if (naxis > 2) {  /* print higher dimensions */
                          fprintf(fout, " ");
                          for (d = naxis - 1; d > 0; d--) {
                            fprintf(fout, "% 2li ", fpixel[d]);
                          }
                        }
                        for (ii = 0; ii < naxes[0]; ii++) {
                          fprintf(fout, format, pixels[ii]);  /* print each value  */
                        }
                        fprintf(fout, "\n");                    /* terminate line */
                      }
                    }
                  }
                }
              }
            }
          }
        }
        free(pixels);
      }
    }
    fits_close_file(fptr, &status);
  }

  pclose(fout);

  if (status) {
    fits_report_error(stderr, status);  /* print any error message */
  }
  return (status);
}
