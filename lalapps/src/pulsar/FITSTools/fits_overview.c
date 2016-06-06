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
  fitsfile *fptr = 0;         /* FITS file pointer, defined in fitsio.h */
  char keyname[FLEN_KEYWORD], colname[FLEN_VALUE], coltype[FLEN_VALUE];
  int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
  int single = 0, hdupos = 0, hdutype = 0, bitpix = 0, naxis = 0, ncols = 0, ii = 0;
  long naxes[10], nrows = 0;

  if (argc != 2) {
    fprintf(stderr, "Usage:  %s filename[ext] \n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "List the structure of a single extension, or, if ext is \n");
    fprintf(stderr, "not given, list the structure of the entire FITS file.  \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note that it may be necessary to enclose the input file\n");
    fprintf(stderr, "name in single quote characters on the Unix command line.\n");
    return (0);
  }

  FILE *fout = popen(PAGER, "w");
  if (fout == NULL) {
    fprintf(stderr, "Could not execute '%s'\n", PAGER);
    return (1);
  }

  if (!fits_open_file(&fptr, argv[1], READONLY, &status)) {
    fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */

    /* List only a single structure if a specific extension was given */
    if (strchr(argv[1], '[') || strchr(argv[1], '+')) {
      single++;
    }

    for (; !status; hdupos++) { /* Main loop for each HDU */
      fits_get_hdu_type(fptr, &hdutype, &status);  /* Get the HDU type */

      fprintf(fout, "\nHDU #%d  ", hdupos);
      if (hdutype == IMAGE_HDU) { /* primary array or image HDU */
        fits_get_img_param(fptr, 10, &bitpix, &naxis, naxes, &status);

        fprintf(fout, "Array:  NAXIS = %d,  BITPIX = %d\n", naxis, bitpix);
        for (ii = 0; ii < naxis; ii++) {
          fprintf(fout, "   NAXIS%d = %ld\n",ii+1, naxes[ii]);
        }
      } else { /* a table HDU */
        fits_get_num_rows(fptr, &nrows, &status);
        fits_get_num_cols(fptr, &ncols, &status);

        if (hdutype == ASCII_TBL) {
          fprintf(fout, "ASCII Table:  ");
        } else {
          fprintf(fout, "Binary Table:  ");
        }

        fprintf(fout, "%d columns x %ld rows\n", ncols, nrows);
        fprintf(fout, " COL NAME             FORMAT\n");

        for (ii = 1; ii <= ncols; ii++) {
          fits_make_keyn("TTYPE", ii, keyname, &status); /* make keyword */
          fits_read_key(fptr, TSTRING, keyname, colname, NULL, &status);
          fits_make_keyn("TFORM", ii, keyname, &status); /* make keyword */
          fits_read_key(fptr, TSTRING, keyname, coltype, NULL, &status);

          fprintf(fout, " %3d %-16s %-16s\n", ii, colname, coltype);
        }
      }

      if (single) {
        break;  /* quit if only listing a single HDU */
      }

      fits_movrel_hdu(fptr, 1, NULL, &status);  /* try move to next ext */
    }

    if (status == END_OF_FILE) {
      status = 0;  /* Reset normal error */
    }
    fits_close_file(fptr, &status);
  }

  pclose(fout);

  if (status) {
    fits_report_error(stderr, status);  /* print any error message */
  }
  return (status);
}
