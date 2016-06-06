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
  fitsfile *infptr = 0, *outfptr = 0;  /* FITS file pointers */
  int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
  int hdutype = 0, hdunum = 0, ii = 0;

  if (argc != 4) {
    fprintf(stderr, "Usage:  %s infile expression outfile\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Copy selected rows from the input table to the output file\n");
    fprintf(stderr, "based on the input boolean expression.  The expression may \n");
    fprintf(stderr, "be a function of the values in other table columns or header \n");
    fprintf(stderr, "keyword values.  If the expression evaluates to 'true' then \n");
    fprintf(stderr, "that row is copied to the output file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example: \n");
    fprintf(stderr, "1. %s intab.fits+1 'counts > 0' outab.fits\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    copy rows that have a positive 'counts' column value\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "2. %s intab.fits+1 'gtifilter()' outab.fits\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    Select rows which have a Time column value that is\n");
    fprintf(stderr, "    within one of the Good Time Intervals (GTI) which are\n");
    fprintf(stderr, "    defined in a separate GTI extension in the same file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "3. %s intab.fits+1 'regfilter(\"pow.reg\")' outab.fits\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    Select rows which have X,Y column coordinates located\n");
    fprintf(stderr, "    within the spatial region defined in the file named\n");
    fprintf(stderr, "    'pow.reg'.  This is an ASCII text file containing a\n");
    fprintf(stderr, "    list of one or more geometric regions such as circle,\n");
    fprintf(stderr, "    rectangle, annulus, etc.\n");
    return (0);
  }
  if (!fits_open_file(&infptr, argv[1], READONLY, &status)) {
    if (fits_get_hdu_type(infptr, &hdutype,&status) ||
        hdutype==IMAGE_HDU) {
      fprintf(stderr, "Error: input HDU is not a table\n");
    } else {

      fits_get_hdu_num(infptr, &hdunum);  /* save current HDU location */

      if (!fits_create_file(&outfptr, argv[3], &status)) {
        /* copy all the HDUs from the input file to the output file */
        for (ii = 1; !status; ii++) {
          if (!fits_movabs_hdu(infptr, ii, NULL, &status)) {
            fits_copy_hdu(infptr, outfptr, 0, &status);
          }
        }

        if (status == END_OF_FILE) {
          status = 0;  /* reset expected error */
        }

        /* move back to initial position in the file */
        fits_movabs_hdu(outfptr, hdunum, NULL, &status);

        /* argv[2] is the expression */
        /* input and output files are the same, so delete rows that */
        /* do not satisfy the expression */
        fits_select_rows(outfptr, outfptr, argv[2], &status);

        fits_close_file(outfptr, &status);  /* Done */
      }
    }
    fits_close_file(infptr, &status);
  }

  if (status) {
    fits_report_error(stderr, status);  /* print any error message */
  }
  return (status);
}
