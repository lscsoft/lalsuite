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
  int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
  int hdutype = 0, hdunum = 0, ii = 0;

  if (argc != 5) {
    fprintf(stderr, "Usage:  %s infile expression colname outfile\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Compute new values for the specified table column using the\n");
    fprintf(stderr, "input arithmetic expression which may be a function of the \n");
    fprintf(stderr, "values in other table columns. The input file is first copied\n");
    fprintf(stderr, "to the output file, then the output file is updated with the\n");
    fprintf(stderr, "new column values.  If the column doesn't already exist,\n");
    fprintf(stderr, "then a new column will be appended to the table.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example: \n");
    fprintf(stderr, "1. %s intab.fits+1 'counts/#exposure' rate outab.fits\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    Calculate the 'rate' column from the values in the\n");
    fprintf(stderr, "    'counts' column and the 'exposure' keyword.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "2. %s intab.fits+1 'sqrt(X**2 + Y**2)' Radius outab.fits\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    Calculate the 'Radius' column from the 'X' and 'Y' cols.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "3. %s intab.fits+1 '(rate{-1}+rate+rate{+1})/3' rate3 outab.fits\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    Calculate the running mean of the rate column by \n");
    fprintf(stderr, "    averaging the values in the previous row, the current\n");
    fprintf(stderr, "    row, and the next row\n");
    return (0);
  }
  if (!fits_open_file(&infptr, argv[1], READONLY, &status)) {
    if (fits_get_hdu_type(infptr, &hdutype,&status) ||
        hdutype==IMAGE_HDU) {
      fprintf(stderr, "Error: input HDU is not a table\n");
    } else {

      fits_get_hdu_num(infptr, &hdunum);  /* save current HDU location */

      if (!fits_create_file(&outfptr, argv[4], &status)) {
        /* copy all the HDUs from the input file to the output file */
        for (ii = 1; !status; ii++) {
          fits_movabs_hdu(infptr, ii, NULL, &status);
          fits_copy_hdu(infptr, outfptr, 0, &status);
        }

        if (status == END_OF_FILE) {
          status = 0;  /* reset expected error */
        }

        /* move back to initial position in the file */
        fits_movabs_hdu(outfptr, hdunum, NULL, &status);

        /* argv[2] is the expression, and argv[3] is the col name */
        fits_calculator(outfptr, argv[2], outfptr, argv[3],
                        NULL, &status);

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
