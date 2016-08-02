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
  fitsfile *fptr = 0;      /* FITS file pointer, defined in fitsio.h */
  char *val = 0, value[1000], nullstr[]="NAN";
  char keyword[FLEN_KEYWORD], colname[1000][FLEN_VALUE];
  int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
  int hdunum = 0, hdutype = 0, ncols = 0, ii = 0, anynul = 0, dispwidth[1000];
  long jj = 0, nrows = 0;

  if (argc != 2) {
    fprintf(stderr, "Usage:  %s filename[ext][col filter][row filter] \n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "List the contents of a FITS table \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  %s tab.fits[GTI]           - list the GTI extension\n", argv[0]);
    fprintf(stderr, "  %s tab.fits[1][#row < 101] - list first 100 rows\n", argv[0]);
    fprintf(stderr, "  %s tab.fits[1][col X;Y]    - list X and Y cols only\n", argv[0]);
    fprintf(stderr, "  %s tab.fits[1][col -PI]    - list all but the PI col\n", argv[0]);
    fprintf(stderr, "  %s tab.fits[1][col -PI][#row < 101]  - combined case\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Display formats can be modified with the TDISPn keywords.\n");
    return (0);
  }

  FILE *fout = popen(PAGER, "w");
  if (fout == NULL) {
    fprintf(stderr, "Could not execute '%s'\n", PAGER);
    return (1);
  }

  if (!fits_open_file(&fptr, argv[1], READONLY, &status)) {
    if (fits_get_hdu_num(fptr, &hdunum) == 1)
      /* This is the primary array;  try to move to the */
      /* first extension and see if it is a table */
    {
      fits_movabs_hdu(fptr, 2, &hdutype, &status);
    } else {
      fits_get_hdu_type(fptr, &hdutype, &status);  /* Get the HDU type */
    }

    if (hdutype == IMAGE_HDU) {
      fprintf(stderr, "Error: this program only displays tables, not images\n");
    } else {
      fits_get_num_rows(fptr, &nrows, &status);
      fits_get_num_cols(fptr, &ncols, &status);

      for (ii = 1; ii <= ncols; ii++) {
        fits_make_keyn("TTYPE", ii, keyword, &status);
        fits_read_key(fptr, TSTRING, keyword, colname[ii], NULL, &status);
        fits_get_col_display_width(fptr, ii, &dispwidth[ii], &status);
        if (dispwidth[ii] < (int)strlen(colname[ii])) {
          dispwidth[ii] = (int)strlen(colname[ii]);
        }
      }

      /* print column names as column headers */
      fprintf(fout, "##\n## ");
      for (ii = 1; ii <= ncols; ii++) {
        fprintf(fout, "%*s ",dispwidth[ii], colname[ii]);
      }
      fprintf(fout, "\n");  /* terminate header line */

      /* print each column, row by row (there are faster ways to do this) */
      val = value;
      for (jj = 1; jj <= nrows && !status; jj++) {
        fprintf(fout, "   ");
        for (ii = 1; ii <= ncols; ii++) {
          /* read value as a string, regardless of intrinsic datatype */
          if (fits_read_col_str(fptr,ii,jj, 1, 1, nullstr,
                                &val, &anynul, &status)) {
            break;  /* jump out of loop on error */
          }

          fprintf(fout, "%*s ",dispwidth[ii], value);
        }
        fprintf(fout, "\n");
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
