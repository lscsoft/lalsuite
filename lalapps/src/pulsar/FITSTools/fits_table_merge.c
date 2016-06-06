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
  int icol = 0, incols = 0, outcols = 0, intype = 0, outtype = 0, check = 1;
  long inrep = 0, outrep = 0, width = 0, inrows = 0, outrows = 0, ii = 0, jj = 0;
  unsigned char *buffer = 0;

  if (argc != 3) {
    fprintf(stderr, "Usage:  %s infile1[ext][filter] outfile[ext]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Merge 2 tables by copying all the rows from the 1st table\n");
    fprintf(stderr, "into the 2nd table.  The  2 tables must have identical\n");
    fprintf(stderr, "structure, with the same number of columns with the same\n");
    fprintf(stderr, "datatypes.  This program modifies the output file in place,\n");
    fprintf(stderr, "rather than creating a whole new output file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "1. %s intab.fit+1 outtab.fit+2\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    merge the table in the 1st extension of intab.fit with\n");
    fprintf(stderr, "    the table in the 2nd extension of outtab.fit.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "2. %s 'intab.fit+1[PI > 45]' outab.fits+2\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    Same as the 1st example, except only rows that have a PI\n");
    fprintf(stderr, "    column value > 45 will be merged into the output table.\n");
    fprintf(stderr, "\n");
    return (0);
  }

  /* open both input and output files and perform validity checks */
  if (fits_open_file(&infptr,  argv[1], READONLY,  &status) ||
      fits_open_file(&outfptr, argv[2], READWRITE, &status)) {
    fprintf(stderr, " Couldn't open both files\n");
  }

  else if (fits_get_hdu_type(infptr,  &intype,  &status) ||
           fits_get_hdu_type(outfptr, &outtype, &status)) {
    fprintf(stderr, "couldn't get the type of HDU for the files\n");
  }

  else if (intype == IMAGE_HDU) {
    fprintf(stderr, "The input HDU is an image, not a table\n");
  }

  else if (outtype == IMAGE_HDU) {
    fprintf(stderr, "The output HDU is an image, not a table\n");
  }

  else if (outtype != intype) {
    fprintf(stderr, "Input and output HDUs are not the same type of table.\n");
  }

  else if (fits_get_num_cols(infptr,  &incols,  &status) ||
           fits_get_num_cols(outfptr, &outcols, &status)) {
    fprintf(stderr, "Couldn't get number of columns in the tables\n");
  }

  else if (incols != outcols) {
    fprintf(stderr, "Input and output HDUs don't have same # of columns.\n");
  }

  else if (fits_read_key(infptr, TLONG, "NAXIS1", &width, NULL, &status)) {
    fprintf(stderr, "Couldn't get width of input table\n");
  }

  else if (!(buffer = (unsigned char *) malloc(width))) {
    fprintf(stderr, "memory allocation error\n");
  }

  else if (fits_get_num_rows(infptr,  &inrows,  &status) ||
           fits_get_num_rows(outfptr, &outrows, &status)) {
    fprintf(stderr, "Couldn't get the number of rows in the tables\n");
  }

  else  {
    /* check that the corresponding columns have the same datatypes */
    for (icol = 1; icol <= incols; icol++) {
      fits_get_coltype(infptr,  icol, &intype,  &inrep,  NULL, &status);
      fits_get_coltype(outfptr, icol, &outtype, &outrep, NULL, &status);
      if (intype != outtype || inrep != outrep) {
        fprintf(stderr, "Column %d is not the same in both tables\n", icol);
        check = 0;
      }
    }

    if (check && !status) {
      /* insert 'inrows' empty rows at the end of the output table */
      fits_insert_rows(outfptr, outrows, inrows, &status);

      for (ii = 1, jj = outrows +1; ii <= inrows; ii++, jj++) {
        /* read row from input and write it to the output table */
        fits_read_tblbytes(infptr,  ii, 1, width, buffer, &status);
        fits_write_tblbytes(outfptr, jj, 1, width, buffer, &status);
        if (status) {
          break;  /* jump out of loop if error occurred */
        }
      }

      /* all done; now free memory and close files */
      fits_close_file(outfptr, &status);
      fits_close_file(infptr,  &status);
    }
  }

  if (buffer) {
    free(buffer);
  }

  if (status) {
    fits_report_error(stderr, status);  /* print any error message */
  }
  return (status);
}
