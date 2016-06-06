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
#include <stdio.h>

#if defined(HAVE_LIBCFITSIO)
#include <fitsio.h>
#else
#error CFITSIO library is not available
#endif

int main(int argc, char *argv[])
{
  fitsfile *infptr = 0, *outfptr = 0;   /* FITS file pointers defined in fitsio.h */
  int status = 0, ii = 1;       /* status must always be initialized = 0  */

  if (argc != 3) {
    fprintf(stderr, "Usage:  %s inputfile outputfile\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Copy an input file to an output file, optionally filtering\n");
    fprintf(stderr, "the file in the process.  This seemingly simple program can\n");
    fprintf(stderr, "apply powerful filters which transform the input file as\n");
    fprintf(stderr, "it is being copied.  Filters may be used to extract a\n");
    fprintf(stderr, "subimage from a larger image, select rows from a table,\n");
    fprintf(stderr, "filter a table with a GTI time extension or a SAO region file,\n");
    fprintf(stderr, "create or delete columns in a table, create an image by\n");
    fprintf(stderr, "binning (histogramming) 2 table columns, and convert IRAF\n");
    fprintf(stderr, "format *.imh or raw binary data files into FITS images.\n");
    fprintf(stderr, "See the CFITSIO User's Guide for a complete description of\n");
    fprintf(stderr, "the Extended File Name filtering syntax.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s in.fit out.fit                   (simple file copy)\n", argv[0]);
    fprintf(stderr, "%s - -                              (stdin to stdout)\n", argv[0]);
    fprintf(stderr, "%s in.fit[11:50,21:60] out.fit      (copy a subimage)\n", argv[0]);
    fprintf(stderr, "%s iniraf.imh out.fit               (IRAF image to FITS)\n", argv[0]);
    fprintf(stderr, "%s in.dat[i512,512] out.fit         (raw array to FITS)\n", argv[0]);
    fprintf(stderr, "%s in.fit[events][pi>35] out.fit    (copy rows with pi>35)\n", argv[0]);
    fprintf(stderr, "%s in.fit[events][bin X,Y] out.fit  (bin an image) \n", argv[0]);
    fprintf(stderr, "%s in.fit[events][col x=.9*y] out.fit        (new x column)\n", argv[0]);
    fprintf(stderr, "%s in.fit[events][gtifilter()] out.fit       (time filter)\n", argv[0]);
    fprintf(stderr, "%s in.fit[2][regfilter(\"pow.reg\")] out.fit (spatial filter)\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Note that it may be necessary to enclose the input file name\n");
    fprintf(stderr, "in single quote characters on the Unix command line.\n");
    return (0);
  }

  /* Open the input file */
  if (!fits_open_file(&infptr, argv[1], READONLY, &status)) {
    /* Create the output file */
    if (!fits_create_file(&outfptr, argv[2], &status)) {
      /* Copy every HDU until we get an error */
      while (!fits_movabs_hdu(infptr, ii++, NULL, &status)) {
        fits_copy_hdu(infptr, outfptr, 0, &status);
      }

      /* Reset status after normal error */
      if (status == END_OF_FILE) {
        status = 0;
      }

      fits_close_file(outfptr,  &status);
    }
    fits_close_file(infptr, &status);
  }

  /* if error occured, print out error message */
  if (status) {
    fits_report_error(stderr, status);
  }
  return (status);
}
