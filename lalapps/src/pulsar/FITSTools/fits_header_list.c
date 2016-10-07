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
  fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
  char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
  int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
  int single = 0, hdupos = 0, nkeys = 0, ii = 0;

  if (argc != 2) {
    fprintf(stderr, "Usage:  %s filename[ext] \n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "List the FITS header keywords in a single extension, or, if \n");
    fprintf(stderr, "ext is not given, list the keywords in all the extensions. \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "   %s file.fits      - list every header in the file \n", argv[0]);
    fprintf(stderr, "   %s file.fits[0]   - list primary array header \n", argv[0]);
    fprintf(stderr, "   %s file.fits[2]   - list header of 2nd extension \n", argv[0]);
    fprintf(stderr, "   %s file.fits+2    - same as above \n", argv[0]);
    fprintf(stderr, "   %s file.fits[GTI] - list header of GTI extension\n", argv[0]);
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

    /* List only a single header if a specific extension was given */
    if (hdupos != 1 || strchr(argv[1], '[')) {
      single = 1;
    }

    for (; !status; hdupos++) { /* Main loop through each extension */
      fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of keywords */

      fprintf(fout, "Header listing for HDU #%d:\n", hdupos);

      for (ii = 1; ii <= nkeys; ii++) { /* Read and print each keywords */

        if (fits_read_record(fptr, ii, card, &status)) {
          break;
        }
        fprintf(fout, "%s\n", card);
      }
      fprintf(fout, "END\n\n");  /* terminate listing with END */

      if (single) {
        break;  /* quit if only listing a single header */
      }

      fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
    }

    if (status == END_OF_FILE) {
      status = 0;  /* Reset after normal error */
    }

    fits_close_file(fptr, &status);
  }

  pclose(fout);

  if (status) {
    fits_report_error(stderr, status);  /* print any error message */
  }
  return (status);
}
