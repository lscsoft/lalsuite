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
  fitsfile *fptr = 0;         /* FITS file pointer, defined in fitsio.h */
  char card[FLEN_CARD], newcard[FLEN_CARD];
  char oldvalue[FLEN_VALUE], comment[FLEN_COMMENT];
  int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
  int keytype = 0;

  if (argc != 4) {
    fprintf(stderr, "Usage:  %s filename[ext] keyword newvalue\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Modify the value of a header keyword.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  %s file.fits dec 30.0 - set DEC = 30.0 \n", argv[0]);
    return (0);
  }

  if (!fits_open_file(&fptr, argv[1], READWRITE, &status)) {
    if (fits_read_card(fptr,argv[2], card, &status)) {
      fprintf(stderr, "Keyword does not exist\n");
      return (1);
    }

    if (fits_get_keyclass(card) == TYP_STRUC_KEY) {
      fprintf(stderr, "Protected keyword cannot be modified\n");
      return (1);
    }

    fits_parse_value(card, oldvalue, comment, &status);

    /* construct template for new keyword */
    strcpy(newcard, argv[2]);     /* copy keyword name */
    strcat(newcard, " = ");       /* '=' value delimiter */
    strcat(newcard, argv[3]);     /* new value */
    if (*comment) {
      strcat(newcard, " / ");  /* comment delimiter */
      strcat(newcard, comment);     /* append the comment */
    }

    /* reformat the keyword string to conform to FITS rules */
    fits_parse_template(newcard, card, &keytype, &status);

    /* overwrite the keyword with the new value */
    fits_update_card(fptr, argv[2], card, &status);

    fits_close_file(fptr, &status);
  }    /* open_file */

  /* if error occured, print out error message */
  if (status) {
    fits_report_error(stderr, status);
  }
  return (status);
}
