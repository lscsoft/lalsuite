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

#if defined(HAVE_LIBCFITSIO)

#include <fitsio.h>

// If fffree() is missing, use free() instead
#if !defined(HAVE_FFFREE)
#undef fits_free_memory
#define fits_free_memory(ptr, status) free(ptr)

// If ffree() is present but not declared, declare it
#elif defined(HAVE_DECL_FFFREE) && !HAVE_DECL_FFFREE
int fffree( void *, int * );
#undef fits_free_memory
#define fits_free_memory fffree

#endif // ffree()

#endif // defined(HAVE_LIBCFITSIO)

int main(int argc, char *argv[])
{
  fitsfile *fptr = 0;         /* FITS file pointer, defined in fitsio.h */
  int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */

  int printhelp = (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0));

  if (printhelp || argc != 3) {
    fprintf(stderr, "Usage:  %s filename[ext] keyword\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Delete a header keyword.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  %s file.fits dec      - delete the DEC keyword\n", argv[0]);
    return (0);
  }

  if (!fits_open_file(&fptr, argv[1], READWRITE, &status)) {
    fits_delete_key(fptr,argv[2], &status);
    fits_close_file(fptr, &status);
  }    /* open_file */

  /* if error occured, print out error message */
  if (status) {
    fits_report_error(stderr, status);
  }
  return (status);
}
