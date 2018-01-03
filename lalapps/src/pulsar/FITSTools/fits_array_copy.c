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
#include <stdlib.h>

#if defined(HAVE_LIBCFITSIO)
#include <fitsio.h>
#else
#error CFITSIO library is not available
#endif

int main(int argc, char *argv[])
{
  fitsfile *infptr = 0, *outfptr = 0;   /* FITS file pointers defined in fitsio.h */
  int status = 0, ii = 1, iteration = 0, single = 0, hdupos = 0;
  int hdutype = 0, bitpix = 0, bytepix = 0, naxis = 0, nkeys = 0, datatype = 0, anynul = 0;
  long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  long first = 0, totpix = 0, npix = 0;
  double *array = 0, bscale = 1.0, bzero = 0.0, nulval = 0.;
  char card[81];

  if (argc != 3) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:  %s inputArray outputArray[compress]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Copy an input array to an output array, optionally compressing\n");
    fprintf(stderr, "or uncompressing the array in the process.  If the [compress]\n");
    fprintf(stderr, "qualifier is appended to the output file name then the input array\n");
    fprintf(stderr, "will be compressed using the tile-compressed format.  In this format,\n");
    fprintf(stderr, "the array is divided into rectangular tiles and each tile of pixels\n");
    fprintf(stderr, "is compressed and stored in a variable-length row of a binary table.\n");
    fprintf(stderr, "If the [compress] qualifier is omitted, and the input array is\n");
    fprintf(stderr, "in tile-compressed format, then the output array will be uncompressed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "If an extension name or number is appended to the input file name, \n");
    fprintf(stderr, "enclosed in square brackets, then only that single extension will be\n");
    fprintf(stderr, "copied to the output file.  Otherwise, every extension in the input file\n");
    fprintf(stderr, "will be processed in turn and copied to the output file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "1)  %s array.fit 'carray.fit[compress]'\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    This compresses the input array using the default parameters, i.e.,\n");
    fprintf(stderr, "    using the Rice compression algorithm and using row by row tiles.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "2)  %s carray.fit array2.fit\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    This uncompress the array created in the first example.\n");
    fprintf(stderr, "    array2.fit should be identical to array.fit if the array\n");
    fprintf(stderr, "    has an integer datatype.  There will be small differences\n");
    fprintf(stderr, "    in the pixel values if it is a floating point array.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "3)  %s array.fit 'carray.fit[compress GZIP 100,100;4]'\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    This compresses the input array using the following parameters:\n");
    fprintf(stderr, "         GZIP compression algorithm;\n");
    fprintf(stderr, "         100 X 100 pixel compression tiles;\n");
    fprintf(stderr, "         noise_bits = 4 (only used with floating point arrays)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The full syntax of the compression qualifier is:\n");
    fprintf(stderr, "    [compress ALGORITHM TDIM1,TDIM2,...; NOISE_BITS]\n");
    fprintf(stderr, "where the allowed ALGORITHM values are Rice, GZIP, PLIO, \n");
    fprintf(stderr, "and TDIMn is the size of the compression tile in each dimension,\n");
    fprintf(stderr, "and NOISE_BITS = 1, 2, 3, or 4 and controls the amount of noise\n");
    fprintf(stderr, "suppression when compressing floating point arrays. \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note that it may be necessary to enclose the file names\n");
    fprintf(stderr, "in single quote characters on the Unix command line.\n");
    return (0);
  }

  /* Open the input file and create output file */
  fits_open_file(&infptr, argv[1], READONLY, &status);
  fits_create_file(&outfptr, argv[2], &status);

  if (status != 0) {
    fits_report_error(stderr, status);
    return (status);
  }

  fits_get_hdu_num(infptr, &hdupos);  /* Get the current HDU position */

  /* Copy only a single HDU if a specific extension was given */
  if (hdupos != 1 || strchr(argv[1], '[')) {
    single = 1;
  }

  for (; !status; hdupos++) { /* Main loop through each extension */

    fits_get_hdu_type(infptr, &hdutype, &status);

    if (hdutype == IMAGE_HDU) {

      /* get array dimensions and total number of pixels in array */
      for (ii = 0; ii < 9; ii++) {
        naxes[ii] = 1;
      }

      fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);

      totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4]
        * naxes[5] * naxes[6] * naxes[7] * naxes[8];
    }

    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0) {

      /* just copy tables and null arrays */
      fits_copy_hdu(infptr, outfptr, 0, &status);

    } else {

      /* Explicitly create new array, to support compression */
      fits_create_img(outfptr, bitpix, naxis, naxes, &status);

      /* copy all the user keywords (not the structural keywords) */
      fits_get_hdrspace(infptr, &nkeys, NULL, &status);

      for (ii = 1; ii <= nkeys; ii++) {
        fits_read_record(infptr, ii, card, &status);
        if (fits_get_keyclass(card) > TYP_CMPRS_KEY) {
          fits_write_record(outfptr, card, &status);
        }
      }

      switch (bitpix) {
      case BYTE_IMG:
        datatype = TBYTE;
        break;
      case SHORT_IMG:
        datatype = TSHORT;
        break;
      case LONG_IMG:
        datatype = TLONG;
        break;
      case FLOAT_IMG:
        datatype = TFLOAT;
        break;
      case DOUBLE_IMG:
        datatype = TDOUBLE;
        break;
      }

      bytepix = abs(bitpix) / 8;

      npix = totpix;
      iteration = 0;

      /* try to allocate memory for the entire array */
      /* use double type to force memory alignment */
      array = (double *) calloc(npix, bytepix);

      /* if allocation failed, divide size by 2 and try again */
      while (!array && iteration < 10)  {
        iteration++;
        npix = npix / 2;
        array = (double *) calloc(npix, bytepix);
      }

      if (!array)  {
        fprintf(stderr, "Memory allocation error\n");
        return (0);
      }

      /* turn off any scaling so that we copy the raw pixel values */
      fits_set_bscale(infptr,  bscale, bzero, &status);
      fits_set_bscale(outfptr, bscale, bzero, &status);

      first = 1;
      while (totpix > 0 && !status) {
        /* read all or part of array then write it back to the output file */
        fits_read_img(infptr, datatype, first, npix,
                      &nulval, array, &anynul, &status);

        fits_write_img(outfptr, datatype, first, npix, array, &status);
        totpix = totpix - npix;
        first  = first  + npix;
      }
      free(array);
    }

    if (single) {
      break;  /* quit if only copying a single HDU */
    }
    fits_movrel_hdu(infptr, 1, NULL, &status);  /* try to move to next HDU */
  }

  if (status == END_OF_FILE) {
    status = 0;  /* Reset after normal error */
  }

  fits_close_file(outfptr,  &status);
  fits_close_file(infptr, &status);

  /* if error occurred, print out error message */
  if (status) {
    fits_report_error(stderr, status);
  }
  return (status);
}
