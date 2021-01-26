/*
 *  Copyright (C) 2021 Karl Wette
 *  Copyright (C) 2004, 2005 Bruce Allen
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/**
 * \author Bruce Allen
 * \file
 * \ingroup lalapps_pulsar_SFTTools
 * \brief
 * Verify that a set of SFT files is valid
 *
 * The exit status will be zero if all SFTs are valid.  The exit status
 * will be non-zero if any of the SFTs was invalid.  grep SFTE
 * SFTReferenceLibrary.h will show the return values.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <lal/LALStdlib.h>
#include <lal/SFTReferenceLibrary.h>
#include <LALAppsVCSInfo.h>

int main(int argc, char** argv) {
  
  if (argc == 2 && (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0)) {
    fprintf(stdout, "%s: %s %s\n", argv[0], lalAppsVCSInfo.vcsId, lalAppsVCSInfo.vcsStatus);
    return EXIT_SUCCESS;
  }

  if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
    fprintf(stdout, "usage:\n");
    fprintf(stdout, "   %s *.sft\n", argv[0]);
    fprintf(stdout, "   ls *.sft | %s\n", argv[0]);
    fprintf(stdout, "   find -name '*.sft' | %s >valid-sfts.txt 2>errors.log\n", argv[0]);
    return EXIT_SUCCESS;
  }

  int errcode = EXIT_SUCCESS;

  if (argc > 1) {

    /* loop over all file names on command line */
    for (int i = 1; i < argc; ++i) {
      /* we on purpose do not call the XLAL version here,
       * so as to have the same stdout printing
       * and return code handling as older versions of this executable.
       */
      if (ValidateSFTFile(argv[i]) != 0) {
        errcode = EXIT_FAILURE;
      } else {
        fprintf(stdout, "%s\n", argv[i]);
        fflush(stdout);
      }
    }

  } else {

    char line[2048];

    /* loop over all file names from standard input */
    while (fgets(line, sizeof(line) - 1, stdin) != NULL) {
      size_t len = strlen(line);
      if (len > 1) {
        line[len - 1] = '\0';
        /* we on purpose do not call the XLAL version here,
         * so as to have the same stdout printing
         * and return code handling as older versions of this executable.
         */
        if (ValidateSFTFile(line) != 0) {
          errcode = EXIT_FAILURE;
        } else {
          fprintf(stdout, "%s\n", line);
          fflush(stdout);
        }
      }
    }

  }

  return errcode;

}
