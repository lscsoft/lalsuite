/*
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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
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
  
  fprintf(stdout, "%s: %s %s\n", argv[0], lalAppsVCSInfo.vcsId, lalAppsVCSInfo.vcsStatus);

  /* loop over all file names on command line */
  for (int i=1; i<argc; i++) {
    /* we on purpose do not call the XLAL version here,
     * so as to have the same stdout printing
     * and return code handling as older versions of this executable.
     */
    int errcode = ValidateSFTFile(argv[i]);
    if (errcode != 0) {
        return errcode;
    }
  }
  return 0;
}
