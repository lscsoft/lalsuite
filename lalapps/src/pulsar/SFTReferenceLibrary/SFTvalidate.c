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
 * \ingroup pulsarApps
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
#include <LALAppsVCSInfo.h>
#include "SFTReferenceLibrary.h"

int main(int argc, char** argv) {
  int i;
  float *data=NULL;
  
  fprintf(stdout, "%s: %s %s\n", argv[0], lalAppsVCSId, lalAppsVCSStatus);

  /* loop over all file names on command line */
  for (i=1; i<argc; i++) {
    FILE *fp;
    int count;

    /* open the file */
    if (!(fp=fopen(argv[i], "r"))) {
      fprintf(stderr,"Unable to open %s", argv[i]);
      if (errno)
	perror(" ");
      return SFTENULLFP;
    }

    /* and read successive SFTs blocks from the file and validate CRC
       checksums */
    for (count=0; 1; count++) {
      struct headertag2 info,lastinfo;
      int err=0, swapendian, move, j;
      
      err=ReadSFTHeader(fp, &info, NULL, &swapendian, 1);

      /* at end of SFT file or merged SFT file blocks */
      if (err==SFTENONE && count)
	break;
      
      /* SFT was invalid: say why */
      if (err) {
	fprintf(stderr, "%s is not a valid SFT. %s\n", argv[i], SFTErrorMessage(err));
	if (errno)
	  perror(NULL);
	return err;
      }

      /* check that various bits of header information are consistent */
      if (count && (err=CheckSFTHeaderConsistency(&lastinfo, &info)))
	{
	  fprintf(stderr, "%s is not a valid SFT. %s\n", argv[i], SFTErrorMessage(err));
	  if (errno)
	    perror(NULL);
	  return err;
	}
      
      /* check that data appears valid */
      data=(float *)realloc((void *)data, info.nsamples*4*2);
      if (!data) {
	errno=SFTENULLPOINTER;
	fprintf(stderr, "ran out of memory at %s. %s\n", argv[i], SFTErrorMessage(err));
	if (errno)
	  perror(NULL);
	return err;
      }

      err=ReadSFTData(fp, data, info.firstfreqindex, info.nsamples, /*comment*/ NULL, /*headerinfo */ NULL);
      if (err) {
	fprintf(stderr, "%s is not a valid SFT. %s\n", argv[i], SFTErrorMessage(err));
	if (errno)
	  perror(NULL);
	return err;
      }

      for (j=0; j<info.nsamples; j++) {
	if (!isfinite(data[2*j]) || !isfinite(data[2*j+1])) {
	  fprintf(stderr, "%s is not a valid SFT (data infinite at freq bin %d)\n", argv[i], j+info.firstfreqindex);
	  return SFTNOTFINITE;
	}
      }

      /* keep copy of header for comparison the next time */
      lastinfo=info;
      
      /* Move forward to next SFT in merged file */
      if (info.version==1)
	move=sizeof(struct headertag1)+info.nsamples*2*sizeof(float);
      else
	move=sizeof(struct headertag2)+info.nsamples*2*sizeof(float)+info.comment_length;
      fseek(fp, move, SEEK_CUR);
    }
    fclose(fp);
  }
  return 0;
}
