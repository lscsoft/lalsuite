/* $Id$ */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "SFTReferenceLibrary.h"

int main(int argc, char **argv) {
  int i;

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

    /* and read successive SFTs blocks from the file and print headers */
    for (count=0; 1; count++) {

      struct headertag2 info,lastinfo;
      int err=0, swapendian, move;
      char *mycomment;

      err=ReadSFTHeader(fp, &info, &mycomment, &swapendian, 1);
      
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
      else {
	printf("File name:            %s\n", argv[i]);
	printf("SFT Version:          %.0f\n", info.version);
	printf("GPS_sec:              %d\n", info.gps_sec);
	printf("GPS_nsec:             %d\n", info.gps_nsec);
	printf("Timebase:             %-16f\n", info.tbase);
	printf("First frequency bin:  %d\n", info.firstfreqindex);
	printf("Number of freq bins:  %d\n", info.nsamples);
	printf("Detector prefix:      %c%c\n", info.detector[0], info.detector[1]);
	printf("64-bit CRC checksum:  %llu\n", info.crc64);
	printf("Comment length bytes: %d\n", info.comment_length);
	printf("Endian order:         %s\n", swapendian?"reversed":"native");
	if (info.comment_length) {
	  printf("Comment:              %s\n", mycomment);
	  free(mycomment);
	}
	printf("\n");
      }
      
      /* check that various bits of header information are consistent */
      if (count && (err=CheckSFTHeaderConsistency(&lastinfo, &info))) {
	fprintf(stderr, "%s is not a valid SFT. %s\n", argv[i], SFTErrorMessage(err));
	return err;
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
