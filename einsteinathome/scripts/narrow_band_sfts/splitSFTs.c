/*
  - splits a SFTv2 into multiple ones containing narrow frequency bands.
  - uses the SFTReferenceLibrary, compile with
  gcc -Wall -O2 splitSFTs.c -o splitSFTs libSFTReferenceLibrary.a

  Author: Bernd Machenschalk
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SFTReferenceLibrary.h"

/* error if value is nonzero */
#define TRY(v,c) { int r; if((r=(v))) { fprintf(stderr,c " (%d)\n", r); exit(-1); } }

int main(int argc, char**argv) {
  unsigned int arg;
  unsigned int start, end, width;
  unsigned int bin;
  struct headertag2 hd;
  FILE *fp;
  char *oldcomment;
  char *cmdline = NULL;
  char *comment = NULL;
  int swap;
  float *data;
  char *outname;
  char *prefix;
  double factor = 1.0;

  /* usage error */
  if(argc < 6) {
    fprintf(stderr, "%s <startbin> <endbin> <sftbins> <outputprefix> <inputfile> ...\n", argv[0]);
    exit(-1);
  } 

  /* get parameters from command line */
  start  = atoi(argv[1]);
  end    = atoi(argv[2]);
  width  = atoi(argv[3]);
  prefix = argv[4];

  /* allocate space for output filename */
  TRY((outname = (char*)malloc(strlen(prefix) + 20)) == NULL,
      "out of memory allocating outname");

  /* record the commandline for the comment */
  for(arg = 0; arg < argc; arg++) {
    cmdline = realloc((void*)cmdline, strlen(cmdline) + strlen(argv[arg]) + 2);
    strcat(cmdline, argv[arg]);
    if(arg == argc - 1)
      strcat(cmdline, "\n");
    else
      strcat(cmdline, " ");
  }

  /* loop over all input files */
  for(arg = 5; arg < argc; arg++) {    

    /* open input SFT */
    TRY((fp = fopen(argv[arg], "r")) == NULL,
	"could not open SFT file for read");
    
    /* read header */
    TRY(ReadSFTHeader(fp, &hd, &oldcomment, &swap, 1),
	"could not read SFT header");
    
    /* allocate space for SFT data */
    TRY((data = (float*)malloc(hd.nsamples * sizeof(float))) == NULL,
	"out of memory allocating data");

    /* allocate space for new comment */
    TRY((comment = (char*)malloc(hd.comment_length + strlen(cmdline) + 1)) == NULL,
	"out of memory allocating comment");

    /* append the commandline of this program to the old comment */
    strcpy(comment,oldcomment);
    strcat(comment,cmdline);

    /* read in complete SFT data */
    TRY(ReadSFTData(fp, data, 0, hd.nsamples, NULL, NULL),
	"could not read SFT data");

    /* apply factor */
    for(bin = 0; bin < hd.nsamples; bin++)
      data[bin] *= factor;

    /* cleanup */
    fclose(fp);

    /* loop over start bins for output SFTs */
    for(bin = start; bin < end; bin += width) {

      /* construct output SFT filename */
      sprintf(outname, "%s%d", prefix, bin);

      /* open for appending */
      TRY((fp = fopen(outname,"a")) == NULL,
	  "could not open SFT for writing");

      /* append the SFT to the "merged" SFT */
      TRY(WriteSFT(fp, hd.gps_sec, hd.gps_nsec, hd.tbase, 
		   bin, width, hd.detector, comment, data),
	  "could not write SFT data");

      /* cleanup */
      fclose(fp);
    }

    /* cleanup */
    free(comment);
    free(data);
  }

  /* cleanup */
  free(outname);
  free(cmdline);

  return(0);
}
