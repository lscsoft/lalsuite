/*
  - splits a SFTv2 into multiple ones containing narrow frequency bands.
  - uses the SFTReferenceLibrary, compile with
  gcc -Wall -O2 splitSFTs.c -o splitSFTs libSFTReferenceLibrary.a

  Author: Bernd Machenschalk
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SFTReferenceLibrary.h"

#define RCSID "$Id: splitSFTs.c,v 1.7 2008/07/15 12:40:16 bema Exp $"

/* rounding for positive numbers!
   taken from SFTfileIO in LALSupport, should be consistent with that */
#define MYROUND(x) ( floor( (x) + 0.5 ) )

/* error if value is nonzero */
#define TRY(v,c) { int r; if((r=(v))) { fprintf(stderr,c " (%d)\n", r); exit(-1); } }

int main(int argc, char**argv) {
  unsigned int arg;
  unsigned int start = 0, end = 0, width = 0;
  unsigned int bin;
  struct headertag2 hd;
  FILE *fp;
  char *oldcomment;
  char *cmdline = NULL;
  char *comment = NULL;
  int swap;
  float *data;
  char *outname;
  char *prefix = "";
  double factor = 1.0;
  double fMin = -1.0, fMax = -1.0, fWidth = -1.0;

  /* help / usage message */
  if((argv[1] == NULL) ||
     (strcmp(argv[1], "-h") == 0) || 
     (strcmp(argv[1], "--help") == 0)) {
    fprintf(stderr,
	    "%s -h\n"
	    "%s [-s <startbin>] [-e <endbin>] [-b <sftbins>]"
	    " [-fs <startfrequency>] [-fs <endfrequency>] [-fs <frequencywidth>]"
	    " [-m <factor>] [-o <outputprefix>] -i <inputfile> ...\n", argv[0], argv[0]);
    exit(0);
  }

  /* record the commandline for the comment */
  cmdline = malloc(strlen(RCSID)+2);
  strcpy(cmdline,RCSID);
  strcat(cmdline, "\n");
  for(arg = 0; arg < argc; arg++) {
    if (strcmp(argv[arg], "-m") == 0) {
      cmdline = realloc((void*)cmdline, strlen(cmdline) + 8);
      strcat(cmdline, "-m xxx ");
      arg++;
    } else {
      cmdline = realloc((void*)cmdline, strlen(cmdline) + strlen(argv[arg]) + 2);
      strcat(cmdline, argv[arg]);
      if(arg == argc - 1)
	strcat(cmdline, "\n");
      else
	strcat(cmdline, " ");
    }
  }

  /* get parameters from command-line */
  for(arg = 1; arg < argc; arg++) {
    if(strcmp(argv[arg], "-s") == 0) {
      start = atoi(argv[++arg]);
    } else if(strcmp(argv[arg], "-e") == 0) {
      end = atoi(argv[++arg]);
    } else if(strcmp(argv[arg], "-b") == 0) {
      width = atoi(argv[++arg]);
    } else if(strcmp(argv[arg], "-fs") == 0) {
      fMin = atof(argv[++arg]);
    } else if(strcmp(argv[arg], "-fe") == 0) {
      fMax = atof(argv[++arg]);
    } else if(strcmp(argv[arg], "-fb") == 0) {
      fWidth = atof(argv[++arg]);
    } else if(strcmp(argv[arg], "-m") == 0) {
      factor = atof(argv[++arg]);
    } else if(strcmp(argv[arg], "-o") == 0) {
      prefix = argv[++arg];
    } else if(strcmp(argv[arg], "-i") == 0) {
      break;
    }
  }

  /* allocate space for output filename */
  TRY((outname = (char*)malloc(strlen(prefix) + 20)) == NULL,
      "out of memory allocating outname");

  /* loop over all input files
     first skip the "-i" option */
  for(arg++; arg < argc; arg++) {    

    /* open input SFT */
    TRY((fp = fopen(argv[arg], "r")) == NULL,
	"could not open SFT file for read");
    
    /* read header */
    TRY(ReadSFTHeader(fp, &hd, &oldcomment, &swap, 1),
	"could not read SFT header");

    /* calculate bins from frequency parameters if they were given */
    /* deltaF = 1.0 / tbase; bins = freq / deltaF => bins = freq * tbase */
    if(fMin >= 0.0)
      start = MYROUND(fMin * hd.tbase);
    if(fMax >= 0.0)
      end   = MYROUND(fMax * hd.tbase);
    if(fWidth >= 0.0)
      width = MYROUND(fWidth * hd.tbase);
    
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
    for(bin = start; bin < end - width + 1 ; bin += width) {

      /* construct output SFT filename */
      sprintf(outname, "%s%d", prefix, bin);

      /* append the SFT to the "merged" SFT with the same name */
      TRY((fp = fopen(outname,"a")) == NULL,
	  "could not open SFT for writing");

      /* write the data */
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
