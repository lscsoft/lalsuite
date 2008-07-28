/*

  This program reads in binary SFTs (v1 and v2) and writes out narrow-banded merged SFTs (v2).

  It links to the SFTReferenceLibrary. To compile, use somehting like
  gcc -Wall -g -O2 splitSFTs.c -o splitSFTs libSFTReferenceLibrary.a -lm

  The frequency bands of the ouput SFTs (first frequency bin of first output SFT,
  last frequency bin of last output SFT, number of bins in each output SFT)
  can be specified either in bins ('-s', '-e', '-b') or Hz ('-fs', '-fe', '-fb')
  (or mixed - if both are given, frequency values take precedence).

  A "mystery factor" can be specified with '-m' option.

  In case of reading v1 SFTs (which don't support detector information in the header) the detector
  needs to be specified on the command-line using '-d' option.

  The name of the output SFTs is created by appending the start bin of the narrow-band SFT to the
  "output prefix" that can be given to the program with '-o'. If an output file already exists,
  the program will append the new SFTs to them, making it possible to construct the final
  narrow-band SFTs by running the program multiple times with different input SFTs. The GPS
  timestamps of the input SFTs need to be in ascending order to get valid merged SFT files.

  The last option on the command-line needs to be '-i', followed by as many input files as you wish
  (or the OS supports - using xargs should be simple with this command-line syntax).

  The program adds its own RCSID and command-line to the comment of the written SFTs,
  a mystery factor should show up as "xxx" there.

  Copyright (C) 2008 Bernd Machenschalk

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SFTReferenceLibrary.h"

#define RCSID "$Id: splitSFTs.c,v 1.27 2008/07/28 17:10:08 bema Exp $"

/* rounding (for positive numbers!)
   taken from SFTfileIO in LALSupport, should be consistent with that */
#define MYROUND(x) ( floor( (x) + 0.5 ) )

#define FALSE 0
#define TRUE (!FALSE)

#define CMT_NONE 0
#define CMT_OLD  1
#define CMT_FULL 2

/* error if value is nonzero */
#define TRY(v,c,errex) { int r; if((r=(v))) { fprintf(stderr,c " (%d @ %d)\n", r, __LINE__); exit(errex); } }
/* for SFT library calls write out the corresponding SFTErrorMessage, too */
#define TRYSFT(v,c) { int r; if((r=(v))) { fprintf(stderr,c " (%s @ %d)\n", SFTErrorMessage(r), __LINE__); exit(r); } }

int main(int argc, char**argv) {
  unsigned int arg;       /* current command-line argument */
  unsigned int bin;       /* current bin */
  struct headertag2 hd;   /* header of input SFT */
  FILE *fp;               /* currently open filepointer */
  char *oldcomment;       /* comment of input SFT */
  char *cmdline = NULL;   /* records command-line to add it to comment */
  char *comment = NULL;   /* comment to be written into output SFT file */
  int swap;               /* do we need to swap bytes? */
  float *data;            /* SFT data */
  char *outname;          /* name of output SFT file */
  char *prefix = "";      /* output filename prefix */
  char *detector = NULL;  /* detector name */
  double factor = 1.0;    /* "mystery" factor */
  int firstfile = TRUE;   /* are we processing the first input SFT file? */
  int add_comment = CMT_FULL; /* add RCSID and full command-line to every SFT file */
  unsigned int start = 0, end = 0, width = 0;     /* start, end and width in bins */
  double fMin = -1.0, fMax = -1.0, fWidth = -1.0; /* start, end and width in Hz */

  /* help / usage message */
  if((argv[1] == NULL) ||
     (strcmp(argv[1], "-h") == 0) || 
     (strcmp(argv[1], "--help") == 0)) {
    fprintf(stderr,
	    "%s -h\n"
	    "%s [-c 0|1|2] [-s <startbin>] [-e <endbin (exclusively)>] [-b <sftbins>]"
	    " [-fs <startfrequency>] [-fe <endfrequency (exclusively)>] [-fb <frequencywidth>]"
	    " [-m <factor>] [-d <detector>] [-o <outputprefix>] -i <inputfile> ...\n",
	    argv[0], argv[0]);
    exit(0);
  }

  /* record RCSID and command-line for the comment */
  TRY((cmdline = (char*)malloc(strlen(RCSID)+2)) == NULL,
      "out of memory allocating cmdline",1);
  strcpy(cmdline,RCSID);
  strcat(cmdline, "\n");
  for(arg = 0; arg < argc; arg++) {
    /* obscure the mystery factor */
    if (strcmp(argv[arg], "-m") == 0) {
      TRY((cmdline = (char*)realloc((void*)cmdline, strlen(cmdline) + 8)) == NULL,
	  "out of memory allocating cmdline",2);
      strcat(cmdline, "-m xxx ");
      arg++;
    } else {
      TRY((cmdline = (char*)realloc((void*)cmdline, strlen(cmdline) + strlen(argv[arg]) + 2)) == NULL,
	  "out of memory allocating cmdline",3);
      strcat(cmdline, argv[arg]);
      if(arg == argc - 1)
	strcat(cmdline, "\n");
      else
	strcat(cmdline, " ");
    }
  }

  /* get parameters from command-line */
  for(arg = 1; arg < argc; arg++) {
    if(strcmp(argv[arg], "-d") == 0) {
      detector = argv[++arg];
    } else if(strcmp(argv[arg], "-s") == 0) {
      start = atoi(argv[++arg]);
    } else if(strcmp(argv[arg], "-c") == 0) {
      add_comment = atoi(argv[++arg]);
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
    } else {
      fprintf(stderr, "unknown option '%s', try '-h' for help\n", argv[arg]);
      exit (-1);
    }
  }

  /* check if there was a "-i" at all */
  TRY(argv[arg] == NULL, "no input files specified",4);
  TRY(strcmp(argv[arg], "-i") != 0, "no input files specified",5);

  /* allocate space for output filename */
  TRY((outname = (char*)malloc(strlen(prefix) + 20)) == NULL,
      "out of memory allocating outname",6);

  /* loop over all input files
     first skip the "-i" option */
  for(arg++; arg < argc; arg++) {    

    /* open input SFT */
    TRY((fp = fopen(argv[arg], "r")) == NULL,
	"could not open SFT file for reading",7);
    
    /* read header */
    TRYSFT(ReadSFTHeader(fp, &hd, &oldcomment, &swap, 1),
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
    TRY((data = (float*)calloc(hd.nsamples, 2*sizeof(float))) == NULL,
	"out of memory allocating data",8);

    /* error if start < hd.firstfreqindex */
    if(start < hd.firstfreqindex) {
      fprintf(stderr,
	      "ERROR: start bin (%d) is smaller than first bin in input SFT (%d)\n",
	      start, hd.firstfreqindex);
      exit(9);
    }

    /* error if end > hd.firstfreqindex + hd.nsamples - 1 */
    if(start + width > end + 1) {
      fprintf(stderr,
	      "ERROR: end bin (%d) is larger than last bin in input SFT (%d)\n",
	      end, hd.firstfreqindex+hd.nsamples - 1);
      exit(10);
    }

    if (add_comment > CMT_OLD) {

      /* allocate space for new comment */
      TRY((comment = (char*)malloc(hd.comment_length + strlen(cmdline) + 1)) == NULL,
	  "out of memory allocating comment",11);
      
      /* append the commandline of this program to the old comment */
      if (oldcomment)
	strcpy(comment,oldcomment);
      else
	*comment = '\0';
      strcat(comment,cmdline);

    } else if (add_comment == CMT_OLD) {

      comment = oldcomment;

    }

    /* get the detector name from SFT header if present there (v2 SFTs),
       or else it needs to have been set on the command-line */
    if(hd.detector)
      detector = hd.detector;

    /* if no detector has been specified, issue an error */
    TRY(detector == NULL, "When reading v1 SFTs a detector needs to be specified with -d",12);

    /* read in complete SFT data */
    TRYSFT(ReadSFTData(fp, data, hd.firstfreqindex, hd.nsamples, NULL, NULL),
	   "could not read SFT data");

    /* apply factor */
    for(bin = 0; bin < 2*hd.nsamples; bin++)
      data[bin] *= factor;

    /* cleanup */
    fclose(fp);

    /* loop over start bins for output SFTs */
    for(bin = start; bin < end; bin += width) {

      int last_index = hd.firstfreqindex + hd.nsamples - 1;
      int last_bin = bin + width - 1;
      int this_width1 = last_bin <= last_index ? width : width - (last_bin - last_index);
      int this_width2 = end - bin + 1;
      int this_width = this_width1 < this_width2 ? this_width1 : this_width2;

      /* construct output SFT filename */
      sprintf(outname, "%s%d", prefix, bin);

      /* append the SFT to the "merged" SFT with the same name */
      TRY((fp = fopen(outname,"a")) == NULL,
	  "could not open SFT for writing",13);

      /* write the data */
      if (firstfile) {
	/* write the comment only to the first SFT of a "block" */
	TRYSFT(WriteSFT(fp, hd.gps_sec, hd.gps_nsec, hd.tbase, 
			bin, this_width, detector, comment,
			data + 2 * (bin - hd.firstfreqindex)),
	       "could not write SFT data");
      } else {
	/* not first SFT => NULL comment, everything else identical */
	TRYSFT(WriteSFT(fp, hd.gps_sec, hd.gps_nsec, hd.tbase, 
			bin, this_width, detector, NULL,
			data + 2 * (bin - hd.firstfreqindex)),
	       "could not write SFT data");
      }

      /* cleanup */
      fclose(fp);

    } /* loop over output SFTs */

    /* cleanup */
    if (add_comment > CMT_OLD)
      free(comment);
    free(data);

    /* next file is not the first file anymore */
    firstfile = FALSE;

  } /* loop over input SFTs */

  /* cleanup */
  free(outname);
    if (add_comment > CMT_OLD)
    free(cmdline);

  return(0);
}
