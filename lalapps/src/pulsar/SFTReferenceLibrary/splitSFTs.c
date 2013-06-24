/*
 *  Copyright (C) 2008, 2010 Bernd Machenschalk, Bruce Allen
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
 * \file
 * \ingroup pulsarApps
 * \author Bernd Machenschalk, Bruce Allen
 *
 * \brief This program reads in binary SFTs (v1 and v2) and writes out narrow-banded merged SFTs (v2).

  This code links to the SFTReferenceLibrary. To compile, use somehting like
  <code>gcc -Wall -g -O2 splitSFTs.c -o splitSFTs libSFTReferenceLibrary.a -lm</code>

  Writen by Bernd Machenschalk for Einstein\@home 2008

  * Revision splitSFTs.c,v 1.41 2008/10/29 16:54:13 was
    reviewed by the LSC CW Review Committee Thur Oct 30, 2008

    Suggested improvements:
    * issue warning at ambiguous frequency values that are so close to the boundary of a bin
      that rounding might end up giving an unintended bin
    * check for consistency of input SFTs (same timebase, ascending timestamps etc., see spec)
      and merged SFTs (check last header of a file we are appending to)

  * Other possible improvements not suggested by the committee
    * keep output files open (if there aren't too many)
    * obscure a mystery factor in command-line record even if given with long option --factor
*/

#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <LALAppsVCSInfo.h>
#include "SFTReferenceLibrary.h"

#define VCSID LALAPPS_VCS_IDENT_ID LALAPPS_VCS_IDENT_STATUS

/** rounding (for positive numbers!)
    taken from SFTfileIO in LALSupport, should be consistent with that */
#define MYROUND(x) ( floor( (x) + 0.5 ) )

#define FALSE 0
#define TRUE (!FALSE)

#define CMT_NONE 0
#define CMT_OLD  1
#define CMT_FULL 2

/* macros for error handling */
/** error if value is nonzero */
#define TRY(v,c,errex) { int r; if((r=(v))) { fprintf(stderr,c " (%d @ %d)\n", r, __LINE__); exit(errex); } }
/** for SFT library calls write out the corresponding SFTErrorMessage, too */
#define TRYSFT(v,c) { int r; if((r=(v))) { fprintf(stderr,c " (%s @ %d)\n", SFTErrorMessage(r), __LINE__); exit(r); } }

typedef struct {
  int resource_units;  /**< number of resource units available for consumption */
  int resource_rate;   /**< this many resource units become available in one second */
  time_t last_checked; /**< time we last checked */
} UNIT_SOURCE;

/** throtteling settings */
UNIT_SOURCE read_bandwidth={0, 0, 0};
UNIT_SOURCE read_open_rate={0, 0, 0};
UNIT_SOURCE write_bandwidth={0, 0, 0};
UNIT_SOURCE write_open_rate={0, 0, 0};

void request_resource(UNIT_SOURCE *, int);

/** request a resurce. Function returns after waiting for throttle time */
void request_resource(UNIT_SOURCE *us, int units) {
  time_t now;
  int seconds;

  /* negative rate indicates no throttling */
  if(us->resource_rate <= 0)
    return;

  us->resource_units -= units;
  while(us->resource_units <= 0) {
    time(&now);
    seconds = now-us->last_checked;
    /* guard against overflow and condor checkpointing */
    if(seconds < 0)
      seconds = 0;
    if(seconds > 10)
      seconds = 10;
    if(seconds == 0)
      sleep(1);
    us->resource_units += seconds*us->resource_rate;
    us->last_checked = now;
  }
}

/** main program */
int main(int argc, char**argv) {
  int arg;                        /* current command-line argument */
  unsigned int bin;               /* current bin */
  struct headertag2 hd;           /* header of input SFT */
  FILE *fp;                       /* currently open filepointer */
  char *oldcomment;               /* comment of input SFT */
  char *cmdline = NULL;           /* records command-line to add it to comment */
  char *comment = NULL;           /* comment to be written into output SFT file */
  int swap;                       /* do we need to swap bytes? */
  float *data;                    /* SFT data */
  char *outname;                  /* name of output SFT file */
  char empty = '\0';
  char *prefix = &empty;          /* output filename prefix */
  char *detector = NULL;          /* detector name */
  double factor = 1.0;            /* "mystery" factor */
  double conversion_factor = 1.0; /* extra factor needed when converting from v1 SFTs */
  int firstfile = TRUE;           /* are we processing the first input SFT file? */
  int allcomments = FALSE;        /* write comment into _every_ SFT in the file */
  int add_comment = CMT_FULL;     /* add VCS ID and full command-line to every SFT file */
  unsigned int start = 0, end = 0;     /* start and end in bins */
  unsigned int width = 0, overlap = 0; /* width and overlap in bins */
  double fMin = -1.0, fMax = -1.0;     /* start and end in Hz */
  double fWidth = -1.0, fOverlap = -1; /* width and overlap in Hz */
  unsigned int nactivesamples;         /* number of bins to actually read in */

  /* initialize throtteling */
  time(&read_bandwidth.last_checked);
  time(&read_open_rate.last_checked);
  time(&write_bandwidth.last_checked);
  time(&write_open_rate.last_checked);

  /* help / usage message */
  if((argv[1] == NULL) ||
     (strcmp(argv[1], "-h") == 0) || 
     (strcmp(argv[1], "--help") == 0)) {
    fprintf(stderr,
	    "%s -h\n"
	    "\n"
	    "  Write this help message\n"
	    "\n"
	    "%s [-c 0|1|2] [-a] [-s <startbin>] [-e <endbin (exclusively)>] [-b <sftbins>]\n"
	    "  [-fs <startfrequency>] [-fe <endfrequency (exclusively)>] [-fb <frequencywidth>]\n"
	    "  [-x <overlap>] [-fx overlap] [-m <factor>] [-d <detector>] [-o <outputprefix>]\n"
	    "  -i <inputfile> ...\n"
	    "\n"
	    "  This program reads in binary SFTs (v1 and v2) and writes out narrow-banded\n"
	    "  merged SFTs (v2).\n"
	    "\n"
	    "  The frequency bands of the ouput SFTs (first frequency bin of first output SFT,\n"
	    "  last frequency bin of last output SFT, number of bins in each output SFT)\n"
	    "  and a possible overlap of the output files can be specified\n"
	    "  either in bins ('-s', '-e', '-b', -'x') or Hz ('-fs', '-fe', '-fb', '-fx')\n"
	    "  (or mixed - if both are given, frequency values take precedence).\n"
	    "\n"
	    "  A 'mystery factor' can be specified with '-m' option.\n"
	    "\n"
	    "  In case of reading v1 SFTs (which don't support detector information in the header)\n"
	    "  the detector needs to be specified on the command-line using '-d' option.\n"
	    "\n"
	    "  The name of the output SFTs is created by appending the start bin of the narrow-band\n"
	    "  SFT to the 'output prefix' that can be given to the program with '-o'. If an output\n"
	    "  file already exists, the program will append the new SFTs to them, making it possible\n"
	    "  to construct the final narrow-band SFTs by running the program multiple times with\n"
	    "  different input SFTs. The GPS timestamps of the input SFTs need to be in ascending\n"
	    "  order to get valid merged SFT files.\n"
	    "\n"
	    "  The '-c' options specifies how to deal with comments - 0 means no comment is written\n"
	    "  at all, 1 means that the comment is taken unmodified from the input SFTs, 2 (default)\n"
	    "  means that the program appends its RCS id and command-line to the comment.\n"
	    "  By default a comment is written only to the first SFT of a merged SFT output 'block'\n"
	    "  (i.e. call to this program). Adding the option '-a' to the command line specifies that\n"
	    "  instead the comment is written into every SFT in the resulting file.\n"
	    "\n"
	    "  The last option on the command-line needs to be '-i', followed by as many input files\n"
	    "  as you wish (or the OS supports - using xargs should be simple with this command-line\n"
	    "  syntax).\n"
	    "\n"
	    "  The program adds its own VCS ID and command-line to the comment of the written SFTs,\n"
	    "  a mystery factor should show up as 'xxx' there.\n",
	    argv[0], argv[0]);
    exit(0);
  }

  /* record VCS ID and command-line for the comment */
  TRY((cmdline = (char*)malloc(strlen(VCSID)+2)) == NULL,
      "out of memory allocating cmdline",1);
  strcpy(cmdline,VCSID);
  strcat(cmdline, "\n");
  for(arg = 0; arg < argc; arg++) {
    if (strcmp(argv[arg], "-m") == 0) {
      /* obscure the mystery factor */
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
    } else if((strcmp(argv[arg], "-c") == 0) ||
	      (strcmp(argv[arg], "--add-comment") == 0)) {
      add_comment = atoi(argv[++arg]);
    } else if((strcmp(argv[arg], "-a") == 0) ||
	      (strcmp(argv[arg], "--all-comments") == 0)) {
      allcomments = TRUE;
    } else if((strcmp(argv[arg], "-s") == 0) ||
	      (strcmp(argv[arg], "--start-bin") == 0)) {
      start = atoi(argv[++arg]);
    } else if((strcmp(argv[arg], "-e") == 0) ||
	      (strcmp(argv[arg], "--end-bin") == 0)) {
      end = atoi(argv[++arg]);
    } else if((strcmp(argv[arg], "-b") == 0) ||
	      (strcmp(argv[arg], "--width") == 0)) {
      width = atoi(argv[++arg]);
    } else if((strcmp(argv[arg], "-x") == 0) ||
	      (strcmp(argv[arg], "--overlap") == 0)) {
      overlap = atoi(argv[++arg]);
    } else if((strcmp(argv[arg], "-fs") == 0) ||
	      (strcmp(argv[arg], "--start-frequency") == 0)) {
      fMin = atof(argv[++arg]);
    } else if((strcmp(argv[arg], "-fe") == 0) ||
	      (strcmp(argv[arg], "--end-frequency") == 0)) {
      fMax = atof(argv[++arg]);
    } else if((strcmp(argv[arg], "-fb") == 0) ||
	      (strcmp(argv[arg], "--frequency-bandwidth") == 0)) {
      fWidth = atof(argv[++arg]);
    } else if((strcmp(argv[arg], "-fx") == 0) ||
	      (strcmp(argv[arg], "--frequency-overlap") == 0)) {
      fOverlap = atof(argv[++arg]);
    } else if((strcmp(argv[arg], "-m") == 0) ||
	      (strcmp(argv[arg], "--factor") == 0)) {
      factor = atof(argv[++arg]);
    } else if((strcmp(argv[arg], "-o") == 0) ||
	      (strcmp(argv[arg], "--output-prefix") == 0)) {
      prefix = argv[++arg];
    } else if((strcmp(argv[arg], "-rb") == 0) ||
	      (strcmp(argv[arg], "--read-bandwidth") == 0)) {
      read_bandwidth.resource_rate=atoi(argv[++arg]);
    } else if((strcmp(argv[arg], "-ror") == 0) ||
	      (strcmp(argv[arg], "--read-open-rate") == 0)) {
      read_open_rate.resource_rate=atoi(argv[++arg]);
    } else if((strcmp(argv[arg], "-wb") == 0) ||
	      (strcmp(argv[arg], "--write-bandwidth") == 0)) {
      write_bandwidth.resource_rate=atoi(argv[++arg]);
    } else if((strcmp(argv[arg], "-wor") == 0) ||
	      (strcmp(argv[arg], "--write-open-rate") == 0)) {
      write_open_rate.resource_rate=atoi(argv[++arg]);
    } else if((strcmp(argv[arg], "-i") == 0) ||
	      (strcmp(argv[arg], "--input-files") == 0)) {
      break;
    } else {
      fprintf(stderr, "unknown option '%s', try '-h' for help\n", argv[arg]);
      exit (-1);
    }
  }

  /* check if there was an input-file option given at all */
  TRY(argv[arg] == NULL, "no input files specified",4);
  TRY((strcmp(argv[arg], "-i") != 0) &&
      (strcmp(argv[arg], "--input-files") != 0),
      "no input files specified",5);

  /* allocate space for output filename */
  TRY((outname = (char*)malloc(strlen(prefix) + 20)) == NULL,
      "out of memory allocating outname",6);

  /* loop over all input SFT files */
  /* first skip the "-i" option */
  for(arg++; arg < argc; arg++) {

    /* open input SFT */
    request_resource(&read_open_rate, 1);
    TRY((fp = fopen(argv[arg], "r")) == NULL,
	"could not open SFT file for reading",7);
    
    /* read header */
    request_resource(&read_bandwidth, 40);
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
    if(fOverlap >= 0.0)
      overlap = MYROUND(fOverlap * hd.tbase);
 
    /* allocate space for SFT data */
    TRY((data = (float*)calloc(hd.nsamples, 2*sizeof(float))) == NULL,
	"out of memory allocating data",8);

    /* error if desired start bin < hd.firstfreqindex */
    if((int)start < hd.firstfreqindex) {
      fprintf(stderr,
	      "ERROR: start bin (%d) is smaller than first bin in input SFT (%d)\n",
	      start, hd.firstfreqindex);
      exit(9);
    }

    /* error if desired end bin > hd.firstfreqindex + hd.nsamples - 1 */
    if(start + width > end + 1) {
      fprintf(stderr,
	      "ERROR: end bin (%d) is larger than last bin in input SFT (%d)\n",
	      end, hd.firstfreqindex+hd.nsamples - 1);
      exit(10);
    }

    /* error if overlap is larger than the width */
    if(overlap >= width) {
      fprintf(stderr,
              "ERROR: overlap (%d) is not smaller than the width (%d)\n",
              overlap, width);
      exit(11);
    }

    /* construct comment for output SFTs */
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

      /* only copied existing comment, no additional space needed */
      comment = oldcomment;

    } /* else (add_comment == CMT_NONE) and (comment == NULL) i.e. no comment at all */

    /* get the detector name from SFT header if present there (in v2 SFTs),
       or else it needs to have been set on the command-line */
    if(hd.detector && *hd.detector)
      detector = hd.detector;

    /* if no detector has been specified, issue an error */
    TRY(!detector || !*detector, "When reading v1 SFTs a detector needs to be specified with -d",12);

    /* calculate number of bins to actually read (from width + overlap) */
    /* add width-overlap samples as lon as they are < the total number og bins to write */
    for(nactivesamples = 0; nactivesamples < end - start; nactivesamples += width - overlap);
    /* add the last overlap */
    nactivesamples += overlap + 1;
    /* if this number is larger than the bins in the input sft, just use the latter */
    if(nactivesamples > hd.nsamples + hd.firstfreqindex - start)
       nactivesamples = hd.nsamples + hd.firstfreqindex - start;

    /* read in SFT bins */
    request_resource(&read_bandwidth, nactivesamples*8);
    TRYSFT(ReadSFTData(fp, data, start, nactivesamples, NULL, NULL),
	   "could not read SFT data");

    /* if reading v1 SFTs set up a factor to be applied for normalization conversion */
    if(hd.version == 1.0) {
      conversion_factor = 0.5 * hd.tbase / hd.nsamples;
    } else {
      conversion_factor = 1.0;
    }

    /* apply mystery factor and possibly normalization factor */
    for(bin = 0; bin < 2 * nactivesamples; bin++)
      data[bin] *= factor * conversion_factor;

    /* close the input sfts */
    fclose(fp);

    /* loop over start bins for output SFTs */
    for(bin = start; bin < end; bin += width - overlap) {
      /* determine the number of bins actually to write from the desired 'width',
	 given that the remaining number of bin may be odd (especially from overlapping)
	 and the bins to write need to be present in the input sft
       */
      int last_input_bin   = hd.firstfreqindex + hd.nsamples - 1;
      int last_output_bin  = bin + width - 1;
      int max_input_width  = last_output_bin <= last_input_bin ? width : width - (last_output_bin - last_input_bin);
      int max_output_width = end - bin + 1;
      int this_width       = max_input_width < max_output_width ? max_input_width : max_output_width;

      /* construct filename for this output SFT */
      sprintf(outname, "%s%d", prefix, bin);

      /* append this SFT to a possible "merged" SFT with the same name */
      request_resource(&write_open_rate, 1);
      TRY((fp = fopen(outname,"a")) == NULL,
	  "could not open SFT for writing",13);

      /* write the data */
      /* write the comment only to the first SFT of a "block", i.e. of a call of this program */
      request_resource(&write_bandwidth, 40 + this_width * 8);
      TRYSFT(WriteSFT(fp, hd.gps_sec, hd.gps_nsec, hd.tbase, 
		      bin, this_width, detector,
		      (firstfile || allcomments) ? comment : NULL,
		      data + 2 * (bin - start)),
	     "could not write SFT data");

      /* close output SFT file */
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
