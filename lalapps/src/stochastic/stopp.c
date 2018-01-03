/*
 * stopp.c - SGWB Standalone Analysis Pipeline
 *         - Post Processing
 *
 * Copyright (C) 2004-2006,2009,2010 Adam Mercer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

/**
 * \file
 * \ingroup lalapps_inspiral
 *
 *
 * <dl>
 * <dt>Name</dt><dd>
 * <tt>lalapps_stopp</tt> --- Stochastic Post Processing.</dd>
 *
 * <dt>Synopsis</dt><dd>
 * <tt>lalapps_stopp</tt> <i>options</i> <i>xml files</i>
 * <tt>--help</tt>
 * <tt>--version</tt>
 * <tt>--verbose</tt>
 * <tt>--cat-only</tt>
 * <tt>--analyse-only</tt>
 * <tt>--text</tt>
 * <tt>--output</tt> <i>FILE</i></dd>
 *
 * <dt>Description</dt><dd>
 * <tt>lalapps_stopp</tt> performs post processing upon output from
 * <tt>lalapps_stochastic</tt>.</dd>
 *
 * <dt>Options</dt><dd>
 * <dl>
 * <dt><tt>--help</tt></dt><dd>
 * Display usage information</dd>
 * <dt><tt>--version</tt></dt><dd>
 * Display version information</dd>
 * <dt><tt>--verbose</tt></dt><dd>
 * Verbose mode</dd>
 * <dt><tt>--cat-only</tt></dt><dd>
 * Only cat XML files together</dd>
 * <dt><tt>--analyse-only</tt></dt><dd>
 * Only combine statistics</dd>
 * <dt><tt>--text</tt></dt><dd>
 * Output file as text</dd>
 * <dt><tt>--output</tt> <i>FILE</i></dt><dd>
 * write output data to <i>FILE</i></dd>
 * </dl></dd>
 *
 * <dt>Example</dt><dd>
 * <tt>lalapps_stopp</tt> is generally run as part of a DAG, as created by
 * the pipeline generation script <tt>lalapps_stochastic_pipe</tt>, however
 * an example usage can be seen below.
 *
 * \code
 * > lalapps_stopp --output S3-H1L1-STOCHASTIC.xml \
 * >   H1L1-STOCHASTIC-753601044-753601242.xml \
 * >   H1L1-STOCHASTIC-753620042-753620352.xml \
 * >   H1L1-STOCHASTIC-753638864-753639462.xml \
 * >   H1L1-STOCHASTIC-753785374-753785707.xml \
 * >   H1L1-STOCHASTIC-753791744-753792342.xml
 * \endcode</dd>
 *
 * <dt>Author</dt><dd>
 * Adam Mercer</dd>
 * </dl>
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <lal/Date.h>
#include <lal/LALgetopt.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLStochasticRead.h>
#include <lal/LIGOMetadataTables.h>

#include <lalapps.h>
#include <LALAppsVCSInfo.h>

/* verbose flag */
extern int vrbflg;


/* cvs info */
#define PROGRAM_NAME "stopp"

#define USAGE \
  "Usage: " PROGRAM_NAME " [options] [xml files]\n"\
  " --help                       display this message\n"\
  " --version                    display version\n"\
  " --verbose                    verbose mode\n"\
  " --cat-only                   only cat xml files\n"\
  " --analyse-only               only combine statistics\n"\
  " --text                       output file as text\n"\
  " --output FILE                write output data to FILE\n"

INT4 main(INT4 argc, CHAR *argv[])
{
  /* status */
  LALStatus status = blank_status;

  /* LALgetopt flags */
  static int text_flag;
  static int cat_flag;
  static int analyse_flag;

  /* counters */
  INT4 i;

  /* combined statistics variables */
  REAL8 numerator = 0;
  REAL8 denominator = 0;
  REAL8 yOpt = 0;
  REAL8 sigmaOpt = 0;

  /* program option variables */
  CHAR *outputFileName = NULL;

  /* xml data structures */
  LIGOLwXMLStream xmlStream;
  INT4 numSegments = 0;
  StochasticTable *stochHead = NULL;
  StochasticTable *thisStoch = NULL;
  MetadataTable outputTable;
  StochasticTable **stochHandle = NULL;

  /* text output file */
  FILE *out;

  /* parse command line arguments */
  while (1)
  {
    /* LALgetopt arguments */
    static struct LALoption long_options[] =
    {
      /* options that set a flag */
      {"verbose", no_argument, &vrbflg, 1},
      {"text", no_argument, &text_flag, 1},
      {"cat-only", no_argument, &cat_flag, 1},
      {"analyse-only", no_argument, &analyse_flag, 1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {"output", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };
    int c;

    /* LALgetopt_long stores the option index here. */
    int option_index = 0;
    size_t LALoptarg_len;

    c = LALgetopt_long_only(argc, argv, "hvo:", long_options, &option_index);

    /* detect the end of the options */
    if (c == - 1)
    {
      /* end of options, break loop */
      break;
    }

    switch (c)
    {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
        {
          break;
        }
        else
        {
          fprintf(stderr, "error parseing option %s with argument %s\n", \
              long_options[option_index].name, LALoptarg);
          exit(1);
        }
        break;

      case 'h':
        fprintf(stdout, USAGE);
        exit(0);
        break;

      case 'v':
        /* display version info and exit */
        fprintf(stdout, "Stochastic Post Processing\n");
        XLALOutputVersionString(stderr,0);
        exit(0);
        break;

      case 'o':
        /* create storage for the output file name */
        LALoptarg_len = strlen(LALoptarg) + 1;
        outputFileName = (CHAR *)calloc(LALoptarg_len, sizeof(CHAR));
        memcpy(outputFileName, LALoptarg, LALoptarg_len);
        break;

      case '?':
        exit(1);
        break;

      default:
        fprintf(stderr, "Unknown error while parsing options\n");
        exit(1);
    }
  }

  /* read in the input data from the rest of the arguments */
  if (LALoptind < argc)
  {
    for (i = LALoptind; i < argc; ++i)
    {
      struct stat infileStatus;

      /* if the named file does not exist, exit with an error */
      if (stat(argv[i], &infileStatus) == -1)
      {
        fprintf(stderr, "Error opening input file \"%s\"\n", argv[i]);
        exit(1);
      }

      if (!stochHead)
      {
        stochHandle = &stochHead;
      }
      else
      {
        stochHandle = &thisStoch->next;
      }

      /* read in the stochastic table */
      numSegments = LALStochasticTableFromLIGOLw(stochHandle, argv[i]);

      if (numSegments < 0)
      {
        fprintf(stderr, "Unable to read stochastic_table from \"%s\"\n", \
            argv[i]);
        exit(1);
      }
      else if (numSegments > 0)
      {
        if (vrbflg)
        {
          fprintf(stdout, "Read in %d segments from file \"%s\"\n", \
              numSegments, argv[i]);
        }

        /* scroll to end of list */
        thisStoch = *stochHandle;
        while (thisStoch->next)
        {
          thisStoch = thisStoch->next;
        }
      }
    }
  }

  if (!cat_flag)
  {
    /* combine statistics */
    for (thisStoch = stochHead; thisStoch; thisStoch = thisStoch->next)
    {
      numerator += thisStoch->cc_stat / (thisStoch->cc_sigma * \
          thisStoch->cc_sigma);
      denominator += 1./(thisStoch->cc_sigma * thisStoch->cc_sigma);
    }
    yOpt = (1./stochHead->duration.gpsSeconds) * (numerator / denominator);
    sigmaOpt = (1./stochHead->duration.gpsSeconds) * (1./sqrt(denominator));

    /* print results */
    fprintf(stdout, "    yOpt = %e\n", yOpt);
    fprintf(stdout, "sigmaOpt = %e\n", sigmaOpt);
  }

  if (!analyse_flag)
  {
    /* output as text file */
    if (text_flag)
    {
      /* open output file */
      if ((out = fopen(outputFileName, "w")) == NULL)
      {
        fprintf(stderr, "Can't open file \"%s\" for output...\n", \
            outputFileName);
        exit(1);
      }

      /* write details of events */
      for (thisStoch = stochHead; thisStoch; thisStoch = thisStoch->next)
      {
        fprintf(out, "%d %e %e\n", thisStoch->start_time.gpsSeconds, \
            thisStoch->cc_stat, thisStoch->cc_sigma);
      }

      /* close output file */
      fclose(out);
    }

    /* output as xml file */
    else
    {
      /* open xml file stream */
      memset(&xmlStream, 0, sizeof(LIGOLwXMLStream));
      LAL_CALL(LALOpenLIGOLwXMLFile(&status, &xmlStream, outputFileName), \
          &status);

      /* write stochastic table */
      if (stochHead)
      {
        outputTable.stochasticTable = stochHead;
        LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlStream, \
              stochastic_table), &status);
        LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlStream, outputTable, \
              stochastic_table), &status);
        LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlStream), &status);
      }

      /* close xml file */
      LAL_CALL(LALCloseLIGOLwXMLFile(&status, &xmlStream), &status);
    }
  }

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  exit(0);
}

/*
 * vim: et
 */
