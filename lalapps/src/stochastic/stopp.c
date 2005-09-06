/*
 * stopp.c - Stochastic Post Processing
 *
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 *
 * $Id$
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <getopt.h>

#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>

#include <lalapps.h>

/* verbose flag */
extern int vrbflg;

NRCSID(STOPPC, "$Id$");
RCSID("$Id$");

/* cvs info */
#define PROGRAM_NAME "stopp"
#define CVS_ID "$Id$"

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

  /* getopt flags */
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
    /* getopt arguments */
    static struct option long_options[] = 
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

    /* getopt_long stores the option index here. */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only(argc, argv, "hvo:", long_options, &option_index);

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
              long_options[option_index].name, optarg);
          exit(1);
        }
        break;

      case 'h':
        fprintf(stdout, USAGE);
        exit(0);
        break;

      case 'v':
        /* display version info and exit */
        fprintf(stdout, "Stochastic Post Processing\n" CVS_ID "\n");
        exit(0);
        break;

      case 'o':
        /* create storage for the output file name */
        optarg_len = strlen(optarg) + 1;
        outputFileName = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(outputFileName, optarg, optarg_len);
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
  if (optind < argc)
  {
    for (i = optind; i < argc; ++i)
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
