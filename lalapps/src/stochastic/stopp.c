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

#include <getopt.h>

#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>

/* cvs info */
#define PROGRAM_NAME "stopp"
#define CVS_ID "$Id$"

#define USAGE \
	"Usage: " PROGRAM_NAME " [options]\n"\
  " --help                       display this message\n"\
  " --version                    display version\n"\
  " --input FILE                 read input data from FILE\n"\
  " --output FILE                write output data to FILE\n"

INT4 main(INT4 argc, CHAR *argv[])
{
	/*  program option variables */
	CHAR *inputFileName = NULL;
	CHAR *outputFileName = NULL;
	FILE *out = NULL;

	/* xml data structures */
	INT4 numSegments = 0;
	StochasticTable *stochHead = NULL;
	StochasticTable *thisStoch = NULL;

	/* parse command line arguments */

	while (1)
	{
		/* getopt arguments */
		static struct option long_options[] = 
		{
			{"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
			{"input", required_argument, 0, 'i'},
			{"output", required_argument, 0, 'o'},
			{0, 0, 0, 0}
		};
		int c;

		/* getopt_long stores the option index here. */
		int option_index = 0;
		size_t optarg_len;

		c = getopt_long_only(argc, argv, "hvi:o:", long_options, &option_index);

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

			case 'i':
				/* create storage for the input file name */
				optarg_len = strlen(optarg) + 1;
				inputFileName = (CHAR *)calloc(optarg_len, sizeof(CHAR));
				memcpy(inputFileName, optarg, optarg_len);
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
				fprintf(stderr, "unknown error while parsing options\n");
				exit(1);
		}   
	}

	if (optind < argc)
	{
		fprintf(stderr, "extraneous command line arguments:\n");
		while (optind < argc)
		{
			fprintf (stderr, "%s\n", argv[optind++]);
		}
		exit(1);
	}

	/* read in the stochastic table */
	numSegments = StochasticTableFromLIGOLw(&stochHead, inputFileName);

	if (numSegments < 0)  
	{
		fprintf(stderr, "error: unable to read stochastic_table from %s\n", \
				inputFileName);
		exit(1);
	}
	else if (numSegments > 0)
	{
		fprintf(stdout, "Read in %d segments from file %s\n\n", numSegments, \
				inputFileName);

		/* open output file */
		if ((out = fopen(outputFileName, "w")) == NULL)
		{
			fprintf(stderr, "Can't open file \"%s\" for output...\n", outputFileName);
			exit(1);
		}

		/* write details of events */
		for(thisStoch=stochHead; thisStoch; thisStoch=thisStoch->next)
		{
			fprintf(out, "%d %e %e\n", thisStoch->start_time.gpsSeconds, \
					thisStoch->cc_stat, thisStoch->cc_sigma);
		}
	}

	/* close output file */
	fclose(out);

	/* free stochastic table */
	while (stochHead)
	{
		thisStoch = stochHead;
		stochHead = stochHead->next;
		LALFree(thisStoch);
	}

	LALCheckMemoryLeaks();
	exit(0);
}

/*
 * vim: et
 */
