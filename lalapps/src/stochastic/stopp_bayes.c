/*
 * stopp_bayes.c - Stochastic Post Processing for the Bayesian Search
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

#include <gsl/gsl_sf_erf.h>

#include <lalapps.h>

/* verbose flag */
extern int vrbflg;

NRCSID(STOPPBAYESC, "$Id$");
RCSID("$Id$");

/* cvs info */
#define PROGRAM_NAME "stopp_bayes"
#define CVS_ID "$Id$"

#define USAGE \
  "Usage: " PROGRAM_NAME " [options] [xml files]\n"\
  " --help                       display this message\n"\
  " --version                    display version\n"\
  " --verbose                    verbose mode\n"\
  " --cat-only                   only cat xml files\n"\
  " --text                       output file as text\n"\
  " --output FILE                write output data to FILE\n"\
  " --confidence LEVEL           set confidence to LEVEL\n"

/* helper functions */

/* function to return the inverse complimentary error function */
static double stopp_erfcinv(double y)
{
  /*
   * based on dierfc() by Takuya OOURA:
   *
   * http://momonga.t.u-tokyo.ac.jp/~ooura/gamerf.html
   * 
   * Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
   * You may use, copy, modify this code for any purpose and
   * without fee. You may distribute this ORIGINAL package.
   */
  double s, t, u, w, x, z;

  z = y;

  if (y > 1)
    z = 2 - y;

  w = 0.916461398268964 - log(z);
  u = sqrt(w);
  s = (log(u) + 0.488826640273108) / w;
  t = 1 / (u + 0.231729200323405);
  x = u * (1 - s * (s * 0.124610454613712 + 0.5)) - \
      ((((-0.0728846765585675 * t + 0.269999308670029) * t + \
      0.150689047360223) * t + 0.116065025341614) * t + \
      0.499999303439796) * t;
  t = 3.97886080735226 / (x + 3.97886080735226);
  u = t - 0.5;
  s = (((((((((0.00112648096188977922 * u + \
      1.05739299623423047e-4) * u - 0.00351287146129100025) * u - \
      7.71708358954120939e-4) * u + 0.00685649426074558612) * u + \
      0.00339721910367775861) * u - 0.011274916933250487) * u - \
      0.0118598117047771104) * u + 0.0142961988697898018) * u + \
      0.0346494207789099922) * u + 0.00220995927012179067;
  s = ((((((((((((s * u - 0.0743424357241784861) * u - \
      0.105872177941595488) * u + 0.0147297938331485121) * u + \
      0.316847638520135944) * u + 0.713657635868730364) * u + \
      1.05375024970847138) * u + 1.21448730779995237) * u + \
      1.16374581931560831) * u + 0.956464974744799006) * u + \
      0.686265948274097816) * u + 0.434397492331430115) * u + \
      0.244044510593190935) * t - \
      z * exp(x * x - 0.120782237635245222);
  x += s * (x * s + 1);

  if (y > 1)
    x = -x;

  return x;
}

INT4 main(INT4 argc, CHAR *argv[])
{
  /* status */
  LALStatus status = blank_status;

  /* getopt flags */
  static int text_flag;
  static int cat_flag;

  /* counters */
  INT4 i;

  /* combined statistics variables */
  REAL8 numerator = 0;
  REAL8 denominator = 0;
  REAL8 yOpt = 0;
  REAL8 sigmaOpt = 0;
  REAL8 confidence = 0.95;
  REAL8 zeta;
  REAL8 upperlimit;

  /* pdf */
  REAL8 exponent;
  REAL8 pdf[100];
  REAL8 min_omega;
  REAL8 max_omega;
  REAL8 omega_i;

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
  FILE *pdf_out;

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
      /* options that don't set a flag */
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {"output", required_argument, 0, 'o'},
      {"confidence", required_argument, 0, 'c'},
      {0, 0, 0, 0}
    };
    int c;

    /* getopt_long stores the option index here. */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only(argc, argv, "hvo:c:", long_options, &option_index);

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
        fprintf(stdout, "Stochastic Post Processing: Bayesian\n" CVS_ID "\n");
        exit(0);
        break;

      case 'o':
        /* create storage for the output file name */
        optarg_len = strlen(optarg) + 1;
        outputFileName = (CHAR *)calloc(optarg_len, sizeof(CHAR));
        memcpy(outputFileName, optarg, optarg_len);
        break;

      case 'c':
        /* confidence level */
        confidence = atof(optarg);
        if ((confidence >= 1) || (confidence <= 0))
        {
          fprintf(stderr, "invalid argument to --%s\n" \
              "confidence must be between 0 and 1, exclusive " \
              "(%.2f specified)\n", long_options[option_index].name, \
              confidence);
          exit(1);
        }
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
  }

  /* calculate upperlimit */
  zeta = yOpt / (sqrt(2) * sigmaOpt);
  upperlimit = yOpt + (sqrt(2) * sigmaOpt * \
      stopp_erfcinv((1 - confidence) * gsl_sf_erfc(-zeta)));
  fprintf(stdout, "upperlimit = %e\n", upperlimit);

  /* calculate pdf */
  min_omega = 0;
  max_omega = yOpt + (3 * sigmaOpt);
  for (i = 0; i < 100; i++)
  {
    omega_i = min_omega + (((i - 1)/99.) * (max_omega - min_omega));
    exponent = ((omega_i - yOpt) / sigmaOpt) * ((omega_i - yOpt) / sigmaOpt);
    pdf[i] = exp(-0.5 * exponent);
  }

  /* save out pdf */
  if ((pdf_out = fopen("pdf.dat", "w")) == NULL)
  {
    fprintf(stderr, "Can't open file for pdf output...\n");
    exit(1);
  }
  for (i = 0; i < 100; i++)
  {
    omega_i = min_omega + (((i - 1)/99.) * (max_omega - min_omega));
    fprintf(pdf_out, "%e %e\n", omega_i, pdf[i]);
  }
  fclose(pdf_out);
  
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
      LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlStream, stochastic_table), \
          &status);
      LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlStream, outputTable, \
            stochastic_table), &status);
      LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlStream), &status);
    }

    /* close xml file */
    LAL_CALL(LALCloseLIGOLwXMLFile(&status, &xmlStream), &status);
  }

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  exit(0);
}

/*
 * vim: et
 */
