/*----------------------------------------------------------------------- 
 * 
 * File Name: findchirp_args.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include "findchirp.h"

extern char *optarg;
extern int optind, opterr, optopt;

extern int vrbflg;

extern UINT4   numPoints;
extern UINT4   numSegments;
extern UINT4   ovrlap;

extern INT4    invSpecTrunc;
extern REAL4   fLow;
extern REAL4   dynRange;

extern INT4    numChisqBins;
extern REAL4   rhosqThreshVec[];
extern REAL4   chisqThreshVec[];


/*-----------------------------------------------------------------------*/

void
findchirp_usage (
    const char *program, 
    int         exitcode
    )
{
  fprintf( stderr, "Usage: %s [options]\n"
      "Options (defauls shown in square brackets):\n"
      "  General:\n"
      "    -h                print this message\n"
      "    -v                verbose\n"
      "    -d debuglevel     LAL status debug level [1]\n"
      "  Data Parameters:\n"
      "    -n numPoints      number of points in a segment [1048576]\n"
      "    -s numSegments    number of data segments [1]\n"
      "    -o ovrlap         overlap beteen of data segments [0]\n"
      "  Data Conditioning:\n"
      "    -i invSpecTrunc   number of points to truncate S^{-1}(f) [262144]\n"
      "    -f fLow           low frequency cutoff for S^{-1}(f) [40.0]\n"
      "    -y dynRange       log_2( dynamicRange ) [69.0]\n"
      "  Filtering:\n"
      "    -b numChisqBins   number of bins for chi squared test [8]\n"
      "    -t rhosqThresh    signal to noise squared threshold [10.0,0.0]\n"
      "    -c chisqThresh    chi squared threshold [5.0,0.0]\n", program);

  exit (exitcode);
}

/*-----------------------------------------------------------------------*/

void
findchirp_parse_options (
    int         argc, 
    char       *argv[]
    )
{
  const char *dbglvl  = NULL;

  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, 
        "hvd:"                          /* general              */
        "n:""s:""o:"                    /* data params          */
        "i:""f:""y:"                    /* data conditioning    */
        "b:""t:""c:"                    /* filtering            */
        );
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'h':
        findchirp_usage (argv[0], 0);
        break;
      case 'd': /* set debuglevel */
        dbglvl = optarg;
        set_debug_level( dbglvl );
        break;
      case 'v': /* set verbosity */
        vrbflg = 1;
        break;

      case 'n': /* set number of points in a segment */
        numPoints = atoi (optarg);
        break;
      case 's': /* set number of segments */
        numSegments = atoi (optarg);
        break;
      case 'o': /* set sampling rate */
        ovrlap = atoi (optarg);
        break;

      case 'i': /* set invSpecTrunc */
        invSpecTrunc = atoi (optarg);
        break;
      case 'f': /* set fLow */
        fLow = (REAL4) atof (optarg);
        break;
      case 'y': /* set dynRange */
        dynRange = (REAL4) atof (optarg);
        break;

      case 'b': /* set numChisqBins */
        numChisqBins = atoi (optarg);
        break;
      case 't': /* set rhosq threshold */
        if ( sscanf( optarg, "%f,%f", 
              rhosqThreshVec, rhosqThreshVec + 1 ) - 2 )
        {
          fprintf( stderr, "error reading signal to noise thresholds\n" );
          exit( 1 );
        }
        break;
      case 'c': /* set rhosq threshold */
        if ( sscanf( optarg, "%f,%f", 
              chisqThreshVec, chisqThreshVec + 1 ) - 2 )
        {
          fprintf( stderr, "error reading chi squared thresholds\n" );
          exit( 1 );
        }
        break;
      default:
        findchirp_usage (argv[0], 1);
    }
  }

  if (optind < argc)
  {
    findchirp_usage (argv[0], 1);
  }

  return;
}

/*-----------------------------------------------------------------------*/

void
findchirp_print_options ( 
    void 
    )
{
  fprintf( stdout, 
      "debug level\t\t\t\t= %d\n\n"
      "number of points per data segment\t= %d\n"
      "number of data segments\t\t\t= %d\n"
      "data segments overlap\t\t\t= %d\n\n"
      "inverse power sectrum truncation\t= %d\n"
      "low frequency cutoff\t\t\t= %e\n"
      "log_2( dynamic range )\t\t\t= %e\n\n"
      "number of chisq bins\t\t\t= %d\n"
      "rhosq thresholds\t\t\t= %e\t%e\n"
      "chisq thresholds\t\t\t= %e\t%e\n",
      lalDebugLevel,
      numPoints,
      numSegments,
      ovrlap,
      invSpecTrunc,
      fLow,
      dynRange,
      numChisqBins,
      rhosqThreshVec[0], rhosqThreshVec[1],
      chisqThreshVec[0], chisqThreshVec[1]
      );

  return;
}


