/*----------------------------------------------------------------------- 
 * 
 * File Name: olapredfcn_args.c
 *
 * Author: John T. Whelan
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include "olapredfcn.h"

extern int lalDebugLevel;
extern BOOLEAN optVerbose;
extern REAL8 optDeltaF;
extern UINT4 optLength;
extern REAL8 optF0;
extern UINT4 optDetector1;
extern UINT4 optDetector2;
extern REAL4 optAzimuth1;
extern REAL4 optAzimuth2;
extern CHAR optFile[LALNameLength];

extern char *optarg;
extern int   optind;

/*-----------------------------------------------------------------------*/

void
olapredfcn_usage (
    const char *program, 
    int         exitcode
    )
{
  INT4 i;

  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h             print this message\n");
  fprintf (stderr, "  -q             quiet: run silently\n");
  fprintf (stderr, "  -v             verbose: print extra information\n");
  fprintf (stderr, "  -d level       set lalDebugLevel to level\n");
  fprintf (stderr, "  -s siteID1     calculate overlap red fcn for site siteID1\n");
  fprintf (stderr, "  -t siteID2       with site siteID2\n");
  for (i=0; i<LALNumCachedDetectors; ++i)
  {
    fprintf (stderr, "                   %d = %s\n",
             i, lalCachedDetectors[i].frDetector.name);
  }
  for (i=0; i<LALNumCachedBars; ++i)
  {
    fprintf (stderr, "                   %d = %s\n",
             LALNumCachedDetectors + i, lalCachedBars[i].name);
  }
  fprintf (stderr, "  -a az     set azimuth of first detector to az degrees East of North\n");
  fprintf (stderr, "              Must be between -360 and 360;\n");
  fprintf (stderr, "              Only valid for bars (ID %d thru %d)\n",
	   LALNumCachedDetectors,
	   LALNumCachedDetectors + LALNumCachedBars - 1);
  fprintf (stderr, "  -b az     set azimuth of second detector to az degrees East of North\n");
  fprintf (stderr, "              Must be between -360 and 360;\n");
  fprintf (stderr, "              Only valid for bars (ID %d thru %d)\n",
	   LALNumCachedDetectors,
	   LALNumCachedDetectors + LALNumCachedBars - 1);
  fprintf (stderr, "  -f f0          set start frequency to f0\n");
  fprintf (stderr, "  -e deltaF      set frequency spacing to deltaF\n");
  fprintf (stderr, "  -n length      set number of points in frequency series to length\n");
  fprintf (stderr, "  -o filename    print overlap reduction function to file filename\n");
  exit (exitcode);
}

/*-----------------------------------------------------------------------*/

void
olapredfcn_parse_options (
    int         argc, 
    char       *argv[]
    )
{
  const char *dbglvl  = NULL;

  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:s:t:a:b:f:e:n:o:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'o': /* specify output file */
        strncpy (optFile, optarg, LALNameLength);
        break;
        
      case 'n': /* specify number of points in frequency series */
        optLength = atoi (optarg);
        break;
        
      case 'e': /* specify frequency resolution */
        optDeltaF = atof (optarg);
        break;
        
      case 'f': /* specify start frequency */
        optF0 = atof (optarg);
        break;

      case 'a': /* specify detector #1 azimuth */
        optAzimuth1 = atof (optarg);
	break;

      case 'b': /* specify detector #2 azimuth */
        optAzimuth2 = atof (optarg);
	break;

      case 's': /* specify detector #1 */
        optDetector1 = atoi (optarg);
	break;

      case 't': /* specify detector #2 */
        optDetector2 = atoi (optarg);
	break;

      case 'd': /* set debug level */
        dbglvl = optarg;
        set_debug_level( dbglvl );
        break;

      case 'v': /* optVerbose */
        optVerbose = OLAPREDFCNH_TRUE;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        freopen ("/dev/null", "w", stderr);
        freopen ("/dev/null", "w", stdout);
        break;

      case 'h':
        olapredfcn_usage(argv[0], 0);
        break;

      default:
        olapredfcn_usage(argv[0], 1);
    }

  }

  if (optind < argc)
  {
    olapredfcn_usage(argv[0], 1);
  }

  return;
}
