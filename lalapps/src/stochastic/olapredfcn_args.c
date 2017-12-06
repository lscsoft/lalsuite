/*
*  Copyright (C) 2007 Jolien Creighton
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

/*----------------------------------------------------------------------- 
 * 
 * File Name: olapredfcn_args.c
 *
 * Author: John T. Whelan
 * 
 * 
 *-----------------------------------------------------------------------
 */

#include "olapredfcn.h"

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
  fprintf (stderr, "  -s siteID1     calculate overlap red fcn for site siteID1\n");
  fprintf (stderr, "  -t siteID2       with site siteID2\n");
  for (i=0; i<LALNumCachedDetectors; ++i)
  {
    fprintf (stderr, "                   %d = %s\n",
             i, lalCachedDetectors[i].frDetector.name);
  }
  fprintf (stderr, "  -a az     set azimuth of first detector to az degrees East of North\n");
  fprintf (stderr, "              Must be between -360 and 360;\n");
  fprintf (stderr, "              Only valid for bars \n");
  fprintf (stderr, "  -b az     set azimuth of second detector to az degrees East of North\n");
  fprintf (stderr, "              Must be between -360 and 360;\n");
  fprintf (stderr, "              Only valid for bars\n");
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
  FILE *fp;

  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvs:t:a:b:f:e:n:o:");
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

      case 'v': /* optVerbose */
        optVerbose = OLAPREDFCNH_TRUE;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        fp = freopen("/dev/null", "w", stderr);
        if ( fp == NULL )
        {
          XLALPrintError( "XLAL Error - %s: unable to open /dev/null\n", __func__);
          XLAL_ERROR_VOID( XLAL_EFAILED );
        }
        fp = freopen("/dev/null", "w", stdout);
        if ( fp == NULL )
        {
          XLALPrintError( "XLAL Error - %s: unable to open /dev/null\n", __func__);
          XLAL_ERROR_VOID( XLAL_EFAILED );
        }
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
