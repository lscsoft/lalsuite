/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
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
 * File Name: olapredfcn.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include "olapredfcn.h"

RCSID("$Id$");

extern char *optarg;
extern int   optind;

BOOLEAN optVerbose  = OLAPREDFCNH_FALSE;
REAL8 optDeltaF     = -1;
UINT4 optLength     = 0;
REAL8 optF0         = 0.0;
UINT4 optDetector1  = LALNumCachedDetectors + LALNumCachedBars;
UINT4 optDetector2  = LALNumCachedDetectors + LALNumCachedBars;
REAL4 optAzimuth1   = OLAPREDFCNH_OOR;
REAL4 optAzimuth2   = OLAPREDFCNH_OOR;
CHAR optFile[LALNameLength] = "";

int main( int argc, char *argv[] )
{
  LALStatus                status = blank_status;
  
  OverlapReductionFunctionParameters   parameters;
  REAL4FrequencySeries     overlap;

  LALDetectorPair     detectors;
  LALFrDetector       barGeom;

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  olapredfcn_parse_options( argc, argv );

  if (optFile[0]) {

    if (optLength == 0) {
      fprintf( stderr, "ERROR: Zero length specified\n");
      return OLAPREDFCNH_EARG;
    }

    if (optDeltaF <= 0) {
      fprintf( stderr, "ERROR: Non-positive frequency spacing specified\n");
      return OLAPREDFCNH_EARG;
    }

    if (optF0 < 0) {
      fprintf( stderr, "ERROR: Negative start frequency specified\n");
      return OLAPREDFCNH_EARG;
    }

    parameters.length = optLength;
    parameters.f0 = optF0;
    parameters.deltaF = optDeltaF;

    if (optDetector1 < LALNumCachedDetectors)
    {
      detectors.detectorOne = lalCachedDetectors[optDetector1];
      if ( optVerbose ) {
	fprintf( stderr, "Using cached detector %s for detector 1\n", 
		 detectors.detectorOne.frDetector.name );
      }
    }
    else if (optDetector1 < LALNumCachedDetectors + LALNumCachedBars) {
      barGeom = lalCachedBars[optDetector1-LALNumCachedDetectors];
      if ( optVerbose ) {
	fprintf( stderr, "Creating detector 1 from known bar detector %s\n", 
		 barGeom.name );
      }
      if ( optAzimuth1 > -360.0 && optAzimuth1 < 360.0 ) {
	if ( optVerbose ) {
	  fprintf( stderr, "Changing azimuth to %.3f degrees East of North\n", 
		   optAzimuth1 );
	}
	barGeom.xArmAzimuthRadians = 
	  optAzimuth1 * LAL_PI_180;
      } 
      else if ( optVerbose ) {
	fprintf( stderr, "Using IGEC azimuth of %.3f degrees East of North\n",
		 (REAL4)
		 ( barGeom.xArmAzimuthRadians * LAL_180_PI )
		 );
      }

      LAL_CALL( 
	 LALCreateDetector( &status, &(detectors.detectorOne),
			    &barGeom, LALDETECTORTYPE_CYLBAR ),
	 &status );
    } else {
      fprintf( stderr, "ERROR: Detector 1 ID %d out of range\n",
	       optDetector1 );
      return OLAPREDFCNH_EARG;
    }

    if (optDetector2 < LALNumCachedDetectors)
    {
      detectors.detectorTwo = lalCachedDetectors[optDetector2];
      if ( optVerbose ) {
	fprintf( stderr, "Using cached detector %s for detector 2\n", 
		 detectors.detectorTwo.frDetector.name );
      }
    }
    else if (optDetector2 < LALNumCachedDetectors + LALNumCachedBars) {
      barGeom = lalCachedBars[optDetector2-LALNumCachedDetectors];
      if ( optVerbose ) {
	fprintf( stderr, "Creating detector 2 from known bar detector %s\n", 
		 barGeom.name );
      }
      if ( optAzimuth2 > -360.0 && optAzimuth2 < 360.0 ) {
	if ( optVerbose ) {
	  fprintf( stderr, "Changing azimuth to %.3f degrees East of North\n", 
		   optAzimuth2 );
	}
	barGeom.xArmAzimuthRadians = 
	  optAzimuth2 * LAL_PI_180;
      } 
      else if ( optVerbose ) {
	fprintf( stderr, "Using IGEC azimuth of %.3f degrees East of North\n",
		 (REAL4)
		 ( barGeom.xArmAzimuthRadians * LAL_180_PI )
		 );
      }

      LAL_CALL( 
	 LALCreateDetector( &status, &(detectors.detectorTwo),
			    &barGeom, LALDETECTORTYPE_CYLBAR ),
	 &status );
    } else {
      fprintf( stderr, "ERROR: Detector 2 ID %d out of range\n",
	       optDetector2 );
      return OLAPREDFCNH_EARG;
    }

    overlap.data = NULL;
    LAL_CALL( 
       LALSCreateVector(&status, &(overlap.data), optLength),
       &status );
    LAL_CALL( 
       LALOverlapReductionFunction(&status, &overlap, &detectors, &parameters),
       &status );
    LAL_CALL( 
       LALPrintFrequencySeries( &overlap, optFile ),
       &status );

 if ( optVerbose ) {
   fprintf(stderr, 
	   "======== Overlap Reduction Function Written to File %s ========\n",
	   optFile);
 }

    LAL_CALL( 
       LALSDestroyVector(&status, &(overlap.data)),
       &status );
    
    LALCheckMemoryLeaks();
    return OLAPREDFCNH_ENOM;
  }
  else 
  {
    fprintf( stderr, "ERROR: No output file specified\n" );
    return OLAPREDFCNH_EARG;
  }
}
