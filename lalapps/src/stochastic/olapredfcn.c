/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
*                2010 Nickolas Fotopoulos
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
 * 
 *-----------------------------------------------------------------------
 */

#include "olapredfcn.h"

extern char *optarg;
extern int   optind;

BOOLEAN optVerbose  = OLAPREDFCNH_FALSE;
REAL8 optDeltaF     = -1;
UINT4 optLength     = 0;
REAL8 optF0         = 0.0;
UINT4 optDetector1  = LALNumCachedDetectors;
UINT4 optDetector2  = LALNumCachedDetectors;
REAL4 optAzimuth1   = OLAPREDFCNH_OOR;
REAL4 optAzimuth2   = OLAPREDFCNH_OOR;
CHAR optFile[LALNameLength] = "";

int main( int argc, char *argv[] )
{
  LALStatus                status = blank_status;
  
  OverlapReductionFunctionParameters   parameters;
  REAL4FrequencySeries     overlap;

  LALDetectorPair     detectors;

  lal_errhandler = LAL_ERR_EXIT;

  olapredfcn_parse_options( argc, argv );

  if ( !optFile[0] ) {
    fprintf( stderr, "ERROR: No output file specified\n" );
    return OLAPREDFCNH_EARG;
  }

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

  if (optDetector1 >= LALNumCachedDetectors) {
    fprintf( stderr, "ERROR: Detector 1 unknown\n");
    return OLAPREDFCNH_EARG;
  }
  if (optDetector2 >= LALNumCachedDetectors) {
    fprintf( stderr, "ERROR: Detector 2 unknown\n");
    return OLAPREDFCNH_EARG;
  }

  if ( (( optAzimuth1 != OLAPREDFCNH_OOR ) &&
        (( optAzimuth1 < -360.0 ) || ( optAzimuth1 > 360.0 ))) ||
       (( optAzimuth2 != OLAPREDFCNH_OOR ) &&
        (( optAzimuth2 < -360.0 ) && ( optAzimuth2 > 360.0 )))) {
    fprintf( stderr, "ERROR: azimuth must be between -360 and 360\n" );
    return OLAPREDFCNH_EARG;
  }

  parameters.length = optLength;
  parameters.f0 = optF0;
  parameters.deltaF = optDeltaF;

  /* set detector from known detectors */
  detectors.detectorOne = lalCachedDetectors[optDetector1];
  detectors.detectorTwo = lalCachedDetectors[optDetector2];
  if ( optVerbose ) {
    fprintf( stderr, "Using cached detector %s for detector 1\n",
      detectors.detectorOne.frDetector.name );
    fprintf( stderr, "Using cached detector %s for detector 2\n",
      detectors.detectorTwo.frDetector.name );
  }

  /* bars are uniquely rotatable */
  if ( (( detectors.detectorOne.type != LALDETECTORTYPE_CYLBAR ) &&
        ( optAzimuth1 != OLAPREDFCNH_OOR )) ||
       (( detectors.detectorTwo.type != LALDETECTORTYPE_CYLBAR ) &&
        ( optAzimuth2 != OLAPREDFCNH_OOR )) ) {
    fprintf( stderr, "ERROR: Can only set azimuth for bar detectors\n" );
    return OLAPREDFCNH_EARG;
  }

  if ( detectors.detectorOne.type == LALDETECTORTYPE_CYLBAR ) {
    if ( optAzimuth1 != OLAPREDFCNH_OOR ) {
      if ( optVerbose ) {
        fprintf( stderr, "Changing azimuth to %.3f degrees East of North\n",
          optAzimuth1 );
      }
      detectors.detectorOne.frDetector.xArmAzimuthRadians =
        optAzimuth1 * LAL_PI_180;
    }
    else if ( optVerbose ) {
      fprintf( stderr, "Using IGEC azimuth of %.3f degrees East of North\n",
        detectors.detectorOne.frDetector.xArmAzimuthRadians * LAL_180_PI );
    }
  }
  if ( detectors.detectorTwo.type == LALDETECTORTYPE_CYLBAR ) {
    if ( optAzimuth2 != OLAPREDFCNH_OOR ) {
      if ( optVerbose ) {
        fprintf( stderr, "Changing azimuth to %.3f degrees East of North\n",
          optAzimuth2 );
      }
      detectors.detectorTwo.frDetector.xArmAzimuthRadians =
        optAzimuth2 * LAL_PI_180;
    }
    else if ( optVerbose ) {
      fprintf( stderr, "Using IGEC azimuth of %.3f degrees East of North\n",
        detectors.detectorTwo.frDetector.xArmAzimuthRadians * LAL_180_PI );
    }
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
