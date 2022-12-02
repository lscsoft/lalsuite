/*
 * Copyright (C) 2006 Reinhard Prix
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/**
 * \author Reinhard Prix
 * \date 2006
 * \file
 * \ingroup lalpulsar_general
 * \brief Header-file defining the API for the LISA-specific functions
 *
 */

#ifndef _LISASPECIFICS_H  /* Double-include protection. */
#define _LISASPECIFICS_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- exported INCLUDES ----------*/
#include <lal/LALDetectors.h>
#include <lal/LALDatatypes.h>
#include <lal/DetectorStates.h>
#include <lal/LALComputeAM.h>

/*---------- exported DEFINES ----------*/

#define LISA_TIME_ORIGIN		700000000; 	/* ad-hoc setting for MLDC 't=0' to GPS for Tue Mar 12 20:26:27 GMT 2002 */

/*---------- exported types ----------*/
/** Translate TDI arm indices to C-indexing */
typedef enum tagLISAarmT {
  LISA_ARM1 = 0,
  LISA_ARM2,
  LISA_ARM3
} LISAarmT;


/**
 * The 'detector tensor' for a GW-detector: symmetric 3x3 matrix, storing only the upper triangle.
 * The coordinate-system is SSB-fixed Cartesian coordinates, in particular EQUATORIAL coords for
 * Earth-based detectors and ECLIPTIC coords for LISA.
 */
typedef struct tagCmplxDetectorTensor
{
  SymmTensor3 re;	/**< tensor holding real-parts of all components */
  SymmTensor3 im;	/**< tensor holding imaginary-parts of all components */
} CmplxDetectorTensor;

/**
 * Convenience container for precomputed pi f L/c  and skyposition vector
 */
typedef struct tagFreqSkypos_t
{
  REAL4 Freq;		/**< signal frequency */
  REAL8 skyposV[3];	/**< unit vector pointing to skyposition of source */
  SymmTensor3 ePlus;	/**< ePlus polarization tensor (skypos-dependent) */
  SymmTensor3 eCross;	/**< eCross polarization tensor (skypos-dependent) */
} FreqSkypos_t;

/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */

/*---------- exported prototypes [API] ----------*/
BOOLEAN XLALisLISAdetector ( const LALDetector *det );
int XLALregisterLISAdetectors ( const CHAR prefixLetter );

int XLALprecomputeLISAarms ( DetectorState *detState );

int XLALgetLISADetectorTensorLWL ( SymmTensor3 *detT, const Detector3Arms detArms, CHAR channelNum );
int XLALgetLISADetectorTensorRAA ( CmplxDetectorTensor *detT, const Detector3Arms detArms, CHAR channelNum, const FreqSkypos_t *freq_skypos );

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
