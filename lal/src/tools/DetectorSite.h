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

/*
 * author John T. Whelan
 * date 2007
 * Dummy header file which includes LALDetectors.h for backwards compatibility
 *
 */

#ifndef _DETECTORSITE_H
#define _DETECTORSITE_H

#include <lal/LALDetectors.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/* Legacy code: should now use data in LALDetectors.h */
enum
{
  LALDetectorIndexLHODIFF = LAL_LHO_4K_DETECTOR,
  LALDetectorIndexLLODIFF = LAL_LLO_4K_DETECTOR,
  LALDetectorIndexVIRGODIFF = LAL_VIRGO_DETECTOR,
  LALDetectorIndexGEO600DIFF = LAL_GEO_600_DETECTOR,
  LALDetectorIndexTAMA300DIFF = LAL_TAMA_300_DETECTOR,
  LALDetectorIndexCIT40DIFF = LAL_CIT_40_DETECTOR,
  LALDetectorIndexE1DIFF = LAL_ET1_DETECTOR,
  LALDetectorIndexE2DIFF = LAL_ET2_DETECTOR,
  LALDetectorIndexE3DIFF = LAL_ET3_DETECTOR,
  LALDetectorIndexKAGRADIFF = LAL_KAGRA_DETECTOR,
  LALNumCachedDetectors = LAL_NUM_DETECTORS
};

#ifdef  __cplusplus
}
#endif

#endif /* _DETECTORSITE_H */
