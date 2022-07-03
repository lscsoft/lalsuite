/*
*  Copyright (C) 2007 Drew Keppel, Duncan Brown, Jolien Creighton, Patrick Brady, Stephen Fairhurst, Thomas Cokelaer
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

/*-----------------------------------------------------------------------
 *
 * File Name: LIGOMetadataUtils.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A.
 * \file
 * \ingroup lalmetaio_general
 * \brief General routines for manipulating LIGO metadatabase tables.
 *
 * ### Description ###
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * LALCalloc, LALMalloc, LALFree.
 *
 * ### Notes ###
 *
 * %% Any relevant notes.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>


void
XLALReturnIFO(
    char                *ifo,
    InterferometerNumber IFONumber
    )

{
  switch( IFONumber )
  {
    case LAL_IFO_G1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "G1");
      break;

    case LAL_IFO_H1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "H1");
      break;

    case LAL_IFO_H2:
      snprintf( ifo, LIGOMETA_IFO_MAX, "H2");
      break;

    case LAL_IFO_L1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "L1");
      break;

    case LAL_IFO_T1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "T1");
      break;

    case LAL_IFO_V1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "V1");
      break;

    default:
      /* Invalid Detector Site */
      snprintf( ifo, LIGOMETA_IFO_MAX, " ");
  }
}



void
XLALReturnDetector(
    LALDetector           *det,
    InterferometerNumber   IFONumber
    )

{
  switch( IFONumber )
  {
    case LAL_IFO_G1:
      *det = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      break;

    case LAL_IFO_H1:
      *det = lalCachedDetectors[LALDetectorIndexLHODIFF];
      break;

    case LAL_IFO_H2:
      *det = lalCachedDetectors[LALDetectorIndexLHODIFF];
      break;

    case LAL_IFO_L1:
      *det = lalCachedDetectors[LALDetectorIndexLLODIFF];
      break;

    case LAL_IFO_T1:
      *det = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      break;

    case LAL_IFO_V1:
      *det = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
      break;

    default:
      /* Invalid Detector Site */
      memset(det, 0, sizeof(LALDetector) );
  }
}
