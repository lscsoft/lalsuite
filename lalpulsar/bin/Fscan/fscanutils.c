/*
*  Copyright (C) 2023 Evan Goetz
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
* \file
* \ingroup lalpulsar_bin_fscan
*/

#include <lal/SFTfileIO.h>
#include <lal/LALStdio.h>

#include "fscanutils.h"

/* Extract a single SFT from an SFTCatalog: the SFT indicated by GPS start time with band f_min to f_max*/
SFTVector *extract_one_sft( const SFTCatalog *full_catalog, const LIGOTimeGPS starttime, const REAL8 f_min, const REAL8 f_max )
{
  // Initialize an SFTCatalog
  SFTCatalog XLAL_INIT_DECL( catalogSlice );

  //Set start time
  //Set end time just 0.01 seconds after the start time. This is sufficiently small to get just one SFT
  LIGOTimeGPS thisSFTendtime = starttime;
  XLAL_CHECK_NULL( XLALGPSAdd( &thisSFTendtime, 0.01 ) != NULL, XLAL_EFUNC );

  // Get the catalog of the single SFT from the full catalog
  XLAL_CHECK_NULL( XLALSFTCatalogTimeslice( &catalogSlice, full_catalog, &starttime, &thisSFTendtime ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Check that we got one SFT
  XLAL_CHECK_NULL( catalogSlice.length == 1, XLAL_EFUNC, "Found no unique SFT starting at %d", starttime.gpsSeconds );

  //Extract the SFT
  SFTVector *sft_vect = NULL;
  XLAL_CHECK_NULL( ( sft_vect = XLALLoadSFTs( &catalogSlice, f_min, f_max ) ) != NULL, XLAL_EFUNC );

  //Check we got only zero or one SFT; no more
  XLAL_CHECK_NULL( sft_vect->length == 1, XLAL_EBADLEN, "SFT in catalog could not cover [%.2f, %.2f) Hz", f_min, f_max );

  return sft_vect;
}
