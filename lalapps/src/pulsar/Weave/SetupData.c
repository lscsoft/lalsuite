//
// Copyright (C) 2016, 2017 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "SetupData.h"

#include <lal/LALInitBarycenter.h>
#include <lal/FITSPulsarIO.h>

///
/// \name Internal routines
///
/// @{

/// @}

///
/// Free contents of setup data
///
void XLALWeaveSetupDataClear(
  WeaveSetupData *setup
  )
{
  if ( setup != NULL ) {
    XLALDestroyStringVector( setup->detectors );
    XLALSegListFree( setup->segments );
    XLALDestroySuperskyMetrics( setup->metrics );
    XLALDestroyEphemerisData( setup->ephemerides );
  }
}

///
/// Write setup data to a FITS file
///
int XLALWeaveSetupDataWrite(
  FITSFile *file,
  const WeaveSetupData *setup
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup->detectors != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup->segments != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup->metrics != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup->ephemerides != NULL, XLAL_EFAULT );

  // Write reference time
  XLAL_CHECK( XLALFITSHeaderWriteGPSTime( file, "date-obs", &setup->ref_time, "reference time" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write list of detector names
  XLAL_CHECK( XLALFITSHeaderWriteStringVector( file, "detect", setup->detectors, "setup detectors" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write ephemerides
  XLAL_CHECK( XLALFITSWriteEphemerisData( file, setup->ephemerides ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write segment list
  XLAL_CHECK( XLALFITSWriteSegmentList( file, "segments", setup->segments, "segment list" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write supersky metrics
  XLAL_CHECK( XLALFITSWriteSuperskyMetrics( file, setup->metrics ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Read setup data from a FITS file
///
int XLALWeaveSetupDataRead(
  FITSFile *file,
  WeaveSetupData *setup
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup != NULL, XLAL_EFAULT );

  // Erase memory
  XLAL_INIT_MEM( *setup );

  // Set coordinate transforms between physical and lattice coordinates
  setup->phys_to_latt = ( WeavePhysicalToLattice ) XLALConvertPhysicalToSuperskyPoint;
  setup->latt_to_phys = ( WeaveLatticeToPhysical ) XLALConvertSuperskyToPhysicalPoint;

  // Read reference time
  XLAL_CHECK( XLALFITSHeaderReadGPSTime( file, "date-obs", &setup->ref_time ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Read list of detector names
  XLAL_CHECK( XLALFITSHeaderReadStringVector( file, "detect", &setup->detectors ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( setup->detectors != NULL, XLAL_EFAULT );

  // Read ephemerides
  XLAL_CHECK( XLALFITSReadEphemerisData( file, &setup->ephemerides ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( setup->ephemerides != NULL, XLAL_EFAULT );

  // Read segment list
  XLAL_CHECK( XLALFITSReadSegmentList( file, "segments", &setup->segments ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( setup->ephemerides != NULL, XLAL_EFAULT );

  // Read supersky metrics
  XLAL_CHECK( XLALFITSReadSuperskyMetrics( file, &setup->metrics ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( setup->metrics != NULL, XLAL_EFAULT );

  return XLAL_SUCCESS;

}

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
