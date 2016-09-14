//
// Copyright (C) 2016 Karl Wette
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

#ifndef _SETUP_DATA_H
#define _SETUP_DATA_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"

#include <lal/LALBarycenter.h>
#include <lal/Segments.h>
#include <lal/SuperskyMetrics.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Setup data which is computed only once for a given search setup
///
typedef struct {
  /// Physical to lattice coordinate transform
  WeavePhysicalToLattice phys_to_latt;
  /// Lattice to physical coordinate transform
  WeaveLatticeToPhysical latt_to_phys;
  /// Reference time at which search is conducted
  LIGOTimeGPS ref_time;
  /// List of detector names for which metrics were computed
  LALStringVector *detectors;
  /// Segment list for which metrics were computed
  LALSegList *segments;
  /// Reduced supersky parameter-space metrics
  SuperskyMetrics *metrics;
  /// Ephemeris data over time-span of segments
  EphemerisData *ephemerides;
} WeaveSetupData;

///
/// \name Routines which handle the setup data
///
/// @{

void XLALWeaveSetupDataClear(
  WeaveSetupData *setup
  );
int XLALWeaveSetupDataWrite(
  FITSFile *file,
  const WeaveSetupData *setup
  );
int XLALWeaveSetupDataRead(
  FITSFile *file,
  WeaveSetupData *setup
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _SETUP_DATA_H
