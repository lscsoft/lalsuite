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

#ifndef _FITSPULSARIO_H
#define _FITSPULSARIO_H

#include <stdint.h>
#include <gsl/gsl_matrix.h>
#include <lal/LALStdlib.h>
#include <lal/FITSFileIO.h>
#include <lal/Segments.h>
#include <lal/LALBarycenter.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// \defgroup FITSPulsarIO_h Header FITSPulsarIO.h
/// \ingroup lalpulsar_general
/// \author Karl Wette
/// \brief Routines for reading/writing pulsar-related data to/from FITS files
///

/// @{

///
/// Write a segment list to a FITS file
///
int XLALFITSWriteSegmentList(
  FITSFile *file,                       ///< [in] FITS file pointer
  const CHAR *name,                     ///< [in] Name of FITS table to write segment list to
  const LALSegList *segments,           ///< [in] Segment list
  const CHAR *comment                   ///< [in] Comment for FITS table
  );

///
/// Read a segment list from a FITS file
///
#ifdef SWIG // SWIG interface directives
SWIGLAL( INOUT_STRUCTS( LALSegList **, segments ) );
#endif
int XLALFITSReadSegmentList(
  FITSFile *file,                       ///< [in] FITS file pointer
  const CHAR *name,                     ///< [in] Name of FITS table to read segment list from
  LALSegList **segments                 ///< [out] Segment list
  );

///
/// Write ephemeris data to a FITS file
///
int XLALFITSWriteEphemerisData(
  FITSFile *file,                       ///< [in] FITS file pointer
  const EphemerisData *ephemerides      ///< [in] Ephemeris data
  );

///
/// Read ephemeris data from a FITS file
///
#ifdef SWIG // SWIG interface directives
SWIGLAL( INOUT_STRUCTS( EphemerisData **, ephemerides ) );
#endif
int XLALFITSReadEphemerisData(
  FITSFile *file,                       ///< [in] FITS file pointer
  EphemerisData **ephemerides           ///< [out] Ephemeris data
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _FITSPULSARIO_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
