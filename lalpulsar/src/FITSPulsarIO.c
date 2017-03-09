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

#include <lal/FITSPulsarIO.h>

///
/// Initialise a FITS table for writing/reading a table of LALSeg entries
///
static int fits_table_init_LALSeg(
  FITSFile *file
  )
{
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_FITS_TABLE_COLUMN_BEGIN( LALSeg );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, start ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, end ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, INT4, id, "numSFTs" ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;
}

int XLALFITSWriteSegmentList(
  FITSFile *file,
  const CHAR *name,
  const LALSegList *segments,
  const CHAR *comment
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( name != NULL, XLAL_EFAULT );
  XLAL_CHECK( segments != NULL, XLAL_EFAULT );
  XLAL_CHECK( comment != NULL, XLAL_EFAULT );

  // Write segment list to a FITS table
  XLAL_CHECK( XLALFITSTableOpenWrite( file, name, comment ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fits_table_init_LALSeg( file ) == XLAL_SUCCESS, XLAL_EFUNC );
  for ( size_t i = 0; i < segments->length; ++i ) {
    XLAL_CHECK( XLALFITSTableWriteRow( file, &segments->segs[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}

int XLALFITSReadSegmentList(
  FITSFile *file,
  const CHAR *name,
  LALSegList **segments
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( name != NULL, XLAL_EFAULT );
  XLAL_CHECK( segments != NULL && *segments == NULL, XLAL_EFAULT );

  // Read segment list from a FITS table
  UINT8 nrows = 0;
  XLAL_CHECK( XLALFITSTableOpenRead( file, name, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fits_table_init_LALSeg( file ) == XLAL_SUCCESS, XLAL_EFUNC );
  *segments = XLALSegListCreate();
  XLAL_CHECK( *segments != NULL, XLAL_EFUNC );
  while ( nrows > 0 ) {
    LALSeg XLAL_INIT_DECL( seg );
    XLAL_CHECK( XLALFITSTableReadRow( file, &seg, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALSegListAppend( *segments, &seg ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
