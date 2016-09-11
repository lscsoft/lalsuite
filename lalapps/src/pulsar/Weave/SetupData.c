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

#include "Weave.h"

#include <lal/LALStdio.h>
#include <lal/LALString.h>
#include <lal/UserInput.h>
#include <lal/GSLHelpers.h>
#include <lal/LALInitBarycenter.h>

static int setup_LALSeg_fits_table_init( FITSFile *file );
static int setup_PosVelAcc_fits_table_init( FITSFile *file );

///
/// Initialise a FITS table for writing/reading a table of PosVelAcc entries
///
int setup_PosVelAcc_fits_table_init(
  FITSFile *file
  )
{
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_FITS_TABLE_COLUMN_BEGIN( PosVelAcc );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, gps ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY( file, REAL8, pos ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY( file, REAL8, vel ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY( file, REAL8, acc ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;
}

///
/// Initialise a FITS table for writing/reading a table of LALSeg entries
///
int setup_LALSeg_fits_table_init(
  FITSFile *file
  )
{
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_FITS_TABLE_COLUMN_BEGIN( LALSeg );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, start ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, end ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, INT4, id ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;
}

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
  {

    // Write Earth ephemerides
    {
      XLAL_CHECK( XLALFITSTableOpenWrite( file, "earth_ephem", "Earth ephemeris" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( setup_PosVelAcc_fits_table_init( file ) == XLAL_SUCCESS, XLAL_EFUNC );
      for ( INT4 i = 0; i < setup->ephemerides->nentriesE; ++i ) {
        XLAL_CHECK( XLALFITSTableWriteRow( file, &setup->ephemerides->ephemE[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
      XLAL_CHECK( XLALFITSHeaderWriteString( file, "filename", setup->ephemerides->ephiles.earthEphemeris, "ephemeris filename" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALFITSHeaderWriteINT4( file, "etype", setup->ephemerides->etype, "ephemeris type" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, "dttable", setup->ephemerides->dtEtable, "spacing in seconds" ) == XLAL_SUCCESS, XLAL_EFUNC );

    }

    // Write Sun ephemerides
    {
      XLAL_CHECK( XLALFITSTableOpenWrite( file, "sun_ephem", "Sun ephemeris" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( setup_PosVelAcc_fits_table_init( file ) == XLAL_SUCCESS, XLAL_EFUNC );
      for ( INT4 i = 0; i < setup->ephemerides->nentriesS; ++i ) {
        XLAL_CHECK( XLALFITSTableWriteRow( file, &setup->ephemerides->ephemS[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
      XLAL_CHECK( XLALFITSHeaderWriteString( file, "filename", setup->ephemerides->ephiles.sunEphemeris, "ephemeris filename" ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALFITSHeaderWriteINT4( file, "etype", setup->ephemerides->etype, "ephemeris type" )  == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALFITSHeaderWriteREAL8( file, "dttable", setup->ephemerides->dtStable, "spacing in seconds" )  == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Write segment list
  {
    XLAL_CHECK( XLALFITSTableOpenWrite( file, "segments", "segment list" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( setup_LALSeg_fits_table_init( file ) == XLAL_SUCCESS, XLAL_EFUNC );
    for ( size_t i = 0; i < setup->segments->length; ++i ) {
      XLAL_CHECK( XLALFITSTableWriteRow( file, &setup->segments->segs[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Write supersky metrics
  {

    // Write coherent metrics
    {
      size_t dims[3] = {
        setup->metrics->coh_rssky_metric[0]->size1,
        setup->metrics->coh_rssky_metric[0]->size2,
        setup->metrics->num_segments
      };
      XLAL_CHECK( XLALFITSArrayOpenWrite( file, "coh_rssky_metric", 3, dims, "coherent supersky metrics" ) == XLAL_SUCCESS, XLAL_EFUNC );
      size_t idx[3];
      for ( idx[2] = 0; idx[2] < dims[2]; ++idx[2] ) {
        XLAL_CHECK( XLALFITSArrayWriteGSLMatrix( file, idx, setup->metrics->coh_rssky_metric[idx[2]] ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }

    // Write coherent metric transform data
    {
      size_t dims[3] = {
        setup->metrics->coh_rssky_transf[0]->size1,
        setup->metrics->coh_rssky_transf[0]->size2,
        setup->metrics->num_segments
      };
      XLAL_CHECK( XLALFITSArrayOpenWrite( file, "coh_rssky_transf", 3, dims, "coherent supersky metric transform data" ) == XLAL_SUCCESS, XLAL_EFUNC );
      size_t idx[3];
      for ( idx[2] = 0; idx[2] < dims[2]; ++idx[2] ) {
        XLAL_CHECK( XLALFITSArrayWriteGSLMatrix( file, idx, setup->metrics->coh_rssky_transf[idx[2]] ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }

    // Write semicoherent metric
    XLAL_CHECK( XLALFITSArrayOpenWrite2( file, "semi_rssky_metric", setup->metrics->semi_rssky_metric->size1, setup->metrics->semi_rssky_metric->size2, "semicoherent supersky metric" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSArrayWriteGSLMatrix( file, NULL, setup->metrics->semi_rssky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Write semicoherent metric transform data
    XLAL_CHECK( XLALFITSArrayOpenWrite2( file, "semi_rssky_transf", setup->metrics->semi_rssky_transf->size1, setup->metrics->semi_rssky_transf->size2, "semicoherent supersky metric transform data" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSArrayWriteGSLMatrix( file, NULL, setup->metrics->semi_rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );

  }

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

  // Read ephemerides
  {

    // Allocate memory
    setup->ephemerides = XLALCalloc( 1, sizeof( *setup->ephemerides ) );
    XLAL_CHECK( setup->ephemerides != NULL, XLAL_ENOMEM );

    // Read Earth ephemerides
    {
      UINT8 nrows = 0;
      XLAL_CHECK( XLALFITSTableOpenRead( file, "earth_ephem", &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( setup_PosVelAcc_fits_table_init( file ) == XLAL_SUCCESS, XLAL_EFUNC );
      setup->ephemerides->nentriesE = nrows;
      setup->ephemerides->ephemE = XLALCalloc( setup->ephemerides->nentriesE, sizeof( setup->ephemerides->ephemE[0] ) );
      XLAL_CHECK( setup->ephemerides->ephemE != NULL, XLAL_EINVAL );
      for ( INT4 i = 0; i < setup->ephemerides->nentriesE; ++i ) {
        XLAL_CHECK( XLALFITSTableReadRow( file, &setup->ephemerides->ephemE[i], NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
      XLAL_CHECK( XLALFITSHeaderReadString( file, "filename", &setup->ephemerides->ephiles.earthEphemeris ) == XLAL_SUCCESS, XLAL_EFUNC );
      INT4 etype = 0;
      XLAL_CHECK( XLALFITSHeaderReadINT4( file, "etype", &etype ) == XLAL_SUCCESS, XLAL_EFUNC );
      setup->ephemerides->etype = etype;
      XLAL_CHECK( XLALFITSHeaderReadREAL8( file, "dttable", &setup->ephemerides->dtEtable ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

    // Read Sun ephemerides
    {
      UINT8 nrows = 0;
      XLAL_CHECK( XLALFITSTableOpenRead( file, "sun_ephem", &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( setup_PosVelAcc_fits_table_init( file ) == XLAL_SUCCESS, XLAL_EFUNC );
      setup->ephemerides->nentriesS = nrows;
      setup->ephemerides->ephemS = XLALCalloc( setup->ephemerides->nentriesS, sizeof( setup->ephemerides->ephemS[0] ) );
      XLAL_CHECK( setup->ephemerides->ephemS != NULL, XLAL_EINVAL );
      for ( INT4 i = 0; i < setup->ephemerides->nentriesS; ++i ) {
        XLAL_CHECK( XLALFITSTableReadRow( file, &setup->ephemerides->ephemS[i], NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
      XLAL_CHECK( XLALFITSHeaderReadString( file, "filename", &setup->ephemerides->ephiles.sunEphemeris ) == XLAL_SUCCESS, XLAL_EFUNC );
      INT4 etype = 0;
      XLAL_CHECK( XLALFITSHeaderReadINT4( file, "etype", &etype ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( setup->ephemerides->etype == ( EphemerisType ) etype, XLAL_EIO );
      XLAL_CHECK( XLALFITSHeaderReadREAL8( file, "dttable", &setup->ephemerides->dtStable ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  }

  // Read segment list
  {
    UINT8 nrows = 0;
    XLAL_CHECK( XLALFITSTableOpenRead( file, "segments", &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( setup_LALSeg_fits_table_init( file ) == XLAL_SUCCESS, XLAL_EFUNC );
    setup->segments = XLALSegListCreate();
    XLAL_CHECK( setup->segments != NULL, XLAL_EFUNC );
    while ( nrows > 0 ) {
      LALSeg XLAL_INIT_DECL( seg );
      XLAL_CHECK( XLALFITSTableReadRow( file, &seg, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALSegListAppend( setup->segments, &seg ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Read supersky metrics
  {

    // Allocate memory
    setup->metrics = XLALCalloc( 1, sizeof( *setup->metrics ) );
    XLAL_CHECK( setup->metrics != NULL, XLAL_ENOMEM );

    // Read coherent metrics
    {
      size_t ndim, dims[FFIO_MAX];
      XLAL_CHECK( XLALFITSArrayOpenRead( file, "coh_rssky_metric", &ndim, dims ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( ndim == 3, XLAL_EIO );
      setup->metrics->num_segments = dims[2];
      setup->metrics->coh_rssky_metric = XLALCalloc( setup->metrics->num_segments, sizeof( setup->metrics->coh_rssky_metric[0] ) );
      XLAL_CHECK( setup->metrics->coh_rssky_metric != NULL, XLAL_ENOMEM );
      size_t idx[3];
      for ( idx[2] = 0; idx[2] < dims[2]; ++idx[2] ) {
        XLAL_CHECK( XLALFITSArrayReadGSLMatrix( file, idx, &setup->metrics->coh_rssky_metric[idx[2]] ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }

    // Read coherent metric transform data
    {
      size_t ndim, dims[FFIO_MAX];
      XLAL_CHECK( XLALFITSArrayOpenRead( file, "coh_rssky_transf", &ndim, dims ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( ndim == 3, XLAL_EIO );
      XLAL_CHECK( dims[2] == setup->metrics->num_segments, XLAL_EIO );
      setup->metrics->coh_rssky_transf = XLALCalloc( setup->metrics->num_segments, sizeof( setup->metrics->coh_rssky_transf[0] ) );
      XLAL_CHECK( setup->metrics->coh_rssky_transf != NULL, XLAL_ENOMEM );
      size_t idx[3];
      for ( idx[2] = 0; idx[2] < dims[2]; ++idx[2] ) {
        XLAL_CHECK( XLALFITSArrayReadGSLMatrix( file, idx, &setup->metrics->coh_rssky_transf[idx[2]] ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }

    // Read semicoherent metric
    XLAL_CHECK( XLALFITSArrayOpenRead2( file, "semi_rssky_metric", NULL, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSArrayReadGSLMatrix( file, NULL, &setup->metrics->semi_rssky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Read semicoherent metric transform data
    XLAL_CHECK( XLALFITSArrayOpenRead2( file, "semi_rssky_transf", NULL, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSArrayReadGSLMatrix( file, NULL, &setup->metrics->semi_rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );

  }

  // Check output
  XLAL_CHECK( setup->detectors != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup->segments != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup->metrics != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup->ephemerides != NULL, XLAL_EFAULT );

  return XLAL_SUCCESS;

}
