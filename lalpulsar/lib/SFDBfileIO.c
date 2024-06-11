/*
 * Copyright (C) 2020 Rodrigo Tenorio
 * Copyright (C) 2020 Pep Covas
 * Copyright (C) 2019, 2020 David Keitel
 * Copyright (C) 2014, 2022 Karl Wette
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

/*---------- includes ----------*/

#include "SFTinternal.h"

/*---------- constants ----------*/

/* detector numbers as defined in Rome SFDBs
 * 0 is Nautilus but we won't support that
 */
typedef enum tagSFDBDetectors {
  SFDB_DET_FIRST = 1,
  SFDB_DET_V1 = 1,
  SFDB_DET_H1 = 2,
  SFDB_DET_L1 = 3,
  SFDB_DET_LAST
} SFDBDetectors;

const char *const SFDB_detector_names[] = {
  [SFDB_DET_V1] = "V1",
  [SFDB_DET_H1] = "H1",
  [SFDB_DET_L1] = "L1"
};

/*---------- internal types ----------*/

/* header contents of SFDBs, many fields unused
 */
typedef struct tagSFDBHeader {
  INT4 det;  // Index of the detector (0=V1, 1=H1, 2=L1)
  INT4 gps_sec;  // GPS time of this SFDB
  INT4 gps_nsec;
  REAL8 tbase;  // Coherent time of the SFDB
  INT4 firstfrind;
  INT4 nsamples;  // Number of frequency bins
  INT4 red;  // Reduction factor
  INT4 typ;
  REAL4 n_flag;
  REAL4 einstein;  // Normalization factor to get back to 10^{-20} units
  REAL8 mjdtime;
  INT4 nfft;  // Number of FFTs in a SFDB
  INT4 wink;
  REAL4 normd;  // Normalization factor
  REAL4 normw;  // Window normalization factor
  REAL8 frinit;  // Initial frequency
  REAL8 tsamplu;  // Sampling time
  REAL8 deltanu;  // Frequency resolution
  REAL8 vx_eq;  // Velocity vector of the detector at this GPS time
  REAL8 vy_eq;
  REAL8 vz_eq;
  REAL8 px_eq;  // Position vector of the detector at this GPS time
  REAL8 py_eq;
  REAL8 pz_eq;
  INT4 n_zeroes;
  REAL8 sat_howmany;
  INT4 lavesp;
} SFDBHeader;

/*---------- internal prototypes ----------*/

static int read_SFDB_header_from_fp( FILE *fp, SFDBHeader *header );
static BOOLEAN CheckIfSFDBInScienceMode( SFDBHeader *SFDBHeader, LALStringVector *detectors, MultiLIGOTimeGPSVector *startingTS, MultiLIGOTimeGPSVector *endingTS );

/*========== function definitions ==========*/

/// \addtogroup SFTfileIO_h
/// @{

/**
 * Return a MultiSFTVector struct from an input set of SFDBs, possibly from more than one detector.
 *
 * An SFDB (Short Fourier DataBase) is the frequency-domain data format created by the Rome group. It has a time-domain cleaning, which is described in \cite Astone2014.
 *
 * In order to only use SFDBs within science segments, it is possible to input files for each detector which have the science segments. Two files for each detector are needed, one with the starting timestamps and the other with the ending timestamps. The format for these files is one timestamp per line. If not needed, the timestamp inputs can be NULL.
 *
 * The returned SFTs in the standard LALSuite format are sorted by increasing GPS-epochs!
 */
MultiSFTVector *
XLALReadSFDB(
  REAL8 f_min,                       // Minimum frequency to be read
  REAL8 f_max,                       // Maximum frequency to be read
  const CHAR *file_pattern,          // String of SFDB files (possibly from more than one detector, separated by a ;)
  const CHAR *timeStampsStarting,    // File(s) containing the starting timestamps of science segments (possibly from more than one detector, separated by a ;)
  const CHAR *timeStampsFinishing    // File(s) containing the finishing timestamps of science segments (possibly from more than one detector, separated by a ;)
)
{

  /* ------ Basic Setup ------ */
  LALStringVector *fnamesSFDB;
  XLAL_CHECK_NULL( ( fnamesSFDB = XLALFindFiles( file_pattern ) ) != NULL, XLAL_EFUNC, "Failed to find filelist matching pattern '%s'.\n\n", file_pattern );
  UINT4 numSFDBFiles = fnamesSFDB->length;

  BOOLEAN flag_timestamps; /* This flag is FALSE if no timestamps are used, and TRUE if timestamps to only load SFDBs within science segments are used */
  LALStringVector *fnamesStartTS = NULL; /* There is one pair (startTS, endTS) for each IFO: #IFOs == startingTS->length == endingTS->length */
  LALStringVector *fnamesEndTS = NULL;
  MultiLIGOTimeGPSVector *startingTS = NULL;
  MultiLIGOTimeGPSVector *endingTS = NULL;

  LALStringVector *detectors = NULL; /* Only used if flag_timestamps == True */

  INT4 Tsft = 0;

  if ( timeStampsStarting && timeStampsFinishing ) {
    flag_timestamps = TRUE;
  } else if ( timeStampsStarting && timeStampsFinishing == NULL ) {
    XLAL_ERROR_NULL( XLAL_EFUNC, "Must give two files with initial and finishing timestamps, missing finishing timestamps\n" );
  } else if ( timeStampsStarting == NULL && timeStampsFinishing ) {
    XLAL_ERROR_NULL( XLAL_EFUNC, "Must give two files with initial and finishing timestamps, missing starting timestamps\n" );
  } else {
    flag_timestamps = FALSE;
  }

  if ( flag_timestamps ) {
    XLAL_CHECK_NULL( ( fnamesStartTS = XLALFindFiles( timeStampsStarting ) ) != NULL, XLAL_EFUNC, "Failed to find filelist matching pattern '%s'.\n\n", timeStampsStarting );
    XLAL_CHECK_NULL( ( fnamesEndTS = XLALFindFiles( timeStampsFinishing ) ) != NULL, XLAL_EFUNC, "Failed to find filelist matching pattern '%s'.\n\n", timeStampsFinishing );
    UINT4 numTSFiles = fnamesStartTS->length;
    XLAL_CHECK_NULL( numTSFiles == fnamesEndTS->length, XLAL_EINVAL );
    XLAL_CHECK_NULL( numTSFiles > 0, XLAL_EINVAL );

    startingTS = XLALReadMultiTimestampsFiles( fnamesStartTS );
    endingTS = XLALReadMultiTimestampsFiles( fnamesEndTS );

    for ( UINT4 X = 0; X < numTSFiles; ++X ) {

      /* Retrieve IFO names following the SFDB convention */
      CHAR *filenameSt = fnamesStartTS->data[X];
      BOOLEAN detectorWasFound = 0;
      for ( UINT4 Y = SFDB_DET_FIRST; Y < SFDB_DET_LAST; Y++ ) {
        if ( strstr( filenameSt, SFDB_detector_names[Y] ) ) {
          XLAL_CHECK_NULL(
            ( detectors = XLALAppendString2Vector( detectors,
                          SFDB_detector_names[Y] )
            ) != NULL,
            XLAL_EFUNC
          );
          detectorWasFound = 1;
        }
      }
      XLAL_CHECK_NULL(
        detectorWasFound,
        XLAL_EINVAL,
        "No matching IFO name was found for time stamp file %s",
        filenameSt
      );

      /* Time stamp stanity check */
      UINT4 numStartStamps = startingTS->data[X]->length;
      UINT4 numEndStamps = endingTS->data[X]->length;
      XLAL_CHECK_NULL(
        numStartStamps == numEndStamps,
        XLAL_EINVAL,
        "Got %u starting and %u finishing timestamps at %s, lengths must be equal.",
        numStartStamps,
        numEndStamps,
        filenameSt
      );

    }/* -- IFO < numTSfiles -- */

  }/* -- if flag_timestamps -- */

  /* ------ First Step: Count SFDB files (and IFOs if required by time stamps) to allocate enough memory ------ */

  /* Y indices loop over the SFDB detector ordering convention (regardless of whether data is actually present or not) */
  UINT4 numSFTsY[SFDB_DET_LAST];
  XLAL_INIT_MEM( numSFTsY );
  for ( UINT4 i = 0; i < numSFDBFiles; i++ ) {
    const CHAR *filename = fnamesSFDB->data[i];
    FILE  *fpPar = NULL;
    XLAL_CHECK_NULL( ( fpPar = fopen( filename, "r" ) ) != NULL, XLAL_EIO, "Failed to open SFDB file '%s' for reading.", filename );
    setvbuf( fpPar, ( CHAR * )NULL, _IOLBF, 0 );

    REAL8 count; /* Index of an SFDBs in the file (a SFDB file can have more than one SFDB) */
    while ( fread( &count, sizeof( REAL8 ), 1, fpPar ) == 1 ) {

      SFDBHeader header;
      XLAL_CHECK_NULL( read_SFDB_header_from_fp( fpPar, &header ) == 0, XLAL_EIO, "Failed to parse SFDB header." );

      /* If time stamps are given, only count science mode SFDBs */
      numSFTsY[header.det] += ( flag_timestamps ) ?
                              CheckIfSFDBInScienceMode( &header, detectors, startingTS, endingTS ) :
                              1;

      Tsft = header.tbase;

      /* Skip number of bytes corresponding to the actual data content */
      INT4 lavespOrRed = 0;
      UINT4 lsps = 0;
      if ( header.lavesp > 0 ) {
        lavespOrRed = header.lavesp;
        lsps = header.lavesp;
      } else {
        lavespOrRed = header.red;
        lsps = header.nsamples / header.red;
      }
      XLAL_CHECK_NULL( fseek( fpPar, lavespOrRed * sizeof( REAL4 ), SEEK_CUR ) == 0, XLAL_EIO );
      XLAL_CHECK_NULL( fseek( fpPar, ( lsps + 2 * header.nsamples ) * sizeof( REAL4 ), SEEK_CUR ) == 0, XLAL_EIO );

    }/* -- while fread(...) == 1 -- */

    fclose( fpPar );

  }/* -- while i < numSFDBFiles -- */

  /* ------ Second Step: Reformat retrieved information  ------ */

  UINT4 numSFTsTotal = 0;
  for ( UINT4 Y = SFDB_DET_FIRST; Y < SFDB_DET_LAST; Y++ ) {
    numSFTsTotal += numSFTsY[Y];
  }
  XLAL_CHECK_NULL( numSFTsTotal > 0, XLAL_EINVAL, "No SFTs found for any detector." );

  /* Prepare a lookup table containing detectors present in the SFDBs */
  INT4 detectorLookupYtoX[SFDB_DET_LAST];
  CHAR detectorNames[SFDB_DET_LAST][3];
  XLAL_INIT_MEM( detectorNames );
  for ( UINT4 Y = 0; Y < SFDB_DET_LAST; Y++ ) {
    detectorLookupYtoX[Y] = -1;
    strncpy( detectorNames[Y], "XX", 3 );
  }

  /* Count IFOs and SFTs per IFO */
  /* To do so, loop over everything, retrieve present */
  /* information and re-format it into an array. */
  UINT4 numIFOs = 0;
  UINT4Vector *auxNumSFTsX = NULL;
  auxNumSFTsX = XLALCreateUINT4Vector( SFDB_DET_LAST );
  for ( UINT4 Y = SFDB_DET_FIRST; Y < SFDB_DET_LAST; Y++ ) {
    if ( numSFTsY[Y] > 0 ) {
      /* numIFOs is used here as an internal loop variable */
      /* of actually present detectors (equivalent to X later) */
      /* and will be equal to the total number at the end of the loop */
      strncpy( detectorNames[numIFOs], SFDB_detector_names[Y], 3 );
      auxNumSFTsX->data[numIFOs] = numSFTsY[Y];
      detectorLookupYtoX[Y] = numIFOs;
      numIFOs += 1;
    }
  }

  /* X indices loop over the actually present detectors */
  UINT4Vector *numSFTsX = NULL;
  numSFTsX =  XLALCreateUINT4Vector( numIFOs );
  XLALPrintInfo( "Number of SFTs we'll load from the SFDBs:\n" );
  for ( UINT4 X = 0; X < numIFOs; X++ ) {
    numSFTsX->data[X] = auxNumSFTsX->data[X];
    XLALPrintInfo( "%s: %d\n", detectorNames[X], numSFTsX->data[X] );
  }

  /* Allocate memory for the SFT structure */
  /* Calling XLALFindCoveringSFTBins() guarantees we use the same bandwidth */
  /* conventions as e.g. XLALCWMakeFakeMultiData() */
  UINT4 firstBinExt, numBinsExt;
  XLAL_CHECK_NULL( XLALFindCoveringSFTBins( &firstBinExt, &numBinsExt, f_min, f_max - f_min, Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );
  MultiSFTVector *outputSFTs = NULL;
  outputSFTs = XLALCreateMultiSFTVector( numBinsExt, numSFTsX );

  /* Clean up auxiliar variables */
  XLALDestroyUINT4Vector( numSFTsX );
  XLALDestroyUINT4Vector( auxNumSFTsX );

  /* ------ Third Step: Fill up SFTs using SFDB data ------ */

  UINT4 numSFTsLoadedInX[numIFOs];
  XLAL_INIT_MEM( numSFTsLoadedInX );
  for ( UINT4 i = 0; i < numSFDBFiles; i++ ) {
    const CHAR *filename = fnamesSFDB->data[i];
    FILE  *fpPar = NULL;
    XLAL_CHECK_NULL( ( fpPar = fopen( filename, "r" ) ) != NULL, XLAL_EIO, "Failed to open SFDB file '%s' for reading.", filename );
    setvbuf( fpPar, ( CHAR * )NULL, _IOLBF, 0 );

    REAL8 count;
    while ( fread( &count, sizeof( REAL8 ), 1, fpPar ) == 1 ) { /* Read SFDBs one by one */

      SFDBHeader header;
      XLAL_CHECK_NULL( read_SFDB_header_from_fp( fpPar, &header ) == 0, XLAL_EIO, "Failed to parse SFDB header." );


      /* Read the actual binary data */
      /* This is needed regardless of this data falling inside science mode */
      /* since fread is what moves the pointer forwards! */
      INT4 lavespOrRed = 0;
      UINT4 lsps = 0;
      REAL4 *buffer1 = NULL,  *buffer2 = NULL, *buffer3 = NULL;
      if ( header.lavesp > 0 ) {
        lavespOrRed = header.lavesp;
        lsps = header.lavesp;
      } else {
        lavespOrRed = header.red;
        lsps = header.nsamples / header.red;
      }

      XLAL_CHECK_NULL( ( buffer1 = XLALCalloc( lavespOrRed, sizeof( REAL4 ) ) ) != NULL, XLAL_ENOMEM );
      XLAL_CHECK_NULL( fread( buffer1, lavespOrRed * sizeof( REAL4 ), 1, fpPar ) == 1, XLAL_EIO );

      XLAL_CHECK_NULL( ( buffer2 = XLALCalloc( lsps, sizeof( REAL4 ) ) ) != NULL, XLAL_ENOMEM );
      XLAL_CHECK_NULL( fread( buffer2, lsps * sizeof( REAL4 ), 1, fpPar ) == 1, XLAL_EIO );

      XLAL_CHECK_NULL( ( buffer3 = XLALCalloc( 2 * header.nsamples, sizeof( REAL4 ) ) ) != NULL, XLAL_ENOMEM );
      XLAL_CHECK_NULL( fread( buffer3, 2 * header.nsamples * sizeof( REAL4 ), 1, fpPar ) == 1, XLAL_EIO );

      if ( flag_timestamps ? CheckIfSFDBInScienceMode( &header, detectors, startingTS, endingTS ) : 1 ) {
        XLAL_CHECK_NULL( detectorLookupYtoX[header.det] >= 0, XLAL_EDOM, "Cannot match detector %d, as read from file, with first run.", header.det );
        UINT4 X = detectorLookupYtoX[header.det];
        numSFTsLoadedInX[X] += 1;

        /* Fill up this SFT */
        SFTtype *thisSFT = NULL;
        thisSFT = outputSFTs->data[X]->data + numSFTsLoadedInX[X] - 1;

        XLALGPSSetREAL8( &( thisSFT->epoch ), header.gps_sec );
        strncpy( thisSFT->name, detectorNames[X], 3 );
        thisSFT->f0 = f_min;
        thisSFT->deltaF = header.deltanu;

        /* Copy to the SFT structure the complex values in the frequency bins (multiplied by some normalization factors in order to agree with the SFT specification) */
        UINT4 thisSFTLength = thisSFT->data->length;
        COMPLEX8 *thisSFTData;
        thisSFTData = thisSFT->data->data;
        for ( UINT4 thisBin = firstBinExt; thisBin < ( thisSFTLength + firstBinExt ); thisBin++ ) {
          *thisSFTData = crectf( buffer3[2 * thisBin], buffer3[2 * thisBin + 1] ) * header.einstein * header.tsamplu * header.normw;
          ++thisSFTData;
        }
      }

      XLALFree( buffer1 );
      XLALFree( buffer2 );
      XLALFree( buffer3 );

    }/* while fread(...) == 1 */

    fclose( fpPar );
  }/* -- i < numSFDBFiles -- */

  XLALDestroyStringVector( fnamesSFDB );
  XLALDestroyStringVector( fnamesStartTS );
  XLALDestroyStringVector( fnamesEndTS );
  XLALDestroyStringVector( detectors );
  XLALDestroyMultiTimestamps( startingTS );
  XLALDestroyMultiTimestamps( endingTS );

  /* Sort the SFDBs from lower GPS timestamp to higher GPS timestamp */
  for ( UINT4 X = 0; X < numIFOs; X++ ) {
    qsort( ( void * )outputSFTs->data[X]->data, outputSFTs->data[X]->length, sizeof( outputSFTs->data[X]->data[0] ), compareSFTepoch );
  }

  return outputSFTs;
} /* XLALReadSFDB() */


static int
read_SFDB_header_from_fp( FILE *fp, SFDBHeader *header )
{

  XLAL_CHECK( fread( &header->det, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO ); // see SFDBDetectors for numbering convention
  XLAL_CHECK( header->det > 0, XLAL_EIO, "Unsupported detector number %d in SFDB.", header->det );
  XLAL_CHECK( header->det < SFDB_DET_LAST, XLAL_EIO, "Unsupported detector number %d in SFDB, highest known number is %d.", header->det, SFDB_DET_LAST - 1 );

  XLAL_CHECK( fread( &header->gps_sec, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO ); // GPS time of this SFDB
  XLAL_CHECK( fread( &header->gps_nsec, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->tbase, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO ); // Coherent time of the SFDB
  XLAL_CHECK( fread( &header->firstfrind, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->nsamples, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO ); // Number of frequency bins
  XLAL_CHECK( fread( &header->red, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO ); // Reductions factor
  XLAL_CHECK( fread( &header->typ, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->n_flag, sizeof( REAL4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->einstein, sizeof( REAL4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->mjdtime, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->nfft, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->wink, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->normd, sizeof( REAL4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->normw, sizeof( REAL4 ), 1, fp ) == 1, XLAL_EIO ); // Normalization factors
  XLAL_CHECK( fread( &header->frinit, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->tsamplu, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO ); // Sampling time
  XLAL_CHECK( fread( &header->deltanu, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO ); // Frequency resolution
  XLAL_CHECK( fread( &header->vx_eq, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->vy_eq, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->vz_eq, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->px_eq, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->py_eq, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->pz_eq, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->n_zeroes, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fread( &header->sat_howmany, sizeof( REAL8 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fseek( fp, 3 * sizeof( REAL8 ), SEEK_CUR ) == 0, XLAL_EIO ); // spare fields
  XLAL_CHECK( fseek( fp, 3 * sizeof( REAL4 ), SEEK_CUR ) == 0, XLAL_EIO ); // spare fields
  XLAL_CHECK( fread( &header->lavesp, sizeof( INT4 ), 1, fp ) == 1, XLAL_EIO );
  XLAL_CHECK( fseek( fp, 2 * sizeof( INT4 ), SEEK_CUR ) == 0, XLAL_EIO ); // more spare spares

  return 0;

} // read_SFDB_header_from_fp


static BOOLEAN
CheckIfSFDBInScienceMode(
  SFDBHeader *header,
  LALStringVector *detectors,
  MultiLIGOTimeGPSVector *startingTS,
  MultiLIGOTimeGPSVector *endingTS
)
{

  /* Get detector index */
  UINT4 detectorIndex = 0;
  while ( strcmp( SFDB_detector_names[header->det],
                  detectors->data[detectorIndex] ) != 0
        ) {
    ++detectorIndex;
  }

  /* Check if SFDB's epoch falls within time stamps */
  /* FIXME: Maybe a binary search is faster? */
  BOOLEAN SFDBWithinScienceMode = 0;
  UINT4 tsIndex = 0;
  REAL8 SFDBEndTime = header->gps_sec + header->tbase;
  while ( header->gps_sec >= startingTS->data[detectorIndex]->data[tsIndex].gpsSeconds ) {
    if ( SFDBEndTime < endingTS->data[detectorIndex]->data[tsIndex].gpsSeconds ) {
      SFDBWithinScienceMode = 1;
      break;
    }
    ++tsIndex;
  }

  return SFDBWithinScienceMode;
}

/// @}
