/*
 * Copyright (C) 2010 Reinhard Prix (XLALified)
*  Copyright (C) 2007 Curt Cutler, Jolien Creighton, Reinhard Prix, Teviet Creighton, Bernd Machenschalk
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

#include <lal/FileIO.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/ConfigFile.h>
#include <lal/LALString.h>

/** \cond DONT_DOXYGEN */

/* ----- defines and macros ---------- */
#define SQ(x) ((x) * (x))
#define NORM3D(x) ( SQ( (x)[0]) + SQ( (x)[1] ) + SQ ( (x)[2] ) )
#define LENGTH3D(x) ( sqrt( NORM3D ( (x) ) ) )

/** \endcond */

/* ----- local type definitions ---------- */

/** Generic ephemeris-vector type, holding one timeseries of pos, vel, acceleration.
 * This is used for the generic ephemeris-reader XLAL-function, at the end the resulting
 * ephemeris-data  will be stored in the 'old type \a EphemerisData for backwards compatibility.
 */
typedef struct
{
  UINT4 length;      	/**< number of ephemeris-data entries */
  REAL8 dt;      	/**< spacing in seconds between consecutive instants in ephemeris table.*/
  PosVelAcc *data;    	/**< array containing pos,vel,acc as extracted from ephem file. Units are sec, 1, 1/sec respectively */
}
EphemerisVector;

/* ----- internal prototypes ---------- */
EphemerisVector *XLALCreateEphemerisVector ( UINT4 length );
void XLALDestroyEphemerisVector ( EphemerisVector *ephemV );

EphemerisVector * XLALReadEphemerisFile ( const CHAR *fname);
int XLALCheckEphemerisRanges ( const EphemerisVector *ephemEarth, REAL8 avg[3], REAL8 range[3] );

/* ----- function definitions ---------- */

/* ========== exported API ========== */

/** An XLAL interface for reading a time correction file containing a table
 * of values for converting between Terrestrial Time TT (or TDT) to either
 * TDB (i.e. the file contains time corrections related to the Einstein delay)
 * or Teph (a time system from Irwin and Fukushima, 1999, closely related to 
 * Coordinate Barycentric Time, TCB) depending on the file.
 *
 * The file contains a header with the GPS start time, GPS end time, time 
 * interval between subsequent entries (seconds), and the number of entries.
 * 
 * The rest of the file contains a list of the time delays (in seconds).
 *
 * The tables for the conversion to TDB and Teph are derived from the ephemeris
 * file TDB.1950.2050 and TIMEEPH_short.te405 within TEMPO2 
 * http://www.atnf.csiro.au/research/pulsar/tempo2/. They are created from the
 * Chebychev polynomials in these files using the conversion in the lalapps code
 * lalapps_create_time_correction_ephemeris
 *
 * \ingroup LALBarycenter_h
 */
TimeCorrectionData *XLALInitTimeCorrections ( const CHAR *timeCorrectionFile /**< File containing Earth's position.  */ ){
  REAL8 *tvec = NULL; /* create time vector */
  LALParsedDataFile *flines = NULL;
  UINT4 numLines = 0, j = 0;
  REAL8 endtime = 0.;
  
  /* check user input consistency */
  if ( !timeCorrectionFile ) {
    XLALPrintError("%s: invalid NULL input timeCorrectionFile=%p\n", __func__, timeCorrectionFile );
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }
  
  /* read in file with XLALParseDataFile to ignore comment header lines */
  if ( XLALParseDataFile ( &flines, timeCorrectionFile ) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  /* prepare output ephemeris struct for returning */
  TimeCorrectionData *tdat;
  if ( ( tdat = XLALCalloc ( 1, sizeof(*tdat) ) ) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( 1, %d ) failed.\n", __func__, sizeof(*tdat) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  numLines = flines->lines->nTokens;
  
  /* read in info from first line (an uncommented header) */
  if ( 4 != sscanf(flines->lines->tokens[0], "%lf %lf %lf %u",
       &tdat->timeCorrStart, &endtime, &tdat->dtTtable, &tdat->nentriesT) ){
    XLALDestroyParsedDataFile ( flines );
    XLALDestroyTimeCorrectionData( tdat );
    XLALPrintError("%s: couldn't parse first line of %s: %d\n", __func__,
                   timeCorrectionFile );
    XLAL_ERROR_NULL ( XLAL_EDOM );
  }
  
  if( numLines != tdat->nentriesT+1 ){
    XLALDestroyParsedDataFile ( flines );
    XLALDestroyTimeCorrectionData( tdat );
    XLALPrintError ("%s: header line inconsistent with number of lines in \
file.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  
  /* allocate memory for table entries */
  if ( (tvec = XLALCalloc( tdat->nentriesT, sizeof(REAL8) )) == NULL ) {
    XLALDestroyParsedDataFile ( flines );
    XLALDestroyTimeCorrectionData( tdat );
    XLALPrintError ("%s: XLALCalloc(%u, sizeof(REAL8))\n", __func__,
                    tdat->nentriesT );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  
  /* read in table data */
  int ret;
  for (j=1; j < numLines; j++){
    ret = sscanf( flines->lines->tokens[j], "%lf", &tvec[j-1] );

    /* check number of scanned items */
    if (ret != 1) {
      XLALFree( tvec );
      XLALDestroyParsedDataFile ( flines );
      XLALDestroyTimeCorrectionData( tdat );
      XLALPrintError("%s: Couldn't parse line %d of %s: %d\n", __func__, j+2, timeCorrectionFile, ret);
      XLAL_ERROR_NULL ( XLAL_EDOM );
    }
  }

  XLALDestroyParsedDataFile ( flines );

  /* set output time delay vector */
  tdat->timeCorrs = tvec;
  
  return tdat;
}

/** Destructor for TimeCorrectionData struct, NULL robust.
 * \ingroup LALBarycenter_h
 */
void XLALDestroyTimeCorrectionData ( TimeCorrectionData *tcd )
{
  if ( !tcd )
    return;

  if ( tcd->timeCorrs )
    XLALFree ( tcd->timeCorrs );

  XLALFree ( tcd );

  return;

} /* XLALDestroyTimeCorrectionData() */

/** XLAL interface to reading ephemeris files 'earth' and 'sun', and return
 * ephemeris-data in old backwards-compatible type \a EphemerisData
 *
 * These ephemeris data files contain arrays
 * of center-of-mass positions for the Earth and Sun, respectively.
 *
 * The tables are derived from the JPL ephemeris.
 *
 * Files tabulate positions for one calendar year
 * (actually, a little more than one year, to deal
 * with overlaps).  The first line of each table summarizes
 * what is in it. Subsequent lines give the time (GPS) and the
 * Earth's position \f$(x,y,z)\f$,
 * velocity \f$(v_x, v_y, v_z)\f$, and acceleration \f$(a_x, a_y, a_z)\f$
 * at that instant.  All in units of seconds; e.g. positions have
 * units of seconds, and accelerations have units 1/sec.
 *
 * \ingroup LALBarycenter_h
 */
EphemerisData *
XLALInitBarycenter ( const CHAR *earthEphemerisFile,         /**< File containing Earth's position.  */
                     const CHAR *sunEphemerisFile            /**< File containing Sun's position. */
                     )
{
  EphemerisType etype;

  /* check user input consistency */
  if ( !earthEphemerisFile || !sunEphemerisFile ) {
    XLALPrintError ("%s: invalid NULL input earthEphemerisFile=%p, sunEphemerisFile=%p\n", __func__, earthEphemerisFile, sunEphemerisFile );
    XLAL_ERROR_NULL (XLAL_EINVAL );
  }

  /* check the ephemeris type from file name for consistency */
  if ( strstr( earthEphemerisFile, "DE200" ) ){
    if ( !strstr( sunEphemerisFile, "DE200" ) ){
      XLALPrintError("%s: %s and %s have inconsistent ephemeris type\n",
                     __func__, earthEphemerisFile, sunEphemerisFile );
      XLAL_ERROR_NULL( XLAL_EINVAL );
    }
    else etype = EPHEM_DE200;
  }
  else if ( strstr( earthEphemerisFile, "DE405" ) ){
    if ( !strstr( sunEphemerisFile, "DE405" ) ){
      XLALPrintError("%s: %s and %s have inconsistent ephemeris type\n",
                     __func__, earthEphemerisFile, sunEphemerisFile );
      XLAL_ERROR_NULL( XLAL_EINVAL );
    }
    else etype = EPHEM_DE405;
  }
  else if ( strstr( earthEphemerisFile, "DE414" ) ){
    if ( !strstr( sunEphemerisFile, "DE414" ) ){
      XLALPrintError("%s: %s and %s have inconsistent ephemeris type\n",
                     __func__, earthEphemerisFile, sunEphemerisFile );
      XLAL_ERROR_NULL( XLAL_EINVAL );
    }
    else etype = EPHEM_DE414;
  }
  else{
    if ( strstr( sunEphemerisFile, "DE200" ) ||
         strstr( sunEphemerisFile, "DE405" ) ||
         strstr( sunEphemerisFile, "DE414" ) ){
      XLALPrintError("%s: %s and %s have inconsistent ephemeris type\n",
                     __func__, earthEphemerisFile, sunEphemerisFile );
      XLAL_ERROR_NULL( XLAL_EINVAL );
    }
    else /* if no ephemeris type extension is found default to DE405 */
      etype = EPHEM_DE405;
  }

  EphemerisVector *ephemV;

  /* ----- read EARTH ephemeris file ---------- */
  if ( ( ephemV = XLALReadEphemerisFile ( earthEphemerisFile )) == NULL ) {
    XLALPrintError ("%s: XLALReadEphemerisFile('%s') failed\n", __func__, earthEphemerisFile );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* typical position, velocity and acceleration and allowed ranged */
  REAL8 avgE[3] = {499.0,  1e-4, 2e-11 };
  REAL8 rangeE[3] = {25.0, 1e-5, 3e-12 };

  if ( XLALCheckEphemerisRanges ( ephemV, avgE, rangeE ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: Earth-ephemeris range error!\n", __func__ );
    XLALDestroyEphemerisVector ( ephemV );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* prepare output ephemeris struct for returning */
  EphemerisData *edat;
  if ( ( edat = XLALCalloc ( 1, sizeof(*edat) ) ) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( 1, %d ) failed.\n", __func__, sizeof(*edat) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  /* store in ephemeris-struct */
  edat->nentriesE = ephemV->length;
  edat->dtEtable  = ephemV->dt;
  edat->ephemE    = ephemV->data;
  edat->etype     = etype;
  XLALFree ( ephemV );	/* don't use 'destroy', as we linked the data into edat! */
  ephemV = NULL;

  /* ----- read SUN ephemeris file ---------- */
  if ( ( ephemV = XLALReadEphemerisFile ( sunEphemerisFile )) == NULL ) {
    XLALPrintError ("%s: XLALReadEphemerisFile('%s') failed\n", __func__, sunEphemerisFile );
    XLALDestroyEphemerisData ( edat );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* typical position, velocity and acceleration and allowed ranged */
  REAL8 avgS[3]   = { 5.5, 5.5e-8, 50.5e-16 };
  REAL8 rangeS[3] = { 4.5, 4.5e-8, 49.5e-16 };

  if ( XLALCheckEphemerisRanges ( ephemV, avgS, rangeS ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: Sun-ephemeris range error!\n", __func__ );
    XLALDestroyEphemerisVector ( ephemV );
    XLALDestroyEphemerisData ( edat );
    XLAL_ERROR_NULL ( XLAL_EDOM );
  }

  /* store in ephemeris-struct */
  edat->nentriesS = ephemV->length;
  edat->dtStable  = ephemV->dt;
  edat->ephemS    = ephemV->data;
  XLALFree ( ephemV );	/* don't use 'destroy', as we linked the data into edat! */
  ephemV = NULL;

  // store *copy* of ephemeris-file names in output structure
  edat->ephiles.earthEphemeris = XLALStringDuplicate( earthEphemerisFile );
  edat->ephiles.sunEphemeris   = XLALStringDuplicate( sunEphemerisFile );

  /* return resulting ephemeris-data */
  return edat;

} /* XLALInitBarycenter() */


/** Destructor for EphemerisData struct, NULL robust.
 * \ingroup LALBarycenter_h
 */
void
XLALDestroyEphemerisData ( EphemerisData *edat )
{
  if ( !edat )
    return;

  if ( edat->ephiles.earthEphemeris )
    XLALFree ( edat->ephiles.earthEphemeris );

  if ( edat->ephiles.sunEphemeris )
    XLALFree ( edat->ephiles.sunEphemeris );

  if ( edat->ephemE )
    XLALFree ( edat->ephemE );

  if ( edat->ephemS )
    XLALFree ( edat->ephemS );

  XLALFree ( edat );

  return;

} /* XLALDestroyEphemerisData() */


/* ========== internal function definitions ========== */

/** simple creator function for EphemerisVector type */
EphemerisVector *
XLALCreateEphemerisVector ( UINT4 length )
{
  EphemerisVector * ret;
  if ( ( ret = XLALCalloc ( 1, sizeof (*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc(1, %d)\n", sizeof (*ret) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( ( ret->data = XLALCalloc ( length, sizeof(*ret->data) ) ) == NULL ) {
    XLALFree ( ret );
    XLALPrintError ("%s: failed to XLALCalloc (%d, %d)\n", __func__, length, sizeof(*ret->data) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  ret->length = length;

  return ret;

} /* XLALCreateEphemerisVector() */

/** Destructor for EphemerisVector, NULL robust.
 */
void
XLALDestroyEphemerisVector ( EphemerisVector *ephemV )
{
  if ( !ephemV )
    return;

  if ( ephemV->data )
    XLALFree ( ephemV->data );

  XLALFree ( ephemV );

  return;

} /* XLALDestroyEphemerisVector() */


/** XLAL function to read ephemeris-data from one file, returning a EphemerisVector.
 * This is a helper-function to XLALInitBarycenter().
 */
EphemerisVector *
XLALReadEphemerisFile ( const CHAR *fname )
{
  UINT4 j = 0, numLines = 0;
  LALParsedDataFile *flines = NULL;

  /* check input consistency */
  if ( !fname ) {
    XLALPrintError ("%s: invalid NULL input\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* read in file with XLALParseDataFile to ignore comment header lines */
  if ( XLALParseDataFile ( &flines, fname ) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  numLines = flines->lines->nTokens;

  INT4 gpsYr; /* gpsYr + leap is the time on the GPS clock
               * at first instant of new year, UTC; equivalently
               * leap is # of leap secs added between Jan.6, 1980 and
               * Jan. 2 of given year
               */
  REAL8 dt;		/* ephemeris-file time-step in seconds */
  UINT4 nEntries;	/* number of ephemeris-file entries */

  /* read first line */
  if ( 3 != sscanf(flines->lines->tokens[0],"%d %le %u\n", &gpsYr, &dt,
     &nEntries)) {
    XLALDestroyParsedDataFile( flines );
    XLALPrintError("%s: couldn't parse first line of %s: %d\n", __func__, fname );
    XLAL_ERROR_NULL ( XLAL_EDOM );
  }

  /* check that number of lines is correct */
  if( nEntries != (numLines - 1)/4 ){
    XLALDestroyParsedDataFile( flines );
    XLALPrintError("%s: inconsistent number of lines in file compared to \
header information\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* prepare output ephemeris vector */
  EphemerisVector *ephemV;
  if ( (ephemV = XLALCreateEphemerisVector ( nEntries )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateEphemerisVector(%d)\n", __func__, nEntries );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  ephemV->dt = dt;

  /* first column in ephemeris-file is gps time--one long integer
   * giving the number of secs that have ticked since start of GPS epoch
   * +  on 1980 Jan. 6 00:00:00 UTC
   */

  /* the ephemeris files are created with each entry spanning 4 lines with the
   * format:
   *  gps\tposX\tposY\n
   *  posZ\tvelX\tvelY\n
   *  velZ\taccX\taccY\n
   *  accZ\n
   ***************************************************************************/

  /* read the remaining lines */
  int ret;
  for (j=0; j < nEntries; j++)
    {
      CHAR *oneline = NULL;

      int i;
      /* concatenate the tokens representing one entry */
      for ( i = 0; i < 4; i++ ){
        oneline = XLALStringAppend( oneline, flines->lines->tokens[4*j+1+i] );
        /* add space between tokens*/
        oneline = XLALStringAppend( oneline, " " );
      }

      ret = sscanf( oneline, "%le %le %le %le %le %le %le %le %le %le",
                    &ephemV->data[j].gps, &ephemV->data[j].pos[0],
                    &ephemV->data[j].pos[1], &ephemV->data[j].pos[2],
                    &ephemV->data[j].vel[0], &ephemV->data[j].vel[1],
                    &ephemV->data[j].vel[2], &ephemV->data[j].acc[0],
                    &ephemV->data[j].acc[1], &ephemV->data[j].acc[2] );

      XLALFree( oneline );

      if (ret != 10) {
        XLALDestroyEphemerisVector ( ephemV );
        XLALDestroyParsedDataFile( flines );
        XLALPrintError("%s: Couldn't parse line %d of %s: %d\n", j+2, fname,
                       ret);
        XLAL_ERROR_NULL ( XLAL_EDOM );
      }

      /* check timestamps */
      if(j == 0)
        {
          if (gpsYr - ephemV->data[j].gps > 3600 * 24 * 365 ) {
            XLALPrintError("%s: Wrong timestamp in line %d of %s: %d/%le\n", __func__, j+2, fname, gpsYr, ephemV->data[j].gps );
            XLALDestroyEphemerisVector ( ephemV );
            XLAL_ERROR_NULL ( XLAL_EDOM );
          }
        }
      else
        {
          if (ephemV->data[j].gps != ephemV->data[j-1].gps + ephemV->dt ) {
            XLALPrintError("%s: Wrong timestamp in line %d of %s: %le/%le\n", __func__, j+2, fname, ephemV->data[j].gps, ephemV->data[j-1].gps + ephemV->dt );
            XLALDestroyEphemerisVector ( ephemV );
            XLAL_ERROR_NULL ( XLAL_EDOM );
          }
        }

    } /* for j < nEntries */

  XLALDestroyParsedDataFile( flines );

  /* return result */
  return ephemV;

} /* XLALReadEphemerisFile() */


/** Function to check rough consistency of ephemeris-data with being an actual
 * 'Earth' ephemeris: ie check position, velocity and acceleration are within
 * reasonable ranges {avg +- range}. where 'avg' and 'range' are 3-D arrays
 * with [0]=position, [1]=velocity and [2]=acceleration
 */
int
XLALCheckEphemerisRanges ( const EphemerisVector *ephemV, REAL8 avg[3], REAL8 range[3] )
{
  /* check input consistency */
  if ( !ephemV ) {
    XLALPrintError ("%s: invalid NULL input \n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  UINT4 numEntries = ephemV->length;
  REAL8 dt = ephemV->dt;

  /* check position, velocity and acceleration */
  UINT4 j;
  REAL8 tjm1 = 0;
  for ( j=0; j < numEntries; j ++ )
    {
      REAL8 length;
      length = LENGTH3D ( ephemV->data[j].pos );
      if ( fabs( avg[0] - length) >  range[0] ) {
        XLALPrintError("%s: position out of range in entry %d: vr=(%le, %le, %le), sqrt{|vr|} = %le [%g +- %g]\n", __func__,
                       j, ephemV->data[j].pos[0], ephemV->data[j].pos[1], ephemV->data[j].pos[2], length, avg[0], range[0] );
        XLAL_ERROR ( XLAL_EDOM );
      }
      length = LENGTH3D ( ephemV->data[j].vel );
      if ( fabs(avg[1] - length) > range[1] ) /* 10% */ {
        XLALPrintError("%s: velocity out of range in entry %d: vv=(%le, %le, %le), sqrt{|vv|} = %le, [%g +- %g]\n", __func__,
                       j, ephemV->data[j].vel[0], ephemV->data[j].vel[1], ephemV->data[j].vel[2], length, avg[1], range[1] );
        XLAL_ERROR ( XLAL_EDOM );
      }
      length = LENGTH3D ( ephemV->data[j].acc );
      if ( fabs(avg[2] - length) > range[2] ) /* 15% */ {
        XLALPrintError("%s: acceleration out of range in entry %d: va=(%le, %le, %le), sqrt{|va|} = %le, [%g +- %g]\n", __func__,
                       j, ephemV->data[j].acc[0], ephemV->data[j].acc[1], ephemV->data[j].acc[2], length, avg[2], range[2] );
        XLAL_ERROR ( XLAL_EDOM );
      }

      /* check timestep */
      if ( j > 0 ) {
        if ( ephemV->data[j].gps - tjm1 != dt ) {
          XLALPrintError ("%s: invalid timestep in entry %d: t_i - t_{i-1} = %g != %g\n", __func__, j, ephemV->data[j].gps - tjm1, dt );
          XLAL_ERROR ( XLAL_EDOM );
        }
      }
      tjm1 = ephemV->data[j].gps;	/* keep track of previous timestamp */

    } /* for j < nEntries */


  /* all seems ok */
  return XLAL_SUCCESS;

} /* XLALCheckEphemerisRanges() */

/* ============================= deprecated LAL interface ============================== */
/** \cond DONTDOXYGEN */
#define ERRMSGLEN 512
CHAR errmsg[ERRMSGLEN];	/* string-buffer for more explicit error-messages */
/** \endcond */

/** \ingroup LALBarycenter_h
 * \brief [DEPRECATED] Reads Earth and Sun ephemeris files. Simple wrapper around XLALInitBarycenter()
 * \deprecated Use XLALInitBarycenter() instead.
 */
void
LALInitBarycenter ( LALStatus *stat,	/**< LAL-status pointer */
                    EphemerisData *edat	/**< [in/out] initialized ephemeris-data */
                    )
{
    INITSTATUS(stat);

    if( edat == NULL )
      ABORT( stat, LALINITBARYCENTERH_EOPEN, "Ephemeris structure is NULL" );

    if( edat->ephiles.earthEphemeris == NULL )
      ABORT( stat, LALINITBARYCENTERH_EOPEN, LALINITBARYCENTERH_MSGEOPEN );

    if( edat->ephiles.sunEphemeris == NULL )
      ABORT( stat, LALINITBARYCENTERH_EOPEN, LALINITBARYCENTERH_MSGEOPEN );

    /* use XLALInitBarycenter */
    EphemerisData *edattmp = XLALInitBarycenter( edat->ephiles.earthEphemeris, edat->ephiles.sunEphemeris );
    if( edattmp == NULL ){
      ABORT( stat, LALINITBARYCENTERH_EOPEN, LALINITBARYCENTERH_MSGEOPEN );
    }

    // We need to be careful about returning this, due to the unfortunate input/output method
    // of this deprecated LAL function: 'edat' is both used as input and output, where only
    // the 'ephiles' entry is supposed to be initialized in the input data, and so we need to
    // preserve that entry:
    EphemerisFilenames tmp;
    memcpy ( &tmp, &edat->ephiles, sizeof(edat->ephiles) );
    // now copy new ephemeris-struct over the input one:
    memcpy ( edat, edattmp, sizeof(*edat) );
    // restore the original 'ephiles' pointer entries
    memcpy ( &edat->ephiles, &tmp, sizeof(edat->ephiles) );
    // free the new 'ephiles' strings (allocated in XLALInitBarycenter())
    XLALFree ( edattmp->ephiles.earthEphemeris );
    XLALFree ( edattmp->ephiles.sunEphemeris );
    // and free the interal 'edattmp' container (but *not* the ephemE/ephemS arrays, which we return in 'edat' !)
    XLALFree ( edattmp );

    /* successful return */
    RETURN(stat);

} /* LALInitBarycenter() */
