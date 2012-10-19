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
  
  /* check user input consistency */
  if ( !timeCorrectionFile ) {
    XLALPrintError("%s: invalid NULL input timeCorrectionFile=%p\n", __func__, timeCorrectionFile );
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }
  
  /* prepare output ephemeris struct for returning */
  TimeCorrectionData *tdat;
  if ( ( tdat = XLALCalloc ( 1, sizeof(*tdat) ) ) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( 1, %d ) failed.\n", __func__, sizeof(*tdat) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  
  /* open ephemeris file */
  FILE *fp;
  if ( (fp = LALOpenDataFile ( timeCorrectionFile )) == NULL ) {
    XLALPrintError ("%s: LALOpenDataFile() failed to open '%s' for reading.\n", __func__, timeCorrectionFile );
    XLAL_ERROR_NULL ( XLAL_ESYS );
  }
  
  REAL8 endtime = 0.;
  
  /* read in first header line */
  if ( 4 != fscanf(fp,"%lf %lf %lf %u\n", &tdat->timeCorrStart, &endtime, &tdat->dtTtable, &tdat->nentriesT)) {
    fclose(fp);
    XLALPrintError("%s: couldn't parse first line of %s: %d\n", __func__, timeCorrectionFile );
    XLAL_ERROR_NULL ( XLAL_EDOM );
  }
  
  /* allocate memory for table entries */
  if ( (tvec = XLALCalloc( tdat->nentriesT, sizeof(REAL8) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc(%u, sizeof(REAL8))\n", __func__, tdat->nentriesT );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  
  /* read in table data */
  UINT4 j;
  int ret;
  for (j=0; j < tdat->nentriesT; j++){
    ret = fscanf( fp, "%lf\n", &tvec[j] );

    /* check number of scanned items */
    if (ret != 1) {
      fclose(fp);
      XLALFree( tvec );
      XLALPrintError("%s: Couldn't parse line %d of %s: %d\n", __func__, j+2, timeCorrectionFile, ret);
      XLAL_ERROR_NULL ( XLAL_EDOM );
    }
    
    /* check we've not hit the end of file before reading in all points */
    if ( feof(fp) && j < (tdat->nentriesT)-1 ){
      fclose(fp);
      XLALFree( tvec );
      XLALPrintError("%s: %s does not contain %u lines!\n", __func__, timeCorrectionFile, tdat->nentriesT);
      XLAL_ERROR_NULL ( XLAL_EDOM );
    }
  }
  
  /* close file */
  fclose(fp);
  
  /* set output time delay vector */
  tdat->timeCorrs = tvec;
  
  return tdat;
}

/** Destructor for TimeCorrectionData struct, NULL robust.
 * \ingroup LALBarycenter_h
 */
void XLALDestroyTimeCorrectionData ( TimeCorrectionData *times )
{
  if ( !times )
    return;

  if ( times->timeCorrs )
    XLALFree ( times->timeCorrs );

  XLALFree ( times );

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
  /* check user input consistency */
  if ( !earthEphemerisFile || !sunEphemerisFile ) {
    XLALPrintError ("%s: invalid NULL input earthEphemerisFile=%p, sunEphemerisFile=%p\n", __func__, earthEphemerisFile, sunEphemerisFile );
    XLAL_ERROR_NULL (XLAL_EINVAL );
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
  /* check input consistency */
  if ( !fname ) {
    XLALPrintError ("%s: invalid NULL input\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* open ephemeris file */
  FILE *fp;
  if ( (fp = LALOpenDataFile ( fname )) == NULL ) {
    XLALPrintError ("%s: LALOpenDataFile() failed to open '%s' for reading.\n", __func__, fname );
    XLAL_ERROR_NULL ( XLAL_ESYS );
  }

  INT4 gpsYr; /* gpsYr + leap is the time on the GPS clock
               * at first instant of new year, UTC; equivalently
               * leap is # of leap secs added between Jan.6, 1980 and
               * Jan. 2 of given year
               */
  REAL8 dt;		/* ephemeris-file time-step in seconds */
  UINT4 nEntries;	/* number of ephemeris-file entries */

  /* read first line */
  if ( 3 != fscanf(fp,"%d %le %u\n", &gpsYr, &dt, &nEntries)) {
    fclose(fp);
    XLALPrintError("%s: couldn't parse first line of %s: %d\n", __func__, fname );
    XLAL_ERROR_NULL ( XLAL_EDOM );
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

  /* read the remaining lines */
  UINT4 j;
  int ret;
  for (j=0; j < nEntries; j++)
    {
      ret = fscanf( fp, "%le %le %le %le %le %le %le %le %le %le\n",
                    &ephemV->data[j].gps,
                    &ephemV->data[j].pos[0], &ephemV->data[j].pos[1], &ephemV->data[j].pos[2],
                    &ephemV->data[j].vel[0], &ephemV->data[j].vel[1], &ephemV->data[j].vel[2],
                    &ephemV->data[j].acc[0], &ephemV->data[j].acc[1], &ephemV->data[j].acc[2]);

      /* check number of scanned items */
      if (ret != 10) {
	fclose(fp);
	XLALDestroyEphemerisVector ( ephemV );
	XLALPrintError("%s: Couldn't parse line %d of %s: %d\n", j+2, fname, ret);
	XLAL_ERROR_NULL ( XLAL_EDOM );
      }

      /* check timestamps */
      if(j == 0)
        {
          if (gpsYr - ephemV->data[j].gps > 3600 * 24 * 365 ) {
            XLALPrintError("%s: Wrong timestamp in line %d of %s: %d/%le\n", __func__, j+2, fname, gpsYr, ephemV->data[j].gps );
            fclose(fp);
            XLALDestroyEphemerisVector ( ephemV );
            XLAL_ERROR_NULL ( XLAL_EDOM );
          }
        }
      else
        {
          if (ephemV->data[j].gps != ephemV->data[j-1].gps + ephemV->dt ) {
            XLALPrintError("%s: Wrong timestamp in line %d of %s: %le/%le\n", __func__, j+2, fname, ephemV->data[j].gps, ephemV->data[j-1].gps + ephemV->dt );
            fclose(fp);
            XLALDestroyEphemerisVector ( ephemV );
            XLAL_ERROR_NULL ( XLAL_EDOM );
          }
        }

    } /* for j < nEntries */

  /* check file-sanity: nothing beyond end of table */
  CHAR dummy;
  if ( fscanf (fp,"%c",&dummy) != EOF) {
    XLALPrintError("%s: Garbage at end of ephemeris file %s\n", __func__, fname );
    fclose(fp);
    XLALDestroyEphemerisVector ( ephemV );
    XLAL_ERROR_NULL ( XLAL_EDOM );
  }

  /* done reading, close ephemeris-file file */
  fclose(fp);

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
 * \brief [DEPRECATED] Reads Earth and Sun ephemeris files.
 *
 * This function fills the contents of \a edat from data
 * read from data files, see \a EphemerisData for the definition of this data-type.
 *
 * The function reads in two data files (specified in the
 * edat->ephiles structure) that contain the position, velocity,
 * and acceleration of the Earth and Sun, respectively, at regular
 * intervals througout the specified year. E.g., for 1998, the two files
 * are <tt>earth98.dat</tt> and <tt>sun98.dat</tt>.  These files are derived
 * from the JPL DE405 ephemeris and are provided by Cutler.  The first
 * line of these files specifies the start time, sampling interval, and
 * number of datapoints stored in each file, which are used to allocate
 * data arrays edat->ephemE and edat->ephemS of appropriate
 * length.  LALInitBarycenter() should be called once near the
 * beginning of the larger analysis package, and the fields
 * edat->ephemE and edat->ephemS should be freed with
 * LALFree() near the end.  See the LALBarycenterTest program for an illustration
 * of how this routine is used.
 *
 * \deprecated Use XLALInitBarycenter() instead.
 */
void
LALInitBarycenter ( LALStatus *stat,	/**< LAL-status pointer */
                    EphemerisData *edat	/**< [in/out] initialized ephemeris-data */
                    )
{

    FILE *fp1, *fp2; /* fp1 is table of Earth location; fp2 is for Sun*/
    CHAR dummy;
    INT4 j; /*dummy index*/
    INT4 gpsYr; /*gpsYr + leap is the time on the GPS clock
                          at first instant of new year, UTC; equivalently
                          leap is # of leap secs added between Jan.6, 1980 and
                          Jan. 2 of given year */
    INT4 ret; /* return value for checking */

    INITSTATUS(stat);
    ATTATCHSTATUSPTR(stat);

    /* open earth file */
    fp1 = LALOpenDataFile(edat->ephiles.earthEphemeris);

    /* check that we could open the file */
    if ( fp1 == NULL ) {
      snprintf (errmsg, ERRMSGLEN, "%s '%s'\n", LALINITBARYCENTERH_MSGEOPEN, edat->ephiles.earthEphemeris);
      errmsg[ERRMSGLEN-1] = '\0';
      ABORT (stat, LALINITBARYCENTERH_EOPEN, errmsg);
    }

    /* read first line */
    ret = fscanf(fp1,"%d %le %d\n", &gpsYr, &edat->dtEtable, &edat->nentriesE);
    if (ret != 3) {
      fclose(fp1);
      XLALPrintError("couldn't parse first line of %s: %d\n", edat->ephiles.earthEphemeris, ret);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }

    /* allocate memory for ephemeris info */
    edat->ephemE  = (PosVelAcc *)LALMalloc(edat->nentriesE*sizeof(PosVelAcc));
    if (edat->ephemE == NULL) {
      fclose(fp1);
      ABORT(stat, LALINITBARYCENTERH_EMEM, LALINITBARYCENTERH_MSGEMEM);
    }

    /* first column in earth.dat or sun.dat is gps time--one long integer
       giving the number of secs that have ticked since start of GPS epoch
       +  on 1980 Jan. 6 00:00:00 UTC
    */

    /* read the remaining lines */
    for (j=0; j < edat->nentriesE; ++j) {
      ret = fscanf(fp1,"%le %le %le %le %le %le %le %le %le %le\n",
		   &edat->ephemE[j].gps,
		   &edat->ephemE[j].pos[0], &edat->ephemE[j].pos[1], &edat->ephemE[j].pos[2],
		   &edat->ephemE[j].vel[0], &edat->ephemE[j].vel[1], &edat->ephemE[j].vel[2],
		   &edat->ephemE[j].acc[0], &edat->ephemE[j].acc[1], &edat->ephemE[j].acc[2]);

      /* check number of scanned items */
      if (ret != 10) {
	fclose(fp1);
	LALFree(edat->ephemE);
	XLALPrintError("Couldn't parse line %d of %s: %d\n", j+2, edat->ephiles.earthEphemeris, ret);
	ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
      }

      /* check timestamps */
      if(j == 0) {
	if (gpsYr - edat->ephemE[j].gps > 3600 * 24 * 365) {
	  XLALPrintError("Wrong timestamp in line %d of %s: %d/%le\n",
			j+2, edat->ephiles.earthEphemeris, gpsYr, edat->ephemE[j].gps);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      } else {
	if (edat->ephemE[j].gps != edat->ephemE[j-1].gps + edat->dtEtable) {
	  XLALPrintError("Wrong timestamp in line %d of %s: %le/%le\n",
			j+2, edat->ephiles.earthEphemeris, edat->ephemE[j].gps, edat->ephemE[j-1].gps + edat->dtEtable);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      }

      /* check position, velocity and acceleration */
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
      {
	REAL8 length;
	length = sqrt(SQR(edat->ephemE[j].pos[0]) + SQR(edat->ephemE[j].pos[1]) + SQR(edat->ephemE[j].pos[2]));
	if ( fabs(499.0 - length) > 25) /* 5% */ {
	  XLALPrintError("earth position out of range in line %d of %s: %le %le %le: %le\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemE[j].pos[0], edat->ephemE[j].pos[1], edat->ephemE[j].pos[2], length);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
	length = sqrt(SQR(edat->ephemE[j].vel[0]) + SQR(edat->ephemE[j].vel[1]) + SQR(edat->ephemE[j].vel[2]));
	if (fabs (1e-4 - length) > 1e-5) /* 10% */ {
	  XLALPrintError("earth velocity out of range in line %d of %s: %le %le %le: %le\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemE[j].vel[0], edat->ephemE[j].vel[1], edat->ephemE[j].vel[2], length);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
	length = sqrt(SQR(edat->ephemE[j].acc[0]) + SQR(edat->ephemE[j].acc[1]) + SQR(edat->ephemE[j].acc[2]));
	if (fabs (2e-11 - length) > 3e-12) /* 15% */ {
	  XLALPrintError("earth acceleration out of range in line %d of %s: %le %le %le: %le\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemE[j].acc[0], edat->ephemE[j].acc[1], edat->ephemE[j].acc[2], length);
	  fclose(fp1);
	  LALFree(edat->ephemE);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      }

      /* debug
      {
	REAL8 length;
	fprintf(stderr,"earth line: %d:  ");
	length = SQR(edat->ephemE[j].pos[0]) + SQR(edat->ephemE[j].pos[1]) + SQR(edat->ephemE[j].pos[2]);
	fprintf(stderr,"pos: %le (%le), ", sqrt(length), length);
	length = SQR(edat->ephemE[j].vel[0]) + SQR(edat->ephemE[j].vel[1]) + SQR(edat->ephemE[j].vel[2]);
	fprintf(stderr,"vel: %le (%le), ", sqrt(length), length);
	length = SQR(edat->ephemE[j].acc[0]) + SQR(edat->ephemE[j].acc[1]) + SQR(edat->ephemE[j].acc[2]);
	fprintf(stderr,"acc: %le (%le)\n", sqrt(length), length);
      }
      */
    }

    if (fscanf(fp1,"%c",&dummy) != EOF) {
      XLALPrintError("Garbage at end of ephemeris file %s\n", edat->ephiles.earthEphemeris);
      fclose(fp1);
      LALFree(edat->ephemE);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }

    /* close earth file */
    fclose(fp1);


    /* open sun file */
    fp2 = LALOpenDataFile(edat->ephiles.sunEphemeris);

    /* check that we could open the file */
    if ( fp2 == NULL ) {
      LALFree(edat->ephemE);
      snprintf (errmsg, ERRMSGLEN, "%s '%s'\n", LALINITBARYCENTERH_MSGEOPEN, edat->ephiles.sunEphemeris);
      errmsg[ERRMSGLEN-1] = 0;
      ABORT (stat, LALINITBARYCENTERH_EOPEN, errmsg);
    }

    /* read first line */
    ret = fscanf(fp2,"%d %le %d\n", &gpsYr, &edat->dtStable, &edat->nentriesS);
    if (ret != 3) {
      LALFree(edat->ephemE);
      fclose(fp2);
      XLALPrintError("Couldn't parse first line of %s: %d\n", edat->ephiles.sunEphemeris, ret);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }

    /* allocate memory for ephemeris info */
    edat->ephemS  = (PosVelAcc *)LALMalloc(edat->nentriesS*sizeof(PosVelAcc));
    if (edat->ephemS == NULL) {
      fclose(fp2);
      LALFree(edat->ephemE);
      ABORT(stat, LALINITBARYCENTERH_EMEM, LALINITBARYCENTERH_MSGEMEM);
    }

    /* read the remaining lines */
    for (j=0; j < edat->nentriesS; ++j) {
      ret = fscanf(fp2,"%le %le %le %le %le %le %le %le %le %le\n",
		   &edat->ephemS[j].gps,
		   &edat->ephemS[j].pos[0], &edat->ephemS[j].pos[1], &edat->ephemS[j].pos[2],
		   &edat->ephemS[j].vel[0], &edat->ephemS[j].vel[1], &edat->ephemS[j].vel[2],
		   &edat->ephemS[j].acc[0], &edat->ephemS[j].acc[1], &edat->ephemS[j].acc[2]);

      /* check number of scanned items */
      if (ret != 10) {
	fclose(fp2);
	LALFree(edat->ephemE);
	LALFree(edat->ephemS);
	XLALPrintError("Couldn't parse line %d of %s: %d\n", j+2, edat->ephiles.sunEphemeris, ret);
	ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
      }

      /* check timestamps */
      if(j == 0) {
	if (gpsYr - edat->ephemS[j].gps > 3600 * 24 * 365) {
	  XLALPrintError("Wrong timestamp in line %d of %s: %d/%le\n",
			j+2, edat->ephiles.sunEphemeris, gpsYr, edat->ephemS[j].gps);
	  fclose(fp2);
	  LALFree(edat->ephemE);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      } else {
	if (edat->ephemS[j].gps != edat->ephemS[j-1].gps + edat->dtStable) {
	  XLALPrintError("Wrong timestamp in line %d of %s: %le/%le\n",
			j+2, edat->ephiles.sunEphemeris, edat->ephemS[j].gps, edat->ephemS[j-1].gps + edat->dtStable);
	  fclose(fp2);
	  LALFree(edat->ephemE);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      }

      /* check position, velocity and acceleration */
      {
	REAL8 length;
	length = sqrt(SQR(edat->ephemS[j].pos[0]) + SQR(edat->ephemS[j].pos[1]) + SQR(edat->ephemS[j].pos[2]));
	if ((1 > length) || (length > 10)) {
	  XLALPrintError("sun position out of range in line %d of %s: %f %f %f: %f\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemS[j].pos[0], edat->ephemS[j].pos[1], edat->ephemS[j].pos[2], length);
	  fclose(fp2);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
	length = sqrt(SQR(edat->ephemS[j].vel[0]) + SQR(edat->ephemS[j].vel[1]) + SQR(edat->ephemS[j].vel[2]));
	if ((1e-8 > length) || (length > 1e-7)) {
	  XLALPrintError("sun velocity out of range in line %d of %s: %f %f %f: %f\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemS[j].vel[0], edat->ephemS[j].vel[1], edat->ephemS[j].vel[2], length);
	  fclose(fp2);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
	length = sqrt(SQR(edat->ephemS[j].acc[0]) + SQR(edat->ephemS[j].acc[1]) + SQR(edat->ephemS[j].acc[2]));
	if ((1e-16 > length) || (length > 1e-14)) {
	  XLALPrintError("sun acceleration out of range in line %d of %s: %f %f %f: %f\n",
			j+2, edat->ephiles.earthEphemeris,
			edat->ephemS[j].acc[0], edat->ephemS[j].acc[1], edat->ephemS[j].acc[2], length);
	  fclose(fp2);
	  LALFree(edat->ephemS);
	  ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
	}
      }
      /* debug
      {
	REAL8 length;
	fprintf(stderr,"sun line: %d:  ");
	length = SQR(edat->ephemS[j].pos[0]) + SQR(edat->ephemS[j].pos[1]) + SQR(edat->ephemS[j].pos[2]);
	fprintf(stderr,"pos: %le (%le), ", sqrt(length), length);
	length = SQR(edat->ephemS[j].vel[0]) + SQR(edat->ephemS[j].vel[1]) + SQR(edat->ephemS[j].vel[2]);
	fprintf(stderr,"vel: %le (%le), ", sqrt(length), length);
	length = SQR(edat->ephemS[j].acc[0]) + SQR(edat->ephemS[j].acc[1]) + SQR(edat->ephemS[j].acc[2]);
	fprintf(stderr,"acc: %le (%le)\n", sqrt(length), length);
      }
      */
    }

    if (fscanf(fp2,"%c",&dummy) != EOF) {
      XLALPrintError("Garbage at end of ephemeris file %s\n", edat->ephiles.sunEphemeris);
      fclose(fp2);
      LALFree(edat->ephemE);
      LALFree(edat->ephemS);
      ABORT(stat, LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
    }

    /* close the file */
    fclose(fp2);

    /* successful return */
    DETATCHSTATUSPTR(stat);
    RETURN(stat);

} /* LALInitBarycenter() */
