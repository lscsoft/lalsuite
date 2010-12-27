/*
*  Copyright (C) 2007 Curt Cutler, David Chin, Jolien Creighton, Reinhard Prix, Teviet Creighton
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

/**
 * \author Curt Cutler
 * \date 2001
 * \file
 * \ingroup LALBarycenter_h
 * \brief  Tests the routine LALBarycenter(): exercises some of the error
 * conditions and makes sure that they work.
 *
 * \heading{Program <tt>LALBarycenterTest.c</tt>}
 *
 * \heading{Usage}
 * \code
 * LALBarycenterTest
 * \endcode
 *
 * \heading{Description}
 *
 * This program demonstrates the use of LALBarycenter.c.
 * The two ephemeris files specified in the EphemerisFilenames
 * structure (e.g., for data taken in 1998, <tt>sun98.dat</tt> and <tt>earth98.dat</tt>)
 * are assumed to be in the same directory as the program as
 * the test program.
 *
 */

#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/DetectorSite.h>
#include <lal/Date.h>

/** \name Error codes */
/*@{*/
#define LALBARYCENTERTESTC_ENOM 0
#define LALBARYCENTERTESTC_EOPEN 1
#define LALBARYCENTERTESTC_EOUTOFRANGEE  4
#define LALBARYCENTERTESTC_EOUTOFRANGES  8
#define LALBARYCENTERTESTC_EBADSOURCEPOS 16
#define LALBARYCENTERTESTC_EEPHFILE 32

#define LALBARYCENTERTESTC_MSGENOM "Nominal exit"
#define LALBARYCENTERTESTC_MSGEOPEN "Error checking failed to catch missing ephemeris file."
#define LALBARYCENTERTESTC_MSGEOUTOFRANGEE "Failed to catch that tgps not in range of earth.dat file"
#define LALBARYCENTERTESTC_MSGEOUTOFRANGES "Failed to catch that tgps not in range of sun.dat file"
#define LALBARYCENTERTESTC_MSGEBADSOURCEPOS "Failed to catch bad source position"
#define LALBARYCENTERTESTC_MSGEEPHFILE "Failed to catch error reading ephemeris file."
/*@}*/

/** \cond DONT_DOXYGEN */
NRCSID(LALBARYCENTERTESTC,"$Id$");

/* ----- internal prototype ---------- */
int compare_ephemeris ( const EphemerisData *edat1, const EphemerisData *edat2 );

/*
  int lalDebugLevel=0;
*/
  INT4 lalDebugLevel=7;

BarycenterInput baryinput;
EmissionTime  emit;
LIGOTimeGPS tGPS;
EarthState earth;

INT4 t2000 = 630720013; /* gps time at Jan 1, 2000 00:00:00 UTC */
INT4 t1998 = 630720013-730*86400-1;/* gps at Jan 1,1998 00:00:00 UTC*/

int
main( void )
{
  static LALStatus stat;

  INT4 i,k; /*dummy indices*/
  EphemerisData *edat = NULL;

  char eEphFileBad[] = "earth47.dat";
  char eEphFile[] = "earth98.dat";
  char sEphFile[] = "sun98.dat";


  REAL8 alpha,delta;  /* RA and DEC (radians) in
			 ICRS realization of J2000 coords.*/

#if 0 /* Parallax is not treated yet. */
  REAL8 dInv; /* 1/(Dist. to Source), in units 1/sec */
#endif

  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));

#define DEBUG 1 /*rem non-zero is TRUE */
#if (DEBUG)

/* Checking response if data files not present */

  (*edat).ephiles.earthEphemeris = eEphFileBad;
  (*edat).ephiles.sunEphemeris = sEphFile;
  LALInitBarycenter(&stat, edat);

  if ( stat.statusCode != LALINITBARYCENTERH_EOPEN)
    {
      printf( "Got error code %d and message '%s', but expected error code %d\n",
          stat.statusCode, stat.statusDescription, LALINITBARYCENTERH_EOPEN);
      return LALBARYCENTERTESTC_EOPEN;
    }

/* Checking response if data files somehow corrupted --to be fixed!

  (*edat).ephiles.earthEphemeris = "earth98.dat";
  (*edat).ephiles.sunEphemeris = "sun98_corrupt.dat";
  LALInitBarycenter(&stat, edat);

      if ( stat.statusCode != LALINITBARYCENTERH_EEPHFILE
        || strcmp(stat.statusDescription, LALINITBARYCENTERH_MSGEEPHFILE) )
    {
      printf( "Got error code %d and message %s\n",
          stat.statusCode, stat.statusDescription );
      printf( "Expected error code %d and message %s\n",
           LALINITBARYCENTERH_EEPHFILE, LALINITBARYCENTERH_MSGEEPHFILE);
      return LALBARYCENTERTESTC_EEPHFILE;
    }
*/
#endif

/*Now inputting kosher ephemeris. files and leap sec, to illustrate
  proper usage. The real, serious TEST of the code is a script written
  by Rejean Dupuis comparing LALBarycenter to TEMPO for thousands
  of source positions and times. */

  (*edat).ephiles.earthEphemeris = eEphFile;
  (*edat).ephiles.sunEphemeris = sEphFile;

  LALInitBarycenter(&stat, edat);
  if ( stat.statusCode ) {
    XLALPrintError ("LALInitBarycenter() failed with code %d\n", stat.statusCode);
    return XLAL_EFAILED;
  }


  /* ===== now test equivalence of new XLALInitBarycenter() function ========== */
  EphemerisData *edat_xlal;
  if ( ( edat_xlal = XLALInitBarycenter ( eEphFile, sEphFile )) == NULL ) {
    XLALPrintError ("Something failed in XLALInitBarycenter(), errno =%d\n", xlalErrno );
    return XLAL_EFAILED;
  }
  if ( compare_ephemeris ( edat, edat_xlal ) != XLAL_SUCCESS ) {
    XLALPrintError ("Equivalence test failed between XLALInitEphemeris() and LALInitEphemeris()\n" );
    return XLAL_EFAILED;
  }
  XLALDestroyEphemerisData ( edat_xlal );

  /* ========================================================================== */


 /* The routines using LALBarycenter package, the code above, leading
    up LALInitBarycenter call, should be near top of main. The idea is
    that ephemeris data is read into RAM once, at the beginning.

    NOTE that the only part of the piece of the LALDetector structure
    baryinput.site that has to be filled in by the driver code is
    the 3-vector: baryinput.site.location[] .

    NOTE that the driver code that calls LALInitBarycenter must
    LALFree(edat->ephemE) and LALFree(edat->ephemS).
    The driver code that calls LALBarycenter must LALFree(edat).
 */


  { /*Now getting coords. for detector. Cached options are:
      LALDetectorIndexLHODIFF, LALDetectorIndexLLODIFF,
      LALDetectorIndexVIRGODIFF, LALDetectorIndexGEO600DIFF,
      LALDetectorIndexTAMA300DIFF,LALDetectorIndexCIT40DIFF */

  LALDetector cachedDetector;
  cachedDetector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  baryinput.site.location[0]=cachedDetector.location[0]/LAL_C_SI;
  baryinput.site.location[1]=cachedDetector.location[1]/LAL_C_SI;
  baryinput.site.location[2]=cachedDetector.location[2]/LAL_C_SI;
  }

#if (DEBUG)
/* Checking error messages when the timestamp is not within the
   1-yr ephemeris files
*/
    tGPS.gpsSeconds = t1998+5.e7;
    tGPS.gpsNanoSeconds = 0;
    LALBarycenterEarth(&stat, &earth, &tGPS, edat);
      if ( stat.statusCode != LALBARYCENTERH_EOUTOFRANGEE
        || strcmp(stat.statusDescription, LALBARYCENTERH_MSGEOUTOFRANGEE) )
    {
      printf( "Got error code %d and message %s\n",
          stat.statusCode, stat.statusDescription );
      printf( "Expected error code %d and message %s\n",
           LALBARYCENTERH_EOUTOFRANGEE, LALBARYCENTERH_MSGEOUTOFRANGEE);
      return LALBARYCENTERTESTC_EOUTOFRANGEE;
    }

/* next try calling for bad choice of RA,DEC (e.g., something
sensible in degrees, but radians)*/

      tGPS.gpsSeconds = t1998+3600;
      tGPS.gpsNanoSeconds = 0;

    LALBarycenterEarth(&stat, &earth, &tGPS, edat);


    baryinput.alpha= 120.e0;
    baryinput.delta=60.e0;
    baryinput.dInv=0.e0;

    LALBarycenter(&stat, &emit, &baryinput, &earth);
      if ( stat.statusCode != LALBARYCENTERH_EBADSOURCEPOS
        || strcmp(stat.statusDescription,LALBARYCENTERH_MSGEBADSOURCEPOS) )
    {
      printf( "Got error code %d and message %s\n",
          stat.statusCode, stat.statusDescription );
      printf( "Expected error code %d and message %s\n",
           LALBARYCENTERH_EBADSOURCEPOS, LALBARYCENTERH_MSGEBADSOURCEPOS);
      return LALBARYCENTERTESTC_EBADSOURCEPOS;
    }

#endif
/* Now running program w/o errors, to illustrate proper use. */

/*First: outer loop over pulse arrival times; LALBarycenterEarth
    called ONCE per arrival time */

  for (i=0;i < 10; i++){

    /*GPS time(sec) =  tGPS.gpsSeconds + 1.e-9*tGPS.gpsNanoSeconds  */

    tGPS.gpsSeconds = t1998;
    tGPS.gpsSeconds +=i*3600*50;
    tGPS.gpsNanoSeconds = 0;

    LALBarycenterEarth(&stat, &earth, &tGPS, edat);
    if ( stat.statusCode ) {
      XLALPrintError ("LALBarycenterEarth() failed with code %d\n", stat.statusCode);
      return XLAL_EFAILED;
    }


/*Next: inner loop over different sky positions, for each arrival time;
     LALBarycenter called ONCE per sky position (or ONCE per detector) */

    for (k=0;k<3;k++){

      alpha=(LAL_PI/12.0)*(14.e0 + 51.e0/60.e0 +
			   +38.56024702e0/3.6e3) + LAL_PI*k/10.e0;
      delta=(LAL_PI/180.e0)*(12.e0+ 19.e0/60.e0
			     +59.1434800e0/3.6e3);

      baryinput.alpha = alpha;
      baryinput.delta = delta;
      baryinput.dInv = 0.e0;

      baryinput.tgps.gpsSeconds = tGPS.gpsSeconds;
      baryinput.tgps.gpsNanoSeconds = tGPS.gpsNanoSeconds;

      LALBarycenter(&stat, &emit, &baryinput, &earth);
      if ( stat.statusCode ) {
        XLALPrintError ("LALBarycenter() failed with code %d\n", stat.statusCode);
        return XLAL_EFAILED;
      }

#if 0
      printf("%d %d %d %25.17e %25.17e\n", k,
	     tGPS.gpsSeconds,  tGPS.gpsNanoSeconds,
	     (emit.deltaT + tGPS.gpsSeconds + tGPS.gpsNanoSeconds*1.e-9),
             emit.tDot);

      printf("%d %d %25.17e\n",
	     emit.te.gpsSeconds, emit.te.gpsNanoSeconds, emit.deltaT);

      printf("%25.17e %25.17e %25.17e\n",
	     emit.rDetector[0],emit.rDetector[1],emit.rDetector[2]);

      printf("%25.17e %25.17e %25.17e\n",
	     emit.vDetector[0],emit.vDetector[1],emit.vDetector[2]);
#endif

    }
  }
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
  LALCheckMemoryLeaks();
  return 0;
}


/** Function to test equivalence between two loaded ephemeris-data structs.
 * Note: we compare everything *except* the ephiles fields, which are actually useless.
 */
int
compare_ephemeris ( const EphemerisData *edat1, const EphemerisData *edat2 )
{
  const char *fn = __func__;

  if ( !edat1 || !edat2 ) {
    XLALPrintError ("%s: invalid NULL input edat1=%p, edat2=%p\n", fn, edat1, edat2 );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  if ( edat1->nentriesE != edat2->nentriesE ) {
    XLALPrintError ("%s: different nentriesE (%d != %d)\n", fn, edat1->nentriesE, edat2->nentriesE );
    XLAL_ERROR ( fn, XLAL_EFAILED );
  }
  if ( edat1->nentriesS != edat2->nentriesS ) {
    XLALPrintError ("%s: different nentriesS (%d != %d)\n", fn, edat1->nentriesS, edat2->nentriesS );
    XLAL_ERROR ( fn, XLAL_EFAILED );
  }
  if ( edat1->dtEtable != edat2->dtEtable ) {
    XLALPrintError ("%s: different dtEtable (%g != %g)\n", fn, edat1->dtEtable, edat2->dtEtable );
    XLAL_ERROR ( fn, XLAL_EFAILED );
  }
  if ( edat1->dtStable != edat2->dtStable ) {
    XLALPrintError ("%s: different dtStable (%g != %g)\n", fn, edat1->dtStable, edat2->dtStable );
    XLAL_ERROR ( fn, XLAL_EFAILED );
  }

  /* compare earth ephemeris data */
  if ( !edat1->ephemE || !edat2->ephemE ) {
    XLALPrintError ("%s: invalid NULL ephemE pointer edat1 (%p), edat2 (%p)\n", fn, edat1->ephemE, edat2->ephemE );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  INT4 i;
  for ( i=0; i < edat1->nentriesE; i ++ )
    {
      if ( edat1->ephemE[i].gps != edat2->ephemE[i].gps ||

           edat1->ephemE[i].pos[0] != edat2->ephemE[i].pos[0] ||
           edat1->ephemE[i].pos[1] != edat2->ephemE[i].pos[1] ||
           edat1->ephemE[i].pos[2] != edat2->ephemE[i].pos[2] ||

           edat1->ephemE[i].vel[0] != edat2->ephemE[i].vel[0] ||
           edat1->ephemE[i].vel[1] != edat2->ephemE[i].vel[1] ||
           edat1->ephemE[i].vel[2] != edat2->ephemE[i].vel[2] ||

           edat1->ephemE[i].acc[0] != edat2->ephemE[i].acc[0] ||
           edat1->ephemE[i].acc[1] != edat2->ephemE[i].acc[1] ||
           edat1->ephemE[i].acc[2] != edat2->ephemE[i].acc[2]
           )
        {
          XLALPrintError ("%s: Inconsistent earth-entry %d:\n", fn, i );
          XLALPrintError ("    edat1 = %g, (%g, %g, %g), (%g, %g, %g), (%g, %g, %g)\n",
                          edat1->ephemE[i].gps,
                          edat1->ephemE[i].pos[0], edat1->ephemE[i].pos[1], edat1->ephemE[i].pos[2],
                          edat1->ephemE[i].vel[0], edat1->ephemE[i].vel[1], edat1->ephemE[i].vel[2],
                          edat1->ephemE[i].acc[0], edat1->ephemE[i].acc[1], edat1->ephemE[i].acc[2]
                          );
          XLALPrintError ("    edat2 = %g, (%g, %g, %g), (%g, %g, %g), (%g, %g, %g)\n",
                          edat2->ephemE[i].gps,
                          edat2->ephemE[i].pos[0], edat2->ephemE[i].pos[1], edat2->ephemE[i].pos[2],
                          edat2->ephemE[i].vel[0], edat2->ephemE[i].vel[1], edat2->ephemE[i].vel[2],
                          edat2->ephemE[i].acc[0], edat2->ephemE[i].acc[1], edat2->ephemE[i].acc[2]
                          );
          XLAL_ERROR ( fn, XLAL_EFAILED );
        } /* if difference in data-set i */

    } /* for i < nentriesE */

  /* compare sun ephemeris data */
  if ( !edat1->ephemS || !edat2->ephemS ) {
    XLALPrintError ("%s: invalid NULL ephemS pointer edat1 (%p), edat2 (%p)\n", fn, edat1->ephemS, edat2->ephemS );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  for ( i=0; i < edat1->nentriesS; i ++ )
    {
      if ( edat1->ephemS[i].gps != edat2->ephemS[i].gps ||

           edat1->ephemS[i].pos[0] != edat2->ephemS[i].pos[0] ||
           edat1->ephemS[i].pos[1] != edat2->ephemS[i].pos[1] ||
           edat1->ephemS[i].pos[2] != edat2->ephemS[i].pos[2] ||

           edat1->ephemS[i].vel[0] != edat2->ephemS[i].vel[0] ||
           edat1->ephemS[i].vel[1] != edat2->ephemS[i].vel[1] ||
           edat1->ephemS[i].vel[2] != edat2->ephemS[i].vel[2] ||

           edat1->ephemS[i].acc[0] != edat2->ephemS[i].acc[0] ||
           edat1->ephemS[i].acc[1] != edat2->ephemS[i].acc[1] ||
           edat1->ephemS[i].acc[2] != edat2->ephemS[i].acc[2]
           )
        {
          XLALPrintError ("%s: Inconsistent sun-entry %d:\n", fn, i );
          XLALPrintError ("    edat1 = %g, (%g, %g, %g), (%g, %g, %g), (%g, %g, %g)\n",
                          edat1->ephemS[i].gps,
                          edat1->ephemS[i].pos[0], edat1->ephemS[i].pos[1], edat1->ephemS[i].pos[2],
                          edat1->ephemS[i].vel[0], edat1->ephemS[i].vel[1], edat1->ephemS[i].vel[2],
                          edat1->ephemS[i].acc[0], edat1->ephemS[i].acc[1], edat1->ephemS[i].acc[2]
                          );
          XLALPrintError ("    edat2 = %g, (%g, %g, %g), (%g, %g, %g), (%g, %g, %g)\n",
                          edat2->ephemS[i].gps,
                          edat2->ephemS[i].pos[0], edat2->ephemS[i].pos[1], edat2->ephemS[i].pos[2],
                          edat2->ephemS[i].vel[0], edat2->ephemS[i].vel[1], edat2->ephemS[i].vel[2],
                          edat2->ephemS[i].acc[0], edat2->ephemS[i].acc[1], edat2->ephemS[i].acc[2]
                          );
          XLAL_ERROR ( fn, XLAL_EFAILED );
        } /* if difference in data-set i */

    } /* for i < nentriesE */

  /* everything seems fine */
  return XLAL_SUCCESS;

} /* compare_ephemeris() */

/** \endcond */
