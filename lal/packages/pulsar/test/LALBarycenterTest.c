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

/****************************** <lalVerbatim file="LALBarycenterTestCV">
Author: Cutler, C.
$Id$
******************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Program \texttt{LALBarycenterTest.c}}
\label{ss:LALBarycenterTest.c}

Tests the routine \verb@LALBarycenter()@.  Exercises some of the error
conditions and makes sure that they work.

\subsubsection*{Usage}
\begin{verbatim}
LALBarycenterTest
\end{verbatim}

\subsubsection*{Description}

This program demonstrates the use of \verb@LALBarycenter.c@.
The two ephemeris files specified in the \verb@EphemerisFilenames@
structure (e.g., for data taken in 1998, \verb@sun98.dat@ and \verb@earth98.dat@)
are assumed to be in the same directory as the program as
the test program.

\subsubsection*{Exit codes}
\input{LALBarycenterTestCE}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALMalloc()
LALFree()
LALBarycenterInit()
LALBarycenterEarth()
LALBarycenter()
LALCheckMemoryLeaks()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALBarycenterTestCV}}

</lalLaTeX> */

#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/DetectorSite.h>
#include <lal/Date.h>

NRCSID(LALBARYCENTERTESTC,"$Id$");

/***************************** <lalErrTable file="LALBarycenterTestCE"> */
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


/***************************** </lalErrTable> */

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
  LALLeapSecFormatAndAcc lsfas = {LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
  INT4 tmpLeap; /* need this because Date pkg defines leap seconds as
                   INT4, while EphemerisData defines it to be INT2. This won't
                   cause problems before, oh, I don't know, the Earth has been
                   destroyed in nuclear holocaust. -- dwchin 2004-02-29 */

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
  (*edat).leap = 12;
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
  (*edat).leap = 12;
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

/* Next give the number of leap secs added from start of GPS epoch to
   tgps. It's perfectly OK to instead give the number of leap
   sec from start of GPS epoch to, say, Jan. 2 in year that contains
   tgps. Currently have to specify leap by hand. This will be
   replaced by a leap sec function being written by D. Chin.
   Use: leap = 11 for 1997, leap = 12 for 1998, leap = 13 for 1999,
   leap = 13 for 2000, leap = 13 for 2001, leap = 13 for 2002.
   Yes, really: the last time it changed was end of 1998, and it's
   not changing at end of 2001.
*/

  (*edat).leap = 12;

  LALInitBarycenter(&stat, edat);
  printf("stat.statusCode = %d\n",stat.statusCode);
  REPORTSTATUS(&stat);

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

    /* addition by dwchin - 2004-02-29 */
    LALLeapSecs(&stat, &tmpLeap, &tGPS, &lsfas);
    edat->leap = (INT2)tmpLeap;

    LALBarycenterEarth(&stat, &earth, &tGPS, edat);
    REPORTSTATUS(&stat);

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
      REPORTSTATUS(&stat);

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
    }
  }
  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
  LALCheckMemoryLeaks();
  return 0;
}




