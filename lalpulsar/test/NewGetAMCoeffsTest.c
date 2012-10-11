/*
 *
 * Copyright (C) 2006 John Whelan, Reinhard Prix
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

/*********************************************************************************/
/** \author John Whelan (based on code by Reinhard Prix)
 * \file
 * \brief Test for LALNewGetAMCoeffs(): compare results to older, well-tested
 * (but less efficient, and harder to understand) function LALComputeAM()
 *
 * This routine currently compares LALComputeAM(), LALGetAMCoeffs() and LALNewGetAMCoeffs() to
 * each other, by computing average and maximal relative differences between the respective
 * a's and b's.
 *
 * Detector and sky-location are picked at random each time, which allows a minimal
 * Monte-Carlo validation by simply running this script repeatedly.
 *
 *********************************************************************************/
#include <config.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <math.h>
#include <sys/times.h>

#include <lal/ComputeFstat.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>
extern char *optarg;

/** \name Error codes */
/*@{*/
#define NEWGETAMCOEFFSTEST_ENORM 	0
#define NEWGETAMCOEFFSTEST_ESUB  	1

#define NEWGETAMCOEFFSTEST_MSGENORM    "Normal exit"
#define NEWGETAMCOEFFSTEST_MSGESUB     "Subroutine failed"
/*@}*/


#define RELERROR(x, y) fabs( 2.0 * ((x) - (y)) / ( (x) + (y) ) )
#define SQ(x) ( (x) * (x) )
#define MYMAX(a,b) ( (a) > (b) ? (a) : (b) )

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, "$Id$", statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              "$Id$", (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( NEWGETAMCOEFFSTEST_ESUB, NEWGETAMCOEFFSTEST_MSGESUB,      	     \
           "Function call \"" #func "\" failed:" );                  \
    return NEWGETAMCOEFFSTEST_ESUB;                                     \
  }                                                                  \
} while (0)


static const LALStatus empty_status;
static const AMCoeffsParams empty_AMCoeffsParams;
static const AMCoeffs empty_AMCoeffs;

extern int lalDebugLevel;

/** Very simple test: pick random skyposition, compute a_i, b_i using
 *  once LALComputeAM() and once LALNewGetAMCoeffs(), and look at the errors
 *  sum_i (a_i - a_i')^2
 */
int main(int argc, char *argv[])
{
  LALStatus status = empty_status;
  int              opt;             /* Command-line option. */

  LIGOTimeGPS startTime = {714180733, 0};
  REAL8 duration = 180000;	/* 50 hours */
  REAL8 Tsft = 1800;		/* assume 30min SFTs */
  LIGOTimeGPSVector *timestamps = NULL;
  DetectorStateSeries *detStates = NULL;
  SkyPosition skypos = empty_SkyPosition;
  EphemerisData edat = empty_EphemerisData;
  BarycenterInput baryinput = empty_BarycenterInput;
  LALDetector *det = NULL;
  AMCoeffs AMold = empty_AMCoeffs, AMnew1 = empty_AMCoeffs, AMnew2 = empty_AMCoeffs;
  REAL8 alpha, delta;
  AMCoeffsParams amParams = empty_AMCoeffsParams;
  EarthState earth;
  UINT4 i;
  REAL8 maxerr01, maxerr02, maxerr12, averr01, averr02, averr12;
  REAL8 tolerance = 1e-2;	/* be generous: allow 1% error */
  struct tms buf;

  const CHAR *sites[]   = {"H1", "L1", "V2", "G1", "T1" };
  REAL8 sinzeta;	/* zeta = IFO opening angle */
  UINT4 pickedSite;
  BOOLEAN ignoreErrors = 0; /* Don't fail if tolerance exceeded */
  UINT4 numChecks = 1; /* Number of times to check */
  char earthEphem[] = "earth00-04.dat";
  char sunEphem[] = "sun00-04.dat";

  /* ----- old testing code to use 9 degree earth rotations ----- */
  /* startTime.gpsSeconds = 714275242;
  duration = 86164;
  Tsft = 2154.1; */

  lalDebugLevel = 0;

  while ((opt = getopt( argc, argv, "n:qv:" )) != -1) {
    switch (opt) {
    case 'v': /* set lalDebugLevel */
      lalDebugLevel = atoi( optarg );
      break;
    case 'q': /* don't fail if tolerance exceeded */
      ignoreErrors = 1;
      break;
    case 'n': /* number of times to check */
      numChecks = atoi( optarg );
      break;
    }
  }

  /* init random-generator */
  srand ( times(&buf) );

  /* ----- init ephemeris ----- */
  edat.ephiles.earthEphemeris = earthEphem;
  edat.ephiles.sunEphemeris = sunEphem;
  SUB ( LALInitBarycenter(&status, &edat), &status);

  /* ----- get timestamps ----- */
  SUB ( LALMakeTimestamps ( &status, &timestamps, startTime, duration, Tsft ), &status );

  /* ----- allocate memory for AM-coeffs ----- */
  AMold.a = XLALCreateREAL4Vector ( timestamps->length );
  AMold.b = XLALCreateREAL4Vector ( timestamps->length );
  AMnew1.a = XLALCreateREAL4Vector ( timestamps->length );
  AMnew1.b = XLALCreateREAL4Vector ( timestamps->length );
  AMnew2.a = XLALCreateREAL4Vector ( timestamps->length );
  AMnew2.b = XLALCreateREAL4Vector ( timestamps->length );

  while ( numChecks-- )
{

  /* ----- pick detector-site at random ----- */
  pickedSite = floor( 5 * (1.0 * rand() / (RAND_MAX + 1.0) ) );  /* int in [0,5) */

  /* NOTE: contrary to ComputeAM() and LALGetAMCoffs(), the new function LALNewGetAMCoeffs()
   * computes 'a * sinzeta' and 'b * sinzeta': for the comparison we therefore need to correct
   * for GEO's opening-angle of 94.33degrees [JKS98]: */
  if ( ! strcmp ( sites[pickedSite], "G1" ) )
    sinzeta = 0.997146;
  else
    sinzeta = 1;

  if ( ( det = XLALGetSiteInfo ( sites[pickedSite] )) == NULL )
    {
      XLALPrintError ("\nCall to XLALGetSiteInfo() has failed for site = '%s'... \n\n",
		     sites[pickedSite]);
      return NEWGETAMCOEFFSTEST_ESUB;
    }

  /* ----- pick skyposition at random ----- */
  alpha = LAL_TWOPI * (1.0 * rand() / ( RAND_MAX + 1.0 ) );  /* uniform in [0, 2pi) */
  delta = LAL_PI_2 - acos ( 1 - 2.0 * rand()/RAND_MAX );	/* sin(delta) uniform in [-1,1] */
  /* ----- old testing code to put source overhead ----- */
  /*  alpha = det->frDetector.vertexLongitudeRadians;
      delta = det->frDetector.vertexLatitudeRadians; */

  /* ===== compute AM-coeffs the 'old way': ===== */
  baryinput.site.location[0] = det->location[0]/LAL_C_SI;
  baryinput.site.location[1] = det->location[1]/LAL_C_SI;
  baryinput.site.location[2] = det->location[2]/LAL_C_SI;
  baryinput.alpha = alpha;
  baryinput.delta = delta;
  baryinput.dInv = 0.e0;

  /* amParams structure to compute a(t) and b(t) */
  amParams.das = (LALDetAndSource *)LALMalloc(sizeof(LALDetAndSource));
  amParams.das->pSource = (LALSource *)LALMalloc(sizeof(LALSource));
  amParams.baryinput = &baryinput;
  amParams.earth = &earth;
  amParams.edat = &edat;
  amParams.das->pDetector = det;
  amParams.das->pSource->equatorialCoords.longitude = alpha;
  amParams.das->pSource->equatorialCoords.latitude = delta;
  amParams.das->pSource->orientation = 0.0;
  amParams.das->pSource->equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  amParams.polAngle = 0;

  SUB (LALComputeAM ( &status, &AMold, timestamps->data, &amParams), &status);

  /* ===== compute AM-coeffs the 'new way' using LALNewGetAMCoeffs() */

  /* ----- get detector-state series ----- */
  SUB ( LALGetDetectorStates (&status, &detStates, timestamps, det, &edat, 0 ), &status );

  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  skypos.longitude = alpha;
  skypos.latitude = delta;

  /* the 'new' and the 'newer' way ... */
  SUB ( LALGetAMCoeffs ( &status, &AMnew1, detStates, skypos ), &status );	/* 'new1' */
  SUB ( LALNewGetAMCoeffs ( &status, &AMnew2, detStates, skypos ), &status );	/* 'new2' */


  /* ===== analyse relative errors ===== */
  maxerr01 = maxerr02 = maxerr12 = 0; /* errors between 0='old', 1='new1', 2='new2' */
  averr01 = averr02 = averr12 = 0;
  for ( i=0; i < timestamps->length; i ++ )
    {
      /*      printf("GPS time: %d s %d ns; GMST in radians: %f\n",
	     detStates->data[i].tGPS.gpsSeconds,
	     detStates->data[i].tGPS.gpsNanoSeconds,
	     fmod(detStates->data[i].earthState.gmstRad,LAL_TWOPI));
	     printf("Old AM coeffs: a=%f, b=%f\nNew AM coeffs: a=%f, b=%f\nNEWER AM coeffs: a=%f b=%f",
	     AMold.a->data[i], AMold.b->data[i],
	     AMnew.a->data[i], AMnew.b->data[i],
	     AMnewer.a->data[i], AMnewer.b->data[i]); */
      REAL8 thisErr;
      /* compare 0-1 */
      thisErr = sqrt( SQ ( AMold.a->data[i] -  AMnew1.a->data[i] ) / AMold.A );
      averr01 += thisErr;
      maxerr01 = MYMAX( thisErr, maxerr01 );
      thisErr = sqrt( SQ ( AMold.b->data[i] -  AMnew1.b->data[i] ) / AMold.B );
      averr01 += thisErr;
      maxerr01 = MYMAX( thisErr, maxerr01 );

      /* compare 0-2 */
      thisErr = sqrt( SQ ( AMold.a->data[i] -  AMnew2.a->data[i]/sinzeta ) / AMold.A );
      averr02 += thisErr;
      maxerr02 = MYMAX( thisErr, maxerr02 );
      thisErr = sqrt( SQ ( AMold.b->data[i] -  AMnew2.b->data[i]/sinzeta ) / AMold.B );
      averr02 += thisErr;
      maxerr02 = MYMAX( thisErr, maxerr02 );

      /* compare 1-2 */
      thisErr = sqrt( SQ ( AMnew1.a->data[i] -  AMnew2.a->data[i]/sinzeta ) / AMold.A );
      averr12 += thisErr;
      maxerr12 = MYMAX( thisErr, maxerr12 );
      thisErr = sqrt( SQ ( AMnew1.b->data[i] -  AMnew2.b->data[i]/sinzeta ) / AMold.B );
      averr12 += thisErr;
      maxerr12 = MYMAX( thisErr, maxerr12 );

    }
  averr01 /= 2.0 * timestamps->length;
  averr02 /= 2.0 * timestamps->length;
  averr12 /= 2.0 * timestamps->length;

  if ( lalDebugLevel )
    {
      printf ("Parameters: IFO = %s, skypos = [%g, %g]\n", sites[pickedSite], alpha, delta );
      printf ("Maximal relative errors: maxerr(0-1) = %g %%, maxerr(0-2) = %g %% maxerr(1-2) = %g %%\n",
	      100.0 * maxerr01, 100.0 * maxerr02, 100.0 * maxerr12 );
      printf ("Average relative errors: averr(0-1)  = %g %%, averr(0-2)  = %g %% averr(1-2)  = %g %%\n",
	      100.0 * averr01, 100.0 * averr02, 100.0 * averr12 );
    }
  else
    printf ("%d %g %g \t %g %g %g \t %g %g %g\n", pickedSite, alpha, delta, averr01, averr02, averr12, maxerr01, maxerr02, maxerr12);

  if ( (averr01 > tolerance) || (averr02 > tolerance) || (averr12 > tolerance)
       || (maxerr01 > tolerance) ||(maxerr02 > tolerance) || (maxerr12 > tolerance) )
    {
      XLALPrintError ("Maximal error-tolerance of %g %% was exceeded!\n", 100.0 * tolerance );
      if (!ignoreErrors)
	return 1;
    }

  if ( lalDebugLevel )
    printf("%d checks left\n", numChecks);

  /* ---- Clean up things that were created in this loop ---- */
  XLALDestroyDetectorStateSeries ( detStates );
  detStates = NULL;
  LALFree ( det );
  LALFree ( amParams.das->pSource );
  LALFree ( amParams.das );

}

  /* ----- free memory ----- */
  XLALDestroyTimestampVector ( timestamps );
  XLALDestroyREAL4Vector ( AMold.a );
  XLALDestroyREAL4Vector ( AMold.b );
  XLALDestroyREAL4Vector ( AMnew1.a );
  XLALDestroyREAL4Vector ( AMnew1.b );
  XLALDestroyREAL4Vector ( AMnew2.a );
  XLALDestroyREAL4Vector ( AMnew2.b );

  LALFree(edat.ephemE);
  LALFree(edat.ephemS);


  LALCheckMemoryLeaks();

  return 0;	/* OK */

} /* main() */


