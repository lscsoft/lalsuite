/*
 *
 * Copyright (C) 2006 Reinhard Prix
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
/** \author Reinhard Prix
 * \file
 * \brief Test for LALGetAMCoeffs(): compare results to older, well-tested
 * (but less efficient, and harder to understand) function  LALComputeAM()
 *
 *********************************************************************************/
#include <math.h>
#include <sys/times.h>

#include <lal/ComputeFstat.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>

NRCSID (GETAMCOEFFSTEST, "$Id$");

/** \name Error codes */
/*@{*/
#define GETAMCOEFFSTEST_ENORM 	0
#define GETAMCOEFFSTEST_ESUB  	1

#define GETAMCOEFFSTEST_MSGENORM    "Normal exit"
#define GETAMCOEFFSTEST_MSGESUB     "Subroutine failed"
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
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, GETAMCOEFFSTEST, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              GETAMCOEFFSTEST, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( GETAMCOEFFSTEST_ESUB, GETAMCOEFFSTEST_MSGESUB,      	     \
           "Function call \"" #func "\" failed:" );                  \
    return GETAMCOEFFSTEST_ESUB;                                     \
  }                                                                  \
} while (0)


static const LALStatus empty_status;
static const AMCoeffsParams empty_AMCoeffsParams;
static const AMCoeffs empty_AMCoeffs;

extern int lalDebugLevel;

/** Very simple test: pick random skyposition, compute a_i, b_i using
 *  once LALComputeAM() and once LALGetAMCoeffs(), and look at the errors
 *  sum_i (a_i - a_i')^2
 */
int main(int argc, char *argv[])
{
  LALStatus status = empty_status;

  LIGOTimeGPS startTime = {714180733, 0};
  REAL8 duration = 180000;	/* 50 hours */
  REAL8 Tsft = 1800;		/* assume 30min SFTs */
  LIGOTimeGPSVector *timestamps = NULL;
  DetectorStateSeries *detStates = NULL;
  SkyPosition skypos = empty_SkyPosition;
  EphemerisData edat = empty_EphemerisData;
  BarycenterInput baryinput = empty_BarycenterInput;
  LALDetector *det = NULL;
  AMCoeffs AMold = empty_AMCoeffs, AMnew = empty_AMCoeffs;
  REAL8 alpha, delta;
  AMCoeffsParams amParams = empty_AMCoeffsParams;
  EarthState earth;
  UINT4 i;
  REAL8 maxerr_a, maxerr_b, averr_a, averr_b;
  REAL8 tolerance = 1e-2;	/* be generous: allow 1% error */
  struct tms buf;

  const CHAR *sites[] = {"H1", "L1", "V2", "G1", "T1" };
  UINT4 pickedSite;
  char earthEphem[] = "earth00-04.dat";
  char sunEphem[] = "sun00-04.dat";

  lalDebugLevel = 0;
  if ( argc == 2 && !strcmp(argv[1], "-v1") )
    lalDebugLevel = 1;


  /* init random-generator */
  srand ( times(&buf) );

  /* ----- init ephemeris ----- */
  edat.ephiles.earthEphemeris = earthEphem;
  edat.ephiles.sunEphemeris = sunEphem;
  edat.leap = 0;
  SUB ( LALInitBarycenter(&status, &edat), &status);

  /* ----- get timestamps ----- */
  SUB ( LALMakeTimestamps ( &status, &timestamps, startTime, duration, Tsft ), &status );

  /* ----- allocate memory for AM-coeffs ----- */
  AMold.a = XLALCreateREAL4Vector ( timestamps->length );
  AMold.b = XLALCreateREAL4Vector ( timestamps->length );
  AMnew.a = XLALCreateREAL4Vector ( timestamps->length );
  AMnew.b = XLALCreateREAL4Vector ( timestamps->length );

  /* ----- pick detector-site at random ----- */
  pickedSite = floor( 5 * (1.0 * rand() / (RAND_MAX + 1.0) ) );  /* int in [0,5) */
  if ( ( det = XLALGetSiteInfo ( sites[pickedSite] )) == NULL )
    {
      LALPrintError ("\nCall to XLALGetSiteInfo() has failed for site = '%s'... \n\n",
		     sites[pickedSite]);
      return GETAMCOEFFSTEST_ESUB;
    }

  /* ----- pick skyposition at random ----- */
  alpha = LAL_TWOPI * (1.0 * rand() / ( RAND_MAX + 1.0 ) );  /* uniform in [0, 2pi) */
  delta = LAL_PI_2 - acos ( 1 - 2.0 * rand()/RAND_MAX );	/* sin(delta) uniform in [-1,1] */

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
  amParams.leapAcc = LALLEAPSEC_STRICT;

  SUB (LALComputeAM ( &status, &AMold, timestamps->data, &amParams), &status);

  /* ===== compute AM-coeffs the 'new way' using LALGetAMCoeffs() */

  /* ----- get detector-state series ----- */
  SUB ( LALGetDetectorStates (&status, &detStates, timestamps, det, &edat, 0 ), &status );

  skypos.system = COORDINATESYSTEM_EQUATORIAL;
  skypos.longitude = alpha;
  skypos.latitude = delta;

  SUB ( LALGetAMCoeffs ( &status, &AMnew, detStates, skypos ), &status );


  /* ===== analyse relative error ===== */
  maxerr_a = maxerr_b = averr_a = averr_b = 0;
  for ( i=0; i < timestamps->length; i ++ )
    {
      REAL8 thisErr;
      thisErr = sqrt( SQ ( AMold.a->data[i] -  AMnew.a->data[i] ) / AMold.A );
      averr_a += thisErr;
      maxerr_a = MYMAX( thisErr, maxerr_a );
      thisErr = sqrt( SQ ( AMold.b->data[i] - AMnew.b->data[i] ) / AMold.B );
      averr_b += thisErr;
      maxerr_b = MYMAX( thisErr, maxerr_b );
    }
  averr_a /= timestamps->length;
  averr_b /= timestamps->length;

  if ( lalDebugLevel )
    {
      printf ("Parameters: IFO = %s, skypos = [%g, %g]\n", sites[pickedSite], alpha, delta );
      printf ("Maximal relative errors: maxerr(a) = %g %%, maxerr(b) = %g %% \n",
	      100.0 * maxerr_a, 100.0 * maxerr_b);
      printf ("Average relative errors: averr(a)  = %g %%, averr(b)  = %g %% \n",
	      100.0 * averr_a, 100.0 * averr_b );
    }
  else
    printf ("%d %g %g %g %g %g %g \n", pickedSite, alpha, delta, averr_a, averr_b, maxerr_a, maxerr_b);

  if ( (averr_a > tolerance) || (averr_b > tolerance) || (maxerr_a > tolerance) ||(maxerr_b > tolerance))
    {
      LALPrintError ("Maximal error-tolerance of %g %% was exceeded!\n", 100.0 * tolerance );
      return 1;
    }

  /* ----- free memory ----- */
  XLALDestroyTimestampVector ( timestamps );
  XLALDestroyREAL4Vector ( AMold.a );
  XLALDestroyREAL4Vector ( AMold.b );
  XLALDestroyREAL4Vector ( AMnew.a );
  XLALDestroyREAL4Vector ( AMnew.b );
  LALFree ( det );
  XLALDestroyDetectorStateSeries ( detStates );
  LALFree ( amParams.das->pSource );
  LALFree ( amParams.das );

  LALFree(edat.ephemE);
  LALFree(edat.ephemS);


  LALCheckMemoryLeaks();

  return 0;	/* OK */

} /* main() */


