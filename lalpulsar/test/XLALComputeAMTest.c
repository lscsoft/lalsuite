/*
 * Copyright (C) 2010 Reinhard Prix
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

#include <config.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <math.h>
#include <sys/times.h>

#include <gsl/gsl_math.h>

#include <lal/ComputeFstat.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>

/** \author Reinhard Prix, John Whelan
 * \file
 * \ingroup LALComputeAM_h
 *
 * \brief Test for XLALComputeAMCoeffs() and XLALComputeMultiAMCoeffs() by
 * comparison with the equivalent LAL functions LALNewGetAMCoeffs() and LALGetMultiAMCoeffs().
 *
 * Note, we run a comparison only for the 2-IFO multiAM functions XLALComputeMultiAMCoeffs()
 * comparing it to LALGetMultiAMCoeffs() [combined with XLALWeightMultiAMCoeffs()],
 * as this excercises the 1-IFO functions as well.
 *
 * Sky-location is picked at random each time, which allows a minimal
 * Monte-Carlo validation by simply running this script repeatedly.
 *
 */

extern char *optarg;

static const LALStatus empty_status;
static const AMCoeffsParams empty_AMCoeffsParams;
static const AMCoeffs empty_AMCoeffs;

extern int lalDebugLevel;

/* ----- internal prototypes ---------- */
int XLALCompareMultiAMCoeffs ( MultiAMCoeffs *multiAM1, MultiAMCoeffs *multiAM2, REAL8 tolerance );


/* ----- function definitions ---------- */
int main(int argc, char *argv[])
{
  int opt;             /* Command-line option. */

  UINT4 numIFOs = 2;
  const CHAR *sites[2] = {"H1", "V1"};

  UINT4 startTime = 714180733;
  UINT4 duration = 180000;	/* 50 hours */
  UINT4 Tsft = 1800;		/* assume 30min SFTs */

  REAL8 tolerance = 2e-6;	/* same algorithm, should be basically identical results */

  char earthEphem[] = DATADIR "earth00-19-DE405.dat.gz";
  char sunEphem[]   = DATADIR "sun00-19-DE405.dat.gz";

  UINT4 numChecks = 1; /* Number of times to check */

  /* read user input */
  lalDebugLevel = 0;

  while ((opt = getopt( argc, argv, "n:qv:" )) != -1) {
    switch (opt) {
    case 'v': /* set lalDebugLevel */
      lalDebugLevel = atoi( optarg );
      break;
    case 'n': /* number of times to check */
      numChecks = atoi( optarg );
      break;
    }
  }

  /* init random-generator */
  struct tms buf;
  UINT4 seed = times(&buf);
  srand ( seed );
  XLALPrintInfo ("%s: seed used = %d\n", __func__, seed );

  /* ----- init ephemeris ----- */
  EphemerisData *edat;
  if ( (edat = XLALInitBarycenter ( earthEphem, sunEphem )) == NULL ) {
    XLALPrintError ("%s: XLALInitBarycenter() failed with xlalErrno = %d\n", __func__, xlalErrno );
    return XLAL_EFAILED;
  }

  /* ----- init detector info ---------- */
  UINT4 X;
  MultiLALDetector *multiDet;
  if ( (multiDet = XLALCreateMultiLALDetector ( numIFOs )) == NULL ) {
    XLALPrintError ("%s: XLALCreateMultiLALDetector(%d) failed with errno=%d\n", __func__, numIFOs, xlalErrno );
    return XLAL_EFAILED;
  }

  for (X=0; X < numIFOs; X ++ )
    {
      LALDetector *site;
      if ( (site = XLALGetSiteInfo ( sites[X] )) == NULL ) {
        XLALPrintError ("%s: Failed to get site-info for detector '%s'\n", __func__, sites[X] );
        return XLAL_EFAILED;
      }
      multiDet->data[X] = (*site); 	/* copy! */
      XLALFree ( site );
    }

  /* ----- init multi-IFO timestamps vector ---------- */
  UINT4 numSteps = (UINT4) ceil ( duration / Tsft );
  MultiLIGOTimeGPSVector * multiTS;
  if ( (multiTS = XLALCalloc ( 1, sizeof(*multiTS))) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  multiTS->length = numIFOs;
  if ( (multiTS->data = XLALCalloc (numIFOs, sizeof(*multiTS->data))) == NULL ) {
    XLAL_ERROR ( XLAL_ENOMEM );
  }
  for ( X=0; X < numIFOs; X ++ )
    {
      if ( (multiTS->data[X] = XLALCreateTimestampVector (numSteps)) == NULL ) {
        XLALPrintError ("%s: XLALCreateTimestampVector(%d) failed.\n", __func__, numSteps );
        return XLAL_EFAILED;
      }
      multiTS->data[X]->deltaT = Tsft;

      UINT4 i;
      for ( i=0; i < numSteps; i ++ )
        {
          UINT4 ti = startTime + i * Tsft;
          multiTS->data[X]->data[i].gpsSeconds = ti;
          multiTS->data[X]->data[i].gpsNanoSeconds = 0;
        } /* for i < numSteps */

    } /* for X < numIFOs */

  /* ---------- compute multi-detector states -------------------- */
  MultiDetectorStateSeries *multiDetStates;
  if ( (multiDetStates = XLALGetMultiDetectorStates ( multiTS, multiDet, edat, 0.5 * Tsft )) == NULL ) {
    XLALPrintError ( "%s: XLALGetMultiDetectorStates() failed.\n", __func__ );
    return XLAL_EFAILED;
  }
  XLALDestroyMultiLALDetector ( multiDet );
  XLALDestroyMultiTimestamps ( multiTS );
  XLALFree(edat->ephemE);
  XLALFree(edat->ephemS);
  XLALFree ( edat );

  /* ========== MAIN LOOP: N-trials of comparisons XLAL <--> LAL multiAM functions ========== */
  while ( numChecks-- )
    {
      LALStatus status = empty_status;

      /* ----- pick skyposition at random ----- */
      SkyPosition skypos;
      skypos.longitude = LAL_TWOPI * (1.0 * rand() / ( RAND_MAX + 1.0 ) );  /* uniform in [0, 2pi) */
      skypos.latitude = LAL_PI_2 - acos ( 1 - 2.0 * rand()/RAND_MAX );	/* sin(delta) uniform in [-1,1] */
      skypos.system = COORDINATESYSTEM_EQUATORIAL;

      MultiNoiseWeights *weights = NULL;	/* for now we only deal with unit-weights case */

      /* ----- compute multiAM using LAL function ----- */
      MultiAMCoeffs *multiAM_LAL  = NULL;
      LALGetMultiAMCoeffs ( &status, &multiAM_LAL, multiDetStates, skypos );
      if ( status.statusCode ) {
        XLALPrintError ("%s: LALGetMultiAMCoeffs() failed with statusCode = %d : %s\n", __func__, status.statusCode, status.statusDescription );
        return XLAL_EFAILED;
      }
      if ( XLALWeightMultiAMCoeffs ( multiAM_LAL, weights ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: XLALWeightMultiAMCoeffs() failed with xlalErrno = %d\n", __func__, xlalErrno );
        return XLAL_EFAILED;
      }


      /* ----- compute multiAM using XLAL function ----- */
      MultiAMCoeffs *multiAM_XLAL;
      if ( ( multiAM_XLAL = XLALComputeMultiAMCoeffs ( multiDetStates, weights, skypos )) == NULL ) {
        XLALPrintError ("%s: XLALComputeMultiAMCoeffs() failed with xlalErrno = %d\n", __func__, xlalErrno );
        return XLAL_EFAILED;
      }

      /* now run comparison */
      if ( XLALCompareMultiAMCoeffs ( multiAM_XLAL, multiAM_LAL, tolerance ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: comparison between multiAM_XLAL and multiAM_LAL failed.\n", __func__ );
        return XLAL_EFAILED;
      }

      /* free memory created inside this loop */
      XLALDestroyMultiAMCoeffs ( multiAM_LAL );
      XLALDestroyMultiAMCoeffs ( multiAM_XLAL );

    } /* for numChecks */

  /* we're done: free memory */
  XLALDestroyMultiDetectorStateSeries ( multiDetStates );

  LALCheckMemoryLeaks();

  printf ("OK. All tests passed successfully\n\n");

  return 0;	/* OK */

} /* main() */

/** Comparison function for two multiAM vectors, return success or failure for given tolerance.
 *
 * we compare avg() and max of |a1_i - a2_i|^2 and |b1_i - b2_i|^2 respectively,
 * and error in |A1 - A2|, |B1 - B2|, |C1 - C2|.
 * These numbers are typically ~ O(1), so we simply compare these absolute errors to the tolerance.
 *
 */
int
XLALCompareMultiAMCoeffs ( MultiAMCoeffs *multiAM1, MultiAMCoeffs *multiAM2, REAL8 tolerance )
{
  /* check input */
  if ( !multiAM1 || !multiAM2 || tolerance <= 0 ) {
    XLALPrintError ("%s: invalid NULL input or non-positive tolerance\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  UINT4 numDet = multiAM1->length;
  if ( numDet != multiAM2->length ) {
    XLALPrintError ("%s: number of detectors differ multiAM1 = %d, multiAM2 = %d\n", __func__, multiAM1->length, multiAM2->length );
    XLAL_ERROR ( XLAL_EFAILED );
  }

  UINT4 X;
  REAL8 maxerr_ab = 0, avgerr_ab = 0;
  UINT4 numTerms = 0;
  for ( X=0; X < numDet; X ++ )
    {
      UINT4 numSteps = multiAM1->data[X]->a->length;
      if ( numSteps != multiAM2->data[X]->a->length ) {
        XLALPrintError ("%s: number of timesteps differ multiAM1[%d]->a = %d, multiAM2[%d]->a = %d\n",__func__, X, multiAM1->data[X]->a->length, X, multiAM2->data[X]->a->length );
        XLAL_ERROR ( XLAL_EFAILED );
      }
      if ( numSteps != multiAM1->data[X]->b->length || numSteps != multiAM2->data[X]->b->length) {
        XLALPrintError ("%s: number of timesteps differ multiAM1[%d]->b = %d, multiAM2[%d]->b = %d\n",__func__, X, multiAM1->data[X]->b->length, X, multiAM2->data[X]->b->length );
        XLAL_ERROR ( XLAL_EFAILED );
      }

      UINT4 i;
      for ( i=0; i < numSteps; i ++ )
        {
          REAL8 err_a = fabs ( multiAM1->data[X]->a->data[i] - multiAM2->data[X]->a->data[i] );
          REAL8 err_b = fabs ( multiAM1->data[X]->b->data[i] - multiAM2->data[X]->b->data[i] );

          if ( err_a > maxerr_ab )
            maxerr_ab = err_a;
          if ( err_b > maxerr_ab )
            maxerr_ab = err_b;

          avgerr_ab += err_a + err_b;
          numTerms += 1;

        } /* for i < numSteps */

    } /* for X < numDet */

  avgerr_ab /= (2.0 * numTerms);

  /* now compute absolute maximal error in AntennaPattern matrix terms */
  AntennaPatternMatrix *Mmunu1 = &multiAM1->Mmunu;
  AntennaPatternMatrix *Mmunu2 = &multiAM2->Mmunu;

  REAL8 maxxerr_Ad = 0, maxxerr_Bd = 0, maxxerr_Cd = 0, maxxerr_Dd = 0;
  REAL8 err;

  err = fabs ( Mmunu1->Ad - Mmunu2->Ad );
  if ( err > maxxerr_Ad )
    maxxerr_Ad = err;
  err = fabs ( Mmunu1->Bd - Mmunu2->Bd );
  if ( err > maxxerr_Bd )
    maxxerr_Bd = err;
  err = fabs ( Mmunu1->Cd - Mmunu2->Cd );
  if ( err > maxxerr_Cd )
    maxxerr_Cd = err;
  err = fabs ( Mmunu1->Dd - Mmunu2->Dd );
  if ( err > maxxerr_Dd )
    maxxerr_Dd = err;

  UINT4 failed = 0;
  /* special treatment for Sinv_Tsft: independent of AM-functions, should agree to within numerics */
  double eps = 1e-15;
  if ( gsl_fcmp ( Mmunu1->Sinv_Tsft, Mmunu2->Sinv_Tsft, eps ) ) {
    XLALPrintError ("%s: Sinv_Tsft differs by more than %g relative error\n", __func__, eps );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: Sinv_Tsft 1 = %g, 2 = %g, 1-2 = %g\n", __func__, Mmunu1->Sinv_Tsft, Mmunu2->Sinv_Tsft, Mmunu1->Sinv_Tsft - Mmunu2->Sinv_Tsft );

  /* ----- compare matrix elements A,B,C -------------------- */
  /* Mmunu = {A,B,C} are sums of N terms a^2, b^2, a*b, each with small numerical errors
   * we therefore need to relax the tolerance from a,b,c by a factor of sqrt(N):
   */
  REAL8 tolerance_Mmunu = tolerance * sqrt ( 1.0 * numTerms );
  XLALPrintError ("%s: tolerances tol(a,b,c)=%g, tol(Mmunu) = %g\n", __func__, tolerance, tolerance_Mmunu );

  if ( maxxerr_Ad > tolerance_Mmunu ) {
    XLALPrintError ("%s: maximal difference in Ad is %g, which exceeds the tolerance %g\n", __func__, maxxerr_Ad, tolerance_Mmunu );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: maxxerr_Ad = %g (< %g)\n", __func__, maxxerr_Ad, tolerance_Mmunu);

  if ( maxxerr_Bd > tolerance_Mmunu ) {
    XLALPrintError ("%s: maximal difference in Bd is %g, which exceeds the tolerance %g\n", __func__, maxxerr_Bd, tolerance_Mmunu );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: maxxerr_Bd = %g (< %g)\n", __func__, maxxerr_Bd, tolerance_Mmunu);

  if ( maxxerr_Cd > tolerance_Mmunu ) {
    XLALPrintError ("%s: maximal difference in Cd is %g, which exceeds the tolerance %g\n", __func__, maxxerr_Cd, tolerance_Mmunu );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: maxxerr_Cd = %g (< %g)\n", __func__, maxxerr_Cd, tolerance_Mmunu);

  /* matrix can be quite ill-conditioned, so the error on Dd is very hard to constrain.
   * in principle for random-parameters this can go arbitrarly badly, so we don't use
   * any constraints in this test to avoid spurious test-failures
   */
  XLALPrintError ("%s: maxxerr_Dd = %g (%g)\n", __func__, maxxerr_Dd, tolerance_Mmunu);

  /* ----- compare individual a,b,c errors -------------------- */
  if ( maxerr_ab > tolerance ) {
    XLALPrintError ("%s: maximal difference in {a, b} coefficients is %g, which exceeds the tolerance %g\n", __func__, maxerr_ab, tolerance );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: maxerr_ab = %g (< %g)\n", __func__, maxerr_ab, tolerance);

  if ( avgerr_ab > tolerance ) {
    XLALPrintError ("%s: average difference in {a, b} coefficients is %g, which exceeds the tolerance %g\n", __func__, avgerr_ab, tolerance );
    failed ++;
  }
  else
    XLALPrintInfo ("%s: avgerr_ab = %g (< %g)\n", __func__, avgerr_ab, tolerance);

  return failed;

} /* XLALCompareMultiAMCoeffs() */
