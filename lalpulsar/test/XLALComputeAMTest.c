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

#include <math.h>
#include <sys/times.h>

#include <gsl/gsl_math.h>

#include <lal/ComputeFstat.h>
#include <lal/LALgetopt.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>
#include <lal/SinCosLUT.h>

/**
 * \author Reinhard Prix, John Whelan
 * \file
 * \ingroup LALComputeAM_h
 *
 * \brief Test for XLALComputeAMCoeffs() and XLALComputeMultiAMCoeffs() by
 * comparison with the old LAL functions old_LALGetAMCoeffs() and old_LALGetMultiAMCoeffs().
 *
 * Note, we run a comparison only for the 2-IFO multiAM functions XLALComputeMultiAMCoeffs()
 * comparing it to old_LALGetMultiAMCoeffs() [combined with XLALWeightMultiAMCoeffs()],
 * as this excercises the 1-IFO functions as well.
 *
 * Sky-location is picked at random each time, which allows a minimal
 * Monte-Carlo validation by simply running this script repeatedly.
 *
 */

/* ----- internal prototypes ---------- */
static int XLALCompareMultiAMCoeffs ( MultiAMCoeffs *multiAM1, MultiAMCoeffs *multiAM2, REAL8 tolerance );
static void old_LALGetMultiAMCoeffs (LALStatus *, MultiAMCoeffs **multiAMcoef, const MultiDetectorStateSeries *multiDetStates, SkyPosition pos );
static void old_LALNewGetAMCoeffs(LALStatus *, AMCoeffs *coeffs, const DetectorStateSeries *DetectorStates, SkyPosition skypos);

#define LALCOMPUTEAMH_ENULL 	   7
#define LALCOMPUTEAMH_EINPUT   	   8
#define LALCOMPUTEAMH_ENONULL 	   9
#define LALCOMPUTEAMH_EMEM   	  10

#define LALCOMPUTEAMH_MSGENULL 	  "Arguments contained an unexpected null pointer"
#define LALCOMPUTEAMH_MSGEINPUT   "Invalid input"
#define LALCOMPUTEAMH_MSGENONULL  "Output pointer is non-NULL"
#define LALCOMPUTEAMH_MSGEMEM     "Out of memory. Bad."

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

  char earthEphem[] = TEST_DATA_DIR "earth00-19-DE405.dat.gz";
  char sunEphem[]   = TEST_DATA_DIR "sun00-19-DE405.dat.gz";

  UINT4 numChecks = 1; /* Number of times to check */

  /* read user input */

  while ((opt = LALgetopt( argc, argv, "n:qv:" )) != -1) {
    switch (opt) {
    case 'v': /* set lalDebugLevel */
      break;
    case 'n': /* number of times to check */
      numChecks = atoi( LALoptarg );
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
  MultiLALDetector multiDet;
  multiDet.length = numIFOs;
  for (X=0; X < numIFOs; X ++ )
    {
      LALDetector *site;
      if ( (site = XLALGetSiteInfo ( sites[X] )) == NULL ) {
        XLALPrintError ("%s: Failed to get site-info for detector '%s'\n", __func__, sites[X] );
        return XLAL_EFAILED;
      }
      multiDet.sites[X] = (*site); 	/* copy! */
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
  if ( (multiDetStates = XLALGetMultiDetectorStates ( multiTS, &multiDet, edat, 0.5 * Tsft )) == NULL ) {
    XLALPrintError ( "%s: XLALGetMultiDetectorStates() failed.\n", __func__ );
    return XLAL_EFAILED;
  }
  XLALDestroyMultiTimestamps ( multiTS );
  XLALDestroyEphemerisData ( edat );

  /* ========== MAIN LOOP: N-trials of comparisons XLAL <--> LAL multiAM functions ========== */
  while ( numChecks-- )
    {
      LALStatus XLAL_INIT_DECL(status);

      /* ----- pick skyposition at random ----- */
      SkyPosition skypos;
      skypos.longitude = LAL_TWOPI * (1.0 * rand() / ( RAND_MAX + 1.0 ) );  /* uniform in [0, 2pi) */
      skypos.latitude = LAL_PI_2 - acos ( 1 - 2.0 * rand()/RAND_MAX );	/* sin(delta) uniform in [-1,1] */
      skypos.system = COORDINATESYSTEM_EQUATORIAL;

      MultiNoiseWeights *weights = NULL;	/* for now we only deal with unit-weights case */

      /* ----- compute multiAM using old LAL function ----- */
      MultiAMCoeffs *multiAM_LAL  = NULL;
      old_LALGetMultiAMCoeffs ( &status, &multiAM_LAL, multiDetStates, skypos );
      if ( status.statusCode ) {
        XLALPrintError ("%s: old_LALGetMultiAMCoeffs() failed with statusCode = %d : %s\n", __func__, status.statusCode, status.statusDescription );
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

/*
 * Comparison function for two multiAM vectors, return success or failure for given tolerance.
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

/*
 * Multi-IFO version of old_LALGetAMCoeffs().
 * Get all antenna-pattern coefficients for all input detector-series.
 *
 * NOTE: contrary to old_LALGetAMCoeffs(), this functions *allocates* the output-vector,
 * use XLALDestroyMultiAMCoeffs() to free this.
 */
void
old_LALGetMultiAMCoeffs (LALStatus *status,			/*< [in/out] LAL status structure pointer */
		     MultiAMCoeffs **multiAMcoef,	/*< [out] AM-coefficients for all input detector-state series */
		     const MultiDetectorStateSeries *multiDetStates, /*< [in] detector-states at timestamps t_i */
		     SkyPosition skypos			/*< source sky-position [in equatorial coords!] */
		     )
{
  UINT4 X, numDetectors;
  MultiAMCoeffs *ret = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT (multiDetStates, status,LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT (multiDetStates->length, status,LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT (multiAMcoef, status,LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT ( *multiAMcoef == NULL, status,LALCOMPUTEAMH_ENONULL, LALCOMPUTEAMH_MSGENONULL);
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status, LALCOMPUTEAMH_EINPUT, LALCOMPUTEAMH_MSGEINPUT );

  numDetectors = multiDetStates->length;

  if ( ( ret = LALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, LALCOMPUTEAMH_EMEM, LALCOMPUTEAMH_MSGEMEM);
  }
  ret->length = numDetectors;
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, LALCOMPUTEAMH_EMEM, LALCOMPUTEAMH_MSGEMEM);
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      AMCoeffs *amcoeX = NULL;
      UINT4 numStepsX = multiDetStates->data[X]->length;

      ret->data[X] = LALCalloc ( 1, sizeof ( *(ret->data[X]) ) );
      amcoeX = ret->data[X];
      amcoeX->a = XLALCreateREAL4Vector ( numStepsX );
      if ( (amcoeX->b = XLALCreateREAL4Vector ( numStepsX )) == NULL ) {
	LALPrintError ("\nOut of memory!\n\n");
	goto failed;
      }

      old_LALNewGetAMCoeffs (status->statusPtr, amcoeX, multiDetStates->data[X], skypos );
      if ( status->statusPtr->statusCode )
	{
	  LALPrintError ( "\nCall to old_LALNewGetAMCoeffs() has failed ... \n\n");
	  goto failed;
	}

    } /* for X < numDetectors */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiAMCoeffs ( ret );
  ABORT ( status, -1, "old_LALGetMultiAMCoeffs() failed" );

 success:
  (*multiAMcoef) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* old_LALGetMultiAMCoeffs() */

/*
 * Compute the 'amplitude coefficients' \f$a(t)\sin\zeta\f$,
 * \f$b(t)\sin\zeta\f$ as defined in \cite JKS98 for a series of
 * timestamps.
 *
 * The input consists of the DetectorState-timeseries, which contains
 * the detector-info and the LMST's corresponding to the different times.
 *
 * In order to allow re-using the output-structure AMCoeffs for subsequent
 * calls, we require the REAL4Vectors a and b to be allocated already and
 * to have the same length as the DetectoStates-timeseries.
 *
 * \note This is an alternative implementation to both LALComputeAM()
 * and old_LALGetAMCoeffs(), which uses the geometrical definition of
 * \f$a\sin\zeta\f$ and \f$b\sin\zeta\f$ as detector response
 * coefficients in a preferred polarization basis.  (It is thereby
 * more general than the JKS expressions and could be used e.g., with
 * the response tensor of a bar detector with no further modification
 * needed.)
 */
void
old_LALNewGetAMCoeffs(LALStatus *status,			/*< [in/out] LAL status structure pointer */
	       AMCoeffs *coeffs,			/*< [out] amplitude-coeffs {a(t_i), b(t_i)} */
	       const DetectorStateSeries *DetectorStates,/*< timeseries of detector states */
	       SkyPosition skypos			/*< {alpha,delta} of the source */
	       )
{
  REAL4 delta, alpha;
  REAL4 sin1delta, cos1delta;
  REAL4 sin1alpha, cos1alpha;

  REAL4 xi1, xi2;
  REAL4 eta1, eta2, eta3;
  REAL4 norm;
  UINT4 i, numSteps;

  INITSTATUS(status);

  /*---------- check input ---------- */
  ASSERT ( DetectorStates, status, LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);

  numSteps = DetectorStates->length;

  /* require the coeffients-vectors to be allocated and consistent with timestamps */
  ASSERT ( coeffs, status, LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT ( coeffs->a && coeffs->b, status, LALCOMPUTEAMH_ENULL, LALCOMPUTEAMH_MSGENULL);
  ASSERT ( (coeffs->a->length == numSteps) && (coeffs->b->length == numSteps), status,
	   LALCOMPUTEAMH_EINPUT,  LALCOMPUTEAMH_MSGEINPUT);

  /* require sky-pos to be in equatorial coordinates */
  ASSERT ( skypos.system == COORDINATESYSTEM_EQUATORIAL, status,
	   SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );

  /*---------- We write components of xi and eta vectors in SSB-fixed coords */
  alpha = skypos.longitude;
  delta = skypos.latitude;

  if( XLALSinCosLUT (&sin1delta, &cos1delta, delta ) != XLAL_SUCCESS )
    ABORT( status->statusPtr, LAL_EXLAL, "XLALSinCosLUT (&sin1delta, &cos1delta, delta ) failed" );

  if( XLALSinCosLUT (&sin1alpha, &cos1alpha, alpha ) != XLAL_SUCCESS )
    ABORT( status->statusPtr, LAL_EXLAL, "XLALSinCosLUT (&sin1alpha, &cos1alpha, alpha ) failed" );

  // see Eq.(17) in CFSv2 notes (version v3):
  // https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=1665&version=3
  xi1 =   sin1alpha;
  xi2 =  -cos1alpha;
  eta1 = -sin1delta * cos1alpha;
  eta2 = -sin1delta * sin1alpha;
  eta3 = cos1delta;

  /*---------- Compute the a(t_i) and b(t_i) ---------- */
  coeffs->A = 0;
  coeffs->B = 0;
  coeffs->C = 0;
  coeffs->D = 0;
  for ( i=0; i < numSteps; i++ )
    {
      REAL4 ai, bi;

      SymmTensor3 *d = &(DetectorStates->data[i].detT);

      ai =    d->d11 * ( xi1 * xi1 - eta1 * eta1 )
	+ 2 * d->d12 * ( xi1*xi2 - eta1*eta2 )
	- 2 * d->d13 *             eta1 * eta3
	+     d->d22 * ( xi2*xi2 - eta2*eta2 )
	- 2 * d->d23 *             eta2 * eta3
	-     d->d33 *             eta3*eta3;

      bi =    d->d11 * 2 * xi1 * eta1
	+ 2 * d->d12 *   ( xi1 * eta2 + xi2 * eta1 )
	+ 2 * d->d13 *     xi1 * eta3
	+     d->d22 * 2 * xi2 * eta2
	+ 2 * d->d23 *     xi2 * eta3;

      /*
      printf("xi = (%f,%f)\n",xi1,xi2);
      printf("eta = (%f,%f,%f)\n",eta1,eta2,eta3);
      printf("d = (%f %f %f\n",d->d11,d->d12,d->d13);
      printf("     %f %f %f\n",d->d12,d->d22,d->d23);
      printf("     %f %f %f)\n",d->d13,d->d23,d->d33);
      */

      coeffs->a->data[i] = ai;
      coeffs->b->data[i] = bi;

      /* sum A, B, C on the fly */
      coeffs->A += ai * ai;
      coeffs->B += bi * bi;
      coeffs->C += ai * bi;

    } /* for i < numSteps */

  /* finish calculation of A,B,C, D */
  norm = 2.0f / numSteps;
  coeffs->A *= norm;
  coeffs->B *= norm;
  coeffs->C *= norm;

  coeffs->D = coeffs->A * coeffs->B - coeffs->C * coeffs->C;

  RETURN(status);

} /* old_LALNewGetAMCoeffs() */
