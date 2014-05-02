/*
 * Copyright (C) 2013 Reinhard Prix
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

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif


#include <lal/ComputeFstat.h>
#include <lal/LALInitBarycenter.h>
#include <lal/FindRoot.h>
#include <lal/UserInput.h>

/**
 * \author Reinhard Prix
 * \file
 * \ingroup ComputeFstat_OldDemodAPI_h
 * \brief Tests for XLALAdd[Multi]BinaryTimes()
 *
 * We simply compare the results to the old+obsolete LAL functions LALGet[Multi]Binarytimes(),
 * which have been moved here, and only serve for this comparison.
 *
 * Sky-location is picked at random each time, which allows a minimal
 * Monte-Carlo validation by simply running this script repeatedly.
 *
 */

// ----- local defines
#define EA_ACC          1E-9                    /* the timing accuracy of LALGetBinaryTimes in seconds */

// ----- macros
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

// ----- global variables
static REAL8 A,B;          /* binary time delay coefficients (need to be global so that the LAL root finding procedure can see them) */

// local types
typedef struct
{
  BOOLEAN help;		/**< Print this help/usage message */
  INT4 randSeed;	/**< allow user to specify random-number seed for reproducible noise-realizations */
} UserInput_t;

/* Type defining the orbital parameters of a binary pulsar */
typedef struct tagBinaryOrbitParams {
  LIGOTimeGPS tp;         /* time of observed periapsis passage (in SSB) */
  REAL8 argp;             /* argument of periapsis (radians) */
  REAL8 asini;            /* projected, normalized orbital semi-major axis (s) */
  REAL8 ecc;              /* orbital eccentricity */
  REAL8 period;           /* orbital period (sec) */
} BinaryOrbitParams;

/* ----- internal prototypes ---------- */
static void LALGetBinarytimes (LALStatus *, SSBtimes *tBinary, const SSBtimes *tSSB, const DetectorStateSeries *DetectorStates, const BinaryOrbitParams *binaryparams, LIGOTimeGPS refTime);
static void LALGetMultiBinarytimes (LALStatus *, MultiSSBtimes **multiBinary, const MultiSSBtimes *multiSSB, const MultiDetectorStateSeries *multiDetStates, const BinaryOrbitParams *binaryparams, LIGOTimeGPS refTime);
static void EccentricAnomoly(LALStatus *status, REAL8 *tr, REAL8 lE, void *x0);

int XLALCompareSSBtimes ( REAL8 *err_DeltaT, REAL8 *err_Tdot, const SSBtimes *t1, const SSBtimes *t2 );
int XLALCompareMultiSSBtimes ( REAL8 *err_DeltaT, REAL8 *err_Tdot, const MultiSSBtimes *m1, const MultiSSBtimes *m2 );


/* ----- function definitions ---------- */
int
main ( int argc, char *argv[] )
{
  LALStatus XLAL_INIT_DECL(status);
  UserInput_t XLAL_INIT_DECL(uvar_s);
  UserInput_t *uvar = &uvar_s;

  struct tms buf;
  uvar->randSeed = times(&buf);

  // ---------- register all our user-variable ----------
  XLALregBOOLUserStruct (  help,                'h', UVAR_HELP    , "Print this help/usage message");
  XLALregINTUserStruct (   randSeed,             's', UVAR_OPTIONAL, "Specify random-number seed for reproducible noise.");

  /* read cmdline & cfgfile  */
  XLAL_CHECK ( XLALUserVarReadAllInput ( argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( uvar->help ) {	/* if help was requested, we're done */
    exit (0);
  }

  srand ( uvar->randSeed );

  REAL8 startTimeREAL8 	= 714180733;
  REAL8 duration 	= 180000;	/* 50 hours */
  REAL8 Tsft 		= 1800;		/* assume 30min SFTs */
  char earthEphem[] 	= TEST_DATA_DIR "earth00-19-DE200.dat.gz";
  char sunEphem[]   	= TEST_DATA_DIR "sun00-19-DE200.dat.gz";

  //REAL8 tolerance = 2e-10;	/* same algorithm, should be basically identical results */

  LIGOTimeGPS startTime, refTime;
  XLALGPSSetREAL8 ( &startTime, startTimeREAL8 );
  refTime = startTime;

  // pick skyposition at random ----- */
  SkyPosition skypos;
  skypos.longitude = LAL_TWOPI * (1.0 * rand() / ( RAND_MAX + 1.0 ) );  // alpha uniform in [0, 2pi)
  skypos.latitude = LAL_PI_2 - acos ( 1 - 2.0 * rand()/RAND_MAX );	// sin(delta) uniform in [-1,1]
  skypos.system = COORDINATESYSTEM_EQUATORIAL;

  // pick binary orbital parameters:
  // somewhat inspired by Sco-X1 parameters from S2-paper (PRD76, 082001 (2007), gr-qc/0605028)
  // but with a more extreme eccentricity, and random argp
  REAL8 argp = LAL_TWOPI * (1.0 * rand() / ( RAND_MAX + 1.0 ) );	// uniform in [0, 2pi)
  BinaryOrbitParams orbit;
  XLALGPSSetREAL8 ( &orbit.tp, 731163327 ); 	// time of observed periapsis passage (in SSB)
  orbit.argp = argp;		// argument of periapsis (radians)
  orbit.asini = 1.44;           // projected, normalized orbital semi-major axis (s) */
  orbit.ecc = 1e-2;             // relatively large value, for better testing
  orbit.period = 68023;		// period (s) : about ~18.9h

  // ----- step 0: prepare test-case input for calling the BinarySSB-functions
  // setup detectors
  const char *sites[3] = { "H1", "L1", "V1" };
  UINT4 numDetectors = sizeof( sites ) / sizeof ( sites[0] );

  MultiLALDetector multiIFO;
  multiIFO.length = numDetectors;
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      LALDetector *det = XLALGetSiteInfo ( sites[X] );
      XLAL_CHECK ( det != NULL, XLAL_EFUNC, "XLALGetSiteInfo ('%s') failed for detector X=%d\n", sites[X], X );
      multiIFO.sites[X] = (*det);	 // struct copy
      XLALFree ( det );
    }

  // load ephemeris
  EphemerisData *edat = XLALInitBarycenter ( earthEphem, sunEphem );
  XLAL_CHECK ( edat != NULL, XLAL_EFUNC, "XLALInitBarycenter('%s','%s') failed\n", earthEphem, sunEphem );

  // setup multi-timeseries
  MultiLIGOTimeGPSVector *multiTS;

  XLAL_CHECK ( (multiTS = XLALCalloc ( 1, sizeof(*multiTS))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (multiTS->data = XLALCalloc (numDetectors, sizeof(*multiTS->data))) != NULL, XLAL_ENOMEM );
  multiTS->length = numDetectors;

  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      multiTS->data[X] = XLALMakeTimestamps ( startTime, duration, Tsft, 0 );
      XLAL_CHECK ( multiTS->data[X] != NULL, XLAL_EFUNC, "XLALMakeTimestamps() failed.\n");
    } /* for X < numIFOs */

  // generate detector-states
  MultiDetectorStateSeries *multiDetStates = XLALGetMultiDetectorStates ( multiTS, &multiIFO, edat, 0 );
  XLAL_CHECK ( multiDetStates != NULL, XLAL_EFUNC, "XLALGetMultiDetectorStates() failed.\n");

  // generate isolated-NS SSB times
  MultiSSBtimes *multiSSBIn = XLALGetMultiSSBtimes ( multiDetStates, skypos, refTime, SSBPREC_RELATIVISTICOPT );
  XLAL_CHECK ( multiSSBIn != NULL, XLAL_EFUNC, "XLALGetMultiSSBtimes() failed.\n");

  // ----- step 1: compute reference-result using old LALGetMultiBinarytimes()
  MultiSSBtimes *multiBinary_ref = NULL;
  LALGetMultiBinarytimes (&status, &(multiBinary_ref), multiSSBIn, multiDetStates, &orbit, refTime );
  XLAL_CHECK ( status.statusCode == 0, XLAL_EFAILED, "LALGetMultiBinarytimes() failed with status = %d : '%s'\n", status.statusCode, status.statusDescription );

  // ----- step 2: compute test-result using new XLALAddMultiBinaryTimes()
  MultiSSBtimes *multiBinary_test = NULL;
  PulsarDopplerParams doppler;
  memset(&doppler, 0, sizeof(doppler));
  doppler.tp = orbit.tp;
  doppler.argp = orbit.argp;
  doppler.asini = orbit.asini;
  doppler.ecc = orbit.ecc;
  doppler.period = orbit.period;
  XLAL_CHECK ( XLALAddMultiBinaryTimes ( &multiBinary_test, multiSSBIn, &doppler ) == XLAL_SUCCESS, XLAL_EFUNC );

  // ----- step 3: compare results
  REAL8 err_DeltaT, err_Tdot;
  REAL8 tolerance = 1e-10;
  int ret = XLALCompareMultiSSBtimes ( &err_DeltaT, &err_Tdot, multiBinary_ref, multiBinary_test );
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALCompareMultiSSBtimes() failed.\n");

  XLALPrintWarning ( "INFO: err(DeltaT) = %g, err(Tdot) = %g\n", err_DeltaT, err_Tdot );

  XLAL_CHECK ( err_DeltaT < tolerance, XLAL_ETOL, "error(DeltaT) = %g exceeds tolerance of %g\n", err_DeltaT, tolerance );
  XLAL_CHECK ( err_Tdot   < tolerance, XLAL_ETOL, "error(Tdot) = %g exceeds tolerance of %g\n", err_Tdot, tolerance );

  // ---- step 4: clean-up memory
  XLALDestroyUserVars();
  XLALDestroyEphemerisData ( edat );
  XLALDestroyMultiSSBtimes ( multiBinary_test );
  XLALDestroyMultiSSBtimes ( multiBinary_ref );
  XLALDestroyMultiSSBtimes ( multiSSBIn );
  XLALDestroyMultiTimestamps ( multiTS );
  XLALDestroyMultiDetectorStateSeries ( multiDetStates );

  // check for memory-leaks
  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()

/**
 * Compare 2 MultiSSBtimes vectors and return the maximal difference in
 * DeltaT in s, and in Tdot (dimensionless)
 */
int
XLALCompareMultiSSBtimes ( REAL8 *err_DeltaT, REAL8 *err_Tdot, const MultiSSBtimes *m1, const MultiSSBtimes *m2 )
{
  XLAL_CHECK ( err_DeltaT != NULL, XLAL_EINVAL );
  XLAL_CHECK ( err_Tdot != NULL, XLAL_EINVAL );
  XLAL_CHECK ( m1 != NULL, XLAL_EINVAL );
  XLAL_CHECK ( m2 != NULL, XLAL_EINVAL );

  XLAL_CHECK ( m1->length == m2->length, XLAL_EINVAL, "Number of detectors differ %d != %d\n", m1->length, m2->length );

  REAL8 max_DeltaT = 0;
  REAL8 max_Tdot = 0;

  UINT4 numDet = m1->length;
  for ( UINT4 X = 0; X < numDet; X ++ )
    {
      REAL8 DeltaT_X, Tdot_X;
      XLAL_CHECK ( XLALCompareSSBtimes ( &DeltaT_X, &Tdot_X, m1->data[X], m2->data[X] ) == XLAL_SUCCESS, XLAL_EFUNC );
      max_DeltaT = fmax ( max_DeltaT, DeltaT_X );
      max_Tdot   = fmax ( max_Tdot, Tdot_X );
    } // for X < numDet

  (*err_DeltaT) = max_DeltaT;
  (*err_Tdot)   = max_Tdot;

  return XLAL_SUCCESS;

} // XLALCompareMultiSSBtimes()

/**
 * Compare 2 SSBtimes vectors and return the maximal difference in
 * DeltaT in s, and in Tdot (dimensionless)
 */
int
XLALCompareSSBtimes ( REAL8 *err_DeltaT, REAL8 *err_Tdot, const SSBtimes *t1, const SSBtimes *t2 )
{
  XLAL_CHECK ( err_DeltaT != NULL, XLAL_EINVAL );
  XLAL_CHECK ( err_Tdot != NULL, XLAL_EINVAL );
  XLAL_CHECK ( t1 != NULL, XLAL_EINVAL );
  XLAL_CHECK ( t2 != NULL, XLAL_EINVAL );
  XLAL_CHECK ( (t1->DeltaT != NULL) && (t1->Tdot != NULL), XLAL_EINVAL );
  XLAL_CHECK ( XLALGPSDiff ( &(t1->refTime), &(t2->refTime) ) == 0, XLAL_ETOL,
               "Different reference-times: %f != %f\n", XLALGPSGetREAL8(&(t1->refTime)), XLALGPSGetREAL8(&(t2->refTime)));

  UINT4 numSteps = t1->DeltaT->length;
  XLAL_CHECK ( t1->Tdot->length == numSteps, XLAL_EINVAL );
  XLAL_CHECK ( t2->DeltaT->length == numSteps, XLAL_EINVAL );
  XLAL_CHECK ( t2->Tdot->length == numSteps, XLAL_EINVAL );

  REAL8 max_DeltaT = 0;
  REAL8 max_Tdot = 0;

  for ( UINT4 i = 0; i < numSteps; i ++ )
    {
      REAL8 DeltaT_i = fabs ( t1->DeltaT->data[i] - t2->DeltaT->data[i] );
      REAL8 Tdot_i   = fabs ( t1->Tdot->data[i]   - t2->Tdot->data[i] );

      max_DeltaT = fmax ( max_DeltaT, DeltaT_i );
      max_Tdot   = fmax ( max_Tdot,   Tdot_i );

    } // for i < numSteps

  (*err_DeltaT) = max_DeltaT;
  (*err_Tdot)   = max_Tdot;

  return XLAL_SUCCESS;

} // XLALCompareSSBtimes()



// ---------- obsolete LAL functions LALGet[Multi]Binarytimes() kept here for comparison purposes

#define COMPUTEFSTATC_ENULL 		1
#define COMPUTEFSTATC_ENONULL 		2
#define COMPUTEFSTATC_EINPUT   		3
#define COMPUTEFSTATC_EMEM   		4
#define COMPUTEFSTATC_EXLAL		5
#define COMPUTEFSTATC_EIEEE		6

#define COMPUTEFSTATC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPUTEFSTATC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPUTEFSTATC_MSGEINPUT   	"Invalid input"
#define COMPUTEFSTATC_MSGEMEM   	"Out of memory. Bad."
#define COMPUTEFSTATC_MSGEXLAL		"XLAL function call failed"
#define COMPUTEFSTATC_MSGEIEEE		"Floating point failure"

/**
 * For a given OrbitalParams, calculate the time-differences
 * \f$\Delta T_\alpha\equiv T(t_\alpha) - T_0\f$, and their
 * derivatives \f$Tdot_\alpha \equiv d T / d t (t_\alpha)\f$.
 *
 * \note The return-vectors \a DeltaT and \a Tdot must be allocated already
 * and have the same length as the input time-series \a DetStates.
 *
 */
void
LALGetBinarytimes (LALStatus *status,				/**< pointer to LALStatus structure */
		   SSBtimes *tBinary,				/**< [out] DeltaT_alpha = T(t_alpha) - T_0; and Tdot(t_alpha) */
		   const SSBtimes *tSSB,			/**< [in] DeltaT_alpha = T(t_alpha) - T_0; and Tdot(t_alpha) */
		   const DetectorStateSeries *DetectorStates,	/**< [in] detector-states at timestamps t_i */
		   const BinaryOrbitParams *binaryparams,	/**< [in] source binary orbit parameters */
		   LIGOTimeGPS refTime				/**< SSB reference-time T_0 of pulsar-parameters */
		   )
{
  UINT4 numSteps, i;
  REAL8 refTimeREAL8;
  REAL8 Porb;           /* binary orbital period */
  REAL8 asini;          /* the projected orbital semimajor axis */
  REAL8 e,ome    ;      /* the eccentricity, one minus eccentricity */
  REAL8 sinw,cosw;      /* the sin and cos of the argument of periapsis */
  REAL8 tSSB_now;       /* the SSB time at the midpoint of each SFT in REAL8 form */
  REAL8 fracorb;        /* the fraction of orbits completed since current SSB time */
  REAL8 E;              /* the eccentric anomoly */
  DFindRootIn input;    /* the input structure for the root finding procedure */
  REAL8 acc;            /* the accuracy in radians of the eccentric anomoly computation */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (DetectorStates, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  numSteps = DetectorStates->length;		/* number of timestamps */

  ASSERT (tSSB, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (tSSB->DeltaT, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (tSSB->Tdot, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (tBinary, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (tBinary->DeltaT, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (tBinary->Tdot, status, COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  ASSERT (tSSB->DeltaT->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT (tSSB->Tdot->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT (tBinary->DeltaT->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);
  ASSERT (tBinary->Tdot->length == numSteps, status, COMPUTEFSTATC_EINPUT, COMPUTEFSTATC_MSGEINPUT);

  /* convenience variables */
  Porb = binaryparams->period;
  e = binaryparams->ecc;
  asini = binaryparams->asini;
  sinw = sin(binaryparams->argp);
  cosw = cos(binaryparams->argp);
  ome = 1.0 - e;
  refTimeREAL8 = GPS2REAL8(refTime);

  /* compute p, q and r coeeficients */
  A  = (LAL_TWOPI/Porb)*cosw*asini*sqrt(1.0-e*e) - e;
  B  = (LAL_TWOPI/Porb)*sinw*asini;
  REAL8 tp = GPS2REAL8(binaryparams->tp);
  REAL8 Tp = tp  + asini*sinw*ome;

  /* Calculate the required accuracy for the root finding procedure in the main loop */
  acc = LAL_TWOPI*(REAL8)EA_ACC/Porb;   /* EA_ACC is defined above and represents the required timing precision in seconds (roughly) */

  /* loop over the SFTs */
  for (i=0; i < numSteps; i++ )
    {

      /* define SSB time for the current SFT midpoint */
      tSSB_now = refTimeREAL8 + (tSSB->DeltaT->data[i]);

      /* define fractional orbit in SSB frame since periapsis (enforce result 0->1) */
      /* the result of fmod uses the dividend sign hence the second procedure */
      REAL8 temp = fmod((tSSB_now - Tp),Porb)/(REAL8)Porb;
      fracorb = temp - (REAL8)floor(temp);
      REAL8 x0 = fracorb * LAL_TWOPI;

      /* compute eccentric anomaly using a root finding procedure */
      input.function = EccentricAnomoly;     /* This is the name of the function we must solve to find E */
      input.xmin = 0.0;                      /* We know that E will be found between 0 and 2PI */
      input.xmax = LAL_TWOPI;
      input.xacc = acc;                      /* The accuracy of the root finding procedure */

      /* expand domain until a root is bracketed */
      LALDBracketRoot(status->statusPtr,&input,&x0);

      /* bisect domain to find eccentric anomoly E corresponding to the SSB time of the midpoint of this SFT */
      LALDBisectionFindRoot(status->statusPtr,&E,&input,&x0);

      /* use our value of E to compute the additional binary time delay */
      tBinary->DeltaT->data[i] = tSSB->DeltaT->data[i] - ( asini*sinw*(cos(E)-e) + asini*cosw*sqrt(1.0-e*e)*sin(E) );

      /* combine with Tdot (dtSSB_by_dtdet) -> dtbin_by_dtdet */
      tBinary->Tdot->data[i] = tSSB->Tdot->data[i] * ( (1.0 - e*cos(E))/(1.0 + A*cos(E) - B*sin(E)) );

    } /* for i < numSteps */


  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetBinarytimes() */

/**
 * For a given set of binary parameters we solve the following function for
 * the eccentric anomoly E
 */
static void EccentricAnomoly(LALStatus *status,
			     REAL8 *tr,
			     REAL8 lE,
			     void *x0
			     )
{
  INITSTATUS(status);
  ASSERT(x0,status, 1, "Null pointer");

  /* this is the function relating the observed time since periapse in the SSB to the true eccentric anomoly E */
  (*tr) = - (*(REAL8 *)x0) + (lE + A*sin(lE) + B*(cos(lE) - 1.0));

  RETURN(status);
}

/**
 * Multi-IFO version of LALGetBinarytimes().
 * Get all binary-timings for all input detector-series.
 *
 */
void
LALGetMultiBinarytimes (LALStatus *status,				/**< pointer to LALStatus structure */
			MultiSSBtimes **multiBinary,			/**< [out] SSB-timings for all input detector-state series */
			const MultiSSBtimes *multiSSB,			/**< [in] SSB-timings for all input detector-state series */
			const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
			const BinaryOrbitParams *binaryparams,		/**< [in] source binary orbit parameters */
			LIGOTimeGPS refTime				/**< SSB reference-time T_0 for SSB-timing */
			)
{
  UINT4 X, numDetectors;
  MultiSSBtimes *ret = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check input */
  ASSERT (multiDetStates, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiDetStates->length, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT (multiSSB, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);
  ASSERT ( *multiBinary == NULL, status,COMPUTEFSTATC_ENONULL, COMPUTEFSTATC_MSGENONULL);
  ASSERT (multiSSB != NULL, status,COMPUTEFSTATC_ENULL, COMPUTEFSTATC_MSGENULL);

  numDetectors = multiDetStates->length;

  if ( ( ret = LALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }
  ret->length = numDetectors;
  if ( ( ret->data = LALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    LALFree ( ret );
    ABORT (status, COMPUTEFSTATC_EMEM, COMPUTEFSTATC_MSGEMEM);
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      SSBtimes *BinarytimesX = NULL;
      UINT4 numStepsX = multiDetStates->data[X]->length;

      ret->data[X] = LALCalloc ( 1, sizeof ( *(ret->data[X]) ) );
      BinarytimesX = ret->data[X];
      BinarytimesX->DeltaT = XLALCreateREAL8Vector ( numStepsX );
      if ( (BinarytimesX->Tdot = XLALCreateREAL8Vector ( numStepsX )) == NULL ) {
	XLALPrintError ("\nOut of memory!\n\n");
	goto failed;
      }
      /* printf("calling  LALGetBinarytimes\n"); */
      LALGetBinarytimes (status->statusPtr, BinarytimesX, multiSSB->data[X], multiDetStates->data[X], binaryparams, refTime);
      /* printf("finished  LALGetBinarytimes\n"); */
      if ( status->statusPtr->statusCode )
	{
	  XLALPrintError ( "\nCall to LALGetBinarytimes() has failed ... \n\n");
	  goto failed;
	}

      // NOTE! It seems the old LAL-function *did not* set the reference time of the output binary-SSB vectors
      // so we 'fix' this here retro-actively, in order to facilitate comparison with the new XLAL function:
      BinarytimesX->refTime = multiSSB->data[X]->refTime;

    } /* for X < numDet */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiSSBtimes ( ret );
  ABORT ( status, -1, "LALGetMultiBinarytimes failed" );

 success:
  (*multiBinary) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetMultiBinarytimes() */
