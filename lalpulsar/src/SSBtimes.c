/*
 *  Copyright (C) 2007 Chris Messenger
 *  Copyright (C) 2006 John T. Whelan, Badri Krishnan
 *  Copyright (C) 2005, 2006, 2007, 2010 Reinhard Prix
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

/*---------- INCLUDES ----------*/
#include <math.h>

#include <lal/SSBtimes.h>
#include <lal/AVFactories.h>

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_roots.h>

/*---------- local DEFINES ----------*/
#define EA_ACC          1E-9                    /* the timing accuracy of LALGetBinaryTimes in seconds */

/*----- Macros ----- */

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

/*---------- Global variables ----------*/

/* empty initializers  */
const SSBtimes empty_SSBtimes;
const MultiSSBtimes empty_MultiSSBtimes;

/*---------- internal prototypes ----------*/

static double gsl_E_solver ( double E, void *p );

struct E_solver_params {
  double p, q, r, fracOrb;
};

static double gsl_E_solver ( REAL8 E, void *par )
{
  struct E_solver_params *params = (struct E_solver_params*) par;
  double p = params->p;
  double q = params->q;
  double r = params->r;
  double x0 = params->fracOrb * LAL_TWOPI;

  /* this is the function relating the observed time since periapse in the SSB to the true eccentric anomoly E */
  double diff = - x0 + ( E + p * sin(E) + q * ( cos(E) - 1.0 ) + r );

  return diff;
} // gsl_E_solver()

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * For a given OrbitalParams, return a newly-allocated SSBtimes structure
 * containing the input SSBtimes with the *additional* time-differences due to
 * the binary orbital motion *added* to it.
 *
 * NOTE: the output vector 'tSSBOut' can be passed either
 *    - unallocated (where it must be (*tSSBOut)==NULL), and it gets allocated here, or
 *    - it can also contain a pre-existing vector, which must then be consistent (same vector lenghts)
 *      than the input SSB-vectors.
 * This is intended to minimize unnecessary repeated memory allocs+frees on successive binary templates.
 *
 * NOTE2: it is allowed to pass the same in- and output-vectors, i.e. (*tSSBOut) = tSSBIn,
 * in which case the input vector will be safely added to.
 */
int
XLALAddBinaryTimes ( SSBtimes **tSSBOut,			/**< [out] SSB timings tSSBIn with binary offsets added, can be NULL */
                     const SSBtimes *tSSBIn,			/**< [in] SSB timings DeltaT_alpha = T(t_alpha) - T_0; and Tdot(t_alpha) */
                     const BinaryOrbitParams *binaryparams	/**< [in] source binary orbit parameters */
                     )
{
  XLAL_CHECK ( tSSBIn != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSB'\n" );
  XLAL_CHECK ( tSSBIn->DeltaT != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSBIn->DeltaT'\n" );
  XLAL_CHECK ( tSSBIn->Tdot != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSBIn->Tdot'\n" );
  XLAL_CHECK ( binaryparams != NULL, XLAL_EINVAL, "Invalid NULL input 'binaryparams'\n");

  UINT4 numSteps = tSSBIn->DeltaT->length;		/* number of timesteps */
  XLAL_CHECK (tSSBIn->Tdot->length == numSteps, XLAL_EINVAL,
                   "Length tSSBIn->DeltaT = %d, while tSSBIn->Tdot = %d\n", numSteps, tSSBIn->Tdot->length );

  SSBtimes *binaryTimes;
  // ----- prepare output timeseries: either allocate or re-use existing
  if ( (*tSSBOut) == NULL )	// creating new output vector
    {
      XLAL_CHECK ( (binaryTimes = XLALDuplicateSSBtimes ( tSSBIn )) != NULL, XLAL_EFUNC );
    }
  else if ( (*tSSBOut) == tSSBIn )	// input==output vector
    {
      binaryTimes = (*tSSBOut);
    }
  else // input vector given, but not identical to output vector
    {
      binaryTimes = (*tSSBOut);
      // need to do a few more sanity checks
      XLAL_CHECK ( binaryTimes->DeltaT->length == numSteps, XLAL_EINVAL,
                   "Length (*tSSBOut)->DeltaT = %d, while tSSBIn->DeltaT = %d\n", binaryTimes->DeltaT->length, numSteps );
      XLAL_CHECK ( binaryTimes->Tdot->length == numSteps, XLAL_EINVAL,
                        "Length tSSBOut->Tdot = %d, while tSSBIn->Tdot = %d\n", binaryTimes->Tdot->length, numSteps );
      // ... and copy the vector contents from the input SSB vector
      binaryTimes->refTime = tSSBIn->refTime;
      memcpy ( binaryTimes->DeltaT->data, tSSBIn->DeltaT->data, numSteps * sizeof(binaryTimes->DeltaT->data[0]) );
      memcpy ( binaryTimes->Tdot->data,   tSSBIn->Tdot->data,   numSteps * sizeof(binaryTimes->Tdot->data[0]) );
    } // re-using input vector

  /* ----- convenience variables */
  REAL8 Porb = binaryparams->period;		/* binary orbital period */
  REAL8 e = binaryparams->ecc;			/* the eccentricity */
  REAL8 ome = 1.0 - e;				/* 1 minus eccentricity */
  REAL8 asini = binaryparams->asini;		/* the projected orbital semimajor axis */
  REAL8 sinw = sin(binaryparams->argp);		/* the sin and cos of the argument of periapsis */
  REAL8 cosw = cos(binaryparams->argp);

  REAL8 refTimeREAL8 = XLALGPSGetREAL8 ( &tSSBIn->refTime );

  /* compute (global) p, q and r coeeficients for root-finder */
  REAL8 p = (LAL_TWOPI/Porb)*cosw*asini*sqrt(1.0-e*e);
  REAL8 q = (LAL_TWOPI/Porb)*sinw*asini;
  REAL8 r = (LAL_TWOPI/Porb)*sinw*asini*ome;

  /* Calculate the required accuracy for the root finding procedure in the main loop */
  REAL8 acc = LAL_TWOPI*(REAL8)EA_ACC/Porb;   /* EA_ACC is defined above and represents the required timing precision in seconds (roughly) */

  /* loop over the SFTs */
  for ( UINT4 i = 0; i < numSteps; i++ )
    {
      /* define SSB time for the current SFT midpoint */
      REAL8 tSSB_now = refTimeREAL8 + tSSBIn->DeltaT->data[i];

      /* define fractional orbit in SSB frame since periapsis (enforce result 0->1) */
      /* the result of fmod uses the dividend sign hence the second procedure */
      REAL8 temp = fmod((tSSB_now - XLALGPSGetREAL8 ( &binaryparams->tp ) ),Porb)/(REAL8)Porb;
      REAL8 fracorb = temp - floor(temp);	        /* the fraction of orbits completed since current SSB time */

      // ---------- FIXME: can we use GSL for this root-finding?
      REAL8 E;              /* the eccentric anomoly */
      {
        const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
        gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
        REAL8 E_lo = 0, E_hi = LAL_TWOPI;
        gsl_function F;
        struct E_solver_params pars = {p, q, r, fracorb};
        F.function = &gsl_E_solver;
        F.params = &pars;

        XLAL_CHECK ( gsl_root_fsolver_set(s, &F, E_lo, E_hi) == 0, XLAL_EFAILED );

        XLALPrintInfo ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "abstol", "err(est)");
        int max_iter = 100;
        int iter = 0;
        int status;
        do
          {
            iter++;
            status = gsl_root_fsolver_iterate(s);
            XLAL_CHECK ( (status == GSL_SUCCESS) || (status == GSL_CONTINUE), XLAL_EFAILED );
            E = gsl_root_fsolver_root(s);
            E_lo = gsl_root_fsolver_x_lower (s);
            E_hi = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval ( E_lo, E_hi, acc, 0 );

            if (status == GSL_SUCCESS) { XLALPrintInfo ("Converged:\n"); }
            XLALPrintInfo ("%5d [%.7f, %.7f] %.7f %+10.7g %10.7g\n", iter, E_lo, E_hi, E, acc, E_hi - E_lo);

          } while ( (status == GSL_CONTINUE) && (iter < max_iter) );

        XLAL_CHECK ( status == GSL_SUCCESS, XLAL_EMAXITER, "Eccentric anomaly: failed to converge within %d iterations\n", max_iter );
        gsl_root_fsolver_free(s);
      }

      /* use our value of E to compute the additional binary time delay */
      binaryTimes->DeltaT->data[i] += - ( asini*sinw*(cos(E)-e) + asini*cosw*sqrt(1.0-e*e)*sin(E) );

      /* combine with Tdot (dtSSB_by_dtdet) -> dtbin_by_dtdet */
      binaryTimes->Tdot->data[i] *= ( (1.0 - e*cos(E))/(1.0 + p*cos(E) - q*sin(E)) );

    } /* for i < numSteps */

  // pass back output SSB timings
  (*tSSBOut) = binaryTimes;

  return XLAL_SUCCESS;

} /* XLALAddBinaryTimes() */

/**
 * Multi-IFO version of XLALAddBinaryTimes().
 * For a given OrbitalParams, return a newly-allocated multiSSBtimes structure
 * containing the input SSBtimes with the *additional* time-differences due to
 * the binary orbital motion *added* to it.
 *
 * NOTE: the output vector 'multiSSBOut' can be passed either
 *    - unallocated (where it must be (*multiSSBOut)==NULL), and it gets allocated here, or
 *    - it can also contain a pre-existing vector, which must then be consistent (same vector lenghts)
 *      than the input SSB-vectors.
 * This is intended to minimize unnecessary repeated memory allocs+frees on successive binary templates.
 *
 * NOTE2: it is allowed to pass the same in- and output-vectors, i.e. (*multiSSBOut) = multiSSBIn,
 * in which case the input vector will be safely added to.
 *
 */
int
XLALAddMultiBinaryTimes ( MultiSSBtimes **multiSSBOut,		/**< [out] output SSB times */
                          const MultiSSBtimes *multiSSBIn,	/**< [in] SSB-timings for all input detector-state series */
                          const BinaryOrbitParams *binaryparams	/**< [in] source binary orbit parameters, NULL = isolated system */
                          )
{
  /* check input */
  XLAL_CHECK ( multiSSBIn != NULL, XLAL_EINVAL, "Invalid NULL input 'multiSSB'\n");
  XLAL_CHECK ( binaryparams != NULL, XLAL_EINVAL, "Invalid NULL input 'binaryparams'\n");

  UINT4 numDetectors = multiSSBIn->length;

  MultiSSBtimes *multiBinaryTimes;
  // ----- prepare output timeseries: either allocate or re-use existing
  if ( (*multiSSBOut) == NULL )	// creating new output vector
    {
      XLAL_CHECK ( (multiBinaryTimes = XLALDuplicateMultiSSBtimes ( multiSSBIn )) != NULL, XLAL_EFUNC );
    }
  else // input vector given
    {
      multiBinaryTimes = (*multiSSBOut);
      XLAL_CHECK ( multiBinaryTimes->length == numDetectors, XLAL_EINVAL,
                   "Inconsistent length (*multiSSBOut)->length = %d, while multiSSBIn->length = %d\n", (*multiSSBOut)->length, numDetectors );
      // we'll leave all other sanity-checks to XLALAddBinaryTimes() calls
    }

  // ----- simply loop over detectors for XLALAddBinaryTimes()
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      int ret = XLALAddBinaryTimes ( &(multiBinaryTimes->data[X]), multiSSBIn->data[X], binaryparams );
      XLAL_CHECK( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALAddBinaryTimes() failed for X=%d\n", X );
    } /* for X < numDet */

  // pass back result-vector
  (*multiSSBOut) = multiBinaryTimes;

  return XLAL_SUCCESS;

} /* XLALAddMultiBinaryTimes() */


/** Duplicate (ie allocate + copy) an input SSBtimes structure.
 * This can be useful for creating a copy before adding binary-orbital corrections in XLALAddBinaryTimes()
 */
SSBtimes *
XLALDuplicateSSBtimes ( const SSBtimes *tSSB )
{
  XLAL_CHECK_NULL ( tSSB != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSB'\n" );

  UINT4 len;
  SSBtimes *ret;
  ret = XLALCalloc ( 1, len = sizeof (*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( 1, %d )\n", len );

  ret->refTime = tSSB->refTime;

  if ( tSSB->DeltaT )
    {
      len = tSSB->DeltaT->length;
      ret->DeltaT = XLALCreateREAL8Vector ( len );
      XLAL_CHECK_NULL ( ret->DeltaT != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(%d) failed\n", len );
      memcpy ( ret->DeltaT->data, tSSB->DeltaT->data, len * sizeof(ret->DeltaT->data[0]) );
    }

  if ( tSSB->Tdot )
    {
      len = tSSB->Tdot->length;
      ret->Tdot = XLALCreateREAL8Vector ( len );
      XLAL_CHECK_NULL ( ret->Tdot != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(%d) failed\n", len );
      memcpy ( ret->Tdot->data, tSSB->Tdot->data, len * sizeof(ret->Tdot->data[0]) );
    }

  return ret;

} /* XLALDuplicateSSBtimes() */


/** Duplicate (ie allocate + copy) an input MultiSSBtimes structure.
 */
MultiSSBtimes *
XLALDuplicateMultiSSBtimes ( const MultiSSBtimes *multiSSB )
{
  XLAL_CHECK_NULL ( multiSSB != NULL, XLAL_EINVAL, "Invalid NULL input 'multiSSB'\n");

  UINT4 len;
  MultiSSBtimes *ret;
  ret = XLALCalloc ( 1, len=sizeof(*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( 1, %d )\n", len );

  UINT4 numDetectors = multiSSB->length;
  ret->length = numDetectors;
  ret->data = XLALCalloc ( numDetectors, len = sizeof(ret->data[0]) );
  XLAL_CHECK_NULL ( ret->data != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( %d, %d )\n", numDetectors, len );

  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      ret->data[X] = XLALDuplicateSSBtimes ( multiSSB->data[X] );
      XLAL_CHECK_NULL( ret->data[X] != NULL, XLAL_EFUNC, "XLALDuplicateSSBtimes() failed for detector X=%d\n", X );
    } // for X < numDetectors

  return ret;

} /* XLALDuplicateMultiSSBtimes() */

/** For a given DetectorStateSeries, calculate the time-differences
 *  \f$\Delta T_\alpha\equiv T(t_\alpha) - T_0\f$, and their
 *  derivatives \f$\dot{T}_\alpha \equiv d T / d t (t_\alpha)\f$.
 *
 *  \note The return-vector is allocated here
 *
 */
SSBtimes *
XLALGetSSBtimes ( const DetectorStateSeries *DetectorStates,	/**< [in] detector-states at timestamps t_i */
                  SkyPosition pos,				/**< source sky-location */
                  LIGOTimeGPS refTime,				/**< SSB reference-time T_0 of pulsar-parameters */
                  SSBprecision precision			/**< relativistic or Newtonian SSB transformation? */
                  )
{
  XLAL_CHECK_NULL ( DetectorStates != NULL, XLAL_EINVAL, "Invalid NULL input 'DetectorStates'\n" );
  XLAL_CHECK_NULL ( precision < SSBPREC_LAST, XLAL_EDOM, "Invalid value precision=%d, allowed are [0, %d]\n", precision, SSBPREC_LAST -1 );
  XLAL_CHECK_NULL ( pos.system == COORDINATESYSTEM_EQUATORIAL, XLAL_EDOM, "Only equatorial coordinate system (=%d) allowed, got %d\n", COORDINATESYSTEM_EQUATORIAL, pos.system );

  UINT4 numSteps = DetectorStates->length;		/* number of timestamps */

  // prepare output SSBtimes struct
  int len;
  SSBtimes *ret = XLALCalloc ( 1, len = sizeof(*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1,%d)\n", len );
  ret->DeltaT = XLALCreateREAL8Vector ( numSteps );
  XLAL_CHECK_NULL ( ret->DeltaT != NULL, XLAL_EFUNC, "ret->DeltaT = XLALCreateREAL8Vector(%d) failed\n", numSteps );
  ret->Tdot = XLALCreateREAL8Vector ( numSteps );
  XLAL_CHECK_NULL ( ret->Tdot != NULL, XLAL_EFUNC, "ret->Tdot = XLALCreateREAL8Vector(%d) failed\n", numSteps );

  /* convenience variables */
  REAL8 alpha = pos.longitude;
  REAL8 delta = pos.latitude;
  REAL8 refTimeREAL8 = XLALGPSGetREAL8 ( &refTime );

  BarycenterInput baryinput = empty_BarycenterInput;

  /*----- now calculate the SSB transformation in the precision required */
  switch (precision)
    {
      REAL8 vn[3];		/* unit-vector pointing to source in Cart. coord. */

    case SSBPREC_NEWTONIAN:	/* use simple vr.vn to calculate time-delay */

      /*----- get the cartesian source unit-vector */
      vn[0] = cos(alpha) * cos(delta);
      vn[1] = sin(alpha) * cos(delta);
      vn[2] = sin(delta);

      for (UINT4 i = 0; i < numSteps; i++ )
	{
	  LIGOTimeGPS *ti = &(DetectorStates->data[i].tGPS);
	  /* DeltaT_alpha */
	  ret->DeltaT->data[i]  = XLALGPSGetREAL8 ( ti );
	  ret->DeltaT->data[i] += SCALAR(vn, DetectorStates->data[i].rDetector);
	  ret->DeltaT->data[i] -= refTimeREAL8;

	  /* Tdot_alpha */
	  ret->Tdot->data[i] = 1.0 + SCALAR(vn, DetectorStates->data[i].vDetector);

	} /* for i < numSteps */

      break;

    case SSBPREC_RELATIVISTIC:	/* use LALBarycenter() to get SSB-times and derivative */

      baryinput.site = DetectorStates->detector;
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;

      baryinput.alpha = alpha;
      baryinput.delta = delta;
      baryinput.dInv = 0;

      for (UINT4 i = 0; i < numSteps; i++ )
	{
	  EmissionTime emit;
	  DetectorState *state = &(DetectorStates->data[i]);

	  baryinput.tgps = state->tGPS;

          if ( XLALBarycenter ( &emit, &baryinput, &(state->earthState) ) != XLAL_SUCCESS )
            XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALBarycenter() failed with xlalErrno = %d\n", xlalErrno );

	  ret->DeltaT->data[i] = XLALGPSGetREAL8 ( &emit.te ) - refTimeREAL8;
	  ret->Tdot->data[i] = emit.tDot;

	} /* for i < numSteps */

      break;

    case SSBPREC_RELATIVISTICOPT:	/* use optimized version XLALBarycenterOpt() */

      baryinput.site = DetectorStates->detector;
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;

      baryinput.alpha = alpha;
      baryinput.delta = delta;
      baryinput.dInv = 0;

      BarycenterBuffer *bBuffer = NULL;
      for ( UINT4 i = 0; i < numSteps; i++ )
        {
          EmissionTime emit;
          DetectorState *state = &(DetectorStates->data[i]);
          baryinput.tgps = state->tGPS;

          if ( XLALBarycenterOpt ( &emit, &baryinput, &(state->earthState), &bBuffer ) != XLAL_SUCCESS )
            XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALBarycenterOpt() failed with xlalErrno = %d\n", xlalErrno );

          ret->DeltaT->data[i] = XLALGPSGetREAL8 ( &emit.te ) - refTimeREAL8;
          ret->Tdot->data[i] = emit.tDot;

        } /* for i < numSteps */
      // free buffer memory
      XLALFree ( bBuffer );

      break;

    default:
      XLAL_ERROR_NULL (XLAL_EFAILED, "\n?? Something went wrong.. this should never be called!\n\n" );
      break;
    } /* switch precision */

  /* finally: store the reference-time used into the output-structure */
  ret->refTime = refTime;

  return ret;

} /* XLALGetSSBtimes() */

/** Multi-IFO version of LALGetSSBtimes().
 * Get all SSB-timings for all input detector-series.
 *
 * NOTE: this functions *allocates* the output-vector,
 * use XLALDestroyMultiSSBtimes() to free this.
 */
MultiSSBtimes *
XLALGetMultiSSBtimes ( const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
                       SkyPosition skypos,		/**< source sky-position [in equatorial coords!] */
                       LIGOTimeGPS refTime,		/**< SSB reference-time T_0 for SSB-timing */
                       SSBprecision precision		/**< use relativistic or Newtonian SSB timing?  */
                       )
{
  /* check input */
  XLAL_CHECK_NULL ( multiDetStates != NULL, XLAL_EINVAL, "Invalid NULL input 'multiDetStates'\n");
  XLAL_CHECK_NULL ( multiDetStates->length > 0, XLAL_EINVAL, "Invalid zero-length 'multiDetStates'\n");

  UINT4 numDetectors = multiDetStates->length;

  // prepare return struct
  int len;
  MultiSSBtimes *ret = XLALCalloc ( 1, len = sizeof( *ret ) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1,%d)\n", len );
  ret->length = numDetectors;
  ret->data = XLALCalloc ( numDetectors, len=sizeof ( *ret->data ) );
  XLAL_CHECK_NULL ( ret->data != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(%d,%d)\n", numDetectors, len );

  // loop over detectors
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      ret->data[X] = XLALGetSSBtimes ( multiDetStates->data[X], skypos, refTime, precision );
      XLAL_CHECK_NULL ( ret->data[X] != NULL, XLAL_EFUNC, "ret->data[%d] = XLALGetSSBtimes() failed with xlalErrno = %d\n", X, xlalErrno );

    } /* for X < numDet */

  return ret;

} /* XLALGetMultiSSBtimes() */

/* ===== Object creation/destruction functions ===== */

/** Destroy a MultiSSBtimes structure.
 * Note, this is "NULL-robust" in the sense that it will not crash
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs
 */
void
XLALDestroyMultiSSBtimes ( MultiSSBtimes *multiSSB )
{
  UINT4 X;
  SSBtimes *tmp;

  if ( ! multiSSB )
    return;

  if ( multiSSB->data )
    {
      for ( X=0; X < multiSSB->length; X ++ )
	{
	  if ( (tmp = multiSSB->data[X]) != NULL )
	    {
	      if ( tmp->DeltaT )
		XLALDestroyREAL8Vector ( tmp->DeltaT );
	      if ( tmp->Tdot )
		XLALDestroyREAL8Vector ( tmp->Tdot );
	      LALFree ( tmp );
	    } /* if multiSSB->data[X] */
	} /* for X < numDetectors */
      LALFree ( multiSSB->data );
    }
  LALFree ( multiSSB );

  return;

} /* XLALDestroyMultiSSBtimes() */
