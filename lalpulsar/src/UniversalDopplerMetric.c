/*
 * Copyright (C) 2012--2015 Karl Wette
 * Copyright (C) 2008, 2009 Reinhard Prix
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
#ifdef _OPENMP
#include <omp.h>
#endif

/* gsl includes */
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_eigen.h>

#include <lal/FlatPulsarMetric.h>
#include <lal/PulsarTimes.h>
#include <lal/ComputeFstat.h>
#include <lal/XLALGSL.h>
#include <lal/Factorial.h>
#include <lal/LogPrintf.h>
#include <lal/MetricUtils.h>

#include <lal/EstimateAmplitudeParams.h>
#include <lal/UniversalDopplerMetric.h>
#include <lal/MetricUtils.h>

#include <lal/GSLHelpers.h>

/**
 * \author Reinhard Prix, Karl Wette
 * \ingroup UniversalDopplerMetric_h
 * \brief Function to compute the full F-statistic metric, including
 * antenna-pattern functions from multi-detector, as derived in \cite Prix07 .
 *
 */

/*---------- SCALING FACTORS ----------*/
/** shortcuts for integer powers */
#define POW2(a)  ( (a) * (a) )
#define POW3(a)  ( (a) * (a) * (a) )
#define POW4(a)  ( (a) * (a) * (a) * (a) )
#define POW5(a)  ( (a) * (a) * (a) * (a) * (a) )

/* These constants define an internal scale used within CW_Phi_i().
 * Distances are normalised by SCALE_R, and times are normalised by SCALE_T.
 */
#define SCALE_R (LAL_AU_SI)
#define SCALE_T (LAL_YRSID_SI)

/*---------- LOOKUP ARRAYS ----------*/

/**
 * Array of symbolic 'names' for various detector-motions
 */
static const struct {
  const DetectorMotionType type;
  const char *const name;
} DetectorMotionNames[] = {

#define DETMOTION_NAMES(orbit_type, spin_orbit_name, nospin_orbit_name) \
  {DETMOTION_SPIN | orbit_type, "spin" spin_orbit_name}, \
  {DETMOTION_SPINZ | orbit_type, "spinz" spin_orbit_name}, \
  {DETMOTION_SPINXY | orbit_type, "spinxy" spin_orbit_name}, \
  {orbit_type, nospin_orbit_name}   /* this must be the last entry in the macro */

  DETMOTION_NAMES(DETMOTION_ORBIT, "+orbit", "orbit"),
  DETMOTION_NAMES(DETMOTION_PTOLEORBIT, "+ptoleorbit", "ptoleorbit"),
  DETMOTION_NAMES(0, "", NULL)   /* last entry provides marker for end of list */

#undef DETMOTION_NAMES
};

/**
 * Array of descriptor structs for each Doppler coordinate name
 */
static const struct {
  const char *const name;	/**< coordinate name */
  const double scale;		/**< multiplicative scaling factor of the coordinate */
  const char *const help;	/**< help string explaining the coordinate's meaning and units */
} DopplerCoordinates[DOPPLERCOORD_LAST] = {

  [DOPPLERCOORD_FREQ]     = {"freq",    SCALE_T,          "Frequency [Units: Hz]."},
  [DOPPLERCOORD_F1DOT]    = {"f1dot",   POW2(SCALE_T),    "First spindown [Units: Hz/s]."},
  [DOPPLERCOORD_F2DOT]    = {"f2dot",   POW3(SCALE_T),    "Second spindown [Units: Hz/s^2]."},
  [DOPPLERCOORD_F3DOT]    = {"f3dot",   POW4(SCALE_T),    "Third spindown [Units: Hz/s^3]."},

  [DOPPLERCOORD_GC_NU0]   = {"gc_nu0",  SCALE_T,          "Global correlation frequency [Units: Hz]."},
  [DOPPLERCOORD_GC_NU1]   = {"gc_nu1",  POW2(SCALE_T),    "Global correlation first spindown [Units: Hz/s]."},
  [DOPPLERCOORD_GC_NU2]   = {"gc_nu2",  POW3(SCALE_T),    "Global correlation second spindown [Units: Hz/s^2]."},
  [DOPPLERCOORD_GC_NU3]   = {"gc_nu3",  POW4(SCALE_T),    "Global correlation third spindown [Units: Hz/s^3]."},

  [DOPPLERCOORD_ALPHA]    = {"alpha",   SCALE_R/LAL_C_SI, "Right ascension [Units: radians]."},
  [DOPPLERCOORD_DELTA]    = {"delta",   SCALE_R/LAL_C_SI, "Declination [Units: radians]."},

  [DOPPLERCOORD_N2X_EQU]  = {"n2x_equ", SCALE_R/LAL_C_SI, "X component of constrained sky position in equatorial coordinates [Units: none]."},
  [DOPPLERCOORD_N2Y_EQU]  = {"n2y_equ", SCALE_R/LAL_C_SI, "Y component of constrained sky position in equatorial coordinates [Units: none]."},

  [DOPPLERCOORD_N2X_ECL]  = {"n2x_ecl", SCALE_R/LAL_C_SI, "X component of constrained sky position in ecliptic coordinates [Units: none]."},
  [DOPPLERCOORD_N2Y_ECL]  = {"n2y_ecl", SCALE_R/LAL_C_SI, "Y component of constrained sky position in ecliptic coordinates [Units: none]."},

  [DOPPLERCOORD_N3X_EQU]  = {"n3x_equ", SCALE_R/LAL_C_SI, "X component of unconstrained super-sky position in equatorial coordinates [Units: none]."},
  [DOPPLERCOORD_N3Y_EQU]  = {"n3y_equ", SCALE_R/LAL_C_SI, "Y component of unconstrained super-sky position in equatorial coordinates [Units: none]."},
  [DOPPLERCOORD_N3Z_EQU]  = {"n3z_equ", SCALE_R/LAL_C_SI, "Z component of unconstrained super-sky position in equatorial coordinates [Units: none]."},

  [DOPPLERCOORD_N3X_ECL]  = {"n3x_ecl", SCALE_R/LAL_C_SI, "X component of unconstrained super-sky position in ecliptic coordinates [Units: none]."},
  [DOPPLERCOORD_N3Y_ECL]  = {"n3y_ecl", SCALE_R/LAL_C_SI, "Y component of unconstrained super-sky position in ecliptic coordinates [Units: none]."},
  [DOPPLERCOORD_N3Z_ECL]  = {"n3z_ecl", SCALE_R/LAL_C_SI, "Z component of unconstrained super-sky position in ecliptic coordinates [Units: none]."},

  [DOPPLERCOORD_N3SX_EQU] = {"n3sx_equ",SCALE_R/LAL_C_SI, "X spin-component of unconstrained super-sky position in equatorial coordinates [Units: none]."},
  [DOPPLERCOORD_N3SY_EQU] = {"n3sy_equ",SCALE_R/LAL_C_SI, "Y spin-component of unconstrained super-sky position in equatorial coordinates [Units: none]."},

  [DOPPLERCOORD_N3OX_ECL] = {"n3ox_ecl",SCALE_R/LAL_C_SI, "X orbit-component of unconstrained super-sky position in ecliptic coordinates [Units: none]."},
  [DOPPLERCOORD_N3OY_ECL] = {"n3oy_ecl",SCALE_R/LAL_C_SI, "Y orbit-component of unconstrained super-sky position in ecliptic coordinates [Units: none]."},
  [DOPPLERCOORD_N3OZ_ECL] = {"n3oz_ecl",SCALE_R/LAL_C_SI, "Z orbit-component of unconstrained super-sky position in ecliptic coordinates [Units: none]."},

  [DOPPLERCOORD_ASINI]    = {"asini",   1,                "Projected semimajor axis of binary orbit in small-eccentricy limit (ELL1 model) [Units: light seconds]."},
  [DOPPLERCOORD_TASC]     = {"tasc",    1,                "Time of ascension (neutron star crosses line of nodes moving away from observer) for binary orbit (ELL1 model) [Units: GPS seconds]."},
  [DOPPLERCOORD_PORB]     = {"porb",    1,                "Period of binary orbit (ELL1 model) [Units: s]."},
  [DOPPLERCOORD_KAPPA]    = {"kappa", 	1,		  "Lagrange parameter 'kappa = ecc * cos(argp)', ('ecc' = eccentricity, 'argp' = argument of periapse) of binary orbit (ELL1 model) [Units: none]."},
  [DOPPLERCOORD_ETA]      = {"eta",     1,                "Lagrange parameter 'eta = ecc * sin(argp) of binary orbit (ELL1 model) [Units: none]."}

};

/*---------- DEFINES ----------*/
#define TRUE  (1==1)
#define FALSE (0==1)

/* get the coordinate ID for a given coordinate index */
#define GET_COORD_ID(coordSys, coord) ( (coord) >= 0 ? (coordSys)->coordIDs[(coord)] : DOPPLERCOORD_NONE )

/* get the coordinate name and scale associated with a given coordinate index */
#define GET_COORD_NAME(coordSys, coord) ( (coord) >= 0 ? DopplerCoordinates[GET_COORD_ID(coordSys, coord)].name : "(none)" )
#define GET_COORD_SCALE(coordSys, coord) ( (coord) >= 0 ? DopplerCoordinates[GET_COORD_ID(coordSys, coord)].scale : 1.0 )

/** copy 3 components of Euklidean vector */
#define COPY_VECT(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define DOT_VECT(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

/** Vector product of two 3-dim vectors in cartesian coords */
#define CROSS_VECT_0(u,v) ((u)[1]*(v)[2] - (u)[2]*(v)[1])
#define CROSS_VECT_1(u,v) ((u)[2]*(v)[0] - (u)[0]*(v)[2])
#define CROSS_VECT_2(u,v) ((u)[0]*(v)[1] - (u)[1]*(v)[0])
#define CROSS_VECT(x,u,v) do { (x)[0] = CROSS_VECT_0(u,v); (x)[1] = CROSS_VECT_1(u,v); (x)[2] = CROSS_VECT_2(u,v); } while (0)

/** Operations on 3-dim vectors  */
#define ZERO_VECT(v) do{ (v)[0] = 0; (v)[1] = 0; (v)[2] = 0; } while(0)
#define NORMSQ_VECT(v) DOT_VECT(v,v)
#define NORM_VECT(v) sqrt(NORMSQ_VECT(v))
#define MULT_VECT(v,lam) do{ (v)[0] *= (lam); (v)[1] *= (lam); (v)[2] *= (lam); } while(0)
#define ADD_VECT(dst,src) do { (dst)[0] += (src)[0]; (dst)[1] += (src)[1]; (dst)[2] += (src)[2]; } while(0)
#define SUB_VECT(dst,src) do { (dst)[0] -= (src)[0]; (dst)[1] -= (src)[1]; (dst)[2] -= (src)[2]; } while(0)

/** Convert 3-D vector from equatorial into ecliptic coordinates */
#define EQU_VECT_TO_ECL(dst,src) do { \
    (dst)[0] = (src)[0]; \
    (dst)[1] = (src)[1]*LAL_COSIEARTH + (src)[2]*LAL_SINIEARTH; \
    (dst)[2] = (src)[2]*LAL_COSIEARTH - (src)[1]*LAL_SINIEARTH; \
  } while(0)

/** Convert 3-D vector from ecliptic into equatorial coordinates */
#define ECL_VECT_TO_EQU(dst,src) do { \
    (dst)[0] = (src)[0]; \
    (dst)[1] = (src)[1]*LAL_COSIEARTH - (src)[2]*LAL_SINIEARTH; \
    (dst)[2] = (src)[2]*LAL_COSIEARTH + (src)[1]*LAL_SINIEARTH; \
  } while(0)

#define SQUARE(x) ((x) * (x))

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYMAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MYMIN(a,b) ( (a) < (b) ? (a) : (b) )

#define RELERR(dx,x) ( (x) == 0 ? (dx) : ( (dx) / (x) ) )

/* highest supported spindown-order */
#define MAX_SPDNORDER 4

/** 5-point derivative formulas (eg see http://math.fullerton.edu/mathews/articles/2003NumericalDiffFormulae.pdf) */
#define DERIV5P_1(pm2,pm1,p0,pp1,pp2,h) ( ( (pm2) - 8.0 * (pm1) + 8.0 * (pp1) - (pp2)) / ( 12.0 * (h) ) )
#define DERIV5P_2(pm2,pm1,p0,pp1,pp2,h) ( (-(pm2) + 16.0 * (pm1) - 30.0 * (p0) + 16.0 * (pp1) - (pp2) ) / ( 12.0 * (h) * (h) ) )
#define DERIV5P_3(pm2,pm1,p0,pp1,pp2,h) ( (-(pm2) + 2.0 * (pm1) - 2.0 * (pp1) + (pp2) ) / ( 2.0 * (h) * (h) * (h) ) )
#define DERIV5P_4(pm2,pm1,p0,pp1,pp2,h) ( ( (pm2) - 4.0 * (pm1) + 6.0 * (p0) - 4.0 * (pp1) + (pp2) ) / ( (h) * (h) * (h) * (h) ) )

/*----- SWITCHES -----*/
/*---------- internal types ----------*/

/**
 * components of antenna-pattern function: q_l = {a(t), b(t)}
 */
typedef enum {
  AMCOMP_NONE   = -1,	/**< no antenna pattern function: (equivalent "a = 1") */
  AMCOMP_A	= 0,	/**< a(t) */
  AMCOMP_B,		/**< b(t) */
  AMCOMP_LAST
} AM_comp_t;

/** parameters for metric-integration */
typedef struct
{
  int errnum;				/**< store XLAL error of any failures within integrator */
  double epsrel;			/**< relative error tolerance for GSL integration routines */
  double epsabs;			/**< absolute error tolerance for GSL integration routines */
  double intT;				/**< length of integration time segments for phase integrals */
  DetectorMotionType detMotionType;	/**< which detector-motion to use in metric integration */
  const DopplerCoordinateSystem *coordSys;/**< Doppler coordinate system in use */
  const gsl_matrix *coordTransf;        /**< coordinate transform to apply to coordinate system */
  int coord1, coord2;			/**< coordinate indexes of the two components of the derivative-product Phi_i_Phi_j to compute*/
  int coord;				/**< coordinate index of the component for single phase-derivative Phi_i compute */
  AM_comp_t amcomp1, amcomp2;		/**< two AM components q_l q_m */
  const PulsarDopplerParams *dopplerPoint;/**< Doppler params to compute metric for */
  REAL8 startTime;			/**< GPS start time of observation */
  REAL8 refTime;			/**< GPS reference time for pulsar parameters */
  REAL8 Tspan;				/**< length of observation time in seconds */
  const LALDetector *site;		/**< detector site to compute metric for */
  const EphemerisData *edat;		/**< ephemeris data */
  vect3Dlist_t *rOrb_n;			/**< list of orbital-radius derivatives at refTime of order n = 0, 1, ... */
  BOOLEAN approxPhase;			/**< use an approximate phase-model, neglecting Roemer delay in spindown coordinates (or orders \>= 1) */
} intparams_t;


/*---------- Global variables ----------*/

/* Some local constants. */
#define rOrb_c  (LAL_AU_SI / LAL_C_SI)
#define rEarth_c  (LAL_REARTH_SI / LAL_C_SI)
#define vOrb_c  (LAL_TWOPI * LAL_AU_SI / LAL_C_SI / LAL_YRSID_SI)

/*---------- internal prototypes ----------*/
static DopplerFstatMetric* XLALComputeFmetricFromAtoms ( const FmetricAtoms_t *atoms, REAL8 cosi, REAL8 psi );
static gsl_matrix* XLALComputeFisherFromAtoms ( const FmetricAtoms_t *atoms, PulsarAmplitudeParams Amp );

static double CW_am1_am2_Phi_i_Phi_j ( double tt, void *params );
static double CW_Phi_i ( double tt, void *params );

static double XLALAverage_am1_am2_Phi_i_Phi_j ( const intparams_t *params, double *relerr_max );
static double XLALCovariance_Phi_ij ( const MultiLALDetector *multiIFO, const MultiNoiseFloor *multiNoiseFloor, const LALSegList *segList,
                                      const intparams_t *params, double *relerr_max );

static UINT4 findHighestGCSpinOrder ( const DopplerCoordinateSystem *coordSys );

/*==================== FUNCTION DEFINITIONS ====================*/


/**
 * Integrate a general quadruple product CW_am1_am2_Phi_i_Phi_j() from 0 to 1.
 * This implements the expression \f$\langle<q_1 q_2 \phi_i \phi_j\rangle\f$
 * for single-IFO average over the observation time.
 *
 * The input parameters correspond to CW_am1_am2_Phi_i_Phi_j()
 */
static double
XLALAverage_am1_am2_Phi_i_Phi_j ( const intparams_t *params, double *relerr_max )
{
  intparams_t par = (*params);	/* struct-copy, as the 'coord' field has to be changeable */
  gsl_function integrand;
  double epsrel = params->epsrel;
  double epsabs = params->epsabs;
  const size_t limit = 512;
  gsl_integration_workspace *wksp = NULL;
  int stat;

  integrand.params = (void*)&par;

  const UINT4 intN = (UINT4) ceil ( params->Tspan / params->intT );
  const REAL8 dT = 1.0 / intN;

  REAL8 res = 0;
  REAL8 abserr2 = 0;

  /* allocate workspace for adaptive integration */
  wksp = gsl_integration_workspace_alloc(limit);
  XLAL_CHECK_REAL8(wksp != NULL, XLAL_EFAULT);

  const double scale12 = GET_COORD_SCALE(par.coordSys, par.coord1) * GET_COORD_SCALE(par.coordSys, par.coord2);

  /* compute <q_1 q_2 phi_i phi_j> as an integral from tt=0 to tt=1 */
  integrand.function = &CW_am1_am2_Phi_i_Phi_j;
  for ( UINT4 n = 0; n < intN; n ++ )
    {
      REAL8 ti = 1.0 * n * dT;
      REAL8 tf = MYMIN( (n+1.0) * dT, 1.0 );
      REAL8 err_n, res_n;

      XLAL_CALLGSL ( stat = gsl_integration_qag (&integrand, ti, tf, epsabs, epsrel, limit, GSL_INTEG_GAUSS61, wksp, &res_n, &err_n) );
      /* NOTE: we don't fail if the requested level of accuracy was not reached, rather we only output a warning
       * and try to estimate the final accuracy of the integral
       */
      if ( stat != 0 )
        {
          XLALPrintWarning ( "\n%s: GSL-integration 'gsl_integration_qag()' of <am1_am2_Phi_i Phi_j> did not reach requested precision!\n", __func__ );
          XLALPrintWarning ("xlalErrno=%i, Segment n=%d, Result = %g, abserr=%g ==> relerr = %.2e > %.2e\n",
                            par.errnum, n, res_n, err_n, fabs(err_n/res_n), epsrel);
          /* XLAL_ERROR_REAL8( XLAL_EFUNC ); */
        }

      res_n *= scale12;
      err_n *= scale12;

      res += res_n;
      abserr2 += SQUARE ( err_n );

    } /* for n < intN */

  REAL8 relerr = RELERR( sqrt(abserr2), fabs(res) );
  if ( relerr_max )
    (*relerr_max) = relerr;

  gsl_integration_workspace_free(wksp);

  return res;

} /* XLALAverage_am1_am2_Phi_i_Phi_j() */


/**
 * For gsl-integration: general quadruple product between two antenna-pattern functions
 * am1, am2 in {a(t),b(t)} and two phase-derivatives phi_i * phi_j,
 * i.e. compute an expression of the form
 * \f$q_1(t) q_2(t) \phi_i(t) \phi_j(t)\f$, where \f$q_i = \{a(t), b(t)\}\f$.
 *
 * NOTE: this can be 'truncated' to any sub-expression by using
 * AMCOMP_NONE for antenna-pattern component and DOPPLERCOORD_NONE for DopplerCoordinate,
 * eg in this way this function can be used to compute \f$a^2(t), b^2(t), a(t) b(t)\f$,
 * or \f$phi_i(t) phi_j(t)\f$.
 */
static double
CW_am1_am2_Phi_i_Phi_j ( double tt, void *params )
{
  intparams_t *par = (intparams_t*) params;
  if (par->errnum) {
    return GSL_NAN;
  }

  REAL8 am1, am2, phi_i, phi_j, ret;

  am1 = am2 = 1.0;	/* default */
  phi_i = phi_j = 1.0;

  /* do we need any antenna-pattern functions in here? */
  if ( (par->amcomp1 != AMCOMP_NONE) || (par->amcomp1 != AMCOMP_NONE) )
    {
      REAL8 ttSI;
      LIGOTimeGPS ttGPS;
      SkyPosition skypos;
      REAL8 ai, bi;

      skypos.system = COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = par->dopplerPoint->Alpha;
      skypos.latitude  = par->dopplerPoint->Delta;

      ttSI = par->startTime + tt * par->Tspan;	/* current GPS time in seconds */
      XLALGPSSetREAL8( &ttGPS, ttSI );

      int errnum;
      XLAL_TRY( XLALComputeAntennaPatternCoeffs ( &ai, &bi, &skypos, &ttGPS, par->site, par->edat ), errnum );
      if ( errnum ) {
	XLALPrintError ( "%s: Call to XLALComputeAntennaPatternCoeffs() failed!\n", __func__);
	GSL_ERROR_VAL( "Failure in CW_am1_am2_Phi_i_Phi_j", GSL_EFAILED, GSL_NAN );
      }

      /* first antenna-pattern component */
      if ( par->amcomp1 == AMCOMP_A )
	am1 = ai;
      else if ( par->amcomp1 == AMCOMP_B )
	am1 = bi;
      /* second antenna-pattern component */
      if ( par->amcomp2 == AMCOMP_A )
	am2 = ai;
      else if ( par->amcomp2 == AMCOMP_B )
	am2 = bi;

    } /* if any antenna-pattern components needed */

  /* first Doppler phase derivative */
  if ( GET_COORD_ID(par->coordSys, par->coord1) != DOPPLERCOORD_NONE )
    {
      par->coord = par->coord1;
      phi_i = CW_Phi_i ( tt, params );
    }
  else
    phi_i = 1.0;

  /* second Doppler phase derivative */
  if ( GET_COORD_ID(par->coordSys, par->coord2) != DOPPLERCOORD_NONE )
    {
      par->coord = par->coord2;
      phi_j = CW_Phi_i ( tt, params );
    }
  else
    phi_j = 1.0;

  ret = am1 * am2 * phi_i * phi_j;

  LogPrintf( LOG_DETAIL,
             "amcomp1=%d amcomp2=%d coord1='%s' coord2='%s' tt=%f am1=%f am2=%f phi_i=%f phi_j=%f ret=%f\n",
             par->amcomp1, par->amcomp2,
             GET_COORD_NAME(par->coordSys, par->coord1), GET_COORD_NAME(par->coordSys, par->coord2),
             tt, am1, am2, phi_i, phi_j, ret );

  return ( ret );

} /* CW_am1_am2_Phi_i_Phi_j() */


REAL8
XLALComputePhaseDerivative ( REAL8 t,                                        ///< time 't' to compute derivative at: (effectively interpreted as SSB time if approxPhase given)
                             const PulsarDopplerParams *dopplerPoint,        ///< phase-evolution ('doppler') parameters to compute phase-derivative at
                             DopplerCoordinateID coordID,                    ///< coord-ID of coordinate wrt which to compute phase derivative
                             const EphemerisData *edat,                      ///< ephemeris data
                             const LALDetector *site,                        ///< optional detector (if NULL: use H1 as default)
                             BOOLEAN includeRoemer			     ///< whether to include Roemer-delay correction 'tSSB = t + dRoemer(t)' or not
                             )
{
  XLAL_CHECK_REAL8 ( dopplerPoint != NULL, XLAL_EINVAL );
  XLAL_CHECK_REAL8 ( edat != NULL, XLAL_EINVAL );
  XLAL_CHECK_REAL8 ( (coordID > DOPPLERCOORD_NONE) && (coordID < DOPPLERCOORD_LAST), XLAL_EINVAL );

  const DopplerCoordinateSystem coordSys = {
    .dim = 1,
    .coordIDs = { coordID },
  };

  REAL8 refTime8 = XLALGPSGetREAL8 ( &(dopplerPoint->refTime) );

 intparams_t params = {
   .detMotionType = DETMOTION_SPIN | DETMOTION_ORBIT, // doesn't matter, dummy value
   .coordSys = &coordSys,
   .dopplerPoint = dopplerPoint,
   .startTime = t,
   .refTime = refTime8,
   .Tspan = 0, // doesn't matter, dummy value
   .site = site ? site : &lalCachedDetectors[LAL_LHO_4K_DETECTOR],	// fall back to H1 if passed as NULL
   .edat = edat,
   .approxPhase = includeRoemer,	// whether to include Roemer-delay correction tSSB = t + dRoemer(t)
 };

 double scale = DopplerCoordinates[coordID].scale;

 return (REAL8) (scale * CW_Phi_i ( 0, &params ) );

}  // XLALComputePhaseDerivativeSSB()

/**
 * Partial derivative of continuous-wave (CW) phase, with respect
 * to Doppler coordinate 'i' := intparams_t->phderiv
 *
 * Time is in 'natural units' of Tspan, i.e. tt is in [0, 1] corresponding
 * to GPS-times in [startTime, startTime + Tspan ]
 *
 */
static double
CW_Phi_i ( double tt, void *params )
{
  intparams_t *par = (intparams_t*) params;
  if (par->errnum) {
    return GSL_NAN;
  }

  vect3D_t nn_equ, nn_ecl;	/* skypos unit vector */
  vect3D_t nDeriv_i;	/* derivative of sky-pos vector wrt i */

  /* positions/velocities at time tt: */
  PosVel3D_t XLAL_INIT_DECL(spin_posvel);
  PosVel3D_t XLAL_INIT_DECL(orbit_posvel);
  PosVel3D_t XLAL_INIT_DECL(posvel);

  /* positions in ecliptic plane */
  vect3D_t XLAL_INIT_DECL(ecl_orbit_pos);
  vect3D_t XLAL_INIT_DECL(ecl_pos);

  /* get skypos-vector */
  const REAL8 cosa = cos(par->dopplerPoint->Alpha);
  const REAL8 sina = sin(par->dopplerPoint->Alpha);
  const REAL8 cosd = cos(par->dopplerPoint->Delta);
  const REAL8 sind = sin(par->dopplerPoint->Delta);

  /* ... in an equatorial coordinate-frame */
  nn_equ[0] = cosd * cosa;
  nn_equ[1] = cosd * sina;
  nn_equ[2] = sind;

  /* and in an ecliptic coordinate-frame */
  EQU_VECT_TO_ECL( nn_ecl, nn_equ );

  /* get current detector position r(t) and velocity v(t) */
  REAL8 ttSI = par->startTime + tt * par->Tspan;	/* current GPS time in seconds */
  LIGOTimeGPS ttGPS;
  XLALGPSSetREAL8( &ttGPS, ttSI );
  if ( XLALDetectorPosVel ( &spin_posvel, &orbit_posvel, &ttGPS, par->site, par->edat, par->detMotionType ) != XLAL_SUCCESS ) {
    par->errnum = xlalErrno;
    XLALPrintError ( "%s: Call to XLALDetectorPosVel() failed!\n", __func__);
    return GSL_NAN;
  }

  /* XLALDetectorPosVel() returns detector positions and velocities from XLALBarycenter(),
     which are in units of LAL_C_SI, so we multiply by LAL_C_SI to get back to SI units.
     We then divide by the internal scale SCALE_R. */
  MULT_VECT(spin_posvel.pos, LAL_C_SI / SCALE_R );
  MULT_VECT(orbit_posvel.pos, LAL_C_SI / SCALE_R );
  MULT_VECT(spin_posvel.vel, LAL_C_SI / SCALE_R );
  MULT_VECT(orbit_posvel.vel, LAL_C_SI / SCALE_R );

  /* calculate detector total motion = spin motion + orbital motion */
  COPY_VECT(posvel.pos, spin_posvel.pos);
  ADD_VECT(posvel.pos, orbit_posvel.pos);
  COPY_VECT(posvel.vel, spin_posvel.vel);
  ADD_VECT(posvel.vel, orbit_posvel.vel);

  /* compute detector positions projected onto ecliptic plane */
  EQU_VECT_TO_ECL(ecl_orbit_pos, orbit_posvel.pos);
  EQU_VECT_TO_ECL(ecl_pos, posvel.pos);

  /* get frequency of Doppler point */
  const REAL8 Freq = par->dopplerPoint->fkdot[0];

  /* get time span, normalised by internal scale SCALE_R */
  const REAL8 Tspan = par->Tspan / SCALE_T;

  /* account for referenceTime != startTime, in units of SCALE_T */
  REAL8 tau0 = ( par->startTime - par->refTime ) / SCALE_T;

  /* barycentric time-delay since reference time, in units of SCALE_T */
  REAL8 tau = tau0 + ( tt * Tspan );

  /* correct for time-delay from SSB to detector (Roemer delay), neglecting relativistic effects */
  if ( !par->approxPhase )
    {
      /* SSB time-delay in internal scaled units */
      REAL8 dTRoemer = DOT_VECT(nn_equ, posvel.pos) * (SCALE_R / LAL_C_SI / SCALE_T);

      tau += dTRoemer;
    }

  /* get 'reduced' detector position of order 'n': r_n(t) [passed in as argument]
   * defined as: r_n(t) = r(t) - dot{r_orb}(tau_ref) tau - 1/2! ddot{r_orb}(tau_re) tau^2 - ....
   */
  vect3D_t rr_ord_Equ, rr_ord_Ecl;
  UINT4 i, n;
  COPY_VECT ( rr_ord_Equ, posvel.pos );		/* 0th order */
  if ( par->rOrb_n )
    {
      /* par->rOrb_n are derivatives with SI time units, so multiply tau by SCALE_T */
      const double tauSI = tau * SCALE_T;
      for (n=0; n < par->rOrb_n->length; n ++ )
        {
          /* par->rOrb_n->data[n][i] is calculated from detector positions given by XLALBarycenter(),
             which are in units of LAL_C_SI, so we multiply by LAL_C_SI to get back to SI units.
             We then divide by the internal scale SCALE_R. */
          REAL8 pre_n = LAL_FACT_INV[n] * pow(tauSI, n) * (LAL_C_SI / SCALE_R);
          for (i=0; i<3; i++)
          {
            rr_ord_Equ[i] -=  pre_n * par->rOrb_n->data[n][i];
          }
        }
    } /* if rOrb_n */
  EQU_VECT_TO_ECL( rr_ord_Ecl, rr_ord_Equ );	  /* convert into ecliptic coordinates */

  // ---------- prepare shortcuts for binary orbital parameters ----------
  // See Leaci, Prix, PRD91, 102003 (2015):  DOI:10.1103/PhysRevD.91.102003
  REAL8 orb_asini = par->dopplerPoint->asini;
  REAL8 orb_Omega = ( LAL_TWOPI / par->dopplerPoint->period );
  REAL8 orb_kappa = par->dopplerPoint->ecc * cos ( par->dopplerPoint->argp );	// Eq.(33)
  REAL8 orb_eta   = par->dopplerPoint->ecc * sin ( par->dopplerPoint->argp );	// Eq.(34)

  REAL8 orb_phase = par->dopplerPoint->argp + orb_Omega * ( ttSI - XLALGPSGetREAL8 ( &(par->dopplerPoint->tp) ) ); // Eq.(35),(36)
  REAL8 sinPsi  = sin ( orb_phase );
  REAL8 cosPsi  = cos ( orb_phase );
  REAL8 sin2Psi = sin ( 2.0 * orb_phase );
  REAL8 cos2Psi = cos ( 2.0 * orb_phase );

  /* now compute the requested (possibly linear combination of) phase derivative(s) */
  REAL8 phase_deriv = 0.0;
  for ( int coord = 0; coord < (int)par->coordSys->dim; ++coord ) {

    /* get the coefficient used to multiply phase derivative term */
    REAL8 coeff = 0.0;
    if ( par->coordTransf != NULL ) {
      coeff = gsl_matrix_get( par->coordTransf, par->coord, coord );
    } else if ( par->coord == coord ) {
      coeff = 1.0;
    }
    if ( coeff == 0.0 ) {
      continue;
    }
    coeff *= GET_COORD_SCALE(par->coordSys, coord) / GET_COORD_SCALE(par->coordSys, par->coord);

    /* compute the phase derivative term */
    REAL8 ret = 0.0;
    const DopplerCoordinateID deriv = GET_COORD_ID(par->coordSys, coord);
    switch ( deriv )
    {

    case DOPPLERCOORD_FREQ:		/**< Frequency [Units: Hz]. */
    case DOPPLERCOORD_GC_NU0:		/**< Global correlation frequency [Units: Hz]. Activates 'reduced' detector position. */
      ret = LAL_TWOPI * tau * LAL_FACT_INV[1];
      break;
    case DOPPLERCOORD_F1DOT:		/**< First spindown [Units: Hz/s]. */
    case DOPPLERCOORD_GC_NU1:		/**< Global correlation first spindown [Units: Hz/s]. Activates 'reduced' detector position. */
      ret = LAL_TWOPI * POW2(tau) * LAL_FACT_INV[2];
      break;
    case DOPPLERCOORD_F2DOT:		/**< Second spindown [Units: Hz/s^2]. */
    case DOPPLERCOORD_GC_NU2:		/**< Global correlation second spindown [Units: Hz/s^2]. Activates 'reduced' detector position. */
      ret = LAL_TWOPI * POW3(tau) * LAL_FACT_INV[3];
      break;
    case DOPPLERCOORD_F3DOT:		/**< Third spindown [Units: Hz/s^3]. */
    case DOPPLERCOORD_GC_NU3:		/**< Global correlation third spindown [Units: Hz/s^3]. Activates 'reduced' detector position. */
      ret = LAL_TWOPI * POW4(tau) * LAL_FACT_INV[4];
      break;

    case DOPPLERCOORD_ALPHA:		/**< Right ascension [Units: radians]. Uses 'reduced' detector position. */
      nDeriv_i[0] = - cosd * sina;
      nDeriv_i[1] =   cosd * cosa;
      nDeriv_i[2] =   0;
      ret = LAL_TWOPI * Freq * DOT_VECT(rr_ord_Equ, nDeriv_i);
      break;
    case DOPPLERCOORD_DELTA:		/**< Declination [Units: radians]. Uses 'reduced' detector position. */
      nDeriv_i[0] = - sind * cosa;
      nDeriv_i[1] = - sind * sina;
      nDeriv_i[2] =   cosd;
      ret = LAL_TWOPI * Freq * DOT_VECT(rr_ord_Equ, nDeriv_i);
      break;

    case DOPPLERCOORD_N2X_EQU:		/**< X component of constrained sky position in equatorial coordinates [Units: none]. Uses 'reduced' detector position. */
      ret = LAL_TWOPI * Freq * ( rr_ord_Equ[0] - (nn_equ[0]/nn_equ[2]) * rr_ord_Equ[2] );
      break;
    case DOPPLERCOORD_N2Y_EQU:		/**< Y component of constrained sky position in equatorial coordinates [Units: none]. Uses 'reduced' detector position. */
      ret = LAL_TWOPI * Freq * ( rr_ord_Equ[1] - (nn_equ[1]/nn_equ[2]) * rr_ord_Equ[2] );
      break;

    case DOPPLERCOORD_N2X_ECL:		/**< X component of constrained sky position in ecliptic coordinates [Units: none]. Uses 'reduced' detector position. */
      ret = LAL_TWOPI * Freq * ( rr_ord_Ecl[0] - (nn_ecl[0]/nn_ecl[2]) * rr_ord_Ecl[2] );
      break;
    case DOPPLERCOORD_N2Y_ECL:		/**< Y component of constrained sky position in ecliptic coordinates [Units: none]. Uses 'reduced' detector position. */
      ret = LAL_TWOPI * Freq * ( rr_ord_Ecl[1] - (nn_ecl[1]/nn_ecl[2]) * rr_ord_Ecl[2] );
      break;

    case DOPPLERCOORD_N3X_EQU:		/**< X component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * posvel.pos[0];
      break;
    case DOPPLERCOORD_N3Y_EQU:		/**< Y component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * posvel.pos[1];
      break;
    case DOPPLERCOORD_N3Z_EQU:		/**< Z component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * posvel.pos[2];
      break;

    case DOPPLERCOORD_N3X_ECL:		/**< X component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * ecl_pos[0];
      break;
    case DOPPLERCOORD_N3Y_ECL:		/**< Y component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * ecl_pos[1];
      break;
    case DOPPLERCOORD_N3Z_ECL:		/**< Z component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * ecl_pos[2];
      break;

    case DOPPLERCOORD_N3SX_EQU:	/**< X spin-component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * spin_posvel.pos[0];
      break;
    case DOPPLERCOORD_N3SY_EQU:	/**< Y spin-component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * spin_posvel.pos[1];
      break;

    case DOPPLERCOORD_N3OX_ECL:	/**< X orbit-component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * ecl_orbit_pos[0];
      break;
    case DOPPLERCOORD_N3OY_ECL:	/**< Y orbit-component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * ecl_orbit_pos[1];
      break;
    case DOPPLERCOORD_N3OZ_ECL:	/**< Z orbit-component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * ecl_orbit_pos[2];
      break;

      // ---------- binary orbital parameters ----------
      // Phase derivates taken from Eq.(39) in Leaci, Prix, PRD91, 102003 (2015):  DOI:10.1103/PhysRevD.91.102003
    case DOPPLERCOORD_ASINI: /**< Projected semimajor axis of binary orbit in small-eccentricy limit (ELL1 model) [Units: light seconds]. */
      ret = - LAL_TWOPI * Freq * ( sinPsi + 0.5 * orb_kappa * sin2Psi - 0.5 * orb_eta * cos2Psi );
      break;
    case DOPPLERCOORD_TASC: /**< Time of ascension (neutron star crosses line of nodes moving away from observer) for binary orbit (ELL1 model) [Units: GPS s]. */
      ret = LAL_TWOPI * Freq * orb_asini * orb_Omega * ( cosPsi + orb_kappa * cos2Psi + orb_eta * sin2Psi );
      break;
    case DOPPLERCOORD_PORB: /**< Period of binary orbit (ELL1 model) [Units: s]. */
      ret = Freq * orb_asini * orb_Omega * orb_phase *  ( cosPsi + orb_kappa * cos2Psi + orb_eta * sin2Psi );
      break;
    case DOPPLERCOORD_KAPPA: /**< Lagrange parameter 'kappa = ecc * cos(argp)', ('ecc' = eccentricity, 'argp' = argument of periapse) of binary orbit (ELL1 model) [Units: none] */
      ret = - LAL_PI * Freq * orb_asini * sin2Psi;
      break;
    case DOPPLERCOORD_ETA: /**< Lagrange parameter 'eta = ecc * sin(argp) of binary orbit (ELL1 model) [Units: none] */
      ret = LAL_PI * Freq * orb_asini * cos2Psi;
      break;

      // ------------------------------------------------
    default:
      par->errnum = XLAL_EINVAL;
      XLALPrintError("%s: Unknown phase-derivative type '%d'\n", __func__, deriv );
      return GSL_NAN;
      break;

    } /* switch deriv */

    phase_deriv += coeff * ret;

  }

  return phase_deriv;

} /* CW_Phi_i() */


/**
 * Given a GPS time and detector, return the current position (and velocity) of the detector.
 */
int
XLALDetectorPosVel ( PosVel3D_t *spin_posvel,		/**< [out] instantaneous sidereal position and velocity vector */
                     PosVel3D_t *orbit_posvel,		/**< [out] instantaneous orbital position and velocity vector */
                     const LIGOTimeGPS *tGPS,		/**< [in] GPS time */
                     const LALDetector *site,		/**< [in] detector info */
                     const EphemerisData *edat,		/**< [in] ephemeris data */
                     DetectorMotionType detMotionType	/**< [in] detector motion type */
                     )
{
  EarthState earth;
  BarycenterInput XLAL_INIT_DECL(baryinput);
  EmissionTime XLAL_INIT_DECL(emit);
  PosVel3D_t Det_wrt_Earth;
  PosVel3D_t PtoleOrbit;
  PosVel3D_t Spin_z, Spin_xy;
  const vect3D_t eZ = {0, -LAL_SINIEARTH, LAL_COSIEARTH};       /* ecliptic z-axis in equatorial coordinates */

  XLAL_CHECK( tGPS, XLAL_EFAULT );
  XLAL_CHECK( site, XLAL_EFAULT );
  XLAL_CHECK( edat, XLAL_EFAULT );

  XLAL_CHECK( detMotionType > 0, XLAL_EINVAL, "Invalid detector motion type '%d'", detMotionType );

  /* ----- find ephemeris-based position of Earth wrt to SSB at this moment */
  XLAL_CHECK( XLALBarycenterEarth( &earth, tGPS, edat ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* ----- find ephemeris-based position of detector wrt to SSB */
  baryinput.tgps = *tGPS;
  baryinput.site = *site;
  baryinput.site.location[0] /= LAL_C_SI; baryinput.site.location[1] /= LAL_C_SI; baryinput.site.location[2] /= LAL_C_SI;
  baryinput.alpha = 0; baryinput.delta = 0; baryinput.dInv = 0;
  XLAL_CHECK( XLALBarycenter ( &emit, &baryinput, &earth ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* ----- determine position-vector of detector wrt center of Earth */
  COPY_VECT(Det_wrt_Earth.pos, emit.rDetector);
  SUB_VECT(Det_wrt_Earth.pos, earth.posNow);
  COPY_VECT(Det_wrt_Earth.vel, emit.vDetector);
  SUB_VECT(Det_wrt_Earth.vel, earth.velNow);

  /* compute ecliptic-z projected spin motion */
  REAL8 pz = DOT_VECT ( Det_wrt_Earth.pos, eZ );
  REAL8 vz = DOT_VECT ( Det_wrt_Earth.vel, eZ );
  COPY_VECT ( Spin_z.pos, eZ );
  MULT_VECT ( Spin_z.pos, pz );
  COPY_VECT ( Spin_z.vel, eZ );
  MULT_VECT ( Spin_z.vel, vz );

  /* compute ecliptic-xy projected spin motion */
  COPY_VECT ( Spin_xy.pos, Det_wrt_Earth.pos );
  SUB_VECT ( Spin_xy.pos, Spin_z.pos );
  COPY_VECT ( Spin_xy.vel, Det_wrt_Earth.vel );
  SUB_VECT ( Spin_xy.vel, Spin_z.vel );

  /* ----- Ptolemaic special case: orbital motion on a circle */
  if ( detMotionType & DETMOTION_PTOLEORBIT ) {
    XLAL_CHECK( XLALPtolemaicPosVel ( &PtoleOrbit, tGPS ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  /* ----- return the requested type of detector motion */
  if ( spin_posvel ) {
    switch ( detMotionType & DETMOTION_MASKSPIN ) {

    case 0:   /**< No spin motion */
      ZERO_VECT(spin_posvel->pos);
      ZERO_VECT(spin_posvel->vel);
      break;

    case DETMOTION_SPIN:   /**< Full spin motion */
      COPY_VECT(spin_posvel->pos, Det_wrt_Earth.pos);
      COPY_VECT(spin_posvel->vel, Det_wrt_Earth.vel);
      break;

    case DETMOTION_SPINZ:   /**< Ecliptic-Z component of spin motion only */
      COPY_VECT ( spin_posvel->pos, Spin_z.pos );
      COPY_VECT ( spin_posvel->vel, Spin_z.vel );
      break;

    case DETMOTION_SPINXY:   /**< Ecliptic-X+Y components of spin motion only */
      COPY_VECT ( spin_posvel->pos, Spin_xy.pos );
      COPY_VECT ( spin_posvel->vel, Spin_xy.vel );
      break;

    default:
      XLAL_ERROR( XLAL_EINVAL, "%s: Illegal spin motion '%d'\n\n", __func__, detMotionType & DETMOTION_MASKSPIN );
      break;
    }
  }
  if ( orbit_posvel ) {
    switch ( detMotionType & DETMOTION_MASKORBIT ) {

    case 0:   /**< No orbital motion */
      ZERO_VECT(orbit_posvel->pos);
      ZERO_VECT(orbit_posvel->vel);
      break;

    case DETMOTION_ORBIT:   /**< Ephemeris-based orbital motion */
      COPY_VECT(orbit_posvel->pos, earth.posNow);
      COPY_VECT(orbit_posvel->vel, earth.velNow);
      break;

    case DETMOTION_PTOLEORBIT:   /**< Ptolemaic (circular) orbital motion */
      COPY_VECT(orbit_posvel->pos, PtoleOrbit.pos);
      COPY_VECT(orbit_posvel->vel, PtoleOrbit.vel);
      break;

    default:
      XLAL_ERROR( XLAL_EINVAL, "%s: Illegal orbit motion '%d'\n\n", __func__, detMotionType & DETMOTION_MASKORBIT );
      break;
    }
  }

  return XLAL_SUCCESS;

} /* XLALDetectorPosVel() */



/**
 * Compute position and velocity assuming a purely "Ptolemaic" orbital motion
 * (i.e. on a circle) around the sun, approximating Earth's orbit
 */
int
XLALPtolemaicPosVel ( PosVel3D_t *posvel,		/**< [out] instantaneous position and velocity vector */
		      const LIGOTimeGPS *tGPS		/**< [in] GPS time */
		      )
{
  PulsarTimesParamStruc XLAL_INIT_DECL(times);
  REAL8 phiOrb;   /* Earth orbital revolution angle, in radians. */
  REAL8 sinOrb, cosOrb;
  LALStatus XLAL_INIT_DECL(status);

  if ( !posvel || !tGPS ) {
    XLALPrintError ( "%s: Illegal NULL pointer passed!\n", __func__);
    XLAL_ERROR( XLAL_EINVAL );
  }

  times.epoch = (*tGPS);	/* get tAutumn */
  LALGetEarthTimes ( &status, &times );
  if ( status.statusCode ) {
    XLALPrintError ( "%s: call to LALGetEarthTimes() failed!\n\n", __func__);
    XLAL_ERROR( XLAL_EFUNC );
  }

  phiOrb = - LAL_TWOPI * times.tAutumn / LAL_YRSID_SI;
  sinOrb = sin(phiOrb);
  cosOrb = cos(phiOrb);

  /* Get instantaneous position. */
  posvel->pos[0] = rOrb_c * cosOrb;
  posvel->pos[1] = rOrb_c * sinOrb * LAL_COSIEARTH;
  posvel->pos[2]=  rOrb_c * sinOrb * LAL_SINIEARTH;

  /* Get instantaneous velocity. */
  posvel->vel[0] = -vOrb_c * sinOrb;
  posvel->vel[1] =  vOrb_c * cosOrb * LAL_COSIEARTH;
  posvel->vel[2] =  vOrb_c * cosOrb * LAL_SINIEARTH;

  return XLAL_SUCCESS;

} /* XLALPtolemaicPosVel() */



/**
 * Compute a pure phase-deriv covariance \f$[\phi_i, \phi_j] = \langle phi_i phi_j\rangle - \langle phi_i\rangle\langle phi_j\rangle\f$
 * which gives a component of the "phase metric".
 *
 * NOTE: for passing unit noise-weights, set MultiNoiseFloor->length=0 (but multiNoiseFloor==NULL is invalid)
 */
static double
XLALCovariance_Phi_ij ( const MultiLALDetector *multiIFO,		//!< [in] detectors to use
                        const MultiNoiseFloor *multiNoiseFloor,	//!< [in] corresponding noise floors for weights, NULL means unit-weights
                        const LALSegList *segList,			//!< [in] segment list
                        const intparams_t *params,			//!< [in] integration parameters
                        double* relerr_max				//!< [in] maximal error for integration
                      )
{
  XLAL_CHECK_REAL8 ( multiIFO != NULL, XLAL_EINVAL );
  UINT4 numDet = multiIFO->length;
  XLAL_CHECK_REAL8 ( numDet > 0, XLAL_EINVAL );

  // either no noise-weights given (multiNoiseFloor->length=0) or same number of detectors
  XLAL_CHECK_REAL8 ( multiNoiseFloor != NULL, XLAL_EINVAL );
  BOOLEAN haveNoiseWeights = (multiNoiseFloor->length > 0);
  XLAL_CHECK_REAL8 ( !haveNoiseWeights || (multiNoiseFloor->length == numDet), XLAL_EINVAL );

  XLAL_CHECK_REAL8 ( segList != NULL, XLAL_EINVAL );
  UINT4 Nseg = segList->length;

  /* sanity-check: don't allow any AM-coeffs being turned on here! */
  if ( params->amcomp1 != AMCOMP_NONE || params->amcomp2 != AMCOMP_NONE ) {
    XLALPrintError ( "%s: Illegal input, amcomp[12] must be set to AMCOMP_NONE!\n", __func__ );
    XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  /* store detector weights and accumulate total weight */
  REAL8 total_weight = 0.0, weights[numDet];
  for (UINT4 X = 0; X < numDet; X ++) {
    weights[X] = haveNoiseWeights ? multiNoiseFloor->sqrtSn[X] : 1.0;
    total_weight += weights[X];
  }
  XLAL_CHECK_REAL8 (total_weight > 0, XLAL_EDOM, "Detectors noise-floors given but all zero!" );

  /* ---------- set up GSL integration ---------- */

  /* constants passed to GSL integration function */
  const double epsrel = params->epsrel;
  const double epsabs = params->epsabs;
  const size_t limit = 512;

  /* array of structs storing input and output info for GSL integration */
  UINT4 intN[Nseg][numDet];
  typedef struct {
    intparams_t par;			/* integration parameters */
    double ti;				/* integration start time */
    double tf;				/* integration end time */
    gsl_integration_workspace *wksp ;	/* GSL integration workspace */
    double av_ij;			/* value of integral <phi_i phi_j> */
    double av_i;			/* value of integral <phi_i> */
    double av_j;			/* value of integral <phi_j> */
    double av_ij_err_sq;		/* squared error in value of integral <phi_i phi_j> */
    double av_i_err_sq;			/* squared error in value of integral <phi_i> */
    double av_j_err_sq;			/* squared error in value of integral <phi_j> */
  } InputOutputInfo;
  InputOutputInfo *intInOut[Nseg][numDet];

  // loop over segments and detectors
  for ( size_t k = 0; k < Nseg; ++k ) {
    for ( size_t X = 0; X < numDet; ++X ) {

      intparams_t par = (*params);   /* struct-copy, as the 'coord' field has to be changeable */

      // set start time and time span for this segment
      const LIGOTimeGPS *startTime = &(segList->segs[k].start);
      const LIGOTimeGPS *endTime   = &(segList->segs[k].end);
      const REAL8 Tspan = XLALGPSDiff( endTime, startTime );
      par.startTime = XLALGPSGetREAL8 ( startTime );
      par.Tspan     = Tspan;

      // set detector for phase integrals
      par.site = &multiIFO->sites[X];

      // calculate how many 'units' integral is split into
      intN[k][X] = (UINT4) ceil ( par.Tspan / par.intT );

      // allocate memory of input/output info
      intInOut[k][X] = XLALCalloc( intN[k][X], sizeof(*intInOut[k][X]) );
      XLAL_CHECK_REAL8( intInOut[k][X] != NULL, XLAL_ENOMEM );

      // initialise input info
      for ( size_t n = 0; n < intN[k][X]; ++n ) {

        InputOutputInfo *io = &intInOut[k][X][n];

        io->par = par;

        const double dT = 1.0 / intN[k][X];
        io->ti = 1.0 * n * dT;
        io->tf = MYMIN( (n+1.0) * dT, 1.0 );

        XLAL_CHECK_REAL8( ( io->wksp = gsl_integration_workspace_alloc(limit) ) != NULL, XLAL_ENOMEM );

      } /* for n < intN */

    } // for X < numDet
  } // for k < Nseg

  /* ---------- perform GSL integration ---------- */

  // turn off GSL error handling
  gsl_error_handler_t *saveGSLErrorHandler;
  saveGSLErrorHandler = gsl_set_error_handler_off();

  // loop over segments, detectors, and integration 'units'
  // together using a single index, for parallelisation
  {
    size_t index_max = 0;
    for ( size_t k = 0; k < Nseg; ++k ) {
      for ( size_t X = 0; X < numDet; ++X ) {
        index_max += intN[k][X];
      }
    }
#pragma omp parallel for schedule(static)
    for ( size_t indx = 0; indx < index_max; ++indx )
      {
        // break single index into 'k', 'X', and 'n'
        size_t k = 0, X = 0, n = 0, index_break = indx + 1;
        for ( k = 0; k < Nseg; ++k ) {
          for ( X = 0; X < numDet; ++X ) {
            if ( index_break > intN[k][X] ) {
              index_break -= intN[k][X];
            } else {
              n = index_break - 1;
              index_break = 0;
              break;
            }
          }
          if ( index_break == 0 ) {
            break;
          }
        }

        InputOutputInfo *io = &intInOut[k][X][n];

        gsl_function integrand;
        integrand.params = (void*)&( io->par );

        int stat;
        double res, abserr;

        const double scale1 = GET_COORD_SCALE(io->par.coordSys, io->par.coord1);
        const double scale2 = GET_COORD_SCALE(io->par.coordSys, io->par.coord2);
        const double scale12 = scale1 * scale2;

        /* compute <phi_i phi_j> */
        integrand.function = &CW_am1_am2_Phi_i_Phi_j;
        stat = gsl_integration_qag (&integrand, io->ti, io->tf, epsabs, epsrel, limit, GSL_INTEG_GAUSS61, io->wksp, &res, &abserr);
        if ( stat != 0 ) {
          XLALPrintWarning ( "\n%s: GSL-integration 'gsl_integration_qag()' of <Phi_i Phi_j> did not reach requested precision!\n", __func__ );
          XLALPrintWarning ( "xlalErrno=%i, seg=%zu, av_ij_n=%g, abserr=%g\n", io->par.errnum, n, res, abserr );
          io->av_ij = GSL_NAN;
        } else {
          io->av_ij = scale12 * res;
          io->av_ij_err_sq = SQUARE( scale12 * abserr);
        }

        /* compute <phi_i> */
        integrand.function = &CW_Phi_i;
        io->par.coord = io->par.coord1;
        stat = gsl_integration_qag (&integrand, io->ti, io->tf, epsabs, epsrel, limit, GSL_INTEG_GAUSS61, io->wksp, &res, &abserr);
        if ( stat != 0 ) {
          XLALPrintWarning ( "\n%s: GSL-integration 'gsl_integration_qag()' of <Phi_i> did not reach requested precision!\n", __func__ );
          XLALPrintWarning ( "xlalErrno=%i, seg=%zu, av_i_n=%g, abserr=%g\n", io->par.errnum, n, res, abserr );
          io->av_i = GSL_NAN;
        } else {
          io->av_i = scale1 * res;
          io->av_i_err_sq = SQUARE( scale1 * abserr);
        }

        /* compute <phi_j> */
        integrand.function = &CW_Phi_i;
        io->par.coord = io->par.coord2;
        stat = gsl_integration_qag (&integrand, io->ti, io->tf, epsabs, epsrel, limit, GSL_INTEG_GAUSS61, io->wksp, &res, &abserr);
        if ( stat != 0 ) {
          XLALPrintWarning ( "\n%s: GSL-integration 'gsl_integration_qag()' of <Phi_j> did not reach requested precision!\n", __func__ );
          XLALPrintWarning ( "xlalErrno=%i, seg=%zu, av_j_n=%g, abserr=%g\n", io->par.errnum, n, res, abserr );
          io->av_j = GSL_NAN;
        } else {
          io->av_j = scale2 * res;
          io->av_j_err_sq = SQUARE( scale2 * abserr);
        }

      } /* for indx < index_max */
  }

  // restore GSL error handling
  gsl_set_error_handler( saveGSLErrorHandler );

  /* ---------- compute final result ---------- */

  double ret = 0, relerr_max_sq = 0;

  // loop over segments
  for (UINT4 k = 0; k < Nseg; k ++) {

    double av_ij = 0, av_ij_err_sq = 0;
    double av_i = 0, av_i_err_sq = 0;
    double av_j = 0, av_j_err_sq = 0;

    // loop over detectors and integration 'units'
    for (UINT4 X = 0; X < numDet; X ++) {
      for ( size_t n = 0; n < intN[k][X]; ++n ) {

        const InputOutputInfo *io = &intInOut[k][X][n];

        /* accumulate <phi_i phi_j> */
        av_ij += weights[X] * io->av_ij;
        av_ij_err_sq += weights[X] * io->av_ij_err_sq;

        /* accumulate <phi_i> */
        av_i += weights[X] * io->av_i;
        av_i_err_sq += weights[X] * io->av_i_err_sq;

        /* accumulate <phi_j> */
        av_j += weights[X] * io->av_j;
        av_j_err_sq += weights[X] * io->av_j_err_sq;

      } /* for n < intN */

    } // for X < numDet

    // normalise by total weight
    av_ij /= total_weight;
    av_i /= total_weight;
    av_j /= total_weight;
    av_ij_err_sq /= total_weight;
    av_i_err_sq /= total_weight;
    av_j_err_sq /= total_weight;

    // work out maximum relative error for this segment
    const double av_ij_relerr = RELERR( sqrt(av_ij_err_sq), fabs(av_ij) );
    const double av_i_relerr = RELERR( sqrt(av_i_err_sq), fabs (av_i) );
    const double av_j_relerr = RELERR( sqrt(av_j_err_sq), fabs (av_j) );
    const double relerr_max_k = MYMAX ( av_ij_relerr, MYMAX ( av_i_relerr, av_j_relerr ) );

    ret += av_ij - av_i * av_j;
    relerr_max_sq += SQUARE( relerr_max_k );

  } // for k < Nseg

  ret /= Nseg;
  relerr_max_sq /= Nseg;

  if ( relerr_max )
    (*relerr_max) = sqrt( relerr_max_sq );

  /* ----- cleanup ----- */

  for ( size_t k = 0; k < Nseg; ++k ) {
    for ( size_t X = 0; X < numDet; ++X ) {
      for ( size_t n = 0; n < intN[k][X]; ++n ) {
        InputOutputInfo *io = &intInOut[k][X][n];
        gsl_integration_workspace_free( io->wksp );
      }
      XLALFree( intInOut[k][X] );
    }
  }

  return ret;

} /* XLALCovariance_Phi_ij() */


/**
 * Calculate an approximate "phase-metric" with the specified parameters.
 *
 * Note: if this function is called with multiple detectors, the phase components
 * are averaged over detectors as well as time. This is a somewhat ad-hoc approach;
 * if you want a more rigorous multi-detector metric you need to use the full
 * Fstat-metric, as computed by XLALComputeDopplerFstatMetric().
 *
 * Return NULL on error.
 */
DopplerPhaseMetric *
XLALComputeDopplerPhaseMetric ( const DopplerMetricParams *metricParams,  	/**< input parameters determining the metric calculation */
                                const EphemerisData *edat			/**< ephemeris data */
                              )
{
  intparams_t XLAL_INIT_DECL(intparams);

  /* ---------- sanity/consistency checks ---------- */
  XLAL_CHECK_NULL ( metricParams != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( edat != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( XLALSegListIsInitialized ( &(metricParams->segmentList) ), XLAL_EINVAL, "Passed un-initialzied segment list 'metricParams->segmentList'\n");
  UINT4 Nseg = metricParams->segmentList.length;

  UINT4 dim = metricParams->coordSys.dim;
  const DopplerCoordinateSystem *coordSys = &(metricParams->coordSys);
  // ----- check that {n2x_equ, n2y_equ} are not used at the equator (delta=0), as metric is undefined there
  BOOLEAN have_n2xy = 0;
  for ( UINT4 i = 0; i < dim; i ++ ) {
    if ( (coordSys->coordIDs[i] == DOPPLERCOORD_N2X_EQU) || ( coordSys->coordIDs[i] == DOPPLERCOORD_N2Y_EQU) ) {
      have_n2xy = 1;
    }
  }
  BOOLEAN at_equator = (metricParams->signalParams.Doppler.Delta == 0);
  XLAL_CHECK_NULL ( !(at_equator && have_n2xy), XLAL_EINVAL, "Can't use 'n2x_equ','n2y_equ' at equator (Delta=0): metric is singular there");

  // ----- allocate output DopplerPhaseMetric() struct
  DopplerPhaseMetric *metric = XLALCalloc( 1, sizeof(*metric) );
  XLAL_CHECK_NULL( metric != NULL, XLAL_ENOMEM );

  // ----- compute segment-list midTime
  LIGOTimeGPS midTime;
  {
    const LIGOTimeGPS *startTime = &(metricParams->segmentList.segs[0].start);
    const LIGOTimeGPS *endTime   = &(metricParams->segmentList.segs[Nseg-1].end);
    const REAL8 Tspan = XLALGPSDiff( endTime, startTime );
    midTime = *startTime;
    XLALGPSAdd( &midTime, 0.5 * Tspan );
  }

  /* ---------- prepare output metric ---------- */
  if ( (metric->g_ij = gsl_matrix_calloc ( dim, dim )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc(%d, %d) failed.\n\n", __func__, dim, dim );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  /* ---------- set up integration parameters ---------- */
  intparams.coordSys = coordSys;
  intparams.edat = edat;
  intparams.refTime   = XLALGPSGetREAL8 ( &midTime );   /* always compute metric at midTime, transform to refTime later */
  intparams.dopplerPoint = &(metricParams->signalParams.Doppler);
  intparams.detMotionType = metricParams->detMotionType;
  intparams.approxPhase = metricParams->approxPhase;
  /* deactivate antenna-patterns for phase-metric */
  intparams.amcomp1 = AMCOMP_NONE;
  intparams.amcomp2 = AMCOMP_NONE;

  /* if using 'global correlation' frequency variables, determine the highest spindown order: */
  UINT4 maxorder = findHighestGCSpinOrder ( coordSys );
  if ( maxorder > 0 ) {

    /* compute rOrb(t) derivatives at reference time */
    if ( (intparams.rOrb_n = XLALComputeOrbitalDerivatives ( maxorder, &intparams.dopplerPoint->refTime, edat )) == NULL ) {
      XLALPrintError ("%s: XLALComputeOrbitalDerivatives() failed.\n", __func__);
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }
    if (lalDebugLevel & LALINFOBIT) {
      /* diagnostic / debug output */
      fprintf( stderr, "%s: rOrb_n(%d) = [ ", __func__, intparams.dopplerPoint->refTime.gpsSeconds );
      for ( UINT4 n = 0; n < intparams.rOrb_n->length; n++ ) {
        fprintf( stderr, "[%g, %g, %g]%s", intparams.rOrb_n->data[n][0], intparams.rOrb_n->data[n][1], intparams.rOrb_n->data[n][2],
                 (n < intparams.rOrb_n->length -1 ) ? ", " : " ]\n" );
      }
    }

  }

  metric->maxrelerr = 0;
  double err = 0;

  /* ========== use numerically-robust method to compute metric ========== */

  /* allocate memory for coordinate transform */
  gsl_matrix *transform = gsl_matrix_alloc(dim, dim);
  XLAL_CHECK_NULL( transform != NULL, XLAL_ENOMEM );
  gsl_matrix_set_identity( transform );
  intparams.coordTransf = transform;

  for ( size_t n = 1; n <= dim; ++n ) {

    /* NOTE: this level of accuracy is only achievable *without* AM-coefficients involved
     * which are computed in REAL4 precision. For the current function this is OK, as this
     * function is only supposed to compute *pure* phase-derivate covariances.
     */
    intparams.epsrel = 1e-6;
    /* we need an abs-cutoff as well, as epsrel can be too restrictive for small integrals */
    intparams.epsabs = 1e-3;
    /* NOTE: this numerical integration still runs into problems when integrating over
     * long durations (~O(23d)), as the integrands are oscillatory functions on order of ~1d
     * and convergence degrades.
     * As a solution, we split the integral into 'intN' units of ~1 day duration, and compute
     * the final integral as a sum over partial integrals.
     * It is VERY important to ensure that 'intT' is not exactly 1 day, since then an integrand
     * with period ~1 day may integrate to zero, which is both slower and MUCH more difficult
     * for the numerical integration functions (since the integrand them becomes small with
     * respect to any integration errors).
     */
    intparams.intT = 0.9 * LAL_DAYSID_SI;

    /* allocate memory for Cholesky decomposition */
    gsl_matrix *cholesky = gsl_matrix_alloc(n, n);
    XLAL_CHECK_NULL( cholesky != NULL, XLAL_ENOMEM );

    /* create views of n-by-n submatrices of metric and coordinate transform */
    gsl_matrix_view g_ij_n = gsl_matrix_submatrix( metric->g_ij, 0, 0, n, n );
    gsl_matrix_view transform_n = gsl_matrix_submatrix( transform, 0, 0, n, n );

    /* try this loop a certain number of times */
    const size_t max_tries = 64;
    size_t tries = 0;
    while ( ++tries <= max_tries ) {

      /* ----- compute last row/column of n-by-n submatrix of metric ----- */
      for ( size_t i = 0; i < n; ++i ) {
        const size_t j = n - 1;

        /* g_ij = [Phi_i, Phi_j] */
        intparams.coord1 = i;
        intparams.coord2 = j;
        REAL8 gg = XLALCovariance_Phi_ij ( &metricParams->multiIFO, &metricParams->multiNoiseFloor, &metricParams->segmentList,
                                           &intparams, &err );
        XLAL_CHECK_NULL( !gsl_isnan(gg), XLAL_EFUNC, "%s: integration of phase metric g_{i=%zu,j=%zu} failed (n=%zu, tries=%zu)", __func__, i, j, n, tries );
        gsl_matrix_set (&g_ij_n.matrix, i, j, gg);
        gsl_matrix_set (&g_ij_n.matrix, j, i, gg);
        metric->maxrelerr = MYMAX ( metric->maxrelerr, err );

      } /* for i < n */

      /* ----- compute L D L^T Cholesky decomposition of metric ----- */
      XLAL_CHECK_NULL( XLALCholeskyLDLTDecompMetric( &cholesky, &g_ij_n.matrix ) == XLAL_SUCCESS, XLAL_EFUNC );
      gsl_vector_view D = gsl_matrix_diagonal( cholesky );
      if ( ( tries > 1 ) && ( lalDebugLevel & LALINFOBIT ) ) {
        /* diagnostic / debug output */
        fprintf(stdout, "%s: n=%zu, try=%zu, Cholesky diagonal =", __func__, n, tries);
        XLALfprintfGSLvector(stdout, "%0.2e", &D.vector);
      }

      /* ----- check that all diagonal elements D are positive after at least 1 try; if so, exit try loop */
      if ( ( tries > 1 ) && ( gsl_vector_min( &D.vector ) > 0.0 ) ) {
        break;
      }

      /* zero out all but last row of L, since we do not want to coordinates before 'n' */
      if ( n > 1 ) {
        gsl_matrix_view cholesky_nm1 = gsl_matrix_submatrix( cholesky, 0, 0, n - 1, n );
        gsl_matrix_set_identity( &cholesky_nm1.matrix );
      }

      /* multiply transform by inverse of L (with unit diagonal), to transform 'n'th coordinates so that metric is diagonal */
      gsl_blas_dtrsm( CblasLeft, CblasLower, CblasNoTrans, CblasUnit, 1.0, cholesky, &transform_n.matrix );

      /* decrease relative error tolerances; don't do this too quickly, since
         too-stringent tolerances may make GSL integration fail to converge */
      intparams.epsrel = intparams.epsrel * 0.9;
      /* decrease absolute error tolerances; don't do this too quickly, since
         too-stringent tolerances may make GSL integration fail to converge,
         and only after a certain number of tries, otherwise small integrals
         fail to converge */
      if ( tries >= 8  ) {
        intparams.epsabs = intparams.epsabs * 0.9;
      }
      /* reduce the length of integration time units, but stop at 900s,
         and ensure that 'intT' does NOT become divisible by 1/day, for
         the same reason given at the initialisation of 'intT' above */
      intparams.intT = MYMAX(900, intparams.intT * 0.9);

    }
    XLAL_CHECK_NULL( tries <= max_tries, XLAL_EMAXITER, "%s: convergence of phase metric failed (n=%zu)", __func__, n );

    /* free memory which is then re-allocated in next loop */
    gsl_matrix_free( cholesky );

  }

  /* apply inverse of transform^T to metric, to get back original coordinates */
  gsl_matrix_transpose( transform );
  XLAL_CHECK_NULL( XLALInverseTransformMetric( &metric->g_ij, transform, metric->g_ij ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* transform phase metric reference time from midTime to refTime */
  const REAL8 Dtau = XLALGPSDiff( &(metricParams->signalParams.Doppler.refTime), &midTime );
  XLAL_CHECK_NULL( XLALChangeMetricReferenceTime( &metric->g_ij, NULL, metric->g_ij, &(metricParams->coordSys), Dtau ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* ----- if requested, project g_ij onto coordinate 'projectCoord' */
  if ( metricParams->projectCoord >= 0 )
    {
      UINT4 projCoord = (UINT4)metricParams->projectCoord;
      if ( XLALProjectMetric ( &metric->g_ij, metric->g_ij, projCoord ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: failed to project g_ij onto coordinate '%d'. errno=%d\n", __func__, projCoord, xlalErrno );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }
    } /* if projectCoordinate >= 0 */

  /* ---- check that metric is positive definite, by checking determinants of submatrices */
  for ( size_t n = 1; n <= dim; ++n ) {
    gsl_matrix_view g_ij_n = gsl_matrix_submatrix( metric->g_ij, 0, 0, n, n );
    const double det_n = XLALMetricDeterminant( &g_ij_n.matrix );
    XLAL_CHECK_NULL( det_n > 0, XLAL_EFAILED, "%s: could not compute a positive-definite phase metric (n=%zu, det_n=%0.3e)", __func__, n, det_n );
  }
  if (lalDebugLevel & LALINFOBIT) {
    /* diagnostic / debug output */
    fprintf(stdout, "%s: phase metric:", __func__);
    XLALfprintfGSLmatrix(stdout, "%0.15e", metric->g_ij);
  }

  /*  attach the metricParams struct as 'meta-info' to the output */
  metric->meta = (*metricParams);

  /* free memory */
  XLALDestroyVect3Dlist ( intparams.rOrb_n );
  gsl_matrix_free(transform);

  return metric;

} /* XLALComputeDopplerPhaseMetric() */


/** Free a DopplerPhaseMetric structure */
void
XLALDestroyDopplerPhaseMetric ( DopplerPhaseMetric *metric )
{
  if ( !metric )
    return;

  gsl_matrix_free ( metric->g_ij );

  XLALFree ( metric );

  return;

} /* XLALDestroyDopplerPhaseMetric() */


/**
 * Calculate the general (single-segment coherent, or multi-segment semi-coherent)
 * *full* (multi-IFO) Fstat-metrix and the Fisher-matrix derived in \cite Prix07 .
 *
 * The semi-coherent metrics \f$g_{ij}\f$ over \f$N\f$ segments are computed according to
 *
 * \f[ \overline{g}_{ij} \equiv \frac{1}{N} \sum_{k=1}^{N} g_{ij,k} \f]
 *
 * where \f$g_{ij,k}\f$ is the coherent single-segment metric of segment k
 *
 * Note: The returned DopplerFstatMetric struct contains the matrices
 * g_ij (the phase metric), gF_ij (the F-metric), gFav_ij (the average F-metric),
 * m1_ij, m2_ij, m3_ij (auxiliary matrices)
 * and Fisher_ab (the full 4+n dimensional Fisher matrix).
 *
 * The returned metric struct also carries the meta-info about
 * the metrics in the field 'DopplerMetricParams meta'.
 *
 * Note2: for backwards-compatibility, we treat params->Nseg==0 equivalent to Nseg==1, ie
 * compute the coherent, single-segment metrics
 *
 * Return NULL on error.
 */
DopplerFstatMetric *
XLALComputeDopplerFstatMetric ( const DopplerMetricParams *metricParams,  	/**< input parameters determining the metric calculation */
                                const EphemerisData *edat			/**< ephemeris data */
                              )
{
  XLAL_CHECK_NULL ( metricParams, XLAL_EINVAL, "Invalid NULL input 'metricParams'\n" );
  XLAL_CHECK_NULL ( XLALSegListIsInitialized ( &(metricParams->segmentList) ), XLAL_EINVAL, "Passed un-initialzied segment list 'metricParams->segmentList'\n");

  UINT4 dim = metricParams->coordSys.dim;
  const DopplerCoordinateSystem *coordSys = &(metricParams->coordSys);
  // ----- check that {n2x_equ, n2y_equ} are not used at the equator (delta=0), as metric is undefined there
  BOOLEAN have_n2xy = 0;
  for ( UINT4 i = 0; i < dim; i ++ ) {
    if ( (coordSys->coordIDs[i] == DOPPLERCOORD_N2X_EQU) || ( coordSys->coordIDs[i] == DOPPLERCOORD_N2Y_EQU) ) {
      have_n2xy = 1;
    }
  }
  BOOLEAN at_equator = (metricParams->signalParams.Doppler.Delta == 0);
  XLAL_CHECK_NULL ( !(at_equator && have_n2xy), XLAL_EINVAL, "Can't use 'n2x_equ','n2y_equ' at equator (Delta=0): metric is singular there");

  UINT4 Nseg = metricParams->segmentList.length;	// number of semi-coherent segments to average metrics over
  XLAL_CHECK_NULL ( Nseg >= 1, XLAL_EINVAL, "Got empty segment list metricParams->segmentList, needs to contain at least 1 segments\n");

  DopplerFstatMetric *metric = NULL;

  LALSegList segList_k;
  LALSeg segment_k;
  XLALSegListInit( &segList_k );	// prepare single-segment list containing segment k
  segList_k.arraySize = 1;
  segList_k.length = 1;
  segList_k.segs = &segment_k;

  DopplerMetricParams metricParams_k = (*metricParams);	// *copy* of input parameters to be used for single-segment coherent metrics
  metricParams_k.segmentList = segList_k;

  for ( UINT4 k = 0; k < Nseg; k ++ )
    {
      DopplerFstatMetric *metric_k;	// per-segment coherent metric
      FmetricAtoms_t *atoms = NULL;

      // setup 1-segment segment-list pointing k-th segment
      metricParams_k.segmentList.segs[0] = metricParams->segmentList.segs[k];

      /* ---------- compute Fmetric 'atoms', ie the averaged <a^2>, <a b Phi_i>, <a^2 Phi_i Phi_j>, etc ---------- */
      if ( (atoms = XLALComputeAtomsForFmetric ( &metricParams_k, edat )) == NULL ) {
        XLALPrintError ("%s: XLALComputeAtomsForFmetric() failed. errno = %d\n\n", __func__, xlalErrno );
        XLAL_ERROR_NULL( XLAL_EFUNC );
      }

      /* ----- compute the F-metric gF_ij and related matrices ---------- */
      REAL8 cosi = metricParams->signalParams.Amp.cosi;
      REAL8 psi  = metricParams->signalParams.Amp.psi;

      if ( (metric_k = XLALComputeFmetricFromAtoms ( atoms, cosi, psi)) == NULL ) {
        XLALPrintError ("%s: XLALComputeFmetricFromAtoms() failed, errno = %d\n\n", __func__, xlalErrno );
        XLALDestroyFmetricAtoms ( atoms );
        XLAL_ERROR_NULL( XLAL_EFUNC );
      }

      /* ----- compute the full 4+n dimensional Fisher matrix ---------- */
      if ( (metric_k->Fisher_ab = XLALComputeFisherFromAtoms ( atoms, metricParams->signalParams.Amp )) == NULL ) {
        XLALPrintError ("%s: XLALComputeFisherFromAtoms() failed. errno = %d\n\n", __func__, xlalErrno );
        XLALDestroyFmetricAtoms ( atoms );
        XLALDestroyDopplerFstatMetric ( metric );
        XLAL_ERROR_NULL( XLAL_EFUNC );
      }

      XLALDestroyFmetricAtoms ( atoms );

      /* ---- add DopplerFstatMetric results over segments ----- */
      if ( metric == NULL ) {
        metric = metric_k;   // just use DopplerFstatMetric from 1st segment
      } else {

        // add matrices
        gsl_matrix_add ( metric->gF_ij,     metric_k->gF_ij );
        gsl_matrix_add ( metric->gFav_ij,   metric_k->gFav_ij );
        gsl_matrix_add ( metric->m1_ij,     metric_k->m1_ij );
        gsl_matrix_add ( metric->m2_ij,     metric_k->m2_ij );
        gsl_matrix_add ( metric->m3_ij,     metric_k->m3_ij );
        gsl_matrix_add ( metric->Fisher_ab, metric_k->Fisher_ab );

        // add errors in quadrature
        metric->maxrelerr = sqrt ( SQUARE(metric->maxrelerr) + SQUARE(metric_k->maxrelerr) );

        // add SNR^2
        metric->rho2 += metric_k->rho2;

        XLALDestroyDopplerFstatMetric ( metric_k );

      }

    } // for k < Nseg

  // semi-coherent metric: <g_ij> = (1/Nseg) sum_{k=1}^Nseg g_{k,ij}
  {
    const double scale = 1.0/Nseg;

    // scale matrices
    gsl_matrix_scale ( metric->gF_ij,     scale );
    gsl_matrix_scale ( metric->gFav_ij,   scale );
    gsl_matrix_scale ( metric->m1_ij,     scale );
    gsl_matrix_scale ( metric->m2_ij,     scale );
    gsl_matrix_scale ( metric->m3_ij,     scale );
    gsl_matrix_scale ( metric->Fisher_ab, scale );

    // scale errors
    metric->maxrelerr *= scale;

    // scale SNR^2
    metric->rho2 *= scale;
  }

  /* ----- if requested, project gF_ij and gFav_ij onto coordinate 'projectCoord' */
  if ( metricParams->projectCoord >= 0 )
    {
      UINT4 projCoord = (UINT4)metricParams->projectCoord;
      if ( XLALProjectMetric ( &metric->gF_ij, metric->gF_ij, projCoord ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: failed to project gF_ij onto coordinate '%d'. errno=%d\n", __func__, projCoord, xlalErrno );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }
      if ( XLALProjectMetric ( &metric->gFav_ij, metric->gFav_ij, projCoord ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: failed to project gFav_ij onto coordinate '%d'. errno=%d\n", __func__, projCoord, xlalErrno );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }
    } /* if projectCoordinate >= 0 */

  /*  attach the metricParams struct as 'meta-info' to the output */
  metric->meta = (*metricParams);

  return metric;

} /* XLALComputeDopplerFstatMetric() */


/** Free a DopplerFstatMetric structure */
void
XLALDestroyDopplerFstatMetric ( DopplerFstatMetric *metric )
{
  if ( !metric )
    return;

  gsl_matrix_free ( metric->gF_ij );
  gsl_matrix_free ( metric->gFav_ij );
  gsl_matrix_free ( metric->m1_ij );
  gsl_matrix_free ( metric->m2_ij );
  gsl_matrix_free ( metric->m3_ij );
  gsl_matrix_free ( metric->Fisher_ab );

  XLALFree ( metric );

  return;

} /* XLALDestroyDopplerFstatMetric() */


/**
 * Function to the compute the FmetricAtoms_t, from which the F-metric and Fisher-matrix can be computed.
 *
 * NOTE: if MultiNoiseFloor.length=0, unit noise weights are assumed.
 */
FmetricAtoms_t*
XLALComputeAtomsForFmetric ( const DopplerMetricParams *metricParams,  	/**< input parameters determining the metric calculation */
			     const EphemerisData *edat			/**< ephemeris data */
			     )
{
  FmetricAtoms_t *ret;		/* return struct */
  intparams_t XLAL_INIT_DECL(intparams);

  UINT4 dim, i=-1, j=-1, X;		/* index counters */
  REAL8 A, B, C;			/* antenna-pattern coefficients (gsl-integrated) */

  const DopplerCoordinateSystem *coordSys;

  REAL8 max_relerr = 0;
  REAL8 relerr_thresh = 1e-2;	/* relatively tolerant integration-threshold */

  /* ---------- sanity/consistency checks ---------- */
  if ( !metricParams || !edat ) {
    XLALPrintError ("\n%s: Illegal NULL pointer passed!\n\n", __func__);
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }
  UINT4 numDet = metricParams->multiIFO.length;
  XLAL_CHECK_NULL ( numDet >= 1, XLAL_EINVAL );
  XLAL_CHECK_NULL ( XLALSegListIsInitialized ( &(metricParams->segmentList) ), XLAL_EINVAL, "Passed un-initialzied segment list 'metricParams->segmentList'\n");
  UINT4 Nseg = metricParams->segmentList.length;
  XLAL_CHECK_NULL ( Nseg == 1, XLAL_EINVAL, "Segment list must only contain Nseg=1 segments, got Nseg=%d", Nseg );

  BOOLEAN haveNoiseWeights = (metricParams->multiNoiseFloor.length > 0);
  XLAL_CHECK_NULL ( !haveNoiseWeights || metricParams->multiNoiseFloor.length == numDet, XLAL_EINVAL );

  LIGOTimeGPS *startTime = &(metricParams->segmentList.segs[0].start);
  LIGOTimeGPS *endTime   = &(metricParams->segmentList.segs[0].end);
  REAL8 Tspan = XLALGPSDiff( endTime, startTime );

  const LIGOTimeGPS *refTime   = &(metricParams->signalParams.Doppler.refTime);

  dim = metricParams->coordSys.dim;	/* shorthand: number of Doppler dimensions */
  coordSys = &(metricParams->coordSys);

  /* ----- create output structure ---------- */
  XLAL_CHECK_NULL ( ( ret = XLALCalloc(1,sizeof(*ret))) != NULL, XLAL_ENOMEM, "%s: XLALCalloc(1,%zu) failed.\n", __func__, sizeof(*ret));
  XLAL_CHECK_NULL ( ( ret->a_a_i = gsl_vector_calloc (dim)) != NULL, XLAL_ENOMEM, "%s: a_a_i = gsl_vector_calloc (%d) failed.\n", __func__, dim );
  XLAL_CHECK_NULL ( ( ret->a_b_i = gsl_vector_calloc (dim)) != NULL, XLAL_ENOMEM, "%s: a_b_i = gsl_vector_calloc (%d) failed.\n", __func__, dim );
  XLAL_CHECK_NULL ( ( ret->b_b_i = gsl_vector_calloc (dim)) != NULL, XLAL_ENOMEM, "%s: b_b_i = gsl_vector_calloc (%d) failed.\n", __func__, dim );
  XLAL_CHECK_NULL ( ( ret->a_a_i_j = gsl_matrix_calloc (dim, dim)) != NULL, XLAL_ENOMEM, "%s: a_a_i_j = gsl_matrix_calloc (%d,%d) failed.\n", __func__, dim, dim );
  XLAL_CHECK_NULL ( ( ret->a_b_i_j = gsl_matrix_calloc (dim, dim)) != NULL, XLAL_ENOMEM, "%s: a_b_i_j = gsl_matrix_calloc (%d,%d) failed.\n", __func__, dim, dim );
  XLAL_CHECK_NULL ( ( ret->b_b_i_j = gsl_matrix_calloc (dim, dim)) != NULL, XLAL_ENOMEM, "%s: b_b_i_j = gsl_matrix_calloc (%d,%d) failed.\n", __func__, dim, dim );

  /* ---------- set up integration parameters ---------- */
  intparams.coordSys = coordSys;
  intparams.detMotionType = metricParams->detMotionType;
  intparams.dopplerPoint = &metricParams->signalParams.Doppler;
  intparams.startTime = XLALGPSGetREAL8 ( startTime );
  intparams.refTime   = XLALGPSGetREAL8 ( refTime );
  intparams.Tspan = Tspan;
  intparams.edat = edat;

  /* NOTE: this level of accuracy should be compatible with AM-coefficients involved
   * which are computed in REAL4 precision. We therefor cannot go lower than this it seems,
   * otherwise the gsl-integration fails to converge in some cases.
   */
  intparams.epsrel = 1e-4;
  intparams.epsabs = 0;
  /* NOTE: this numerical integration runs into problems when integrating over
   * several days (~O(5d)), as the integrands are oscillatory functions on order of ~1/4d
   * and convergence degrades.
   * As a solution, we split the integral into 'intN' units of 1/4 day duration, and compute
   * the final integral as a sum over partial integrals.
   */
  intparams.intT = 0.25 * LAL_DAYSID_SI;

  /* if using 'global correlation' frequency variables, determine the highest spindown order: */
  UINT4 maxorder = findHighestGCSpinOrder ( coordSys );
  if ( maxorder > 0 ) {

    /* compute rOrb(t) derivatives at reference time */
    if ( (intparams.rOrb_n = XLALComputeOrbitalDerivatives ( maxorder, &intparams.dopplerPoint->refTime, edat )) == NULL ) {
      XLALPrintError ("%s: XLALComputeOrbitalDerivatives() failed.\n", __func__);
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }

  }

  /* ----- integrate antenna-pattern coefficients A, B, C */
  REAL8 sum_weights = 0;
  A = B = C = 0;
  for ( X = 0; X < numDet; X ++ )
    {
      REAL8 weight = haveNoiseWeights ? metricParams->multiNoiseFloor.sqrtSn[X] : 1.0;
      sum_weights += weight;
      REAL8 av, relerr;
      intparams.site = &(metricParams->multiIFO.sites[X]);

      intparams.coord1 = -1;
      intparams.coord2 = -1;

      /* A = < a^2 > (67)*/
      intparams.amcomp1 = AMCOMP_A;
      intparams.amcomp2 = AMCOMP_A;
      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
      max_relerr = MYMAX ( max_relerr, relerr );
      if ( xlalErrno ) goto failed;
      A += weight * av;

      /* B = < b^2 > (67) */
      intparams.amcomp1 = AMCOMP_B;
      intparams.amcomp2 = AMCOMP_B;
      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
      max_relerr = MYMAX ( max_relerr, relerr );
      if ( xlalErrno ) goto failed;
      B += weight * av;

      /* C = < a b > (67) */
      intparams.amcomp1 = AMCOMP_A;
      intparams.amcomp2 = AMCOMP_B;
      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
      max_relerr = MYMAX ( max_relerr, relerr );
      if ( xlalErrno ) goto failed;
      C += weight * av;

    } /* for X < numDetectors */

  REAL8 norm_weight = 1.0 / sum_weights;

  A *= norm_weight;
  B *= norm_weight;
  C *= norm_weight;

  ret->a_a = A;
  ret->b_b = B;
  ret->a_b = C;

  /* ---------- compute components of the phase-metric ---------- */
  for ( i=0; i < dim; i ++ )
    {
      for ( j = 0; j <= i; j ++ )
	{
	  REAL8 a_a_i_j, b_b_i_j, a_b_i_j;
	  REAL8 a_a_i, b_b_i, a_b_i;
	  REAL8 a_a_j, b_b_j, a_b_j;

	  a_a_i_j = b_b_i_j = a_b_i_j = 0;
	  a_a_i = b_b_i = a_b_i = 0;
	  a_a_j = b_b_j = a_b_j = 0;

	  for ( X = 0; X < numDet; X ++ )
	    {
	      REAL8 weight = haveNoiseWeights ? metricParams->multiNoiseFloor.sqrtSn[X] : 1.0;
	      REAL8 av, relerr;
	      intparams.site = &(metricParams->multiIFO.sites[X]);

	      /* ------------------------------ */
	      intparams.coord1 = i;
	      intparams.coord2 = j;

	      /* <a^2 Phi_i Phi_j> */
	      intparams.amcomp1 = AMCOMP_A;
	      intparams.amcomp2 = AMCOMP_A;
	      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
              max_relerr = MYMAX ( max_relerr, relerr );
	      if ( xlalErrno ) goto failed;
	      a_a_i_j += weight * av;

	      /* <b^2 Phi_i Phi_j> */
	      intparams.amcomp1 = AMCOMP_B;
	      intparams.amcomp2 = AMCOMP_B;
	      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
              max_relerr = MYMAX ( max_relerr, relerr );
	      if ( xlalErrno ) goto failed;
	      b_b_i_j += weight * av;

	      /* <a b Phi_i Phi_j> */
	      intparams.amcomp1 = AMCOMP_A;
	      intparams.amcomp2 = AMCOMP_B;
	      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
              max_relerr = MYMAX ( max_relerr, relerr );
	      if ( xlalErrno ) goto failed;
	      a_b_i_j += weight * av;

	      /* ------------------------------ */
	      intparams.coord1 = i;
	      intparams.coord2 = -1;

	      /* <a^2 Phi_i> */
	      intparams.amcomp1 = AMCOMP_A;
	      intparams.amcomp2 = AMCOMP_A;
	      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
              max_relerr = MYMAX ( max_relerr, relerr );
	      if ( xlalErrno ) goto failed;
	      a_a_i += weight * av;

	      /* <b^2 Phi_i> */
	      intparams.amcomp1 = AMCOMP_B;
	      intparams.amcomp2 = AMCOMP_B;
	      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
              max_relerr = MYMAX ( max_relerr, relerr );
	      if ( xlalErrno ) goto failed;
	      b_b_i += weight * av;

	      /* <a b Phi_i> */
	      intparams.amcomp1 = AMCOMP_A;
	      intparams.amcomp2 = AMCOMP_B;
	      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
              max_relerr = MYMAX ( max_relerr, relerr );
	      if ( xlalErrno ) goto failed;
	      a_b_i += weight * av;

	      /* ------------------------------ */
	      intparams.coord1 = -1;
	      intparams.coord2 = j;

	      /* <a^2 Phi_j> */
	      intparams.amcomp1 = AMCOMP_A;
	      intparams.amcomp2 = AMCOMP_A;
	      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
              max_relerr = MYMAX ( max_relerr, relerr );
	      if ( xlalErrno ) goto failed;
	      a_a_j += weight * av;

	      /* <b^2 Phi_j> */
	      intparams.amcomp1 = AMCOMP_B;
	      intparams.amcomp2 = AMCOMP_B;
	      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
              max_relerr = MYMAX ( max_relerr, relerr );
	      if ( xlalErrno ) goto failed;
	      b_b_j += weight * av;

	      /* <a b Phi_j> */
	      intparams.amcomp1 = AMCOMP_A;
	      intparams.amcomp2 = AMCOMP_B;
	      av = XLALAverage_am1_am2_Phi_i_Phi_j ( &intparams, &relerr );
              max_relerr = MYMAX ( max_relerr, relerr );
	      if ( xlalErrno ) goto failed;
	      a_b_j += weight * av;

	    } /* for X < numDetectors */

	  gsl_vector_set (ret->a_a_i, i, a_a_i * norm_weight );
	  gsl_vector_set (ret->a_b_i, i, a_b_i * norm_weight );
	  gsl_vector_set (ret->b_b_i, i, b_b_i * norm_weight );


	  gsl_matrix_set (ret->a_a_i_j, i, j, a_a_i_j * norm_weight);
	  gsl_matrix_set (ret->a_a_i_j, j, i, a_a_i_j * norm_weight);

	  gsl_matrix_set (ret->a_b_i_j, i, j, a_b_i_j * norm_weight);
	  gsl_matrix_set (ret->a_b_i_j, j, i, a_b_i_j * norm_weight);

	  gsl_matrix_set (ret->b_b_i_j, i, j, b_b_i_j * norm_weight);
	  gsl_matrix_set (ret->b_b_i_j, j, i, b_b_i_j * norm_weight);

	} /* for j <= i */

    } /* for i < dim */

  /* return error-estimate */
  ret->maxrelerr = max_relerr;

  /* free memory */
  XLALDestroyVect3Dlist ( intparams.rOrb_n );

  /* FIXME: should probably not be hardcoded */
  if ( max_relerr > relerr_thresh )
    {
      XLALPrintError ("Maximal relative F-metric error too high: %.2e > %.2e\n", max_relerr, relerr_thresh );
      XLALDestroyFmetricAtoms ( ret );
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }
  else
    XLALPrintInfo ("\nMaximal relative error in F-metric: %.2e\n", max_relerr );

  return ret;

 failed:
  XLALDestroyFmetricAtoms ( ret );
  XLALPrintError ( "%s: XLALAverage_am1_am2_Phi_i_Phi_j() FAILED with errno = %d: am1 = %d, am2 = %d, coord1 = '%s', coord2 = '%s'\n",
		   __func__, xlalErrno, intparams.amcomp1, intparams.amcomp2,
                   GET_COORD_NAME(intparams.coordSys, intparams.coord1), GET_COORD_NAME(intparams.coordSys, intparams.coord2)
                 );
  XLAL_ERROR_NULL( XLAL_EFUNC );

} /* XLALComputeAtomsForFmetric() */


/**
 * Parse a detector-motion type string into the corresponding enum-number,
 */
int
XLALParseDetectorMotionString ( const CHAR *detMotionString )
{

  XLAL_CHECK( detMotionString, XLAL_EINVAL );

  for ( int i = 0; DetectorMotionNames[i].type > 0; ++i ) {
    if ( strcmp ( detMotionString, DetectorMotionNames[i].name ) != 0 )
      continue;
    return DetectorMotionNames[i].type;   /* found the right entry */
  }

  XLAL_ERROR ( XLAL_EINVAL, "Could not parse '%s' into a valid detector-motion type!", detMotionString );

} /* XLALParseDetectorMotionString() */


/**
 * Provide a pointer to a static string containing the DopplerCoordinate-name
 * cooresponding to the enum DopplerCoordinateID
 */
const CHAR *
XLALDetectorMotionName ( DetectorMotionType detMotionType )
{

  for ( int i = 0; DetectorMotionNames[i].type > 0; ++i ) {
    if ( DetectorMotionNames[i].type != detMotionType )
      continue;
    return DetectorMotionNames[i].name;   /* found the right entry */
  }

  XLAL_ERROR_NULL ( XLAL_EINVAL, "Could not parse '%d' into a valid detector-motion name!", detMotionType );

} /* XLALDetectorMotionName() */



/**
 * Parse a DopplerCoordinate-name into the corresponding DopplerCoordinateID
 */
int
XLALParseDopplerCoordinateString ( const CHAR *coordName )
{
  int i;

  if ( !coordName )
    XLAL_ERROR ( XLAL_EINVAL );

  for ( i=0; i < DOPPLERCOORD_LAST; i ++ )
    {
      if ( DopplerCoordinates[i].name && strcmp ( coordName, DopplerCoordinates[i].name ) )
	continue;
      return i;	/* found the right entry */
    }

  XLALPrintError ("\nCould not parse '%s' into a valid coordinate-ID!\n\n", coordName );
  XLAL_ERROR ( XLAL_EINVAL );

} /* XLALParseDopplerCoordinateString() */

/**
 * Given a LALStringVector of coordinate-names, parse them into a
 * 'DopplerCoordinateSystem', namely a list of coordinate-IDs
 */
int
XLALDopplerCoordinateNames2System ( DopplerCoordinateSystem *coordSys,	/**< [out] pointer to coord-system struct */
				    const LALStringVector *coordNames 	/**< [in] list of coordinate names */
				    )
{
  UINT4 i;

  if ( !coordSys || !coordNames )
    XLAL_ERROR ( XLAL_EINVAL );

  coordSys->dim = coordNames->length;
  for ( i=0; i < coordNames->length; i++ )
    {
      coordSys->coordIDs[i] = (DopplerCoordinateID)XLALParseDopplerCoordinateString ( coordNames->data[i] );
      if ( xlalErrno )
	XLAL_ERROR ( XLAL_EFUNC );
    }

  return XLAL_SUCCESS;

} /* XLALDopplerCoordinateNames2System() */



/**
 * Given a coordinate ID 'coordID', return its dimension within the given coordinate system 'coordSys',
 * or return -1 if 'coordID' is not found
 */
int XLALFindDopplerCoordinateInSystem ( const DopplerCoordinateSystem *coordSys, const DopplerCoordinateID coordID )
{
  for ( int i = 0; i < ((int)coordSys->dim); ++i )
    {
      if ( coordSys->coordIDs[i] == coordID )
        {
          return i;
        }
    }
  return -1;
}



/**
 * Provide a pointer to a static string containing the DopplerCoordinate-name
 * cooresponding to the enum DopplerCoordinateID
 */
const CHAR *
XLALDopplerCoordinateName ( DopplerCoordinateID coordID )
{
  if ( coordID >= DOPPLERCOORD_LAST ) {
    XLAL_ERROR_NULL ( XLAL_EINVAL, "coordID '%d' outside valid range [0, %d]\n\n", coordID, DOPPLERCOORD_LAST - 1 );
  }

  if ( !DopplerCoordinates[coordID].name ) {
    XLAL_ERROR_NULL ( XLAL_EINVAL, "coordID '%d' has no associated name\n\n", coordID );
  }

  return ( DopplerCoordinates[coordID].name );

} /* XLALDopplerCoordinateName() */


/**
 * Provide a pointer to a static string containing the a descriptive
 * 'help-string' describing the coordinate DopplerCoordinateID
 */
const CHAR *
XLALDopplerCoordinateHelp ( DopplerCoordinateID coordID )
{
  if ( coordID >= DOPPLERCOORD_LAST ) {
    XLAL_ERROR_NULL ( XLAL_EINVAL, "coordID '%d' outside valid range [0, %d]\n\n", coordID, DOPPLERCOORD_LAST - 1 );
  }

  if ( !DopplerCoordinates[coordID].help ) {
    XLAL_ERROR_NULL ( XLAL_EINVAL, "coordID '%d' has no associated help text\n\n", coordID );
  }

  return ( DopplerCoordinates[coordID].help );

} /* XLALDopplerCoordinateHelp() */

/**
 * Return a string (allocated here) containing a full name - helpstring listing
 * for all doppler-coordinates DopplerCoordinateID allowed by UniversalDopplerMetric.c
 */
CHAR *
XLALDopplerCoordinateHelpAll ( void )
{
  CHAR *helpstr;
  const CHAR *name;
  const CHAR *help;
  UINT4 i, len, maxlen = 0;
  CHAR buf[512];
  CHAR fmt[512];

#define HEADER "Doppler-coordinate names and explanations:\n--------------------------------------------------\n"
  if ( (helpstr = XLALCalloc ( strlen(HEADER)+1, sizeof(CHAR) )) == NULL ) {
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  strcpy ( helpstr, HEADER );
  len = strlen(helpstr);

  /* get maximal field-length of coordinate names */
  for ( i = 0; i < DOPPLERCOORD_LAST; i ++ )
    {
      if ( (name = XLALDopplerCoordinateName (  (DopplerCoordinateID) i )) == NULL ) {
	XLAL_ERROR_NULL ( XLAL_EINVAL );
      }
      maxlen = MYMAX ( maxlen, strlen(name) );
      sprintf ( fmt, "%%-%ds: %%s\n", maxlen + 2 );
    }

  /* assemble help-lines */
  for ( i = 0; i < DOPPLERCOORD_LAST; i ++ )
    {
      if ( (name = XLALDopplerCoordinateName ( (DopplerCoordinateID) i )) == NULL ) {
	XLAL_ERROR_NULL ( XLAL_EINVAL );
      }
      if ( (help = XLALDopplerCoordinateHelp ( (DopplerCoordinateID) i )) == NULL ) {
	XLAL_ERROR_NULL ( XLAL_EINVAL );
      }

      snprintf ( buf, sizeof(buf) - 1, fmt, name, help );
      len += strlen ( buf ) + 1;
      if ( (helpstr = XLALRealloc ( helpstr, len )) == NULL ) {
	XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      helpstr = strcat ( helpstr, buf );

    } /* for i < DOPPLERCOORD_LAST */

  return ( helpstr );

} /* XLALDopplerCoordinateHelpAll() */

/**
 * Free a FmetricAtoms_t structure, allowing any pointers to be NULL
 */
void
XLALDestroyFmetricAtoms ( FmetricAtoms_t *atoms )
{
  if ( !atoms )
    return;

  if ( atoms->a_a_i ) gsl_vector_free ( atoms->a_a_i );
  if ( atoms->a_b_i ) gsl_vector_free ( atoms->a_b_i );
  if ( atoms->b_b_i ) gsl_vector_free ( atoms->b_b_i );

  if ( atoms->a_a_i_j ) gsl_matrix_free ( atoms->a_a_i_j );
  if ( atoms->a_b_i_j ) gsl_matrix_free ( atoms->a_b_i_j );
  if ( atoms->b_b_i_j ) gsl_matrix_free ( atoms->b_b_i_j );

  XLALFree ( atoms );

  return;

} /* XLALDestroyFmetricAtoms() */


/**
 * Compute the 'F-metric' gF_ij (and also gFav_ij, m1_ij, m2_ij, m3_ij)
 * from the given FmetricAtoms and the signal amplitude parameters.
 *
 */
static DopplerFstatMetric*
XLALComputeFmetricFromAtoms ( const FmetricAtoms_t *atoms, REAL8 cosi, REAL8 psi )
{
  DopplerFstatMetric *metric;		/* output matrix */

  UINT4 dim, i, j;			/* Doppler index counters */
  REAL8 A, B, C, D;			/* 'standard' antenna-pattern coefficients (gsl-integrated, though) */
  REAL8 alpha1, alpha2, alpha3, eta2, cos2psi, sin2psi;

  if ( !atoms ) {
    XLALPrintError ("%s: illegal NULL input.\n\n", __func__ );
    XLAL_ERROR_NULL (  XLAL_EINVAL );
  }

  if ( !atoms->a_a_i || !atoms->a_b_i || !atoms->b_b_i || !atoms->a_a_i_j || !atoms->a_b_i_j || !atoms->b_b_i_j ) {
    XLALPrintError ("%s: input Fmetric-atoms not fully allocated.\n\n", __func__ );
    XLAL_ERROR_NULL (  XLAL_EINVAL );
  }

  dim = atoms->a_a_i->size;

  /* allocate output metric structure */
  if ( (metric = XLALCalloc ( 1, sizeof(*metric) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( 1, %zu) failed.\n\n", __func__, sizeof(*metric) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  metric->gF_ij = gsl_matrix_calloc ( dim, dim );
  metric->gFav_ij = gsl_matrix_calloc ( dim, dim );
  metric->m1_ij = gsl_matrix_calloc ( dim, dim );
  metric->m2_ij = gsl_matrix_calloc ( dim, dim );
  metric->m3_ij = gsl_matrix_calloc ( dim, dim );

  if ( !metric->gF_ij || !metric->gFav_ij || !metric->m1_ij || !metric->m2_ij || !metric->m3_ij ) {
    XLALPrintError ("%s: failed to gsl_matrix_calloc(%d,%d) for gF_ij, gFav_ij, m1_ij, m2_ij, m3_ij\n\n", __func__, dim, dim );
    XLALDestroyDopplerFstatMetric ( metric );
    XLAL_ERROR_NULL (  XLAL_ENOMEM );
  }

  A = atoms->a_a;
  B = atoms->b_b;
  C = atoms->a_b;

  D = A * B - C * C;	/* determinant of [A, C; C, B] */

  /* get amplitude-parameter factors alpha_1, alpha_2, alpha_3 */
  eta2 = SQUARE ( cosi );
  cos2psi = cos ( 2.0 * psi );
  sin2psi = sin ( 2.0 * psi );
  alpha1 = 0.25 * SQUARE ( 1.0 + eta2 ) * SQUARE ( cos2psi ) + eta2 * SQUARE ( sin2psi );
  alpha2 = 0.25 * SQUARE ( 1.0 + eta2 ) * SQUARE ( sin2psi ) + eta2 * SQUARE ( cos2psi );
  alpha3 = 0.25 * SQUARE ( 1.0 - eta2 ) * sin2psi * cos2psi;

  metric->rho2 = alpha1 * A + alpha2 * B + 2.0 * alpha3 * C;

  metric->maxrelerr = atoms->maxrelerr;

  /* ---------- compute components of the metric ---------- */
  for ( i=0; i < dim; i ++ )
    {
      REAL8 a_a_i, b_b_i, a_b_i;

      a_a_i = gsl_vector_get ( atoms->a_a_i, i );
      a_b_i = gsl_vector_get ( atoms->a_b_i, i );
      b_b_i = gsl_vector_get ( atoms->b_b_i, i );

      for ( j = 0; j <= i; j ++ )
	{
	  REAL8 a_a_i_j, b_b_i_j, a_b_i_j;
	  REAL8 a_a_j, b_b_j, a_b_j;

	  REAL8 P1_ij, P2_ij, P3_ij;	/* ingredients to m_r_ij */
	  REAL8 Q1_ij, Q2_ij, Q3_ij;	/* ingredients to m_r_ij */
	  REAL8 gg;

	  a_a_j = gsl_vector_get ( atoms->a_a_i, j );
	  a_b_j = gsl_vector_get ( atoms->a_b_i, j );
	  b_b_j = gsl_vector_get ( atoms->b_b_i, j );

	  a_a_i_j = gsl_matrix_get ( atoms->a_a_i_j, i, j );
	  a_b_i_j = gsl_matrix_get ( atoms->a_b_i_j, i, j );
	  b_b_i_j = gsl_matrix_get ( atoms->b_b_i_j, i, j );

	  /* trivial assignments, see Eq.(76) in \cite Prix07 */
	  P1_ij = a_a_i_j;
	  P2_ij = b_b_i_j;
	  P3_ij = a_b_i_j;

	  /* bit more involved, see Eq.(80)-(82) in \cite Prix07 [includes *explicit* index-symmetrization!!] */
	  Q1_ij = A * a_b_i * a_b_j + B * a_a_i * a_a_j - C * ( a_a_i * a_b_j + a_a_j * a_b_i );	/* (80) symmetrized */
	  Q1_ij /= D;

	  Q2_ij = A * b_b_i * b_b_j + B * a_b_i * a_b_j - C * ( a_b_i * b_b_j + a_b_j * b_b_i );	/* (81) symmetrized */
	  Q2_ij /= D;

	  Q3_ij = 0.5 * A * ( a_b_i * b_b_j + a_b_j * b_b_i )
	    + 0.5 * B * ( a_b_i * a_a_j + a_b_j * a_a_i )
	    - 0.5 * C * ( b_b_i * a_a_j + b_b_j * a_a_i + 2.0 * a_b_i * a_b_j );	/* (83) symmetrized */
	  Q3_ij /= D;

	  /* put the pieces together to compute m1_ij, m2_ij and m3_ij according to (85) */
	  gg = P1_ij  - Q1_ij;
	  gsl_matrix_set (metric->m1_ij, i, j, gg);
	  gsl_matrix_set (metric->m1_ij, j, i, gg);

	  gg = P2_ij  - Q2_ij;
	  gsl_matrix_set (metric->m2_ij, i, j, gg);
	  gsl_matrix_set (metric->m2_ij, j, i, gg);

	  gg = P3_ij  - Q3_ij;
	  gsl_matrix_set (metric->m3_ij, i, j, gg);
	  gsl_matrix_set (metric->m3_ij, j, i, gg);


	  /* assemble the 'full' F-stat metric from these ingredients, see Eq.(87) */
	  gg = alpha1 * gsl_matrix_get (metric->m1_ij, i, j ) + alpha2 * gsl_matrix_get (metric->m2_ij, i, j )
	    + 2.0 * alpha3 * gsl_matrix_get (metric->m3_ij, i, j );
	  gg /= metric->rho2;

	  gsl_matrix_set (metric->gF_ij, i, j, gg);
	  gsl_matrix_set (metric->gF_ij, j, i, gg);

	  /* compute 'average' F-stat metric as given by Eq.(93) */
	  gg = B * gsl_matrix_get (metric->m1_ij, i, j ) + A * gsl_matrix_get (metric->m2_ij, i, j )
	    - 2.0 * C * gsl_matrix_get (metric->m3_ij, i, j );
	  gg /= 2.0 * D;

	  gsl_matrix_set (metric->gFav_ij, i, j, gg);
	  gsl_matrix_set (metric->gFav_ij, j, i, gg);

	} /* for j <= i */

    } /* for i < dim */

  return metric;

} /* XLALComputeFmetricFromAtoms() */


/**
 * Function to compute *full* 4+n dimensional Fisher matric for the
 * full CW parameter-space of Amplitude + Doppler parameters !
 */
static gsl_matrix*
XLALComputeFisherFromAtoms ( const FmetricAtoms_t *atoms, PulsarAmplitudeParams Amp )
{
  gsl_matrix *fisher = NULL;	/* output matrix */

  UINT4 dimDoppler, dimFull, i, j;
  REAL8 al1, al2, al3;

  /* check input consistency */
  if ( !atoms ) {
    XLALPrintError ("%s: illegal NULL input.\n\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  if ( !atoms->a_a_i || !atoms->a_b_i || !atoms->b_b_i || !atoms->a_a_i_j || !atoms->a_b_i_j || !atoms->b_b_i_j ) {
    XLALPrintError ("%s: input Fmetric-atoms not fully allocated.\n\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  REAL8 A1,A2,A3,A4;
  {
    PulsarAmplitudeVect Amu;
    if ( XLALAmplitudeParams2Vect ( Amu, Amp ) != XLAL_SUCCESS ) {
      XLALPrintError ( "%s: XLALAmplitudeParams2Vect() failed with errno = %d\n\n", __func__, xlalErrno );
      XLAL_ERROR_NULL ( XLAL_EFUNC );
    }

    A1 = Amu[0];
    A2 = Amu[1];
    A3 = Amu[2];
    A4 = Amu[3];
  }

  dimDoppler = atoms->a_a_i->size;
  dimFull = 4 + dimDoppler;	/* 4 amplitude params + n Doppler params */

  if ( (fisher = gsl_matrix_calloc ( dimFull, dimFull )) == NULL ) {
    XLALPrintError ("%s: gsl_matric_calloc(%d,%d) failed.\n\n", __func__, dimFull, dimFull );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  /* ----- set pure Amplitude block 4x4: M_mu_nu ---------- */
  {
    REAL8 A, B, C;
    A = atoms->a_a;
    B = atoms->b_b;
    C = atoms->a_b;

    gsl_matrix_set ( fisher, 0, 0, A );
    gsl_matrix_set ( fisher, 2, 2, A );

    gsl_matrix_set ( fisher, 1, 1, B );
    gsl_matrix_set ( fisher, 3, 3, B );

    gsl_matrix_set ( fisher, 0, 1, C );
    gsl_matrix_set ( fisher, 1, 0, C );
    gsl_matrix_set ( fisher, 2, 3, C );
    gsl_matrix_set ( fisher, 3, 2, C );
  } /* amplitude-param block M_mu_nu */


  /* ----- set Doppler block (4+i,4+j) ---------- */
  al1 = SQUARE(A1) + SQUARE(A3);
  al2 = A1 * A2 + A3 * A4;
  al3 = SQUARE(A2) + SQUARE(A4);

  for ( i=0; i < dimDoppler; i ++ )
    {
      for ( j=0; j <= i; j ++ )
	{
	  REAL8 gg, a_a_i_j, a_b_i_j, b_b_i_j;

	  a_a_i_j = gsl_matrix_get(atoms->a_a_i_j, i, j);
	  a_b_i_j = gsl_matrix_get(atoms->a_b_i_j, i, j);
	  b_b_i_j = gsl_matrix_get(atoms->b_b_i_j, i, j);

	  gg = al1 * a_a_i_j + 2.0 * al2 * a_b_i_j + al3 * b_b_i_j;

	  gsl_matrix_set ( fisher, 4 + i, 4 + j, gg );
	  gsl_matrix_set ( fisher, 4 + j, 4 + i, gg );

	} /* for j <= i */
    } /* for i < dimDoppler */


  /* ----- compute mixed Amplitude-Doppler block ( non-square ) */
  for ( i=0; i < dimDoppler; i ++ )
    {
      REAL8 a_a_i, a_b_i, b_b_i;
      REAL8 AR[4];
      UINT4 mu;

      a_a_i = gsl_vector_get ( atoms->a_a_i, i );
      a_b_i = gsl_vector_get ( atoms->a_b_i, i );
      b_b_i = gsl_vector_get ( atoms->b_b_i, i );

      AR[0] =  A3 * a_a_i + A4 * a_b_i;
      AR[1] =  A3 * a_b_i + A4 * b_b_i;
      AR[2] = -A1 * a_a_i - A2 * a_b_i;
      AR[3] = -A1 * a_b_i - A2 * b_b_i;

      for ( mu = 0; mu < 4; mu ++ )
	{
	  gsl_matrix_set ( fisher, mu, 4 + i, AR[mu] );
	  gsl_matrix_set ( fisher, 4 + i, mu, AR[mu] );
	}

    } /* for i < dimDoppler */


  return fisher;

} /* XLALComputeFisherFromAtoms() */


/** Convert 3-D vector from equatorial into ecliptic coordinates */
void
XLALequatorialVect2ecliptic ( vect3D_t out, const vect3D_t in )
{
  static mat33_t rotEqu2Ecl = { { 1.0,        0,       0 },
                                { 0.0,  LAL_COSIEARTH, LAL_SINIEARTH },
                                { 0.0, -LAL_SINIEARTH, LAL_COSIEARTH } };

  XLALmatrix33_in_vect3 ( out, rotEqu2Ecl, in );

} /* XLALequatorialVect2ecliptic() */

/** Convert 3-D vector from ecliptic into equatorial coordinates */
void
XLALeclipticVect2equatorial ( vect3D_t out, const vect3D_t in )
{
  static mat33_t rotEcl2Equ =  { { 1.0,        0,       0 },
                                 { 0.0,  LAL_COSIEARTH, -LAL_SINIEARTH },
                                 { 0.0,  LAL_SINIEARTH,  LAL_COSIEARTH } };

  XLALmatrix33_in_vect3 ( out, rotEcl2Equ, in );

} /* XLALeclipticVect2equatorial() */

/** compute matrix product mat . vect */
void
XLALmatrix33_in_vect3 ( vect3D_t out, mat33_t mat, const vect3D_t in )
{

  UINT4 i,j;
  for ( i=0; i < 3; i ++ )
    {
      out[i] = 0;
      for ( j=0; j < 3; j ++ )
        {
          out[i] += mat[i][j] * in[j];
        }
    }

} /* XLALmatrix33_in_vect3() */

/**
 * Compute time-derivatives up to 'maxorder' of the Earths' orbital position vector
 * \f$r_{\mathrm{orb}}(t)\f$.
 *
 * Algorithm: using 5-point differentiation expressions on r_orb(t) returned from LALBarycenterEarth().
 *
 * Returns a vector of derivatives \f$\frac{d^n\,r_{\mathrm{orb}}}{d\,t^n}\f$ at the given
 * GPS time. Note, the return vector includes the zeroth-order derivative, so we return
 * (maxorder + 1) derivatives: n = 0 ... maxorder
 *
 */
vect3Dlist_t *
XLALComputeOrbitalDerivatives ( UINT4 maxorder,			/**< [in] highest derivative-order to compute */
                                const LIGOTimeGPS *tGPS,	/**< [in] GPS time at which to compute the derivatives */
                                const EphemerisData *edat	/**< [in] ephemeris data */
                                )
{
  EarthState earth;
  LIGOTimeGPS ti;
  REAL8 h = 0.5 * 86400.0;	/* finite-differencing step-size for rOrb. Before CAREFUL before changing this! */
  vect3D_t r0m2h, r0mh, r0, r0_h, r0_2h;
  vect3Dlist_t *ret = NULL;

  /* check input consistency */
  if ( maxorder > MAX_SPDNORDER ) {
    XLALPrintError ("%s: maxorder = %d too large, currently supports only up to maxorder = %d.\n", __func__, maxorder, MAX_SPDNORDER );
    XLAL_ERROR_NULL ( XLAL_EDOM );
  }

  if ( !tGPS ) {
    XLALPrintError ("%s: invalid NULL pointer received for 'tGPS'.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( !edat ) {
    XLALPrintError ("%s: invalid NULL pointer received for 'edat'.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }


  /* ----- find Earth's position at the 5 points: t0 -2h, t0-h, t0, t0+h, t0 + 2h ----- */

  /* t = t0 */
  ti = (*tGPS);
  if ( XLALBarycenterEarth( &earth, &ti, edat ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: call to XLALBarycenterEarth() failed!\n\n", __func__);
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  COPY_VECT ( r0, earth.posNow );

  /* t = t0 - h*/
  ti.gpsSeconds = (*tGPS).gpsSeconds - h;
  if ( XLALBarycenterEarth( &earth, &ti, edat ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: call to XLALBarycenterEarth() failed!\n\n", __func__);
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  COPY_VECT ( r0mh, earth.posNow );

  /* t = t0 - 2h*/
  ti.gpsSeconds = (*tGPS).gpsSeconds - 2 * h;
  if ( XLALBarycenterEarth( &earth, &ti, edat ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: call to XLALBarycenterEarth() failed!\n\n", __func__);
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  COPY_VECT ( r0m2h, earth.posNow );

  /* t = t0 + h*/
  ti.gpsSeconds = (*tGPS).gpsSeconds + h;
  if ( XLALBarycenterEarth( &earth, &ti, edat ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: call to XLALBarycenterEarth() failed!\n\n", __func__);
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  COPY_VECT ( r0_h, earth.posNow );

  /* t = t0 + 2h*/
  ti.gpsSeconds = (*tGPS).gpsSeconds + 2 * h;
  if ( XLALBarycenterEarth( &earth, &ti, edat ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: call to XLALBarycenterEarth() failed!\n\n", __func__);
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  COPY_VECT ( r0_2h, earth.posNow );

  /* use these 5 points to estimate derivatives */
  UINT4 i;
  vect3D_t rdn[MAX_SPDNORDER+1];
  COPY_VECT ( rdn[0], r0 );	// 0th order is imply r_orb(t0)
  for ( i=0; i < 3; i ++ )
    {
      rdn[1][i] = DERIV5P_1(r0m2h[i], r0mh[i], r0[i], r0_h[i], r0_2h[i], h );
      rdn[2][i] = DERIV5P_2(r0m2h[i], r0mh[i], r0[i], r0_h[i], r0_2h[i], h );
      rdn[3][i] = DERIV5P_3(r0m2h[i], r0mh[i], r0[i], r0_h[i], r0_2h[i], h );
      rdn[4][i] = DERIV5P_4(r0m2h[i], r0mh[i], r0[i], r0_h[i], r0_2h[i], h );
    } /* for i < 3 */

  /* allocate return list */
  if ( (ret = XLALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc(1,%zu)\n", __func__, sizeof(*ret) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  ret->length = maxorder + 1;
  if ( (ret->data = XLALCalloc ( ret->length, sizeof(*ret->data) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc(%d,%zu)\n", __func__, ret->length, sizeof(*ret->data) );
    XLALFree ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  UINT4 n;
  for ( n=0; n <= maxorder; n ++ )
    {
      COPY_VECT ( ret->data[n], rdn[n] );
    }

  return ret;

} /* XLALComputeOrbitalDerivatives() */

void
XLALDestroyVect3Dlist ( vect3Dlist_t *list )
{
  if ( !list )
    return;

  if ( list->data )
    XLALFree ( list->data );

  XLALFree ( list );

  return;

} /* XLALDestroyVect3Dlist() */

/**
 * Return the highest 'global-correlation' spindown order found in this coordinate system.
 * Counting nu0 = order1, nu1 = order2, nu2 = order3, ...,
 * order = 0 therefore means there are no GC spin coordinates at all
 */
static UINT4
findHighestGCSpinOrder ( const DopplerCoordinateSystem *coordSys )
{

  UINT4 maxorder = 0;
  UINT4 i;

  for ( i=0; i < coordSys->dim; i ++ )
    {
      UINT4 order = 0;
      if ( coordSys->coordIDs[i] ==  DOPPLERCOORD_GC_NU0 ) order = 1;
      if ( coordSys->coordIDs[i] ==  DOPPLERCOORD_GC_NU1 ) order = 2;
      if ( coordSys->coordIDs[i] ==  DOPPLERCOORD_GC_NU2 ) order = 3;
      if ( coordSys->coordIDs[i] ==  DOPPLERCOORD_GC_NU3 ) order = 4;
      maxorder = MYMAX ( maxorder, order );
    }

  return maxorder;
} /*  findHighestSpinOrder() */
