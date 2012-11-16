/*
 * Copyright (C) 2012 Karl Wette
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

/**
 * \file
 *
 * \author Reinhard Prix, Karl Wette
 * \ingroup PulsarMetric
 * \brief Function to compute the full F-statistic metric, including
 *  antenna-pattern functions from multi-detector, as derived in \ref Prix07.
 *
 *
 */

/*---------- INCLUDES ----------*/
#include <math.h>

/* gsl includes */
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>

#include <lal/FlatPulsarMetric.h>
#include <lal/PulsarTimes.h>
#include <lal/ComputeFstat.h>
#include <lal/XLALGSL.h>
#include <lal/Factorial.h>

#include <lal/UniversalDopplerMetric.h>

/*---------- SCALING FACTORS ----------*/

/** shortcuts for integer powers */
#define POW2(a)  ( (a) * (a) )
#define POW3(a)  ( (a) * (a) * (a) )
#define POW4(a)  ( (a) * (a) * (a) * (a) )
#define POW5(a)  ( (a) * (a) * (a) * (a) * (a) )

/* These constants define an internal scale used within CWPhaseDeriv_i().
 * Distances are normalised by SCALE_R, and times are normalised by SCALE_T.
 */
#define SCALE_R (LAL_AU_SI)
#define SCALE_T (LAL_YRSID_SI)

/*---------- LOOKUP ARRAYS ----------*/

/** Array of symbolic 'names' for various detector-motions
 */
static const CHAR *const DetectorMotionNames[DETMOTION_LAST] = {
  [DETMOTION_SPIN_ORBIT] = "spin+orbit",
  [DETMOTION_ORBIT] = "orbit",
  [DETMOTION_SPIN] = "spin",

  [DETMOTION_SPIN_PTOLEORBIT] = "spin+ptoleorbit",
  [DETMOTION_PTOLEORBIT] = "ptoleorbit",

  [DETMOTION_ORBIT_SPINZ] = "orbit+spin_Z",
  [DETMOTION_ORBIT_SPINXY] = "orbit+spin_XY",
};

/** Array of descriptor structs for each Doppler coordinate name
 */
const struct {
  const char *const name;	/**< coordinate name */
  const double scale;		/**< multiplicative scaling factor of the coordinate */
  const char *const help;	/**< help string explaining the coordinate's meaning and units */
} DopplerCoordinates[DOPPLERCOORD_LAST] = {

  [DOPPLERCOORD_FREQ] = {"freq", SCALE_T, "Frequency [Units: Hz]."},
  [DOPPLERCOORD_F1DOT] = {"f1dot", POW2(SCALE_T), "First spindown [Units: Hz/s]."},
  [DOPPLERCOORD_F2DOT] = {"f2dot", POW3(SCALE_T), "Second spindown [Units: Hz/s^2]."},
  [DOPPLERCOORD_F3DOT] = {"f3dot", POW4(SCALE_T), "Third spindown [Units: Hz/s^3]."},

  [DOPPLERCOORD_GC_NU0] = {"gc_nu0", SCALE_T, "Global correlation frequency [Units: Hz]. Activates 'reduced' detector position."},
  [DOPPLERCOORD_GC_NU1] = {"gc_nu1", POW2(SCALE_T), "Global correlation first spindown [Units: Hz/s]. Activates 'reduced' detector position."},
  [DOPPLERCOORD_GC_NU2] = {"gc_nu2", POW3(SCALE_T), "Global correlation second spindown [Units: Hz/s^2]. Activates 'reduced' detector position."},
  [DOPPLERCOORD_GC_NU3] = {"gc_nu3", POW4(SCALE_T), "Global correlation third spindown [Units: Hz/s^3]. Activates 'reduced' detector position."},

  [DOPPLERCOORD_ALPHA] = {"alpha", SCALE_R/LAL_C_SI, "Right ascension [Units: radians]. Uses 'reduced' detector position."},
  [DOPPLERCOORD_DELTA] = {"delta", SCALE_R/LAL_C_SI, "Declination [Units: radians]. Uses 'reduced' detector position."},

  [DOPPLERCOORD_N2X_EQU] = {"n2x_equ", SCALE_R/LAL_C_SI, "X component of contrained sky position in equatorial coordinates [Units: none]. Uses 'reduced' detector position."},
  [DOPPLERCOORD_N2Y_EQU] = {"n2y_equ", SCALE_R/LAL_C_SI, "Y component of contrained sky position in equatorial coordinates [Units: none]. Uses 'reduced' detector position."},

  [DOPPLERCOORD_N2X_ECL] = {"n2x_ecl", SCALE_R/LAL_C_SI, "X component of contrained sky position in ecliptic coordinates [Units: none]. Uses 'reduced' detector position."},
  [DOPPLERCOORD_N2Y_ECL] = {"n2y_ecl", SCALE_R/LAL_C_SI, "Y component of contrained sky position in ecliptic coordinates [Units: none]. Uses 'reduced' detector position."},

  [DOPPLERCOORD_N3X_EQU] = {"n3x_equ", SCALE_R/LAL_C_SI, "X component of unconstrained super-sky position in equatorial coordinates [Units: none]."},
  [DOPPLERCOORD_N3Y_EQU] = {"n3y_equ", SCALE_R/LAL_C_SI, "Y component of unconstrained super-sky position in equatorial coordinates [Units: none]."},
  [DOPPLERCOORD_N3Z_EQU] = {"n3z_equ", SCALE_R/LAL_C_SI, "Z component of unconstrained super-sky position in equatorial coordinates [Units: none]."},

  [DOPPLERCOORD_N3X_ECL] = {"n3x_ecl", SCALE_R/LAL_C_SI, "X component of unconstrained super-sky position in ecliptic coordinates [Units: none]."},
  [DOPPLERCOORD_N3Y_ECL] = {"n3y_ecl", SCALE_R/LAL_C_SI, "Y component of unconstrained super-sky position in ecliptic coordinates [Units: none]."},
  [DOPPLERCOORD_N3Z_ECL] = {"n3z_ecl", SCALE_R/LAL_C_SI, "Z component of unconstrained super-sky position in ecliptic coordinates [Units: none]."},

  [DOPPLERCOORD_N3SX_EQU] = {"n3sx_equ", SCALE_R/LAL_C_SI, "X spin-component of unconstrained super-sky position in equatorial coordinates [Units: none]."},
  [DOPPLERCOORD_N3SY_EQU] = {"n3sy_equ", SCALE_R/LAL_C_SI, "Y spin-component of unconstrained super-sky position in equatorial coordinates [Units: none]."},

  [DOPPLERCOORD_N3OX_ECL] = {"n3ox_ecl", SCALE_R/LAL_C_SI, "X orbit-component of unconstrained super-sky position in equatorial coordinates [Units: none]."},
  [DOPPLERCOORD_N3OY_ECL] = {"n3oy_ecl", SCALE_R/LAL_C_SI, "Y orbit-component of unconstrained super-sky position in equatorial coordinates [Units: none]."},

};

/*---------- DEFINES ----------*/
#define TRUE  (1==1)
#define FALSE (0==1)

/* get the scale associated with Doppler coordinate ID */
#define GET_SCALE(ID) ( (ID) == DOPPLERCOORD_NONE ? 1.0 : DopplerCoordinates[ID].scale )

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

/** components of antenna-pattern function: q_l = {a(t), b(t)}
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
  DetectorMotionType detMotionType;	/**< which detector-motion to use in metric integration */
  DopplerCoordinateID deriv1, deriv2;	/**< the two components of the derivative-product Phi_i_Phi_j to compute*/
  DopplerCoordinateID deriv;		/**< component for single phase-derivative Phi_i compute */
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


/*---------- empty initializers ---------- */
static const LALStatus empty_status;
static const EmissionTime empty_EmissionTime;
static const intparams_t empty_intparams;
static const PulsarTimesParamStruc empty_PulsarTimesParamStruc;

const PosVel3D_t empty_PosVel3D_t;
const DopplerMetricParams empty_DopplerMetricParams;
const DopplerCoordinateSystem empty_DopplerCoordinateSystem;
const MultiDetectorInfo empty_MultiDetectorInfo;

/*---------- Global variables ----------*/

BOOLEAN outputIntegrand = 0;

/* Some local constants. */
#define rOrb_c  (LAL_AU_SI / LAL_C_SI)
#define vOrb_c  (LAL_TWOPI * LAL_AU_SI / LAL_C_SI / LAL_YRSID_SI)

/*---------- internal prototypes ----------*/
DopplerMetric* XLALComputeFmetricFromAtoms ( const FmetricAtoms_t *atoms, REAL8 cosi, REAL8 psi );
gsl_matrix* XLALComputeFisherFromAtoms ( const FmetricAtoms_t *atoms, PulsarAmplitudeParams Amp );

double CW_am1_am2_Phi_i_Phi_j ( double tt, void *params );
double CWPhaseDeriv_i ( double tt, void *params );

double XLALAverage_am1_am2_Phi_i_Phi_j ( const intparams_t *params, double *relerr_max );
double CWPhase_cov_Phi_ij ( const MultiDetectorInfo *detInfo, const intparams_t *params, double *relerr_max );

int XLALPtolemaicPosVel ( PosVel3D_t *posvel, const LIGOTimeGPS *tGPS );
gsl_matrix *XLALProjectMetric ( const gsl_matrix * g_ij, const UINT4 c );

void equatorialVect2ecliptic ( vect3D_t out, const vect3D_t in );
void eclipticVect2equatorial ( vect3D_t out, const vect3D_t in );
void matrix33_in_vect3 ( vect3D_t out, mat33_t mat, const vect3D_t in );

UINT4 findHighestGCSpinOrder ( const DopplerCoordinateSystem *coordSys );

/*==================== FUNCTION DEFINITIONS ====================*/


/** Integrate a general quadruple product CW_am1_am2_Phi_i_Phi_j() from 0 to 1.
 * This implements the expression \f$\langle<q_1 q_2 \phi_i \phi_j\rangle\f$
 * for single-IFO average over the observation time.
 *
 * The input parameters correspond to CW_am1_am2_Phi_i_Phi_j()
 */
double
XLALAverage_am1_am2_Phi_i_Phi_j ( const intparams_t *params, double *relerr_max )
{
  intparams_t par = (*params);	/* struct-copy, as the 'deriv' field has to be changeable */
  gsl_function integrand;
  double epsrel = 1e-4;
  /* NOTE: this level of accuracy should be compatible with AM-coefficients involved
   * which are computed in REAL4 precision. We therefor cannot go lower than this it seems,
   * otherwise the gsl-integration fails to converge in some cases.
   */
  double epsabs = 0;
  const size_t limit = 64;
  gsl_integration_workspace *wksp = NULL;
  int stat;

  integrand.params = (void*)&par;

  /* compute <q_1 q_2 phi_i phi_j> as an integral from tt=0 to tt=1 */

  /* NOTE: this numerical integration runs into problems when integrating over
   * several days (~O(5d)), as the integrands are oscillatory functions on order of ~1/4d
   * and convergence degrades.
   * As a solution, we split the integral into N segments of 1/4 day duration, and compute
   * the final integral as a sum over partial integrals
   */
  REAL8 Tseg = 0.25 * LAL_DAYSID_SI;
  UINT4 Nseg = (UINT4) ceil ( params->Tspan / Tseg );
  UINT4 n;
  REAL8 dT = 1.0 / Nseg;

  REAL8 res = 0;
  REAL8 abserr2 = 0;

  /* allocate workspace for adaptive integration */
  wksp = gsl_integration_workspace_alloc(limit);
  XLAL_CHECK_REAL8(wksp != NULL, XLAL_EFAULT);

  const double scale12 = GET_SCALE(par.deriv1) * GET_SCALE(par.deriv2);

  integrand.function = &CW_am1_am2_Phi_i_Phi_j;
  for (n=0; n < Nseg; n ++ )
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
          XLALPrintWarning ("Segment n=%d, Result = %g, abserr=%g ==> relerr = %.2e > %.2e\n", n, res_n, err_n, fabs(err_n/res_n), epsrel);
          /* XLAL_ERROR_REAL8( XLAL_EFUNC ); */
        }

      res_n *= scale12;
      err_n *= scale12;

      res += res_n;
      abserr2 += SQUARE ( err_n );

    } /* for i < Nseg */

  REAL8 relerr = RELERR( sqrt(abserr2), fabs(res) );
  if ( relerr_max )
    (*relerr_max) = relerr;

  gsl_integration_workspace_free(wksp);

  return res;

} /* XLALAverage_am1_am2_Phi_i_Phi_j() */


/** For gsl-integration: general quadruple product between two antenna-pattern functions
 * am1, am2 in {a(t),b(t)} and two phase-derivatives phi_i * phi_j,
 * i.e. compute an expression of the form
 * \f$q_1(t) q_2(t) \phi_i(t) \phi_j(t)\f$, where \f$q_i = \{a(t), b(t)\}\f$.
 *
 * NOTE: this can be 'truncated' to any sub-expression by using
 * AMCOMP_NONE for antenna-pattern component and DOPPLERCOORD_NONE for DopplerCoordinate,
 * eg in this way this function can be used to compute \f$a^2(t), b^2(t), a(t) b(t)\f$,
 * or \f$phi_i(t) phi_j(t)\f$.
 */
double
CW_am1_am2_Phi_i_Phi_j ( double tt, void *params )
{
  intparams_t *par = (intparams_t*) params;

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
  if ( par->deriv1 != DOPPLERCOORD_NONE )
    {
      par->deriv = par->deriv1;
      phi_i = CWPhaseDeriv_i ( tt, params );
    }
  else
    phi_i = 1.0;

  /* second Doppler phase derivative */
  if ( par->deriv2 != DOPPLERCOORD_NONE )
    {
      par->deriv = par->deriv2;
      phi_j = CWPhaseDeriv_i ( tt, params );
    }
  else
    phi_j = 1.0;

  ret = am1 * am2 * phi_i * phi_j;

  if ( outputIntegrand )
    {
      printf ( "%d %d %d %d   %f  %f  %f  %f  %f   %f\n",
               par->amcomp1, par->amcomp2, par->deriv1, par->deriv2,
               tt, am1, am2, phi_i, phi_j, ret );
    }


  return ( ret );

} /* CW_am1_am2_Phi_i_Phi_j() */



/** Partial derivative of continuous-wave (CW) phase, with respect
 * to Doppler coordinate 'i' := intparams_t->phderiv
 *
 * Time is in 'natural units' of Tspan, i.e. tt is in [0, 1] corresponding
 * to GPS-times in [startTime, startTime + Tspan ]
 *
 */
double
CWPhaseDeriv_i ( double tt, void *params )
{
  REAL8 ret = 0;
  intparams_t *par = (intparams_t*) params;
  vect3D_t nn_equ, nn_ecl;	/* skypos unit vector */
  vect3D_t nDeriv_i;	/* derivative of sky-pos vector wrt i */

  /* positions/velocities at time tt: */
  PosVel3D_t spin_posvel = empty_PosVel3D_t;
  PosVel3D_t orbit_posvel = empty_PosVel3D_t;
  PosVel3D_t posvel = empty_PosVel3D_t;

  /* orbit position in ecliptic plane */
  vect3D_t ecl_pos = empty_vect3D_t;
  vect3D_t ecl_orbit_pos = empty_vect3D_t;

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
  equatorialVect2ecliptic ( nn_ecl, nn_equ );

  /* get current detector position r(t) and velocity v(t) */
  REAL8 ttSI = par->startTime + tt * par->Tspan;	/* current GPS time in seconds */
  LIGOTimeGPS ttGPS;
  XLALGPSSetREAL8( &ttGPS, ttSI );
  int errnum;
  XLAL_TRY( XLALDetectorPosVel ( &spin_posvel, &orbit_posvel, &ttGPS, par->site, par->edat, par->detMotionType ), errnum );
  if ( errnum ) {
    XLALPrintError ( "%s: Call to XLALDetectorPosVel() failed!\n", __func__);
    GSL_ERROR_VAL( "Failure in CWPhaseDeriv_i", GSL_EFAILED, GSL_NAN );
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

  /* compute orbital detector positions projected onto ecliptic plane */
  equatorialVect2ecliptic(ecl_pos, posvel.pos);
  equatorialVect2ecliptic(ecl_orbit_pos, orbit_posvel.pos);

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
  equatorialVect2ecliptic ( rr_ord_Ecl, rr_ord_Equ );	  /* convert into ecliptic coordinates */

  /* now compute the requested phase derivative */
  switch ( par->deriv )
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

    case DOPPLERCOORD_N2X_EQU:		/**< X component of contrained sky position in equatorial coordinates [Units: none]. Uses 'reduced' detector position. */
      ret = LAL_TWOPI * Freq * ( rr_ord_Equ[0] - (nn_equ[0]/nn_equ[2]) * rr_ord_Equ[2] );
      break;
    case DOPPLERCOORD_N2Y_EQU:		/**< Y component of contrained sky position in equatorial coordinates [Units: none]. Uses 'reduced' detector position. */
      ret = LAL_TWOPI * Freq * ( rr_ord_Equ[1] - (nn_equ[1]/nn_equ[2]) * rr_ord_Equ[2] );
      break;

    case DOPPLERCOORD_N2X_ECL:		/**< X component of contrained sky position in ecliptic coordinates [Units: none]. Uses 'reduced' detector position. */
      ret = LAL_TWOPI * Freq * ( rr_ord_Ecl[0] - (nn_ecl[0]/nn_ecl[2]) * rr_ord_Ecl[2] );
      break;
    case DOPPLERCOORD_N2Y_ECL:		/**< Y component of contrained sky position in ecliptic coordinates [Units: none]. Uses 'reduced' detector position. */
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

    case DOPPLERCOORD_N3OX_ECL:	/**< X orbit-component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * ecl_orbit_pos[0];
      break;
    case DOPPLERCOORD_N3OY_ECL:	/**< Y orbit-component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
      ret = LAL_TWOPI * Freq * ecl_orbit_pos[1];
      break;

    default:
      XLALPrintError("%s: Unknown phase-derivative type '%d'\n", __func__, par->deriv );
      GSL_ERROR_VAL( "Failure in CWPhaseDeriv_i", GSL_EFAILED, GSL_NAN );
      break;

    } /* switch par->deriv */

  return ret;

} /* CWPhaseDeriv_i() */


/** Given a GPS time and detector, return the current position (and velocity) of the detector.
 *
 * NOTE: the 'special' flag allows to simulate artifically truncated
 * detector motions, such as pure orbital motion.
 *
 */
int
XLALDetectorPosVel ( PosVel3D_t *spin_posvel,	/**< [out] instantaneous sidereal position and velocity vector */
                     PosVel3D_t *orbit_posvel,	/**< [out] instantaneous orbital position and velocity vector */
		     const LIGOTimeGPS *tGPS,	/**< [in] GPS time */
		     const LALDetector *site,	/**< [in] detector info */
		     const EphemerisData *edat,	/**< [in] ephemeris data */
		     DetectorMotionType special	/**< [in] detector motion type */
		     )
{
  EarthState earth;
  BarycenterInput baryinput = empty_BarycenterInput;
  EmissionTime emit = empty_EmissionTime;
  PosVel3D_t Det_wrt_Earth;
  PosVel3D_t PtoleOrbit;
  PosVel3D_t Spin_z, Spin_xy;
  REAL8 eZ[3];

  if ( !tGPS || !site || !edat ) {
    XLALPrintError ( "%s: Illegal NULL pointer passed!\n", __func__);
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* ----- find ephemeris-based position of Earth wrt to SSB at this moment */
  if ( XLALBarycenterEarth( &earth, tGPS, edat ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: call to XLALBarycenterEarth() failed!\n\n", __func__);
    XLAL_ERROR( XLAL_EFUNC );
  }
  /* ----- find ephemeris-based position of detector wrt to SSB */
  baryinput.tgps = *tGPS;
  baryinput.site = *site;
  baryinput.site.location[0] /= LAL_C_SI; baryinput.site.location[1] /= LAL_C_SI; baryinput.site.location[2] /= LAL_C_SI;
  baryinput.alpha = 0; baryinput.delta = 0; baryinput.dInv = 0;
  if ( XLALBarycenter ( &emit, &baryinput, &earth ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: call to XLALBarycenter() failed!\n\n", __func__);
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* ----- determine position-vector of detector wrt center of Earth */
  COPY_VECT(Det_wrt_Earth.pos, emit.rDetector);
  SUB_VECT(Det_wrt_Earth.pos, earth.posNow);

  COPY_VECT(Det_wrt_Earth.vel, emit.vDetector);
  SUB_VECT(Det_wrt_Earth.vel, earth.velNow);

  eZ[0] = 0; eZ[1] = -LAL_SINIEARTH; eZ[2] = LAL_COSIEARTH; 	/* ecliptic z-axis in equatorial coordinates */
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
  if ( (special == DETMOTION_SPIN_PTOLEORBIT) || (special == DETMOTION_PTOLEORBIT) )
    {
      if ( XLALPtolemaicPosVel ( &PtoleOrbit, tGPS ) ) {
	XLALPrintError ( "%s: call to XLALPtolemaicPosVel() failed!\n\n", __func__);
	XLAL_ERROR( XLAL_EFUNC );
      }
    }

  /* ----- return the requested type of detector motion */
  switch ( special )
    {
      /* full detector-motion: ephemeris-orbital + Earth-spin */
    case DETMOTION_SPIN_ORBIT:
      if ( spin_posvel ) {
        COPY_VECT(spin_posvel->pos, Det_wrt_Earth.pos);
        COPY_VECT(spin_posvel->vel, Det_wrt_Earth.vel);
      }

      if ( orbit_posvel ) {
        COPY_VECT(orbit_posvel->pos, earth.posNow);
        COPY_VECT(orbit_posvel->vel, earth.velNow);
      }
      break;

      /* full ephemeris orbital detector-motion, neglecting Earth-spin */
    case DETMOTION_ORBIT:
      if ( spin_posvel ) {
        ZERO_VECT(spin_posvel->pos);
        ZERO_VECT(spin_posvel->vel);
      }

      if ( orbit_posvel ) {
        COPY_VECT(orbit_posvel->pos, earth.posNow);
        COPY_VECT(orbit_posvel->vel, earth.velNow);
      }
      break;

      /* detector-motion including only Earth-spin, no orbital motion */
    case DETMOTION_SPIN:
      if ( spin_posvel ) {
        COPY_VECT(spin_posvel->pos, Det_wrt_Earth.pos);
        COPY_VECT(spin_posvel->vel, Det_wrt_Earth.vel);
      }

      if ( orbit_posvel ) {
        ZERO_VECT(orbit_posvel->pos);
        ZERO_VECT(orbit_posvel->vel);
      }
      break;

      /* pure orbital detector motion, using "Ptolemaic" (ie. circular) approximation */
    case DETMOTION_PTOLEORBIT:
      if ( spin_posvel ) {
        ZERO_VECT(spin_posvel->pos);
        ZERO_VECT(spin_posvel->vel);
      }

      if ( orbit_posvel ) {
        COPY_VECT(orbit_posvel->pos, PtoleOrbit.pos);
        COPY_VECT(orbit_posvel->vel, PtoleOrbit.vel);
      }
      break;

      /* Ptolemaic-orbital motion, plus Earth spin */
    case DETMOTION_SPIN_PTOLEORBIT:
      if ( spin_posvel ) {
        COPY_VECT(spin_posvel->pos, Det_wrt_Earth.pos);
        COPY_VECT(spin_posvel->vel, Det_wrt_Earth.vel);
      }

      if ( orbit_posvel ) {
        COPY_VECT(orbit_posvel->pos, PtoleOrbit.pos);
        COPY_VECT(orbit_posvel->vel, PtoleOrbit.vel);
      }
      /*
      printf ("\nPtole = [ %f, %f, %f ], Ephem = [%f, %f, %f]\n",
	      posvel->pos[0], posvel->pos[1], posvel->pos[2],
	      emit.rDetector[0], emit.rDetector[1], emit.rDetector[2] );
      */
      break;

      /**< orbital motion plus *only* z-component of Earth spin-motion wrt to ecliptic plane */
    case DETMOTION_ORBIT_SPINZ:
      if ( spin_posvel ) {
        COPY_VECT ( spin_posvel->pos, Spin_z.pos );
        COPY_VECT ( spin_posvel->vel, Spin_z.vel );
      }

      if ( orbit_posvel ) {
        COPY_VECT(orbit_posvel->pos, earth.posNow);
        COPY_VECT(orbit_posvel->vel, earth.velNow);
      }

      break;

      /**< orbital motion plus *only* x+y component of Earth spin-motion in the ecliptic */
    case DETMOTION_ORBIT_SPINXY:
      if ( spin_posvel ) {
        COPY_VECT ( spin_posvel->pos, Spin_xy.pos );
        COPY_VECT ( spin_posvel->vel, Spin_xy.vel );
      }

      if ( orbit_posvel ) {
        COPY_VECT(orbit_posvel->pos, earth.posNow);
        COPY_VECT(orbit_posvel->vel, earth.velNow);
      }

      break;

    default:
      XLALPrintError("\n%s: Illegal 'special' value passed: '%d'\n\n", __func__, special );
      XLAL_ERROR( XLAL_EINVAL );
      break;
    } /* switch(special) */


#if 0
  /* debug output */
  printf ("%.6f  %.16g  %.16g  %.16g    %.16g  %.16g %.16g \n",
          GPS2REAL8((*tGPS)),
          posvel->pos[0], posvel->pos[1], posvel->pos[2],
          posvel->vel[0], posvel->vel[1], posvel->vel[2] );

#endif

  return XLAL_SUCCESS;

} /* XLALDetectorPosVel() */



/** Compute position and velocity assuming a purely "Ptolemaic" orbital motion
 * (i.e. on a circle) around the sun, approximating Earth's orbit
 */
int
XLALPtolemaicPosVel ( PosVel3D_t *posvel,		/**< [out] instantaneous position and velocity vector */
		      const LIGOTimeGPS *tGPS		/**< [in] GPS time */
		      )
{
  PulsarTimesParamStruc times = empty_PulsarTimesParamStruc;
  REAL8 phiOrb;   /* Earth orbital revolution angle, in radians. */
  REAL8 sinOrb, cosOrb;
  LALStatus status = empty_status;

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



/** Compute a pure phase-deriv covariance \f$[\phi_i, \phi_j] = \langle phi_i phi_j\rangle - \langle phi_i\rangle\langle phi_j\rangle\f$
 * which gives a component of the "phase metric"
 */
double
CWPhase_cov_Phi_ij ( const MultiDetectorInfo *detInfo, const intparams_t *params, double* relerr_max )
{
  gsl_function integrand;

  intparams_t par = (*params);	/* struct-copy, as the 'deriv' field has to be changeable */
  int stat;
  double ret;

  /* sanity-check: don't allow any AM-coeffs being turned on here! */
  if ( par.amcomp1 != AMCOMP_NONE || par.amcomp2 != AMCOMP_NONE ) {
    XLALPrintError ( "%s: Illegal input, amcomp[12] must be set to AMCOMP_NONE!\n", __func__ );
    XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  integrand.params = (void*)&par;


  double epsrel = 1e-6;
  /* NOTE: this level of accuracy is only achievable *without* AM-coefficients involved
   * which are computed in REAL4 precision. For the current function this is OK, as this
   * function is only supposed to compute *pure* phase-derivate covariances.
   */

  /* NOTE: this numerical integration still runs into problems when integrating over
   * long durations (~O(23d)), as the integrands are oscillatory functions on order of ~1d
   * and convergence degrades.
   * As a solution, we split the integral into N segments of 1 day duration, and compute
   * the final integral as a sum over partial integrals
   */
  REAL8 Tseg = LAL_DAYSID_SI;
  UINT4 Nseg = (UINT4) ceil ( params->Tspan / Tseg );
  UINT4 n;
  REAL8 dT = 1.0 / Nseg;

  double epsabs = 1e-3; 	/* we need an abs-cutoff as well, as epsrel can be too restrictive for small integrals */
  double abserr, maxrelerr = 0;
  const size_t limit = 64;
  gsl_integration_workspace *wksp = NULL;
  double av_ij = 0, av_i = 0, av_j = 0;
  double av_ij_err = 0, av_i_err = 0, av_j_err = 0;

  /* allocate workspace for adaptive integration */
  wksp = gsl_integration_workspace_alloc(limit);
  XLAL_CHECK_REAL8(wksp != NULL, XLAL_EFAULT);

  const double scale1 = GET_SCALE(par.deriv1);
  const double scale2 = GET_SCALE(par.deriv2);
  const double scale12 = scale1 * scale2;

  // loop over detectors
  REAL8 total_weight = 0;
  for (UINT4 det = 0; det < detInfo->length; ++det) {

    // set detector for phase integrals
    par.site = &detInfo->sites[det];

    // accumulate detector weights
    const REAL8 weight = detInfo->detWeights[det];
    total_weight += weight;

    for (n=0; n < Nseg; n ++ ) {

      REAL8 ti = 1.0 * n * dT;
      REAL8 tf = MYMIN( (n+1.0) * dT, 1.0 );
      double res_n;

      /* compute <phi_i phi_j> */
      integrand.function = &CW_am1_am2_Phi_i_Phi_j;
      XLAL_CALLGSL ( stat = gsl_integration_qag (&integrand, ti, tf, epsabs, epsrel, limit, GSL_INTEG_GAUSS61, wksp, &res_n, &abserr) );
      if ( outputIntegrand ) printf ("\n");
      if ( stat != 0 ) {
        XLALPrintError ( "\n%s: GSL-integration 'gsl_integration_qag()' of <Phi_i Phi_j> failed! seg=%d, av_ij_n=%g, abserr=%g\n",
                         __func__, n, res_n, abserr);
        XLAL_ERROR_REAL8( XLAL_EFUNC );
      }
      res_n *= scale12;
      abserr *= scale12;
      av_ij += weight * res_n;
      av_ij_err += weight * SQUARE (abserr);

      /* compute <phi_i> */
      integrand.function = &CWPhaseDeriv_i;
      par.deriv = par.deriv1;
      XLAL_CALLGSL ( stat = gsl_integration_qag (&integrand, ti, tf, epsabs, epsrel, limit, GSL_INTEG_GAUSS61, wksp, &res_n, &abserr) );
      if ( stat != 0 ) {
        XLALPrintError ( "\n%s: GSL-integration 'gsl_integration_qag()' of <Phi_i> failed! seg=%d, av_i_n=%g, abserr=%g\n",
                         __func__, n, res_n, abserr);
        XLAL_ERROR_REAL8( XLAL_EFUNC );
      }
      res_n *= scale1;
      abserr *= scale1;
      av_i += weight * res_n;
      av_i_err += weight * SQUARE (abserr);

      /* compute <phi_j> */
      integrand.function = &CWPhaseDeriv_i;
      par.deriv = par.deriv2;
      XLAL_CALLGSL ( stat = gsl_integration_qag (&integrand, ti, tf, epsabs, epsrel, limit, GSL_INTEG_GAUSS61, wksp, &res_n, &abserr) );
      if ( stat != 0 ) {
        XLALPrintError ( "\n%s: GSL-integration 'gsl_integration_qag()' of <Phi_j> failed! seg=%d, av_j_n=%g, abserr=%g\n",
                         __func__, n, res_n, abserr);
        XLAL_ERROR_REAL8( XLAL_EFUNC );
      }
      res_n *= scale2;
      abserr *= scale2;
      av_j += weight * res_n;
      av_j_err += weight * SQUARE (abserr);

    } /* for i < Nseg */

  } // for det < detInfo->length

  // raise error if no detector weights were given
  if (total_weight == 0) {
    XLAL_ERROR( XLAL_EDOM, "Detectors weights are all zero!" );
  }

  // normalise by total weight
  const REAL8 inv_total_weight = 1.0 / total_weight;
  av_ij *= inv_total_weight;
  av_i *= inv_total_weight;
  av_j *= inv_total_weight;
  av_ij_err *= inv_total_weight;
  av_i_err *= inv_total_weight;
  av_j_err *= inv_total_weight;

  av_ij_err = RELERR( sqrt(av_ij_err), fabs(av_ij) );
  av_i_err = RELERR( sqrt(av_i_err), fabs (av_i) );
  av_j_err = RELERR( sqrt(av_j_err), fabs (av_j) );

  maxrelerr = MYMAX ( av_ij_err, av_i_err );
  maxrelerr = MYMAX ( maxrelerr, av_j_err );

  if ( relerr_max )
    (*relerr_max) = maxrelerr;

  ret = av_ij - av_i * av_j;

  gsl_integration_workspace_free(wksp);

  return ret;	/* return covariance */

} /* CWPhase_cov_Phi_ij() */


/** Calculate an approximate "phase-metric" with the specified parameters.
 *
 * Note: if this function is called with multiple detectors, the phase components
 * are averaged over detectors as well as time. This is a somewhat ad-hoc approach;
 * if you want a more rigorous multi-detector metric you need to use the full
 * Fstat-metric, as computed by XLALDopplerFstatMetric().
 *
 * Return NULL on error.
 */
gsl_matrix *
XLALDopplerPhaseMetric ( const DopplerMetricParams *metricParams,  	/**< input parameters determining the metric calculation */
			 const EphemerisData *edat,			/**< ephemeris data */
                         double *relerr_max				/**< [in] maximal relative error in integration */
			 )
{
  gsl_matrix *g_ij = NULL;
  intparams_t intparams = empty_intparams;
  UINT4 i, j;
  REAL8 gg;
  UINT4 dim;
  const LIGOTimeGPS *refTime, *startTime;
  const DopplerCoordinateSystem *coordSys;

  /* ---------- sanity/consistency checks ---------- */
  if ( !metricParams || !edat ) {
    XLALPrintError ("\n%s: Illegal NULL pointer passed.\n\n", __func__);
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }

  startTime = &(metricParams->startTime);
  refTime   = &(metricParams->signalParams.Doppler.refTime);

  dim = metricParams->coordSys.dim;
  coordSys = &(metricParams->coordSys);

  /* ---------- prepare output metric ---------- */
  if ( (g_ij = gsl_matrix_calloc ( dim, dim )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc(%d, %d) failed.\n\n", __func__, dim, dim );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  /* ---------- set up integration parameters ---------- */
  intparams.edat = edat;
  intparams.startTime = XLALGPSGetREAL8 ( startTime );
  intparams.refTime   = XLALGPSGetREAL8 ( refTime );
  intparams.Tspan     = metricParams->Tspan;
  intparams.dopplerPoint = &(metricParams->signalParams.Doppler);
  intparams.detMotionType = metricParams->detMotionType;
  intparams.site = NULL;
  intparams.approxPhase = metricParams->approxPhase;
  /* deactivate antenna-patterns for phase-metric */
  intparams.amcomp1 = AMCOMP_NONE;
  intparams.amcomp2 = AMCOMP_NONE;

  /* if using 'global correlation' frequency variables, determine the highest spindown order: */
  UINT4 maxorder = findHighestGCSpinOrder ( coordSys );

  /* compute rOrb(t) derivatives at reference time */
  if ( (intparams.rOrb_n = XLALComputeOrbitalDerivatives ( maxorder, &intparams.dopplerPoint->refTime, edat )) == NULL ) {
    XLALPrintError ("%s: XLALComputeOrbitalDerivatives() failed.\n", __func__);
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }

#if 0
  /* diagnostic / debug output */
  UINT4 n;
  printf (" rOrb_n(%d) = [ ", intparams.dopplerPoint->refTime.gpsSeconds );
  for ( n=0; n < intparams.rOrb_n->length; n ++)
    printf ("[%g, %g, %g]%s", intparams.rOrb_n->data[n][0], intparams.rOrb_n->data[n][1], intparams.rOrb_n->data[n][2],
            (n < intparams.rOrb_n->length -1 ) ? ", " : " ]\n" );
#endif

  /* ---------- compute components of the phase-metric ---------- */
  double maxrelerr = 0, err;
  for ( i=0; i < dim; i ++ )
    {
      for ( j = 0; j <= i; j ++ )
	{
	  /* g_ij */
	  intparams.deriv1 = coordSys->coordIDs[i];
	  intparams.deriv2 = coordSys->coordIDs[j];
	  gg = CWPhase_cov_Phi_ij ( &metricParams->detInfo, &intparams, &err );	/* [Phi_i, Phi_j] */
          maxrelerr = MYMAX ( maxrelerr, err );
	  if ( xlalErrno ) {
	    XLALPrintError ("\n%s: Integration of g_ij (i=%d, j=%d) failed. errno = %d\n", __func__, i, j, xlalErrno );
            xlalErrno = 0;
            BOOLEAN sav = outputIntegrand;
            outputIntegrand = 1;
            gg = CWPhase_cov_Phi_ij ( &metricParams->detInfo, &intparams, &maxrelerr );	/* [Phi_i, Phi_j] */
            outputIntegrand = sav;

	    XLAL_ERROR_NULL( XLAL_EFUNC );
	  }
	  gsl_matrix_set (g_ij, i, j, gg);
	  gsl_matrix_set (g_ij, j, i, gg);

	} /* for j <= i */

    } /* for i < dim */

  if ( relerr_max )
    (*relerr_max) = maxrelerr;

  /* free memory */
  XLALDestroyVect3Dlist ( intparams.rOrb_n );

  return g_ij;

} /* XLALDopplerPhaseMetric() */


/** Calculate the phase-metric, the *full* (multi-IFO) Fstat-metrix
 *  and the Fisher-matrix derived in \ref Prix07.
 *
 * Note: The returned DopplerMetric struct contains the matrices
 * g_ij (the phase metric), gF_ij (the F-metric), gFav_ij (the average F-metric),
 * m1_ij, m2_ij, m3_ij (auxiliary matrices)
 * and Fisher_ab (the full 4+n dimensional Fisher matrix).
 *
 * The returned metric struct also carries the meta-info about
 * the metrics in the field 'DopplerMetricParams meta'.
 *
 *
 * Return NULL on error.
 */
DopplerMetric *
XLALDopplerFstatMetric ( const DopplerMetricParams *metricParams,  	/**< input parameters determining the metric calculation */
			 const EphemerisData *edat			/**< ephemeris data */
			 )
{
  DopplerMetric *metric = NULL;
  REAL8 cosi, psi;
  double relerr;
  gsl_matrix *tmp;

  /* ---------- sanity/consistency checks ---------- */
  if ( !metricParams || !edat ) {
    XLALPrintError ("%s: Illegal NULL pointer passed!\n\n", __func__);
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }

  if ( metricParams->metricType >= METRIC_TYPE_LAST ) {
    XLALPrintError ("%s: Invalid value '%d' for metricType received. Must be within [%d,%d]!\n\n",
                    __func__, metricParams->metricType, 0, METRIC_TYPE_LAST - 1);
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }

  /* if we're asked to compute F-metric only (1) or F-metric + phase-metric (2) */
  if ( metricParams->metricType == METRIC_TYPE_FSTAT || metricParams->metricType == METRIC_TYPE_ALL )
    {
      FmetricAtoms_t *atoms = NULL;

      /* ---------- compute Fmetric 'atoms', ie the averaged <a^2>, <a b Phi_i>, <a^2 Phi_i Phi_j>, etc ---------- */
      if ( (atoms = XLALComputeAtomsForFmetric ( metricParams, edat )) == NULL ) {
        XLALPrintError ("%s: XLALComputeAtomsForFmetric() failed. errno = %d\n\n", __func__, xlalErrno );
        XLAL_ERROR_NULL( XLAL_EFUNC );
      }

      /* ----- compute the F-metric gF_ij and related matrices ---------- */
      cosi = metricParams->signalParams.Amp.cosi;
      psi  = metricParams->signalParams.Amp.psi;

      if ( (metric = XLALComputeFmetricFromAtoms ( atoms, cosi, psi)) == NULL ) {
        XLALPrintError ("%s: XLALComputeFmetricFromAtoms() failed, errno = %d\n\n", __func__, xlalErrno );
        XLALDestroyFmetricAtoms ( atoms );
        XLAL_ERROR_NULL( XLAL_EFUNC );
      }

      /* ----- compute the full 4+n dimensional Fisher matrix ---------- */
      if ( (metric->Fisher_ab = XLALComputeFisherFromAtoms ( atoms, metricParams->signalParams.Amp )) == NULL ) {
        XLALPrintError ("%s: XLALComputeFisherFromAtoms() failed. errno = %d\n\n", xlalErrno );
        XLALDestroyFmetricAtoms ( atoms );
        XLALDestroyDopplerMetric ( metric );
        XLAL_ERROR_NULL( XLAL_EFUNC );
      }

      XLALDestroyFmetricAtoms ( atoms );
    } /* if compute F-metric */

  if ( metricParams->metricType == METRIC_TYPE_PHASE || metricParams->metricType == METRIC_TYPE_ALL )
    {
      /* if return-container 'metric' hasn't been allocated already earlier, we do it now */
      if (!metric ) {
        if ( (metric = XLALCalloc ( 1, sizeof(*metric) )) == NULL ) {
          XLALPrintError ("%s: XLALCalloc ( 1, %d) failed.\n\n", sizeof(*metric) );
          XLAL_ERROR_NULL ( XLAL_ENOMEM );
        }
      }

      /* ----- compute the standard phase-metric g_ij ---------- */
      if ( (metric->g_ij = XLALDopplerPhaseMetric ( metricParams, edat, &relerr )) == NULL ) {
        XLALPrintError ("%s: XLALDopplerPhaseMetric() failed, errno = %d.\n\n", __func__, xlalErrno );
        XLALDestroyDopplerMetric ( metric );
        XLAL_ERROR_NULL( XLAL_EFUNC );
      }
      metric->maxrelerr_gPh = relerr;

    } /* if compute phase-metric */

  /* ----- if requested, project gF_ij, gFav_ij and g_ij onto coordinate 'projectCoord' */
  if ( metricParams->projectCoord >= 0 )
    {
      UINT4 projCoord = (UINT4)metricParams->projectCoord;

      /* gF_ij */
      if ( metric->gF_ij )
        {
          if ( (tmp = XLALProjectMetric ( metric->gF_ij, projCoord )) == NULL ) {
            XLALPrintError ("%s: failed to project gF_ij onto coordinate '%d'. errno=%d\n", __func__, projCoord, xlalErrno );
            XLAL_ERROR_NULL ( XLAL_EFUNC );
          }
          gsl_matrix_free ( metric->gF_ij );
          metric->gF_ij = tmp;
        }

      /* gFav_ij */
      if ( metric->gFav_ij )
        {
          if ( (tmp = XLALProjectMetric ( metric->gFav_ij, projCoord )) == NULL ) {
            XLALPrintError ("%s: failed to project gFav_ij onto coordinate '%d'. errno=%d\n", __func__, projCoord, xlalErrno );
            XLAL_ERROR_NULL ( XLAL_EFUNC );
          }
          gsl_matrix_free ( metric->gFav_ij );
          metric->gFav_ij = tmp;
        }

      /* phase-metric g_ij */
      if ( metric->g_ij )
        {
          if ( (tmp = XLALProjectMetric ( metric->g_ij, projCoord )) == NULL ) {
            XLALPrintError ("%s: failed to project g_ij onto coordinate '%d'. errno=%d\n", __func__, projCoord, xlalErrno );
            XLAL_ERROR_NULL ( XLAL_EFUNC );
          }
          gsl_matrix_free ( metric->g_ij );
          metric->g_ij = tmp;
        }

    } /* if projectCoordinate >= 0 */


  /*  attach the metricParams struct as 'meta-info' to the output */
  metric->meta = (*metricParams);

  return metric;

} /* XLALDopplerFstatMetric() */

/** Function to the compute the FmetricAtoms_t, from which the F-metric and Fisher-matrix can be computed.
 */
FmetricAtoms_t*
XLALComputeAtomsForFmetric ( const DopplerMetricParams *metricParams,  	/**< input parameters determining the metric calculation */
			     const EphemerisData *edat			/**< ephemeris data */
			     )
{
  FmetricAtoms_t *ret;		/* return struct */
  intparams_t intparams = empty_intparams;

  UINT4 dim, numDet, i=-1, j=-1, X;		/* index counters */
  REAL8 A, B, C;			/* antenna-pattern coefficients (gsl-integrated) */

  const LIGOTimeGPS *refTime, *startTime;
  const DopplerCoordinateSystem *coordSys;

  REAL8 max_relerr = 0;
  REAL8 relerr_thresh = 1e-2;	/* relatively tolerant integration-threshold */

  /* ---------- sanity/consistency checks ---------- */
  if ( !metricParams || !edat ) {
    XLALPrintError ("\n%s: Illegal NULL pointer passed!\n\n", __func__);
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }

  startTime = &(metricParams->startTime);
  refTime   = &(metricParams->signalParams.Doppler.refTime);

  dim = metricParams->coordSys.dim;	/* shorthand: number of Doppler dimensions */
  numDet = metricParams->detInfo.length;
  coordSys = &(metricParams->coordSys);

  /* ----- create output structure ---------- */
  if ( (ret = XLALCreateFmetricAtoms ( dim )) == NULL ) {
    XLALPrintError ("%s: call to XLALCreateFmetricAtoms (%s) failed. errno = %d\n\n", __func__, dim, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* ---------- set up integration parameters ---------- */
  intparams.detMotionType = metricParams->detMotionType;
  intparams.dopplerPoint = &metricParams->signalParams.Doppler;
  intparams.startTime = XLALGPSGetREAL8 ( startTime );
  intparams.refTime   = XLALGPSGetREAL8 ( refTime );
  intparams.Tspan = metricParams->Tspan;
  intparams.edat = edat;

  /* if using 'global correlation' frequency variables, determine the highest spindown order: */
  UINT4 maxorder = findHighestGCSpinOrder ( coordSys );

  /* compute rOrb(t) derivatives at reference time */
  if ( (intparams.rOrb_n = XLALComputeOrbitalDerivatives ( maxorder, &intparams.dopplerPoint->refTime, edat )) == NULL ) {
    XLALPrintError ("%s: XLALComputeOrbitalDerivatives() failed.\n", __func__);
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }

  /* ----- integrate antenna-pattern coefficients A, B, C */
  A = B = C = 0;
  for ( X = 0; X < numDet; X ++ )
    {
      REAL8 weight = metricParams->detInfo.detWeights[X];
      REAL8 av, relerr;
      intparams.site = &(metricParams->detInfo.sites[X]);

      intparams.deriv1 = DOPPLERCOORD_NONE;
      intparams.deriv2 = DOPPLERCOORD_NONE;

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
	      REAL8 weight = metricParams->detInfo.detWeights[X];
	      REAL8 av, relerr;
	      intparams.site = &(metricParams->detInfo.sites[X]);

	      /* ------------------------------ */
	      intparams.deriv1 = coordSys->coordIDs[i];
	      intparams.deriv2 = coordSys->coordIDs[j];

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
	      intparams.deriv1 = coordSys->coordIDs[i];
	      intparams.deriv2 = DOPPLERCOORD_NONE;

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
	      intparams.deriv1 = DOPPLERCOORD_NONE;
	      intparams.deriv2 = coordSys->coordIDs[j];

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

	  gsl_vector_set (ret->a_a_i, i, a_a_i);
	  gsl_vector_set (ret->a_b_i, i, a_b_i);
	  gsl_vector_set (ret->b_b_i, i, b_b_i);


	  gsl_matrix_set (ret->a_a_i_j, i, j, a_a_i_j);
	  gsl_matrix_set (ret->a_a_i_j, j, i, a_a_i_j);

	  gsl_matrix_set (ret->a_b_i_j, i, j, a_b_i_j);
	  gsl_matrix_set (ret->a_b_i_j, j, i, a_b_i_j);

	  gsl_matrix_set (ret->b_b_i_j, i, j, b_b_i_j);
	  gsl_matrix_set (ret->b_b_i_j, j, i, b_b_i_j);

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
  XLALPrintError ( "%s: XLALAverage_am1_am2_Phi_i_Phi_j() FAILED with errno = %d: am1 = %d, am2 = %d, i = %d : '%s', j = %d : '%s'\n",
		   __func__, xlalErrno, intparams.amcomp1, intparams.amcomp2,
                   i, XLALDopplerCoordinateName(intparams.deriv1),
                   j, XLALDopplerCoordinateName(intparams.deriv2) );
  XLAL_ERROR_NULL( XLAL_EFUNC );

} /* XLALComputeAtomsForFmetric() */


/** Free a DopplerMetric structure */
void
XLALDestroyDopplerMetric ( DopplerMetric *metric )
{
  if ( !metric )
    return;

  if ( metric->g_ij ) 	gsl_matrix_free ( metric->g_ij );
  if ( metric->gFav_ij) gsl_matrix_free ( metric->gFav_ij );
  if ( metric->m1_ij ) 	gsl_matrix_free ( metric->m1_ij );
  if ( metric->m2_ij ) 	gsl_matrix_free ( metric->m2_ij );
  if ( metric->m3_ij ) 	gsl_matrix_free ( metric->m3_ij );
  if ( metric->Fisher_ab ) gsl_matrix_free ( metric->Fisher_ab );

  XLALFree ( metric );

  return;

} /* XLALDestroyDopplerMetric() */


/** Parse a detector-motion type string into the corresponding enum-number,
 */
int
XLALParseDetectorMotionString ( const CHAR *detMotionString )
{
  int i;

  if ( ! detMotionString ) {
    XLAL_ERROR ( XLAL_EINVAL );
  }

  for ( i=0; i < DETMOTION_LAST; i ++ )
    {
      if ( DetectorMotionNames[i] && strcmp ( detMotionString, DetectorMotionNames[i] ) )
	continue;
      return i;	/* found the right entry */
    }

  XLALPrintError ("\nCould not parse '%s' into a valid detector-motion type!\n\n", detMotionString );
  XLAL_ERROR ( XLAL_EINVAL );

} /* XLALParseDetectorMotionString() */


/** Provide a pointer to a static string containing the DopplerCoordinate-name
 * cooresponding to the enum DopplerCoordinateID
 */
const CHAR *
XLALDetectorMotionName ( DetectorMotionType detType )
{
  if ( detType >= DETMOTION_LAST ) {
    XLAL_ERROR_NULL ( XLAL_EINVAL, "detector-motion type '%d' outside valid range [0, %d]\n\n", detType, DETMOTION_LAST - 1 );
  }

  if ( !DetectorMotionNames[detType] ) {
    XLAL_ERROR_NULL ( XLAL_EINVAL, "detector-motion type '%d' has no associated name\n\n", detType );
  }

  return ( DetectorMotionNames[detType] );

} /* XLALDetectorMotionName() */



/** Parse a DopplerCoordinate-name into the corresponding DopplerCoordinateID
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

/** Given a LALStringVector of coordinate-names, parse them into a
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



/** Provide a pointer to a static string containing the DopplerCoordinate-name
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


/** Provide a pointer to a static string containing the a descriptive
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

/** Return a string (allocated here) containing a full name - helpstring listing
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


/** Parse string-vectors (typically input by user) of detector-names
 * and relative noise-weights, and return a MultiDetectorInfo struct.
 *
 * NOTE: you can pass detWeights == NULL, corresponding to equal-sensitivity detectors,
 * ie. all noise-weights equal.
 *
 * NOTE: the input noise-weights dont have to be normalized, but the
 * returned noise-weights will be properly normalized, i.e. \f$\sum_{i=1}^N w_i = 1\f$.
 *
 * Return  0 == OK, nonzero == ERROR
 */
int
XLALParseMultiDetectorInfo ( MultiDetectorInfo *detInfo,	/**< [out] parsed detector-info struct */
			     const LALStringVector *detNames,	/**< [in] list of detector names */
			     const LALStringVector *detWeights	/**< [in] list of (strings) with detector weights (NULL if all 1) */
			     )
{
  UINT4 X, numDet;
  REAL8 totalWeight;

  if ( !detInfo || !detNames || (detNames->length == 0) ) {
    XLALPrintError ("\n%s: Illegal NULL pointer input\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  numDet = detNames->length;
  if ( detWeights && (detWeights->length != numDet ) ) {
    XLALPrintError ("\n%s: Illegal input: number of noise-weights must agree with number of detectors\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* initialize empty return struct */
  memset ( detInfo, 0, sizeof(*detInfo) );

  detInfo->length = numDet;

  totalWeight = 0;
  /* parse input strings and fill detInfo */
  for ( X = 0; X < numDet; X ++ )
    {
      LALDetector *ifo;
      /* first parse detector name */
      if ( ( ifo = XLALGetSiteInfo ( detNames->data[X] ) ) == NULL ) {
	XLALPrintError ("%s: Failed to get site-info for detector '%s'\n", __func__, detNames->data[X] );
	XLAL_ERROR ( XLAL_EINVAL );
      }
      detInfo->sites[X] = (*ifo);
      XLALFree ( ifo );

      /* parse noise weights if any */
      if ( detWeights )
	{
	  if ( 1 != sscanf ( detWeights->data[X], "%lf", &(detInfo->detWeights[X]) ) )
	    {
	      XLALPrintError ("%s: Failed to parse noise-weight '%s' into float.\n", __func__, detWeights->data[X] );
	      XLAL_ERROR ( XLAL_EINVAL );
	    }
	} /* if detWeights */
      else
	detInfo->detWeights[X] = 1;

      totalWeight += detInfo->detWeights[X];

    } /* for X < numDet */

  /* normalized noise-weights to sum weights = 1 */
  for ( X = 0; X < numDet; X ++ )
    detInfo->detWeights[X] /= totalWeight;

  return XLAL_SUCCESS;

} /* XLALParseMultiDetectorInfo() */


/** Free a FmetricAtoms_t structure, allowing any pointers to be NULL
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



/** Allocate an FmetricAtoms_t structure for given number of dimension.
 */
FmetricAtoms_t*
XLALCreateFmetricAtoms ( UINT4 dim )
{
  FmetricAtoms_t *ret;		/* output structure */

  if ( ( ret = XLALCalloc(1,sizeof(*ret))) == NULL ) {
    XLALPrintError ( "%s: XLALCalloc(1,%s) failed.\n", __func__, sizeof(*ret));
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( ( ret->a_a_i = gsl_vector_calloc (dim)) == NULL ) {
    XLALPrintError ( "%s: a_a_i = gsl_vector_calloc (%d) failed.\n", __func__, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( ( ret->a_b_i = gsl_vector_calloc (dim)) == NULL ) {
    XLALPrintError ( "%s: a_b_i = gsl_vector_calloc (%d) failed.\n", __func__, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( ( ret->b_b_i = gsl_vector_calloc (dim)) == NULL ) {
    XLALPrintError ( "%s: b_b_i = gsl_vector_calloc (%d) failed.\n", __func__, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }


  if ( ( ret->a_a_i_j = gsl_matrix_calloc (dim, dim)) == NULL ) {
    XLALPrintError ( "%s: a_a_i_j = gsl_matrix_calloc (%d,%d) failed.\n", __func__, dim, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( ( ret->a_b_i_j = gsl_matrix_calloc (dim, dim)) == NULL ) {
    XLALPrintError ( "%s: a_b_i_j = gsl_matrix_calloc (%d,%d) failed.\n", __func__, dim, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( ( ret->b_b_i_j = gsl_matrix_calloc (dim, dim)) == NULL ) {
    XLALPrintError ( "%s: b_b_i_j = gsl_matrix_calloc (%d,%d) failed.\n", __func__, dim, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  return ret;

} /* XLALCreateFmetricAtoms() */


/** Compute the 'F-metric' gF_ij (and also gFav_ij, m1_ij, m2_ij, m3_ij)
 * from the given FmetricAtoms and the signal amplitude parameters.
 *
 */
DopplerMetric*
XLALComputeFmetricFromAtoms ( const FmetricAtoms_t *atoms, REAL8 cosi, REAL8 psi )
{
  DopplerMetric *metric;		/* output matrix */

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
    XLALPrintError ("%s: XLALCalloc ( 1, %d) failed.\n\n", sizeof(*metric) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  metric->gF_ij = gsl_matrix_calloc ( dim, dim );
  metric->gFav_ij = gsl_matrix_calloc ( dim, dim );
  metric->m1_ij = gsl_matrix_calloc ( dim, dim );
  metric->m2_ij = gsl_matrix_calloc ( dim, dim );
  metric->m3_ij = gsl_matrix_calloc ( dim, dim );

  if ( !metric->gF_ij || !metric->gFav_ij || !metric->m1_ij || !metric->m2_ij || !metric->m3_ij ) {
    XLALPrintError ("%s: failed to gsl_matrix_calloc(%d,%d) for gF_ij, gFav_ij, m1_ij, m2_ij, m3_ij\n\n", __func__, dim, dim );
    XLALDestroyDopplerMetric ( metric );
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

  metric->maxrelerr_gF = atoms->maxrelerr;

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

	  /* trivial assignments, see Eq.(76) in \ref Prix07 */
	  P1_ij = a_a_i_j;
	  P2_ij = b_b_i_j;
	  P3_ij = a_b_i_j;

	  /* bit more involved, see Eq.(80)-(82) in \ref Prix07 [includes *explicit* index-symmetrization!!] */
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


/** Function to compute *full* 4+n dimensional Fisher matric for the
 *  full CW parameter-space of Amplitude + Doppler parameters !
 */
gsl_matrix*
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


/** Calculate the projected metric onto the subspace orthogonal to coordinate-axis 'c', namely
 * ret_ij = g_ij - ( g_ic * g_jc / g_cc ) , where c is the value of the projected coordinate
 * The output-matrix is allocate here
 *
 * return 0 = OK, -1 on error.
 */
gsl_matrix *
XLALProjectMetric ( const gsl_matrix * g_ij, const UINT4 c )
{
  UINT4 i,j, dim1, dim2;
  gsl_matrix *ret_ij;

  if ( !g_ij ) {
    XLALPrintError ("%s: invalid NULL input 'g_ij'.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  dim1 = g_ij->size1;
  dim2 = g_ij->size2;

  if ( dim1 != dim2 ) {
    XLALPrintError ( "%s: input matrix g_ij must be square! (got %d x %d)\n", __func__, dim1, dim2 );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  if ( (ret_ij = gsl_matrix_alloc ( dim1, dim2 )) == NULL ) {
    XLALPrintError ("%s: failed to gsl_matrix_alloc(%d, %d)\n", __func__, dim1, dim2 );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  for ( i=0; i < dim1; i++)
    {
    for ( j=0; j < dim2; j++ )
      {
        if ( i==c || j==c )
          {
            gsl_matrix_set ( ret_ij, i, j, 0.0 );
          }
        else
          {
            double proj = gsl_matrix_get(g_ij, i, j) - (gsl_matrix_get(g_ij, i, c) * gsl_matrix_get(g_ij, j, c) / gsl_matrix_get(g_ij, c, c));
            gsl_matrix_set ( ret_ij, i, j, proj );
          }
      } /* for j < dim2 */

    } /* for i < dim1 */

  return ret_ij;

} /* XLALProjectMetric() */

/** Convert 3-D vector from equatorial into ecliptic coordinates */
void
equatorialVect2ecliptic ( vect3D_t out, const vect3D_t in )
{
  static mat33_t rotEqu2Ecl = { { 1.0,        0,       0 },
                                { 0.0,  LAL_COSIEARTH, LAL_SINIEARTH },
                                { 0.0, -LAL_SINIEARTH, LAL_COSIEARTH } };

  matrix33_in_vect3 ( out, rotEqu2Ecl, in );

} /* equatorialVect2ecliptic() */

/** Convert 3-D vector from ecliptic into equatorial coordinates */
void
eclipticVect2equatorial ( vect3D_t out, const vect3D_t in )
{
  static mat33_t rotEcl2Equ =  { { 1.0,        0,       0 },
                                 { 0.0,  LAL_COSIEARTH, -LAL_SINIEARTH },
                                 { 0.0,  LAL_SINIEARTH,  LAL_COSIEARTH } };

  matrix33_in_vect3 ( out, rotEcl2Equ, in );

} /* eclipticVect2equatorial() */

/** compute matrix product mat . vect */
void
matrix33_in_vect3 ( vect3D_t out, mat33_t mat, const vect3D_t in )
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

} /* matrix33_in_vect3() */

/** Compute time-derivatives up to 'maxorder' of the Earths' orbital position vector
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
    XLALPrintError ("%s: failed to XLALCalloc(1,%d)\n", __func__, sizeof(*ret) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  ret->length = maxorder + 1;
  if ( (ret->data = XLALCalloc ( ret->length, sizeof(*ret->data) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc(%d,%d)\n", __func__, ret->length, sizeof(*ret->data) );
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

/** Return the highest 'global-correlation' spindown order found in this coordinate system.
    Counting nu0 = order1, nu1 = order2, nu2 = order3, ...,
    order = 0 therefore means there are no GC spin coordinates at all
*/
UINT4
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


/**
 * Return a metric in "naturalized" coordinates.
 * Frequency coordinates of spindown order \f$s\f$ are scaled by
 * \f[ \frac{2\pi}{(s+1)!} \left(\frac{T}{2}\right)^{s+1} \f]
 * where \f$T\f$ is the observation time-span.
 * Sky coordinates are scaled by
 * \f[ \frac{2\pi \bar{f} R_{ES}}{c} \f]
 * where \f$\bar{f}\f$ is a fiducial frequency and
 * \f$R_{ES}\f$ the mean Earth--Sun distance.
 *
 * Returns NULL on error, otherwise a new matrix is allocated.
 */
gsl_matrix* XLALNaturalizeMetric(
  const gsl_matrix* g_ij,			/**< [in] Input metric */
  const DopplerMetricParams *metricParams	/**< [in] Input parameters used to calculate g_ij */
  )
{

  /* Check input */
  XLAL_CHECK_NULL( g_ij, XLAL_EINVAL );
  XLAL_CHECK_NULL( g_ij->size1 == g_ij->size2, XLAL_EINVAL, "Input matrix g_ij must be square! (got %d x %d)\n", g_ij->size1, g_ij->size2 );

  /* Compute naturalization scale */
  double nat_scale[g_ij->size1];
  for (size_t i = 0; i < g_ij->size1; ++i) {
    const DopplerCoordinateID coordID = metricParams->coordSys.coordIDs[i];
    const double Freq = metricParams->signalParams.Doppler.fkdot[0];
    const double T = metricParams->Tspan;
    double scale;
    switch (coordID) {
    case DOPPLERCOORD_NONE:
      scale = 1;
      break;

    case DOPPLERCOORD_FREQ:
    case DOPPLERCOORD_GC_NU0:
      scale = LAL_TWOPI * LAL_FACT_INV[1] * (0.5 * T);
      break;

    case DOPPLERCOORD_F1DOT:
    case DOPPLERCOORD_GC_NU1:
      scale = LAL_TWOPI * LAL_FACT_INV[2] * POW2(0.5 * T);
      break;

    case DOPPLERCOORD_F2DOT:
    case DOPPLERCOORD_GC_NU2:
      scale = LAL_TWOPI * LAL_FACT_INV[3] * POW3(0.5 * T);
      break;

    case DOPPLERCOORD_F3DOT:
    case DOPPLERCOORD_GC_NU3:
      scale = LAL_TWOPI * LAL_FACT_INV[4] * POW4(0.5 * T);
      break;

    case DOPPLERCOORD_ALPHA:
    case DOPPLERCOORD_DELTA:
    case DOPPLERCOORD_N2X_EQU:
    case DOPPLERCOORD_N2Y_EQU:
    case DOPPLERCOORD_N2X_ECL:
    case DOPPLERCOORD_N2Y_ECL:
    case DOPPLERCOORD_N3X_EQU:
    case DOPPLERCOORD_N3Y_EQU:
    case DOPPLERCOORD_N3Z_EQU:
    case DOPPLERCOORD_N3X_ECL:
    case DOPPLERCOORD_N3Y_ECL:
    case DOPPLERCOORD_N3Z_ECL:
    case DOPPLERCOORD_N3SX_EQU:
    case DOPPLERCOORD_N3SY_EQU:
    case DOPPLERCOORD_N3OX_ECL:
    case DOPPLERCOORD_N3OY_ECL:
      scale = LAL_TWOPI * Freq * rOrb_c;
      break;

    default:
      XLAL_ERROR_NULL(XLAL_EINVAL, "Unknown phase-derivative type '%d'\n", coordID );
    }
    nat_scale[i] = scale;
  }

  /* Allocate return matrix */
  gsl_matrix* ret_ij = gsl_matrix_alloc ( g_ij->size1, g_ij->size2 );
  XLAL_CHECK_NULL( ret_ij, XLAL_ENOMEM );

  /* Rescale metric to naturalized coordinates */
  for (size_t i = 0; i < g_ij->size1; ++i) {
    for (size_t j = 0; j < g_ij->size2; ++j) {
      double tmp = gsl_matrix_get(g_ij, i, j);
      tmp /= nat_scale[i] * nat_scale[j];
      gsl_matrix_set(ret_ij, i, j, tmp);
    }
  }

  return ret_ij;

} /* XLALNaturalizeMetric() */


/** "DiagNormalize" a metric matrix.
 * DiagNormalization means normalize metric by its diagonal, namely apply the transformation
 * G_ij = g_ij /sqrt(g_ii * g_jj), to all elements, resulting in lower
 * condition number and unit diagonal elements.
 *
 * return NULL on error, otherwise new matrix is allocated here.
 */
gsl_matrix *
XLALDiagNormalizeMetric ( const gsl_matrix * g_ij )
{
  const char *fn = __func__;
  UINT4 i,j, dim1, dim2;
  gsl_matrix *ret_ij;

  if ( !g_ij ) {
    XLALPrintError ("%s: invalid NULL input 'g_ij'.\n", fn );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  dim1 = g_ij->size1;
  dim2 = g_ij->size2;

  if ( dim1 != dim2 ) {
    XLALPrintError ( "%s: input matrix g_ij must be square! (got %d x %d)\n", fn, dim1, dim2 );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  if ( (ret_ij = gsl_matrix_alloc ( dim1, dim2 )) == NULL ) {
    XLALPrintError ("%s: failed to gsl_matrix_alloc(%d, %d)\n", fn, dim1, dim2 );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  for ( i=0; i < dim1; i++)
    {
    for ( j=0; j < dim2; j++ )
      {
        if ( i == j )
          {
            gsl_matrix_set ( ret_ij, i, j, 1.0 );	/* use exact result on diagonal */
          }
        else
          {
            double gtmp_ii = gsl_matrix_get(g_ij, i, i);
            double gtmp_jj = gsl_matrix_get(g_ij, j, j);

            if ( (gtmp_ii <= 0) || (gtmp_jj <= 0 ) ) {
              XLALPrintError ("%f: DiagNormalize not defined for non-positive diagonal elements! i=%d, j=%d, g_ii=%g, g_jj=%g\n", i, j, gtmp_ii, gtmp_jj );
              XLAL_ERROR_NULL ( XLAL_EDOM );
            }

            double gtmp_ij = gsl_matrix_get(g_ij, i, j);
            double new_ij = gtmp_ij / sqrt ( gtmp_ii * gtmp_jj );

            gsl_matrix_set ( ret_ij, i, j, new_ij );
          }
      } /* for j < dim2 */

    } /* for i < dim1 */

  return ret_ij;

} /* XLALDiagNormalizeMetric() */
