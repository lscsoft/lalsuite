/*
 * Copyright (C) 2008 Reinhard Prix
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
 * \author Reinhard Prix
 * \ingroup PulsarMetric
 * \brief Function to compute the full F-statistic metric, including
 *  antenna-pattern functions from multi-detector, as derived in \ref Prix07.
 *
 * Revision: $Id$
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

#define IN_UNIVERSALDOPPLERMETRICC
#include <lal/UniversalDopplerMetric.h>

/*---------- DEFINES ----------*/
#define TRUE  (1==1)
#define FALSE (0==1)

/** copy 3 components of Euklidean vector */
#define COPY_VECT(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])
#define MULT_VECT(lam, v) do{ (v)[0] *= (lam); (v)[1] *= (lam); (v)[2] *= (lam); } while(0)
#define ADD_VECT(dst,src) do { (dst)[0] += (src)[0]; (dst)[1] += (src)[1]; (dst)[2] += (src)[2]; } while(0)
#define SUB_VECT(dst,src) do { (dst)[0] -= (src)[0]; (dst)[1] -= (src)[1]; (dst)[2] -= (src)[2]; } while(0)

#define SQUARE(x) ((x) * (x))

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYMAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MYMIN(a,b) ( (a) < (b) ? (a) : (b) )

/** 5-point derivative formulas */
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
} intparams_t;


/*---------- empty initializers ---------- */
static LALStatus empty_status;
static EmissionTime empty_EmissionTime;
static intparams_t empty_intparams;
static PosVel3D_t empty_PosVel3D_t;
static const PulsarTimesParamStruc empty_PulsarTimesParamStruc;

DopplerMetricParams empty_DopplerMetricParams;
DopplerCoordinateSystem empty_DopplerCoordinateSystem;
MultiDetectorInfo empty_MultiDetectorInfo;

/*---------- Global variables ----------*/
NRCSID( UNIVERSALDOPPLERMETRICC, "$Id$");

BOOLEAN outputIntegrand = 0;

/* Some local constants. */
#define rOrb_c  (LAL_AU_SI / LAL_C_SI)
#define vOrb_c  (LAL_TWOPI * LAL_AU_SI / LAL_C_SI / LAL_YRSID_SI)
// sin,cos(LAL_IEARTH);
#define cosiEcl 0.917482062157619
#define siniEcl 0.397777155727931

/*---------- internal prototypes ----------*/
DopplerMetric* XLALComputeFmetricFromAtoms ( const FmetricAtoms_t *atoms, REAL8 cosi, REAL8 psi );
gsl_matrix* XLALComputeFisherFromAtoms ( const FmetricAtoms_t *atoms, const PulsarAmplitudeParams *Amp );

double CW_am1_am2_Phi_i_Phi_j ( double tt, void *params );
double CWPhaseDeriv_i ( double tt, void *params );

double XLALAverage_am1_am2_Phi_i_Phi_j ( const intparams_t *params, double *relerr_max );
double CWPhase_cov_Phi_ij ( const intparams_t *params, double *relerr_max );

int XLALPtolemaicPosVel ( PosVel3D_t *posvel, const LIGOTimeGPS *tGPS );
gsl_matrix *XLALProjectMetric ( const gsl_matrix * g_ij, const UINT4 c );

int equatorialVect2ecliptic ( vect3D_t *out, vect3D_t * const in );
int eclipticVect2equatorial ( vect3D_t *out, vect3D_t * const in );
int matrix33_in_vect3 ( vect3D_t *out, mat33_t * mat, vect3D_t * const in );

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
  const CHAR *fn = "XLALAverage_am1_am2_Phi_i_Phi_j()";

  intparams_t par = (*params);	/* struct-copy, as the 'deriv' field has to be changeable */
  gsl_function integrand;
  double epsrel = 1e-4;
  /* NOTE: this level of accuracy should be compatible with AM-coefficients involved
   * which are computed in REAL4 precision. We therefor cannot go lower than this it seems,
   * otherwise the gsl-integration fails to converge in some cases.
   */
  double epsabs = 0;
  size_t neval;
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

  integrand.function = &CW_am1_am2_Phi_i_Phi_j;
  for (n=0; n < Nseg; n ++ )
    {
      REAL8 ti = 1.0 * n * dT;
      REAL8 tf = MYMIN( (n+1.0) * dT, 1.0 );
      REAL8 err_n, res_n;

      XLAL_CALLGSL ( stat = gsl_integration_qng (&integrand, ti, tf, epsabs, epsrel, &res_n, &err_n, &neval) );
      /* NOTE: we don't fail if the requested level of accuracy was not reached, rather we only output a warning
       * and try to estimate the final accuracy of the integral
       */
      if ( stat != 0 )
        {
          XLALPrintWarning ( "\n%s: GSL-integration 'gsl_integration_qng()' of <am1_am2_Phi_i Phi_j> did not reach requested precision!\n", fn );
          XLALPrintWarning ("Segment n=%d, neval=%d: Result = %g, abserr=%g ==> relerr = %.2e > %.2e\n", n, neval, res_n, err_n, fabs(err_n/res_n), epsrel);
          /* XLAL_ERROR_REAL8( fn, XLAL_EFUNC ); */
        }

      res += res_n;
      abserr2 += SQUARE ( err_n );

    } /* for i < Nseg */

  REAL8 relerr = sqrt(abserr2) / fabs(res);
  if ( relerr_max )
    (*relerr_max) = relerr;

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
  const CHAR *fn = "CW_am1_am2_Phi_i_Phi_j()";
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

      if ( XLALComputeAntennaPatternCoeffs ( &ai, &bi, &skypos, &ttGPS, par->site, par->edat ) ) {
	XLALPrintError ( "%s: Call to XLALComputeAntennaPatternCoeffs() failed!\n", fn);
	XLAL_ERROR( fn, XLAL_EFUNC );
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
  const CHAR *fn = "CWPhaseDeriv_i()";
  REAL8 ret;
  intparams_t *par = (intparams_t*) params;
  vect3D_t nn_equ, nn_ecl;	/* skypos unit vector */
  vect3D_t nDeriv_i;	/* derivative of sky-pos vector wrt i */
  PosVel3D_t posvel = empty_PosVel3D_t;
  vect3D_t detpos_ecl, detpos_equ;

  REAL8 Freq = par->dopplerPoint->fkdot[0];
  REAL8 Tspan = par->Tspan;
  static REAL8 kfactinv[] = { 1.0, 1.0/1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0 };	/* 1/k! */

  /* get skypos-vector */
  REAL8 cosa = cos(par->dopplerPoint->Alpha);
  REAL8 sina = sin(par->dopplerPoint->Alpha);
  REAL8 cosd = cos(par->dopplerPoint->Delta);
  REAL8 sind = sin(par->dopplerPoint->Delta);
  /* ... in an equatorial coordinate-frame */
  nn_equ[0] = cosd * cosa;
  nn_equ[1] = cosd * sina;
  nn_equ[2] = sind;

  /* and in an ecliptic coordinate-frame */
  equatorialVect2ecliptic ( &nn_ecl, (vect3D_t * const )&nn_equ );

  if ( abs(nn_ecl[2]) < 1e-6 )	/* avoid singularity at ecliptic equator */
    nn_ecl[2] = 1e-6;

  /*
  vect3D_t nn_ecl0, nn_equ0;
  printf ("equatorialVect2ecliptic(): nn_ecl0 = { %g, %g, %g}\n", nn_ecl0[0], nn_ecl0[1], nn_ecl0[2] );
  nn_ecl[0] = nn_equ[0];
  nn_ecl[1] = cosiEcl * nn_equ[1] + siniEcl * nn_equ[2];
  nn_ecl[2] = cosiEcl * nn_equ[2] - siniEcl * nn_equ[1];
  printf ("manual conversion        : nn_ecl  = { %g, %g, %g}\n", nn_ecl[0], nn_ecl[1], nn_ecl[2] );
  eclipticVect2equatorial ( nn_equ0, nn_ecl0 );
  printf ("eclipticVect2equatorial(): nn_equ0 = { %g, %g, %g}\n", nn_equ0[0], nn_equ0[1], nn_equ0[2] );
  printf ("original input           : nn_equ  = { %g, %g, %g}\n", nn_equ[0], nn_equ[1], nn_equ[2] );
  exit(0);
  */

  /* get current detector position r(t) */
  REAL8 ttSI = par->startTime + tt * Tspan;	/* current GPS time in seconds */
  LIGOTimeGPS ttGPS;
  XLALGPSSetREAL8( &ttGPS, ttSI );

  if ( XLALDetectorPosVel ( &posvel, &ttGPS, par->site, par->edat, par->detMotionType ) ) {
    XLALPrintError ( "%s: Call to XLALDetectorPosVel() failed!\n", fn);
    XLAL_ERROR( fn, XLAL_EFUNC );
  }
  COPY_VECT ( detpos_equ, posvel.pos );


  /* get 'reduced' detector position of order 'n': r_n(t),
   * defined as: r_n(t) = r(t) - dot{r_orb}(tau_ref) - 1/2! ddot{r_orb}(tau_re) - ....
   */
  vect3D_t rr_ord;
  COPY_VECT ( rr_ord, posvel.pos );


  /* convert detector position in ecliptic coordinates */
  equatorialVect2ecliptic ( &detpos_ecl, &detpos_equ );


  /* account for referenceTime != startTime */
  REAL8 tau0 = ( par->startTime - par->refTime ) / Tspan;

  /* correct for time-delay from SSB to detector, neglecting relativistic effects */
  REAL8 dTRoemerSI = SCALAR(nn_equ, detpos_equ );
  REAL8 dTRoemer = dTRoemerSI / Tspan;		/* SSB time-delay in 'natural units' */

  /* barycentric time-delay since reference time, measured in units of Tspan */
  REAL8 tau = tau0 + tt + dTRoemer;

  REAL8 nNat = Freq * Tspan * 1e-4;	/* 'natural sky-units': Freq * Tspan * V/c */

  switch ( par->deriv )
    {
      /* ----- sky derivatives ----- */
    case DOPPLERCOORD_ALPHA_RAD:				/* longitude/right-ascension/Alpha in radians */
    case DOPPLERCOORD_ALPHA_NAT:				/* longitude/right-ascension/Alpha in natural units */
      nDeriv_i[0] = - cosd * sina;
      nDeriv_i[1] =   cosd * cosa;
      nDeriv_i[2] =   0;

      ret = LAL_TWOPI * Freq * SCALAR(posvel.pos, nDeriv_i);	/* dPhi/dAlpha = 2 pi f (r/c) . (dn/dAlpha) */
      if ( par->deriv == DOPPLERCOORD_ALPHA_NAT )
        ret *= nNat;
      break;

    case DOPPLERCOORD_DELTA_RAD:				/* latitude/declination/Delta in radians */
    case DOPPLERCOORD_DELTA_NAT:				/* latitude/declination/Delta in natural units */
      nDeriv_i[0] = - sind * cosa;
      nDeriv_i[1] = - sind * sina;
      nDeriv_i[2] =   cosd;

      ret = LAL_TWOPI * Freq * SCALAR(posvel.pos, nDeriv_i);	/* dPhi/dDelta = 2 pi f (r/c) . (dn/dDelta) */
      if ( par->deriv == DOPPLERCOORD_DELTA_NAT )
        ret *= nNat;
      break;

    case DOPPLERCOORD_NECL_X_NAT:
      ret = ( detpos_ecl[0] - (nn_ecl[0]/nn_ecl[2]) * detpos_ecl[2] ) / rOrb_c;
      break;
    case DOPPLERCOORD_NECL_Y_NAT:
      ret = ( detpos_ecl[1] - (nn_ecl[1]/nn_ecl[2]) * detpos_ecl[2] ) / rOrb_c;
      break;

      /* experimental 'global correlation' sky coordinate, if holding {nu, nu1, ...} fixed.
       * Note: derivatives wrt to nu, nu1, ... are equivalent to f, fdot, ...
       */
    case DOPPLERCOORD_NEQU_X_GC:
      break;

    case DOPPLERCOORD_NEQU_Y_GC:
      break;

      /* experimental: unconstrained skypos vector n3 */
    case DOPPLERCOORD_N3X:
      ret = detpos_ecl[0] / rOrb_c;
      break;
    case DOPPLERCOORD_N3Y:
      ret = detpos_ecl[1] / rOrb_c;
      break;
    case DOPPLERCOORD_N3Z:
      ret = detpos_ecl[2] / rOrb_c;
      break;

    case DOPPLERCOORD_NEQU_X_NAT:
      ret = ( detpos_equ[0] - (nn_equ[0]/nn_equ[2]) * detpos_equ[2] ) / rOrb_c;
      break;
    case DOPPLERCOORD_NEQU_Y_NAT:
      ret = ( detpos_equ[1] - (nn_equ[1]/nn_equ[2]) * detpos_equ[2] ) / rOrb_c;
      break;



      /* ----- frequency derivatives SI-units ----- */
    case DOPPLERCOORD_FREQ_SI:
    case DOPPLERCOORD_FREQ_NAT:				/* om0 = 2pi f T */
      ret = tau;					/* in natural units: dPhi/dom0 = tau */
      if ( par->deriv == DOPPLERCOORD_FREQ_SI )
        ret *= LAL_TWOPI * Tspan * kfactinv[1];		/* dPhi/dFreq = 2 * pi * tSSB_i */
      break;

    case DOPPLERCOORD_F1DOT_SI:
    case DOPPLERCOORD_F1DOT_NAT:			/* om1 = 2pi f/2! T^2 */
      ret = tau * tau;					/* in natural units: dPhi/dom1 = tau^2 */
      if ( par->deriv == DOPPLERCOORD_F1DOT_SI )
        ret *= LAL_TWOPI * (Tspan*Tspan) * kfactinv[2];/* dPhi/df1dot = 2pi * (tSSB_i)^2/2! */
      break;

    case DOPPLERCOORD_F2DOT_SI:
    case DOPPLERCOORD_F2DOT_NAT:			/* om2 = 2pi f/3! T^3 */
      ret =  tau * tau * tau;				/* in natural units: dPhi/dom2 = tau^3 */
      if ( par->deriv == DOPPLERCOORD_F2DOT_SI )
        ret *= LAL_TWOPI * (Tspan*Tspan*Tspan) * kfactinv[3];/* dPhi/f2dot = 2pi * (tSSB_i)^3/3! */
      break;

    case DOPPLERCOORD_F3DOT_SI:
    case DOPPLERCOORD_F3DOT_NAT:			/* om3 = 2pi f/4! T^4 */
      ret = tau * tau * tau * tau;			/* in natural units: dPhi/dom3 = tau^4 */
      if ( par->deriv == DOPPLERCOORD_F3DOT_SI )
        ret *= LAL_TWOPI * (Tspan*Tspan*Tspan*Tspan) * kfactinv[4];/* dPhi/df3dot = 2pi * (tSSB_i)^4/4! */
      break;

    default:
      XLALPrintError("%s: Unknown phase-derivative type '%d'\n", fn, par->deriv );
      XLAL_ERROR( fn, XLAL_EINVAL );
      break;
    } /* switch phderiv */

  return ret;

} /* CWPhaseDeriv_i() */


/** Given a GPS time and detector, return the current position (and velocity) of the detector.
 *
 * NOTE: the 'special' flag allows to simulate artifically truncated
 * detector motions, such as pure orbital motion.
 *
 */
int
XLALDetectorPosVel ( PosVel3D_t *posvel,	/**< [out] instantaneous position and velocity vector */
		     const LIGOTimeGPS *tGPS,	/**< [in] GPS time */
		     const LALDetector *site,	/**< [in] detector info */
		     const EphemerisData *edat,	/**< [in] ephemeris data */
		     DetectorMotionType special	/**< [in] 'special' flag: 0 = full motion, 1 = pure orbital, 2 = pure spin */
		     )
{
  const CHAR *fn = "XLALDetectorPosition()";
  EarthState earth;
  BarycenterInput baryinput = empty_BarycenterInput;
  LALStatus status = empty_status;
  EmissionTime emit = empty_EmissionTime;
  PosVel3D_t Det_wrt_Earth;
  PosVel3D_t PtoleOrbit;
  PosVel3D_t Spin_z, Spin_xy;
  REAL8 eZ[3];

  if ( !posvel || !tGPS || !site || !edat ) {
    XLALPrintError ( "%s: Illegal NULL pointer passed!\n", fn);
    XLAL_ERROR( fn, XLAL_EINVAL );
  }

  /* ----- find ephemeris-based position of Earth wrt to SSB at this moment */
  LALBarycenterEarth( &status, &earth, tGPS, edat );
  if ( status.statusCode != 0 ) {
    XLALPrintError ( "%s: call to LALBarycenterEarth() failed!\n\n", fn);
    XLAL_ERROR( fn, XLAL_EFUNC );
  }
  /* ----- find ephemeris-based position of detector wrt to SSB */
  baryinput.tgps = *tGPS;
  baryinput.site = *site;
  baryinput.site.location[0] /= LAL_C_SI; baryinput.site.location[1] /= LAL_C_SI; baryinput.site.location[2] /= LAL_C_SI;
  baryinput.alpha = 0; baryinput.delta = 0; baryinput.dInv = 0;
  status = empty_status;
  LALBarycenter ( &status, &emit, &baryinput, &earth );
  if ( status.statusCode != 0 ) {
    XLALPrintError ( "%s: call to LALBarycenter() failed!\n\n", fn);
    XLAL_ERROR( fn, XLAL_EFAILED );
  }

  /* ----- determine position-vector of detector wrt center of Earth */
  Det_wrt_Earth.pos[0] = emit.rDetector[0] - earth.posNow[0];
  Det_wrt_Earth.pos[1] = emit.rDetector[1] - earth.posNow[1];
  Det_wrt_Earth.pos[2] = emit.rDetector[2] - earth.posNow[2];

  Det_wrt_Earth.vel[0] = emit.vDetector[0] - earth.velNow[0];
  Det_wrt_Earth.vel[1] = emit.vDetector[1] - earth.velNow[1];
  Det_wrt_Earth.vel[2] = emit.vDetector[2] - earth.velNow[2];

  eZ[0] = 0; eZ[1] = -siniEcl; eZ[2] = cosiEcl; 	/* ecliptic z-axis in equatorial coordinates */
  /* compute ecliptic-z projected spin motion */
  REAL8 pz = SCALAR ( Det_wrt_Earth.pos, eZ );
  REAL8 vz = SCALAR ( Det_wrt_Earth.vel, eZ );

  COPY_VECT ( Spin_z.pos, eZ );
  MULT_VECT ( pz, Spin_z.pos );

  COPY_VECT ( Spin_z.vel, eZ );
  MULT_VECT ( vz, Spin_z.vel );

  /* compute ecliptic-xy projected spin motion */
  COPY_VECT ( Spin_xy.pos, Det_wrt_Earth.pos );
  SUB_VECT ( Spin_xy.pos, Spin_z.pos );

  COPY_VECT ( Spin_xy.vel, Det_wrt_Earth.vel );
  SUB_VECT ( Spin_xy.vel, Spin_z.vel );

  /* ----- Ptolemaic special case: orbital motion on a circle */
  if ( (special == DETMOTION_SPIN_PTOLEORBIT) || (special == DETMOTION_PTOLEORBIT) )
    {
      if ( XLALPtolemaicPosVel ( &PtoleOrbit, tGPS ) ) {
	XLALPrintError ( "%s: call to XLALPtolemaicPosVel() failed!\n\n", fn);
	XLAL_ERROR( fn, XLAL_EFUNC );
      }
    }

  /* ----- return the requested type of detector motion */
  switch ( special )
    {
      /* full detector-motion: ephemeris-orbital + Earth-spin */
    case DETMOTION_SPIN_ORBIT:
      COPY_VECT(posvel->pos, emit.rDetector);
      COPY_VECT(posvel->vel, emit.vDetector);
      break;

      /* full ephemeris orbital detector-motion, neglecting Earth-spin */
    case DETMOTION_ORBIT:
      COPY_VECT(posvel->pos, earth.posNow);
      COPY_VECT(posvel->vel, earth.velNow);
      break;

      /* detector-motion including only Earth-spin, no orbital motion */
    case DETMOTION_SPIN:
      COPY_VECT(posvel->pos,Det_wrt_Earth.pos);
      COPY_VECT(posvel->vel,Det_wrt_Earth.vel);
      break;

      /* pure orbital detector motion, using "Ptolemaic" (ie. circular) approximation */
    case DETMOTION_PTOLEORBIT:
      COPY_VECT(posvel->pos,PtoleOrbit.pos);
      COPY_VECT(posvel->vel,PtoleOrbit.vel);
      break;

      /* Ptolemaic-orbital motion, plus Earth spin */
    case DETMOTION_SPIN_PTOLEORBIT:
      COPY_VECT(posvel->pos,PtoleOrbit.pos);
      COPY_VECT(posvel->vel,PtoleOrbit.vel);
      posvel->pos[0] += Det_wrt_Earth.pos[0];
      posvel->pos[1] += Det_wrt_Earth.pos[1];
      posvel->pos[2] += Det_wrt_Earth.pos[2];

      posvel->vel[0] += Det_wrt_Earth.vel[0];
      posvel->vel[1] += Det_wrt_Earth.vel[1];
      posvel->vel[2] += Det_wrt_Earth.vel[2];
      /*
      printf ("\nPtole = [ %f, %f, %f ], Ephem = [%f, %f, %f]\n",
	      posvel->pos[0], posvel->pos[1], posvel->pos[2],
	      emit.rDetector[0], emit.rDetector[1], emit.rDetector[2] );
      */
      break;

      /**< orbital motion plus *only* z-component of Earth spin-motion wrt to ecliptic plane */
    case DETMOTION_ORBIT_SPINZ:

      COPY_VECT(posvel->pos, earth.posNow);
      COPY_VECT(posvel->vel, earth.velNow);

      ADD_VECT ( posvel->pos, Spin_z.pos );
      ADD_VECT ( posvel->vel, Spin_z.vel );

      break;

      /**< orbital motion plus *only* x+y component of Earth spin-motion in the ecliptic */
    case DETMOTION_ORBIT_SPINXY:

      COPY_VECT(posvel->pos, earth.posNow);
      COPY_VECT(posvel->vel, earth.velNow);

      ADD_VECT ( posvel->pos, Spin_xy.pos );
      ADD_VECT ( posvel->vel, Spin_xy.vel );

      break;

    default:
      XLALPrintError("\n%s: Illegal 'special' value passed: '%d'\n\n", fn, special );
      XLAL_ERROR( fn, XLAL_EINVAL );
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

} /* XLALDetectorPosition() */



/** Compute position and velocity assuming a purely "Ptolemaic" orbital motion
 * (i.e. on a circle) around the sun, approximating Earth's orbit
 */
int
XLALPtolemaicPosVel ( PosVel3D_t *posvel,		/**< [out] instantaneous position and velocity vector */
		      const LIGOTimeGPS *tGPS		/**< [in] GPS time */
		      )
{
  const CHAR *fn = "XLALPtolemaicPosVel()";
  PulsarTimesParamStruc times = empty_PulsarTimesParamStruc;
  REAL8 phiOrb;   /* Earth orbital revolution angle, in radians. */
  REAL8 sinOrb, cosOrb;
  LALStatus status = empty_status;

  if ( !posvel || !tGPS ) {
    XLALPrintError ( "%s: Illegal NULL pointer passed!\n", fn);
    XLAL_ERROR( fn, XLAL_EINVAL );
  }

  times.epoch = (*tGPS);	/* get tAutumn */
  LALGetEarthTimes ( &status, &times );
  if ( status.statusCode ) {
    XLALPrintError ( "%s: call to LALGetEarthTimes() failed!\n\n", fn);
    XLAL_ERROR( fn, XLAL_EFUNC );
  }

  phiOrb = - LAL_TWOPI * times.tAutumn / LAL_YRSID_SI;
  sinOrb = sin(phiOrb);
  cosOrb = cos(phiOrb);

  /* Get instantaneous position. */
  posvel->pos[0] = rOrb_c * cosOrb;
  posvel->pos[1] = rOrb_c * sinOrb * cosiEcl;
  posvel->pos[2]=  rOrb_c * sinOrb * siniEcl;

  /* Get instantaneous velocity. */
  posvel->vel[0] = -vOrb_c * sinOrb;
  posvel->vel[1] =  vOrb_c * cosOrb * cosiEcl;
  posvel->vel[2] =  vOrb_c * cosOrb * siniEcl;

  return XLAL_SUCCESS;

} /* XLALPtolemaicPosVel() */



/** Compute a pure phase-deriv covariance [phi_i, phi_j] = <phi_i phi_j> - <phi_i><phi_j>
 * which gives a component of the "phase metric"
 */
double
CWPhase_cov_Phi_ij ( const intparams_t *params, double* relerr_max )
{
  const CHAR *fn = "CWPhase_cov_Phi_ij()";
  gsl_function integrand;

  intparams_t par = (*params);	/* struct-copy, as the 'deriv' field has to be changeable */
  int stat;

  /* sanity-check: don't allow any AM-coeffs being turned on here! */
  if ( par.amcomp1 != AMCOMP_NONE || par.amcomp2 != AMCOMP_NONE ) {
    XLALPrintError ( "%s: Illegal input, amcomp[12] must be set to AMCOMP_NONE!\n", fn );
    XLAL_ERROR_REAL8( fn, XLAL_EINVAL );
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
  size_t neval;
  double av_ij = 0, av_i = 0, av_j = 0;
  double av_ij_err = 0, av_i_err = 0, av_j_err = 0;

  for (n=0; n < Nseg; n ++ )
    {
      REAL8 ti = 1.0 * n * dT;
      REAL8 tf = MYMIN( (n+1.0) * dT, 1.0 );
      double res_n;

      /* compute <phi_i phi_j> */
      integrand.function = &CW_am1_am2_Phi_i_Phi_j;
      XLAL_CALLGSL ( stat = gsl_integration_qng (&integrand, ti, tf, epsabs, epsrel, &res_n, &abserr, &neval) );
      if ( outputIntegrand ) printf ("\n");
      if ( stat != 0 ) {
        XLALPrintError ( "\n%s: GSL-integration 'gsl_integration_qng()' of <Phi_i Phi_j> failed! seg=%d, av_ij_n=%g, abserr=%g, neval=%d\n",
                         fn, n, res_n, abserr, neval);
        XLAL_ERROR_REAL8( fn, XLAL_EFUNC );
      }
      av_ij += res_n;
      av_ij_err += SQUARE (abserr);


      /* compute <phi_i> */
      integrand.function = &CWPhaseDeriv_i;
      par.deriv = par.deriv1;
      XLAL_CALLGSL ( stat = gsl_integration_qng (&integrand, ti, tf, epsabs, epsrel, &res_n, &abserr, &neval) );
      if ( stat != 0 ) {
        XLALPrintError ( "\n%s: GSL-integration 'gsl_integration_qng()' of <Phi_i> failed! seg=%d, av_i_n=%g, abserr=%g, neval=%d\n",
                         fn, n, res_n, abserr, neval);
        XLAL_ERROR_REAL8( fn, XLAL_EFUNC );
      }
      av_i += res_n;
      av_i_err += SQUARE (abserr);

      /* compute <phi_j> */
      integrand.function = &CWPhaseDeriv_i;
      par.deriv = par.deriv2;
      XLAL_CALLGSL ( stat = gsl_integration_qng (&integrand, ti, tf, epsabs, epsrel, &res_n, &abserr, &neval) );
      if ( stat != 0 ) {
        XLALPrintError ( "\n%s: GSL-integration 'gsl_integration_qng()' of <Phi_j> failed! seg=%d, av_j_n=%g, abserr=%g, neval=%d\n",
                         fn, n, res_n, abserr, neval);
        XLAL_ERROR_REAL8( fn, XLAL_EFUNC );
      }
      av_j += res_n;
      av_j_err += SQUARE (abserr);

    } /* for i < Nseg */

  av_ij_err = sqrt(av_ij_err) / fabs(av_ij);
  av_i_err = sqrt(av_i_err) / fabs (av_i);
  av_j_err = sqrt(av_j_err) / fabs (av_j);

  maxrelerr = MYMAX ( av_ij_err, av_i_err );
  maxrelerr = MYMAX ( maxrelerr, av_j_err );

  if ( relerr_max )
    (*relerr_max) = maxrelerr;

  return ( av_ij - av_i * av_j );	/* return covariance */

} /* CWPhase_cov_Phi_ij() */


/** Calculate an approximate "phase-metric" with the specified parameters.
 *
 * The phase metric can only be computed for a single detector, if you want a
 * multi-detector metric you need to use the full Fstat-metric, as computed
 * by XLALDopplerFstatMetric().
 *
 * Note: if this function is called with multiple detectors, we compute the
 * phase metric using the *first* detector in the list!
 *
 *
 * Return NULL on error.
 */
gsl_matrix *
XLALDopplerPhaseMetric ( const DopplerMetricParams *metricParams,  	/**< input parameters determining the metric calculation */
			 const EphemerisData *edat,			/**< ephemeris data */
                         double *relerr_max
			 )
{
  const CHAR *fn = "XLALDopplerPhaseMetric()";
  gsl_matrix *g_ij = NULL;
  intparams_t intparams = empty_intparams;
  UINT4 i, j;
  REAL8 gg;
  const LALDetector *ifo;
  UINT4 dim;
  const LIGOTimeGPS *refTime, *startTime;
  const DopplerCoordinateSystem *coordSys;

  /* ---------- sanity/consistency checks ---------- */
  if ( !metricParams || !edat ) {
    XLALPrintError ("\n%s: Illegal NULL pointer passed.\n\n", fn);
    XLAL_ERROR_NULL( fn, XLAL_EINVAL );
  }

  startTime = &(metricParams->startTime);
  refTime   = &(metricParams->signalParams.Doppler.refTime);

  dim = metricParams->coordSys.dim;
  coordSys = &(metricParams->coordSys);

  /* ---------- prepare output metric ---------- */
  if ( (g_ij = gsl_matrix_calloc ( dim, dim )) == NULL ) {
    XLALPrintError ("%s: gsl_matrix_calloc(%d, %d) failed.\n\n", fn, dim, dim );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  /* always use first dector in list */
  ifo = &(metricParams->detInfo.sites[0]);

  /* ---------- set up integration parameters ---------- */
  intparams.edat = edat;
  intparams.startTime = XLALGPSGetREAL8 ( startTime );
  intparams.refTime   = XLALGPSGetREAL8 ( refTime );
  intparams.Tspan     = metricParams->Tspan;
  intparams.dopplerPoint = &(metricParams->signalParams.Doppler);
  intparams.detMotionType = metricParams->detMotionType;
  intparams.site = ifo;
  /* deactivate antenna-patterns for phase-metric */
  intparams.amcomp1 = AMCOMP_NONE;
  intparams.amcomp2 = AMCOMP_NONE;

  /* ---------- compute components of the phase-metric ---------- */
  double maxrelerr = 0, err;
  for ( i=0; i < dim; i ++ )
    {
      for ( j = 0; j <= i; j ++ )
	{
	  /* g_ij */
	  intparams.deriv1 = coordSys->coordIDs[i];
	  intparams.deriv2 = coordSys->coordIDs[j];
	  gg = CWPhase_cov_Phi_ij ( &intparams, &err );	/* [Phi_i, Phi_j] */
          maxrelerr = MYMAX ( maxrelerr, err );
	  if ( xlalErrno ) {
	    XLALPrintError ("\n%s: Integration of g_ij (i=%d, j=%d) failed. errno = %d\n", fn, i, j, xlalErrno );
            xlalErrno = 0;
            BOOLEAN sav = outputIntegrand;
            outputIntegrand = 1;
            gg = CWPhase_cov_Phi_ij ( &intparams, &maxrelerr );	/* [Phi_i, Phi_j] */
            outputIntegrand = sav;

	    XLAL_ERROR_NULL( fn, XLAL_EFUNC );
	  }
	  gsl_matrix_set (g_ij, i, j, gg);
	  gsl_matrix_set (g_ij, j, i, gg);

	} /* for j <= i */

    } /* for i < dim */

  if ( relerr_max )
    (*relerr_max) = maxrelerr;

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
  const CHAR *fn = "XLALDopplerFstatMetric()";
  DopplerMetric *metric = NULL;
  REAL8 cosi, psi;
  double relerr;
  gsl_matrix *tmp;

  /* ---------- sanity/consistency checks ---------- */
  if ( !metricParams || !edat ) {
    XLALPrintError ("%s: Illegal NULL pointer passed!\n\n", fn);
    XLAL_ERROR_NULL( fn, XLAL_EINVAL );
  }

  if ( metricParams->metricType >= METRIC_TYPE_LAST ) {
    XLALPrintError ("%s: Invalid value '%d' for metricType received. Must be within [%d,%d]!\n\n",
                    metricParams->metricType, 0, METRIC_TYPE_LAST - 1, fn);
    XLAL_ERROR_NULL( fn, XLAL_EINVAL );
  }

  /* if we're asked to compute F-metric only (1) or F-metric + phase-metric (2) */
  if ( metricParams->metricType == METRIC_TYPE_FSTAT || metricParams->metricType == METRIC_TYPE_ALL )
    {
      FmetricAtoms_t *atoms = NULL;

      /* ---------- compute Fmetric 'atoms', ie the averaged <a^2>, <a b Phi_i>, <a^2 Phi_i Phi_j>, etc ---------- */
      if ( (atoms = XLALComputeAtomsForFmetric ( metricParams, edat )) == NULL ) {
        XLALPrintError ("%s: XLALComputeAtomsForFmetric() failed. errno = %d\n\n", fn, xlalErrno );
        XLAL_ERROR_NULL( fn, XLAL_EFUNC );
      }

      /* ----- compute the F-metric gF_ij and related matrices ---------- */
      cosi = metricParams->signalParams.Amp.cosi;
      psi  = metricParams->signalParams.Amp.psi;

      if ( (metric = XLALComputeFmetricFromAtoms ( atoms, cosi, psi)) == NULL ) {
        XLALPrintError ("%s: XLALComputeFmetricFromAtoms() failed, errno = %d\n\n", fn, xlalErrno );
        XLALDestroyFmetricAtoms ( atoms );
        XLAL_ERROR_NULL( fn, XLAL_EFUNC );
      }

      /* ----- compute the full 4+n dimensional Fisher matrix ---------- */
      if ( (metric->Fisher_ab = XLALComputeFisherFromAtoms ( atoms, &metricParams->signalParams.Amp )) == NULL ) {
        XLALPrintError ("%s: XLALComputeFisherFromAtoms() failed. errno = %d\n\n", xlalErrno );
        XLALDestroyFmetricAtoms ( atoms );
        XLALDestroyDopplerMetric ( metric );
        XLAL_ERROR_NULL( fn, XLAL_EFUNC );
      }

      XLALDestroyFmetricAtoms ( atoms );
    } /* if compute F-metric */

  if ( metricParams->metricType == METRIC_TYPE_PHASE || metricParams->metricType == METRIC_TYPE_ALL )
    {
      /* if return-container 'metric' hasn't been allocated already earlier, we do it now */
      if (!metric ) {
        if ( (metric = XLALCalloc ( 1, sizeof(*metric) )) == NULL ) {
          XLALPrintError ("%s: XLALCalloc ( 1, %d) failed.\n\n", sizeof(*metric) );
          XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
        }
      }

      /* ----- compute the standard phase-metric g_ij ---------- */
      if ( (metric->g_ij = XLALDopplerPhaseMetric ( metricParams, edat, &relerr )) == NULL ) {
        XLALPrintError ("%s: XLALDopplerPhaseMetric() failed, errno = %d.\n\n", fn, xlalErrno );
        XLALDestroyDopplerMetric ( metric );
        XLAL_ERROR_NULL( fn, XLAL_EFUNC );
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
            XLALPrintError ("%s: failed to project gF_ij onto coordinate '%d'. errno=%d\n", fn, projCoord, xlalErrno );
            XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
          }
          gsl_matrix_free ( metric->gF_ij );
          metric->gF_ij = tmp;
        }

      /* gFav_ij */
      if ( metric->gFav_ij )
        {
          if ( (tmp = XLALProjectMetric ( metric->gFav_ij, projCoord )) == NULL ) {
            XLALPrintError ("%s: failed to project gFav_ij onto coordinate '%d'. errno=%d\n", fn, projCoord, xlalErrno );
            XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
          }
          gsl_matrix_free ( metric->gFav_ij );
          metric->gFav_ij = tmp;
        }

      /* phase-metric g_ij */
      if ( metric->g_ij )
        {
          if ( (tmp = XLALProjectMetric ( metric->g_ij, projCoord )) == NULL ) {
            XLALPrintError ("%s: failed to project g_ij onto coordinate '%d'. errno=%d\n", fn, projCoord, xlalErrno );
            XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
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
  const CHAR *fn = "XLALComputeAtomsForFmetric()";

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
    XLALPrintError ("\n%s: Illegal NULL pointer passed!\n\n", fn);
    XLAL_ERROR_NULL( fn, XLAL_EINVAL );
  }

  startTime = &(metricParams->startTime);
  refTime   = &(metricParams->signalParams.Doppler.refTime);

  dim = metricParams->coordSys.dim;	/* shorthand: number of Doppler dimensions */
  numDet = metricParams->detInfo.length;
  coordSys = &(metricParams->coordSys);

  /* ----- create output structure ---------- */
  if ( (ret = XLALCreateFmetricAtoms ( dim )) == NULL ) {
    XLALPrintError ("%s: call to XLALCreateFmetricAtoms (%s) failed. errno = %d\n\n", fn, dim, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  /* ---------- set up integration parameters ---------- */
  intparams.detMotionType = metricParams->detMotionType;
  intparams.dopplerPoint = &metricParams->signalParams.Doppler;
  intparams.startTime = XLALGPSGetREAL8 ( startTime );
  intparams.refTime   = XLALGPSGetREAL8 ( refTime );
  intparams.Tspan = metricParams->Tspan;
  intparams.edat = edat;

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

  /* FIXME: should probably not be hardcoded */
  if ( max_relerr > relerr_thresh )
    {
      XLALPrintError ("Maximal relative F-metric error too high: %.2e > %.2e\n", max_relerr, relerr_thresh );
      XLALDestroyFmetricAtoms ( ret );
      XLAL_ERROR_NULL( fn, XLAL_EFUNC );
    }
  else
    XLALPrintInfo ("\nMaximal relative error in F-metric: %.2e\n", max_relerr );

  return ret;

 failed:
  XLALDestroyFmetricAtoms ( ret );
  XLALPrintError ( "%s: XLALAverage_am1_am2_Phi_i_Phi_j() FAILED with errno = %d: am1 = %d, am2 = %d, i = %d : '%s', j = %d : '%s'\n",
		   fn, xlalErrno, intparams.amcomp1,intparams.amcomp2, i, DopplerCoordinateNames[i], j, DopplerCoordinateNames[j] );
  XLAL_ERROR_NULL( fn, XLAL_EFUNC );

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
DetectorMotionType
XLALParseDetectorMotionString ( const CHAR *detMotionString )
{
  const CHAR *fn = "XLALParseDetectorMotionString()";
  UINT4 i;

  if ( ! detMotionString ) {
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  for ( i=0; i < DETMOTION_LAST; i ++ )
    {
      if ( strcmp ( detMotionString, DetectorMotionNames[i] ) )
	continue;
      return i;	/* found the right entry */
    }

  XLALPrintError ("\nCould not parse '%s' into a valid detector-motion type!\n\n", detMotionString );
  XLAL_ERROR ( fn, XLAL_EINVAL );

} /* XLALParseDetectorMotionString() */


/** Provide a pointer to a static string containing the DopplerCoordinate-name
 * cooresponding to the enum DopplerCoordinateID
 */
const CHAR *
XLALDetectorMotionName ( DetectorMotionType detType )
{
  const CHAR *fn = "XLALDetectorMotionName()";

  if ( detType >= DETMOTION_LAST ) {
    XLALPrintError ( "%s: detector-motion type '%d' outside valid range [0, %d]\n\n", detType, DETMOTION_LAST - 1 );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  return ( DetectorMotionNames[detType] );

} /* XLALDetectorMotionName() */



/** Parse a DopplerCoordinate-name into the corresponding DopplerCoordinateID
 */
DopplerCoordinateID
XLALParseDopplerCoordinateString ( const CHAR *coordName )
{
  const CHAR *fn = "XLALParseDopplerCoordinateString()";
  UINT4 i;

  if ( !coordName )
    XLAL_ERROR ( fn, XLAL_EINVAL );

  for ( i=0; i < DOPPLERCOORD_LAST; i ++ )
    {
      if ( strcmp ( coordName, DopplerCoordinateNames[i] ) )
	continue;
      return i;	/* found the right entry */
    }

  XLALPrintError ("\nCould not parse '%s' into a valid coordinate-ID!\n\n", coordName );
  XLAL_ERROR ( fn, XLAL_EINVAL );

} /* XLALParseDopplerCoordinateString() */

/** Given a LALStringVector of coordinate-names, parse them into a
 * 'DopplerCoordinateSystem', namely a list of coordinate-IDs
 */
int
XLALDopplerCoordinateNames2System ( DopplerCoordinateSystem *coordSys,	/**< [out] pointer to coord-system struct */
				    const LALStringVector *coordNames 	/**< [in] list of coordinate names */
				    )
{
  const CHAR *fn = "XLALDopplerCoordinateNames2System()";
  UINT4 i;

  if ( !coordSys || !coordNames )
    XLAL_ERROR ( fn, XLAL_EINVAL );

  coordSys->dim = coordNames->length;
  for ( i=0; i < coordNames->length; i++ )
    {
      coordSys->coordIDs[i] = XLALParseDopplerCoordinateString ( coordNames->data[i] );
      if ( xlalErrno )
	XLAL_ERROR ( fn, XLAL_EFUNC );
    }

  return XLAL_SUCCESS;

} /* XLALDopplerCoordinateNames2System() */



/** Provide a pointer to a static string containing the DopplerCoordinate-name
 * cooresponding to the enum DopplerCoordinateID
 */
const CHAR *
XLALDopplerCoordinateName ( DopplerCoordinateID coordID )
{
  const CHAR *fn = "XLALDopplerCoordinateName()";

  if ( coordID >= DOPPLERCOORD_LAST ) {
    XLALPrintError ( "%s: coordID '%d' outside valid range [0, %d]\n\n", coordID, DOPPLERCOORD_LAST - 1 );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  return ( DopplerCoordinateNames[coordID] );

} /* XLALDopplerCoordinateName() */


/** Provide a pointer to a static string containing the a descriptive
 * 'help-string' describing the coordinate DopplerCoordinateID
 */
const CHAR *
XLALDopplerCoordinateHelp ( DopplerCoordinateID coordID )
{
  const CHAR *fn = "XLALDopplerCoordinateHelp()";

  if ( coordID >= DOPPLERCOORD_LAST ) {
    XLALPrintError ( "%s: coordID '%d' outside valid range [0, %d]\n\n", coordID, DOPPLERCOORD_LAST - 1 );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  return ( DopplerCoordinateNamesHelp[coordID] );

} /* XLALDopplerCoordinateHelp() */

/** Return a string (allocated here) containing a full name - helpstring listing
 * for all doppler-coordinates DopplerCoordinateID allowed by UniversalDopplerMetric.c
 */
CHAR *
XLALDopplerCoordinateHelpAll ( void )
{
  const CHAR *fn = "XLALDopplerCoordinateHelpAll()";
  CHAR *helpstr;
  const CHAR *name;
  const CHAR *help;
  UINT4 i, len, maxlen = 0;
  CHAR buf[512];
  CHAR fmt[512];

#define HEADER "Doppler-coordinate names and explanations:\n--------------------------------------------------\n"
  if ( (helpstr = XLALCalloc ( strlen(HEADER)+1, sizeof(CHAR) )) == NULL ) {
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  strcpy ( helpstr, HEADER );
  len = strlen(helpstr);

  /* get maximal field-length of coordinate names */
  for ( i = 0; i < DOPPLERCOORD_LAST; i ++ )
    {
      if ( (name = XLALDopplerCoordinateName ( i )) == NULL ) {
	XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
      }
      maxlen = MYMAX ( maxlen, strlen(name) );
      sprintf ( fmt, "%%-%ds: %%s\n", maxlen + 2 );
    }

  /* assemble help-lines */
  for ( i = 0; i < DOPPLERCOORD_LAST; i ++ )
    {
      if ( (name = XLALDopplerCoordinateName ( i )) == NULL ) {
	XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
      }
      if ( (help = XLALDopplerCoordinateHelp ( i )) == NULL ) {
	XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
      }

      snprintf ( buf, sizeof(buf) - 1, fmt, name, help );
      len += strlen ( buf ) + 1;
      if ( (helpstr = XLALRealloc ( helpstr, len )) == NULL ) {
	XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
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
XLALParseMultiDetectorInfo ( MultiDetectorInfo *detInfo,	/** [out] parsed detector-info struct */
			     const LALStringVector *detNames,	/**< [in] list of detector names */
			     const LALStringVector *detWeights	/**< [in] list of (strings) with detector weights (NULL if all 1) */
			     )
{
  const CHAR *fn = "XLALParseMultiDetectorInfo()";

  UINT4 X, numDet;
  REAL8 totalWeight;

  if ( !detInfo || !detNames || (detNames->length == 0) ) {
    XLALPrintError ("\n%s: Illegal NULL pointer input\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  numDet = detNames->length;
  if ( detWeights && (detWeights->length != numDet ) ) {
    XLALPrintError ("\n%s: Illegal input: number of noise-weights must agree with number of detectors\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
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
	XLALPrintError ("%s: Failed to get site-info for detector '%s'\n", fn, detNames->data[X] );
	XLAL_ERROR ( fn, XLAL_EINVAL );
      }
      detInfo->sites[X] = (*ifo);
      XLALFree ( ifo );

      /* parse noise weights if any */
      if ( detWeights )
	{
	  if ( 1 != sscanf ( detWeights->data[X], "%lf", &(detInfo->detWeights[X]) ) )
	    {
	      XLALPrintError ("%s: Failed to parse noise-weight '%s' into float.\n", fn, detWeights->data[X] );
	      XLAL_ERROR ( fn, XLAL_EINVAL );
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
  const CHAR *fn = "XLALCreateFmetricAtoms()";

  FmetricAtoms_t *ret;		/* output structure */

  if ( ( ret = XLALCalloc(1,sizeof(*ret))) == NULL ) {
    XLALPrintError ( "%s: XLALCalloc(1,%s) failed.\n", fn, sizeof(*ret));
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  if ( ( ret->a_a_i = gsl_vector_calloc (dim)) == NULL ) {
    XLALPrintError ( "%s: a_a_i = gsl_vector_calloc (%d) failed.\n", fn, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  if ( ( ret->a_b_i = gsl_vector_calloc (dim)) == NULL ) {
    XLALPrintError ( "%s: a_b_i = gsl_vector_calloc (%d) failed.\n", fn, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  if ( ( ret->b_b_i = gsl_vector_calloc (dim)) == NULL ) {
    XLALPrintError ( "%s: b_b_i = gsl_vector_calloc (%d) failed.\n", fn, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }


  if ( ( ret->a_a_i_j = gsl_matrix_calloc (dim, dim)) == NULL ) {
    XLALPrintError ( "%s: a_a_i_j = gsl_matrix_calloc (%d,%d) failed.\n", fn, dim, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  if ( ( ret->a_b_i_j = gsl_matrix_calloc (dim, dim)) == NULL ) {
    XLALPrintError ( "%s: a_b_i_j = gsl_matrix_calloc (%d,%d) failed.\n", fn, dim, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  if ( ( ret->b_b_i_j = gsl_matrix_calloc (dim, dim)) == NULL ) {
    XLALPrintError ( "%s: b_b_i_j = gsl_matrix_calloc (%d,%d) failed.\n", fn, dim, dim );
    XLALDestroyFmetricAtoms ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  return ret;

} /* XLALCreateFmetricAtoms() */


/** Convert amplitude-params from 'physical' coordinates {h0, cosi, psi, phi0} into
 * 'canonical' coordinates A^mu = {A1, A2, A3, A4}. The equations are found in
 * \ref JKS98 or \ref Prix07 Eq.(2).
 */
int
XLALAmplitudeParams2Vect ( PulsarAmplitudeVect *Amu, 		/**< [out] canonical amplitude coordinates A^mu = {A1, A2, A3, A4} */
			   const PulsarAmplitudeParams *Amp	/**< [in] 'physical' amplitude params {h0, cosi, psi, phi0} */
			   )
{
  const CHAR *fn = "XLALAmplitudeParams2Vect()";

  REAL8 Aplus, Across;
  REAL8 cosphi, sinphi, cos2psi, sin2psi;

  if ( !Amu || !Amp ) {
    XLALPrintError ("%s: illegal NULL pointer passed.\n\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  Aplus  = 0.5 * Amp->h0 * ( 1.0 + SQUARE(Amp->cosi) );
  Across = Amp->h0 * Amp->cosi;

  cosphi  = cos ( Amp->phi0 );
  sinphi  = sin ( Amp->phi0 );
  cos2psi = cos ( 2.0 * Amp->psi );
  sin2psi = sin ( 2.0 * Amp->psi );

  Amu->A1 =  Aplus * cosphi * cos2psi - Across * sinphi * sin2psi;
  Amu->A2 =  Aplus * cosphi * sin2psi + Across * sinphi * cos2psi;
  Amu->A3 = -Aplus * sinphi * cos2psi - Across * cosphi * sin2psi;
  Amu->A4 = -Aplus * sinphi * sin2psi + Across * cosphi * cos2psi;

  return XLAL_SUCCESS;

} /* XLALAmplitudeParams2Vect() */


/** Compute the 'F-metric' gF_ij (and also gFav_ij, m1_ij, m2_ij, m3_ij)
 * from the given FmetricAtoms and the signal amplitude parameters.
 *
 */
DopplerMetric*
XLALComputeFmetricFromAtoms ( const FmetricAtoms_t *atoms, REAL8 cosi, REAL8 psi )
{
  const CHAR *fn = "XLALComputeFmetricFromAtoms()";

  DopplerMetric *metric;		/* output matrix */

  UINT4 dim, i, j;			/* Doppler index counters */
  REAL8 A, B, C, D;			/* 'standard' antenna-pattern coefficients (gsl-integrated, though) */
  REAL8 alpha1, alpha2, alpha3, eta2, cos2psi, sin2psi;

  if ( !atoms ) {
    XLALPrintError ("%s: illegal NULL input.\n\n", fn );
    XLAL_ERROR_NULL (fn, XLAL_EINVAL );
  }

  if ( !atoms->a_a_i || !atoms->a_b_i || !atoms->b_b_i || !atoms->a_a_i_j || !atoms->a_b_i_j || !atoms->b_b_i_j ) {
    XLALPrintError ("%s: input Fmetric-atoms not fully allocated.\n\n", fn );
    XLAL_ERROR_NULL (fn, XLAL_EINVAL );
  }

  dim = atoms->a_a_i->size;

  /* allocate output metric structure */
  if ( (metric = XLALCalloc ( 1, sizeof(*metric) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc ( 1, %d) failed.\n\n", sizeof(*metric) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  metric->gF_ij = gsl_matrix_calloc ( dim, dim );
  metric->gFav_ij = gsl_matrix_calloc ( dim, dim );
  metric->m1_ij = gsl_matrix_calloc ( dim, dim );
  metric->m2_ij = gsl_matrix_calloc ( dim, dim );
  metric->m3_ij = gsl_matrix_calloc ( dim, dim );

  if ( !metric->gF_ij || !metric->gFav_ij || !metric->m1_ij || !metric->m2_ij || !metric->m3_ij ) {
    XLALPrintError ("%s: failed to gsl_matrix_calloc(%d,%d) for gF_ij, gFav_ij, m1_ij, m2_ij, m3_ij\n\n", fn, dim, dim );
    XLALDestroyDopplerMetric ( metric );
    XLAL_ERROR_NULL (fn, XLAL_ENOMEM );
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
XLALComputeFisherFromAtoms ( const FmetricAtoms_t *atoms, const PulsarAmplitudeParams *Amp )
{
  const CHAR *fn = "XLALComputeFisherFromAtoms()";
  gsl_matrix *fisher = NULL;	/* output matrix */

  UINT4 dimDoppler, dimFull, i, j;
  PulsarAmplitudeVect Amu;
  REAL8 al1, al2, al3;

  /* check input consistency */
  if ( !atoms || !Amp ) {
    XLALPrintError ("%s: illegal NULL input.\n\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  if ( !atoms->a_a_i || !atoms->a_b_i || !atoms->b_b_i || !atoms->a_a_i_j || !atoms->a_b_i_j || !atoms->b_b_i_j ) {
    XLALPrintError ("%s: input Fmetric-atoms not fully allocated.\n\n", fn );
    XLAL_ERROR_NULL (fn, XLAL_EINVAL );
  }

  if ( XLALAmplitudeParams2Vect ( &Amu, Amp ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: XLALAmplitudeParams2Vect() failed with errno = %d\n\n", fn, xlalErrno );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }

  dimDoppler = atoms->a_a_i->size;
  dimFull = 4 + dimDoppler;	/* 4 amplitude params + n Doppler params */

  if ( (fisher = gsl_matrix_calloc ( dimFull, dimFull )) == NULL ) {
    XLALPrintError ("%s: gsl_matric_calloc(%d,%d) failed.\n\n", fn, dimFull, dimFull );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
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
  al1 = SQUARE(Amu.A1) + SQUARE(Amu.A3);
  al2 = Amu.A1 * Amu.A2 + Amu.A3 * Amu.A4;
  al3 = SQUARE(Amu.A2) + SQUARE(Amu.A4);

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

      AR[0] =  Amu.A3 * a_a_i + Amu.A4 * a_b_i;
      AR[1] =  Amu.A3 * a_b_i + Amu.A4 * b_b_i;
      AR[2] = -Amu.A1 * a_a_i - Amu.A2 * a_b_i;
      AR[3] = -Amu.A1 * a_b_i - Amu.A2 * b_b_i;

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
  const char *fn = __func__;
  UINT4 i,j, dim1, dim2;
  gsl_matrix *ret_ij;

  if ( !g_ij ) {
    XLALPrintError ("%s: invalid NULL input 'g_ij'.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  dim1 = g_ij->size1;
  dim2 = g_ij->size2;

  if ( dim1 != dim2 ) {
    XLALPrintError ( "%s: input matrix g_ij must be square! (got %d x %d)\n", fn, dim1, dim2 );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  if ( (ret_ij = gsl_matrix_alloc ( dim1, dim2 )) == NULL ) {
    XLALPrintError ("%s: failed to gsl_matrix_alloc(%d, %d)\n", fn, dim1, dim2 );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
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

/** Convert 3-D vector from equatorial into ecliptic coordinates
 * return: 0 = OK, -1 = ERROR
 */
int
equatorialVect2ecliptic ( vect3D_t *out, vect3D_t * const in )
{
  static mat33_t rotEqu2Ecl = { { 1.0,        0,       0 },
                                { 0.0,  cosiEcl, siniEcl },
                                { 0.0, -siniEcl, cosiEcl } };
  if (!out || !in )
    return -1;

  return matrix33_in_vect3 ( out, &rotEqu2Ecl, in );

} /* equatorialVect2ecliptic() */

/** Convert 3-D vector from ecliptic into equatorial coordinates
 * return: 0 = OK, -1 = ERROR
 */
int
eclipticVect2equatorial ( vect3D_t *out, vect3D_t * const in )
{
  static mat33_t rotEcl2Equ =  { { 1.0,        0,       0 },
                                 { 0.0,  cosiEcl, -siniEcl },
                                 { 0.0,  siniEcl,  cosiEcl } };

  if (!out || !in )
    return -1;

  return matrix33_in_vect3 ( out, &rotEcl2Equ, in );

} /* eclipticVect2equatorial() */

/** compute matrix product mat . vect
 * return: 0 = OK, -1 = ERROR
 */
int
matrix33_in_vect3 ( vect3D_t *out, mat33_t * mat, vect3D_t * const in )
{
  if ( !out || !mat || !in )
    return -1;

  UINT4 i,j;
  for ( i=0; i < 3; i ++ )
    {
      (*out)[i] = 0;
      for ( j=0; j < 3; j ++ )
        {
          (*out)[i] += (*mat)[i][j] * (*in)[j];
        }
    }

  return 0;

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
  const char *fn = __func__;

  EarthState earth;
  LALStatus status;
  LIGOTimeGPS ti;
  INT4 h = 86400;	/* finite-differencing step-size: for rOrb, 1days seems like a reasonable choice */
  vect3D_t r0m2h, r0mh, r0, r0_h, r0_2h;
  vect3Dlist_t *ret = NULL;

  /* check input consistency */
#define MAX_MAXORDER 4
  if ( maxorder > MAX_MAXORDER ) {
    XLALPrintError ("%s: maxorder = %d too large, currently supports only up to maxorder = %d.\n", fn, maxorder, MAX_MAXORDER );
    XLAL_ERROR_NULL ( fn, XLAL_EDOM );
  }

  if ( !tGPS ) {
    XLALPrintError ("%s: invalid NULL pointer received for 'tGPS'.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }
  if ( !edat ) {
    XLALPrintError ("%s: invalid NULL pointer received for 'edat'.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }


  /* ----- find Earth's position at the 5 points: t0 -2h, t0-h, t0, t0+h, t0 + 2h ----- */

  /* t = t0 */
  ti = (*tGPS);
  status = empty_status;
  LALBarycenterEarth( &status, &earth, tGPS, edat );
  if ( status.statusCode != 0 ) {
    XLALPrintError ( "%s: call to LALBarycenterEarth() failed!\n\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  COPY_VECT ( r0, earth.posNow );

  /* t = t0 - h*/
  ti.gpsSeconds = (*tGPS).gpsSeconds - h;
  status = empty_status;
  LALBarycenterEarth( &status, &earth, tGPS, edat );
  if ( status.statusCode != 0 ) {
    XLALPrintError ( "%s: call to LALBarycenterEarth() failed!\n\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  COPY_VECT ( r0mh, earth.posNow );

  /* t = t0 - 2h*/
  ti.gpsSeconds = (*tGPS).gpsSeconds - 2 * h;
  status = empty_status;
  LALBarycenterEarth( &status, &earth, tGPS, edat );
  if ( status.statusCode != 0 ) {
    XLALPrintError ( "%s: call to LALBarycenterEarth() failed!\n\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  COPY_VECT ( r0m2h, earth.posNow );

  /* t = t0 + h*/
  ti.gpsSeconds = (*tGPS).gpsSeconds + h;
  status = empty_status;
  LALBarycenterEarth( &status, &earth, tGPS, edat );
  if ( status.statusCode != 0 ) {
    XLALPrintError ( "%s: call to LALBarycenterEarth() failed!\n\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  COPY_VECT ( r0_h, earth.posNow );

  /* t = t0 + 2h*/
  ti.gpsSeconds = (*tGPS).gpsSeconds + 2 * h;
  status = empty_status;
  LALBarycenterEarth( &status, &earth, tGPS, edat );
  if ( status.statusCode != 0 ) {
    XLALPrintError ( "%s: call to LALBarycenterEarth() failed!\n\n", fn);
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  COPY_VECT ( r0_2h, earth.posNow );

  /* use these 5 points to estimate derivatives */
  UINT4 i;
  vect3D_t rdn[MAX_MAXORDER+1];
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
    XLALPrintError ("%s: failed to XLALCalloc(1,%d)\n", fn, sizeof(*ret) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  if ( (ret->data = XLALCalloc ( maxorder + 1, sizeof(*ret->data) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc(%d,%d)\n", fn, maxorder + 1, sizeof(*ret->data) );
    XLALFree ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  UINT4 n;
  for ( n=0; n <= maxorder; n ++ ) {
    COPY_VECT ( ret->data[n], rdn[n] );
  }

  return ret;

} /* XLALComputeOrbitalDerivatives() */

