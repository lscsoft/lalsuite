/*
 * Copyright (C) 2005, 2006 Reinhard Prix
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
 * \author Reinhard Prix
 * \date 2005, 2006
 * \file 
 * \ingroup PulsarMetric
 * \brief Calculate flat approximation to the pulsar-metric.
 *
 * $Id$
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

#include <lal/DetectorStates.h>
#include <lal/PulsarTimes.h>

#include "FlatPulsarMetric.h"

NRCSID( FLATPULSARMETRICC, "$Id$");

/*---------- DEFINES ----------*/
#define TRUE (1==1)
#define FALSE (1==0)

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- empty initializers ---------- */

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/
REAL8 XLALcov ( const REAL8Vector *at, const REAL8Vector *bt );

/*==================== FUNCTION DEFINITIONS ====================*/

/** Compute a flat approximation for the continuous-wave metric (by neglecting the z-motion of 
 * the detector in ecliptic coordinates. 
 * 
 * gij has to be an allocated symmetric matrix of dimension \a dim: the order of coordinates
 * is \f$[ \omega_0, \tilde{k}_x, \tilde{k}_y, \omega_1, \omega_2, ... ]\f$, 
 * but \a dim must be at least 3 and maximally 6 (Freq + 2 sky + 3 spin-downs)
 * The dimensionless coordinates are defined as
 * \f$\omega_s \equiv 2\pi \, f^{(s)}\, T^{s+1}\f$ in terms of the observation time \f$T\f$,
 * and \f$\tilde{k}_l \equiv - 2\pi \bar{f} \hat{n} R_\mathrm{orb} / c\f$.
 * 
 * NOTE: the different detectors are currently combined simply by unit-weight averaging.
 */
int
XLALFlatMetricCW ( gsl_matrix *gij, 				/**< [out] metric */
		      const MultiDetectorStateSeries *mDetStates,/**< [in] detector-state series */
		      LIGOTimeGPS refTime			/**< [in] reference time for spin-parameters */
		      )
{
  UINT4 dim, numSpins, s, sp;
  UINT4 X, numDet;
  UINT4 i, l, numSteps;
  REAL8Vector *rx, *ry;		/* position time-series in ecliptic coordinates */
  REAL8Vector **ti;		/* time-series of deltaT^(s+1)/(Tspan (s+1)!) */

  REAL8 coseps, sineps;		/* inclination equatorial wrt ecliptic coords (eps ~23.4deg) */
  REAL8 Tspan;
  REAL8 gg;

  if ( !gij || !mDetStates || ( gij->size1 != gij->size2  ) )
    return -1;

  dim = gij->size1;
  if ( dim < 3 || dim > 6 )
    {
      LALPrintError ("\nMetric dimension must be 3 <= dim <= 6!\n\n");
      return -1;
    }
  numSpins = dim - 2;

  switch ( mDetStates->data[0]->system )
    {
    case COORDINATESYSTEM_ECLIPTIC:
      coseps = 1;
      sineps = 0;
      break;
    case COORDINATESYSTEM_EQUATORIAL:
      sineps = sin ( LAL_IEARTH );
      coseps = cos ( LAL_IEARTH );
      break;
    default:
      LALPrintError("\nIllegal coordinate system '%d' for XLALFlatMetricCW(): only EQUATORIAL or ECLIPTIC allowed!\n\n");
      return -1;
    } /* switch (system) */

  numDet = mDetStates->length;    
  numSteps = 0;
  for ( X=0; X < numDet; X ++ )
    numSteps += mDetStates->data[X]->length;

  Tspan = mDetStates->Tspan;
  
  if ( (rx = XLALCreateREAL8Vector ( numSteps ) ) == NULL )
    goto failed;
  if ( (ry = XLALCreateREAL8Vector ( numSteps ) ) == NULL )
    goto failed;
  ti = LALCalloc ( numSpins, sizeof(REAL8Vector*) );
  for ( s=0; s < numSpins; s ++ )
    if ( ( ti[s] = XLALCreateREAL8Vector ( numSteps )) == NULL )
      goto failed;


  /* determine the 2 + numSpins time-series {rx(t), ry(t), deltaT^{s+1}/(s+1)!} */
  l = 0;
  for ( X=0; X < numDet; X ++ )
    {
      const DetectorStateSeries *detStates = mDetStates->data[X];
      UINT4 numStepsX = detStates->length;
      for ( i=0; i < numStepsX; i ++ )
	{
	  /* detector position-series in ecliptic coordinates, normalized by 1AU */
	  rx->data[l] = (LAL_C_SI / LAL_AU_SI) * detStates->data[i].rDetector[0] ;
	  ry->data[i] = (LAL_C_SI / LAL_AU_SI) * ( coseps * detStates->data[i].rDetector[1]  + sineps * detStates->data[i].rDetector[2] );
	  /* spins */
	  ti[0]->data[l] = XLALGPSDiff ( &refTime, &(detStates->data[i].tGPS) )  / Tspan;
	  for ( s=1; s < numSpins; s ++ )
	    ti[s]->data[l] = ti[s-1]->data[l] * ti[0]->data[l] / (s + 1.0);

	  l ++ ;	/* counter over all time-steps */
	} /* for i < numStepsX */
    } /* for X < numDet */

  /* gxx */
  gg = XLALcov ( rx, rx );
  gsl_matrix_set (gij, 0, 0, gg);
  /* gyy */
  gg = XLALcov ( ry, ry );
  gsl_matrix_set (gij, 1, 1, gg);
  /* gxy */
  gg = XLALcov ( rx, ry );
  gsl_matrix_set (gij, 0, 1, gg);
  gsl_matrix_set (gij, 1, 0, gg);
  
  /* spins */
  for ( s=0; s < numSpins; s ++ )
    {
      gg = XLALcov ( rx, ti[s] );
      gsl_matrix_set (gij, 0, s+2, gg);
      gsl_matrix_set (gij, s+2, 0, gg);

      gg = XLALcov ( ry, ti[s] );
      gsl_matrix_set (gij, 1, s+2, gg);
      gsl_matrix_set (gij, s+2, 1, gg);

      for ( sp=s; sp < numSpins; s ++ )
	{
	  gg = XLALcov ( ti[s], ti[sp] );
	  gsl_matrix_set (gij, s+2,  sp+2, gg);
	  gsl_matrix_set (gij, sp+2, s+2, gg);
	} /* for s' in [s, numSpins) */

    } /* for s < numSpins */

  /* get metric from {kx, ky, w0, w1, w2,..} into 'canonical' order {w0, kx, ky, w1, w2, ... } */
  gsl_matrix_swap_rows (gij, 0, 2);	/* swap w0 <--> kx */
  gsl_matrix_swap_rows (gij, 1, 2);	/* swap kx <--> ky */
  gsl_matrix_swap_columns (gij, 0, 2);	/* swap w0 <--> kx */
  gsl_matrix_swap_columns (gij, 1, 2);	/* swap kx <--> ky */


  /* Ok */
  if ( rx ) XLALDestroyREAL8Vector ( rx);
  if ( ry ) XLALDestroyREAL8Vector ( ry);
  for ( s=0; s < numSpins; s ++ )
    if ( ti[s] )
      XLALDestroyREAL8Vector ( ti[s] );

  return 0;

 failed:
  if ( rx ) XLALDestroyREAL8Vector ( rx);
  if ( ry ) XLALDestroyREAL8Vector ( ry);
  for ( s=0; s < numSpins; s ++ )
    if ( ti[s] )
      XLALDestroyREAL8Vector ( ti[s] );
  
  return -1;
  

} /* XLALFlatMetricCW() */


/** Compute the 'covariance' between two time-series a(t), and b(t), namely
 * \f$\mathrm{cov}(a,b) = \langle a\, b\rangle - \langle a \rangle \langle b \rangle\f$, 
 * where \f$ \langle a \rangle = (1/T) \int_0^T a(t) \, d t\f$.
 */
REAL8
XLALcov ( const REAL8Vector *at, const REAL8Vector *bt )
{
  UINT4 i, N;
  REAL8 av_a, av_b, av_ab;
  const REAL8 *this_a, *this_b;

  if ( !at || !bt || at->length != bt->length ) {
    LALPrintError("\nXLALcov(): NULL input or time-series of different lengths!\n\n");
    XLAL_ERROR_REAL8 ( "XLALcov()", XLAL_EINVAL );
  }

  N = at->length;

  av_ab = av_a = av_b = 0;
  this_a = &at->data[0];
  this_b = &bt->data[0];

  for ( i=0; i <  N; i ++ )
    {
      REAL8 a = *this_a;
      REAL8 b = *this_b;

      av_ab += a * b;
      av_a  += a;
      av_b  += b;

      this_a ++;
      this_b ++;
    } /* for i < N */
  
  av_ab *= (1.0 / N);
  av_a  *= (1.0 / N);
  av_b  *= (1.0 / N);

  return ( av_ab - av_a * av_b );

} /* XLALcov() */



/* [OBSOLETE?] ================================================================================*/

/** [OBSOLETE?] Compute the constant-coefficient approximate pulsar-metric.
 *  The returned metric symmetric matrix \f$g_{\alpha \beta}\f$ 
 * is encoded in a 1D vector \f$v_l\f$ in the same
 *  way as in LALPulsarMetric(), namely: 
 * \f$g_{\alpha \beta} = v_l\f$, where for \f$\alpha<=\beta\f$ 
 * the vector-index \f$l\f$ is \f$l = \alpha + \beta \,(\beta + 1)/2\f$.
 *
 * The (dimensionless) parameter-space coordinates (and their order) are:
 * \f$\{ \kappa^X, \kappa^Y, \varpi_0, \varpi_1, \varpi_2, ...\}\f$, defined as
 * \f[ \kappa^i \equiv R_{ES} {2\pi \over c} f \, n^i\,, \f]
 * \f[ \varpi_s \equiv T^{s+1}\, 2\pi f^{(s)}\, \f] 
 * where \f$R_{ES} = 1\,\textrm{AU} \sim 1.5\times10^{11}\f$ m is the orbital radius,
 * \f$f\f$ is the frequency, \f$n^i\f$ is the unit-vector pointing to a sky-location, 
 * \f$T\f$ is the observation time and \f$f^{(s)}\f$ is the s-th time-derivative
 * of the frequency.
 * 
 * The two cartesian coordinate-directions \f$X, Y\f$ are referring to the 
 * <em>ecliptic</em> cartesian reference-frame.
 *
 */
void
LALFlatPulsarMetric ( LALStatus *status,
		      REAL8Vector **metric,	/**< [out] constant pulsar-metric */
		      LIGOTimeGPS startTime,	/**< start time of observation */
		      REAL8 duration,		/**< duration of observation in seconds */
		      const LALDetector *site 	/**< [in] detector location */
		      )
{
  REAL8Vector *ret = NULL;		/* return final metric */

  PulsarTimesParamStruc zero_phases; 	/* Needed to calculate initial phases */
  REAL8 lat, lon;			/* Detector latitude and longitude */

  UINT4 spinOrder = 3;			/* number of spin-parameters to include */
  UINT4 dim = 2 + spinOrder;		/* parameter-space dimension : 2 sky + spin-params */
  UINT4 s0, s1; 			/* spin-parameter indices */

  /* parameters */
  REAL8 phi_o_0, phi_s_0;       	/* Initial phases of orbital- and spin- motion */
  REAL8 phi_o_1, phi_s_1;       	/* Final phases of orbital- and spin- motion */
  REAL8 delta_phi_o, delta_phi_s;	/* phase-differences (final - intial) */
  REAL8 omega_o, omega_s;		/* orbital- and spin- rotation rates */
  REAL8 REX, REY;			/* detector spin-motion amplitudes in units of 1AU */

  /* intermediate results */
  REAL8 rX_av, rY_av;			/* averages of rX(t) and rY(t) over observation-time */
  REAL8 rX_2_av, rY_2_av;		/* averages of rX(t)^2 and rY(t)^2 over observation-time */
  REAL8 rX_rY_av;			/* average  of rX(t)*rY(t) over observation-time */
  REAL8 cosLambda;			/* lambda = detector's latitude */
  REAL8 cosEps;				/* Eps = ecliptic inclination of earth-spin (ca 23.4deg) */

  /*----------*/
  INITSTATUS( status, "LALFlatPulsarMetric", FLATPULSARMETRICC );
  ATTATCHSTATUSPTR (status);

  /* get detector's position */
  lon = site->frDetector.vertexLongitudeRadians;
  lat = site->frDetector.vertexLatitudeRadians;

  /* allocate memory for resulting metric components */
  TRY ( LALDCreateVector( status->statusPtr, &ret, dim * (dim+1) / 2 ), status);

  
  /* ---------- pure spin-spin components g_{s0 s1} ---------- */
  {
    UINT4 s0fact = 1, s1fact = 1;	/* factorials of s0 and s1 */

    for ( s0 = 0; s0 < spinOrder; s0 ++ )
      {
	s0fact *= (s0 ? s0 : 1);
	s1fact = 1;
	for ( s1 = s0; s1 < spinOrder; s1 ++ )
	  {
	    UINT4 tmp;
	    REAL8 gs0s1;

	    s1fact *= (s1 ? s1 : 1);
	    
	    tmp = s0fact * s1fact * (s0 + 2) * (s1 + 2) * (s0 + s1 + 3);
	    gs0s1 = 1.0 / tmp;

	    ret->data[ PMETRIC_INDEX (s0+2, s1+2) ] = gs0s1;
	
	  } /* for s0 <= s1 < spinOrder */
      }

  } /* spin-spin components */

  /* ---------- pure sky-sky components ---------- */

  /* amplitudes of detector spin-motion in ecliptic plane, in units of 1AU */
  cosLambda = cos(lat);
  cosEps = cos(LAL_IEARTH);
  REX = LAL_REARTH_SI * cosLambda / LAL_AU_SI;
  REY = REX * cosEps;

  /* Angular velocities: */
  omega_s = LAL_TWOPI / LAL_DAYSID_SI;
  omega_o = LAL_TWOPI / LAL_YRSID_SI;

  /* Calculation of phases of spin and orbit at start: */
  zero_phases.epoch = startTime;
  TRY (LALGetEarthTimes( status->statusPtr, &zero_phases), status);
  phi_o_0 = LAL_TWOPI - zero_phases.tAutumn/LAL_YRSID_SI*LAL_TWOPI;
  phi_s_0 = LAL_TWOPI - zero_phases.tMidnight/LAL_DAYSID_SI*LAL_TWOPI + lon;

  /* phases at end of observation-time */
  delta_phi_o = omega_o * duration;
  delta_phi_s = omega_s * duration;
  phi_o_1 = phi_o_0 + delta_phi_o;
  phi_s_1 = phi_s_0 + delta_phi_s;

  /* rX, rY are in units of 1AU ! */
  rX_av  = 
    1.0 * ( sin(phi_o_1) - sin(phi_o_0) ) / delta_phi_o + 
    REX * ( sin(phi_s_1) - sin(phi_s_0) ) / delta_phi_s;

  rY_av = 
    1.0 * (-cos(phi_o_1) + cos(phi_o_0) ) / delta_phi_o + 
    REY * (-cos(phi_s_1) + cos(phi_s_0) ) / delta_phi_s;

  /* calculate <rX^2>, <rY^2>, and <rX rY> */
  {
    /* intermediate values for making calculation of <ri rj>  more readable */
    REAL8 cos_o_2_av = 0.5 + (0.25 / delta_phi_o) * ( sin(2.0*phi_o_1) - sin(2.0*phi_o_0) );
    REAL8 cos_s_2_av = 0.5 + (0.25 / delta_phi_s) * ( sin(2.0*phi_s_1) - sin(2.0*phi_s_0) );

    REAL8 sin_o_2_av = 0.5 - (0.25 / delta_phi_o) * ( sin(2.0*phi_o_1) - sin(2.0*phi_o_0) );
    REAL8 sin_s_2_av = 0.5 - (0.25 / delta_phi_s) * ( sin(2.0*phi_s_1) - sin(2.0*phi_s_0) );

    REAL8 cos_o_cos_s_av = 
      ( 0.5 / (delta_phi_s + delta_phi_o) ) * ( sin(phi_s_1 + phi_o_1) - sin(phi_s_0 + phi_o_0) ) + 
      ( 0.5 / (delta_phi_s - delta_phi_o) ) * ( sin(phi_s_1 - phi_o_1) - sin(phi_s_0 - phi_o_0) );

    REAL8 sin_o_sin_s_av = 
      (-0.5 / (delta_phi_s + delta_phi_o) ) * ( sin(phi_s_1 + phi_o_1) - sin(phi_s_0 + phi_o_0) ) +
      ( 0.5 / (delta_phi_s - delta_phi_o) ) * ( sin(phi_s_1 - phi_o_1) - sin(phi_s_0 - phi_o_0) );

    REAL8 cos_o_sin_o_av = (-0.25 / delta_phi_o) * ( cos(2.0*phi_o_1) - cos(2.0*phi_o_0) );
    REAL8 cos_s_sin_s_av = (-0.25 / delta_phi_s) * ( cos(2.0*phi_s_1) - cos(2.0*phi_s_0) );
    
    REAL8 cos_o_sin_s_av = 
      (-0.5 / (delta_phi_s + delta_phi_o) ) * ( cos(phi_s_1 + phi_o_1) - cos(phi_s_0 + phi_o_0) ) +
      (-0.5 / (delta_phi_s - delta_phi_o) ) * ( cos(phi_s_1 - phi_o_1) - cos(phi_s_0 - phi_o_0) );

    /* simply exchange s <--> o for this next one: */
    REAL8 cos_s_sin_o_av = 
      (-0.5 / (delta_phi_o + delta_phi_s) ) * ( cos(phi_o_1 + phi_s_1) - cos(phi_o_0 + phi_s_0) ) +
      (-0.5 / (delta_phi_o - delta_phi_s) ) * ( cos(phi_o_1 - phi_s_1) - cos(phi_o_0 - phi_s_0) );


    /* therefore we can now write: */
    rX_2_av = cos_o_2_av  + REX*REX * cos_s_2_av     + 2.0 * REX * cos_o_cos_s_av;
    rY_2_av = sin_o_2_av  + REY*REY * sin_s_2_av     + 2.0 * REY * sin_o_sin_s_av;

    rX_rY_av= cos_o_sin_o_av + REX*REY*cos_s_sin_s_av + REX*cos_s_sin_o_av + REY*cos_o_sin_s_av;
  }

  /* ==> sky-sky metric components: */
  ret->data[ PMETRIC_INDEX(0,0) ] = rX_2_av  - rX_av * rX_av;	/* cov(rX, rX) */
  ret->data[ PMETRIC_INDEX(1,1) ] = rY_2_av  - rY_av * rY_av;	/* cov(rY, rY) */
  ret->data[ PMETRIC_INDEX(0,1) ] = rX_rY_av - rX_av * rY_av;	/* cos(rX, rY) */

  /* ---------- mixed sky-spin components ---------- */
  { /* < t^1 r_i > */
    REAL8 t_1_cos_o_av = (1.0 / pow(delta_phi_o,2)) * 
      ( delta_phi_o * sin(phi_o_1) + cos(phi_o_1) - cos(phi_o_0) );
    REAL8 t_1_cos_s_av = (1.0 / pow(delta_phi_s,2)) * 
      ( delta_phi_s * sin(phi_s_1) + cos(phi_s_1) - cos(phi_s_0) ); /* o -> s */

    REAL8 t_1_sin_o_av = (1.0 / pow(delta_phi_o,2)) * 
      (-delta_phi_o * cos(phi_o_1) + sin(phi_o_1) - sin(phi_o_0) );
    REAL8 t_1_sin_s_av = (1.0 / pow(delta_phi_s,2)) * 
      (-delta_phi_s * cos(phi_s_1) + sin(phi_s_1) - sin(phi_s_0) );

    REAL8 t_1_rX_av = t_1_cos_o_av + REX * t_1_cos_s_av;
    REAL8 t_1_rY_av = t_1_sin_o_av + REY * t_1_sin_s_av;

    REAL8 t_1_av = 1.0 / 2.0;

    ret->data[ PMETRIC_INDEX(0,2) ] = t_1_rX_av - t_1_av * rX_av;   /* g_w0_X */
    ret->data[ PMETRIC_INDEX(1,2) ] = t_1_rY_av - t_1_av * rY_av;   /* g_w0_Y */
  }
  {  /* < t^2 r_i > */
    REAL8 t_2_cos_o_av = (1.0 / pow(delta_phi_o,3)) * 
      ( (pow(delta_phi_o,2) - 2.0) * sin(phi_o_1) + 2.0*delta_phi_o * cos(phi_o_1) + 2.0*sin(phi_o_0));
    REAL8 t_2_cos_s_av = (1.0 / pow(delta_phi_s,3)) * 
      ( (pow(delta_phi_s,2) - 2.0) * sin(phi_s_1) + 2.0*delta_phi_s * cos(phi_s_1) + 2.0*sin(phi_s_0));
    
    REAL8 t_2_sin_o_av = (1.0 / pow(delta_phi_o,3)) * 
      (-(pow(delta_phi_o,2) - 2.0) * cos(phi_o_1) + 2.0*delta_phi_o * sin(phi_o_1) - 2.0*cos(phi_o_0));
    REAL8 t_2_sin_s_av = (1.0 / pow(delta_phi_s,3)) * 
      (-(pow(delta_phi_s,2) - 2.0) * cos(phi_s_1) + 2.0*delta_phi_s * sin(phi_s_1) - 2.0*cos(phi_s_0));


    REAL8 t_2_rX_av = t_2_cos_o_av + REX * t_2_cos_s_av;
    REAL8 t_2_rY_av = t_2_sin_o_av + REY * t_2_sin_s_av;

    REAL8 t_2_av = 1.0 / 3.0;

    ret->data[ PMETRIC_INDEX(0,3) ] = t_2_rX_av - t_2_av * rX_av;   /* g_w1_X */
    ret->data[ PMETRIC_INDEX(1,3) ] = t_2_rY_av - t_2_av * rY_av;   /* g_w1_Y */
  }
  { /* < t^3 r_i > */
    REAL8 t_3_cos_o_av = (1.0 / pow(delta_phi_o,4)) * 
      ( (pow(delta_phi_o,3) - 6.0*delta_phi_o) * sin(phi_o_1) 
	+ (3.0*pow(delta_phi_o,2)-6.0) * cos(phi_o_1) + 6.0*cos(phi_o_0) );
    REAL8 t_3_cos_s_av = (1.0 / pow(delta_phi_s,4)) * 
      ( (pow(delta_phi_s,3) - 6.0*delta_phi_s) * sin(phi_s_1) 
	+ (3.0*pow(delta_phi_s,2)-6.0) * cos(phi_s_1) + 6.0*cos(phi_s_0) );

    REAL8 t_3_sin_o_av = (1.0 / pow(delta_phi_o,4)) * 
      (-(pow(delta_phi_o,3) - 6.0*delta_phi_o) * cos(phi_o_1) 
	+ (3.0*pow(delta_phi_o,2)-6.0) * sin(phi_o_1) + 6.0*sin(phi_o_0) );
    REAL8 t_3_sin_s_av = (1.0 / pow(delta_phi_s,4)) * 
      (-(pow(delta_phi_s,3) - 6.0*delta_phi_s) * cos(phi_s_1) 
	+ (3.0*pow(delta_phi_s,2)-6.0) * sin(phi_s_1) + 6.0*sin(phi_s_0) );
    
    REAL8 t_3_rX_av = t_3_cos_o_av + REX * t_3_cos_s_av;
    REAL8 t_3_rY_av = t_3_sin_o_av + REY * t_3_sin_s_av;

    REAL8 t_3_av = 1.0 / 4.0;

    ret->data[ PMETRIC_INDEX(0,4) ] = t_3_rX_av - t_3_av * rX_av;   /* g_w2_X */
    ret->data[ PMETRIC_INDEX(1,4) ] = t_3_rY_av - t_3_av * rY_av;   /* g_w2_Y */
  }

  
  (*metric) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALFlatPulsarMetric() */


