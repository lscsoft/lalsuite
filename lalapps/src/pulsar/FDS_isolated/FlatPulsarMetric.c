/*
 * Copyright (C) 2005 Reinhard Prix
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
 * \date 2005
 * \file 
 * \ingroup flatPulsarMetric
 * \brief Calculate the flat approximation to the pulsar-metric.
 *
 * $Id$
 *
 */

/*---------- INCLUDES ----------*/
#include <math.h>

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

/*==================== FUNCTION DEFINITIONS ====================*/

/** Compute the constant-coefficient approximate pulsar-metric.
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
  REAL8 phi_o_i;          /* Phase of Earth's orbit at startTime */
  REAL8 phi_s_i;          /* Phase of Earth's spin at startTime */
  PulsarTimesParamStruc zero_phases; /* Needed to calculate initial phases */
  REAL8 lat, lon;	/* Detector latitude and longitude */
  REAL8Vector *ret = NULL;	/* return final metric */

  UINT4 spdnOrder = 3;	/* number of spindowns to include: now set to 3 */
  UINT4 dim = 2 + spdnOrder;	/* parameter-space dimension : 2 sky + spindowns */
  UINT4 s0, s1; 	/* spin-indices */
  
  INITSTATUS( status, "LALFlatPulsarMetric", FLATPULSARMETRICC );
  ATTATCHSTATUSPTR (status);

  /* get detector's position */
  lon = site->frDetector.vertexLongitudeRadians;
  lat = site->frDetector.vertexLatitudeRadians;

  /* Calculation of phases of spin and orbit at start: */
  zero_phases.epoch = startTime;
  TRY (LALGetEarthTimes( status->statusPtr, &zero_phases), status);
  phi_o_i = LAL_TWOPI - zero_phases.tAutumn/LAL_YRSID_SI*LAL_TWOPI;
  phi_s_i = LAL_TWOPI - zero_phases.tMidnight/LAL_DAYSID_SI*LAL_TWOPI + lon;

  /* allocate memory for resulting metric components */
  TRY ( LALDCreateVector( status->statusPtr, &ret, dim * (dim+1) / 2 ), status);

  
  /* ----- pure spin-spin components g_{s0 s1} ----- */
  {
    UINT4 s0fact = 1, s1fact = 1;	/* factorials of s0 and s1 */

    for ( s0 = 0; s0 < spdnOrder; s0 ++ )
      {
	s0fact *= (s0 ? s0 : 1);
	s1fact = 1;
	for ( s1 = s0; s1 < spdnOrder; s1 ++ )
	  {
	    UINT4 tmp;
	    REAL8 gs0s1;

	    s1fact *= (s1 ? s1 : 1);
	    
	    tmp = s0fact * s1fact * (s0 + 2) * (s1 + 2) * (s0 + s1 + 3);
	    gs0s1 = 1.0 / tmp;

	    ret->data[ PMETRIC_INDEX (s0+2, s1+2) ] = gs0s1;
	
	  } /* for s0 <= s1 < spdnOrder */
      }

  } /* spin-spin components */
  /* ----- pure sky-components ----- */

  /* ----- mixed sky-spin components ----- */



  
  (*metric) = ret;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALFlatPulsarMetric() */
