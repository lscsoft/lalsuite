/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Teviet Creighton, John Whelan
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

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SkyCoordinates.h>

#define LAL_DGALCORE_SI (2.62e20) /* Galactic core distance (metres) */

/**
 * \author Creighton, T. D.
 *
 * \brief Computes the input parameters for a PPN inspiral.
 *
 * ### Description ###
 *
 * This function takes a Galactic location and pair of masses from
 * <tt>*input</tt> and uses them to set the \c PPNParamStruc fields
 * <tt>output->position</tt>, <tt>output->mTot</tt>, <tt>output->eta</tt>, and
 * <tt>output->d</tt>.  The fields <tt>output->psi</tt>, <tt>output->inc</tt>,
 * and <tt>output->phi</tt> are set randomly to reflect a uniform
 * distribution in solid angle (that is, cosine of inclination is uniform
 * between \f$-1\f$ and 1, other angles are uniform between 0 and \f$2\pi\f$).
 * The routine uses the random sequence specified by <tt>*params</tt> when
 * given, but if <tt>*params</tt>=\c NULL a new sequence is started
 * internally using the current execution time as a seed. The field
 * <tt>input->geocentEndTime</tt> is ignored by this routine.
 *
 * The other \c PPNParamStruc input fields are not touched by this
 * routine, and must be specified externally before generating a waveform
 * with this structure.
 *
 * ### Algorithm ###
 *
 * Galactocentric Galactic axial coordinates \f$\rho\f$, \f$z\f$, and \f$l_G\f$ are
 * transformed to geocentric Galactic Cartesian coordinates:
 * \f{eqnarray}{
 * x_e & = & R_e + \rho\cos l_G \;,\\
 * y_e & = & \rho\sin l_G \;,\\
 * z_e & = & z \;,
 * \f}
 * where
 * \f[
 * R_e \approx 8.5\,\mathrm{kpc}
 * \f]
 * is the distance to the Galactic core (this constant will probably
 * migrate into \ref LALConstants.h eventually).  These are converted
 * to geocentric Galactic spherical coordinates:
 * \f{eqnarray}{
 * d & = & \sqrt{x_e^2 + y_e^2 + z_e^2} \;,\\
 * b & = & \arcsin\left(\frac{z_e}{d_e}\right) \;,\\
 * l & = & \arctan\!2(y_e,x_e) \;.
 * \f}
 * In the calculation of \f$d\f$ we factor out the leading order term from
 * the square root to avoid inadvertent overflow, and check for underflow
 * in case the location lies on top of the Earth.  The angular
 * coordinates are then transformed to equatorial celestial coordinates
 * \f$\alpha\f$ and \f$\delta\f$ using the routines in \ref SkyCoordinates.h.
 *
 */
void
LALGetInspiralParams( LALStatus                  *stat,
		      PPNParamStruc              *output,
		      GalacticInspiralParamStruc *input,
		      RandomParams               *params )
{
  REAL4 x, y, z;  /* geocentric Galactic Cartesian coordinates */
  REAL4 max, d;   /* maximum of x, y, and z, and normalized distance */
  REAL4 psi, phi, inc; /* polarization, phase, and inclination angles */
  REAL4 mTot;          /* total binary mass */
  SkyPosition direction; /* direction to the source */
  RandomParams *localParams = NULL; /* local random parameters pointer */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures exist. */
  ASSERT( output, stat, GENERATEPPNINSPIRALH_ENUL,
	  GENERATEPPNINSPIRALH_MSGENUL );
  ASSERT( input, stat, GENERATEPPNINSPIRALH_ENUL,
	  GENERATEPPNINSPIRALH_MSGENUL );

  /* Compute total mass. */
  mTot = input->m1 + input->m2;
  if ( mTot == 0.0 ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EMBAD,
	   GENERATEPPNINSPIRALH_MSGEMBAD );
  }

  /* Compute Galactic geocentric Cartesian coordinates. */
  x = LAL_DGALCORE_SI + input->rho*1e3*LAL_PC_SI*cos( input->lGal );
  y = input->rho*1e3*LAL_PC_SI*sin( input->lGal );
  z = input->z*1e3*LAL_PC_SI;

  /* Compute Galactic geocentric spherical coordinates. */
  max = fabs( x );
  if ( fabs( y ) > max )
    max = fabs( y );
  if ( fabs( z ) > max )
    max = fabs( z );
  if ( max == 0.0 ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EDBAD,
	   GENERATEPPNINSPIRALH_MSGEDBAD );
  }
  x /= max;
  y /= max;
  z /= max;
  d = sqrt( x*x + y*y + z*z );
  direction.latitude = asin( z/d );
  direction.longitude = atan2( y, x );
  direction.system = COORDINATESYSTEM_GALACTIC;

  /* Compute equatorial coordinates. */
  TRY( LALGalacticToEquatorial( stat->statusPtr, &direction,
				&direction ), stat );
  output->position = direction;
  output->d = max*d;

  /* If we haven't been given a random sequence, generate one. */
  if ( params )
    localParams = params;
  else {
    TRY( LALCreateRandomParams( stat->statusPtr, &localParams, 0 ),
	 stat );
  }

  /* Compute random inclination and polarization angle. */
  LALUniformDeviate( stat->statusPtr, &psi, localParams );
  if ( params )
    CHECKSTATUSPTR( stat );
  else
    BEGINFAIL( stat )
      TRY( LALDestroyRandomParams( stat->statusPtr, &localParams ),
	   stat );
      ENDFAIL( stat );
  LALUniformDeviate( stat->statusPtr, &phi, localParams );
  if ( params )
    CHECKSTATUSPTR( stat );
  else
    BEGINFAIL( stat )
      TRY( LALDestroyRandomParams( stat->statusPtr, &localParams ),
	   stat );
    ENDFAIL( stat );
  LALUniformDeviate( stat->statusPtr, &inc, localParams );
  if ( params )
    CHECKSTATUSPTR( stat );
  else
    BEGINFAIL( stat )
      TRY( LALDestroyRandomParams( stat->statusPtr, &localParams ),
	   stat );
    ENDFAIL( stat );
  output->psi = LAL_TWOPI*psi;
  output->phi = LAL_TWOPI*phi;
  inc = 2.0*inc - 1.0;
  output->inc = acos( inc );

  /* Set output masses. */
  output->mTot = mTot;
  output->eta = (input->m1/mTot)*(input->m2/mTot);

  /* Clean up and exit. */
  if ( !params ) {
    TRY( LALDestroyRandomParams( stat->statusPtr, &localParams ),
	 stat );
  }
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
