/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Teviet Creighton, John Whelan
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
#include <lal/SkyCoordinates.h>

#define LAL_ALPHAGAL (3.366032942)
#define LAL_DELTAGAL (0.473477302)
#define LAL_LGAL     (0.576)

/**
 * \author Creighton, T. D.
 * \addtogroup SkyCoordinates_c
 * \brief Automatically converts among sky coordinate systems.
 *
 * The function <tt>LALConvertSkyCoordinates()</tt> transforms the contents
 * of <tt>*input</tt> to the system
 * specified in <tt>*params</tt>, storing the result in <tt>*output</tt>
 * (which may point to the same object as <tt>*input</tt> for an in-place
 * transformation).  The routine makes calls to the functions in
 * \ref CelestialCoordinates.c and \ref TerrestrialCoordinates.c as
 * required; the <tt>*params</tt> object must store any data fields
 * required by these functions, or an error will occur.
 *
 * The function <tt>LALNormalizeSkyPosition()</tt> "normalizes" any given
 * (spherical) sky-position (in radians), which means it projects the
 * angles into \f$[0, 2\pi) \times [-\pi/2, \pi/2]\f$ if they lie outside.
 *
 * \heading{Algorithm}
 *
 * <tt>LALConvertSkyCoordinates()</tt> is structured as a simple loop over
 * transformations, each
 * of which moves the output sky position one step closer to the desired
 * final coordinates system.  The usual "flow" of the algorithm is:
 *
 * \image html  SkyCoordinates_conversions.png
 * \image latex SkyCoordinates_conversions.eps
 *
 * although one can also convert directly between equatorial and horizon
 * coordinate systems if <tt>params->zenith</tt> is given in equatorial
 * coordinates (i.e.\ if its longitudinal coordinate is the local mean
 * sidereal time rather than the geographic longitude of the observer).
 * This leads to the only error checking done within this function: when
 * transforming to horizon coordinates, it checks that
 * <tt>params->zenith</tt> is either in sky-fixed equatorial or Earth-fixed
 * geographic coordinates.  Other than this, error checking is left to
 * the secondary function call; if a parameter is absent or poorly
 * formatted, the called function will return an error.
 *
 */
/*@{*/

/** \see See documentation in \ref SkyCoordinates_c */
void
LALConvertSkyCoordinates( LALStatus        *stat,
			  SkyPosition      *output,
			  SkyPosition      *input,
			  ConvertSkyParams *params )
{
  SkyPosition temp; /* temporary sky position (duh) */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( params, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Start looping! */
  temp = *input;
  while ( temp.system != params->system ) {

    if ( temp.system == COORDINATESYSTEM_HORIZON )
      /* Only one possible direction. */
      TRY( LALHorizonToSystem( stat->statusPtr, &temp, &temp,
			       params->zenith ), stat );

    else if ( temp.system == COORDINATESYSTEM_GEOGRAPHIC ) {
      /* Two possible directions.  But, if headed towards the horizon
         coordinate system, make sure we can get there! */
      if ( params->system == COORDINATESYSTEM_HORIZON ) {
	if ( params->zenith ) {
	  if ( params->zenith->system == COORDINATESYSTEM_GEOGRAPHIC )
	    TRY( LALSystemToHorizon( stat->statusPtr, &temp, &temp,
				     params->zenith ), stat );
	  else if ( params->zenith->system == COORDINATESYSTEM_EQUATORIAL )
	    TRY( LALGeographicToEquatorial( stat->statusPtr, &temp, &temp,
					    params->gpsTime ), stat );
	  else
	    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
	} else
	  ABORT( stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
      } else
	TRY( LALGeographicToEquatorial( stat->statusPtr, &temp, &temp,
					params->gpsTime ), stat );
    }

    else if ( temp.system == COORDINATESYSTEM_EQUATORIAL ) {
      /* Up to four possible directions, depending on
         params->zenith->system. */
      if ( ( params->system == COORDINATESYSTEM_HORIZON ) &&
	   ( params->zenith ) &&
	   ( params->zenith->system == COORDINATESYSTEM_EQUATORIAL ) )
	TRY( LALSystemToHorizon( stat->statusPtr, &temp, &temp,
				 params->zenith ), stat );
      else if ( params->system == COORDINATESYSTEM_ECLIPTIC )
	TRY( LALEquatorialToEcliptic( stat->statusPtr, &temp, &temp ),
	     stat );
      else if ( params->system == COORDINATESYSTEM_GALACTIC )
	TRY( LALEquatorialToGalactic( stat->statusPtr, &temp, &temp ),
	     stat );
      else
	TRY( LALEquatorialToGeographic( stat->statusPtr, &temp, &temp,
					params->gpsTime ), stat );
    }

    else if ( temp.system == COORDINATESYSTEM_ECLIPTIC ) {
      /* Only one possible direction. */
      TRY( LALEclipticToEquatorial( stat->statusPtr, &temp, &temp ),
	   stat );
    }

    else if ( temp.system == COORDINATESYSTEM_GALACTIC ) {
      /* Only one possible direction. */
      TRY( LALGalacticToEquatorial( stat->statusPtr, &temp, &temp ),
	   stat );
    }

    else
      ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* We've gotten to the correct coordinate system.  Set the output
     and return. */
  *output = temp;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );

} /* LALConvertSkyCoordinates */



/**
 * \deprecated Use XLALNormalizeSkyPosition() instead.
 */
void
LALNormalizeSkyPosition (LALStatus *stat,		/**< pointer to LALStatus structure */
			 SkyPosition *posOut, 		/**< [out] normalized sky-position */
			 const SkyPosition *posIn)	/**< [in] general sky-position */
{
  SkyPosition tmp;	/* allow posOut == posIn */

  INITSTATUS(stat);

  ASSERT (posIn, stat, SKYCOORDINATESH_ENUL ,  SKYCOORDINATESH_MSGENUL );
  ASSERT (posOut, stat, SKYCOORDINATESH_ENUL ,  SKYCOORDINATESH_MSGENUL );

  tmp = *posIn;

  XLALNormalizeSkyPosition ( &tmp.longitude, &tmp.latitude );

  *posOut = tmp;

  RETURN (stat);

} /* LALNormalizeSkyPosition() */


/**
 * If sky-position is not in the canonical range
 * \f$(\alpha,\delta)\in [0,2\pi]\times[-\pi/2, \pi/2]\f$, normalize it
 * by mapping it into this coordinate-interval.
 * Based on Alicia's function with some additional "unwinding" added.
 */
void
XLALNormalizeSkyPosition ( double *restrict longitude,   /**< [in,out] sky-position longitude to normalize*/
                           double *restrict latitude     /**< [in,out] sky-position latitude to normalize*/
                           )
{

  /* FIRST STEP: completely "unwind" positions, i.e. make sure that
   * [0 <= alpha < 2pi] and [-pi < delta <= pi] */
  /* normalize longitude */
  *longitude -= floor(*longitude / LAL_TWOPI) * LAL_TWOPI;

  /* pre-normalize (unwind) latitude */
  *latitude += LAL_PI;
  *latitude -= floor(*latitude / LAL_TWOPI) * LAL_TWOPI;
  *latitude -= LAL_PI;

  /* SECOND STEP: get latitude into canonical interval [-pi/2 <= delta <= pi/2 ] */
  /* this requires also a change in longitude by adding/subtracting PI */
  if (*latitude > LAL_PI_2)
    {
      *latitude = LAL_PI - *latitude;
      if (*longitude < LAL_PI)
	{
	  *longitude += LAL_PI;
	}
      else
	{
	  *longitude -= LAL_PI;
	}
    }

  if (*latitude < -LAL_PI_2)
    {
      *latitude = -LAL_PI - *latitude;
      if (*longitude < LAL_PI)
	{
	  *longitude += LAL_PI;
	}
      else
	{
	  *longitude -= LAL_PI;
	}
    }

  return;

} /* XLALNormalizeSkyPosition() */

/*@}*/
