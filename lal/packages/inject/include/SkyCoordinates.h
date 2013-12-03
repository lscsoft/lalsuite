/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Teviet Creighton
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

#ifndef _SKYCOORDINATES_H
#define _SKYCOORDINATES_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup SkyCoordinates_h
 * @{
 */

/** \name Error codes */
/*@{*/
#define SKYCOORDINATESH_ENUL  1	/**< Unexpected null pointer in arguments */
#define SKYCOORDINATESH_ESYS  2	/**< Wrong coordinate system in input */
#define SKYCOORDINATESH_EZERO 3	/**< Angular coordinates undefined at origin */
#define SKYCOORDINATESH_ESING 4	/**< Point is inside singular ellipsoid */
/*@}*/

/** \cond DONT_DOXYGEN */
#define SKYCOORDINATESH_MSGENUL  "Unexpected null pointer in arguments"
#define SKYCOORDINATESH_MSGESYS  "Wrong coordinate system in input"
#define SKYCOORDINATESH_MSGEZERO "Angular coordinates undefined at origin"
#define SKYCOORDINATESH_MSGESING "Point is inside singular ellipsoid"
/** \endcond */


/*---------- exported types ---------- */

/**
 * This enumerated type is used to identify data as being in one of the
 * coordinate systems discussed in \ref SkyCoordinates_h.
 */
typedef enum {
  COORDINATESYSTEM_HORIZON,	/**< A horizon coordinate system. */
  COORDINATESYSTEM_GEOGRAPHIC,	/**< The Earth-fixed geographic coordinate system. */
  COORDINATESYSTEM_EQUATORIAL,	/**< The sky-fixed equatorial coordinate system. */
  COORDINATESYSTEM_ECLIPTIC,	/**< The ecliptic coordinate system. */
  COORDINATESYSTEM_GALACTIC	/**< The galactic coordinate system. */
} CoordinateSystem;

/**
 * This structure stores the two spherical coordinates of a sky position;
 * ie a generic latitude and longitude; the structure is not defined
 * specific to a particular coordinate system, but maintains a tag
 * indicating which coordinate system it is expressed in.
 */
typedef struct tagSkyPosition {
  REAL8 longitude;		/**< The longitudinal coordinate (in radians), as defined above.*/
  REAL8 latitude;		/**< The latitudinal coordinate (in radians), as defined above. */
  CoordinateSystem system; 	/**< The coordinate system in which latitude/longitude are expressed. */
} SkyPosition;

/**
 * This structure stores the location of a point on (or near) the surface
 * of the Earth in both geodetic and geocentric coordinates, as described
 * in TerrestrialCoordinates.c .
 */
typedef struct tagEarthPosition {
  SkyPosition geodetic; 	/**< The geographic coordinates of the
				 * upward vertical direction from the point; that is, the point's
				 * <em>geodetic</em> latitude and longitude. */

  REAL8 elevation;		/**< The vertical distance of the point above the reference ellipsoid,
				 * in metres.*/

  REAL8 x, y, z;		/**< The Earth-fixed geocentric Cartesian coordinates of the point,
				 *in metres.*/

  REAL8 radius;			/**< The distance of the point from the geocentre, in metres. */

  SkyPosition geocentric;	/**<  The geographic coordinates of the direction from the centre
				 * of the Earth through the point; that is, the point's
				 * <em>geocentric</em> latitude and longitude.*/
} EarthPosition;


/**
 * This structure stores parameters for the function <tt>LALConvertSkyPosition()</tt>.
 */
typedef struct tagConvertSkyParams {
  CoordinateSystem system;	/**<  The coordinate system to which one is transforming. */

  SkyPosition *zenith;		/**< The position of the zenith of the horizon coordinate system;
				 * may be <tt>NULL</tt> if one is neither converting to nor from
				 * a horizon system. */

  LIGOTimeGPS *gpsTime;		/**< The GPS time for conversions between Earth-fixed and
				 * sky-fixed coordinates; may be <tt>NULL</tt> if no such conversion
				 * is required (or if one is transforming to or from horizon
				 * coordinates and <tt>*zenith</tt> is given in the sky-fixed
				 * equatorial system). */
} ConvertSkyParams;

/*@}*/

/* ---------- Function prototypes ---------- */

void
LALGalacticToEquatorial( LALStatus   *,
			 SkyPosition *output,
			 SkyPosition *input );

void
LALEquatorialToGalactic( LALStatus   *,
			 SkyPosition *output,
			 SkyPosition *input );

void
LALEclipticToEquatorial( LALStatus   *,
			 SkyPosition *output,
			 SkyPosition *input );

void
LALEquatorialToEcliptic( LALStatus   *,
			 SkyPosition *output,
			 SkyPosition *input );

void
LALGeographicToEquatorial( LALStatus   *,
			   SkyPosition *output,
			   SkyPosition *input,
			   LIGOTimeGPS *gpsTime );

void
LALEquatorialToGeographic( LALStatus   *,
			   SkyPosition *output,
			   SkyPosition *input,
			   LIGOTimeGPS *gpsTime );

void
LALSystemToHorizon( LALStatus   *,
		    SkyPosition *output,
		    SkyPosition *input,
		    const SkyPosition *zenith );

void
LALHorizonToSystem( LALStatus   *,
		    SkyPosition *output,
		    SkyPosition *input,
		    const SkyPosition *zenith );

void
LALGeodeticToGeocentric( LALStatus *, EarthPosition *location );

void
LALGeocentricToGeodetic( LALStatus *, EarthPosition *location );

void
LALConvertSkyCoordinates( LALStatus        *,
			  SkyPosition      *output,
			  SkyPosition      *input,
			  ConvertSkyParams *params );

void LALNormalizeSkyPosition (LALStatus *status, SkyPosition *posOut, const SkyPosition *posIn);

#ifdef SWIG /* SWIG interface directives */
SWIGLAL(INOUT_SCALARS(double*, longitude, latitude));
#endif

void XLALNormalizeSkyPosition ( double *restrict longitude, double *restrict latitude );

#ifdef  __cplusplus
}
#endif

#endif /* _SKYCOORDINATES_H */
