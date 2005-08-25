/** \file
 *  \ingroup SkyCoordinates
 *  \author Creighton, T. D.
 *  \date 2002
 *  \brief This header covers routines to perform coordinate transformations
 *   among the various spherical coordinate systems used in astronomy.
 */

#ifndef _SKYCOORDINATES_H
#define _SKYCOORDINATES_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( SKYCOORDINATESH, "$Id$" );

/** \name Error codes */
/*@{*/
#define SKYCOORDINATESH_ENUL  1
#define SKYCOORDINATESH_ESYS  2
#define SKYCOORDINATESH_EZERO 3
#define SKYCOORDINATESH_ESING 4

#define SKYCOORDINATESH_MSGENUL  "Unexpected null pointer in arguments"
#define SKYCOORDINATESH_MSGESYS  "Wrong coordinate system in input"
#define SKYCOORDINATESH_MSGEZERO "Angular coordinates undefined at origin"
#define SKYCOORDINATESH_MSGESING "Point is inside singular ellipsoid"
/*@}*/

/*---------- exported types ---------- */

/** This enumerated type is used to identify data as being in one of the
 *  coordinate systems discussed in \ref SkyCoordinates.  */
typedef enum {
  COORDINATESYSTEM_HORIZON,	/**< A horizon coordinate system. */
  COORDINATESYSTEM_GEOGRAPHIC,	/**< The Earth-fixed geographic coordinate system. */
  COORDINATESYSTEM_EQUATORIAL,	/**< The sky-fixed equatorial coordinate system. */
  COORDINATESYSTEM_ECLIPTIC,	/**< The ecliptic coordinate system. */
  COORDINATESYSTEM_GALACTIC	/**< The galactic coordinate system. */
} CoordinateSystem;

/** This structure stores the two spherical coordinates of a sky position;
 * ie a generic latitude and longitude; the structure is not defined
 * specific to a particular coordinate system, but maintains a tag
 * indicating which coordinate system it is expressed in.
 */
typedef struct tagSkyPosition {
  REAL8 longitude;		/**< The longitudinal coordinate (in radians), as defined above.*/
  REAL8 latitude;		/**< The latitudinal coordinate (in radians), as defined above. */
  CoordinateSystem system; 	/**< The coordinate system in which latitude/longitude are expressed. */
} SkyPosition;

/** This structure stores the location of a point on (or near) the surface
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


/** This structure stores parameters for the function <tt>LALConvertSkyPosition()</tt>. 
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

#ifdef  __cplusplus
}
#endif

#endif /* _SKYCOORDINATES_H */
