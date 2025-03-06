/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef _LALDETECTORS_H
#define _LALDETECTORS_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALDetectors_h Header LALDetectors.h
 * \ingroup lal_tools
 * \author J. T. Whelan and J. D. E. Creighton
 *
 * \brief This header defines structures to hold the basic data describing a gravitational wave detector.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LALDetectors.h>
 * \endcode
 *
 * According to the common frame format specification
 * \cite LIGOVIRGO_2000 the geometry of an interferometric
 * detector will be stored in a \c FrDetector structure, specifying
 * the location of the detector vertex and the orientation of its arms in
 * geodetic coördinates suited to geographical surveying.  Resonant
 * bars and other sorts of detectors, if they write their data to frames,
 * are expected to fill this structure with their location and
 * orientation in some way suited to the detector type.
 *
 * For most data analysis tasks, however, any gravitational wave detector
 * can be described by its location in an Earth-fixed rotating reference
 * frame, as well as a <em>response tensor</em> \f$d^{ab}\f$, constant in the
 * same frame, which defines the "strain" \f$h\f$ measured by the detector in
 * terms of the metric perturbation \f$h_{ab}\f$ as
 * \f{equation}{
 * h = h_{ab} \, d^{ab}
 * \ .
 * \f}
 *
 * This header defines a \c LALFrDetector structure which contains
 * essentially the same information as the \c FrDetector structure,
 * as well as a \c LALDetector structure which contains the
 * Cartesian coördinates of the detector along with the components of
 * the response tensor \f$d^{ab}\f$ in the same coördinate system.
 *
 * ### The Geodetic Coördinate System ###
 *
 * Geodetic coördinates are spheroidal coördinates
 * based on the WGS-84 Earth Model, which is an
 * oblate spheroid with equatorial radius \f$a=6.378137\times
 * 10^6\,\textrm{m}\f$ and polar radius \f$b=6.356752314\times
 * 10^6\,\textrm{m}\f$.  Any point in space can be located according to its
 * longitude, latitude, and elevation.  The \e longitude \f$\lambda\f$
 * is the angle between the half-plane bounded by the symmetry axis of
 * the reference ellipsoid containing the point in question and the
 * half-plane plane containing the Prime Meridian; it is measured in
 * radians, increases to the East, and ranges from
 * \f$-\pi\f$ to \f$\pi\f$.  The \e latitude \f$\beta\f$ is the
 * angle between the ray which is normal to the ellipsoid and passes
 * through the point in question and the equatorial plane; it is measured
 * in radians, increases to the North, and ranges
 * from \f$-\pi/2\f$ to \f$\pi/2\f$.  The \e elevation \f$h\f$ is the
 * signed distance along this ray from the reference ellipsoid to the
 * point in question.  This coördinate system is described in more
 * detail in \cite Althouse_1999 .
 *
 * ### Altitude and Azimuth Angles ###
 *
 * The \c LALFrDetector structure stores the directions along the
 * two arms of an interferometer in an altitude/azimuth representation
 * with respect to the local tangent plane to the reference ellipsoid,
 * known as the local horizontal.  The altitude \f${\mathcal{A}}\f$ is the angle the
 * direction vector makes with the horizontal, \f${\mathcal{A}} > 0\f$ meaning above
 * horizontal, \f${\mathcal{A}} < 0\f$ below.  The azimuth angle \f$\zeta\f$ is found by
 * projecting the direction onto the local horizontal plane, then
 * measuring the angle clockwise from North to this projected direction.
 *
 * ### The Cartesian Coördinate System ###
 *
 * The position vector and response tensor contained in the
 * \c LALDetector structure are defined in
 * a simple orthonormal coördinate system with its origin at
 * the center of the earth, an \f$x^1\f$ axis which pierces the Earth's
 * surface at the intersection of the equator and the prime meridian, an
 * \f$x^2\f$ axis which pierces the earth's surface at \f$\pi/2\f$ radians East
 * longitude on the equator, and an \f$x^3\f$ axis which pierces the Earth's
 * surface at the North Pole.  The coördinates \f$x^1\f$, \f$x^2\f$, \f$x^3\f$
 * correspond to the Earth-fixed coördinates \f$X_E\f$, \f$Y_E\f$, \f$Z_E\f$
 * defined in \cite Althouse_1999 , respectively.
 *
 * The relationship between geodetic and Cartesian coördinates is
 * given by
 * \f{align}{
 * \label{tools_e_cart1}
 * x^1 &=\left(
 * \frac{a^2}{\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}}
 * + h
 * \right) \cos\beta\cos\lambda \\
 * \label{tools_e_cart2}
 * x^2 &=\left(
 * \frac{a^2}{\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}}
 * + h
 * \right) \cos\beta\sin\lambda \\
 * \label{tools_e_cart3}
 * x^3 &=\left(
 * \frac{b^2}{\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}}
 * + h
 * \right) \sin\beta \\
 * \f}
 *
 * ### Cached Detectors ###
 *
 * In practice, we will often be
 * working with fixed unchanging site geometry, e.g., for the LIGO
 * interferometers; to avoid constantly reconstructing the corresponding
 * \c LALDetectors, we should define some constant
 * \c LALDetectors describing them.  Those are stored in a constant
 * array of \c LALDetector structures known as
 * \c lalCachedDetectors, which is declared \c extern in this
 * header and defined in \ref CreateDetector_c.
 *
 * The <tt>LALCreateDetector()</tt> routine will first look through the
 * \c lalCachedDetectors array for a \c LALDetector structure
 * with matching \c type and <tt>frDetector.name</tt> fields; if it
 * finds one, it returns a copy of that; if not, it creates one.
 *
 * For example, the \c LALDetector representing LIGO Hanford 4km (H1) in
 * differential mode is <tt>lalCachedDetectors[LAL_LHO_4K_DETECTOR]</tt>.
 *
 */
/** @{ */

/** \name Error Codes */
/** @{ */
#define LALDETECTORSH_ENULLP        1	/**< Null pointer */
#define LALDETECTORSH_ETYPE         2	/**< Unsupported detector type */
/** @} */

/** \cond DONT_DOXYGEN */
#define LALDETECTORSH_MSGENULLP     "Null pointer"
#define LALDETECTORSH_MSGETYPE      "Unsupported detector type"

#define LALDETECTORSH_PRINTF        0
/** \endcond */

/** Enumeration of Detectors: follows order of DQ bit assignments */
typedef enum tagLALDetectorEnum {
	LAL_TAMA_300_DETECTOR	=	0,
	LAL_VIRGO_CITF_DETECTOR	=	1,
	LAL_VIRGO_DETECTOR	=	2,
	LAL_GEO_600_DETECTOR	=	3,
	LAL_LHO_2K_DETECTOR	=	4,
	LAL_LHO_4K_DETECTOR	=	5,
	LAL_LLO_4K_DETECTOR	=	6,
	LAL_CIT_40_DETECTOR	=	7,
	LAL_ALLEGRO_DETECTOR	=	8,
	LAL_AURIGA_DETECTOR	=	9,
	LAL_EXPLORER_DETECTOR	=	10,
	LAL_NIOBE_DETECTOR	=	11,
	LAL_NAUTILUS_DETECTOR	=	12,
	LAL_ACIGA_DETECTOR	=	13,
	LAL_KAGRA_DETECTOR	=	14,
	LAL_LIO_4K_DETECTOR	=	15,
	LAL_ET1_DETECTOR	=	16,
	LAL_ET2_DETECTOR	=	17,
	LAL_ET3_DETECTOR	=	18,
	LAL_ET0_DETECTOR	=	19,
	LAL_NUM_DETECTORS	=	20
}
LALDetectorEnum;

/** \name Detector DQ bit assignments (2 bits per detector) */
/** @{ */
#define LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(d)	(LAL_INT8_C(1) << 2 * (d))

#define LAL_TAMA_300_DETECTOR_BIT	LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_TAMA_300_DETECTOR)
#define LAL_VIRGO_CITF_DETECTOR_BIT	LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_VIRGO_CITF_DETECTOR)
#define LAL_VIRGO_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_VIRGO_DETECTOR)
#define LAL_GEO_600_DETECTOR_BIT	LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_GEO_600_DETECTOR)
#define LAL_LHO_2K_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_LHO_2K_DETECTOR)
#define LAL_LHO_4K_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_LHO_4K_DETECTOR)
#define LAL_LLO_4K_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_LLO_4K_DETECTOR)
#define LAL_CIT_40_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_CIT_40_DETECTOR)
#define LAL_ALLEGRO_DETECTOR_BIT	LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_ALLEGRO_DETECTOR)
#define LAL_AURIGA_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_AURIGA_DETECTOR)
#define LAL_NIOBE_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_NIOBE_DETECTOR)
#define LAL_NAUTILUS_DETECTOR_BIT	LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_NAUTILUS_DETECTOR)
#define LAL_ACIGA_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_ACIGA_DETECTOR)
#define LAL_KAGRA_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_KAGRA_DETECTOR)
#define LAL_LIO_4K_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_LIO_4K_DETECTOR)
#define LAL_ET1_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_ET1_DETECTOR)
#define LAL_ET2_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_ET2_DETECTOR)
#define LAL_ET3_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_ET3_DETECTOR)
#define LAL_ET0_DETECTOR_BIT		LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(LAL_ET0_DETECTOR)
/** @} */

#ifdef SWIG /* SWIG interface only; SWIG cannot parse the _DETECTOR_BIT constants */
%inline %{
INT8 XLALDetectorDQBitFromDetectorEnum(const LALDetectorEnum d)
{
	return LAL_DETECTOR_DQ_BIT_FROM_DETECTOR_ENUM(d);
}
%}
#endif /* SWIG */


/**
 * Detector type, which determines how the detector response is determined.
 * Since data from bars as well as interferometers can be written to
 * frames, we need an additional piece of information to interpret the
 * site geometry data specified in the \c LALFrDetector
 * structure; for instance, is the x arm really the x arm or is it the
 * long axis of a bar?  The \c LALDetectorType enumeration
 * provides a way to keep track of that.
 *
 */
typedef enum tagLALDetectorType {
	LALDETECTORTYPE_ABSENT,		/**< No FrDetector associated with this detector */
	LALDETECTORTYPE_IFODIFF,	/**< IFO in differential mode */
	LALDETECTORTYPE_IFOXARM,	/**< IFO in one-armed mode (X arm) */
	LALDETECTORTYPE_IFOYARM,	/**< IFO in one-armed mode (Y arm) */
	LALDETECTORTYPE_IFOCOMM,	/**< IFO in common mode */
	LALDETECTORTYPE_CYLBAR		/**< Cylindrical bar */
}
LALDetectorType;


/**
 * Detector frame data structure
 * Structure to contain the data that appears in a FrDetector structure
 * in frame data.
 */
typedef struct tagLALFrDetector
{
	CHAR	name[LALNameLength];	/**< A unique identifying string */
	CHAR	prefix[3];		/**< Two-letter prefix for detector's channel names */
	REAL8	vertexLongitudeRadians;	/**< The geodetic longitude \f$\lambda\f$ of the vertex in radians */
	REAL8	vertexLatitudeRadians;	/**< The geodetic latitude \f$\beta\f$ of the vertex in radians */
	REAL4	vertexElevation;	/**< The height of the vertex above the reference ellipsoid in meters */
	REAL4	xArmAltitudeRadians;	/**< The angle \f${\mathcal{A}}_X\f$ up from the local tangent plane of the reference ellipsoid to the X arm (or bar's cylidrical axis) in radians */
	REAL4	xArmAzimuthRadians;	/**< The angle \f$\zeta_X\f$ clockwise from North to the projection of the X arm (or bar's cylidrical axis) into the local tangent plane of the reference ellipsoid in radians */
	REAL4	yArmAltitudeRadians;	/**< The angle \f${\mathcal{A}}_Y\f$ up from the local tangent plane of the reference ellipsoid to the Y arm in radians (unused for bars: set it to zero) */
	REAL4	yArmAzimuthRadians;	/**< The angle \f$\zeta_Y\f$ clockwise from North to the projection of the Y arm into the local tangent plane of the reference ellipsoid in radians (unused for bars: set it to zero) */
	REAL4	xArmMidpoint;		/**< The distance to the midpoint of the X arm in meters (unused for bars: set it to zero) */
	REAL4	yArmMidpoint;		/**< The distance to the midpoint of the Y arm in meters (unused for bars: set it to zero) */
}
LALFrDetector;


/**
 * Detector structure
 *
 * Structure to contain detector data in the format most easily used
 * by the LAL routines.
 */
typedef struct tagLALDetector
{
	REAL8		location[3];	/**< The three components, in an Earth-fixed Cartesian coordinate system, of the position vector from the center of the Earth to the detector in meters */
	REAL4		response[3][3];	/**< The Earth-fixed Cartesian components of the detector's response tensor \f$d^{ab}\f$ */
	LALDetectorType	type;		/**< The type of the detector (e.g., IFO in differential mode, cylindrical bar, etc.) */
	LALFrDetector	frDetector;	/**< The original LALFrDetector structure from which this was created */
}
LALDetector;


/** Pre-existing detectors. */
extern const LALDetector lalCachedDetectors[LAL_NUM_DETECTORS];

/** @} */


/* Routine to create a LALDetector. */
LALDetector * XLALCreateDetector( LALDetector *detector, const LALFrDetector *frDetector, LALDetectorType type );



/* Interferometric Detectors */

/**
 * \defgroup DetectorConstants Detector Constants
 * \ingroup LALDetectors_h
 * \brief Constants describing various gravitational wave detectors
 *
 * The \ref LALDetectors_h also defines numerical constants that describe the location and
 * geometry of several operating gravitational wave detectors.
 * These detectors are both resonant mass (bar) detectors and interferometric
 * detectors.
 * <ul>
 * <li> Data for the resonant mass detectors is taken from:
 * http://igec.lnl.infn.it/cgi-bin/browser.pl?Level=0,3,1
 * and \cite FinnLazzarini_2001
 *
 * <li> Data for LIGO detectors is taken from \cite Althouse_1999
 *
 * <li> Data for the VIRGO detector is provided by Benoit Mours.
 *
 * <li> Data for the GEO detector is taken from:
 *
 * http://www.geo600.uni-hannover.de/geo600/project/location.html
 *
 * <li> Data for the TAMA detector is provided by Masa-Katsu Fujimoto
 *
 * <li> Data for the Caltech detector is taken from \cite Allen_1996
 *
 * <li> Data for the KAGRA detector was provided by Yousuke Itoh
 * and is derived from data from
 *
 * > Yoshio Saito, "KAGRA location", KAGRA Technical Document JGW-G1503824
 * > http://gwdoc.icrr.u-tokyo.ac.jp/cgi-bin/DocDB/ShowDocument?docid=3824
 *
 * </ul>
 *
 * See the technical document \cite ABCCRW_2001 for details.
 *
 * Data in this file (e.g., angle conventions etc.) is intended
 * to conform to the conventions of the Frame format specification \cite LIGOVIRGO_2000
 *
 */
/** @{ */

/**
 * \name TAMA 300m Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * TAMA 300m Interferometric Detector.
 */
/** @{ */
#define LAL_TAMA_300_DETECTOR_NAME               	"TAMA_300"	/**< TAMA_300 detector name string */
#define LAL_TAMA_300_DETECTOR_PREFIX             	"T1"	/**< TAMA_300 detector prefix string */
#define LAL_TAMA_300_DETECTOR_LONGITUDE_RAD      	2.43536359469	/**< TAMA_300 vertex longitude (rad) */
#define LAL_TAMA_300_DETECTOR_LATITUDE_RAD       	0.62267336022	/**< TAMA_300 vertex latitude (rad) */
#define LAL_TAMA_300_DETECTOR_ELEVATION_SI       	90	/**< TAMA_300 vertex elevation (m) */
#define LAL_TAMA_300_DETECTOR_ARM_X_AZIMUTH_RAD  	4.71238898038	/**< TAMA_300 x arm azimuth (rad) */
#define LAL_TAMA_300_DETECTOR_ARM_Y_AZIMUTH_RAD  	3.14159265359	/**< TAMA_300 y arm azimuth (rad) */
#define LAL_TAMA_300_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< TAMA_300 x arm altitude (rad) */
#define LAL_TAMA_300_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< TAMA_300 y arm altitude (rad) */
#define LAL_TAMA_300_DETECTOR_ARM_X_MIDPOINT_SI  	150.00000000000	/**< TAMA_300 x arm midpoint (m) */
#define LAL_TAMA_300_DETECTOR_ARM_Y_MIDPOINT_SI  	150.00000000000	/**< TAMA_300 y arm midpoint (m) */
#define LAL_TAMA_300_VERTEX_LOCATION_X_SI        	-3.94640899111e+06	/**< TAMA_300 x-component of vertex location in Earth-centered frame (m) */
#define LAL_TAMA_300_VERTEX_LOCATION_Y_SI        	3.36625902802e+06	/**< TAMA_300 y-component of vertex location in Earth-centered frame (m) */
#define LAL_TAMA_300_VERTEX_LOCATION_Z_SI        	3.69915069233e+06	/**< TAMA_300 z-component of vertex location in Earth-centered frame (m) */
#define LAL_TAMA_300_ARM_X_DIRECTION_X           	0.64896940530	/**< TAMA_300 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_X_DIRECTION_Y           	0.76081450498	/**< TAMA_300 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_X_DIRECTION_Z           	-0.00000000000	/**< TAMA_300 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_Y_DIRECTION_X           	-0.44371376921	/**< TAMA_300 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_Y_DIRECTION_Y           	0.37848471479	/**< TAMA_300 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_TAMA_300_ARM_Y_DIRECTION_Z           	-0.81232223390	/**< TAMA_300 z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */

/**
 * \name VIRGO_CITF Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * VIRGO_CITF Interferometric Detector.  FIXME: the armlength is a stub.
 * FIXME: this IFO is called V0 rather than to avoid name clash.
 */
/** @{ */
#define LAL_VIRGO_CITF_DETECTOR_NAME               	"VIRGO_CITF"	/**< VIRGO_CITF detector name string */
#define LAL_VIRGO_CITF_DETECTOR_PREFIX             	"V0"	/**< VIRGO_CITF detector prefix string */
#define LAL_VIRGO_CITF_DETECTOR_LONGITUDE_RAD      	0.18333805213	/**< VIRGO_CITF vertex longitude (rad) */
#define LAL_VIRGO_CITF_DETECTOR_LATITUDE_RAD       	0.76151183984	/**< VIRGO_CITF vertex latitude (rad) */
#define LAL_VIRGO_CITF_DETECTOR_ELEVATION_SI       	51.884	/**< VIRGO_CITF vertex elevation (m) */
#define LAL_VIRGO_CITF_DETECTOR_ARM_X_AZIMUTH_RAD  	0.33916285222	/**< VIRGO_CITF x arm azimuth (rad) */
#define LAL_VIRGO_CITF_DETECTOR_ARM_Y_AZIMUTH_RAD  	5.05155183261	/**< VIRGO_CITF y arm azimuth (rad) */
#define LAL_VIRGO_CITF_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< VIRGO_CITF x arm altitude (rad) */
#define LAL_VIRGO_CITF_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< VIRGO_CITF y arm altitude (rad) */
#define LAL_VIRGO_CITF_DETECTOR_ARM_X_MIDPOINT_SI  	0.00000000000	/**< VIRGO_CITF x arm midpoint (m) */
#define LAL_VIRGO_CITF_DETECTOR_ARM_Y_MIDPOINT_SI  	0.00000000000	/**< VIRGO_CITF y arm midpoint (m) */
#define LAL_VIRGO_CITF_VERTEX_LOCATION_X_SI        	4.54637409900e+06	/**< VIRGO_CITF x-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_CITF_VERTEX_LOCATION_Y_SI        	8.42989697626e+05	/**< VIRGO_CITF y-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_CITF_VERTEX_LOCATION_Z_SI        	4.37857696241e+06	/**< VIRGO_CITF z-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_CITF_ARM_X_DIRECTION_X           	-0.70045821479	/**< VIRGO_CITF x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_CITF_ARM_X_DIRECTION_Y           	0.20848948619	/**< VIRGO_CITF y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_CITF_ARM_X_DIRECTION_Z           	0.68256166277	/**< VIRGO_CITF z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_CITF_ARM_Y_DIRECTION_X           	-0.05379255368	/**< VIRGO_CITF x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_VIRGO_CITF_ARM_Y_DIRECTION_Y           	-0.96908180549	/**< VIRGO_CITF y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_VIRGO_CITF_ARM_Y_DIRECTION_Z           	0.24080451708	/**< VIRGO_CITF z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */

/**
 * \name VIRGO 3km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * VIRGO 3km Interferometric Detector.
 */
/** @{ */
#define LAL_VIRGO_DETECTOR_NAME               	"VIRGO"	/**< VIRGO detector name string */
#define LAL_VIRGO_DETECTOR_PREFIX             	"V1"	/**< VIRGO detector prefix string */
#define LAL_VIRGO_DETECTOR_LONGITUDE_RAD      	0.18333805213	/**< VIRGO vertex longitude (rad) */
#define LAL_VIRGO_DETECTOR_LATITUDE_RAD       	0.76151183984	/**< VIRGO vertex latitude (rad) */
#define LAL_VIRGO_DETECTOR_ELEVATION_SI       	51.884	/**< VIRGO vertex elevation (m) */
#define LAL_VIRGO_DETECTOR_ARM_X_AZIMUTH_RAD  	0.33916285222	/**< VIRGO x arm azimuth (rad) */
#define LAL_VIRGO_DETECTOR_ARM_Y_AZIMUTH_RAD  	5.05155183261	/**< VIRGO y arm azimuth (rad) */
#define LAL_VIRGO_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< VIRGO x arm altitude (rad) */
#define LAL_VIRGO_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< VIRGO y arm altitude (rad) */
#define LAL_VIRGO_DETECTOR_ARM_X_MIDPOINT_SI  	1500.00000000000	/**< VIRGO x arm midpoint (m) */
#define LAL_VIRGO_DETECTOR_ARM_Y_MIDPOINT_SI  	1500.00000000000	/**< VIRGO y arm midpoint (m) */
#define LAL_VIRGO_VERTEX_LOCATION_X_SI        	4.54637409900e+06	/**< VIRGO x-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_VERTEX_LOCATION_Y_SI        	8.42989697626e+05	/**< VIRGO y-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_VERTEX_LOCATION_Z_SI        	4.37857696241e+06	/**< VIRGO z-component of vertex location in Earth-centered frame (m) */
#define LAL_VIRGO_ARM_X_DIRECTION_X           	-0.70045821479	/**< VIRGO x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_ARM_X_DIRECTION_Y           	0.20848948619	/**< VIRGO y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_ARM_X_DIRECTION_Z           	0.68256166277	/**< VIRGO z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_VIRGO_ARM_Y_DIRECTION_X           	-0.05379255368	/**< VIRGO x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_VIRGO_ARM_Y_DIRECTION_Y           	-0.96908180549	/**< VIRGO y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_VIRGO_ARM_Y_DIRECTION_Z           	0.24080451708	/**< VIRGO z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */


/**
 * \name GEO 600m Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * GEO 600m Interferometric Detector.
 */
/** @{ */
#define LAL_GEO_600_DETECTOR_NAME               	"GEO_600"	/**< GEO_600 detector name string */
#define LAL_GEO_600_DETECTOR_PREFIX             	"G1"	/**< GEO_600 detector prefix string */
#define LAL_GEO_600_DETECTOR_LONGITUDE_RAD      	0.17116780435	/**< GEO_600 vertex longitude (rad) */
#define LAL_GEO_600_DETECTOR_LATITUDE_RAD       	0.91184982752	/**< GEO_600 vertex latitude (rad) */
#define LAL_GEO_600_DETECTOR_ELEVATION_SI       	114.425	/**< GEO_600 vertex elevation (m) */
#define LAL_GEO_600_DETECTOR_ARM_X_AZIMUTH_RAD  	1.19360100484	/**< GEO_600 x arm azimuth (rad) */
#define LAL_GEO_600_DETECTOR_ARM_Y_AZIMUTH_RAD  	5.83039279401	/**< GEO_600 y arm azimuth (rad) */
#define LAL_GEO_600_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< GEO_600 x arm altitude (rad) */
#define LAL_GEO_600_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< GEO_600 y arm altitude (rad) */
#define LAL_GEO_600_DETECTOR_ARM_X_MIDPOINT_SI  	300.00000000000	/**< GEO_600 x arm midpoint (m) */
#define LAL_GEO_600_DETECTOR_ARM_Y_MIDPOINT_SI  	300.00000000000	/**< GEO_600 y arm midpoint (m) */
#define LAL_GEO_600_VERTEX_LOCATION_X_SI        	3.85630994926e+06	/**< GEO_600 x-component of vertex location in Earth-centered frame (m) */
#define LAL_GEO_600_VERTEX_LOCATION_Y_SI        	6.66598956317e+05	/**< GEO_600 y-component of vertex location in Earth-centered frame (m) */
#define LAL_GEO_600_VERTEX_LOCATION_Z_SI        	5.01964141725e+06	/**< GEO_600 z-component of vertex location in Earth-centered frame (m) */
#define LAL_GEO_600_ARM_X_DIRECTION_X           	-0.44530676905	/**< GEO_600 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_GEO_600_ARM_X_DIRECTION_Y           	0.86651354130	/**< GEO_600 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_GEO_600_ARM_X_DIRECTION_Z           	0.22551311312	/**< GEO_600 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_GEO_600_ARM_Y_DIRECTION_X           	-0.62605756776	/**< GEO_600 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_GEO_600_ARM_Y_DIRECTION_Y           	-0.55218609524	/**< GEO_600 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_GEO_600_ARM_Y_DIRECTION_Z           	0.55058372486	/**< GEO_600 z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */


/**
 * \name LIGO Hanford Observatory 2km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Hanford Observatory 2km Interferometric Detector.
 */
/** @{ */
#define LAL_LHO_2K_DETECTOR_NAME               	"LHO_2k"	/**< LHO_2k detector name string */
#define LAL_LHO_2K_DETECTOR_PREFIX             	"H2"	/**< LHO_2k detector prefix string */
#define LAL_LHO_2K_DETECTOR_LONGITUDE_RAD      	-2.08405676917	/**< LHO_2k vertex longitude (rad) */
#define LAL_LHO_2K_DETECTOR_LATITUDE_RAD       	0.81079526383	/**< LHO_2k vertex latitude (rad) */
#define LAL_LHO_2K_DETECTOR_ELEVATION_SI       	142.554	/**< LHO_2k vertex elevation (m) */
#define LAL_LHO_2K_DETECTOR_ARM_X_AZIMUTH_RAD  	5.65487724844	/**< LHO_2k x arm azimuth (rad) */
#define LAL_LHO_2K_DETECTOR_ARM_Y_AZIMUTH_RAD  	4.08408092164	/**< LHO_2k y arm azimuth (rad) */
#define LAL_LHO_2K_DETECTOR_ARM_X_ALTITUDE_RAD 	-0.00061950000	/**< LHO_2k x arm altitude (rad) */
#define LAL_LHO_2K_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00001250000	/**< LHO_2k y arm altitude (rad) */
#define LAL_LHO_2K_DETECTOR_ARM_X_MIDPOINT_SI  	1004.50000000000	/**< LHO_2k x arm midpoint (m) */
#define LAL_LHO_2K_DETECTOR_ARM_Y_MIDPOINT_SI  	1004.50000000000	/**< LHO_2k y arm midpoint (m) */
#define LAL_LHO_2K_VERTEX_LOCATION_X_SI        	-2.16141492636e+06	/**< LHO_2k x-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_2K_VERTEX_LOCATION_Y_SI        	-3.83469517889e+06	/**< LHO_2k y-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_2K_VERTEX_LOCATION_Z_SI        	4.60035022664e+06	/**< LHO_2k z-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_2K_ARM_X_DIRECTION_X           	-0.22389266154	/**< LHO_2k x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_X_DIRECTION_Y           	0.79983062746	/**< LHO_2k y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_X_DIRECTION_Z           	0.55690487831	/**< LHO_2k z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_Y_DIRECTION_X           	-0.91397818574	/**< LHO_2k x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_Y_DIRECTION_Y           	0.02609403989	/**< LHO_2k y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_2K_ARM_Y_DIRECTION_Z           	-0.40492342125	/**< LHO_2k z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */


/**
 * \name LIGO Hanford Observatory 4km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Hanford Observatory 4km Interferometric Detector.
 */
/** @{ */
#define LAL_LHO_4K_DETECTOR_NAME               	"LHO_4k"	/**< LHO_4k detector name string */
#define LAL_LHO_4K_DETECTOR_PREFIX             	"H1"	/**< LHO_4k detector prefix string */
#define LAL_LHO_4K_DETECTOR_LONGITUDE_RAD      	-2.08405676917	/**< LHO_4k vertex longitude (rad) */
#define LAL_LHO_4K_DETECTOR_LATITUDE_RAD       	0.81079526383	/**< LHO_4k vertex latitude (rad) */
#define LAL_LHO_4K_DETECTOR_ELEVATION_SI       	142.554	/**< LHO_4k vertex elevation (m) */
#define LAL_LHO_4K_DETECTOR_ARM_X_AZIMUTH_RAD  	5.65487724844	/**< LHO_4k x arm azimuth (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_Y_AZIMUTH_RAD  	4.08408092164	/**< LHO_4k y arm azimuth (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_X_ALTITUDE_RAD 	-0.00061950000	/**< LHO_4k x arm altitude (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00001250000	/**< LHO_4k y arm altitude (rad) */
#define LAL_LHO_4K_DETECTOR_ARM_X_MIDPOINT_SI  	1997.54200000000	/**< LHO_4k x arm midpoint (m) */
#define LAL_LHO_4K_DETECTOR_ARM_Y_MIDPOINT_SI  	1997.52200000000	/**< LHO_4k y arm midpoint (m) */
#define LAL_LHO_4K_VERTEX_LOCATION_X_SI        	-2.16141492636e+06	/**< LHO_4k x-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_4K_VERTEX_LOCATION_Y_SI        	-3.83469517889e+06	/**< LHO_4k y-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_4K_VERTEX_LOCATION_Z_SI        	4.60035022664e+06	/**< LHO_4k z-component of vertex location in Earth-centered frame (m) */
#define LAL_LHO_4K_ARM_X_DIRECTION_X           	-0.22389266154	/**< LHO_4k x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_X_DIRECTION_Y           	0.79983062746	/**< LHO_4k y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_X_DIRECTION_Z           	0.55690487831	/**< LHO_4k z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_Y_DIRECTION_X           	-0.91397818574	/**< LHO_4k x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_Y_DIRECTION_Y           	0.02609403989	/**< LHO_4k y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LHO_4K_ARM_Y_DIRECTION_Z           	-0.40492342125	/**< LHO_4k z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */


/**
 * \name LIGO Livingston Observatory 4km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Livingston Observatory 4km Interferometric Detector.
 */
/** @{ */
#define LAL_LLO_4K_DETECTOR_NAME               	"LLO_4k"	/**< LLO_4k detector name string */
#define LAL_LLO_4K_DETECTOR_PREFIX             	"L1"	/**< LLO_4k detector prefix string */
#define LAL_LLO_4K_DETECTOR_LONGITUDE_RAD      	-1.58430937078	/**< LLO_4k vertex longitude (rad) */
#define LAL_LLO_4K_DETECTOR_LATITUDE_RAD       	0.53342313506	/**< LLO_4k vertex latitude (rad) */
#define LAL_LLO_4K_DETECTOR_ELEVATION_SI       	-6.574	/**< LLO_4k vertex elevation (m) */
#define LAL_LLO_4K_DETECTOR_ARM_X_AZIMUTH_RAD  	4.40317772346	/**< LLO_4k x arm azimuth (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_Y_AZIMUTH_RAD  	2.83238139666	/**< LLO_4k y arm azimuth (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_X_ALTITUDE_RAD 	-0.00031210000	/**< LLO_4k x arm altitude (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_Y_ALTITUDE_RAD 	-0.00061070000	/**< LLO_4k y arm altitude (rad) */
#define LAL_LLO_4K_DETECTOR_ARM_X_MIDPOINT_SI  	1997.57500000000	/**< LLO_4k x arm midpoint (m) */
#define LAL_LLO_4K_DETECTOR_ARM_Y_MIDPOINT_SI  	1997.57500000000	/**< LLO_4k y arm midpoint (m) */
#define LAL_LLO_4K_VERTEX_LOCATION_X_SI        	-7.42760447238e+04	/**< LLO_4k x-component of vertex location in Earth-centered frame (m) */
#define LAL_LLO_4K_VERTEX_LOCATION_Y_SI        	-5.49628371971e+06	/**< LLO_4k y-component of vertex location in Earth-centered frame (m) */
#define LAL_LLO_4K_VERTEX_LOCATION_Z_SI        	3.22425701744e+06	/**< LLO_4k z-component of vertex location in Earth-centered frame (m) */
#define LAL_LLO_4K_ARM_X_DIRECTION_X           	-0.95457412153	/**< LLO_4k x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_X_DIRECTION_Y           	-0.14158077340	/**< LLO_4k y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_X_DIRECTION_Z           	-0.26218911324	/**< LLO_4k z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_Y_DIRECTION_X           	0.29774156894	/**< LLO_4k x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_Y_DIRECTION_Y           	-0.48791033647	/**< LLO_4k y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LLO_4K_ARM_Y_DIRECTION_Z           	-0.82054461286	/**< LLO_4k z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */



/**
 * \name LIGO India 4km Interferometric Detector constants
 * @warning These numbers are subject to change.
 * The following constants describe hypothetical location and geometry
 * of the LIGO India 4km Interferometric Detector that have been used
 * in several studies with LALInference. Note that these data do not
 * represent an actual prospective site.
 */
/** @{ */
#define LAL_LIO_4K_DETECTOR_NAME                 "LIO_4k" /**< LIO_4K detector name string */
#define LAL_LIO_4K_DETECTOR_PREFIX               "I1"    /**< LIO_4K detector prefix string */
#define LAL_LIO_4K_DETECTOR_LONGITUDE_RAD        1.34444215058   /**< LIO_4K vertex longitude (rad; equal to 77°02') */
#define LAL_LIO_4K_DETECTOR_LATITUDE_RAD         0.34231676739   /**< LIO_4K vertex latitude (rad; equal to 19°37') */
#define LAL_LIO_4K_DETECTOR_ELEVATION_SI         440.0  /**< LIO_4K vertex elevation (m) */
#define LAL_LIO_4K_DETECTOR_ARM_X_AZIMUTH_RAD    5.80120119264   /**< LIO_4K x arm azimuth (rad) */
#define LAL_LIO_4K_DETECTOR_ARM_Y_AZIMUTH_RAD    4.23039066080   /**< LIO_4K y arm azimuth (rad) */
#define LAL_LIO_4K_DETECTOR_ARM_X_ALTITUDE_RAD   0.0   /**< LIO_4K x arm altitude (rad) */
#define LAL_LIO_4K_DETECTOR_ARM_Y_ALTITUDE_RAD   0.0   /**< LIO_4K y arm altitude (rad) */
#define LAL_LIO_4K_DETECTOR_ARM_X_MIDPOINT_SI    2000.00000000000        /**< LIO_4K x arm midpoint (m) */
#define LAL_LIO_4K_DETECTOR_ARM_Y_MIDPOINT_SI    2000.00000000000        /**< LIO_4K y arm midpoint (m) */
#define LAL_LIO_4K_VERTEX_LOCATION_X_SI          1.34897115479e+06       /**< LIO_4K x-component of vertex location in Earth-centered frame (m) */
#define LAL_LIO_4K_VERTEX_LOCATION_Y_SI          5.85742826577e+06       /**< LIO_4K y-component of vertex location in Earth-centered frame (m) */
#define LAL_LIO_4K_VERTEX_LOCATION_Z_SI          2.12756925209e+06       /**< LIO_4K z-component of vertex location in Earth-centered frame (m) */
#define LAL_LIO_4K_ARM_X_DIRECTION_X            0.38496278183  /**< LIO_4K x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_X_DIRECTION_Y             -0.39387275094   /**< LIO_4K y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_X_DIRECTION_Z            0.83466634811 /**< LIO_4K z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_Y_DIRECTION_X             0.89838844906  /**< LIO_4K x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_Y_DIRECTION_Y            -0.04722636126   /**< LIO_4K y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_LIO_4K_ARM_Y_DIRECTION_Z             -0.43665531647  /**< LIO_4K z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */


/**
 * \name Caltech 40m Prototype Detector constants
 * The following constants describe the location and geometry of the
 * Caltech 40m Prototype Detector.
 */
/** @{ */
#define LAL_CIT_40_DETECTOR_NAME               	"CIT_40"	/**< CIT_40 detector name string */
#define LAL_CIT_40_DETECTOR_PREFIX             	"C1"	/**< CIT_40 detector prefix string */
#define LAL_CIT_40_DETECTOR_LONGITUDE_RAD      	-2.06175744538	/**< CIT_40 vertex longitude (rad) */
#define LAL_CIT_40_DETECTOR_LATITUDE_RAD       	0.59637900541	/**< CIT_40 vertex latitude (rad) */
#define LAL_CIT_40_DETECTOR_ELEVATION_SI       	0	/**< CIT_40 vertex elevation (m) */
#define LAL_CIT_40_DETECTOR_ARM_X_AZIMUTH_RAD  	3.14159265359	/**< CIT_40 x arm azimuth (rad) */
#define LAL_CIT_40_DETECTOR_ARM_Y_AZIMUTH_RAD  	1.57079632679	/**< CIT_40 y arm azimuth (rad) */
#define LAL_CIT_40_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< CIT_40 x arm altitude (rad) */
#define LAL_CIT_40_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< CIT_40 y arm altitude (rad) */
#define LAL_CIT_40_DETECTOR_ARM_X_MIDPOINT_SI  	19.12500000000	/**< CIT_40 x arm midpoint (m) */
#define LAL_CIT_40_DETECTOR_ARM_Y_MIDPOINT_SI  	19.12500000000	/**< CIT_40 y arm midpoint (m) */
#define LAL_CIT_40_VERTEX_LOCATION_X_SI        	-2.49064958347e+06	/**< CIT_40 x-component of vertex location in Earth-centered frame (m) */
#define LAL_CIT_40_VERTEX_LOCATION_Y_SI        	-4.65869968211e+06	/**< CIT_40 y-component of vertex location in Earth-centered frame (m) */
#define LAL_CIT_40_VERTEX_LOCATION_Z_SI        	3.56206411403e+06	/**< CIT_40 z-component of vertex location in Earth-centered frame (m) */
#define LAL_CIT_40_ARM_X_DIRECTION_X           	-0.26480331633	/**< CIT_40 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_CIT_40_ARM_X_DIRECTION_Y           	-0.49530818538	/**< CIT_40 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_CIT_40_ARM_X_DIRECTION_Z           	-0.82737476706	/**< CIT_40 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_CIT_40_ARM_Y_DIRECTION_X           	0.88188012386	/**< CIT_40 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_CIT_40_ARM_Y_DIRECTION_Y           	-0.47147369718	/**< CIT_40 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_CIT_40_ARM_Y_DIRECTION_Z           	0.00000000000	/**< CIT_40 z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */


/**
 * \name Einstein Telescop 10km Interferometric Detector constants
 * The following constants describe the locations and geometrys of the
 * three 10km Interferometric Detectors for the planned third generation
 * Einstein Telescop detector as well as the theoretical null stream.
 * See T1400308 
 */
/** @{ */
#define LAL_ET1_DETECTOR_NAME                  	"ET1_T1400308"	/**< ET1 detector name string */
#define LAL_ET1_DETECTOR_PREFIX                	"E1"	/**< ET1 detector prefix string */
#define LAL_ET1_DETECTOR_LONGITUDE_RAD         	0.18333805213	/**< ET1 vertex longitude (rad) */
#define LAL_ET1_DETECTOR_LATITUDE_RAD          	0.76151183984	/**< ET1 vertex latitude (rad) */
#define LAL_ET1_DETECTOR_ELEVATION_SI          	51.884	/**< ET1 vertex elevation (m) */
#define LAL_ET1_DETECTOR_ARM_X_AZIMUTH_RAD     	0.33916285222	/**< ET1 x arm azimuth (rad) */
#define LAL_ET1_DETECTOR_ARM_Y_AZIMUTH_RAD     	5.57515060820	/**< ET1 y arm azimuth (rad) */
#define LAL_ET1_DETECTOR_ARM_X_ALTITUDE_RAD    	0.00000000000	/**< ET1 x arm altitude (rad) */
#define LAL_ET1_DETECTOR_ARM_Y_ALTITUDE_RAD    	0.00000000000	/**< ET1 y arm altitude (rad) */
#define LAL_ET1_DETECTOR_ARM_X_MIDPOINT_SI     	5000.00000000000	/**< ET1 x arm midpoint (m) */
#define LAL_ET1_DETECTOR_ARM_Y_MIDPOINT_SI     	5000.00000000000	/**< ET1 y arm midpoint (m) */
#define LAL_ET1_VERTEX_LOCATION_X_SI           	4.54637409900e+06	/**< ET1 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ET1_VERTEX_LOCATION_Y_SI           	8.42989697626e+05	/**< ET1 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ET1_VERTEX_LOCATION_Z_SI           	4.37857696241e+06	/**< ET1 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ET1_ARM_X_DIRECTION_X              	-0.70045821479	/**< ET1 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET1_ARM_X_DIRECTION_Y              	0.20848948619	/**< ET1 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET1_ARM_X_DIRECTION_Z              	0.68256166277	/**< ET1 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET1_ARM_Y_DIRECTION_X              	-0.39681482542	/**< ET1 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET1_ARM_Y_DIRECTION_Y              	-0.73500471881	/**< ET1 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET1_ARM_Y_DIRECTION_Z              	0.54982366052	/**< ET1 z-component of unit vector pointing along y arm in Earth-centered frame */


#define LAL_ET2_DETECTOR_NAME                  	"ET2_T1400308"	/**< ET2 detector name string */
#define LAL_ET2_DETECTOR_PREFIX                	"E2"	/**< ET2 detector prefix string */
#define LAL_ET2_DETECTOR_LONGITUDE_RAD         	0.18405858870	/**< ET2 vertex longitude (rad) */
#define LAL_ET2_DETECTOR_LATITUDE_RAD          	0.76299307990	/**< ET2 vertex latitude (rad) */
#define LAL_ET2_DETECTOR_ELEVATION_SI          	59.735	/**< ET2 vertex elevation (m) */
#define LAL_ET2_DETECTOR_ARM_X_AZIMUTH_RAD     	4.52845115854	/**< ET2 x arm azimuth (rad) */
#define LAL_ET2_DETECTOR_ARM_Y_AZIMUTH_RAD     	3.48125307555	/**< ET2 y arm azimuth (rad) */
#define LAL_ET2_DETECTOR_ARM_X_ALTITUDE_RAD    	-0.00078362156	/**< ET2 x arm altitude (rad) */
#define LAL_ET2_DETECTOR_ARM_Y_ALTITUDE_RAD    	-0.00157024452	/**< ET2 y arm altitude (rad) */
#define LAL_ET2_DETECTOR_ARM_X_MIDPOINT_SI     	5000.00000000000	/**< ET2 x arm midpoint (m) */
#define LAL_ET2_DETECTOR_ARM_Y_MIDPOINT_SI     	5000.00000000000	/**< ET2 y arm midpoint (m) */
#define LAL_ET2_VERTEX_LOCATION_X_SI           	4.53936951685e+06	/**< ET2 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ET2_VERTEX_LOCATION_Y_SI           	8.45074592488e+05	/**< ET2 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ET2_VERTEX_LOCATION_Z_SI           	4.38540257904e+06	/**< ET2 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ET2_ARM_X_DIRECTION_X              	0.30364338937	/**< ET2 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET2_ARM_X_DIRECTION_Y              	-0.94349420500	/**< ET2 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET2_ARM_X_DIRECTION_Z              	-0.13273800225	/**< ET2 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET2_ARM_Y_DIRECTION_X              	0.70045821479	/**< ET2 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET2_ARM_Y_DIRECTION_Y              	-0.20848948619	/**< ET2 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET2_ARM_Y_DIRECTION_Z              	-0.68256166277	/**< ET2 z-component of unit vector pointing along y arm in Earth-centered frame */


#define LAL_ET3_DETECTOR_NAME                  	"ET3_T1400308"	/**< ET3 detector name string */
#define LAL_ET3_DETECTOR_PREFIX                	"E3"	/**< ET3 detector prefix string */
#define LAL_ET3_DETECTOR_LONGITUDE_RAD         	0.18192996730	/**< ET3 vertex longitude (rad) */
#define LAL_ET3_DETECTOR_LATITUDE_RAD          	0.76270463257	/**< ET3 vertex latitude (rad) */
#define LAL_ET3_DETECTOR_ELEVATION_SI          	59.727	/**< ET3 vertex elevation (m) */
#define LAL_ET3_DETECTOR_ARM_X_AZIMUTH_RAD     	2.43258574281	/**< ET3 x arm azimuth (rad) */
#define LAL_ET3_DETECTOR_ARM_Y_AZIMUTH_RAD     	1.38538766217	/**< ET3 y arm azimuth (rad) */
#define LAL_ET3_DETECTOR_ARM_X_ALTITUDE_RAD    	-0.00156852107	/**< ET3 x arm altitude (rad) */
#define LAL_ET3_DETECTOR_ARM_Y_ALTITUDE_RAD    	-0.00078189811	/**< ET3 y arm altitude (rad) */
#define LAL_ET3_DETECTOR_ARM_X_MIDPOINT_SI     	5000.00000000000	/**< ET3 x arm midpoint (m) */
#define LAL_ET3_DETECTOR_ARM_Y_MIDPOINT_SI     	5000.00000000000	/**< ET3 y arm midpoint (m) */
#define LAL_ET3_VERTEX_LOCATION_X_SI           	4.54240595075e+06	/**< ET3 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ET3_VERTEX_LOCATION_Y_SI           	8.35639650438e+05	/**< ET3 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ET3_VERTEX_LOCATION_Z_SI           	4.38407519902e+06	/**< ET3 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ET3_ARM_X_DIRECTION_X              	0.39681482542	/**< ET3 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET3_ARM_X_DIRECTION_Y              	0.73500471881	/**< ET3 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET3_ARM_X_DIRECTION_Z              	-0.54982366052	/**< ET3 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET3_ARM_Y_DIRECTION_X              	-0.30364338937	/**< ET3 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET3_ARM_Y_DIRECTION_Y              	0.94349420500	/**< ET3 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET3_ARM_Y_DIRECTION_Z              	0.13273800225	/**< ET3 z-component of unit vector pointing along y arm in Earth-centered frame */


#define LAL_ET0_DETECTOR_NAME                  	"ET0_T1400308"	/**< ET0 detector name string */
#define LAL_ET0_DETECTOR_PREFIX                	"E0"	/**< ET0 detector prefix string */
#define LAL_ET0_DETECTOR_LONGITUDE_RAD         	0.18192996730	/**< ET0 vertex longitude (rad) */
#define LAL_ET0_DETECTOR_LATITUDE_RAD          	0.76270463257	/**< ET0 vertex latitude (rad) */
#define LAL_ET0_DETECTOR_ELEVATION_SI          	59.727	/**< ET0 vertex elevation (m) */
#define LAL_ET0_DETECTOR_ARM_X_AZIMUTH_RAD     	0.00000000000	/**< ET0 x arm azimuth (rad) */
#define LAL_ET0_DETECTOR_ARM_Y_AZIMUTH_RAD     	0.00000000000	/**< ET0 y arm azimuth (rad) */
#define LAL_ET0_DETECTOR_ARM_X_ALTITUDE_RAD    	0.00000000000	/**< ET0 x arm altitude (rad) */
#define LAL_ET0_DETECTOR_ARM_Y_ALTITUDE_RAD    	0.00000000000	/**< ET0 y arm altitude (rad) */
#define LAL_ET0_DETECTOR_ARM_X_MIDPOINT_SI     	0.00000000000	/**< ET0 x arm midpoint (m) */
#define LAL_ET0_DETECTOR_ARM_Y_MIDPOINT_SI     	0.00000000000	/**< ET0 y arm midpoint (m) */
#define LAL_ET0_VERTEX_LOCATION_X_SI           	4.54240595075e+06	/**< ET0 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ET0_VERTEX_LOCATION_Y_SI           	8.35639650438e+05	/**< ET0 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ET0_VERTEX_LOCATION_Z_SI           	4.38407519902e+06	/**< ET0 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ET0_ARM_X_DIRECTION_X              	0.00000000000	/**< ET0 x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET0_ARM_X_DIRECTION_Y              	0.00000000000	/**< ET0 y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET0_ARM_X_DIRECTION_Z              	0.00000000000	/**< ET0 z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ET0_ARM_Y_DIRECTION_X              	0.00000000000	/**< ET0 x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET0_ARM_Y_DIRECTION_Y              	0.00000000000	/**< ET0 y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ET0_ARM_Y_DIRECTION_Z              	0.00000000000	/**< ET0 z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */

/**
 * \name KAGRA Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * KAGRA Interferometric Detector.
 * \sa
 * > Yoshio Saito, "KAGRA location", KAGRA Technical Document JGW-G1503824
 * > http://gwdoc.icrr.u-tokyo.ac.jp/cgi-bin/DocDB/ShowDocument?docid=3824
 */
/** @{ */
#define LAL_KAGRA_DETECTOR_NAME               	"KAGRA"	/**< KAGRA detector name string */
#define LAL_KAGRA_DETECTOR_PREFIX             	"K1"	/**< KAGRA detector prefix string */
#define LAL_KAGRA_DETECTOR_LONGITUDE_RAD      	2.396441015	/**< KAGRA vertex longitude (rad) */
#define LAL_KAGRA_DETECTOR_LATITUDE_RAD       	0.6355068497	/**< KAGRA vertex latitude (rad) */
#define LAL_KAGRA_DETECTOR_ELEVATION_SI       	414.181	/**< KAGRA vertex elevation (m) */
#define LAL_KAGRA_DETECTOR_ARM_X_AZIMUTH_RAD  	1.054113	/**< KAGRA x arm azimuth (rad) */
#define LAL_KAGRA_DETECTOR_ARM_Y_AZIMUTH_RAD  	-0.5166798	/**< KAGRA y arm azimuth (rad) */
#define LAL_KAGRA_DETECTOR_ARM_X_ALTITUDE_RAD 	0.0031414	/**< KAGRA x arm altitude (rad) */
#define LAL_KAGRA_DETECTOR_ARM_Y_ALTITUDE_RAD 	-0.0036270	/**< KAGRA y arm altitude (rad) */
#define LAL_KAGRA_DETECTOR_ARM_X_MIDPOINT_SI  	1513.2535	/**< KAGRA x arm midpoint (m) */
#define LAL_KAGRA_DETECTOR_ARM_Y_MIDPOINT_SI  	1511.611	/**< KAGRA y arm midpoint (m) */
#define LAL_KAGRA_VERTEX_LOCATION_X_SI        	-3777336.024	/**< KAGRA x-component of vertex location in Earth-centered frame (m) */
#define LAL_KAGRA_VERTEX_LOCATION_Y_SI        	3484898.411	/**< KAGRA y-component of vertex location in Earth-centered frame (m) */
#define LAL_KAGRA_VERTEX_LOCATION_Z_SI        	3765313.697	/**< KAGRA z-component of vertex location in Earth-centered frame (m) */
#define LAL_KAGRA_ARM_X_DIRECTION_X           	-0.3759040	/**< KAGRA x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_KAGRA_ARM_X_DIRECTION_Y           	-0.8361583	/**< KAGRA y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_KAGRA_ARM_X_DIRECTION_Z           	0.3994189	/**< KAGRA z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_KAGRA_ARM_Y_DIRECTION_X           	0.7164378	/**< KAGRA x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_KAGRA_ARM_Y_DIRECTION_Y           	0.01114076	/**< KAGRA y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_KAGRA_ARM_Y_DIRECTION_Z           	0.6975620	/**< KAGRA z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */


/**
 * \name ACIGA Interferometric Detector constants (not implemented)
 * The following constants are stubs for the location and geometry of the
 * ACIGA Interferometric Detector.
 */
/** @{ */
#define LAL_ACIGA_DETECTOR_NAME               	"ACIGA"	/**< ACIGA detector name string */
#define LAL_ACIGA_DETECTOR_PREFIX             	"U1"	/**< ACIGA detector prefix string */
#define LAL_ACIGA_DETECTOR_LONGITUDE_RAD      	0.0	/**< ACIGA vertex longitude (rad) */
#define LAL_ACIGA_DETECTOR_LATITUDE_RAD       	0.0	/**< ACIGA vertex latitude (rad) */
#define LAL_ACIGA_DETECTOR_ELEVATION_SI       	0.0	/**< ACIGA vertex elevation (m) */
#define LAL_ACIGA_DETECTOR_ARM_X_AZIMUTH_RAD  	0.0	/**< ACIGA x arm azimuth (rad) */
#define LAL_ACIGA_DETECTOR_ARM_Y_AZIMUTH_RAD  	0.0	/**< ACIGA y arm azimuth (rad) */
#define LAL_ACIGA_DETECTOR_ARM_X_ALTITUDE_RAD 	0.0	/**< ACIGA x arm altitude (rad) */
#define LAL_ACIGA_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.0	/**< ACIGA y arm altitude (rad) */
#define LAL_ACIGA_DETECTOR_ARM_X_MIDPOINT_SI  	0.0	/**< ACIGA x arm midpoint (m) */
#define LAL_ACIGA_DETECTOR_ARM_Y_MIDPOINT_SI  	0.0	/**< ACIGA y arm midpoint (m) */
#define LAL_ACIGA_VERTEX_LOCATION_X_SI        	0.0	/**< ACIGA x-component of vertex location in Earth-centered frame (m) */
#define LAL_ACIGA_VERTEX_LOCATION_Y_SI        	0.0	/**< ACIGA y-component of vertex location in Earth-centered frame (m) */
#define LAL_ACIGA_VERTEX_LOCATION_Z_SI        	0.0	/**< ACIGA z-component of vertex location in Earth-centered frame (m) */
#define LAL_ACIGA_ARM_X_DIRECTION_X           	0.0	/**< ACIGA x-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ACIGA_ARM_X_DIRECTION_Y           	0.0	/**< ACIGA y-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ACIGA_ARM_X_DIRECTION_Z           	0.0	/**< ACIGA z-component of unit vector pointing along x arm in Earth-centered frame */
#define LAL_ACIGA_ARM_Y_DIRECTION_X           	0.0	/**< ACIGA x-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ACIGA_ARM_Y_DIRECTION_Y           	0.0	/**< ACIGA y-component of unit vector pointing along y arm in Earth-centered frame */
#define LAL_ACIGA_ARM_Y_DIRECTION_Z           	0.0	/**< ACIGA z-component of unit vector pointing along y arm in Earth-centered frame */
/** @} */


/* Resonant Mass (Bar) Detectors */


/**
 * \name ALLEGRO Resonant Mass Detector with 320 degree azimuth "IGEC axis" constants
 * The following constants describe the location and geometry of the
 * ALLEGRO Resonant Mass Detector with 320 degree azimuth "IGEC axis".
 */
/** @{ */
#define LAL_ALLEGRO_320_DETECTOR_NAME               	"ALLEGRO_320"	/**< ALLEGRO_320 detector name string */
#define LAL_ALLEGRO_320_DETECTOR_PREFIX             	"A1"	/**< ALLEGRO_320 detector prefix string */
#define LAL_ALLEGRO_320_DETECTOR_LONGITUDE_RAD      	-1.59137068496	/**< ALLEGRO_320 vertex longitude (rad) */
#define LAL_ALLEGRO_320_DETECTOR_LATITUDE_RAD       	0.53079879206	/**< ALLEGRO_320 vertex latitude (rad) */
#define LAL_ALLEGRO_320_DETECTOR_ELEVATION_SI       	0	/**< ALLEGRO_320 vertex elevation (m) */
#define LAL_ALLEGRO_320_DETECTOR_ARM_X_AZIMUTH_RAD  	-0.69813170080	/**< ALLEGRO_320 x arm azimuth (rad) */
#define LAL_ALLEGRO_320_DETECTOR_ARM_Y_AZIMUTH_RAD  	0.00000000000	/**< ALLEGRO_320 y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_ALLEGRO_320_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< ALLEGRO_320 x arm altitude (rad) */
#define LAL_ALLEGRO_320_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< ALLEGRO_320 y arm altitude (rad) UNUSED FOR BARS */
#define LAL_ALLEGRO_320_DETECTOR_ARM_X_MIDPOINT_SI  	0.00000000000	/**< ALLEGRO_320 x arm midpoint (m) UNUSED FOR BARS */
#define LAL_ALLEGRO_320_DETECTOR_ARM_Y_MIDPOINT_SI  	0.00000000000	/**< ALLEGRO_320 y arm midpoint (m) UNUSED FOR BARS */
#define LAL_ALLEGRO_320_VERTEX_LOCATION_X_SI        	-1.13258964140e+05	/**< ALLEGRO_320 x-component of vertex location in Earth-centered frame (m) */
#define LAL_ALLEGRO_320_VERTEX_LOCATION_Y_SI        	-5.50408337391e+06	/**< ALLEGRO_320 y-component of vertex location in Earth-centered frame (m) */
#define LAL_ALLEGRO_320_VERTEX_LOCATION_Z_SI        	3.20989567981e+06	/**< ALLEGRO_320 z-component of vertex location in Earth-centered frame (m) */
#define LAL_ALLEGRO_320_AXIS_DIRECTION_X            	-0.63467362345	/**< ALLEGRO_320 x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_ALLEGRO_320_AXIS_DIRECTION_Y            	0.40093077976	/**< ALLEGRO_320 y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_ALLEGRO_320_AXIS_DIRECTION_Z            	0.66063901000	/**< ALLEGRO_320 z-component of unit vector pointing along axis in Earth-centered frame */
/** @} */

/**
 * \name AURIGA Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * AURIGA Resonant Mass Detector.
 */
/** @{ */
#define LAL_AURIGA_DETECTOR_NAME               	"AURIGA"	/**< AURIGA detector name string */
#define LAL_AURIGA_DETECTOR_PREFIX             	"O1"	/**< AURIGA detector prefix string */
#define LAL_AURIGA_DETECTOR_LONGITUDE_RAD      	0.20853775679	/**< AURIGA vertex longitude (rad) */
#define LAL_AURIGA_DETECTOR_LATITUDE_RAD       	0.79156499342	/**< AURIGA vertex latitude (rad) */
#define LAL_AURIGA_DETECTOR_ELEVATION_SI       	0	/**< AURIGA vertex elevation (m) */
#define LAL_AURIGA_DETECTOR_ARM_X_AZIMUTH_RAD  	0.76794487088	/**< AURIGA x arm azimuth (rad) */
#define LAL_AURIGA_DETECTOR_ARM_Y_AZIMUTH_RAD  	0.00000000000	/**< AURIGA y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_AURIGA_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< AURIGA x arm altitude (rad) */
#define LAL_AURIGA_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< AURIGA y arm altitude (rad) UNUSED FOR BARS */
#define LAL_AURIGA_DETECTOR_ARM_X_MIDPOINT_SI  	0.00000000000	/**< AURIGA x arm midpoint (m) UNUSED FOR BARS */
#define LAL_AURIGA_DETECTOR_ARM_Y_MIDPOINT_SI  	0.00000000000	/**< AURIGA y arm midpoint (m) UNUSED FOR BARS */
#define LAL_AURIGA_VERTEX_LOCATION_X_SI        	4.39246733007e+06	/**< AURIGA x-component of vertex location in Earth-centered frame (m) */
#define LAL_AURIGA_VERTEX_LOCATION_Y_SI        	9.29508666967e+05	/**< AURIGA y-component of vertex location in Earth-centered frame (m) */
#define LAL_AURIGA_VERTEX_LOCATION_Z_SI        	4.51502913071e+06	/**< AURIGA z-component of vertex location in Earth-centered frame (m) */
#define LAL_AURIGA_AXIS_DIRECTION_X            	-0.64450412225	/**< AURIGA x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_AURIGA_AXIS_DIRECTION_Y            	0.57365538956	/**< AURIGA y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_AURIGA_AXIS_DIRECTION_Z            	0.50550364038	/**< AURIGA z-component of unit vector pointing along axis in Earth-centered frame */
/** @} */

/**
 * \name EXPLORER Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * EXPLORER Resonant Mass Detector.
 */
/** @{ */
#define LAL_EXPLORER_DETECTOR_NAME               	"EXPLORER"	/**< EXPLORER detector name string */
#define LAL_EXPLORER_DETECTOR_PREFIX             	"X1"	        /**< EXPLORER detector prefix string */
#define LAL_EXPLORER_DETECTOR_LONGITUDE_RAD      	0.10821041362	/**< EXPLORER vertex longitude (rad) */
#define LAL_EXPLORER_DETECTOR_LATITUDE_RAD       	0.81070543755	/**< EXPLORER vertex latitude (rad) */
#define LAL_EXPLORER_DETECTOR_ELEVATION_SI       	0	/**< EXPLORER vertex elevation (m) */
#define LAL_EXPLORER_DETECTOR_ARM_X_AZIMUTH_RAD  	0.68067840828	/**< EXPLORER x arm azimuth (rad) */
#define LAL_EXPLORER_DETECTOR_ARM_Y_AZIMUTH_RAD  	0.00000000000	/**< EXPLORER y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_EXPLORER_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< EXPLORER x arm altitude (rad) */
#define LAL_EXPLORER_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< EXPLORER y arm altitude (rad) UNUSED FOR BARS */
#define LAL_EXPLORER_DETECTOR_ARM_X_MIDPOINT_SI  	0.00000000000	/**< EXPLORER x arm midpoint (m) UNUSED FOR BARS */
#define LAL_EXPLORER_DETECTOR_ARM_Y_MIDPOINT_SI  	0.00000000000	/**< EXPLORER y arm midpoint (m) UNUSED FOR BARS */
#define LAL_EXPLORER_VERTEX_LOCATION_X_SI        	4.37645395452e+06	/**< EXPLORER x-component of vertex location in Earth-centered frame (m) */
#define LAL_EXPLORER_VERTEX_LOCATION_Y_SI        	4.75435044067e+05	/**< EXPLORER y-component of vertex location in Earth-centered frame (m) */
#define LAL_EXPLORER_VERTEX_LOCATION_Z_SI        	4.59985274450e+06	/**< EXPLORER z-component of vertex location in Earth-centered frame (m) */
#define LAL_EXPLORER_AXIS_DIRECTION_X            	-0.62792641437	/**< EXPLORER x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_EXPLORER_AXIS_DIRECTION_Y            	0.56480832712	/**< EXPLORER y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_EXPLORER_AXIS_DIRECTION_Z            	0.53544371484	/**< EXPLORER z-component of unit vector pointing along axis in Earth-centered frame */
/** @} */

/**
 * \name Nautilus Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * Nautilus Resonant Mass Detector.
 */
/** @{ */
#define LAL_NAUTILUS_DETECTOR_NAME               	"Nautilus"	/**< Nautilus detector name string */
#define LAL_NAUTILUS_DETECTOR_PREFIX             	"N1"	/**< Nautilus detector prefix string */
#define LAL_NAUTILUS_DETECTOR_LONGITUDE_RAD      	0.22117684946	/**< Nautilus vertex longitude (rad) */
#define LAL_NAUTILUS_DETECTOR_LATITUDE_RAD       	0.72996456710	/**< Nautilus vertex latitude (rad) */
#define LAL_NAUTILUS_DETECTOR_ELEVATION_SI       	0	/**< Nautilus vertex elevation (m) */
#define LAL_NAUTILUS_DETECTOR_ARM_X_AZIMUTH_RAD  	0.76794487088	/**< Nautilus x arm azimuth (rad) */
#define LAL_NAUTILUS_DETECTOR_ARM_Y_AZIMUTH_RAD  	0.00000000000	/**< Nautilus y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_NAUTILUS_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< Nautilus x arm altitude (rad) */
#define LAL_NAUTILUS_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< Nautilus y arm altitude (rad) UNUSED FOR BARS */
#define LAL_NAUTILUS_DETECTOR_ARM_X_MIDPOINT_SI  	0.00000000000	/**< Nautilus x arm midpoint (m) UNUSED FOR BARS */
#define LAL_NAUTILUS_DETECTOR_ARM_Y_MIDPOINT_SI  	0.00000000000	/**< Nautilus y arm midpoint (m) UNUSED FOR BARS */
#define LAL_NAUTILUS_VERTEX_LOCATION_X_SI        	4.64410999868e+06	/**< Nautilus x-component of vertex location in Earth-centered frame (m) */
#define LAL_NAUTILUS_VERTEX_LOCATION_Y_SI        	1.04425342477e+06	/**< Nautilus y-component of vertex location in Earth-centered frame (m) */
#define LAL_NAUTILUS_VERTEX_LOCATION_Z_SI        	4.23104713307e+06	/**< Nautilus z-component of vertex location in Earth-centered frame (m) */
#define LAL_NAUTILUS_AXIS_DIRECTION_X            	-0.62039441384	/**< Nautilus x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_NAUTILUS_AXIS_DIRECTION_Y            	0.57250373141	/**< Nautilus y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_NAUTILUS_AXIS_DIRECTION_Z            	0.53605060283	/**< Nautilus z-component of unit vector pointing along axis in Earth-centered frame */
/** @} */

/**
 * \name NIOBE Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * NIOBE Resonant Mass Detector.
 */
/** @{ */
#define LAL_NIOBE_DETECTOR_NAME               	"NIOBE"	/**< NIOBE detector name string */
#define LAL_NIOBE_DETECTOR_PREFIX             	"B1"	/**< NIOBE detector prefix string */
#define LAL_NIOBE_DETECTOR_LONGITUDE_RAD      	2.02138216202	/**< NIOBE vertex longitude (rad) */
#define LAL_NIOBE_DETECTOR_LATITUDE_RAD       	-0.55734180780	/**< NIOBE vertex latitude (rad) */
#define LAL_NIOBE_DETECTOR_ELEVATION_SI       	0	/**< NIOBE vertex elevation (m) */
#define LAL_NIOBE_DETECTOR_ARM_X_AZIMUTH_RAD  	0.00000000000	/**< NIOBE x arm azimuth (rad) */
#define LAL_NIOBE_DETECTOR_ARM_Y_AZIMUTH_RAD  	0.00000000000	/**< NIOBE y arm azimuth (rad) UNUSED FOR BARS */
#define LAL_NIOBE_DETECTOR_ARM_X_ALTITUDE_RAD 	0.00000000000	/**< NIOBE x arm altitude (rad) */
#define LAL_NIOBE_DETECTOR_ARM_Y_ALTITUDE_RAD 	0.00000000000	/**< NIOBE y arm altitude (rad) UNUSED FOR BARS */
#define LAL_NIOBE_DETECTOR_ARM_X_MIDPOINT_SI  	0.00000000000	/**< NIOBE x arm midpoint (m) UNUSED FOR BARS */
#define LAL_NIOBE_DETECTOR_ARM_Y_MIDPOINT_SI  	0.00000000000	/**< NIOBE y arm midpoint (m) UNUSED FOR BARS */
#define LAL_NIOBE_VERTEX_LOCATION_X_SI        	-2.35948871453e+06	/**< NIOBE x-component of vertex location in Earth-centered frame (m) */
#define LAL_NIOBE_VERTEX_LOCATION_Y_SI        	4.87721571259e+06	/**< NIOBE y-component of vertex location in Earth-centered frame (m) */
#define LAL_NIOBE_VERTEX_LOCATION_Z_SI        	-3.35416003274e+06	/**< NIOBE z-component of vertex location in Earth-centered frame (m) */
#define LAL_NIOBE_AXIS_DIRECTION_X            	-0.23034623759	/**< NIOBE x-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_NIOBE_AXIS_DIRECTION_Y            	0.47614056486	/**< NIOBE y-component of unit vector pointing along axis in Earth-centered frame */
#define LAL_NIOBE_AXIS_DIRECTION_Z            	0.84866411101	/**< NIOBE z-component of unit vector pointing along axis in Earth-centered frame */
/** @} */

/** @} */ /* end: DetectorConstants */

#ifdef __cplusplus
}
#endif

#endif /* _LALDETECTORS_H */
