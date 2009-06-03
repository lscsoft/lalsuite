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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/********************************* <lalVerbatim file="LALDetectorsHV">
Author: J. T. Whelan <john.whelan@ligo.org>
$Id$
********************************** </lalVerbatim> */

/********************************* <lalLaTeX>

\section{Header \texttt{LALDetectors.h}}
\label{tools:s:LALDetectors.h}


This header defines structures to hold the basic data describing
a gravitational wave detector.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALDetectors.h>
\end{verbatim}

According to the common frame format specification
\cite{tools:LIGOVIRGO:2000} the geometry of an interferometric
detector will be stored in a \texttt{FrDetector} structure, specifying
the location of the detector vertex and the orientation of its arms in
geodetic co\"{o}rdinates suited to geographical surveying.  Resonant
bars and other sorts of detectors, if they write their data to frames,
are expected to fill this structure with their location and
orientation in some way suited to the detector type.

For most data analysis tasks, however, any gravitational wave detector
can be described by its location in an Earth-fixed rotating reference
frame, as well as a \textit{response tensor} $d^{ab}$, constant in the
same frame, which defines the ``strain'' $h$ measured by the detector in
terms of the metric perturbation $h_{ab}$ as
\begin{equation}
h = h_{ab} \, d^{ab}
\ .
\end{equation}

This header defines a \texttt{LALFrDetector} structure which contains
essentially the same information as the \texttt{FrDetector} structure,
as well as a \texttt{LALDetector} structure which contains the
Cartesian co\"{o}rdinates of the detector along with the components of
the response tensor $d^{ab}$ in the same co\"{o}rdinate system.

\subsubsection*{The Geodetic Co\"{o}rdinate System}

Geodetic co\"{o}rdinates are spheroidal co\"{o}rdinates
based on the WGS-84 Earth Model, which is an
oblate spheroid with equatorial radius $a=6.378137\times
10^6\,\textrm{m}$ and polar radius $b=6.356752314\times
10^6\,\textrm{m}$.  Any point in space can be located according to its
longitude, latitude, and elevation.  The \textit{longitude} $\lambda$
is the angle between the half-plane bounded by the symmetry axis of
the reference ellipsoid containing the point in question and the
half-plane plane containing the Prime Meridian; it is measured in
radians, increases to the East, and ranges from
$-\pi$ to $\pi$.  The \textit{latitude} $\beta$ is the
angle between the ray which is normal to the ellipsoid and passes
through the point in question and the equatorial plane; it is measured
in radians, increases to the North, and ranges
from $-\pi/2$ to $\pi/2$.  The \textit{elevation} $h$ is the
signed distance along this ray from the reference ellipsoid to the
point in question.  This co\"{o}rdinate system is described in more
detail in \cite{tools:Althouse:1999}.

\subsubsection*{Altitude and Azimuth Angles}

The \texttt{LALFrDetector} structure stores the directions along the
two arms of an interferometer in an altitude/azimuth representation
with respect to the local tangent plane to the reference ellipsoid,
known as the local horizontal.  The altitude ${\mathcal{A}}$ is the angle the
direction vector makes with the horizontal, ${\mathcal{A}} > 0$ meaning above
horizontal, ${\mathcal{A}} < 0$ below.  The azimuth angle $\zeta$ is found by
projecting the direction onto the local horizontal plane, then
measuring the angle clockwise from North to this projected direction.

\subsubsection*{The Cartesian Co\"{o}rdinate System}

The position vector and response tensor contained in the
\texttt{LALDetector} structure are defined in
 a simple orthonormal co\"{o}rdinate system with its origin at
the center of the earth, an $x^1$ axis which pierces the Earth's
surface at the intersection of the equator and the prime meridian, an
$x^2$ axis which pierces the earth's surface at $\pi/2$ radians East
longitude on the equator, and an $x^3$ axis which pierces the Earth's
surface at the North Pole.  The co\"{o}rdinates $x^1$, $x^2$, $x^3$
correspond to the Earth-fixed co\"{o}rdinates $X_E$, $Y_E$, $Z_E$
defined in \cite{tools:Althouse:1999}, respectively.

The relationship between geodetic and Cartesian co\"{o}rdinates is
given by
\begin{eqnarray}
\label{tools:e:cart1}
x^1&=&\left(
          \frac{a^2}{\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}}
          + h
       \right) \cos\beta\cos\lambda             \\
\label{tools:e:cart2}
x^2&=&\left(
          \frac{a^2}{\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}}
          + h
       \right) \cos\beta\sin\lambda             \\
\label{tools:e:cart3}
x^3&=&\left(
          \frac{b^2}{\sqrt{a^2\cos^2\beta+b^2\sin^2\beta}}
          + h
       \right) \sin\beta
\end{eqnarray}

\subsection*{Error conditions}
\input{LALDetectorsHE}

********************************** </lalLaTeX> */

/**************************************** <lalLaTeX file="LALDetectorsHB">

\bibitem{tools:LIGOVIRGO:2000}
LIGO Data Group and VIRGO Data Acquisition Group, ``Specification of a
Common Data Frame Format for Interferometric Gravitational Wave
Detectors (IGWD)'', LIGO Technical Note
\href{http://www.ligo.caltech.edu/docs/T/T970130-D.pdf}{LIGO-T970130}

\bibitem{tools:Althouse:1999}
William Althouse, Larry Jones, and Albert Lazzarini, ``Determination
of Global and Local Coordinate Axes for the LIGO Sites'', LIGO
Technical Note
\href{http://www.ligo.caltech.edu/docs/T/T980044-10.pdf}{LIGO-T980044}

% \bibitem{tools:Lazzarini:1995}
% Albert Lazzarini, ``Derivation
% of Global and Local Coordinate Axes for the LIGO Sites'', LIGO
% Technical Note LIGO-T950004

******************************************************* </lalLaTeX> */


#ifndef _LALDETECTORS_H
#define _LALDETECTORS_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( LALDETECTORSH, "$Id$" );

/** \file
 * \ingroup tools
 * \author J. T. Whelan and J. D. E. Creighton
 * \date $Id$
 * \brief Provides constants describing various gravitational wave detectors
 *
 * This header defines numerical constants that describe the location and
 * geometry of several operating gravitational wave detectors.
 * These detectors are both resonant mass (bar) detectors and interferometric
 * detectors.  Data for the resonant mass detectors is taken from:
 *
 *      http://igec.lnl.infn.it/cgi-bin/browser.pl?Level=0,3,1
 *
 * and
 *
 *      L. S. Finn and A. Lazzarini, Phys. Rev. D 64, 082002 (2001)
 *
 * Data for LIGO detectors is taken from:
 *
 *      William Althouse, Larry Jones, Albert Lazzarini (1999)
 *      "Determination of Global and Local Coordinate Axes for the LIGO Sites"
 *      LIGO-T980044-08-E
 *
 * Data for the VIRGO detector is provided by Benoit Mours.
 *
 * Data for the GEO detector is taken from:
 *
 *      http://www.geo600.uni-hannover.de/geo600/project/location.html
 *
 * Data for the TAMA detector is provided by Masa-Katsu Fujimoto
 *
 * Data for the Caltech detector is taken from:
 *
 *      B. Allen, "Gravitational Wave Detector Sites," gr-qc/9607075 (1996).
 *
 * See the technical document
 *
 *      Warren Anderson, Patrick Brady, David Chin, Jolien Creighton,
 *      Keith Riles, and John Whelan
 *      "Beam Pattern Response Functions and Times of Arrival
 *      for Earthbound Interferometer"
 *      LIGO-T010110-00-Z
 *      http://www.lsc-group.phys.uwm.edu/daswg/docs/technical/T010110.pdf
 *
 * for details.
 *
 * Data in this file (e.g., angle conventions etc.) is intended
 * to conform to the conventions of the Frame format specification:
 *
 *      LIGO Data and Computing Group and Virgo Data Acquisition Group
 *      Specification of a Common Data Frame Format for
 *      Interferometric Gravitational Wave Detectors
 *      (IGWD)
 *      LIGO-T970130-F-E and VIRGO-SPE-LAP-5400-102 (Version 6)
 *      http://www.ligo.caltech.edu/docs/T/T970130-F.pdf
 *
 */


/******************************** <lalErrTable file="LALDetectorsHE"> */

#define LALDETECTORSH_ENULLP        1
#define LALDETECTORSH_ETYPE         2

#define LALDETECTORSH_MSGENULLP     "Null pointer"
#define LALDETECTORSH_MSGETYPE      "Unsupported detector type"

/************************************ </lalErrTable> */

#define LALDETECTORSH_PRINTF        0

/********************************* <lalLaTeX>

\subsubsection*{The \texttt{LALDetectorType} enumeration}
\idx[Type]{LALDetectorType}
\idx[Constant]{LALDETECTORTYPE\_ABSENT}
\idx[Constant]{LALDETECTORTYPE\_IFODIFF}
\idx[Constant]{LALDETECTORTYPE\_IFOXARM}
\idx[Constant]{LALDETECTORTYPE\_IFOYARM}
\idx[Constant]{LALDETECTORTYPE\_IFOCOMM}
\idx[Constant]{LALDETECTORTYPE\_CYLBAR}

Since data from bars as well as interferometers can be written to
   frames, we need an additional piece of information to interpret the
   site geometry data specified in the \texttt{LALFrDetector}
   structure; for instance, is the x arm really the x arm or is it the
   long axis of a bar?  The \texttt{LALDetectorType} enumeration
   provides a way to keep track of that.

The possible values are (each value is prefaced by
\texttt{LALDETECTORTYPE\_}):
\begin{description}
  \item[\texttt{LALDETECTORTYPE\_ABSENT}] No \texttt{FrDetector}
	associated with the structure
  \item[\texttt{LALDETECTORTYPE\_IFODIFF}] Interferometer
	in differential mode
  \item[\texttt{LALDETECTORTYPE\_IFOXARM}] Interferometer
	in one-armed mode (X arm)
  \item[\texttt{LALDETECTORTYPE\_IFOYARM}] Interferometer
	in one-armed mode (Y arm)
  \item[\texttt{LALDETECTORTYPE\_IFOCOMM}] Interferometer in common mode
  \item[\texttt{LALDETECTORTYPE\_CYLBAR}]  Cylindrical bar
\end{description}

\subsubsection*{Cached Detectors}

\idx[Variable]{lalCachedDetectors[]} In practice, we will often be
working with fixed unchanging site geometry, e.g., for the LIGO
interferometers; to avoid constantly reconstructing the corresponding
\texttt{LALDetector}s, we should define some constant
\texttt{LALDetector}s describing them.  Those are stored in a constant
array of \texttt{LALDetector} structures known as
\texttt{lalCachedDetectors}, which is declared \texttt{extern} in this
header and defined in \texttt{CreateDetector.c} (see
Sec.~\ref{tools:ss:CreateDetector.c}).

The \texttt{LALCreateDetector()} routine will first look through the
\texttt{lalCachedDetectors} array for a \texttt{LALDetector} structure
with matching \texttt{type} and \texttt{frDetector.name} fields; if it
finds one, it returns a copy of that; if not, it creates one.

The header \texttt{LALDetectors.h} also defines an enumeration of the
indices of the known detectors:
\idx[Constant]{LAL\_TAMA\_300\_DETECTOR}
\idx[Constant]{LAL\_VIRGO\_DETECTOR}
\idx[Constant]{LAL\_GEO\_600\_DETECTOR}
\idx[Constant]{LAL\_LHO\_2K\_DETECTOR}
\idx[Constant]{LAL\_LHO\_4K\_DETECTOR}
\idx[Constant]{LAL\_LLO\_4K\_DETECTOR}
\idx[Constant]{LAL\_CIT\_40\_DETECTOR}
\idx[Constant]{LAL\_ALLEGRO\_DETECTOR}
\idx[Constant]{LAL\_AURIGA\_DETECTOR}
\idx[Constant]{LAL\_EXPLORER\_DETECTOR}
\idx[Constant]{LAL\_NIOBE\_DETECTOR}
\idx[Constant]{LAL\_NAUTILUS\_DETECTOR}
\idx[Constant]{LAL\_NUM\_DETECTORS}
********************************** </lalLaTeX> */

  /********************************* <lalLaTeX>
For example, the \texttt{LALDetector} representing LIGO Hanford 4km (H1) in
differential mode is
\texttt{lalCachedDetectors[LAL\_LHO\_4K\_DETECTOR]}.

\subsection*{Structures}

********************************** </lalLaTeX> */

/********************************* <lalLaTeX>

\subsubsection*{Structure \texttt{LALFrDetector}}
\idx[Type]{LALFrDetector}

The \texttt{LALFrDetector} structure holds site geometry information
in the same format as the \texttt{FrDetector} structure defined in the
frames spec. \cite{tools:LIGOVIRGO:2000}  The fields are:
\begin{description}
  \item[\texttt{CHAR name[LALNameLength]}] A unique identifying string.
  \item[\texttt{CHAR prefix[3]}] Two-letter prefix for detector names.
  \item[\texttt{REAL8 vertexLongitudeRadians}] The geodetic longitude
$\lambda$ of the vertex, in radians.
  \item[\texttt{REAL8 vertexLatitudeRadians}] The geodetic latitude
$\beta$ of the vertex, in radians.
  \item[\texttt{REAL4 vertexElevation}] The height of the vertex above
  the reference ellipsoid, in meters.
  \item[\texttt{REAL4 xArmAltitudeRadians}]  The angle ${\mathcal{A}}_X$ up from the
  local tangent plane of the reference ellipsoid to the X arm, in radians.
  \item[\texttt{REAL4 xArmAzimuthRadians}] The angle $\zeta_X$ clockwise
  from North to the projection of the X arm into the local tangent plane of
  the reference ellipsoid, in radians.
  \item[\texttt{REAL4 yArmAltitudeRadians}]  The angle ${\mathcal{A}}_Y$ up from the
  local tangent plane of the reference ellipsoid to the Y arm, in radians.
  \item[\texttt{REAL4 yArmAzimuthRadians}] The angle $\zeta_Y$ clockwise
  from North to the projection of the Y arm into the local tangent plane of
  the reference ellipsoid, in radians.
  \item[\texttt{REAL4 xArmMidpoint}] The distance to the midpoint of the X arm in meters (unused for bars: set it to zero).
  \item[\texttt{REAL4 yArmMidpoint}] The distance to the midpoint of the Y arm in meters (unused for bars: set it to zero).
\end{description}

\subsubsection*{Structure \texttt{LALDetector}}
\idx[Type]{LALDetector}

The \texttt{LALDetector} structure is intended to be the way that detector
   geometry information is passed to LAL routines.
This structure describes a detector geometry in a way independent of
the type of detector.  The fields are:
\begin{description}
  \item[\texttt{REAL8 location[3]}]  The three components, in an
  Earth-fixed Cartesian co\"{o}rdinate system, of the
  position vector from the center of the Earth to the detector,
  in meters.
  \item[\texttt{REAL4 response[3][3]}] The Earth-fixed Cartesian components
 of the detector's response tensor   $d^{ab}$.
  \item[\texttt{LALDetectorType type}] The type of detector (e.g., IFO in
  differential mode, cylindrical bar, etc.)
  \item[\texttt{LALFrDetector frDetector}] The original
  \texttt{LALFrDetector} structure from which this was created.
\end{description}

\vfill{\footnotesize\input{LALDetectorsHV}}
\newpage\input{CreateDetectorC}
\newpage\input{DetectorSiteTestC}

********************************** </lalLaTeX> */



/** Enumeration of Detectors: follows order of DQ bit assignments */
enum {
	LAL_TAMA_300_DETECTOR	=	0,
	LAL_VIRGO_DETECTOR	=	1,
	LAL_GEO_600_DETECTOR	=	2,
	LAL_LHO_2K_DETECTOR	=	3,
	LAL_LHO_4K_DETECTOR	=	4,
	LAL_LLO_4K_DETECTOR	=	5,
	LAL_CIT_40_DETECTOR	=	6,
	LAL_ALLEGRO_DETECTOR	=	7,
	LAL_AURIGA_DETECTOR	=	8,
	LAL_EXPLORER_DETECTOR	=	9,
	LAL_NIOBE_DETECTOR	=	10,
	LAL_NAUTILUS_DETECTOR	=	11,
	LAL_NUM_DETECTORS	=	12
};

/**< Detector DQ bit assignments (2 bits per detector) */
enum {
	LAL_TAMA_300_DETECTOR_BIT	=	1 << 2 * LAL_TAMA_300_DETECTOR,
	LAL_VIRGO_DETECTOR_BIT   	=	1 << 2 * LAL_VIRGO_DETECTOR,
	LAL_GEO_600_DETECTOR_BIT 	=	1 << 2 * LAL_GEO_600_DETECTOR,
	LAL_LHO_2K_DETECTOR_BIT  	=	1 << 2 * LAL_LHO_2K_DETECTOR,
	LAL_LHO_4K_DETECTOR_BIT  	=	1 << 2 * LAL_LHO_4K_DETECTOR,
	LAL_LLO_4K_DETECTOR_BIT  	=	1 << 2 * LAL_LLO_4K_DETECTOR,
	LAL_CIT_40_DETECTOR_BIT  	=	1 << 2 * LAL_CIT_40_DETECTOR,
	LAL_ALLEGRO_DETECTOR_BIT 	=	1 << 2 * LAL_ALLEGRO_DETECTOR,
	LAL_AURIGA_DETECTOR_BIT  	=	1 << 2 * LAL_AURIGA_DETECTOR,
	LAL_NIOBE_DETECTOR_BIT   	=	1 << 2 * LAL_NIOBE_DETECTOR,
	LAL_NAUTILUS_DETECTOR_BIT	=	1 << 2 * LAL_NAUTILUS_DETECTOR
};


/** Detector type
 *
 * The type of detector.  This determines how the detector response
 * is determined.
 */
typedef enum {
	LALDETECTORTYPE_ABSENT,	/**< No FrDetector associated with this detector */
	LALDETECTORTYPE_IFODIFF,	/**< IFO in differential mode */
	LALDETECTORTYPE_IFOXARM,	/**< IFO in one-armed mode (X arm) */
	LALDETECTORTYPE_IFOYARM,	/**< IFO in one-armed mode (Y arm) */
	LALDETECTORTYPE_IFOCOMM,	/**< IFO in common mode */
	LALDETECTORTYPE_CYLBAR	/**< Cylindrical bar */
}
LALDetectorType;


/** Detector frame data structure
 *
 * Structure to contain the data that appears in a FrDetector structure
 * in frame data.
 */
typedef struct tagLALFrDetector
{
	CHAR	name[LALNameLength];	/**< A unique identifying string. */
	CHAR	prefix[3];		/**< Two-letter prefix for detector's channel names. */
	REAL8	vertexLongitudeRadians;	/**< The geodetic longitude \f$\lambda\f$ of the vertex in radians. */
	REAL8	vertexLatitudeRadians;	/**< The geodetic latitude \f$\beta\f$ of the vertex in radians. */
	REAL4	vertexElevation;	/**< The height of the vertex above the reference ellipsoid in meters. */
	REAL4	xArmAltitudeRadians;	/**< The angle \f${\mathcal{A}}_X\f$ up from the local tangent plane of the reference ellipsoid to the X arm (or bar's cylidrical axis) in radians. */
	REAL4	xArmAzimuthRadians;	/**< The angle \f$\zeta_X\f$ clockwise from North to the projection of the X arm (or bar's cylidrical axis) into the local tangent plane of the reference ellipsoid in radians. */
	REAL4	yArmAltitudeRadians;	/**< The angle \f${\mathcal{A}}_Y\f$ up from the local tangent plane of the reference ellipsoid to the Y arm in radians (unused for bars: set it to zero). */
	REAL4	yArmAzimuthRadians;	/**< The angle \f$\zeta_Y\f$ clockwise from North to the projection of the Y arm into the local tangent plane of the reference ellipsoid in radians (unused for bars: set it to zero). */
	REAL4	xArmMidpoint;	/**< The distance to the midpoint of the X arm in meters (unused for bars: set it to zero). */
	REAL4	yArmMidpoint;	/**< The distance to the midpoint of the Y arm in meters (unused for bars: set it to zero). */
}
LALFrDetector;


/** Detector structure
 *
 * Structure to contain detector data in the format most easily used
 * by the LAL routines.
 */
typedef struct tagLALDetector
{
	REAL8		location[3];	/**< The three components, in an Earth-fixed Cartesian coordinate system, of the position vector from the center of the Earth to the detector in meters. */
	REAL4		response[3][3];	/**< The Earth-fixed Cartesian components of the detector's response tensor \f$d^{ab}\f$. */
	LALDetectorType	type;		/**< The type of the detector (e.g., IFO in differential mode, cylindrical bar, etc.). */
	LALFrDetector	frDetector;	/**< The original LALFrDetector structure from which this was created. */
}
LALDetector;


/** Pre-existing detectors. */
extern const LALDetector lalCachedDetectors[LAL_NUM_DETECTORS];



/** Routine to create a LALDetector. */
LALDetector * XLALCreateDetector( LALDetector *detector, const LALFrDetector *frDetector, LALDetectorType type );
void LALCreateDetector( LALStatus *status, LALDetector *output, const LALFrDetector *input, const LALDetectorType type );



/* Interferometric Detectors */


/** \name TAMA 300m Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * TAMA 300m Interferometric Detector. */
/*@{*/
#define LAL_TAMA_300_DETECTOR_NAME               	"TAMA_300"	/**< TAMA_300 detector name string */
#define LAL_TAMA_300_DETECTOR_PREFIX             	"T1"	/**< TAMA_300 detector prefix string */
#define LAL_TAMA_300_DETECTOR_LONGITUDE_RAD      	2.43536359469	/**< TAMA_300 vertex longitude (rad) */
#define LAL_TAMA_300_DETECTOR_LATITUDE_RAD       	0.62267336022	/**< TAMA_300 vertex latitude (rad) */
#define LAL_TAMA_300_DETECTOR_ELEVATION_SI       	90	/**< TAMA_300 vertex elevation (rad) */
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
/*@}*/


/** \name VIRGO 3km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * VIRGO 3km Interferometric Detector. */
/*@{*/
#define LAL_VIRGO_DETECTOR_NAME               	"VIRGO"	/**< VIRGO detector name string */
#define LAL_VIRGO_DETECTOR_PREFIX             	"V2"	/**< VIRGO detector prefix string */
#define LAL_VIRGO_DETECTOR_LONGITUDE_RAD      	0.18333805213	/**< VIRGO vertex longitude (rad) */
#define LAL_VIRGO_DETECTOR_LATITUDE_RAD       	0.76151183984	/**< VIRGO vertex latitude (rad) */
#define LAL_VIRGO_DETECTOR_ELEVATION_SI       	51.884	/**< VIRGO vertex elevation (rad) */
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
/*@}*/


/** \name GEO 600m Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * GEO 600m Interferometric Detector. */
/*@{*/
#define LAL_GEO_600_DETECTOR_NAME               	"GEO_600"	/**< GEO_600 detector name string */
#define LAL_GEO_600_DETECTOR_PREFIX             	"G1"	/**< GEO_600 detector prefix string */
#define LAL_GEO_600_DETECTOR_LONGITUDE_RAD      	0.17116780435	/**< GEO_600 vertex longitude (rad) */
#define LAL_GEO_600_DETECTOR_LATITUDE_RAD       	0.91184982752	/**< GEO_600 vertex latitude (rad) */
#define LAL_GEO_600_DETECTOR_ELEVATION_SI       	114.425	/**< GEO_600 vertex elevation (rad) */
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
/*@}*/


/** \name LIGO Hanford Observatory 2km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Hanford Observatory 2km Interferometric Detector. */
/*@{*/
#define LAL_LHO_2K_DETECTOR_NAME               	"LHO_2k"	/**< LHO_2k detector name string */
#define LAL_LHO_2K_DETECTOR_PREFIX             	"H2"	/**< LHO_2k detector prefix string */
#define LAL_LHO_2K_DETECTOR_LONGITUDE_RAD      	-2.08405676917	/**< LHO_2k vertex longitude (rad) */
#define LAL_LHO_2K_DETECTOR_LATITUDE_RAD       	0.81079526383	/**< LHO_2k vertex latitude (rad) */
#define LAL_LHO_2K_DETECTOR_ELEVATION_SI       	142.554	/**< LHO_2k vertex elevation (rad) */
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
/*@}*/


/** \name LIGO Hanford Observatory 4km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Hanford Observatory 4km Interferometric Detector. */
/*@{*/
#define LAL_LHO_4K_DETECTOR_NAME               	"LHO_4k"	/**< LHO_4k detector name string */
#define LAL_LHO_4K_DETECTOR_PREFIX             	"H1"	/**< LHO_4k detector prefix string */
#define LAL_LHO_4K_DETECTOR_LONGITUDE_RAD      	-2.08405676917	/**< LHO_4k vertex longitude (rad) */
#define LAL_LHO_4K_DETECTOR_LATITUDE_RAD       	0.81079526383	/**< LHO_4k vertex latitude (rad) */
#define LAL_LHO_4K_DETECTOR_ELEVATION_SI       	142.554	/**< LHO_4k vertex elevation (rad) */
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
/*@}*/


/** \name LIGO Livingston Observatory 4km Interferometric Detector constants
 * The following constants describe the location and geometry of the
 * LIGO Livingston Observatory 4km Interferometric Detector. */
/*@{*/
#define LAL_LLO_4K_DETECTOR_NAME               	"LLO_4k"	/**< LLO_4k detector name string */
#define LAL_LLO_4K_DETECTOR_PREFIX             	"L1"	/**< LLO_4k detector prefix string */
#define LAL_LLO_4K_DETECTOR_LONGITUDE_RAD      	-1.58430937078	/**< LLO_4k vertex longitude (rad) */
#define LAL_LLO_4K_DETECTOR_LATITUDE_RAD       	0.53342313506	/**< LLO_4k vertex latitude (rad) */
#define LAL_LLO_4K_DETECTOR_ELEVATION_SI       	-6.574	/**< LLO_4k vertex elevation (rad) */
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
/*@}*/


/** \name Caltech 40m Prototype Detector constants
 * The following constants describe the location and geometry of the
 * Caltech 40m Prototype Detector. */
/*@{*/
#define LAL_CIT_40_DETECTOR_NAME               	"CIT_40"	/**< CIT_40 detector name string */
#define LAL_CIT_40_DETECTOR_PREFIX             	"P1"	/**< CIT_40 detector prefix string */
#define LAL_CIT_40_DETECTOR_LONGITUDE_RAD      	-2.06175744538	/**< CIT_40 vertex longitude (rad) */
#define LAL_CIT_40_DETECTOR_LATITUDE_RAD       	0.59637900541	/**< CIT_40 vertex latitude (rad) */
#define LAL_CIT_40_DETECTOR_ELEVATION_SI       	0	/**< CIT_40 vertex elevation (rad) */
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
/*@}*/


/* Resonant Mass (Bar) Detectors */


/** \name ALLEGRO Resonant Mass Detector with 320 degree azimuth "IGEC axis" constants
 * The following constants describe the location and geometry of the
 * ALLEGRO Resonant Mass Detector with 320 degree azimuth "IGEC axis". */
/*@{*/
#define LAL_ALLEGRO_320_DETECTOR_NAME               	"ALLEGRO_320"	/**< ALLEGRO_320 detector name string */
#define LAL_ALLEGRO_320_DETECTOR_PREFIX             	"A1"	/**< ALLEGRO_320 detector prefix string */
#define LAL_ALLEGRO_320_DETECTOR_LONGITUDE_RAD      	-1.59137068496	/**< ALLEGRO_320 vertex longitude (rad) */
#define LAL_ALLEGRO_320_DETECTOR_LATITUDE_RAD       	0.53079879206	/**< ALLEGRO_320 vertex latitude (rad) */
#define LAL_ALLEGRO_320_DETECTOR_ELEVATION_SI       	0	/**< ALLEGRO_320 vertex elevation (rad) */
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


/** \name AURIGA Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * AURIGA Resonant Mass Detector. */
/*@{*/
#define LAL_AURIGA_DETECTOR_NAME               	"AURIGA"	/**< AURIGA detector name string */
#define LAL_AURIGA_DETECTOR_PREFIX             	"O1"	/**< AURIGA detector prefix string */
#define LAL_AURIGA_DETECTOR_LONGITUDE_RAD      	0.20853775679	/**< AURIGA vertex longitude (rad) */
#define LAL_AURIGA_DETECTOR_LATITUDE_RAD       	0.79156499342	/**< AURIGA vertex latitude (rad) */
#define LAL_AURIGA_DETECTOR_ELEVATION_SI       	0	/**< AURIGA vertex elevation (rad) */
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


/** \name EXPLORER Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * EXPLORER Resonant Mass Detector. */
/*@{*/
#define LAL_EXPLORER_DETECTOR_NAME               	"EXPLORER"	/**< EXPLORER detector name string */
#define LAL_EXPLORER_DETECTOR_PREFIX             	"E1"	/**< EXPLORER detector prefix string */
#define LAL_EXPLORER_DETECTOR_LONGITUDE_RAD      	0.10821041362	/**< EXPLORER vertex longitude (rad) */
#define LAL_EXPLORER_DETECTOR_LATITUDE_RAD       	0.81070543755	/**< EXPLORER vertex latitude (rad) */
#define LAL_EXPLORER_DETECTOR_ELEVATION_SI       	0	/**< EXPLORER vertex elevation (rad) */
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


/** \name Nautilus Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * Nautilus Resonant Mass Detector. */
/*@{*/
#define LAL_NAUTILUS_DETECTOR_NAME               	"Nautilus"	/**< Nautilus detector name string */
#define LAL_NAUTILUS_DETECTOR_PREFIX             	"N1"	/**< Nautilus detector prefix string */
#define LAL_NAUTILUS_DETECTOR_LONGITUDE_RAD      	0.22117684946	/**< Nautilus vertex longitude (rad) */
#define LAL_NAUTILUS_DETECTOR_LATITUDE_RAD       	0.72996456710	/**< Nautilus vertex latitude (rad) */
#define LAL_NAUTILUS_DETECTOR_ELEVATION_SI       	0	/**< Nautilus vertex elevation (rad) */
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


/** \name NIOBE Resonant Mass Detector constants
 * The following constants describe the location and geometry of the
 * NIOBE Resonant Mass Detector. */
/*@{*/
#define LAL_NIOBE_DETECTOR_NAME               	"NIOBE"	/**< NIOBE detector name string */
#define LAL_NIOBE_DETECTOR_PREFIX             	"B1"	/**< NIOBE detector prefix string */
#define LAL_NIOBE_DETECTOR_LONGITUDE_RAD      	2.02138216202	/**< NIOBE vertex longitude (rad) */
#define LAL_NIOBE_DETECTOR_LATITUDE_RAD       	-0.55734180780	/**< NIOBE vertex latitude (rad) */
#define LAL_NIOBE_DETECTOR_ELEVATION_SI       	0	/**< NIOBE vertex elevation (rad) */
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

#ifdef __cplusplus
}
#endif

#endif /* _LALDETECTORS_H */
