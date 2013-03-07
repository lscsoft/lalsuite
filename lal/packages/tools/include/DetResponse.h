/*
*  Copyright (C) 2007 David Chin, Jolien Creighton, Kipp Cannon, Teviet Creighton
*  Copyright (C) 2012 Matthew Pitkin
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

#ifndef _DETRESPONSE_H
#define _DETRESPONSE_H

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <lal/Date.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** \addtogroup DetResponse_h
    \author David Chin <dwchin@umich.edu>, Kipp Cannon <kipp@gravity.phys.uwm.edu>

    \brief Provides routines to compute gravitational wave detector response to
    polarized planar gravitational wave originating from a given source,
    detected at a given time.

\heading{Synopsis}

\code
#include <lal/DetResponse.h>
\endcode

\heading{Description}

These routines compute the antenna beam pattern for all supported detector
types.  <tt>XLALComputeDetAMResponse()</tt> computes the response at one
instance in time, and <tt>XLALComputeDetAMResponseSeries()</tt> computes a
vector of response for some length of time.

\heading{Algorithm}

This code is a translation of the algorithm in the Maple worksheet by
Anderson, <em>et al.</em> [\ref Anderson_2000].  We compute the \f$h\f$-tensors for
\f$+\f$- and \f$\times\f$-polarized in the Earth-fixed frame, and then contract
them (take the scalar product) with the detector response tensors as
described in the \ref LALDetectors_h section of the \c tools package.

\ref LALDetectors_h provides predefined
\c LALDetector structures representing most current detectors,
including LIGO (Hanford and Livingston), and GEO.

\heading{Notes}

For examples of usage, please see the test programs in the \c test directory.

*/
/*@{*/

/** \name Error Codes */
/*@{*/
#define DETRESPONSEH_ENULLINPUT  1		/**< Input is NULL */
#define DETRESPONSEH_ENULLOUTPUT 2		/**< Output is NULL */
#define DETRESPONSEH_ESRCNOTEQUATORIAL 3	/**< Source coordinates not in Equatorial system */
/*@}*/

/** \cond DONT_DOXYGEN */
#define DETRESPONSEH_MSGENULLINPUT "Input is NULL"
#define DETRESPONSEH_MSGENULLOUTPUT "Output is NULL"
#define DETRESPONSEH_MSGESRCNOTEQUATORIAL "Source coordinates not in Equatorial system"
/** \endcond */


/** This structure contains gravitational wave source position (in Equatorial
 * coördinates), and orientation angle.
 * The orientation is measured counter-clockwise with respect to the "line of ascending nodes",
 * i.e. counter-clockwise with respect to a line perpendicular to the
 * source's meridian and extending westwards.  For a source in the Northern
 * celestial hemisphere, and an observer in the Northern hemisphere standing
 * such that they are facing South, this angle is measured counter-clockwise
 * from a 3 o'clock position (pointing West) at the source.  The polarization
 * convention is chosen such that if the source orientation were zero, the
 * source would be a pure \f$+\f$-polarized source.
 */
typedef struct
tagLALSource
{
  CHAR         name[LALNameLength];  /**< name of source, eg catalog number */
  SkyPosition  equatorialCoords;     /**< equatorial coordinates of source, in decimal RADIANS */
  REAL8        orientation;          /**< Orientation angle (\f$\psi\f$) of source:
                                      * counter-clockwise angle \f$x\f$-axis makes with a line perpendicular to
                                      * meridian of source in Westward direction (i.e. North of West),
                                      * in decimal radians.
                                      */
}
LALSource;

/** This structure aggregates a pointer to a \c LALDetector and a
 * \c LALSource.  Its sole function is to allow the user to pass
 * detector and source parameters to the functions
 * LALComputeDetAMResponse() and
 * LALComputeDetAMResponseSeries().
 */
typedef struct
tagLALDetAndSource
{
  const LALDetector  *pDetector;/**< Pointer to ::LALDetector object containing information about the detector */
  LALSource    *pSource;	/**< Pointer to ::LALSource object containing information about the source */
}
LALDetAndSource;

/** This structure encapsulates the detector AM (beam pattern) coefficients for
 * one source at one instance in time.
 */
typedef struct
tagLALDetAMResponse
{
  REAL4 plus;	/**< Detector response to \f$+\f$-polarized gravitational radiation  */
  REAL4 cross;	/**< Detector response to \f$\times\f$-polarized gravitational radiation */
  REAL4 scalar;	/**< Detector response to scalar gravitational radiation (NB: ignored at present -- scalar response computation not yet implemented) */
}
LALDetAMResponse;

/** This structure aggregates together three ::REAL4TimeSeries objects containing
 * time series of detector AM response.
 */
typedef struct
tagLALDetAMResponseSeries
{
  REAL4TimeSeries *pPlus;	/**< timeseries of detector response to \f$+\f$-polarized gravitational radiation */
  REAL4TimeSeries *pCross;	/**< timeseries of detector response to \f$\times\f$-polarized gravitational radiation */
  REAL4TimeSeries *pScalar;	/**< timeseries of detector response to scalar gravitational radiation (NB: not yet implemented.) */
}
LALDetAMResponseSeries;


/** This structure encapsulates time and sampling information for computing a
 * ::LALDetAMResponseSeries. Its fields correspond to some fields of the
 * TimeSeries structures for easy conversion.
 */
typedef struct
tagLALTimeIntervalAndNSample
{
  LIGOTimeGPS     epoch;	/**< The start time \f$t_0\f$ of the time series */
  REAL8           deltaT;	/**< The sampling interval \f$\Delta t\f$, in seconds */
  UINT4           nSample;	/**< The total number of samples to be computed */
}
LALTimeIntervalAndNSample;

/*
 * Function prototypes
 */

void
LALComputeDetAMResponse( LALStatus             *status,
                         LALDetAMResponse      *pResponse,
                         const LALDetAndSource *pDetAndSrc,
                         const LIGOTimeGPS     *gps);

void XLALComputeDetAMResponse(
	double *fplus,
	double *fcross,
	const REAL4 D[3][3],
	const double ra,
	const double dec,
	const double psi,
	const double gmst
);


void XLALComputeDetAMResponseExtraModes(
  double *fplus,
  double *fcross,
  double *fb,
  double *fl,
  double *fx,
  double *fy,
  REAL4 D[3][3],
  const double ra,
  const double dec,
  const double psi,
  const double gmst
);

/*
 * Gives a time series of the detector's response to plus and cross
 * polarization
 */
void
LALComputeDetAMResponseSeries( LALStatus                      *status,
                               LALDetAMResponseSeries         *pResponseSeries,
                               const LALDetAndSource          *pDetAndSource,
                               const LALTimeIntervalAndNSample *pTimeInfo);

int XLALComputeDetAMResponseSeries(
	REAL4TimeSeries **fplus,
	REAL4TimeSeries **fcross,
	const REAL4 D[3][3],
	const double ra,
	const double dec,
	const double psi,
	const LIGOTimeGPS *start,
	const double deltaT,
	const int n
);

int XLALComputeDetAMResponseExtraModesSeries(
  REAL4TimeSeries **fplus,
  REAL4TimeSeries **fcross,
  REAL4TimeSeries **fb,
  REAL4TimeSeries **fl,
  REAL4TimeSeries **fx,
  REAL4TimeSeries **fy,
  REAL4 D[3][3],
  const double ra,
  const double dec,
  const double psi,
  const LIGOTimeGPS *start,
  const double deltaT,
  const int n  
);

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* !defined _DETRESPONSE_H */
