//
// Copyright (C) 2014 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#ifndef _SUPERSKYMETRICS_H
#define _SUPERSKYMETRICS_H

#include <gsl/gsl_matrix.h>
#include <lal/LALStdlib.h>
#include <lal/UniversalDopplerMetric.h>
#include <lal/LatticeTiling.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// \file
/// \author Karl Wette
/// \brief Functions for using the super-sky and reduced super-sky metrics of \cite WettePrix2013a
///

///
/// Coordinate systems associated with the super-sky metrics
///
typedef enum {
  SSC_PHYSICAL,					///< Physical: right ascension, declination, frequency and spindowns
  SSC_SUPER_SKY,				///< Super-sky: 3-dimensional sky, spindowns and frequency
  SSC_REDUCED_SUPER_SKY,			///< Reduced super-sky: 2-dimensional sky, reduced spindowns and frequency
  SSC_MAX
} SuperSkyCoordinates;

///
/// Compute the expanded super-sky metric, which separates spin and orbital sky components.
///
int XLALExpandedSuperSkyMetric(
  gsl_matrix **essky_metric,			///< [out] Pointer to allocated expanded super-sky metric
  const size_t spindowns,			///< [in] Number of frequency spindown coordinates
  const LIGOTimeGPS* ref_time,			///< [in] Reference time for the metric
  const LALSegList* segments,			///< [in] List of segments to average metric over
  const double fiducial_freq,			///< [in] Fiducial frequency for sky-position coordinates
  const MultiLALDetector* detectors,		///< [in] List of detector to average metric over
  const MultiNoiseFloor* detector_weights,	///< [in] Weights used to combine single-detector metrics (default: unit weights)
  const DetectorMotionType detector_motion,	///< [in] Which detector motion to use
  const EphemerisData* ephemerides		///< [in] Earth/Sun ephemerides
  );

///
/// Compute the (untransformed) super-sky metric in equatorial coordinates from the expanded super-sky metric.
///
int XLALSuperSkyMetric(
  gsl_matrix **ssky_metric,			///< [out] Pointer to allocated super-sky metric
  const gsl_matrix* essky_metric		///< [in] Input expanded super-sky metric
  );

///
/// Compute the reduced super-sky metric and coordinate transform data from the expanded super-sky metric.
///
int XLALReducedSuperSkyMetric(
  gsl_matrix **rssky_metric,			///< [out] Pointer to allocated reduced super-sky metric
  gsl_matrix **rssky_transf,			///< [out] Pointer to allocated coordinate transform data
  const gsl_matrix* essky_metric		///< [in] Input expanded super-sky metric
  );

///
/// Convert a series of points between super-sky coordinate systems
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(gsl_matrix**, out_points));
#endif
int XLALConvertSuperSkyCoordinates(
  const SuperSkyCoordinates out,		///< [in] Coordinate system of the output points
  gsl_matrix** out_points,			///< [in/out] Matrix whose columns are the output points
  const SuperSkyCoordinates in,			///< [in] Coordinate system of the input points
  const gsl_matrix* in_points,			///> [in] Matrix whose columns are the input points
  const gsl_matrix* rssky_transf		///< [in] Reduced super-sky coordinate transform data
  );

///
/// Convert a single point from physical to super-sky coordinates
///
int XLALConvertPhysicalToSuperSky(
  const SuperSkyCoordinates out,		///< [in] Coordinate system of the output point
  gsl_vector* out_point,			///< [in/out] Output point in super-sky coordinates
  const PulsarDopplerParams* in_phys,		///< [in] Input point in physical coordinates
  const gsl_matrix* rssky_transf		///< [in] Reduced super-sky coordinate transform data
  );

///
/// Convert a single point from super-sky to physical coordinates
///
int XLALConvertSuperSkyToPhysical(
  PulsarDopplerParams* out_phys,		///< [in/out] Output point in physical coordinates
  const SuperSkyCoordinates in,			///< [in] Coordinate system of the input point
  const gsl_vector* in_point,			///< [in] Input point in super-sky coordinates
  const gsl_matrix* rssky_transf		///< [in] Reduced super-sky coordinate transform data
  );

///
/// Set all-sky parameter-space bounds on a lattice tiling using the reduced super-sky metric
///
int XLALSetLatticeReducedSuperSkyBounds(
  LatticeTiling* tiling				///< [in] Tiling state
  );

///
/// Set a sky point parameter-space bound on a lattice tiling using the reduced super-sky metric
///
int XLALSetLatticeReducedSuperSkyPointBounds(
  LatticeTiling* tiling,			///< [in] Tiling state
  const gsl_matrix* rssky_transf,		///< [in] Reduced super-sky coordinate transform data
  const double alpha,				///< [in] Sky point right ascension
  const double delta				///< [in] Sky point declination
  );

///
/// Set lattice tiling parameter-space bounds on the (untransformed) super-sky frequency/spindowns \f$f^{(s)}\f$
///
int XLALSetLatticeSuperSkySpinBound(
  LatticeTiling* tiling,			///< [in] Tiling state
  const gsl_matrix* rssky_transf,		///< [in] Reduced super-sky coordinate transform data
  const size_t s,				///< [in] Spindown order; 0=frequency, 1=first spindown, etc.
  const double bound1,				///< [in] First bound on frequency/spindown
  const double bound2,				///< [in] Second bound on frequency/spindown
  bool extra_padding				///< [in] Use extra padding to cover sky position mismatch
  );

///
/// Set lattice tiling parameter-space bounds on the reduced super-sky frequency/spindowns \f$\nu^{(s)}\f$,
/// which are related to the super-sky frequency/spindowns by \f$\nu^{(s)} = f^{(s)} + \vec\Delta^s \cdot \vec n\f$
///
int XLALSetLatticeReducedSuperSkySpinBound(
  LatticeTiling* tiling,			///< [in] Tiling state
  const gsl_matrix* rssky_transf,		///< [in] Reduced super-sky coordinate transform data
  const size_t s,				///< [in] Spindown order; 0=frequency, 1=first spindown, etc.
  const double bound1,				///< [in] First bound on frequency/spindown
  const double bound2				///< [in] Second bound on frequency/spindown
  );

#ifdef __cplusplus
}
#endif

#endif
