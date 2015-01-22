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
/// \defgroup SuperSkyMetrics_h Header SuperSkyMetrics.h
/// \ingroup lalpulsar_metric
/// \author Karl Wette
/// \brief Functions which compute the super-sky and reduced super-sky metrics
/// of \cite WettePrix2013a .
///
/// @{
///

///
/// Coordinate systems associated with the super-sky metrics.
///
typedef enum {
  /// Physical: right ascension, declination, frequency and spindowns.
  SSC_PHYSICAL,
  /// Super-sky: 3-dimensional sky, spindowns and frequency.
  SSC_SUPER_SKY,
  /// Reduced super-sky: 2-dimensional sky, reduced spindowns and frequency.
  SSC_REDUCED_SUPER_SKY,
  /// \cond DONT_DOXYGEN
  SSC_MAX
  /// \endcond
} SuperSkyCoordinates;

///
/// Compute the expanded super-sky metric, which separates spin and orbital sky
/// components.
///
int XLALExpandedSuperSkyMetric(
  /// [out] Pointer to allocated expanded super-sky metric.
  gsl_matrix **essky_metric,
  /// [in] Number of frequency spindown coordinates.
  const size_t spindowns,
  /// [in] Reference time for the metric.
  const LIGOTimeGPS* ref_time,
  /// [in] List of segments to average metric over.
  const LALSegList* segments,
  /// [in] Fiducial frequency for sky-position coordinates.
  const double fiducial_freq,
  /// [in] List of detector to average metric over.
  const MultiLALDetector* detectors,
  /// [in] Weights used to combine single-detector metrics (default: unit
  /// weights).
  const MultiNoiseFloor* detector_weights,
  /// [in] Which detector motion to use.
  const DetectorMotionType detector_motion,
  /// [in] Earth/Sun ephemerides.
  const EphemerisData* ephemerides
  );

///
/// Compute the (untransformed) super-sky metric in equatorial coordinates from
/// the expanded super-sky metric.
///
int XLALSuperSkyMetric(
  /// [out] Pointer to allocated super-sky metric.
  gsl_matrix **ssky_metric,
  /// [in] Input expanded super-sky metric.
  const gsl_matrix* essky_metric
  );

///
/// Compute the reduced super-sky metric and coordinate transform data from the
/// expanded super-sky metric.
///
int XLALReducedSuperSkyMetric(
  /// [out] Pointer to allocated reduced super-sky metric.
  gsl_matrix **rssky_metric,
  /// [out] Pointer to allocated coordinate transform data.
  gsl_matrix **rssky_transf,
  /// [in] Input expanded super-sky metric.
  const gsl_matrix* essky_metric
  );

///
/// Convert a series of points between super-sky coordinate systems.
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(gsl_matrix**, out_points));
#endif
int XLALConvertSuperSkyCoordinates(
  /// [in] Coordinate system of the output points.
  const SuperSkyCoordinates out,
  /// [in/out] Matrix whose columns are the output points.
  gsl_matrix** out_points,
  /// [in] Coordinate system of the input points.
  const SuperSkyCoordinates in,
  /// [in] Matrix whose columns are the input points.
  const gsl_matrix* in_points,
  /// [in] Reduced super-sky coordinate transform data.
  const gsl_matrix* rssky_transf
  );

///
/// Convert a single point from physical to super-sky coordinates.
///
int XLALConvertPhysicalToSuperSky(
  /// [in] Coordinate system of the output point.
  const SuperSkyCoordinates out,
  /// [in/out] Output point in super-sky coordinates.
  gsl_vector* out_point,
  /// [in] Input point in physical coordinates.
  const PulsarDopplerParams* in_phys,
  /// [in] Reduced super-sky coordinate transform data.
  const gsl_matrix* rssky_transf
  );

///
/// Convert a single point from super-sky to physical coordinates.
///
int XLALConvertSuperSkyToPhysical(
  /// [in/out] Output point in physical coordinates.
  PulsarDopplerParams* out_phys,
  /// [in] Coordinate system of the input point.
  const SuperSkyCoordinates in,
  /// [in] Input point in super-sky coordinates.
  const gsl_vector* in_point,
  /// [in] Reduced super-sky coordinate transform data.
  const gsl_matrix* rssky_transf
  );

///
/// Set all-sky parameter-space bounds on a lattice tiling using the reduced
/// super-sky metric.
///
int XLALSetLatticeTilingReducedSuperSkyBounds(
  /// [in] Lattice tiling parameter space.
  LatticeTilingSpace* space
  );

///
/// Set a sky point parameter-space bound on a lattice tiling using the reduced
/// super-sky metric.
///
int XLALSetLatticeTilingReducedSuperSkyPointBounds(
  /// [in] Lattice tiling parameter space.
  LatticeTilingSpace* space,
  /// [in] Reduced super-sky coordinate transform data.
  const gsl_matrix* rssky_transf,
  /// [in] Sky point right ascension.
  const double alpha,
  /// [in] Sky point declination.
  const double delta
  );

///
/// Set lattice tiling parameter-space bounds on the physical
/// frequency/spindowns \f$f^{(s)}\f$.
///
int XLALSetLatticeTilingPhysicalSpinBound(
  /// [in] Lattice tiling parameter space.
  LatticeTilingSpace* space,
  /// [in] Reduced super-sky coordinate transform data.
  const gsl_matrix* rssky_transf,
  /// [in] Spindown order; 0=frequency, 1=first spindown, etc.
  const size_t s,
  /// [in] First bound on frequency/spindown.
  const double bound1,
  /// [in] Second bound on frequency/spindown.
  const double bound2
  );

///
/// Set lattice tiling parameter-space bounds on the reduced super-sky
/// frequency/spindowns \f$\nu^{(s)}\f$, which are related to the super-sky
/// frequency/spindowns by \f$\nu^{(s)} = f^{(s)} + \vec\Delta^s \cdot \vec
/// n\f$.
///
int XLALSetLatticeTilingReducedSuperSkySpinBound(
  /// [in] Lattice tiling parameter space.
  LatticeTilingSpace* space,
  /// [in] Reduced super-sky coordinate transform data.
  const gsl_matrix* rssky_transf,
  /// [in] Spindown order; 0=frequency, 1=first spindown, etc.
  const size_t s,
  /// [in] First bound on frequency/spindown.
  const double bound1,
  /// [in] Second bound on frequency/spindown.
  const double bound2
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif
