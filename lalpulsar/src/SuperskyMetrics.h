//
// Copyright (C) 2014, 2015 Karl Wette
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
/// \defgroup SuperskyMetrics_h Header SuperskyMetrics.h
/// \ingroup lalpulsar_metric
/// \author Karl Wette
/// \brief Compute the supersky metrics and coordinate transforms of \cite WettePrix2013a .
///
/// @{
///

///
/// Coordinate systems associated with the supersky metrics.
///
typedef enum {
  SC_PHYS,					///< Physical coordinates, in order: right ascension, declination, frequency, spindowns
  SC_USSKY,					///< Unrestricted supersky (equatorial) coordinates, in order: 3-dimensional sky, frequency, spindowns
  SC_RSSKY,					///< Reduced supersky coordinates, in order: 2-dimensional sky, reduced spindowns, <b>frequency last</b>
  SC_MAX
} SuperskyCoordinates;

///
/// Computed supersky metrics, returned by XLALComputeSuperskyMetrics().
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(ARRAY_MULTIPLE_LENGTHS(tagSuperskyMetrics, num_segments));
#endif // SWIG
typedef struct tagSuperskyMetrics {
  size_t num_segments;				///< Number of segments
  double fiducial_freq;				///< Fiducial frequency at which metrics were calculated

#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(SuperskyMetrics, gsl_matrix*, ussky_metric_seg, size_t, num_segments));
#endif // SWIG
  gsl_matrix **ussky_metric_seg;		///< Unrestricted supersky metric (3-dimensional sky), for each segment
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(SuperskyMetrics, gsl_matrix*, rssky_metric_seg, size_t, num_segments));
#endif // SWIG
  gsl_matrix **rssky_metric_seg;		///< Reduced supersky metric (2-dimensional sky), for each segment
#ifdef SWIG // SWIG interface directives
  SWIGLAL(ARRAY_1D(SuperskyMetrics, gsl_matrix*, rssky_transf_seg, size_t, num_segments));
#endif // SWIG
  gsl_matrix **rssky_transf_seg;		///< Coordinate transform data of reduced supersky metric, for each segment

  gsl_matrix *ussky_metric_avg;			///< Unrestricted supersky metric (3-dimensional sky), averaged over segments
  gsl_matrix *rssky_metric_avg;			///< Reduced supersky metric (2-dimensional sky), averaged over segments
  gsl_matrix *rssky_transf_avg;			///< Coordinate transform data of reduced supersky metric, averaged over segments

} SuperskyMetrics;

///
/// Compute the various supersky metrics, which are returned in a #SuperskyMetrics struct.
///
SuperskyMetrics *XLALComputeSuperskyMetrics(
  const size_t spindowns,			///< [in] Number of frequency+spindown coordinates
  const LIGOTimeGPS *ref_time,			///< [in] Reference time for the metrics
  const LALSegList *segments,			///< [in] List of segments to average metrics over
  const double fiducial_freq,			///< [in] Fiducial frequency for sky-position coordinates
  const MultiLALDetector *detectors,		///< [in] List of detectors to average metrics over
  const MultiNoiseFloor *detector_weights,	///< [in] Weights used to combine single-detector metrics (default: unit weights)
  const DetectorMotionType detector_motion,	///< [in] Which detector motion to use
  const EphemerisData *ephemerides		///< [in] Earth/Sun ephemerides
  );

///
/// Destroy a #SuperskyMetrics struct.
///
void XLALDestroySuperskyMetrics(
  SuperskyMetrics *metrics			///< [in] Supersky metrics struct
  );

///
/// Scale all supersky metrics to match the given fiducial frequency.
///
int XLALScaleSuperskyMetricsFiducialFreq(
  SuperskyMetrics *metrics,			///< [in] Supersky metrics struct
  const double fiducial_freq			///< [in] Fiducial frequency for sky-position coordinates
  );

///
/// Project and rescale the averaged reduced supersky metrics in the frequency dimension, such that
/// all reduced supersky metrics have the same frequency spacing for the given maximum mismatches.
///
int XLALEqualizeReducedSuperskyMetricsFreqSpacing(
  SuperskyMetrics *metrics,			///< [in] Supersky metrics struct
  const double coh_max_mismatch,		///< [in] Maximum coherent mismatch, for per-segment metrics
  const double semi_max_mismatch		///< [in] Maximum semicoherent mismatch, for averaged metrics
  );

///
/// Convert a series of points between supersky coordinate systems.
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(gsl_matrix **, out_points));
#endif
int XLALConvertSuperskyCoordinates(
  const SuperskyCoordinates out,		///< [in] Coordinate system of the output points
  gsl_matrix **out_points,			///< [in/out] Matrix whose columns are the output points
  const SuperskyCoordinates in,			///< [in] Coordinate system of the input points
  const gsl_matrix *in_points,			///< [in] Matrix whose columns are the input points
  const gsl_matrix *rssky_transf		///< [in] Reduced supersky coordinate transform data
  );

///
/// Convert a single point from physical to supersky coordinates.
///
int XLALConvertPhysicalToSupersky(
  const SuperskyCoordinates out,		///< [in] Coordinate system of the output point
  gsl_vector *out_point,			///< [in/out] Output point in supersky coordinates
  const PulsarDopplerParams *in_phys,		///< [in] Input point in physical coordinates
  const gsl_matrix *rssky_transf,		///< [in] Reduced supersky coordinate transform data
  const LIGOTimeGPS *ref_time			///< [in] Reference time for the coordinate transform data
  );

///
/// Convert a single point from supersky to physical coordinates.
///
int XLALConvertSuperskyToPhysical(
  PulsarDopplerParams *out_phys,		///< [in/out] Output point in physical coordinates
  const SuperskyCoordinates in,			///< [in] Coordinate system of the input point
  const gsl_vector *in_point,			///< [in] Input point in supersky coordinates
  const gsl_matrix *rssky_transf,		///< [in] Reduced supersky coordinate transform data
  const LIGOTimeGPS *ref_time			///< [in] Reference time for the coordinate transform data
  );

#ifdef SWIG // SWIG interface directives
SWIGLAL(COPYINOUT_ARRAYS(gsl_matrix, rssky_metric, rssky_transf));
#endif // SWIG

///
/// Set parameter-space bounds on the physical sky position \f$(\alpha, \delta)\f$ for a lattice
/// tiling using the reduced supersky metric. The metric and coordinate transform data must be supplied,
/// since they will be transformed such that the physical sky region maps to a region in the reduced
/// supersky coordinates \f$(n_a,n_b)\f$ which may be covered by the lattice tiling.
///
int XLALSetSuperskyLatticeTilingPhysicalSkyBounds(
  LatticeTiling *tiling,			///< [in] Lattice tiling
  gsl_matrix *rssky_metric,			///< [in] Reduced supersky metric
  gsl_matrix *rssky_transf,			///< [in] Coordinate transform data of reduced supersky metric
  const double alpha1,				///< [in] First bound on sky position right ascension
  const double alpha2,				///< [in] Second bound on sky position right ascension
  const double delta1,				///< [in] First bound on sky position declination
  const double delta2				///< [in] Second bound on sky position declination
  );

///
/// Set parameter-space bounds on an equal-area fraction of the physical sky \f$(\alpha, \delta)\f$
/// for a lattice tiling using the reduced supersky metric. The metric and coordinate transform data
/// must be supplied, since they will be transformed such that the physical sky patch maps to a region
/// in the reduced supersky coordinates \f$(n_a,n_b)\f$ which may be covered by the lattice tiling.
///
int XLALSetSuperskyLatticeTilingPhysicalSkyPatch(
  LatticeTiling *tiling,			///< [in] Lattice tiling
  gsl_matrix *rssky_metric,			///< [in] Reduced supersky metric
  gsl_matrix *rssky_transf,			///< [in] Coordinate transform data of reduced supersky metric
  const UINT4 patch_count,			///< [in] Number of equal-area patches to divide sky into
  const UINT4 patch_index			///< [in] Index of the patch for which to set bounds
  );

#ifdef SWIG // SWIG interface directives
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(gsl_matrix, rssky_metric, rssky_transf));
#endif // SWIG

///
/// Set parameter-space bounds on the physical frequency/spindowns \f$f^{(s)}\f$ for a lattice
/// tiling using the reduced supersky metric.
///
int XLALSetSuperskyLatticeTilingPhysicalSpinBound(
  LatticeTiling *tiling,			///< [in] Lattice tiling
  const gsl_matrix *rssky_transf,		///< [in] Reduced supersky coordinate transform data
  const size_t s,				///< [in] Spindown order; 0=frequency, 1=first spindown, etc.
  const double bound1,				///< [in] First bound on frequency/spindown
  const double bound2				///< [in] Second bound on frequency/spindown
  );

///
/// Set parameter-space bounds on the reduced supersky frequency/spindown coordinates \f$\nu^{(s)}\f$
/// for a lattice tiling using the reduced supersky metric. These coordinates are related to the
/// physical frequency/spindowns by \f$\nu^{(s)} = f^{(s)} + \vec\Delta^s \cdot \vec n\f$.
///
int XLALSetSuperskyLatticeTilingCoordinateSpinBound(
  LatticeTiling *tiling,			///< [in] Lattice tiling.
  const gsl_matrix *rssky_transf,		///< [in] Reduced supersky coordinate transform data
  const size_t s,				///< [in] Spindown order; 0=frequency, 1=first spindown, etc.
  const double bound1,				///< [in] First bound on frequency/spindown
  const double bound2				///< [in] Second bound on frequency/spindown
  );

///
/// Fill a PulsarSpinRange with the physical frequency/spindown ranges covered by a reduced supersky
/// lattice tiling.
///
int XLALSuperskyLatticePulsarSpinRange(
  PulsarSpinRange *spin_range,			///< [in,out] Physical frequency/spindown range
  LatticeTiling *tiling,			///< [in] Lattice tiling
  const gsl_matrix *rssky_transf,		///< [in] Reduced supersky coordinate transform data
  const LIGOTimeGPS *ref_time			///< [in] Reference time for the coordinate transform data
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _SUPERSKYMETRICS_H
