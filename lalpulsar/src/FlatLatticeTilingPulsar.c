//
// Copyright (C) 2007, 2008, 2012 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

#include <math.h>
#include <gsl/gsl_math.h>

#include <lal/FlatLatticeTilingPulsar.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static void SuperSkyNZBound(
  const size_t dimension,
  const gsl_vector_uint* bound UNUSED,
  const gsl_vector* point,
  const void* data,
  const gsl_vector* incr UNUSED,
  const gsl_vector* bbox UNUSED,
  gsl_vector* lower,
  gsl_vector* upper UNUSED,
  double* lower_pad UNUSED,
  double* upper_pad UNUSED
  )
{

  // Get bounds data
  const FLSSNZ type = *((const FLSSNZ*)data);

  // Calculate nz from nx and ny
  const double nx = gsl_vector_get(point, dimension - 2);
  const double ny = gsl_vector_get(point, dimension - 1);
  const double nxysqr = nx * nx + ny * ny;
  const double nz = (nxysqr < 1.0) ? sqrt(1.0 - nxysqr) : 0.0;

  // Set singular bounds on nz
  if (nz > 0.0) {
    switch (type) {
    case FLSSNZ_LOWER:
      gsl_vector_set(lower, 0, -nz);
      break;
    case FLSSNZ_PLANE:
      gsl_vector_set(lower, 0, 0.0);
      break;
    case FLSSNZ_UPPER:
      gsl_vector_set(lower, 0, nz);
      break;
    case FLSSNZ_SPHERE:
      gsl_vector_set(lower, 0, -nz);
      gsl_vector_set(lower, 1, nz);
      break;
    default:
      return;
    }
  } else {
    gsl_vector_set(lower, 0, 0.0);
  }

};

int XLALSetFlatLatticeSuperSkyNZBound(
  FlatLatticeTiling* tiling,
  const size_t nz_dimension,
  const FLSSNZ type
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(type < FLSSNZ_LAST, XLAL_EINVAL);

  // Allocate and set bounds data
  FLSSNZ* info = XLALCalloc(1, sizeof(FLSSNZ));
  XLAL_CHECK(info != NULL, XLAL_ENOMEM);
  *info = type;

  // Set parameter space bound
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, nz_dimension, true, SuperSkyNZBound, (void*)info) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

typedef struct {
  size_t nx_dim;
  double offset[3];
  double bounds[2];
} FnDotConstantBoundInfo;

static void FnDotConstantBound(
  const size_t dimension UNUSED,
  const gsl_vector_uint* bound UNUSED,
  const gsl_vector* point UNUSED,
  const void* data,
  const gsl_vector* incr UNUSED,
  const gsl_vector* bbox UNUSED,
  gsl_vector* lower,
  gsl_vector* upper,
  double* lower_pad UNUSED,
  double* upper_pad UNUSED
  )
{

  // Get bounds data
  const FnDotConstantBoundInfo* info = (const FnDotConstantBoundInfo*)data;

  // Calculate sky offset
  double offset = 0.0;
  for (size_t i = 0; i < 3; ++i) {
    offset += info->offset[i] * gsl_vector_get(point, info->nx_dim + i);
  }

  // Set constant lower and upper bounds on offset frequency/spindowns
  gsl_vector_set(lower, 0, info->bounds[0] + offset);
  if (upper) {
    gsl_vector_set(upper, 0, info->bounds[1] + offset);
  }

}

int XLALSetFlatLatticeFnDotConstantBound(
  FlatLatticeTiling* tiling,
  const size_t nx_dimension,
  const double offset[3],
  size_t dimension,
  double bound1,
  double bound2
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(nx_dimension + 3 <= dimension, XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound1), XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound2), XLAL_EINVAL);

  // Allocate and set bounds data
  FnDotConstantBoundInfo* info = XLALCalloc(1, sizeof(FnDotConstantBoundInfo));
  XLAL_CHECK(info != NULL, XLAL_ENOMEM);
  info->nx_dim = nx_dimension;
  for (size_t i = 0; i < 3; ++i) {
    if (offset) {
      XLAL_CHECK(isfinite(info->offset[i]), XLAL_EINVAL);
      info->offset[i] = offset[i];
    } else {
      info->offset[i] = 0;
    }
  }
  info->bounds[0] = GSL_MIN(bound1, bound2);
  info->bounds[1] = GSL_MAX(bound1, bound2);

  // Set parameter space bound
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, dimension, bound1 == bound2, FnDotConstantBound, (void*)info) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

typedef struct {
  size_t freq_dim;
  double lower;
  double upper;
} F1DotAgeBrakingBoundInfo;

static void F1DotAgeBrakingBound(
  const size_t dimension UNUSED,
  const gsl_vector_uint* bound UNUSED,
  const gsl_vector* point,
  const void* data,
  const gsl_vector* incr UNUSED,
  const gsl_vector* bbox UNUSED,
  gsl_vector* lower,
  gsl_vector* upper,
  double* lower_pad UNUSED,
  double* upper_pad UNUSED
  )
{

  // Get bounds info
  const F1DotAgeBrakingBoundInfo* info = (const F1DotAgeBrakingBoundInfo*)data;

  // Get current value of frequency
  const double freq = gsl_vector_get(point, info->freq_dim);

  // Set first spindown bounds
  gsl_vector_set(lower, 0, info->lower * freq);
  if (upper != NULL) {
    gsl_vector_set(upper, 0, info->upper * freq);
  }

}

int XLALSetFlatLatticeF1DotAgeBrakingBound(
  FlatLatticeTiling* tiling,
  const size_t freq_dimension,
  const size_t f1dot_dimension,
  const double age,
  const double min_braking,
  const double max_braking
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(freq_dimension < f1dot_dimension, XLAL_EINVAL);
  XLAL_CHECK(age > 0.0, XLAL_EINVAL);
  XLAL_CHECK(min_braking > 1.0, XLAL_EINVAL);
  XLAL_CHECK(max_braking > 1.0, XLAL_EINVAL);
  XLAL_CHECK(min_braking <= max_braking, XLAL_EINVAL);

  // Allocate memory
  F1DotAgeBrakingBoundInfo* info = XLALCalloc(1, sizeof(F1DotAgeBrakingBoundInfo));
  XLAL_CHECK(info != NULL, XLAL_ENOMEM);

  // Set bounds info
  info->freq_dim = freq_dimension;
  info->lower = -1.0 / ((min_braking - 1.0) * age);
  info->upper = -1.0 / ((max_braking - 1.0) * age);

  // Set parameter space bound
  const bool is_singular = (min_braking == max_braking);
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, f1dot_dimension, is_singular, F1DotAgeBrakingBound, (void*)info) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

typedef struct {
  size_t freq_dim;
  size_t f1dot_dim;
  double lower;
  double upper;
} F2DotBrakingBoundInfo;

static void F2DotBrakingBound(
  const size_t dimension UNUSED,
  const gsl_vector_uint* bound UNUSED,
  const gsl_vector* point,
  const void* data,
  const gsl_vector* incr UNUSED,
  const gsl_vector* bbox UNUSED,
  gsl_vector* lower,
  gsl_vector* upper,
  double* lower_pad UNUSED,
  double* upper_pad UNUSED
  )
{

  // Get bounds info
  const F2DotBrakingBoundInfo* info = (const F2DotBrakingBoundInfo*)data;

  // Get current values of frequency and first spindown
  const double freq = gsl_vector_get(point, info->freq_dim);
  const double f1dot = gsl_vector_get(point, info->f1dot_dim);

  // Set first spindown bounds
  const double x = f1dot * f1dot / freq;
  gsl_vector_set(lower, 0, info->lower * x);
  if (upper != NULL) {
    gsl_vector_set(upper, 0, info->upper * x);
  }

}

int XLALSetFlatLatticeF2DotBrakingBound(
  FlatLatticeTiling* tiling,
  const size_t freq_dimension,
  const size_t f1dot_dimension,
  const size_t f2dot_dimension,
  const double min_braking,
  const double max_braking
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(freq_dimension < f1dot_dimension, XLAL_EINVAL);
  XLAL_CHECK(f1dot_dimension < f2dot_dimension, XLAL_EINVAL);
  XLAL_CHECK(min_braking > 0.0, XLAL_EINVAL);
  XLAL_CHECK(max_braking > 0.0, XLAL_EINVAL);
  XLAL_CHECK(min_braking <= max_braking, XLAL_EINVAL);

  // Allocate memory
  F2DotBrakingBoundInfo* info = XLALCalloc(1, sizeof(F2DotBrakingBoundInfo));
  XLAL_CHECK(info != NULL, XLAL_ENOMEM);

  // Set bounds info
  info->freq_dim = freq_dimension;
  info->f1dot_dim = f1dot_dimension;
  info->lower = min_braking;
  info->upper = max_braking;

  // Set parameter space bound
  const bool is_singular = (min_braking == max_braking);
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, f2dot_dimension, is_singular, F2DotBrakingBound, (void*)info) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}
