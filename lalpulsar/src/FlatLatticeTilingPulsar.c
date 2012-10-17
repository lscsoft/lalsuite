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

#include <lal/FlatLatticeTilingPulsar.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

typedef struct {
  size_t freq_dim;
  double lower;
  double upper;
} F1DotAgeBrakingBoundInfo;

static void F1DotAgeBrakingBound(
  const gsl_vector_uint* bound UNUSED,
  const gsl_vector* point,
  const void* data,
  gsl_vector* lower,
  gsl_vector* upper
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
  const gsl_vector_uint* bound UNUSED,
  const gsl_vector* point,
  const void* data,
  gsl_vector* lower,
  gsl_vector* upper
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
