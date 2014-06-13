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

#include <stdio.h>
#include <inttypes.h>
#include <math.h>

#include <lal/LatticeTiling.h>
#include <lal/LALStdlib.h>
#include <lal/XLALError.h>

static int CheckLatticeTiling(
  const size_t n,
  const LatticeType lattice,
  const size_t total_ref
  )
{

  // Create a lattice tiling
  fprintf(stderr, "Number of dimensions: %zu\n", n);
  LatticeTiling* tiling = XLALCreateLatticeTiling(n);
  XLAL_CHECK(tiling != NULL, XLAL_ENOMEM);

  // Add bounds
  for (size_t i = 0; i < n; ++i) {
    XLAL_CHECK(XLALSetLatticeConstantBound(tiling, i, 0.0, pow(100.0, 1.0/n)) == XLAL_SUCCESS, XLAL_EFUNC);
  }

  // Set metric to the Lehmer matrix
  gsl_matrix* metric = gsl_matrix_alloc(n, n);
  XLAL_CHECK(metric != NULL, XLAL_ENOMEM);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      const double ii = i+1, jj = j+1;
      gsl_matrix_set(metric, i, j, jj >= ii ? ii/jj : jj/ii);
    }
  }

  // Set lattice and metric
  fprintf(stderr, "Lattice type: %u\n", lattice);
  XLAL_CHECK(XLALSetLatticeTypeAndMetric(tiling, lattice, metric, 0.3) == XLAL_SUCCESS, XLAL_EFUNC);

  // Count number of templates
  size_t total = XLALCountLatticePoints(tiling);
  fprintf(stderr, "Number of lattice points: %zu\n", total);
  XLAL_CHECK(total == total_ref, XLAL_EFUNC, "ERROR: total = %zu != %zu = total_ref", total, total_ref);

  // Get all templates
  gsl_matrix* templates = gsl_matrix_alloc(n, total);
  XLAL_CHECK(templates != NULL, XLAL_ENOMEM);
  for (size_t i = 0; i < total; ++i) {
    gsl_vector_view point = gsl_matrix_column(templates, i);
    XLALNextLatticePoint(tiling, &point.vector);
    XLAL_CHECK(XLALGetLatticePointCount(tiling) == i + 1, XLAL_EFAILED);
  }
  XLAL_CHECK(XLALNextLatticePoint(tiling, NULL) < 0, XLAL_EFAILED);
  XLAL_CHECK(XLALRestartLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);

  // Get nearest point to each template; should be template itself
  gsl_matrix* nearest = gsl_matrix_alloc(n, total);
  XLAL_CHECK(nearest != NULL, XLAL_ENOMEM);
  UINT8Vector* indices = XLALCreateUINT8Vector(total);
  XLAL_CHECK(indices != NULL, XLAL_ENOMEM);
  XLAL_CHECK(XLALBuildLatticeIndexLookup(tiling) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK(XLALNearestLatticePoints(tiling, templates, nearest, indices) == XLAL_SUCCESS, XLAL_EFUNC);
  size_t failed = 0;
  for (size_t i = 0; i < total; ++i) {
    if (indices->data[i] != i) {
      ++failed;
      fprintf(stderr, "warning: indices->data[i] = %"PRIu64" != %zu\n", indices->data[i], i);
    }
  }
  if (failed >= 5) {
    XLAL_ERROR(XLAL_EFAILED, "ERROR: number of failed index lookups = %zu >= 10", failed);
  } else if (failed > 0) {
    fprintf(stderr, "warning: number of failed index lookups = %zu > 0", failed);
  }

  // Cleanup
  XLALDestroyLatticeTiling(tiling);
  gsl_matrix_free(metric);
  gsl_matrix_free(templates);
  XLALDestroyUINT8Vector(indices);
  LALCheckMemoryLeaks();
  fprintf(stderr, "\n");

  return XLAL_SUCCESS;

}

int main(void) {

  // Check tiling in 1 dimension
  XLAL_CHECK(CheckLatticeTiling(1, LATTICE_TYPE_CUBIC,    93) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK(CheckLatticeTiling(1, LATTICE_TYPE_ANSTAR,   93) == XLAL_SUCCESS, XLAL_EFAILED);

  // Check tiling in 2 dimensions
  XLAL_CHECK(CheckLatticeTiling(2, LATTICE_TYPE_CUBIC,   189) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK(CheckLatticeTiling(2, LATTICE_TYPE_ANSTAR,  121) == XLAL_SUCCESS, XLAL_EFAILED);

  // Check tiling in 3 dimensions
  XLAL_CHECK(CheckLatticeTiling(3, LATTICE_TYPE_CUBIC,   623) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK(CheckLatticeTiling(3, LATTICE_TYPE_ANSTAR,  295) == XLAL_SUCCESS, XLAL_EFAILED);

  // Check tiling in 4 dimensions
  XLAL_CHECK(CheckLatticeTiling(4, LATTICE_TYPE_CUBIC,  2417) == XLAL_SUCCESS, XLAL_EFAILED);
  XLAL_CHECK(CheckLatticeTiling(4, LATTICE_TYPE_ANSTAR, 1031) == XLAL_SUCCESS, XLAL_EFAILED);

  return 0;

}
