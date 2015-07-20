/*
 *  Copyright (C) 2014 Michael Puerrer, John Veitch
 *  Reduced Order Model for SEOBNR
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

/**
 * \author Michael Puerrer
 *
 * \file
 *
 * \brief Auxiliary functions for SEOBNRv1/v2 reduced order modeling codes in
 * LALSimIMRSEOBNRv2ROMDoubleSpin.c, LALSimIMRSEOBNRv1ROMDoubleSpin.c,
 * LALSimIMRSEOBNRv2ROMEffectiveSpin.c, LALSimIMRSEOBNRv1ROMEffectiveSpin.c,
 * LALSimIMRSEOBNRv2ChirpTime.c.
 *
 * Here I collect common auxiliary functions pertaining to reading
 * data stored in gsl binary vectors and matrices,
 * parameter space interpolation with B-splines,
 * a custom gsl error handler and adjustment of nearby parameter values.
 */

#include <stdio.h>
#include <lal/XLALError.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>


static void err_handler(const char *reason, const char *file, int line, int gsl_errno);
static int read_vector(const char dir[], const char fname[], gsl_vector *v);
UNUSED static int read_matrix(const char dir[], const char fname[], gsl_matrix *m);

UNUSED static REAL8 Interpolate_Coefficent_Tensor(
  gsl_vector *v,
  REAL8 eta,
  REAL8 chi1,
  REAL8 chi2,
  int ncy,
  int ncz,
  gsl_bspline_workspace *bwx,
  gsl_bspline_workspace *bwy,
  gsl_bspline_workspace *bwz
);

UNUSED static REAL8 Interpolate_Coefficent_Matrix(
  gsl_vector *v,
  REAL8 eta,
  REAL8 chi,
  int ncx,
  int ncy,
  gsl_bspline_workspace *bwx,
  gsl_bspline_workspace *bwy
);

static bool approximately_equal(REAL8 x, REAL8 y, REAL8 epsilon);
static void nudge(REAL8 *x, REAL8 X, REAL8 epsilon);


// Definitions

static void err_handler(const char *reason, const char *file, int line, int gsl_errno) {
  XLALPrintError("gsl: %s:%d: %s - %d\n", file, line, reason, gsl_errno);
}

// Helper functions to read gsl_vector and gsl_matrix data with error checking
static int read_vector(const char dir[], const char fname[], gsl_vector *v) {
  size_t size = strlen(dir) + strlen(fname) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, fname);
  
  FILE *f = fopen(path, "rb");
  if (!f) 
    XLAL_ERROR(XLAL_EIO, "Could not find ROM data file at path `%s'", path);
  int ret = gsl_vector_fread(f, v);
  if (ret != 0)
     XLAL_ERROR(XLAL_EIO, "Error reading data from `%s'", path);
  fclose(f);
  
  XLAL_PRINT_INFO("Sucessfully read data file `%s'", path);
  XLALFree(path);
  return(XLAL_SUCCESS);
}

static int read_matrix(const char dir[], const char fname[], gsl_matrix *m) {
  size_t size = strlen(dir) + strlen(fname) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, fname);

  FILE *f = fopen(path, "rb");
  if (!f) 
    XLAL_ERROR(XLAL_EIO, "Could not find ROM data file at path `%s'", path);
  int ret = gsl_matrix_fread(f, m);
  if (ret != 0)
     XLAL_ERROR(XLAL_EIO, "Error reading data from `%s'", path);
  fclose(f);
  
  XLAL_PRINT_INFO("Sucessfully read data file `%s'", path);
  XLALFree(path);
  return(XLAL_SUCCESS);
}

// Helper function to perform tensor product spline interpolation with gsl
// The gsl_vector v contains the ncx x ncy x ncz dimensional coefficient tensor in vector form
// that should be interpolated and evaluated at position (eta,chi1,chi2).
static REAL8 Interpolate_Coefficent_Tensor(
  gsl_vector *v,
  REAL8 eta,
  REAL8 chi1,
  REAL8 chi2,
  int ncy,
  int ncz,
  gsl_bspline_workspace *bwx,
  gsl_bspline_workspace *bwy,
  gsl_bspline_workspace *bwz
) {
  // Store nonzero cubic (order k=4) B-spline basis functions in the eta and chi directions.
  gsl_vector *Bx4 = gsl_vector_alloc(4);
  gsl_vector *By4 = gsl_vector_alloc(4);
  gsl_vector *Bz4 = gsl_vector_alloc(4);

  size_t isx, isy, isz; // first non-zero spline
  size_t iex, iey, iez; // last non-zero spline
  // Evaluate all potentially nonzero cubic B-spline basis functions for
  // positions (eta,chi) and store them in the vectors Bx4, By4, Bz4.
  // Since the B-splines are of compact support we only need to store a small
  // number of basis functions to avoid computing terms that would be zero anyway.
  // https://www.gnu.org/software/gsl/manual/html_node/Overview-of-B_002dsplines.html#Overview-of-B_002dsplines
  gsl_bspline_eval_nonzero(eta,  Bx4, &isx, &iex, bwx);
  gsl_bspline_eval_nonzero(chi1, By4, &isy, &iey, bwy);
  gsl_bspline_eval_nonzero(chi2, Bz4, &isz, &iez, bwz);

  // Now compute coefficient at desired parameters (q,chi1,chi2)
  // from C(eta,chi1,chi2) = c_ijk * Beta_i * Bchi1_j * Bchi2_k
  // while summing over indices i,j,k where the B-splines are nonzero.
  // Note: in the 2D case we were able to use gsl_matrix c = gsl_matrix_view_vector(&v, ncx, ncy).matrix
  // to convert vector view of the coefficient matrix to a matrix view.
  // However, since tensors are not supported in gsl, we have to do the indexing explicitly.
  double sum = 0;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      for (int k=0; k<4; k++) {
        int ii = isx + i;
        int jj = isy + j;
        int kk = isz + k;
        double cijk = gsl_vector_get(v, (ii*ncy + jj)*ncz + kk);
        sum += cijk * gsl_vector_get(Bx4, i) * gsl_vector_get(By4, j) * gsl_vector_get(Bz4, k);
      }

  gsl_vector_free(Bx4);
  gsl_vector_free(By4);
  gsl_vector_free(Bz4);

  return sum;
}

// Helper function to perform tensor product spline interpolation with gsl
// The gsl_vector v contains the ncx x ncy dimensional coefficient matrix in vector form
// that should be interpolated and evaluated at position (eta,chi).
static REAL8 Interpolate_Coefficent_Matrix(
  gsl_vector *v,
  REAL8 eta,
  REAL8 chi,
  int ncx,
  int ncy,
  gsl_bspline_workspace *bwx,
  gsl_bspline_workspace *bwy
) {
  gsl_matrix c = gsl_matrix_view_vector(v, ncx, ncy).matrix;   // Convert coefficient matrix from vector view to matrix view c_ij.

  // Store nonzero cubic (order k=4) B-spline basis functions in the eta and chi directions.
  gsl_vector *Bx4 = gsl_vector_alloc(4);
  gsl_vector *By4 = gsl_vector_alloc(4);

  REAL8 sum = 0;
  size_t isx, isy; // first non-zero spline
  size_t iex, iey; // last non-zero spline
  // Evaluate all potentially nonzero cubic B-spline basis functions for positions (eta,chi) and stores them in the vectors Bx4, By4.
  // Since the B-splines are of compact support we only need to store a small number of basis functions
  // to avoid computing terms that would be zero anyway.
  // https://www.gnu.org/software/gsl/manual/html_node/Overview-of-B_002dsplines.html#Overview-of-B_002dsplines
  gsl_bspline_eval_nonzero(eta, Bx4, &isx, &iex, bwx);
  gsl_bspline_eval_nonzero(chi, By4, &isy, &iey, bwy);

  // Now compute coefficient at desired parameters (eta,chi) from C(eta,chi) = c_ij * Beta_i * Bchi_j
  // summing over indices i,j where the B-splines are nonzero.
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      sum += gsl_matrix_get(&c, isx + i, isy + j) * gsl_vector_get(Bx4, i) * gsl_vector_get(By4, j);

  gsl_vector_free(Bx4);
  gsl_vector_free(By4);

  return sum;
}

// This function determines whether x and y are approximately equal to a relative accuracy epsilon.
// Note that x and y are compared to relative accuracy, so this function is not suitable for testing whether a value is approximately zero.
static bool approximately_equal(REAL8 x, REAL8 y, REAL8 epsilon) {
  return !gsl_fcmp(x, y, epsilon);
}

// If x and X are approximately equal to relative accuracy epsilon then set x = X.
// If X = 0 then use an absolute comparison.
static void nudge(REAL8 *x, REAL8 X, REAL8 epsilon) {
  if (X != 0.0) {
    if (approximately_equal(*x, X, epsilon)) {
      XLAL_PRINT_INFO("Nudging value %.15g to %.15g\n", *x, X);
      *x = X;
    }
  }
  else {
    if (fabs(*x - X) < epsilon)
      *x = X;
  }
}
