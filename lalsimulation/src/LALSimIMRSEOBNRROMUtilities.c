/*
 *  Copyright (C) 2014, 2015 Michael Puerrer, John Veitch
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
 * LALSimIMRSEOBNRv2ChirpTime.c, LALSimIMRSEOBNRv2ROMDoubleSpinHI.c.
 *
 * Here I collect common auxiliary functions pertaining to reading
 * data stored in gsl binary vectors and matrices, HDF5 files using the LAL interface,
 * parameter space interpolation with B-splines, fitting to a cubic,
 * a custom gsl error handler and adjustment of nearby parameter values.
 */

#include <stdio.h>
#include <lal/XLALError.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>

#ifdef LAL_HDF5_ENABLED
#include <lal/H5FileIO.h>
#endif

UNUSED static void err_handler(const char *reason, const char *file, int line, int gsl_errno);
UNUSED static int read_vector(const char dir[], const char fname[], gsl_vector *v);
UNUSED static int read_matrix(const char dir[], const char fname[], gsl_matrix *m);

#ifdef LAL_HDF5_ENABLED
UNUSED static int CheckVectorFromHDF5(LALH5File *file, const char name[], const double *v, size_t n);
UNUSED static int ReadHDF5RealVectorDataset(LALH5File *file, const char *name, gsl_vector **data);
UNUSED static int ReadHDF5RealMatrixDataset(LALH5File *file, const char *name, gsl_matrix **data);
UNUSED static void PrintInfoStringAttribute(LALH5File *file, const char attribute[]);
UNUSED static int ROM_check_version_number(LALH5File *file, 	INT4 version_major_in, INT4 version_minor_in, INT4 version_micro_in);
#endif

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

UNUSED static gsl_vector *Fit_cubic(const gsl_vector *xi, const gsl_vector *yi);

UNUSED static bool approximately_equal(REAL8 x, REAL8 y, REAL8 epsilon);
UNUSED static void nudge(REAL8 *x, REAL8 X, REAL8 epsilon);

UNUSED static double SEOBNRROM_Ringdown_Mf_From_Mtot_q(
  const double Mtot_sec,
  const double q,
  const double chi1,
  const double chi2,
  Approximant apx
);
UNUSED static double SEOBNRROM_Ringdown_Mf_From_Mtot_Eta(
  const double Mtot_sec,
  const double eta,
  const double chi1,
  const double chi2,
  Approximant apx
);


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

#ifdef LAL_HDF5_ENABLED
static int CheckVectorFromHDF5(LALH5File *file, const char name[], const double *v, size_t n) {
  gsl_vector *temp = NULL;
  ReadHDF5RealVectorDataset(file, name, &temp);

  if (temp->size != n)
    XLAL_ERROR(XLAL_EIO, "Length of data %s disagrees", name);

  for (size_t i=0; i<temp->size; i++)
    if (gsl_vector_get(temp, i) != v[i])
      XLAL_ERROR(XLAL_EIO, "Data %s disagrees", name);
  gsl_vector_free(temp);
  return XLAL_SUCCESS;
}

static int ReadHDF5RealVectorDataset(LALH5File *file, const char *name, gsl_vector **data) {
	LALH5Dataset *dset;
	UINT4Vector *dimLength;
	size_t n;

	if (file == NULL || name == NULL || data == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	dset = XLALH5DatasetRead(file, name);
	if (dset == NULL)
		XLAL_ERROR(XLAL_EFUNC);

	if (XLALH5DatasetQueryType(dset) != LAL_D_TYPE_CODE) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_ETYPE, "Dataset `%s' is wrong type", name);
	}

	dimLength = XLALH5DatasetQueryDims(dset);
	if (dimLength == NULL) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC);
	}
	if (dimLength->length != 1) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EDIMS, "Dataset `%s' must be 1-dimensional", name);
	}

	n = dimLength->data[0];
	XLALDestroyUINT4Vector(dimLength);

	if (*data == NULL) {
		*data = gsl_vector_alloc(n);
		if (*data == NULL) {
			XLALH5DatasetFree(dset);
			XLAL_ERROR(XLAL_ENOMEM, "gsl_vector_alloc(%zu) failed", n);
		}
	}
	else if ((*data)->size != n) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EINVAL, "Expected gsl_vector `%s' of size %zu", name, n);
	}

  // Now read the data
	if (XLALH5DatasetQueryData((*data)->data, dset) < 0) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC);
	}

	XLALH5DatasetFree(dset);
	return 0;
}

static int ReadHDF5RealMatrixDataset(LALH5File *file, const char *name, gsl_matrix **data) {
	LALH5Dataset *dset;
	UINT4Vector *dimLength;
	size_t n1, n2;

	if (file == NULL || name == NULL || data == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	dset = XLALH5DatasetRead(file, name);
	if (dset == NULL)
		XLAL_ERROR(XLAL_EFUNC);

	if (XLALH5DatasetQueryType(dset) != LAL_D_TYPE_CODE) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_ETYPE, "Dataset `%s' is wrong type", name);
	}

	dimLength = XLALH5DatasetQueryDims(dset);
	if (dimLength == NULL) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC);
	}
	if (dimLength->length != 2) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EDIMS, "Dataset `%s' must be 2-dimensional", name);
	}

	n1 = dimLength->data[0];
	n2 = dimLength->data[1];
	XLALDestroyUINT4Vector(dimLength);

	if (*data == NULL) {
		*data = gsl_matrix_alloc(n1, n2);
		if (*data == NULL) {
			XLALH5DatasetFree(dset);
			XLAL_ERROR(XLAL_ENOMEM, "gsl_matrix_alloc(%zu, %zu) failed", n1, n2);
		}
	}
	else if ((*data)->size1 != n1 || (*data)->size2 != n2) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EINVAL, "Expected gsl_matrix `%s' of size %zu x %zu", name, n1, n2);
	}

  // Now read the data
	if (XLALH5DatasetQueryData((*data)->data, dset) < 0) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EFUNC);
	}

	XLALH5DatasetFree(dset);
	return 0;
}

static void PrintInfoStringAttribute(LALH5File *file, const char attribute[]) {
  int len = XLALH5FileQueryStringAttributeValue(NULL, 0, file, attribute) + 1;
  char *str = XLALMalloc(len);
  XLALH5FileQueryStringAttributeValue(str, len, file, attribute);
  XLALPrintInfo("%s\n", str);
  LALFree(str);
}

static int ROM_check_version_number(LALH5File *file, 	INT4 version_major_in, INT4 version_minor_in, INT4 version_micro_in) {
  INT4 version_major;
  INT4 version_minor;
  INT4 version_micro;

  XLALH5FileQueryScalarAttributeValue(&version_major, file, "version_major");
  XLALH5FileQueryScalarAttributeValue(&version_minor, file, "version_minor");
  XLALH5FileQueryScalarAttributeValue(&version_micro, file, "version_micro");

  if ((version_major_in != version_major) || (version_minor_in != version_minor) || (version_micro_in != version_micro)) {
    XLAL_ERROR(XLAL_EIO, "Expected ROM data version %d.%d.%d, but got version %d.%d.%d.",
    version_major_in, version_minor_in, version_micro_in, version_major, version_minor, version_micro);
  }
  else {
    XLALPrintInfo("Reading ROM data version %d.%d.%d.\n", version_major, version_minor, version_micro);
    return XLAL_SUCCESS;
  }
}
#endif

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

// Returns fitting coefficients for cubic y = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3
static gsl_vector *Fit_cubic(const gsl_vector *xi, const gsl_vector *yi) {
  const int n = xi->size; // how many data points are we fitting
  const int p = 4;        // order of fitting polynomial

  gsl_matrix *X = gsl_matrix_alloc(n, p); // design matrix
  gsl_vector *y = gsl_vector_alloc(n);    // data vector
  gsl_vector *c = gsl_vector_alloc(p);    // coefficient vector
  gsl_matrix *cov = gsl_matrix_alloc(p, p);
  double chisq;

  for (int i=0; i < n; i++) {
      double xval = gsl_vector_get(xi, i);
      double xval2 = xval*xval;
      double xval3 = xval2*xval;
      gsl_matrix_set(X, i, 0, 1.0);
      gsl_matrix_set(X, i, 1, xval);
      gsl_matrix_set(X, i, 2, xval2);
      gsl_matrix_set(X, i, 3, xval3);

      double yval = gsl_vector_get(yi, i);
      gsl_vector_set (y, i, yval);
    }

  gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, p);
  gsl_multifit_linear(X, y, c, cov, &chisq, work); // perform linear least-squares fit

  // clean up
  gsl_multifit_linear_free(work);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);

  return c;
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

/* compute ringdown frequency for (2,2) mode from total mass and mass ratio */
static double SEOBNRROM_Ringdown_Mf_From_Mtot_q(
  const double Mtot_sec,
  const double q,
  const double chi1,
  const double chi2,
  Approximant apx
)
{
  double Mtot_SI = Mtot_sec / LAL_MTSUN_SI * LAL_MSUN_SI;
  double m1_SI = Mtot_SI * q / (1. + q);
  double m2_SI = Mtot_SI * 1. / (1. + q);
  return XLALSimInspiralGetFinalFreq(m1_SI, m2_SI, 0, 0, chi1, 0, 0, chi2, apx) * Mtot_sec;
}

/* compute ringdown frequency for (2,2) mode from total mass and eta */
static double SEOBNRROM_Ringdown_Mf_From_Mtot_Eta(
  const double Mtot_sec,
  const double eta,
  const double chi1,
  const double chi2,
  Approximant apx
)
{
  double q = (1. + sqrt(1. - 4. * eta) - 2. * eta) / (2. * eta);
  return SEOBNRROM_Ringdown_Mf_From_Mtot_q(Mtot_sec, q, chi1, chi2, apx);
}
