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
#include <stdarg.h> 
#include <lal/XLALError.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_fit.h>
#include <LALSimBlackHoleRingdown.h>

#ifdef LAL_HDF5_ENABLED
#include <lal/H5FileIO.h>
#endif

UNUSED static int read_vector(const char dir[], const char fname[], gsl_vector *v);
UNUSED static int read_matrix(const char dir[], const char fname[], gsl_matrix *m);
/* SEOBNRv4HM_ROM functions */
UNUSED static char* concatenate_strings(int count, ...);
UNUSED static REAL8 sigmoid(REAL8 x);
UNUSED static UINT4 blend(gsl_vector * freqs, gsl_vector * out_fun, REAL8 freq_1, REAL8 freq_2);
UNUSED static UINT4 blend_functions(gsl_vector * freqs_out, gsl_vector * out_fun, gsl_vector * freq_in_1, gsl_vector * fun_in_1, gsl_vector * freq_in_2, gsl_vector * fun_in_2, REAL8 freq_1, REAL8 freq_2);
UNUSED static UINT4 compute_i_max_LF_i_min_HF(INT8 * i_max_LF, INT8 * i_min_LF,gsl_vector * freqs_in_1,gsl_vector * freqs_in_2,REAL8 freq_1);
UNUSED static REAL8 Get_omegaQNM_SEOBNRv4(REAL8 q, REAL8 chi1z, REAL8 chi2z, UINT4 l, UINT4 m);
UNUSED static UINT4 unwrap_phase(gsl_vector* phaseout,gsl_vector* phasein);
UNUSED static UINT8 compute_i_at_f(gsl_vector * freq_array, REAL8 freq);
UNUSED static UINT4 align_wfs_window(gsl_vector* f_array_1, gsl_vector* f_array_2, gsl_vector* phase_1, gsl_vector* phase_2, REAL8 f_align_start, REAL8 f_align_end);

#ifdef LAL_HDF5_ENABLED
UNUSED static int CheckVectorFromHDF5(LALH5File *file, const char name[], const double *v, size_t n);
UNUSED static int ReadHDF5RealVectorDataset(LALH5File *file, const char *name, gsl_vector **data);
UNUSED static int ReadHDF5RealMatrixDataset(LALH5File *file, const char *name, gsl_matrix **data);
UNUSED static int ReadHDF5LongVectorDataset(LALH5File *file, const char *name, gsl_vector_long **data);
UNUSED static int ReadHDF5LongMatrixDataset(LALH5File *file, const char *name, gsl_matrix_long **data);
UNUSED static void PrintInfoStringAttribute(LALH5File *file, const char attribute[]);
UNUSED static int ROM_check_version_number(LALH5File *file, INT4 version_major_in, INT4 version_minor_in, INT4 version_micro_in);
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

static int ReadHDF5LongVectorDataset(LALH5File *file, const char *name, gsl_vector_long **data) {
	LALH5Dataset *dset;
	UINT4Vector *dimLength;
	size_t n;

	if (file == NULL || name == NULL || data == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	dset = XLALH5DatasetRead(file, name);
	if (dset == NULL)
		XLAL_ERROR(XLAL_EFUNC);

	if (XLALH5DatasetQueryType(dset) != LAL_I8_TYPE_CODE) {
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
		*data = gsl_vector_long_alloc(n);
		if (*data == NULL) {
			XLALH5DatasetFree(dset);
			XLAL_ERROR(XLAL_ENOMEM, "gsl_vector_long_alloc(%zu) failed", n);
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

static int ReadHDF5LongMatrixDataset(LALH5File *file, const char *name, gsl_matrix_long **data) {
	LALH5Dataset *dset;
	UINT4Vector *dimLength;
	size_t n1, n2;

	if (file == NULL || name == NULL || data == NULL)
		XLAL_ERROR(XLAL_EFAULT);

	dset = XLALH5DatasetRead(file, name);
	if (dset == NULL)
		XLAL_ERROR(XLAL_EFUNC);

	if (XLALH5DatasetQueryType(dset) != LAL_I8_TYPE_CODE) {
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
		*data = gsl_matrix_long_alloc(n1, n2);
		if (*data == NULL) {
			XLALH5DatasetFree(dset);
			XLAL_ERROR(XLAL_ENOMEM, "gsl_matrix_long_alloc(%zu, %zu) failed", n1, n2);
		}
	}
	else if ((*data)->size1 != n1 || (*data)->size2 != n2) {
		XLALH5DatasetFree(dset);
		XLAL_ERROR(XLAL_EINVAL, "Expected gsl_matrix_long `%s' of size %zu x %zu", name, n1, n2);
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
  LALH5Generic gfile = {.file = file};
  int len = XLALH5AttributeQueryStringValue(NULL, 0, gfile, attribute) + 1;
  char *str = XLALMalloc(len);
  XLALH5FileQueryStringAttributeValue(str, len, file, attribute);
  XLALPrintInfo("%s:\n%s\n", attribute, str);
  LALFree(str);
}

static int ROM_check_version_number(LALH5File *file, 	INT4 version_major_in, INT4 version_minor_in, INT4 version_micro_in) {
  INT4 version_major;
  INT4 version_minor;
  INT4 version_micro;

  LALH5Generic gfile = {.file = file};
  XLALH5AttributeQueryScalarValue(&version_major, gfile, "version_major");
  XLALH5AttributeQueryScalarValue(&version_minor, gfile, "version_minor");
  XLALH5AttributeQueryScalarValue(&version_micro, gfile, "version_micro");

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

/* SEOBNRv4HM_ROM functions */


// Helper function to concatenate strings
char* concatenate_strings(int count, ...)
{
    va_list ap;
    int i;

    // Find required length to store merged string
    int len = 1; // room for NULL
    va_start(ap, count);
    for(i=0 ; i<count ; i++)
        len += strlen(va_arg(ap, char*));
    va_end(ap);

    // Allocate memory to concat strings
    char *merged = XLALCalloc(sizeof(char),len);
    int null_pos = 0;

    // Actually concatenate strings
    va_start(ap, count);
    for(i=0 ; i<count ; i++)
    {
        char *s = va_arg(ap, char*);
        strcpy(merged+null_pos, s);
        null_pos += strlen(s);
    }
    va_end(ap);

    return merged;
}

/* Sigmoid function for hybridization */
REAL8 sigmoid(REAL8 x)
{
     REAL8 exp_value;
     REAL8 return_value;

     /*** Exponential calculation ***/
     exp_value = exp(-x);

     /*** Final sigmoid value ***/
     return_value = 1 / (1 + exp_value);

     return return_value;
}

// Compute activation function in a window [freq_1, freq_2]
UINT4 blend(gsl_vector * freqs, gsl_vector * out_fun, REAL8 freq_1, REAL8 freq_2){
  for(unsigned int i = 0; i < freqs->size; i++){
    if(freqs->data[i] <= freq_1){
      out_fun->data[i] = 0.;
    }
    else if((freqs->data[i]>freq_1) && (freqs->data[i]<freq_2)){
      out_fun->data[i] = sigmoid(- (freq_2 - freq_1)/(freqs->data[i]-freq_1) - (freq_2 - freq_1)/(freqs->data[i]-freq_2));
    }
    else{
      out_fun->data[i] = 1.;
    }
  }
  return XLAL_SUCCESS;
}

// Function to compute the indices associated to [freq_1] in two arrays
UINT4 compute_i_max_LF_i_min_HF(INT8 * i_max_LF, INT8 * i_min_HF,gsl_vector * freqs_in_1,gsl_vector * freqs_in_2, REAL8 freq_1){

  for(unsigned int i = 0; i < freqs_in_1->size; i++){
    if(freqs_in_1->data[i]<freq_1){
      *i_max_LF = i;
    }
  }
  for(unsigned int i = freqs_in_2->size -1; i > 0; i--){
    if(freqs_in_2->data[i]>=freq_1){
      *i_min_HF = i;
    }
  }
return XLAL_SUCCESS;
}

// Function to compute the index associated to [freq_1] in an array
UINT8 compute_i_at_f(gsl_vector * freq_array, REAL8 freq){
  UINT8 index = 0;
  REAL8 diff = fabs(freq-freq_array->data[index]);
  for(unsigned int i = 1; i < freq_array->size; i++){
    if(fabs(freq-freq_array->data[i])< diff){
      diff = fabs(freq-freq_array->data[i]);
      index = i;
    }
  }
  return index;
}

/* Function to blend LF and HF ROM functions in a window [freq_1, freq_2] */
UINT4 blend_functions(gsl_vector * freqs_out, gsl_vector * out_fun, gsl_vector * freq_in_1, gsl_vector * fun_in_1, gsl_vector * freq_in_2, gsl_vector * fun_in_2, REAL8 freq_1, REAL8 freq_2){

  gsl_vector *bld_2 = gsl_vector_calloc(out_fun->size);

  /* Generate the window function */
  UNUSED UINT4 ret = blend(freqs_out, bld_2,freq_1,freq_2);

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, freq_in_1->size);
  gsl_spline_init (spline, freq_in_1->data, fun_in_1->data, freq_in_1->size);

  unsigned int i = 0;
  REAL8 inter_func;
  /* Build the contribution to out_fun from fun_in_1 in the range [f_ini, freq_2] */
  /* where f_ini is the starting frequency of the array freq_in_1 */
  /* for freq > freq_2 the contribution of fun_in_1 to out_fun is null since 1 - bld_2 = 0 there*/
  /* for freq < freq_1, bld_2 = 0, so I could have another easier while in the range [freq_1,freq_2]*/

  while(freqs_out->data[i]<=freq_2){
    /* Interpolate the original function on the new grid */
    inter_func = gsl_spline_eval (spline, freqs_out->data[i], acc);
    /* Multiply the original function by the window function */
    out_fun->data[i] = inter_func*(1.-bld_2->data[i]);
    i++;
  }
  gsl_spline_free (spline);

  /* Build the contribution to out_fun from fun_in_2 in the range [freq_2, f_end] */
  /* where f_end is the ending frequency of the array freq_in_2 */
  /* for freq > freq_2, bld_2 = 1*/
  spline = gsl_spline_alloc (gsl_interp_cspline, freq_in_2->size);
  gsl_spline_init (spline, freq_in_2->data, fun_in_2->data, freq_in_2->size);
  i = freqs_out->size-1;
  while(freqs_out->data[i]>freq_2){
    inter_func = gsl_spline_eval (spline, freqs_out->data[i], acc);
    out_fun->data[i] = inter_func;
    i--;
  }
  
  /* Build the contribution to out_fun from fun_in_2 in the range [freq_1, freq_2] */
  /* where f_end is the ending frequency of the array freq_in_2 */
  while((freqs_out->data[i]>=freq_1)&&(freqs_out->data[i]<=freq_2)){
    inter_func = gsl_spline_eval (spline, freqs_out->data[i], acc);
    out_fun->data[i] += inter_func*bld_2->data[i];
    /* Note that the += is important here */
    i--;
  }

  /* Cleanup */
  gsl_vector_free(bld_2);
  gsl_spline_free (spline);
  gsl_interp_accel_free(acc);

  return XLAL_SUCCESS;
}

// Function to align two phases in a window [f_align_start,f_align_end]
UINT4 align_wfs_window(gsl_vector* f_array_1, gsl_vector* f_array_2, gsl_vector* phase_1, gsl_vector* phase_2, REAL8 f_align_start, REAL8 f_align_end){

    // Search the indices corresponding to [f_align_start,f_align_end] in f_array_2
    INT8 i_align_start = compute_i_at_f(f_array_2, f_align_start);
    INT8 i_align_end = compute_i_at_f(f_array_2, f_align_end);
    // Interpolate phase_1 phase on phase_2 frequency grid between [f_align_start,f_align_end] and compute phase difference
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, phase_1->size);
    gsl_spline_init (spline, f_array_1->data, phase_1->data, phase_1->size);
    gsl_vector* delta_phi = gsl_vector_alloc(i_align_end-i_align_start);
    gsl_vector* freq_delta_phi = gsl_vector_alloc(i_align_end-i_align_start);
    for(unsigned int i = 0; i <freq_delta_phi->size; i++){
      freq_delta_phi->data[i] = f_array_2->data[i+i_align_start];
      delta_phi->data[i] = gsl_spline_eval(spline, f_array_2->data[i+i_align_start], acc) - phase_2->data[i+i_align_start];
    }
    REAL8 c0, c1, cov00, cov01, cov11, chisq;
    gsl_fit_linear (freq_delta_phi->data, 1, delta_phi->data, 1, freq_delta_phi->size, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
    
    // Remove the linear fit from phase_1
    for(unsigned int i = 0; i < f_array_1->size; i++){
      phase_1->data[i] = phase_1->data[i] - (c1*f_array_1->data[i] + c0);
    }
    /* Cleanup */
    gsl_vector_free(delta_phi);
    gsl_vector_free(freq_delta_phi);
    gsl_spline_free (spline);
    gsl_interp_accel_free(acc);

    return XLAL_SUCCESS;

}

// Function to compute the QNM frequency
REAL8 Get_omegaQNM_SEOBNRv4(REAL8 q, REAL8 chi1z, REAL8 chi2z, UINT4 l, UINT4 m){
    // Total mass M is not important here, we return the QNM frequencies in geometric units
    REAL8 M = 100.; 
    REAL8 Ms = M * LAL_MTSUN_SI;
    REAL8 m1 = M * q/(1+q);
    REAL8 m2 = M * 1/(1+q);
    Approximant SpinAlignedEOBapproximant = SEOBNRv4;
    COMPLEX16Vector modefreqVec;
    COMPLEX16 modeFreq;
    modefreqVec.length = 1;
    modefreqVec.data = &modeFreq;
    REAL8 spin1[3] = {0., 0., chi1z};
    REAL8 spin2[3] = {0., 0., chi2z}; 
    // XLALSimIMREOBGenerateQNMFreqV2 is returning QNM frequencies in SI units, we multiply by Ms to convert it in geometric units 
    UNUSED UINT4 ret = XLALSimIMREOBGenerateQNMFreqV2(&modefreqVec, m1, m2, spin1, spin2, l, m, 1, SpinAlignedEOBapproximant);
    return Ms * creal(modefreqVec.data[0]);
}

/* Function to unwrap phases mod 2pi  - acts on a REAL8Vector representing the phase */
/* FIXME: for long vectors there are small differences with the numpy unwrap function - to be checked */
static UINT4 unwrap_phase(
  gsl_vector* phaseout,    /* Output: unwrapped phase vector, already allocated - can be the same as phasein, in which case unwrapping is done in place */
  gsl_vector* phasein)     /* Input: phase vector */
{
  int N = phasein->size;
  double* p = phasein->data;
  double* pmod = (double*) malloc(sizeof(double) * N);
  int* jumps = (int*) malloc(sizeof(int) * N);
  int* cumul = (int*) malloc(sizeof(int) * N);

  /* Compute phase mod 2pi (shifted to be between -pi and pi) */
  for(int i=0; i<N; i++) {
    pmod[i] = p[i] - floor((p[i] + LAL_PI) / (2*LAL_PI))*(2*LAL_PI);
  }

  /* Identify jumps */
  jumps[0] = 0;
  double d = 0.;
  for(int i=1; i<N; i++) {
    d = pmod[i] - pmod[i-1];
    if(d<-LAL_PI) jumps[i] = 1;
    else if(d>LAL_PI) jumps[i] = -1;
    else jumps[i] = 0;
  }

  /* Cumulative of the jump sequence */
  int c = 0;
  cumul[0] = 0;
  for(int i=1; i<N; i++) {
    c += jumps[i];
    cumul[i] = c;
  }

  /* Correct mod 2pi phase series by the number of 2pi factor given by the cumulative of the jumps */
  double* pout = phaseout->data;
  for(int i=0; i<N; i++) {
    pout[i] = pmod[i] + 2*LAL_PI*cumul[i];
  }

  /* Cleanup */
  free(pmod);
  free(jumps);
  free(cumul);

	return XLAL_SUCCESS;
}