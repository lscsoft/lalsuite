/*
 *  LALInferenceCreateROQ.c: Reduced order quadrature basis and interpolant generation
 *
 *  Copyright (C) 2014, 2016 Matthew Pitkin, Rory Smith
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

#include <lal/LALInferenceGenerateROQ.h>

#ifndef _OPENMP
#define omp ignore
#endif


/* internal function definitions */

/* define function to project model vector onto the training set of models */
void project_onto_basis(gsl_vector *weight,
                        gsl_matrix *RB,
                        gsl_matrix *TS,
                        gsl_matrix *projections,
                        INT4 idx);

void complex_project_onto_basis(gsl_vector *weight,
                                gsl_matrix_complex *RB,
                                gsl_matrix_complex *TS,
                                gsl_matrix_complex *projections,
                                INT4 idx);

/* the dot product of two real vectors scaled by a weighting factor */
REAL8 weighted_dot_product(gsl_vector *weight, gsl_vector *a, gsl_vector *b);

/* the dot product of two complex vectors scaled by a weighting factor */
gsl_complex complex_weighted_dot_product(gsl_vector *weight, gsl_vector_complex *a, gsl_vector_complex *b);

void normalise(gsl_vector *weight, gsl_vector *a);
void complex_normalise(gsl_vector *weight, gsl_vector_complex *a);

void normalise_training_set(gsl_vector *weight, gsl_matrix *TS);
void complex_normalise_training_set(gsl_vector *weight, gsl_matrix_complex *TS);

/* get the B_matrix */
REAL8Array *B_matrix(gsl_matrix *V, gsl_matrix *RB);
COMPLEX16Array *complex_B_matrix(gsl_matrix_complex *V, gsl_matrix_complex *RB);

/* find the index of the absolute maximum value for a complex vector */
int complex_vector_maxabs_index( gsl_vector_complex *c );


/** \brief Function to project the training set onto a given basis vector
 *
 * This is an internal function to be used by \c LALInferenceCreateREAL8OrthonormalBasis
 *
 * @param[in] weight The normalisation weight(s) for the training set waveforms (e.g. time or frequency step(s) size)
 * @param[in] RB The reduced basis set
 * @param[in] TS The training set of waveforms
 * @param[in] projections The set of projections (this is updated in this function)
 * @param[in] idx The index of the reduced basis vector that the training set will project onto
 */
void project_onto_basis(gsl_vector *weight,
                        gsl_matrix *RB,
                        gsl_matrix *TS,
                        gsl_matrix *projections,
                        INT4 idx){
  size_t row = 0;
  gsl_vector_view basis;

  XLAL_CALLGSL( basis = gsl_matrix_row(RB, idx) );

  for ( row=0; row < TS->size1; row++ ){
    double prod;
    gsl_vector_view proj, model;
    gsl_vector *basisscale;

    XLAL_CALLGSL( proj = gsl_matrix_row(projections, row) );
    XLAL_CALLGSL( basisscale = gsl_vector_calloc(TS->size2) );

    XLAL_CALLGSL( model = gsl_matrix_row(TS, row) ); /* get model from training set */

    prod = weighted_dot_product(weight, &basis.vector, &model.vector);

    XLAL_CALLGSL( gsl_vector_memcpy(basisscale, &basis.vector) );
    XLAL_CALLGSL( gsl_vector_scale(basisscale, prod) );
    XLAL_CALLGSL( gsl_vector_add(&proj.vector, basisscale) );
    XLAL_CALLGSL( gsl_vector_free(basisscale) );
  }
}


/** \brief Function to project the complex training set onto a given basis vector
 *
 * This is an internal function to be used by \c LALInferenceCreateCOMPLEX16OrthonormalBasis
 *
 * @param[in] weight The normalisation weight(s) for the training set waveforms (e.g. time or frequency step(s) size)
 * @param[in] RB The reduced basis set
 * @param[in] TS The training set of waveforms
 * @param[in] projections The set of projections (this is updated in this function)
 * @param[in] idx The index of the reduced basis vector that the training set will project onto
 */
void complex_project_onto_basis(gsl_vector *weight,
                                gsl_matrix_complex *RB,
                                gsl_matrix_complex *TS,
                                gsl_matrix_complex *projections,
                                INT4 idx){
  size_t row = 0;
  gsl_vector_complex_view basis;

  XLAL_CALLGSL( basis = gsl_matrix_complex_row(RB, idx) );

  for ( row=0; row < TS->size1; row++ ){
    gsl_complex cprod;
    gsl_vector_complex_view proj, model;
    gsl_vector_complex *basisscale;

    XLAL_CALLGSL( proj = gsl_matrix_complex_row(projections, row) );
    XLAL_CALLGSL( basisscale = gsl_vector_complex_calloc(TS->size2) );

    XLAL_CALLGSL( model = gsl_matrix_complex_row(TS, row) ); /* get model from training set */

    cprod = complex_weighted_dot_product(weight, &basis.vector, &model.vector);

    XLAL_CALLGSL( gsl_vector_complex_memcpy(basisscale, &basis.vector) );
    XLAL_CALLGSL( gsl_vector_complex_scale(basisscale, cprod) );
    XLAL_CALLGSL( gsl_vector_complex_add(&proj.vector, basisscale) );
    XLAL_CALLGSL( gsl_vector_complex_free(basisscale) );
  }
}


/** \brief The dot product of two real vectors scaled by a given weight factor
 *
 * @param[in] weight A (set of) scaling factor(s) for the dot product
 * @param[in] a The first vector
 * @param[in] b The second vector
 *
 * @return The real dot product of the two vectors
 */
REAL8 weighted_dot_product(gsl_vector *weight, gsl_vector *a, gsl_vector *b){
  REAL8 dp;
  gsl_vector *weighted;

  XLAL_CHECK_REAL8( a->size == b->size, XLAL_EFUNC, "Size of input vectors are not the same.");

  XLAL_CALLGSL( weighted = gsl_vector_calloc(a->size) );
  XLAL_CALLGSL( gsl_vector_memcpy(weighted, a) );

  /* multiply vector by weight */
  if ( weight->size == 1 ){ /* just a single weight to scale with */
    XLAL_CALLGSL( gsl_vector_scale(weighted, gsl_vector_get(weight, 0)) );
  }
  else if ( weight->size == a->size ){
    XLAL_CALLGSL( gsl_vector_mul(weighted, weight) );
  }
  else{
    XLAL_ERROR_REAL8( XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors." );
  }

  /* get dot product */
  XLAL_CALLGSL( gsl_blas_ddot(weighted, b, &dp) );

  XLAL_CALLGSL( gsl_vector_free(weighted) );

  return dp;
}


/** \brief The dot product of two complex vectors scaled by a given weight factor
 *
 * The dot product is produced using the complex conjugate of the first vector.
 *
 * @param[in] weight A real scaling factor for the dot product
 * @param[in] a The first complex vector
 * @param[in] b The second complex vector
 *
 * @return The absolute value of the complex dot product of the two vectors
 */
gsl_complex complex_weighted_dot_product(gsl_vector *weight, gsl_vector_complex *a, gsl_vector_complex *b){
  gsl_complex dp;
  gsl_vector_complex *weighted;

  if ( a->size != b->size ){ XLAL_PRINT_ERROR( "Size of input vectors are not the same." ); }

  XLAL_CALLGSL( weighted = gsl_vector_complex_calloc(a->size) );
  XLAL_CALLGSL( gsl_vector_complex_memcpy(weighted, a) );

  /* multiply vector by weight */
  if ( weight->size == 1 ){ /* just a single weight to scale with */
    XLAL_CALLGSL( gsl_blas_zdscal(gsl_vector_get(weight, 0), weighted) );
  }
  else if ( weight->size == a->size ){
    gsl_vector_view rview, iview;

    XLAL_CALLGSL( rview = gsl_vector_complex_real(weighted) );
    XLAL_CALLGSL( iview = gsl_vector_complex_imag(weighted) );

    XLAL_CALLGSL( gsl_vector_mul(&rview.vector, weight) );
    XLAL_CALLGSL( gsl_vector_mul(&iview.vector, weight) );
  }
  else{
    XLAL_PRINT_ERROR( "Vector of weights must either contain a single value, or be the same length as the other input vectors." );
  }

  /* get dot product (taking the complex conjugate of the first vector) */
  XLAL_CALLGSL( gsl_blas_zdotc(weighted, b, &dp) );
  XLAL_CALLGSL( gsl_vector_complex_free(weighted) );

  return dp;
}


/** \brief Normalise a real vector with a given weighting
 *
 * @param[in] weight The weighting(s) in the normalisation (e.g. time of frequency step(s) between points)
 * @param[in] a The vector to be normalise (this will be change by the function to return the
 * normalised vector.
 */
void normalise(gsl_vector *weight, gsl_vector *a){
  double norm;

  if ( weight->size == 1 ){
    XLAL_CALLGSL( norm = gsl_blas_dnrm2(a) ); /* use GSL normalisation calculation function */
    XLAL_CALLGSL( gsl_vector_scale(a, 1./(norm*sqrt(gsl_vector_get(weight, 0)))) );
  }
  else if ( weight->size == a->size ){
    norm = 1./sqrt(weighted_dot_product(weight, a, a));
    XLAL_CALLGSL( gsl_vector_scale(a, norm) );
  }
  else{
    XLAL_ERROR_VOID( XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors." );
  }
}


/** \brief Normalise a complex vector with a given (real) weighting
 *
 * @param[in] weight The weighting(s) in the normalisation (e.g. time of frequency step(s) between points)
 * @param[in] a The vector to be normalise (this will be change by the function to return the
 * normalised vector.
 */
void complex_normalise(gsl_vector *weight, gsl_vector_complex *a){
  double norm;

  if ( weight->size == 1 ){
    XLAL_CALLGSL( norm = gsl_blas_dznrm2(a) ); /* use GSL normalisation calculation function */
    XLAL_CALLGSL( gsl_blas_zdscal(1./(norm*sqrt(gsl_vector_get(weight, 0))), a) );
  }
  else if ( weight->size == a->size ){
    norm = 1./sqrt(gsl_complex_abs(complex_weighted_dot_product(weight, a, a)));
    XLAL_CALLGSL( gsl_blas_zdscal(norm, a) );
  }
  else{
    XLAL_ERROR_VOID( XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors." );
  }
}


/** \brief Normalise the set of training waveforms
 *
 * This function will normalise a set of training waveforms. This will be used within the
 * \c LALInferenceCreateREAL8OrthonormalBasis function.
 *
 * @param[in] weight The e.g. time/frequency step in the waveforms used to normalise the waveforms
 * @param[in] TS The training set to be normalised.
 */
void normalise_training_set(gsl_vector *weight, gsl_matrix *TS){
  gsl_vector_view rowview;
  size_t i = 0;
  for ( i=0; i<TS->size1; i++ ){
    XLAL_CALLGSL( rowview = gsl_matrix_row(TS, i) );
    normalise(weight, &rowview.vector);
  }
}


/** \brief Normalise the set of complex training waveforms
 *
 * This function will normalise a set of complex training waveforms. This will be used within the
 * \c LALInferenceCreateCOMPLEX16OrthonormalBasis function.
 *
 * @param[in] weight The e.g. time/frequency step in the waveforms used to normalise the waveforms
 * @param[in] TS The training set to be normalised.
 */
void complex_normalise_training_set(gsl_vector *weight, gsl_matrix_complex *TS){
  gsl_vector_complex_view rowview;
  size_t i = 0;

  for ( i=0; i<TS->size1; i++ ){
    XLAL_CALLGSL( rowview = gsl_matrix_complex_row(TS, i) );
    complex_normalise(weight, &rowview.vector);
  }
}


/** \brief Get the interpolant of a reduced basis set
 *
 * This is used internally by \c LALInferenceCreateREALROQInterpolant when
 * iteratively calculating the interpolation matrix.
 *
 * @param[in] V The matrix containing the basis vector points at the current interpolation nodes
 * @param[in] RB The set of basis vectors
 *
 * @return The interpolant matrix
 */
REAL8Array *B_matrix(gsl_matrix *V, gsl_matrix *RB){
  /* get inverse of V */
  size_t n = V->size1;
  gsl_matrix *invV, *LU;
  gsl_permutation *p;
  gsl_matrix_view subRB;
  REAL8Array *B = NULL;
  UINT4Vector *dims = NULL;

  XLAL_CALLGSL( invV = gsl_matrix_alloc(n, n) );
  int signum;

  /* use LU decomposition to get inverse */
  XLAL_CALLGSL( LU = gsl_matrix_alloc(n, n) );
  XLAL_CALLGSL( gsl_matrix_memcpy(LU, V) );

  XLAL_CALLGSL( p = gsl_permutation_alloc(n) );
  XLAL_CALLGSL( gsl_linalg_LU_decomp(LU, p, &signum) );
  XLAL_CALLGSL( gsl_linalg_LU_invert(LU, p, invV) );
  XLAL_CALLGSL( gsl_permutation_free(p) );
  XLAL_CALLGSL( gsl_matrix_free(LU) );

  /* get B matrix */
  dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = n;
  dims->data[1] = RB->size2;
  B = XLALCreateREAL8Array( dims );
  XLALDestroyUINT4Vector( dims );
  gsl_matrix_view Bview;
  Bview = gsl_matrix_view_array(B->data, n, RB->size2);
  XLAL_CALLGSL( subRB = gsl_matrix_submatrix(RB, 0, 0, n, RB->size2) );
  XLAL_CALLGSL( gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, invV, &subRB.matrix, 0., &Bview.matrix) );

  XLAL_CALLGSL( gsl_matrix_free(invV) );

  return B;
}


/** \brief Get the interpolant of a complex reduced basis set
 *
 * This is used internally by \c LALInferenceCreateCOMPLEXROQInterpolant when
 * iteratively calculating the interpolation matrix.
 *
 * @param[in] V The matrix containing the basis vector points at the current interpolation nodes
 * @param[in] RB The set of basis vectors
 *
 * @return The interpolant matrix
 */
COMPLEX16Array *complex_B_matrix(gsl_matrix_complex *V, gsl_matrix_complex *RB){
  /* get inverse of V */
  size_t n = V->size1;
  gsl_matrix_complex *invV, *LU;
  gsl_permutation *p;
  gsl_matrix_complex_view subRB;
  gsl_complex scale1, scale0;
  COMPLEX16Array *B = NULL;
  UINT4Vector *dims = NULL;

  GSL_SET_COMPLEX(&scale1, 1., 0.);
  GSL_SET_COMPLEX(&scale0, 0., 0.);

  XLAL_CALLGSL( invV = gsl_matrix_complex_alloc(n, n) );
  int signum;

  /* use LU decomposition to get inverse */
  XLAL_CALLGSL( LU = gsl_matrix_complex_alloc(n, n) );
  XLAL_CALLGSL( gsl_matrix_complex_memcpy(LU, V) );

  XLAL_CALLGSL( p = gsl_permutation_alloc(n) );
  XLAL_CALLGSL( gsl_linalg_complex_LU_decomp(LU, p, &signum) );
  XLAL_CALLGSL( gsl_linalg_complex_LU_invert(LU, p, invV) );
  XLAL_CALLGSL( gsl_permutation_free(p) );
  XLAL_CALLGSL( gsl_matrix_complex_free(LU) );

  /* get B matrix */
  dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = n;
  dims->data[1] = RB->size2;
  B = XLALCreateCOMPLEX16Array( dims );
  XLALDestroyUINT4Vector( dims );
  gsl_matrix_complex_view Bview;
  Bview = gsl_matrix_complex_view_array((double *)B->data, n, RB->size2);
  XLAL_CALLGSL( subRB = gsl_matrix_complex_submatrix(RB, 0, 0, n, RB->size2) );
  XLAL_CALLGSL( gsl_blas_zgemm(CblasTrans, CblasNoTrans, scale1, invV, &subRB.matrix, scale0, &Bview.matrix) );

  XLAL_CALLGSL( gsl_matrix_complex_free(invV) );

  return B;
}

/** \brief Get the index of the maximum absolute value for a complex vector
 *
 * @param[in] c A complex vector
 *
 * @return The index of the maximum absolute value of that vector
 */
int complex_vector_maxabs_index( gsl_vector_complex *c ){
  double maxv = -INFINITY, absval = 0.;
  int idx = 0;
  size_t i = 0;

  for ( i=0; i<c->size; i++ ){
    XLAL_CALLGSL( absval = gsl_complex_abs(gsl_vector_complex_get(c, i)) );

    if ( absval > maxv ){
      maxv = absval;
      idx = (int)i;
    }
  }

  return idx;
}


/* main functions */

/**
 * \brief Create a orthonormal basis set from a training set of real waveforms
 *
 * Given a \c gsl_matrix containing a training set of real waveforms (where the waveforms
 * are created at time or frequency steps seperated by \c delta) an orthonormal basis
 * will be generated using the greedy binning Algorithm 1 of \cite FGHKT2014 . The stopping
 * criteria for the algorithm is controlled by the \c tolerance value, which defined the
 * maximum residual between the current basis set (at a given iteration) and the training
 * set (and example tolerance is \f$10^{-12}\f$). In this function the training set will be
 * normalised, so the input \c TS will be modified.
 *
 * If the \c RBin value is \c NULL then a new reduced basis will be formed from the given
 * training set. However, if \c RBin already contains a previously produced basis, then this
 * basis will be enriched with bases if possible using the new training set.  <b>NOTE</b>: when
 * using  small tolerances enriching the basis in this way can lead to numerical precision issues,
 * so in general you should use \c LALInferenceEnrichREAL8Basis for enrichment.
 *
 * @param[out] RBin A \c REAL8Array to return the reduced basis.
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] tolerance The tolerance used as a stopping criteria for the basis generation.
 * @param[in] TS A \c REAL8Array matrix containing the complex training set, where the number of waveforms in the
 * training set is given by the rows and the waveform points by the columns.
 *
 * @return A \c REAL8 with the maximum projection error for the final reduced basis.
 *
 * \sa LALInferenceEnrichREAL8Basis
 */
REAL8 LALInferenceGenerateREAL8OrthonormalBasis(REAL8Array **RBin,
                                                REAL8Vector *delta,
                                                REAL8 tolerance,
                                                REAL8Array *TS){
  XLAL_CHECK_REAL8( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK_REAL8( TS != NULL, XLAL_EFUNC, "Training set array is NULL!" );
  XLAL_CHECK_REAL8( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );
  XLAL_CHECK_REAL8( TS->dimLength->length == 2, XLAL_EFUNC, "Training set array must have only two dimensions" );

  REAL8Array *RB = NULL;
  UINT4Vector *dims = NULL;
  gsl_matrix_view RBview;

  gsl_matrix *projections; /* projections of the basis onto the training set */
  gsl_matrix *residual;
  gsl_vector *projection_errors;

  REAL8 sigma = 1., prevsigma = 1.;
  size_t dlength = TS->dimLength->data[1], nts = TS->dimLength->data[0];
  size_t mindex = 0, k=0;
  UINT4 idx = 0, nprev = 0, j = 0;

  /* normalise the training set */
  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  gsl_matrix_view TSview;
  XLAL_CALLGSL( TSview = gsl_matrix_view_array(TS->data, nts, dlength) );
  normalise_training_set(&deltaview.vector, &TSview.matrix);

  XLAL_CALLGSL( projections = gsl_matrix_calloc(nts, dlength) );

  /* if a NULL basis has been passed then create reduced basis (initially just one model vector in length) */
  dims = XLALCreateUINT4Vector( 2 );
  if ( *RBin == NULL ){
    gsl_vector *firstrow;
    dims->data[0] = 1; /* one row */
    dims->data[1] = dlength;
    RB = XLALCreateREAL8Array( dims );

    XLAL_CALLGSL( RBview = gsl_matrix_view_array(RB->data, 1, dlength) );
    XLAL_CALLGSL( firstrow = gsl_vector_calloc(dlength) );
    XLAL_CALLGSL( gsl_matrix_get_row(firstrow, &TSview.matrix, 0) );
    XLAL_CALLGSL( gsl_matrix_set_row(&RBview.matrix, 0, firstrow) );
    XLAL_CALLGSL( gsl_vector_free(firstrow) );
    *RBin = RB;
    nprev = 1;
  }
  else{ /* if the reduced basis already exists then view it as a gsl_matrix */
    RB = *RBin;
    dims->data[0] = RB->dimLength->data[0];
    dims->data[1] = RB->dimLength->data[1];
    XLAL_CHECK_REAL8( dims->data[1] == dlength, XLAL_EFUNC, "Training set and reduced basis must have the same basis lengths" );
    XLAL_CALLGSL( RBview = gsl_matrix_view_array(RB->data, dims->data[0], dims->data[1]) );
    nprev = dims->data[0];
  }

  for ( j = 0; j < nprev; j++ ){ project_onto_basis( &deltaview.vector, &RBview.matrix, &TSview.matrix, projections, j ); }

  XLAL_CALLGSL( projection_errors = gsl_vector_calloc(nts) );
  XLAL_CALLGSL( residual = gsl_matrix_calloc(nts, dlength) );

  /* create reduced basis set using greedy binning Algorithm 1 of http://arxiv.org/abs/1308.3565 */
  while ( 1 ){
    gsl_vector *next_basis;
    prevsigma = sigma; /* check that sigma is decreasing */

    if ( idx > 0 ){ project_onto_basis(&deltaview.vector, &RBview.matrix, &TSview.matrix, projections, idx+nprev-1); }

    XLAL_CALLGSL( gsl_matrix_memcpy(residual, &TSview.matrix) ); /* copy training set into residual */

    /* get residuals by subtracting projections from training set */
    XLAL_CALLGSL( gsl_matrix_sub(residual, projections) );

    /* get projection errors */
    for( k=0; k < nts; k++ ){
      REAL8 err;
      gsl_vector_view resrow;

      XLAL_CALLGSL( resrow = gsl_matrix_row(residual, k) );

      err = weighted_dot_product(&deltaview.vector, &resrow.vector, &resrow.vector);

      XLAL_CALLGSL( gsl_vector_set(projection_errors, k, err) );
    }

    sigma = fabs(gsl_vector_max(projection_errors));

    //if ( sigma > 1e-5 ){ fprintf(stderr, "%.12lf\t%d\n", sigma, idx); }
    //else { fprintf(stderr, "%.12le\t%d\n", sigma, idx); }

    idx++;

    /* break point */
    if ( sigma < tolerance || idx == (UINT4)nts || sigma > prevsigma ) {
      if ( idx == (UINT4)nts ){
        XLAL_PRINT_WARNING( "Not enough training models (%zu) to produce orthonormal basis given the tolerance of %le\n", nts, tolerance );
      }
      if ( sigma > prevsigma ){
        XLAL_PRINT_WARNING( "Numerical precision issues have cause the residuals to increase (%12.le -> %.12le). Stop adding more bases with this training set.\n", prevsigma, sigma );
        /* remove previously added reduced basis */
        dims->data[0] = dims->data[0]-1;
        RB = XLALResizeREAL8Array( RB, dims );
      }
      break;
    }

    /* get index of training set with the largest projection errors */
    XLAL_CALLGSL( mindex = gsl_vector_max_index( projection_errors ) );

    XLAL_CALLGSL( next_basis = gsl_vector_calloc(dlength) );
    XLAL_CALLGSL( gsl_matrix_get_row(next_basis, residual, mindex) );

    /* normalise vector */
    normalise(&deltaview.vector, next_basis);

    /* expand reduced basis */
    dims->data[0] = idx+nprev; /* add row */
    RB = XLALResizeREAL8Array( RB, dims );

    /* add on next basis */
    XLAL_CALLGSL( RBview = gsl_matrix_view_array(RB->data, idx+nprev, dlength) );
    XLAL_CALLGSL( gsl_matrix_set_row(&RBview.matrix, idx+nprev-1, next_basis) );

    XLAL_CALLGSL( gsl_vector_free(next_basis) );
  }

  /* free memory */
  XLAL_CALLGSL( gsl_matrix_free(projections) );
  XLAL_CALLGSL( gsl_vector_free(projection_errors) );
  XLAL_CALLGSL( gsl_matrix_free(residual) );
  XLALDestroyUINT4Vector( dims );

  return sigma;
}


/**
 * \brief Create a orthonormal basis set from a training set of complex waveforms
 *
 * Given a \c gsl_matrix containing a training set \c TS of complex waveforms (where the waveforms
 * are created at time or frequency steps seperated by \c delta) an orthonormal basis
 * will be generated using the greedy binning Algorithm 1 of \cite FGHKT2014 . The stopping
 * criteria for the algorithm is controlled by the \c tolerance value, which is defined the
 * maximum residual between the current basis set (at a given iteration) and the training
 * set (and example tolerance is \f$10^{-12}\f$). In this function the training set will be
 * normalised, so the input \c TS will be modified.
 *
 * If the \c RBin value is \c NULL then a new reduced basis will be formed from the given
 * training set. However, if \c RBin already contains a previously produced basis, then this
 * basis will be enriched with bases if possible using the new training set. <b>NOTE</b>: when
 * using  small tolerances enriching the basis in this way can lead to numerical precision issues,
 * so in general you should use \c LALInferenceEnrichCOMPLEX16Basis for enrichment.
 *
 * Note that in this function we have to cast the \c COMPLEX16 array as a double to use
 * \c gsl_matrix_view_array, which assume that the data is passed as a double array with
 * memory laid out so that adjacent double memory blocks hold the corresponding real and
 * imaginary parts.
 *
 * @param[out] RBin A \c COMPLEX16Array to return the reduced basis.
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] tolerance The tolerance used as a stopping criteria for the basis generation.
 * @param[in] TS A \c COMPLEX16Array matrix containing the complex training set, where the number of waveforms in the
 * training set is given by the rows and the waveform points by the columns.
 *
 * @return A \c REAL8 with the maximum projection error for the final reduced basis.
 *
 * \sa LALInferenceEnrichCOMPLEX16Basis
 */
REAL8 LALInferenceGenerateCOMPLEX16OrthonormalBasis(COMPLEX16Array **RBin,
                                                    REAL8Vector *delta,
                                                    REAL8 tolerance,
                                                    COMPLEX16Array *TS){
  XLAL_CHECK_REAL8( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK_REAL8( TS != NULL, XLAL_EFUNC, "Training set array is NULL!" );
  XLAL_CHECK_REAL8( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );
  XLAL_CHECK_REAL8( TS->dimLength->length == 2, XLAL_EFUNC, "Training set array must have only two dimensions" );

  COMPLEX16Array *RB = NULL;
  UINT4Vector *dims = NULL;
  gsl_matrix_complex_view RBview;

  gsl_matrix_complex *projections; /* projections of the basis onto the training set */
  gsl_matrix_complex *residual;
  gsl_vector *projection_errors;

  REAL8 sigma = 1., prevsigma = 1.;
  size_t dlength = TS->dimLength->data[1], nts = TS->dimLength->data[0];
  size_t mindex = 0, k=0;
  UINT4 idx = 0, nprev = 0, j = 0;

  /* normalise the training set */
  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  gsl_matrix_complex_view TSview;
  XLAL_CALLGSL( TSview = gsl_matrix_complex_view_array((double *)TS->data, nts, dlength) );
  complex_normalise_training_set(&deltaview.vector, &TSview.matrix);

  XLAL_CALLGSL( projections = gsl_matrix_complex_calloc(nts, dlength) );

  /* if a NULL basis has been passed then create reduced basis (initially just one model vector in length) */
  dims = XLALCreateUINT4Vector( 2 );
  if ( *RBin == NULL ){
    gsl_vector_complex *firstrow;
    dims->data[0] = 1; /* one row */
    dims->data[1] = dlength;
    RB = XLALCreateCOMPLEX16Array( dims );

    XLAL_CALLGSL( RBview = gsl_matrix_complex_view_array((double*)RB->data, 1, dlength) );
    XLAL_CALLGSL( firstrow = gsl_vector_complex_calloc(dlength) );
    XLAL_CALLGSL( gsl_matrix_complex_get_row(firstrow, &TSview.matrix, 0) );
    XLAL_CALLGSL( gsl_matrix_complex_set_row(&RBview.matrix, 0, firstrow) );
    XLAL_CALLGSL( gsl_vector_complex_free(firstrow) );
    *RBin = RB;
    nprev = 1;
  }
  else{ /* if the reduced basis already exists then view it as a gsl_matrix */
    RB = *RBin;
    dims->data[0] = RB->dimLength->data[0];
    dims->data[1] = RB->dimLength->data[1];
    XLAL_CHECK_REAL8( dims->data[1] == dlength, XLAL_EFUNC, "Training set and reduced basis must have the same basis lengths" );
    XLAL_CALLGSL( RBview = gsl_matrix_complex_view_array((double*)RB->data, dims->data[0], dims->data[1]) );
    nprev = dims->data[0];
  }

  for ( j = 0; j < nprev; j++ ){ complex_project_onto_basis( &deltaview.vector, &RBview.matrix, &TSview.matrix, projections, j ); }

  XLAL_CALLGSL( projection_errors = gsl_vector_calloc(nts) );
  XLAL_CALLGSL( residual = gsl_matrix_complex_calloc(nts, dlength) );

  /* create reduced basis set using greedy binning Algorithm 1 of http://arxiv.org/abs/1308.3565 */
  while ( 1 ){
    gsl_vector_complex *next_basis;
    prevsigma = sigma; /* check that sigma is decreasing */

    if ( idx > 0 ){ complex_project_onto_basis(&deltaview.vector, &RBview.matrix, &TSview.matrix, projections, idx+nprev-1); }

    XLAL_CALLGSL( gsl_matrix_complex_memcpy(residual, &TSview.matrix) ); /* copy training set into residual */

    /* get residuals by subtracting projections from training set */
    XLAL_CALLGSL( gsl_matrix_complex_sub(residual, projections) );

    /* get projection errors */
    for( k=0; k < nts; k++ ){
      gsl_complex err;
      gsl_vector_complex_view resrow;

      XLAL_CALLGSL( resrow = gsl_matrix_complex_row(residual, k) );

      err = complex_weighted_dot_product(&deltaview.vector, &resrow.vector, &resrow.vector);

      XLAL_CALLGSL( gsl_vector_set(projection_errors, k, GSL_REAL(err)) );
    }

    sigma = fabs(gsl_vector_max(projection_errors));

    //if ( sigma > 1e-5 ){ fprintf(stderr, "%.12lf\t%d\n", sigma, idx); }
    //else { fprintf(stderr, "%.12le\t%d\n", sigma, idx); }

    idx++;

    /* break point */
    if ( sigma < tolerance || idx == (UINT4)nts || sigma > prevsigma ) {
      if ( idx == (UINT4)nts ){
        XLAL_PRINT_WARNING( "Not enough training models (%zu) to produce orthonormal basis given the tolerance of %le\n", nts, tolerance );
      }
      if ( sigma > prevsigma ){
        XLAL_PRINT_WARNING( "Numerical precision issues have cause the residuals to increase (%12.le -> %.12le). Stop adding more bases with this training set.\n", prevsigma, sigma );
        /* remove previously added reduced basis */
        dims->data[0] = dims->data[0]-1;
        RB = XLALResizeCOMPLEX16Array( RB, dims );
      }
      break;
    }

    /* get index of training set with the largest projection errors */
    XLAL_CALLGSL( mindex = gsl_vector_max_index( projection_errors ) );

    XLAL_CALLGSL( next_basis = gsl_vector_complex_calloc(dlength) );
    XLAL_CALLGSL( gsl_matrix_complex_get_row( next_basis, residual, mindex ) );

    /* normalise vector */
    complex_normalise(&deltaview.vector, next_basis);

    /* expand reduced basis */
    dims->data[0] = idx+nprev; /* add row */
    RB = XLALResizeCOMPLEX16Array( RB, dims );

    /* add on next basis */
    XLAL_CALLGSL( RBview = gsl_matrix_complex_view_array((double*)RB->data, idx+nprev, dlength) );
    XLAL_CALLGSL( gsl_matrix_complex_set_row(&RBview.matrix, idx+nprev-1, next_basis) );

    XLAL_CALLGSL( gsl_vector_complex_free(next_basis) );
  }

  /* free memory */
  XLAL_CALLGSL( gsl_matrix_complex_free(projections) );
  XLAL_CALLGSL( gsl_vector_free(projection_errors) );
  XLAL_CALLGSL( gsl_matrix_complex_free(residual) );
  XLALDestroyUINT4Vector( dims );

  return sigma;
}


/**
 * \brief Test the real reduced basis against another set of waveforms
 *
 * This function projects a set of waveforms onto the reduced basis and
 * checks that the residuals are within a given tolerance
 *
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] tolerance The allowed residual tolerence for the test waveforms
 * @param[in] RB The reduced basis set
 * @param[in] testmodels The set of waveform models to project onto the basis
 *
 * @return Returns \c XLAL_SUCCESS if all test waveforms meet the tolerance
 */
INT4 LALInferenceTestREAL8OrthonormalBasis(REAL8Vector *delta,
                                           REAL8 tolerance,
                                           REAL8Array *RB,
                                           REAL8Array *testmodels){
  XLAL_CHECK( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK( RB != NULL, XLAL_EFUNC, "Reduced basis set array is NULL!" );
  XLAL_CHECK( RB->dimLength->length == 2, XLAL_EFUNC, "Reduced basis set array must have only two dimensions" );
  XLAL_CHECK( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );
  XLAL_CHECK( testmodels != NULL, XLAL_EFUNC, "Set of test models is NULL!" );
  XLAL_CHECK( testmodels->dimLength->length == 2, XLAL_EFUNC, "Test set does not have 2 dimensions!" );

  size_t dlength = testmodels->dimLength->data[1], nts = testmodels->dimLength->data[0];
  size_t j = 0, k = 0;

  /* normalise the test set */
  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  gsl_matrix_view testmodelsview;
  XLAL_CALLGSL( testmodelsview = gsl_matrix_view_array(testmodels->data, nts, dlength) );
  normalise_training_set(&deltaview.vector, &testmodelsview.matrix);

  gsl_matrix_view RBview;
  RBview = gsl_matrix_view_array(RB->data, RB->dimLength->data[0], dlength);

  /* get projection errors for each test model */
  for ( k = 0; k < nts; k++ ){
    REAL8 projerr = 0.;
    gsl_vector *proj;
    gsl_vector_view testrow;
    XLAL_CALLGSL( proj = gsl_vector_calloc(dlength) ); /* allocated to zero */
    XLAL_CALLGSL( testrow = gsl_matrix_row(&testmodelsview.matrix, k) );

    /* loop over reduced basis getting the projection coefficients */
    for ( j = 0; j < RBview.matrix.size1; j++ ){
      REAL8 projcoeff;
      gsl_vector *RBrow;
      XLAL_CALLGSL( RBrow = gsl_vector_alloc( RBview.matrix.size2 ) );
      XLAL_CALLGSL( gsl_matrix_get_row( RBrow, &RBview.matrix, j ) );

      /* get dot product of reduced basis vector with test model */
      projcoeff = weighted_dot_product(&deltaview.vector, RBrow, &testrow.vector);

      XLAL_CALLGSL( gsl_vector_scale( RBrow, projcoeff ) );
      XLAL_CALLGSL( gsl_vector_add( proj, RBrow ) );
      XLAL_CALLGSL( gsl_vector_free( RBrow ) );
    }

    /* get residual */
    XLAL_CALLGSL( gsl_vector_sub( proj, &testrow.vector ) );
    projerr = weighted_dot_product(&deltaview.vector, proj, proj);

    /* check projection error against tolerance */
    if ( fabs(projerr) > tolerance ) { return XLAL_FAILURE; }

    XLAL_CALLGSL( gsl_vector_free( proj ) );
  }

  return XLAL_SUCCESS;
}


/**
 * \brief Expand the real training waveforms with ones from a set of new training waveforms
 *
 * Expands an original set of training waveforms with waveforms from a new set a training waveforms
 * if, when projected onto the reduced basis, those new waveforms are greater than the given tolerance.
 * The reduced basis will then be recomputed using the expanded training set.
 *
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] tolerance The allowed residual tolerence for the test waveforms
 * @param[in] RB The reduced basis set
 * @param[in] testmodels The set of waveform models used to produce the reduced basis
 * @param[in] testmodelsnew A new set of waveforms to project onto the current reduced basis
 *
 * @return Returns the maximum projection error of the new basis set
 */
REAL8 LALInferenceEnrichREAL8Basis(REAL8Vector *delta,
                                   REAL8 tolerance,
                                   REAL8Array **RB,
                                   REAL8Array **testmodels,
                                   REAL8Array *testmodelsnew){
  XLAL_CHECK( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );

  size_t dlength = testmodelsnew->dimLength->data[1], nts = testmodelsnew->dimLength->data[0];
  size_t j = 0, k = 0;

  REAL8 maxprojerr = 0.;

  /* normalise the test set */
  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  gsl_matrix_view testmodelsview;
  XLAL_CALLGSL( testmodelsview = gsl_matrix_view_array(testmodelsnew->data, nts, dlength) );
  normalise_training_set(&deltaview.vector, &testmodelsview.matrix);

  REAL8Array *RBin = NULL;
  RBin = *RB;
  REAL8Array *tm = NULL;
  tm = *testmodels;

  gsl_matrix_view RBview;
  RBview = gsl_matrix_view_array(RBin->data, RBin->dimLength->data[0], dlength);

  XLAL_CHECK_REAL8( dlength == tm->dimLength->data[1], XLAL_EFUNC, "New training set contains waveforms of different length to the original set!" );

  UINT4Vector *dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = tm->dimLength->data[0];
  dims->data[1] = dlength;
  size_t ntsorig = tm->dimLength->data[0];

  /* get projection errors for each test model */
  for ( k = 0; k < nts; k++ ){
    REAL8 projerr = 0.;
    gsl_vector *proj;
    gsl_vector_view testrow;
    XLAL_CALLGSL( proj = gsl_vector_calloc(dlength) ); /* allocated to zero */
    XLAL_CALLGSL( testrow = gsl_matrix_row(&testmodelsview.matrix, k) );

    /* loop over reduced basis getting the projection coefficients */
    for ( j = 0; j < RBview.matrix.size1; j++ ){
      REAL8 projcoeff;
      gsl_vector *RBrow;
      XLAL_CALLGSL( RBrow = gsl_vector_alloc( RBview.matrix.size2 ) );
      XLAL_CALLGSL( gsl_matrix_get_row( RBrow, &RBview.matrix, j ) );

      /* get dot product of reduced basis vector with test model */
      projcoeff = weighted_dot_product(&deltaview.vector, RBrow, &testrow.vector);

      XLAL_CALLGSL( gsl_vector_scale( RBrow, projcoeff ) );
      XLAL_CALLGSL( gsl_vector_add( proj, RBrow ) );
      XLAL_CALLGSL( gsl_vector_free( RBrow ) );
    }

    /* get residual */
    XLAL_CALLGSL( gsl_vector_sub( proj, &testrow.vector ) );
    projerr = weighted_dot_product(&deltaview.vector, proj, proj);

    /* check projection error against tolerance and add to training waveforms if required */
    if ( fabs(projerr) > tolerance ) {
      gsl_matrix_view tmview;
      dims->data[0] = tm->dimLength->data[0] + 1;
      tm = XLALResizeREAL8Array( tm, dims );
      tmview = gsl_matrix_view_array(tm->data, tm->dimLength->data[0], tm->dimLength->data[1]);
      gsl_matrix_set_row(&tmview.matrix, tm->dimLength->data[0]-1, &testrow.vector);
    }
    if ( fabs(projerr) > maxprojerr ){ maxprojerr = fabs(projerr); }

    XLAL_CALLGSL( gsl_vector_free( proj ) );
  }

  XLALDestroyUINT4Vector( dims );

  /* recalculate the reduced basis using the updated training set */
  if ( tm->dimLength->data[0] > ntsorig ){
    XLALDestroyREAL8Array( *RB ); /* remove current basis, so it is rebuilt using the entire new training set */
    *RB = NULL;
    maxprojerr = LALInferenceGenerateREAL8OrthonormalBasis(RB, delta, tolerance, tm);
  }

  return maxprojerr;
}


/**
 * \brief Test the complex reduced basis against another set of waveforms
 *
 * This function projects a set of waveforms onto the reduced basis and
 * checks that the residuals are within a given tolerance
 *
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] tolerance The allowed residual tolerence for the test waveforms
 * @param[in] RB The reduced basis set
 * @param[in] testmodels The set of waveform models to project onto the basis
 *
 * @return Returns \c XLAL_SUCCESS if all test waveforms meet the tolerance
 */
INT4 LALInferenceTestCOMPLEX16OrthonormalBasis(REAL8Vector *delta,
                                               REAL8 tolerance,
                                               COMPLEX16Array *RB,
                                               COMPLEX16Array *testmodels){
  XLAL_CHECK( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );
  XLAL_CHECK( RB != NULL, XLAL_EFUNC, "Reduced basis is NULL!" );
  XLAL_CHECK( RB->dimLength->length == 2, XLAL_EFUNC, "Reduced basis set does not have 2 dimensions!" );
  XLAL_CHECK( testmodels != NULL, XLAL_EFUNC, "Set of test models is NULL!" );
  XLAL_CHECK( testmodels->dimLength->length == 2, XLAL_EFUNC, "Test set does not have 2 dimensions!" );

  size_t dlength = testmodels->dimLength->data[1], nts = testmodels->dimLength->data[0];
  size_t j = 0, k = 0;

  /* normalise the test set */
  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  gsl_matrix_complex_view testmodelsview;
  XLAL_CALLGSL( testmodelsview = gsl_matrix_complex_view_array((double *)testmodels->data, nts, dlength) );
  complex_normalise_training_set(&deltaview.vector, &testmodelsview.matrix);

  gsl_matrix_complex_view RBview;
  RBview = gsl_matrix_complex_view_array((double *)RB->data, RB->dimLength->data[0], dlength);

  /* get projection errors for each test model */
  for ( k = 0; k < nts; k++ ){
    gsl_complex projerr;
    gsl_vector_complex *proj;
    gsl_vector_complex_view testrow;
    XLAL_CALLGSL( proj = gsl_vector_complex_calloc(dlength) ); /* allocated to zero */
    XLAL_CALLGSL( testrow = gsl_matrix_complex_row(&testmodelsview.matrix, k) );

    /* loop over reduced basis getting the projection coefficients */
    for ( j = 0; j < RBview.matrix.size1; j++ ){
      gsl_complex projcoeff;
      gsl_vector_complex *RBrow;
      XLAL_CALLGSL( RBrow = gsl_vector_complex_alloc( RBview.matrix.size2 ) );
      XLAL_CALLGSL( gsl_matrix_complex_get_row( RBrow, &RBview.matrix, j ) );

      /* get dot product of reduced basis vector with test model */
      projcoeff = complex_weighted_dot_product(&deltaview.vector, RBrow, &testrow.vector);

      XLAL_CALLGSL( gsl_vector_complex_scale( RBrow, projcoeff ) );
      XLAL_CALLGSL( gsl_vector_complex_add( proj, RBrow ) );
      XLAL_CALLGSL( gsl_vector_complex_free( RBrow ) );
    }

    /* get residual */
    XLAL_CALLGSL( gsl_vector_complex_sub( proj, &testrow.vector ) );
    projerr = complex_weighted_dot_product(&deltaview.vector, proj, proj);

    //fprintf(stderr, "%zu: Projected error = %le ", k+1, GSL_REAL(projerr));

    /* check projection error against tolerance */
    if ( GSL_REAL(projerr) > tolerance ) { return XLAL_FAILURE; }

    XLAL_CALLGSL( gsl_vector_complex_free( proj ) );
  }

  return XLAL_SUCCESS;
}


/**
 * \brief Expand the complex training waveforms with ones from a set of new training waveforms
 *
 * Expands an original set of training waveforms with waveforms from a new set a training waveforms
 * if, when projected onto the reduced basis, those new waveforms are greater than the given tolerance.
 * The reduced basis will then be recomputed using the expanded training set.
 *
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] tolerance The allowed residual tolerence for the test waveforms
 * @param[in] RB The reduced basis set
 * @param[in] testmodels The set of waveform models used to produce the reduced basis
 * @param[in] testmodelsnew A new set of waveforms to project onto the current reduced basis
 *
 * @return Returns the maximum projection error of the new basis set
 */
REAL8 LALInferenceEnrichCOMPLEX16Basis(REAL8Vector *delta,
                                       REAL8 tolerance,
                                       COMPLEX16Array **RB,
                                       COMPLEX16Array **testmodels,
                                       COMPLEX16Array *testmodelsnew){
  XLAL_CHECK( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );

  size_t dlength = testmodelsnew->dimLength->data[1], nts = testmodelsnew->dimLength->data[0];
  size_t j = 0, k = 0;

  REAL8 maxprojerr = 0.;

  /* normalise the test set */
  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  gsl_matrix_complex_view testmodelsview;
  XLAL_CALLGSL( testmodelsview = gsl_matrix_complex_view_array((double *)testmodelsnew->data, nts, dlength) );
  complex_normalise_training_set(&deltaview.vector, &testmodelsview.matrix);

  COMPLEX16Array *RBin = NULL;
  RBin = *RB;
  COMPLEX16Array *tm = NULL;
  tm = *testmodels;

  gsl_matrix_complex_view RBview;
  RBview = gsl_matrix_complex_view_array((double *)RBin->data, RBin->dimLength->data[0], dlength);

  XLAL_CHECK_REAL8( dlength == tm->dimLength->data[1], XLAL_EFUNC, "New training set contains waveforms of different length to the original set!" );

  UINT4Vector *dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = tm->dimLength->data[0];
  dims->data[1] = dlength;
  size_t ntsorig = tm->dimLength->data[0];

  /* get projection errors for each test model */
  for ( k = 0; k < nts; k++ ){
    gsl_complex projerr;
    gsl_vector_complex *proj;
    gsl_vector_complex_view testrow;
    XLAL_CALLGSL( proj = gsl_vector_complex_calloc(dlength) ); /* allocated to zero */
    XLAL_CALLGSL( testrow = gsl_matrix_complex_row(&testmodelsview.matrix, k) );

    /* loop over reduced basis getting the projection coefficients */
    for ( j = 0; j < RBview.matrix.size1; j++ ){
      gsl_complex projcoeff;
      gsl_vector_complex *RBrow;
      XLAL_CALLGSL( RBrow = gsl_vector_complex_alloc( RBview.matrix.size2 ) );
      XLAL_CALLGSL( gsl_matrix_complex_get_row( RBrow, &RBview.matrix, j ) );

      /* get dot product of reduced basis vector with test model */
      projcoeff = complex_weighted_dot_product(&deltaview.vector, RBrow, &testrow.vector);

      XLAL_CALLGSL( gsl_vector_complex_scale( RBrow, projcoeff ) );
      XLAL_CALLGSL( gsl_vector_complex_add( proj, RBrow ) );
      XLAL_CALLGSL( gsl_vector_complex_free( RBrow ) );
    }

    /* get residual */
    XLAL_CALLGSL( gsl_vector_complex_sub( proj, &testrow.vector ) );
    projerr = complex_weighted_dot_product(&deltaview.vector, proj, proj);

    //fprintf(stderr, "%zu: Projected error = %le\n", k+1, GSL_REAL(projerr));

    /* check projection error against tolerance */
    if ( GSL_REAL(projerr) > tolerance ) {
      gsl_matrix_complex_view tmview;
      dims->data[0] = tm->dimLength->data[0] + 1;
      tm = XLALResizeCOMPLEX16Array( tm, dims );
      tmview = gsl_matrix_complex_view_array((double *)tm->data, tm->dimLength->data[0], tm->dimLength->data[1]);
      gsl_matrix_complex_set_row(&tmview.matrix, tm->dimLength->data[0]-1, &testrow.vector);
    }
    if ( GSL_REAL(projerr) > maxprojerr ){ maxprojerr = GSL_REAL(projerr); }

    XLAL_CALLGSL( gsl_vector_complex_free( proj ) );
  }

  XLALDestroyUINT4Vector( dims );

  /* recalculate the reduced basis using the updated training set */
  if ( tm->dimLength->data[0] > ntsorig ){
    XLALDestroyCOMPLEX16Array( *RB ); /* remove current basis, so it is rebuilt using the entire new training set */
    *RB = NULL;
    maxprojerr = LALInferenceGenerateCOMPLEX16OrthonormalBasis(RB, delta, tolerance, tm);
  }

  return maxprojerr;
}


/**
 * \brief Create a real empirical interpolant from a set of orthonormal basis functions
 *
 * Given a real \c REAL8Array containing a set of orthonormal basis functions generate an
 * empirical intopolant, and set of interpolation points, using Algorithm 2 of
 * \cite FGHKT2014 .
 *
 * @param[in] RB The set of basis functions
 *
 * @return A \c LALInferenceREALROQInterpolant structure containing the interpolant and its nodes
 */
LALInferenceREALROQInterpolant *LALInferenceGenerateREALROQInterpolant(REAL8Array *RB){
  XLAL_CHECK_NULL( RB != NULL, XLAL_EFUNC, "Reduced basis array is NULL!" );

  size_t RBsize = RB->dimLength->data[0];  /* reduced basis size (no. of reduced bases) */
  size_t dlength = RB->dimLength->data[1]; /* length of each base */
  size_t i=1, j=0, k=0;
  REAL8 *V = XLALMalloc(sizeof(REAL8));
  gsl_matrix_view Vview;

  LALInferenceREALROQInterpolant *interp = XLALMalloc(sizeof(LALInferenceREALROQInterpolant));
  int idmax = 0, newidx = 0;

  /* get index of maximum absolute value of first basis */
  gsl_matrix_view RBview;
  RBview = gsl_matrix_view_array( RB->data, RBsize, dlength );
  gsl_vector_view firstbasis = gsl_matrix_row(&RBview.matrix, 0);
  XLAL_CALLGSL( idmax = (int)gsl_blas_idamax(&firstbasis.vector) ); /* function gets index of maximum absolute value */

  interp->nodes = XLALMalloc(RBsize*sizeof(UINT4));
  interp->nodes[0] = idmax;

  for ( i=1; i<RBsize; i++ ){
    gsl_vector *interpolant, *subbasis;
    gsl_vector_view subview;

    Vview = gsl_matrix_view_array(V, i, i);

    for ( j=0; j<i; j++ ){
      for ( k=0; k<i; k++ ){
        XLAL_CALLGSL( gsl_matrix_set(&Vview.matrix, k, j, gsl_matrix_get(&RBview.matrix, j, interp->nodes[k])) );
      }
    }

    /* get B matrix */
    REAL8Array *B = B_matrix(&Vview.matrix, &RBview.matrix);
    gsl_matrix_view Bview;
    Bview = gsl_matrix_view_array(B->data, B->dimLength->data[0], B->dimLength->data[1]);

    /* make empirical interpolant of basis */
    XLAL_CALLGSL( interpolant = gsl_vector_calloc(dlength) );
    XLAL_CALLGSL( subbasis = gsl_vector_calloc(i) );
    XLAL_CALLGSL( subview = gsl_matrix_row(&RBview.matrix, i) );

    for ( k=0; k<i; k++ ){
      XLAL_CALLGSL( gsl_vector_set(subbasis, k, gsl_vector_get(&subview.vector, interp->nodes[k])) );
    }

    XLAL_CALLGSL( gsl_blas_dgemv(CblasTrans, 1.0, &Bview.matrix, subbasis, 0., interpolant) );

    /* get residuals of interpolant */
    XLAL_CALLGSL( gsl_vector_sub(interpolant, &subview.vector) );

    XLAL_CALLGSL( newidx = (int)gsl_blas_idamax(interpolant) );

    interp->nodes[i] = newidx;

    XLAL_CALLGSL( gsl_vector_free(subbasis) );
    XLALDestroyREAL8Array( B );
    XLAL_CALLGSL( gsl_vector_free(interpolant) );

    /* reallocate memory for V */
    V = XLALRealloc(V, (i+1)*(i+1)*sizeof(REAL8));
  }

  /* get final B vector with all the indices */
  Vview = gsl_matrix_view_array(V, RBsize, RBsize);
  for( j=0; j<RBsize; j++ ){
    for( k=0; k<RBsize; k++ ){
      XLAL_CALLGSL( gsl_matrix_set(&Vview.matrix, k, j, gsl_matrix_get(&RBview.matrix, j, interp->nodes[k])) );
    }
  }

  /* allocate memory for intpolant array */
  interp->B = B_matrix(&Vview.matrix, &RBview.matrix);

  XLALFree(V);

  return interp;
}


/**
 * \brief Create a complex empirical interpolant from a set of orthonormal basis functions
 *
 * Given a complex \c gsl_matrix_complex containing a set of orthonormal basis functions generate an
 * empirical intopolant, and set of interpolation points, using Algorithm 2 of
 * \cite FGHKT2014 .
 *
 * @param[in] RB The set of basis functions
 *
 * @return A \c LALInferenceCOMPLEXROQInterpolant structure containing the interpolant and its nodes
 */
LALInferenceCOMPLEXROQInterpolant *LALInferenceGenerateCOMPLEXROQInterpolant(COMPLEX16Array *RB){
  XLAL_CHECK_NULL( RB != NULL, XLAL_EFUNC, "Reduced basis array is NULL!" );

  size_t RBsize = RB->dimLength->data[0]; /* reduced basis size (no. of reduced bases) */
  size_t dlength = RB->dimLength->data[1]; /* length of each base */
  size_t i=1, j=0, k=0;
  REAL8 *V = XLALMalloc(sizeof(COMPLEX16));
  gsl_matrix_complex_view Vview;

  gsl_complex scale1, scale0;
  GSL_SET_COMPLEX(&scale1, 1., 0.);
  GSL_SET_COMPLEX(&scale0, 0., 0.);

  LALInferenceCOMPLEXROQInterpolant *interp = XLALMalloc(sizeof(LALInferenceCOMPLEXROQInterpolant));
  int idmax = 0, newidx = 0;

  /* get index of maximum absolute value of first basis */
  gsl_matrix_complex_view RBview;
  RBview = gsl_matrix_complex_view_array((double *)RB->data, RBsize, dlength);
  gsl_vector_complex_view firstbasis = gsl_matrix_complex_row(&RBview.matrix, 0);
  idmax = complex_vector_maxabs_index(&firstbasis.vector);

  interp->nodes = XLALMalloc(RBsize*sizeof(UINT4));
  interp->nodes[0] = idmax;

  for ( i=1; i<RBsize; i++ ){
    gsl_vector_complex *interpolant, *subbasis;
    gsl_vector_complex_view subview;

    Vview = gsl_matrix_complex_view_array(V, i, i);

    for ( j=0; j<i; j++ ){
      for ( k=0; k<i; k++ ){
        XLAL_CALLGSL( gsl_matrix_complex_set(&Vview.matrix, k, j, gsl_matrix_complex_get(&RBview.matrix, j, interp->nodes[k])) );
      }
    }

    /* get B matrix */
    COMPLEX16Array *B = complex_B_matrix(&Vview.matrix, &RBview.matrix);
    gsl_matrix_complex_view Bview;
    Bview = gsl_matrix_complex_view_array((double *)B->data, B->dimLength->data[0], B->dimLength->data[1]);

    /* make empirical interpolant of basis */
    XLAL_CALLGSL( interpolant = gsl_vector_complex_calloc(dlength) );
    XLAL_CALLGSL( subbasis = gsl_vector_complex_calloc(i) );
    XLAL_CALLGSL( subview = gsl_matrix_complex_row(&RBview.matrix, i) );

    for ( k=0; k<i; k++ ){
      XLAL_CALLGSL( gsl_vector_complex_set(subbasis, k, gsl_vector_complex_get(&subview.vector, interp->nodes[k])) );
    }

    XLAL_CALLGSL( gsl_blas_zgemv(CblasTrans, scale1, &Bview.matrix, subbasis, scale0, interpolant) );

    /* get residuals of interpolant */
    XLAL_CALLGSL( gsl_vector_complex_sub(interpolant, &subview.vector) );

    newidx = complex_vector_maxabs_index(interpolant);

    interp->nodes[i] = newidx;

    XLAL_CALLGSL( gsl_vector_complex_free(subbasis) );
    XLALDestroyCOMPLEX16Array( B );
    XLAL_CALLGSL( gsl_vector_complex_free(interpolant) );

    /* reallocate memory for V */
    V = XLALRealloc(V, (i+1)*(i+1)*sizeof(COMPLEX16));
  }

  /* get final B vector with all the indices */
  Vview = gsl_matrix_complex_view_array((double*)V, RBsize, RBsize);
  for( j=0; j<RBsize; j++ ){
    for( k=0; k<RBsize; k++ ){
      XLAL_CALLGSL( gsl_matrix_complex_set(&Vview.matrix, k, j, gsl_matrix_complex_get(&RBview.matrix, j, interp->nodes[k])) );
    }
  }

  /* allocate memory for intpolant array */
  interp->B = complex_B_matrix(&Vview.matrix, &RBview.matrix);

  XLALFree(V);

  return interp;
}


/** \brief Create the weights for the ROQ interpolant for the linear data and model dot product <d|h>
 *
 * @param[in] B The interpolant matrix
 * @param[in] data The real data vector
 * @param[in] vars A vector of data noise variance values (or a single value) to weight the "weights"
 *
 * @return The vector of weights
 */
REAL8Vector *LALInferenceGenerateREAL8LinearWeights(REAL8Array *B, REAL8Vector *data, REAL8Vector *vars){
  REAL8Vector *weights = NULL;
  gsl_vector *datacopy;

  /* get view of the interpolant matrix */
  gsl_matrix_view Bview;
  XLAL_CALLGSL( Bview = gsl_matrix_view_array(B->data, B->dimLength->data[0], B->dimLength->data[1]) );

  /* get view of the data and variances */
  gsl_vector_view dataview, varsview;
  XLAL_CALLGSL( dataview = gsl_vector_view_array(data->data, data->length) );
  XLAL_CALLGSL( varsview = gsl_vector_view_array(vars->data, vars->length) );

  XLAL_CHECK_NULL( vars->length == 1 || vars->length == B->dimLength->data[1], XLAL_EFUNC, "Vector of variance values is the wrong size" );

  XLAL_CALLGSL( datacopy = gsl_vector_alloc(B->dimLength->data[1]) );
  XLAL_CALLGSL( gsl_vector_memcpy(datacopy, &dataview.vector) );

  if ( vars->length == 1 ){ XLAL_CALLGSL( gsl_vector_scale( datacopy, 1./vars->data[0]) ); }
  else{ XLAL_CALLGSL( gsl_vector_div( datacopy, &varsview.vector ) ); }

  /* create weights */
  weights = XLALCreateREAL8Vector( B->dimLength->data[0] );
  gsl_vector_view weightsview;
  XLAL_CALLGSL( weightsview = gsl_vector_view_array(weights->data, weights->length) );
  XLAL_CALLGSL( gsl_blas_dgemv(CblasNoTrans, 1.0, &Bview.matrix, datacopy, 0., &weightsview.vector) );
  XLAL_CALLGSL( gsl_vector_free( datacopy ) );

  return weights;
}


/** \brief Create the weights for the ROQ interpolant for the model quadratic model term real(<h|h>).
 *
 * @param[in] B The interpolant matrix
 * @param[in] vars A vector of data noise variance values (or a single value) to weight the "weights"
 *
 * @return The vector of weights
 */
REAL8Vector *LALInferenceGenerateQuadraticWeights(REAL8Array *B, REAL8Vector *vars){
  REAL8Vector *weights = NULL;
  gsl_vector *vecones = gsl_vector_alloc(B->dimLength->data[1]);
  XLAL_CALLGSL( gsl_vector_set_all(vecones, 1.0) ); /* set all to 1 */

  /* get view of the interpolant matrix */
  gsl_matrix_view Bview;
  XLAL_CALLGSL( Bview = gsl_matrix_view_array(B->data, B->dimLength->data[0], B->dimLength->data[1]) );

  /* get view of the data and variances */
  gsl_vector_view varsview;
  XLAL_CALLGSL( varsview = gsl_vector_view_array(vars->data, vars->length) );

  XLAL_CHECK_NULL( vars->length == 1 || vars->length == B->dimLength->data[1], XLAL_EFUNC, "Vector of variance values is the wrong size" );

  if ( vars->length == 1 ){ XLAL_CALLGSL( gsl_vector_scale( vecones, 1./vars->data[0]) ); }
  else{ XLAL_CALLGSL( gsl_vector_div( vecones, &varsview.vector ) ); }

  /* create weights */
  weights = XLALCreateREAL8Vector( B->dimLength->data[0] );
  gsl_vector_view weightsview;
  XLAL_CALLGSL( weightsview = gsl_vector_view_array(weights->data, weights->length) );
  XLAL_CALLGSL( gsl_blas_dgemv(CblasNoTrans, 1.0, &Bview.matrix, vecones, 0., &weightsview.vector) );
  XLAL_CALLGSL( gsl_vector_free( vecones ) );

  return weights;
}


/** \brief Create the weights for the ROQ interpolant for the complex data and model dot product <d|h>
 *
 * @param[in] B The interpolant matrix
 * @param[in] data The complex data vector
 * @param[in] vars A vector of data noise variance values (or a single value) to weight the "weights"
 *
 * @return The vector of weights
 */
COMPLEX16Vector *LALInferenceGenerateCOMPLEX16LinearWeights(COMPLEX16Array *B, COMPLEX16Vector *data, REAL8Vector *vars){
  COMPLEX16Vector *weights = NULL;
  gsl_vector_complex *conjdata;
  gsl_complex scale1, scale0, cconj;
  size_t i = 0;

  /* get view of the interpolant matrix */
  gsl_matrix_complex_view Bview;
  XLAL_CALLGSL( Bview = gsl_matrix_complex_view_array((double *)B->data, B->dimLength->data[0], B->dimLength->data[1]) );

  XLAL_CHECK_NULL( vars->length == 1 || vars->length == B->dimLength->data[1], XLAL_EFUNC, "Vector of variance values is the wrong size" );

  XLAL_CALLGSL( conjdata = gsl_vector_complex_alloc(B->dimLength->data[1]) );

  /* get conjugate of data and scale it */
  for ( i=0; i<conjdata->size; i++ ){
    gsl_complex cdata;
    GSL_SET_COMPLEX(&cdata, creal(data->data[i]), cimag(data->data[i]));
    XLAL_CALLGSL( cconj = gsl_complex_conjugate(cdata) );
    if ( vars->length == 1 ) { XLAL_CALLGSL( cconj = gsl_complex_div_real( cconj, vars->data[0] ) ); }
    else { XLAL_CALLGSL( cconj = gsl_complex_div_real( cconj, vars->data[i] ) ); }
    XLAL_CALLGSL( gsl_vector_complex_set(conjdata, i, cconj) );
  }

  /* create weights */
  weights = XLALCreateCOMPLEX16Vector( B->dimLength->data[0] );
  gsl_vector_complex_view weightsview;
  XLAL_CALLGSL( weightsview = gsl_vector_complex_view_array((double *)weights->data, weights->length) );
  GSL_SET_COMPLEX(&scale1, 1., 0.);
  GSL_SET_COMPLEX(&scale0, 0., 0.);
  XLAL_CALLGSL( gsl_blas_zgemv(CblasNoTrans, scale1, &Bview.matrix, conjdata, scale0, &weightsview.vector) );
  XLAL_CALLGSL( gsl_vector_complex_free(conjdata) );

  return weights;
}


/** \brief Calculate the dot product of two vectors using the ROQ iterpolant
 *
 * This function calculates the dot product of two real vectors using the ROQ
 * interpolant. This required the interpolant weights computed with \c LALInferenceCreateREAL8DataModelWeights
 * and the waveform model defined at the interpolation node.
 *
 * @param[in] weights The ROQ interpolation weights
 * @param[in] model The waveform model (or quadratic model) defined at the interpolation points
 *
 * @return The dot product \c weights and \c model
 */
REAL8 LALInferenceROQREAL8DotProduct(REAL8Vector *weights, REAL8Vector *model){
  XLAL_CHECK_REAL8( weights->length == model->length, XLAL_EFUNC, "Two vectors must have the same lengths." );

  REAL8 d_dot_m = 0.;
  gsl_vector_view weightsview, modelview;
  XLAL_CALLGSL( weightsview = gsl_vector_view_array(weights->data, weights->length) );
  XLAL_CALLGSL( modelview = gsl_vector_view_array(model->data, model->length) );
  XLAL_CALLGSL( gsl_blas_ddot(&weightsview.vector, &modelview.vector, &d_dot_m) );

  return d_dot_m;
}

/** \brief Calculate the dot product of two complex vectors using the ROQ iterpolant
 *
 * This function calculates the dot product of two complex vectors using the ROQ
 * interpolant. This required the interpolant weights computed with \c LALInferenceCreateCOMPLEX16DataModelWeights
 * and the waveform model defined at the interolation node.
 *
 * @param[in] weights The ROQ interpolation weights
 * @param[in] model The waveform model defined at the interpolation points
 *
 * @return The dot product of the data and model
 */
COMPLEX16 LALInferenceROQCOMPLEX16DotProduct(COMPLEX16Vector *weights, COMPLEX16Vector *model){
  XLAL_CHECK_REAL8( weights->length == model->length, XLAL_EFUNC, "Two vectors must have the same lengths." );

  gsl_complex d_dot_m;
  gsl_vector_complex_view weightsview, modelview;
  XLAL_CALLGSL( weightsview = gsl_vector_complex_view_array((double *)weights->data, weights->length) );
  XLAL_CALLGSL( modelview = gsl_vector_complex_view_array((double *)model->data, model->length) );
  XLAL_CALLGSL( gsl_blas_zdotu(&weightsview.vector, &modelview.vector, &d_dot_m) );

  return GSL_REAL(d_dot_m) + I*GSL_IMAG(d_dot_m);
}


/** \brief Free memory for a \c LALInferenceREALROQInterpolant
 *
 * @param[in] a A pointer to a  \c LALInferenceREALROQInterpolant
 */
void LALInferenceRemoveREALROQInterpolant( LALInferenceREALROQInterpolant *a ){
  if ( a == NULL ) { return; }

  if ( a->B != NULL ){ XLALDestroyREAL8Array( a->B ); }
  if ( a->nodes != NULL ){ XLALFree( a->nodes ); }
  XLALFree( a );
  a = NULL;
}


/** \brief Free memory for a \c LALInferenceCOMPLEXROQInterpolant
 *
 * @param[in] a A pointer to a  \c LALInferenceCOMPLEXROQInterpolant
 */
void LALInferenceRemoveCOMPLEXROQInterpolant( LALInferenceCOMPLEXROQInterpolant *a ){
  if ( a == NULL ) { return; }

  if ( a->B != NULL ){ XLALDestroyCOMPLEX16Array( a->B ); }
  if ( a->nodes != NULL ){ XLALFree( a->nodes ); }
  XLALFree( a );
  a = NULL;
}
