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
#include <lal/XLALGSL.h>

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
REAL8 weighted_dot_product(const gsl_vector *weight, gsl_vector *a, gsl_vector *b);

/* the dot product of two complex vectors scaled by a weighting factor */
gsl_complex complex_weighted_dot_product(const gsl_vector *weight, gsl_vector_complex *a, gsl_vector_complex *b);

void normalise(const gsl_vector *weight, gsl_vector *a);
void complex_normalise(const gsl_vector *weight, gsl_vector_complex *a);

void normalise_training_set(gsl_vector *weight, gsl_matrix *TS);
void complex_normalise_training_set(gsl_vector *weight, gsl_matrix_complex *TS);

double normalisation(const gsl_vector *weight, gsl_vector *a);
double complex_normalisation(const gsl_vector *weight, gsl_vector_complex *a);

void modified_gs_complex(gsl_vector_complex *ru,
                 gsl_vector_complex *ortho_basis,
                 const gsl_matrix_complex *RB_space,
                 const gsl_vector *wQuad,
                 const int dim_RB);

void iterated_modified_gm_complex(gsl_vector_complex *ru,
                                  gsl_vector_complex *ortho_basis,
                                  const gsl_matrix_complex *RB_space,
                                  const gsl_vector *wQuad,
                                  const int dim_RB);

void modified_gs(gsl_vector *ru,
                 gsl_vector *ortho_basis,
                 const gsl_matrix *RB_space,
                 const gsl_vector *wQuad,
                 const int dim_RB);

void iterated_modified_gm(gsl_vector *ru,
                          gsl_vector *ortho_basis,
                          const gsl_matrix *RB_space,
                          const gsl_vector *wQuad,
                          const int dim_RB);

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
REAL8 weighted_dot_product(const gsl_vector *weight, gsl_vector *a, gsl_vector *b){
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
gsl_complex complex_weighted_dot_product(const gsl_vector *weight, gsl_vector_complex *a, gsl_vector_complex *b){
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
void normalise(const gsl_vector *weight, gsl_vector *a){
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
 * @param[in] a The vector to be normalised (this will be change by the function to return the
 * normalised vector).
 */
void complex_normalise(const gsl_vector *weight, gsl_vector_complex *a){
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


/** \brief Get normalisation constant for a real vector
 *
 * @param[in] weight The weighting(s) in the normalisation (e.g. time of frequency step(s) between points)
 * @param[in] a The vector to calculate the normalisation constant for
 *
 * @return The normalisation
 */
double normalisation(const gsl_vector *weight, gsl_vector *a){
  double norm;

  if ( weight->size == 1 ){
    XLAL_CALLGSL( norm = gsl_blas_dnrm2(a) ); /* use GSL normalisation calculation function */
    norm *= sqrt(gsl_vector_get(weight, 0));
  }
  else if ( weight->size == a->size ){
    norm = sqrt(fabs(weighted_dot_product(weight, a, a)));
  }
  else{
    XLAL_ERROR_REAL8( XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors." );
  }

  return norm;
}


/** \brief Get normalisation constant for a complex vector
 *
 * @param[in] weight The weighting(s) in the normalisation (e.g. time of frequency step(s) between points)
 * @param[in] a The vector to calculate the normalisation constant for
 *
 * @return The normalisation
 */
double complex_normalisation(const gsl_vector *weight, gsl_vector_complex *a){
  double norm;

  if ( weight->size == 1 ){
    XLAL_CALLGSL( norm = gsl_blas_dznrm2(a) ); /* use GSL normalisation calculation function */
    norm *= sqrt(gsl_vector_get(weight, 0));
  }
  else if ( weight->size == a->size ){
    norm = sqrt(gsl_complex_abs(complex_weighted_dot_product(weight, a, a)));
  }
  else{
    XLAL_ERROR_REAL8( XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors." );
  }

  return norm;
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
  COMPLEX16Array *B = NULL;
  UINT4Vector *dims = NULL;

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
  XLAL_CALLGSL( gsl_blas_zgemm(CblasTrans, CblasNoTrans, GSL_COMPLEX_ONE, invV, &subRB.matrix, GSL_COMPLEX_ZERO, &Bview.matrix) );

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


/** \brief Modified Gram-Schmidt algorithm for complex data
 *
 * A modified Gram-Schmidt algorithm taken from the
 * <a href="https://bitbucket.org/sfield83/greedycpp">greedycpp</a> code.
 *
 * @param[out] ru A 1-by-(\c dim_RB + 1) complex vector returning a slice of the R matrix 
 * @param[out] ortho_basis A complex vector returning the orthogonal basis vector
 * @param[in] RB_space A complex matrix containing the existing set of orthonormal basis vectors
 * @param[in] wQuad A vector of weights for the inner products
 * @param[in] dim_RB The number of elements in the current \c RB_space
 */
void modified_gs_complex(gsl_vector_complex *ru,
                         gsl_vector_complex *ortho_basis,
                         const gsl_matrix_complex *RB_space,
                         const gsl_vector *wQuad,
                         const int dim_RB){
  INT4 quad_num = RB_space->size2;
  gsl_complex L2_proj;
  gsl_vector_complex *basis;

  basis = gsl_vector_complex_alloc(quad_num);

  /* if not done, R matrix fills up below diagonal with .5 instead of 0 */
  gsl_vector_complex_set_zero(ru); 

  for(INT4 i = 0; i < dim_RB; i++){
    gsl_matrix_complex_get_row(basis, RB_space, i);

    /* ortho_basis = ortho_basis - L2_proj*basis; */
    L2_proj = complex_weighted_dot_product(wQuad, basis, ortho_basis);
    gsl_vector_complex_set(ru, i, L2_proj);
    gsl_vector_complex_scale(basis, L2_proj); /* basis <- basis*L2_proj */
    gsl_vector_complex_sub(ortho_basis, basis); /* ortho <- ortho - basis */
  }

  double nrm = complex_normalisation(wQuad, ortho_basis);
  gsl_complex nrmc;
  GSL_SET_COMPLEX(&nrmc, nrm, 0.0);
  gsl_vector_complex_set(ru, dim_RB, nrmc);

  complex_normalise(wQuad, ortho_basis);

  gsl_vector_complex_free(basis);
}


/** \brief Modified Gram-Schmidt algorithm for real data
 *
 * A modified Gram-Schmidt algorithm taken from the
 * <a href="https://bitbucket.org/sfield83/greedycpp">greedycpp</a> code.
 *
 * @param[out] ru A 1-by-(\c dim_RB + 1) vector returning a slice of the R matrix 
 * @param[out] ortho_basis A vector returning the orthogonal basis vector
 * @param[in] RB_space A matrix containing the existing set of orthonormal basis vectors
 * @param[in] wQuad A vector of weights for the inner products
 * @param[in] dim_RB The number of elements in the current \c RB_space
 */
void modified_gs(gsl_vector *ru,
                 gsl_vector *ortho_basis,
                 const gsl_matrix *RB_space,
                 const gsl_vector *wQuad,
                 const int dim_RB){
  INT4 quad_num = RB_space->size2;
  REAL8 L2_proj;
  gsl_vector *basis;

  basis = gsl_vector_alloc(quad_num);

  /* if not done, R matrix fills up below diagonal with .5 instead of 0 */
  gsl_vector_set_zero(ru); 

  for(INT4 i = 0; i < dim_RB; i++){
    gsl_matrix_get_row(basis, RB_space, i);

    /* ortho_basis = ortho_basis - L2_proj*basis; */
    L2_proj = weighted_dot_product(wQuad, basis, ortho_basis);
    gsl_vector_set(ru, i, L2_proj);
    gsl_vector_scale(basis, L2_proj); /* basis <- basis*L2_proj */
    gsl_vector_sub(ortho_basis, basis); /* ortho <- ortho - basis */
  }

  double nrm = normalisation(wQuad, ortho_basis);
  gsl_vector_set(ru, dim_RB, nrm);

  normalise(wQuad, ortho_basis);

  gsl_vector_free(basis);
}


/** \brief Iterated modified Gram-Schmidt algorithm for complex data
 *
 * An iterated modified Gram-Schmidt algorithm taken from the
 * <a href="https://bitbucket.org/sfield83/greedycpp">greedycpp</a> code. This uses the method
 * given in \cite Hoffmann1989
 *
 * @param[out] ru A complex 1-by-(\c dim_RB + 1) vector returning a slice of the R matrix 
 * @param[out] ortho_basis A complex vector returning the orthogonal basis vector
 * @param[in] RB_space A complex matrix containing the existing set of orthonormal basis vectors
 * @param[in] wQuad A vector of weights for the inner products
 * @param[in] dim_RB The number of elements in the current \c RB_space
 */
void iterated_modified_gm_complex(gsl_vector_complex *ru,
                                  gsl_vector_complex *ortho_basis,
                                  const gsl_matrix_complex *RB_space,
                                  const gsl_vector *wQuad,
                                  const int dim_RB){
  REAL8 ortho_condition = .5; /* hard coded IMGS stopping condition */

  INT4 quad_num = RB_space->size2;
  INT4 r_size = ru->size;
  REAL8 nrm_prev = complex_normalisation(wQuad, ortho_basis);
  UINT4 flag = 0, iter = 0;
  gsl_vector_complex *e,*r_last;
  REAL8 nrm_current = 0.;
  gsl_complex nrmc_current;

  /* allocate memory */
  e = gsl_vector_complex_alloc(quad_num);
  r_last = gsl_vector_complex_alloc(r_size);

  gsl_vector_complex_memcpy(e,ortho_basis);
  complex_normalise(wQuad, e);

  /* if not done, R matrix fills up below diagonal with .5 instead of 0 */
  gsl_vector_complex_set_zero(ru);  

  while( !flag ){
    gsl_vector_complex_memcpy(ortho_basis, e);
    gsl_vector_complex_set_zero(r_last);

    modified_gs_complex(r_last, ortho_basis, RB_space, wQuad, dim_RB);

    gsl_vector_complex_add(ru, r_last);
    nrmc_current = gsl_vector_complex_get(r_last, dim_RB);
    nrm_current = GSL_REAL(nrmc_current);

    gsl_vector_complex_scale(ortho_basis, nrmc_current);

    if( nrm_current/nrm_prev <= ortho_condition ) {
      nrm_prev = nrm_current;
      iter = iter + 1;
      gsl_vector_complex_memcpy(e, ortho_basis);
    }
    else{ flag = 1; }

    nrm_current = complex_normalisation(wQuad, ortho_basis);
    GSL_SET_COMPLEX(&nrmc_current, nrm_current, 0.0);
    gsl_vector_complex_set(ru, dim_RB, nrmc_current);

    complex_normalise(wQuad, ortho_basis);
  }

  gsl_vector_complex_free(e);
  gsl_vector_complex_free(r_last);
}


/** \brief Iterated modified Gram-Schmidt algorithm for real data
 *
 * An iterated modified Gram-Schmidt algorithm taken from the
 * <a href="https://bitbucket.org/sfield83/greedycpp">greedycpp</a> code. This uses the method
 * given in \cite Hoffmann1989
 *
 * @param[out] ru A 1-by-(\c dim_RB + 1) vector returning a slice of the R matrix 
 * @param[out] ortho_basis A vector returning the orthogonal basis vector
 * @param[in] RB_space A matrix containing the existing set of orthonormal basis vectors
 * @param[in] wQuad A vector of weights for the inner products
 * @param[in] dim_RB The number of elements in the current \c RB_space
 */
void iterated_modified_gm(gsl_vector *ru,
                          gsl_vector *ortho_basis,
                          const gsl_matrix *RB_space,
                          const gsl_vector *wQuad,
                          const int dim_RB){
  REAL8 ortho_condition = .5; /* hard coded IMGS stopping condition */

  INT4 quad_num = RB_space->size2;
  INT4 r_size = ru->size;
  REAL8 nrm_prev = normalisation(wQuad, ortho_basis);
  UINT4 flag = 0, iter = 0;
  gsl_vector *e,*r_last;
  REAL8 nrm_current = 0.;

  /* allocate memory */
  e = gsl_vector_alloc(quad_num);
  r_last = gsl_vector_alloc(r_size);

  gsl_vector_memcpy(e,ortho_basis);
  normalise(wQuad, e);

  /* if not done, R matrix fills up below diagonal with .5 instead of 0 */
  gsl_vector_set_zero(ru);  

  while( !flag ){
    gsl_vector_memcpy(ortho_basis, e);
    gsl_vector_set_zero(r_last);

    modified_gs(r_last, ortho_basis, RB_space, wQuad, dim_RB);

    gsl_vector_add(ru, r_last);
    nrm_current = gsl_vector_get(r_last, dim_RB);

    gsl_vector_scale(ortho_basis, nrm_current);

    if( nrm_current/nrm_prev <= ortho_condition ) {
      nrm_prev = nrm_current;
      iter = iter + 1;
      gsl_vector_memcpy(e, ortho_basis);
    }
    else{ flag = 1; }

    nrm_current = normalisation(wQuad, ortho_basis);
    gsl_vector_set(ru, dim_RB, nrm_current);

    normalise(wQuad, ortho_basis);
  }

  gsl_vector_free(e);
  gsl_vector_free(r_last);
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
 * @param[in] TS A \c REAL8Array matrix containing the complex training set, where the number of
 * waveforms in the training set is given by the rows and the waveform points by the columns. This
 * will be modified by this function, as the training set will be normalised.
 * @param[out] greedypoints A \c UINT4Vector to return the indices of the training set rows that
 * have been used to form the reduced basis.
 *
 * @return A \c REAL8 with the maximum projection error for the final reduced basis.
 *
 * \sa LALInferenceEnrichREAL8Basis
 */
REAL8 LALInferenceGenerateREAL8OrthonormalBasis(REAL8Array **RBin,
                                                const REAL8Vector *delta,
                                                REAL8 tolerance,
                                                REAL8Array **TS,
                                                UINT4Vector **greedypoints){
  REAL8Array *ts = NULL;
  ts = *TS; // pointer to training set

  size_t cols = ts->dimLength->data[1], rows = ts->dimLength->data[0];

  /* view of training set */
  gsl_matrix_view TSview;
  XLAL_CALLGSL( TSview = gsl_matrix_view_array((double *)ts->data, rows, cols) );

  size_t max_RB = rows;

  //REAL8 greedy_err[max_RB];           /* approximate error */
  UINT4Vector *gpts = NULL;
  gpts = XLALCreateUINT4Vector(max_RB); /* selected greedy points (row selection) */
  *greedypoints = gpts;

  REAL8 worst_err;          /* errors in greedy sweep */
  UINT4 worst_app = 0;      /* worst error stored */
  REAL8 tmpc;               /* worst error temp */

  gsl_vector *ts_el, *last_rb, *ortho_basis, *ru;
  gsl_matrix *R_matrix;
  REAL8 A_row_norms2[rows];              // || A(i,:) ||^2
  REAL8 projection_norms2[rows];
  REAL8 errors[rows];                    // approximation errors at i^{th} sweep

  REAL8Array *RB = NULL;
  UINT4Vector *dims = NULL;
  gsl_matrix_view RBview;

  /* this memory should be freed here */
  ts_el         = gsl_vector_alloc(cols);
  last_rb       = gsl_vector_alloc(cols);
  ortho_basis   = gsl_vector_alloc(cols);
  ru            = gsl_vector_alloc(max_RB);

  //project_coeff = gsl_matrix_alloc(max_RB,rows);
  REAL8 projection_coeff;
  R_matrix = gsl_matrix_alloc(max_RB, max_RB);

  /* initialise projection norms with zeros */
  for(size_t i=0; i<rows; ++i){ projection_norms2[i] = 0; }

  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );

  /* normalise the training set */
  normalise_training_set(&deltaview.vector, &TSview.matrix);

  /* compute norm of each training space element */
  for(size_t i=0; i<rows; ++i) {
    gsl_matrix_get_row(ts_el, &TSview.matrix, i);
    A_row_norms2[i] = normalisation(&deltaview.vector, ts_el);
  }

  /* initialize algorithm with first training set value */
  dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = 1; /* one row */
  dims->data[1] = cols;
  RB = XLALCreateREAL8Array( dims );
  *RBin = RB;

  XLAL_CALLGSL( RBview = gsl_matrix_view_array((double*)RB->data, 1, cols) );
  gsl_matrix_get_row(ts_el, &TSview.matrix, 0);
  gsl_matrix_set_row(&RBview.matrix, 0, ts_el);
  tmpc = 1.;
  gsl_matrix_set(R_matrix, 0, 0, tmpc);  // assumes normalized solutions

  gpts->data[0] = 0;
  UINT4 dim_RB          = 1;
  //greedy_err[0]    = 1.0;

  /* loop to find reduced basis */
  while( 1 ){
    gsl_matrix_get_row(last_rb, &RBview.matrix, dim_RB-1); /* previous basis */

    /* Compute overlaps of pieces of training set with rb_new */
    for(size_t i = 0; i < rows; i++){
      gsl_matrix_get_row(ts_el, &TSview.matrix, i);
      projection_coeff = weighted_dot_product(&deltaview.vector, last_rb, ts_el);
      projection_norms2[i] += (projection_coeff*projection_coeff);
      errors[i] = A_row_norms2[i] - projection_norms2[i];
    }

    /* find worst represented training set element, and add to basis */
    worst_err = 0.0;
    for(size_t i = 0; i < rows; i++) {
      if(worst_err < errors[i]) {
        worst_err = errors[i];
        worst_app = i;
      }
    }

    // This block exists for cases when the reduced basis is being used during enrichment.
    // It exists to make sure that all new training elements get added to the reduced basis
    // even if their errors are very small.
    if ( tolerance == 0. ){
      // check if worst_app is already in gpts
      UINT4 idxexists = 0;
      for(size_t i = 0; i < dim_RB; i++) {
        if ( worst_app == gpts->data[i] ){
          idxexists = 1;
          break;
        }
      }

      // if it is, then just find the first index that is not in gpts already
      if ( idxexists ){
        UINT4 newidxexists = 0;
        for (size_t i = 0; i < rows; i++){
          newidxexists = 0;
          for (size_t j = 0; j < dim_RB; j++){
            if ( i == gpts->data[j] ){
              newidxexists = 1;
              break;
            }
          }
          if ( !newidxexists ){
            worst_app = i;
            break;
          }
        }
      }
    }
    //////////// end check ////////////

    gpts->data[dim_RB] = worst_app;

    /* add worst approximated solution to basis set */
    gsl_matrix_get_row(ortho_basis, &TSview.matrix, worst_app);
    iterated_modified_gm(ru, ortho_basis, &RBview.matrix, &deltaview.vector, dim_RB); /* use IMGS */

    /* check normalisation of generated orthogonal basis is not NaN (cause by a new orthogonal basis
      having zero residual with the current basis) - if this is the case do not add the new basis. */
    if ( gsl_isnan(gsl_vector_get(ru, dim_RB)) ){ break; }

    /* add to reduced basis */
    dims->data[0] = dim_RB+1; /* add row */
    RB = XLALResizeREAL8Array( RB, dims );

    /* add on next basis */
    XLAL_CALLGSL( RBview = gsl_matrix_view_array((double*)RB->data, dim_RB+1, cols) );

    gsl_matrix_set_row(&RBview.matrix, dim_RB, ortho_basis);
    gsl_matrix_set_row(R_matrix, dim_RB, ru);

    ++dim_RB;

    /* decide if another greedy sweep is needed */
    if( (dim_RB == max_RB) || (worst_err < tolerance) || (rows == dim_RB) ){ break; }
  }

  gpts = XLALResizeUINT4Vector( gpts, dim_RB );

  XLALDestroyUINT4Vector(dims);
  gsl_vector_free(ts_el);
  gsl_vector_free(last_rb);
  gsl_vector_free(ortho_basis);
  gsl_vector_free(ru);
  gsl_matrix_free(R_matrix);

  return worst_err;
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
 * @param[in] TS A \c COMPLEX16Array matrix containing the complex training set, where the number
 * of waveforms in the training set is given by the rows and the waveform points by the columns.
 * This will be modified by this function, as the training set will be normalised.
 * @param[out] greedypoints A \c UINT4Vector to return the indices of the training set rows that
 * have been used to form the reduced basis.
 *
 * @return A \c REAL8 with the maximum projection error for the final reduced basis.
 *
 * \sa LALInferenceEnrichCOMPLEX16Basis
 */
REAL8 LALInferenceGenerateCOMPLEX16OrthonormalBasis(COMPLEX16Array **RBin,
                                                    const REAL8Vector *delta,
                                                    REAL8 tolerance,
                                                    COMPLEX16Array **TS,
                                                    UINT4Vector **greedypoints){
  COMPLEX16Array *ts = NULL;
  ts = *TS; // pointer to training set

  size_t cols = ts->dimLength->data[1], rows = ts->dimLength->data[0];

  /* view of training set */
  gsl_matrix_complex_view TSview;
  XLAL_CALLGSL( TSview = gsl_matrix_complex_view_array((double *)ts->data, rows, cols) );

  size_t max_RB = rows;

  UINT4Vector *gpts = NULL;
  gpts = XLALCreateUINT4Vector(max_RB); /* selected greedy points (row selection) */
  *greedypoints = gpts;

  REAL8 worst_err;          /* errors in greedy sweep */
  UINT4 worst_app = 0;      /* worst error stored */
  gsl_complex tmpc;         /* worst error temp */

  gsl_vector_complex *ts_el, *last_rb, *ortho_basis, *ru;
  gsl_matrix_complex *R_matrix;
  REAL8 A_row_norms2[rows];              // || A(i,:) ||^2
  REAL8 projection_norms2[rows];
  REAL8 errors[rows];                    // approximation errors at i^{th} sweep

  COMPLEX16Array *RB = NULL;
  UINT4Vector *dims = NULL;
  gsl_matrix_complex_view RBview;
  
  /* this memory should be freed here */
  ts_el         = gsl_vector_complex_alloc(cols);
  last_rb       = gsl_vector_complex_alloc(cols);
  ortho_basis   = gsl_vector_complex_alloc(cols);
  ru            = gsl_vector_complex_alloc(max_RB);

  //project_coeff = gsl_matrix_complex_alloc(max_RB,rows);
  gsl_complex projection_coeff;
  R_matrix = gsl_matrix_complex_alloc(max_RB, max_RB);

  /* initialise projection norms with zeros */
  for(size_t i=0; i<rows; ++i){ projection_norms2[i] = 0; }

  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  
  /* normalise the training set */
  complex_normalise_training_set(&deltaview.vector, &TSview.matrix);

  /* compute norm of each training space element */
  for(size_t i=0; i<rows; ++i) {
    gsl_matrix_complex_get_row(ts_el, &TSview.matrix, i);
    A_row_norms2[i] = complex_normalisation(&deltaview.vector, ts_el);
  }

  /* initialize algorithm with first training set value */
  dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = 1; /* one row */
  dims->data[1] = cols;
  RB = XLALCreateCOMPLEX16Array( dims );
  *RBin = RB;

  XLAL_CALLGSL( RBview = gsl_matrix_complex_view_array((double*)RB->data, 1, cols) );
  gsl_matrix_complex_get_row(ts_el, &TSview.matrix, 0);
  gsl_matrix_complex_set_row(&RBview.matrix, 0, ts_el);
  GSL_SET_COMPLEX(&tmpc, 1.0, 0.0);
  gsl_matrix_complex_set(R_matrix, 0, 0, tmpc);  // assumes normalized solutions

  gpts->data[0] = 0;
  UINT4 dim_RB          = 1;

  /* loop to find reduced basis */
  while( 1 ){
    gsl_matrix_complex_get_row(last_rb, &RBview.matrix, dim_RB-1); /* previous basis */

    /* Compute overlaps of pieces of training set with rb_new */
    for(size_t i = 0; i < rows; i++){
      gsl_matrix_complex_get_row(ts_el, &TSview.matrix, i);
      projection_coeff = complex_weighted_dot_product(&deltaview.vector, last_rb, ts_el);
      projection_norms2[i] += (projection_coeff.dat[0]*projection_coeff.dat[0] + projection_coeff.dat[1]*projection_coeff.dat[1]);
      errors[i] = A_row_norms2[i] - projection_norms2[i];
    }

    /* find worst represented training set element, and add to basis */
    worst_err = 0.0;
    for(size_t i = 0; i < rows; i++) {
      if(worst_err < errors[i]) {
        worst_err = errors[i];
        worst_app = i;
      }
    }

    // This block exists for cases when the reduced basis is being used during enrichment.
    // It exists to make sure that all new training elements get added to the reduced basis
    // even if their errors are very small.
    if ( tolerance == 0. ){
      // check if worst_app is already in gpts
      UINT4 idxexists = 0;
      for(size_t i = 0; i < dim_RB; i++) {
        if ( worst_app == gpts->data[i] ){
          idxexists = 1;
          break;
        }
      }

      // if it is, then just find the first index that is not in gpts already
      if ( idxexists ){
        UINT4 newidxexists = 0;
        for (size_t i = 0; i < rows; i++){
          newidxexists = 0;
          for (size_t j = 0; j < dim_RB; j++){
            if ( i == gpts->data[j] ){
              newidxexists = 1;
              break;
            }
          }
          if ( !newidxexists ){
            worst_app = i;
            break;
          }
        }
      }
    }
    //////////// end check ////////////

    gpts->data[dim_RB] = worst_app;

    /* add worst approximated solution to basis set */
    gsl_matrix_complex_get_row(ortho_basis, &TSview.matrix, worst_app);
    iterated_modified_gm_complex(ru, ortho_basis, &RBview.matrix, &deltaview.vector, dim_RB); /* use IMGS */

    /* check normalisation of generated orthogonal basis is not NaN (cause by a new orthogonal basis
      having zero residual with the current basis) - if this is the case do not add the new basis. */
    if ( gsl_isnan(GSL_REAL(gsl_vector_complex_get(ru, dim_RB))) ){ break; }

    /* add to reduced basis */
    dims->data[0] = dim_RB+1; /* add row */
    RB = XLALResizeCOMPLEX16Array( RB, dims );

    /* add on next basis */
    XLAL_CALLGSL( RBview = gsl_matrix_complex_view_array((double*)RB->data, dim_RB+1, cols) );

    gsl_matrix_complex_set_row(&RBview.matrix, dim_RB, ortho_basis);
    gsl_matrix_complex_set_row(R_matrix, dim_RB, ru);

    ++dim_RB;

    /* decide if another greedy sweep is needed */
    if( (dim_RB == max_RB) || (worst_err < tolerance) || (rows == dim_RB) ){ break; }
  }

  gpts = XLALResizeUINT4Vector( gpts, dim_RB );

  XLALDestroyUINT4Vector(dims);
  gsl_vector_complex_free(ts_el);
  gsl_vector_complex_free(last_rb);
  gsl_vector_complex_free(ortho_basis);
  gsl_vector_complex_free(ru);
  gsl_matrix_complex_free(R_matrix);

  return worst_err;
}


/**
 * \brief Validate the real reduced basis against another set of waveforms
 *
 * This function projects a set of waveforms onto the reduced basis and
 * checks that the residuals are within a given tolerance. It returns
 * the projection errors.
 *
 * Note that the projection error returned are the square of the residual
 * errors, as is used as the criterion for adding new bases in the reduced
 * basis generation function \c LALInferenceGenerateREAL8OrthonormalBasis.
 * This is different to the \c validation.cpp code in <a href="https://bitbucket.org/sfield83/greedycpp">greedycpp</a>,
 * which returns the square root of the value we return.
 *
 * @param[out] projerr The projection errors for each test waveform.
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] RB The reduced basis set
 * @param[in] testmodels The set of waveform models to project onto the basis (these will
 * be changed by this function, as the waveforms will get normalised).
 *
 * \sa LALInferenceTestREAL8OrthonormalBasis
 */
void LALInferenceValidateREAL8OrthonormalBasis(REAL8Vector **projerr,
                                               const REAL8Vector *delta,
                                               const REAL8Array *RB,
                                               REAL8Array **testmodels){
  XLAL_CHECK_VOID( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK_VOID( RB != NULL, XLAL_EFUNC, "Reduced basis set array is NULL!" );
  XLAL_CHECK_VOID( RB->dimLength->length == 2, XLAL_EFUNC, "Reduced basis set array must have only two dimensions" );

  REAL8Array *tm = NULL;
  tm = *testmodels;

  size_t dlength = RB->dimLength->data[1], nts = tm->dimLength->data[0];
  size_t k = 0;
  gsl_vector *r_tmp, *testrow;

  /* normalise the test set */
  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  gsl_matrix_view testmodelsview;
  XLAL_CALLGSL( testmodelsview = gsl_matrix_view_array(tm->data, nts, dlength) );
  normalise_training_set(&deltaview.vector, &testmodelsview.matrix);

  gsl_matrix_view RBview;
  RBview = gsl_matrix_view_array(RB->data, RB->dimLength->data[0], dlength);

  XLAL_CALLGSL( testrow = gsl_vector_alloc(dlength) );
  XLAL_CALLGSL( r_tmp = gsl_vector_alloc(RB->dimLength->data[0]) );

  REAL8Vector *pe = NULL;
  pe = XLALCreateREAL8Vector( nts );
  *projerr = pe;

  for ( k = 0; k < nts; k++ ){
    pe->data[k] = 0.;
    REAL8 r_tmp_nrm = 0.;

    XLAL_CALLGSL( gsl_matrix_get_row(testrow, &testmodelsview.matrix, k) );

    REAL8 nrm = normalisation(&deltaview.vector, testrow); // normalisation (should be 1 as test models are normalised)

    // un-scale testrow
    if ( delta->length == 1 ){
      XLAL_CALLGSL( gsl_vector_scale(testrow, delta->data[0]) );
    }
    else if ( testrow->size == deltaview.vector.size ){
      XLAL_CALLGSL( gsl_vector_mul(testrow, &deltaview.vector) );
    }
    else{
      XLAL_ERROR_VOID(XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors.");
    }

    // get projections
    XLAL_CALLGSL( gsl_blas_dgemv( CblasNoTrans, 1., &RBview.matrix, testrow, 0., r_tmp ) );

    r_tmp_nrm = gsl_blas_dnrm2( r_tmp );
    pe->data[k] = nrm - r_tmp_nrm*r_tmp_nrm;

    if ( pe->data[k] < 0. ) { pe->data[k] = 1.0e-16; } // floating point error can trigger this
  }

  XLAL_CALLGSL( gsl_vector_free( testrow ) );
  XLAL_CALLGSL( gsl_vector_free( r_tmp ) );
}


/**
 * \brief Test the reduced basis against another set of waveforms
 *
 * This function projects a set of waveforms onto the reduced basis and
 * checks that the residuals are within a given tolerance
 *
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] tolerance The allowed (squared) residual tolerence for the test waveforms
 * @param[in] RB The reduced basis set
 * @param[in] testmodels The set of waveform models to project onto the basis
 *
 * @return Returns \c XLAL_SUCCESS if all test waveforms meet the tolerance
 *
 * \sa LALInferenceValidateREAL8OrthonormalBasis
 */
INT4 LALInferenceTestREAL8OrthonormalBasis(const REAL8Vector *delta,
                                           REAL8 tolerance,
                                           const REAL8Array *RB,
                                           REAL8Array **testmodels){
  XLAL_CHECK( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK( RB != NULL, XLAL_EFUNC, "Reduced basis set array is NULL!" );
  XLAL_CHECK( RB->dimLength->length == 2, XLAL_EFUNC, "Reduced basis set array must have only two dimensions" );
  XLAL_CHECK( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );

  REAL8Array *tm = NULL;
  tm = *testmodels;

  size_t nts = tm->dimLength->data[0], k = 0;

  REAL8Vector *projerr = NULL;

  LALInferenceValidateREAL8OrthonormalBasis(&projerr, delta, RB, testmodels);

  for ( k = 0; k < nts; k++ ){
    /* check projection error against tolerance */
    if ( projerr->data[k] > tolerance ) {
      XLALDestroyREAL8Vector( projerr );
      return XLAL_FAILURE;
    }
  }

  XLALDestroyREAL8Vector( projerr );

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
 * @param[in] greedypoints The rows in \c testmodels that formed the reduced basis
 * @param[in] testmodels The set of waveform models used to produce the reduced basis (already normalised)
 * @param[in] testmodelsnew A new set of waveforms to project onto the current reduced basis (these will
 * be changed by this function, as the waveforms will get normalised when pass to
 * \c LALInferenceValidateREAL8OrthonormalBasis, and editted to contain the new set of training waveforms
 * from which the enriched basis was formed).
 *
 * @return Returns the maximum projection error of the new basis set
 */
REAL8 LALInferenceEnrichREAL8Basis(const REAL8Vector *delta,
                                   REAL8 tolerance,
                                   REAL8Array **RB,
                                   UINT4Vector **greedypoints,
                                   const REAL8Array *testmodels,
                                   REAL8Array **testmodelsnew){
  XLAL_CHECK( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );

  size_t dlength = testmodels->dimLength->data[1]; // length of training vectors
  size_t k = 0;

  REAL8 maxprojerr = 0.;

  REAL8Array *RBin = NULL;
  RBin = *RB;
  UINT4Vector *gp = NULL;
  gp = *greedypoints;
  REAL8Array *tm = NULL;
  tm = *testmodelsnew;

  size_t nts = tm->dimLength->data[0]; // number of new training vectors

  XLAL_CHECK_REAL8( dlength == tm->dimLength->data[1], XLAL_EFUNC, "New training set contains waveforms of different length to the original set!" );
  
  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  
  // get projection errors of new test model set against the reduced basis
  REAL8Vector *projerr = NULL;
  LALInferenceValidateREAL8OrthonormalBasis( &projerr, delta, *RB, testmodelsnew );

  // check for errors outside the tolerance value, and reduced new training set accordingly
  size_t keepcounter = 0;
  gsl_matrix_view tmview;
  tmview = gsl_matrix_view_array(tm->data, nts, dlength);
  for ( k = 0; k < projerr->length; k++ ){
    if ( projerr->data[k] > tolerance ){
      if ( keepcounter != k ){
        gsl_vector_view row;
        row = gsl_matrix_row( &tmview.matrix, k );

        // un-scale row (so that it works with the reduced basis generation function below)
        if ( delta->length == 1 ){
          XLAL_CALLGSL( gsl_vector_scale(&row.vector, delta->data[0]) );
        }
        else if ( row.vector.size == delta->length ){
          XLAL_CALLGSL( gsl_vector_mul(&row.vector, &deltaview.vector) );
        }
        else{
          XLAL_ERROR_REAL8(XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors.");
        }
        
        gsl_matrix_set_row( &tmview.matrix, keepcounter, &row.vector ); // copy row to earlier in the matrix
      }
      keepcounter++;
    }
  }

  XLALDestroyREAL8Vector( projerr );

  if ( keepcounter == 0 ){ return maxprojerr; } // nothing was above the tolerance

  // add vectors used to produce the original reduced basis
  UINT4Vector *dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = RBin->dimLength->data[0] + keepcounter;
  dims->data[1] = dlength;
  tm = XLALResizeREAL8Array( tm, dims );
  tmview = gsl_matrix_view_array(tm->data, dims->data[0], dlength);

  gsl_matrix_view otmview;
  otmview = gsl_matrix_view_array(testmodels->data, testmodels->dimLength->data[0], dlength); // view of original training set
  for ( k = 0; k < RBin->dimLength->data[0]; k++ ){
    gsl_vector_view row;
    row = gsl_matrix_row( &otmview.matrix, gp->data[k] );

    // un-scale row (so that it works with the reduced basis generation function below)
    if ( delta->length == 1 ){
      XLAL_CALLGSL( gsl_vector_scale(&row.vector, delta->data[0]) );
    }
    else if ( row.vector.size == delta->length ){
      XLAL_CALLGSL( gsl_vector_mul(&row.vector, &deltaview.vector) );
    }
    else{
      XLAL_ERROR_REAL8(XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors.");
    }

    gsl_matrix_set_row( &tmview.matrix, keepcounter + k, &row.vector );
  }

  /* recalculate the reduced basis using the updated training set */
  XLALDestroyREAL8Array( *RB ); /* remove current basis, so it is rebuilt using the entire new training set */
  *RB = NULL;
  XLALDestroyUINT4Vector( *greedypoints );
  *greedypoints = NULL;
  maxprojerr = LALInferenceGenerateREAL8OrthonormalBasis(RB, delta, 0.0, &tm, greedypoints); // set tolerance to 0, so that points are always added

  XLALDestroyUINT4Vector( dims );

  return maxprojerr;
}


/**
 * \brief Validate the complex reduced basis against another set of waveforms
 *
 * This function projects a set of waveforms onto the reduced basis and
 * checks that the residuals are within a given tolerance. It returns
 * the projection errors.
 *
 * Note that the projection error returned are the square of the residual
 * errors, as is used as the criterion for adding new bases in the reduced
 * basis generation function \c LALInferenceGenerateCOMPLEX16OrthonormalBasis.
 * This is different to the \c validation.cpp code in <a href="https://bitbucket.org/sfield83/greedycpp">greedycpp</a>,
 * which returns the square root of the value we return.
 *
 * @param[out] projerr The projection errors (square of the residual) for each test waveform.
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] RB The reduced basis set
 * @param[in] testmodels The set of waveform models to project onto the basis (these will
 * be changed by this function, as the waveforms will get normalised).
 *
 * \sa LALInferenceTestCOMPLEX16OrthonormalBasis
 */
void LALInferenceValidateCOMPLEX16OrthonormalBasis(REAL8Vector **projerr,
                                                   const REAL8Vector *delta,
                                                   const COMPLEX16Array *RB,
                                                   COMPLEX16Array **testmodels){
  XLAL_CHECK_VOID( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK_VOID( RB != NULL, XLAL_EFUNC, "Reduced basis is NULL!" );
  XLAL_CHECK_VOID( RB->dimLength->length == 2, XLAL_EFUNC, "Reduced basis set does not have 2 dimensions!" );

  COMPLEX16Array *tm = NULL;
  tm = *testmodels;

  size_t dlength = RB->dimLength->data[1], nts = tm->dimLength->data[0];
  size_t k = 0, j = 0;
  gsl_vector_complex *r_tmp, *testrow;

  /* normalise the test set */
  gsl_vector_view deltaview;
  XLAL_CALLGSL( deltaview = gsl_vector_view_array(delta->data, delta->length) );
  gsl_matrix_complex_view testmodelsview;
  XLAL_CALLGSL( testmodelsview = gsl_matrix_complex_view_array((double *)tm->data, nts, dlength) );
  complex_normalise_training_set(&deltaview.vector, &testmodelsview.matrix);

  gsl_matrix_complex_view RBview;
  RBview = gsl_matrix_complex_view_array((double *)RB->data, RB->dimLength->data[0], dlength);

  XLAL_CALLGSL( testrow = gsl_vector_complex_alloc(dlength) );
  XLAL_CALLGSL( r_tmp = gsl_vector_complex_alloc(RB->dimLength->data[0]) );

  REAL8Vector *pe = NULL;
  pe = XLALCreateREAL8Vector( nts );
  *projerr = pe;

  /* get projection errors for each test model */
  for ( k = 0; k < nts; k++ ){
    pe->data[k] = 0.;
    REAL8 r_tmp_nrm = 0.;

    XLAL_CALLGSL( gsl_matrix_complex_get_row(testrow, &testmodelsview.matrix, k) );

    REAL8 nrm = complex_normalisation(&deltaview.vector, testrow); // normalisation (should be 1 as test models are normalised)

    // un-scale testrow (gsl_vector_complex_mul doesn't work when delta is a real vector)
    if ( delta->length == 1 || testrow->size == deltaview.vector.size ){
      for ( j = 0; j < dlength; j++ ){
        gsl_complex ctmp;
        XLAL_CALLGSL( ctmp = gsl_vector_complex_get( testrow, j ) );
        if ( delta->length == 1 ){
          XLAL_CALLGSL( ctmp = gsl_complex_mul_real( ctmp, delta->data[0] ) );
        }
        else{
          XLAL_CALLGSL( ctmp = gsl_complex_mul_real( ctmp, delta->data[j] ) );
        }
        XLAL_CALLGSL( gsl_vector_complex_set( testrow, j, ctmp ) );
      }
    }
    else{
      XLAL_ERROR_VOID( XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors." );
    }
    
    // get complex conjugate
    for ( j = 0; j < testrow->size; j++ ){ testrow->data[2*j+1] = -testrow->data[2*j+1]; }

    // get projections
    XLAL_CALLGSL( gsl_blas_zgemv( CblasNoTrans, GSL_COMPLEX_ONE, &RBview.matrix, testrow, GSL_COMPLEX_ZERO, r_tmp ) );

    r_tmp_nrm = gsl_blas_dznrm2(r_tmp);
    pe->data[k] = nrm - r_tmp_nrm*r_tmp_nrm;

    if ( pe->data[k] < 0. ) { pe->data[k] = 1.0e-16; } // floating point error can trigger this
  }

  XLAL_CALLGSL( gsl_vector_complex_free( testrow ) );
  XLAL_CALLGSL( gsl_vector_complex_free( r_tmp ) );
}


/**
 * \brief Test the reduced basis against another set of waveforms
 *
 * This function projects a set of waveforms onto the reduced basis and
 * checks that the residuals are within a given tolerance
 *
 * @param[in] delta The time/frequency step(s) in the training set used to normalise the models.
 * This can be a vector containing just one value.
 * @param[in] tolerance The allowed residual (squared) tolerence for the test waveforms
 * @param[in] RB The reduced basis set
 * @param[in] testmodels The set of waveform models to project onto the basis
 *
 * @return Returns \c XLAL_SUCCESS if all test waveforms meet the tolerance
 *
 * \sa LALInferenceValidateCOMPLEX16OrthonormalBasis
 */
INT4 LALInferenceTestCOMPLEX16OrthonormalBasis(const REAL8Vector *delta,
                                               REAL8 tolerance,
                                               const COMPLEX16Array *RB,
                                               COMPLEX16Array **testmodels){
  XLAL_CHECK( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK( RB != NULL, XLAL_EFUNC, "Reduced basis is NULL!" );
  XLAL_CHECK( RB->dimLength->length == 2, XLAL_EFUNC, "Reduced basis set does not have 2 dimensions!" );
  XLAL_CHECK( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );

  COMPLEX16Array *tm = NULL;
  tm = *testmodels;

  size_t nts = tm->dimLength->data[0], k = 0;

  REAL8Vector *projerr = NULL;

  LALInferenceValidateCOMPLEX16OrthonormalBasis(&projerr, delta, RB, testmodels);

  for ( k = 0; k < nts; k++ ){
    /* check projection error against tolerance */
    if ( projerr->data[k] > tolerance ) {
      XLALDestroyREAL8Vector( projerr );
      return XLAL_FAILURE;
    }
  }

  XLALDestroyREAL8Vector( projerr );

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
 * @param[in] greedypoints The rows in \c testmodels that formed the reduced basis
 * @param[in] testmodels The set of waveform models used to produce the reduced basis
 * @param[in] testmodelsnew A new set of waveforms to project onto the current reduced basis (these will
 * be changed by this function, as the waveforms will get normalised when pass to
 * \c LALInferenceValidateCOMPLEX16OrthonormalBasis, and editted to contain the new set of training waveforms
 * from which the enriched basis was formed)
 *
 * @return Returns the maximum projection error of the new basis set
 */
REAL8 LALInferenceEnrichCOMPLEX16Basis(const REAL8Vector *delta,
                                       REAL8 tolerance,
                                       COMPLEX16Array **RB,
                                       UINT4Vector **greedypoints,
                                       const COMPLEX16Array *testmodels,
                                       COMPLEX16Array **testmodelsnew){
  XLAL_CHECK( delta != NULL, XLAL_EFUNC, "Vector of 'delta' values is NULL!" );
  XLAL_CHECK( tolerance > 0, XLAL_EFUNC, "Tolerance is less than, or equal to, zero!" );

  size_t dlength = testmodels->dimLength->data[1]; // length of training vectors
  size_t k = 0, j = 0;

  REAL8 maxprojerr = 0.;

  COMPLEX16Array *RBin = NULL;
  RBin = *RB;
  UINT4Vector *gp = NULL;
  gp = *greedypoints;
  COMPLEX16Array *tm = NULL;
  tm = *testmodelsnew;

  size_t nts = tm->dimLength->data[0]; // number of new training vectors

  XLAL_CHECK_REAL8( dlength == tm->dimLength->data[1], XLAL_EFUNC, "New training set contains waveforms of different length to the original set!" );

  // get projection errors of new test model set against the reduced basis
  REAL8Vector *projerr = NULL;
  LALInferenceValidateCOMPLEX16OrthonormalBasis( &projerr, delta, *RB, testmodelsnew );

  // check for errors outside the tolerance value, and reduced new training set accordingly
  size_t keepcounter = 0;
  gsl_matrix_complex_view tmview;
  tmview = gsl_matrix_complex_view_array((double *)tm->data, nts, dlength);
  for ( k = 0; k < projerr->length; k++ ){
    if ( projerr->data[k] > tolerance ){
      if ( keepcounter != k ){
        gsl_vector_complex_view row;
        row = gsl_matrix_complex_row( &tmview.matrix, k );

        // un-scale row (so that this works with the RB generation function later)
        if ( delta->length == 1 || row.vector.size == delta->length ){
          for ( j = 0; j < row.vector.size; j++ ){
            gsl_complex ctmp;
            XLAL_CALLGSL( ctmp = gsl_vector_complex_get( &row.vector, j ) );
            if ( delta->length == 1 ){
              XLAL_CALLGSL( ctmp = gsl_complex_mul_real( ctmp, delta->data[0] ) );
            }
            else{
              XLAL_CALLGSL( ctmp = gsl_complex_mul_real( ctmp, delta->data[j] ) );
            }
            XLAL_CALLGSL( gsl_vector_complex_set( &row.vector, j, ctmp ) );
          }
        }
        else{
          XLAL_ERROR_REAL8( XLAL_EFUNC, "Vector of weights must either contain a single value, or be the same length as the other input vectors." );
        }
        gsl_matrix_complex_set_row( &tmview.matrix, keepcounter, &row.vector ); // copy row to earlier in the matrix
      }
      keepcounter++;
    }
  }

  XLALDestroyREAL8Vector( projerr );

  if ( keepcounter == 0 ){ return maxprojerr; } // nothing was above the tolerance

  // add vectors used to produce the orgial reduced basis
  UINT4Vector *dims = XLALCreateUINT4Vector( 2 );
  dims->data[0] = RBin->dimLength->data[0] + keepcounter;
  dims->data[1] = dlength;
  tm = XLALResizeCOMPLEX16Array( tm, dims );
  tmview = gsl_matrix_complex_view_array((double *)tm->data, dims->data[0], dlength);

  gsl_matrix_complex_view otmview;
  otmview = gsl_matrix_complex_view_array((double *)testmodels->data, nts, dlength); // view of original training set
  for ( k = 0; k < RBin->dimLength->data[0]; k++ ){
    gsl_vector_complex_view row;
    row = gsl_matrix_complex_row( &otmview.matrix, gp->data[k] );

    // un-scale row (so that this works with the RB generation function later)
    for ( j = 0; j < row.vector.size; j++ ){
      gsl_complex ctmp;
      XLAL_CALLGSL( ctmp = gsl_vector_complex_get( &row.vector, j ) );
      if ( delta->length == 1 ){
        XLAL_CALLGSL( ctmp = gsl_complex_mul_real( ctmp, delta->data[0] ) );
      }
      else{
        XLAL_CALLGSL( ctmp = gsl_complex_mul_real( ctmp, delta->data[j] ) );
      }
      XLAL_CALLGSL( gsl_vector_complex_set( &row.vector, j, ctmp ) );
    }

    gsl_matrix_complex_set_row( &tmview.matrix, keepcounter + k, &row.vector );
  }

  /* recalculate the reduced basis using the updated training set */
  XLALDestroyCOMPLEX16Array( *RB ); /* remove current basis, so it is rebuilt using the entire new training set */
  *RB = NULL;
  XLALDestroyUINT4Vector( *greedypoints );
  *greedypoints = NULL;
  maxprojerr = LALInferenceGenerateCOMPLEX16OrthonormalBasis(RB, delta, 0.0, &tm, greedypoints); // set tolerance to 0, so that points are always added

  XLALDestroyUINT4Vector( dims );

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

    XLAL_CALLGSL( gsl_blas_zgemv(CblasTrans, GSL_COMPLEX_ONE, &Bview.matrix, subbasis, GSL_COMPLEX_ZERO, interpolant) );

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
  gsl_complex cconj;
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
  XLAL_CALLGSL( gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, &Bview.matrix, conjdata, GSL_COMPLEX_ZERO, &weightsview.vector) );
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
