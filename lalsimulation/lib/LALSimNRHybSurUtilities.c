/*  Copyright (C) 2018 Vijay Varma
 *  Utility functions that are useful for evaluating NRHybrid surrogates.
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/**
 * \author Vijay Varma
 *
 * \file
 *
 * \brief Utilities needed for aligned-spin NR-hybrid surrogate models.
 *
 * Called from:
 * LALSimIMRNRHybSur3dq8.h
 */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


#include "LALSimNRHybSurUtilities.h"

#include <assert.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex_math.h>

#include <lal/SeqFactories.h>
#include <lal/H5FileIO.h>

#include <lal/XLALError.h>
#include <lal/LALConstants.h>
#include <lal/LALSimIMR.h>

#include "LALSimBlackHoleRingdown.h"
#include "LALSimIMRSEOBNRROMUtilities.c"


//*************************************************************************/
//************************* function definitions **************************/
//*************************************************************************/




//********************* H5 wrapper functions ************************/

/**
 * Reads a REAL8 value from a H5 file/group.
 */
int ReadHDF5DoubleDataset(
    REAL8 *param,      /**< Output: REAL8 value from H5 file/group. */
    LALH5File *sub,     /**< H5 file or group. */
    const char *name    /**< Name of REAL8 dataset within file/group. */
) {
    REAL8Vector *data = XLALH5FileReadREAL8Vector(sub, name);

    // Expecting a single double
    if(data == NULL || data->length != 1) {
        XLALDestroyREAL8Vector(data);
        XLAL_ERROR(XLAL_EFUNC,
            "Failed to load `%s' scalar dataset\n", name);
    }

    *param = data->data[0];
    XLALDestroyREAL8Vector(data);

    return XLAL_SUCCESS;
}


/**
 * Reads an INT8 value from a H5 file/group.
 */
int ReadHDF5IntDataset(
    int *param,         /**< Output: int value from H5 file/group. */
    LALH5File *sub,     /**< H5 file or group. */
    const char *name    /**< Name of int dataset within file/group. */
) {
    INT8Vector *data = XLALH5FileReadINT8Vector(sub, name);

    // Expecting a single int
    if(data == NULL || data->length != 1) {
        XLALDestroyINT8Vector(data);
        XLAL_ERROR(XLAL_EFUNC,
            "Failed to load `%s' scalar dataset\n", name);
    }

    *param = (int)data->data[0];
    XLALDestroyINT8Vector(data);

    return XLAL_SUCCESS;
}





//********************* Surrogate loading functions ************************/
// These should only be called once, when initializing the surrogate

/**
 * Loads a single waveform data piece from a H5 file.
 */
static int NRHybSur_LoadDataPiece(
    DataPiece **data_piece,   /**< Output: Waveform data piece. *data_piece
                                should be NULL. Space will be allocated. */
    LALH5File *file,          /**< Opened HDF5 file. */
    const char *sub_grp_name  /**< H5 group name. */
) {

    if (data_piece == NULL || *data_piece != NULL) {
        XLAL_ERROR(XLAL_EFAULT, "*data_piece should be NULL");
    }
    if (file == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "file should not be NULL");
    }

    // Open h5 group
    LALH5File *sub;
    sub = XLALH5GroupOpen(file, sub_grp_name);
    *data_piece = XLALMalloc(sizeof(DataPiece));

    gsl_matrix *ei_basis = NULL;
    int ret = ReadHDF5RealMatrixDataset(sub, "ei_basis", &ei_basis);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load ei_basis.");
    }

    (*data_piece)->ei_basis = ei_basis;

    // Get number of empirical time nodes
    int n_nodes;
    ret = ReadHDF5IntDataset(&n_nodes, sub, "n_nodes");
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load n_nodes.");
    }

    (*data_piece)->n_nodes = n_nodes;
    (*data_piece)->fit_data = XLALMalloc( n_nodes * sizeof(NRHybSurFitData *) );

    // Load fit data for each empirical time node
    const size_t str_size = 20;
    char *node_name = XLALMalloc(str_size);
    UNUSED size_t nwritten;
    for (int i=0; i<n_nodes; i++) {

        nwritten = snprintf(node_name, str_size, "node_num_%d", i);
        assert(nwritten < str_size);
        LALH5File *node_function = XLALH5GroupOpen(sub, node_name);
        NRHybSurFitData *fit_data = XLALMalloc(sizeof(NRHybSurFitData));

        GPRHyperParams *hyperparams = XLALMalloc(sizeof(GPRHyperParams));

        // Load scalars needed for fit
        ret = ReadHDF5DoubleDataset(&(hyperparams->constant_value),
                node_function, "constant_value");
        if (ret != XLAL_SUCCESS) {
            XLALFree(node_name);
            XLAL_ERROR(XLAL_EFUNC, "Failed to load constant_value.");
        }

        ret = ReadHDF5DoubleDataset(&(hyperparams->y_train_mean),
                node_function, "y_train_mean");
        if (ret != XLAL_SUCCESS) {
            XLALFree(node_name);
            XLAL_ERROR(XLAL_EFUNC, "Failed to load y_train_mean.");
        }

        ret = ReadHDF5DoubleDataset(&(fit_data->data_mean),
                node_function, "data_mean");
        if (ret != XLAL_SUCCESS) {
            XLALFree(node_name);
            XLAL_ERROR(XLAL_EFUNC, "Failed to load data_mean.");
        }

        ret = ReadHDF5DoubleDataset(&(fit_data->data_std),
                node_function, "data_std");
        if (ret != XLAL_SUCCESS) {
            XLALFree(node_name);
            XLAL_ERROR(XLAL_EFUNC, "Failed to load data_std.");
        }

        ret = ReadHDF5DoubleDataset(&(fit_data->lin_intercept),
                node_function, "lin_intercept");
        if (ret != XLAL_SUCCESS) {
            XLALFree(node_name);
            XLAL_ERROR(XLAL_EFUNC, "Failed to load lin_intercept.");
        }


        // Load arrays needed for fit
        hyperparams->length_scale = NULL;
        ret = ReadHDF5RealVectorDataset(node_function, "length_scale",
                &(hyperparams->length_scale));
        if (ret != XLAL_SUCCESS) {
            XLALFree(node_name);
            XLAL_ERROR(XLAL_EFUNC, "Failed to load length_scale.");
        }

        hyperparams->alpha = NULL;
        ret = ReadHDF5RealVectorDataset(node_function, "alpha",
                &(hyperparams->alpha));
        if (ret != XLAL_SUCCESS) {
            XLALFree(node_name);
            XLAL_ERROR(XLAL_EFUNC, "Failed to load alpha.");
        }

        fit_data->lin_coef = NULL;
        ret = ReadHDF5RealVectorDataset(node_function, "lin_coef",
                &(fit_data->lin_coef));
        if (ret != XLAL_SUCCESS) {
            XLALFree(node_name);
            XLAL_ERROR(XLAL_EFUNC, "Failed to load lin_coef.");
        }

        XLALH5FileClose(node_function);

        fit_data->hyperparams = hyperparams;
        (*data_piece)->fit_data[i] = fit_data;
    }

    XLALFree(node_name);
    XLALH5FileClose(sub);

    return ret;
}

/**
 * Loads all data pieces of a single waveform mode.
 */
static int NRHybSur_LoadSingleModeData(
    ModeDataPieces **mode_data_pieces, /**< Output: Waveform data pieces of a
                                        given mode. Space will be allocated to
                                        **mode_data_pieces. */
    LALH5File *file,                  /**< Opened HDF5 file. */
    const int mode_idx,               /**< Index corresponding to the mode. */
    const gsl_matrix_long *mode_list  /**< List of all modes. */
) {

    if (mode_data_pieces == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "mode_data_pieces should not be NULL");
    }
    if (file == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "file should not be NULL");
    }

    int ret;
    const UINT4 ell = gsl_matrix_long_get(mode_list, mode_idx, 0);
    const UINT4 m = gsl_matrix_long_get(mode_list, mode_idx, 1);

    ModeDataPieces *data_pieces
        = XLALMalloc(sizeof(*mode_data_pieces[mode_idx]));

    data_pieces->ell = ell;
    data_pieces->m = m;

    // At most two of these data pieces are needed for each mode,
    // so let's set them to NULL by default. This will also raise
    // an error if the wrong data piece is used for a mode
    data_pieces->ampl_data_piece = NULL;
    data_pieces->phase_res_data_piece = NULL;
    data_pieces->coorb_re_data_piece = NULL;
    data_pieces->coorb_im_data_piece = NULL;

    const size_t str_size = 20;
    char *sub_grp_name = XLALMalloc(str_size);

    UNUSED size_t nwritten;
    if (ell == 2 && m ==2){

        // Amplitude of 22 mode
        nwritten = snprintf(sub_grp_name, str_size, "l%u_m%u/amp", ell, m);
        assert(nwritten < str_size);
        ret = NRHybSur_LoadDataPiece(&(data_pieces)->ampl_data_piece,
                file, sub_grp_name);
        if(ret != XLAL_SUCCESS) {
            XLALFree(sub_grp_name);
            XLAL_ERROR(XLAL_EFUNC,
                "Failed to load `%s' data piece", sub_grp_name);
        }

        // phase of 22 mode
        nwritten = snprintf(sub_grp_name, str_size, "l%u_m%u/phase", ell, m);
        assert(nwritten < str_size);
        ret = NRHybSur_LoadDataPiece(&(data_pieces)->phase_res_data_piece,
                file, sub_grp_name);
        if(ret != XLAL_SUCCESS) {
            XLALFree(sub_grp_name);
            XLAL_ERROR(XLAL_EFUNC,
                "Failed to load `%s' data piece", sub_grp_name);
        }

    } else {
        // For m=0, l=even, the imaginary part is zero, so we only need to
        // load the real part. But when m!=0, we still want to load the real
        // part.
        if (m != 0 || ell % 2 == 0) {
            // Real part of coorbital frame mode
            nwritten = snprintf(sub_grp_name, str_size, "l%u_m%u/re", ell, m);
            assert(nwritten < str_size);
            ret = NRHybSur_LoadDataPiece(&(data_pieces)->coorb_re_data_piece,
                    file, sub_grp_name);
            if(ret != XLAL_SUCCESS) {
                XLALFree(sub_grp_name);
                XLAL_ERROR(XLAL_EFUNC,
                    "Failed to load `%s' data piece", sub_grp_name);
            }
        }

        // Similarly, for m=0, l=odd, the real part is zero, so we only need
        // to load the imaginary part. But when m!=0, we still want to load
        // the imaginary part.
        if (m != 0 || ell % 2 == 1) {
            // Imaginary part of coorbital frame mode
            nwritten = snprintf(sub_grp_name, str_size, "l%u_m%u/im", ell, m);
            assert(nwritten < str_size);
            ret = NRHybSur_LoadDataPiece(&(data_pieces)->coorb_im_data_piece,
                    file, sub_grp_name);
            if(ret != XLAL_SUCCESS) {
                XLALFree(sub_grp_name);
                XLAL_ERROR(XLAL_EFUNC,
                    "Failed to load `%s' data piece", sub_grp_name);
            }
        }
    }

    mode_data_pieces[mode_idx] = data_pieces;
    XLALFree(sub_grp_name);
    return ret;
}


/**
 * Initialize a NRHybSurData structure from an open H5 file.
 * This will typically only be called once.
 * For example from NRHybSur3dq8_Init_LALDATA.
 *
 * The HDF5 file should have the following structure:
 *
 * - Top level: \n
 * 'GPR_X_train': h5 dataset, shape=(N,dim) where N=size of training dataset. \n
 * 'domain': h5 dataset, size = number of time samples, can be sparse. \n
 * 'TaylorT3_t_ref': h5 dataset, INT8, used in TaylorT3 piece. \n
 * 'mode_list': h5 dataset, INT8, shape=(M,2), where M=num. of modes. Each
 *      element is (ell,m) of that mode. Only m>=0 modes are modeled. First
 *      mode should always be (2,2). \n
 * 'phaseAlignIdx: h5 dataset, INT8. Needed in TaylorT3 piece. \n
 * 'l2_m0': h5 group, surrogate data needed for this mode. \n
 * 'l2_m1': \n
 * 'l2_m2': \n
 *
 *  - Subgroups for each mode: \n
 *    - 'l2_m2': \n
 *      'amp': h5 group, surrogate data needed for amplitude data piece. \n
 *      'phase' : h5 group, surrogate data needed for phase data piece. \n
 *
 *    - 'l2_m1' and all other modes: \n
 *      're': h5 group, data needed for real part of coorbital frame
 *          waveform data piece. \n
 *      'im': h5 group, data needed for imag part of coorbital frame
 *          waveform data piece. \n
 *    - Subgroups for each data piece: \n
 *      'ei_basis': h5 dataset, shape=(n_nodes, L), where L=len(domain). \n
 *      'n_nodes': h5 dataset, INT8. Num. of empirical time nodes for this data
 *           piece. \n
 *      'node_num_0': h5 group, fit data for this node. \n
 *      'node_num_1': \n
 *           - Subgroups for each node: \n
 *      'alpha': h5 dataset, size=N, where N=size of training dataset. \n
 *      'constant_value': h5 dataset, REAL8, constant value in kernel. \n
 *      'data_mean': h5 dataset, REAL8. Mean of fit data. \n
 *      'data_std': h5 dataset, REAL8. Standard deviation of fit data. \n
 *      'length_scale': h5 dataset, size=dom, where dim=dimensions of the model.
 *           Length scale parameters in the kernel. \n
 *      'lin_coef':  h5 dataset, size=dim. Coefs of of linear fit. \n
 *      'lin_intercept': h5 dataset, REAL8. Intercept of linear fit. \n
 *      'L': h5 dataset, shape=(N,N), where N=size of training dataset. Not used
 *           right now, needed to evaluate error estimate. \n
 *      'noise_level': h5 dataset, REAL8. Noise parameter, not needed for
 *           evaluating fit. \n
 *      'y_train_mean': h5 dataset, REAL8. Mean value before doing GPR fit,
 *           usually zero. We already remove the mean and store it in
 *           data_mean. \n
 */
int NRHybSur_Init(
    NRHybSurData *NR_hybsur_data, /**< Output: Struct to save surrogate data. */
    LALH5File *file               /**< Opened HDF5 file. */
) {

    if (NR_hybsur_data == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "NR_hybsur_data should not be NULL");
    }
    if (file == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "file should not be NULL");
    }

    if (NR_hybsur_data->setup) {
        XLAL_ERROR(XLAL_FAILURE,
            "Model was already initialized. Ignoring.");
    }

    gsl_vector *domain = NULL;
    int ret = ReadHDF5RealVectorDataset(file, "domain", &domain);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load domain.");
    }
    NR_hybsur_data->domain = domain;

    gsl_matrix *x_train = NULL;
    ret = ReadHDF5RealMatrixDataset(file, "GPR_X_train", &x_train);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load GPR_X_train.");
    }
    NR_hybsur_data->x_train = x_train;

    // dimension of the surrogate
    NR_hybsur_data->params_dim = x_train->size2;

    gsl_matrix_long *mode_list = NULL;
    ret = ReadHDF5LongMatrixDataset(file, "mode_list", &mode_list);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load mode_list.");
    }
    NR_hybsur_data->mode_list = mode_list;

    UINT4 num_modes_modeled = mode_list->size1;
    NR_hybsur_data->num_modes_modeled = num_modes_modeled;

    ModeDataPieces **mode_data_pieces
        = XLALMalloc((NR_hybsur_data->num_modes_modeled)
            * sizeof(*mode_data_pieces));

    // These are needed for the TaylorT3 term
    int phaseAlignIdx;
    ret = ReadHDF5IntDataset(&phaseAlignIdx, file, "phaseAlignIdx");
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load phaseAlignIdx.");
    }

    REAL8 TaylorT3_t_ref;
    ret = ReadHDF5DoubleDataset(&TaylorT3_t_ref, file, "TaylorT3_t_ref");
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load TaylorT3_t_ref.");
    }

    // This stores a precomputed vector needed in computing the TaylorT3 term
    gsl_vector *TaylorT3_factor_without_eta = gsl_vector_alloc(domain->size);
    if (TaylorT3_factor_without_eta == NULL){
        XLAL_ERROR(XLAL_ENOMEM, "gsl_vector_alloc failed.");
    }

    for (UINT4 i=0; i<domain->size; i++){
        const REAL8 tVal = gsl_vector_get(domain, i);
        const REAL8 theta_without_eta = pow((TaylorT3_t_ref - tVal)/5., -1./8);
        gsl_vector_set(TaylorT3_factor_without_eta, i,
            -2./pow(theta_without_eta, 5));
    }
    NR_hybsur_data->TaylorT3_factor_without_eta = TaylorT3_factor_without_eta;
    NR_hybsur_data->phaseAlignIdx = phaseAlignIdx;

    for (UINT4 mode_idx = 0; mode_idx < num_modes_modeled; mode_idx++) {
        ret = NRHybSur_LoadSingleModeData(mode_data_pieces, file,
                mode_idx, mode_list);
        if (ret != XLAL_SUCCESS) {
            XLAL_ERROR(XLAL_EFUNC, "Failed to load mode data.");
        }
    }

    NR_hybsur_data->mode_data_pieces = mode_data_pieces;

    if (ret == XLAL_SUCCESS){
        NR_hybsur_data->setup = 1;
    }

    return ret;
}






//******************* Surrogate evaluation functions **********************/

/**
 * The GPR Kernel evaluated at a single pair of points.
 *
 * We follow sklearn closely. We use the kernel given in
 * Eq.(S3) of arxiv:1809.09125 :
 * \f[
 * K(x, x') = \sigma_k^2 exp(-1/2 \sum_i^D (x^{i} - x'^{i})^2/\sigma_i^2)
 * \f]
 * where D is the dimension of the model.
 *
 * Note that Eq.(S3) also includes the WhiteKernel, which has a noise
 * parameter, but we don't need that here since we only need to evaluate \f$
 * K_{x* x} \f$ when evaluating the fits (See the mean value in Eq.(S2) of same
 * paper). The other term we need is alpha = \f$ K_{x x}^{-1} {\bf f}\f$, which
 * involves the WhiteKernel, but is precomputed offline. alpha is a vector of
 * size N, where N is the number of cases in the training data set.
 */
static REAL8 kernel(
    const gsl_vector *x1,                 /**< Parameter space point 1. */
    const gsl_vector *x2,                 /**< Parameter space point 2. */
    const GPRHyperParams *hyperparams,    /**< GPR hyperparameters. */
    gsl_vector *dummy_worker        /**< Dummy worker array for computations. */
    )
{
    const gsl_vector ls = *hyperparams->length_scale;

    XLAL_CHECK(
        (x1->size == x2->size) && (x1->size == ls.size)
        && (x1->size == dummy_worker->size), XLAL_EDIMS,
        "Mismatch in size of x1, x2, dummy_worker, ls: %zu, %zu, %zu, %zu.\n",
        x1->size, x2->size, dummy_worker->size, ls.size);

    gsl_vector_memcpy(dummy_worker, x1);
    gsl_vector_sub(dummy_worker, x2);
    gsl_vector_div(dummy_worker, &ls); // (x1 - x2) / ls
    const REAL8 r = gsl_blas_dnrm2(dummy_worker);

    // RBF kernel
    return hyperparams->constant_value * exp(-r*r/2.0);
}


/**
 * Evaluate a GPR fit. See Eq.(S2) of arxiv:1809.09125.
 */
static REAL8 gp_predict(
    const gsl_vector *xst,      /**< Point \f$ x_* \f$ where fit will be
                                evaluated. */
    const GPRHyperParams *hyperparams,  /**< Hyperparams for the GPR kernel. */
    const gsl_matrix *x_train,          /**< Training set points. */
    gsl_vector *dummy_worker    /**< Dummy worker array for computations. */
    )
{
    // Evaluate vector K_*
    const UINT4 n = x_train->size1;
    gsl_vector *Kst = gsl_vector_alloc(n);
    if (Kst == NULL){
        XLAL_ERROR(XLAL_ENOMEM, "gsl_vector_alloc failed.");
    }
    for (UINT4 i=0; i < n; i++) {
        const gsl_vector x = gsl_matrix_const_row(x_train, i).vector;
        const REAL8 ker = kernel(xst, &x, hyperparams, dummy_worker);
        gsl_vector_set(Kst, i, ker);
    }

    // Evaluate y_*
    REAL8 res = 0;
    gsl_blas_ddot(Kst, hyperparams->alpha, &res);
    gsl_vector_free(Kst);

    return res + hyperparams->y_train_mean;
}



/**
 * Evaluate a NRHybSur fit.
 */
REAL8 NRHybSur_eval_fit(
    const NRHybSurFitData *fit_data,   /**< Data for fit. */
    const gsl_vector *fit_params, /**< Parameter space point to evaluate the fit
                                at. size=D, the dimension of the model. */
    const gsl_matrix *x_train,    /**< Training set points. */
    gsl_vector *dummy_worker      /**< Dummy worker array for computations. */
) {

    // Get GPR evaluation
    REAL8 fit_val = gp_predict(fit_params, fit_data->hyperparams,
            x_train, dummy_worker);

    // The data was zero-meaned and normalized before constructing the fit,
    // now do the inverse of that.
    fit_val = fit_val * fit_data->data_std + fit_data->data_mean;

    // A linear fit was removed first, now add that back
    gsl_vector_memcpy(dummy_worker, fit_data->lin_coef);
    gsl_vector_mul(dummy_worker, fit_params);
    for (UINT4 i=0; i < dummy_worker->size; i++) {
        fit_val += gsl_vector_get(dummy_worker, i);
    }
    fit_val += fit_data->lin_intercept;

    return fit_val;
}

/**
 * Evaluate a single NRHybSur waveform data piece.
 */
static int NRHybSur_eval_data_piece(
    gsl_vector **result, /**< Output: Should have already been assigned
                         space.*/
    const gsl_vector *fit_params, /**< Parameter space point to evaluate the
                                 fit at. size=D, the dimension of the model. */
    const DataPiece *data_piece,  /**< The waveform data piece to evaluate */
    const gsl_matrix *x_train,        /**< Training set points. */
    gsl_vector *dummy_worker    /**< Dummy worker array for computations. */
) {

    gsl_vector *nodes = gsl_vector_alloc(data_piece->n_nodes);
    if (nodes == NULL){
        XLAL_ERROR(XLAL_ENOMEM, "gsl_vector_alloc failed.");
    }
    for (int i=0; i < data_piece->n_nodes; i++) {
        const REAL8 fit_val = NRHybSur_eval_fit(data_piece->fit_data[i],
            fit_params, x_train, dummy_worker);
        gsl_vector_set(nodes, i, fit_val);
    }

    // Evaluate the empirical interpolant
    gsl_blas_dgemv(CblasTrans, 1.0, data_piece->ei_basis, nodes, 0.0, *result);
    gsl_vector_free(nodes);

    return XLAL_SUCCESS;
}


/**
 * Do cubic spline interpolation using a gsl_interp_cspline.
 * Results differ slightly from scipy.interpolate.InterpolatedUnivariateSpline
 * due to different boundary conditions.
 *
 * gwsurrogate by default will use gsl_interp_cspline, and results should agree
 * exactly. However, gwsurrogate also has the option to use scipy interpolation,
 * in which case the waveform values at the start and end can be slightly
 * different. gwsurrogate will only use the scipy interpolation if for some
 * reasong the gsl build fails.
 */
static gsl_vector *spline_array_interp(
    const gsl_vector *xout,   /**< The vector of points onto which we want to
                          interpolate. */
    const gsl_vector *x,      /**< The x values of the data to interpolate. */
    const gsl_vector *y       /**< The y values of the data to interpolate. */
) {
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, x->size);
    gsl_spline_init(spline, x->data, y->data, x->size);

    gsl_vector *res = gsl_vector_alloc(xout->size);
    REAL8 tmp;
    for (UINT4 i=0; i<xout->size; i++) {
        tmp = gsl_spline_eval(spline, gsl_vector_get(xout, i), acc);
        gsl_vector_set(res, i, tmp);
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return res;
}

/**
 * Computes omega_22 by forward differences.
 */
static REAL8 compute_omega_22(
    const gsl_vector *times,   /**< Time array. */
    const gsl_vector *phi_22,  /**< Phase of (2,2) mode. */
    const int idx /**< Index to compute omega_22 at. */
) {
    const REAL8 t_next =  gsl_vector_get(times, idx+1);
    const REAL8 t =  gsl_vector_get(times, idx);
    const REAL8 phase_next =  gsl_vector_get(phi_22, idx+1);
    const REAL8 phase =  gsl_vector_get(phi_22, idx);

    // Compute omegaM_22 by forward differences
    const REAL8 omegaM_22 = (phase_next - phase)/(t_next - t);

    return omegaM_22;
}


/**
 * Finds closest index such that omegaM_22[index] = omegaM_22_val.
 */
static int search_omega_22(
    int *omega_idx,      /**< Output: closest index such that
                     omegaM_22[index] = omegaM_22_val. */
    const gsl_vector *times,   /**< Time array. */
    const gsl_vector *phi_22,  /**< Phase of (2,2) mode. */
    const REAL8 omegaM_22_val /**< Desired frequency of (2,2) mode. */
) {
    REAL8 omegaM_22 = -10;
    REAL8 omegaM_22_prev = -10;
    int idx = -1;
    while (omegaM_22 < omegaM_22_val){
        // save omegaM_22 of previous index and increase index
        omegaM_22_prev = omegaM_22;
        idx++;

        // Compute omegaM_22 by forward differences
        omegaM_22 = compute_omega_22(times, phi_22, idx);

        if (gsl_vector_get(times, idx) > 0){
            XLAL_ERROR(XLAL_EINVAL,
                "fMin or fRef is larger than the peak frequency=%.7f,"
                " reduce it. Note that this is in code units, radians/M.",
                omegaM_22);
        }
    }

    // At this point idx corresponds to first index with
    // omegaM_22 > omegaM_22_val. If idx-1 is closer to omegaM_22_val,
    // we pick that instead.
    if (fabs(omegaM_22_prev - omegaM_22_val) < fabs(omegaM_22-omegaM_22_val)){
        idx -= 1;
    }

    // output idx
    *omega_idx = idx;

    return XLAL_SUCCESS;
}

/**
 * Evaluate the phase of the (2,2) mode on the sprase surrogate domain.
 *
 * The surrogate actually models the phase residual \f$ \phi^{res}_{22} \f$
 * defined in Eq.(44) of arxiv:1812.07865. Here we first evaluate that and then
 * add the 0PN TaylorT3 phase to get the (2,2) mode phase.
 */
static int NRHybSur_eval_phase_22_sparse(
    gsl_vector **phi_22_sparse,/**< Output: (2,2) mode phase on sparse
                                    domain.*/
    const REAL8 eta,             /**< Symmetric mass ratio. */
    const gsl_vector *fit_params, /**< Parameter space point to evaluate the
                                fit at. size=D, the dimension of the model. */
    gsl_vector *dummy_dp,       /**< Dummy vector to store phase evaluation. */
    const gsl_matrix *x_train,  /**< Training set points. */
    gsl_vector *dummy_worker,   /**< Dummy worker array for computations. */
    const NRHybSurData *NR_hybsur_data /**< Loaded surrogate data. */
){

    // The first mode in mode_list is always the (2,2) mode
    const ModeDataPieces *data_pieces = NR_hybsur_data->mode_data_pieces[0];

    // sanity check
    if (data_pieces->ell != 2 || data_pieces->m != 2){
        XLAL_ERROR(XLAL_EINVAL, "Expected first mode to be (2,2)");
    }

    // Get phi_22_residual, this is phi_22 with the 0PN phase subtracted.
    // The 0PN contribution will be added below.
    int ret = NRHybSur_eval_data_piece(&dummy_dp, fit_params,
            data_pieces->phase_res_data_piece, x_train, dummy_worker);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to evaluate phase_res_data_piece.");
    }

    gsl_vector_memcpy(*phi_22_sparse, dummy_dp);


    // compute 0PN TaylorT3 phase
    gsl_vector *phi_22_T3 = gsl_vector_alloc((*phi_22_sparse)->size);
    gsl_vector_memcpy(phi_22_T3, NR_hybsur_data->TaylorT3_factor_without_eta);
    gsl_vector_scale(phi_22_T3, 1./pow(eta, 3./8.));

    // set phi_22_T3 to zero at phaseAlignIdx
    const int phaseAlignIdx = NR_hybsur_data->phaseAlignIdx;
    gsl_vector_add_constant(phi_22_T3,
        -gsl_vector_get(phi_22_T3, phaseAlignIdx));

    // Add 0PN TaylorT3 phase to get the (2,2) mode phase
    gsl_vector_add(*phi_22_sparse, phi_22_T3);

    gsl_vector_free(phi_22_T3);

    return ret;
}

/**
 * Upsamples sparse \f$ \phi_{22} \f$ such that time step size is deltaTOverM,
 * and initial frequency is roughly omegaM_22_min. The initial frequency will
 * be refined later on.
 */
static int NRHybSur_upsample_phi_22(
    gsl_vector **phi_22_dense, /**< Output: Densely sampled (2,2) mode phase.*/
    gsl_vector **times_dense,  /**< Output: Dense time array. */
    const REAL8 deltaTOverM,  /**< Time step in M. */
    const REAL8 omegaM_22_min, /**< Start frequency of (2,2) mode in rad/M. */
    const gsl_vector *phi_22_sparse, /**< (2,2) mode phase on sparse domain. */
    const gsl_vector *domain         /**< Sparse surrogate domain. */
){

    // Sanity checks
    const REAL8 min_allowed_omegaM_22
        = compute_omega_22(domain, phi_22_sparse, 0);
    if (omegaM_22_min < min_allowed_omegaM_22){
        XLAL_ERROR(XLAL_EINVAL,
            "fMin=%.7f is lesser than the minimum allowed value=%.7f."
            " Note that these are in code units, radians/M.",
            omegaM_22_min, min_allowed_omegaM_22);
    }

    // Find init_idx such that omegaM_22 ~ omegaM_22_min, using the
    // sparse data. This will be refined below.
    int init_idx;
    int ret = search_omega_22(&init_idx, domain, phi_22_sparse, omegaM_22_min);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed fMin search.\n");
    }

    // Go a few indices below to ensure omegaM_22_min is included.
    // This is sometimes required because the search based on the sparse
    // phi_22_sparse above can be inaccurate.
    init_idx -= 5;

    // But if going out of bounds, just choose the first index.
    if (init_idx < 0){
        init_idx = 0;
    }

    // Truncated sparse phi_22 and domain
    gsl_vector *domain_trunc
        = gsl_vector_alloc(phi_22_sparse->size - init_idx);
    gsl_vector *phi_22_sparse_trunc = gsl_vector_alloc(domain_trunc->size);
    for (UINT4 j=0; j < domain_trunc->size; j++) {
        gsl_vector_set(domain_trunc, j, gsl_vector_get(domain, j+init_idx));
        gsl_vector_set(phi_22_sparse_trunc, j,
                gsl_vector_get(phi_22_sparse, j+init_idx));
    }

    // initial time (this is temporary, as we will truncate further to refine
    // initial frequency)
    const REAL8 t0 = gsl_vector_get(domain_trunc, 0);

    // final time
    const REAL8 tf = gsl_vector_get(domain_trunc, domain_trunc->size-1);

    // Setup times_dense, this a dense time array with the required
    // step size. Later on, we will exclude frequencies < omegaM_22_min.
    const int num_times = (int) ceil((tf - t0)/deltaTOverM);
    *times_dense = gsl_vector_alloc(num_times);
    for (int j=0; j < num_times; j++) {
        gsl_vector_set(*times_dense, j, t0 + deltaTOverM*j);
    }

    // interpolate phi_22_sparse_trunc onto times_dense
    *phi_22_dense = spline_array_interp(*times_dense, domain_trunc,
        phi_22_sparse_trunc);

    gsl_vector_free(phi_22_sparse_trunc);
    gsl_vector_free(domain_trunc);

    return ret;
}


/**
 * Evaluate the phase of the (2,2) mode
 *
 * The surrogate actually models the phase residual \f$ \phi^{res}_{22} \f$
 * defined in Eq.(44) of arxiv:1812.07865. Here we first evaluate that and then
 * add the 0PN TaylorT3 phase to get the (2,2) mode phase.
 *
 * Sets the orbital phase to phiRef at the reference frequency omegaM_22_ref.
 * The orbital phase is obtained as phi_22/2, so this leaves a pi ambiguity.
 * But the surrogate data is already aligned such that the heavier BH is on
 * the +ve x-axis at t=-1000M. See Sec.VI.A.4 of arxiv:1812.07865, the resolves
 * the pi ambiguity. This means that the after the realignment, the orbital
 * phase at reference frequency omegaM_22_ref is phiRef, or the heavier BH is
 * at azimuthal angle = phiRef from the +ve x-axis.
 *
 * Only uses data at (2,2) mode frequencies >= omegaM_22_min. This determines
 * the start time. The start time, along with the step size deltaTOverM, is used
 * to determine the output_times. Uses cubic spline interpolation to
 * interpolate from the surrogate's time array to output_times.
 */
int NRHybSur_eval_phase_22(
    gsl_vector **phi_22,        /**< Output: (2,2) mode phase. */
    gsl_vector **output_times,  /**< Output: Time array. */
    const REAL8 eta,                 /**< Symmetric mass ratio. */
    const gsl_vector *fit_params, /**< Parameter space point to evaluate the fit
                                at. size=D, the dimension of the model. */
    const REAL8 omegaM_22_min, /**< Start frequency of (2,2) mode in rad/M. */
    const REAL8 deltaTOverM,   /**< Time step in M. */
    const REAL8 phiRef,        /**< Orbital phase at reference frequency. */
    const REAL8 omegaM_22_ref, /**< Reference freq of (2,2) mode in rad/M. */
    gsl_vector *dummy_dp,       /**< Dummy vector to store phase evaluation. */
    const gsl_matrix *x_train,  /**< Training set points. */
    gsl_vector *dummy_worker,   /**< Dummy worker array for computations. */
    const NRHybSurData *NR_hybsur_data  /**< Loaded surrogate data. */
) {

    if (omegaM_22_ref + 1e-13 < omegaM_22_min){
        XLAL_ERROR(XLAL_EINVAL, "fRef cannot be lesser than fMin.");
    }

    // Evaluate phi_22 on sparse surrogate domain
    const gsl_vector *domain = NR_hybsur_data->domain;
    gsl_vector *phi_22_sparse = gsl_vector_alloc(domain->size);
    int ret = NRHybSur_eval_phase_22_sparse(&phi_22_sparse, eta, fit_params,
            dummy_dp, x_train, dummy_worker, NR_hybsur_data);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed phi_22 sparse evaluation.\n");
    }

    // Upsample to time step deltaTOverM, but restrict to frequencies
    // approximately >= omegaM_22_min.
    gsl_vector *phi_22_dense = NULL;
    gsl_vector *times_dense = NULL;
    ret = NRHybSur_upsample_phi_22(&phi_22_dense, &times_dense, deltaTOverM,
        omegaM_22_min, phi_22_sparse, domain);
    gsl_vector_free(phi_22_sparse);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed phi_22 upsampling.\n");
    }

    // Now refine the start frequency. Find start_idx such that
    // omegaM_22 = omegaM_22_min, and drop everyting below
    int start_idx;
    ret = search_omega_22(&start_idx, times_dense, phi_22_dense,
            omegaM_22_min);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed fMin search.\n");
    }

    *output_times = gsl_vector_alloc(times_dense->size - start_idx);
    *phi_22 = gsl_vector_alloc((*output_times)->size);
    for (UINT4 i=0; i < (*output_times)->size; i++){
        gsl_vector_set(*phi_22, i, gsl_vector_get(phi_22_dense, i+start_idx));
        gsl_vector_set(*output_times, i,
                gsl_vector_get(times_dense, i+start_idx));
    }
    gsl_vector_free(phi_22_dense);
    gsl_vector_free(times_dense);

    // Set reference orbital phase at omegaM_22_ref

    // If omegaM_22_ref is the same as omegaM_22_min, ref_idx should be 0
    int ref_idx = 0;

    // Else, find ref_idx such that omegaM_22[ref_idx] = omegaM_22_ref
    if (fabs(omegaM_22_ref - omegaM_22_min)/omegaM_22_min > 1e-13){
        ret = search_omega_22(&ref_idx, *output_times, *phi_22, omegaM_22_ref);
        if (ret != XLAL_SUCCESS) {
            XLAL_ERROR(XLAL_EFUNC, "Failed fRef search.\n");
        }
    }

    // set (2,2) phase to 2*phiRef at ref_idx
    gsl_vector_add_constant(*phi_22,
        -gsl_vector_get(*phi_22,ref_idx)+2*phiRef);

    return ret;
}


/**
 * Evaluate waveform data pieces of a single mode.
 *
 * For (2,2) mode we model the amplitude and phase, but the phase is
 * evaluated using NRHybSur_eval_phase_22, since it is required for all
 * modes to transform from the coorbital frame to the inertial frame.
 * For all other modes we evaluate the real and imaginary parts of the coorbital
 * frame strain, defined in Eq.(39) of arxiv:1812.07865.
 *
 * Only the required data pieces for a given mode will be evaluated.
 */
int NRHybSur_eval_mode_data_pieces(
    EvaluatedDataPieces **this_mode_eval_dp, /**< Output: evaluated waveform
                                            data pieces of a single mode. */
    const UINT4 ell,                 /**< \f$\ell\f$ index of mode. */
    const UINT4 m,                   /**< m index of mode. */
    const ModeDataPieces *data_pieces,/**< Surrogate data pieces of this mode.*/
    const gsl_vector *output_times,   /**< Time vector to evaluate at. */
    const gsl_vector *fit_params, /**< Parameter space point to evaluate the fit
                                at. size=D, the dimension of the model. */
    gsl_vector *dummy_dp, /**< Dummy vector to store phase evaluation. */
    const gsl_matrix *x_train,        /**< Training set points. */
    gsl_vector *dummy_worker,    /**< Dummy worker array for computations. */
    const NRHybSurData *NR_hybsur_data  /**< Loaded surrogate data. */
) {

    int ret;
    const gsl_vector *domain = NR_hybsur_data->domain;
    (*this_mode_eval_dp)->ell = ell;
    (*this_mode_eval_dp)->m = m;

    if (ell == 2 && m ==2){

        // The phase was already evaluated so, only evaluate the
        // amplitude
        ret = NRHybSur_eval_data_piece(&dummy_dp, fit_params,
                data_pieces->ampl_data_piece, x_train, dummy_worker);
        if (ret != XLAL_SUCCESS) {
            XLAL_ERROR(XLAL_EFUNC, "Failed (2,2) mode amplitude evaluation.\n");
        }
        (*this_mode_eval_dp)->ampl_eval
            = spline_array_interp(output_times, domain, dummy_dp);

    } else {
        // For m=0, l=even, the imaginary part is zero, so we only need to
        // evaluate the real part. But when m!=0, we still want to evaluate
        // the real part.
        if (m != 0 || ell % 2 == 0) {

            // evaluate real part of coorbital frame mode
            ret = NRHybSur_eval_data_piece(&dummy_dp, fit_params,
                    data_pieces->coorb_re_data_piece, x_train, dummy_worker);
            if (ret != XLAL_SUCCESS) {
                XLAL_ERROR(XLAL_EFUNC, "Failed (%u,%u) mode real part evaluation.\n",
                    ell, m);
            }
            // interpolate on to output_times
            (*this_mode_eval_dp)->coorb_re_eval
                = spline_array_interp(output_times, domain, dummy_dp);
        }

        // For m=0, l=odd, the imaginary part is zero, so we only need to
        // evaluate the imaginary part. But when m!=0, we still want to
        // evaluate the imaginary part.
        if (m != 0 || ell % 2 == 1) {

            // evaluate imaginary part of coorbital frame mode
            ret = NRHybSur_eval_data_piece(&dummy_dp, fit_params,
                    data_pieces->coorb_im_data_piece, x_train, dummy_worker);
            if (ret != XLAL_SUCCESS) {
                XLAL_ERROR(XLAL_EFUNC, "Failed (%u,%u) mode imag part evaluation.\n",
                    ell, m);
            }

            // interpolate on to output_times
            (*this_mode_eval_dp)->coorb_im_eval
                = spline_array_interp(output_times, domain, dummy_dp);
        }
    }

    return ret;
}


/**
 * Destroy phi_22 and an EvaluatedDataPieces structure.
 *
 * Free all associated memory.
 */
void NRHybSur_DestroyEvaluatedDataPieces(
    gsl_vector *phi_22,     /**< \f$\phi_{22}\f$ data piece. */
    EvaluatedDataPieces **evaluated_mode_dps, /**< All other data pieces. */
    const UINT4 num_modes_incl  /**< Number of models included. */
) {

    // The phase of the (2,2) mode is always evaluated, free it.
    gsl_vector_free(phi_22);

    // Free all other data pieces
    for (UINT4 mode_idx = 0; mode_idx < num_modes_incl; mode_idx++){
        EvaluatedDataPieces *this_mode_eval_dp = evaluated_mode_dps[mode_idx];
        const UINT4 ell = this_mode_eval_dp->ell;
        const UINT4 m = this_mode_eval_dp->m;

        if (ell == 2 && m ==2){
            // For (2,2) mode we only have the amplitude, the phase is
            // evaluated separately
            gsl_vector_free(this_mode_eval_dp->ampl_eval);
        } else {
            // for m=0, l=odd, the real part is zero and not evaluated. For all
            // other cases free the real part
            if (m != 0 || ell % 2 == 0) {
                gsl_vector_free(this_mode_eval_dp->coorb_re_eval);
            }

            // for m=0, l=even, the imaginary part is zero and not evaluated.
            // For all other cases free the imaginary part
            if (m != 0 || ell % 2 == 1) {
                gsl_vector_free(this_mode_eval_dp->coorb_im_eval);
            }
        }
        XLALFree(this_mode_eval_dp);
    }
    XLALFree(evaluated_mode_dps);

}


/**
 * Sanity check (warning only, not error) that the sample rate is high enough
 * to capture the ringdown frequencies, by ensuring Nyquist frequency is
 * greater than the QNM frequency of the (max_ell,max_ell,0) mode, where
 * max_ell is the maximum ell index among the included modes.
 */
int NRHybSur_sanity_check_sample_rate(
    REAL8 deltaT,                   /**< Sampling interval (s). */
    REAL8 m1,                       /**< Mass of Bh1 (kg). */
    REAL8 m2,                       /**< Mass of Bh2 (kg). */
    REAL8 chi1z,                    /**< Dimensionless spin of Bh1. */
    REAL8 chi2z,                    /**< Dimensionless spin of Bh2. */
    UINT4 max_ell                   /**< Max ell index included. */
){
    /* Ringdown freq used to check the sample rate */
    COMPLEX16Vector modefreqVec;
    COMPLEX16 modeFreq;

    /* Before calculating everything else, check sample freq is high enough */
    modefreqVec.length = 1;
    modefreqVec.data = &modeFreq;

    // m=ell mode should have the highest frequency
    UINT4 mode_highest_freqL = max_ell;
    UINT4 mode_highest_freqM = max_ell;
    UINT4 num_overtones = 1; // One overtone should be good enough for a test
    Approximant SpinAlignedEOBapproximant = SEOBNRv4;

    const REAL8 spin1[3] = {0, 0, chi1z};
    const REAL8 spin2[3] = {0, 0, chi2z};

    // NOTE: XLALSimIMREOBGenerateQNMFreqV2 expects masses in Solar mass
    int ret = XLALSimIMREOBGenerateQNMFreqV2(&modefreqVec, m1/LAL_MSUN_SI,
            m2/LAL_MSUN_SI, spin1, spin2, mode_highest_freqL,
            mode_highest_freqM, num_overtones, SpinAlignedEOBapproximant);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "XLALSimIMREOBGenerateQNMFreqV2 failed");
    }

    const REAL8 nyquist_freq = 1./2./deltaT;
    const REAL8 ringdown_freq = creal(modeFreq)/(2.*LAL_PI);

    if (nyquist_freq < ringdown_freq){
        XLAL_PRINT_WARNING("Nyquist frequency=%.7f Hz is lesser than the QNM"
            " frequency=%.7f Hz of the (%u,%u,0) mode. Consider reducing time"
            " step.", nyquist_freq, ringdown_freq, max_ell, max_ell);
    }

    return XLAL_SUCCESS;
}


/**
 * Activates all modes of an NRHybSur model.  For NRHybSur3dq8 that is \f$ \ell
 * \leq 4, m \geq 0 \f$, and (5,5), but not (4,1) or (4,0).
 */
void NRHybSur_set_default_modes(
    LALValue *ModeArray,                /**< Output: Container for modes. */
    const NRHybSurData *NR_hybsur_data  /**< Loaded surrogate data. */
    )
{
    // list of available modes
    const gsl_matrix_long *mode_list = NR_hybsur_data->mode_list;
    // number of available modes
    const UINT4 num_modes_modeled = NR_hybsur_data->num_modes_modeled;
    UINT4 ell, m;
    for (UINT4 idx = 0; idx < num_modes_modeled; idx++){
        ell = gsl_matrix_long_get(mode_list, idx, 0);
        m = gsl_matrix_long_get(mode_list, idx, 1);
        XLALSimInspiralModeArrayActivateMode(ModeArray, ell, m);
    }
}

/**
 * Sanity checks on ModeArray.
 *
 * Will raise an error if an unavailable mode is requested. Note that we only
 * look for m>=0 modes in ModeArray, and will ignore m<0 modes even if present.
 * The m<0 modes automatically get added when evaluting the waveform.
 */
int NRHybSur_check_mode_array(
    UINT4 *num_modes_incl,   /**< Output: Number of modes to include. */
    UINT4 *max_ell,          /**< Output: Max ell index included. */
    LALValue *ModeArray,                /**< Container for modes. */
    const NRHybSurData *NR_hybsur_data  /**< Loaded surrogate data. */
    )
{
    INT4 modeAvailable = 0;
    UINT4 ell, m;

    // number of available modes
    const UINT4 num_modes_modeled = NR_hybsur_data->num_modes_modeled;

    const gsl_matrix_long *mode_list = NR_hybsur_data->mode_list;

    *num_modes_incl = 0;
    *max_ell = 2;
    for (UINT4 ELL = 2; ELL <= LAL_SIM_L_MAX_MODE_ARRAY; ELL++) {
        for (UINT4 EMM = 0; EMM <= ELL; EMM++) {

            modeAvailable = 0;
            if (XLALSimInspiralModeArrayIsModeActive(ModeArray,ELL,EMM) == 1) {

                for (UINT4 idx = 0; idx < num_modes_modeled; idx++){

                    ell = gsl_matrix_long_get(mode_list, idx, 0);
                    m = gsl_matrix_long_get(mode_list, idx, 1);

                    // Check if the (ELL, EMM) mode is included
                    if ((ell == ELL) && (m == EMM)) {
                        modeAvailable=1;
                        *num_modes_incl += 1;
                        if (ell > *max_ell) {
                            *max_ell = ell;
                        }
                    }
                }

                if (modeAvailable != 1) {
                    XLAL_ERROR(XLAL_EINVAL,
                        "Mode (%d,%d) is not available.",ELL,EMM);
                }
            }
        }
    }

    if (*num_modes_incl == 0) {
        XLAL_ERROR(XLAL_EINVAL, "Zero available modes requested.");
    }

    return XLAL_SUCCESS;
}
