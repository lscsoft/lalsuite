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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALSimIMR.h>


//*************************************************************************/
//************************* struct definitions ****************************/
//*************************************************************************/

/**
 * Data used in a single GPR fit.
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
typedef struct tagGPRHyperParams {
    REAL8 constant_value;    /**< \f$ \sigma_k^2 \f$ in kernel. */
    REAL8 y_train_mean;      /**< Mean value before GPR fit, usually zero. */
    gsl_vector *length_scale; /**< \f$ \sigma_i \f$ in kernel. */
    gsl_vector *alpha;        /**< Precomputed \f$ K_{x x}^{-1} {\bf f}\f$. */
} GPRHyperParams;

/**
 * Data used in a single NRHybSur fit.
 *
 * Described in supplemental materials of arxiv:1809.09125.
 * We first subtract a linear fit to the data. Then we subtract the mean of
 * the resulting values. Then we normalize by dividing by the standard
 * deviation of the resulting values. Finally, we construct a GPR fit for the
 * resulting values. Below, D = dimension of the model.
 */
typedef struct tagNRHybSurFitData {
    REAL8 data_mean;            /**< Mean of data after linear fit. */
    REAL8 data_std;             /**< Standard deviation after removing mean. */
    REAL8 lin_intercept;        /**< Constant intercept from linear fit. */
    gsl_vector *lin_coef;        /**< Linear coefs from linear fit, size=D. */
    GPRHyperParams *hyperparams; /**< Hyperparameters from GPR fit. */
} NRHybSurFitData;

/**
 * NRHybSur data for a single waveform data piece.
 */
typedef struct tagDataPiece {
    int n_nodes;                /**< Number of empirical nodes. */
    NRHybSurFitData **fit_data; /**< NRHybSurFitData at each empirical node. */
    gsl_matrix *ei_basis;       /**< Empirical interpolation matrix. */
} DataPiece;

/**
 * NRHybSur data pieces of a single mode.
 *
 * The different data pieces are described in Sec.VI.B of arxiv:1812.07865.
 * For (2,2) mode we model the amplitude \f$ A_{22} \f$ and phase residual
 * \f$ \phi^{res}_{22} \f$ defined in Eq.(44) of arxiv:1812.07865.
 * For all other modes we model the real and imaginary parts of the coorbital
 * frame strain, defined in Eq.(39) of arxiv:1812.07865.
 * Only the required data pieces for a given mode will be loaded.
 */
typedef struct tagModeDataPieces {
    UINT4 ell;                         /**< Mode \f$ \ell \f$ value. */
    UINT4 m;                           /**< Mode m value. */
    DataPiece *ampl_data_piece;         /**< Amplitude data piece. */
    DataPiece *phase_res_data_piece;    /**< Phase residual data piece. */
    DataPiece *coorb_re_data_piece;/**< Coorbital frame real part data piece. */
    DataPiece *coorb_im_data_piece;/**< Coorbital frame imag part data piece. */
} ModeDataPieces;

/**
 * NRHybSur surrogate data for all modes, to be loaded from a h5 file.
 */
typedef struct tagNRHybSurData {
    UINT4 setup;         /**< Indicates if NRHybSur has been initialized */
    UINT4 num_modes_modeled;   /**< Number of modeled modes. */
    int phaseAlignIdx;          /**< Alignment index for phase data piece. */
    REAL8 params_dim;          /**< Dimensions of the model. */
    gsl_vector *domain;         /**< Common time samples for the surrogate. */
    gsl_vector *TaylorT3_factor_without_eta; /**< Precomputed quantity for use
                        in evaluating the 0PN TaylorT3 phase contribution. */
    gsl_matrix_long *mode_list; /**< List of modeled modes, first element is
                                always (2,2). */
    gsl_matrix *x_train; /**< Training set parameters, needed for GPR fits. */
    ModeDataPieces **mode_data_pieces; /**< Data pieces of all modes, same
                                order as mode_list. */
} NRHybSurData;

/**
 * NRHybSur evaluated data for a single mode
 *
 * For the (2,2) mode only the amplitude is stored in this struct. The phase
 * is computed separately as it is needed for all modes to transform from the
 * coorbital frame to the inertial frame.
 * For all other modes the real and imaginary parts of the strain in the
 * coorbital frame are stored.
 * Only the required data pieces for a given mode will be evaluated.
 */
typedef struct tagEvaluatedDataPieces {
    UINT4 ell;                 /**< Mode \f$ \ell \f$ value. */
    UINT4 m;                   /**< Mode m value. */
    gsl_vector *ampl_eval;      /**< Amplitude data piece evaluation. */
    gsl_vector *coorb_re_eval;  /**< Coorbital frame real part evaluation. */
    gsl_vector *coorb_im_eval;  /**< Coorbital frame imag part evaluation. */
} EvaluatedDataPieces;


//*************************************************************************/
//************************* function declarations *************************/
//*************************************************************************/

int ReadHDF5DoubleDataset(
    REAL8 *param,
    LALH5File *sub,
    const char *name
);

int ReadHDF5IntDataset(
    int *param,
    LALH5File *sub,
    const char *name
);

int NRHybSur_Init(
    NRHybSurData *data,
    LALH5File *file
    );

REAL8 NRHybSur_eval_fit(
    const NRHybSurFitData *fit_data,
    const gsl_vector *fit_params,
    const gsl_matrix *x_train,
    gsl_vector *dummy_worker
);

int NRHybSur_eval_phase_22(
    gsl_vector **phi_22,
    gsl_vector **output_times,
    const REAL8 eta,
    const gsl_vector *fit_params,
    const REAL8 omegaM_22_min,
    const REAL8 deltaTOverM,
    const REAL8 phiRef,
    const REAL8 omegaM_22_ref,
    gsl_vector *dummy_dp,
    const gsl_matrix *x_train,
    gsl_vector *dummy_worker,
    const NRHybSurData *NR_hybsur_data
);

int NRHybSur_eval_mode_data_pieces(
    EvaluatedDataPieces **this_mode_eval_dp,
    const UINT4 ell,
    const UINT4 m,
    const ModeDataPieces *data_pieces,
    const gsl_vector *output_times,
    const gsl_vector *fit_params,
    gsl_vector *dummy_dp,
    const gsl_matrix *x_train,
    gsl_vector *dummy_worker,
    const NRHybSurData *NR_hybsur_data
);

void NRHybSur_DestroyEvaluatedDataPieces(
    gsl_vector *phi_22,
    EvaluatedDataPieces **evaluated_mode_dps,
    const UINT4 num_modes_incl
);

int NRHybSur_sanity_check_sample_rate(
    REAL8 deltaT,
    REAL8 m1,
    REAL8 m2,
    REAL8 chi1z,
    REAL8 chi2z,
    UINT4 max_ell
);

void NRHybSur_set_default_modes(
    LALValue *ModeArray,
    const NRHybSurData *NR_hybsur_data
    );

int NRHybSur_check_mode_array(
    UINT4 *num_modes_incl,
    UINT4 *max_ell,
    LALValue *ModeArray,
    const NRHybSurData *NR_hybsur_data
    );
