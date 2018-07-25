/* *  Copyright (C) 2017 Jonathan Blackman
 *  NRSur7dq2 NR surrogate model.
 *  Paper: https://arxiv.org/abs/1705.07089
 *  Based on the python implementation found at:
 *  https://www.black-holes.org/data/surrogates/index.html
 *  which uses the same NRSur7dq2.hdf5 data file.
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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <alloca.h>
#include <string.h>
#include <libgen.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex_math.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimSphHarmMode.h>
#include <lal/LALSimIMR.h>
#include <lal/H5FileIO.h>

#include "LALSimIMRNRSur7dq2.h"
#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif


#ifdef LAL_PTHREAD_LOCK
static pthread_once_t NRSur7dq2_is_initialized = PTHREAD_ONCE_INIT;
#endif


/**
 * Global surrogate data.
 * This data will be loaded at most once. Any executable which calls NRSur7dq2_Init_LALDATA
 * directly or by calling any XLAL NRSur7dq2 function will have a memory leak according
 * to valgrind, because we never free this memory.
 */
static NRSur7dq2Data __lalsim_NRSur7dq2_data;

/**
 * @addtogroup LALSimNRSur7dq2_c
 * @{
 *
 * @name Routines for NR surrogate model "NRSur7dq2"
 * @{
 *
 * @author Jonathan Blackman
 *
 * @brief C code for NRSur7dq2 NR surrogate waveform model.
 *
 * This is a fully precessing time domain model including all subdominant modes up to ell=4.
 * See Blackman et al \cite Blackman:q27d for details.
 * Any studies that use this waveform model should include a reference to that paper.
 * Using this model requires the file lalsuite-extra/data/lalsimulation/NRSur7dq2.h5
 * Make sure your $LAL_DATA_PATH points to lalsuite-extra/data/lalsimulation/.
 * The lalsuite-extra commit hash at the time of review was 77613e7f5f01d5ea11829ded5677783cafc0d298
 *
 * @note The range of validity of the model is:
 * * Mass ratios 1 <= q <= 2
 * * Spin magnitudes |chi_i| <= 0.8
 * * Total time before merger <= 4500M, which in practice leads to a parameter-dependent lower bound for fmin.
 *
 * @note Additional notes:
 * * Passing in a non-trivial ModeArray controls which co-orbital frame modes are evaluated.
 * * A python version of this model can be installed with "pip install NRSur7dq2".
 * * This lalsimulation implementation has been verified to agree with version 1.0.3 up to
 * very small differences due to slightly differing time interpolation methods.
 * * Note that for conventions to agree with ChooseTDWaveform (and XLALSimInspiralNRSur7dq2Polarizations),
 * you must pass use_lalsimulation_conventions=True when calling the pip version of NRSur7dq2.
 *
 * @review NRSur7dq2 model and routines reviewed by Sebastian Khan, Harald Pfeiffer, Geraint Pratten, and Michael Pürrer.
 * Reviewees were Jonathan Blackman, Scott Field, and Vijay Varma.
 * The review page can be found at https://git.ligo.org/waveforms/reviews/nrsur/wikis/home
 */


/***********************************************************************************/
/****************************** Function Definitions *******************************/
/***********************************************************************************/




/**
 * This needs to be called once, before __lalsim_NRSur7dq2_data is used.
 * It finds the hdf5 data file with the NRSur7dq2 data and calls NRSur7dq2_Init.
 */
static void NRSur7dq2_Init_LALDATA(void) {
    if (NRSur7dq2_IsSetup()) return;

    char *path = XLALFileResolvePathLong(NRSUR7DQ2_DATAFILE, PKG_DATA_DIR);
    if (path==NULL)
        XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", NRSUR7DQ2_DATAFILE);
    char *dir = dirname(path);
    size_t size = strlen(dir) + strlen(NRSUR7DQ2_DATAFILE) + 2;
    char *file_path = XLALMalloc(size);
    snprintf(file_path, size, "%s/%s", dir, NRSUR7DQ2_DATAFILE);

    LALH5File *file = XLALH5FileOpen(file_path, "r");

    int ret = NRSur7dq2_Init(&__lalsim_NRSur7dq2_data, file);

    if (ret != XLAL_SUCCESS)
        XLAL_ERROR_VOID(XLAL_FAILURE, "Failure loading data from %s\n", file_path);

    XLALFree(path);
    XLALFree(file_path);
}

/**
 * Initialize a NRSur7dq2Data structure from an open hdf5 file.
 * This will typically only be called once, from NRSur7dq2_Init_LALDATA.
 */
static int NRSur7dq2_Init(NRSur7dq2Data *data, LALH5File *file) {
    size_t i;

    if (data->setup) {
        XLALPrintError("You tried to setup the NRSur7dq2 model that was already initialized. Ignoring.\n");
        return XLAL_FAILURE;
    }


    // Get the dynamics time nodes
    gsl_vector *t_ds_with_halves = NULL;
    ReadHDF5RealVectorDataset(file, "t_ds", &t_ds_with_halves);
    gsl_vector *t_ds = gsl_vector_alloc(t_ds_with_halves->size - 3);
    gsl_vector *t_ds_half_times = gsl_vector_alloc(3);
    for (i=0; i < 3; i++) {
        gsl_vector_set(t_ds, i, gsl_vector_get(t_ds_with_halves, 2*i));
        gsl_vector_set(t_ds_half_times, i, gsl_vector_get(t_ds_with_halves, 2*i + 1));
    }
    for (i=3; i < t_ds->size; i++) {
        gsl_vector_set(t_ds, i, gsl_vector_get(t_ds_with_halves, i+3));
    }
    gsl_vector_free(t_ds_with_halves);
    data->t_ds = t_ds;
    data->t_ds_half_times = t_ds_half_times;

    // Load the dynamics time node data. The nodes are stored in HDF5 groups "ds_node_%d"%(i)
    // where i=0, 1, 2, ...
    // These INCLUDE the half-nodes at t_1/2, t_3/2, and t_5/2, so the (node time, index) pairs are
    // (t_0, 0), (t_1/2, 1), (t_1, 2), (t_3/2, 3), (t_2, 4), (t_5/2, 5), (t_3, 6), (t_4, 7), ...
    // (t_n, n+3) for n >= 3
    DynamicsNodeFitData **ds_node_data = XLALMalloc( (t_ds->size) * sizeof(*ds_node_data));
    DynamicsNodeFitData **ds_half_node_data = XLALMalloc( 3 * sizeof(*ds_node_data) );
    for (i=0; i < (t_ds->size); i++) ds_node_data[i] = NULL;
    for (i=0; i < 3; i++) ds_half_node_data[i] = NULL;
    LALH5File *sub;
    char *sub_name = XLALMalloc(15); // Should be enough for j < 1000000
    int j;
    for (i=0; i < (t_ds->size); i++) {
        if (i < 3) {j = 2*i;} else {j = i+3;}
        snprintf(sub_name, 15, "ds_node_%d", j);
        sub = XLALH5GroupOpen(file, sub_name);
        NRSur7dq2_LoadDynamicsNode(ds_node_data, sub, i);

        if (i < 3) {
            snprintf(sub_name, 15, "ds_node_%d", j+1);
            sub = XLALH5GroupOpen(file, sub_name);
            NRSur7dq2_LoadDynamicsNode(ds_half_node_data, sub, i);
        }
    }
    XLALFree(sub_name);
    data->ds_node_data = ds_node_data;
    data->ds_half_node_data = ds_half_node_data;

    // Get the coorbital time array
    gsl_vector *t_coorb = NULL;
    ReadHDF5RealVectorDataset(file, "t_coorb", &t_coorb);
    data->t_coorb = t_coorb;

    // Load coorbital waveform surrogate data
    WaveformFixedEllModeData **coorbital_mode_data = XLALMalloc( (NRSUR7DQ2_LMAX - 1) * sizeof(*coorbital_mode_data) );
    for (int ell_idx=0; ell_idx < NRSUR7DQ2_LMAX-1; ell_idx++) {
        NRSur7dq2_LoadCoorbitalEllModes(coorbital_mode_data, file, ell_idx);
    }
    data->coorbital_mode_data = coorbital_mode_data;

    XLAL_PRINT_INFO("Successfully loaded NRSur7dq2 data!");
    data->LMax = NRSUR7DQ2_LMAX;
    data->setup = 1;

    return XLAL_SUCCESS;
}

/**
 * Loads the data for a single dynamics node into a DynamicsNodeFitData struct.
 * This is only called during the initialization of the surrogate data through NRSur7dq2_Init.
 */
static void NRSur7dq2_LoadDynamicsNode(
    DynamicsNodeFitData **ds_node_data, /**< Entry i should be NULL; Will malloc space and load data into it. */
    LALH5File *sub,                     /**< Subgroup containing data for dynamics node i. */
    int i                               /**< Dynamics node index. */
) {
    ds_node_data[i] = XLALMalloc( sizeof(*ds_node_data[i]) );

    // omega
    FitData *omega_data = XLALMalloc(sizeof(FitData));
    omega_data->coefs = NULL;
    omega_data->basisFunctionOrders = NULL;
    ReadHDF5RealVectorDataset(sub, "omega_coefs", &(omega_data->coefs));
    ReadHDF5LongMatrixDataset(sub, "omega_bfOrders", &(omega_data->basisFunctionOrders));
    omega_data->n_coefs = omega_data->coefs->size;
    ds_node_data[i]->omega_data = omega_data;

    // omega_copr
    VectorFitData *omega_copr_data = XLALMalloc(sizeof(VectorFitData));
    omega_copr_data->coefs = NULL;
    omega_copr_data->basisFunctionOrders = NULL;
    omega_copr_data->componentIndices = NULL;
    ReadHDF5RealVectorDataset(sub, "omega_orb_coefs", &(omega_copr_data->coefs));
    ReadHDF5LongMatrixDataset(sub, "omega_orb_bfOrders", &(omega_copr_data->basisFunctionOrders));
    ReadHDF5LongVectorDataset(sub, "omega_orb_bVecIndices", &(omega_copr_data->componentIndices));
    omega_copr_data->n_coefs = omega_copr_data->coefs->size;
    omega_copr_data->vec_dim = 2;
    ds_node_data[i]->omega_copr_data = omega_copr_data;

    // chiA_dot
    VectorFitData *chiA_dot_data = XLALMalloc(sizeof(VectorFitData));
    chiA_dot_data->coefs = NULL;
    chiA_dot_data->basisFunctionOrders = NULL;
    chiA_dot_data->componentIndices = NULL;
    ReadHDF5RealVectorDataset(sub, "chiA_coefs", &(chiA_dot_data->coefs));
    ReadHDF5LongMatrixDataset(sub, "chiA_bfOrders", &(chiA_dot_data->basisFunctionOrders));
    ReadHDF5LongVectorDataset(sub, "chiA_bVecIndices", &(chiA_dot_data->componentIndices));
    chiA_dot_data->n_coefs = chiA_dot_data->coefs->size;
    chiA_dot_data->vec_dim = 3;
    ds_node_data[i]->chiA_dot_data = chiA_dot_data;

    // chiB_dot
    // One chiB_dot node has 0 coefficients, and ReadHDF5RealVectorDataset fails.
    VectorFitData *chiB_dot_data = XLALMalloc(sizeof(VectorFitData));
    chiB_dot_data->coefs = NULL;
    chiB_dot_data->basisFunctionOrders = NULL;
    chiB_dot_data->componentIndices = NULL;

    UINT4Vector *dimLength;
    size_t n;
    LALH5Dataset *dset;
    dset = XLALH5DatasetRead(sub, "chiB_coefs");
    dimLength = XLALH5DatasetQueryDims(dset);
    n = dimLength->data[0];
    if (n==0) {
        chiB_dot_data->n_coefs = 0;
    } else {
        ReadHDF5RealVectorDataset(sub, "chiB_coefs", &(chiB_dot_data->coefs));
        ReadHDF5LongMatrixDataset(sub, "chiB_bfOrders", &(chiB_dot_data->basisFunctionOrders));
        ReadHDF5LongVectorDataset(sub, "chiB_bVecIndices", &(chiB_dot_data->componentIndices));
        chiB_dot_data->n_coefs = chiB_dot_data->coefs->size;
    }
    chiB_dot_data->vec_dim = 3;
    ds_node_data[i]->chiB_dot_data = chiB_dot_data;
}

/**
 * Load the WaveformFixedEllModeData from file for a single value of ell.
 * This is only called during the initialization of the surrogate data through NRSur7dq2_Init.
 */
static void NRSur7dq2_LoadCoorbitalEllModes(
    WaveformFixedEllModeData **coorbital_mode_data, /**< Entry i should be NULL; will malloc space and load data into it.*/
    LALH5File *file, /**< The open NRSur7dq2.hdf5 file */
    int i /**< The index of coorbital_mode_data. Equivalently, ell-2. */
) {
    WaveformFixedEllModeData *mode_data = XLALMalloc( sizeof(*coorbital_mode_data[i]) );
    mode_data->ell = i+2;

    LALH5File *sub;
    int str_size = 30; // Enough for L with 15 digits...
    char *sub_name = XLALMalloc(str_size);

    // Real part of m=0 mode
    snprintf(sub_name, str_size, "hCoorb_%d_0_real", i+2);
    sub = XLALH5GroupOpen(file, sub_name);
    NRSur7dq2_LoadWaveformDataPiece(sub, &(mode_data->m0_real_data), false);

    // Imag part of m=0 mode
    snprintf(sub_name, str_size, "hCoorb_%d_0_imag", i+2);
    sub = XLALH5GroupOpen(file, sub_name);
    NRSur7dq2_LoadWaveformDataPiece(sub, &(mode_data->m0_imag_data), false);

    // NOTE:
    // In the paper https://arxiv.org/abs/1705.07089, Eq. 16 uses
    // h_\pm^{\ell, m} = \frac{1}{2} \left( h_\mathrm{coorb}^{\ell, m} \pm h_\mathrm{coorb}^{\ell, -m\, *} \right)
    // and we often denote h_\pm^{\ell, m} == X_\pm^{\ell, m} in this file and other documentation,
    // but when actually building the surrogate the (\ell, -m) mode was taken as the reference mode, so we store
    // Y_\pm^{\ell, m} = \frac{1}{2} \left( h_\mathrm{coorb}^{\ell, -m} \pm h_\mathrm{coorb}^{\ell, m\, *} \right)
    // Note that we still have everything we want:
    // X_\pm^{\ell, m} = \pm Y_\pm^{\ell, m\, *}
    // or:
    // Re[X_+] = Re[Y_+],
    // Im[X_+] = -Im[Y_+],
    // Re[X_-] = -Re[Y_-],
    // Im[X_-] = Im[Y_-]
    // We work with X rather than Y in this file, so we need the minus signs when loading from Y.
    mode_data->X_real_plus_data = XLALMalloc( (i+2) * sizeof(WaveformDataPiece *) );
    mode_data->X_real_minus_data = XLALMalloc( (i+2) * sizeof(WaveformDataPiece *) );
    mode_data->X_imag_plus_data = XLALMalloc( (i+2) * sizeof(WaveformDataPiece *) );
    mode_data->X_imag_minus_data = XLALMalloc( (i+2) * sizeof(WaveformDataPiece *) );
    for (int m=1; m<=(i+2); m++) {
        snprintf(sub_name, str_size, "hCoorb_%d_%d_Re+", i+2, m);
        sub = XLALH5GroupOpen(file, sub_name);
        NRSur7dq2_LoadWaveformDataPiece(sub, &(mode_data->X_real_plus_data[m-1]), false);
        snprintf(sub_name, str_size, "hCoorb_%d_%d_Re-", i+2, m);
        sub = XLALH5GroupOpen(file, sub_name);
        NRSur7dq2_LoadWaveformDataPiece(sub, &(mode_data->X_real_minus_data[m-1]), true);
        snprintf(sub_name, str_size, "hCoorb_%d_%d_Im+", i+2, m);
        sub = XLALH5GroupOpen(file, sub_name);
        NRSur7dq2_LoadWaveformDataPiece(sub, &(mode_data->X_imag_plus_data[m-1]), true);
        snprintf(sub_name, str_size, "hCoorb_%d_%d_Im-", i+2, m);
        sub = XLALH5GroupOpen(file, sub_name);
        NRSur7dq2_LoadWaveformDataPiece(sub, &(mode_data->X_imag_minus_data[m-1]), false);
    }
    XLALFree(sub_name);
    coorbital_mode_data[i] = mode_data;
}

/**
 * Loads a single NRSur7dq2 coorbital waveform data piece from file into a WaveformDataPiece.
 * This is only called during the initialization of the surrogate data through NRSur7dq2_Init.
 */
static void NRSur7dq2_LoadWaveformDataPiece(
    LALH5File *sub,             /**< HDF5 group containing data for this waveform data piece */
    WaveformDataPiece **data,   /**< Output - *data should be NULL. Space will be allocated. */
    bool invert_sign            /**< If true, multiply the empirical interpolation matrix by -1. */
) {
    *data = XLALMalloc(sizeof(WaveformDataPiece));

    gsl_matrix *EI_basis = NULL;
    ReadHDF5RealMatrixDataset(sub, "EIBasis", &EI_basis);
    if (invert_sign) {
        gsl_matrix_scale(EI_basis, -1);
    }
    (*data)->empirical_interpolant_basis = EI_basis;

    gsl_vector_long *node_indices = NULL;
    ReadHDF5LongVectorDataset(sub, "nodeIndices", &node_indices);
    (*data)->empirical_node_indices = node_indices;

    int n_nodes = (*data)->empirical_node_indices->size;
    (*data)->n_nodes = n_nodes;
    (*data)->fit_data = XLALMalloc( n_nodes * sizeof(FitData *) );

    LALH5File *nodeModelers = XLALH5GroupOpen(sub, "nodeModelers");
    int str_size = 20; // Enough for L with 11 digits...
    char *sub_name = XLALMalloc(str_size);
    for (int i=0; i<n_nodes; i++) {
        FitData *node_data = XLALMalloc(sizeof(FitData));
        node_data->coefs = NULL;
        node_data->basisFunctionOrders = NULL;
        snprintf(sub_name, str_size, "coefs_%d", i);
        ReadHDF5RealVectorDataset(nodeModelers, sub_name, &(node_data->coefs));
        snprintf(sub_name, str_size, "bfOrders_%d", i);
        ReadHDF5LongMatrixDataset(nodeModelers, sub_name, &(node_data->basisFunctionOrders));
        node_data->n_coefs = node_data->coefs->size;
        (*data)->fit_data[i] = node_data;
    }
}

/**
 * Helper function which returns whether or not the global NRSur7dq2 surrogate data has been initialized.
 */
static bool NRSur7dq2_IsSetup(void) {
    if(__lalsim_NRSur7dq2_data.setup)   return true;
    else return false;
}

/**
 * Helper function for integer powers
 */
double ipow(double base, int exponent) {
    if (exponent == 0) return 1.0;
    double res = base;
    while (exponent > 1) {
        res = res*base;
        exponent -= 1;
    }
    return res;
}


/*
 * Evaluate a NRSur7dq2 scalar fit.
 * The fit result is given by
 *      \sum_{i=1}^{n} c_i * \prod_{j=1}^7 B_j(k_{i, j}; x_j)
 * where i runs over fit coefficients, j runs over the 7 dimensional parameter
 * space, and B_j is a basis function, taking an integer order k_{i, j} and
 * the parameter component x_j. For this surrogate, B_j are monomials in the spin
 * components, and monomials in an affine transformation of the mass ratio.
 */
double NRSur7dq2_eval_fit(
    FitData *data,  /**< Data for fit */
    double *x       /**< size 7, giving mass ratio q, and dimensionless spin components */
) {
    double x_powers[22]; // 3 per spin component, 4 for mass ratio
    double res = 0.0;
    double prod;
    int i, j;

    // The fits were constructed using this rather than using q directly
    double q_fit = NRSUR7DQ2_Q_FIT_OFFSET + NRSUR7DQ2_Q_FIT_SLOPE*x[0];

    // Compute powers of components of x
    for (i=0; i<22; i++) {
        if (i%7==0) {
            x_powers[i] = ipow(q_fit, i/7);
        } else {
            x_powers[i] = ipow(x[i%7], i/7);
        }
    }

    // Sum up fit terms
    for (i=0; i < data->n_coefs; i++) {
        // Initialize with q basis function:
        prod = x_powers[7 * gsl_matrix_long_get(data->basisFunctionOrders, i, 0)];
        // Multiply with spin basis functions:
        for (j=1; j<7; j++) {
            prod *= x_powers[7 * gsl_matrix_long_get(data->basisFunctionOrders, i, j) + j];
        }
        res += gsl_vector_get(data->coefs, i) * prod;
    }

    return res;
}

/*
 * This is very similar to NRSur7dq2_eval_fit except that the result is a
 * 2d or 3d vector instead of a scalar. Each fit coefficient now applies to
 * just a single component of the result.
 */
static void NRSur7dq2_eval_vector_fit(
    double *res,            /**< Result */
    VectorFitData *data,    /**< Data for fit */
    double *x               /**< size 7, giving mass ratio q, and dimensionless spin components */
) {
    double x_powers[22]; // 3 per spin component, 4 for mass ratio
    double prod;
    int i, j;

    // Initialize the result
    for (i=0; i < data->vec_dim; i++) {
        res[i] = 0.0;
    }

    // The fits were constructed using this rather than using q directly
    double q_fit = NRSUR7DQ2_Q_FIT_OFFSET + NRSUR7DQ2_Q_FIT_SLOPE*x[0];

    // Compute powers of components of x
    for (i=0; i<22; i++) {
        if (i%7==0) {
            x_powers[i] = ipow(q_fit, i/7);
        } else {
            x_powers[i] = ipow(x[i%7], i/7);
        }
    }

    // Sum up fit terms
    for (i=0; i < data->n_coefs; i++) {
        // Initialize with q basis function:
        prod = x_powers[7 * gsl_matrix_long_get(data->basisFunctionOrders, i, 0)];
        // Multiply with spin basis functions:
        for (j=1; j<7; j++) {
            prod *= x_powers[7 * gsl_matrix_long_get(data->basisFunctionOrders, i, j) + j];
        }
        res[gsl_vector_long_get(data->componentIndices, i)] += gsl_vector_get(data->coefs, i) * prod;
    }
}

/* During the ODE integration, the norm of the spins will change due to
 * integration errors and fit modeling errors. Keep them normalized.
 * Normalizes in-place
 */
static void NRSur7dq2_normalize_y(
    double chiANorm,    /**< ||vec{chi}_A|| */
    double chiBNorm,    /**< ||vec{chi}_B|| */
    double *y           /**< [q0, qx, qy, qz, phi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz] */
) {
    double nQ, nA, nB, sum;
    int i;

    // Compute current norms
    sum = 0.0;
    for (i=0; i<4; i++) {
        sum += y[i]*y[i];
    }
    nQ = sqrt(sum);

    sum = 0.0;
    for (i=5; i<8; i++) {
        sum += y[i]*y[i];
    }
    nA = sqrt(sum);

    sum = 0.0;
    for (i=8; i<11; i++) {
        sum += y[i]*y[i];
    }
    nB = sqrt(sum);

    // Compute normalized output
    for (i=0; i<4; i++) {
        y[i] = y[i] / nQ;
    }
    y[4] = y[4];
    for (i=5; i<8; i++) {
        y[i] = y[i] * chiANorm / nA;
    }
    for (i=8; i<11; i++) {
        y[i] = y[i] * chiBNorm / nB;
    }
}

/*
 * Similar to normalize_y, but components given individually as arrays with n samples
 */
static void NRSur7dq2_normalize_results(
    double normA,       /**< ||vec{chi}_A|| */
    double normB,       /**< ||vec{chi}_B|| */
    gsl_vector **quat,  /**< The four quaternion time-dependent components */
    gsl_vector **chiA,  /**< Time-dependent components of chiA */
    gsl_vector **chiB   /**< Time-dependent components of chiB */
) {
    double nA, nB, nQ;
    int i;
    int n = quat[0]->size;
    double *chiAx = chiA[0]->data;
    double *chiAy = chiA[1]->data;
    double *chiAz = chiA[2]->data;
    double *chiBx = chiB[0]->data;
    double *chiBy = chiB[1]->data;
    double *chiBz = chiB[2]->data;
    double *q0 = quat[0]->data;
    double *qx = quat[1]->data;
    double *qy = quat[2]->data;
    double *qz = quat[3]->data;

    if (normA > 0.0) {
        for (i=0; i<n; i++) {
            nA = sqrt(chiAx[i]*chiAx[i] + chiAy[i]*chiAy[i] + chiAz[i]*chiAz[i]);
            chiAx[i] *= normA / nA;
            chiAy[i] *= normA / nA;
            chiAz[i] *= normA / nA;
        }
    }

    if (normB > 0.0) {
        for (i=0; i<n; i++) {
            nB = sqrt(chiBx[i]*chiBx[i] + chiBy[i]*chiBy[i] + chiBz[i]*chiBz[i]);
            chiBx[i] *= normB / nB;
            chiBy[i] *= normB / nB;
            chiBz[i] *= normB / nB;
        }
    }

    for (i=0; i<n; i++) {
        nQ = sqrt(q0[i]*q0[i] + qx[i]*qx[i] + qy[i]*qy[i] + qz[i]*qz[i]);
        q0[i] /= nQ;
        qx[i] /= nQ;
        qy[i] /= nQ;
        qz[i] /= nQ;
    }
}

/**
 * Transforms chiA and chiB from the coprecessing frame to the coorbital frame
 * using the orbital phase.
 */
static void NRSur7dq2_rotate_spins(
    gsl_vector **chiA,  /**< 3 time-dependent components of chiA in the coorbital frame */
    gsl_vector **chiB,  /**< 3 time-dependent components of chiB in the coorbital frame */
    gsl_vector *phi_vec /**< The time-dependent orbital phase */
) {
    int i;
    int n = phi_vec->size;
    double sp, cp, tmp;
    double *phi = phi_vec->data;

    double *chix = chiA[0]->data;
    double *chiy = chiA[1]->data;
    for (i=0; i<n; i++) {
        cp = cos(phi[i]);
        sp = sin(phi[i]);
        tmp = chix[i];
        chix[i] = tmp*cp + chiy[i]*sp;
        chiy[i] = -1*tmp*sp + chiy[i]*cp;
    }

    chix = chiB[0]->data;
    chiy = chiB[1]->data;
    for (i=0; i<n; i++) {
        cp = cos(phi[i]);
        sp = sin(phi[i]);
        tmp = chix[i];
        chix[i] = tmp*cp + chiy[i]*sp;
        chiy[i] = -1*tmp*sp + chiy[i]*cp;
    }
}

/*
 * When integrating the dynamics ODE, we need to evaluate fits.
 * This helper function computes the fit inputs from the current ODE solution.
 * The spin components of y are in the coprecessing frame, but the spin
 * components of x are in the coorbital frame.
 */
static void NRSur7dq2_ds_fit_x(
    double *x,  /**< Result, length 7 */
    double q,   /**< Mass ratio */
    double *y   /**< [q0, qx, qy, qz, phi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz] */
) {
    double sp = sin(y[4]);
    double cp = cos(y[4]);

    // q
    x[0] = q;

    // chiA, transformed to the coorbital frame
    x[1] = y[5]*cp + y[6]*sp;
    x[2] = -1*y[5]*sp + y[6]*cp;
    x[3] = y[7];

    // chiB, transformed to the coorbital frame
    x[4] = y[8]*cp + y[9]*sp;
    x[5] = -1*y[8]*sp + y[9]*cp;
    x[6] = y[10];
}

/*
 * After evaluating fits, we have many time derivatives, with components in
 * the coorbital frame. dydt has components in the coprecessing frame, so we
 * transform them, and we also compute the time derivative of the unit
 * quaternions using the coprecessing components of Omega
 */
static void NRSur7dq2_assemble_dydt(
    double *dydt,           /**< Result, length 11 */
    double *y,              /**< ODE solution at the current time step */
    double *Omega_coorb_xy, /**< a form of time derivative of the coprecessing frame */
    double omega,           /**< orbital angular frequency in the coprecessing frame */
    double *chiA_dot,       /**< chiA time derivative */
    double *chiB_dot        /**< chiB time derivative */
) {
    double omega_quat_x, omega_quat_y;

    double cp = cos(y[4]);
    double sp = sin(y[4]);

    // Quaternion derivative:
    // Omega = 2 * quat^{-1} * dqdt -> dqdt = 0.5 * quat * omegaQuat, where
    // omegaQuat = [0, Omega_copr_x, Omega_copr_y, 0]
    omega_quat_x = Omega_coorb_xy[0]*cp - Omega_coorb_xy[1]*sp;
    omega_quat_y = Omega_coorb_xy[0]*sp + Omega_coorb_xy[1]*cp;
    dydt[0] = (-0.5)*y[1]*omega_quat_x - 0.5*y[2]*omega_quat_y;
    dydt[1] = (-0.5)*y[3]*omega_quat_y + 0.5*y[0]*omega_quat_x;
    dydt[2] = 0.5*y[3]*omega_quat_x + 0.5*y[0]*omega_quat_y;
    dydt[3] = 0.5*y[1]*omega_quat_y - 0.5*y[2]*omega_quat_x;

    // Orbital phase derivative
    dydt[4] = omega;

    // Spin derivatives
    dydt[5] = chiA_dot[0]*cp - chiA_dot[1]*sp;
    dydt[6] = chiA_dot[0]*sp + chiA_dot[1]*cp;
    dydt[7] = chiA_dot[2];
    dydt[8] = chiB_dot[0]*cp - chiB_dot[1]*sp;
    dydt[9] = chiB_dot[0]*sp + chiB_dot[1]*cp;
    dydt[10] = chiB_dot[2];
}

/**
 * Cubic interpolation of 4 data points
 * This gives a much closer result to scipy.interpolate.InterpolatedUnivariateSpline than using gsl_interp_cspline
 * (see comment in spline_array_interp)
 */
static double cubic_interp(
    double xout,    /**< The target x value */
    double *x,      /**< The x values of the points to interpolate. Length 4, must be increasing. */
    double *y       /**< The y values of the points to interpolate. Length 4. */
) {
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *interpolant = gsl_spline_alloc(gsl_interp_polynomial, 4);
    gsl_spline_init(interpolant, x, y, 4);
    double res = gsl_spline_eval(interpolant, xout, acc);
    gsl_spline_free(interpolant);
    gsl_interp_accel_free(acc);
    return res;
}

/**
 * Do cubic spline interpolation using a gsl_interp_cspline.
 * Results differ from scipy.interpolate.InterpolatedUnivariateSpline due to different boundary conditions.
 * This difference leads to small differences between this implementation of NRSur7dq2 and the python
 * implementation, especially near the start and end of the waveform.
 */
static gsl_vector *spline_array_interp(
    gsl_vector *xout,   /**< The vector of points onto which we want to interpolate. */
    gsl_vector *x,      /**< The x values of the data to interpolate. */
    gsl_vector *y       /**< The y values of the data to interpolate. */
) {
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, x->size);
    gsl_spline_init(spline, x->data, y->data, x->size);

    gsl_vector *res = gsl_vector_alloc(xout->size);
    double tmp;
    for (size_t i=0; i<xout->size; i++) {
        tmp = gsl_spline_eval(spline, gsl_vector_get(xout, i), acc);
        gsl_vector_set(res, i, tmp);
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return res;
}

/**
 * Computes the orbital angular frequency at a dynamics node.
 */
static double NRSur7dq2_get_omega(
    size_t node_index,  /**< The index of the dynamics node. */
    double q,           /**< The mass ratio. */
    double *y0          /**< The value of the ODE state y = [q0, qx, qy, qz, orbphase, chiAx, chiAy, chiAz,
                                                            chiBx, chiBy, chiBz] */
){
    double x[7];
    NRSur7dq2_ds_fit_x(x, q, y0);
    FitData *data = (&__lalsim_NRSur7dq2_data)->ds_node_data[node_index]->omega_data;
    double omega = NRSur7dq2_eval_fit(data, x);
    return omega;
}

/**
 * Computes a reference time from a reference orbital angular frequency.
 */
static double NRSur7dq2_get_t_ref(
    double omega_ref,   /**< reference orbital angular frequency */
    double q,           /**< mass ratio */
    double *chiA0,      /**< chiA at reference point */
    double *chiB0,      /**< chiB at reference point */
    double *q_ref,      /**< coprecessing frame quaternion at reference point */
    double phi_ref      /**< orbital phase at reference point */
) {
    if (fabs(omega_ref) < 1.e-10) {
        XLAL_PRINT_WARNING("WARNING: Treating omega_ref = 0 as a flag to use t_ref = t_0 = -4500M");
        return -4499.99999999; // The first time node is ever so slightly larger than -4500; this avoids out-of-range issues
    }

    // omega_ref <= 0.201 guarantees omega is smaller than the final and merger omegas for all masses and spins.
    if (omega_ref > 0.201) XLAL_ERROR_REAL8(XLAL_FAILURE, "Reference frequency omega_ref=%0.4f > 0.2, too large!\n", omega_ref);

    double y0[11];
    int i;
    for (i=0; i<4; i++) y0[i] = q_ref[i];
    y0[4] = phi_ref;
    for (i=0; i<3; i++) {
        y0[i+5] = chiA0[i];
        y0[i+8] = chiB0[i];
    }

    double omega_min = NRSur7dq2_get_omega(0, q, y0);
    if (omega_ref < omega_min) {
        XLAL_ERROR_REAL8(XLAL_FAILURE, "Got omega_ref=%0.4f smaller than omega0=%0.4f for this configuration!", omega_ref, omega_min);
    }

    // i0=0 is a lower bound; find the first index where omega > omega_ref, and the previous index will have omega <= omega_ref.
    size_t imax = 1;
    double omega_max = NRSur7dq2_get_omega(imax, q, y0);
    gsl_vector *t_ds = (&__lalsim_NRSur7dq2_data)->t_ds;
    while ((omega_max <= omega_ref) && (imax < t_ds->size)) {
        imax += 1;
        omega_min = omega_max;
        omega_max = NRSur7dq2_get_omega(imax, q, y0);
    }

    if (omega_max <= omega_ref) {
        XLAL_ERROR_REAL8(XLAL_FAILURE, "Tried all nodes and still have omega=%0.4f <= omega_ref=%0.4f\n", omega_min, omega_ref);
    }

    // Do a linear interpolation between omega_min and omega_max
    double t_min = gsl_vector_get(t_ds, imax-1);
    double t_max = gsl_vector_get(t_ds, imax);
    double t_ref = (t_min * (omega_max - omega_ref) + t_max * (omega_ref - omega_min)) / (omega_max - omega_min);
    return t_ref;
}

/**
 * Compute dydt at a given dynamics node, where y is the numerical solution to the dynamics ODE.
 */
static void NRSur7dq2_get_time_deriv_from_index(
    double *dydt,   /**< Output: dy/dt evaluated at the ODE time node with index i0. Must have space for 11 entries. */
    int i0,         /**< Time node index. i0=-1, -2, and -3 are used for time nodes 1/2, 3/2, and 5/2 respectively. */
    double q,       /**< Mass ratio */
    double *y       /**< Current ODE state: [q0, qx, qy, qz, orbphase, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz] */
) {
    // Setup fit variables
    double x[7];
    NRSur7dq2_ds_fit_x(x, q, y);

    // Get fit data
    DynamicsNodeFitData *ds_node;
    if (i0 >= 0) {
        ds_node = (&__lalsim_NRSur7dq2_data)->ds_node_data[i0];
    } else {
        ds_node = (&__lalsim_NRSur7dq2_data)->ds_half_node_data[-1*i0 - 1];
    }

    // Evaluate fits
    double omega, Omega_coorb_xy[2], chiA_dot[3], chiB_dot[3];
    omega = NRSur7dq2_eval_fit(ds_node->omega_data, x);
    NRSur7dq2_eval_vector_fit(Omega_coorb_xy, ds_node->omega_copr_data, x);
    NRSur7dq2_eval_vector_fit(chiA_dot, ds_node->chiA_dot_data, x);
    NRSur7dq2_eval_vector_fit(chiB_dot, ds_node->chiB_dot_data, x);
    NRSur7dq2_assemble_dydt(dydt, y, Omega_coorb_xy, omega, chiA_dot, chiB_dot);
}

/**
 * Compute dydt at any time by evaluating dydt at 4 nearby dynamics nodes and
 * using cubic spline interpolation to evaluate at the desired time.
 */
static void NRSur7dq2_get_time_deriv(
    double *dydt,   /**< Output: dy/dt evaluated at time t. Must have space for 11 entries. */
    double t,       /**< Time at which the ODE should be evaluated. */
    double q,       /**< Mass ratio */
    double *y       /**< Current ODE state */
) {
    // Make sure we are within the valid range
    gsl_vector *t_ds = (&__lalsim_NRSur7dq2_data)->t_ds;
    int i1 = 0;
    int imax = t_ds->size - 1;
    double t1 = gsl_vector_get(t_ds, i1);
    double tmax = gsl_vector_get(t_ds, imax);
    if (t < t1 || t > tmax) {
        XLAL_ERROR_VOID(XLAL_FAILURE, "Tried to get dydt at t=%0.12f, not in valid range [%0.12f, %0.12f]\n", t, t1, tmax);
    }

    double t2 = gsl_vector_get(t_ds, i1+1);
    // Put t2 slightly larger than t
    while (t2 < t) {
        i1 += 1;
        t1 = t2;
        t2 = gsl_vector_get(t_ds, i1+1);
    }

    // Do cubic spline interpolation using 4 data points, cenetered if possible
    int i0 = i1-1;
    if (i0 < 0) i0 = 0;
    if (i0 > imax-3) i0 = imax-3;
    double times[4], derivs[4], dydt0[11], dydt1[11], dydt2[11], dydt3[11];
    int j;
    for (j=0; j<4; j++) times[j] = gsl_vector_get(t_ds, i0+j);
    NRSur7dq2_get_time_deriv_from_index(dydt0, i0, q, y);
    NRSur7dq2_get_time_deriv_from_index(dydt1, i0+1, q, y);
    NRSur7dq2_get_time_deriv_from_index(dydt2, i0+2, q, y);
    NRSur7dq2_get_time_deriv_from_index(dydt3, i0+3, q, y);

    for (j=0; j<11; j++) {
        derivs[0] = dydt0[j];
        derivs[1] = dydt1[j];
        derivs[2] = dydt2[j];
        derivs[3] = dydt3[j];
        dydt[j] = cubic_interp(t, times, derivs);
    }

}

/**
 * Initialize the dynamics ODE at a dynamics node.
 * Given t_ref, finds the nearest dynamics node and integrates the initial conditions at t_ref
 * a tiny bit to obtain the ODE state at the dynamics node.
 */
static int NRSur7dq2_initialize_at_dynamics_node(
    double *dynamics_data,  /**< ODE output */
    double t_ref,           /**< reference time. t_ds[i0] will be close to t_ref. */
    double q,               /**< mass ratio */
    double *chiA0,          /**< chiA at t_ref. */
    double *chiB0,          /**< chiB at t_ref. */
    double phi_ref,         /**< orbital phase at t_ref. */
    double *q_ref,          /**< quaternion at t_ref. */
    double normA,           /**< |chiA| */
    double normB            /**< |chiB| */
) {
    int imin, imax, i0, j;
    double tmin, tmax, t0, dt, y_ref[11], dydt_ref[11], *node_data;
    gsl_vector *t_ds;

    // First find i0, the closest dynamics time node to time t_ref, with a binary search
    t_ds = (&__lalsim_NRSur7dq2_data)->t_ds;
    imin = 0;
    imax = t_ds->size - 1;
    tmin = gsl_vector_get(t_ds, imin);
    tmax = gsl_vector_get(t_ds, imax);
    while (imax - imin > 1) {
        i0 = (imax + imin)/2;
        t0 = gsl_vector_get(t_ds, i0);
        if (t0 > t_ref) {
            imax = i0;
            tmax = t0;
        } else {
            imin = i0;
            tmin = t0;
        }
    }

    // Now we have tmin <= t_ref < tmax, and imax = imin+1.
    // Step towards the closest of tmin, tmax.
    if (fabs(t_ref - tmin) < fabs(tmax - t_ref)) {
        i0 = imin;
        dt = tmin - t_ref;
    } else {
        i0 = imax;
        dt = tmax - t_ref;
    }

    // Setup y at t_ref
    for (j=0; j<4; j++) y_ref[j] = q_ref[j];
    y_ref[4] = phi_ref;
    for (j=0; j<3; j++) {
        y_ref[5+j] = chiA0[j];
        y_ref[8+j] = chiB0[j];
    }

    // Compute dydt at t_ref
    NRSur7dq2_get_time_deriv(dydt_ref, t_ref, q, y_ref);

    // modify y_ref to be y at t0 = t_ds[i0]
    for (j=0; j<11; j++) y_ref[j] += dt * dydt_ref[j];
    NRSur7dq2_normalize_y(normA, normB, y_ref);

    // transfer results to ODE output data structures
    node_data = dynamics_data + i0*11; // Point to the output at i0
    for (j=0; j<11; j++) node_data[j] = y_ref[j];

    return i0;
}

/**
 * Initializes the AB4 ODE system from the surrogate start time by taking 3 RK4 steps,
 * making use of the three additional half-time-step nodes for the RK4 substeps.
 * This is the recommended way to initialize the AB4 ODE system - the additional nodes
 * are there to increase accuracy during initialization.
 * The drawback is that we are forced to accept fRef to be the
 * surrogate start frequency, which will depend on masses and spins.
 * The ODE system is for the vector of 11 quantities:
 * [q0, qx, qy, qz, varphi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
 * where (q0, qx, qy, qz) is the coprecessing frame quaternion, varphi is the
 * orbital phase, and the spin components are in the coprecessing frame.
 */
static void NRSur7dq2_initialize_RK4_with_half_nodes(
    double *dynamics_data,  /**< A pointer to the start of the ODE output; the first 11 entries
                                 (for node 0) should have already been set.*/
    double *time_steps,     /**< Output: first three time steps. Should have size 3. */
    double *dydt0,          /**< Output: dydt at node 0. Should have size 11. */
    double *dydt1,          /**< Output: dydt at node 1. Should have size 11. */
    double *dydt2,          /**< Output: dydt at node 2. Should have size 11. */
    double *dydt3,          /**< Output: dydt at node 3. Should have size 11. */
    double normA,           /**< |chiA| */
    double normB,           /**< |chiB| */
    double q                /**< mass ratio */
) {
    gsl_vector *t_ds = (&__lalsim_NRSur7dq2_data)->t_ds;
    double t1, t2;
    int i, j;
    double k2[11], k3[11], k4[11]; // dydt for all but the first RK4 substep (which will be one of dydt0, dydt1, dydt2)
    double *node_data;
    double y_tmp[11];
    double *dydt[3];
    dydt[0] = dydt0;
    dydt[1] = dydt1;
    dydt[2] = dydt2;

    // Setup time steps
    t1 = gsl_vector_get(t_ds, 0);
    for (i=0; i<3; i++) {
        t2 = gsl_vector_get(t_ds, i+1);
        time_steps[i] = t2 - t1;
        t1 = t2;
    }

    // Three steps of RK4.
    for (i=0; i<3; i++ ) {

        // Initial substep
        node_data = dynamics_data + 11*i;
        NRSur7dq2_get_time_deriv_from_index(dydt[i], i, q, node_data);

        // Next substep: evaluate at the half node (by using -1, -2, or -3 for i=0, 1, or 2)
        for (j=0; j<11; j++) {
            y_tmp[j] = node_data[j] + 0.5*time_steps[i]*dydt[i][j];
        }
        NRSur7dq2_get_time_deriv_from_index(k2, -1-i, q, y_tmp);

        // Next substep: also at the half node, but update y_tmp
        for (j=0; j<11; j++) {
            y_tmp[j] = node_data[j] + 0.5*time_steps[i]*k2[j];
        }
        NRSur7dq2_get_time_deriv_from_index(k3, -1-i, q, y_tmp);

        // Final substep: evaluate at the next node
        for (j=0; j<11; j++) {
            y_tmp[j] = node_data[j] + time_steps[i]*k3[j];
        }
        NRSur7dq2_get_time_deriv_from_index(k4, i+1, q, y_tmp);

        // Compute the RK4 expression for the next node, normalize, and store the data
        for (j=0; j<11; j++) {
            y_tmp[j] = node_data[j] + (time_steps[i]/6.0)*(dydt[i][j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        NRSur7dq2_normalize_y(normA, normB, y_tmp);
        node_data = node_data + 11;
        for (j=0; j<11; j++) {
            node_data[j] = y_tmp[j];
        }
    }

    // Finally, we need to compute dydt3;
    node_data = dynamics_data + 33;
    NRSur7dq2_get_time_deriv_from_index(dydt3, 3, q, node_data);
}

/**
 * Initializes the AB4 ODE system from a single arbitrary dynamics node by taking 3 RK4 steps.
 * Compared to NRSur7dq2_initialize_RK4_with_half_nodes, this may be slightly less accurate
 * and slightly more time consuming as we must interpolate many dynamics node fit evaluations,
 * but this is more flexible as we can initialize from any node instead of just node 0.
 */
static int NRSur7dq2_initialize_RK4(
    double *dynamics_data,  /**< A pointer to the start of the ODE output.
                                 Entries with i0*11 <= i < (i0+1)*11 should have already be computed. */
    double *time_steps,     /**< Output: first three time steps. Should have size 3. */
    double *dydt0,          /**< Output: dydt at node i0 + 0. Should have size 11. */
    double *dydt1,          /**< Output: dydt at node i0 + 1. Should have size 11. */
    double *dydt2,          /**< Output: dydt at node i0 + 2. Should have size 11. */
    double *dydt3,          /**< Output: dydt at node i0 + 3. Should have size 11. */
    double normA,           /**< |chiA| */
    double normB,           /**< |chiB| */
    double q,               /**< mass ratio */
    int i0                  /**< the node that is already initialized */
) {
    gsl_vector *t_ds = (&__lalsim_NRSur7dq2_data)->t_ds;
    double t1, t2;
    int i, j;
    double k2[11], k3[11], k4[11]; // dydt for all but the first RK4 substep (which will be one of dydt0, dydt1, dydt2, dydt3)
    double *node_data;
    double y_tmp[11];
    double *dydt[4];
    dydt[0] = dydt0;
    dydt[1] = dydt1;
    dydt[2] = dydt2;
    dydt[3] = dydt3;
    int i_start;

    // If possible, do 3 steps backwards with RK4. Otherwise, do 3 steps forwards.
    if (i0 > 2) {
        i_start = i0 - 3;

        // Setup time steps
        t1 = gsl_vector_get(t_ds, i_start);
        for (i=0; i<3; i++) {
            t2 = gsl_vector_get(t_ds, i_start+i+1);
            time_steps[i] = t2 - t1;
            t1 = t2;
        }

        // Three steps of RK4 BACKWARDS, so include a minus sign beside all dydts.
        for (i=0; i<3; i++ ) {

            // Initial substep
            node_data = dynamics_data + 11*(i0-i);
            NRSur7dq2_get_time_deriv_from_index(dydt[3-i], i0-i, q, node_data);

            // Next substep: evaluate halfway to the next timestep
            t1 = 0.5*(gsl_vector_get(t_ds, i0-i) + gsl_vector_get(t_ds, i0-i-1));
            for (j=0; j<11; j++) {
                y_tmp[j] = node_data[j] - 0.5*time_steps[2-i]*dydt[3-i][j];
            }
            NRSur7dq2_get_time_deriv(k2, t1, q, y_tmp);

            // Next substep: also halfway, but update y_tmp
            for (j=0; j<11; j++) {
                y_tmp[j] = node_data[j] - 0.5*time_steps[2-i]*k2[j];
            }
            NRSur7dq2_get_time_deriv(k3, t1, q, y_tmp);

            // Final substep: evaluate at the next node
            for (j=0; j<11; j++) {
                y_tmp[j] = node_data[j] - time_steps[2-i]*k3[j];
            }
            NRSur7dq2_get_time_deriv_from_index(k4, i0-i-1, q, y_tmp);

            // Compute the RK4 expression for the next node, normalize, and store the data
            for (j=0; j<11; j++) {
                y_tmp[j] = node_data[j] - (time_steps[2-i]/6.0)*(dydt[3-i][j] + 2*k2[j] + 2*k3[j] + k4[j]);
            }
            NRSur7dq2_normalize_y(normA, normB, y_tmp);
            node_data = node_data - 11;
            for (j=0; j<11; j++) {
                node_data[j] = y_tmp[j];
            }
        }
        // Finally, we need to compute dydt0;
        node_data = dynamics_data + (i0 - 3)*11;
        NRSur7dq2_get_time_deriv_from_index(dydt0, i0-3, q, node_data);
    } else {
        i_start = i0;

        // Setup time steps
        t1 = gsl_vector_get(t_ds, i0);
        for (i=0; i<3; i++) {
            t2 = gsl_vector_get(t_ds, i0+i+1);
            time_steps[i] = t2 - t1;
            t1 = t2;
        }

        // Three steps of RK4.
        for (i=0; i<3; i++ ) {

            // Initial substep
            node_data = dynamics_data + 11*(i+i0);
            NRSur7dq2_get_time_deriv_from_index(dydt[i], i0+i, q, node_data);

            // Next substep: evaluate halfway to the next timestep
            t1 = 0.5*(gsl_vector_get(t_ds, i+i0) + gsl_vector_get(t_ds, i+i0+1));
            for (j=0; j<11; j++) {
                y_tmp[j] = node_data[j] + 0.5*time_steps[i]*dydt[i][j];
            }
            NRSur7dq2_get_time_deriv(k2, t1, q, y_tmp);

            // Next substep: also halfway, but update y_tmp
            for (j=0; j<11; j++) {
                y_tmp[j] = node_data[j] + 0.5*time_steps[i]*k2[j];
            }
            NRSur7dq2_get_time_deriv(k3, t1, q, y_tmp);

            // Final substep: evaluate at the next node
            for (j=0; j<11; j++) {
                y_tmp[j] = node_data[j] + time_steps[i]*k3[j];
            }
            NRSur7dq2_get_time_deriv_from_index(k4, i0+i+1, q, y_tmp);

            // Compute the RK4 expression for the next node, normalize, and store the data
            for (j=0; j<11; j++) {
                y_tmp[j] = node_data[j] + (time_steps[i]/6.0)*(dydt[i][j] + 2*k2[j] + 2*k3[j] + k4[j]);
            }
            NRSur7dq2_normalize_y(normA, normB, y_tmp);
            node_data = node_data + 11;
            for (j=0; j<11; j++) {
                node_data[j] = y_tmp[j];
            }
        }
        // Finally, we need to compute dydt3;
        node_data = dynamics_data + (i0 + 3)*11;
        NRSur7dq2_get_time_deriv_from_index(dydt3, i0+3, q, node_data);
    }

    return i_start;
}

/**
 * Integrates the AB4 ODE system in time forwards, and backwards if needed.
 * The system should have already been initialized at 4 adjacent dynamics nodes.
 * Output is a flattened array where entries i*11 <= j < (i+1)*11 are
 * [q0, qx, qy, qz, varphi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz] at dynamics node i,
 * where (q0, qx, qy, qz) is the coprecessing frame quaternion, varphi is the
 * orbital phase, and the spin components are in the coprecessing frame.
 */
static void NRSur7dq2_integrate_AB4(
    double *dynamics_data,  /**< ODE output */
    double *time_steps,     /**< The first three time steps beginning at i_start. */
    double *dydt0,          /**< dydt at node i_start */
    double *dydt1,          /**< dydt at node i_start + 1 */
    double *dydt2,          /**< dydt at node i_start + 2 */
    double *dydt3,          /**< dydt at node i_start + 3 */
    double normA,           /**< |chiA| */
    double normB,           /**< |chiB| */
    double q,               /**< mass ratio */
    int i_start             /**< nodes i_start through i_start+3 are already initialized */
) {
    // At each step, we will use the 4 most recent nodes where we have the solution, and integrate up to the first new node.
    double dt1, dt2, dt3, dt4; // Time steps between the 4 known nodes, as well as up to the first new node (dt4)
    double k1[11], k2[11], k3[11], k4[11]; // dydt(y(t); t) evaluated at the 4 known nodes
    double ynext[11]; // Temporary solution
    double *node_data; // Pointer for writing output
    int i, j;

    gsl_vector *t_ds = (&__lalsim_NRSur7dq2_data)->t_ds;

    // Initialize to integrate forwards
    dt1 = time_steps[0];
    dt2 = time_steps[1];
    dt3 = time_steps[2];
    for (j=0; j<11; j++) {
        k1[j] = dydt0[j];
        k2[j] = dydt1[j];
        k3[j] = dydt2[j];
        k4[j] = dydt3[j];
    }

    // Integrate forwards
    node_data = dynamics_data + 11*(i_start + 3); // Point to latest known solution
    for (i=i_start+4; i< (int)t_ds->size; i++) { // i indexes the output

        // Compute dy = dydt*dt, write it in ynext
        dt4 = gsl_vector_get(t_ds, i) - gsl_vector_get(t_ds, i-1);
        NRSur_ab4_dy(ynext, k1, k2, k3, k4, dt1, dt2, dt3, dt4, 11);

        // Add the latest known node, to obtain y at the next node
        for (j=0; j<11; j++) ynext[j] += node_data[j];

        // Normalize and write output
        NRSur7dq2_normalize_y(normA, normB, ynext);
        node_data += 11;
        for (j=0; j<11; j++) node_data[j] = ynext[j];

        // Setup for next step
        if (i < (int) t_ds->size - 1) {
            dt1 = dt2;
            dt2 = dt3;
            dt3 = dt4;
            for (j=0; j<11; j++) {
                k1[j] = k2[j];
                k2[j] = k3[j];
                k3[j] = k4[j];
            }
            NRSur7dq2_get_time_deriv_from_index(k4, i, q, node_data);
        }
    }

    // Initialize to integrate backwards
    dt1 = -1 * time_steps[2];
    dt2 = -1 * time_steps[1];
    dt3 = -1 * time_steps[0];
    for (j=0; j<11; j++) {
        k1[j] = dydt3[j];
        k2[j] = dydt2[j];
        k3[j] = dydt1[j];
        k4[j] = dydt0[j];
    }

    // Integrate backwards
    node_data = dynamics_data + 11*(i_start); // Point to earliest known solution
    for (i=i_start-1; i>=0; i--) { // i indexes the output

        // Compute dy = dydt*dt, write it in ynext
        dt4 = gsl_vector_get(t_ds, i) - gsl_vector_get(t_ds, i+1);
        NRSur_ab4_dy(ynext, k1, k2, k3, k4, dt1, dt2, dt3, dt4, 11);

        // Add the earliest known node, to obtain y at the previous
        for (j=0; j<11; j++) ynext[j] += node_data[j];

        // Normalize and write output
        NRSur7dq2_normalize_y(normA, normB, ynext);
        node_data -= 11;
        for (j=0; j<11; j++) node_data[j] = ynext[j];

        // Setup for next step
        if (i > 0) {
            dt1 = dt2;
            dt2 = dt3;
            dt3 = dt4;
            for (j=0; j<11; j++) {
                k1[j] = k2[j];
                k2[j] = k3[j];
                k3[j] = k4[j];
            }
            NRSur7dq2_get_time_deriv_from_index(k4, i, q, node_data);
        }
    }
}

/**
 * Evaluates a single NRSur7dq2 coorbital waveoform data piece.
 * The dynamics ODE must have already been solved, since this requires the
 * spins evaluated at all of the empirical nodes for this waveform data piece.
 */
static void NRSur7dq2_eval_data_piece(
    gsl_vector *result, /**< Output: Should have already been assigned space */
    double q,           /**< Mass ratio */
    gsl_vector **chiA,  /**< 3 gsl_vector *s, one for each (coorbital) component */
    gsl_vector **chiB,  /**< similar to chiA */
    WaveformDataPiece *data /**< The data piece to evaluate */
) {
    gsl_vector *nodes = gsl_vector_alloc(data->n_nodes);
    double x[7];
    int i, j, node_index;

    // Evaluate the fits at the empirical nodes, using the spins at the empirical node times
    x[0] = q;
    for (i=0; i<data->n_nodes; i++) {
        node_index = gsl_vector_long_get(data->empirical_node_indices, i);
        for (j=0; j<3; j++) {
            x[1+j] = gsl_vector_get(chiA[j], node_index);
            x[4+j] = gsl_vector_get(chiB[j], node_index);
        }
        gsl_vector_set(nodes, i, NRSur7dq2_eval_fit(data->fit_data[i], x));
    }

    // Evaluate the empirical interpolant
    gsl_blas_dgemv(CblasTrans, 1.0, data->empirical_interpolant_basis, nodes, 0.0, result);

    gsl_vector_free(nodes);
}

/************************ Main Waveform Generation Routines ***********/

/**
 * This is the main NRSur7dq2 dynamics surrogate integration routine.
 * Given omega_ref and the system parameters at omega_ref, we find the corresponding t_ref
 * and first (if needed) take a tiny time step and evolve the system to the nearest
 * dynamics node.
 * We then initialize the AB4 ODE system at 4 consecutive dynamics nodes by taking 3 RK4 steps.
 * Finally, we integrate forwards (and backwards if needed) to obtain the solution at all
 * dynamics nodes.
 * Output is a flattened array where entries i*11 <= j < (i+1)*11 are
 * [q0, qx, qy, qz, varphi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz] at dynamics node i,
 * where (q0, qx, qy, qz) is the coprecessing frame quaternion, varphi is the
 * orbital phase, and the spin components are in the coprecessing frame.
 */
static int NRSur7dq2_IntegrateDynamics(
    double *dynamics_data,  /**< Output: length (n_dynamics_nodes * 11) */
    double q,               /**< Mass ratio mA / mB */
    double *chiA0,          /**< chiA at the reference point */
    double *chiB0,          /**< chiB at the reference point */
    double omega_ref,       /**< orbital angular frequency at the reference point */
    double phi_ref,         /**< orbital phase at the reference point */
    double *q_ref           /**< coprecessing quaterion at the reference point */
) {
    double normA = sqrt(chiA0[0]*chiA0[0] + chiA0[1]*chiA0[1] + chiA0[2]*chiA0[2]);
    double normB = sqrt(chiB0[0]*chiB0[0] + chiB0[1]*chiB0[1] + chiB0[2]*chiB0[2]);

    // Sanity checks and warnings
    if (q < NRSUR7DQ2_Q_MIN) XLAL_ERROR_REAL8(XLAL_FAILURE, "Invalid mass ratio q = %0.4f < 1\n", q);
    if (q > NRSUR7DQ2_Q_MAX) XLAL_ERROR_REAL8(XLAL_FAILURE, "Invalid mass ratio q = %0.4f > 2\n", q);
    if (normA > 1.001) XLAL_ERROR_REAL8(XLAL_FAILURE, "Invalid spin magnitude |chiA| = %0.4f > 1\n", normA);
    if (normB > 1.001) XLAL_ERROR_REAL8(XLAL_FAILURE, "Invalid spin magnitude |chiB| = %0.4f > 1\n", normB);
    if (normA > NRSUR7DQ2_CHI_MAX) XLAL_PRINT_WARNING("Extrapolating to spin magnitude |chiA| = %0.4f > 0.8\n", normA);
    if (normB > NRSUR7DQ2_CHI_MAX) XLAL_PRINT_WARNING("Extrapolating to spin magnitude |chiB| = %0.4f > 0.8\n", normB);

    double t_ref = NRSur7dq2_get_t_ref(omega_ref, q, chiA0, chiB0, q_ref, phi_ref);
    if (XLAL_IS_REAL8_FAIL_NAN(t_ref)) XLAL_ERROR_REAL8(XLAL_FAILURE, "Failed to determine t_ref");

    // Initialize the ODE system by stepping to a dynamics node indexed by i0
    int i0 = NRSur7dq2_initialize_at_dynamics_node(dynamics_data, t_ref, q, chiA0, chiB0, phi_ref, q_ref, normA, normB);

    // To initialize the RK4 ODE integration scheme, we need y to be evaluated at 4 consecutive nodes. Currently we have 1.
    // The original method assumed we always use t_ref=-4500M and took 3 steps using RK4, making use of the time nodes 1/2, 3/2, and 5/2.
    // If i0=0, we can use that scheme and obtain y at nodes 0, 1, 2, and 3. Otherwise, obtain 4 consecutive nodes starting at i_start,
    // where i_start = max(0, i0 - 3).
    int i_start;
    double time_steps[3];
    double dydt0[11], dydt1[11], dydt2[11], dydt3[11];

    if (i0 == 0) {
        NRSur7dq2_initialize_RK4_with_half_nodes(dynamics_data, time_steps, dydt0, dydt1, dydt2, dydt3, normA, normB, q);
        i_start = 0;
    } else {
        i_start = NRSur7dq2_initialize_RK4(dynamics_data, time_steps, dydt0, dydt1, dydt2, dydt3, normA, normB, q, i0);
    }

    // Now that we have 4 consecutive evaluated nodes beginning at i_start, we can integrate from i_start backwards to 0
    // with AB4 as well as from i_start+3 to the end
    NRSur7dq2_integrate_AB4(dynamics_data, time_steps, dydt0, dydt1, dydt2, dydt3, normA, normB, q, i_start);
    return XLAL_SUCCESS;
}


/* This is the core function of the NRSur7dq2 model.
 * It evaluates the model, and returns waveform modes in the inertial frame, sampled on t_coorb.
 * When using a custom ModeArray the user must explicitly supply all modes, -ve and +ve m-modes.
 * Note that the co-orbital frame modes are specified in ModeArray, not the inertial frame modes.
 */
static void NRSur7dq2_core(
    MultiModalWaveform **h, /**< Output. Dimensionless waveform modes sampled on t_coorb */
    double q,               /**< Mass ratio mA / mB */
    double *chiA0,          /**< chiA at the reference point */
    double *chiB0,          /**< chiB at the reference point */
    double omega_ref,       /**< orbital angular frequency at the reference point */
    double phi_ref,         /**< orbital phase at the reference point */
    double *q_ref,          /**< coprecessing quaterion at the reference point */
    LALValue* ModeArray     /**< Container for the ell and m co-orbital modes to generate */
) {
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&NRSur7dq2_is_initialized, NRSur7dq2_Init_LALDATA);
#else
    NRSur7dq2_Init_LALDATA();
#endif

    // Make sure we didn't request any unavailable modes
    int ell, m;
    for (ell=NRSUR7DQ2_LMAX+1; ell <= LAL_SIM_L_MAX_MODE_ARRAY; ell++) {
        for (m=-ell; m<=ell; m++) {
            if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, m) == 1) {
                XLAL_ERROR_VOID(XLAL_FAILURE, "A requested mode has ell larger than available");
            }
        }
    }

    // time arrays
    gsl_vector *t_ds = (&__lalsim_NRSur7dq2_data)->t_ds;
    gsl_vector *t_coorb = (&__lalsim_NRSur7dq2_data)->t_coorb;
    int n_ds = t_ds->size;
    int n_coorb = t_coorb->size;

    // Get dynamics
    double *dynamics_data = XLALCalloc(n_ds * 11, sizeof(double));

    int ret = NRSur7dq2_IntegrateDynamics(dynamics_data, q, chiA0, chiB0, omega_ref, phi_ref, q_ref);
    if(ret != XLAL_SUCCESS) XLAL_ERROR_VOID(XLAL_FAILURE, "Failed to integrate dynamics");

    // Put output into appropriate vectors
    int i, j;
    gsl_vector *dynamics[11];
    for (j=0; j<11; j++) {
        dynamics[j] = gsl_vector_alloc(n_ds);
    }
    for (i=0; i<n_ds; i++) {
        for (j=0; j<11; j++) {
            gsl_vector_set(dynamics[j], i, dynamics_data[11*i + j]);
        }
    }

    XLALFree(dynamics_data);

    // Interpolate onto the coorbital time grid
    gsl_vector *quat_coorb[4], *chiA_coorb[3], *chiB_coorb[3];
    quat_coorb[0] = spline_array_interp(t_coorb, t_ds, dynamics[0]);
    quat_coorb[1] = spline_array_interp(t_coorb, t_ds, dynamics[1]);
    quat_coorb[2] = spline_array_interp(t_coorb, t_ds, dynamics[2]);
    quat_coorb[3] = spline_array_interp(t_coorb, t_ds, dynamics[3]);
    gsl_vector *phi_coorb = spline_array_interp(t_coorb, t_ds, dynamics[4]);
    chiA_coorb[0] = spline_array_interp(t_coorb, t_ds, dynamics[5]);
    chiA_coorb[1] = spline_array_interp(t_coorb, t_ds, dynamics[6]);
    chiA_coorb[2] = spline_array_interp(t_coorb, t_ds, dynamics[7]);
    chiB_coorb[0] = spline_array_interp(t_coorb, t_ds, dynamics[8]);
    chiB_coorb[1] = spline_array_interp(t_coorb, t_ds, dynamics[9]);
    chiB_coorb[2] = spline_array_interp(t_coorb, t_ds, dynamics[10]);

    for (j=0; j<11; j++) {
        gsl_vector_free(dynamics[j]);
    }

    // Normalize spins after interpolation
    double normA = sqrt(chiA0[0]*chiA0[0] + chiA0[1]*chiA0[1] + chiA0[2]*chiA0[2]);
    double normB = sqrt(chiB0[0]*chiB0[0] + chiB0[1]*chiB0[1] + chiB0[2]*chiB0[2]);
    NRSur7dq2_normalize_results(normA, normB, quat_coorb, chiA_coorb, chiB_coorb);

    // Transform spins from coprecessing frame to coorbital frame for use in coorbital waveform surrogate
    NRSur7dq2_rotate_spins(chiA_coorb, chiB_coorb, phi_coorb);

    // Evaluate the coorbital waveform surrogate
    MultiModalWaveform *h_coorb = NULL;
    MultiModalWaveform_Init(&h_coorb, NRSUR7DQ2_LMAX, n_coorb);
    gsl_vector *data_piece_eval = gsl_vector_alloc(n_coorb);
    WaveformDataPiece *data_piece_data;
    int i0; // for indexing the (ell, m=0) mode, such that the (ell, m) mode is index (i0 + m).
    WaveformFixedEllModeData *ell_data;
    for (ell=2; ell<=NRSUR7DQ2_LMAX; ell++) {
        ell_data = (&__lalsim_NRSur7dq2_data)->coorbital_mode_data[ell - 2];
        i0 = ell*(ell+1) - 4;

        // m=0
        if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, 0) != 1) {
          XLAL_PRINT_INFO("SKIPPING (%i, 0) co-orbital mode", ell);
        }
        else {
            XLAL_PRINT_INFO("Generating (%i, 0) co-orbital mode", ell);
            data_piece_data = ell_data->m0_real_data;
            NRSur7dq2_eval_data_piece(data_piece_eval, q, chiA_coorb, chiB_coorb, data_piece_data);
            gsl_vector_add(h_coorb->modes_real_part[i0], data_piece_eval);

            data_piece_data = ell_data->m0_imag_data;
            NRSur7dq2_eval_data_piece(data_piece_eval, q, chiA_coorb, chiB_coorb, data_piece_data);
            gsl_vector_add(h_coorb->modes_imag_part[i0], data_piece_eval);
        }

        // Other modes
        for (m=1; m<=ell; m++) {
            if ((XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, m) != 1) &&
                (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, -m) != 1)) {
                XLAL_PRINT_INFO("SKIPPING (%i, %i) co-orbital mode", ell, m);
                continue;
            }
            XLAL_PRINT_INFO("Generating (%i, +/- %i) co-orbital modes", ell, m);

            // h^{ell, m} = X_plus + X_minus
            // h^{ell, -m} = (X_plus - X_minus)* <- complex conjugate

            // Re[X_plus] gets added to both Re[h^{ell, m}] and Re[h^{ell, -m}]
            data_piece_data = ell_data->X_real_plus_data[m-1];
            NRSur7dq2_eval_data_piece(data_piece_eval, q, chiA_coorb, chiB_coorb, data_piece_data);
            gsl_vector_add(h_coorb->modes_real_part[i0+m], data_piece_eval);
            gsl_vector_add(h_coorb->modes_real_part[i0-m], data_piece_eval);

            // Re[X_minus] gets added to Re[h^{ell, m}] and subtracted from Re[h^{ell, -m}]
            data_piece_data = ell_data->X_real_minus_data[m-1];
            NRSur7dq2_eval_data_piece(data_piece_eval, q, chiA_coorb, chiB_coorb, data_piece_data);
            gsl_vector_add(h_coorb->modes_real_part[i0+m], data_piece_eval);
            gsl_vector_sub(h_coorb->modes_real_part[i0-m], data_piece_eval);

            // Im[X_plus] gets added to Re[h^{ell, m}] and subtracted from Re[h^{ell, -m}]
            data_piece_data = ell_data->X_imag_plus_data[m-1];
            NRSur7dq2_eval_data_piece(data_piece_eval, q, chiA_coorb, chiB_coorb, data_piece_data);
            gsl_vector_add(h_coorb->modes_imag_part[i0+m], data_piece_eval);
            gsl_vector_sub(h_coorb->modes_imag_part[i0-m], data_piece_eval);

            // Im[X_minus] gets added to both Re[h^{ell, m}] and Re[h^{ell, -m}]
            data_piece_data = ell_data->X_imag_minus_data[m-1];
            NRSur7dq2_eval_data_piece(data_piece_eval, q, chiA_coorb, chiB_coorb, data_piece_data);
            gsl_vector_add(h_coorb->modes_imag_part[i0+m], data_piece_eval);
            gsl_vector_add(h_coorb->modes_imag_part[i0-m], data_piece_eval);
        }
    }

    // Rotate to the inertial frame, write results in h
    MultiModalWaveform_Init(h, NRSUR7DQ2_LMAX, n_coorb);
    TransformModesCoorbitalToInertial(*h, h_coorb, quat_coorb, phi_coorb);

    // Cleanup
    MultiModalWaveform_Destroy(h_coorb);
    for (i=0; i<3; i++) {
        gsl_vector_free(chiA_coorb[i]);
        gsl_vector_free(chiB_coorb[i]);
        gsl_vector_free(quat_coorb[i]);
    }
    gsl_vector_free(quat_coorb[3]);
    gsl_vector_free(phi_coorb);
    gsl_vector_free(data_piece_eval);
}

/**
 * Computes the starting frequency of the NRSur7dq2 waveform approximant when
 * evaluated using fRef=0 (which uses the entire surrogate waveform starting
 * 4500M before the peak amplitude).
 */
double XLALSimInspiralNRSur7dq2StartFrequency(
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 s1x,                      /**< initial value of S1x */
        REAL8 s1y,                      /**< initial value of S1y */
        REAL8 s1z,                      /**< initial value of S1z */
        REAL8 s2x,                      /**< initial value of S2x */
        REAL8 s2y,                      /**< initial value of S2y */
        REAL8 s2z                      /**< initial value of S2z */
) {
#ifdef LAL_PTHREAD_LOCK
    (void) pthread_once(&NRSur7dq2_is_initialized, NRSur7dq2_Init_LALDATA);
#else
    NRSur7dq2_Init_LALDATA();
#endif
    if (m1 < m2) {
        // Now switch the labeling; NRSur7dq2 assumes m1 >= m2.
        REAL8 tmp = m1;
        m1 = m2;
        m2 = tmp;
        tmp = s1x;
        s1x = s2x;
        s2x = tmp;
        tmp = s1y;
        s1y = s2y;
        s2y = tmp;
        tmp = s1z;
        s1z = s2z;
        s2z = tmp;
    }

    double Mtot = (m1 + m2) / LAL_MSUN_SI;
    double Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */
    double q = m1 / m2;
    double y0[11];
    y0[0] = 1.0; // Scalar component of reference quaternion
    int i;
    for (i=1; i<4; i++) y0[i] = 0.0; // Vector components of quaternion
    y0[4] = 0.0; // Reference orbital phase - doesn't affect frequency
    y0[5] = s1x;
    y0[6] = s1y;
    y0[7] = s1z;
    y0[8] = s2x;
    y0[9] = s2y;
    y0[10] = s2z;
    double omega_dimless_start = NRSur7dq2_get_omega(0, q, y0);
    // GW freq is roughly twice the orbital freq
    double f_start_Hz = omega_dimless_start / (LAL_PI * Mtot_sec);
    return f_start_Hz;
}


/**
 * This function evaluates the NRSur7dq2 surrogate model and sums over all ell <= 4 modes to
 * obtain the + and x polarizations.
 * The system is initialized at a time where the orbital angular frequency of the waveform
 * coorbital frame (Eq 10. of https://arxiv.org/abs/1705.07089) is pi * fRef.
 * At this reference point, the system is initialized with the z-axis along the
 * principal eigenvector of the angular momentum operator acting on the waveform (the system is
 * instantaneously aligned with the coprecessing frame) and with an orbital phase of phiRef
 * where this orbital phase is the one in Eq 6 of https://arxiv.org/abs/1705.07089 (...and the
 * coorbital frame).
 * This essentially places the larger black hole on the +ive x-axis, with the smaller black hole
 * on the -ive x-axis, and the initial orbital angular momentum is in the +z direction.
 * Finally, the modes are summed up with the specified inclination angle and an azimuthal angle of 0.
 * When using a custom ModeArray the user must explicitly supply all modes, -ve and +ve m-modes.
 * Note that these are modes in the co-orbital frame, not the inertial frame. */
int XLALSimInspiralNRSur7dq2Polarizations(
        REAL8TimeSeries **hplus,        /**< OUTPUT h_+ vector */
        REAL8TimeSeries **hcross,       /**< OUTPUT h_x vector */
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 inclination,              /**< inclination angle */
        REAL8 deltaT,                   /**< sampling interval (s) */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 distance,                 /**< distance of source (m) */
        REAL8 fMin,                     /**< start GW frequency (Hz) */
        REAL8 fRef,                     /**< reference GW frequency (Hz) */
        REAL8 s1x,                      /**< initial value of S1x */
        REAL8 s1y,                      /**< initial value of S1y */
        REAL8 s1z,                      /**< initial value of S1z */
        REAL8 s2x,                      /**< initial value of S2x */
        REAL8 s2y,                      /**< initial value of S2y */
        REAL8 s2z,                      /**< initial value of S2z */
        LALValue* ModeArray             /**< Container for the ell and m co-orbital modes to generate. To generate all available modes pass NULL */
) {
    int ell, m;
    if ( ModeArray == NULL ) {
        // Use all available modes
        ModeArray = XLALSimInspiralCreateModeArray();
        for (ell=2; ell <= NRSUR7DQ2_LMAX; ell++) {
            XLALSimInspiralModeArrayActivateAllModesAtL(ModeArray, ell);
        }
    } // Otherwise, use the specified modes.

    // BH A is defined to be the one with the larger mass, BH B with the smaller mass.
    if (m1 < m2) {
        // Switch the physical locations of the black holes
        phiRef += LAL_PI;
        // Now switch the labeling
        double tmp = m1;
        m1 = m2;
        m2 = tmp;
        tmp = s1x;
        s1x = s2x;
        s2x = tmp;
        tmp = s1y;
        s1y = s2y;
        s2y = tmp;
        tmp = s1z;
        s1z = s2z;
        s2z = tmp;
    }

    // Parameters
    double massA = m1 / LAL_MSUN_SI;
    double massB = m2 / LAL_MSUN_SI;
    double Mtot = massA+massB;
    double Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */
    double q = massA / massB;
    double chiA0[3], chiB0[3];
    /* The spin components are given in the Lalsimulation source frame
     * (see See Harald Pfeiffer, T18002260-v1 for a diagram). The surrogate frame has the same z
     * but has its x along the line of ascending nodes, so we must rotate the (x, y) spin components
     * by phiRef */
    chiA0[0] = s1x * cos(phiRef) - s1y * sin(phiRef);
    chiA0[1] = s1y * cos(phiRef) + s1x * sin(phiRef);
    chiA0[2] = s1z;
    chiB0[0] = s2x * cos(phiRef) - s2y * sin(phiRef);
    chiB0[1] = s2y * cos(phiRef) + s2x * sin(phiRef);
    chiB0[2] = s2z;
    // Orbital angular frequency = 0.5*(wave angular frequency) = pi*(wave frequency)
    double omegaRef_dimless = (Mtot_sec * fRef) * LAL_PI;
    // Take the inertial frame to be aligned with the coprecessing frame at the reference point:
    double q_ref[4];
    q_ref[0] = 1.0;
    q_ref[1] = 0.0;
    q_ref[2] = 0.0;
    q_ref[3] = 0.0;

    // Evaluate the model modes
    MultiModalWaveform *h_inertial_modes = NULL;
    NRSur7dq2_core(&h_inertial_modes, q, chiA0, chiB0, omegaRef_dimless, phiRef, q_ref, ModeArray);
    if (!h_inertial_modes) {
        return XLAL_FAILURE;
    }

    // Setup the time series for h_+ and h_x evaluated on the surrogate time array
    size_t length = h_inertial_modes->n_times;
    gsl_vector *hplus_model_times = gsl_vector_calloc(length);
    gsl_vector *hcross_model_times = gsl_vector_calloc(length);

    // Sum up the modes
    COMPLEX16 swsh_coef;// = XLALSpinWeightedSphericalHarmonic(spheretheta, spherephi, -2, swsh_L, swsh_m);
    double c_re, c_im;// = creal(swsh_coef);
    double hmre, hmim;
    int i, j;
    bool skip_ell_modes;
    i=0;
    for (ell=2; ell <= NRSUR7DQ2_LMAX; ell++) {
        /* If *any* co-orbital frame mode of this ell is non-zero we need to sum over all modes of this ell,
         * as the frame rotation will have mixed the modes. */
        skip_ell_modes = true;
        for (m=-ell; m<=ell; m++) {
            if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, m) == 1) skip_ell_modes = false;
        }
        if (skip_ell_modes) {
            XLAL_PRINT_INFO("SKIPPING ell=%i modes", ell);
            i += 2*ell + 1;
            continue;
        }
        for (m=-ell; m<=ell; m++) {
            XLAL_PRINT_INFO("Interpolating (%i, %i) mode", ell, m);
            /* See Harald Pfeiffer, T18002260-v1 for frame diagram. The surrogate frame has its z
             * direction aligned with the orbital angular momentum, which agrees with the lowercase
             * source frame z in the diagram. The surrogate x direction, however, points along the
             * line of ascending nodes (the funny omega with circles on the ends). The uppercase Z,
             * which is the direction along which we want to evaluate the waveform, is always in the
             * surrogate frame (y, z) plane. Z is rotated from z towards the +ive y surrogate axis,
             * so we should always evaluate at (inclination, pi/2). */
            swsh_coef = XLALSpinWeightedSphericalHarmonic(inclination, 0.5 * LAL_PI, -2, ell, m);
            c_re = creal(swsh_coef);
            c_im = cimag(swsh_coef);
            for (j=0; j < (int)length; j++) {
                hmre = gsl_vector_get(h_inertial_modes->modes_real_part[i], j);
                hmim = gsl_vector_get(h_inertial_modes->modes_imag_part[i], j);
                hplus_model_times->data[j] += (c_re * hmre - c_im * hmim);
                hcross_model_times->data[j] -= (c_im * hmre + c_re * hmim); // - sign since h = h_+ - i*h_x
            }
            i += 1;
        }
    }

    // Scale the amplitude based on the distance
    double a0 = Mtot_sec * LAL_C_SI / distance;
    gsl_vector_scale(hplus_model_times, a0);
    gsl_vector_scale(hcross_model_times, a0);

    // Setup output times
    gsl_vector *model_times = (&__lalsim_NRSur7dq2_data)->t_coorb;
    double deltaT_dimless = deltaT / Mtot_sec;
    double t0 = gsl_vector_get(model_times, 0);
    double start_freq = XLALSimInspiralNRSur7dq2StartFrequency(m1, m2, s1x, s2y, s1z, s2x, s2y, s2z);
    if (fMin >= start_freq) {
        double omegaMin_dimless = (Mtot_sec * fMin) * LAL_PI;
        t0 = NRSur7dq2_get_t_ref(omegaMin_dimless, q, chiA0, chiB0, q_ref, phiRef);
    } else if (fMin > 0) {
        XLAL_ERROR_REAL8(XLAL_FAILURE, "fMin should be 0 or >= %0.8f for this configuration, got %0.8f", start_freq, fMin);
    }
    double tf = gsl_vector_get(model_times, length-1);
    int nt = (int) ceil((tf - t0) / deltaT_dimless);
    gsl_vector *output_times = gsl_vector_alloc(nt);
    for (j=0; j<nt; j++) {
        gsl_vector_set(output_times, j, t0 + deltaT_dimless*j);
    }

    // Interpolate onto output times
    double t;
    LIGOTimeGPS epoch = LIGOTIMEGPSZERO;
    XLALGPSAdd( &epoch, Mtot_sec * t0);
    *hplus = XLALCreateREAL8TimeSeries("hp: TD waveform", &epoch, 0.0, deltaT, &lalStrainUnit, nt);
    *hcross = XLALCreateREAL8TimeSeries("hc: TD waveform", &epoch, 0.0, deltaT, &lalStrainUnit, nt);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spl_hplus = gsl_spline_alloc(gsl_interp_cspline, length);
    gsl_spline *spl_hcross = gsl_spline_alloc(gsl_interp_cspline, length);
    gsl_spline_init(spl_hplus, model_times->data, hplus_model_times->data, length);
    gsl_spline_init(spl_hcross, model_times->data, hcross_model_times->data, length);
    for (j=0; j<nt; j++) {
        t = gsl_vector_get(output_times, j);
        (*hplus)->data->data[j] = gsl_spline_eval(spl_hplus, t, acc);
        (*hcross)->data->data[j] = gsl_spline_eval(spl_hcross, t, acc);
    }

    // Cleanup
    gsl_vector_free(hplus_model_times);
    gsl_vector_free(hcross_model_times);
    gsl_vector_free(output_times);
    gsl_interp_accel_free(acc);
    gsl_spline_free(spl_hplus);
    gsl_spline_free(spl_hcross);
    MultiModalWaveform_Destroy(h_inertial_modes);

    return XLAL_SUCCESS;
}

/**
 * This function evaluates the NRSur7dq2 surrogate model and returns the inertial frame modes
 * in the form of a SphHarmTimeSeries. The system is initialized at a time where the orbital
 * angular frequency of the waveform coorbital frame (Eq 10. of https://arxiv.org/abs/1705.07089)
 * is pi * fRef. At this reference point, the system is initialized with the z-axis along the
 * principal eigenvector of the angular momentum operator acting on the waveform (the system is
 * instantaneously aligned with the coprecessing frame) and with an orbital phase of phiRef
 * where this orbital phase is the one in Eq 6 of https://arxiv.org/abs/1705.07089 (...and the
 * coorbital frame).
 * This essentially places the larger black hole on the +ive x-axis, with the smaller black hole
 * on the -ive x-axis, and the initial orbital angular momentum is in the +z direction.
 * When using a custom ModeArray the user must explicitly supply all modes, -ve and +ve m-modes.
 * Note that these are modes in the co-orbital frame, not the inertial frame. */

SphHarmTimeSeries *XLALSimInspiralNRSur7dq2Modes(
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 deltaT,                   /**< sampling interval (s) */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 S1x,                      /**< x-component of the dimensionless spin of object 1 */
        REAL8 S1y,                      /**< y-component of the dimensionless spin of object 1 */
        REAL8 S1z,                      /**< z-component of the dimensionless spin of object 1 */
        REAL8 S2x,                      /**< x-component of the dimensionless spin of object 2 */
        REAL8 S2y,                      /**< y-component of the dimensionless spin of object 2 */
        REAL8 S2z,                      /**< z-component of the dimensionless spin of object 2 */
        REAL8 fMin,                     /**< start GW frequency (Hz) */
        REAL8 fRef,                     /**< reference GW frequency (Hz) */
        REAL8 distance,                 /**< distance of source (m) */
        LALValue* ModeArray             /**< Container for the ell and m modes to generate. To generate all available modes pass NULL */
) {
    SphHarmTimeSeries *hlms = NULL;

    int ell, m;
    if ( ModeArray == NULL ) {
        // Use all available modes
        ModeArray = XLALSimInspiralCreateModeArray();
        for (ell=2; ell <= NRSUR7DQ2_LMAX; ell++) {
            XLALSimInspiralModeArrayActivateAllModesAtL(ModeArray, ell);
        }
    } // Otherwise, use the specified modes.

    double chiA0[3], chiB0[3];
    chiA0[0] = S1x;
    chiA0[1] = S1y;
    chiA0[2] = S1z;
    chiB0[0] = S2x;
    chiB0[1] = S2y;
    chiB0[2] = S2z;

    // In NRSur7dq2, BH A is defined to be the one with the larger mass, BH B with the smaller mass.
    // TODO: Do we need an appropriate rotation by pi after evaluating to agree with LAL conventions?
    int i;
    if (m2 > m1) {
        double tmp = m1;
        m1 = m2;
        m2 = tmp;
        phiRef += LAL_PI;
        for (i=0; i<3; i++) {
            tmp = chiA0[i];
            chiA0[i] = chiB0[i];
            chiB0[i] = tmp;
        }
    }

    // Parameters
    double massA = m1 / LAL_MSUN_SI;
    double massB = m2 / LAL_MSUN_SI;
    double Mtot = massA+massB;
    double Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */
    double q = massA / massB;
    double omegaRef_dimless = (Mtot_sec * fRef) * LAL_PI; // Extra factor of 0.5 to get orbital frequency
    double q_ref[4]; // For now, just use identity; coprecessing frame is aligned at the reference point.
    q_ref[0] = 1.0;
    q_ref[1] = 0.0;
    q_ref[2] = 0.0;
    q_ref[3] = 0.0;

    // Evaluate the model modes
    MultiModalWaveform *h_inertial = NULL;
    NRSur7dq2_core(&h_inertial, q, chiA0, chiB0, omegaRef_dimless, phiRef, q_ref, ModeArray);
    if (!h_inertial) {
        XLAL_PRINT_INFO("NRSur7dq2_core failed!");
        return hlms;
    }

    // Scale the amplitude based on the distance
    double a0 = Mtot_sec * LAL_C_SI / distance;

    // Setup dimensionless output times
    gsl_vector *model_times = (&__lalsim_NRSur7dq2_data)->t_coorb;
    double deltaT_dimless = deltaT / Mtot_sec;
    double t0 = gsl_vector_get(model_times, 0);
    double start_freq = XLALSimInspiralNRSur7dq2StartFrequency(m1, m2, S1x, S2y, S1z, S2x, S2y, S2z);
    if (fMin >= start_freq) {
        double omegaMin_dimless = (Mtot_sec * fMin) * LAL_PI;
        t0 = NRSur7dq2_get_t_ref(omegaMin_dimless, q, chiA0, chiB0, q_ref, phiRef);
    } else if (fMin > 0) {
        XLAL_ERROR_NULL(XLAL_FAILURE, "fMin should be 0 or >= %0.8f for this configuration, got %0.8f", start_freq, fMin);
    }
    size_t length = model_times->size;
    double tf = gsl_vector_get(model_times, length-1);
    int nt = (int) ceil((tf - t0) / deltaT_dimless);
    gsl_vector *output_times = gsl_vector_alloc(nt);
    int j;
    for (j=0; j<nt; j++) {
        gsl_vector_set(output_times, j, t0 + deltaT_dimless*j);
    }

    // Setup interpolation onto dimensionless output times
    double t;
    LIGOTimeGPS epoch = LIGOTIMEGPSZERO;
    XLALGPSAdd( &epoch, Mtot_sec * t0);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    // Create LAL modes
    COMPLEX16TimeSeries *tmp_mode;
    i=0;
    bool skip_ell_modes;
    char mode_name[32];
    for (ell=2; ell <= NRSUR7DQ2_LMAX; ell++) {
        /* If *any* co-orbital frame mode of this ell is non-zero we need to sum over all modes of this ell,
         * as the frame rotation will have mixed the modes. */
        skip_ell_modes = true;
        for (m=-ell; m<=ell; m++) {
            if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, m) == 1) skip_ell_modes = false;
        }
        if (skip_ell_modes) {
            XLAL_PRINT_INFO("SKIPPING ell=%i modes", ell);
            i += 2*ell + 1;
            continue;
        }
        for (m=-ell; m<=ell; m++) {
            XLAL_PRINT_INFO("Interpolating (%i, %i) mode", ell, m);

            gsl_vector_scale(h_inertial->modes_real_part[i], a0);
            gsl_vector_scale(h_inertial->modes_imag_part[i], a0);
            snprintf(mode_name, sizeof(mode_name), "(%d, %d) mode", ell, m);
            tmp_mode = XLALCreateCOMPLEX16TimeSeries(mode_name, &epoch, 0.0,
                    deltaT, &lalStrainUnit, nt);
            gsl_spline *spl_re = gsl_spline_alloc(gsl_interp_cspline, length);
            gsl_spline *spl_im = gsl_spline_alloc(gsl_interp_cspline, length);
            gsl_spline_init(spl_re, model_times->data, h_inertial->modes_real_part[i]->data, length);
            gsl_spline_init(spl_im, model_times->data, h_inertial->modes_imag_part[i]->data, length);
            for (j=0; j<nt; j++) {
                t = gsl_vector_get(output_times, j);
                tmp_mode->data->data[j] = gsl_spline_eval(spl_re, t, acc);
                tmp_mode->data->data[j] += I * gsl_spline_eval(spl_im, t, acc);
            }
            gsl_spline_free(spl_re);
            gsl_spline_free(spl_im);
            hlms = XLALSphHarmTimeSeriesAddMode(hlms, tmp_mode, ell, m);
            XLALDestroyCOMPLEX16TimeSeries(tmp_mode);
            i += 1;
        }
    }

    // Cleanup
    MultiModalWaveform_Destroy(h_inertial);
    gsl_vector_free(output_times);
    gsl_interp_accel_free(acc);

    return hlms;
}

/** @} */
/** @} */
