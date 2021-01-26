/*  Copyright (C) 2019 Vijay Varma
 *  Utils for evaluates NR surrogate remnant fits for the mass, spin and recoil
 *  kick of the final BH left behind in a BBH merger.
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
 * \brief Utils for NR surrogates for remnant BH mass, spin and recoil kick.
 */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include "LALSimNRSurRemnantUtils.h"

#include <libgen.h>

#include <assert.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex_math.h>

#include <lal/SeqFactories.h>
#include <lal/FileIO.h>
#include <lal/H5FileIO.h>

#include <lal/XLALError.h>
#include <lal/LALConstants.h>
#include <lal/LALSimIMR.h>

#include "LALSimBlackHoleRingdown.h"
#include "LALSimIMRSEOBNRROMUtilities.c"

//*************************************************************************/
//************************* function definitions **************************/
//*************************************************************************/

//********************* Surrogate loading functions ************************/
// These should only be called once, when initializing the surrogate

/**
 * Loads H5 file for a NRSurRemnant model.
 */
void NRSurRemnant_LoadH5File(
    LALH5File **file,   /**< Output: Returns the opened H5 file. */
    const char* NRSurRemnant_DATAFILE   /**< H5 file name to load. */
) {

    char *path = XLALFileResolvePathLong(NRSurRemnant_DATAFILE, PKG_DATA_DIR);
    if (path==NULL) {
        XLAL_ERROR_VOID(XLAL_ENOENT,
            "Unable to find data file %s in $LAL_DATA_PATH\n",
            NRSurRemnant_DATAFILE);
    }

    char *dir = dirname(path);
    const UINT4 size = strlen(dir) + strlen(NRSurRemnant_DATAFILE) + 2;
    char *file_path = XLALMalloc(size);
    snprintf(file_path, size, "%s/%s", dir, NRSurRemnant_DATAFILE);

    *file = XLALH5FileOpen(file_path, "r");
    if (*file==NULL) {
        XLAL_ERROR_VOID(XLAL_EIO,
            "Unable to load data file %s in $LAL_DATA_PATH."
            " File may be corrupted.\n", NRSurRemnant_DATAFILE);
    }
    XLALFree(path);
    XLALFree(file_path);
}


/**
 * Loads a single NRSurRemnant GPR fit, as described in the supplementary
 * materials of arxiv:1809.09125.
 */
int NRSurRemnant_LoadScalarFit(
    ScalarFitData **fit_data, /**< Output: Fit data. *fit_data
                                should be NULL. Space will be allocated. */
    LALH5File *file,          /**< Opened H5 file. */
    const char *grp_name      /**< H5 group name. */
) {

    if (fit_data == NULL || *fit_data != NULL) {
        XLAL_ERROR(XLAL_EFAULT, "*fit_data should be NULL");
    }
    if (file == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "file should not be NULL");
    }

    // Open h5 group
    LALH5File *sub;
    sub = XLALH5GroupOpen(file, grp_name);
    *fit_data = XLALMalloc(sizeof(ScalarFitData));

    // Load fit data
    GPRHyperParams *hyperparams = XLALMalloc(sizeof(GPRHyperParams));

    //// Load scalars needed for fit

    // constant_value = \sigma_k^2 in Eq. S3 of arxiv:1809.09125
    int ret = ReadHDF5DoubleDataset(&(hyperparams->constant_value),
            sub, "constant_value");
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load constant_value.");
    }

    // Mean value before doing GPR fit, usually zero. We already remove the
    // mean and store it in data_mean.
    ret = ReadHDF5DoubleDataset(&(hyperparams->y_train_mean),
            sub, "y_train_mean");
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load y_train_mean.");
    }

    // Mean of fit data.
    ret = ReadHDF5DoubleDataset(&((*fit_data)->data_mean),
            sub, "data_mean");
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load data_mean.");
    }

    // Standard deviation of fit data.
    ret = ReadHDF5DoubleDataset(&((*fit_data)->data_std),
            sub, "data_std");
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load data_std.");
    }

    // Intercept of linear fit.
    ret = ReadHDF5DoubleDataset(&((*fit_data)->lin_intercept),
            sub, "lin_intercept");
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load lin_intercept.");
    }


    //// Load arrays needed for fit

    // Length scale parameters in the kernel.
    // Array of size=dom, where dim=dimensions of the model.
    hyperparams->length_scale = NULL;
    ret = ReadHDF5RealVectorDataset(sub, "length_scale",
            &(hyperparams->length_scale));
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load length_scale.");
    }

    // alpha = \f$ K_{x x}^{-1} {\bf f}\f$, which is used in the mean value
    // in Eq. S2 of arxiv:1809.09125. This can be precomputed as it does not
    // depend on the x_{*} values.
    // Array of size=N, where N=size of training dataset.
    hyperparams->alpha = NULL;
    ret = ReadHDF5RealVectorDataset(sub, "alpha",
            &(hyperparams->alpha));
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load alpha.");
    }

    // Coefs of of linear fit.
    (*fit_data)->lin_coef = NULL;
    ret = ReadHDF5RealVectorDataset(sub, "lin_coef", &((*fit_data)->lin_coef));
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load lin_coef.");
    }

    (*fit_data)->hyperparams = hyperparams;


    XLALH5FileClose(sub);
    return ret;
}


/**
 * Loads a vector of NRSurRemnant GPR fits
 */
int NRSurRemnant_LoadVectorFit(
    VectorFitData **vector_data, /**< Output: Vector of fit data. *vector_data
                                 should be NULL. Space will be allocated. */
    UINT4 vec_dim,               /**< Length of the vector */
    LALH5File *file,             /**< Opened H5 file. */
    const char *grp_name         /**< H5 group name. */
) {

    if (vector_data == NULL || *vector_data != NULL) {
        XLAL_ERROR(XLAL_EFAULT, "*vector_data should be NULL");
    }

    *vector_data = XLALMalloc(sizeof(VectorFitData));
    (*vector_data)->fit_data = XLALMalloc( vec_dim * sizeof(ScalarFitData *) );

    const size_t str_size = 20;
    char *sub_grp_name = XLALMalloc(str_size);
    UNUSED size_t nwritten;

    // For each component load a Scalar fit
    int ret = XLAL_SUCCESS;
    for (UINT4 i=0; i<vec_dim; i++){

        nwritten = snprintf(sub_grp_name, str_size, "%s/comp_%u", grp_name, i);
        assert(nwritten < str_size);

        ScalarFitData *fit_data = NULL;
        ret = NRSurRemnant_LoadScalarFit(&fit_data, file, sub_grp_name);
        (*vector_data)->fit_data[i] = fit_data;
    }
    (*vector_data)->vec_dim = vec_dim;

    return ret;
}

/**
 * Initializes fit data for a precessing NRSurRemnant.
 *
 * The data includes a ScalarFitData for the final mass, and a VectorFitData
 * for the final spin and kick 3-vectors.
 */
int PrecessingNRSurRemnant_Init(
    PrecessingRemnantFitData *sur_data, /**< Output: Loaded surrogate data. */
    LALH5File *file                     /**< Opened H5 file. */
) {

    if (sur_data == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "sur_data should not be NULL");
    }
    if (file == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "file should not be NULL");
    }

    if (sur_data->setup) {
        XLAL_ERROR(XLAL_FAILURE,
            "Model was already initialized. Exiting.");
    }

    // x_train contains the training parameters and are common for all fits
    gsl_matrix *x_train = NULL;
    int ret = ReadHDF5RealMatrixDataset(file, "GPR_X_train", &x_train);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load GPR_X_train.");
    }
    sur_data->x_train = x_train;

    // dimension of the surrogate
    sur_data->params_dim = x_train->size2;

    // final mass
    sur_data->mf_data = NULL;
    ret = NRSurRemnant_LoadScalarFit(&(sur_data->mf_data), file, "mf");

    // final spin vector
    sur_data->chif_data = NULL;
    NRSurRemnant_LoadVectorFit(&(sur_data->chif_data), 3, file, "chif");

    // recoil kick vector
    sur_data->vf_data = NULL;
    NRSurRemnant_LoadVectorFit(&(sur_data->vf_data), 3, file, "vf");

    if (ret == XLAL_SUCCESS){
        sur_data->setup = 1;
    }

    return ret;
}


/**
 * Initializes fit data for an aligned-spin NRSurRemnant.
 *
 * The data includes a ScalarFitData for the final mass, z-component of the
 * final spin, and x and y components of the recoil kick. The other components
 * of the spin and kick are zero due to the symmetries of aligned-spin systems
 * and are not modelled by the surrogate.
 */
int AlignedSpinNRSurRemnant_Init(
    AlignedSpinRemnantFitData *sur_data, /**< Output: Loaded surrogate data. */
    LALH5File *file                      /**< Opened H5 file. */
) {

    if (sur_data == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "sur_data should not be NULL");
    }
    if (file == NULL) {
        XLAL_ERROR(XLAL_EFAULT, "file should not be NULL");
    }

    if (sur_data->setup) {
        XLAL_ERROR(XLAL_FAILURE,
            "Model was already initialized. Exiting.");
    }

    // x_train contains the training parameters and are common for all fits
    gsl_matrix *x_train = NULL;
    int ret = ReadHDF5RealMatrixDataset(file, "GPR_X_train", &x_train);
    if (ret != XLAL_SUCCESS) {
        XLAL_ERROR(XLAL_EFUNC, "Failed to load GPR_X_train.");
    }
    sur_data->x_train = x_train;

    // dimension of the surrogate
    sur_data->params_dim = x_train->size2;

    // final mass
    sur_data->mf_data = NULL;
    ret = NRSurRemnant_LoadScalarFit(&(sur_data->mf_data), file, "mf");

    // z-component of final spin (other components are zero)
    sur_data->chifz_data = NULL;
    ret = NRSurRemnant_LoadScalarFit(&(sur_data->chifz_data), file, "chifz");

    // x-component of recoil kick
    sur_data->vfx_data = NULL;
    ret = NRSurRemnant_LoadScalarFit(&(sur_data->vfx_data), file, "vfx");

    // y-component of recoil kick, z-component is zero
    sur_data->vfy_data = NULL;
    ret = NRSurRemnant_LoadScalarFit(&(sur_data->vfy_data), file, "vfy");

    if (ret == XLAL_SUCCESS){
        sur_data->setup = 1;
    }

    return ret;
}
