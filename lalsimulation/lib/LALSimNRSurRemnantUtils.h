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

#include "LALSimNRHybSurUtilities.h"




/*
 * Custom Warning function because XLAL_PRINT_WARNING is not enabled by
 * default.
 */
#define print_warning(...) \
    if (lalDebugLevel & LALERRORBIT) \
        printf("Warning - %s (%s:%d): ", __func__, __FILE__, __LINE__); \
        printf(__VA_ARGS__);

//*************************************************************************/
//************************* struct definitions ****************************/
//*************************************************************************/

// These GPR fits need the same struct as NRHybSurFitData, but renaming it
// here to avoid confusion
#define ScalarFitData NRHybSurFitData

/**
 * Vector of fit data, where each component is a ScalarFitData
 */
typedef struct tagVectorFitData {
    UINT4 vec_dim;              /**< Length of the vector. */
    ScalarFitData **fit_data;   /**< One ScalarFitData per component */
} VectorFitData;

/**
 * NRSurRemnant GPR fit data for the mass, spin, and recoil kick for
 * generically precessing BBHs.
 *
 * The final mass is included as a sclar, while the final spin and kick are
 * included as 3-vectors.
 */
typedef struct tagPrecessingRemnantFitData {
    UINT4 setup;         /**< Indicates if NRSurRemnant has been initialized */
    UINT4 params_dim;    /**< Dimensions of the model. */
    ScalarFitData *mf_data;     /**< Fit data for final mass. */
    VectorFitData *chif_data;   /**< Fit data for final spin. */
    VectorFitData *vf_data;     /**< Fit data for recoil kick. */
    gsl_matrix *x_train; /**< Training set parameters, needed for GPR fits. */
} PrecessingRemnantFitData;

/**
 * NRSurRemnant GPR fit data for the mass, spin, and recoil kick for
 * aligned-spin BBHs.
 *
 * The final mass, z-component of the final spin, and x and y components of the
 * recoil kick are included as scalars. The other components of the spin and
 * kick are zero due to the symmetries of aligned-spin systems and are not
 * modelled by the surrogate.
 */
typedef struct tagAlignedSpinRemnantFitData {
    UINT4 setup;         /**< Indicates if NRSurRemnant has been initialized */
    UINT4 params_dim;    /**< Dimensions of the model. */
    ScalarFitData *mf_data;     /**< Fit data for final mass. */
    ScalarFitData *chifz_data;     /**< Fit data for final mass. */
    ScalarFitData *vfx_data;     /**< Fit data for final mass. */
    ScalarFitData *vfy_data;     /**< Fit data for final mass. */
    gsl_matrix *x_train; /**< Training set parameters, needed for GPR fits. */
} AlignedSpinRemnantFitData;

//*************************************************************************/
//************************* function declarations *************************/
//*************************************************************************/

void NRSurRemnant_LoadH5File(
    LALH5File **file,
    const char* NRSurRemnant_DATAFILE
);

int NRSurRemnant_LoadScalarFit(
    ScalarFitData **fit_data,
    LALH5File *file,
    const char *grp_name
);

int NRSurRemnant_LoadVectorFit(
    VectorFitData **vector_data,
    UINT4 vec_dim,
    LALH5File *file,
    const char *grp_name
);

int PrecessingNRSurRemnant_Init(
    PrecessingRemnantFitData *sur_data,
    LALH5File *file
);

int AlignedSpinNRSurRemnant_Init(
    AlignedSpinRemnantFitData *sur_data,
    LALH5File *file
);
