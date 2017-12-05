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

#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif

// These are to rescale the mass ratio fit range from [0.99, 2.01] to [-1, 1].
static const double NRSUR7DQ2_Q_FIT_OFFSET = -2.9411764705882355; // = (0.5*(2.01 - 0.99)) * (2 / (2.01 - 0.99))
static const double NRSUR7DQ2_Q_FIT_SLOPE = 1.9607843137254901; // = 2 / (2.01 - 0.99)

// Allow a tiny bit of extrapolation in mass ratio and spin
static const double NRSUR7DQ2_Q_MIN = 0.999;
static const double NRSUR7DQ2_Q_MAX = 2.001;
static const double NRSUR7DQ2_CHI_MAX = 0.801;

static const int NRSUR7DQ2_LMAX = 4;

// Surrogate model data, in LAL_DATA_PATH. File available in lalsuite-extra or at
// https://www.black-holes.org/surrogates
static const char NRSUR7DQ2_DATAFILE[] = "NRSur7dq2.h5";

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t NRSur7dq2_is_initialized = PTHREAD_ONCE_INIT;
#endif

/***********************************************************************************/
/****************************** Type definitions ***********************************/
/***********************************************************************************/

// Data used in a single fit
typedef struct tagFitData {
    gsl_matrix_long *basisFunctionOrders; // matrix of (n_coefs x 7) basis function orders
    gsl_vector *coefs;
    int n_coefs;
} FitData;

// Data used in vector fit
typedef struct tagVectorFitData {
    gsl_matrix_long *basisFunctionOrders;
    gsl_vector *coefs;
    gsl_vector_long *componentIndices;
    int n_coefs;
    int vec_dim;
} VectorFitData;

typedef struct tagDynamicsNodeFitData {
    FitData *omega_data;
    VectorFitData *omega_copr_data;
    VectorFitData *chiA_dot_data;
    VectorFitData *chiB_dot_data;
} DynamicsNodeFitData;

typedef struct tagWaveformDataPiece {
    int n_nodes;
    FitData **fit_data;
    gsl_matrix *empirical_interpolant_basis;
    gsl_vector_long *empirical_node_indices;
} WaveformDataPiece;

typedef struct tagWaveformFixedEllModeData {
    int ell;
    WaveformDataPiece *m0_real_data;
    WaveformDataPiece *m0_imag_data;
    // X_\pm^{l, m} = \frac{1}{2}( h^{l, m} \pm h^{l, -m} )
    WaveformDataPiece **X_real_plus_data; // One for each 1 <= m <= ell
    WaveformDataPiece **X_real_minus_data;
    WaveformDataPiece **X_imag_plus_data; // One for each 1 <= m <= ell
    WaveformDataPiece **X_imag_minus_data;
} WaveformFixedEllModeData;

typedef struct tagNRSur7dq2Data {
    UINT4 setup;
    int LMax;
    gsl_vector *t_ds;
    gsl_vector *t_ds_half_times; // t_1/2, t_3/2, and t_5/2 used to start up integration
    gsl_vector *t_coorb;
    DynamicsNodeFitData **ds_node_data;
    DynamicsNodeFitData **ds_half_node_data; // 3 half nodes at t_ds_half_times
    WaveformFixedEllModeData **coorbital_mode_data; // One for each 2 <= ell <= LMax
} NRSur7dq2Data;

typedef struct tagMultiModalWaveform {
    int n_modes;
    int n_times;
    int *lvals;
    int *mvals;
    gsl_vector **modes_real_part;
    gsl_vector **modes_imag_part;
} MultiModalWaveform;

typedef struct tagWignerDMatrices {
    int LMax;
    int n_entries;
    gsl_vector **real_part;
    gsl_vector **imag_part;
} WignerDMatrices;

typedef struct tagRealPowers {
    int LMax;
    int n_entries;
    int n_times;
    gsl_vector **powers;
} RealPowers;

typedef struct tagComplexPowers {
    int LMax;
    int n_entries;
    int n_times;
    gsl_vector **real_part;
    gsl_vector **imag_part;
} ComplexPowers;

/**** Global surrogate data ****/
static NRSur7dq2Data __lalsim_NRSur7dq2_data;


/***********************************************************************************/
/****************************** Function declarations*******************************/
/***********************************************************************************/
static void NRSur7dq2_MultiModalWaveform_Init(MultiModalWaveform **wave, int LMax, int n_times);
static void NRSur7dq2_MultiModalWaveform_Destroy(MultiModalWaveform *wave);
static void NRSur7dq2_WignerDMatrices_Init(WignerDMatrices **matrices, int n_times, int LMax);
static void NRSur7dq2_WignerDMatrices_Destroy(WignerDMatrices *matrices);
static int NRSur7dq2_WignerDMatrix_Index(int ell, int m, int mp);
static void NRSur7dq2_WignerDMatrices_Compute(WignerDMatrices *matrices, gsl_vector **quat);
static void NRSur7dq2_ComplexPowers_Init(ComplexPowers **cp, int LMax, int n_times);
static void NRSur7dq2_ComplexPowers_Destroy(ComplexPowers *cp);
static void NRSur7dq2_RealPowers_Init(RealPowers **rp, int LMax, int n_times);
static void NRSur7dq2_RealPowers_Destroy(RealPowers *rp);
static void NRSur7dq2_ComplexPowers_Compute(ComplexPowers *cp, gsl_vector *x, gsl_vector *y);
static void NRSur7dq2_RealPowers_Compute(RealPowers *rp, gsl_vector *x);
static void NRSur7dq2_Init_LALDATA(void);
static int NRSur7dq2_Init(NRSur7dq2Data *data, LALH5File *file);
static void NRSur7dq2_LoadDynamicsNode(DynamicsNodeFitData **ds_node_data, LALH5File *sub, int i);
static void NRSur7dq2_LoadCoorbitalEllModes(WaveformFixedEllModeData **coorbital_mode_data, LALH5File *file, int i);
static void NRSur7dq2_LoadWaveformDataPiece(LALH5File *sub, WaveformDataPiece **data, bool invert_sign);
static bool NRSur7dq2_IsSetup(void);
static double ipow(double base, int exponent); // integer powers

static void complex_vector_mult(gsl_vector *x1, gsl_vector *y1, gsl_vector *x2, gsl_vector *y2, gsl_vector *tmp1, gsl_vector *tmp2);
static double NRSur7dq2_eval_fit(FitData *data, double *x);

static void NRSur7dq2_eval_vector_fit(
    double *res, // Result
    VectorFitData *data, // Data for fit
    double *x // size 7, giving mass ratio q, and dimensionless spin components
);

static void NRSur7dq2_normalize_y(
    double chiANorm,
    double chiBNorm,
    double *y
); // Helper to keep spin magnitude constant

static void NRSur7dq2_normalize_results(
    double normA,
    double normB,
    gsl_vector **quat,
    gsl_vector **chiA,
    gsl_vector **chiB
);

static void NRSur7dq2_rotate_spins(gsl_vector **chiA, gsl_vector **chiB, gsl_vector *phi, bool backwards);

static void NRSur7dq2_ds_fit_x(
    double *x, //Result, length 7
    double q, // Mass ratio
    double *y // [q0, qx, qy, qz, phi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
);

static void NRSur7dq2_assemble_dydt(
    double *dydt,           // Result, length 11
    double *y,              // ODE solution at the current time step
    double *Omega_coorb_xy, // a form of time derivative of the coprecessing frame
    double omega,           // orbital angular frequency in the coprecessing frame
    double *chiA_dot,       // chiA time derivative
    double *chiB_dot        // chiB time derivative
);

static void NRSur7dq2_ab4_dy(
    double *dy,     // Result, length 11
    double *k1,     // dy/dt evaluated at node 1
    double *k2,     // dy/dt evaluated at node 2
    double *k3,     // dy/dt evaluated at node 3
    double *k4,     // dy/dt evaluated at node 4
    double dt1,     // t_2 - t_1
    double dt2,     // t_3 - t_2
    double dt3,     // t_4 - t_3
    double dt4      // t_5 - t_4
);

static double factorial(int n);
static double factorial_ratio(int n, int k);
static double binomial(int n, int k);
static double wigner_coef(int ell, int mp, int m);

static double cubic_interp(double xout, double *x, double *y);
static gsl_vector *spline_array_interp(gsl_vector *xout, gsl_vector *x, gsl_vector *y);

static double NRSur7dq2_get_omega(size_t node_index, double q, double *y0);
static double NRSur7dq2_get_t_ref(double omega_ref, double q, double *chiA0, double *chiB0, double *q_ref, double phi_ref);

static void NRSur7dq2_get_time_deriv_from_index(
    double *dydt,       // Output: dy/dt evaluated at the ODE time node with index i0. Must have space for 11 entries.
    int i0,             // Time node index. i0=-1, -2, and -3 are used for time nodes 1/2, 3/2, and 5/2 respectively.
    double q,           // Mass ratio
    double *y           // Current ODE state: [q0, qx, qy, qz, orbphase, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
);

static void NRSur7dq2_get_time_deriv(
    double *dtdy,       // Output: dy/dt evaluated at time t. Must have space for 11 entries.
    double t,           // Time at which the ODE should be evaluated.
    double q,           // Mass ratio
    double *y           // Current ODE state
);

static int NRSur7dq2_initialize_at_dynamics_node(
    double *dynamics_data,  // ODE output
    double t_ref,           // reference time. t_ds[i0] will be close to t_ref.
    double q,               // mass ratio
    double *chiA0,          // chiA at t_ref.
    double *chiB0,          // chiB at t_ref.
    double phi_ref,         // orbital phase at t_ref.
    double *q_ref,          // quaternion at t_ref.
    double normA,           // |chiA|
    double normB            // |chiB|
);

static void NRSur7dq2_initialize_RK4_with_half_nodes(
    double *dynamics_data,  // ODE output
    double *time_steps,     // Output: first three time steps. Should have size 3.
    double *dydt0,             // Output: dydt at node 0. Should have size 11.
    double *dydt1,             // Output: dydt at node 1. Should have size 11.
    double *dydt2,             // Output: dydt at node 2. Should have size 11.
    double *dydt3,             // Output: dydt at node 3. Should have size 11.
    double normA,           // |chiA|
    double normB,           // |chiB|
    double q                // mass ratio
);

static int NRSur7dq2_initialize_RK4(
    double *dynamics_data,  // ODE output
    double *time_steps,     // Output: first three time steps. Should have size 3.
    double *dydt0,             // Output: dydt at node i0 + 0. Should have size 11.
    double *dydt1,             // Output: dydt at node i0 + 1. Should have size 11.
    double *dydt2,             // Output: dydt at node i0 + 2. Should have size 11.
    double *dydt3,             // Output: dydt at node i0 + 3. Should have size 11.
    double normA,           // |chiA|
    double normB,           // |chiB|
    double q,               // mass ratio
    int i0                  // the node that is already initialized
);

static void NRSur7dq2_integrate_AB4(
    double *dynamics_data,  // ODE output
    double *time_steps,     // The first three time steps beginning at i_start.
    double *dydt0,          // dydt at node i_start
    double *dydt1,          // dydt at node i_start + 1
    double *dydt2,          // dydt at node i_start + 2
    double *dydt3,          // dydt at node i_start + 3
    double normA,           // |chiA|
    double normB,           // |chiB|
    double q,               // mass ratio
    int i_start             // nodes i_start through i_start+3 are already initialized
);

static int NRSur7dq2_IntegrateDynamics(
    double *dynamics_data,  // Output: length (n * 11), where entries 11*m <= i < 11*(m+1) are
                            // [q0, qx, qy, qz, varphi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
    double q,               // Mass ratio mA / mB
    double *chiA0,          // chiA at the reference point
    double *chiB0,          // chiB at the reference point
    double omega_ref,       // orbital angular frequency at the reference point
    double phi_ref,         // orbital phase at the reference point
    double *q_ref           // coprecessing quaterion at the reference point
);

static void NRSur7dq2_eval_data_piece(
    gsl_vector *result, // Output: Should have already been assigned space
    double q,           // Mass ratio
    gsl_vector **chiA,  // 3 gsl_vector *s, one for each (coorbital) component
    gsl_vector **chiB,  // similar to chiA
    WaveformDataPiece *data // The data piece to evaluate
);

static void NRSur7dq2_TransformModes(
    MultiModalWaveform *h,          // Output. Dimensionless waveform modes sampled on t_coorb. Should be initialized.
    MultiModalWaveform *h_coorb,    // Coorbital frame waveform modes.
    gsl_vector **quat,              // should be *quat[4], one for each component. Coprecessing frame quaternions.
    gsl_vector *orbphase            // Orbital phase, used to obtain the coprecessing modes.
);

static void NRSur7dq2_core(
    MultiModalWaveform **h, // Output. Dimensionless waveform modes sampled on t_coorb
    double q,               // Mass ratio mA / mB
    double *chiA0,          // chiA at the reference point
    double *chiB0,          // chiB at the reference point
    double omega_ref,       // orbital angular frequency at the reference point
    double phi_ref,         // orbital phase at the reference point
    double *q_ref,          // coprecessing quaterion at the reference point
    int LMax                // Maximum ell mode to evaluate (NRSur7dq2 contains ell=2, 3, 4 modes)
);

/***********************************************************************************/
/****************************** Function Definitions *******************************/
/***********************************************************************************/
static void NRSur7dq2_MultiModalWaveform_Init(
    MultiModalWaveform **wave,
    int LMax,
    int n_times
) {
    if (!wave) exit(1);
    if (LMax < 2) XLAL_ERROR_VOID(XLAL_FAILURE, "Got LMax=%d < 2!\n", LMax);
    if (*wave) NRSur7dq2_MultiModalWaveform_Destroy(*wave);
    (*wave) = XLALCalloc(1, sizeof(MultiModalWaveform));

    int n_modes = LMax*(LMax+2) - 3;

    (*wave)->n_modes = n_modes;
    (*wave)->lvals = XLALCalloc(n_modes, sizeof(int));
    (*wave)->mvals = XLALCalloc(n_modes, sizeof(int));
    (*wave)->n_times = n_times;
    gsl_vector **modes_real_part = malloc(n_modes * sizeof(*modes_real_part));
    gsl_vector **modes_imag_part = malloc(n_modes * sizeof(*modes_imag_part));
    (*wave)->modes_real_part = modes_real_part;
    (*wave)->modes_imag_part = modes_imag_part;

    int ell=2;
    int m=-2;
    for (int i=0; i<n_modes; i++) {
        (*wave)->lvals[i] = ell;
        (*wave)->mvals[i] = m;
        (*wave)->modes_real_part[i] = gsl_vector_calloc(n_times);
        (*wave)->modes_imag_part[i] = gsl_vector_calloc(n_times);
        m += 1;
        if (m > ell) {
            ell += 1;
            m = -1 * ell;
        }
    }
}

static void NRSur7dq2_MultiModalWaveform_Destroy(MultiModalWaveform *wave) {
    if (!wave) return;
    for (int i=0; i<wave->n_modes; i++) {
        if (wave->modes_real_part[i]) gsl_vector_free(wave->modes_real_part[i]);
        if (wave->modes_imag_part[i]) gsl_vector_free(wave->modes_imag_part[i]);
    }
    free(wave->modes_real_part);
    free(wave->modes_imag_part);
    XLALFree(wave->lvals);
    XLALFree(wave->mvals);
    XLALFree(wave);
}

static void NRSur7dq2_WignerDMatrices_Init(
    WignerDMatrices **matrices,
    int n_times,
    int LMax
) {
    if (!matrices) exit(1);
    if (LMax < 2) XLAL_ERROR_VOID(XLAL_FAILURE, "Got LMax=%d < 2!\n", LMax);
    if (*matrices) NRSur7dq2_WignerDMatrices_Destroy(*matrices);
    (*matrices) = XLALCalloc(1, sizeof(WignerDMatrices));

    int n_entries = 0;
    for (int ell=2; ell<=LMax; ell++) {
        n_entries += (2*ell+1) * (2*ell+1);
    }

    (*matrices)->LMax = LMax;
    (*matrices)->n_entries = n_entries;

    gsl_vector **real_part = malloc(n_entries * sizeof(*real_part));
    gsl_vector **imag_part = malloc(n_entries * sizeof(*imag_part));
    (*matrices)->real_part = real_part;
    (*matrices)->imag_part = imag_part;

    for (int i=0; i<n_entries; i++) {
        (*matrices)->real_part[i] = gsl_vector_calloc(n_times);
        (*matrices)->imag_part[i] = gsl_vector_calloc(n_times);
    }
}

static void NRSur7dq2_WignerDMatrices_Destroy(WignerDMatrices *matrices) {
    if (!matrices) return;
    for (int i=0; i<matrices->n_entries; i++) {
        if (matrices->real_part[i]) gsl_vector_free(matrices->real_part[i]);
        if (matrices->imag_part[i]) gsl_vector_free(matrices->imag_part[i]);
    }
    free(matrices->real_part);
    free(matrices->imag_part);
    XLALFree(matrices);
}

static int NRSur7dq2_WignerDMatrix_Index(int ell, int m, int mp) {
    int i0 = (ell*(ell*ell*4 - 1))/3 - 10; // Start of the (m, mp) matrix, which has size (2*ell + 1) X (2*ell + 1)
    int res = i0 + (2*ell+1)*(ell+m) + (ell+mp);
    return res;
}

static void NRSur7dq2_ComplexPowers_Init(ComplexPowers **cp, int LMax, int n_times) {
    // include powers from -2*LMax to 2*LMax inclusive
    if (!cp) exit(1);
    if (*cp) NRSur7dq2_ComplexPowers_Destroy(*cp);
    (*cp) = XLALCalloc(1, sizeof(ComplexPowers));

    int n_entries = 4*LMax+1;
    (*cp)->LMax = LMax;
    (*cp)->n_entries = n_entries;

    gsl_vector **real_part = malloc(n_entries * sizeof(*real_part));
    gsl_vector **imag_part = malloc(n_entries * sizeof(*imag_part));
    (*cp)->real_part = real_part;
    (*cp)->imag_part = imag_part;

    for (int i=0; i<n_entries; i++) {
        (*cp)->real_part[i] = gsl_vector_calloc(n_times);
        (*cp)->imag_part[i] = gsl_vector_calloc(n_times);
    }
}

static void NRSur7dq2_ComplexPowers_Destroy(ComplexPowers *cp) {
    if (!cp) return;
    for (int i=0; i<cp->n_entries; i++) {
        if (cp->real_part[i]) gsl_vector_free(cp->real_part[i]);
        if (cp->imag_part[i]) gsl_vector_free(cp->imag_part[i]);
    }
    free(cp->real_part);
    free(cp->imag_part);
    XLALFree(cp);
}

static void NRSur7dq2_RealPowers_Init(RealPowers **rp, int LMax, int n_times) {
    // include powers from 0 to 2*LMax inclusive
    if (!rp) exit(1);
    if (*rp) NRSur7dq2_RealPowers_Destroy(*rp);
    (*rp) = XLALCalloc(1, sizeof(RealPowers));

    int n_entries = 2*LMax+1;
    (*rp)->LMax = LMax;
    (*rp)->n_entries = n_entries;

    gsl_vector **powers = malloc(n_entries * sizeof(*powers));
    (*rp)->powers = powers;

    for (int i=0; i<n_entries; i++) {
        (*rp)->powers[i] = gsl_vector_calloc(n_times);
    }
}

static void NRSur7dq2_RealPowers_Destroy(RealPowers *rp) {
    if (!rp) return;
    for (int i=0; i<rp->n_entries; i++) {
        if (rp->powers[i]) gsl_vector_free(rp->powers[i]);
    }
    free(rp->powers);
    XLALFree(rp);
}

static void NRSur7dq2_ComplexPowers_Compute(ComplexPowers *cp, gsl_vector *x, gsl_vector *y) {
    int i, j, power;

    // z = x + iy
    gsl_vector *tmp = gsl_vector_calloc(x->size);
    gsl_vector *z_mag_sqr = gsl_vector_calloc(x->size);

    // Set zero'th power
    i = 2*cp->LMax;
    gsl_vector_add_constant(cp->real_part[i], 1.0);

    // Set first power
    i += 1;
    gsl_vector_add(cp->real_part[i], x);
    gsl_vector_add(cp->imag_part[i], y);

    // Compute positive powers
    for (power=2; power <= 2*cp->LMax; power++) {
        // Compute z^n = z^{n-1} * z. Currently, i indexes z^{n-1}.
        // Re[z^{n-1]] * Re[z]
        gsl_vector_add(tmp, cp->real_part[i]);
        gsl_vector_mul(tmp, x);
        gsl_vector_add(cp->real_part[i+1], tmp);
        gsl_vector_set_zero(tmp);

        // Re[z^{n-1]] * Im[z]
        gsl_vector_add(tmp, cp->real_part[i]);
        gsl_vector_mul(tmp, y);
        gsl_vector_add(cp->imag_part[i+1], tmp);
        gsl_vector_set_zero(tmp);

        // Im[z^{n-1]] * Re[z]
        gsl_vector_add(tmp, cp->imag_part[i]);
        gsl_vector_mul(tmp, x);
        gsl_vector_add(cp->imag_part[i+1], tmp);
        gsl_vector_set_zero(tmp);

        // Im[z^{n-1]] * Im[z]
        gsl_vector_add(tmp, cp->imag_part[i]);
        gsl_vector_mul(tmp, y);
        gsl_vector_sub(cp->real_part[i+1], tmp);
        gsl_vector_set_zero(tmp);

        i += 1;
    }

    // Compute z^{-n} = (z^n)* / |(z^n)^2|
    for (power=1; power <= 2*cp->LMax; power++) {
        i = 2*cp->LMax + power;
        j = 2*cp->LMax - power;

        // Compute |(z^n)|^2
        gsl_vector_add(tmp, cp->real_part[i]);
        gsl_vector_mul(tmp, tmp);
        gsl_vector_add(z_mag_sqr, tmp);
        gsl_vector_set_zero(tmp);
        gsl_vector_add(tmp, cp->imag_part[i]);
        gsl_vector_mul(tmp, tmp);
        gsl_vector_add(z_mag_sqr, tmp);
        gsl_vector_set_zero(tmp);

        // Set z^{-n}
        gsl_vector_add(cp->real_part[j], cp->real_part[i]);
        gsl_vector_div(cp->real_part[j], z_mag_sqr);
        gsl_vector_sub(cp->imag_part[j], cp->imag_part[i]);
        gsl_vector_div(cp->imag_part[j], z_mag_sqr);

        gsl_vector_set_zero(z_mag_sqr);
    }

    gsl_vector_free(tmp);
    gsl_vector_free(z_mag_sqr);

}

static void NRSur7dq2_RealPowers_Compute(RealPowers *rp, gsl_vector *x) {
    int power;

    gsl_vector *tmp = gsl_vector_calloc(x->size);

    // Set zero'th power
    gsl_vector_add_constant(rp->powers[0], 1.0);

    // Set first power
    gsl_vector_add(rp->powers[1], x);

    // Compute positive powers
    for (power=2; power <= 2*rp->LMax; power++) {
        // Compute x^n = x^{n-1} * x.
        gsl_vector_add(tmp, rp->powers[power-1]);
        gsl_vector_mul(tmp, x);
        gsl_vector_add(rp->powers[power], tmp);
        gsl_vector_set_zero(tmp);
    }

    gsl_vector_free(tmp);
}

static void NRSur7dq2_WignerDMatrices_Compute(WignerDMatrices *matrices, gsl_vector **quat) {
    // We compute them all at once because a lot of the work is just processing quat.
    // Parts of this function are adapted from GWFrames:
    // https://github.com/moble/GWFrames
    // written by Michael Boyle, based on his paper:
    // http://arxiv.org/abs/1302.2919
    // although we are working with the conjugate of quat (replace q with [q[0], -1*q[1], -1*q[2], -1*q[3]]

    int n_times = quat[0]->size;
    int i, j, ell, m, mp, rho_min, rho_max, rho;
    double tmp_re, tmp_im, tmp_abs_sqr, coef, c;
    double eps_sqr = 1.0e-24;

    // ra = q[0] - I*q[3], and rb = -q[2] - I*q[1]
    // It's possible that |ra| = 0 or |rb| = 0, making it unsafe to divide by them.
    // We keep track of which indices have this problem, use a safe-but-garbage value of quat=[0.5, 0.5, 0.5, 0.5] at those indices,
    // and then fix them afterwards.
    gsl_vector_int *index_types = gsl_vector_int_calloc(n_times);
    gsl_vector *safe_ra_mag_sqr = gsl_vector_calloc(n_times);
    gsl_vector *safe_rb_mag_sqr = gsl_vector_calloc(n_times);
    gsl_vector *safe_ra_real = gsl_vector_calloc(n_times);
    gsl_vector *safe_ra_imag = gsl_vector_calloc(n_times);
    gsl_vector *safe_rb_real = gsl_vector_calloc(n_times);
    gsl_vector *safe_rb_imag = gsl_vector_calloc(n_times);
    gsl_vector_add(safe_ra_real, quat[0]);
    gsl_vector_sub(safe_ra_imag, quat[3]);
    gsl_vector_sub(safe_rb_real, quat[2]);
    gsl_vector_sub(safe_rb_imag, quat[1]);
    for (i=0; i<n_times; i++) {
        tmp_re = gsl_vector_get(quat[2], i);
        tmp_im = gsl_vector_get(quat[1], i);
        tmp_abs_sqr = tmp_re*tmp_re + tmp_im*tmp_im;
        if (tmp_abs_sqr < eps_sqr) {
            gsl_vector_int_set(index_types, i, 2); // Use 0 for normal, 1 for ra small, 2 for ra ok but rb small
            gsl_vector_set(safe_rb_mag_sqr, i, 0.5);
            gsl_vector_set(safe_rb_real, i, 0.25);
            gsl_vector_set(safe_rb_imag, i, 0.25);
        } else {
            gsl_vector_set(safe_rb_mag_sqr, i, tmp_abs_sqr);
        }

        tmp_re = gsl_vector_get(quat[0], i);
        tmp_im = gsl_vector_get(quat[3], i);
        tmp_abs_sqr = tmp_re*tmp_re + tmp_im*tmp_im;
        if (tmp_abs_sqr < eps_sqr) {
            gsl_vector_int_set(index_types, i, 1); // Use 0 for normal, 1 for ra small, 2 for ra ok but rb small
            gsl_vector_set(safe_ra_mag_sqr, i, 0.5);
            gsl_vector_set(safe_ra_real, i, 0.25);
            gsl_vector_set(safe_ra_imag, i, 0.25);
        } else {
            gsl_vector_set(safe_ra_mag_sqr, i, tmp_abs_sqr);
        }
    }

    // We will need various integer powers of ra, rb, |ra^2|, and |(rb/ra)^2| that will be used in many WignerD terms
    gsl_vector *safe_abs_rb_over_ra_sqr = gsl_vector_calloc(n_times);
    gsl_vector_add(safe_abs_rb_over_ra_sqr, safe_rb_mag_sqr);
    gsl_vector_div(safe_abs_rb_over_ra_sqr, safe_ra_mag_sqr);

    ComplexPowers *ra_powers = NULL;
    ComplexPowers *rb_powers = NULL;
    RealPowers *abs_ra_sqr_powers = NULL;
    RealPowers *abs_rb_over_ra_sqr_powers = NULL;

    NRSur7dq2_ComplexPowers_Init(&ra_powers, matrices->LMax, n_times);
    NRSur7dq2_ComplexPowers_Init(&rb_powers, matrices->LMax, n_times);
    NRSur7dq2_RealPowers_Init(&abs_ra_sqr_powers, matrices->LMax, n_times);
    NRSur7dq2_RealPowers_Init(&abs_rb_over_ra_sqr_powers, matrices->LMax, n_times);


    NRSur7dq2_ComplexPowers_Compute(ra_powers, safe_ra_real, safe_ra_imag);
    NRSur7dq2_ComplexPowers_Compute(rb_powers, safe_rb_real, safe_rb_imag);
    NRSur7dq2_RealPowers_Compute(abs_ra_sqr_powers, safe_ra_mag_sqr);
    NRSur7dq2_RealPowers_Compute(abs_rb_over_ra_sqr_powers, safe_abs_rb_over_ra_sqr);


    // Compute the matrices. We don't use safe_* anymore, so we can use them as temporary results
    for (ell=2; ell <= matrices->LMax; ell++) {
        for (m = -1*ell; m <= ell; m++) {
            for (mp = -1*ell; mp <= ell; mp++) {
                i = NRSur7dq2_WignerDMatrix_Index(ell, m, mp);
                coef = wigner_coef(ell, mp, m);
                gsl_vector_set_zero(safe_ra_mag_sqr);

                // factor = coef * (ra^(m+mp)) * (rb^(m-mp)) * (|ra^2|^(ell-m))
                // The result is
                //      factor * \sum_{rho} c(rho) * (|(rb/ra)^2|^rho)
                //  where c(rho) = (-1)^rho * binom(ell+mp, rho) * binom(ell-mp, ell-rho-m)

                // Store factor in matrices, using safe_rb_real and safe_rb_imag as temporary results
                gsl_vector_add(matrices->real_part[i], ra_powers->real_part[2*matrices->LMax + m + mp]);
                gsl_vector_add(matrices->imag_part[i], ra_powers->imag_part[2*matrices->LMax + m + mp]);
                complex_vector_mult(matrices->real_part[i], matrices->imag_part[i],
                                    rb_powers->real_part[2*matrices->LMax + m - mp], rb_powers->imag_part[2*matrices->LMax + m - mp],
                                    safe_rb_real, safe_rb_imag);
                gsl_vector_scale(matrices->real_part[i], coef);
                gsl_vector_scale(matrices->imag_part[i], coef);
                gsl_vector_mul(matrices->real_part[i], abs_ra_sqr_powers->powers[ell-m]);
                gsl_vector_mul(matrices->imag_part[i], abs_ra_sqr_powers->powers[ell-m]);

                // Compute the sum, storing it in safe_ra_mag_sqr
                rho_min = (mp > m) ? (mp - m) : 0;
                rho_max = (mp > -1*m) ? (ell-m) : (ell+mp);
                for (rho=rho_min; rho <= rho_max; rho++) {
                    c = binomial(ell+mp, rho) * binomial(ell-mp, ell-rho-m);
                    if ((rho%2)==1) {
                        c *= -1;
                    }
                    // Store the temporary term in safe_rb_mag_sqr
                    gsl_vector_set_zero(safe_rb_mag_sqr);
                    gsl_vector_add(safe_rb_mag_sqr, abs_rb_over_ra_sqr_powers->powers[rho]);
                    gsl_vector_scale(safe_rb_mag_sqr, c);
                    gsl_vector_add(safe_ra_mag_sqr, safe_rb_mag_sqr);
                }

                // multiply by the (real) sum to get the result
                gsl_vector_mul(matrices->real_part[i], safe_ra_mag_sqr);
                gsl_vector_mul(matrices->imag_part[i], safe_ra_mag_sqr);

            }
        }
    }

    double zx, zy;
    int k;
    // Now fix the values at bad indices
    for (j=0; j<n_times; j++) {
        rho = gsl_vector_int_get(index_types, j);
        if (rho != 0) {
            for (ell=2; ell <= matrices->LMax; ell++) {
                for (m=-1*ell; m <= ell; m++) {
                    for (mp = -1*ell; mp <= ell; mp++) {
                        i = NRSur7dq2_WignerDMatrix_Index(ell, m, mp);
                        zx = 0.0;
                        zy = 0.0;
                        if (rho==1 && (mp == -m)) {
                            // z = (-1)^(ell+m+1) * rb^(2*m)
                            tmp_re = -1 * gsl_vector_get(quat[2], j);
                            tmp_im = -1 * gsl_vector_get(quat[1], j);
                            zx = 1.0;
                            zy = 0.0;
                            for (k=0; k < 2*m; k++) {
                                c = zx * tmp_im;
                                zx = zx*tmp_re - zy*tmp_im;
                                zy = zy*tmp_re + c;
                            }
                            if ((ell+m)%2 == 0) {
                                zx *= -1;
                                zy *= -1;
                            }
                            gsl_vector_set(matrices->real_part[i], j, zx);
                            gsl_vector_set(matrices->imag_part[i], j, zy);
                        } else if (rho==2 && (mp == m)) {
                            // z = ra^(2*m)
                            tmp_re = gsl_vector_get(quat[0], j);
                            tmp_im = -1 * gsl_vector_get(quat[3], j);
                            zx = 1.0;
                            zy = 0.0;
                            for (k=0; k < 2*m; k++) {
                                c = zx * tmp_im;
                                zx = zx*tmp_re - zy*tmp_im;
                                zy = zy*tmp_re + c;
                            }
                            gsl_vector_set(matrices->real_part[i], j, zx);
                            gsl_vector_set(matrices->imag_part[i], j, zy);
                        }
                        gsl_vector_set(matrices->real_part[i], j, zx);
                        gsl_vector_set(matrices->imag_part[i], j, zy);
                    }
                }
            }
        }
    }

    // Cleanup
    gsl_vector_int_free(index_types);
    gsl_vector_free(safe_ra_mag_sqr);
    gsl_vector_free(safe_rb_mag_sqr);
    gsl_vector_free(safe_ra_real);
    gsl_vector_free(safe_ra_imag);
    gsl_vector_free(safe_rb_real);
    gsl_vector_free(safe_rb_imag);
    gsl_vector_free(safe_abs_rb_over_ra_sqr);
    NRSur7dq2_ComplexPowers_Destroy(ra_powers);
    NRSur7dq2_ComplexPowers_Destroy(rb_powers);
    NRSur7dq2_RealPowers_Destroy(abs_ra_sqr_powers);
    NRSur7dq2_RealPowers_Destroy(abs_rb_over_ra_sqr_powers);
}

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
    DynamicsNodeFitData **ds_node_data = malloc( (t_ds->size) * sizeof(*ds_node_data));
    DynamicsNodeFitData **ds_half_node_data = malloc( 3 * sizeof(*ds_node_data) );
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
    WaveformFixedEllModeData **coorbital_mode_data = malloc( (NRSUR7DQ2_LMAX - 1) * sizeof(*coorbital_mode_data) );
    for (int ell_idx=0; ell_idx < NRSUR7DQ2_LMAX-1; ell_idx++) {
        NRSur7dq2_LoadCoorbitalEllModes(coorbital_mode_data, file, ell_idx);
    }
    data->coorbital_mode_data = coorbital_mode_data;

    printf("Successfully loaded NRSur7dq2 data!\n");
    data->LMax = NRSUR7DQ2_LMAX;
    data->setup = 1;

    return XLAL_SUCCESS;
}

static void NRSur7dq2_LoadDynamicsNode(
    DynamicsNodeFitData **ds_node_data,     // Entry i should be NULL; Will malloc space and load data
    LALH5File *sub,                         // Subgroup containing data for dynamics node i
    int i                                   // Dynamics node index. Note that since ) {
) {
    ds_node_data[i] = malloc( sizeof(*ds_node_data[i]) );

    // omega
    FitData *omega_data = malloc(sizeof(FitData));
    omega_data->coefs = NULL;
    omega_data->basisFunctionOrders = NULL;
    ReadHDF5RealVectorDataset(sub, "omega_coefs", &(omega_data->coefs));
    ReadHDF5LongMatrixDataset(sub, "omega_bfOrders", &(omega_data->basisFunctionOrders));
    omega_data->n_coefs = omega_data->coefs->size;
    ds_node_data[i]->omega_data = omega_data;

    // omega_copr
    VectorFitData *omega_copr_data = malloc(sizeof(VectorFitData));
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
    VectorFitData *chiA_dot_data = malloc(sizeof(VectorFitData));
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
    VectorFitData *chiB_dot_data = malloc(sizeof(VectorFitData));
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

static void NRSur7dq2_LoadCoorbitalEllModes(
    WaveformFixedEllModeData **coorbital_mode_data, // All coorbital waveform modes; Space will be allocated and data will be loaded for entry i.
    LALH5File *file, // NRSur7dq2.hdf5
    int i // The index of coorbital_mode_data, also ell-2.
) {
    WaveformFixedEllModeData *mode_data = malloc( sizeof(*coorbital_mode_data[i]) );
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
    mode_data->X_real_plus_data = malloc( (i+2) * sizeof(WaveformDataPiece *) );
    mode_data->X_real_minus_data = malloc( (i+2) * sizeof(WaveformDataPiece *) );
    mode_data->X_imag_plus_data = malloc( (i+2) * sizeof(WaveformDataPiece *) );
    mode_data->X_imag_minus_data = malloc( (i+2) * sizeof(WaveformDataPiece *) );
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

static void NRSur7dq2_LoadWaveformDataPiece(
    LALH5File *sub, // HDF5 group containing data for this waveform data piece
    WaveformDataPiece **data, // Output
    bool invert_sign
) {
    *data = malloc(sizeof(WaveformDataPiece));

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
    (*data)->fit_data = malloc( n_nodes * sizeof(FitData *) );

    LALH5File *nodeModelers = XLALH5GroupOpen(sub, "nodeModelers");
    int str_size = 20; // Enough for L with 11 digits...
    char *sub_name = XLALMalloc(str_size);
    for (int i=0; i<n_nodes; i++) {
        FitData *node_data = malloc(sizeof(FitData));
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

static bool NRSur7dq2_IsSetup(void) {
    if(__lalsim_NRSur7dq2_data.setup)   return true;
    else return false;
}

double ipow(double base, int exponent) {
    if (exponent == 0) return 1.0;
    double res = base;
    while (exponent > 1) {
        res = res*base;
        exponent -= 1;
    }
    return res;
}

// This is just a helper to save space.
// It's probably better to use complex vectors but I'm almost done and don't want to rewrite everything...
// Multiplies (x1 + I*y1) * (x2 + I*y2), storing the real part in x1 and the imaginary part in x2.
// Uses tmp1 and tmp2 for temporary results; x2 and y2 will be left unchanged but tmpx and tmpy are modified.
static void complex_vector_mult(gsl_vector *x1, gsl_vector *y1, gsl_vector *x2, gsl_vector *y2, gsl_vector *tmpx, gsl_vector *tmpy) {
    gsl_vector_set_zero(tmpx);
    gsl_vector_set_zero(tmpy);

    // Store the temporary results x1*y2 and y1*y2
    gsl_vector_add(tmpx, y1);
    gsl_vector_mul(tmpx, y2);
    gsl_vector_add(tmpy, x1);
    gsl_vector_mul(tmpy, y2);

    // Now we can safely transform x1->x1*x2 and y1->y1*x2
    gsl_vector_mul(x1, x2);
    gsl_vector_mul(y1, x2);

    // Now we add in the temporary results
    gsl_vector_sub(x1, tmpx);
    gsl_vector_add(y1, tmpy);
}

/*
 * The fit result is given by
 *      \sum_{i=1}^{n} c_i * \prod_{j=1}^7 B_j(k_{i, j}; x_j)
 * where i runs over fit coefficients, j runs over the 7 dimensional parameter
 * space, and B_j is a basis function, taking an integer order k_{i, j} and
 * the parameter component x_j. For this surrogate, B_j are monomials in the spin
 * components, and monomials in an affine transformation of the mass ratio.
 */
double NRSur7dq2_eval_fit(
    FitData *data, // Data for fit
    double *x // size 7, giving mass ratio q, and dimensionless spin components
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
        prod = x_powers[7 * gsl_matrix_long_get(data->basisFunctionOrders, i, 0)]; // Initialize with q basis function
        for (j=1; j<7; j++) { // Multiply with spin basis functions
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
    double *res, // Result
    VectorFitData *data, // Data for fit
    double *x // size 7, giving mass ratio q, and dimensionless spin components
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
        prod = x_powers[7 * gsl_matrix_long_get(data->basisFunctionOrders, i, 0)]; // Initialize with q basis function
        for (j=1; j<7; j++) { // Multiply with spin basis functions
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
    double chiANorm, // ||\vec{\chi}_A||
    double chiBNorm, // ||\vec{\chi}_B||
    double *y // [q0, qx, qy, qz, phi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
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
    double normA,
    double normB,
    gsl_vector **quat,
    gsl_vector **chiA,
    gsl_vector **chiB
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

static void NRSur7dq2_rotate_spins(gsl_vector **chiA, gsl_vector **chiB, gsl_vector *phi_vec, bool backwards) {

    int i;
    int n = phi_vec->size;
    double sp, cp, tmp;
    double *phi = phi_vec->data;
    int sgn = 1;
    if (backwards) sgn = -1;

    double *chix = chiA[0]->data;
    double *chiy = chiA[1]->data;
    for (i=0; i<n; i++) {
        cp = cos(phi[i]);
        sp = sin(sgn * phi[i]);
        tmp = chix[i];
        chix[i] = tmp*cp + chiy[i]*sp;
        chiy[i] = -1*tmp*sp + chiy[i]*cp;
    }

    chix = chiB[0]->data;
    chiy = chiB[1]->data;
    for (i=0; i<n; i++) {
        cp = cos(phi[i]);
        sp = sin(sgn * phi[i]);
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
    double *x, // Result, length 7
    double q, // Mass ratio
    double *y // [q0, qx, qy, qz, phi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
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
    double *dydt,           // Result, length 11
    double *y,              // ODE solution at the current time step
    double *Omega_coorb_xy, // a form of time derivative of the coprecessing frame
    double omega,           // orbital angular frequency in the coprecessing frame
    double *chiA_dot,       // chiA time derivative
    double *chiB_dot        // chiB time derivative
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

/*
 * Given dy/dt evaluated as a function of y and t at 4 time nodes, compute
 * dy = y_5 - y_4 using a 4th-order accurate Adams-Bashforth scheme.
 */
static void NRSur7dq2_ab4_dy(
    double *dy,     // Result, length 11
    double *k1,     // dy/dt evaluated at node 1
    double *k2,     // dy/dt evaluated at node 2
    double *k3,     // dy/dt evaluated at node 3
    double *k4,     // dy/dt evaluated at node 4
    double dt1,     // t_2 - t_1
    double dt2,     // t_3 - t_2
    double dt3,     // t_4 - t_3
    double dt4      // t_5 - t_4
) {

    double dt12, dt23, dt123, D1, D2, D3, B41, B42, B43, B4, C41, C42, C43, C4;
    double A, B, C, D;
    int i;

    // Various time intervals
    dt12 = dt1 + dt2;
    dt123 = dt12 + dt3;
    dt23 = dt2 + dt3;

    // Denomenators and coefficients
    D1 = dt1 * dt12 * dt123;
    D2 = dt1 * dt2 * dt23;
    D3 = dt2 * dt12 * dt3;

    B41 = dt3 * dt23 / D1;
    B42 = -1 * dt3 * dt123 / D2;
    B43 = dt23 * dt123 / D3;
    B4 = B41 + B42 + B43;

    C41 = (dt23 + dt3) / D1;
    C42 = -1 * (dt123 + dt3) / D2;
    C43 = (dt123 + dt23) / D3;
    C4 = C41 + C42 + C43;

    // Polynomial coefficients and result
    for (i=0; i<11; i++) {
        A = k4[i];
        B = k4[i]*B4 - k1[i]*B41 - k2[i]*B42 - k3[i]*B43;
        C = k4[i]*C4 - k1[i]*C41 - k2[i]*C42 - k3[i]*C43;
        D = (k4[i]-k1[i])/D1 - (k4[i]-k2[i])/D2 + (k4[i]-k3[i])/D3;
        dy[i] = dt4 * (A + dt4 * (0.5*B + dt4*( C/3.0 + dt4*0.25*D)));
    }
}

static double factorial(int n) {
    if (n <= 0) return 1.0;
    return factorial(n-1) * n;
}

static double factorial_ratio(int n, int k) {
    if (n <= k) return 1.0;
    return factorial_ratio(n-1, k) * n;
}

static double binomial(int n, int k) {
    return factorial_ratio(n, k) / factorial(n-k);
}

static double wigner_coef(int ell, int mp, int m) {
    return sqrt(factorial(ell+m) * factorial(ell-m) / (factorial(ell+mp) * factorial(ell-mp)));
}

static double cubic_interp(double xout, double *x, double *y) {
    // This gives a much closer result to scipy.interpolate.InterpolatedUnivariateSpline than using gsl_interp_cspline
    // (see comment below in spline_array_interp)
    // x and y should have length 4, with x increasing.
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *interpolant = gsl_spline_alloc(gsl_interp_polynomial, 4);
    gsl_spline_init(interpolant, x, y, 4);
    double res = gsl_spline_eval(interpolant, xout, acc);
    gsl_spline_free(interpolant);
    gsl_interp_accel_free(acc);
    return res;
}

static gsl_vector *spline_array_interp(gsl_vector *xout, gsl_vector *x, gsl_vector *y) {
    // Results differ from scipy.interpolate.InterpolatedUnivariateSpline due to different boundary conditions.
    // We should really implement fast not-a-knot 1d spline interpolation and use it everywhere instead of gsl.
    // This difference leads to differences between this implementation of NRSur7dq2 and the python implementation.
    // As far as I'm aware, it's the only difference, and everything else agrees to within machine precision.
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

static double NRSur7dq2_get_omega(
    size_t node_index,
    double q,
    double *y0
){
    double x[7];
    NRSur7dq2_ds_fit_x(x, q, y0);
    FitData *data = (&__lalsim_NRSur7dq2_data)->ds_node_data[node_index]->omega_data;
    double omega = NRSur7dq2_eval_fit(data, x);
    return omega;
}

static double NRSur7dq2_get_t_ref(
    double omega_ref,   // reference orbital angular frequency
    double q,           // mass ratio
    double *chiA0,      // chiA at reference point
    double *chiB0,      // chiB at reference point
    double *q_ref,      // coprecessing frame quaternion at reference point
    double phi_ref      // orbital phase at reference point
) {
    if (fabs(omega_ref) < 1.e-10) {
        XLAL_PRINT_WARNING("WARNING: Treating omega_ref = 0 as a flag to use t_ref = t_0 = -4500M");
        return -4499.99999999; // The first time node is ever so slightly larger than -4500; this avoids out-of-range issues
    }

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

static void NRSur7dq2_get_time_deriv_from_index(
    double *dydt,       // Output: dy/dt evaluated at the ODE time node with index i0. Must have space for 11 entries.
    int i0,             // Time node index. i0=-1, -2, and -3 are used for time nodes 1/2, 3/2, and 5/2 respectively.
    double q,           // Mass ratio
    double *y           // Current ODE state: [q0, qx, qy, qz, orbphase, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
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

static void NRSur7dq2_get_time_deriv(
    double *dydt,       // Output: dy/dt evaluated at time t. Must have space for 11 entries.
    double t,           // Time at which the ODE should be evaluated.
    double q,           // Mass ratio
    double *y           // Current ODE state
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

static int NRSur7dq2_initialize_at_dynamics_node(
    double *dynamics_data,  // ODE output
    double t_ref,           // reference time. t_ds[i0] will be close to t_ref.
    double q,               // mass ratio
    double *chiA0,          // chiA at t_ref.
    double *chiB0,          // chiB at t_ref.
    double phi_ref,         // orbital phase at t_ref.
    double *q_ref,          // quaternion at t_ref.
    double normA,           // |chiA|
    double normB            // |chiB|
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

static void NRSur7dq2_initialize_RK4_with_half_nodes(
    double *dynamics_data,  // ODE output
    double *time_steps,     // Output: first three time steps. Should have size 3.
    double *dydt0,             // Output: dydt at node 0. Should have size 11.
    double *dydt1,             // Output: dydt at node 1. Should have size 11.
    double *dydt2,             // Output: dydt at node 2. Should have size 11.
    double *dydt3,             // Output: dydt at node 3. Should have size 11.
    double normA,           // |chiA|
    double normB,           // |chiB|
    double q                // mass ratio
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

static int NRSur7dq2_initialize_RK4(
    double *dynamics_data,  // ODE output - quaternions
    double *time_steps,     // Output: first three time steps. Should have size 3.
    double *dydt0,             // Output: dydt at node i0 + 0. Should have size 11.
    double *dydt1,             // Output: dydt at node i0 + 1. Should have size 11.
    double *dydt2,             // Output: dydt at node i0 + 2. Should have size 11.
    double *dydt3,             // Output: dydt at node i0 + 3. Should have size 11.
    double normA,           // |chiA|
    double normB,           // |chiB|
    double q,               // mass ratio
    int i0                  // the node that is already initialized
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
                y_tmp[j] = node_data[j] - (time_steps[i]/6.0)*(dydt[3-i][j] + 2*k2[j] + 2*k3[j] + k4[j]);
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

static void NRSur7dq2_integrate_AB4(
    double *dynamics_data,  // ODE output
    double *time_steps,     // The first three time steps beginning at i_start.
    double *dydt0,          // dydt at node i_start
    double *dydt1,          // dydt at node i_start + 1
    double *dydt2,          // dydt at node i_start + 2
    double *dydt3,          // dydt at node i_start + 3
    double normA,           // |chiA|
    double normB,           // |chiB|
    double q,               // mass ratio
    int i_start             // nodes i_start through i_start+3 are already initialized
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
        NRSur7dq2_ab4_dy(ynext, k1, k2, k3, k4, dt1, dt2, dt3, dt4);

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
        NRSur7dq2_ab4_dy(ynext, k1, k2, k3, k4, dt1, dt2, dt3, dt4);

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

static void NRSur7dq2_eval_data_piece(
    gsl_vector *result, // Output: Should have already been assigned space
    double q,           // Mass ratio
    gsl_vector **chiA,  // 3 gsl_vector *s, one for each (coorbital) component
    gsl_vector **chiB,  // similar to chiA
    WaveformDataPiece *data // The data piece to evaluate
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

static int NRSur7dq2_IntegrateDynamics(
    double *dynamics_data,  // Output: length (n * 11), where entries 11*m <= i < 11*(m+1) are
                            // [q0, qx, qy, qz, varphi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
    double q,               // Mass ratio mA / mB
    double *chiA0,          // chiA at the reference point
    double *chiB0,          // chiB at the reference point
    double omega_ref,       // orbital angular frequency at the reference point
    double phi_ref,         // orbital phase at the reference point
    double *q_ref           // coprecessing quaterion at the reference point
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

static void NRSur7dq2_TransformModes(
    MultiModalWaveform *h,          // Output. Dimensionless waveform modes sampled on t_coorb. Should be initialized.
    MultiModalWaveform *h_coorb,    // Coorbital frame waveform modes.
    gsl_vector **quat,              // should be *quat[4], one for each component. Coprecessing frame quaternions.
    gsl_vector *orbphase            // Orbital phase, used to obtain the coprecessing modes.
) {
    int i, j, ell, m, mp;
    int n_times = h_coorb->n_times;
    int lmax = h_coorb->lvals[h_coorb->n_modes -1];
    int lmin;
    gsl_vector *cosmphi = gsl_vector_alloc(n_times);
    gsl_vector *sinmphi = gsl_vector_alloc(n_times);
    gsl_vector *tmp_vec = gsl_vector_alloc(n_times);

    // First transform to the coprecessing frame:
    // h^{\ell, m}_\mathrm{copr} = e^{-i m \varphi} h^{\ell, m}_\mathrm{coorb}.
    // Loop over m first, since fixed m modes transform the same way.
    MultiModalWaveform *h_copr = NULL;
    NRSur7dq2_MultiModalWaveform_Init(&h_copr, lmax, n_times);
    for (m = -1*lmax; m <= lmax; m++) {
        if (m==0) {
            // No transformation
            for (ell = 2; ell <= lmax; ell++) {
                i = ell*(ell+1) - 4;
                gsl_vector_add(h_copr->modes_real_part[i], h_coorb->modes_real_part[i]);
                gsl_vector_add(h_copr->modes_imag_part[i], h_coorb->modes_imag_part[i]);
            }
        } else {
            for (j=0; j<n_times; j++) {
                gsl_vector_set(cosmphi, j, cos(m * gsl_vector_get(orbphase, j)));
                gsl_vector_set(sinmphi, j, sin(m * gsl_vector_get(orbphase, j)));
            }
            lmin = abs(m);
            if (lmin < 2) {
                lmin = 2;
            }
            for (ell = lmin; ell<=lmax; ell++) {
                i = ell*(ell+1) - 4 + m;
                // compute and add multiplicative combinations of {real, imag} and {cosmphi, sinmphi}
                // real * cosmphi
                gsl_vector_set_zero(tmp_vec);
                gsl_vector_add(tmp_vec, h_coorb->modes_real_part[i]);
                gsl_vector_mul(tmp_vec, cosmphi);
                gsl_vector_add(h_copr->modes_real_part[i], tmp_vec);
                // -1 * real * I * sinmphi
                gsl_vector_set_zero(tmp_vec);
                gsl_vector_add(tmp_vec, h_coorb->modes_real_part[i]);
                gsl_vector_mul(tmp_vec, sinmphi);
                gsl_vector_sub(h_copr->modes_imag_part[i], tmp_vec);
                // I * imag * cosmphi
                gsl_vector_set_zero(tmp_vec);
                gsl_vector_add(tmp_vec, h_coorb->modes_imag_part[i]);
                gsl_vector_mul(tmp_vec, cosmphi);
                gsl_vector_add(h_copr->modes_imag_part[i], tmp_vec);
                // imag * sinmphi
                gsl_vector_set_zero(tmp_vec);
                gsl_vector_add(tmp_vec, h_coorb->modes_imag_part[i]);
                gsl_vector_mul(tmp_vec, sinmphi);
                gsl_vector_add(h_copr->modes_real_part[i], tmp_vec);
            }
        }
    }

    // Now we transform the coprecessing modes to the intertial frame, using the quaternions.
    WignerDMatrices *matrices = NULL;
    NRSur7dq2_WignerDMatrices_Init(&matrices, n_times, lmax);
    NRSur7dq2_WignerDMatrices_Compute(matrices, quat);
    int matrix_index;
    for (ell = 2; ell <= lmax; ell++) {
        for (m= -1 * ell; m <= ell; m++) {
            i = ell*(ell+1) - 4 + m;
            for (mp = -1 * ell; mp <= ell; mp++) {
                j = ell*(ell+1) - 4 + mp;
                matrix_index = NRSur7dq2_WignerDMatrix_Index(ell, m, mp);
                // Re[h] * Re[D]
                gsl_vector_set_zero(tmp_vec);
                gsl_vector_add(tmp_vec, h_copr->modes_real_part[j]);
                gsl_vector_mul(tmp_vec, matrices->real_part[matrix_index]);
                gsl_vector_add(h->modes_real_part[i], tmp_vec);

                // Re[h] * Im[D]
                gsl_vector_set_zero(tmp_vec);
                gsl_vector_add(tmp_vec, h_copr->modes_real_part[j]);
                gsl_vector_mul(tmp_vec, matrices->imag_part[matrix_index]);
                gsl_vector_add(h->modes_imag_part[i], tmp_vec);

                // Im[h] * Re[D]
                gsl_vector_set_zero(tmp_vec);
                gsl_vector_add(tmp_vec, h_copr->modes_imag_part[j]);
                gsl_vector_mul(tmp_vec, matrices->real_part[matrix_index]);
                gsl_vector_add(h->modes_imag_part[i], tmp_vec);

                // Im[h] * Im[D]
                gsl_vector_set_zero(tmp_vec);
                gsl_vector_add(tmp_vec, h_copr->modes_imag_part[j]);
                gsl_vector_mul(tmp_vec, matrices->imag_part[matrix_index]);
                gsl_vector_sub(h->modes_real_part[i], tmp_vec);
            }
        }
    }
    NRSur7dq2_WignerDMatrices_Destroy(matrices);
    NRSur7dq2_MultiModalWaveform_Destroy(h_copr);
    gsl_vector_free(cosmphi);
    gsl_vector_free(sinmphi);
    gsl_vector_free(tmp_vec);
}

/* This is the core function of the NRSur7dq2 model.
 * It evaluates the model, and returns waveform modes in the inertial frame, sampled on
 * t_coorb.
 */
static void NRSur7dq2_core(
    MultiModalWaveform **h, // Output. Dimensionless waveform modes sampled on t_coorb
    double q,               // Mass ratio mA / mB
    double *chiA0,          // chiA at the reference point
    double *chiB0,          // chiB at the reference point
    double omega_ref,       // orbital angular frequency at the reference point
    double phi_ref,         // orbital phase at the reference point
    double *q_ref,          // coprecessing quaterion at the reference point
    int LMax                // Maximum ell mode to evaluate (NRSur7dq2 contains ell=2, 3, 4 modes)
) {
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&NRSur7dq2_is_initialized, NRSur7dq2_Init_LALDATA);
#else
    NRSur7dq2_Init_LALDATA();
#endif

    // time arrays
    gsl_vector *t_ds = (&__lalsim_NRSur7dq2_data)->t_ds;
    gsl_vector *t_coorb = (&__lalsim_NRSur7dq2_data)->t_coorb;
    int n_ds = t_ds->size;
    int n_coorb = t_coorb->size;

    // Get dynamics
    double *dynamics_data = XLALCalloc(n_ds * 11, sizeof(double));

    int i, j;
    int ret = NRSur7dq2_IntegrateDynamics(dynamics_data, q, chiA0, chiB0, omega_ref, phi_ref, q_ref);
    if(ret != XLAL_SUCCESS) XLAL_ERROR_VOID(XLAL_FAILURE, "Failed to integrate dynamics");

    // Put output into appropriate vectors
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
    NRSur7dq2_rotate_spins(chiA_coorb, chiB_coorb, phi_coorb, false);

    // Evaluate the coorbital waveform surrogate
    MultiModalWaveform *h_coorb = NULL;
    NRSur7dq2_MultiModalWaveform_Init(&h_coorb, LMax, n_coorb);
    gsl_vector *data_piece_eval = gsl_vector_alloc(n_coorb);
    WaveformDataPiece *data_piece_data;
    int ell, m;
    int i0; // for indexing the (ell, m=0) mode, such that the (ell, m) mode is index (i0 + m).
    WaveformFixedEllModeData *ell_data;
    for (ell=2; ell<=LMax; ell++) {
        ell_data = (&__lalsim_NRSur7dq2_data)->coorbital_mode_data[ell - 2];
        i0 = ell*(ell+1) - 4;

        // m=0
        data_piece_data = ell_data->m0_real_data;
        NRSur7dq2_eval_data_piece(data_piece_eval, q, chiA_coorb, chiB_coorb, data_piece_data);
        gsl_vector_add(h_coorb->modes_real_part[i0], data_piece_eval);

        data_piece_data = ell_data->m0_imag_data;
        NRSur7dq2_eval_data_piece(data_piece_eval, q, chiA_coorb, chiB_coorb, data_piece_data);
        gsl_vector_add(h_coorb->modes_imag_part[i0], data_piece_eval);

        // Other modes
        for (m=1; m<=ell; m++) {
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
    NRSur7dq2_MultiModalWaveform_Init(h, LMax, n_coorb);
    NRSur7dq2_TransformModes(*h, h_coorb, quat_coorb, phi_coorb);

    // Cleanup
    NRSur7dq2_MultiModalWaveform_Destroy(h_coorb);
    for (i=0; i<3; i++) {
        gsl_vector_free(chiA_coorb[i]);
        gsl_vector_free(chiB_coorb[i]);
        gsl_vector_free(quat_coorb[i]);
    }
    gsl_vector_free(quat_coorb[3]);
    gsl_vector_free(phi_coorb);
    gsl_vector_free(data_piece_eval);
}

/** This function evaluates the NRSur7dq2 surrogate model and sums over all ell <= 4 modes to obtain the + and x polarizations.*/
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
        REAL8 s2z                      /**< initial value of S2z */
) {
    if (fabs(fMin) > 1.e-6) XLAL_PRINT_WARNING("NRSur7dq2 ignores fMin. Set fMin=0 to ignore this warning.");

    // TODO: Could lmax, as well as the frame choice at fRef, be obtained from the LALparams dict?
    // If so, how? lmax would be quite useful to have, both to study the effect of lmax and since
    // lmax=2 is more than twice as fast as lmax=4 (and for low SNR lmax=2 is sufficient).

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
    chiA0[0] = s1x;
    chiA0[1] = s1y;
    chiA0[2] = s1z;
    chiB0[0] = s2x;
    chiB0[1] = s2y;
    chiB0[2] = s2z;
    // Orbital angular frequency = 0.5*(wave angular frequency) = pi*(wave frequency)
    double omegaRef_dimless = (Mtot_sec * fRef) * LAL_PI;
    // Take the inertial frame to be aligned with the coprecessing frame at the reference point:
    double q_ref[4];
    q_ref[0] = 1.0;
    q_ref[1] = 0.0;
    q_ref[2] = 0.0;
    q_ref[3] = 0.0;
    int lmax=4; // Use all available modes.

    // Evaluate the model modes
    MultiModalWaveform *h_inertial_modes = NULL;
    NRSur7dq2_core(&h_inertial_modes, q, chiA0, chiB0, omegaRef_dimless, phiRef, q_ref, lmax);
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
    int ell, m, i, j;
    i=0;
    for (ell=2; ell<=lmax; ell++) {
        for (m=-ell; m<=ell; m++) {
            // phi=0. phiRef controls the *orbital phase*.
            // To effectively obtain phi != 0, take phiRef -> phiRef - phi, and rotate the (x + iy) spin components
            // in the complex plane by e^{-i*phi} for both chi1 and chi2.
            swsh_coef = XLALSpinWeightedSphericalHarmonic(inclination, 0.0, -2, ell, m);
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
    double tf = gsl_vector_get(model_times, length-1);
    int nt = (int) ceil((tf - t0) / deltaT_dimless);
    gsl_vector *output_times = gsl_vector_alloc(nt);
    for (j=0; j<nt; j++) {
        gsl_vector_set(output_times, j, t0 + deltaT_dimless*j);
    }

    // Interpolate onto output times
    double t;
    LIGOTimeGPS epoch = LIGOTIMEGPSZERO; // Dummy time
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
    NRSur7dq2_MultiModalWaveform_Destroy(h_inertial_modes);

    return XLAL_SUCCESS;
}

/** This function evaluates the NRSur7dq2 surrogate model and returns the inertial frame modes.*/
SphHarmTimeSeries *XLALSimInspiralNRSur7dq2Modes(
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 deltaT,                   /**< sampling interval (s) */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 fMin,                     /**< start GW frequency (Hz) */
        REAL8 fRef,                     /**< reference GW frequency (Hz) */
        REAL8 distance,                 /**< distance of source (m) */
        int lmax                        /**< Evaluates (l, m) modes with l <= lmax. The model contains modes up to l=4. */
) {
    if (fabs(fMin) > 1.e-6) XLAL_PRINT_WARNING("NRSur7dq2 ignores fMin. Set fMin=0 to ignore this warning.");

    // TODO: Implement spins somehow?!?!
    SphHarmTimeSeries *hlms = NULL;

    if (lmax > NRSUR7DQ2_LMAX) {
        XLAL_PRINT_WARNING("Got lmax=%d > NRSUR7DQ2_LMAX=%d, decreasing lmax\n", lmax, NRSUR7DQ2_LMAX);
        lmax = NRSUR7DQ2_LMAX;
    }

    double chiA0[3], chiB0[3];
    chiA0[0] = 0.0;
    chiA0[1] = 0.0;
    chiA0[2] = 0.0;
    chiB0[0] = 0.0;
    chiB0[1] = 0.0;
    chiB0[2] = 0.0;

    // BH A is defined to be the one with the larger mass, BH B with the smaller mass.
    if (m2 > m1) {
        double tmp = m1;
        m1 = m2;
        m2 = tmp;
        phiRef += LAL_PI;
        // TODO: If spins are implemented, also switch spin labels here.
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
    NRSur7dq2_core(&h_inertial, q, chiA0, chiB0, omegaRef_dimless, phiRef, q_ref, lmax);
    if (!h_inertial) {
        printf("NRSur7dq2_core failed!\n");
        return hlms;
    }

    // Scale the amplitude based on the distance
    double a0 = Mtot_sec * LAL_C_SI / distance;

    // Setup dimensionless output times
    gsl_vector *model_times = (&__lalsim_NRSur7dq2_data)->t_coorb;
    size_t length = model_times->size;
    double deltaT_dimless = deltaT / Mtot_sec;
    double t0 = gsl_vector_get(model_times, 0);
    double tf = gsl_vector_get(model_times, length-1);
    int nt = (int) ceil((tf - t0) / deltaT_dimless);
    gsl_vector *output_times = gsl_vector_alloc(nt);
    int j;
    for (j=0; j<nt; j++) {
        gsl_vector_set(output_times, j, t0 + deltaT_dimless*j);
    }

    // Setup interpolation onto dimensionless output times
    double t;
    LIGOTimeGPS epoch = LIGOTIMEGPSZERO; // Dummy time
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    // Sum over modes
    COMPLEX16TimeSeries *tmp_mode;
    int ell, m;
    int i=0;
    char mode_name[32];
    for (ell=2; ell<=lmax; ell++) {
        for (m=-ell; m<=ell; m++) {
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
    NRSur7dq2_MultiModalWaveform_Destroy(h_inertial);
    gsl_vector_free(output_times);
    gsl_interp_accel_free(acc);

    return hlms;
}
