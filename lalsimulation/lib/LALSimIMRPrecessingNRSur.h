/* *  Copyright (C) 2017 Jonathan Blackman, Vijay Varma
 *  NRSur7dq2 and NRSur7dq4 NR surrogate models.
 *  Papers: https://arxiv.org/abs/1705.07089, https://arxiv.org/abs/1905.09300.
 *  Based on the python implementation found at:
 *  https://www.black-holes.org/data/surrogates/index.html
 *  which uses the same hdf5 data file.
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

#include "LALSimNRSurrogateUtilities.c"


/****************************** Constants ***********************************/

///// NRSur7dq2 constants

// These are to rescale the mass ratio fit range from [0.99, 2.01] to [-1, 1].
static const double NRSUR7DQ2_Q_FIT_OFFSET = -2.9411764705882355; // = (0.5*(2.01 - 0.99)) * (2 / (2.01 - 0.99))
static const double NRSUR7DQ2_Q_FIT_SLOPE = 1.9607843137254901; // = 2 / (2.01 - 0.99)

// Allow extrapolation in mass ratio and spin beyond these values but throw
// a warning
static const double NRSUR7DQ2_Q_MAX_WARN = 2.01;
static const double NRSUR7DQ2_CHI_MAX_WARN = 0.81;

// Do not allow extrapolation beyond these, unless unlimited_extrapolation=1
static const double NRSUR7DQ2_Q_MAX = 3.01;
static const double NRSUR7DQ2_CHI_MAX = 1;

static const double NRSUR7DQ2_START_TIME = -4499.99999999; // The first time node is ever so slightly larger than -4500; this avoids out-of-range issues;


///// NRSur7dq4 constants

// These are to rescale log(q) fit range from [-0.01, np.log(4.01)] to [-1, 1].
// slope = 2./(np.log(4.01) + 0.01)
static const double NRSUR7DQ4_Q_FIT_SLOPE = 1.4298059216576398  ;
// offset = np.log(4.01)*1.4298059216576398 - 1
static const double NRSUR7DQ4_Q_FIT_OFFSET = -0.9857019407834238;

// Allow extrapolation in mass ratio and spin beyond these values but throw
// a warning
static const double NRSUR7DQ4_Q_MAX_WARN = 4.01;
static const double NRSUR7DQ4_CHI_MAX_WARN = 0.81;

// Do not allow extrapolation beyond these, unless unlimited_extrapolation=1
static const double NRSUR7DQ4_Q_MAX = 6.01;
static const double NRSUR7DQ4_CHI_MAX = 1;

static const double NRSUR7DQ4_START_TIME = -4299.99999999; // The first time node is ever so slightly larger than -4300; this avoids out-of-range issues;





static const int NRSUR_LMAX = 4;

// Surrogate model data, in LAL_DATA_PATH. File available in lalsuite-extra or at
// https://www.black-holes.org/surrogates
static const char NRSUR7DQ2_DATAFILE[] = "NRSur7dq2.h5";
static const char NRSUR7DQ4_DATAFILE[] = "NRSur7dq4.h5";

/***********************************************************************************/
/****************************** Type definitions ***********************************/
/***********************************************************************************/

/**
 * Data used in a single scalar fit
 */
typedef struct tagFitData {
    gsl_matrix_long *basisFunctionOrders;   /**< matrix of (n_coefs x 7) basis function orders
                        giving the polynomial order in f(q), chiA components, and chiB components. */
    gsl_vector *coefs;                      /**< coefficient vector of length n_coefs */
    int n_coefs;                            /**< Number of coefficients in the fit */
} FitData;

/**
 * Data used in a single vector fit
 * NOTE: basisFunctionOrders, coefs, componentIndices, and n_coefs are only
 * used by NRSur7dq2. While fit_data is only used by NRSur7dq4.
 */
typedef struct tagVectorFitData {
    gsl_matrix_long *basisFunctionOrders;   /**< matrix of (n_coefs x 7) basis function orders
                        giving the polynomial order in f(q), chiA components, and chiB components. */
    gsl_vector *coefs;                      /**< coefficient vector of length n_coefs */
    gsl_vector_long *componentIndices;      /**< Each fit coefficient applies to a single component
                                                 of the vector; this gives the component indices. */
    int n_coefs;                            /**< Number of coefficients in the fit */
    int vec_dim;                            /**< Dimension of the vector */
    FitData **fit_data;            /**< Vector of FitData */
} VectorFitData;

/**
 * Data for a single dynamics node
 */
typedef struct tagDynamicsNodeFitData {
    FitData *omega_data;            /**< A fit to the orbital angular frequency */
    VectorFitData *omega_copr_data; /**< A 2d vector fit for the x and y components of
                                         Omega^{coorb}(t) in arxiv 1705.07089 */
    VectorFitData *chiA_dot_data;   /**< A 3d vector fit for the coorbital components of the
                                         time derivative of chiA taken in the coprecessing frame */
    VectorFitData *chiB_dot_data;   /**< A 3d vector fit for the coorbital components of the
                                         time derivative of chiB taken in the coprecessing frame */
} DynamicsNodeFitData;

/**
 * Data for a single waveform data piece.
 * Waveform data pieces are real time-dependent components of the waveform in the coorbital frame.
 * See the beginning of section IV of https://arxiv.org/abs/1705.07089
 */
typedef struct tagWaveformDataPiece {
    int n_nodes;                                /**< Number of empirical nodes */
    FitData **fit_data;                         /**< FitData at each empirical node */
    gsl_matrix *empirical_interpolant_basis;    /**< The empirical interpolation matrix */
    gsl_vector_long *empirical_node_indices;    /**< The empirical node indices */
} WaveformDataPiece;

/**
 * All WaveformDataPieces needed to evaluate all modes with a fixed value of ell.
 * For m=0 modes we model the real and imaginary parts separately.
 * For m>0, we model the real and imaginary parts of
 * X_pm^{ell, m} = frac{1}{2}( h^{ell, m} pm h^{ell, -m} )
 * where this h is in the coorbital frame.
 * (see https://arxiv.org/abs/1705.07089)
 */
typedef struct tagWaveformFixedEllModeData {
    int ell;                                /**< The fixed value of ell */
    WaveformDataPiece *m0_real_data;        /**< The real (ell, 0) mode data piece */
    WaveformDataPiece *m0_imag_data;        /**< The imag (ell, 0) mode data piece */
    WaveformDataPiece **X_real_plus_data;   /**< One Re[X_+] for each 1 <= m <= ell */
    WaveformDataPiece **X_real_minus_data;  /**< One Re[X_-] for each 1 <= m <= ell */
    WaveformDataPiece **X_imag_plus_data;   /**< One Im[X_+] for each 1 <= m <= ell */
    WaveformDataPiece **X_imag_minus_data;  /**< One Im[X_-] for each 1 <= m <= ell */
} WaveformFixedEllModeData;

/**
 * All data needed by the full surrogate model
 */
typedef struct tagPrecessingNRSurData {
    UINT4 setup; /**< Indicates if this has been initialized */
    int LMax; /**< Maximum ell mode that will ever be evaluated */
    gsl_vector *t_ds; /** Vector of the dynamics surrogate node times, not including half times */
    gsl_vector *t_ds_half_times; /**< t_1/2, t_3/2, and t_5/2 used to start up integration*/
    gsl_vector *t_coorb; /**< Vector of the coorbital surrogate output times. */
    DynamicsNodeFitData **ds_node_data; /** A DynamicsNodeFitData for each time in t_ds.*/
    DynamicsNodeFitData **ds_half_node_data; /** A DynamicsNodeFitData for each time in t_ds_half_times. */
    WaveformFixedEllModeData **coorbital_mode_data; /** One for each 2 <= ell <= LMax */
    UINT4 PrecessingNRSurVersion;   /**< 0 for NRSur7dq2, 1 for NRSur7dq4 */
} PrecessingNRSurData;


/***********************************************************************************/
/****************************** Function declarations*******************************/
/***********************************************************************************/
static void NRSur7dq2_Init_LALDATA(void);
static void NRSur7dq4_Init_LALDATA(void);
static int PrecessingNRSur_Init(PrecessingNRSurData *data, LALH5File *file, UINT4 PrecessingNRSurVersion);
static void PrecessingNRSur_LoadFitData(FitData **fit_data, LALH5File *sub, const char *name);
static void NRSur7dq4_LoadVectorFitData(VectorFitData **vector_fit_data, LALH5File *sub, const char *name, const size_t size);
static void PrecessingNRSur_LoadDynamicsNode(DynamicsNodeFitData **ds_node_data, LALH5File *sub, int i, UINT4 PrecessingNRSurVersion);
static void PrecessingNRSur_LoadCoorbitalEllModes(WaveformFixedEllModeData **coorbital_mode_data, LALH5File *file, int i);
static void PrecessingNRSur_LoadWaveformDataPiece(LALH5File *sub, WaveformDataPiece **data, bool invert_sign);
static bool NRSur7dq2_IsSetup(void);
static bool NRSur7dq4_IsSetup(void);
static double ipow(double base, int exponent); // integer powers

static double NRSur7dq2_eval_fit(FitData *data, double *x);

static void NRSur7dq2_eval_vector_fit(
    double *res, // Result
    VectorFitData *data, // Data for fit
    double *x // size 7, giving mass ratio q, and dimensionless spin components
);

static int NRSur7dq4_effective_spins(REAL8 *chiHat, REAL8 *chi_a,
        const double q, const double chi1z, const double chi2z);
static double NRSur7dq4_eval_fit(FitData *data, double *x);

static void NRSur7dq4_eval_vector_fit(
    double *res, // Result
    VectorFitData *data, // Data for fit
    double *x // size 7, giving mass ratio q, and dimensionless spin components
);

double PrecessingNRSur_eval_fit(FitData *data, double *x, PrecessingNRSurData *__sur_data);

static void PrecessingNRSur_eval_vector_fit(double *res, VectorFitData *data, double *x, PrecessingNRSurData *__sur_data);

static void PrecessingNRSur_normalize_y(
    double chiANorm,
    double chiBNorm,
    double *y
); // Helper to keep spin magnitude constant

static void PrecessingNRSur_normalize_results(
    double normA,
    double normB,
    gsl_vector **quat,
    gsl_vector **chiA,
    gsl_vector **chiB
);

static void PrecessingNRSur_rotate_spins(gsl_vector **chiA, gsl_vector **chiB, gsl_vector *phi);

static void PrecessingNRSur_ds_fit_x(
    double *x, //Result, length 7
    double q, // Mass ratio
    double *y // [q0, qx, qy, qz, phi, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
);

static void PrecessingNRSur_assemble_dydt(
    double *dydt,           // Result, length 11
    double *y,              // ODE solution at the current time step
    double *Omega_coorb_xy, // a form of time derivative of the coprecessing frame
    double omega,           // orbital angular frequency in the coprecessing frame
    double *chiA_dot,       // chiA time derivative
    double *chiB_dot        // chiB time derivative
);


static double cubic_interp(double xout, double *x, double *y);
static gsl_vector *spline_array_interp(gsl_vector *xout, gsl_vector *x, gsl_vector *y);

static double PrecessingNRSur_get_omega(size_t node_index, double q, double *y0, PrecessingNRSurData *__sur_data);

static REAL8 PrecessingNRSur_get_t_ref(
    REAL8 omega_ref,
    REAL8 q,
    REAL8 *chiA0,
    REAL8 *chiB0,
    REAL8 *init_quat,
    REAL8 init_orbphase,
    PrecessingNRSurData *__sur_data
);

static void PrecessingNRSur_get_time_deriv_from_index(
    double *dydt,       // Output: dy/dt evaluated at the ODE time node with index i0. Must have space for 11 entries.
    int i0,             // Time node index. i0=-1, -2, and -3 are used for time nodes 1/2, 3/2, and 5/2 respectively.
    double q,           // Mass ratio
    double *y,           // Current ODE state: [q0, qx, qy, qz, orbphase, chiAx, chiAy, chiAz, chiBx, chiBy, chiBz]
    PrecessingNRSurData *__sur_data
);

static void PrecessingNRSur_get_time_deriv(
    double *dtdy,       // Output: dy/dt evaluated at time t. Must have space for 11 entries.
    double t,           // Time at which the ODE should be evaluated.
    double q,           // Mass ratio
    double *y,           // Current ODE state
    PrecessingNRSurData *__sur_data
);

static int PrecessingNRSur_initialize_at_dynamics_node(
    double *dynamics_data,  // ODE output
    double t_ref,           // reference time. t_ds[i0] will be close to t_ref.
    double q,               // mass ratio
    double *chiA0,          // chiA at t_ref.
    double *chiB0,          // chiB at t_ref.
    double init_orbphase,         // orbital phase at t_ref.
    double *init_quat,          // quaternion at t_ref.
    double normA,           // |chiA|
    double normB,            // |chiB|
    PrecessingNRSurData *__sur_data
);

static void PrecessingNRSur_initialize_RK4_with_half_nodes(
    double *dynamics_data,  // ODE output
    double *time_steps,     // Output: first three time steps. Should have size 3.
    double *dydt0,             // Output: dydt at node 0. Should have size 11.
    double *dydt1,             // Output: dydt at node 1. Should have size 11.
    double *dydt2,             // Output: dydt at node 2. Should have size 11.
    double *dydt3,             // Output: dydt at node 3. Should have size 11.
    double normA,           // |chiA|
    double normB,           // |chiB|
    double q,                // mass ratio
    PrecessingNRSurData *__sur_data
);

static int PrecessingNRSur_initialize_RK4(
    double *dynamics_data,  // ODE output
    double *time_steps,     // Output: first three time steps. Should have size 3.
    double *dydt0,             // Output: dydt at node i0 + 0. Should have size 11.
    double *dydt1,             // Output: dydt at node i0 + 1. Should have size 11.
    double *dydt2,             // Output: dydt at node i0 + 2. Should have size 11.
    double *dydt3,             // Output: dydt at node i0 + 3. Should have size 11.
    double normA,           // |chiA|
    double normB,           // |chiB|
    double q,               // mass ratio
    int i0,                  // the node that is already initialized
    PrecessingNRSurData *__sur_data
);

static void PrecessingNRSur_integrate_AB4(
    double *dynamics_data,  // ODE output
    double *time_steps,     // The first three time steps beginning at i_start.
    double *dydt0,          // dydt at node i_start
    double *dydt1,          // dydt at node i_start + 1
    double *dydt2,          // dydt at node i_start + 2
    double *dydt3,          // dydt at node i_start + 3
    double normA,           // |chiA|
    double normB,           // |chiB|
    double q,               // mass ratio
    int i_start,             // nodes i_start through i_start+3 are already initialized
    PrecessingNRSurData *__sur_data
);

static int PrecessingNRSur_IntegrateDynamics(
    double *dynamics_data,
    double q,
    double *chiA0,
    double *chiB0,
    double omega_ref,
    double init_orbphase,
    double *init_quat,
    LALDict* LALparams,
    UINT4 PrecessingNRSurVersion
);

static void PrecessingNRSur_eval_data_piece(
    gsl_vector *result,
    double q,
    gsl_vector **chiA,
    gsl_vector **chiB,
    WaveformDataPiece *data,
    PrecessingNRSurData *__sur_data
);

static PrecessingNRSurData* PrecessingNRSur_LoadData(Approximant approximant);

static PrecessingNRSurData* PrecessingNRSur_core(
    MultiModalWaveform **h,
    double q,
    double *chiA0,
    double *chiB0,
    double omega_ref,
    double init_orbphase,
    double *init_quat,
    LALValue* ModeArray,
    LALDict* LALparams,
    Approximant approximant
);

static double PrecessingNRSur_StartFrequency(
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 s1x,                      /**< initial value of S1x */
        REAL8 s1y,                      /**< initial value of S1y */
        REAL8 s1z,                      /**< initial value of S1z */
        REAL8 s2x,                      /**< initial value of S2x */
        REAL8 s2y,                      /**< initial value of S2y */
        REAL8 s2z,                      /**< initial value of S2z */
        PrecessingNRSurData *__sur_data    /**< 0 for NRSur7dq2, 1 for NRSur7dq4 */
);


static bool PrecessingNRSur_switch_labels_if_needed(
        REAL8 *m1,                       /**< mass of companion 1 (kg) */
        REAL8 *m2,                       /**< mass of companion 2 (kg) */
        REAL8 *s1x,                      /**< initial value of S1x */
        REAL8 *s1y,                      /**< initial value of S1y */
        REAL8 *s1z,                      /**< initial value of S1z */
        REAL8 *s2x,                      /**< initial value of S2x */
        REAL8 *s2y,                      /**< initial value of S2y */
        REAL8 *s2z                       /**< initial value of S2z */
);
