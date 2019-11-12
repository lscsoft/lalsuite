/*  Copyright (C) 2018 Jonathan Blackman
 *  Utility functions that are useful for evaluating NR surrogates.
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

#include <math.h>
#include <lal/XLALError.h>
#include <lal/Units.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Structure for a multi-modal waveform. Can store any set of modes.
 */
typedef struct tagMultiModalWaveform {
    REAL8 phiRef;                   /**< orbital phase at reference pt. */
    int n_modes;                    /**< Number of modes */
    int n_times;                    /**< Number of time samples in each mode */
    int *lvals;                     /**< The ell values of each mode (n_modes entries) */
    int *mvals;                     /**< The m values of each mode (n_modes entries) */
    gsl_vector **modes_real_part;   /**< The real part of each mode - n_modes vectors of length n_times */
    gsl_vector **modes_imag_part;   /**< The imaginary part of each mode - n_modes vectors of length n_times */
} MultiModalWaveform;

/**
 * Helper structure for computing WignerDMatrices, which require x(t)^n for many values of n.
 */
typedef struct tagRealPowers {
    int LMax;               /**< Maximum L; powers computed depend on LMax. */
    int n_entries;          /**< Number of entries in **powers. */
    int n_times;            /**< Length of each entry. */
    gsl_vector **powers;    /**< The data; x^n. */
} RealPowers;

/**
 * Helper structure for computing WignerDMatrices, which require z(t)^n for many values of n.
 */
typedef struct tagComplexPowers {
    int LMax;               /**< Maximum L; powers computed depend on LMax. */
    int n_entries;          /**< Number of entries in **powers. */
    int n_times;            /**< Length of each entry. */
    gsl_vector **real_part; /**< The real part of z^n. */
    gsl_vector **imag_part; /**< The imag part of z^n. */
} ComplexPowers;

typedef struct tagWignerDMatrices {
    int LMax;                       /**< Includes (ell, m) modes with 2 <= ell <= LMax. */
    int n_entries;                  /**< Number of (ell, m, m') entries. */
    int n_times;                    /**< Number of time samples per entry. */
    gsl_vector **real_part;         /**< The real part of each entry. */
    gsl_vector **imag_part;         /**< The imaginary part of each entry. */
} WignerDMatrices;

UNUSED static void MultiModalWaveform_Init(MultiModalWaveform **wave, int LMax, int n_times);
UNUSED static void MultiModalWaveform_Destroy(MultiModalWaveform *wave);
UNUSED static void WignerDMatrices_Init(WignerDMatrices **matrices, int n_times, int LMax);
UNUSED static void WignerDMatrices_Destroy(WignerDMatrices *matrices);
UNUSED static int WignerDMatrix_Index(int ell, int m, int mp);
UNUSED static void WignerDMatrices_Compute(WignerDMatrices *matrices, gsl_vector **quat);
UNUSED static void RealPowers_Init(RealPowers **rp, int LMax, int n_times);
UNUSED static void RealPowers_Destroy(RealPowers *rp);
UNUSED static void RealPowers_Compute(RealPowers *rp, gsl_vector *x);
UNUSED static void ComplexPowers_Init(ComplexPowers **cp, int LMax, int n_times);
UNUSED static void ComplexPowers_Destroy(ComplexPowers *cp);
UNUSED static void ComplexPowers_Compute(ComplexPowers *cp, gsl_vector *x, gsl_vector *y);

UNUSED static double factorial(int n);
UNUSED static double factorial_ratio(int n, int k);
UNUSED static double binomial(int n, int k);
UNUSED static double wigner_coef(int ell, int mp, int m);

UNUSED static void complex_vector_mult(gsl_vector *x1, gsl_vector *y1, gsl_vector *x2, gsl_vector *y2,
                                       gsl_vector *tmp1, gsl_vector *tmp2);

UNUSED static void NRSur_ab4_dy(
    double *dy,     // Result, length 11
    double *k1,     // dy/dt evaluated at node 1
    double *k2,     // dy/dt evaluated at node 2
    double *k3,     // dy/dt evaluated at node 3
    double *k4,     // dy/dt evaluated at node 4
    double dt12,    // t_2 - t_1
    double dt23,    // t_3 - t_2
    double dt34,    // t_4 - t_3
    double dt45,    // t_5 - t_4
    int dim         // Vector dimension of dy and the k vectors.
);

/**
 * Function to transform waveform modes from the coorbital frame
 * to the inertial frame. First transforms to the coprecessing frame
 * by rotating about the z-axis by the orbital phase, then rotates
 * to the inertial frame using the coprecessing frame quaternions.
 */
UNUSED static void TransformModesCoorbitalToInertial(
    MultiModalWaveform *h,          /**< Output. Dimensionless waveform modes. Should be initialized already. */
    MultiModalWaveform *h_coorb,    /**< Coorbital frame waveform modes. */
    gsl_vector **quat,              /**< Coprecessing frame quaternions. */
    gsl_vector *orbphase            /**< Orbital phase. */
);

UNUSED static int quatInv(double *qInv, double *q);
UNUSED static int multiplyQuats(double *prod, double *q1, double *q2);
UNUSED static int quaternionTransformVector(
        double *new_vec, double *quat, double *vec);
UNUSED static int transformTimeDependentVector(gsl_vector **vec, gsl_vector **quat);

UNUSED static int get_dimless_omega(
        REAL8 *omegaMin_dimless,
        REAL8 *omegaRef_dimless,
        const REAL8 fMin,
        const REAL8 fRef,
        const REAL8 Mtot_sec
);
