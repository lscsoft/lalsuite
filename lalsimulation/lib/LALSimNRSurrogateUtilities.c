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

#include "LALSimNRSurrogateUtilities.h"
#include <lal/LALConstants.h>


/****************** Init and Destroy functions ****************/

/**
 * Initialize a MultiModalWaveforme.  We also have the usual frame data, which is initially the
m structure
 */
static void MultiModalWaveform_Init(
    MultiModalWaveform **wave, /**< Output; the MultiModalWaveform. */
    int LMax,                  /**< (ell, m) modes with 2 <= ell <= LMax will be created. */
    int n_times                /**< Number of time samples per mode. */
) {
    if (!wave) exit(1);
    if (LMax < 2) XLAL_ERROR_VOID(XLAL_FAILURE, "Got LMax=%d < 2!\n", LMax);
    if (*wave) MultiModalWaveform_Destroy(*wave);
    (*wave) = XLALCalloc(1, sizeof(MultiModalWaveform));

    int n_modes = LMax*(LMax+2) - 3;

    (*wave)->n_modes = n_modes;
    (*wave)->lvals = XLALCalloc(n_modes, sizeof(int));
    (*wave)->mvals = XLALCalloc(n_modes, sizeof(int));
    (*wave)->n_times = n_times;
    gsl_vector **modes_real_part = XLALMalloc(n_modes * sizeof(*modes_real_part));
    gsl_vector **modes_imag_part = XLALMalloc(n_modes * sizeof(*modes_imag_part));
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

/**
 * Destroy a MultiModalWaveform structure; free all associated memory.
 */
static void MultiModalWaveform_Destroy(MultiModalWaveform *wave) {
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

/**
 * Initialize a ComplexPowers structure
 */
static void ComplexPowers_Init(
    ComplexPowers **cp, /**< The ComplexPowers structure to be initialized.*/
    int LMax,           /**< The maximum ell that will be needed by WignerDMatrices.*/
    int n_times         /**< The number of time samples.*/
) {
    // include powers from -2*LMax to 2*LMax inclusive
    if (!cp) exit(1);
    if (*cp) ComplexPowers_Destroy(*cp);
    (*cp) = XLALCalloc(1, sizeof(ComplexPowers));

    int n_entries = 4*LMax+1;
    (*cp)->LMax = LMax;
    (*cp)->n_entries = n_entries;

    gsl_vector **real_part = XLALMalloc(n_entries * sizeof(*real_part));
    gsl_vector **imag_part = XLALMalloc(n_entries * sizeof(*imag_part));
    (*cp)->real_part = real_part;
    (*cp)->imag_part = imag_part;

    for (int i=0; i<n_entries; i++) {
        (*cp)->real_part[i] = gsl_vector_calloc(n_times);
        (*cp)->imag_part[i] = gsl_vector_calloc(n_times);
    }
}

/**
 * Destroy a ComplexPowers structure; free all associated memory.
 */
static void ComplexPowers_Destroy(ComplexPowers *cp) {
    if (!cp) return;
    for (int i=0; i<cp->n_entries; i++) {
        if (cp->real_part[i]) gsl_vector_free(cp->real_part[i]);
        if (cp->imag_part[i]) gsl_vector_free(cp->imag_part[i]);
    }
    free(cp->real_part);
    free(cp->imag_part);
    XLALFree(cp);
}

/**
 * Initialize a RealPowers structure
 */
static void RealPowers_Init(
    RealPowers **rp, /**< The RealPowers structure to be initialized.*/
    int LMax,           /**< The maximum ell that will be needed by WignerDMatrices.*/
    int n_times         /**< The number of time samples.*/
) {
    // include powers from 0 to 2*LMax inclusive
    if (!rp) exit(1);
    if (*rp) RealPowers_Destroy(*rp);
    (*rp) = XLALCalloc(1, sizeof(RealPowers));

    int n_entries = 2*LMax+1;
    (*rp)->LMax = LMax;
    (*rp)->n_entries = n_entries;

    gsl_vector **powers = XLALMalloc(n_entries * sizeof(*powers));
    (*rp)->powers = powers;

    for (int i=0; i<n_entries; i++) {
        (*rp)->powers[i] = gsl_vector_calloc(n_times);
    }
}

/**
 * Destroy a RealPowers structure; free all associated memory.
 */
static void RealPowers_Destroy(RealPowers *rp) {
    if (!rp) return;
    for (int i=0; i<rp->n_entries; i++) {
        if (rp->powers[i]) gsl_vector_free(rp->powers[i]);
    }
    free(rp->powers);
    XLALFree(rp);
}

/**
 * Initialize a WignerDMatrices structure
 */
static void WignerDMatrices_Init(
    WignerDMatrices **matrices, /**< Output; the WignerDMatrices.*/
    int n_times,                /**< The number of time samples for each (ell, m, m') component.*/
    int LMax                    /**< The maximum ell to generate. */
) {
    if (!matrices) exit(1);
    if (LMax < 2) XLAL_ERROR_VOID(XLAL_FAILURE, "Got LMax=%d < 2!\n", LMax);
    if (*matrices) WignerDMatrices_Destroy(*matrices);
    (*matrices) = XLALCalloc(1, sizeof(WignerDMatrices));

    int n_entries = 0;
    for (int ell=2; ell<=LMax; ell++) {
        n_entries += (2*ell+1) * (2*ell+1);
    }

    (*matrices)->LMax = LMax;
    (*matrices)->n_entries = n_entries;
    (*matrices)->n_times = n_times;

    gsl_vector **real_part = XLALMalloc(n_entries * sizeof(*real_part));
    gsl_vector **imag_part = XLALMalloc(n_entries * sizeof(*imag_part));
    (*matrices)->real_part = real_part;
    (*matrices)->imag_part = imag_part;

    for (int i=0; i<n_entries; i++) {
        (*matrices)->real_part[i] = gsl_vector_calloc(n_times);
        (*matrices)->imag_part[i] = gsl_vector_calloc(n_times);
    }
}

/**
 * Destroy a WignerDMatrices structure; free all associated memory.
 */
static void WignerDMatrices_Destroy(WignerDMatrices *matrices) {
    if (!matrices) return;
    for (int i=0; i<matrices->n_entries; i++) {
        if (matrices->real_part[i]) gsl_vector_free(matrices->real_part[i]);
        if (matrices->imag_part[i]) gsl_vector_free(matrices->imag_part[i]);
    }
    free(matrices->real_part);
    free(matrices->imag_part);
    XLALFree(matrices);
}



/**************** Utility functions **************************/

/**
 * Computes the index of the (ell, m, mp) component in a WignerDMatices structure
 */
static int WignerDMatrix_Index(int ell, int m, int mp) {
    // Start of the (m, mp) matrix, which has size (2*ell + 1) X (2*ell + 1)
    int i0 = (ell*(ell*ell*4 - 1))/3 - 10;
    int res = i0 + (2*ell+1)*(ell+m) + (ell+mp);
    return res;
}

static REAL8 factorial(int n) {
    if (n <= 0) return 1.0;
    return factorial(n-1) * n;
}

static REAL8 factorial_ratio(int n, int k) {
    if (n <= k) return 1.0;
    return factorial_ratio(n-1, k) * n;
}

static REAL8 binomial(int n, int k) {
    return factorial_ratio(n, k) / factorial(n-k);
}

static REAL8 wigner_coef(int ell, int mp, int m) {
    return sqrt(factorial(ell+m) * factorial(ell-m) / (factorial(ell+mp) * factorial(ell-mp)));
}

/**
 * Multiplies z1(t) * z2(t) = (x1(t) + i*y1(t)) * (x2(t) + i*y2(t)),
 * storing the result in x1 and y1.
 */
static void complex_vector_mult(
    gsl_vector *x1,     /**< Real part of z1. */
    gsl_vector *y1,     /**< Imag part of z1. */
    gsl_vector *x2,     /**< Real part of z2. */
    gsl_vector *y2,     /**< Imag part of z2. */
    gsl_vector *tmpx,   /**< A vector for storing temporary results. */
    gsl_vector *tmpy    /**< A second vector for storing temporary results. */
) {
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


/*************** Heavy lifting functions ********************/
/*
 * Given dy/dt evaluated as a function of y and t at 4 time nodes, compute
 * dy = y_5 - y_4 using a 4th-order accurate Adams-Bashforth scheme.
 * See https://dcc.ligo.org/T1800069 for details.
 */
static void NRSur_ab4_dy(
    REAL8 *dy,     /**< Result */
    REAL8 *k1,     /**< dy/dt evaluated at node 1 */
    REAL8 *k2,     /**< dy/dt evaluated at node 2 */
    REAL8 *k3,     /**< dy/dt evaluated at node 3 */
    REAL8 *k4,     /**< dy/dt evaluated at node 4 */
    REAL8 t12,     /**< t_2 - t_1 */
    REAL8 t23,     /**< t_3 - t_2 */
    REAL8 t34,     /**< t_4 - t_3 */
    REAL8 t45,     /**< t_5 - t_4 */
    int dim         /**< Vector dimension. Length of dy, k1, k2, k3, k4. */
) {

    // These are defined in https://dcc.ligo.org/T1800069, except we rename Dj to denomj.
    REAL8 t13, t14, t24, denom1, denom2, denom3;
    REAL8 A1, A2, A3, A4, B1, B2, B3, B4, C1, C2, C3, C4; // Matrix coefficients in the M matrix.
    REAL8 dx_vec[4]; // Integral of the p vector from t4 to t5; [t45, t45^2/2, t45^3/3, t45^4/4]

    // Temp variable for the coefficients multiplying 1, x, x^2, and x^3 (k * M in the DCC document)
    REAL8 tmp_coef;

    int i;

    // Compute time intervals
    t13 = t12 + t23;
    t14 = t13 + t34;
    t24 = t23 + t34;

    // Compute denomenators
    denom1 = t12 * t13 * t14;
    denom2 = t12 * t23 * t24;
    denom3 = t13 * t23 * t34;

    // Compute dx_vec
    dx_vec[0] = t45;
    dx_vec[1] = t45 * t45 / 2.0;
    dx_vec[2] = t45 * t45 * t45 / 3.0;
    dx_vec[3] = dx_vec[1] * dx_vec[1];

    // Column 1, which goes with dx_vec[1]:
    A1 = -1.0 * t24 * t34 / denom1;
    A2 = t14 * t34 / denom2;
    A3 = -1.0 * t14 * t24 / denom3;
    A4 = -1 * (A1 + A2 + A3);

    // Column 2, which goes with dx_vec[2]:
    B1 = -1.0 * (t24 + t34) / denom1;
    B2 = (t14 + t34) / denom2;
    B3 = -1.0 * (t14 + t24) / denom3;
    B4 = -1 * (B1 + B2 + B3);

    // Column 3, which goes with dx_vec[3]:
    C1 = -1.0 / denom1;
    C2 = 1.0 / denom2;
    C3 = -1.0 / denom3;
    C4 = -1 * (C1 + C2 + C3);

    // Compute results
    for (i=0; i<dim; i++) {
        tmp_coef = k4[i] * dx_vec[0];
        tmp_coef += (A1 * k1[i] + A2 * k2[i] + A3 * k3[i] + A4 * k4[i]) * dx_vec[1];
        tmp_coef += (B1 * k1[i] + B2 * k2[i] + B3 * k3[i] + B4 * k4[i]) * dx_vec[2];
        tmp_coef += (C1 * k1[i] + C2 * k2[i] + C3 * k3[i] + C4 * k4[i]) * dx_vec[3];
        dy[i] = tmp_coef;
    }
}



/**
 * Computes all powers of z(t) = x(t) + i*y(t) needed for WignerDMatrices.
 */
static void ComplexPowers_Compute(ComplexPowers *cp, gsl_vector *x, gsl_vector *y) {
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

/**
 * Computes all powers of x(t) needed for WignerDMatrices.
 */
static void RealPowers_Compute(RealPowers *rp, gsl_vector *x) {
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


/**
 * Transforms modes from the coorbital frame to the inertial frame.
 * Note that this is a coorbital frame determined entirely from the waveform
 * modes and *not* from the trajectories; it is a rotation of the coprecessing
 * frame about the z-axis by an orbital phase computed from the waveform modes.
 * See eqns 6-9 and the surrounding text of https://arxiv.org/abs/1705.07089
 */
static void TransformModesCoorbitalToInertial(
    MultiModalWaveform *h,          /**< Output. Dimensionless waveform modes.
                                         Should be initialized already. */
    MultiModalWaveform *h_coorb,    /**< Coorbital frame waveform modes. */
    gsl_vector **quat,              /**< Coprecessing frame quaternions - 4 vectors. */
    gsl_vector *orbphase            /**< Orbital phase. */
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
    MultiModalWaveform_Init(&h_copr, lmax, n_times);
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
                // compute and add combinations of {real, imag} and {cosmphi, sinmphi}
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
    WignerDMatrices_Init(&matrices, n_times, lmax);
    WignerDMatrices_Compute(matrices, quat);
    int matrix_index;
    for (ell = 2; ell <= lmax; ell++) {
        for (m= -1 * ell; m <= ell; m++) {
            i = ell*(ell+1) - 4 + m;
            for (mp = -1 * ell; mp <= ell; mp++) {
                j = ell*(ell+1) - 4 + mp;
                matrix_index = WignerDMatrix_Index(ell, m, mp);
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
    WignerDMatrices_Destroy(matrices);
    MultiModalWaveform_Destroy(h_copr);
    gsl_vector_free(cosmphi);
    gsl_vector_free(sinmphi);
    gsl_vector_free(tmp_vec);
}

/**
 * Computes all the Wigner D Matrices.
 * We compute them all at once because a lot of the work is just processing quat.
 * Parts of this function are adapted from GWFrames:
 * https://github.com/moble/GWFrames
 * written by Michael Boyle, based on his paper: http://arxiv.org/abs/1302.2919
 * although we are working with the conjugate of quat.
 */
static void WignerDMatrices_Compute(
    WignerDMatrices *matrices, /**< The initialized WignerDMatrices which will contain the output.*/
    gsl_vector **quat          /**< The 4 quaternion components representing the coorbital frame.*/
) {

    int n_times = quat[0]->size;
    int i, j, ell, m, mp, rho_min, rho_max, rho;
    REAL8 tmp_re, tmp_im, tmp_abs_sqr, coef, c;
    REAL8 eps_sqr = 1.0e-24;

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

    ComplexPowers_Init(&ra_powers, matrices->LMax, n_times);
    ComplexPowers_Init(&rb_powers, matrices->LMax, n_times);
    RealPowers_Init(&abs_ra_sqr_powers, matrices->LMax, n_times);
    RealPowers_Init(&abs_rb_over_ra_sqr_powers, matrices->LMax, n_times);


    ComplexPowers_Compute(ra_powers, safe_ra_real, safe_ra_imag);
    ComplexPowers_Compute(rb_powers, safe_rb_real, safe_rb_imag);
    RealPowers_Compute(abs_ra_sqr_powers, safe_ra_mag_sqr);
    RealPowers_Compute(abs_rb_over_ra_sqr_powers, safe_abs_rb_over_ra_sqr);


    // Compute the matrices. We don't use safe_* anymore, so we can use them as temporary results
    for (ell=2; ell <= matrices->LMax; ell++) {
        for (m = -1*ell; m <= ell; m++) {
            for (mp = -1*ell; mp <= ell; mp++) {
                i = WignerDMatrix_Index(ell, m, mp);
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

    REAL8 zx, zy;
    int k;
    // Now fix the values at bad indices
    for (j=0; j<n_times; j++) {
        rho = gsl_vector_int_get(index_types, j);
        if (rho != 0) {
            for (ell=2; ell <= matrices->LMax; ell++) {
                for (m=-1*ell; m <= ell; m++) {
                    for (mp = -1*ell; mp <= ell; mp++) {
                        i = WignerDMatrix_Index(ell, m, mp);
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
    ComplexPowers_Destroy(ra_powers);
    ComplexPowers_Destroy(rb_powers);
    RealPowers_Destroy(abs_ra_sqr_powers);
    RealPowers_Destroy(abs_rb_over_ra_sqr_powers);
}


/**
 * Computes inverse, qInv, of a quaternion q such that q*qInv = 1.
 */
UNUSED static int quatInv(
    REAL8 *qInv,   /**< Output: Inverse of q. Should be initialized with size=4 */
    REAL8 *q       /**< Input quaternion. Should have size=4*/
) {
    const REAL8 normSqr =  q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    qInv[0] = q[0]/normSqr;
    for (int i=1; i<4; i++) {
        qInv[i] = -q[i]/normSqr;
    }

    return XLAL_SUCCESS;
}

/**
 * Multiplies two quaternions.
 */
UNUSED static int multiplyQuats(
    REAL8 *prod,   /**< Output: Product of q1 and q2. Should be initialized with size=4. */
    REAL8 *q1,     /**< First quaternion. Should have size=4. */
    REAL8 *q2      /**< Second quaternion. Should have size=4. */
) {
    prod[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
    prod[1] = q1[2]*q2[3] - q2[2]*q1[3] + q1[0]*q2[1] + q2[0]*q1[1];
    prod[2] = q1[3]*q2[1] - q2[3]*q1[1] + q1[0]*q2[2] + q2[0]*q1[2];
    prod[3] = q1[1]*q2[2] - q2[1]*q1[2] + q1[0]*q2[3] + q2[0]*q1[3];
    return XLAL_SUCCESS;
}


/**
 * Given a coprecessing frame quaternion (quat), transforms a 3-vector (vec)
 * from the coprecessing frame to the inertial frame.
 * Reference: Eq.(A6) of https://arxiv.org/pdf/1302.2919.pdf
 */
UNUSED static int quaternionTransformVector(
    REAL8 *new_vec,  /**< Output: Transformed vector in inertial frame. Should be initialized with size=3 */
    REAL8 *quat,     /**< Coprecessing frame quaternion. Should have size=4.*/
    REAL8 *vec       /**< Input vector in coprecessing frame. Should have size=3.*/
) {

    // quaternion representation of vec, such that the scalar component is zero
    REAL8 vec_quat[4] = {0, vec[0], vec[1], vec[2]};

    // q_prod = vec_quat * inv(quat)
    REAL8 q_prod[4];
    REAL8 qInv[4];
    quatInv(qInv, quat);
    multiplyQuats(q_prod, vec_quat, qInv);

    // quaternion representation of output vec, such that the
    // scalar component is zero
    REAL8 new_vec_quat[4];

    // new_vec_quat = quat * vec_quat * inv(quat)
    multiplyQuats(new_vec_quat, quat, q_prod);
    new_vec[0] = new_vec_quat[1];
    new_vec[1] = new_vec_quat[2];
    new_vec[2] = new_vec_quat[3];

    return XLAL_SUCCESS;
}

/**
 * Transforms a vector from the coprecessing frame to the inertial frame
 * using the coprecessing frame quaternions.
 */
UNUSED static int transformTimeDependentVector(
    gsl_vector **vec,  /**< Input and Output: 3 time-dependent components of vec in the coprecessing frame */
    gsl_vector **quat  /**< 4 time-dependent components of coprecessing frame quaternions. */
) {

    REAL8 tmp_vec[3];
    REAL8 tmp_new_vec[3];
    REAL8 tmp_quat[4];

    REAL8 *x = vec[0]->data;
    REAL8 *y = vec[1]->data;
    REAL8 *z = vec[2]->data;

    REAL8 *q0 = quat[0]->data;
    REAL8 *q1 = quat[1]->data;
    REAL8 *q2 = quat[2]->data;
    REAL8 *q3 = quat[3]->data;

    int n = quat[0]->size;
    for (int i=0; i<n; i++) {

        tmp_vec[0] = x[i];
        tmp_vec[1] = y[i];
        tmp_vec[2] = z[i];

        tmp_quat[0] = q0[i];
        tmp_quat[1] = q1[i];
        tmp_quat[2] = q2[i];
        tmp_quat[3] = q3[i];

        // update x,y,z, which also updates vec, with the transformed vector
        quaternionTransformVector(tmp_new_vec, tmp_quat, tmp_vec);
        x[i] = tmp_new_vec[0];
        y[i] = tmp_new_vec[1];
        z[i] = tmp_new_vec[2];
    }

    return XLAL_SUCCESS;
}


/**
 * Wrapper to get dimensionless omegaMin/omegaRef from fMin/fRef.
 */
UNUSED static int get_dimless_omega(
        REAL8 *omegaMin_dimless,  /**< Output. omegaMin in units rad/M. */
        REAL8 *omegaRef_dimless,  /**< Output. omegaRef in units rad/M. */
        const REAL8 fMin,         /**< fMin in Hertz. */
        const REAL8 fRef,         /**< fRef in Hertz. */
        const REAL8 Mtot_sec      /**< Total mass in seconds. */
) {

    // Orbital angular frequency = 0.5*(wave angular frequency) = pi*(wave frequency)
    *omegaMin_dimless = (Mtot_sec * fMin) * LAL_PI;

    if (fRef == 0) {
        // If fRef is 0, set it to fMin
        *omegaRef_dimless = *omegaMin_dimless;
    }
    else {
        *omegaRef_dimless = (Mtot_sec * fRef) * LAL_PI;
    }

    if (*omegaRef_dimless + 1e-13 < *omegaMin_dimless){
        XLAL_ERROR(XLAL_EINVAL, "fRef cannot be lesser than fMin.");
    }

    return XLAL_SUCCESS;
}
