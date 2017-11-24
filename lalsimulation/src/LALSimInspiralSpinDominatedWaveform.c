/*
 * Copyright (C) 2014 M. Tapai
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
#include <complex.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>

#define LAL_SDW_ABSOLUTE_TOLERANCE 1.e-12
#define LAL_SDW_RELATIVE_TOLERANCE 1.e-12
#define LAL_SDW_NUM_VARIABLES 3
#define LAL_SDW_MAX_PN_PARAM 0.8

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif



static const REAL8 G_CP2 = LAL_G_SI / LAL_C_SI / LAL_C_SI;

enum {
    PNDEF = -1 ,PN00 = 0, PN05 = 1, PN10 = 2, PN15 =3, PN20 =4, PN25 =5, PN3 =6, PN_ORDER =7,
};

typedef enum {
    PHI, OMEGA, PSI, PHASE = PSI,
} TIME_DEPENDANT;

typedef enum {
    PLUS_ = 0, MINUS = 1, CROSS_ = 1, PLUS_MINUS_DIM = 2, PLUS_CROSS_DIM = 2,
} COMPONENTS;

typedef enum {
    PN00DIM = 2, PN05DIM = 11, PN10DIM = 15, PN15DIM = 17,
} COEFFICIENT_DIMENSIONS;

typedef enum {
    TRIGONOMETRIC_POWER = 5, AMPCOEFF_DIM = 11, OMEGA_POWER_DIM = 6, PHI_PSI_DIM = 6,
} CONSTANTS;

#define vectorProd(lhs, rhs, denominator, result)                    \
    result[X] = (lhs[Y] * rhs[Z] - lhs[Z] * rhs[Y]) / denominator;    \
    result[Y] = (lhs[Z] * rhs[X] - lhs[X] * rhs[Z]) / denominator;    \
    result[Z] = (lhs[X] * rhs[Y] - lhs[Y] * rhs[X]) / denominator;

/**
 * Structure containing the prefered variabloes for Spin-Dominated waveforms.
 */
typedef struct tagLALSDWaveformParams {
    REAL8 totalmass; //total mass of the binary
    REAL8 nu;  // mass ratio, <1
    REAL8 chi1; // chi1 dimensionless spin parameter
    REAL8 dist; // distance to the source
    REAL8 kappa1; // angle between L and S1
    REAL8 beta1;  // angle between J and S1
    REAL8 theta; // angle between J and N
    REAL8 eps; // PN paramter
    REAL8 xi; // second small parameter
    REAL8 omega;
    int pnamp;
    int pnphase;
    REAL8 ccoeff00pn[4];
    REAL8 ccoeff05pn[22];
    REAL8 ccoeff10pn[30];
    REAL8 ccoeff15pn[34];
    REAL8 prevdomega;
    REAL8 polarizationangle;
} LALSDWaveformParams;

static INT4 XLALSpinDominatedWaveformStoppingTest(UNUSED REAL8 t, const REAL8 values[], REAL8 dvalues[],
UNUSED void *mparams);

static INT4 XLALSpinDominatedWaveformDerivatives(UNUSED REAL8 t, const REAL8 values[], REAL8 dvalues[], void *mparams);

int XLALSpinDominatedWaveformBuild(LALSDWaveformParams *params, REAL8 expr[], REAL8TimeSeries **hplus,
        REAL8TimeSeries **hcross, int idx);

int XLALSpinDominatedWaveformConstantCoefficients(LALSDWaveformParams * params);

/**
 * Function for allocating memory for a matrix
 */
static REAL8 *XLALDmatrix(INT8 nrh, INT8 nch) {
    INT8 size = (nrh) * (nch) * sizeof(REAL8);
    REAL8 *ptr = (REAL8 *) LALMalloc(size);
    if (ptr != NULL) {
        return ptr;
    }
    printf("malloc error");
    return NULL;
}

/**
 * Function for freeing memory for a matrix
 */
static void XLALFreeDmatrix(REAL8 *m) {
    LALFree(m);
}

/**
 * Function for calculating the constant coefficients of Spin-Dominated waveforms
 * See tables 1 to 5 in the appendix of Arxiv:1209.1722
 */
int XLALSpinDominatedWaveformConstantCoefficients(LALSDWaveformParams * params) {

    int i, j;
    REAL8 *acoeff00pn, *b0coeff0pn, *d0coeff0pn, *acoeff0_5pn, *b0coeff0_5pn, *d0coeff0_5pn, *acoeff1pn, *b0coeff1pn,
            *d0coeff1pn, *b1coeff1pn, *d1coeff1pn, *acoeff1_5pn, *b0coeff1_5pn, *d0coeff1_5pn, *b1coeff1_5pn,
            *d1coeff1_5pn;

    REAL8 skp[TRIGONOMETRIC_POWER], ckp[TRIGONOMETRIC_POWER], stp[TRIGONOMETRIC_POWER], ctp[TRIGONOMETRIC_POWER];
    skp[0] = ckp[0] = stp[0] = ctp[0] = 1.;
    skp[1] = sin(params->kappa1);
    ckp[1] = cos(params->kappa1);
    stp[1] = sin(params->theta);
    ctp[1] = cos(params->theta);
    for (i = 2; i < TRIGONOMETRIC_POWER; ++i) {
        skp[i] = skp[1] * skp[i - 1];
        ckp[i] = ckp[1] * ckp[i - 1];
        stp[i] = stp[1] * stp[i - 1];
        ctp[i] = ctp[1] * ctp[i - 1];
    }
    REAL8 k[PLUS_MINUS_DIM] = { ckp[1] - 1., ckp[1] + 1. };
    REAL8 c2t = cos(2. * params->theta);


    acoeff00pn = XLALDmatrix(PN00DIM, PLUS_MINUS_DIM);
    b0coeff0pn = XLALDmatrix(PN00DIM, PLUS_MINUS_DIM);
    d0coeff0pn = XLALDmatrix(PN00DIM, PLUS_MINUS_DIM);
    acoeff0_5pn = XLALDmatrix(PN05DIM, PLUS_MINUS_DIM);
    b0coeff0_5pn = XLALDmatrix(PN05DIM, PLUS_MINUS_DIM);
    d0coeff0_5pn = XLALDmatrix(PN05DIM, PLUS_MINUS_DIM);
    acoeff1pn = XLALDmatrix(PN10DIM, PLUS_MINUS_DIM);
    b0coeff1pn = XLALDmatrix(PN10DIM, PLUS_MINUS_DIM);
    d0coeff1pn = XLALDmatrix(PN10DIM, PLUS_MINUS_DIM);
    b1coeff1pn = XLALDmatrix(PN10DIM, PLUS_MINUS_DIM);
    d1coeff1pn = XLALDmatrix(PN10DIM, PLUS_MINUS_DIM);
    acoeff1_5pn = XLALDmatrix(PN15DIM, PLUS_MINUS_DIM);
    b0coeff1_5pn = XLALDmatrix(PN15DIM, PLUS_MINUS_DIM);
    d0coeff1_5pn = XLALDmatrix(PN15DIM, PLUS_MINUS_DIM);
    b1coeff1_5pn = XLALDmatrix(PN15DIM, PLUS_MINUS_DIM);
    d1coeff1_5pn = XLALDmatrix(PN15DIM, PLUS_MINUS_DIM);
    /*
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%   0 PN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     */
    acoeff00pn[0 + PN00DIM * PLUS_] = -2. * k[PLUS_];
    acoeff00pn[0 + PN00DIM * MINUS] = +2. * k[MINUS];
    acoeff00pn[1 + PN00DIM * PLUS_] = -k[PLUS_];
    acoeff00pn[1 + PN00DIM * MINUS] = +k[MINUS];
    b0coeff0pn[0 + PN00DIM * PLUS_] = -1.;
    b0coeff0pn[0 + PN00DIM * MINUS] = -1.;
    b0coeff0pn[1 + PN00DIM * PLUS_] = -2.;
    b0coeff0pn[1 + PN00DIM * MINUS] = -2.;
    d0coeff0pn[0 + PN00DIM * PLUS_] = +0.;
    d0coeff0pn[0 + PN00DIM * MINUS] = +0.;
    d0coeff0pn[1 + PN00DIM * PLUS_] = +0.;
    d0coeff0pn[1 + PN00DIM * MINUS] = +0.;
    for (i = 0; i <= 1; i++) {
        for (j = 0; j <= 1; j++) {
            params->ccoeff00pn[i + PN00DIM * j] = acoeff00pn[i + PN00DIM * j]
                    + skp[2] * (b0coeff0pn[i + PN00DIM * j] + ckp[1] * d0coeff0pn[i + PN00DIM * j]);
        }
    }
    /*
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%   0.5 PN   %%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     */
    acoeff0_5pn[0 + PN05DIM * PLUS_] = +4. * k[PLUS_] * (6. - stp[2]);
    acoeff0_5pn[0 + PN05DIM * MINUS] = +4. * k[MINUS] * (6. - stp[2]);
    acoeff0_5pn[1 + PN05DIM * PLUS_] = +4. * k[PLUS_];
    acoeff0_5pn[1 + PN05DIM * MINUS] = +4. * k[MINUS];
    acoeff0_5pn[2 + PN05DIM * PLUS_] = +2. * k[PLUS_] * (6. - stp[2]);
    acoeff0_5pn[2 + PN05DIM * MINUS] = -2. * k[MINUS] * (6. - stp[2]);
    acoeff0_5pn[3 + PN05DIM * PLUS_] = 12. * k[PLUS_];
    acoeff0_5pn[3 + PN05DIM * MINUS] = 12. * k[MINUS];
    acoeff0_5pn[4 + PN05DIM * PLUS_] = +2. * k[PLUS_] * (2. * stp[2] - 3.);
    acoeff0_5pn[4 + PN05DIM * MINUS] = -2. * k[MINUS] * (2. * stp[2] - 3.);
    acoeff0_5pn[5 + PN05DIM * PLUS_] = -2. * k[PLUS_];
    acoeff0_5pn[5 + PN05DIM * MINUS] = -2. * k[MINUS];
    acoeff0_5pn[6 + PN05DIM * PLUS_] = -2. * params->ccoeff00pn[1 + PN00DIM * PLUS_] * (6. - stp[2]);
    acoeff0_5pn[6 + PN05DIM * MINUS] = +2. * params->ccoeff00pn[1 + PN00DIM * MINUS] * (6. - stp[2]);
    acoeff0_5pn[7 + PN05DIM * PLUS_] = +2. * k[PLUS_];
    acoeff0_5pn[7 + PN05DIM * MINUS] = -2. * k[MINUS];
    acoeff0_5pn[8 + PN05DIM * PLUS_] = (44. - 34. * stp[2] + 2. * (5. * stp[2] - 46.) * ckp[1]);
    acoeff0_5pn[8 + PN05DIM * MINUS] = (44. - 34. * stp[2] - 2. * (5. * stp[2] - 46.) * ckp[1]);
    acoeff0_5pn[9 + PN05DIM * PLUS_] = 22. + 46. * ckp[1];
    acoeff0_5pn[9 + PN05DIM * MINUS] = 22. - 46. * ckp[1];
    acoeff0_5pn[10 + PN05DIM * PLUS_] = -2. * k[PLUS_] * (3. - 2. * stp[2]);
    acoeff0_5pn[10 + PN05DIM * MINUS] = -2. * k[MINUS] * (3. - 2. * stp[2]);
    b0coeff0_5pn[0 + PN05DIM * PLUS_] = +(46. - 5. * stp[2]);
    b0coeff0_5pn[0 + PN05DIM * MINUS] = -(46. - 5. * stp[2]);
    b0coeff0_5pn[1 + PN05DIM * PLUS_] = +3.;
    b0coeff0_5pn[1 + PN05DIM * MINUS] = -3.;
    b0coeff0_5pn[2 + PN05DIM * PLUS_] = (2. - 3. * stp[2]);
    b0coeff0_5pn[2 + PN05DIM * MINUS] = (2. - 3. * stp[2]);
    b0coeff0_5pn[3 + PN05DIM * PLUS_] = +23.;
    b0coeff0_5pn[3 + PN05DIM * MINUS] = -23.;
    b0coeff0_5pn[4 + PN05DIM * PLUS_] = -c2t;
    b0coeff0_5pn[4 + PN05DIM * MINUS] = -c2t;
    b0coeff0_5pn[5 + PN05DIM * PLUS_] = -4.;
    b0coeff0_5pn[5 + PN05DIM * MINUS] = +4.;
    b0coeff0_5pn[6 + PN05DIM * PLUS_] = +0.;
    b0coeff0_5pn[6 + PN05DIM * MINUS] = +0.;
    b0coeff0_5pn[7 + PN05DIM * PLUS_] = +3.;
    b0coeff0_5pn[7 + PN05DIM * MINUS] = +3.;
    b0coeff0_5pn[8 + PN05DIM * PLUS_] = -15. * (2. - 3. * stp[2]);
    b0coeff0_5pn[8 + PN05DIM * MINUS] = -15. * (2. - 3. * stp[2]);
    b0coeff0_5pn[9 + PN05DIM * PLUS_] = -15.;
    b0coeff0_5pn[9 + PN05DIM * MINUS] = -15.;
    b0coeff0_5pn[10 + PN05DIM * PLUS_] = -4. * (3. - 2. * stp[2]);
    b0coeff0_5pn[10 + PN05DIM * MINUS] = +4. * (3. - 2. * stp[2]);
    d0coeff0_5pn[0 + PN05DIM * PLUS_] = +5. * (3. * stp[2] - 2.);
    d0coeff0_5pn[0 + PN05DIM * MINUS] = +5. * (3. * stp[2] - 2.);
    d0coeff0_5pn[1 + PN05DIM * PLUS_] = -1.;
    d0coeff0_5pn[1 + PN05DIM * MINUS] = -1.;
    d0coeff0_5pn[2 + PN05DIM * PLUS_] = +0.;
    d0coeff0_5pn[2 + PN05DIM * MINUS] = +0.;
    d0coeff0_5pn[3 + PN05DIM * PLUS_] = -5.;
    d0coeff0_5pn[3 + PN05DIM * MINUS] = -5.;
    d0coeff0_5pn[4 + PN05DIM * PLUS_] = +0.;
    d0coeff0_5pn[4 + PN05DIM * MINUS] = +0.;
    d0coeff0_5pn[5 + PN05DIM * PLUS_] = +3.;
    d0coeff0_5pn[5 + PN05DIM * MINUS] = +3.;
    d0coeff0_5pn[6 + PN05DIM * PLUS_] = -3. * (2. - 3. * stp[2]);
    d0coeff0_5pn[6 + PN05DIM * MINUS] = -3. * (2. - 3. * stp[2]);
    d0coeff0_5pn[7 + PN05DIM * PLUS_] = +0.;
    d0coeff0_5pn[7 + PN05DIM * MINUS] = +0.;
    d0coeff0_5pn[8 + PN05DIM * PLUS_] = +0.;
    d0coeff0_5pn[8 + PN05DIM * MINUS] = +0.;
    d0coeff0_5pn[9 + PN05DIM * PLUS_] = +0.;
    d0coeff0_5pn[9 + PN05DIM * MINUS] = +0.;
    d0coeff0_5pn[10 + PN05DIM * PLUS_] = +3. * c2t;
    d0coeff0_5pn[10 + PN05DIM * MINUS] = +3. * c2t;
    for (i = 0; i <= 10; i++) {
        for (j = 0; j <= 1; j++) {
            params->ccoeff05pn[i + PN05DIM * j] = acoeff0_5pn[i + PN05DIM * j]
                    + skp[2] * (b0coeff0_5pn[i + PN05DIM * j] + ckp[1] * d0coeff0_5pn[i + PN05DIM * j]);
        }
    }
    /*
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%   1 PN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     */
    acoeff1pn[0 + PN10DIM * PLUS_] = +8. * k[PLUS_];
    acoeff1pn[0 + PN10DIM * MINUS] = -8. * k[MINUS];
    acoeff1pn[1 + PN10DIM * PLUS_] = +6. * k[PLUS_] * (stp[2] + 5.);
    acoeff1pn[1 + PN10DIM * MINUS] = +6. * k[MINUS] * (stp[2] + 5.);
    acoeff1pn[2 + PN10DIM * PLUS_] = +2. * k[PLUS_] * (4. - stp[2]);
    acoeff1pn[2 + PN10DIM * MINUS] = +2. * k[MINUS] * (4. - stp[2]);
    acoeff1pn[3 + PN10DIM * PLUS_] = +2. * k[PLUS_] * (2. * stp[4] + 11. * stp[2] - 38.);
    acoeff1pn[3 + PN10DIM * MINUS] = -2. * k[MINUS] * (2. * stp[4] + 11. * stp[2] - 38.);
    acoeff1pn[4 + PN10DIM * PLUS_] = +6. * k[PLUS_] * (3. * stp[2] + 5.);
    acoeff1pn[4 + PN10DIM * MINUS] = +6. * k[MINUS] * (3. * stp[2] + 5.);
    acoeff1pn[5 + PN10DIM * PLUS_] = +2. * k[PLUS_] * (4. * stp[2] + 19.);
    acoeff1pn[5 + PN10DIM * MINUS] = -2. * k[MINUS] * (4. * stp[2] + 19.);
    acoeff1pn[6 + PN10DIM * PLUS_] = -2. * k[PLUS_] * (3. * stp[2] - 4.);
    acoeff1pn[6 + PN10DIM * MINUS] = -2. * k[MINUS] * (3. * stp[2] - 4.);
    acoeff1pn[7 + PN10DIM * PLUS_] = +2. * k[PLUS_] * (4. - stp[2]);
    acoeff1pn[7 + PN10DIM * MINUS] = -2. * k[MINUS] * (4. - stp[2]);
    acoeff1pn[8 + PN10DIM * PLUS_] = +6. * k[PLUS_] * (5. + stp[2]);
    acoeff1pn[8 + PN10DIM * MINUS] = -6. * k[MINUS] * (5. + stp[2]);
    acoeff1pn[9 + PN10DIM * PLUS_] = -4. * k[PLUS_];
    acoeff1pn[9 + PN10DIM * MINUS] = +4. * k[MINUS];
    acoeff1pn[10 + PN10DIM * PLUS_] = k[PLUS_] * (22. + 29. * stp[2] - 16. * stp[4]);
    acoeff1pn[10 + PN10DIM * MINUS] = k[MINUS] * (22. + 29. * stp[2] - 16. * stp[4]);
    acoeff1pn[11 + PN10DIM * PLUS_] = +2. * k[PLUS_];
    acoeff1pn[11 + PN10DIM * MINUS] = +2. * k[MINUS];
    acoeff1pn[12 + PN10DIM * PLUS_] = +6. * k[PLUS_] * (3. * stp[2] + 5.);
    acoeff1pn[12 + PN10DIM * MINUS] = -6. * k[MINUS] * (3. * stp[2] + 5.);
    acoeff1pn[13 + PN10DIM * PLUS_] = -k[PLUS_] * (20. * stp[2] + 11.);
    acoeff1pn[13 + PN10DIM * MINUS] = -k[MINUS] * (20. * stp[2] + 11.);
    acoeff1pn[14 + PN10DIM * PLUS_] = -2. * k[PLUS_] * (3. * stp[2] - 4.);
    acoeff1pn[14 + PN10DIM * MINUS] = +2. * k[MINUS] * (3. * stp[2] - 4.);
    b0coeff1pn[0 + PN10DIM * PLUS_] = +8.;
    b0coeff1pn[0 + PN10DIM * MINUS] = +8.;
    b0coeff1pn[1 + PN10DIM * PLUS_] = -18. + 7. * stp[2];
    b0coeff1pn[1 + PN10DIM * MINUS] = +18. - 7. * stp[2];
    b0coeff1pn[2 + PN10DIM * PLUS_] = -3. * stp[2] + 6.;
    b0coeff1pn[2 + PN10DIM * MINUS] = +3. * stp[2] - 6.;
    b0coeff1pn[3 + PN10DIM * PLUS_] = (-22. - 29. * stp[2] + 16. * stp[4]);
    b0coeff1pn[3 + PN10DIM * MINUS] = (-22. - 29. * stp[2] + 16. * stp[4]);
    b0coeff1pn[4 + PN10DIM * PLUS_] = +26. * stp[2] - 18.;
    b0coeff1pn[4 + PN10DIM * MINUS] = -26. * stp[2] + 18.;
    b0coeff1pn[5 + PN10DIM * PLUS_] = +11. + 20. * stp[2];
    b0coeff1pn[5 + PN10DIM * MINUS] = +11. + 20. * stp[2];
    b0coeff1pn[6 + PN10DIM * PLUS_] = -6. * stp[2] + 6.;
    b0coeff1pn[6 + PN10DIM * MINUS] = +6. * stp[2] - 6.;
    b0coeff1pn[7 + PN10DIM * PLUS_] = +2. * (11. - 5. * stp[2]);
    b0coeff1pn[7 + PN10DIM * MINUS] = +2. * (11. - 5. * stp[2]);
    b0coeff1pn[8 + PN10DIM * PLUS_] = +6. * (7. + 9. * stp[2]);
    b0coeff1pn[8 + PN10DIM * MINUS] = +6. * (7. + 9. * stp[2]);
    b0coeff1pn[9 + PN10DIM * PLUS_] = -11.;
    b0coeff1pn[9 + PN10DIM * MINUS] = -11.;
    b0coeff1pn[10 + PN10DIM * PLUS_] = -3. * (8. - 20. * stp[2] + 7. * stp[4]);
    b0coeff1pn[10 + PN10DIM * MINUS] = +3. * (8. - 20. * stp[2] + 7. * stp[4]);
    b0coeff1pn[11 + PN10DIM * PLUS_] = +3.;
    b0coeff1pn[11 + PN10DIM * MINUS] = -3.;
    b0coeff1pn[12 + PN10DIM * PLUS_] = +3. * (19. * stp[2] + 14.);
    b0coeff1pn[12 + PN10DIM * MINUS] = +3. * (19. * stp[2] + 14.);
    b0coeff1pn[13 + PN10DIM * PLUS_] = +12. * c2t;
    b0coeff1pn[13 + PN10DIM * MINUS] = -12. * c2t;
    b0coeff1pn[14 + PN10DIM * PLUS_] = (22. - 21. * stp[2]);
    b0coeff1pn[14 + PN10DIM * MINUS] = (22. - 21. * stp[2]);
    d0coeff1pn[0 + PN10DIM * PLUS_] = -4.;
    d0coeff1pn[0 + PN10DIM * MINUS] = +4.;
    d0coeff1pn[1 + PN10DIM * PLUS_] = (6. - 14. * stp[2]);
    d0coeff1pn[1 + PN10DIM * MINUS] = (6. - 14. * stp[2]);
    d0coeff1pn[2 + PN10DIM * PLUS_] = +2. * (stp[2] - 1.);
    d0coeff1pn[2 + PN10DIM * MINUS] = +2. * (stp[2] - 1.);
    d0coeff1pn[3 + PN10DIM * PLUS_] = -2. * (8. - 20. * stp[2] + 7. * stp[4]);
    d0coeff1pn[3 + PN10DIM * MINUS] = +2. * (8. - 20. * stp[2] + 7. * stp[4]);
    d0coeff1pn[4 + PN10DIM * PLUS_] = (6. - 7. * stp[2]);
    d0coeff1pn[4 + PN10DIM * MINUS] = (6. - 7. * stp[2]);
    d0coeff1pn[5 + PN10DIM * PLUS_] = (-16. * stp[2] + 8.);
    d0coeff1pn[5 + PN10DIM * MINUS] = (+16. * stp[2] - 8.);
    d0coeff1pn[6 + PN10DIM * PLUS_] = (3. * stp[2] - 2.);
    d0coeff1pn[6 + PN10DIM * MINUS] = (3. * stp[2] - 2.);
    d0coeff1pn[7 + PN10DIM * PLUS_] = +9. * (stp[2] - 2.);
    d0coeff1pn[7 + PN10DIM * MINUS] = -9. * (stp[2] - 2.);
    d0coeff1pn[8 + PN10DIM * PLUS_] = +3. * (18. - 7. * stp[2]);
    d0coeff1pn[8 + PN10DIM * MINUS] = -3. * (18. - 7. * stp[2]);
    d0coeff1pn[9 + PN10DIM * PLUS_] = +9.;
    d0coeff1pn[9 + PN10DIM * MINUS] = -9.;
    d0coeff1pn[10 + PN10DIM * PLUS_] = +4. * (2. - 8. * stp[2] + 7. * stp[4]);
    d0coeff1pn[10 + PN10DIM * MINUS] = +4. * (2. - 8. * stp[2] + 7. * stp[4]);
    d0coeff1pn[11 + PN10DIM * PLUS_] = -2.;
    d0coeff1pn[11 + PN10DIM * MINUS] = -2.;
    d0coeff1pn[12 + PN10DIM * PLUS_] = +6. * (9. - 13. * stp[2]);
    d0coeff1pn[12 + PN10DIM * MINUS] = -6. * (9. - 13. * stp[2]);
    d0coeff1pn[13 + PN10DIM * PLUS_] = +2. * (7. * stp[2] - 2.);
    d0coeff1pn[13 + PN10DIM * MINUS] = +2. * (7. * stp[2] - 2.);
    d0coeff1pn[14 + PN10DIM * PLUS_] = -18. * ctp[2];
    d0coeff1pn[14 + PN10DIM * MINUS] = +18. * ctp[2];
    b1coeff1pn[0 + PN10DIM * PLUS_] = -1.;
    b1coeff1pn[0 + PN10DIM * MINUS] = -1.;
    b1coeff1pn[1 + PN10DIM * PLUS_] = +0.;
    b1coeff1pn[1 + PN10DIM * MINUS] = +0.;
    b1coeff1pn[2 + PN10DIM * PLUS_] = +0.;
    b1coeff1pn[2 + PN10DIM * MINUS] = +0.;
    b1coeff1pn[3 + PN10DIM * PLUS_] = -2. * (2. - 8. * stp[2] + 7. * stp[4]);
    b1coeff1pn[3 + PN10DIM * MINUS] = -2. * (2. - 8. * stp[2] + 7. * stp[4]);
    b1coeff1pn[4 + PN10DIM * PLUS_] = +0.;
    b1coeff1pn[4 + PN10DIM * MINUS] = +0.;
    b1coeff1pn[5 + PN10DIM * PLUS_] = (2. - 7. * stp[2]);
    b1coeff1pn[5 + PN10DIM * MINUS] = (2. - 7. * stp[2]);
    b1coeff1pn[6 + PN10DIM * PLUS_] = +0.;
    b1coeff1pn[6 + PN10DIM * MINUS] = +0.;
    b1coeff1pn[7 + PN10DIM * PLUS_] = -8. * ctp[2];
    b1coeff1pn[7 + PN10DIM * MINUS] = -8. * ctp[2];
    b1coeff1pn[8 + PN10DIM * PLUS_] = +8. * (3. - 7. * stp[2]);
    b1coeff1pn[8 + PN10DIM * MINUS] = +8. * (3. - 7. * stp[2]);
    b1coeff1pn[9 + PN10DIM * PLUS_] = +4.;
    b1coeff1pn[9 + PN10DIM * MINUS] = +4.;
    b1coeff1pn[10 + PN10DIM * PLUS_] = +0.;
    b1coeff1pn[10 + PN10DIM * MINUS] = +0.;
    b1coeff1pn[11 + PN10DIM * PLUS_] = +0.;
    b1coeff1pn[11 + PN10DIM * MINUS] = +0.;
    b1coeff1pn[12 + PN10DIM * PLUS_] = -4. * (7. * stp[2] - 6.);
    b1coeff1pn[12 + PN10DIM * MINUS] = -4. * (7. * stp[2] - 6.);
    b1coeff1pn[13 + PN10DIM * PLUS_] = +0.;
    b1coeff1pn[13 + PN10DIM * MINUS] = +0.;
    b1coeff1pn[14 + PN10DIM * PLUS_] = -4. * (2. - 3. * stp[2]);
    b1coeff1pn[14 + PN10DIM * MINUS] = -4. * (2. - 3. * stp[2]);
    for (i = 0; i < 15; i++) {
        d1coeff1pn[i + PN10DIM * PLUS_] = +0.;
        d1coeff1pn[i + PN10DIM * MINUS] = +0.;
    }
    for (i = 0; i <= 14; i++) {
        for (j = 0; j <= 1; j++) {
            params->ccoeff10pn[i + PN10DIM * j] = acoeff1pn[i + PN10DIM * j]
                    + skp[2] * (b0coeff1pn[i + PN10DIM * j] + ckp[1] * d0coeff1pn[i + PN10DIM * j])
                    + skp[4] * (b1coeff1pn[i + PN10DIM * j] + ckp[1] * d1coeff1pn[i + PN10DIM * j]);
        }
    }
    /*
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%   1.5 PN   %%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     */
    acoeff1_5pn[0 + PN15DIM * PLUS_] = +4. * k[PLUS_] * (stp[2] - 6.);
    acoeff1_5pn[0 + PN15DIM * MINUS] = -4. * k[MINUS] * (stp[2] - 6.);
    acoeff1_5pn[1 + PN15DIM * PLUS_] = +4. * k[PLUS_] * (stp[4] + 42. * stp[2] - 166.);
    acoeff1_5pn[1 + PN15DIM * MINUS] = -4. * k[MINUS] * (stp[4] + 42. * stp[2] - 166.);
    acoeff1_5pn[2 + PN15DIM * PLUS_] = +16. * k[PLUS_];
    acoeff1_5pn[2 + PN15DIM * MINUS] = +16. * k[MINUS];
    acoeff1_5pn[3 + PN15DIM * PLUS_] = +8. * k[PLUS_] * (stp[4] + 8. * stp[2] - 28.);
    acoeff1_5pn[3 + PN15DIM * MINUS] = +8. * k[MINUS] * (stp[4] + 8. * stp[2] - 28.);
    acoeff1_5pn[4 + PN15DIM * PLUS_] = +8. * k[PLUS_] * (-332. + 94. * stp[2] + stp[4]);
    acoeff1_5pn[4 + PN15DIM * MINUS] = +8. * k[MINUS] * (-332. + 94. * stp[2] + stp[4]);
    acoeff1_5pn[5 + PN15DIM * PLUS_] = +8. * k[PLUS_] * (38. - 42. * stp[2] - 9. * stp[4]);
    acoeff1_5pn[5 + PN15DIM * MINUS] = -8. * k[MINUS] * (38. - 42. * stp[2] - 9. * stp[4]);
    acoeff1_5pn[6 + PN15DIM * PLUS_] = -16. * k[PLUS_] * (152. - 46. * stp[2] - 9. * stp[4]);
    acoeff1_5pn[6 + PN15DIM * MINUS] = -16. * k[MINUS] * (152. - 46. * stp[2] - 9. * stp[4]);
    acoeff1_5pn[7 + PN15DIM * PLUS_] = +24. * k[PLUS_] * (3. * stp[2] - 10.);
    acoeff1_5pn[7 + PN15DIM * MINUS] = -24. * k[MINUS] * (3. * stp[2] - 10.);
    acoeff1_5pn[8 + PN15DIM * PLUS_] = -8. * k[PLUS_] * (160. - 204. * stp[2] - 63. * stp[4]);
    acoeff1_5pn[8 + PN15DIM * MINUS] = -8. * k[MINUS] * (160. - 204. * stp[2] - 63. * stp[4]);
    acoeff1_5pn[9 + PN15DIM * PLUS_] = +4. * k[PLUS_] * (3. - 2. * stp[2]);
    acoeff1_5pn[9 + PN15DIM * MINUS] = -4. * k[MINUS] * (3. - 2. * stp[2]);
    acoeff1_5pn[10 + PN15DIM * PLUS_] = -8. * k[PLUS_] * (14. + 3. * stp[2]);
    acoeff1_5pn[10 + PN15DIM * MINUS] = -8. * k[MINUS] * (14. + 3. * stp[2]);
    acoeff1_5pn[11 + PN15DIM * PLUS_] = -16. * k[PLUS_] * (15. * stp[2] + 76.);
    acoeff1_5pn[11 + PN15DIM * MINUS] = -16. * k[MINUS] * (15. * stp[2] + 76.);
    acoeff1_5pn[12 + PN15DIM * PLUS_] = -8. * k[PLUS_] * (5. * stp[2] + 166.);
    acoeff1_5pn[12 + PN15DIM * MINUS] = -8. * k[MINUS] * (5. * stp[2] + 166.);
    acoeff1_5pn[13 + PN15DIM * PLUS_] = -8. * k[PLUS_] * (80. + 63. * stp[2]);
    acoeff1_5pn[13 + PN15DIM * MINUS] = -8. * k[MINUS] * (80. + 63. * stp[2]);
    acoeff1_5pn[14 + PN15DIM * PLUS_] = +4. * k[PLUS_] * (166. - 125. * stp[2] - 8. * stp[4]);
    acoeff1_5pn[14 + PN15DIM * MINUS] = -4. * k[MINUS] * (166. - 125. * stp[2] - 8. * stp[4]);
    acoeff1_5pn[15 + PN15DIM * PLUS_] = -8. * k[PLUS_] * (38. - 61. * stp[2] - 24. * stp[4]);
    acoeff1_5pn[15 + PN15DIM * MINUS] = +8. * k[MINUS] * (38. - 61. * stp[2] - 24. * stp[4]);
    acoeff1_5pn[16 + PN15DIM * PLUS_] = +8. * k[PLUS_] * (5. - 4. * stp[2]);
    acoeff1_5pn[16 + PN15DIM * MINUS] = -8. * k[MINUS] * (5. - 4. * stp[2]);
    b0coeff1_5pn[0 + PN15DIM * PLUS_] = (5. * stp[2] - 6.);
    b0coeff1_5pn[0 + PN15DIM * MINUS] = (5. * stp[2] - 6.);
    b0coeff1_5pn[1 + PN15DIM * PLUS_] = (18. * stp[4] + 252. * stp[2] - 188.);
    b0coeff1_5pn[1 + PN15DIM * MINUS] = (18. * stp[4] + 252. * stp[2] - 188.);
    b0coeff1_5pn[2 + PN15DIM * PLUS_] = +20.;
    b0coeff1_5pn[2 + PN15DIM * MINUS] = -20.;
    b0coeff1_5pn[3 + PN15DIM * PLUS_] = (+9. * stp[4] - 90. * stp[2] + 56.);
    b0coeff1_5pn[3 + PN15DIM * MINUS] = (-9. * stp[4] + 90. * stp[2] - 56.);
    b0coeff1_5pn[4 + PN15DIM * PLUS_] = -4. * (1184. - 172. * stp[2] - 7. * stp[4]);
    b0coeff1_5pn[4 + PN15DIM * MINUS] = +4. * (1184. - 172. * stp[2] - 7. * stp[4]);
    b0coeff1_5pn[5 + PN15DIM * PLUS_] = +2. * (46. + 48. * stp[2] - 99. * stp[4]);
    b0coeff1_5pn[5 + PN15DIM * MINUS] = +2. * (46. + 48. * stp[2] - 99. * stp[4]);
    b0coeff1_5pn[6 + PN15DIM * PLUS_] = -12. * (72. + 110. * stp[2] - 63. * stp[4]);
    b0coeff1_5pn[6 + PN15DIM * MINUS] = +12. * (72. + 110. * stp[2] - 63. * stp[4]);
    b0coeff1_5pn[7 + PN15DIM * PLUS_] = +144. * (stp[2] - 2.);
    b0coeff1_5pn[7 + PN15DIM * MINUS] = +144. * (stp[2] - 2.);
    b0coeff1_5pn[8 + PN15DIM * PLUS_] = -3. * (-204. + 406. * stp[2] - 189. * stp[4]);
    b0coeff1_5pn[8 + PN15DIM * MINUS] = +3. * (-204. + 406. * stp[2] - 189. * stp[4]);
    b0coeff1_5pn[9 + PN15DIM * PLUS_] = +3. - 4. * stp[2];
    b0coeff1_5pn[9 + PN15DIM * MINUS] = +3. - 4. * stp[2];
    b0coeff1_5pn[10 + PN15DIM * PLUS_] = +28. - 31. * stp[2];
    b0coeff1_5pn[10 + PN15DIM * MINUS] = -28. + 31. * stp[2];
    b0coeff1_5pn[11 + PN15DIM * PLUS_] = -432. - 876. * stp[2];
    b0coeff1_5pn[11 + PN15DIM * MINUS] = +432. + 876. * stp[2];
    b0coeff1_5pn[12 + PN15DIM * PLUS_] = -4. * (71. * stp[2] + 592.);
    b0coeff1_5pn[12 + PN15DIM * MINUS] = +4. * (71. * stp[2] + 592.);
    b0coeff1_5pn[13 + PN15DIM * PLUS_] = +306. - 651. * stp[2];
    b0coeff1_5pn[13 + PN15DIM * MINUS] = -306. + 651. * stp[2];
    b0coeff1_5pn[14 + PN15DIM * PLUS_] = +2. * (94. - 173. * stp[2] - 24. * stp[4]);
    b0coeff1_5pn[14 + PN15DIM * MINUS] = +2. * (94. - 173. * stp[2] - 24. * stp[4]);
    b0coeff1_5pn[15 + PN15DIM * PLUS_] = -2. * (46. + 25. * stp[2] - 180. * stp[4]);
    b0coeff1_5pn[15 + PN15DIM * MINUS] = -2. * (46. + 25. * stp[2] - 180. * stp[4]);
    b0coeff1_5pn[16 + PN15DIM * PLUS_] = +48. * ctp[2];
    b0coeff1_5pn[16 + PN15DIM * MINUS] = +48. * ctp[2];
    d0coeff1_5pn[0 + PN15DIM * PLUS_] = +0.;
    d0coeff1_5pn[0 + PN15DIM * MINUS] = +0.;
    d0coeff1_5pn[1 + PN15DIM * PLUS_] = (-6. * stp[4] + 72. * stp[2] - 20.);
    d0coeff1_5pn[1 + PN15DIM * MINUS] = (+6. * stp[4] - 72. * stp[2] + 20.);
    d0coeff1_5pn[2 + PN15DIM * PLUS_] = -12.;
    d0coeff1_5pn[2 + PN15DIM * MINUS] = -12.;
    d0coeff1_5pn[3 + PN15DIM * PLUS_] = (-15. * stp[4] + 22. * stp[2] - 8.);
    d0coeff1_5pn[3 + PN15DIM * MINUS] = (-15. * stp[4] + 22. * stp[2] - 8.);
    d0coeff1_5pn[4 + PN15DIM * PLUS_] = (1920. - 2832. * stp[2] - 84. * stp[4]);
    d0coeff1_5pn[4 + PN15DIM * MINUS] = (1920. - 2832. * stp[2] - 84. * stp[4]);
    d0coeff1_5pn[5 + PN15DIM * PLUS_] = +6. * (10. - 44. * stp[2] + 27. * stp[4]);
    d0coeff1_5pn[5 + PN15DIM * MINUS] = -6. * (10. - 44. * stp[2] + 27. * stp[4]);
    d0coeff1_5pn[6 + PN15DIM * PLUS_] = -4. * (88. - 422. * stp[2] + 171. * stp[4]);
    d0coeff1_5pn[6 + PN15DIM * MINUS] = -4. * (88. - 422. * stp[2] + 171. * stp[4]);
    d0coeff1_5pn[7 + PN15DIM * PLUS_] = +12. * (14. - 9. * stp[2]);
    d0coeff1_5pn[7 + PN15DIM * MINUS] = -12. * (14. - 9. * stp[2]);
    d0coeff1_5pn[8 + PN15DIM * PLUS_] = -9. * (28. - 126. * stp[2] + 105. * stp[4]);
    d0coeff1_5pn[8 + PN15DIM * MINUS] = -9. * (28. - 126. * stp[2] + 105. * stp[4]);
    d0coeff1_5pn[9 + PN15DIM * PLUS_] = +0.;
    d0coeff1_5pn[9 + PN15DIM * MINUS] = +0.;
    d0coeff1_5pn[10 + PN15DIM * PLUS_] = (9. * stp[2] - 4.);
    d0coeff1_5pn[10 + PN15DIM * MINUS] = (9. * stp[2] - 4.);
    d0coeff1_5pn[11 + PN15DIM * PLUS_] = (-176. + 756. * stp[2]);
    d0coeff1_5pn[11 + PN15DIM * MINUS] = (-176. + 756. * stp[2]);
    d0coeff1_5pn[12 + PN15DIM * PLUS_] = +12. * (7. * stp[2] + 80.);
    d0coeff1_5pn[12 + PN15DIM * MINUS] = +12. * (7. * stp[2] + 80.);
    d0coeff1_5pn[13 + PN15DIM * PLUS_] = (-126. + 189. * stp[2]);
    d0coeff1_5pn[13 + PN15DIM * MINUS] = (-126. + 189. * stp[2]);
    d0coeff1_5pn[14 + PN15DIM * PLUS_] = +2. * (10. - 41. * stp[2] + 36. * stp[4]);
    d0coeff1_5pn[14 + PN15DIM * MINUS] = -2. * (10. - 41. * stp[2] + 36. * stp[4]);
    d0coeff1_5pn[15 + PN15DIM * PLUS_] = -6. * (10. - 49. * stp[2] + 44. * stp[4]);
    d0coeff1_5pn[15 + PN15DIM * MINUS] = +6. * (10. - 49. * stp[2] + 44. * stp[4]);
    d0coeff1_5pn[16 + PN15DIM * PLUS_] = -4. * (7. - 8. * stp[2]);
    d0coeff1_5pn[16 + PN15DIM * MINUS] = +4. * (7. - 8. * stp[2]);
    b1coeff1_5pn[0 + PN15DIM * PLUS_] = +0.;
    b1coeff1_5pn[0 + PN15DIM * MINUS] = +0.;
    b1coeff1_5pn[1 + PN15DIM * PLUS_] = (-15. * stp[4] + 12. * stp[2] - 2.);
    b1coeff1_5pn[1 + PN15DIM * MINUS] = (-15. * stp[4] + 12. * stp[2] - 2.);
    b1coeff1_5pn[2 + PN15DIM * PLUS_] = -5.;
    b1coeff1_5pn[2 + PN15DIM * MINUS] = +5.;
    b1coeff1_5pn[3 + PN15DIM * PLUS_] = +0.;
    b1coeff1_5pn[3 + PN15DIM * MINUS] = +0.;
    b1coeff1_5pn[4 + PN15DIM * PLUS_] = -(236. - 294. * stp[2] + 21. * stp[4]);
    b1coeff1_5pn[4 + PN15DIM * MINUS] = +(236. - 294. * stp[2] + 21. * stp[4]);
    b1coeff1_5pn[5 + PN15DIM * PLUS_] = +3. * (6. - 36. * stp[2] + 45. * stp[4]);
    b1coeff1_5pn[5 + PN15DIM * MINUS] = +3. * (6. - 36. * stp[2] + 45. * stp[4]);
    b1coeff1_5pn[6 + PN15DIM * PLUS_] = -3. * (232. - 510. * stp[2] + 243. * stp[4]);
    b1coeff1_5pn[6 + PN15DIM * MINUS] = +3. * (232. - 510. * stp[2] + 243. * stp[4]);
    b1coeff1_5pn[7 + PN15DIM * PLUS_] = +9. * (6. - 5. * stp[2]);
    b1coeff1_5pn[7 + PN15DIM * MINUS] = +9. * (6. - 5. * stp[2]);
    b1coeff1_5pn[8 + PN15DIM * PLUS_] = +0.;
    b1coeff1_5pn[8 + PN15DIM * MINUS] = +0.;
    b1coeff1_5pn[9 + PN15DIM * PLUS_] = +0.;
    b1coeff1_5pn[9 + PN15DIM * MINUS] = +0.;
    b1coeff1_5pn[10 + PN15DIM * PLUS_] = +0.;
    b1coeff1_5pn[10 + PN15DIM * MINUS] = +0.;
    b1coeff1_5pn[11 + PN15DIM * PLUS_] = (-348. + 591. * stp[2]);
    b1coeff1_5pn[11 + PN15DIM * MINUS] = (+348. - 591. * stp[2]);
    b1coeff1_5pn[12 + PN15DIM * PLUS_] = (+273. * stp[2] - 118.);
    b1coeff1_5pn[12 + PN15DIM * MINUS] = (-273. * stp[2] + 118.);
    b1coeff1_5pn[13 + PN15DIM * PLUS_] = +0.;
    b1coeff1_5pn[13 + PN15DIM * MINUS] = +0.;
    b1coeff1_5pn[14 + PN15DIM * PLUS_] = (2. - 13. * stp[2] + 12. * stp[4]);
    b1coeff1_5pn[14 + PN15DIM * MINUS] = (2. - 13. * stp[2] + 12. * stp[4]);
    b1coeff1_5pn[15 + PN15DIM * PLUS_] = -9. * (2. - 13. * stp[2] + 12. * stp[4]);
    b1coeff1_5pn[15 + PN15DIM * MINUS] = -9. * (2. - 13. * stp[2] + 12. * stp[4]);
    b1coeff1_5pn[16 + PN15DIM * PLUS_] = -3. * (3. - 4. * stp[2]);
    b1coeff1_5pn[16 + PN15DIM * MINUS] = -3. * (3. - 4. * stp[2]);
    d1coeff1_5pn[0 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[0 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[1 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[1 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[2 + PN15DIM * PLUS_] = +1.;
    d1coeff1_5pn[2 + PN15DIM * MINUS] = +1.;
    d1coeff1_5pn[3 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[3 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[4 + PN15DIM * PLUS_] = (28. - 126. * stp[2] + 105. * stp[4]);
    d1coeff1_5pn[4 + PN15DIM * MINUS] = (28. - 126. * stp[2] + 105. * stp[4]);
    d1coeff1_5pn[5 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[5 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[6 + PN15DIM * PLUS_] = +27. * (8. - 22. * stp[2] + 15. * stp[4]);
    d1coeff1_5pn[6 + PN15DIM * MINUS] = +27. * (8. - 22. * stp[2] + 15. * stp[4]);
    d1coeff1_5pn[7 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[7 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[8 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[8 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[9 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[9 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[10 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[10 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[11 + PN15DIM * PLUS_] = (-243. * stp[2] + 108.);
    d1coeff1_5pn[11 + PN15DIM * MINUS] = (-243. * stp[2] + 108.);
    d1coeff1_5pn[12 + PN15DIM * PLUS_] = +7. * (2. - 3. * stp[2]);
    d1coeff1_5pn[12 + PN15DIM * MINUS] = +7. * (2. - 3. * stp[2]);
    d1coeff1_5pn[13 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[13 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[14 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[14 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[15 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[15 + PN15DIM * MINUS] = +0.;
    d1coeff1_5pn[16 + PN15DIM * PLUS_] = +0.;
    d1coeff1_5pn[16 + PN15DIM * MINUS] = +0.;
    for (i = 0; i <= 16; i++) {
        for (j = 0; j <= 1; j++) {
            params->ccoeff15pn[i + PN15DIM * j] = acoeff1_5pn[i + PN15DIM * j]
                    + skp[2] * (b0coeff1_5pn[i + PN15DIM * j] + ckp[1] * d0coeff1_5pn[i + PN15DIM * j])
                    + skp[4] * (b1coeff1_5pn[i + PN15DIM * j] + ckp[1] * d1coeff1_5pn[i + PN15DIM * j]);
        }
    }
    XLALFreeDmatrix(acoeff00pn);
    XLALFreeDmatrix(b0coeff0pn);
    XLALFreeDmatrix(d0coeff0pn);
    XLALFreeDmatrix(acoeff0_5pn);
    XLALFreeDmatrix(b0coeff0_5pn);
    XLALFreeDmatrix(d0coeff0_5pn);
    XLALFreeDmatrix(acoeff1pn);
    XLALFreeDmatrix(b0coeff1pn);
    XLALFreeDmatrix(d0coeff1pn);
    XLALFreeDmatrix(b1coeff1pn);
    XLALFreeDmatrix(d1coeff1pn);
    XLALFreeDmatrix(acoeff1_5pn);
    XLALFreeDmatrix(b0coeff1_5pn);
    XLALFreeDmatrix(d0coeff1_5pn);
    XLALFreeDmatrix(b1coeff1_5pn);
    XLALFreeDmatrix(d1coeff1_5pn);
    return XLAL_SUCCESS;
}

/**
 * Function building the wavefrom from the calculated parameters at a given time
 * For the formulae see the appendix of Arxiv:1209.1722
 */
int XLALSpinDominatedWaveformBuild(LALSDWaveformParams *params, /**< The SDW parameters */
REAL8 expr[], /**< The 3 time dependent variables of the waveform at the time indexed by idx */
REAL8TimeSeries **hplus, /**< +-polarization waveform */
REAL8TimeSeries **hcross, /**< x-polarization waveform */
int idx /**< index of the time dependent variables */) {
    enum {
        PN00_A, PN00_B, PN05_A, PN05_B, PN10_A, PN10_B, PN10_C, PN10_D, PN15_A, PN15_B, PN15_C,
    };
    REAL8 omega0 = 1;
    REAL8 chi1_1 = params->chi1;
    REAL8 skp[TRIGONOMETRIC_POWER], ckp[TRIGONOMETRIC_POWER], stp[TRIGONOMETRIC_POWER], ctp[TRIGONOMETRIC_POWER];
    skp[0] = ckp[0] = stp[0] = ctp[0] = 1.;
    skp[1] = sin(params->kappa1);
    ckp[1] = cos(params->kappa1);
    stp[1] = sin(params->theta);
    ctp[1] = cos(params->theta);
    for (int i = 2; i < TRIGONOMETRIC_POWER; ++i) {
        skp[i] = skp[1] * skp[i - 1];
        ckp[i] = ckp[1] * ckp[i - 1];
        stp[i] = stp[1] * stp[i - 1];
        ctp[i] = ctp[1] * ctp[i - 1];
    }
    REAL8 s2t = sin(2. * params->theta);
    REAL8 c2t = cos(2. * params->theta);
    REAL8 k[PLUS_MINUS_DIM] = { ckp[1] - 1., ckp[1] + 1. };
    REAL8 s2k1 = sin(2. * params->kappa1);
    REAL8 c2k1 = cos(2. * params->kappa1);

    REAL8 vP[OMEGA_POWER_DIM];
    vP[0] = 1.;
    vP[1] = cbrt(params->totalmass * expr[OMEGA] * G_CP2 / LAL_C_SI);
    for (int i = 2; i < OMEGA_POWER_DIM; ++i) {
        vP[i] = vP[1] * vP[i - 1];
    }
    params->eps = vP[2];
    REAL8 eta = params->nu / (1. + params->nu) / (1. + params->nu);
    REAL8 eps_corr[PN_ORDER][PN_ORDER];
    eps_corr[PN00][PN10] = params->eps
            * (1. + vP[2] * (1. - eta / 3.) );
    eps_corr[PN05][PN10] = eps_corr[PN00][PN10];
    eps_corr[PN00][PN15] = params->eps
            * (1. + vP[2] * (1. - eta / 3.) + vP[3] * (eta + (2. / 3.) / (eta / params->nu)));
    REAL8 epssqrt = sqrt(params->eps);
    params->xi = params->nu / epssqrt; // second small parameter

    REAL8 phin[PHI_PSI_DIM], psin[PHI_PSI_DIM];
    phin[0] = psin[0] = 1.;
    for (int i = 1; i < PHI_PSI_DIM; ++i) {
        phin[i] = i * expr[PHI];
        psin[i] = i * expr[PSI];

    }
    REAL8 *waveampcoeffs = XLALDmatrix(AMPCOEFF_DIM, PLUS_MINUS_DIM);
    switch (params->pnamp) {
    case (PNDEF): //default value -1 contains all corrections
    case (PN15): {
        waveampcoeffs[PN15_C + AMPCOEFF_DIM * PLUS_] = LAL_PI_2 * (    //
                6. * skp[2] * stp[2] * cos(psin[2]) + (    //
                        (stp[2] - 2.) * (    //
                                /*  */cos(phin[2] + psin[2]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                                /**/+ cos(phin[2] - psin[2]) * params->ccoeff00pn[0 + PN00DIM * MINUS]    //
                                ) - 2. * skp[1] * s2t * (    //
                                /*  */sin(phin[1] + psin[2]) * k[PLUS_]    //
                                /**/+ sin(phin[1] - psin[2]) * k[MINUS]    //
                                )//
                        )//
                ) + 3. * log(eps_corr[PN00][PN15] / omega0) * (    //
                4. * ctp[1] * skp[1] * stp[1] * (    //
                        -k[MINUS] * cos(phin[1] - psin[2])
                        + k[PLUS_] * cos(phin[1] + psin[2])    //
                        )//
                + sin(phin[2] - psin[2]) * (2. * ckp[1] - skp[2] + 2.) * (-stp[2] + 2.)  //
                + sin(phin[2] + psin[2]) * (2. * ckp[1] + skp[2] - 2.) //
                + 6. * skp[2] * sin(psin[2]) * stp[2] //
                );
        waveampcoeffs[PN15_C + AMPCOEFF_DIM * MINUS] = LAL_PI * ( //
                ctp[1] * ( //
                        /*  */sin(phin[2] + psin[2]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                        /**/+ sin(phin[2] - psin[2]) * params->ccoeff00pn[0 + PN00DIM * MINUS] //
                        ) - 2. * skp[1] * stp[1] * ( //
                        /*  */cos(phin[1] + psin[2]) * k[PLUS_]
                        /**/+ cos(phin[1] - psin[2]) * k[MINUS] //
                        )//
                ) + 6. * log(eps_corr[PN00][PN15] / omega0) * ( //
                (2. * ckp[1] - skp[2] + 2.) * ctp[1] * cos(phin[2] - psin[2])
                + (2. * ckp[1] + skp[2] - 2.) * ctp[1] * cos(phin[2] + psin[2])
                + 2. * k[MINUS] * skp[1] * stp[1] * sin(phin[1] - psin[2])
                - 2. * k[PLUS_] * skp[1] * stp[1] * sin(phin[1] + psin[2]) //
                );

        waveampcoeffs[PN15_B + AMPCOEFF_DIM * PLUS_] = chi1_1 / 2. * (4. * skp[1] * ( //
                /**/ckp[1] * skp[1] * cos(phin[2]) - c2k1 * ctp[1] * stp[1] * sin(phin[1])
                /**/+ ckp[1] * skp[1] * stp[2] * ( //
                        /*  */6. * sin(psin[1]) * sin(psin[1]) - 2. + sin(phin[1]) * sin(phin[1]) //
                        )//
                ) + ( //
                2. * ctp[1] * skp[1] * stp[1] * ( //
                        /*  */(-3. * k[PLUS_] - 4. * skp[2]) * sin(phin[1] + psin[2])
                        /**/+ (+3. * k[MINUS] - 4. * skp[2]) * sin(phin[1] - psin[2]) //
                        ) + (stp[2] - 2.) * ( //
                        /*  */(-2. * k[PLUS_] + (2. * ckp[1] - 3.) * skp[2]) * cos(phin[2] + psin[2])
                        /**/+ (-2. * k[MINUS] + (2. * ckp[1] + 3.) * skp[2]) * cos(phin[2] - psin[2]) //
                        )//
                ));
        waveampcoeffs[PN15_B + AMPCOEFF_DIM * MINUS] = chi1_1 * ( //
                -2. * cos(phin[1]) * skp[1] * (stp[1] * c2k1 + ctp[1] * s2k1 * sin(phin[1])) + ( //
                        ctp[1] * ( //
                                /*  */(-2. * k[PLUS_] + (2. * ckp[1] - 3.) * skp[2]) * sin(phin[2] + psin[2])
                                /**/+ (-2. * k[MINUS] + (2. * ckp[1] + 3.) * skp[2]) * sin(phin[2] - psin[2]) //
                                ) + skp[1] * stp[1] * ( //
                                /*  */(-3. * k[PLUS_] - 4. * skp[2]) * cos(phin[1] + psin[2])
                                /**/+ (+3. * k[MINUS] - 4. * skp[2]) * cos(phin[1] - psin[2]) //
                                )//
                        )//
                );
        waveampcoeffs[PN15_A + AMPCOEFF_DIM * PLUS_] = 1. / 12288. * (12. * ctp[1] * skp[1] * stp[2] * ( //
                cos(psin[3]) * ( //
                        /**/1701. * (2. - 3. * stp[2]) * skp[4] + 72. * skp[2] * (63. * stp[2] + 178.) //
                        ) + cos(psin[1]) * ( //
                        /**/-14. * (2. - 3. * stp[2]) * skp[4] //
                        /**/- 8. * skp[2] * (7. * stp[2] + 162.) + 16. * (stp[2] + 66.) //
                        ) - 4375. * (2. - 3. * stp[2]) * skp[4] * cos(psin[5]) //
                ) + ( //
                2. * (stp[2] - 2.) * skp[4] * stp[3] * ( //
                        /*  */sin(phin[5] + psin[1]) * k[PLUS_]
                        /**/+ sin(phin[5] - psin[1]) * k[MINUS] //
                        ) + 4. * ctp[1] * skp[3] * stp[2] * ( //
                        /*  */cos(phin[4] + psin[1]) * params->ccoeff15pn[0 + PN15DIM * PLUS_]
                        /**/+ cos(phin[4] - psin[1]) * params->ccoeff15pn[0 + PN15DIM * MINUS] //
                        ) + 16. * ctp[1] * skp[1] * ( //
                        /*  */cos(phin[2] + psin[1]) * params->ccoeff15pn[1 + PN15DIM * PLUS_]
                        /**/+ cos(phin[2] - psin[1]) * params->ccoeff15pn[1 + PN15DIM * MINUS] //
                        ) + 1250. * skp[4] * stp[1] * (105. * stp[4] - 126. * stp[2] + 28.) * ( //
                        /*  */sin(phin[1] + psin[5]) * k[PLUS_]
                        /**/+ sin(phin[1] - psin[5]) * k[MINUS] //
                        ) + 625. * (stp[2] - 2.) * stp[3] * ( //
                        /*  */sin(phin[5] + psin[5]) * params->ccoeff15pn[2 + PN15DIM * PLUS_]
                        /**/+ sin(phin[5] - psin[5]) * params->ccoeff15pn[2 + PN15DIM * MINUS] //
                        ) + 6. * skp[2] * stp[1] * ( //
                        /*  */sin(phin[3] + psin[1]) * params->ccoeff15pn[3 + PN15DIM * PLUS_]
                        /**/+ sin(phin[3] - psin[1]) * params->ccoeff15pn[3 + PN15DIM * MINUS] //
                        ) + 243. * (stp[2] - 2.) * skp[2] * stp[3] * ( //
                        /*  */sin(phin[5] + psin[3]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                        /**/+ sin(phin[5] - psin[3]) * params->ccoeff05pn[1 + PN05DIM * MINUS] //
                        ) + 4. * stp[1] * ( //
                        /*  */sin(phin[1] + psin[1]) * params->ccoeff15pn[4 + PN15DIM * PLUS_]
                        /**/+ sin(phin[1] - psin[1]) * params->ccoeff15pn[4 + PN15DIM * MINUS] //
                        ) + 5000. * ctp[1] * skp[3] * (15. * stp[4] - 12. * stp[2] + 2.) * ( //
                        /*  */cos(phin[2] + psin[5]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                        /**/+ cos(phin[2] - psin[5]) * params->ccoeff00pn[0 + PN00DIM * MINUS] //
                        ) - 1250. * ctp[1] * skp[1] * stp[2] * (5. * stp[2] - 6.) * ( //
                        /*  */cos(phin[4] + psin[5]) * params->ccoeff10pn[0 + PN10DIM * PLUS_]
                        /**/+ cos(phin[4] - psin[5]) * params->ccoeff10pn[0 + PN10DIM * MINUS] //
                        ) + 1875. * skp[2] * stp[1] * (8. - 22. * stp[2] + 15. * stp[4]) * ( //
                        /*  */sin(phin[3] + psin[5]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                        /**/+ sin(phin[3] - psin[5]) * params->ccoeff05pn[1 + PN05DIM * MINUS] //
                        ) + 216. * ctp[1] * skp[1] * ( //
                        /*  */cos(phin[2] + psin[3]) * params->ccoeff15pn[5 + PN15DIM * PLUS_]
                        /**/+ cos(phin[2] - psin[3]) * params->ccoeff15pn[5 + PN15DIM * MINUS] //
                        ) + 27. * stp[1] * ( //
                        /*  */sin(phin[3] + psin[3]) * params->ccoeff15pn[6 + PN15DIM * PLUS_]
                        /**/+ sin(phin[3] - psin[3]) * params->ccoeff15pn[6 + PN15DIM * MINUS] //
                        ) + 54. * ctp[1] * skp[1] * stp[2] * ( //
                        /*  */cos(phin[4] + psin[3]) * params->ccoeff15pn[7 + PN15DIM * PLUS_]
                        /**/+ cos(phin[4] - psin[3]) * params->ccoeff15pn[7 + PN15DIM * MINUS] //
                        ) + 54. * skp[2] * stp[1] * ( //
                        /*  */sin(phin[1] + psin[3]) * params->ccoeff15pn[8 + PN15DIM * PLUS_]
                        /**/+ sin(phin[1] - psin[3]) * params->ccoeff15pn[8 + PN15DIM * MINUS] //
                        )//
                ));
        waveampcoeffs[PN15_A + AMPCOEFF_DIM * MINUS] = 1. / 6144. * (192. * ckp[1] * skp[1] * stp[2] * ( //
                /*        */sin(psin[1]) * (64. - skp[2] * (7. * stp[2] - 6.) + 4. * stp[2])
                /**/+ 27. * sin(psin[3]) * skp[2] * (7. * stp[2] - 6.) //
                ) + ( //
                4. * skp[3] * stp[2] * ( //
                        /*  */sin(phin[4] + psin[1]) * params->ccoeff15pn[9 + PN15DIM * PLUS_]
                        /**/+ sin(phin[4] - psin[1]) * params->ccoeff15pn[9 + PN15DIM * MINUS] //
                        ) - 2. * ctp[1] * skp[4] * stp[3] * ( //
                        /*  */cos(phin[5] + psin[1]) * k[PLUS_]
                        /**/+ cos(phin[5] - psin[1]) * k[MINUS] //
                        ) - 243. * ctp[1] * skp[2] * stp[3] * ( //
                        /*  */cos(phin[5] + psin[3]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                        /**/+ cos(phin[5] - psin[3]) * params->ccoeff05pn[1 + PN05DIM * MINUS] //
                        ) - 625. * ctp[1] * stp[3] * ( //
                        /*  */cos(phin[5] + psin[5]) * params->ccoeff15pn[2 + PN15DIM * PLUS_]
                        /**/+ cos(phin[5] - psin[5]) * params->ccoeff15pn[2 + PN15DIM * MINUS] //
                        ) + 3. * skp[2] * s2t * ( //
                        /*  */cos(phin[3] + psin[1]) * params->ccoeff15pn[10 + PN15DIM * PLUS_]
                        /**/+ cos(phin[3] - psin[1]) * params->ccoeff15pn[10 + PN15DIM * MINUS] //
                        ) + 27. * ctp[1] * stp[1] * ( //
                        /*  */cos(phin[3] + psin[3]) * params->ccoeff15pn[11 + PN15DIM * PLUS_]
                        /**/+ cos(phin[3] - psin[3]) * params->ccoeff15pn[11 + PN15DIM * MINUS] //
                        ) + 1875. * ctp[1] * skp[2] * stp[1] * (4. - 9. * stp[2]) * ( //
                        /*  */cos(phin[3] + psin[5]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                        /**/+ cos(phin[3] - psin[5]) * params->ccoeff05pn[1 + PN05DIM * MINUS] //
                        ) + 2. * s2t * ( //
                        /*  */cos(phin[1] + psin[1]) * params->ccoeff15pn[12 + PN15DIM * PLUS_]
                        /**/+ cos(phin[1] - psin[1]) * params->ccoeff15pn[12 + PN15DIM * MINUS] //
                        ) + 27. * skp[2] * s2t * ( //
                        /*  */cos(phin[1] + psin[3]) * params->ccoeff15pn[13 + PN15DIM * PLUS_]
                        /**/+ cos(phin[1] - psin[3]) * params->ccoeff15pn[13 + PN15DIM * MINUS] //
                        ) + 8. * skp[1] * ( //
                        /*  */sin(phin[2] + psin[1]) * params->ccoeff15pn[14 + PN15DIM * PLUS_]
                        /**/+ sin(phin[2] - psin[1]) * params->ccoeff15pn[14 + PN15DIM * MINUS] //
                        ) + 4375. * (2. - 3. * stp[2]) * skp[4] * s2t * ( //
                        /*  */cos(phin[1] + psin[5]) * k[PLUS_]
                        /**/+ cos(phin[1] - psin[5]) * k[MINUS] //
                        ) + 108. * skp[1] * ( //
                        /*  */sin(phin[2] + psin[3]) * params->ccoeff15pn[15 + PN15DIM * PLUS_]
                        /**/+ sin(phin[2] - psin[3]) * params->ccoeff15pn[15 + PN15DIM * MINUS] //
                        ) + 162. * skp[1] * stp[2] * ( //
                        /*  */sin(phin[4] + psin[3]) * params->ccoeff15pn[16 + PN15DIM * PLUS_]
                        /**/+ sin(phin[4] - psin[3]) * params->ccoeff15pn[16 + PN15DIM * MINUS] //
                        ) - 2500. * skp[3] * (2. - 13. * stp[2] + 12. * stp[4]) * ( //
                        /*  */sin(phin[2] + psin[5]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                        /**/+ sin(phin[2] - psin[5]) * params->ccoeff00pn[0 + PN00DIM * MINUS] //
                        ) - 1250. * skp[1] * stp[2] * (3. - 4. * stp[2]) * ( //
                        /*  */sin(phin[4] + psin[5]) * params->ccoeff10pn[0 + PN10DIM * PLUS_]
                        /**/+ sin(phin[4] - psin[5]) * params->ccoeff10pn[0 + PN10DIM * MINUS] //
                        )//
                ));
    }
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN10): {
        REAL8 phi1_1 = LAL_PI_2;
        waveampcoeffs[PN10_D + AMPCOEFF_DIM * PLUS_] = chi1_1 / 2. * (    //
                (    //
                /* */-k[PLUS_] * ctp[1] * sin(phin[2] + psin[1] - phi1_1)
                /**/+ k[MINUS] * ctp[1] * sin(phin[2] - psin[1] - phi1_1)    //
                ) + skp[1] * stp[1] * (    //
                        /*  */sin(phin[1] + psin[1]) - sin(phin[1] - psin[1])
                        /**/+ cos(phin[1] + psin[1] - phi1_1) - cos(phin[1] - psin[1] - phi1_1)    //
                        )//
                );
        waveampcoeffs[PN10_D + AMPCOEFF_DIM * MINUS] = chi1_1 / 4. * (-stp[2] * (    //
                2. * (sin(phi1_1) + 2.) * ckp[1] * sin(psin[1]) + 2. * cos(phi1_1) * cos(psin[1])    //
                ) + (    //
                2. * ctp[1] * skp[1] * stp[1] * (    //
                        /**/cos(phin[1] + psin[1]) - cos(phin[1] - psin[1])
                        /**/- sin(phin[1] + psin[1] - phi1_1) + sin(phin[1] - psin[1] - phi1_1)    //
                        ) + (stp[2] - 2.) * (    //
                        /*  */k[PLUS_] * cos(phin[2] + psin[1] - phi1_1)
                        /**/- k[MINUS] * cos(phin[2] - psin[1] - phi1_1) //
                        )//
                ));
        waveampcoeffs[PN10_C + AMPCOEFF_DIM * PLUS_] = chi1_1 / 2. * stp[1]
                * (k[PLUS_] * sin(phin[1] + psin[1]) - k[MINUS] * sin(phin[1] - psin[1]));
        waveampcoeffs[PN10_C + AMPCOEFF_DIM * MINUS] = chi1_1 / 2. * stp[1]
                * (2. * skp[1] * sin(psin[1]) * stp[1]
                        + ctp[1] * (k[PLUS_] * cos(phin[1] + psin[1]) - k[MINUS] * cos(phin[1] - psin[1])));
        waveampcoeffs[PN10_B + AMPCOEFF_DIM * PLUS_] = 1. / 24. * ( //
                ckp[1] * skp[1] * stp[2] * ( //
                        4. * (15. * stp[2] + 51.) * cos(psin[2])
                                + 20. * skp[2] * (7. * stp[2] - 6.) * (4. * cos(psin[4]) - cos(psin[2])) //
                        ) + ( //
                        s2t * ( //
                                /*  */sin(phin[3] + psin[2]) * params->ccoeff10pn[7 + PN10DIM * PLUS_]
                                /**/+ sin(phin[3] - psin[2]) * params->ccoeff10pn[7 + PN10DIM * MINUS] //
                                ) + s2t * ( //
                                /*  */sin(phin[1] + psin[2]) * params->ccoeff10pn[8 + PN10DIM * PLUS_]
                                /**/+ sin(phin[1] - psin[2]) * params->ccoeff10pn[8 + PN10DIM * MINUS] //
                                ) + ctp[2] * s2t * 8. * ( //
                                /*  */sin(phin[3] + psin[4]) * params->ccoeff10pn[9 + PN10DIM * PLUS_]
                                /**/+ sin(phin[3] - psin[4]) * params->ccoeff10pn[9 + PN10DIM * MINUS] //
                                ) + 2. * skp[1] * ( //
                                /*  */cos(phin[2] + psin[2]) * params->ccoeff10pn[10 + PN10DIM * PLUS_]
                                /**/+ cos(phin[2] - psin[2]) * params->ccoeff10pn[10 + PN10DIM * MINUS] //
                                ) + 8. * skp[2] * s2t * (7. * stp[2] - 3.) * ( //
                                /*  */sin(phin[1] + psin[4]) * (+3. * k[PLUS_] + 4. * skp[2])
                                /**/+ sin(phin[1] - psin[4]) * (-3. * k[MINUS] + 4. * skp[2]) //
                                ) + 16. * skp[1] * (2. - 8. * stp[2] + 7. * stp[4]) * ( //
                                /*  */cos(phin[2] + psin[4]) * params->ccoeff10pn[11 + PN10DIM * PLUS_]
                                /**/+ cos(phin[2] - psin[4]) * params->ccoeff10pn[11 + PN10DIM * MINUS] //
                                ) + skp[1] * stp[2] * (stp[2] - 2.) * ( //
                                /*  */cos(phin[4] + psin[2]) * params->ccoeff10pn[11 + PN10DIM * PLUS_]
                                /**/+ cos(phin[4] - psin[2]) * params->ccoeff10pn[11 + PN10DIM * MINUS] - 8. * ( //
                                        /*  */cos(phin[4] + psin[4]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                                        /**/+ cos(phin[4] - psin[4]) * params->ccoeff05pn[1 + PN05DIM * MINUS] //
                                        )//
                                )//
                        )//
                );
        waveampcoeffs[PN10_B + AMPCOEFF_DIM * MINUS] = 1. / 12. * ( //
                30. * ctp[1] * skp[1] * sin(psin[2]) * stp[2] * (3. * skp[2] - 2.) + (stp[1] * ( //
                        /*  */cos(phin[1] + psin[2]) * params->ccoeff10pn[12 + PN10DIM * PLUS_]
                        /**/+ cos(phin[1] - psin[2]) * params->ccoeff10pn[12 + PN10DIM * MINUS] //
                        ) + 2. * ctp[1] * skp[1] * ( //
                        ( //
                        /*  */sin(phin[2] + psin[2]) * params->ccoeff10pn[13 + PN10DIM * PLUS_]
                        /**/+ sin(phin[2] - psin[2]) * params->ccoeff10pn[13 + PN10DIM * MINUS] //
                        ) + 4. * (7. * stp[2] - 2.) * ( //
                                /*  */sin(phin[2] + psin[4]) * params->ccoeff10pn[11 + PN10DIM * PLUS_]
                                /**/+ sin(phin[2] - psin[4]) * params->ccoeff10pn[11 + PN10DIM * MINUS] //
                                )//
                        ) + ctp[1] * skp[1] * stp[2] * ( //
                        ( //
                        /*  */sin(phin[4] + psin[2]) * params->ccoeff10pn[11 + PN10DIM * PLUS_]
                        /**/+ sin(phin[4] - psin[2]) * params->ccoeff10pn[11 + PN10DIM * MINUS] //
                        ) - 8. * ( //
                                /*  */sin(phin[4] + psin[4]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                                /**/+ sin(phin[4] - psin[4]) * params->ccoeff05pn[1 + PN05DIM * MINUS] //
                                )//
                        ) + stp[1] * ( //
                        ( //
                        /*  */cos(phin[3] + psin[2]) * params->ccoeff10pn[14 + PN10DIM * PLUS_]
                        /**/+ cos(phin[3] - psin[2]) * params->ccoeff10pn[14 + PN10DIM * MINUS] //
                        ) + 4. * (2. - 3. * stp[2]) * ( //
                                /*  */cos(phin[3] + psin[4]) * params->ccoeff10pn[9 + PN10DIM * PLUS_]
                                /**/+ cos(phin[3] - psin[4]) * params->ccoeff10pn[9 + PN10DIM * MINUS] //
                                )//
                        ) + 4. * skp[2] * stp[1] * (6. - 7. * stp[2]) * ( //
                        /*  */(-3. * k[PLUS_] - 4. * skp[2]) * cos(phin[1] + psin[4])
                        /**/+ (+3. * k[MINUS] - 4. * skp[2]) * cos(phin[1] - psin[4]) //
                        ))//
                );
        waveampcoeffs[PN10_A + AMPCOEFF_DIM * PLUS_] = 1. / 48. * (2. * skp[2] * stp[2] * ( //
                /*  */5. * skp[2] * (7. * stp[2] - 6.) * (cos(psin[2]) - 4. * cos(psin[4]))
                /**/- 2. * (15. * stp[2] + 51.) * cos(psin[2]) //
                ) + ( //
                /*  */16. * skp[3] * s2t * (7. * stp[2] - 3.) * ( //
                        /*  */sin(phin[1] + psin[4]) * k[PLUS_]
                        /**/+ sin(phin[1] - psin[4]) * k[MINUS] //
                        ) - (stp[2] - 2.) * skp[2] * stp[2] * ( //
                        /*  */cos(phin[4] + psin[2]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                        /**/+ cos(phin[4] - psin[2]) * params->ccoeff00pn[0 + PN00DIM * MINUS] //
                        ) + 4. * (stp[2] - 2.) * stp[2] * ( //
                        /*  */cos(phin[4] + psin[4]) * params->ccoeff10pn[0 + PN10DIM * PLUS_]
                        /**/+ cos(phin[4] - psin[4]) * params->ccoeff10pn[0 + PN10DIM * MINUS] //
                        ) + 2. * skp[1] * s2t * ( //
                        ( //
                        /*  */sin(phin[1] + psin[2]) * params->ccoeff10pn[1 + PN10DIM * PLUS_]
                        /**/+ sin(phin[1] - psin[2]) * params->ccoeff10pn[1 + PN10DIM * MINUS] //
                        ) + (
                        /*  */sin(phin[3] + psin[2]) * params->ccoeff10pn[2 + PN10DIM * PLUS_]
                        /**/+ sin(phin[3] - psin[2]) * params->ccoeff10pn[2 + PN10DIM * MINUS] //
                        ) - 8. * ctp[2] * (
                        /*  */sin(phin[3] + psin[4]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                        /**/+ sin(phin[3] - psin[4]) * params->ccoeff05pn[1 + PN05DIM * MINUS] //
                        )//
                        ) + 2. * ( //
                        /*  */cos(phin[2] + psin[2]) * params->ccoeff10pn[3 + PN10DIM * PLUS_]
                        /**/+ cos(phin[2] - psin[2]) * params->ccoeff10pn[3 + PN10DIM * MINUS] //
                        ) - 16. * skp[2] * (7. * stp[4] - 2. * (4. * stp[2] - 1.)) * ( //
                        /*  */cos(phin[2] + psin[4]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                        /**/+ cos(phin[2] - psin[4]) * params->ccoeff00pn[0 + PN00DIM * MINUS] //
                        )//
                ));
        waveampcoeffs[PN10_A + AMPCOEFF_DIM * MINUS] = 1. / 24. * ( //
                60. * ckp[1] * ctp[1] * skp[2] * sin(psin[2]) * stp[2] + ( //
                        8. * skp[3] * stp[1] * (7. * stp[2] - 6.) * ( //
                                /*  */cos(phin[1] + psin[4]) * k[PLUS_]
                                /**/+ cos(phin[1] - psin[4]) * k[MINUS] //
                                ) - ctp[1] * skp[2] * stp[2] * ( //
                                /*  */sin(phin[4] + psin[2]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                                /**/+ sin(phin[4] - psin[2]) * params->ccoeff00pn[0 + PN00DIM * MINUS] //
                                ) + 2. * skp[1] * stp[1] * ( //
                                /*  */cos(phin[1] + psin[2]) * params->ccoeff10pn[4 + PN10DIM * PLUS_]
                                /**/+ cos(phin[1] - psin[2]) * params->ccoeff10pn[4 + PN10DIM * MINUS] //
                                ) + 4. * ctp[1] * stp[2] * ( //
                                /*  */sin(phin[4] + psin[4]) * params->ccoeff10pn[0 + PN10DIM * PLUS_]
                                /**/+ sin(phin[4] - psin[4]) * params->ccoeff10pn[0 + PN10DIM * MINUS] //
                                ) + 2. * ctp[1] * ( //
                                /*  */sin(phin[2] + psin[2]) * params->ccoeff10pn[5 + PN10DIM * PLUS_]
                                /**/+ sin(phin[2] - psin[2]) * params->ccoeff10pn[5 + PN10DIM * MINUS] //
                                ) - 8. * ctp[1] * skp[2] * (7. * stp[2] - 2.) * ( //
                                /*  */sin(phin[2] + psin[4]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                                /**/+ sin(phin[2] - psin[4]) * params->ccoeff00pn[0 + PN00DIM * MINUS] //
                                ) - 8. * (2. - 3. * stp[2]) * skp[1] * stp[1] * ( //
                                /*  */cos(phin[3] + psin[4]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                                /**/+ cos(phin[3] - psin[4]) * params->ccoeff05pn[1 + PN05DIM * MINUS] //
                                ) + 2. * skp[1] * stp[1] * ( //
                                /*  */cos(phin[3] + psin[2]) * params->ccoeff10pn[6 + PN10DIM * PLUS_]
                                /**/+ cos(phin[3] - psin[2]) * params->ccoeff10pn[6 + PN10DIM * MINUS] //
                                )//
                        )//
                );
    }
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN05): {
        waveampcoeffs[PN05_B + AMPCOEFF_DIM * PLUS_] = 1. / 64. * (    //
                4. * ckp[1] * ctp[1] * stp[2] * (135. * cos(psin[3]) * skp[2] + cos(psin[1]) * (4. - 15. * skp[2])    //
                ) + (    //
                        2. * ctp[1] * (    //
                                9. * (2. - 3. * stp[2]) * (    //
                                        /*  */cos(phin[2] + psin[3]) * params->ccoeff05pn[5 + PN05DIM * PLUS_]
                                        /**/+ cos(phin[2] - psin[3]) * params->ccoeff05pn[5 + PN05DIM * MINUS]    //
                                        )//
                                /*    */+ cos(phin[2] + psin[1]) * params->ccoeff05pn[6 + PN05DIM * PLUS_]
                                        + cos(phin[2] - psin[1]) * params->ccoeff05pn[6 + PN05DIM * MINUS]    //
                                ) - skp[1] * stp[1] * (stp[2] - 2.) * (    //
                                27. * (    //
                                        /*  */sin(phin[3] + psin[3]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                                        /**/+ sin(phin[3] - psin[3]) * params->ccoeff00pn[0 + PN00DIM * MINUS] //
                                        )//
                                /**/+ sin(phin[3] + psin[1]) * params->ccoeff05pn[7 + PN05DIM * PLUS_]
                                /**/+ sin(phin[3] - psin[1]) * params->ccoeff05pn[7 + PN05DIM * MINUS]    //
                                ) + skp[1] * stp[1] * (    //
                                45. * (2. - 3. * stp[2]) * (    //
                                        /*  */sin(phin[1] + psin[3]) * params->ccoeff05pn[7 + PN05DIM * PLUS_]
                                        /**/+ sin(phin[1] - psin[3]) * params->ccoeff05pn[7 + PN05DIM * MINUS] //
                                        )//
                                /**/+ sin(phin[1] + psin[1]) * params->ccoeff05pn[8 + PN05DIM * PLUS_]
                                /**/+ sin(phin[1] - psin[1]) * params->ccoeff05pn[8 + PN05DIM * MINUS] //
                                ))//
                );
        waveampcoeffs[PN05_B + AMPCOEFF_DIM * MINUS] = 1. / 32. * (32. * sin(psin[1]) * stp[2] * c2k1 + (    //
                ctp[1] * skp[1] * stp[1] * (27. * (    //
                        /*  */cos(phin[3] + psin[3]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                        /**/+ cos(phin[3] - psin[3]) * params->ccoeff00pn[0 + PN00DIM * MINUS]    //
                        ) + (    //
                /**/cos(phin[3] + psin[1]) + 45. * cos(phin[1] + psin[3])    //
                ) * params->ccoeff05pn[7 + PN05DIM * PLUS_] + (    //
                /**/cos(phin[3] - psin[1]) + 45. * cos(phin[1] - psin[3])    //
                ) * params->ccoeff05pn[7 + PN05DIM * MINUS] + (    //
                        /*  */params->ccoeff05pn[9 + PN05DIM * PLUS_] * cos(phin[1] + psin[1])
                        /**/+ params->ccoeff05pn[9 + PN05DIM * MINUS] * cos(phin[1] - psin[1])    //
                        )) - (18. * c2t) * (    //
                        /*  */sin(phin[2] + psin[3]) * params->ccoeff05pn[5 + PN05DIM * PLUS_]
                        /**/+ sin(phin[2] - psin[3]) * params->ccoeff05pn[5 + PN05DIM * MINUS]    //
                        ) - 2. * (    //
                        /*  */sin(phin[2] + psin[1]) * params->ccoeff05pn[10 + PN05DIM * PLUS_]
                        /**/+ sin(phin[2] - psin[1]) * params->ccoeff05pn[10 + PN05DIM * MINUS]    //
                        )//
                ));
        waveampcoeffs[PN05_A + AMPCOEFF_DIM * PLUS_] = 1. / 64. * (4. * ctp[1] * skp[1] * stp[2] * (    //
                /**/-45. * skp[2] * cos(psin[3]) + cos(psin[1]) * (5. * skp[2] - 4.)    //
                ) - skp[2] * stp[1] * ((stp[2] - 2.) * (    //
                /*  */sin(phin[3] + psin[1]) * k[PLUS_]
                /**/+ sin(phin[3] - psin[1]) * k[MINUS]    //
                ) - 45. * (2. - 3. * stp[2]) * (    //
                /*  */sin(phin[1] + psin[3]) * k[PLUS_]
                /**/+ sin(phin[1] - psin[3]) * k[MINUS]    //
                )) + stp[1] * (    //
                /*  */sin(phin[1] + psin[1]) * params->ccoeff05pn[0 + PN05DIM * PLUS_]
                /**/+ sin(phin[1] - psin[1]) * params->ccoeff05pn[0 + PN05DIM * MINUS]    //
                - 9. * (stp[2] - 2.) * ( //
                        /*  */sin(phin[3] + psin[3]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                        /**/+ sin(phin[3] - psin[3]) * params->ccoeff05pn[1 + PN05DIM * MINUS] //
                        )//
                ) + 2. * ctp[1] * skp[1] * ( //
                /*  */cos(phin[2] + psin[1]) * params->ccoeff05pn[2 + PN05DIM * PLUS_]
                /**/+ cos(phin[2] - psin[1]) * params->ccoeff05pn[2 + PN05DIM * MINUS] //
                + 9. * (2. - 3. * stp[2]) * ( //
                        /*  */cos(phin[2] + psin[3]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                        /**/+ cos(phin[2] - psin[3]) * params->ccoeff00pn[0 + PN00DIM * MINUS] //
                        )//
                ));
        waveampcoeffs[PN05_A + AMPCOEFF_DIM * MINUS] = 1. / 32. * (    //
                -16. * s2k1 * stp[2] * sin(psin[1]) + ctp[1] * stp[1] * skp[2] * (    //
                        /*  */cos(phin[3] + psin[1]) * k[PLUS_]
                        /**/+ cos(phin[3] - psin[1]) * k[MINUS]    //
                        + 45. * (    //
                                /*  */cos(phin[1] + psin[3]) * k[PLUS_]
                                /**/+ cos(phin[1] - psin[3]) * k[MINUS]    //
                                )//
                        ) + 0.5 * s2t * (    //
                        /*  */cos(phin[1] + psin[1]) * params->ccoeff05pn[3 + PN05DIM * PLUS_]
                        /**/+ cos(phin[1] - psin[1]) * params->ccoeff05pn[3 + PN05DIM * MINUS] //
                        + 9. * (    //
                                /*  */cos(phin[3] + psin[3]) * params->ccoeff05pn[1 + PN05DIM * PLUS_]
                                /**/+ cos(phin[3] - psin[3]) * params->ccoeff05pn[1 + PN05DIM * MINUS]    //
                                )//
                        ) + 2. * skp[1] * (    //
                        /*  */sin(phin[2] + psin[1]) * params->ccoeff05pn[4 + PN05DIM * PLUS_]
                        /**/+ sin(phin[2] - psin[1]) * params->ccoeff05pn[4 + PN05DIM * MINUS] //
                        - 9. * c2t * (    //
                                /*  */sin(phin[2] + psin[3]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                                /**/+ sin(phin[2] - psin[3]) * params->ccoeff00pn[0 + PN00DIM * MINUS]    //
                                )//
                        )//
                );
    }
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN00): {
        waveampcoeffs[PN00_B + AMPCOEFF_DIM * PLUS_] = 0.5 * (s2t * (    //
                /*  */sin(phin[1] + psin[2]) * params->ccoeff00pn[1 + PN00DIM * PLUS_]
                /**/+ sin(phin[1] - psin[2]) * params->ccoeff00pn[1 + PN00DIM * MINUS]    //
                ) + skp[1] * (stp[2] - 2.) * (    //
                /*  */cos(phin[2] + psin[2]) * k[PLUS_]
                /**/+ cos(phin[2] - psin[2]) * k[MINUS]    //
                ) - 3. * s2k1 * stp[2] * cos(psin[2]));
        waveampcoeffs[PN00_B + AMPCOEFF_DIM * MINUS] = ctp[1] * skp[1] * (    //
                /*  */sin(phin[2] + psin[2]) * k[PLUS_]
                /**/+ sin(phin[2] - psin[2]) * k[MINUS]    //
                ) + stp[1] * (    //
                /*  */cos(phin[1] + psin[2]) * params->ccoeff00pn[1 + PN00DIM * PLUS_]
                /**/+ cos(phin[1] - psin[2]) * params->ccoeff00pn[1 + PN00DIM * MINUS]);

        waveampcoeffs[PN00_A + AMPCOEFF_DIM * PLUS_] = 0.25 * ((stp[2] - 2.) * (    //
                /*  */cos(phin[2] + psin[2]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                /**/+ cos(phin[2] - psin[2]) * params->ccoeff00pn[0 + PN00DIM * MINUS]    //
                ) - 2. * skp[1] * s2t * (    //
                /*  */sin(phin[1] + psin[2]) * k[PLUS_]
                /**/+ sin(phin[1] - psin[2]) * k[MINUS]    //
                ) + 6. * skp[2] * stp[2] * cos(psin[2]));
        waveampcoeffs[PN00_A + AMPCOEFF_DIM * MINUS] = 0.5 * (ctp[1] * (    //
                /*  */sin(phin[2] + psin[2]) * params->ccoeff00pn[0 + PN00DIM * PLUS_]
                /**/+ sin(phin[2] - psin[2]) * params->ccoeff00pn[0 + PN00DIM * MINUS]    //
                ) - 2. * stp[1] * skp[1] * (    //
                /*  */cos(phin[1] + psin[2]) * k[PLUS_]
                /**/+ cos(phin[1] - psin[2]) * k[MINUS]));

    }
    }


    REAL8 h[PLUS_CROSS_DIM] = { 0., 0. };
    switch (params->pnamp) {
    case (PNDEF): //default value -1 contains all corrections
    case (PN15):
        // Highest order, only leading order of eps(omega) is needed.
        for (int i = PLUS_; i < PLUS_CROSS_DIM; ++i) {
            h[i] += params->eps * epssqrt// eps_corr[PN00][PN15] * sqrt(eps_corr[PN00][PN15])
                    * (waveampcoeffs[PN15_A + AMPCOEFF_DIM * i] + waveampcoeffs[PN15_B + AMPCOEFF_DIM * i]
                            + waveampcoeffs[PN15_C + AMPCOEFF_DIM * i]);
        }
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN10):
        // Since highest order is 1.5 PN and there is no 0.5 PN order correction to eps(omega), leading order eps is enough.
        for (int i = PLUS_; i < PLUS_CROSS_DIM; ++i) {
            h[i] += params->eps // eps_corr[PN00][PN10]
                    * (4. * params->xi * waveampcoeffs[PN05_A + AMPCOEFF_DIM * i]
                    + waveampcoeffs[PN10_A + AMPCOEFF_DIM * i] + waveampcoeffs[PN10_C + AMPCOEFF_DIM * i]
                            + params->beta1 * waveampcoeffs[PN10_B + AMPCOEFF_DIM * i]
                            + params->beta1 * waveampcoeffs[PN10_D + AMPCOEFF_DIM * i]);
        }
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN05):
        //The 0.5 PN correction needs to include 1 PN correction of eps(omega) the amplitude is taken to 1.5 PN order
        for (int i = PLUS_; i < PLUS_CROSS_DIM; ++i) {
            REAL8 temp = waveampcoeffs[PN05_A + AMPCOEFF_DIM * i]
                    + params->beta1 * waveampcoeffs[PN05_B + AMPCOEFF_DIM * i]
                    - 2. * params->xi * waveampcoeffs[PN00_A + AMPCOEFF_DIM * i];
            if (params->pnamp == PN15) {
                h[i] += 1. / (params->eps * epssqrt) * (eps_corr[PN05][PN10] * sqrt(eps_corr[PN05][PN10])) * epssqrt * temp;
            } else {
                h[i] +=epssqrt * temp;
            }
        }
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN00):
        // If the amplitude is taken to 1 PN order, the eps(omega) needs to include the 1 PN correction, if amplitude is taken to 1.5 PN order, that eps(omega) needs to include corrections up to 1.5 PN order
        for (int i = PLUS_; i < PLUS_CROSS_DIM; ++i) {
            REAL8 temp = waveampcoeffs[PN00_A + AMPCOEFF_DIM * i]
                    + params->beta1 * waveampcoeffs[PN00_B + AMPCOEFF_DIM * i];
            if (params->pnamp == PN15) {
                h[i] += 1. / (params->eps * epssqrt) * eps_corr[PN00][PN15] * sqrt(eps_corr[PN00][PN15]) * temp;
            } else if (params->pnamp == PN10) {
                h[i] += 1. / (params->eps * epssqrt) * eps_corr[PN00][PN10] * sqrt(eps_corr[PN00][PN10]) * temp;
            } else {
                h[i] += /*                                                                               */temp;
            }
        }
        break;
    }
    REAL8 ampcoeff = 2. * G_CP2 * params->totalmass * params->eps * epssqrt * params->xi / params->dist;
    // rotating h+ and hx into the LALSimulation precessing radiation frame
    (*hplus)->data->data[idx] = ampcoeff * (h[PLUS_]*cos(params->polarizationangle)-h[MINUS]*sin(params->polarizationangle));
    (*hcross)->data->data[idx] = ampcoeff * (h[MINUS]*cos(params->polarizationangle)+h[PLUS_]*sin(params->polarizationangle));
    XLALFreeDmatrix(waveampcoeffs);
    return XLAL_SUCCESS;
}

/**
 * Interface routine, calculating the prefered variables for the Spin-dominated waveforms
 */
int XLALSimInspiralSpinDominatedWaveformInterfaceTD(REAL8TimeSeries **hplus, /**< +-polarization waveform */
REAL8TimeSeries **hcross, /**< x-polarization waveform */
REAL8 deltaT, /**< sampling interval (s) */
REAL8 m1, /**< mass of companion 1 (kg) */
REAL8 m2, /**< mass of companion 2 (kg) */
REAL8 fStart, /**< start GW frequency (Hz) */
REAL8 fRef, /**< end GW frequency (Hz) */
REAL8 D, /**< distance of source (m) */
REAL8 s1x, /**< initial value of S1x */
REAL8 s1y, /**< initial value of S1y */
REAL8 s1z, /**< initial value of S1z */
REAL8 lnhatx, /**< initial value of LNhatx */
REAL8 lnhaty, /**< initial value of LNhaty */
REAL8 lnhatz, /**< initial value of LNhatz */
REAL8 incl, /**< inclination, angle between L and line of sight N */
int phaseO, /**< twice PN phase order */
int amplitudeO, /**< twice PN amplitude order */
REAL8 phiRef /**< Reference phase at the Reference Frequency */
) {
    enum {
        X, Y, Z, DIM_XYZ,
    };
    REAL8 totalmass, nu, chi1, beta1, kappa1, totalJ, S1, omega, eta, romega, v, LN, theta, polarizationangle;
    REAL8 phin0;
    totalmass = m1 + m2;
    if (m1 > m2) {
        nu = m2 / m1;
    } else {
        nu = m1 / m2;
    }
    if (LAL_SDW_MAX_PN_PARAM < 100. * nu * nu) {
        XLALPrintError(
                "XLAL Error: Spin-dominated waveforms error: Please make sure that the total mass is higher than 45 solar mass, and mass ratio is lower than 0.03125. Also above 130 solar mass be aware that high starting frequency may result in termination right after start, due to high value of the pn parameter. \n");
        XLAL_ERROR(XLAL_EDOM);
    } //too high mass ratio for the waveform, abort
    omega = fStart * LAL_PI;
    eta = nu / (1. + nu) / (1. + nu);
    chi1 = sqrt(s1x * s1x + s1y * s1y + s1z * s1z);
    if (chi1 < 0.5) {
        XLALPrintError(
                "XLAL Error: Spin-dominated waveforms error: Please make sure that the dimensionless spin parameter is higher than 0.5 \n");
        XLAL_ERROR(XLAL_EDOM);
    }
    REAL8 Lnh[DIM_XYZ] = { lnhatx, lnhaty, lnhatz };
    REAL8 LNHdotS1 = (Lnh[X] * s1x + Lnh[Y] * s1y + Lnh[Z] * s1z) / chi1;
    if (LNHdotS1 - 1.0 > 0 && LNHdotS1 - 1.0 < 1.0e-12)
        {kappa1 = acos(1.0);
    } else {
        kappa1 = acos((Lnh[X] * s1x + Lnh[Y] * s1y + Lnh[Z] * s1z) / chi1);
        }

    // Calculate the orbital angular momentum, up to 1.5 PN, with SO corrections
    v = cbrt(G_CP2 * totalmass * omega / LAL_C_SI);
    romega = G_CP2 * totalmass / v / v * (1.);
    LN = eta * totalmass * romega * romega * omega;
    // Calculate Spin magnitude, and the total angular momentum J
    S1 = chi1 * LAL_G_SI / LAL_C_SI * totalmass * totalmass * eta / nu;

    REAL8 JxN[DIM_XYZ], JxL[DIM_XYZ], Nvec[DIM_XYZ];
    REAL8 J[DIM_XYZ] = { LN * Lnh[X] + S1 * s1x / chi1, LN * Lnh[Y] + S1 * s1y / chi1, LN * Lnh[Z] + S1 * s1z / chi1 };
    totalJ = sqrt(J[X] * J[X] + J[Y] * J[Y] + J[Z] * J[Z]);
    // calculate the remaining angles
    if (kappa1 < 1e-7) {
        kappa1 = 0.;
        phin0 = LAL_PI_2;
        beta1 = 0.;
        polarizationangle=LAL_PI;
        theta = incl;
    } else if (kappa1 - LAL_PI > 0 && kappa1 - LAL_PI < 1.0e-12) {
        kappa1 = LAL_PI;
        phin0 = LAL_PI_2;
        beta1 = LAL_PI;
        polarizationangle=LAL_PI;
        theta = incl;
    } else {
        REAL8 JdotS = (J[X] * s1x + J[Y] * s1y + J[Z] * s1z) / totalJ / chi1 ;
        if (JdotS - 1.0 > 0 && JdotS - 1.0 < 1.0e-12){
            beta1 = acos(1.0);
        } else{
        beta1 = acos((J[X] * s1x + J[Y] * s1y + J[Z] * s1z) / totalJ / chi1);
        }
        // Line of sight vectorProd
        Nvec[X]=0.;
        Nvec[Y]=cos(incl);
        Nvec[Z]=sin(incl);
        theta = acos((J[X] * Nvec[X] + J[Y] * Nvec[Y] + J[Z] * Nvec[Z]) / totalJ );



        vectorProd(J, Lnh, totalJ, JxL);
        vectorProd(J, Nvec, totalJ, JxN);
        REAL8 JxNxN[DIM_XYZ], NxLxN[DIM_XYZ], LxN[DIM_XYZ], NLNxJNN[DIM_XYZ];

        REAL8 JxNxNamp, NxLxNamp, polarizationanglesign, JxLamp;
        vectorProd(Lnh, Nvec, 1.0, LxN);
        vectorProd(JxN, Nvec, 1.0, JxNxN);
        vectorProd(Nvec, LxN, 1.0, NxLxN);

        JxLamp = sqrt(JxL[X]*JxL[X]+JxL[Y]*JxL[Y]+JxL[Z]*JxL[Z]);
        JxNxNamp=sqrt(JxNxN[X]*JxNxN[X]+JxNxN[Y]*JxNxN[Y]+JxNxN[Z]*JxNxN[Z]);
        NxLxNamp=sqrt(NxLxN[X]*NxLxN[X]+NxLxN[Y]*NxLxN[Y]+NxLxN[Z]*NxLxN[Z]);

        REAL8 JxNxNdotNxLxN = (JxNxN[X]*NxLxN[X]+JxNxN[Y]*NxLxN[Y]+JxNxN[Z]*NxLxN[Z])/JxNxNamp/NxLxNamp ;

        if (JxNxNdotNxLxN - 1.0 > 0 && JxNxNdotNxLxN - 1.0 < 1.0e-12){
            polarizationangle = acos(1.0);
        } else if (JxNxNdotNxLxN + 1.0 < 0 && fabs(JxNxNdotNxLxN + 1.0) < 1.0e-12){
            polarizationangle = acos(-1.0);
        }
        else {
        polarizationangle = acos((JxNxN[X]*NxLxN[X]+JxNxN[Y]*NxLxN[Y]+JxNxN[Z]*NxLxN[Z])/JxNxNamp/NxLxNamp);
        }

        vectorProd(NxLxN, JxNxN, JxNxNamp*NxLxNamp, NLNxJNN);

        if (beta1 < 1e-7){
            polarizationangle =LAL_PI;
        } else{
        polarizationanglesign = NLNxJNN[X]*Nvec[X]+NLNxJNN[Y]*Nvec[Y]+NLNxJNN[Z]*Nvec[Z];

        if (polarizationanglesign < 0.) {
            polarizationangle *= -1.0;
        }
        }

        REAL8 JxNxJ[DIM_XYZ];
        REAL8 alpha0;
        vectorProd(JxN, J, totalJ, JxNxJ);

        REAL8 JxNxJamp;
        JxNxJamp = sqrt(JxNxJ[X]*JxNxJ[X]+JxNxJ[Y]*JxNxJ[Y]+JxNxJ[Z]*JxNxJ[Z]);

        alpha0 = acos((JxNxJ[X] * JxL[X] + JxNxJ[Y] * JxL[Y] + JxNxJ[Z] * JxL[Z])/JxNxJamp/JxLamp);

        phin0 =  alpha0;

        REAL8 JNJxJL[DIM_XYZ];
        vectorProd(JxNxJ, JxL, JxNxJamp*JxLamp, JNJxJL);
        REAL8 JNJxJLsign = JNJxJL[X]*J[X]+JNJxJL[Y]*J[Y]+JNJxJL[Z]*J[Z];

        if (JNJxJLsign > 0.) {
            phin0 *= -1.0;
        }
        phiRef += LAL_PI_2;
    }
    // calling the SDW driver with the prefered variables
    int n = XLALSimInspiralSpinDominatedWaveformDriver(hplus, hcross, totalmass, nu, chi1, D, kappa1, beta1, theta,
            fStart, fRef, phaseO, amplitudeO, deltaT, phiRef, phin0, polarizationangle);
    return n;
}

/**
 * Function calculating the Spin-Dominated waveforms
 * This waveform is an inspiral only, 1 spin, precessing waveform.
 * For the formulae see the appendix of Arxiv:1209.1722
 */
int XLALSimInspiralSpinDominatedWaveformDriver(REAL8TimeSeries **hplus, /**< +-polarization waveform */
REAL8TimeSeries **hcross, /**< x-polarization waveform */
REAL8 totalmass, /**< total mass of the binary */
REAL8 nu, /**< mass ratio */
REAL8 chi1, /**< dimensionless spin paramter */
REAL8 D, /**< Distance to the source */
REAL8 kappa1, /**< Angle span by S_1 and L */
REAL8 beta1, /**< Angle span by J and S_1 */
REAL8 theta, /**< Angle span by the line of sight and J */
REAL8 fStart, /**< Starting gravitational wave frequency*/
REAL8 fRef, /**< Ending gravitational wave frequency*/
int phaseO, /**< twice PN phase order */
int amplitudeO, /**< twice PN amplitude order */
REAL8 deltaT, /**< Sampling time interval */
REAL8 phiRef, /**< Reference phase at the Reference Frequency */
REAL8 phin0, /**< Starting value of the phi_n parameter */
REAL8 polarizationangle /**< Angle to rotate the radiaton frame to the default LALSimulation radiation frame */
) {
    int idx;
    int n;
    unsigned int i;
    REAL8 phiShift;
    LIGOTimeGPS tStart = LIGOTIMEGPSZERO;
    /* check inputs for sanity */
    if (*hplus) {
        XLAL_ERROR(XLAL_EFAULT);
    }
    if (*hcross) {
        XLAL_ERROR(XLAL_EFAULT);
    }
    if (deltaT <= 0) {
        XLAL_ERROR(XLAL_EDOM);
    }
    if (totalmass < 0) {
        XLAL_ERROR(XLAL_EDOM);
    }
    if (fStart <= 0) {
        XLAL_ERROR(XLAL_EDOM);
    }
    /* set up the integrator*/
    LALAdaptiveRungeKuttaIntegrator *integrator = XLALAdaptiveRungeKutta4Init(LAL_SDW_NUM_VARIABLES,
            XLALSpinDominatedWaveformDerivatives, XLALSpinDominatedWaveformStoppingTest, LAL_SDW_ABSOLUTE_TOLERANCE,
            LAL_SDW_RELATIVE_TOLERANCE);
    if (!integrator) {
        XLALPrintError("XLAL Error - %s: Cannot allocate integrator\n", __func__);
        XLAL_ERROR(XLAL_EFUNC);
    }
    /* stop the integration only when the test is true */
    integrator->stopontestonly = 1;
    LALSDWaveformParams params;
    params.totalmass = totalmass;
    params.nu = nu;
    params.chi1 = chi1;
    params.dist = D;
    params.kappa1 = kappa1;
    params.beta1 = beta1;
    params.theta = theta;
    params.eps = 0.;
    params.xi = 0.;
    params.pnamp = amplitudeO;
    params.pnphase = phaseO;
    params.prevdomega = 0.;
    params.polarizationangle = polarizationangle;
    n = XLALSpinDominatedWaveformConstantCoefficients(&params);
    if (n < 0) {
        XLAL_ERROR(XLAL_EFUNC);
    }
    REAL8 yin[LAL_SDW_NUM_VARIABLES];
    yin[PHI] = phin0;
    yin[OMEGA] = fStart * LAL_PI;
    yin[PSI] = 0.;
    REAL8Array *yout;
    // estimating the length of the waveform
    REAL8 length = 5. / 256. * pow(fStart * LAL_PI, -8. / 3.) * (1. + params.nu) * (1. + params.nu) / params.nu
            * pow(G_CP2 * params.totalmass / LAL_C_SI, -5. / 3.);
    INT4 intLen = XLALAdaptiveRungeKutta4Hermite(integrator, (void *) &params, yin, 0.0, length, deltaT, &yout);
    UNUSED INT4 intReturn = integrator->returncode;
    XLALAdaptiveRungeKuttaFree(integrator);
    REAL8TimeSeries *phin = XLALCreateREAL8TimeSeries("PHI_N", &tStart, 0., deltaT, &lalDimensionlessUnit, intLen);
    REAL8TimeSeries *omega = XLALCreateREAL8TimeSeries("OMEGA", &tStart, 0., deltaT, &lalDimensionlessUnit, intLen);
    REAL8TimeSeries *psi = XLALCreateREAL8TimeSeries("ORBITAL_PHASE", &tStart, 0., deltaT, &lalDimensionlessUnit,
            intLen);
    for (idx = 0; idx < intLen; idx++) {
        phin->data->data[idx] = yout->data[intLen + idx];
        omega->data->data[idx] = yout->data[2 * intLen + idx];
        psi->data->data[idx] = yout->data[3 * intLen + idx];
    }

    if (fRef == 0) {
        phiShift = phiRef - psi->data->data[0];
        for (i = 0; i < psi->data->length; i++) {
            psi->data->data[i] += phiShift;
        }
    } else if (fRef == fStart) {
        phiShift = phiRef - psi->data->data[0];
        for (i = 0; i < psi->data->length; i++) {
            psi->data->data[i] += phiShift;
        }
    } else {
        XLALPrintError(
                "XLAL Error: Spin-dominated waveforms error: Please set the reference frequency as the starting frequency, Setting 0 will default to the starting frequency. \n");
        XLAL_ERROR(XLAL_EDOM);
    }
    if ((*hplus) && (*hcross)) {
        if ((*hplus)->data->length != (*hcross)->data->length) {
            XLALPrintError("***  h+ and hx differ in length\n");
            XLAL_ERROR(XLAL_EFAILED);
        } else {
            if ((int) (*hplus)->data->length < intLen) {
                XLALPrintError("*** ERROR: h+ and hx too short\n");
                XLAL_ERROR(XLAL_EFAILED);
            } else {
                XLALGPSAdd(&((*hplus)->epoch), -intLen * deltaT);
                XLALGPSAdd(&((*hcross)->epoch), -intLen * deltaT);
            }
        }
    } else {
        XLALGPSAdd(&tStart, -intLen * deltaT);
        *hplus = XLALCreateREAL8TimeSeries("H+", &tStart, 0.0, deltaT, &lalStrainUnit, intLen);
        *hcross = XLALCreateREAL8TimeSeries("Hx", &tStart, 0.0, deltaT, &lalStrainUnit, intLen);
        if (*hplus == NULL || *hcross == NULL) {
            XLAL_ERROR(XLAL_ENOMEM);
        }
    }


    REAL8 expr[3];
    for (idx = 0; idx < intLen; idx++) {
        expr[PHI] = phin->data->data[idx];
        expr[OMEGA] = omega->data->data[idx];
        expr[PSI] = psi->data->data[idx];

    n = XLALSpinDominatedWaveformBuild(&params, expr, hplus, hcross, idx);
        if (n < 0) {
            XLAL_ERROR(XLAL_EFUNC);
        }

    }

    XLALDestroyREAL8Array(yout);
    XLALDestroyREAL8TimeSeries(phin);
    XLALDestroyREAL8TimeSeries(omega);
    XLALDestroyREAL8TimeSeries(psi);
    return intLen;
}

/**
 * Function calculating the derivatives of the three time dependent variables of the Spin-Dominated waveforms (SDW)
 * The first paramter is phi_n, Eq 27 of Arxiv:1005.5330, taken for 1 spin case, and integrated over an orbital period.
 * The second parameter is omega, the derivative is taken from Arxiv: astro-ph/0504538, up to 2 PN orders with 1 spin. (In order to stay consistent with SDW)
 * The thirs parameter is the phase.
 */
static INT4 XLALSpinDominatedWaveformDerivatives(UNUSED REAL8 t, const REAL8 values[], REAL8 dvalues[], void *mparams) {
    LALSDWaveformParams *params = (LALSDWaveformParams *) mparams;
    // parameters required for the time derivatives
    REAL8 vP[OMEGA_POWER_DIM];
    vP[0] = 1.;
    vP[1] = cbrt(LAL_G_SI * params->totalmass * values[OMEGA] / LAL_C_SI / LAL_C_SI / LAL_C_SI);
    for (int i = 2; i < OMEGA_POWER_DIM; ++i) {
        vP[i] = vP[1] * vP[i - 1];
    }

    params->eps = vP[2];
    REAL8 epsP3 = params->eps * params->eps * params->eps;
    params->xi = params->nu / sqrt(params->eps);
    REAL8 eta = params->nu / (1. + params->nu) / (1. + params->nu);
    REAL8 phasecoeff = 96. / 5. * eta * vP[5] * values[OMEGA] * values[OMEGA];
    REAL8 sinKappa1 = sin(params->kappa1);
    REAL8 cosKappa1 = cos(params->kappa1);
    // Calculating the derivatives

    dvalues[PHI] = 0;
    if (params->kappa1 != 0 ) {
        switch (params->pnphase) {
        case (PNDEF): //default value -1 contains all corrections
        case (PN20):
            dvalues[PHI] += +3. / 2. / LAL_G_SI / params->totalmass / sinKappa1 * params->chi1 * params->chi1 * LAL_C_SI
                    * LAL_C_SI * LAL_C_SI * epsP3 * sqrt(params->eps)
                    * (1. - 2. * params->xi * sqrt(params->eps)) * (sinKappa1 * cosKappa1
                            + params->beta1 * cosKappa1 * cosKappa1);
#if __GNUC__ >= 7
            __attribute__ ((fallthrough)); /* no break */
#endif
        case (PN15):
            dvalues[PHI] += params->chi1 * LAL_C_SI * LAL_C_SI * LAL_C_SI * epsP3 / 2. / LAL_G_SI / params->totalmass
                    * ( 5. * sqrt(params->eps) * params->xi - 4.)*(1. + cosKappa1 * params->beta1 / sinKappa1 );
#if __GNUC__ >= 7
            __attribute__ ((fallthrough)); /* no break */
#endif
        case (PN10):
            /* no break */
        case (PN05):
            /* no break */
        case (PN00):
            break;
        }
    }

    dvalues[OMEGA] = 0;
    switch (params->pnphase) {
    case (PNDEF): //default value -1 contains all corrections
    case (PN20):
        dvalues[OMEGA] += phasecoeff * vP[PN20] * (34103. / 18144. + 13661. / 2016. * eta + 59. / 18. * eta * eta);
        // QM and Self-Spin components taken together
        dvalues[OMEGA] += phasecoeff * vP[PN20]
                * (5. / 2. * params->chi1 * params->chi1 * eta / params->nu * (3. * cosKappa1 * cosKappa1 - 1.)
                        + 1. / 96. * params->chi1 * params->chi1 * eta / params->nu * (6. + sinKappa1 * sinKappa1));
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN15):
        dvalues[OMEGA] += phasecoeff * vP[PN15] * (4. * LAL_PI);
        // SO component
        dvalues[OMEGA] += phasecoeff * vP[PN15]
                * (-1 / 12. * cosKappa1 * params->chi1 * (113. * eta / params->nu + 75. * eta));
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN10):
        dvalues[OMEGA] += -phasecoeff * vP[PN10] * (743. / 336. + 11. / 4. * eta);
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN05):
        /* no break */
    case (PN00):
        dvalues[OMEGA] += phasecoeff * 1.;
        break;
    }
    dvalues[PSI] = values[OMEGA]-dvalues[PHI]*cos(params->kappa1-params->beta1);

    return GSL_SUCCESS;
} /* end of XLALSpinDominatedWaveformDerivatives */

/**
 * Stopping test for the Spin-Dominated waveforms. Using MECO, or the desired ending frequency. The maximum value of the PN parameter is set to 0.8.
 */
static INT4 XLALSpinDominatedWaveformStoppingTest(UNUSED REAL8 t, const REAL8 values[], REAL8 dvalues[],
UNUSED void *mparams) {
    LALSDWaveformParams *params = (LALSDWaveformParams *) mparams;
    REAL8 vP[OMEGA_POWER_DIM];
    vP[0] = 1.;
    vP[1] = cbrt(G_CP2 * params->totalmass * values[OMEGA] / LAL_C_SI);
    for (int i = 2; i < OMEGA_POWER_DIM; ++i) {
        vP[i] = vP[1] * vP[i - 1];
    }
    REAL8 cosKappa1 = cos(params->kappa1);
    REAL8 eta = params->nu / (1. + params->nu) / (1. + params->nu);
    REAL8 omega = values[OMEGA];
    REAL8 domega = dvalues[OMEGA];
    REAL8 d2omega = dvalues[OMEGA] - params->prevdomega;
    params->prevdomega = dvalues[OMEGA];
    REAL8 mecotest;
    mecotest = 0.;

    switch (params->pnphase) {
    case (PNDEF): //default value -1 contains all corrections
    case (PN20):
        mecotest += +6. * vP[PN20]
                * (1. / 8. * (-27. + 19. * eta - eta * eta / 3.)
                        - (3. * cosKappa1 * cosKappa1 - 1.) / 2. * params->chi1 * params->chi1 * eta / params->nu);
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN15):
        mecotest += +5. * vP[PN15] * (8. / 3. * eta / params->nu + 2. * eta) * cosKappa1 * params->chi1;
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN10):
        mecotest += -4. * vP[PN10] * (3. + eta / 3.) / 4.;
#if __GNUC__ >= 7
        __attribute__ ((fallthrough)); /* no break */
#endif
    case (PN05):
        /* no break */
    case (PN00):
        mecotest += +2.;
        break;
    }

    /* Check value of SEOBNRv2 stopping frequency

    REAL8 mcheck1 = params->totalmass / (1. + params->nu);
    REAL8 mcheck2 = params->totalmass - mcheck1;
    REAL8 spincheckz1 = params->chi1 * cos(params->kappa1);
    REAL8 spincheckz2 = 0.0;
    REAL8 seobnr_stop_freq = XLALSimIMRSpinAlignedEOBPeakFrequency(mcheck1, mcheck2, spincheckz1, spincheckz2, 2);
*/

    if (mecotest < 0) {
        XLALPrintWarning("** LALSimInspiralSDW WARNING **: MECO reached\n");
        return -1;

    }  /*else if (omega > seobnr_stop_freq*LAL_PI) {
        XLALPrintWarning("** LALSimInspiralSDW WARNING **: SEOBNR stopping frequency reached\n");
        return -1;

    } */else if (isnan(omega)) {
        XLALPrintWarning("** LALSimInspiralSDW WARNING **: omega is NAN\n");
        return -1;
    } else if (vP[1] >= 1.) {
        XLALPrintWarning("** LALSimInspiralSDW WARNING **: PN parameter is too large\n");
        return -1;
    } else if (domega < 0.0) {
        XLALPrintWarning("** LALSimInspiralSDW WARNING **: domega < 0\n");
        return -1;
    } else if (d2omega <= 0.) {
        XLALPrintWarning("** LALSimInspiralSDW WARNING **: d2omega < 0\n");
        return -1;
    }
    return GSL_SUCCESS;
} /* End of XLALSpinDominatedWaveformStoppingTest */
