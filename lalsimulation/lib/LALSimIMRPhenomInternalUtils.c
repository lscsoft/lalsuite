/*
 *  Copyright (C) 2017 Sebastian Khan
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

/**
 * \author Sebastian Khan
 *
 * \file
 *
 * \brief Internal (not SWIG'd) Auxiliary functions for phenomenological model development
 *
 * Helper functions for phenomenological waveform models
 * Cannot be used through python SWIG wrapping.
 * NOTE: The convention for naming functions in there is to use
 * the prefix 'PhenomInternal_'
 * This is to avoid the use of defining functions as static and
 * therefore not able to be used by other modules by including
 * LALSimIMRPhenomInternalUtils.h.
 */

#include "LALSimIMRPhenomInternalUtils.h"

// /**
//  * Example how to write an internal phenom function
//  */
// void PhenomInternal_UtilsTest(){
//     printf("Hello! I am the PhenomInternal_UtilsTest function\n");
// }

// This function determines whether x and y are approximately equal to a relative accuracy epsilon.
// Note that x and y are compared to relative accuracy, so this function is not suitable for testing whether a value is approximately zero.

/**
 * This function determines whether x and y are approximately equal to a
 * relative accuracy epsilon.
 * Note that x and y are compared to relative accuracy, so this function is not
 * suitable for testing whether a value is approximately zero.
 */
bool PhenomInternal_approx_equal(REAL8 x, REAL8 y, REAL8 epsilon)
{
    return !gsl_fcmp(x, y, epsilon);
}

/**
 * If x and X are approximately equal to relative accuracy epsilon
 * then set x = X.
 * If X = 0 then use an absolute comparison.
 */
void PhenomInternal_nudge(REAL8 *x, REAL8 X, REAL8 epsilon)
{
    if (X != 0.0)
    {
        if (PhenomInternal_approx_equal(*x, X, epsilon))
        {
            XLAL_PRINT_INFO("Nudging value %.15g to %.15g\n", *x, X);
            *x = X;
        }
    }
    else
    {
        if (fabs(*x - X) < epsilon)
            *x = X;
    }
}

/**
 * Return the closest higher power of 2
 */
size_t PhenomInternal_NextPow2(const size_t n)
{
    // use pow here, not bit-wise shift, as the latter seems to run against an upper cutoff long before SIZE_MAX, at least on some platforms
    return (size_t)pow(2, ceil(log2(n)));
}

/**
 * Given m1 with aligned-spin chi1z and m2 with aligned-spin chi2z.
 * Enforce that m1 >= m2 and swap spins accordingly.
 * Enforce that the primary object (heavier) is indexed by 1.
 * To be used with aligned-spin waveform models.
 * There is another function for precessing waveform models.
 * This is currently not used
 */
int PhenomInternal_AlignedSpinEnforcePrimaryIsm1(
    REAL8 *m1,    /**< [out] mass of body 1 */
    REAL8 *m2,    /**< [out] mass of body 2 */
    REAL8 *chi1z, /**< [out] aligned-spin component of body 1 */
    REAL8 *chi2z  /**< [out] aligned-spin component of body 2 */
)
{
    REAL8 chi1z_tmp, chi2z_tmp, m1_tmp, m2_tmp;
    if (*m1 > *m2)
    {
        chi1z_tmp = *chi1z;
        chi2z_tmp = *chi2z;
        m1_tmp = *m1;
        m2_tmp = *m2;
    }
    else
    { /* swap spins and masses */
        chi1z_tmp = *chi2z;
        chi2z_tmp = *chi1z;
        m1_tmp = *m2;
        m2_tmp = *m1;
    }
    *m1 = m1_tmp;
    *m2 = m2_tmp;
    *chi1z = chi1z_tmp;
    *chi2z = chi2z_tmp;

    if (*m1 < *m2)
        XLAL_ERROR(XLAL_EDOM, "XLAL_ERROR in EnforcePrimaryIsm1. When trying\
 to enfore that m1 should be the larger mass.\
 After trying to enforce this m1 = %f and m2 = %f\n", *m1, *m2);

    return XLAL_SUCCESS;
}

/**
 * Given m1 with spins (chi1x, chi1y, chi1z) and m2 with spins (chi2x,chi2y,chi2z).
 * Enforce that m1 >= m2 and swap spins accordingly.
 * Enforce that the primary object (heavier) is indexed by 1.
 * To be used with precessing-spin waveform models.
 */
int PhenomInternal_PrecessingSpinEnforcePrimaryIsm1(
    REAL8 *m1,    /**< [out] mass of body 1 */
    REAL8 *m2,    /**< [out] mass of body 2 */
    REAL8 *chi1x, /**< [out] x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi1y, /**< [out] y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi1z, /**< [out] z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi2x, /**< [out] x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi2y, /**< [out] y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi2z  /**< [out] z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
)
{
    REAL8 m1_tmp, m2_tmp;
    REAL8 chi1x_tmp, chi1y_tmp, chi1z_tmp;
    REAL8 chi2x_tmp, chi2y_tmp, chi2z_tmp;
    if (*m1 > *m2)
    {
        chi1x_tmp = *chi1x;
        chi1y_tmp = *chi1y;
        chi1z_tmp = *chi1z;

        chi2x_tmp = *chi2x;
        chi2y_tmp = *chi2y;
        chi2z_tmp = *chi2z;

        m1_tmp = *m1;
        m2_tmp = *m2;
    }
    else
    { /* swap spins and masses */
        chi1x_tmp = *chi2x;
        chi1y_tmp = *chi2y;
        chi1z_tmp = *chi2z;

        chi2x_tmp = *chi1x;
        chi2y_tmp = *chi1y;
        chi2z_tmp = *chi1z;

        m1_tmp = *m2;
        m2_tmp = *m1;
    }
    *m1 = m1_tmp;
    *m2 = m2_tmp;
    *chi1x = chi1x_tmp;
    *chi1y = chi1y_tmp;
    *chi1z = chi1z_tmp;

    *chi2x = chi2x_tmp;
    *chi2y = chi2y_tmp;
    *chi2z = chi2z_tmp;

    if (*m1 < *m2)
        XLAL_ERROR(XLAL_EDOM, "XLAL_ERROR in EnforcePrimaryIsm1. When trying\
 to enfore that m1 should be the larger mass.\
 After trying to enforce this m1 = %f and m2 = %f\n",
                   *m1, *m2);

    return XLAL_SUCCESS;
}

/**
 * atan2 wrapper that returns 0 when both magnitudes of x and y are
 * below tol, otherwise it returns atan2(x, y)
 */
REAL8 PhenomInternal_atan2tol(REAL8 a, REAL8 b, REAL8 tol)
{
  REAL8 c;
  if (fabs(a) < tol && fabs(b) < tol)
    c = 0.;
  else
    c = atan2(a, b);
  return c;
}


/**
 * function to convert from 3d cartesian components to polar angles and vector magnitude.
 * https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates
 */
void PhenomInternal_ComputeIMRPhenomPv3CartesianToPolar(
    REAL8 *polar,
    REAL8 *azimuthal,
    REAL8 *magnitude,
    REAL8 x,
    REAL8 y,
    REAL8 z
)
{
    /* TODO: check that this is the correct convention */
    *magnitude = sqrt( x*x + y*y + z*z );
    if (PhenomInternal_approx_equal(*magnitude, 0, 1e-9)){
        *polar = 0.;
        *azimuthal = 0.;
    } else {
        *polar = acos( z / *magnitude );
        *azimuthal = PhenomInternal_atan2tol( y, x, PhenomInternal_MAX_TOL_ATAN );
    }
}


/**
 * Wrapper function for XLALOrbitalAngMom3PNSpinning.
 * Used to compute the orbital angular momentum at 3PN order
 * at a single frequency.
 * We assume that Lhat = (0,0,1)
 */
double PhenomInternal_OrbAngMom3PN(
    const double f_orb_hz,   /**< Orbtial frequency (Hz)  */
    const double m1,         /**< mass of primary in SI (kg) */
    const double m2,         /**< mass of secondary in SI (kg) */
    const double s1x,        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const double s1y,        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const double s1z,        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    const double s2x,        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const double s2y,        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const double s2z,        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const double f_0,        /**< Reference Gravitational Wave frequency (Hz) */
    const int ExpansionOrder /**< Keep terms upto ExpansionOrder in precession angles phi_z and zeta (1,2,3,4,5 or -1 for all orders) */
)
{

    const double mul = 1.0; /* Cosine of Polar angle of orbital angular momentum */
    const double phl = 0.0; /* Azimuthal angle of orbital angular momentum */
    double mu1 = 0.0;
    double ph1 = 0.0;
    double ch1 = 0.0;
    double mu2 = 0.0;
    double ph2 = 0.0;
    double ch2 = 0.0;

    REAL8Sequence *L_norm_3PN_Seq = XLALCreateREAL8Sequence(1);
    REAL8Sequence *freqs_seq = XLALCreateREAL8Sequence(1);

    L_norm_3PN_Seq->data[0] = 0.;
    freqs_seq->data[0] = f_orb_hz;

    PhenomInternal_ComputeIMRPhenomPv3CartesianToPolar(&mu1, &ph1, &ch1, s1x, s1y, s1z);
    PhenomInternal_ComputeIMRPhenomPv3CartesianToPolar(&mu2, &ph2, &ch2, s2x, s2y, s2z);

    int retcode = XLALOrbitalAngMom3PNSpinning(
        L_norm_3PN_Seq, freqs_seq,
        m1, m2,
        mul, phl,
        cos(mu1), ph1, ch1,
        cos(mu2), ph2, ch2,
        f_0, ExpansionOrder);
    XLAL_CHECK(retcode == XLAL_SUCCESS, XLAL_EFUNC, "XLALOrbitalAngMom3PNSpinning failed.");

    double L_norm_3PN = L_norm_3PN_Seq->data[0];

    XLALDestroyREAL8Sequence(L_norm_3PN_Seq);
    XLALDestroyREAL8Sequence(freqs_seq);

    return L_norm_3PN;
}

/**
 * Formula to predict the total radiated energy. Equation 3.7 and 3.8 arXiv:1508.07250
 * Input parameter s defined around Equation 3.7 and 3.8.
 */
static double EradRational0815_s(double eta, double s)
{
    double eta2 = eta * eta;
    double eta3 = eta2 * eta;

    return (eta * (0.055974469826360077 + 0.5809510763115132 * eta - 0.9606726679372312 * eta2 + 3.352411249771192 * eta3) *
            (1. + (-0.0030302335878845507 - 2.0066110851351073 * eta + 7.7050567802399215 * eta2) * s)) /
           (1. + (-0.6714403054720589 - 1.4756929437702908 * eta + 7.304676214885011 * eta2) * s);
}

/**
 * Wrapper function for EradRational0815_s.
 */
double PhenomInternal_EradRational0815(double eta, double chi1, double chi2)
{
    // Convention m1 >= m2
    double Seta = sqrt(1.0 - 4.0 * eta);
    double m1 = 0.5 * (1.0 + Seta);
    double m2 = 0.5 * (1.0 - Seta);
    double m1s = m1 * m1;
    double m2s = m2 * m2;
    // arXiv:1508.07250
    double s = (m1s * chi1 + m2s * chi2) / (m1s + m2s);

    return EradRational0815_s(eta, s);
}

/**
 * helper function to multiple hlm with Ylm.
 * Adapted from LALSimIMREOBNRv2HMROMUtilities.c
 */
int PhenomInternal_IMRPhenomHMFDAddMode(
    COMPLEX16FrequencySeries *hptilde,
    COMPLEX16FrequencySeries *hctilde,
    COMPLEX16FrequencySeries *hlmtilde,
    REAL8 theta,
    REAL8 phi,
    INT4 l,
    INT4 m,
    INT4 sym)
{
    COMPLEX16 Y;
    UINT4 j;
    COMPLEX16 hlm; /* helper variable that contain a single point of hlmtilde */

    INT4 minus1l; /* (-1)^l */
    if (l % 2)
        minus1l = -1;
    else
        minus1l = 1;
    if (sym)
    { /* Equatorial symmetry: add in -m mode */
        Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
        COMPLEX16 Ymstar = conj(XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m));
        COMPLEX16 factorp = 0.5 * (Y + minus1l * Ymstar);
        COMPLEX16 factorc = -I * 0.5 * (Y - minus1l * Ymstar);
        for (j = 0; j < hlmtilde->data->length; ++j)
        {
            hlm = (hlmtilde->data->data[j]);
            hptilde->data->data[j] += factorp * hlm;
            hctilde->data->data[j] += factorc * hlm;
        }
    }
    else
    { /* not adding in the -m mode */
        Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
        COMPLEX16 factorp = 0.5 * Y;
        COMPLEX16 factorc = -I * factorp;
        for (j = 0; j < hlmtilde->data->length; ++j)
        {
            hlm = (hlmtilde->data->data[j]);
            hptilde->data->data[j] += factorp * hlm;
            hctilde->data->data[j] += factorc * hlm;
        }
    }

    return XLAL_SUCCESS;
}
