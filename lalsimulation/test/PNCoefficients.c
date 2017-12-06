/*
 *  Copyright (C) 2014 Andrew Lundgren
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

#include <stdlib.h>
#include <math.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include "LALSimInspiralPNCoefficients.c"

#define EPSILON 1.e-10

static int compare(
    REAL8 val1,
    REAL8 val2,
    int v_order,
    int log_order)
{
    if (fabs(val1 - val2) > EPSILON)
    {
        if (log_order == 0) { fprintf(stderr, "FAILED at %.1f PN order: %.10f versus %.10f\n", 0.5*v_order, val1, val2); }
        else { fprintf(stderr, "FAILED at %.1f PN order, in log^%u term: %.10f versus %.10f\n", 0.5*v_order, log_order, val1, val2); }
        return 1;
    }
    else
    {
        return 0;
    }
}

static int compare_dtdv(
    PNPhasingSeries *dtdv1,
    PNPhasingSeries *dtdv2)
{
    int ret = 0;

    ret += compare(dtdv1->v[0], dtdv2->v[0], 0, 0);
    ret += compare(dtdv1->v[2], dtdv2->v[2], 2, 0);
    ret += compare(dtdv1->v[3], dtdv2->v[3], 3, 0);
    ret += compare(dtdv1->v[4], dtdv2->v[4], 4, 0);
    ret += compare(dtdv1->v[5], dtdv2->v[5], 5, 0);
    ret += compare(dtdv1->v[6], dtdv2->v[6], 6, 0);
    ret += compare(dtdv1->vlogv[6], dtdv2->vlogv[6], 6, 1);
    ret += compare(dtdv1->v[7], dtdv2->v[7], 7, 0);

    return ret;
}

/* *
 * The factor to go from dtdv to phasing.
 * Derived from the stationary phase approximation,
 * Psi(f) = 2 Phi_orb - 2 Pi f t,
 * where f = v^3 / Pi,
 * Phi_orb = int v^3 dt/dv dv,
 *       t = int dt / dv dv
 * */
static REAL8 pfac(int n) { return 40./(n-5.)/(n-8.); }

/* The factor to go from a dt/dv log term to non-log phasing */
static REAL8 logpfac(int n)
{
    REAL8 temp = (n-5.)*(n-8.);
    return 40.*(13.-2.*n)/temp/temp;
}

static int compare_phasing_to_dtdv(
    PNPhasingSeries *phasing,
    PNPhasingSeries *dtdv)
{
    int ret = 0;

    ret += compare(phasing->v[2], pfac(2)*dtdv->v[2], 2, 0);
    ret += compare(phasing->v[3], pfac(3)*dtdv->v[3], 3, 0);
    ret += compare(phasing->v[4], pfac(4)*dtdv->v[4], 4, 0);
    /* The 2.5 pN term is anomalous - it integrates to a log term */
    ret += compare(phasing->v[5], (-40./9.)*dtdv->v[5], 5, 0);
    ret += compare(phasing->vlogv[5], (-40./3.)*dtdv->v[5], 5, 1);
    /* The 3 pN term has an extra piece coming from the log */
    ret += compare(phasing->v[6], pfac(6)*dtdv->v[6] + logpfac(6)*dtdv->vlogv[6], 6, 0);
    ret += compare(phasing->vlogv[6], pfac(6)*dtdv->vlogv[6], 6, 1);
    ret += compare(phasing->v[7], pfac(7)*dtdv->v[7], 7, 0);

    return ret;
}

/* Helper function to calculate sum of flux[j]*dtdv[k-j] from 1 to k-1 */
static REAL8 sum(
    REAL8 *arr1,
    REAL8 *arr2,
    int k)
{
    REAL8 accum = 0.;
    for (int j=1; j < k; j++)
    {
        accum += arr1[j]*arr2[k-j];
    }
    return accum;
}

/* The TaylorT2/F2 waveform is defined by
 * dt/dv = trunc[ -(dE/dv) / flux(v) ]
 * where trunc[ ] means truncate the expression at the desired PN order.
 * There's a recursive expression for dt/dv at each order. Write
 * dt/dv = T_0 ( 1 + \sum_{i=1}^n T_i v^i )
 * e(v) = E_0 ( 1 + \sum_{i=1}^n E_i v^i )
 * flux(v) = F_0 ( 1 + \sum_{i=1}^n F_i v^i )
 *
 * Then
 *
 * T_i = (1+i/2)*E_i - F_i - \sum_{j=1}^{i-1} F_j T_{i-j}
 *
 * This doesn't handle terms involving log(v), so we'll treat them separately.
 * Be careful with the dE/dv, which causes the log and non-log terms to mix.
 * Also, because spin-spin PN terms are only known to leading order, we have
 * to be careful with spin terms. We first add in the spin-orbit terms, but
 * only keeping terms linear in spin. Then we add in the spin-spin terms at 2 PN.
 */

static void dtdv_from_energy_flux(
    PNPhasingSeries *dtdv,
    const REAL8 m1M,
    const REAL8 chi1,
    const REAL8 chi2,
    const REAL8 qm_def1,
    const REAL8 qm_def2
    )
{
    REAL8 m2M = 1.-m1M;
    REAL8 eta = m1M*m2M;
    /* Spins use the wacky lal convention */
    REAL8 S1L = m1M*m1M*chi1;
    REAL8 S2L = m2M*m2M*chi2;

    REAL8 energy[9];
    REAL8 flux[9];

    energy[1] = 0.;
    energy[2] = XLALSimInspiralPNEnergy_2PNCoeff(eta);
    energy[3] = XLALSimInspiralPNEnergy_3PNSOCoeff(m1M)*S1L + XLALSimInspiralPNEnergy_3PNSOCoeff(m2M)*S2L;
    energy[4] = XLALSimInspiralPNEnergy_4PNCoeff(eta);
    energy[5] = XLALSimInspiralPNEnergy_5PNSOCoeff(m1M)*S1L + XLALSimInspiralPNEnergy_5PNSOCoeff(m2M)*S2L;
    energy[6] = XLALSimInspiralPNEnergy_6PNCoeff(eta);
    energy[7] = XLALSimInspiralPNEnergy_7PNSOCoeff(m1M)*S1L + XLALSimInspiralPNEnergy_7PNSOCoeff(m2M)*S2L;

    flux[1] = 0.;
    flux[2] = XLALSimInspiralPNFlux_2PNCoeff(eta);
    flux[3] = XLALSimInspiralPNFlux_3PNCoeff(eta);
    flux[4] = XLALSimInspiralPNFlux_4PNCoeff(eta);
    flux[5] = XLALSimInspiralPNFlux_5PNCoeff(eta);
    flux[6] = XLALSimInspiralPNFlux_6PNCoeff(eta);
    flux[7] = XLALSimInspiralPNFlux_7PNCoeff(eta);
    /* Add the spin-orbit fluxes */
    flux[3] += XLALSimInspiralPNFlux_3PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_3PNSOCoeff(m2M)*S2L;
    flux[5] += XLALSimInspiralPNFlux_5PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_5PNSOCoeff(m2M)*S2L;
    flux[6] += XLALSimInspiralPNFlux_6PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_6PNSOCoeff(m2M)*S2L;
    flux[7] += XLALSimInspiralPNFlux_7PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_7PNSOCoeff(m2M)*S2L;

    /* FIXME: Change this when the log convention is fixed in flux
     Needed because the flux term currently multiplies log(16v^2) */
    REAL8 flux6l = 2.*XLALSimInspiralPNFlux_6PNLogCoeff(eta);
    flux[6] += XLALSimInspiralPNFlux_6PNLogCoeff(eta)*log(16.);

    memset(dtdv, 0, sizeof(PNPhasingSeries));
    dtdv->v[0] = -2.*XLALSimInspiralPNEnergy_0PNCoeff(eta) / XLALSimInspiralPNFlux_0PNCoeff(eta);
    dtdv->v[1] = 0.; /* there's no 0.5 PN term */
    for (int i = 2; i < 8; i++)
    {
        dtdv->v[i] = (1.+0.5*i)*energy[i] - flux[i] - sum(flux, dtdv->v, i);
    }
    dtdv->vlogv[6] = -flux6l;

    /* Remove the (spin-orbit)^2 term from the 3 PN term, because spin^2 is incomplete above 2 PN */
    REAL8 energy_so3 = XLALSimInspiralPNEnergy_3PNSOCoeff(m1M)*S1L + XLALSimInspiralPNEnergy_3PNSOCoeff(m2M)*S2L;
    REAL8 flux_so3 = XLALSimInspiralPNFlux_3PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_3PNSOCoeff(m2M)*S2L;
    dtdv->v[6] -= -flux_so3 * ( (5./2.)*energy_so3 - flux_so3 );

    /* Calculate the leading-order spin-spin terms separately
     * FIXME: For now, only do aligned spins
     */
    REAL8 energy_ss4 = XLALSimInspiralPNEnergy_4PNS1S2LCoeff(eta)*S1L*S2L;
    energy_ss4 += XLALSimInspiralPNEnergy_4PNS1S2Coeff(eta)*S1L*S2L;
    energy_ss4 += qm_def1*XLALSimInspiralPNEnergy_4PNQM2SCoeff(m1M)*S1L*S1L;
    energy_ss4 += qm_def1*XLALSimInspiralPNEnergy_4PNQM2SOCoeff(m1M)*S1L*S1L;
    energy_ss4 += qm_def2*XLALSimInspiralPNEnergy_4PNQM2SCoeff(m2M)*S2L*S2L;
    energy_ss4 += qm_def2*XLALSimInspiralPNEnergy_4PNQM2SOCoeff(m2M)*S2L*S2L;

    REAL8 flux_ss4 = XLALSimInspiralPNFlux_4PNS1S2LCoeff(eta)*S1L*S2L;
    flux_ss4 += XLALSimInspiralPNFlux_4PNS1S2Coeff(eta)*S1L*S2L;

    flux_ss4 += qm_def1*XLALSimInspiralPNFlux_4PNQM2SCoeff(m1M)*S1L*S1L;
    flux_ss4 += qm_def1*XLALSimInspiralPNFlux_4PNQM2SOCoeff(m1M)*S1L*S1L;
    flux_ss4 += XLALSimInspiralPNFlux_4PNSelf2SCoeff(m1M)*S1L*S1L;
    flux_ss4 += XLALSimInspiralPNFlux_4PNSelf2SOCoeff(m1M)*S1L*S1L;

    flux_ss4 += qm_def2*XLALSimInspiralPNFlux_4PNQM2SCoeff(m2M)*S2L*S2L;
    flux_ss4 += qm_def2*XLALSimInspiralPNFlux_4PNQM2SOCoeff(m2M)*S2L*S2L;
    flux_ss4 += XLALSimInspiralPNFlux_4PNSelf2SCoeff(m2M)*S2L*S2L;
    flux_ss4 += XLALSimInspiralPNFlux_4PNSelf2SOCoeff(m2M)*S2L*S2L;

    dtdv->v[4] += 3.*energy_ss4 - flux_ss4;

    return;
}

static void dtdv_from_pncoefficients(
    PNPhasingSeries *dtdv,
    const REAL8 m1M,
    const REAL8 chi1,
    const REAL8 chi2,
    const REAL8 qm_def1,
    const REAL8 qm_def2
    )
{
    REAL8 m2M = 1.-m1M;
    REAL8 eta = m1M*m2M;
    /* Spins use the wacky lal convention */
    REAL8 S1L = m1M*m1M*chi1;
    REAL8 S2L = m2M*m2M*chi2;

    memset(dtdv, 0, sizeof(PNPhasingSeries));
    dtdv->v[0] = XLALSimInspiralTaylorT2dtdv_0PNCoeff(eta);
    dtdv->v[1] = 0.;
    dtdv->v[2] = XLALSimInspiralTaylorT2dtdv_2PNCoeff(eta);
    dtdv->v[3] = XLALSimInspiralTaylorT2dtdv_3PNCoeff(eta);
    dtdv->v[4] = XLALSimInspiralTaylorT2dtdv_4PNCoeff(eta);
    dtdv->v[5] = XLALSimInspiralTaylorT2dtdv_5PNCoeff(eta);
    dtdv->v[6] = XLALSimInspiralTaylorT2dtdv_6PNCoeff(eta);
    dtdv->vlogv[6] = XLALSimInspiralTaylorT2dtdv_6PNLogCoeff(eta);
    dtdv->v[7] = XLALSimInspiralTaylorT2dtdv_7PNCoeff(eta);

    dtdv->v[3] += XLALSimInspiralTaylorT2dtdv_3PNSOCoeff(m1M)*S1L + XLALSimInspiralTaylorT2dtdv_3PNSOCoeff(m2M)*S2L;
    dtdv->v[5] += XLALSimInspiralTaylorT2dtdv_5PNSOCoeff(m1M)*S1L + XLALSimInspiralTaylorT2dtdv_5PNSOCoeff(m2M)*S2L;
    dtdv->v[6] += XLALSimInspiralTaylorT2dtdv_6PNSOCoeff(m1M)*S1L + XLALSimInspiralTaylorT2dtdv_6PNSOCoeff(m2M)*S2L;
    dtdv->v[7] += XLALSimInspiralTaylorT2dtdv_7PNSOCoeff(m1M)*S1L + XLALSimInspiralTaylorT2dtdv_7PNSOCoeff(m2M)*S2L;

    dtdv->v[4] += XLALSimInspiralTaylorT2dtdv_4PNS1S2Coeff(eta)*S1L*S2L;
    dtdv->v[4] += XLALSimInspiralTaylorT2dtdv_4PNS1S2LCoeff(eta)*S1L*S2L;
    dtdv->v[4] += XLALSimInspiralTaylorT2dtdv_4PNSSCoeff(m1M, qm_def1)*S1L*S1L;
    dtdv->v[4] += XLALSimInspiralTaylorT2dtdv_4PNSSLCoeff(m1M, qm_def1)*S1L*S1L;
    dtdv->v[4] += XLALSimInspiralTaylorT2dtdv_4PNSSCoeff(m2M, qm_def2)*S2L*S2L;
    dtdv->v[4] += XLALSimInspiralTaylorT2dtdv_4PNSSLCoeff(m2M, qm_def2)*S2L*S2L;

    return;
}

static int test_consistency(
    const REAL8 m1M,
    const REAL8 chi1,
    const REAL8 chi2,
    const REAL8 qm_def1,
    const REAL8 qm_def2
    )
{
    REAL8 m2M = 1.-m1M;
    REAL8 eta = m1M*m2M;

    int ret = 0;

    fprintf(stdout, "\n=== Testing eta=%.4f, chi1=%.4f, chi2=%.4f ===\n", eta, chi1, chi2);

    PNPhasingSeries dtdv_ef;
    dtdv_from_energy_flux(&dtdv_ef, m1M, chi1, chi2, qm_def1, qm_def2);

    PNPhasingSeries dtdv_pn;
    dtdv_from_pncoefficients(&dtdv_pn, m1M, chi1, chi2, qm_def1, qm_def2);

    PNPhasingSeries phasing;
    XLALSimInspiralPNPhasing_F2(&phasing, m1M,m2M, chi1, chi2,\
                                chi1*chi1, chi2*chi2, chi1*chi2,\
                                qm_def1, qm_def2, 7);

    /* Divide the phasing by the leading-order term */
    REAL8 phase0 = phasing.v[0];
    for (int i = 1; i < PN_PHASING_SERIES_MAX_ORDER; i++)
    {
        phasing.v[i] /= phase0;
        phasing.vlogv[i] /= phase0;
        phasing.vlogvsq[i] /= phase0;
    }

    fprintf(stdout, "Testing dtdv consistency with energy and flux.\n");
    ret += compare_dtdv(&dtdv_pn, &dtdv_ef);

    fprintf(stdout, "Testing phasing consistency with dtdv.\n");
    ret += compare(phasing.v[0], 3./20.*dtdv_pn.v[0], 0, 0);
    ret += compare_phasing_to_dtdv(&phasing, &dtdv_pn);

    fprintf(stdout, "Testing phasing consistency with energy and flux.\n");
    ret += compare(phasing.v[0], 3./20.*dtdv_ef.v[0], 0, 0);
    ret += compare_phasing_to_dtdv(&phasing, &dtdv_ef);

    return ret;
}

/* Testing tidal coefficients. Since they are symmetric with respect to both objects
 * it is sufficient to test only one non-zero coefficient (lambda2 below).  */

static int test_tidal(
    const REAL8 m2M,
    const REAL8 lambda2
    )
{
    REAL8 m1M = 1.L-m2M;
    REAL8 eta = m1M*m2M;

    REAL8 chi1 = m1M;
    REAL8 chi2 = m2M;

    REAL8 energy2 = XLALSimInspiralPNEnergy_2PNCoeff(eta);
    REAL8 flux2 = XLALSimInspiralPNFlux_2PNCoeff(eta);
    REAL8 energy10 = XLALSimInspiralPNEnergy_10PNTidalCoeff(chi1, chi2, lambda2);
    REAL8 flux10 = XLALSimInspiralPNFlux_10PNTidalCoeff(chi2, lambda2);
    REAL8 energy12 = XLALSimInspiralPNEnergy_12PNTidalCoeff(chi1, chi2, lambda2);
    REAL8 flux12 = XLALSimInspiralPNFlux_12PNTidalCoeff(chi2, lambda2);

    REAL8 dtdv2 = 2.L*energy2 - flux2;
    REAL8 dtdv10 = 6.L*energy10 - flux10;
    REAL8 dtdv12 = (7.L*energy12 - flux12) - flux2*dtdv10 - flux10*dtdv2;

    REAL8 calc_phasing10 = 4.L*dtdv10;
    REAL8 calc_phasing12 = (10.L/7.L)*dtdv12;

    REAL8 phasing10 = XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(chi2, lambda2);
    REAL8 phasing12 = XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(chi2, lambda2);

    int ret = 0;
    ret += compare(calc_phasing10, phasing10, 10, 0);
    ret += compare(calc_phasing12, phasing12, 12, 0);

    return ret;
}


int main (int argc, char **argv)
{
    /* Ignore unused parameters. */
    (void)argc;
    (void)argv;

    int ret = 0;

    ret += test_consistency(0.5, 0., 0., 1., 4.5);
    ret += test_consistency(0.9, 0., 0., 1., 4.5);
    ret += test_consistency(0.01, 0., 0., 1., 4.5);

    ret += test_consistency(0.5, 1., 0., 1., 4.5);
    ret += test_consistency(0.9, 1., 0., 1., 4.5);
    ret += test_consistency(0.01, 1., 0., 1., 4.5);

    ret += test_consistency(0.5, 0., -0.4, 1., 4.5);
    ret += test_consistency(0.9, 0., -0.4, 1., 4.5);
    ret += test_consistency(0.01, 0., -0.4, 1., 4.5);

    ret += test_consistency(0.5, 0.9, -0.9, 2., 2.);
    ret += test_consistency(0.9, 0.9, -0.9, 3., 3.);
    ret += test_consistency(0.01, 0.9, -0.9, 4., 4.);

    fprintf(stdout, "Testing tidal terms.\n");
    ret += test_tidal(0.5, 0.);
    ret += test_tidal(0.5, 10.);
    ret += test_tidal(0.5, 9241.);
    ret += test_tidal(0.1, 0.);
    ret += test_tidal(0.1, 10.);
    ret += test_tidal(0.86, 9241.);
    ret += test_tidal(0.86, 0.);
    ret += test_tidal(0.86, 10.);
    ret += test_tidal(0.1, 9241.);

    if (ret == 0)
    {
        fprintf(stdout, "\nAll PNCoefficients tests passed.\n");
    }
    else
    {
        fprintf(stderr, "\nFAILURE: %u PNCoefficients comparisons incorrect.\n", ret);
    }

    return ret;
}
