/*
 *  Copyright (C) 2014 Andrew Lundgren, 2017 Riccardo Sturani
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

#include <stdlib.h>
#include <math.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include <lal/LALSimInspiralTestGRParams.h>
#include "LALSimInspiralPNCoefficients.c"

#define EPSILON 1.e-11

static int compare_value(
    REAL8 val1,
    REAL8 val2)
{
    if (fabs(val1 - val2) > EPSILON)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

static int compare(
    REAL8 val1,
    REAL8 val2,
    int v_order,
    int log_order)
{
    if (fabs(val1 - val2) > EPSILON)
    {
        if (log_order == 0) { fprintf(stderr, "FAILED at %.1f PN order: %.11f versus %.11f\n", 0.5*v_order, val1, val2); }
        else { fprintf(stderr, "FAILED at %.1f PN order, in log^%u term: %.11f versus %.11f\n", 0.5*v_order, log_order, val1, val2); }
        return 1;
    }
    else
    {
        return 0;
    }
}

static int compare_pnseries(
    PNPhasingSeries *s1,
    PNPhasingSeries *s2)
{
    int ret = 0;

    ret += compare(s1->v[0], s2->v[0], 0, 0);
    ret += compare(s1->v[2], s2->v[2], 2, 0);
    ret += compare(s1->v[3], s2->v[3], 3, 0);
    ret += compare(s1->v[4], s2->v[4], 4, 0);
    ret += compare(s1->v[5], s2->v[5], 5, 0);
    ret += compare(s1->v[6], s2->v[6], 6, 0);
    ret += compare(s1->vlogv[6], s2->vlogv[6], 6, 1);
    ret += compare(s1->v[7], s2->v[7], 7, 0);

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

/* Helper function to calculate sum of (1+j/2)energy[j]*wdot[k-j] from 1 to k-1 */
static REAL8 sumE(
    REAL8 *ene,
    REAL8 *wdot,
    int k)
{
    REAL8 accum = 0.;
    for (int j=1; j < k; j++)
    {
      accum += (1.+((double) j)/2.)*ene[j]*wdot[k-j];
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
    /* Check is performed for aligned spin only*/
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

    REAL8 flux6l = XLALSimInspiralPNFlux_6PNLogCoeff(eta);

    memset(dtdv, 0, sizeof(PNPhasingSeries));
    dtdv->v[0] = -2.*XLALSimInspiralPNEnergy_0PNCoeff(eta) / XLALSimInspiralPNFlux_0PNCoeff(eta);
    dtdv->v[1] = 0.; /* there's no 0.5 PN term */
    for (int i = 2; i < 8; i++)
    {
      dtdv->v[i] = (1.+0.5*((double)i))*energy[i] - flux[i] - sum(flux, dtdv->v, i);
    }
    dtdv->vlogv[6] = -flux6l;

    /* Calculate the leading-order spin-spin terms separately
     * FIXME: For now, only do aligned spins
     */
    REAL8 energy_ss4 = XLALSimInspiralPNEnergy_4PNS1S2Coeff(eta)*S1L*S2L;
    energy_ss4 += qm_def1*XLALSimInspiralPNEnergy_4PNQMS1S1Coeff(m1M)*S1L*S1L;
    energy_ss4 += qm_def2*XLALSimInspiralPNEnergy_4PNQMS1S1Coeff(m2M)*S2L*S2L;

    REAL8 flux_ss4 = XLALSimInspiralPNFlux_4PNS1S2Coeff(eta)*S1L*S2L;
    flux_ss4 += XLALSimInspiralPNFlux_4PNS1S1Coeff(m1M)*S1L*S1L;
    flux_ss4 += XLALSimInspiralPNFlux_4PNS1S1Coeff(m2M)*S2L*S2L;
    flux_ss4 += qm_def1*XLALSimInspiralPNFlux_4PNQMS1S1Coeff(m1M)*S1L*S1L;
    flux_ss4 += qm_def2*XLALSimInspiralPNFlux_4PNQMS1S1Coeff(m2M)*S2L*S2L;

    dtdv->v[4] += 3.*energy_ss4 - flux_ss4;

    REAL8 energy_ss6 = XLALSimInspiralPNEnergy_6PNS1S2Coeff(eta)*S1L*S2L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1OS2OCoeff(eta)*S1L*S2L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1S1Coeff(m1M)*S1L*S1L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1OS1OCoeff(m1M)*S1L*S1L;
    energy_ss6 += qm_def1*XLALSimInspiralPNEnergy_6PNQMS1S1Coeff(m1M)*S1L*S1L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1S1Coeff(m2M)*S2L*S2L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1OS1OCoeff(m2M)*S2L*S2L;
    energy_ss6 += qm_def2*XLALSimInspiralPNEnergy_6PNQMS1S1Coeff(m2M)*S2L*S2L;

    REAL8 flux_ss6 = XLALSimInspiralPNFlux_6PNS1S2Coeff(eta)*S1L*S2L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1OS2OCoeff(eta)*S1L*S2L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1S1Coeff(m1M)*S1L*S1L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1OS1OCoeff(m1M)*S1L*S1L;
    flux_ss6 += qm_def1*XLALSimInspiralPNFlux_6PNQMS1S1Coeff(m1M)*S1L*S1L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1S1Coeff(m2M)*S2L*S2L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1OS1OCoeff(m2M)*S2L*S2L;
    flux_ss6 += qm_def2*XLALSimInspiralPNFlux_6PNQMS1S1Coeff(m2M)*S2L*S2L;

    dtdv->v[6] += 4.*energy_ss6 - flux_ss6 -3.*flux[2]*energy_ss4 - 2.*flux_ss4*energy[2] + 2.*flux[2]*flux_ss4;

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

    //Test for spin aligned only
    dtdv->v[4] += XLALSimInspiralTaylorT2dtdv_4PNS1S2Coeff(eta)*S1L*S2L;
    dtdv->v[4] += (XLALSimInspiralTaylorT2dtdv_4PNS1S1Coeff(m1M)+qm_def1*XLALSimInspiralTaylorT2dtdv_4PNQMS1S1Coeff(m1M))*S1L*S1L;
    dtdv->v[4] += (XLALSimInspiralTaylorT2dtdv_4PNS1S1Coeff(m2M)+qm_def2*XLALSimInspiralTaylorT2dtdv_4PNQMS1S1Coeff(m2M))*S2L*S2L;

    dtdv->v[6] += XLALSimInspiralTaylorT2dtdv_6PNS1S2Coeff(eta)*S1L*S2L;
    dtdv->v[6] += XLALSimInspiralTaylorT2dtdv_6PNS1OS2OCoeff(eta)*S1L*S2L;
    dtdv->v[6] += (XLALSimInspiralTaylorT2dtdv_6PNS1S1Coeff(m1M)+XLALSimInspiralTaylorT2dtdv_6PNS1OS1OCoeff(m1M)+qm_def1*XLALSimInspiralTaylorT2dtdv_6PNQMS1S1Coeff(m1M))*S1L*S1L;
    dtdv->v[6] += (XLALSimInspiralTaylorT2dtdv_6PNS1S1Coeff(m2M)+XLALSimInspiralTaylorT2dtdv_6PNS1OS1OCoeff(m2M)+qm_def2*XLALSimInspiralTaylorT2dtdv_6PNQMS1S1Coeff(m2M))*S2L*S2L;

    return;
}

static int test_average(const REAL8 m1M)
{
  const REAL8 m2M=1.-m1M;
  const REAL8 eta=m1M*m2M;

  int ret=0;

  printf("Testing consistency between averaged and instantaneous coefficients (2 and 3PN spin^2)\n");

  ret+= compare_value(XLALSimInspiralPNEnergy_4PNS1S2Coeff(eta),XLALSimInspiralPNEnergy_4PNS1S2CoeffAvg(eta)+XLALSimInspiralPNEnergy_4PNS1OS2OCoeffAvg(eta));
  ret+= compare_value(-0.5*XLALSimInspiralPNEnergy_4PNS1nS2nCoeff(eta),XLALSimInspiralPNEnergy_4PNS1OS2OCoeffAvg(eta));

  ret+=compare_value(XLALSimInspiralPNEnergy_4PNQMS1S1Coeff(m1M),XLALSimInspiralPNEnergy_4PNQMS1S1CoeffAvg(m1M)+XLALSimInspiralPNEnergy_4PNQMS1OS1OCoeffAvg(m1M));
  ret+=compare_value(-0.5*XLALSimInspiralPNEnergy_4PNQMS1nS1nCoeff(m1M),XLALSimInspiralPNEnergy_4PNQMS1OS1OCoeffAvg(m1M));

  ret+=compare_value(XLALSimInspiralPNFlux_4PNS1S2Coeff(eta),XLALSimInspiralPNFlux_4PNS1S2CoeffAvg(eta)+XLALSimInspiralPNFlux_4PNS1OS2OCoeffAvg(eta));
  ret+=compare_value(-0.5*(XLALSimInspiralPNFlux_4PNS1nS2nCoeff(eta)+XLALSimInspiralPNFlux_4PNS1vS2vCoeff(eta)),XLALSimInspiralPNFlux_4PNS1OS2OCoeffAvg(eta));

  ret+=compare_value(XLALSimInspiralPNFlux_4PNS1S1Coeff(m1M),XLALSimInspiralPNFlux_4PNS1S1CoeffAvg(m1M)+XLALSimInspiralPNFlux_4PNS1OS1OCoeffAvg(m1M));
  ret+=compare_value(-0.5*XLALSimInspiralPNFlux_4PNS1vS1vCoeff(m1M),XLALSimInspiralPNFlux_4PNS1OS1OCoeffAvg(m1M));

  ret+=compare_value(XLALSimInspiralPNFlux_4PNQMS1S1Coeff(m1M),XLALSimInspiralPNFlux_4PNQMS1S1CoeffAvg(m1M)+XLALSimInspiralPNFlux_4PNQMS1OS1OCoeffAvg(m1M));
  ret+=compare_value(-0.5*(XLALSimInspiralPNFlux_4PNQMS1nS1nCoeff(m1M)+XLALSimInspiralPNFlux_4PNQMS1vS1vCoeff(m1M)),XLALSimInspiralPNFlux_4PNQMS1OS1OCoeffAvg(m1M));

  ret+=compare_value(XLALSimInspiralTaylorT4wdot_4PNS1S2Coeff(eta),XLALSimInspiralTaylorT4wdot_4PNS1S2CoeffAvg(eta)+XLALSimInspiralTaylorT4wdot_4PNS1OS2OCoeffAvg(eta));

  ret+=compare_value(XLALSimInspiralTaylorT4wdot_4PNS1S1Coeff(m1M),XLALSimInspiralTaylorT4wdot_4PNS1S1CoeffAvg(m1M)+XLALSimInspiralTaylorT4wdot_4PNS1OS1OCoeffAvg(m1M));

  ret+=compare_value(XLALSimInspiralTaylorT4wdot_4PNQMS1S1Coeff(m1M),XLALSimInspiralTaylorT4wdot_4PNQMS1S1CoeffAvg(m1M)+XLALSimInspiralTaylorT4wdot_4PNQMS1OS1OCoeffAvg(m1M));

  ret+=compare_value(XLALSimInspiralTaylorT2dtdv_4PNS1S2Coeff(eta),XLALSimInspiralTaylorT2dtdv_4PNS1S2CoeffAvg(eta)+XLALSimInspiralTaylorT2dtdv_4PNS1OS2OCoeffAvg(eta));

  ret+=compare_value(XLALSimInspiralTaylorT2dtdv_4PNS1S1Coeff(m1M),XLALSimInspiralTaylorT2dtdv_4PNS1S1CoeffAvg(m1M)+XLALSimInspiralTaylorT2dtdv_4PNS1OS1OCoeffAvg(m1M));

  ret+=compare_value(XLALSimInspiralTaylorT2dtdv_4PNQMS1S1Coeff(m1M),XLALSimInspiralTaylorT2dtdv_4PNQMS1S1CoeffAvg(m1M)+XLALSimInspiralTaylorT2dtdv_4PNQMS1OS1OCoeffAvg(m1M));

  return ret;

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

    fprintf(stdout, "\n=== Testing eta=%.4f, chi1=%.4f, chi2=%.4f, qm1=%.4f, qm2=%.4f ===\n", eta, chi1, chi2, qm_def1, qm_def2);

    PNPhasingSeries dtdv_ef;
    dtdv_from_energy_flux(&dtdv_ef, m1M, chi1, chi2, qm_def1, qm_def2);

    PNPhasingSeries dtdv_pn;
    dtdv_from_pncoefficients(&dtdv_pn, m1M, chi1, chi2, qm_def1, qm_def2);

    LALDict *extraParams=XLALCreateDict();
    XLALSimInspiralWaveformParamsInsertdQuadMon1(extraParams,qm_def1-1.);
    XLALSimInspiralWaveformParamsInsertdQuadMon2(extraParams,qm_def2-1.);
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,7);
    PNPhasingSeries phasing;
    XLALSimInspiralPNPhasing_F2(&phasing, m1M,m2M, chi1, chi2,\
                                chi1*chi1, chi2*chi2, chi1*chi2,\
                                extraParams);
    XLALDestroyDict(extraParams);

    /* Divide the phasing by the leading-order term */
    REAL8 phase0 = phasing.v[0];
    for (int i = 1; i < PN_PHASING_SERIES_MAX_ORDER; i++)
    {
        phasing.v[i] /= phase0;
        phasing.vlogv[i] /= phase0;
        phasing.vlogvsq[i] /= phase0;
    }

    fprintf(stdout, "Testing dtdv consistency with energy and flux.\n");
    ret += compare_pnseries(&dtdv_pn, &dtdv_ef);

    fprintf(stdout, "Testing phasing consistency with dtdv.\n");
    ret += compare(phasing.v[0], 3./20.*dtdv_pn.v[0], 0, 0);
    ret += compare_phasing_to_dtdv(&phasing, &dtdv_pn);

    fprintf(stdout, "Testing phasing consistency with energy and flux.\n");
    ret += compare(phasing.v[0], 3./20.*dtdv_ef.v[0], 0, 0);
    ret += compare_phasing_to_dtdv(&phasing, &dtdv_ef);

    return ret;
}

/* The TaylorT4 waveform is defined by
 * dw/dt = trunc[ -flux(v)/dE/dw(v) ]
 * where trunc[ ] means truncate the expression at the desired PN order.
 * There's a recursive expression for dw/dt at each order. Write
 * dw/dt   = W_0 ( 1 + \sum_{i=1}^n W_i v^i )
 * e(v)    = E_0 ( 1 + \sum_{i=1}^n E_i v^i )
 * flux(v) = F_0 ( 1 + \sum_{i=1}^n F_i v^i )
 *
 * Then
 *
 * W_i = F_i - (1+i/2)E_i - \sum_{j=1}^{i-1} (1+i/2) E_j T_{i-j}
 *
 * This doesn't handle terms involving log(v), so we'll treat them separately.
 * Be careful with the dE/dv, which causes the log and non-log terms to mix.
 * We first add in the spin-orbit terms, but
 * only keeping terms linear in spin. Then we add in the spin-spin terms.
 */

static void T4wdot_from_energy_flux(
    PNPhasingSeries *wdot,
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
    energy[3] = 0.;
    energy[4] = XLALSimInspiralPNEnergy_4PNCoeff(eta);
    energy[5] = 0.;
    energy[6] = XLALSimInspiralPNEnergy_6PNCoeff(eta);
    energy[7] = 0.;
    //We do not add the non-spin part of the 4PN energy as the 4PN flux is unknown, apart from the SL terms
    energy[8] = 0.;

    flux[1] = 0.;
    flux[2] = XLALSimInspiralPNFlux_2PNCoeff(eta);
    flux[3] = XLALSimInspiralPNFlux_3PNCoeff(eta);
    flux[4] = XLALSimInspiralPNFlux_4PNCoeff(eta);
    flux[5] = XLALSimInspiralPNFlux_5PNCoeff(eta);
    flux[6] = XLALSimInspiralPNFlux_6PNCoeff(eta);
    flux[7] = XLALSimInspiralPNFlux_7PNCoeff(eta);
    flux[8] = 0.;

    /* Add the spin-orbit contributions */
    REAL8 energy_so3=XLALSimInspiralPNEnergy_3PNSOCoeff(m1M)*S1L + XLALSimInspiralPNEnergy_3PNSOCoeff(m2M)*S2L;
    energy[3] += energy_so3;
    REAL8 energy_so5=XLALSimInspiralPNEnergy_5PNSOCoeff(m1M)*S1L + XLALSimInspiralPNEnergy_5PNSOCoeff(m2M)*S2L;
    energy[5] += energy_so5;
    energy[7] +=XLALSimInspiralPNEnergy_7PNSOCoeff(m1M)*S1L + XLALSimInspiralPNEnergy_7PNSOCoeff(m2M)*S2L;

    flux[3] += XLALSimInspiralPNFlux_3PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_3PNSOCoeff(m2M)*S2L;
    flux[5] += XLALSimInspiralPNFlux_5PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_5PNSOCoeff(m2M)*S2L;
    REAL8 flux_so6=XLALSimInspiralPNFlux_6PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_6PNSOCoeff(m2M)*S2L;
    flux[6] += flux_so6;
    flux[7] += XLALSimInspiralPNFlux_7PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_7PNSOCoeff(m2M)*S2L;
    flux[8] += XLALSimInspiralPNFlux_8PNSOCoeff(m1M)*S1L + XLALSimInspiralPNFlux_8PNSOCoeff(m2M)*S2L;

    memset(wdot, 0, sizeof(PNPhasingSeries));
    wdot->v[0] = -XLALSimInspiralPNFlux_0PNCoeff(eta)/(2./3.*XLALSimInspiralPNEnergy_0PNCoeff(eta));
    wdot->v[1] = 0.; /* there's no 0.5 PN term */
    for (int i = 2; i <8; i++)
    {
        wdot->v[i] = flux[i] - (1.+0.5*i)*energy[i] - sumE(energy, wdot->v, i);
    }
    wdot->vlogv[6]= XLALSimInspiralPNFlux_6PNLogCoeff(eta);

    // The 8PN SO term is checked by hand
    wdot->v[8]=flux[8]-5.*energy[8] - flux_so6*2.*XLALSimInspiralPNEnergy_2PNCoeff(eta) - XLALSimInspiralPNFlux_5PNCoeff(eta)*5./2.*energy_so3 + XLALSimInspiralPNFlux_3PNCoeff(eta)*(-7./2.*energy_so5 + 10.*XLALSimInspiralPNEnergy_2PNCoeff(eta)*energy_so3);

    /* Calculate the leading-order spin-spin terms separately */
    REAL8 energy_ss4 = XLALSimInspiralPNEnergy_4PNS1S2Coeff(eta)*S1L*S2L;
    energy_ss4 += qm_def1*XLALSimInspiralPNEnergy_4PNQMS1S1Coeff(m1M)*S1L*S1L;
    energy_ss4 += qm_def2*XLALSimInspiralPNEnergy_4PNQMS1S1Coeff(m2M)*S2L*S2L;

    REAL8 flux_ss4 = XLALSimInspiralPNFlux_4PNS1S2Coeff(eta)*S1L*S2L;
    flux_ss4 += XLALSimInspiralPNFlux_4PNS1S1Coeff(m1M)*S1L*S1L;
    flux_ss4 += XLALSimInspiralPNFlux_4PNS1S1Coeff(m2M)*S2L*S2L;
    flux_ss4 += qm_def1*XLALSimInspiralPNFlux_4PNQMS1S1Coeff(m1M)*S1L*S1L;
    flux_ss4 += qm_def2*XLALSimInspiralPNFlux_4PNQMS1S1Coeff(m2M)*S2L*S2L;

    wdot->v[4] += flux_ss4 -3.*energy_ss4;

    REAL8 energy_ss6 = XLALSimInspiralPNEnergy_6PNS1S2Coeff(eta)*S1L*S2L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1OS2OCoeff(eta)*S1L*S2L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1S1Coeff(m1M)*S1L*S1L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1OS1OCoeff(m1M)*S1L*S1L;
    energy_ss6 += qm_def1*XLALSimInspiralPNEnergy_6PNQMS1S1Coeff(m1M)*S1L*S1L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1S1Coeff(m2M)*S2L*S2L;
    energy_ss6 += XLALSimInspiralPNEnergy_6PNS1OS1OCoeff(m2M)*S2L*S2L;
    energy_ss6 += qm_def2*XLALSimInspiralPNEnergy_6PNQMS1S1Coeff(m2M)*S2L*S2L;

    REAL8 flux_ss6 = XLALSimInspiralPNFlux_6PNS1S2Coeff(eta)*S1L*S2L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1OS2OCoeff(eta)*S1L*S2L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1S1Coeff(m1M)*S1L*S1L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1OS1OCoeff(m1M)*S1L*S1L;
    flux_ss6 += qm_def1*XLALSimInspiralPNFlux_6PNQMS1S1Coeff(m1M)*S1L*S1L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1S1Coeff(m2M)*S2L*S2L;
    flux_ss6 += XLALSimInspiralPNFlux_6PNS1OS1OCoeff(m2M)*S2L*S2L;
    flux_ss6 += qm_def2*XLALSimInspiralPNFlux_6PNQMS1S1Coeff(m2M)*S2L*S2L;

    wdot->v[6] += flux_ss6 - 4.*energy_ss6 -3.*flux[2]*energy_ss4 - 2.*flux_ss4*energy[2] + 12.*energy[2]*energy_ss4;

    return;
}

static void T4wdot_from_pncoefficients(
    PNPhasingSeries *wdot,
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

    memset(wdot, 0, sizeof(PNPhasingSeries));
    wdot->v[0] = XLALSimInspiralTaylorT4wdot_0PNCoeff(eta);
    wdot->v[1] = 0.;
    wdot->v[2] = XLALSimInspiralTaylorT4wdot_2PNCoeff(eta);
    wdot->v[3] = XLALSimInspiralTaylorT4wdot_3PNCoeff(eta);
    wdot->v[4] = XLALSimInspiralTaylorT4wdot_4PNCoeff(eta);
    wdot->v[5] = XLALSimInspiralTaylorT4wdot_5PNCoeff(eta);
    wdot->v[6] = XLALSimInspiralTaylorT4wdot_6PNCoeff(eta);
    wdot->vlogv[6] = XLALSimInspiralTaylorT4wdot_6PNLogCoeff(eta);
    wdot->v[7] = XLALSimInspiralTaylorT4wdot_7PNCoeff(eta);
    wdot->v[8] = 0.;

    wdot->v[3] += XLALSimInspiralTaylorT4wdot_3PNSOCoeff(m1M)*S1L + XLALSimInspiralTaylorT4wdot_3PNSOCoeff(m2M)*S2L;
    wdot->v[5] += XLALSimInspiralTaylorT4wdot_5PNSOCoeff(m1M)*S1L + XLALSimInspiralTaylorT4wdot_5PNSOCoeff(m2M)*S2L;
    wdot->v[6] += XLALSimInspiralTaylorT4wdot_6PNSOCoeff(m1M)*S1L + XLALSimInspiralTaylorT4wdot_6PNSOCoeff(m2M)*S2L;
    wdot->v[7] += XLALSimInspiralTaylorT4wdot_7PNSOCoeff(m1M)*S1L + XLALSimInspiralTaylorT4wdot_7PNSOCoeff(m2M)*S2L;
    wdot->v[8] += XLALSimInspiralTaylorT4wdot_8PNSOCoeff(m1M)*S1L + XLALSimInspiralTaylorT4wdot_8PNSOCoeff(m2M)*S2L;

    wdot->v[4] += XLALSimInspiralTaylorT4wdot_4PNS1S2Coeff(eta)*S1L*S2L;
    wdot->v[4] += XLALSimInspiralTaylorT4wdot_4PNS1S1Coeff(m1M)*S1L*S1L;
    wdot->v[4] += XLALSimInspiralTaylorT4wdot_4PNS1S1Coeff(m2M)*S2L*S2L;
    wdot->v[4] += qm_def1*XLALSimInspiralTaylorT4wdot_4PNQMS1S1Coeff(m1M)*S1L*S1L;
    wdot->v[4] += qm_def2*XLALSimInspiralTaylorT4wdot_4PNQMS1S1Coeff(m2M)*S2L*S2L;

    wdot->v[6] += XLALSimInspiralTaylorT4wdot_6PNS1S2Coeff(eta)*S1L*S2L;
    wdot->v[6] += XLALSimInspiralTaylorT4wdot_6PNS1OS2OCoeff(eta)*S1L*S2L;
    wdot->v[6] += XLALSimInspiralTaylorT4wdot_6PNS1S1Coeff(m1M)*S1L*S1L;
    wdot->v[6] += XLALSimInspiralTaylorT4wdot_6PNS1OS1OCoeff(m1M)*S1L*S1L;
    wdot->v[6] += XLALSimInspiralTaylorT4wdot_6PNS1S1Coeff(m2M)*S2L*S2L;
    wdot->v[6] += XLALSimInspiralTaylorT4wdot_6PNS1OS1OCoeff(m2M)*S2L*S2L;
    wdot->v[6] += qm_def1*XLALSimInspiralTaylorT4wdot_6PNQMS1S1Coeff(m1M)*S1L*S1L;
    wdot->v[6] += qm_def2*XLALSimInspiralTaylorT4wdot_6PNQMS1S1Coeff(m2M)*S2L*S2L;

    return;
}

static int test_consistency_T4(
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

    fprintf(stdout, "\n=== Testing T4  eta=%.4f, chi1=%.4f, chi2=%.4f, qm1=%.4f, qm2=%.4f ===\n", eta, chi1, chi2, qm_def1, qm_def2);

    PNPhasingSeries wdot_ef;
    T4wdot_from_energy_flux(&wdot_ef, m1M, chi1, chi2, qm_def1, qm_def2);

    PNPhasingSeries wdot_pn;
    T4wdot_from_pncoefficients(&wdot_pn, m1M, chi1, chi2, qm_def1, qm_def2);

    fprintf(stdout, "Testing wdot consistency with energy and flux.\n");
    ret += compare_pnseries(&wdot_pn, &wdot_ef);

    return ret;

}

/* Testing tidal coefficients. Since they are symmetric with respect to both objects
 * it is sufficient to test only one non-zero coefficient.  */

static int test_tidal_F2(
    const REAL8 m2M
    )
{
    REAL8 m1M = 1.L-m2M;
    REAL8 eta = m1M*m2M;

    REAL8 energy2 = XLALSimInspiralPNEnergy_2PNCoeff(eta);
    REAL8 flux2 = XLALSimInspiralPNFlux_2PNCoeff(eta);
    REAL8 energy10 = XLALSimInspiralPNEnergy_10PNTidalCoeff(m2M);
    REAL8 flux10 = XLALSimInspiralPNFlux_10PNTidalCoeff(m2M);
    REAL8 energy12 = XLALSimInspiralPNEnergy_12PNTidalCoeff(m2M);
    REAL8 flux12 = XLALSimInspiralPNFlux_12PNTidalCoeff(m2M);

    REAL8 dtdv2 = 2.L*energy2 - flux2;
    REAL8 dtdv10 = 6.L*energy10 - flux10;
    REAL8 dtdv12 = (7.L*energy12 - flux12) - flux2*dtdv10 - flux10*dtdv2;

    REAL8 calc_phasing10 = 4.L*dtdv10;
    REAL8 calc_phasing12 = (10.L/7.L)*dtdv12;

    REAL8 phasing10 = XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m2M);
    REAL8 phasing12 = XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m2M);

    int ret = 0;
    ret += compare(calc_phasing10, phasing10, 10, 0);
    ret += compare(calc_phasing12, phasing12, 12, 0);

    return ret;
}

static int test_tidal_T4(
    const REAL8 m2M
    )
{
    REAL8 m1M = 1.L-m2M;
    REAL8 eta = m1M*m2M;

    REAL8 energy2 = XLALSimInspiralPNEnergy_2PNCoeff(eta);
    REAL8 flux2 = XLALSimInspiralPNFlux_2PNCoeff(eta);
    REAL8 energy10 = XLALSimInspiralPNEnergy_10PNTidalCoeff(m2M);
    REAL8 flux10 = XLALSimInspiralPNFlux_10PNTidalCoeff(m2M);
    REAL8 energy12 = XLALSimInspiralPNEnergy_12PNTidalCoeff(m2M);
    REAL8 flux12 = XLALSimInspiralPNFlux_12PNTidalCoeff(m2M);

    REAL8 dvdt2 = flux2 - 2.L*energy2;
    REAL8 dvdt10 = flux10 - 6.L*energy10;
    REAL8 dvdt12 = (flux12 -7.L*energy12) - 2.L*energy2*dvdt10 - 6.L*energy10*dvdt2;

    REAL8 phasing10 = XLALSimInspiralTaylorT4wdot_10PNTidalCoeff(m2M);
    REAL8 phasing12 = XLALSimInspiralTaylorT4wdot_12PNTidalCoeff(m2M);

    int ret = 0;
    ret += compare(dvdt10, phasing10, 10, 0);
    ret += compare(dvdt12, phasing12, 12, 0);

    return ret;
}

static int test_tidal_T2(
    const REAL8 m2M
    )
{
    REAL8 m1M = 1.L-m2M;
    REAL8 eta = m1M*m2M;

    REAL8 energy2 = XLALSimInspiralPNEnergy_2PNCoeff(eta);
    REAL8 flux2 = XLALSimInspiralPNFlux_2PNCoeff(eta);
    REAL8 energy10 = XLALSimInspiralPNEnergy_10PNTidalCoeff(m2M);
    REAL8 flux10 = XLALSimInspiralPNFlux_10PNTidalCoeff(m2M);
    REAL8 energy12 = XLALSimInspiralPNEnergy_12PNTidalCoeff(m2M);
    REAL8 flux12 = XLALSimInspiralPNFlux_12PNTidalCoeff(m2M);

    REAL8 dtdv2 = 2.L*energy2 - flux2;
    REAL8 dtdv10 = 6.L*energy10 - flux10;
    REAL8 dtdv12 = (7.L*energy12 - flux12) - flux2*dtdv10 - flux10*dtdv2;

    REAL8 phasing10 = XLALSimInspiralTaylorT2dtdv_10PNTidalCoeff(m2M);
    REAL8 phasing12 = XLALSimInspiralTaylorT2dtdv_12PNTidalCoeff(m2M);

    int ret = 0;
    ret += compare(dtdv10, phasing10, 10, 0);
    ret += compare(dtdv12, phasing12, 12, 0);

    return ret;
}

/* The dL PN coefficients are defined in LALSimInspiralPNCoefficients.c.
 */

static void dL_from_pncoefficients(
    PNPhasingSeries *dL1,  /* Coefficients of \epsilon_{ijk}S1_jL_k */
    PNPhasingSeries *dL2,  /* Coefficients of \epsilon_{ijk}S2_jL_k */
    const REAL8 m1M
    )
{
    /* Check is performed for aligned spin only*/
    REAL8 m2M = 1.-m1M;
    /* Spins use the LAL convention */

    memset(dL1, 0, sizeof(PNPhasingSeries));
    memset(dL2, 0, sizeof(PNPhasingSeries));
    dL1->v[0] = 0.;   dL2->v[0] = 0.;
    dL1->v[1] = 0.;   dL2->v[1] = 0.;
    dL1->v[2] = 0.;   dL2->v[2] = 0.;
    dL1->v[3] = XLALSimInspiralLDot_3PNSOCoeff(m1M);
    dL2->v[3] = XLALSimInspiralLDot_3PNSOCoeff(m2M);
    dL1->v[4] = 0.;
    dL2->v[4] = 0.;
    dL1->v[5] = XLALSimInspiralLDot_5PNSOCoeff(m1M);
    dL2->v[5] = XLALSimInspiralLDot_5PNSOCoeff(m2M);
    dL1->v[6] = 0.;
    dL2->v[6] = 0.;
    dL1->v[7] = XLALSimInspiralLDot_7PNSOCoeff(m1M);
    dL2->v[7] = XLALSimInspiralLDot_7PNSOCoeff(m2M);

    return;
}

/* The dL PN coefficients are defined via the dS ones.
 */

static void dL_from_dSpncoefficients(
    PNPhasingSeries *dL1,  /* Coefficients of \epsilon_{ijk}S1_jL_k */
    PNPhasingSeries *dL2,  /* Coefficients of \epsilon_{ijk}S2_jL_k */
    const REAL8 m1M
    )
{
    /* Check is performed for aligned spin only*/
    REAL8 m2M = 1.-m1M;
    REAL8 eta = m1M*m2M;
    /* Spins use the LAL convention */

    memset(dL1, 0, sizeof(PNPhasingSeries));
    memset(dL2, 0, sizeof(PNPhasingSeries));
    dL1->v[0] = 0.;   dL2->v[0] = 0.;
    dL1->v[1] = 0.;   dL2->v[1] = 0.;
    dL1->v[2] = 0.;   dL2->v[2] = 0.;
    dL1->v[3] = XLALSimInspiralSpinDot_3PNCoeff(m1M)/eta;
    dL2->v[3] = XLALSimInspiralSpinDot_3PNCoeff(m2M)/eta;
    dL1->v[4] = 0.;
    dL2->v[4] = 0.;
    dL1->v[5] = XLALSimInspiralSpinDot_5PNCoeff(m1M)/eta;
    dL2->v[5] = XLALSimInspiralSpinDot_5PNCoeff(m2M)/eta;
    dL1->v[6] = 0.;
    dL2->v[6] = 0.;
    dL1->v[7] = XLALSimInspiralSpinDot_7PNCoeff(m1M)/eta;
    dL2->v[7] = XLALSimInspiralSpinDot_7PNCoeff(m2M)/eta;

    return;
}

static int test_consistency_dL(
    const REAL8 m1M)
{

    REAL8 m2M = 1.-m1M;
    REAL8 eta = m1M*m2M;

    int ret = 0;
    int idx;

    fprintf(stdout, "\n=== Testing dL  eta=%.4f ===\n", eta);

    PNPhasingSeries dL1_dL, dL2_dL;
    dL_from_pncoefficients(&dL1_dL, &dL2_dL, m1M);

    PNPhasingSeries dL1_dS, dL2_dS;
    dL_from_dSpncoefficients(&dL1_dS, &dL2_dS, m1M);

    fprintf(stdout, "Testing dL consistency with dS.\n");
    for (idx=0;idx<9;idx++) {
      ret += compare(dL1_dL.v[idx], dL1_dS.v[idx],idx,0);
      ret += compare(dL2_dL.v[idx], dL2_dS.v[idx],idx,0);
    }

    return ret;

}

int main (int argc, char **argv)
{
    /* Ignore unused parameters. */
    (void)argc;
    (void)argv;

    int ret = 0;

    ret += test_consistency(0.5, 0., 0., 0., 0.);
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

    ret += test_consistency_T4(0.5, 0., 0., 1., 4.5);
    ret += test_consistency_T4(0.9, 0., 0., 1., 4.5);
    ret += test_consistency_T4(0.01, 0., 0., 1., 4.5);

    ret += test_consistency_T4(0.5, 1., 0., 1., 4.5);
    ret += test_consistency_T4(0.9, 1., 0., 1., 4.5);
    ret += test_consistency_T4(0.01, 1., 0., 1., 4.5);

    ret += test_consistency_T4(0.5, 0., -0.4, 1., 4.5);
    ret += test_consistency_T4(0.9, 0., -0.4, 1., 4.5);
    ret += test_consistency_T4(0.01, 0., -0.4, 1., 4.5);

    ret += test_consistency_T4(0.5, 0.9, 0.9, 0., 0.);
    ret += test_consistency_T4(0.9, 0.9, -0.9, 3., 2.);
    ret += test_consistency_T4(0.01, 0.9, 0.9, 4., 4.);

    fprintf(stdout, "Testing tidal terms.\n");
    for (UINT4 idx=1;idx<=9;idx++) {
      ret += test_tidal_F2(0.1*((REAL8)idx));
      ret += test_tidal_T2(0.1*((REAL8)idx));
      ret += test_tidal_T4(0.1*((REAL8)idx));
      ret += test_average(0.1*((REAL8)idx));
    }

    for (UINT4 idx=1;idx<=5;idx++) {
      ret += test_consistency_dL(0.1*((REAL8)idx));
    }

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
