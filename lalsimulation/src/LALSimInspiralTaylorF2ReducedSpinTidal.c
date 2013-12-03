/*
*  Copyright (C) 2011 P. Ajith, N. Fotopoulos, J. Read
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

#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>

/**
 * Generate the "reduced-spin templates" proposed in http://arxiv.org/abs/1107.1267
 * Add the tidal phase terms from http://arxiv.org/abs/1101.1673 (Eqs. 3.9, 3.10)
 * The chi parameter should be determined from XLALSimInspiralTaylorF2ReducedSpinComputeChi.
 */
int XLALSimInspiralTaylorF2ReducedSpinTidal(
    COMPLEX16FrequencySeries **htilde,   /**< FD waveform */
    const REAL8 phic,                /**< orbital coalescence phase (rad) */
    const REAL8 deltaF,              /**< frequency resolution (Hz) */
    const REAL8 m1_SI,               /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,               /**< mass of companion 2 (kg) */
    const REAL8 chi,                 /**< dimensionless aligned-spin param */
    const REAL8 lam1,                /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
    const REAL8 lam2,                /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
    const REAL8 fStart,              /**< start GW frequency (Hz) */
    const REAL8 fEnd,                /**< highest GW frequency (Hz) of waveform generation - if 0, end at Schwarzschild ISCO */
    const REAL8 r,                   /**< distance of source (m) */
    const INT4 phaseO,               /**< twice PN phase order */
    const INT4 ampO                  /**< twice PN amplitude order */
    ) {
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 xi1 = m1 / m;
    const REAL8 xi2 = m2 / m;
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 mSevenBySix = -7./6.;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    REAL8 v0 = cbrt(piM * fStart);
    REAL8 shft, amp0, f_max;
    REAL8 psiNewt, psi2, psi3, psi4, psi5, psi6, psi6L, psi7, psi3S, psi4S, psi5S, psi10T1, psi10T2, psi10, psi12T1, psi12T2, psi12;
    REAL8 alpha2, alpha3, alpha4, alpha5, alpha6, alpha6L, alpha7, alpha3S, alpha4S, alpha5S;
    size_t i, n, iStart;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tStart = {0, 0};

    /* check inputs for sanity */
    if (*htilde) XLAL_ERROR(XLAL_EFAULT);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fabs(chi) > 1) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);
    if (ampO > 7) XLAL_ERROR(XLAL_EDOM); /* only implemented to pN 3.5 */
    if (phaseO > 7) XLAL_ERROR(XLAL_EDOM); /* only implemented to pN 3.5 */

    /* allocate htilde */
    if ( fEnd == 0. ) // End at ISCO
        f_max = fISCO;
    else // End at user-specified freq.
        f_max = fEnd;
    n = (size_t) (f_max / deltaF + 1);
    XLALGPSAdd(&tStart, -1 / deltaF);  /* coalesce at t=0 */
    *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tStart, 0.0, deltaF, &lalStrainUnit, n);
    if (!(*htilde)) XLAL_ERROR(XLAL_EFUNC);
    memset((*htilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);

    /* extrinsic parameters */
    amp0 = -pow(m_sec, 5./6.) * sqrt(5. * eta / 24.) / (cbrt(LAL_PI * LAL_PI) * r / LAL_C_SI);
    shft = LAL_TWOPI * (tStart.gpsSeconds + 1e-9 * tStart.gpsNanoSeconds);

    /* spin terms in the amplitude and phase (in terms of the reduced
     * spin parameter */
    psi3S = 113.*chi/3.;
    psi4S = 63845.*(-81. + 4.*eta)*chi*chi/(8.*pow(-113. + 76.*eta, 2.));
    psi5S = -565.*(-146597. + 135856.*eta + 17136.*eta*eta)*chi/(2268.*(-113. + 76.*eta));

    alpha3S = (113.*chi)/24.;
    alpha4S = (12769.*chi*chi*(-81. + 4.*eta))/(32.*pow(-113. + 76.*eta,2));
    alpha5S = (-113.*chi*(502429. - 591368.*eta + 1680*eta*eta))/(16128.*(-113 + 76*eta));

    /* tidal terms in the phase */
    psi10T2 = -24./xi2 * (1. + 11. * xi1) * lam2 * xi2*xi2*xi2*xi2*xi2;
    psi10T1 = -24./xi1 * (1. + 11. * xi2) * lam1 * xi1*xi1*xi1*xi1*xi1;
    psi12T2 = -5./28./xi2 * (3179. - 919.* xi2 - 2286.* xi2*xi2 + 260.* xi2*xi2*xi2)* lam2 * xi2*xi2*xi2*xi2*xi2;
    psi12T1 = -5./28./xi1 * (3179. - 919.* xi1 - 2286.* xi1*xi1 + 260.* xi1*xi1*xi1)* lam1 * xi1*xi1*xi1*xi1*xi1;

    /* coefficients of the phase at PN orders from 0 to 3.5PN */
    psiNewt = 3./(128.*eta);
    psi2 = 3715./756. + 55.*eta/9.;
    psi3 = psi3S - 16.*LAL_PI;
    psi4 = 15293365./508032. + 27145.*eta/504. + 3085.*eta*eta/72. + psi4S;
    psi5 = (38645.*LAL_PI/756. - 65.*LAL_PI*eta/9. + psi5S);
    psi6 = 11583231236531./4694215680. - (640.*LAL_PI*LAL_PI)/3. - (6848.*LAL_GAMMA)/21.
             + (-5162.983708047263 + 2255.*LAL_PI*LAL_PI/12.)*eta
             + (76055.*eta*eta)/1728. - (127825.*eta*eta*eta)/1296.;
    psi6L = -6848./21.;
    psi7 = (77096675.*LAL_PI)/254016. + (378515.*LAL_PI*eta)/1512.
             - (74045.*LAL_PI*eta*eta)/756.;
    psi10 = psi10T1+psi10T2;
    psi12 = psi12T1+psi12T2;

    /* amplitude coefficients */
    alpha2 = 1.1056547619047619 + (11*eta)/8.;
    alpha3 = -LAL_TWOPI + alpha3S;
    alpha4 = 0.8939214212884228 + (18913*eta)/16128. + (1379*eta*eta)/1152. + alpha4S;
    alpha5 = (-4757*LAL_PI)/1344. + (57*eta*LAL_PI)/16. + alpha5S;
    alpha6 = -58.601030974347324 + (3526813753*eta)/2.7869184e7 -
                (1041557*eta*eta)/258048. + (67999*eta*eta*eta)/82944. +
                (10*LAL_PI*LAL_PI)/3. - (451*eta*LAL_PI*LAL_PI)/96.;
    alpha6L = 856/105.;
    alpha7 = (-5111593*LAL_PI)/2.709504e6 - (72221*eta*LAL_PI)/24192. -
                (1349*eta*eta*LAL_PI)/24192.;

    /* select the terms according to the PN order chosen */
    switch (ampO) {
        case 0:
        case 1:
            alpha2 = 0.;
        case 2:
            alpha3 = 0.;
        case 3:
            alpha4 = 0.;
        case 4:
            alpha5 = 0.;
        case 5:
            alpha6 = 0.;
            alpha6L = 0.;
        case 6:
            alpha7 = 0.;
        default:
            break;
    }

    switch (phaseO) {
        case 0:
        case 1:
            psi2 = 0.;
        case 2:
            psi3 = 0.;
        case 3:
            psi4 = 0.;
        case 4:
            psi5 = 0.;
        case 5:
            psi6 = 0.;
            psi6L = 0.;
        case 6:
            psi7 = 0.;
        default:
            break;
    }

    /* Fill with non-zero vals from fStart to lesser of fEnd, fISCO */
    iStart = (size_t) ceil(fStart / deltaF);
    data = (*htilde)->data->data;
    const REAL8 logv0=log(v0);
    const REAL8 log4=log(4.0);
    
    for (i = iStart; i < n; i++) {
        /* fourier frequency corresponding to this bin */
        const REAL8 f = i * deltaF;
        const REAL8 v3 = piM*f;

        /* PN expansion parameter */
        REAL8 Psi, amp;
        const REAL8 v = cbrt(v3);
	const REAL8 logv=log(v);
        const REAL8 v2 = v*v;
	const REAL8 v4 = v3*v;
	const REAL8 v5 = v4*v;
	const REAL8 v6 = v3*v3;
	const REAL8 v7 = v6*v;
        const REAL8 v10 = v5*v5;
	const REAL8 v12 = v6*v6;

        /* compute the phase and amplitude */
        Psi = psiNewt / v5 * (1.
            + psi2 * v2 + psi3 * v3 + psi4 * v4
            + psi5 * v5 * (1. + 3. * (logv - logv0))
            + (psi6 + psi6L * (log4 + logv)) * v6 + psi7 * v7)
            + psi10 * v10 + psi12 * v12;

        amp = amp0 * pow(f, mSevenBySix) * (1.
            + alpha2 * v2 + alpha3 * v3 + alpha4 * v4 + alpha5 * v5
            + (alpha6 + alpha6L * (LAL_GAMMA + (log4 + logv))) * v6
            + alpha7 * v7);

        data[i] = amp * cos(Psi + shft * f - 2.*phic - LAL_PI_4) - I * (amp * sin(Psi + shft * f - 2.*phic - LAL_PI_4));
    }

    return XLAL_SUCCESS;
}
