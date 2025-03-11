/*
 * * Copyright (C) 2017 Tim Dietrich, Sebastiano Bernuzzi, Nathan Johnson-McDaniel,
 * * Shasvath J Kapadia, Francesco Pannarale and Sebastian Khan, Michael Puerrer.
 * *
 * *  This program is free software; you can redistribute it and/or modify
 * *  it under the terms of the GNU General Public License as published by
 * *  the Free Software Foundation; either version 2 of the License, or
 * *  (at your option) any later version.
 * *
 * *  This program is distributed in the hope that it will be useful,
 * *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 * *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * *  GNU General Public License for more details.
 * *
 * *  You should have received a copy of the GNU General Public License
 * *  along with with program; see the file COPYING. If not, write to the
 * *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * *  MA  02110-1301  USA
 * */

static int EnforcePrimaryMassIsm1(REAL8 *m1, REAL8 *m2, REAL8 *lambda1, REAL8 *lambda2);

static int EnforcePrimaryMassIsm1_v3(REAL8 *m1, REAL8 *m2, REAL8 *lambda1, REAL8 *lambda2, REAL8 *chi1, REAL8 *chi2);


static double SimNRTunedTidesFDTidalPhase(
    const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
    const REAL8 Xa, /**< Mass of companion 1 divided by total mass */
    const REAL8 Xb, /**< Mass of companion 2 divided by total mass */
    const REAL8 mtot, /**< total mass (Msun) */
    const REAL8 kappa2T /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
);

static REAL8 PlanckTaper(const REAL8 t, const REAL8 t1, const REAL8 t2);

static REAL8 SimNRTunedTidesFDTidalAmplitude(
    const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
    const REAL8 mtot, /**< Total mass in solar masses */
    const REAL8 kappa2T /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
    );

static double SimNRTunedTidesFDTidalPhase_v2(
    const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
    const REAL8 Xa, /**< Mass of companion 1 divided by total mass */
    const REAL8 Xb, /**< Mass of companion 2 divided by total mass */
    const REAL8 mtot, /**< total mass (Msun) */
    const REAL8 kappa2T /**< tidal coupling constant. Eq. 2 in arXiv:1706.02969 */
    );

static double SimNRTunedTidesFDTidalPhase_v3(
    const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
    const REAL8 mtot, /**< total mass (Msun) */
    const REAL8 NRTidalv3_coeffs[20], /**< NRTidalv3 coefficients*/
    const REAL8 PN_coeffs[10] /**< 7.5 PN coefficients to be used as constraints*/
    );

static double SimNRTunedTidesFDTidalPhase_PN(
    const REAL8 fHz, /**< Gravitational wave frequency (Hz) */
    const REAL8 Xa, /**< Mass of companion 1 divided by total mass */
    const REAL8 mtot, /**< total mass (Msun) */
    const REAL8 lambda1, /**< dimensionless tidal deformability of companion 1*/
    const REAL8 lambda2, /**< dimensionless tidal deformability of companion 2*/
    const REAL8 PN_coeffs[10] /**< 7.5 PN coefficients*/
    );
