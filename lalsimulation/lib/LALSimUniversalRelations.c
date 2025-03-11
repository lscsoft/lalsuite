/*
 *  Copyright (C) 2017 Andrea Taracchini
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

#include <lal/LALSimUniversalRelations.h>

/**< Generic form of universal relation */
REAL8 XLALSimUniversalRelation( REAL8 x, REAL8 coeffs[] ) {
    return coeffs[0] + coeffs[1] * x + coeffs[2] * x * x + coeffs[3] * x * x * x + coeffs[4] * x *x * x * x;
}

/**< Eq. (60) with coeffs from 1st row of Table I of  https://arxiv.org/pdf/1311.0872.pdf */
/*    Gives the dimensionless l=3 tidal deformability: lambda3bar = 2/15 k3 C^7
 where k3 is the l=3 Love number and C is the compactness. It is a
 function the dimensionless l=2 tidal deformability: lambda2bar = 2/3 k2 C^5.
 Compared to NR for 1 <= lambda2bar <= 3000
 */
REAL8 XLALSimUniversalRelationlambda3TidalVSlambda2Tidal(
                                                                REAL8 lambda2bar /**< l=2 dimensionless tidal defomability */
)
{
    REAL8 coeffs[] = {-1.15, 1.18, 2.51e-2, -1.31e-3, 2.52e-5};
    REAL8 lnx;
    if ( lambda2bar < 0. ) {
        XLAL_ERROR (XLAL_EFUNC);
    }
    else if ( 0. <= lambda2bar && lambda2bar  < 0.01 ) {
        /* This is a function fitted to the universal relation in the range
         0.00001 <= lambda2hat <= 0.01 with the requirements that it goes to 0 at
         lambda2hat=0 and is exactly equal to the universal relation at lambda2hat=0.01 */
        return 0.4406491912035266*lambda2bar - 34.63232296075433*lambda2bar*lambda2bar
        + 1762.112913125107*lambda2bar*lambda2bar*lambda2bar;
    }
    else {
        lnx = log( lambda2bar );
    }
    REAL8 lny = XLALSimUniversalRelation( lnx, coeffs );
    return exp(lny);
}

/**< Eq. (3.5) with coeffs from 1st column of Table I of https://arxiv.org/pdf/1408.3789.pdf */
/*    Gives the l=2 f-mode frequency M_{NS}omega_{02}
 as a function the dimensionless l=2 tidal deformability: lambda2bar = 2/3 k2 C^5
 where k2 is the l=2 Love number and C is the compactness.
 Compared to NR for 0 <= log(lambda2bar) <= 9, that is
 1 <= lambda2bar <= 8100
 */
REAL8 XLALSimUniversalRelationomega02TidalVSlambda2Tidal(
                                                                REAL8 lambda2bar /**< l=2 dimensionless tidal defomability */
)
{
    REAL8 coeffs[] = {1.82e-1, -6.836e-3, -4.196e-3, 5.215e-4, -1.857e-5};
    REAL8 lnx;
    if ( lambda2bar < 0. ) {
        XLAL_ERROR (XLAL_EFUNC);
    }
    else if ( 0. <= lambda2bar && lambda2bar < 1. ) {
        lnx = 0.;
    }
    else if ( 1. <= lambda2bar && lambda2bar < exp(9.) ) {
        lnx = log( lambda2bar );
    }
    else {
        lnx = 9.;
    }
    return XLALSimUniversalRelation( lnx, coeffs );
}

/**< Eq. (3.5) with coeffs from 2nd column of Table I of https://arxiv.org/pdf/1408.3789.pdf */
/*    Gives the l=3 f-mode frequency M_{NS}omega_{03}
 as a function the dimensionless l=3 tidal deformability: lambda3bar = 2/15 k3 C^5
 where k3 is the l=3 Love number and C is the compactness.
 Compared to NR for -1 <= log(lambda3bar) <= 10, that is
 0.37 <= lambda3bar <= 20000
 */
REAL8 XLALSimUniversalRelationomega03TidalVSlambda3Tidal(
                                                                REAL8 lambda3bar /**< l=3 dimensionless tidal defomability */
)
{
    REAL8 coeffs[] = {2.245e-1, -1.5e-2, -1.412e-3, 1.832e-4, -5.561e-6};
    REAL8 lnx;
    if ( lambda3bar < 0. ) {
        XLAL_ERROR (XLAL_EFUNC);
    }
    else if ( 0. <= lambda3bar && lambda3bar < exp(-1.) ) {
        lnx = -1.;
    }
    else if ( exp(-1.) <= lambda3bar && lambda3bar < exp(10.) ) {
        lnx = log( lambda3bar );
    }
    else {
        lnx = 10.;
    }
    return XLALSimUniversalRelation( lnx, coeffs );
}

/**< Eq. (15) with coeffs from third row of Table I of https://arxiv.org/pdf/1608.02582.pdf (Yagi-Yunes) */
/* Gives the spin-induced quadrupole coefficient as a function the dimensionless l=2
 tidal deformability: lambda2bar = 2/3 k2 C^5.
 This coefficient quadparam relates the spin-induced quadrupole to the square of the spin, according to Q = -quadparam*chi^2*m^3*.
 It takes the value 1 for BH, and can reach ~10 for NS.
 The notation is Qbar in Yagi-Yunes. In the PN literature, the notation is often kappa (e.g. in
 https://arxiv.org/pdf/1501.01529.pdf). In https://arxiv.org/pdf/gr-qc/9709032.pdf the notation is a.
 The Yagi-Yunes fit does not cover the BH limit, where lambda2bar->0 and kappa->1.
 We extend it with a polynomial below lambda2bar=1. so that the function and its two
 first derivatives are smooth at the junction, while enforcing the BH limit.
 */
REAL8 XLALSimUniversalRelationQuadMonVSlambda2Tidal(
                                                    REAL8 lambda2bar /**< l=2 dimensionless tidal deformability */
)
{
    REAL8 coeffs[] = {0.1940, 0.09163, 0.04812, -4.283e-3, 1.245e-4};
    REAL8 lnx;
    if ( lambda2bar < 0. ) {
        XLAL_ERROR (XLAL_EFUNC);
    }
    else if ( 0. <= lambda2bar && lambda2bar  < 1. ) {
        /* Extension of the fit in the range lambda2bar=[0,1.] so that
        the BH limit is enforced, lambda2bar->0 gives quadparam->1. and
        the junction with the universal relation is smooth, of class C2  */
        return 1. + lambda2bar*(0.427688866723244 + lambda2bar*(-0.324336526985068 + lambda2bar*0.1107439432180572));
    }
    else {
        lnx = log( lambda2bar );
    }
    REAL8 lny = XLALSimUniversalRelation( lnx, coeffs );
    return exp(lny);
}

/* Quasi universal relation between spin-induced quadrupole and
 *   spin-induced octupole moment based on Yagi & Yunes arxiv:1608.02582;
 *    Table 2 of the review also given explicitly in NRTidalv2 paper https://arxiv.org/abs/1905.06011
*/

REAL8 XLALSimUniversalRelationSpinInducedOctupoleVSSpinInducedQuadrupole(
                                                                        REAL8 qm_def /**< spin-induced quadrupole moment */
)
{
    REAL8 coeffs[] = {0.003131, 2.071, -0.7152, 0.2458, -0.03309};
    REAL8 lnx;
    lnx = log( qm_def );

    REAL8 lny = 1;
    lny = XLALSimUniversalRelation( lnx, coeffs );
    return exp(lny);

}
