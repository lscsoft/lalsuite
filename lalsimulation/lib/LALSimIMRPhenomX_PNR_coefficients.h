#ifndef _LALSIM_IMR_PHENOMX_PNR_COEFFICIENTS_H
#define _LALSIM_IMR_PHENOMX_PNR_COEFFICIENTS_H
/*
 * Copyright (C) 2022 Cardiff University
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

/**
 * \author Eleanor Hamilton, Sebastian Khan, Jonathan E. Thompson
 *
 */

#ifdef __cplusplus
extern "C"
{
#endif

#include <lal/LALDatatypes.h>

#include <math.h>

    REAL8 IMRPhenomX_PNR_evaluate_coefficient_array(REAL8 coeff_array[4][4][5], REAL8 eta, REAL8 chi, REAL8 costheta);

    REAL8 IMRPhenomX_PNR_alpha_A1_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);
    REAL8 IMRPhenomX_PNR_alpha_A2_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);
    REAL8 IMRPhenomX_PNR_alpha_A3_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);
    REAL8 IMRPhenomX_PNR_alpha_A4_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);

    REAL8 IMRPhenomX_PNR_beta_B0_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);
    REAL8 IMRPhenomX_PNR_beta_B1_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);
    REAL8 IMRPhenomX_PNR_beta_B2_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);
    REAL8 IMRPhenomX_PNR_beta_B3_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);
    REAL8 IMRPhenomX_PNR_beta_B4_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);
    REAL8 IMRPhenomX_PNR_beta_B5_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);

    REAL8 IMRPhenomX_PNR_beta_Bf_coefficient(REAL8 eta, REAL8 chi, REAL8 costheta);

#ifdef __cplusplus
}
#endif

#endif /*_LALSIM_IMR_PHENOMX_PNR_COEFFICIENTS_H*/