#ifndef _LALSIM_IMR_PHENOMNSBH_H
#define _LALSIM_IMR_PHENOMNSBH_H

#ifdef __cplusplus
extern "C"
{
#endif

/*
 * Copyright (C) 2019 Jonathan Thompson, Edward Fauchon-Jones, Sebastian Khan
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

#include <math.h>
#include <lal/LALSimIMR.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>


#include "LALSimIMRPhenomInternalUtils.h"
#include "LALSimIMRPhenomUtils.h"
#include "LALSimIMRPhenomC_internals.h"

typedef enum tagMergerType
{
    DISRUPTIVE,
    MILDLY_DISRUPTIVE_NO_TORUS_REMNANT,
    MILDLY_DISRUPTIVE_TORUS_REMNANT,
    NON_DISRUPTIVE,
    BBH
} MergerType;

typedef struct tagBBHPhenomNSBHParams
{
    REAL8 lambda;
    REAL8 m_sec;
    REAL8 gamma_correction;
    REAL8 delta_2_prime;
    REAL8 sigma;
    REAL8 sigma_tide;
    REAL8 d_param;
    REAL8 epsilon_ins;
    REAL8 epsilon_tide;
    REAL8 f0_tilde_PN;
    REAL8 f0_tilde_PM;
    REAL8 f0_tilde_RD;
    REAL8 chif;
    REAL8 f_RD;
    REAL8 q_factor;
    REAL8 f_tide;
    MergerType merger_type;
    REAL8 Mtorus;
    REAL8 C;
    REAL8 final_mass;
} BBHPhenomNSBHParams;

static int IMRPhenomNSBH_Core(
    COMPLEX16FrequencySeries **htilde,
    REAL8 phiRef,
    REAL8 fRef,
    REAL8 distance,
    REAL8 mBH_SI,
    REAL8 mNS_SI,
    REAL8 chi_BH,
    REAL8 chi_NS,
    LALDict *extraParams,
    const REAL8Sequence *freqs_in,
    REAL8 deltaF
);

static BBHPhenomNSBHParams *ComputeIMRPhenomNSBHParams(
    const REAL8 m1,
    const REAL8 m2,
    const REAL8 chi,
    const REAL8 lambda,
    const BBHPhenomCParams *params);

static REAL8 PhenomNSBHAmplitudeOneFrequency(const REAL8 Mf, const BBHPhenomCParams *params, const BBHPhenomNSBHParams *params_NSBH);



#ifdef __cplusplus
}
#endif

#endif
