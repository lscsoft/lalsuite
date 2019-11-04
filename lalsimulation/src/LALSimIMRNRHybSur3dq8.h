/*
 * Copyright (C) 2018 Vijay Varma
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
 * \author Vijay Varma
 *
 * \file
 *
 * \brief C code for NRHybSur3dq8 waveform model, an NR-hybrid surrogate model
 * for aligned-spin BBH.
 *
 * The binary data file is available at https://dcc.ligo.org/LIGO-T1900034.
 * Get the lalsuite-extra repo or put the data into a location in your
 * LAL_DATA_PATH.
 *
 * **Paper**: https://arxiv.org/abs/1812.07865
 *
 * **Parameter ranges**:
 *
 *   q = [1, 10.1] and \f$\chi_{1z}, \chi_{2z}\f$ = [-0.81, 0.81]
 *   or
 *   q = [1, 9.1] and \f$\chi_{1z}, \chi_{2z}\f$ = [-0.91, 0.91]
 *
 *   modes: \f$ \ell \leq 4, m \geq 0 \f$, and (5,5), but not (4,1) or (4,0).
 *   m<0 modes are determined from the m \f$\geq0\f$ modes.
 *
 *   \f$M \geq 2.25 M_{\odot} \f$, for fstart=20Hz, for all modes.
 *
 * **Training parameter ranges**:
 *
 *   q = [1, 8]
 *
 *   \f$\chi_{1z}, \chi_{2z}\f$ = [-0.8, 0.8]
 *
 *   But extrapolates reasonably to the above mass ratios and spins.
 */


#ifndef _LALSIM_INSPIRAL_NRHybSur3dq8_H
#define _LALSIM_INSPIRAL_NRHybSur3dq8_H
#endif


#include <stdbool.h>

#include <gsl/gsl_vector.h>

#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALSimIMR.h>

#include "LALSimNRHybSurUtilities.h"



// https://dcc.ligo.org/LIGO-T1900034
static const char NRHybSur3dq8_DATAFILE[] = "NRHybSur3dq8_lal.h5";


//*************************************************************************/
//************************* function declarations *************************/
//*************************************************************************/

static bool NRHybSur3dq8_IsSetup(void);

static void NRHybSur3dq8_Init_LALDATA(void);

int NRHybSur3dq8_fitParams(
    gsl_vector* fit_params,
    const REAL8 q,
    const REAL8 chi1z,
    const REAL8 chi2z
);


int NRHybSur3dq8_core(
    gsl_vector **phi_22,
    EvaluatedDataPieces **evaluated_mode_dps,
    LIGOTimeGPS *epoch,
    const REAL8 deltaTOverM,
    const REAL8 fMin,
    const REAL8 fRef,
    REAL8 q,
    const REAL8 Mtot_sec,
    REAL8 chi1z,
    REAL8 chi2z,
    LALValue* ModeArray,
    LALDict* LALparams
);
