/*  Copyright (C) 2019 Vijay Varma
 *  Evaluates NRSur7dq4Remnant model for remnant BH mass, spin and recoil kick
 *  for generically precessing BBH.
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
 * \author Vijay Varma
 *
 * \file
 *
 * \brief NRSur7dq4Remnant model for remnant BH mass, spin and recoil kick for
 * generically precessing BBH.
 *
 * The binary data file (NRSur7dq4Remnant_v1.0.h5) is available at:
 * https://git.ligo.org/waveforms/software/lalsuite-waveform-data
 * and on Zenodo at: https://zenodo.org/records/14999310.
 * Get the lalsuite-waveform-data repo or put the data into a location in your
 * LAL_DATA_PATH.
 * The data is also available on CIT at /home/lalsimulation_data and via CVMFS
 * at /cvmfs/shared.storage.igwn.org/igwn/shared/auxiliary/obs_sci/cbc/waveform/lalsimulation_data
 *
 * **Paper**: https://arxiv.org/abs/1905.09300,
 * https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.1.033015
 *
 * **Parameter ranges**:
 *
 *   q = [1, 6.01]
 *
 *   \f$|\chi_{1}|, |\chi_{2}| \leq 1 \f$
 *
 * **Training parameter ranges**:
 *
 *   q = [1, 4.01]
 *
 *   \f$|\chi_{1}|, |\chi_{2}| \leq 0.81 \f$
 *
 *   But extrapolates reasonably to the above mass ratios and spins. However,
 *   if a guarantee of accuracy is required, this model should be used within
 *   the training parameter range.
 */

// https://git.ligo.org/waveforms/software/lalsuite-waveform-data, should be placed in $LAL_DATA_PATH.
// Also on Zenodo at: https://zenodo.org/records/14999310.
// The data is also available on CIT at /home/lalsimulation_data and via CVMFS
// at /cvmfs/shared.storage.igwn.org/igwn/shared/auxiliary/obs_sci/cbc/waveform/lalsimulation_data
static const char NRSur7dq4Remnant_DATAFILE[] = "NRSur7dq4Remnant_v1.0.h5";

//*************************************************************************/
//************************* function declarations *************************/
//*************************************************************************/

static bool NRSur7dq4Remnant_IsSetup(void);

static void NRSur7dq4Remnant_Init_LALDATA(void);

static int NRSur7dq4Remnant_fitParams(
    gsl_vector* fit_params,
    const REAL8 q,
    const REAL8 chiAx,
    const REAL8 chiAy,
    const REAL8 chiAz,
    const REAL8 chiBx,
    const REAL8 chiBy,
    const REAL8 chiBz,
    LALDict* LALparams
);
