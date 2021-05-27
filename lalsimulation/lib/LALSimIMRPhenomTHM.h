#ifndef _LALSIM_IMR_PHENOMT_H
#define _LALSIM_IMR_PHENOMT_H

/*
 * Copyright (C) 2020 Hector Estelles
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
 * \author Hector Estelles
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

/*+++++++++++++++++++++++ Internal functions *****************/

/* Function for generating a specific lm mode */
int LALSimIMRPhenomTHM_OneMode(
	COMPLEX16TimeSeries **hlm,         /**< UNDOCUMENTED */
	IMRPhenomTWaveformStruct *pWF,     /**< UNDOCUMENTED */
	IMRPhenomTPhase22Struct *pPhase,   /**< UNDOCUMENTED */
	REAL8Sequence *phase22,            /**< Values of the 22 phase for the waveform time array */
	REAL8Sequence *xorb,               /**< UNDOCUMENTED */
	UINT4 ell,                         /**< UNDOCUMENTED */
	UINT4 emm                          /**< UNDOCUMENTED */
        );

/* Function for generating all modes in a requested modelist */
int LALSimIMRPhenomTHM_Modes(
  SphHarmTimeSeries **hlms,   /**< [out] Time domain modes */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1L,                /**< Dimensionless aligned spin of companion 1 */
  REAL8 chi2L,                /**< Dimensionless aligned spin of companion 2 */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams,            /**< LAL dictionary containing accessory parameters */
  UINT4 only22              /* Internal flag for mode array check */
  );

/* Function to adapt user requested mode list to modes included in the model */
LALDict *IMRPhenomTHM_setup_mode_array(LALDict *lalParams);
int check_input_mode_array_THM(LALDict *lalParams);
int check_input_mode_array_22_THM(LALDict *lalParams);


#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMT_H */