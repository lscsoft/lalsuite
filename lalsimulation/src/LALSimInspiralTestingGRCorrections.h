/*
 *  Copyright (C) 2017 Walter Del Pozzo
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

#ifndef _LALSIMINSPIRALTESTINGGRCORRECTIONS_H
#define _LALSIMINSPIRALTESTINGGRCORRECTIONS_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>

int XLALSimInspiralTestingGRCorrections(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */                                   
                                        const REAL8 m1_SI,
                                        const REAL8 m2_SI,
                                        const REAL8 chi1z,
                                        const REAL8 chi2z,
                                        const REAL8 f_low,
                                        const REAL8 f_ref,
                                        const REAL8 lambda1,
                                        const REAL8 lambda2,
					const REAL8 f_window_div_f_Peak,     /** Frequency at which to attach non-GR and GR waveforms, inputted as a fraction of f_Peak (should be between 0 and 1) */
					const REAL8 NCyclesStep,                /** Number of GW cycles over which to taper the non-GR phase correction */
                                        const LALSimInspiralTestGRParam *pnCorrections    /**< input linked list of testing gr parameters */
);

void XLALSimInspiralPNCorrections(PNPhasingSeries *pfa, const REAL8 m1, const REAL8 m2, const REAL8 chi1L, const REAL8 chi2L, const REAL8 chi1sq, const REAL8 chi2sq, const REAL8 chi1dotchi2, const REAL8 qm_def1, const REAL8 qm_def2, const LALSimInspiralTestGRParam *pnCorrections);

int XLALSimInspiralPhaseCorrectionsPhasing(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                           const REAL8Sequence *freqs,
                                           const UINT4 iStart,
                                           const UINT4 iRef,
                                           const UINT4 iEnd,
                                           const UINT4 iPeak,
                                           PNPhasingSeries pfa,
                                           const REAL8 mtot,
                                           const REAL8 eta,
					   const REAL8 NCyclesStep);  /** Choose number of GW cycles over which to taper the non-GR phase correction */

int XLALSimInspiralTestingGRCorrectionsWithDS(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                        const REAL8 distance,
                                        const REAL8 m1_SI,
                                        const REAL8 m2_SI,
                                        const REAL8 chi1z,
                                        const REAL8 chi2z,
                                        const REAL8 f_low,
                                        const REAL8 f_ref,
					const REAL8 f_window_div_f_Peak,     /** Frequency at which to attach non-GR to GR waveforms, inputted as a fraction of f_Peak (should be between 0 and 1) */
					const REAL8 NCyclesStep,             /** Number of GW cycles over which to taper the non-GR phase correction to GR waveform */
					const REAL8 f_DS,                    /** Frequency at which to attach GR to non-GR waveforms, inputted in Hz */
					const REAL8 NCyclesDS,               /** Number of GW cycles over which to taper the GR waveform to non-GR phase corrections */                 
                                        const LALSimInspiralTestGRParam *pnCorrections    /**< input linked list of testing gr parameters */
);

int XLALSimInspiralPhaseCorrectionsPhasingWithDS(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                           const REAL8 distance,
                                           const REAL8Sequence *freqs,
                                           const UINT4 iStart,
                                           const UINT4 iEnd,
                                           PNPhasingSeries pfa,
                                           const REAL8 mtot,
                                           const REAL8 eta,
			  		   const REAL8 f_ref, /** this must be in seconds **/
					   const REAL8 NCyclesStep, /** Choose number of GW cycles over which to taper the non-GR phase correction */
					   const REAL8 f_DS,  /** Frequency in which dynamical scalarization turns on */
					   const REAL8 NCyclesDS); /** Choose number of GW cycles over which to taper the GR phase correction */

REAL8 PNPhase(REAL8 f, PNPhasingSeries pfa, const REAL8 mtot);
REAL8 PNPhaseDerivative(REAL8 f, PNPhasingSeries pfa, const REAL8 mtot);
REAL8 PNPhaseSecondDerivative(REAL8 f, PNPhasingSeries pfa, const REAL8 mtot);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRALTESTINGGRCORRECTIONS_H */
