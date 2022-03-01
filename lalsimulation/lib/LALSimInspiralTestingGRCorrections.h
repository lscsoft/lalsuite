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


/* Main function for waveform generation. Accepts a GR waveform htilde, binary parameters, testing-GR parameters, and paremeters 
 * that determine how the correction will be tapered. Calls PNCorrections to compute non-GR corrections to phase, then
 * PhaseCorrectionsPhasing to smoothly taper the correction to the baseline GR waveform.
 */
int XLALSimInspiralTestingGRCorrections(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                        const UINT4 l,
                                        const UINT4 m,                                   
                                        const REAL8 m1_SI,
                                        const REAL8 m2_SI,
                                        const REAL8 chi1z,
                                        const REAL8 chi2z,
                                        const REAL8 f_low,
                                        const REAL8 f_ref,
					const REAL8 f_window_div_f_Peak,     /** Frequency at which to attach non-GR and GR waveforms, inputted as a fraction of f_Peak (should be between 0 and 1) */
					const REAL8 NCyclesStep,                /** Number of GW cycles over which to taper the non-GR phase correction */
                                        LALDict *LALpars    /**< input linked list of testing gr parameters */
);

/* Accepts binary parameters and testing-GR parameters, computes the PN phase corrections specified by the testing-GR parameters and stores in PNPhasingSeries pfa */
void XLALSimInspiralPNCorrections(PNPhasingSeries *pfa, const REAL8 m1, const REAL8 m2, const REAL8 chi1L, const REAL8 chi2L, const REAL8 chi1sq, const REAL8 chi2sq, const REAL8 chi1dotchi2, const REAL8 qm_def1, const REAL8 qm_def2, LALDict *LALpars);

/* Accepts GR baseline waveform htilde, PN phase corrections pfa and parameters that determine how the correction will be tapered.
 * Tapers the phase correction by multiplying the second derivative of the correction w.r.t. frequency by a Heaviside function
 * and then integrating back to recover the phase. Finally, adds the phaes correction to the waveform and stores in htilde.
 */
int XLALSimInspiralPhaseCorrectionsPhasing(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                           const REAL8Sequence *freqs,
                                           const UINT4 m,
                                           const UINT4 iStart,
                                           const UINT4 iRef,
                                           const UINT4 iEnd,
                                           const UINT4 iPeak,
                                           PNPhasingSeries pfa,
                                           const REAL8 mtot,
                                           const REAL8 eta,
					   const REAL8 NCyclesStep);  /** Choose number of GW cycles over which to taper the non-GR phase correction */

/* Main function for waveform generation for signals with dynamical scalarization. Accepts a GR waveform htilde, binary parameters,
 * testing-GR parameters, and paremeters that determine how the corrections will be tapered. Calls PNCorrections to compute 
 * non-GR corrections to phase, then PhaseCorrectionsPhasing to smoothly taper the correction to the baseline GR waveform.
 */

REAL8 PNPhase(REAL8 f, UINT4 m, PNPhasingSeries pfa, const REAL8 mtot); /* Returns phase computed from PN coefficients pfa at frequency f*/
REAL8 PNPhaseDerivative(REAL8 f, UINT4 m, PNPhasingSeries pfa, const REAL8 mtot); /* Returns derivative of phase w.r.t. frequency computed from PN coefficients pfa at frequency f*/
REAL8 PNPhaseSecondDerivative(REAL8 f, UINT4 m, PNPhasingSeries pfa, const REAL8 mtot); /* Returns second derivative of phase w.r.t. frequency computed from PN coefficients pfa at frequency f*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRALTESTINGGRCORRECTIONS_H */
