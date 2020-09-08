/*
 * Copyright (C) 2019 Marta Colleoni, Cecilio Garcia Quiros
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
//
//  LALSimIMRPhenomXHM_internals.h
//
//
//  Created by Marta on 06/02/2019.
//


#ifndef _LALSIM_IMR_PHENOMXHM_INTERNALS_H
#define _LALSIM_IMR_PHENOMXHM_INTERNALS_H


#ifdef __cplusplus
extern "C" {
  #endif
  #ifdef __GNUC__
  #define UNUSED __attribute__ ((unused))
  #else
  #define UNUSED
  #endif


  #include <stdlib.h>
  #include <stdio.h>
  #include <math.h>
  #include <complex.h>

  #include <lal/LALStdlib.h>
  #include <lal/LALConstants.h>
  #include <lal/Date.h>
  #include <lal/FrequencySeries.h>
  #include <lal/Units.h>

  #include "LALSimIMRPhenomXHM_structs.h"
  #include "LALSimIMRPhenomX_internals.h"

  #include <lal/LALSimIMR.h>

  #define NMAX_INSPIRAL_COEFFICIENTS 13

  //  You should not declare static functions here, since this file is included in other files apart form the source one.

  /*********** Useful Powers of pi **************/
  extern IMRPhenomX_UsefulPowers powers_of_lalpiHM;

  /**************** QNMs and mixing coefficients ************** */
  void IMRPhenomXHM_Initialize_QNMs(QNMFits *qnmsFits);
  void IMRPhenomXHM_Initialize_MixingCoeffs(IMRPhenomXHMWaveformStruct *wf,IMRPhenomXWaveformStruct *wf22);

  /****************** Initialization of higher-mode waveform struct ******************/
  void IMRPhenomXHM_SetHMWaveformVariables(int ell, int emm, IMRPhenomXHMWaveformStruct *wf, IMRPhenomXWaveformStruct *wf22, QNMFits *qnms,LALDict *LALParams);


  /*************** AMPLITUDE ****************/

  /***************** set pointers to par-space fits *******************/
  void IMRPhenomXHM_FillAmpFitsArray(IMRPhenomXHMAmpCoefficients *pAmp);

  /***************** Cutting frequencies for amplitude ******************************/
  double IMRPhenomXHM_Amplitude_fcutInsp(IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22);
  double IMRPhenomXHM_Amplitude_fcutRD(IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22);

  /***************** Amplitude coefficients ****************/
  void IMRPhenomXHM_GetAmplitudeCoefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22);
  void IMRPhenomXHM_GetPNAmplitudeCoefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXWaveformStruct *pWF22);
  void Get21PNAmplitudeCoefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXWaveformStruct *pWF22);

  /**************** Amplitude reconstruction ******************/
  // Functions that return suitable ansatz at a given frequency
  double IMRPhenomXHM_Amplitude_noModeMixing(double f, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMWaveformStruct *pWF);
  double IMRPhenomXHM_Amplitude_ModeMixing(double f, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF,  IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22);
  double IMRPhenomXHM_Amplitude_ModeMixingRecycle(double f, IMRPhenomX_UsefulPowers *powers_of_f, COMPLEX16 wf22, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF);



  /***************** PHASE *******************/

  /***************** set pointers to par-space fits *******************/
  void IMRPhenomXHM_FillPhaseFitsArray(IMRPhenomXHMPhaseCoefficients *pPhase);
  void IMRPhenomXHM_Intermediate_CollocPtsFreqs(IMRPhenomXHMPhaseCoefficients *pPhase,IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22);

  /***************** intermediate region reconstruction utilities *******************/
  void EquidistantNodes(double nodes[], double fmin, double fmax, int npts);
  double GetfcutInsp(IMRPhenomXWaveformStruct *pWF22, IMRPhenomXHMWaveformStruct *pWFHM);

  /***************** spheroidal->spherical harmonic conversion  *******************/
  double complex SpheroidalToSpherical(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXHMAmpCoefficients *pAmplm, IMRPhenomXHMPhaseCoefficients *pPhaselm, IMRPhenomXHMWaveformStruct *pWFlm, IMRPhenomXWaveformStruct *pWF22);
  double complex SpheroidalToSphericalRecycle(double ff, IMRPhenomX_UsefulPowers *powers_of_f, COMPLEX16 wf22, IMRPhenomXHMAmpCoefficients *pAmplm, IMRPhenomXHMPhaseCoefficients *pPhaselm, IMRPhenomXHMWaveformStruct *pWFlm);

  /***************** spheroidal-harmonic ringdown reconstruction ********************/
  void IMRPhenomXHM_Ringdown_CollocPtsFreqs(IMRPhenomXHMPhaseCoefficients *pPhase,IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22);
  void GetSpheroidalCoefficients(IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22);

  /***************** phase coefficients *****************/
  int IMRPhenomXHM_PN21AmpSign (double ff,IMRPhenomXWaveformStruct *wf22);
  void IMRPhenomXHM_GetPhaseCoefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22, LALDict *lalParams);

  /**************** Phase reconstruction ******************/
  // functions that return suitable ansatz at a given frequency
  double IMRPhenomXHM_Phase_noModeMixing(double f, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF,IMRPhenomXWaveformStruct *pWF22);
  double IMRPhenomXHM_Phase_ModeMixing(double f, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF,IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22);
  double IMRPhenomXHM_Phase_ModeMixingRecycle(double f, IMRPhenomX_UsefulPowers *powers_of_f, COMPLEX16 wf22, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF);


  // Debugging function
  int ParametersToFile(
    IMRPhenomXWaveformStruct *pWF,           /**< Wf structure for the 22 mode*/
    IMRPhenomXHMWaveformStruct *pWFHM,         /**< Wf structure for the lm mode*/
    IMRPhenomXHMAmpCoefficients *pAmp,         /**< Coefficients struct of the lm Amplitude */
    UNUSED IMRPhenomXHMPhaseCoefficients *pPhase      /**< Coefficients struct of the lm Phase */
  );

#ifdef __cplusplus
}
#endif


#endif /* LALSimIMRPhenomXHM_internals_h */
