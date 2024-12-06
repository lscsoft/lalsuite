#ifndef _LALSIMINSPIRALWAVEFORMPARAMS_H
#define _LALSIMINSPIRALWAVEFORMPARAMS_H

#include <lal/LALDatatypes.h>
#include <lal/LALDict.h>
#include <lal/LALSimInspiralWaveformFlags.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/* Parameters of the New Waveforms Interface */
int XLALSimInspiralCheckKnownREAL8Key(const char* key);
int XLALSimInspiralCheckDeterminationOfMasses(LALDict *params);

int XLALSimInspiralWaveformParamsInsertMass1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertMass2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTotalMass(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertMassRatio(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSymMassRatio(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertChirpMass(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertMassDifference(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertReducedMass(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin1x(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin2x(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin1y(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin2y(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin1z(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin2z(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin1norm(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin2norm(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin1tilt(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin2tilt(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin1phi(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertSpin2phi(LALDict *params, REAL8 value);

int XLALSimInspiralWaveformParamsInsertDeltaF(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDeltaT(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertF22Ref(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertRefPhase(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertF22Start(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertFMax(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDistance(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertInclination(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertLongAscNodes(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertEccentricity(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertMeanPerAno(LALDict *params, REAL8 value);

int XLALSimInspiralWaveformParamsInsertLmax(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertModesChoice(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertFrameAxis(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertSideband(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertNumRelData(LALDict *params, const char * value);

int XLALSimInspiralWaveformParamsInsertModeArray(LALDict *params,  LALValue *value);
int XLALSimInspiralWaveformParamsInsertModeArrayJframe(LALDict *params,  LALValue *value);
int XLALSimInspiralWaveformParamsInsertModeArrayFromModeString(LALDict *params, const char *modestr);
int XLALSimInspiralWaveformParamsInsertModeArrayJframeFromModeString(LALDict *params, const char *modestr);

int XLALSimInspiralWaveformParamsInsertPNPhaseOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPNEccentricityOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPNSpinOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPNTidalOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertGETides(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertGMTides(LALDict *params, INT4 value);

int XLALSimInspiralWaveformParamsInsertTidalLambda1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalLambda2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalOctupolarLambda2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalHexadecapolarLambda1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalHexadecapolarLambda2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalQuadrupolarFMode1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalQuadrupolarFMode2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalOctupolarFMode1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalOctupolarFMode2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertLscorr(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertFinalFreq(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertOnlyFinal(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertdQuadMon1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertdQuadMon2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertRedshift(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertEccentricityFreq(LALDict *params, REAL8 value);

/* IMRPhenomX Parameters */
int XLALSimInspiralWaveformParamsInsertPhenomXInspiralPhaseVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXInspiralAmpVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXIntermediatePhaseVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXIntermediateAmpVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXRingdownPhaseVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXRingdownAmpVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPrecVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXReturnCoPrec(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPExpansionOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPConvention(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPFinalSpinMod(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPTransPrecessionMethod(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXTidalFlag(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPSpinTaylorVersion(LALDict *params, const char * value);
int XLALSimInspiralWaveformParamsInsertPhenomXPSpinTaylorCoarseFactor(LALDict *params, INT4 value);

/* IMRPhenomXHM Parameters */
int XLALSimInspiralWaveformParamsInsertPhenomXHMReleaseVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMInspiralPhaseVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMIntermediatePhaseVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMRingdownPhaseVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMInspiralAmpVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMIntermediateAmpVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMRingdownAmpVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMInspiralAmpFitsVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMIntermediateAmpFitsVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMRingdownAmpFitsVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMInspiralAmpFreqsVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMIntermediateAmpFreqsVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMRingdownAmpFreqsVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMPhaseRef21(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMThresholdMband(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXHMAmpInterpolMB(LALDict *params, INT4 value);

/* IMRPhenomTHM Parameters */
int XLALSimInspiralWaveformParamsInsertPhenomTHMInspiralVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomTPHMMergerVersion(LALDict *params, INT4 value);

/* IMRPhenomX_PNR Parameters */
int XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedAngles(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedCoprec(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedCoprec33(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPNRUseInputCoprecDeviations(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPNRForceXHMAlignment(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPNRInterpTolerance(LALDict *params, REAL8 value);

/* IMRPhenomX_PNR anti-symmetric waveform Parameters */
int XLALSimInspiralWaveformParamsInsertPhenomXAntisymmetricWaveform(LALDict *params, INT4 value);


/* IMRPhenomXPHM Parameters */
int XLALSimInspiralWaveformParamsInsertPhenomXPHMMBandVersion(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPHMUseModes(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPHMModesL0Frame(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPHMPrecModes(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPhenomXPHMTwistPhenomHM(LALDict *params, INT4 value);

/* IMRPhenomXCP Parameters */
int XLALSimInspiralWaveformParamsInsertPhenomXCPMU1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPMU2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPMU3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPMU4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPNU0(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPNU4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPNU5(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPNU6(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPZETA1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPZETA2(LALDict *params, REAL8 value);
/* l=3, m=3 */
int XLALSimInspiralWaveformParamsInsertPhenomXCPMU1l3m3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPMU2l3m3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPMU3l3m3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPMU4l3m3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPNU0l3m3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPNU4l3m3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPNU5l3m3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPNU6l3m3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPZETA1l3m3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertPhenomXCPZETA2l3m3(LALDict *params, REAL8 value);
/* Toggle output of XAS phase for debugging purposes */
int XLALSimInspiralWaveformParamsInsertPhenomXOnlyReturnPhase(LALDict *params, INT4 value);

int XLALSimInspiralWaveformParamsInsertNonGRPhi1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRPhi2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRPhi3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRPhi4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChiMinus2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChiMinus1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi0(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi3S(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi3NS(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi4S(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi4NS(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi5(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi5S(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi5NS(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi5L(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi5LS(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi5LNS(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi6(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi6S(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi6NS(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi6L(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi7(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi7S(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi7NS(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDXi1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDXi2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDXi3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDXi4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDXi5(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDXi6(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDSigma1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDSigma2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDSigma3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDSigma4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDAlpha1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDAlpha2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDAlpha3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDAlpha4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDAlpha5(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDBeta1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDBeta2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDBeta3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRBetaPPE(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE0(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRBetaPPE0(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRBetaPPE1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRBetaPPE2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRBetaPPE3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRBetaPPE4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE5(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRBetaPPE5(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE6(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRBetaPPE6(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE7(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRBetaPPE7(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertEnableLIV(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertNonGRLIVLogLambdaEff(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRLIVASign(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRLIVAlpha(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChikappaS(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChikappaA(LALDict *params, REAL8 value);

/* PhenomX TGR Coefficients */
int XLALSimInspiralWaveformParamsInsertNonGRDC1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDC2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDC4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDCL(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDB1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDB2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDB3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDB4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRParameterization(LALDict *params, INT4 value);


/* NLTides parameters */
/* used within LALSimInspiralTaylorF2NLTides.c */
int XLALSimInspiralWaveformParamsInsertNLTidesA1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesN1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesF1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesA2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesN2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesF2(LALDict *params, REAL8 value);

/* SEOBNRv4P */
INT4 XLALSimInspiralWaveformParamsInsertEOBChooseNumOrAnalHamDer(LALDict *params, INT4 value);
INT4 XLALSimInspiralWaveformParamsInsertEOBEllMaxForNyquistCheck(LALDict *params, INT4 value);


/* new interface */
REAL8 XLALSimInspiralWaveformParamsLookupMass1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupMass2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTotalMass(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupMassRatio(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSymMassRatio(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupChirpMass(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupMassDifference(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupReducedMass(LALDict *params);

REAL8 XLALSimInspiralGetMassRatioFromSymMassRatio(REAL8 sym_mass_ratio);
REAL8 XLALSimInspiralGetMassRatioFromChirpMassComponentMass1(REAL8 chirp_mass, REAL8 component_mass);
REAL8 XLALSimInspiralGetMassRatioFromChirpMassComponentMass2(REAL8 chirp_mass, REAL8 component_mass);

REAL8 XLALSimInspiralGetCartesianSpinXFromPolar(REAL8 spin_norm, REAL8 spin_tilt, REAL8 spin_phi);
REAL8 XLALSimInspiralGetCartesianSpinYFromPolar(REAL8 spin_norm, REAL8 spin_tilt, REAL8 spin_phi);
REAL8 XLALSimInspiralGetCartesianSpinZFromPolar(REAL8 spin_norm, REAL8 spin_tilt);

REAL8 XLALSimInspiralGetPolarSpin_normFromCartesian(REAL8 spinx, REAL8 spiny, REAL8 spinz);
REAL8 XLALSimInspiralGetPolarSpin_tiltFromCartesian(REAL8 spinx, REAL8 spiny, REAL8 spinz);
REAL8 XLALSimInspiralGetPolarSpin_phiFromCartesian(REAL8 spiny, REAL8 spinz);


REAL8 XLALSimInspiralWaveformParamsLookupSpin1x(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin1y(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin1z(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin2x(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin2y(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin2z(LALDict *params);

REAL8 XLALSimInspiralWaveformParamsLookupSpin1norm(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin1tilt(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin1phi(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin2norm(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin2tilt(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupSpin2phi(LALDict *params);

REAL8 XLALSimInspiralWaveformParamsLookupDeltaF(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDeltaT(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupF22Ref(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupRefPhase(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupF22Start(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupFMax(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDistance(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupInclination(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupLongAscNodes(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupEccentricity(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupMeanPerAno(LALDict *params);

INT4 XLALSimInspiralWaveformParamsLookupLmax(LALDict *params);

int XLALSimInspiralWaveformParamsInsertDOmega220(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDTau220(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDOmega210(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDTau210(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDOmega330(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDTau330(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDOmega440(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDTau440(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDOmega550(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertDTau550(LALDict *params, REAL8 value);

INT4 XLALSimInspiralWaveformParamsLookupModesChoice(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupFrameAxis(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupSideband(LALDict *params);
const char * XLALSimInspiralWaveformParamsLookupNumRelData(LALDict *params);

LALValue* XLALSimInspiralWaveformParamsLookupModeArray(LALDict *params);
LALValue* XLALSimInspiralWaveformParamsLookupModeArrayJframe(LALDict *params);

INT4 XLALSimInspiralWaveformParamsLookupPNPhaseOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPNAmplitudeOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPNEccentricityOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPNSpinOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupGETides(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupGMTides(LALDict *params);

REAL8 XLALSimInspiralWaveformParamsLookupTidalLambda1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTidalLambda2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTidalOctupolarLambda1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTidalOctupolarLambda2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTidalHexadecapolarLambda1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTidalHexadecapolarLambda2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTidalQuadrupolarFMode1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTidalQuadrupolarFMode2(LALDict *params);;
REAL8 XLALSimInspiralWaveformParamsLookupTidalOctupolarFMode1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTidalOctupolarFMode2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupdQuadMon1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupdQuadMon2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupRedshift(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupEccentricityFreq(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupLscorr(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupFinalFreq(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupOnlyFinal(LALDict *params);

/* IMRPhenomX Parameters */
INT4 XLALSimInspiralWaveformParamsLookupPhenomXInspiralPhaseVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXInspiralAmpVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXIntermediatePhaseVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXIntermediateAmpVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXRingdownPhaseVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXRingdownAmpVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXReturnCoPrec(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPExpansionOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPConvention(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPFinalSpinMod(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXTidalFlag(LALDict *params);
const char * XLALSimInspiralWaveformParamsLookupPhenomXPSpinTaylorVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPSpinTaylorCoarseFactor(LALDict *params);

/* IMRPhenomXHM Parameters */
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMReleaseVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMInspiralPhaseVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMIntermediatePhaseVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMRingdownPhaseVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMInspiralAmpVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMIntermediateAmpVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMRingdownAmpVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMInspiralAmpFitsVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMIntermediateAmpFitsVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMRingdownAmpFitsVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMInspiralAmpFreqsVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMIntermediateAmpFreqsVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMRingdownAmpFreqsVersion(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXHMPhaseRef21(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXHMThresholdMband(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXHMAmpInterpolMB(LALDict *params);

/* IMRPhenomTHM Parameters */
INT4 XLALSimInspiralWaveformParamsLookupPhenomTHMInspiralVersion(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomTPHMMergerVersion(LALDict *params);

/* IMRPhenomX_PNR Parameters */
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedCoprec(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedCoprec33(LALDict *params);
// Option to only be used when actively tuning PNR Coprec relative to XHM wherein the non-precessing final spin is used
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPNRUseInputCoprecDeviations(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPNRForceXHMAlignment(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXPNRInterpTolerance(LALDict *params);

INT4 XLALSimInspiralWaveformParamsLookupPhenomXAntisymmetricWaveform(LALDict *params);


/* IMRPhenomXPHM Parameters */
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPHMMBandVersion(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXPHMThresholdMband(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPHMUseModes(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPHMModesL0Frame(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPHMPrecModes(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPHMTwistPhenomHM(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPhenomXPTransPrecessionMethod(LALDict *params);

/* IMRPhenomXCP Parameters */
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPMU1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPMU2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPMU3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPMU4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPNU0(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPNU4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPNU5(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPNU6(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPZETA1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPZETA2(LALDict *params);
/* l=3, m=3 */
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPMU1l3m3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPMU2l3m3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPMU3l3m3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPMU4l3m3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPNU0l3m3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPNU4l3m3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPNU5l3m3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPNU6l3m3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPZETA1l3m3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupPhenomXCPZETA2l3m3(LALDict *params);
/* Toggle output of XAS phase for debugging purposes */
int XLALSimInspiralWaveformParamsLookupPhenomXOnlyReturnPhase(LALDict *params);

/* NonGR Parameters */
REAL8 XLALSimInspiralWaveformParamsLookupNonGRPhi1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRPhi2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRPhi3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRPhi4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChiMinus2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChiMinus1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi0(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi3S(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi3NS(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi4S(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi4NS(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi5(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi5S(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi5NS(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi5L(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi5LS(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi5LNS(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi6(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi6S(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi6NS(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi6L(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi7(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi7S(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi7NS(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDXi1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDXi2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDXi3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDXi4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDXi5(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDXi6(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDSigma1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDSigma2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDSigma3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDSigma4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDAlpha1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDAlpha2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDAlpha3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDAlpha4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDAlpha5(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDBeta1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDBeta2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDBeta3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRAlphaPPE(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRBetaPPE(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRAlphaPPE0(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRBetaPPE0(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRAlphaPPE1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRBetaPPE1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRAlphaPPE2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRBetaPPE2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRAlphaPPE3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRBetaPPE3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRAlphaPPE4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRBetaPPE4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRAlphaPPE5(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRBetaPPE5(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRAlphaPPE6(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRBetaPPE6(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRAlphaPPE7(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRBetaPPE7(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupEnableLIV(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRLIVLogLambdaEff(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRLIVASign(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRLIVAlpha(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChikappaS(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChikappaA(LALDict *params);

/* PhenomX TGR Coefficients */
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDC1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDC2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDC4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDCL(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDB1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDB2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDB3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDB4(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupNonGRParameterization(LALDict *params);


/* NLTides parameters */
/* used within LALSimInspiralTaylorF2NLTides.c */
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesA1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesN1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesF1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesA2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesN2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesF2(LALDict *params);

/* SEOBNRv4P */
INT4 XLALSimInspiralWaveformParamsLookupEOBChooseNumOrAnalHamDer(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupEOBEllMaxForNyquistCheck(LALDict *params);

int XLALSimInspiralWaveformParamsLmaxIsDefault(LALDict *params);

REAL8 XLALSimInspiralWaveformParamsLookupDOmega220(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDTau220(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDOmega210(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDTau210(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDOmega330(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDTau330(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDOmega440(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDTau440(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDOmega550(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupDTau550(LALDict *params);

int XLALSimInspiralWaveformParamsModesChoiceIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsFrameAxisIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsSidebandIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNumRelDataIsDefault(LALDict *params);

int XLALSimInspiralWaveformParamsModeArrayIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsModeArrayJframeIsDefault(LALDict *params);

int XLALSimInspiralWaveformParamsPNPhaseOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPNAmplitudeOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPNEccentricityOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPNTidalOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsGETidesIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsGMTidesIsDefault(LALDict *params);

int XLALSimInspiralWaveformParamsTidalLambda1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalLambda2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalOctupolarLambda1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalOctupolarLambda2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalHexadecapolarLambda1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalHexadecapolarLambda2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalQuadrupolarFMode1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalQuadrupolarFMode2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalOctupolarFMode1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalOctupolarFMode2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsdQuadMon1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsdQuadMon2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsRedshiftIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsEccentricityFreqIsDefault(LALDict *params);

/* IMRPhenomX Parameters */
int XLALSimInspiralWaveformParamsPhenomXInspiralPhaseVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXInspiralAmpVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXIntermediatePhaseVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXIntermediateAmpVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXRingdownPhaseVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXRingdownAmpVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPrecVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXReturnCoPrecIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPExpansionOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPConventionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPFinalSpinModIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPTransPrecessionMethodIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXTidalFlagIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPSpinTaylorVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPSpinTaylorCoarseFactorIsDefault(LALDict *params);

/* IMRPhenomXHM Parameters */
int XLALSimInspiralWaveformParamsPhenomXHMReleaseVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMInspiralPhaseVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMIntermediatePhaseVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMRingdownPhaseVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMInspiralAmpVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMIntermediateAmpVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMRingdownAmpVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMInspiralAmpFitsVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMIntermediateAmpFitsVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMRingdownAmpFitsVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMInspiralAmpFreqsVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMIntermediateAmpFreqsVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMRingdownAmpFreqsVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMPhaseRef21IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMThresholdMbandIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXHMAmpInterpolMBIsDefault(LALDict *params);

/* IMRPhenomXPHM Parameters */
int XLALSimInspiralWaveformParamsPhenomXPHMMBandVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPHMThresholdMbandIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPHMUseModesIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPHMModesL0FrameIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPHMPrecModesIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPHMTwistPhenomHMIsDefault(LALDict *params);

/* IMRPhenomTHM Parameters */
int XLALSimInspiralWaveformParamsPhenomTHMInspiralVersionIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomTPHMMergerVersionIsDefault(LALDict *params);

/* IMRPhenomX_PNR Parameters */
int XLALSimInspiralWaveformParamsPhenomXPNRUseTunedAnglesIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPNRUseTunedCoprecIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPNRUseTunedCoprec33IsDefault(LALDict *params);
// Option to only be used when actively tuning PNR Coprec relative to XHM wherein the non-precessing final spin is used
int XLALSimInspiralWaveformParamsPhenomXPNRUseInputCoprecDeviationsIsDefault(LALDict *params);
// Dev option for forcing 22 phase derivative inspiral values to align with XHM at a low ref frequency
int XLALSimInspiralWaveformParamsPhenomXPNRForceXHMAlignmentIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXPNRInterpToleranceIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXAntisymmetricWaveformIsDefault(LALDict *params);


/* IMRPhenomXCP Parameters */
int XLALSimInspiralWaveformParamsPhenomXCPMU1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPMU2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPMU3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPMU4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPNU0IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPNU4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPNU5IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPNU6IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPZETA1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPZETA2IsDefault(LALDict *params);
/* l=3, m=3 */
int XLALSimInspiralWaveformParamsPhenomXCPMU1l3m3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPMU2l3m3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPMU3l3m3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPMU4l3m3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPNU0l3m3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPNU4l3m3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPNU5l3m3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPNU6l3m3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPZETA1l3m3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPhenomXCPZETA2l3m3IsDefault(LALDict *params);
/* Toggle output of XAS phase for debugging purposes */
int XLALSimInspiralWaveformParamsPhenomXOnlyReturnPhaseIsDefault(LALDict *params);

/* NonGR Parameters */
int XLALSimInspiralWaveformParamsNonGRPhi1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRPhi2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRPhi3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRPhi4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChiMinus2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChiMinus1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi0IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi3SIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi3NSIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi4SIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi4NSIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi5IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi5SIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi5NSIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi5LIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi5LSIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi5LNSIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi6IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi6SIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi6NSIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi6LIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi7IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi7SIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi7NSIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDXi1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDXi2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDXi3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDXi4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDXi5IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDXi6IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDSigma1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDSigma2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDSigma3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDSigma4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDAlpha1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDAlpha2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDAlpha3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDAlpha4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDAlpha5IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDBeta1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDBeta2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDBeta3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRAlphaPPEIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRBetaPPEIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRAlphaPPE0IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRBetaPPE0IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRAlphaPPE1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRBetaPPE1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRAlphaPPE2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRBetaPPE2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRAlphaPPE3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRBetaPPE3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRAlphaPPE4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRBetaPPE4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRAlphaPPE5IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRBetaPPE5IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRAlphaPPE6IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRBetaPPE6IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRAlphaPPE7IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRBetaPPE7IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsEnableLIVIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRLIVLogLambdaEffIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRLIVASignIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRLIVAlphaIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChikappaSIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChikappaAIsDefault(LALDict *params);

/* PhenomX TGR Coefficients */
int XLALSimInspiralWaveformParamsNonGRDC1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDC2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDC4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDCLIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDB1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDB2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDB3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDB4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRParameterizationIsDefault(LALDict *params);


/* SEOBNRv4P */
INT4 XLALSimInspiralWaveformParamsEOBChooseNumOrAnalHamDerIsDefault(LALDict *params);
INT4 XLALSimInspiralWaveformParamsEOBEllMaxForNyquistCheckIsDefault(LALDict *params);

LALDict* XLALSimInspiralParamsDict(const REAL8 m1, const REAL8 m2, const REAL8 S1x, const REAL8 S1y, const REAL8 S1z, const REAL8 S2x, const REAL8 S2y, const REAL8 S2z, const REAL8 distance, const REAL8 inclination, const REAL8 phiRef, const REAL8 longAscNodes, const REAL8 eccentricity, const REAL8 f_ref, LALDict *LALparams);

int XLALSimInspiralWaveformParamsDOmega220IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsDTau220IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsDOmega210IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsDTau210IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsDOmega330IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsDTau330IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsDOmega440IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsDTau440IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsDOmega550IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsDTau550IsDefault(LALDict *params);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRALWAVEFORMPARAMS_H */
