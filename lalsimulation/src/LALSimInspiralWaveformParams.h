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

int XLALSimInspiralWaveformParamsInsertModesChoice(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertFrameAxis(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertSideband(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertNumRelData(LALDict *params, const char * value);

int XLALSimInspiralWaveformParamsInsertLmax(LALDict *params, INT4 value);

int XLALSimInspiralWaveformParamsInsertPNPhaseOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPNEccentricityOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPNSpinOrder(LALDict *params, INT4 value);
int XLALSimInspiralWaveformParamsInsertPNTidalOrder(LALDict *params, INT4 value);

int XLALSimInspiralWaveformParamsInsertTidalLambda1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertTidalLambda2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertdQuadMon1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertdQuadMon2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertRedshift(LALDict *params, REAL8 value);

int XLALSimInspiralWaveformParamsInsertNonGRPhi1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRPhi2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRPhi3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRPhi4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi0(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi3(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi4(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi5(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi5L(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi6(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi6L(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNonGRDChi7(LALDict *params, REAL8 value);
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

/* NLTides parameters */
/* used within LALSimInspiralTaylorF2NLTides.c */
int XLALSimInspiralWaveformParamsInsertNLTidesA1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesN1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesF1(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesA2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesN2(LALDict *params, REAL8 value);
int XLALSimInspiralWaveformParamsInsertNLTidesF2(LALDict *params, REAL8 value);

INT4 XLALSimInspiralWaveformParamsLookupModesChoice(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupFrameAxis(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupSideband(LALDict *params);
const char * XLALSimInspiralWaveformParamsLookupNumRelData(LALDict *params);

INT4 XLALSimInspiralWaveformParamsLookupLmax(LALDict *params);

INT4 XLALSimInspiralWaveformParamsLookupPNPhaseOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPNAmplitudeOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPNEccentricityOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPNSpinOrder(LALDict *params);
INT4 XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALDict *params);

REAL8 XLALSimInspiralWaveformParamsLookupTidalLambda1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupTidalLambda2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupdQuadMon1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupdQuadMon2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupRedshift(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRPhi1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRPhi2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRPhi3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRPhi4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi0(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi3(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi4(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi5(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi5L(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi6(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi6L(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNonGRDChi7(LALDict *params);
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

/* NLTides parameters */
/* used within LALSimInspiralTaylorF2NLTides.c */
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesA1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesN1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesF1(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesA2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesN2(LALDict *params);
REAL8 XLALSimInspiralWaveformParamsLookupNLTidesF2(LALDict *params);

int XLALSimInspiralWaveformParamsModesChoiceIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsFrameAxisIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsSidebandIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNumRelDataIsDefault(LALDict *params);

int XLALSimInspiralWaveformParamsLmaxIsDefault(LALDict *params);

int XLALSimInspiralWaveformParamsPNPhaseOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPNAmplitudeOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPNEccentricityOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsPNTidalOrderIsDefault(LALDict *params);

int XLALSimInspiralWaveformParamsTidalLambda1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsTidalLambda2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsdQuadMon1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsdQuadMon2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsRedshiftIsDefault(LALDict *params);

int XLALSimInspiralWaveformParamsNonGRPhi1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRPhi2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRPhi3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRPhi4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi0IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi1IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi2IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi3IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi4IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi5IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi5LIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi6IsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi6LIsDefault(LALDict *params);
int XLALSimInspiralWaveformParamsNonGRDChi7IsDefault(LALDict *params);
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

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRALWAVEFORMPARAMS_H */
