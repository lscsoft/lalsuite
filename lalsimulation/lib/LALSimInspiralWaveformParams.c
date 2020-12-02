#include <lal/LALStdio.h>
#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralWaveformParams.h>

#if 1 /* generate definitions for source */

#define DEFINE_INSERT_FUNC(NAME, TYPE, KEY, DEFAULT) \
	int XLALSimInspiralWaveformParamsInsert ## NAME(LALDict *params, TYPE value) \
	{ \
		return XLALDictInsert ## TYPE ## Value(params, KEY, value); \
	}

#define DEFINE_LOOKUP_FUNC(NAME, TYPE, KEY, DEFAULT) \
	TYPE XLALSimInspiralWaveformParamsLookup ## NAME(LALDict *params) \
	{ \
		TYPE value = DEFAULT; \
		if (params && XLALDictContains(params, KEY)) \
			value = XLALDictLookup ## TYPE ## Value(params, KEY); \
		return value; \
	}

#define DEFINE_ISDEFAULT_FUNC(NAME, TYPE, KEY, DEFAULT) \
	int XLALSimInspiralWaveformParams ## NAME ## IsDefault(LALDict *params) \
	{ \
		return XLALSimInspiralWaveformParamsLookup ## NAME(params) == DEFAULT; \
	}

#else /* generate prototypes for header */

#define DEFINE_INSERT_FUNC(NAME, TYPE, KEY, DEFAULT) \
	int XLALSimInspiralWaveformParamsInsert ## NAME(LALDict *params, TYPE value);

#define DEFINE_LOOKUP_FUNC(NAME, TYPE, KEY, DEFAULT) \
	TYPE XLALSimInspiralWaveformParamsLookup ## NAME(LALDict *params);

#define DEFINE_ISDEFAULT_FUNC(NAME, TYPE, KEY, DEFAULT) \
	int XLALSimInspiralWaveformParams ## NAME ## IsDefault(LALDict *params);

#endif

/* "String" is function names becomes type "const char *" */
#ifdef String
#undef String
#endif
#define String const char *

/*
 * Note: missing one type of data for SpinTaylorF2:
 * DEFINE_INSERT_FUNC(PNSideband, INT4, "sideband", 0)
 */

/* INSERT FUNCTIONS */

DEFINE_INSERT_FUNC(ModesChoice, INT4, "modes", LAL_SIM_INSPIRAL_MODES_CHOICE_ALL)
DEFINE_INSERT_FUNC(FrameAxis, INT4, "axis", LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L)
DEFINE_INSERT_FUNC(Sideband, INT4, "sideband", 0)
DEFINE_INSERT_FUNC(NumRelData, String, "numreldata", NULL)

int XLALSimInspiralWaveformParamsInsertModeArray(LALDict *params, LALValue *value)
{
	return XLALDictInsertValue(params, "ModeArray", value);
}

DEFINE_INSERT_FUNC(PNPhaseOrder, INT4, "phaseO", -1)
DEFINE_INSERT_FUNC(PNAmplitudeOrder, INT4, "ampO", -1)
DEFINE_INSERT_FUNC(PNEccentricityOrder, INT4, "eccO", -1)
DEFINE_INSERT_FUNC(PNSpinOrder, INT4, "spinO", -1)
DEFINE_INSERT_FUNC(PNTidalOrder, INT4, "tideO", -1)
DEFINE_INSERT_FUNC(GETides, INT4, "GEtideO", 0)
DEFINE_INSERT_FUNC(GMTides, INT4, "GMtideO", 0)

DEFINE_INSERT_FUNC(TidalLambda1, REAL8, "lambda1", 0)
DEFINE_INSERT_FUNC(TidalLambda2, REAL8, "lambda2", 0)
DEFINE_INSERT_FUNC(TidalOctupolarLambda1, REAL8, "TidalOctupolarLambda1", 0)
DEFINE_INSERT_FUNC(TidalOctupolarLambda2, REAL8, "TidalOctupolarLambda2", 0)
DEFINE_INSERT_FUNC(TidalHexadecapolarLambda1, REAL8, "TidalHexadecapolarLambda1", 0)
DEFINE_INSERT_FUNC(TidalHexadecapolarLambda2, REAL8, "TidalHexadecapolarLambda2", 0)
DEFINE_INSERT_FUNC(TidalQuadrupolarFMode1, REAL8, "TidalQuadrupolarFMode1", 0)
DEFINE_INSERT_FUNC(TidalQuadrupolarFMode2, REAL8, "TidalQuadrupolarFMode2", 0)
DEFINE_INSERT_FUNC(TidalOctupolarFMode1, REAL8, "TidalOctupolarFMode1", 0)
DEFINE_INSERT_FUNC(TidalOctupolarFMode2, REAL8, "TidalOctupolarFMode2", 0)
/* Note: some approximants like SEOBNRv2T/SEOBNRv4T will by default compute dQuadMon1, dQuadMon2 */
/* from TidalLambda1, TidalLambda2 using universal relations, rather than using the default value 0 */
DEFINE_INSERT_FUNC(dQuadMon1, REAL8, "dQuadMon1", 0)
DEFINE_INSERT_FUNC(dQuadMon2, REAL8, "dQuadMon2", 0)
DEFINE_INSERT_FUNC(Redshift, REAL8, "redshift", 0)
DEFINE_INSERT_FUNC(EccentricityFreq, REAL8, "f_ecc", LAL_DEFAULT_F_ECC)
DEFINE_INSERT_FUNC(Lscorr, INT4, "lscorr", 0)

DEFINE_INSERT_FUNC(NonGRPhi1, REAL8, "phi1", 0)
DEFINE_INSERT_FUNC(NonGRPhi2, REAL8, "phi2", 0)
DEFINE_INSERT_FUNC(NonGRPhi3, REAL8, "phi3", 0)
DEFINE_INSERT_FUNC(NonGRPhi4, REAL8, "phi4", 0)
DEFINE_INSERT_FUNC(NonGRDChi0, REAL8, "dchi0", 0)
DEFINE_INSERT_FUNC(NonGRDChi1, REAL8, "dchi1", 0)
DEFINE_INSERT_FUNC(NonGRDChi2, REAL8, "dchi2", 0)
DEFINE_INSERT_FUNC(NonGRDChi3, REAL8, "dchi3", 0)
DEFINE_INSERT_FUNC(NonGRDChi4, REAL8, "dchi4", 0)
DEFINE_INSERT_FUNC(NonGRDChi5, REAL8, "dchi5", 0)
DEFINE_INSERT_FUNC(NonGRDChi5L, REAL8, "dchi5l", 0)
DEFINE_INSERT_FUNC(NonGRDChi6, REAL8, "dchi6", 0)
DEFINE_INSERT_FUNC(NonGRDChi6L, REAL8, "dchi6l", 0)
DEFINE_INSERT_FUNC(NonGRDChi7, REAL8, "dchi7", 0)
DEFINE_INSERT_FUNC(NonGRDXi1, REAL8, "dxi1", 0)
DEFINE_INSERT_FUNC(NonGRDXi2, REAL8, "dxi2", 0)
DEFINE_INSERT_FUNC(NonGRDXi3, REAL8, "dxi3", 0)
DEFINE_INSERT_FUNC(NonGRDXi4, REAL8, "dxi4", 0)
DEFINE_INSERT_FUNC(NonGRDXi5, REAL8, "dxi5", 0)
DEFINE_INSERT_FUNC(NonGRDXi6, REAL8, "dxi6", 0)
DEFINE_INSERT_FUNC(NonGRDSigma1, REAL8, "dsigma1", 0)
DEFINE_INSERT_FUNC(NonGRDSigma2, REAL8, "dsigma2", 0)
DEFINE_INSERT_FUNC(NonGRDSigma3, REAL8, "dsigma3", 0)
DEFINE_INSERT_FUNC(NonGRDSigma4, REAL8, "dsigma4", 0)
DEFINE_INSERT_FUNC(NonGRDAlpha1, REAL8, "dalpha1", 0)
DEFINE_INSERT_FUNC(NonGRDAlpha2, REAL8, "dalpha2", 0)
DEFINE_INSERT_FUNC(NonGRDAlpha3, REAL8, "dalpha3", 0)
DEFINE_INSERT_FUNC(NonGRDAlpha4, REAL8, "dalpha4", 0)
DEFINE_INSERT_FUNC(NonGRDAlpha5, REAL8, "dalpha5", 0)
DEFINE_INSERT_FUNC(NonGRDBeta1, REAL8, "dbeta1", 0)
DEFINE_INSERT_FUNC(NonGRDBeta2, REAL8, "dbeta2", 0)
DEFINE_INSERT_FUNC(NonGRDBeta3, REAL8, "dbeta3", 0)
DEFINE_INSERT_FUNC(NonGRAlphaPPE, REAL8, "alphaPPE", 0)
DEFINE_INSERT_FUNC(NonGRBetaPPE, REAL8, "betaPPE", 0)
DEFINE_INSERT_FUNC(NonGRAlphaPPE0, REAL8, "alphaPPE0", 0)
DEFINE_INSERT_FUNC(NonGRBetaPPE0, REAL8, "betaPPE0", 0)
DEFINE_INSERT_FUNC(NonGRAlphaPPE1, REAL8, "alphaPPE1", 0)
DEFINE_INSERT_FUNC(NonGRBetaPPE1, REAL8, "betaPPE1", 0)
DEFINE_INSERT_FUNC(NonGRAlphaPPE2, REAL8, "alphaPPE2", 0)
DEFINE_INSERT_FUNC(NonGRBetaPPE2, REAL8, "betaPPE2", 0)
DEFINE_INSERT_FUNC(NonGRAlphaPPE3, REAL8, "alphaPPE3", 0)
DEFINE_INSERT_FUNC(NonGRBetaPPE3, REAL8, "betaPPE3", 0)
DEFINE_INSERT_FUNC(NonGRAlphaPPE4, REAL8, "alphaPPE4", 0)
DEFINE_INSERT_FUNC(NonGRBetaPPE4, REAL8, "betaPPE4", 0)
DEFINE_INSERT_FUNC(NonGRAlphaPPE5, REAL8, "alphaPPE5", 0)
DEFINE_INSERT_FUNC(NonGRBetaPPE5, REAL8, "betaPPE5", 0)
DEFINE_INSERT_FUNC(NonGRAlphaPPE6, REAL8, "alphaPPE6", 0)
DEFINE_INSERT_FUNC(NonGRBetaPPE6, REAL8, "betaPPE6", 0)
DEFINE_INSERT_FUNC(NonGRAlphaPPE7, REAL8, "alphaPPE7", 0)
DEFINE_INSERT_FUNC(NonGRBetaPPE7, REAL8, "betaPPE7", 0)
DEFINE_INSERT_FUNC(EnableLIV, INT4, "liv", 0)
DEFINE_INSERT_FUNC(NonGRLIVLogLambdaEff, REAL8, "log10lambda_eff", 100)
DEFINE_INSERT_FUNC(NonGRLIVASign, REAL8, "LIV_A_sign", 1)
DEFINE_INSERT_FUNC(NonGRLIVAlpha, REAL8, "nonGR_alpha", 0)

/* NLTides parameters */
/* used within LALSimInspiralTaylorF2NLTides.c */
DEFINE_INSERT_FUNC(NLTidesA1, REAL8, "nlTidesA1", 0)
DEFINE_INSERT_FUNC(NLTidesN1, REAL8, "nlTidesN1", 0)
DEFINE_INSERT_FUNC(NLTidesF1, REAL8, "nlTidesF1", 0)
DEFINE_INSERT_FUNC(NLTidesA2, REAL8, "nlTidesA2", 0)
DEFINE_INSERT_FUNC(NLTidesN2, REAL8, "nlTidesN2", 0)
DEFINE_INSERT_FUNC(NLTidesF2, REAL8, "nlTidesF2", 0)

/* SEOBNRv4P */
DEFINE_INSERT_FUNC(EOBChooseNumOrAnalHamDer, INT4, "EOBChooseNumOrAnalHamDer", 1)
DEFINE_INSERT_FUNC(EOBEllMaxForNyquistCheck, INT4, "EOBEllMaxForNyquistCheck", 5)


/* IMRPhenomX Parameters */
DEFINE_INSERT_FUNC(PhenomXInspiralPhaseVersion, INT4, "InsPhaseVersion", 104)
DEFINE_INSERT_FUNC(PhenomXInspiralAmpVersion, INT4, "InsAmpVersion", 103)
DEFINE_INSERT_FUNC(PhenomXIntermediatePhaseVersion, INT4, "IntPhaseVersion", 105)
DEFINE_INSERT_FUNC(PhenomXIntermediateAmpVersion, INT4, "IntAmpVersion", 104)
DEFINE_INSERT_FUNC(PhenomXRingdownPhaseVersion, INT4, "RDPhaseVersion", 105)
DEFINE_INSERT_FUNC(PhenomXRingdownAmpVersion, INT4, "RDAmpVersion", 103)
DEFINE_INSERT_FUNC(PhenomXPrecVersion, INT4, "PrecVersion", 223)
DEFINE_INSERT_FUNC(PhenomXPExpansionOrder, INT4, "ExpansionOrder", 5)
DEFINE_INSERT_FUNC(PhenomXPConvention, INT4, "Convention", 1)
DEFINE_INSERT_FUNC(PhenomXPFinalSpinMod, INT4, "FinalSpinMod", 3)
DEFINE_INSERT_FUNC(PhenomXPTransPrecessionMethod, INT4, "TransPrecessionMethod", 1)

/* IMRPhenomXHM Parameters */
DEFINE_INSERT_FUNC(PhenomXHMInspiralPhaseVersion, INT4, "InsPhaseHMVersion", 122019)
DEFINE_INSERT_FUNC(PhenomXHMIntermediatePhaseVersion, INT4, "IntPhaseHMVersion", 122019)
DEFINE_INSERT_FUNC(PhenomXHMRingdownPhaseVersion, INT4, "RDPhaseHMVersion", 122019)
DEFINE_INSERT_FUNC(PhenomXHMInspiralAmpVersion, INT4, "InsAmpHMVersion", 3)
DEFINE_INSERT_FUNC(PhenomXHMIntermediateAmpVersion, INT4, "IntAmpHMVersion", 2)
DEFINE_INSERT_FUNC(PhenomXHMRingdownAmpVersion, INT4, "RDAmpHMVersion", 0)
DEFINE_INSERT_FUNC(PhenomXHMInspiralAmpFitsVersion, INT4, "InsAmpFitsVersion", 122018)
DEFINE_INSERT_FUNC(PhenomXHMIntermediateAmpFitsVersion, INT4, "IntAmpFitsVersion", 122018)
DEFINE_INSERT_FUNC(PhenomXHMRingdownAmpFitsVersion, INT4, "RDAmpFitsVersion", 122018)
DEFINE_INSERT_FUNC(PhenomXHMPhaseRef21, REAL8, "PhaseRef21", 0.)
DEFINE_INSERT_FUNC(PhenomXHMThresholdMband, REAL8, "ThresholdMband", 0.001)
DEFINE_INSERT_FUNC(PhenomXHMAmpInterpolMB, INT4, "AmpInterpol", 1)

/* IMRPhenomXPHM Parameters */
DEFINE_INSERT_FUNC(PhenomXPHMMBandVersion, INT4, "MBandPrecVersion", 0)
DEFINE_INSERT_FUNC(PhenomXPHMThresholdMband, REAL8, "PrecThresholdMband", 0.001)
DEFINE_INSERT_FUNC(PhenomXPHMUseModes, INT4, "UseModes", 0)
DEFINE_INSERT_FUNC(PhenomXPHMModesL0Frame, INT4, "ModesL0Frame", 0)
DEFINE_INSERT_FUNC(PhenomXPHMPrecModes, INT4, "PrecModes", 0)
DEFINE_INSERT_FUNC(PhenomXPHMTwistPhenomHM, INT4, "TwistPhenomHM", 0)

/* LOOKUP FUNCTIONS */

DEFINE_LOOKUP_FUNC(ModesChoice, INT4, "modes", LAL_SIM_INSPIRAL_MODES_CHOICE_ALL)
DEFINE_LOOKUP_FUNC(FrameAxis, INT4, "axis", LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L)
DEFINE_LOOKUP_FUNC(Sideband, INT4, "sideband", 0)
DEFINE_LOOKUP_FUNC(NumRelData, String, "numreldata", NULL)

LALValue* XLALSimInspiralWaveformParamsLookupModeArray(LALDict *params)
{
	/* Initialise and set Default to NULL */
	LALValue * value = NULL;
	if (params && XLALDictContains(params, "ModeArray"))
	{
		LALDictEntry * entry = XLALDictLookup(params, "ModeArray");
		value = XLALValueDuplicate(XLALDictEntryGetValue(entry));
	}
	return value;
}

DEFINE_LOOKUP_FUNC(PNPhaseOrder, INT4, "phaseO", -1)
DEFINE_LOOKUP_FUNC(PNAmplitudeOrder, INT4, "ampO", -1)
DEFINE_LOOKUP_FUNC(PNEccentricityOrder, INT4, "eccO", -1)
DEFINE_LOOKUP_FUNC(PNSpinOrder, INT4, "spinO", -1)
DEFINE_LOOKUP_FUNC(PNTidalOrder, INT4, "tideO", -1)
DEFINE_LOOKUP_FUNC(GETides, INT4, "GEtideO", 0)
DEFINE_LOOKUP_FUNC(GMTides, INT4, "GMtideO", 0)

DEFINE_LOOKUP_FUNC(TidalLambda1, REAL8, "lambda1", 0)
DEFINE_LOOKUP_FUNC(TidalLambda2, REAL8, "lambda2", 0)
DEFINE_LOOKUP_FUNC(TidalOctupolarLambda1, REAL8, "TidalOctupolarLambda1", 0)
DEFINE_LOOKUP_FUNC(TidalOctupolarLambda2, REAL8, "TidalOctupolarLambda2", 0)
DEFINE_LOOKUP_FUNC(TidalHexadecapolarLambda1, REAL8, "TidalHexadecapolarLambda1", 0)
DEFINE_LOOKUP_FUNC(TidalHexadecapolarLambda2, REAL8, "TidalHexadecapolarLambda2", 0)
DEFINE_LOOKUP_FUNC(TidalQuadrupolarFMode1, REAL8, "TidalQuadrupolarFMode1", 0)
DEFINE_LOOKUP_FUNC(TidalQuadrupolarFMode2, REAL8, "TidalQuadrupolarFMode2", 0)
DEFINE_LOOKUP_FUNC(TidalOctupolarFMode1, REAL8, "TidalOctupolarFMode1", 0)
DEFINE_LOOKUP_FUNC(TidalOctupolarFMode2, REAL8, "TidalOctupolarFMode2", 0)
/* Note: some approximants like SEOBNRv2T/SEOBNRv4T will by default compute dQuadMon1, dQuadMon2 */
/* from TidalLambda1, TidalLambda2 using universal relations, rather than using the default value 0 */
DEFINE_LOOKUP_FUNC(dQuadMon1, REAL8, "dQuadMon1", 0)
DEFINE_LOOKUP_FUNC(dQuadMon2, REAL8, "dQuadMon2", 0)
DEFINE_LOOKUP_FUNC(Redshift, REAL8, "redshift", 0)
DEFINE_LOOKUP_FUNC(EccentricityFreq, REAL8, "f_ecc", LAL_DEFAULT_F_ECC)
DEFINE_LOOKUP_FUNC(Lscorr, INT4, "lscorr", 0)

DEFINE_LOOKUP_FUNC(NonGRPhi1, REAL8, "phi1", 0)
DEFINE_LOOKUP_FUNC(NonGRPhi2, REAL8, "phi2", 0)
DEFINE_LOOKUP_FUNC(NonGRPhi3, REAL8, "phi3", 0)
DEFINE_LOOKUP_FUNC(NonGRPhi4, REAL8, "phi4", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi0, REAL8, "dchi0", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi1, REAL8, "dchi1", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi2, REAL8, "dchi2", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi3, REAL8, "dchi3", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi4, REAL8, "dchi4", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi5, REAL8, "dchi5", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi5L, REAL8, "dchi5l", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi6, REAL8, "dchi6", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi6L, REAL8, "dchi6l", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi7, REAL8, "dchi7", 0)
DEFINE_LOOKUP_FUNC(NonGRDXi1, REAL8, "dxi1", 0)
DEFINE_LOOKUP_FUNC(NonGRDXi2, REAL8, "dxi2", 0)
DEFINE_LOOKUP_FUNC(NonGRDXi3, REAL8, "dxi3", 0)
DEFINE_LOOKUP_FUNC(NonGRDXi4, REAL8, "dxi4", 0)
DEFINE_LOOKUP_FUNC(NonGRDXi5, REAL8, "dxi5", 0)
DEFINE_LOOKUP_FUNC(NonGRDXi6, REAL8, "dxi6", 0)
DEFINE_LOOKUP_FUNC(NonGRDSigma1, REAL8, "dsigma1", 0)
DEFINE_LOOKUP_FUNC(NonGRDSigma2, REAL8, "dsigma2", 0)
DEFINE_LOOKUP_FUNC(NonGRDSigma3, REAL8, "dsigma3", 0)
DEFINE_LOOKUP_FUNC(NonGRDSigma4, REAL8, "dsigma4", 0)
DEFINE_LOOKUP_FUNC(NonGRDAlpha1, REAL8, "dalpha1", 0)
DEFINE_LOOKUP_FUNC(NonGRDAlpha2, REAL8, "dalpha2", 0)
DEFINE_LOOKUP_FUNC(NonGRDAlpha3, REAL8, "dalpha3", 0)
DEFINE_LOOKUP_FUNC(NonGRDAlpha4, REAL8, "dalpha4", 0)
DEFINE_LOOKUP_FUNC(NonGRDAlpha5, REAL8, "dalpha5", 0)
DEFINE_LOOKUP_FUNC(NonGRDBeta1, REAL8, "dbeta1", 0)
DEFINE_LOOKUP_FUNC(NonGRDBeta2, REAL8, "dbeta2", 0)
DEFINE_LOOKUP_FUNC(NonGRDBeta3, REAL8, "dbeta3", 0)
DEFINE_LOOKUP_FUNC(NonGRAlphaPPE, REAL8, "alphaPPE", 0)
DEFINE_LOOKUP_FUNC(NonGRBetaPPE, REAL8, "betaPPE", 0)
DEFINE_LOOKUP_FUNC(NonGRAlphaPPE0, REAL8, "alphaPPE0", 0)
DEFINE_LOOKUP_FUNC(NonGRBetaPPE0, REAL8, "betaPPE0", 0)
DEFINE_LOOKUP_FUNC(NonGRAlphaPPE1, REAL8, "alphaPPE1", 0)
DEFINE_LOOKUP_FUNC(NonGRBetaPPE1, REAL8, "betaPPE1", 0)
DEFINE_LOOKUP_FUNC(NonGRAlphaPPE2, REAL8, "alphaPPE2", 0)
DEFINE_LOOKUP_FUNC(NonGRBetaPPE2, REAL8, "betaPPE2", 0)
DEFINE_LOOKUP_FUNC(NonGRAlphaPPE3, REAL8, "alphaPPE3", 0)
DEFINE_LOOKUP_FUNC(NonGRBetaPPE3, REAL8, "betaPPE3", 0)
DEFINE_LOOKUP_FUNC(NonGRAlphaPPE4, REAL8, "alphaPPE4", 0)
DEFINE_LOOKUP_FUNC(NonGRBetaPPE4, REAL8, "betaPPE4", 0)
DEFINE_LOOKUP_FUNC(NonGRAlphaPPE5, REAL8, "alphaPPE5", 0)
DEFINE_LOOKUP_FUNC(NonGRBetaPPE5, REAL8, "betaPPE5", 0)
DEFINE_LOOKUP_FUNC(NonGRAlphaPPE6, REAL8, "alphaPPE6", 0)
DEFINE_LOOKUP_FUNC(NonGRBetaPPE6, REAL8, "betaPPE6", 0)
DEFINE_LOOKUP_FUNC(NonGRAlphaPPE7, REAL8, "alphaPPE7", 0)
DEFINE_LOOKUP_FUNC(NonGRBetaPPE7, REAL8, "betaPPE7", 0)
DEFINE_LOOKUP_FUNC(EnableLIV, INT4, "liv", 0)
DEFINE_LOOKUP_FUNC(NonGRLIVLogLambdaEff, REAL8, "log10lambda_eff", 100)
DEFINE_LOOKUP_FUNC(NonGRLIVASign, REAL8, "LIV_A_sign", 1)
DEFINE_LOOKUP_FUNC(NonGRLIVAlpha, REAL8, "nonGR_alpha", 0)

/* NLTides parameters */
/* used within LALSimInspiralTaylorF2NLTides.c */
DEFINE_LOOKUP_FUNC(NLTidesA1, REAL8, "nlTidesA1", 0)
DEFINE_LOOKUP_FUNC(NLTidesN1, REAL8, "nlTidesN1", 0)
DEFINE_LOOKUP_FUNC(NLTidesF1, REAL8, "nlTidesF1", 0)
DEFINE_LOOKUP_FUNC(NLTidesA2, REAL8, "nlTidesA2", 0)
DEFINE_LOOKUP_FUNC(NLTidesN2, REAL8, "nlTidesN2", 0)
DEFINE_LOOKUP_FUNC(NLTidesF2, REAL8, "nlTidesF2", 0)
/* SEOBNRv4P */
DEFINE_LOOKUP_FUNC(EOBChooseNumOrAnalHamDer, INT4, "EOBChooseNumOrAnalHamDer", 1)
DEFINE_LOOKUP_FUNC(EOBEllMaxForNyquistCheck, INT4, "EOBEllMaxForNyquistCheck", 5)

/* IMRPhenomX Parameters */
DEFINE_LOOKUP_FUNC(PhenomXInspiralPhaseVersion, INT4, "InsPhaseVersion", 104)
DEFINE_LOOKUP_FUNC(PhenomXInspiralAmpVersion, INT4, "InsAmpVersion", 103)
DEFINE_LOOKUP_FUNC(PhenomXIntermediatePhaseVersion, INT4, "IntPhaseVersion", 105)
DEFINE_LOOKUP_FUNC(PhenomXIntermediateAmpVersion, INT4, "IntAmpVersion", 104)
DEFINE_LOOKUP_FUNC(PhenomXRingdownPhaseVersion, INT4, "RDPhaseVersion", 105)
DEFINE_LOOKUP_FUNC(PhenomXRingdownAmpVersion, INT4, "RDAmpVersion", 103)
DEFINE_LOOKUP_FUNC(PhenomXPrecVersion, INT4, "PrecVersion", 223)
DEFINE_LOOKUP_FUNC(PhenomXPExpansionOrder, INT4, "ExpansionOrder", 5)
DEFINE_LOOKUP_FUNC(PhenomXPConvention, INT4, "Convention", 1)
DEFINE_LOOKUP_FUNC(PhenomXPFinalSpinMod, INT4, "FinalSpinMod", 3)
DEFINE_LOOKUP_FUNC(PhenomXPTransPrecessionMethod, INT4, "TransPrecessionMethod", 1)

/* IMRPhenomXHM Parameters */
DEFINE_LOOKUP_FUNC(PhenomXHMInspiralPhaseVersion, INT4, "InsPhaseHMVersion", 122019)
DEFINE_LOOKUP_FUNC(PhenomXHMIntermediatePhaseVersion, INT4, "IntPhaseHMVersion", 122019)
DEFINE_LOOKUP_FUNC(PhenomXHMRingdownPhaseVersion, INT4, "RDPhaseHMVersion", 122019)
DEFINE_LOOKUP_FUNC(PhenomXHMInspiralAmpVersion, INT4, "InsAmpHMVersion", 3)
DEFINE_LOOKUP_FUNC(PhenomXHMIntermediateAmpVersion, INT4, "IntAmpHMVersion", 2)
DEFINE_LOOKUP_FUNC(PhenomXHMRingdownAmpVersion, INT4, "RDAmpHMVersion", 0)
DEFINE_LOOKUP_FUNC(PhenomXHMInspiralAmpFitsVersion, INT4, "InsAmpFitsVersion", 122018)
DEFINE_LOOKUP_FUNC(PhenomXHMIntermediateAmpFitsVersion, INT4, "IntAmpFitsVersion", 122018)
DEFINE_LOOKUP_FUNC(PhenomXHMRingdownAmpFitsVersion, INT4, "RDAmpFitsVersion", 122018)
DEFINE_LOOKUP_FUNC(PhenomXHMPhaseRef21, REAL8, "PhaseRef21", 0.)
DEFINE_LOOKUP_FUNC(PhenomXHMThresholdMband, REAL8, "ThresholdMband", 0.001)
DEFINE_LOOKUP_FUNC(PhenomXHMAmpInterpolMB, INT4, "AmpInterpol", 1)

/* IMRPhenomXPHM */
DEFINE_LOOKUP_FUNC(PhenomXPHMMBandVersion, INT4, "MBandPrecVersion", 0)
DEFINE_LOOKUP_FUNC(PhenomXPHMThresholdMband, REAL8, "PrecThresholdMband", 0.001)
DEFINE_LOOKUP_FUNC(PhenomXPHMUseModes, INT4, "UseModes", 0)
DEFINE_LOOKUP_FUNC(PhenomXPHMModesL0Frame, INT4, "ModesL0Frame", 0)
DEFINE_LOOKUP_FUNC(PhenomXPHMPrecModes, INT4, "PrecModes", 0)
DEFINE_LOOKUP_FUNC(PhenomXPHMTwistPhenomHM, INT4, "TwistPhenomHM", 0)

/* ISDEFAULT FUNCTIONS */

DEFINE_ISDEFAULT_FUNC(ModesChoice, INT4, "modes", LAL_SIM_INSPIRAL_MODES_CHOICE_ALL)
DEFINE_ISDEFAULT_FUNC(FrameAxis, INT4, "axis", LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L)
DEFINE_ISDEFAULT_FUNC(Sideband, INT4, "sideband", 0)
DEFINE_ISDEFAULT_FUNC(NumRelData, String, "numreldata", NULL)

int XLALSimInspiralWaveformParamsModeArrayIsDefault(LALDict *params)
{
	return XLALSimInspiralWaveformParamsLookupModeArray(params) == NULL;
}

DEFINE_ISDEFAULT_FUNC(PNPhaseOrder, INT4, "phaseO", -1)
DEFINE_ISDEFAULT_FUNC(PNAmplitudeOrder, INT4, "ampO", -1)
DEFINE_ISDEFAULT_FUNC(PNEccentricityOrder, INT4, "eccO", -1)
DEFINE_ISDEFAULT_FUNC(PNSpinOrder, INT4, "spinO", -1)
DEFINE_ISDEFAULT_FUNC(PNTidalOrder, INT4, "tideO", -1)
DEFINE_ISDEFAULT_FUNC(GETides, INT4, "GEtideO", LAL_SIM_INSPIRAL_GETIDES_GSF23)
DEFINE_ISDEFAULT_FUNC(GMTides, INT4, "GMtideO", LAL_SIM_INSPIRAL_GMTIDES_PN)

DEFINE_ISDEFAULT_FUNC(TidalLambda1, REAL8, "lambda1", 0)
DEFINE_ISDEFAULT_FUNC(TidalLambda2, REAL8, "lambda2", 0)
DEFINE_ISDEFAULT_FUNC(TidalOctupolarLambda1, REAL8, "TidalOctupolarLambda1", 0)
DEFINE_ISDEFAULT_FUNC(TidalOctupolarLambda2, REAL8, "TidalOctupolarLambda2", 0)
DEFINE_ISDEFAULT_FUNC(TidalHexadecapolarLambda1, REAL8, "TidalHexadecapolarLambda1", 0)
DEFINE_ISDEFAULT_FUNC(TidalHexadecapolarLambda2, REAL8, "TidalHexadecapolarLambda2", 0)
DEFINE_ISDEFAULT_FUNC(TidalQuadrupolarFMode1, REAL8, "TidalQuadrupolarFMode1", 0)
DEFINE_ISDEFAULT_FUNC(TidalQuadrupolarFMode2, REAL8, "TidalQuadrupolarFMode2", 0)
DEFINE_ISDEFAULT_FUNC(TidalOctupolarFMode1, REAL8, "TidalOctupolarFMode1", 0)
DEFINE_ISDEFAULT_FUNC(TidalOctupolarFMode2, REAL8, "TidalOctupolarFMode2", 0)
/* Note: some approximants like SEOBNRv2T/SEOBNRv4T will by default compute dQuadMon1, dQuadMon2 */
/* from TidalLambda1, TidalLambda2 using universal relations, rather than using the default value 0 */
DEFINE_ISDEFAULT_FUNC(dQuadMon1, REAL8, "dQuadMon1", 0)
DEFINE_ISDEFAULT_FUNC(dQuadMon2, REAL8, "dQuadMon2", 0)
DEFINE_ISDEFAULT_FUNC(Redshift, REAL8, "redshift", 0)
DEFINE_ISDEFAULT_FUNC(EccentricityFreq, REAL8, "f_ecc", LAL_DEFAULT_F_ECC)

DEFINE_ISDEFAULT_FUNC(NonGRPhi1, REAL8, "phi1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRPhi2, REAL8, "phi2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRPhi3, REAL8, "phi3", 0)
DEFINE_ISDEFAULT_FUNC(NonGRPhi4, REAL8, "phi4", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi0, REAL8, "dchi0", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi1, REAL8, "dchi1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi2, REAL8, "dchi2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi3, REAL8, "dchi3", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi4, REAL8, "dchi4", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi5, REAL8, "dchi5", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi5L, REAL8, "dchi5l", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi6, REAL8, "dchi6", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi6L, REAL8, "dchi6l", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi7, REAL8, "dchi7", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDXi1, REAL8, "dxi1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDXi2, REAL8, "dxi2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDXi3, REAL8, "dxi3", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDXi4, REAL8, "dxi4", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDXi5, REAL8, "dxi5", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDXi6, REAL8, "dxi6", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDSigma1, REAL8, "dsigma1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDSigma2, REAL8, "dsigma2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDSigma3, REAL8, "dsigma3", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDSigma4, REAL8, "dsigma4", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDAlpha1, REAL8, "dalpha1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDAlpha2, REAL8, "dalpha2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDAlpha3, REAL8, "dalpha3", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDAlpha4, REAL8, "dalpha4", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDAlpha5, REAL8, "dalpha5", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDBeta1, REAL8, "dbeta1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDBeta2, REAL8, "dbeta2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDBeta3, REAL8, "dbeta3", 0)
DEFINE_ISDEFAULT_FUNC(NonGRAlphaPPE, REAL8, "alphaPPE", 0)
DEFINE_ISDEFAULT_FUNC(NonGRBetaPPE, REAL8, "betaPPE", 0)
DEFINE_ISDEFAULT_FUNC(NonGRAlphaPPE0, REAL8, "alphaPPE0", 0)
DEFINE_ISDEFAULT_FUNC(NonGRBetaPPE0, REAL8, "betaPPE0", 0)
DEFINE_ISDEFAULT_FUNC(NonGRAlphaPPE1, REAL8, "alphaPPE1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRBetaPPE1, REAL8, "betaPPE1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRAlphaPPE2, REAL8, "alphaPPE2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRBetaPPE2, REAL8, "betaPPE2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRAlphaPPE3, REAL8, "alphaPPE3", 0)
DEFINE_ISDEFAULT_FUNC(NonGRBetaPPE3, REAL8, "betaPPE3", 0)
DEFINE_ISDEFAULT_FUNC(NonGRAlphaPPE4, REAL8, "alphaPPE4", 0)
DEFINE_ISDEFAULT_FUNC(NonGRBetaPPE4, REAL8, "betaPPE4", 0)
DEFINE_ISDEFAULT_FUNC(NonGRAlphaPPE5, REAL8, "alphaPPE5", 0)
DEFINE_ISDEFAULT_FUNC(NonGRBetaPPE5, REAL8, "betaPPE5", 0)
DEFINE_ISDEFAULT_FUNC(NonGRAlphaPPE6, REAL8, "alphaPPE6", 0)
DEFINE_ISDEFAULT_FUNC(NonGRBetaPPE6, REAL8, "betaPPE6", 0)
DEFINE_ISDEFAULT_FUNC(NonGRAlphaPPE7, REAL8, "alphaPPE7", 0)
DEFINE_ISDEFAULT_FUNC(NonGRBetaPPE7, REAL8, "betaPPE7", 0)
DEFINE_ISDEFAULT_FUNC(EnableLIV, INT4, "liv", 0)
DEFINE_ISDEFAULT_FUNC(NonGRLIVLogLambdaEff, REAL8, "log10lambda_eff", 100)
DEFINE_ISDEFAULT_FUNC(NonGRLIVASign, REAL8, "LIV_A_sign", 1)
DEFINE_ISDEFAULT_FUNC(NonGRLIVAlpha, REAL8, "nonGR_alpha", 0)
/* SEOBNRv4P */
DEFINE_ISDEFAULT_FUNC(EOBChooseNumOrAnalHamDer, INT4, "EOBChooseNumOrAnalHamDer", 1)
DEFINE_ISDEFAULT_FUNC(EOBEllMaxForNyquistCheck, INT4, "EOBEllMaxForNyquistCheck", 5)

/* IMRPhenomX Parameters */
DEFINE_ISDEFAULT_FUNC(PhenomXInspiralPhaseVersion, INT4, "InsPhaseVersion", 104)
DEFINE_ISDEFAULT_FUNC(PhenomXInspiralAmpVersion, INT4, "InsAmpVersion", 103)
DEFINE_ISDEFAULT_FUNC(PhenomXIntermediatePhaseVersion, INT4, "IntPhaseVersion", 105)
DEFINE_ISDEFAULT_FUNC(PhenomXIntermediateAmpVersion, INT4, "IntAmpVersion", 104)
DEFINE_ISDEFAULT_FUNC(PhenomXRingdownPhaseVersion, INT4, "RDPhaseVersion", 105)
DEFINE_ISDEFAULT_FUNC(PhenomXRingdownAmpVersion, INT4, "RDAmpVersion", 103)
DEFINE_ISDEFAULT_FUNC(PhenomXPrecVersion, INT4, "PrecVersion", 223)
DEFINE_ISDEFAULT_FUNC(PhenomXPExpansionOrder, INT4, "ExpansionOrder", 5)
DEFINE_ISDEFAULT_FUNC(PhenomXPConvention, INT4, "Convention", 1)
DEFINE_ISDEFAULT_FUNC(PhenomXPFinalSpinMod, INT4, "FinalSpinMod", 3)
DEFINE_ISDEFAULT_FUNC(PhenomXPTransPrecessionMethod, INT4, "TransPrecessionMethod", 1)

/* IMRPhenomXHM Parameters */
DEFINE_ISDEFAULT_FUNC(PhenomXHMInspiralPhaseVersion, INT4, "InsPhaseHMVersion", 122019)
DEFINE_ISDEFAULT_FUNC(PhenomXHMIntermediatePhaseVersion, INT4, "IntPhaseHMVersion", 122019)
DEFINE_ISDEFAULT_FUNC(PhenomXHMRingdownPhaseVersion, INT4, "RDPhaseHMVersion", 122019)
DEFINE_ISDEFAULT_FUNC(PhenomXHMInspiralAmpVersion, INT4, "InsAmpHMVersion", 3)
DEFINE_ISDEFAULT_FUNC(PhenomXHMIntermediateAmpVersion, INT4, "IntAmpHMVersion", 2)
DEFINE_ISDEFAULT_FUNC(PhenomXHMRingdownAmpVersion, INT4, "RDAmpHMVersion", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXHMInspiralAmpFitsVersion, INT4, "InsAmpFitsVersion", 122018)
DEFINE_ISDEFAULT_FUNC(PhenomXHMIntermediateAmpFitsVersion, INT4, "IntAmpFitsVersion", 122018)
DEFINE_ISDEFAULT_FUNC(PhenomXHMRingdownAmpFitsVersion, INT4, "RDAmpFitsVersion", 122018)
DEFINE_ISDEFAULT_FUNC(PhenomXHMPhaseRef21, REAL8, "PhaseRef21", 0.)
DEFINE_ISDEFAULT_FUNC(PhenomXHMThresholdMband, REAL8, "ThresholdMband", 0.001)
DEFINE_ISDEFAULT_FUNC(PhenomXHMAmpInterpolMB, INT4, "AmpInterpol", 1)

/* IMRPhenomXPHM */
DEFINE_ISDEFAULT_FUNC(PhenomXPHMMBandVersion, INT4, "MBandPrecVersion", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMThresholdMband, REAL8, "PrecThresholdMband", 0.001)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMUseModes, INT4, "UseModes", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMModesL0Frame, INT4, "ModesL0Frame", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMPrecModes, INT4, "PrecModes", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMTwistPhenomHM, INT4, "TwistPhenomHM", 0)

#undef String
