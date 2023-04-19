#include <lal/LALStdio.h>
#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralWaveformFlags.h>
#include <lal/LALSimInspiralWaveformParams.h>
#include <math.h>
#include <gsl/gsl_poly.h>
#include "LALSimInspiralWaveformParams_common.c"

/* Warning message for unreviewed code.
   The XLAL warning messages are suppressed by default.
	 Here we temporarily change the lalDebugLevel to print the warning message and
	 inmediately afterwards we reset the lalDebugLevel to its original value.
	 In this way we avoid showing any other unwanted warning messages.
*/
#define UNREVIEWED_CODE_WARNING \
  int debug_level = XLALGetDebugLevel(); \
	XLALClobberDebugLevel(2); \
	XLAL_PRINT_WARNING("This code is currently UNREVIEWED, use with caution!"); \
	XLALClobberDebugLevel(debug_level);

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

/* New Waveforms Interface */

/** Check if the key belong to the know waveform params defined in LALSimInspiralWaveformParams_common.c */
int XLALSimInspiralCheckKnownREAL8Key(const char* key){
	size_t num_par = XLAL_NUM_ELEM(lalSimInspiralREAL8WaveformParams);
	for(size_t i = 0; i < num_par; i++) {
		if (strcmp(lalSimInspiralREAL8WaveformParams[i].name, key) == 0){
			return 1;
		}
 	}
	return 0;
}

/**
 * Check if the mass paramters inserted in the LALDict allow to determine the two mass components mass1, mass2.
 * It accepts only two mass parameters and at least one must be dimensionful.
 */
int XLALSimInspiralCheckDeterminationOfMasses(LALDict *params){
	UNREVIEWED_CODE_WARNING

	UINT2 dim_number =  0; /* dimensionful-mass counter */
	UINT2 nodim_number = 0; /* dimensionless-mass counter */
	UINT2 sym_number = 0; /* symmetric masses counter */
	const char *dimensionful_masses[6] = {"mass1", "mass2", "total_mass", \
	 "chirp_mass", "mass_difference", "reduced_mass"};
	const char *dimensionless_masses[2] = {"mass_ratio", "sym_mass_ratio"};
	const char *symetric_masses[6] = {"mass1", "mass2", "total_mass", "chirp_mass", "sym_mass_ratio", "reduced_mass"};

	for (size_t j = 0; j < sizeof(dimensionful_masses)/sizeof(*dimensionful_masses); ++j){
    if (XLALDictContains(params, dimensionful_masses[j]) == 1){
        dim_number += 1;
		}
	}
	for (size_t j = 0; j < sizeof(dimensionless_masses)/sizeof(*dimensionless_masses); ++j){
    if (XLALDictContains(params, dimensionless_masses[j]) == 1){
        nodim_number += 1;
		}
	}
	for (size_t j = 0; j < sizeof(symetric_masses)/sizeof(*symetric_masses); ++j){
    if (XLALDictContains(params, symetric_masses[j]) == 1){
        sym_number += 1;
		}
	}
	if(XLALDictContains(params, "mass1") && XLALDictContains(params, "mass2")){
		sym_number = 0;
	}

	if ((dim_number == 2 && nodim_number == 0) || (dim_number == 1 && nodim_number == 1)){
		if(sym_number == 2){
			XLAL_PRINT_WARNING("The larger object cannot be determined, assuming m1 >= m2.");
		}
		return XLAL_SUCCESS;
	}
	else if ((dim_number == 1 && nodim_number == 0) || dim_number == 0){
		XLAL_ERROR(XLAL_FAILURE, "Mass parameters are underspecified. Please include" \
		" one dimensionless and one dimensionful mass parameters, or two dimensionful masses.");
	}
	else{
		XLAL_ERROR(XLAL_FAILURE, "Mass parameters are overspecified. Please include" \
		" one dimensionless and one dimensionful mass parameters, or two dimensionful masses.");
	}
}


/* INSERT FUNCTIONS */
DEFINE_INSERT_FUNC(Mass1, REAL8, "mass1", 0)
DEFINE_INSERT_FUNC(Mass2, REAL8, "mass2", 0)
DEFINE_INSERT_FUNC(TotalMass, REAL8, "total_mass", 0)
DEFINE_INSERT_FUNC(MassRatio, REAL8, "mass_ratio", 0)
DEFINE_INSERT_FUNC(SymMassRatio, REAL8, "sym_mass_ratio", 0)
DEFINE_INSERT_FUNC(ChirpMass, REAL8, "chirp_mass", 0)
DEFINE_INSERT_FUNC(MassDifference, REAL8, "mass_difference", 0)
DEFINE_INSERT_FUNC(ReducedMass, REAL8, "reduced_mass", 0)
DEFINE_INSERT_FUNC(Spin1x, REAL8, "spin1x", 0)
DEFINE_INSERT_FUNC(Spin2x, REAL8, "spin2x", 0)
DEFINE_INSERT_FUNC(Spin1y, REAL8, "spin1y", 0)
DEFINE_INSERT_FUNC(Spin2y, REAL8, "spin2y", 0)
DEFINE_INSERT_FUNC(Spin1z, REAL8, "spin1z", 0)
DEFINE_INSERT_FUNC(Spin2z, REAL8, "spin2z", 0)
DEFINE_INSERT_FUNC(Spin1norm, REAL8, "spin1_norm", 0)
DEFINE_INSERT_FUNC(Spin2norm, REAL8, "spin2_norm", 0)
DEFINE_INSERT_FUNC(Spin1tilt, REAL8, "spin1_tilt", 0)
DEFINE_INSERT_FUNC(Spin2tilt, REAL8, "spin2_tilt", 0)
DEFINE_INSERT_FUNC(Spin1phi, REAL8, "spin1_phi", 0)
DEFINE_INSERT_FUNC(Spin2phi, REAL8, "spin2_phi", 0)
DEFINE_INSERT_FUNC(DeltaF, REAL8, "deltaF", 0)
DEFINE_INSERT_FUNC(DeltaT, REAL8, "deltaT", 0)
DEFINE_INSERT_FUNC(F22Ref, REAL8, "f22_ref", 0)
DEFINE_INSERT_FUNC(RefPhase, REAL8, "phi_ref", 0)
DEFINE_INSERT_FUNC(F22Start, REAL8, "f22_start", 0)
DEFINE_INSERT_FUNC(FMax, REAL8, "f_max", 0)
DEFINE_INSERT_FUNC(Distance, REAL8, "distance", 0)
DEFINE_INSERT_FUNC(Inclination, REAL8, "inclination", 0)
DEFINE_INSERT_FUNC(LongAscNodes, REAL8, "longAscNodes", 0)
DEFINE_INSERT_FUNC(Eccentricity, REAL8, "eccentricity", 0)
DEFINE_INSERT_FUNC(MeanPerAno, REAL8, "meanPerAno", 0)

DEFINE_INSERT_FUNC(Lmax, INT4, "lmax", 0)
DEFINE_INSERT_FUNC(ModesChoice, INT4, "modes", LAL_SIM_INSPIRAL_MODES_CHOICE_ALL)
DEFINE_INSERT_FUNC(FrameAxis, INT4, "axis", LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L)
DEFINE_INSERT_FUNC(Sideband, INT4, "sideband", 0)
DEFINE_INSERT_FUNC(NumRelData, String, "numreldata", NULL)

int XLALSimInspiralWaveformParamsInsertModeArray(LALDict *params, LALValue *value)
{
	return XLALDictInsertValue(params, "ModeArray", value);
}

int XLALSimInspiralWaveformParamsInsertModeArrayJframe(LALDict *params, LALValue *value)
{
	return XLALDictInsertValue(params, "ModeArrayJframe", value);
}

int XLALSimInspiralWaveformParamsInsertModeArrayFromModeString(LALDict *params, const char *modestr)
{
	UNREVIEWED_CODE_WARNING
        LALValue *modes = XLALSimInspiralModeArrayFromModeString(modestr);
        XLAL_CHECK(modes, XLAL_EFUNC);
	return XLALSimInspiralWaveformParamsInsertModeArray(params, modes);
}

int XLALSimInspiralWaveformParamsInsertModeArrayJframeFromModeString(LALDict *params, const char *modestr)
{
	UNREVIEWED_CODE_WARNING
        LALValue *modes = XLALSimInspiralModeArrayFromModeString(modestr);
        XLAL_CHECK(modes, XLAL_EFUNC);
	return XLALSimInspiralWaveformParamsInsertModeArrayJframe(params, modes);
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
DEFINE_INSERT_FUNC(FinalFreq, REAL8, "fend", 0)
DEFINE_INSERT_FUNC(OnlyFinal, INT4, "OnlyFinal", 0)

DEFINE_INSERT_FUNC(NonGRPhi1, REAL8, "phi1", 0)
DEFINE_INSERT_FUNC(NonGRPhi2, REAL8, "phi2", 0)
DEFINE_INSERT_FUNC(NonGRPhi3, REAL8, "phi3", 0)
DEFINE_INSERT_FUNC(NonGRPhi4, REAL8, "phi4", 0)
DEFINE_INSERT_FUNC(NonGRDChiMinus2, REAL8, "dchiMinus2", 0)
DEFINE_INSERT_FUNC(NonGRDChiMinus1, REAL8, "dchiMinus1", 0)
DEFINE_INSERT_FUNC(NonGRDChi0, REAL8, "dchi0", 0)
DEFINE_INSERT_FUNC(NonGRDChi1, REAL8, "dchi1", 0)
DEFINE_INSERT_FUNC(NonGRDChi2, REAL8, "dchi2", 0)
DEFINE_INSERT_FUNC(NonGRDChi3, REAL8, "dchi3", 0)
DEFINE_INSERT_FUNC(NonGRDChi3S, REAL8, "dchi3S", 0)
DEFINE_INSERT_FUNC(NonGRDChi3NS, REAL8, "dchi3NS", 0)
DEFINE_INSERT_FUNC(NonGRDChi4, REAL8, "dchi4", 0)
DEFINE_INSERT_FUNC(NonGRDChi4S, REAL8, "dchi4S", 0)
DEFINE_INSERT_FUNC(NonGRDChi4NS, REAL8, "dchi4NS", 0)
DEFINE_INSERT_FUNC(NonGRDChi5, REAL8, "dchi5", 0)
DEFINE_INSERT_FUNC(NonGRDChi5S, REAL8, "dchi5S", 0)
DEFINE_INSERT_FUNC(NonGRDChi5NS, REAL8, "dchi5NS", 0)
DEFINE_INSERT_FUNC(NonGRDChi5L, REAL8, "dchi5l", 0)
DEFINE_INSERT_FUNC(NonGRDChi5LS, REAL8, "dchi5lS", 0)
DEFINE_INSERT_FUNC(NonGRDChi5LNS, REAL8, "dchi5lNS", 0)
DEFINE_INSERT_FUNC(NonGRDChi6, REAL8, "dchi6", 0)
DEFINE_INSERT_FUNC(NonGRDChi6S, REAL8, "dchi6S", 0)
DEFINE_INSERT_FUNC(NonGRDChi6NS, REAL8, "dchi6NS", 0)
DEFINE_INSERT_FUNC(NonGRDChi6L, REAL8, "dchi6l", 0)
DEFINE_INSERT_FUNC(NonGRDChi7, REAL8, "dchi7", 0)
DEFINE_INSERT_FUNC(NonGRDChi7S, REAL8, "dchi7S", 0)
DEFINE_INSERT_FUNC(NonGRDChi7NS, REAL8, "dchi7NS", 0)
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
DEFINE_INSERT_FUNC(NonGRDChikappaS, REAL8, "dchikappaS", 0)
DEFINE_INSERT_FUNC(NonGRDChikappaA, REAL8, "dchikappaA", 0)

/* NLTides parameters */
/* used within LALSimInspiralTaylorF2NLTides.c */
DEFINE_INSERT_FUNC(NLTidesA1, REAL8, "nlTidesA1", 0)
DEFINE_INSERT_FUNC(NLTidesN1, REAL8, "nlTidesN1", 0)
DEFINE_INSERT_FUNC(NLTidesF1, REAL8, "nlTidesF1", 0)
DEFINE_INSERT_FUNC(NLTidesA2, REAL8, "nlTidesA2", 0)
DEFINE_INSERT_FUNC(NLTidesN2, REAL8, "nlTidesN2", 0)
DEFINE_INSERT_FUNC(NLTidesF2, REAL8, "nlTidesF2", 0)
DEFINE_INSERT_FUNC(DOmega220, REAL8, "domega220", 0)
DEFINE_INSERT_FUNC(DTau220, REAL8, "dtau220", 0)
DEFINE_INSERT_FUNC(DOmega210, REAL8, "domega210", 0)
DEFINE_INSERT_FUNC(DTau210, REAL8, "dtau210", 0)
DEFINE_INSERT_FUNC(DOmega330, REAL8, "domega330", 0)
DEFINE_INSERT_FUNC(DTau330, REAL8, "dtau330", 0)
DEFINE_INSERT_FUNC(DOmega440, REAL8, "domega440", 0)
DEFINE_INSERT_FUNC(DTau440, REAL8, "dtau440", 0)
DEFINE_INSERT_FUNC(DOmega550, REAL8, "domega550", 0)
DEFINE_INSERT_FUNC(DTau550, REAL8, "dtau550", 0)

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
DEFINE_INSERT_FUNC(PhenomXPrecVersion, INT4, "PrecVersion", 300)
DEFINE_INSERT_FUNC(PhenomXReturnCoPrec, INT4, "ReturnCoPrec", 0)
DEFINE_INSERT_FUNC(PhenomXPExpansionOrder, INT4, "ExpansionOrder", 5)
DEFINE_INSERT_FUNC(PhenomXPConvention, INT4, "Convention", 1)
DEFINE_INSERT_FUNC(PhenomXPFinalSpinMod, INT4, "FinalSpinMod", 4)
DEFINE_INSERT_FUNC(PhenomXPTransPrecessionMethod, INT4, "TransPrecessionMethod", 1)
DEFINE_INSERT_FUNC(PhenomXPSpinTaylorVersion, String, "SpinTaylorVersion", NULL)
DEFINE_INSERT_FUNC(PhenomXPSpinTaylorCoarseFactor, INT4, "SpinTaylorCoarseFactor",10);

/* IMRPhenomXAS_NRTidalvX */
DEFINE_INSERT_FUNC(PhenomXTidalFlag, INT4, "PhenXTidal", 0)

/* IMRPhenomXHM Parameters */
DEFINE_INSERT_FUNC(PhenomXHMReleaseVersion, INT4, "PhenomXHMReleaseVersion", 122022)
DEFINE_INSERT_FUNC(PhenomXHMInspiralPhaseVersion, INT4, "InsPhaseHMVersion", 122019)
DEFINE_INSERT_FUNC(PhenomXHMIntermediatePhaseVersion, INT4, "IntPhaseHMVersion", 122019)
DEFINE_INSERT_FUNC(PhenomXHMRingdownPhaseVersion, INT4, "RDPhaseHMVersion", 122019)
DEFINE_INSERT_FUNC(PhenomXHMInspiralAmpVersion, INT4, "InsAmpHMVersion", 3)
DEFINE_INSERT_FUNC(PhenomXHMIntermediateAmpVersion, INT4, "IntAmpHMVersion", 2)
DEFINE_INSERT_FUNC(PhenomXHMRingdownAmpVersion, INT4, "RDAmpHMVersion", 0)
DEFINE_INSERT_FUNC(PhenomXHMInspiralAmpFitsVersion, INT4, "InsAmpFitsVersion", 122018)
DEFINE_INSERT_FUNC(PhenomXHMIntermediateAmpFitsVersion, INT4, "IntAmpFitsVersion", 122018)
DEFINE_INSERT_FUNC(PhenomXHMRingdownAmpFitsVersion, INT4, "RDAmpFitsVersion", 122018)
DEFINE_INSERT_FUNC(PhenomXHMInspiralAmpFreqsVersion, INT4, "InsAmpFreqsVersion", 122018)
DEFINE_INSERT_FUNC(PhenomXHMIntermediateAmpFreqsVersion, INT4, "IntAmpFreqsVersion", 122018)
DEFINE_INSERT_FUNC(PhenomXHMRingdownAmpFreqsVersion, INT4, "RDAmpFreqsVersion", 122018)
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

/* IMRPhenomTHM Parameters */
DEFINE_INSERT_FUNC(PhenomTHMInspiralVersion, INT4, "InspiralVersion", 0)
DEFINE_INSERT_FUNC(PhenomTPHMMergerVersion, INT4, "MergerVersion", 1)

/* IMRPhenomX_PNR Parameters */
DEFINE_INSERT_FUNC(PhenomXPNRUseTunedAngles, INT4, "PNRUseTunedAngles", 0)
DEFINE_INSERT_FUNC(PhenomXPNRUseTunedCoprec, INT4, "PNRUseTunedCoprec", 0)
DEFINE_INSERT_FUNC(PhenomXPNRUseTunedCoprec33, INT4, "PNRUseTunedCoprec33", 0)
// Option to only be used when actively tuning PNR Coprec relative to XHM wherein the non-precessing final spin is used
DEFINE_INSERT_FUNC(PhenomXPNRUseInputCoprecDeviations, INT4, "PNRUseInputCoprecDeviations", 0)
// Dev option for forcing 22 phase derivative inspiral values to align with XHM at a low ref frequency
DEFINE_INSERT_FUNC(PhenomXPNRForceXHMAlignment, INT4, "PNRForceXHMAlignment", 0)
/* Toggle output of XAS phase for debugging purposes */
DEFINE_INSERT_FUNC(PhenomXOnlyReturnPhase, INT4, "PhenomXOnlyReturnPhase", 0)
DEFINE_INSERT_FUNC(PhenomXPNRInterpTolerance, REAL8, "PNRInterpTolerance", 0.01)

/* IMRPhenomX_PNR_Asymmetry Parameters */
DEFINE_INSERT_FUNC(PhenomXAntisymmetricWaveform, INT4, "AntisymmetricWaveform", 0)


/* IMRPhenomXCP Parameters */
DEFINE_INSERT_FUNC(PhenomXCPMU1, REAL8, "MU1", 0)
DEFINE_INSERT_FUNC(PhenomXCPMU2, REAL8, "MU2", 0)
DEFINE_INSERT_FUNC(PhenomXCPMU3, REAL8, "MU3", 0)
DEFINE_INSERT_FUNC(PhenomXCPMU4, REAL8, "MU4", 0)
DEFINE_INSERT_FUNC(PhenomXCPNU0, REAL8, "NU0", 0)
DEFINE_INSERT_FUNC(PhenomXCPNU4, REAL8, "NU4", 0)
DEFINE_INSERT_FUNC(PhenomXCPNU5, REAL8, "NU5", 0)
DEFINE_INSERT_FUNC(PhenomXCPNU6, REAL8, "NU6", 0)
DEFINE_INSERT_FUNC(PhenomXCPZETA1, REAL8, "ZETA1", 0)
DEFINE_INSERT_FUNC(PhenomXCPZETA2, REAL8, "ZETA2", 0)
/* l=3, m=3 */
DEFINE_INSERT_FUNC(PhenomXCPMU1l3m3, REAL8, "MU1l3m3", 0)
DEFINE_INSERT_FUNC(PhenomXCPMU2l3m3, REAL8, "MU2l3m3", 0)
DEFINE_INSERT_FUNC(PhenomXCPMU3l3m3, REAL8, "MU3l3m3", 0)
DEFINE_INSERT_FUNC(PhenomXCPMU4l3m3, REAL8, "MU4l3m3", 0)
DEFINE_INSERT_FUNC(PhenomXCPNU0l3m3, REAL8, "NU0l3m3", 0)
DEFINE_INSERT_FUNC(PhenomXCPNU4l3m3, REAL8, "NU4l3m3", 0)
DEFINE_INSERT_FUNC(PhenomXCPNU5l3m3, REAL8, "NU5l3m3", 0)
DEFINE_INSERT_FUNC(PhenomXCPNU6l3m3, REAL8, "NU6l3m3", 0)
DEFINE_INSERT_FUNC(PhenomXCPZETA1l3m3, REAL8, "ZETA1l3m3", 0)
DEFINE_INSERT_FUNC(PhenomXCPZETA2l3m3, REAL8, "ZETA2l3m3", 0)


/* FLEXIBLE INPUT PARAMETERS FUNCTIONS */

/* Auxiliar mass arguments transformation functions */

REAL8 XLALSimInspiralGetMassRatioFromSymMassRatio(REAL8 sym_mass_ratio){
	UNREVIEWED_CODE_WARNING

		REAL8 mass_ratio;
		REAL8 x;

		if (sym_mass_ratio <= 0.25){
			x = 1.0 - 4.0 * sym_mass_ratio;
			mass_ratio = 0.5 * (1. - pow(x, 0.5)) / sym_mass_ratio - 1.0;
		}
		else{
				XLAL_ERROR(XLAL_EINVAL, "Invalid value of symmetric mass ratio given");
		}
		return mass_ratio;
}

REAL8 XLALSimInspiralGetMassRatioFromChirpMassComponentMass2(REAL8 chirp_mass, REAL8 component_mass){
	UNREVIEWED_CODE_WARNING
	return 1./XLALSimInspiralGetMassRatioFromChirpMassComponentMass1(chirp_mass, component_mass);
}

REAL8 XLALSimInspiralGetMassRatioFromChirpMassComponentMass1(REAL8 chirp_mass, REAL8 component_mass){
		UNREVIEWED_CODE_WARNING
		REAL8 c;
		REAL8 x;
		REAL8 mass_ratio;

		c = pow((chirp_mass / component_mass), 5);
		x = 1.5 * pow((3.0/c), 0.5);  // Implement trigonometric and hyperbolic solutions of cubic equation
		if (x < 1.0){
			mass_ratio = 3.0 * cos(acos(x) / 3.0) / x;
		}
		else{
			mass_ratio = 3.0 * cosh(acosh(x) / 3.0) / x;
		}
		return mass_ratio;
}

/* MASS parameters LOOKUP functions */

/**
 * Compute mass1 from any possible combination of 2 mass parameters inserted in the LALDict.
 * If the combination does not allow to distinguish the largest object then it assumes m1 > m2.
 * mass_ratio is defined as q = m2/m1 and mass_difference as m1 - m2.
 */
REAL8 XLALSimInspiralWaveformParamsLookupMass1(LALDict *params){

	REAL8 mass1 = 0;
	REAL8 mass2;
	REAL8 total_mass;
	REAL8 mass_ratio;
	REAL8 sym_mass_ratio;
	REAL8 mass_difference;
	REAL8 chirp_mass;
	REAL8 reduced_mass;
	REAL8 x;


		if (XLALDictContains(params, "mass1") == 1){
			mass1 = XLALDictLookupREAL8Value(params, "mass1");
			XLAL_CHECK(mass1 > 0, XLAL_EDOM, "mass1 must be positive");
			return mass1;
		}

		UNREVIEWED_CODE_WARNING

		INT4 status = XLALSimInspiralCheckDeterminationOfMasses(params);
		XLAL_CHECK(status==XLAL_SUCCESS, status, "Mass1 cannot be determined");

		if (XLALDictContains(params, "mass2") == 1){
			mass2 = XLALDictLookupREAL8Value(params, "mass2");
			XLAL_CHECK(mass2 > 0, XLAL_EDOM, "mass2 must be positive");
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				mass1 = mass2 / mass_ratio;
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				mass1 = mass2 / mass_ratio;
			}
			else if (XLALDictContains(params, "mass_difference") == 1){
				mass_difference = XLALSimInspiralWaveformParamsLookupMassDifference(params);
				mass1 = mass2 + mass_difference;
			}
			else if (XLALDictContains(params, "total_mass") == 1){
				total_mass = XLALSimInspiralWaveformParamsLookupTotalMass(params);
				mass1 = total_mass - mass2;
			}
			else if (XLALDictContains(params, "reduced_mass") == 1){
				reduced_mass = XLALSimInspiralWaveformParamsLookupReducedMass(params);
				mass1 = reduced_mass * mass2 / (mass2 - reduced_mass);
			}
			else if (XLALDictContains(params, "chirp_mass") == 1){
				chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromChirpMassComponentMass2(chirp_mass, mass2);
				mass1 = mass2 / mass_ratio;
			}
		}
		else if (XLALDictContains(params, "total_mass") == 1){
			total_mass = XLALSimInspiralWaveformParamsLookupTotalMass(params);
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				mass1 = total_mass / (1. + mass_ratio);
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				mass1 = total_mass / (1. + mass_ratio);
			}
			else if (XLALDictContains(params, "mass_difference") == 1){
				mass_difference = XLALSimInspiralWaveformParamsLookupMassDifference(params);
				mass1 = 0.5 * (total_mass + mass_difference);
			}
			else if (XLALDictContains(params, "reduced_mass") == 1){
				reduced_mass = XLALSimInspiralWaveformParamsLookupReducedMass(params);
				if (total_mass < 4.0 * reduced_mass){
					XLAL_ERROR(XLAL_EINVAL, "Invalid combination of total mass and reduced mass given");
				}
				x = total_mass * (total_mass - 4.0 * reduced_mass);
				mass_difference = pow(x, 0.5);
				mass1 = total_mass - 0.5 * (total_mass - mass_difference);
			}
			else if (XLALDictContains(params, "chirp_mass") == 1){
				chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
				sym_mass_ratio = pow((chirp_mass / total_mass), 5.0/3.0);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				mass1 = total_mass / (1.0 + mass_ratio);
			}
		}
		else if (XLALDictContains(params, "reduced_mass") == 1){
			reduced_mass = XLALSimInspiralWaveformParamsLookupReducedMass(params);
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				mass1 = (1. + mass_ratio) * reduced_mass / mass_ratio;
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				mass1 = (1. + mass_ratio) * reduced_mass / mass_ratio;
			}
			else if (XLALDictContains(params, "mass_difference") == 1){
				mass_difference = XLALSimInspiralWaveformParamsLookupMassDifference(params);
				x = 4.0 * pow(reduced_mass,2) + pow(mass_difference,2);
				mass1 = reduced_mass + 0.5 * mass_difference + 0.5 * pow(x, 0.5);
			}
			else if (XLALDictContains(params, "chirp_mass") == 1){
				chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
				total_mass = pow(pow(chirp_mass, 5) / pow(reduced_mass, 3), 0.5);
				x = total_mass * (total_mass - 4.0 * reduced_mass);
				if (x >= 0){
					mass_difference = pow(x, 0.5);
					mass1 = total_mass - 0.5 * (total_mass - mass_difference);
				}
				else {
					XLAL_ERROR(XLAL_FAILURE, "Invalid combination of reduced mass and chirp mass given");
				}
			}
		}
		else if (XLALDictContains(params, "mass_difference") == 1){
			mass_difference = XLALSimInspiralWaveformParamsLookupMassDifference(params);
			XLAL_CHECK(fabs(mass_difference) > 0, XLAL_EDOM, "Mass difference cannot be zero if it is combined with a dimensionless mass parameter.");
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				mass1 = mass_difference / (1.0 - mass_ratio);
				XLAL_CHECK(mass1 > 0, XLAL_EDOM, "Inconsistent values of mass_difference and mass_ratio.");
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				if (mass_difference < 0){
					mass_ratio = 1./mass_ratio;
				}
				mass1 = mass_difference / (1.0 - mass_ratio);
			}
			else if (XLALDictContains(params, "chirp_mass") == 1){
				chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
				/* We will solve the equation
				*	 m2^6 + 3 mdiff m2^5 + 3 mdiff^2 m2^4 + mdiff^3 m2^3 - 2 chirpmass^5 m2 - mdiff chripmass^5 = 0
				*  which is obtained by substituting m1 = mdiff + m2 into the chirpmass formula.
				*  There will be 6 possible solutions that can be complex. We need to pick those that are real and positive, and then the one that also gives m1 positive.
				*/
				double chirp_mass_5 = pow(chirp_mass, 5);
				double coefficients[7] = { -mass_difference * chirp_mass_5, -2 * chirp_mass_5, 0, pow(mass_difference, 3), 3 * pow(mass_difference, 2), 3 * mass_difference, 1 };
				double complex m2[6];

				gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(7);
				gsl_poly_complex_solve(coefficients, 7, w, (double *)m2);
				gsl_poly_complex_workspace_free(w);

				// Pick real and positive solution
				for (UINT2 i = 0; i < 6; i++)
		    {
					if(cimag(m2[i]) == 0 && creal(m2[i])>0){
							mass1 = creal(m2[i]) + mass_difference;
							if (mass1 > 0) break;
					}
		    }
				XLAL_CHECK(mass1>0, XLAL_FAILURE, "Could not find solution for mass1\n");
			}
		}
		else if (XLALDictContains(params, "chirp_mass") == 1){
			chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				sym_mass_ratio = mass_ratio / pow((1.0 + mass_ratio), 2);
				total_mass = chirp_mass / pow(sym_mass_ratio, 3.0/5.0);
				mass1 = total_mass / (1.0 + mass_ratio);
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				total_mass = chirp_mass / pow(sym_mass_ratio, 3.0/5.0);
				mass1 = total_mass / (1.0 + mass_ratio);
			}
		}

		return mass1;

}

/**
 * Compute mass2 from any possible combination of 2 mass parameters inserted in the LALDict.
 * If the combination does not allow to distinguish the largest object then it assumes m1 > m2.
 * mass_ratio is defined as q = m2/m1 and mass_difference as m1 - m2.
 */
REAL8 XLALSimInspiralWaveformParamsLookupMass2(LALDict *params){

	REAL8 mass2 = 0;
	REAL8 mass1;
	REAL8 total_mass;
	REAL8 mass_ratio;
	REAL8 sym_mass_ratio;
	REAL8 mass_difference;
	REAL8 chirp_mass;
	REAL8 reduced_mass;
	REAL8 x;

		if (XLALDictContains(params, "mass2") == 1){
			mass2 = XLALDictLookupREAL8Value(params, "mass2");
      return mass2;
		}

		UNREVIEWED_CODE_WARNING

		INT4 status = XLALSimInspiralCheckDeterminationOfMasses(params);
		XLAL_CHECK(status==XLAL_SUCCESS, status, "Mass2 cannot be determined");

		if (XLALDictContains(params, "mass1") == 1){
			mass1 = XLALDictLookupREAL8Value(params, "mass1");
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				mass2 = mass1 * mass_ratio;
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				mass2 = mass1 * mass_ratio;
			}
			else if (XLALDictContains(params, "mass_difference") == 1){
				mass_difference = XLALSimInspiralWaveformParamsLookupMassDifference(params);
				mass2 = mass1 - mass_difference;
			}
			else if (XLALDictContains(params, "total_mass") == 1){
				total_mass = XLALSimInspiralWaveformParamsLookupTotalMass(params);
				mass2 = total_mass - mass1;
			}
			else if (XLALDictContains(params, "reduced_mass") == 1){
				reduced_mass = XLALSimInspiralWaveformParamsLookupReducedMass(params);
				mass2 = reduced_mass * mass1 / (mass1 - reduced_mass);
			}
			else if (XLALDictContains(params, "chirp_mass") == 1){
				chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromChirpMassComponentMass1(chirp_mass, mass1);
				mass2 = mass1 * mass_ratio;
			}
		}

		else if (XLALDictContains(params, "total_mass") == 1){
			total_mass = XLALSimInspiralWaveformParamsLookupTotalMass(params);
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				mass2 = total_mass / (1. + mass_ratio) * mass_ratio;
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				mass2 = total_mass / (1. + mass_ratio) * mass_ratio;
			}
			else if (XLALDictContains(params, "mass_difference") == 1){
				mass_difference = XLALSimInspiralWaveformParamsLookupMassDifference(params);
				mass2 = 0.5 * (total_mass - mass_difference);
			}
			else if (XLALDictContains(params, "reduced_mass") == 1){
				reduced_mass = XLALSimInspiralWaveformParamsLookupReducedMass(params);
				if (total_mass < 4.0 * reduced_mass){
					XLAL_ERROR(XLAL_EINVAL, "Invalid combination of total mass and reduced mass given");
				}
				x = total_mass * (total_mass - 4.0 * reduced_mass);
				mass_difference = pow(x, 0.5);
				mass2 = 0.5 * (total_mass - mass_difference);
			}
			else if (XLALDictContains(params, "chirp_mass") == 1){
				chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
				sym_mass_ratio = pow((chirp_mass / total_mass), 5.0/3.0);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				mass2 = total_mass / (1.0 + mass_ratio) * mass_ratio;
			}
		}

		else if (XLALDictContains(params, "reduced_mass") == 1){
			reduced_mass = XLALSimInspiralWaveformParamsLookupReducedMass(params);
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				mass2 = (1. + mass_ratio) * reduced_mass;
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				mass2 = (1. + mass_ratio) * reduced_mass;
			}
			else if (XLALDictContains(params, "mass_difference") == 1){
				mass_difference = XLALSimInspiralWaveformParamsLookupMassDifference(params);
				x = 4.0 * pow(reduced_mass,2) + pow(mass_difference,2);
				mass2 = reduced_mass + 0.5 * mass_difference + 0.5 * pow(x, 0.5) - mass_difference;
			}
			else if (XLALDictContains(params, "chirp_mass") == 1){
				chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
				total_mass = pow(pow(chirp_mass, 5) / pow(reduced_mass, 3), 0.5);
				x = total_mass * (total_mass - 4.0 * reduced_mass);
				if (x >= 0){
					mass_difference = pow(x, 0.5);
					mass2 = 0.5 * (total_mass - mass_difference);
				}
				else {
					XLAL_ERROR(XLAL_FAILURE, "Invalid combination of reduced mass and chirp mass given");
				}
			}
		}

		else if (XLALDictContains(params, "mass_difference") == 1){
			mass_difference = XLALSimInspiralWaveformParamsLookupMassDifference(params);
			XLAL_CHECK(fabs(mass_difference) > 0, XLAL_EDOM, "Mass difference cannot be zero if it is combined with a dimensionless mass parameter.");
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				mass2 = mass_difference / (1.0 - mass_ratio) * mass_ratio;
				XLAL_CHECK(mass2 > 0, XLAL_EDOM, "Inconsistent values of mass_difference and mass_ratio.");
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				if (mass_difference < 0){
					mass_ratio = 1./mass_ratio;
				}
				mass2 = mass_ratio * mass_difference / (1.0 - mass_ratio);
			}
			else if (XLALDictContains(params, "chirp_mass") == 1){
				chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
				/* We will solve the equation
				*	 m2^6 + 3 mdiff m2^5 + 3 mdiff^2 m2^4 + mdiff^3 m2^3 - 2 chirpmass^5 m2 - mdiff chripmass^5 = 0
				*  which is obtained by substituting m1 = mdiff + m2 into the chirpmass formula.
				*  There will be 6 possible solutions that can be complex. We need to pick those that are real and positive, and then the one that also gives m1 positive.
				*/
				double chirp_mass_5 = pow(chirp_mass, 5);
				double coefficients[7] = { -mass_difference * chirp_mass_5, -2 * chirp_mass_5, 0, pow(mass_difference, 3), 3 * pow(mass_difference, 2), 3 * mass_difference, 1 };
				double complex m2[6];

				gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(7);
				gsl_poly_complex_solve(coefficients, 7, w, (double *)m2);
				gsl_poly_complex_workspace_free(w);

				// Pick real and positive solution
				for (UINT2 i = 0; i < 6; i++)
		    {
					if(cimag(m2[i]) == 0 && creal(m2[i])>0){
							mass1 = creal(m2[i]) + mass_difference;
							if (mass1 > 0){
								mass2 = creal(m2[i]);
								break;
							}
					}
		    }
				XLAL_CHECK(mass2 > 0, XLAL_FAILURE, "Could not find solution for mass2\n");
			}
		}
		else if (XLALDictContains(params, "chirp_mass") == 1){
			chirp_mass = XLALSimInspiralWaveformParamsLookupChirpMass(params);
			if (XLALDictContains(params, "mass_ratio") == 1){
				mass_ratio = XLALSimInspiralWaveformParamsLookupMassRatio(params);
				sym_mass_ratio = mass_ratio / pow((1.0 + mass_ratio), 2);
				total_mass = chirp_mass / pow(sym_mass_ratio, 3.0/5.0);
				mass2 = total_mass / (1.0 + mass_ratio) * mass_ratio;
			}
			else if (XLALDictContains(params, "sym_mass_ratio") == 1){
				sym_mass_ratio = XLALSimInspiralWaveformParamsLookupSymMassRatio(params);
				mass_ratio = XLALSimInspiralGetMassRatioFromSymMassRatio(sym_mass_ratio);
				total_mass = chirp_mass / pow(sym_mass_ratio, 3.0/5.0);
				mass2 = total_mass / (1.0 + mass_ratio) * mass_ratio;
			}
		}

		return mass2;

}


REAL8 XLALSimInspiralWaveformParamsLookupTotalMass(LALDict *params){
	UNREVIEWED_CODE_WARNING
	REAL8 mass1;
	REAL8 mass2;
	REAL8 total_mass;

	if (XLALDictContains(params, "total_mass") == 1){
		total_mass = XLALDictLookupREAL8Value(params, "total_mass");
		XLAL_CHECK(total_mass > 0, XLAL_EDOM, "total_mass must be positive");
	}
	else {
	mass1 = XLALSimInspiralWaveformParamsLookupMass1(params);
	mass2 = XLALSimInspiralWaveformParamsLookupMass2(params);
	total_mass = mass1 + mass2;
	}
	return total_mass;
}

REAL8 XLALSimInspiralWaveformParamsLookupMassRatio(LALDict *params){
	UNREVIEWED_CODE_WARNING
	REAL8 mass1;
	REAL8 mass2;
	REAL8 mass_ratio;

	if (XLALDictContains(params, "mass_ratio") == 1){
		mass_ratio = XLALDictLookupREAL8Value(params, "mass_ratio");
		XLAL_CHECK(mass_ratio > 0, XLAL_EDOM, "mass_ratio must be positive");
	}
	else {
	mass1 = XLALSimInspiralWaveformParamsLookupMass1(params);
	mass2 = XLALSimInspiralWaveformParamsLookupMass2(params);
	mass_ratio = mass2 / mass1;
	}
	return mass_ratio;
}

REAL8 XLALSimInspiralWaveformParamsLookupSymMassRatio(LALDict *params){
  UNREVIEWED_CODE_WARNING
	REAL8 mass1;
	REAL8 mass2;
	REAL8 sym_mass_ratio;

	if (XLALDictContains(params, "sym_mass_ratio") == 1){
		sym_mass_ratio = XLALDictLookupREAL8Value(params, "sym_mass_ratio");
		XLAL_CHECK(sym_mass_ratio > 0 && sym_mass_ratio <= 0.25, XLAL_EDOM, "sym_mass_ratio must be between (0, 0.25]");
	}
	else {
	mass1 = XLALSimInspiralWaveformParamsLookupMass1(params);
	mass2 = XLALSimInspiralWaveformParamsLookupMass2(params);
	sym_mass_ratio = mass1 * mass2 / pow(mass1 + mass2, 2);
	}
	return sym_mass_ratio;
}

REAL8 XLALSimInspiralWaveformParamsLookupChirpMass(LALDict *params){
	UNREVIEWED_CODE_WARNING
	REAL8 mass1;
	REAL8 mass2;
	REAL8 chirp_mass;

	if (XLALDictContains(params, "chirp_mass") == 1){
		chirp_mass = XLALDictLookupREAL8Value(params, "chirp_mass");
		XLAL_CHECK(chirp_mass > 0, XLAL_EDOM, "chirp_mass must be positive");
	}
	else {
	mass1 = XLALSimInspiralWaveformParamsLookupMass1(params);
	mass2 = XLALSimInspiralWaveformParamsLookupMass2(params);
	chirp_mass = pow(mass1 * mass2, 3.0 / 5.0) / pow(mass1 + mass2, 1.0 / 5.0);
	}
	return chirp_mass;
}

REAL8 XLALSimInspiralWaveformParamsLookupMassDifference(LALDict *params){
	UNREVIEWED_CODE_WARNING
	REAL8 mass1;
	REAL8 mass2;
	REAL8 mass_difference;

	if (XLALDictContains(params, "mass_difference") == 1){
		mass_difference = XLALDictLookupREAL8Value(params, "mass_difference");
	}
	else {
	mass1 = XLALSimInspiralWaveformParamsLookupMass1(params);
	mass2 = XLALSimInspiralWaveformParamsLookupMass2(params);
	mass_difference = mass1 - mass2;
	}
	return mass_difference;
}

REAL8 XLALSimInspiralWaveformParamsLookupReducedMass(LALDict *params){
	UNREVIEWED_CODE_WARNING
	REAL8 mass1;
	REAL8 mass2;
	REAL8 reduced_mass;

	if (XLALDictContains(params, "reduced_mass") == 1){
		reduced_mass = XLALDictLookupREAL8Value(params, "reduced_mass");
		XLAL_CHECK(reduced_mass > 0, XLAL_EDOM, "reduced_mass must be positive");
	}
	else {
	mass1 = XLALSimInspiralWaveformParamsLookupMass1(params);
	mass2 = XLALSimInspiralWaveformParamsLookupMass2(params);
	reduced_mass = mass1 * mass2 / (mass1 + mass2);
	}
	return reduced_mass;
}

/*  Auxiliar spin transformations functions
 *  In polar components, the tilt angle is the angles between the Z axis and the spin vector.
 *  The phi angle is the angle between the X axis and the projection of the spin vector in the X-Y plane.
 */

REAL8 XLALSimInspiralGetCartesianSpinXFromPolar(REAL8 spin_norm, REAL8 spin_tilt, REAL8 spin_phi){
		UNREVIEWED_CODE_WARNING
		REAL8 spinx;
		spinx = spin_norm * sin(spin_tilt) * cos(spin_phi);
		return spinx;
}


REAL8 XLALSimInspiralGetCartesianSpinYFromPolar(REAL8 spin_norm, REAL8 spin_tilt, REAL8 spin_phi){
		UNREVIEWED_CODE_WARNING
		REAL8 spiny;
		spiny = spin_norm * sin(spin_tilt) * sin(spin_phi);
		return spiny;
}

REAL8 XLALSimInspiralGetCartesianSpinZFromPolar(REAL8 spin_norm, REAL8 spin_tilt){
		UNREVIEWED_CODE_WARNING
		REAL8 spinz;
		spinz = spin_norm * cos(spin_tilt);
		return spinz;
}

REAL8 XLALSimInspiralGetPolarSpin_normFromCartesian(REAL8 spinx, REAL8 spiny, REAL8 spinz){
		UNREVIEWED_CODE_WARNING
		REAL8 spin_norm;
		spin_norm = sqrt(spinx* spinx+ spiny* spiny+ spinz * spinz);
		return spin_norm;
}

REAL8 XLALSimInspiralGetPolarSpin_tiltFromCartesian(REAL8 spinx, REAL8 spiny, REAL8 spinz){
		UNREVIEWED_CODE_WARNING
		REAL8 spin_tilt;
		spin_tilt = acos(spinz / sqrt(spinx* spinx+ spiny* spiny+ spinz * spinz));
		return spin_tilt;
}

REAL8 XLALSimInspiralGetPolarSpin_phiFromCartesian(REAL8 spiny, REAL8 spinz){
		UNREVIEWED_CODE_WARNING
		REAL8 spin_phi;
		spin_phi = atan(spiny / spinz);
		return spin_phi;
}

/* SPIN parameters LOOKUP functions 
 *  Read cartesian and polar spin components from LALDict
 *  These functions do not check if the spins are overdetermined.
 *  For one spin you cannot mix cartesian and polar components. The three components must be or cartesian or polar.
 */

REAL8 XLALSimInspiralWaveformParamsLookupSpin1x(LALDict *params){
				REAL8 spin1x = -1;
        REAL8 spin1_norm;
        REAL8 spin1_tilt;
        REAL8 spin1_phi;
				if (XLALDictContains(params, "spin1x") == 1){
					spin1x = XLALDictLookupREAL8Value(params, "spin1x");
					return spin1x;
				}
				UNREVIEWED_CODE_WARNING
				if ((XLALDictContains(params, "spin1_norm") == 1) && (XLALDictContains(params, "spin1_tilt") == 1) && (XLALDictContains(params, "spin1_phi") == 1)){
    			spin1_norm = XLALDictLookupREAL8Value(params,"spin1_norm");
    			spin1_tilt =  XLALDictLookupREAL8Value(params,"spin1_tilt");
          spin1_phi = XLALDictLookupREAL8Value(params,"spin1_phi");
          spin1x = XLALSimInspiralGetCartesianSpinXFromPolar(spin1_norm, spin1_tilt, spin1_phi);
				}
        else {
					spin1x = 0;
					XLAL_PRINT_WARNING("Not enough information provided to determine spin1x. Assuming zero as a default value. \n");
    		}
        return spin1x;
	}

REAL8 XLALSimInspiralWaveformParamsLookupSpin1y(LALDict *params){
        REAL8 spin1y = -1;
        REAL8 spin1_norm;
        REAL8 spin1_tilt;
        REAL8 spin1_phi;
        if (XLALDictContains(params, "spin1y") == 1){
          spin1y = XLALDictLookupREAL8Value(params, "spin1y");
					return spin1y;
        }
				UNREVIEWED_CODE_WARNING
        if  ((XLALDictContains(params, "spin1_norm") == 1) && (XLALDictContains(params, "spin1_tilt") == 1) \
				&& (XLALDictContains(params, "spin1_phi") == 1)){
          spin1_norm = XLALDictLookupREAL8Value(params,"spin1_norm");
          spin1_tilt =  XLALDictLookupREAL8Value(params,"spin1_tilt");
          spin1_phi = XLALDictLookupREAL8Value(params,"spin1_phi");
          spin1y = XLALSimInspiralGetCartesianSpinYFromPolar(spin1_norm, spin1_tilt, spin1_phi);}
        else {
					spin1y = 0;
					XLAL_PRINT_WARNING("Not enough information provided to determine spin1y. Assuming zero as a default value. \n");
    		}
        return spin1y;
}

REAL8 XLALSimInspiralWaveformParamsLookupSpin1z(LALDict *params){
        REAL8 spin1z = -1;
        REAL8 spin1_norm;
        REAL8 spin1_tilt;
        if (XLALDictContains(params, "spin1z") == 1){
          spin1z = XLALDictLookupREAL8Value(params, "spin1z");
					return spin1z;
        }
				UNREVIEWED_CODE_WARNING
				if  ((XLALDictContains(params, "spin1_norm") == 1) && (XLALDictContains(params, "spin1_tilt") == 1)){
          spin1_norm = XLALDictLookupREAL8Value(params,"spin1_norm");
          spin1_tilt = XLALDictLookupREAL8Value(params,"spin1_tilt");
          spin1z = XLALSimInspiralGetCartesianSpinZFromPolar(spin1_norm, spin1_tilt);}
        else {
					spin1z = 0;
					XLAL_PRINT_WARNING("Not enough information provided to determine spin1z. Assuming zero as a default value. \n");
    		}
        return spin1z;
      }

REAL8 XLALSimInspiralWaveformParamsLookupSpin2x(LALDict *params){
        REAL8 spin2x = -1;
        REAL8 spin2_norm;
        REAL8 spin2_tilt;
        REAL8 spin2_phi;
        if (XLALDictContains(params, "spin2x") == 1){
          spin2x = XLALDictLookupREAL8Value(params, "spin2x");
					return spin2x;
        }
				UNREVIEWED_CODE_WARNING
				if ((XLALDictContains(params, "spin2_norm") == 1) && (XLALDictContains(params, "spin2_tilt") == 1) \
				&& (XLALDictContains(params, "spin2_phi") == 1)){
          spin2_norm = XLALDictLookupREAL8Value(params,"spin2_norm");
          spin2_tilt =  XLALDictLookupREAL8Value(params,"spin2_tilt");
          spin2_phi = XLALDictLookupREAL8Value(params,"spin2_phi");
          spin2x = XLALSimInspiralGetCartesianSpinXFromPolar(spin2_norm, spin2_tilt, spin2_phi);}
        else {
					spin2x = 0;
					XLAL_PRINT_WARNING("Not enough information provided to determine spin2x. Assuming zero as a default value. \n");
    		}
        return spin2x;
      }

REAL8 XLALSimInspiralWaveformParamsLookupSpin2y(LALDict *params){
        REAL8 spin2y = -1;
        REAL8 spin2_norm;
        REAL8 spin2_tilt;
        REAL8 spin2_phi;
        if (XLALDictContains(params, "spin2y") == 1){
                spin2y = XLALDictLookupREAL8Value(params, "spin2y");
								return spin2y;
        }
				UNREVIEWED_CODE_WARNING
        if ((XLALDictContains(params, "spin2_norm") == 1) && (XLALDictContains(params, "spin2_tilt") == 1) \
				 && (XLALDictContains(params, "spin2_phi") == 1)){
          spin2_norm = XLALDictLookupREAL8Value(params,"spin2_norm");
          spin2_tilt =  XLALDictLookupREAL8Value(params,"spin2_tilt");
          spin2_phi = XLALDictLookupREAL8Value(params,"spin2_phi");
          spin2y = XLALSimInspiralGetCartesianSpinYFromPolar(spin2_norm, spin2_tilt, spin2_phi);}
        else {
					spin2y = 0;
					XLAL_PRINT_WARNING("Not enough information provided to determine spin2y. Assuming zero as a default value. \n");
    		}
        return spin2y;
      }

REAL8 XLALSimInspiralWaveformParamsLookupSpin2z(LALDict *params){
        REAL8 spin2z = -1;
        REAL8 spin2_norm;
        REAL8 spin2_tilt;
        if (XLALDictContains(params, "spin2z") == 1){
          spin2z = XLALDictLookupREAL8Value(params, "spin2z");
					return spin2z;
        }
				UNREVIEWED_CODE_WARNING
        if ((XLALDictContains(params, "spin2_norm") == 1) && (XLALDictContains(params, "spin2_tilt") == 1)){
          spin2_norm = XLALDictLookupREAL8Value(params,"spin2_norm");
          spin2_tilt = XLALDictLookupREAL8Value(params,"spin2_tilt");
          spin2z = XLALSimInspiralGetCartesianSpinZFromPolar(spin2_norm, spin2_tilt);}
        else {
					spin2z = 0;
					XLAL_PRINT_WARNING("Not enough information provided to determine spin2z. Assuming zero as a default value. \n");
    		}
        return spin2z;
      }


REAL8 XLALSimInspiralWaveformParamsLookupSpin1norm(LALDict *params){
	UNREVIEWED_CODE_WARNING
        REAL8 spin1_norm = -1;
        REAL8 spin1x;
        REAL8 spin1y;
        REAL8 spin1z;
        if (XLALDictContains(params, "spin1_norm") == 1){
                spin1_norm = XLALDictLookupREAL8Value(params, "spin1_norm");
        }
        else if  ((XLALDictContains(params, "spin1x") == 1) && (XLALDictContains(params, "spin1y") == 1) && (XLALDictContains(params, "spin1z") == 1))
               {
               spin1x = XLALDictLookupREAL8Value(params,"spin1x");
               spin1y =  XLALDictLookupREAL8Value(params,"spin1y");
               spin1z = XLALDictLookupREAL8Value(params,"spin1z");
               spin1_norm = XLALSimInspiralGetPolarSpin_normFromCartesian(spin1x, spin1y, spin1z);}
        else {
        XLAL_ERROR_REAL8(XLAL_FAILURE, "Not enough information provided for spin1_norm calculation\n");
    }
        return spin1_norm;
        }

REAL8 XLALSimInspiralWaveformParamsLookupSpin1tilt(LALDict *params){
	UNREVIEWED_CODE_WARNING
        REAL8 spin1_tilt = -1;
        REAL8 spin1x;
        REAL8 spin1y;
        REAL8 spin1z;
        if (XLALDictContains(params, "spin1_tilt") == 1){
                spin1_tilt = XLALDictLookupREAL8Value(params, "spin1_tilt");
        }
        else if  ((XLALDictContains(params, "spin1x") == 1) && (XLALDictContains(params, "spin1y") == 1) && (XLALDictContains(params, "spin1z") == 1))
               {
               spin1x = XLALDictLookupREAL8Value(params,"spin1x");
               spin1y =  XLALDictLookupREAL8Value(params,"spin1y");
               spin1z = XLALDictLookupREAL8Value(params,"spin1z");
               spin1_tilt = XLALSimInspiralGetPolarSpin_tiltFromCartesian(spin1x, spin1y, spin1z);}
        else {
        XLAL_ERROR_REAL8(XLAL_FAILURE, "Not enough information provided for spin1_tilt calculation\n");
    }
        return spin1_tilt;
        }


REAL8 XLALSimInspiralWaveformParamsLookupSpin1phi(LALDict *params){
	UNREVIEWED_CODE_WARNING
        REAL8 spin1_phi = -1;
        REAL8 spin1y;
        REAL8 spin1z;
        if (XLALDictContains(params, "spin1_phi") == 1){
                spin1_phi = XLALDictLookupREAL8Value(params, "spin1_phi");
        }
        else if  ((XLALDictContains(params, "spin1y") == 1) && (XLALDictContains(params, "spin1z") == 1))
               {
               spin1y =  XLALDictLookupREAL8Value(params,"spin1y");
               spin1z = XLALDictLookupREAL8Value(params,"spin1z");
               spin1_phi = XLALSimInspiralGetPolarSpin_phiFromCartesian(spin1y, spin1z);}
        else {
        XLAL_ERROR_REAL8(XLAL_FAILURE, "Not enough information provided for spin1_phi calculation\n");
    }
        return spin1_phi;
        }


REAL8 XLALSimInspiralWaveformParamsLookupSpin2norm(LALDict *params){
	UNREVIEWED_CODE_WARNING
        REAL8 spin2_norm = -1;
        REAL8 spin2x;
        REAL8 spin2y;
        REAL8 spin2z;
        if (XLALDictContains(params, "spin2_norm") == 1){
                spin2_norm = XLALDictLookupREAL8Value(params, "spin2_norm");
        }
        else if  ((XLALDictContains(params, "spin2x") == 1) && (XLALDictContains(params, "spin2y") == 1) && (XLALDictContains(params, "spin2z") == 1))
               {
               spin2x = XLALDictLookupREAL8Value(params,"spin2x");
               spin2y =  XLALDictLookupREAL8Value(params,"spin2y");
               spin2z = XLALDictLookupREAL8Value(params,"spin2z");
               spin2_norm = XLALSimInspiralGetPolarSpin_normFromCartesian(spin2x, spin2y, spin2z);}
        else {
        XLAL_ERROR_REAL8(XLAL_FAILURE, "Not enough information provided for spin2_norm calculation\n");
    }
        return spin2_norm;
        }

REAL8 XLALSimInspiralWaveformParamsLookupSpin2tilt(LALDict *params){
	UNREVIEWED_CODE_WARNING
        REAL8 spin2_tilt = -1;
        REAL8 spin2x;
        REAL8 spin2y;
        REAL8 spin2z;
        if (XLALDictContains(params, "spin2_tilt") == 1){
                spin2_tilt = XLALDictLookupREAL8Value(params, "spin2_tilt");
        }
        else if  ((XLALDictContains(params, "spin2x") == 1) && (XLALDictContains(params, "spin2y") == 1) && (XLALDictContains(params, "spin2z") == 1))
               {
               spin2x = XLALDictLookupREAL8Value(params,"spin2x");
               spin2y =  XLALDictLookupREAL8Value(params,"spin2y");
               spin2z = XLALDictLookupREAL8Value(params,"spin2z");
               spin2_tilt = XLALSimInspiralGetPolarSpin_tiltFromCartesian(spin2x, spin2y, spin2z);}
        else {
        XLAL_ERROR_REAL8(XLAL_FAILURE, "Not enough information provided for spin2_tilt calculation\n");
    }
        return spin2_tilt;
        }

REAL8 XLALSimInspiralWaveformParamsLookupSpin2phi(LALDict *params){
	UNREVIEWED_CODE_WARNING
        REAL8 spin2_phi = -1;
        REAL8 spin2y;
        REAL8 spin2z;
        if (XLALDictContains(params, "spin2_phi") == 1){
                spin2_phi = XLALDictLookupREAL8Value(params, "spin2_phi");
        }
        else if  ((XLALDictContains(params, "spin2y") == 1) && (XLALDictContains(params, "spin2z") == 1))
               {
               spin2y =  XLALDictLookupREAL8Value(params,"spin2y");
               spin2z = XLALDictLookupREAL8Value(params,"spin2z");
               spin2_phi = XLALSimInspiralGetPolarSpin_phiFromCartesian(spin2y, spin2z);}
        else {
        XLAL_ERROR_REAL8(XLAL_FAILURE, "Not enough information provided for spin2_phi calculation\n");
    }
        return spin2_phi;
}

DEFINE_LOOKUP_FUNC(DeltaF, REAL8, "deltaF", 0)
DEFINE_LOOKUP_FUNC(DeltaT, REAL8, "deltaT", 0)
DEFINE_LOOKUP_FUNC(F22Ref, REAL8, "f22_ref", 0)
DEFINE_LOOKUP_FUNC(RefPhase, REAL8, "phi_ref", 0)
DEFINE_LOOKUP_FUNC(F22Start, REAL8, "f22_start", 0)
DEFINE_LOOKUP_FUNC(FMax, REAL8, "f_max", 0)
DEFINE_LOOKUP_FUNC(Distance, REAL8, "distance", 0)
DEFINE_LOOKUP_FUNC(Inclination, REAL8, "inclination", 0)
DEFINE_LOOKUP_FUNC(LongAscNodes, REAL8, "longAscNodes", 0)
DEFINE_LOOKUP_FUNC(Eccentricity, REAL8, "eccentricity", 0)
DEFINE_LOOKUP_FUNC(MeanPerAno, REAL8, "meanPerAno", 0)

DEFINE_LOOKUP_FUNC(Lmax, INT4, "lmax", 0)
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

LALValue* XLALSimInspiralWaveformParamsLookupModeArrayJframe(LALDict *params)
{
	/* Initialise and set Default to NULL */
	LALValue * value = NULL;
	if (params && XLALDictContains(params, "ModeArrayJframe"))
	{
		LALDictEntry * entry = XLALDictLookup(params, "ModeArrayJframe");
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
DEFINE_LOOKUP_FUNC(FinalFreq, REAL8, "fend", 0)
DEFINE_LOOKUP_FUNC(OnlyFinal, INT4, "OnlyFinal", 0)

DEFINE_LOOKUP_FUNC(NonGRPhi1, REAL8, "phi1", 0)
DEFINE_LOOKUP_FUNC(NonGRPhi2, REAL8, "phi2", 0)
DEFINE_LOOKUP_FUNC(NonGRPhi3, REAL8, "phi3", 0)
DEFINE_LOOKUP_FUNC(NonGRPhi4, REAL8, "phi4", 0)
DEFINE_LOOKUP_FUNC(NonGRDChiMinus2, REAL8, "dchiMinus2", 0)
DEFINE_LOOKUP_FUNC(NonGRDChiMinus1, REAL8, "dchiMinus1", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi0, REAL8, "dchi0", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi1, REAL8, "dchi1", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi2, REAL8, "dchi2", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi3, REAL8, "dchi3", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi3S, REAL8, "dchi3S", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi3NS, REAL8, "dchi3NS", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi4, REAL8, "dchi4", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi4S, REAL8, "dchi4S", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi4NS, REAL8, "dchi4NS", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi5, REAL8, "dchi5", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi5S, REAL8, "dchi5S", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi5NS, REAL8, "dchi5NS", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi5L, REAL8, "dchi5l", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi5LS, REAL8, "dchi5lS", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi5LNS, REAL8, "dchi5lNS", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi6, REAL8, "dchi6", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi6S, REAL8, "dchi6S", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi6NS, REAL8, "dchi6NS", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi6L, REAL8, "dchi6l", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi7, REAL8, "dchi7", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi7S, REAL8, "dchi7S", 0)
DEFINE_LOOKUP_FUNC(NonGRDChi7NS, REAL8, "dchi7NS", 0)
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
DEFINE_LOOKUP_FUNC(NonGRDChikappaS, REAL8, "dchikappaS", 0)
DEFINE_LOOKUP_FUNC(NonGRDChikappaA, REAL8, "dchikappaA", 0)

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
DEFINE_LOOKUP_FUNC(PhenomXPrecVersion, INT4, "PrecVersion", 300)
DEFINE_LOOKUP_FUNC(PhenomXReturnCoPrec, INT4, "ReturnCoPrec", 0)
DEFINE_LOOKUP_FUNC(PhenomXPExpansionOrder, INT4, "ExpansionOrder", 5)
DEFINE_LOOKUP_FUNC(PhenomXPConvention, INT4, "Convention", 1)
DEFINE_LOOKUP_FUNC(PhenomXPFinalSpinMod, INT4, "FinalSpinMod", 4)
DEFINE_LOOKUP_FUNC(PhenomXPTransPrecessionMethod, INT4, "TransPrecessionMethod", 1)
DEFINE_LOOKUP_FUNC(PhenomXPSpinTaylorVersion, String, "SpinTaylorVersion", NULL)
DEFINE_LOOKUP_FUNC(PhenomXPSpinTaylorCoarseFactor, INT4, "SpinTaylorCoarseFactor",10);

/* IMRPhenomX_NRTidalvX Parameters */
DEFINE_LOOKUP_FUNC(PhenomXTidalFlag, INT4, "PhenXTidal", 0)

/* IMRPhenomXHM Parameters */
DEFINE_LOOKUP_FUNC(PhenomXHMReleaseVersion, INT4, "PhenomXHMReleaseVersion", 122022)
DEFINE_LOOKUP_FUNC(PhenomXHMInspiralPhaseVersion, INT4, "InsPhaseHMVersion", 122019)
DEFINE_LOOKUP_FUNC(PhenomXHMIntermediatePhaseVersion, INT4, "IntPhaseHMVersion", 122019)
DEFINE_LOOKUP_FUNC(PhenomXHMRingdownPhaseVersion, INT4, "RDPhaseHMVersion", 122019)
DEFINE_LOOKUP_FUNC(PhenomXHMInspiralAmpVersion, INT4, "InsAmpHMVersion", 3)
DEFINE_LOOKUP_FUNC(PhenomXHMIntermediateAmpVersion, INT4, "IntAmpHMVersion", 2)
DEFINE_LOOKUP_FUNC(PhenomXHMRingdownAmpVersion, INT4, "RDAmpHMVersion", 0)
DEFINE_LOOKUP_FUNC(PhenomXHMInspiralAmpFitsVersion, INT4, "InsAmpFitsVersion", 122018)
DEFINE_LOOKUP_FUNC(PhenomXHMIntermediateAmpFitsVersion, INT4, "IntAmpFitsVersion", 122018)
DEFINE_LOOKUP_FUNC(PhenomXHMRingdownAmpFitsVersion, INT4, "RDAmpFitsVersion", 122018)
DEFINE_LOOKUP_FUNC(PhenomXHMInspiralAmpFreqsVersion, INT4, "InsAmpFreqsVersion", 122018)
DEFINE_LOOKUP_FUNC(PhenomXHMIntermediateAmpFreqsVersion, INT4, "IntAmpFreqsVersion", 122018)
DEFINE_LOOKUP_FUNC(PhenomXHMRingdownAmpFreqsVersion, INT4, "RDAmpFreqsVersion", 122018)
DEFINE_LOOKUP_FUNC(PhenomXHMPhaseRef21, REAL8, "PhaseRef21", 0.)
DEFINE_LOOKUP_FUNC(PhenomXHMThresholdMband, REAL8, "ThresholdMband", 0.001)
DEFINE_LOOKUP_FUNC(PhenomXHMAmpInterpolMB, INT4, "AmpInterpol", 1)
DEFINE_LOOKUP_FUNC(DOmega220, REAL8, "domega220", 0)
DEFINE_LOOKUP_FUNC(DTau220, REAL8, "dtau220", 0)
DEFINE_LOOKUP_FUNC(DOmega210, REAL8, "domega210", 0)
DEFINE_LOOKUP_FUNC(DTau210, REAL8, "dtau210", 0)
DEFINE_LOOKUP_FUNC(DOmega330, REAL8, "domega330", 0)
DEFINE_LOOKUP_FUNC(DTau330, REAL8, "dtau330", 0)
DEFINE_LOOKUP_FUNC(DOmega440, REAL8, "domega440", 0)
DEFINE_LOOKUP_FUNC(DTau440, REAL8, "dtau440", 0)
DEFINE_LOOKUP_FUNC(DOmega550, REAL8, "domega550", 0)
DEFINE_LOOKUP_FUNC(DTau550, REAL8, "dtau550", 0)

/* IMRPhenomXPHM */
DEFINE_LOOKUP_FUNC(PhenomXPHMMBandVersion, INT4, "MBandPrecVersion", 0)
DEFINE_LOOKUP_FUNC(PhenomXPHMThresholdMband, REAL8, "PrecThresholdMband", 0.001)
DEFINE_LOOKUP_FUNC(PhenomXPHMUseModes, INT4, "UseModes", 0)
DEFINE_LOOKUP_FUNC(PhenomXPHMModesL0Frame, INT4, "ModesL0Frame", 0)
DEFINE_LOOKUP_FUNC(PhenomXPHMPrecModes, INT4, "PrecModes", 0)
DEFINE_LOOKUP_FUNC(PhenomXPHMTwistPhenomHM, INT4, "TwistPhenomHM", 0)

/* IMRPhenomTHM Parameters */
DEFINE_LOOKUP_FUNC(PhenomTHMInspiralVersion, INT4, "InspiralVersion", 0)
DEFINE_LOOKUP_FUNC(PhenomTPHMMergerVersion, INT4, "MergerVersion", 1)

/* IMRPhenomX_PNR Parameters */
DEFINE_LOOKUP_FUNC(PhenomXPNRUseTunedAngles, INT4, "PNRUseTunedAngles", 0)
DEFINE_LOOKUP_FUNC(PhenomXPNRUseTunedCoprec, INT4, "PNRUseTunedCoprec", 0)
DEFINE_LOOKUP_FUNC(PhenomXPNRUseTunedCoprec33, INT4, "PNRUseTunedCoprec33", 0)
// Option to only be used when actively tuning PNR Coprec relative to XHM wherein the non-precessing final spin is used
DEFINE_LOOKUP_FUNC(PhenomXPNRUseInputCoprecDeviations, INT4, "PNRUseInputCoprecDeviations", 0)
// Dev option for forcing 22 phase derivative inspiral values to align with XHM at a low ref frequency
DEFINE_LOOKUP_FUNC(PhenomXPNRForceXHMAlignment, INT4, "PNRForceXHMAlignment", 0)
/* Toggle output of XAS phase for debugging purposes */
DEFINE_LOOKUP_FUNC(PhenomXOnlyReturnPhase, INT4, "PhenomXOnlyReturnPhase", 0)
DEFINE_LOOKUP_FUNC(PhenomXPNRInterpTolerance, REAL8, "PNRInterpTolerance", 0.01)

/* IMRPhenomX_PNR_Asymmetry Parameters */
DEFINE_LOOKUP_FUNC(PhenomXAntisymmetricWaveform, INT4, "AntisymmetricWaveform", 0)

/* IMRPhenomXCP Parameters */
DEFINE_LOOKUP_FUNC(PhenomXCPMU1, REAL8, "MU1", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPMU2, REAL8, "MU2", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPMU3, REAL8, "MU3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPMU4, REAL8, "MU4", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPNU0, REAL8, "NU0", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPNU4, REAL8, "NU4", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPNU5, REAL8, "NU5", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPNU6, REAL8, "NU6", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPZETA1, REAL8, "ZETA1", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPZETA2, REAL8, "ZETA2", 0)
/* l=3, m=3 */
DEFINE_LOOKUP_FUNC(PhenomXCPMU1l3m3, REAL8, "MU1l3m3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPMU2l3m3, REAL8, "MU2l3m3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPMU3l3m3, REAL8, "MU3l3m3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPMU4l3m3, REAL8, "MU4l3m3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPNU0l3m3, REAL8, "NU0l3m3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPNU4l3m3, REAL8, "NU4l3m3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPNU5l3m3, REAL8, "NU5l3m3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPNU6l3m3, REAL8, "NU6l3m3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPZETA1l3m3, REAL8, "ZETA1l3m3", 0)
DEFINE_LOOKUP_FUNC(PhenomXCPZETA2l3m3, REAL8, "ZETA2l3m3", 0)

/* ISDEFAULT FUNCTIONS */

DEFINE_ISDEFAULT_FUNC(Lmax, INT4, "lmax", 0)
DEFINE_ISDEFAULT_FUNC(ModesChoice, INT4, "modes", LAL_SIM_INSPIRAL_MODES_CHOICE_ALL)
DEFINE_ISDEFAULT_FUNC(FrameAxis, INT4, "axis", LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L)
DEFINE_ISDEFAULT_FUNC(Sideband, INT4, "sideband", 0)
DEFINE_ISDEFAULT_FUNC(NumRelData, String, "numreldata", NULL)

int XLALSimInspiralWaveformParamsModeArrayIsDefault(LALDict *params)
{
	return XLALSimInspiralWaveformParamsLookupModeArray(params) == NULL;
}

int XLALSimInspiralWaveformParamsModeArrayJframeIsDefault(LALDict *params)
{
	return XLALSimInspiralWaveformParamsLookupModeArrayJframe(params) == NULL;
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
DEFINE_ISDEFAULT_FUNC(NonGRDChiMinus2, REAL8, "dchiMinus2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChiMinus1, REAL8, "dchiMinus1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi0, REAL8, "dchi0", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi1, REAL8, "dchi1", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi2, REAL8, "dchi2", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi3, REAL8, "dchi3", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi3S, REAL8, "dchi3S", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi3NS, REAL8, "dchi3NS", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi4, REAL8, "dchi4", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi4S, REAL8, "dchi4S", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi4NS, REAL8, "dchi4NS", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi5, REAL8, "dchi5", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi5S, REAL8, "dchi5S", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi5NS, REAL8, "dchi5NS", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi5L, REAL8, "dchi5l", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi5LS, REAL8, "dchi5lS", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi5LNS, REAL8, "dchi5lNS", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi6, REAL8, "dchi6", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi6S, REAL8, "dchi6S", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi6NS, REAL8, "dchi6NS", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi6L, REAL8, "dchi6l", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi7, REAL8, "dchi7", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi7S, REAL8, "dchi7S", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChi7NS, REAL8, "dchi7NS", 0)
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
DEFINE_ISDEFAULT_FUNC(NonGRDChikappaS, REAL8, "dchikappaS", 0)
DEFINE_ISDEFAULT_FUNC(NonGRDChikappaA, REAL8, "dchikappaA", 0)
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
DEFINE_ISDEFAULT_FUNC(PhenomXPrecVersion, INT4, "PrecVersion", 300)
DEFINE_ISDEFAULT_FUNC(PhenomXReturnCoPrec, INT4, "ReturnCoPrec", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPExpansionOrder, INT4, "ExpansionOrder", 5)
DEFINE_ISDEFAULT_FUNC(PhenomXPConvention, INT4, "Convention", 1)
DEFINE_ISDEFAULT_FUNC(PhenomXPFinalSpinMod, INT4, "FinalSpinMod", 4)
DEFINE_ISDEFAULT_FUNC(PhenomXPTransPrecessionMethod, INT4, "TransPrecessionMethod", 1)
DEFINE_ISDEFAULT_FUNC(PhenomXPSpinTaylorVersion, String, "SpinTaylorVersion", NULL)
DEFINE_ISDEFAULT_FUNC(PhenomXPSpinTaylorCoarseFactor, INT4, "SpinTaylorCoarseFactor",10);

/* IMRPhenomX_NRTidal Parameters */
DEFINE_ISDEFAULT_FUNC(PhenomXTidalFlag, INT4, "PhenXTidal", 0)

/* IMRPhenomXHM Parameters */
DEFINE_ISDEFAULT_FUNC(PhenomXHMReleaseVersion, INT4, "PhenomXHMReleaseVersion", 122022)
DEFINE_ISDEFAULT_FUNC(PhenomXHMInspiralPhaseVersion, INT4, "InsPhaseHMVersion", 122019)
DEFINE_ISDEFAULT_FUNC(PhenomXHMIntermediatePhaseVersion, INT4, "IntPhaseHMVersion", 122019)
DEFINE_ISDEFAULT_FUNC(PhenomXHMRingdownPhaseVersion, INT4, "RDPhaseHMVersion", 122019)
DEFINE_ISDEFAULT_FUNC(PhenomXHMInspiralAmpVersion, INT4, "InsAmpHMVersion", 3)
DEFINE_ISDEFAULT_FUNC(PhenomXHMIntermediateAmpVersion, INT4, "IntAmpHMVersion", 2)
DEFINE_ISDEFAULT_FUNC(PhenomXHMRingdownAmpVersion, INT4, "RDAmpHMVersion", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXHMInspiralAmpFitsVersion, INT4, "InsAmpFitsVersion", 122018)
DEFINE_ISDEFAULT_FUNC(PhenomXHMIntermediateAmpFitsVersion, INT4, "IntAmpFitsVersion", 122018)
DEFINE_ISDEFAULT_FUNC(PhenomXHMRingdownAmpFitsVersion, INT4, "RDAmpFitsVersion", 122018)
DEFINE_ISDEFAULT_FUNC(PhenomXHMInspiralAmpFreqsVersion, INT4, "InsAmpFreqsVersion", 122018)
DEFINE_ISDEFAULT_FUNC(PhenomXHMIntermediateAmpFreqsVersion, INT4, "IntAmpFreqsVersion", 122018)
DEFINE_ISDEFAULT_FUNC(PhenomXHMRingdownAmpFreqsVersion, INT4, "RDAmpFreqsVersion", 122018)
DEFINE_ISDEFAULT_FUNC(PhenomXHMPhaseRef21, REAL8, "PhaseRef21", 0.)
DEFINE_ISDEFAULT_FUNC(PhenomXHMThresholdMband, REAL8, "ThresholdMband", 0.001)
DEFINE_ISDEFAULT_FUNC(PhenomXHMAmpInterpolMB, INT4, "AmpInterpol", 1)
DEFINE_ISDEFAULT_FUNC(DOmega220, REAL8, "domega220", 0)
DEFINE_ISDEFAULT_FUNC(DTau220, REAL8, "dtau220", 0)
DEFINE_ISDEFAULT_FUNC(DOmega210, REAL8, "domega210", 0)
DEFINE_ISDEFAULT_FUNC(DTau210, REAL8, "dtau210", 0)
DEFINE_ISDEFAULT_FUNC(DOmega330, REAL8, "domega330", 0)
DEFINE_ISDEFAULT_FUNC(DTau330, REAL8, "dtau330", 0)
DEFINE_ISDEFAULT_FUNC(DOmega440, REAL8, "domega440", 0)
DEFINE_ISDEFAULT_FUNC(DTau440, REAL8, "dtau440", 0)
DEFINE_ISDEFAULT_FUNC(DOmega550, REAL8, "domega550", 0)
DEFINE_ISDEFAULT_FUNC(DTau550, REAL8, "dtau550", 0)

/* IMRPhenomXPHM */
DEFINE_ISDEFAULT_FUNC(PhenomXPHMMBandVersion, INT4, "MBandPrecVersion", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMThresholdMband, REAL8, "PrecThresholdMband", 0.001)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMUseModes, INT4, "UseModes", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMModesL0Frame, INT4, "ModesL0Frame", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMPrecModes, INT4, "PrecModes", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPHMTwistPhenomHM, INT4, "TwistPhenomHM", 0)

/* IMRPhenomTHM Parameters */
DEFINE_ISDEFAULT_FUNC(PhenomTHMInspiralVersion, INT4, "InspiralVersion", 0)
DEFINE_ISDEFAULT_FUNC(PhenomTPHMMergerVersion, INT4, "MergerVersion", 1)

/* IMRPhenomX_PNR Parameters */
DEFINE_ISDEFAULT_FUNC(PhenomXPNRUseTunedAngles, INT4, "PNRUseTunedAngles", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPNRUseTunedCoprec, INT4, "PNRUseTunedCoprec", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPNRUseTunedCoprec33, INT4, "PNRUseTunedCoprec33", 0)
// Option to only be used when actively tuning PNR Coprec relative to XHM wherein the non-precessing final spin is used
DEFINE_ISDEFAULT_FUNC(PhenomXPNRUseInputCoprecDeviations, INT4, "PNRUseInputCoprecDeviations", 0)
// Dev option for forcing 22 phase derivative inspiral values to align with XHM at a low ref frequency
DEFINE_ISDEFAULT_FUNC(PhenomXPNRForceXHMAlignment, INT4, "PNRForceXHMAlignment", 0)
/* Toggle output of XAS phase for debugging purposes */
DEFINE_ISDEFAULT_FUNC(PhenomXOnlyReturnPhase, INT4, "PhenomXOnlyReturnPhase", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXPNRInterpTolerance, REAL8, "PNRInterpTolerance", 0.01)

/* IMRPhenomX_PNR_Asymmetry Parameters */
DEFINE_ISDEFAULT_FUNC(PhenomXAntisymmetricWaveform, INT4, "AntisymmetricWaveform", 0)

/* IMRPhenomXCP Parameters */
DEFINE_ISDEFAULT_FUNC(PhenomXCPMU1, REAL8, "MU1", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPMU2, REAL8, "MU2", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPMU3, REAL8, "MU3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPMU4, REAL8, "MU4", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPNU0, REAL8, "NU0", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPNU4, REAL8, "NU4", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPNU5, REAL8, "NU5", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPNU6, REAL8, "NU6", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPZETA1, REAL8, "ZETA1", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPZETA2, REAL8, "ZETA2", 0)
/* l=3, m=3 */
DEFINE_ISDEFAULT_FUNC(PhenomXCPMU1l3m3, REAL8, "MU1l3m3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPMU2l3m3, REAL8, "MU2l3m3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPMU3l3m3, REAL8, "MU3l3m3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPMU4l3m3, REAL8, "MU4l3m3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPNU0l3m3, REAL8, "NU0l3m3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPNU4l3m3, REAL8, "NU4l3m3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPNU5l3m3, REAL8, "NU5l3m3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPNU6l3m3, REAL8, "NU6l3m3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPZETA1l3m3, REAL8, "ZETA1l3m3", 0)
DEFINE_ISDEFAULT_FUNC(PhenomXCPZETA2l3m3, REAL8, "ZETA2l3m3", 0)


/** Insert common waveform parameters into LALDict */
LALDict* XLALSimInspiralParamsDict(const REAL8 m1, const REAL8 m2, const REAL8 S1x, const REAL8 S1y, const REAL8 S1z, const REAL8 S2x, const REAL8 S2y, const REAL8 S2z, const REAL8 distance, const REAL8 inclination, const REAL8 phiRef, const REAL8 longAscNodes, const REAL8 eccentricity, const REAL8 f_ref, LALDict *LALparams)

{
       LALDict* params = XLALDictDuplicate(LALparams);

   {   XLALSimInspiralWaveformParamsInsertMass1(params, m1);
       XLALSimInspiralWaveformParamsInsertMass2(params, m2);
       XLALSimInspiralWaveformParamsInsertSpin1x(params, S1x);
       XLALSimInspiralWaveformParamsInsertSpin1y(params, S1y);
       XLALSimInspiralWaveformParamsInsertSpin1z(params, S1z);
       XLALSimInspiralWaveformParamsInsertSpin2x(params, S2x);
       XLALSimInspiralWaveformParamsInsertSpin2y(params, S2y);
       XLALSimInspiralWaveformParamsInsertSpin2z(params, S2z);

       XLALSimInspiralWaveformParamsInsertDistance(params, distance);
       XLALSimInspiralWaveformParamsInsertInclination(params, inclination);
       XLALSimInspiralWaveformParamsInsertRefPhase(params, phiRef);
       XLALSimInspiralWaveformParamsInsertLongAscNodes(params, longAscNodes);
       XLALSimInspiralWaveformParamsInsertEccentricity(params, eccentricity);
       XLALSimInspiralWaveformParamsInsertF22Ref(params, f_ref);

}



    return params;
}

#undef String
