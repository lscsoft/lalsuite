/*
 * Macro procedure for aborting if non-default LALSimInspiralWaveformFlags
 * struct was provided, but that approximant does not use the struct
 * and only has a single default use case.
 *
 * The ChooseWaveform functions will fail in such a case, so the user does not
 * think they are including features that are unavailable.
 *
 * All of the macros below will destroy the LALSimInspiralWaveformFlags struct,
 * print a specific warning and raise a general XLAL_ERROR for invalid argument.
 */
#define ABORT_NONDEFAULT_WAVEFORM_FLAGS(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Same as above macro, but returns a null pointer rather than XLAL_FAILURE int
 */
#define ABORT_NONDEFAULT_WAVEFORM_FLAGS_NULL(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralWaveformFlags given, but this approximant does not support this case.\n", __func__);\
XLAL_ERROR_NULL(XLAL_EINVAL);\
} while (0)

#define ABORT_NONDEFAULT_LALDICT_FLAGS(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-default flags given, but this approximant does not support this case.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Same as above macro, but returns a null pointer rather than XLAL_FAILURE int
 */
#define ABORT_NONDEFAULT_LALDICT_FLAGS_NULL(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-default flags given, but this approximant does not support this case.\n", __func__);\
XLAL_ERROR_NULL(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-zero spins
 * given to a non-spinning approximant
 */
#define ABORT_NONZERO_SPINS_OLD(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-zero spins were given, but this is a non-spinning approximant.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-zero spins
 * given to a non-spinning approximant
 */
#define ABORT_NONZERO_SPINS(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-zero spins were given, but this is a non-spinning approximant.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-zero transverse spin
 * components given to a non-precessing approximant
 */
#define ABORT_NONZERO_TRANSVERSE_SPINS_OLD(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-zero transverse spins were given, but this is a non-precessing approximant.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-zero transverse spin
 * components given to a non-precessing approximant
 */
#define ABORT_NONZERO_TRANSVERSE_SPINS(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-zero transverse spins were given, but this is a non-precessing approximant.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-zero tidal parameters
 * given to an approximant with no tidal corrections
 */
#define ABORT_NONZERO_TIDES_OLD(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-zero tidal parameters
 * given to an approximant with no tidal corrections
 */
#define ABORT_NONZERO_TIDES(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-zero tidal parameters were given, but this is approximant doe not have tidal corrections.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralSpinOrder is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_SPIN_ORDER_OLD(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralSpinOrder is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_SPIN_ORDER(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralSpinOrder provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralTidalOrder is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_TIDAL_ORDER_OLD(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralTidalOrder provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralTidalOrder is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_TIDAL_ORDER(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralTidalOrder provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralFrameAxis is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_FRAME_AXIS_OLD(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Same as above macro, but returns a null pointer rather than XLAL_FAILURE int
 */
#define ABORT_NONDEFAULT_FRAME_AXIS_NULL_OLD(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR_NULL(XLAL_EINVAL);\
} while (0)

#define ABORT_NONDEFAULT_FRAME_AXIS(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Same as above macro, but returns a null pointer rather than XLAL_FAILURE int
 */
#define ABORT_NONDEFAULT_FRAME_AXIS_NULL(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralFrameAxis provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR_NULL(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-default value of
 * LALSimInspiralModesChoice is given for an approximant
 * which does not use that flag
 */
#define ABORT_NONDEFAULT_MODES_CHOICE_OLD(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

#define ABORT_NONDEFAULT_MODES_CHOICE(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR(XLAL_EINVAL);\
} while (0)

/*
 * Same as above macro, but returns a null pointer rather than XLAL_FAILURE int
 */
#define ABORT_NONDEFAULT_MODES_CHOICE_NULL_OLD(waveFlags)\
do {\
XLALSimInspiralDestroyWaveformFlags(waveFlags);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR_NULL(XLAL_EINVAL);\
} while (0)

/*
 * Same as above macro, but returns a null pointer rather than XLAL_FAILURE int
 */
#define ABORT_NONDEFAULT_MODES_CHOICE_NULL(LALparams)\
do {\
XLALDestroyDict(LALparams);\
XLALPrintError("XLAL Error - %s: Non-default LALSimInspiralModesChoice provided, but this approximant does not use that flag.\n", __func__);\
XLAL_ERROR_NULL(XLAL_EINVAL);\
} while (0)

/*
 * Macro procedure for aborting if non-zero spin2
 * components given to central-spin-only approximant
 */
#define ABORT_NONZERO_SPIN2_OLD(waveFlags)\
	do {\
	XLALSimInspiralDestroyWaveformFlags(waveFlags);\
	XLALPrintError("XLAL Error - %s: Non-zero CO spin given, but this approximant does not support this case.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/*
 * Macro procedure for aborting if non-zero spin2
 * components given to central-spin-only approximant
 */
#define ABORT_NONZERO_SPIN2(LALparams)\
	do {\
	XLALDestroyDict(LALparams);\
	XLALPrintError("XLAL Error - %s: Non-zero CO spin given, but this approximant does not support this case.\n", __func__);\
	XLAL_ERROR(XLAL_EINVAL);\
	} while (0)

/* Internal utility macro to check all spin components are zero
   returns 1 if all spins zero, otherwise returns 0 */
#define checkSpinsZero(s1x, s1y, s1z, s2x, s2y, s2z) \
    (((s1x) != 0. || (s1y) != 0. || (s1z) != 0. || (s2x) != 0. || (s2y) != 0. || (s2z) != 0.) ? 0 : 1)

/* Internal utility macro to check that the second body's spin components are zero.
   Returns 1 if all components are zero, otherwise returns 0 */
#define checkCOSpinZero(s2x, s2y, s2z) \
    (((s2x) != 0. || (s2y) != 0. || (s2z) != 0.) ? 0 : 1)

/* Internal utility macro to check transverse spins are zero
   returns 1 if x and y components of spins are zero, otherwise returns 0 */
#define checkTransverseSpinsZero(s1x, s1y, s2x, s2y) \
    (((s1x) != 0. || (s1y) != 0. || (s2x) != 0. || (s2y) != 0. ) ? 0 : 1)

/* Internal utility macro to check aligned spins very close to equal
   returns 1 if z components of spins are very close to equal, otherwise returns 0 */
#define checkAlignedSpinsEqual(s1z, s2z) \
    ((fabs((s1z) - (s2z)) > 1e-6) ? 0 : 1)

/* Internal utility macro to check tidal parameters are zero
   returns 1 if both tidal parameters zero, otherwise returns 0 */
#define checkTidesZero(lambda1, lambda2) \
    (((lambda1) != 0. || (lambda2) != 0. ) ? 0 : 1)
