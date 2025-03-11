/* Copyright (C) 2012 Evan Ochsner
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

#include <stdio.h>
#include <lal/LALString.h>
#include <lal/LALValue.h>
#include <lal/LALStdlib.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include <lal/LALSimInspiralWaveformFlags.h>
#include <lal/LALSimInspiralWaveformParams.h>

/**
 * Struct containing several enumerated flags that control specialized behavior
 * for some waveform approximants.
 *
 * Users: Access this struct only through the Create/Destroy/Set/Get/IsDefault
 * functions declared in this file.
 *
 * Developers: Do not add anything but enumerated flags to this struct. Avoid
 * adding extra flags whenever possible.
 * DEPRECATED, use LALDict instead.
 */
struct tagLALSimInspiralWaveformFlags
{
    LALSimInspiralSpinOrder spinO; /**< PN order of spin effects */
    LALSimInspiralTidalOrder tideO; /**< PN order of spin effects */
    LALSimInspiralFrameAxis axisChoice; /**< Flag to set frame z-axis convention */
    LALSimInspiralModesChoice modesChoice; /**< Flag to control which modes are included in IMR models */
    char numreldata[FILENAME_MAX]; /**< Location of NR data file for NR waveforms */
};

/**
 * @addtogroup LALSimInspiralWaveformFlags_c
 * @brief Routines to manipulate inspiral waveform flags structures.
 * @{
 */

/**
 * Create a new LALSimInspiralWaveformFlags struct
 * with all flags set to their default values.
 *
 * Remember to destroy the struct when you are done with it.
 */
LALSimInspiralWaveformFlags *XLALSimInspiralCreateWaveformFlags(void)
{
    LALSimInspiralWaveformFlags *waveFlags;
    /* Allocate memory for the waveform flags */
    waveFlags = XLALMalloc( sizeof(*waveFlags) );
    if( !waveFlags )
    {
        XLALFree(waveFlags);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }

    /* Set all flags to their default values */
    XLALSimInspiralSetSpinOrder(waveFlags,
            LAL_SIM_INSPIRAL_SPIN_ORDER_DEFAULT);
    XLALSimInspiralSetTidalOrder(waveFlags,
            LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT);
    XLALSimInspiralSetFrameAxis(waveFlags,
            LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT);
    XLALSimInspiralSetModesChoice(waveFlags,
            LAL_SIM_INSPIRAL_MODES_CHOICE_DEFAULT);

    return waveFlags;
}

/**
 * Destroy a LALSimInspiralWaveformFlags struct.
 */
void XLALSimInspiralDestroyWaveformFlags(
        LALSimInspiralWaveformFlags *waveFlags
        )
{
    if( waveFlags )
        XLALFree(waveFlags);
    return;
}

/**
 * Returns true if waveFlags is non-NULL and all of its fields have default
 * value; returns false otherwise.
 */
bool XLALSimInspiralWaveformParamsFlagsAreDefault(LALDict *params)
{
    /* Check every field of WaveformFlags, each returns 1/0 for true/false.
     * Return true iff waveFlags is non-NULL and all checks are true. */
  return ( XLALSimInspiralWaveformParamsPNSpinOrderIsDefault(params) &&
	   XLALSimInspiralWaveformParamsPNTidalOrderIsDefault(params) &&
	   XLALSimInspiralWaveformParamsFrameAxisIsDefault(params) &&
	   XLALSimInspiralWaveformParamsModesChoiceIsDefault(params));
}

/**
 * Returns true if waveFlags is non-NULL and all of its fields have default
 * value; returns false otherwise.
 */
bool XLALSimInspiralWaveformFlagsIsDefaultOLD(
        LALSimInspiralWaveformFlags *waveFlags
        )
{
    /* Check every field of WaveformFlags, each returns 1/0 for true/false.
     * Return true iff waveFlags is non-NULL and all checks are true. */
    return !waveFlags || (
        XLALSimInspiralSpinOrderIsDefault(waveFlags->spinO) &&
        XLALSimInspiralTidalOrderIsDefault(waveFlags->tideO) &&
        XLALSimInspiralFrameAxisIsDefault(waveFlags->axisChoice) &&
        XLALSimInspiralModesChoiceIsDefault(waveFlags->modesChoice));
}

/**
 * Checks if all flags in two LALSimInspiralWaveformFlags structs are equal.
 * Returns true if all flags are equal. Returns false if one or more differ.
 */
bool XLALSimInspiralWaveformFlagsEqualOLD(
        LALSimInspiralWaveformFlags *waveFlags1,
        LALSimInspiralWaveformFlags *waveFlags2
        )
{
    LALSimInspiralSpinOrder spinO1, spinO2;
    LALSimInspiralTidalOrder tideO1, tideO2;
    LALSimInspiralFrameAxis axisChoice1, axisChoice2;
    LALSimInspiralModesChoice modesChoice1, modesChoice2;

    spinO1 = XLALSimInspiralGetSpinOrder(waveFlags1);
    spinO2 = XLALSimInspiralGetSpinOrder(waveFlags2);
    tideO1 = XLALSimInspiralGetTidalOrder(waveFlags1);
    tideO2 = XLALSimInspiralGetTidalOrder(waveFlags2);
    axisChoice1 = XLALSimInspiralGetFrameAxis(waveFlags1);
    axisChoice2 = XLALSimInspiralGetFrameAxis(waveFlags2);
    modesChoice1 = XLALSimInspiralGetModesChoice(waveFlags1);
    modesChoice2 = XLALSimInspiralGetModesChoice(waveFlags2);

    return ( (spinO1==spinO2) && (tideO1==tideO2) && (axisChoice1==axisChoice2)
            && (modesChoice1==modesChoice2) );
}

/**
 * Checks if all flags in two LALSimInspiralWaveformFlags structs are equal.
 * Returns true if all flags are equal. Returns false if one or more differ.
 */
bool XLALSimInspiralWaveformFlagsEqual(
        LALDict *LALpars1,
        LALDict *LALpars2
        )
{
    LALSimInspiralSpinOrder spinO1, spinO2;
    LALSimInspiralTidalOrder tideO1, tideO2;
    LALSimInspiralFrameAxis axisChoice1, axisChoice2;
    LALSimInspiralModesChoice modesChoice1, modesChoice2;

    spinO1 = XLALSimInspiralWaveformParamsLookupPNSpinOrder(LALpars1);
    spinO2 = XLALSimInspiralWaveformParamsLookupPNSpinOrder(LALpars2);
    tideO1 = XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALpars1);
    tideO2 = XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALpars2);
    axisChoice1 = XLALSimInspiralWaveformParamsLookupFrameAxis(LALpars1);
    axisChoice2 = XLALSimInspiralWaveformParamsLookupFrameAxis(LALpars2);
    modesChoice1 = XLALSimInspiralWaveformParamsLookupModesChoice(LALpars1);
    modesChoice2 = XLALSimInspiralWaveformParamsLookupModesChoice(LALpars2);

    return ( (spinO1==spinO2) && (tideO1==tideO2) && (axisChoice1==axisChoice2)
            && (modesChoice1==modesChoice2) );
}

/**
 * Set the LALSimInspiralSpinOrder within a LALSimInspiralWaveformFlags struct
 */
void XLALSimInspiralSetSpinOrder(
        LALSimInspiralWaveformFlags *waveFlags, /**< Struct whose flag will be set */

        LALSimInspiralSpinOrder spinO /**< value to set flag to */
        )
{
    waveFlags->spinO = spinO;
    return;
}

/**
 * Get the LALSimInspiralSpinOrder within a LALSimInspiralWaveformFlags struct,
 * or LAL_SIM_INSPIRAL_SPIN_ORDER_DEFAULT if waveFlags is NULL
 */
LALSimInspiralSpinOrder XLALSimInspiralGetSpinOrder(
        LALSimInspiralWaveformFlags *waveFlags
        )
{
    if ( waveFlags )
        return waveFlags->spinO;
    else
        return LAL_SIM_INSPIRAL_SPIN_ORDER_DEFAULT;
}

/**
 * Returns true if LALSimInspiralSpinOrder has default value
 * returns false otherwise
 */
bool XLALSimInspiralSpinOrderIsDefault(
        LALSimInspiralSpinOrder spinO
        )
{
    if( spinO == LAL_SIM_INSPIRAL_SPIN_ORDER_DEFAULT )
        return true;
    else
        return false;
}

/**
 * Set the LALSimInspiralTidalOrder within a LALSimInspiralWaveformFlags struct
 */
void XLALSimInspiralSetTidalOrder(
        LALSimInspiralWaveformFlags *waveFlags, /**< Struct whose flag will be set */

        LALSimInspiralTidalOrder tideO /**< value to set flag to */
        )
{
    waveFlags->tideO = tideO;
    return;
}

/**
 * Get the LALSimInspiralTidalOrder within a LALSimInspiralWaveformFlags struct,
 * or LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT if waveFlags is NULL
 */
LALSimInspiralTidalOrder XLALSimInspiralGetTidalOrder(
        LALSimInspiralWaveformFlags *waveFlags
        )
{
    if ( waveFlags )
        return waveFlags->tideO;
    else
        return LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT;
}

/**
 * Returns true if LALSimInspiralTidalOrder has default value
 * returns false otherwise
 */
bool XLALSimInspiralTidalOrderIsDefault(
        LALSimInspiralTidalOrder tideO
        )
{
    if( tideO == LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT )
        return true;
    else
        return false;
}

/**
 * Set the LALSimInspiralFrameAxis within a LALSimInspiralWaveformFlags struct
 */
void XLALSimInspiralSetFrameAxis(
        LALSimInspiralWaveformFlags *waveFlags, /**< Struct whose flag will be set */
        LALSimInspiralFrameAxis axisChoice /**< value to set flag to */
        )
{
    waveFlags->axisChoice = axisChoice;
    return;
}

/**
 * Get the LALSimInspiralFrameAxis within a LALSimInspiralWaveformFlags struct,
 * or LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT if waveFlags is NULL
 */
LALSimInspiralFrameAxis XLALSimInspiralGetFrameAxis(
        LALSimInspiralWaveformFlags *waveFlags
        )
{
    if ( waveFlags )
        return waveFlags->axisChoice;
    else
        return LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT;
}

/**
 * Returns true if LALSimInspiralFrameAxis has default value
 * returns false otherwise
 */
bool XLALSimInspiralFrameAxisIsDefault(
        LALSimInspiralFrameAxis axisChoice
        )
{
    if( axisChoice == LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT )
        return true;
    else
        return false;
}

/**
 * Set the LALSimInspiralModesChoice within a LALSimInspiralWaveformFlags struct
 */
void XLALSimInspiralSetModesChoice(
        LALSimInspiralWaveformFlags *waveFlags, /**< Struct whose flag will be set */
        LALSimInspiralModesChoice modesChoice /**< value to set flag to */
        )
{
    waveFlags->modesChoice = modesChoice;
    return;
}

/**
 * Get the LALSimInspiralModesChoice within a LALSimInspiralWaveformFlags struct,
 * or LAL_SIM_INSPIRAL_MODES_CHOICE_DEFAULT if waveFlags is NULL
 */
LALSimInspiralModesChoice XLALSimInspiralGetModesChoice(
        LALSimInspiralWaveformFlags *waveFlags
        )
{
    if ( waveFlags )
        return waveFlags->modesChoice;
    else
        return LAL_SIM_INSPIRAL_MODES_CHOICE_DEFAULT;
}

/**
 * Returns true if LALSimInspiralModesChoice has default value
 * returns false otherwise
 */
bool XLALSimInspiralModesChoiceIsDefault(
        LALSimInspiralModesChoice modesChoice
        )
{
    if( modesChoice == LAL_SIM_INSPIRAL_MODES_CHOICE_DEFAULT )
        return true;
    else
        return false;
}

/**
 * Set the numreldata string within a LALSimInspiralWaveformFlags struct
 */
void XLALSimInspiralSetNumrelDataOLD(
        LALSimInspiralWaveformFlags *waveFlags, /**< Struct whose value will be set */
        const char* numreldata /**< value to set numreldata to */
        )
{
    XLALStringCopy(waveFlags->numreldata, numreldata, sizeof(waveFlags->numreldata));
    return;
}

/**
 * Returns a deepcopy of the pointer of the numeraldata attribute of the
 * waveFlags structure. If this is NULL then NULL will be returned.
 * The returned value is independent of the waveFlags structure and will
 * need to be LALFree-d.
 */
char* XLALSimInspiralGetNumrelDataOLD(
        LALSimInspiralWaveformFlags *waveFlags
        )
{
    char *ret_string;
    if ( waveFlags )
    {
        ret_string = XLALMalloc(FILENAME_MAX * sizeof(char));
        XLALStringCopy(ret_string, waveFlags->numreldata, sizeof(waveFlags->numreldata));
        return ret_string;
    }
    else
    {
        return NULL;
    }
}


static char empty_modes[((LAL_SIM_L_MAX_MODE_ARRAY + 1) * (LAL_SIM_L_MAX_MODE_ARRAY + 1)) / CHAR_BIT + 2] = { '\0' };

/**
 * Create a LALValue pointer to store the mode array.
 */

LALValue * XLALSimInspiralCreateModeArray(void)
{
	return XLALCreateValue(empty_modes, sizeof(empty_modes), LAL_CHAR_TYPE_CODE);
}

LALValue * XLALSimInspiralModeArrayActivateMode(LALValue *modes, unsigned l, int m)
{
	char *data;
	unsigned bit = l * l + l + m;
	unsigned byte = bit / CHAR_BIT;
	bit %= CHAR_BIT;

	/* sanity checks on l and m */
	XLAL_CHECK_NULL(l <= LAL_SIM_L_MAX_MODE_ARRAY, XLAL_EINVAL, "Invalid value of l=%u must not be greater than %u", l, LAL_SIM_L_MAX_MODE_ARRAY);
	XLAL_CHECK_NULL((unsigned)abs(m) <= l, XLAL_EINVAL, "Invalid value of m=%d for l=%u", m, l);

	/* sanity checks on modes */
	data = (char *)(intptr_t)XLALValueGetString(modes);
	XLAL_CHECK_NULL(data, XLAL_EFUNC);
	XLAL_CHECK_NULL(XLALValueGetSize(modes) == sizeof(empty_modes), XLAL_EINVAL, "Invalid data size for modes");

	data[byte] |= (1 << bit);
	return modes;
}

LALValue * XLALSimInspiralModeArrayDeactivateMode(LALValue *modes, unsigned l, int m)
{
	char *data;
	unsigned bit = l * l + l + m;
	unsigned byte = bit / CHAR_BIT;
	bit %= CHAR_BIT;

	/* sanity checks on l and m */
	XLAL_CHECK_NULL(l <= LAL_SIM_L_MAX_MODE_ARRAY, XLAL_EINVAL, "Invalid value of l=%u must not be greater than %u", l, LAL_SIM_L_MAX_MODE_ARRAY);
	XLAL_CHECK_NULL((unsigned)abs(m) <= l, XLAL_EINVAL, "Invalid value of m=%d for l=%u", m, l);

	/* sanity checks on modes */
	data = (char *)(intptr_t)XLALValueGetString(modes);
	XLAL_CHECK_NULL(data, XLAL_EFUNC);
	XLAL_CHECK_NULL(XLALValueGetSize(modes) == sizeof(empty_modes), XLAL_EINVAL, "Invalid data size for modes");

	data[byte] &= ~(1 << bit);
	return modes;
}

LALValue * XLALSimInspiralModeArrayActivateAllModes(LALValue *modes)
{
	char *data;
	data = (char *)(intptr_t)XLALValueGetString(modes);
	XLAL_CHECK_NULL(data, XLAL_EFUNC);
	XLAL_CHECK_NULL(XLALValueGetSize(modes) == sizeof(empty_modes), XLAL_EINVAL, "Invalid data size for modes");
	memset(data, ~0, sizeof(empty_modes) - 1);

    /* Deactivate the unphysical modes: (l,m) = ((0,0), (1,-1), (1,0), (1,1)) */
    XLALSimInspiralModeArrayDeactivateMode(modes, 0, 0);
    XLALSimInspiralModeArrayDeactivateMode(modes, 1, -1);
    XLALSimInspiralModeArrayDeactivateMode(modes, 1, 0);
    XLALSimInspiralModeArrayDeactivateMode(modes, 1, 1);

	return modes;
}

LALValue * XLALSimInspiralModeArrayDeactivateAllModes(LALValue *modes)
{
	char *data;
	data = (char *)(intptr_t)XLALValueGetString(modes);
	XLAL_CHECK_NULL(data, XLAL_EFUNC);
	XLAL_CHECK_NULL(XLALValueGetSize(modes) == sizeof(empty_modes), XLAL_EINVAL, "Invalid data size for modes");
	memset(data, 0, sizeof(empty_modes) - 1);
	return modes;
}

int XLALSimInspiralModeArrayIsModeActive(LALValue *modes, unsigned l, int m)
{
	const char *data;
	unsigned bit = l * l + l + m;
	unsigned byte = bit / CHAR_BIT;
	// unsigned bit;
	// unsigned byte;
	// positionOfModeInString(&bit, &byte, l, m);
	bit %= CHAR_BIT;

	/* sanity checks on l and m */
	XLAL_CHECK(l <= LAL_SIM_L_MAX_MODE_ARRAY, XLAL_EINVAL, "Invalid value of l=%u must not be greater than %u", l, LAL_SIM_L_MAX_MODE_ARRAY);
	XLAL_CHECK((unsigned)abs(m) <= l, XLAL_EINVAL, "Invalid value of m=%d for l=%u", m, l);

	/* sanity checks on modes */
	data = XLALValueGetString(modes);
	XLAL_CHECK(data, XLAL_EFUNC);
	XLAL_CHECK(XLALValueGetSize(modes) == sizeof(empty_modes), XLAL_EINVAL, "Invalid data size for modes");

	return (data[byte] & (1 << bit)) != 0;
}


LALValue * XLALSimInspiralModeArrayActivateAllModesAtL(LALValue *modes, unsigned l)
{
	for(int m =-l; m <= (int) l; ++m)
	{
		XLALSimInspiralModeArrayActivateMode(modes, l, m);
	}
	return modes;
}

LALValue * XLALSimInspiralModeArrayDeactivateAllModesAtL(LALValue *modes, unsigned l)
{
	for(int m =-l; m <= (int) l; ++m)
	{
		XLALSimInspiralModeArrayDeactivateMode(modes, l, m);
	}
	return modes;
}

int XLALSimInspiralModeArrayPrintModes(LALValue *modes)
{
	int l;
	for (l = 0; l <= LAL_SIM_L_MAX_MODE_ARRAY; ++l) {
		int m;
		for (m = -l; m <= l; ++m)
			printf("(%u,%+d) : %d\n", l, m, XLALSimInspiralModeArrayIsModeActive(modes, l, m));
		printf("\n");
	}
	return 0;
}

INT2Sequence *XLALSimInspiralModeArrayReadModes(LALValue *modes)
{
	INT2Sequence *seqmodes = XLALCreateINT2Sequence(4*LAL_SIM_L_MAX_MODE_ARRAY+2);
	int nmodes = 0;
	for (int l = 0; l <= LAL_SIM_L_MAX_MODE_ARRAY; l++) {
		for (int m = -l; m <= l; m++)
			 if(XLALSimInspiralModeArrayIsModeActive(modes, l, m)==1)
			 {
				 seqmodes->data[2*nmodes] = l;
				 seqmodes->data[2*nmodes+1] = m;
				 nmodes++;
			 }
	}
	seqmodes = XLALShrinkINT2Sequence(seqmodes, 0, 2*nmodes);
	return seqmodes;
}

char * XLALSimInspiralModeArrayToModeString(LALValue *modes)
{
    char *s = NULL;
    int n = 0;
    if ((s = XLALStringAppend(s, "[")) == NULL)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    for (int l = 0; l <= LAL_SIM_L_MAX_MODE_ARRAY; ++l)
        for (int m = -l; m <= l; ++m)
            if (XLALSimInspiralModeArrayIsModeActive(modes, l, m))
                if ((s = XLALStringAppendFmt(s, "%s(%u,%+d)", n++ ? "," : "", l, m)) == NULL)
                    XLAL_ERROR_NULL(XLAL_EFUNC);
    if ((s = XLALStringAppend(s, "]")) == NULL)
        XLAL_ERROR_NULL(XLAL_EFUNC);
    return s;
}

LALValue * XLALSimInspiralModeArrayFromModeString(const char *modestr)
{
    LALValue *modes = NULL;
    const char *s = modestr;
    int inlist = 0;
    int intup = 0;
    int done = 0;
    unsigned l = UINT_MAX; // invalid value
    int m = INT_MAX; // invalid value
    char c;

    XLAL_CHECK_NULL(modestr, XLAL_EFAULT);

    modes = XLALSimInspiralCreateModeArray();
    XLAL_CHECK_NULL(modes, XLAL_ENOMEM);

    while ((c = *s++)) {
        char *endp;

        if (isspace(c))
            continue;

        XLAL_CHECK_FAIL(!done, XLAL_EINVAL, "Malformed mode string \"%s\": trailing characters found", modestr);

        switch (c) {

        case '[':
            XLAL_CHECK_FAIL(!inlist, XLAL_EINVAL, "Malformed mode string \"%s\": cannot have a list inside a list", modestr);
            inlist = 1;
            break;

        case ']':
            XLAL_CHECK_FAIL(inlist, XLAL_EINVAL, "Malformed mode string \"%s\": end of list when not in list", modestr);
            inlist = 0;
            done = 1;
            break;

        case '(':
            XLAL_CHECK_FAIL(inlist, XLAL_EINVAL, "Malformed mode string \"%s\": tuple found outside of list", modestr);
            XLAL_CHECK_FAIL(!intup, XLAL_EINVAL, "Malformed mode string \"%s\": tuple found within tuple", modestr);
            intup = 1;
            l = strtoul(s, &endp, 0);
            XLAL_CHECK_FAIL(s != endp, XLAL_EINVAL, "Malformed mode string \"%s\": could not convert unsigned integer", modestr);
            s = endp;
            break;

        case ')':
            XLAL_CHECK_FAIL(intup, XLAL_EINVAL, "Malformed mode string \"%s\": end of tuple when not in tuple", modestr);
            XLAL_CHECK_FAIL(inlist, XLAL_EINVAL, "Malformed mode string \"%s\": tuple found outside of list", modestr);
            intup = 0;
            break;

        case ',':
            XLAL_CHECK_FAIL(inlist, XLAL_EINVAL, "Malformed mode string \"%s\": separater found when not in list", modestr);
            if (intup) {
                m = strtol(s, &endp, 0);
                XLAL_CHECK_FAIL(s != endp, XLAL_EINVAL, "Malformed mode string \"%s\": could not convert signed integer", modestr);
                XLAL_CHECK_FAIL(XLALSimInspiralModeArrayActivateMode(modes, l, m), XLAL_EFUNC);
                s = endp;
            }
            break;

        default:
            XLAL_ERROR_FAIL(XLAL_EINVAL, "Invalid character '%c' in mode string \"%s\"", c, modestr);
        }
    }

    if (done)
        return modes;

XLAL_FAIL:
    XLALDestroyValue(modes);
    return NULL;
}

/** @} */
