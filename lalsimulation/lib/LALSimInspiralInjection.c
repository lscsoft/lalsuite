/*
*  Copyright (C) 2023 Jolien Creighton
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
#include <string.h>
#include <sys/types.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDetectors.h>
#include <lal/LALString.h>
#include <lal/LALDict.h>
#include <lal/LALDictSequence.h>
#include <lal/Date.h>
#include <lal/Sort.h>
#include <lal/AVFactories.h>
#include <lal/StringVector.h>
#include <lal/LALSimulation.h>
#include <lal/LALSimInspiral.h>
#include <lal/H5FileIO.h>
#include <lal/LALSimInspiralInjection.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* translates from the injection file parameter names
 * https://git.ligo.org/waveforms/o4-injection-file/-/blob/main/ESSENTIAL_SPEC.md * to the waveform interface parameter names
 * https://git.ligo.org/waveforms/new-waveforms-interface/-/blob/master/parameter-names.md
 */
static const char * translate_key(const char *key, int inverse)
{
    static const char *map[][2] = {
        { "mass1_det", "mass1" },
        { "mass2_det", "mass2" },
        { "d_lum", "distance" },
        { "f22_ref_spin", "f22_ref" },
        { "coa_phase", "phi_ref" },
        { "cbc_model", "approximant" },
    };
    int to = inverse ? 0 : 1;
    int from = inverse ? 1 : 0;
    for (size_t i = 0; i < XLAL_NUM_ELEM(map); ++i)
        if (strcmp(key, map[i][from]) == 0)
            return map[i][to];
    return key;
}

static double si_scale_factor(const char *key)
{
    if (strcmp(key, "mass1") == 0)
        return LAL_MSUN_SI;
    else if (strcmp(key, "mass2") == 0)
        return LAL_MSUN_SI;
    else if (strcmp(key, "distance") == 0)
        return 1e6 * LAL_PC_SI;
    return 1.0;
}

LALDictSequence * XLALSimInspiralInjectionSequenceFromH5File(const char *fname)
{
    LALDictSequence *injseq = NULL;
    LALStringVector *strvec = NULL;
    INT4Vector *intvec = NULL;
    REAL8Vector *dblvec = NULL;
    LALH5File *file = NULL;
    LALH5File *group = NULL;
    LALH5Dataset *dset = NULL;
    char *name = NULL;
    size_t npar;

    XLAL_CHECK_NULL(fname, XLAL_EFAULT);
    XLAL_CHECK_NULL(*fname, XLAL_EINVAL);

    file = XLALH5FileOpen(fname, "r");
    XLAL_CHECK_FAIL(file, XLAL_EFUNC);
    group = XLALH5GroupOpen(file, "cbc_waveform_params");
    XLAL_CHECK_FAIL(group, XLAL_EFUNC);
    npar = XLALH5FileQueryNDatasets(group);
    XLAL_CHECK_FAIL((ssize_t)npar >= 0, XLAL_EFUNC);
    XLAL_CHECK_FAIL(npar > 0, XLAL_ESIZE, "No datasets found");

    for (size_t pos = 0; pos < npar; ++pos) {
        const char *key;
        int namelen;
        LALTYPECODE type;
        size_t ninj;
        double scale;

        namelen = XLALH5FileQueryDatasetName(NULL, 0, group, pos);
        XLAL_CHECK_FAIL(namelen >= 0, XLAL_EFUNC);
        XLAL_CHECK_FAIL(namelen > 0, XLAL_ENAME);
        name = LALMalloc(namelen + 1);
        XLAL_CHECK_FAIL(name, XLAL_ENOMEM);
        XLAL_CHECK_FAIL(XLALH5FileQueryDatasetName(name, namelen + 1, group, pos) == namelen, XLAL_EFUNC);
	dset = XLALH5DatasetRead(group, name);
        XLAL_CHECK_FAIL(dset, XLAL_EFUNC);
        type = XLALH5DatasetQueryType(dset);
        XLAL_CHECK_FAIL((int)type >= 0, XLAL_EFUNC);
        ninj = XLALH5DatasetQueryNPoints(dset);
        XLAL_CHECK_FAIL((ssize_t)ninj >= 0, XLAL_EFUNC);
        if (injseq == NULL) {
            injseq = XLALCreateDictSequence(ninj);
            XLAL_CHECK_FAIL(injseq, XLAL_EFUNC);
        } else if (ninj != injseq->length)
            XLAL_ERROR_FAIL(XLAL_EBADLEN, "Datasets have different number of points");

        key = strrchr(name, '/');
        key = key ? key + 1: name;
        key = translate_key(key, 0);

        switch (type) {
        case LAL_CHAR_TYPE_CODE:
           strvec = XLALH5DatasetReadStringVector(dset);
           XLAL_CHECK_FAIL(strvec, XLAL_EFUNC);
           for (size_t i = 0; i < strvec->length; ++i) {
               int retval;
               if (!strlen(strvec->data[i]))
                   continue;
               /* FIXME: Handle ModeArray when new waveform interface has been merged */
               //if (strcmp(key, "ModeArray") == 0) {
               //   retval = 0; /* handle mode array as a special case */
               //} else {
                   retval = XLALDictInsertStringValue(injseq->data[i], key, strvec->data[i]);
               //}
               XLAL_CHECK_FAIL(retval == 0, XLAL_EFUNC);
           }
           XLALDestroyStringVector(strvec);
           strvec = NULL;
           break;
        case LAL_I4_TYPE_CODE:
           intvec = XLALH5DatasetReadINT4Vector(dset);
           XLAL_CHECK_FAIL(intvec, XLAL_EFUNC);
           for (size_t i = 0; i < intvec->length; ++i) {
               int retval;
               retval = XLALDictInsertINT4Value(injseq->data[i], key, intvec->data[i]);
               XLAL_CHECK_FAIL(retval == 0, XLAL_EFUNC);
           }
           XLALDestroyINT4Vector(intvec);
           intvec = NULL;
           break;
        case LAL_D_TYPE_CODE:
           scale = si_scale_factor(key);
           dblvec = XLALH5DatasetReadREAL8Vector(dset);
           XLAL_CHECK_FAIL(dblvec, XLAL_EFUNC);
           for (size_t i = 0; i < dblvec->length; ++i) {
               int retval;
               if (!isfinite(dblvec->data[i]))
                   continue;
               retval = XLALDictInsertREAL8Value(injseq->data[i], key, scale * dblvec->data[i]);
               XLAL_CHECK_FAIL(retval == 0, XLAL_EFUNC);
           }
           XLALDestroyREAL8Vector(dblvec);
           dblvec = NULL;
           break;
        default:
           XLAL_ERROR_FAIL(XLAL_ETYPE, "Unsupported data type for dataset \"%s\"", name);
           break;
	}

        XLALH5DatasetFree(dset);
        XLALFree(name);
        dset = NULL;
        name = NULL;
    }

    XLALH5FileClose(group);
    XLALH5FileClose(file);

    return injseq;

XLAL_FAIL:

    XLALFree(name);
    XLALH5DatasetFree(dset);
    XLALH5FileClose(group);
    XLALH5FileClose(file);
    XLALDestroyREAL8Vector(dblvec);
    XLALDestroyINT4Vector(intvec);
    XLALDestroyStringVector(strvec);
    XLALDestroyDictSequence(injseq);
    return NULL;
}


static void XLALSimInspiralDictInsertParameterType(char *key, LALValue *value, void *thunk)
{
    LALDict *types = (LALDict *)thunk;
    XLALDictInsertINT4Value(types, key, XLALValueGetType(value));
    return;
}

int XLALSimInspiralInjectionSequenceToH5File(const LALDictSequence *injseq, const char *fname)
{
    LALH5File *file = NULL;
    LALH5File *group = NULL;
    LALStringVector *strvec = NULL;
    INT4Vector *intvec = NULL;
    REAL8Vector *dblvec = NULL;
    LALDict *types = NULL;
    LALDictIter iter;
    LALDictEntry *entry;

    file = XLALH5FileOpen(fname, "w");
    XLAL_CHECK_FAIL(file, XLAL_EFUNC);
    group = XLALH5GroupOpen(file, "cbc_waveform_params");
    XLAL_CHECK_FAIL(group, XLAL_EFUNC);

    /* get the types of all parameter keys */
    types = XLALCreateDict();
    for (size_t i = 0; i < injseq->length; ++i)
        XLALDictForeach(injseq->data[i], XLALSimInspiralDictInsertParameterType, types);

    /* loop over parameters */
    XLALDictIterInit(&iter, types);
    while ((entry = XLALDictIterNext(&iter))) {
        const char *key = XLALDictEntryGetKey(entry);
        const char *new = translate_key(key, 1);
        const LALValue *value = XLALDictEntryGetValue(entry);
        LALTYPECODE type = XLALValueGetINT4(value);
        LALH5Dataset *dset;
        double scale;
        switch (type) {
        case LAL_CHAR_TYPE_CODE:
            strvec = XLALCreateEmptyStringVector(injseq->length);
            for (size_t i = 0; i < injseq->length; ++i)
                strvec->data[i] = XLALStringDuplicate(XLALDictContains(injseq->data[i], key) ? XLALDictLookupStringValue(injseq->data[i], key) : "");
            dset = XLALH5DatasetAllocStringVector(group, new, strvec);
            XLAL_CHECK_FAIL(dset, XLAL_EFUNC);
            XLALDestroyStringVector(strvec);
            strvec = NULL;
            break;
        case LAL_I4_TYPE_CODE:
            intvec = XLALCreateINT4Vector(injseq->length);
            for (size_t i = 0; i < injseq->length; ++i)
                intvec->data[i] = XLALDictLookupINT4Value(injseq->data[i], key);
            dset = XLALH5DatasetAllocINT4Vector(group, new, intvec);
            XLAL_CHECK_FAIL(dset, XLAL_EFUNC);
            XLALDestroyINT4Vector(intvec);
            intvec = NULL;
            break;
        case LAL_D_TYPE_CODE:
            dblvec = XLALCreateREAL8Vector(injseq->length);
            scale = si_scale_factor(key);
            for (size_t i = 0; i < injseq->length; ++i)
                dblvec->data[i] = XLALDictContains(injseq->data[i], key) ? XLALDictLookupREAL8Value(injseq->data[i], key) / scale : NAN;
            dset = XLALH5DatasetAllocREAL8Vector(group, new, dblvec);
            XLAL_CHECK_FAIL(dset, XLAL_EFUNC);
            XLALDestroyREAL8Vector(dblvec);
            dblvec = NULL;
            break;
        default:
           XLAL_ERROR_FAIL(XLAL_ETYPE, "Unsupported data type for dataset \"%s\"", key);
           break;
        }
        XLALH5DatasetFree(dset);
        dset = NULL;
    }

    XLALDestroyDict(types);
    XLALH5FileClose(group);
    XLALH5FileClose(file);
    return XLAL_SUCCESS;

XLAL_FAIL:
    XLALDestroyStringVector(strvec);
    XLALDestroyINT4Vector(intvec);
    XLALDestroyREAL8Vector(dblvec);
    XLALDestroyDict(types);
    XLALH5FileClose(group);
    XLALH5FileClose(file);
    return XLAL_FAILURE;
}

LIGOTimeGPS * XLALSimInspiralInjectionEndTime(LIGOTimeGPS *epoch, LALDict *injparams)
{
    REAL8 mass1;
    REAL8 mass2;
    REAL8 spin1z;
    REAL8 spin2z;
    REAL8 spin_final;
    REAL8 t_ring;
    REAL8 t_co_gps;
    REAL8 t_co_gps_add = 0.0;

    XLAL_CHECK_NULL(epoch, XLAL_EFAULT);
    XLAL_CHECK_NULL(injparams, XLAL_EFAULT);

    t_co_gps = XLALDictLookupREAL8Value(injparams, "t_co_gps");
    XLAL_CHECK_VAL(0, !XLAL_IS_REAL8_FAIL_NAN(t_co_gps), XLAL_EFUNC);

    if (XLALDictContains(injparams, "t_co_gps_add")) {
        t_co_gps_add = XLALDictLookupREAL8Value(injparams, "t_co_gps_add");
        XLAL_CHECK_VAL(0, !XLAL_IS_REAL8_FAIL_NAN(t_co_gps_add), XLAL_EFUNC);
    }

    XLALGPSSetREAL8(epoch, t_co_gps);
    XLALGPSAdd(epoch, t_co_gps_add);

    mass1 = XLALDictLookupREAL8Value(injparams, "mass1");
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(mass1), XLAL_EFUNC);

    mass2 = XLALDictLookupREAL8Value(injparams, "mass2");
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(mass2), XLAL_EFUNC);

    spin1z = XLALDictLookupREAL8Value(injparams, "spin1z");
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(spin1z), XLAL_EFUNC);

    spin2z = XLALDictLookupREAL8Value(injparams, "spin2z");
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(spin2z), XLAL_EFUNC);

    spin_final = XLALSimInspiralFinalBlackHoleSpinBound(spin1z, spin2z);
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(spin_final), XLAL_EFUNC);

    t_ring = XLALSimInspiralRingdownTimeBound(mass1 + mass2, spin_final);
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(t_ring), XLAL_EFUNC);
 
    return XLALGPSAdd(epoch, t_ring);
}

LIGOTimeGPS * XLALSimInspiralInjectionStartTime(LIGOTimeGPS *epoch, LALDict *injparams)
{
    REAL8 f22_start;
    REAL8 mass1;
    REAL8 mass2;
    REAL8 spin1z;
    REAL8 spin2z;
    REAL8 t_chirp;
    REAL8 t_merge;

    epoch = XLALSimInspiralInjectionEndTime(epoch, injparams);
    XLAL_CHECK_NULL(epoch, XLAL_EFUNC);

    f22_start = XLALDictLookupREAL8Value(injparams, "f22_start");
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(f22_start), XLAL_EFUNC);

    mass1 = XLALDictLookupREAL8Value(injparams, "mass1");
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(mass1), XLAL_EFUNC);

    mass2 = XLALDictLookupREAL8Value(injparams, "mass2");
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(mass2), XLAL_EFUNC);

    spin1z = XLALDictLookupREAL8Value(injparams, "spin1z");
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(spin1z), XLAL_EFUNC);

    spin2z = XLALDictLookupREAL8Value(injparams, "spin2z");
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(spin2z), XLAL_EFUNC);

    t_chirp = XLALSimInspiralChirpTimeBound(f22_start, mass1, mass2, spin1z, spin2z);
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(t_chirp), XLAL_EFUNC);

    t_merge = XLALSimInspiralMergeTimeBound(mass1, mass2);
    XLAL_CHECK_NULL(!XLAL_IS_REAL8_FAIL_NAN(t_merge), XLAL_EFUNC);

    return XLALGPSAdd(epoch, -(t_chirp + t_merge));
}

static int XLALSimInspiralParamsCompareEndTime(void *thunk UNUSED, const void *a, const void *b)
{
    /* note: discarding const qualifier is safe in this context */
    LALDict *dict1 = *(LALDict **)(uintptr_t)a;
    LALDict *dict2 = *(LALDict **)(uintptr_t)b;
    LIGOTimeGPS t1;
    LIGOTimeGPS t2;
    XLALSimInspiralInjectionEndTime(&t1, dict1);
    XLALSimInspiralInjectionEndTime(&t2, dict2);
    return XLALGPSCmp(&t1, &t2);
}

static int XLALSimInspiralParamsCompareEndTimeToGPSTime(void *thunk UNUSED, const void *a, const void *b)
{
    /* note: discarding const qualifier is safe in this context */
    LALDict *dict = *(LALDict **)(uintptr_t)a;
    LIGOTimeGPS *t2 = (LIGOTimeGPS *)(uintptr_t)b;
    LIGOTimeGPS t1;
    XLALSimInspiralInjectionEndTime(&t1, dict);
    return XLALGPSCmp(&t1, t2);
}

static int XLALSimInspiralParamsCompareStartTime(void *thunk UNUSED, const void *a, const void *b)
{
    /* note: discarding const qualifier is safe in this context */
    LALDict *dict1 = *(LALDict **)(uintptr_t)a;
    LALDict *dict2 = *(LALDict **)(uintptr_t)b;
    LIGOTimeGPS t1;
    LIGOTimeGPS t2;
    XLALSimInspiralInjectionStartTime(&t1, dict1);
    XLALSimInspiralInjectionStartTime(&t2, dict2);
    return XLALGPSCmp(&t1, &t2);
}

static int XLALSimInspiralParamsCompareStartTimeToGPSTime(void *thunk UNUSED, const void *a, const void *b)
{
    /* note: discarding const qualifier is safe in this context */
    LALDict *dict = *(LALDict **)(uintptr_t)a;
    LIGOTimeGPS *t2 = (LIGOTimeGPS *)(uintptr_t)b;
    LIGOTimeGPS t1;
    XLALSimInspiralInjectionStartTime(&t1, dict);
    return XLALGPSCmp(&t1, t2);
}

int XLALSimInspiralInjectionSequenceIsEndTimeOrdered(LALDictSequence *injseq)
{
    int errnum;
    int retval;
    XLAL_CHECK(injseq, XLAL_EFAULT);
    errnum = XLALClearErrno(); /* clear xlalErrno and preserve value */
    retval = XLALIsSorted(injseq->data, injseq->length, sizeof(*injseq->data), NULL, XLALSimInspiralParamsCompareEndTime);
    if (retval < 0 || XLALGetBaseErrno()) /* something failed */
        XLAL_ERROR(XLAL_EFUNC);
    XLALSetErrno(errnum); /* restore previous value of xlalErrno */
    return retval;
}

int XLALSimInspiralInjectionSequenceIsStartTimeOrdered(LALDictSequence *injseq)
{
    int errnum;
    int retval;
    XLAL_CHECK(injseq, XLAL_EFAULT);
    errnum = XLALClearErrno(); /* clear xlalErrno and preserve value */
    retval = XLALIsSorted(injseq->data, injseq->length, sizeof(*injseq->data), NULL, XLALSimInspiralParamsCompareStartTime);
    if (retval < 0 || XLALGetBaseErrno()) /* something failed */
        XLAL_ERROR(XLAL_EFUNC);
    XLALSetErrno(errnum); /* restore previous value of xlalErrno */
    return retval;
}

int XLALSimInspiralInjectionSequenceOrderByEndTime(LALDictSequence *injseq)
{
    int retval;
    /* it is likely that it is already time ordered, so check first */
    retval = XLALSimInspiralInjectionSequenceIsEndTimeOrdered(injseq);
    XLAL_CHECK(retval >= 0, XLAL_EFUNC);
    if (!retval) { /* sort it */
        int errnum;
        errnum = XLALClearErrno(); /* clear xlalErrno and preserve value */
        retval = XLALMergeSort(injseq->data, injseq->length, sizeof(*injseq->data), NULL, XLALSimInspiralParamsCompareEndTime);
        if (retval < 0 || XLALGetBaseErrno()) /* something failed */
            XLAL_ERROR(XLAL_EFUNC);
        XLALSetErrno(errnum); /* restore previous value of xlalErrno */
    }
    return 0;
}

int XLALSimInspiralInjectionSequenceOrderByStartTime(LALDictSequence *injseq)
{
    int retval;
    /* it is likely that it is already time ordered, so check first */
    retval = XLALSimInspiralInjectionSequenceIsStartTimeOrdered(injseq);
    XLAL_CHECK(retval >= 0, XLAL_EFUNC);
    if (!retval) { /* sort it */
        int errnum;
        errnum = XLALClearErrno(); /* clear xlalErrno and preserve value */
        retval = XLALMergeSort(injseq->data, injseq->length, sizeof(*injseq->data), NULL, XLALSimInspiralParamsCompareStartTime);
        if (retval < 0 || XLALGetBaseErrno()) /* something failed */
            XLAL_ERROR(XLAL_EFUNC);
        XLALSetErrno(errnum); /* restore previous value of xlalErrno */
    }
    return 0;
}

LALDictSequence * XLALSimInspiralInjectionSequenceInInterval(const LALDictSequence *injseq, const LIGOTimeGPS *start, const LIGOTimeGPS *end)
{
    LALDictSequence *new = NULL;
    LALDictSequence *tmp;
    ssize_t i;
    int retval;

    XLAL_CHECK_NULL(injseq, XLAL_EFAULT);

    /* copy sequence */
    new = XLALCopyDictSequence(injseq);
    XLAL_CHECK_NULL(new, XLAL_EFUNC);

    /* keep only injections where injection end time > start time */

    /* sort by injection end time */
    retval = XLALSimInspiralInjectionSequenceOrderByEndTime(new);
    XLAL_CHECK_FAIL(retval == XLAL_SUCCESS, XLAL_EFUNC);

    /* find index of first injection ending after start of interval */
    i = XLALSearchSorted(start, new->data, new->length, sizeof(*new->data), NULL, XLALSimInspiralParamsCompareEndTimeToGPSTime, +1);
    XLAL_CHECK_FAIL(i >= 0, XLAL_EFUNC);
    tmp = XLALResizeDictSequence(new, i, new->length - i);
    XLAL_CHECK_FAIL(tmp, XLAL_EFUNC);
    new = tmp;

    /* sort by injection start time */
    retval = XLALSimInspiralInjectionSequenceOrderByStartTime(new);
    XLAL_CHECK_FAIL(retval == XLAL_SUCCESS, XLAL_EFUNC);

    /* find index of last injection starting before end of interval */
    i = XLALSearchSorted(end, new->data, new->length, sizeof(*new->data), NULL, XLALSimInspiralParamsCompareStartTimeToGPSTime, -1);
    XLAL_CHECK_FAIL(i >= 0, XLAL_EFUNC);
    tmp = XLALResizeDictSequence(new, 0, i);
    XLAL_CHECK_FAIL(tmp, XLAL_EFUNC);
    new = tmp;

    /* note: returned injection sequence is ordered by injection start time */
    return new;

XLAL_FAIL:
    XLALDestroyDictSequence(new);
    return NULL;
}

int XLALSimInspiralInjectionTDWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LALDict *injparams, REAL8 deltaT)
{
    LALSimInspiralGenerator *generator = NULL;
    LALDict *wfmparams;
    Approximant approx;
    char *approximant;
    int err;
    int ret = XLAL_FAILURE;

    XLAL_CHECK(hplus, XLAL_EFAULT);
    XLAL_CHECK(hcross, XLAL_EFAULT);
    XLAL_CHECK(injparams, XLAL_EFAULT);
    XLAL_CHECK(*hplus == NULL, XLAL_EINVAL);
    XLAL_CHECK(*hcross == NULL, XLAL_EINVAL);

    wfmparams = XLALDictDuplicate(injparams);
    XLAL_CHECK(wfmparams, XLAL_EFUNC);

    /* remove non-waveform parameters */
    err = XLALDictContains(wfmparams, "ra") ? XLALDictRemove(wfmparams, "ra"): 0;
    XLAL_CHECK(err == 0, XLAL_EFUNC);
    err = XLALDictContains(wfmparams, "dec") ? XLALDictRemove(wfmparams, "dec"): 0;
    XLAL_CHECK(err == 0, XLAL_EFUNC);
    err = XLALDictContains(wfmparams, "polarization") ? XLALDictRemove(wfmparams, "polarization"): 0;
    XLAL_CHECK(err == 0, XLAL_EFUNC);
    err = XLALDictContains(wfmparams, "t_co_gps") ? XLALDictRemove(wfmparams, "t_co_gps"): 0;
    XLAL_CHECK(err == 0, XLAL_EFUNC);
    err = XLALDictContains(wfmparams, "t_co_gps_add") ? XLALDictRemove(wfmparams, "t_co_gps_add"): 0;
    XLAL_CHECK(err == 0, XLAL_EFUNC);

    /* pop approximant string and convert to Approximant */
    approximant = XLALDictPopStringValue(wfmparams, "approximant");
    XLAL_CHECK_FAIL(approximant, XLAL_EFUNC);
    approx = XLALSimInspiralGetApproximantFromString(approximant);
    LALFree(approximant);

    /* insert deltaT */
    ret = XLALSimInspiralWaveformParamsInsertDeltaT(wfmparams, deltaT);
    XLAL_CHECK_FAIL(!(ret < 0), XLAL_EFUNC);

    /* create generator FIXME? allow parameters to be added here? */
    generator = XLALSimInspiralChooseGenerator(approx, NULL);
    XLAL_CHECK_FAIL(generator, XLAL_EFUNC);

    /* add conditioning for approximant */
    ret = XLALSimInspiralGeneratorAddConditioningForApproximant(generator, approx);
    XLAL_CHECK_FAIL(!(ret < 0), XLAL_EFUNC);

    /* TODO translate ModeString into ModeArray */

    /* generate the waveform */
    ret = XLALSimInspiralGenerateTDWaveform(hplus, hcross, wfmparams, generator);
    XLAL_CHECK_FAIL(!(ret < 0), XLAL_EFUNC);

XLAL_FAIL:
    XLALDestroySimInspiralGenerator(generator);
    XLALDestroyDict(wfmparams);
    return ret;
}

REAL8TimeSeries * XLALSimInspiralInjectionStrain(LALDict *injparams, REAL8 deltaT, const LALDetector *detector)
{
    REAL8TimeSeries *hplus = NULL;
    REAL8TimeSeries *hcross = NULL;
    REAL8TimeSeries *strain = NULL;
    REAL8 t_co_gps_add = 0.0;
    REAL8 t_co_gps;
    REAL8 ra;
    REAL8 dec;
    REAL8 polarization;
    int ret;

    t_co_gps = XLALDictLookupREAL8Value(injparams, "t_co_gps");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(t_co_gps), XLAL_EFUNC);

    if (XLALDictContains(injparams, "t_co_gps_add")) {
        t_co_gps = XLALDictLookupREAL8Value(injparams, "t_co_gps_add");
        XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(t_co_gps_add), XLAL_EFUNC);
    }

    ra = XLALDictLookupREAL8Value(injparams, "ra");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(ra), XLAL_EFUNC);

    dec = XLALDictLookupREAL8Value(injparams, "dec");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(dec), XLAL_EFUNC);

    polarization = XLALDictLookupREAL8Value(injparams, "polarization");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(polarization), XLAL_EFUNC);

    ret = XLALSimInspiralInjectionTDWaveform(&hplus, &hcross, injparams, deltaT);
    XLAL_CHECK_FAIL(!(ret < 0), XLAL_EFUNC);

    XLALGPSAdd(&hplus->epoch, t_co_gps);
    XLALGPSAdd(&hplus->epoch, t_co_gps_add);
    XLALGPSAdd(&hcross->epoch, t_co_gps);
    XLALGPSAdd(&hcross->epoch, t_co_gps_add);

    strain = XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, ra, dec, polarization, detector);
    XLAL_CHECK_FAIL(strain, XLAL_EFUNC);

XLAL_FAIL:
    XLALDestroyREAL8TimeSeries(hcross);
    XLALDestroyREAL8TimeSeries(hplus);
    return strain;
}
