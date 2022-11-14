#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDetectors.h>
#include <lal/LALString.h>
#include <lal/LALDict.h>
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


void XLALSimInspiralDestroyInjectionSequence(LALSimInspiralInjectionSequence *sequence)
{
    if (sequence) {
        if (sequence->data) {
            for (size_t i = 0; i < sequence->length; ++i)
                XLALDestroyDict(sequence->data[i]);
            LALFree(sequence->data);
        }
        LALFree(sequence);
    }
    return;
}

LALSimInspiralInjectionSequence * XLALSimInspiralCreateInjectionSequence(size_t length)
{
    LALSimInspiralInjectionSequence *sequence;
    sequence = LALMalloc(sizeof(*sequence));
    XLAL_CHECK_NULL(sequence, XLAL_ENOMEM);
    sequence->length = length;
    sequence->data = LALCalloc(length, sizeof(*sequence->data));
    if (sequence->data == NULL) {
        LALFree(sequence);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    for (size_t i = 0; i < length; ++i) {
        sequence->data[i] = XLALCreateDict();
        if (sequence->data[i] == NULL) {
            XLALSimInspiralDestroyInjectionSequence(sequence);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }
    }
    return sequence;
}

LALSimInspiralInjectionSequence * XLALSimInspiralCutInjectionSequence(const LALSimInspiralInjectionSequence *sequence, size_t first, size_t length)
{
    LALSimInspiralInjectionSequence *new;
    XLAL_CHECK_NULL(sequence, XLAL_EFAULT);
    XLAL_CHECK_NULL(first + length <= sequence->length, XLAL_EBADLEN);
    new = XLALSimInspiralCreateInjectionSequence(length);
    XLAL_CHECK_NULL(new, XLAL_EFUNC);
    for (size_t i = 0; i < length; ++i) {
        int retval = XLALDictUpdate(new->data[i], sequence->data[i + first]);
        if (retval < 0) {
            XLALSimInspiralDestroyInjectionSequence(new);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }
    }
    return new;
}

LALSimInspiralInjectionSequence * XLALSimInspiralCopyInjectionSequence(const LALSimInspiralInjectionSequence *sequence)
{
    XLAL_CHECK_NULL(sequence, XLAL_EFAULT);
    return XLALSimInspiralCutInjectionSequence(sequence, 0, sequence->length);
}

void XLALSimInspiralShiftInjectionSequence(LALSimInspiralInjectionSequence *sequence, ssize_t count)
{
    size_t abscount = count > 0 ? count : -count;

    XLAL_CHECK_VOID(sequence, XLAL_EFAULT);
    XLAL_CHECK_VOID(sequence->data || sequence->length == 0, XLAL_EINVAL);

    if (!count)
        return;

    /* shift memory if abs(count) < sequence->length */
    if (sequence->length > abscount) {
        char *a = (char *)sequence->data + (count < 0 ? 0 : count);
        char *b = (char *)sequence->data - (count < 0 ? count : 0);
        memmove(a, b, (sequence->length - abscount) * sizeof(*sequence->data));
    }

    /* clear unshifted memory */
    if (abscount >= sequence->length) {
        for (size_t i = 0; i < sequence->length; ++i)
            XLALClearDict(sequence->data[i]);
    } else if (count > 0) {
        for (size_t i = 0; i < abscount; ++i)
            XLALClearDict(sequence->data[i]);
    } else {
        for (size_t i = sequence->length - abscount; i < abscount; ++i)
            XLALClearDict(sequence->data[i]);
    }

    return;
}

LALSimInspiralInjectionSequence * XLALSimInspiralResizeInjectionSequence(LALSimInspiralInjectionSequence *sequence, ssize_t first, size_t length)
{
    LALDict **data;
    XLAL_CHECK_NULL(sequence, XLAL_EFAULT);
    XLAL_CHECK_NULL(sequence->data || sequence->length == 0, XLAL_EINVAL);

    if (length > sequence->length) { /* need to increase memory */
        data = XLALRealloc(sequence->data, length * sizeof(*sequence->data));
        XLAL_CHECK_NULL(data, XLAL_EFUNC);
        for (size_t i = sequence->length; i < length - sequence->length; ++i)
            data[i] = XLALCreateDict();
        sequence->data = data;
        sequence->length = length;
        XLALSimInspiralShiftInjectionSequence(sequence, -first);
    } else if (length > 0) { /* need to decrease memory */
        XLALSimInspiralShiftInjectionSequence(sequence, -first);
        for (size_t i = length; i < sequence->length; ++i)
            XLALDestroyDict(sequence->data[i]);
        data = XLALRealloc(sequence->data, length * sizeof(*sequence->data));
        XLAL_CHECK_NULL(data, XLAL_EFUNC);
        sequence->data = data;
        sequence->length = length;
    } else { /* length == 0: need to release all memory */
        for (size_t i = 0; i < sequence->length; ++i)
            XLALDestroyDict(sequence->data[i]);
        XLALFree(sequence->data);
        sequence->data = NULL;
        sequence->length = 0;
    }

    return sequence;
}

/* translates from the injection file parameter names
 * https://git.ligo.org/waveforms/o4-injection-file/-/blob/main/ESSENTIAL_SPEC.md * to the waveform interface parameter names
 * https://git.ligo.org/waveforms/new-waveforms-interface/-/blob/master/parameter-names.md
 */
static const char * translate_key(const char *key, int inverse)
{
    const char *map[][2] = {
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

LALSimInspiralInjectionSequence * XLALSimInspiralInjectionSequenceFromH5File(const char *fname)
{
    LALSimInspiralInjectionSequence *sequence = NULL;
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
        if (sequence == NULL) {
            sequence = XLALSimInspiralCreateInjectionSequence(ninj);
            XLAL_CHECK_FAIL(sequence, XLAL_EFUNC);
        } else if (ninj != sequence->length)
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
                   retval = XLALDictInsertStringValue(sequence->data[i], key, strvec->data[i]);
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
               retval = XLALDictInsertINT4Value(sequence->data[i], key, intvec->data[i]);
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
               retval = XLALDictInsertREAL8Value(sequence->data[i], key, scale * dblvec->data[i]);
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

    return sequence;

XLAL_FAIL:

    XLALFree(name);
    XLALH5DatasetFree(dset);
    XLALH5FileClose(group);
    XLALH5FileClose(file);
    XLALDestroyREAL8Vector(dblvec);
    XLALDestroyINT4Vector(intvec);
    XLALDestroyStringVector(strvec);
    XLALSimInspiralDestroyInjectionSequence(sequence);
    return NULL;
}


static void XLALSimInspiralDictInsertParameterType(char *key, LALValue *value, void *thunk)
{
    LALDict *types = (LALDict *)thunk;
    XLALDictInsertINT4Value(types, key, XLALValueGetType(value));
    return;
}

int XLALSimInspiralInjectionSequenceToH5File(const LALSimInspiralInjectionSequence *sequence, const char *fname)
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
    for (size_t i = 0; i < sequence->length; ++i)
        XLALDictForeach(sequence->data[i], XLALSimInspiralDictInsertParameterType, types);

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
            strvec = XLALCreateEmptyStringVector(sequence->length);
            for (size_t i = 0; i < sequence->length; ++i)
                strvec->data[i] = XLALStringDuplicate(XLALDictContains(sequence->data[i], key) ? XLALDictLookupStringValue(sequence->data[i], key) : "");
            dset = XLALH5DatasetAllocStringVector(group, new, strvec);
            XLAL_CHECK_FAIL(dset, XLAL_EFUNC);
            XLALDestroyStringVector(strvec);
            strvec = NULL;
            break;
        case LAL_I4_TYPE_CODE:
            intvec = XLALCreateINT4Vector(sequence->length);
            for (size_t i = 0; i < sequence->length; ++i)
                intvec->data[i] = XLALDictLookupINT4Value(sequence->data[i], key);
            dset = XLALH5DatasetAllocINT4Vector(group, new, intvec);
            XLAL_CHECK_FAIL(dset, XLAL_EFUNC);
            XLALDestroyINT4Vector(intvec);
            intvec = NULL;
            break;
        case LAL_D_TYPE_CODE:
            dblvec = XLALCreateREAL8Vector(sequence->length);
            scale = si_scale_factor(key);
            for (size_t i = 0; i < sequence->length; ++i)
                dblvec->data[i] = XLALDictContains(sequence->data[i], key) ? XLALDictLookupREAL8Value(sequence->data[i], key) / scale : nan(NULL);
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

int XLALSimInspiralInjectionSequenceIsEndTimeOrdered(LALSimInspiralInjectionSequence *sequence)
{
    int errnum;
    int retval;
    XLAL_CHECK(sequence, XLAL_EFAULT);
    errnum = XLALClearErrno(); /* clear xlalErrno and preserve value */
    retval = XLALIsSorted(sequence->data, sequence->length, sizeof(*sequence->data), NULL, XLALSimInspiralParamsCompareEndTime);
    if (retval < 0 || XLALGetBaseErrno()) /* something failed */
        XLAL_ERROR(XLAL_EFUNC);
    XLALSetErrno(errnum); /* restore previous value of xlalErrno */
    return retval;
}

int XLALSimInspiralInjectionSequenceIsStartTimeOrdered(LALSimInspiralInjectionSequence *sequence)
{
    int errnum;
    int retval;
    XLAL_CHECK(sequence, XLAL_EFAULT);
    errnum = XLALClearErrno(); /* clear xlalErrno and preserve value */
    retval = XLALIsSorted(sequence->data, sequence->length, sizeof(*sequence->data), NULL, XLALSimInspiralParamsCompareStartTime);
    if (retval < 0 || XLALGetBaseErrno()) /* something failed */
        XLAL_ERROR(XLAL_EFUNC);
    XLALSetErrno(errnum); /* restore previous value of xlalErrno */
    return retval;
}

int XLALSimInspiralInjectionSequenceOrderByEndTime(LALSimInspiralInjectionSequence *sequence)
{
    int retval;
    /* it is likely that it is already time ordered, so check first */
    retval = XLALSimInspiralInjectionSequenceIsEndTimeOrdered(sequence);
    XLAL_CHECK(retval >= 0, XLAL_EFUNC);
    if (!retval) { /* sort it */
        int errnum;
        errnum = XLALClearErrno(); /* clear xlalErrno and preserve value */
        retval = XLALMergeSort(sequence->data, sequence->length, sizeof(*sequence->data), NULL, XLALSimInspiralParamsCompareEndTime);
        if (retval < 0 || XLALGetBaseErrno()) /* something failed */
            XLAL_ERROR(XLAL_EFUNC);
        XLALSetErrno(errnum); /* restore previous value of xlalErrno */
    }
    return 0;
}

int XLALSimInspiralInjectionSequenceOrderByStartTime(LALSimInspiralInjectionSequence *sequence)
{
    int retval;
    /* it is likely that it is already time ordered, so check first */
    retval = XLALSimInspiralInjectionSequenceIsStartTimeOrdered(sequence);
    XLAL_CHECK(retval >= 0, XLAL_EFUNC);
    if (!retval) { /* sort it */
        int errnum;
        errnum = XLALClearErrno(); /* clear xlalErrno and preserve value */
        retval = XLALMergeSort(sequence->data, sequence->length, sizeof(*sequence->data), NULL, XLALSimInspiralParamsCompareStartTime);
        if (retval < 0 || XLALGetBaseErrno()) /* something failed */
            XLAL_ERROR(XLAL_EFUNC);
        XLALSetErrno(errnum); /* restore previous value of xlalErrno */
    }
    return 0;
}

LALSimInspiralInjectionSequence * XLALSimInspiralInjectionSequenceInInterval(const LALSimInspiralInjectionSequence *sequence, const LIGOTimeGPS *start, const LIGOTimeGPS *end)
{
    LALSimInspiralInjectionSequence *new = NULL;
    LALSimInspiralInjectionSequence *tmp;
    ssize_t i;
    int retval;

    XLAL_CHECK_NULL(sequence, XLAL_EFAULT);

    /* copy sequence */
    new = XLALSimInspiralCopyInjectionSequence(sequence);
    XLAL_CHECK_NULL(new, XLAL_EFUNC);

    /* keep only injections where injection end time > start time */

    /* sort by injection end time */
    retval = XLALSimInspiralInjectionSequenceOrderByEndTime(new);
    XLAL_CHECK_FAIL(retval == XLAL_SUCCESS, XLAL_EFUNC);

    /* find index of first injection ending after start of interval */
    i = XLALSearchSorted(start, new->data, new->length, sizeof(*new->data), NULL, XLALSimInspiralParamsCompareEndTimeToGPSTime, +1);
    XLAL_CHECK_FAIL(i >= 0, XLAL_EFUNC);
    tmp = XLALSimInspiralResizeInjectionSequence(new, i, new->length - i);
    XLAL_CHECK_FAIL(tmp, XLAL_EFUNC);
    new = tmp;

    /* sort by injection start time */
    retval = XLALSimInspiralInjectionSequenceOrderByStartTime(new);
    XLAL_CHECK_FAIL(retval == XLAL_SUCCESS, XLAL_EFUNC);

    /* find index of last injection starting before end of interval */
    i = XLALSearchSorted(end, new->data, new->length, sizeof(*new->data), NULL, XLALSimInspiralParamsCompareStartTimeToGPSTime, -1);
    XLAL_CHECK_FAIL(i >= 0, XLAL_EFUNC);
    tmp = XLALSimInspiralResizeInjectionSequence(new, 0, i);
    XLAL_CHECK_FAIL(tmp, XLAL_EFUNC);
    new = tmp;

    /* note: returned injection sequence is ordered by injection start time */
    return new;

XLAL_FAIL:
    XLALSimInspiralDestroyInjectionSequence(new);
    return NULL;
}

int XLALSimInspiralInjectionTDWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LALDict *injparams, REAL8 deltaT)
{
    /* TODO: change to new waveform interface */

    LALDict *wfmparams;
    Approximant approx;
    REAL8 mass1;
    REAL8 mass2;
    REAL8 spin1x;
    REAL8 spin1y;
    REAL8 spin1z;
    REAL8 spin2x;
    REAL8 spin2y;
    REAL8 spin2z;
    REAL8 inclination;
    REAL8 phi_ref;
    REAL8 distance;
    REAL8 f22_ref;
    REAL8 f22_start;
    REAL8 longAscNodes;
    REAL8 eccentricity;
    REAL8 meanPerAno;
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

    /* pop required parameters */
    mass1 = XLALDictPopREAL8Value(wfmparams, "mass1");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(mass1), XLAL_EFUNC);
    mass2 = XLALDictPopREAL8Value(wfmparams, "mass2");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(mass2), XLAL_EFUNC);
    spin1x = XLALDictPopREAL8Value(wfmparams, "spin1x");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(spin1x), XLAL_EFUNC);
    spin1y = XLALDictPopREAL8Value(wfmparams, "spin1y");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(spin1y), XLAL_EFUNC);
    spin1z = XLALDictPopREAL8Value(wfmparams, "spin1z");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(spin1z), XLAL_EFUNC);
    spin2x = XLALDictPopREAL8Value(wfmparams, "spin2x");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(spin2x), XLAL_EFUNC);
    spin2y = XLALDictPopREAL8Value(wfmparams, "spin2y");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(spin2y), XLAL_EFUNC);
    spin2z = XLALDictPopREAL8Value(wfmparams, "spin2z");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(spin2z), XLAL_EFUNC);
    distance = XLALDictPopREAL8Value(wfmparams, "distance");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(distance), XLAL_EFUNC);
    inclination = XLALDictPopREAL8Value(wfmparams, "inclination");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(inclination), XLAL_EFUNC);
    phi_ref = XLALDictPopREAL8Value(wfmparams, "phi_ref");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(phi_ref), XLAL_EFUNC);
    f22_ref = XLALDictPopREAL8Value(wfmparams, "f22_ref");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(f22_ref), XLAL_EFUNC);
    f22_start = XLALDictPopREAL8Value(wfmparams, "f22_start");
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(f22_start), XLAL_EFUNC);
    approximant = XLALDictPopStringValue(wfmparams, "approximant");
    XLAL_CHECK_FAIL(approximant, XLAL_EFUNC);
    approx = XLALSimInspiralGetApproximantFromString(approximant);
    LALFree(approximant);

    /* pop optional parameters */
    longAscNodes = XLALDictContains(wfmparams, "longAscNodes") ? XLALDictPopREAL8Value(wfmparams, "longAscNodes"): 0.0;
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(longAscNodes), XLAL_EFUNC);
    eccentricity = XLALDictContains(wfmparams, "eccentricity") ? XLALDictPopREAL8Value(wfmparams, "eccentricity"): 0.0;
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(eccentricity), XLAL_EFUNC);
    meanPerAno = XLALDictContains(wfmparams, "meanPerAno") ? XLALDictPopREAL8Value(wfmparams, "meanPerAno"): 0.0;
    XLAL_CHECK_FAIL(!XLAL_IS_REAL8_FAIL_NAN(meanPerAno), XLAL_EFUNC);
  
    /* TODO translate ModeString into ModeArray */

    /* generate the waveform */
    ret = XLALSimInspiralTD(hplus, hcross, mass1, mass2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, distance, inclination, phi_ref, longAscNodes, eccentricity, meanPerAno, deltaT, f22_start, f22_ref, wfmparams, approx);
    XLAL_CHECK_FAIL(!(ret < 0), XLAL_EFUNC);

XLAL_FAIL:
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
