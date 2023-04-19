#include "config.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#ifndef HAVE_PYTHON3_EMBED

#include <lal/LALStdlib.h>
#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include "LALSimInspiralGenerator_private.h"

static int initialize(UNUSED LALSimInspiralGenerator *myself, UNUSED LALDict *params)
{
    XLAL_ERROR(XLAL_EINVAL, "Python generator not implemented");
}

const LALSimInspiralGenerator lalPythonGeneratorTemplate = {
    .name = "ExternalPython",
    .initialize = initialize,
    .finalize = NULL,
    .generate_td_waveform = NULL,
    .generate_fd_waveform = NULL,
    .generate_td_modes = NULL,
    .generate_fd_modes = NULL,
    .internal_data = NULL
};

#else /* HAVE_PYTHON_EMBED */

/* Python.h requires the assert macro, so don't poison it */
#undef LAL_STRICT_DEFS_ENABLED
#include <lal/LALStddef.h>

#include <lal/LALStdlib.h>
#include <lal/LALDict.h>
#include <lal/LALSimInspiral.h>
#include "LALSimInspiralGenerator_private.h"

#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALString.h>
#include <lal/LALValue.h>
#include <lal/Sequence.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimSphHarmSeries.h>
#include "LALSimInspiralWaveformParams_common.c"

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif
#include <numpy/arrayobject.h>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

/* Create laldict where the units of all the parameters in LALSimInspiralWavofrmParams_common.c are inserted as strings */
static LALDict *create_waveform_parameter_units_dict(void)
{

    LALDict *dict;
    dict = XLALCreateDict();
    if (dict == NULL)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    for (size_t i = 0; i < XLAL_NUM_ELEM(lalSimInspiralREAL8WaveformParams); ++i)
        if (XLALDictInsertStringValue(dict, lalSimInspiralREAL8WaveformParams[i].name,
				lalSimInspiralREAL8WaveformParams[i].unit)) {
            XLALDestroyDict(dict);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }

    return dict;
}

/* emit error messages corresponding to python exceptions */
static int python_error(void)
{
    if (PyErr_Occurred()) {
        XLAL_PRINT_ERROR("Python exception raised");
        PyErr_Print();
        return 1;
    }
    return 0;
}


/* ROUTINES TO ENCODE A LALDICT AS A PYTHON DICT */

/* Helper struct containing a python dictionary, an astropy quantity and a laldict containing the SI units of the waveform parameters */
struct add_to_dict_params {
    PyObject *dict;
    PyObject *quantity;
    LALDict *waveform_parameter_units;
};

/* Set a new value for a particular key in a input python dictionary (thunk.dict) */
static void add_to_dict(char *key, LALValue *value, void *thunk)
{
    struct add_to_dict_params *add_to_dict_params = thunk;
    LALDict *units = add_to_dict_params->waveform_parameter_units;
    PyObject *quantity = add_to_dict_params->quantity;
    PyObject *dict = add_to_dict_params->dict;  /* python dict */
    PyObject *args = NULL;
    PyObject *val = NULL; /* python value */
    /* Transform LALValue to python value */
    switch (XLALValueGetType(value)) {
    case LAL_CHAR_TYPE_CODE:
        val = PyUnicode_FromString(XLALValueGetString(value));
        XLAL_CHECK_FAIL(val && !PyErr_Occurred(), XLAL_EFAILED, "Failed to create Python value for parameter %s", key);
        break;
    case LAL_I4_TYPE_CODE:
        val = PyLong_FromLong(XLALValueGetINT4(value));
        XLAL_CHECK_FAIL(val && !PyErr_Occurred(), XLAL_EFAILED, "Failed to create Python value for parameter %s", key);
        break;
    case LAL_D_TYPE_CODE:
        val = PyFloat_FromDouble(XLALValueGetREAL8(value));
        XLAL_CHECK_FAIL(val && !PyErr_Occurred(), XLAL_EFAILED, "Failed to create Python value for parameter %s", key);
        /* Use an astropy quantity. The input LALValue is assumed to be in SI units, those are read from thunk.waveform_parameter_units */
        if (quantity && units) {
            const char *unit = ""; /* dimensionless */
            if (XLALDictContains(units, key))
                unit = XLALDictLookupStringValue(units, key);
            XLAL_CHECK_FAIL(unit, XLAL_EFUNC);
            args = Py_BuildValue("Os", val, unit);
            XLAL_CHECK_FAIL(args, XLAL_EFAILED, "Failed to build args for parameter %s", key);
            Py_CLEAR(val);
            val = PyObject_CallObject(quantity, args);
            XLAL_CHECK_FAIL(val, XLAL_EFAILED, "Failed to create Quantity for parameter %s with unit %s", key, unit);
            Py_CLEAR(args);
        }
        break;
    default:
        XLAL_ERROR_FAIL(XLAL_ETYPE, "Invalid type for parameter %s", key);
        break;
    }
    /* Insert value in python dictionary */
    XLAL_CHECK_FAIL(PyDict_SetItemString(dict, key, val) == 0, XLAL_EFAILED, "Failed to insert parameter %s in Python dict", key);
    Py_DECREF(val);
    return;

XLAL_FAIL:
    python_error();
    Py_XDECREF(args);
    Py_XDECREF(val);
    return;
}

/* Transform a laldictionary into a python dictionary */
static PyObject *PyDict_FromLALDict(LALDict *ldict, PyObject *quantity, LALDict *waveform_parameter_units)
{
    struct add_to_dict_params add_to_dict_params;
    PyObject *dict;
    int errnum;
    XLAL_CHECK_NULL(ldict, XLAL_EFAULT);
    dict = PyDict_New();
    XLAL_CHECK_NULL(dict, XLAL_EFAILED, "Could not create Python dict");
    add_to_dict_params.dict = dict;
    add_to_dict_params.quantity = quantity;
    add_to_dict_params.waveform_parameter_units = waveform_parameter_units;
    errnum = XLALClearErrno(); /* clear xlalErrno and preserve value */
    /* Loop over all the keys in the laldict and insert the values into the python dictionary */
    XLALDictForeach(ldict, add_to_dict, &add_to_dict_params);
    if (XLALGetBaseErrno()) { /* something failed... */
        Py_DECREF(dict);
        XLAL_ERROR_NULL(XLAL_EFUNC);
    }
    XLALSetErrno(errnum); /* restore previous value of xlalErrno */
    return dict;
}


/* ROUTINES TO CONVERT NUMPY ARRAYS AS C ARRAYS */

/* 
 * Store data from PyArrayObject in buf of size bufsz, returns number of bytes.
 * If buf is NULL then return the number of bytes in the array.
 * Routine fails if buf is not NULL and bufsz is not equal to the number of
 * bytes in the array; if buf needs to be allocated, call the routine once
 * with buf == NULL to determine the amount of memory to allocate to buf.
 */
static size_t PyArray_AsBytes(void *buf, size_t bufsz, PyArrayObject *arr)
{
    size_t nbytes;
    XLAL_CHECK(arr, XLAL_EFAULT);
    nbytes = PyArray_NBYTES(arr);
    if (buf) {
        XLAL_CHECK(nbytes == bufsz, XLAL_ESIZE, "Inconsisent number of bytes");
        memcpy(buf, PyArray_DATA(arr), bufsz);
    }
    return nbytes;
}

/* Macro to define helper functions that return the number of bytes needed to allocate a numpy array of a specific type */
#define DEFINE_PYARRAY_ASTYPE(LALTYPE, NPYTYPE) \
    static size_t PyArray_As ## LALTYPE (LALTYPE *data, size_t length, PyObject *o) \
    { \
        PyArrayObject *arr; \
        size_t nbytes; \
        XLAL_CHECK(o, XLAL_EFAULT); \
        XLAL_CHECK(PyArray_Check(o), XLAL_EINVAL, "Python object is not a ndarray"); \
        arr = (PyArrayObject *)PyArray_FROM_OF(o, NPY_ARRAY_IN_ARRAY); \
        if (python_error()) \
            XLAL_ERROR(XLAL_EFAILED); \
        XLAL_CHECK(PyArray_TYPE(arr) == NPYTYPE, XLAL_ETYPE, "Ndarray has wrong dtype"); \
        nbytes = PyArray_AsBytes(data, length * sizeof(*data), arr); \
        if ((ssize_t)(nbytes) < 0) \
            XLAL_ERROR(XLAL_EFUNC); \
        return nbytes / sizeof(LALTYPE); \
    }

/* Only a few of these are used... */
// DEFINE_PYARRAY_ASTYPE(CHAR, NPY_BYTE)
// DEFINE_PYARRAY_ASTYPE(INT2, NPY_INT16)
// DEFINE_PYARRAY_ASTYPE(INT4, NPY_INT32)
// DEFINE_PYARRAY_ASTYPE(INT8, NPY_INT64)
// DEFINE_PYARRAY_ASTYPE(UCHAR, NPY_UBYTE)
// DEFINE_PYARRAY_ASTYPE(UINT2, NPY_UINT16)
// DEFINE_PYARRAY_ASTYPE(UINT4, NPY_UINT32)
// DEFINE_PYARRAY_ASTYPE(UINT8, NPY_UINT64)
// DEFINE_PYARRAY_ASTYPE(REAL4, NPY_FLOAT)
DEFINE_PYARRAY_ASTYPE(REAL8, NPY_DOUBLE)
// DEFINE_PYARRAY_ASTYPE(COMPLEX8, NPY_CFLOAT)
DEFINE_PYARRAY_ASTYPE(COMPLEX16, NPY_CDOUBLE)


/* ROUTINES TO CONVERT NUMPY ARRAYS AS LAL SEQUENCES */

#define DEFINE_PYARRAY_ASSEQTYPE(LALTYPE) \
    static LALTYPE ## Sequence * PyArray_As ## LALTYPE ## Sequence (PyObject *o) \
    { \
        LALTYPE ## Sequence *seq; \
        size_t length; \
        XLAL_CHECK_NULL(o, XLAL_EFAULT); \
        if ((ssize_t)(length = PyArray_As ## LALTYPE(NULL, 0, o)) < 0) \
            XLAL_ERROR_NULL(XLAL_EFUNC); \
        if ((seq = XLALCreate ## LALTYPE ## Sequence(length)) == NULL) \
            XLAL_ERROR_NULL(XLAL_EFUNC); \
        if (PyArray_As ## LALTYPE(seq->data, seq->length, o) != length) { \
            XLALDestroy ## LALTYPE ## Sequence(seq); \
            XLAL_ERROR_NULL(XLAL_ESIZE, "Failed to read ndarray"); \
        } \
        return seq; \
    }

/* Only a few of these are used... */
// DEFINE_PYARRAY_ASSEQTYPE(CHAR)
// DEFINE_PYARRAY_ASSEQTYPE(INT2)
// DEFINE_PYARRAY_ASSEQTYPE(INT4)
// DEFINE_PYARRAY_ASSEQTYPE(INT8)
// DEFINE_PYARRAY_ASSEQTYPE(UCHAR)
// DEFINE_PYARRAY_ASSEQTYPE(UINT2)
// DEFINE_PYARRAY_ASSEQTYPE(UINT4)
// DEFINE_PYARRAY_ASSEQTYPE(UINT8)
// DEFINE_PYARRAY_ASSEQTYPE(REAL4)
DEFINE_PYARRAY_ASSEQTYPE(REAL8)
// DEFINE_PYARRAY_ASSEQTYPE(COMPLEX8)
DEFINE_PYARRAY_ASSEQTYPE(COMPLEX16)


/* ROUTINE TO CONVERT ASTROPY TIME AND UNIT AS LAL */

static int AstropyTime_AsLIGOTimeGPS(LIGOTimeGPS *lt, PyObject *at)
{
    PyObject *gpstime = NULL;
    PyObject *result = NULL;
    double sec, nan;

    XLAL_CHECK(lt, XLAL_EFAULT);
    XLAL_CHECK(at, XLAL_EFAULT);

    /* get the gps time at full precision */
    gpstime = PyObject_CallMethod(at, "to_value", "ss", "gps", "decimal");
    XLAL_CHECK_FAIL(gpstime, XLAL_EFAILED, "Could not call method .to_value()");

    /* break it into seconds and nanoseconds */
    result = PyObject_CallMethod(gpstime, "__divmod__", "i", 1);
    XLAL_CHECK_FAIL(result, XLAL_EFAILED, "Could not call method .__divmod__(1)");
    Py_CLEAR(gpstime);

    /* result must be a tuple of two objects */
    XLAL_CHECK_FAIL(PyTuple_Check(result), XLAL_EFAILED, "Invalid result of method .__divmod__(1)");
    XLAL_CHECK_FAIL(PyTuple_Size(result) == 2, XLAL_EFAILED, "Invalid result of method .__divmod__(1)");
    sec = PyFloat_AsDouble(PyTuple_GetItem(result, 0));
    if (PyErr_Occurred())
        XLAL_ERROR_FAIL(XLAL_EFAILED, "Could not convert GPS seconds to double");
    nan = PyFloat_AsDouble(PyTuple_GetItem(result, 1));
    if (PyErr_Occurred())
        XLAL_ERROR_FAIL(XLAL_EFAILED, "Could not convert GPS nanoseconds to double");
    Py_CLEAR(result);

    XLALGPSSet(lt, sec, 1e9 * nan);
    return 0;

XLAL_FAIL:
    python_error();
    Py_XDECREF(gpstime);
    Py_XDECREF(result);
    return (int)XLAL_FAILURE;
}

static double AstropyUnit_AsLALUnit(LALUnit *lu, PyObject *au)
{
    PyObject *si = NULL;
    PyObject *scale = NULL;
    PyObject *bases = NULL;
    PyObject *powers = NULL;
    double scale_factor;

    XLAL_CHECK_REAL8(lu, XLAL_EFAULT);
    XLAL_CHECK_REAL8(au, XLAL_EFAULT);

    /* Get the astropy value in SI units */
    si = PyObject_GetAttrString(au, "si");
    XLAL_CHECK_FAIL(si, XLAL_EFAILED, "Could not get attribute .si");

    /* Get the 3 attributes of the CompositeUnit astropy class (scale, bases, powers) */
    scale = PyObject_GetAttrString(si, "scale");
    XLAL_CHECK_FAIL(scale, XLAL_EFAILED, "Could not get attribute .scale");
    XLAL_CHECK_FAIL(PyFloat_Check(scale), XLAL_EFAILED, "Attribute .scale is not a float");
    scale_factor = PyFloat_AsDouble(scale);
    if (PyErr_Occurred())
        XLAL_ERROR_FAIL(XLAL_EFAILED, "Could not convert attribute .scale to double");
    Py_CLEAR(scale);

    bases = PyObject_GetAttrString(si, "bases");
    XLAL_CHECK_FAIL(bases, XLAL_EFAILED, "Could not get attribute .bases");

    powers = PyObject_GetAttrString(si, "powers");
    XLAL_CHECK_FAIL(powers, XLAL_EFAILED, "Could not get attribute .powers");

    Py_CLEAR(si);

    if (!PyList_Check(bases) || !PyList_Check(powers))
        XLAL_ERROR_FAIL(XLAL_EFAILED, "Invalid astropy unit data");
    if (PyList_Size(bases) != PyList_Size(powers))
        XLAL_ERROR_FAIL(XLAL_EFAILED, "Internal data inconsistency");

    /* The bases attribute contains a sequence of units this unit is composed of.
     * We loop over each of them and transform them into a LALUnit.
     * The powers item tell us the power of the corresponding base unit (this is a general convention to describe units) and to account for
     * float powers they are decomposed into a fraction of integer numbers (numerator and denominator) which are passed then to the LALUnit struct.
     */
    memset(lu, 0, sizeof(*lu));
    for (Py_ssize_t i = 0; i < PyList_Size(bases); ++i) {
        PyObject *base = PyList_GetItem(bases, i);  // borrowed reference
        PyObject *power = PyList_GetItem(powers, i);  // borrowed reference
        PyObject *intratio = NULL;
        PyObject *basename = NULL;
        int numerator;
        int denominator;

        intratio = PyObject_CallMethod(power, "as_integer_ratio", NULL);
        if (intratio == NULL || !PyTuple_Check(intratio)) {
            Py_XDECREF(intratio);
            XLAL_ERROR_FAIL(XLAL_EFAILED, "Failed method .as_integer_ratio");
        }
        if (!PyArg_ParseTuple(intratio, "ii", &numerator, &denominator)) {
            Py_DECREF(intratio);
            XLAL_ERROR_FAIL(XLAL_EFAILED, "Failed to unpack integer ratio");
        }
        Py_DECREF(intratio);

        basename = PyObject_GetAttrString(base, "name");
        if (basename == NULL || !PyUnicode_Check(basename)) {
            Py_XDECREF(basename);
            XLAL_ERROR_FAIL(XLAL_EFAILED, "Failed get attribute .name");
        }

        if (PyUnicode_CompareWithASCIIString(basename, "m") == 0) {
            lu->unitNumerator[LALUnitIndexMeter] = numerator;
            lu->unitDenominatorMinusOne[LALUnitIndexMeter] = denominator - 1;
        } else if (PyUnicode_CompareWithASCIIString(basename, "kg") == 0) {
            lu->unitNumerator[LALUnitIndexKiloGram] = numerator;
            lu->unitDenominatorMinusOne[LALUnitIndexKiloGram] = denominator - 1;
        } else if (PyUnicode_CompareWithASCIIString(basename, "s") == 0) {
            lu->unitNumerator[LALUnitIndexSecond] = numerator;
            lu->unitDenominatorMinusOne[LALUnitIndexSecond] = denominator - 1;
        } else if (PyUnicode_CompareWithASCIIString(basename, "A") == 0) {
            lu->unitNumerator[LALUnitIndexAmpere] = numerator;
            lu->unitDenominatorMinusOne[LALUnitIndexAmpere] = denominator - 1;
        } else if (PyUnicode_CompareWithASCIIString(basename, "K") == 0) {
            lu->unitNumerator[LALUnitIndexKelvin] = numerator;
            lu->unitDenominatorMinusOne[LALUnitIndexKelvin] = denominator - 1;
        } else if (PyUnicode_CompareWithASCIIString(basename, "strain") == 0) {
            lu->unitNumerator[LALUnitIndexStrain] = numerator;
            lu->unitDenominatorMinusOne[LALUnitIndexStrain] = denominator - 1;
        } else if (PyUnicode_CompareWithASCIIString(basename, "ct") == 0) {
            lu->unitNumerator[LALUnitIndexADCCount] = numerator;
            lu->unitDenominatorMinusOne[LALUnitIndexADCCount] = denominator - 1;
        } else if (PyUnicode_CompareWithASCIIString(basename, "rad") == 0) {
            ; /* this is just dimensionless in LAL */
        } else {
            Py_DECREF(basename);
            XLAL_ERROR_FAIL(XLAL_EFAILED, "Astropy unit has no associated LAL unit");
        }

        Py_DECREF(basename);
    }

    Py_DECREF(bases);
    Py_DECREF(powers);
    return scale_factor;

XLAL_FAIL:
    python_error();
    Py_XDECREF(si);
    Py_XDECREF(bases);
    Py_XDECREF(powers);
    return XLAL_REAL8_FAIL_NAN;
}


/* ROUTINES TO CONVERT GWPY TIME/FREQUENCY SERIES AS LAL SERIES */

#define DEFINE_GWPY_TIMESERIES_TO_LAL_TIMESERIES(TYPE) \
    static TYPE ## TimeSeries *GWpyTimeSeries_As ## TYPE ## TimeSeries(PyObject *gwpyser) \
    { \
        TYPE ## TimeSeries *series; \
        double scale; \
        PyObject *name = NULL; \
        PyObject *bytes = NULL; \
        PyObject *epoch = NULL; \
        PyObject *dt = NULL; \
        PyObject *value = NULL; \
        PyObject *units = NULL; \
        PyObject *si = NULL; \
     \
        series = XLALCalloc(1, sizeof(*series)); \
        XLAL_CHECK_NULL(series, XLAL_ENOMEM); \
     \
        /* extract the metadata */ \
     \
        name = PyObject_GetAttrString(gwpyser, "name"); \
        XLAL_CHECK_FAIL(name, XLAL_EFAILED, "Could not get attribute .name"); \
        XLAL_CHECK_FAIL(PyUnicode_Check(name), XLAL_EFAILED, "Attribute .name is not a string"); \
        bytes = PyUnicode_AsASCIIString(name); \
        XLAL_CHECK_FAIL(PyBytes_Check(bytes), XLAL_EFAILED, "Could not decode attribute .name as ASCII string"); \
        Py_CLEAR(name); \
        XLALStringCopy(series->name, PyBytes_AsString(bytes), sizeof(series->name)); \
        Py_CLEAR(bytes); \
     \
        epoch = PyObject_GetAttrString(gwpyser, "epoch"); \
        XLAL_CHECK_FAIL(epoch, XLAL_EFAILED, "Could not get attribute .epoch"); \
        if (AstropyTime_AsLIGOTimeGPS(&series->epoch, epoch) < 0) \
            XLAL_ERROR_FAIL(XLAL_EFUNC); \
        Py_CLEAR(epoch); \
     \
        dt = PyObject_GetAttrString(gwpyser, "dt"); \
        XLAL_CHECK_FAIL(dt, XLAL_EFAILED, "Could not get attribute .dt"); \
        si = PyObject_GetAttrString(dt, "si"); \
        XLAL_CHECK_FAIL(si, XLAL_EFAILED, "Could not get attribute .si"); \
        Py_CLEAR(dt); \
        value = PyObject_GetAttrString(si, "value"); \
        XLAL_CHECK_FAIL(value, XLAL_EFAILED, "Could not get attribute .value"); \
        Py_CLEAR(si); \
        XLAL_CHECK_FAIL(PyFloat_Check(value), XLAL_EFAILED, "Attribute .value is not a float"); \
        series->deltaT = PyFloat_AsDouble(value); \
        if (PyErr_Occurred()) \
            XLAL_ERROR_FAIL(XLAL_EFAILED, "Could not convert attribute .value to double"); \
        Py_CLEAR(value); \
     \
        /* extract the data */ \
        value = PyObject_GetAttrString(gwpyser, "value"); \
        XLAL_CHECK_FAIL(value, XLAL_EFAILED, "Could not get attribute .value"); \
        units = PyObject_GetAttrString(gwpyser, "unit"); \
        XLAL_CHECK_FAIL(units, XLAL_EFAILED, "Could not get attribute .unit"); \
        scale = AstropyUnit_AsLALUnit(&series->sampleUnits, units); \
        if (XLAL_IS_REAL8_FAIL_NAN(scale)) \
            XLAL_ERROR_FAIL(XLAL_EFUNC); \
        Py_CLEAR(units); \
        series->data = PyArray_As ## TYPE ## Sequence(value); \
        XLAL_CHECK_FAIL(series->data, XLAL_EFUNC); \
        Py_CLEAR(value); \
        if (scale != 1.0) \
            for (ssize_t i = 0; i < series->data->length; ++i) \
                series->data->data[i] *= scale; \
     \
        return series; \
     \
    XLAL_FAIL: \
        python_error(); \
        Py_XDECREF(name); \
        Py_XDECREF(bytes); \
        Py_XDECREF(epoch); \
        Py_XDECREF(dt); \
        Py_XDECREF(value); \
        Py_XDECREF(units); \
        Py_XDECREF(si); \
        XLALDestroy ## TYPE ## TimeSeries(series); \
        return NULL; \
    }

#define DEFINE_GWPY_FREQUENCYSERIES_TO_LAL_FREQUENCYSERIES(TYPE) \
    static TYPE ## FrequencySeries *GWpyFrequencySeries_As ## TYPE ## FrequencySeries(PyObject *gwpyser) \
    { \
        TYPE ## FrequencySeries *series; \
        double scale; \
        PyObject *name = NULL; \
        PyObject *bytes = NULL; \
        PyObject *epoch = NULL; \
        PyObject *f0 = NULL; \
        PyObject *df = NULL; \
        PyObject *value = NULL; \
        PyObject *units = NULL; \
        PyObject *si = NULL; \
     \
        series = XLALCalloc(1, sizeof(*series)); \
        XLAL_CHECK_NULL(series, XLAL_ENOMEM); \
     \
        /* extract the metadata */ \
     \
        name = PyObject_GetAttrString(gwpyser, "name"); \
        XLAL_CHECK_FAIL(name, XLAL_EFAILED, "Could not get attribute .name"); \
        XLAL_CHECK_FAIL(PyUnicode_Check(name), XLAL_EFAILED, "Attribute .name is not a string"); \
        bytes = PyUnicode_AsASCIIString(name); \
        XLAL_CHECK_FAIL(PyBytes_Check(bytes), XLAL_EFAILED, "Could not decode attribute .name as ASCII string"); \
        Py_CLEAR(name); \
        XLALStringCopy(series->name, PyBytes_AsString(bytes), sizeof(series->name)); \
        Py_CLEAR(bytes); \
     \
        epoch = PyObject_GetAttrString(gwpyser, "epoch"); \
        XLAL_CHECK_FAIL(epoch, XLAL_EFAILED, "Could not get attribute .epoch"); \
        if (AstropyTime_AsLIGOTimeGPS(&series->epoch, epoch) < 0) \
            XLAL_ERROR_FAIL(XLAL_EFUNC); \
        Py_CLEAR(epoch); \
     \
        f0 = PyObject_GetAttrString(gwpyser, "f0"); \
        XLAL_CHECK_FAIL(f0, XLAL_EFAILED, "Could not get attribute .f0"); \
        si = PyObject_GetAttrString(f0, "si"); \
        XLAL_CHECK_FAIL(si, XLAL_EFAILED, "Could not get attribute .si"); \
        Py_CLEAR(f0); \
        value = PyObject_GetAttrString(si, "value"); \
        XLAL_CHECK_FAIL(value, XLAL_EFAILED, "Could not get attribute .value"); \
        Py_CLEAR(si); \
        XLAL_CHECK_FAIL(PyFloat_Check(value), XLAL_EFAILED, "Attribute .value is not a float"); \
        series->f0 = PyFloat_AsDouble(value); \
        if (PyErr_Occurred()) \
            XLAL_ERROR_FAIL(XLAL_EFAILED, "Could not convert attribute .value to double"); \
        Py_CLEAR(value); \
     \
        df = PyObject_GetAttrString(gwpyser, "df"); \
        XLAL_CHECK_FAIL(df, XLAL_EFAILED, "Could not get attribute .df"); \
        si = PyObject_GetAttrString(df, "si"); \
        XLAL_CHECK_FAIL(si, XLAL_EFAILED, "Could not get attribute .si"); \
        Py_CLEAR(df); \
        value = PyObject_GetAttrString(si, "value"); \
        XLAL_CHECK_FAIL(value, XLAL_EFAILED, "Could not get attribute .value"); \
        Py_CLEAR(si); \
        XLAL_CHECK_FAIL(PyFloat_Check(value), XLAL_EFAILED, "Attribute .value is not a float"); \
        series->deltaF = PyFloat_AsDouble(value); \
        if (PyErr_Occurred()) \
            XLAL_ERROR_FAIL(XLAL_EFAILED, "Could not convert attribute .value to double"); \
        Py_CLEAR(value); \
     \
        /* extract the data */ \
        value = PyObject_GetAttrString(gwpyser, "value"); \
        XLAL_CHECK_FAIL(value, XLAL_EFAILED, "Could not get attribute .value"); \
        units = PyObject_GetAttrString(gwpyser, "unit"); \
        XLAL_CHECK_FAIL(units, XLAL_EFAILED, "Could not get attribute .unit"); \
        scale = AstropyUnit_AsLALUnit(&series->sampleUnits, units); \
        if (XLAL_IS_REAL8_FAIL_NAN(scale)) \
            XLAL_ERROR_FAIL(XLAL_EFUNC); \
        Py_CLEAR(units); \
        series->data = PyArray_As ## TYPE ## Sequence(value); \
        XLAL_CHECK_FAIL(series->data, XLAL_EFUNC); \
        Py_CLEAR(value); \
        if (scale != 1.0) \
            for (ssize_t i = 0; i < series->data->length; ++i) \
                series->data->data[i] *= scale; \
     \
        return series; \
     \
    XLAL_FAIL: \
        python_error(); \
        Py_XDECREF(name); \
        Py_XDECREF(bytes); \
        Py_XDECREF(epoch); \
        Py_XDECREF(f0); \
        Py_XDECREF(df); \
        Py_XDECREF(value); \
        Py_XDECREF(units); \
        Py_XDECREF(si); \
        XLALDestroy ## TYPE ## FrequencySeries(series); \
        return NULL; \
    }

DEFINE_GWPY_TIMESERIES_TO_LAL_TIMESERIES(REAL8)
DEFINE_GWPY_TIMESERIES_TO_LAL_TIMESERIES(COMPLEX16)
DEFINE_GWPY_FREQUENCYSERIES_TO_LAL_FREQUENCYSERIES(COMPLEX16)


/* STANDARD GENERATOR METHODS */

static int generate_td_modes(
    SphHarmTimeSeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *myself
);

static int generate_td_waveform(
    REAL8TimeSeries **hplus,
    REAL8TimeSeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *myself
);

static int generate_fd_modes(
    SphHarmFrequencySeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *myself
);

static int generate_fd_waveform(
    COMPLEX16FrequencySeries **hplus,
    COMPLEX16FrequencySeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *myself
);

/* Helper struct to store the python quantity and objects used for waveform generation */
struct internal_data {
    PyObject *module;
    PyObject *object;
    PyObject *instance;
    PyObject *generate_td_modes;
    PyObject *generate_fd_modes;
    PyObject *generate_td_waveform;
    PyObject *generate_fd_waveform;
    PyObject *quantity;
    LALDict *waveform_parameter_units;
};

/* Free internal_data struct */
static void free_internal_data(struct internal_data *internal_data)
{
    if (internal_data) {
        Py_XDECREF(internal_data->module);
        Py_XDECREF(internal_data->object);
        Py_XDECREF(internal_data->instance);
        Py_XDECREF(internal_data->generate_td_modes);
        Py_XDECREF(internal_data->generate_fd_modes);
        Py_XDECREF(internal_data->generate_td_waveform);
        Py_XDECREF(internal_data->generate_fd_waveform);
        Py_XDECREF(internal_data->quantity);
        XLALDestroyDict(internal_data->waveform_parameter_units);
        XLALFree(internal_data);
    }
}

/* Initialize python generator.
 * The path to the python code that actually generates the model is stored in internal_data.module and .object
 */
static int initialize(LALSimInspiralGenerator *myself, LALDict *params)
{
    struct internal_data *internal_data = NULL;
    const char *module_name;
    const char *object_name;
    LALDict *genparams = NULL;
    PyObject *args = NULL;
    PyObject *kwargs = NULL;
    PyObject *fromlist = NULL;
    PyObject *result = NULL;

    XLAL_CHECK(myself, XLAL_EFAULT);
    XLAL_CHECK(params, XLAL_EINVAL);

    module_name = XLALDictLookupStringValue(params, "module");
    object_name = XLALDictLookupStringValue(params, "object");
    XLAL_CHECK(module_name, XLAL_EINVAL, "No parameter \"module\" found");
    XLAL_CHECK(object_name, XLAL_EINVAL, "No parameter \"object\" found");

    Py_Initialize();
    if (_import_array() < 0) {
        PyErr_Print();
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        XLAL_ERROR(XLAL_EFAILED, "Failed to import numpy.core.multiarray");
    }

    internal_data = XLALCalloc(1, sizeof(struct internal_data));
    XLAL_CHECK(internal_data, XLAL_ENOMEM);

    internal_data->waveform_parameter_units = create_waveform_parameter_units_dict();
    XLAL_CHECK_FAIL(internal_data->waveform_parameter_units, XLAL_EFUNC);

    /* from astropy.units import Quantity */
    fromlist = Py_BuildValue("(s)", "Quantity");
    XLAL_CHECK_FAIL(fromlist, XLAL_EFAILED, "Could not build fromlist");
    result = PyImport_ImportModuleEx("astropy.units", NULL, NULL, fromlist);
    XLAL_CHECK_FAIL(result, XLAL_EFAILED, "Could not load module astropy.units");
    Py_CLEAR(fromlist);

    internal_data->quantity = PyObject_GetAttrString(result, "Quantity");
    XLAL_CHECK_FAIL(internal_data->quantity, XLAL_EFAILED, "Could not load Quantity from module astropy.units");
    Py_CLEAR(result);
    
    /* import specified module */
    internal_data->module = PyImport_ImportModule(module_name);
    XLAL_CHECK_FAIL(internal_data->module, XLAL_EFAILED, "Could not find module %s", module_name);

    /* load specified object from module */
    internal_data->object = PyObject_GetAttrString(internal_data->module, object_name);
    XLAL_CHECK_FAIL(internal_data->object, XLAL_EFAILED, "Could not find object %s in module %s", object_name, module_name);

    /* parameters to pass off to the python generator's .__init__ */
    genparams = XLALDictDuplicate(params);
    XLALDictRemove(genparams, "module");
    XLALDictRemove(genparams, "object");
    kwargs = PyDict_FromLALDict(genparams, NULL, NULL);
    XLAL_CHECK_FAIL(kwargs, XLAL_EFUNC);
    XLALDestroyDict(genparams);
    genparams = NULL;
    args = PyTuple_New(0); /* empty args */
    XLAL_CHECK_FAIL(args, XLAL_EFAILED, "Could not create tuple");

    /* create an instance of the specified object from module */
    internal_data->instance = PyObject_Call(internal_data->object, args, kwargs);
    XLAL_CHECK_FAIL(internal_data->instance, XLAL_EFAILED, "Could not create instance of object %s in module %s", object_name, module_name);
    Py_CLEAR(args);
    Py_CLEAR(kwargs);
    
    /* figure out what methods the object supports */

    /* generate_td_waveform */
    if (PyObject_HasAttrString(internal_data->object, "generate_td_waveform")) {
        internal_data->generate_td_waveform = PyObject_GetAttrString(internal_data->object, "generate_td_waveform");
        XLAL_CHECK_FAIL(internal_data->generate_td_waveform, XLAL_EFAILED, "Failed to get method generate_td_waveform");
    } else
        internal_data->generate_td_waveform = NULL;

    /* generate_fd_waveform */
    if (PyObject_HasAttrString(internal_data->object, "generate_fd_waveform")) {
        internal_data->generate_fd_waveform = PyObject_GetAttrString(internal_data->object, "generate_fd_waveform");
        XLAL_CHECK_FAIL(internal_data->generate_fd_waveform, XLAL_EFAILED, "Failed to get method generate_fd_waveform");
    } else
        internal_data->generate_fd_waveform = NULL;

    /* generate_td_modes */
    if (PyObject_HasAttrString(internal_data->object, "generate_td_modes")) {
        internal_data->generate_td_modes = PyObject_GetAttrString(internal_data->object, "generate_td_modes");
        XLAL_CHECK_FAIL(internal_data->generate_td_modes, XLAL_EFAILED, "Failed to get method generate_td_modes");
    } else
        internal_data->generate_td_modes = NULL;

    /* generate_fd_modes */
    if (PyObject_HasAttrString(internal_data->object, "generate_fd_modes")) {
        internal_data->generate_fd_modes = PyObject_GetAttrString(internal_data->object, "generate_fd_modes");
        XLAL_CHECK_FAIL(internal_data->generate_fd_modes, XLAL_EFAILED, "Failed to get method generate_fd_modes");
    } else
        internal_data->generate_fd_modes = NULL;

    /* make sure that at least one method is defined */
    XLAL_CHECK_FAIL(internal_data->generate_td_waveform || internal_data->generate_fd_waveform || internal_data->generate_td_modes || internal_data->generate_fd_modes, XLAL_EFAILED, "Object %s in module %s must provide at least one following methods: generate_td_waveform, generate_fd_waveform, generate_td_modes, generate_fd_modes", object_name, module_name);

    /* attach those methods that are found to myself */
    if (internal_data->generate_td_waveform)
        myself->generate_td_waveform = generate_td_waveform;
    if (internal_data->generate_fd_waveform)
        myself->generate_fd_waveform = generate_fd_waveform;
    if (internal_data->generate_td_modes)
        myself->generate_td_modes = generate_td_modes;
    if (internal_data->generate_fd_modes)
        myself->generate_fd_modes = generate_fd_modes;

    myself->internal_data = internal_data;
    return 0;

XLAL_FAIL:
    python_error();
    free_internal_data(internal_data);
    XLALDestroyDict(genparams);
    Py_XDECREF(args);
    Py_XDECREF(kwargs);
    Py_XDECREF(result);
    Py_XDECREF(fromlist);
    return (int)XLAL_FAILURE;
}

/* Free memory */
static int finalize(LALSimInspiralGenerator *myself)
{
    XLAL_CHECK(myself, XLAL_EFAULT); 
    XLAL_CHECK(myself->internal_data, XLAL_EINVAL); 
    free_internal_data(myself->internal_data);
    myself->generate_td_waveform = NULL;
    myself->generate_fd_waveform = NULL;
    myself->generate_td_modes = NULL;
    myself->generate_fd_modes = NULL;
    myself->internal_data = NULL;
    /* don't Py_Finalize() since other things might be using the interpreter */
    // Py_Finalize();
    return 0;
}

/* Generate time-domain polarizations for python generator */
static int generate_td_waveform(
    REAL8TimeSeries **hplus,
    REAL8TimeSeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *myself
)
{
    struct internal_data *internal_data;
    PyObject *method;
    PyObject *args = NULL;
    PyObject *kwargs = NULL;
    PyObject *result = NULL;
    PyObject *hp; /* borrowed reference */
    PyObject *hc; /* borrowed reference */

    /* sanity check arguments */
    XLAL_CHECK(hplus, XLAL_EFAULT);
    XLAL_CHECK(hcross, XLAL_EFAULT);
    XLAL_CHECK(params, XLAL_EFAULT);
    XLAL_CHECK(myself, XLAL_EFAULT);
    XLAL_CHECK(*hplus == NULL, XLAL_EINVAL, "Argument hplus must be a pointer to NULL");
    XLAL_CHECK(*hcross == NULL, XLAL_EINVAL, "Argument hcross must be a pointer to NULL");

    internal_data = myself->internal_data;
    XLAL_CHECK(internal_data, XLAL_EINVAL);

    method = internal_data->generate_td_waveform;
    XLAL_CHECK(method, XLAL_EFAILED, "Object does not provide method %s", __func__);

    /* call method */
    args = PyTuple_Pack(1, internal_data->instance); /* self */
    XLAL_CHECK_FAIL(args, XLAL_EFAILED, "Failed to create tuple");
    kwargs = PyDict_FromLALDict(params, internal_data->quantity, internal_data->waveform_parameter_units);
    XLAL_CHECK_FAIL(kwargs, XLAL_EFUNC);
    result = PyObject_Call(method, args, kwargs);
    XLAL_CHECK_FAIL(result, XLAL_EFUNC, "Failed to call method %s", __func__);
    Py_CLEAR(kwargs);
    Py_CLEAR(args);

    /* unpack result */
    /* should be a tuple of size 2 */
    XLAL_CHECK_FAIL(PyTuple_Check(result), XLAL_EFAILED, "Invalid value returned by method %s", __func__);
    XLAL_CHECK_FAIL(PyTuple_Size(result) == 2, XLAL_EFAILED, "Invalid value returned by method %s", __func__);
    hp = PyTuple_GetItem(result, 0);
    hc = PyTuple_GetItem(result, 1);
    XLAL_CHECK_FAIL(hp && hc, XLAL_EFAILED, "Invalid value returned by method %s", __func__);

    /* convert results */
    *hplus = GWpyTimeSeries_AsREAL8TimeSeries(hp);
    XLAL_CHECK_FAIL(*hplus, XLAL_EFUNC);
    *hcross = GWpyTimeSeries_AsREAL8TimeSeries(hc);
    XLAL_CHECK_FAIL(*hcross, XLAL_EFUNC);

    Py_CLEAR(result);
    Py_CLEAR(hp);
    Py_CLEAR(hc);
    return 0;

XLAL_FAIL:
    python_error();
    XLALDestroyREAL8TimeSeries(*hcross);
    XLALDestroyREAL8TimeSeries(*hplus);
    Py_XDECREF(args);
    Py_XDECREF(kwargs);
    Py_XDECREF(result);
    return (int)XLAL_FAILURE;
}

/* Generate time-domain modes for python generator */
static int generate_td_modes(
    SphHarmTimeSeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *myself
)
{
    struct internal_data *internal_data;
    PyObject *method;
    PyObject *args = NULL;
    PyObject *kwargs = NULL;
    PyObject *result = NULL;
    PyObject *key; /* borrowed reference */
    PyObject *val; /* borrowed reference */
    Py_ssize_t pos = 0;

    /* sanity check arguments */
    XLAL_CHECK(hlm, XLAL_EFAULT);
    XLAL_CHECK(params, XLAL_EFAULT);
    XLAL_CHECK(myself, XLAL_EFAULT);
    XLAL_CHECK(*hlm == NULL, XLAL_EINVAL, "Argument hplus must be a pointer to NULL");

    internal_data = myself->internal_data;
    XLAL_CHECK(internal_data, XLAL_EINVAL);

    method = internal_data->generate_td_modes;
    XLAL_CHECK(method, XLAL_EFAILED, "Object does not provide method %s", __func__);

    /* call method */
    args = PyTuple_Pack(1, internal_data->instance); /* self */
    XLAL_CHECK_FAIL(args, XLAL_EFAILED, "Failed to create tuple");
    kwargs = PyDict_FromLALDict(params, internal_data->quantity, internal_data->waveform_parameter_units);
    XLAL_CHECK_FAIL(kwargs, XLAL_EFUNC);
    result = PyObject_Call(method, args, kwargs);
    XLAL_CHECK_FAIL(result, XLAL_EFUNC, "Failed to call method %s", __func__);
    Py_CLEAR(kwargs);
    Py_CLEAR(args);

    /* unpack result */
    /* make sure it is a dict */
    XLAL_CHECK_FAIL(PyDict_Check(result), XLAL_EFAILED, "Invalid value returned by method %s", __func__);

    /* loop over items in dict */
    while (PyDict_Next(result, &pos, &key, &val)) {
        SphHarmTimeSeries *this;
        PyObject *l; /* borrowed reference */
        PyObject *m; /* borrowed reference */

        /* prepend this mode */
        this = XLALCalloc(1, sizeof(*this));
        XLAL_CHECK_FAIL(this, XLAL_ENOMEM);
        this->next = *hlm;
        *hlm = this;

        /* key should be a tuple of size 2 */
        XLAL_CHECK_FAIL(PyTuple_Check(key) && PyTuple_Size(key) == 2, XLAL_EFAILED, "Invalid value returned by method %s", __func__);
        l = PyTuple_GetItem(key, 0);
        m = PyTuple_GetItem(key, 1);
        XLAL_CHECK_FAIL(l && m && PyLong_Check(l) && PyLong_Check(m), XLAL_EFAILED, "Invalid value returned by method %s", __func__);

        /* convert results */
        this->l = PyLong_AsUnsignedLong(l);
        XLAL_CHECK_FAIL(!PyErr_Occurred(), XLAL_EFAILED, "Invalid value returned by method %s", __func__);
        this->m = PyLong_AsLong(m);
        XLAL_CHECK_FAIL(!PyErr_Occurred(), XLAL_EFAILED, "Invalid value returned by method %s", __func__);
        this->mode = GWpyTimeSeries_AsCOMPLEX16TimeSeries(val);
        XLAL_CHECK_FAIL(this->mode, XLAL_EFUNC);
    }

    Py_CLEAR(result);
    return 0;

XLAL_FAIL:
    python_error();
    XLALDestroySphHarmTimeSeries(*hlm);
    *hlm = NULL;
    Py_XDECREF(args);
    Py_XDECREF(kwargs);
    Py_XDECREF(result);
    return (int)XLAL_FAILURE;
}

/* Generate Fourier-domain polarizations for python generator */
static int generate_fd_waveform(
    COMPLEX16FrequencySeries **hplus,
    COMPLEX16FrequencySeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *myself
)
{
    struct internal_data *internal_data;
    PyObject *method;
    PyObject *args = NULL;
    PyObject *kwargs = NULL;
    PyObject *result = NULL;
    PyObject *hp; /* borrowed reference */
    PyObject *hc; /* borrowed reference */

    /* sanity check arguments */
    XLAL_CHECK(hplus, XLAL_EFAULT);
    XLAL_CHECK(hcross, XLAL_EFAULT);
    XLAL_CHECK(params, XLAL_EFAULT);
    XLAL_CHECK(myself, XLAL_EFAULT);
    XLAL_CHECK(*hplus == NULL, XLAL_EINVAL, "Argument hplus must be a pointer to NULL");
    XLAL_CHECK(*hcross == NULL, XLAL_EINVAL, "Argument hcross must be a pointer to NULL");

    internal_data = myself->internal_data;
    XLAL_CHECK(internal_data, XLAL_EINVAL);

    method = internal_data->generate_fd_waveform;
    XLAL_CHECK(method, XLAL_EFAILED, "Object does not provide method %s", __func__);

    /* call method */
    args = PyTuple_Pack(1, internal_data->instance); /* self */
    XLAL_CHECK_FAIL(args, XLAL_EFAILED, "Failed to create tuple");
    kwargs = PyDict_FromLALDict(params, internal_data->quantity, internal_data->waveform_parameter_units);
    XLAL_CHECK_FAIL(kwargs, XLAL_EFUNC);
    result = PyObject_Call(method, args, kwargs);
    XLAL_CHECK_FAIL(result, XLAL_EFUNC, "Failed to call method %s", __func__);
    Py_CLEAR(kwargs);
    Py_CLEAR(args);

    /* unpack result */
    /* should be a tuple of size 2 */
    XLAL_CHECK_FAIL(PyTuple_Check(result), XLAL_EFAILED, "Invalid value returned by method %s", __func__);
    XLAL_CHECK_FAIL(PyTuple_Size(result) == 2, XLAL_EFAILED, "Invalid value returned by method %s", __func__);
    hp = PyTuple_GetItem(result, 0);
    hc = PyTuple_GetItem(result, 1);
    XLAL_CHECK_FAIL(hp && hc, XLAL_EFAILED, "Invalid value returned by method %s", __func__);

    /* convert results */
    *hplus = GWpyFrequencySeries_AsCOMPLEX16FrequencySeries(hp);
    XLAL_CHECK_FAIL(*hplus, XLAL_EFUNC);
    *hcross = GWpyFrequencySeries_AsCOMPLEX16FrequencySeries(hc);
    XLAL_CHECK_FAIL(*hcross, XLAL_EFUNC);


    Py_CLEAR(result);
    Py_CLEAR(hp);
    Py_CLEAR(hc);
    return 0;

XLAL_FAIL:
    python_error();
    XLALDestroyCOMPLEX16FrequencySeries(*hcross);
    XLALDestroyCOMPLEX16FrequencySeries(*hplus);
    Py_XDECREF(args);
    Py_XDECREF(kwargs);
    Py_XDECREF(result);
    return (int)XLAL_FAILURE;
}

/* Generate Fourier-domain modes for python generator */
static int generate_fd_modes(
    SphHarmFrequencySeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *myself
)
{
    struct internal_data *internal_data;
    PyObject *method;
    PyObject *args = NULL;
    PyObject *kwargs = NULL;
    PyObject *result = NULL;
    PyObject *key; /* borrowed reference */
    PyObject *val; /* borrowed reference */
    Py_ssize_t pos = 0;

    /* sanity check arguments */
    XLAL_CHECK(hlm, XLAL_EFAULT);
    XLAL_CHECK(params, XLAL_EFAULT);
    XLAL_CHECK(myself, XLAL_EFAULT);
    XLAL_CHECK(*hlm == NULL, XLAL_EINVAL, "Argument hplus must be a pointer to NULL");

    internal_data = myself->internal_data;
    XLAL_CHECK(internal_data, XLAL_EINVAL);

    method = internal_data->generate_fd_modes;
    XLAL_CHECK(method, XLAL_EFAILED, "Object does not provide method %s", __func__);

    /* call method */
    args = PyTuple_Pack(1, internal_data->instance); /* self */
    XLAL_CHECK_FAIL(args, XLAL_EFAILED, "Failed to create tuple");
    kwargs = PyDict_FromLALDict(params, internal_data->quantity, internal_data->waveform_parameter_units);
    XLAL_CHECK_FAIL(kwargs, XLAL_EFUNC);
    result = PyObject_Call(method, args, kwargs);
    XLAL_CHECK_FAIL(result, XLAL_EFUNC, "Failed to call method %s", __func__);
    Py_CLEAR(kwargs);
    Py_CLEAR(args);

    /* unpack result */
    /* make sure it is a dict */
    XLAL_CHECK_FAIL(PyDict_Check(result), XLAL_EFAILED, "Invalid value returned by method %s", __func__);

    /* loop over items in dict */
    while (PyDict_Next(result, &pos, &key, &val)) {
        SphHarmFrequencySeries *this;
        PyObject *l; /* borrowed reference */
        PyObject *m; /* borrowed reference */

        /* prepend this mode */
        this = XLALCalloc(1, sizeof(*this));
        XLAL_CHECK_FAIL(this, XLAL_ENOMEM);
        this->next = *hlm;
        *hlm = this;

        /* key should be a tuple of size 2 */
        XLAL_CHECK_FAIL(PyTuple_Check(key) && PyTuple_Size(key) == 2, XLAL_EFAILED, "Invalid value returned by method %s", __func__);
        l = PyTuple_GetItem(key, 0);
        m = PyTuple_GetItem(key, 1);
        XLAL_CHECK_FAIL(l && m && PyLong_Check(l) && PyLong_Check(m), XLAL_EFAILED, "Invalid value returned by method %s", __func__);

        /* convert results */
        this->l = PyLong_AsUnsignedLong(l);
        XLAL_CHECK_FAIL(!PyErr_Occurred(), XLAL_EFAILED, "Invalid value returned by method %s", __func__);
        this->m = PyLong_AsLong(m);
        XLAL_CHECK_FAIL(!PyErr_Occurred(), XLAL_EFAILED, "Invalid value returned by method %s", __func__);
        this->mode = GWpyFrequencySeries_AsCOMPLEX16FrequencySeries(val);
        XLAL_CHECK_FAIL(this->mode, XLAL_EFUNC);
    }

    Py_CLEAR(result);
    return 0;

XLAL_FAIL:
    python_error();
    XLALDestroySphHarmFrequencySeries(*hlm);
    *hlm = NULL;
    Py_XDECREF(args);
    Py_XDECREF(kwargs);
    Py_XDECREF(result);
    return (int)XLAL_FAILURE;
}

/* Define methods supported by the python generator.
 * For any external python model, the approximant name must be ExternalPython.
 * The path to the python code that actually generates the model is stored in internal_data.module and .object
 */
const LALSimInspiralGenerator lalPythonGeneratorTemplate = {
    .name = "ExternalPython",
    .initialize = initialize,
    .finalize = finalize,
    .generate_td_waveform = generate_td_waveform,
    .generate_fd_waveform = generate_fd_waveform,
    .generate_td_modes = generate_td_modes,
    .generate_fd_modes = generate_fd_modes,
    .internal_data = NULL
};

#endif /* HAVE_PYTHON_EMBED */
