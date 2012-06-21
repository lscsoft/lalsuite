/*
 *
 * Copyright (C) 2008-9 Antony Searle
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#ifndef SKYMAP_H
#define SKYMAP_H

/* Special cases for Python bindings of functions defined in this header.
 * FIXME: add Octave section, or modify function declarations to work more
 * smoothly with stock swiglal typemaps. */
#if defined(SWIG) && defined(SWIGPYTHON)
// special typemaps for XLALSkymapPlanConstruct
%typemap(in, noblock=1) (int n, int *detectors) (int i, PyObject *seq) {
    seq = PySequence_Fast($input, "expected detectors to be a sequence");
    if (!seq) SWIG_fail;

    $1 = PySequence_Fast_GET_SIZE(seq);
    $2 = ($2_type) calloc($1, sizeof($*2_type));
    if (!$2)
    {
        Py_DECREF(seq);
        SWIG_exception(SWIG_MemoryError,
            "failed to allocate temporary memory");
    }

    for (i = 0; i < $1; i ++)
    {
        long longValue = PyInt_AsLong(PySequence_Fast_GET_ITEM(seq, i));
        if (longValue == -1 && PyErr_Occurred())
        {
            Py_DECREF(seq);
            SWIG_fail;
        }
        $2[i] = longValue;
    }
    Py_DECREF(seq);
}
%typemap(freearg, noblock=1) (int n, int *detectors) {
    free($2);
}

// special typemaps for XLALSkymapKernelConstruct
// and XLALSkymapUncertainKernelConstruct
%typemap(in, noblock=1) double *wSw (int i, int n, PyObject *seq) {
    seq = PySequence_Fast($input, "expected direction to be a sequence");
    if (!seq) SWIG_fail;

    n = PySequence_Fast_GET_SIZE(seq);

    $1 = ($1_type) calloc(n, sizeof($*1_type));
    if (!$1)
    {
        Py_DECREF(seq);
        SWIG_exception(SWIG_MemoryError,
            "failed to allocate temporary memory");
    }

    for (i = 0; i < n; i ++)
    {
        double doubleValue = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
        if (doubleValue == -1 && PyErr_Occurred())
        {
            Py_DECREF(seq);
            SWIG_fail;
        }
        $1[i] = doubleValue;
    }

    Py_DECREF(seq);
}
%typemap(freearg, noblock=1) double *wSw {
    free($1);
}
%apply double *wSw { double *error };

// special typemaps for XLALSkymapApply
%typemap(in, noblock=1) double **xSw (PyObject *seq, int n, int i) {
    seq = PySequence_Fast($input, "expected xSw to be a sequence");
    if (!seq) SWIG_fail;

    n = PySequence_Fast_GET_SIZE(seq);

    $1 = ($1_type) calloc(n, sizeof($*1_type));
    if (!$1)
    {
        Py_DECREF(seq);
        SWIG_exception(SWIG_MemoryError,
            "failed to allocate temporary memory");
    }

    for (i = 0; i < n; i ++)
    {
        PyObject *arr;
        int len;
        void *data;

        arr = PyArray_CheckFromAny(PySequence_Fast_GET_ITEM(seq, i),
            PyArray_DescrFromType(NPY_DOUBLE), 1, 1,
            NPY_IN_ARRAY | NPY_NOTSWAPPED, NULL);
        if (!arr)
        {
            Py_DECREF(seq);
            SWIG_fail;
        }

        len = PyArray_SIZE(arr);
        data = malloc(len * sizeof(double));
        if (!data)
        {
            Py_DECREF(arr);
            Py_DECREF(seq);
            SWIG_exception(SWIG_MemoryError,
                "failed to allocate temporary memory");
        }
        memcpy(data, PyArray_DATA(arr), len * sizeof(double));
        Py_DECREF(arr);
        $1[i] = ($*1_type) data;
    }

    Py_DECREF(seq);
}
%typemap(argout, noblock=1) double **xSw {
}
%typemap(freearg, noblock=1) double **xSw (int i) {
    for (i = 0; i < n$argnum; i ++)
        free($1[i]);
    free($1);
}
%typemap(in, noblock=1, numinputs=0) double *logPosterior (double ret) {
    $1 = &ret;
}
%typemap(argout, noblock=1) double *logPosterior {
    $result = PyFloat_FromDouble(*$1);
}
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <lal/LALDetectors.h>

// Manipulate log-represented values without overflowing

double XLALSkymapLogSumExp(double a, double b);
double XLALSkymapLogTotalExp(double* begin, double* end);

// Lightweight coordinate transformations

void XLALSkymapCartesianFromSpherical(double a[3], double b[2]);
void XLALSkymapSphericalFromCartesian(double a[2], double b[3]);

// Cubic interpolation

double XLALSkymapInterpolate(double t, double* x);

// Largest network of interest (needed to allocate storage)

#define XLALSKYMAP_N 5

// Use detector names from lal/packages/tools/include/LALDetectors.h :
//	LAL_TAMA_300_DETECTOR	=	0
//	LAL_VIRGO_DETECTOR      =	1
//	LAL_GEO_600_DETECTOR	=	2
//	LAL_LHO_2K_DETECTOR     =	3
//	LAL_LHO_4K_DETECTOR     =	4
//	LAL_LLO_4K_DETECTOR     =	5

// Struct to store basic properties of the analysis: the sample rate, and
// the detectors involved (described by their LAL ID numbers)

typedef struct tagXLALSkymapPlanType
{
    int sampleFrequency;
    int n;
    LALDetector site[XLALSKYMAP_N];
} XLALSkymapPlanType;

void XLALSkymapPlanConstruct(
    int sampleFrequency,
    int n,
    int* detectors,
    XLALSkymapPlanType* plan
    );

// Struct to store reuseable pre-computed quantities for a specific
// direction, set of detectors, and sample rate

typedef struct tagXLALSkymapDirectionPropertiesType
{
    double f[XLALSKYMAP_N][2];
    double delay[XLALSKYMAP_N];
} XLALSkymapDirectionPropertiesType;

void XLALSkymapDirectionPropertiesConstruct(
    XLALSkymapPlanType* plan,
    double direction[2],
    XLALSkymapDirectionPropertiesType* properties
    );

// Struct to store reuseable pre-computed kernel for a specific direction,
// power spectra, and sample rate

typedef struct tagXLALSkymapKernelType
{
    double k[XLALSKYMAP_N][XLALSKYMAP_N];
    double logNormalization;
} XLALSkymapKernelType;

void XLALSkymapKernelConstruct(
    XLALSkymapPlanType* plan,
    XLALSkymapDirectionPropertiesType* properties,
    double* wSw,
    XLALSkymapKernelType* kernel
    );

void XLALSkymapUncertainKernelConstruct(
    XLALSkymapPlanType* plan,
    XLALSkymapDirectionPropertiesType* properties,
    double* wSw,
    double* error,
    XLALSkymapKernelType* kernel
    );

// Compute the Bayesian marginalization integral for the specified system

void XLALSkymapApply(
    XLALSkymapPlanType* plan,
    XLALSkymapDirectionPropertiesType* properties,
    XLALSkymapKernelType* kernel,
    double** xSw,
    double tau,
    double* logPosterior
    );

#ifdef __cplusplus
}
#endif

#endif // SKYMAP_H

