/*
 * Copyright (C) 2013  Leo Singer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */

#define NPY_NO_DEPRECATED_API 7

#include <Python.h>
#include <numpy/arrayobject.h>
#include <chealpix.h>
#include <gsl/gsl_errno.h>
#include <lal/bayestar_sky_map.h>


/**
 * Premalloced objects.
 */


typedef struct {
    PyObject_HEAD
    void *data;
} premalloced_object;


static void premalloced_dealloc(premalloced_object *self)
{
    free(self->data);
    self->ob_type->tp_free((PyObject *) self);
}


static PyTypeObject premalloced_type = {
    PyObject_HEAD_INIT(NULL)
    .tp_name = "premalloced",
    .tp_basicsize = sizeof(premalloced_object),
    .tp_dealloc = (destructor)premalloced_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = "Pre-malloc'd memory"
};


static PyObject *premalloced_new(void *data)
{
    premalloced_object *obj = PyObject_New(premalloced_object, &premalloced_type);
    if (obj)
        obj->data = data;
    else
        free(data);
    return (PyObject *) obj;
}


static void
my_gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
    PyObject *exception_type;
    switch (gsl_errno)
    {
        case GSL_EINVAL:
            exception_type = PyExc_ValueError;
            break;
        case GSL_ENOMEM:
            exception_type = PyExc_MemoryError;
            break;
        default:
            exception_type = PyExc_ArithmeticError;
            break;
    }
    PyErr_Format(exception_type, "%s:%d: %s\n", file, line, reason);
}


static PyObject *sky_map_tdoa(PyObject *module, PyObject *args, PyObject *kwargs)
{
    long i;
    Py_ssize_t n;
    long nside = -1;
    long npix;
    long nifos = 0;
    double gmst;
    PyObject *toas_obj, *w_toas_obj, *locations_obj;

    PyArrayObject *toas_npy = NULL, *w_toas_npy = NULL, **locations_npy = NULL;

    double *toas;
    double *w_toas;
    const double **locations = NULL;

    npy_intp dims[1];
    PyArrayObject  *out = NULL, *ret = NULL;
    PyObject *premalloced = NULL;
    double *P;
    gsl_error_handler_t *old_handler;

    /* Names of arguments */
    static const char *keywords[] = {"gmst", "toas",
        "w_toas", "locations", "nside", NULL};

    /* Silence warning about unused parameter. */
    (void)module;

    /* Parse arguments */
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dOOO|l", keywords,
        &gmst, &toas_obj, &w_toas_obj, &locations_obj, &nside))
        goto fail;

    if (nside == -1)
    {
        npix = -1;
    } else {
        npix = nside2npix(nside);
        if (npix < 0)
        {
            PyErr_SetString(PyExc_ValueError, "nside must be a power of 2");
            goto fail;
        }
    }

    toas_npy = (PyArrayObject *) PyArray_ContiguousFromAny(toas_obj, NPY_DOUBLE, 1, 1);
    if (!toas_npy) goto fail;
    nifos = PyArray_DIM(toas_npy, 0);
    toas = PyArray_DATA(toas_npy);

    w_toas_npy = (PyArrayObject *) PyArray_ContiguousFromAny(w_toas_obj, NPY_DOUBLE, 1, 1);
    if (!w_toas_npy) goto fail;
    if (PyArray_DIM(w_toas_npy, 0) != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and w_toas must have the same length");
        goto fail;
    }
    w_toas = PyArray_DATA(w_toas_npy);

    locations_npy = malloc(nifos * sizeof(PyObject *));
    if (!locations_npy)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
        locations_npy[i] = NULL;
    locations = malloc(nifos * sizeof(double *));
    if (!locations)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }

    n = PySequence_Length(locations_obj);
    if (n < 0) goto fail;
    if (n != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and locations must have the same length");
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
    {
        PyObject *obj = PySequence_GetItem(locations_obj, i);
        if (!obj) goto fail;
        locations_npy[i] = (PyArrayObject *) PyArray_ContiguousFromAny(obj, NPY_DOUBLE, 1, 1);
        Py_XDECREF(obj);
        if (!locations_npy[i]) goto fail;
        if (PyArray_DIM(locations_npy[i], 0) != 3)
        {
            PyErr_SetString(PyExc_ValueError, "expected every element of locations to be a vector of length 3");
            goto fail;
        }
        locations[i] = PyArray_DATA(locations_npy[i]);
    }

    old_handler = gsl_set_error_handler(my_gsl_error);
    P = bayestar_sky_map_tdoa(&npix, gmst, nifos, locations, toas, w_toas);
    gsl_set_error_handler(old_handler);

    if (!P)
        goto fail;
    premalloced = premalloced_new(P);
    if (!premalloced)
        goto fail;
    dims[0] = npix;
    out = (PyArrayObject *) PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, P);
    if (!out)
        goto fail;
    Py_INCREF(premalloced);
#ifdef PyArray_BASE
    /* FIXME: PyArray_BASE changed from a macro to a getter function in
     * Numpy 1.7. When we drop Numpy 1.6 support, remove this #ifdef block. */
    PyArray_BASE(out) = premalloced;
#else
    if (PyArray_SetBaseObject(out, premalloced))
        goto fail;
#endif

    ret = out;
    out = NULL;
fail:
    Py_XDECREF(toas_npy);
    Py_XDECREF(w_toas_npy);
    if (locations_npy)
        for (i = 0; i < nifos; i ++)
            Py_XDECREF(locations_npy[i]);
    free(locations_npy);
    free(locations);
    Py_XDECREF(premalloced);
    Py_XDECREF(out);
    return (PyObject *) ret;
};


static PyObject *sky_map_tdoa_snr(PyObject *module, PyObject *args, PyObject *kwargs)
{
    long i;
    Py_ssize_t n;
    long nside = -1;
    long npix;
    long nifos = 0;
    double gmst;
    PyObject *toas_obj, *snrs_obj, *w_toas_obj, *responses_obj,
        *locations_obj, *horizons_obj;

    PyArrayObject *toas_npy = NULL, *snrs_npy = NULL, *w_toas_npy = NULL, **responses_npy = NULL, **locations_npy = NULL, *horizons_npy = NULL;

    double *toas;
    double *snrs;
    double *w_toas;
    const float **responses = NULL;
    const double **locations = NULL;
    double *horizons;

    double min_distance, max_distance;
    int prior_distance_power;

    npy_intp dims[1];
    PyArrayObject *out = NULL, *ret = NULL;
    PyObject *premalloced = NULL;
    double *P;
    gsl_error_handler_t *old_handler;

    /* Names of arguments */
    static const char *keywords[] = {"gmst", "toas", "snrs",
        "w_toas", "responses", "locations", "horizons",
        "min_distance", "max_distance", "prior_distance_power", "nside", NULL};

    /* Silence warning about unused parameter. */
    (void)module;

    /* Parse arguments */
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dOOOOOOddi|l", keywords,
        &gmst, &toas_obj, &snrs_obj, &w_toas_obj,
        &responses_obj, &locations_obj, &horizons_obj,
        &min_distance, &max_distance, &prior_distance_power, &nside)) goto fail;

    if (nside == -1)
    {
        npix = -1;
    } else {
        npix = nside2npix(nside);
        if (npix < 0)
        {
            PyErr_SetString(PyExc_ValueError, "nside must be a power of 2");
            goto fail;
        }
    }

    toas_npy = (PyArrayObject *) PyArray_ContiguousFromAny(toas_obj, NPY_DOUBLE, 1, 1);
    if (!toas_npy) goto fail;
    nifos = PyArray_DIM(toas_npy, 0);
    toas = PyArray_DATA(toas_npy);

    snrs_npy = (PyArrayObject *) PyArray_ContiguousFromAny(snrs_obj, NPY_DOUBLE, 1, 1);
    if (!snrs_npy) goto fail;
    if (PyArray_DIM(snrs_npy, 0) != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and snrs must have the same length");
        goto fail;
    }
    snrs = PyArray_DATA(snrs_npy);

    w_toas_npy = (PyArrayObject *) PyArray_ContiguousFromAny(w_toas_obj, NPY_DOUBLE, 1, 1);
    if (!w_toas_npy) goto fail;
    if (PyArray_DIM(w_toas_npy, 0) != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and w_toas must have the same length");
        goto fail;
    }
    w_toas = PyArray_DATA(w_toas_npy);

    responses_npy = malloc(nifos * sizeof(PyObject *));
    if (!responses_npy)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
        responses_npy[i] = NULL;
    responses = malloc(nifos * sizeof(float *));
    if (!responses)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }

    n = PySequence_Length(responses_obj);
    if (n < 0) goto fail;
    if (n != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and responses must have the same length");
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
    {
        PyObject *obj = PySequence_GetItem(responses_obj, i);
        if (!obj) goto fail;
        responses_npy[i] = (PyArrayObject *) PyArray_ContiguousFromAny(obj, NPY_FLOAT, 2, 2);
        Py_XDECREF(obj);
        if (!responses_npy[i]) goto fail;
        if (PyArray_DIM(responses_npy[i], 0) != 3 || PyArray_DIM(responses_npy[i], 1) != 3)
        {
            PyErr_SetString(PyExc_ValueError, "expected every element of responses to be a 3x3 matrix");
            goto fail;
        }
        responses[i] = PyArray_DATA(responses_npy[i]);
    }

    locations_npy = malloc(nifos * sizeof(PyObject *));
    if (!locations_npy)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
        locations_npy[i] = NULL;
    locations = malloc(nifos * sizeof(double *));
    if (!locations)
    {
        PyErr_SetNone(PyExc_MemoryError);
        goto fail;
    }

    n = PySequence_Length(locations_obj);
    if (n < 0) goto fail;
    if (n != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and locations must have the same length");
        goto fail;
    }
    for (i = 0; i < nifos; i ++)
    {
        PyObject *obj = PySequence_GetItem(locations_obj, i);
        if (!obj) goto fail;
        locations_npy[i] = (PyArrayObject *) PyArray_ContiguousFromAny(obj, NPY_DOUBLE, 1, 1);
        Py_XDECREF(obj);
        if (!locations_npy[i]) goto fail;
        if (PyArray_DIM(locations_npy[i], 0) != 3)
        {
            PyErr_SetString(PyExc_ValueError, "expected every element of locations to be a vector of length 3");
            goto fail;
        }
        locations[i] = PyArray_DATA(locations_npy[i]);
    }

    horizons_npy = (PyArrayObject *) PyArray_ContiguousFromAny(horizons_obj, NPY_DOUBLE, 1, 1);
    if (!horizons_npy) goto fail;
    if (PyArray_DIM(horizons_npy, 0) != nifos)
    {
        PyErr_SetString(PyExc_ValueError, "toas and horizons must have the same length");
        goto fail;
    }
    horizons = PyArray_DATA(horizons_npy);

    old_handler = gsl_set_error_handler(my_gsl_error);
    P = bayestar_sky_map_tdoa_snr(&npix, gmst, nifos, responses, locations, toas, snrs, w_toas, horizons, min_distance, max_distance, prior_distance_power);
    gsl_set_error_handler(old_handler);

    if (!P)
        goto fail;
    premalloced = premalloced_new(P);
    if (!premalloced)
        goto fail;
    dims[0] = npix;
    out = (PyArrayObject *) PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, P);
    if (!out)
        goto fail;
    Py_INCREF(premalloced);
#ifdef PyArray_BASE
    /* FIXME: PyArray_BASE changed from a macro to a getter function in
     * Numpy 1.7. When we drop Numpy 1.6 support, remove this #ifdef block. */
    PyArray_BASE(out) = premalloced;
#else
    if (PyArray_SetBaseObject(out, premalloced))
        goto fail;
#endif

    ret = out;
    out = NULL;
fail:
    Py_XDECREF(toas_npy);
    Py_XDECREF(snrs_npy);
    Py_XDECREF(w_toas_npy);
    if (responses_npy)
        for (i = 0; i < nifos; i ++)
            Py_XDECREF(responses_npy[i]);
    free(responses_npy);
    free(responses);
    if (locations_npy)
        for (i = 0; i < nifos; i ++)
            Py_XDECREF(locations_npy[i]);
    free(locations_npy);
    free(locations);
    Py_XDECREF(horizons_npy);
    Py_XDECREF(premalloced);
    Py_XDECREF(out);
    return (PyObject *) ret;
};


static PyMethodDef methods[] = {
    {"tdoa", (PyCFunction)sky_map_tdoa, METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {"tdoa_snr", (PyCFunction)sky_map_tdoa_snr, METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
initsky_map(void) {
    premalloced_type.tp_new = PyType_GenericNew;
    if (PyType_Ready(&premalloced_type) < 0)
        return;

    Py_InitModule("sky_map", methods);
    import_array();
}
