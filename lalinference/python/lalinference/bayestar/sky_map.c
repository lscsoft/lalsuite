/*
 * Copyright (C) 2013-2016  Leo Singer
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

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>
#include <chealpix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_nan.h>
#include <lal/bayestar_sky_map.h>
#include <lal/bayestar_distance.h>
#include <assert.h>


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
    Py_TYPE(self)->tp_free((PyObject *) self);
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


static PyObject *premalloced_npy_double_array(double *data, int ndims, npy_intp *dims)
{
    PyObject *premalloced = premalloced_new(data);
    if (!premalloced)
        return NULL;

    PyArrayObject *out = (PyArrayObject *)
        PyArray_SimpleNewFromData(ndims, dims, NPY_DOUBLE, data);
    if (!out)
    {
        Py_DECREF(premalloced);
        return NULL;
    }

#ifdef PyArray_BASE
    /* FIXME: PyArray_BASE changed from a macro to a getter function in
     * Numpy 1.7. When we drop Numpy 1.6 support, remove this #ifdef block. */
    PyArray_BASE(out) = premalloced;
#else
    if (PyArray_SetBaseObject(out, premalloced))
    {
        Py_DECREF(out);
        return NULL;
    }
#endif

    return (PyObject *)out;
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


#define INPUT_LIST_OF_ARRAYS(NAME, NPYTYPE, DEPTH, CHECK) \
{ \
    const Py_ssize_t n = PySequence_Length(NAME##_obj); \
    if (n < 0) \
        goto fail; \
    else if (n != nifos) { \
        PyErr_SetString(PyExc_ValueError, #NAME \
        " appears to be the wrong length for the number of detectors"); \
        goto fail; \
    } \
    for (unsigned int iifo = 0; iifo < nifos; iifo ++) \
    { \
        PyObject *obj = PySequence_ITEM(NAME##_obj, iifo); \
        if (!obj) goto fail; \
        PyArrayObject *npy = (PyArrayObject *) \
            PyArray_ContiguousFromAny(obj, NPYTYPE, DEPTH, DEPTH); \
        Py_XDECREF(obj); \
        if (!npy) goto fail; \
        NAME##_npy[iifo] = npy; \
        { CHECK } \
        NAME[iifo] = PyArray_DATA(npy); \
    } \
}

#define FREE_INPUT_LIST_OF_ARRAYS(NAME) \
{ \
    for (unsigned int iifo = 0; iifo < nifos; iifo ++) \
        Py_XDECREF(NAME##_npy[iifo]); \
}

#define INPUT_VECTOR_NIFOS(CTYPE, NAME, NPYTYPE) \
    NAME##_npy = (PyArrayObject *) \
        PyArray_ContiguousFromAny(NAME##_obj, NPYTYPE, 1, 1); \
    if (!NAME##_npy) goto fail; \
    if (PyArray_DIM(NAME##_npy, 0) != nifos) \
    { \
        PyErr_SetString(PyExc_ValueError, #NAME \
            " appears to be the wrong length for the number of detectors"); \
        goto fail; \
    } \
    const CTYPE *NAME = PyArray_DATA(NAME##_npy);

#define INPUT_VECTOR_DOUBLE_NIFOS(NAME) \
    INPUT_VECTOR_NIFOS(double, NAME, NPY_DOUBLE)


static PyObject *sky_map_toa_phoa_snr(
    PyObject *NPY_UNUSED(module), PyObject *args, PyObject *kwargs)
{
    /* Input arguments */
    long nside = -1;
    long npix;
    double min_distance;
    double max_distance;
    int prior_distance_power;
    double gmst;
    unsigned int nifos;
    unsigned long nsamples = 0;
    double sample_rate;
    PyObject *acors_obj;
    PyObject *responses_obj;
    PyObject *locations_obj;
    PyObject *horizons_obj;
    PyObject *toas_obj;
    PyObject *phoas_obj;
    PyObject *snrs_obj;

    /* Names of arguments */
    static const char *keywords[] = {"min_distance", "max_distance",
        "prior_distance_power", "gmst", "sample_rate", "acors", "responses",
        "locations", "horizons", "toas", "phoas", "snrs", "nside", NULL};

    /* Parse arguments */
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ddiddOOOOOOO|l",
        keywords, &min_distance, &max_distance, &prior_distance_power, &gmst,
        &sample_rate, &acors_obj, &responses_obj, &locations_obj,
        &horizons_obj, &toas_obj, &phoas_obj, &snrs_obj, &nside)) return NULL;

    /* Determine HEALPix resolution, if specified */
    if (nside == -1)
    {
        npix = -1;
    } else {
        npix = nside2npix(nside);
        if (npix == -1)
        {
            PyErr_SetString(PyExc_ValueError, "nside must be a power of 2");
            return NULL;
        }
    }

    /* Determine number of detectors */
    {
        Py_ssize_t n = PySequence_Length(acors_obj);
        if (n < 0) return NULL;
        nifos = n;
    }

    /* Return value */
    PyObject *out = NULL;

    /* Numpy array objects */
    PyArrayObject *acors_npy[nifos], *responses_npy[nifos],
        *locations_npy[nifos], *horizons_npy = NULL, *toas_npy = NULL,
        *phoas_npy = NULL, *snrs_npy = NULL;
    memset(acors_npy, 0, sizeof(acors_npy));
    memset(responses_npy, 0, sizeof(responses_npy));
    memset(locations_npy, 0, sizeof(locations_npy));

    /* Arrays of pointers for inputs with multiple dimensions */
    const double complex *acors[nifos];
    const float (*responses[nifos])[3];
    const double *locations[nifos];

    /* Gather C-aligned arrays from Numpy types */
    INPUT_LIST_OF_ARRAYS(acors, NPY_CDOUBLE, 1,
        npy_intp dim = PyArray_DIM(npy, 0);
        if (iifo == 0)
            nsamples = dim;
        else if ((unsigned long)dim != nsamples)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of acors to be vectors of the same length");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(responses, NPY_FLOAT, 2,
        if (PyArray_DIM(npy, 0) != 3 || PyArray_DIM(npy, 1) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of responses to be 3x3 arrays");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(locations, NPY_DOUBLE, 1,
        if (PyArray_DIM(npy, 0) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of locations to be vectors of length 3");
            goto fail;
        }
    )
    INPUT_VECTOR_DOUBLE_NIFOS(horizons)
    INPUT_VECTOR_DOUBLE_NIFOS(toas)
    INPUT_VECTOR_DOUBLE_NIFOS(phoas)
    INPUT_VECTOR_DOUBLE_NIFOS(snrs)

    /* Call function */
    gsl_error_handler_t *old_handler = gsl_set_error_handler(my_gsl_error);
    double (*ret)[4] = bayestar_sky_map_toa_phoa_snr(&npix, min_distance,
        max_distance, prior_distance_power, gmst, nifos, nsamples, sample_rate,
        acors, responses, locations, horizons, toas, phoas, snrs);
    gsl_set_error_handler(old_handler);

    /* Prepare output object */
    if (ret)
    {
        npy_intp dims[] = {npix, 4};
        out = premalloced_npy_double_array(ret, 2, dims);
    }

fail: /* Cleanup */
    FREE_INPUT_LIST_OF_ARRAYS(acors)
    FREE_INPUT_LIST_OF_ARRAYS(responses)
    FREE_INPUT_LIST_OF_ARRAYS(locations)
    Py_XDECREF(horizons_npy);
    Py_XDECREF(toas_npy);
    Py_XDECREF(phoas_npy);
    Py_XDECREF(snrs_npy);
    return out;
};


static PyObject *log_likelihood_toa_snr(
    PyObject *NPY_UNUSED(module), PyObject *args, PyObject *kwargs)
{
    /* Input arguments */
    double ra;
    double sin_dec;
    double distance;
    double u;
    double twopsi;
    double t;
    double gmst;
    unsigned int nifos;
    unsigned long nsamples = 0;
    double sample_rate;
    PyObject *acors_obj;
    PyObject *responses_obj;
    PyObject *locations_obj;
    PyObject *horizons_obj;
    PyObject *toas_obj;
    PyObject *snrs_obj;

    /* Names of arguments */
    static const char *keywords[] = {"params", "gmst", "sample_rate", "acors",
        "responses", "locations", "horizons", "toas", "snrs", NULL};

    /* Parse arguments */
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "(dddddd)ddOOOOOO",
        keywords, &ra, &sin_dec, &distance, &u, &twopsi, &t, &gmst,
        &sample_rate, &acors_obj, &responses_obj, &locations_obj, &horizons_obj,
        &toas_obj, &snrs_obj)) return NULL;

    /* Determine number of detectors */
    {
        Py_ssize_t n = PySequence_Length(acors_obj);
        if (n < 0) return NULL;
        nifos = n;
    }

    /* Return value */
    PyObject *out = NULL;

    /* Numpy array objects */
    PyArrayObject *acors_npy[nifos], *responses_npy[nifos],
        *locations_npy[nifos], *horizons_npy = NULL, *toas_npy = NULL,
        *snrs_npy = NULL;
    memset(acors_npy, 0, sizeof(acors_npy));
    memset(responses_npy, 0, sizeof(responses_npy));
    memset(locations_npy, 0, sizeof(locations_npy));

    /* Arrays of pointers for inputs with multiple dimensions */
    const double complex *acors[nifos];
    const float (*responses[nifos])[3];
    const double *locations[nifos];

    /* Gather C-aligned arrays from Numpy types */
    INPUT_LIST_OF_ARRAYS(acors, NPY_CDOUBLE, 1,
        npy_intp dim = PyArray_DIM(npy, 0);
        if (iifo == 0)
            nsamples = dim;
        else if ((unsigned long)dim != nsamples)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of acors to be vectors of the same length");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(responses, NPY_FLOAT, 2,
        if (PyArray_DIM(npy, 0) != 3 || PyArray_DIM(npy, 1) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of responses to be 3x3 arrays");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(locations, NPY_DOUBLE, 1,
        if (PyArray_DIM(npy, 0) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of locations to be vectors of length 3");
            goto fail;
        }
    )
    INPUT_VECTOR_DOUBLE_NIFOS(horizons)
    INPUT_VECTOR_DOUBLE_NIFOS(toas)
    INPUT_VECTOR_DOUBLE_NIFOS(snrs)

    /* Call function */
    const double ret = bayestar_log_likelihood_toa_snr(ra, sin_dec,
        distance, u, twopsi, t, gmst, nifos, nsamples, sample_rate, acors,
        responses, locations, horizons, toas, snrs);

    /* Prepare output object */
    out = PyFloat_FromDouble(ret);

fail: /* Cleanup */
    FREE_INPUT_LIST_OF_ARRAYS(acors)
    FREE_INPUT_LIST_OF_ARRAYS(responses)
    FREE_INPUT_LIST_OF_ARRAYS(locations)
    Py_XDECREF(horizons_npy);
    Py_XDECREF(toas_npy);
    Py_XDECREF(snrs_npy);
    return out;
};


static PyObject *log_likelihood_toa_phoa_snr(
    PyObject *NPY_UNUSED(module), PyObject *args, PyObject *kwargs)
{
    /* Input arguments */
    double ra;
    double sin_dec;
    double distance;
    double u;
    double twopsi;
    double t;
    double gmst;
    unsigned int nifos;
    unsigned long nsamples = 0;
    double sample_rate;
    PyObject *acors_obj;
    PyObject *responses_obj;
    PyObject *locations_obj;
    PyObject *horizons_obj;
    PyObject *toas_obj;
    PyObject *phoas_obj;
    PyObject *snrs_obj;

    /* Names of arguments */
    static const char *keywords[] = {"params", "gmst", "sample_rate", "acors",
        "responses", "locations", "horizons", "toas", "phoas", "snrs", NULL};

    /* Parse arguments */
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "(dddddd)ddOOOOOOO",
        keywords, &ra, &sin_dec, &distance, &u, &twopsi, &t, &gmst,
        &sample_rate, &acors_obj, &responses_obj, &locations_obj, &horizons_obj,
        &toas_obj, &phoas_obj, &snrs_obj)) return NULL;

    /* Determine number of detectors */
    {
        Py_ssize_t n = PySequence_Length(acors_obj);
        if (n < 0) return NULL;
        nifos = n;
    }

    /* Return value */
    PyObject *out = NULL;

    /* Numpy array objects */
    PyArrayObject *acors_npy[nifos], *responses_npy[nifos],
        *locations_npy[nifos], *horizons_npy = NULL, *toas_npy = NULL,
        *phoas_npy = NULL, *snrs_npy = NULL;
    memset(acors_npy, 0, sizeof(acors_npy));
    memset(responses_npy, 0, sizeof(responses_npy));
    memset(locations_npy, 0, sizeof(locations_npy));

    /* Arrays of pointers for inputs with multiple dimensions */
    const double complex *acors[nifos];
    const float (*responses[nifos])[3];
    const double *locations[nifos];

    /* Gather C-aligned arrays from Numpy types */
    INPUT_LIST_OF_ARRAYS(acors, NPY_CDOUBLE, 1,
        npy_intp dim = PyArray_DIM(npy, 0);
        if (iifo == 0)
            nsamples = dim;
        else if ((unsigned long)dim != nsamples)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of acors to be vectors of the same length");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(responses, NPY_FLOAT, 2,
        if (PyArray_DIM(npy, 0) != 3 || PyArray_DIM(npy, 1) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of responses to be 3x3 arrays");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(locations, NPY_DOUBLE, 1,
        if (PyArray_DIM(npy, 0) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of locations to be vectors of length 3");
            goto fail;
        }
    )
    INPUT_VECTOR_DOUBLE_NIFOS(horizons)
    INPUT_VECTOR_DOUBLE_NIFOS(toas)
    INPUT_VECTOR_DOUBLE_NIFOS(phoas)
    INPUT_VECTOR_DOUBLE_NIFOS(snrs)

    /* Call function */
    const double ret = bayestar_log_likelihood_toa_phoa_snr(ra, sin_dec,
        distance, u, twopsi, t, gmst, nifos, nsamples, sample_rate, acors,
        responses, locations, horizons, toas, phoas, snrs);

    /* Prepare output object */
    out = PyFloat_FromDouble(ret);

fail: /* Cleanup */
    FREE_INPUT_LIST_OF_ARRAYS(acors)
    FREE_INPUT_LIST_OF_ARRAYS(responses)
    FREE_INPUT_LIST_OF_ARRAYS(locations)
    Py_XDECREF(horizons_npy);
    Py_XDECREF(toas_npy);
    Py_XDECREF(phoas_npy);
    Py_XDECREF(snrs_npy);
    return out;
};


static void distance_moments_to_parameters_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        bayestar_distance_moments_to_parameters(
            *(double *) &args[0][i * steps[0]],
            *(double *) &args[1][i * steps[1]],
             (double *) &args[2][i * steps[2]],
             (double *) &args[3][i * steps[3]],
             (double *) &args[4][i * steps[4]]);
    }
}


static void distance_parameters_to_moments_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        bayestar_distance_parameters_to_moments(
            *(double *) &args[0][i * steps[0]],
            *(double *) &args[1][i * steps[1]],
             (double *) &args[2][i * steps[2]],
             (double *) &args[3][i * steps[3]],
             (double *) &args[4][i * steps[4]]);
    }
}


static const PyUFuncGenericFunction distance_moments_to_parameters_loops[] = {
    distance_moments_to_parameters_loop};
static const PyUFuncGenericFunction distance_parameters_to_moments_loops[] = {
    distance_parameters_to_moments_loop};
static const char distance_ufunc_types[] = {
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};
static const void *distance_ufunc_data[] = {NULL};


static void volume_render_kernel_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];
    const long nside = npix2nside(dimensions[2]);

    /* Assert that rotation matrix is 3x3 */
    assert(dimensions[1] == 3);

    /* FIXME: Check that array arguments are stored contiguously */

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        *(double *) &args[11][i * steps[11]] = bayestar_volume_render_kernel(
            *(double *)   &args[0][i * steps[0]],
            *(double *)   &args[1][i * steps[1]],
            *(double *)   &args[2][i * steps[2]],
            *(int *)      &args[3][i * steps[3]],
            *(int *)      &args[4][i * steps[4]],
             (double *)   &args[5][i * steps[5]],
            nside,
            *(npy_bool *) &args[6][i * steps[6]],
             (double *)   &args[7][i * steps[7]],
             (double *)   &args[8][i * steps[8]],
             (double *)   &args[9][i * steps[9]],
             (double *)   &args[10][i * steps[10]]);
    }
}


static const PyUFuncGenericFunction volume_render_kernel_loops[] = {
    volume_render_kernel_loop};
static const char volume_render_kernel_ufunc_types[] = {
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_INT, NPY_INT, NPY_DOUBLE, NPY_BOOL, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};
static const void *volume_render_kernel_ufunc_data[] = {NULL};


static void marginal_distance_distribution_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];
    const long npix = dimensions[1];

    /* Assert that array arguments are stored contiguously */
    assert(steps[6] == sizeof(double));

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        *(double *) &args[5][i * steps[5]] =
            bayestar_marginal_distance_distribution(
            *(double *) &args[0][i * steps[0]], npix,
             (double *) &args[1][i * steps[1]],
             (double *) &args[2][i * steps[2]],
             (double *) &args[3][i * steps[3]],
             (double *) &args[4][i * steps[4]]);
    }
}


static const PyUFuncGenericFunction marginal_distance_distribution_loops[] = {
    marginal_distance_distribution_loop};
static const char marginal_distance_distribution_ufunc_types[] = {
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};
static const void *marginal_distance_distribution_ufunc_data[] = {NULL};


static PyObject *test(
    PyObject *NPY_UNUSED(module), PyObject *NPY_UNUSED(arg))
{
    return PyLong_FromLong(bayestar_test());
}


static PyMethodDef methods[] = {
    {"toa_phoa_snr", (PyCFunction)sky_map_toa_phoa_snr,
        METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {"log_likelihood_toa_snr", (PyCFunction)log_likelihood_toa_snr,
        METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {"log_likelihood_toa_phoa_snr", (PyCFunction)log_likelihood_toa_phoa_snr,
        METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {"test", (PyCFunction)test,
        METH_NOARGS, "fill me in"},
    {NULL, NULL, 0, NULL}
};


typedef struct {
    PyObject_HEAD
    log_radial_integrator *integrator;
} LogRadialIntegrator;


static void LogRadialIntegrator_dealloc(LogRadialIntegrator *self)
{
    log_radial_integrator_free(self->integrator);
    Py_TYPE(self)->tp_free((PyObject *)self);
}


static int LogRadialIntegrator_init(
    LogRadialIntegrator *self, PyObject *args, PyObject *kwargs)
{
    static const char *keywords[] = {"r1", "r2", "k", "pmax", "size", NULL};
    self->integrator = NULL;
    double r1, r2, pmax;
    int k, size;

    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "ddidi", keywords, &r1, &r2, &k, &pmax, &size))
        return -1;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(my_gsl_error);
    self->integrator = log_radial_integrator_init(r1, r2, k, pmax, size);
    gsl_set_error_handler(old_handler);

    if (self->integrator)
        return 0;
    else
        return -1;
}


static PyObject *LogRadialIntegrator_call(
    LogRadialIntegrator *self, PyObject *args, PyObject *kwargs)
{
    static const char *keywords[] = {"p", "b", NULL};
    double p, b, result;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dd", keywords, &p, &b))
        return NULL;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(my_gsl_error);
    result = log_radial_integrator_eval(self->integrator, p, b);
    gsl_set_error_handler(old_handler);

    if (PyErr_Occurred())
        return NULL;

    return PyFloat_FromDouble(result);
}


static PyTypeObject LogRadialIntegrator_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "sky_map.LogRadialIntegrator",          /*tp_name*/
    sizeof(LogRadialIntegrator),            /*tp_basicsize*/
    0,                                      /*tp_itemsize*/
    (destructor)LogRadialIntegrator_dealloc,/*tp_dealloc*/
    0,                                      /*tp_print*/
    0,                                      /*tp_getattr*/
    0,                                      /*tp_setattr*/
    0,                                      /*tp_compare*/
    0,                                      /*tp_repr*/
    0,                                      /*tp_as_number*/
    0,                                      /*tp_as_sequence*/
    0,                                      /*tp_as_mapping*/
    0,                                      /*tp_hash */
    (ternaryfunc)LogRadialIntegrator_call,  /*tp_call*/
    0,                                      /*tp_str*/
    0,                                      /*tp_getattro*/
    0,                                      /*tp_setattro*/
    0,                                      /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,                     /*tp_flags*/
    "fill me in",                           /* tp_doc */
    0,                                      /* tp_traverse */
    0,                                      /* tp_clear */
    0,                                      /* tp_richcompare */
    0,                                      /* tp_weaklistoffset */
    0,                                      /* tp_iter */
    0,                                      /* tp_iternext */
    0,                                      /* tp_methods */
    0,                                      /* tp_members */
    0,                                      /* tp_getset */
    0,                                      /* tp_base */
    0,                                      /* tp_dict */
    0,                                      /* tp_descr_get */
    0,                                      /* tp_descr_set */
    0,                                      /* tp_dictoffset */
    (initproc)LogRadialIntegrator_init,     /* tp_init */
};


static const char modulename[] = "sky_map";


#if PY_MAJOR_VERSION >= 3
static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    modulename, NULL, -1, methods
};
#endif


#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initsky_map(void); /* Silence -Wmissing-prototypes */
PyMODINIT_FUNC initsky_map(void)
#else
PyMODINIT_FUNC PyInit_sky_map(void); /* Silence -Wmissing-prototypes */
PyMODINIT_FUNC PyInit_sky_map(void)
#endif
{
    PyObject *module;
    import_array();
    import_umath();

    premalloced_type.tp_new = PyType_GenericNew;
    if (PyType_Ready(&premalloced_type) < 0)
    {
#if PY_MAJOR_VERSION < 3
        return;
#else
        return NULL;
#endif
    }

    LogRadialIntegrator_type.tp_new = PyType_GenericNew;
    if (PyType_Ready(&LogRadialIntegrator_type) < 0)
    {
#if PY_MAJOR_VERSION < 3
        return;
#else
        return NULL;
#endif
    }

#if PY_MAJOR_VERSION < 3
    module = Py_InitModule(modulename, methods);
#else
    module = PyModule_Create(&moduledef);
#endif

    PyModule_AddObject(
        module, "distance_moments_to_parameters", PyUFunc_FromFuncAndData(
            distance_moments_to_parameters_loops, distance_ufunc_data,
            distance_ufunc_types, 1, 2, 3, PyUFunc_None,
            "distance_moments_to_parameters", NULL, 0));

    PyModule_AddObject(
        module, "distance_parameters_to_moments", PyUFunc_FromFuncAndData(
            distance_parameters_to_moments_loops, distance_ufunc_data,
            distance_ufunc_types, 1, 2, 3, PyUFunc_None,
            "distance_parameters_to_moments", NULL, 0));

    PyModule_AddObject(
        module, "volume_render_kernel", PyUFunc_FromFuncAndDataAndSignature(
            volume_render_kernel_loops, volume_render_kernel_ufunc_data,
            volume_render_kernel_ufunc_types, 1, 11, 1, PyUFunc_None,
            "volume_render_kernel", NULL, 0, "(),(),(),(),(),(i,i),(),(n),(n),(n),(n)->()"));

    PyModule_AddObject(
        module, "marginal_distance_distribution", PyUFunc_FromFuncAndDataAndSignature(
            marginal_distance_distribution_loops, marginal_distance_distribution_ufunc_data,
            marginal_distance_distribution_ufunc_types, 1, 5, 1, PyUFunc_None,
            "marginal_distance_distribution", NULL, 0, "(),(n),(n),(n),(n)->()"));

    Py_INCREF((PyObject *)&LogRadialIntegrator_type);
    PyModule_AddObject(
        module, "LogRadialIntegrator", (PyObject *)&LogRadialIntegrator_type);

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
