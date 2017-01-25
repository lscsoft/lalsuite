/*
 * Copyright (C) 2017  Leo Singer
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

#include "config.h"
#include <lal/bayestar_moc.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>
#include "six.h"


static void capsule_free(PyObject *self)
{
    free(PyCapsule_GetPointer(self, NULL));
}


static PyObject *rasterize(PyObject *NPY_UNUSED(module), PyObject *arg)
{
    PyArrayObject *arr = (PyArrayObject *) PyArray_FromAny(
        arg, NULL, 1, 1, NPY_ARRAY_CARRAY_RO, NULL);
    PyObject *uniq_key = PyUnicode_FromString("UNIQ");
    PyObject *new_fields = PyDict_New();
    PyObject *capsule = NULL;
    PyArrayObject *ret = NULL;

    if (!arr || !uniq_key || !new_fields)
        goto done;

    PyObject *fields = PyArray_DESCR(arr)->fields;
    if (!fields)
    {
        PyErr_SetString(PyExc_ValueError, "expected record array");
        goto done;
    }

    PyObject *uniq_field = PyDict_GetItem(fields, uniq_key);
    if (!uniq_field)
    {
        PyErr_SetString(PyExc_ValueError, "array does not have UNIQ column");
        goto done;
    }

    PyObject *uniq_dtype;
    int uniq_offset;
    if (!PyArg_ParseTuple(uniq_field, "Oi", &uniq_dtype, &uniq_offset))
        return NULL;

    if (!PyArray_DescrCheck(uniq_dtype))
    {
        PyErr_SetString(PyExc_ValueError, "not a dtype");
        goto done;
    }

    int uniq_typenum = ((PyArray_Descr *) uniq_dtype)->type_num;
    if (!PyArray_EquivTypenums(uniq_typenum, NPY_UINT64))
    {
        PyErr_SetString(PyExc_ValueError, "'uniq' field must be uint64");
        goto done;
    }

    if (uniq_offset != 0)
    {
        PyErr_SetString(PyExc_ValueError, "'uniq' field must be at offset 0");
        goto done;
    }

    const void *pixels = PyArray_DATA(arr);
    const size_t offset = ((PyArray_Descr *) uniq_dtype)->elsize;
    const size_t itemsize = PyArray_ITEMSIZE(arr) - offset;
    const size_t len = PyArray_SIZE(arr);
    size_t npix;

    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(fields, &pos, &key, &value))
    {
        PyObject *new_descr;
        int new_offset;

        if (PyObject_RichCompareBool(uniq_key, key, Py_EQ))
            continue;

        if (!PyArg_ParseTuple(value, "Oi", &new_descr, &new_offset))
            goto done;

        new_offset -= offset;
        PyObject *new_value = Py_BuildValue("Oi", new_descr, new_offset);
        if (!new_value)
            goto done;
        if (PyDict_SetItem(new_fields, key, new_value))
            goto done;
    }

    void *out;
    Py_BEGIN_ALLOW_THREADS
    out = moc_rasterize64(pixels, offset, itemsize, len, &npix);
    Py_END_ALLOW_THREADS
    if (!out)
        goto done;

    /* Prepare output object */
    capsule = PyCapsule_New(out, NULL, capsule_free);
    if (!capsule)
        goto done;

    PyArray_Descr *descr;
    if (PyArray_DescrConverter(new_fields, &descr) != NPY_SUCCEED)
        goto done;

    npy_intp dims[] = {npix};
    ret = (PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, NULL,
        out, NPY_ARRAY_DEFAULT, NULL);
    if (!ret)
        goto done;

    if (PyArray_SetBaseObject(ret, capsule))
    {
        Py_DECREF(ret);
        ret = NULL;
    } else {
        capsule = NULL;
    }

done:
    Py_XDECREF(arr);
    Py_XDECREF(uniq_key);
    Py_XDECREF(new_fields);
    Py_XDECREF(capsule);
    return ret;
}


static void nest2uniq_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        *(int64_t *) &args[2][i * steps[2]] = nest2uniq64(
        *(int8_t *)  &args[0][i * steps[0]],
        *(int64_t *) &args[1][i * steps[1]]);
    }
}


static void uniq2nest_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        *(int64_t *) &args[2][i * steps[2]] = *(int64_t *) &args[0][i * steps[0]];
        *(int8_t *)  &args[1][i * steps[1]] = uniq2nest64(
         (int64_t *) &args[2][i * steps[2]]);
    }
}


static void uniq2order_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        *(int8_t *)  &args[1][i * steps[1]] = uniq2order64(
        *(int64_t *) &args[0][i * steps[0]]);
    }
}


static void uniq2pixarea_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        *(double *)  &args[1][i * steps[1]] = uniq2pixarea64(
        *(int64_t *) &args[0][i * steps[0]]);
    }
}


static const PyUFuncGenericFunction
    nest2uniq_loops[] = {nest2uniq_loop},
    uniq2nest_loops[] = {uniq2nest_loop},
    uniq2order_loops[] = {uniq2order_loop},
    uniq2pixarea_loops[] = {uniq2pixarea_loop};

static const char nest2uniq_types[] = {NPY_INT8, NPY_UINT64, NPY_UINT64},
                  uniq2nest_types[] = {NPY_UINT64, NPY_INT8, NPY_UINT64},
                  uniq2order_types[] = {NPY_UINT64, NPY_INT8},
                  uniq2pixarea_types[] = {NPY_UINT64, NPY_DOUBLE};

static const void *no_ufunc_data[] = {NULL};

static const char modulename[] = "_moc";


static PyMethodDef methods[] = {
    {"rasterize", rasterize, METH_O, "fill me in"},
    {NULL, NULL, 0, NULL}
};


static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    modulename, NULL, -1, methods,
    NULL, NULL, NULL, NULL
};


PyMODINIT_FUNC PyInit__moc(void); /* Silence -Wmissing-prototypes */
PyMODINIT_FUNC PyInit__moc(void)
{
    PyObject *module = NULL;

    gsl_set_error_handler_off();
    import_array();
    import_umath();

    module = PyModule_Create(&moduledef);
    if (!module)
        goto done;

    PyModule_AddObject(
        module, "nest2uniq", PyUFunc_FromFuncAndData(
            nest2uniq_loops, no_ufunc_data,
            nest2uniq_types, 1, 2, 1, PyUFunc_None,
            "nest2uniq", NULL, 0));

    PyModule_AddObject(
        module, "uniq2nest", PyUFunc_FromFuncAndData(
            uniq2nest_loops, no_ufunc_data,
            uniq2nest_types, 1, 1, 2, PyUFunc_None,
            "uniq2nest", NULL, 0));

    PyModule_AddObject(
        module, "uniq2order", PyUFunc_FromFuncAndData(
            uniq2order_loops, no_ufunc_data,
            uniq2order_types, 1, 1, 1, PyUFunc_None,
            "uniq2order", NULL, 0));

    PyModule_AddObject(
        module, "uniq2pixarea", PyUFunc_FromFuncAndData(
            uniq2pixarea_loops, no_ufunc_data,
            uniq2pixarea_types, 1, 1, 1, PyUFunc_None,
            "uniq2pixarea", NULL, 0));

done:
    return module;
}


SIX_COMPAT_MODULE(_moc)
