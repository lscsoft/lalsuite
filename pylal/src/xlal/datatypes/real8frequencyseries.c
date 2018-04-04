/*
 * Copyright (C) 2009  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <Python.h>
#include <numpy/arrayobject.h>
#include <lal/LALDatatypes.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/Units.h>
#include <lalunit.h>
#include <ligotimegps.h>
#include <real8frequencyseries.h>


#define MODULE_NAME PYLAL_REAL8FREQUENCYSERIES_MODULE_NAME


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


/*
 * Methods
 */


static PyObject *__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	pylal_REAL8FrequencySeries *obj;
	LIGOTimeGPS zero = {0, 0};

	obj = (pylal_REAL8FrequencySeries *) PyType_GenericNew(type, args, kwds);
	if(!obj)
		return NULL;
	obj->series = XLALCreateREAL8FrequencySeries(NULL, &zero, 0.0, 0.0, &lalDimensionlessUnit, 0);
	obj->owner = NULL;
	return (PyObject *) obj;
}


static int __init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	PyObject *key, *value;
	Py_ssize_t pos = 0;

	if(kwds)
		while(PyDict_Next(kwds, &pos, &key, &value))
			if(PyObject_SetAttr(self, key, value) < 0)
				return -1;
	return 0;
}


static void __del__(PyObject *self)
{
	pylal_REAL8FrequencySeries *obj = (pylal_REAL8FrequencySeries *) self;
	if(obj->owner)
		Py_DECREF(obj->owner);
	else
		/* we are the owner */
		XLALDestroyREAL8FrequencySeries(obj->series);
	self->ob_type->tp_free(self);
}


static PyObject *__getattro__(PyObject *self, PyObject *attr_name)
{
	const char *name = PyString_AsString(attr_name);
	pylal_REAL8FrequencySeries *obj = (pylal_REAL8FrequencySeries *) self;

	if(!strcmp(name, "name"))
		return PyString_FromString(obj->series->name);
	if(!strcmp(name, "epoch"))
		return pylal_LIGOTimeGPS_new(obj->series->epoch);
	if(!strcmp(name, "f0"))
		return PyFloat_FromDouble(obj->series->f0);
	if(!strcmp(name, "deltaF"))
		return PyFloat_FromDouble(obj->series->deltaF);
	if(!strcmp(name, "sampleUnits"))
		return pylal_LALUnit_new(0, obj->series->sampleUnits);
	if(!strcmp(name, "data")) {
		npy_intp dims[] = {obj->series->data->length};
		PyObject *array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, obj->series->data->data);
		PyObject *copy;
		if(!array)
			return NULL;
		/* incref self to prevent data from disappearing while
		 * array is still in use, and tell numpy to decref self
		 * when the array is deallocated */
		Py_INCREF(self);
		PyArray_BASE(array) = self;
		copy = PyArray_NewCopy((PyArrayObject *) array, NPY_ANYORDER);
		Py_DECREF(array);
		return copy;
	}
	PyErr_SetString(PyExc_AttributeError, name);
	return NULL;
}


static int __setattro__(PyObject *self, PyObject *attr_name, PyObject *value)
{
	const char *name = PyString_AsString(attr_name);
	pylal_REAL8FrequencySeries *obj = (pylal_REAL8FrequencySeries *) self;

	if(!strcmp(name, "name")) {
		const char *s = PyString_AsString(value);
		if(PyErr_Occurred())
			return -1;
		if(strlen(s) >= sizeof(obj->series->name)) {
			PyErr_Format(PyExc_ValueError, "name too long \"%s\"", s);
			return -1;
		}
		strcpy(obj->series->name, s);
		return 0;
	}
	if(!strcmp(name, "epoch")) {
		if(!PyObject_TypeCheck(value, &pylal_LIGOTimeGPS_Type)) {
			PyErr_SetObject(PyExc_TypeError, value);
			return -1;
		}
		obj->series->epoch = ((pylal_LIGOTimeGPS *) value)->gps;
		return 0;
	}
	if(!strcmp(name, "f0")) {
		double f0 = PyFloat_AsDouble(value);
		if(PyErr_Occurred())
			return -1;
		obj->series->f0 = f0;
		return 0;
	}
	if(!strcmp(name, "deltaF")) {
		double deltaF = PyFloat_AsDouble(value);
		if(PyErr_Occurred())
			return -1;
		obj->series->deltaF = deltaF;
		return 0;
	}
	if(!strcmp(name, "sampleUnits")) {
		if(!PyObject_TypeCheck(value, &pylal_LALUnit_Type)) {
			PyErr_SetObject(PyExc_TypeError, value);
			return -1;
		}
		obj->series->sampleUnits = ((pylal_LALUnit *) value)->unit;
		return 0;
	}
	if(!strcmp(name, "data")) {
		int n;
		/* require array of double precision floats */
		if(!PyArray_Check(value) || (PyArray_TYPE(value) != NPY_DOUBLE)) {
			PyErr_SetObject(PyExc_TypeError, value);
			return -1;
		}
		/* require exactly 1 dimension */
		if(((PyArrayObject *) value)->nd != 1) {
			PyErr_SetObject(PyExc_ValueError, value);
			return -1;
		}
		n = PyArray_DIM(value, 0);
		if(n != obj->series->data->length)
			obj->series->data = XLALResizeREAL8Sequence(obj->series->data, 0, n);
		memcpy(obj->series->data->data, PyArray_GETPTR1(value, 0), n * sizeof(*obj->series->data->data));
		return 0;
	}
	PyErr_SetString(PyExc_AttributeError, name);
	return -1;
}


/*
 * Type
 */


PyTypeObject pylal_real8frequencyseries_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_REAL8FrequencySeries),
	.tp_dealloc = __del__,
	.tp_doc = "REAL8FrequencySeries structure",
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_getattro = __getattro__,
	.tp_setattro = __setattro__,
	.tp_name = MODULE_NAME ".REAL8FrequencySeries",
	.tp_new = __new__,
	.tp_init = __init__
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


PyMODINIT_FUNC initreal8frequencyseries(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "Wrapper for LAL's REAL8FrequencySeries type.");

	import_array();
	pylal_lalunit_import();
	pylal_ligotimegps_import();

	/*
	 * REAL8FrequencySeries
	 */

	_pylal_REAL8FrequencySeries_Type = &pylal_real8frequencyseries_type;
	if(PyType_Ready(&pylal_REAL8FrequencySeries_Type) < 0)
		return;
	Py_INCREF(&pylal_REAL8FrequencySeries_Type);
	PyModule_AddObject(module, "REAL8FrequencySeries", (PyObject *) &pylal_REAL8FrequencySeries_Type);
}
