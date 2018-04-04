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


#ifndef _PYLAL_XLAL_DATATYPES_REAL8FFTPLAN_H_
#define _PYLAL_XLAL_DATATYPES_REAL8FFTPLAN_H_


#include <Python.h>
#include <lal/RealFFT.h>


#define PYLAL_REAL8FFTPLAN_MODULE_NAME "pylal.xlal.datatypes.real8fftplan"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_REAL8FFTPlan_Type = NULL;
#define pylal_REAL8FFTPlan_Type (*_pylal_REAL8FFTPlan_Type)


typedef struct {
	PyObject_HEAD
	PyObject *owner;
	REAL8FFTPlan *plan;
} pylal_REAL8FFTPlan;


static PyObject *pylal_real8fftplan_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_REAL8FFTPLAN_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("REAL8FFTPlan");
	_pylal_REAL8FFTPlan_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_REAL8FFTPlan_Type);
	Py_DECREF(name);

	return module;
}


static PyObject *pylal_REAL8FFTPlan_new(REAL8FFTPlan *plan, PyObject *owner)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_REAL8FFTPlan *obj = (pylal_REAL8FFTPlan *) PyType_GenericNew(&pylal_REAL8FFTPlan_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj) {
		if(!owner)
			XLALDestroyREAL8FFTPlan(plan);
		return NULL;
	}
	if(owner)
		Py_INCREF(owner);
	obj->owner = owner;
	XLALDestroyREAL8FFTPlan(obj->plan);
	obj->plan = plan;
	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_REAL8FFTPLAN_H_ */
