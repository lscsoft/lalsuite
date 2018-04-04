/*
 * Copyright (C) 2011  Kipp Cannon
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


#ifndef _PYLAL_XLAL_DATATYPES_SNGLBURST_H_
#define _PYLAL_XLAL_DATATYPES_SNGLBURST_H_


#include <Python.h>
#include <lal/LIGOMetadataTables.h>


#define PYLAL_SNGLBURST_MODULE_NAME "pylal.xlal.datatypes.snglburst"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_SnglBurst_Type = NULL;
#define pylal_SnglBurst_Type (*_pylal_SnglBurst_Type)


typedef struct {
	PyObject_HEAD
	SnglBurst sngl_burst;
} pylal_SnglBurst;


static PyObject *pylal_snglburst_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_SNGLBURST_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("SnglBurst");
	_pylal_SnglBurst_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_SnglBurst_Type);
	Py_DECREF(name);

	return module;
}


static PyObject *pylal_SnglBurst_new(const SnglBurst *row)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_SnglBurst *obj = (pylal_SnglBurst *) PyType_GenericNew(&pylal_SnglBurst_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj)
		return NULL;

	obj->sngl_burst = *row;

	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_SNGLBURST_H_ */
