/*
 * Copyright (C) 2010  Kipp Cannon
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


#ifndef _PYLAL_XLAL_DATATYPES_SIMBURST_H_
#define _PYLAL_XLAL_DATATYPES_SIMBURST_H_


#include <Python.h>
#include <lal/LIGOMetadataTables.h>


#define PYLAL_SIMBURST_MODULE_NAME "pylal.xlal.datatypes.simburst"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_SimBurst_Type = NULL;
#define pylal_SimBurst_Type (*_pylal_SimBurst_Type)


typedef struct {
	PyObject_HEAD
	SimBurst sim_burst;
} pylal_SimBurst;


static PyObject *pylal_simburst_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_SIMBURST_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("SimBurst");
	_pylal_SimBurst_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_SimBurst_Type);
	Py_DECREF(name);

	return module;
}


static PyObject *pylal_SimBurst_new(const SimBurst *row)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_SimBurst *obj = (pylal_SimBurst *) PyType_GenericNew(&pylal_SimBurst_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj)
		return NULL;

	obj->sim_burst = *row;

	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_SIMBURST_H_ */
