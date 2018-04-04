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


#ifndef _PYLAL_XLAL_DATATYPES_LALUNIT_H_
#define _PYLAL_XLAL_DATATYPES_LALUNIT_H_


#include <Python.h>
#include <lal/LALDatatypes.h>


#define PYLAL_LALUNIT_MODULE_NAME "pylal.xlal.datatypes.lalunit"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_LALUnit_Type = NULL;
#define pylal_LALUnit_Type (*_pylal_LALUnit_Type)


typedef struct {
	PyObject_HEAD
	LALUnit unit;
} pylal_LALUnit;


extern PyObject *pylal_LALUnitMeter;
extern PyObject *pylal_LALUnitKiloGram;
extern PyObject *pylal_LALUnitSecond;
extern PyObject *pylal_LALUnitAmpere;
extern PyObject *pylal_LALUnitKelvin;
extern PyObject *pylal_LALUnitStrain;
extern PyObject *pylal_LALUnitADCCount;


static PyObject *pylal_lalunit_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_LALUNIT_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("LALUnit");
	_pylal_LALUnit_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_LALUnit_Type);
	Py_DECREF(name);

	return module;
}


static PyObject *pylal_LALUnit_new(int power_of_ten, LALUnit unit)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_LALUnit *obj = (pylal_LALUnit *) PyType_GenericNew(&pylal_LALUnit_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	obj->unit = unit;
	obj->unit.powerOfTen += power_of_ten;
	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_LALUNIT_H_ */
