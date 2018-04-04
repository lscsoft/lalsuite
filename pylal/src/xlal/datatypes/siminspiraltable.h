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


#ifndef _PYLAL_XLAL_DATATYPES_SIMINSPIRALTABLE_H_
#define _PYLAL_XLAL_DATATYPES_SIMINSPIRALTABLE_H_


#include <Python.h>
#include <lal/LIGOMetadataTables.h>


#define PYLAL_SIMINSPIRALTABLE_MODULE_NAME "pylal.xlal.datatypes.siminspiraltable"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_SimInspiralTable_Type = NULL;
#define pylal_SimInspiralTable_Type (*_pylal_SimInspiralTable_Type)


typedef struct {
	PyObject_HEAD
	SimInspiralTable sim_inspiral;
	/* FIXME:  these should be incorporated into the LAL structure */
	long process_id_i;
	EventIDColumn event_id;
} pylal_SimInspiralTable;


static PyObject *pylal_siminspiraltable_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_SIMINSPIRALTABLE_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("SimInspiralTable");
	_pylal_SimInspiralTable_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_SimInspiralTable_Type);
	Py_DECREF(name);

	return module;
}


static PyObject *pylal_SimInspiralTable_new(const SimInspiralTable *row)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_SimInspiralTable *obj = (pylal_SimInspiralTable *) PyType_GenericNew(&pylal_SimInspiralTable_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj)
		return NULL;

	obj->sim_inspiral = *row;
	obj->sim_inspiral.event_id = &obj->event_id;
	obj->process_id_i = 0;
	if(row->event_id)
		obj->event_id.id = row->event_id->id;
	else
		obj->event_id.id = 0;

	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_SIMINSPIRALTABLE_H_ */
