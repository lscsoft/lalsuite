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


#ifndef _PYLAL_XLAL_DATATYPES_SNGLINSPIRALTABLE_H_
#define _PYLAL_XLAL_DATATYPES_SNGLINSPIRALTABLE_H_


#include <Python.h>
#include <lal/LIGOMetadataTables.h>


#define PYLAL_SNGLINSPIRALTABLE_MODULE_NAME "pylal.xlal.datatypes.snglinspiraltable"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_SnglInspiralTable_Type = NULL;
#define pylal_SnglInspiralTable_Type (*_pylal_SnglInspiralTable_Type)


typedef struct {
	PyObject_HEAD
	SnglInspiralTable sngl_inspiral;
	/* FIXME:  these should be incorporated into the LAL structure */
	long process_id_i;
	EventIDColumn event_id;
} pylal_SnglInspiralTable;


static PyObject *pylal_snglinspiraltable_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_SNGLINSPIRALTABLE_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("SnglInspiralTable");
	_pylal_SnglInspiralTable_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_SnglInspiralTable_Type);
	Py_DECREF(name);

	return module;
}


static PyObject *pylal_SnglInspiralTable_new(const SnglInspiralTable *row)
{
	PyObject *empty_tuple = PyTuple_New(0);
	pylal_SnglInspiralTable *obj = (pylal_SnglInspiralTable *) PyType_GenericNew(&pylal_SnglInspiralTable_Type, empty_tuple, NULL);
	Py_DECREF(empty_tuple);
	if(!obj)
		return NULL;

	obj->sngl_inspiral = *row;
	obj->sngl_inspiral.event_id = &obj->event_id;
	obj->process_id_i = 0;
	if(row->event_id)
		obj->event_id.id = row->event_id->id;
	else
		obj->event_id.id = 0;

	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_SNGLINSPIRALTABLE_H_ */
