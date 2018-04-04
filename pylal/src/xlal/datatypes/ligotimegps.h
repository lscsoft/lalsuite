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


#ifndef _PYLAL_XLAL_DATATYPES_LIGOTIMEGPS_H_
#define _PYLAL_XLAL_DATATYPES_LIGOTIMEGPS_H_


#include <Python.h>
#include <lal/LALDatatypes.h>
#include <lal/Date.h>


#define PYLAL_LIGOTIMEGPS_MODULE_NAME "pylal.xlal.datatypes.ligotimegps"


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyTypeObject *_pylal_LIGOTimeGPS_Type = NULL;
#define pylal_LIGOTimeGPS_Type (*_pylal_LIGOTimeGPS_Type)


typedef struct {
	PyObject_HEAD
	LIGOTimeGPS gps;
} pylal_LIGOTimeGPS;


static PyObject *pylal_ligotimegps_import(void)
{
	PyObject *name = PyString_FromString(PYLAL_LIGOTIMEGPS_MODULE_NAME);
	PyObject *module = PyImport_Import(name);
	Py_DECREF(name);

	name = PyString_FromString("LIGOTimeGPS");
	_pylal_LIGOTimeGPS_Type = (PyTypeObject *) PyDict_GetItem(PyModule_GetDict(module), name);
	Py_INCREF(&pylal_LIGOTimeGPS_Type);
	Py_DECREF(name);

	return module;
}


static PyObject *pylal_LIGOTimeGPS_new(LIGOTimeGPS gps)
{
	pylal_LIGOTimeGPS *obj = (pylal_LIGOTimeGPS *) _PyObject_New(&pylal_LIGOTimeGPS_Type);

	XLALGPSSet(&obj->gps, gps.gpsSeconds, gps.gpsNanoSeconds);

	return (PyObject *) obj;
}


#endif /* _PYLAL_XLAL_DATATYPES_LIGOTIMEGPS_H_ */
