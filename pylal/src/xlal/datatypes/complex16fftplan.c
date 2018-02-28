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
#include <lal/ComplexFFT.h>
#include <complex16fftplan.h>


#define MODULE_NAME PYLAL_COMPLEX16FFTPLAN_MODULE_NAME


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


static PyObject *__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	unsigned size;
	int fwdflg;
	int measurelvl;
	pylal_COMPLEX16FFTPlan *obj;

	if(!PyArg_ParseTuple(args, "Iii", &size, &fwdflg, &measurelvl))
		return NULL;

	obj = (pylal_COMPLEX16FFTPlan *) PyType_GenericNew(type, args, kwds);
	if(!obj)
		return NULL;
	obj->owner = NULL;
	obj->plan = XLALCreateCOMPLEX16FFTPlan(size, fwdflg, measurelvl);
	/* FIXME:  check for failure of XLALCreateCOMPLEX16FFTPlan() */
	return (PyObject *) obj;
}


static void __del__(PyObject *self)
{
	pylal_COMPLEX16FFTPlan *obj = (pylal_COMPLEX16FFTPlan *) self;

	if(obj->owner)
		Py_DECREF(obj->owner);
	else
		/* we are the owner */
		XLALDestroyCOMPLEX16FFTPlan(obj->plan);

	self->ob_type->tp_free(self);
}


static PyTypeObject pylal_complex16fftplan_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_COMPLEX16FFTPlan),
	.tp_dealloc = __del__,
	.tp_doc = "COMPLEX16FFTPlan structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES,
	.tp_name = MODULE_NAME ".COMPLEX16FFTPlan",
	.tp_new = __new__
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


PyMODINIT_FUNC initcomplex16fftplan(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "Wrapper for LAL's COMPLEX16FFTPlan type.");

	/* COMPLEX16FFTPlan */
	_pylal_COMPLEX16FFTPlan_Type = &pylal_complex16fftplan_type;
	if(PyType_Ready(&pylal_COMPLEX16FFTPlan_Type) < 0)
		return;
	Py_INCREF(&pylal_COMPLEX16FFTPlan_Type);
	PyModule_AddObject(module, "COMPLEX16FFTPlan", (PyObject *) &pylal_COMPLEX16FFTPlan_Type);
}
