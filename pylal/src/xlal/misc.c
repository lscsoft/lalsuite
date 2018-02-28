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
 *                              Helper Utilities
 *
 * ============================================================================
 */


#include <Python.h>
#include <lal/XLALError.h>
#include <misc.h>


/*
 * ============================================================================
 *
 *                         Attribute Get/Set Helpers
 *
 * ============================================================================
 */


PyObject *pylal_inline_string_get(PyObject *obj, void *data)
{
	const struct pylal_inline_string_description *desc = data;
	char *s = (void *) obj + desc->offset;

	if(strlen(s) >= desc->length) {
		/* something's wrong, obj probably isn't a valid address */
	}

	return PyString_FromString(s);
}


int pylal_inline_string_set(PyObject *obj, PyObject *val, void *data)
{
	const struct pylal_inline_string_description *desc = data;
	char *v = PyString_AsString(val);
	char *s = (void *) obj + desc->offset;

	if(!v)
		return -1;
	if(strlen(v) >= desc->length) {
		PyErr_Format(PyExc_ValueError, "string too long \'%s\'", v);
		return -1;
	}

	strncpy(s, v, desc->length - 1);
	s[desc->length - 1] = '\0';

	return 0;
}


int pylal_ilwdchar_id_set(PyObject *obj, PyObject *val, void *data)
{
	const struct pylal_ilwdchar_id_description *desc = data;
	void *location = (void *) obj + desc->offset;
	long i = PyInt_AsLong(val);

	if(PyErr_Occurred())
		return -1;

	if((PyObject *) val->ob_type != *desc->id_type) {
		PyErr_SetObject(PyExc_TypeError, val);
		return -1;
	}

	*(long *) location = i;

	return 0;
}


PyObject *pylal_ilwdchar_id_get(PyObject *obj, void *data)
{
	const struct pylal_ilwdchar_id_description *desc = data;
	void *location = (void *) obj + desc->offset;

	return PyObject_CallFunction(*desc->id_type, "l", *(long *) location);
}


/*
 * ============================================================================
 *
 *                               Error Handling
 *
 * ============================================================================
 */


PyObject *pylal_exception_from_errno(enum XLALErrorValue code, const char **msg)
{
	if(msg)
		*msg = XLALErrorString(code);

	switch(code) {
	case XLAL_SUCCESS:
		return NULL;

	case XLAL_EIO:
		return PyExc_IOError;

	case XLAL_ENOMEM:
	case XLAL_EFAULT:
		return PyExc_MemoryError;

	case XLAL_EINVAL:
	case XLAL_EDOM:
	case XLAL_ERANGE:
	case XLAL_EBADLEN:
	case XLAL_ESIZE:
	case XLAL_EDIMS:
		return PyExc_ValueError;

	case XLAL_ETYPE:
		return PyExc_TypeError;

	case XLAL_ETIME:
	case XLAL_EFREQ:
	case XLAL_EUNIT:
	case XLAL_ENAME:
	case XLAL_EDATA:
		return PyExc_ValueError;

	case XLAL_ESYS:
		return PyExc_SystemError;

	case XLAL_EERR:
		return PyExc_RuntimeError;

	case XLAL_EFPDIV0:
		return PyExc_ZeroDivisionError;

	case XLAL_EFPINVAL:
	case XLAL_EFPOVRFLW:
	case XLAL_EFPUNDFLW:
	case XLAL_EFPINEXCT:
	case XLAL_EMAXITER:
	case XLAL_EDIVERGE:
	case XLAL_ESING:
	case XLAL_ETOL:
	case XLAL_ELOSS:
		return PyExc_FloatingPointError;

	default:
		return PyExc_Exception;
	}
}


void pylal_set_exception_from_xlalerrno(void)
{
	const char *msg;
	PyObject *exc = pylal_exception_from_errno(XLALGetBaseErrno(), &msg);
	PyErr_SetString(exc, msg);
	XLALClearErrno();
}


/*
 * ============================================================================
 *
 *                            ILWD_CHAR Utilities
 *
 * ============================================================================
 */


PyObject *pylal_get_ilwdchar_class(const char *table_name, const char *column_name)
{
	PyObject *module_name;
	PyObject *module;
	PyObject *func;
	PyObject *class;

	module_name = PyString_FromString("glue.ligolw.ilwd");
	if(!module_name)
		return NULL;

	module = PyImport_Import(module_name);
	Py_DECREF(module_name);
	if(!module)
		return NULL;

	func = PyMapping_GetItemString(PyModule_GetDict(module), "get_ilwdchar_class");
	Py_DECREF(module);
	if(!func)
		return NULL;

	class = PyObject_CallFunction(func, "ss", table_name, column_name);
	Py_DECREF(func);

	return class;
}
