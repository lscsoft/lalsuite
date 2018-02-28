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


#ifndef _PYLAL_XLAL_MISC_H_
#define _PYLAL_XLAL_MISC_H_


#include <Python.h>
#include <lal/XLALError.h>


/*
 * ============================================================================
 *
 *                         Attribute Get/Set Helpers
 *
 * ============================================================================
 */


struct pylal_inline_string_description {
#if PY_VERSION_HEX < 0x02050000
	int offset;
	int length;
#else
	Py_ssize_t offset;
	Py_ssize_t length;
#endif
};


PyObject *pylal_inline_string_get(PyObject *obj, void *data);
int pylal_inline_string_set(PyObject *obj, PyObject *val, void *data);


struct pylal_ilwdchar_id_description {
#if PY_VERSION_HEX < 0x02050000
	int offset;
#else
	Py_ssize_t offset;
#endif
	PyObject **id_type;
};


int pylal_ilwdchar_id_set(PyObject *obj, PyObject *val, void *data);
PyObject *pylal_ilwdchar_id_get(PyObject *obj, void *data);


/*
 * ============================================================================
 *
 *                            XLAL Error Handling
 *
 * ============================================================================
 */


PyObject *pylal_exception_from_errno(enum XLALErrorValue code, const char **msg);
void pylal_set_exception_from_xlalerrno(void);


/*
 * ============================================================================
 *
 *                            ILWD_CHAR Utilities
 *
 * ============================================================================
 */


PyObject *pylal_get_ilwdchar_class(const char *table_name, const char *column_name);


#endif /* _PYLAL_XLAL_MISC_H_ */
