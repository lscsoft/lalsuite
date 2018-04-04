/*
 *
 * Copyright (C) 2007,2009,2011-2013  Kipp C. Cannon
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
 *                   glue.ligolw._ilwd Extension Module
 *
 * ============================================================================
 */


#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <wchar.h>
#include <ilwd.h>


#define MODULE_NAME "glue.ligolw._ilwd"


/*
 * ============================================================================
 *
 *                               ilwdchar Type
 *
 * ============================================================================
 */


/*
 * "table_name" and "column_name" as Python objects for use in retrieving
 * attributes.
 */


static PyObject *table_name;
static PyObject *column_name;


/*
 * Some notes.  LIGO Light Weight XML "ilwd:char" IDs are strings of the
 * form "table:column:integer", for example "process:process_id:10".  Large
 * complex documents can have many millions of these strings, and their
 * storage represents a significant RAM burden.  At the same time, however,
 * while there can be millions of ID strings in use there might be only a
 * small number (e.g. 10 or fewer) ID prefixes in use (the table name and
 * column name part).  This C extension module implements a class that
 * tries to reduce the storage requirements of these ID strings by storing
 * only one copy of each prefix.  The book keeping associated with the
 * construction of new subclasses on the fly, and caching
 * previously-constructed classes is handled by Python code that can be
 * found elsewhere.  The code here is only the underlying C engine to make
 * the manipulation of these objects quick.
 *
 * The objects, like strings, are (considered to be) immutable.  For this
 * reason references to them can be shared, and so initialization is
 * performed in the __new__() method instead of the __init__() method as is
 * normal.  This allows references to existing objects to be returned
 * instead of unconditionally constructing new instances.
 */


/*
 * Methods
 */


static PyObject *ligolw_ilwdchar___add__(PyObject *self, PyObject *other)
{
	long delta = PyInt_AsLong(other);
	PyObject *new;

	if(PyErr_Occurred())
		/* argument is not int-like --> type error */
		return NULL;

	if(!delta) {
		/* no change in value */
		Py_INCREF(self);
		return self;
	}

	new = PyType_GenericNew(self->ob_type, NULL, NULL);

	if(new)
		((ligolw_ilwdchar *) new)->i = ((ligolw_ilwdchar *) self)->i + delta;

	return new;
}


static long ligolw_ilwdchar___hash__(PyObject *self)
{
	PyObject *tbl = PyObject_GetAttr(self, table_name);
	PyObject *col = PyObject_GetAttr(self, column_name);
	long hash;

	if(tbl && col) {
		hash = PyObject_Hash(tbl) ^ PyObject_Hash(col) ^ ((ligolw_ilwdchar *) self)->i;
		if(hash == -1)
			/* -1 is reserved for error conditions */
			hash = -2;
	} else
		hash = -1;

	Py_XDECREF(tbl);
	Py_XDECREF(col);

	return hash;
}


static PyObject *ligolw_ilwdchar___int__(PyObject *self)
{
	return PyInt_FromLong(((ligolw_ilwdchar *) self)->i);
}


static PyObject *ligolw_ilwdchar___new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	/* call the generic __new__() */
	PyObject *new = PyType_GenericNew(type, NULL, NULL);
	PyObject *obj;

	if(!new)
		return NULL;

	/* initialize to default value */
	((ligolw_ilwdchar *) new)->i = 0;

	if(PyArg_ParseTuple(args, "U", &obj)) {
		/* we've been passed a unicode string, see if we can parse
		 * it */
		Py_ssize_t len = PyUnicode_GetSize(obj);
		wchar_t s[len + 1];
		int converted_len = -1;
		wchar_t *tbl_name = NULL, *col_name = NULL;

		PyUnicode_AsWideChar((PyUnicodeObject *) obj, s, len);
		s[len] = 0;

		/* can we parse it as an ilwd:char string? */
		swscanf(s, L" %ml[^:]:%ml[^:]:%zu %n", &tbl_name, &col_name, &((ligolw_ilwdchar *) new)->i, &converted_len);
		if(converted_len < len) {
			/* nope, how 'bout just an int? */
			converted_len = -1;
			swscanf(s, L" %zu %n", &((ligolw_ilwdchar *) new)->i, &converted_len);
			if(converted_len < len) {
				/* nope */
#if PY_MAJOR_VERSION < 3
				PyObject *str = PyUnicode_AsEncodedString(obj, NULL, NULL);
				PyErr_Format(PyExc_TypeError, "invalid literal for ilwdchar(): '%s'", PyString_AS_STRING(str));
				Py_DECREF(str);
#else
				PyErr_Format(PyExc_ValueError, "invalid literal for ilwdchar(): '%u'", obj);
#endif
				Py_DECREF(new);
				new = NULL;
			}
		} else {
			/* yes, it's an ilwd:char string, so confirm
			 * that the table and column names are
			 * correct */
			PyObject *tbl, *col;
			PyObject *tbl_attr, *col_attr;

			tbl = PyUnicode_FromWideChar(tbl_name, wcslen(tbl_name));
			col = PyUnicode_FromWideChar(col_name, wcslen(col_name));
			tbl_attr = PyObject_GetAttr(new, table_name);
			col_attr = PyObject_GetAttr(new, column_name);

			if(!tbl || !col || !tbl_attr || !col_attr ||
			PyUnicode_Compare(tbl_attr, tbl) ||
			PyUnicode_Compare(col_attr, col)) {
				/* mismatch */
#if PY_MAJOR_VERSION < 3
				PyObject *str = PyUnicode_AsEncodedString(obj, NULL, NULL);
				PyErr_Format(PyExc_TypeError, "ilwdchar type mismatch: '%s'", PyString_AS_STRING(str));
				Py_DECREF(str);
#else
				PyErr_Format(PyExc_TypeError, "ilwdchar type mismatch: '%u'", obj);
#endif
				Py_DECREF(new);
				new = NULL;
			}

			Py_XDECREF(tbl_attr);
			Py_XDECREF(col_attr);
			Py_XDECREF(tbl);
			Py_XDECREF(col);
		}

		free(tbl_name);
		free(col_name);
	} else if(PyArg_ParseTuple(args, "|l", &((ligolw_ilwdchar *) new)->i)) {
		/* we were passed nothing or an int:  i has either been set
		 * from the int, or we'll use the default value of 0.
		 * clear the error from the earlier ParseTuple failures */
		PyErr_Clear();
	} else if(PyArg_ParseTuple(args, "O!", type, &obj)) {
		/* we've been passed an instance of our own type, just
		 * incref and return it.  clear the error from the earlier
		 * ParseTuple failures. */
		PyErr_Clear();
		Py_DECREF(new);
		new = obj;
		Py_INCREF(new);
	} else {
		/* we weren't passed a unicode or an int or an ilwdchar
		 * instance:  if args contains exactly 1 object then this
		 * is a type error, otherwise it's a wrong number of
		 * arguments error */
		Py_DECREF(new);
		if(PyArg_ParseTuple(args, "O", &new))
			PyErr_SetObject(PyExc_TypeError, new);
		new = NULL;
	}

	return new;
}


static PyObject *ligolw_ilwdchar___richcompare__(PyObject *self, PyObject *other, int op)
{
	PyObject *tbl_s, *col_s;
	PyObject *tbl_o, *col_o;
	PyObject *result;

	switch(PyObject_IsInstance(other, (PyObject *) &ligolw_ilwdchar_Type)) {
	case -1:
		/* function call failed, it will have set the exception */
		return NULL;

	case 0:
		/* type mismatch.  if other is a string or unicode raise
		 * TypeError to prevent bugs:  people will assume that the
		 * comparison is doing a compare-by-value.  FIXME:  try
		 * coercing strings/unicodes to ilwdchars for comparison */
#if PY_MAJOR_VERSION < 3
		if(PyString_Check(other) || PyUnicode_Check(other)) {
#else
		if(PyUnicode_Check(other)) {
#endif
			PyErr_SetObject(PyExc_TypeError, other);
			return NULL;
		}
		Py_INCREF(Py_NotImplemented);
		return Py_NotImplemented;

	case 1:
		/* type agreement */
		break;

	default:
		/* ?? */
		PyErr_BadArgument();
		return NULL;
	}

	tbl_s = PyObject_GetAttr(self, table_name);
	col_s = PyObject_GetAttr(self, column_name);
	tbl_o = PyObject_GetAttr(other, table_name);
	col_o = PyObject_GetAttr(other, column_name);

	if(!(tbl_s && col_s && tbl_o && col_o))
		result = NULL;
	else {
		long r;

		/* compare by table name */
		r = PyUnicode_Compare(tbl_s, tbl_o);
		if(!r)
			/* break ties by comparing by column name */
			r = PyUnicode_Compare(col_s, col_o);
		if(!r)
			/* break ties by comparing by row ID */
			r = ((ligolw_ilwdchar *) self)->i - ((ligolw_ilwdchar *) other)->i;

		switch(op) {
		case Py_LT:
			result = (r < 0) ? Py_True : Py_False;
			break;

		case Py_LE:
			result = (r <= 0) ? Py_True : Py_False;
			break;

		case Py_EQ:
			result = (r == 0) ? Py_True : Py_False;
			break;

		case Py_NE:
			result = (r != 0) ? Py_True : Py_False;
			break;

		case Py_GT:
			result = (r > 0) ? Py_True : Py_False;
			break;

		case Py_GE:
			result = (r >= 0) ? Py_True : Py_False;
			break;

		default:
			PyErr_BadArgument();
			result = NULL;
		}
	}

	Py_XDECREF(tbl_s);
	Py_XDECREF(col_s);
	Py_XDECREF(tbl_o);
	Py_XDECREF(col_o);

	Py_XINCREF(result);

	return result;
}


static PyObject *ligolw_ilwdchar___str__(PyObject *self)
{
	ligolw_ilwdchar *ilwd = (ligolw_ilwdchar *) self;
	PyObject *tbl = PyObject_GetAttr(self, table_name);
	PyObject *col = PyObject_GetAttr(self, column_name);
	PyObject *result;

	if(!tbl || !col) {
		result = NULL;
	} else if(!PyUnicode_Check(tbl)) {
		PyErr_SetObject(PyExc_TypeError, tbl);
		result = NULL;
	} else if(!PyUnicode_Check(col)) {
		PyErr_SetObject(PyExc_TypeError, col);
		result = NULL;
	} else {
		/* 23 = 20 characters for a long int (2^63 == 19 digits,
		 * plus a possible "-" sign) + 2 ":" characters + a null
		 * terminator */
		Py_ssize_t tbl_len = PyUnicode_GetSize(tbl);
		Py_ssize_t col_len = PyUnicode_GetSize(col);
		wchar_t buff[tbl_len + col_len + 23];
		PyUnicode_AsWideChar((PyUnicodeObject *) tbl, buff, tbl_len);
		buff[tbl_len] = L':';
		PyUnicode_AsWideChar((PyUnicodeObject *) col, buff + tbl_len + 1, col_len);
		buff[tbl_len + 1 + col_len] = L':';
		swprintf(&buff[tbl_len + 1 + col_len], tbl_len + col_len + 23, L":%ld", ilwd->i);
		result = PyUnicode_FromWideChar(buff, wcslen(buff));
	}

	Py_XDECREF(tbl);
	Py_XDECREF(col);

	return result;
}


static PyObject *ligolw_ilwdchar___sub__(PyObject *self, PyObject *other)
{
	long delta = PyInt_AsLong(other);
	PyObject *new;

	if(PyErr_Occurred()) {
		/* can't be converted to int, maybe it's an ilwd:char of
		 * the same type as us */
		if(other->ob_type != self->ob_type)
			/* nope --> type error (already set from the int
			 * conversion failure) */
			return NULL;

		/* yes it is, return the ID difference as an int */
		PyErr_Clear();
		return PyInt_FromLong(((ligolw_ilwdchar *) self)->i - ((ligolw_ilwdchar *) other)->i);
	}

	if(!delta) {
		/* no change in value */
		Py_INCREF(self);
		return self;
	}

	new = PyType_GenericNew(self->ob_type, NULL, NULL);

	if(new)
		((ligolw_ilwdchar *) new)->i = ((ligolw_ilwdchar *) self)->i - delta;

	return new;
}


/*
 * The presence of this method allows ilwdchar sub-classes to be inserted
 * directly into SQLite databases as strings. See
 * http://www.python.org/dev/peps/pep-0246 for more information.
 *
 * FIXME: GvR has rejected that PEP, so this mechanism is obsolete.  Be
 * prepared to fix this, replacing it with whatever replaces it.
 *
 * FIXME: The return should be inside an "if protocol is
 * sqlite3.PrepareProtocol:" conditional, but that would require importing
 * sqlite3 which would break this module on FC4 boxes, and I'm not going to
 * spend time fixing something that's obsolete anyway.
 */


static PyObject *ligolw_ilwdchar___conform__(PyObject *self, PyObject *protocol)
{
	return PyObject_Unicode(self);
}


/*
 * Type
 */


static struct PyMethodDef methods[] = {
	{"__conform__", ligolw_ilwdchar___conform__, METH_O, "See http://www.python.org/dev/peps/pep-0246 for more information."},
	{NULL,}
};


PyTypeObject ligolw_ilwdchar_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(ligolw_ilwdchar),
	.tp_name = MODULE_NAME ".ilwdchar",
	.tp_doc =
"RAM-efficient row ID parent class.  This is only useful when subclassed in\n" \
"order to provide specific values of the class attributes \"table_name\"\n" \
"and \"column_name\".  The attributes must be unicode objects.\n" \
"\n" \
"Example:\n" \
"\n" \
">>> class ID(ilwdchar):\n" \
"...     __slots__ = ()\n" \
"...     table_name = u\"table_a\"\n" \
"...     column_name = u\"column_b\"\n" \
"... \n" \
">>> x = ID(10)\n" \
">>> print x\n" \
"table_a:column_b:10\n" \
">>> x = ID(u\" 10 \")	# ignores whitespace\n" \
">>> print x\n" \
"table_a:column_b:10\n" \
">>> print x + 35\n" \
"table_a:column_b:45\n" \
">>> y = ID(u\" table_a:column_b:10 \")	# ignores whitespace\n" \
">>> print x - y\n" \
"table_a:column_b:0\n" \
">>> x == y\n" \
"True\n" \
">>> x is y\n" \
"False\n" \
">>> len(set([x, y]))\n" \
"1\n" \
"\n" \
"Note that the two instances have the same hash value and compare as equal,\n" \
"and so only one of them remains in the set although they are not the same\n" \
"object.",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_hash = ligolw_ilwdchar___hash__,
	.tp_richcompare = ligolw_ilwdchar___richcompare__,
	.tp_str = ligolw_ilwdchar___str__,
	.tp_as_number = &(PyNumberMethods) {
		.nb_add = ligolw_ilwdchar___add__,
		.nb_int = ligolw_ilwdchar___int__,
		.nb_subtract = ligolw_ilwdchar___sub__,
	},
	.tp_methods = methods,
	.tp_new = ligolw_ilwdchar___new__,
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


void init_ilwd(void)
{
	/*
	 * Create the module.
	 */

	PyObject *module = Py_InitModule3(MODULE_NAME, NULL,
"C extension module providing the ilwdchar parent class for row ID classes."
	);

	/*
	 * Add the ilwdchar class.
	 */

	if(PyType_Ready(&ligolw_ilwdchar_Type) < 0)
		return;
	Py_INCREF(&ligolw_ilwdchar_Type);
	PyModule_AddObject(module, "ilwdchar", (PyObject *) &ligolw_ilwdchar_Type);

	/*
	 * Pre-construct table_name and column_name strings for attribute
	 * look-up
	 */

	table_name = PyUnicode_FromString("table_name");
	column_name = PyUnicode_FromString("column_name");
}
