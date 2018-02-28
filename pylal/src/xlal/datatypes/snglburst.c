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


#include <Python.h>
#include <structmember.h>
#include <lal/LIGOMetadataTables.h>
#include <misc.h>
#include <snglburst.h>


#define MODULE_NAME PYLAL_SNGLBURST_MODULE_NAME


/*
 * ============================================================================
 *
 *                                    Type
 *
 * ============================================================================
 */


/*
 * Cached ID types
 */


static PyObject *sngl_burst_event_id_type = NULL;
static PyObject *process_id_type = NULL;


/*
 * Member access
 */


static struct PyMemberDef members[] = {
	{"start_time", T_INT, offsetof(pylal_SnglBurst, sngl_burst.start_time.gpsSeconds), 0, "start_time"},
	{"start_time_ns", T_INT, offsetof(pylal_SnglBurst, sngl_burst.start_time.gpsNanoSeconds), 0, "start_time_ns"},
	{"peak_time", T_INT, offsetof(pylal_SnglBurst, sngl_burst.peak_time.gpsSeconds), 0, "peak_time"},
	{"peak_time_ns", T_INT, offsetof(pylal_SnglBurst, sngl_burst.peak_time.gpsNanoSeconds), 0, "peak_time_ns"},
	{"duration", T_FLOAT, offsetof(pylal_SnglBurst, sngl_burst.duration), 0, "duration"},
	{"central_freq", T_FLOAT, offsetof(pylal_SnglBurst, sngl_burst.central_freq), 0, "central_freq"},
	{"bandwidth", T_FLOAT, offsetof(pylal_SnglBurst, sngl_burst.bandwidth), 0, "bandwidth"},
	{"amplitude", T_FLOAT, offsetof(pylal_SnglBurst, sngl_burst.amplitude), 0, "amplitude"},
	{"snr", T_FLOAT, offsetof(pylal_SnglBurst, sngl_burst.snr), 0, "snr"},
	{"confidence", T_FLOAT, offsetof(pylal_SnglBurst, sngl_burst.confidence), 0, "confidence"},
	{"chisq", T_DOUBLE, offsetof(pylal_SnglBurst, sngl_burst.chisq), 0, "chisq"},
	{"chisq_dof", T_DOUBLE, offsetof(pylal_SnglBurst, sngl_burst.chisq_dof), 0, "chisq_dof"},
	{NULL,}
};


static struct PyGetSetDef getset[] = {
	{"ifo", pylal_inline_string_get, pylal_inline_string_set, "ifo", &(struct pylal_inline_string_description) {offsetof(pylal_SnglBurst, sngl_burst.ifo), LIGOMETA_IFO_MAX}},
	{"search", pylal_inline_string_get, pylal_inline_string_set, "search", &(struct pylal_inline_string_description) {offsetof(pylal_SnglBurst, sngl_burst.search), LIGOMETA_SEARCH_MAX}},
	{"channel", pylal_inline_string_get, pylal_inline_string_set, "channel", &(struct pylal_inline_string_description) {offsetof(pylal_SnglBurst, sngl_burst.channel), LIGOMETA_CHANNEL_MAX}},
	{"process_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "process_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SnglBurst, sngl_burst.process_id), &process_id_type}},
	{"event_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "event_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SnglBurst, sngl_burst.event_id), &sngl_burst_event_id_type}},
	{NULL,}
};


static Py_ssize_t getreadbuffer(PyObject *self, Py_ssize_t segment, void **ptrptr)
{
	if(segment) {
		PyErr_SetString(PyExc_SystemError, "bad segment");
		return -1;
	}
	*ptrptr = &((pylal_SnglBurst*)self)->sngl_burst;
	return sizeof(((pylal_SnglBurst*)self)->sngl_burst);
}


static Py_ssize_t getsegcount(PyObject *self, Py_ssize_t *lenp)
{
	if(lenp)
		*lenp = sizeof(((pylal_SnglBurst*)self)->sngl_burst);
	return 1;
}


static PyBufferProcs as_buffer = {
	.bf_getreadbuffer = getreadbuffer,
	.bf_getsegcount = getsegcount,
	.bf_getwritebuffer = NULL,
	.bf_getcharbuffer = NULL
};


/*
 * Methods
 */


static PyObject *__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	/* call the generic __new__() */
	pylal_SnglBurst *new = (pylal_SnglBurst *) PyType_GenericNew(type, args, kwds);

	if(!new)
		return NULL;

	/* done */
	return (PyObject *) new;
}


static PyObject *from_buffer(PyObject *cls, PyObject *args)
{
	const SnglBurst *data;
	Py_ssize_t length;
	unsigned i;
	PyObject *result;

	if(!PyArg_ParseTuple(args, "s#", (const char **) &data, &length))
		return NULL;

	if(length % sizeof(SnglBurst)) {
		PyErr_SetString(PyExc_ValueError, "buffer size is not an integer multiple of SnglBurst struct size");
		return NULL;
	}
	length /= sizeof(SnglBurst);

	result = PyTuple_New(length);
	if(!result)
		return NULL;
	for(i = 0; i < length; i++) {
		PyObject *item = PyType_GenericNew((PyTypeObject *) cls, NULL, NULL);
		if(!item) {
			Py_DECREF(result);
			return NULL;
		}
		/* memcpy sngl_burst row */
		((pylal_SnglBurst*)item)->sngl_burst = *data++;

		PyTuple_SET_ITEM(result, i, item);
	}

	return result;
}


static struct PyMethodDef methods[] = {
	{"from_buffer", from_buffer, METH_VARARGS | METH_CLASS, "Construct a tuple of SnglBurst objects from a buffer object.  The buffer is interpreted as a C array of SnglBurst structures."},
	{NULL, }
};


/*
 * Type
 */


PyTypeObject pylal_snglburst_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_SnglBurst),
	.tp_doc = "LAL's SnglBurst structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = members,
	.tp_methods = methods,
	.tp_getset = getset,
	.tp_as_buffer = &as_buffer,
	.tp_name = MODULE_NAME ".SnglBurst",
	.tp_new = __new__,
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


PyMODINIT_FUNC initsnglburst(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "Wrapper for LAL's SnglBurst type.");

	/* Cached ID types */
	sngl_burst_event_id_type = pylal_get_ilwdchar_class("sngl_burst", "event_id");
	process_id_type = pylal_get_ilwdchar_class("process", "process_id");

	/* SnglBurst */
	_pylal_SnglBurst_Type = &pylal_snglburst_type;
	if(PyType_Ready(&pylal_SnglBurst_Type) < 0)
		return;
	Py_INCREF(&pylal_SnglBurst_Type);
	PyModule_AddObject(module, "SnglBurst", (PyObject *) &pylal_SnglBurst_Type);
}
