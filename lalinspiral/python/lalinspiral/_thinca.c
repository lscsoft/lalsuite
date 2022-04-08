/*
 * Copyright (C) 2018  Kipp C. Cannon
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
 * how does this work?  this is a "get_coincs" class for the gstlal
 * inspiral search to use with the snglcoinc machinery.  a get_coincs
 * object is a callable object initialized with a time-ordered sequence of
 * triggers from one detector.  subsequently, when called with a single
 * trigger (from some other detector) it returns the in-order sequence of
 * triggers from its initialization set that are coincident with the given
 * trigger.  this provides the double-coincidences upon which the rest of
 * the coincidence machinery is built.
 *
 * for the gstlal inspiral search, two triggers are coincident if they were
 * obtained from the same template and are within some time window of each
 * other.  to find coincidences quickly, the initial sequence of triggers
 * is organized into a data structure like:
 *
 * array of event lists:
 *	template ID 0: [event, event, event, ...]
 *	template ID 1: [event, event, event, ...]
 *	...
 *
 * the event lists are sorted by template ID, and the events in each list
 * are time-ordered (which is ensured by providing them to the
 * initialization code in time order).  the events coincident with a
 * trigger are found by retrieving the event list matching that trigger's
 * ID with a bisection search, then the start and stop of the range of
 * times coincident with the trigger are computed and the sequence of
 * events that fall within that time range is returned by performing a
 * bisection search to identify the first such event, and then a brute
 * force search to find the last one (measured to be faster than a second
 * bisection search).
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


/* Ignore warnings in Python API itself */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"
#include <Python.h>
#pragma GCC diagnostic pop

#include <math.h>
#include <stdlib.h>

#include <lal/LALAtomicDatatypes.h>


/*
 * ============================================================================
 *
 *                               Internal Code
 *
 * ============================================================================
 */


/* retrieve the template ID of an event.  in modern versions of the
 * pipeline any integer may be used as a template ID, including negative
 * values, so error status must be checked in the calling code using
 * PyErr_Occurred() instead of by the return value alone.  mirroring the
 * behaviour of PyLong_AsLong(), we return -1 on error so that calling
 * codes can test for that special value before calling PyErr_Occurred() to
 * avoid the cost of that additional test when not needed */


static long get_template_id(PyObject *event)
{
	PyObject *template_id = PyObject_GetAttrString(event, "template_id");
	long val;

	if(!template_id)
		return -1;

	val = PyLong_AsLong(template_id);
	Py_DECREF(template_id);
	/* PyLong_AsLong()'s error code is the same as ours, so we don't do
	 * error checking here.  the calling code is responsible */

	return val;
}


/* get the end time of an event as integer nanoseconds from the GPS epoch.
 * returns -1 on error, but use PyErr_Occured() to confirm that an error
 * has occured, and that -1 is not actually the return value. */


static INT8 get_end_ns(PyObject *event)
{
	INT8 val = -1;
	PyObject *end = PyObject_GetAttrString(event, "end");
	if(end) {
		PyObject *ns = PyObject_CallMethod(end, "ns", NULL);
		Py_DECREF(end);
		if(ns) {
			val = PyLong_AsLong(ns);
			Py_DECREF(ns);
		}
	}
	return val;
}


/* compute the range of times, [start, stop), for which triggers from other
 * detectors are coincident with this one. */


static int compute_start_stop(INT8 *start_ns, INT8 *stop_ns, PyObject *event, PyObject *offset, PyObject* coinc_window)
{
	INT8 end_ns = get_end_ns(event);
	double offset_ns = PyFloat_AsDouble(offset);
	double coinc_window_ns = PyFloat_AsDouble(coinc_window);

	if((end_ns == -1 || offset_ns == -1. || coinc_window_ns == -1.) && PyErr_Occurred())
		return -1;
	offset_ns *= 1e9;
	coinc_window_ns *= 1e9;

	*start_ns = end_ns + lround(offset_ns - coinc_window_ns);
	*stop_ns = end_ns + lround(offset_ns + coinc_window_ns);

	return 0;
}


/*
 * event table
 */


struct event {
	INT8 end_ns;
	PyObject *event;
};


struct event_sequence {
	long template_id;
	struct event *events;
	int length;
};

/* comparison operator used with bsearch() to find an event list
 * corresponding to a given template ID */

static int event_sequence_get_cmp(const void *key, const void *obj)
{
	return (long) key - ((const struct event_sequence *) obj)->template_id;
}

/* comparison operator used with qsort() to sort event lists by template
 * ID */

static int event_sequence_sort_cmp(const void *a, const void *b)
{
	return ((const struct event_sequence *) a)->template_id - ((const struct event_sequence *) b)->template_id;
}

/* bsearch() wrapper to retrieve an event list by template ID */

static struct event_sequence *event_sequence_get(struct event_sequence *event_sequences, int n, long template_id)
{
	return bsearch((void *) template_id, event_sequences, n, sizeof(*event_sequences), event_sequence_get_cmp);
}

/* insert a time-ordered event into the array of event list arrays */

static int event_sequence_insert(struct event_sequence **event_sequences, int *n, PyObject *event)
{
	long template_id = get_template_id(event);
	INT8 end_ns = get_end_ns(event);
	struct event_sequence *event_sequence;
	int need_sort = 0;

	if((template_id == -1 || end_ns == -1) && PyErr_Occurred())
		return -1;

	/* get the event list for this template ID or create a new empty
	 * one if there isn't one yet.  the new list is appeneded, and the
	 * array of event lists put into the correct order later.  */

	event_sequence = event_sequence_get(*event_sequences, *n, template_id);
	if(!event_sequence) {
		struct event_sequence *new = realloc(*event_sequences, (*n + 1) * sizeof(**event_sequences));
		if(!new) {
			PyErr_SetString(PyExc_MemoryError, "realloc() failed");
			return -1;
		}
		*event_sequences = new;
		event_sequence = &new[*n];
		(*n)++;
		event_sequence->template_id = template_id;
		event_sequence->events = NULL;
		event_sequence->length = 0;
		need_sort = 1;
	}

	/* make room for the new event.  allocate space for triggers 128 at
	 * a time */

	if(!(event_sequence->length & 0x7f)) {
		struct event *new = realloc(event_sequence->events, ((event_sequence->length | 0x7f) + 1) * sizeof(*event_sequence->events));
		if(!new) {
			PyErr_SetString(PyExc_MemoryError, "realloc() failed");
			return -1;
		}
		event_sequence->events = new;
	}

	/* record the event in the event list */

	Py_INCREF(event);
	event_sequence->events[event_sequence->length].end_ns = end_ns;
	event_sequence->events[event_sequence->length].event = event;
	event_sequence->length++;

	/* if a new event list was created, put the array of event lists in
	 * order by template ID */

	if(need_sort)
		qsort(*event_sequences, *n, sizeof(**event_sequences), event_sequence_sort_cmp);

	return 0;
}


static void event_sequence_free(struct event_sequence *event_sequences, int n)
{
	while(n--) {
		struct event_sequence *event_sequence = &event_sequences[n];
		while(event_sequence->length--)
			Py_DECREF(event_sequence->events[event_sequence->length].event);
		free(event_sequence->events);
		event_sequence->events = NULL;
	}

	free(event_sequences);
}


static Py_ssize_t bisect_left(struct event *items, Py_ssize_t hi, INT8 end_ns)
{
	Py_ssize_t lo = 0;

	while(lo < hi) {
		Py_ssize_t mid = (lo + hi) / 2;
		if(items[mid].end_ns < end_ns)
			lo = mid + 1;
		else
			/* items[mid].end_ns >= end_ns */
			hi = mid;
	}

	return lo;
}

static PyObject *event_sequence_extract(struct event_sequence *event_sequence, INT8 start_ns, INT8 stop_ns)
{
	Py_ssize_t lo = bisect_left(event_sequence->events, event_sequence->length, start_ns);
	Py_ssize_t hi;
	Py_ssize_t i;
	PyObject *result;

	if(lo < 0)
		return NULL;

	for(hi = lo; event_sequence->events[hi].end_ns < stop_ns && hi < event_sequence->length; hi++);

	result = PyTuple_New(hi - lo);
	if(!result)
		return NULL;
	for(i = 0; lo < hi; i++, lo++) {
		PyObject *item = event_sequence->events[lo].event;
		Py_INCREF(item);
		PyTuple_SET_ITEM(result, i, item);
	}

	return result;
}


/*
 * ============================================================================
 *
 *                              get_coincs Class
 *
 * ============================================================================
 */


/*
 * Structure
 */


struct get_coincs {
	PyObject_HEAD

	struct event_sequence *event_sequences;
	int n_sequences;
};


/*
 * __del__() method
 */


static void get_coincs__del__(PyObject *self)
{
	struct get_coincs *get_coincs = (struct get_coincs *) self;

	event_sequence_free(get_coincs->event_sequences, get_coincs->n_sequences);
	get_coincs->event_sequences = NULL;
	get_coincs->n_sequences = 0;

	self->ob_type->tp_free(self);
}


/*
 * __init__() method
 */


static int get_coincs__init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	struct get_coincs *get_coincs = (struct get_coincs *) self;
	PyObject *events;
	PyObject **item, **last;

	if(kwds) {
		PyErr_SetString(PyExc_ValueError, "unexpected keyword arguments");
		return -1;
	}
	if(!PyArg_ParseTuple(args, "O", &events))
		return -1;

	if(get_coincs->event_sequences)
		event_sequence_free(get_coincs->event_sequences, get_coincs->n_sequences);
	get_coincs->event_sequences = NULL;
	get_coincs->n_sequences = 0;

	item = PySequence_Fast_ITEMS(events);
	if(!item)
		return -1;
	last = item + PySequence_Length(events);
	for(; item < last; item++)
		if(event_sequence_insert(&get_coincs->event_sequences, &get_coincs->n_sequences, *item) < 0) {
			event_sequence_free(get_coincs->event_sequences, get_coincs->n_sequences);
			get_coincs->event_sequences = NULL;
			get_coincs->n_sequences = 0;
			return -1;
		}

	return 0;
}


/*
 * __call__() method
 */


static PyObject *get_coincs__call__(PyObject *self, PyObject *args, PyObject *kwds)
{
	struct get_coincs *get_coincs = (struct get_coincs *) self;
	PyObject *event_a;
	PyObject *offset_a;
	PyObject *coinc_window;
	long template_id;
	struct event_sequence *event_sequence;
	INT8 start_ns, stop_ns;

	if(kwds) {
		PyErr_SetString(PyExc_ValueError, "unexpected keyword arguments");
		return NULL;
	}
	if(!PyArg_ParseTuple(args, "OOO", &event_a, &offset_a, &coinc_window))
		return NULL;

	template_id = get_template_id(event_a);
	if(template_id == -1 && PyErr_Occurred())
		return NULL;
	event_sequence = event_sequence_get(get_coincs->event_sequences, get_coincs->n_sequences, template_id);
	if(!(event_sequence && event_sequence->length)) {
		/* no events */
		Py_INCREF(Py_None);
		return Py_None;
	}

	if(compute_start_stop(&start_ns, &stop_ns, event_a, offset_a, coinc_window) < 0)
		return NULL;

	/* return the events matching event_a's .template_id and having end
	 * times in the range t - coinc_window <= end < t + coinc_window
	 * where t = event_a.end + offset_a */

	return event_sequence_extract(event_sequence, start_ns, stop_ns);
}


/*
 * Type information
 */


static PyTypeObject get_coincs_Type = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = MODULE_NAME ".get_coincs",
	.tp_basicsize = sizeof(struct get_coincs),
	.tp_dealloc = get_coincs__del__,
	.tp_call = get_coincs__call__,
	.tp_doc = "",
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_init = get_coincs__init__,
	.tp_new = PyType_GenericNew,
};


/*
 * ============================================================================
 *
 *                                Entry Point
 *
 * ============================================================================
 */


PyMODINIT_FUNC PyInit__thinca(void);	/* silence warning */
PyMODINIT_FUNC PyInit__thinca(void)
{
	static struct PyModuleDef modef = {
		PyModuleDef_HEAD_INIT,
		.m_name = MODULE_NAME,
		.m_doc = "",
		.m_size = -1,
		.m_methods = NULL,
	};
	PyObject *module = PyModule_Create(&modef);

	if(PyType_Ready(&get_coincs_Type) < 0) {
		Py_DECREF(module);
		module = NULL;
		goto done;
	}
	Py_INCREF((PyObject *) &get_coincs_Type);
	if(PyModule_AddObject(module, "get_coincs", (PyObject *) &get_coincs_Type) < 0) {
		Py_DECREF(module);
		module = NULL;
		goto done;
	}

done:
	return module;
}
