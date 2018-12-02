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
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <Python.h>
#include <stdlib.h>


/*
 * ============================================================================
 *
 *                               Internal Code
 *
 * ============================================================================
 */


static long get_template_id(PyObject *event)
{
	PyObject *template_id = PyObject_GetAttrString(event, "template_id");
	long val;

	if(!template_id)
		return -1;

#if PY_MAJOR_VERSION < 3
	val = PyInt_AsLong(template_id);
#else
	val = PyLong_AsLong(template_id);
#endif
	Py_DECREF(template_id);
	if(val < 0 && PyErr_Occurred())
		return -1;

	return val;
}


static int compute_start_stop(PyObject **start, PyObject **stop, PyObject *event, PyObject *offset, PyObject* coinc_window)
{
	PyObject *end = PyObject_GetAttrString(event, "end");
	if(!end)
		return -1;

	{
	PyObject *tmp = PyNumber_InPlaceAdd(end, offset);
	Py_DECREF(end);
	if(!tmp)
		return -1;
	end = tmp;
	}

	*start = PyNumber_Subtract(end, coinc_window);
	*stop = PyNumber_Add(end, coinc_window);
	Py_DECREF(end);
	if(!*start || !*stop) {
		Py_XDECREF(*start);
		Py_XDECREF(*stop);
		return -1;
	}

	return 0;
}


static Py_ssize_t bisect_left(PyObject **items, Py_ssize_t hi, PyObject *val)
{
	Py_ssize_t lo = 0;

	while(lo < hi) {
		Py_ssize_t mid = (lo + hi) / 2;
		int result = PyObject_RichCompareBool(items[mid], val, Py_LT);
		if(result > 0)
			/* items[mid] < val */
			lo = mid + 1;
		else if(result == 0)
			/* items[mid] >= val */
			hi = mid;
		else
			/* error */
			return -1;
	}

	return lo;
}


/*
 * event table
 */


struct event_sequence {
	long template_id;
	PyObject **events;
	int length;
};


static int event_sequence_get_cmp(const void *key, const void *obj)
{
	return (long) key - ((const struct event_sequence *) obj)->template_id;
}


static int event_sequence_sort_cmp(const void *a, const void *b)
{
	return ((const struct event_sequence *) a)->template_id - ((const struct event_sequence *) b)->template_id;
}


static struct event_sequence *event_sequence_get(struct event_sequence *event_sequences, int n, long template_id)
{
	return bsearch((void *) template_id, event_sequences, n, sizeof(*event_sequences), event_sequence_get_cmp);
}


/* this consumes a reference for event */
static int event_sequence_insert(struct event_sequence **event_sequences, int *n, PyObject *event)
{
	long template_id = get_template_id(event);
	struct event_sequence *event_sequence;
	PyObject **new_events;
	int need_sort = 0;

	if(template_id < 0)
		return -1;

	event_sequence = event_sequence_get(*event_sequences, *n, template_id);
	if(!event_sequence) {
		struct event_sequence *new = realloc(*event_sequences, (*n + 1) * sizeof(**event_sequences));
		if(!new) {
			Py_DECREF(event);
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

	new_events = realloc(event_sequence->events, (event_sequence->length + 1) * sizeof(*event_sequence->events));
	if(!new_events) {
		Py_DECREF(event);
		PyErr_SetString(PyExc_MemoryError, "realloc() failed");
		return -1;
	}
	event_sequence->events = new_events;
	event_sequence->events[event_sequence->length] = event;
	event_sequence->length++;

	if(need_sort)
		qsort(*event_sequences, *n, sizeof(**event_sequences), event_sequence_sort_cmp);

	return 0;
}


static void event_sequence_free(struct event_sequence *event_sequences, int n)
{
	while(n--) {
		struct event_sequence *event_sequence = &event_sequences[n];
		while(event_sequence->length--)
			Py_DECREF(event_sequence->events[event_sequence->length]);
		free(event_sequence->events);
		event_sequence->events = NULL;
	}

	free(event_sequences);
}


static PyObject *event_sequence_extract(struct event_sequence *event_sequence, PyObject *start, PyObject *stop)
{
	Py_ssize_t lo = bisect_left(event_sequence->events, event_sequence->length, start);
	Py_ssize_t hi;
	Py_ssize_t i;
	PyObject *result;

	Py_DECREF(start);

	if(lo < 0) {
		Py_DECREF(stop);
		return NULL;
	}

	for(hi = lo; hi < event_sequence->length; hi++) {
		int LT = PyObject_RichCompareBool(event_sequence->events[hi], stop, Py_LT);
		if(LT > 0)
			/* events[hi] < stop */
			continue;
		else if(LT == 0)
			/* events[hi] >= stop */
			break;
		/* error */
		Py_DECREF(stop);
		return NULL;
	}
	Py_DECREF(stop);

	result = PyTuple_New(hi - lo);
	if(!result)
		return NULL;
	for(i = 0; lo < hi; i++, lo++) {
		PyObject *item = event_sequence->events[lo];
		Py_INCREF(item);
		PyTuple_SET_ITEM(result, i, item);
	}

	return result;
}


/*
 * ============================================================================
 *
 *                    Rate Estimation --- Posterior Class
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
	for(; item < last; item++) {
		long template_id = get_template_id(*item);
		if(template_id < 0) {
			event_sequence_free(get_coincs->event_sequences, get_coincs->n_sequences);
			return -1;
		}
		Py_INCREF(*item);
		if(event_sequence_insert(&get_coincs->event_sequences, &get_coincs->n_sequences, *item) < 0) {
			event_sequence_free(get_coincs->event_sequences, get_coincs->n_sequences);
			return -1;
		}
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
	PyObject *start, *stop;

	if(kwds) {
		PyErr_SetString(PyExc_ValueError, "unexpected keyword arguments");
		return NULL;
	}
	if(!PyArg_ParseTuple(args, "OOO", &event_a, &offset_a, &coinc_window))
		return NULL;

	template_id = get_template_id(event_a);
	if(template_id < 0)
		return NULL;
	event_sequence = event_sequence_get(get_coincs->event_sequences, get_coincs->n_sequences, template_id);
	if(!(event_sequence && event_sequence->length)) {
		/* no events */
		Py_INCREF(Py_None);
		return Py_None;
	}

	if(compute_start_stop(&start, &stop, event_a, offset_a, coinc_window) < 0)
		return NULL;

	/* return the events matching event_a's .template_id and having end
	 * times in the range t - coinc_window <= end < t + coinc_window
	 * where t = event_a.end + offset_a */
	return event_sequence_extract(event_sequence, start, stop);
}


/*
 * Type information
 */


static PyTypeObject get_coincs_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(struct get_coincs),
	.tp_call = get_coincs__call__,
	.tp_dealloc = get_coincs__del__,
	.tp_doc = "",
#if PY_MAJOR_VERSION < 3
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES,
#else
	.tp_flags = Py_TPFLAGS_DEFAULT,
#endif
	.tp_init = get_coincs__init__,
	.tp_name = MODULE_NAME ".get_coincs",
	.tp_new = PyType_GenericNew,
};


/*
 * ============================================================================
 *
 *                                Entry Point
 *
 * ============================================================================
 */


void init_thinca(void);	/* silence warning */
void init_thinca(void)
{
#if PY_MAJOR_VERSION < 3
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "");
#else
	static struct PyModuleDef modef = {
		PyModuleDef_HEAD_INIT,
		MODULE_NAME,
		"",
		-1,
		NULL
	};
	PyObject *module = PyModule_Create(&modef);
#endif

	if(PyType_Ready(&get_coincs_Type) < 0)
		return;
	Py_INCREF((PyObject *) &get_coincs_Type);
	PyModule_AddObject(module, "get_coincs", (PyObject *) &get_coincs_Type);
}
