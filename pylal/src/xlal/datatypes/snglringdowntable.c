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


#include <Python.h>
#include <structmember.h>
#include <lal/LIGOMetadataTables.h>
#include <misc.h>
#include <snglringdowntable.h>


#define MODULE_NAME PYLAL_SNGLRINGDOWNTABLE_MODULE_NAME


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


static PyObject *sngl_ringdown_event_id_type = NULL;
static PyObject *process_id_type = NULL;


/*
 * Member access
 */


static struct PyMemberDef members[] = {
	{"start_time", T_INT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.start_time.gpsSeconds), 0, "start_time"},
	{"start_time_ns", T_INT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.start_time.gpsNanoSeconds), 0, "start_time_ns"},
	{"start_time_gmst", T_DOUBLE, offsetof(pylal_SnglRingdownTable, sngl_ringdown.start_time_gmst), 0, "start_time_gmst"},
	{"frequency", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.frequency), 0, "frequency"},
	{"quality", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.quality), 0, "quality"},
	{"phase", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.phase), 0, "phase"},
	{"mass", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.mass), 0, "mass"},
	{"spin", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.spin), 0, "spin"},
	{"epsilon", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.epsilon), 0, "epsilon"},
	{"num_clust_trigs", T_INT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.num_clust_trigs), 0, "num_clust_trigs"},
	{"ds2_H1H2", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.ds2_H1H2), 0, "ds2_H1H2"},
	{"ds2_H1L1", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.ds2_H1L1), 0, "ds2_H1L1"},
	{"ds2_H1V1", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.ds2_H1L1), 0, "ds2_H1V1"},
	{"ds2_H2L1", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.ds2_H2L1), 0, "ds2_H2L1"},
	{"ds2_H2V1", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.ds2_H2L1), 0, "ds2_H2V1"},
	{"ds2_L1V1", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.ds2_H2L1), 0, "ds2_L1V1"},
	{"amplitude", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.amplitude), 0, "amplitude"},
	{"snr", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.snr), 0, "snr"},
	{"eff_dist", T_FLOAT, offsetof(pylal_SnglRingdownTable, sngl_ringdown.eff_dist), 0, "eff_dist"},
	{"sigma_sq", T_DOUBLE, offsetof(pylal_SnglRingdownTable, sngl_ringdown.sigma_sq), 0, "sigma_sq"},
	{NULL,}
};


static struct PyGetSetDef getset[] = {
	{"ifo", pylal_inline_string_get, pylal_inline_string_set, "ifo", &(struct pylal_inline_string_description) {offsetof(pylal_SnglRingdownTable, sngl_ringdown.ifo), LIGOMETA_IFO_MAX}},
	{"channel", pylal_inline_string_get, pylal_inline_string_set, "channel", &(struct pylal_inline_string_description) {offsetof(pylal_SnglRingdownTable, sngl_ringdown.channel), LIGOMETA_CHANNEL_MAX}},
	{"process_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "process_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SnglRingdownTable, process_id_i), &process_id_type}},
	{"event_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "event_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SnglRingdownTable, event_id.id), &sngl_ringdown_event_id_type}},
	{NULL,}
};


/*
 * Methods
 */


static PyObject *__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	/* call the generic __new__() */
	pylal_SnglRingdownTable *new = (pylal_SnglRingdownTable *) PyType_GenericNew(type, args, kwds);

	if(!new)
		return NULL;

	/* link the event_id pointer in the sngl_ringdown table structure
	 * to the event_id structure */
	new->sngl_ringdown.event_id = &new->event_id;

	new->process_id_i = 0;
	new->event_id.id = 0;

	/* done */
	return (PyObject *) new;
}


/*
 * Type
 */


PyTypeObject pylal_snglringdowntable_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_SnglRingdownTable),
	.tp_doc = "LAL's SnglRingdownTable structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = members,
	.tp_getset = getset,
	.tp_name = MODULE_NAME ".SnglRingdownTable",
	.tp_new = __new__,
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef functions[] = {
	{NULL,}
};


PyMODINIT_FUNC initsnglringdowntable(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, functions, "Wrapper for LAL's SnglRingdownTable type.");

	/* Cached ID types */
	process_id_type = pylal_get_ilwdchar_class("process", "process_id");
	sngl_ringdown_event_id_type = pylal_get_ilwdchar_class("sngl_ringdown", "event_id");

	/* SnglRingdownTable */
	_pylal_SnglRingdownTable_Type = &pylal_snglringdowntable_type;
	if(PyType_Ready(&pylal_SnglRingdownTable_Type) < 0)
		return;
	Py_INCREF(&pylal_SnglRingdownTable_Type);
	PyModule_AddObject(module, "SnglRingdownTable", (PyObject *) &pylal_SnglRingdownTable_Type);
}
