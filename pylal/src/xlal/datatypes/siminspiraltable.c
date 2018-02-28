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
#include <siminspiraltable.h>


#define MODULE_NAME PYLAL_SIMINSPIRALTABLE_MODULE_NAME


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


static PyObject *sim_inspiral_simulation_id_type = NULL;
static PyObject *process_id_type = NULL;


/*
 * Member access
 */


static struct PyMemberDef members[] = {
	{"geocent_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.geocent_end_time.gpsSeconds), 0, "geocent_end_time"},
	{"geocent_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.geocent_end_time.gpsNanoSeconds), 0, "geocent_end_time_ns"},
	{"h_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.h_end_time.gpsSeconds), 0, "h_end_time"},
	{"h_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.h_end_time.gpsNanoSeconds), 0, "h_end_time_ns"},
	{"l_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.l_end_time.gpsSeconds), 0, "l_end_time"},
	{"l_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.l_end_time.gpsNanoSeconds), 0, "l_end_time_ns"},
	{"g_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.g_end_time.gpsSeconds), 0, "g_end_time"},
	{"g_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.g_end_time.gpsNanoSeconds), 0, "g_end_time_ns"},
	{"t_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.t_end_time.gpsSeconds), 0, "t_end_time"},
	{"t_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.t_end_time.gpsNanoSeconds), 0, "t_end_time_ns"},
	{"v_end_time", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.v_end_time.gpsSeconds), 0, "v_end_time"},
	{"v_end_time_ns", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.v_end_time.gpsNanoSeconds), 0, "v_end_time_ns"},
	{"end_time_gmst", T_DOUBLE, offsetof(pylal_SimInspiralTable, sim_inspiral.end_time_gmst), 0, "end_time_gmst"},
	{"mass1", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.mass1), 0, "mass1"},
	{"mass2", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.mass2), 0, "mass2"},
	{"mchirp", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.mchirp), 0, "mchirp"},
	{"eta", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eta), 0, "eta"},
	{"distance", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.distance), 0, "distance"},
	{"longitude", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.longitude), 0, "longitude"},
	{"latitude", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.latitude), 0, "latitude"},
	{"inclination", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.inclination), 0, "inclination"},
	{"coa_phase", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.coa_phase), 0, "coa_phase"},
	{"polarization", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.polarization), 0, "polarization"},
	{"psi0", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.psi0), 0, "psi0"},
	{"psi3", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.psi3), 0, "psi3"},
	{"alpha", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha), 0, "alpha"},
	{"alpha1", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha1), 0, "alpha1"},
	{"alpha2", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha2), 0, "alpha2"},
	{"alpha3", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha3), 0, "alpha3"},
	{"alpha4", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha4), 0, "alpha4"},
	{"alpha5", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha5), 0, "alpha5"},
	{"alpha6", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.alpha6), 0, "alpha6"},
	{"beta", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.beta), 0, "beta"},
	{"spin1x", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin1x), 0, "spin1x"},
	{"spin1y", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin1y), 0, "spin1y"},
	{"spin1z", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin1z), 0, "spin1z"},
	{"spin2x", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin2x), 0, "spin2x"},
	{"spin2y", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin2y), 0, "spin2y"},
	{"spin2z", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.spin2z), 0, "spin2z"},
	{"theta0", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.theta0), 0, "theta0"},
	{"phi0", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.phi0), 0, "phi0"},
	{"f_lower", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.f_lower), 0, "f_lower"},
	{"f_final", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.f_final), 0, "f_final"},
	{"eff_dist_h", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_h), 0, "eff_dist_h"},
	{"eff_dist_l", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_l), 0, "eff_dist_l"},
	{"eff_dist_g", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_t), 0, "eff_dist_t"},
	{"eff_dist_t", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_g), 0, "eff_dist_g"},
	{"eff_dist_v", T_FLOAT, offsetof(pylal_SimInspiralTable, sim_inspiral.eff_dist_v), 0, "eff_dist_v"},
	{"numrel_mode_min", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.numrel_mode_min), 0, "numrel_mode_min"},
	{"numrel_mode_max", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.numrel_mode_max), 0, "numrel_mode_max"},
	{"amp_order", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.amp_order), 0, "amp_order"},
	{"bandpass", T_INT, offsetof(pylal_SimInspiralTable, sim_inspiral.bandpass), 0, "bandpass"},
	{NULL,}
};


static struct PyGetSetDef getset[] = {
	{"waveform", pylal_inline_string_get, pylal_inline_string_set, "waveform", &(struct pylal_inline_string_description) {offsetof(pylal_SimInspiralTable, sim_inspiral.waveform), LIGOMETA_WAVEFORM_MAX}},
	{"source", pylal_inline_string_get, pylal_inline_string_set, "source", &(struct pylal_inline_string_description) {offsetof(pylal_SimInspiralTable, sim_inspiral.source), LIGOMETA_SOURCE_MAX}},
	{"numrel_data", pylal_inline_string_get, pylal_inline_string_set, "numrel_data", &(struct pylal_inline_string_description) {offsetof(pylal_SimInspiralTable, sim_inspiral.numrel_data), LIGOMETA_STRING_MAX}},
	{"taper", pylal_inline_string_get, pylal_inline_string_set, "taper", &(struct pylal_inline_string_description) {offsetof(pylal_SimInspiralTable, sim_inspiral.taper), LIGOMETA_INSPIRALTAPER_MAX}},
	{"process_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "process_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SimInspiralTable, process_id_i), &process_id_type}},
	{"simulation_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "simulation_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SimInspiralTable, event_id.id), &sim_inspiral_simulation_id_type}},
	{NULL,}
};


static Py_ssize_t getreadbuffer(PyObject *self, Py_ssize_t segment, void **ptrptr)
{
	if(segment) {
		PyErr_SetString(PyExc_SystemError, "bad segment");
		return -1;
	}
	*ptrptr = &((pylal_SimInspiralTable*)self)->sim_inspiral;
	return sizeof(((pylal_SimInspiralTable*)self)->sim_inspiral);
}


static Py_ssize_t getsegcount(PyObject *self, Py_ssize_t *lenp)
{
	if(lenp)
		*lenp = sizeof(((pylal_SimInspiralTable*)self)->sim_inspiral);
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
	pylal_SimInspiralTable *new = (pylal_SimInspiralTable *) PyType_GenericNew(type, args, kwds);

	if(!new)
		return NULL;

	/* link the event_id pointer in the sngl_inspiral table structure
	 * to the event_id structure */
	new->sim_inspiral.event_id = &new->event_id;

	new->process_id_i = 0;
	new->event_id.id = 0;

	/* done */
	return (PyObject *) new;
}


static PyObject *from_buffer(PyObject *cls, PyObject *args)
{
	const SimInspiralTable *data;
	Py_ssize_t length;
	unsigned i;
	PyObject *result;

	if(!PyArg_ParseTuple(args, "s#", (const char **) &data, &length))
		return NULL;

	if(length % sizeof(SimInspiralTable)) {
		PyErr_SetString(PyExc_ValueError, "buffer size is not an integer multiple of SimInspiralTable struct size");
		return NULL;
	}
	length /= sizeof(SimInspiralTable);

	result = PyTuple_New(length);
	if(!result)
		return NULL;
	for(i = 0; i < length; i++) {
		PyObject *item = PyType_GenericNew((PyTypeObject *) cls, NULL, NULL);
		if(!item) {
			Py_DECREF(result);
			return NULL;
		}
		/* memcpy sim_inspiral row */
		((pylal_SimInspiralTable*)item)->sim_inspiral = *data++;
		/* repoint event_id to event_id structure */
		((pylal_SimInspiralTable*)item)->sim_inspiral.event_id = &((pylal_SimInspiralTable*)item)->event_id;

		PyTuple_SET_ITEM(result, i, item);
	}

	return result;
}


static struct PyMethodDef methods[] = {
	{"from_buffer", from_buffer, METH_VARARGS | METH_CLASS, "Construct a tuple of SimInspiralTable objects from a buffer object.  The buffer is interpreted as a C array of SimInspiralTable structures."},
	{NULL,}
};


/*
 * Type
 */


PyTypeObject pylal_siminspiraltable_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_SimInspiralTable),
	.tp_doc = "LAL's SimInspiralTable structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = members,
	.tp_methods = methods,
	.tp_getset = getset,
	.tp_as_buffer = &as_buffer,
	.tp_name = MODULE_NAME ".SimInspiralTable",
	.tp_new = __new__,
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


PyMODINIT_FUNC initsiminspiraltable(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "Wrapper for LAL's SimInspiralTable type.");

	/* Cached ID types */
	sim_inspiral_simulation_id_type = pylal_get_ilwdchar_class("sim_inspiral", "simulation_id");
	process_id_type = pylal_get_ilwdchar_class("process", "process_id");

	/* SimInspiralTable */
	_pylal_SimInspiralTable_Type = &pylal_siminspiraltable_type;
	if(PyType_Ready(&pylal_SimInspiralTable_Type) < 0)
		return;
	Py_INCREF(&pylal_SimInspiralTable_Type);
	PyModule_AddObject(module, "SimInspiralTable", (PyObject *) &pylal_SimInspiralTable_Type);
}
