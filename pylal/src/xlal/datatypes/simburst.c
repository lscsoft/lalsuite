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
#include <simburst.h>


#define MODULE_NAME PYLAL_SIMBURST_MODULE_NAME


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


static PyObject *sim_burst_simulation_id_type = NULL;
static PyObject *time_slide_id_type = NULL;
static PyObject *process_id_type = NULL;


/*
 * Member access
 */


static struct PyMemberDef members[] = {
	{"ra", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.ra), 0, "ra"},
	{"dec", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.dec), 0, "dec"},
	{"psi", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.psi), 0, "psi"},
	{"time_geocent_gps", T_INT, offsetof(pylal_SimBurst, sim_burst.time_geocent_gps.gpsSeconds), 0, "time_geocent_gps"},
	{"time_geocent_gps_ns", T_INT, offsetof(pylal_SimBurst, sim_burst.time_geocent_gps.gpsNanoSeconds), 0, "time_geocent_gps_ns"},
	{"time_geocent_gmst", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.time_geocent_gmst), 0, "time_geocent_gmst"},
	{"duration", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.duration), 0, "duration"},
	{"frequency", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.frequency), 0, "frequency"},
	{"bandwidth", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.bandwidth), 0, "bandwidth"},
	{"q", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.q), 0, "q"},
	{"pol_ellipse_angle", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.pol_ellipse_angle), 0, "pol_ellipse_angle"},
	{"pol_ellipse_e", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.pol_ellipse_e), 0, "pol_ellipse_e"},
	{"amplitude", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.amplitude), 0, "amplitude"},
	{"hrss", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.hrss), 0, "hrss"},
	{"egw_over_rsquared", T_DOUBLE, offsetof(pylal_SimBurst, sim_burst.egw_over_rsquared), 0, "egw_over_rsquared"},
	{"waveform_number", T_INT, offsetof(pylal_SimBurst, sim_burst.waveform_number), 0, "waveform_number"},
	{NULL,}
};


static struct PyGetSetDef getset[] = {
	{"process_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "process_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SimBurst, sim_burst.process_id), &process_id_type}},
	{"waveform", pylal_inline_string_get, pylal_inline_string_set, "waveform", &(struct pylal_inline_string_description) {offsetof(pylal_SimBurst, sim_burst.waveform), LIGOMETA_WAVEFORM_MAX}},
	{"time_slide_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "time_slide_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SimBurst, sim_burst.time_slide_id), &time_slide_id_type}},
	{"simulation_id", pylal_ilwdchar_id_get, pylal_ilwdchar_id_set, "simulation_id", &(struct pylal_ilwdchar_id_description) {offsetof(pylal_SimBurst, sim_burst.simulation_id), &sim_burst_simulation_id_type}},
	{NULL,}
};


static Py_ssize_t getreadbuffer(PyObject *self, Py_ssize_t segment, void **ptrptr)
{
	if(segment) {
		PyErr_SetString(PyExc_SystemError, "bad segment");
		return -1;
	}
	*ptrptr = &((pylal_SimBurst*)self)->sim_burst;
	return sizeof(((pylal_SimBurst*)self)->sim_burst);
}


static Py_ssize_t getsegcount(PyObject *self, Py_ssize_t *lenp)
{
	if(lenp)
		*lenp = sizeof(((pylal_SimBurst*)self)->sim_burst);
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
	pylal_SimBurst *new = (pylal_SimBurst *) PyType_GenericNew(type, args, kwds);

	if(!new)
		return NULL;

	/* done */
	return (PyObject *) new;
}


static PyObject *from_buffer(PyObject *cls, PyObject *args)
{
	const SimBurst *data;
	Py_ssize_t length;
	unsigned i;
	PyObject *result;

	if(!PyArg_ParseTuple(args, "s#", (const char **) &data, &length))
		return NULL;

	if(length % sizeof(SimBurst)) {
		PyErr_SetString(PyExc_ValueError, "buffer size is not an integer multiple of SimBurst struct size");
		return NULL;
	}
	length /= sizeof(SimBurst);

	result = PyTuple_New(length);
	if(!result)
		return NULL;
	for(i = 0; i < length; i++) {
		PyObject *item = PyType_GenericNew((PyTypeObject *) cls, NULL, NULL);
		if(!item) {
			Py_DECREF(result);
			return NULL;
		}
		/* memcpy sim_burst row */
		((pylal_SimBurst*)item)->sim_burst = *data++;

		PyTuple_SET_ITEM(result, i, item);
	}

	return result;
}


static struct PyMethodDef methods[] = {
	{"from_buffer", from_buffer, METH_VARARGS | METH_CLASS, "Construct a tuple of SimBurst objects from a buffer object.  The buffer is interpreted as a C array of SimBurst structures."},
	{NULL, }
};


/*
 * Type
 */


PyTypeObject pylal_simburst_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_SimBurst),
	.tp_doc = "LAL's SimBurst structure",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	.tp_members = members,
	.tp_methods = methods,
	.tp_getset = getset,
	.tp_as_buffer = &as_buffer,
	.tp_name = MODULE_NAME ".SimBurst",
	.tp_new = __new__,
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


PyMODINIT_FUNC initsimburst(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL, "Wrapper for LAL's SimBurst type.");

	/* Cached ID types */
	sim_burst_simulation_id_type = pylal_get_ilwdchar_class("sim_burst", "simulation_id");
	time_slide_id_type = pylal_get_ilwdchar_class("time_slide", "time_slide_id");
	process_id_type = pylal_get_ilwdchar_class("process", "process_id");

	/* SimBurst */
	_pylal_SimBurst_Type = &pylal_simburst_type;
	if(PyType_Ready(&pylal_SimBurst_Type) < 0)
		return;
	Py_INCREF(&pylal_SimBurst_Type);
	PyModule_AddObject(module, "SimBurst", (PyObject *) &pylal_SimBurst_Type);
}
