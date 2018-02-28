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
 *                   Python Wrapper For LAL's LALUnit Type
 *
 * ============================================================================
 */


#include <math.h>
#include <Python.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lalunit.h>


#define MODULE_NAME PYLAL_LALUNIT_MODULE_NAME


/*
 * ============================================================================
 *
 *                                LALUnit Type
 *
 * ============================================================================
 */


/*
 * Cached instances
 */


PyObject *pylal_LALUnitDimensionless = NULL;
PyObject *pylal_LALUnitMeter = NULL;
PyObject *pylal_LALUnitKiloGram = NULL;
PyObject *pylal_LALUnitSecond = NULL;
PyObject *pylal_LALUnitAmpere = NULL;
PyObject *pylal_LALUnitKelvin = NULL;
PyObject *pylal_LALUnitStrain = NULL;
PyObject *pylal_LALUnitADCCount = NULL;


/*
 * Methods
 */


static PyObject *__new__(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
	const char *s = NULL;
	pylal_LALUnit *new;

	if(!PyArg_ParseTuple(args, "|s", &s))
		return NULL;
	new = (pylal_LALUnit *) PyType_GenericNew(type, args, kwds);
	if(!new)
		return NULL;
	if(!s)
		new->unit = lalDimensionlessUnit;
	else if(!XLALParseUnitString(&new->unit, s)) {
		PyErr_Format(PyExc_ValueError, "cannot parse \"%s\"", s);
		Py_DECREF(new);
		return NULL;
	}
	return (PyObject *) new;
}


static PyObject *__repr__(PyObject *self)
{
	static char s[50];
	XLALUnitAsString(s, sizeof(s), &((pylal_LALUnit *) self)->unit);
	return PyString_FromFormat("LALUnit(\"%s\")", s);
}


static PyObject *__str__(PyObject *self)
{
	static char s[50];
	XLALUnitAsString(s, sizeof(s), &((pylal_LALUnit *) self)->unit);
	return PyString_FromString(s);
}


static int __cmp__(PyObject *self, PyObject *other)
{
	return XLALUnitCompare(&((pylal_LALUnit *) self)->unit, &((pylal_LALUnit *) other)->unit);
}


static long __hash__(PyObject *self)
{
	long hash = ((pylal_LALUnit *) self)->unit.powerOfTen;
	int i;
	for(i = 0; i < LALNumUnits; i++)
		hash ^= (((pylal_LALUnit *) self)->unit.unitNumerator[i] + ((pylal_LALUnit *) self)->unit.unitDenominatorMinusOne[i]) << ((8 * sizeof(hash) / LALNumUnits) * i);
	/* -1 reserved for errors */
	return hash == -1 ? -2 : hash;
}


static PyObject *__mul__(PyObject *self, PyObject *other)
{
	LALUnit new;
	return pylal_LALUnit_new(0, *XLALUnitMultiply(&new, &((pylal_LALUnit *) self)->unit, &((pylal_LALUnit *) other)->unit));
}


static PyObject *__div__(PyObject *self, PyObject *other)
{
	LALUnit new;
	return pylal_LALUnit_new(0, *XLALUnitDivide(&new, &((pylal_LALUnit *) self)->unit, &((pylal_LALUnit *) other)->unit));
}


static PyObject *__pow__(PyObject *self, PyObject *other, PyObject *modulo)
{
	/* FIXME */
	PyErr_SetNone(PyExc_NotImplementedError);
	return NULL;
}


static PyObject *__float__(PyObject *self)
{
	if(!XLALUnitIsDimensionless(&((pylal_LALUnit *) self)->unit)) {
		PyErr_SetString(PyExc_ValueError, "not dimensionless");
		return NULL;
	}
	return PyFloat_FromDouble(pow(10, ((pylal_LALUnit *) self)->unit.powerOfTen));
}


static PyObject *__invert__(PyObject *self)
{
	LALUnit new;
	return pylal_LALUnit_new(0, *XLALUnitInvert(&new, &((pylal_LALUnit *) self)->unit));
}


static int __coerce__(PyObject **p1, PyObject **p2)
{
	if(!PyObject_TypeCheck(*p1, &pylal_LALUnit_Type)) {
		PyObject *o;
		double power_of_ten;
		if(!PyNumber_Check(*p1))
			return -1;
		o = PyNumber_Float(*p1);
		if(!o)
			return -1;
		power_of_ten = log10(PyFloat_AsDouble(o));
		Py_DECREF(o);
		if(floor(power_of_ten) != power_of_ten) {
			PyErr_SetObject(PyExc_ValueError, *p1);
			return -1;
		}
		*p1 = pylal_LALUnit_new(floor(power_of_ten), lalDimensionlessUnit);
	} else
		Py_INCREF(*p1);
	if(!PyObject_TypeCheck(*p2, &pylal_LALUnit_Type)) {
		PyObject *o;
		double power_of_ten;
		if(!PyNumber_Check(*p2))
			return -1;
		o = PyNumber_Float(*p2);
		if(!o)
			return -1;
		power_of_ten = log10(PyFloat_AsDouble(o));
		Py_DECREF(o);
		if(floor(power_of_ten) != power_of_ten) {
			PyErr_SetObject(PyExc_ValueError, *p2);
			return -1;
		}
		*p2 = pylal_LALUnit_new(floor(power_of_ten), lalDimensionlessUnit);
	} else
		Py_INCREF(*p2);
	return 0;
}


/*
 * Type
 */


static PyNumberMethods as_number = {
	.nb_multiply = __mul__,
	.nb_divide = __div__,
	.nb_power = __pow__,
	.nb_invert = __invert__,
	.nb_coerce = __coerce__,
	.nb_float = __float__
};


static PyTypeObject pylal_lalunit_type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(pylal_LALUnit),
	.tp_doc =
"LALUnit structure.  This is an immutable type (it can be hashed, put into\n" \
"sets, used as a dictionary key, etc.).\n" \
"\n" \
"Example:\n" \
"\n" \
">>> x = LALUnit(\"m\")\n" \
">>> y = LALUnit(\"s\")\n" \
">>> 1000 * x / y\n" \
"LALUnit(\"10^3 m s^-1\")\n",
	.tp_flags = Py_TPFLAGS_DEFAULT,
	.tp_name = MODULE_NAME ".LALUnit",
	.tp_new = __new__,
	.tp_compare = __cmp__,
	.tp_repr = __repr__,
	.tp_as_number = &as_number,
	.tp_hash = __hash__,
	.tp_str = __str__
};


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


PyMODINIT_FUNC initlalunit(void)
{
	PyObject *module = Py_InitModule3(MODULE_NAME, NULL,
"Wrapper for LAL's LALUnit type.  Several pre-defined unit constants are also provided.");

	/*
	 * LALUnit
	 */

	_pylal_LALUnit_Type = &pylal_lalunit_type;
	if(PyType_Ready(&pylal_LALUnit_Type) < 0)
		return;
	Py_INCREF(&pylal_LALUnit_Type);
	PyModule_AddObject(module, "LALUnit", (PyObject *) &pylal_LALUnit_Type);

	/*
	 * cached instances
	 */

	pylal_LALUnitDimensionless = pylal_LALUnit_new(0, lalDimensionlessUnit);
	PyModule_AddObject(module, "lalDimensionlessUnit", (PyObject *) pylal_LALUnitDimensionless);
	pylal_LALUnitMeter = pylal_LALUnit_new(0, lalMeterUnit);
	PyModule_AddObject(module, "lalMeterUnit", (PyObject *) pylal_LALUnitMeter);
	pylal_LALUnitKiloGram = pylal_LALUnit_new(0, lalKiloGramUnit);
	PyModule_AddObject(module, "lalKiloGramUnit", (PyObject *) pylal_LALUnitKiloGram);
	pylal_LALUnitSecond = pylal_LALUnit_new(0, lalSecondUnit);
	PyModule_AddObject(module, "lalSecondUnit", (PyObject *) pylal_LALUnitSecond);
	pylal_LALUnitAmpere = pylal_LALUnit_new(0, lalAmpereUnit);
	PyModule_AddObject(module, "lalAmpereUnit", (PyObject *) pylal_LALUnitAmpere);
	pylal_LALUnitKelvin = pylal_LALUnit_new(0, lalKelvinUnit);
	PyModule_AddObject(module, "lalKelvinUnit", (PyObject *) pylal_LALUnitKelvin);
	pylal_LALUnitStrain = pylal_LALUnit_new(0, lalStrainUnit);
	PyModule_AddObject(module, "lalStrainUnit", (PyObject *) pylal_LALUnitStrain);
	pylal_LALUnitADCCount = pylal_LALUnit_new(0, lalADCCountUnit);
	PyModule_AddObject(module, "lalADCCountUnit", (PyObject *) pylal_LALUnitADCCount);
}
