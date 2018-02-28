/*
 * Copyright (C) 2006  Kipp Cannon
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
 *                   Python Wrapper For LAL's Date Package
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <stdlib.h>
#include <time.h>
#include <numpy/arrayobject.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <misc.h>
#include <datatypes/ligotimegps.h>


#define MODULE_NAME "pylal.xlal.date"


/*
 * ============================================================================
 *
 *                              Time Conversion
 *
 * ============================================================================
 */


/*
 * Convert a struct tm of the kind Python uses to/from a struct tm of the
 * kind the C library uses.
 *
 * Python:
 *	1900 = year 1900
 *	January = month 1
 *	Monday = week day 0
 *	January 1st = year day 1
 *
 * C:
 *	1900 = year 0
 *	January = month 0
 *	Monday = week day 1
 *	January 1st = year day 0
 */


static struct tm *struct_tm_python_to_c(struct tm *tm)
{
	tm->tm_year -= 1900;
	tm->tm_mon -= 1;
	tm->tm_wday = (tm->tm_wday + 8) % 7;
	tm->tm_yday -= 1;

	return tm;
}


static struct tm *struct_tm_c_to_python(struct tm *tm)
{
	tm->tm_year += 1900;
	tm->tm_mon += 1;
	tm->tm_wday = (tm->tm_wday + 6) % 7;
	tm->tm_yday += 1;

	return tm;
}


/*
 * GPS to/from UTC
 */

PyDoc_STRVAR(pylal_XLALGPSToUTC__doc__,
"Returns a time structure for a given GPS time (for integer second only!)\n"
"Example:\n"
"\n"
">>> from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS\n"
">>> gps = LIGOTimeGPS(969953934)\n"
">>> XLALGPSToUTC(gps)\n"
"(2010, 10, 1, 7, 38, 39, 4, 274, 0)\n"
"\n"
"The inverse operation is:\n"
"\n"
">>> import time\n"
">>> tm = time.struct_time((2010, 10, 1, 7, 38, 39, 4, 274, 0))\n"
">>> XLALUTCToGPS(tm)\n"
"LIGOTimeGPS(969953934,0)\n");

static PyObject *pylal_XLALGPSToUTC(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *gps;
	struct tm utc;

	/* LIGOTimeGPS */
	if(!PyArg_ParseTuple(args, "O!:XLALGPSToUTC", &pylal_LIGOTimeGPS_Type, &gps))
		return NULL;
	if(gps->gps.gpsNanoSeconds) {
		PyErr_SetString(PyExc_TypeError, "cannot convert non-integer seconds");
		return NULL;
	}

	XLALGPSToUTC(&utc, gps->gps.gpsSeconds);
	struct_tm_c_to_python(&utc);

	/* time.struct_time */
	return Py_BuildValue("(iiiiiiiii)", utc.tm_year, utc.tm_mon, utc.tm_mday, utc.tm_hour, utc.tm_min, utc.tm_sec, utc.tm_wday, utc.tm_yday, utc.tm_isdst);
}

PyDoc_STRVAR(pylal_XLALUTCToGPS__doc__,
"Returns the GPS time for a specified UTC time structure.\n"
"Example:\n"
"\n"
">>> import time\n"
">>> tm = time.gmtime()\n"
">>> XLALUTCToGPS(tm)\n"
"LIGOTimeGPS(971226989,0)\n"
"\n"
"The inverse operation is:\n"
"\n"
">>> from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS\n"
">>> gps = LIGOTimeGPS(971226989)\n"
">>> XLALGPSToUTC(gps)\n"
"(2010, 10, 16, 1, 16, 14, 5, 289, 0)\n"
);

static PyObject *pylal_XLALUTCToGPS(PyObject *self, PyObject *args)
{
	LIGOTimeGPS gps;
	struct tm utc;

	/* time.struct_time */
	if(!PyArg_ParseTuple(args, "(iiiiiiiii):XLALUTCToGPS", &utc.tm_year, &utc.tm_mon, &utc.tm_mday, &utc.tm_hour, &utc.tm_min, &utc.tm_sec, &utc.tm_wday, &utc.tm_yday, &utc.tm_isdst))
		return NULL;

	XLALGPSSet(&gps, XLALUTCToGPS(struct_tm_python_to_c(&utc)), 0);

	/* LIGOTimeGPS */
	return pylal_LIGOTimeGPS_new(gps);
}


/*
 * Sidereal time
 */

PyDoc_STRVAR(pylal_XLALGreenwichSiderealTime__doc__,
"Returns the Greenwich Sidereal Time IN RADIANS corresponding to a specified\n"
"GPS time.\n"
"\n"
"Apparent sidereal time is computed by providing the equation of equinoxes in\n"
"units of seconds. For mean sidereal time, set this parameter to 0.\n"
"\n"
"This function returns the sidereal time in radians measured from the Julian\n"
"epoch (current J2000). The result is NOT modulo 2 pi.\n"
"\n"
"Inspired by the function sidereal_time() in the NOVAS-C library, version\n"
"2.0.1, which is dated December 10th, 1999, and carries the following references:\n"
"\n"
"Aoki, et al. (1982) Astronomy and Astrophysics 105, 359-361. Kaplan, G. H.\n"
"\"NOVAS: Naval Observatory Vector Astrometry Subroutines\"; USNO internal\n"
"document dated 20 Oct 1988; revised 15 Mar 1990.\n"
"\n"
"See http://aa.usno.navy.mil/software/novas for more information.\n"
"\n"
"Example:\n"
"\n"
">>> from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS\n"
">>> gps = LIGOTimeGPS(969953934)\n"
">>> XLALGreenwichSiderealTime(gps, 0)\n"
"24739.075161218552\n");

static PyObject *pylal_XLALGreenwichSiderealTime(PyObject *self, PyObject *args)
{
	pylal_LIGOTimeGPS *gps;
	double equation_of_equinoxes;

	/* LIGOTimeGPS, float */
	if(!PyArg_ParseTuple(args, "O!d:XLALGreenwichSiderealTime", &pylal_LIGOTimeGPS_Type, &gps, &equation_of_equinoxes))
		return NULL;

	/* float */
	return PyFloat_FromDouble(XLALGreenwichSiderealTime(&gps->gps, equation_of_equinoxes));
}


/*
 * ============================================================================
 *
 *                             Propogation Delay
 *
 * ============================================================================
 */

PyDoc_STRVAR(pylal_XLALArrivalTimeDiff__doc__,
"Calculate the light-travel time difference (in units of nanoseconds) between\n"
"two locations, given a sky location and a time.\n"
"\n"
"Example: The arrival time difference from a signal in Andromeda on 1 Feb\n"
"2010 at 15:23:11. \n"
"\n"
">>> import time\n"
">>> tm = time.strptime(\"1 Feb 2007 15:23:11\",\"%d %b %Y %H:%M:%S\")\n"
">>> gps = XLALUTCToGPS(tm)\n"
">>> from pylal.xlal import tools\n"
">>> location_L1 = tools.cached_detector[\"LLO_4k\"].location\n"
">>> location_H1 = tools.cached_detector[\"LHO_4k\"].location\n"
">>> ra = 0.1864\n"
">>> de = 0.702\n"
">>> XLALArrivalTimeDiff(location_L1, location_H1, ra, de, gps)\n"
"-0.0016739778381789972\n");

static PyObject *pylal_XLALArrivalTimeDiff(PyObject *self, PyObject *args)
{
	PyObject *pos1, *pos2;
	double ra, dec;
	pylal_LIGOTimeGPS *gps;
	PyObject *result;

	/* 3-element list, 3-element list, float, float, LIGOTimeGPS */
	if(!PyArg_ParseTuple(args, "OOddO!:XLALArrivalTimeDiff", &pos1, &pos2, &ra, &dec, &pylal_LIGOTimeGPS_Type, &gps))
		return NULL;
	pos1 = PyArray_FromAny(pos1, PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_CONTIGUOUS, NULL);
	if(!pos1) {
		PyErr_SetString(PyExc_TypeError, "XLALArrivalTimeDiff() unable to convert argument 1 to array");
		return NULL;
	}
	pos2 = PyArray_FromAny(pos2, PyArray_DescrFromType(NPY_FLOAT64), 1, 1, NPY_CONTIGUOUS, NULL);
	if(!pos2) {
		PyErr_SetString(PyExc_TypeError, "XLALArrivalTimeDiff() unable to convert argument 2 to array");
		Py_DECREF(pos1);
		return NULL;
	}
	if(PyArray_DIM(pos1, 0) != 3) {
		PyErr_SetString(PyExc_TypeError, "XLALArrivalTimeDiff() argument 1 must have length 3");
		Py_DECREF(pos1);
		Py_DECREF(pos2);
	}
	if(PyArray_DIM(pos2, 0) != 3) {
		PyErr_SetString(PyExc_TypeError, "XLALArrivalTimeDiff() argument 2 must have length 3");
		Py_DECREF(pos1);
		Py_DECREF(pos2);
	}

	result = PyFloat_FromDouble(XLALArrivalTimeDiff(PyArray_DATA(pos1), PyArray_DATA(pos2), ra, dec, &gps->gps));
	Py_DECREF(pos1);
	Py_DECREF(pos2);

	/* float */
	return result;
}


/*
 * ============================================================================
 *
 *                            Module Registration
 *
 * ============================================================================
 */


static struct PyMethodDef module_methods[] = {
	{"XLALArrivalTimeDiff", pylal_XLALArrivalTimeDiff, METH_VARARGS, pylal_XLALArrivalTimeDiff__doc__},
	{"XLALGPSToUTC", pylal_XLALGPSToUTC, METH_VARARGS, pylal_XLALGPSToUTC__doc__},
	{"XLALGreenwichSiderealTime", pylal_XLALGreenwichSiderealTime, METH_VARARGS, pylal_XLALGreenwichSiderealTime__doc__},
	{"XLALUTCToGPS", pylal_XLALUTCToGPS, METH_VARARGS, pylal_XLALUTCToGPS__doc__},
	{NULL,}
};


PyMODINIT_FUNC initdate(void)
{
	/* commented out to silence warning */
	/*PyObject *module = */Py_InitModule3(MODULE_NAME, module_methods, "Wrapper for LAL's date package.");

	import_array();
	pylal_ligotimegps_import();
}
