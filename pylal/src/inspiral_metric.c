/*
    Copyright (C) 2012  Melissa Frei, Nickolas Fotopoulos

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <Python.h>

#include <string.h>
#include <lal/LALInspiralBank.h> /* where prototypes live and InspiralMomentsEtc and InspiralMetric */
#include <lal/LALInspiral.h> /* where PNOrder definition lives */
#include <lal/FrequencySeries.h> /* where the REAL8FrequencySeries lives */
#include <lal/LALDatatypes.h>

#include <misc.h>
#include <real8frequencyseries.h>

#define MODULE_NAME "pylal.inspiral_metric"

static PyObject *
inspiral_metric_getGammas(PyObject *self, PyObject *args)
{
    double fLower;
    double fCutoff;
    pylal_REAL8FrequencySeries *psd = NULL;
    int order; // Order must be between 0 and 7!
    double t0;
    double t3;

    InspiralMomentsEtc moments;
    InspiralMetric metric;

    memset(&moments, 0, sizeof(moments));
    memset(&metric, 0, sizeof(metric));

    if (!PyArg_ParseTuple(args, "ddiddO!", &fLower, &fCutoff, &order, &t0, &t3, &pylal_REAL8FrequencySeries_Type, &psd))
        return NULL;

    if (XLALGetInspiralMoments(&moments, fLower, fCutoff, psd->series)) {
        pylal_set_exception_from_xlalerrno();
        return NULL;
    }

    if (XLALInspiralComputeMetric(&metric, &moments, fLower, order, t0, t3)) {
        pylal_set_exception_from_xlalerrno();
        return NULL;
    }

    return Py_BuildValue("(ffffffffff)", metric.Gamma[0], metric.Gamma[1], metric.Gamma[2], metric.Gamma[3], metric.Gamma[4], metric.Gamma[5], metric.Gamma[6], metric.Gamma[7], metric.Gamma[8], metric.Gamma[9]);

}

static PyMethodDef Inspiral_MetricMethods[] = {
    {"compute_metric",  inspiral_metric_getGammas, METH_VARARGS, "Return the Gammas for given values of fLower, fCutoff, order, t0, t3, and psd. Gamma0 through Gamma5 are the upper diagonal values of the 2PN (tc, tau0, tau3) metric. Gamma6 through Gamma 9 were once used for various PTF metrics, but are return as 0 to match the SnglInspiralTable structure. Moments and metric are computed internally, so you can't reuse moments; this was done as a first cut to avoid wrapping the InspiralMomentsEtc structure."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initinspiral_metric(void)
{
    (void) Py_InitModule3("pylal.inspiral_metric", Inspiral_MetricMethods, "Compute (tc, tau0, tau3) metric coefficients");
    pylal_real8frequencyseries_import();
}
