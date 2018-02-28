#include <Python.h>
#include <string.h>
#include <stdlib.h>
#include <lal/XLALError.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/CoincInspiralEllipsoid.h>
#include <lal/Date.h>
#include <lal/EllipsoidOverlapTools.h>

static void GetAttrInPlaceString(char *dest, int n, PyObject *obj, char *attr)
{
    PyObject *val = PyObject_GetAttrString(obj, attr);
    PyObject *str;

    memset(dest, 0, n);

    if(!val)
        return;

    /* this extra step is done because the attribute is typically not a
     * string object.  it's either a unicode object or some kind of other
     * thing that is string-like but not actually a string (like an
     * ilwdchar subclass instance). */
    str = PyObject_Str(val);
    Py_DECREF(val);

    if(!str)
        return;

    strncpy(dest, PyString_AsString(str), n - 1);

    Py_DECREF(str);
}

static double GetAttrFloat(PyObject *obj, char *attr)
{
    PyObject *val = PyObject_GetAttrString(obj, attr);
    double result;

    if(!val)
        return XLAL_REAL8_FAIL_NAN;

    result = val == Py_None ? XLAL_REAL8_FAIL_NAN : PyFloat_AsDouble(val);

    Py_DECREF(val);

    return result;
}

static long GetAttrLong(PyObject *obj, char *attr)
{
    PyObject *val = PyObject_GetAttrString(obj, attr);
    long result;

    if(!val)
        return 0;

    result = val == Py_None ? 0 : PyInt_AsLong(val);

    Py_DECREF(val);

    return result;
}

static long long GetAttrLongLong(PyObject *obj, char *attr)
{
    PyObject *val = PyObject_GetAttrString(obj, attr);
    long long result;

    if(!val)
        return 0;

    result = val == Py_None ? 0 : PyLong_AsLongLong(val);

    Py_DECREF(val);

    return result;
}


SnglInspiralTable *PySnglInspiral2CSnglInspiral(PyObject *row) {
    /* Convert a Python SnglInspiral (row) to a C SnglInspiralTable.
    Used in function PyCalculateEThincaParameter and
      PyThincaParameterForInjection. */

    SnglInspiralTable *event; /* Return value */

    /* allocate new memory for row */
    event = calloc(1, sizeof(*event));
    event->event_id = calloc(1, sizeof(*event->event_id));

    /* copy to C SnglInspiral row */
    GetAttrInPlaceString(event->ifo, LIGOMETA_IFO_MAX, row, "ifo");
    GetAttrInPlaceString(event->search, LIGOMETA_SEARCH_MAX, row, "search");
    GetAttrInPlaceString(event->channel, LIGOMETA_CHANNEL_MAX, row, "channel");

    event->end_time.gpsSeconds = GetAttrLong(row, "end_time");
    event->end_time.gpsNanoSeconds = GetAttrLong(row, "end_time_ns");
    event->end_time_gmst = GetAttrFloat(row, "end_time_gmst");
    event->impulse_time.gpsSeconds = GetAttrLong(row, "impulse_time");
    event->impulse_time.gpsNanoSeconds = GetAttrLong(row, "impulse_time_ns");
    event->template_duration = GetAttrFloat(row, "template_duration");
    event->event_duration = GetAttrFloat(row, "event_duration");
    event->amplitude = GetAttrFloat(row, "amplitude");
    event->eff_distance = GetAttrFloat(row, "eff_distance");
    event->coa_phase = GetAttrFloat(row, "coa_phase");
    event->mass1 = GetAttrFloat(row, "mass1");
    event->mass2 = GetAttrFloat(row, "mass2");
    event->mchirp = GetAttrFloat(row, "mchirp");
    event->mtotal = GetAttrFloat(row, "mtotal");
    event->eta = GetAttrFloat(row, "eta");
    event->kappa = GetAttrFloat(row, "kappa");
    event->chi = GetAttrFloat(row, "chi");
    event->tau0 = GetAttrFloat(row, "tau0");
    event->tau2 = GetAttrFloat(row, "tau2");
    event->tau3 = GetAttrFloat(row, "tau3");
    event->tau4 = GetAttrFloat(row, "tau4");
    event->tau5 = GetAttrFloat(row, "tau5");
    event->ttotal = GetAttrFloat(row, "ttotal");
    event->psi0 = GetAttrFloat(row, "psi0");
    event->psi3 = GetAttrFloat(row, "psi3");
    event->alpha = GetAttrFloat(row, "alpha");
    event->alpha1 = GetAttrFloat(row, "alpha1");
    event->alpha2 = GetAttrFloat(row, "alpha2");
    event->alpha3 = GetAttrFloat(row, "alpha3");
    event->alpha4 = GetAttrFloat(row, "alpha4");
    event->alpha5 = GetAttrFloat(row, "alpha5");
    event->alpha6 = GetAttrFloat(row, "alpha6");
    event->beta = GetAttrFloat(row, "beta");
    event->f_final = GetAttrFloat(row, "f_final");
    event->snr = GetAttrFloat(row, "snr");
    event->chisq = GetAttrFloat(row, "chisq");
    event->chisq_dof = GetAttrLong(row, "chisq_dof");
    event->sigmasq = GetAttrFloat(row, "sigmasq");
    event->rsqveto_duration = GetAttrFloat(row, "rsqveto_duration");

    event->Gamma[0] = GetAttrFloat(row, "Gamma0");
    event->Gamma[1] = GetAttrFloat(row, "Gamma1");
    event->Gamma[2] = GetAttrFloat(row, "Gamma2");
    event->Gamma[3] = GetAttrFloat(row, "Gamma3");
    event->Gamma[4] = GetAttrFloat(row, "Gamma4");
    event->Gamma[5] = GetAttrFloat(row, "Gamma5");
    event->Gamma[6] = GetAttrFloat(row, "Gamma6");
    event->Gamma[7] = GetAttrFloat(row, "Gamma7");
    event->Gamma[8] = GetAttrFloat(row, "Gamma8");
    event->Gamma[9] = GetAttrFloat(row, "Gamma9");

    event->spin1x = GetAttrFloat(row, "spin1x");
    event->spin1y = GetAttrFloat(row, "spin1y");
    event->spin1z = GetAttrFloat(row, "spin1z");
    event->spin2x = GetAttrFloat(row, "spin2x");
    event->spin2y = GetAttrFloat(row, "spin2y");
    event->spin2z = GetAttrFloat(row, "spin2z");

    event->event_id->id = GetAttrLongLong(row, "event_id");

    if(PyErr_Occurred()) {
        free(event->event_id);
        free(event);
        return NULL;
    }

    return event;
}

SimInspiralTable *PySimInspiral2CSimInspiral(PyObject *row) {
    /* Convert a Python SimInspiral (row) to a C SimInspiralTable.
    Used in function PyThincaParameterForInjection. */

    SimInspiralTable *event; /* Return value */

    /* allocate new memory for row */
    event = calloc(1, sizeof(*event));
    event->event_id = calloc(1, sizeof(*event->event_id));

    /* copy to C SimInspiral row */
    GetAttrInPlaceString(event->waveform, LIGOMETA_WAVEFORM_MAX, row, "waveform");

    event->geocent_end_time.gpsSeconds = GetAttrLong(row, "geocent_end_time");
    event->geocent_end_time.gpsNanoSeconds = GetAttrLong(row, "geocent_end_time_ns");
    event->h_end_time.gpsSeconds = GetAttrLong(row, "h_end_time");
    event->h_end_time.gpsNanoSeconds = GetAttrLong(row, "h_end_time_ns");
    event->l_end_time.gpsSeconds = GetAttrLong(row, "l_end_time");
    event->l_end_time.gpsNanoSeconds = GetAttrLong(row, "l_end_time_ns");
    event->g_end_time.gpsSeconds = GetAttrLong(row, "g_end_time");
    event->g_end_time.gpsNanoSeconds = GetAttrLong(row, "g_end_time_ns");
    event->t_end_time.gpsSeconds = GetAttrLong(row, "t_end_time");
    event->t_end_time.gpsNanoSeconds = GetAttrLong(row, "t_end_time_ns");
    event->v_end_time.gpsSeconds = GetAttrLong(row, "v_end_time");
    event->v_end_time.gpsNanoSeconds = GetAttrLong(row, "v_end_time_ns");

    event->end_time_gmst = GetAttrFloat(row, "end_time_gmst");

    GetAttrInPlaceString(event->source, LIGOMETA_SOURCE_MAX, row, "source");

    event->mass1 = GetAttrFloat(row, "mass1");
    event->mass2 = GetAttrFloat(row, "mass2");
    event->eta = GetAttrFloat(row, "eta");
    event->distance = GetAttrFloat(row, "distance");
    event->longitude = GetAttrFloat(row, "longitude");
    event->latitude = GetAttrFloat(row, "latitude");
    event->inclination = GetAttrFloat(row, "inclination");
    event->coa_phase = GetAttrFloat(row, "coa_phase");
    event->polarization = GetAttrFloat(row, "polarization");
    event->psi0 = GetAttrFloat(row, "psi0");
    event->psi3 = GetAttrFloat(row, "psi3");
    event->alpha = GetAttrFloat(row, "alpha");
    event->alpha1 = GetAttrFloat(row, "alpha1");
    event->alpha2 = GetAttrFloat(row, "alpha2");
    event->alpha3 = GetAttrFloat(row, "alpha3");
    event->alpha4 = GetAttrFloat(row, "alpha4");
    event->alpha5 = GetAttrFloat(row, "alpha5");
    event->alpha6 = GetAttrFloat(row, "alpha6");
    event->beta = GetAttrFloat(row, "beta");
    event->spin1x = GetAttrFloat(row, "spin1x");
    event->spin1y = GetAttrFloat(row, "spin1y");
    event->spin1z = GetAttrFloat(row, "spin1z");
    event->spin2x = GetAttrFloat(row, "spin2x");
    event->spin2y = GetAttrFloat(row, "spin2y");
    event->spin2z = GetAttrFloat(row, "spin2z");
    event->theta0 = GetAttrFloat(row, "theta0");
    event->phi0 = GetAttrFloat(row, "phi0");
    event->f_lower = GetAttrFloat(row, "f_lower");
    event->f_final = GetAttrFloat(row, "f_final");
    event->mchirp = GetAttrFloat(row, "mchirp");
    event->eff_dist_h = GetAttrFloat(row, "eff_dist_h");
    event->eff_dist_l = GetAttrFloat(row, "eff_dist_l");
    event->eff_dist_g = GetAttrFloat(row, "eff_dist_g");
    event->eff_dist_t = GetAttrFloat(row, "eff_dist_t");
    event->eff_dist_v = GetAttrFloat(row, "eff_dist_v");

    GetAttrInPlaceString(event->event_id->textId, LIGOMETA_UNIQUE_MAX, row, "simulation_id");

    if(PyErr_Occurred()) {
        free(event->event_id);
        free(event);
        return NULL;
    }

    return event;
}

static PyObject *PyCalculateEThincaParameter(PyObject *self, PyObject *args) {
    /* Take two Python SnglInspiral values (rows of SnglInspiralTable) and
    call XLALCalculateEThincaParameter on their contents. */

    double result;
    PyObject *py_row1, *py_row2;
    SnglInspiralTable *c_row1, *c_row2;
    InspiralAccuracyList accuracyParams;

    if(!PyArg_ParseTuple(args, "OO", &py_row1, &py_row2))
        return NULL;

    /* Get rows into a format suitable for the LAL call */
    c_row1 = PySnglInspiral2CSnglInspiral(py_row1);
    if(!c_row1)
        return NULL;
    c_row2 = PySnglInspiral2CSnglInspiral(py_row2);
    if(!c_row2) {
        free(c_row1->event_id);
        free(c_row1);
        return NULL;
    }

    memset(&accuracyParams, 0, sizeof(accuracyParams));
    XLALPopulateAccuracyParams(&accuracyParams);

    /* This is the main call */
    result = XLALCalculateEThincaParameter(c_row1, c_row2, &accuracyParams);

    /* Free temporary memory */
    free(c_row1->event_id);
    free(c_row1);
    free(c_row2->event_id);
    free(c_row2);

    if(XLAL_IS_REAL8_FAIL_NAN(result)) {
        /* convert XLAL exception to Python exception */
        XLALClearErrno();
        PyErr_SetString(PyExc_ValueError, "SnglInspiral triggers are not coincident.");
        return NULL;
    }

    return PyFloat_FromDouble(result);
}

static PyObject *PyCalculateEThincaParameterExt(PyObject *self, PyObject *args) {
    /* Take two Python SnglInspiral values (rows of SnglInspiralTable) and
    call XLALCalculateEThincaParameter on their contents. */

    double result;
    double ra_deg, dec_deg;
    long gps;
    LIGOTimeGPS gpstime;
    PyObject    *py_row1, *py_row2;
    SnglInspiralTable *c_row1, *c_row2;
    InspiralAccuracyList accuracyParams;
    
    if(!PyArg_ParseTuple(args, "OOldd", &py_row1, &py_row2, &gps, &ra_deg, &dec_deg))
        return NULL;
    
    /* check the values */
    if (ra_deg<0 || ra_deg > 360) 
    {
      XLALPrintError("Right ascension value outside [0; 360]. Value given: %f\n", ra_deg);
      return NULL;
    }
    if (dec_deg<-90 || dec_deg>90) 
    {
      XLALPrintError("Declination value outside [-90; 90]. Value given: %f\n", dec_deg);
      return NULL;
    }

    /* Get rows into a format suitable for the LAL call */
    c_row1 = PySnglInspiral2CSnglInspiral(py_row1);
    if(!c_row1)
        return NULL;
    c_row2 = PySnglInspiral2CSnglInspiral(py_row2);
    if(!c_row2) {
        free(c_row1->event_id);
        free(c_row1);
        return NULL;
    }

    XLALGPSSet(&gpstime, gps, 0);
    memset(&accuracyParams, 0, sizeof(accuracyParams));
    XLALPopulateAccuracyParamsExt( &accuracyParams, &gpstime, ra_deg, dec_deg );

    /* This is the main call */    
    result = XLALCalculateEThincaParameter(c_row1, c_row2, &accuracyParams);

    /* Free temporary memory */
    free(c_row1->event_id);
    free(c_row1);
    free(c_row2->event_id);
    free(c_row2);

    if(XLAL_IS_REAL8_FAIL_NAN(result)) {
        /* convert XLAL exception to Python exception */
        XLALClearErrno();
        PyErr_SetString(PyExc_ValueError, "SnglInspiral triggers are not coincident.");
        return NULL;
    }

    return PyFloat_FromDouble(result);
}

static PyObject *PyEThincaParameterForInjection(PyObject *self, PyObject *args) {
    /* Take a Python SimInspiral value and a Python SnglInspiral value
    (rows of SimInspiralTable and SnglInspiralTable, respectively) and
    call XLALEThincaParameterForInjection on their contents. */
    
    double result;
    PyObject *py_row1, *py_row2;
    SimInspiralTable *c_row1;
    SnglInspiralTable *c_row2;
    
    if(!PyArg_ParseTuple(args, "OO", &py_row1, &py_row2))
        return NULL;
    
    /* Get rows into a format suitable for the LAL call */
    c_row1 = PySimInspiral2CSimInspiral(py_row1);
    if(!c_row1)
        return NULL;
    c_row2 = PySnglInspiral2CSnglInspiral(py_row2);
    if(!c_row2) {
        free(c_row1->event_id);
        free(c_row1);
        return NULL;
    }
    
    /* This is the main call */
    result = XLALEThincaParameterForInjection(c_row1, c_row2);
    
    /* Free temporary memory */
    free(c_row1->event_id);
    free(c_row1);
    free(c_row2->event_id);
    free(c_row2);
    
    return PyFloat_FromDouble(result);
}


static struct PyMethodDef tools_methods[] = {
    {"XLALCalculateEThincaParameter", PyCalculateEThincaParameter,
     METH_VARARGS,
     "XLALCalculateEThincaParameter(SnglInspiral1, SnglInspiral2)\n"
     "\n"
     "Takes two SnglInspiral objects (rows of a SnglInspiralTable) and\n"
     "calculates the overlap factor between them."},
    {"XLALCalculateEThincaParameterExt", PyCalculateEThincaParameterExt,
     METH_VARARGS,
     "XLALCalculateEThincaParameterExt(SnglInspiral1, SnglInspiral2, gps, ra_deg, dec_deg)\n"
     "\n"
     "Takes two SnglInspiral objects (rows of a SnglInspiralTable) and\n"
     "calculates the overlap factor between them, using the known time delay \n"
     "between the IFO's at a given time for a given sky location (given in degrees)."},
     {"XLALEThincaParameterForInjection", PyEThincaParameterForInjection,
      METH_VARARGS, \
      "XLALEThincaParameterForInjection(SimInspiral, SnglInspiral)\n"
      "\n"
      "Takes a SimInspiral and a SnglInspiral object (rows of\n"
      "SimInspiralTable and SnglInspiralTable, respectively) and\n"
      "calculates the ethinca parameter required to put the SimInspiral\n"},   
    {NULL, NULL, 0}
};

void inittools (void) {
    (void) Py_InitModule("pylal.tools", tools_methods);
}
