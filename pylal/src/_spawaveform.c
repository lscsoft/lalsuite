#include <Python.h>

/* standard includes */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>

/* LAL Includes */

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/LALInspiral.h>

/* GSL includes */
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <numpy/arrayobject.h>


/* static functions used by the python wrappers */
static int SPAWaveform (double mass1, double mass2, int order, double deltaF, double deltaT, double fLower, double fFinal, int numPoints, complex double *expPsi);
static double chirp_time (double m1, double m2, double fLower, int order,double chi);
static double chirp_time_between_f1_and_f2(double m1, double m2, double fLower, double fUpper, int order, double chi);
static double schwarz_isco(double m1, double m2);
static double bkl_isco(double m1, double m2);
static double light_ring(double m1, double m2);
static int IMRSPAWaveform(double mass1, double mass2, double spin1,  double spin2, double deltaF, double fLower, int numPoints, complex double *hOfF);
static int SPAWaveformReduceSpin (double mass1, double mass2, double chi, int order, double startTime, double phi0, double deltaF, double fLower, double fFinal, int numPoints, complex double *hOfF);
static int IMRSPAWaveformFromChi(double mass1, double mass2, double chi, double deltaF, double fLower, int numPoints, complex double *hOfF);
static int GenericSPAWaveform (double *psis, double *psils, int order, double deltaF, double deltaT, double fLower, double fFinal, int numPoints, complex double *expPsi);
static double imr_merger(double m1, double m2, double chi);
static double imr_ring(double m1, double m2, double chi);
static double imr_fcut(double m1, double m2, double chi);
static double compute_chi(double m1, double m2, double spin1, double spin2);

/* doc string for help() */
const char SPADocstring[] =
"This module wraps SPA inspiral waveform generation and some useful functions "
"related to them.\n\n"
"EXAMPLE CODE: (you could cut and paste this into the interpreter)\n"
"\n"
"from pylal import spawaveform\n"
"import numpy\n"
"import pylab\n"
"import scipy\n"
"import time\n"
"\n"
"# some parameters\n"
"m1 = 25.0\n"
"m2 = 25.0\n"
"order = 7\n"
"fLower = 30.0\n"
"fFinal = spawaveform.ffinal(m1,m2,'light_ring')\n"
"# find next highest powers of 2\n"
"sr = 2**numpy.ceil(numpy.log2(fFinal*2))\n"
"dur = 2**numpy.ceil(numpy.log2(spawaveform.chirptime(m1,m2,order,fLower)))\n"
"deltaF = 1.0 / dur\n"
"deltaT = 1.0 / sr\n"
"\n"
"pylab.figure(1)\n"
"for j, endfreq in enumerate(['schwarz_isco','bkl_isco','light_ring']):\n"
"\tfFinal = spawaveform.ffinal(m1,m2,endfreq)\n"
"\tprint endfreq, fFinal\n"
"\t# time stamp vector\n"
"\ttimestamp = numpy.arange(dur * sr) / sr\n"
"\t# initialise data for the template\n"
"\tz = numpy.empty(sr * dur, 'complex128')\n"
"\t# make a spawaveform\n"
"\tspawaveform.waveform(m1, m2, order, deltaF, deltaT, fLower, fFinal, z)\n"
"\tz = scipy.ifft(z)\n"
"\t# plot it\n"
"\tpylab.subplot(3,1,j+1)\n"
"\tpylab.plot(timestamp, numpy.real(z))\n"
"\tpylab.hold(1)\n"
"\tpylab.plot(timestamp, numpy.imag(z),'r')\n"
"\tpylab.legend([endfreq + ' hc(t)', endfreq + ' hs(t)'],loc='lower left')\n"
"\tpylab.xlim([dur - spawaveform.chirptime(m1,m2,order,40.0,fFinal), dur])\n"
"pylab.hold(0)\n"
"pylab.show()\n";

/* FIXME someday do better error handling? */
/* static PyObject* SPAWaveformError;      */

/* Function to return the final frequency of spa waveforms */
static PyObject *PyFFinal(PyObject *self, PyObject *args)
	{
	double mass1, mass2;
	const char *s = NULL;
	if(!PyArg_ParseTuple(args, "dd|s", &mass1, &mass2, &s)) return NULL;
	/* Default is schwarz isco */
	if ( !s || !strcmp(s, "schwarz_isco") ) return Py_BuildValue("d", schwarz_isco(mass1,mass2));
	if ( s && !strcmp(s, "bkl_isco") ) return Py_BuildValue("d", bkl_isco(mass1,mass2));
	if ( s && !strcmp(s, "light_ring") ) return Py_BuildValue("d", light_ring(mass1,mass2));
	PyErr_SetString(PyExc_ValueError, "Unrecognized ending frequency, must be schwarz_isco | bkl_isco | light_ring");
	return NULL;
	}

/* Function to return the characteristic frequencies of IMR waveforms */
static PyObject *PyIMRFFinal(PyObject *self, PyObject *args)
	{
	double mass1, mass2, chi;
	const char *s = NULL;
	if(!PyArg_ParseTuple(args, "ddd|s", &mass1, &mass2, &chi, &s)) return NULL;
	/* Default is schwarz isco */
	if ( s && !strcmp(s, "merger") ) return Py_BuildValue("d",imr_merger(mass1,mass2,chi));
	if (!s || !strcmp(s, "fcut") ) return Py_BuildValue("d", imr_fcut(mass1,mass2,chi));
	if ( s && !strcmp(s, "ringdown") ) return Py_BuildValue("d", imr_ring(mass1,mass2,chi));
	PyErr_SetString(PyExc_ValueError, "Unrecognized ending frequency, must be merger | ringdown | fcut");
	return NULL;
	}

/* Function to compute the "chi" parameter (mass weighted combined spin) */
static PyObject *PyComputeChi(PyObject *self, PyObject *args)
	{
	double mass1, mass2, spin1, spin2;
	if(!PyArg_ParseTuple(args, "dddd", &mass1, &mass2, &spin1, &spin2)) return NULL;
	return Py_BuildValue("d",compute_chi(mass1, mass2, spin1, spin2));
	}

/* Function to compute the frequency domain SPA waveform */
static PyObject *PySPAWaveform(PyObject *self, PyObject *args)
	{
	/* Generate a SPA (frequency domain) waveform at a given PN order */
	PyObject *arg9, *py_spa_array;
	double mass1, mass2, deltaF, deltaT, fLower, fFinal;
	int order;
	npy_intp *dims = NULL;
	complex double *data = NULL;
	double spin1 = -100.0;
	double spin2 = -100.0;
	double chi = 0.0;

	/* FIXME properly handle references */
	if(!PyArg_ParseTuple(args, "ddiddddO|dd", &mass1, &mass2, &order, &deltaF, &deltaT, &fLower, &fFinal, &arg9, &spin1, &spin2)) return NULL;
	/* this gets a contiguous memory numpy array */
        py_spa_array = PyArray_FROM_OTF(arg9, NPY_CDOUBLE, NPY_IN_ARRAY);
	if (py_spa_array == NULL) return NULL;
	/* Actually call the SPA waveform C function */
	/* FIXME no checking of the array dimensions, this could be done in a python wrapper */
	dims = PyArray_DIMS(py_spa_array);
	data = PyArray_DATA(py_spa_array);
	if (spin1 == -100.0 && spin2 == -100.0)
		SPAWaveform(mass1, mass2, order, deltaF, deltaT, fLower, fFinal, dims[0], data);
	if (spin1 != -100.0 && spin2 == -100.0) {
		chi = spin1; // only one spin argument is interpreted as chi
		SPAWaveformReduceSpin(mass1, mass2, chi, order, 0.0, 0.0, deltaF, fLower, fFinal, dims[0], data);
		}
	if (spin1 != -100.0 && spin2 != -100.0)	{
		chi = compute_chi(mass1, mass2, spin1, spin2);
		SPAWaveformReduceSpin(mass1, mass2, chi, order, 0.0, 0.0, deltaF, fLower, fFinal, dims[0], data);
		}
	Py_DECREF(py_spa_array);
        Py_INCREF(Py_None);
        return Py_None;
	}

/* Function to compute the frequency domain IMR waveform */
static PyObject *PyIMRSPAWaveform(PyObject *self, PyObject *args)
	{
	/* Generate a SPA (frequency domain) waveform at a given PN order */
	PyObject *arg9, *py_spa_array;
	double mass1, mass2, deltaF, fLower;
	npy_intp *dims = NULL;
	complex double *data = NULL;
	/*FIXME get rid of this hack to handle optional spin arguments */
	double spin1 = -100.0;
	double spin2 = -100.0;
	/* FIXME properly handle references */
	if (!PyArg_ParseTuple(args, "ddddO|dd", &mass1, &mass2, &deltaF, &fLower, &arg9, &spin1, &spin2)) return NULL;
	/* Check for no spin case */
	if (spin1 == -100.0 && spin2 == -100.0) spin1 = spin2 = 0.0;
	/* this gets a contiguous memory numpy array */
        py_spa_array = PyArray_FROM_OTF(arg9, NPY_CDOUBLE, NPY_IN_ARRAY);
	if (py_spa_array == NULL) return NULL;
	/* Actually call the SPA waveform C function */
	/* FIXME no checking of the array dimensions, this could be done in a python wrapper */
	dims = PyArray_DIMS(py_spa_array);
	data = PyArray_DATA(py_spa_array);
	/* depending on the number of arguments given call a different function */
	if (spin1 != -100.0 && spin2 == -100.0) IMRSPAWaveformFromChi(mass1, mass2, spin1, deltaF, fLower, dims[0], data);
	else IMRSPAWaveform(mass1, mass2, spin1, spin2, deltaF, fLower, dims[0], data);
	Py_DECREF(py_spa_array);
        Py_INCREF(Py_None);
        return Py_None;
	}

/* Function to compute the frequency domain generic inspiral waveform */
static PyObject *PyGenericSPAWaveform(PyObject *self, PyObject *args)
	{
	/* Generate a generic SPA (frequency domain) waveform of given PN order */
	PyObject *arg1, *py_psi_array;
	PyObject *arg2, *py_psil_array;
	PyObject *arg8, *py_spa_array;
	double deltaF, deltaT, fLower, fFinal;
	int order;
	npy_intp *psidims = NULL;
	double *psis = NULL;
	npy_intp *psildims = NULL;
	double *psils = NULL;
	npy_intp *dims = NULL;
	complex double *data = NULL;

	/* FIXME properly handle references */
	if(!PyArg_ParseTuple(args, "OOiddddO", &arg1, &arg2, &order, &deltaF, &deltaT, &fLower, &fFinal, &arg8)) return NULL;
	/* this gets a contiguous memory numpy arrays */
        py_psi_array = PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_IN_ARRAY);
	if (py_psi_array == NULL) return NULL;
        py_psil_array = PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_IN_ARRAY);
	if (py_psil_array == NULL) return NULL;
        py_spa_array = PyArray_FROM_OTF(arg8, NPY_CDOUBLE, NPY_IN_ARRAY);
	if (py_spa_array == NULL) return NULL;

	/* get PN coefficients and check they match the order requested */
	psidims = PyArray_DIMS(py_psi_array);
	psis = PyArray_DATA(py_psi_array);
	psildims = PyArray_DIMS(py_psil_array);
	psils = PyArray_DATA(py_psil_array);
	if (order != psidims[0]-1) return NULL;
	if (order != psildims[0]-1) return NULL;

	/* Actually call the SPA waveform C function */
	/* FIXME no checking of the array dimensions, this could be done in a python wrapper */
	dims = PyArray_DIMS(py_spa_array);
	data = PyArray_DATA(py_spa_array);
	GenericSPAWaveform(psis, psils, order, deltaF, deltaT, fLower, fFinal, dims[0], data);
	Py_DECREF(py_psi_array);
	Py_DECREF(py_psil_array);
	Py_DECREF(py_spa_array);
        Py_INCREF(Py_None);
        return Py_None;
	}

/* Function to wrap GSLs SVD */
static PyObject *PySVD(PyObject *self, PyObject *args, PyObject *keywds)
        {

        /*
	 * input A (m x n matrix)
	 * output U,S,V
	 */

	/* input array */
	PyObject *a, *A;

	/* output data  */
	PyObject *U, *V, *S, *out;

	/* views of output data */
	gsl_matrix_view gU, gV;
	gsl_vector_view gS;

	/* array sizes */
	npy_intp *Adims = NULL;
	npy_intp Vdims[] = {0.0, 0.0};
	npy_intp Sdims[] = {0.0};

	/* native C-datatype representation of numpy array data */
	double *cA = NULL;
	double *cU = NULL;
	double *cV = NULL;
	double *cS = NULL;

	/* workspace data */
	gsl_vector *gW;
	gsl_matrix *gX;

	/* list of available keywords */
	char *kwlist[] = {"array","inplace","mod",NULL};

	/* keyword argument vars */
	int *modGolubReinsch = 0;
	int *inplace = 0;

	/*
	 * end declarations
	 */

	/* Read in input array, represent in a,A,cA,Adims */
	if(!PyArg_ParseTupleAndKeywords(args, keywds, "O|ii", kwlist, &a, &inplace, &modGolubReinsch)) return NULL;

	A = PyArray_FROM_OTF(a, NPY_DOUBLE, NPY_IN_ARRAY);
	Adims = PyArray_DIMS(A);
	cA = PyArray_DATA(A);

	/* Don't support M < N */
	if ( Adims[0] < Adims[1] ) return NULL;

	/* allocate new ouput matrices */
	if ( inplace )
	{
	    U = A;
	    cU = cA;
	}
	else
	{
	    U = PyArray_SimpleNew(2, Adims, NPY_DOUBLE);
	    cU = PyArray_DATA(U);
	    memcpy( cU, cA, Adims[0]*Adims[1]*sizeof(double) );
	}

	Vdims[0] = Vdims[1] = Adims[1];
	V = PyArray_SimpleNew(2, Vdims, NPY_DOUBLE);
	cV = PyArray_DATA(V);

	Sdims[0] = Adims[1];
	S = PyArray_SimpleNew(1, Sdims, NPY_DOUBLE);
	cS = PyArray_DATA(S);

	/* Allocate workspace gsl matrix */
	gW = gsl_vector_calloc(Adims[1]);

	/* Get gsl matrix views of the numpy array data */
	/* U will be overwritten */
	gU = gsl_matrix_view_array(cU, Adims[0], Adims[1]);
	gS = gsl_vector_view_array(cS, Sdims[0]);
	gV = gsl_matrix_view_array(cV, Vdims[0], Vdims[1]);

	/* PERFORM THE SVD! */
	if ( modGolubReinsch )
	{
	    gX = gsl_matrix_calloc(Adims[1],Adims[1]);
	    gsl_linalg_SV_decomp_mod (&(gU.matrix),gX, &(gV.matrix), &(gS.vector), gW);
	    gsl_matrix_free(gX);
	}
	else
	{
	    gsl_linalg_SV_decomp (&(gU.matrix), &(gV.matrix), &(gS.vector), gW);
	}
	gsl_vector_free(gW);

	/* take the transpose to be consistent with scipys svd */
	gsl_matrix_transpose(&(gV.matrix));
	out = Py_BuildValue("OOO", U, S, V);

	/* FIXME HELP I don't know what to do about references... */
	/*Py_DECREF(A); */
	/*Py_DECREF(V); */
	/*Py_DECREF(S); */
        /*Py_INCREF(Py_None); */

	return out;
	}

/* Function to compute chirp time */
static PyObject *PyChirpTime(PyObject *self, PyObject *args)
	{
	/* Generate a SPA (frequency domain) waveform at a given PN order */
	double mass1, mass2, fLower, fFinal, time, chi;
	int order;
	fFinal = 0;
	chi = 0.;
	/* FIXME properly handle references */
	if (!PyArg_ParseTuple(args, "ddid|dd", &mass1, &mass2, &order, &fLower, &fFinal, &chi)) return NULL;
	if (fFinal)
		time = chirp_time_between_f1_and_f2(mass1, mass2, fLower, fFinal, order, chi);
	else
		time = chirp_time(mass1, mass2, fLower, order, chi);
	return Py_BuildValue("d", time);
	}

static PyObject *PyIIR(PyObject *self, PyObject *args)
{
	PyObject *amp, *phase;
	PyObject *amp_array, *phase_array;
	double eps, alpha, beta, padding;
	REAL8Vector amp_real8, phase_real8;
	COMPLEX16Vector *a1 =NULL;
	COMPLEX16Vector *b0 = NULL;
	INT4Vector *delay =NULL;
	PyObject *a1_pyob, *b0_pyob, *delay_pyob;
	npy_intp *amp_arraydims = NULL;
	npy_intp *phase_arraydims = NULL;
	npy_intp a1_length[] = {0};
	npy_intp b0_length[] = {0};
	npy_intp delay_length[] = {0};
	PyObject *out;

	if (!PyArg_ParseTuple(args, "OOdddd", &amp, &phase, &eps, &alpha, &beta, &padding)) return NULL;
	amp_array = PyArray_FROM_OTF(amp, NPY_DOUBLE, NPY_IN_ARRAY);
	phase_array = PyArray_FROM_OTF(phase, NPY_DOUBLE, NPY_IN_ARRAY);

	amp_arraydims = PyArray_DIMS(amp_array);
	amp_real8.length = amp_arraydims[0];
	amp_real8.data = PyArray_DATA(amp_array);
	phase_arraydims = PyArray_DIMS(phase_array);
	phase_real8.length = phase_arraydims[0];
	phase_real8.data = PyArray_DATA(phase_array);
	XLALInspiralGenerateIIRSet(&amp_real8, &phase_real8, eps, alpha, beta, padding, &a1, &b0, &delay);
	a1_length[0] = a1->length;
	b0_length[0] = b0->length;
	delay_length[0] = delay->length;
	a1_pyob = PyArray_SimpleNewFromData(1, a1_length, NPY_CDOUBLE, (void *) a1->data);
	b0_pyob = PyArray_SimpleNewFromData(1, b0_length, NPY_CDOUBLE, (void *) b0->data);
	delay_pyob = PyArray_SimpleNewFromData(1, delay_length, NPY_INT, (void *) delay->data);
	out = Py_BuildValue("OOO", a1_pyob, b0_pyob, delay_pyob);
	Py_DECREF(amp_array);
	Py_DECREF(phase_array);
	return out;
}

static PyObject *PyIIRResponse(PyObject *self, PyObject *args)
{
	PyObject *a1, *b0, *delay;
	PyObject *a1_array, *b0_array, *delay_array;
	int N;
	COMPLEX16Vector a1_complex16, b0_complex16;
	INT4Vector delay_int4;
	COMPLEX16Vector *resp =NULL;
	PyObject *resp_pyob;
	npy_intp *a1_arraydims = NULL;
	npy_intp *b0_arraydims = NULL;
	npy_intp *delay_arraydims = NULL;
	npy_intp resp_length[] = {0};

	if (!PyArg_ParseTuple(args, "iOOO", &N, &a1, &b0, &delay)) return NULL;
	a1_array = PyArray_FROM_OTF(a1, NPY_CDOUBLE, NPY_IN_ARRAY);
	b0_array = PyArray_FROM_OTF(b0, NPY_CDOUBLE, NPY_IN_ARRAY);
	delay_array = PyArray_FROM_OTF(delay, NPY_INT, NPY_IN_ARRAY);
	a1_arraydims = PyArray_DIMS(a1_array);
	a1_complex16.length = a1_arraydims[0];
	a1_complex16.data = PyArray_DATA(a1_array);
	b0_arraydims = PyArray_DIMS(b0_array);
	b0_complex16.length = b0_arraydims[0];
	b0_complex16.data = PyArray_DATA(b0_array);
	delay_arraydims = PyArray_DIMS(delay_array);
	delay_int4.length = delay_arraydims[0];
	delay_int4.data = PyArray_DATA(delay_array);

	resp = XLALCreateCOMPLEX16Vector(N);

	XLALInspiralIIRSetResponse(&a1_complex16, &b0_complex16, &delay_int4, resp);
	resp_length[0] = resp->length;
	resp_pyob = PyArray_SimpleNew(1, resp_length, NPY_CDOUBLE);
	memcpy(PyArray_DATA(resp_pyob), resp->data, resp->length * sizeof(*resp->data));
		
	XLALDestroyCOMPLEX16Vector(resp);
	Py_DECREF(a1_array);
	Py_DECREF(b0_array);
	Py_DECREF(delay_array);

	return resp_pyob;
}

static PyObject *PyIIRInnerProduct(PyObject *self, PyObject *args)
{
	PyObject *a1, *b0, *delay, *psd;
	REAL8 ip = 0.0;
	PyObject *a1_array, *b0_array, *delay_array, *psd_array;
	COMPLEX16Vector a1_complex16, b0_complex16;
	INT4Vector delay_int4;
	REAL8Vector psd_real8;
	npy_intp *a1_arraydims = NULL;
	npy_intp *b0_arraydims = NULL;
	npy_intp *delay_arraydims = NULL;
	npy_intp *psd_arraydims = NULL;

	if (!PyArg_ParseTuple(args, "OOOO", &a1, &b0, &delay, &psd)) return NULL;
	a1_array = PyArray_FROM_OTF(a1, NPY_CDOUBLE, NPY_IN_ARRAY);
	b0_array = PyArray_FROM_OTF(b0, NPY_CDOUBLE, NPY_IN_ARRAY);
	delay_array = PyArray_FROM_OTF(delay, NPY_INT, NPY_IN_ARRAY);
	psd_array = PyArray_FROM_OTF(psd, NPY_DOUBLE, NPY_IN_ARRAY);
	a1_arraydims = PyArray_DIMS(a1_array);
	a1_complex16.length = a1_arraydims[0];
	a1_complex16.data = PyArray_DATA(a1_array);
	b0_arraydims = PyArray_DIMS(b0_array);
	b0_complex16.length = b0_arraydims[0];
	b0_complex16.data = PyArray_DATA(b0_array);
	delay_arraydims = PyArray_DIMS(delay_array);
	delay_int4.length = delay_arraydims[0];
	delay_int4.data = PyArray_DATA(delay_array);
	psd_arraydims = PyArray_DIMS(psd_array);
	psd_real8.length = psd_arraydims[0];
	psd_real8.data = PyArray_DATA(psd_array);
	XLALInspiralCalculateIIRSetInnerProduct(&a1_complex16, &b0_complex16, &delay_int4, &psd_real8, &ip);
	Py_DECREF(a1_array);
	Py_DECREF(b0_array);
	Py_DECREF(delay_array);
	Py_DECREF(psd_array);

	return Py_BuildValue("d", ip);
}

/* Structure defining the functions of this module and doc strings etc... */
static struct PyMethodDef methods[] = {
	{"waveform", PySPAWaveform, METH_VARARGS,
	 "This function produces a frequency domain waveform at a "
	 "specified mass1, mass2 and PN order.\n\n"
	 "waveform(m1, m2, order, deltaF, deltaT, fLower, fFinal, signalArray)\n\n"
	 "You can produce a spin aligned waveform by doing\n\n"
	 "waveform(m1, m2, order, deltaF, deltaT, fLower, fFinal, signalArray, spin1, spin2)"
	 "Or you can produce a spin aligned waveform by doing\n\n"
	 "waveform(m1, m2, order, deltaF, deltaT, fLower, fFinal, signalArray, chi)"
	},
	{"genericwaveform", PyGenericSPAWaveform, METH_VARARGS,
	 "This function produces a frequency domain waveform at a "
	 "specified PN order using user defined PN coefficients.\n\n"
	 "genericwaveform(psis, psils, order, deltaF, deltaT, fLower, fFinal, signalArray)"
	},
	{"imrwaveform", PyIMRSPAWaveform, METH_VARARGS,
	 "This function produces a frequency domain IMR waveform at a "
	 "specified mass1, mass2 by calling \n\n"
	 "imrwaveform(m1, m2, deltaF, fLower, signalArray)\n\n"
	 "This function is overlouded to produce a frequency domain IMR waveform at a "
	 "specified mass1, mass2, chi by calling \n\n"
	 "imrwaveform(m1, m2, deltaF, fLower, signalArray, chi)\n\n"
	 "This function is overlouded to produce a frequency domain IMR waveform at a "
	 "specified mass1, mass2, z component spin1, z component spin2 by calling \n\n"
	 "imrwaveform(m1, m2, deltaF, fLower, signalArray, spin1, spin2)\n\n"
	},
	{"chirptime", PyChirpTime, METH_VARARGS,
	 "This function calculates the SPA chirptime at a specified mass1, mass2 "
	 "and PN order between two frequencies.  If the second frequency is omitted "
	 "it is assumed to be infinite. To include spin (chi) you can place it as " 
	 "the last argument, but you must give an fFinal.\n\n"
	 "chirptime(m1, m2, order, fLower, [fFinal, chi])\n\n"
	},
	{"ffinal", PyFFinal, METH_VARARGS,
	 "This function calculates the ending frequency specified by "
	 "mass1 and mass2.\n\n"
	 "ffinal(m1, m2, ['schwarz_isco'|'bkl_isco'|'light_ring'])\n\n"
	},
	{"imrffinal", PyIMRFFinal, METH_VARARGS,
	 "This function calculates the ending frequency specified by "
	 "mass1 and mass2 and chi for the imr waveforms. The default is "
	 "to return the cutoff frequency unless you specify merger or ringdown\n\n"
	 "ffinal(m1, m2, chi, ['merger'|'ringdown'])\n\n"
	},
	{"computechi", PyComputeChi, METH_VARARGS,
	 "This function calculates the mass weighted spin parameter chi\n\n"
	 "computechi(m1, m2, spin1, spin2)\n\n"
	},
	{"svd", (PyCFunction) PySVD, METH_KEYWORDS,
	 "This function calculates the singular value decomposition of a matrix\n"
	 "via the GSL implementation of the Golub-Reinsch algorithm.  The default\n"
	 "GSL function used is\n"
	 "gsl_linalg_SV_decomp (gsl_matrix * A, gsl_matrix * V, gsl_vector * S, gsl_vector * work)\n\n"
	 "USAGE:\n\tU, S, V = svd(A,mod=False,inplace=False)\n\n"
	 "A is an MxN numpy array, V is an Nx1 numpy array and V is an MxM numpy array.\n"
	 "The case M<N is not supported by GSL and this wrapping does not add support for this.\n"
	 "If mod=True, then the modified Golub-Reinsch algorithm is used (implemented in GSL\n"
	 "as gsl_linalg_SV_decomp_mod). This algorithm uses slightly more memory but is expected\n"
	 "to out-perform the standard Golub-Reinsch the limit M>>N.\n"
	 "If the inplace=True, the input array A is overwritten by U, instead of the default\n"
	 "behavior which is to allocate new space for U and preserve A.  Both variables continue\n"
	 "to exist but point to the same data! USE THIS OPTION WITH CARE!\n\n"
	 "EXAMPLE:\n\tfrom pylal import spawaveform\n"
	 "\timport numpy\n"
	 "\tA = numpy.random.randn(4,3)\n"
	 "\tprint A\n"
	 "\tU,S,V = spawaveform.svd(A)\n"
	 "\tB = U * S\n"
	 "\tAprime = numpy.dot(B,numpy.transpose(V))\n"
	 "\tprint Aprime\n\n"
	},
	{"iir", PyIIR, METH_VARARGS,
	 "This function calculates the a set of single pole IIR filters corresponding to a complex\n"
	 "time series\n\n"
	 "ffinal(amplitude, phase, epsilion, alpha, beta, a1, b0, delay\n\n"
	},
	{"iirresponse", PyIIRResponse, METH_VARARGS,
	 "This function produces a truncated Impulse response function for a given delay of IIR filters with a1, b0 and delays.\n\n"
	 "iirresponse(length_of_impulse_response, a1_set, b0_set, delay_set\n\n"
	},
	{"iirinnerproduct", PyIIRInnerProduct, METH_VARARGS,
	 "This function outputs the inner product of the sum of iir responses\n\n"
	 "iirinnerproduct(a1_set, b0_set, delay_set, psd\n\n"
	},
	{NULL, NULL, 0, NULL}
	};

/* The init function for this module */
void init_spawaveform(void)
	{
	(void) Py_InitModule3("pylal._spawaveform", methods, SPADocstring);
	import_array();
	/* FIXME someday handle errors
	 * SVMError = PyErr_NewException("_spawaveform.SPAWaveformError", NULL, NULL);
	 * Py_INCREF(SPAWaveformError);
	 * PyModule_AddObject(m, "SPAWaveformError", SPAWaveformError);
         */
	}

/*****************************************************************************/
/* The remainder of this code defines the static functions that the python   */
/* functions will use to compute various quantities.  They are not exposed   */
/* outside of this file directly but can be called from python via           */
/* the documentation described when doing help() on this module.             */
/* A lot of this code is in lal in some form or must be moved to lal.  Here  */
/* the functin prototypes have been vastly simplified (use of native c types */
/* and double precision)                                                     */
/*****************************************************************************/

static int SPAWaveformReduceSpin (double mass1, double mass2, double chi, 
        int order, double startTime, double phi0, double deltaF,
        double fLower, double fFinal, int numPoints, complex double *hOfF) {

	double m = mass1 + mass2;
	double eta = mass1 * mass2 / m / m;
    
    double psi2 = 0., psi3 = 0., psi4 = 0., psi5 = 0., psi6 = 0., psi6L = 0., psi7 = 0.; 
    double psi3S = 0., psi4S = 0., psi5S = 0., psi0; 
    double alpha2 = 0., alpha3 = 0., alpha4 = 0., alpha5 = 0., alpha6 = 0., alpha6L = 0.;
    double alpha7 = 0., alpha3S = 0., alpha4S = 0., alpha5S = 0.; 
    double f, v, v2, v3, v4, v5, v6, v7, Psi, amp, shft, amp0, d_eff; 
    int k, kmin, kmax; 

    double mSevenBySix = -7./6.;
    double piM = LAL_PI*m*LAL_MTSUN_SI;
    double oneByThree = 1./3.;
    double piBy4 = LAL_PI/4.;

    /************************************************************************/
    /* spin terms in the ampl & phase in terms of the 'reduced-spin' param. */
    /************************************************************************/
    psi3S = 113.*chi/3.;
    psi4S = 63845.*(-81. + 4.*eta)*chi*chi/(8.*pow(-113. + 76.*eta, 2.));  
    psi5S = -565.*(-146597. + 135856.*eta + 17136.*eta*eta)*chi/(2268.*(-113. + 76.*eta)); 

    alpha3S = (113*chi)/24.; 
    alpha4S = (12769*pow(chi,2)*(-81 + 4*eta))/(32.*pow(-113 + 76*eta,2)); 
    alpha5S = (-113*chi*(502429 - 591368*eta + 1680*pow(eta,2)))/(16128.*(-113 + 76*eta)); 

    /* coefficients of the phase at PN orders from 0 to 3.5PN */
    psi0 = 3./(128.*eta);

    /************************************************************************/
    /* set the amplitude and phase coefficients according to the PN order   */
    /************************************************************************/
    switch (order) {
        case 7: 
            psi7 = (77096675.*LAL_PI)/254016. + (378515.*LAL_PI*eta)/1512.  
                     - (74045.*LAL_PI*eta*eta)/756.;
            alpha7 = (-5111593*LAL_PI)/2.709504e6 - (72221*eta*LAL_PI)/24192. - 
                        (1349*pow(eta,2)*LAL_PI)/24192.; 
        case 6:
            psi6 = 11583231236531./4694215680. - (640.*LAL_PI*LAL_PI)/3. - (6848.*LAL_GAMMA)/21. 
                     + (-5162.983708047263 + 2255.*LAL_PI*LAL_PI/12.)*eta 
                     + (76055.*eta*eta)/1728. - (127825.*eta*eta*eta)/1296.;
            psi6L = -6848./21.;
            alpha6 = -58.601030974347324 + (3526813753*eta)/2.7869184e7 - 
                        (1041557*pow(eta,2))/258048. + (67999*pow(eta,3))/82944. + 
                        (10*pow(LAL_PI,2))/3. - (451*eta*pow(LAL_PI,2))/96.; 
            alpha6L = 856/105.; 
        case 5:
            psi5 = (38645.*LAL_PI/756. - 65.*LAL_PI*eta/9. + psi5S);
            alpha5 = (-4757*LAL_PI)/1344. + (57*eta*LAL_PI)/16. + alpha5S; 
        case 4:
            psi4 = 15293365./508032. + 27145.*eta/504. + 3085.*eta*eta/72. + psi4S;
            alpha4 = 0.8939214212884228 + (18913*eta)/16128. + (1379*pow(eta,2))/1152. + alpha4S; 
        case 3:
            psi3 = psi3S - 16.*LAL_PI;
            alpha3 = -2*LAL_PI + alpha3S; 
        case 2:
            psi2 = 3715./756. + 55.*eta/9.;
            alpha2 = 1.1056547619047619 + (11*eta)/8.; 
        default:
            break;
    }

    /* compute the amplitude assuming effective dist. of 1 Mpc */
    d_eff = 1e6*LAL_PC_SI/LAL_C_SI;  /*1 Mpc in seconds */
    amp0 = sqrt(5.*eta/24.)*pow(m*LAL_MTSUN_SI, 5./6.)/(d_eff*pow(LAL_PI, 2./3.));

    shft = 2.*LAL_PI *startTime;

    /* zero outout */    
    memset (hOfF, 0, numPoints * sizeof (complex double));

	kmin = fLower / deltaF > 1 ? fLower / deltaF : 1;
	kmax = fFinal / deltaF < numPoints  ? fFinal / deltaF : numPoints ;

    /************************************************************************/
    /*          now generate the waveform at all frequency bins             */
    /************************************************************************/
    for (k = kmin; k < kmax; k++) {

        /* fourier frequency corresponding to this bin */
      	f = k * deltaF;
        v = pow(piM*f, oneByThree);

        v2 = v*v;   v3 = v2*v;  v4 = v3*v;  v5 = v4*v;  v6 = v5*v;  v7 = v6*v;

        /* compute the phase and amplitude */
        if ((f < fLower) || (f > fFinal)) {
            amp = 0.;
            Psi = 0.;
        }
        else {

            Psi = psi0*pow(v, -5.)*(1. 
                    + psi2*v2 + psi3*v3 + psi4*v4 
                    + psi5*v5*(1.+3.*log(v)) 
                    + (psi6 + psi6L*log(4.*v))*v6 + psi7*v7); 

            amp = amp0*pow(f, mSevenBySix)*(1. 
                    + alpha2*v2 + alpha3*v3 + alpha4*v4 
                    + alpha5*v5 + (alpha6 + alpha6L*(LAL_GAMMA+log(4.*v)) )*v6 
                    + alpha7*v7); 

        }

        /* generate the waveform */
       	hOfF[k] = amp * (cos(Psi+shft*f+phi0+piBy4) - I*sin(Psi+shft*f+phi0+piBy4)); 

    }    

	return 0;
}

/* FIXME make this function exist in LAL and have the LAL SPA waveform generator call it? */
static int SPAWaveform (double mass1, double mass2, int order, double deltaF, double deltaT, double fLower, double fFinal, int numPoints,  complex double *expPsi)
	{
	double m = mass1 + mass2;
	double eta = mass1 * mass2 / m / m;
	double mchirp = m * pow(eta, 3.0 / 5.0);

	double x1 = pow (LAL_PI * m * LAL_MTSUN_SI * deltaF, -1.0 / 3.0);
	double psi = 0.0;
	double psi0 = 0.0;
	int k = 0;
	int kmin = fLower / deltaF > 1 ? fLower / deltaF : 1;
	int kmax = fFinal / deltaF < numPoints / 2 ? fFinal / deltaF : numPoints / 2;

	const double cannonDist = 1.0; /* Mpc */
	double distNorm = LAL_MRSUN_SI / (cannonDist * 1.0e6 * LAL_PC_SI);
	/* from FINDCHIRP paper arXiv:gr-qc/0509116v2 */
	double tNorm = sqrt(5.0 / 24.0 / LAL_PI) * pow(LAL_PI * LAL_MTSUN_SI, -1.0 / 6.0) * pow(mchirp, 5.0 / 6.0);

	/* pn constants */
	double c0, c10, c15, c20, c25, c25Log, c30, c30Log, c35, c40P;
	double x;

	/* chebychev coefficents for expansion of sin and cos */
	const double s2 = -0.16605;
	const double s4 = 0.00761;
	const double c2 = -0.49670;
	const double c4 = 0.03705;

	complex double value;

	/* template norm */
	tNorm *= distNorm;

	/* zero output */
	memset (expPsi, 0, numPoints * sizeof (complex double));

	/* Initialize all PN phase coeffs to zero. */
	c0 = c10 = c15 = c20 = c25 = c25Log = c30 = c30Log = c35 = c40P = 0.;

	/* Switch on PN order, set the appropriate phase coeffs for that order */
	switch (order)
		{
		case 8:
			c40P = 3923.0;
		case 7:
			c35 = LAL_PI * (77096675.0 / 254016.0 + eta * 378515.0 / 1512.0 - eta * eta * 74045.0 / 756.0);
		case 6: 
			c30 = 11583231236531.0 / 4694215680.0 - LAL_GAMMA * 6848.0 / 21.0 - LAL_PI * LAL_PI * 640.0 / 3.0 + eta * (LAL_PI * LAL_PI * 2255.0 / 12.0 - 15737765635.0 / 3048192.0) + eta * eta * 76055.0 / 1728.0 - eta * eta * eta * 127825.0 / 1296.0 - 6848.0 * log (4.0) / 21.0;
			c30Log = -6848.0 / 21.0;
		case 5:
			c25 = LAL_PI * 38645.0 / 756.0 - LAL_PI * eta * 65.0 / 9.0;
			c25Log = 3 * c25;
		case 4:
			c20 = 15293365.0 / 508032.0 + eta * (27145.0 / 504.0 + eta * 3085.0 / 72.0);
			c15 = -16 * LAL_PI;
			c10 = 3715.0 / 756.0 + eta * 55.0 / 9.0;
			c0 = 3.0 / (eta * 128.0);
			break;
		default:
			fprintf (stderr, "unknown PN order (use 8 for 4PN 7 for 3.5 PN ... 4 for 2PN ) nothing above 4PN or below 2PN is recognized");
			break;
		}

	/* x1 */
	x1 = pow (LAL_PI * m * LAL_MTSUN_SI * deltaF, -1.0 / 3.0);
	x = x1 * pow ((double) kmin, -1.0 / 3.0);

	psi = c0 * (x * (c20 + x * (c15 + x * (c10 + x * x))) + c25 - c25Log * log (x) + (1.0 / x) * (c30 - c30Log * log (x) + (1.0 / x) * (c35 - (1.0 / x) * c40P * log (x))));
	psi0 = -2 * LAL_PI * (floor (0.5 * psi / LAL_PI));
	
	/* Chirp Time */
	/* This formula works for any PN order, because */
	/* higher order coeffs will be set to zero.     */
	for (k = kmin; k < kmax; ++k)
		{
		double x = x1 * pow ((double) k, -1.0 / 3.0);
		double psi = c0 * (x * (c20 + x * (c15 + x * (c10 + x * x))) + c25 - c25Log * log (x) + (1.0 / x) * (c30 - c30Log * log (x) + (1.0 / x) * (c35 - (1.0 / x) * c40P * log (x))));

		double psi1 = psi + psi0;
		double psi2;

		/* range reduction of psi1 */
		while (psi1 < -LAL_PI)
			{
			psi1 += 2 * LAL_PI;
			psi0 += 2 * LAL_PI;
			}
		while (psi1 > LAL_PI)
			{
			psi1 -= 2 * LAL_PI;
			psi0 -= 2 * LAL_PI;
			}

		/* compute approximate sine and cosine of psi1 */
		if (psi1 < -LAL_PI / 2)
			{
			psi1 = -LAL_PI - psi1;
			psi2 = psi1 * psi1;
			/* XXX minus sign added because of new sign convention for fft */
			/* FIXME minus sign put back because it makes a reverse chirp with scipy's ifft */
			value = psi1 * (1 + psi2 * (s2 + psi2 * s4)) +  I * (0. - 1. - psi2 * (c2 + psi2 * c4));
			expPsi[k]  = value;
			}
		else if (psi1 > LAL_PI / 2)
			{
			psi1 = LAL_PI - psi1;
			psi2 = psi1 * psi1;
			/* XXX minus sign added because of new sign convention for fft */
			/* FIXME minus sign put back because it makes a reverse chirp with scipy's ifft */
			value = psi1 * (1 + psi2 * (s2 + psi2 * s4)) + I * (0. - 1. - psi2 * (c2 + psi2 * c4));
			expPsi[k] = value;
			}
		else
			{
			psi2 = psi1 * psi1;
			/* XXX minus sign added because of new sign convention for fft */
			/* FIXME minus sign put back because it makes a reverse chirp with scipy's ifft */
			value = psi1 * (1 + psi2 * (s2 + psi2 * s4)) + I * (1. + psi2 * (c2 + psi2 * c4));
			expPsi[k] = value;
			}
		/* put in the first order amplitude factor */
		expPsi[k] *= pow(k*deltaF, -7.0 / 6.0) * tNorm;
		}
	return 0;
	}

/* FIXME make this function exist in LAL and have the LAL SPA waveform generator call it? */
static int GenericSPAWaveform (double *psis, double *psils, int order, double deltaF, double deltaT, double fLower, double fFinal, int numPoints,  complex double *expPsi)
	{
	int i,k = 0;
	int kmin = fLower / deltaF > 1 ? fLower / deltaF : 1;
	int kmax = fFinal / deltaF < numPoints / 2 ? fFinal / deltaF : numPoints / 2;

	double x1, x, logx, psi, dpsi, psi0, psi1, psi2;

	/* chebychev coefficents for expansion of sin and cos */
	const double s2 = -0.16605;
	const double s4 = 0.00761;
	const double c2 = -0.49670;
	const double c4 = 0.03705;

	complex double value;

	/* zero output */
	memset (expPsi, 0, numPoints * sizeof (complex double));

	/* f1 */
	x1 = pow ((double) deltaF, 1.0 / 3.0);
	x = x1 * pow ((double) kmin, 1.0 / 3.0);
	logx = log(x);

	psi = 0;
	for (i = order; i >= 0; i--)
		{
		dpsi = psis[i] + psils[i]*3.*logx;
		psi = dpsi + x*psi;
		}
	psi *= pow(x, -5.);
	psi0 = -2 * LAL_PI * (floor (0.5 * psi / LAL_PI));
	
	for (k = kmin; k < kmax; ++k)
		{
		x = x1 * pow ((double) k, 1.0 / 3.0);
		logx = log(x);

		psi = 0;
		for (i = order; i >= 0; i--)
			{
			dpsi = psis[i] + psils[i]*3.*logx;
			psi = dpsi + x*psi;
			}
		psi *= pow(x, -5);

		psi1 = psi + psi0;

		/* range reduction of psi1 */
		while (psi1 < -LAL_PI)
			{
			psi1 += 2 * LAL_PI;
			psi0 += 2 * LAL_PI;
			}
		while (psi1 > LAL_PI)
			{
			psi1 -= 2 * LAL_PI;
			psi0 -= 2 * LAL_PI;
			}

		/* compute approximate sine and cosine of psi1 */
		if (psi1 < -LAL_PI / 2)
			{
			psi1 = -LAL_PI - psi1;
			psi2 = psi1 * psi1;
			/* XXX minus sign added because of new sign convention for fft */
			/* FIXME minus sign put back because it makes a reverse chirp with scipy's ifft */
			value = psi1 * (1 + psi2 * (s2 + psi2 * s4)) +  I * (0. - 1. - psi2 * (c2 + psi2 * c4));
			expPsi[k]  = value;
			}
		else if (psi1 > LAL_PI / 2)
			{
			psi1 = LAL_PI - psi1;
			psi2 = psi1 * psi1;
			/* XXX minus sign added because of new sign convention for fft */
			/* FIXME minus sign put back because it makes a reverse chirp with scipy's ifft */
			value = psi1 * (1 + psi2 * (s2 + psi2 * s4)) + I * (0. - 1. - psi2 * (c2 + psi2 * c4));
			expPsi[k] = value;
			}
		else
			{
			psi2 = psi1 * psi1;
			/* XXX minus sign added because of new sign convention for fft */
			/* FIXME minus sign put back because it makes a reverse chirp with scipy's ifft */
			value = psi1 * (1 + psi2 * (s2 + psi2 * s4)) + I * (1. + psi2 * (c2 + psi2 * c4));
			expPsi[k] = value;
			}
		/* put in the first order amplitude factor */
		expPsi[k] *= pow(k*deltaF, -7.0 / 6.0);
		}
	return 0;
	}

/* A function to compute the time of a spa waveform starting at a given low
 * frequency and going to infinity
 */

static double chirp_time (double m1, double m2, double fLower, int order, double chi)
	{

	/* variables used to compute chirp time */
	double c0T, c2T, c3T, c4T, c5T, c6T, c6LogT, c7T;
	double xT, x2T, x3T, x4T, x5T, x6T, x7T, x8T;
	double m = m1 + m2;
	double eta = m1 * m2 / m / m;

	c0T = c2T = c3T = c4T = c5T = c6T = c6LogT = c7T = 0.;

	/* Switch on PN order, set the chirp time coeffs for that order */
	switch (order)
		{
		case 8:
		case 7:
			c7T = LAL_PI * (14809.0 * eta * eta - 75703.0 * eta / 756.0 - 15419335.0 / 127008.0);
		case 6:
			c6T = LAL_GAMMA * 6848.0 / 105.0 - 10052469856691.0 / 23471078400.0 + LAL_PI * LAL_PI * 128.0 / 3.0 + eta * (3147553127.0 / 3048192.0 - LAL_PI * LAL_PI * 451.0 / 12.0) - eta * eta * 15211.0 / 1728.0 + eta * eta * eta * 25565.0 / 1296.0 + log (4.0) * 6848.0 / 105.0;
     			c6LogT = 6848.0 / 105.0;
		case 5:
			c5T = 13.0 * LAL_PI * eta / 3.0 - 7729.0 / 252.0 - (0.4*565.*(-146597. + 135856.*eta + 17136.*eta*eta)*chi/(2268.*(-113. + 76.*eta))); // last term is 0 if chi is 0;
		case 4:
			c4T = 3058673.0 / 508032.0 + eta * (5429.0 / 504.0 + eta * 617.0 / 72.0) + (0.4*63845.*(-81. + 4.*eta)*chi*chi/(8.*pow(-113. + 76.*eta, 2.))); // last term is 0 if chi is 0;
			c3T = -32.0 * LAL_PI / 5.0 + (0.4*113.*chi/3.); // last term is 0 if chi is 0;
			c2T = 743.0 / 252.0 + eta * 11.0 / 3.0;
			c0T = 5.0 * m * LAL_MTSUN_SI / (256.0 * eta);	
			break;
		default:
			fprintf (stderr, "ERROR!!!\n");
			break;
		}

	/* This is the PN parameter v evaluated at the lower freq. cutoff */
	xT = pow (LAL_PI * m * LAL_MTSUN_SI * fLower, 1.0 / 3.0);
	x2T = xT * xT;
	x3T = xT * x2T;
	x4T = x2T * x2T;
	x5T = x2T * x3T;
	x6T = x3T * x3T;
	x7T = x3T * x4T;
	x8T = x4T * x4T;

	/* Computes the chirp time as tC = t(v_low)    */
	/* tC = t(v_low) - t(v_upper) would be more    */
	/* correct, but the difference is negligble.   */

	/* This formula works for any PN order, because */
	/* higher order coeffs will be set to zero.     */

	return c0T * (1 + c2T * x2T + c3T * x3T + c4T * x4T + c5T * x5T + (c6T + c6LogT * log (xT)) * x6T + c7T * x7T) / x8T;
}

/* A convenience function to compute the time between two frequencies */
static double chirp_time_between_f1_and_f2(double m1, double m2, double fLower, double fUpper, int order, double chi)
	{
	return chirp_time(m1,m2,fLower,order,chi) - chirp_time(m1,m2,fUpper,order,chi);
	}


static double schwarz_isco(double m1, double m2)
	{
	double m = LAL_MTSUN_SI * (m1 + m2);
	return 1.0 / pow(6.0, 1.5) / m / LAL_PI; 
	}

static double light_ring(double m1, double m2)
	{
	double m = LAL_MTSUN_SI * (m1 + m2);
	return 1.0 / pow(3.0, 1.5) / m / LAL_PI; 
	}

static double bkl_isco(double m1, double m2)
	{
	double q;
	q = (m1 < m2) ? (m1/m2) : (m2/m1);
	return (0.8 * q * q * q - 2.6 * q * q + 2.8 * q + 1.0 ) * schwarz_isco(m1,m2);
	}
	
static double compute_chi(double m1, double m2, double spin1, double spin2) {
	double totalMass = m1+m2;
	double eta = m1*m2/pow(totalMass,2.);
	double delta = sqrt(1.-4.*eta);
	/* spin parameter used for the search */
	double chi = 0.5*(spin1*(1.+delta) + spin2*(1.-delta));
	return chi;
	}

static double imr_merger(double m1, double m2, double chi) {
	double totalMass = m1+m2;
	double eta = m1*m2/pow(totalMass,2.);
	double piM = totalMass*LAL_PI*LAL_MTSUN_SI;
	double fMerg =  1. - 4.4547*pow(1.-chi,0.217) + 3.521*pow(1.-chi,0.26) +
		6.4365e-01*eta + 8.2696e-01*eta*chi + -2.7063e-01*eta*pow(chi,2.) +
		-5.8218e-02*pow(eta,2.) + -3.9346e+00*pow(eta,2.)*chi +
		-7.0916e+00*pow(eta,3.);
    
	fMerg /= piM;
	return fMerg;
	}

static double imr_ring(double m1, double m2, double chi){
	double totalMass = m1+m2;
	double eta = m1*m2/pow(totalMass,2.);
	double piM = totalMass*LAL_PI*LAL_MTSUN_SI;
	double fRing = (1. - 0.63*pow(1.-chi,0.3))/2. +
		1.4690e-01*eta + -1.2281e-01*eta*chi + -2.6091e-02*eta*pow(chi,2.) +
		-2.4900e-02*pow(eta,2.) + 1.7013e-01*pow(eta,2.)*chi +
		2.3252e+00*pow(eta,3.);
	fRing /= piM;
	return fRing;
	}

static double imr_fcut(double m1, double m2, double chi) {
	double totalMass = m1+m2;
	double eta = m1*m2/pow(totalMass,2.);
	double piM = totalMass*LAL_PI*LAL_MTSUN_SI;
	double fCut = 3.2361e-01 + 4.8935e-02*chi + 1.3463e-02*pow(chi,2.) +
		-1.3313e-01*eta + -8.1719e-02*eta*chi + 1.4512e-01*eta*pow(chi,2.) +
		-2.7140e-01*pow(eta,2.) + 1.2788e-01*pow(eta,2.)*chi +
		4.9220e+00*pow(eta,3.);
	fCut /= piM;
	return fCut;
	}

/* FIXME a hack fit to unwrap the imr waveforms and provide a time correction*/
static double _imrdur(double m1, double m2, double chi)
	{
	/* FIT done to + 1500 M */
	double c1 = 1.39 - 0.88 * chi;
	double c2 = -4.6 - 1.2 * pow(chi + 0.7, 3);
	double c3 = 400.0 + 500.0 * pow(chi + 1.0, 3);
	double eta = m1*m2 / (m1+m2) / (m1+m2);
	double out = c1 * exp(c2 * eta + c3 * pow(eta,6)) * ((m1+m2) / 50.0 ) - (1500 - fabs(chi) * 200) * (m1+m2) * LAL_MTSUN_SI;
	return out;
	}

int IMRSPAWaveform(double mass1, double mass2, double spin1, double spin2, double deltaF, double fLower, int numPoints,  complex double *hOfF) {
	double chi = compute_chi(mass1, mass2, spin1, spin2);
	return IMRSPAWaveformFromChi(mass1, mass2, chi, deltaF, fLower, numPoints, hOfF);
	}

/****************************************************************************/
/* Ajith's code *************************************************************/
/* FIXME make white space style similar	*************************************/
/****************************************************************************/
int IMRSPAWaveformFromChi(double mass1, double mass2, double chi, double deltaF, double fLower, int numPoints, complex double *hOfF) {

    double totalMass, piM, eta;
    double psi0, psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8, fMerg, fRing, fCut, sigma;
    double f, shft, amp0, ampEff, psiEff, fNorm;
    double v, alpha2, alpha3, w1, vMerg, epsilon_1, epsilon_2, w2, vRing;
    double startPhase = 0., startTime, distance, Lorentzian;
    int k, kmin, kmax;

    /* calculate the total mass, symmetric mass ratio and asymmetric       */
    /* mass ratio                                                          */
    totalMass = mass1+mass2;
    eta = mass1*mass2/pow(totalMass,2.);
    piM = totalMass*LAL_PI*LAL_MTSUN_SI;
    
    /***********************************************************************/
    /*              compute the phenomenological parameters                */
    /***********************************************************************/

    psi0 = 3./(128.*eta);

    psi2 = 3715./756. +
            -9.2091e+02*eta + 4.9213e+02*eta*chi + 1.3503e+02*eta*pow(chi,2.) +
            6.7419e+03*pow(eta,2.) + -1.0534e+03*pow(eta,2.)*chi +
            -1.3397e+04*pow(eta,3.);

    psi3 = -16.*LAL_PI + 113.*chi/3. +
            1.7022e+04*eta + -9.5659e+03*eta*chi + -2.1821e+03*eta*pow(chi,2.) +
            -1.2137e+05*pow(eta,2.) + 2.0752e+04*pow(eta,2.)*chi +
            2.3859e+05*pow(eta,3.);

    psi4 = 15293365./508032. - 405.*pow(chi,2.)/8. +
            -1.2544e+05*eta + 7.5066e+04*eta*chi + 1.3382e+04*eta*pow(chi,2.) +
            8.7354e+05*pow(eta,2.) + -1.6573e+05*pow(eta,2.)*chi +
            -1.6936e+06*pow(eta,3.);

    psi6 = -8.8977e+05*eta + 6.3102e+05*eta*chi + 5.0676e+04*eta*pow(chi,2.) +
            5.9808e+06*pow(eta,2.) + -1.4148e+06*pow(eta,2.)*chi +
            -1.1280e+07*pow(eta,3.);

    psi7 = 8.6960e+05*eta + -6.7098e+05*eta*chi + -3.0082e+04*eta*pow(chi,2.) +
            -5.8379e+06*pow(eta,2.) + 1.5145e+06*pow(eta,2.)*chi +
            1.0891e+07*pow(eta,3.);

    psi8 = -3.6600e+05*eta + 3.0670e+05*eta*chi + 6.3176e+02*eta*pow(chi,2.) +
            2.4265e+06*pow(eta,2.) + -7.2180e+05*pow(eta,2.)*chi + 
            -4.5524e+06*pow(eta,3.);

    fMerg =  1. - 4.4547*pow(1.-chi,0.217) + 3.521*pow(1.-chi,0.26) +
            6.4365e-01*eta + 8.2696e-01*eta*chi + -2.7063e-01*eta*pow(chi,2.) +
            -5.8218e-02*pow(eta,2.) + -3.9346e+00*pow(eta,2.)*chi +
            -7.0916e+00*pow(eta,3.);

    fRing = (1. - 0.63*pow(1.-chi,0.3))/2. +
            1.4690e-01*eta + -1.2281e-01*eta*chi + -2.6091e-02*eta*pow(chi,2.) +
            -2.4900e-02*pow(eta,2.) + 1.7013e-01*pow(eta,2.)*chi +
            2.3252e+00*pow(eta,3.);

    sigma = (1. - 0.63*pow(1.-chi,0.3))*pow(1.-chi,0.45)/4. +
            -4.0979e-01*eta + -3.5226e-02*eta*chi + 1.0082e-01*eta*pow(chi,2.) +
            1.8286e+00*pow(eta,2.) + -2.0169e-02*pow(eta,2.)*chi +
            -2.8698e+00*pow(eta,3.);

    fCut = 3.2361e-01 + 4.8935e-02*chi + 1.3463e-02*pow(chi,2.) +
            -1.3313e-01*eta + -8.1719e-02*eta*chi + 1.4512e-01*eta*pow(chi,2.) +
            -2.7140e-01*pow(eta,2.) + 1.2788e-01*pow(eta,2.)*chi +
            4.9220e+00*pow(eta,3.);

    psi1    = 0.;
    psi5    = 0.;

    /* appropriately scale with mass */
    fCut  /= piM;
    fMerg /= piM;
    fRing /= piM;
    sigma /= piM;

	kmin = fLower / deltaF > 1 ? fLower / deltaF : 1;
	kmax = fCut / deltaF < numPoints / 2 ? fCut / deltaF : numPoints / 2;

    /*********************************************************************/
    /*      set other parameters required for the waveform generation    */
    /*********************************************************************/
    startTime = -50.*totalMass*LAL_MTSUN_SI;
    shft = 2.*LAL_PI *startTime;
    distance = 1e6*LAL_PC_SI;  /*1 Mpc in meters */
        
    /* Now compute the amplitude.  NOTE the distance is assumed to          */
    /* be in meters. This is, in principle, inconsistent with the LAL       */
    /* documentation (inspiral package). But this seems to be the           */
    /* convention employed in the injection codes                           */
    
    amp0 = pow(LAL_MTSUN_SI*totalMass, 5./6.)*pow(fMerg,-7./6.)/pow(LAL_PI,2./3.);
    amp0 *= pow(5.*eta/24., 1./2.)/(distance/LAL_C_SI);

    /* PN correctiosn to the frequency domain amplitude of the (2,2) mode   */
    alpha2   = -323./224. + 451.*eta/168.;
    alpha3   = (27./8. - 11.*eta/6.)*chi;

    /* spin-dependant corrections to the merger amplitude                   */
    epsilon_1 =  1.4547*chi - 1.8897;
    epsilon_2 = -1.8153*chi + 1.6557;

    /* normalisation constant of the inspiral amplitude                     */
    vMerg = pow(LAL_PI*totalMass*LAL_MTSUN_SI*fMerg, 1./3.);
    vRing = pow(LAL_PI*totalMass*LAL_MTSUN_SI*fRing, 1./3.);

    w1 = 1. + alpha2*pow(vMerg,2.) + alpha3*pow(vMerg,3.);
    w1 = w1/(1. + epsilon_1*vMerg + epsilon_2*vMerg*vMerg);
    w2 = w1*(LAL_PI*sigma/2.)*pow(fRing/fMerg, -2./3.)*(1. + epsilon_1*vRing
            + epsilon_2*vRing*vRing);

    /* zero output */
    memset (hOfF, 0, numPoints * sizeof (complex double));
    ampEff = 0.;
    psiEff = 0.;

    /************************************************************************/
    /*          now generate the waveform at all frequency bins             */
    /************************************************************************/
    for (k = kmin; k < kmax; k++) {

        /* fourier frequency corresponding to this bin                      */
      	f = k * deltaF;
        fNorm = f/fMerg;

        /* PN expansion parameter                                           */
        v = pow(LAL_PI*totalMass*LAL_MTSUN_SI*f, 1./3.);

    	/* compute the amplitude                                            */
        if (f <= fMerg) {
            ampEff = pow(fNorm, -7./6.)*(1. + alpha2*pow(v,2.) + alpha3*pow(v,3.));
        }
        else if ((f > fMerg) & (f <= fRing)) {
            ampEff = w1*pow(fNorm, -2./3.)*(1. + epsilon_1*v + epsilon_2*v*v);
        }
        else if (f > fRing) {
            Lorentzian =  sigma / (2.*LAL_PI * (pow(f-fRing, 2.) + sigma*sigma/4.0));
            ampEff = w2*Lorentzian;
        }

        /* now compute the phase                                             */
        psiEff =  shft*f + startPhase
                    + 3./(128.*eta*pow(v,5.))*(1 + psi2*pow(v, 2.)
                    + psi3*pow(v, 3.) + psi4*pow(v, 4.)
                    + psi5*pow(v, 5.) + psi6*pow(v, 6.)
                    + psi7*pow(v, 7.) + psi8*pow(v, 8.));

        /* generate the waveform                                             */
        hOfF[k] = amp0*ampEff * (cos(psiEff) - I * sin(psiEff)); 

    }

    return 0;

}

