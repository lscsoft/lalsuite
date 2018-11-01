/*
*  Copyright (C) 2007 Xavier Siemens
*  Copyright (C) 2010 Andrew Mergl
*  Copyright (C) 2018 Daichi Tsuna
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*********************************************************************************/
/*         Cosmic string burst rate computation code for small/large loops       */
/*                                                                               */
/*                  Xavier Siemens, Jolien Creighton, Irit Maor                  */
/*                                                                               */
/*                         UWM/Caltech - September 2006                          */
/*********************************************************************************/
/*Modified June 2010 by Andrew Mergl for use with Python*/
/* Modified November 2018 by Daichi Tsuna to include large loop scenarios */
/* Updated to be able to test 3 models of large loop distributions */

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <lal/cs_cosmo.h>
#include <lal/cs_lambda_cosmo.h>
#include "six.h"

#define CUSPS_PER_LOOP 1.0		/* c */
#define LOOP_RAD_POWER 50.0		/* Gamma */

#define H0 LAMBDA_H_0


/*****************************************************************************/
/*int finddRdz(double Gmu, double alpha, double f, double Gamma, int Namp, double *zofA, double *dRdz)
 * Find dR/dz given Gmu, alpha, f, and Gamma. This is a C function that was taken from cs_gamma.c
 * and modified so that it could be called from Python. The global variables
 * that the function used to use are now all passed through arguments.
 *
 * Arguments:
 * Gmu, alpha, gamma: parameters calculated by the main program. See technical
 *  document for what they represent cosmologically.
 * f: the frequency that is passed to the main program with --frequency opt
 * Namp: The size of the data arrays. This is set when opening the data file
 * *zofA, *dRdz: 1D arrays of length Namp. See technical document for what
 *  they represent cosmologically.
 */
/*****************************************************************************/
static PyObject *cs_gamma_finddRdz(PyObject *self, PyObject *args)
{
  PyArrayObject *Numpy_zofA;
  PyObject *Numpy_dRdz;
  double Gmu, alpha, f, Gamma, *zofA, *dRdz;
  long int Namp;
  cs_cosmo_functions_t cosmofns;
  int j;
  (void)self;	/* silence unused parameter warning */

  if (!PyArg_ParseTuple(args, "ddddO!", &Gmu, &alpha, &f, &Gamma, &PyArray_Type, &Numpy_zofA))
    return NULL;

  Numpy_zofA = PyArray_GETCONTIGUOUS(Numpy_zofA);
  if(!Numpy_zofA)
    return NULL;
  Namp = PyArray_DIM(Numpy_zofA, 0);
  zofA = PyArray_DATA(Numpy_zofA);

  {
  npy_intp dims[1] = {Namp};
  Numpy_dRdz = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  }
  dRdz = PyArray_DATA((PyArrayObject *) Numpy_dRdz);

  cosmofns = XLALCSCosmoFunctions( zofA, Namp);

  for ( j = 0; j < Namp; j++ )
    {
      /*double theta = pow((1+cosmofns.z[j]) * f * alpha * cosmofns.phit[j] / H0, -1.0/3.0);

      if (theta > 1.0)
          dRdz[j] = 0.0;
      else*/
          dRdz[j] = 0.5 * H0 * pow(f/H0,-2.0/3.0) * pow(alpha, -5.0/3.0) / (Gamma*Gmu) * pow(cosmofns.phit[j],-14.0/3.0) * cosmofns.phiV[j] * pow(1+cosmofns.z[j],-5.0/3.0);
      if(gsl_isnan(dRdz[j])) {
        Py_DECREF(Numpy_dRdz);
        Numpy_dRdz = NULL;
        break;
      }
    }

  XLALCSCosmoFunctionsFree( cosmofns );
  Py_DECREF(Numpy_zofA);

  return Numpy_dRdz;
}
/*****************************************************************************/
/*int findzofA(double Gmu, double alpha, int Namp, double *zofA, double *amp)
 *Find z(A) given Gmu and alpha. This function was taken from cs_gamma.c and
 * modified so that it could be called from Python. The global variables it
 * used to use are now passed to it.
 *
 * Arugments:
 * Gmu, alpha: values calculated by the main program. See S4 technical
 *  documentation for their cosmological meanings.
 * Namp: The length of the data arrays. This is set when reading in the data
 * zofA, amp: 1D data arrays of length Namp. See S4 technical documentation
 *  for their cosmological meanings.
 */
/*****************************************************************************/
static PyObject *cs_gamma_findzofA(PyObject *self, PyObject *args)
{
  PyArrayObject *Numpy_amp;
  PyObject *Numpy_zofA;
  double Gmu, alpha, *zofA, *amp;
  unsigned long int Namp;
  (void)self;	/* silence unused parameter warning */

  double z_min = 1e-20, z_max = 1e10;
  double dlnz = 0.05;
  unsigned numz = floor( (log(z_max) - log(z_min)) / dlnz );
  unsigned long int i;
  cs_cosmo_functions_t cosmofns;
  double *fz,*z;
  double a;
  gsl_interp *zofa_interp;
  gsl_interp_accel *acc_zofa = gsl_interp_accel_alloc();

  if (!PyArg_ParseTuple(args, "ddO!", &Gmu, &alpha, &PyArray_Type, &Numpy_amp))
    return NULL;

  Numpy_amp = PyArray_GETCONTIGUOUS(Numpy_amp);
  if(!Numpy_amp)
    return NULL;
  Namp = PyArray_DIM(Numpy_amp, 0);
  amp = PyArray_DATA(Numpy_amp);

  {
  npy_intp dims[1] = {Namp};
  Numpy_zofA = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  }
  zofA = PyArray_DATA((PyArrayObject *) Numpy_zofA);

  cosmofns = XLALCSCosmoFunctionsAlloc( z_min, dlnz, numz );

  zofa_interp = gsl_interp_alloc (gsl_interp_linear, cosmofns.n);

  fz = calloc( cosmofns.n, sizeof( *fz ) );
  z = calloc( cosmofns.n, sizeof( *z ) );

  /* first compute the function that relates A and z */
  /* invert order; b/c fz is a monotonically decreasing func of z */
  for ( i = cosmofns.n ; i > 0; i-- )
    {
      unsigned long int j = cosmofns.n - i;
      z[j] = cosmofns.z[i-1];
      fz[j] = pow(cosmofns.phit[i-1], 2.0/3.0) * pow(1+z[j], -1.0/3.0) / cosmofns.phiA[i-1];
    }

  gsl_interp_init (zofa_interp, fz, z, cosmofns.n);

  /* now compute the amplitudes (suitably multiplied) that are equal to fz for some z*/
  for ( i = 0; i < Namp; i++ )
    {
      a = amp[i] * pow(H0,-1.0/3.0) * pow(alpha,-2.0/3.0) / Gmu;
      /* evaluate z(fz) at fz=a */
      zofA[i] = gsl_interp_eval (zofa_interp, fz, z, a, acc_zofa );
      if(gsl_isnan(zofA[i])) {
        Py_DECREF(Numpy_zofA);
        Numpy_zofA = NULL;
        break;
      }
    }

  XLALCSCosmoFunctionsFree( cosmofns );
  Py_DECREF(Numpy_amp);
  free(fz);
  free(z);
  gsl_interp_free (zofa_interp);
  gsl_interp_accel_free(acc_zofa);

  return Numpy_zofA;
}


/*******************************************************************************/
/* nu(double Gmu, double Gamma, double *amp, double *z, char model)
 * Find loop distribution givem Gmu, Gamma, amplitude, z, and a particular large loop model.
 * The loop distribution nu is calculated as a function of amp and z, and is
 * returned as an array having a size of Namp * Nz.
 */
/*****************************************************************************/
static double nu(double Gmu, double Gamma, double A, double z, double phit, double phiA, const char *model)
{
	double l = pow( A / Gmu / H0 * pow(1.0+z,1.0/3.0)* phiA, 3.0/2.0);
	double alpha, nuR, nuM;
	double P_R=1.60, P_M=1.41, a_index_R=1.0/2.0, a_index_M=2.0/3.0;
	double Upsilon = 10.0;
	double chi_R=1.0-P_R/2.0;
	double chi_M=1.0-P_M/2.0;
	double gamma_c_R, gamma_c_M;
	double t=phit/H0;
	double cuspdist;
	double teq=8.122570474611143e+11;
	double crateR, crateRadStragglers, crateM;
	double crate;

	if (strcmp(model,"Siemens06") == 0) {
		alpha=1e-1;
		nuR=0.4*15*sqrt(alpha);
		nuM=0.12*4;
		/* Radiation era loops */
		crateR = 0.0;
		if( (l < alpha*t) && (t < teq) )
			crateR = nuR * pow(t,-1.5) * pow( l + Gamma*Gmu/H0 * phit, -2.5 );

		/* Radiation stragglers */
		crateRadStragglers = 0.0;
		if ( (l <  alpha*teq-Gamma*Gmu*(t-teq) ) && ( t > teq ) )
			crateRadStragglers = nuR*pow(teq, 0.5)/t/t* pow( l + Gamma*Gmu/H0 * phit, -2.5);

		/* matter era loops */
		crateM = 0.0;
		if( (l < alpha*t) && (t > teq) && (l >  alpha*teq-Gamma*Gmu*(t-teq)) )
			crateM = nuM / t / t / ( l + Gamma*Gmu/H0 * phit) / ( l + Gamma*Gmu/H0 * phit);

		crate = crateR + crateRadStragglers + crateM;
		} else if (strcmp(model,"Blanco-Pillado14") == 0) {
		alpha = 0.1;
		nuR = 0.18;
		/* Radiation era loops */
		crateR = 0.0;
		if( (l < alpha*t) && (t < teq) )
			crateR = nuR * pow(t,-1.5) * pow( l + Gamma*Gmu/H0 * phit, -2.5 );

		/* Radiation stragglers */
		crateRadStragglers = 0.0;
		if ( (l <  alpha*teq-Gamma*Gmu*(t-teq) ) && ( t > teq ) )
			crateRadStragglers = nuR*pow(teq, 0.5)/t/t* pow( l + Gamma*Gmu/H0 * phit, -2.5);

		/* matter era loops */
		crateM = 0.0;
		if( (l < 0.18*t) && (t > teq) && (l >  alpha*teq-Gamma*Gmu*(t-teq)) )
			crateM = (0.27-0.45*pow(l/t,0.31)) / t / t / ( l + Gamma*Gmu/H0 * phit) / ( l + Gamma*Gmu/H0 * phit);

		crate = crateR + crateRadStragglers + crateM;
	} else if (strcmp(model,"Ringeval07") == 0) {
		nuR = 0.21 * pow(1.0-a_index_R, 3.0-P_R);
		nuM = 0.09 * pow(1.0-a_index_M, 3.0-P_M);
		crateR = 0.0;
		crateM = 0.0;
		gamma_c_R = Upsilon * pow(Gmu,1.0+2.0*chi_R);
		gamma_c_M = Upsilon * pow(Gmu,1.0+2.0*chi_M);
		/* Radiation era loops */
		if( (l > Gamma * Gmu * t) && (l <= t/(1.0-a_index_R)) )
			crateR = pow(t,-4.0) * nuR / pow(l/t + Gamma*Gmu, P_R+1);

		if( (l > gamma_c_R * t) && (l <= Gamma * Gmu * t) )
			crateR = pow(t,-4.0) * nuR * (3.0*a_index_R-2.0*chi_R-1.0) / (2.0-2.0*chi_R) / (Gamma*Gmu) / pow(l/t,P_R);

		if(l <= gamma_c_R * t)
			crateR = pow(t,-4.0) * nuR * (3.0*a_index_R-2.0*chi_R-1.0) / (2.0-2.0*chi_R) / (Gamma*Gmu) / pow(gamma_c_R,P_R);
		/* Matter era loops */
		if( (l > Gamma * Gmu * t) && (l <= t/(1.0-a_index_M)) )
			crateM = pow(t,-4.0) * nuM / pow(l/t + Gamma*Gmu, P_M+1);

		if( (l > gamma_c_M * t) && (l <= Gamma * Gmu * t) )
			crateM = pow(t,-4.0) * nuM * (3.0*a_index_M-2.0*chi_M-1.0) / (2.0-2.0*chi_M) / (Gamma*Gmu) / pow(l/t,P_M);

		if(l <= gamma_c_M * t)
			crateM = pow(t,-4.0) * nuM * (3.0*a_index_M-2.0*chi_M-1.0) / (2.0-2.0*chi_M) / (Gamma*Gmu) / pow(gamma_c_M,P_M);

		crate = crateR + crateM;
	}

	cuspdist = 3.0/A * crate;

	return cuspdist;
}

/*******************************************************************************/
/* finddRdzdA(double Gmu, double f, double Gamma, double *amp, double *z, double *nu)
 * Find d^2R/dzdA givem Gmu, f, Gamma, amplitude, z, and the loop distribution.
 * Using the equations formulated in e.g. arXiv:1712.01168 Appendix B, d^2R/dzdA is calculated
 * as a function of amp and z, and is returned as an array having a size of Namp * Nz.
 */
/*****************************************************************************/
static PyObject *cs_gamma_finddRdzdA(PyObject *self, PyObject *args)
{
	PyArrayObject *Numpy_amp, *Numpy_z;
	PyObject *Numpy_dRdzdA;
	double Gmu, f, Gamma, *amp, *z;
	long int Namp, Nz;
	int i,j;
	char *model;
	cs_cosmo_functions_t cosmofns;
	(void)self;   /* silence unused parameter warning */

	if (!PyArg_ParseTuple(args, "dddO!O!s", &Gmu, &f, &Gamma, &PyArray_Type, &Numpy_amp, &PyArray_Type, &Numpy_z, &model))
		return NULL;

	Numpy_amp = PyArray_GETCONTIGUOUS(Numpy_amp);
	Numpy_z = PyArray_GETCONTIGUOUS(Numpy_z);
	if(!Numpy_amp || !Numpy_z){
		Py_DECREF(Numpy_amp);
		Py_DECREF(Numpy_z);
		return NULL;
	}
	Namp = PyArray_DIM(Numpy_amp, 0);
	amp = PyArray_DATA(Numpy_amp);
	Nz = PyArray_DIM(Numpy_z, 0);
	z = PyArray_DATA(Numpy_z);

	{
	npy_intp dims[] = {Namp,Nz};
	Numpy_dRdzdA = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
	}

	cosmofns = XLALCSCosmoFunctions(z, Nz);

	/* loop over amplitudes */
	for ( i = 0; i < Namp; i++ ){
		double A = amp[i];
		/* loop over redshifts */
		for ( j = 0; j < Nz; j++ ){
			double phit = cosmofns.phit[j];
			double phiA = cosmofns.phiA[j];
			double phiV = cosmofns.phiV[j];
			double l = pow ( A / Gmu / H0 * pow(1+z[j],1.0/3.0)* phiA, 3.0/2.0);
			double theta = pow(f*(1+z[j])*l, -1.0/3.0);
			double dRdzdA;

			if(theta > 1.0){
				dRdzdA = 0.0;
			} else {
				double Delta = 0.25*theta*theta;
				dRdzdA = pow(H0,-3.0) * phiV / (1.0+z[j]) * nu(Gmu, Gamma, A, z[j], phit, phiA, model) * Delta;
				if(gsl_isnan(dRdzdA)) {
					Py_DECREF(Numpy_dRdzdA);
					Numpy_dRdzdA = NULL;
					XLALCSCosmoFunctionsFree( cosmofns );
					Py_DECREF(Numpy_amp);
					Py_DECREF(Numpy_z);
					return NULL;
				}
			}
			*(double *) PyArray_GETPTR2((PyArrayObject *) Numpy_dRdzdA, i, j) = dRdzdA;
		}
	}

	XLALCSCosmoFunctionsFree( cosmofns );
	Py_DECREF(Numpy_amp);
	Py_DECREF(Numpy_z);

	return Numpy_dRdzdA;
}


/*******************************************************************************/

//List of functions available to the Python module.
static PyMethodDef cs_gammaMethods[] = {
  {"findzofA", cs_gamma_findzofA, METH_VARARGS,
  "Function to find z(A). From cs_gamma.c; modified to be called from Python."},
  {"finddRdz", cs_gamma_finddRdz, METH_VARARGS,
  "Function to find dR/dz. From cs_gamma.c; modified to be called from Python."},
  {"finddRdzdA", cs_gamma_finddRdzdA, METH_VARARGS,
  "Function to find d^2R/dzdA. From cs_gammaLargeLoops.c; modified to be called from Python."},
  {NULL, NULL, 0, NULL}
};

static PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "cs_gamma", NULL, -1, cs_gammaMethods,
  NULL, NULL, NULL, NULL
};

//They Python module initialization function.
PyMODINIT_FUNC PyInit_cs_gamma(void);	/* silence no-previous-prototype warning */
PyMODINIT_FUNC
PyInit_cs_gamma(void)
{
  import_array();
  return PyModule_Create(&moduledef);
}


SIX_COMPAT_MODULE(cs_gamma)
