/*
*  Copyright (C) 2007 Xavier Siemens
*  Copyright (C) 2010 Andrew Mergl
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
/*            Cosmic string burst rate computation code for small loops          */
/*                                                                               */
/*                  Xavier Siemens, Jolien Creighton, Irit Maor                  */
/*                                                                               */
/*                         UWM/Caltech - September 2006                          */
/*********************************************************************************/
/*Modified June 2010 by Andrew Mergl for use with Python*/

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <lal/cs_cosmo.h>
#include <lal/cs_lambda_cosmo.h>

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
  int Namp;
  cs_cosmo_functions_t cosmofns;
  int j;

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
  dRdz = PyArray_DATA(Numpy_dRdz);

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
  int Namp;

  double z_min = 1e-20, z_max = 1e10;
  double dlnz = 0.05;
  unsigned numz = floor( (log(z_max) - log(z_min)) / dlnz );
  int i;
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
  zofA = PyArray_DATA(Numpy_zofA);

  cosmofns = XLALCSCosmoFunctionsAlloc( z_min, dlnz, numz );

  zofa_interp = gsl_interp_alloc (gsl_interp_linear, cosmofns.n);

  fz = calloc( cosmofns.n, sizeof( *fz ) ); 
  z = calloc( cosmofns.n, sizeof( *z ) ); 
  
  /* first compute the function that relates A and z */
  /* invert order; b/c fz is a monotonically decreasing func of z */
  for ( i = cosmofns.n-1 ; i >= 0; i-- )
    {
      int j = cosmofns.n-1 - i;
      z[j] = cosmofns.z[i];
      fz[j] = pow(cosmofns.phit[i], 2.0/3.0) * pow(1+z[j], -1.0/3.0) / cosmofns.phiA[i];
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

//List of functions available to the Python module.
static PyMethodDef cs_gammaMethods[] = {
  {"findzofA", cs_gamma_findzofA, METH_VARARGS,
  "Function to find z(A). From cs_gamma.c; modified to be called from Python."},
  {"finddRdz", cs_gamma_finddRdz, METH_VARARGS,
  "Function to find dR/dz. From cs_gamma.c; modified to be called from Python."},
  {NULL, NULL, 0, NULL}
};

//They Python module initialization function.
PyMODINIT_FUNC
initcs_gamma(void)
{
  Py_InitModule("cs_gamma", cs_gammaMethods);
  import_array();
}
