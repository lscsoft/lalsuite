/*
 *  Copyright (C) 2012 Nickolas Fotopoulos, Leo Singer, Alexander Dietz,
 *   Justin Wagner
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


#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <lal/DetResponse.h>
#include <lal/LALConstants.h>
#include <lal/LALDetectors.h>
#include <time.h>
#include <ctype.h>

// #ifdef _OPENMP
// #include <omp.h>
// #endif

static ssize_t str2network(LALDetector network[], size_t net_size, char *str);
static double rchisq_2(gsl_rng *rng, double lambda);

static void simulate(size_t *num_found, LALDetector *network, size_t network_size, double* BNS_horizon_dist, double jet_semiangle_deg, double m1_min, double m1_max, double m2_min, double m2_max, double rho_thresh, double *distances, size_t dist_bins, size_t samples){
  double beam_fac=1-cos(jet_semiangle_deg)*M_PI/180.0;
  static double twopi=2 * M_PI;
  const double BNS_mchirp_fac= 1.4 * 1.4 / cbrt(2.8);  /* mc^{5/3} */
  const double rho_thresh_2=rho_thresh*rho_thresh;
  memset(num_found, 0, dist_bins * sizeof(*num_found));

// #pragma omp parallel
  {
    /*Create and seed random number generator*/
    gsl_rng *rng=gsl_rng_alloc(gsl_rng_mt19937);
// #ifdef _OPENMP
//     gsl_rng_set(rng, omp_get_thread_num()*time(NULL));//Seed with thread number.
// #else
    gsl_rng_set(rng, time(NULL));
// #endif

    /* create per-thread workspace */
    double *lambda_const = NULL;
    size_t *threadTotals = NULL;
    lambda_const = malloc(2 * network_size * sizeof(*lambda_const));
    if (!lambda_const) goto done;
    threadTotals = calloc(dist_bins, sizeof(*threadTotals));
    if (!threadTotals) goto done;

    size_t j, k, l;

// #pragma omp for schedule(static)
    for(k=0; k<samples; ++k){
      /* draw orientations and determine network response */
      const double lon=twopi*gsl_rng_uniform(rng);
      const double lat=M_PI_2 - acos(2*gsl_rng_uniform(rng) - 1);
      const double zeta=twopi*gsl_rng_uniform(rng);

      /* draw cos iota */
      const double cosiota=1-beam_fac*gsl_rng_uniform(rng);//beam_fac determines max iota.
      const double cosiotasq=cosiota*cosiota;
      const double iotafac=0.25*(1+cosiotasq)*(1+cosiotasq);

      /* draw masses */
      const double m1=m1_min+gsl_rng_uniform(rng)*(m1_max-m1_min);
      const double m2=m2_min+gsl_rng_uniform(rng)*(m2_max-m2_min);
      /* amplitude^2 scales as mchirp^{5/3} */
      const double mass_correction = m1 * m2 / (cbrt(m1 + m2) * BNS_mchirp_fac);

      /* gather distance-independent amplitude factors */
      for(l=network_size; l--;) {
        double fplus = 0, fcross = 0;
        XLALComputeDetAMResponse(&fplus, &fcross, (const REAL4 (*)[3])network[l].response, lon, lat, zeta, 0.0);
        const double S2 = (fplus*fplus)*iotafac+(fcross*fcross)*cosiotasq;
        lambda_const[l] = (64 - 2) * BNS_horizon_dist[l] * BNS_horizon_dist[l] * S2 * mass_correction;
      }

      /* determine detections */
      for(j=dist_bins; j--;){
        int successes=0;
        const double dist_sq = distances[j] * distances[j];
        for(l=network_size; l--;){
          /* determine non-centrality parameter lambda (~= expected SNR^2) */
          const double lambda = lambda_const[l] / dist_sq;
          /* draw actual SNR^2 value from distribution over noise */
          const double rand = rchisq_2(rng,lambda);
          /* Require two or more detectors to detect a signal
             for a valid detection to be added to the total. */
          successes+=rand>rho_thresh_2;
        }
        threadTotals[j]+=(successes>=2);
      }
    }

// #pragma omp critical
    for(j = 0; j < dist_bins; ++j) {
      num_found[j] += threadTotals[j];
    }

done:
    gsl_rng_free(rng);
    free(threadTotals);
    free(lambda_const);
  }
}

static double rchisq_2(gsl_rng *rng, double lambda) {
  /* Generates a random value from a degree 2
     chi_squared with non-centrality parameter of lambda*/
  double a=sqrt(0.5*lambda);
  const double temp=gsl_ran_gaussian(rng, 1.0)+a;
  const double temp2=gsl_ran_gaussian(rng, 1.0)+a;
  return (temp*temp)+(temp2*temp2);
}

static ssize_t str2network(LALDetector network[LAL_NUM_DETECTORS], size_t net_size, char *str){
  /*Convert string like "HLVK" to an array of LALDetectors.
    Return the size of the network.*/
  ssize_t k=0;
  while (k < net_size) {
    if (str[k]=='H') {
      network[k++]=lalCachedDetectors[LAL_LHO_4K_DETECTOR];
    }
    else if (str[k]=='L') {
      network[k++]=lalCachedDetectors[LAL_LLO_4K_DETECTOR];
    }
    else if (str[k]=='V') {
      network[k++]=lalCachedDetectors[LAL_VIRGO_DETECTOR];
    }
    else if (str[k]=='K') {
      /* numbers from private communication with Koji Arai */
      LALDetector *detector=network+k;
      LALFrDetector *frDetector=&(detector->frDetector);
      strncpy(frDetector->name, "KAGRA", LALNameLength);
      strncpy(frDetector->prefix, "K1", 3);
      frDetector->vertexLatitudeRadians=2.396511595913414;
      frDetector->vertexLongitudeRadians=0.6354743806511354;
      frDetector->vertexElevation=372.0;
      frDetector->xArmAltitudeRadians=0.0;
      frDetector->xArmAzimuthRadians=1.076693615555302;
      frDetector->yArmAltitudeRadians=0.0;
      frDetector->yArmAzimuthRadians=5.789082595939991;
      detector=XLALCreateDetector(network+k, frDetector, LALDETECTORTYPE_IFODIFF);
      if (!detector) {
        fprintf(stderr, "Failed to create KAGRA detector\n");
        return -1;
      }
      k++;
    }
    else if (str[k] == 'I') {
      /* numbers from Schutz 2011 network FOMs */
      LALDetector *detector=network+k;
      LALFrDetector *frDetector=&(detector->frDetector);
      strncpy(frDetector->name, "Indigo", LALNameLength);
      strncpy(frDetector->prefix, "I1", 3);
      frDetector->vertexLatitudeRadians=1.3098647554849334;
      frDetector->vertexLongitudeRadians=0.33329486135237268;
      frDetector->vertexElevation=0.0;
      frDetector->xArmAltitudeRadians=0.0;
      frDetector->xArmAzimuthRadians=3.9269908169872414;
      frDetector->yArmAltitudeRadians=0.0;
      frDetector->yArmAzimuthRadians=5.497787143782138;
      detector=XLALCreateDetector(network+k, frDetector, LALDETECTORTYPE_IFODIFF);
      if (!detector) {
        fprintf(stderr, "Failed to create Indigo detector\n");
        return -1;
      }
      k++;
    }
    else {
      fprintf(stderr, "Unrecognized site: %c\n", str[k]);
      return -1;
    }
  }
  return k;
}

static PyObject *pylalsimulate(PyObject *self, PyObject *args){
  PyObject *py_network, *distances, *out_array = NULL;
  double m1_min, m1_max;
  double m2_min, m2_max;
  double rho_thresh, jet_semiangle_deg;
  Py_ssize_t samples, network_length, i;
  int allocated_distances = 0;
  char *net_names = NULL;
  double *net_hor_dists = NULL;
  LALDetector *network = NULL;

  if(!PyArg_ParseTuple(args, "OdddddOdn", &py_network, &jet_semiangle_deg, &m1_min, &m1_max, &m2_min, &m2_max, &distances, &rho_thresh, &samples)) return NULL;

  /* make sure we have a distance array */
  if (!PyArray_Check(distances) || (PyArray_TYPE(distances) != NPY_DOUBLE)
      || !PyArray_ISCARRAY_RO(distances)) {
      distances = PyArray_FROMANY(distances, NPY_DOUBLE, 1, 1, NPY_CARRAY_RO);
      allocated_distances = 1;
  }

  /* build network */
  if (!PySequence_Check(py_network)) {
    PyErr_SetString(PyExc_ValueError, "require network to be a sequence");
    goto done;
  }
  network_length = PySequence_Length(py_network);
  net_names = malloc(network_length * sizeof(*net_names));
  net_hor_dists = malloc(network_length * sizeof(*net_hor_dists));
  for(i=0; i<network_length; i++){
    PyObject *name_dist_tup = PySequence_ITEM(py_network, i);
    if (!name_dist_tup || !PySequence_Check(name_dist_tup) || PySequence_Length(name_dist_tup) != 2) {
      Py_XDECREF(name_dist_tup);
      PyErr_SetString(PyExc_ValueError, "network must be composed of 2-tuples");
      goto done;
    }
    PyObject *name = PySequence_ITEM(name_dist_tup, 0);
    if (!name || !PyString_Check(name) || PyString_Size(name) != 1) {
      Py_XDECREF(name);
      Py_DECREF(name_dist_tup);
      PyErr_SetString(PyExc_ValueError, "network name must be a single-character string");
      goto done;
    }
    net_names[i] = toupper(PyString_AsString(name)[0]);
    Py_DECREF(name);
    PyObject *dist = PySequence_ITEM(name_dist_tup, 1);
    if (!dist || !PyNumber_Check(dist)) {
      Py_XDECREF(dist);
      Py_DECREF(name_dist_tup);
      PyErr_SetString(PyExc_ValueError, "network distance must be a number");
      goto done;
    }
    PyObject *dist_float = PyNumber_Float(dist);
    if (!dist_float) {
      Py_DECREF(dist);
      Py_DECREF(name_dist_tup);
      PyErr_SetString(PyExc_ValueError, "could not cast network distance to float");
      goto done;
    }
    net_hor_dists[i] = PyFloat_AsDouble(dist);
    Py_DECREF(dist_float);
    Py_DECREF(dist);
    Py_DECREF(name_dist_tup);
  }

  /* map detector names to detectors */
  network = malloc(network_length * sizeof(*network));
  if (str2network(network, network_length, net_names) != network_length){
    PyErr_SetString(PyExc_ValueError, "Failed to build network");
    goto done;
  }

  /* do simulation */
  out_array = PyArray_SimpleNew(1, PyArray_DIMS(distances), NPY_UINTP);
  simulate(PyArray_DATA(out_array), network, network_length, net_hor_dists, jet_semiangle_deg, m1_min, m1_max, m2_min, m2_max, rho_thresh, PyArray_DATA(distances), PyArray_SIZE(distances), samples);

done:

  if (allocated_distances) { Py_DECREF(distances); }

  free(network);
  free(net_names);
  free(net_hor_dists);

  return out_array;

}

const char docstring[] =
  "Simulates the detection pipeline of a detector network using\n"
  "requiring a response above the given SNR threshold in two or more\n"
  "detectors. Returns the number of signals found in each distance bin.\n"
  "========================\n"
  "       Parameters\n"
  "========================\n"
  "Network: Sequence of (Site letter, Horizon distance (Mpc) for 1.4, 1.4 Msun binary) tuples.\n"
  "jet_semiangle_deg: Maximum jet angle iota (degrees)\n"
  "m1_min: Min mass for m1 (Msun).\n"
  "m1_max: Max mass for m1 (Msun).\n"
  "m2_min: Min mass for m2 (Msun).\n"
  "m2_max: Max mass for m2 (Msun).\n"
  "distances: Distances (Mpc) at which to evaluate efficiency.\n"
  "rho_thresh: SNR threshold (required in two or more detectors).\n"
  "N: Number of simulations.\n";

static struct PyMethodDef methods[] = {
  {"simulate", pylalsimulate, METH_VARARGS, docstring},
  {NULL, NULL, 0, NULL}/* Sentinel */
};

void initcbc_network_efficiency(void) {
  (void) Py_InitModule3("pylal.cbc_network_efficiency", methods, docstring);
  import_array();
}
