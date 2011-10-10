// SWIG interface to the LALSimulation library

// common SWIG interface code
%include <lal/swiglal-common.i>

// import LAL definitions
%import <lal/swig-lal.i>

/* All of the functions in lalsimulation that return int use it as an error
 * code, but the swiglal-common.i already registers an error handler that
 * picks up XLAL_ERROR calls and raises Python exceptions.
 *
 * Aside: It's possible to return an error code but forget to call XLAL_ERROR;
 * this miss that situation.
 */
%typemap(out) int ""

/*
 * Every occurrence of the arguments
 *   REAL8TimeSeries **hplus
 * and
 *   REAL8TimeSeries **hcross
 * represents an output argument.  Apply the double-pointer typemap from
 * swiglal-common.i.
 */
%apply SWIGTYPE **OUTPUT {REAL8TimeSeries **hplus, REAL8TimeSeries **hcross};

/*
 * Every occurrence of the arguments
 *   COMPLEX16FrequencySeries **htilde
 * represents an output argument.  Apply the double-pointer typemap from
 * swiglal-common.i.
 */
%apply SWIGTYPE **OUTPUT {COMPLEX16FrequencySeries **htilde};

// allocate these structs with LAL memory
SWIGLAL_STRUCT_LALMEM(expnCoeffsdEnergyFlux);
SWIGLAL_STRUCT_LALMEM(TofVIntegrandIn);
SWIGLAL_STRUCT_LALMEM(TofVIn);

// LALSimulation header files which form the SWIG interface
