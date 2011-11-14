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

// LALSimulation header files which form the SWIG interface
%apply REAL8 *OUTPUT { REAL8 *gamma00, REAL8 *gamma01, REAL8 *gamma02, REAL8 *gamma11, REAL8 *gamma12, REAL8 *gamma22 };
