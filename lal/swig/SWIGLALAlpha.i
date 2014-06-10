//
// Copyright (C) 2011--2014 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

// Header containing SWIG code which must appear *before* the LAL headers.
// Author: Karl Wette

////////// GSL error handling //////////

// Set custom GSL error handler which raises an XLAL error (instead of aborting).
%header %{
#include <gsl/gsl_errno.h>
static void swiglal_gsl_error_handler(const char *reason, const char *file, int line, int errnum) {
  XLALPrintError("GSL function failed: %s (errnum=%i)\n", reason, errnum);
  XLALError("<GSL function>", file, line, XLAL_EFAILED);
}
%}
%init %{
gsl_set_error_handler(swiglal_gsl_error_handler);
%}

////////// GSL vectors and matrices //////////

// Include GSL headers.
%header %{
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
%}

// This macro create wrapping structs for GSL vectors and matrices.
%define %swig_lal_gsl_vector_matrix(TYPE, NAME)

// GSL vector of type NAME.
typedef struct {
  %extend {
    gsl_vector##NAME(const size_t n) {
      return gsl_vector##NAME##_calloc(n);
    }
    gsl_vector##NAME(gsl_vector##NAME *v0) {
      gsl_vector##NAME *v = gsl_vector##NAME##_alloc(v0->size);
      gsl_vector##NAME##_memcpy(v, v0);
      return v;
    }
    ~gsl_vector##NAME() {
      %swiglal_struct_call_dtor(gsl_vector##NAME##_free, $self);
    }
  }
  %swiglal_array_dynamic_size(size_t, size);
  %swiglal_array_dynamic_1D(gsl_vector##NAME, TYPE, size_t, data, arg1->size, arg1->stride);
} gsl_vector##NAME;

// GSL matrix of type NAME.
typedef struct {
  %extend {
    gsl_matrix##NAME(const size_t n1, const size_t n2) {
      return gsl_matrix##NAME##_calloc(n1, n2);
    }
    gsl_matrix##NAME(gsl_matrix##NAME *m0) {
      gsl_matrix##NAME *m = gsl_matrix##NAME##_alloc(m0->size1, m0->size2);
      gsl_matrix##NAME##_memcpy(m, m0);
      return m;
    }
    ~gsl_matrix##NAME() {
      %swiglal_struct_call_dtor(gsl_matrix##NAME##_free, $self);
    }
  }
  %swiglal_array_dynamic_size(size_t, size1);
  %swiglal_array_dynamic_size(size_t, size2);
  %swiglal_array_dynamic_2D(gsl_matrix##NAME, TYPE, size_t, data, arg1->size1, arg1->size2, arg1->tda, 1);
} gsl_matrix##NAME;

%enddef // %swig_lal_gsl_vector_matrix

// GSL integer vectors and matrices.
%swig_lal_gsl_vector_matrix(short, _short);
%swig_lal_gsl_vector_matrix(unsigned short, _ushort);
%swig_lal_gsl_vector_matrix(int, _int);
%swig_lal_gsl_vector_matrix(unsigned int, _uint);
%swig_lal_gsl_vector_matrix(long, _long);
%swig_lal_gsl_vector_matrix(unsigned long, _ulong);

// GSL real and complex vectors and matrices.
%swig_lal_gsl_vector_matrix(float, _float);
%swig_lal_gsl_vector_matrix(double, ); // GSL double vec./mat. has no typename suffix.
%swig_lal_gsl_vector_matrix(gsl_complex_float, _complex_float);
%swig_lal_gsl_vector_matrix(gsl_complex, _complex);

////////// LAL error handling //////////

// Replace default LAL raise/abort hooks with custom handlers which
// translate them into an XLAL error, so that they will be caught
// by the SWIG %exception handler.
%header %{

#include <signal.h>

// Print the supplied error message, then raise an XLAL error.
static int swig_lal_LALRaise(int sig, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  (void) vfprintf(stderr, fmt, ap);
  va_end(ap);
  (void) fprintf(stderr, "LALRaise: %s\n", strsignal(sig));
  XLALSetErrno(XLAL_EFAILED);
  return 0;
}

// Print the supplied error message, then raise an XLAL error.
static void swig_lal_LALAbort(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  (void) vfprintf(stderr, fmt, ap);
  va_end(ap);
  XLALSetErrno(XLAL_EFAILED);
}

%} // %header
%init %{
lalRaiseHook = swig_lal_LALRaise;
lalAbortHook = swig_lal_LALAbort;
%}

////////// Specialised typemaps for LIGOTimeGPS //////////

// Specialised input typemaps for LIGOTimeGPS structs.
// Accepts a SWIG-wrapped LIGOTimeGPS or a double as input.
%fragment("swiglal_specialised_tagLIGOTimeGPS", "header", fragment=SWIG_AsVal_frag(double)) {
  int swiglal_specialised_tagLIGOTimeGPS(SWIG_Object in, LIGOTimeGPS *out) {
    double val = 0;
    int res = SWIG_AsVal(double)(in, &val);
    if (!SWIG_IsOK(res)) {
      return res;
    }
    XLALGPSSetREAL8(out, val);
    return SWIG_OK;
  }
}
%swiglal_specialised_typemaps(tagLIGOTimeGPS, "swiglal_specialised_tagLIGOTimeGPS");

// Local Variables:
// mode: c
// End:
