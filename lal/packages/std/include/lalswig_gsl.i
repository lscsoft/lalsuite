//
//  Copyright (C) 2011, 2012 Karl Wette
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with with program; see the file COPYING. If not, write to the
//  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//  MA  02111-1307  USA
//

// SWIG interface to GSL structures
// Author: Karl Wette

// Only in SWIG interface.
#ifdef SWIG

// Exclude from preprocessing interface.
#ifndef SWIGXML

// Include basic GSL headers.
%header %{
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
%}

////////// Error handling //////////

// Set custom GSL error handler which raises an XLAL error (instead of aborting).
%header %{
  void swiglal_gsl_error_handler(const char *reason, const char *file, int line, int errnum) {
    XLALPrintError("GSL function failed: %s (errnum=%i)\n", reason, errnum);
    XLALError("<GSL function>", file, line, XLAL_EFAILED);
  }
%}
%init %{
  gsl_set_error_handler(swiglal_gsl_error_handler);
%}

////////// GSL vectors and matrices //////////

// This macro create wrapping structs for GSL vectors and matrices.
%define %lalswig_gsl_vector_matrix(TYPE, NAME)

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
        %swiglal_call_dtor(gsl_vector##NAME##_free, $self);
      }
    }
    %swiglal_array_dynamic_1D(TYPE, size_t, data, size, arg1->stride);
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
        %swiglal_call_dtor(gsl_matrix##NAME##_free, $self);
      }
    }
    %swiglal_array_dynamic_2D(TYPE, size_t, data, size1, size2, arg1->tda, 1);
  } gsl_matrix##NAME;

%enddef // %lalswig_gsl_vector_matrix

// GSL integer vectors and matrices.
%lalswig_gsl_vector_matrix(short, _short);
%lalswig_gsl_vector_matrix(unsigned short, _ushort);
%lalswig_gsl_vector_matrix(int, _int);
%lalswig_gsl_vector_matrix(unsigned int, _uint);
%lalswig_gsl_vector_matrix(long, _long);
%lalswig_gsl_vector_matrix(unsigned long, _ulong);

// GSL real and complex vectors and matrices.
%lalswig_gsl_vector_matrix(float, _float);
%lalswig_gsl_vector_matrix(double, ); // GSL double vec./mat. has no typename suffix.
%lalswig_gsl_vector_matrix(gsl_complex_float, _complex_float);
%lalswig_gsl_vector_matrix(gsl_complex, _complex);

////////// GSL random number generators //////////

// GSL random number generator
%header %{
typedef gsl_rng GSL_rng;
%}
typedef struct {
  %extend {

    // Constructor
    GSL_rng(const char* name, unsigned long int seed) {

      // Check input
      XLAL_CHECK_NULL(name != NULL, XLAL_EFAULT, "Generator name must be non-NULL");

      // Read environment variables for default generators
      gsl_rng_env_setup();

      // Find generator
      const gsl_rng_type* T = NULL;
      if (strcmp(name, "default") == 0) {
        T = gsl_rng_default;
      } else {
        const gsl_rng_type **types = gsl_rng_types_setup();
        for (const gsl_rng_type **t = types; *t != NULL; ++t) {
          if (strcmp(name, (*t)->name) == 0) {
            T = *t;
            break;
          }
        }
      }
      XLAL_CHECK_NULL(T != NULL, XLAL_EINVAL, "Could not find generator named '%s'", name);

      // Create generator and set seed
      gsl_rng* rng = gsl_rng_alloc(T);
      gsl_rng_set(rng, seed);

      return rng;

    }

    // Copy constructor
    GSL_rng(const GSL_rng* src) {

      // Check input
      XLAL_CHECK_NULL(src != NULL, XLAL_EFAULT, "Generator must be non-NULL");

      // Clone generator
      return gsl_rng_clone(src);

    }

    // Destructor
    ~GSL_rng() {
      %swiglal_call_dtor(gsl_rng_free, $self);
    }

    // Properties and methods
    void set_seed(unsigned long int seed) {
      gsl_rng_set($self, seed);
    }
    unsigned long int get_value() {
      return gsl_rng_get($self);
    }
    double uniform() {
      return gsl_rng_uniform($self);
    }
    double uniform_pos() {
      return gsl_rng_uniform_pos($self);
    }
    unsigned long int uniform_int(unsigned long int n) {
      return gsl_rng_uniform_int($self, n);
    }
    const char* name() {
      return gsl_rng_name($self);
    }
    unsigned long int max_value() {
      return gsl_rng_max($self);
    }
    unsigned long int min_value() {
      return gsl_rng_min($self);
    }

  }
} GSL_rng;

#endif // !SWIGXML

#endif // SWIG
