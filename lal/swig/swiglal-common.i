//
//  Copyright (C) 2011 Karl Wette
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

// common SWIG interface code
// Author: Karl Wette, 2011

/////////////// Basic definitions and includes ///////////////

// Ensure that NDEBUG is defined iff SWIGLAL_NDEBUG is defined
%begin %{
  #undef NDEBUG
  #ifdef SWIGLAL_NDEBUG
    #define NDEBUG
  #endif
%}

// SWIG version used to generate wrapping code
%inline %{
  const int swig_version_hex = SWIGVERSION;
%}

// Was the wrapping code generated in debugging mode?
#ifdef SWIGLAL_NDEBUG
  %inline %{
    const bool swiglal_debug = false;
  %}
#else
  %inline %{
    const bool swiglal_debug = true;
  %}
#endif

// Ignore GCC warnings about unitialised or unused variables.
// These warnings may be produced by some SWIG-generated code.
%begin %{
  #ifdef __GNUC__
    #pragma GCC diagnostic ignored "-Wuninitialized"
    #pragma GCC diagnostic ignored "-Wunused-variable"
  #endif
%}

// Include basic C++ and LAL headers in wrapping code.
%header %{
  #include <cstdlib>
  #include <cstring>
  #include <iostream>
  #include <string>
  #include <sstream>
  #include <lal/XLALError.h>
  #include <lal/LALMalloc.h>
%}

// Allow SWIG wrapping code can raise exceptions.
%include <exception.i>

// Tell SWIG about C99 integer typedefs,
// on which LAL integer types are based.
%include <stdint.i>

// Generate copy constructors for all structs.
%copyctor;

// Turn on auto-documentation of functions.
%feature("autodoc", 1);

// Enable keyword argument globally, if
// supported by the target scripting language.
%feature("kwargs", 1);

// Suppress warnings about not being able to use
// keywords arguments in various circumstances.
#pragma SWIG nowarn=SWIGWARN_PARSE_KEYWORD
#pragma SWIG nowarn=SWIGWARN_LANG_VARARGS_KEYWORD
#pragma SWIG nowarn=SWIGWARN_LANG_OVERLOAD_KEYWORD

// This macro is used to both %include a LAL header into the
// SWIG interface, so that SWIG will generate wrappings for
// its contents, and also #include it in the wrapping code.
%define swiglal_include(HEADER)
  %include <HEADER>
  %header %{
    #include <HEADER>
  %}
%enddef

// Remove LAL RCS ID macros from SWIG interface.
#define NRCSID(name,id)
#define RCSID(id)

// So that SWIG knows about basic LAL datatypes.
%header %{
  #include <lal/LALAtomicDatatypes.h>
  #include <lal/LALComplex.h>
%}

// So that SWIG wrapping code knows about basic GSL types.
%header %{
  #include <gsl/gsl_complex_math.h>
  #include <gsl/gsl_vector.h>
  #include <gsl/gsl_matrix.h>
  // GSL doesn't provide a constructor function for
  // gsl_complex_float, so we provide one here.
  SWIGINTERN gsl_complex_float gsl_complex_float_rect(float x, float y) {
    gsl_complex_float z;
    GSL_SET_COMPLEX(&z, x, y);
    return z;
  }
%}

// Function which tests whether the pointer passed to it is non-zero.
// This function does the right thing if it is passed an actual pointer,
// or the name of a statically-allocated array (which is implicitly
// converted into an address), and it will also correctly fail to compile
// if its argument is neither of these possibilities.
%header %{
  template<class T> SWIGINTERNINLINE bool swiglal_check_ptr(T *p) {
    return p != NULL;
  }
%}

// Wrapper function for SWIG_exception which allows MESSAGE to
// contain streaming (<<) operators; this allows error messages
// to easily include non-string data.
%define swiglal_exception(CODE, MESSAGE)
  std::stringstream msg;
  msg << MESSAGE;
  SWIG_exception(CODE, msg.str().c_str());
%enddef

// Include SWIG interface code for specific scripting languages.
#ifdef SWIGOCTAVE
%include <lal/swiglal-octave.i>
#endif
#ifdef SWIGPYTHON
%include <lal/swiglal-python.i>
#endif

/////////////// Memory allocation ///////////////

// SWIG generates default constructor and destructors for wrapped structs.
// These contructors/destructors always use 'malloc'/'free' (C mode) or
// 'new'/'delete' (C++ mode); see e.g. Source/Swig/cwrap.c in SWIG 1.3.40.
// It is not possible to specify alternative memory allocators through any
// options to SWIG. However, we would like all LAL structs to be allocated
// with LAL memory allocation functions, so that the LAL memory debugging
// system can be used inside a scripting language. For example, if a LAL
// struct allocated by SWIG is passed to a LAL function which tries to
// de-allocate it, the LAL memory debugging system will fail unless the
// LAL struct was originally allocated by LAL.
//
// The solution is to use a C++ feature: a class C can supply their own
// 'new'/'delete' operators, which will then be called by 'new C()' to
// allocate the class. We do the same for all LAL structs by adding the
// symbol 'SWIGLAL_STRUCT_LALALLOC' to their definitions. This symbol
// defines a 'new' operator, which allocates the struct using XLALCalloc
// (so that the struct is zero-initialised), and a 'delete' operator
// which calls XLALFree to de-allocate the struct.
//
// Note that, if a LAL struct does not contain 'SWIGLAL_STRUCT_LALALLOC',
// it will continue to be use the default 'new'/'delete' operators. This
// should break anything as long as LAL memory debugging is not used,
// either by compiling with --disable-debug or by setting lalDebugLevel=0.

// Remove SWIGLAL_STRUCT_LALALLOC from SWIG interface
#define SWIGLAL_STRUCT_LALALLOC()

// Define SWIGLAL_STRUCT_LALALLOC in SWIG wrapping code
%header %{
  #define SWIGLAL_STRUCT_LALALLOC() \
    void* operator new (size_t n) throw() { \
      return XLALCalloc(1, n); \
    } \
    void operator delete(void* p) { \
      XLALFree(p); \
    }
%}
