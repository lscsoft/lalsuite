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

// Common SWIG interface code.
// Author: Karl Wette

////////// General SWIG directives //////////

// Ensure that all LAL library wrapping modules share type information.
%begin %{
#define SWIG_TYPE_TABLE swiglal
%}

// Include SWIG headers.
%include <exception.i>
%include <typemaps.i>

// Suppress some SWIG warnings.
#pragma SWIG nowarn=SWIGWARN_PARSE_KEYWORD
#pragma SWIG nowarn=SWIGWARN_LANG_VARARGS_KEYWORD
#pragma SWIG nowarn=SWIGWARN_LANG_OVERLOAD_KEYWORD

// Turn on auto-documentation of functions.
%feature("autodoc", 1);

// Enable keyword arguments for functions.
%feature("kwargs", 1);

////////// Public macros //////////

// Public macros (i.e. those used in headers) are contained inside
// the SWIGLAL() macro, so that they can be easily removed from the
// proprocessing interface. The SWIGLAL_CLEAR() macro is used to
// clear any effects of a public macro. Calls to SWIGLAL_CLEAR()
// are generated from the preprocessing interface.
#define SWIGLAL(...) %swiglal_public_##__VA_ARGS__
#define SWIGLAL_CLEAR(...) %swiglal_public_clear_##__VA_ARGS__

////////// Utility macros and typemaps //////////

// Define a SWIG preprocessor symbol NAME1_NAME2
%define %swiglal_define2(NAME1, NAME2)
%define NAME1##_##NAME2 %enddef
%enddef

// The macro %swiglal_map() maps a one-argument MACRO(X) onto a list
// of arguments (which may be empty). Based on SWIG's %formacro().
%define %_swiglal_map(MACRO, X, ...)
#if #X != ""
MACRO(X);
%_swiglal_map(MACRO, __VA_ARGS__);
#endif
%enddef
%define %swiglal_map(MACRO, ...)
%_swiglal_map(MACRO, __VA_ARGS__, );
%enddef

// The macro %swiglal_map_a() maps a two-argument MACRO(A, X) onto a list
// of arguments (which may be empty), with a common first argument A.
%define %_swiglal_map_a(MACRO, A, X, ...)
#if #X != ""
MACRO(A, X);
%_swiglal_map_a(MACRO, A, __VA_ARGS__);
#endif
%enddef
%define %swiglal_map_a(MACRO, A, ...)
%_swiglal_map_a(MACRO, A, __VA_ARGS__, );
%enddef

// The macro %swiglal_map_ab() maps a three-argument MACRO(A, B, X) onto a list
// of arguments (which may be empty), with common first arguments A and B.
%define %_swiglal_map_ab(MACRO, A, B, X, ...)
#if #X != ""
MACRO(A, B, X);
%_swiglal_map_ab(MACRO, A, B, __VA_ARGS__);
#endif
%enddef
%define %swiglal_map_ab(MACRO, A, B, ...)
%_swiglal_map_ab(MACRO, A, B, __VA_ARGS__, );
%enddef

// The macro %swiglal_map_abc() maps a four-argument MACRO(A, B, C, X) onto a list
// of arguments (which may be empty), with common first arguments A, B, and C.
%define %_swiglal_map_abc(MACRO, A, B, C, X, ...)
#if #X != ""
MACRO(A, B, C, X);
%_swiglal_map_abc(MACRO, A, B, C, __VA_ARGS__);
#endif
%enddef
%define %swiglal_map_abc(MACRO, A, B, C, ...)
%_swiglal_map_abc(MACRO, A, B, C, __VA_ARGS__, );
%enddef

// Apply and clear SWIG typemaps.
%define %swiglal_apply(TYPEMAP, TYPE, NAME)
%apply TYPEMAP { TYPE NAME };
%enddef
%define %swiglal_clear(TYPE, NAME)
%clear TYPE NAME;
%enddef

// Apply a SWIG feature.
%define %swiglal_feature(FEATURE, VALUE, NAME)
%feature(FEATURE, VALUE) NAME;
%enddef
%define %swiglal_feature_nspace(FEATURE, VALUE, NSPACE, NAME)
%feature(FEATURE, VALUE) NSPACE::NAME;
%enddef

// Suppress a SWIG warning.
%define %swiglal_warnfilter(WARNING, NAME)
%warnfilter(WARNING) NAME;
%enddef
%define %swiglal_warnfilter_nspace(WARNING, NSPACE, NAME)
%warnfilter(WARNING) NSPACE::NAME;
%enddef

// Macros for allocating/copying new instances and arrays
// Analogous to SWIG macros but using XLAL memory functions
#define %swiglal_new_instance(TYPE...) \
  %reinterpret_cast(XLALMalloc(sizeof(TYPE)), TYPE*)
#define %swiglal_new_copy(VAL, TYPE...) \
  %reinterpret_cast(memcpy(%swiglal_new_instance(TYPE), &VAL, sizeof(TYPE)), TYPE*)
#define %swiglal_new_array(SIZE, TYPE...) \
  %reinterpret_cast(XLALMalloc((SIZE)*sizeof(TYPE)), TYPE*)
#define %swiglal_new_copy_array(PTR, SIZE, TYPE...) \
  %reinterpret_cast(memcpy(%swiglal_new_array(SIZE, TYPE), PTR, sizeof(TYPE)*(SIZE)), TYPE*)

////////// General interface code //////////

// Include C99/C++ headers.
%header %{
#ifdef __cplusplus
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <complex>
#else
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#endif
%}

// Remove LAL macros for certain keywords.
#define _LAL_INLINE_
#define _LAL_RESTRICT_

// If necessary, redefine 'restrict' keyword for C++ code.
%header %{
#ifdef __cplusplus
# ifndef restrict
#  if defined __GNUC__
#   define restrict __restrict__
#  else
#   define restrict
#  endif
# endif
#endif
%}

// Include SWIG configuration header generated from 'config.h'.
%header %{
#include <swig_config.h>
%}

// Include LAL headers.
%header %{
#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/XLALError.h>
#include <lal/Date.h>
%}

// Version of SWIG used to generate wrapping code.
%constant int swig_version = SWIG_VERSION;

// Constructors for GSL complex numbers, if required.
%header %{
#ifdef SWIGLAL_HAVE_LIBGSL
#include <gsl/gsl_complex_math.h>   // provides gsl_complex_rect()
SWIGINTERNINLINE gsl_complex_float gsl_complex_float_rect(float x, float y) {
  gsl_complex_float z;
  GSL_SET_COMPLEX(&z, x, y);
  return z;
}
#endif
%}

// Convert XLAL/LAL errors into native scripting-language exceptions:
//  - XLAL: Before performing any action, clear the XLAL error number.
//    Check it after the action is performed; if it is non-zero, raise
//    a native scripting-language exception with the appropriate XLAL
//    error message.
//  - LAL: For any function which has a LALStatus* argument, create a
//    blank LALStatus struct and pass its pointer to the LAL function.
//    After the function is called, check the LALStatus statusCode;
//    if it is non-zero, raise a generic XLAL error (to set the XLAL
//    error number), then a native scripting-language exception. The
//    swiglal_check_LALStatus (C) preprocessor symbol determines if
//    the LAL error handling code is invoked; by default this symbol
//    should be undefined, i.e. functions are XLAL functions by default.
%header %{
static const LALStatus swiglal_empty_LALStatus = {0, NULL, NULL, NULL, NULL, 0, NULL, 0};
#undef swiglal_check_LALStatus
%}
%typemap(in, noblock=1, numinputs=0) LALStatus* {
  LALStatus lalstatus = swiglal_empty_LALStatus;
  $1 = &lalstatus;
%#define swiglal_check_LALStatus
}
%exception %{
  XLALClearErrno();
  $action
#ifdef swiglal_check_LALStatus
  if (lalstatus.statusCode) {
    XLALSetErrno(XLAL_EFAILED);
    SWIG_exception(SWIG_RuntimeError, lalstatus.statusDescription);
  }
#else
  if (xlalErrno) {
    SWIG_exception(SWIG_RuntimeError, XLALErrorString(xlalErrno));
  }
#endif
#undef swiglal_check_LALStatus
%}

////////// General fragments //////////

// Empty fragment, for fragment-generating macros.
%fragment("swiglal_empty_frag", "header") {}

// Wrappers around SWIG's pointer to/from SWIG-wrapped scripting language object functions.
// swiglal_from_SWIGTYPE() simply returns a SWIG-wrapped object containing the input pointer.
// swiglal_from_SWIGTYPE() extracts a pointer from a SWIG-wrapped object, then struct-copies
// the pointer to the supplied output pointer.
%fragment("swiglal_from_SWIGTYPE", "header") {
  SWIGINTERNINLINE SWIG_Object swiglal_from_SWIGTYPE(SWIG_Object self, void *ptr, bool isptr, swig_type_info *tinfo, int tflags) {
    return SWIG_NewPointerObj(isptr ? *((void**)ptr) : ptr, tinfo, tflags);
  }
}
%fragment("swiglal_as_SWIGTYPE", "header") {
  SWIGINTERN int swiglal_as_SWIGTYPE(SWIG_Object self, SWIG_Object obj, void *ptr, size_t len, bool isptr, swig_type_info *tinfo, int tflags) {
    void *vptr = NULL;
    int ecode = SWIG_ConvertPtr(obj, &vptr, tinfo, tflags);
    if (!SWIG_IsOK(ecode)) {
      return ecode;
    }
    memcpy(ptr, isptr ? &vptr : vptr, len);
    return ecode;
  }
}

////////// Generate interface //////////

// The SWIG interface file is generated by 'generate_swiglal_iface.py', which
// extracts lists of interface symbols (functions, structs, typedefs to structs)
// and headers from the preprocessing interface, and generates calls to the
// following macros:

// Process an interface function NAME: rename it to RENAME, and set it to
// always return SWIG-owned wrapping objects (unless the function is being
// ignored). If RETN_1STARG is true, define 'swiglal_return_1starg_##NAME'.
%define %swiglal_process_function(NAME, RENAME, RETN_1STARG)
%rename(#RENAME) NAME;
#if #RENAME != "$ignore"
%feature("new", "1") NAME;
#if RETN_1STARG
%header {#define swiglal_return_1starg_##NAME}
#endif
#endif
%enddef

// Process a typedef to an interface struct TAGNAME: rename it to RENAME.
%define %swiglal_process_tdstruct(TAGNAME, RENAME)
%rename(#RENAME) TAGNAME;
%enddef

// Do not generate any default (copy) contructors or destructors.
%nodefaultctor;
%nocopyctor;
%nodefaultdtor;

// Call the destructor DTORFUNC, but first call swiglal_release_parent()
// to check whether PTR own its own memory (and release any parents).
%define %swiglal_call_dtor(DTORFUNC, PTR)
if (swiglal_release_parent(PTR)) {
  XLALClearErrno();
  (void)DTORFUNC(PTR);
  XLALClearErrno();
}
%enddef

// Generate constructors and destructors for a struct NAME.
%define %swiglal_generate_struct_cdtor(NAME, TAGNAME, OPAQUE, DTORFUNC)

  // If this is an opaque struct:
  #if OPAQUE

    // Create an empty struct to represent the opaque struct,
    // so that SWIG has something to attach the destructor to.
    // No constructors are generated, since it is assumed that
    // the struct will be created internally by some function
    // and returned.
    struct TAGNAME {
    };

  #else

    // Generate constructors, using XLALCalloc() to allocate memory.
    %extend TAGNAME {
      TAGNAME() {
        NAME* self = %reinterpret_cast(XLALCalloc(1, sizeof(NAME)), NAME*);
        return self;
      }
      TAGNAME(NAME* src) {
        NAME* self = %reinterpret_cast(XLALCalloc(1, sizeof(NAME)), NAME*);
        memcpy(self, src, sizeof(NAME));
        return self;
      }
    }

  #endif

  // Generate destructor, using either the destructor function DTORFUNC,
  // or else XLALFree() if DTORFUNC is undefined, to destroy memory.
  #if #DTORFUNC == ""
    %extend TAGNAME {
      ~TAGNAME() {
        %swiglal_call_dtor(XLALFree, $self);
      }
    }
  #else
    %extend TAGNAME {
      ~TAGNAME() {
        %swiglal_call_dtor(%arg(DTORFUNC), $self);
      }
    }
  #endif

%enddef // %swiglal_generate_struct_cdtor()

////////// Fragments and typemaps for arrays //////////

// Map fixed-array types to special variables of their elements,
// e.g. $typemap(swiglal_fixarr_ltype, const int[][]) returns "int".
// Fixed-array types are assumed never to be arrays of pointers.
%typemap(swiglal_fixarr_ltype) SWIGTYPE "$ltype";
%typemap(swiglal_fixarr_ltype) SWIGTYPE[ANY] "$typemap(swiglal_fixarr_ltype, $*type)";
%typemap(swiglal_fixarr_ltype) SWIGTYPE[ANY][ANY] "$typemap(swiglal_fixarr_ltype, $*type)";
%typemap(swiglal_fixarr_tinfo) SWIGTYPE "$&descriptor";
%typemap(swiglal_fixarr_tinfo) SWIGTYPE[ANY] "$typemap(swiglal_fixarr_tinfo, $*type)";
%typemap(swiglal_fixarr_tinfo) SWIGTYPE[ANY][ANY] "$typemap(swiglal_fixarr_tinfo, $*type)";

// The conversion of C arrays to/from scripting-language arrays are performed
// by the following functions:
//  - %swiglal_array_copyin...() copies a scripting-language array into a C array.
//  - %swiglal_array_copyout...() copies a C array into a scripting-language array.
//  - %swiglal_array_viewout...() wraps a C array inside a scripting-language array,
//    if this is supported by the target scripting language.
// These functions accept a subset of the following arguments:
//  - SWIG_Object obj: input scripting-language array.
//  - SWIG_Object parent: SWIG-wrapped object containing parent struct.
//  - void* ptr: pointer to C array.
//  - const size_t esize: size of one C array element, in bytes.
//  - const size_t ndims: number of C array dimensions.
//  - const size_t dims[]: length of each C array dimension.
//  - const size_t strides[]: strides of C array dimension, in number of elements.
//  - swig_type_info *tinfo: SWIG type info of the C array element datatype.
//  - const int tflags: SWIG type conversion flags.
// Return values are:
//  - %swiglal_array_copyin...(): int: SWIG error code.
//  - %swiglal_array_copyout...() SWIG_Object: output scripting-language array.
//  - %swiglal_array_viewout...() SWIG_Object: output scripting-language array view.

// Names of array conversion functions for array type ACFTYPE.
#define %swiglal_array_copyin_func(ACFTYPE) swiglal_array_copyin_##ACFTYPE
#define %swiglal_array_copyout_func(ACFTYPE) swiglal_array_copyout_##ACFTYPE
#define %swiglal_array_viewout_func(ACFTYPE) swiglal_array_viewout_##ACFTYPE

// Names of fragments containing the conversion functions for ACFTYPE.
#define %swiglal_array_copyin_frag(ACFTYPE) "swiglal_array_copyin_" %str(ACFTYPE)
#define %swiglal_array_copyout_frag(ACFTYPE) "swiglal_array_copyout_" %str(ACFTYPE)
#define %swiglal_array_viewout_frag(ACFTYPE) "swiglal_array_viewout_" %str(ACFTYPE)

// The %swiglal_array_type() macro maps TYPEs of C arrays to an ACFTYPE of the
// appropriate array conversion functions, using a special typemap. The typemap also
// ensures that fragments containing the required conversion functions are included.
%define %_swiglal_array_type(ACFTYPE, TYPE)
%fragment("swiglal_array_frags_" %str(ACFTYPE), "header",
          fragment=%swiglal_array_copyin_frag(ACFTYPE),
          fragment=%swiglal_array_copyout_frag(ACFTYPE),
          fragment=%swiglal_array_viewout_frag(ACFTYPE)) {};
%typemap(swiglal_array_typeid, fragment="swiglal_array_frags_" %str(ACFTYPE)) TYPE* %str(ACFTYPE);
%typemap(swiglal_array_typeid, fragment="swiglal_array_frags_" %str(ACFTYPE)) const TYPE* %str(ACFTYPE);
%typemap(swiglal_array_typeid, fragment="swiglal_array_frags_" %str(ACFTYPE)) TYPE[ANY] %str(ACFTYPE);
%typemap(swiglal_array_typeid, fragment="swiglal_array_frags_" %str(ACFTYPE)) TYPE[ANY][ANY] %str(ACFTYPE);
%typemap(swiglal_array_typeid, fragment="swiglal_array_frags_" %str(ACFTYPE)) const TYPE[ANY] %str(ACFTYPE);
%typemap(swiglal_array_typeid, fragment="swiglal_array_frags_" %str(ACFTYPE)) const TYPE[ANY][ANY] %str(ACFTYPE);
%enddef
%define %swiglal_array_type(ACFTYPE, ...)
%swiglal_map_a(%_swiglal_array_type, ACFTYPE, __VA_ARGS__);
%enddef

// Map C array TYPEs to array conversion function ACFTYPEs.
%swiglal_array_type(SWIGTYPE, SWIGTYPE);
%swiglal_array_type(LALchar, char*);
%swiglal_array_type(int8_t, char, signed char, int8_t);
%swiglal_array_type(uint8_t, unsigned char, uint8_t);
%swiglal_array_type(int16_t, short, int16_t);
%swiglal_array_type(uint16_t, unsigned short, uint16_t);
%swiglal_array_type(int32_t, int, int32_t, enum SWIGTYPE);
%swiglal_array_type(uint32_t, unsigned int, uint32_t);
%swiglal_array_type(int64_t, long long, int64_t);
%swiglal_array_type(uint64_t, unsigned long long, uint64_t);
%swiglal_array_type(float, float);
%swiglal_array_type(double, double);
%swiglal_array_type(gsl_complex_float, gsl_complex_float);
%swiglal_array_type(gsl_complex, gsl_complex);
%swiglal_array_type(COMPLEX8, COMPLEX8);
%swiglal_array_type(COMPLEX16, COMPLEX16);

// On modern systems, 'long' could be either 32- or 64-bit...
%swiglal_array_type(int32or64_t, long);
%swiglal_array_type(uint32or64_t, unsigned long);
%define %swiglal_array_int32or64_frags(FRAG, FUNC, INT)
%fragment(FRAG(INT##32or64_t), "header") {
%#if LONG_MAX > INT_MAX
%#define FUNC(INT##32or64_t) FUNC(INT##64_t)
%#else
%#define FUNC(INT##32or64_t) FUNC(INT##32_t)
%#endif
}
%enddef
%swiglal_array_int32or64_frags(%swiglal_array_copyin_frag, %swiglal_array_copyin_func, int);
%swiglal_array_int32or64_frags(%swiglal_array_copyout_frag, %swiglal_array_copyout_func, int);
%swiglal_array_int32or64_frags(%swiglal_array_viewout_frag, %swiglal_array_viewout_func, int);
%swiglal_array_int32or64_frags(%swiglal_array_copyin_frag, %swiglal_array_copyin_func, uint);
%swiglal_array_int32or64_frags(%swiglal_array_copyout_frag, %swiglal_array_copyout_func, uint);
%swiglal_array_int32or64_frags(%swiglal_array_viewout_frag, %swiglal_array_viewout_func, uint);

// Call the appropriate conversion function for C TYPE arrays.
#define %swiglal_array_copyin(TYPE) %swiglal_array_copyin_func($typemap(swiglal_array_typeid, TYPE))
#define %swiglal_array_copyout(TYPE) %swiglal_array_copyout_func($typemap(swiglal_array_typeid, TYPE))
#define %swiglal_array_viewout(TYPE) %swiglal_array_viewout_func($typemap(swiglal_array_typeid, TYPE))

//
// Typemaps which convert to/from fixed-size arrays.
//

// Input typemaps for functions and structs:
%typemap(in, noblock=1) SWIGTYPE[ANY], SWIGTYPE INOUT[ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0];
  $1 = &temp$argnum[0];
  {
    const size_t dims[] = {$1_dim0};
    const size_t strides[] = {1};
    // swiglal_array_typeid input type: $1_type
    int ecode = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr($1),
                                                      sizeof($1[0]), 1, dims, strides,
                                                      false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                                      $disown | %convertptr_flags);
    if (!SWIG_IsOK(ecode)) {
      %argument_fail(ecode, "$type", $symname, $argnum);
    }
  }
}
%typemap(in, noblock=1) SWIGTYPE[ANY][ANY], SWIGTYPE INOUT[ANY][ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0][$1_dim1];
  $1 = &temp$argnum[0];
  {
    const size_t dims[] = {$1_dim0, $1_dim1};
    const size_t strides[] = {$1_dim1, 1};
    // swiglal_array_typeid input type: $1_type
    int ecode = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr($1),
                                               sizeof($1[0][0]), 2, dims, strides,
                                               false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                               $disown | %convertptr_flags);
    if (!SWIG_IsOK(ecode)) {
      %argument_fail(ecode, "$type", $symname, $argnum);
    }
  }
}

// Input typemaps for global variables:
%typemap(varin) SWIGTYPE[ANY] {
  const size_t dims[] = {$1_dim0};
  const size_t strides[] = {1};
  // swiglal_array_typeid input type: $1_type
  int ecode = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr($1),
                                             sizeof($1[0]), 1, dims, strides,
                                             false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                             %convertptr_flags);
  if (!SWIG_IsOK(ecode)) {
    %variable_fail(ecode, "$type", $symname);
  }
}
%typemap(varin) SWIGTYPE[ANY][ANY] {
  const size_t dims[] = {$1_dim0, $1_dim1};
  const size_t strides[] = {$1_dim1, 1};
  // swiglal_array_typeid input type: $1_type
  int ecode = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr($1),
                                             sizeof($1[0][0]), 2, dims, strides,
                                             false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                             %convertptr_flags);
  if (!SWIG_IsOK(ecode)) {
    %variable_fail(ecode, "$type", $symname);
  }
}

// Output typemaps for functions and structs:
%typemap(out) SWIGTYPE[ANY] {
  const size_t dims[] = {$1_dim0};
  const size_t strides[] = {1};
  // swiglal_array_typeid input type: $1_type
%#if $owner & SWIG_POINTER_OWN
  %set_output(%swiglal_array_copyout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                              sizeof($1[0]), 1, dims, strides,
                                              false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                              $owner | %newpointer_flags));
%#else
  %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                              sizeof($1[0]), 1, dims, strides,
                                              false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                              $owner | %newpointer_flags));
%#endif
}
%typemap(out) SWIGTYPE[ANY][ANY] {
  const size_t dims[] = {$1_dim0, $1_dim1};
  const size_t strides[] = {$1_dim1, 1};
  // swiglal_array_typeid input type: $1_type
%#if $owner & SWIG_POINTER_OWN
  %set_output(%swiglal_array_copyout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                              sizeof($1[0][0]), 2, dims, strides,
                                              false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                              $owner | %newpointer_flags));
%#else
  %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                              sizeof($1[0][0]), 2, dims, strides,
                                              false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                              $owner | %newpointer_flags));
%#endif
}

// Output typemaps for global variables:
%typemap(varout) SWIGTYPE[ANY] {
  const size_t dims[] = {$1_dim0};
  const size_t strides[] = {1};
  // swiglal_array_typeid input type: $1_type
  %set_output(%swiglal_array_viewout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                              sizeof($1[0]), 1, dims, strides,
                                              false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                              %newpointer_flags));
}
%typemap(varout) SWIGTYPE[ANY][ANY] {
  const size_t dims[] = {$1_dim0, $1_dim1};
  const size_t strides[] = {$1_dim1, 1};
  // swiglal_array_typeid input type: $1_type
  %set_output(%swiglal_array_viewout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                              sizeof($1[0][0]), 2, dims, strides,
                                              false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                              %newpointer_flags));
}

// Argument-output typemaps for functions:
%typemap(in, noblock=1, numinputs=0) SWIGTYPE OUTPUT[ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0];
  $1 = &temp$argnum[0];
}
%typemap(argout) SWIGTYPE OUTPUT[ANY], SWIGTYPE INOUT[ANY] {
  const size_t dims[] = {$1_dim0};
  const size_t strides[] = {1};
  // swiglal_array_typeid input type: $1_type
  %append_output(%swiglal_array_copyout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                                 sizeof($1[0]), 1, dims, strides,
                                                 false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                                 SWIG_POINTER_OWN | %newpointer_flags));
}
%typemap(in, noblock=1, numinputs=0) SWIGTYPE OUTPUT[ANY][ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0][$1_dim1];
  $1 = &temp$argnum[0];
}
%typemap(argout) SWIGTYPE OUTPUT[ANY][ANY], SWIGTYPE INOUT[ANY][ANY] {
  const size_t dims[] = {$1_dim0, $1_dim1};
  const size_t strides[] = {$1_dim1, 1};
  // swiglal_array_typeid input type: $1_type
  %append_output(%swiglal_array_copyout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                                 sizeof($1[0][0]), 2, dims, strides,
                                                 false, $typemap(swiglal_fixarr_tinfo, $1_type),
                                                 SWIG_POINTER_OWN | %newpointer_flags));
}

// Public macros to make fixed nD arrays:
// * output-only arguments: SWIGLAL(OUTPUT_ARRAY_nD(TYPE, ...))
%define %swiglal_public_OUTPUT_ARRAY_1D(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE OUTPUT[ANY], TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_OUTPUT_ARRAY_1D(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_OUTPUT_ARRAY_2D(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE OUTPUT[ANY][ANY], TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_OUTPUT_ARRAY_2D(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
// * input-output arguments: SWIGLAL(INOUT_ARRAY_nD(TYPE, ...))
%define %swiglal_public_INOUT_ARRAY_1D(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE INOUT[ANY], TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_INOUT_ARRAY_1D(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_INOUT_ARRAY_2D(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE INOUT[ANY][ANY], TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_INOUT_ARRAY_2D(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef

// Get the correct descriptor for a dynamic array element:
// Always return a pointer description, even for non-pointer types, and
// determine whether array is an array of pointers or of data blocks.
%typemap(swiglal_dynarr_isptr) SWIGTYPE  "false";
%typemap(swiglal_dynarr_tinfo) SWIGTYPE  "$&descriptor";
%typemap(swiglal_dynarr_isptr) SWIGTYPE* "true";
%typemap(swiglal_dynarr_tinfo) SWIGTYPE* "$descriptor";

// Create immutable members for accessing the array's dimensions.
// NI is the name of the dimension member, and SIZET is its type.
%define %swiglal_array_dynamic_size(SIZET, NI)
  %feature("action") NI {
    result = %static_cast(arg1->NI, SIZET);
  }
  %extend {
    const SIZET NI;
  }
  %feature("action", "") NI;
%enddef // %swiglal_array_dynamic_size()

// Check that array strides are non-zero, otherwise fail.
%define %swiglal_array_dynamic_check_strides(NAME, DATA, I)
  if (strides[I-1] == 0) {
    SWIG_exception_fail(SWIG_IndexError, "Stride of dimension "#I" of "#NAME"."#DATA" is zero");
  }
%enddef // %swiglal_array_dynamic_check_strides()

// The %swiglal_array_dynamic_<n>D() macros create typemaps which convert
// <n>-D dynamically-allocated arrays in structs NAME. The macros must be
// added inside the definition of the struct, before the struct members
// comprising the array are defined. The DATA and N{I,J} members give
// the array data and dimensions, TYPE and SIZET give their respective
// types. The S{I,J} give the strides of the array, in number of elements.
// If the sizes or strides are members of the struct, 'arg1->' should be
// used to access the struct itself.
// 1-D arrays:
%define %swiglal_array_dynamic_1D(NAME, TYPE, SIZET, DATA, NI, SI)

  // Typemaps which convert to/from the dynamically-allocated array.
  %typemap(in, noblock=1) TYPE* DATA {
    if (arg1) {
      const size_t dims[] = {NI};
      const size_t strides[] = {SI};
      %swiglal_array_dynamic_check_strides(NAME, DATA, 1);
      $1 = %reinterpret_cast(arg1->DATA, TYPE*);
      // swiglal_array_typeid input type: $1_type
      int ecode = %swiglal_array_copyin($1_type)(swiglal_self(), $input, %as_voidptr($1),
                                                 sizeof(TYPE), 1, dims, strides,
                                                 $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                 $disown | %convertptr_flags);
      if (!SWIG_IsOK(ecode)) {
        %argument_fail(ecode, "$type", $symname, $argnum);
      }
    }
  }
  %typemap(out, noblock=1) TYPE* DATA {
    if (arg1) {
      const size_t dims[] = {NI};
      const size_t strides[] = {SI};
      %swiglal_array_dynamic_check_strides(NAME, DATA, 1);
      $1 = %reinterpret_cast(arg1->DATA, TYPE*);
      // swiglal_array_typeid input type: $1_type
      %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                                  sizeof(TYPE), 1, dims, strides,
                                                  $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                  $owner | %newpointer_flags));
    }
  }

  // Clear unneeded typemaps and features.
  %typemap(memberin, noblock=1) TYPE* DATA "";
  %typemap(argout, noblock=1) TYPE* DATA "";
  %typemap(freearg, noblock=1) TYPE* DATA "";
  %feature("action") DATA "";
  %feature("except") DATA "";

  // Create the array's data member.
  %extend {
    TYPE *DATA;
  }

  // Restore modified typemaps and features.
  %feature("action", "") DATA;
  %feature("except", "") DATA;
  %clear TYPE* DATA;

%enddef // %swiglal_array_dynamic_1D()
// a 1-D array of pointers to 1-D arrays:
// * NI refer to the array of pointers, NJ/SJ refer to the arrays pointed to
%define %swiglal_array_dynamic_1d_ptr_1d(NAME, TYPE, SIZET, DATA, NI, NJ, SJ)

  // Typemaps which convert from the dynamically-allocated array, indexed by 'arg2'
  %typemap(out, noblock=1) TYPE* DATA {
    if (arg1) {
      const size_t dims[] = {NJ};
      const size_t strides[] = {SJ};
      %swiglal_array_dynamic_check_strides(NAME, DATA, 1);
      if (!(0 <= arg2 && arg2 < NI)) {
        SWIG_exception_fail(SWIG_IndexError, "Index to "#NAME"."#DATA" is outside of range [0,"#NI"]");
      }
      $1 = %reinterpret_cast(arg1->DATA, TYPE*)[arg2];
      // swiglal_array_typeid input type: $1_type
      %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                                  sizeof(TYPE), 1, dims, strides,
                                                  $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                  $owner | %newpointer_flags));
    }
  }

  // Clear unneeded typemaps and features.
  %typemap(memberin, noblock=1) TYPE* DATA "";
  %typemap(argout, noblock=1) TYPE* DATA "";
  %typemap(freearg, noblock=1) TYPE* DATA "";
  %feature("action") DATA "";
  %feature("except") DATA "";

  // Create the array's indexed data member.
  %extend {
    TYPE *DATA(const SIZET index);
  }

  // Restore modified typemaps and features.
  %feature("action", "") DATA;
  %feature("except", "") DATA;
  %clear TYPE* DATA;

%enddef // %swiglal_array_dynamic_1d_ptr_1d()
// 2-D arrays:
%define %swiglal_array_dynamic_2D(NAME, TYPE, SIZET, DATA, NI, NJ, SI, SJ)

  // Typemaps which convert to/from the dynamically-allocated array.
  %typemap(in, noblock=1) TYPE* DATA {
    if (arg1) {
      const size_t dims[] = {NI, NJ};
      const size_t strides[] = {SI, SJ};
      %swiglal_array_dynamic_check_strides(NAME, DATA, 1);
      %swiglal_array_dynamic_check_strides(NAME, DATA, 2);
      $1 = %reinterpret_cast(arg1->DATA, TYPE*);
      // swiglal_array_typeid input type: $1_type
      int ecode = %swiglal_array_copyin($1_type)(swiglal_self(), $input, %as_voidptr($1),
                                                 sizeof(TYPE), 2, dims, strides,
                                                 $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                 $disown | %convertptr_flags);
      if (!SWIG_IsOK(ecode)) {
        %argument_fail(ecode, "$type", $symname, $argnum);
      }
    }
  }
  %typemap(out, noblock=1) TYPE* DATA {
    if (arg1) {
      const size_t dims[] = {NI, NJ};
      const size_t strides[] = {SI, SJ};
      %swiglal_array_dynamic_check_strides(NAME, DATA, 1);
      %swiglal_array_dynamic_check_strides(NAME, DATA, 2);
      $1 = %reinterpret_cast(arg1->DATA, TYPE*);
      // swiglal_array_typeid input type: $1_type
      %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                                  sizeof(TYPE), 2, dims, strides,
                                                  $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                  $owner | %newpointer_flags));
    }
  }

  // Clear unneeded typemaps and features.
  %typemap(memberin, noblock=1) TYPE* DATA "";
  %typemap(argout, noblock=1) TYPE* DATA "";
  %typemap(freearg, noblock=1) TYPE* DATA "";
  %feature("action") DATA "";
  %feature("except") DATA "";

  // Create the array's data member.
  %extend {
    TYPE *DATA;
  }

  // Restore modified typemaps and features.
  %feature("action", "") DATA;
  %feature("except", "") DATA;
  %clear TYPE* DATA;

%enddef // %swiglal_array_dynamic_2D()

// These macros should be called from within the definitions of
// LAL structs NAME containing dynamically-allocated arrays.
// 1-D arrays, e.g:
//   SIZET NI;
//   TYPE* DATA;
%define %swiglal_public_ARRAY_1D(NAME, TYPE, DATA, SIZET, NI)
  %swiglal_array_dynamic_size(SIZET, NI);
  %swiglal_array_dynamic_1D(NAME, TYPE, SIZET, DATA, arg1->NI, 1);
  %ignore DATA;
  %ignore NI;
%enddef
#define %swiglal_public_clear_ARRAY_1D(NAME, TYPE, DATA, SIZET, NI)
// a 1-D array of pointers to 1-D arrays, e.g.:
//   SIZET NI;
//   SIZET NJ;
//   ATYPE* DATA[(NI members)];
%define %swiglal_public_ARRAY_1D_PTR_1D(NAME, TYPE, DATA, SIZET, NI, NJ)
  %swiglal_array_dynamic_size(SIZET, NI);
  %swiglal_array_dynamic_size(SIZET, NJ);
  %swiglal_array_dynamic_1d_ptr_1d(NAME, TYPE, SIZET, DATA, arg1->NI, arg1->NJ, 1);
  %ignore DATA;
  %ignore NJ;
%enddef
#define %swiglal_public_clear_ARRAY_1D_PTR_1D(NAME, TYPE, DATA, SIZET, NI, NJ)
// 2-D arrays of fixed-length arrays, e.g:
//   typedef ETYPE[NJ] ATYPE;
//   SIZET NI;
//   ATYPE* DATA;
%define %swiglal_public_ARRAY_2D_FIXED(NAME, ETYPE, ATYPE, DATA, SIZET, NI)
  %swiglal_array_dynamic_size(SIZET, NI);
  %swiglal_array_dynamic_2D(NAME, ETYPE, SIZET, DATA, arg1->NI, (sizeof(ATYPE)/sizeof(ETYPE)), (sizeof(ATYPE)/sizeof(ETYPE)), 1);
  %ignore DATA;
  %ignore NI;
%enddef
#define %swiglal_public_clear_ARRAY_2D_FIXED(NAME, ETYPE, ATYPE, DATA, SIZET, NI)
// 2-D arrays, e.g:
//   SIZET NI, NJ;
//   TYPE* DATA;
%define %swiglal_public_ARRAY_2D(NAME, TYPE, DATA, SIZET, NI, NJ)
  %swiglal_array_dynamic_size(SIZET, NI);
  %swiglal_array_dynamic_size(SIZET, NJ);
  %swiglal_array_dynamic_2D(NAME, TYPE, SIZET, DATA, arg1->NI, arg1->NJ, arg1->NJ, 1);
  %ignore DATA;
  %ignore NI;
  %ignore NJ;
%enddef
#define %swiglal_public_clear_ARRAY_2D(NAME, TYPE, DATA, SIZET, NI, NJ)

// If multiple arrays in the same struct use the same length parameter, this
// directive is required before the struct definition to suppress warnings.
%define %swiglal_public_ARRAY_MULTIPLE_LENGTHS(TAGNAME, ...)
%swiglal_map_ab(%swiglal_warnfilter_nspace, SWIGWARN_PARSE_REDEFINED, TAGNAME, __VA_ARGS__);
%enddef
#define %swiglal_public_clear_ARRAY_MULTIPLE_LENGTHS(TAGNAME, ...)

////////// Include scripting-language-specific interface headers //////////

#ifdef SWIGOCTAVE
%include <lal/swiglal_octave.i>
#endif
#ifdef SWIGPYTHON
%include <lal/swiglal_python.i>
#endif

////////// General typemaps and macros //////////

// The SWIGLAL(RETURN_VOID(TYPE,...)) public macro can be used to ensure
// that the return value of a function is always ignored.
%define %swiglal_public_RETURN_VOID(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE SWIGLAL_RETURN_VOID, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_RETURN_VOID(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
%typemap(out, noblock=1) SWIGTYPE SWIGLAL_RETURN_VOID {
  %set_output(VOID_Object);
}

// The SWIGLAL(RETURN_VALUE(TYPE,...)) public macro can be used to ensure
// that the return value of a function is not ignored, if the return value
// has previously been ignored in the generated wrappings.
%define %swiglal_public_RETURN_VALUE(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
#define %swiglal_public_clear_RETURN_VALUE(TYPE, ...)

// The SWIGLAL(DISABLE_EXCEPTIONS(...)) public macro is useful for
// functions which manipulate XLAL error codes, which thus require
// XLAL exception handling to be disabled.
%define %swiglal_public_DISABLE_EXCEPTIONS(...)
%swiglal_map_ab(%swiglal_feature, "except", "$action", __VA_ARGS__);
%enddef
#define %swiglal_public_clear_DISABLE_EXCEPTIONS(...)

// The SWIGLAL(FUNCTION_POINTER(...)) macro can be used to create
// a function pointer constant, for functions which need to be used
// as callback functions.
%define %swiglal_public_FUNCTION_POINTER(...)
%swiglal_map_ab(%swiglal_feature, "callback", "%sPtr", __VA_ARGS__);
%enddef
#define %swiglal_public_clear_FUNCTION_POINTER(...)

// The SWIGLAL(IMMUTABLE_MEMBERS(TAGNAME, ...)) macro can be used to make
// the listed members of the struct TAGNAME immutable.
%define %swiglal_public_IMMUTABLE_MEMBERS(TAGNAME, ...)
%swiglal_map_abc(%swiglal_feature_nspace, "immutable", "1", TAGNAME, __VA_ARGS__);
%enddef
#define %swiglal_public_clear_IMMUTABLE_MEMBERS(...)

// Typemap for functions which return 'int'. If these functions also return
// other output arguments (via 'argout' typemaps), the 'int' return value is
// ignored. This is because 'int' is very commonly used to return an XLAL
// error code, which will be converted into a native scripting-language
// exception, and so the error code itself is not needed directly. To avoid
// having to unpack the error code when collecting the other output arguments,
// therefore, it is ignored in the wrappings. Functions which fit this criteria
// but do return a useful 'int' can use SWIGLAL(RETURN_VALUE(int, ...)) to
// disable this behaviour.
//
// For functions, since %feature("new") is set, the 'out' typemap will have $owner=1,
// and the 'newfree' typemap is also applied. The 'out' typemap ignores the 'int'
// return value by setting the output argument list to VOID_Object; the wrapping
// function them proceeds to add other output arguments to the list, if any. After
// this, the 'newfree' typemap is triggered, which appends the 'int' return if the
// output argument list is empty, using the scripting-language-specific macro
// swiglal_append_output_if_empty(). For structs, $owner=0, so the int return is
// set straight away, and the 'newfree' typemap is never applied.
%typemap(out, noblock=1, fragment=SWIG_From_frag(int)) int {
%#if $owner
  %set_output(VOID_Object);
%#else
  %set_output(SWIG_From(int)($1));
%#endif
}
%typemap(newfree, noblock=1, fragment=SWIG_From_frag(int)) int {
  swiglal_append_output_if_empty(SWIG_From(int)($1));
}

// Typemaps for empty arguments. These typemaps are useful when no input from the
// scripting language is required, and an empty struct needs to be supplied to
// the C function. The SWIGLAL(EMPTY_ARGUMENT(TYPE, ...)) macro applies the typemap which
// supplies a static struct, while the SWIGLAL(NEW_EMPTY_ARGUMENT(TYPE, ...)) macro
// applies the typemap which supplies a dynamically-allocated struct. These typemaps
// may cause there to be no SWIG-wrapped object for the first argument; if so,
// 'swiglal_no_1starg' is defined for the duration of the wrapping function.
%define %swiglal_public_EMPTY_ARGUMENT(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE* SWIGLAL_EMPTY_ARGUMENT, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_EMPTY_ARGUMENT(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
%typemap(in, noblock=1, numinputs=0) SWIGTYPE* SWIGLAL_EMPTY_ARGUMENT ($*ltype emptyarg) {
  memset(&emptyarg, 0, sizeof($*type));
  $1 = &emptyarg;
%#if $argnum == 1
%#define swiglal_no_1starg
%#endif
}
%typemap(freearg, noblock=1) SWIGTYPE* SWIGLAL_EMPTY_ARGUMENT {
%#undef swiglal_no_1starg
}
%define %swiglal_public_NEW_EMPTY_ARGUMENT(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE* SWIGLAL_NEW_EMPTY_ARGUMENT, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_NEW_EMPTY_ARGUMENT(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
%typemap(in, noblock=1, numinputs=0) SWIGTYPE* SWIGLAL_NEW_EMPTY_ARGUMENT {
  $1 = %swiglal_new_instance($*type);
%#if $argnum == 1
%#define swiglal_no_1starg
%#endif
}
%typemap(freearg, noblock=1) SWIGTYPE* SWIGLAL_NEW_EMPTY_ARGUMENT {
%#undef swiglal_no_1starg
}

// SWIG conversion functions for C99 integer types.
// These are mapped to the corresponding basic C types,
// conversion functions for which are supplied by SWIG.
%define %swiglal_numeric_typedef(CHECKCODE, BASETYPE, TYPE)
%numeric_type_from(TYPE, BASETYPE);
%numeric_type_asval(TYPE, BASETYPE, "swiglal_empty_frag", false);
%typemaps_primitive(%checkcode(CHECKCODE), TYPE);
%enddef
%swiglal_numeric_typedef(INT8, signed char, int8_t);
%swiglal_numeric_typedef(UINT8, unsigned char, uint8_t);
%swiglal_numeric_typedef(INT16, short, int16_t);
%swiglal_numeric_typedef(UINT16, unsigned short, uint16_t);
%swiglal_numeric_typedef(INT32, int, int32_t);
%swiglal_numeric_typedef(UINT32, unsigned int, uint32_t);
%swiglal_numeric_typedef(INT64, long long, int64_t);
%swiglal_numeric_typedef(UINT64, unsigned long long, uint64_t);

// Fragments and typemaps for the LAL BOOLEAN type. The fragments re-use
// existing scriping-language conversion functions for the C/C++ boolean type.
// Appropriate typemaps are then generated by %typemaps_asvalfromn().
%fragment(SWIG_From_frag(BOOLEAN), "header", fragment=SWIG_From_frag(bool)) {
  SWIGINTERNINLINE SWIG_Object SWIG_From_dec(BOOLEAN)(BOOLEAN value) {
    return SWIG_From(bool)(value ? true : false);
  }
}
%fragment(SWIG_AsVal_frag(BOOLEAN), "header", fragment=SWIG_AsVal_frag(bool)) {
  SWIGINTERN int SWIG_AsVal_dec(BOOLEAN)(SWIG_Object obj, BOOLEAN *val) {
    bool v;
    int ecode = SWIG_AsVal(bool)(obj, val ? &v : 0);
    if (!SWIG_IsOK(ecode)) {
      return SWIG_TypeError;
    }
    if (val) {
      *val = v ? 1 : 0;
    }
    return ecode;
  }
}
%typemaps_primitive(%checkcode(BOOL), BOOLEAN);

// Fragments and typemaps for LAL strings, which should be (de)allocated
// using LAL memory functions. The fragments re-use existing scriping-language
// conversion functions for ordinary char* strings. Appropriate typemaps are
// then generated by %typemaps_string_alloc(), with custom memory allocators.
%fragment("SWIG_FromLALcharPtrAndSize", "header", fragment="SWIG_FromCharPtrAndSize") {
  SWIGINTERNINLINE SWIG_Object SWIG_FromLALcharPtrAndSize(const char *str, size_t size) {
    return SWIG_FromCharPtrAndSize(str, size);
  }
}
%fragment("SWIG_AsLALcharPtrAndSize", "header", fragment="SWIG_AsCharPtrAndSize") {
  SWIGINTERN int SWIG_AsLALcharPtrAndSize(SWIG_Object obj, char **pstr, size_t *psize, int *alloc) {
    char *slstr = 0;
    size_t slsize = 0;
    int slalloc = 0;
    // Get pointer to scripting-language string 'slstr' and size 'slsize'.
    // The 'slalloc' argument indicates whether a new string was allocated.
    int ecode = SWIG_AsCharPtrAndSize(obj, &slstr, &slsize, &slalloc);
    if (!SWIG_IsOK(ecode)) {
      return SWIG_TypeError;
    }
    // Return the string, if needed.
    if (pstr) {
      // Free the LAL string if it is already allocated.
      if (*pstr) {
        XLALFree(*pstr);
      }
      if (alloc) {
        // Copy the scripting-language string into a LAL-managed memory string.
        *pstr = %swiglal_new_copy_array(slstr, slsize, char);
        *alloc = SWIG_NEWOBJ;
      }
      else {
        return SWIG_TypeError;
      }
    }
    // Return the size (length+1) of the string, if needed.
    if (psize) {
      *psize = slsize;
    }
    // Free the scripting-language string, if it was allocated.
    if (slalloc == SWIG_NEWOBJ) {
      %delete_array(slstr);
    }
    return ecode;
  }
}
%fragment("SWIG_AsNewLALcharPtr", "header", fragment="SWIG_AsLALcharPtrAndSize") {
  SWIGINTERN int SWIG_AsNewLALcharPtr(SWIG_Object obj, char **pstr) {
    int alloc = 0;
    return SWIG_AsLALcharPtrAndSize(obj, pstr, 0, &alloc);
  }
}
#if SWIG_VERSION >= 0x030000
%typemaps_string_alloc(%checkcode(STRING), %checkcode(char), char, LALchar,
                       SWIG_AsLALcharPtrAndSize, SWIG_FromLALcharPtrAndSize,
                       strlen, SWIG_strnlen, %swiglal_new_copy_array, XLALFree,
                       "<limits.h>", CHAR_MIN, CHAR_MAX);
#else
%typemaps_string_alloc(%checkcode(STRING), %checkcode(char), char, LALchar,
                       SWIG_AsLALcharPtrAndSize, SWIG_FromLALcharPtrAndSize,
                       strlen, %swiglal_new_copy_array, XLALFree,
                       "<limits.h>", CHAR_MIN, CHAR_MAX);
#endif

// Typemaps for string pointers.  By default, treat arguments of type char**
// as output-only arguments, which do not require a scripting-language input
// argument, and return their results in the output argument list. Also
// supply an INOUT typemap for input-output arguments, which allows a
// scripting-language input string to be supplied. The INOUT typemaps can be
// applied as needed using the SWIGLAL(INOUT_STRINGS(...)) macro.
%typemap(in, noblock=1, numinputs=0) char ** (char *str = NULL, int alloc = 0) {
  $1 = %reinterpret_cast(&str, $ltype);
  alloc = 0;
}
%typemap(in, noblock=1, fragment="SWIG_AsLALcharPtrAndSize") char ** INOUT (char *str = NULL, int alloc = 0, int ecode = 0) {
  ecode = SWIG_AsLALcharPtr($input, &str, &alloc);
  if (!SWIG_IsOK(ecode)) {
    %argument_fail(ecode, "$type", $symname, $argnum);
  }
  $1 = %reinterpret_cast(&str, $ltype);
}
%typemap(argout, noblock=1) char ** {
  %append_output(SWIG_FromLALcharPtr(str$argnum));
}
%typemap(freearg, match="in") char ** {
  if (SWIG_IsNewObj(alloc$argnum)) {
    XLALFree(str$argnum);
  }
}
%define %swiglal_public_INOUT_STRINGS(...)
%swiglal_map_ab(%swiglal_apply, char ** INOUT, char **, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_INOUT_STRINGS(...)
%swiglal_map_a(%swiglal_clear, char **, __VA_ARGS__);
%enddef

// Do not try to free const char* return arguments.
%typemap(newfree,noblock=1) const char* "";

// Input typemap for pointer-to-const SWIGTYPEs. This typemap is identical to the
// standard SWIGTYPE pointer typemap, except $disown is commented out. This prevents
// SWIG transferring ownership of SWIG-wrapped objects when assigning to pointer-to-const
// members of structs. In this case, it is assumed that the struct does not want to take
// ownership of the pointer, since it cannot free it (since it is a pointer-to-const).
%typemap(in,noblock=1) const SWIGTYPE * (void *argp = 0, int res = 0) {
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  $1 = %reinterpret_cast(argp, $ltype);
}
%typemap(freearg) const SWIGTYPE * "";

// Typemaps for output SWIGTYPEs. This typemaps will match either the SWIG-wrapped
// return argument from functions (which will have the SWIG_POINTER_OWN bit set
// in $owner) or return a member of a struct through a 'get' functions (in which
// case SWIG_POINTER_OWN will not be set). They require the following macros:
//
// The macro %swiglal_store_parent() is called to store a reference to the struct
// containing the member being accessed, in order to prevent it from being destroyed
// as long as the SWIG-wrapped member object is in scope. The return object is then
// always created with SWIG_POINTER_OWN, so that its destructor will always be called.
%define %swiglal_store_parent(PTR, OWNER, SELF)
%#if !(OWNER & SWIG_POINTER_OWN)
  if (%as_voidptr(PTR) != NULL) {
    swiglal_store_parent(%as_voidptr(PTR), SELF);
  }
%#endif
%enddef
//
// The macro %swiglal_set_output() sets the output of the wrapping function. If the
// (pointer) return type of the function is the same as its first argument, then
// 'swiglal_return_1starg_##NAME' is defined. Unless 'swiglal_no_1starg' is defined
// (in which case the first argument is being handled by e.g. the EMPTY_ARGUMENT
// typemap), the macro compares the pointer of the return value (result) to that of
// the first argument (arg1). If they're equal, the SWIG-wrapped function will return
// a *reference* to the SWIG object wrapping the first argument, i.e. the same object
// with its reference count incremented. That way, if the return value is assigned to
// a different scripting-language variable than the first argument, the underlying
// C struct will not be destroyed until both scripting-language variables are cleared.
// If the pointers are not equal, or one pointer is NULL, the macro return a SWIG
// object wrapping the new C struct.
%define %swiglal_set_output(NAME, OBJ)
%#if defined(swiglal_return_1starg_##NAME) && !defined(swiglal_no_1starg)
  if (result != NULL && result == arg1) {
    %set_output(swiglal_get_reference(swiglal_1starg()));
  } else {
    %set_output(OBJ);
  }
%#else
  %set_output(OBJ);
%#endif
%enddef
//
// Typemaps:
%typemap(out,noblock=1) SWIGTYPE *, SWIGTYPE &, SWIGTYPE[] {
  %swiglal_store_parent($1, $owner, swiglal_self());
  %swiglal_set_output($1_name, SWIG_NewPointerObj(%as_voidptr($1), $descriptor, ($owner | %newpointer_flags) | SWIG_POINTER_OWN));
}
%typemap(out,noblock=1) const SWIGTYPE *, const SWIGTYPE &, const SWIGTYPE[] {
  %set_output(SWIG_NewPointerObj(%as_voidptr($1), $descriptor, ($owner | %newpointer_flags) & ~SWIG_POINTER_OWN));
}
%typemap(out, noblock=1) SWIGTYPE *const& {
  %swiglal_store_parent(*$1, $owner, swiglal_self());
  %swiglal_set_output($1_name, SWIG_NewPointerObj(%as_voidptr(*$1), $*descriptor, ($owner | %newpointer_flags) | SWIG_POINTER_OWN));
}
%typemap(out, noblock=1) const SWIGTYPE *const& {
  %set_output(SWIG_NewPointerObj(%as_voidptr(*$1), $*descriptor, ($owner | %newpointer_flags) & ~SWIG_POINTER_OWN));
}
%typemap(out, noblock=1) SWIGTYPE (void* copy = NULL) {
  copy = %swiglal_new_copy($1, $ltype);
  %swiglal_store_parent(copy, SWIG_POINTER_OWN, swiglal_self());
  %set_output(SWIG_NewPointerObj(copy, $&descriptor, (%newpointer_flags) | SWIG_POINTER_OWN));
}
%typemap(varout, noblock=1) SWIGTYPE *, SWIGTYPE [] {
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr($1), $descriptor, (%newpointer_flags) & ~SWIG_POINTER_OWN));
}
%typemap(varout, noblock=1) SWIGTYPE & {
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(&$1), $descriptor, (%newpointer_flags) & ~SWIG_POINTER_OWN));
}
%typemap(varout, noblock=1) SWIGTYPE {
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(&$1), $&descriptor, (%newpointer_flags) & ~SWIG_POINTER_OWN));
}

// The SWIGLAL(RETURNS_PROPERTY(...)) macro is used when a function returns an object whose
// memory is owned by the object supplied as the first argument to the function.
// Typically this occurs when the function is returning some property of its first
// argument. The macro applies a typemap which calles swiglal_store_parent() to store
// a reference to the first argument as the 'parent' of the return argument, so that
// the parent will not be destroyed as long as the return value is in scope.
%define %swiglal_public_RETURNS_PROPERTY(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE* SWIGLAL_RETURNS_PROPERTY, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_RETURNS_PROPERTY(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
%typemap(out,noblock=1) SWIGTYPE* SWIGLAL_RETURNS_PROPERTY {
%#ifndef swiglal_no_1starg
  %swiglal_store_parent($1, 0, swiglal_1starg());
%#endif
  %set_output(SWIG_NewPointerObj(%as_voidptr($1), $descriptor, ($owner | %newpointer_flags) | SWIG_POINTER_OWN));
}

// The SWIGLAL(ACQUIRES_OWNERSHIP(...)) macro indicates that a function will acquire ownership
// of a particular argument, e.g. by storing that argument in some container, and that therefore
// the SWIG object wrapping that argument should no longer own its memory.
%define %swiglal_public_ACQUIRES_OWNERSHIP(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE* DISOWN, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_ACQUIRES_OWNERSHIP(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef

// Typemaps for pointers to primitive scalars. These are treated as output-only
// arguments by default, by globally applying the SWIG OUTPUT typemaps.
%apply int* OUTPUT { enum SWIGTYPE* };
%apply short* OUTPUT { short* };
%apply unsigned short* OUTPUT { unsigned short* };
%apply int* OUTPUT { int* };
%apply unsigned int* OUTPUT { unsigned int* };
%apply long* OUTPUT { long* };
%apply unsigned long* OUTPUT { unsigned long* };
%apply long long* OUTPUT { long long* };
%apply unsigned long long* OUTPUT { unsigned long long* };
%apply float* OUTPUT { float* };
%apply double* OUTPUT { double* };
%apply int8_t* OUTPUT { int8_t* };
%apply uint8_t* OUTPUT { uint8_t* };
%apply int16_t* OUTPUT { int16_t* };
%apply uint16_t* OUTPUT { uint16_t* };
%apply int32_t* OUTPUT { int32_t* };
%apply uint32_t* OUTPUT { uint32_t* };
%apply int64_t* OUTPUT { int64_t* };
%apply uint64_t* OUTPUT { uint64_t* };
%apply gsl_complex_float* OUTPUT { gsl_complex_float* };
%apply gsl_complex* OUTPUT { gsl_complex* };
%apply COMPLEX8* OUTPUT { COMPLEX8* };
%apply COMPLEX16* OUTPUT { COMPLEX16* };
// The INOUT typemaps can be applied instead using the macro SWIGLAL(INOUT_SCALARS(TYPE, ...)).
%define %swiglal_public_INOUT_SCALARS(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, TYPE INOUT, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_INOUT_SCALARS(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
// The INPUT typemaps can be applied instead using the macro SWIGLAL(INPUT_SCALARS(TYPE, ...)).
%typemap(argout, noblock=1) SWIGTYPE* SWIGLAL_INPUT_SCALAR "";
%define %swiglal_public_INPUT_SCALARS(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, TYPE INPUT, TYPE, __VA_ARGS__);
%swiglal_map_ab(%swiglal_apply, SWIGTYPE* SWIGLAL_INPUT_SCALAR, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_INPUT_SCALARS(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef

// Typemaps for double pointers. By default, treat arguments of type TYPE**
// as output-only arguments, which do not require a scripting-language input
// argument, and return their results in the output argument list. Also supply
// an INOUT typemap for input-output arguments, which allows a scripting-language
// input argument to be supplied. The INOUT typemaps can be applied as needed
// using the SWIGLAL(INOUT_STRUCTS(TYPE, ...)) macro.
%typemap(in, noblock=1, numinputs=0) SWIGTYPE ** (void *argp = NULL, int owner = 0) {
  $1 = %reinterpret_cast(&argp, $ltype);
  owner = (argp == NULL) ? SWIG_POINTER_OWN : 0;
}
%typemap(in, noblock=1) SWIGTYPE ** INOUT (void  *argp = NULL, int owner = 0, int ecode = 0) {
  ecode = SWIG_ConvertPtr($input, &argp, $*descriptor, ($disown | %convertptr_flags) | SWIG_POINTER_DISOWN);
  if (!SWIG_IsOK(ecode)) {
    %argument_fail(ecode, "$type", $symname, $argnum);
  }
  $1 = %reinterpret_cast(&argp, $ltype);
  owner = (argp == NULL) ? SWIG_POINTER_OWN : 0;
}
%typemap(argout, noblock=1) SWIGTYPE ** {
  %append_output(SWIG_NewPointerObj(%as_voidptr(*$1), $*descriptor, owner$argnum | %newpointer_flags));
}
%typemap(freearg) SWIGTYPE ** "";
%define %swiglal_public_INOUT_STRUCTS(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE ** INOUT, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_INOUT_STRUCTS(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef

// Make the wrapping of printf-style LAL functions a little safer, as suggested in
// the SWIG 2.0 documentation (section 13.5). These functions should now be safely
// able to print any string, so long as the format string is named "format" or "fmt".
%typemap(in, fragment="SWIG_AsLALcharPtrAndSize") (const char *SWIGLAL_PRINTF_FORMAT, ...)
(char fmt[] = "%s", char *str = 0, int alloc = 0)
{
  $1 = fmt;
  int ecode = SWIG_AsLALcharPtr($input, &str, &alloc);
  if (!SWIG_IsOK(ecode)) {
    %argument_fail(ecode, "$type", $symname, $argnum);
  }
  $2 = (void *) str;
}
%typemap(freearg, match="in") (const char *format, ...) {
  if (SWIG_IsNewObj(alloc$argnum)) {
    XLALFree(str$argnum);
  }
}
%apply (const char *SWIGLAL_PRINTF_FORMAT, ...) {
  (const char *format, ...), (const char *fmt, ...)
};

// Specialised input typemap for C file pointers. Generally it is not possible to
// convert scripting-language file objects into FILE*, since the scripting language
// may not provide access to the FILE*, or even be using FILE* internally for I/O.
// The FILE* will therefore have to be supplied from another SWIG-wrapped C function.
// For convenience, however, we allow the user to pass integers 0, 1, or 2 in place
// of a FILE*, as an instruction to use standard input, output, or error respectively.
%typemap(in, noblock=1, fragment=SWIG_AsVal_frag(int)) FILE* (void *argp = 0, int res = 0) {
  res = SWIG_ConvertPtr($input, &argp, $descriptor, $disown | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    int val = 0;
    res = SWIG_AsVal(int)($input, &val);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    } else {
      switch (val) {
      case 0:
        $1 = stdin;
        break;
      case 1:
        $1 = stdout;
        break;
      case 2:
        $1 = stderr;
        break;
      default:
        %argument_fail(SWIG_ValueError, "$type", $symname, $argnum);
      }
    }
  } else {
    $1 = %reinterpret_cast(argp, $ltype);
  }
}
%typemap(freearg) FILE* "";

// This macro supports functions which require a variable-length
// list of arguments of type TYPE, i.e. a list of strings. It
// generates SWIG "compact" default arguments, i.e. only one
// wrapping function where all missing arguments are assigned ENDVALUE,
// generates 11 additional optional arguments of type TYPE, and
// creates a contract ensuring that the last argument is always
// ENDVALUE, so that the argument list is terminated by ENDVALUE.
%define %swiglal_public_VARIABLE_ARGUMENT_LIST(FUNCTION, TYPE, ENDVALUE)
%feature("kwargs", 0) FUNCTION;
%feature("compactdefaultargs") FUNCTION;
%varargs(11, TYPE arg = ENDVALUE) FUNCTION;
%contract FUNCTION {
require:
  arg11 == ENDVALUE;
}
%enddef
#define %swiglal_public_clear_VARIABLE_ARGUMENT_LIST(FUNCTION, TYPE, ENDVALUE)

// This macro can be used to support structs which are not
// declared in LALSuite. It treats the struct as opaque,
// and attaches a destructor function to it.
%define %swiglal_public_EXTERNAL_STRUCT(NAME, DTORFUNC)
typedef struct {} NAME;
%ignore DTORFUNC;
%extend NAME {
  ~NAME() {
    (void)DTORFUNC($self);
  }
}
%enddef
#define %swiglal_public_clear_EXTERNAL_STRUCT(NAME, DTORFUNC)

// Local Variables:
// mode: c
// End:
