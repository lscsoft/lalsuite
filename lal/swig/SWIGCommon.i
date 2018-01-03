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

///
/// \defgroup SWIGCommon_i Interface SWIGCommon.i
/// \ingroup lal_swig
/// \brief Common SWIG interface code.
/// \author Karl Wette
///

///
/// # General SWIG directives
///

///
/// Ensure that all LAL library bindings share type information.
///
%begin %{
#define SWIG_TYPE_TABLE swiglal
%}

///
/// Include SWIG headers.
///
%include <exception.i>
%include <typemaps.i>

///
/// Suppress some SWIG warnings.
///
#pragma SWIG nowarn=SWIGWARN_PARSE_KEYWORD
#pragma SWIG nowarn=SWIGWARN_LANG_VARARGS_KEYWORD
#pragma SWIG nowarn=SWIGWARN_LANG_OVERLOAD_KEYWORD

///
/// Turn on auto-documentation of functions.
///
%feature("autodoc", 1);

///
/// Enable keyword arguments for functions.
///
%feature("kwargs", 1);

///
/// # Public macros
///

///
/// Public macros (i.e. those used in headers) are contained inside the \b SWIGLAL() macro, so that
/// they can be easily removed from the proprocessing interface. The \b SWIGLAL_CLEAR() macro is
/// used to clear any effects of a public macro. Calls to \b SWIGLAL_CLEAR() are generated from the
/// preprocessing interface.
///
#define SWIGLAL(...) %swiglal_public_##__VA_ARGS__
#define SWIGLAL_CLEAR(...) %swiglal_public_clear_##__VA_ARGS__

///
/// # Utility macros
///

///
/// The macro \b %swiglal_define2 defines a SWIG preprocessor symbol \c NAME1_NAME2.
///
%define %swiglal_define2(NAME1, NAME2)
%define NAME1##_##NAME2 %enddef
%enddef

///
/// The macro \b %swiglal_map() maps a one-argument \c MACRO(X) onto a list of arguments (which may
/// be empty). Based on SWIG's \c %formacro().
///
%define %_swiglal_map(MACRO, X, ...)
#if #X != ""
MACRO(X);
%_swiglal_map(MACRO, __VA_ARGS__);
#endif
%enddef
%define %swiglal_map(MACRO, ...)
%_swiglal_map(MACRO, __VA_ARGS__, );
%enddef

///
/// The macro \b %swiglal_map_a() maps a two-argument \c MACRO(A,X) onto a list of arguments (which
/// may be empty), with a common first argument \c A.
///
%define %_swiglal_map_a(MACRO, A, X, ...)
#if #X != ""
MACRO(A, X);
%_swiglal_map_a(MACRO, A, __VA_ARGS__);
#endif
%enddef
%define %swiglal_map_a(MACRO, A, ...)
%_swiglal_map_a(MACRO, A, __VA_ARGS__, );
%enddef

///
/// The macro \b %swiglal_map_ab() maps a three-argument \c MACRO(A,B,X) onto a list of arguments
/// (which may be empty), with common first arguments \c A and \c B.
///
%define %_swiglal_map_ab(MACRO, A, B, X, ...)
#if #X != ""
MACRO(A, B, X);
%_swiglal_map_ab(MACRO, A, B, __VA_ARGS__);
#endif
%enddef
%define %swiglal_map_ab(MACRO, A, B, ...)
%_swiglal_map_ab(MACRO, A, B, __VA_ARGS__, );
%enddef

///
/// The macro \b %swiglal_map_abc() maps a four-argument \c MACRO(A,B,C,X) onto a list of arguments
/// (which may be empty), with common first arguments \c A, \c B, and \c C.
///
%define %_swiglal_map_abc(MACRO, A, B, C, X, ...)
#if #X != ""
MACRO(A, B, C, X);
%_swiglal_map_abc(MACRO, A, B, C, __VA_ARGS__);
#endif
%enddef
%define %swiglal_map_abc(MACRO, A, B, C, ...)
%_swiglal_map_abc(MACRO, A, B, C, __VA_ARGS__, );
%enddef

///
/// The macros \b %swiglal_apply() and \b %swiglal_clear() apply and clear SWIG typemaps.
///
%define %swiglal_apply(TYPEMAP, TYPE, NAME)
%apply TYPEMAP { TYPE NAME };
%enddef
%define %swiglal_clear(TYPE, NAME)
%clear TYPE NAME;
%enddef

///
/// The macros \b %swiglal_feature and \b %swiglal_feature_nspace apply SWIG features.
///
%define %swiglal_feature(FEATURE, VALUE, NAME)
%feature(FEATURE, VALUE) NAME;
%enddef
%define %swiglal_feature_nspace(FEATURE, VALUE, NSPACE, NAME)
%feature(FEATURE, VALUE) NSPACE::NAME;
%enddef

///
/// The macros \b %swiglal_warnfilter and \b %swiglal_warnfilter_nspace suppress SWIG warnings.
///
%define %swiglal_warnfilter(WARNING, NAME)
%warnfilter(WARNING) NAME;
%enddef
%define %swiglal_warnfilter_nspace(WARNING, NSPACE, NAME)
%warnfilter(WARNING) NSPACE::NAME;
%enddef

///
/// The macro \b %swiglal_ignore_nspace ignores a symbol.
///
%define %swiglal_ignore_nspace(NSPACE, NAME)
%ignore NSPACE::NAME;
%enddef

///
/// These macros allocate and/or copy new instances and arrays.  They are analogous to
/// similarly-named SWIG macros but use XLALCalloc().
///
#define %swiglal_new_instance(TYPE...) %reinterpret_cast(XLALCalloc(1, sizeof(TYPE)), TYPE*)
#define %swiglal_new_copy(VAL, TYPE...) %reinterpret_cast(memcpy(%swiglal_new_instance(TYPE), &(VAL), sizeof(TYPE)), TYPE*)
#define %swiglal_new_array(SIZE, TYPE...) %reinterpret_cast(XLALCalloc(SIZE, sizeof(TYPE)), TYPE*)
#define %swiglal_new_copy_array(PTR, SIZE, TYPE...) %reinterpret_cast(memcpy(%swiglal_new_array(SIZE, TYPE), PTR, sizeof(TYPE)*(SIZE)), TYPE*)

///
/// # General interface code
///

///
/// Include C99/C++ headers.
///
%header %{
#ifdef __cplusplus
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <complex>
#include <csignal>
#else
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <signal.h>
#endif
%}

///
/// Include SWIG configuration header generated from \c config.h, but don't generate wrappers for
/// any of its symbols
///
%rename("$ignore", regexmatch$name="^SWIGLAL_") "";
%include <swiglal_config.h>
%header %{
#include <swiglal_config.h>
%}

///
/// Include GSL headers.
///
#ifdef SWIGLAL_HAVE_LIBGSL
%header %{
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
%}
#endif

///
/// Add missing GSL constructor for single-precision complex numbers.
///
#ifdef SWIGLAL_HAVE_LIBGSL
%header %{
SWIGINTERNINLINE gsl_complex_float gsl_complex_float_rect(float x, float y) {
  gsl_complex_float z;
  GSL_SET_COMPLEX(&z, x, y);
  return z;
}
%}
#endif

///
/// Include LAL standard definitions header \c LALStddef.h: this contains macro definitions which
/// are required for parsing LAL headers.
///
%include <lal/LALStddef.h>

///
/// Include LAL headers.
///
%header %{
#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/XLALError.h>
#include <lal/Date.h>
%}

///
/// Version of SWIG used to generate wrapping code.
///
%constant int swig_version = SWIG_VERSION;

///
/// Convert XLAL/LAL errors into native scripting-language exceptions:
/// <dl>
///    <dt>XLAL</dt><dd>Before performing any action, clear the XLAL error number.  Check it after
///    the action is performed; if it is non-zero, raise a native scripting-language exception with
///    the appropriate XLAL error message.</dd>
///    <dt>LAL</dt><dd>For any function which has a LALStatus* argument, create a blank LALStatus
///    struct and pass its pointer to the LAL function.  After the function is called, check the
///    LALStatus statusCode; if it is non-zero, raise a generic XLAL error (to set the XLAL error
///    number), then a native scripting-language exception. The \c swiglal_check_LALStatus (C)
///    preprocessor symbol determines if the LAL error handling code is invoked; by default this
///    symbol should be undefined, i.e. functions are XLAL functions by default.</dd>
/// </dl>
///
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

///
/// # General fragments
///

///
/// Empty fragment, for fragment-generating macros.
///
%fragment("swiglal_empty_frag", "header") {}

///
/// Wrappers around SWIG's pointer to/from SWIG-wrapped scripting language object functions.
/// \b swiglal_from_<i>SWIGTYPE</i>() simply returns a SWIG-wrapped object containing the input
/// pointer.  \b swiglal_from_<i>SWIGTYPE</i>() extracts a pointer from a SWIG-wrapped object, then
/// struct-copies the pointer to the supplied output pointer.
///
%fragment("swiglal_from_SWIGTYPE", "header") {
  SWIGINTERNINLINE SWIG_Object swiglal_from_SWIGTYPE(SWIG_Object self, void *ptr, bool isptr, swig_type_info *tinfo, int tflags) {
    return SWIG_NewPointerObj(isptr ? *((void**)ptr) : ptr, tinfo, tflags);
  }
}
%fragment("swiglal_as_SWIGTYPE", "header") {
  SWIGINTERN int swiglal_as_SWIGTYPE(SWIG_Object self, SWIG_Object obj, void *ptr, size_t len, bool isptr, swig_type_info *tinfo, int tflags) {
    void *vptr = NULL;
    int res = SWIG_ConvertPtr(obj, &vptr, tinfo, tflags);
    if (!SWIG_IsOK(res)) {
      return res;
    }
    memcpy(ptr, isptr ? &vptr : vptr, len);
    return res;
  }
}

///
/// # Create constructors and destructors for structs
///

///
/// Do not generate any default (copy) contructors or destructors.
///
%nodefaultctor;
%nocopyctor;
%nodefaultdtor;

///
/// Call the destructor \c DTORFUNC, but first call \b swiglal_release_parent() to check whether
/// \c PTR own its own memory (and release any parents).
///
%define %swiglal_struct_call_dtor(DTORFUNC, PTR)
if (swiglal_release_parent(PTR)) {
  XLALClearErrno();
  (void)DTORFUNC(PTR);
  XLALClearErrno();
}
%enddef

///
/// Extend a <tt>struct TAGNAME</tt>.
%define %swiglal_struct_extend(TAGNAME, OPAQUE, DTORFUNC)
/// <ul><li>

/// If this is an opaque struct, create an empty struct to represent the opaque struct, so that SWIG
/// has something to attach the destructor to.  No constructors are generated, since it is assumed
/// that the struct will have a creator function. Otherwise, if this is not an opaque struct,
/// generate basic constructor and (shallow) copy constructor, using XLALCalloc() to allocate memory.
#if OPAQUE
struct TAGNAME {
};
#else
%extend TAGNAME {
  TAGNAME() {
    return %swiglal_new_instance(struct TAGNAME);
  }
  TAGNAME(const struct TAGNAME *src) {
    return %swiglal_new_copy(*src, struct TAGNAME);
  }
}
#endif

/// </li><li>

/// Create destructor, using either the destructor function \c DTORFUNC, or else XLALFree() if
/// \c DTORFUNC is undefined, to destroy memory.
#if #DTORFUNC == ""
%extend TAGNAME {
  ~TAGNAME() {
    %swiglal_struct_call_dtor(XLALFree, $self);
  }
}
#else
%extend TAGNAME {
  ~TAGNAME() {
    %swiglal_struct_call_dtor(%arg(DTORFUNC), $self);
  }
}
#endif

/// </li><li>

/// Create shallow copy function __copy__() for the use of Python's copy.copy() function. It is
/// always defined but will fail for opaque structs, which cannot be copied.
#if !OPAQUE
%extend TAGNAME {
  struct TAGNAME *__copy__() {
    return %swiglal_new_copy(*$self, struct TAGNAME);
  }
}
#else
%extend TAGNAME {
  struct TAGNAME *__copy__() {
    XLALSetErrno(XLAL_ENOSYS); /* Silently signal an error to wrapper function */
    return NULL;
  }
}
#endif

/// </li><li>

/// Create deep copy function __deepcopy__() for the use of Python's copy.deepcopy() function. It is
/// always defined but will fail for opaque structs, which cannot be copied, and for structs with a
/// destructor, which presumably cannot be trivially copied with memcpy().
#if !OPAQUE && #DTORFUNC == ""
%extend TAGNAME {
  %typemap(in, noblock=1) const void *memo "";
  struct TAGNAME *__deepcopy__(const void *memo) {
    return %swiglal_new_copy(*$self, struct TAGNAME);
  }
  %clear const void *memo;
}
#else
%extend TAGNAME {
  %typemap(in, noblock=1) const void *memo "";
  struct TAGNAME *__deepcopy__(const void *memo) {
    XLALSetErrno(XLAL_ENOSYS); /* Silently signal an error to wrapper function */
    return NULL;
  }
  %clear const void *memo;
}
#endif

/// </li></ul>
%enddef

///
/// # Fragments and typemaps for arrays
///

///
/// Map fixed-array types to special variables of their elements, e.g.
/// <tt>$typemap(swiglal_fixarr_ltype, const int[][])</tt> returns \c int.
///
%typemap(swiglal_fixarr_ltype) SWIGTYPE "$ltype";
%typemap(swiglal_fixarr_ltype) SWIGTYPE[ANY] "$typemap(swiglal_fixarr_ltype, $*type)";
%typemap(swiglal_fixarr_ltype) SWIGTYPE[ANY][ANY] "$typemap(swiglal_fixarr_ltype, $*type)";

%typemap(swiglal_fixarr_tinfo) SWIGTYPE "$&descriptor";
%typemap(swiglal_fixarr_tinfo) SWIGTYPE * "$descriptor";
%typemap(swiglal_fixarr_tinfo) SWIGTYPE *const "$descriptor";
%typemap(swiglal_fixarr_tinfo) SWIGTYPE[ANY] "$typemap(swiglal_fixarr_tinfo, $*type)";
%typemap(swiglal_fixarr_tinfo) SWIGTYPE[ANY][ANY] "$typemap(swiglal_fixarr_tinfo, $*type)";

%typemap(swiglal_fixarr_isptr) SWIGTYPE "false";
%typemap(swiglal_fixarr_isptr) SWIGTYPE * "true";
%typemap(swiglal_fixarr_isptr) SWIGTYPE *const "true";
%typemap(swiglal_fixarr_isptr) SWIGTYPE[ANY] "$typemap(swiglal_fixarr_isptr, $*type)";
%typemap(swiglal_fixarr_isptr) SWIGTYPE[ANY][ANY] "$typemap(swiglal_fixarr_isptr, $*type)";

///
/// The conversion of C arrays to/from scripting-language arrays are performed by the following
/// functions:
/// <ul>
///    <li>\b %swiglal_array_copyin...() copies a scripting-language array into a C array.</li>
///    <li>\b %swiglal_array_copyout...() copies a C array into a scripting-language array.</li>
///    <li>\b %swiglal_array_viewin...() tries to view a scripting-language array as a C array.</li>
///    <li>\b %swiglal_array_viewout...() wraps a C array inside a scripting-language array, if this
///    is supported by the target scripting language.</li>
/// </ul>
///
/// These functions accept a subset of the following arguments:
/// <ul>
///    <li><tt>SWIG_Object obj</tt>: input scripting-language array.</li>
///    <li><tt>SWIG_Object parent</tt>: SWIG-wrapped object containing parent struct.</li>
///    <li><tt>void* ptr</tt>: pointer to C array.</li>
///    <li><tt>const size_t esize</tt>: size of one C array element, in bytes.</li>
///    <li><tt>const size_t ndims</tt>: number of C array dimensions.</li>
///    <li><tt>const size_t dims[]</tt>: length of each C array dimension.</li>
///    <li><tt>const size_t strides[]</tt>: strides of C array dimension, in number of elements.</li>
///    <li><tt>swig_type_info *tinfo</tt>: SWIG type info of the C array element datatype.</li>
///    <li><tt>const int tflags</tt>: SWIG type conversion flags.</li>
/// </ul>
///
/// Return values are:
/// <ul>
///    <li>\b %swiglal_array_copyin...(): SWIG error code as an \c int.</li>
///    <li>\b %swiglal_array_copyout...(): output scripting-language array as a \c SWIG_Object.</li>
///    <li>\b %swiglal_array_viewin...(): SWIG error code as an \c int.</li>
///    <li>\b %swiglal_array_viewout...(): output scripting-language array view as a \c SWIG_Object.</li>
/// </ul>
///

///
/// Names of array conversion functions for array type \c ACFTYPE.
///
#define %swiglal_array_copyin_func(ACFTYPE) swiglal_array_copyin_##ACFTYPE
#define %swiglal_array_copyout_func(ACFTYPE) swiglal_array_copyout_##ACFTYPE
#define %swiglal_array_viewin_func(ACFTYPE) swiglal_array_viewin_##ACFTYPE
#define %swiglal_array_viewout_func(ACFTYPE) swiglal_array_viewout_##ACFTYPE

///
/// Names of fragments containing the conversion functions for \c ACFTYPE.
///
#define %swiglal_array_copyin_frag(ACFTYPE) "swiglal_array_copyin_" %str(ACFTYPE)
#define %swiglal_array_copyout_frag(ACFTYPE) "swiglal_array_copyout_" %str(ACFTYPE)
#define %swiglal_array_viewin_frag(ACFTYPE) "swiglal_array_viewin_" %str(ACFTYPE)
#define %swiglal_array_viewout_frag(ACFTYPE) "swiglal_array_viewout_" %str(ACFTYPE)

///
/// The \b %swiglal_array_type() macro maps <tt>TYPE</tt>s of C arrays to an \c ACFTYPE of the
/// appropriate array conversion functions, using a special typemap. The typemap also ensures that
/// fragments containing the required conversion functions are included.
///
%define %_swiglal_array_type(ACFTYPE, TYPE)
%fragment("swiglal_array_frags_" %str(ACFTYPE), "header",
          fragment=%swiglal_array_copyin_frag(ACFTYPE),
          fragment=%swiglal_array_copyout_frag(ACFTYPE),
          fragment=%swiglal_array_viewin_frag(ACFTYPE),
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

///
/// Map C array <tt>TYPE</tt>s to array conversion function <tt>ACFTYPE</tt>s.
///
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

///
/// On modern systems, \c long could be either 32- or 64-bit...
///
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
%swiglal_array_int32or64_frags(%swiglal_array_viewin_frag, %swiglal_array_viewin_func, int);
%swiglal_array_int32or64_frags(%swiglal_array_viewout_frag, %swiglal_array_viewout_func, int);
%swiglal_array_int32or64_frags(%swiglal_array_copyin_frag, %swiglal_array_copyin_func, uint);
%swiglal_array_int32or64_frags(%swiglal_array_copyout_frag, %swiglal_array_copyout_func, uint);
%swiglal_array_int32or64_frags(%swiglal_array_viewin_frag, %swiglal_array_viewin_func, uint);
%swiglal_array_int32or64_frags(%swiglal_array_viewout_frag, %swiglal_array_viewout_func, uint);

///
/// Call the appropriate conversion function for C \c TYPE arrays.
///
#define %swiglal_array_copyin(TYPE) %swiglal_array_copyin_func($typemap(swiglal_array_typeid, TYPE))
#define %swiglal_array_copyout(TYPE) %swiglal_array_copyout_func($typemap(swiglal_array_typeid, TYPE))
#define %swiglal_array_viewin(TYPE) %swiglal_array_viewin_func($typemap(swiglal_array_typeid, TYPE))
#define %swiglal_array_viewout(TYPE) %swiglal_array_viewout_func($typemap(swiglal_array_typeid, TYPE))

///
/// ## Typemaps which convert to/from fixed-size arrays.
///

///
/// Type checkers for overloaded functions:
///
%typecheck(SWIG_TYPECHECK_POINTER) SWIGTYPE[ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0];
  const size_t dims[] = {$1_dim0};
  const size_t strides[] = {1};
  /* swiglal_array_typeid input type: $1_type */
  int res = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr(&temp$argnum[0]),
                                           sizeof(temp$argnum[0]), 1, dims, strides,
                                           $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                           %convertptr_flags);
  $1 = SWIG_CheckState(res);
}
%typecheck(SWIG_TYPECHECK_POINTER) SWIGTYPE[ANY][ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0][$1_dim1];
  const size_t dims[] = {$1_dim0, $1_dim1};
  const size_t strides[] = {$1_dim1, 1};
  /* swiglal_array_typeid input type: $1_type */
  int res = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr(&temp$argnum[0]),
                                           sizeof(temp$argnum[0][0]), 2, dims, strides,
                                           $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                           %convertptr_flags);
  $1 = SWIG_CheckState(res);
}

///
/// Input typemaps for functions and structs:
///
%typemap(in, noblock=1) SWIGTYPE[ANY], SWIGTYPE INOUT[ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0];
  $1 = &temp$argnum[0];
  {
    const size_t dims[] = {$1_dim0};
    const size_t strides[] = {1};
    /* swiglal_array_typeid input type: $1_type */
    int res = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr($1),
                                             sizeof($1[0]), 1, dims, strides,
                                             $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                             $disown | %convertptr_flags);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
}
%typemap(in, noblock=1) SWIGTYPE[ANY][ANY], SWIGTYPE INOUT[ANY][ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0][$1_dim1];
  $1 = &temp$argnum[0];
  {
    const size_t dims[] = {$1_dim0, $1_dim1};
    const size_t strides[] = {$1_dim1, 1};
    /* swiglal_array_typeid input type: $1_type */
    int res = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr($1),
                                             sizeof($1[0][0]), 2, dims, strides,
                                             $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                             $disown | %convertptr_flags);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
}

///
/// Input typemaps for global variables:
///
%typemap(varin) SWIGTYPE[ANY] {
  const size_t dims[] = {$1_dim0};
  const size_t strides[] = {1};
  /* swiglal_array_typeid input type: $1_type */
  int res = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr($1),
                                           sizeof($1[0]), 1, dims, strides,
                                           $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                           %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    %variable_fail(res, "$type", $symname);
  }
}
%typemap(varin) SWIGTYPE[ANY][ANY] {
  const size_t dims[] = {$1_dim0, $1_dim1};
  const size_t strides[] = {$1_dim1, 1};
  /* swiglal_array_typeid input type: $1_type */
  int res = %swiglal_array_copyin($1_type)(swiglal_no_self(), $input, %as_voidptr($1),
                                           sizeof($1[0][0]), 2, dims, strides,
                                           $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                           %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    %variable_fail(res, "$type", $symname);
  }
}

///
/// Output typemaps for functions and structs:
///
%typemap(out) SWIGTYPE[ANY] {
  const size_t dims[] = {$1_dim0};
  const size_t strides[] = {1};
  /* swiglal_array_typeid input type: $1_type */
%#if $owner & SWIG_POINTER_OWN
  %set_output(%swiglal_array_copyout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                              sizeof($1[0]), 1, dims, strides,
                                              $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                              $owner | %newpointer_flags));
%#else
  %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                              sizeof($1[0]), 1, dims, strides,
                                              $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                              $owner | %newpointer_flags));
%#endif
}
%typemap(out) SWIGTYPE[ANY][ANY] {
  const size_t dims[] = {$1_dim0, $1_dim1};
  const size_t strides[] = {$1_dim1, 1};
  /* swiglal_array_typeid input type: $1_type */
%#if $owner & SWIG_POINTER_OWN
  %set_output(%swiglal_array_copyout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                              sizeof($1[0][0]), 2, dims, strides,
                                              $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                              $owner | %newpointer_flags));
%#else
  %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                              sizeof($1[0][0]), 2, dims, strides,
                                              $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                              $owner | %newpointer_flags));
%#endif
}

///
/// Output typemaps for global variables:
///
%typemap(varout) SWIGTYPE[ANY] {
  const size_t dims[] = {$1_dim0};
  const size_t strides[] = {1};
  /* swiglal_array_typeid input type: $1_type */
  %set_output(%swiglal_array_viewout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                              sizeof($1[0]), 1, dims, strides,
                                              $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                              %newpointer_flags));
}
%typemap(varout) SWIGTYPE[ANY][ANY] {
  const size_t dims[] = {$1_dim0, $1_dim1};
  const size_t strides[] = {$1_dim1, 1};
  /* swiglal_array_typeid input type: $1_type */
  %set_output(%swiglal_array_viewout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                              sizeof($1[0][0]), 2, dims, strides,
                                              $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                              %newpointer_flags));
}

///
/// Argument-output typemaps for functions:
///
%typemap(in, noblock=1, numinputs=0) SWIGTYPE OUTPUT[ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0];
  $1 = &temp$argnum[0];
}
%typemap(argout) SWIGTYPE OUTPUT[ANY], SWIGTYPE INOUT[ANY] {
  const size_t dims[] = {$1_dim0};
  const size_t strides[] = {1};
  /* swiglal_array_typeid input type: $1_type */
  %append_output(%swiglal_array_copyout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                                 sizeof($1[0]), 1, dims, strides,
                                                 $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                                 SWIG_POINTER_OWN | %newpointer_flags));
}
%typemap(in, noblock=1, numinputs=0) SWIGTYPE OUTPUT[ANY][ANY] {
  $typemap(swiglal_fixarr_ltype, $1_type) temp$argnum[$1_dim0][$1_dim1];
  $1 = &temp$argnum[0];
}
%typemap(argout) SWIGTYPE OUTPUT[ANY][ANY], SWIGTYPE INOUT[ANY][ANY] {
  const size_t dims[] = {$1_dim0, $1_dim1};
  const size_t strides[] = {$1_dim1, 1};
  /* swiglal_array_typeid input type: $1_type */
  %append_output(%swiglal_array_copyout($1_type)(swiglal_no_self(), %as_voidptr($1),
                                                 sizeof($1[0][0]), 2, dims, strides,
                                                 $typemap(swiglal_fixarr_isptr, $1_type), $typemap(swiglal_fixarr_tinfo, $1_type),
                                                 SWIG_POINTER_OWN | %newpointer_flags));
}

///
/// Public macros to make fixed <i>n</i>D arrays:
/// <ul>
///    <li>output-only arguments: <b>SWIGLAL(OUTPUT_ARRAY_<i>n</i>D(TYPE, ...))</b></li>
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
///    <li>input-output arguments: <b>SWIGLAL(INOUT_ARRAY_<i>n</i>D(TYPE, ...))</b></li>
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
/// </ul>

///
/// ## Typemaps which convert to/from dynamically-sized arrays.
///

///
/// Get the correct descriptor for a dynamic array element.  Always return a pointer description,
/// even for non-pointer types, and determine whether array is an array of pointers or of data
/// blocks.
///
%typemap(swiglal_dynarr_tinfo) SWIGTYPE  "$&descriptor";
%typemap(swiglal_dynarr_tinfo) SWIGTYPE* "$descriptor";

%typemap(swiglal_dynarr_isptr) SWIGTYPE  "false";
%typemap(swiglal_dynarr_isptr) SWIGTYPE* "true";

///
/// Create immutable members for accessing the array's dimensions.  \c NI is the name of the
/// dimension member, and \c SIZET is its type.
///
%define %swiglal_array_dynamic_size(SIZET, NI)
%feature("action") NI {
  result = %static_cast(arg1->NI, SIZET);
}
%extend {
  const SIZET NI;
}
%feature("action", "") NI;
%enddef

///
/// Check that array strides are non-zero, otherwise fail.
///
%define %swiglal_array_dynamic_check_strides(NAME, DATA, I)
if (strides[I-1] == 0) {
  SWIG_exception_fail(SWIG_IndexError, "Stride of dimension "#I" of "#NAME"."#DATA" is zero");
}
%enddef

///
/// The \b %swiglal_array_dynamic_<i>n</i>D() macros create typemaps which convert <i>n</i>D
/// dynamically-allocated arrays in structs \c NAME. The macros must be added inside the definition
/// of the struct, before the struct members comprising the array are defined. The \c DATA and \c
/// N{I,J} members give the array data and dimensions, \c TYPE and \c SIZET give their respective
/// types. The \c S{I,J} give the strides of the array, in number of elements.  If the sizes or
/// strides are members of the struct, <tt>arg1-\></tt> should be used to access the struct itself.
///
/// <ul><li>
/// 1-D arrays:
%define %swiglal_array_dynamic_1D(NAME, TYPE, SIZET, DATA, NI, SI)
/// <ul><li>

/// Typemaps which convert to/from the dynamically-allocated array.
%typemap(in, noblock=1) TYPE* DATA {
  if (arg1) {
    const size_t dims[] = {NI};
    const size_t strides[] = {SI};
    %swiglal_array_dynamic_check_strides(NAME, DATA, 1);
    $1 = %reinterpret_cast(arg1->DATA, TYPE*);
    /* swiglal_array_typeid input type: $1_type */
    int res = %swiglal_array_copyin($1_type)(swiglal_self(), $input, %as_voidptr($1),
                                             sizeof(TYPE), 1, dims, strides,
                                             $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                             $disown | %convertptr_flags);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
}
%typemap(out, noblock=1) TYPE* DATA {
  if (arg1) {
    const size_t dims[] = {NI};
    const size_t strides[] = {SI};
    %swiglal_array_dynamic_check_strides(NAME, DATA, 1);
    $1 = %reinterpret_cast(arg1->DATA, TYPE*);
    /* swiglal_array_typeid input type: $1_type */
    %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                                sizeof(TYPE), 1, dims, strides,
                                                $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                $owner | %newpointer_flags));
  }
}

/// </li><li>

/// Clear unneeded typemaps and features.
%typemap(memberin, noblock=1) TYPE* DATA "";
%typemap(argout, noblock=1) TYPE* DATA "";
%typemap(freearg, noblock=1) TYPE* DATA "";
%feature("action") DATA "";
%feature("except") DATA "";

/// </li><li>

/// Create the array's data member.
%extend {
  TYPE *DATA;
}

/// </li><li>

/// Restore modified typemaps and features.
%feature("action", "") DATA;
%feature("except", "") DATA;
%clear TYPE* DATA;

/// </li></ul>
%enddef

/// </li><li>
/// a 1-D array of pointers to 1-D arrays; \c NI refer to the array of pointers, \c NJ/SJ refer to
/// the arrays pointed to:
%define %swiglal_array_dynamic_1d_ptr_1d(NAME, TYPE, SIZET, DATA, NI, NJ, SJ)
/// <ul><li>

/// Typemaps which convert from the dynamically-allocated array, indexed by \c arg2
%typemap(out, noblock=1) TYPE* DATA {
  if (arg1) {
    const size_t dims[] = {NJ};
    const size_t strides[] = {SJ};
    %swiglal_array_dynamic_check_strides(NAME, DATA, 1);
    if (((uint64_t)arg2) >= ((uint64_t)NI)) {
      SWIG_exception_fail(SWIG_IndexError, "Index to "#NAME"."#DATA" is outside of range [0,"#NI"]");
    }
    $1 = %reinterpret_cast(arg1->DATA, TYPE*)[arg2];
    /* swiglal_array_typeid input type: $1_type */
    %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                                sizeof(TYPE), 1, dims, strides,
                                                $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                $owner | %newpointer_flags));
  }
}

/// </li><li>

/// Clear unneeded typemaps and features.
%typemap(memberin, noblock=1) TYPE* DATA "";
%typemap(argout, noblock=1) TYPE* DATA "";
%typemap(freearg, noblock=1) TYPE* DATA "";
%feature("action") DATA "";
%feature("except") DATA "";

/// </li><li>

/// Create the array's indexed data member.
%extend {
  TYPE *DATA(const SIZET index);
}

/// </li><li>

/// Restore modified typemaps and features.
%feature("action", "") DATA;
%feature("except", "") DATA;
%clear TYPE* DATA;

/// </li></ul>
%enddef

/// </li><li>
/// 2-D arrays:
%define %swiglal_array_dynamic_2D(NAME, TYPE, SIZET, DATA, NI, NJ, SI, SJ)
/// <ul><li>

/// Typemaps which convert to/from the dynamically-allocated array.
%typemap(in, noblock=1) TYPE* DATA {
  if (arg1) {
    const size_t dims[] = {NI, NJ};
    const size_t strides[] = {SI, SJ};
    %swiglal_array_dynamic_check_strides(NAME, DATA, 1);
    %swiglal_array_dynamic_check_strides(NAME, DATA, 2);
    $1 = %reinterpret_cast(arg1->DATA, TYPE*);
    /* swiglal_array_typeid input type: $1_type */
    int res = %swiglal_array_copyin($1_type)(swiglal_self(), $input, %as_voidptr($1),
                                             sizeof(TYPE), 2, dims, strides,
                                             $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                             $disown | %convertptr_flags);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
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
    /* swiglal_array_typeid input type: $1_type */
    %set_output(%swiglal_array_viewout($1_type)(swiglal_self(), %as_voidptr($1),
                                                sizeof(TYPE), 2, dims, strides,
                                                $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                $owner | %newpointer_flags));
  }
}

/// </li><li>

/// Clear unneeded typemaps and features.
%typemap(memberin, noblock=1) TYPE* DATA "";
%typemap(argout, noblock=1) TYPE* DATA "";
%typemap(freearg, noblock=1) TYPE* DATA "";
%feature("action") DATA "";
%feature("except") DATA "";

/// </li><li>

/// Create the array's data member.
%extend {
  TYPE *DATA;
}

/// </li><li>

/// Restore modified typemaps and features.
%feature("action", "") DATA;
%feature("except", "") DATA;
%clear TYPE* DATA;

/// </li></ul>
%enddef

/// </li></ul>
///

///
/// The \b %swiglal_array_struct_<i>n</i>D() macros create typemaps which attempt to view a
/// scripting-language array as a C array <tt>struct NAME</tt>. If the input is not already a
/// SWIG-wrapped object wrapping a <tt>struct NAME</tt>, an input view is attempted using \b
/// %swiglal_array_viewin...().
///
/// The main concern with passing scripting-language memory to C code is that it might try to
/// re-allocate or free the memory, which would certainly lead to undefined behaviour. Therefore, by
/// default a view is attempted only for pointers to <tt>const NAME*</tt>, since it is reasonable to
/// assume the called C code will not try to re-allocate or free constant memory. When it is known
/// that the called C function will not try to re-allocate or free a particular argument, the
/// <b>SWIGLAL(VIEWIN_ARRAYS(NAME, ...))</b> macro can be used to apply the typemap to pointers to
/// (non-<tt>const</tt>) \c NAME*. Alternatively, the <b>SWIGLAL(COPYINOUT_ARRAYS(NAME, ...))</b>
/// macro treats \c NAME* as an input-output argument, and makes an internal copy of it if necessary.
///
/// <ul><li>
/// 1-D arrays:
%define %swiglal_array_struct_1D(NAME, TYPE, SIZET, DATA, NI)
/// <ul><li>

/// Typemap which attempts to view pointers to <tt>const NAME*</tt>.
%typemap(in, noblock=1) const NAME* (void *argp = 0, int res = 0, NAME temp, void *temp_data = 0) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    typedef struct { SIZET NI; TYPE* DATA; } sizchk_t;
    if (!($disown) && sizeof(sizchk_t) == sizeof(NAME)) {
      size_t numel = 0;
      size_t dims[] = {0};
      temp.DATA = NULL;
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&temp.DATA),
                                         sizeof(TYPE), 1, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (!SWIG_IsOK(res)) {
        temp_data = temp.DATA = %reinterpret_cast(XLALMalloc(numel * sizeof(TYPE)), TYPE*);
        size_t strides[] = {1};
        res = %swiglal_array_copyin(TYPE*)(swiglal_no_self(), $input, %as_voidptr(temp.DATA),
                                           sizeof(TYPE), 1, dims, strides,
                                           $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                           $disown | %convertptr_flags);
        if (!SWIG_IsOK(res)) {
          %argument_fail(res, "$type", $symname, $argnum);
        } else {
          temp.NI = %reinterpret_cast(dims[0], SIZET);
          argp = &temp;
        }
      } else {
        temp.NI = %reinterpret_cast(dims[0], SIZET);
        argp = &temp;
      }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}
%typemap(freearg, match="in", noblock=1) const NAME* %{
  if (temp_data$argnum) {
    XLALFree(temp_data$argnum);
  }
%}

/// </li><li>

/// Typemap which attempts to view pointers to non-<tt>const</tt> \c NAME*.
%typemap(in, noblock=1) NAME* SWIGLAL_VIEWIN_ARRAY (void *argp = 0, int res = 0, NAME temp) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    typedef struct { SIZET NI; TYPE* DATA; } sizchk_t;
    if (!($disown) && sizeof(sizchk_t) == sizeof(NAME)) {
      size_t numel = 0;
      size_t dims[] = {0};
      temp.DATA = NULL;
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&temp.DATA),
                                         sizeof(TYPE), 1, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (!SWIG_IsOK(res)) {
        %argument_fail(res, "$type", $symname, $argnum);
      } else {
        temp.NI = %reinterpret_cast(dims[0], SIZET);
        argp = &temp;
      }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}

/// </li><li>

/// Typemap which treats pointers to non-<tt>const</tt> \c NAME* as input-output arguments.
/// The type of the output argument should always match that of the input argument, so:
/// - If the input argument is a SWIG-wrapped \c NAME*, just unwrap it and return a reference.
/// - If the input argument is a native scripting-language array, make an internal copy of it,
///   use the copy, and return a native scripting-language array copy of the internal copy.
%typemap(in, noblock=1) NAME* SWIGLAL_COPYINOUT_ARRAY (void *argp = 0, int res = 0, NAME temp, SWIG_Object input_ref, void *temp_data = 0) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    typedef struct { SIZET NI; TYPE* DATA; } sizchk_t;
    if (!($disown) && sizeof(sizchk_t) == sizeof(NAME)) {
      size_t numel = 0;
      size_t dims[] = {0};
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&temp_data),
                                         sizeof(TYPE), 1, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (numel > 0) {
        temp_data = temp.DATA = %reinterpret_cast(XLALMalloc(numel * sizeof(TYPE)), TYPE*);
        size_t strides[] = {1};
        res = %swiglal_array_copyin(TYPE*)(swiglal_no_self(), $input, %as_voidptr(temp_data),
                                           sizeof(TYPE), 1, dims, strides,
                                           $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                           $disown | %convertptr_flags);
        if (!SWIG_IsOK(res)) {
          %argument_fail(res, "$type", $symname, $argnum);
        } else {
          temp.NI = %reinterpret_cast(dims[0], SIZET);
          argp = &temp;
        }
      } else {
        %argument_fail(res, "$type", $symname, $argnum);
      }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  } else {
    input_ref = $input;
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}
%typemap(argout, match="in", noblock=1) NAME* SWIGLAL_COPYINOUT_ARRAY %{
  if (temp_data$argnum) {
    const size_t dims[] = {temp$argnum.NI};
    const size_t strides[] = {1};
    /* swiglal_array_typeid input type: TYPE* */
    %append_output(%swiglal_array_copyout(TYPE*)(swiglal_no_self(), %as_voidptr(temp_data$argnum),
                                                 sizeof(TYPE), 1, dims, strides,
                                                 $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                 SWIG_POINTER_OWN | %newpointer_flags));
  } else {
    %append_output(swiglal_get_reference(input_ref$argnum));
  }
%}
%typemap(freearg, match="in", noblock=1) NAME* SWIGLAL_COPYINOUT_ARRAY %{
  if (temp_data$argnum) {
    XLALFree(temp_data$argnum);
  }
%}

/// </li></ul>
%enddef

/// </li><li>
/// 2-D arrays:
%define %swiglal_array_struct_2D(NAME, TYPE, SIZET, DATA, NI, NJ)
/// <ul><li>

/// Typemap which attempts to view pointers to <tt>const NAME*</tt>.
%typemap(in, noblock=1) const NAME* (void *argp = 0, int res = 0, NAME temp, void *temp_data = 0) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    typedef struct { SIZET NI; SIZET NJ; TYPE* DATA; } sizchk_t;
    if (!($disown) && sizeof(sizchk_t) == sizeof(NAME)) {
      size_t numel = 0;
      size_t dims[] = {0, 0};
      temp.DATA = NULL;
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&temp.DATA),
                                         sizeof(TYPE), 2, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (!SWIG_IsOK(res)) {
        temp_data = temp.DATA = %reinterpret_cast(XLALMalloc(numel * sizeof(TYPE)), TYPE*);
        size_t strides[] = {dims[1], 1};
        res = %swiglal_array_copyin(TYPE*)(swiglal_no_self(), $input, %as_voidptr(temp.DATA),
                                           sizeof(TYPE), 2, dims, strides,
                                           $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                           $disown | %convertptr_flags);
        if (!SWIG_IsOK(res)) {
          %argument_fail(res, "$type", $symname, $argnum);
        } else {
          temp.NI = %reinterpret_cast(dims[0], SIZET);
          temp.NJ = %reinterpret_cast(dims[1], SIZET);
          argp = &temp;
        }
      } else {
        temp.NI = %reinterpret_cast(dims[0], SIZET);
        temp.NJ = %reinterpret_cast(dims[1], SIZET);
        argp = &temp;
      }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}
%typemap(freearg, match="in", noblock=1) const NAME* %{
  if (temp_data$argnum) {
    XLALFree(temp_data$argnum);
  }
%}

/// </li><li>

/// Typemap which attempts to view pointers to non-<tt>const</tt> \c NAME*.
%typemap(in, noblock=1) NAME* SWIGLAL_VIEWIN_ARRAY (void *argp = 0, int res = 0, NAME temp) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    typedef struct { SIZET NI; SIZET NJ; TYPE* DATA; } sizchk_t;
    if (!($disown) && sizeof(sizchk_t) == sizeof(NAME)) {
      size_t numel = 0;
      size_t dims[] = {0, 0};
      temp.DATA = NULL;
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&temp.DATA),
                                         sizeof(TYPE), 2, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (!SWIG_IsOK(res)) {
        %argument_fail(res, "$type", $symname, $argnum);
      } else {
        temp.NI = %reinterpret_cast(dims[0], SIZET);
        temp.NJ = %reinterpret_cast(dims[1], SIZET);
        argp = &temp;
      }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}

/// </li><li>

/// Typemap which treats pointers to non-<tt>const</tt> \c NAME* as input-output arguments.
/// The type of the output argument should always match that of the input argument, so:
/// - If the input argument is a SWIG-wrapped \c NAME*, just unwrap it and return a reference.
/// - If the input argument is a native scripting-language array, make an internal copy of it,
///   use the copy, and return a native scripting-language array copy of the internal copy.
%typemap(in, noblock=1) NAME* SWIGLAL_COPYINOUT_ARRAY (void *argp = 0, int res = 0, NAME temp, SWIG_Object input_ref, void *temp_data = 0) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    typedef struct { SIZET NI; SIZET NJ; TYPE* DATA; } sizchk_t;
    if (!($disown) && sizeof(sizchk_t) == sizeof(NAME)) {
      size_t numel = 0;
      size_t dims[] = {0, 0};
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&temp_data),
                                         sizeof(TYPE), 2, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (numel > 0) {
        temp_data = temp.DATA = %reinterpret_cast(XLALMalloc(numel * sizeof(TYPE)), TYPE*);
        size_t strides[] = {dims[1], 1};
        res = %swiglal_array_copyin(TYPE*)(swiglal_no_self(), $input, %as_voidptr(temp_data),
                                           sizeof(TYPE), 2, dims, strides,
                                           $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                           $disown | %convertptr_flags);
        if (!SWIG_IsOK(res)) {
          %argument_fail(res, "$type", $symname, $argnum);
        } else {
          temp.NI = %reinterpret_cast(dims[0], SIZET);
          temp.NJ = %reinterpret_cast(dims[1], SIZET);
          argp = &temp;
        }
      } else {
        %argument_fail(res, "$type", $symname, $argnum);
      }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  } else {
    input_ref = $input;
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}
%typemap(argout, match="in", noblock=1) NAME* SWIGLAL_COPYINOUT_ARRAY %{
  if (temp_data$argnum) {
    const size_t dims[] = {temp$argnum.NI, temp$argnum.NJ};
    size_t strides[] = {dims[1], 1};
    /* swiglal_array_typeid input type: TYPE* */
    %append_output(%swiglal_array_copyout(TYPE*)(swiglal_no_self(), %as_voidptr(temp_data$argnum),
                                                 sizeof(TYPE), 2, dims, strides,
                                                 $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                 SWIG_POINTER_OWN | %newpointer_flags));
  } else {
    %append_output(swiglal_get_reference(input_ref$argnum));
  }
%}
%typemap(freearg, match="in", noblock=1) NAME* SWIGLAL_COPYINOUT_ARRAY %{
  if (temp_data$argnum) {
    XLALFree(temp_data$argnum);
  }
%}

/// </li></ul>
%enddef

/// </li></ul>
///

///
/// The <b>SWIGLAL(VIEWIN_ARRAYS(NAME, ...))</b> macro can be used to apply the above input view
/// typemaps to pointers to (non-<tt>const</tt>) \c NAME*.
///
%define %swiglal_public_VIEWIN_ARRAYS(NAME, ...)
%swiglal_map_ab(%swiglal_apply, NAME* SWIGLAL_VIEWIN_ARRAY, NAME*, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_VIEWIN_ARRAYS(NAME, ...)
%swiglal_map_a(%swiglal_clear, NAME*, __VA_ARGS__);
%enddef

///
/// The <b>SWIGLAL(COPYINOUT_ARRAYS(NAME, ...))</b> treats \c NAME* as an input-output argument.
/// The type of the output argument should always match that of the input argument, i.e. either
/// a SWIG-wrapped \c NAME* struct, or a native scripting language array.
///
%define %swiglal_public_COPYINOUT_ARRAYS(NAME, ...)
%swiglal_map_ab(%swiglal_apply, NAME* SWIGLAL_COPYINOUT_ARRAY, NAME*, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_COPYINOUT_ARRAYS(NAME, ...)
%swiglal_map_a(%swiglal_clear, NAME*, __VA_ARGS__);
%enddef

///
/// The <b>SWIGLAL(ARRAY_<i>n</i>D...(NAME, ...))</b> macros should be called from within the
/// definitions of LAL structs NAME containing dynamically-allocated arrays.
/// <ul><li>
/// 1-D arrays, e.g <tt>SIZET NI; TYPE* DATA;</tt>:
%define %swiglal_public_ARRAY_1D(NAME, TYPE, DATA, SIZET, NI)
%swiglal_array_dynamic_size(SIZET, NI);
%swiglal_array_dynamic_1D(NAME, TYPE, SIZET, DATA, arg1->NI, 1);
%ignore DATA;
%ignore NI;
%swiglal_array_struct_1D(NAME, TYPE, SIZET, DATA, NI);
%enddef
#define %swiglal_public_clear_ARRAY_1D(NAME, TYPE, DATA, SIZET, NI)
/// </li><li>
/// a 1-D array of pointers to 1-D arrays, e.g. <tt>SIZET NI; SIZET NJ; ATYPE* DATA[(NI
/// members)];</tt>
%define %swiglal_public_ARRAY_1D_PTR_1D(NAME, TYPE, DATA, SIZET, NI, NJ)
%swiglal_array_dynamic_size(SIZET, NI);
%swiglal_array_dynamic_size(SIZET, NJ);
%swiglal_array_dynamic_1d_ptr_1d(NAME, TYPE, SIZET, DATA, arg1->NI, arg1->NJ, 1);
%ignore DATA;
%ignore NJ;
%enddef
#define %swiglal_public_clear_ARRAY_1D_PTR_1D(NAME, TYPE, DATA, SIZET, NI, NJ)
/// </li><li>
/// 2-D arrays of fixed-length arrays, e.g <tt>typedef ETYPE[NJ] ATYPE; SIZET NI; ATYPE* DATA;</tt>:
%define %swiglal_public_ARRAY_2D_FIXED(NAME, ETYPE, ATYPE, DATA, SIZET, NI)
%swiglal_array_dynamic_size(SIZET, NI);
%swiglal_array_dynamic_2D(NAME, ETYPE, SIZET, DATA, arg1->NI, (sizeof(ATYPE)/sizeof(ETYPE)), (sizeof(ATYPE)/sizeof(ETYPE)), 1);
%ignore DATA;
%ignore NI;
%enddef
#define %swiglal_public_clear_ARRAY_2D_FIXED(NAME, ETYPE, ATYPE, DATA, SIZET, NI)
/// </li><li>
/// 2-D arrays, e.g <tt>SIZET NI, NJ; TYPE* DATA;</tt>:
%define %swiglal_public_ARRAY_2D(NAME, TYPE, DATA, SIZET, NI, NJ)
%swiglal_array_dynamic_size(SIZET, NI);
%swiglal_array_dynamic_size(SIZET, NJ);
%swiglal_array_dynamic_2D(NAME, TYPE, SIZET, DATA, arg1->NI, arg1->NJ, arg1->NJ, 1);
%ignore DATA;
%ignore NI;
%ignore NJ;
%swiglal_array_struct_2D(NAME, TYPE, SIZET, DATA, NI, NJ);
%enddef
#define %swiglal_public_clear_ARRAY_2D(NAME, TYPE, DATA, SIZET, NI, NJ)
/// </li></ul>
///

///
/// If multiple arrays in the same struct use the same length parameter, the
/// <b>SWIGLAL(ARRAY_MULTIPLE_LENGTHS(TAGNAME, ...))</b> macro is required
/// before the struct definition to suppress warnings.
///
%define %swiglal_public_ARRAY_MULTIPLE_LENGTHS(TAGNAME, ...)
%swiglal_map_ab(%swiglal_warnfilter_nspace, SWIGWARN_PARSE_REDEFINED, TAGNAME, __VA_ARGS__);
%enddef
#define %swiglal_public_clear_ARRAY_MULTIPLE_LENGTHS(TAGNAME, ...)

///
/// # Include scripting-language-specific interface headers
///
#if defined(SWIGOCTAVE)
%include <lal/SWIGOctave.i>
#elif defined(SWIGPYTHON)
%include <lal/SWIGPython.i>
#endif

///
/// # General typemaps
///

///
/// Typemap for functions which return \c int. If these functions also return other output arguments
/// (via \c argout typemaps), the \c int return value is ignored. This is because \c int is very
/// commonly used to return an XLAL error code, which will be converted into a native
/// scripting-language exception, and so the error code itself is not needed directly. To avoid
/// having to unpack the error code when collecting the other output arguments, therefore, it is
/// ignored in the wrappings. Functions which fit this criteria but do return a useful \c int can
/// use <b>SWIGLAL(RETURN_VALUE(int, ...))</b> to disable this behaviour.
///
/// For functions, since \c %feature("new") is set, the \c out typemap will have <tt>$owner=1</tt>,
/// and the \c newfree typemap is also applied. The \c out typemap ignores the \c int return value
/// by setting the output argument list to \c VOID_Object; the wrapping function them proceeds to
/// add other output arguments to the list, if any. After this, the \c newfree typemap is triggered,
/// which appends the \c int return if the output argument list is empty, using the
/// scripting-language-specific macro \b swiglal_append_output_if_empty(). For structs,
/// <tt>$owner=0</tt>, so the int return is set straight away, and the \c newfree typemap is never
/// applied.
///
%typemap(out, noblock=1, fragment=SWIG_From_frag(int)) int SWIGLAL_MAYBE_RETURN_INT {
%#if $owner
  %set_output(VOID_Object);
%#else
  %set_output(SWIG_From(int)($1));
%#endif
}
%typemap(newfree, noblock=1, fragment=SWIG_From_frag(int)) int SWIGLAL_MAYBE_RETURN_INT {
  swiglal_append_output_if_empty(SWIG_From(int)($1));
}

///
/// Typemaps for empty arguments. These typemaps are useful when no input from the scripting
/// language is required, and an empty struct needs to be supplied to the C function. The
/// <b>SWIGLAL(EMPTY_ARGUMENT(TYPE, ...))</b> macro applies the typemap which supplies a static
/// struct, while the <b>SWIGLAL(NEW_EMPTY_ARGUMENT(TYPE, ...))</b> macro applies the typemap which
/// supplies a dynamically-allocated struct. These typemaps may cause there to be no SWIG-wrapped
/// object for the first argument; if so, \c swiglal_no_1starg is defined for the duration of the
/// wrapping function.
///
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

///
/// SWIG conversion functions for C99 integer types.  These are mapped to the corresponding basic C
/// types, conversion functions for which are supplied by SWIG.
///
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

///
/// Fragments and typemaps for the LAL ::BOOLEAN type. The fragments re-use existing
/// scriping-language conversion functions for the C/C++ boolean type.  Appropriate
/// typemaps are then generated by \c %typemaps_asvalfromn().
///
%fragment(SWIG_From_frag(BOOLEAN), "header", fragment=SWIG_From_frag(bool)) {
  SWIGINTERNINLINE SWIG_Object SWIG_From_dec(BOOLEAN)(BOOLEAN value) {
    return SWIG_From(bool)(value ? true : false);
  }
}
%fragment(SWIG_AsVal_frag(BOOLEAN), "header", fragment=SWIG_AsVal_frag(bool)) {
  SWIGINTERN int SWIG_AsVal_dec(BOOLEAN)(SWIG_Object obj, BOOLEAN *val) {
    bool v;
    int res = SWIG_AsVal(bool)(obj, val ? &v : 0);
    if (!SWIG_IsOK(res)) {
      return SWIG_TypeError;
    }
    if (val) {
      *val = v ? 1 : 0;
    }
    return res;
  }
}
%typemaps_primitive(%checkcode(BOOL), BOOLEAN);

///
/// Fragments and typemaps for LAL strings, which should be (de)allocated using XLALMalloc() and
/// XLALFree(). The fragments re-use existing scriping-language conversion functions for ordinary
/// \c char* strings. Appropriate typemaps are then generated by \c %typemaps_string_alloc(), with
/// custom memory allocators.
///
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
    /* Get pointer to scripting-language string 'slstr' and size 'slsize'. */
    /* The 'slalloc' argument indicates whether a new string was allocated. */
    int res = SWIG_AsCharPtrAndSize(obj, &slstr, &slsize, &slalloc);
    if (!SWIG_IsOK(res)) {
      return SWIG_TypeError;
    }
    /* Return the string, if needed. */
    if (pstr) {
      /* Free the LAL string if it is already allocated. */
      if (*pstr) {
        XLALFree(*pstr);
      }
      if (alloc) {
        /* Copy the scripting-language string into a LAL-managed memory string. */
        *pstr = %swiglal_new_copy_array(slstr, slsize, char);
        *alloc = SWIG_NEWOBJ;
      }
      else {
        return SWIG_TypeError;
      }
    }
    /* Return the size (length+1) of the string, if needed. */
    if (psize) {
      *psize = slsize;
    }
    /* Free the scripting-language string, if it was allocated. */
    if (slalloc == SWIG_NEWOBJ) {
      %delete_array(slstr);
    }
    return res;
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

///
/// Typemaps for string pointers.  By default, treat arguments of type \c char** as output-only
/// arguments, which do not require a scripting-language input argument, and return their results in
/// the output argument list. Also supply an \c INOUT typemap for input-output arguments, which
/// allows a scripting-language input string to be supplied. The \c INOUT typemaps can be applied as
/// needed using the <b>SWIGLAL(INOUT_STRINGS(...))</b> macro.
///
%typemap(in, noblock=1, numinputs=0) char ** (char *str = NULL, int alloc = 0) {
  $1 = %reinterpret_cast(&str, $ltype);
  alloc = 0;
}
%typemap(in, noblock=1, fragment="SWIG_AsLALcharPtrAndSize") char ** INOUT (char *str = NULL, int alloc = 0, int res = 0) {
  res = SWIG_AsLALcharPtr($input, &str, &alloc);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
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

///
/// Do not try to free <tt>const char*</tt> return arguments.
///
%typemap(newfree, noblock=1) const char* "";

///
/// Input typemap for pointer-to-<tt>const</tt> <tt>SWIGTYPE</tt>s. This typemap is identical to the
/// standard \c SWIGTYPE pointer typemap, except <tt>$disown</tt> is commented out. This prevents
/// SWIG transferring ownership of SWIG-wrapped objects when assigning to pointer-to-<tt>const</tt>
/// members of structs. In this case, it is assumed that the struct does not want to take ownership
/// of the pointer, since it cannot free it (since it is a pointer-to-<tt>const</tt>).
///
%typemap(in, noblock=1) const SWIGTYPE * (void *argp = 0, int res = 0) {
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }
  $1 = %reinterpret_cast(argp, $ltype);
}
%typemap(freearg) const SWIGTYPE * "";

///
/// Typemaps for output <tt>SWIGTYPE</tt>s. This typemaps will match either the SWIG-wrapped return
/// argument from functions (which will have the \c SWIG_POINTER_OWN bit set in <tt>$owner</tt>) or
/// return a member of a struct through a \c get functions (in which case \c SWIG_POINTER_OWN will
/// not be set). They require the following macros:
///
/// The macro \b %swiglal_store_parent() is called to store a reference to the struct containing the
/// member being accessed, in order to prevent it from being destroyed as long as the SWIG-wrapped
/// member object is in scope. The return object is then always created with \c SWIG_POINTER_OWN, so
/// that its destructor will always be called.
///
%define %swiglal_store_parent(PTR, OWNER, SELF)
%#if !(OWNER & SWIG_POINTER_OWN)
  if (%as_voidptr(PTR) != NULL) {
    swiglal_store_parent(%as_voidptr(PTR), SELF);
  }
%#endif
%enddef
///
/// The macro \b %swiglal_set_output() sets the output of the wrapping function. If the (pointer)
/// return type of the function is the same as its first argument, then
/// <tt>swiglal_return_1starg_\#\#<i>NAME</i></tt> is defined. Unless \c swiglal_no_1starg is
/// defined (in which case the first argument is being handled by e.g. the \c EMPTY_ARGUMENT
/// typemap), the macro compares the pointer of the return value (\c result) to that of the first
/// argument (\c arg1). If they're equal, the SWIG-wrapped function will return a \e reference to
/// the SWIG object wrapping the first argument, i.e. the same object with its reference count
/// incremented. That way, if the return value is assigned to a different scripting-language
/// variable than the first argument, the underlying C struct will not be destroyed until both
/// scripting-language variables are cleared.  If the pointers are not equal, or one pointer is \c
/// NULL, the macro return a SWIG object wrapping the new C struct.
///
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
%typemap(out, noblock=1) SWIGTYPE *, SWIGTYPE &, SWIGTYPE[] {
  %swiglal_store_parent($1, $owner, swiglal_self());
  %swiglal_set_output($1_name, SWIG_NewPointerObj(%as_voidptr($1), $descriptor, ($owner | %newpointer_flags) | SWIG_POINTER_OWN));
}
%typemap(out, noblock=1) const SWIGTYPE *, const SWIGTYPE &, const SWIGTYPE[] {
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

///
/// Typemaps for pointers to primitive scalars. These are treated as output-only arguments by
/// default, by globally applying the SWIG \c OUTPUT typemaps.
///
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
///
/// The INOUT typemaps can be applied instead using the macro <b>SWIGLAL(INOUT_SCALARS(TYPE,
/// ...))</b>.
///
%define %swiglal_public_INOUT_SCALARS(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, TYPE INOUT, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_INOUT_SCALARS(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
///
/// The INPUT typemaps can be applied instead using the macro <b>SWIGLAL(INPUT_SCALARS(TYPE,
/// ...))</b>.
///
%typemap(argout, noblock=1) SWIGTYPE* SWIGLAL_INPUT_SCALAR "";
%define %swiglal_public_INPUT_SCALARS(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, TYPE INPUT, TYPE, __VA_ARGS__);
%swiglal_map_ab(%swiglal_apply, SWIGTYPE* SWIGLAL_INPUT_SCALAR, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_INPUT_SCALARS(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef

///
/// Typemaps for double pointers. By default, treat arguments of type \c TYPE** as output-only
/// arguments, which do not require a scripting-language input argument, and return their results in
/// the output argument list. Also supply an \c INOUT typemap for input-output arguments, which
/// allows a scripting-language input argument to be supplied. The INOUT typemaps can be applied as
/// needed using the <b>SWIGLAL(INOUT_STRUCTS(TYPE, ...))</b> macro.
///
/// The typemaps apply the following convention for \c NULL pointers. For an argument \c arg of type
/// \c TYPE**, where <b>SWIGLAL(INOUT_STRUCTS(TYPE, arg))</b> has been applied:
/// <ul>
///    <li>if \c arg is passed a scripting language \c empty value (e.g. <tt>[]</tt> in Octave, or
///    <tt>None</tt> in Python), this is interpreted as a \c NULL pointer of type \c TYPE**, e.g.
///    <tt>TYPE** arg = NULL;</tt>.</li>
///    <li>if \c arg is passed the integer value \c 0, this is interpreted as a valid pointer of
///    type TYPE** to a NULL pointer of type TYPE*, e.g.  <tt>TYPE* ptr = NULL; TYPE** arg =
///    \&ptr;</tt>.</li>
/// </ul>
///
%typemap(in, noblock=1, numinputs=0) SWIGTYPE ** (void *argp = NULL, int owner = 0) {
  $1 = %reinterpret_cast(&argp, $ltype);
  owner = SWIG_POINTER_OWN;
}
%typemap(in, noblock=1, fragment=SWIG_AsVal_frag(int)) SWIGTYPE ** INOUT (void  *argp = NULL, int owner = 0, int res = 0) {
  res = SWIG_ConvertPtr($input, &argp, $*descriptor, ($disown | %convertptr_flags) | SWIG_POINTER_DISOWN);
  if (!SWIG_IsOK(res)) {
    int val = 0;
    res = SWIG_AsVal(int)($input, &val);
    if (!SWIG_IsOK(res) || val != 0) {
      %argument_fail(res, "$type", $symname, $argnum);
    } else {
      argp = NULL;
      $1 = %reinterpret_cast(&argp, $ltype);
      owner = SWIG_POINTER_OWN;
    }
  } else {
    if (argp == NULL) {
      $1 = NULL;
      owner = 0;
    } else {
      $1 = %reinterpret_cast(&argp, $ltype);
      owner = 0;
    }
  }
}
%typemap(argout, noblock=1) SWIGTYPE ** {
  %append_output(SWIG_NewPointerObj($1 != NULL ? %as_voidptr(*$1) : NULL, $*descriptor, owner$argnum | %newpointer_flags));
}
%typemap(freearg) SWIGTYPE ** "";
%define %swiglal_public_INOUT_STRUCTS(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE ** INOUT, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_INOUT_STRUCTS(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef

///
/// Make the wrapping of <tt>printf()</tt>-style LAL functions a little safer, as suggested in the
/// SWIG 2.0 documentation (section 13.5). These functions should now be safely able to print any
/// string, so long as the format string is named \c format or \c fmt.
///
%typemap(in, fragment="SWIG_AsLALcharPtrAndSize") (const char *SWIGLAL_PRINTF_FORMAT, ...)
(char fmt[] = "%s", char *str = 0, int alloc = 0)
{
  $1 = fmt;
  int res = SWIG_AsLALcharPtr($input, &str, &alloc);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
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

///
/// Specialised input typemap for C file pointers. Generally it is not possible to convert
/// scripting-language file objects into \c FILE*, since the scripting language may not provide
/// access to the \c FILE*, or even be using \c FILE* internally for I/O.  The \c FILE* will
/// therefore have to be supplied from another SWIG-wrapped C function.  For convenience, however,
/// we allow the user to pass integers 0, 1, or 2 in place of a \c FILE*, as an instruction to use
/// standard input, output, or error respectively.
///
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

///
/// Specialised input typemaps for a given <tt>struct TAGNAME</tt>. If conversion from a
/// SWIG-wrapped scripting language object \c INPUT fails, call the function defined in \c FRAGMENT,
/// <b>swiglal_specialised_<i>TAGNAME</i>(INPUT, OUTPUT)</b>, to try to convert \c INPUT from some
/// other value, and if successful store it in <tt>struct TAGNAME OUTPUT</tt>.  Otherwise, raise a
/// SWIG error. Separate typemaps are needed for <tt>struct TAGNAME</tt> and pointers to <tt>struct
/// TAGNAME</tt>. Typecheck typemaps are used by overloaded functions, e.g. constructors.
///
%define %swiglal_specialised_typemaps(TAGNAME, FRAGMENT)
%typemap(in, noblock=1, fragment=FRAGMENT)
  struct TAGNAME (void *argp = 0, int res = 0),
  const struct TAGNAME (void *argp = 0, int res = 0)
{
  res = SWIG_ConvertPtr($input, &argp, $&descriptor, $disown | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    res = swiglal_specialised_##TAGNAME($input, &$1);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  } else {
    if (!argp) {
      %argument_nullref("$type", $symname, $argnum);
    } else {
      $&ltype temp = %reinterpret_cast(argp, $&ltype);
      $1 = *temp;
      if (SWIG_IsNewObj(res)) {
        %delete(temp);
      }
    }
  }
}
%typemap(freearg) struct TAGNAME, const struct TAGNAME "";
%typemap(in, noblock=1, fragment=FRAGMENT)
  struct TAGNAME* (struct TAGNAME tmp, void *argp = 0, int res = 0),
  const struct TAGNAME* (struct TAGNAME tmp, void *argp = 0, int res = 0)
{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, $disown | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    res = swiglal_specialised_##TAGNAME($input, &tmp);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    } else {
      $1 = %reinterpret_cast(&tmp, $ltype);
    }
  } else {
    $1 = %reinterpret_cast(argp, $ltype);
  }
}
%typemap(freearg) struct TAGNAME*, const struct TAGNAME* "";
%typemap(typecheck, fragment=FRAGMENT, precedence=SWIG_TYPECHECK_SWIGOBJECT) struct TAGNAME, const struct TAGNAME {
  void *argp = 0;
  int res = SWIG_ConvertPtr($input, &argp, $&descriptor, 0);
  $1 = SWIG_CheckState(res);
  if (!$1) {
    struct TAGNAME tmp;
    res = swiglal_specialised_##TAGNAME($input, &tmp);
    $1 = SWIG_CheckState(res);
  }
}
%typemap(typecheck, fragment=FRAGMENT, precedence=SWIG_TYPECHECK_SWIGOBJECT) struct TAGNAME*, const struct TAGNAME* {
  void *argp = 0;
  int res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
  $1 = SWIG_CheckState(res);
  if (!$1) {
    struct TAGNAME tmp;
    res = swiglal_specialised_##TAGNAME($input, &tmp);
    $1 = SWIG_CheckState(res);
  }
}
%enddef

///
/// # General macros
///

///
/// The <b>SWIGLAL(RETURN_VOID(TYPE,...))</b> public macro can be used to ensure that the return
/// value of a function is always ignored.
///
%define %swiglal_public_RETURN_VOID(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE SWIGLAL_RETURN_VOID, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_RETURN_VOID(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
%typemap(out, noblock=1) SWIGTYPE SWIGLAL_RETURN_VOID {
  %set_output(VOID_Object);
}

///
/// The <b>SWIGLAL(RETURN_VALUE(TYPE,...))</b> public macro can be used to ensure that the return
/// value of a function is not ignored, if the return value has previously been ignored in the
/// generated wrappings.
///
%define %swiglal_public_RETURN_VALUE(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
#define %swiglal_public_clear_RETURN_VALUE(TYPE, ...)

///
/// The <b>SWIGLAL(DISABLE_EXCEPTIONS(...))</b> public macro is useful for functions which
/// manipulate XLAL error codes, which thus require XLAL exception handling to be disabled.
///
%define %swiglal_public_DISABLE_EXCEPTIONS(...)
%swiglal_map_ab(%swiglal_feature, "except", "$action", __VA_ARGS__);
%enddef
#define %swiglal_public_clear_DISABLE_EXCEPTIONS(...)

///
/// The <b>SWIGLAL(FUNCTION_POINTER(...))</b> macro can be used to create a function pointer
/// constant, for functions which need to be used as callback functions.
///
%define %swiglal_public_FUNCTION_POINTER(...)
%swiglal_map_ab(%swiglal_feature, "callback", "%sPtr", __VA_ARGS__);
%enddef
#define %swiglal_public_clear_FUNCTION_POINTER(...)

///
/// The <b>SWIGLAL(VARIABLE_ARGUMENT_LIST(...))</b> macro supports functions which require a
/// variable-length list of arguments of type \c TYPE, i.e. a list of strings. It generates SWIG
/// \e compact default arguments, i.e. only one wrapping function where all missing arguments are
/// assigned \e ENDVALUE, generates 11 additional optional arguments of type \e TYPE, and creates a
/// contract ensuring that the last argument is always \c ENDVALUE, so that the argument list is
/// terminated by \c ENDVALUE.
///
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

///
/// The <b>SWIGLAL(OWNED_BY_1ST_ARG(...))</b> macro is used when a function returns an object whose
/// memory is owned by the object supplied as the first argument to the function.  Typically this
/// occurs when the function is returning some property of its first argument. The macro applies a
/// typemap which calles \c swiglal_store_parent() to store a reference to the first argument as the
/// \c parent of the return argument, so that the parent will not be destroyed as long as the return
/// value is in scope.
///
%define %swiglal_public_OWNED_BY_1ST_ARG(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE* SWIGLAL_OWNED_BY_1ST_ARG, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_OWNED_BY_1ST_ARG(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef
%typemap(out, noblock=1) SWIGTYPE* SWIGLAL_OWNED_BY_1ST_ARG {
%#ifndef swiglal_no_1starg
  %swiglal_store_parent($1, 0, swiglal_1starg());
%#endif
  %set_output(SWIG_NewPointerObj(%as_voidptr($1), $descriptor, ($owner | %newpointer_flags) | SWIG_POINTER_OWN));
}

///
/// The <b>SWIGLAL(OWNS_THIS_ARG(...))</b> macro indicates that a function will acquire
/// ownership of a particular argument, e.g. by storing that argument in some container, and that
/// therefore the SWIG object wrapping that argument should no longer own its memory.
///
%define %swiglal_public_OWNS_THIS_ARG(TYPE, ...)
%swiglal_map_ab(%swiglal_apply, SWIGTYPE* DISOWN, TYPE, __VA_ARGS__);
%enddef
%define %swiglal_public_clear_OWNS_THIS_ARG(TYPE, ...)
%swiglal_map_a(%swiglal_clear, TYPE, __VA_ARGS__);
%enddef

///
/// The <b>SWIGLAL(EXTERNAL_STRUCT(...))</b> macro can be used to support structs which are not
/// declared in LALSuite. It treats the struct as opaque, and attaches a destructor function to it.
///
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

///
/// The <b>SWIGLAL(CAST_STRUCT_TO(...))</b> macro adds to the containing struct methods which cast
/// it to each of the given list of types. For example:
/// \code
/// typedef struct tagBase { ... } Base;
/// typedef struct tagDerived {
///   SWIGLAL(CAST_STRUCT_TO(Base));
///    ...
/// } Derived;
/// \endcode
/// adds a method <i>Base* cast2Base()</i> method to Derived. Obviously this should be a valid cast
/// for the given types! The SWIG-wrapped object returned by the <i>cast2...()</i> methods will
/// remain in scope as long as the struct that was cast from, by using a typemap similar to that of
/// the SWIGLAL(OWNED_BY_1ST_ARG(...)) macro.
///
%typemap(out, noblock=1) SWIGTYPE* SWIGLAL_RETURNS_SELF {
%#ifndef swiglal_no_1starg
  %swiglal_store_parent($1, 0, swiglal_self());
%#endif
  %set_output(SWIG_NewPointerObj(%as_voidptr($1), $descriptor, ($owner | %newpointer_flags) | SWIG_POINTER_OWN));
}
%define %swiglal_cast_struct(TOTYPE)
%extend {
  %apply SWIGTYPE* SWIGLAL_RETURNS_SELF { TOTYPE* cast2##TOTYPE() };
  TOTYPE* cast2##TOTYPE() {
    return (TOTYPE*) $self;
  }
}
%enddef
%define %swiglal_public_CAST_STRUCT_TO(...)
%swiglal_map(%swiglal_cast_struct, __VA_ARGS__);
%enddef
#define %swiglal_public_clear_CAST_STRUCT_TO(...)

///
/// The <b>SWIGLAL(IMMUTABLE_MEMBERS(TAGNAME, ...))</b> macro can be used to make the listed members
/// of the struct TAGNAME immutable.
///
%define %swiglal_public_IMMUTABLE_MEMBERS(TAGNAME, ...)
%swiglal_map_abc(%swiglal_feature_nspace, "immutable", "1", TAGNAME, __VA_ARGS__);
%enddef
#define %swiglal_public_clear_IMMUTABLE_MEMBERS(...)

///
/// The <b>SWIGLAL(IGNORE_MEMBERS(TAGNAME, ...))</b> macro can be used to ignore the listed members
/// of the struct TAGNAME.
///
%define %swiglal_public_IGNORE_MEMBERS(TAGNAME, ...)
%swiglal_map_a(%swiglal_ignore_nspace, TAGNAME, __VA_ARGS__);
%enddef
#define %swiglal_public_clear_IGNORE_MEMBERS(...)

// Local Variables:
// mode: c
// End:
