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
/// \defgroup SWIGLALAlpha_i Interface SWIGLALAlpha.i
/// \ingroup lal_swig
/// \brief SWIG code which must appear \e before the LAL headers.
/// \author Karl Wette
///

///
/// # Error handling
///

///
/// Custom GSL/LAL error handlers which raise XLAL errors, so that they will be caught by the SWIG
/// \c %exception handler (instead of aborting, which will crash the user's scripting language
/// session).
%header %{
/// <ul><li>

/// Print the supplied error message, then raise an XLAL error.
static void swig_lal_gsl_error_handler(const char *reason, const char *file, int line, int errnum) {
  XLALPrintError("GSL function failed: %s (errnum=%i)\n", reason, errnum);
  XLALError("<GSL function>", file, line, XLAL_EFAILED);
}

/// </li><li>

/// Print the supplied error message, then raise an XLAL error.
static int swig_lal_raise_hook(int sig, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  (void) vfprintf(stderr, fmt, ap);
  va_end(ap);
  (void) fprintf(stderr, "LALRaise: %s\n", strsignal(sig));
  XLALSetErrno(XLAL_EFAILED);
  return 0;
}

/// </li><li>

/// Print the supplied error message, then raise an XLAL error.
static void swig_lal_abort_hook(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  (void) vfprintf(stderr, fmt, ap);
  va_end(ap);
  XLALSetErrno(XLAL_EFAILED);
}

/// </li></ul>
%}
///

///
/// Replace default GSL/LAL error handler with nice custom handlers, and ensure default XLAL error
/// handler is used.
///
%inline %{
void swig_set_nice_error_handlers(void) {
  gsl_set_error_handler(swig_lal_gsl_error_handler);
  lalRaiseHook = swig_lal_raise_hook;
  lalAbortHook = swig_lal_abort_hook;
  XLALSetErrorHandler(XLALDefaultErrorHandler);
}
%}

///
/// Use \c abort() error handlers in GSL/LAL/XLAL, which can be useful when running scripting
/// language interpreter under a debugger.
///
%inline %{
void swig_set_nasty_error_handlers(void) {
  fprintf(stderr, "*** WARNING: GSL/LAL/XLAL functions will now abort() on error ***\n");
  gsl_set_error_handler(NULL);
  lalRaiseHook = LALRaise;
  lalAbortHook = LALAbort;
  XLALSetErrorHandler(XLALAbortErrorHandler);
}
%}

///
/// Use nice custom handlers by default.
///
%init %{
swig_set_nice_error_handlers();
%}

///
/// # GSL vectors and matrices
///

///
/// This macro create wrapping structs for GSL vectors and matrices.
%define %swig_lal_gsl_vector_matrix(BASETYPE, TYPE, NAME)
/// <ul><li>

/// GSL vector of type \c NAME.
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

/// </li><li>

/// Typemap which attempts to view pointers to \c const GSL vector.
%typemap(in, noblock=1) const gsl_vector##NAME* (void *argp = 0, int res = 0, gsl_vector##NAME##_view temp, void *temp_data = 0) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    if (!($disown)) {
      size_t numel = 0;
      size_t dims[] = {0};
      void *data = NULL;
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&data),
                                         sizeof(TYPE), 1, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (!SWIG_IsOK(res)) {
        temp_data = data = XLALMalloc(numel * sizeof(TYPE));
        size_t strides[] = {1};
        res = %swiglal_array_copyin(TYPE*)(swiglal_no_self(), $input, %as_voidptr(data),
                                           sizeof(TYPE), 1, dims, strides,
                                           $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                           $disown | %convertptr_flags);
        if (!SWIG_IsOK(res)) {
          %argument_fail(res, "$type", $symname, $argnum);
        } else {
          temp = gsl_vector##NAME##_view_array(%reinterpret_cast(data, BASETYPE*), dims[0]);
          argp = &temp.vector;
        }
      } else {
        temp = gsl_vector##NAME##_view_array(%reinterpret_cast(data, BASETYPE*), dims[0]);
        argp = &temp.vector;
    }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}
%typemap(freearg, match="in", noblock=1) const gsl_vector##NAME* %{
  if (temp_data$argnum) {
    XLALFree(temp_data$argnum);
  }
%}

/// </li><li>

/// Typemap which attempts to view pointers to non-\c const GSL vector.
%typemap(in, noblock=1) gsl_vector##NAME* SWIGLAL_VIEWIN_ARRAY (void *argp = 0, int res = 0, gsl_vector##NAME##_view temp) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    if (!($disown)) {
      size_t numel = 0;
      size_t dims[] = {0};
      void *data = NULL;
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&data),
                                         sizeof(TYPE), 1, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (!SWIG_IsOK(res)) {
        %argument_fail(res, "$type", $symname, $argnum);
      } else {
        temp = gsl_vector##NAME##_view_array(%reinterpret_cast(data, BASETYPE*), dims[0]);
        argp = &temp.vector;
      }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}

/// </li><li>

/// Typemap which treats pointers to non-\c const GSL vector as input-output arguments.
/// The type of the output argument should always match that of the input argument, so:
/// - If the input argument is a SWIG-wrapped \c NAME*, just unwrap it and return a reference.
/// - If the input argument is a native scripting-language array, make an internal copy of it,
///   use the copy, and return a native scripting-language array copy of the internal copy.
%typemap(in, noblock=1) gsl_vector##NAME* SWIGLAL_COPYINOUT_ARRAY (void *argp = 0, int res = 0, gsl_vector##NAME##_view temp, SWIG_Object input_ref, void *temp_data = 0) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    if (!($disown)) {
      size_t numel = 0;
      size_t dims[] = {0};
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&temp_data),
                                         sizeof(TYPE), 1, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (numel > 0) {
        temp_data = %reinterpret_cast(XLALMalloc(numel * sizeof(TYPE)), TYPE*);
        size_t strides[] = {1};
        res = %swiglal_array_copyin(TYPE*)(swiglal_no_self(), $input, %as_voidptr(temp_data),
                                           sizeof(TYPE), 1, dims, strides,
                                           $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                           $disown | %convertptr_flags);
        if (!SWIG_IsOK(res)) {
          %argument_fail(res, "$type", $symname, $argnum);
        } else {
          temp = gsl_vector##NAME##_view_array(%reinterpret_cast(temp_data, BASETYPE*), dims[0]);
          argp = &temp.vector;
        }
      } else {
        %argument_fail(res, "$type", $symname, $argnum);
      }
    }
  } else {
    input_ref = $input;
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}
%typemap(argout, match="in", noblock=1) gsl_vector##NAME* SWIGLAL_COPYINOUT_ARRAY %{
  if (temp_data$argnum) {
    const size_t dims[] = {temp$argnum.vector.size};
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
%typemap(freearg, match="in", noblock=1) gsl_vector##NAME* SWIGLAL_COPYINOUT_ARRAY %{
  if (temp_data$argnum) {
    XLALFree(temp_data$argnum);
  }
%}

/// </li><li>

/// GSL matrix of type \c NAME.
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

/// </li><li>

/// Typemap which attempts to view pointers to \c const GSL matrix.
%typemap(in, noblock=1) const gsl_matrix##NAME* (void *argp = 0, int res = 0, gsl_matrix##NAME##_view temp, void *temp_data = 0) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    if (!($disown)) {
      size_t numel = 0;
      size_t dims[] = {0, 0};
      void *data = NULL;
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&data),
                                         sizeof(TYPE), 2, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (!SWIG_IsOK(res)) {
        temp_data = data = XLALMalloc(numel * sizeof(TYPE));
        size_t strides[] = {dims[1], 1};
        res = %swiglal_array_copyin(TYPE*)(swiglal_no_self(), $input, %as_voidptr(data),
                                           sizeof(TYPE), 2, dims, strides,
                                           $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                           $disown | %convertptr_flags);
        if (!SWIG_IsOK(res)) {
          %argument_fail(res, "$type", $symname, $argnum);
        } else {
          temp = gsl_matrix##NAME##_view_array(%reinterpret_cast(data, BASETYPE*), dims[0], dims[1]);
          argp = &temp.matrix;
        }
      } else {
        temp = gsl_matrix##NAME##_view_array(%reinterpret_cast(data, BASETYPE*), dims[0], dims[1]);
        argp = &temp.matrix;
      }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}
%typemap(freearg, match="in", noblock=1) const gsl_matrix##NAME* %{
  if (temp_data$argnum) {
    XLALFree(temp_data$argnum);
  }
%}

/// </li><li>

/// Typemap which attempts to view pointers to non-\c const GSL matrix.
%typemap(in, noblock=1) gsl_matrix##NAME* SWIGLAL_VIEWIN_ARRAY (void *argp = 0, int res = 0, gsl_matrix##NAME##_view temp) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    if (!($disown)) {
      size_t numel = 0;
      size_t dims[] = {0, 0};
      void *data = NULL;
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&data),
                                         sizeof(TYPE), 2, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (!SWIG_IsOK(res)) {
        %argument_fail(res, "$type", $symname, $argnum);
      } else {
        temp = gsl_matrix##NAME##_view_array(%reinterpret_cast(data, BASETYPE*), dims[0], dims[1]);
        argp = &temp.matrix;
      }
    } else {
      %argument_fail(res, "$type", $symname, $argnum);
    }
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}

/// </li><li>

/// Typemap which treats pointers to non-\c const GSL matrix as input-output arguments.
/// The type of the output argument should always match that of the input argument, so:
/// - If the input argument is a SWIG-wrapped \c NAME*, just unwrap it and return a reference.
/// - If the input argument is a native scripting-language array, make an internal copy of it,
///   use the copy, and return a native scripting-language array copy of the internal copy.
%typemap(in, noblock=1) gsl_matrix##NAME* SWIGLAL_COPYINOUT_ARRAY (void *argp = 0, int res = 0, gsl_matrix##NAME##_view temp, SWIG_Object input_ref, void *temp_data = 0) %{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, 0 /*$disown*/ | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    if (!($disown)) {
      size_t numel = 0;
      size_t dims[] = {0, 0};
      /* swiglal_array_typeid input type: TYPE* */
      res = %swiglal_array_viewin(TYPE*)(swiglal_no_self(), $input, %as_voidptrptr(&temp_data),
                                         sizeof(TYPE), 2, &numel, dims,
                                         $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                         $disown | %convertptr_flags);
      if (numel > 0) {
        temp_data = %reinterpret_cast(XLALMalloc(numel * sizeof(TYPE)), TYPE*);
        size_t strides[] = {dims[1], 1};
        res = %swiglal_array_copyin(TYPE*)(swiglal_no_self(), $input, %as_voidptr(temp_data),
                                           sizeof(TYPE), 2, dims, strides,
                                           $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                           $disown | %convertptr_flags);
        if (!SWIG_IsOK(res)) {
          %argument_fail(res, "$type", $symname, $argnum);
        } else {
          temp = gsl_matrix##NAME##_view_array(%reinterpret_cast(temp_data, BASETYPE*), dims[0], dims[1]);
          argp = &temp.matrix;
        }
      } else {
        %argument_fail(res, "$type", $symname, $argnum);
      }
    }
  } else {
    input_ref = $input;
  }
  $1 = %reinterpret_cast(argp, $ltype);
%}
%typemap(argout, match="in", noblock=1) gsl_matrix##NAME* SWIGLAL_COPYINOUT_ARRAY %{
  if (temp_data$argnum) {
    const size_t dims[] = {temp$argnum.matrix.size1, temp$argnum.matrix.size2};
    const size_t strides[] = {dims[1], 1};
    /* swiglal_array_typeid input type: TYPE* */
    %append_output(%swiglal_array_copyout(TYPE*)(swiglal_no_self(), %as_voidptr(temp_data$argnum),
                                                 sizeof(TYPE), 2, dims, strides,
                                                 $typemap(swiglal_dynarr_isptr, TYPE), $typemap(swiglal_dynarr_tinfo, TYPE),
                                                 SWIG_POINTER_OWN | %newpointer_flags));
  } else {
    %append_output(swiglal_get_reference(input_ref$argnum));
  }
%}
%typemap(freearg, match="in", noblock=1) gsl_matrix##NAME* SWIGLAL_COPYINOUT_ARRAY %{
  if (temp_data$argnum) {
    XLALFree(temp_data$argnum);
  }
%}

/// </li></ul>
%enddef
///

///
/// GSL integer vectors and matrices.
///
%swig_lal_gsl_vector_matrix(short, short, _short);
%swig_lal_gsl_vector_matrix(unsigned short, unsigned short, _ushort);
%swig_lal_gsl_vector_matrix(int, int, _int);
%swig_lal_gsl_vector_matrix(unsigned int, unsigned int, _uint);
%swig_lal_gsl_vector_matrix(long, long, _long);
%swig_lal_gsl_vector_matrix(unsigned long, unsigned long, _ulong);

///
/// GSL real and complex vectors and matrices.
///
%swig_lal_gsl_vector_matrix(float, float, _float);
%swig_lal_gsl_vector_matrix(double, double, ); // GSL double vector/matrix has no typename suffix.
%swig_lal_gsl_vector_matrix(float, gsl_complex_float, _complex_float);
%swig_lal_gsl_vector_matrix(double, gsl_complex, _complex);

///
/// # Specialised typemaps for ::LIGOTimeGPS
///

///
/// Specialised input typemaps for ::LIGOTimeGPS structs.  Accepts a SWIG-wrapped ::LIGOTimeGPS or a
/// double as input; in Python, also accepts any object with .gpsSeconds and .gpsNanoSeconds attributes.
///
%fragment("swiglal_specialised_tagLIGOTimeGPS", "header", fragment=SWIG_AsVal_frag(double), fragment=SWIG_AsVal_frag(int32_t)) {
  int swiglal_specialised_tagLIGOTimeGPS(SWIG_Object in, LIGOTimeGPS *out) {
    double val = 0;
    int res = SWIG_AsVal(double)(in, &val);
    if (!SWIG_IsOK(res)) {
#ifdef SWIGPYTHON
      if (PyObject_HasAttrString(in, "gpsSeconds") && PyObject_HasAttrString(in, "gpsNanoSeconds")) {
        int32_t gpsSeconds = 0, gpsNanoSeconds = 0;
        res = SWIG_AsVal(int32_t)(PyObject_GetAttrString(in, "gpsSeconds"), &gpsSeconds);
        if (!SWIG_IsOK(res)) {
          return res;
        }
        res = SWIG_AsVal(int32_t)(PyObject_GetAttrString(in, "gpsNanoSeconds"), &gpsNanoSeconds);
        if (!SWIG_IsOK(res)) {
          return res;
        }
        XLALGPSSet(out, gpsSeconds, gpsNanoSeconds);
        return SWIG_OK;
      }
#endif
      return res;
    }
    XLALGPSSetREAL8(out, val);
    return SWIG_OK;
  }
}
%swiglal_specialised_typemaps(tagLIGOTimeGPS, "swiglal_specialised_tagLIGOTimeGPS");

///
/// # Specialised typemaps for ::LALUnit
///

///
/// Specialised input typemaps for ::LALUnit structs.  Accepts a SWIG-wrapped ::LALUnit or
/// power-of-10 double as input.
///
%fragment("swiglal_specialised_tagLALUnit", "header", fragment="SWIG_AsCharPtr", fragment=SWIG_AsVal_frag(double)) {
  int swiglal_specialised_tagLALUnit(SWIG_Object in, LALUnit *out) {
    char *str = 0;
    int alloc = 0;
    int res = SWIG_AsCharPtr(in, &str, &alloc);
    if (SWIG_IsOK(res)) {
      if (XLALParseUnitString(out, str) == NULL) {
        res = SWIG_ValueError;
      }
    }
    if (alloc == SWIG_NEWOBJ) {
      %delete_array(str);
    }
    if (SWIG_IsOK(res)) {
      return res;
    }
    double val = 0;
    res = SWIG_AsVal(double)(in, &val);
    if (!SWIG_IsOK(res)) {
      return res;
    }
    if (val <= 0) {
      return SWIG_ValueError;
    }
    double pow10 = 0;
    if (modf(log10(val), &pow10) != 0) {
      return SWIG_ValueError;
    }
    if (pow10 < INT16_MIN || INT16_MAX < pow10) {
      return SWIG_ValueError;
    }
    *out = lalDimensionlessUnit;
    out->powerOfTen = (INT2)pow10;
    return SWIG_OK;
  }
}
%swiglal_specialised_typemaps(tagLALUnit, "swiglal_specialised_tagLALUnit");

// Local Variables:
// mode: c
// End:
