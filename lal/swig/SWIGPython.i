//
// Copyright (C) 2011--2014, 2022 Karl Wette
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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301  USA
//

// SWIG interface code specific to Python.
// Author: Karl Wette

//
// General SWIG directives and interface code
//

// In SWIG Python modules, everything is namespaced, so it makes sense to rename symbols to remove
// superfluous C-API prefixes.
#define SWIGLAL_MODULE_RENAME_CONSTANTS
#define SWIGLAL_MODULE_RENAME_FUNCTIONS
#define SWIGLAL_MODULE_RENAME_TDSTRUCTS
#define SWIGLAL_MODULE_RENAME_VARIABLES

// Include SWIG Python headers.
%include <pycomplex.swg>

// Include NumPy headers in wrapping code, and ensure that NumPy array module is loaded along with
// this module.
%header %{
#include <numpy/arrayobject.h>
%}
%init %{
import_array();
%}

// Evaluates true if a PyObject is not empty, false otherwise.
%header %{
#define swiglal_not_empty(v)  ((v) != NULL)
%}

// Name of PyObject containing the SWIG wrapping of the struct whose members are being accessed.
%header %{
#define swiglal_self()    (self)
#define swiglal_no_self() (NULL)
%}

// Name of PyObject containing the SWIG wrapping of the first argument to a function.
%header %{
#define swiglal_1starg()  (obj0)
%}

// Return a reference to the supplied PyObject; increment its reference count, then return it.
%header %{
SWIGINTERNINLINE PyObject* swiglal_get_reference(PyObject* v) { Py_XINCREF(v); return v; }
%}

// Append an argument to the output argument list of an Python SWIG-wrapped function, if the list is
// empty.
%header %{
#define swiglal_append_output_if_empty(v) if (resultobj == Py_None) resultobj = SWIG_Python_AppendOutput(resultobj, v)
%}

// Evaluates true if a PyObject represents a null pointer, false otherwise.
%header %{
#define swiglal_null_ptr(v)  ((v) == Py_None)
%}

// Python-specific function for standard output/error redirection
%header %{
SWIGINTERN int swiglal_output_stdouterr(void) {

  // Flush and rewind temporary files
  fflush(swiglal_tmp_stdout);
  rewind(swiglal_tmp_stdout);
  fflush(swiglal_tmp_stderr);
  rewind(swiglal_tmp_stderr);

  // Write standard output
  {
    char buf[512];
    while (fgets(buf, sizeof(buf), swiglal_tmp_stdout) != NULL) {
      PySys_WriteStdout("%s", buf);
    }
  }

  // Write standard error
  {
    char buf[512];
    while (fgets(buf, sizeof(buf), swiglal_tmp_stderr) != NULL) {
      PySys_WriteStderr("%s", buf);
    }
  }

  // Close temporary files
  fclose(swiglal_tmp_stdout);
  fclose(swiglal_tmp_stderr);

  return 1;

}
%}

//
// SWIG directives for operators
//

// These macros apply the correct python:slot directives to map Python __operator__ functions (which
// may be defined in %extend) to the correct PyTypeObject slots.

// Unary operators which do not return a new object.
%define %swiglal_py_ury_op(NAME, FUNCTYPE, SLOT)
%pythonmaybecall *::__##NAME##__;
%feature("python:slot", #SLOT, functype=#FUNCTYPE) *::__##NAME##__;
%feature("kwargs", 0) *::__##NAME##__;
%enddef
%swiglal_py_ury_op(float, unaryfunc, nb_float);
%swiglal_py_ury_op(hash, hashfunc, tp_hash);
%swiglal_py_ury_op(int, unaryfunc, nb_int);
%swiglal_py_ury_op(long, unaryfunc, nb_long);
%swiglal_py_ury_op(nonzero, inquiry, nb_nonzero);
%swiglal_py_ury_op(repr, reprfunc, tp_repr);
%swiglal_py_ury_op(str, reprfunc, tp_str);

// Unary operators which return a new object, and thus require %newobject to be set.
%define %swiglal_py_urn_op(NAME, FUNCTYPE, SLOT)
%newobject *::__##NAME##__;
%pythonmaybecall *::__##NAME##__;
%feature("python:slot", #SLOT, functype=#FUNCTYPE) *::__##NAME##__;
%feature("kwargs", 0) *::__##NAME##__;
%enddef
%swiglal_py_urn_op(abs, unaryfunc, nb_absolute);
%swiglal_py_urn_op(neg, unaryfunc, nb_negative);
%swiglal_py_urn_op(pos, unaryfunc, nb_positive);

// Binary operators, which always must return a new object, and thus require %newobject to be
// set. The SWIG Python module with -builtin does not support reverse operators, so they are removed
// from the interface.
%define %swiglal_py_bin_op(NAME, FUNCTYPE, SLOT)
%newobject *::__##NAME##__;
%pythonmaybecall *::__##NAME##__;
%feature("python:slot", #SLOT, functype=#FUNCTYPE) *::__##NAME##__;
%feature("kwargs", 0) *::__##NAME##__;
%ignore *::__r##NAME##__;
%enddef
%swiglal_py_bin_op(add, binaryfunc, nb_add);
%swiglal_py_bin_op(and, binaryfunc, nb_and);
%swiglal_py_bin_op(div, binaryfunc, nb_divide);
%swiglal_py_bin_op(lshift, binaryfunc, nb_lshift);
%swiglal_py_bin_op(mod, binaryfunc, nb_remainder);
%swiglal_py_bin_op(mul, binaryfunc, nb_multiply);
%swiglal_py_bin_op(or, binaryfunc, nb_or);
%swiglal_py_bin_op(pow, ternaryfunc, nb_power);
%swiglal_py_bin_op(rshift, binaryfunc, nb_rshift);
%swiglal_py_bin_op(sub, binaryfunc, nb_subtract);
%swiglal_py_bin_op(xor, binaryfunc, nb_xor);

// Python __pow__() operator must accept 3 arguments, but we do not use the 3rd.
%typemap(in) void* SWIGLAL_OP_POW_3RDARG {
  if ($input != Py_None) {
    %argument_fail(SWIG_TypeError, "$type", $symname, $argnum);
  }
}

// SWIG's SWIGPY_HASHFUNC_CLOSURE() macro requires the return of a __hash__()
// function to be of PyLong type, which is not guaranteed by SWIG_from_long()
// (which may return a PyInt), so use this custom typemap to guarantee this.
%typemap(out, noblock=1) long __hash__ {
  %set_output(PyLong_FromLong($1));
}

// Comparison operators.
%typemap(in, numinputs=0, noblock=1) int SWIGLAL_CMP_OP_RETN_HACK "";
%define %swiglal_py_cmp_op(NAME, COMPTYPE)
%pythonmaybecall *::__##NAME##__;
%feature("python:compare", #COMPTYPE) *::__##NAME##__;
%feature("kwargs", 0) *::__##NAME##__;
%feature("new", 1) *::__##NAME##__;
%typemap(out, noblock=1, fragment=SWIG_From_frag(bool)) bool __##NAME##__ {
  return SWIG_From_bool($1);
}
%typemap(freearg, noblock=1) int SWIGLAL_CMP_OP_RETN_HACK {
  PyErr_Clear();
  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}
%enddef
%swiglal_py_cmp_op(eq, Py_EQ);
%swiglal_py_cmp_op(ge, Py_GE);
%swiglal_py_cmp_op(gt, Py_GT);
%swiglal_py_cmp_op(le, Py_LE);
%swiglal_py_cmp_op(lt, Py_LT);
%swiglal_py_cmp_op(ne, Py_NE);

//
// Python-specific extensions to structs
//

// Mark all SWIG wrappings of C structs as "non-dynamic", so that they do not
// allow arbitrary attributes to be set. This prevents bugs due to trivial
// typos/misspellings of C struct fields; without %pythonnondynamic, if a C
// struct has a field "foo" but one misspells it as "Foo" in assigning the field
// in Python, SWIG will silently store "Foo" somewhere in the Python wrapper but
// will not modify the underlying C struct field "foo".
%pythonnondynamic;

// Extend a struct TAGNAME.
%define %swiglal_struct_extend_specific(TAGNAME, OPAQUE, DTORFUNC)

// Create shallow copy function __copy__() for the use of Python's copy.copy() function. It is
// always defined but will fail for opaque structs, which cannot be copied.
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

// Create deep copy function __deepcopy__() for the use of Python's copy.deepcopy() function. It is
// always defined but will fail for opaque structs, which cannot be copied, and for structs with a
// destructor, which presumably cannot be trivially copied with memcpy().
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

%enddef // %swiglal_struct_extend_specific

//
// General fragments, typemaps, and macros
//

// SWIG conversion fragments and typemaps for GSL complex numbers.
%swig_cplxflt_convn(gsl_complex_float, gsl_complex_float_rect, GSL_REAL, GSL_IMAG);
%swig_cplxdbl_convn(gsl_complex, gsl_complex_rect, GSL_REAL, GSL_IMAG);
%typemaps_primitive(%checkcode(CPLXFLT), gsl_complex_float);
%typemaps_primitive(%checkcode(CPLXDBL), gsl_complex);

// SWIG conversion fragments and typemaps for LAL complex numbers.
%swig_cplxflt_convn(COMPLEX8, crectf, crealf, cimagf);
%swig_cplxdbl_convn(COMPLEX16, crect, creal, cimag);
%typemaps_primitive(%checkcode(CPLXFLT), COMPLEX8);
%typemaps_primitive(%checkcode(CPLXDBL), COMPLEX16);

// Handle NumPy fixed-width integer/float types
// - Since all SWIG integer type conversions ultimately use either SWIG_AsVal_long() or
//   SWIG_AsVal_unsigned_SS_long(), and all SWIG floating-point type conversions ultimately use
//   SWIG_AsVal_float() or SWIG_AsVal_double(), it is most straightforward to replace those
//   functions with custom versions using C preprocessor macros
// - SWIGLAL maps complex floating-point types to COMPLEX{8|16} via %swig_cplx{flt|dbl}_convn()
%fragment(SWIG_AsVal_frag(long));
%fragment(SWIG_AsVal_frag(unsigned long));
%fragment(SWIG_AsVal_frag(float));
%fragment(SWIG_AsVal_frag(double));
%fragment(SWIG_AsVal_frag(COMPLEX8));
%fragment(SWIG_AsVal_frag(COMPLEX16));
%header %{
SWIGINTERN int swiglal_SWIG_AsVal_long(PyObject *obj, long* val) {
  if (PyArray_IsScalar(obj, Integer)) {
    /* handle NumPy signed integer types */
    if (val) {
      PyArray_Descr *longDescr = PyArray_DescrFromType(NPY_LONG);
      PyArray_CastScalarToCtype(obj, (void*)val, longDescr);
      Py_DECREF(longDescr);
    }
    return SWIG_OK;
  }
  /* fall back to SWIG default behaviour */
  return SWIG_AsVal_long(obj, val);
}
SWIGINTERN int swiglal_SWIG_AsVal_unsigned_SS_long(PyObject *obj, unsigned long *val) {
  if (PyArray_IsScalar(obj, Integer)) {
    /* handle NumPy unsigned integer types */
    if (val) {
      PyArray_Descr *ulongDescr = PyArray_DescrFromType(NPY_ULONG);
      PyArray_CastScalarToCtype(obj, (void*)val, ulongDescr);
      Py_DECREF(ulongDescr);
    }
    return SWIG_OK;
  }
  /* fall back to SWIG default behaviour */
  return SWIG_AsVal_unsigned_SS_long(obj, val);
}
SWIGINTERN int swiglal_SWIG_AsVal_float(PyObject *obj, float* val) {
  if (PyArray_IsScalar(obj, Integer) || PyArray_IsScalar(obj, Floating)) {
    /* handle NumPy signed integer types */
    if (val) {
      PyArray_Descr *floatDescr = PyArray_DescrFromType(NPY_FLOAT);
      PyArray_CastScalarToCtype(obj, (void*)val, floatDescr);
      Py_DECREF(floatDescr);
    }
    return SWIG_OK;
  }
  /* fall back to SWIG default behaviour */
  return SWIG_AsVal_float(obj, val);
}
SWIGINTERN int swiglal_SWIG_AsVal_double(PyObject *obj, double* val) {
  if (PyArray_IsScalar(obj, Integer) || PyArray_IsScalar(obj, Floating)) {
    /* handle NumPy signed integer types */
    if (val) {
      PyArray_Descr *doubleDescr = PyArray_DescrFromType(NPY_DOUBLE);
      PyArray_CastScalarToCtype(obj, (void*)val, doubleDescr);
      Py_DECREF(doubleDescr);
    }
    return SWIG_OK;
  }
  /* fall back to SWIG default behaviour */
  return SWIG_AsVal_double(obj, val);
}
SWIGINTERN int swiglal_SWIG_AsVal_COMPLEX8(PyObject *obj, COMPLEX8* val) {
  if (PyArray_IsScalar(obj, Integer) || PyArray_IsScalar(obj, Floating) || PyArray_IsScalar(obj, ComplexFloating)) {
    /* handle NumPy signed integer types */
    if (val) {
      PyArray_Descr *floatComplexDescr = PyArray_DescrFromType(NPY_COMPLEX64);
      PyArray_CastScalarToCtype(obj, (void*)val, floatComplexDescr);
      Py_DECREF(floatComplexDescr);
    }
    return SWIG_OK;
  }
  /* fall back to SWIG default behaviour */
  return SWIG_AsVal_COMPLEX8(obj, val);
}
SWIGINTERN int swiglal_SWIG_AsVal_COMPLEX16(PyObject *obj, COMPLEX16* val) {
  if (PyArray_IsScalar(obj, Integer) || PyArray_IsScalar(obj, Floating) || PyArray_IsScalar(obj, ComplexFloating)) {
    /* handle NumPy signed integer types */
    if (val) {
      PyArray_Descr *doubleComplexDescr = PyArray_DescrFromType(NPY_COMPLEX128);
      PyArray_CastScalarToCtype(obj, (void*)val, doubleComplexDescr);
      Py_DECREF(doubleComplexDescr);
    }
    return SWIG_OK;
  }
  /* fall back to SWIG default behaviour */
  return SWIG_AsVal_COMPLEX16(obj, val);
}
#define SWIG_AsVal_long(obj, val) swiglal_SWIG_AsVal_long(obj, val)
#define SWIG_AsVal_unsigned_SS_long(obj, val) swiglal_SWIG_AsVal_unsigned_SS_long(obj, val)
#define SWIG_AsVal_float(obj, val) swiglal_SWIG_AsVal_float(obj, val)
#define SWIG_AsVal_double(obj, val) swiglal_SWIG_AsVal_double(obj, val)
#define SWIG_AsVal_COMPLEX8(obj, val) swiglal_SWIG_AsVal_COMPLEX8(obj, val)
#define SWIG_AsVal_COMPLEX16(obj, val) swiglal_SWIG_AsVal_COMPLEX16(obj, val)
%}

// Typemaps which convert to/from the C broken-down date/time struct.
%typemap(in) struct tm* (struct tm temptm) {

  // Set 'tm' struct to zero
  memset(&temptm, 0, sizeof(temptm));

  if ($input != NULL && $input != Py_None) {

    // Check that the $input PyObject is a sequence of 9 integer elements
    if (!PySequence_Check($input)) {
      %argument_fail(SWIG_ValueError, "$type (not a sequence)", $symname, $argnum);
    }
    if (PySequence_Size($input) != 9) {
      %argument_fail(SWIG_ValueError, "$type (must have 9 elements)", $symname, $argnum);
    }
    PyObject *seq = PySequence_Fast($input, "$type (not a sequence)");
    temptm.tm_year  = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 0)), int);
    temptm.tm_mon   = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 1)), int);
    temptm.tm_mday  = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 2)), int);
    temptm.tm_hour  = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 3)), int);
    temptm.tm_min   = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 4)), int);
    temptm.tm_sec   = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 5)), int);
    temptm.tm_wday  = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 6)), int);
    temptm.tm_yday  = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 7)), int);
    temptm.tm_isdst = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 8)), int);
    Py_CLEAR(seq);
    if (PyErr_Occurred()) {  // Catch any errors while converting items to integers
      SWIG_fail;
    }

    // Convert Python date ranges to 'tm' struct date ranges
    // Python:  1900 = year 1900, January = month 1, Monday = week day 0,
    // January 1st = year day 1
    // C:  1900 = year 0, January = month 0, Monday = week day 1, January
    // 1st = year day 0
    temptm.tm_year -= 1900;
    temptm.tm_mon  -= 1;
    temptm.tm_wday  = (temptm.tm_wday + 8) % 7;
    temptm.tm_yday -= 1;
  }

  $1 = &temptm;

}
%typemap(freearg) struct tm* "";
%typemap(out) struct tm* {

  // Convert 'tm' struct date ranges to Python date ranges
  // Python:  1900 = year 1900, January = month 1, Monday = week day 0,
  // January 1st = year day 1
  // C:  1900 = year 0, January = month 0, Monday = week day 1, January 1st
  // = year day 0
  $1->tm_year += 1900;
  $1->tm_mon  += 1;
  $1->tm_wday  = ($1->tm_wday + 6) % 7;
  $1->tm_yday += 1;

  // Build a 9-element tuple (Python struct_time is immutable)
  $result = Py_BuildValue("(iiiiiiiii)",
                          $1->tm_year, $1->tm_mon, $1->tm_mday,
                          $1->tm_hour, $1->tm_min, $1->tm_sec,
                          $1->tm_wday, $1->tm_yday, $1->tm_isdst);

}

//
// Interface code to track object parents
//

// Interface code which tracks the parent structs of SWIG-wrapped struct members, so that the parent
// struct is not destroyed as long as a SWIG-wrapped object containing any of its members exists.
%header %{

// Internal map from member pointers to PyObjects containing the member parent struct, as well as an
// internal reference count of how many SWIG-wrapped member objects are extant.
static PyObject *parent_map = NULL;

// Store a reference to the parent of ptr in the internal map. If there is already such a reference,
// increment the internal reference count instead.
SWIGINTERN void swiglal_store_parent(void* ptr, PyObject* parent) {
  PyObject *pyerr_type = NULL, *pyerr_value = NULL, *pyerr_traceback = NULL;
  PyErr_Fetch(&pyerr_type, &pyerr_value, &pyerr_traceback);
  int ecode;
  assert(ptr != NULL);
  assert(parent != NULL);
  PyObject* key = PyLong_FromVoidPtr(ptr);
  assert(key != NULL);
  PyObject* parent_tuple = PyDict_GetItem(parent_map, key);
  if (parent_tuple == NULL) {
    const long ref_count = 1;
    parent_tuple = Py_BuildValue("Ol", parent, ref_count);
    assert(parent_tuple != NULL);
    ecode = PyDict_SetItem(parent_map, key, parent_tuple);
    assert(ecode == 0);
    Py_CLEAR(parent_tuple);
  }
  else {
    Py_INCREF(parent_tuple);
    PyObject* stored_parent = NULL;
    long ref_count = 0;
    ecode = PyArg_ParseTuple(parent_tuple, "Ol", &stored_parent, &ref_count);
    assert(ecode);
    ++ref_count;
    Py_INCREF(stored_parent);
    Py_CLEAR(parent_tuple);
    parent_tuple = Py_BuildValue("Nl", stored_parent, ref_count);
    assert(parent_tuple != NULL);
    ecode = PyDict_SetItem(parent_map, key, parent_tuple);
    assert(ecode == 0);
    Py_CLEAR(parent_tuple);
  }
  Py_CLEAR(key);
  assert(PyErr_Occurred() == NULL);
  PyErr_Restore(pyerr_type, pyerr_value, pyerr_traceback);
}

// Check if ptr stored a reference to a parent struct. If there is no parent object, then ptr
// *really* owns its memory, and it's okay for it to destroy it (so return true). Otherwise,
// decrement the internal reference count, erase the parent map entry if it reaches zero, and
// return false to prevent any destructors being called.
SWIGINTERN bool swiglal_release_parent(void *ptr) {
  PyObject *pyerr_type = NULL, *pyerr_value = NULL, *pyerr_traceback = NULL;
  PyErr_Fetch(&pyerr_type, &pyerr_value, &pyerr_traceback);
  int ecode;
  bool retn = true;
  assert(ptr != NULL);
  PyObject* key = PyLong_FromVoidPtr(ptr);
  assert(key != NULL);
  PyObject* parent_tuple = PyDict_GetItem(parent_map, key);
  if (parent_tuple != NULL) {
    Py_INCREF(parent_tuple);
    retn = false;
    PyObject* stored_parent = NULL;
    long ref_count = 0;
    ecode = PyArg_ParseTuple(parent_tuple, "Ol", &stored_parent, &ref_count);
    assert(ecode);
    Py_INCREF(stored_parent);
    Py_CLEAR(parent_tuple);
    if (--ref_count == 0) {
      ecode = PyDict_DelItem(parent_map, key);
      assert(ecode == 0);
    }
    else {
      parent_tuple = Py_BuildValue("Ol", stored_parent, ref_count);
      ecode = PyDict_SetItem(parent_map, key, parent_tuple);
      assert(ecode == 0);
      Py_CLEAR(parent_tuple);
    }
    Py_CLEAR(stored_parent);
  }
  Py_CLEAR(key);
  assert(PyErr_Occurred() == NULL);
  PyErr_Restore(pyerr_type, pyerr_value, pyerr_traceback);
  return retn;
}

%} // %header
%init %{

// Get a pointer to the internal parent map. Look for an attribute 'parent_map' of an internal
// module 'swiglal_runtime_data'; if it does not exist, create a new map and assign the module
// attribute, otherwise store the attribute's value. In this way each wrapping module gets a pointer
// to the same map.
{
  const char *const module_name = "swiglal_runtime_data";
  const char *const parent_map_name = "parent_map";
#if PY_VERSION_HEX >= 0x03000000
  PyObject* module = PyImport_AddModule(module_name);
#else
  PyObject* module = Py_InitModule(module_name, NULL);
#endif
  assert(module != NULL);
  if (PyObject_HasAttrString(module, parent_map_name)) {
    parent_map = PyObject_GetAttrString(module, parent_map_name);
  }
  else {
    parent_map = PyDict_New();
    PyObject_SetAttrString(module, parent_map_name, parent_map);
  }
  assert(parent_map != NULL);
  Py_INCREF(parent_map);
}

%} // %init

//
// Fragments and typemaps for arrays
//

// This section implements array conversion functions for basic C array types, and custom NumPy
// array descriptors for viewing C arrays of object, e.g. structs.

// Fragment defining helper functions for the array conversion functions.
%fragment("swiglal_py_array_helpers", "header") {

  // Compute the scalar index of the C array element, and return a pointer to the element itself.
  void* swiglal_py_get_element_ptr(void* ptr,
                                   const size_t esize,
                                   const size_t ndims,
                                   const size_t strides[],
                                   npy_intp idx[])
  {
    size_t elemidx = 0;
    for (size_t j = 0; j < ndims; ++j) {
      elemidx += ((size_t)idx[j]) * strides[j];
    }
    return %reinterpret_cast(%reinterpret_cast(ptr, char*) + elemidx*esize, void*);
  }

  // Increment the NumPy array index in row-major order, to match the ordering of the C array.
  void swiglal_py_increment_idx(const size_t ndims,
                                const size_t dims[],
                                npy_intp idx[])
  {
    for (int j = ((int)ndims) - 1; j >= 0; --j) {
      if (++idx[j] < ((npy_intp)dims[j])) {
        break;
      }
      idx[j] = 0;
    }
  }

 } // fragment swiglal_py_array_helpers

// Fragment defining helper functions for the NumPy object-view array descriptors.
%fragment("swiglal_py_array_objview", "header") {

  // Struct which associates a SWIG type descriptor with two NumPy array descriptors, one for arrays
  // of data blocks (_noptr), and one for arrays of pointers (_isptr).
  typedef struct {
    swig_type_info* tinfo;
    PyArray_Descr* descr_noptr;
    PyArray_Descr* descr_isptr;
  } swiglal_py_array_type_pair;

  // Static array of SWIG type/NumPy array descriptor pairs. This array should always be long enough
  // to accommodate all possible swig_type_info*, since they are always members of the
  // SWIG-generated global array swig_types[].  This array in turn is always one longer than the
  // total number of types, so there should always be a sentinal NULL element at the end.
  static swiglal_py_array_type_pair swiglal_py_array_types[sizeof(swig_types) / sizeof(swig_types[0])];

  // This function maps a SWIG type descriptor to a NumPy array descriptor, or returns the first
  // NULL element if a mapping doesn't exist yet.
  SWIGINTERN PyArray_Descr** swiglal_py_array_descr_from_tinfo(const bool isptr, swig_type_info* tinfo) {
    size_t i = 0;
    while (swiglal_py_array_types[i].tinfo != NULL && swiglal_py_array_types[i].tinfo != tinfo)
      ++i;
    if (swiglal_py_array_types[i].tinfo == NULL)
      swiglal_py_array_types[i].tinfo = tinfo;
    return isptr ? &swiglal_py_array_types[i].descr_isptr : &swiglal_py_array_types[i].descr_noptr;
  }

  // This function maps a NumPy array descriptor to a SWIG type descriptor, or returns NULL element
  // if a mapping doesn't exist.
  SWIGINTERN void swiglal_py_array_tinfo_from_descr(bool *isptr, swig_type_info** tinfo, PyArray_Descr* descr) {
    size_t i = 0;
    while ( ( swiglal_py_array_types[i].descr_noptr != NULL  || swiglal_py_array_types[i].descr_isptr != NULL  ) &&
            ( swiglal_py_array_types[i].descr_noptr != descr && swiglal_py_array_types[i].descr_isptr != descr ) )
      ++i;
    *isptr = (swiglal_py_array_types[i].descr_isptr == descr);
    *tinfo = swiglal_py_array_types[i].tinfo;
  }

  // Array of NumPy types that a NumPy object-view array can be safely cast to.
  static int swiglal_py_array_objview_copyswap_cancastto[2] = {NPY_OBJECT, NPY_NOTYPE};

  // NumPy array descriptor function for copying/byte-swapping an array element.
  static void swiglal_py_array_objview_copyswap(void* dst, void* src, int swap, void* arr) {

    // Check input.
    assert(arr != NULL);
    PyArrayObject* nparr = (PyArrayObject*)arr;
    assert(PyArray_DESCR(nparr) != NULL);

    // Copy array element.
    if (src != NULL) {
      memcpy(dst, src, PyArray_DESCR(nparr)->elsize);
    }

    // Byte-swap array element, if required.
    if (swap) {
      const size_t n = PyArray_DESCR(nparr)->elsize / 2;
      char *a, *b, c;
      a = (char *)dst;
      b = a + (PyArray_DESCR(nparr)->elsize-1);
      for (size_t i = 0; i < n; i++) {
        c = *a;
        *a++ = *b;
        *b-- = c;
      }
    }

  }

} // fragment swiglal_py_array_objview

// Name of fragment containing a NumPy object-view array descriptor for type ACFTYPE.
#define %swiglal_py_array_objview_frag(ACFTYPE) "swiglal_py_array_objview_" %str(ACFTYPE)

// Name of fragment containing NumPy object-view array descriptor initialisation code for type
// ACFTYPE.
#define %swiglal_py_array_objview_init_frag(ACFTYPE) "swiglal_py_array_objview_init_" %str(ACFTYPE)

// Macro which generates fragments containing a NumPy object-view array descriptor for type ACFTYPE.
//  - IN/OUTFRAG are names of fragments required by the in/out conversion functions IN/OUTCALL.
%define %swiglal_py_array_objview(ACFTYPE, INFRAG, OUTFRAG, INCALL, OUTCALL)

// Fragment containing NumPy object-view array descriptor initialisation code for type ACFTYPE.
%fragment(%swiglal_py_array_objview_init_frag(ACFTYPE), "init") {
  swiglal_py_array_objview_##ACFTYPE##_arrfuncs.cast[NPY_OBJECT] =
    (PyArray_VectorUnaryFunc*)swiglal_py_array_objview_##ACFTYPE##_cast_to_object;
}

// Fragment containing a NumPy object-view array descriptor for type ACFTYPE.
%fragment(%swiglal_py_array_objview_frag(ACFTYPE), "header",
          fragment="swiglal_py_array_objview",
          fragment=%swiglal_py_array_objview_init_frag(ACFTYPE),
          fragment=INFRAG, fragment=OUTFRAG)
{

  // NumPy array descriptor function which gets an element from the viewed array.
  static PyObject* swiglal_py_array_objview_##ACFTYPE##_getitem(void* elemptr, void* arr) {

    // Check input.
    assert(elemptr != NULL);
    assert(arr != NULL);
    PyArrayObject* nparr = (PyArrayObject*)arr;
    assert(PyArray_DESCR(nparr) != NULL);

    // Look up the SWIG type descriptor for this array.
    bool isptr;
    swig_type_info* tinfo = NULL;
    swiglal_py_array_tinfo_from_descr(&isptr, &tinfo, PyArray_DESCR(nparr));
    assert(tinfo != NULL);

    // Get the Python object wrapping the C array element.
    const bool copyobj = false;
    const size_t esize = PyArray_DESCR(nparr)->elsize;
    const int tflags = 0;
    PyObject* parent = PyArray_BASE(nparr);
    return OUTCALL;

  }

  // NumPy array descriptor function which assigns an element in the viewed array.
  static int swiglal_py_array_objview_##ACFTYPE##_setitem(PyObject* objelem, void* elemptr, void* arr) {

    // Check input.
    assert(elemptr != NULL);
    assert(arr != NULL);
    PyArrayObject* nparr = (PyArrayObject*)arr;
    assert(PyArray_DESCR(nparr) != NULL);

    // Look up the SWIG type descriptor for this array.
    bool isptr;
    swig_type_info* tinfo = NULL;
    swiglal_py_array_tinfo_from_descr(&isptr, &tinfo, PyArray_DESCR(nparr));
    assert(tinfo != NULL);

    // When assigning Python objects to a C array of pointers, assume the struct
    // who owns the C array takes ownership of the memory of the C array element.
    // The Python object wrapping the C array element should therefore disown the
    // underlying memory.
    // When assigning Python objects to a C array of data blocks, however, the C
    // array just struct-copies the object rather than taking ownership of its
    // pointer, and so the Python object should not be disowned so that it can
    // be garbage-collected later.
    const int tflags = isptr ? SWIG_POINTER_DISOWN : 0;

    // Set the C array element to the supplied Python object.
    const size_t esize = PyArray_DESCR(nparr)->elsize;
    PyObject* parent = PyArray_BASE(nparr);
    int elemalloc = 0;
    int *pelemalloc = &elemalloc;
    int res = INCALL;
    if (!SWIG_IsOK(res)) {
      SWIG_Error(res, "failure in swiglal_py_array_objview_" #ACFTYPE "_setitem()");
      return -1;
    }
    return 0;

  }

  // NumPy array descriptor function which casts elements of the viewed array to NPY_OBJECTs.
  static void swiglal_py_array_objview_##ACFTYPE##_cast_to_object(void *from, void *to, npy_intp n, void *fromarr, void *toarr) {

    // Check input.
    assert(fromarr != NULL);
    PyArrayObject* npfromarr = (PyArrayObject*)fromarr;
    assert(PyArray_DESCR(npfromarr) != NULL);
    assert(toarr != NULL);
    PyArrayObject* nptoarr = (PyArrayObject*)toarr;
    assert(PyArray_DESCR(nptoarr) != NULL);

    // 'toarr' should be an array of pointers to PyObjects.
    assert(PyArray_DESCR(nptoarr)->elsize == sizeof(PyObject*));

    // Loop over 'n' elements, and assign each element of 'toarr' the Python object wrapping the
    // corresponding element of 'fromarr'.
    char* fromelem = (void*)from;
    PyObject** toelem = (PyObject**)to;
    while (--n >= 0) {
      *toelem = swiglal_py_array_objview_##ACFTYPE##_getitem(fromelem, fromarr);
      fromelem += PyArray_DESCR(npfromarr)->elsize;
      ++toelem;
    }

  }

  // NumPy array descriptor function table for type ACFTYPE.
  static PyArray_ArrFuncs swiglal_py_array_objview_##ACFTYPE##_arrfuncs = {
    {(PyArray_VectorUnaryFunc*)NULL},   // cast
    (PyArray_GetItemFunc*)&swiglal_py_array_objview_##ACFTYPE##_getitem,   // getitem
    (PyArray_SetItemFunc*)&swiglal_py_array_objview_##ACFTYPE##_setitem,   // setitem
    (PyArray_CopySwapNFunc*)NULL,   // copyswapn
    (PyArray_CopySwapFunc*)&swiglal_py_array_objview_copyswap,   // copyswap
    (PyArray_CompareFunc*)NULL,   // compare
    (PyArray_ArgFunc*)NULL,   // argmax
    (PyArray_DotFunc*)NULL,   // dotfunc
    (PyArray_ScanFunc*)NULL,   // scanfunc
    (PyArray_FromStrFunc*)NULL,   // fromstr
    (PyArray_NonzeroFunc*)NULL,   // nonzero
    (PyArray_FillFunc*)NULL,   // fill
    (PyArray_FillWithScalarFunc*)NULL,   // fillwithscalar
    {(PyArray_SortFunc*)NULL},   // sort
    {(PyArray_ArgSortFunc*)NULL},   // argsort
    (PyObject*)NULL,   // castdict
    (PyArray_ScalarKindFunc*)NULL,   // scalarkind
    (int**)NULL,   // cancastscalarkindto
    (int*)swiglal_py_array_objview_copyswap_cancastto,   // cancastto
    (PyArray_FastClipFunc*)NULL,   // fastclip
    (PyArray_FastPutmaskFunc*)NULL,   // fastputmask
    (PyArray_FastTakeFunc*)NULL,   // fasttake
  };

  // This function returns the NumPy array descriptor appropriate for the supplied SWIG type
  // descriptor. If no array descriptor exists, it creates one from the array descriptor for type
  // ACFTYPE.
  SWIGINTERN PyArray_Descr* swiglal_py_array_objview_##ACFTYPE##_descr(const bool isptr, swig_type_info* tinfo, const int esize) {

    // Lookup existing NumPy array descriptor for SWIG type descriptor.
    PyArray_Descr* *pdescr = swiglal_py_array_descr_from_tinfo(isptr, tinfo);

    // Create NumPy array descriptor if none yet exists.
    if (*pdescr == NULL) {
      *pdescr = PyArray_DescrNewFromType(NPY_VOID);
      if (*pdescr == NULL) {
        return NULL;
      }
      (*pdescr)->typeobj = SwigPyObject_type();
      (*pdescr)->byteorder = '=';
      (*pdescr)->flags = NPY_LIST_PICKLE | NPY_NEEDS_INIT | NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM;
      (*pdescr)->type_num = 0;
      (*pdescr)->elsize = esize;
      (*pdescr)->alignment = 1;
      (*pdescr)->subarray = NULL;
      (*pdescr)->names = NULL;
      (*pdescr)->fields = NULL;
      (*pdescr)->f = &swiglal_py_array_objview_##ACFTYPE##_arrfuncs;

      if (PyArray_RegisterDataType(*pdescr) < 0) {
        return NULL;
      }
    }

    // PyArray_NewFromDescr appears to steal a reference to the descriptor passed to it, so a
    // reference count increment is needed here.
    Py_INCREF(*pdescr);

    return *pdescr;

  }

} // %swiglal_py_array_objview_frag(ACFTYPE)

%enddef // %swiglal_py_array_objview

// Macro which generates fragments which define ACFTYPE-specific array view classes and conversion
// functions:
//  - IN/OUTFRAG are names of fragments required by the in/out conversion functions IN/OUTCALL.
//  - VIEWFRAG is the name of a fragment needed for array views.
//  - NPYTYPE/NPYDESCR is the appropriate NumPy array typenum/descriptor.
%define %swiglal_py_array_frags(ACFTYPE, INFRAG, OUTFRAG, INCALL, OUTCALL, VIEWFRAG, NPYTYPE, NPYDESCR)

// Input copy conversion fragment for arrays of type ACFTYPE.
%fragment(%swiglal_array_copyin_frag(ACFTYPE), "header",
          fragment="swiglal_py_array_helpers", fragment=INFRAG)
{
  SWIGINTERN int %swiglal_array_copyin_func(ACFTYPE)(PyObject* parent,
                                                     PyObject* obj,
                                                     void* ptr,
                                                     int *pelemalloc,
                                                     const size_t esize,
                                                     const size_t ndims,
                                                     const size_t dims[],
                                                     const size_t strides[],
                                                     const bool isptr,
                                                     swig_type_info *tinfo,
                                                     const int tflags)
  {
    PyArrayObject* nparr = NULL;
    int res = 0;
    npy_intp idx[ndims];

    // Check that C array pointer is valid.
    if (ptr == NULL) {
      return SWIG_MemoryError;
    }

    // Convert the input Python object to a NumPy array.
    if (PyArray_Converter(obj, (PyObject**)&nparr) != NPY_SUCCEED) {
      return SWIG_ValueError;
    }

    // Check that NumPy array dimensions are consistent with C array dimensions.
    if (((size_t)PyArray_NDIM(nparr)) != ndims) {
      res = SWIG_ValueError;
      goto end;
    }
    size_t nelem = 1;
    for (size_t i = 0; i < ndims; ++i) {
      if (((size_t)PyArray_DIM(nparr, i)) != dims[i]) {
        res = SWIG_ValueError;
        goto end;
      }
      nelem *= dims[i];
    }

    // Iterate over all elements in the C array.
    memset(idx, 0, ndims*sizeof(npy_intp));
    for (size_t i = 0; i < nelem; ++i) {

      // Get a pointer to the element of the C array.
      void* elemptr = swiglal_py_get_element_ptr(ptr, esize, ndims, strides, idx);

      // Copy the NumPy array element to the C array.
      PyObject* objelem = PyArray_GETITEM(nparr, PyArray_GetPtr((PyArrayObject*)nparr, idx));
      res = INCALL;
      if (!SWIG_IsOK(res)) {
        goto end;
      }
      Py_CLEAR(objelem);

      // Increment the NumPy array index.
      swiglal_py_increment_idx(ndims, dims, idx);

    }

    res = SWIG_OK;

  end:
    Py_CLEAR(nparr);
    return res;

  }
}

// Output copy conversion fragment for arrays of type ACFTYPE.
%fragment(%swiglal_array_copyout_frag(ACFTYPE), "header",
          fragment="swiglal_py_array_helpers", fragment=OUTFRAG)
{
  SWIGINTERN PyObject* %swiglal_array_copyout_func(ACFTYPE)(PyObject* parent,
                                                            void* ptr,
                                                            const size_t esize,
                                                            const size_t ndims,
                                                            const size_t dims[],
                                                            const size_t strides[],
                                                            const bool isptr,
                                                            swig_type_info *tinfo,
                                                            const int tflags)
  {
    PyArrayObject* nparr = NULL;
    npy_intp objdims[ndims];
    npy_intp idx[ndims];

    // Check that C array pointer is valid.
    if (ptr == NULL) {
      goto fail;
    }

    // Copy C array dimensions.
    size_t nelem = 1;
    for (size_t i = 0; i < ndims; ++i) {
      objdims[i] = dims[i];
      nelem *= dims[i];
    }

    // Create new NumPy array.
    nparr = (PyArrayObject*)PyArray_EMPTY(ndims, objdims, NPYTYPE, 0);
    if (nparr == NULL) {
      goto fail;
    }

    // Iterate over all elements in the C array.
    memset(idx, 0, ndims*sizeof(npy_intp));
    for (size_t i = 0; i < nelem; ++i) {

      // Get a pointer to the element of the C array.
      void* elemptr = swiglal_py_get_element_ptr(ptr, esize, ndims, strides, idx);

      // Copy the C array element to the NumPy array.
      const bool copyobj = true;
      PyObject* objelem = OUTCALL;
      PyArray_SETITEM(nparr, PyArray_GetPtr((PyArrayObject*)nparr, idx), objelem);
      Py_CLEAR(objelem);

      // Increment the NumPy array index.
      swiglal_py_increment_idx(ndims, dims, idx);

    }

    return (PyObject*)nparr;

  fail:
    Py_CLEAR(nparr);
    Py_INCREF(Py_None);
    return Py_None;

  }
}

// Input view conversion fragment for arrays of type ACFTYPE.
%fragment(%swiglal_array_viewin_frag(ACFTYPE), "header",
          fragment="swiglal_py_array_helpers", fragment=INFRAG)
{
  SWIGINTERN int %swiglal_array_viewin_func(ACFTYPE)(PyObject* parent,
                                                     PyObject* obj,
                                                     void** ptr,
                                                     const size_t esize,
                                                     const size_t ndims,
                                                     size_t dims[],
                                                     const bool isptr,
                                                     swig_type_info *tinfo,
                                                     const int tflags)
  {
    PyArrayObject* nparr = NULL;
    int res = 0;

    // Check that C array pointer is valid.
    if (ptr == NULL) {
      return SWIG_MemoryError;
    }

    // Convert the input Python object to a NumPy array.
    if (PyArray_Converter(obj, (PyObject**)&nparr) != NPY_SUCCEED) {
      return SWIG_ValueError;
    }

    // Check that 'nparr' has the correct number of dimensions.
    if (((size_t)PyArray_NDIM(nparr)) != ndims) {
      res = SWIG_ValueError;
      goto end;
    }

    // Return dimensions of Python array.
    for (size_t i = 0; i < ndims; ++i) {
      dims[i] = PyArray_DIM(nparr, i);
    }

    // Cannot view an object which is not a NumPy array.
    if (!PyArray_Check(obj)) {
      res = SWIG_TypeError;
      goto end;
    }

    // Cannot view an array of pointers.
    if (isptr) {
      res = SWIG_TypeError;
      goto end;
    }

    // Cannot view an array of objects.
    if (NPYTYPE == NPY_OBJECT) {
      res = SWIG_TypeError;
      goto end;
    }

    // Cannot view an array which is not in C-array order.
    if (!PyArray_ISCARRAY(nparr)) {
      res = SWIG_TypeError;
      goto end;
    }

    // Check that 'nparr' is of the correct type.
    if (PyArray_TYPE(nparr) != NPYTYPE) {
      res = SWIG_TypeError;
      goto end;
    }

    // Check that the elements of 'nparr' have the correct size.
    if (((size_t)PyArray_ITEMSIZE(nparr)) != esize) {
      res = SWIG_TypeError;
      goto end;
    }

    // Get pointer to Python array data.
    *ptr = PyArray_DATA(nparr);
    if (*ptr == NULL) {
      res = SWIG_ValueError;
      goto end;
    }

    res = SWIG_OK;

  end:
    Py_CLEAR(nparr);
    return res;

  }
}

// Output view conversion fragment for arrays of type ACFTYPE.
%fragment(%swiglal_array_viewout_frag(ACFTYPE), "header",
          fragment="swiglal_py_array_helpers", fragment=VIEWFRAG, fragment=OUTFRAG)
{
  SWIGINTERN PyObject* %swiglal_array_viewout_func(ACFTYPE)(PyObject* parent,
                                                            void* ptr,
                                                            const size_t esize,
                                                            const size_t ndims,
                                                            const size_t dims[],
                                                            const size_t strides[],
                                                            const bool isptr,
                                                            swig_type_info *tinfo,
                                                            const int tflags)
  {
    PyArrayObject* nparr = NULL;
    npy_intp objdims[ndims];
    npy_intp objstrides[ndims];

    // Check that C array pointer is valid.
    if (ptr == NULL) {
      goto fail;
    }

    // Copy C array dimensions and strides.
    for (size_t i = 0; i < ndims; ++i) {
      objdims[i] = dims[i];
      objstrides[i] = strides[i] * esize;
    }

    // Create a new NumPy array view.
    PyArray_Descr* descr = NPYDESCR;
    if (descr == NULL) {
      goto fail;
    }
    nparr = (PyArrayObject*)PyArray_NewFromDescr(&PyArray_Type, descr, ndims, objdims, objstrides, ptr, NPY_ARRAY_WRITEABLE, NULL);
    if (nparr == NULL) {
      goto fail;
    }

    // Set the NumPy array view parent, if given.
    if (parent) {
      Py_INCREF(parent);
      PyArray_SetBaseObject(nparr, parent);
    }

    return (PyObject*)nparr;

  fail:
    Py_CLEAR(nparr);
    Py_INCREF(Py_None);
    return Py_None;

  }
}

%enddef // %swiglal_py_array_frags

// Macro which generates array conversion function fragments to/from Python arrays for object
// arrays, which require additional code for views.
%define %swiglal_py_array_objview_frags(ACFTYPE, INFRAG, OUTFRAG, INCALL, OUTCALL)
%swiglal_py_array_objview(ACFTYPE, INFRAG, OUTFRAG, INCALL, OUTCALL);
%swiglal_py_array_frags(ACFTYPE, INFRAG, OUTFRAG, INCALL, OUTCALL,
                        %swiglal_py_array_objview_frag(ACFTYPE),
                        NPY_OBJECT, %arg(swiglal_py_array_objview_##ACFTYPE##_descr(isptr, tinfo, esize)));
%enddef

// Array conversion fragments for generic arrays, e.g. SWIG-wrapped types.
%swiglal_py_array_objview_frags(SWIGTYPE, "swiglal_as_SWIGTYPE", "swiglal_from_SWIGTYPE",
                                %arg(swiglal_as_SWIGTYPE(parent, objelem, elemptr, esize, isptr, tinfo, tflags)),
                                %arg(swiglal_from_SWIGTYPE(parent, copyobj, elemptr, esize, isptr, tinfo, tflags)));

// Array conversion fragments for arrays of LAL strings.
%swiglal_py_array_objview_frags(LALchar, "SWIG_AsLALcharPtrAndSize", "SWIG_FromLALcharPtr",
                                %arg(SWIG_AsLALcharPtrAndSize(objelem, %reinterpret_cast(elemptr, char**), 0, pelemalloc)),
                                %arg(SWIG_FromLALcharPtr(*%reinterpret_cast(elemptr, char**))));

// Macro which generates array conversion function fragments to/from Python arrays for real/fragment
// TYPEs which use SWIG_AsVal/From fragments.
%define %swiglal_py_array_asvalfrom_frags(TYPE, NPYTYPE)
%swiglal_py_array_frags(TYPE, SWIG_AsVal_frag(TYPE), SWIG_From_frag(TYPE),
                        %arg(SWIG_AsVal(TYPE)(objelem, %reinterpret_cast(elemptr, TYPE*))),
                        %arg(SWIG_From(TYPE)(*%reinterpret_cast(elemptr, TYPE*))),
                        "swiglal_empty_frag", NPYTYPE, PyArray_DescrFromType(NPYTYPE));
%enddef

// Array conversion fragments for integer arrays.
%swiglal_py_array_asvalfrom_frags(int8_t, NPY_INT8);
%swiglal_py_array_asvalfrom_frags(uint8_t, NPY_UINT8);
%swiglal_py_array_asvalfrom_frags(int16_t, NPY_INT16);
%swiglal_py_array_asvalfrom_frags(uint16_t, NPY_UINT16);
%swiglal_py_array_asvalfrom_frags(int32_t, NPY_INT32);
%swiglal_py_array_asvalfrom_frags(uint32_t, NPY_UINT32);
%swiglal_py_array_asvalfrom_frags(int64_t, NPY_INT64);
%swiglal_py_array_asvalfrom_frags(uint64_t, NPY_UINT64);

// Array conversion fragments for floating-precision real arrays.
%swiglal_py_array_asvalfrom_frags(float, NPY_FLOAT);
%swiglal_py_array_asvalfrom_frags(double, NPY_DOUBLE);

// Array conversion fragments for floating-precision complex arrays.
%swiglal_py_array_asvalfrom_frags(gsl_complex_float, NPY_CFLOAT);
%swiglal_py_array_asvalfrom_frags(gsl_complex, NPY_CDOUBLE);
%swiglal_py_array_asvalfrom_frags(COMPLEX8, NPY_CFLOAT);
%swiglal_py_array_asvalfrom_frags(COMPLEX16, NPY_CDOUBLE);

// Local Variables:
// mode: c
// End:
