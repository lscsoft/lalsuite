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

// SWIG interface code specific to Python.
// Author: Karl Wette

////////// General SWIG directives and interface code //////////

// Include SWIG Python headers.
%include <pycomplex.swg>

// Include NumPy headers in wrapping code, and ensure that
// NumPy array module is loaded along with this module.
%header %{
#include <numpy/arrayobject.h>
%}
%init %{
  import_array();
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

// Append an argument to the output argument list of an Python SWIG-wrapped function, if the list is empty.
%header %{
#define swiglal_append_output_if_empty(v) \
  if (PySequence_Length(resultobj)) resultobj = SWIG_Python_AppendOutput(resultobj, v)
%}

////////// SWIG directives for operators //////////

// These macros apply the correct python:slot directives
// to map Python __operator__ functions (which may be
// defined in %extend) to the correct PyTypeObject slots.

// Unary operators which do not return a new object.
%define %swiglal_py_ury_op(NAME, FUNCTYPE, SLOT)
%pythonmaybecall *::__##NAME##__;
%feature("python:slot", #SLOT, functype=#FUNCTYPE) *::__##NAME##__;
%enddef
%swiglal_py_ury_op(float, unaryfunc, nb_float);
%swiglal_py_ury_op(hash, hashfunc, tp_hash);
%swiglal_py_ury_op(int, unaryfunc, nb_int);
%swiglal_py_ury_op(long, unaryfunc, nb_long);
%swiglal_py_ury_op(nonzero, inquiry, nb_nonzero);
%swiglal_py_ury_op(repr, reprfunc, tp_repr);
%swiglal_py_ury_op(str, reprfunc, tp_str);

// Unary operators which return a new object, and thus
// require %newobject to be set.
%define %swiglal_py_urn_op(NAME, FUNCTYPE, SLOT)
%newobject *::__##NAME##__;
%pythonmaybecall *::__##NAME##__;
%feature("python:slot", #SLOT, functype=#FUNCTYPE) *::__##NAME##__;
%enddef
%swiglal_py_urn_op(abs, unaryfunc, nb_absolute);
%swiglal_py_urn_op(neg, unaryfunc, nb_negative);
%swiglal_py_urn_op(pos, unaryfunc, nb_positive);

// Binary operators, which always must return a new object,
// and thus require %newobject to be set. The SWIG Python
// module with -builtin does not support reverse operators,
// so they are removed from the interface.
%define %swiglal_py_bin_op(NAME, FUNCTYPE, SLOT)
%newobject *::__##NAME##__;
%pythonmaybecall *::__##NAME##__;
%feature("python:slot", #SLOT, functype=#FUNCTYPE) *::__##NAME##__;
%ignore *::__r##NAME##__;
%enddef
%swiglal_py_bin_op(add, binaryfunc, nb_add);
%swiglal_py_bin_op(and, binaryfunc, nb_and);
%swiglal_py_bin_op(div, binaryfunc, nb_divide);
%swiglal_py_bin_op(lshift, binaryfunc, nb_lshift);
%swiglal_py_bin_op(mod, binaryfunc, nb_remainder);
%swiglal_py_bin_op(mul, binaryfunc, nb_multiply);
%swiglal_py_bin_op(or, binaryfunc, nb_or);
%swiglal_py_bin_op(rshift, binaryfunc, nb_rshift);
%swiglal_py_bin_op(sub, binaryfunc, nb_subtract);
%swiglal_py_bin_op(xor, binaryfunc, nb_xor);

// Comparison operators.
%define %swiglal_py_cmp_op(NAME, COMPTYPE)
%pythonmaybecall *::__##NAME##__;
%feature("python:compare", #COMPTYPE) *::__##NAME##__;
%enddef
%swiglal_py_cmp_op(eq, Py_EQ);
%swiglal_py_cmp_op(ge, Py_GE);
%swiglal_py_cmp_op(gt, Py_GT);
%swiglal_py_cmp_op(le, Py_LE);
%swiglal_py_cmp_op(lt, Py_LT);
%swiglal_py_cmp_op(ne, Py_NE);

////////// General fragments, typemaps, and macros //////////

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

// Typemaps which convert to/from the C broken-down date/time struct.
// Uses some code from PyLAL: Copyright (C) 2006 Kipp Cannon
%typemap(in) struct tm* (struct tm temptm) {

  // Set 'tm' struct to zero
  memset(&temptm, 0, sizeof(temptm));

  if ($input != NULL && $input != Py_None) {

    // Check that the $input PyObject is a sequence of either 6 or 9 integer elements; this
    // should also handle a time.struct_time (assuming its attributes are in the order below).
    // Note that the 7th ('tm_wday') and 8th ('tm_yday') elements are ignored; see below.
    if (!PySequence_Check($input)) {
      %argument_fail(SWIG_ValueError, "$type (not a sequence)", $symname, $argnum);
    }
    if (PySequence_Size($input) != 6 && PySequence_Size($input) != 9) {
      %argument_fail(SWIG_ValueError, "$type (must have 6 or 9 elements)", $symname, $argnum);
    }
    PyObject *seq = PySequence_Fast($input, "$type (not a sequence)");
    temptm.tm_year  = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 0)), int);
    temptm.tm_mon   = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 1)), int);
    temptm.tm_mday  = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 2)), int);
    temptm.tm_hour  = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 3)), int);
    temptm.tm_min   = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 4)), int);
    temptm.tm_sec   = %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 5)), int);
    temptm.tm_isdst = PySequence_Size($input) > 8 ?
      %static_cast(PyInt_AsLong(PySequence_Fast_GET_ITEM($input, 8)), int) : -1;
    Py_CLEAR(seq);
    if (PyErr_Occurred())   // Catch any errors while converting items to integers
      SWIG_fail;

    // Convert Python date ranges to 'tm' struct date ranges
    temptm.tm_year -= 1900;   // 'tm' struct years start from 1900
    temptm.tm_mon  -= 1;      // 'tm' struct months start from 0

    // Fill in values for 'tm_wday' and 'tm_yday', and normalise member ranges
    int errnum = 0;
    XLAL_TRY( XLALFillBrokenDownTime(&temptm), errnum );
    if (errnum != XLAL_SUCCESS) {
      %argument_fail(SWIG_ValueError, "$type (invalid date/time)", $symname, $argnum);
    }

  }

  $1 = &temptm;

}
%typemap(freearg) struct tm* "";
%typemap(out) struct tm* {

  // Convert 'tm' struct date ranges to Python date ranges
  $1->tm_year += 1900;                    // Python stores 4-digit years
  $1->tm_mon  += 1;                       // Python months start from 1
  $1->tm_wday  = ($1->tm_wday + 6) % 7;   // Python week days start from 0=Monday
  $1->tm_yday += 1;                       // Python year days start from 1

  // Build a 9-element list
  $result = Py_BuildValue("[iiiiiiiii]",
                          $1->tm_year, $1->tm_mon, $1->tm_mday,
                          $1->tm_hour, $1->tm_min, $1->tm_sec,
                          $1->tm_wday, $1->tm_yday, $1->tm_isdst);

}

////////// Interface code to track object parents //////////

// Interface code which tracks the parent structs of SWIG-wrapped struct members,
// so that the parent struct is not destroyed as long as a SWIG-wrapped object
// containing any of its members exists.
%header %{

  // Internal map from member pointers to octave_values containing the
  // member parent struct, as well as an internal reference count of
  // how many SWIG-wrapped member objects are extant.
  static PyObject *parent_map = NULL;

  // Store a reference to the parent of ptr in the internal map.
  // If there is already such a reference, increment the internal
  // reference count instead.
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

  // Check if ptr stored a reference to a parent struct. If there is
  // no parent object, then ptr *really* owns its memory, and it's okay
  // for it to destroy it (so return true). Otherwise, decrement the
  // internal reference count, erase the parent map entry if it reaches
  // zero, and return false to prevent any destructors being called.
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


%}
%init %{

  // Get a pointer to the internal parent map. Look for an attribute
  // 'parent_map' of an internal module 'swiglal_runtime_data'; if it
  // does not exist, create a new map and assign the module attribute,
  // otherwise store the attribute's value. In this way each wrapping
  // module gets a pointer to the same map.
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

%}

////////// Fragments and typemaps for arrays //////////

// This section implements array conversion functions for basic C array types,
// and custom NumPy array descriptors for viewing C arrays of object, e.g. structs.

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
    for (int j = 0; j < ndims; ++j) {
      elemidx += idx[j] * strides[j];
    }
    return %reinterpret_cast(%reinterpret_cast(ptr, char*) + elemidx*esize, void*);
  }

  // Increment the NumPy array index in row-major order, to match the ordering of the C array.
  void swiglal_py_increment_idx(const size_t ndims,
                                const size_t dims[],
                                npy_intp idx[])
  {
    for (int j = ndims-1; j >= 0; --j) {
      if (++idx[j] < dims[j]) {
        break;
      }
      idx[j] = 0;
    }
  }

 } // fragment swiglal_py_array_helpers

// Fragment defining helper functions for the NumPy object-view array descriptors.
%fragment("swiglal_py_array_objview", "header") {

  // Struct which associates a SWIG type descriptor with two NumPy array descriptors,
  // one for arrays of data blocks (_noptr), and one for arrays of pointers (_isptr).
  typedef struct {
    swig_type_info* tinfo;
    PyArray_Descr* descr_noptr;
    PyArray_Descr* descr_isptr;
  } swiglal_py_array_type_pair;

  // Static array of SWIG type/NumPy array descriptor pairs. This array should
  // always be long enough to accommodate all possible swig_type_info*, since
  // they are always members of the SWIG-generated global array swig_types[].
  // This array in turn is always one longer than the total number of types,
  // so there should always be a sentinal NULL element at the end.
  static swiglal_py_array_type_pair swiglal_py_array_types[sizeof(swig_types) / sizeof(swig_types[0])];

  // This function maps a SWIG type descriptor to a NumPy array descriptor,
  // or returns the first NULL element if a mapping doesn't exist yet.
  SWIGINTERN PyArray_Descr** swiglal_py_array_descr_from_tinfo(const bool isptr, swig_type_info* tinfo) {
    size_t i = 0;
    while (swiglal_py_array_types[i].tinfo != NULL && swiglal_py_array_types[i].tinfo != tinfo)
      ++i;
    if (swiglal_py_array_types[i].tinfo == NULL)
      swiglal_py_array_types[i].tinfo = tinfo;
    return isptr ? &swiglal_py_array_types[i].descr_isptr : &swiglal_py_array_types[i].descr_noptr;
  }

  // This function maps a NumPy array descriptor to a SWIG type descriptor,
  // or returns NULL element if a mapping doesn't exist.
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
    assert(nparr->descr != NULL);

    // Copy array element.
    if (src != NULL) {
      memcpy(dst, src, nparr->descr->elsize);
    }

    // Byte-swap array element, if required.
    if (swap) {
      const size_t n = nparr->descr->elsize / 2;
      char *a, *b, c;
      a = (char *)dst;
      b = a + (nparr->descr->elsize-1);
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

// Name of fragment containing NumPy object-view array descriptor initialisation code for type ACFTYPE.
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
    assert(nparr->descr != NULL);

    // Look up the SWIG type descriptor for this array.
    bool isptr;
    swig_type_info* tinfo = NULL;
    swiglal_py_array_tinfo_from_descr(&isptr, &tinfo, nparr->descr);
    assert(tinfo != NULL);

    // Get the Python object wrapping the C array element.
    const int tflags = 0;
    PyObject* parent = nparr->base;
    return OUTCALL;

  }

  // NumPy array descriptor function which assigns an element in the viewed array.
  static int swiglal_py_array_objview_##ACFTYPE##_setitem(PyObject* objelem, void* elemptr, void* arr) {

    // Check input.
    assert(elemptr != NULL);
    assert(arr != NULL);
    PyArrayObject* nparr = (PyArrayObject*)arr;
    assert(nparr->descr != NULL);

    // Look up the SWIG type descriptor for this array.
    bool isptr;
    swig_type_info* tinfo = NULL;
    swiglal_py_array_tinfo_from_descr(&isptr, &tinfo, nparr->descr);
    assert(tinfo != NULL);

    // Set the C array element to the supplied Python object.
    const int tflags = 0;
    const size_t esize = nparr->descr->elsize;
    PyObject* parent = nparr->base;
    int ecode = INCALL;
    if (!SWIG_IsOK(ecode)) {
      SWIG_Error(ecode, "failure in swiglal_py_array_objview_" #ACFTYPE "_setitem()");
      return -1;
    }
    return 0;

  }

  // NumPy array descriptor function which casts elements of the viewed array to NPY_OBJECTs.
  static void swiglal_py_array_objview_##ACFTYPE##_cast_to_object(void *from, void *to, npy_intp n, void *fromarr, void *toarr) {

    // Check input.
    assert(fromarr != NULL);
    PyArrayObject* npfromarr = (PyArrayObject*)fromarr;
    assert(npfromarr->descr != NULL);
    assert(toarr != NULL);
    PyArrayObject* nptoarr = (PyArrayObject*)toarr;
    assert(nptoarr->descr != NULL);

    // toarr should be an array of pointers to PyObjects.
    assert(nptoarr->descr->elsize == sizeof(PyObject*));

    // Loop over n elements, and assign each element of toarr
    // the Python object wrapping the corresponding element of fromarr.
    char* fromelem = (void*)from;
    PyObject** toelem = (PyObject**)to;
    while (--n >= 0) {
      *toelem = swiglal_py_array_objview_##ACFTYPE##_getitem(fromelem, fromarr);
      fromelem += npfromarr->descr->elsize;
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

  // NumPy array descriptor function for type ACFTYPE.
  static PyArray_Descr swiglal_py_array_objview_##ACFTYPE##_arrdescr = {
    PyObject_HEAD_INIT(NULL)
    (PyTypeObject*) NULL,   // typeobj
    NPY_VOIDLTR,   // kind
    NPY_VOIDLTR,   // type
    '=',   // byteorder
    NPY_LIST_PICKLE | NPY_USE_GETITEM | NPY_USE_SETITEM
    | NPY_ITEM_IS_POINTER | NPY_NEEDS_INIT | NPY_NEEDS_PYAPI,   // hasobject
    0,   // type_num
    0,   // elsize
    0,   // alignment
    (PyArray_ArrayDescr*)NULL,   // subarray
    (PyObject*)NULL,   // fields
    (PyObject*)NULL,   // names
    &swiglal_py_array_objview_##ACFTYPE##_arrfuncs,   // f
  };

  // This function returns the NumPy array descriptor appropriate for the
  // supplied SWIG type descriptor. If no array descriptor exists, it creates
  // one from the array descriptor for type ACFTYPE.
  SWIGINTERN PyArray_Descr* swiglal_py_array_objview_##ACFTYPE##_descr(const bool isptr, swig_type_info* tinfo, const int esize) {

    // Lookup existing NumPy array descriptor for SWIG type descriptor.
    PyArray_Descr* *pdescr = swiglal_py_array_descr_from_tinfo(isptr, tinfo);

    // Create NumPy array descriptor if none yet exists.
    if (*pdescr == NULL) {
      *pdescr = PyArray_DescrNew(&swiglal_py_array_objview_##ACFTYPE##_arrdescr);
      if (*pdescr == NULL) {
        return NULL;
      }
      (*pdescr)->typeobj = SwigPyObject_type();
      (*pdescr)->elsize = esize;
      (*pdescr)->alignment = 1;
      if (PyArray_RegisterDataType(*pdescr) < 0) {
        return NULL;
      }
    }

    // PyArray_NewFromDescr appears to steal a reference to the descriptor
    // passed to it, so a reference count increment is needed here.
    Py_INCREF(*pdescr);

    return *pdescr;

  }

} // %swiglal_py_array_objview_frag(ACFTYPE)

%enddef // %swiglal_py_array_objview

// Macro which generates fragments which define ACFTYPE-specific array view classes and conversion functions:
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
                                                       const size_t esize,
                                                       const size_t ndims,
                                                       const size_t dims[],
                                                       const size_t strides[],
                                                       const bool isptr,
                                                       swig_type_info *tinfo,
                                                       const int tflags)
    {
      int ecode = 0;
      npy_intp idx[ndims];

      // Check that C array is non-NULL.
      if (ptr == NULL) {
        return SWIG_MemoryError;
      }

      // Convert the input Python object to a NumPy array.
      PyObject* nparr = NULL;
      if (PyArray_Converter(obj, &nparr) != NPY_SUCCEED) {
        ecode = SWIG_ValueError;
        goto end;
      }

      // Check that NumPy array dimensions are consistent with C array dimensions.
      if (PyArray_NDIM(nparr) != ndims) {
        ecode = SWIG_ValueError;
        goto end;
      }
      size_t nelem = 1;
      for (int i = 0; i < ndims; ++i) {
        if (PyArray_DIM(nparr, i) != dims[i]) {
          ecode = SWIG_ValueError;
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
        ecode = INCALL;
        if (!SWIG_IsOK(ecode)) {
          goto end;
        }
        Py_CLEAR(objelem);

        // Increment the NumPy array index.
        swiglal_py_increment_idx(ndims, dims, idx);

      }

      ecode = SWIG_OK;

    end:
      Py_CLEAR(nparr);
      return ecode;

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
      PyObject* nparr = NULL;
      npy_intp objdims[ndims];
      npy_intp idx[ndims];

      // Check that C array is non-NULL.
      if (ptr == NULL) {
        goto fail;
      }

      // Copy C array dimensions.
      size_t nelem = 1;
      for (int i = 0; i < ndims; ++i) {
        objdims[i] = dims[i];
        nelem *= dims[i];
      }

      // Create new NumPy array.
      nparr = PyArray_EMPTY(ndims, objdims, NPYTYPE, 0);
      if (nparr == NULL) {
        goto fail;
      }

      // Iterate over all elements in the C array.
      memset(idx, 0, ndims*sizeof(npy_intp));
      for (size_t i = 0; i < nelem; ++i) {

        // Get a pointer to the element of the C array.
        void* elemptr = swiglal_py_get_element_ptr(ptr, esize, ndims, strides, idx);

        // Copy the C array element to the NumPy array.
        PyObject* objelem = OUTCALL;
        PyArray_SETITEM(nparr, PyArray_GetPtr((PyArrayObject*)nparr, idx), objelem);

        // Increment the NumPy array index.
        swiglal_py_increment_idx(ndims, dims, idx);

      }

      return nparr;

    fail:
      Py_CLEAR(nparr);
      Py_INCREF(Py_None);
      return Py_None;

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
      PyObject* nparr = NULL;
      npy_intp objdims[ndims];
      npy_intp objstrides[ndims];

      // Copy C array dimensions and strides.
      for (int i = 0; i < ndims; ++i) {
        objdims[i] = dims[i];
        objstrides[i] = strides[i] * esize;
      }

      // Create a new NumPy array view.
      PyArray_Descr* descr = NPYDESCR;
      if (descr == NULL) {
        goto fail;
      }
      nparr = PyArray_NewFromDescr(&PyArray_Type, descr, ndims, objdims, objstrides, ptr, NPY_WRITEABLE, NULL);
      if (nparr == NULL) {
        goto fail;
      }

      // Set the NumPy array view parent, if given.
      if (parent) {
        Py_INCREF(parent);
        ((PyArrayObject*)nparr)->base = parent;
      }

      return nparr;

    fail:
      Py_CLEAR(nparr);
      Py_INCREF(Py_None);
      return Py_None;

    }
  }

%enddef // %swiglal_py_array_frags

// Macro which generates array conversion function fragments to/from Python
// arrays for object arrays, which require additional code for views.
%define %swiglal_py_array_objview_frags(ACFTYPE, INFRAG, OUTFRAG, INCALL, OUTCALL)
%swiglal_py_array_objview(ACFTYPE, INFRAG, OUTFRAG, INCALL, OUTCALL);
%swiglal_py_array_frags(ACFTYPE, INFRAG, OUTFRAG, INCALL, OUTCALL,
                        %swiglal_py_array_objview_frag(ACFTYPE),
                        NPY_OBJECT, %arg(swiglal_py_array_objview_##ACFTYPE##_descr(isptr, tinfo, esize)));
%enddef

// Array conversion fragments for generic arrays, e.g. SWIG-wrapped types.
%swiglal_py_array_objview_frags(SWIGTYPE, "swiglal_as_SWIGTYPE", "swiglal_from_SWIGTYPE",
                                %arg(swiglal_as_SWIGTYPE(parent, objelem, elemptr, esize, isptr, tinfo, tflags)),
                                %arg(swiglal_from_SWIGTYPE(parent, elemptr, isptr, tinfo, tflags)));

// Array conversion fragments for arrays of LAL strings.
%swiglal_py_array_objview_frags(LALchar, "SWIG_AsNewLALcharPtr", "SWIG_FromLALcharPtr",
                                %arg(SWIG_AsNewLALcharPtr(objelem, %reinterpret_cast(elemptr, char**))),
                                %arg(SWIG_FromLALcharPtr(*%reinterpret_cast(elemptr, char**))));

// Macro which generates array conversion function fragments to/from Python
// arrays for real/fragment TYPEs which use SWIG_AsVal/From fragments.
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
