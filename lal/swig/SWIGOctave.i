//
// Copyright (C) 2011--2017, 2022 Karl Wette
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

// SWIG interface code specific to Octave.
// Author: Karl Wette

//
// General SWIG directives and interface code
//

// In SWIG Octave modules, only variables are namespaced, everything else
// is inserted in the global symbol table, so we rename only variables
#define SWIGLAL_MODULE_RENAME_VARIABLES

// Include SWIG Octave headers.
%include <octcomplex.swg>

// Include Octave headers and set versioning macros.
%header %{
extern "C++" {
#include <octave/ov-cell.h>
#include <octave/ov-int-traits.h>
#include <octave/ov-flt-re-mat.h>
#include <octave/ov-re-mat.h>
#include <octave/ov-flt-cx-mat.h>
#include <octave/ov-cx-mat.h>
#include <octave/Array-util.h>
}
%}

// Evaluates true if an octave_value is not empty, false otherwise.
%header %{
#define swiglal_not_empty(v)  (!(v).is_empty())
%}

// Name of octave_value containing the SWIG wrapping of the struct whose members are being accessed.
%header %{
#define swiglal_self()    (args.length() > 0 ? args(0) : octave_value())
#define swiglal_no_self() (octave_value())
%}

// Name of octave_value containing the SWIG wrapping of the first argument to a function.
%header %{
#define swiglal_1starg()  (args.length() > 0 ? args(0) : octave_value())
%}

// Return a reference to the supplied octave_value; since Octave handles reference counting, just return it.
%header %{
#define swiglal_get_reference(v) (v)
%}

// Append an argument to the output argument list of an Octave SWIG-wrapped function, if the list is empty.
%header %{
#define swiglal_append_output_if_empty(v) if (_outp->length() == 0) _outp = SWIG_Octave_AppendOutput(_outp, v)
%}

// Evaluates true if an octave_value represents a null pointer, false otherwise.
%header %{
#define swiglal_null_ptr(v)  (!(v).is_string() && (v).is_matrix_type() && (v).rows() == 0 && (v).columns() == 0)
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
    octave_value_list args = feval("stdout", octave_value_list(), 1);
    if (args.length() < 1) {
      return 0;
    }
    args.resize(3);
    args(1) = octave_value(std::string("%s"));
    char buf[512];
    while (fgets(buf, sizeof(buf), swiglal_tmp_stdout) != NULL) {
      args(2) = octave_value(std::string(buf));
      feval("fprintf", args, 0);
    }
  }

  // Write standard error
  {
    octave_value_list args = feval("stderr", octave_value_list(), 1);
    if (args.length() < 1) {
      return 0;
    }
    args.resize(3);
    args(1) = octave_value(std::string("%s"));
    char buf[512];
    while (fgets(buf, sizeof(buf), swiglal_tmp_stderr) != NULL) {
      args(2) = octave_value(std::string(buf));
      feval("fprintf", args, 0);
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

// Unary operators which return a new object, and thus require %newobject to be set.
%define %swiglal_oct_urn_op(NAME, OCTNAME)
%rename(__##OCTNAME##__) *::__##NAME##__;
%newobject *::__##NAME##__;
%enddef
%swiglal_oct_urn_op(abs, abs);
%swiglal_oct_urn_op(neg, uminus);
%swiglal_oct_urn_op(pos, uplus);

// Binary operators, which always must return a new object, and thus require %newobject to be set.
%define %swiglal_oct_bin_op(NAME)
%newobject *::__##NAME##__;
%newobject *::__r##NAME##__;
%enddef
%swiglal_oct_bin_op(add);
%swiglal_oct_bin_op(and);
%swiglal_oct_bin_op(div);
%swiglal_oct_bin_op(lshift);
%swiglal_oct_bin_op(mod);
%swiglal_oct_bin_op(mul);
%swiglal_oct_bin_op(or);
%swiglal_oct_bin_op(pow);
%swiglal_oct_bin_op(rshift);
%swiglal_oct_bin_op(sub);
%swiglal_oct_bin_op(xor);

// Octave __pow__() operator takes 2 arguments, so we ignore the 3rd.
%typemap(in, numinputs=0) void* SWIGLAL_OP_POW_3RDARG "";

// Comparison operators.
%typemap(in, numinputs=0, noblock=1) int SWIGLAL_CMP_OP_RETN_HACK "";

//
// Octave-specific extensions to structs
//

// Extend a struct TAGNAME.
%define %swiglal_struct_extend_specific(TAGNAME, OPAQUE, DTORFUNC)

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

// Typemaps which convert to/from the C broken-down date/time struct.
%typemap(in, fragment=SWIG_AsVal_frag(double)) struct tm* (struct tm temptm) {

  // Convert '$input' to an Octave date number 'datenum'
  octave_value_list datenum_args;
  if ($input.is_string()) {
    datenum_args.append($input);
  } else {
    dim_vector dims = $input.dims();
    if (dims.length() == 2 && dims.num_ones() == 1 && 3 <= dims.numel() && dims.numel() <= 6) {
      RowVector datevec = $input.row_vector_value();
      for (int i = 0; i < datevec.numel(); ++i) {
        datenum_args.append(octave_value(datevec(i)));
      }
    }
  }
  octave_value_list retn;
  if (datenum_args.length() > 0) {
    retn = feval("datenum", datenum_args, 1);
  }
  if (retn.length() == 0) {
    %argument_fail(SWIG_ValueError, "$type", $symname, $argnum);
  }
  double datenum = 0;
  int res = SWIG_AsVal(double)(retn(0), &datenum);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "$type", $symname, $argnum);
  }

  // Convert 'datenum' to a C time_t 't', using the following:
  //   1970-01-01 00:00:00 +0000: 'datenum'=719529, 't'=0
  //   1970-01-02 00:00:00 +0000: 'datenum'=719530, 't'=86400
  time_t t = (time_t) llround((datenum - 719529) * 86400);

  // Convert 't' to a C struct tm 'temptm'
  memset(&temptm, 0, sizeof(temptm));
  if (gmtime_r(&t, &temptm) == NULL) {
    %argument_fail(SWIG_RuntimeError, "$type", $symname, $argnum);
  }

  $1 = &temptm;

}
%typemap(freearg) struct tm* "";
%typemap(out) struct tm* {

  // Create an Octave date vector 'datevec'
  RowVector datevec(6);

  // Assign members 'datevec', converting 'tm' struct date ranges to Octave date ranges
  datevec(0) = $1->tm_year + 1900;   // Octave stores 4-digit years
  datevec(1) = $1->tm_mon  + 1;      // Octave months start from 1
  datevec(2) = $1->tm_mday;
  datevec(3) = $1->tm_hour;
  datevec(4) = $1->tm_min;
  datevec(5) = $1->tm_sec;

  $result = octave_value(datevec);

}

//
// Interface code to track object parents
//

// Interface code which tracks the parent structs of SWIG-wrapped struct members, so that the parent
// struct is not destroyed as long as a SWIG-wrapped object containing any of its members exists.
%header %{

#include <map>

// Internal map from member pointers to octave_values containing the member parent struct, as well
// as an internal reference count of how many SWIG-wrapped member objects are extant.
typedef std::pair<octave_value, int> swiglal_oct_parent;
typedef std::pair<void*, swiglal_oct_parent> swiglal_oct_parent_pair;
typedef std::map<void*, swiglal_oct_parent> swiglal_oct_parent_map;
static swiglal_oct_parent_map* parent_map = 0;

// Store a reference to the parent of ptr in the internal map. If there is already such a reference,
// increment the internal reference count instead.
SWIGINTERN void swiglal_store_parent(void* ptr, octave_value parent) {
  assert(ptr);
  assert(parent.is_defined());
  swiglal_oct_parent_map::iterator i = parent_map->find(ptr);
  if (i == parent_map->end()) {
    parent_map->insert(swiglal_oct_parent_pair(ptr, swiglal_oct_parent(parent, 1)));
  }
  else {
    ++i->second.second;
  }
}

// Check if ptr stored a reference to a parent struct. If there is no parent object, then ptr
// *really* owns its memory, and it's okay for it to destroy it (so return true). Otherwise,
// decrement the internal reference count, erase the parent map entry if it reaches zero, and return
// false to prevent any destructors being called.
SWIGINTERN bool swiglal_release_parent(void *ptr) {
  bool retn = true;
  assert(ptr);
  swiglal_oct_parent_map::iterator i = parent_map->find(ptr);
  if (i != parent_map->end()) {
    retn = false;
    if (--i->second.second == 0) {
      parent_map->erase(i);
    }
  }
  return retn;
}

%} // %header
%init %{

// Get a pointer to the internal parent map. Look for a global variable ' __SWIGLAL_parent_map__',
// then try to extract a pointer to the parent map from a SWIG packed object. If the packed object
// does not yet exist, create a new map, pack its pointer in a SWIG packed object, then assign it to
// global variable. Thus each binding gets a pointer to the same map.
{
  const char* const parent_map_name = "__SWIGLAL_parent_map__";
  octave_value ov = SWIG_Octave_GetGlobalValue(parent_map_name);
  int res = SWIG_ConvertPacked(ov, &parent_map, sizeof(parent_map), 0);
  if (!SWIG_IsOK(res)) {
    parent_map = new swiglal_oct_parent_map();
    ov = SWIG_NewPackedObj(&parent_map, sizeof(parent_map), 0);
    SWIG_Octave_SetGlobalValue(parent_map_name, ov);
  }
  assert(parent_map);
}

%} // %init

//
// Fragments and typemaps for arrays
//

// This section implements a series of array view classes, through which arbitrary C array data can
// be viewed as native Octave matrices, etc.

// Name of array view class for array conversion type ACFTYPE.
#define %swiglal_oct_array_view_class(ACFTYPE) swiglal_oct_array_view_##ACFTYPE

// Name of helper class for array view class for array conversion type ACFTYPE.
#define %swiglal_oct_array_view_helper_class(ACFTYPE) swiglal_oct_array_view_helper_##ACFTYPE

// Name of base array view template for array conversion type ACFTYPE.
#define %swiglal_oct_array_view_tmpl(ACFTYPE) swiglal_oct_array_view<%swiglal_oct_array_view_helper_class(ACFTYPE) >

// String denoting octave_value type of array view class for array conversion type ACFTYPE.
#define %swiglal_oct_array_view_ovtype(ACFTYPE) "swiglal_oct_array_view_" %str(ACFTYPE)

// Name of fragment containing array view class for array conversion type ACFTYPE.
#define %swiglal_oct_array_view_frag(ACFTYPE) "swiglal_oct_array_view_" %str(ACFTYPE)

// Name of fragment containing initialisation code for array view class for array conversion type
// ACFTYPE.
#define %swiglal_oct_array_view_init_frag(ACFTYPE) "swiglal_oct_array_view_init_" %str(ACFTYPE)

// Fragment defining a base array view template, where all functionality is implemented, and from
// which ACFTYPE-specific array view classes inherit. The template argument is a helper class which
// supplied ACFTYPE-specific functions and types.
%fragment("swiglal_oct_array_view", "header") {

  extern "C++" {

    template<class HELPER>
    class swiglal_oct_array_view : public octave_base_value {

    private:

      // Instance of the corresponding octave_base_value-derived class of the array view class, for
      // consulting as to various properties.
      const typename HELPER::OVClass sloav_class;

      // Keep a reference to the SWIG-wrapped struct containing the C array being viewed, to prevent
      // it being destroyed if the struct goes out of scope by the array view remains.
      const octave_value sloav_parent;

      // Parameters of the C array data being viewed, and associated SWIG type information.
      void *const sloav_ptr;
      const size_t sloav_esize;
      const size_t sloav_ndims;
      const dim_vector sloav_dims;
      const dim_vector sloav_strides;
      const bool sloav_isptr;
      swig_type_info *const sloav_tinfo;
      const int sloav_tflags;

      // Construct an Octave dim_vector from a C array.
      static dim_vector sloav_make_dim_vector(const size_t n, const size_t v[]) {
        dim_vector dv(n, 1);
        for (size_t i = 0; i < n; ++i) {
          dv(i) = v[i];
        }
        return dv;
      }

      // Numeric conversion function for converting an array view into an Octave array.
      static octave_base_value* sloav_numeric_conversion_function(const octave_base_value& v) {
        const swiglal_oct_array_view& oav = dynamic_cast<const swiglal_oct_array_view&>(v);
        octave_value ov = oav.sloav_array_out();
        return new typename HELPER::OVClass(HELPER::ovvalue(ov));
      }

      // Compute the scalar index of the C array element, and return a pointer to the element itself.
      void* sloav_get_element_ptr(Array<octave_idx_type>& idx) const {
        size_t elemidx = 0;
        for (size_t j = 0; j < sloav_ndims; ++j) {
          elemidx += idx(j) * sloav_strides(j);
        }
        return reinterpret_cast<void*>(reinterpret_cast<char*>(sloav_ptr) + elemidx*sloav_esize);
      }

      // Increment the Octave array index in row-major order, to match the ordering of the C array.
      void sloav_increment_idx(Array<octave_idx_type>& idx) const {
        for (int j = sloav_ndims-1; j >= 0; --j) {
          if (++idx(j) < sloav_dims(j)) {
            break;
          }
          idx(j) = 0;
        }
      }

    public:

      virtual ~swiglal_oct_array_view()
      { }

      swiglal_oct_array_view()
        : octave_base_value(), sloav_class(typename HELPER::OVClass()),
          sloav_parent(), sloav_ptr(0), sloav_esize(0), sloav_ndims(0),
          sloav_dims(), sloav_strides(),
          sloav_isptr(false), sloav_tinfo(0), sloav_tflags(0)
      { }

      swiglal_oct_array_view(const octave_value& parent,
                             void* ptr,
                             const size_t esize,
                             const size_t ndims,
                             const size_t dims[],
                             const size_t strides[],
                             const bool isptr,
                             swig_type_info* tinfo,
                             const int tflags)
        : octave_base_value(), sloav_class(typename HELPER::OVClass()),
          sloav_parent(parent), sloav_ptr(ptr), sloav_esize(esize), sloav_ndims(ndims),
          sloav_dims(sloav_make_dim_vector(ndims, dims)),
          sloav_strides(sloav_make_dim_vector(ndims, strides)),
          sloav_isptr(isptr), sloav_tinfo(tinfo), sloav_tflags(tflags)
      { }

      // Copy the Octave array obj to the C array.
      int sloav_array_in(octave_value& obj, int *pelemalloc, const int tflags) {

        // Check that C array pointer is valid.
        if (!sloav_ptr) {
          return SWIG_MemoryError;
        }

        // Check that Octave array dimensions are consistent with C array dimensions.  1-D arrays
        // are a special case, since Octave arrays are always at least 2-dimensional, so need to
        // check that one of those dimensions is singular.
        dim_vector objdims = obj.dims();
        if (sloav_ndims == 1) {
          if (objdims.length() > 2 || objdims.num_ones() == 0 || objdims.numel() != sloav_dims(0)) {
            return SWIG_ValueError;
          }
        }
        else if (objdims != sloav_dims) {
          return SWIG_ValueError;
        }

        // Iterate over all elements in the C array.
        Array<octave_idx_type> idx(dim_vector(1, sloav_ndims), 0);
        std::list<octave_value_list> objidx(1);
        for (int i = 0; i < objdims.numel(); ++i) {

          // Get the scalar index of the Octave array element, and the element itself.
          objidx.front()(0) = get_scalar_idx(idx, objdims) + 1;
          octave_value objelem = obj.subsref(obj.is_cell() ? "{" : "(", objidx);

          // Copy the Octave array element to the C array.
          int res = HELPER::incall(sloav_parent, objelem, sloav_get_element_ptr(idx), pelemalloc, sloav_esize, sloav_isptr, sloav_tinfo, sloav_tflags | tflags);
          if (!SWIG_IsOK(res)) {
            return res;
          }

          // Increment the Octave array index.
          sloav_increment_idx(idx);

        }

        return SWIG_OK;

      }

      // Copy the C array to the returned Octave array.
      octave_value sloav_array_out(const bool copyobj = false) const {

        // Check that C array pointer is valid.
        if (!sloav_ptr) {
          return octave_value();
        }

        // Create a new Octave array.
        dim_vector objdims = sloav_dims;
        typename HELPER::OVType objval(objdims);
        octave_value obj(objval);

        // Iterate over all elements in the C array.
        Array<octave_idx_type> idx(dim_vector(1, sloav_ndims), 0);
        std::list<octave_value_list> objidx(1);
        for (int i = 0; i < objdims.numel(); ++i) {

          // Get the scalar index of the Octave array element.
          objidx.front()(0) = get_scalar_idx(idx, objdims) + 1;

          // Copy the C array element to the Octave array.
          octave_value objelem = HELPER::outcall(sloav_parent, copyobj, sloav_get_element_ptr(idx), sloav_esize, sloav_isptr, sloav_tinfo, sloav_tflags);
          obj = obj.subsasgn(obj.is_cell() ? "{" : "(", objidx, objelem);

          // Increment the Octave array index.
          sloav_increment_idx(idx);

        }

        return obj;

      }

      // The following methods override virtual methods in octave_base_value.

      // Array dimensions.
      dim_vector dims() const {
        return sloav_dims;
      }

      // Return an Octave array.
      octave_value full_value() const {
        return sloav_array_out();
      }

      // Return an empty Octave array.
      octave_base_value* empty_clone() const {
        return sloav_array_out().empty_clone();
      }

      // Return the numeric conversion function.
      octave_base_value::type_conv_info numeric_conversion_function() const {
        return octave_base_value::type_conv_info(sloav_numeric_conversion_function, sloav_class.type_id());
      }

      // Do indexing without resizing.
      octave_value do_index_op(const octave_value_list& idx, bool resize_ok) {
        return sloav_array_out().do_index_op(idx, false);
      }
      octave_value_list do_multi_index_op(int nargout, const octave_value_list& idx) {
        return sloav_array_out().do_multi_index_op(nargout, idx);
      }

      // Do subscripting on Octave array.
      octave_value subsref(const std::string& type, const std::list<octave_value_list>& idx) {
        return sloav_array_out().subsref(type, idx);
      }
      octave_value_list subsref(const std::string& type, const std::list<octave_value_list>& idx, int nargout) {
        return sloav_array_out().subsref(type, idx, nargout);
      }

      // Do subscript assignment, and copy result back to C array.
      octave_value subsasgn(const std::string& type, const std::list<octave_value_list>& idx, const octave_value& rhs) {
        octave_value obj = sloav_array_out().subsasgn(type, idx, rhs);
        int elemalloc = 0;
        // When assigning Octave objects to a C array of pointers, assume the struct
        // who owns the C array takes ownership of the memory of the C array element.
        // The Octave object wrapping the C array element should therefore disown the
        // underlying memory.
        // When assigning Octave objects to a C array of data blocks, however, the C
        // array just struct-copies the object rather than taking ownership of its
        // pointer, and so the Octave object should not be disowned so that it can
        // be garbage-collected later.
        int res = sloav_array_in(obj, &elemalloc, sloav_isptr ? SWIG_POINTER_DISOWN : 0);
        if (!SWIG_IsOK(res)) {
          std::string n = type_name();
          std::string e = SWIG_ErrorType(res).string_value();
          error("failed to perform indexed assignment for %s type: %s", n.c_str(), e.c_str());
          return octave_value();
        }
        count++;
        return octave_value(this);
      }

      // Save and load from ASCII.
      bool save_ascii (std::ostream& os) {
        return sloav_array_out().save_ascii(os);
      }
      bool load_ascii(std::istream& is) {
        octave_value obj = sloav_array_out();
        int elemalloc = 0;
        return obj.load_ascii(is) && SWIG_IsOK(sloav_array_in(obj, &elemalloc, 0));
      }

      // Save and load from binary.
      bool save_binary(std::ostream& os, bool& save_as_floats) {
        return sloav_array_out().save_binary(os, save_as_floats);
      }
      bool load_binary(std::istream& is, bool swap, oct_mach_info::float_format fmt) {
        octave_value obj = sloav_array_out();
        int elemalloc = 0;
        return obj.load_binary(is, swap, fmt) && SWIG_IsOK(sloav_array_in(obj, &elemalloc, 0));
      }

      // Save and load from HDF5.
%#ifdef HAVE_HDF5
%#if SWIG_OCTAVE_PREREQ(4,0,0)
      bool save_hdf5(octave_hdf5_id loc_id, const char *name, bool save_as_floats) {
        return sloav_array_out().save_hdf5(loc_id, name, save_as_floats);
      }
%#else
      bool save_hdf5(hid_t loc_id, const char *name, bool save_as_floats) {
        return sloav_array_out().save_hdf5(loc_id, name, save_as_floats);
      }
%#endif
%#if SWIG_OCTAVE_PREREQ(4,0,0)
      bool load_hdf5(octave_hdf5_id loc_id, const char *name) {
        octave_value obj = sloav_array_out();
        int elemalloc = 0;
        return obj.load_hdf5(loc_id, name) && SWIG_IsOK(sloav_array_in(obj, &elemalloc, 0));
      }
%#elif SWIG_OCTAVE_PREREQ(3,3,52)
      bool load_hdf5(hid_t loc_id, const char *name) {
        octave_value obj = sloav_array_out();
        int elemalloc = 0;
        return obj.load_hdf5(loc_id, name) && SWIG_IsOK(sloav_array_in(obj, &elemalloc, 0));
      }
%#else
      bool load_hdf5(hid_t loc_id, const char *name, bool have_h5giterate_bug) {
        octave_value obj = sloav_array_out();
        int elemalloc = 0;
        return obj.load_hdf5(loc_id, name, have_h5giterate_bug) && SWIG_IsOK(sloav_array_in(obj, &elemalloc, 0));
      }
%#endif
%#endif

      // Print array.
%#if SWIG_OCTAVE_PREREQ(4,0,0)
      void print(std::ostream &os, bool pr_as_read_syntax = false)
%#else
      void print(std::ostream &os, bool pr_as_read_syntax = false) const
%#endif
      {
        return sloav_array_out().print(os, pr_as_read_syntax);
      }

      // The following methods override virtual const-methods in octave_base_value.  These methods
      // are mapped to the equivalent method of the octave_base_value-derived class of the array
      // view class.

#define SLOAV_OBV_METH_FROM_CLASS_0(N, R) R N() const { return sloav_class.N(); }
      SLOAV_OBV_METH_FROM_CLASS_0(is_all_va_args, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_bool_matrix, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_bool_scalar, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_bool_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_builtin_function, bool);
%#if SWIG_OCTAVE_PREREQ(4,4,0)
      SLOAV_OBV_METH_FROM_CLASS_0(iscell, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(iscellstr, bool);
%#else
      SLOAV_OBV_METH_FROM_CLASS_0(is_cell, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_cellstr, bool);
%#endif
      SLOAV_OBV_METH_FROM_CLASS_0(is_char_matrix, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_classdef_meta, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_classdef_object, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_classdef_superclass_ref, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_complex_matrix, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_complex_scalar, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_complex_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_constant, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_cs_list, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_defined, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_diag_matrix, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_dld_function, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_double_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_float_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_function, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_function_handle, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_inline_function, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_int16_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_int32_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_int64_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_int8_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_integer_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_list, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_magic_colon, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_map, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_matrix_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_mex_function, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_null_value, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_numeric_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_object, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_perm_matrix, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_range, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_real_matrix, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_real_nd_array, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_real_scalar, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_real_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_scalar_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_single_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_sparse_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_sq_string, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_string, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_true, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_uint16_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_uint32_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_uint64_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_uint8_type, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_user_code, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_user_function, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(is_user_script, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(isjava, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(isobject, bool);
      SLOAV_OBV_METH_FROM_CLASS_0(isstruct, bool);
#undef SLOAV_OBV_METH_FROM_CLASS_0

      // The following methods override virtual const-methods in octave_base_value.  These methods
      // are mapped to the equivalent method of the Octave array.

#define SLOAV_OBV_METH_FROM_ARRAY_0(N, R) R N() const { return sloav_array_out().N(); }
#define SLOAV_OBV_METH_FROM_ARRAY_1(N, R, A) R N(A a) const { return sloav_array_out().N(a); }
#define SLOAV_OBV_METH_FROM_ARRAY_2(N, R, A, B) R N(A a, B b) const { return sloav_array_out().N(a, b); }
#define SLOAV_OBV_METH_FROM_ARRAY_3(N, R, A, B, C) R N(A a, B b, C c) const { return sloav_array_out().N(a, b, c); }
#define SLOAV_OBV_METH_FROM_ARRAY_4(N, R, A, B, C, D) R N(A a, B b, C c, D d) const { return sloav_array_out().N(a, b, c, d); }
#define SLOAV_OBV_METH_FROM_ARRAY_5(N, R, A, B, C, D, E) R N(A a, B b, C c, D d, E e) const { return sloav_array_out().N(a, b, c, d, e); }
%#if SWIG_OCTAVE_PREREQ(4,2,0)
      SLOAV_OBV_METH_FROM_ARRAY_0(as_double, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(as_single, octave_value);
%#endif
%#if SWIG_OCTAVE_PREREQ(3,3,52)
      SLOAV_OBV_METH_FROM_ARRAY_0(map_value, octave_map);
%#else
      SLOAV_OBV_METH_FROM_ARRAY_0(map_value, Octave_map);
%#endif
      SLOAV_OBV_METH_FROM_ARRAY_0(abs, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(acos, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(acosh, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(angle, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(arg, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(asin, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(asinh, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(atan, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(atanh, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(byte_size, size_t);
      SLOAV_OBV_METH_FROM_ARRAY_0(ceil, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(cell_value, Cell);
      SLOAV_OBV_METH_FROM_ARRAY_0(cellstr_value, Array<std::string>);
      SLOAV_OBV_METH_FROM_ARRAY_0(conj, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(cos, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(cosh, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(erf, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(erfc, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(exp, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(expm1, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(finite, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(fix, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(floor, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(gamma, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(imag, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(index_vector, idx_vector);
      SLOAV_OBV_METH_FROM_ARRAY_0(int16_array_value, int16NDArray);
      SLOAV_OBV_METH_FROM_ARRAY_0(int16_scalar_value, octave_int16);
      SLOAV_OBV_METH_FROM_ARRAY_0(int32_array_value, int32NDArray);
      SLOAV_OBV_METH_FROM_ARRAY_0(int32_scalar_value, octave_int32);
      SLOAV_OBV_METH_FROM_ARRAY_0(int64_array_value, int64NDArray);
      SLOAV_OBV_METH_FROM_ARRAY_0(int64_scalar_value, octave_int64);
      SLOAV_OBV_METH_FROM_ARRAY_0(int8_array_value, int8NDArray);
      SLOAV_OBV_METH_FROM_ARRAY_0(int8_scalar_value, octave_int8);
      SLOAV_OBV_METH_FROM_ARRAY_0(isinf, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(isna, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(isnan, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(lgamma, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(log, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(log10, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(log1p, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(log2, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(matrix_type, MatrixType);
      SLOAV_OBV_METH_FROM_ARRAY_0(nnz, octave_idx_type);
      SLOAV_OBV_METH_FROM_ARRAY_0(nzmax, octave_idx_type);
      SLOAV_OBV_METH_FROM_ARRAY_0(perm_matrix_value, PermMatrix);
      SLOAV_OBV_METH_FROM_ARRAY_0(range_value, Range);
      SLOAV_OBV_METH_FROM_ARRAY_0(real, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(round, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(roundb, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(signum, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(sin, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(sinh, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(sqrt, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(tan, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(tanh, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(uint16_array_value, uint16NDArray);
      SLOAV_OBV_METH_FROM_ARRAY_0(uint16_scalar_value, octave_uint16);
      SLOAV_OBV_METH_FROM_ARRAY_0(uint32_array_value, uint32NDArray);
      SLOAV_OBV_METH_FROM_ARRAY_0(uint32_scalar_value, octave_uint32);
      SLOAV_OBV_METH_FROM_ARRAY_0(uint64_array_value, uint64NDArray);
      SLOAV_OBV_METH_FROM_ARRAY_0(uint64_scalar_value, octave_uint64);
      SLOAV_OBV_METH_FROM_ARRAY_0(uint8_array_value, uint8NDArray);
      SLOAV_OBV_METH_FROM_ARRAY_0(uint8_scalar_value, octave_uint8);
      SLOAV_OBV_METH_FROM_ARRAY_0(xisalnum, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xisalpha, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xisascii, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xiscntrl, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xisdigit, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xisgraph, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xislower, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xisprint, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xispunct, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xisspace, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xisupper, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xisxdigit, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xtoascii, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xtolower, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_0(xtoupper, octave_value);
      SLOAV_OBV_METH_FROM_ARRAY_1(all, octave_value, int);
      SLOAV_OBV_METH_FROM_ARRAY_1(all_strings, string_vector, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(any, octave_value, int);
      SLOAV_OBV_METH_FROM_ARRAY_1(array_value, NDArray, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(bool_array_value, boolNDArray, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(bool_matrix_value, boolMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(bool_value, bool, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(char_array_value, charNDArray, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(char_matrix_value, charMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(complex_array_value, ComplexNDArray, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(complex_diag_matrix_value, ComplexDiagMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(complex_matrix_value, ComplexMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(complex_value, Complex, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(diag, octave_value, octave_idx_type);
      SLOAV_OBV_METH_FROM_ARRAY_1(diag_matrix_value, DiagMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(double_value, double, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(float_array_value, FloatNDArray, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(float_complex_array_value, FloatComplexNDArray, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(float_complex_diag_matrix_value, FloatComplexDiagMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(float_complex_matrix_value, FloatComplexMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(float_complex_value, FloatComplex, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(float_diag_matrix_value, FloatDiagMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(float_matrix_value, FloatMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(float_value, float, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(is_sorted, sortmode, sortmode);
      SLOAV_OBV_METH_FROM_ARRAY_1(is_sorted_rows, sortmode, sortmode);
      SLOAV_OBV_METH_FROM_ARRAY_1(matrix_value, Matrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(reshape, octave_value, const dim_vector&);
      SLOAV_OBV_METH_FROM_ARRAY_1(sort_rows_ids, Array<octave_idx_type>, sortmode);
      SLOAV_OBV_METH_FROM_ARRAY_1(sparse_bool_matrix_value, SparseBoolMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(sparse_complex_matrix_value, SparseComplexMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(sparse_matrix_value, SparseMatrix, bool);
      SLOAV_OBV_METH_FROM_ARRAY_1(string_value, std::string, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(int_value, int, bool, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(long_value, long int, bool, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(permute, octave_value, const Array<int>&, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(print_info, void, std::ostream&, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(print_raw, void, std::ostream&, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(resize, octave_value, const dim_vector&, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(short_value, short int, bool, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(sort, octave_value, octave_idx_type, sortmode);
      SLOAV_OBV_METH_FROM_ARRAY_2(uint_value, unsigned int, bool, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(ulong_value, unsigned long int, bool, bool);
      SLOAV_OBV_METH_FROM_ARRAY_2(ushort_value, unsigned short int, bool, bool);
      SLOAV_OBV_METH_FROM_ARRAY_3(convert_to_str_internal, octave_value, bool, bool, char);
      SLOAV_OBV_METH_FROM_ARRAY_3(sort, octave_value, Array<octave_idx_type>&, octave_idx_type, sortmode);
      SLOAV_OBV_METH_FROM_ARRAY_5(write, int, octave_stream&, int, oct_data_conv::data_type, int, oct_mach_info::float_format);
#undef SLOAV_OBV_METH_FROM_ARRAY_0
#undef SLOAV_OBV_METH_FROM_ARRAY_1
#undef SLOAV_OBV_METH_FROM_ARRAY_2
#undef SLOAV_OBV_METH_FROM_ARRAY_3
#undef SLOAV_OBV_METH_FROM_ARRAY_4
#undef SLOAV_OBV_METH_FROM_ARRAY_5

    }; // class swiglal_oct_array_view

  } // extern "C++"

} // fragment swiglal_oct_array_view

// Macro which generates fragments which define ACFTYPE-specific array view classes and conversion
// functions:
//  - IN/OUTFRAG are names of fragments required by the in/out conversion functions IN/OUTCALL.
//  - OVCLASS is the octave_base_value-derived class of the array view class.
//  - OVTYPE is the type of octave_value array value.
//  - OVVALUE() is the method of octave_value which returns an OVTYPE.
//  - ISOVTYPEEXPR returns true if the octave_value 'obj' is of type OVTYPE.
%define %swiglal_oct_array_frags(ACFTYPE, INFRAG, OUTFRAG, INCALL, OUTCALL, OVCLASS, OVTYPE, OVVALUE, ISOVTYPEEXPR)

// Register the ACFTYPE-specific array view class as an Octave type.
%fragment(%swiglal_oct_array_view_init_frag(ACFTYPE), "init") %{
#if SWIG_OCTAVE_PREREQ(4,4,0)
  {
    octave::type_info& typeinfo = octave::interpreter::the_interpreter()->get_type_info();
    string_vector types = typeinfo.installed_type_names();
    bool register_octave_array_view_class = true;
    bool register_octave_swig_packed = true;
    for (int i = 0; i < types.numel(); ++i) {
      if (types(i) == %swiglal_oct_array_view_class(ACFTYPE)::static_type_name()) {
        register_octave_array_view_class = false;
      }
    }
    if (register_octave_array_view_class) {
      %swiglal_oct_array_view_class(ACFTYPE)::register_type();
    }
  }
#else
  %swiglal_oct_array_view_class(ACFTYPE)::register_type();
#endif
%}

// ACFTYPE-specific array view class fragment.
%fragment(%swiglal_oct_array_view_frag(ACFTYPE), "header",
          fragment=%swiglal_oct_array_view_init_frag(ACFTYPE),
          fragment="swiglal_oct_array_view", fragment=INFRAG, fragment=OUTFRAG)
%{

  extern "C++" {

    // Helper class which supplies ACFTYPE-specific types and functions to the base template
    // swiglal_oct_array_view<>.
    class %swiglal_oct_array_view_helper_class(ACFTYPE) {

    public:

      // Octave array-related types.
      typedef OVCLASS OVClass;
      typedef OVTYPE OVType;

      // Extract an OVType from an octave_value.
      static OVType ovvalue(octave_value& ov) {
        return ov.OVVALUE();
      }

      // Convert the octave_value objelem to an array element stored at elemptr.
      static int incall(const octave_value& parent, octave_value& objelem, void *elemptr, int *pelemalloc, const size_t esize, const bool isptr, swig_type_info *const tinfo, const int tflags)
      {
        return INCALL;
      }

      // Convert the array element stored at elemptr to an octave_value.
      static octave_value outcall(const octave_value& parent, const bool copyobj, void *elemptr, const size_t esize, const bool isptr, swig_type_info *const tinfo, const int tflags) {
        return OUTCALL;
      }

    };

    // ACFTYPE-specific array view class fragment.
    class OCTINTERP_API %swiglal_oct_array_view_class(ACFTYPE) : public %swiglal_oct_array_view_tmpl(ACFTYPE) {

    private:

      // Octave type constructs.
#if !SWIG_OCTAVE_PREREQ(4,0,0)
      DECLARE_OCTAVE_ALLOCATOR;
#endif
      DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA;

    protected:

    public:

      virtual ~%swiglal_oct_array_view_class(ACFTYPE)()
      { }

      %swiglal_oct_array_view_class(ACFTYPE)()
        : %swiglal_oct_array_view_tmpl(ACFTYPE)()
      { }

      %swiglal_oct_array_view_class(ACFTYPE)(const octave_value& parent,
                                             void* ptr,
                                             const size_t esize,
                                             const size_t ndims,
                                             const size_t dims[],
                                             const size_t strides[],
                                             const bool isptr,
                                             swig_type_info* tinfo,
                                             const int tflags)
        : %swiglal_oct_array_view_tmpl(ACFTYPE)(parent, ptr, esize, ndims, dims, strides, isptr, tinfo, tflags)
      { }

    };

    // Octave type constructs.
#if !SWIG_OCTAVE_PREREQ(4,0,0)
    DEFINE_OCTAVE_ALLOCATOR(%swiglal_oct_array_view_class(ACFTYPE));
#endif
    DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(%swiglal_oct_array_view_class(ACFTYPE),
                                        %swiglal_oct_array_view_ovtype(ACFTYPE),
                                        OVCLASS::static_class_name());

  } // extern "C++"

%} // %swiglal_oct_array_view_frag()

// Input copy conversion fragment for arrays of type ACFTYPE.
%fragment(%swiglal_array_copyin_frag(ACFTYPE), "header",
          fragment=%swiglal_oct_array_view_frag(ACFTYPE))
%{
  SWIGINTERN int %swiglal_array_copyin_func(ACFTYPE)(const octave_value& parent,
                                                     octave_value obj,
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
    // Create a local array view, then use its sloav_array_in() member to copy the input Octave
    // array to the viewed C array.
    %swiglal_oct_array_view_class(ACFTYPE) arrview(parent, ptr, esize, ndims, dims, strides, isptr, tinfo, tflags);
    return arrview.sloav_array_in(obj, pelemalloc, 0);
  }
%}

// Output copy conversion fragment for arrays of type ACFTYPE.
%fragment(%swiglal_array_copyout_frag(ACFTYPE), "header",
          fragment=%swiglal_oct_array_view_frag(ACFTYPE))
%{
  SWIGINTERN octave_value %swiglal_array_copyout_func(ACFTYPE)(const octave_value& parent,
                                                               void* ptr,
                                                               const size_t esize,
                                                               const size_t ndims,
                                                               const size_t dims[],
                                                               const size_t strides[],
                                                               const bool isptr,
                                                               swig_type_info *tinfo,
                                                               const int tflags)
  {
    // Create a local array view, then use its sloav_array_out() member to copy the viewed C array
    // to the output Octave array.
    %swiglal_oct_array_view_class(ACFTYPE) arrview(parent, ptr, esize, ndims, dims, strides, isptr, tinfo, tflags);
    return arrview.sloav_array_out(true);
  }
%}

// Input view conversion fragment for arrays of type ACFTYPE.
%fragment(%swiglal_array_viewin_frag(ACFTYPE), "header",
          fragment=%swiglal_oct_array_view_frag(ACFTYPE))
%{
  SWIGINTERN int %swiglal_array_viewin_func(ACFTYPE)(const octave_value& parent,
                                                     octave_value obj,
                                                     void** ptr,
                                                     const size_t esize,
                                                     const size_t ndims,
                                                     size_t dims[],
                                                     const bool isptr,
                                                     swig_type_info *tinfo,
                                                     const int tflags)
  {

    // Check that C array pointer is valid.
    if (ptr == NULL) {
      return SWIG_MemoryError;
    }

    // Convert 'obj' to the desired OVTYPE, store in 'val'.
    OVTYPE val = obj.OVVALUE();

    // Check that 'val' has the correct number of dimensions.
    // - 1-D arrays are a special case, since Octave arrays are always at least
    //   2-dimensional, so need to check that one of those dimensions is singular.
    dim_vector valdims = val.dims();
    if (ndims == 1) {
      if (valdims.length() > 2 || valdims.num_ones() == 0) {
        return SWIG_ValueError;
      }
    }
    else if (val.ndims() != %reinterpret_cast(ndims, octave_idx_type)) {
      return SWIG_ValueError;
    }

    // Return dimensions of Octave array.
    if (ndims == 1) {
      dims[0] = val.numel();
    } else {
      for (size_t i = 0; i < ndims; ++i) {
        dims[i] = valdims(i);
      }
    }

    // Cannot view an array of pointers.
    if (isptr) {
      return SWIG_TypeError;
    }

    // Cannot view an array which is not itself an Octave view of a C array.
    if (obj.type_name().find("swiglal_oct_array_view_") == 0) {
      return SWIG_TypeError;
    }

    // Since Octave stores arrays in column-major order, we can only view 1-D arrays.
    if (ndims != 1) {
      return SWIG_ValueError;
    }

    // Check that 'obj' is of the correct type.
    if (!(ISOVTYPEEXPR)) {
      return SWIG_TypeError;
    }

    // Check that the elements of 'val' have the correct size.
    if (val.byte_size() != val.numel() * esize) {
      return SWIG_TypeError;
    }

    // Get pointer to Octave array data, a highly complicated and dodgy process!  Usually
    // mex_get_data() does the job, apart from complex arrays where that creates a copy ...
    // in which case try data() and try to detect copying ...
    if (obj.is_complex_type() && !obj.is_scalar_type()) {
      if (obj.is_double_type()) {
        Complex c;
        {
          const ComplexNDArray t = obj.complex_array_value();
          *ptr = %reinterpret_cast(t.data(), void*);
          c = *((Complex*) *ptr);
        }
        if (c != *((Complex*) *ptr)) {
          return SWIG_UnknownError;
        }
      } else if (obj.is_single_type()) {
        FloatComplex c;
        {
          const FloatComplexNDArray t = obj.float_complex_array_value();
          *ptr = %reinterpret_cast(t.data(), void*);
          c = *((FloatComplex*) *ptr);
        }
        if (c != *((FloatComplex*) *ptr)) {
          return SWIG_UnknownError;
        }
      } else {
        return SWIG_TypeError;
      }
    } else {
      *ptr = %reinterpret_cast(obj.mex_get_data(), void*);
    }
    if (!*ptr) {
      return SWIG_ValueError;
    }

    return SWIG_OK;

  }
%}

// Output view conversion fragment for arrays of type ACFTYPE.
%fragment(%swiglal_array_viewout_frag(ACFTYPE), "header",
          fragment=%swiglal_oct_array_view_frag(ACFTYPE))
%{
  SWIGINTERN octave_value %swiglal_array_viewout_func(ACFTYPE)(const octave_value& parent,
                                                               void* ptr,
                                                               const size_t esize,
                                                               const size_t ndims,
                                                               const size_t dims[],
                                                               const size_t strides[],
                                                               const bool isptr,
                                                               swig_type_info *tinfo,
                                                               const int tflags)
  {
    // Return an Octave array view of the C array.
    octave_base_value *objval = new %swiglal_oct_array_view_class(ACFTYPE)(parent, ptr, esize, ndims, dims, strides, isptr, tinfo, tflags);
    return octave_value(objval);
  }
%}

%enddef // %swiglal_oct_array_frags

// Array conversion fragments for generic arrays, e.g. SWIG-wrapped types.  Note that input views
// are not supported, and so ISOVTYPEEXPR is 'false'.
%swiglal_oct_array_frags(SWIGTYPE, "swiglal_as_SWIGTYPE", "swiglal_from_SWIGTYPE",
                         %arg(swiglal_as_SWIGTYPE(parent, objelem, elemptr, esize, isptr, tinfo, tflags)),
                         %arg(swiglal_from_SWIGTYPE(parent, copyobj, elemptr, esize, isptr, tinfo, tflags)),
                         octave_cell, Cell, cell_value, false);

// Array conversion fragments for arrays of LAL strings.  Note that input views are not supported,
// and so ISOVTYPEEXPR is 'false'.
%swiglal_oct_array_frags(LALchar, "SWIG_AsLALcharPtrAndSize", "SWIG_FromLALcharPtr",
                         %arg(SWIG_AsLALcharPtrAndSize(objelem, %reinterpret_cast(elemptr, char**), 0, pelemalloc)),
                         %arg(SWIG_FromLALcharPtr(*%reinterpret_cast(elemptr, char**))),
                         octave_cell, Cell, cell_value, false);

// Macro which generates array conversion function fragments to/from Octave arrays for real/fragment
// TYPEs which use SWIG_AsVal/From fragments.
%define %swiglal_oct_array_asvalfrom_frags(TYPE, OVCLASS, OVTYPE, OVVALUE, ISOVTYPEEXPR)
%swiglal_oct_array_frags(TYPE, SWIG_AsVal_frag(TYPE), SWIG_From_frag(TYPE),
                         %arg(SWIG_AsVal(TYPE)(objelem, %reinterpret_cast(elemptr, TYPE*))),
                         %arg(SWIG_From(TYPE)(*%reinterpret_cast(elemptr, TYPE*))),
                         OVCLASS, OVTYPE, OVVALUE, ISOVTYPEEXPR);
%enddef

// Array conversion fragments for integer arrays.
%swiglal_oct_array_asvalfrom_frags(int8_t, octave_int8_matrix, intNDArray<octave_int<int8_t> >, int8_array_value, obj.is_int8_type());
%swiglal_oct_array_asvalfrom_frags(uint8_t, octave_uint8_matrix, intNDArray<octave_int<uint8_t> >, uint8_array_value, obj.is_uint8_type());
%swiglal_oct_array_asvalfrom_frags(int16_t, octave_int16_matrix, intNDArray<octave_int<int16_t> >, int16_array_value, obj.is_int16_type());
%swiglal_oct_array_asvalfrom_frags(uint16_t, octave_uint16_matrix, intNDArray<octave_int<uint16_t> >, uint16_array_value, obj.is_uint16_type());
%swiglal_oct_array_asvalfrom_frags(int32_t, octave_int32_matrix, intNDArray<octave_int<int32_t> >, int32_array_value, obj.is_int32_type());
%swiglal_oct_array_asvalfrom_frags(uint32_t, octave_uint32_matrix, intNDArray<octave_int<uint32_t> >, uint32_array_value, obj.is_uint32_type());
%swiglal_oct_array_asvalfrom_frags(int64_t, octave_int64_matrix, intNDArray<octave_int<int64_t> >, int64_array_value, obj.is_int64_type());
%swiglal_oct_array_asvalfrom_frags(uint64_t, octave_uint64_matrix, intNDArray<octave_int<uint64_t> >, uint64_array_value, obj.is_uint64_type());

// Array conversion fragments for floating-precision real arrays.
%swiglal_oct_array_asvalfrom_frags(float, octave_float_matrix, FloatMatrix, float_matrix_value, obj.is_real_type() && obj.is_single_type());
%swiglal_oct_array_asvalfrom_frags(double, octave_matrix, Matrix, matrix_value, obj.is_real_type() && obj.is_double_type());

// Array conversion fragments for floating-precision complex arrays.
%swiglal_oct_array_asvalfrom_frags(gsl_complex_float, octave_float_complex_matrix, FloatComplexMatrix, float_complex_matrix_value, obj.is_complex_type() && obj.is_single_type());
%swiglal_oct_array_asvalfrom_frags(gsl_complex, octave_complex_matrix, ComplexMatrix, complex_matrix_value, obj.is_complex_type() && obj.is_double_type());
%swiglal_oct_array_asvalfrom_frags(COMPLEX8, octave_float_complex_matrix, FloatComplexMatrix, float_complex_matrix_value, obj.is_complex_type() && obj.is_single_type());
%swiglal_oct_array_asvalfrom_frags(COMPLEX16, octave_complex_matrix, ComplexMatrix, complex_matrix_value, obj.is_complex_type() && obj.is_double_type());

// Local Variables:
// mode: c
// End:
