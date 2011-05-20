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

// SWIG interface code specific to Octave
// Author: Karl Wette, 2011

// This function is declared in the Octave runtime but never defined
// (as of SWIG 2.0.4). This empty definition is to prevents warnings.
%runtime %{
  void Swig::swig_register_director(octave_swig_type*, void*, Swig::Director*) {}
%}

// octcomplex.swg was broken for single-precision complex numbers; fixed
// in SWIG version 2.0.4. For older SWIG we provide the fixed version here.
#if SWIG_VERSION < 0x020004
/////////////// license Lib/octave/octcomplex.swg ///////////////
// SWIG is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version. See the LICENSE-GPL file for
// the full terms of the GNU General Public license version 3.
//
// Portions of SWIG are also licensed under the terms of the licenses
// in the file LICENSE-UNIVERSITIES. You must observe the terms of
// these licenses, as well as the terms of the GNU General Public License,
// when you distribute SWIG.
//
// The SWIG library and examples, under the Lib and Examples top level
// directories, are distributed under the following terms:
//
//   You may copy, modify, distribute, and make derivative works based on
//   this software, in source code or object code form, without
//   restriction. If you distribute the software to others, you may do
//   so according to the terms of your choice. This software is offered as
//   is, without warranty of any kind.
//
// See the COPYRIGHT file for a list of contributors to SWIG and their
// copyright notices.
/////////////// begin Lib/octave/octcomplex.swg ///////////////
/*
  Defines the As/From conversors for double/float complex, you need to
  provide complex Type, the Name you want to use in the conversors,
  the complex Constructor method, and the Real and Imag complex
  accesor methods.
  See the std_complex.i and ccomplex.i for concrete examples.
*/
/* the common from conversor */
%define %swig_fromcplx_conv(Type, OctConstructor, Real, Imag)
     %fragment(SWIG_From_frag(Type),"header")
{
  SWIGINTERNINLINE octave_value
    SWIG_From(Type)(const Type& c)
    {
      return octave_value(OctConstructor(Real(c), Imag(c)));
    }
}
%enddef
// the double case
%define %swig_cplxdbl_conv(Type, Constructor, Real, Imag)
     %fragment(SWIG_AsVal_frag(Type),"header",
               fragment=SWIG_AsVal_frag(double))
{
  SWIGINTERN int
    SWIG_AsVal(Type) (const octave_value& ov, Type* val)
    {
      if (ov.is_complex_scalar()) {
        if (val) {
          Complex c(ov.complex_value());
          *val=Constructor(c.real(),c.imag());
        }
        return SWIG_OK;
      } else {
        double d;
        int res = SWIG_AddCast(SWIG_AsVal(double)(ov, &d));
        if (SWIG_IsOK(res)) {
          if (val)
            *val = Constructor(d, 0.0);
          return res;
        }
      }
      return SWIG_TypeError;
    }
}
%swig_fromcplx_conv(Type, Complex, Real, Imag);
%enddef
// the float case
%define %swig_cplxflt_conv(Type, Constructor, Real, Imag)
     %fragment(SWIG_AsVal_frag(Type),"header",
               fragment=SWIG_AsVal_frag(float)) {
  SWIGINTERN int
    SWIG_AsVal(Type) (const octave_value& ov, Type* val)
    {
      if (ov.is_complex_scalar()) {
        if (val) {
          Complex c(ov.complex_value());
          double re = c.real();
          double im = c.imag();
          if ((-FLT_MAX <= re && re <= FLT_MAX) && (-FLT_MAX <= im && im <= FLT_MAX)) {
            if (val)
              *val = Constructor(%numeric_cast(re, float),
                                 %numeric_cast(im, float));
            return SWIG_OK;
          } else
            return SWIG_OverflowError;
        }
      } else {
        float d;
        int res = SWIG_AddCast(SWIG_AsVal(float)(ov, &d));
        if (SWIG_IsOK(res)) {
          if (val)
            *val = Constructor(d, 0.0);
          return res;
        }
      }
      return SWIG_TypeError;
    }
}
%swig_fromcplx_conv(Type, FloatComplex, Real, Imag);
%enddef
#define %swig_cplxflt_convn(Type, Constructor, Real, Imag) \
%swig_cplxflt_conv(Type, Constructor, Real, Imag)
#define %swig_cplxdbl_convn(Type, Constructor, Real, Imag) \
%swig_cplxdbl_conv(Type, Constructor, Real, Imag)
/////////////// end Lib/octave/octcomplex.swg ///////////////
#else
%include <octcomplex.swg>
#endif

// Provide SWIG UTL conversion functions SWIG_From and SWIG_AsVal
// between Octave and LAL/GSL complex numbers.
%swig_cplxflt_convn(gsl_complex_float, gsl_complex_float_rect, GSL_REAL, GSL_IMAG);
%swig_cplxdbl_convn(gsl_complex, gsl_complex_rect, GSL_REAL, GSL_IMAG);
%swig_cplxflt_convn(COMPLEX8, XLALCOMPLEX8Rect, LAL_REAL, LAL_IMAG);
%swig_cplxdbl_convn(COMPLEX16, XLALCOMPLEX16Rect, LAL_REAL, LAL_IMAG);

// Return true, since an octave_value is always a valid object.
#define swiglal_object_valid(OBJ)   true

// Do nothing, since an octave_value is automatically destroyed.
#define swiglal_object_free(OBJ)   /*nothing*/

// Functions for manipulating vectors in Octave.
%header %{

  // Return whether an octave_value is a non-zero-length vector,
  // i.e. whether it has 2 dimensions (all octave_value array
  // objects have at least 2 dimensions), one of which  must be
  // of size 1, and the other must be of non-zero size.
  SWIGINTERN bool swiglal_is_vector(const octave_value& v) {
    return v.ndims() == 2 &&
      ((v.rows() == 1 && v.columns() >= 1) ||
       (v.rows() >= 1 && v.columns() == 1));
  }

  // Return the length of an octave_value vector.
  SWIGINTERN size_t swiglal_vector_length(const octave_value& v) {
    return v.length();
  }

  // Get the (i)th element of an octave_value vector.
  SWIGINTERN octave_value swiglal_vector_get(octave_value v, const size_t i) {
    std::list<octave_value_list> idx(1);
    idx.front()(0) = i + 1;
    return v.subsref(v.is_cell() ? "{" : "(", idx);
  }

  // Set the (i)th element of an octave_value vector to vi.
  SWIGINTERN bool swiglal_vector_set(octave_value& v, const size_t i, octave_value vi) {
    std::list<octave_value_list> idx(1);
    idx.front()(0) = i + 1;
    v = v.subsasgn(v.is_cell() ? "{" : "(", idx, vi);
    return true;
  }

  // Create a new octave_value vector for a general type.
  // We return a one-dimensional cell array, so that it can
  // contain arbitrary data such as swig_type-wrapped pointers.
  template<class TYPE > SWIGINTERN octave_value swiglal_new_vector(const size_t n) {
    return octave_value(Cell(dim_vector(1, n)));
  }

%}

// Functions for manipulating matrices in Octave.
%header %{

  // Return whether an octave_value is a non-zero-size matrix,
  // i.e. whether it has 2 dimensions, each of non-zero size.
  SWIGINTERN bool swiglal_is_matrix(const octave_value& v) {
    return v.ndims() == 2 &&
      v.rows() >= 1 && v.columns() >= 1;
  }

  // Return the number of rows of an octave_value matrix.
  SWIGINTERN size_t swiglal_matrix_rows(const octave_value& v) {
    return v.rows();
  }

  // Return the number of columns of an octave_value matrix.
  SWIGINTERN size_t swiglal_matrix_cols(const octave_value& v) {
    return v.columns();
  }

  // Get the (i,j)th element of an octave_value matrix.
  SWIGINTERN octave_value swiglal_matrix_get(octave_value v, const size_t i, const size_t j) {
    std::list<octave_value_list> idx(1);
    idx.front()(0) = i + 1;
    idx.front()(1) = j + 1;
    return v.subsref(v.is_cell() ? "{" : "(", idx);
  }

  // Set the (i,j)th element of an octave_value matrix to vij.
  SWIGINTERN bool swiglal_matrix_set(octave_value& v, const size_t i, const size_t j, octave_value vij) {
    std::list<octave_value_list> idx(1);
    idx.front()(0) = i + 1;
    idx.front()(1) = j + 1;
    v = v.subsasgn(v.is_cell() ? "{" : "(", idx, vij);
    return true;
  }

  // Create a new octave_value matrix for a general type.
  // We return a two-dimensional cell array, so that it can
  // contain arbitrary data such as swig_type-wrapped pointers.
  template<class TYPE > SWIGINTERN octave_value swiglal_new_matrix(const size_t ni, const size_t nj) {
    return octave_value(Cell(dim_vector(ni, nj)));
  }

%}

// Create new octave_value vectors and matrices of a specific TYPE.
%define swiglal_new_oct_vecmat(TYPE, CTOR)
  %header %{
    template<> octave_value swiglal_new_vector<TYPE >(const size_t n) {
      return octave_value(CTOR(dim_vector(1, n)));
    }
    template<> octave_value swiglal_new_matrix<TYPE >(const size_t ni, const size_t nj) {
      return octave_value(CTOR(dim_vector(ni, nj)));
    }
  %}
%enddef

// Create new octave_value vectors and matrices for integer types
%define swiglal_new_oct_int_vecmat(TYPE)
  swiglal_new_oct_vecmat(TYPE, intNDArray<octave_int<TYPE > >)
%enddef
swiglal_new_oct_int_vecmat(int8_t);
swiglal_new_oct_int_vecmat(uint8_t);
swiglal_new_oct_int_vecmat(int16_t);
swiglal_new_oct_int_vecmat(uint16_t);
swiglal_new_oct_int_vecmat(int32_t);
swiglal_new_oct_int_vecmat(uint32_t);
swiglal_new_oct_int_vecmat(int64_t);
swiglal_new_oct_int_vecmat(uint64_t);

// Create new octave_value vectors and matrices for real and complex types
swiglal_new_oct_vecmat(float, FloatMatrix);
swiglal_new_oct_vecmat(double, Matrix);
swiglal_new_oct_vecmat(gsl_complex_float, FloatComplexMatrix);
swiglal_new_oct_vecmat(gsl_complex, ComplexMatrix);
swiglal_new_oct_vecmat(COMPLEX8, FloatComplexMatrix);
swiglal_new_oct_vecmat(COMPLEX16, ComplexMatrix);
