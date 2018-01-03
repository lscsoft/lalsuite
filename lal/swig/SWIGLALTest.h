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

// Code for SWIG tests of the LAL bindings.
// Author: Karl Wette

#ifndef _SWIGLALTEST_H
#define _SWIGLALTEST_H

#ifdef SWIGLAL_HAVE_LIBGSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#endif // SWIGLAL_HAVE_LIBGSL

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif

// Test various combinations of 1D and 2D fixed arrays
// with structs, struct/enum type, and global variables.
typedef enum tagswig_lal_test_enum {
  swig_lal_test_enum_a,
  swig_lal_test_enum_b,
  swig_lal_test_enum_c
} swig_lal_test_enum;
typedef struct tagswig_lal_test_struct {
  int n;
  INT4 i;
  REAL4 f;
  CHAR str[10];
  INT4 vec[3];
  INT4 mat[2][3];
  swig_lal_test_enum evec[3];
} swig_lal_test_struct;
extern const swig_lal_test_struct swig_lal_test_struct_const;
extern swig_lal_test_struct swig_lal_test_struct_vector[3];
extern swig_lal_test_struct swig_lal_test_struct_matrix[2][3];
extern swig_lal_test_enum swig_lal_test_enum_vector[3];
extern swig_lal_test_enum swig_lal_test_enum_matrix[2][3];
extern INT4 swig_lal_test_empty_INT4_vector[0];
extern INT4 swig_lal_test_INT4_vector[3];
extern INT4 swig_lal_test_INT4_matrix[2][3];
extern const INT4 swig_lal_test_INT4_const_vector[3];
extern const INT4 swig_lal_test_INT4_const_matrix[2][3];
extern REAL8 swig_lal_test_REAL8_vector[3];
extern REAL8 swig_lal_test_REAL8_matrix[2][3];
extern COMPLEX8 swig_lal_test_COMPLEX8_vector[3];
extern COMPLEX8 swig_lal_test_COMPLEX8_matrix[2][3];

// Test fixed and dynamic arrays typemaps.
void swig_lal_test_copyin_array1(REAL8 INPUT[3], REAL8 scale, REAL8 OUTPUT[3]);
void swig_lal_test_copyin_array2(INT4 INPUT[3][2], INT4 scale, INT4 OUTPUT[3][2]);

// Test input views of array structs.
BOOLEAN swig_lal_test_viewin_REAL4Vector(REAL4Vector* copyout, const REAL4Vector* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(REAL4Vector, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_REAL4Vector(REAL4Vector* viewout, REAL4Vector* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(REAL4Vector, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(REAL4Vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_REAL4Vector(REAL4Vector* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(REAL4Vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_REAL8Vector(REAL8Vector* copyout, const REAL8Vector* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(REAL8Vector, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_REAL8Vector(REAL8Vector* viewout, REAL8Vector* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(REAL8Vector, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(REAL8Vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_REAL8Vector(REAL8Vector* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(REAL8Vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_COMPLEX8Vector(COMPLEX8Vector* copyout, const COMPLEX8Vector* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(COMPLEX8Vector, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_COMPLEX8Vector(COMPLEX8Vector* viewout, COMPLEX8Vector* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(COMPLEX8Vector, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(COMPLEX8Vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_COMPLEX8Vector(COMPLEX8Vector* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(COMPLEX8Vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_COMPLEX16Vector(COMPLEX16Vector* copyout, const COMPLEX16Vector* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(COMPLEX16Vector, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_COMPLEX16Vector(COMPLEX16Vector* viewout, COMPLEX16Vector* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(COMPLEX16Vector, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(COMPLEX16Vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_COMPLEX16Vector(COMPLEX16Vector* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(COMPLEX16Vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_REAL4VectorSequence(REAL4VectorSequence* copyout, const REAL4VectorSequence* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(REAL4VectorSequence, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_REAL4VectorSequence(REAL4VectorSequence* viewout, REAL4VectorSequence* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(REAL4VectorSequence, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(REAL4VectorSequence, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_REAL4VectorSequence(REAL4VectorSequence* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(REAL4VectorSequence, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_REAL8VectorSequence(REAL8VectorSequence* copyout, const REAL8VectorSequence* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(REAL8VectorSequence, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_REAL8VectorSequence(REAL8VectorSequence* viewout, REAL8VectorSequence* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(REAL8VectorSequence, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(REAL8VectorSequence, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_REAL8VectorSequence(REAL8VectorSequence* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(REAL8VectorSequence, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_COMPLEX8VectorSequence(COMPLEX8VectorSequence* copyout, const COMPLEX8VectorSequence* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(COMPLEX8VectorSequence, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_COMPLEX8VectorSequence(COMPLEX8VectorSequence* viewout, COMPLEX8VectorSequence* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(COMPLEX8VectorSequence, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(COMPLEX8VectorSequence, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_COMPLEX8VectorSequence(COMPLEX8VectorSequence* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(COMPLEX8VectorSequence, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_COMPLEX16VectorSequence(COMPLEX16VectorSequence* copyout, const COMPLEX16VectorSequence* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(COMPLEX16VectorSequence, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_COMPLEX16VectorSequence(COMPLEX16VectorSequence* viewout, COMPLEX16VectorSequence* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(COMPLEX16VectorSequence, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(COMPLEX16VectorSequence, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_COMPLEX16VectorSequence(COMPLEX16VectorSequence* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(COMPLEX16VectorSequence, copyinout));
#endif // SWIG
#ifdef SWIGLAL_HAVE_LIBGSL
BOOLEAN swig_lal_test_viewin_gsl_vector_float(gsl_vector_float* copyout, const gsl_vector_float* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(gsl_vector_float, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_gsl_vector_float(gsl_vector_float* viewout, gsl_vector_float* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(gsl_vector_float, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(gsl_vector_float, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_gsl_vector_float(gsl_vector_float* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(gsl_vector_float, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_gsl_vector(gsl_vector* copyout, const gsl_vector* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(gsl_vector, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_gsl_vector(gsl_vector* viewout, gsl_vector* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(gsl_vector, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(gsl_vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_gsl_vector(gsl_vector* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(gsl_vector, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_gsl_vector_complex_float(gsl_vector_complex_float* copyout, const gsl_vector_complex_float* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(gsl_vector_complex_float, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_gsl_vector_complex_float(gsl_vector_complex_float* viewout, gsl_vector_complex_float* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(gsl_vector_complex_float, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(gsl_vector_complex_float, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_gsl_vector_complex_float(gsl_vector_complex_float* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(gsl_vector_complex_float, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_gsl_vector_complex(gsl_vector_complex* copyout, const gsl_vector_complex* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(gsl_vector_complex, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_gsl_vector_complex(gsl_vector_complex* viewout, gsl_vector_complex* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(gsl_vector_complex, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(gsl_vector_complex, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_gsl_vector_complex(gsl_vector_complex* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(gsl_vector_complex, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_gsl_matrix_float(gsl_matrix_float* copyout, const gsl_matrix_float* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(gsl_matrix_float, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_gsl_matrix_float(gsl_matrix_float* viewout, gsl_matrix_float* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(gsl_matrix_float, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(gsl_matrix_float, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_gsl_matrix_float(gsl_matrix_float* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(gsl_matrix_float, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_gsl_matrix(gsl_matrix* copyout, const gsl_matrix* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(gsl_matrix, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_gsl_matrix(gsl_matrix* viewout, gsl_matrix* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(gsl_matrix, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(gsl_matrix, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_gsl_matrix(gsl_matrix* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(gsl_matrix, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_gsl_matrix_complex_float(gsl_matrix_complex_float* copyout, const gsl_matrix_complex_float* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(gsl_matrix_complex_float, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_gsl_matrix_complex_float(gsl_matrix_complex_float* viewout, gsl_matrix_complex_float* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(gsl_matrix_complex_float, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(gsl_matrix_complex_float, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_gsl_matrix_complex_float(gsl_matrix_complex_float* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(gsl_matrix_complex_float, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_viewin_gsl_matrix_complex(gsl_matrix_complex* copyout, const gsl_matrix_complex* viewin);
#ifdef SWIG
SWIGLAL(VIEWIN_ARRAYS(gsl_matrix_complex, viewin, viewout));
#endif // SWIG
BOOLEAN swig_lal_test_viewinout_gsl_matrix_complex(gsl_matrix_complex* viewout, gsl_matrix_complex* viewin);
#ifdef SWIG
SWIGLAL_CLEAR(VIEWIN_ARRAYS(gsl_matrix_complex, viewin, viewout));
#endif // SWIG
#ifdef SWIG
SWIGLAL(COPYINOUT_ARRAYS(gsl_matrix_complex, copyinout));
#endif // SWIG
BOOLEAN swig_lal_test_copyinout_gsl_matrix_complex(gsl_matrix_complex* copyinout);
#ifdef SWIG
SWIGLAL_CLEAR(COPYINOUT_ARRAYS(gsl_matrix_complex, copyinout));
#endif // SWIG
#endif // SWIGLAL_HAVE_LIBGSL

// Test dynamic array of pointer access.
typedef struct tagswig_lal_test_arrayofdata {
#ifdef SWIG
  SWIGLAL(ARRAY_1D(swig_lal_test_arrayofdata, INT4, data, UINT4, length));
#endif // SWIG
  UINT4 length;
  INT4 *data;
} swig_lal_test_arrayofdata;
typedef struct tagswig_lal_test_arrayofptrs {
#ifdef SWIG
  SWIGLAL(ARRAY_1D(swig_lal_test_arrayofptrs, swig_lal_test_arrayofdata*, data, UINT4, length));
#endif // SWIG
  UINT4 length;
  swig_lal_test_arrayofdata **data;
} swig_lal_test_arrayofptrs;
swig_lal_test_arrayofptrs* swig_lal_test_Create_arrayofptrs(UINT4 length);
void swig_lal_test_Destroy_arrayofptrs(swig_lal_test_arrayofptrs* ap);

// Test LIGOTimeGPS operations.
typedef struct tagswig_lal_test_gps {
  LIGOTimeGPS t;
} swig_lal_test_gps;
REAL8 swig_lal_test_noptrgps(const LIGOTimeGPS gps);

#ifdef __cplusplus
}
#endif

#endif // _SWIGLALTEST_H
