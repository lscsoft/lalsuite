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

// Test various combinations of 1D and 2D fixed arrays
// with structs, struct/enum type, and global variables.
typedef enum {
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
static const swig_lal_test_struct swig_lal_test_struct_const = {
  2, 3, 5.7, "abcde", {1,2,3}, {{4,5,6},{7,8,9}},
  {swig_lal_test_enum_a,swig_lal_test_enum_b,swig_lal_test_enum_c}
};
static swig_lal_test_struct swig_lal_test_struct_vector[3];
static swig_lal_test_struct swig_lal_test_struct_matrix[2][3];
static swig_lal_test_enum swig_lal_test_enum_vector[3];
static swig_lal_test_enum swig_lal_test_enum_matrix[2][3];
static INT4 swig_lal_test_empty_INT4_vector[0];
static INT4 swig_lal_test_INT4_vector[3];
static INT4 swig_lal_test_INT4_matrix[2][3];
static const INT4 swig_lal_test_INT4_const_vector[3] = {1, 2, 4};
static const INT4 swig_lal_test_INT4_const_matrix[2][3] = {{1, 2, 4}, {2, 4, 8}};
static REAL8 swig_lal_test_REAL8_vector[3];
static REAL8 swig_lal_test_REAL8_matrix[2][3];
static COMPLEX8 swig_lal_test_COMPLEX8_vector[3];
static COMPLEX8 swig_lal_test_COMPLEX8_matrix[2][3];

// Test fixed and dynamic arrays typemaps.
void swig_lal_test_copyin_array1(REAL8 INPUT[3], REAL8 scale, REAL8 OUTPUT[3]) {
  for (int i = 0; i < 3; ++i) {
    OUTPUT[i] = scale * INPUT[i];
  }
}
void swig_lal_test_copyin_array2(INT4 INPUT[3][2], INT4 scale, INT4 OUTPUT[3][2]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      OUTPUT[i][j] = scale * INPUT[i][j];
    }
  }
}

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
swig_lal_test_arrayofptrs* swig_lal_test_Create_arrayofptrs(UINT4 length) {
  swig_lal_test_arrayofptrs* ap = (swig_lal_test_arrayofptrs*)XLALMalloc(sizeof(swig_lal_test_arrayofptrs));
  XLAL_CHECK_NULL(ap != NULL, XLAL_ENOMEM);
  ap->length = length;
  ap->data = (swig_lal_test_arrayofdata**)XLALCalloc(ap->length, sizeof(swig_lal_test_arrayofdata*));
  XLAL_CHECK_NULL(ap->data != NULL, XLAL_ENOMEM);
  for (UINT4 i = 0; i < ap->length; ++i) {
    ap->data[i] = (swig_lal_test_arrayofdata*)XLALMalloc(sizeof(swig_lal_test_arrayofdata));
    XLAL_CHECK_NULL(ap->data[i] != NULL, XLAL_ENOMEM);
    ap->data[i]->length = 2*length;
    ap->data[i]->data = (INT4*)XLALCalloc(ap->data[i]->length, sizeof(INT4));
    XLAL_CHECK_NULL(ap->data[i]->data != NULL, XLAL_ENOMEM);
    for (UINT4 j = 0; j < ap->data[i]->length; ++j) {
      ap->data[i]->data[j] = 42*length*i + j;
    }
  }
  return ap;
}
void swig_lal_test_Destroy_arrayofptrs(swig_lal_test_arrayofptrs* ap) {
  if (ap) {
    if (ap->data) {
      for (UINT4 i = 0; i < ap->length; ++i) {
        if (ap->data[i]) {
          if (ap->data[i]->data) {
            XLALFree(ap->data[i]->data);
          }
          XLALFree(ap->data[i]);
        }
      }
      XLALFree(ap->data);
    }
    XLALFree(ap);
  }
}

// Test LIGOTimeGPS operations.
typedef struct tagswig_lal_test_gps {
  LIGOTimeGPS t;
} swig_lal_test_gps;
REAL8 swig_lal_test_noptrgps(const LIGOTimeGPS gps) {
  return XLALGPSGetREAL8(&gps);
}
