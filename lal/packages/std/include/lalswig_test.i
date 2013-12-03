// Tests of SWIG interface code
// Author: Karl Wette, 2011, 2012

// Test various combinations of 1D and 2D fixed arrays with
// structs, struct/enum type, and global variables (LAL only).
typedef enum {
  lalswig_test_enum_a,
  lalswig_test_enum_b,
  lalswig_test_enum_c
} lalswig_test_enum;
typedef struct taglalswig_test_struct {
  int n;
  INT4 i;
  REAL4 f;
  CHAR str[10];
  INT4 vec[3];
  INT4 mat[2][3];
  lalswig_test_enum evec[3];
} lalswig_test_struct;
static const lalswig_test_struct lalswig_test_struct_const = {
  2, 3, 5.7, "abcde", {1,2,3}, {{4,5,6},{7,8,9}},
  {lalswig_test_enum_a,lalswig_test_enum_b,lalswig_test_enum_c}
};
static lalswig_test_struct lalswig_test_struct_vector[3];
static lalswig_test_struct lalswig_test_struct_matrix[2][3];
static lalswig_test_enum lalswig_test_enum_vector[3];
static lalswig_test_enum lalswig_test_enum_matrix[2][3];
static INT4 lalswig_test_empty_INT4_vector[0];
static INT4 lalswig_test_INT4_vector[3];
static INT4 lalswig_test_INT4_matrix[2][3];
static const INT4 lalswig_test_INT4_const_vector[3] = {1, 2, 4};
static const INT4 lalswig_test_INT4_const_matrix[2][3] = {{1, 2, 4}, {2, 4, 8}};
static REAL8 lalswig_test_REAL8_vector[3];
static REAL8 lalswig_test_REAL8_matrix[2][3];
static COMPLEX8 lalswig_test_COMPLEX8_vector[3];
static COMPLEX8 lalswig_test_COMPLEX8_matrix[2][3];

// Test dynamic array of pointer access
typedef struct taglalswig_test_arrayofdata {
#ifdef SWIG
  SWIGLAL(ARRAY_1D(lalswig_test_arrayofdata, INT4, data, UINT4, length));
#endif // SWIG
  UINT4 length;
  INT4 *data;
} lalswig_test_arrayofdata;
typedef struct taglalswig_test_arrayofptrs {
#ifdef SWIG
  SWIGLAL(ARRAY_1D(lalswig_test_arrayofptrs, lalswig_test_arrayofdata*, data, UINT4, length));
#endif // SWIG
  UINT4 length;
  lalswig_test_arrayofdata **data;
} lalswig_test_arrayofptrs;
static lalswig_test_arrayofptrs* lalswig_test_Create_arrayofptrs(UINT4);
static void lalswig_test_Destroy_arrayofptrs(lalswig_test_arrayofptrs*);
#ifdef SWIG
%header %{
  static lalswig_test_arrayofptrs* lalswig_test_Create_arrayofptrs(UINT4 length) {
    lalswig_test_arrayofptrs* ap = (lalswig_test_arrayofptrs*)XLALMalloc(sizeof(lalswig_test_arrayofptrs));
    XLAL_CHECK_NULL(ap != NULL, XLAL_ENOMEM);
    ap->length = length;
    ap->data = (lalswig_test_arrayofdata**)XLALCalloc(ap->length, sizeof(lalswig_test_arrayofdata*));
    XLAL_CHECK_NULL(ap->data != NULL, XLAL_ENOMEM);
    for (UINT4 i = 0; i < ap->length; ++i) {
      ap->data[i] = (lalswig_test_arrayofdata*)XLALMalloc(sizeof(lalswig_test_arrayofdata));
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
  static void lalswig_test_Destroy_arrayofptrs(lalswig_test_arrayofptrs* ap) {
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
%}
#endif // SWIG

// Test LIGOTimeGPS operations.
typedef struct taglalswig_test_gps {
  LIGOTimeGPS t;
} lalswig_test_gps;
static REAL8 lalswig_test_noptrgps(const LIGOTimeGPS gps);
#ifdef SWIG
%header %{
  static REAL8 lalswig_test_noptrgps(const LIGOTimeGPS gps) {
    return XLALGPSGetREAL8(&gps);
  }
%}
#endif // SWIG
