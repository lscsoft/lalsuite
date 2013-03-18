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
  INT4 i;
  REAL4 f;
  CHAR str[10];
  INT4 vec[3];
  INT4 mat[2][3];
  lalswig_test_enum evec[3];
} lalswig_test_struct;
const lalswig_test_struct lalswig_test_struct_const = {
  3, 5.7, "abcde", {1,2,3}, {{4,5,6},{7,8,9}},
  {lalswig_test_enum_a,lalswig_test_enum_b,lalswig_test_enum_c}
};
lalswig_test_struct lalswig_test_struct_vector[3];
lalswig_test_struct lalswig_test_struct_matrix[2][3];
lalswig_test_enum lalswig_test_enum_vector[3];
lalswig_test_enum lalswig_test_enum_matrix[2][3];
INT4 lalswig_test_INT4_vector[3];
INT4 lalswig_test_INT4_matrix[2][3];
const INT4 lalswig_test_INT4_const_vector[3] = {1, 2, 4};
const INT4 lalswig_test_INT4_const_matrix[2][3] = {{1, 2, 4}, {2, 4, 8}};
REAL8 lalswig_test_REAL8_vector[3];
REAL8 lalswig_test_REAL8_matrix[2][3];
COMPLEX8 lalswig_test_COMPLEX8_vector[3];
COMPLEX8 lalswig_test_COMPLEX8_matrix[2][3];
