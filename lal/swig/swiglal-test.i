// tests of SWIG interface code
// Author: Karl Wette, 2011

// Test variables/vectors/matrices for arrays of struct types.
// Use to test that *_getel/*_setel accessors work correctly.
%inline %{
  struct swiglal_test_struct_type {
    int a;
    float b;
    char c[10];
  };
  const swiglal_test_struct_type swiglal_test_struct_const = {3, 5.7, "abcde"};
  swiglal_test_struct_type swiglal_test_struct_vector[1];
  swiglal_test_struct_type swiglal_test_struct_matrix[1][1];
%}
SWIGLAL_GLOBAL_FIXED_1DARRAY_ELEM(swiglal_test_struct_type, swiglal_test_struct_vector);
SWIGLAL_GLOBAL_FIXED_2DARRAY_ELEM(swiglal_test_struct_type, swiglal_test_struct_matrix);

// Not all of the static vector/matrix member/global variable/enum
// combinations appear to exist in LAL code, so the following definitions
// exist to ensure all possible combinations are tested. They also test
// the SWIGLAL*_ELEM macros, which provides *_getel and *_setel methods.
%inline %{
  enum swiglal_test_static_enum {
    swiglal_test_static_enum_a
  };
  struct swiglal_test_static_struct {
    int vector[3];
    int matrix[2][3];
    swiglal_test_static_enum enum_vector[3];
    swiglal_test_static_enum enum_matrix[2][3];
  };
  int swiglal_test_static_vector[3];
  int swiglal_test_static_matrix[2][3];
  swiglal_test_static_enum swiglal_test_static_enum_vector[3];
  swiglal_test_static_enum swiglal_test_static_enum_matrix[2][3];
  const int swiglal_test_static_const_vector[3] = {1, 2, 4};
  const int swiglal_test_static_const_matrix[2][3] = {{1, 2, 4}, {2, 4, 8}};
%}
SWIGLAL_FIXED_1DARRAY_ELEM(int, vector,
                           swiglal_test_static_struct);
SWIGLAL_FIXED_2DARRAY_ELEM(int, matrix,
                           swiglal_test_static_struct);
SWIGLAL_FIXED_1DARRAY_ELEM(swiglal_test_static_enum, enum_vector,
                           swiglal_test_static_struct);
SWIGLAL_FIXED_2DARRAY_ELEM(swiglal_test_static_enum, enum_matrix,
                           swiglal_test_static_struct);
SWIGLAL_GLOBAL_FIXED_1DARRAY_ELEM(int, swiglal_test_static_vector);
SWIGLAL_GLOBAL_FIXED_2DARRAY_ELEM(int, swiglal_test_static_matrix);
SWIGLAL_GLOBAL_FIXED_1DARRAY_ELEM(swiglal_test_static_enum, swiglal_test_static_enum_vector);
SWIGLAL_GLOBAL_FIXED_2DARRAY_ELEM(swiglal_test_static_enum, swiglal_test_static_enum_matrix);
SWIGLAL_GLOBAL_CONST_FIXED_1DARRAY_ELEM(int, swiglal_test_static_const_vector);
SWIGLAL_GLOBAL_CONST_FIXED_2DARRAY_ELEM(int, swiglal_test_static_const_matrix);
