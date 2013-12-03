// Tests of SWIG interface code
// Author: Karl Wette, 2011, 2012

// Include LAL test code.
#include <lal/lalswig_test.i>

// Test object parent tracking between modules
typedef struct taglalmetaioswig_test_parent_map {
  lalswig_test_struct s;
} lalmetaioswig_test_parent_map_struct;
lalmetaioswig_test_parent_map_struct lalmetaioswig_test_parent_map;
