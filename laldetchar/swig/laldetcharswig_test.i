// Tests of SWIG interface code
// Author: Karl Wette, 2011, 2012

// Include LAL test code.
#include <lal/lalswig_test.i>

// Test object parent tracking between modules
typedef struct taglaldetcharswig_test_parent_map {
  lalswig_test_struct s;
} laldetcharswig_test_parent_map_struct;
laldetcharswig_test_parent_map_struct laldetcharswig_test_parent_map;
