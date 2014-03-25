// Tests of SWIG interface code
// Author: Karl Wette, 2011, 2012

// Include LAL test code.
#include <lal/lalswig_test.i>

// Test object parent tracking between modules
typedef struct tagswig_lalframe_test_parent_map {
  swig_lal_test_struct s;
} swig_lalframe_test_parent_map_struct;
swig_lalframe_test_parent_map_struct swig_lalframe_test_parent_map;
