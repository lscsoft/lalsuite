// Tests of SWIG interface code
// Author: Karl Wette, 2011, 2012

// Include LAL test code.
#include <lal/lalswig_test.i>

// Test object parent tracking between modules
typedef struct tagswig_lalxml_test_parent_map {
  swig_lal_test_struct s;
} swig_lalxml_test_parent_map_struct;
swig_lalxml_test_parent_map_struct swig_lalxml_test_parent_map;
