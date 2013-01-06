// Tests of SWIG interface code
// Author: Karl Wette, 2011, 2012

// Include LAL test code.
#include <lal/lalswig_test.i>

// Test object parent tracking between modules
typedef struct taglalinspiralswig_test_parent_map {
  lalswig_test_struct s;
} lalinspiralswig_test_parent_map_struct;
lalinspiralswig_test_parent_map_struct lalinspiralswig_test_parent_map;
