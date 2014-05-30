# Check SWIG Python bindings for LALInference
# Author: Karl Wette, 2011--2014

# check module load
print("checking module load ...")
import lal
import lalinference
from lalinference import globalvar as lalinferenceglobalvar
from lal import globalvar as lalglobalvar
print("PASSED module load")

# check object parent tracking
print("checking object parent tracking ...")
a = lalinference.swig_lalinference_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalinferenceglobalvar.swig_lalinference_test_parent_map.s
    lalinferenceglobalvar.swig_lalinference_test_parent_map.s = lalglobalvar.swig_lal_test_struct_const
del c
del b
del a
lal.CheckMemoryLeaks()
print("PASSED object parent tracking")

# passed all tests!
print("PASSED all tests")
