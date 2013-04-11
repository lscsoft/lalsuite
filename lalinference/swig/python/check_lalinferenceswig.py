# Check SWIG Python module wrapping lalinference
# Author: Karl Wette, 2011, 2012

# check module load
import lal
import lalinference
from lalinference import cvar as lalinferencecvar
from lal import cvar as lalcvar
lalcvar.lalDebugLevel = 1
print("passed module load")

# check object parent tracking
a = lalinference.lalinferenceswig_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalinferencecvar.lalinferenceswig_test_parent_map.s
    lalinferencecvar.lalinferenceswig_test_parent_map.s = lalcvar.lalswig_test_struct_const
del a, b, c
lal.CheckMemoryLeaks()
print("passed object parent tracking")

# passed all tests!
print("PASSED all tests")
