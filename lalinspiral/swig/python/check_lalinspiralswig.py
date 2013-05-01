# Check SWIG Python module wrapping lalinspiral
# Author: Karl Wette, 2011, 2012

# check module load
import lal
import lalinspiral
from lalinspiral import cvar as lalinspiralcvar
from lal import cvar as lalcvar
lalcvar.lalDebugLevel = lal.LALERROR | lal.LALMEMTRACE
print("passed module load")

# check object parent tracking
a = lalinspiral.lalinspiralswig_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalinspiralcvar.lalinspiralswig_test_parent_map.s
    lalinspiralcvar.lalinspiralswig_test_parent_map.s = lalcvar.lalswig_test_struct_const
del a, b, c
lal.CheckMemoryLeaks()
print("passed object parent tracking")

# passed all tests!
print("PASSED all tests")
