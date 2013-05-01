# Check SWIG Python module wrapping lalburst
# Author: Karl Wette, 2011, 2012

# check module load
import lal
import lalburst
from lalburst import cvar as lalburstcvar
from lal import cvar as lalcvar
lalcvar.lalDebugLevel = lal.LALERROR | lal.LALMEMTRACE
print("passed module load")

# check object parent tracking
a = lalburst.lalburstswig_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalburstcvar.lalburstswig_test_parent_map.s
    lalburstcvar.lalburstswig_test_parent_map.s = lalcvar.lalswig_test_struct_const
del a, b, c
lal.CheckMemoryLeaks()
print("passed object parent tracking")

# passed all tests!
print("PASSED all tests")
