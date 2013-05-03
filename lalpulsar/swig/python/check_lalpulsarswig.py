# Check SWIG Python module wrapping lalpulsar
# Author: Karl Wette, 2011, 2012

# check module load
import lal
import lalpulsar
from lalpulsar import cvar as lalpulsarcvar
from lal import cvar as lalcvar
lalcvar.lalDebugLevel = lal.LALERROR | lal.LALMEMDBG
print("passed module load")

# check object parent tracking
a = lalpulsar.lalpulsarswig_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalpulsarcvar.lalpulsarswig_test_parent_map.s
    lalpulsarcvar.lalpulsarswig_test_parent_map.s = lalcvar.lalswig_test_struct_const
del a, b, c
lal.CheckMemoryLeaks()
print("passed object parent tracking")

# passed all tests!
print("PASSED all tests")
