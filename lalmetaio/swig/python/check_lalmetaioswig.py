# Check SWIG Python module wrapping lalmetaio
# Author: Karl Wette, 2011, 2012

# check module load
import lal
import lalmetaio
from lalmetaio import cvar as lalmetaiocvar
from lal import cvar as lalcvar
lalcvar.lalDebugLevel = 1
print("passed module load")

# check object parent tracking
a = lalmetaio.lalmetaioswig_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalmetaiocvar.lalmetaioswig_test_parent_map.s
    lalmetaiocvar.lalmetaioswig_test_parent_map.s = lalcvar.lalswig_test_struct_const
del a, b, c
lal.CheckMemoryLeaks()
print("passed object parent tracking")

# passed all tests!
print("PASSED all tests")
