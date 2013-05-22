# Check SWIG Python module wrapping lalstochastic
# Author: Karl Wette, 2011, 2012

# check module load
import lal
import lalstochastic
from lalstochastic import cvar as lalstochasticcvar
from lal import cvar as lalcvar
print("passed module load")

# check object parent tracking
a = lalstochastic.lalstochasticswig_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalstochasticcvar.lalstochasticswig_test_parent_map.s
    lalstochasticcvar.lalstochasticswig_test_parent_map.s = lalcvar.lalswig_test_struct_const
del a, b, c
lal.CheckMemoryLeaks()
print("passed object parent tracking")

# passed all tests!
print("PASSED all tests")
