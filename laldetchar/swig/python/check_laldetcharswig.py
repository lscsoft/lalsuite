# Check SWIG Python module wrapping laldetchar
# Author: Karl Wette, 2011, 2012

# check module load
print("checking module load ...")
import lal
import laldetchar
from laldetchar import cvar as laldetcharcvar
from lal import cvar as lalcvar
print("PASSED module load")

# check object parent tracking
print("checking object parent tracking ...")
a = laldetchar.laldetcharswig_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = laldetcharcvar.laldetcharswig_test_parent_map.s
    laldetcharcvar.laldetcharswig_test_parent_map.s = lalcvar.lalswig_test_struct_const
del a, b, c
lal.CheckMemoryLeaks()
print("PASSED object parent tracking")

# passed all tests!
print("PASSED all tests")
