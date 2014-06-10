# Check SWIG Python bindings for LALDetChar
# Author: Karl Wette, 2011--2014

# check module load
print("checking module load ...")
import lal
import laldetchar
from laldetchar import cvar as laldetcharcvar
from lal import cvar as lalcvar
print("PASSED module load")

# check object parent tracking
print("checking object parent tracking ...")
a = laldetchar.swig_laldetchar_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = laldetcharcvar.swig_laldetchar_test_parent_map.s
    laldetcharcvar.swig_laldetchar_test_parent_map.s = lalcvar.swig_lal_test_struct_const
del c
del b
del a
lal.CheckMemoryLeaks()
print("PASSED object parent tracking")

# passed all tests!
print("PASSED all tests")
