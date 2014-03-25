# Check SWIG Python module wrapping lalxml
# Author: Karl Wette, 2011, 2012

# check module load
print("checking module load ...")
import lal
import lalxml
from lalxml import cvar as lalxmlcvar
from lal import cvar as lalcvar
print("PASSED module load")

# check object parent tracking
print("checking object parent tracking ...")
a = lalxml.swig_lalxml_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalxmlcvar.swig_lalxml_test_parent_map.s
    lalxmlcvar.swig_lalxml_test_parent_map.s = lalcvar.swig_lal_test_struct_const
del a, b, c
lal.CheckMemoryLeaks()
print("PASSED object parent tracking")

# passed all tests!
print("PASSED all tests")
