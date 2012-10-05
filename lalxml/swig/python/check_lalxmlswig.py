# Check SWIG Python module wrapping lalxml
# Author: Karl Wette, 2011, 2012

# check module load
import lal
import lalxml
from lalxml import cvar as lalxmlcvar
from lal import cvar as lalcvar
lalcvar.lalDebugLevel = 1
print("passed module load")

# check object parent tracking
a = lalxml.lalxmlswig_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalxmlcvar.lalxmlswig_test_parent_map.s
    lalxmlcvar.lalxmlswig_test_parent_map.s = lalcvar.lalswig_test_struct_const
del a, b, c
lal.CheckMemoryLeaks()
print("passed object parent tracking")

# passed all tests!
print("PASSED all tests")
