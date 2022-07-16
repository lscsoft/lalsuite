# Check SWIG Python bindings for LALBurst
# Author: Karl Wette, 2011--2014

# check module load
print("checking module load ...")
import lal
import lalburst
from lalburst import globalvar as lalburstglobalvar
from lal import globalvar as lalglobalvar
print("PASSED module load")

# set error handlers
def set_nice_error_handlers():
    lal.swig_set_nice_error_handlers()
def set_default_error_handlers():
    lal.swig_set_nasty_error_handlers()
set_default_error_handlers()

# check object parent tracking
print("checking object parent tracking ...")
a = lalburst.swig_lalburst_test_parent_map_struct()
for i in range(0, 7):
    b = a.s
    c = lalburstglobalvar.swig_lalburst_test_parent_map.s
    lalburstglobalvar.swig_lalburst_test_parent_map.s = lalglobalvar.swig_lal_test_struct_const
del c
del b
del a
lal.CheckMemoryLeaks()
print("PASSED object parent tracking")

# passed all tests!
print("PASSED all tests")
