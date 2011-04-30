# check SWIG Python module wrapping the LAL library
# Author: Karl Wette, 2011

import os
import datetime

def msg(str):
    print(os.path.basename(__file__) + ": " + str)
class error(Exception):
    def __init__(self, str):
        self.str = str
    def __str__(self):
        return str

# check memory allocation
if not cvar.swiglal_debug:
    msg("skipping memory allocation")
else:
    LALCheckMemoryLeaks()
    mem1 = LALDetector()
    mem2 = LALStringVector()
    mem3 = COMPLEX8Vector()
    mem4 = XLALCreateREAL8Vector(3)
    msg("*** below should be an error message from LALCheckMemoryLeaks() ***")
    try:
        LALCheckMemoryLeaks()
        raise error("expected exception")
    except:
        pass
    msg("*** above should be an error message from LALCheckMemoryLeaks() ***")
    del mem1, mem2, mem3
    XLALDestroyREAL8Vector(mem4)
    LALCheckMemoryLeaks()
    msg("passed memory allocation")

# check complex number conversions
assert(XLALCOMPLEX8Add(complex(1, 3), complex(2, -5)) == complex(3, -2))
assert(XLALCOMPLEX8Sub(complex(4, 2), complex(10, 5)) == complex(-6, -3))
assert(XLALCOMPLEX16Mul(complex(10, -7), complex(30, 17)) == complex(419, -40))
assert(XLALCOMPLEX16Div(complex(111.75, -120.50), complex(5, 12)) == complex(-5.25, -11.5))
msg("passed complex number conversions")

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
