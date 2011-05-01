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

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
