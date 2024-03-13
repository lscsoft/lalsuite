# Check SWIG Python bindings for LALInference
# Author: Karl Wette, 2011--2014

import os
import sys
import inspect
import gc

import pytest

# check module load
print("checking module load ...", file=sys.stderr)
import lal
import lalinference
from lalinference import globalvar as lalinferenceglobalvar
from lal import globalvar as lalglobalvar

print("PASSED module load", file=sys.stderr)


# -- configure error handling

# set error handlers
def set_nice_error_handlers():
    lal.swig_set_nice_error_handlers()


def set_default_error_handlers():
    if "NASTY_ERROR_HANDLERS" in os.environ:
        lal.swig_set_nasty_error_handlers()
    else:
        lal.swig_set_nice_error_handlers()


set_default_error_handlers()


# -- check for memory leaks


def check_memory_leaks():
    # pytest's rewrite of assert() can keep references
    # to SWIGLAL objects around; find and clear them
    frame = inspect.currentframe()
    try:
        for v in frame.f_back.f_locals:
            if v.startswith("@py_assert"):
                frame.f_back.f_locals[v] = None
    finally:
        del frame

    # garbage collector should free all SWIGLAL objects
    gc.collect()

    # check that all LAL memory has been freed
    lal.CheckMemoryLeaks()


# -- tests


def test_object_parent_tracking():
    """check object parent tracking"""

    print("checking object parent tracking ...", file=sys.stderr)
    a = lalinference.swig_lalinference_test_parent_map_struct()
    for i in range(0, 7):
        b = a.s
        c = lalinferenceglobalvar.swig_lalinference_test_parent_map.s
        lalinferenceglobalvar.swig_lalinference_test_parent_map.s = (
            lalglobalvar.swig_lal_test_struct_const
        )
    del c
    del b
    del a
    check_memory_leaks()
    print("PASSED object parent tracking", file=sys.stderr)


if __name__ == "__main__":
    args = sys.argv[1:] or [
        "-v",
        "-rs",
        "--junit-xml=junit-SWIGTestLALInferencePython.xml",
    ]
    sys.exit(pytest.main(args=[__file__] + args))
