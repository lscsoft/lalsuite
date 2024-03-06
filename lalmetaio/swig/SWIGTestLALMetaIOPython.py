# Check SWIG Python bindings for LALMetaIO
# Author: Karl Wette, 2011--2014

import os
import sys
import pytest

# check module load
print("checking module load ...", file=sys.stderr)
import lal
import lalmetaio
from lalmetaio import globalvar as lalmetaioglobalvar
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


# -- tests


def test_object_parent_tracking():
    """check object parent tracking"""

    print("checking object parent tracking ...", file=sys.stderr)
    a = lalmetaio.swig_lalmetaio_test_parent_map_struct()
    for i in range(0, 7):
        b = a.s
        c = lalmetaioglobalvar.swig_lalmetaio_test_parent_map.s
        lalmetaioglobalvar.swig_lalmetaio_test_parent_map.s = (
            lalglobalvar.swig_lal_test_struct_const
        )
    del c
    del b
    del a
    locals()  # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED object parent tracking", file=sys.stderr)


if __name__ == "__main__":
    args = sys.argv[1:] or [
        "-v",
        "-rs",
        "--junit-xml=junit-SWIGTestLALMetaIOPython.xml",
    ]
    sys.exit(pytest.main(args=[__file__] + args))
