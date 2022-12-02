# Check SWIG Python bindings for LALPulsar
# Author: Karl Wette, 2011--2014

import os
import sys
import pytest

# check module load
print("checking module load ...", file=sys.stderr)
import lal
import lalpulsar
from lalpulsar import globalvar as lalpulsarglobalvar
from lal import globalvar as lalglobalvar
print("PASSED module load", file=sys.stderr)

# set error handlers
def set_nice_error_handlers():
    lal.swig_set_nice_error_handlers()
def set_default_error_handlers():
    if 'NASTY_ERROR_HANDLERS' in os.environ:
        lal.swig_set_nasty_error_handlers()
    else:
        lal.swig_set_nice_error_handlers()
set_default_error_handlers()

def test_object_parent_tracking():
    """check object parent tracking
    """

    print("checking object parent tracking ...", file=sys.stderr)
    a = lalpulsar.swig_lalpulsar_test_parent_map_struct()
    for i in range(0, 7):
        b = a.s
        c = lalpulsarglobalvar.swig_lalpulsar_test_parent_map.s
        lalpulsarglobalvar.swig_lalpulsar_test_parent_map.s = lalglobalvar.swig_lal_test_struct_const
    del c
    del b
    del a
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED object parent tracking", file=sys.stderr)

def test_array_element_assignment():
    """check array element assignment
    """

    print("checking array element assignment ...", file=sys.stderr)
    mts = lalpulsar.CreateMultiLIGOTimeGPSVector(2)
    ts0 = lalpulsar.CreateTimestampVector(3)
    for i in range(ts0.length):
        ts0.data[i] = lal.LIGOTimeGPS(900000000 + i)
    mts.data[0] = ts0
    lal.swig_set_nasty_error_handlers()
    del mts
    del ts0
    set_default_error_handlers()
    lal.CheckMemoryLeaks()
    print("PASSED array element assignment", file=sys.stderr)

def test_array_element_parent_tracking():
    """check array element parent tracking
    """

    print("checking array element parent tracking", file=sys.stderr)
    assert lalpulsar.MakeTimestamps(0, 100, 1, 0).data[-1] == lal.LIGOTimeGPS(99)
    assert lalpulsar.MakeMultiTimestamps(0, 100, 1, 0, 3).data[-1].data[-1] == lal.LIGOTimeGPS(99)
    lal.CheckMemoryLeaks()
    print("PASSED array element parent tracking", file=sys.stderr)

def test_CWMFDataParams_usage():
    """check CWMFDataParams usage
    """

    print("checking CWMFDataParams usage", file=sys.stderr)
    def make_CWMFDataParams():
        cwmf = lalpulsar.CWMFDataParams()
        cwmf.multiTimestamps = lalpulsar.MakeMultiTimestamps(0, 10, 1, 0, 3)
        return cwmf
    assert make_CWMFDataParams().multiTimestamps.data[-1].data[-1] == lal.LIGOTimeGPS(9)
    lal.CheckMemoryLeaks()
    print("PASSED CWMFDataParams usage", file=sys.stderr)

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-SWIGTestLALPulsarPython.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
