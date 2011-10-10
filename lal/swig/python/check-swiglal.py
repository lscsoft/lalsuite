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

# check module load
import swiglal
from swiglal import *
swiglal.cvar.lalDebugLevel = 1
msg("passed module load")

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

# check string conversions
strs = ["a", "bc", "def"]
sv = XLALCreateStringVector(*strs)
assert(sv.length == 3)
assert((sv.data == strs).all())
strs[0] = "ghijk"
sv.data_setel(0, strs[0])
strs.append("lmnopq")
XLALAppendString2Vector(sv, strs[3])
assert(sv.length == 4)
for i in range(0,4):
    assert(sv.data_getel(i) == strs[i])
XLALDestroyStringVector(sv)
msg("passed string conversions")

# check vector/matrix struct type accessors
if not cvar.swiglal_debug:
    msg("skipping vector/matrix struct type accessors")
else:
    swiglal_test_struct_vector_setel(0, cvar.swiglal_test_struct_const)
    assert(swiglal_test_struct_vector_getel(0).a == cvar.swiglal_test_struct_const.a)
    assert(swiglal_test_struct_vector_getel(0).b == cvar.swiglal_test_struct_const.b)
    assert(swiglal_test_struct_vector_getel(0).c == cvar.swiglal_test_struct_const.c)
    swiglal_test_struct_matrix_setel(0, 0, cvar.swiglal_test_struct_const)
    assert(swiglal_test_struct_matrix_getel(0, 0).a == cvar.swiglal_test_struct_const.a)
    assert(swiglal_test_struct_matrix_getel(0, 0).b == cvar.swiglal_test_struct_const.b)
    assert(swiglal_test_struct_matrix_getel(0, 0).c == cvar.swiglal_test_struct_const.c)
    msg("passed vector/matrix struct type accessors")

# check static vector/matrix conversions
if not cvar.swiglal_debug:
    msg("skipping static vector/matrix conversions")
else:
    sts = swiglal_test_static_struct()
    assert(len(sts.vector) == 3)
    assert(len(sts.enum_vector) == 3)
    assert(sts.matrix.shape == (2, 3))
    assert(sts.enum_matrix.shape == (2, 3))
    sts.vector = [3, 2, 1]
    assert((sts.vector == [3, 2, 1]).all())
    sts.matrix = [[4, 5, 6], (9, 8, 7)]
    try:
        sts.matrix = [[1.1, 2.3, 4.5], [6.5, 4.3, 2.1]]
        raise error("expected exception")
    except:
        pass
    assert((sts.matrix == [[4, 5, 6], [9, 8, 7]]).all())
    for i in range(0,3):
        sts.enum_vector_setel(i, 2*i + 3)
        assert(sts.enum_vector_getel(i) == (2*i + 3))
    del sts
    assert(not cvar.swiglal_test_static_vector.any())
    assert(not cvar.swiglal_test_static_matrix.any())
    assert(not cvar.swiglal_test_static_enum_vector.any())
    assert(not cvar.swiglal_test_static_enum_matrix.any())
    swiglal_test_static_vector_setel(0, 10)
    assert(swiglal_test_static_vector_getel(0) == 10)
    swiglal_test_static_matrix_setel(0, 0, 11)
    assert(swiglal_test_static_matrix_getel(0, 0) == 11)
    cvar.swiglal_test_static_vector = cvar.swiglal_test_static_const_vector
    assert((cvar.swiglal_test_static_vector == [1, 2, 4]).all())
    assert(swiglal_test_static_const_vector_getel(2) == 4)
    cvar.swiglal_test_static_matrix = cvar.swiglal_test_static_const_matrix
    assert((cvar.swiglal_test_static_matrix == [[1, 2, 4], [2, 4, 8]]).all())
    assert(swiglal_test_static_const_matrix_getel(1, 2) == 8)
    try:
        swiglal_test_static_const_vector_getel(20)
        raise error("expected exception")
    except:
        pass
    msg("passed static vector/matrix conversions")

# check dynamic vector/matrix conversions
def check_dynamic_vector_matrix(iv, ivl, rv, rvl, cm, cms1, cms2):
    assert(ivl == 5)
    iv.data = [1, 3, 2, 4, 3]
    assert((iv.data == [1, 3, 2, 4, 3]).all())
    iv.data_setel(3, 7)
    assert(iv.data_getel(3) == 7)
    assert(rvl == 5)
    rv.data = [1.2, 3.4, 2.6, 4.8, 3.5]
    assert((rv.data == [1.2, 3.4, 2.6, 4.8, 3.5]).all())
    rv.data_setel(rvl - 1, 7.5)
    assert(rv.data_getel(rvl - 1) == 7.5)
    try:
        rv.data_setel(rvl, 99.9)
        raise error("expected exception")
    except:
        pass
    try:
        iv.data = rv.data
        raise error("expected exception")
    except:
        pass
    rv.data = iv.data
    assert((rv.data == iv.data).all())
    assert(cms1 == 4)
    assert(cms2 == 6)
    for i in range(0,cms1):
        for j in range(0,cms2):
            cm.data_setel(i, j, complex(i / 4.0, j / 2.0))
    assert(cm.data_getel(2, 3) == complex(0.5, 1.5))
    assert(cm.data_getel(3, 2) == complex(0.75, 1.0))
    try:
        iv.data_setel(0, cm.data_getel(2, 3))
        raise error("expected exception")
    except:
        pass
    try:
        rv.data_setel(0, cm.data_getel(3, 2))
        raise error("expected exception")
    except:
        pass
# check LAL vector and matrix datatypes
iv = XLALCreateINT4Vector(5)
rv = XLALCreateREAL8Vector(5)
cm = XLALCreateCOMPLEX8VectorSequence(4, 6)
check_dynamic_vector_matrix(iv, iv.length, rv, rv.length,
                            cm, cm.length, cm.vectorLength)
XLALDestroyINT4Vector(iv)
XLALDestroyREAL8Vector(rv)
XLALDestroyCOMPLEX8VectorSequence(cm)
LALCheckMemoryLeaks()
msg("passed dynamic vector/matrix conversions (LAL)")
# check GSL vectors and matrices
iv = gsl_vector_int(5)
rv = gsl_vector(5)
cm = gsl_matrix_complex_float(4, 6)
check_dynamic_vector_matrix(iv, iv.size, rv, rv.size,
                            cm, cm.size1, cm.size2)
del iv, rv, cm
msg("passed dynamic vector/matrix conversions (GSL)")

# check 'tm' struct conversions
gps = 989168284
utc = [2011, 5, 11, 16, 57, 49, 2, 131, 0]
assert(XLALGPSToUTC(None, gps) == utc)
assert(XLALUTCToGPS(utc) == gps)
assert(XLALUTCToGPS(utc[0:6]) == gps)
utc[6] = utc[7] = 0
for i in [-1, 0, 1]:
    utc[8] = i
    assert(XLALUTCToGPS(utc) == gps)
utcd = utc
for i in range(0,10):
    utcd[2] = utc[2] + i
    utcd = XLALGPSToUTC(None, XLALUTCToGPS(utcd))
    dt = datetime.datetime(*utcd[0:6])
    assert(utcd[6] == dt.weekday())
msg("passed 'tm' struct conversions")

# check LIGOTimeGPS operations
t0 = LIGOTimeGPS()
assert(t0 == 0 and isinstance(t0, LIGOTimeGPS))
t1 = LIGOTimeGPS(10.5)
t2 = LIGOTimeGPS(10, 500000000)
assert(not t0 and t1 and t2)
assert(t1 == t2 and isinstance(t1, LIGOTimeGPS))
+t1
-t2
assert(t1 == t2 and t1 >= t2 and t2 >= t1)
assert(abs(-t1) == t1)
assert(int(t1) == 10 and int(t1) == int(t2))
assert(float(t1) == 10.5)
assert(t1 + 3.5 == 14 and isinstance(t1 + 3.5, LIGOTimeGPS))
assert(3.5 + t1 == 14 and isinstance(3.5 + t1, LIGOTimeGPS))
t2 -= 5.5
assert(t2 == 5 and isinstance(t2, LIGOTimeGPS))
assert(t2 + 5.5 >= t1 and t2 + 3 != t2)
assert(t2 - 5 == t0 and isinstance(t2 - 5, LIGOTimeGPS))
assert(t1 * 3 == 31.5 and isinstance(t1 * 3, LIGOTimeGPS))
assert(3 * t1 == 31.5 and isinstance(3 * t1, LIGOTimeGPS))
assert(t2 / 2.5 == 2 and isinstance(t2 / 2.5, LIGOTimeGPS))
assert(21 / t1  == 2 and isinstance(21 / t1 , LIGOTimeGPS))
assert(t1 + t2 == 15.5 and isinstance(t1 + t2, LIGOTimeGPS))
assert(t1 - t2 == 5.5 and isinstance(t1 - t2, LIGOTimeGPS))
assert(t1 * t2 == 52.5 and isinstance(t1 * t2, LIGOTimeGPS))
assert(t2 * t1 == 52.5 and isinstance(t2 * t1, LIGOTimeGPS))
assert(t1 / t2 == 2.1 and isinstance(t1 / t2, LIGOTimeGPS))
assert(t1 % t2 == 0.5 and isinstance(t1 % t2, LIGOTimeGPS))
assert(t1 // t2 == 2.0 and isinstance(t1 // t2, LIGOTimeGPS))
assert(t1 > t2 and t2 < t1 and t1 >= t2 and t2 <= t1)
class tc:
    seconds = t1.gpsSeconds
    nanoseconds = t1.gpsNanoSeconds
assert(t1 == tc)
assert(t2 + tc == 15.5 and isinstance(t2 + tc, LIGOTimeGPS))
assert(t2 - tc == -5.5 and isinstance(t2 - tc, LIGOTimeGPS))
assert(t2 * tc == 52.5 and isinstance(t2 * tc, LIGOTimeGPS))
assert(t2 < tc and t2 <= tc)
t1 += 812345667.75
assert(str(t1) == "812345678.250000000")
assert(LIGOTimeGPS(str(t1)) == t1)
assert(long(t1) == 812345678)
assert(t1.ns() == 812345678250000000L)
assert(hash(t1) == 1049484238)
msg("passed LIGOTimeGPS operations")

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
