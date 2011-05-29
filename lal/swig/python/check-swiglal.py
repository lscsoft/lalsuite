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

# check static vector/matrix conversions
if not cvar.swiglal_debug:
    msg("skipping static vector/matrix conversions")
else:
    sts = swiglal_static_test_struct()
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
    assert(not cvar.swiglal_static_test_vector.any())
    assert(not cvar.swiglal_static_test_matrix.any())
    assert(not cvar.swiglal_static_test_enum_vector.any())
    assert(not cvar.swiglal_static_test_enum_matrix.any())
    cvar.swiglal_static_test_vector = cvar.swiglal_static_test_const_vector
    assert((cvar.swiglal_static_test_vector == [1, 2, 4]).all())
    assert(swiglal_static_test_const_vector_getel(2) == 4)
    try:
        swiglal_static_test_const_vector_getel(20)
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
iv = gsl_vector_int_calloc(5)
rv = gsl_vector_calloc(5)
cm = gsl_matrix_complex_float_calloc(4, 6)
check_dynamic_vector_matrix(iv, iv.size, rv, rv.size,
                            cm, cm.size1, cm.size2)
gsl_vector_int_free(iv)
gsl_vector_free(rv)
gsl_matrix_complex_float_free(cm)
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

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
