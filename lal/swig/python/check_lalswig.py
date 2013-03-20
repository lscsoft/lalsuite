# Check SWIG Python module wrapping lal
# Author: Karl Wette, 2011, 2012

import datetime
import numpy

# check module load
import lal
from lal import cvar as lalcvar
lalcvar.lalDebugLevel = 1
print("passed module load")

# check memory allocation
lal.CheckMemoryLeaks()
mem1 = lal.Detector()
mem2 = lal.CreateCOMPLEX8Vector(5)
mem3 = lal.CreateREAL8Vector(3)
mem4 = lal.CreateREAL4TimeSeries("test", lal.LIGOTimeGPS(0), 100, 0.1, lalcvar.lalDimensionlessUnit, 10)
print("*** below should be an error message from CheckMemoryLeaks() ***")
try:
    lal.CheckMemoryLeaks()
    raise Exception("expected exception")
except:
    pass
print("*** above should be an error message from CheckMemoryLeaks() ***")
del mem1, mem2, mem3, mem4
lal.CheckMemoryLeaks()
print("passed memory allocation")

# check string conversions
strs = ["a", "bc", "def"]
sv = lal.CreateStringVector(*strs)
assert(sv.length == 3)
assert((sv.data.astype(numpy.object) == strs).all())
strs[0] = "ghijk"
sv.data[0] = strs[0]
strs.append("lmnopq")
lal.AppendString2Vector(sv, strs[3])
assert(sv.length == 4)
for i in range(0, 4):
    assert(sv.data[i] == strs[i])
del sv
lal.CheckMemoryLeaks()
print("passed string conversions")

## check static vector/matrix conversions
lalcvar.lalswig_test_struct_vector[0] = lalcvar.lalswig_test_struct_const
assert(lalcvar.lalswig_test_struct_vector[0].i == lalcvar.lalswig_test_struct_const.i)
assert(lalcvar.lalswig_test_struct_vector[0].f == lalcvar.lalswig_test_struct_const.f)
assert(lalcvar.lalswig_test_struct_vector[0].str == lalcvar.lalswig_test_struct_const.str)
assert((lalcvar.lalswig_test_struct_vector[0].vec == lalcvar.lalswig_test_struct_const.vec).all())
lalcvar.lalswig_test_struct_matrix[0, 0] = lalcvar.lalswig_test_struct_const
assert(lalcvar.lalswig_test_struct_matrix[0, 0].i == lalcvar.lalswig_test_struct_const.i)
assert(lalcvar.lalswig_test_struct_matrix[0, 0].f == lalcvar.lalswig_test_struct_const.f)
assert(lalcvar.lalswig_test_struct_matrix[0, 0].str == lalcvar.lalswig_test_struct_const.str)
assert((lalcvar.lalswig_test_struct_matrix[0, 0].vec == lalcvar.lalswig_test_struct_const.vec).all())
sts = lal.lalswig_test_struct()
assert(len(sts.vec) == 3)
assert(len(sts.evec) == 3)
assert(sts.mat.shape == (2, 3))
sts.vec = [3, 2, 1]
assert((sts.vec == [3, 2, 1]).all())
sts.mat = [[4, 5, 6], (9, 8, 7)]
try:
    sts.mat = [[1.1, 2.3, 4.5], [6.5, 4.3, 2.1]]
    raise Exception("expected exception")
except:
    pass
assert((sts.mat == [[4, 5, 6], [9, 8, 7]]).all())
for i in range(0, 3):
    sts.evec[i] = 2*i + 3
    assert(sts.evec[i] == (2*i + 3))
del sts
assert(not lalcvar.lalswig_test_enum_vector.any())
assert(not lalcvar.lalswig_test_enum_matrix.any())
assert(not lalcvar.lalswig_test_INT4_vector.any())
assert(not lalcvar.lalswig_test_INT4_matrix.any())
assert(not lalcvar.lalswig_test_REAL8_vector.any())
assert(not lalcvar.lalswig_test_REAL8_matrix.any())
assert(not lalcvar.lalswig_test_COMPLEX8_vector.any())
assert(not lalcvar.lalswig_test_COMPLEX8_matrix.any())
lalcvar.lalswig_test_INT4_vector[0] = 10
assert(lalcvar.lalswig_test_INT4_vector[0] == 10)
lalcvar.lalswig_test_INT4_matrix[0, 0] = 11
assert(lalcvar.lalswig_test_INT4_matrix[0, 0] == 11)
lalcvar.lalswig_test_INT4_vector = lalcvar.lalswig_test_INT4_const_vector
assert((lalcvar.lalswig_test_INT4_vector == [1, 2, 4]).all())
assert(lalcvar.lalswig_test_INT4_const_vector[2] == 4)
lalcvar.lalswig_test_INT4_matrix = lalcvar.lalswig_test_INT4_const_matrix
assert((lalcvar.lalswig_test_INT4_matrix == [[1, 2, 4], [2, 4, 8]]).all())
assert(lalcvar.lalswig_test_INT4_const_matrix[1, 2] == 8)
try:
    lalcvar.lalswig_test_INT4_const_vector(20)
    raise Exception("expected exception")
except:
    pass
lalcvar.lalswig_test_REAL8_vector[0] = 3.4
assert(lalcvar.lalswig_test_REAL8_vector[0] == 3.4)
lalcvar.lalswig_test_REAL8_matrix[0, 0] = 5.6
assert(lalcvar.lalswig_test_REAL8_matrix[0, 0] == 5.6)
lalcvar.lalswig_test_COMPLEX8_vector[0] = complex(3.5, 4.75)
assert(lalcvar.lalswig_test_COMPLEX8_vector[0] == complex(3.5, 4.75))
lalcvar.lalswig_test_COMPLEX8_matrix[0, 0] = complex(5.5, 6.25)
assert(lalcvar.lalswig_test_COMPLEX8_matrix[0, 0] == complex(5.5, 6.25))
print("passed static vector/matrix conversions")

# check dynamic vector/matrix conversions
def check_dynamic_vector_matrix(iv, ivl, rv, rvl, cm, cms1, cms2):
    assert(ivl == 5)
    iv.data = [1, 3, 2, 4, 3]
    assert((iv.data == [1, 3, 2, 4, 3]).all())
    iv.data[3] = 7
    assert(iv.data[3] == 7)
    assert(rvl == 5)
    rv.data = [1.2, 3.4, 2.6, 4.8, 3.5]
    assert((rv.data == [1.2, 3.4, 2.6, 4.8, 3.5]).all())
    rv.data[rvl - 1] = 7.5
    assert(rv.data[rvl - 1] == 7.5)
    try:
        rv.data_setel(rvl, 99.9)
        raise Exception("expected exception")
    except:
        pass
    try:
        iv.data = rv.data
        raise Exception("expected exception")
    except:
        pass
    rv.data = iv.data
    assert((rv.data == iv.data).all())
    assert(cms1 == 4)
    assert(cms2 == 6)
    for i in range(0, cms1):
        for j in range(0, cms2):
            cm.data[i, j] = complex(i / 4.0, j / 2.0)
    assert(cm.data[2, 3] == complex(0.5, 1.5))
    assert(cm.data[3, 2] == complex(0.75, 1.0))
    try:
        iv.data_setel(0, cm.data_getel(2, 3))
        raise Exception("expected exception")
    except:
        pass
    try:
        iv.data[0] = cm.data[2, 3]
        raise Exception("expected exception")
    except:
        pass
    try:
        rv.data_setel(0, cm.data_getel(3, 2))
        raise Exception("expected exception")
    except:
        pass
    try:
        rv.data[0] = cm.data[3, 2]
        raise Exception("expected exception")
    except:
        pass
# check LAL vector and matrix datatypes
iv = lal.CreateINT4Vector(5)
rv = lal.CreateREAL8Vector(5)
cm = lal.CreateCOMPLEX8VectorSequence(4, 6)
check_dynamic_vector_matrix(iv, iv.length, rv, rv.length,
                            cm, cm.length, cm.vectorLength)
del iv, rv, cm
lal.CheckMemoryLeaks()
print("passed dynamic vector/matrix conversions (LAL)")
# check GSL vectors and matrices
iv = lal.gsl_vector_int(5)
rv = lal.gsl_vector(5)
cm = lal.gsl_matrix_complex_float(4, 6)
check_dynamic_vector_matrix(iv, iv.size, rv, rv.size,
                            cm, cm.size1, cm.size2)
del iv, rv, cm
print("passed dynamic vector/matrix conversions (GSL)")

# check 'tm' struct conversions
gps = 989168284
utc = [2011, 5, 11, 16, 57, 49, 2, 131, 0]
assert(lal.GPSToUTC(gps) == utc)
assert(lal.UTCToGPS(utc) == gps)
assert(lal.UTCToGPS(utc[0:6]) == gps)
utc[6] = utc[7] = 0
for i in [-1, 0, 1]:
    utc[8] = i
    assert(lal.UTCToGPS(utc) == gps)
utcd = utc
for i in range(0, 10):
    utcd[2] = utc[2] + i
    utcd = lal.GPSToUTC(lal.UTCToGPS(utcd))
    dt = datetime.datetime(*utcd[0:6])
    assert(utcd[6] == dt.weekday())
lal.CheckMemoryLeaks()
print("passed 'tm' struct conversions")

# check LIGOTimeGPS operations
from lal import LIGOTimeGPS
t0 = LIGOTimeGPS()
assert(t0 == 0 and isinstance(t0, LIGOTimeGPS))
assert(t0 != None and not t0 is None)
t1 = LIGOTimeGPS(10.5)
t2 = LIGOTimeGPS(10, 500000000)
assert(not t0 and t1 and t2)
assert(t1 == t2 and isinstance(t1, LIGOTimeGPS))
t3 = +t1
t3 = -t2
assert(t1 == t2 and t1 >= t2 and t2 >= t1)
assert(abs(-t1) == t1)
assert(float(t1) == 10.5)
assert(t1 + 3.5 == 14 and isinstance(t1 + 3.5, LIGOTimeGPS))
t2 -= 5.5
assert(t2 == 5 and isinstance(t2, LIGOTimeGPS))
assert(t2 + 5.5 >= t1 and t2 + 3 != t2)
assert(t2 - 5 == t0 and isinstance(t2 - 5, LIGOTimeGPS))
assert(t1 * 3 == 31.5 and isinstance(t1 * 3, LIGOTimeGPS))
assert(t2 / 2.5 == 2 and isinstance(t2 / 2.5, LIGOTimeGPS))
assert(t1 + t2 == 15.5 and isinstance(t1 + t2, LIGOTimeGPS))
assert(t1 - t2 == 5.5 and isinstance(t1 - t2, LIGOTimeGPS))
assert(t1 * t2 == 52.5 and isinstance(t1 * t2, LIGOTimeGPS))
assert(t2 * t1 == 52.5 and isinstance(t2 * t1, LIGOTimeGPS))
assert(t1 / t2 == 2.1 and isinstance(t1 / t2, LIGOTimeGPS))
assert(t1 % t2 == 0.5 and isinstance(t1 % t2, LIGOTimeGPS))
assert(t1 > t2 and t2 < t1 and t1 >= t2 and t2 <= t1)
assert(LIGOTimeGPS(333333333,333333333) == LIGOTimeGPS(1000000000) / 3)
assert(LIGOTimeGPS(666666666,666666667) == LIGOTimeGPS(2000000000) / 3)
t1 += 812345667.75
assert(str(t1) == "812345678.250000000")
assert(LIGOTimeGPS(repr(t1)) == t1)
assert(long(t1) == 812345678)
assert(lal.GPSToINT8NS(t1) == 812345678250000000L)
assert(hash(t1) == 1049484238)
t4struct = lal.lalswig_test_gps()
t4struct.t = 1234.5
assert(t4struct.t == 1234.5)
del t0, t1, t2, t3, t4struct
lal.CheckMemoryLeaks()
print("passed LIGOTimeGPS operations")

# check object parent tracking
a = lal.gsl_vector(3)
a.data = [1.1, 2.2, 3.3]
b = a.data
assert(not b.flags['OWNDATA'])
assert((b == [1.1, 2.2, 3.3]).all())
del a
assert((b == [1.1, 2.2, 3.3]).all())
ts = lal.CreateREAL8TimeSeries("test", lal.LIGOTimeGPS(0), 0, 0.1, lalcvar.lalDimensionlessUnit, 10)
ts.data.data = range(0, 10)
for i in range(0, 7):
    v = ts.data
assert((v.data == range(0, 10)).all())
del ts
assert((v.data == range(0, 10)).all())
del v
lal.CheckMemoryLeaks()
print("passed object parent tracking")

# passed all tests!
print("PASSED all tests")
