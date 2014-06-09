# Check SWIG Python bindings for LAL
# Author: Karl Wette, 2011--2014

import sys
import warnings
import datetime
import numpy
expected_exception = False

# turn NumPy's ComplexWarning into an error, if available
if hasattr(numpy, "ComplexWarning"):
    warnings.simplefilter("error", numpy.ComplexWarning)

# check module load
print("checking module load ...")
import lal
from lal import globalvar as lalglobalvar
print("PASSED module load")

# check memory allocation
print("checking memory allocation ...")
if not lal.NoDebug:
    lal.CheckMemoryLeaks()
    mem1 = lal.Detector()
    mem2 = lal.CreateCOMPLEX8Vector(5)
    mem3 = lal.CreateREAL8Vector(3)
    mem4 = lal.CreateREAL4TimeSeries("test", lal.LIGOTimeGPS(0), 100, 0.1, lal.DimensionlessUnit, 10)
    print("*** below should be an error message from CheckMemoryLeaks() ***")
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        lal.CheckMemoryLeaks()
        expected_exception = True
    except:
        pass
    assert(not expected_exception)
    print("*** above should be an error message from CheckMemoryLeaks() ***")
    del mem1
    del mem2
    del mem3
    del mem4
    lal.CheckMemoryLeaks()
    print("PASSED memory allocation")
else:
    print("skipped memory allocation")

# check object parent tracking
print("checking object parent tracking ...")
a = lal.gsl_vector(3)
a.data = [1.1, 2.2, 3.3]
b = a.data
assert(not b.flags['OWNDATA'])
assert((b == [1.1, 2.2, 3.3]).all())
del a
assert((b == [1.1, 2.2, 3.3]).all())
ts = lal.CreateREAL8TimeSeries("test", lal.LIGOTimeGPS(0), 0, 0.1, lal.DimensionlessUnit, 10)
ts.data.data = range(0, 10)
for i in range(0, 7):
    v = ts.data
assert((v.data == range(0, 10)).all())
del ts
assert((v.data == range(0, 10)).all())
del v
lal.CheckMemoryLeaks()
print("PASSED object parent tracking")

# check equal return/first argument type handling
print("checking equal return/first argument type handling")
sv = lal.CreateStringVector("1")
assert(sv.length == 1)
lal.AppendString2Vector(sv, "2")
assert(sv.length == 2)
sv = lal.AppendString2Vector(sv, "3")
assert(sv.length == 3)
sv2 = lal.AppendString2Vector(sv, "4")
assert(sv.length == 4)
assert(sv2.length == 4)
assert(sv == sv2)
del sv
del sv2
lal.CheckMemoryLeaks()
ts = lal.CreateREAL8TimeSeries("ts", 800000000, 100, 0.1, lal.HertzUnit, 10)
assert(ts.data.length == 10)
lal.ResizeREAL8TimeSeries(ts, 0, 20)
assert(ts.data.length == 20)
ts = lal.ResizeREAL8TimeSeries(ts, 0, 30)
assert(ts.data.length == 30)
ts2 = lal.ResizeREAL8TimeSeries(ts, 0, 40)
assert(ts.data.length == 40)
assert(ts2.data.length == 40)
assert(ts == ts2)
del ts
del ts2
lal.CheckMemoryLeaks()
print("PASSED equal return/first argument type handling")

# check string conversions
print("checking string conversions ...")
strs = ["a", "bc", "def"]
sv = lal.CreateStringVector(*strs)
assert(sv.length == 3)
assert((sv.data.astype(numpy.object) == strs).all())
strs[0] = "ghijk"
sv.data[0] = strs[0]
strs.append("lmnopq")
sv = lal.AppendString2Vector(sv, strs[3])
assert(sv.length == 4)
for i in range(0, 4):
    assert(sv.data[i] == strs[i])
del sv
lal.CheckMemoryLeaks()
print("PASSED string conversions")

# check static vector/matrix conversions
print("checking static vector/matrix conversions ...")
lalglobalvar.swig_lal_test_struct_vector[0] = lalglobalvar.swig_lal_test_struct_const
assert(lalglobalvar.swig_lal_test_struct_vector[0].n == lalglobalvar.swig_lal_test_struct_const.n)
assert(lalglobalvar.swig_lal_test_struct_vector[0].i == lalglobalvar.swig_lal_test_struct_const.i)
assert(lalglobalvar.swig_lal_test_struct_vector[0].f == lalglobalvar.swig_lal_test_struct_const.f)
assert(lalglobalvar.swig_lal_test_struct_vector[0].str == lalglobalvar.swig_lal_test_struct_const.str)
assert((lalglobalvar.swig_lal_test_struct_vector[0].vec == lalglobalvar.swig_lal_test_struct_const.vec).all())
lalglobalvar.swig_lal_test_struct_matrix[0, 0] = lalglobalvar.swig_lal_test_struct_const
assert(lalglobalvar.swig_lal_test_struct_matrix[0, 0].n == lalglobalvar.swig_lal_test_struct_const.n)
assert(lalglobalvar.swig_lal_test_struct_matrix[0, 0].i == lalglobalvar.swig_lal_test_struct_const.i)
assert(lalglobalvar.swig_lal_test_struct_matrix[0, 0].f == lalglobalvar.swig_lal_test_struct_const.f)
assert(lalglobalvar.swig_lal_test_struct_matrix[0, 0].str == lalglobalvar.swig_lal_test_struct_const.str)
assert((lalglobalvar.swig_lal_test_struct_matrix[0, 0].vec == lalglobalvar.swig_lal_test_struct_const.vec).all())
sts = lal.swig_lal_test_struct()
assert(len(sts.vec) == 3)
assert(len(sts.evec) == 3)
assert(sts.mat.shape == (2, 3))
sts.vec = [3, 2, 1]
assert((sts.vec == [3, 2, 1]).all())
sts.mat = [[4, 5, 6], (9, 8, 7)]
try:
    sts.mat = [[1.1, 2.3, 4.5], [6.5, 4.3, 2.1]]
    expected_exception = True
except:
    pass
assert(not expected_exception)
assert((sts.mat == [[4, 5, 6], [9, 8, 7]]).all())
for i in range(0, 3):
    sts.evec[i] = 2*i + 3
    assert(sts.evec[i] == (2*i + 3))
del sts
assert(not lalglobalvar.swig_lal_test_enum_vector.any())
assert(not lalglobalvar.swig_lal_test_enum_matrix.any())
assert(len(lalglobalvar.swig_lal_test_empty_INT4_vector) == 0)
assert(not lalglobalvar.swig_lal_test_INT4_vector.any())
assert(not lalglobalvar.swig_lal_test_INT4_matrix.any())
assert(not lalglobalvar.swig_lal_test_REAL8_vector.any())
assert(not lalglobalvar.swig_lal_test_REAL8_matrix.any())
assert(not lalglobalvar.swig_lal_test_COMPLEX8_vector.any())
assert(not lalglobalvar.swig_lal_test_COMPLEX8_matrix.any())
lalglobalvar.swig_lal_test_INT4_vector[0] = 10
assert(lalglobalvar.swig_lal_test_INT4_vector[0] == 10)
lalglobalvar.swig_lal_test_INT4_matrix[0, 0] = 11
assert(lalglobalvar.swig_lal_test_INT4_matrix[0, 0] == 11)
lalglobalvar.swig_lal_test_INT4_vector = lalglobalvar.swig_lal_test_INT4_const_vector
assert((lalglobalvar.swig_lal_test_INT4_vector == [1, 2, 4]).all())
assert(lalglobalvar.swig_lal_test_INT4_const_vector[2] == 4)
lalglobalvar.swig_lal_test_INT4_matrix = lalglobalvar.swig_lal_test_INT4_const_matrix
assert((lalglobalvar.swig_lal_test_INT4_matrix == [[1, 2, 4], [2, 4, 8]]).all())
assert(lalglobalvar.swig_lal_test_INT4_const_matrix[1, 2] == 8)
try:
    lalglobalvar.swig_lal_test_INT4_const_vector(20)
    expected_exception = True
except:
    pass
assert(not expected_exception)
lalglobalvar.swig_lal_test_REAL8_vector[0] = 3.4
assert(lalglobalvar.swig_lal_test_REAL8_vector[0] == 3.4)
lalglobalvar.swig_lal_test_REAL8_matrix[0, 0] = 5.6
assert(lalglobalvar.swig_lal_test_REAL8_matrix[0, 0] == 5.6)
lalglobalvar.swig_lal_test_COMPLEX8_vector[0] = complex(3.5, 4.75)
assert(lalglobalvar.swig_lal_test_COMPLEX8_vector[0] == complex(3.5, 4.75))
lalglobalvar.swig_lal_test_COMPLEX8_matrix[0, 0] = complex(5.5, 6.25)
assert(lalglobalvar.swig_lal_test_COMPLEX8_matrix[0, 0] == complex(5.5, 6.25))
print("PASSED static vector/matrix conversions")

# check dynamic vector/matrix conversions
print("checking dynamic vector/matrix conversions ...")
def check_dynamic_vector_matrix(iv, ivl, rv, rvl, cm, cms1, cms2):
    expected_exception = False
    iv.data = numpy.zeros(ivl, dtype=iv.data.dtype)
    rv.data = numpy.zeros(rvl, dtype=rv.data.dtype)
    cm.data = numpy.zeros((cms1, cms2), dtype=cm.data.dtype)
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
        rv.data[rvl] = 99.9
        expected_exception = True
    except:
        pass
    assert(not expected_exception)
    try:
        iv.data = rv.data
        expected_exception = True
    except:
        pass
    assert(not expected_exception)
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
        iv.data[0] = cm.data[2, 3]
        if not hasattr(numpy, "ComplexWarning"):
            raise Exception("NumPy %s does not have ComplexWarning" % numpy.__version__)
        expected_exception = True
    except:
        pass
    assert(not expected_exception)
    try:
        rv.data[0] = cm.data[3, 2]
        if not hasattr(numpy, "ComplexWarning"):
            raise Exception("NumPy %s does not have ComplexWarning" % numpy.__version__)
        expected_exception = True
    except:
        pass
    assert(not expected_exception)
# check LAL vector and matrix datatypes
iv = lal.CreateINT4Vector(5)
rv = lal.CreateREAL8Vector(5)
cm = lal.CreateCOMPLEX8VectorSequence(4, 6)
check_dynamic_vector_matrix(iv, iv.length, rv, rv.length,
                            cm, cm.length, cm.vectorLength)
del iv
del rv
del cm
rv0 = lal.CreateREAL8Vector(0)
assert(rv0.length == 0)
assert(len(rv0.data) == 0)
del rv0
rv1 = lal.CreateREAL8Vector(1)
rv1.data[0] = 1
del rv1
lal.CheckMemoryLeaks()
print("PASSED dynamic vector/matrix conversions (LAL)")
# check GSL vectors and matrices
iv = lal.gsl_vector_int(5)
rv = lal.gsl_vector(5)
cm = lal.gsl_matrix_complex_float(4, 6)
check_dynamic_vector_matrix(iv, iv.size, rv, rv.size,
                            cm, cm.size1, cm.size2)
del iv
del rv
del cm
rv1 = lal.gsl_vector(1)
rv1.data[0] = 1
del rv1
print("PASSED dynamic vector/matrix conversions (GSL)")

# check fixed and dynamic arrays typemaps
print("checking fixed and dynamic arrays typemaps ...")
a1in = numpy.array([1.2, 3.5, 7.9], dtype=numpy.double)
a1out = a1in * 2.5
assert((lal.swig_lal_test_copyin_array1(a1in, 2.5) == a1out).all())
a2in = numpy.array([[3,2], [7,6], [12,10]], dtype=numpy.int32)
a2out = a2in * 15
assert((lal.swig_lal_test_copyin_array2(a2in, 15) == a2out).all())
try:
    lal.swig_lal_test_copyin_array1(numpy.array([0,0,0,0], dtype=numpy.double), 0)
    expected_exception = True
except:
    pass
assert(not expected_exception)
try:
    lal.swig_lal_test_copyin_array2(numpy.array([[1.2,3.4],[0,0],[0,0]], dtype=numpy.double), 0)
    expected_exception = True
except:
    pass
assert(not expected_exception)
print("PASSED fixed and dynamic arrays typemaps")

## check input views of array structs
print("checking input views of array structs ...")
r4dat = numpy.array([1.2, 2.3, 3.4, 4.5, 5.6], dtype=numpy.float32)
r8dat = numpy.array([3.4, 4.5, 5.6, 6.7, 7.8, 8.9], dtype=numpy.float64)
c8dat = numpy.array(numpy.vectorize(complex)(r4dat, 8 + r4dat), dtype=numpy.complex64)
c16dat = numpy.array(numpy.vectorize(complex)(r8dat, 16 + r8dat), dtype=numpy.complex128)
r4 = lal.CreateREAL4Vector(len(r4dat))
r4.data = r4dat
r4out = lal.CreateREAL4Vector(len(r4dat))
r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
assert(lal.swig_lal_test_viewin_REAL4Vector(r4out, r4))
assert((r4out.data == r4.data).all())
r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
assert(lal.swig_lal_test_viewin_REAL4Vector(r4out, r4dat))
assert((r4out.data == r4dat).all())
r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
assert(lal.swig_lal_test_viewinout_REAL4Vector(r4out, r4))
assert((2 * r4out.data == r4.data).all())
r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
assert(lal.swig_lal_test_viewinout_REAL4Vector(r4out, r4dat))
assert((2 * r4out.data == r4dat).all())
del r4
del r4out
del r4dat
lal.CheckMemoryLeaks()
r8 = lal.CreateREAL8Vector(len(r8dat))
r8.data = r8dat
r8out = lal.CreateREAL8Vector(len(r8dat))
r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
assert(lal.swig_lal_test_viewin_REAL8Vector(r8out, r8))
assert((r8out.data == r8.data).all())
r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
assert(lal.swig_lal_test_viewin_REAL8Vector(r8out, r8dat))
assert((r8out.data == r8dat).all())
r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
assert(lal.swig_lal_test_viewinout_REAL8Vector(r8out, r8))
assert((2 * r8out.data == r8.data).all())
r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
assert(lal.swig_lal_test_viewinout_REAL8Vector(r8out, r8dat))
assert((2 * r8out.data == r8dat).all())
del r8
del r8out
del r8dat
lal.CheckMemoryLeaks()
c8 = lal.CreateCOMPLEX8Vector(len(c8dat))
c8.data = c8dat
c8out = lal.CreateCOMPLEX8Vector(len(c8dat))
c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
assert(lal.swig_lal_test_viewin_COMPLEX8Vector(c8out, c8))
assert((c8out.data == c8.data).all())
c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
assert(lal.swig_lal_test_viewin_COMPLEX8Vector(c8out, c8dat))
assert((c8out.data == c8dat).all())
c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
assert(lal.swig_lal_test_viewinout_COMPLEX8Vector(c8out, c8))
assert((2 * c8out.data == c8.data).all())
c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
assert(lal.swig_lal_test_viewinout_COMPLEX8Vector(c8out, c8dat))
assert((2 * c8out.data == c8dat).all())
del c8
del c8out
del c8dat
lal.CheckMemoryLeaks()
c16 = lal.CreateCOMPLEX16Vector(len(c16dat))
c16.data = c16dat
c16out = lal.CreateCOMPLEX16Vector(len(c16dat))
c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
assert(lal.swig_lal_test_viewin_COMPLEX16Vector(c16out, c16))
assert((c16out.data == c16.data).all())
c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
assert(lal.swig_lal_test_viewin_COMPLEX16Vector(c16out, c16dat))
assert((c16out.data == c16dat).all())
c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
assert(lal.swig_lal_test_viewinout_COMPLEX16Vector(c16out, c16))
assert((2 * c16out.data == c16.data).all())
c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
assert(lal.swig_lal_test_viewinout_COMPLEX16Vector(c16out, c16dat))
assert((2 * c16out.data == c16dat).all())
del c16
del c16out
del c16dat
lal.CheckMemoryLeaks()
r4dat = numpy.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]], dtype=numpy.float32)
r8dat = numpy.array([[3.4, 4.5], [5.6, 6.7], [7.8, 8.9]], dtype=numpy.float64)
c8dat = numpy.array(numpy.vectorize(complex)(r4dat, 8 + r4dat), dtype=numpy.complex64)
c16dat = numpy.array(numpy.vectorize(complex)(r8dat, 16 + r8dat), dtype=numpy.complex128)
r4 = lal.CreateREAL4VectorSequence(r4dat.shape[0], r4dat.shape[1])
r4.data = r4dat
r4out = lal.CreateREAL4VectorSequence(r4dat.shape[0], r4dat.shape[1])
r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
assert(lal.swig_lal_test_viewin_REAL4VectorSequence(r4out, r4))
assert((r4out.data == r4.data).all())
r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
assert(lal.swig_lal_test_viewin_REAL4VectorSequence(r4out, r4dat))
assert((r4out.data == r4dat).all())
r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
assert(lal.swig_lal_test_viewinout_REAL4VectorSequence(r4out, r4))
assert((2 * r4out.data == r4.data).all())
r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
assert(lal.swig_lal_test_viewinout_REAL4VectorSequence(r4out, r4dat))
assert((2 * r4out.data == r4dat).all())
del r4
del r4out
del r4dat
lal.CheckMemoryLeaks()
r8 = lal.CreateREAL8VectorSequence(r8dat.shape[0], r8dat.shape[1])
r8.data = r8dat
r8out = lal.CreateREAL8VectorSequence(r8dat.shape[0], r8dat.shape[1])
r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
assert(lal.swig_lal_test_viewin_REAL8VectorSequence(r8out, r8))
assert((r8out.data == r8.data).all())
r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
assert(lal.swig_lal_test_viewin_REAL8VectorSequence(r8out, r8dat))
assert((r8out.data == r8dat).all())
r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
assert(lal.swig_lal_test_viewinout_REAL8VectorSequence(r8out, r8))
assert((2 * r8out.data == r8.data).all())
r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
assert(lal.swig_lal_test_viewinout_REAL8VectorSequence(r8out, r8dat))
assert((2 * r8out.data == r8dat).all())
del r8
del r8out
del r8dat
lal.CheckMemoryLeaks()
c8 = lal.CreateCOMPLEX8VectorSequence(c8dat.shape[0], c8dat.shape[1])
c8.data = c8dat
c8out = lal.CreateCOMPLEX8VectorSequence(c8dat.shape[0], c8dat.shape[1])
c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
assert(lal.swig_lal_test_viewin_COMPLEX8VectorSequence(c8out, c8))
assert((c8out.data == c8.data).all())
c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
assert(lal.swig_lal_test_viewin_COMPLEX8VectorSequence(c8out, c8dat))
assert((c8out.data == c8dat).all())
c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
assert(lal.swig_lal_test_viewinout_COMPLEX8VectorSequence(c8out, c8))
assert((2 * c8out.data == c8.data).all())
c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
assert(lal.swig_lal_test_viewinout_COMPLEX8VectorSequence(c8out, c8dat))
assert((2 * c8out.data == c8dat).all())
del c8
del c8out
del c8dat
lal.CheckMemoryLeaks()
c16 = lal.CreateCOMPLEX16VectorSequence(c16dat.shape[0], c16dat.shape[1])
c16.data = c16dat
c16out = lal.CreateCOMPLEX16VectorSequence(c16dat.shape[0], c16dat.shape[1])
c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
assert(lal.swig_lal_test_viewin_COMPLEX16VectorSequence(c16out, c16))
assert((c16out.data == c16.data).all())
c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
assert(lal.swig_lal_test_viewin_COMPLEX16VectorSequence(c16out, c16dat))
assert((c16out.data == c16dat).all())
c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
assert(lal.swig_lal_test_viewinout_COMPLEX16VectorSequence(c16out, c16))
assert((2 * c16out.data == c16.data).all())
c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
assert(lal.swig_lal_test_viewinout_COMPLEX16VectorSequence(c16out, c16dat))
assert((2 * c16out.data == c16dat).all())
del c16
del c16out
del c16dat
lal.CheckMemoryLeaks()
print("PASSED input views of array structs (LAL)")
vfdat = numpy.array([1.2, 2.3, 3.4, 4.5, 5.6], dtype=numpy.float32)
vddat = numpy.array([3.4, 4.5, 5.6, 6.7, 7.8, 8.9], dtype=numpy.float64)
vcfdat = numpy.array(numpy.vectorize(complex)(vfdat, 8 + vfdat), dtype=numpy.complex64)
vcddat = numpy.array(numpy.vectorize(complex)(vddat, 16 + vddat), dtype=numpy.complex128)
vf = lal.gsl_vector_float(len(vfdat))
vf.data = vfdat
vfout = lal.gsl_vector_float(len(vfdat))
vfout.data = numpy.zeros(numpy.shape(vfdat), dtype=vfdat.dtype)
assert(lal.swig_lal_test_viewin_gsl_vector_float(vfout, vf))
assert((vfout.data == vf.data).all())
vfout.data = numpy.zeros(numpy.shape(vfdat), dtype=vfdat.dtype)
assert(lal.swig_lal_test_viewin_gsl_vector_float(vfout, vfdat))
assert((vfout.data == vfdat).all())
vfout.data = numpy.zeros(numpy.shape(vfdat), dtype=vfdat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_vector_float(vfout, vf))
assert((2 * vfout.data == vf.data).all())
vfout.data = numpy.zeros(numpy.shape(vfdat), dtype=vfdat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_vector_float(vfout, vfdat))
assert((2 * vfout.data == vfdat).all())
del vf
del vfout
del vfdat
lal.CheckMemoryLeaks()
vd = lal.gsl_vector(len(vddat))
vd.data = vddat
vdout = lal.gsl_vector(len(vddat))
vdout.data = numpy.zeros(numpy.shape(vddat), dtype=vddat.dtype)
assert(lal.swig_lal_test_viewin_gsl_vector(vdout, vd))
assert((vdout.data == vd.data).all())
vdout.data = numpy.zeros(numpy.shape(vddat), dtype=vddat.dtype)
assert(lal.swig_lal_test_viewin_gsl_vector(vdout, vddat))
assert((vdout.data == vddat).all())
vdout.data = numpy.zeros(numpy.shape(vddat), dtype=vddat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_vector(vdout, vd))
assert((2 * vdout.data == vd.data).all())
vdout.data = numpy.zeros(numpy.shape(vddat), dtype=vddat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_vector(vdout, vddat))
assert((2 * vdout.data == vddat).all())
del vd
del vdout
del vddat
lal.CheckMemoryLeaks()
vcf = lal.gsl_vector_complex_float(len(vcfdat))
vcf.data = vcfdat
vcfout = lal.gsl_vector_complex_float(len(vcfdat))
vcfout.data = numpy.zeros(numpy.shape(vcfdat), dtype=vcfdat.dtype)
assert(lal.swig_lal_test_viewin_gsl_vector_complex_float(vcfout, vcf))
assert((vcfout.data == vcf.data).all())
vcfout.data = numpy.zeros(numpy.shape(vcfdat), dtype=vcfdat.dtype)
assert(lal.swig_lal_test_viewin_gsl_vector_complex_float(vcfout, vcfdat))
assert((vcfout.data == vcfdat).all())
vcfout.data = numpy.zeros(numpy.shape(vcfdat), dtype=vcfdat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_vector_complex_float(vcfout, vcf))
assert((2 * vcfout.data == vcf.data).all())
vcfout.data = numpy.zeros(numpy.shape(vcfdat), dtype=vcfdat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_vector_complex_float(vcfout, vcfdat))
assert((2 * vcfout.data == vcfdat).all())
del vcf
del vcfout
del vcfdat
lal.CheckMemoryLeaks()
vcd = lal.gsl_vector_complex(len(vcddat))
vcd.data = vcddat
vcdout = lal.gsl_vector_complex(len(vcddat))
vcdout.data = numpy.zeros(numpy.shape(vcddat), dtype=vcddat.dtype)
assert(lal.swig_lal_test_viewin_gsl_vector_complex(vcdout, vcd))
assert((vcdout.data == vcd.data).all())
vcdout.data = numpy.zeros(numpy.shape(vcddat), dtype=vcddat.dtype)
assert(lal.swig_lal_test_viewin_gsl_vector_complex(vcdout, vcddat))
assert((vcdout.data == vcddat).all())
vcdout.data = numpy.zeros(numpy.shape(vcddat), dtype=vcddat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_vector_complex(vcdout, vcd))
assert((2 * vcdout.data == vcd.data).all())
vcdout.data = numpy.zeros(numpy.shape(vcddat), dtype=vcddat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_vector_complex(vcdout, vcddat))
assert((2 * vcdout.data == vcddat).all())
del vcd
del vcdout
del vcddat
lal.CheckMemoryLeaks()
mfdat = numpy.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]], dtype=numpy.float32)
mddat = numpy.array([[3.4, 4.5], [5.6, 6.7], [7.8, 8.9]], dtype=numpy.float64)
mcfdat = numpy.array(numpy.vectorize(complex)(mfdat, 8 + mfdat), dtype=numpy.complex64)
mcddat = numpy.array(numpy.vectorize(complex)(mddat, 16 + mddat), dtype=numpy.complex128)
mf = lal.gsl_matrix_float(mfdat.shape[0], mfdat.shape[1])
mf.data = mfdat
mfout = lal.gsl_matrix_float(mfdat.shape[0], mfdat.shape[1])
mfout.data = numpy.zeros(numpy.shape(mfdat), dtype=mfdat.dtype)
assert(lal.swig_lal_test_viewin_gsl_matrix_float(mfout, mf))
assert((mfout.data == mf.data).all())
mfout.data = numpy.zeros(numpy.shape(mfdat), dtype=mfdat.dtype)
assert(lal.swig_lal_test_viewin_gsl_matrix_float(mfout, mfdat))
assert((mfout.data == mfdat).all())
mfout.data = numpy.zeros(numpy.shape(mfdat), dtype=mfdat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_matrix_float(mfout, mf))
assert((2 * mfout.data == mf.data).all())
mfout.data = numpy.zeros(numpy.shape(mfdat), dtype=mfdat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_matrix_float(mfout, mfdat))
assert((2 * mfout.data == mfdat).all())
del mf
del mfout
del mfdat
lal.CheckMemoryLeaks()
md = lal.gsl_matrix(mddat.shape[0], mddat.shape[1])
md.data = mddat
mdout = lal.gsl_matrix(mddat.shape[0], mddat.shape[1])
mdout.data = numpy.zeros(numpy.shape(mddat), dtype=mddat.dtype)
assert(lal.swig_lal_test_viewin_gsl_matrix(mdout, md))
assert((mdout.data == md.data).all())
mdout.data = numpy.zeros(numpy.shape(mddat), dtype=mddat.dtype)
assert(lal.swig_lal_test_viewin_gsl_matrix(mdout, mddat))
assert((mdout.data == mddat).all())
mdout.data = numpy.zeros(numpy.shape(mddat), dtype=mddat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_matrix(mdout, md))
assert((2 * mdout.data == md.data).all())
mdout.data = numpy.zeros(numpy.shape(mddat), dtype=mddat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_matrix(mdout, mddat))
assert((2 * mdout.data == mddat).all())
del md
del mdout
del mddat
lal.CheckMemoryLeaks()
mcf = lal.gsl_matrix_complex_float(mcfdat.shape[0], mcfdat.shape[1])
mcf.data = mcfdat
mcfout = lal.gsl_matrix_complex_float(mcfdat.shape[0], mcfdat.shape[1])
mcfout.data = numpy.zeros(numpy.shape(mcfdat), dtype=mcfdat.dtype)
assert(lal.swig_lal_test_viewin_gsl_matrix_complex_float(mcfout, mcf))
assert((mcfout.data == mcf.data).all())
mcfout.data = numpy.zeros(numpy.shape(mcfdat), dtype=mcfdat.dtype)
assert(lal.swig_lal_test_viewin_gsl_matrix_complex_float(mcfout, mcfdat))
assert((mcfout.data == mcfdat).all())
mcfout.data = numpy.zeros(numpy.shape(mcfdat), dtype=mcfdat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_matrix_complex_float(mcfout, mcf))
assert((2 * mcfout.data == mcf.data).all())
mcfout.data = numpy.zeros(numpy.shape(mcfdat), dtype=mcfdat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_matrix_complex_float(mcfout, mcfdat))
assert((2 * mcfout.data == mcfdat).all())
del mcf
del mcfout
del mcfdat
lal.CheckMemoryLeaks()
mcd = lal.gsl_matrix_complex(mcddat.shape[0], mcddat.shape[1])
mcd.data = mcddat
mcdout = lal.gsl_matrix_complex(mcddat.shape[0], mcddat.shape[1])
mcdout.data = numpy.zeros(numpy.shape(mcddat), dtype=mcddat.dtype)
assert(lal.swig_lal_test_viewin_gsl_matrix_complex(mcdout, mcd))
assert((mcdout.data == mcd.data).all())
mcdout.data = numpy.zeros(numpy.shape(mcddat), dtype=mcddat.dtype)
assert(lal.swig_lal_test_viewin_gsl_matrix_complex(mcdout, mcddat))
assert((mcdout.data == mcddat).all())
mcdout.data = numpy.zeros(numpy.shape(mcddat), dtype=mcddat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_matrix_complex(mcdout, mcd))
assert((2 * mcdout.data == mcd.data).all())
mcdout.data = numpy.zeros(numpy.shape(mcddat), dtype=mcddat.dtype)
assert(lal.swig_lal_test_viewinout_gsl_matrix_complex(mcdout, mcddat))
assert((2 * mcdout.data == mcddat).all())
del mcd
del mcdout
del mcddat
lal.CheckMemoryLeaks()
print("PASSED input views of array structs (GSL)")
def check_input_view_type_safety(f, a, b, expect_exception):
    expected_exception = False
    if expect_exception:
        try:
            f(a, b)
            expected_exception = True
        except:
            pass
        assert(not expected_exception)
        try:
            f(b, a)
            expected_exception = True
        except:
            pass
        assert(not expected_exception)
    else:
        f(a, b)
        f(b, a)
r4 = numpy.zeros(10, dtype=numpy.float32)
r8 = numpy.zeros(10, dtype=numpy.float64)
c8 = numpy.zeros(10, dtype=numpy.complex64)
c16 = numpy.zeros(10, dtype=numpy.complex128)
check_input_view_type_safety(lal.swig_lal_test_viewinout_REAL4Vector, r4, r4, False)
check_input_view_type_safety(lal.swig_lal_test_viewinout_REAL4Vector, r4, r8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_REAL4Vector, r4, c8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_REAL4Vector, r4, c16, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_REAL8Vector, r8, r4, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_REAL8Vector, r8, r8, False)
check_input_view_type_safety(lal.swig_lal_test_viewinout_REAL8Vector, r8, c8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_REAL8Vector, r8, c16, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_COMPLEX8Vector, c8, r4, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_COMPLEX8Vector, c8, r8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_COMPLEX8Vector, c8, c8, False)
check_input_view_type_safety(lal.swig_lal_test_viewinout_COMPLEX8Vector, c8, c16, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_COMPLEX16Vector, c16, r4, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_COMPLEX16Vector, c16, r8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_COMPLEX16Vector, c16, c8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_COMPLEX16Vector, c16, c16, False)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_float, r4, r4, False)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_float, r4, r8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_float, r4, c8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_float, r4, c16, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector, r8, r4, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector, r8, r8, False)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector, r8, c8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector, r8, c16, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_complex_float, c8, r4, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_complex_float, c8, r8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_complex_float, c8, c8, False)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_complex_float, c8, c16, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_complex, c16, r4, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_complex, c16, r8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_complex, c16, c8, True)
check_input_view_type_safety(lal.swig_lal_test_viewinout_gsl_vector_complex, c16, c16, False)
del r4
del r8
del c8
del c16
lal.CheckMemoryLeaks()
print("PASSED input views of array structs (type safety)")

# check FFT functions with input views
print("check FFT functions with input views ...")
r4in = numpy.array(range(0, 32), dtype=numpy.float32)
r8in = numpy.array(range(0, 64), dtype=numpy.float64)
c8in = numpy.array(numpy.vectorize(complex)(8 + r4in, r4in), dtype=numpy.complex64)
c16in = numpy.array(numpy.vectorize(complex)(16 + r8in, r8in), dtype=numpy.complex128)
c8inv = lal.CreateCOMPLEX8Vector(len(c8in))
c8inv.data = c8in
c8outv = lal.CreateCOMPLEX8Vector(len(c8in))
plan = lal.CreateForwardCOMPLEX8FFTPlan(len(c8in), 0)
lal.COMPLEX8VectorFFT(c8outv, c8inv, plan)
c8out = numpy.zeros(numpy.shape(c8outv.data), dtype=c8outv.data.dtype)
lal.COMPLEX8VectorFFT(c8out, c8in, plan)
assert((c8out == c8outv.data).all())
del c8inv
del c8outv
del plan
lal.CheckMemoryLeaks()
c16inv = lal.CreateCOMPLEX16Vector(len(c16in))
c16inv.data = c16in
c16outv = lal.CreateCOMPLEX16Vector(len(c16in))
plan = lal.CreateForwardCOMPLEX16FFTPlan(len(c16in), 0)
lal.COMPLEX16VectorFFT(c16outv, c16inv, plan)
c16out = numpy.zeros(numpy.shape(c16outv.data), dtype=c16outv.data.dtype)
lal.COMPLEX16VectorFFT(c16out, c16in, plan)
assert((c16out == c16outv.data).all())
del c16inv
del c16outv
del plan
lal.CheckMemoryLeaks()
r4inv = lal.CreateREAL4Vector(len(r4in))
r4inv.data = r4in
c8outv = lal.CreateCOMPLEX8Vector(len(r4in)/2 + 1)
plan = lal.CreateForwardREAL4FFTPlan(len(r4in), 0)
lal.REAL4ForwardFFT(c8outv, r4inv, plan)
c8out = numpy.zeros(numpy.shape(c8outv.data), dtype=c8outv.data.dtype)
lal.REAL4ForwardFFT(c8out, r4in, plan)
assert((c8out == c8outv.data).all())
del r4inv
del c8outv
del plan
lal.CheckMemoryLeaks()
c8inv = lal.CreateCOMPLEX8Vector(len(c8in))
c8inv.data = c8in
r4outv = lal.CreateREAL4Vector((len(c8in)-1)*2)
plan = lal.CreateReverseREAL4FFTPlan((len(c8in)-1)*2, 0)
lal.REAL4ReverseFFT(r4outv, c8inv, plan)
r4out = numpy.zeros(numpy.shape(r4outv.data), dtype=r4outv.data.dtype)
lal.REAL4ReverseFFT(r4out, c8in, plan)
assert((r4out == r4outv.data).all())
del c8inv
del r4outv
del plan
lal.CheckMemoryLeaks()
r8inv = lal.CreateREAL8Vector(len(r8in))
r8inv.data = r8in
c16outv = lal.CreateCOMPLEX16Vector(len(r8in)/2 + 1)
plan = lal.CreateForwardREAL8FFTPlan(len(r8in), 0)
lal.REAL8ForwardFFT(c16outv, r8inv, plan)
c16out = numpy.zeros(numpy.shape(c16outv.data), dtype=c16outv.data.dtype)
lal.REAL8ForwardFFT(c16out, r8in, plan)
assert((c16out == c16outv.data).all())
del r8inv
del c16outv
del plan
lal.CheckMemoryLeaks()
c16inv = lal.CreateCOMPLEX16Vector(len(c16in))
c16inv.data = c16in
r8outv = lal.CreateREAL8Vector((len(c16in)-1)*2)
plan = lal.CreateReverseREAL8FFTPlan((len(c16in)-1)*2, 0)
lal.REAL8ReverseFFT(r8outv, c16inv, plan)
r8out = numpy.zeros(numpy.shape(r8outv.data), dtype=r8outv.data.dtype)
lal.REAL8ReverseFFT(r8out, c16in, plan)
assert((r8out == r8outv.data).all())
del c16inv
del r8outv
del plan
lal.CheckMemoryLeaks()
r4inv = lal.CreateREAL4Vector(len(r4in))
r4inv.data = r4in
r4outv = lal.CreateREAL4Vector(len(r4in))
plan = lal.CreateForwardREAL4FFTPlan(len(r4in), 0)
lal.REAL4VectorFFT(r4outv, r4inv, plan)
r4out = numpy.zeros(numpy.shape(r4outv.data), dtype=r4outv.data.dtype)
lal.REAL4VectorFFT(r4out, r4in, plan)
assert((r4out == r4outv.data).all())
del r4inv
del r4outv
del plan
lal.CheckMemoryLeaks()
r8inv = lal.CreateREAL8Vector(len(r8in))
r8inv.data = r8in
r8outv = lal.CreateREAL8Vector(len(r8in))
plan = lal.CreateForwardREAL8FFTPlan(len(r8in), 0)
lal.REAL8VectorFFT(r8outv, r8inv, plan)
r8out = numpy.zeros(numpy.shape(r8outv.data), dtype=r8outv.data.dtype)
lal.REAL8VectorFFT(r8out, r8in, plan)
assert((r8out == r8outv.data).all())
del r8inv
del r8outv
del plan
lal.CheckMemoryLeaks()
r4inv = lal.CreateREAL4Vector(len(r4in))
r4inv.data = r4in
r4outv = lal.CreateREAL4Vector(len(r4in)/2 + 1)
plan = lal.CreateForwardREAL4FFTPlan(len(r4in), 0)
lal.REAL4PowerSpectrum(r4outv, r4inv, plan)
r4out = numpy.zeros(numpy.shape(r4outv.data), dtype=r4outv.data.dtype)
lal.REAL4PowerSpectrum(r4out, r4in, plan)
assert((r4out == r4outv.data).all())
del r4inv
del r4outv
del plan
lal.CheckMemoryLeaks()
r8inv = lal.CreateREAL8Vector(len(r8in))
r8inv.data = r8in
r8outv = lal.CreateREAL8Vector(len(r8in)/2 + 1)
plan = lal.CreateForwardREAL8FFTPlan(len(r8in), 0)
lal.REAL8PowerSpectrum(r8outv, r8inv, plan)
r8out = numpy.zeros(numpy.shape(r8outv.data), dtype=r8outv.data.dtype)
lal.REAL8PowerSpectrum(r8out, r8in, plan)
assert((r8out == r8outv.data).all())
del r8inv
del r8outv
del plan
lal.CheckMemoryLeaks()
print("PASSED FFT functions with input views ...")

# check dynamic array of pointers access
print("checking dynamic array of pointers access ...")
ap = lal.swig_lal_test_Create_arrayofptrs(3)
assert(ap.length == 3)
for i in range(0, ap.length):
    assert(ap.data[i].length == 6)
    for j in range(0, ap.data[i].length):
        assert(ap.data[i].data[j] == 42*ap.length*i + j)
del ap
lal.CheckMemoryLeaks()
print("PASSED dynamic array of pointers access")

# check 'tm' struct conversions
print("checking 'tm' struct conversions ...")
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
print("PASSED 'tm' struct conversions")

# check LIGOTimeGPS operations
print("checking LIGOTimeGPS operations ...")
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
assert(3.5 + t1 == 14 and isinstance(3.5 + t1, LIGOTimeGPS))
t2 -= 5.5
assert(t2 == 5 and isinstance(t2, LIGOTimeGPS))
assert(t2 + 5.5 >= t1 and t2 + 3 != t2)
assert(t2 - 5 == t0 and isinstance(t2 - 5, LIGOTimeGPS))
assert(t1 * 3 == 31.5 and isinstance(t1 * 3, LIGOTimeGPS))
assert(3 * t1 == 31.5 and isinstance(3 * t1, LIGOTimeGPS))
assert(t2 / 2.5 == 2 and isinstance(t2 / 2.5, LIGOTimeGPS))
assert(21 / t1  == 2 and isinstance(21 / t1, LIGOTimeGPS))
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
assert(t1.ns() == 812345678250000000L)
assert(hash(t1) == 1049484238)
t4struct = lal.swig_lal_test_gps()
t4struct.t = 1234.5
assert(t4struct.t == 1234.5)
t5 = LIGOTimeGPS("1000")
assert(t5 == 1000)
try:
    t5 = LIGOTimeGPS("abc1000")
    expected_exception = True
except:
    pass
assert(not expected_exception)
try:
    t5 = LIGOTimeGPS("1000abc")
    expected_exception = True
except:
    pass
assert(not expected_exception)
assert(lal.swig_lal_test_noptrgps(LIGOTimeGPS(1234.5)) == lal.swig_lal_test_noptrgps(1234.5))
del t0
del t1
del t2
del t3
del t4struct
del t5
lal.CheckMemoryLeaks()
print("PASSED LIGOTimeGPS operations")

# check LALUnit operations
print("checking LALUnit operations ...")
u1 = lal.Unit("kg m s^-2")
assert(u1 == lal.NewtonUnit and isinstance(u1, lal.Unit))
assert(str(u1) == "m kg s^-2")
u2 = lal.MeterUnit * lal.KiloGramUnit / lal.SecondUnit ** 2
assert(u1 == u2 and isinstance(u2, lal.Unit));
u2 = lal.MeterUnit**(1,2) * lal.KiloGramUnit**(1,2) * lal.SecondUnit ** -1
assert(u1**(1,2) == u2 and isinstance(u2, lal.Unit));
try:
    lal.SecondUnit ** (1,0);
    expected_exception = True
except:
    pass
assert(not expected_exception)
u1 *= lal.MeterUnit
assert(u1 == lal.JouleUnit and isinstance(u1, lal.Unit))
assert(repr(u1) == "m^2 kg s^-2")
u1 /= lal.SecondUnit
assert(u1 == lal.WattUnit and isinstance(u1, lal.Unit))
assert(u1 == "m^2 kg s^-3")
u1 *= 1000
assert(u1 == lal.KiloUnit * lal.WattUnit)
assert(u1 == 1000 * lal.WattUnit)
assert(u1 == lal.WattUnit * 1000)
assert(u1 == lal.MegaUnit / 1000 * lal.WattUnit)
assert(int(u1) == 1000)
u1 /= 10000
assert(u1 == 100 * lal.MilliUnit * lal.WattUnit)
try:
    u1 *= 1.234
    expected_exception = True
except:
    pass
assert(not expected_exception)
assert(u1.norm() == u1)
del u1
del u2
lal.CheckMemoryLeaks()
print("PASSED LALUnit operations")

# passed all tests!
print("PASSED all tests")
