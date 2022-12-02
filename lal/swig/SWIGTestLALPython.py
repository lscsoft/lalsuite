# Check SWIG Python bindings for LAL
# Author: Karl Wette, 2011--2014

import os
import sys
import pytest

import warnings
import datetime
import pickle
import numpy

# return if 'x' has both value 'v' and type 't'
def is_value_and_type(x, v, t):
    return x == v and type(x) is t

# turn NumPy's ComplexWarning into an error, if available
if hasattr(numpy, "ComplexWarning"):
    warnings.simplefilter("error", numpy.ComplexWarning)

# check module load
print("checking module load ...", file=sys.stderr)
import lal
from lal import globalvar as lalglobalvar
lal_c_si = lal.C_SI
lal_180_pi = lal.LAL_180_PI
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

def test_memory_allocation():
    """check memory allocation
    """

    if not lal.MEMORY_FUNCTIONS_DISABLED:
        pytest.skip("LAL was built with MEMORY_FUNCTIONS_DISABLED")

    print("checking memory allocation ...", file=sys.stderr)
    lal.CheckMemoryLeaks()
    mem1 = lal.Detector()
    mem2 = lal.CreateCOMPLEX8Vector(5)
    mem3 = lal.CreateREAL8Vector(3)
    mem4 = lal.CreateREAL4TimeSeries("test", lal.LIGOTimeGPS(0), 100, 0.1, lal.DimensionlessUnit, 10)
    print("*** below should be an error message from CheckMemoryLeaks() ***", file=sys.stderr)
    set_nice_error_handlers()
    expected_exception = False
    try:
        lal.CheckMemoryLeaks()
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    print("*** above should be an error message from CheckMemoryLeaks() ***", file=sys.stderr)
    del mem1
    del mem2
    del mem3
    del mem4
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED memory allocation", file=sys.stderr)

def test_object_parent_tracking():
    """check object parent tracking
    """

    print("checking object parent tracking ...", file=sys.stderr)
    a = lal.gsl_vector(3)
    a.data = [1.1, 2.2, 3.3]
    b = a.data
    assert not b.flags['OWNDATA']
    assert (b == [1.1, 2.2, 3.3]).all()
    del a
    assert (b == [1.1, 2.2, 3.3]).all()
    ts = lal.CreateREAL8TimeSeries("test", lal.LIGOTimeGPS(0), 0, 0.1, lal.DimensionlessUnit, 10)
    ts.data.data = list(range(0, 10))
    for i in range(0, 7):
        v = ts.data
    assert (v.data == list(range(0, 10))).all()
    del ts
    assert (v.data == list(range(0, 10))).all()
    del v
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED object parent tracking", file=sys.stderr)

def test_equal_return_first_argument_type_handling():
    """check equal return/first argument type handling
    """

    print("checking equal return/first argument type handling", file=sys.stderr)
    sv = lal.CreateStringVector("1")
    assert sv.length == 1
    lal.AppendString2Vector(sv, "2")
    assert sv.length == 2
    sv = lal.AppendString2Vector(sv, "3")
    assert sv.length == 3
    sv2 = lal.AppendString2Vector(sv, "4")
    assert sv.length == 4
    assert sv2.length == 4
    assert sv == sv2
    del sv
    del sv2
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    ts = lal.CreateREAL8TimeSeries("ts", 800000000, 100, 0.1, lal.HertzUnit, 10)
    assert ts.data.length == 10
    lal.ResizeREAL8TimeSeries(ts, 0, 20)
    assert ts.data.length == 20
    ts = lal.ResizeREAL8TimeSeries(ts, 0, 30)
    assert ts.data.length == 30
    ts2 = lal.ResizeREAL8TimeSeries(ts, 0, 40)
    assert ts.data.length == 40
    assert ts2.data.length == 40
    assert ts == ts2
    del ts
    del ts2
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED equal return/first argument type handling", file=sys.stderr)

def test_string_conversions():
    """check string conversions
    """

    print("checking string conversions ...", file=sys.stderr)
    strs = ["a", "bc", "def"]
    sv = lal.CreateStringVector(*strs)
    assert sv.length == 3
    assert (sv.data.astype(numpy.object) == strs).all()
    strs[0] = "ghijk"
    sv.data[0] = strs[0]
    strs.append("lmnopq")
    sv = lal.AppendString2Vector(sv, strs[3])
    assert sv.length == 4
    for i in range(0, 4):
        assert sv.data[i] == strs[i]
    del sv
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED string conversions", file=sys.stderr)

def test_static_vector_matrix_conversions():
    """check static vector/matrix conversions
    """

    print("checking static vector/matrix conversions ...", file=sys.stderr)
    lalglobalvar.swig_lal_test_struct_vector[0] = lalglobalvar.swig_lal_test_struct_const
    assert lalglobalvar.swig_lal_test_struct_vector[0].n == lalglobalvar.swig_lal_test_struct_const.n
    assert lalglobalvar.swig_lal_test_struct_vector[0].i == lalglobalvar.swig_lal_test_struct_const.i
    assert lalglobalvar.swig_lal_test_struct_vector[0].f == lalglobalvar.swig_lal_test_struct_const.f
    assert lalglobalvar.swig_lal_test_struct_vector[0].str == lalglobalvar.swig_lal_test_struct_const.str
    assert (lalglobalvar.swig_lal_test_struct_vector[0].vec == lalglobalvar.swig_lal_test_struct_const.vec).all()
    lalglobalvar.swig_lal_test_struct_matrix[0, 0] = lalglobalvar.swig_lal_test_struct_const
    assert lalglobalvar.swig_lal_test_struct_matrix[0, 0].n == lalglobalvar.swig_lal_test_struct_const.n
    assert lalglobalvar.swig_lal_test_struct_matrix[0, 0].i == lalglobalvar.swig_lal_test_struct_const.i
    assert lalglobalvar.swig_lal_test_struct_matrix[0, 0].f == lalglobalvar.swig_lal_test_struct_const.f
    assert lalglobalvar.swig_lal_test_struct_matrix[0, 0].str == lalglobalvar.swig_lal_test_struct_const.str
    assert (lalglobalvar.swig_lal_test_struct_matrix[0, 0].vec == lalglobalvar.swig_lal_test_struct_const.vec).all()
    sts = lal.swig_lal_test_struct()
    assert len(sts.vec) == 3
    assert len(sts.evec) == 3
    assert sts.mat.shape == (2, 3)
    sts.vec = [3, 2, 1]
    assert (sts.vec == [3, 2, 1]).all()
    sts.mat = [[4, 5, 6], (9, 8, 7)]
    set_nice_error_handlers()
    expected_exception = False
    try:
        sts.mat = [[1.1, 2.3, 4.5], [6.5, 4.3, 2.1]]
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    assert (sts.mat == [[4, 5, 6], [9, 8, 7]]).all()
    for i in range(0, 3):
        sts.evec[i] = 2*i + 3
        assert sts.evec[i] == (2*i + 3)
    del sts
    assert not lalglobalvar.swig_lal_test_enum_vector.any()
    assert not lalglobalvar.swig_lal_test_enum_matrix.any()
    assert len(lalglobalvar.swig_lal_test_empty_INT4_vector) == 0
    assert not lalglobalvar.swig_lal_test_INT4_vector.any()
    assert not lalglobalvar.swig_lal_test_INT4_matrix.any()
    assert not lalglobalvar.swig_lal_test_REAL8_vector.any()
    assert not lalglobalvar.swig_lal_test_REAL8_matrix.any()
    assert not lalglobalvar.swig_lal_test_COMPLEX8_vector.any()
    assert not lalglobalvar.swig_lal_test_COMPLEX8_matrix.any()
    lalglobalvar.swig_lal_test_INT4_vector[0] = 10
    assert lalglobalvar.swig_lal_test_INT4_vector[0] == 10
    lalglobalvar.swig_lal_test_INT4_matrix[0, 0] = 11
    assert lalglobalvar.swig_lal_test_INT4_matrix[0, 0] == 11
    lalglobalvar.swig_lal_test_INT4_vector = lalglobalvar.swig_lal_test_INT4_const_vector
    assert (lalglobalvar.swig_lal_test_INT4_vector == [1, 2, 4]).all()
    assert lalglobalvar.swig_lal_test_INT4_const_vector[2] == 4
    lalglobalvar.swig_lal_test_INT4_matrix = lalglobalvar.swig_lal_test_INT4_const_matrix
    assert (lalglobalvar.swig_lal_test_INT4_matrix == [[1, 2, 4], [2, 4, 8]]).all()
    assert lalglobalvar.swig_lal_test_INT4_const_matrix[1, 2] == 8
    set_nice_error_handlers()
    expected_exception = False
    try:
        lalglobalvar.swig_lal_test_INT4_const_vector(20)
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    lalglobalvar.swig_lal_test_REAL8_vector[0] = 3.4
    assert lalglobalvar.swig_lal_test_REAL8_vector[0] == 3.4
    lalglobalvar.swig_lal_test_REAL8_matrix[0, 0] = 5.6
    assert lalglobalvar.swig_lal_test_REAL8_matrix[0, 0] == 5.6
    lalglobalvar.swig_lal_test_COMPLEX8_vector[0] = complex(3.5, 4.75)
    assert lalglobalvar.swig_lal_test_COMPLEX8_vector[0] == complex(3.5, 4.75)
    lalglobalvar.swig_lal_test_COMPLEX8_matrix[0, 0] = complex(5.5, 6.25)
    assert lalglobalvar.swig_lal_test_COMPLEX8_matrix[0, 0] == complex(5.5, 6.25)
    print("PASSED static vector/matrix conversions", file=sys.stderr)

def test_dynamic_vector_matrix_conversions():
    """check dynamic vector/matrix conversions
    """

    print("checking dynamic vector/matrix conversions ...", file=sys.stderr)
    def check_dynamic_vector_matrix(iv, ivl, rv, rvl, cm, cms1, cms2):
        expected_exception = False
        iv.data = numpy.zeros(ivl, dtype=iv.data.dtype)
        rv.data = numpy.zeros(rvl, dtype=rv.data.dtype)
        cm.data = numpy.zeros((cms1, cms2), dtype=cm.data.dtype)
        assert ivl == 5
        iv.data = [1, 3, 2, 4, 3]
        assert (iv.data == [1, 3, 2, 4, 3]).all()
        iv.data[3] = 7
        assert iv.data[3] == 7
        assert rvl == 5
        rv.data = [1.2, 3.4, 2.6, 4.8, 3.5]
        assert (rv.data == [1.2, 3.4, 2.6, 4.8, 3.5]).all()
        rv.data[rvl - 1] = 7.5
        assert rv.data[rvl - 1] == 7.5
        set_nice_error_handlers()
        expected_exception = False
        try:
            rv.data[rvl] = 99.9
            expected_exception = True
        except:
            pass
        set_default_error_handlers()
        assert not expected_exception
        set_nice_error_handlers()
        expected_exception = False
        try:
            iv.data = rv.data
            expected_exception = True
        except:
            pass
        set_default_error_handlers()
        assert not expected_exception
        rv.data = iv.data
        assert (rv.data == iv.data).all()
        assert cms1 == 4
        assert cms2 == 6
        for i in range(0, cms1):
            for j in range(0, cms2):
                cm.data[i, j] = complex(i / 4.0, j / 2.0)
        assert cm.data[2, 3] == complex(0.5, 1.5)
        assert cm.data[3, 2] == complex(0.75, 1.0)
        set_nice_error_handlers()
        expected_exception = False
        try:
            iv.data[0] = cm.data[2, 3]
            if not hasattr(numpy, "ComplexWarning"):
                raise Exception("NumPy %s does not have ComplexWarning" % numpy.__version__)
            expected_exception = True
        except:
            pass
        set_default_error_handlers()
        assert not expected_exception
        set_nice_error_handlers()
        expected_exception = False
        try:
            rv.data[0] = cm.data[3, 2]
            if not hasattr(numpy, "ComplexWarning"):
                raise Exception("NumPy %s does not have ComplexWarning" % numpy.__version__)
            expected_exception = True
        except:
            pass
        set_default_error_handlers()
        assert not expected_exception
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
    assert rv0.length == 0
    assert rv0.data is None
    del rv0
    rv1 = lal.CreateREAL8Vector(1)
    rv1.data[0] = 1
    del rv1
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED dynamic vector/matrix conversions (LAL)", file=sys.stderr)
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
    print("PASSED dynamic vector/matrix conversions (GSL)", file=sys.stderr)

def test_fixed_and_dynamic_arrays_typemaps():
    """check fixed and dynamic arrays typemaps
    """

    print("checking fixed and dynamic arrays typemaps ...", file=sys.stderr)
    a1in = numpy.array([1.2, 3.5, 7.9], dtype=numpy.double)
    a1out = a1in * 2.5
    assert (lal.swig_lal_test_copyin_array1(a1in, 2.5) == a1out).all()
    a2in = numpy.array([[3,2], [7,6], [12,10]], dtype=numpy.int32)
    a2out = a2in * 15
    assert (lal.swig_lal_test_copyin_array2(a2in, 15) == a2out).all()
    a3in = numpy.array([lal.LIGOTimeGPS(1234.5), lal.LIGOTimeGPS(678.9)])
    a3out = a3in * 3
    assert (lal.swig_lal_test_copyin_array3(a3in, 3) == a3out).all()
    set_nice_error_handlers()
    expected_exception = False
    try:
        lal.swig_lal_test_copyin_array1(numpy.array([0,0,0,0], dtype=numpy.double), 0)
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    set_nice_error_handlers()
    expected_exception = False
    try:
        lal.swig_lal_test_copyin_array2(numpy.array([[1.2,3.4],[0,0],[0,0]], dtype=numpy.double), 0)
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    del a3in
    del a3out
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED fixed and dynamic arrays typemaps", file=sys.stderr)

def test_input_views_of_string_array_structs():
    """check input views of string array structs
    """

    print("checking input views of string array structs ...", file=sys.stderr)
    svdat = ["a", "bc", "def"]
    sv = lal.CreateEmptyStringVector(len(svdat))
    sv.data = svdat;
    svout = lal.CreateEmptyStringVector(len(svdat));
    svout.data = [""] * len(svdat)
    assert lal.swig_lal_test_viewin_LALStringVector(svout, sv)
    assert all(map(lambda x, y: x == y, svout.data, sv.data))
    svout.data = [""] * len(svdat)
    assert lal.swig_lal_test_viewin_LALStringVector(svout, svdat)
    assert all(map(lambda x, y: x == y, svout.data, svdat))
    sv.data = svdat
    assert lal.swig_lal_test_copyinout_LALStringVector(sv)
    assert all(map(lambda x, y: x == y.upper(), sv.data, svdat))
    sv.data = svdat
    retn, sv = lal.swig_lal_test_copyinout_LALStringVector(sv)
    assert retn
    assert all(map(lambda x, y: x == y.upper(), sv.data, svdat))
    sv = svdat
    retn, sv = lal.swig_lal_test_copyinout_LALStringVector(sv)
    assert retn
    assert all(map(lambda x, y: x == y.upper(), sv, svdat))
    del sv
    del svout
    del svdat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED input views of string array structs", file=sys.stderr)

def test_input_views_of_numeric_array_structs():
    """check input views of numeric array structs
    """

    print("checking input views of numeric array structs ...", file=sys.stderr)
    r4dat = numpy.array([1.2, 2.3, 3.4, 4.5, 5.6], dtype=numpy.float32)
    r8dat = numpy.array([3.4, 4.5, 5.6, 6.7, 7.8, 8.9], dtype=numpy.float64)
    c8dat = numpy.array(numpy.vectorize(complex)(r4dat, 8 + r4dat), dtype=numpy.complex64)
    c16dat = numpy.array(numpy.vectorize(complex)(r8dat, 16 + r8dat), dtype=numpy.complex128)
    r4 = lal.CreateREAL4Vector(len(r4dat))
    r4.data = r4dat
    r4out = lal.CreateREAL4Vector(len(r4dat))
    r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
    assert lal.swig_lal_test_viewin_REAL4Vector(r4out, r4)
    assert (r4out.data == r4.data).all()
    r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
    assert lal.swig_lal_test_viewin_REAL4Vector(r4out, r4dat)
    assert (r4out.data == r4dat).all()
    r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
    assert lal.swig_lal_test_viewinout_REAL4Vector(r4out, r4)
    assert (2 * r4out.data == r4.data).all()
    r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
    assert lal.swig_lal_test_viewinout_REAL4Vector(r4out, r4dat)
    assert (2 * r4out.data == r4dat).all()
    r4.data = r4dat
    assert lal.swig_lal_test_copyinout_REAL4Vector(r4)
    assert (r4.data == 3 * r4dat).all()
    r4.data = r4dat
    retn, r4 = lal.swig_lal_test_copyinout_REAL4Vector(r4)
    assert retn
    assert (r4.data == 3 * r4dat).all()
    r4 = r4dat
    retn, r4 = lal.swig_lal_test_copyinout_REAL4Vector(r4)
    assert retn
    assert (r4 == 3 * r4dat).all()
    del r4
    del r4out
    del r4dat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    r8 = lal.CreateREAL8Vector(len(r8dat))
    r8.data = r8dat
    r8out = lal.CreateREAL8Vector(len(r8dat))
    r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
    assert lal.swig_lal_test_viewin_REAL8Vector(r8out, r8)
    assert (r8out.data == r8.data).all()
    r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
    assert lal.swig_lal_test_viewin_REAL8Vector(r8out, r8dat)
    assert (r8out.data == r8dat).all()
    r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
    assert lal.swig_lal_test_viewinout_REAL8Vector(r8out, r8)
    assert (2 * r8out.data == r8.data).all()
    r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
    assert lal.swig_lal_test_viewinout_REAL8Vector(r8out, r8dat)
    assert (2 * r8out.data == r8dat).all()
    r8.data = r8dat
    assert lal.swig_lal_test_copyinout_REAL8Vector(r8)
    assert (r8.data == 3 * r8dat).all()
    r8.data = r8dat
    retn, r8 = lal.swig_lal_test_copyinout_REAL8Vector(r8)
    assert retn
    assert (r8.data == 3 * r8dat).all()
    r8 = r8dat
    retn, r8 = lal.swig_lal_test_copyinout_REAL8Vector(r8)
    assert retn
    assert (r8 == 3 * r8dat).all()
    del r8
    del r8out
    del r8dat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    c8 = lal.CreateCOMPLEX8Vector(len(c8dat))
    c8.data = c8dat
    c8out = lal.CreateCOMPLEX8Vector(len(c8dat))
    c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
    assert lal.swig_lal_test_viewin_COMPLEX8Vector(c8out, c8)
    assert (c8out.data == c8.data).all()
    c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
    assert lal.swig_lal_test_viewin_COMPLEX8Vector(c8out, c8dat)
    assert (c8out.data == c8dat).all()
    c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
    assert lal.swig_lal_test_viewinout_COMPLEX8Vector(c8out, c8)
    assert (2 * c8out.data == c8.data).all()
    c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
    assert lal.swig_lal_test_viewinout_COMPLEX8Vector(c8out, c8dat)
    assert (2 * c8out.data == c8dat).all()
    c8.data = c8dat
    assert lal.swig_lal_test_copyinout_COMPLEX8Vector(c8)
    assert (c8.data == 3 * c8dat).all()
    c8.data = c8dat
    retn, c8 = lal.swig_lal_test_copyinout_COMPLEX8Vector(c8)
    assert retn
    assert (c8.data == 3 * c8dat).all()
    c8 = c8dat
    retn, c8 = lal.swig_lal_test_copyinout_COMPLEX8Vector(c8)
    assert retn
    assert (c8 == 3 * c8dat).all()
    del c8
    del c8out
    del c8dat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    c16 = lal.CreateCOMPLEX16Vector(len(c16dat))
    c16.data = c16dat
    c16out = lal.CreateCOMPLEX16Vector(len(c16dat))
    c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
    assert lal.swig_lal_test_viewin_COMPLEX16Vector(c16out, c16)
    assert (c16out.data == c16.data).all()
    c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
    assert lal.swig_lal_test_viewin_COMPLEX16Vector(c16out, c16dat)
    assert (c16out.data == c16dat).all()
    c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
    assert lal.swig_lal_test_viewinout_COMPLEX16Vector(c16out, c16)
    assert (2 * c16out.data == c16.data).all()
    c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
    assert lal.swig_lal_test_viewinout_COMPLEX16Vector(c16out, c16dat)
    assert (2 * c16out.data == c16dat).all()
    c16.data = c16dat
    assert lal.swig_lal_test_copyinout_COMPLEX16Vector(c16)
    assert (c16.data == 3 * c16dat).all()
    c16.data = c16dat
    retn, c16 = lal.swig_lal_test_copyinout_COMPLEX16Vector(c16)
    assert retn
    assert (c16.data == 3 * c16dat).all()
    c16 = c16dat
    retn, c16 = lal.swig_lal_test_copyinout_COMPLEX16Vector(c16)
    assert retn
    assert (c16 == 3 * c16dat).all()
    del c16
    del c16out
    del c16dat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    r4dat = numpy.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]], dtype=numpy.float32)
    r8dat = numpy.array([[3.4, 4.5], [5.6, 6.7], [7.8, 8.9]], dtype=numpy.float64)
    c8dat = numpy.array(numpy.vectorize(complex)(r4dat, 8 + r4dat), dtype=numpy.complex64)
    c16dat = numpy.array(numpy.vectorize(complex)(r8dat, 16 + r8dat), dtype=numpy.complex128)
    r4 = lal.CreateREAL4VectorSequence(r4dat.shape[0], r4dat.shape[1])
    r4.data = r4dat
    r4out = lal.CreateREAL4VectorSequence(r4dat.shape[0], r4dat.shape[1])
    r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
    assert lal.swig_lal_test_viewin_REAL4VectorSequence(r4out, r4)
    assert (r4out.data == r4.data).all()
    r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
    assert lal.swig_lal_test_viewin_REAL4VectorSequence(r4out, r4dat)
    assert (r4out.data == r4dat).all()
    r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
    assert lal.swig_lal_test_viewinout_REAL4VectorSequence(r4out, r4)
    assert (2 * r4out.data == r4.data).all()
    r4out.data = numpy.zeros(numpy.shape(r4dat), dtype=r4dat.dtype)
    assert lal.swig_lal_test_viewinout_REAL4VectorSequence(r4out, r4dat)
    assert (2 * r4out.data == r4dat).all()
    r4.data = r4dat
    assert lal.swig_lal_test_copyinout_REAL4VectorSequence(r4)
    assert (r4.data == 3 * r4dat).all()
    r4.data = r4dat
    retn, r4 = lal.swig_lal_test_copyinout_REAL4VectorSequence(r4)
    assert retn
    assert (r4.data == 3 * r4dat).all()
    r4 = r4dat
    retn, r4 = lal.swig_lal_test_copyinout_REAL4VectorSequence(r4)
    assert retn
    assert (r4 == 3 * r4dat).all()
    del r4
    del r4out
    del r4dat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    r8 = lal.CreateREAL8VectorSequence(r8dat.shape[0], r8dat.shape[1])
    r8.data = r8dat
    r8out = lal.CreateREAL8VectorSequence(r8dat.shape[0], r8dat.shape[1])
    r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
    assert lal.swig_lal_test_viewin_REAL8VectorSequence(r8out, r8)
    assert (r8out.data == r8.data).all()
    r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
    assert lal.swig_lal_test_viewin_REAL8VectorSequence(r8out, r8dat)
    assert (r8out.data == r8dat).all()
    r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
    assert lal.swig_lal_test_viewinout_REAL8VectorSequence(r8out, r8)
    assert (2 * r8out.data == r8.data).all()
    r8out.data = numpy.zeros(numpy.shape(r8dat), dtype=r8dat.dtype)
    assert lal.swig_lal_test_viewinout_REAL8VectorSequence(r8out, r8dat)
    assert (2 * r8out.data == r8dat).all()
    r8.data = r8dat
    assert lal.swig_lal_test_copyinout_REAL8VectorSequence(r8)
    assert (r8.data == 3 * r8dat).all()
    r8.data = r8dat
    retn, r8 = lal.swig_lal_test_copyinout_REAL8VectorSequence(r8)
    assert retn
    assert (r8.data == 3 * r8dat).all()
    r8 = r8dat
    retn, r8 = lal.swig_lal_test_copyinout_REAL8VectorSequence(r8)
    assert retn
    assert (r8 == 3 * r8dat).all()
    del r8
    del r8out
    del r8dat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    c8 = lal.CreateCOMPLEX8VectorSequence(c8dat.shape[0], c8dat.shape[1])
    c8.data = c8dat
    c8out = lal.CreateCOMPLEX8VectorSequence(c8dat.shape[0], c8dat.shape[1])
    c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
    assert lal.swig_lal_test_viewin_COMPLEX8VectorSequence(c8out, c8)
    assert (c8out.data == c8.data).all()
    c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
    assert lal.swig_lal_test_viewin_COMPLEX8VectorSequence(c8out, c8dat)
    assert (c8out.data == c8dat).all()
    c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
    assert lal.swig_lal_test_viewinout_COMPLEX8VectorSequence(c8out, c8)
    assert (2 * c8out.data == c8.data).all()
    c8out.data = numpy.zeros(numpy.shape(c8dat), dtype=c8dat.dtype)
    assert lal.swig_lal_test_viewinout_COMPLEX8VectorSequence(c8out, c8dat)
    assert (2 * c8out.data == c8dat).all()
    c8.data = c8dat
    assert lal.swig_lal_test_copyinout_COMPLEX8VectorSequence(c8)
    assert (c8.data == 3 * c8dat).all()
    c8.data = c8dat
    retn, c8 = lal.swig_lal_test_copyinout_COMPLEX8VectorSequence(c8)
    assert retn
    assert (c8.data == 3 * c8dat).all()
    c8 = c8dat
    retn, c8 = lal.swig_lal_test_copyinout_COMPLEX8VectorSequence(c8)
    assert retn
    assert (c8 == 3 * c8dat).all()
    del c8
    del c8out
    del c8dat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    c16 = lal.CreateCOMPLEX16VectorSequence(c16dat.shape[0], c16dat.shape[1])
    c16.data = c16dat
    c16out = lal.CreateCOMPLEX16VectorSequence(c16dat.shape[0], c16dat.shape[1])
    c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
    assert lal.swig_lal_test_viewin_COMPLEX16VectorSequence(c16out, c16)
    assert (c16out.data == c16.data).all()
    c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
    assert lal.swig_lal_test_viewin_COMPLEX16VectorSequence(c16out, c16dat)
    assert (c16out.data == c16dat).all()
    c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
    assert lal.swig_lal_test_viewinout_COMPLEX16VectorSequence(c16out, c16)
    assert (2 * c16out.data == c16.data).all()
    c16out.data = numpy.zeros(numpy.shape(c16dat), dtype=c16dat.dtype)
    assert lal.swig_lal_test_viewinout_COMPLEX16VectorSequence(c16out, c16dat)
    assert (2 * c16out.data == c16dat).all()
    c16.data = c16dat
    assert lal.swig_lal_test_copyinout_COMPLEX16VectorSequence(c16)
    assert (c16.data == 3 * c16dat).all()
    c16.data = c16dat
    retn, c16 = lal.swig_lal_test_copyinout_COMPLEX16VectorSequence(c16)
    assert retn
    assert (c16.data == 3 * c16dat).all()
    c16 = c16dat
    retn, c16 = lal.swig_lal_test_copyinout_COMPLEX16VectorSequence(c16)
    assert retn
    assert (c16 == 3 * c16dat).all()
    del c16
    del c16out
    del c16dat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED input views of numeric array structs (LAL)", file=sys.stderr)
    vfdat = numpy.array([1.2, 2.3, 3.4, 4.5, 5.6], dtype=numpy.float32)
    vddat = numpy.array([3.4, 4.5, 5.6, 6.7, 7.8, 8.9], dtype=numpy.float64)
    vcfdat = numpy.array(numpy.vectorize(complex)(vfdat, 8 + vfdat), dtype=numpy.complex64)
    vcddat = numpy.array(numpy.vectorize(complex)(vddat, 16 + vddat), dtype=numpy.complex128)
    vf = lal.gsl_vector_float(len(vfdat))
    vf.data = vfdat
    vfout = lal.gsl_vector_float(len(vfdat))
    vfout.data = numpy.zeros(numpy.shape(vfdat), dtype=vfdat.dtype)
    assert lal.swig_lal_test_viewin_gsl_vector_float(vfout, vf)
    assert (vfout.data == vf.data).all()
    vfout.data = numpy.zeros(numpy.shape(vfdat), dtype=vfdat.dtype)
    assert lal.swig_lal_test_viewin_gsl_vector_float(vfout, vfdat)
    assert (vfout.data == vfdat).all()
    vfout.data = numpy.zeros(numpy.shape(vfdat), dtype=vfdat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_vector_float(vfout, vf)
    assert (2 * vfout.data == vf.data).all()
    vfout.data = numpy.zeros(numpy.shape(vfdat), dtype=vfdat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_vector_float(vfout, vfdat)
    assert (2 * vfout.data == vfdat).all()
    vf.data = vfdat
    assert lal.swig_lal_test_copyinout_gsl_vector_float(vf)
    assert (vf.data == 3 * vfdat).all()
    vf.data = vfdat
    retn, vf = lal.swig_lal_test_copyinout_gsl_vector_float(vf)
    assert retn
    assert (vf.data == 3 * vfdat).all()
    vf = vfdat
    retn, vf = lal.swig_lal_test_copyinout_gsl_vector_float(vf)
    assert retn
    assert (vf == 3 * vfdat).all()
    del vf
    del vfout
    del vfdat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    vd = lal.gsl_vector(len(vddat))
    vd.data = vddat
    vdout = lal.gsl_vector(len(vddat))
    vdout.data = numpy.zeros(numpy.shape(vddat), dtype=vddat.dtype)
    assert lal.swig_lal_test_viewin_gsl_vector(vdout, vd)
    assert (vdout.data == vd.data).all()
    vdout.data = numpy.zeros(numpy.shape(vddat), dtype=vddat.dtype)
    assert lal.swig_lal_test_viewin_gsl_vector(vdout, vddat)
    assert (vdout.data == vddat).all()
    vdout.data = numpy.zeros(numpy.shape(vddat), dtype=vddat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_vector(vdout, vd)
    assert (2 * vdout.data == vd.data).all()
    vdout.data = numpy.zeros(numpy.shape(vddat), dtype=vddat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_vector(vdout, vddat)
    assert (2 * vdout.data == vddat).all()
    vd.data = vddat
    assert lal.swig_lal_test_copyinout_gsl_vector(vd)
    assert (vd.data == 3 * vddat).all()
    vd.data = vddat
    retn, vd = lal.swig_lal_test_copyinout_gsl_vector(vd)
    assert retn
    assert (vd.data == 3 * vddat).all()
    vd = vddat
    retn, vd = lal.swig_lal_test_copyinout_gsl_vector(vd)
    assert retn
    assert (vd == 3 * vddat).all()
    del vd
    del vdout
    del vddat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    vcf = lal.gsl_vector_complex_float(len(vcfdat))
    vcf.data = vcfdat
    vcfout = lal.gsl_vector_complex_float(len(vcfdat))
    vcfout.data = numpy.zeros(numpy.shape(vcfdat), dtype=vcfdat.dtype)
    assert lal.swig_lal_test_viewin_gsl_vector_complex_float(vcfout, vcf)
    assert (vcfout.data == vcf.data).all()
    vcfout.data = numpy.zeros(numpy.shape(vcfdat), dtype=vcfdat.dtype)
    assert lal.swig_lal_test_viewin_gsl_vector_complex_float(vcfout, vcfdat)
    assert (vcfout.data == vcfdat).all()
    vcfout.data = numpy.zeros(numpy.shape(vcfdat), dtype=vcfdat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_vector_complex_float(vcfout, vcf)
    assert (2 * vcfout.data == vcf.data).all()
    vcfout.data = numpy.zeros(numpy.shape(vcfdat), dtype=vcfdat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_vector_complex_float(vcfout, vcfdat)
    assert (2 * vcfout.data == vcfdat).all()
    vcf.data = vcfdat
    assert lal.swig_lal_test_copyinout_gsl_vector_complex_float(vcf)
    assert (vcf.data == 3 * vcfdat).all()
    vcf.data = vcfdat
    retn, vcf = lal.swig_lal_test_copyinout_gsl_vector_complex_float(vcf)
    assert retn
    assert (vcf.data == 3 * vcfdat).all()
    vcf = vcfdat
    retn, vcf = lal.swig_lal_test_copyinout_gsl_vector_complex_float(vcf)
    assert retn
    assert (vcf == 3 * vcfdat).all()
    del vcf
    del vcfout
    del vcfdat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    vcd = lal.gsl_vector_complex(len(vcddat))
    vcd.data = vcddat
    vcdout = lal.gsl_vector_complex(len(vcddat))
    vcdout.data = numpy.zeros(numpy.shape(vcddat), dtype=vcddat.dtype)
    assert lal.swig_lal_test_viewin_gsl_vector_complex(vcdout, vcd)
    assert (vcdout.data == vcd.data).all()
    vcdout.data = numpy.zeros(numpy.shape(vcddat), dtype=vcddat.dtype)
    assert lal.swig_lal_test_viewin_gsl_vector_complex(vcdout, vcddat)
    assert (vcdout.data == vcddat).all()
    vcdout.data = numpy.zeros(numpy.shape(vcddat), dtype=vcddat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_vector_complex(vcdout, vcd)
    assert (2 * vcdout.data == vcd.data).all()
    vcdout.data = numpy.zeros(numpy.shape(vcddat), dtype=vcddat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_vector_complex(vcdout, vcddat)
    assert (2 * vcdout.data == vcddat).all()
    vcd.data = vcddat
    assert lal.swig_lal_test_copyinout_gsl_vector_complex(vcd)
    assert (vcd.data == 3 * vcddat).all()
    vcd.data = vcddat
    retn, vcd = lal.swig_lal_test_copyinout_gsl_vector_complex(vcd)
    assert retn
    assert (vcd.data == 3 * vcddat).all()
    vcd = vcddat
    retn, vcd = lal.swig_lal_test_copyinout_gsl_vector_complex(vcd)
    assert retn
    assert (vcd == 3 * vcddat).all()
    del vcd
    del vcdout
    del vcddat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    mfdat = numpy.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]], dtype=numpy.float32)
    mddat = numpy.array([[3.4, 4.5], [5.6, 6.7], [7.8, 8.9]], dtype=numpy.float64)
    mcfdat = numpy.array(numpy.vectorize(complex)(mfdat, 8 + mfdat), dtype=numpy.complex64)
    mcddat = numpy.array(numpy.vectorize(complex)(mddat, 16 + mddat), dtype=numpy.complex128)
    mf = lal.gsl_matrix_float(mfdat.shape[0], mfdat.shape[1])
    mf.data = mfdat
    mfout = lal.gsl_matrix_float(mfdat.shape[0], mfdat.shape[1])
    mfout.data = numpy.zeros(numpy.shape(mfdat), dtype=mfdat.dtype)
    assert lal.swig_lal_test_viewin_gsl_matrix_float(mfout, mf)
    assert (mfout.data == mf.data).all()
    mfout.data = numpy.zeros(numpy.shape(mfdat), dtype=mfdat.dtype)
    assert lal.swig_lal_test_viewin_gsl_matrix_float(mfout, mfdat)
    assert (mfout.data == mfdat).all()
    mfout.data = numpy.zeros(numpy.shape(mfdat), dtype=mfdat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_matrix_float(mfout, mf)
    assert (2 * mfout.data == mf.data).all()
    mfout.data = numpy.zeros(numpy.shape(mfdat), dtype=mfdat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_matrix_float(mfout, mfdat)
    assert (2 * mfout.data == mfdat).all()
    mf.data = mfdat
    assert lal.swig_lal_test_copyinout_gsl_matrix_float(mf)
    assert (mf.data == 3 * mfdat).all()
    mf.data = mfdat
    retn, mf = lal.swig_lal_test_copyinout_gsl_matrix_float(mf)
    assert retn
    assert (mf.data == 3 * mfdat).all()
    mf = mfdat
    retn, mf = lal.swig_lal_test_copyinout_gsl_matrix_float(mf)
    assert retn
    assert (mf == 3 * mfdat).all()
    del mf
    del mfout
    del mfdat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    md = lal.gsl_matrix(mddat.shape[0], mddat.shape[1])
    md.data = mddat
    mdout = lal.gsl_matrix(mddat.shape[0], mddat.shape[1])
    mdout.data = numpy.zeros(numpy.shape(mddat), dtype=mddat.dtype)
    assert lal.swig_lal_test_viewin_gsl_matrix(mdout, md)
    assert (mdout.data == md.data).all()
    mdout.data = numpy.zeros(numpy.shape(mddat), dtype=mddat.dtype)
    assert lal.swig_lal_test_viewin_gsl_matrix(mdout, mddat)
    assert (mdout.data == mddat).all()
    mdout.data = numpy.zeros(numpy.shape(mddat), dtype=mddat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_matrix(mdout, md)
    assert (2 * mdout.data == md.data).all()
    mdout.data = numpy.zeros(numpy.shape(mddat), dtype=mddat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_matrix(mdout, mddat)
    assert (2 * mdout.data == mddat).all()
    md.data = mddat
    assert lal.swig_lal_test_copyinout_gsl_matrix(md)
    assert (md.data == 3 * mddat).all()
    md.data = mddat
    retn, md = lal.swig_lal_test_copyinout_gsl_matrix(md)
    assert retn
    assert (md.data == 3 * mddat).all()
    md = mddat
    retn, md = lal.swig_lal_test_copyinout_gsl_matrix(md)
    assert retn
    assert (md == 3 * mddat).all()
    del md
    del mdout
    del mddat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    mcf = lal.gsl_matrix_complex_float(mcfdat.shape[0], mcfdat.shape[1])
    mcf.data = mcfdat
    mcfout = lal.gsl_matrix_complex_float(mcfdat.shape[0], mcfdat.shape[1])
    mcfout.data = numpy.zeros(numpy.shape(mcfdat), dtype=mcfdat.dtype)
    assert lal.swig_lal_test_viewin_gsl_matrix_complex_float(mcfout, mcf)
    assert (mcfout.data == mcf.data).all()
    mcfout.data = numpy.zeros(numpy.shape(mcfdat), dtype=mcfdat.dtype)
    assert lal.swig_lal_test_viewin_gsl_matrix_complex_float(mcfout, mcfdat)
    assert (mcfout.data == mcfdat).all()
    mcfout.data = numpy.zeros(numpy.shape(mcfdat), dtype=mcfdat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_matrix_complex_float(mcfout, mcf)
    assert (2 * mcfout.data == mcf.data).all()
    mcfout.data = numpy.zeros(numpy.shape(mcfdat), dtype=mcfdat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_matrix_complex_float(mcfout, mcfdat)
    assert (2 * mcfout.data == mcfdat).all()
    mcf.data = mcfdat
    assert lal.swig_lal_test_copyinout_gsl_matrix_complex_float(mcf)
    assert (mcf.data == 3 * mcfdat).all()
    mcf.data = mcfdat
    retn, mcf = lal.swig_lal_test_copyinout_gsl_matrix_complex_float(mcf)
    assert retn
    assert (mcf.data == 3 * mcfdat).all()
    mcf = mcfdat
    retn, mcf = lal.swig_lal_test_copyinout_gsl_matrix_complex_float(mcf)
    assert retn
    assert (mcf == 3 * mcfdat).all()
    del mcf
    del mcfout
    del mcfdat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    mcd = lal.gsl_matrix_complex(mcddat.shape[0], mcddat.shape[1])
    mcd.data = mcddat
    mcdout = lal.gsl_matrix_complex(mcddat.shape[0], mcddat.shape[1])
    mcdout.data = numpy.zeros(numpy.shape(mcddat), dtype=mcddat.dtype)
    assert lal.swig_lal_test_viewin_gsl_matrix_complex(mcdout, mcd)
    assert (mcdout.data == mcd.data).all()
    mcdout.data = numpy.zeros(numpy.shape(mcddat), dtype=mcddat.dtype)
    assert lal.swig_lal_test_viewin_gsl_matrix_complex(mcdout, mcddat)
    assert (mcdout.data == mcddat).all()
    mcdout.data = numpy.zeros(numpy.shape(mcddat), dtype=mcddat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_matrix_complex(mcdout, mcd)
    assert (2 * mcdout.data == mcd.data).all()
    mcdout.data = numpy.zeros(numpy.shape(mcddat), dtype=mcddat.dtype)
    assert lal.swig_lal_test_viewinout_gsl_matrix_complex(mcdout, mcddat)
    assert (2 * mcdout.data == mcddat).all()
    mcd.data = mcddat
    assert lal.swig_lal_test_copyinout_gsl_matrix_complex(mcd)
    assert (mcd.data == 3 * mcddat).all()
    mcd.data = mcddat
    retn, mcd = lal.swig_lal_test_copyinout_gsl_matrix_complex(mcd)
    assert retn
    assert (mcd.data == 3 * mcddat).all()
    mcd = mcddat
    retn, mcd = lal.swig_lal_test_copyinout_gsl_matrix_complex(mcd)
    assert retn
    assert (mcd == 3 * mcddat).all()
    del mcd
    del mcdout
    del mcddat
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED input views of numeric array structs (GSL)", file=sys.stderr)
    def check_input_view_type_safety(f, a, b, expect_exception):
        expected_exception = False
        if expect_exception:
            set_nice_error_handlers()
            expected_exception = False
            try:
                f(a, b)
                expected_exception = True
            except:
                pass
            set_default_error_handlers()
            assert not expected_exception
            set_nice_error_handlers()
            expected_exception = False
            try:
                f(b, a)
                expected_exception = True
            except:
                pass
            set_default_error_handlers()
            assert not expected_exception
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
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED input views of numeric array structs (type safety)", file=sys.stderr)

def test_FFT_functions_with_input_views():
    """check FFT functions with input views
    """

    print("check FFT functions with input views ...", file=sys.stderr)
    r4in = numpy.array(list(range(0, 32)), dtype=numpy.float32)
    r8in = numpy.array(list(range(0, 64)), dtype=numpy.float64)
    c8in = numpy.array(numpy.vectorize(complex)(8 + r4in, r4in), dtype=numpy.complex64)
    c16in = numpy.array(numpy.vectorize(complex)(16 + r8in, r8in), dtype=numpy.complex128)
    c8inv = lal.CreateCOMPLEX8Vector(len(c8in))
    c8inv.data = c8in
    c8outv = lal.CreateCOMPLEX8Vector(len(c8in))
    plan = lal.CreateForwardCOMPLEX8FFTPlan(len(c8in), 0)
    lal.COMPLEX8VectorFFT(c8outv, c8inv, plan)
    c8out = numpy.zeros(numpy.shape(c8outv.data), dtype=c8outv.data.dtype)
    lal.COMPLEX8VectorFFT(c8out, c8in, plan)
    assert (c8out == c8outv.data).all()
    del c8inv
    del c8outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    c16inv = lal.CreateCOMPLEX16Vector(len(c16in))
    c16inv.data = c16in
    c16outv = lal.CreateCOMPLEX16Vector(len(c16in))
    plan = lal.CreateForwardCOMPLEX16FFTPlan(len(c16in), 0)
    lal.COMPLEX16VectorFFT(c16outv, c16inv, plan)
    c16out = numpy.zeros(numpy.shape(c16outv.data), dtype=c16outv.data.dtype)
    lal.COMPLEX16VectorFFT(c16out, c16in, plan)
    assert (c16out == c16outv.data).all()
    del c16inv
    del c16outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    r4inv = lal.CreateREAL4Vector(len(r4in))
    r4inv.data = r4in
    c8outv = lal.CreateCOMPLEX8Vector(len(r4in)//2 + 1)
    plan = lal.CreateForwardREAL4FFTPlan(len(r4in), 0)
    lal.REAL4ForwardFFT(c8outv, r4inv, plan)
    c8out = numpy.zeros(numpy.shape(c8outv.data), dtype=c8outv.data.dtype)
    lal.REAL4ForwardFFT(c8out, r4in, plan)
    assert (c8out == c8outv.data).all()
    del r4inv
    del c8outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    c8inv = lal.CreateCOMPLEX8Vector(len(c8in))
    c8inv.data = c8in
    r4outv = lal.CreateREAL4Vector((len(c8in)-1)*2)
    plan = lal.CreateReverseREAL4FFTPlan((len(c8in)-1)*2, 0)
    lal.REAL4ReverseFFT(r4outv, c8inv, plan)
    r4out = numpy.zeros(numpy.shape(r4outv.data), dtype=r4outv.data.dtype)
    lal.REAL4ReverseFFT(r4out, c8in, plan)
    assert (r4out == r4outv.data).all()
    del c8inv
    del r4outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    r8inv = lal.CreateREAL8Vector(len(r8in))
    r8inv.data = r8in
    c16outv = lal.CreateCOMPLEX16Vector(len(r8in)//2 + 1)
    plan = lal.CreateForwardREAL8FFTPlan(len(r8in), 0)
    lal.REAL8ForwardFFT(c16outv, r8inv, plan)
    c16out = numpy.zeros(numpy.shape(c16outv.data), dtype=c16outv.data.dtype)
    lal.REAL8ForwardFFT(c16out, r8in, plan)
    assert (c16out == c16outv.data).all()
    del r8inv
    del c16outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    c16inv = lal.CreateCOMPLEX16Vector(len(c16in))
    c16inv.data = c16in
    r8outv = lal.CreateREAL8Vector((len(c16in)-1)*2)
    plan = lal.CreateReverseREAL8FFTPlan((len(c16in)-1)*2, 0)
    lal.REAL8ReverseFFT(r8outv, c16inv, plan)
    r8out = numpy.zeros(numpy.shape(r8outv.data), dtype=r8outv.data.dtype)
    lal.REAL8ReverseFFT(r8out, c16in, plan)
    assert (r8out == r8outv.data).all()
    del c16inv
    del r8outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    r4inv = lal.CreateREAL4Vector(len(r4in))
    r4inv.data = r4in
    r4outv = lal.CreateREAL4Vector(len(r4in))
    plan = lal.CreateForwardREAL4FFTPlan(len(r4in), 0)
    lal.REAL4VectorFFT(r4outv, r4inv, plan)
    r4out = numpy.zeros(numpy.shape(r4outv.data), dtype=r4outv.data.dtype)
    lal.REAL4VectorFFT(r4out, r4in, plan)
    assert (r4out == r4outv.data).all()
    del r4inv
    del r4outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    r8inv = lal.CreateREAL8Vector(len(r8in))
    r8inv.data = r8in
    r8outv = lal.CreateREAL8Vector(len(r8in))
    plan = lal.CreateForwardREAL8FFTPlan(len(r8in), 0)
    lal.REAL8VectorFFT(r8outv, r8inv, plan)
    r8out = numpy.zeros(numpy.shape(r8outv.data), dtype=r8outv.data.dtype)
    lal.REAL8VectorFFT(r8out, r8in, plan)
    assert (r8out == r8outv.data).all()
    del r8inv
    del r8outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    r4inv = lal.CreateREAL4Vector(len(r4in))
    r4inv.data = r4in
    r4outv = lal.CreateREAL4Vector(len(r4in)//2 + 1)
    plan = lal.CreateForwardREAL4FFTPlan(len(r4in), 0)
    lal.REAL4PowerSpectrum(r4outv, r4inv, plan)
    r4out = numpy.zeros(numpy.shape(r4outv.data), dtype=r4outv.data.dtype)
    lal.REAL4PowerSpectrum(r4out, r4in, plan)
    assert (r4out == r4outv.data).all()
    del r4inv
    del r4outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    r8inv = lal.CreateREAL8Vector(len(r8in))
    r8inv.data = r8in
    r8outv = lal.CreateREAL8Vector(len(r8in)//2 + 1)
    plan = lal.CreateForwardREAL8FFTPlan(len(r8in), 0)
    lal.REAL8PowerSpectrum(r8outv, r8inv, plan)
    r8out = numpy.zeros(numpy.shape(r8outv.data), dtype=r8outv.data.dtype)
    lal.REAL8PowerSpectrum(r8out, r8in, plan)
    assert (r8out == r8outv.data).all()
    del r8inv
    del r8outv
    del plan
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED FFT functions with input views ...", file=sys.stderr)

def test_dynamic_array_of_pointers_access():
    """check dynamic array of pointers access
    """

    print("checking dynamic array of pointers access ...", file=sys.stderr)
    ap = lal.swig_lal_test_Create_arrayofptrs(3)
    assert ap.length == 3
    for i in range(0, ap.length):
        assert ap.data[i].length == 6
        for j in range(0, ap.data[i].length):
            assert ap.data[i].data[j] == 42*ap.length*i + j
    del ap
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED dynamic array of pointers access", file=sys.stderr)

def test_typemaps_for_strings_and_double_pointers():
    """check typemaps for strings and double pointers
    """

    print("checking typemaps for strings and double pointers ...", file=sys.stderr)
    sts = lal.swig_lal_test_struct()
    ptr_ptr, ptr_null_ptr, null_ptr_ptr = lal.swig_lal_test_typemaps_string_ptrptr("abcde", "", None, sts, 0, None)
    assert ptr_ptr == sts
    assert ptr_null_ptr is not None
    assert null_ptr_ptr is None
    del sts
    del ptr_ptr
    del ptr_null_ptr
    del null_ptr_ptr
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    ptr_ptr = 0
    for i in range(1, 10):
        ptr_ptr = lal.swig_lal_test_typemaps_ptrptr(ptr_ptr)
        assert ptr_ptr is not None
        assert ptr_ptr.n == i
    del ptr_ptr
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    ptr_ptr_list = [0]
    for i in range(1, 10):
        ptr_ptr_list.append(lal.swig_lal_test_typemaps_ptrptr(ptr_ptr_list[-1]))
        assert ptr_ptr_list[-1] is not None
        assert ptr_ptr_list[-1].n == i
    while len(ptr_ptr_list) > 0:
        assert ptr_ptr_list[-1] is not None
        assert ptr_ptr_list[-1].n == i
        del ptr_ptr_list[0]
    del ptr_ptr_list
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED typemaps for strings and double pointers", file=sys.stderr)

def test_tm_struct_conversions():
    """check 'tm' struct conversions
    """

    print("checking 'tm' struct conversions ...", file=sys.stderr)
    gps0 = 989168284
    utc0 = [2011, 5, 11, 16, 57, 49, 2, 131, 0]
    assert lal.GPSToUTC(gps0) == tuple(utc0)
    assert lal.UTCToGPS(utc0) == gps0
    for i in range(0, 10):
        gps = gps0 + i * 86400
        utc = list(utc0)
        utc[2] = utc[2] + i
        utc[6] = (utc[6] + i) % 7
        utc[7] = utc[7] + i
        utc[8] = -1 + (i % 3)
        assert lal.GPSToUTC(gps)[0:8] == tuple(utc[0:8])
        assert lal.UTCToGPS(utc) == gps
        utc = lal.GPSToUTC(lal.UTCToGPS(utc))
        dt = datetime.datetime(*utc[0:6])
        assert utc[6] == dt.weekday()
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED 'tm' struct conversions", file=sys.stderr)

def test_LIGOTimeGPS_operations():
    """check LIGOTimeGPS operations
    """

    print("checking LIGOTimeGPS operations ...", file=sys.stderr)
    from lal import LIGOTimeGPS
    t0 = LIGOTimeGPS()
    assert type(LIGOTimeGPS(t0)) is LIGOTimeGPS
    assert is_value_and_type(t0, 0, LIGOTimeGPS)
    assert t0 != None and not t0 is None
    t1 = LIGOTimeGPS(10.5)
    t2 = LIGOTimeGPS(10, 500000000)
    assert not t0 and t1 and t2
    assert is_value_and_type(t1, t2, LIGOTimeGPS)
    t3 = +t1
    t3 = -t2
    assert t1 == t2 and t1 >= t2 and t2 >= t1
    assert abs(-t1) == t1
    assert float(t1) == 10.5
    assert is_value_and_type(t1 + 3.5, 14, LIGOTimeGPS)
    assert is_value_and_type(3.5 + t1, 14, LIGOTimeGPS)
    t2 -= 5.5
    assert is_value_and_type(t2, 5, LIGOTimeGPS)
    assert t2 + 5.5 >= t1 and t2 + 3 != t2
    assert is_value_and_type(t2 - 5, t0, LIGOTimeGPS)
    assert is_value_and_type(t1 * 3, 31.5, LIGOTimeGPS)
    assert is_value_and_type(3 * t1, 31.5, LIGOTimeGPS)
    assert is_value_and_type(t2 / 2.5, 2, LIGOTimeGPS)
    assert is_value_and_type(21 / t1, 2, LIGOTimeGPS)
    assert is_value_and_type(t1 + t2, 15.5, LIGOTimeGPS)
    assert is_value_and_type(t1 - t2, 5.5, LIGOTimeGPS)
    assert is_value_and_type(t1 * t2, 52.5, LIGOTimeGPS)
    assert is_value_and_type(t2 * t1, 52.5, LIGOTimeGPS)
    assert is_value_and_type(t1 / t2, 2.1, LIGOTimeGPS)
    assert is_value_and_type(t1 % t2, 0.5, LIGOTimeGPS)
    assert t1 > t2 and t2 < t1 and t1 >= t2 and t2 <= t1
    assert LIGOTimeGPS(333333333,333333333) == LIGOTimeGPS(1000000000) / 3
    assert LIGOTimeGPS(666666666,666666667) == LIGOTimeGPS(2000000000) / 3
    assert LIGOTimeGPS("-62997760.825036067") == LIGOTimeGPS("-47044285.062262587") - LIGOTimeGPS("15953475.76277348")
    assert LIGOTimeGPS("-6542354.389038577") == LIGOTimeGPS("-914984.929117316") * 7.1502318572066237
    assert LIGOTimeGPS("-6542354.389038577") == 7.1502318572066237 * LIGOTimeGPS("-914984.929117316")
    assert LIGOTimeGPS("-127965.770535834") == LIGOTimeGPS("-914984.929117316") / 7.1502318572066237
    t1 += 812345667.75
    assert str(t1) == "812345678.25"
    assert type(eval(repr(t1))) is type(t1)
    assert eval(repr(t1)) == t1
    assert(LIGOTimeGPS(1100000000).asutcstr() == "Fri, 14 Nov 2014 11:33:04 +0000")   # lalapps_tconvert -u -R
    assert LIGOTimeGPS(1100000000, 100).asutcstr() == "Fri, 14 Nov 2014 11:33:04.0000001 +0000"
    assert LIGOTimeGPS(0, 0).asutcstr() == "Sun, 06 Jan 1980 00:00:00 +0000"
    assert LIGOTimeGPS(-1, 0).asutcstr() == "Sat, 05 Jan 1980 23:59:59 +0000"
    assert LIGOTimeGPS(0, -1).asutcstr() == "Sat, 05 Jan 1980 23:59:59.999999999 +0000"
    assert int(t1) == 812345678
    assert t1.ns() == 812345678250000000
    assert hash(t1) == 1049484238
    t4struct = lal.swig_lal_test_gps()
    t4struct.t = 1234.5
    assert t4struct.t == 1234.5
    t5 = LIGOTimeGPS("1000")
    assert t5 == 1000
    print("*** below should be error messages from LIGOTimeGPS constructor ***", file=sys.stderr)
    set_nice_error_handlers()
    expected_exception = False
    try:
        t5 = LIGOTimeGPS("abc1000")
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    set_nice_error_handlers()
    expected_exception = False
    try:
        t5 = LIGOTimeGPS("1000abc")
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    print("*** above should be error messages from LIGOTimeGPS constructor ***", file=sys.stderr)
    assert lal.swig_lal_test_noptrgps(LIGOTimeGPS(1234.5)) == lal.swig_lal_test_noptrgps(1234.5)
    print("*** below should be error messages from LIGOTimeGPS constructor ***", file=sys.stderr)
    set_nice_error_handlers()
    expected_exception = False
    try:
        LIGOTimeGPS(None)
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    print("*** above should be error messages from LIGOTimeGPS constructor ***", file=sys.stderr)
    del t0
    del t1
    del t2
    del t3
    del t4struct
    del t5
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED LIGOTimeGPS operations", file=sys.stderr)
    print("checking LIGOTimeGPS operations (Python specific) ...", file=sys.stderr)
    class my_gps_class:
        def __init__(self, s, ns):
            self.gpsSeconds = s
            self.gpsNanoSeconds = ns
    tmy = my_gps_class(987, 654321)
    tsw = LIGOTimeGPS(tmy)
    assert type(tsw) is LIGOTimeGPS
    assert tsw.gpsSeconds == 987
    assert tsw.gpsNanoSeconds == 654321
    assert is_value_and_type(tsw + tmy, tsw + tsw, LIGOTimeGPS)
    assert is_value_and_type(tmy + tsw, tsw + tsw, LIGOTimeGPS)
    assert is_value_and_type(tsw - tmy, tsw - tsw, LIGOTimeGPS)
    assert is_value_and_type(tmy - tsw, tsw - tsw, LIGOTimeGPS)
    assert is_value_and_type(tsw * tmy, tsw * tsw, LIGOTimeGPS)
    assert is_value_and_type(tmy * tsw, tsw * tsw, LIGOTimeGPS)
    assert is_value_and_type(tsw / tmy, tsw / tsw, LIGOTimeGPS)
    assert is_value_and_type(tmy / tsw, tsw / tsw, LIGOTimeGPS)
    assert lal.swig_lal_test_noptrgps(tmy) == lal.swig_lal_test_noptrgps(tsw)
    del tsw
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED LIGOTimeGPS operations (Python specific)", file=sys.stderr)

def test_LALUnit_operations():
    """check LALUnit operations
    """

    print("checking LALUnit operations ...", file=sys.stderr)
    u1 = lal.Unit("kg m s^-2")
    assert type(lal.Unit(u1)) is lal.Unit
    assert is_value_and_type(u1, lal.NewtonUnit, lal.Unit)
    assert str(u1) == "m kg s^-2"
    u2 = lal.MeterUnit * lal.KiloGramUnit / lal.SecondUnit ** 2
    assert is_value_and_type(u2, u1, lal.Unit)
    u2 = lal.MeterUnit**(1,2) * lal.KiloGramUnit**(1,2) * lal.SecondUnit ** -1
    assert is_value_and_type(u2, u1**(1,2), lal.Unit)
    set_nice_error_handlers()
    expected_exception = False
    try:
        lal.SecondUnit ** (1,0)
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    u1 *= lal.MeterUnit
    assert is_value_and_type(u1, lal.JouleUnit, lal.Unit)
    assert repr(u1) == "m^2 kg s^-2"
    u1 /= lal.SecondUnit
    assert is_value_and_type(u1, lal.WattUnit, lal.Unit)
    assert u1 == "m^2 kg s^-3"
    u1 *= 1000
    assert u1 == lal.KiloUnit * lal.WattUnit
    assert u1 == 1000 * lal.WattUnit
    assert u1 == lal.WattUnit * 1000
    assert u1 == lal.MegaUnit / 1000 * lal.WattUnit
    assert int(u1) == 1000
    u1 /= 10000
    assert u1 == 100 * lal.MilliUnit * lal.WattUnit
    set_nice_error_handlers()
    expected_exception = False
    try:
        u1 *= 1.234
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    assert u1.norm() == u1
    del u1
    del u2
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED LALUnit operations", file=sys.stderr)

def test_Python_non_dynamic_structs():
    """check Python non-dynamic structs
    """

    print("checking Python non-dynamic structs", file=sys.stderr)
    sts = lal.swig_lal_test_struct()
    sts.n = 1
    sts.Alpha = 1.234
    set_nice_error_handlers()
    expected_exception = False
    try:
        sts.N = 1
        sts.alpha = 1.234
        sts.zzzzz = "does not exist"
        expected_exception = True
    except:
        pass
    set_default_error_handlers()
    assert not expected_exception
    del sts
    t = lal.LIGOTimeGPS(1234)
    t.event = "this happened"
    del t
    locals()   # update locals() to remove lingering references
    lal.CheckMemoryLeaks()
    print("PASSED Python non-dynamic structs", file=sys.stderr)

def test_Python_pickling():
    """check Python pickling
    """

    print("checking Python pickling ...", file=sys.stderr)
    for datatype in ['INT2', 'INT4', 'INT8', 'UINT2', 'UINT4', 'UINT8',
                     'REAL4', 'REAL8', 'COMPLEX8', 'COMPLEX16']:

        creator = getattr(lal, 'Create{}FrequencySeries'.format(datatype))
        a = creator(
            'foobar', lal.LIGOTimeGPS(1e9), 0, 1, lal.StrainUnit, 1024)
        a.data.data = numpy.arange(1024)
        pickled = pickle.dumps(a, 0)
        b = pickle.loads(pickled)
        assert type(a) == type(b)
        assert a.name == b.name
        assert a.epoch == b.epoch
        assert a.f0 == b.f0
        assert a.deltaF == b.deltaF
        assert a.sampleUnits == b.sampleUnits
        assert (a.data.data == b.data.data).all()

        creator = getattr(lal, 'Create{}TimeSeries'.format(datatype))
        a = creator(
            'foobar', lal.LIGOTimeGPS(1e9), 0, 1, lal.StrainUnit, 1024)
        a.data.data = numpy.arange(1024)
        pickled = pickle.dumps(a, 0)
        b = pickle.loads(pickled)
        assert type(a) == type(b)
        assert a.name == b.name
        assert a.epoch == b.epoch
        assert a.f0 == b.f0
        assert a.deltaT == b.deltaT
        assert a.sampleUnits == b.sampleUnits
        assert (a.data.data == b.data.data).all()
    print("PASSED Python pickling", file=sys.stderr)

def test_Python_dict_to_LALDict_typemap():
    """check Python dict to LALDict typemap
    """

    print("checking Python Python dict to LALDict typemap ...", file=sys.stderr)
    pydict = {
        "str": "A string value",
        "2-byte-unsigned:UINT2": 32767,
        "4-byte-unsigned:UINT4": 123456,
        "8-byte-unsigned:UINT8": 9223372036854775807,
        "2-byte-signed:INT2": -32768,
        "4-byte-signed:INT4": -123456,
        "8-byte-signed:INT8": 9223372036854775807,
        "single:REAL4": 987e6,
        "double:REAL8": -543e-21,
        "single complex:COMPLEX8": complex(987e6, -123e4),
        "double complex:COMPLEX16": complex(-543e-21, 345e43)
    }
    laldict = lal.CreateDict()
    lal.DictInsertStringValue(laldict, "str", pydict["str"])
    lal.DictInsertUINT2Value(laldict, "2-byte-unsigned", pydict["2-byte-unsigned:UINT2"])
    lal.DictInsertUINT4Value(laldict, "4-byte-unsigned", pydict["4-byte-unsigned:UINT4"])
    lal.DictInsertUINT8Value(laldict, "8-byte-unsigned", pydict["8-byte-unsigned:UINT8"])
    lal.DictInsertINT2Value(laldict, "2-byte-signed", pydict["2-byte-signed:INT2"])
    lal.DictInsertINT4Value(laldict, "4-byte-signed", pydict["4-byte-signed:INT4"])
    lal.DictInsertINT8Value(laldict, "8-byte-signed", pydict["8-byte-signed:INT8"])
    lal.DictInsertREAL4Value(laldict, "single", pydict["single:REAL4"])
    lal.DictInsertREAL8Value(laldict, "double", pydict["double:REAL8"])
    lal.DictInsertCOMPLEX8Value(laldict, "single complex", pydict["single complex:COMPLEX8"])
    lal.DictInsertCOMPLEX16Value(laldict, "double complex", pydict["double complex:COMPLEX16"])
    lal.swig_lal_test_pydict_to_laldict(laldict)
    lal.swig_lal_test_pydict_to_laldict(pydict)
    print("PASSED Python dict to LALDict typemap", file=sys.stderr)

def test_Python_conversion_of_NumPy_fixed_width_integer_float_types():
    """check Python conversion of NumPy fixed-width integer/float types
    """

    print("checking Python conversion of NumPy fixed-width integer/float types", file=sys.stderr)
    assert lal.swig_lal_test_numpy_int_types(10, numpy.int16(20), numpy.int32(30), numpy.int64(-40)) == 20
    assert lal.swig_lal_test_numpy_int_types(numpy.int8(10), numpy.int16(20), numpy.int32(30), numpy.int64(-40)) == 20
    assert lal.swig_lal_test_numpy_int_types(numpy.int16(10), numpy.int16(20), numpy.int32(30), numpy.int64(-40)) == 20
    assert lal.swig_lal_test_numpy_int_types(numpy.int32(10), numpy.int16(20), numpy.int32(30), numpy.int64(-40)) == 20
    assert lal.swig_lal_test_numpy_int_types(numpy.int64(10), numpy.int16(20), numpy.int32(30), numpy.int64(-40)) == 20
    assert lal.swig_lal_test_numpy_uint_types(10, numpy.uint16(20), numpy.uint32(30), numpy.uint64(40)) == 100
    assert lal.swig_lal_test_numpy_uint_types(numpy.uint8(10), numpy.uint16(20), numpy.uint32(30), numpy.uint64(40)) == 100
    assert lal.swig_lal_test_numpy_uint_types(numpy.uint16(10), numpy.uint16(20), numpy.uint32(30), numpy.uint64(40)) == 100
    assert lal.swig_lal_test_numpy_uint_types(numpy.uint32(10), numpy.uint16(20), numpy.uint32(30), numpy.uint64(40)) == 100
    assert lal.swig_lal_test_numpy_uint_types(numpy.uint64(10), numpy.uint16(20), numpy.uint32(30), numpy.uint64(40)) == 100
    assert lal.swig_lal_test_numpy_flt_types(numpy.int8(10), numpy.int16(20), numpy.int32(30), numpy.int64(-40)) == 20
    assert lal.swig_lal_test_numpy_flt_types(numpy.uint8(10), numpy.uint16(20), numpy.uint32(30), numpy.uint64(40)) == 100
    assert lal.swig_lal_test_numpy_flt_types(numpy.float_(10), numpy.float16(20), numpy.float32(30), numpy.float64(40)) == 100
    assert lal.swig_lal_test_numpy_cpx_types(numpy.int8(10), numpy.int16(20), numpy.int32(30), numpy.int64(-40)) == 20
    assert lal.swig_lal_test_numpy_cpx_types(numpy.uint8(10), numpy.uint16(20), numpy.uint32(30), numpy.uint64(40)) == 100
    assert lal.swig_lal_test_numpy_cpx_types(numpy.float_(10), numpy.float16(20), numpy.float32(30), numpy.float64(40)) == 100
    assert lal.swig_lal_test_numpy_cpx_types(numpy.complex_(10), numpy.complex64(20), numpy.complex64(30), numpy.complex128(40)) == 100
    assert int(lal.LIGOTimeGPS(numpy.int8(123))) == 123
    assert int(lal.LIGOTimeGPS(numpy.int16(123))) == 123
    assert int(lal.LIGOTimeGPS(numpy.int32(123))) == 123
    assert int(lal.LIGOTimeGPS(numpy.int64(123))) == 123
    print("PASSED Python conversion of NumPy fixed-width integer/float types", file=sys.stderr)

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-SWIGTestLALPython.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
