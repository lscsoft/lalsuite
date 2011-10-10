## check SWIG Octave module wrapping the LAL library
## Author: Karl Wette, 2011

function msg(str)
  disp([program_name, ": ", str]);
endfunction

## check module load
swiglal;
swiglal.cvar.lalDebugLevel = 1;
msg("passed module load");

## check memory allocation
if !cvar.swiglal_debug
  msg("skipping memory allocation");
else
  LALCheckMemoryLeaks();
  mem1 = new_LALDetector();
  mem2 = new_LALStringVector();
  mem3 = XLALCreateCOMPLEX8Vector(5);
  mem4 = XLALCreateREAL8Vector(3);
  msg("*** below should be an error message from LALCheckMemoryLeaks() ***");
  try
    LALCheckMemoryLeaks();
    error("expected exception");
  end_try_catch
  msg("*** above should be an error message from LALCheckMemoryLeaks() ***");
  clear mem1 mem2 mem3 mem4;
  LALCheckMemoryLeaks();
  msg("passed memory allocation");
endif

## check complex number conversions
assert(XLALCOMPLEX8Add(complex(1, 3), complex(2, -5)) == complex(3, -2));
assert(XLALCOMPLEX8Sub(complex(4, 2), complex(10, 5)) == complex(-6, -3));
assert(XLALCOMPLEX16Mul(complex(10, -7), complex(30, 17)) == complex(419, -40));
assert(XLALCOMPLEX16Div(complex(111.75, -120.50), complex(5, 12)) == complex(-5.25, -11.5));
msg("passed complex number conversions");

## check string conversions
strs = {"a", "bc", "def"};
sv = XLALCreateStringVector(strs{:});
assert(sv.length == 3);
assert(all(strcmp(sv.data, strs)));
strs{1} = "ghijk";
sv.data_setel(0, strs{1});
strs{end+1} = "lmnopq";
XLALAppendString2Vector(sv, strs{4});
assert(sv.length == 4);
for i = 1:4
  assert(strcmp(sv.data_getel(i-1), strs{i}));
endfor
XLALDestroyStringVector(sv);
msg("passed string conversions");

## check vector/matrix struct type accessors
if !cvar.swiglal_debug
  msg("skipping vector/matrix struct type accessors");
else
  swiglal_test_struct_vector_setel(0, cvar.swiglal_test_struct_const);
  assert(swiglal_test_struct_vector_getel(0).a == cvar.swiglal_test_struct_const.a);
  assert(swiglal_test_struct_vector_getel(0).b == cvar.swiglal_test_struct_const.b);
  assert(strcmp(swiglal_test_struct_vector_getel(0).c, cvar.swiglal_test_struct_const.c));
  swiglal_test_struct_matrix_setel(0, 0, cvar.swiglal_test_struct_const);
  assert(swiglal_test_struct_matrix_getel(0, 0).a == cvar.swiglal_test_struct_const.a);
  assert(swiglal_test_struct_matrix_getel(0, 0).b == cvar.swiglal_test_struct_const.b);
  assert(strcmp(swiglal_test_struct_matrix_getel(0, 0).c, cvar.swiglal_test_struct_const.c));
  msg("passed vector/matrix struct type accessors");
endif

## check static vector/matrix conversions
if !cvar.swiglal_debug
  msg("skipping static vector/matrix conversions");
else
  sts = new_swiglal_test_static_struct();
  assert(length(sts.vector) == 3);
  assert(length(sts.enum_vector) == 3);
  assert(all(size(sts.matrix) == [2, 3]));
  assert(all(size(sts.enum_matrix) == [2, 3]));
  sts.vector = [3, 2, 1];
  assert(all(sts.vector == [3, 2, 1]));
  sts.matrix = [4, 5, 6; 9, 8, 7];
  try
    sts.matrix = [1.1, 2.3, 4.5; 6.5, 4.3, 2.1];
    error("expected exception");
  end_try_catch
  assert(all(all(sts.matrix == [4, 5, 6; 9, 8, 7])));
  for i = 1:3
    sts.enum_vector_setel(i-1, 2*i + 3);
    assert(sts.enum_vector_getel(i-1) == (2*i + 3));
  endfor
  clear sts;
  assert(!any(cvar.swiglal_test_static_vector));
  assert(!any(cvar.swiglal_test_static_matrix(:)));
  assert(!any(cvar.swiglal_test_static_enum_vector));
  assert(!any(cvar.swiglal_test_static_enum_matrix(:)));
  swiglal_test_static_vector_setel(0, 10);
  assert(swiglal_test_static_vector_getel(0) == 10);
  swiglal_test_static_matrix_setel(0, 0, 11);
  assert(swiglal_test_static_matrix_getel(0, 0) == 11);
  cvar.swiglal_test_static_vector = cvar.swiglal_test_static_const_vector;
  assert(all(cvar.swiglal_test_static_vector == [1, 2, 4]));
  assert(swiglal_test_static_const_vector_getel(2) == 4);
  cvar.swiglal_test_static_matrix = cvar.swiglal_test_static_const_matrix;
  assert(all(cvar.swiglal_test_static_matrix == [[1, 2, 4]; [2, 4, 8]]));
  assert(swiglal_test_static_const_matrix_getel(1, 2) == 8);
  try
    swiglal_test_static_const_vector_getel(20);
    error("expected exception");
  end_try_catch
  msg("passed static vector/matrix conversions");
endif

## check dynamic vector/matrix conversions
function check_dynamic_vector_matrix(iv, ivl, rv, rvl, cm, cms1, cms2)
  assert(ivl == 5);
  iv.data = [1, 3, 2, 4, 3];
  assert(all(iv.data == [1, 3, 2, 4, 3]));
  iv.data_setel(3, 7);
  assert(iv.data_getel(3) == 7);
  assert(rvl == 5);
  rv.data = [1.2, 3.4, 2.6, 4.8, 3.5];
  assert(all(rv.data == [1.2, 3.4, 2.6, 4.8, 3.5]));
  rv.data_setel(rvl - 1, 7.5);
  assert(rv.data_getel(rvl - 1) == 7.5);
  try
    rv.data_setel(rvl, 99.9);
    error("expected exception");
  end_try_catch
  try
    iv.data = rv.data;
    error("expected exception");
  end_try_catch
  rv.data = iv.data;
  assert(all(rv.data == iv.data));
  assert(cms1 == 4);
  assert(cms2 == 6);
  for i = 0:cms1-1
    for j = 0:cms2-1
      cm.data_setel(i, j, complex(i / 4.0, j / 2.0));
    endfor
  endfor
  assert(cm.data_getel(2, 3) == complex(0.5, 1.5));
  assert(cm.data_getel(3, 2) == complex(0.75, 1.0));
  try
    iv.data_setel(0, cm.data_getel(2, 3));
    error("expected exception");
  end_try_catch
  try
    rv.data_setel(0, cm.data_getel(3, 2));
    error("expected exception");
  end_try_catch
endfunction
## check LAL vector and matrix datatypes
iv = XLALCreateINT4Vector(5);
rv = XLALCreateREAL8Vector(5);
cm = XLALCreateCOMPLEX8VectorSequence(4, 6);
check_dynamic_vector_matrix(iv, iv.length, rv, rv.length,
                            cm, cm.length, cm.vectorLength);
clear iv rv cm;
LALCheckMemoryLeaks();
msg("passed dynamic vector/matrix conversions (LAL)");
## check GSL vectors and matrices
iv = new_gsl_vector_int(5);
rv = new_gsl_vector(5);
cm = new_gsl_matrix_complex_float(4, 6);
check_dynamic_vector_matrix(iv, iv.size, rv, rv.size,
                            cm, cm.size1, cm.size2);
clear iv rv cm;
msg("passed dynamic vector/matrix conversions (GSL)");

## check 'tm' struct conversions
gps = 989168284;
utc = [2011, 5, 11, 16, 57, 49, 4, 131, 0];
assert(all(XLALGPSToUTC([], gps) == utc));
assert(XLALUTCToGPS(utc) == gps);
assert(XLALUTCToGPS(utc(1:6)) == gps);
utc(7) = utc(8) = 0;
for i = [-1, 0, 1]
  utc(9) = i;
  assert(XLALUTCToGPS(utc) == gps);
endfor
utcd = utc;
for i = 0:9
  utcd(3) = utc(3) + i;
  utcd = XLALGPSToUTC([], XLALUTCToGPS(utcd));
  assert(utcd(7) == weekday(datenum(utcd(1:6))));
endfor
msg("passed 'tm' struct conversions");

## check LIGOTimeGPS operations
if swiglal.cvar.swig_version_hex < 0x020004
  global op_any_add_LIGOTimeGPS = inline("LIGOTimeGPS___radd__(y,x)", "x", "y")
  global op_any_sub_LIGOTimeGPS = inline("LIGOTimeGPS___rsub__(y,x)", "x", "y")
  global op_any_mul_LIGOTimeGPS = inline("LIGOTimeGPS___rmul__(y,x)", "x", "y")
  global op_any_div_LIGOTimeGPS = inline("LIGOTimeGPS___rdiv__(y,x)", "x", "y")
endif
t0 = new_LIGOTimeGPS();
assert(t0 == 0 && strcmp(swig_type(t0), "LIGOTimeGPS"));
t1 = new_LIGOTimeGPS(10.5);
t2 = new_LIGOTimeGPS(10, 500000000);
assert(t1 == t2 && strcmp(swig_type(t1), "LIGOTimeGPS"));
+t1;
-t2;
assert(t1 == t2 && t1 >= t2 && t2 >= t1);
assert(t1 + 3.5 == 14 && strcmp(swig_type(t1 + 3.5), "LIGOTimeGPS"));
assert(3.5 + t1 == 14 && strcmp(swig_type(3.5 + t1), "LIGOTimeGPS"));
t2 -= 5.5;
assert(t2 == 5 && strcmp(swig_type(t2), "LIGOTimeGPS"));
assert(t2 + 5.5 >= t1 && t2 + 3 != t2);
assert(t2 - 5 == t0 && strcmp(swig_type(t2 - 5), "LIGOTimeGPS"));
assert(t1 * 3 == 31.5 && strcmp(swig_type(t1 * 3), "LIGOTimeGPS"));
assert(3 * t1 == 31.5 && strcmp(swig_type(3 * t1), "LIGOTimeGPS"));
assert(t2 / 2.5 == 2 && strcmp(swig_type(t2 / 2.5), "LIGOTimeGPS"));
assert(21 / t1  == 2 && strcmp(swig_type(21 / t1 ), "LIGOTimeGPS"));
assert(t1 + t2 == 15.5 && strcmp(swig_type(t1 + t2), "LIGOTimeGPS"));
assert(t1 - t2 == 5.5 && strcmp(swig_type(t1 - t2), "LIGOTimeGPS"));
assert(t1 * t2 == 52.5 && strcmp(swig_type(t1 * t2), "LIGOTimeGPS"));
assert(t2 * t1 == 52.5 && strcmp(swig_type(t2 * t1), "LIGOTimeGPS"));
assert(t1 / t2 == 2.1 && strcmp(swig_type(t1 / t2), "LIGOTimeGPS"));
assert(t1 > t2 && t2 < t1 && t1 >= t2 && t2 <= t1)
t1 += 812345667.75;
assert(strcmp(t1.__str__(), "812345678.250000000"));
assert(new_LIGOTimeGPS(t1.__str__()) == t1);
assert(t1.ns() == 812345678250000000);
msg("passed LIGOTimeGPS operations");

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
