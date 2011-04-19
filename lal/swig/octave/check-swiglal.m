## check SWIG Octave module wrapping the LAL library
## Author: Karl Wette, 2011

function msg(str)
  disp([program_name, ": ", str]);
endfunction

## check memory allocation
if !cvar.swiglal_debug
  msg("skipping memory allocation");
else
  LALCheckMemoryLeaks();
  mem1 = new_LALDetector();
  mem2 = new_LALStringVector();
  mem3 = new_COMPLEX8Vector();
  mem4 = XLALCreateREAL8Vector(3);
  msg("*** below should be an error message from LALCheckMemoryLeaks() ***");
  try
    LALCheckMemoryLeaks();
    error("expected exception");
  end_try_catch
  msg("*** above should be an error message from LALCheckMemoryLeaks() ***");
  clear mem1 mem2 mem3;
  XLALDestroyREAL8Vector(mem4);
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

## check static vector/matrix conversions
if !cvar.swiglal_debug
  msg("skipping static vector/matrix conversions");
else
  sts = new_swiglal_static_test_struct();
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
  assert(!any(cvar.swiglal_static_test_vector));
  assert(!any(cvar.swiglal_static_test_matrix(:)));
  assert(!any(cvar.swiglal_static_test_enum_vector));
  assert(!any(cvar.swiglal_static_test_enum_matrix(:)));
  cvar.swiglal_static_test_vector = cvar.swiglal_static_test_const_vector;
  assert(all(cvar.swiglal_static_test_vector == [1, 2, 4]));
  assert(swiglal_static_test_const_vector_getel(2) == 4);
  try
    swiglal_static_test_const_vector_getel(20);
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
XLALDestroyINT4Vector(iv);
XLALDestroyREAL8Vector(rv);
XLALDestroyCOMPLEX8VectorSequence(cm);
LALCheckMemoryLeaks();
msg("passed dynamic vector/matrix conversions (LAL)");

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
