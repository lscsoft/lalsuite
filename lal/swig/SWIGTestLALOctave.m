## Check SWIG Octave bindings for LAL
## Author: Karl Wette, 2011--2014

page_screen_output(0);
crash_dumps_octave_core(0);

# assertion helper functions
function assert_LIGOTimeGPS(t1, t2)
  if strcmp(typeinfo(t2), "swig_ref")
    assert(char(t1), char(t2));
  else
    assert(double(t1), t2);
  endif
endfunction
function assert_LALUnit(t1, t2)
  assert(char(t1), char(t2));
endfunction

## check module load
disp("checking module load ...");
lal;
assert(exist("lal"));
lal_c_si = LAL_C_SI;
lal_180_pi = LAL_180_PI;
disp("PASSED module load");

# set error handlers
function set_nice_error_handlers()
  swig_set_nice_error_handlers();
endfunction
function set_default_error_handlers()
  lal;
  if length(getenv("NASTY_ERROR_HANDLERS")) > 0
    swig_set_nasty_error_handlers();
  else
    swig_set_nice_error_handlers();
  endif
endfunction
set_default_error_handlers();

## check memory allocation
disp("checking memory allocation ...");
if !LAL_MEMORY_FUNCTIONS_DISABLED
  LALCheckMemoryLeaks();
  mem1 = new_LALDetector();
  mem2 = XLALCreateCOMPLEX8Vector(5);
  mem3 = XLALCreateREAL8Vector(3);
  mem4 = XLALCreateREAL4TimeSeries("test", LIGOTimeGPS(0), 100, 0.1, lal.DimensionlessUnit, 10);
  disp("*** below should be an error message from CheckMemoryLeaks() ***");
  set_nice_error_handlers();
  expected_exception = 0;
  try
    LALCheckMemoryLeaks();
    expected_exception = 1;
  end_try_catch
  set_default_error_handlers();
  assert(!expected_exception);
  disp("*** above should be an error message from CheckMemoryLeaks() ***");
  clear mem1;
  clear mem2;
  clear mem3;
  clear mem4;
  clear ans;
  LALCheckMemoryLeaks();
  disp("PASSED memory allocation");
else
  disp("skipped memory allocation");
endif

## check object parent tracking
disp("checking object parent tracking ...");
a = new_gsl_vector(3);
a.data = [1.1; 2.2; 3.3];
b = a.data;
assert(typeinfo(b), "swiglal_oct_array_view_double");
assert(b(:), [1.1; 2.2; 3.3]);
clear a;
clear ans;
assert(b(:), [1.1; 2.2; 3.3]);
ts = XLALCreateREAL8TimeSeries("test", LIGOTimeGPS(0), 0, 0.1, lal.DimensionlessUnit, 10);
ts.data.data = (1:10)';
for i = 1:7
  v = ts.data;
endfor
assert(v.data(:), (1:10)');
clear ts;
clear ans;
assert(v.data(:), (1:10)');
clear v;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED object parent tracking");

## check equal return/first argument type handling
disp("checking equal return/first argument type handling");
sv = XLALCreateStringVector("1");
assert(sv.length, 1);
XLALAppendString2Vector(sv, "2");
assert(sv.length, 2);
sv = XLALAppendString2Vector(sv, "3");
assert(sv.length, 3);
sv2 = XLALAppendString2Vector(sv, "4");
assert(sv.length, 4);
assert(sv2.length, 4);
assert(swig_this(sv), swig_this(sv2));
clear sv;
clear sv2;
clear ans;
LALCheckMemoryLeaks();
ts = XLALCreateREAL8TimeSeries("ts", 800000000, 100, 0.1, lal.HertzUnit, 10);
assert(ts.data.length == 10)
XLALResizeREAL8TimeSeries(ts, 0, 20);
assert(ts.data.length == 20)
ts = XLALResizeREAL8TimeSeries(ts, 0, 30);
assert(ts.data.length == 30)
ts2 = XLALResizeREAL8TimeSeries(ts, 0, 40);
assert(ts.data.length, 40);
assert(ts2.data.length, 40);
assert(swig_this(ts), swig_this(ts2));
clear ts;
clear ts2;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED equal return/first argument type handling");

## check string conversions
disp("checking string conversions ...");
strs = {"a"; "bc"; "def"};
sv = XLALCreateStringVector(strs{:});
assert(sv.length, 3);
assert(all(strcmp(sv.data, strs)));
strs{1} = "ghijk";
sv.data{1} = strs{1};
strs{end+1} = "lmnopq";
sv = XLALAppendString2Vector(sv, strs{4});
assert(sv.length, 4);
for i = 1:4
  assert(sv.data{i}, strs{i});
endfor
clear sv;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED string conversions");

## check static vector/matrix conversions
disp("checking static vector/matrix conversions ...");
lal.swig_lal_test_struct_vector{1} = lal.swig_lal_test_struct_const;
assert(lal.swig_lal_test_struct_vector{1}.n, lal.swig_lal_test_struct_const.n);
assert(lal.swig_lal_test_struct_vector{1}.i, lal.swig_lal_test_struct_const.i);
assert(lal.swig_lal_test_struct_vector{1}.f, lal.swig_lal_test_struct_const.f);
assert(lal.swig_lal_test_struct_vector{1}.str, lal.swig_lal_test_struct_const.str);
assert(lal.swig_lal_test_struct_vector{1}.vec(:), lal.swig_lal_test_struct_const.vec(:));
lal.swig_lal_test_struct_matrix{1, 1} = lal.swig_lal_test_struct_const;
assert(lal.swig_lal_test_struct_matrix{1, 1}.n, lal.swig_lal_test_struct_const.n);
assert(lal.swig_lal_test_struct_matrix{1, 1}.i, lal.swig_lal_test_struct_const.i);
assert(lal.swig_lal_test_struct_matrix{1, 1}.f, lal.swig_lal_test_struct_const.f);
assert(lal.swig_lal_test_struct_matrix{1, 1}.str, lal.swig_lal_test_struct_const.str);
assert(lal.swig_lal_test_struct_matrix{1, 1}.vec(:), lal.swig_lal_test_struct_const.vec(:));
sts = new_swig_lal_test_struct();
assert(length(sts.vec), 3);
assert(length(sts.evec), 3);
assert(size(sts.mat), [2, 3]);
sts.vec = [3; 2; 1];
assert(sts.vec(:), int32([3; 2; 1]));
sts.mat = [4, 5, 6; 9, 8, 7];
set_nice_error_handlers();
expected_exception = 0;
try
  sts.mat = [1.1, 2.3, 4.5; 6.5, 4.3, 2.1];
  expected_exception = 1;
end_try_catch
set_default_error_handlers();
assert(!expected_exception);
assert(all(sts.mat, [4, 5, 6; 9, 8, 7]));
for i = 1:3
  sts.evec(i) = 2*i + 3;
  assert(sts.evec(i), int32(2*i + 3));
endfor
clear sts;
clear ans;
assert(!any(lal.swig_lal_test_enum_vector));
assert(!any(lal.swig_lal_test_enum_matrix(:)));
assert(length(lal.swig_lal_test_empty_INT4_vector), 0);
assert(!any(lal.swig_lal_test_INT4_vector));
assert(!any(lal.swig_lal_test_INT4_matrix(:)));
assert(!any(lal.swig_lal_test_REAL8_vector));
assert(!any(lal.swig_lal_test_REAL8_matrix(:)));
assert(!any(lal.swig_lal_test_COMPLEX8_vector));
assert(!any(lal.swig_lal_test_COMPLEX8_matrix(:)));
lal.swig_lal_test_INT4_vector(1) = 10;
assert(lal.swig_lal_test_INT4_vector(1), int32(10));
lal.swig_lal_test_INT4_matrix(1, 1) = 11;
assert(lal.swig_lal_test_INT4_matrix(1, 1), int32(11));
lal.swig_lal_test_INT4_vector = lal.swig_lal_test_INT4_const_vector;
assert(lal.swig_lal_test_INT4_vector(:), int32([1; 2; 4]));
assert(lal.swig_lal_test_INT4_const_vector(3), int32(4));
lal.swig_lal_test_INT4_matrix = lal.swig_lal_test_INT4_const_matrix;
assert(lal.swig_lal_test_INT4_matrix(:,:), int32([[1, 2, 4]; [2, 4, 8]]));
assert(lal.swig_lal_test_INT4_const_matrix(2, 3), int32(8));
set_nice_error_handlers();
expected_exception = 0;
try
  lal.swig_lal_test_INT4_const_vector(20);
  expected_exception = 1;
end_try_catch
set_default_error_handlers();
assert(!expected_exception);
lal.swig_lal_test_REAL8_vector(1) = 3.4;
assert(lal.swig_lal_test_REAL8_vector(1), 3.4);
lal.swig_lal_test_REAL8_matrix(1, 1) = 5.6;
assert(lal.swig_lal_test_REAL8_matrix(1, 1), 5.6);
lal.swig_lal_test_COMPLEX8_vector(1) = complex(3.5, 4.75);
assert(lal.swig_lal_test_COMPLEX8_vector(1), complex(single(3.5), single(4.75)));
lal.swig_lal_test_COMPLEX8_matrix(1, 1) = complex(5.5, 6.25);
assert(lal.swig_lal_test_COMPLEX8_matrix(1, 1), complex(single(5.5), single(6.25)));
disp("PASSED static vector/matrix conversions");

## check dynamic vector/matrix conversions
disp("checking dynamic vector/matrix conversions ...");
function check_dynamic_vector_matrix(iv, ivl, rv, rvl, cm, cms1, cms2)
  expected_exception = 0;
  iv.data = zeros(ivl, 1);
  rv.data = zeros(rvl, 1);
  cm.data = zeros(cms1, cms2);
  assert(ivl, 5);
  iv.data = [1; 3; 2; 4; 3];
  assert(iv.data(:), int32([1; 3; 2; 4; 3]));
  iv.data(4) = 7;
  assert(iv.data(4), int32(7));
  assert(rvl, 5);
  rv.data = [1.2; 3.4; 2.6; 4.8; 3.5];
  assert(rv.data(:), [1.2; 3.4; 2.6; 4.8; 3.5]);
  rv.data(rvl) = 7.5;
  assert(rv.data(rvl), 7.5);
  set_nice_error_handlers();
  expected_exception = 0;
  try
    rv.data(rvl + 1) = 99.9;
    expected_exception = 1;
  end_try_catch
  set_default_error_handlers();
  assert(!expected_exception);
  set_nice_error_handlers();
  expected_exception = 0;
  try
    iv.data = rv.data;
    expected_exception = 1;
  end_try_catch
  set_default_error_handlers();
  assert(!expected_exception);
  rv.data = iv.data;
  assert(double(rv.data(:)), double(iv.data(:)));
  assert(cms1, 4);
  assert(cms2, 6);
  for i = 1:cms1
    for j = 1:cms2
      cm.data(i, j) = complex(i / 4.0, j / 2.0);
    endfor
  endfor
  assert(double(cm.data(2, 3)), complex(0.5, 1.5));
  assert(double(cm.data(3, 2)), complex(0.75, 1.0));
  set_nice_error_handlers();
  expected_exception = 0;
  try
    iv.data(0) = cm.data(2, 3);
    expected_exception = 1;
  end_try_catch
  set_default_error_handlers();
  assert(!expected_exception);
  set_nice_error_handlers();
  expected_exception = 0;
  try
    rv.data(0) = cm.data(3, 2);
    expected_exception = 1;
  end_try_catch
  set_default_error_handlers();
  assert(!expected_exception);
endfunction
## check LAL vector and matrix datatypes
iv = XLALCreateINT4Vector(5);
rv = XLALCreateREAL8Vector(5);
cm = XLALCreateCOMPLEX8VectorSequence(4, 6);
check_dynamic_vector_matrix(iv, iv.length, rv, rv.length,
                            cm, cm.length, cm.vectorLength);
clear iv;
clear rv;
clear cm;
clear ans;
rv0 = XLALCreateREAL8Vector(0);
assert(rv0.length, 0);
assert(length(rv0.data), 0);
clear rv0;
clear ans;
rv1 = XLALCreateREAL8Vector(1);
rv1.data(1) = 1;
clear rv1;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED dynamic vector/matrix conversions (LAL)");
## check GSL vectors and matrices
iv = new_gsl_vector_int(5);
rv = new_gsl_vector(5);
cm = new_gsl_matrix_complex_float(4, 6);
check_dynamic_vector_matrix(iv, iv.size, rv, rv.size,
                            cm, cm.size1, cm.size2);
clear iv;
clear rv;
clear cm;
clear ans;
rv1 = new_gsl_vector(1);
rv1.data(1) = 1;
clear rv1;
clear ans;
disp("PASSED dynamic vector/matrix conversions (GSL)");

## check fixed and dynamic arrays typemaps
disp("checking fixed and dynamic arrays typemaps ...");
a1in = double([1.2; 3.5; 7.9]);
a1out = a1in * 2.5;
assert(swig_lal_test_copyin_array1(a1in, 2.5), a1out);
a2in = int32([3,2; 7,6; 12,10]);
a2out = a2in * 15;
assert(swig_lal_test_copyin_array2(a2in, 15), a2out);
a3in = {lal.LIGOTimeGPS(1234.5); lal.LIGOTimeGPS(678.9)};
a3out = {a3in{1} * 3; a3in{2} * 3};
assert(all(cellfun(@(x, y) x == y, lal.swig_lal_test_copyin_array3(a3in, 3), a3out)));
set_nice_error_handlers();
expected_exception = 0;
try
  swig_lal_test_viewin_array1([0,0,0,0], 0);
  expected_exception = 1;
end_try_catch
set_default_error_handlers();
assert(!expected_exception);
set_nice_error_handlers();
expected_exception = 0;
try
  swig_lal_test_viewin_array2([1.2,3.4; 0,0; 0,0], 0);
  expected_exception = 1;
end_try_catch
set_default_error_handlers();
assert(!expected_exception);
clear a3in;
clear a3out;
LALCheckMemoryLeaks();
disp("PASSED fixed and dynamic arrays typemaps")

## check input views of string array structs
disp("checking input views of string array structs ...");
svdat = {"a"; "bc"; "def"};
sv = XLALCreateEmptyStringVector(length(svdat));
sv.data = svdat;
svout = XLALCreateEmptyStringVector(length(svdat));
[svout.data{1:length(svdat)}] = deal("");
assert(swig_lal_test_viewin_LALStringVector(svout, sv));
assert(all(cellfun(@(x, y) strcmp(x, y), svout.data, sv.data)));
[svout.data{1:length(svdat)}] = deal("");
assert(swig_lal_test_viewin_LALStringVector(svout, svdat));
assert(all(cellfun(@(x, y) strcmp(x, y), svout.data, svdat)));
sv.data = svdat;
assert(swig_lal_test_copyinout_LALStringVector(sv));
assert(all(cellfun(@(x, y) strcmp(x, toupper(y)), sv.data, svdat)));
sv.data = svdat;
[retn, sv] = swig_lal_test_copyinout_LALStringVector(sv);
assert(retn);
assert(all(cellfun(@(x, y) strcmp(x, toupper(y)), sv.data, svdat)));
sv = svdat;
[retn, sv] = swig_lal_test_copyinout_LALStringVector(sv);
assert(retn);
assert(all(cellfun(@(x, y) strcmp(x, toupper(y)), sv, svdat)));
clear sv;
clear svout;
clear svdat;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED input views of string array structs");

## check input views of numeric array structs
disp("checking input views of numeric array structs ...");
r4dat = single([1.2; 2.3; 3.4; 4.5; 5.6]);
r8dat = double([3.4; 4.5; 5.6; 6.7; 7.8; 8.9]);
c8dat = complex(r4dat, 8 + r4dat);
c16dat = complex(r8dat, 16 + r8dat);
r4 = XLALCreateREAL4Vector(length(r4dat));
r4.data = r4dat;
r4out = XLALCreateREAL4Vector(length(r4dat));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewin_REAL4Vector(r4out, r4));
assert(r4out.data(:), r4.data(:));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewin_REAL4Vector(r4out, r4dat));
assert(r4out.data(:), r4dat);
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewinout_REAL4Vector(r4out, r4));
assert(2 * r4out.data, r4.data(:));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewinout_REAL4Vector(r4out, r4dat));
assert(2 * r4out.data, r4dat);
r4.data = r4dat;
assert(swig_lal_test_copyinout_REAL4Vector(r4));
assert(r4.data(:), 3 * r4dat);
r4.data = r4dat;
[retn, r4] = swig_lal_test_copyinout_REAL4Vector(r4);
assert(retn);
assert(r4.data(:), 3 * r4dat);
r4 = r4dat;
[retn, r4] = swig_lal_test_copyinout_REAL4Vector(r4);
assert(retn);
assert(r4, 3 * r4dat);
clear r4;
clear r4out;
clear r4dat;
clear ans;
LALCheckMemoryLeaks();
r8 = XLALCreateREAL8Vector(length(r8dat));
r8.data = r8dat;
r8out = XLALCreateREAL8Vector(length(r8dat));
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewin_REAL8Vector(r8out, r8));
assert(r8out.data(:), r8.data(:));
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewin_REAL8Vector(r8out, r8dat));
assert(r8out.data(:), r8dat);
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewinout_REAL8Vector(r8out, r8));
assert(2 * r8out.data, r8.data(:));
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewinout_REAL8Vector(r8out, r8dat));
assert(2 * r8out.data, r8dat);
r8.data = r8dat;
assert(swig_lal_test_copyinout_REAL8Vector(r8));
assert(r8.data(:), 3 * r8dat);
r8.data = r8dat;
[retn, r8] = swig_lal_test_copyinout_REAL8Vector(r8);
assert(retn);
assert(r8.data(:), 3 * r8dat);
r8 = r8dat;
[retn, r8] = swig_lal_test_copyinout_REAL8Vector(r8);
assert(retn);
assert(r8, 3 * r8dat);
clear r8;
clear r8out;
clear r8dat;
clear ans;
LALCheckMemoryLeaks();
c8 = XLALCreateCOMPLEX8Vector(length(c8dat));
c8.data = c8dat;
c8out = XLALCreateCOMPLEX8Vector(length(c8dat));
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewin_COMPLEX8Vector(c8out, c8));
assert(c8out.data(:), c8.data(:));
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewin_COMPLEX8Vector(c8out, c8dat));
assert(c8out.data(:), c8dat);
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewinout_COMPLEX8Vector(c8out, c8));
assert(2 * c8out.data, c8.data(:));
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewinout_COMPLEX8Vector(c8out, c8dat));
assert(2 * c8out.data, c8dat);
c8.data = c8dat;
assert(swig_lal_test_copyinout_COMPLEX8Vector(c8));
assert(c8.data(:), 3 * c8dat);
c8.data = c8dat;
[retn, c8] = swig_lal_test_copyinout_COMPLEX8Vector(c8);
assert(retn);
assert(c8.data(:), 3 * c8dat);
c8 = c8dat;
[retn, c8] = swig_lal_test_copyinout_COMPLEX8Vector(c8);
assert(retn);
assert(c8, 3 * c8dat);
clear c8;
clear c8out;
clear c8dat;
clear ans;
LALCheckMemoryLeaks();
c16 = XLALCreateCOMPLEX16Vector(length(c16dat));
c16.data = c16dat;
c16out = XLALCreateCOMPLEX16Vector(length(c16dat));
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewin_COMPLEX16Vector(c16out, c16));
assert(c16out.data(:), c16.data(:));
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewin_COMPLEX16Vector(c16out, c16dat));
assert(c16out.data(:), c16dat);
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewinout_COMPLEX16Vector(c16out, c16));
assert(2 * c16out.data, c16.data(:));
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewinout_COMPLEX16Vector(c16out, c16dat));
assert(2 * c16out.data, c16dat);
c16.data = c16dat;
assert(swig_lal_test_copyinout_COMPLEX16Vector(c16));
assert(c16.data(:), 3 * c16dat);
c16.data = c16dat;
[retn, c16] = swig_lal_test_copyinout_COMPLEX16Vector(c16);
assert(retn);
assert(c16.data(:), 3 * c16dat);
c16 = c16dat;
[retn, c16] = swig_lal_test_copyinout_COMPLEX16Vector(c16);
assert(retn);
assert(c16, 3 * c16dat);
clear c16;
clear c16out;
clear c16dat;
clear ans;
LALCheckMemoryLeaks();
r4dat = single([1.2, 2.3, 3.4; 4.5, 5.6, 6.7]);
r8dat = double([3.4, 4.5; 5.6, 6.7; 7.8, 8.9]);
c8dat = complex(r4dat, 8 + r4dat);
c16dat = complex(r8dat, 16 + r8dat);
r4 = XLALCreateREAL4VectorSequence(size(r4dat,1), size(r4dat,2));
r4.data = r4dat;
r4out = XLALCreateREAL4VectorSequence(size(r4dat,1), size(r4dat,2));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewin_REAL4VectorSequence(r4out, r4));
assert(r4out.data(:,:), r4.data(:,:));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewin_REAL4VectorSequence(r4out, r4dat));
assert(r4out.data(:,:), r4dat);
r4.data = r4dat;
assert(swig_lal_test_copyinout_REAL4VectorSequence(r4));
assert(r4.data(:,:), 3 * r4dat);
r4.data = r4dat;
[retn, r4] = swig_lal_test_copyinout_REAL4VectorSequence(r4);
assert(retn);
assert(r4.data(:,:), 3 * r4dat);
r4 = r4dat;
[retn, r4] = swig_lal_test_copyinout_REAL4VectorSequence(r4);
assert(retn);
assert(r4, 3 * r4dat);
clear r4;
clear r4out;
clear r4dat;
clear ans;
LALCheckMemoryLeaks();
r8 = XLALCreateREAL8VectorSequence(size(r8dat,1), size(r8dat,2));
r8.data = r8dat;
r8out = XLALCreateREAL8VectorSequence(size(r8dat,1), size(r8dat,2));
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewin_REAL8VectorSequence(r8out, r8));
assert(r8out.data(:,:), r8.data(:,:));
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewin_REAL8VectorSequence(r8out, r8dat));
assert(r8out.data(:,:), r8dat);
r8.data = r8dat;
assert(swig_lal_test_copyinout_REAL8VectorSequence(r8));
assert(r8.data(:,:), 3 * r8dat);
r8.data = r8dat;
[retn, r8] = swig_lal_test_copyinout_REAL8VectorSequence(r8);
assert(retn);
assert(r8.data(:,:), 3 * r8dat);
r8 = r8dat;
[retn, r8] = swig_lal_test_copyinout_REAL8VectorSequence(r8);
assert(retn);
assert(r8, 3 * r8dat);
clear r8;
clear r8out;
clear r8dat;
clear ans;
LALCheckMemoryLeaks();
c8 = XLALCreateCOMPLEX8VectorSequence(size(c8dat,1), size(c8dat,2));
c8.data = c8dat;
c8out = XLALCreateCOMPLEX8VectorSequence(size(c8dat,1), size(c8dat,2));
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewin_COMPLEX8VectorSequence(c8out, c8));
assert(c8out.data(:,:), c8.data(:,:));
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewin_COMPLEX8VectorSequence(c8out, c8dat));
assert(c8out.data(:,:), c8dat);
c8.data = c8dat;
assert(swig_lal_test_copyinout_COMPLEX8VectorSequence(c8));
assert(c8.data(:,:), 3 * c8dat);
c8.data = c8dat;
[retn, c8] = swig_lal_test_copyinout_COMPLEX8VectorSequence(c8);
assert(retn);
assert(c8.data(:,:), 3 * c8dat);
c8 = c8dat;
[retn, c8] = swig_lal_test_copyinout_COMPLEX8VectorSequence(c8);
assert(retn);
assert(c8, 3 * c8dat);
clear c8;
clear c8out;
clear c8dat;
clear ans;
LALCheckMemoryLeaks();
c16 = XLALCreateCOMPLEX16VectorSequence(size(c16dat,1), size(c16dat,2));
c16.data = c16dat;
c16out = XLALCreateCOMPLEX16VectorSequence(size(c16dat,1), size(c16dat,2));
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewin_COMPLEX16VectorSequence(c16out, c16));
assert(c16out.data(:,:), c16.data(:,:));
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewin_COMPLEX16VectorSequence(c16out, c16dat));
assert(c16out.data(:,:), c16dat);
c16.data = c16dat;
assert(swig_lal_test_copyinout_COMPLEX16VectorSequence(c16));
assert(c16.data(:,:), 3 * c16dat);
c16.data = c16dat;
[retn, c16] = swig_lal_test_copyinout_COMPLEX16VectorSequence(c16);
assert(retn);
assert(c16.data(:,:), 3 * c16dat);
c16 = c16dat;
[retn, c16] = swig_lal_test_copyinout_COMPLEX16VectorSequence(c16);
assert(retn);
assert(c16, 3 * c16dat);
clear c16;
clear c16out;
clear c16dat;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED input views of numeric array structs (LAL)");
vfdat = single([1.2; 2.3; 3.4; 4.5; 5.6]);
vddat = double([3.4; 4.5; 5.6; 6.7; 7.8; 8.9]);
vcfdat = complex(vfdat, 8 + vfdat);
vcddat = complex(vddat, 16 + vddat);
vf = new_gsl_vector_float(length(vfdat));
vf.data = vfdat;
vfout = new_gsl_vector_float(length(vfdat));
vfout.data = zeros(size(vfdat));
assert(swig_lal_test_viewin_gsl_vector_float(vfout, vf));
assert(vfout.data(:), vf.data(:));
vfout.data = zeros(size(vfdat));
assert(swig_lal_test_viewin_gsl_vector_float(vfout, vfdat));
assert(vfout.data(:), vfdat);
vfout.data = zeros(size(vfdat));
assert(swig_lal_test_viewinout_gsl_vector_float(vfout, vf));
assert(2 * vfout.data, vf.data(:));
vfout.data = zeros(size(vfdat));
assert(swig_lal_test_viewinout_gsl_vector_float(vfout, vfdat));
assert(2 * vfout.data, vfdat);
vf.data = vfdat;
assert(swig_lal_test_copyinout_gsl_vector_float(vf));
assert(vf.data(:), 3 * vfdat);
vf.data = vfdat;
[retn, vf] = swig_lal_test_copyinout_gsl_vector_float(vf);
assert(retn);
assert(vf.data(:), 3 * vfdat);
vf = vfdat;
[retn, vf] = swig_lal_test_copyinout_gsl_vector_float(vf);
assert(retn);
assert(vf, 3 * vfdat);
clear vf;
clear vfout;
clear vfdat;
clear ans;
LALCheckMemoryLeaks();
vd = new_gsl_vector(length(vddat));
vd.data = vddat;
vdout = new_gsl_vector(length(vddat));
vdout.data = zeros(size(vddat));
assert(swig_lal_test_viewin_gsl_vector(vdout, vd));
assert(vdout.data(:), vd.data(:));
vdout.data = zeros(size(vddat));
assert(swig_lal_test_viewin_gsl_vector(vdout, vddat));
assert(vdout.data(:), vddat);
vdout.data = zeros(size(vddat));
assert(swig_lal_test_viewinout_gsl_vector(vdout, vd));
assert(2 * vdout.data, vd.data(:));
vdout.data = zeros(size(vddat));
assert(swig_lal_test_viewinout_gsl_vector(vdout, vddat));
assert(2 * vdout.data, vddat);
vd.data = vddat;
assert(swig_lal_test_copyinout_gsl_vector(vd));
assert(vd.data(:), 3 * vddat);
vd.data = vddat;
[retn, vd] = swig_lal_test_copyinout_gsl_vector(vd);
assert(retn);
assert(vd.data(:), 3 * vddat);
vd = vddat;
[retn, vd] = swig_lal_test_copyinout_gsl_vector(vd);
assert(retn);
assert(vd, 3 * vddat);
clear vd;
clear vdout;
clear vddat;
clear ans;
LALCheckMemoryLeaks();
vcf = new_gsl_vector_complex_float(length(vcfdat));
vcf.data = vcfdat;
vcfout = new_gsl_vector_complex_float(length(vcfdat));
vcfout.data = zeros(size(vcfdat));
assert(swig_lal_test_viewin_gsl_vector_complex_float(vcfout, vcf));
assert(vcfout.data(:), vcf.data(:));
vcfout.data = zeros(size(vcfdat));
assert(swig_lal_test_viewin_gsl_vector_complex_float(vcfout, vcfdat));
assert(vcfout.data(:), vcfdat);
vcfout.data = zeros(size(vcfdat));
assert(swig_lal_test_viewinout_gsl_vector_complex_float(vcfout, vcf));
assert(2 * vcfout.data, vcf.data(:));
vcfout.data = zeros(size(vcfdat));
assert(swig_lal_test_viewinout_gsl_vector_complex_float(vcfout, vcfdat));
assert(2 * vcfout.data, vcfdat);
vcf.data = vcfdat;
assert(swig_lal_test_copyinout_gsl_vector_complex_float(vcf));
assert(vcf.data(:), 3 * vcfdat);
vcf.data = vcfdat;
[retn, vcf] = swig_lal_test_copyinout_gsl_vector_complex_float(vcf);
assert(retn);
assert(vcf.data(:), 3 * vcfdat);
vcf = vcfdat;
[retn, vcf] = swig_lal_test_copyinout_gsl_vector_complex_float(vcf);
assert(retn);
assert(vcf, 3 * vcfdat);
clear vcf;
clear vcfout;
clear vcfdat;
clear ans;
LALCheckMemoryLeaks();
vcd = new_gsl_vector_complex(length(vcddat));
vcd.data = vcddat;
vcdout = new_gsl_vector_complex(length(vcddat));
vcdout.data = zeros(size(vcddat));
assert(swig_lal_test_viewin_gsl_vector_complex(vcdout, vcd));
assert(vcdout.data(:), vcd.data(:));
vcdout.data = zeros(size(vcddat));
assert(swig_lal_test_viewin_gsl_vector_complex(vcdout, vcddat));
assert(vcdout.data(:), vcddat);
vcdout.data = zeros(size(vcddat));
assert(swig_lal_test_viewinout_gsl_vector_complex(vcdout, vcd));
assert(2 * vcdout.data, vcd.data(:));
vcdout.data = zeros(size(vcddat));
assert(swig_lal_test_viewinout_gsl_vector_complex(vcdout, vcddat));
assert(2 * vcdout.data, vcddat);
vcd.data = vcddat;
assert(swig_lal_test_copyinout_gsl_vector_complex(vcd));
assert(vcd.data(:), 3 * vcddat);
vcd.data = vcddat;
[retn, vcd] = swig_lal_test_copyinout_gsl_vector_complex(vcd);
assert(retn);
assert(vcd.data(:), 3 * vcddat);
vcd = vcddat;
[retn, vcd] = swig_lal_test_copyinout_gsl_vector_complex(vcd);
assert(retn);
assert(vcd, 3 * vcddat);
clear vcd;
clear vcdout;
clear vcddat;
clear ans;
LALCheckMemoryLeaks();
mfdat = single([1.2, 2.3, 3.4; 4.5, 5.6, 6.7]);
mddat = double([3.4, 4.5; 5.6, 6.7; 7.8, 8.9]);
mcfdat = complex(mfdat, 8 + mfdat);
mcddat = complex(mddat, 16 + mddat);
mf = new_gsl_matrix_float(size(mfdat,1), size(mfdat,2));
mf.data = mfdat;
mfout = new_gsl_matrix_float(size(mfdat,1), size(mfdat,2));
mfout.data = zeros(size(mfdat));
assert(swig_lal_test_viewin_gsl_matrix_float(mfout, mf));
assert(mfout.data(:,:), mf.data(:,:));
mfout.data = zeros(size(mfdat));
assert(swig_lal_test_viewin_gsl_matrix_float(mfout, mfdat));
assert(mfout.data(:,:), mfdat);
mf.data = mfdat;
assert(swig_lal_test_copyinout_gsl_matrix_float(mf));
assert(mf.data(:,:), 3 * mfdat);
mf.data = mfdat;
[retn, mf] = swig_lal_test_copyinout_gsl_matrix_float(mf);
assert(retn);
assert(mf.data(:,:), 3 * mfdat);
mf = mfdat;
[retn, mf] = swig_lal_test_copyinout_gsl_matrix_float(mf);
assert(retn);
assert(mf, 3 * mfdat);
clear mf;
clear mfout;
clear mfdat;
clear ans;
LALCheckMemoryLeaks();
md = new_gsl_matrix(size(mddat,1), size(mddat,2));
md.data = mddat;
mdout = new_gsl_matrix(size(mddat,1), size(mddat,2));
mdout.data = zeros(size(mddat));
assert(swig_lal_test_viewin_gsl_matrix(mdout, md));
assert(mdout.data(:,:), md.data(:,:));
mdout.data = zeros(size(mddat));
assert(swig_lal_test_viewin_gsl_matrix(mdout, mddat));
assert(mdout.data(:,:), mddat);
md.data = mddat;
assert(swig_lal_test_copyinout_gsl_matrix(md));
assert(md.data(:,:), 3 * mddat);
md.data = mddat;
[retn, md] = swig_lal_test_copyinout_gsl_matrix(md);
assert(retn);
assert(md.data(:,:), 3 * mddat);
md = mddat;
[retn, md] = swig_lal_test_copyinout_gsl_matrix(md);
assert(retn);
assert(md, 3 * mddat);
clear md;
clear mdout;
clear mddat;
clear ans;
LALCheckMemoryLeaks();
mcf = new_gsl_matrix_complex_float(size(mcfdat,1), size(mcfdat,2));
mcf.data = mcfdat;
mcfout = new_gsl_matrix_complex_float(size(mcfdat,1), size(mcfdat,2));
mcfout.data = zeros(size(mcfdat));
assert(swig_lal_test_viewin_gsl_matrix_complex_float(mcfout, mcf));
assert(mcfout.data(:,:), mcf.data(:,:));
mcfout.data = zeros(size(mcfdat));
assert(swig_lal_test_viewin_gsl_matrix_complex_float(mcfout, mcfdat));
assert(mcfout.data(:,:), mcfdat);
mcf.data = mcfdat;
assert(swig_lal_test_copyinout_gsl_matrix_complex_float(mcf));
assert(mcf.data(:,:), 3 * mcfdat);
mcf.data = mcfdat;
[retn, mcf] = swig_lal_test_copyinout_gsl_matrix_complex_float(mcf);
assert(retn);
assert(mcf.data(:,:), 3 * mcfdat);
mcf = mcfdat;
[retn, mcf] = swig_lal_test_copyinout_gsl_matrix_complex_float(mcf);
assert(retn);
assert(mcf, 3 * mcfdat);
clear mcf;
clear mcfout;
clear mcfdat;
clear ans;
LALCheckMemoryLeaks();
mcd = new_gsl_matrix_complex(size(mcddat,1), size(mcddat,2));
mcd.data = mcddat;
mcdout = new_gsl_matrix_complex(size(mcddat,1), size(mcddat,2));
mcdout.data = zeros(size(mcddat));
assert(swig_lal_test_viewin_gsl_matrix_complex(mcdout, mcd));
assert(mcdout.data(:,:), mcd.data(:,:));
mcdout.data = zeros(size(mcddat));
assert(swig_lal_test_viewin_gsl_matrix_complex(mcdout, mcddat));
assert(mcdout.data(:,:), mcddat);
mcd.data = mcddat;
assert(swig_lal_test_copyinout_gsl_matrix_complex(mcd));
assert(mcd.data(:,:), 3 * mcddat);
mcd.data = mcddat;
[retn, mcd] = swig_lal_test_copyinout_gsl_matrix_complex(mcd);
assert(retn);
assert(mcd.data(:,:), 3 * mcddat);
mcd = mcddat;
[retn, mcd] = swig_lal_test_copyinout_gsl_matrix_complex(mcd);
assert(retn);
assert(mcd, 3 * mcddat);
clear mcd;
clear mcdout;
clear mcddat;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED input views of numeric array structs (GSL)");
function check_input_view_type_safety(f, a, b, expect_exception)
  expected_exception = 0;
  if expect_exception
    set_nice_error_handlers();
    expected_exception = 0;
    try
      f(a, b);
      expected_exception = 1;
    end_try_catch
    set_default_error_handlers();
    assert(!expected_exception);
    set_nice_error_handlers();
    expected_exception = 0;
    try
      f(b, a);
      expected_exception = 1;
    end_try_catch
    set_default_error_handlers();
    assert(!expected_exception);
  else
    f(a, b);
    f(b, a);
  endif
endfunction
r4 = single(zeros(10, 1));
r8 = zeros(10, 1);
c8 = complex(r4);
c16 = complex(r8);
check_input_view_type_safety(@swig_lal_test_viewinout_REAL4Vector, r4, r4, false);
check_input_view_type_safety(@swig_lal_test_viewinout_REAL4Vector, r4, r8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_REAL4Vector, r4, c8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_REAL4Vector, r4, c16, true);
check_input_view_type_safety(@swig_lal_test_viewinout_REAL8Vector, r8, r4, true);
check_input_view_type_safety(@swig_lal_test_viewinout_REAL8Vector, r8, r8, false);
check_input_view_type_safety(@swig_lal_test_viewinout_REAL8Vector, r8, c8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_REAL8Vector, r8, c16, true);
check_input_view_type_safety(@swig_lal_test_viewinout_COMPLEX8Vector, c8, r4, true);
check_input_view_type_safety(@swig_lal_test_viewinout_COMPLEX8Vector, c8, r8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_COMPLEX8Vector, c8, c8, false);
check_input_view_type_safety(@swig_lal_test_viewinout_COMPLEX8Vector, c8, c16, true);
check_input_view_type_safety(@swig_lal_test_viewinout_COMPLEX16Vector, c16, r4, true);
check_input_view_type_safety(@swig_lal_test_viewinout_COMPLEX16Vector, c16, r8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_COMPLEX16Vector, c16, c8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_COMPLEX16Vector, c16, c16, false);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_float, r4, r4, false);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_float, r4, r8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_float, r4, c8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_float, r4, c16, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector, r8, r4, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector, r8, r8, false);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector, r8, c8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector, r8, c16, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_complex_float, c8, r4, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_complex_float, c8, r8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_complex_float, c8, c8, false);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_complex_float, c8, c16, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_complex, c16, r4, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_complex, c16, r8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_complex, c16, c8, true);
check_input_view_type_safety(@swig_lal_test_viewinout_gsl_vector_complex, c16, c16, false);
clear r4;
clear r8;
clear c8;
clear c16;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED input views of numeric array structs (type safety)");

## check FFT functions with input views
disp("check FFT functions with input views ...");
r4in = single(0:31)';
r8in = double(0:63)';
c8in = complex(8 + r4in, r4in);
c16in = complex(16 + r8in, r8in);
c8inv = XLALCreateCOMPLEX8Vector(length(c8in));
c8inv.data = c8in;
c8outv = XLALCreateCOMPLEX8Vector(length(c8in));
plan = XLALCreateForwardCOMPLEX8FFTPlan(length(c8in), 0);
XLALCOMPLEX8VectorFFT(c8outv, c8inv, plan);
c8out = complex(single(zeros(size(c8outv.data))));
XLALCOMPLEX8VectorFFT(c8out, c8in, plan);
assert(c8out, c8outv.data(:));
clear c8inv;
clear c8outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
c16inv = XLALCreateCOMPLEX16Vector(length(c16in));
c16inv.data = c16in;
c16outv = XLALCreateCOMPLEX16Vector(length(c16in));
plan = XLALCreateForwardCOMPLEX16FFTPlan(length(c16in), 0);
XLALCOMPLEX16VectorFFT(c16outv, c16inv, plan);
c16out = complex(double(zeros(size(c16outv.data))));
XLALCOMPLEX16VectorFFT(c16out, c16in, plan);
assert(c16out, c16outv.data(:));
clear c16inv;
clear c16outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
r4inv = XLALCreateREAL4Vector(length(r4in));
r4inv.data = r4in;
c8outv = XLALCreateCOMPLEX8Vector(length(r4in)/2 + 1);
plan = XLALCreateForwardREAL4FFTPlan(length(r4in), 0);
XLALREAL4ForwardFFT(c8outv, r4inv, plan);
c8out = complex(single(zeros(size(c8outv.data))));
XLALREAL4ForwardFFT(c8out, r4in, plan);
assert(c8out, c8outv.data(:));
clear r4inv;
clear c8outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
c8inv = XLALCreateCOMPLEX8Vector(length(c8in));
c8inv.data = c8in;
r4outv = XLALCreateREAL4Vector((length(c8in)-1)*2);
plan = XLALCreateReverseREAL4FFTPlan((length(c8in)-1)*2, 0);
XLALREAL4ReverseFFT(r4outv, c8inv, plan);
r4out = single(zeros(size(r4outv.data)));
XLALREAL4ReverseFFT(r4out, c8in, plan);
assert(r4out, r4outv.data(:));
clear c8inv;
clear r4outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
r8inv = XLALCreateREAL8Vector(length(r8in));
r8inv.data = r8in;
c16outv = XLALCreateCOMPLEX16Vector(length(r8in)/2 + 1);
plan = XLALCreateForwardREAL8FFTPlan(length(r8in), 0);
XLALREAL8ForwardFFT(c16outv, r8inv, plan);
c16out = complex(double(zeros(size(c16outv.data))));
XLALREAL8ForwardFFT(c16out, r8in, plan);
assert(c16out, c16outv.data(:));
clear r8inv;
clear c16outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
c16inv = XLALCreateCOMPLEX16Vector(length(c16in));
c16inv.data = c16in;
r8outv = XLALCreateREAL8Vector((length(c16in)-1)*2);
plan = XLALCreateReverseREAL8FFTPlan((length(c16in)-1)*2, 0);
XLALREAL8ReverseFFT(r8outv, c16inv, plan);
r8out = double(zeros(size(r8outv.data)));
XLALREAL8ReverseFFT(r8out, c16in, plan);
assert(r8out, r8outv.data(:));
clear c16inv;
clear r8outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
r4inv = XLALCreateREAL4Vector(length(r4in));
r4inv.data = r4in;
r4outv = XLALCreateREAL4Vector(length(r4in));
plan = XLALCreateForwardREAL4FFTPlan(length(r4in), 0);
XLALREAL4VectorFFT(r4outv, r4inv, plan);
r4out = single(zeros(size(r4outv.data)));
XLALREAL4VectorFFT(r4out, r4in, plan);
assert(r4out, r4outv.data(:));
clear r4inv;
clear r4outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
r8inv = XLALCreateREAL8Vector(length(r8in));
r8inv.data = r8in;
r8outv = XLALCreateREAL8Vector(length(r8in));
plan = XLALCreateForwardREAL8FFTPlan(length(r8in), 0);
XLALREAL8VectorFFT(r8outv, r8inv, plan);
r8out = double(zeros(size(r8outv.data)));
XLALREAL8VectorFFT(r8out, r8in, plan);
assert(r8out, r8outv.data(:));
clear r8inv;
clear r8outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
r4inv = XLALCreateREAL4Vector(length(r4in));
r4inv.data = r4in;
r4outv = XLALCreateREAL4Vector(length(r4in)/2 + 1);
plan = XLALCreateForwardREAL4FFTPlan(length(r4in), 0);
XLALREAL4PowerSpectrum(r4outv, r4inv, plan);
r4out = single(zeros(size(r4outv.data)));
XLALREAL4PowerSpectrum(r4out, r4in, plan);
assert(r4out, r4outv.data(:));
clear r4inv;
clear r4outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
r8inv = XLALCreateREAL8Vector(length(r8in));
r8inv.data = r8in;
r8outv = XLALCreateREAL8Vector(length(r8in)/2 + 1);
plan = XLALCreateForwardREAL8FFTPlan(length(r8in), 0);
XLALREAL8PowerSpectrum(r8outv, r8inv, plan);
r8out = double(zeros(size(r8outv.data)));
XLALREAL8PowerSpectrum(r8out, r8in, plan);
assert(r8out, r8outv.data(:));
clear r8inv;
clear r8outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED FFT functions with input views ...");

## check dynamic array of pointers access
disp("checking dynamic array of pointers access ...");
ap = swig_lal_test_Create_arrayofptrs(3);
assert(ap.length, 3);
for i = 1:ap.length
  assert(ap.data{i}.length, 6);
  for j = 1:ap.data{i}.length
    assert(ap.data{i}.data(j), int32(42*ap.length*(i-1) + (j-1)));
  endfor
endfor
clear ap;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED dynamic array of pointers access");

## check typemaps for strings and double pointers
disp("checking typemaps for strings and double pointers ...");
sts = new_swig_lal_test_struct();
[ptr_ptr, ptr_null_ptr, null_ptr_ptr] = swig_lal_test_typemaps_string_ptrptr("abcde", "", [], sts, 0, []);
assert(swig_this(ptr_ptr), swig_this(sts));
assert(swig_this(ptr_null_ptr) != 0);
assert(swig_this(null_ptr_ptr) == 0);
clear sts;
clear ptr_ptr;
clear ptr_null_ptr;
clear null_ptr_ptr;
LALCheckMemoryLeaks();
ptr_ptr = 0;
for i = 1:9
  ptr_ptr = swig_lal_test_typemaps_ptrptr(ptr_ptr);
  assert(swig_this(ptr_ptr) != 0);
  assert(ptr_ptr.n, i);
endfor
clear ptr_ptr;
ptr_ptr_list = {0};
for i = 1:9
  ptr_ptr_list{end+1} = swig_lal_test_typemaps_ptrptr(ptr_ptr_list{end});
  assert(swig_this(ptr_ptr_list{end}) != 0);
  assert(ptr_ptr_list{end}.n, i);
endfor
while length(ptr_ptr_list) > 0
  assert(swig_this(ptr_ptr_list{end}) != 0);
  assert(ptr_ptr_list{end}.n, i);
  ptr_ptr_list = ptr_ptr_list(2:end);
endwhile
clear ptr_ptr_list;
LALCheckMemoryLeaks();
disp("PASSED typemaps for strings and double pointers");

## check 'tm' struct conversions
disp("checking 'tm' struct conversions ...");
gps0 = 989168284;
utc0 = [2011, 5, 11, 16, 57, 49];
assert(XLALGPSToUTC(gps0), utc0);
assert(XLALUTCToGPS(utc0), gps0);
for i = 0:9
  gps = gps0 + i * 86400;
  utc = utc0;
  utc(3) = utc(3) + i;
  assert(XLALGPSToUTC(gps), utc);
  assert(XLALUTCToGPS(utc), gps);
endfor
LALCheckMemoryLeaks();
disp("PASSED 'tm' struct conversions");

## check LIGOTimeGPS operations
disp("checking LIGOTimeGPS operations ...");
t0 = new_LIGOTimeGPS();
assert(swig_type(LIGOTimeGPS(t0)), "LIGOTimeGPS");
assert(t0 == 0);
assert(swig_type(t0), "LIGOTimeGPS");
t1 = new_LIGOTimeGPS(10.5);
t2 = new_LIGOTimeGPS(10, 500000000);
assert_LIGOTimeGPS(t1, t2);
assert(swig_type(t1), "LIGOTimeGPS");
t3 = +t1;
t3 = -t2;
assert_LIGOTimeGPS(t1, t2);
assert(t1 >= t2 && t2 >= t1);
assert_LIGOTimeGPS(t1 + 3.5, 14);
assert(swig_type(t1 + 3.5), "LIGOTimeGPS");
assert_LIGOTimeGPS(3.5 + t1, 14);
assert(swig_type(3.5 + t1), "LIGOTimeGPS");
t2 -= 5.5;
assert_LIGOTimeGPS(t2, 5);
assert(swig_type(t2), "LIGOTimeGPS");
assert(t2 + 5.5 >= t1);
assert(t2 + 3 != t2);
assert_LIGOTimeGPS(t2 - 5, t0);
assert(swig_type(t2 - 5), "LIGOTimeGPS");
assert_LIGOTimeGPS(t1 * 3, 31.5);
assert(swig_type(t1 * 3), "LIGOTimeGPS");
assert_LIGOTimeGPS(3 * t1, 31.5);
assert(swig_type(3 * t1), "LIGOTimeGPS");
assert_LIGOTimeGPS(t2 / 2.5, 2);
assert(swig_type(t2 / 2.5), "LIGOTimeGPS");
assert_LIGOTimeGPS(21 / t1, 2);
assert(swig_type(21 / t1), "LIGOTimeGPS");
assert_LIGOTimeGPS(t1 + t2, 15.5);
assert(swig_type(t1 + t2), "LIGOTimeGPS");
assert_LIGOTimeGPS(t1 - t2, 5.5);
assert(swig_type(t1 - t2), "LIGOTimeGPS");
assert_LIGOTimeGPS(t1 * t2, 52.5);
assert(swig_type(t1 * t2), "LIGOTimeGPS");
assert_LIGOTimeGPS(t2 * t1, 52.5);
assert(swig_type(t2 * t1), "LIGOTimeGPS");
assert_LIGOTimeGPS(t1 / t2, 2.1);
assert(swig_type(t1 / t2), "LIGOTimeGPS");
assert(t1 > t2);
assert(t2 < t1 && t1 >= t2 && t2 <= t1);
assert_LIGOTimeGPS(LIGOTimeGPS(333333333,333333333), LIGOTimeGPS(1000000000) / 3);
assert_LIGOTimeGPS(LIGOTimeGPS(666666666,666666667), LIGOTimeGPS(2000000000) / 3);
assert_LIGOTimeGPS(LIGOTimeGPS("-62997760.825036067"), LIGOTimeGPS("-47044285.062262587") - LIGOTimeGPS("15953475.76277348"));
assert_LIGOTimeGPS(LIGOTimeGPS("-6542354.389038577"), LIGOTimeGPS("-914984.929117316") * 7.1502318572066237);
assert_LIGOTimeGPS(LIGOTimeGPS("-6542354.389038577"), 7.1502318572066237 * LIGOTimeGPS("-914984.929117316"));
assert_LIGOTimeGPS(LIGOTimeGPS("-127965.770535834"), LIGOTimeGPS("-914984.929117316") / 7.1502318572066237);
t1 += 812345667.75;
assert(t1.__str__(), "812345678.25");
assert_LIGOTimeGPS(new_LIGOTimeGPS(t1.__str__()), t1);
assert(LIGOTimeGPS(1100000000).asutcstr(), "Fri, 14 Nov 2014 11:33:04 +0000");   # lalapps_tconvert -u -R
assert(LIGOTimeGPS(1100000000, 100).asutcstr(), "Fri, 14 Nov 2014 11:33:04.0000001 +0000");
assert(LIGOTimeGPS(0, 0).asutcstr(), "Sun, 06 Jan 1980 00:00:00 +0000");
assert(LIGOTimeGPS(-1, 0).asutcstr(), "Sat, 05 Jan 1980 23:59:59 +0000");
assert(LIGOTimeGPS(0, -1).asutcstr(), "Sat, 05 Jan 1980 23:59:59.999999999 +0000");
assert(t1.ns(), int64(812345678250000000));
t4struct = new_swig_lal_test_gps;
t4struct.t = 1234.5;
assert_LIGOTimeGPS(t4struct.t, 1234.5);
t5 = LIGOTimeGPS("1000");
assert_LIGOTimeGPS(t5, 1000);
disp("*** below should be error messages from LIGOTimeGPS constructor ***");
set_nice_error_handlers();
expected_exception = 0;
try
  t5 = LIGOTimeGPS("abc1000");
  expected_exception = 1;
end_try_catch
set_default_error_handlers();
assert(!expected_exception);
set_nice_error_handlers();
expected_exception = 0;
try
  t5 = LIGOTimeGPS("1000abc");
  expected_exception = 1;
end_try_catch
set_default_error_handlers();
assert(!expected_exception);
disp("*** above should be error messages from LIGOTimeGPS constructor ***");
assert(swig_lal_test_noptrgps(LIGOTimeGPS(1234.5)) == swig_lal_test_noptrgps(1234.5))
disp("*** below should be error messages from LIGOTimeGPS constructor ***");
set_nice_error_handlers();
expected_exception = 0;
try
  LIGOTimeGPS([]);
  expected_exception = 1;
end_try_catch
set_default_error_handlers();
assert(!expected_exception);
disp("*** above should be error messages from LIGOTimeGPS constructor ***");
clear t0;
clear t1;
clear t2;
clear t3;
clear t4struct;
clear t5;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED LIGOTimeGPS operations");

## check LALUnit operations
disp("checking LALUnit operations ...");
u1 = LALUnit("kg m s^-2");
assert(swig_type(LALUnit(u1)), "LALUnit");
assert_LALUnit(u1, lal.NewtonUnit);
assert(swig_type(u1), "LALUnit");
assert(u1.__str__(), "m kg s^-2");
u2 = lal.MeterUnit * lal.KiloGramUnit / lal.SecondUnit ^ 2;
assert_LALUnit(u1, u2);
assert(swig_type(u2), "LALUnit");
u2 = lal.MeterUnit^[1,2] * lal.KiloGramUnit^[1,2] * lal.SecondUnit ^ -1;
assert_LALUnit(u1^[1,2], u2);
assert(swig_type(u2), "LALUnit");
set_nice_error_handlers();
expected_exception = 0;
try
  lal.SecondUnit ^ [1,0];
  expected_exception = 1;
end_try_catch
set_default_error_handlers();
assert(!expected_exception);
u1 *= lal.MeterUnit;
assert_LALUnit(u1, lal.JouleUnit);
assert(swig_type(u1), "LALUnit");
assert(u1.__repr__(), "m^2 kg s^-2");
u1 /= lal.SecondUnit;
assert_LALUnit(u1, lal.WattUnit);
assert(swig_type(u1), "LALUnit");
assert_LALUnit(u1, "m^2 kg s^-3");
u1 *= 1000;
assert_LALUnit(u1, lal.KiloUnit * lal.WattUnit);
assert_LALUnit(u1, 1000 * lal.WattUnit);
assert_LALUnit(u1, lal.WattUnit * 1000);
assert_LALUnit(u1, lal.MegaUnit / 1000 * lal.WattUnit);
assert_LALUnit(u1.__int__(), 1000);
u1 /= 10000;
assert_LALUnit(u1, 100 * lal.MilliUnit * lal.WattUnit);
set_nice_error_handlers();
expected_exception = 0;
try
  u1 *= 1.234;
  expected_exception = 1;
end_try_catch
set_default_error_handlers();
assert(!expected_exception);
assert_LALUnit(u1.norm(), u1);
clear u1;
clear u2;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED LALUnit operations");

## passed all tests!
disp("PASSED all tests");
if exist("swig_exit")
   swig_exit;
endif
