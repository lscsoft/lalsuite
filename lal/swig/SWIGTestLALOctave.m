## Check SWIG Octave bindings for LAL
## Author: Karl Wette, 2011--2014

crash_dumps_octave_core(0);
expected_exception = 0;

## check module load
disp("checking module load ...");
lal;
assert(exist("lal", "var"));
disp("PASSED module load");

## check memory allocation
disp("checking memory allocation ...");
if !lal.NoDebug
  LALCheckMemoryLeaks();
  mem1 = new_LALDetector();
  mem2 = XLALCreateCOMPLEX8Vector(5);
  mem3 = XLALCreateREAL8Vector(3);
  mem4 = XLALCreateREAL4TimeSeries("test", LIGOTimeGPS(0), 100, 0.1, lal.DimensionlessUnit, 10);
  disp("*** below should be an error message from CheckMemoryLeaks() ***");
  try
    LALCheckMemoryLeaks();
    expected_exception = 1;
  end_try_catch
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
assert(strcmp(typeinfo(b), "swiglal_oct_array_view_double"));
assert(all(b == [1.1; 2.2; 3.3]));
clear a;
clear ans;
assert(all(b == [1.1; 2.2; 3.3]));
ts = XLALCreateREAL8TimeSeries("test", LIGOTimeGPS(0), 0, 0.1, lal.DimensionlessUnit, 10);
ts.data.data = (1:10)';
for i = 1:7
  v = ts.data;
endfor
assert(v.data == (1:10)');
clear ts;
clear ans;
assert(v.data == (1:10)');
clear v;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED object parent tracking");

## check equal return/first argument type handling
disp("checking equal return/first argument type handling");
sv = XLALCreateStringVector("1");
assert(sv.length == 1);
XLALAppendString2Vector(sv, "2");
assert(sv.length == 2);
sv = XLALAppendString2Vector(sv, "3");
assert(sv.length == 3);
sv2 = XLALAppendString2Vector(sv, "4");
assert(sv.length == 4);
assert(sv2.length == 4);
assert(swig_this(sv) == swig_this(sv2));
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
assert(ts.data.length == 40);
assert(ts2.data.length == 40);
assert(swig_this(ts) == swig_this(ts2));
clear ts;
clear ts2;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED equal return/first argument type handling");

## check string conversions
disp("checking string conversions ...");
strs = {"a"; "bc"; "def"};
sv = XLALCreateStringVector(strs{:});
assert(sv.length == 3);
assert(all(strcmp(sv.data, strs)));
strs{1} = "ghijk";
sv.data{1} = strs{1};
strs{end+1} = "lmnopq";
sv = XLALAppendString2Vector(sv, strs{4});
assert(sv.length == 4);
for i = 1:4
  assert(strcmp(sv.data{i}, strs{i}));
endfor
clear sv;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED string conversions");

## check static vector/matrix conversions
disp("checking static vector/matrix conversions ...");
lal.swig_lal_test_struct_vector{1} = lal.swig_lal_test_struct_const;
assert(lal.swig_lal_test_struct_vector{1}.n == lal.swig_lal_test_struct_const.n);
assert(lal.swig_lal_test_struct_vector{1}.i == lal.swig_lal_test_struct_const.i);
assert(lal.swig_lal_test_struct_vector{1}.f == lal.swig_lal_test_struct_const.f);
assert(strcmp(lal.swig_lal_test_struct_vector{1}.str, lal.swig_lal_test_struct_const.str));
assert(all(lal.swig_lal_test_struct_vector{1}.vec == lal.swig_lal_test_struct_const.vec));
lal.swig_lal_test_struct_matrix{1, 1} = lal.swig_lal_test_struct_const;
assert(lal.swig_lal_test_struct_matrix{1, 1}.n == lal.swig_lal_test_struct_const.n);
assert(lal.swig_lal_test_struct_matrix{1, 1}.i == lal.swig_lal_test_struct_const.i);
assert(lal.swig_lal_test_struct_matrix{1, 1}.f == lal.swig_lal_test_struct_const.f);
assert(strcmp(lal.swig_lal_test_struct_matrix{1, 1}.str, lal.swig_lal_test_struct_const.str));
assert(all(lal.swig_lal_test_struct_matrix{1, 1}.vec == lal.swig_lal_test_struct_const.vec));
sts = new_swig_lal_test_struct();
assert(length(sts.vec) == 3);
assert(length(sts.evec) == 3);
assert(all(size(sts.mat) == [2, 3]));
sts.vec = [3; 2; 1];
assert(all(sts.vec == [3; 2; 1]));
sts.mat = [4, 5, 6; 9, 8, 7];
try
  sts.mat = [1.1, 2.3, 4.5; 6.5, 4.3, 2.1];
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
assert(all(all(sts.mat == [4, 5, 6; 9, 8, 7])));
for i = 1:3
  sts.evec(i) = 2*i + 3;
  assert(sts.evec(i) == (2*i + 3));
endfor
clear sts;
clear ans;
assert(!any(lal.swig_lal_test_enum_vector));
assert(!any(lal.swig_lal_test_enum_matrix(:)));
assert(length(lal.swig_lal_test_empty_INT4_vector) == 0);
assert(!any(lal.swig_lal_test_INT4_vector));
assert(!any(lal.swig_lal_test_INT4_matrix(:)));
assert(!any(lal.swig_lal_test_REAL8_vector));
assert(!any(lal.swig_lal_test_REAL8_matrix(:)));
assert(!any(lal.swig_lal_test_COMPLEX8_vector));
assert(!any(lal.swig_lal_test_COMPLEX8_matrix(:)));
lal.swig_lal_test_INT4_vector(1) = 10;
assert(lal.swig_lal_test_INT4_vector(1) == 10);
lal.swig_lal_test_INT4_matrix(1, 1) = 11;
assert(lal.swig_lal_test_INT4_matrix(1, 1) == 11);
lal.swig_lal_test_INT4_vector = lal.swig_lal_test_INT4_const_vector;
assert(all(lal.swig_lal_test_INT4_vector == [1; 2; 4]));
assert(lal.swig_lal_test_INT4_const_vector(3) == 4);
lal.swig_lal_test_INT4_matrix = lal.swig_lal_test_INT4_const_matrix;
assert(all(lal.swig_lal_test_INT4_matrix == [[1, 2, 4]; [2, 4, 8]]));
assert(lal.swig_lal_test_INT4_const_matrix(2, 3) == 8);
try
  lal.swig_lal_test_INT4_const_vector(20);
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
lal.swig_lal_test_REAL8_vector(1) = 3.4;
assert(lal.swig_lal_test_REAL8_vector(1) == 3.4);
lal.swig_lal_test_REAL8_matrix(1, 1) = 5.6;
assert(lal.swig_lal_test_REAL8_matrix(1, 1) == 5.6);
lal.swig_lal_test_COMPLEX8_vector(1) = complex(3.5, 4.75);
assert(lal.swig_lal_test_COMPLEX8_vector(1) == complex(3.5, 4.75));
lal.swig_lal_test_COMPLEX8_matrix(1, 1) = complex(5.5, 6.25);
assert(lal.swig_lal_test_COMPLEX8_matrix(1, 1) == complex(5.5, 6.25));
disp("PASSED static vector/matrix conversions");

## check dynamic vector/matrix conversions
disp("checking dynamic vector/matrix conversions ...");
function check_dynamic_vector_matrix(iv, ivl, rv, rvl, cm, cms1, cms2)
  expected_exception = 0;
  iv.data = zeros(ivl, 1);
  rv.data = zeros(rvl, 1);
  cm.data = zeros(cms1, cms2);
  assert(ivl == 5);
  iv.data = [1; 3; 2; 4; 3];
  assert(all(iv.data == [1; 3; 2; 4; 3]));
  iv.data(4) = 7;
  assert(iv.data(4) == 7);
  assert(rvl == 5);
  rv.data = [1.2; 3.4; 2.6; 4.8; 3.5];
  assert(all(rv.data == [1.2; 3.4; 2.6; 4.8; 3.5]));
  rv.data(rvl) = 7.5;
  assert(rv.data(rvl) == 7.5);
  try
    rv.data(rvl + 1) = 99.9;
    expected_exception = 1;
  end_try_catch
  assert(!expected_exception);
  try
    iv.data = rv.data;
    expected_exception = 1;
  end_try_catch
  assert(!expected_exception);
  rv.data = iv.data;
  assert(all(rv.data == iv.data));
  assert(cms1 == 4);
  assert(cms2 == 6);
  for i = 1:cms1
    for j = 1:cms2
      cm.data(i, j) = complex(i / 4.0, j / 2.0);
    endfor
  endfor
  assert(cm.data(2, 3) == complex(0.5, 1.5));
  assert(cm.data(3, 2) == complex(0.75, 1.0));
  try
    iv.data(0) = cm.data(2, 3);
    expected_exception = 1;
  end_try_catch
  assert(!expected_exception);
  try
    rv.data(0) = cm.data(3, 2);
    expected_exception = 1;
  end_try_catch
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
assert(rv0.length == 0);
assert(length(rv0.data) == 0);
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
assert(all(swig_lal_test_copyin_array1(a1in, 2.5) == a1out));
a2in = int32([3,2; 7,6; 12,10]);
a2out = a2in * 15;
assert(all(swig_lal_test_copyin_array2(a2in, 15) == a2out));
try
  swig_lal_test_viewin_array1([0,0,0,0], 0);
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
try
  swig_lal_test_viewin_array2([1.2,3.4; 0,0; 0,0], 0);
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
disp("PASSED fixed and dynamic arrays typemaps")

## check input views of array structs
disp("checking input views of array structs ...");
r4dat = single([1.2; 2.3; 3.4; 4.5; 5.6]);
r8dat = double([3.4; 4.5; 5.6; 6.7; 7.8; 8.9]);
c8dat = complex(r4dat, 8 + r4dat);
c16dat = complex(r8dat, 16 + r8dat);
r4 = XLALCreateREAL4Vector(length(r4dat));
r4.data = r4dat;
r4out = XLALCreateREAL4Vector(length(r4dat));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewin_REAL4Vector(r4out, r4));
assert(all(r4out.data == r4.data));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewin_REAL4Vector(r4out, r4dat));
assert(all(r4out.data == r4dat));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewinout_REAL4Vector(r4out, r4));
assert(all(2 * r4out.data == r4.data));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewinout_REAL4Vector(r4out, r4dat));
assert(all(2 * r4out.data == r4dat));
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
assert(all(r8out.data == r8.data));
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewin_REAL8Vector(r8out, r8dat));
assert(all(r8out.data == r8dat));
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewinout_REAL8Vector(r8out, r8));
assert(all(2 * r8out.data == r8.data));
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewinout_REAL8Vector(r8out, r8dat));
assert(all(2 * r8out.data == r8dat));
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
assert(all(c8out.data == c8.data));
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewin_COMPLEX8Vector(c8out, c8dat));
assert(all(c8out.data == c8dat));
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewinout_COMPLEX8Vector(c8out, c8));
assert(all(2 * c8out.data == c8.data));
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewinout_COMPLEX8Vector(c8out, c8dat));
assert(all(2 * c8out.data == c8dat));
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
assert(all(c16out.data == c16.data));
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewin_COMPLEX16Vector(c16out, c16dat));
assert(all(c16out.data == c16dat));
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewinout_COMPLEX16Vector(c16out, c16));
assert(all(2 * c16out.data == c16.data));
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewinout_COMPLEX16Vector(c16out, c16dat));
assert(all(2 * c16out.data == c16dat));
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
assert(all(r4out.data == r4.data));
r4out.data = zeros(size(r4dat));
assert(swig_lal_test_viewin_REAL4VectorSequence(r4out, r4dat));
assert(all(r4out.data == r4dat));
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
assert(all(r8out.data == r8.data));
r8out.data = zeros(size(r8dat));
assert(swig_lal_test_viewin_REAL8VectorSequence(r8out, r8dat));
assert(all(r8out.data == r8dat));
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
assert(all(c8out.data == c8.data));
c8out.data = zeros(size(c8dat));
assert(swig_lal_test_viewin_COMPLEX8VectorSequence(c8out, c8dat));
assert(all(c8out.data == c8dat));
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
assert(all(c16out.data == c16.data));
c16out.data = zeros(size(c16dat));
assert(swig_lal_test_viewin_COMPLEX16VectorSequence(c16out, c16dat));
assert(all(c16out.data == c16dat));
clear c16;
clear c16out;
clear c16dat;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED input views of array structs (LAL)");
vfdat = single([1.2; 2.3; 3.4; 4.5; 5.6]);
vddat = double([3.4; 4.5; 5.6; 6.7; 7.8; 8.9]);
vcfdat = complex(vfdat, 8 + vfdat);
vcddat = complex(vddat, 16 + vddat);
vf = new_gsl_vector_float(length(vfdat));
vf.data = vfdat;
vfout = new_gsl_vector_float(length(vfdat));
vfout.data = zeros(size(vfdat));
assert(swig_lal_test_viewin_gsl_vector_float(vfout, vf));
assert(all(vfout.data == vf.data));
vfout.data = zeros(size(vfdat));
assert(swig_lal_test_viewin_gsl_vector_float(vfout, vfdat));
assert(all(vfout.data == vfdat));
vfout.data = zeros(size(vfdat));
assert(swig_lal_test_viewinout_gsl_vector_float(vfout, vf));
assert(all(2 * vfout.data == vf.data));
vfout.data = zeros(size(vfdat));
assert(swig_lal_test_viewinout_gsl_vector_float(vfout, vfdat));
assert(all(2 * vfout.data == vfdat));
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
assert(all(vdout.data == vd.data));
vdout.data = zeros(size(vddat));
assert(swig_lal_test_viewin_gsl_vector(vdout, vddat));
assert(all(vdout.data == vddat));
vdout.data = zeros(size(vddat));
assert(swig_lal_test_viewinout_gsl_vector(vdout, vd));
assert(all(2 * vdout.data == vd.data));
vdout.data = zeros(size(vddat));
assert(swig_lal_test_viewinout_gsl_vector(vdout, vddat));
assert(all(2 * vdout.data == vddat));
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
assert(all(vcfout.data == vcf.data));
vcfout.data = zeros(size(vcfdat));
assert(swig_lal_test_viewin_gsl_vector_complex_float(vcfout, vcfdat));
assert(all(vcfout.data == vcfdat));
vcfout.data = zeros(size(vcfdat));
assert(swig_lal_test_viewinout_gsl_vector_complex_float(vcfout, vcf));
assert(all(2 * vcfout.data == vcf.data));
vcfout.data = zeros(size(vcfdat));
assert(swig_lal_test_viewinout_gsl_vector_complex_float(vcfout, vcfdat));
assert(all(2 * vcfout.data == vcfdat));
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
assert(all(vcdout.data == vcd.data));
vcdout.data = zeros(size(vcddat));
assert(swig_lal_test_viewin_gsl_vector_complex(vcdout, vcddat));
assert(all(vcdout.data == vcddat));
vcdout.data = zeros(size(vcddat));
assert(swig_lal_test_viewinout_gsl_vector_complex(vcdout, vcd));
assert(all(2 * vcdout.data == vcd.data));
vcdout.data = zeros(size(vcddat));
assert(swig_lal_test_viewinout_gsl_vector_complex(vcdout, vcddat));
assert(all(2 * vcdout.data == vcddat));
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
assert(all(mfout.data == mf.data));
mfout.data = zeros(size(mfdat));
assert(swig_lal_test_viewin_gsl_matrix_float(mfout, mfdat));
assert(all(mfout.data == mfdat));
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
assert(all(mdout.data == md.data));
mdout.data = zeros(size(mddat));
assert(swig_lal_test_viewin_gsl_matrix(mdout, mddat));
assert(all(mdout.data == mddat));
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
assert(all(mcfout.data == mcf.data));
mcfout.data = zeros(size(mcfdat));
assert(swig_lal_test_viewin_gsl_matrix_complex_float(mcfout, mcfdat));
assert(all(mcfout.data == mcfdat));
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
assert(all(mcdout.data == mcd.data));
mcdout.data = zeros(size(mcddat));
assert(swig_lal_test_viewin_gsl_matrix_complex(mcdout, mcddat));
assert(all(mcdout.data == mcddat));
clear mcd;
clear mcdout;
clear mcddat;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED input views of array structs (GSL)");
function check_input_view_type_safety(f, a, b, expect_exception)
  expected_exception = 0;
  if expect_exception
    try
      f(a, b);
      expected_exception = 1;
    end_try_catch
    assert(!expected_exception);
    try
      f(b, a);
      expected_exception = 1;
    end_try_catch
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
disp("PASSED input views of array structs (type safety)");

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
assert(all(c8out == c8outv.data));
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
assert(all(c16out == c16outv.data));
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
assert(all(c8out == c8outv.data));
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
assert(all(r4out == r4outv.data));
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
assert(all(c16out == c16outv.data));
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
assert(all(r8out == r8outv.data));
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
assert(all(r4out == r4outv.data));
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
assert(all(r8out == r8outv.data));
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
assert(all(r4out == r4outv.data));
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
assert(all(r8out == r8outv.data));
clear r8inv;
clear r8outv;
clear plan;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED FFT functions with input views ...");

## check dynamic array of pointers access
disp("checking dynamic array of pointers access ...");
ap = swig_lal_test_Create_arrayofptrs(3);
assert(ap.length == 3);
for i = 1:ap.length
  assert(ap.data{i}.length == 6);
  for j = 1:ap.data{i}.length
    assert(ap.data{i}.data(j) == 42*ap.length*(i-1) + (j-1));
  endfor
endfor
clear ap;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED dynamic array of pointers access");

## check 'tm' struct conversions
disp("checking 'tm' struct conversions ...");
gps = 989168284;
utc = [2011, 5, 11, 16, 57, 49, 4, 131, 0];
assert(all(XLALGPSToUTC(gps) == utc));
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
  utcd = XLALGPSToUTC(XLALUTCToGPS(utcd));
  assert(utcd(7) == weekday(datenum(utcd(1:6))));
endfor
LALCheckMemoryLeaks();
disp("PASSED 'tm' struct conversions");

## check LIGOTimeGPS operations
disp("checking LIGOTimeGPS operations ...");
t0 = new_LIGOTimeGPS();
assert(t0 == 0 && strcmp(swig_type(t0), "LIGOTimeGPS"));
t1 = new_LIGOTimeGPS(10.5);
t2 = new_LIGOTimeGPS(10, 500000000);
assert(t1 == t2 && strcmp(swig_type(t1), "LIGOTimeGPS"));
t3 = +t1;
t3 = -t2;
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
assert(21 / t1  == 2 && strcmp(swig_type(21 / t1), "LIGOTimeGPS"));
assert(t1 + t2 == 15.5 && strcmp(swig_type(t1 + t2), "LIGOTimeGPS"));
assert(t1 - t2 == 5.5 && strcmp(swig_type(t1 - t2), "LIGOTimeGPS"));
assert(t1 * t2 == 52.5 && strcmp(swig_type(t1 * t2), "LIGOTimeGPS"));
assert(t2 * t1 == 52.5 && strcmp(swig_type(t2 * t1), "LIGOTimeGPS"));
assert(t1 / t2 == 2.1 && strcmp(swig_type(t1 / t2), "LIGOTimeGPS"));
assert(t1 > t2 && t2 < t1 && t1 >= t2 && t2 <= t1)
assert(LIGOTimeGPS(333333333,333333333) == LIGOTimeGPS(1000000000) / 3)
assert(LIGOTimeGPS(666666666,666666667) == LIGOTimeGPS(2000000000) / 3)
t1 += 812345667.75;
assert(strcmp(t1.__str__(), "812345678.250000000"));
assert(new_LIGOTimeGPS(t1.__repr__()) == t1);
assert(t1.ns() == 812345678250000000);
t4struct = new_swig_lal_test_gps;
t4struct.t = 1234.5;
assert(t4struct.t == 1234.5);
t5 = LIGOTimeGPS("1000");
assert(t5 == 1000);
try
  t5 = LIGOTimeGPS("abc1000");
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
try
  t5 = LIGOTimeGPS("1000abc");
  expected_exception = 1;
end_try_catch
assert(!expected_exception);
assert(swig_lal_test_noptrgps(LIGOTimeGPS(1234.5)) == swig_lal_test_noptrgps(1234.5))
clear t0;
clear t1;
clear t2;
clear t3;
clear t4struct;
clear t5;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED LIGOTimeGPS operations");

## passed all tests!
disp("PASSED all tests");
