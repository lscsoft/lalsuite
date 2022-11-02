## Check SWIG Octave bindings for LALPulsar
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

## check module load
disp("checking module load ...");
lalpulsar;
assert(exist("lalpulsar"));
lal;
assert(exist("lal"));
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

## check object parent tracking
disp("checking object parent tracking ...");
a = lalpulsar.new_swig_lalpulsar_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalpulsar.swig_lalpulsar_test_parent_map.s;
  lalpulsar.swig_lalpulsar_test_parent_map.s = lal.swig_lal_test_struct_const;
endfor
clear c;
clear b;
clear a;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED object parent tracking");

# check array element assignment
disp("checking array element assignment ...");
mts = XLALCreateMultiLIGOTimeGPSVector(2);
ts0 = XLALCreateTimestampVector(3);
for i = 1:ts0.length
  ts0.data(i) = LIGOTimeGPS(900000000 + i - 1);
endfor
mts.data(1) = ts0;
lal.swig_set_nasty_error_handlers();
clear mts;
clear ts0;
set_default_error_handlers();
LALCheckMemoryLeaks();
disp("PASSED array element assignment");

# check array element parent tracking
disp("checking array element parent tracking");
assert_LIGOTimeGPS(XLALMakeTimestamps(0, 100, 1, 0).data{end}, LIGOTimeGPS(99));
assert_LIGOTimeGPS(XLALMakeMultiTimestamps(0, 100, 1, 0, 3).data{end}.data{end}, LIGOTimeGPS(99));
LALCheckMemoryLeaks();
disp("PASSED array element parent tracking");

# check CWMFDataParams usage
disp("checking CWMFDataParams usage");
function cwmf = make_CWMFDataParams()
  lalpulsar;
  cwmf = CWMFDataParams();
  cwmf.multiTimestamps = XLALMakeMultiTimestamps(0, 10, 1, 0, 3);
endfunction
assert_LIGOTimeGPS(make_CWMFDataParams().multiTimestamps.data{end}.data{end}, LIGOTimeGPS(9));
LALCheckMemoryLeaks();
disp("PASSED CWMFDataParams usage");

## passed all tests!
disp("PASSED all tests");
if exist("swig_exit")
   swig_exit;
endif
