## Check SWIG Octave bindings for LALPulsar
## Author: Karl Wette, 2011--2014

page_screen_output(0);
crash_dumps_octave_core(0);

## check module load
disp("checking module load ...");
lalpulsar;
assert(exist("lalpulsar"));
lal;
assert(exist("lal"));
disp("PASSED module load");

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

# check multi-vector element assignment
disp("checking multi-vector element assignment ...");
mts = XLALCreateMultiLIGOTimeGPSVector(2);
ts0 = XLALCreateTimestampVector(3);
mts.data(1) = ts0;
lal.swig_set_nasty_error_handlers();
clear mts;
clear ts0;
lal.swig_set_nice_error_handlers();
disp("PASSED multi-vector element assignment");

## passed all tests!
disp("PASSED all tests");
if exist("swig_exit")
   swig_exit;
endif
