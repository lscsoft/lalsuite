## Check SWIG Octave bindings for LALMetaIO
## Author: Karl Wette, 2011--2014

crash_dumps_octave_core(0);

## check module load
disp("checking module load ...");
lalmetaio;
assert(exist("lalmetaio", "var"));
lal;
assert(exist("lal", "var"));
disp("PASSED module load");

## check object parent tracking
disp("checking object parent tracking ...");
a = lalmetaio.new_swig_lalmetaio_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalmetaio.swig_lalmetaio_test_parent_map.s;
  lalmetaio.swig_lalmetaio_test_parent_map.s = lal.swig_lal_test_struct_const;
endfor
clear c;
clear b;
clear a;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
