## Check SWIG Octave bindings for LALPulsar
## Author: Karl Wette, 2011--2014

crash_dumps_octave_core(0);

## check module load
disp("checking module load ...");
lalpulsar;
assert(exist("lalpulsar", "var"));
assert(exist("lalpulsarcvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
disp("PASSED module load");

## check object parent tracking
disp("checking object parent tracking ...");
a = lalpulsar.new_swig_lalpulsar_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalpulsarcvar.swig_lalpulsar_test_parent_map.s;
  lalpulsarcvar.swig_lalpulsar_test_parent_map.s = lalcvar.swig_lal_test_struct_const;
endfor
clear c;
clear b;
clear a;
clear ans;
CheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
