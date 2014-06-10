## Check SWIG Octave bindings for LALXML
## Author: Karl Wette, 2011--2014

crash_dumps_octave_core(0);

## check module load
disp("checking module load ...");
lalxml;
assert(exist("lalxml", "var"));
assert(exist("lalxmlcvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
disp("PASSED module load");

## check object parent tracking
disp("checking object parent tracking ...");
a = lalxml.new_swig_lalxml_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalxmlcvar.swig_lalxml_test_parent_map.s;
  lalxmlcvar.swig_lalxml_test_parent_map.s = lalcvar.swig_lal_test_struct_const;
endfor
clear c;
clear b;
clear a;
clear ans;
CheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
