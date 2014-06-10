## Check SWIG Octave bindings for LALDetChar
## Author: Karl Wette, 2011--2014

crash_dumps_octave_core(0);

## check module load
disp("checking module load ...");
laldetchar;
assert(exist("laldetchar", "var"));
assert(exist("laldetcharcvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
disp("PASSED module load");

## check object parent tracking
disp("checking object parent tracking ...");
a = laldetchar.new_swig_laldetchar_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = laldetcharcvar.swig_laldetchar_test_parent_map.s;
  laldetcharcvar.swig_laldetchar_test_parent_map.s = lalcvar.swig_lal_test_struct_const;
endfor
clear c;
clear b;
clear a;
clear ans;
CheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
