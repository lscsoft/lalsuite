## Check SWIG Octave module wrapping laldetchar
## Author: Karl Wette, 2011, 2012

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
a = laldetchar.new_laldetcharswig_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = laldetcharcvar.laldetcharswig_test_parent_map.s;
  laldetcharcvar.laldetcharswig_test_parent_map.s = lalcvar.lalswig_test_struct_const;
endfor
clear ans a b c;
CheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
