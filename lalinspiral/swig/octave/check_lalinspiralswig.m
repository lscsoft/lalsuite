## Check SWIG Octave module wrapping lalinspiral
## Author: Karl Wette, 2011, 2012

## check module load
disp("checking module load ...");
lalinspiral;
assert(exist("lalinspiral", "var"));
assert(exist("lalinspiralcvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
disp("PASSED module load");

## check object parent tracking
disp("checking object parent tracking ...");
a = lalinspiral.new_lalinspiralswig_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalinspiralcvar.lalinspiralswig_test_parent_map.s;
  lalinspiralcvar.lalinspiralswig_test_parent_map.s = lalcvar.lalswig_test_struct_const;
endfor
clear ans a b c;
CheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
