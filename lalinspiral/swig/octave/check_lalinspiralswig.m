## Check SWIG Octave module wrapping lalinspiral
## Author: Karl Wette, 2011, 2012

## check module load
lalinspiral;
assert(exist("lalinspiral", "var"));
assert(exist("lalinspiralcvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
lalcvar.lalDebugLevel = bitor(LALERROR, LALMEMDBG);
disp("passed module load");

# check object parent tracking
a = lalinspiral.new_lalinspiralswig_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalinspiralcvar.lalinspiralswig_test_parent_map.s;
  lalinspiralcvar.lalinspiralswig_test_parent_map.s = lalcvar.lalswig_test_struct_const;
endfor
clear a b c;
CheckMemoryLeaks();
disp("passed object parent tracking");

## passed all tests!
disp("PASSED all tests");
