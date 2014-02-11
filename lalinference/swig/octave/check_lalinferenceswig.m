## Check SWIG Octave module wrapping lalinference
## Author: Karl Wette, 2011, 2012

## check module load
disp("checking module load ...");
lalinference;
assert(exist("lalinference", "var"));
assert(exist("lalinferencecvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
disp("PASSED module load");

## check object parent tracking
disp("checking object parent tracking ...");
a = lalinference.new_swig_lalinference_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalinferencecvar.swig_lalinference_test_parent_map.s;
  lalinferencecvar.swig_lalinference_test_parent_map.s = lalcvar.swig_lal_test_struct_const;
endfor
clear ans a b c;
CheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
