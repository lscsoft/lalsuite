## Check SWIG Octave module wrapping lalinference
## Author: Karl Wette, 2011, 2012

## check module load
lalinference;
assert(exist("lalinference", "var"));
assert(exist("lalinferencecvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
lalcvar.lalDebugLevel = bitor(LALERROR, LALMEMTRACE);
disp("passed module load");

# check object parent tracking
a = lalinference.new_lalinferenceswig_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalinferencecvar.lalinferenceswig_test_parent_map.s;
  lalinferencecvar.lalinferenceswig_test_parent_map.s = lalcvar.lalswig_test_struct_const;
endfor
clear a b c;
CheckMemoryLeaks();
disp("passed object parent tracking");

## passed all tests!
disp("PASSED all tests");
