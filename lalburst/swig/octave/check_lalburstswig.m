## Check SWIG Octave module wrapping lalburst
## Author: Karl Wette, 2011, 2012

## check module load
lalburst;
assert(exist("lalburst", "var"));
assert(exist("lalburstcvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
lalcvar.lalDebugLevel = 1;
disp("passed module load");

# check object parent tracking
a = lalburst.new_lalburstswig_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalburstcvar.lalburstswig_test_parent_map.s;
  lalburstcvar.lalburstswig_test_parent_map.s = lalcvar.lalswig_test_struct_const;
endfor
clear a b c;
CheckMemoryLeaks();
disp("passed object parent tracking");

## passed all tests!
disp("PASSED all tests");
