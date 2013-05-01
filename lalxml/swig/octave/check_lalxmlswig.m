## Check SWIG Octave module wrapping lalxml
## Author: Karl Wette, 2011, 2012

## check module load
lalxml;
assert(exist("lalxml", "var"));
assert(exist("lalxmlcvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
lalcvar.lalDebugLevel = bitor(LALERROR, LALMEMTRACE);
disp("passed module load");

# check object parent tracking
a = lalxml.new_lalxmlswig_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalxmlcvar.lalxmlswig_test_parent_map.s;
  lalxmlcvar.lalxmlswig_test_parent_map.s = lalcvar.lalswig_test_struct_const;
endfor
clear a b c;
CheckMemoryLeaks();
disp("passed object parent tracking");

## passed all tests!
disp("PASSED all tests");
