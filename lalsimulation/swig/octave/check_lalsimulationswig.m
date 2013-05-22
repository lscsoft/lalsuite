## Check SWIG Octave module wrapping lalsimulation
## Author: Karl Wette, 2011, 2012

## check module load
lalsimulation;
assert(exist("lalsimulation", "var"));
assert(exist("lalsimulationcvar", "var"));
lal;
assert(exist("lal", "var"));
assert(exist("lalcvar", "var"));
disp("passed module load");

# check object parent tracking
a = lalsimulation.new_lalsimulationswig_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalsimulationcvar.lalsimulationswig_test_parent_map.s;
  lalsimulationcvar.lalsimulationswig_test_parent_map.s = lalcvar.lalswig_test_struct_const;
endfor
clear a b c;
CheckMemoryLeaks();
disp("passed object parent tracking");

## passed all tests!
disp("PASSED all tests");
