## Check SWIG Octave bindings for LALSimulation
## Author: Karl Wette, 2011--2014

crash_dumps_octave_core(0);

## check module load
disp("checking module load ...");
lalsimulation;
assert(exist("lalsimulation", "var"));
lal;
assert(exist("lal", "var"));
disp("PASSED module load");

## check object parent tracking
disp("checking object parent tracking ...");
a = lalsimulation.new_swig_lalsimulation_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalsimulation.swig_lalsimulation_test_parent_map.s;
  lalsimulation.swig_lalsimulation_test_parent_map.s = lal.swig_lal_test_struct_const;
endfor
clear c;
clear b;
clear a;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
