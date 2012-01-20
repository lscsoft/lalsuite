## check SWIG Octave module wrapping the LALSimulation library
## Author: Karl Wette, 2011

function msg(str)
  disp([program_name, ": ", str]);
endfunction

## check module load
swiglalsimulation;
msg("passed module load");

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
