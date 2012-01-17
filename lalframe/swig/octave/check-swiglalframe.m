## check SWIG Octave module wrapping the LALFrame library
## Author: Karl Wette, 2011

function msg(str)
  disp([program_name, ": ", str]);
endfunction

## check module load
swiglalframe;
msg("passed module load");

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
