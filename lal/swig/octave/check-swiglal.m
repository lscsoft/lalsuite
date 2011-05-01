## check SWIG Octave module wrapping the LAL library
## Author: Karl Wette, 2011

function msg(str)
  disp([program_name, ": ", str]);
endfunction

## check memory allocation
if !cvar.swiglal_debug
  msg("skipping memory allocation");
else
  LALCheckMemoryLeaks();
  mem1 = new_LALDetector();
  mem2 = new_LALStringVector();
  mem3 = new_COMPLEX8Vector();
  mem4 = XLALCreateREAL8Vector(3);
  msg("*** below should be an error message from LALCheckMemoryLeaks() ***");
  try
    LALCheckMemoryLeaks();
    error("expected exception");
  end_try_catch
  msg("*** above should be an error message from LALCheckMemoryLeaks() ***");
  clear mem1 mem2 mem3;
  XLALDestroyREAL8Vector(mem4);
  LALCheckMemoryLeaks();
  msg("passed memory allocation");
endif

## passed all tests!
msg("================");
msg("PASSED all tests");
msg("================");
