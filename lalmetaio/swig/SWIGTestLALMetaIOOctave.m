## Check SWIG Octave bindings for LALMetaIO
## Author: Karl Wette, 2011--2014

page_screen_output(0);
crash_dumps_octave_core(0);

## check module load
disp("checking module load ...");
lalmetaio;
assert(exist("lalmetaio"));
lal;
assert(exist("lal"));
disp("PASSED module load");

# set error handlers
function set_nice_error_handlers()
  swig_set_nice_error_handlers();
endfunction
function set_default_error_handlers()
  lal;
  if length(getenv("NASTY_ERROR_HANDLERS")) > 0
    swig_set_nasty_error_handlers();
  else
    swig_set_nice_error_handlers();
  endif
endfunction
set_default_error_handlers();

## check object parent tracking
disp("checking object parent tracking ...");
a = lalmetaio.new_swig_lalmetaio_test_parent_map_struct();
for i = 1:7
  b = a.s;
  c = lalmetaio.swig_lalmetaio_test_parent_map.s;
  lalmetaio.swig_lalmetaio_test_parent_map.s = lal.swig_lal_test_struct_const;
endfor
clear c;
clear b;
clear a;
clear ans;
LALCheckMemoryLeaks();
disp("PASSED object parent tracking");

## passed all tests!
disp("PASSED all tests");
if exist("swig_exit")
   swig_exit;
endif
