## make Octave extension from SWIG wrapping code
## Author: Karl Wette, 2011

## this is a script
1;

## get value of environment variable,
## or raise error if it is not defined.
function str = getdefenv(name)
  str = getenv(name);
  if isempty(str)
    error("Required environment variable '%s' is not defined!", name);
  endif
endfunction

## macro definitions
defflags = strcat("-D", strsplit(getdefenv("swig_defines"), " ", true));

## turn off optimisation unless SWIGLAL_NDEBUG is defined
if isempty(strfind([" ", getdefenv("swig_defines"), " "], " SWIGLAL_NDEBUG "))
  setenv("CFLAGS", [octave_config_info("CFLAGS"), " -O0"]);
  setenv("CXXFLAGS", [octave_config_info("CXXFLAGS"), " -O0"]);
endif

## include directories
inclflags = strcat("-I", strsplit(getdefenv("swig_inclpath"), " ", true));

## compile-time library directories
libflags = strcat("-L", strsplit(getdefenv("swig_libpath"), " ", true));

## run-time library directory
rtlibdir = strcat("-Wl,-rpath=", getdefenv("swig_libdir"));

## libraries to link against
libs = strcat("-l", strsplit(getdefenv("swig_libs"), " ", true));

## source file
srcfile = getdefenv("swig_wrapfile");

## get path of mkoctfile executable
mkoctfilebin = fullfile(
                        octave_config_info("bindir"),
                        sprintf("mkoctfile-%s", OCTAVE_VERSION)
                        );

## arguments to mkoctfile
cmdargs = {"-v", "-Wall", \
           defflags{:}, inclflags{:}, \
           libflags{:}, libs{:}, rtlibdir, \
           srcfile};

## build command line
cmd = strcat("'", mkoctfilebin, "'");
for i = 1:length(cmdargs)
  cmd = strcat(cmd, " '", cmdargs{i}, "'");
endfor

## execute mkoctfile
status = system(cmd);

## return status
exit(status);
