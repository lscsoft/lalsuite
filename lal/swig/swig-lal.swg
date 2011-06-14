// SWIG interface to the LAL library

// common SWIG interface code
%include <lal/swiglal-common.i>

// SWIG interface to GSL
%include <lal/swiglal-gsl.i>

// tests of SWIG interface code
#if defined(SWIGLAL_IS_LAL) && !defined(SWIGLAL_NDEBUG)
%include "swiglal-test.i"
#endif

// LAL header files which form the SWIG interface
swiglal_include(lal/LALVCSInfo.h)
