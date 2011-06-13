// SWIG interface to the LALPulsar library

// common SWIG interface code
%include <lal/swiglal-common.i>

// import LAL definitions
%import <lal/swig-lal.i>

// import LALXML definitions
#ifdef SWIGLAL_USE_LALXML
%import <lal/swig-lalxml.i>
#endif

// LALPulsar header files which form the SWIG interface
