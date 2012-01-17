// SWIG interface to the LALXML library

// common SWIG interface code
%include <lal/swiglal-common.i>

// import LAL definitions
%import <lal/swig-lal.i>

// allocate these structs with LAL memory
SWIGLAL_STRUCT_LALMEM(XML_NAMESPACE);
SWIGLAL_STRUCT_LALMEM(XML_NAMESPACE_VECTOR);
SWIGLAL_STRUCT_LALMEM(VOTField);
SWIGLAL_STRUCT_LALMEM(VOTFieldVector);

// LALXML header files which form the SWIG interface
