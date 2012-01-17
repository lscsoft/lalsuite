// SWIG interface to the LALFrame library

// common SWIG interface code
%include <lal/swiglal-common.i>

// import LAL definitions
%import <lal/swig-lal.i>

// allocate these structs with LAL memory
SWIGLAL_STRUCT_LALMEM(FrStat);
SWIGLAL_STRUCT_LALMEM(FrCache);
SWIGLAL_STRUCT_LALMEM(FrCacheSieve);
SWIGLAL_STRUCT_LALMEM(FrameData);
SWIGLAL_STRUCT_LALMEM(FrFileInfo);
SWIGLAL_STRUCT_LALMEM(FrStream);
SWIGLAL_STRUCT_LALMEM(FrPos);
SWIGLAL_STRUCT_LALMEM(FrChanIn);
SWIGLAL_STRUCT_LALMEM(FrOutPar);

// LALFrame header files which form the SWIG interface
