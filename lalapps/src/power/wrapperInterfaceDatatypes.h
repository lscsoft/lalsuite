/* $Id$ */

#ifndef WRAPPER_INTERFACE_DATATYPES_H
#define WRAPPER_INTERFACE_DATATYPES_H

/* MPI Header File */
#include <mpi.h>   

/* Local Header Files */
#include <lal/LALAtomicDatatypes.h>
   
   
#ifdef __cplusplus
extern "C"{
#endif   

   
/* Domain type */  
typedef enum
{ 
   timeD,
   freqD,
   bothD,
   databaseD
} domain;

   
/* Datatype */
typedef enum
{
   boolean_lu,
   char_s,
   char_u,
   int_2s,
   int_2u,
   int_4s,
   int_4u,
   int_8s,
   int_8u,
   real_4,
   real_8,
   complex_8,
   complex_16,
   char_s_ptr,
   char_u_ptr
}datatype;
   

/* Pointer type checking */   
typedef union
{
   BOOLEAN*   boolean;
   CHAR*      chars;
   UCHAR*     charu;
   INT2*      int2s;
   UINT2*     int2u;
   INT4*      int4s;
   UINT4*     int4u;
   INT8*      int8s;
   UINT8*     int8u;   
   REAL4*     real4;
   REAL8*     real8;
   COMPLEX8*  complex8;
   COMPLEX16* complex16;   
}dataPointer;


#define maxDetectorName 256

/* Detector geometry */
typedef struct
{
  CHAR  name[maxDetectorName];
  REAL8 longitude;    /* of vertex, decimal degrees */
  REAL8 latitude;     /* of vertex, decimal degrees */
  REAL4 elevation;    /* of vertex, meters above WGS-84 ellipsoid */
  REAL4 armXazimuth;  /* Radians North of East */
  REAL4 armYazimuth;  /* Radians North of East */
  REAL4 armXaltitude; /* Radians above tangent to ellipsoid */
  REAL4 armYaltitude; /* Radians above tangent to ellipsoid */
  INT4  localTime;    /* Offset from GMT in seconds */
  UINT4 dataQuality;
} detGeom;
   

/* Time domain interval */   
typedef struct 
{
   UINT8 numberSamples;
   UINT4 startSec;
   UINT4 startNan;
   UINT4 stopSec;
   UINT4 stopNan;   
   REAL8 timeStepSize;
   REAL8 baseFreq;
   REAL8 phase;
   UINT4 numberChannels;
   UINT4 numberDetectors;   
   CHAR** channelName;
   detGeom* geometry;
}gpsTimeInterval;   
   

/* Frequency domain interval */   
typedef struct
{
   UINT8 numberSamples;
   UINT4 gpsStartTimeSec;
   UINT4 gpsStartTimeNan;   
   UINT4 gpsStopTimeSec;
   UINT4 gpsStopTimeNan;      
   REAL8 startFreq;   
   REAL8 stopFreq;      
   REAL8 freqStepSize;
   REAL8 baseFreq;
   UINT4 numberChannels;
   UINT4 numberDetectors;
   CHAR** channelName;   
   detGeom* geometry;   
}frequencyInterval;
   

/* Time & Frequency domain interval */   
typedef struct
{
   UINT8 numberSamples;
   UINT4 gpsStartTimeSec;
   UINT4 gpsStartTimeNan;    
   UINT4 gpsStopTimeSec;
   UINT4 gpsStopTimeNan;       
   REAL8 startFreq;   
   REAL8 stopFreq; 
   REAL8 timeStepSize;
   REAL8 freqStepSize;
   UINT4 numberChannels;
   UINT4 numberDetectors;
   CHAR** channelName;   
   detGeom* geometry;   
}timeFreqInterval;

   
/* Maximum length of SQL query supported by db2: 4K+16   */
#define maxSQL 4112      

   
/* Database domain interval */   
typedef struct
{
   UINT8 numberRows;
   CHAR sqlQuery[ maxSQL ];
}databaseInterval;
   
  
/* Interval domain */   
typedef union
{
   gpsTimeInterval   dTime;
   frequencyInterval dFreq;
   timeFreqInterval  dBoth;
   databaseInterval  dDatabase;
}interval;


#define maxHistoryName  64   
#define maxHistoryUnits 64   

   
/* History linked list */   
typedef struct dcHistoryTag
{
   struct dcHistoryTag* previous;
   CHAR name[maxHistoryName];
   CHAR units[maxHistoryUnits];
   datatype type;
   UINT4 numberValues;
   dataPointer value;
   struct dcHistoryTag* next;
}dcHistory;

#define maxMultiDimName    256
#define maxMultiDimUnits   256   
#define maxMultiDimComment 256   
   
/* Multidimensional data */   
typedef struct 
{
   CHAR name[maxMultiDimName];
   CHAR units[maxMultiDimUnits];
   CHAR comment[maxMultiDimComment];   
   domain space;
   datatype type;
   interval range;
   UINT4 numberDimensions;
   UINT4* dimensions;
   dcHistory* history;
   dataPointer data;
}multiDimData;

   
#define maxStateName 64
   
   
/* State vector */   
typedef struct stateVectorTag
{
   struct stateVectorTag* previous;
   CHAR stateName[maxStateName];
   multiDimData* store;
   struct stateVectorTag* next;
}stateVector;
   

/* Input data structure */   
typedef struct
{
   UINT4 numberSequences;
   stateVector* states;
   multiDimData* sequences;
}inPut;
   

/* Astrophysical/instrumental search catagories */   
typedef enum
{
   binaryInspiral,
   ringDown,
   periodic,
   burst,
   stocastic,
   timeFreq,
   instrumental,
   protoType,
   experimental
}catagory;

   
#define dbNameLimit 19
   

/* Database structure */   
typedef struct dataBaseTag
{
   struct dataBaseTag* previous;
   CHAR tableName[dbNameLimit];
   CHAR columnName[dbNameLimit];
   datatype type;
   UINT4 numberRows;
   dataPointer rows;
   UINT4* rowDimensions;   /* Non zero only for char_s_ptr or char_u_ptr */
   struct dataBaseTag* next;
}dataBase;

   
/* Output data structure */   
typedef struct 
{
   INT8 indexNumber;
   catagory search;
   BOOLEAN significant;
   stateVector* states;
   dataBase* results;
   multiDimData* optional;
}outPut;
   

/* Structures specific to the search functions */

typedef struct
{
   INT4    argc;
   CHAR**  argv;
   INT8    startTime;     /* seconds since January 1, 1970 */
   INT8    dataDuration;  /* length of data rounded up to nearest second */
   REAL4   realtimeRatio; /* from wrapperAPI command line */
   INT4    nodeClass;     /* 0 if searchSlave, 1 if searchMaster */
}InitParams;
   
   
typedef struct 
{
   INT4    add;       /* number of nodes added or subtracted */
   BOOLEAN mpiAPIio;  /* command from mpiAPI: false - exit, true - continue */
}MPIapiAction;

   
typedef struct
{
   MPI_Comm*     comm;   /* wrapper slave COMM_WORLD */
   MPIapiAction* action; /* instruction from mpiAPI to search code */
}SearchParams;
   

typedef struct 
{
   INT4    numOutput;
   outPut* result;
   REAL4   fracRemaining; /* fraction of search remaining */
   BOOLEAN notFinished;   /* false indicates that applySearch is finished */
}SearchOutput;
   
   
#ifdef __cplusplus
}
#endif   
   
   
#endif   
