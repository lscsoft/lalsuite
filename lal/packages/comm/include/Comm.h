/*----------------------------------------------------------------------- 
 * 
 * File Name: Comm.h
 *
 * Author: Allen, B. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _COMM_H
#define _COMM_H

#include "LALDatatypes.h"
#include "AVFactories.h"
#include "mpi.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (COMMH, "$Id$");

#define COMM_ENULL 1
#define COMM_ESIZE 2
#define COMM_ESTRL 4
#define COMM_EMPIE 8
#define COMM_ESENV 16
#define COMM_ESYSC 32

#define COMM_MSGENULL "Null pointer"
#define COMM_MSGESIZE "Invalid size"
#define COMM_MSGESTRL "String too long"
#define COMM_MSGEMPIE "MPI error"
#define COMM_MSGESENV "Couldn't set environment variable"
#define COMM_MSGESYSC "Error executing system command"

/* Structure for identifying processors */
typedef struct
tagMPIId
{
  INT4 numProcs;
  INT4 myId;
  INT4 nameLen;
  CHAR procName[MPI_MAX_PROCESSOR_NAME];
}
MPIId;

typedef struct
tagMPIDebugParams
{
  CHAR *debugger;
  CHAR *progName;
  INT4  delay;
  INT4  myId;
}
MPIDebugParams;

typedef enum
{
  MPIDone,
  MPIMisc,
  MPIErr,
  MPIMsg,
  MPICHARVector,
  MPICHARVectorData,
  MPII2Vector,
  MPII2VectorData,
  MPII4Vector,
  MPII4VectorData,
  MPII8Vector,
  MPII8VectorData,
  MPISVector,
  MPISVectorData,
  MPIDVector,
  MPIDVectorData,
  MPICVector,
  MPICVectorData,
  MPIZVector,
  MPIZVectorData,
  MPII2TimeSeries,
  MPII4TimeSeries,
  MPII8TimeSeries,
  MPISTimeSeries,
  MPIDTimeSeries,
  MPICTimeSeries,
  MPIZTimeSeries,
  MPII2FrequencySeries,
  MPII4FrequencySeries,
  MPII8FrequencySeries,
  MPISFrequencySeries,
  MPIDFrequencySeries,
  MPICFrequencySeries,
  MPIZFrequencySeries
}
MPIMsgCode;

typedef struct
tagMPIMessage
{
  INT4 msg;
  INT4 send;
  INT4 source;
}
MPIMessage;


void
LALMPIExportEnvironment (
    LALStatus     *status,
    const CHAR *env,
    INT4        myId
    );

void
LALMPIDebug (
    LALStatus         *status,
    MPIDebugParams *params
    );

void
LALMPIKillScript (
    LALStatus  *status,
    MPIId   *id
    );


void
LALMPISendMsg (
    LALStatus     *status,
    MPIMessage *msg,
    INT4        dest
    );

void
LALMPIRecvMsg (
    LALStatus     *status,
    MPIMessage *msg
    );

void
LALMPISendCHARVector (
    LALStatus     *status,
    CHARVector *vector,
    INT4        dest
    );

void
LALMPIRecvCHARVector (
    LALStatus     *status,
    CHARVector *vector,
    INT4        source
    );

void
LALMPISendINT2Vector (
    LALStatus     *status,
    INT2Vector *vector,
    INT4        dest
    );

void
LALMPIRecvINT2Vector (
    LALStatus     *status,
    INT2Vector *vector,
    INT4        source
    );

void
LALMPISendREAL4Vector (
    LALStatus      *status,
    REAL4Vector *vector,
    INT4         dest
    );

void
LALMPIRecvREAL4Vector (
    LALStatus      *status,
    REAL4Vector *vector,
    INT4         source
    );

void
LALMPISendCOMPLEX8Vector (
    LALStatus         *status,
    COMPLEX8Vector *vector,
    INT4            dest
    );

void
LALMPIRecvCOMPLEX8Vector (
    LALStatus         *status,
    COMPLEX8Vector *vector,
    INT4            source
    );

void
LALMPISendINT2TimeSeries (
    LALStatus         *status,
    INT2TimeSeries *series,
    INT4            dest
    );

void
LALMPIRecvINT2TimeSeries (
    LALStatus         *status,
    INT2TimeSeries *series,
    INT4            source
    );

void
LALMPISendREAL4TimeSeries (
    LALStatus          *status,
    REAL4TimeSeries *series,
    INT4             dest
    );

void
LALMPIRecvREAL4TimeSeries (
    LALStatus          *status,
    REAL4TimeSeries *series,
    INT4             source
    );

void
LALMPISendCOMPLEX8TimeSeries (
    LALStatus             *status,
    COMPLEX8TimeSeries *series,
    INT4                dest
    );

void
LALMPIRecvCOMPLEX8TimeSeries (
    LALStatus             *status,
    COMPLEX8TimeSeries *series,
    INT4                source
    );

void
LALMPISendREAL4FrequencySeries (
    LALStatus               *status,
    REAL4FrequencySeries *series,
    INT4                  dest
    );

void
LALMPIRecvREAL4FrequencySeries (
    LALStatus               *status,
    REAL4FrequencySeries *series,
    INT4                  source
    );

void
LALMPISendCOMPLEX8FrequencySeries (
    LALStatus                  *status,
    COMPLEX8FrequencySeries *series,
    INT4                     dest
    );

void
LALMPIRecvCOMPLEX8FrequencySeries (
    LALStatus                  *status,
    COMPLEX8FrequencySeries *series,
    INT4                     source
    );


#ifdef  __cplusplus
}
#endif

#endif /* _COMM_H */
