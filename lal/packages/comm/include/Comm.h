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
MPIExportEnvironment (
    Status     *status,
    const CHAR *env,
    INT4        myId
    );

void
MPIDebug (
    Status         *status,
    MPIDebugParams *params
    );

void
MPIKillScript (
    Status  *status,
    MPIId   *id
    );


void
MPISendMsg (
    Status     *status,
    MPIMessage *msg,
    INT4        dest
    );

void
MPIRecvMsg (
    Status     *status,
    MPIMessage *msg
    );

void
MPISendCHARVector (
    Status     *status,
    CHARVector *vector,
    INT4        dest
    );

void
MPIRecvCHARVector (
    Status     *status,
    CHARVector *vector,
    INT4        source
    );

void
MPISendINT2Vector (
    Status     *status,
    INT2Vector *vector,
    INT4        dest
    );

void
MPIRecvINT2Vector (
    Status     *status,
    INT2Vector *vector,
    INT4        source
    );

void
MPISendREAL4Vector (
    Status      *status,
    REAL4Vector *vector,
    INT4         dest
    );

void
MPIRecvREAL4Vector (
    Status      *status,
    REAL4Vector *vector,
    INT4         source
    );

void
MPISendCOMPLEX8Vector (
    Status         *status,
    COMPLEX8Vector *vector,
    INT4            dest
    );

void
MPIRecvCOMPLEX8Vector (
    Status         *status,
    COMPLEX8Vector *vector,
    INT4            source
    );

void
MPISendINT2TimeSeries (
    Status         *status,
    INT2TimeSeries *series,
    INT4            dest
    );

void
MPIRecvINT2TimeSeries (
    Status         *status,
    INT2TimeSeries *series,
    INT4            source
    );

void
MPISendREAL4TimeSeries (
    Status          *status,
    REAL4TimeSeries *series,
    INT4             dest
    );

void
MPIRecvREAL4TimeSeries (
    Status          *status,
    REAL4TimeSeries *series,
    INT4             source
    );

void
MPISendCOMPLEX8TimeSeries (
    Status             *status,
    COMPLEX8TimeSeries *series,
    INT4                dest
    );

void
MPIRecvCOMPLEX8TimeSeries (
    Status             *status,
    COMPLEX8TimeSeries *series,
    INT4                source
    );

void
MPISendREAL4FrequencySeries (
    Status               *status,
    REAL4FrequencySeries *series,
    INT4                  dest
    );

void
MPIRecvREAL4FrequencySeries (
    Status               *status,
    REAL4FrequencySeries *series,
    INT4                  source
    );

void
MPISendCOMPLEX8FrequencySeries (
    Status                  *status,
    COMPLEX8FrequencySeries *series,
    INT4                     dest
    );

void
MPIRecvCOMPLEX8FrequencySeries (
    Status                  *status,
    COMPLEX8FrequencySeries *series,
    INT4                     source
    );


#ifdef  __cplusplus
}
#endif

#endif /* _COMM_H */
