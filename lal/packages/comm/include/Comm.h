/*----------------------------------------------------------------------- 
 * 
 * File Name: Comm.h
 *
 * Author: Allen, B., Brown D. A. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _COMM_H
#define _COMM_H

#include <mpi.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

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
  MPI_Comm mpiComm;
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
    INT4        myId,
    MPI_Comm    mpiComm
    );

void
LALMPIDebug (
    LALStatus         *status,
    MPIDebugParams *params
    );

void
LALMPIKillScript (
    LALStatus  *status,
    MPIId   *id,
    MPI_Comm    mpiComm
    );


void
LALMPISendMsg (
    LALStatus     *status,
    MPIMessage *msg,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvMsg (
    LALStatus     *status,
    MPIMessage *msg,
    MPI_Comm    mpiComm
    );

void
LALMPISendCHARVector (
    LALStatus     *status,
    CHARVector *vector,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvCHARVector (
    LALStatus     *status,
    CHARVector *vector,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT2Vector (
    LALStatus     *status,
    INT2Vector *vector,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT2Vector (
    LALStatus     *status,
    INT2Vector *vector,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT4Vector (
    LALStatus     *status,
    INT4Vector *vector,
    INT4        dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT4Vector (
    LALStatus     *status,
    INT4Vector *vector,
    INT4        source,
    MPI_Comm    mpiComm
    );

void
LALMPISendREAL4Vector (
    LALStatus      *status,
    REAL4Vector *vector,
    INT4         dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvREAL4Vector (
    LALStatus      *status,
    REAL4Vector *vector,
    INT4         source,
    MPI_Comm    mpiComm
    );

void
LALMPISendCOMPLEX8Vector (
    LALStatus         *status,
    COMPLEX8Vector *vector,
    INT4            dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvCOMPLEX8Vector (
    LALStatus         *status,
    COMPLEX8Vector *vector,
    INT4            source,
    MPI_Comm    mpiComm
    );

void
LALMPISendINT2TimeSeries (
    LALStatus         *status,
    INT2TimeSeries *series,
    INT4            dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvINT2TimeSeries (
    LALStatus         *status,
    INT2TimeSeries *series,
    INT4            source,
    MPI_Comm    mpiComm
    );

void
LALMPISendREAL4TimeSeries (
    LALStatus          *status,
    REAL4TimeSeries *series,
    INT4             dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvREAL4TimeSeries (
    LALStatus          *status,
    REAL4TimeSeries *series,
    INT4             source,
    MPI_Comm    mpiComm
    );

void
LALMPISendCOMPLEX8TimeSeries (
    LALStatus             *status,
    COMPLEX8TimeSeries *series,
    INT4                dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvCOMPLEX8TimeSeries (
    LALStatus             *status,
    COMPLEX8TimeSeries *series,
    INT4                source,
    MPI_Comm    mpiComm
    );

void
LALMPISendREAL4FrequencySeries (
    LALStatus               *status,
    REAL4FrequencySeries *series,
    INT4                  dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvREAL4FrequencySeries (
    LALStatus               *status,
    REAL4FrequencySeries *series,
    INT4                  source,
    MPI_Comm    mpiComm
    );

void
LALMPISendCOMPLEX8FrequencySeries (
    LALStatus                  *status,
    COMPLEX8FrequencySeries *series,
    INT4                     dest,
    MPI_Comm    mpiComm
    );

void
LALMPIRecvCOMPLEX8FrequencySeries (
    LALStatus                  *status,
    COMPLEX8FrequencySeries *series,
    INT4                     source,
    MPI_Comm    mpiComm
    );


#ifdef  __cplusplus
}
#endif

#endif /* _COMM_H */
