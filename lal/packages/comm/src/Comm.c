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

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "Comm.h"
#include "LALStdlib.h"

NRCSID (COMMC, "$Id$");

void
LALMPIExportEnvironment (
    LALStatus     *status,
    const CHAR *env,
    INT4        myId
    )
{
  CHAR command[256];
  INT4 code;

  INITSTATUS (status, "LALMPIExportEnvironment", COMMC);

  /* if master, get environment variable */
  if (myId == 0)
  {
    CHAR *var = NULL;
    INT4  len;

    var = getenv (env);
    ASSERT (var, status, COMM_ENULL, COMM_MSGENULL);

    /* calculate length to pass */
    len = strlen(env) + strlen(var) + 9;
    ASSERT (len < sizeof(command), status, COMM_ESTRL, COMM_MSGESTRL);

    /* copy variable into string */
    sprintf (command, "export %s=%s", env, var);

    /* broadcast it */
    code = MPI_Bcast (command, sizeof(command), MPI_CHAR, 0, MPI_COMM_WORLD);
    ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  }
  else
  {
    /* get environment variable string */
    code = MPI_Bcast (command, sizeof(command), MPI_CHAR, 0, MPI_COMM_WORLD);
    ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

    /* set environment variable */
    code = system (command);
    ASSERT (code == 0, status, COMM_ESENV, COMM_MSGESENV);
  }

  RETURN (status);
}


/* be sure to have exported DISPLAY first! */
void
LALMPIDebug (
    LALStatus         *status,
    MPIDebugParams *params
    )
{
  INT4  code;

  INITSTATUS (status, "LALMPIDebug", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (params, status, COMM_ENULL, COMM_MSGENULL);

  /* startup debugger */
  if (params->debugger)
  {
    CHAR cmdstr[256];
    INT4 procId;
    INT4 len;
    
    ASSERT (params->progName, status, COMM_ENULL, COMM_MSGENULL);

    len = strlen(params->debugger) + strlen(params->progName) + 16;
    ASSERT (len < sizeof(cmdstr), status, COMM_ESTRL, COMM_MSGESTRL);

    procId = (INT4) getpid();

    sprintf (cmdstr, "%s %s %d &", params->debugger, params->progName, procId);

    code = system (cmdstr);
    ASSERT (code == 0, status, COMM_ESYSC, COMM_MSGESYSC);
  }

  /* time to delay startup on each node */
  if (params->delay > 0)
  {
    sleep (params->delay);
  }

  /* synchronize */
  code = MPI_Barrier (MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALMPIKillScript (
    LALStatus  *status,
    MPIId   *id
    )
{
  MPI_Status mpiStatus;
  INT4       procId;
  INT4       code;

  INITSTATUS (status, "LALMPIKillScript", COMMC);

  procId = (INT4) getpid();

  if (id->myId == 0) /* I am the master */
  {
    INT4  slave;
    FILE *fp;

    fp = fopen ("killscript", "w");

    fprintf (fp, "#!/bin/sh\n");
    fprintf (fp, "rsh %s kill -9 %d\n", id->procName, procId);
    fflush (fp);

    for (slave = 1; slave < id->numProcs; ++slave)
    {
      CHAR procName[MPI_MAX_PROCESSOR_NAME];

      code = MPI_Recv (&procId, sizeof(INT4), MPI_BYTE, slave,
                       MPIMisc, MPI_COMM_WORLD, &mpiStatus);
      ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
      ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
              COMM_EMPIE, COMM_MSGEMPIE);

      code = MPI_Recv (procName, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, slave,
                       MPIMisc, MPI_COMM_WORLD, &mpiStatus);
      ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
      ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
              COMM_EMPIE, COMM_MSGEMPIE);

      fprintf (fp, "rsh %s kill -9 %d\n", procName, procId);
      fflush (fp);
    }

    fclose (fp);
  }
  else /* I am a slave */
  {
    code = MPI_Send (&procId, sizeof(INT4), MPI_BYTE, 0,
                     MPIMisc, MPI_COMM_WORLD);
    ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

    code = MPI_Send (id->procName, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, 0,
                     MPIMisc, MPI_COMM_WORLD);
    ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  }

  chmod ("killscript", S_IRWXU);

  RETURN (status);
}


void
LALMPISendMsg (
    LALStatus     *status,
    MPIMessage *msg,
    INT4        dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendMsg", COMMC);

  size = sizeof(MPIMessage);
  code = MPI_Send (msg, size, MPI_BYTE, dest, MPIMsg, MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}


void
LALMPIRecvMsg (
    LALStatus     *status,
    MPIMessage *msg
    )
{
  MPI_Status mpiStatus;
  INT4       code;
  INT4       size;

  INITSTATUS (status, "LALMPIRecvMsg", COMMC);

  size = sizeof(MPIMessage);
  code = MPI_Recv (msg, size, MPI_BYTE, MPI_ANY_SOURCE, MPIMsg,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}



/*
 *
 * Routines to send and receive vectors
 *
 */

/* Send/Recv CHARVector: */

void
LALMPISendCHARVector (
    LALStatus     *status,
    CHARVector *vector,
    INT4        dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendCHARVector", COMMC);

  ASSERT (vector, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* send the vector structure */
  size = sizeof(CHARVector);
  code = MPI_Send (vector, size, MPI_BYTE, dest, MPICHARVector,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  /* send the vector data */
  size = vector->length*sizeof(CHAR);
  code = MPI_Send (vector->data, size, MPI_BYTE, dest, MPICHARVectorData,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}


void
LALMPIRecvCHARVector (
    LALStatus     *status,
    CHARVector *vector,
    INT4        source
    )
{
  MPI_Status mpiStatus;
  CHARVector tmpVec;
  INT4       code;
  INT4       size;

  INITSTATUS (status, "LALMPIRecvCHARVector", COMMC);

  ASSERT (vector, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* receive temporary vector structure to check */
  size = sizeof(CHARVector);
  code = MPI_Recv (&tmpVec, size, MPI_BYTE, source, MPICHARVector,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  /* make sure vector lengths agree */
  ASSERT (vector->length == tmpVec.length, status,
          16, "Exchange size mismatch");

  /* receive the vector data */
  size = vector->length*sizeof(CHAR);
  code = MPI_Recv (vector->data, size, MPI_BYTE, source, MPICHARVectorData,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}


/* Send/Recv INT2Vector: */


void
LALMPISendINT2Vector (
    LALStatus     *status,
    INT2Vector *vector,
    INT4        dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendINT2Vector", COMMC);

  ASSERT (vector, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* send the vector structure */
  size = sizeof(INT2Vector);
  code = MPI_Send (vector, size, MPI_BYTE, dest, MPII2Vector, MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  /* send the vector data */
  size = vector->length*sizeof(INT2);
  code = MPI_Send (vector->data, size, MPI_BYTE, dest, MPII2VectorData,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}


void
LALMPIRecvINT2Vector (
    LALStatus     *status,
    INT2Vector *vector,
    INT4        source
    )
{
  MPI_Status mpiStatus;
  INT2Vector tmpVec;
  INT4       code;
  INT4       size;

  INITSTATUS (status, "LALMPIRecvINT2Vector", COMMC);

  ASSERT (vector, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* receive temporary vector structure to check */
  size = sizeof(INT2Vector);
  code = MPI_Recv (&tmpVec, size, MPI_BYTE, source, MPII2Vector,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  /* make sure vector lengths agree */
  ASSERT (vector->length == tmpVec.length, status,
          16, "Exchange size mismatch");

  /* receive the vector data */
  size = vector->length*sizeof(INT2);
  code = MPI_Recv (vector->data, size, MPI_BYTE, source, MPII2VectorData,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}


/* Send/Recv REAL4Vector: */


void
LALMPISendREAL4Vector (
    LALStatus      *status,
    REAL4Vector *vector,
    INT4         dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendREAL4Vector", COMMC);

  ASSERT (vector, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* send the vector structure */
  size = sizeof(REAL4Vector);
  code = MPI_Send (vector, size, MPI_BYTE, dest, MPISVector, MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  /* send the vector data */
  size = vector->length*sizeof(REAL4);
  code = MPI_Send (vector->data, size, MPI_BYTE, dest, MPISVectorData,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}


void
LALMPIRecvREAL4Vector (
    LALStatus      *status,
    REAL4Vector *vector,
    INT4         source
    )
{
  MPI_Status  mpiStatus;
  REAL4Vector tmpVec;
  INT4        code;
  INT4        size;

  INITSTATUS (status, "LALMPIRecvREAL4Vector", COMMC);

  ASSERT (vector, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* receive temporary vector structure to check */
  size = sizeof(REAL4Vector);
  code = MPI_Recv (&tmpVec, size, MPI_BYTE, source, MPISVector,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  /* make sure vector lengths agree */
  ASSERT (vector->length == tmpVec.length, status,
          16, "Exchange size mismatch");

  /* receive the vector data */
  size = vector->length*sizeof(REAL4);
  code = MPI_Recv (vector->data, size, MPI_BYTE, source, MPISVectorData,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}


/* Send/Recv COMPLEX8Vector: */


void
LALMPISendCOMPLEX8Vector (
    LALStatus         *status,
    COMPLEX8Vector *vector,
    INT4            dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendCOMPLEX8Vector", COMMC);

  ASSERT (vector, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* send the vector structure */
  size = sizeof(COMPLEX8Vector);
  code = MPI_Send (vector, size, MPI_BYTE, dest, MPICVector, MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  /* send the vector data */
  size = vector->length*sizeof(COMPLEX8);
  code = MPI_Send (vector->data, size, MPI_BYTE, dest, MPICVectorData,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}


void
LALMPIRecvCOMPLEX8Vector (
    LALStatus         *status,
    COMPLEX8Vector *vector,
    INT4            source
    )
{
  MPI_Status     mpiStatus;
  COMPLEX8Vector tmpVec;
  INT4           code;
  INT4           size;

  INITSTATUS (status, "LALMPIRecvCOMPLEX8Vector", COMMC);

  ASSERT (vector, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (vector->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* receive temporary vector structure to check */
  size = sizeof(COMPLEX8Vector);
  code = MPI_Recv (&tmpVec, size, MPI_BYTE, source, MPICVector,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  /* make sure vector lengths agree */
  ASSERT (vector->length == tmpVec.length, status,
          16, "Exchange size mismatch");

  /* receive the vector data */
  size = vector->length*sizeof(COMPLEX8);
  code = MPI_Recv (vector->data, size, MPI_BYTE, source, MPICVectorData,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  RETURN (status);
}


/*
 *
 * Routines to send and receive time series
 *
 */


/* Send/Recv INT2TimeSeries: */


void
LALMPISendINT2TimeSeries (
    LALStatus          *status,
    INT2TimeSeries *series,
    INT4             dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendINT2TimeSeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* send the series structure */
  size = sizeof(INT2TimeSeries);
  code = MPI_Send (series, size, MPI_BYTE, dest, MPII2TimeSeries,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  /* not sure how to send name --- ignore it! */

  /* send sampleUnits vector */
  LALMPISendCHARVector (status->statusPtr, series->sampleUnits, dest);
  CHECKSTATUSPTR (status);

  /* send data vector */
  LALMPISendINT2Vector (status->statusPtr, series->data, dest);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALMPIRecvINT2TimeSeries (
    LALStatus         *status,
    INT2TimeSeries *series,
    INT4            source
    )
{
  MPI_Status     mpiStatus;
  INT2TimeSeries tmpSer;
  INT4            code;
  INT4            size;

  INITSTATUS (status, "LALMPIRecvINT2TimeSeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* receive temporary time series structure */
  size = sizeof(INT2TimeSeries);
  code = MPI_Recv (&tmpSer, size, MPI_BYTE, source, MPII2TimeSeries,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  /* set the appropriate fields */
  series->epoch  = tmpSer.epoch;
  series->f0     = tmpSer.f0;
  series->deltaT = tmpSer.deltaT;

  /* don't know what to do with name... */
  series->name   = "anonymous";

  /* receive sampleUnits vector */
  LALMPIRecvCHARVector (status->statusPtr, series->sampleUnits, source);
  CHECKSTATUSPTR (status);

  /* receive data vector */
  LALMPIRecvINT2Vector (status->statusPtr, series->data, source);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/* Send/Recv REAL4TimeSeries: */


void
LALMPISendREAL4TimeSeries (
    LALStatus          *status,
    REAL4TimeSeries *series,
    INT4             dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendREAL4TimeSeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* send the series structure */
  size = sizeof(REAL4TimeSeries);
  code = MPI_Send (series, size, MPI_BYTE, dest, MPISTimeSeries,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  /* not sure how to send name --- ignore it! */

  /* send sampleUnits vector */
  LALMPISendCHARVector (status->statusPtr, series->sampleUnits, dest);
  CHECKSTATUSPTR (status);

  /* send data vector */
  LALMPISendREAL4Vector (status->statusPtr, series->data, dest);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALMPIRecvREAL4TimeSeries (
    LALStatus          *status,
    REAL4TimeSeries *series,
    INT4             source
    )
{
  MPI_Status      mpiStatus;
  REAL4TimeSeries tmpSer;
  INT4            code;
  INT4            size;

  INITSTATUS (status, "LALMPIRecvREAL4TimeSeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* receive temporary time series structure */
  size = sizeof(REAL4TimeSeries);
  code = MPI_Recv (&tmpSer, size, MPI_BYTE, source, MPISTimeSeries,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status, COMM_EMPIE,
          COMM_MSGEMPIE);

  /* set the appropriate fields */
  series->epoch  = tmpSer.epoch;
  series->f0     = tmpSer.f0;
  series->deltaT = tmpSer.deltaT;

  /* don't know what to do with name... */
  series->name   = "anonymous";

  /* receive sampleUnits vector */
  LALMPIRecvCHARVector (status->statusPtr, series->sampleUnits, source);
  CHECKSTATUSPTR (status);

  /* receive data vector */
  LALMPIRecvREAL4Vector (status->statusPtr, series->data, source);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* Send/Recv COMPLEX8TimeSeries: */


void
LALMPISendCOMPLEX8TimeSeries (
    LALStatus             *status,
    COMPLEX8TimeSeries *series,
    INT4                dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendCOMPLEX8TimeSeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* send the series structure */
  size = sizeof(COMPLEX8TimeSeries);
  code = MPI_Send (series, size, MPI_BYTE, dest, MPICTimeSeries,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  /* not sure how to send name --- ignore it! */

  /* send sampleUnits vector */
  LALMPISendCHARVector (status->statusPtr, series->sampleUnits, dest);
  CHECKSTATUSPTR (status);

  /* send data vector */
  LALMPISendCOMPLEX8Vector (status->statusPtr, series->data, dest);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALMPIRecvCOMPLEX8TimeSeries (
    LALStatus             *status,
    COMPLEX8TimeSeries *series,
    INT4                source
    )
{
  MPI_Status         mpiStatus;
  COMPLEX8TimeSeries tmpSer;
  INT4               code;
  INT4               size;

  INITSTATUS (status, "LALMPIRecvCOMPLEX8TimeSeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* receive temporary time series structure */
  size = sizeof(COMPLEX8TimeSeries);
  code = MPI_Recv (&tmpSer, size, MPI_BYTE, source, MPICTimeSeries,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  /* set the appropriate fields */
  series->epoch  = tmpSer.epoch;
  series->f0     = tmpSer.f0;
  series->deltaT = tmpSer.deltaT;

  /* don't know what to do with name... */
  series->name   = "anonymous";

  /* receive sampleUnits vector */
  LALMPIRecvCHARVector (status->statusPtr, series->sampleUnits, source);
  CHECKSTATUSPTR (status);

  /* receive data vector */
  LALMPIRecvCOMPLEX8Vector (status->statusPtr, series->data, source);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/*
 *
 * Routines to send and receive frequency series
 *
 */


/* Send/Recv REAL4FrequencySeries: */


void
LALMPISendREAL4FrequencySeries (
    LALStatus               *status,
    REAL4FrequencySeries *series,
    INT4                  dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendREAL4FrequencySeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* send the series structure */
  size = sizeof(REAL4FrequencySeries);
  code = MPI_Send (series, size, MPI_BYTE, dest, MPISFrequencySeries,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  /* not sure how to send name --- ignore it! */

  /* send sampleUnits vector */
  LALMPISendCHARVector (status->statusPtr, series->sampleUnits, dest);
  CHECKSTATUSPTR (status);

  /* send data vector */
  LALMPISendREAL4Vector (status->statusPtr, series->data, dest);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALMPIRecvREAL4FrequencySeries (
    LALStatus               *status,
    REAL4FrequencySeries *series,
    INT4                  source
    )
{
  MPI_Status           mpiStatus;
  REAL4FrequencySeries tmpSer;
  INT4                 code;
  INT4                 size;

  INITSTATUS (status, "LALMPIRecvREAL4FrequencySeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* receive temporary time series structure */
  size = sizeof(REAL4FrequencySeries);
  code = MPI_Recv (&tmpSer, size, MPI_BYTE, source, MPISFrequencySeries,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  /* set the appropriate fields */
  series->epoch  = tmpSer.epoch;
  series->f0     = tmpSer.f0;
  series->deltaF = tmpSer.deltaF;

  /* don't know what to do with name... */
  series->name   = "anonymous";

  /* receive sampleUnits vector */
  LALMPIRecvCHARVector (status->statusPtr, series->sampleUnits, source);
  CHECKSTATUSPTR (status);

  /* receive data vector */
  LALMPIRecvREAL4Vector (status->statusPtr, series->data, source);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* Send/Recv COMPLEX8FrequencySeries: */


void
LALMPISendCOMPLEX8FrequencySeries (
    LALStatus                  *status,
    COMPLEX8FrequencySeries *series,
    INT4                     dest
    )
{
  INT4 code;
  INT4 size;

  INITSTATUS (status, "LALMPISendCOMPLEX8FrequencySeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* send the series structure */
  size = sizeof(COMPLEX8FrequencySeries);
  code = MPI_Send (series, size, MPI_BYTE, dest, MPICFrequencySeries,
                   MPI_COMM_WORLD);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);

  /* not sure how to send name --- ignore it! */

  /* send sampleUnits vector */
  LALMPISendCHARVector (status->statusPtr, series->sampleUnits, dest);
  CHECKSTATUSPTR (status);

  /* send data vector */
  LALMPISendCOMPLEX8Vector (status->statusPtr, series->data, dest);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALMPIRecvCOMPLEX8FrequencySeries (
    LALStatus                  *status,
    COMPLEX8FrequencySeries *series,
    INT4                     source
    )
{
  MPI_Status              mpiStatus;
  COMPLEX8FrequencySeries tmpSer;
  INT4                    code;
  INT4                    size;

  INITSTATUS (status, "LALMPIRecvCOMPLEX8FrequencySeries", COMMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (series, status, COMM_ENULL, COMM_MSGENULL);
  /* ASSERT (series->name, status, 1, "Null pointer"); */
  ASSERT (series->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->data->length, status, COMM_ESIZE, COMM_MSGESIZE);
  ASSERT (series->sampleUnits, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->data, status, COMM_ENULL, COMM_MSGENULL);
  ASSERT (series->sampleUnits->length, status, COMM_ESIZE, COMM_MSGESIZE);

  /* receive temporary time series structure */
  size = sizeof(COMPLEX8FrequencySeries);
  code = MPI_Recv (&tmpSer, size, MPI_BYTE, source, MPICFrequencySeries,
                   MPI_COMM_WORLD, &mpiStatus);
  ASSERT (code == MPI_SUCCESS, status, COMM_EMPIE, COMM_MSGEMPIE);
  ASSERT (mpiStatus.MPI_ERROR == MPI_SUCCESS, status,
          COMM_EMPIE, COMM_MSGEMPIE);

  /* set the appropriate fields */
  series->epoch  = tmpSer.epoch;
  series->f0     = tmpSer.f0;
  series->deltaF = tmpSer.deltaF;

  /* don't know what to do with name... */
  series->name   = "anonymous";

  /* receive sampleUnits vector */
  LALMPIRecvCHARVector (status->statusPtr, series->sampleUnits, source);
  CHECKSTATUSPTR (status);

  /* receive data vector */
  LALMPIRecvCOMPLEX8Vector (status->statusPtr, series->data, source);
  CHECKSTATUSPTR (status);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
