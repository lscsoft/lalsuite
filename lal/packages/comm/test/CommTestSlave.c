/*----------------------------------------------------------------------- 
 * 
 * File Name: CommTestSlave.c
 *
 * Author: Allen, B. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Comm.h>

NRCSID (COMMTESTSLAVEC, "$Id$");

/* global variables */
#include "CommTestGlobal.h"

void
Slave (LALStatus *status, MPIId id)
{
  MPIMessage   message;
  REAL4Vector *vector = NULL;
  INT4         i;

  INITSTATUS (status, "Slave", COMMTESTSLAVEC);
  ATTATCHSTATUSPTR (status);

  printf ("Slave %d starting up\n", id.myId);

  LALSCreateVector (status->statusPtr, &vector, numPoints);
  CHECKSTATUSPTR (status);

  printf ("Slave %d sending message code %d to master\n", id.myId, MPISVector);
  message.msg    = MPISVector;
  message.send   = 1;
  message.source = id.myId;
  LALMPISendMsg (status->statusPtr, &message, 0, MPI_COMM_WORLD);
  CHECKSTATUSPTR (status);

  for (i = 0; i < (INT4) vector->length; ++i)
  {
    vector->data[i] = i % 5 - 2;
  }

  printf ("Slave %d sending REAL4Vector to master\n", id.myId);
  LALMPISendREAL4Vector (status->statusPtr, vector, 0, MPI_COMM_WORLD);
  CHECKSTATUSPTR (status);

  LALSDestroyVector (status->statusPtr, &vector);
  CHECKSTATUSPTR (status);

  printf ("Slave %d shutting down\n", id.myId);
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
