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

#ifndef _STDIO_H
#include <stdio.h>
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _COMM_H
#include "Comm.h"
#ifndef _COMM_H
#define _COMM_H
#endif
#endif

NRCSID (COMMTESTSLAVEC, "$Id$");

/* global variables */
#include "CommTestGlobal.h"

void
Slave (Status *status, MPIId id)
{
  MPIMessage   message;
  REAL4Vector *vector = NULL;
  INT4         i;

  INITSTATUS (status, COMMTESTSLAVEC);
  ATTATCHSTATUSPTR (status);

  printf ("Slave %d starting up\n", id.myId);

  SCreateVector (status->statusPtr, &vector, numPoints);
  CHECKSTATUSPTR (status);

  printf ("Slave %d sending message code %d to master\n", id.myId, MPISVector);
  message.msg    = MPISVector;
  message.send   = 1;
  message.source = id.myId;
  MPISendMsg (status->statusPtr, &message, 0);
  CHECKSTATUSPTR (status);

  for (i = 0; i < vector->length; ++i)
  {
    vector->data[i] = i % 5 - 2;
  }

  printf ("Slave %d sending REAL4Vector to master\n", id.myId);
  MPISendREAL4Vector (status->statusPtr, vector, 0);
  CHECKSTATUSPTR (status);

  SDestroyVector (status->statusPtr, &vector);
  CHECKSTATUSPTR (status);

  printf ("Slave %d shutting down\n", id.myId);
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
