/*----------------------------------------------------------------------- 
 * 
 * File Name: CommTestMaster.c
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

NRCSID (COMMTESTMASTERC, "$Id$");

/* global variables */
#include "CommTestGlobal.h"

void
Master (LALStatus *status, MPIId id)
{
  MPIMessage   message;
  REAL4Vector *vector = NULL;
  INT4         numProcs;
  UINT4        i;

  INITSTATUS (status, "Master", COMMTESTMASTERC);
  ATTATCHSTATUSPTR (status);

  printf ("Master starting up\n");

  LALSCreateVector (status->statusPtr, &vector, numPoints);
  CHECKSTATUSPTR (status);

  numProcs = id.numProcs;

  while (--numProcs > 0)
  {
    LALMPIRecvMsg (status->statusPtr, &message, MPI_COMM_WORLD);
    CHECKSTATUSPTR (status);

    printf ("Master received message code %d from slave %d\n",
            message.msg, message.source);

    switch (message.msg)
    {
      case MPISVector:

        LALMPIRecvREAL4Vector (status->statusPtr, vector, 
            message.source, MPI_COMM_WORLD);
        CHECKSTATUSPTR (status);

        ASSERT (message.send == 1, status, 99, "Master expects to receive");

        printf ("Master received REAL4Vector from slave %d\n",
                message.source);

        for (i = 0; i < vector->length; ++i)
        {
          printf ("%f\t%f\n", vector->data[i], (REAL4)(i % 5 - 2));
        }

        break;

      default:

        printf ("Master: unexpected message code!\n");
        ++numProcs;
    }

  }

  LALSDestroyVector (status->statusPtr, &vector);
  CHECKSTATUSPTR (status);

  printf ("Master shutting down\n");

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
