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

NRCSID (COMMTESTMASTERC, "$Id$");

/* global variables */
#include "CommTestGlobal.h"

void
Master (Status *status, MPIId id)
{
  MPIMessage   message;
  REAL4Vector *vector = NULL;
  INT4         numProcs;
  INT4         i;

  INITSTATUS (status, COMMTESTMASTERC);
  ATTATCHSTATUSPTR (status);

  printf ("Master starting up\n");

  SCreateVector (status->statusPtr, &vector, numPoints);
  CHECKSTATUSPTR (status);

  numProcs = id.numProcs;

  while (--numProcs > 0)
  {
    MPIRecvMsg (status->statusPtr, &message);
    CHECKSTATUSPTR (status);

    printf ("Master received message code %d from slave %d\n",
            message.msg, message.source);

    switch (message.msg)
    {
      case MPISVector:

        MPIRecvREAL4Vector (status->statusPtr, vector, message.source);
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

  SDestroyVector (status->statusPtr, &vector);
  CHECKSTATUSPTR (status);

  printf ("Master shutting down\n");

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
