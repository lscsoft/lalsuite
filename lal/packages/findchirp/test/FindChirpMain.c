/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpMain.c
 *
 * Author: Allen, B., Brown, D., and Creighton, J. D. E.
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

#ifndef _MPI_INCLUDE
#include "mpi.h"
#ifndef _MPI_INCLUDE
#define _MPI_INCLUDE
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

NRCSID (MAIN, "$Id$");

void
Master (Status *, MPIId);

void
Slave (Status *, MPIId);

/* global variables */
#define FINDCHIRPGLOBAL_INIT
#include "FindChirpGlobal.h"
#undef FINDCHIRPGLOBAL_INIT

int debuglevel = 1;

int
main (int argc, char *argv[])
{
  static Status  stat;
  MPIDebugParams debug;
  MPIId          id;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &id.numProcs);
  MPI_Comm_rank (MPI_COMM_WORLD, &id.myId);
  MPI_Get_processor_name (id.procName, &id.nameLen);

  /* write kill script */
  MPIKillScript (&stat, &id);

  /* debugging example */
  if (0)
  {
    MPIExportEnvironment (&stat, "DISPLAY", id.myId);
    debug.debugger = id.myId == 0 ? "ddd" : NULL ; /* attatch ddd to master */
    debug.progName = argv[0];
    debug.delay    = 15;
    debug.myId     = id.myId;
    MPIDebug (&stat, &debug);
  }

  if (id.myId == 0)
  {
    Master (&stat, id);
    REPORTSTATUS (&stat);
  }
  else
  {
    char fname[64];
    sprintf (fname, "Slave%03d.err", id.myId);
    freopen (fname, "w", stderr);
    sprintf (fname, "Slave%03d.out", id.myId);
    freopen (fname, "w", stdout);
    Slave (&stat, id);
    REPORTSTATUS (&stat);
  }

  MPI_Finalize ();

  LALCheckMemoryLeaks ();
  return 0;
}

