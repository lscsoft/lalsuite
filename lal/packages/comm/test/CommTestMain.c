/*----------------------------------------------------------------------- 
 * 
 * File Name: CommTestMain.c
 *
 * Author: Allen, B. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <mpi.h>
#include <lal/LALStdlib.h>
#include <lal/Comm.h>

NRCSID (MAIN, "$Id$");

/* global variables */
#define COMMTESTGLOBAL_INIT
#include "CommTestGlobal.h"
#undef COMMTESTGLOBAL_INIT

void
Master (LALStatus *, MPIId);

void
Slave (LALStatus *, MPIId);

int lalDebugLevel = 1;

int
main (int argc, char *argv[])
{
  static LALStatus  stat;
  MPIDebugParams debug;
  MPIId          id;

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &id.numProcs);
  MPI_Comm_rank (MPI_COMM_WORLD, &id.myId);
  MPI_Get_processor_name (id.procName, &id.nameLen);

  /* write kill script */
  LALMPIKillScript (&stat, &id, MPI_COMM_WORLD);

  /* debugging example */
  if (0)
  {
    LALMPIExportEnvironment (&stat, "DISPLAY", id.myId, MPI_COMM_WORLD);
    debug.debugger = id.myId == 0 ? "ddd" : NULL ; /* attatch ddd to master */
    debug.progName = argv[0];
    debug.delay    = 15;
    debug.myId     = id.myId;
    debug.mpiComm  = MPI_COMM_WORLD;
    LALMPIDebug (&stat, &debug);
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

