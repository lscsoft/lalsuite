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

#include <stdio.h>
#include <mpi.h>
#include <lal/LALStdlib.h>
#include <lal/Comm.h>

NRCSID (MAIN, "$Id$");

/*
 *
 * removed test programe until new engine is complete
 *
 */

int lalDebugLevel = 0;

int
main (int argc, char *argv[])
{
  MPI_Init( &argc, &argv );
  MPI_Finalize();

  return 0;
}

