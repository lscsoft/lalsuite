/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChripTmpltTest.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/FindChirpEngine.h>

NRCSID (MAIN, "$Id$");

int lalDebugLevel = 1;

static
void PrintInspiralBank (
    LALStatus           *status,
    InspiralTemplate    *head,
    FILE                *fp
    );

int
main (int argc, char *argv[])
{
  static LALStatus  stat;

  InspiralCoarseBankIn          *bankIn = NULL;
  InspiralTemplate              *tmplt = NULL;

  LALFindChirpCreateInspiralBank( &stat, bankIn, &tmplt );
  REPORTSTATUS( &stat );
  if ( stat.statusCode ) return stat.statusCode;

  LALCheckMemoryLeaks ();

  return 0;
}

static
void
PrintInspiralBank (
    LALStatus           *status,
    InspiralTemplate    *head,
    FILE                *fp
    )
{
  InspiralTemplate    *current;

  INITSTATUS( status, "PrintInspiralBank", MAIN );
  ATTATCHSTATUSPTR( status );

  ASSERT( head, status, 1, "Null Pointer" );


  /*
   *
   * print the inspiral template parameter bank
   *
   */


  for ( current = head; current; current = current->next)
  {
    fprintf( fp, 
        "%p\t%d\t%f\t%f\t%f\t%f\t%f\t%p\t%p\n",
        current, current->number, 
        current->mass1, current->mass2, current->totalMass,
        current->eta, current->mu,
        current->next, current->fine );
    fflush( fp );
    if ( current->fine )
    {
      PrintInspiralBank( status->statusPtr, current->fine, fp );
      CHECKSTATUSPTR( status );
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
