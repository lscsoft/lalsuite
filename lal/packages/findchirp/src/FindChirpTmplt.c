/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpTmplt.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/FindChirpEngine.h>

NRCSID( FINDCHIRPTMPLTC, "$Id$" );


void
LALFindChirpCreateInspiralBank (
    LALStatus                  *status,
    InspiralCoarseBankIn       *coarseBankIn,
    InspiralTemplate          **head
    )
{

  INITSTATUS( status, "LALFindChirpCreateInspiralBank", FINDCHIRPTMPLTC );
  ATTATCHSTATUSPTR( status );


  DETATCHSTATUSPTR( status );
  RETURN( status );
}


void
LALFindChirpDestroyInspiralBank (
    LALStatus           *status,
    InspiralTemplate   **head
    )
{
  InspiralTemplate    *current;

  INITSTATUS( status, "LALFindChirpDestroyInspiralBank", FINDCHIRPTMPLTC );
  ATTATCHSTATUSPTR( status );

  ASSERT( *head, status, FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );

  /* destroy the inspiral template parameter bank */
  while ( *head )
  {
    current = *head;
    if ( current->fine )
    {
      LALFindChirpDestroyInspiralBank( status->statusPtr, &(current->fine) );
      CHECKSTATUSPTR( status );
    }
    *head = (*head)->next;
    LALFree( current );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

