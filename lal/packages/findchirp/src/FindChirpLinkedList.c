/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpLinkedList.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/FindChirpEngine.h>

NRCSID (FINDCHIRPLINKEDLISTC, "$Id$");

void
LALFindChirpCreateTmpltNode (
    LALStatus                  *status,
    InspiralTemplate           *tmplt,
    InspiralTemplateNode      **tmpltNode
    )
{
  InspiralTemplateNode         *current = NULL;

  INITSTATUS( status, "LALFindChirpCreateTmpltNode", FINDCHIRPLINKEDLISTC );
  ATTATCHSTATUSPTR( status );

  ASSERT( tmplt, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );


  /*
   *
   * create a new template node after the current one
   *
   */


  /* store the address of the current node */
  current = *tmpltNode;

  /* create memory for the template node */
  *tmpltNode = (InspiralTemplateNode *) 
    LALCalloc( 1, sizeof(InspiralTemplateNode) );
  if ( !tmpltNode )
  {
    ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
  }

  (*tmpltNode)->prev     = NULL;
  (*tmpltNode)->next     = NULL;
  (*tmpltNode)->tmpltPtr = tmplt;

  /* link the list */
  if ( current ) 
  {
    if ( current->next )
    {
      (*tmpltNode)->next = current->next;
      current->next->prev = *tmpltNode;
    }
    (*tmpltNode)->prev = current;
    current->next = *tmpltNode;
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}  


void
LALFindChirpDestroyTmpltNode ( 
    LALStatus                  *status,
    InspiralTemplateNode      **tmpltNode
    )
{
  InspiralTemplateNode  *prev = NULL;
  InspiralTemplateNode  *next = NULL;

  INITSTATUS( status, "FindChirpDestroyTmpltNode", FINDCHIRPLINKEDLISTC );
  ATTATCHSTATUSPTR( status );

  ASSERT( tmpltNode, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );


  /*
   *
   * destroy the node pointed at and return a pointer to the node before
   *
   */

  
  /* store the previous and next nodes */
  prev = (*tmpltNode)->prev;
  next = (*tmpltNode)->next;

  LALFree( *tmpltNode );
  *tmpltNode = NULL;

  /* relink the list */
  if ( next ) next->prev = prev;
  if ( prev ) 
  {
    prev->next = next;
    *tmpltNode = prev;
  }
  else
  {
    *tmpltNode = next;
  }
  
  DETATCHSTATUSPTR( status );
  RETURN( status );
}  
