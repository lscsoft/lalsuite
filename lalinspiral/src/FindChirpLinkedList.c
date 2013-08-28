/*
*  Copyright (C) 2007 Duncan Brown
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpLinkedList.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown D. A.
 * \file
 * \ingroup FindChirp_h
 *
 * \brief This module provides memory management functions for creating and
 * destroying linked lists of inspiral template nodes for flat and heirarchical
 * search management.
 *
 * It is often convienient to deal with the inspiral templates as a doubly linked
 * list.
 *
 * ### Description ###
 *
 * The function <tt>LALFindChirpCreateTmpltNode()</tt> adds the inspiral template
 * parameter structure pointed to by \c tmplt to the linked list of
 * template nodes \c tmpltNode. On entry \c tmpltNode should be set
 * to memory address of the last node of the current linked list (or NULL if it
 * is a new linked list) and on exit \c tmpltNode is set to the memory
 * address of the last node in the linked list.
 *
 * The function <tt>LALFindChirpDestroyTmpltNode()</tt> removed the node pointed
 * to by \c tmpltNode from the doubly linked list. On exit
 * \c tmpltNode is set to the address of the previous node in the list for
 * removal of a node in the middle or at the end of the list. If the first node
 * is removed \c tmpltNode is set to the address of the new first node.
 *
 * ### Uses ###
 *
 * \code
 * LALCalloc()
 * LALFree()
 * \endcode
 *
 */


#include <lal/FindChirp.h>

void
LALFindChirpCreateTmpltNode (
    LALStatus                  *status,
    InspiralTemplate           *tmplt,
    InspiralTemplateNode      **tmpltNode
    )

{
  InspiralTemplateNode         *current = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( tmplt, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );


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
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
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

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( tmpltNode, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );


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
