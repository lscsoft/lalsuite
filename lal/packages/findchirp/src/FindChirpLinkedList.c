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

#if 0
<lalVerbatim file="FindChirpLinkedListCV">
Author: Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpLinkedList.c}}
\label{ss:FindChirpLinkedList.c}

It is often convienient to deal with the inspiral templates as a doubly linked
list.  This module provides memory management functions for creating and
destroying linked lists of inspiral template nodes for flat and heirarchical
search management.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpLinkedListCP}
\idx{LALFindChirpCreateTmpltNode()}
\idx{LALFindChirpDestroyTmpltNode()}

\subsubsection*{Description}

The function \texttt{LALFindChirpCreateTmpltNode()} adds the inspiral template
parameter structure pointed to by \texttt{tmplt} to the linked list of
template nodes \texttt{tmpltNode}. On entry \texttt{tmpltNode} should be set
to memory address of the last node of the current linked list (or NULL if it
is a new linked list) and on exit \texttt{tmpltNode} is set to the memory
address of the last node in the linked list.

The function \texttt{LALFindChirpDestroyTmpltNode()} removed the node pointed
to by \texttt{tmpltNode} from the doubly linked list. On exit
\texttt{tmpltNode} is set to the address of the previous node in the list for
removal of a node in the middle or at the end of the list. If the first node
is removed \texttt{tmpltNode} is set to the address of the new first node.

\subsubsection*{Algorithm}

None.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpLinkedListCV}}
</lalLaTeX>
#endif


#include <lal/FindChirp.h>

NRCSID (FINDCHIRPLINKEDLISTC, "$Id$");

/* <lalVerbatim file="FindChirpLinkedListCP"> */
void
LALFindChirpCreateTmpltNode (
    LALStatus                  *status,
    InspiralTemplate           *tmplt,
    InspiralTemplateNode      **tmpltNode
    )
/* </lalVerbatim> */
{
  InspiralTemplateNode         *current = NULL;

  INITSTATUS( status, "LALFindChirpCreateTmpltNode", FINDCHIRPLINKEDLISTC );
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


/* <lalVerbatim file="FindChirpLinkedListCP"> */
void
LALFindChirpDestroyTmpltNode ( 
    LALStatus                  *status,
    InspiralTemplateNode      **tmpltNode
    )
/* </lalVerbatim> */
{
  InspiralTemplateNode  *prev = NULL;
  InspiralTemplateNode  *next = NULL;

  INITSTATUS( status, "FindChirpDestroyTmpltNode", FINDCHIRPLINKEDLISTC );
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
