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

Memory management functions for creating and destroying linked
lists of inspiral template nodes for flat and heirarchical search management.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpLinkedListCP}
\idx{LALFindChirpCreateTmpltNode()}
\idx{LALFindChirpDestroyTmpltNode()}

\subsubsection*{Description}

Placeholder.

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
