/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpEngine.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpEngineHV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{FindChirpEngine.h}}
\label{s:FindChirpEngine.h}

Provides functions to drive the filtering functions in the package
\texttt{findchirp}. 

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpEngine.h>
\end{verbatim}

\input{FindChirpEngineHDoc}

</lalLaTeX>
#endif

#ifndef _FINDCHIRPENGINEH_H
#define _FINDCHIRPENGINEH_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <lal/LALStdlib.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPENGINEHH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPENGINEH_ENULL 1
#define FINDCHIRPENGINEH_ENNUL 2
#define FINDCHIRPENGINEH_EALOC 3
#define FINDCHIRPENGINEH_ELVEL 4
#define FINDCHIRPENGINEH_MSGENULL "Null pointer"
#define FINDCHIRPENGINEH_MSGENNUL "Non-null pointer"
#define FINDCHIRPENGINEH_MSGEALOC "Memory allocation error"
#define FINDCHIRPENGINEH_MSGELVEL "Invalid heriarchical template bank level"
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif
/* --- structure for managing a list of inspiral templates --------------- */
/* <lalVerbatim file="FindChirpEngineHInspiralTemplateNode"> */
typedef struct
tagInspiralTemplateNode
{
  struct tagInspiralTemplateNode       *next;
  struct tagInspiralTemplateNode       *prev;
  InspiralTemplate                     *tmpltPtr;
}
InspiralTemplateNode;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{InspiralTemplateNode}}
\idx[Type]{InspiralTemplateNode}

\input{FindChirpEngineHInspiralTemplateNode}

\noindent This structure provides a method of constucting doubly linked
lists of \texttt{InspiralTemplate} structures. The fields are:

\begin{description}
\item[\texttt{struct tagInspiralTemplateNode *next}] The next structure in
the linked list.

\item[\texttt{struct tagInspiralTemplateNode *prev}] The previous structure in
the linked list.

\item[\texttt{InspiralTemplate *tmpltPtr}] A pointer to an \texttt{InspiralTemplate} structure.
\end{description}
</lalLaTeX>
#endif

void
LALFindChirpCreateTmpltNode (
    LALStatus                  *status,
    InspiralTemplate           *tmplt,
    InspiralTemplateNode      **tmpltNode
    );

void
LALFindChirpDestroyTmpltNode ( 
    LALStatus                  *status,
    InspiralTemplateNode      **tmpltNode
    );

#ifdef  __cplusplus
}
#endif

#endif /* _FINDCHIRPENGINEH_H */
