/****************************** <lalVerbatim file="CreateZPGFilterCV">
Author: Creighton, T. D.
$Id$
******************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{CreateZPGFilter.c}}
\label{ss:CreateZPGFilter.c}

Creates ZPG filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CreateZPGFilterCP}
\index{\verb&CreateCOMPLEX8ZPGFilter()&}
\index{\verb&CreateCOMPLEX16ZPGFilter()&}

\subsubsection*{Description}

These functions create an object \verb@**output@, of type
\verb@COMPLEX8ZPGFilter@ or \verb@COMPLEX16ZPGFilter@, having
\verb@numZeros@ zeros and \verb@numPoles@ poles.  The values of those
zeros and poles are not set by these routines (in general they will
start out as garbage).  The handle passed into the functions must be a
valid handle (i.e.\ \verb@output@$\neq$\verb@NULL@), but must not
point to an existing object (\i.e.\ \verb@*output@=\verb@NULL@).

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()
CCreateVector()
ZCreateVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CreateZPGFilterCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include "AVFactories.h"
#include "ZPGFilter.h"

NRCSID(CREATEZPGFILTERC,"$Id$");


/* <lalVerbatim file="CreateZPGFilterCP"> */
void CreateCOMPLEX8ZPGFilter(Status            *stat,
			     COMPLEX8ZPGFilter **output,
			     INT4              numZeros,
			     INT4              numPoles)
{ /* </lalVerbatim> */
  INITSTATUS(stat,"CreateCOMPLEX8ZPGFilter",CREATEZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);
  ASSERT(!*output,stat,ZPGFILTER_EOUT,ZPGFILTER_MSGEOUT);

  /* Make sure that numZeros and numPoles are non-negative. */
  ASSERT(numZeros>=0,stat,ZPGFILTER_EBAD,ZPGFILTER_MSGEBAD);
  ASSERT(numPoles>=0,stat,ZPGFILTER_EBAD,ZPGFILTER_MSGEBAD);

  /* Create the output structure. */
  *output=(COMPLEX8ZPGFilter *)LALMalloc(sizeof(COMPLEX8ZPGFilter));
  ASSERT(*output,stat,ZPGFILTER_EMEM,ZPGFILTER_MSGEMEM);
  memset(*output,0,sizeof(COMPLEX8ZPGFilter));

  /* Allocate the data fields.  If the number of poles or zeros is 0,
     the corresponding field(s) should remain null. */
  if(numZeros>0)
    TRY(CCreateVector(stat->statusPtr,&((*output)->zeros),numZeros),
	stat);
  if(numPoles>0)
    TRY(CCreateVector(stat->statusPtr,&((*output)->poles),numPoles),
	stat);

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="CreateZPGFilterCP"> */
void CreateCOMPLEX16ZPGFilter(Status             *stat,
			      COMPLEX16ZPGFilter **output,
			      INT4               numZeros,
			      INT4               numPoles)
{ /* </lalVerbatim> */
  INITSTATUS(stat,"CreateCOMPLEX16ZPGFilter",CREATEZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);
  ASSERT(!*output,stat,ZPGFILTER_EOUT,ZPGFILTER_MSGEOUT);

  /* Make sure that numZeros and numPoles are non-negative. */
  ASSERT(numZeros>=0,stat,ZPGFILTER_EBAD,ZPGFILTER_MSGEBAD);
  ASSERT(numPoles>=0,stat,ZPGFILTER_EBAD,ZPGFILTER_MSGEBAD);

  /* Create the output structure. */
  *output=(COMPLEX16ZPGFilter *)LALMalloc(sizeof(COMPLEX16ZPGFilter));
  ASSERT(*output,stat,ZPGFILTER_EMEM,ZPGFILTER_MSGEMEM);
  memset(*output,0,sizeof(COMPLEX16ZPGFilter));

  /* Allocate the data fields.  If the number of poles or zeros is 0,
     the corresponding field(s) should remain null. */
  if(numZeros>0)
    TRY(ZCreateVector(stat->statusPtr,&((*output)->zeros),numZeros),
	stat);
  if(numPoles>0)
    TRY(ZCreateVector(stat->statusPtr,&((*output)->poles),numPoles),
	stat);

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
