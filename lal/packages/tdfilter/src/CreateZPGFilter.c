/****************************** <lalVerbatim file="CreateZPGFilterCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{CreateZPGFilter.c}}
\label{ss:CreateZPGFilter.c}

Creates ZPG filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CreateZPGFilterCP}
\index{\texttt{LALCreateCOMPLEX8ZPGFilter()}}
\index{\texttt{LALCreateCOMPLEX16ZPGFilter()}}

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
LALMalloc()                     LALFree()
LALCCreateVector()              LALCDestroyVector()
LALZCreateVector()              LALZDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CreateZPGFilterCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ZPGFilter.h>

NRCSID(CREATEZPGFILTERC,"$Id$");


/* <lalVerbatim file="CreateZPGFilterCP"> */
void
LALCreateCOMPLEX8ZPGFilter( LALStatus         *stat,
			    COMPLEX8ZPGFilter **output,
			    INT4              numZeros,
			    INT4              numPoles )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALCreateCOMPLEX8ZPGFilter",CREATEZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(!*output,stat,ZPGFILTERH_EOUT,ZPGFILTERH_MSGEOUT);

  /* Make sure that numZeros and numPoles are non-negative. */
  ASSERT(numZeros>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);
  ASSERT(numPoles>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);

  /* Create the output structure. */
  *output=(COMPLEX8ZPGFilter *)LALMalloc(sizeof(COMPLEX8ZPGFilter));
  if ( !(*output) ) {
    ABORT(stat,ZPGFILTERH_EMEM,ZPGFILTERH_MSGEMEM);
  }
  memset(*output,0,sizeof(COMPLEX8ZPGFilter));

  /* Allocate the data fields.  If the number of poles or zeros is 0,
     the corresponding field(s) should remain null. */
  if(numZeros>0){
    LALCCreateVector(stat->statusPtr,&((*output)->zeros),numZeros);
    BEGINFAIL(stat) {
      LALFree(*output);
      *output=NULL;
    } ENDFAIL(stat);
  }
  if(numPoles>0){
    LALCCreateVector(stat->statusPtr,&((*output)->poles),numPoles);
    BEGINFAIL(stat) {
      TRY(LALCDestroyVector(stat->statusPtr,&((*output)->zeros)),
	  stat);
      LALFree(*output);
      *output=NULL;
    } ENDFAIL(stat);
  }

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="CreateZPGFilterCP"> */
void
LALCreateCOMPLEX16ZPGFilter( LALStatus          *stat,
			     COMPLEX16ZPGFilter **output,
			     INT4               numZeros,
			     INT4               numPoles )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALCreateCOMPLEX16ZPGFilter",CREATEZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(!*output,stat,ZPGFILTERH_EOUT,ZPGFILTERH_MSGEOUT);

  /* Make sure that numZeros and numPoles are non-negative. */
  ASSERT(numZeros>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);
  ASSERT(numPoles>=0,stat,ZPGFILTERH_EBAD,ZPGFILTERH_MSGEBAD);

  /* Create the output structure. */
  *output=(COMPLEX16ZPGFilter *)LALMalloc(sizeof(COMPLEX16ZPGFilter));
  if ( !(*output) ) {
    ABORT(stat,ZPGFILTERH_EMEM,ZPGFILTERH_MSGEMEM);
  }
  memset(*output,0,sizeof(COMPLEX16ZPGFilter));

  /* Allocate the data fields.  If the number of poles or zeros is 0,
     the corresponding field(s) should remain null. */
  if(numZeros>0){
    LALZCreateVector(stat->statusPtr,&((*output)->zeros),numZeros);
    BEGINFAIL(stat) {
      LALFree(*output);
      *output=NULL;
    } ENDFAIL(stat);
  }
  if(numPoles>0){
    LALZCreateVector(stat->statusPtr,&((*output)->poles),numPoles);
    BEGINFAIL(stat) {
      TRY(LALZDestroyVector(stat->statusPtr,&((*output)->zeros)),
	  stat);
      LALFree(*output);
      *output=NULL;
    } ENDFAIL(stat);
  }

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
