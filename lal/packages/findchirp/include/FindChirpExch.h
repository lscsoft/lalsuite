/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpExch.h
 *
 * Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpExchHV">
Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
$Id$
</lalVerbatim> 

<lalLaTeX>

\section{Header \texttt{FindChirpExch.h}}
\label{s:FindChirpExch.h}

Provides routines for MPI communication of various \texttt{findchirp} specific
data types.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpExch.h>
\end{verbatim}

</lalLaTeX>
#endif

#ifndef _FINDCHIRPEXCHH_H
#define _FINDCHIRPEXCHH_H

#include <lal/LALDatatypes.h>
#include <lal/Comm.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (FINDCHIRPEXCHH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPEXCHH_ENULL 1
#define FINDCHIRPEXCHH_ENNUL 2
#define FINDCHIRPEXCHH_ENOBJ 4
#define FINDCHIRPEXCHH_EHAND 8
#define FINDCHIRPEXCHH_EMPIE 16
#define FINDCHIRPEXCHH_MSGENULL "Null pointer"
#define FINDCHIRPEXCHH_MSGENNUL "Non-null pointer"
#define FINDCHIRPEXCHH_MSGENOBJ "Invalid number of objects"
#define FINDCHIRPEXCHH_MSGEHAND "Wrong handshake"
#define FINDCHIRPEXCHH_MSGEMPIE "Problem exchanging event list"
/* </lalErrTable> */

#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpExchHV}}
</lalLaTeX> 
#endif

#if 0
<lalLaTeX>
\newpage\input{FindChirpExchC}
</lalLaTeX>
#endif

void
LALExchangeDataSegment (
    LALStatus      *status,
    DataSegment *segment,
    ExchParams  *exchParams
    );

void
LALExchangeInspiralTemplate (
    LALStatus           *status,
    InspiralTemplate *tmplt,
    ExchParams       *exchParams
    );

void
LALExchangeInspiralEvent (
    LALStatus        *status,
    InspiralEvent *event,
    ExchParams    *exchParams
    );

void
LALExchangeInspiralEventList (
    LALStatus     *status,
    InspiralEvent **eventHead,
    ExchParams    *exchParams
    );

void
LALExchangeTemplateBank (
    LALStatus         *status,
    InspiralTemplate **tmpltHead,
    ExchParams        *exchParms
                 );


#ifdef  __cplusplus
}
#endif

#endif /* _FINDCHIRPEXCHH_H */
