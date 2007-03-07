/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpPTFTemplate.c
 *
 * Author: Brown, D. A., and Fazi, D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpPTFTemplateCV">
Author: Brown, D. A., and Fazi, D.
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpPTFTemplate.c}}
\label{ss:FindChirpPTFTemplate.c}

Provides functions to create physical template family templates in a
form that can be used by the \texttt{FindChirpPTFFilter()} function.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpPTFTemplateCP}
\idx{LALFindChirpPTFTemplate()}

The function \texttt{LALFindChirpPTFTemplate()} creates a physical template
family template as described by the algorithm below.

\subsubsection*{Algorithm}

Blah.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpPTFTemplateCV}}
</lalLaTeX> 
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

NRCSID(FINDCHIRPPTFTEMPLATEC, "$Id$");

/* <lalVerbatim file="FindChirpPTFTemplateCP"> */
void
LALFindChirpPTFTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALFindChirpPTFTemplate", FINDCHIRPPTFTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  /* check that the output structures exist */
  ASSERT( fcTmplt, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
#if 0
  ASSERT( fcTmplt->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( fcTmplt->data->data, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
#endif

  /* check that the parameter structure exists */
  ASSERT( params, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->xfacVec, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->xfacVec->data, status, 
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the timestep is positive */
  ASSERT( params->deltaT > 0, status, 
      FINDCHIRPSPH_EDELT, FINDCHIRPSPH_MSGEDELT );

  /* check that the input exists */
  ASSERT( tmplt, status, FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  if ( params->approximant != FindChirpSP )
  {
    ABORT( status, FINDCHIRPSPH_EMAPX, FINDCHIRPSPH_MSGEMAPX );
  }
  LALInfo( status, "Generating template using FindChirpSP" );


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

