/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpTemplate.c
 *
 * Author: Brown D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpTemplateCV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpTemplate.c}}
\label{ss:FindChirpTemplat.c}

Provides functions to initialize template creation routines.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpTemplateCP}
\idx{LALFindChirpTemplateInit()}
\idx{LALFindChirpTemplateFinalize()}

The function \texttt{LALFindChirpTemplateInit()} takes as input the address
of a structure of type \texttt{FindChirpInitParams} containing the correct
values to intialize a search. It creates a structure of type
\texttt{FindChirpTmpltParams} as described above and returns its address.

The function \texttt{LALFindChirpTemplateFinalize()} takes as the address
of a structure of type \texttt{FindChirpTmpltParams} destroys this 
structure and sets the address to NULL.

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

\vfill{\footnotesize\input{FindChirpTemplateCV}}
</lalLaTeX> 
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>


NRCSID (FINDCHIRPTEMPLATEC, "$Id$");

/* <lalVerbatim file="FindChirpTemplateCP"> */
void
LALFindChirpTemplateInit (
    LALStatus                  *status,
    FindChirpTmpltParams      **output,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
{
  UINT4                         k;
  FindChirpTmpltParams         *outputPtr;
  REAL4                        *xfac = NULL;
  const REAL4                   exponent = -1.0/3.0;

  INITSTATUS( status, "LALFindChirpTemplateInit", FINDCHIRPTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( output, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );  
  ASSERT( !*output, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  
  /* make sure that the number of points is positive */
  ASSERT( params->numPoints > 0, status, 
      FINDCHIRPH_ENUMZ, FINDCHIRPH_MSGENUMZ );


  /*
   *
   * create tmplt generation parameters structure
   *
   */


  /* create the output structure */
  outputPtr = *output = (FindChirpTmpltParams *)
    LALCalloc( 1, sizeof(FindChirpTmpltParams) );
  if ( ! outputPtr )
  {
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* store the waveform approximant */
  outputPtr->approximant = params->approximant;

  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case GeneratePPN:
    case PadeT1:
    case EOB:
      /* time domain waveforms use xfac to store the time domain waveform */
      LALCreateVector( status->statusPtr, &(outputPtr->xfacVec), 
          params->numPoints );
      BEGINFAIL( status )
      {
        LALFree( outputPtr );
        *output = NULL;
      }
      ENDFAIL( status );
      
      memset( outputPtr->xfacVec->data, 0, 
          outputPtr->xfacVec->length * sizeof(REAL4) );
      
      /* create an fft plan for the time domain waveform */
      LALCreateForwardRealFFTPlan( status->statusPtr, &(outputPtr->fwdPlan), 
          params->numPoints, 0 );
      BEGINFAIL( status )
      {
        TRY( LALDestroyVector( status->statusPtr, &(outputPtr->xfacVec) ), 
            status );

        LALFree( outputPtr );
        *output = NULL;
      }
      ENDFAIL( status );
      break;
      
    case TaylorF2:
    case BCV:
    case BCVSpin:
      /* freq domain waveforms need xfac vector containing k^(-7/6) */
      LALCreateVector( status->statusPtr, &(outputPtr->xfacVec), 
          params->numPoints/2 + 1 );
      BEGINFAIL( status )
      {
        LALFree( outputPtr );
        *output = NULL;
      }
      ENDFAIL( status );

      xfac = outputPtr->xfacVec->data;

      xfac[0] = 0;
      for (k = 1; k < outputPtr->xfacVec->length; ++k) 
        xfac[k] = pow( (REAL4) k, exponent );
      break;

    default:
      /* unknown approximant type */
      LALFree( outputPtr );
      *output = NULL;
      ABORT( status, FINDCHIRPH_EUAPX, FINDCHIRPH_MSGEUAPX );
      break;
  }
      
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN (status);
}



/* <lalVerbatim file="FindChirpTemplateCP"> */
void
LALFindChirpTemplateFinalize (
    LALStatus                  *status,
    FindChirpTmpltParams      **output
    )
/* </lalVerbatim> */
{
  FindChirpTmpltParams         *outputPtr;

  INITSTATUS( status, "LALFindChirpTemplateFinalize", FINDCHIRPTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure handle is non-null and points to a non-null pointer */
  ASSERT( output, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( *output, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  

  /*
   *
   * destroy tmplt generation parameters structure
   *
   */


  /* local pointer to output */
  outputPtr = *output;

  /* destroy the fft plan if it exists */
  if ( outputPtr->fwdPlan )
  {
    LALDestroyRealFFTPlan( status->statusPtr, &(outputPtr->fwdPlan) );
    CHECKSTATUSPTR( status );
  }

  /* destroy the vector used to store part of the template */
  LALDestroyVector( status->statusPtr, &(outputPtr->xfacVec) );
  CHECKSTATUSPTR( status );

  /* free the structure */
  LALFree( outputPtr );
  *output = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
