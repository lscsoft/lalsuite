/**** <lalVerbatim file="RingSearchInitCV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/String.h>
#include <lal/StringSearch.h>

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{RingSearchInit.c}}
 *
 * Routines to initialize and finalize a ring search.
 *
 * \subsubsection*{Prototypes}
 * \input{RingSearchInitCP}
 * \idx{LALRingSearchInit()}
 * \idx{LALRingSearchFini()}
 * 
 * \subsubsection*{Description}
 *
 * The routine \verb+LALRingSearchInit()+ creates and initialized a
 * \verb+RingSearchParams+ structure that will be used for a given ring
 * search.  Some of the required memory in this structure is allocated and
 * the ring template bank is generated.  As input, it requires parameters
 * specified by the argument vector \verb+argv+ (with \verb+argc+ arguments).
 *
 * The allowed arguments are:
 * \begin{description}
 * \item[\texttt{-segsz} \textit{segsz}] Set the number of points in each
 *   segment of data to \textit{segsz} [no default --- manditory].
 * \item[\texttt{-scale} \textit{scale}] Set the dynamical range scale to
 *   \textit{scale} [default is unity].
 * \item[\texttt{-speclen} \textit{speclen}] Set the number of points in the
 *   inverse spectrum truncation to \textit{speclen}
 *   [default is zero --- no inverse spectrum truncation].
 * \item[\texttt{-flow} \textit{flow}] Set the low frequency cutoff to
 *   \textit{flow} Hz [default is zero].
 * \item[\texttt{-fmin} \textit{fmin}] Set the minimum ring frequency of the
 *   bank to \textit{fmin} Hz [no default --- manditory].
 * \item[\texttt{-fmax} \textit{fmax}] Set the maximum ring frequency of the
 *   bank to \textit{fmax} Hz [no default --- manditory].
 * \item[\texttt{-qmin} \textit{qmin}] Set the minimum ring quality factor
 *   of the bank to \textit{qmin} [no default --- manditory].
 * \item[\texttt{-qmax} \textit{qmax}] Set the maximum ring quality factor
 *   of the bank to \textit{qmax} [no default --- manditory].
 * \item[\texttt{-maxmm} \textit{maxmm}] Set the maximum mismatch of templates
 *   in the bank to \textit{maxmm} [no default --- manditory].
 * \item[\texttt{-thresh} \textit{thresh}] Set the snr threshold for events
 *   to \textit{thresh} [no default --- manditory].
 * \end{description}
 *
 * The zeroth argument \verb+argv[0]+ is ignored and \verb+argv[argc]+ should
 * be \verb+NULL+.  For example:
 * \begin{verbatim}
 * const char *argv[] = { "filterparams", "-segsz", "65536", "-speclen", "4096",
 *    "-flow", "40", "-fmin", "150", "-fmax", "200", "-qmin", "2",
 *    "-qmax", "10", "-maxmm", "0.1", "-thresh", "6", "-scale", "2000", NULL };
 * int argc = sizeof( argv ) / sizeof( *argv ) - 1;
 * \end{verbatim}
 *
 * The function \verb+LALRingSearchFini()+ cleans up all memory allocated in
 * the parameter structure.
 *
 * \vfill{\footnotesize\input{RingCV}}
 *
 **** </lalLaTeX> */


NRCSID( STRINGSEARCHINITC, "$Id$" );

/* <lalVerbatim file="RingSearchInitCP"> */
void LALStringSearchInit(
    LALStatus         *status,
    StringSearchParams **searchParams,
    const CHAR       **argv,
    INT4               argc
    )
{ /* </lalVerbatim> */
  StringSearchParams      *params;
  /* StringTemplateBankInput  bankin; */

  INITSTATUS( status, "StringSearchInit", STRINGSEARCHINITC );
  ATTATCHSTATUSPTR( status );

  ASSERT( argc < 1 ||  argv, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );
  ASSERT( argc < 1 || *argv, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );
  ASSERT( searchParams, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );
  ASSERT( ! *searchParams, status, STRINGSEARCHH_ENNUL, STRINGSEARCHH_MSGENNUL );

  params = *searchParams = LALCalloc( 1, sizeof( *params ) );
  if ( ! params )
  {
    ABORT( status, STRINGSEARCHH_EALOC, STRINGSEARCHH_MSGEALOC );
  }

  params->numSlaves      = -1;
  params->myProcNumber   = -1;
  params->dynRangeFac    =  1;
  params->maximizeEvents = -1; /* natural amount: duration of filter */

  while ( --argc > 0 )
  {
    ++argv;
    if ( strstr( *argv, "-segsz" ) )
    {
      params->segmentSize = atoi( *++argv );
      --argc;
    }
    else if ( strstr( *argv, "-scale" ) )
    {
      params->dynRangeFac = atof( *++argv );
      --argc;
    }
    else if ( strstr( *argv, "-speclen" ) )
    {
      params->invSpecTrunc = atoi( *++argv );
      --argc;
    }
    else if ( strstr( *argv, "-flow" ) )
    {
      params->lowFrequency = atof( *++argv );
      --argc;
    }
    else if ( strstr( *argv, "-expt" ) )
    {
      params->freqPower = atof( *++argv);
      --argc;
    }
    else if ( strstr( *argv, "-fmin" ) )
    {
      params->minFrequency = atof( *++argv );
      --argc;
    }
    else if ( strstr( *argv, "-fmax" ) )
    {
      params->maxFrequency = atof( *++argv );
      --argc;
    }
    else if ( strstr( *argv, "-maxmm" ) )
    {
      params->maxMismatch = atof( *++argv );
      --argc;
    }
    else if ( strstr( *argv, "-thresh" ) )
    {
      params->threshold = atof( *++argv );
      --argc;
    }
    else
    {
      ABORT( status, STRINGSEARCHH_EIOPT, STRINGSEARCHH_MSGEIOPT );
    }
  }

  if ( ! params->segmentSize )
  {
    ABORT( status, STRINGSEARCHH_ESIZE, STRINGSEARCHH_MSGESIZE );
  }

  if ( ! ( params->lowFrequency > 0 ) )
  {
    ABORT( status, STRINGSEARCHH_EFLOW, STRINGSEARCHH_MSGEFLOW );
  }

  if ( ! ( params->minFrequency > 0
        && params->maxFrequency > params->minFrequency ) )
  {
    ABORT( status, STRINGSEARCHH_EFREQ, STRINGSEARCHH_MSGEFREQ );
  }


  /* memory for invSpectrum */

  params->invSpectrum = LALCalloc( 1, sizeof( *params->invSpectrum ) );
  if ( ! params->invSpectrum )
  {
    ABORT( status, STRINGSEARCHH_EALOC, STRINGSEARCHH_MSGEALOC );
  }

  LALSCreateVector( status->statusPtr, &params->invSpectrum->data,
      params->segmentSize / 2 + 1 );
  CHECKSTATUSPTR( status );

  LALCreateForwardRealFFTPlan( status->statusPtr, &params->forwardPlan,
      params->segmentSize, 0 );
  CHECKSTATUSPTR( status );
  
  LALCreateReverseRealFFTPlan( status->statusPtr, &params->reversePlan,
      params->segmentSize, 0 );
  CHECKSTATUSPTR( status );

/*  bankin.minQuality   = params->minQuality;
  bankin.maxQuality   = params->maxQuality;
  bankin.minFrequency = params->minFrequency;
  bankin.maxFrequency = params->maxFrequency;
  bankin.maxMismatch  = params->maxMismatch;
  LALCreateStringTemplateBank( status->statusPtr, &params->templateBank,
      &bankin );
  CHECKSTATUSPTR( status );
*/
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="RingSearchInitCP"> */
void
LALStringSearchFini(
    LALStatus         *status,
    StringSearchParams **searchParams
    )
{ /* </lalVerbatim> */
  StringSearchParams *params;

  INITSTATUS( status, "StringSearchFini", STRINGSEARCHINITC );
  ATTATCHSTATUSPTR( status );

  ASSERT(  searchParams, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );
  ASSERT( *searchParams, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );

  params = *searchParams;

  while ( params->numResults )
  {
    --params->numResults;
    if ( params->result[params->numResults].data )
    {
      LALSDestroyVector( status->statusPtr,
          &params->result[params->numResults].data );
      CHECKSTATUSPTR( status );
    }
  }
  if ( params->result )
  {
    LALFree( params->result );
    params->result = NULL;
  }

  while ( params->numSegments )
  {
    --params->numSegments;
    if ( params->dataSegment[params->numSegments].data )
    {
      LALCDestroyVector( status->statusPtr,
          &params->dataSegment[params->numSegments].data );
      CHECKSTATUSPTR( status );
    }
  }
  if ( params->dataSegment )
  {
    LALFree( params->dataSegment );
    params->dataSegment = NULL;
  }

  if ( params->templateBank )
  {
    LALDestroyStringTemplateBank( status->statusPtr, &params->templateBank );
    CHECKSTATUSPTR( status );
  }

  if ( params->reversePlan )
  {
    LALDestroyRealFFTPlan( status->statusPtr, &params->reversePlan );
    CHECKSTATUSPTR( status );
  }

  if ( params->forwardPlan )
  {
    LALDestroyRealFFTPlan( status->statusPtr, &params->forwardPlan );
    CHECKSTATUSPTR( status );
  }
  
  if ( params->invSpectrum )
  {
    if ( params->invSpectrum->data )
    {
      LALSDestroyVector( status->statusPtr, &params->invSpectrum->data );
      CHECKSTATUSPTR( status );
    }
    LALFree( params->invSpectrum );
    params->invSpectrum = NULL;
  }

  LALFree( params );
  params = NULL;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
