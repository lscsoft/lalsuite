/**** <lalVerbatim file="RingSearchInitCV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/Ring.h>
#include <lal/RingSearch.h>

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
 * \item[\texttt{-fhighpass} \textit{fhighpass}] Highpass filter the data at
 *   frequency \textit{flow} Hz (negative to disable) [default is zero].
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


NRCSID( RINGSEARCHINITC, "$Id$" );

/* <lalVerbatim file="RingSearchInitCP"> */
void LALRingSearchInit(
    LALStatus         *status,
    RingSearchParams **searchParams,
    CHAR             **argv,
    INT4               argc
    )
{ /* </lalVerbatim> */
  RingSearchParams      *params;
  RingTemplateBankInput  bankin;

  INITSTATUS( status, "LALRingSearchInit", RINGSEARCHINITC );
  ATTATCHSTATUSPTR( status );

  ASSERT( argc < 1 ||  argv, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( argc < 1 || *argv, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( searchParams, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( ! *searchParams, status, RINGSEARCHH_ENNUL, RINGSEARCHH_MSGENNUL );

  params = *searchParams = LALCalloc( 1, sizeof( *params ) );
  if ( ! params )
  {
    ABORT( status, RINGSEARCHH_EALOC, RINGSEARCHH_MSGEALOC );
  }

  params->numSlaves      = -1;
  params->myProcNumber   = -1;
  params->dynRangeFac    =  1;
  params->maximizeEvents = -1; /* natural amount: duration of filter */
  params->avgSpecMeth    = useMedian; /* default method */

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
    else if ( strstr( *argv, "-avgspec" ) )
    {
      char *meth = *++argv;
      --argc;
      if ( strstr( meth, "unity" ) )
        params->avgSpecMeth = useUnity;
      else if ( strstr( meth, "mean" ) )
        params->avgSpecMeth = useMean;
      else
        params->avgSpecMeth = useMedian; /* default */
    }
    else if ( strstr( *argv, "-flow" ) )
    {
      params->lowFrequency = atof( *++argv );
      --argc;
    }
    else if ( strstr( *argv, "-fhighpass" ) )
    {
      params->highpassFrequency = atof( *++argv );
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
    else if ( strstr( *argv, "-qmin" ) )
    {
      params->minQuality = atof( *++argv );
      --argc;
    }
    else if ( strstr( *argv, "-qmax" ) )
    {
      params->maxQuality = atof( *++argv );
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
      ABORT( status, RINGSEARCHH_EIOPT, RINGSEARCHH_MSGEIOPT );
    }
  }

  if ( ! params->segmentSize )
  {
    ABORT( status, RINGSEARCHH_ESIZE, RINGSEARCHH_MSGESIZE );
  }

  if ( ! ( params->lowFrequency > 0 ) )
  {
    ABORT( status, RINGSEARCHH_EFLOW, RINGSEARCHH_MSGEFLOW );
  }

  if ( ! ( params->minFrequency > 0
        && params->maxFrequency > params->minFrequency ) )
  {
    ABORT( status, RINGSEARCHH_EFREQ, RINGSEARCHH_MSGEFREQ );
  }

  if ( ! ( params->minQuality > 0
        && params->maxQuality > params->minQuality ) )
  {
    ABORT( status, RINGSEARCHH_EQUAL, RINGSEARCHH_MSGEQUAL );
  }


  /* memory for invSpectrum */

  params->invSpectrum = LALCalloc( 1, sizeof( *params->invSpectrum ) );
  if ( ! params->invSpectrum )
  {
    ABORT( status, RINGSEARCHH_EALOC, RINGSEARCHH_MSGEALOC );
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

  bankin.minQuality   = params->minQuality;
  bankin.maxQuality   = params->maxQuality;
  bankin.minFrequency = params->minFrequency;
  bankin.maxFrequency = params->maxFrequency;
  bankin.maxMismatch  = params->maxMismatch;
  LALCreateRingTemplateBank( status->statusPtr, &params->templateBank,
      &bankin );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="RingSearchInitCP"> */
void
LALRingSearchFini(
    LALStatus         *status,
    RingSearchParams **searchParams
    )
{ /* </lalVerbatim> */
  RingSearchParams *params;

  INITSTATUS( status, "LALRingSearchFini", RINGSEARCHINITC );
  ATTATCHSTATUSPTR( status );

  ASSERT(  searchParams, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( *searchParams, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );

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
    LALDestroyRingTemplateBank( status->statusPtr, &params->templateBank );
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
