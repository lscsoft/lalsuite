/************************** <lalVerbatim file="GenerateBurstCV">
Author: Brady, P. B.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Module \texttt{GenerateBurst.c}}
\label{ss:GenerateBurst.c}

Computes one of the standard burst waveforms with specified $h_{rss}$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GenerateBurstCP}
\idx{LALGenerateBurst()}

\subsubsection*{Description}

This function computes one of the following burst waveforms:
\begin{description}
\item[Sine-Gaussian]:  a linearly polarized sine-Gaussian with the specified
frequency and decay constant.
\end{description}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                   LALFree()
LALSCreateVectorSequence()    LALSDestroyVectorSequence()
LALSCreateVector()            LALSDestroyVector()
LALDCreateVector()            LALDDestroyVector()
LALSnprintf()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GenerateBurstCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateBurst.h>
#include <lal/LIGOMetadataTables.h>

NRCSID( GENERATEBURSTC, "$Id$" );

/* <lalVerbatim file="GenerateBurstCP"> */
void
LALGenerateBurst( 
    LALStatus          *stat, 
    CoherentGW         *output, 
    BurstParamStruc    *params 
    )
{ /* </lalVerbatim> */
  UINT4 n, i;          /* number of and index over samples */
  REAL8 t, dt;         /* time, interval */
  REAL8 t0, tau, gtime;  /* central time, decay time, gaussian time */
  REAL8 f0, phi0;      /* initial phase and frequency */
  REAL8 twopif0;       /* 2*pi*f0 */
  REAL8 f;             /* current value of frequency */
  REAL4 hrss;          /* root sum square strain for burst */
  REAL4 df = 0.0;      /* maximum difference between f */
  REAL8 phi;           /* current value of phase */
  REAL4 *fData;        /* pointer to frequency data */
  REAL8 *phiData;      /* pointer to phase data */
  REAL4 *aData;        /* pointer to frequency data */

  INITSTATUS( stat, "LALGenerateBurst", GENERATEBURSTC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter and output structures exist. */
  ASSERT( params, stat, GENERATEBURSTH_ENUL,
	  GENERATEBURSTH_MSGENUL );
  ASSERT( output, stat, GENERATEBURSTH_ENUL,
	  GENERATEBURSTH_MSGENUL );

  /* Make sure output fields don't exist. */
  ASSERT( !( output->a ), stat, GENERATEBURSTH_EOUT,
	  GENERATEBURSTH_MSGEOUT );
  ASSERT( !( output->f ), stat, GENERATEBURSTH_EOUT,
	  GENERATEBURSTH_MSGEOUT );
  ASSERT( !( output->phi ), stat, GENERATEBURSTH_EOUT,
	  GENERATEBURSTH_MSGEOUT );
  ASSERT( !( output->shift ), stat, GENERATEBURSTH_EOUT,
	  GENERATEBURSTH_MSGEOUT );

  /* Set up some other constants, to avoid repeated dereferencing. */
  n = params->length;
  dt = params->deltaT;
  f0 = params->f0;
  twopif0 = f0*LAL_TWOPI;

  /* Allocate output structures. */
  if ( ( output->a = (REAL4TimeVectorSeries *)
	 LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  memset( output->a, 0, sizeof(REAL4TimeVectorSeries) );
  if ( ( output->f = (REAL4TimeSeries *)
	 LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  memset( output->f, 0, sizeof(REAL4TimeSeries) );
  if ( ( output->phi = (REAL8TimeSeries *)
	 LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    LALFree( output->f ); output->f = NULL;
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  memset( output->phi, 0, sizeof(REAL8TimeSeries) );

  /* Set output structure metadata fields. */
  output->position = params->position;
  output->psi = params->psi;
  output->a->epoch = output->f->epoch = output->phi->epoch
    = params->epoch;
  output->a->deltaT = params->deltaT;
  output->f->deltaT = output->phi->deltaT = params->deltaT;
  output->a->sampleUnits = lalStrainUnit;
  output->f->sampleUnits = lalHertzUnit;
  output->phi->sampleUnits = lalDimensionlessUnit;
  LALSnprintf( output->a->name, LALNameLength, "Burst amplitudes" );
  LALSnprintf( output->a->name, LALNameLength, "Burst frequency" );
  LALSnprintf( output->a->name, LALNameLength, "Burst phase" );

  /* Allocate phase and frequency arrays. */
  LALSCreateVector( stat->statusPtr, &( output->f->data ), n );
  BEGINFAIL( stat ) {
    LALFree( output->a );   output->a = NULL;
    LALFree( output->f );   output->f = NULL;
    LALFree( output->phi ); output->phi = NULL;
  } ENDFAIL( stat );
  LALDCreateVector( stat->statusPtr, &( output->phi->data ), n );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
	 stat );
    LALFree( output->a );   output->a = NULL;
    LALFree( output->f );   output->f = NULL;
    LALFree( output->phi ); output->phi = NULL;
  } ENDFAIL( stat );

  /* Allocate amplitude array. */
  {
    CreateVectorSequenceIn in; /* input to create output->a */
    in.length = 2;
    in.vectorLength = n;
    LALSCreateVectorSequence( stat->statusPtr, &(output->a->data), &in );
    BEGINFAIL( stat ) {
      TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
	   stat );
      TRY( LALDDestroyVector( stat->statusPtr, &( output->phi->data ) ),
	   stat );
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
  }

  /* Fill frequency and phase arrays. */
  fData = output->f->data->data;
  phiData = output->phi->data->data;
  aData = output->a->data->data; 

  /* this depends on the waveform type */
  switch ( params->burstType )
  {
    /* sine-Gaussian burst */
    case sineGaussian:
      for ( i = 0; i < n; i++ ) {
        t = i*dt;
        gtime = (t-t0)/tau;
        *(fData++) = f0;
        *(phiData++) = twopif0 * t;
        *(aData++) = hrss * exp( - gtime * gtime );
        *(aData++) = 0.0;
      }
    /* Gaussian burst */
    case Gaussian:
      for ( i = 0; i < n; i++ ) {
        t = i*dt;
        gtime = (t-t0)/tau;
        *(fData++) = 0.0;
        *(phiData++) = 0.0;
        *(aData++) = hrss * exp( - gtime * gtime );
        *(aData++) = 0.0;
      }
  }

  /* Set output field and return. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


void
LALBurstInjectSignals( 
    LALStatus               *stat, 
    REAL4TimeSeries         *series, 
    SimBurstTable           *injections,
    COMPLEX8FrequencySeries *resp
    )
{
  UINT4              k;
  INT8               chanStartTime;
  DetectorResponse   detector;
  COMPLEX8Vector     *unity = NULL;
  CoherentGW         waveform;
  BurstParamStruc    burstParam;
  REAL4TimeSeries    signal;

  INITSTATUS( stat, "LALBurstInjectSignals", GENERATEBURSTC );
  ATTATCHSTATUSPTR( stat );


  /*  set up parameters */
  LALGPStoINT8( stat->statusPtr, &chanStartTime, &(series->epoch) );
  CHECKSTATUSPTR( stat );

  /* 
   *compute the transfer function 
   */

  /* allocate memory and copy the parameters describing the freq series */
  memset( &detector, 0, sizeof( DetectorResponse ) );
  detector.transfer = (COMPLEX8FrequencySeries *)
    LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
  if ( ! detector.transfer ) 
  {
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  memcpy( &(detector.transfer->epoch), &(resp->epoch),
      sizeof(LIGOTimeGPS) );
  detector.transfer->f0 = resp->f0;
  detector.transfer->deltaF = resp->deltaF;

  detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );
  /* set the detector site */
  switch ( series->name[0] )
  {
    case 'H':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexLHODIFF];
      LALWarning( stat, "computing waveform for Hanford." );
      break;
    case 'L':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
      LALWarning( stat, "computing waveform for Livingston." );
      break;
    default:
      LALFree( detector.site );
      detector.site = NULL;
      LALWarning( stat, "Unknown detector site, computing plus mode "
          "waveform with no time delay" );
      break;
  }

  /* set up units for the transfer function */
  {
    RAT4 negOne = { -1, 0 };
    LALUnit unit;
    LALUnitPair pair;
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &lalStrainUnit;
    LALUnitRaise( stat->statusPtr, &unit, pair.unitTwo, &negOne );
    CHECKSTATUSPTR( stat );
    pair.unitTwo = &unit;
    LALUnitMultiply( stat->statusPtr, &(detector.transfer->sampleUnits),
        &pair );
    CHECKSTATUSPTR( stat );
  }

  /* invert the response function to get the transfer function */
  LALCCreateVector( stat->statusPtr, &( detector.transfer->data ),
      resp->data->length );
  CHECKSTATUSPTR( stat );

  LALCCreateVector( stat->statusPtr, &unity, resp->data->length );
  CHECKSTATUSPTR( stat );
  for ( k = 0; k < resp->data->length; ++k ) 
  {
    unity->data[k].re = 1.0;
    unity->data[k].im = 0.0;
  }

  LALCCVectorDivide( stat->statusPtr, detector.transfer->data, unity,
      resp->data );
  CHECKSTATUSPTR( stat );

  LALCDestroyVector( stat->statusPtr, &unity );
  CHECKSTATUSPTR( stat );

  /*
   * inject into time series 
   */

  burstParam.position.longitude  = injections->longitude;;
  burstParam.position.latitude   = injections->latitude;;
  burstParam.position.system     = COORDINATESYSTEM_EQUATORIAL;
  burstParam.psi = LAL_PI/2.0;
  burstParam.epoch.gpsSeconds =  injections->geocent_start_time.gpsSeconds;
  burstParam.epoch.gpsNanoSeconds =  injections->geocent_start_time.gpsNanoSeconds;
  burstParam.deltaT = series->deltaT;
  burstParam.length = series->data->length;
  burstParam.hrss = injections->hrss;
  burstParam.burstType = sineGaussian;
  burstParam.f0 = (REAL8)injections->freq;
  burstParam.tau  = (REAL8)injections->tau;
  
 /*generate waveform*/

  memset( &waveform, 0, sizeof(CoherentGW) );

  LALGenerateBurst( &stat, &waveform, &burstParam );

  /*set the parameters for the to be injected signal */
  signal.deltaT = series->deltaT;
   if ( ( signal.f0 = series->f0 ) != 0 )
    {
      ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
    }
    signal.sampleUnits = lalADCCountUnit;

   /* simulate the detectors response  */
    LALSCreateVector( stat->statusPtr, &(signal.data), 
        (UINT4)burstParam.length );
    CHECKSTATUSPTR( stat );

    LALSimulateCoherentGW( stat->statusPtr, 
        &signal, &waveform, &detector );
    CHECKSTATUSPTR( stat );

    /* inject the signal into the data channel */
    LALSSInjectTimeSeries( stat->statusPtr, series, &signal );
    CHECKSTATUSPTR( stat );

   /* destroy the signal */
    LALSDestroyVector( stat->statusPtr, &(signal.data) );
    CHECKSTATUSPTR( stat );

  LALCDestroyVector( stat->statusPtr, &( detector.transfer->data ) );
  CHECKSTATUSPTR( stat );

  if ( detector.site ) LALFree( detector.site );
  LALFree( detector.transfer );

  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
