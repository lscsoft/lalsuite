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
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/VectorOps.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateBurst.h>
#include <lal/LIGOMetadataTables.h>

NRCSID( GENERATEBURSTC, "$Id$" );

/* <lalVerbatim file="GenerateBurstCP"> */
void
LALGenerateBurst( 
    LALStatus          *stat, 
    CoherentGW         *output,
    SimBurstTable      *simBurst,
    BurstParamStruc    *params 
    )
{ /* </lalVerbatim> */
  UINT4 n, i;          /* number of and index over samples */
  REAL8 t, dt, duration;         /* time, interval */
  REAL8 t0, tau, gtime;  /* central time, decay time, gaussian time */
  REAL8 f0, phi0;      /* initial phase and frequency */
  REAL8 twopif0;       /* 2*pi*f0 */
  REAL8 f;             /* current value of frequency */
  REAL4 hpeak;         /* peak strain for burst */
  REAL4 df = 0.0;      /* maximum difference between f */
  REAL8 phi;           /* current value of phase */
  REAL4 *fData;        /* pointer to frequency data */
  REAL8 *phiData;      /* pointer to phase data */
  REAL4 *aData;        /* pointer to frequency data */
  LIGOTimeGPS startTime;  /* start time of injection */
  LALTimeInterval dummyInterval;

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
  duration = (REAL8)(simBurst->dtplus + simBurst->dtminus);
  dt = params->deltaT;
  if ( ( n = (INT4) (2.0 * duration / dt) ) == 0 )
  {
    ABORT(stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }

  /* start time of data is peak time duration */
  TRY( LALFloatToInterval( stat->statusPtr, &dummyInterval, &duration ), stat );
  TRY( LALDecrementGPS( stat->statusPtr, &startTime, 
        &(simBurst->geocent_peak_time), &dummyInterval), stat);

  /* Generic burst parameters */
  hpeak = simBurst->hpeak;
  tau = (REAL8)simBurst->tau;
  f0 = (REAL8)simBurst->freq;
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
  output->position.longitude = simBurst->longitude;
  output->position.latitude = simBurst->latitude;
  output->position.system = params->system;
  output->psi = simBurst->polarization;
  output->a->epoch = output->f->epoch = output->phi->epoch = startTime;
  output->a->deltaT = params->deltaT;
  output->f->deltaT = output->phi->deltaT = params->deltaT;
  output->a->sampleUnits = lalStrainUnit;
  output->f->sampleUnits = lalHertzUnit;
  output->phi->sampleUnits = lalDimensionlessUnit;
  LALSnprintf( output->a->name, LALNameLength, "Burst amplitudes" );
  LALSnprintf( output->f->name, LALNameLength, "Burst frequency" );
  LALSnprintf( output->phi->name, LALNameLength, "Burst phase" );

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
    in.length = n;
    in.vectorLength = 2;
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
  if( !( strcmp( simBurst->waveform, "SineGaussian" ) ) )
  {
    /* find the peak time as a REAL8 relative to start of segment */
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );

    /* construct the signal */
    for ( i = 0; i < n; i++ ) {
      t = i*dt;
      gtime = (t-t0)/tau;
      *(fData++) = f0;
      *(phiData++) = twopif0 * (t-t0);
      *(aData++) = hpeak * exp( - gtime * gtime );
      *(aData++) = 0.0;
    }
  }
  else if ( !( strcmp( simBurst->waveform, "Gaussian" ) ) )
  {
    /* find the peak time as a REAL8 relative to start of segment */
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );

    /* construct the signal */
    for ( i = 0; i < n; i++ ) {
      t = i*dt;
      gtime = (t-t0)/tau;
      *(fData++) = 0.0;
      *(phiData++) = 0.0;
      *(aData++) = hpeak * exp( - gtime * gtime );
      *(aData++) = 0.0;
    }
  }
  else if ( !( strcmp( simBurst->waveform, "Ringdown" ) ) )
  {
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );
    for ( i = 0; i < n; i++ )
    {
      t = i * dt;
      gtime = ( t - t0 ) / tau;
      *fData++   = f0;
      *phiData++ = twopif0 * ( t - t0 );
      if ( gtime > 0 )
        *aData++ = hpeak * exp( - gtime );
      else
        *aData++ = 0;
      *aData++   = 0;
    }
  }
  else if ( !( strcmp( simBurst->waveform, "Ringup" ) ) )
  {
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );
    for ( i = 0; i < n; i++ )
    {
      t = i * dt;
      gtime = ( t - t0 ) / tau;
      *fData++   = f0;
      *phiData++ = twopif0 * ( t - t0 );
      if ( gtime < 0 )
        *aData++ = hpeak * exp( gtime );
      else
        *aData++ = 0;
      *aData++   = 0;
    }
  }

  else if ( !( strcmp( simBurst->waveform, "StringCusp" ) ) )
  {
    TRY( LALDeltaGPS( stat->statusPtr, &dummyInterval, 
          &(simBurst->geocent_peak_time), &startTime ), stat );
    TRY( LALIntervalToFloat( stat->statusPtr, &t0, 
          &dummyInterval ), stat );
    for ( i = 0; i < n; i++ )
    {
      t = i * dt;
      gtime = fabs( t - t0 );
      *fData++   = 0.0;
      *phiData++ = 0.0;
      *aData++ = hpeak * (pow ( pow(duration,2)+pow(1/f0,2),1./6.) -
			  pow ( pow(gtime,2)   +pow(1/f0,2),1./6.) ); 
      *aData++   = 0;
    }
  }
  else
  {
    ABORT( stat, GENERATEBURSTH_ETYP, GENERATEBURSTH_MSGETYP );
  }

  /* Set output field and return. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}





/* <lalVerbatim file="GenerateBurstCP"> */
void
LALBurstInjectSignals( 
    LALStatus               *stat, 
    REAL4TimeSeries         *series, 
    SimBurstTable           *injections,
    COMPLEX8FrequencySeries *resp
    )
/* </lalVerbatim> */
{
  UINT4              k;
  INT4               injStartTime;
  INT4               injStopTime;
  DetectorResponse   detector;
  COMPLEX8Vector    *unity = NULL;
  CoherentGW         waveform;
  BurstParamStruc    burstParam;
  REAL4TimeSeries    signal;
  SimBurstTable     *simBurst=NULL;
  LALDetector       *tmpDetector=NULL,*nullDetector=NULL;

  INITSTATUS( stat, "LALBurstInjectSignals", GENERATEBURSTC );
  ATTATCHSTATUSPTR( stat );

  /* set up start and end of injection zone TODO: fix this hardwired 10 */
  injStartTime = series->epoch.gpsSeconds - 10;
  injStopTime = series->epoch.gpsSeconds + 10 + (INT4)(series->data->length
      * series->deltaT);

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

  tmpDetector = detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );
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
      tmpDetector = NULL;
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

  /* Set up a time series to hold signal in ADC counts */
  signal.deltaT = series->deltaT;
  if ( ( signal.f0 = series->f0 ) != 0 )
  {
    ABORT( stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM );
  }
  signal.sampleUnits = lalADCCountUnit;

  signal.data=NULL;
  LALSCreateVector( stat->statusPtr, &(signal.data), 
      series->data->length );
  CHECKSTATUSPTR( stat );

  /* loop over list of waveforms and inject into data stream */
  simBurst = injections;
  while ( simBurst )
  {
    /* only do the work if the burst is in injection zone */
    if( (injStartTime - simBurst->geocent_peak_time.gpsSeconds) *
        (injStopTime - simBurst->geocent_peak_time.gpsSeconds) > 0 )
    {
      simBurst = simBurst->next;
      continue;
    }

    /* set the burt params */
    burstParam.deltaT = series->deltaT;
    if( !( strcmp( simBurst->coordinates, "HORIZON" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_HORIZON;
    }
    else if ( !( strcmp( simBurst->coordinates, "ZENITH" ) ) )
    {
      /* set coordinate system for completeness */
      burstParam.system = COORDINATESYSTEM_EQUATORIAL;
      detector.site = NULL;
    }
    else if ( !( strcmp( simBurst->coordinates, "GEOGRAPHIC" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_GEOGRAPHIC;
    }
    else if ( !( strcmp( simBurst->coordinates, "EQUATORIAL" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_EQUATORIAL;
    }
    else if ( !( strcmp( simBurst->coordinates, "ECLIPTIC" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_ECLIPTIC;
    }
    else if ( !( strcmp( simBurst->coordinates, "GALACTIC" ) ) )
    {
      burstParam.system = COORDINATESYSTEM_GALACTIC;
    }
    else
      burstParam.system = COORDINATESYSTEM_EQUATORIAL;

    /* generate the burst */
    memset( &waveform, 0, sizeof(CoherentGW) );
    LALGenerateBurst( stat->statusPtr, &waveform, simBurst, &burstParam );
    CHECKSTATUSPTR( stat );

    /* must set the epoch of signal since it's used by coherent GW */
    signal.epoch = waveform.a->epoch;
    memset( signal.data->data, 0, signal.data->length * sizeof(REAL4) );
    
    /* convert this into an ADC signal */
    LALSimulateCoherentGW( stat->statusPtr, 
        &signal, &waveform, &detector );
    CHECKSTATUSPTR( stat );

    /* inject the signal into the data channel */
    LALSSInjectTimeSeries( stat->statusPtr, series, &signal );
    CHECKSTATUSPTR( stat );

    /* free memory in coherent GW structure.  TODO:  fix this */
    LALSDestroyVectorSequence( stat->statusPtr, &( waveform.a->data ) );
    CHECKSTATUSPTR( stat );
    LALSDestroyVector( stat->statusPtr, &( waveform.f->data ) );
    CHECKSTATUSPTR( stat );
    LALDDestroyVector( stat->statusPtr, &( waveform.phi->data ) );
    CHECKSTATUSPTR( stat );
    LALFree( waveform.a );   waveform.a = NULL;
    LALFree( waveform.f );   waveform.f = NULL;
    LALFree( waveform.phi );  waveform.phi = NULL;

    /* reset the detector site information in case it changed */
    detector.site = tmpDetector;

    /* move on to next one */
    simBurst = simBurst->next;
  }

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
