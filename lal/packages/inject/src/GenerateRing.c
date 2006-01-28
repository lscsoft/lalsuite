/************************** 
Author: Goggin, L. M., based on GenerateBurst.c by Brady, P. B.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Module \texttt{GenerateRing.c}}
\label{ss:GenerateRing.c}

Computes the ringdown waveform with specified $h_{rss}$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GenerateRingCP}
\idx{LALGenerateRing()}

\subsubsection*{Description}

This function the following burst waveforms:
\begin{description}
\item[Sine-Gaussian]:  exponentially decaying sinusoid with specified frequency and decay constant.
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

\vfill{\footnotesize\input{GenerateRingCV}}

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
#include <lal/GenerateRing.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>


NRCSID( GENERATERINGC, "$Id$" );

/* <lalVerbatim file="GenerateRingCP"> */
void
LALGenerateRing( 
    LALStatus          *stat, 
    CoherentGW         *output,
    REAL4TimeSeries    *series,
    SimRingdownTable   *simRingdown,
    RingParamStruc     *params
    )
    
{ /* </lalVerbatim> */
  UINT4 n, i;      /* number of and index over samples */
  REAL8 t, dt;         /* time, interval */
  REAL8 t0, gtime ;    /* central time, decay time */
  REAL8 f0, quality;   /* frequency and quality factor */
  REAL8 twopif0;       /* 2*pi*f0 */
  REAL4 h0;            /* peak strain for ringdown */
  REAL4 *fData;        /* pointer to frequency data */
  REAL8 *phiData;      /* pointer to phase data */
  REAL4 *aData;        /* pointer to frequency data */
  LIGOTimeGPS startTime;  /* start time of injection */
  REAL4TimeSeries signal; /* start time of block that injection is injected into */
  UINT4 N_point_block; /* number of data points in a block */
  UINT4 start_i, bb;   /* index at which to start injection, debugging const */
  INT8 geoc_tns;       /* geocentric_start_time of the injection in ns */
  INT8 block_tns;      /* start time of block in ns */
  INT8 deltaTns;        /* deltaT in ns */
  LALTimeInterval dummyInterval;
/* REAL8  deltaTns; */
  
  INITSTATUS( stat, "LALGenerateRing", GENERATERINGC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter and output structures exist. */
  ASSERT( params, stat, GENERATERINGH_ENUL,
	  GENERATERINGH_MSGENUL );
  ASSERT( output, stat, GENERATERINGH_ENUL,
	  GENERATERINGH_MSGENUL );

  /* Make sure output fields don't exist. */
  ASSERT( !( output->a ), stat, GENERATERINGH_EOUT,
	  GENERATERINGH_MSGEOUT );
  ASSERT( !( output->f ), stat, GENERATERINGH_EOUT,
	  GENERATERINGH_MSGEOUT );
  ASSERT( !( output->phi ), stat, GENERATERINGH_EOUT,
	  GENERATERINGH_MSGEOUT );
  ASSERT( !( output->shift ), stat, GENERATERINGH_EOUT,
	  GENERATERINGH_MSGEOUT );

  /* Set up some other constants, to avoid repeated dereferencing. */
  dt = params->deltaT; 
  startTime = simRingdown->geocent_start_time;
  N_point_block = series->data->length;
  
  /* Generic ring parameters */
  h0 = simRingdown->h0;
  quality = (REAL8)simRingdown->quality;
  f0 = (REAL8)simRingdown->frequency;
  twopif0 = f0*LAL_TWOPI;
  
  /* fprintf( stderr, "n = %d\n", n );*/

  /* Allocate output structures. */
  if ( ( output->a = (REAL4TimeVectorSeries *)
	 LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
    ABORT( stat, GENERATERINGH_EMEM, GENERATERINGH_MSGEMEM );
  }
  memset( output->a, 0, sizeof(REAL4TimeVectorSeries) );
  if ( ( output->f = (REAL4TimeSeries *)
	 LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    ABORT( stat, GENERATERINGH_EMEM, GENERATERINGH_MSGEMEM );
  }
  memset( output->f, 0, sizeof(REAL4TimeSeries) );
  if ( ( output->phi = (REAL8TimeSeries *)
	 LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    LALFree( output->f ); output->f = NULL;
    ABORT( stat, GENERATERINGH_EMEM, GENERATERINGH_MSGEMEM );
  }
  memset( output->phi, 0, sizeof(REAL8TimeSeries) );

  /* Set output structure metadata fields. */
  output->position.longitude = simRingdown->longitude; 
  output->position.latitude = simRingdown->latitude;
  output->position.system = params->system;
  output->psi = simRingdown->polarization;
   /* set epoch of output time series to that of the block */
  output->a->epoch = output->f->epoch = output->phi->epoch = series->epoch;
  output->a->deltaT = params->deltaT; 
  output->f->deltaT = output->phi->deltaT = params->deltaT; 
  output->a->sampleUnits = lalStrainUnit;
  output->f->sampleUnits = lalHertzUnit;
  output->phi->sampleUnits = lalDimensionlessUnit;
  LALSnprintf( output->a->name, LALNameLength, "Ring amplitudes" );
  LALSnprintf( output->f->name, LALNameLength, "Ring frequency" );
  LALSnprintf( output->phi->name, LALNameLength, "Ring phase" );

  /* Allocate phase and frequency arrays. */
  LALSCreateVector( stat->statusPtr, &( output->f->data ), N_point_block );
  BEGINFAIL( stat ) {
    LALFree( output->a );   output->a = NULL;
    LALFree( output->f );   output->f = NULL;
    LALFree( output->phi ); output->phi = NULL;
  } ENDFAIL( stat );
  
  LALDCreateVector( stat->statusPtr, &( output->phi->data ), N_point_block );
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
    in.length = N_point_block;
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
  
  /*  set arrays to zero */
  memset( output->f->data->data, 0, sizeof( REAL4 ) *  output->f->data->length );
  memset( output->phi->data->data, 0, sizeof( REAL8 ) * output->phi->data->length );
  memset( output->a->data->data, 0, sizeof( REAL4 ) * 
      output->a->data->length * output->a->data->vectorLength );
  
    
  geoc_tns = XLALGPStoINT8( &startTime );
  block_tns =  XLALGPStoINT8( &series->epoch );

  
  
  /* Find starting index, three options, 3rd way gives exactly the same result
   * for start_time in the sngl_ringdown table for found inject F as before any 
   * changes were made to the code.   */
  
  /* 1. This is what we coded on Wednesday, start_time in the sgnl_r table now
   * differs by 3ms*/
  deltaTns = (INT8) floor( 0.5 + 1e9 * series->deltaT ); 
  start_i =( geoc_tns - block_tns ) / deltaTns;   /* but the result of this does not 
                                                     necessarily give an integer, right? */
  
  /* 2. Perhaps its better to leave the rounding off until the end. (Change
   * declaration of deltaTns to REAL8)  */
/* deltaTns = 1e9 * series->deltaT;*/
/* start_i =(INT8) floor( 0.5 + ( ( geoc_tns - block_tns ) / deltaTns ));*/

  /* 3. this way seems to be best */
/* start_i=floor(0.5 + ( ( geoc_tns-block_tns )*16384/1e9)); */


  
   
  /* Fill frequency and phase arrays starting at time of injection NOT start
   * of block. */
  fData = output->f->data->data + start_i;
  phiData = output->phi->data->data + start_i;
  aData = output->a->data->data + start_i + start_i; 
  
  if ( !( strcmp( simRingdown->waveform, "Ringdown" ) ) )
  {
    
    for ( i = start_i; i < N_point_block; i++ )
    {
      /* bb = i - start_i; */
      t = (i - start_i) * dt;
      /* if ( bb <20 ) */
      /*  fprintf( stdout, "i = %d, bb = %d, t = %22.16e,",i,bb,t); */
      gtime = twopif0 / 2 / quality * t ;
      *(fData++)   = f0;
      *(phiData++) = twopif0 * t;
      /* if ( bb <20 ) */
      /*  fprintf( stdout, " phi = %22.16e,",*(phiData)); */
      *(aData++) = h0 * ( 1.0 + pow( cos( simRingdown->inclination ), 2 ) ) * 
        exp( - gtime );
      /* if ( bb <20 ) */
      /*  fprintf( stdout, " Aplus = %e,",*(aData)); */
      *(aData++) = h0* 2.0 * cos( simRingdown->inclination ) * exp( - gtime );
      /* if ( bb <20 ) */
      /*  fprintf( stdout, " Across = %e \n",*(aData)); */
      
    }
  }
  else
  {
    ABORT( stat, GENERATERINGH_ETYP, GENERATERINGH_MSGETYP );
  }

  /* Set output field and return. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="GenerateRingCP"> */
void
LALRingInjectSignals( 
    LALStatus               *stat, 
    REAL4TimeSeries         *series, 
    SimRingdownTable        *injections,
    COMPLEX8FrequencySeries *resp,
    INT4                     calType
    )
/* </lalVerbatim> */
{
  UINT4              k;
  INT4               injStartTime;
  INT4               injStopTime;
  DetectorResponse   detector;
  COMPLEX8Vector    *unity = NULL;
  CoherentGW         waveform;
  RingParamStruc     ringParam;
  REAL4TimeSeries    signal;
  SimRingdownTable  *simRingdown=NULL;
  LALDetector       *tmpDetector=NULL /*,*nullDetector=NULL*/;
  COMPLEX8FrequencySeries    *transfer = NULL;

  INITSTATUS( stat, "LALRingInjectSignals", GENERATERINGC );
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
  transfer = (COMPLEX8FrequencySeries *)
    LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
  if ( ! transfer ) 
  {
    ABORT( stat, GENERATERINGH_EMEM, GENERATERINGH_MSGEMEM );
  }
  memcpy( &(transfer->epoch), &(resp->epoch),
      sizeof(LIGOTimeGPS) );
  transfer->f0 = resp->f0;
  transfer->deltaF = resp->deltaF;

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
    LALUnitMultiply( stat->statusPtr, &(transfer->sampleUnits),
        &pair );
    CHECKSTATUSPTR( stat );
  }

  /* invert the response function to get the transfer function */
  LALCCreateVector( stat->statusPtr, &( transfer->data ),
      resp->data->length );
  CHECKSTATUSPTR( stat );

  LALCCreateVector( stat->statusPtr, &unity, resp->data->length );
  CHECKSTATUSPTR( stat );
  for ( k = 0; k < resp->data->length; ++k ) 
  {
    unity->data[k].re = 1.0;
    unity->data[k].im = 0.0;
  }

  LALCCVectorDivide( stat->statusPtr, transfer->data, unity,
      resp->data );
  CHECKSTATUSPTR( stat );

  LALCDestroyVector( stat->statusPtr, &unity );
  CHECKSTATUSPTR( stat );

  /* Set up a time series to hold signal in ADC counts */
  signal.deltaT = series->deltaT;
  if ( ( signal.f0 = series->f0 ) != 0 )
  {
    ABORT( stat, GENERATERINGH_EMEM, GENERATERINGH_MSGEMEM );
  }
  signal.sampleUnits = lalADCCountUnit;

  signal.data=NULL;
  LALSCreateVector( stat->statusPtr, &(signal.data), 
      series->data->length );
  CHECKSTATUSPTR( stat );

  /* loop over list of waveforms and inject into data stream */
  simRingdown = injections;
  while ( simRingdown )
  {
    /* only do the work if the ring is in injection zone */
    if( (injStartTime - simRingdown->geocent_start_time.gpsSeconds) *
        (injStopTime - simRingdown->geocent_start_time.gpsSeconds) > 0 )
    {
      simRingdown = simRingdown->next;
      continue;
    }

    /* set the ring params */
    ringParam.deltaT = series->deltaT;
    if( !( strcmp( simRingdown->coordinates, "HORIZON" ) ) )
    {
      ringParam.system = COORDINATESYSTEM_HORIZON;
    }
    else if ( !( strcmp( simRingdown->coordinates, "ZENITH" ) ) )
    {
      /* set coordinate system for completeness */
      ringParam.system = COORDINATESYSTEM_EQUATORIAL;
      detector.site = NULL;
    }
    else if ( !( strcmp( simRingdown->coordinates, "GEOGRAPHIC" ) ) )
    {
     ringParam.system = COORDINATESYSTEM_GEOGRAPHIC;
    }
    else if ( !( strcmp( simRingdown->coordinates, "EQUATORIAL" ) ) )
    {
      ringParam.system = COORDINATESYSTEM_EQUATORIAL;
    }
    else if ( !( strcmp( simRingdown->coordinates, "ECLIPTIC" ) ) )
    {
      ringParam.system = COORDINATESYSTEM_ECLIPTIC;
    }
    else if ( !( strcmp( simRingdown->coordinates, "GALACTIC" ) ) )
    {
      ringParam.system = COORDINATESYSTEM_GALACTIC;
    }
    else
      ringParam.system = COORDINATESYSTEM_EQUATORIAL;

    /* generate the ring */
    memset( &waveform, 0, sizeof(CoherentGW) );
    LALGenerateRing( stat->statusPtr, &waveform, series, simRingdown, &ringParam );
    CHECKSTATUSPTR( stat );

 
    /* print the waveform to a file */
    if ( 1 )
      {
        FILE *fp;
        char fname[512];
        UINT4 jj, kplus, kcross;
        LALSnprintf( fname, sizeof(fname) / sizeof(*fname), 
            "waveform-%d-%d-%s.txt", 
            simRingdown->geocent_start_time.gpsSeconds,
            simRingdown->geocent_start_time.gpsNanoSeconds,
            simRingdown->waveform );
        fp = fopen( fname, "w" );
         
        for( jj = 0, kplus = 0, kcross = 1; jj < waveform.phi->data->length; 
            ++jj, kplus += 2, kcross +=2 )
          {
            fprintf(fp, "%e %e %le %e\n", 
                waveform.a->data->data[kplus], 
                waveform.a->data->data[kcross], 
                waveform.phi->data->data[jj], 
                waveform.f->data->data[jj]);
            }
        fclose( fp );     
        }
    /* end */

    

    /* must set the epoch of signal since it's used by coherent GW */
    signal.epoch = waveform.a->epoch;
    memset( signal.data->data, 0, signal.data->length * sizeof(REAL4) );

    /* decide which way to calibrate the data; defaul to old way */
    if( calType )
      detector.transfer=NULL;
    else
      detector.transfer=transfer;
    
    /* convert this into an ADC signal */
    LALSimulateCoherentGW( stat->statusPtr, 
        &signal, &waveform, &detector );
    CHECKSTATUSPTR( stat );

    /* if calibration using RespFilt */
    if( calType == 1 )
      XLALRespFilt(&signal, transfer);

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
    simRingdown = simRingdown->next;
  }

  /* destroy the signal */
  LALSDestroyVector( stat->statusPtr, &(signal.data) );
  CHECKSTATUSPTR( stat );

  LALCDestroyVector( stat->statusPtr, &( transfer->data ) );
  CHECKSTATUSPTR( stat );

  if ( detector.site ) LALFree( detector.site );
  LALFree( transfer );

  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
