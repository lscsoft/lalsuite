ifelse(TYPE,`COMPLEX16',`define(`FRTYPE',`FR_VECT_16C')')
ifelse(TYPE,`COMPLEX8',`define(`FRTYPE',`FR_VECT_8C')')
ifelse(TYPE,`REAL8',`define(`FRTYPE',`FR_VECT_8R')')
ifelse(TYPE,`REAL4',`define(`FRTYPE',`FR_VECT_4R')')
ifelse(TYPE,`INT8',`define(`FRTYPE',`FR_VECT_8S')')
ifelse(TYPE,`INT4',`define(`FRTYPE',`FR_VECT_4S')')
ifelse(TYPE,`INT2',`define(`FRTYPE',`FR_VECT_2S')')

ifelse(TYPE,`COMPLEX16',`define(`FRDATA',`dataD')')
ifelse(TYPE,`COMPLEX8',`define(`FRDATA',`dataF')')
ifelse(TYPE,`REAL8',`define(`FRDATA',`dataD')')
ifelse(TYPE,`REAL4',`define(`FRDATA',`dataF')')
ifelse(TYPE,`INT8',`define(`FRDATA',`dataL')')
ifelse(TYPE,`INT4',`define(`FRDATA',`dataI')')
ifelse(TYPE,`INT2',`define(`FRDATA',`dataS')')

define(`STYPE',`format(`%sTimeSeries',TYPE)')
define(`FSTYPE',`format(`%sFrequencySeries',TYPE)')
define(`FUNC',`format(`LALFrGet%s',STYPE)')
define(`FUNCM',`format(`LALFrGet%sMetadata',STYPE)')
define(`FSFUNC',`format(`LALFrGet%s',FSTYPE)')

/* <lalVerbatim file="FrameSeriesCP"> */
void
FSFUNC (
    LALStatus		*status,
    FSTYPE	*series,
    FrChanIn		*chanin,
    FrStream		*stream
    )
{ /* </lalVerbatim> */
  struct FrVect	*vect;
  INITSTATUS( status, "FSFUNC", FRAMESERIESC );  

  ASSERT( series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( ! series->data, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  if ( stream->state & LAL_FR_END )
  {
    ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
  }

  strncpy( series->name, chanin->name, sizeof( series->name ) );
  vect = loadFrVect( stream, series->name );
  if ( ! vect || ! vect->data )
  {
    ABORT( status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN );
  }
  if ( vect->type != FRTYPE )
  {
    ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );
  }

#if defined FR_VERS && FR_VERS >= 5000
  series->epoch.gpsSeconds     = floor( vect->GTime );
  series->epoch.gpsNanoSeconds = floor( 1e9 * ( vect->GTime - floor( vect->GTime ) ) );
#else
  series->epoch.gpsSeconds     = vect->GTimeS;
  series->epoch.gpsNanoSeconds = vect->GTimeN;
#endif
  series->deltaF = vect->dx[0];
  series->f0 = 0; /* FIXME: should get correct value... */
  series->sampleUnits = lalADCCountUnit;
  series->data = LALCalloc( 1, sizeof( *series->data ) );
  if ( ! series->data )
  {
    ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
  }
  series->data->length = vect->nData;
  series->data->data = LALMalloc( series->data->length * sizeof( *series->data->data ) );
  if ( ! series->data->data )
  {
    ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
  }
  memcpy( series->data->data, vect->FRDATA, series->data->length * sizeof( *series->data->data ) );

  FrVectFree(vect);

  vect=NULL;

  RETURN( status );
}


/* <lalVerbatim file="FrameSeriesCP"> */
void
FUNC (
    LALStatus		*status,
    STYPE	*series,
    FrChanIn		*chanin,
    FrStream		*stream
    )
{ /* </lalVerbatim> */
  const REAL8    fuzz = 0.1 / 16384.0; /* smallest discernable unit of time */
  struct FrVect	*vect;
  UINT4		 need;
  UINT4		 noff;
  UINT4          mult;
  UINT4		 ncpy;
  TYPE 		*dest;
  INT8		 tnow;
  INT8		 tbeg;
  INT8           tend;
  REAL8          rate;
  INT4           gap = 0;

  INITSTATUS( status, "FUNC", FRAMESERIESC );  

  ASSERT( series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  if ( stream->state & LAL_FR_END )
  {
    ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
  }

  strncpy( series->name, chanin->name, sizeof( series->name ) );
  vect = loadFrVect( stream, series->name );
  if ( ! vect || ! vect->data )
  {
    ABORT( status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN );
  }
  if ( vect->type != FRTYPE )
  {
    FrVectFree(vect);	
    ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );
  }

  tnow = EPOCH_TO_I8TIME( stream->epoch );
#if defined FR_VERS && FR_VERS >= 5000
  tbeg = 1e9 * vect->GTime;
#else
  tbeg = SECNAN_TO_I8TIME( vect->GTimeS, vect->GTimeN );
#endif
  if ( tnow + 1000 < tbeg )  /* added 1000 ns to account for double precision */
  {
    FrVectFree(vect);
    ABORT( status, FRAMESTREAMH_ETIME, FRAMESTREAMH_MSGETIME );
  }

  /* compute number of points offset very carefully:
   * if current time is within fuzz of a sample, get that sample;
   * otherwise get the sample just after the requested time */
  rate = vect->dx[0] ? 1.0 / vect->dx[0] : 0.0;
  noff = ceil( ( 1e-9 * ( tnow - tbeg ) - fuzz ) * rate ); 

  /* adjust current time to be exactly the first sample
   * (rounded to nearest nanosecond) */
  tnow = tbeg + floor( 1e9 * noff * vect->dx[0] + 0.5 );


  SET_EPOCH( &series->epoch, tnow );
  series->deltaT = vect->dx[0];
  series->sampleUnits = lalADCCountUnit;

  if ( ! series->data ) /* no data requested: return now */
  {
    FrVectFree(vect);
    RETURN( status );
  }
  ASSERT( series->data->data, status, FRAMESTREAMH_ENULL,
      FRAMESTREAMH_MSGENULL );
  ASSERT( series->data->length > 0, status, FRAMESTREAMH_ESIZE,
      FRAMESTREAMH_MSGESIZE );

  ATTATCHSTATUSPTR( status );

  /* mult is two if output series is complex */
  mult = sizeof( *series->data->data ) / sizeof( *vect->FRDATA );
  dest = series->data->data;
  need = series->data->length;
  if ( noff > vect->nData )
  {
    FrVectFree(vect);
    ABORT( status, FRAMESTREAMH_ETIME, FRAMESTREAMH_MSGETIME );
  }

  /* number of points to copy */
  ncpy = ( vect->nData - noff < need ) ? ( vect->nData - noff ) : need;
  memcpy( dest, vect->FRDATA + noff * mult, ncpy * sizeof( *series->data->data ) );

  FrVectFree(vect);
  vect=NULL;

  dest += ncpy;
  need -= ncpy;


  /* if still data remaining */
  while ( need )
  {
    LALFrNext( status->statusPtr, stream );
    BEGINFAIL( status )
    {
      if(vect) FrVectFree(vect);	
      memset( dest, 0, need * sizeof( *series->data->data ) );
    }
    ENDFAIL( status );
    if ( stream->state & LAL_FR_END )
    {
      if(vect) FrVectFree(vect);	
      memset( dest, 0, need * sizeof( *series->data->data ) );
      ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
    }

    /* load more data */
    vect = loadFrVect( stream, series->name );
    if ( ! vect || ! vect->data )
    {
      memset( dest, 0, need * sizeof( *series->data->data ) );
      ABORT( status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN );
    }
    if ( vect->type != FRTYPE )
    {
      FrVectFree(vect);
      memset( dest, 0, need * sizeof( *series->data->data ) );
      ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );
    }

    if ( stream->state & LAL_FR_GAP ) /* gap in data */
    {
      dest = series->data->data;
      need = series->data->length;
#if defined FR_VERS && FR_VERS >= 5000
      series->epoch.gpsSeconds = floor( vect->GTime );
      series->epoch.gpsNanoSeconds = (INT8)( 1e9 * vect->GTime ) % (INT8)1000000000;
#else
      series->epoch.gpsSeconds = vect->GTimeS;
      series->epoch.gpsNanoSeconds = vect->GTimeN;
#endif
      gap = 1; /* need to record this because next FrNext will erase! */
    }

    /* copy data */
    ncpy = vect->nData < need ? vect->nData : need;
    memcpy( dest, vect->FRDATA, ncpy * sizeof( *series->data->data ) );


  FrVectFree(vect);
  vect=NULL;


    dest += ncpy;
    need -= ncpy;
  }

  /* update stream start time very carefully:
   * start time must be the exact time of the next sample, rounded to the
   * nearest nanosecond */
  SET_EPOCH( &stream->epoch, EPOCH_TO_I8TIME( series->epoch )
      + (INT8)floor( 1e9 * series->data->length * series->deltaT + 0.5 ) );

  /* check to see that we are still within current frame */
  tend  = SECNAN_TO_I8TIME( stream->file->toc->GTimeS[stream->pos], stream->file->toc->GTimeN[stream->pos] );
  tend += (INT8)floor( 1e9 * stream->file->toc->dt[stream->pos] );
  if ( tend <= EPOCH_TO_I8TIME( stream->epoch ) )
  {
    LIGOTimeGPS keep;
    keep = stream->epoch;
    /* advance a frame */
    /* failure is benign, so we return results */
    TRY( LALFrNext( status->statusPtr, stream ), status );
    if ( ! stream->state & LAL_FR_GAP )
    {
      stream->epoch = keep;
    }
  }

  if ( gap ) /* there was a gap in the data */
  {
    stream->state |= LAL_FR_GAP;
  }

  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }

  /* remove this: the error will be reported on the *next* call! */
  /*
  if ( stream->state & LAL_FR_END )
  {
    ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
  }
  */

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameSeriesCP"> */
void
FUNCM (
    LALStatus		*status,
    STYPE	*series,
    FrChanIn		*chanin,
    FrStream		*stream
    )
{ /* </lalVerbatim> */
  void *sequence;
                                                                                
  INITSTATUS (status, "FUNCM", FRAMESERIESC);
  ATTACHSTATUSPTR (status);

  ASSERT (series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
  ASSERT (stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL);
                                                                                
  /* save the sequence address, then wipe the series structure */
  sequence = series->data;
  memset (series, 0, sizeof(*series));

  /* call FUNC to populate the series' metadata */
  FUNC (status, series, chanin, stream);
  CHECKSTATUSPTR (status);

  /* restore the sequence address */
  series->data = sequence;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
