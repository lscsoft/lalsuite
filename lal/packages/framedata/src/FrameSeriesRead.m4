ifelse(TYPE,`REAL8',`define(`FRTYPE',`FR_VECT_8R')')
ifelse(TYPE,`REAL4',`define(`FRTYPE',`FR_VECT_4R')')
ifelse(TYPE,`INT8',`define(`FRTYPE',`FR_VECT_8S')')
ifelse(TYPE,`INT4',`define(`FRTYPE',`FR_VECT_4S')')
ifelse(TYPE,`INT2',`define(`FRTYPE',`FR_VECT_2S')')

ifelse(TYPE,`REAL8',`define(`FRDATA',`dataD')')
ifelse(TYPE,`REAL4',`define(`FRDATA',`dataF')')
ifelse(TYPE,`INT8',`define(`FRDATA',`dataL')')
ifelse(TYPE,`INT4',`define(`FRDATA',`dataI')')
ifelse(TYPE,`INT2',`define(`FRDATA',`dataS')')

define(`STYPE',`format(`%sTimeSeries',TYPE)')
define(`FUNC',`format(`LALFrGet%s',STYPE)')

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
  UINT4		 ncpy;
  TYPE 		*dest;
  INT8		 tnow;
  INT8		 tbeg;
  REAL8          rate;

  INITSTATUS( status, "FUNC", FRAMESERIESC );  

  ASSERT( series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( stream->err )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  if ( stream->end )
  {
    ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
  }

  strncpy( series->name, chanin->name, sizeof( series->name ) );
  vect = loadFrVect( stream->frame, series->name, chanin->type );
  if ( ! vect || ! vect->data )
  {
    ABORT( status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN );
  }
  if ( vect->type != FRTYPE )
  {
    ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );
  }

  tnow = EPOCH_TO_I8TIME( stream->epoch );
  tbeg = SECNAN_TO_I8TIME( vect->GTimeS, vect->GTimeN );
  if ( tnow < tbeg )
  {
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
    RETURN( status );
  }
  ASSERT( series->data->data, status, FRAMESTREAMH_ENULL,
      FRAMESTREAMH_MSGENULL );
  ASSERT( series->data->length > 0, status, FRAMESTREAMH_ESIZE,
      FRAMESTREAMH_MSGESIZE );

  ATTATCHSTATUSPTR( status );

  dest = series->data->data;
  need = series->data->length;
  if ( noff > vect->nData )
  {
    ABORT( status, FRAMESTREAMH_ETIME, FRAMESTREAMH_MSGETIME );
  }

  /* number of points to copy */
  ncpy = ( vect->nData - noff < need ) ? vect->nData - noff : need;
  memcpy( dest, vect->FRDATA + noff, ncpy * sizeof( *series->data->data ) );
  dest += ncpy;
  need -= ncpy;

  /* if still data remaining */
  while ( need )
  {
    LALFrNext( status->statusPtr, stream );
    BEGINFAIL( status )
    {
      memset( dest, 0, need * sizeof( *series->data->data ) );
    }
    ENDFAIL( status );
    if ( stream->end )
    {
      memset( dest, 0, need * sizeof( *series->data->data ) );
      ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
    }

    /* load more data */
    vect = loadFrVect( stream->frame, series->name, chanin->type );
    if ( ! vect || ! vect->data )
    {
      memset( dest, 0, need * sizeof( *series->data->data ) );
      ABORT( status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN );
    }
    if ( vect->type != FRTYPE )
    {
      memset( dest, 0, need * sizeof( *series->data->data ) );
      ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );
    }

    if ( stream->gap ) /* gap in data */
    {
      dest = series->data->data;
      need = series->data->length;
      series->epoch.gpsSeconds = vect->GTimeS;
      series->epoch.gpsNanoSeconds = vect->GTimeN;
    }

    /* copy data */
    ncpy = vect->nData < need ? vect->nData : need;
    memcpy( dest, vect->FRDATA, ncpy * sizeof( *series->data->data ) );
    dest += ncpy;
    need -= ncpy;
  }

  /* update stream start time very carefully:
   * start time must be the exact time of the next sample, rounded to the
   * nearest nanosecond */
  SET_EPOCH( &stream->epoch, EPOCH_TO_I8TIME( series->epoch )
      + (INT8)floor( 1e9 * series->data->length * series->deltaT + 0.5 ) );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
