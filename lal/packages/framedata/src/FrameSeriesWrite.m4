ifelse(TYPE,`COMPLEX16',`define(`FRTYPE',`FR_VECT_16C')')
ifelse(TYPE,`COMPLEX8',`define(`FRTYPE',`FR_VECT_8C')')
ifelse(TYPE,`REAL8',`define(`FRTYPE',`FR_VECT_8R')')
ifelse(TYPE,`REAL4',`define(`FRTYPE',`FR_VECT_4R')')
ifelse(TYPE,`INT8',`define(`FRTYPE',`FR_VECT_8S')')
ifelse(TYPE,`INT4',`define(`FRTYPE',`FR_VECT_4S')')
ifelse(TYPE,`INT2',`define(`FRTYPE',`FR_VECT_2S')')

ifelse(TYPE,`COMPLEX16',`define(`FRDATA',`dataD')')
ifelse(TYPE,`COMPLEX8',`define(`FRDATA',`dataD')')
ifelse(TYPE,`REAL8',`define(`FRDATA',`dataD')')
ifelse(TYPE,`REAL4',`define(`FRDATA',`dataF')')
ifelse(TYPE,`INT8',`define(`FRDATA',`dataL')')
ifelse(TYPE,`INT4',`define(`FRDATA',`dataI')')
ifelse(TYPE,`INT2',`define(`FRDATA',`dataS')')

define(`STYPE',`format(`%sTimeSeries',TYPE)')
define(`FUNC',`format(`LALFrWrite%s',STYPE)')

/* <lalVerbatim file="FrameSeriesCP"> */
void
FUNC (
    LALStatus		*status,
    STYPE 	*series,
    FrOutPar		*params
    )
{ /* </lalVerbatim> */
  TYPE 	*data;
  CHAR   seconds[] = "s";
  CHAR   comment[] = "Created by FUNC $Id$";
  CHAR   source[256];
  CHAR   fname[256];
  CHAR   units[LALUnitNameSize];
  CHARVector vnits;
  struct FrFile *frfile;
  UINT4 nframes;
  INT8 t;
  INT8 tend;
  INT4 dt;

  INITSTATUS( status, "FUNC", FRAMESERIESC );  
  ASSERT( series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( params, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ATTATCHSTATUSPTR( status );

  vnits.length = sizeof( units );
  vnits.data = units;

  strncpy( source, params->source ? params->source : "F", sizeof( source ) );
  
  TRY( LALUnitAsString( status->statusPtr, &vnits, &series->sampleUnits ),
      status );

  t     = EPOCH_TO_I8TIME( series->epoch );
  tend  = t + (INT8)floor( 1e9 * series->data->length * series->deltaT );
  dt    = tend / (INT8)1000000000;
  dt   -= series->epoch.gpsSeconds;
  dt    = dt < 1 ? 1 : dt; /* must be at least one */

  LALSnprintf( fname, sizeof( fname ), "%s-%s-%d-%d.gwf", source,
      params->description ? params->description : "UNKNOWN",
      series->epoch.gpsSeconds, dt );
  frfile = FrFileONew( fname, 0 );

  data = series->data->data;
  nframes = params->nframes;
  while ( nframes-- > 0 )
  {
    UINT4 ncpy;
    struct FrameH *frame;
    struct FrVect *vect;
    if ( nframes )
    {
      ncpy = series->data->length / params->nframes;
    }
    else
    {
      ncpy = series->data->length - ( data - series->data->data );
    }
    frame = FrameHNew( source );
    frame->run = params->run;
    frame->frame = params->frame++;
    frame->GTimeS = t / (INT8)1000000000;
    frame->GTimeN = t % (INT8)1000000000;
    frame->dt = ncpy * series->deltaT;
#if !defined FR_VERS || FR_VERS < 5000
    frame->localTime = 0;
#endif

    /*** CHECK FSHIFT ***/
    vect = makeFrVect1D( frame, params->type, series->name, comment, seconds,
        units, FRTYPE, series->deltaT ? 1.0 / series->deltaT : 0.0,
        series->f0, series->deltaT, ncpy );
    if ( ! vect )
    {
      ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
    }
    memcpy( vect->FRDATA, data, ncpy * sizeof( *series->data->data ) );

    FrameWrite( frame, frfile );
    data += ncpy;
    t += 1e9 * ncpy * series->deltaT;
  }

  FrFileOEnd( frfile );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
