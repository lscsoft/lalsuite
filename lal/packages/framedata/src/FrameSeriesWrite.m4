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
define(`FSTYPE',`format(`%sFrequencySeries',TYPE)')
define(`FUNC',`format(`LALFrWrite%s',STYPE)')
define(`FSFUNC',`format(`LALFrWrite%s',FSTYPE)')
define(`XFUNC',`format(`XLALFrWrite%s',STYPE)')
define(`XFSFUNC',`format(`XLALFrWrite%s',FSTYPE)')


int XFUNC ( STYPE *series, int frnum )
{
  static const char func[] = "XFUNC";
  char fname[FILENAME_MAX];
  char tmpfname[FILENAME_MAX];
  char comment[] = "XFUNC $Id$";
  char seconds[] = "s";
  char units[LALUnitTextSize];
  const int run = 0; /* this routine always sets run number to zero */
  char site = 0;
  char ifo[3] = "\0\0\0";
  char description[FILENAME_MAX] = "";
  char *s;
  int detector;
  int detectorFlag = 0;
  REAL8 duration;
  INT4 t0;
  INT4 dt;
  struct FrProcData *proc;
  struct FrVect *vect;
  struct FrameH *frame;
  struct FrFile *frfile;
  int c;

  if ( ! series )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! series->data || ! series->data->data )
    XLAL_ERROR( func, XLAL_EINVAL );

  /* if channel name is of the form Xn:... then Xn is the ifo name */
  if ( isupper( series->name[0] ) )
  {
    site = series->name[0];
    if ( isdigit( series->name[1] ) && series->name[2] == ':' )
    {
      memcpy( ifo, series->name, 2 );
      /* figure out what detector this is */
      for ( detector = 0; detector < LAL_NUM_DETECTORS; ++detector )
        if ( ! strcmp( ifo, lalCachedDetectors[detector].frDetector.prefix ) )
        {
          detectorFlag = 1 << 2 * detector;
          break;
        }
    }
  }

  /* if site was not set, do something sensible with it */
  if ( ! site )
    site = 'F'; /* fake site */

  /* transform channel name to be a valid description field of file name */
  memset( description, 0, sizeof( description ) );
  strncpy( description, series->name, sizeof( description) - 1 );
  s = description;
  while ( *s )
  {
    if ( ! isalnum( *s ) )
      *s = '_'; /* replace invalid description character with underscore */
    ++s;
  }

  /* compute duration of the data and the stop and end times of frame file */
  duration = series->data->length * series->deltaT;
  t0 = series->epoch.gpsSeconds;
  dt = (INT4)ceil( XLALGPSGetREAL8( &series->epoch ) + duration ) - t0;
  if ( t0 < 0 || dt < 1 )
    XLAL_ERROR( func, XLAL_ETIME );

  /* construct file name */
  c = snprintf( fname, sizeof( fname ), "%c-%s-%d-%d.gwf", site, description, t0, dt );
  if ( c < 0 || c > (int)sizeof(fname) - 2 )
    XLAL_ERROR( func, XLAL_ENAME );
  c = snprintf( tmpfname, sizeof( tmpfname ), "%s.tmp", fname );
  if ( c < 0 || c > (int)sizeof(tmpfname) - 2 )
    XLAL_ERROR( func, XLAL_ENAME );

  /* get sample unit string */
  if (NULL == XLALUnitAsString( units, sizeof( units ), &series->sampleUnits ))
    XLAL_ERROR( func, XLAL_EFUNC );

  /* construct vector and copy data */
  vect = FrVectNew1D( series->name, FRTYPE, series->data->length, series->deltaT, seconds, units );
  if ( ! vect )
    XLAL_ERROR( func, XLAL_EERR );
  vect->startX[0] = 0.0;
  memcpy( vect->data, series->data->data, series->data->length * sizeof( *series->data->data ) );

  frame = XLALFrameNew( &series->epoch, duration, ifo, run, frnum, detectorFlag );
  if ( ! frame )
  {
    FrVectFree( vect );
    XLAL_ERROR( func, XLAL_EFUNC );
  }

  proc = FrProcDataNewV( frame, vect );
  if ( ! proc )
  {
    FrameFree( frame );
    FrVectFree( vect );
    XLAL_ERROR( func, XLAL_EERR );
  }
  FrStrCpy( &proc->comment, comment );
  proc->timeOffset = 0;
  proc->type       = 1;
  proc->subType    = 0;
  proc->tRange     = duration;
  proc->fShift     = series->f0;
  proc->phase      = 0.0;
  proc->BW         = 0.0;

  /* write first to tmpfile then rename it */
  frfile = FrFileONew( tmpfname, 0 );
  if ( ! frfile )
  {
    FrameFree( frame ); /* this frees proc and vect */
    XLAL_ERROR( func, XLAL_EIO );
  }
  if ( FR_OK != FrameWrite( frame, frfile ) )
  {
    FrameFree( frame ); /* this frees proc and vect */
    FrFileFree( frfile );
    XLAL_ERROR( func, XLAL_EERR );
  }

  FrFileOEnd( frfile );
  FrameFree( frame ); /* this frees proc and vect */

  /* now rename tmpfile */
  if ( rename( tmpfname, fname ) < 0 )
    XLAL_ERROR( func, XLAL_ESYS );

  return 0;
}



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
  CHAR   source[FILENAME_MAX];
  CHAR   fname[FILENAME_MAX];
  CHAR   tmpfname[FILENAME_MAX];
  CHAR   units[LALUnitTextSize];
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

  snprintf( tmpfname, sizeof( tmpfname ), "%s-%s-%d-%d.gwf.tmp", source,
      params->description ? params->description : "UNKNOWN",
      series->epoch.gpsSeconds, dt );
  snprintf( fname, sizeof( fname ), "%s-%s-%d-%d.gwf", source,
      params->description ? params->description : "UNKNOWN",
      series->epoch.gpsSeconds, dt );
  frfile = FrFileONew( tmpfname, 0 );

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

      FrVectFree(vect);
      vect=NULL;

      ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
    }
    memcpy( vect->FRDATA, data, ncpy * sizeof( *series->data->data ) );

    FrameWrite( frame, frfile );
    data += ncpy;
    t += 1e9 * ncpy * series->deltaT;

    FrVectFree(vect);
    vect=NULL;

  }

  FrFileOEnd( frfile );
  rename( tmpfname, fname );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}




/* <lalVerbatim file="FrameSeriesCP"> */
void
FSFUNC (
    LALStatus		*status,
    FSTYPE 	*series,
    FrOutPar		*params,
    INT4                 subtype
    )
{ /* </lalVerbatim> */
  TYPE 	*data;
  CHAR   hertz[] = "Hz";
  CHAR   comment[] = "Created by FSFUNC $Id$";
  CHAR   source[FILENAME_MAX];
  CHAR   fname[FILENAME_MAX];
  CHAR   tmpfname[FILENAME_MAX];
  CHAR   units[LALUnitTextSize];
  CHARVector vnits;
  struct FrFile *frfile;
  UINT4 nframes;
  REAL4 deltaT;
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

  deltaT = series->deltaF ? 1.0 / series->deltaF : 1;
  t     = EPOCH_TO_I8TIME( series->epoch );
  tend  = t + (INT8)floor( 1e9 * deltaT );
  dt    = tend / (INT8)1000000000;
  dt   -= series->epoch.gpsSeconds;
  dt    = dt < 1 ? 1 : dt; /* must be at least one */

  snprintf( tmpfname, sizeof( tmpfname ), "%s-%s-%d-%d.gwf.tmp", source,
      params->description ? params->description : "UNKNOWN",
      series->epoch.gpsSeconds, dt );
  snprintf( fname, sizeof( fname ), "%s-%s-%d-%d.gwf", source,
      params->description ? params->description : "UNKNOWN",
      series->epoch.gpsSeconds, dt );
  frfile = FrFileONew( tmpfname, 0 );

  data = series->data->data;
  /* nframes = params->nframes; */
  nframes = 1; /* params->nframes is ignored */
  {
    UINT4 ncpy;
    struct FrameH     *frame;
    struct FrVect     *vect;
    struct FrProcData *proc;
    ncpy = series->data->length;
    frame = FrameHNew( source );
    frame->run = params->run;
    frame->frame = params->frame++;
    frame->GTimeS = t / (INT8)1000000000;
    frame->GTimeN = t % (INT8)1000000000;
    frame->dt = deltaT;
#if !defined FR_VERS || FR_VERS < 5000
    frame->localTime = 0;
#endif

    /* NOTE: always make a proc channel! */
    proc = calloc( 1, sizeof( *proc ) );
    if ( ! proc )
    {
      ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
    }
    vect = FrVectNew1D( series->name, FRTYPE, series->data->length,
        series->deltaF, hertz, units );
    if ( ! vect )
    {
      free( proc );

      FrVectFree(vect);
      vect=NULL;

      ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
    }
    vect->startX[0] = series->f0;

    FrStrCpy( &proc->name, series->name );
    FrStrCpy( &proc->comment, comment );
    proc->next = frame->procData;
    frame->procData = proc;
    proc->classe = FrProcDataDef();
    proc->type = 2;
    proc->data = vect;
    proc->subType = subtype;
    proc->tRange = deltaT;
    proc->fRange = series->data->length * series->deltaF;
    /* frequency shift is in vector */

    memcpy( vect->FRDATA, data, ncpy * sizeof( *series->data->data ) );

    FrameWrite( frame, frfile );

    FrVectFree(vect);
    vect=NULL;

  }

  FrFileOEnd( frfile );
  rename( tmpfname, fname );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
