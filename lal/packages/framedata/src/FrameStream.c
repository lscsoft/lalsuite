/**** <lalVerbatim file="FrameStreamCV">
 * Author: Jolien D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{FrameStream.c}}
 *
 * These are the low-level functions for manupulating a frame stream.
 *
 * \subsubsection*{Prototypes}
 * \input{FrameStreamCP}
 * \idx{FrOpen}
 * \idx{FrClose}
 * \idx{FrEnd}
 * \idx{FrNext}
 * \idx{FrRewind}
 * \idx{FrSeek}
 * \idx{FrTell}
 * \idx{FrGetPos}
 * \idx{FrSetPos}
 * \idx{FrGetINT2TimeSeries}
 * \idx{FrGetREAL4TimeSeries}
 * \idx{FrWriteREAL4TimeSeries}
 *
 * \subsubsection*{Description}
 *
 * Many of these routines perform functions that are similar to standard C
 * file stream manipulation routines.  The names have been chosen to also be
 * similar to the standard C routines.
 *
 * The routines \texttt{LALFrOpen()} and \texttt{LALFrClose()} are used to open
 * and close a frame stream.  The stream is created by \texttt{LALFrOpen()},
 * and must be a pointer to \texttt{NULL} before it is opened.  It must have
 * been created prior to calling \texttt{LALFrClose()}, and after this call,
 * the stream will be a pointer to \texttt{NULL}.  The routine
 * \texttt{LALFrOpen()} requires the user to specify the directory name of the
 * frame files and the head names.  If the directory is \texttt{NULL}, the
 * routine uses the current director (\texttt{.}).  The head names specifies
 * which files are the wanted files in the specified directory.  Wildcards are
 * allowed.  For example, to get LLO frames only, the head names could be set
 * to \texttt{L-*.F}.  If the head name is \texttt{NULL}, the default value
 * \texttt{*.F} is used.
 *
 * The routine \texttt{LALFrEnd()} determines if the end-of-frame-data flag for
 * the data stream has been set.
 *
 * The routine \texttt{LALFrNext()} advances the frame stream to the
 * beginning of the next frame.
 *
 * The routine \texttt{LALFrRewind()} rewinds the frame stream to the first
 * frame.
 *
 * The routine \texttt{LALFrSeek()} sets the frame stream to a specified time,
 * or the earliest time after the specified time if that time is not available
 * (e.g., if it is before the beginning of the frame stream or if it is in a
 * gap in the frame data).  The routine \texttt{LALFrTell()} returns the
 * current time within the frame stream.
 *
 * The routine \texttt{LALFrGetPos()} returns a structure containing the
 * current frame stream position.  The frame stream can later be restored to
 * this position using \texttt{LALFrSetPos()}.
 *
 * The routines
 * \texttt{LALFrGet}$\langle\mbox{datatype}\rangle$\texttt{TimeSeries()}
 * search the frame for a specified channel.  If the time series supplied has
 * data storage allocated, then the specified amount of data is filled from
 * the frame stream.  If no space has been allocated (so that the data field
 * is \texttt{NULL}), then only the channel information is returned in the
 * time series (e.g., the start time of the next data and the time step size).
 *
 * The routines
 * \texttt{LALFrWrite}$\langle\mbox{datatype}\rangle$\texttt{TimeSeries()}
 * outputs a given time series as a new frame file with the filename
 * $\langle\mbox{prefix}\rangle$\texttt{-}$\langle\mbox{GPS-seconds}\rangle$%
 * \texttt{-}$\langle\mbox{number-of-frames}\rangle$\texttt{-}%
 * $\langle\mbox{seconds-per-frame}\rangle$\texttt{.F} (or
 * $\langle\mbox{prefix}\rangle$\texttt{-}$\langle\mbox{gpstime}\rangle$%
 * \texttt{.F} if there is only one one-second frame) where prefix is the
 * specified frame prefix, GPS-seconds is the start time of the data in
 * seconds since 0h UTC 6 Jan 1980 (GPS time origin), number-of-frames is the
 * specified number of frames in the file, and seconds-per-frame is the
 * computed seconds per frame in the file (given the amount of data for output).
 *
 * \vfill{\footnotesize\input{FrameStreamCV}}
 *
 **** </lalLaTeX> */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <fnmatch.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <FrameL.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/FrameStream.h>

NRCSID( FRAMESTREAMC, "$Id$" );


#define SECNAN_TO_I8TIME( sec, nan ) \
  ((INT8)1000000000*(INT8)(sec)+(INT8)(nan))
/* Dangerous!!!: */
#define EPOCH_TO_I8TIME( epoch ) \
  SECNAN_TO_I8TIME( (epoch).gpsSeconds, (epoch).gpsNanoSeconds )
#define SET_EPOCH( pepoch, i8time ) \
  do { INT8 t=(i8time); LIGOTimeGPS *pe=(pepoch); \
    pe->gpsSeconds=t/(INT8)1000000000; pe->gpsNanoSeconds=t%(INT8)1000000000; \
  } while( 0 )


struct
tagFrStream
{
  CHAR          **filelist;
  UINT4           filenum;
  struct FrFile  *frfile;
  UINT4           frnum;
  struct FrameH  *frame;
  LIGOTimeGPS     epoch;
  INT4            end;
  INT4            err;
  INT4            gap;
};


static int strsort( const void *p1, const void *p2 )
{
  return strcmp( *((char * const *)p1), *((char * const *)p2) );
}

static int list_files( char ***flist, const char *dirname, const char *pattern )
{
  DIR           *dir;
  struct dirent *ent;
  int            nfile;
  size_t         dirnamelen;

  pattern = pattern ? pattern : "*";
  dirname = dirname ? dirname : ".";
  dirnamelen = strlen( dirname );
  if ( ! ( dir = opendir( dirname ) ) )
    return -1;

  nfile = 0;
  while ( ( ent = readdir( dir ) ) )
    if ( ! fnmatch( pattern, ent->d_name, 0 ) )
    {
      size_t size = strlen( ent->d_name ) + dirnamelen + 2;
      struct stat buf;
      char *fname;
      fname = LALMalloc( size );
      LALSnprintf( fname, size, "%s/%s", dirname, ent->d_name );
      stat( fname, &buf );
      if ( S_ISREG( buf.st_mode ) )
        ++nfile;
      LALFree( fname );
    }

  if ( ! nfile )
  {
    if ( closedir( dir ) )
      return -1;
    return 0;
  }

  rewinddir( dir );
  *flist = LALCalloc( nfile + 1, sizeof( **flist ) );
  nfile = 0;
  while ( ( ent = readdir( dir ) ) )
    if ( ! fnmatch( pattern, ent->d_name, 0 ) )
    {
      size_t size = strlen( ent->d_name ) + dirnamelen + 2;
      struct stat buf;
      char *fname;
      fname = LALMalloc( size );
      LALSnprintf( fname, size, "%s/%s", dirname, ent->d_name );
      stat( fname, &buf );
      if ( S_ISREG( buf.st_mode ) )
        (*flist)[nfile++] = fname;
      else
        LALFree( fname );
    }

  qsort( *flist, nfile, sizeof( **flist ), strsort );

  if ( closedir( dir ) )
    return -1;
  return nfile;
}  


/*
 *
 * The following routines are designed to manipulate an input frame stream
 * as much like standard C input streams as possible.  Hence the names are
 * based on the standard C stream manipulation functions.
 *
 */


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrOpen(
    LALStatus    *status,
    FrStream    **stream,
    const CHAR   *dirname,
    const CHAR   *pattern
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALFrOpen", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( ! *stream, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL );

  *stream = LALCalloc( 1, sizeof( **stream ) );
  if ( ! *stream )
  {
    ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
  }

  if ( 1 > list_files( &(*stream)->filelist, dirname,
        pattern ? pattern : "*.F" ) )
  {
    LALFree( *stream );
    *stream = NULL;
    ABORT( status, FRAMESTREAMH_EFILE, FRAMESTREAMH_MSGEFILE );
  }

  FrLibIni( NULL, stderr, 0 );

  (*stream)->frfile = FrFileINew( (*stream)->filelist[0] );
  if ( ! (*stream)->frfile )
  {
    CHAR **tmp = (*stream)->filelist;
    while ( *tmp )
      LALFree( *tmp++ );
    LALFree( (*stream)->filelist );
    LALFree( *stream );
    *stream = NULL;
    ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
  }

  (*stream)->frame = FrameRead( (*stream)->frfile );
  if ( ! (*stream)->frame )
  {
    CHAR **tmp = (*stream)->filelist;
    FrFileOEnd( (*stream)->frfile );
    while ( *tmp )
      LALFree( *tmp++ );
    LALFree( (*stream)->filelist );
    LALFree( *stream );
    *stream = NULL;
    ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
  }

  (*stream)->epoch.gpsSeconds     = (*stream)->frame->GTimeS;
  (*stream)->epoch.gpsNanoSeconds = (*stream)->frame->GTimeN;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrClose(
    LALStatus  *status,
    FrStream  **stream
    )
{ /* </lalVerbatim> */
  CHAR **tmp;
  INITSTATUS( status, "LALFrClose", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( *stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  FrameFree( (*stream)->frame );
  FrFileIEnd( (*stream)->frfile );
  tmp = (*stream)->filelist;
  while ( *tmp )
    LALFree( *tmp++ );
  LALFree( *stream );
  *stream = NULL;
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrEnd(
    LALStatus *status,
    INT4      *end,
    FrStream  *stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALFrEnd", FRAMESTREAMC );  
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( end, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  *end = stream->end;
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrNext(
    LALStatus *status,
    FrStream  *stream
    )
{ /* </lalVerbatim> */
  /* timing accuracy: tenth of a sample interval for a 16kHz fast channel */
  const INT8 tacc = (INT8)floor( 0.1 * 1e9 / 16384.0 );
  INT8 texp;
  INT8 tact;
  INITSTATUS( status, "LALNextFrame", FRAMESTREAMC );  
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( stream->err )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  if ( stream->end )
  {
    ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
  }
  stream->gap = 0;

  /* compute expected start time of the next frame */
  texp  = SECNAN_TO_I8TIME( stream->frame->GTimeS, stream->frame->GTimeN );
  texp += (INT8)floor( 1e9 * stream->frame->dt );

  FrameFree( stream->frame );
  stream->frame = FrameRead( stream->frfile );
  ++stream->frnum;

  if ( ! stream->frame ) /* no more frames in file */
  {
    FrFileOEnd( stream->frfile );
    stream->frfile = NULL;
    stream->frnum = 0;
    ++stream->filenum;

    /* no more frame files: close stream */
    if ( ! stream->filelist[stream->filenum] )
    {
      stream->end = EOF;
      RETURN( status );
    }

    stream->frfile = FrFileINew( stream->filelist[stream->filenum] );
    if ( ! stream->frfile )
    {
      stream->err = FRAMESTREAMH_EOPEN; /* error: couldn't read frame file! */
      ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
    }

    stream->frame = FrameRead( stream->frfile );
    if ( ! stream->frame )
    {
      stream->err = FRAMESTREAMH_EREAD; /* error: empty frame file! */
      ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
    }
  }

  /* compute actual start time of this new frame */
  tact = SECNAN_TO_I8TIME( stream->frame->GTimeS, stream->frame->GTimeN );
  stream->epoch.gpsSeconds     = stream->frame->GTimeS;
  stream->epoch.gpsNanoSeconds = stream->frame->GTimeN;

  if ( abs( texp - tact ) > tacc ) /* there is a gap */
  {
    LALInfo( status, "Gap in frame data." );
    stream->gap = 1;
  }

  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrRewind( 
    LALStatus *status,
    FrStream  *stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALFrRewind", FRAMESTREAMC );  
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( stream->frame )
  {
    FrameFree( stream->frame );
  }
  if ( stream->frfile )
  {
    FrFileOEnd( stream->frfile );
  }
  stream->filenum = 0;
  stream->frnum   = 0;
  stream->end     = 0;
  stream->err     = 0;
  stream->gap     = 0;
  stream->frfile  = FrFileINew( stream->filelist[0] );
  if ( ! stream->frfile )
  {
    stream->err = FRAMESTREAMH_EOPEN; /* error: couldn't read frame file! */
    ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
  }
  stream->frame = FrameRead( stream->frfile );
  if ( ! stream->frame )
  {
    stream->err = FRAMESTREAMH_EREAD; /* error: empty frame file! */
    ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
  }
  stream->epoch.gpsSeconds     = stream->frame->GTimeS;
  stream->epoch.gpsNanoSeconds = stream->frame->GTimeN;

  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrSeek(
    LALStatus   *status,
    LIGOTimeGPS *epoch,
    FrStream    *stream
    )
{ /* </lalVerbatim> */
  INT8 twant;
  INT8 tstart;
  INITSTATUS( status, "LALFrSeek", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( epoch, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  twant  = EPOCH_TO_I8TIME( *epoch );
  tstart = SECNAN_TO_I8TIME( stream->frame->GTimeS, stream->frame->GTimeN );

  if ( twant < tstart )
  {
    TRY( LALFrRewind( status->statusPtr, stream ), status );
    tstart = SECNAN_TO_I8TIME( stream->frame->GTimeS, stream->frame->GTimeN );
    if ( twant < tstart )
    {
      LALWarning( status, "Requested time before first frame." );
      stream->gap = 1; /* requested time is just too early */
    }
    else
    {
      TRY( LALFrSeek( status->statusPtr, epoch, stream ), status );
      DETATCHSTATUSPTR( status );
      RETURN( status );
    }
  }
  else
  {
    INT8 tstop = tstart + (INT8)floor( 1e9 * stream->frame->dt );
    while ( twant > tstop )
    {
      TRY( LALFrNext( status->statusPtr, stream ), status );
      if ( stream->gap ) /* location may be in gap */
      {
        tstart = SECNAN_TO_I8TIME( stream->frame->GTimeS, stream->frame->GTimeN );
        if ( twant > tstart ) /* location not in gap: clear gap code */
        {
          stream->gap = 0;
        }
        else
        {
          LALWarning( status, "Requested time in a gap in the frame data." );
        }
      }
      tstop  = SECNAN_TO_I8TIME( stream->frame->GTimeS, stream->frame->GTimeN );
      tstop += (INT8)floor( 1e9 * stream->frame->dt );
    }
  }

  if ( ! stream->gap ) /* set current time to requested time */
  {
    stream->epoch = *epoch;
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrTell(
    LALStatus   *status,
    LIGOTimeGPS *epoch,
    FrStream    *stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALFrTell", FRAMESTREAMC );  
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( epoch, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  *epoch = stream->epoch;
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrGetPos(
    LALStatus *status,
    FrPos     *position,
    FrStream  *stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALFrGetPos", FRAMESTREAMC );  
  ASSERT( position, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  position->epoch   = stream->epoch;
  position->filenum = stream->filenum;
  position->frnum   = stream->frnum;
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrSetPos(
    LALStatus *status,
    FrPos     *position,
    FrStream  *stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALFrSetPos", FRAMESTREAMC );  
  ASSERT( position, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  stream->epoch = position->epoch;
  if ( stream->filenum == position->filenum
      && stream->frnum == position->frnum )
  {
    RETURN( status );
  }

  if ( stream->frame )
  {
    FrameFree( stream->frame );
    stream->frame = NULL;
  }
  if ( stream->frfile )
  {
    FrFileOEnd( stream->frfile );
    stream->frfile = NULL;
  }
  stream->filenum = position->filenum;
  stream->frfile = FrFileINew( stream->filelist[stream->filenum] );
  if ( ! stream->frfile )
  {
    stream->err = FRAMESTREAMH_EOPEN;
    ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
  }

  for ( stream->frnum = 0; stream->frnum <= position->frnum; ++stream->frnum )
  {
    stream->frame = FrameRead( stream->frfile );
    if ( ! stream->frame )
    {
      stream->err = FRAMESTREAMH_EREAD;
      ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
    }
  }
  stream->frnum = position->frnum;

  stream->end = 0;
  stream->err = 0;
  stream->gap = 0;
  RETURN( status );
}



/*
 *
 * The following routines are designed to read and write time series data.
 *
 */


/*
 *
 * Routine to load a FrVect associated with a given channel from a given frame.
 *
 */
static struct FrVect *loadFrVect( struct FrameH *frame, char *name, int chtype )
{
  struct FrVect *vect = NULL;
  if ( chtype == ProcDataChannel )
  {
    struct FrProcData *proc = FrProcDataFind( frame, name );
    vect = proc ? proc->data : NULL;
  }
  else if ( chtype == ADCDataChannel )
  {
    struct FrAdcData *adc = FrAdcDataFind( frame, name );
    vect = adc ? adc->data : NULL;
  }
  else if ( chtype == SimDataChannel )
  {
    struct FrSimData *sim = FrSimDataFind( frame, name );
    vect = sim ? sim->data : NULL;
  }
  return vect;
}


/*
 *
 * Routine to make a 1D FrVect structure and attach it to a frame as the
 * appropriate channel type.
 *
 */
static struct FrVect *makeFrVect1D( struct FrameH *frame, int chtype,
    char *name, char *comment, char *unitx, char *unity, int datatype,
    double rate, double fshift, double dx, unsigned int npts )
{
  FrVect *vect = FrVectNew1D( name, datatype, npts, dx, unitx, unity );
  if ( ! vect ) return NULL;
  if ( chtype == ProcDataChannel )
  {
    struct FrProcData *proc = calloc( 1, sizeof( *proc ) );
    if ( ! proc )
    {
      FrVectFree( vect );
      return NULL;
    }
    proc->classe = FrProcDataDef();
    proc->sampleRate = rate;
    proc->fShift = fshift;
    proc->data = vect;
    proc->next = frame->procData;
    frame->procData = proc;
    if ( ! FrStrCpy( &proc->name, name ) ) return NULL;
    if ( ! FrStrCpy( &proc->comment, comment ) ) return NULL;
  }
  else if ( chtype == ADCDataChannel )
  {
    struct FrAdcData *adc = calloc( 1, sizeof( *adc ) );
    struct FrRawData *raw = frame->rawData ? frame->rawData :
      FrRawDataNew( frame );
    if ( ! adc || ! raw )
    {
      FrVectFree( vect );
      adc ? free( adc ) : 0;
      return NULL;
    }
    adc->classe = FrAdcDataDef();
    adc->sampleRate = rate;
    adc->fShift = fshift;
    adc->data = vect;
    adc->next = raw->firstAdc;
    raw->firstAdc = adc;
    if ( ! FrStrCpy( &adc->name, name ) ) return NULL;
    if ( ! FrStrCpy( &adc->comment, comment ) ) return NULL;
    if ( ! FrStrCpy( &adc->units, unity ) ) return NULL;
  }
  else if ( chtype == SimDataChannel )
  {
    struct FrSimData *sim = calloc( 1, sizeof( *sim ) );
    if ( ! sim )
    {
      FrVectFree( vect );
      return NULL;
    }
    sim->classe = FrSimDataDef();
    sim->sampleRate = rate;
    sim->data = vect;
    sim->next = frame->simData;
    frame->simData = sim;
    if ( ! FrStrCpy( &sim->name, name ) ) return NULL;
    if ( ! FrStrCpy( &sim->comment, comment ) ) return NULL;
  }
  else
  {
    FrVectFree( vect );
    return NULL;
  }
  return vect;
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrGetINT2TimeSeries(
    LALStatus      *status,
    INT2TimeSeries *series,
    FrChanIn       *chanin,
    FrStream       *stream
    )
{ /* </lalVerbatim> */
  struct FrVect *vect;
  UINT4  need;
  UINT4  noff;
  UINT4  ncpy;
  INT2  *dest;
  INT8   tnow;
  INT8   tbeg;

  INITSTATUS( status, "LALFrGetINT2ADC", FRAMESTREAMC );  
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
  if ( vect->type != FR_VECT_2S )
  {
    ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );
  }

  tnow = EPOCH_TO_I8TIME( stream->epoch );
  tbeg = SECNAN_TO_I8TIME( vect->GTimeS, vect->GTimeN );
  if ( tnow < tbeg )
  {
    ABORT( status, FRAMESTREAMH_ETIME, FRAMESTREAMH_MSGETIME );
  }

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
  noff = floor( 1e-9 * ( tnow - tbeg )
      * ( vect->dx[0] ? 1.0 / vect->dx[0] : 0.0 ) );
  if ( noff > vect->nData )
  {
    ABORT( status, FRAMESTREAMH_ETIME, FRAMESTREAMH_MSGETIME );
  }

  /* number of points to copy */
  ncpy = ( vect->nData - noff < need ) ? vect->nData - noff : need;
  memcpy( dest, vect->dataS + noff, ncpy * sizeof( *series->data->data ) );
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
    if ( vect->type != FR_VECT_2S )
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
    memcpy( dest, vect->dataS, ncpy * sizeof( *series->data->data ) );
    dest += ncpy;
    need -= ncpy;
  }

  /* update stream start time */
  SET_EPOCH( &stream->epoch, EPOCH_TO_I8TIME( series->epoch )
      + (INT8)floor( 1e9 * series->data->length * series->deltaT + 0.5 ) );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrGetREAL4TimeSeries(
    LALStatus       *status,
    REAL4TimeSeries *series,
    FrChanIn        *chanin,
    FrStream        *stream
    )
{ /* </lalVerbatim> */
  struct FrVect *vect;
  UINT4  need;
  UINT4  noff;
  UINT4  ncpy;
  REAL4 *dest;
  INT8   tnow;
  INT8   tbeg;

  INITSTATUS( status, "LALFrGetREAL4ADC", FRAMESTREAMC );  
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
  if ( vect->type != FR_VECT_4R )
  {
    ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );
  }

  tnow = EPOCH_TO_I8TIME( stream->epoch );
  tbeg = SECNAN_TO_I8TIME( vect->GTimeS, vect->GTimeN );
  if ( tnow < tbeg )
  {
    ABORT( status, FRAMESTREAMH_ETIME, FRAMESTREAMH_MSGETIME );
  }

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
  noff = floor( 1e-9 * ( tnow - tbeg )
      * ( vect->dx[0] ? 1.0 / vect->dx[0] : 0.0 ) );
  if ( noff > vect->nData )
  {
    ABORT( status, FRAMESTREAMH_ETIME, FRAMESTREAMH_MSGETIME );
  }

  /* number of points to copy */
  ncpy = ( vect->nData - noff < need ) ? vect->nData - noff : need;
  memcpy( dest, vect->dataF + noff, ncpy * sizeof( *series->data->data ) );
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
    if ( vect->type != FR_VECT_4R )
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
    memcpy( dest, vect->dataF, ncpy * sizeof( *series->data->data ) );
    dest += ncpy;
    need -= ncpy;
  }

  /* update stream start time */
  SET_EPOCH( &stream->epoch, EPOCH_TO_I8TIME( series->epoch )
      + (INT8)floor( 1e9 * series->data->length * series->deltaT + 0.5 ) );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrWriteREAL4TimeSeries(
    LALStatus       *status,
    REAL4TimeSeries *series,
    FrOutPar        *params
    )
{ /* </lalVerbatim> */
  REAL4 duration;
  REAL4 *data;
  CHAR seconds[] = "s";
  CHAR comment[] = "Created by LALFrWriteREAL4TimeSeries $Id$";
  CHAR prefix[256];
  CHAR fname[256];
  CHAR units[LALUnitNameSize];
  CHARVector vnits = { LALUnitNameSize, units };
  struct FrFile *frfile;
  UINT4 nframes;
  INT8 t;

  INITSTATUS( status, "LALFrWriteREAL4TimeSeries", FRAMESTREAMC );  
  ASSERT( series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( params, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ATTATCHSTATUSPTR( status );

  strncpy( prefix, params->prefix ? params->prefix : "F", sizeof( prefix ) );
  
  TRY( LALUnitAsString( status->statusPtr, &vnits, &series->sampleUnits ),
      status );

  duration = series->deltaT * series->data->length;
  if ( fabs( duration - 1.0 ) < series->deltaT && params->nframes == 1 )
  {
    /* within one sample of being a one-second frame and one frame per file */
    sprintf( fname, "%s-%u.F", prefix, series->epoch.gpsSeconds );
  }
  else
  {
    sprintf( fname, "%s-%u-%u-%g.F", prefix, series->epoch.gpsSeconds,
        params->nframes, duration / params->nframes );
  }
  frfile = FrFileONew( fname, 0 );

  t = EPOCH_TO_I8TIME( series->epoch );
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
    frame = FrameHNew( prefix );
    frame->run = params->run;
    frame->frame = params->frame++;
    frame->GTimeS = t / (INT8)1000000000;
    frame->GTimeN = t % (INT8)1000000000;
    frame->dt = ncpy * series->deltaT;
    frame->localTime = 0;

    /*** CHECK FSHIFT ***/
    vect = makeFrVect1D( frame, params->type, series->name, comment, seconds,
        units, FR_VECT_4R, series->deltaT ? 1.0 / series->deltaT : 0.0,
        series->f0, series->deltaT, ncpy );
    if ( ! vect )
    {
      ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
    }
    memcpy( vect->dataF, data, ncpy * sizeof( *series->data->data ) );

    FrameWrite( frame, frfile );
    data += ncpy;
    t += 1e9 * ncpy * series->deltaT;
  }

  FrFileOEnd( frfile );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
