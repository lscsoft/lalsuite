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
 * \idx{LALOpenFrameStream()}
 * \idx{LALCloseFrameStream()}
 * \idx{LALSFrameReadADCTimeSeries()}
 * \idx{LALI2FrameReadADCTimeSeries()}
 * \idx{LALNextFrame()}
 * \idx{LALFrameStreamError()}
 * \idx{LALFrameStreamGetPos()}
 * \idx{LALFrameStreamSetPos()}
 *
 * \subsubsection*{Description}
 *
 * The routines \texttt{LALOpenFrameStream()} and
 * \texttt{LALCloseFrameStream()} are used to open and close a frame stream.
 * The stream is created by \texttt{LALOpenFrameStream()}, and must be a
 * pointer to \texttt{NULL} before it is opened.  It must have been created
 * prior to calling \texttt{LALCloseFrameStream()}, and after this call, the
 * stream will be a pointer to \texttt{NULL}.  The routine
 * \texttt{LALOpenFrameStream()} requires the user to specify the directory
 * name of the frame files and the head names.  If the directory is
 * \texttt{NULL}, the routine uses the current director (\texttt{.}).  The
 * head names specifies which files are the wanted files in the specified
 * directory.  Wildcards are allowed.  For example, to get LLO frames only,
 * the head names could be set to \texttt{L-*.F}.  If the head name is
 * \texttt{NULL}, the default value \texttt{*.F} is used.
 *
 * The routines
 * \texttt{LAL}$\langle\mbox{datatype}\rangle$\texttt{FrameReadADCTimeSeries()}
 * extract the specified ADC channel (the user must input the channel name)
 * and creates a time series (of the specified datatype) to store the data.  If
 * the datatype of the calling routine differs from the datatype of the stored
 * data, the routine fails --- it does not try to convert the data.  The data
 * stream is \emph{not} advanced by this routine (so that several channels
 * can be extracted from the same frame).
 *
 * The routine \texttt{LALNextFrame()} advances the frame stream to the next
 * frame (possibly advancing to a new frame file too).  An error code is
 * returned (as output) containing the state of stream.  If the error code
 * is zero, no errors have been detected.  If the error code is \texttt{EOF}
 * then an error has occurred, usually because the end of frames has been
 * reached.  In this case, any subsequent call to \texttt{LALNextFrame()} or
 * one of the frame reading routines will fail.  Also, if a gap in time is
 * detected between the new frame and the previous frame, the error code is
 * set to 1, but this is not stored in the state of the frame stream and can
 * not be determined using \texttt{LALFrameStreamError()} (below).
 *
 * The routine \texttt{LALFrameStreamError()} return the value of the frame
 * stream error status, which is 0 if there is no error or \texttt{EOF} if
 * an error has occured (usually meaning that the end of the frame files
 * has been reached).
 *
 * The routines \texttt{LALFrameStreamGetPos()} and
 * \texttt{LALFrameStreamSetPos()} are used to get a record of the frame stream
 * state and to restore the frame stream to the recorded state respectively.
 * This is useful if you want to repeatedly read several frames: you can save
 * the starting position, read several frames, and then restore the stream to
 * this starting position.
 *
 * \vfill{\footnotesize\input{FrameStreamCV}}
 *
 **** </lalLaTeX> */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <FrameL.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/FrameStream.h>

NRCSID( FRAMESTREAMC, "$Id$" );

typedef struct
tagFrameFiles
{
  UINT4          nfiles;
  CHAR         **fnames;
}
FrameFiles;

struct
tagFrameStream
{
  FrameFiles    *filelist;
  UINT4          filenum;
  struct FrFile *frfile;
  UINT4          frnum;
  struct FrameH *frame;
  INT4           error;
};

static void
GetFrameFiles(
    LALStatus   *status,
    FrameFiles **frfiles,
    const CHAR  *dirname,
    const CHAR  *headname
    )
{
  char cmd[1024];
  char tmp[1024];
  UINT4 file;
  FILE *fp;

  INITSTATUS( status, "LALGetFrameFiles", FRAMESTREAMC );  
  ASSERT( frfiles, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( ! *frfiles, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL );

  *frfiles = LALCalloc( 1, sizeof( **frfiles ) );
  if ( ! *frfiles )
  {
    ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
  }

  LALSnprintf( cmd, sizeof( cmd ), "ls %s/%s 2>/dev/null",
      dirname ? dirname : ".", headname ? headname : "*.F" );

  /*
   *
   * Count number of frame files
   *
   */

  (*frfiles)->nfiles = 0;
  fp = popen( cmd, "r" );
  if ( ! fp )
  {
    ABORT( status, FRAMESTREAMH_EPIPE, FRAMESTREAMH_MSGEPIPE );
  }
  while ( fgets( tmp, sizeof( tmp ), fp ) )
    ++(*frfiles)->nfiles;
  pclose( fp );
  if ( ! (*frfiles)->nfiles )
  {
    ABORT( status, FRAMESTREAMH_EFILE, FRAMESTREAMH_MSGEFILE );
  }

  /* null-terminated list of file names */
  (*frfiles)->fnames =
    LALCalloc( (*frfiles)->nfiles + 1, sizeof( *(*frfiles)->fnames ) );
  if ( ! (*frfiles)->fnames )
  {
    ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
  }


  /*
   *
   * Get the file names
   *
   */

  file = 0;
  fp = popen( cmd, "r" );
  if ( ! fp )
  {
    ABORT( status, FRAMESTREAMH_EPIPE, FRAMESTREAMH_MSGEPIPE );
  }
  while ( fgets( tmp, sizeof( tmp ), fp ) )
  {
    (*frfiles)->fnames[file] = LALMalloc( strlen( tmp ) );
    if ( ! (*frfiles)->fnames[file] )
    {
      ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
    }
    sscanf( tmp, "%s\n", (*frfiles)->fnames[file] );
    ++file;
  }
  pclose( fp );

  RETURN( status );
}


static void
DestroyFrameFiles(
    LALStatus   *status,
    FrameFiles **frfiles
    )
{
  UINT4 file;
  INITSTATUS( status, "LALDestroyFrameFiles", FRAMESTREAMC );  
  ASSERT( frfiles, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( *frfiles, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  for ( file = 0; file < (*frfiles)->nfiles; ++file )
  {
    LALFree( (*frfiles)->fnames[file] );
  }
  LALFree( (*frfiles)->fnames );
  LALFree( *frfiles );
  *frfiles = NULL;

  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALOpenFrameStream(
    LALStatus    *status,
    FrameStream **stream,
    const CHAR   *dirname,
    const CHAR   *headname
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALOpenFrameStream", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( ! *stream, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL );

  FrLibIni( NULL, stderr, 0 );

  *stream = LALCalloc( 1, sizeof( **stream ) );
  if ( ! *stream )
  {
    ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
  }

  GetFrameFiles( status->statusPtr, &(*stream)->filelist, dirname, headname );
  BEGINFAIL( status )
  {
    LALFree( *stream );
    *stream = NULL;
  }
  ENDFAIL( status );

  (*stream)->frfile = FrFileINew( (*stream)->filelist->fnames[0] );
  if ( ! (*stream)->frfile )
  {
    TRY( DestroyFrameFiles( status->statusPtr, &(*stream)->filelist ), status );
    LALFree( *stream );
    *stream = NULL;
    ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
  }
  (*stream)->frame = FrameRead( (*stream)->frfile );
  if ( ! (*stream)->frame )
  {
    FrFileOEnd( (*stream)->frfile );
    TRY( DestroyFrameFiles( status->statusPtr, &(*stream)->filelist ), status );
    LALFree( *stream );
    *stream = NULL;
    ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALCloseFrameStream(
    LALStatus    *status,
    FrameStream **stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALCloseFrameStream", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( *stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( (*stream)->frame )
  {
    FrameFree( (*stream)->frame );
  }

  if ( (*stream)->frfile )
  {
    FrFileOEnd( (*stream)->frfile );
  }

  TRY( DestroyFrameFiles( status->statusPtr, &(*stream)->filelist ), status );

  LALFree( *stream );
  *stream = NULL;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALSFrameReadADCTimeSeries(
    LALStatus        *status,
    REAL4TimeSeries **series,
    const CHAR       *channel,
    FrameStream      *stream
    )
{ /* </lalVerbatim> */
  char chan[128];
  struct FrAdcData *adc;
  size_t size;

  INITSTATUS( status, "LALSFrameReadADCTimeSeries", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( channel, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( ! *series, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL );

  if ( stream->error == EOF )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }

  strncpy( chan, channel, sizeof( chan ) );
  adc = FrAdcDataFind( stream->frame, chan );
  if ( ! adc )
  {
    ABORT( status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN );
  }
  if ( adc->data->type != FR_VECT_4R )
  {
    ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );
  }

  *series = LALCalloc( 1, sizeof( **series ) );
  if ( ! *series )
  {
    ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
  }

  strncpy( (*series)->name, channel, sizeof( (*series)->name ) );
  (*series)->epoch.gpsSeconds = stream->frame->GTimeS;
  (*series)->epoch.gpsNanoSeconds = stream->frame->GTimeN;
  (*series)->deltaT = ( adc->sampleRate == 0.0 ? 0.0 : 1.0 / adc->sampleRate );
  (*series)->f0 = adc->fShift; /* IS THIS CORRECT??? */
  
  LALSCreateVector( status->statusPtr, &(*series)->data, adc->data->nData );
  BEGINFAIL( status )
  {
    LALFree( *series );
    *series = NULL;
  }
  ENDFAIL( status );
  size = adc->data->nData * sizeof( *(*series)->data->data );
  memcpy( (*series)->data->data, adc->data->dataF, size );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALI2FrameReadADCTimeSeries(
    LALStatus       *status,
    INT2TimeSeries **series,
    const CHAR      *channel,
    FrameStream     *stream
    )
{ /* </lalVerbatim> */
  char chan[128];
  struct FrAdcData *adc;
  size_t size;

  INITSTATUS( status, "LALI2FrameReadADCTimeSeries", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( channel, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( series, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( ! *series, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL );

  if ( stream->error == EOF )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }

  strncpy( chan, channel, sizeof( chan ) );
  adc = FrAdcDataFind( stream->frame, chan );
  if ( ! adc )
  {
    ABORT( status, FRAMESTREAMH_ECHAN, FRAMESTREAMH_MSGECHAN );
  }
  if ( adc->data->type != FR_VECT_2S )
  {
    ABORT( status, FRAMESTREAMH_ETYPE, FRAMESTREAMH_MSGETYPE );
  }

  *series = LALCalloc( 1, sizeof( **series ) );
  if ( ! *series )
  {
    ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
  }

  strncpy( (*series)->name, channel, sizeof( (*series)->name ) );
  (*series)->epoch.gpsSeconds = stream->frame->GTimeS;
  (*series)->epoch.gpsNanoSeconds = stream->frame->GTimeN;
  (*series)->deltaT = ( adc->sampleRate == 0.0 ? 0.0 : 1.0 / adc->sampleRate );
  (*series)->f0 = adc->fShift; /* IS THIS CORRECT??? */
  
  LALI2CreateVector( status->statusPtr, &(*series)->data, adc->data->nData );
  BEGINFAIL( status )
  {
    LALFree( *series );
    *series = NULL;
  }
  ENDFAIL( status );
  size = adc->data->nData * sizeof( *(*series)->data->data );
  memcpy( (*series)->data->data, adc->data->dataS, size );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALNextFrame(
    LALStatus   *status,
    INT4        *error,
    FrameStream *stream
    )
{ /* </lalVerbatim> */
  /* timing accuracy: tenth of a sample interval for a 16kHz fast channel */
  const INT8 tacc = (INT8)floor( 0.1 * 1e9 / 16384.0 );
  INT8 texp;
  INT8 tact;
  INITSTATUS( status, "LALNextFrame", FRAMESTREAMC );  
  ASSERT( error, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  *error = stream->error;
  if ( stream->error == EOF )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }

  /* compute expected start time of the next frame */
  texp  = (INT8)stream->frame->GTimeS * (INT8)1000000000;
  texp += (INT8)stream->frame->GTimeN;
  texp += (INT8)floor( 1e9 * stream->frame->dt );

  if ( stream->frame )
  {
    FrameFree( stream->frame );
    stream->frame = NULL;
  }

  stream->frame = FrameRead( stream->frfile );
  ++stream->frnum;

  if ( ! stream->frame ) /* no more frames in file */
  {
    FrFileOEnd( stream->frfile );
    stream->frfile = NULL;
    stream->frnum = 0;

    /* no more frame files: close stream */
    if ( ++stream->filenum >= stream->filelist->nfiles )
    {
      stream->error = *error = EOF;
      RETURN( status );
    }

    stream->frfile = FrFileINew( stream->filelist->fnames[stream->filenum] );
    if ( ! stream->frfile )
    {
      stream->error = *error = EOF;
      ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
    }

    stream->frame = FrameRead( stream->frfile );
    if ( ! stream->frame )
    {
      stream->error = *error = EOF;
      ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
    }
  }

  /* compute actual start time of this new frame */
  tact  = (INT8)stream->frame->GTimeS * (INT8)1000000000;
  tact += (INT8)stream->frame->GTimeN;

  if ( abs( texp - tact ) > tacc ) /* there is a gap */
  {
    *error = 1;
  }

  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrameStreamError(
    LALStatus   *status,
    INT4        *error,
    FrameStream *stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALFrameStreamError", FRAMESTREAMC );  
  ASSERT( error, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  *error = stream->error;
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrameStreamGetPos(
    LALStatus      *status,
    FrameStreamPos *position,
    FrameStream    *stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALFrameStreamGetPos", FRAMESTREAMC );  
  ASSERT( position, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  position->gpstime.gpsSeconds = stream->frame->GTimeS;
  position->gpstime.gpsNanoSeconds = stream->frame->GTimeN;
  position->filenum = stream->filenum;
  position->frnum = stream->frnum;
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrameStreamSetPos(
    LALStatus      *status,
    FrameStreamPos *position,
    FrameStream    *stream
    )
{ /* </lalVerbatim> */
  UINT4 fr;
  INITSTATUS( status, "LALFrameStreamSetPos", FRAMESTREAMC );  
  ASSERT( position, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  stream->filenum = position->filenum;
  stream->frnum = position->frnum;
  stream->frame = NULL;
  stream->frfile = FrFileINew( stream->filelist->fnames[stream->filenum] );
  if ( ! stream->frfile )
  {
    stream->error = EOF;
    ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
  }
  for ( fr = 0; fr <= stream->frnum; ++fr )
  {
    stream->frame = FrameRead( stream->frfile );
    if ( ! stream->frame )
    {
      stream->error = EOF;
      ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
    }
  }
  stream->error = 0;
  RETURN( status );
}
