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
 * \idx{LALFrOpen}
 * \idx{LALFrCacheOpen}
 * \idx{LALFrClose}
 * \idx{LALFrEnd}
 * \idx{LALFrNext}
 * \idx{LALFrRewind}
 * \idx{LALFrSeek}
 * \idx{LALFrTell}
 * \idx{LALFrGetPos}
 * \idx{LALFrSetPos}
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
 * to \texttt{L-*.gwf}.  If the head name is \texttt{NULL}, the default value
 * \texttt{*.gwf} is used.  The routine \texttt{LALFrCacheOpen()} is like
 * \texttt{LALFrOpen()} except that the list of frame files is taken from a
 * frame file cache.  [In fact, \texttt{LALFrOpen()} simply uses
 * \texttt{LALFrCacheGenerate()} and \texttt{LALFrCacheOpen()} to create the
 * stream.]
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
 *
 * \vfill{\footnotesize\input{FrameStreamCV}}
 *
 **** </lalLaTeX> */

#include <FrameL.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <FrameStreamDef.h>

NRCSID( FRAMESTREAMC, "$Id$" );


/*
 *
 * These functions are for internal use.
 *
 */


/* open a specified frame file URL for input */
static struct FrFile *URLFrFileINew( FrFileInfo *file )
{
  FrFile *frfile = NULL;
  if ( file && file->url )
  {
    if ( ! strncmp( file->url, "file:/", 6 ) )
    {
      frfile = FrFileINew( file->url + 5 );
    }
    else
    {
      if ( lalDebugLevel & LALWARNING )
      {
        LALPrintError( "Warning: function %s, file %s, line %d, %s\n"
            "        %s\n", "URLFrFileINew", __FILE__, __LINE__, FRAMESTREAMC,
            "Unknown URL Type" );
      }
    }
  }
  return frfile;
}

/* compare frame file infos for sorting by time */
static int FileListSortCompare( const void *p1, const void *p2 )
{
  const FrFileInfo *file1 = p1;
  const FrFileInfo *file2 = p2;
  int ans = 0;
  if ( file1->t0 < file2->t0 )
    ans = -1;
  else if ( file1->t0 > file2->t0 )
    ans = 1;
  return ans;
}

/* compare frame file infos for selecting a time */
static int FileListSelectCompare( const void *p1, const void *p2 )
{
  const FrFileInfo *key  = p1;
  const FrFileInfo *file = p2;
  int ans = 0;
  if ( key->t0 < file->t0 )
    ans = -1;
  else if ( key->t0 > file->t0 + file->dt )
    ans = 1;
  return ans;
}

/* destroy a frame file list */
static void
LALDestroyFrFileList(
    LALStatus   *status,
    FrFileInfo **list
    )
{
  FrFileInfo *file;
  INITSTATUS( status, "LALDestroyFrFileList", FRAMESTREAMC );  
  file = *list;
  while ( file->ind >= 0 && file->url )
  {
    LALFree( file->url );
    ++file;
  }
  LALFree( *list );
  *list = NULL;
  RETURN( status );
}

/* create a frame file list from a cache */
static void
LALCreateFrFileList(
    LALStatus   *status,
    FrFileInfo **output,
    FrCache     *cache
    )
{
  FrFileInfo *list;
  UINT4 i;
  INITSTATUS( status, "LALCreateFrFileList", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  list = *output = LALCalloc( cache->numFrameFiles + 1, sizeof( *list ) );
  list[cache->numFrameFiles].ind = -1; /* indicates end */
  for ( i = 0; i < cache->numFrameFiles; ++i )
  {
    FrStat *file = cache->frameFiles + i;
    list[i].url = LALMalloc( strlen( file->url ) + 1 );
    if ( ! list[i].url )
    {
      TRY( LALDestroyFrFileList( status->statusPtr, output ), status );
      ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
    }
    strcpy( list[i].url, file->url );
    if ( file->startTime > 0 && file->duration > 0 )
    {
      list[i].t0 = file->startTime;
      list[i].dt = file->duration;
    }
    else
    {
      struct FrFile *frfile;
      struct FrameH *frame;
      double endTime;
      frfile = URLFrFileINew( list + i );
      if ( ! frfile )
      {
        TRY( LALDestroyFrFileList( status->statusPtr, output ), status );
        ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
      }
      frame = FrameRead( frfile );
      if ( ! frame )
      {
        TRY( LALDestroyFrFileList( status->statusPtr, output ), status );
        ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
      }
      list[i].t0 = frame->GTimeS;
      do
      {
        endTime  = frame->GTimeS + 1e-9 * frame->GTimeN;
        endTime += frame->dt;
        FrameFree( frame );
      }
      while ( ( frame = FrameRead( frfile ) ) );
      list[i].dt = ceil( endTime ) - list[i].t0;
      FrFileIEnd( frfile );
    }
  }

  qsort( list, cache->numFrameFiles, sizeof( *list ), FileListSortCompare );
  for ( i = 0; i < cache->numFrameFiles; ++i )
    list[i].ind = i;

  DETATCHSTATUSPTR( status );
  RETURN( status );
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
LALFrCacheOpen(
    LALStatus  *status,
    FrStream  **output,
    FrCache    *cache
    )
{ /* </lalVerbatim> */
  FrStream *stream;

  INITSTATUS( status, "LALFrCacheOpen", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( cache, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( output, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( ! *output, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL );

  stream = *output = LALCalloc( 1, sizeof( **output ) );
  if ( ! stream )
  {
    ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
  }

  stream->numfiles = cache->numFrameFiles;
  LALCreateFrFileList( status->statusPtr, &stream->filelist, cache );
  BEGINFAIL( status )
  {
    LALFree( *output );
    *output = NULL;
  }
  ENDFAIL( status );

  FrLibIni( NULL, stderr, 0 );

  stream->frfile = URLFrFileINew( stream->filelist );
  if ( ! stream->frfile )
  {
    TRY( LALDestroyFrFileList( status->statusPtr, &stream->filelist ), status );
    LALFree( *output );
    *output = NULL;
    ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
  }

  stream->frame = FrameRead( stream->frfile );
  if ( ! stream->frame )
  {
    FrFileIEnd( stream->frfile );
    TRY( LALDestroyFrFileList( status->statusPtr, &stream->filelist ), status );
    LALFree( *output );
    *output = NULL;
    ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
  }

  stream->epoch.gpsSeconds     = stream->frame->GTimeS;
  stream->epoch.gpsNanoSeconds = stream->frame->GTimeN;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrOpen(
    LALStatus    *status,
    FrStream    **stream,
    const CHAR   *dirname,
    const CHAR   *pattern
    )
{ /* </lalVerbatim> */
  FrCache *cache = NULL;

  INITSTATUS( status, "LALFrOpen", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( ! *stream, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL );

  LALFrCacheGenerate( status->statusPtr, &cache, dirname, pattern );
  CHECKSTATUSPTR( status );

  LALFrCacheOpen( status->statusPtr, stream, cache );
  BEGINFAIL( status )
  {
    TRY( LALDestroyFrCache( status->statusPtr, &cache ), status );
  }
  ENDFAIL( status );

  LALDestroyFrCache( status->statusPtr, &cache );
  BEGINFAIL( status )
  {
    TRY( LALFrClose( status->statusPtr, stream ), status );
  }
  ENDFAIL( status );

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
  INITSTATUS( status, "LALFrClose", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( *stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  FrameFree( (*stream)->frame );
  FrFileIEnd( (*stream)->frfile );
  TRY( LALDestroyFrFileList( status->statusPtr, &(*stream)->filelist ), status );
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

  if ( ! stream->frame ) /* no more frames in file */
  {
    FrFileIEnd( stream->frfile );
    stream->frfile = NULL;
    ++stream->filenum;

    /* no more frame files: close stream */
    if ( stream->filenum >= stream->numfiles )
    {
      stream->end = EOF;
      RETURN( status );
    }

    stream->frfile = URLFrFileINew( stream->filelist + stream->filenum );
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
    FrFileIEnd( stream->frfile );
  }
  stream->filenum = 0;
  stream->end     = 0;
  stream->err     = 0;
  stream->gap     = 0;
  stream->frfile  = URLFrFileINew( stream->filelist );
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
  FrFileInfo  key;
  FrFileInfo *file;
  INT4 overlap = 0;

  INITSTATUS( status, "LALFrSeek", FRAMESTREAMC );  
  ATTATCHSTATUSPTR( status );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( epoch, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( stream->frame )
  {
    FrameFree( stream->frame );
    stream->frame = NULL;
  }
  if ( stream->frfile )
  {
    FrFileIEnd( stream->frfile );
    stream->frfile = NULL;
  }
  stream->end = 0;
  stream->err = 0;
  stream->gap = 0;

  file = stream->filelist;
  if ( epoch->gpsSeconds < file->t0 )
  {
    TRY( LALFrRewind( status->statusPtr, stream ), status );
    LALWarning( status, "Requested time before first frame." );
    stream->gap = 1; /* requested time is just too early */
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

  file = stream->filelist + stream->numfiles - 1;
  if ( epoch->gpsSeconds > file->t0 + file->dt )
  {
    LALWarning( status, "Requested time after last frame." );
    stream->filenum = stream->numfiles;
    stream->end = EOF;
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

  key.t0 = epoch->gpsSeconds;
  file = bsearch( &key, stream->filelist, stream->numfiles,
      sizeof( *stream->filelist ), FileListSelectCompare );
  if ( file ) /* check for possible overlap */
  {
    FrFileInfo *file2;
    file2 = stream->filelist + file->ind - 1;
    if ( file->ind > 0 && epoch->gpsSeconds < file2->t0 + file2->dt )
    { /* choose earlier file */
      overlap = 1;
      file = file2;
    }
    file2 = stream->filelist + file->ind + 1;
    if ( file->ind + 2 < (int)stream->numfiles && file2->t0 < epoch->gpsSeconds )
    {
      overlap = 1;
    }
  }
  else /* no match found: must be in a gap */
  {
    UINT4 i;
    stream->gap = 1;
    LALWarning( status, "Requested time in gap in frame data." );
    /* scan frame list in order to find earliest file after time requested */
    for ( i = 0; i < stream->numfiles; ++i )
    {
      file = stream->filelist + i;
      if ( file->t0 > epoch->gpsSeconds )
        break;
    }
  }

  stream->filenum = file->ind;
  stream->frfile  = URLFrFileINew( stream->filelist + stream->filenum );
  if ( ! stream->frfile )
  {
    stream->err = FRAMESTREAMH_EOPEN;
    ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
  }

  if ( stream->gap ) /* just open first frame in file */
  {
    stream->frame = FrameRead( stream->frfile );
  }
  else
  {
    double twant = epoch->gpsSeconds + 1e-9 * epoch->gpsNanoSeconds;
    stream->frame = FrameReadT( stream->frfile, twant );
    if ( overlap && ! stream->frame ) /* could be in the other frame file */
    {
      ++stream->filenum;
      stream->frfile = URLFrFileINew( stream->filelist + stream->filenum );
      if ( ! stream->frfile )
      {
        stream->err = FRAMESTREAMH_EOPEN;
        ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
      }
      stream->frame = FrameReadT( stream->frfile, twant );
    }
    /* FIXME: probably can't handle gaps within a frame file */
  }
  if ( ! stream->frame )
  {
    stream->err = FRAMESTREAMH_EREAD;
    ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
  }

  if ( stream->gap )
  {
    stream->epoch.gpsSeconds     = stream->frame->GTimeS;
    stream->epoch.gpsNanoSeconds = stream->frame->GTimeN;
  }
  else
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
  position->frame   = stream->frame->frame;
  position->run     = stream->frame->run;
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
      && stream->frame->frame == position->frame
      && stream->frame->run == position->run )
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
    FrFileIEnd( stream->frfile );
    stream->frfile = NULL;
  }
  stream->filenum = position->filenum;
  stream->frfile  = URLFrFileINew( stream->filelist + stream->filenum );
  if ( ! stream->frfile )
  {
    stream->err = FRAMESTREAMH_EOPEN;
    ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
  }

  stream->frame = FrameReadN( stream->frfile, position->run, position->frame );
  if ( ! stream->frame )
  {
    stream->err = FRAMESTREAMH_EREAD;
    ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
  }

  stream->end = 0;
  stream->err = 0;
  stream->gap = 0;
  RETURN( status );
}
