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

NRCSID( FRAMESTREAMC, "$Id$" );


/*
 *
 * These functions are for internal use.
 *
 */



#include <errno.h>
#include <unistd.h>
#include <sys/param.h>
#define STR( x ) #x
#define XSTR( x ) STR( x )
#define MAXPROTOCOLLEN 16
static struct FrFile *URLFrFileINew( FrFileInfo *file )
{
  struct FrFile *frfile = NULL;
  char prot[MAXPROTOCOLLEN + 1] = "";
  char host[MAXHOSTNAMELEN + 1] = "";
  char path[MAXPATHLEN + 1]     = "";
  int n;

  /* get protocol, hostname, and path */
  if ( ! file || ! file->url )
  {
    errno = EFAULT;
    return NULL;
  }
  n = sscanf( file->url, "%" XSTR( MAXPROTOCOLLEN ) "[^:]://%"
      XSTR( MAXHOSTNAMELEN ) "[^/]%" XSTR( MAXPATHLEN ) "s", prot, host, path );
  if ( n != 3 ) /* perhaps the hostname has been omitted */
  {
    n = sscanf( file->url, "%" XSTR( MAXPROTOCOLLEN ) "[^:]://%"
        XSTR( MAXPATHLEN ) "s", prot, path );
    if ( n == 2 )
      strcpy( host, "localhost" );
    else
    {
      strncpy( path, file->url, MAXPATHLEN );
      strcpy( prot, "none" );
    }
  }

  /* process various protocols */
  if ( ! strcmp( prot, "none" ) )
  { /* assume a file on the localhost */
    /* TODO: should check for leading ~ */
    frfile = FrFileINew( path );
  }
  else if ( ! strcmp( prot, "file" ) )
  {
    if ( strcmp( host, "localhost" ) )
    { /* make sure the host *is* localhost */
      char localhost[MAXHOSTNAMELEN + 1];
      gethostname( localhost, MAXHOSTNAMELEN );
      if ( strcmp( host, localhost ) )
      { /* nope */
        fprintf( stderr, "Can not read files from remote hosts.\n" );
        errno = EINVAL;
        return NULL;
      }
    }
    frfile = FrFileINew( path );
  }
  else
  {
    fprintf( stderr, "Unsupported protocol %s.", prot );
    errno = EINVAL;
    return NULL;
  }

  return frfile;
}

static void free_flist( FrFileInfo *flist )
{
  FrFileInfo *p = flist;
  if ( ! p )
    return;
  while ( p->ind > -1 )
  {
    if ( p->url )
      LALFree( p->url );
    memset( p, 0, sizeof( *p ) );
    ++p;
  }
  LALFree( flist );
  return;
}

#include <glob.h>
static FrStream *fr_open( const char *pattern )
{
  char host[MAXHOSTNAMELEN + 1] = "localhost";
  glob_t g;
  int i;
  FrStream *stream;

  if ( ! pattern )
    pattern = "*.gwf";
  if ( gethostname( host, MAXHOSTNAMELEN ) )
    return NULL;
  g.gl_offs = 0;
  if ( glob( pattern, GLOB_DOOFFS, NULL, &g ) )
    return NULL;

  stream = LALCalloc( 1, sizeof( *stream ) );
  if ( ! stream )
  {
    errno = ENOMEM;
    return NULL;
  }
  stream->nfile = g.gl_matchc;
  stream->flist = LALCalloc( stream->nfile + 1, sizeof( *stream->flist ) );
  if ( ! stream->flist )
  {
    errno = ENOMEM;
    LALFree( stream );
    return NULL;
  }

  for ( i = 0; i < g.gl_matchc; ++i )
  {
    int t0 = 0;
    int dt = 0;
    char *path = g.gl_pathv[i];
    char *base = strrchr( path, '/' );
    stream->flist[i].ind = i;
    if ( ! base )
      base = path;
    if ( *path == '/' ) /* absolute path */
    {
      stream->flist[i].url = LALMalloc( sizeof( "file://" )
          + strlen( host ) + strlen( path ) );
      if ( ! stream->flist[i].url )
      {
        errno = ENOMEM;
        free_flist( stream->flist );
        LALFree( stream );
        return NULL;
      }
      sprintf( stream->flist[i].url, "file://%s%s", host, path );
    }
    else
    {
      char cwd[MAXPATHLEN + 1] = "./";
      getcwd( cwd, MAXPATHLEN );
      stream->flist[i].url = LALMalloc( sizeof( "file://" )
          + strlen( host ) + strlen( cwd ) + strlen( path ) );
      if ( ! stream->flist[i].url )
      {
        errno = ENOMEM;
        free_flist( stream->flist );
        LALFree( stream );
        return NULL;
      }
      sprintf( stream->flist[i].url, "file://%s%s/%s", host, cwd, path );
    }
    /* get t0 and dt from file name */
    if ( 2 == sscanf( base ? base : path, "%*[^-]-%*[^-]-%d-%d.gwf", &t0, &dt )
        && t0 > 0 && dt > 0 )
    {
      stream->flist[i].t0 = t0;
      stream->flist[i].dt = dt;
    }
    else /* have to open the file to determine t0 and dt */
    {
      stream->file = URLFrFileINew( stream->flist + i );
      if ( ! stream->file )
      {
        free_flist( stream->flist );
        LALFree( stream );
        return NULL;
      }
      if ( FrTOCReadFull( stream->file ) == NULL )
      {
        FrFileIEnd( stream->file );
        free_flist( stream->flist );
        LALFree( stream );
        return NULL;
      }
      t0 = floor( stream->file->toc->GTimeS[0] );
      dt = ceil( (double)stream->file->toc->GTimeS[0]
            + 1e-9 * (double)stream->file->toc->GTimeN[0]
            + stream->file->toc->dt[0] ) - t0;
      stream->flist[i].t0 = t0;
      stream->flist[i].dt = dt;
    }
  }
  stream->flist[stream->nfile].ind = -1;
  stream->file = URLFrFileINew( stream->flist );
  if ( ! stream->file )
  {
    free_flist( stream->flist );
    LALFree( stream );
    return NULL;
  }
  if ( FrTOCReadFull( stream->file ) == NULL )
  {
    FrFileIEnd( stream->file );
    free_flist( stream->flist );
    LALFree( stream );
    return NULL;
  }
  stream->epoch.gpsSeconds     = stream->file->toc->GTimeS[0];
  stream->epoch.gpsNanoSeconds = stream->file->toc->GTimeN[0];
  return stream;
}

static int fr_close( FrStream *stream )
{
  if ( ! stream )
  {
    errno = EFAULT;
    return EOF;
  }
  FrFileIEnd( stream->file );
  free_flist( stream->flist );
  memset( stream, 0, sizeof( *stream ) );
  return 0;
}

static int fr_state( FrStream *stream )
{
  if ( ! stream )
  {
    errno = EFAULT;
    return EOF;
  }
  return stream->state;
}

static void fr_clearerr( FrStream *stream )
{
  if ( ! stream )
    errno = EFAULT;
  else
    stream->state = LAL_FR_OK;
  return;
}

static int fr_rewind( FrStream *stream )
{
  if ( ! stream )
  {
    errno = EFAULT;
    return EOF;
  }
  if ( stream->file )
  {
    FrFileIEnd( stream->file );
    stream->file = NULL;
  }
  stream->pos   = 0;
  stream->fnum  = 0;
  stream->file  = URLFrFileINew( stream->flist );
  stream->state = LAL_FR_OK;
  if ( ! stream->file )
  {
    stream->state |= LAL_FR_ERR | LAL_FR_URL;
    return EOF;
  }
  if ( FrTOCReadFull( stream->file ) == NULL )
  {
    FrFileIEnd( stream->file );
    stream->file   = NULL;
    stream->state |= LAL_FR_ERR | LAL_FR_TOC;
    return EOF;
  }
  stream->epoch.gpsSeconds     = stream->file->toc->GTimeS[0];
  stream->epoch.gpsNanoSeconds = stream->file->toc->GTimeN[0];
  return 0;
}

static int fr_next( FrStream *stream )
{
  /* timing accuracy: tenth of a sample interval for a 16kHz fast channel */
  const INT8 tacc = (INT8)floor( 0.1 * 1e9 / 16384.0 );
  INT8 tnow = 0;
  INT8 texp = 0;
  INT8 tact;

  if ( ! stream )
  {
    errno = EFAULT;
    return EOF;
  }

  if ( stream->state & LAL_FR_END )
    return EOF;

  /* turn off gap bit */
  stream->state &= ~LAL_FR_GAP;

  /* FIXME: assume that stream->file is open */
  if ( stream->file && stream->file->toc )
  {
    if ( stream->pos < stream->file->toc->nFrame )
    {
      tnow  = (INT8)1000000000 * (INT8)stream->file->toc->GTimeS[stream->pos];
      tnow += stream->file->toc->GTimeN[stream->pos];
      texp  = tnow;
      texp += (INT8)floor( 1e9 * stream->file->toc->dt[stream->pos] );
      ++stream->pos;
    }
    if ( stream->pos >= stream->file->toc->nFrame )
    {
      FrFileIEnd( stream->file );
      stream->file = NULL;
      stream->pos = 0;
      ++stream->fnum;
    }
  }

  /* open a new file if necessary */
  if ( ! stream->file )
  {
    stream->pos = 0;
    if ( stream->fnum >= stream->nfile )
    {
      stream->state |= LAL_FR_END;
      return EOF;
    }
    stream->pos  = 0;
    stream->file = URLFrFileINew( stream->flist + stream->fnum );
    if ( ! stream->file )
    {
      stream->state |= LAL_FR_ERR | LAL_FR_URL;
      return EOF;
    }
  }

  /* open TOC if necessary */
  if ( ! stream->file->toc )
  {
    if ( FrTOCReadFull( stream->file ) == NULL )
    {
      FrFileIEnd( stream->file );
      stream->file   = NULL;
      stream->state |= LAL_FR_ERR | LAL_FR_TOC;
      return EOF;
    }
  }

  /* compute actual start time of this new frame */
  tact  = (INT8)1000000000 * (INT8)stream->file->toc->GTimeS[stream->pos];
  tact += (INT8)stream->file->toc->GTimeN[stream->pos];
  stream->epoch.gpsSeconds     = stream->file->toc->GTimeS[stream->pos];
  stream->epoch.gpsNanoSeconds = stream->file->toc->GTimeN[stream->pos];

  if ( abs( texp - tact ) > tacc ) /* there is a gap */
    stream->state |= LAL_FR_GAP;

  return 0;
}

/* compare routine for binary search bsearch to locate the first file
 * containing wanted time, or the first file after a gap if the wanted time
 * occurs during a gap */
static int flist_tcompare( const void *key, const void *ptr )
{
  const FrFileInfo *file = ptr; /* this file */
  const FrFileInfo *prev = file - 1; /* the previous file */
  double twant = *((double*)key);
  if ( twant < file->t0 )
  {
    /* check if previous file was before wanted time */
    if ( file->ind > 0 && twant > prev->t0 + prev->dt )
      return 0; /* gap during wanted time and this is first file after gap */
    return -1; /* this file is after wanted time */
  }
  if ( twant > file->t0 + file->dt )
    return 1; /* this file is before wanted time */
  if ( file->ind > 0 && twant < prev->t0 + prev->dt )
    return -1; /* the previous file contains wanted time too */
  return 0;
}

static int fr_seek( FrStream *stream, LIGOTimeGPS *epoch )
{
  FrFileInfo *file1;
  FrFileInfo *file2;
  double twant;
  UINT4 i;

  if ( ! stream || ! epoch )
  {
    errno = EFAULT;
    return EOF;
  }

  if ( stream->file )
  {
    FrFileIEnd( stream->file );
    stream->file = NULL;
    stream->pos = 0;
  }

  /* clear EOF or GAP states; preserve ERR state */
  if ( stream->state & LAL_FR_ERR )
    stream->state = LAL_FR_ERR;
  else
    stream->state = LAL_FR_OK;

  /* if epoch is before first file */
  file1 = stream->flist;
  if ( epoch->gpsSeconds < file1->t0 )
  {
    errno = EINVAL;
    return fr_rewind( stream );
    stream->state |= LAL_FR_GAP;
  }

  /* if epoch is after last file */
  file2 = stream->flist + stream->nfile - 1;
  if ( epoch->gpsSeconds > file2->t0 + file2->dt )
  {
    errno = EINVAL;
    stream->fnum   = stream->nfile;
    stream->state |= LAL_FR_END;
    return EOF;
  }

  /* search for the correct file in the list */
  twant = epoch->gpsSeconds + 1e-9 * epoch->gpsNanoSeconds;
  file1 = bsearch( &twant, stream->flist, stream->nfile,
      sizeof( *stream->flist ), flist_tcompare );
  if ( ! file1 ) /* no matching file: loop from beginning (should NOT happen) */
    file1 = stream->flist;

  /* find first file after time */
  for ( i = file1->ind; i < stream->nfile; ++i )
  {
    file1 = stream->flist + i;
    if ( epoch->gpsSeconds >= file1->t0
        && epoch->gpsSeconds < file1->t0 + file1->dt )
    {
      double tbeg;
      double tend;
      stream->fnum = i;
      stream->file = URLFrFileINew( file1 );
      if ( ! stream->file )
      {
        stream->state |= LAL_FR_ERR | LAL_FR_URL;
        return EOF;
      }
      if ( FrTOCReadFull( stream->file ) == NULL )
      {
        FrFileIEnd( stream->file );
        stream->file   = NULL;
        stream->state |= LAL_FR_ERR | LAL_FR_TOC;
        return EOF;
      }
      for ( stream->pos = 0; stream->pos < stream->file->toc->nFrame;
          ++stream->pos )
      {
        tbeg  = stream->file->toc->GTimeS[stream->pos];
        tbeg += 1e9 * stream->file->toc->GTimeN[stream->pos];
        tend  = tbeg + stream->file->toc->dt[stream->pos];
        if ( twant >= tbeg && twant < tend ) /* this is the frame */
        {
          goto found;
        }
        if ( twant < tbeg ) /* detect a gap */
        {
          stream->state |= LAL_FR_GAP;
          goto found;
        }
      }
      /* not in this frame file: close it */
      FrFileIEnd( stream->file );
      stream->file = NULL;
      stream->pos = 0;
    }
    if ( epoch->gpsSeconds < file1->t0 ) /* detect a gap */
    {
      stream->state |= LAL_FR_GAP;
      stream->fnum   = i;
      stream->file   = URLFrFileINew( file1 );
      stream->pos    = 0;
      if ( ! stream->file )
      {
        stream->state |= LAL_FR_ERR | LAL_FR_URL;
        return EOF;
      }
      if ( FrTOCReadFull( stream->file ) == NULL )
      {
        FrFileIEnd( stream->file );
        stream->file   = NULL;
        stream->state |= LAL_FR_ERR | LAL_FR_TOC;
        return EOF;
      }
      goto found;
    }
  }
found:

  /* set time of stream */
  if ( stream->state & LAL_FR_GAP )
  {
    stream->epoch.gpsSeconds     = stream->file->toc->GTimeS[stream->pos];
    stream->epoch.gpsNanoSeconds = stream->file->toc->GTimeN[stream->pos];
  }
  else
  {
    stream->epoch = *epoch;
  }

  return 0;
}

static int fr_tell( FrStream *stream, LIGOTimeGPS *epoch )
{
  if ( ! stream || ! epoch )
  {
    errno = EFAULT;
    return EOF;
  }
  *epoch = stream->epoch;
  return 0;
}

static int fr_getpos( FrStream *stream, FrPos *position )
{
  if ( ! stream || ! position )
  {
    errno = EFAULT;
    return EOF;
  }
  position->epoch = stream->epoch;
  position->fnum  = stream->fnum;
  position->pos   = stream->pos;
  return 0;
}

static int fr_setpos( FrStream *stream, FrPos *position )
{
  if ( ! stream || ! position )
  {
    errno = EFAULT;
    return EOF;
  }
  /* clear EOF or GAP states; preserve ERR state */
  if ( stream->state & LAL_FR_ERR )
    stream->state = LAL_FR_ERR;
  else
    stream->state = LAL_FR_OK;
  stream->epoch = position->epoch;
  if ( stream->fnum != position->fnum )
  {
    if ( stream->file )
    {
      FrFileIEnd( stream->file );
      stream->file = NULL;
    }
    if ( position->fnum >= stream->fnum )
    {
      errno = EINVAL;
      stream->fnum  = stream->nfile;
      stream->state |= LAL_FR_END;
      return EOF;
    }
    stream->fnum = position->fnum;
    stream->file = URLFrFileINew( stream->flist + stream->fnum );
    if ( ! stream->file )
    {
      stream->state |= LAL_FR_ERR | LAL_FR_URL;
      return EOF;
    }
  }
  if ( ! stream->file->toc )
  {
    if ( FrTOCReadFull( stream->file ) == NULL )
    {
      FrFileIEnd( stream->file );
      stream->file   = NULL;
      stream->state |= LAL_FR_ERR | LAL_FR_TOC;
      return EOF;
    }
  }
  stream->pos = position->pos;
  if ( stream->pos > stream->file->toc->nFrame )
  {
    errno = EINVAL;
    stream->state |= LAL_FR_ERR;
    return EOF;
  }
  return 0;
}


/* compare routine for qsort to sort flists in increasing order of start time */
static int flist_compare( const void *p, const void *q )
{
  const FrFileInfo *file1 = p;
  const FrFileInfo *file2 = q;
  if ( file1->t0 < file2->t0 )
    return -1;
  if ( file1->t0 > file2->t0 )
    return 1;
  return 0;
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
  list = *output = LALCalloc( cache->numFrameFiles + 1, sizeof( *list ) );
  list[cache->numFrameFiles].ind = -1; /* indicates end */
  for ( i = 0; i < cache->numFrameFiles; ++i )
  {
    FrStat *file = cache->frameFiles + i;
    list[i].url = LALMalloc( strlen( file->url ) + 1 );
    if ( ! list[i].url )
    {
      free_flist( *output );
      LALFree( *output );
      *output = NULL;
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
      frfile = URLFrFileINew( list + i );
      if ( ! frfile )
      {
        free_flist( *output );
        LALFree( *output );
        *output = NULL;
        ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
      }
      if ( FrTOCReadFull( frfile ) == NULL )
      {
        FrFileIEnd( frfile );
        free_flist( *output );
        LALFree( *output );
        ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
      }
      /* TODO: loop over frames */
      list[i].t0 = floor( frfile->toc->GTimeS[0] );
      list[i].dt = ceil( (double)frfile->toc->GTimeS[0]
          + 1e-9 * (double)frfile->toc->GTimeN[0]
          + frfile->toc->dt[0] ) - list[i].t0;
      FrFileIEnd( frfile );
    }
  }

  qsort( list, cache->numFrameFiles, sizeof( *list ), flist_compare );
  for ( i = 0; i < cache->numFrameFiles; ++i )
    list[i].ind = i;

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

  stream->nfile = cache->numFrameFiles;
  LALCreateFrFileList( status->statusPtr, &stream->flist, cache );
  BEGINFAIL( status )
  {
    LALFree( *output );
    *output = NULL;
  }
  ENDFAIL( status );

  FrLibIni( NULL, stderr, 0 );

  stream->file = URLFrFileINew( stream->flist );
  if ( ! stream->file )
  {
    free_flist( stream->flist );
    LALFree( *output );
    *output = NULL;
    ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
  }

  if ( FrTOCReadFull( stream->file ) == NULL )
  {
    FrFileIEnd( stream->file );
    free_flist( stream->flist );
    LALFree( *output );
    *output = NULL;
    ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
  }

  stream->epoch.gpsSeconds     = stream->file->toc->GTimeS[0];
  stream->epoch.gpsNanoSeconds = stream->file->toc->GTimeN[0];

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
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( *stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  fr_close( *stream );
  LALFree( *stream );
  *stream = NULL;
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
  *end = fr_state( stream ) & LAL_FR_END;
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
  if ( fr_rewind( stream ) )
  {
    if ( stream->state & LAL_FR_URL ) /* problem was in opening a file */
    {
      ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
    }
    if ( stream->state & LAL_FR_TOC ) /* problem was in reading a file */
    {
      ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
    }
  }
  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrNext(
    LALStatus   *status,
    FrStream    *stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALNextFrame", FRAMESTREAMC );  
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  if ( stream->state & LAL_FR_END )
  {
    ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
  }

  if ( fr_next( stream ) )
  {
    if ( stream->state & LAL_FR_ERR )
    {
      if ( stream->state & LAL_FR_URL ) /* must have failed to open a file */
      {
        ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
      }
      if ( stream->state & LAL_FR_TOC ) /* must have failed to read a file */
      {
        ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
      }
    }
  }

  if ( stream->state & LAL_FR_GAP )
  {
    LALInfo( status, "Gap in frame data." );
  }

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
  INITSTATUS( status, "LALFrSeek", FRAMESTREAMC );  
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( epoch, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }

  if ( fr_seek( stream, epoch ) )
  {
    if ( stream->state & LAL_FR_ERR )
    {
      if ( stream->state & LAL_FR_URL ) /* must have failed to open a file */
      {
        ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
      }
      if ( stream->state & LAL_FR_TOC ) /* must have failed to read a file */
      {
        ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
      }
    }
    else
    {
      if ( stream->state & LAL_FR_GAP ) /* too early */
      {
        LALWarning( status, "Requested time before first frame." );
        RETURN( status );
      }
      if ( stream->state & LAL_FR_END ) /* too late */
      {
        LALWarning( status, "Requested time after last frame." );
        RETURN( status );
      }
    }
  }
  
  if ( stream->state & LAL_FR_GAP )
  {
    LALWarning( status, "Requested time in gap in frame data." );
  }
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
  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  fr_tell( stream, epoch );
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
  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  fr_getpos( stream, position );
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
  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  errno = 0;
  if ( fr_setpos( stream, position ) )
  {
    if ( stream->state & LAL_FR_ERR )
    {
      if ( stream->state & LAL_FR_URL ) /* must have failed to open a file */
      {
        ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
      }
      if ( stream->state & LAL_FR_TOC ) /* must have failed to read a file */
      {
        ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
      }
    }
  }
  RETURN( status );
}
