/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Kipp Cannon
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
 * The routine \texttt{LALFrSetMode()} is used to change the operating mode
 * of a frame stream, which determines how the routines try to accomodate
 * gaps in data and requests for times when there is no data (e.g., before
 * the beginning of the data, after the end of the data, or in some missing
 * data).  The default mode, which is given the value
 * \verb+LAL_FR_DEFAULT_MODE+, prints warnings if a time requested
 * corresponds to a time when there is no data (but then skips to the first
 * avaliable data) and prints an info message when a gap in the data occurs
 * (but then skips beyond the gap).  This default mode is equal to the
 * combination
 * \verb+LAL_FR_VERBOSE_MODE | LAL_FR_IGNOREGAP_MODE | LAL_FR_IGNORETIME_MODE+
 * where \verb+LAL_FR_VERBOSE_MODE+ is equal to the combination
 * \verb+LAL_FR_TIMEWARN_MODE | LAL_FR_GAPINFO_MODE+.  Use
 * \verb+LAL_FR_VERBOSE_MODE+ to print out warnings when requesting times
 * with no data and print out an info message when a gap in the data is
 * encountered.  Unless the mode is supplemented with
 * \verb+LAL_FR_IGNOREGAP_MODE+, gaps encountered in the data will cause
 * a routine to exit with a non-zero status code; similarly,
 * \verb+LAL_FR_IGNORETIME_MODE+ prevents routines from failing if a time
 * when there is not data is requested.  Set the mode to
 * \verb+LAL_FR_SILENT_MODE+ to suppress the warning and info messages but
 * still cause routines to fail when data is not available.
 * Note: the default value \verb+LAL_FR_DEFAULT_MODE+ is assumed initially,
 * but this is not necessarily the recommended mode --- it is adopted for
 * compatibility reasons.
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

#include <config.h>
#include <unistd.h>
#ifndef HAVE_GETHOSTNAME_PROTOTYPE
int gethostname(char *name, int len);
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/Date.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALFrameIO.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>

NRCSID( FRAMESTREAMC, "$Id$" );


/*
 *
 * These functions are for internal use.
 *
 */



#include <errno.h>


/*
 *
 * Routine to Open a Frame File URL.
 *
 */


static struct FrFile *URLFrFileINew( FrFileInfo *file, int dochecksum )
{
  struct FrFile *frfile;
  if ( ! file )
    return NULL;
  frfile = XLALFrOpenURL( file->url );
  if ( frfile && dochecksum ) {
    int code;
    code = XLALFrFileCheckSum( frfile );
    switch ( code ) {
      case 0:
        break;
      case -1:
        XLALPrintError( "XLAL Error: wrong checksum in file %s\n", file->url );
        FrFileIEnd( frfile );
        return NULL;
      case 1:
        XLALPrintWarning( "XLAL Warning: missing checksum in file %s\n", file->url );
        break;
      default: /* unhandled case */
        FrFileIEnd( frfile );
        XLAL_ERROR_NULL( "URLFrFileINew", XLAL_EERR );
    }
  }
  return frfile;
}



/*
 *
 * Local Routines for Manipulating Frame File Lists.
 *
 */


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

/* compare routine to identify duplicate entries: two files with the same
 * start time and duration and the same basename (but possibly different
 * urls) ... the purpose is to take two entries such as:
 *
 *   file://localhost/archive/frames/S5/strain-LX/LHO/H-H2_RDS_C02_LX-8157/H-H2_RDS_C02_LX-815750954-128.gwf
 *   file://localhost/data/node134/frames/S5/strain-LX/LHO/H-H2_RDS_C02_LX-8157/H-H2_RDS_C02_LX-815750954-128.gwf
 *
 * and get rid of one of them.  */
static int flist_compare_uniq( const void *p, const void *q )
{
  const FrFileInfo *file1 = p;
  const FrFileInfo *file2 = q;
  char *fname1;
  char *fname2;
  int c;
  if ( file1->t0 < file2->t0 )
    return -1;
  if ( file1->t0 > file2->t0 )
    return 1;
  fname1 = strrchr( file1->url, '/' );
  fname1 = fname1 ? fname1 + 1 : file1->url;
  fname2 = strrchr( file2->url, '/' );
  fname2 = fname2 ? fname2 + 1 : file2->url;
  c = strcmp( fname1, fname2 );
  return c;
}


/* compare routine for binary search bsearch to locate the first file
 * containing wanted time, or the first file after a gap if the wanted time
 * occurs during a gap */
static int flist_tcompare( const void *key, const void *ptr )
{
  const FrFileInfo *file = ptr; /* this file */
  const FrFileInfo *prev = file - 1; /* the previous file */
  double twant = *((const double*)key);
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


static void free_flist( FrFileInfo *flist )
{
  FrFileInfo *p = flist;
  if ( ! p )
    return;
  while ( p->ind > -2 )
  {
    if ( p->url )
      LALFree( p->url );
    memset( p, 0, sizeof( *p ) );
    ++p;
  }
  LALFree( flist );
  return;
}


/* create a frame file list from a cache */
static UINT4 create_flist( FrFileInfo **plist, FrCache *cache )
{
  FrFileInfo *list;
  INT4 i, j, n;
  n = cache->numFrameFiles;
  list = *plist = LALCalloc( n + 1, sizeof( *list ) );
  list[n].ind = -2; /* indicates end */
  for ( i = 0; i < n; ++i )
  {
    FrStat *file = cache->frameFiles + i;
    list[i].url = LALMalloc( strlen( file->url ) + 1 );
    if ( ! list[i].url )
    {
      free_flist( list );
      LALFree( list );
      *plist = NULL;
      return 0;
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
      double t1;
      frfile = URLFrFileINew( list + i, 0 );
      if ( ! frfile )
      {
        free_flist( list );
        LALFree( list );
	*plist = NULL;
        return 0;
      }
      if ( FrTOCReadFull( frfile ) == NULL )
      {
        XLALPrintError( "XLAL Error: could not read frame TOC from URL %s\n",
            list[i].url );
        free_flist( list );
        LALFree( list );
	*plist = NULL;
        return 0;
      }
      /* TODO: loop over frames */
      list[i].t0 = floor( frfile->toc->GTimeS[0] );
      t1  = (double)frfile->toc->GTimeS[frfile->toc->nFrame - 1];
      t1 += (double)frfile->toc->GTimeN[frfile->toc->nFrame - 1] * 1e-9;
      t1 += (double)frfile->toc->dt[frfile->toc->nFrame - 1];
      list[i].dt = ceil( t1 ) - list[i].t0;
      FrFileIEnd( frfile );
    }
  }

  qsort( list, cache->numFrameFiles, sizeof( *list ), flist_compare );

  /* index lines and look for duplicate entries (and get rid of them) */
  i = 0;
  for ( j = 0; j < n; ++j )
  {
    FrFileInfo swap;
    swap = list[i];
    list[i] = list[j];
    list[j] = swap;
    if ( j+1 == n || flist_compare_uniq( list + i, list + j + 1 ) )
    {
      list[i].ind = i;
      ++i;
    }
    else if ( j+1 < n )
    {
      XLALPrintWarning( "XLAL Warning: duplicate entry in file list\n"
          "\tkeep:    %s\n\tdiscard: %s\n", list[j+1].url, list[i].url );
    }
  }
  n = i;
  for ( ; i < (INT4)cache->numFrameFiles; ++i )
      list[i].ind = -1; /* beyond the last file in the list */
  list[cache->numFrameFiles].ind = -2; /* last item in list */

  return n;
}


/*
 *
 * XLAL Routines.
 *
 */


FrStream * XLALFrCacheOpen( FrCache *cache )
{
  static const char func[] = "XLALFrCacheOpen";
  FrStream *stream;

  if ( ! cache )
    XLAL_ERROR_NULL( func, XLAL_EFAULT );

  stream = LALCalloc( 1, sizeof( *stream ) );
  if ( ! stream )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  stream->mode = LAL_FR_DEFAULT_MODE;
  stream->nfile = create_flist( &stream->flist, cache );
  if ( ! stream->flist )
  {
    LALFree( stream );
    XLAL_ERROR_NULL( func, XLAL_EIO ); /* assume it was an I/O error */
  }

  stream->file = URLFrFileINew( stream->flist, 0 );
  if ( ! stream->file )
  {
    free_flist( stream->flist );
    LALFree( stream );
    XLAL_ERROR_NULL( func, XLAL_EIO );
  }

  if ( FrTOCReadFull( stream->file ) == NULL )
  {
    XLALPrintError( "XLAL Error - %s: could not read frame TOC from URL %s\n",
        func, stream->flist->url );
    FrFileIEnd( stream->file );
    free_flist( stream->flist );
    LALFree( stream );
    XLAL_ERROR_NULL( func, XLAL_EIO );
  }

  stream->epoch.gpsSeconds     = stream->file->toc->GTimeS[0];
  stream->epoch.gpsNanoSeconds = stream->file->toc->GTimeN[0];

  return stream;
}


FrStream * XLALFrOpen( const char *dirname, const char *pattern )
{
  static const char func[] = "XLALFrOpen";
  FrStream *stream;
  FrCache *cache;

  cache = XLALFrGenerateCache( dirname, pattern );
  if ( ! cache )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  stream = XLALFrCacheOpen( cache );
  if ( ! stream )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  XLALFrDestroyCache( cache );
  return stream;
}


int XLALFrClose( FrStream *stream )
{
  static const char func[] = "XLALFrClose";
  if ( ! stream )
    XLAL_ERROR( func, XLAL_EFAULT );
  FrFileIEnd( stream->file );
  free_flist( stream->flist );
  LALFree( stream );
  return 0;
}


int XLALFrSetMode( FrStream *stream, int mode )
{
  static const char func[] = "XLALFrSetMode";
  if ( ! stream )
    XLAL_ERROR( func, XLAL_EFAULT );
  stream->mode = mode;
  /* if checksum mode is turned on, do checksum on current file */
  if ( (mode & LAL_FR_CHECKSUM_MODE) && (stream->file) )
    return XLALFrFileCheckSum( stream->file );
  return 0;
}


int XLALFrState( FrStream *stream )
{
  static const char func[] = "XLALFrState";
  if ( ! stream )
    XLAL_ERROR( func, XLAL_EFAULT );
  return stream->state;
}


int XLALFrClearErr( FrStream *stream )
{
  static const char func[] = "XLALFrClearErr";
  if ( ! stream )
    XLAL_ERROR( func, XLAL_EFAULT );
  stream->state = LAL_FR_OK;
  return 0;
}


int XLALFrRewind( FrStream *stream )
{
  static const char func[] = "XLALFrRewind";
  if ( ! stream )
    XLAL_ERROR( func, XLAL_EFAULT );

  if ( stream->file )
  {
    FrFileIEnd( stream->file );
    stream->file = NULL;
  }
  stream->pos   = 0;
  stream->fnum  = 0;
  stream->file  = URLFrFileINew( stream->flist, stream->mode & LAL_FR_CHECKSUM_MODE );
  stream->state = LAL_FR_OK;
  if ( ! stream->file )
  {
    stream->state |= LAL_FR_ERR | LAL_FR_URL;
    XLAL_ERROR( func, XLAL_EIO );
  }
  if ( FrTOCReadFull( stream->file ) == NULL )
  {
    XLALPrintError( "XLAL Error - %s: could not read frame TOC from URL %s\n",
        func, stream->flist->url );
    FrFileIEnd( stream->file );
    stream->file   = NULL;
    stream->state |= LAL_FR_ERR | LAL_FR_TOC;
    XLAL_ERROR( func, XLAL_EIO );
  }
  stream->epoch.gpsSeconds     = stream->file->toc->GTimeS[0];
  stream->epoch.gpsNanoSeconds = stream->file->toc->GTimeN[0];
  return 0;
}


int XLALFrNext( FrStream *stream )
{
  static const char func[] = "XLALFrNext";

  /* timing accuracy: tenth of a sample interval for a 16kHz fast channel */
  const INT8 tacc = (INT8)floor( 0.1 * 1e9 / 16384.0 );
  const char *url1;
  const char *url2;
  int pos1;
  int pos2;
  INT8 tnow = 0;
  INT8 texp = 0;
  INT8 tact;

  if ( ! stream )
    XLAL_ERROR( func, XLAL_EFAULT );

  if ( stream->state & LAL_FR_END )
    return 1; /* end code */

  /* turn off gap bit */
  stream->state &= ~LAL_FR_GAP;

  url2 = url1 = stream->flist[stream->fnum].url;
  pos2 = pos1 = stream->pos;
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
    pos2 = stream->pos;
  }

  /* open a new file if necessary */
  if ( ! stream->file )
  {
    stream->pos = 0;
    if ( stream->fnum >= stream->nfile )
    {
      stream->state |= LAL_FR_END;
      return 1;
    }
    stream->pos  = 0;
    stream->file = URLFrFileINew( stream->flist + stream->fnum, stream->mode & LAL_FR_CHECKSUM_MODE );
    url2 = stream->flist[stream->fnum].url;
    if ( ! stream->file )
    {
      stream->state |= LAL_FR_ERR | LAL_FR_URL;
      XLAL_ERROR( func, XLAL_EIO );
    }
    pos2 = stream->pos;
  }

  /* open TOC if necessary */
  if ( ! stream->file->toc )
  {
    if ( FrTOCReadFull( stream->file ) == NULL )
    {
      XLALPrintError( "XLAL Error - %s: could not read frame TOC from URL %s\n",
          func, stream->flist[stream->fnum].url );
      FrFileIEnd( stream->file );
      stream->file   = NULL;
      stream->state |= LAL_FR_ERR | LAL_FR_TOC;
      XLAL_ERROR( func, XLAL_EIO );
    }
  }

  /* compute actual start time of this new frame */
  tact  = (INT8)1000000000 * (INT8)stream->file->toc->GTimeS[stream->pos];
  tact += (INT8)stream->file->toc->GTimeN[stream->pos];
  stream->epoch.gpsSeconds     = stream->file->toc->GTimeS[stream->pos];
  stream->epoch.gpsNanoSeconds = stream->file->toc->GTimeN[stream->pos];

  if ( abs( texp - tact ) > tacc ) /* there is a gap */
  {
    stream->state |= LAL_FR_GAP;
    if ( stream->mode & LAL_FR_GAPINFO_MODE )
    {
      XLALPrintInfo( "XLAL Info - %s: gap in frame data between times "
          "%.6f and %.6f\n", func, 1e-9 * texp, 1e-9 * tact );
    }
    if ( ! ( stream->mode & LAL_FR_IGNOREGAP_MODE ) )
    {
      XLALPrintError( "XLAL Error - %s: gap in frame data\n", func );
      XLALPrintError( "XLAL Error - %s: time %.6f "
          "is end of frame %d of file URL %s\n", func, 1e-9*texp, pos1, url1 );
      XLALPrintError( "XLAL Error - %s: time %.6f "
          "is start of frame %d of file URL %s\n", func, 1e-9*tact, pos2, url2);
      XLAL_ERROR( func, XLAL_ETIME );
    }
    return 2; /* gap code */
  }

  return 0;
}


int XLALFrSeek( FrStream *stream, const LIGOTimeGPS *epoch )
{
  static const char func[] = "XLALFrSeek";
  FrFileInfo *fileinfo;
  double twant = XLALGPSGetREAL8( epoch );

  /* close file if one is open */
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

  /* is epoch before first file? */
  if ( epoch->gpsSeconds < stream->flist->t0 )
  {
    XLALFrRewind( stream );
    stream->state |= LAL_FR_GAP;
    /* is this reported as an error? */
    if ( ! ( stream->mode & LAL_FR_IGNORETIME_MODE ) )
    {
      /* FIXME:  if this is an error, should the stream state say so? */
      /*stream->state |= LAL_FR_ERR;*/
      XLAL_ERROR( func, XLAL_ETIME );
    }
    if ( stream->mode & LAL_FR_TIMEWARN_MODE )
      XLALPrintWarning( "XLAL Warning - %s: "
          "requested time %d before first frame\n", func, epoch->gpsSeconds );
    return 1; /* before first file code */
  }

  /* search for the correct file in the list */
  fileinfo = bsearch( &twant, stream->flist, stream->nfile,
      sizeof( *stream->flist ), flist_tcompare );
  if ( ! fileinfo )
    /* bsearch() failed:  brute-force search from beginning (should only
     * happen if requested time is past end of cache, but we'll figure that
     * out below) */
    fileinfo = stream->flist;

  /* find first file after time */
  while ( 1 )
  {
    if ( fileinfo - stream->flist >= (ptrdiff_t)stream->nfile )
    {
      /* end of file list, time not yet found */
      stream->fnum   = stream->nfile;
      stream->state |= LAL_FR_END;
      /* is this reported as an error? */
      if ( ! ( stream->mode & LAL_FR_IGNORETIME_MODE ) )
      {
        /* FIXME:  if this is an error, should the stream state say so? */
        /*stream->state |= LAL_FR_ERR;*/
        XLAL_ERROR( func, XLAL_ETIME );
      }
      if ( stream->mode & LAL_FR_TIMEWARN_MODE )
        XLALPrintWarning( "XLAL Warning - %s: "
            "requested time %d after last frame\n", func, epoch->gpsSeconds );
      return 2; /* after last file code */
    }
    if ( epoch->gpsSeconds >= fileinfo->t0
        && epoch->gpsSeconds < fileinfo->t0 + fileinfo->dt )
    {
      /* epoch matches this cache entry, check file contents */
      stream->fnum = fileinfo - stream->flist;
      stream->file = URLFrFileINew( fileinfo, stream->mode & LAL_FR_CHECKSUM_MODE );
      if ( ! stream->file )
      {
        stream->state |= LAL_FR_ERR | LAL_FR_URL;
        XLAL_ERROR( func, XLAL_EIO );
      }
      if ( FrTOCReadFull( stream->file ) == NULL )
      {
        XLALPrintError(
            "XLAL Error - %s: could not read frame TOC from URL %s\n",
            func, fileinfo->url );
        FrFileIEnd( stream->file );
        stream->file   = NULL;
        stream->state |= LAL_FR_ERR | LAL_FR_TOC;
        XLAL_ERROR( func, XLAL_EIO );
      }
      /* loop over frames */
      for ( stream->pos = 0; stream->pos < stream->file->toc->nFrame;
          ++stream->pos )
      {
        LIGOTimeGPS tbeg;
        XLALGPSSet( &tbeg, stream->file->toc->GTimeS[stream->pos], stream->file->toc->GTimeN[stream->pos] );
        if ( XLALGPSCmp( epoch, &tbeg ) >= 0 && XLALGPSDiff( epoch, &tbeg ) < stream->file->toc->dt[stream->pos] ) /* this is the frame */
        {
          goto found;
        }
        if ( XLALGPSCmp( epoch, &tbeg ) < 0 ) /* detect a gap */
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
    else if ( epoch->gpsSeconds < fileinfo->t0 ) /* detect a gap */
    {
      stream->state |= LAL_FR_GAP;
      stream->fnum   = fileinfo - stream->flist;
      stream->file   = URLFrFileINew( fileinfo, stream->mode & LAL_FR_CHECKSUM_MODE );
      stream->pos    = 0;
      if ( ! stream->file )
      {
        stream->state |= LAL_FR_ERR | LAL_FR_URL;
        XLAL_ERROR( func, XLAL_EIO );
      }
      if ( FrTOCReadFull( stream->file ) == NULL )
      {
        XLALPrintError(
            "XLAL Error - %s: could not read frame TOC from URL %s\n",
            func, fileinfo->url );
        FrFileIEnd( stream->file );
        stream->file   = NULL;
        stream->state |= LAL_FR_ERR | LAL_FR_TOC;
        XLAL_ERROR( func, XLAL_EIO );
      }
      goto found;
    }
    /* advance to next cache entry */
    fileinfo++;
  }
found:

  /* set time of stream */
  if ( stream->state & LAL_FR_GAP )
  {
    XLALGPSSet( &stream->epoch, stream->file->toc->GTimeS[stream->pos], stream->file->toc->GTimeN[stream->pos] );
    if ( stream->mode & LAL_FR_TIMEWARN_MODE )
      XLALPrintWarning( "XLAL Warning - %s: "
          "requested time %.6f in gap in frame data\n", func, twant );
    if ( ! ( stream->mode & LAL_FR_IGNORETIME_MODE ) )
      XLAL_ERROR( func, XLAL_ETIME );
    return 3; /* in a gap code */
  }

  stream->epoch = *epoch;
  return 0;
}


int XLALFrTell( LIGOTimeGPS *epoch, FrStream *stream )
{
  static const char func[] = "XLALFrTell";
  if ( ! epoch || ! stream )
    XLAL_ERROR( func, XLAL_EFAULT );
  *epoch = stream->epoch;
  return 0;
}


int XLALFrGetpos( FrPos *position, FrStream *stream )
{
  static const char func[] = "XLALFrGetpos";
  if ( ! position || ! stream )
    XLAL_ERROR( func, XLAL_EFAULT );
  position->epoch = stream->epoch;
  position->fnum  = stream->fnum;
  position->pos   = stream->pos;
  return 0;
}


int XLALFrSetpos( FrStream *stream, FrPos *position )
{
  static const char func[] = "XLALFrSetpos";
  if ( ! stream || ! position )
    XLAL_ERROR( func, XLAL_EFAULT );
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
      stream->fnum  = stream->nfile;
      stream->state |= LAL_FR_END;
      XLAL_ERROR( func, XLAL_EINVAL );
    }
    stream->fnum = position->fnum;
    stream->file = URLFrFileINew( stream->flist + stream->fnum, stream->mode & LAL_FR_CHECKSUM_MODE );
    if ( ! stream->file )
    {
      stream->state |= LAL_FR_ERR | LAL_FR_URL;
      XLAL_ERROR( func, XLAL_EIO );
    }
  }
  if ( ! stream->file->toc )
  {
    if ( FrTOCReadFull( stream->file ) == NULL )
    {
      XLALPrintError( "XLAL Error - %s: could not read frame TOC from URL %s\n",
          func, stream->flist[stream->fnum].url );
      FrFileIEnd( stream->file );
      stream->file   = NULL;
      stream->state |= LAL_FR_ERR | LAL_FR_TOC;
      XLAL_ERROR( func, XLAL_EIO );
    }
  }
  stream->pos = position->pos;
  if ( stream->pos > stream->file->toc->nFrame )
  {
    stream->state |= LAL_FR_ERR;
    XLAL_ERROR( func, XLAL_EINVAL );
  }
  return 0;
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
  ASSERT( cache, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( output, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( ! *output, status, FRAMESTREAMH_ENNUL, FRAMESTREAMH_MSGENNUL );

  stream = *output = XLALFrCacheOpen( cache );
  if ( ! stream )
  {
    int errnum = xlalErrno;
    XLALClearErrno();
    switch ( errnum )
    {
      case XLAL_ENOMEM:
        ABORT( status, FRAMESTREAMH_EALOC, FRAMESTREAMH_MSGEALOC );
      case XLAL_EIO:
        ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
      default:
        ABORTXLAL( status );
    }
  }

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
  XLALFrClose( *stream );
  *stream = NULL;
  RETURN( status );
}

/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrSetMode(
    LALStatus *status,
    INT4       mode,
    FrStream  *stream
    )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALFrSetMode", FRAMESTREAMC );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  stream->mode = mode;
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
  *end = XLALFrState( stream ) & LAL_FR_END;
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
  if ( XLALFrRewind( stream ) )
  {
    XLALClearErrno();
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
  CHAR frErrMsg[1024];
  int code;
  INITSTATUS( status, "LALFrNext", FRAMESTREAMC );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );

  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }
  if ( stream->state & LAL_FR_END )
  {
    ABORT( status, FRAMESTREAMH_EDONE, FRAMESTREAMH_MSGEDONE );
  }

  code = XLALFrNext( stream );
  if ( code < 0 )
  {
    XLALClearErrno();
    if ( stream->state & LAL_FR_ERR )
    {
      if ( stream->state & LAL_FR_URL ) /* must have failed to open a file */
      {
        snprintf( frErrMsg, sizeof(frErrMsg)/sizeof(*frErrMsg),
            "Could not open URL %s\n", (stream->flist+stream->fnum)->url );
        LALError( status, frErrMsg );
        ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
      }
      if ( stream->state & LAL_FR_TOC ) /* must have failed to read a file */
      {
        snprintf( frErrMsg, sizeof(frErrMsg)/sizeof(*frErrMsg),
            "Could not read TOC from %s\n", (stream->flist+stream->fnum)->url );
        LALError( status, frErrMsg );
        ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
      }
    }
    else /* must be a gap error */
    {
      ABORT( status, FRAMESTREAMH_EDGAP, FRAMESTREAMH_MSGEDGAP );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="FrameStreamCP"> */
void
LALFrSeek(
    LALStatus         *status,
    const LIGOTimeGPS *epoch,
    FrStream          *stream
    )
{ /* </lalVerbatim> */
  CHAR frErrMsg[1024];
  int code;
  INITSTATUS( status, "LALFrSeek", FRAMESTREAMC );
  ASSERT( stream, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  ASSERT( epoch, status, FRAMESTREAMH_ENULL, FRAMESTREAMH_MSGENULL );
  if ( stream->state & LAL_FR_ERR )
  {
    ABORT( status, FRAMESTREAMH_ERROR, FRAMESTREAMH_MSGERROR );
  }

  code = XLALFrSeek( stream, epoch );
  if ( code < 0 )
  {
    XLALClearErrno();
    if ( stream->state & LAL_FR_ERR ) /* a file error */
    {
      if ( stream->state & LAL_FR_URL ) /* must have failed to open a file */
      {
        snprintf( frErrMsg, sizeof(frErrMsg)/sizeof(*frErrMsg),
            "Could not open URL %s\n", (stream->flist+stream->fnum)->url );
        LALError( status, frErrMsg );
        ABORT( status, FRAMESTREAMH_EOPEN, FRAMESTREAMH_MSGEOPEN );
      }
      if ( stream->state & LAL_FR_TOC ) /* must have failed to read a file */
      {
        snprintf( frErrMsg, sizeof(frErrMsg)/sizeof(*frErrMsg),
            "Could not read TOC from %s\n", (stream->flist+stream->fnum)->url );
        LALError( status, frErrMsg );
        ABORT( status, FRAMESTREAMH_EREAD, FRAMESTREAMH_MSGEREAD );
      }
    }
    else /* must be too early, too late, or in a gap */
    {
      ABORT( status, FRAMESTREAMH_ETREQ, FRAMESTREAMH_MSGETREQ );
    }
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
  XLALFrTell( epoch, stream );
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
  XLALFrGetpos( position, stream );
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
  if ( XLALFrSetpos( stream, position ) )
  {
    XLALClearErrno();
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
