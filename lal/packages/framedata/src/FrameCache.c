/**** <lalVerbatim file="FrameCacheCV">
 * Author: Jolien D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{FrameCache.c}}
 *
 * Routines for importing, exporting, generating, and manipulating frame
 * catalogs and cache structures.
 *
 * \subsubsection*{Prototypes}
 * \input{FrameCacheCP}
 * \idx{LALFrCacheImport}
 * \idx{LALFrCacheExport}
 * \idx{LALFrCacheSieve}
 * \idx{LALFrCacheGenerate}
 *
 * \subsubsection*{Description}
 *
 * A frame catalogue file has several pieces of metadata, including:
 * \begin{description}
 * \item[source] the source identifier, often a combination of upper-case
 *     letters representing the detector sites (e.g., `H' for Hanford, `L'
 *     for Livingston) from which the frame data was generated.
 * \item[description] the description identifier, often a single upper-case
 *     letter describing what kind of frame data is present (e.g., `R' for
 *     raw data, `M' for minute-trend data).
 * \item[GPS-time] the GPS time in seconds (rounded down) of the start of the
 *     data contained in the frame file.
 * \item[duration] the difference between the GPS time in seconds (rounded up)
 *     of the end of the data contained in the frame file and the GPS time in
 *     seconds (rounded down) of the start of the data contained in the frame
 *     file.
 * \item[URL] the URL of the frame data.  If the protocol (i.e.,
 *     \verb+file://+ or \verb+http://+) is missing then it is assumed that
 *     this is the absolute or relative path of a frame file.  File URLs
 *     must have a fully-qualified domain name or the word \verb+localhost+
 *     if on the localhost (or else localhost is assumed if absent).
 *     Examples:
 *     \begin{verbatim}
 *     file://machine.university.edu/usr/share/lal/F-TEST-600000000-60.gwf
 *     file://localhost/usr/share/lal/F-TEST-600000000-60.gwf
 *     file:///usr/share/lal/F-TEST-600000000-60.gwf
 *     /usr/share/lal/F-TEST-600000000-60.gwf
 *     F-TEST-600000000-60.gwf
 *     \end{verbatim}
 * \end{description}
 * Other types of metadata, such as the md5 checksum and the file size, are not
 * used by these LAL routines.  The format of the catalogues is
 * \begin{verbatim}
 * source description GPS-time duration URL (ignored additional fields)
 * \end{verbatim}
 * for example:
 * \begin{verbatim}
 * F TEST 600000000 60 file://localhost/usr/share/lal/F-TEST-600000000-60.gwf
 * F TEST 600000060 60 file://localhost/usr/share/lal/F-TEST-600000060-60.gwf
 * F TEST 600000120 60 file://localhost/usr/share/lal/F-TEST-600000120-60.gwf
 * \end{verbatim}
 * If any data is missing or unknown, it is represented in the catalogue by a
 * single hyphen (\texttt{-}).
 *
 * The routine \texttt{LALFrCacheImport} reads in a specified frame catalogue
 * file and creates a frame cache.  The routine \texttt{LALFrCacheExport}
 * exports a frame cache as a frame catalogue file.  The routine
 * \texttt{LALFrCacheGenerate} scans a colon-delimited list of directory paths
 * (or \texttt{.} if \texttt{NULL}) for files that match a given glob pattern
 * (default is \texttt{*.gwf} if \texttt{NULL}), and uses these to generate a
 * frame cache (with the metadata extracted from the file names).  The routine
 * \texttt{LALFrCacheSieve} applies regular expression filters to various
 * pieces of metadata to distill a frame cache into a sorted sub-cache of frame
 * files of interst.  The routine \texttt{LALDestroyFrCache} destroys a frame
 * cache.
 *
 * \vfill{\footnotesize\input{FrameCacheCV}}
 *
 **** </lalLaTeX> */




#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/param.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <regex.h>
#include <glob.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/FrameCache.h>

NRCSID( FRAMECACHEC, "$Id$" );

#define TMPSTRLEN 15
#define XSTR( x ) #x
#define STR( x ) XSTR( x )

/* <lalVerbatim file="FrameCacheCP"> */
void LALFrCacheImport(
    LALStatus   *status,
    FrCache    **output,
    const CHAR  *fname
    )
{ /* </lalVerbatim> */
  UINT4 numLines = 0;
  CHAR  line[1024];
  FILE *fp;
  FrCache *cache;

  INITSTATUS( status, "LALFrCacheImport", FRAMECACHEC );
  ATTATCHSTATUSPTR( status );
  ASSERT( output, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );
  ASSERT( ! *output, status, FRAMECACHEH_ENNUL, FRAMECACHEH_MSGENNUL );
  ASSERT( fname, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );

  fp = LALFopen( fname, "r" );
  if ( ! fp )
  {
    ABORT( status, FRAMECACHEH_EIEIO, FRAMECACHEH_MSGEIEIO );
  }
  numLines = 0;
  while ( fgets( line, sizeof( line ), fp ) )
  {
    if ( strlen( line ) > sizeof( line ) - 2 )
    {
      ABORT( status, FRAMECACHEH_ELINE, FRAMECACHEH_MSGELINE );
    }
    ++numLines;
  }
  if ( ! numLines )
  {
    ABORT( status, FRAMECACHEH_EIEIO, FRAMECACHEH_MSGEIEIO );
  }

  cache = *output = LALCalloc( 1, sizeof( **output ) );
  if ( ! cache )
  {
    ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
  }
  cache->numFrameFiles = numLines;
  cache->frameFiles = LALCalloc( numLines, sizeof( *cache->frameFiles ) );
  if ( ! cache->frameFiles )
  {
    LALFree( *output );
    *output = NULL;
    ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
  }

  rewind( fp );
  numLines = 0;
  while ( fgets( line, sizeof( line ), fp ) )
  {
    FrStat *file = cache->frameFiles + numLines++;
    char t0[TMPSTRLEN + 1];
    char dt[TMPSTRLEN + 1];
    int src0, src1;
    int dsc0, dsc1;
    int url0, url1;
    int c;
    sscanf( line, "%n%*[a-zA-Z0-9_+#-]%n %n%*[a-zA-Z0-9_+#-]%n %*"
        STR( TMPSTRLEN ) "s %*" STR( TMPSTRLEN ) "s %n%*s%n",
        &src0, &src1, &dsc0, &dsc1, &url0, &url1 );
    file->source = LALMalloc( src1 - src0 + 1 );
    if ( ! file->source )
    {
      TRY( LALDestroyFrCache( status->statusPtr, output ), status );
      ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
    }
    file->description = LALMalloc( dsc1 - dsc0 + 1 );
    if ( ! file->description )
    {
      TRY( LALDestroyFrCache( status->statusPtr, output ), status );
      ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
    }
    file->url = LALMalloc( url1 - url0 + 1 );
    if ( ! file->url )
    {
      TRY( LALDestroyFrCache( status->statusPtr, output ), status );
      ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
    }
    c = sscanf( line, "%[a-zA-Z0-9_+#-] %[a-zA-Z0-9_+#-] %" STR( TMPSTRLEN )
        "s %" STR( TMPSTRLEN ) "s %s", file->source, file->description,
        t0, dt, file->url );
    if ( strchr( file->source, '-' ) )
    {
      LALFree( file->source );
      file->source = NULL;
    }
    if ( strchr( file->description, '-' ) )
    {
      LALFree( file->description );
      file->description = NULL;
    }
    if ( *file->url == '-' )
    {
      LALFree( file->url );
      file->url = NULL;
    }
    if ( strchr( t0, '-' ) || 0 == ( file->startTime = atoi( t0 ) ) )
    {
      file->startTime = 0; /* invalid value */
    }
    if ( strchr( dt, '-' ) || 0 == ( file->duration = atoi( dt ) ) )
    {
      file->duration = 0; /* invalid value */
    }
  }

  LALFclose( fp );
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="FrameCacheCP"> */
void LALFrCacheExport(
    LALStatus  *status,
    FrCache    *cache,
    const CHAR *fname
    )
{ /* </lalVerbatim> */
  UINT4 i;
  FILE *fp;

  INITSTATUS( status, "LALFrCacheExport", FRAMECACHEC );
  ASSERT( cache, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );
  ASSERT( fname, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );

  fp = LALFopen( fname, "w" );
  if ( ! fp )
  {
    ABORT( status, FRAMECACHEH_EIEIO, FRAMECACHEH_MSGEIEIO );
  }
  for ( i = 0; i < cache->numFrameFiles; ++i )
  {
    FrStat *file = cache->frameFiles + i;
    char t0[TMPSTRLEN + 1];
    char dt[TMPSTRLEN + 1];
    int c;
    if ( file->startTime > 0 )
    {
      LALSnprintf( t0, sizeof( t0 ), "%d", file->startTime );
    }
    else
    {
      strncpy( t0, "-", sizeof( t0 ) );
    }
    if ( file->duration > 0 )
    {
      LALSnprintf( dt, sizeof( dt ), "%d", file->duration );
    }
    else
    {
      strncpy( dt, "-", sizeof( dt ) );
    }
    c = fprintf( fp, "%s %s %s %s %s\n", file->source ? file->source : "-",
        file->description ? file->description : "-", t0, dt,
        file->url ? file->url : "-" );
    if ( c < 1 )
    {
      ABORT( status, FRAMECACHEH_EIEIO, FRAMECACHEH_MSGEIEIO );
    }
  }
  LALFclose( fp );
  RETURN( status );
}

/* <lalVerbatim file="FrameCacheCP"> */
void
LALDestroyFrCache(
    LALStatus  *status,
    FrCache   **cache
    )
{ /* </lalVerbatim> */
  UINT4 i;
  INITSTATUS( status, "LALDestroyFrCache", FRAMECACHEC );
  ASSERT( cache, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );
  ASSERT( *cache, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );

  for ( i = 0; i < (*cache)->numFrameFiles; ++i )
  {
    FrStat *file = (*cache)->frameFiles + i;
    if ( file->source )
    {
      LALFree( file->source );
    }
    if ( file->description )
    {
      LALFree( file->description );
    }
    if ( file->url )
    {
      LALFree( file->url );
    }
    memset( file, 0, sizeof( *file ) );
  }

  if ( (*cache)->numFrameFiles ) 
    LALFree( (*cache)->frameFiles );
  LALFree( *cache );
  *cache = NULL;

  RETURN( status );
}

static int FrStatCompare( const void *p1, const void *p2 )
{
  const FrStat *file1 = p1;
  const FrStat *file2 = p2;
  int ans = 0;

  /* see if any source/description is missing; put knowns before unknowns */
  if ( ! file1->source || ! file1->description
      || ! file2->source || ! file2->description )
  {
    if ( file1->source && ! file2->source )
      return -1;
    if ( file2->source && ! file1->source )
      return 1;
    if ( file1->description && ! file2->description )
      return -1;
    if ( file2->description && ! file1->description )
      return 1;
  }

  /* sort by source or description first */
  if ( file1->source && file2->source
      && ( ans = strcmp( file1->source, file2->source ) ) )
    return ans;
  else if ( file1->description && file2->description
      && ( ans = strcmp( file1->description, file2->description ) ) )
    return ans;

  /* then sort by start time */
  if ( file1->startTime < file2->startTime )
    ans = -1;
  else if ( file1->startTime > file2->startTime )
    ans = 1;
  else /* finally, sort by url */
    if ( file1->url )
      if ( file2->url )
        ans = strcmp( file1->url, file2->url );
      else
        ans = -1;
    else if ( file2->url )
      ans = 1;
    else
      ans = 0;

  return ans;
}

/* <lalVerbatim file="FrameCacheCP"> */
void
LALFrCacheSieve(
    LALStatus     *status,
    FrCache      **output,
    FrCache       *input,
    FrCacheSieve  *params
    )
{ /* </lalVerbatim> */
  FrCache *cache;

  regex_t srcReg;
  regex_t dscReg;
  regex_t urlReg;

  UINT4 n;
  UINT4 i;

  INITSTATUS( status, "LALFrCacheSieve", FRAMECACHEC );
  ATTATCHSTATUSPTR( status );
  ASSERT( output, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );
  ASSERT( ! *output, status, FRAMECACHEH_ENNUL, FRAMECACHEH_MSGENNUL );
  ASSERT( input, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );
  ASSERT( params, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );

  if ( params->srcRegEx )
    regcomp( &srcReg, params->srcRegEx, REG_NOSUB );
  if ( params->dscRegEx )
    regcomp( &dscReg, params->dscRegEx, REG_NOSUB );
  if ( params->urlRegEx )
    regcomp( &urlReg, params->urlRegEx, REG_NOSUB );

  n = 0;
  for ( i = 0; i < input->numFrameFiles; ++i )
  {
    FrStat *file = input->frameFiles + i;
    if ( params->earliestTime > 0 )
      if ( params->earliestTime > file->startTime + file->duration )
        continue; /* file is too early */
    if ( params->latestTime > 0 )
      if ( file->startTime <= 0 || params->latestTime < file->startTime )
        continue; /* file is too late */
    if ( params->srcRegEx )
      if ( ! file->source || regexec( &srcReg, file->source, 0, NULL, 0 ) )
        continue; /* source doesn't match regex */
    if ( params->dscRegEx )
      if ( ! file->description
          || regexec( &dscReg, file->description, 0, NULL, 0 ) )
      continue; /* description doesn't match regex */
    if ( params->urlRegEx )
      if ( ! file->source || regexec( &urlReg, file->source, 0, NULL, 0 ) )
        continue; /* url doesn't match regex */
    ++n;
  }

  cache = *output = LALCalloc( 1, sizeof( **output ) );
  if ( ! cache )
  {
    ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
  }
  cache->numFrameFiles = n;
  if ( ! cache->numFrameFiles )
  {
    cache->frameFiles = NULL;
  }
  else
  {
    cache->frameFiles = LALCalloc( n, sizeof( *cache->frameFiles ) );
    if ( ! cache->frameFiles )
    {
      LALFree( *output );
      *output = NULL;
      ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
    }

    n = 0;
    for ( i = 0; i < input->numFrameFiles; ++i )
    {
      FrStat *file = input->frameFiles + i;
      if ( params->earliestTime > 0 )
        if ( params->earliestTime > file->startTime + file->duration )
          continue; /* file is too early */
      if ( params->latestTime > 0 )
        if ( file->startTime <= 0 || params->latestTime < file->startTime )
          continue; /* file is too late */
      if ( params->srcRegEx )
        if ( ! file->source || regexec( &srcReg, file->source, 0, NULL, 0 ) )
          continue; /* source doesn't match regex */
      if ( params->dscRegEx )
        if ( ! file->description
            || regexec( &dscReg, file->description, 0, NULL, 0 ) )
          continue; /* description doesn't match regex */
      if ( params->urlRegEx )
        if ( ! file->source || regexec( &urlReg, file->source, 0, NULL, 0 ) )
          continue; /* url doesn't match regex */

      /* copy frame file stat */
      cache->frameFiles[n] = *file;
      if ( file->source )
      {
        size_t size = strlen( file->source ) + 1;
        cache->frameFiles[n].source = LALMalloc( size );
        if ( ! cache->frameFiles[n].source )
        {
          TRY( LALDestroyFrCache( status->statusPtr, output ), status );
          ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
        }
        memcpy( cache->frameFiles[n].source, file->source, size );
      }
      if ( file->description )
      {
        size_t size = strlen( file->description ) + 1;
        cache->frameFiles[n].description = LALMalloc( size );
        if ( ! cache->frameFiles[n].description )
        {
          TRY( LALDestroyFrCache( status->statusPtr, output ), status );
          ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
        }
        memcpy( cache->frameFiles[n].description, file->description, size );
      }
      if ( file->url )
      {
        size_t size = strlen( file->url ) + 1;
        cache->frameFiles[n].url = LALMalloc( size );
        if ( ! cache->frameFiles[n].url )
        {
          TRY( LALDestroyFrCache( status->statusPtr, output ), status );
          ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
        }
        memcpy( cache->frameFiles[n].url, file->url, size );
      }
      ++n;
    }

    qsort( cache->frameFiles, cache->numFrameFiles, 
        sizeof( *cache->frameFiles ), FrStatCompare );
  }

  if ( params->srcRegEx ) regfree( &srcReg );
  if ( params->dscRegEx ) regfree( &dscReg );
  if ( params->urlRegEx ) regfree( &urlReg );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="FrameCacheCP"> */
void
LALFrCacheGenerate(
    LALStatus   *status,
    FrCache    **output,
    const CHAR  *dirstr,
    const CHAR  *fnptrn
    )
{ /* </lalVerbatim> */
  FrCache *cache;
  glob_t g;
  int globflags = 0;
  int i;

  INITSTATUS( status, "LALFrCacheGenerate", FRAMECACHEC );
  ATTATCHSTATUSPTR( status );
  ASSERT( output, status, FRAMECACHEH_ENULL, FRAMECACHEH_MSGENULL );
  ASSERT( ! *output, status, FRAMECACHEH_ENNUL, FRAMECACHEH_MSGENNUL );

  fnptrn = fnptrn ? fnptrn : "*.gwf";
  dirstr = dirstr ? dirstr : ".";

  if ( fnptrn[0] && 
      ( fnptrn[0] == '/'
        || ( fnptrn[0] == '.' && fnptrn[1] &&
          ( fnptrn[1] == '/' || ( fnptrn[1] == '.' && fnptrn[2] == '/' ) ) )
      )
     )
  {
    glob( fnptrn, globflags, NULL, &g );
  }
  else /* prepend path from dirname */
  {
    CHAR  path[MAXPATHLEN];
    CHAR  dirname[MAXPATHLEN];
    CHAR *nextdir;
    strncpy( dirname, dirstr,  sizeof( dirname ) - 1 );
    do
    {
      nextdir = strchr( dirname, ':' );
      if ( nextdir )
        *nextdir++ = 0;
      LALSnprintf( path, sizeof( path ) - 1, "%s/%s", 
          *dirname ? dirname : ".", fnptrn );
      glob( path, globflags, NULL, &g );
      globflags |= GLOB_APPEND;
    }
    while ( nextdir );
  }

  if ( ! g.gl_pathc )
  {
    globfree( &g );
    ABORT( status, FRAMECACHEH_EPATH, FRAMECACHEH_MSGEPATH );
  }

  cache = *output = LALCalloc( 1, sizeof( **output ) );
  if ( ! cache )
  {
    globfree( &g );
    ABORT( status, FRAMECACHEH_EPATH, FRAMECACHEH_MSGEPATH );
  }
  cache->numFrameFiles = g.gl_pathc;
  cache->frameFiles = LALCalloc( g.gl_pathc, sizeof( *cache->frameFiles ) );
  if ( ! cache->frameFiles )
  {
    globfree( &g );
    LALFree( *output );
    *output = NULL;
    ABORT( status, FRAMECACHEH_EPATH, FRAMECACHEH_MSGEPATH );
  }

  /* copy file names */
  for ( i = 0; i < g.gl_pathc; ++i )
  {
    FrStat *file = cache->frameFiles + i;
    char src[MAXPATHLEN];
    char dsc[MAXPATHLEN];
    int t0 = 0;
    int dt = 0;
    char *path = g.gl_pathv[i];
    char *base = strrchr( path, '/' );
    int c;
    if ( ! base )
      base = path;
    else
      ++base; /* one past the final '/' */
    if ( *path == '/' ) /* absolute path */
    {
      size_t urlsz = strlen( path ) + sizeof( "file://localhost" );
      file->url = LALMalloc( urlsz );
      if ( ! file->url )
      {
        globfree( &g );
        TRY( LALDestroyFrCache( status->statusPtr, output ), status );
        ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
      }
      LALSnprintf( file->url, urlsz, "file://localhost%s", path );
    }
    else /* relative path */
    {
      size_t urlsz = strlen( path ) + sizeof( "file://localhost" );
      CHAR cwd[MAXPATHLEN];
      getcwd( cwd, MAXPATHLEN - 1 );
      urlsz += strlen( cwd ) + 1;
      file->url = LALMalloc( urlsz );
      if ( ! file->url )
      {
        globfree( &g );
        TRY( LALDestroyFrCache( status->statusPtr, output ), status );
        ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
      }
      LALSnprintf( file->url, urlsz, "file://localhost%s/%s", cwd, path );
    }

    /* extract src, dsc, t0, and dt from file name */
    c = sscanf( base ? base : path, "%[a-zA-Z0-9_+#]-%[a-zA-Z0-9_+#]-%d-%d",
        src, dsc, &t0, &dt );
    if ( c == 4 ) /* expected format */
    {
      file->source = LALMalloc( strlen( src ) + 1 );
      if ( ! file->source )
      {
        globfree( &g );
        TRY( LALDestroyFrCache( status->statusPtr, output ), status );
        ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
      }
      strcpy( file->source, src );
      file->description = LALMalloc( strlen( dsc ) + 1 );
      if ( ! file->description )
      {
        globfree( &g );
        TRY( LALDestroyFrCache( status->statusPtr, output ), status );
        ABORT( status, FRAMECACHEH_EALOC, FRAMECACHEH_MSGEALOC );
      }
      strcpy( file->description, dsc );
      file->startTime = t0;
      file->duration = dt;
    }
  }

  globfree( &g );

  qsort( cache->frameFiles, cache->numFrameFiles, sizeof( *cache->frameFiles ),
      FrStatCompare );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
