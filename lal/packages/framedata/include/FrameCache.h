/**** <lalVerbatim file="FrameCacheHV">
 * Author: Jolien D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{FrameCache.h}}
 *
 * Routines for manipulating a cache of available frame files.
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/FrameCache.h>
 * \end{verbatim}
 *
 * Frame file catalogues contain URLs for available frame files along with
 * various pieces of metadata.  These routines allow frame catalogues to be
 * imported, exported, or generated, and stored in a frame cache structure,
 * which can be manipulated.
 *
 **** </lalLaTeX> */

#ifndef _FRAMECACHE_H
#define _FRAMECACHE_H

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

#include <lal/LALDatatypes.h>

NRCSID( FRAMECACHEH, "$Id$" );

/**** <lalLaTeX>
 *
 * \subsection*{Error conditions}
 *
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define FRAMECACHEH_ENULL  1
#define FRAMECACHEH_ENNUL  2
#define FRAMECACHEH_EALOC  4
#define FRAMECACHEH_EIEIO  8
#define FRAMECACHEH_ELINE 16
#define FRAMECACHEH_EPATH 32
#define FRAMECACHEH_ENFRM 64

#define FRAMECACHEH_MSGENULL "Null pointer"
#define FRAMECACHEH_MSGENNUL "Non-null pointer"
#define FRAMECACHEH_MSGEALOC "Memory allocation error"
#define FRAMECACHEH_MSGEIEIO "Import/export I/O error"
#define FRAMECACHEH_MSGELINE "Input line too long"
#define FRAMECACHEH_MSGEPATH "Unable to glob frame files to build cache"
#define FRAMECACHEH_MSGENFRM "No frame files"
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 * \idx[Type]{FrStat}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagFrStat
{
  CHAR *source;
  CHAR *description;
  INT4  startTime;
  INT4  duration;
  CHAR *url;
}
FrStat;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure contains a frame file status.  The fields are:
 * \begin{description}
 * \item[\texttt{source}] the source detector(s) of the data in the frame
 *     file, or other identifier.
 * \item[\texttt{description}] the description of the type of data contained
 *     in the frame file, or other identifier.
 * \item[\texttt{startTime}] the GPS time of the second equal to
 *     (or just before) the start of the data contained in the frame file.
 * \item[\texttt{duration}] the number of seconds between \texttt{startTime}
 *     and the GPS time of the second equal to (or just after) the end of the
 *     data contained in the frame file.
 * \item[\texttt{url}] the URL of the frame file.
 * \end{description}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagFrCache
{
  UINT4   numFrameFiles;
  FrStat *frameFiles;
}
FrCache;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure contains a list of all frame files available.  The fields are:
 * \begin{description}
 * \item[\texttt{numFrameFiles}] the total number of frame files in the list.
 * \item[\texttt{frameFiles}] array of frame file status descriptors.
 * \end{description}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagFrCacheSieve
{
  const CHAR *srcRegEx;
  const CHAR *dscRegEx;
  const CHAR *urlRegEx;
  INT4 earliestTime;
  INT4 latestTime;
}
FrCacheSieve;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure contains parameters to use to extract those frame files
 * of interest from a cache.  The parameters include regular expressions and
 * time ranges.  The fields are:
 * \begin{description}
 * \item[\texttt{srcRegEx}] regular expression to use in selecting frame files
 *     with a specified source identifier.  (Not used if \texttt{NULL}.)
 * \item[\texttt{dscRegEx}] regular expression to use in selecting frame files
 *     with a specified description identifier.  (Not used if \texttt{NULL}.)
 * \item[\texttt{urlRegEx}] regular expression to use in selecting frame files
 *     with a specified URL.  (Not used if \texttt{NULL}.)
 * \item[\texttt{earliestTime}] earliest time (GPS seconds) of frame files of
 *     interest.  (Not used if zero or less.)
 * \item[\texttt{latestTime}] latest time (GPS seconds) of frame files of
 *     interest.  (Not used if zero or less.)
 * \end{description}
 *
 * \vfill{\footnotesize\input{FrameCacheHV}}
 * \newpage\input{FrameCacheC}
 *
 **** </lalLaTeX> */

FrCache * XLALFrImportCache( const char *fname );
FrCache * XLALFrSieveCache( FrCache *input, FrCacheSieve *params );
FrCache * XLALFrGenerateCache( const CHAR *dirstr, const CHAR *fnptrn );
int XLALFrExportCache( FrCache *cache, const CHAR *fname );
void XLALFrDestroyCache( FrCache *cache );

void LALFrCacheImport(
    LALStatus   *status,
    FrCache    **output,
    const CHAR  *fname
    );

void LALFrCacheExport(
    LALStatus  *status,
    FrCache    *cache,
    const CHAR *fname
    );

void LALDestroyFrCache(
    LALStatus  *status,
    FrCache   **cache
    );

void
LALFrCacheSieve(
    LALStatus     *status,
    FrCache      **output,
    FrCache       *input,
    FrCacheSieve  *params
    );

void
LALFrCacheGenerate(
    LALStatus   *status,
    FrCache    **output,
    const CHAR  *dirstr,
    const CHAR  *fnptrn
    );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _FRAMECACHE_H */
