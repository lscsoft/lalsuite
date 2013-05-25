/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton
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
#ifndef _FRAMECACHE_H
#define _FRAMECACHE_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#include <lal/LALDatatypes.h>

/**
 * \defgroup FrameCache_h Header FrameCache.h
 * \ingroup pkg_framedata
 *
 * \author Jolien D. E. Creighton
 *
 * Routines for manipulating a cache of available frame files.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/FrameCache.h>
 * \endcode
 *
 * Frame file catalogues contain URLs for available frame files along with
 * various pieces of metadata.  These routines allow frame catalogues to be
 * imported, exported, or generated, and stored in a frame cache structure,
 * which can be manipulated.
 *
*/
/*@{*/

/**\name Error Codes */
/*@{*/
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
/*@}*/

/**
 *
 * This structure contains a frame file status.  The fields are:
 * <dl>
 * <dt>source</dt><dd> the source detector(s) of the data in the frame
 *     file, or other identifier.
 * </dd><dt>description</dt><dd> the description of the type of data contained
 *     in the frame file, or other identifier.
 * </dd><dt>startTime</dt><dd> the GPS time of the second equal to
 *     (or just before) the start of the data contained in the frame file.
 * </dd><dt>duration</dt><dd> the number of seconds between \c startTime
 *     and the GPS time of the second equal to (or just after) the end of the
 *     data contained in the frame file.
 * </dd><dt>url</dt><dd> the URL of the frame file.
 * </dd></dl>
 *
*/
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

/**
 *
 * This structure contains a list of all frame files available.  The fields are:
 * <dl>
 * <dt>numFrameFiles</dt><dd> the total number of frame files in the list.
 * </dd><dt>frameFiles</dt><dd> array of frame file status descriptors.
 * </dd></dl>
 *
*/
typedef struct
tagFrCache
{
  UINT4   numFrameFiles;
  FrStat *frameFiles;
}
FrCache;

/**
 *
 * This structure contains parameters to use to extract those frame files
 * of interest from a cache.  The parameters include regular expressions and
 * time ranges.  The fields are:
 * <dl>
 * <dt>srcRegEx</dt><dd> regular expression to use in selecting frame files
 *     with a specified source identifier.  (Not used if \c NULL.) </dd>
 * <dt>dscRegEx</dt><dd> regular expression to use in selecting frame files
 *     with a specified description identifier.  (Not used if \c NULL.)</dd>
 * <dt>urlRegEx</dt><dd> regular expression to use in selecting frame files
 *     with a specified URL.  (Not used if \c NULL.)</dd>
 * <dt>earliestTime</dt><dd> earliest time (GPS seconds) of frame files of
 *     interest.  (Not used if zero or less.)</dd>
 * <dt>latestTime</dt><dd> latest time (GPS seconds) of frame files of
 *     interest.  (Not used if zero or less.)</dd>
 * </dl>
 *
 *
 *
 *
*/
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(STRUCT_IMMUTABLE(tagFrCacheSieve, srcRegEx, dscRegEx, urlRegEx));
#endif /* SWIG */
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
LALFrSieveCache(
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

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _FRAMECACHE_H */
