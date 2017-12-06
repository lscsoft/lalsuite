/*
*  Copyright (C) 2007 Jolien Creighton
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

/**
 * \file
 * \ingroup support
 * \author Creighton, J. D. E.
 * \date 2007
 * \brief This header covers routines to create and manipulate LALCache
 * structures and to read LAL cache files.
 *
 */
#ifndef _LALCACHE_H_
#define _LALCACHE_H_

#include <stdio.h>
#include <lal/LALDatatypes.h>
#include <lal/FileIO.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}       /* to match preceeding brace */
#endif
struct tagLALCacheEntry;
struct tagLALCache;

/** An entry in a LAL cache. */
typedef struct tagLALCacheEntry {
    CHAR *src;          /**< File source field */
    CHAR *dsc;          /**< File description field */
    INT4 t0;            /**< GPS time (seconds) of beginning of data in file */
    INT4 dt;            /**< Duration (seconds) of data in file */
    CHAR *url;          /**< URL of file */
} LALCacheEntry;

/** The LALCache structure is an array of entries. */
typedef struct tagLALCache {
    UINT4 length;
    LALCacheEntry *list;
} LALCache;

/** Creates a LALCache structure. */
LALCache *XLALCreateCache(UINT4 length);

/** Destroys a LALCache structure. */
void XLALDestroyCache(LALCache * cache);

/** Duplicates a LALCache structure. */
LALCache *XLALCacheDuplicate(const LALCache * cache);

/** Returns a new LALCache structure that is the merge of two. */
LALCache *XLALCacheMerge(const LALCache * cache1, const LALCache * cache2);

/** Reads a LAL cache file and produces a LALCache structure. */
LALCache *XLALCacheFileRead(LALFILE * fp);

/** Reads a LAL cache file and produces a LALCache structure. */
LALCache *XLALCacheImport(const char *fname);

/**
 * Globs a directory and construct LALCache from matching entries.
 * \param [in] dirstr Colon-delimited list of directories.
 * \param [in] fnptrn Glob pattern for matching files.
 * \returns LALCache structure.
 */
LALCache *XLALCacheGlob(const char *dirstr, const char *fnptrn);

/** Writes a LALCache structure to output LALFILE. */
int XLALCacheFileWrite(LALFILE * fp, const LALCache * cache);

/** Exports a LALCache structure to an output LAL cache file. */
int XLALCacheExport(const LALCache * cache, const char *filename);

/** Sorts entries in a LALCache structure. */
int XLALCacheSort(LALCache * cache);

/**
 * Prunes duplicate entries keeping the second one; cache is reduced in
 * length if there are.  Entries are duplicates if their metadata are
 * the same (even if the urls are different
 */
int XLALCacheUniq(LALCache * cache);

/**
 * Selects only matching entries in a LALCache structure -- other entries
 * are deleted from the LALCache structure.
 * \param cache *UNDOCUMENTED*
 * \param t0 Remove entries ending before t0 (0 to disable).
 * \param t1 Remove entries ending after t1 (0 to disable).
 * \param srcregex Regular expression to match src field (NULL to disable).
 * \param dscregex Regular expression to match dsc field (NULL to disable).
 * \param urlregex Regular expression to match url field (NULL to disable).
 */
int XLALCacheSieve(LALCache * cache, INT4 t0, INT4 t1,
                   const char *srcregex, const char *dscregex,
                   const char *urlregex);

/**
 * Finds the first entry that contains the requested time, or the first entry
 * after the time if the time is in a gap or before the first entry.  Returns
 * NULL if the time is after the last entry.
 */
LALCacheEntry *XLALCacheEntrySeek(const LALCache * cache, double t);


/** Open a file identified by an entry in a LALCache structure. */
LALFILE *XLALCacheEntryOpen(const LALCacheEntry * entry);

#if 0
{       /* to match succeding brace */
#endif
#if defined(__cplusplus)
}
#endif
#endif /* _LALCACHE_H_ */
