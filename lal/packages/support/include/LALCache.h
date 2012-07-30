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

/** \file
 * \ingroup support
 * \author Creighton, J. D. E.
 * \date 2007
 * \brief This header covers routines to create and manipulate LALCache
 * structures and to read LAL cache files.
 *
 */
#ifndef _LALCACHE_H_
#define _LALCACHE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <lal/LALDatatypes.h>
#include <lal/FileIO.h>

/** An entry in a LAL cache */
typedef struct tagLALCacheEntry {
        CHAR *src; /**< File source field */
        CHAR *dsc; /**< File description field */
        INT4 t0;   /**< GPS time (seconds) of beginning of data in file */
        INT4 dt;   /**< Duration (seconds) of data in file */
        CHAR *url; /**< URL of file */
} LALCacheEntry;

/** The LALCache structure is an array of entries */
typedef struct tagLALCache {
        UINT4 length;
        LALCacheEntry *list;
} LALCache;

/** Creates a LALCache structure */
LALCache * XLALCreateCache( UINT4 length );

/** Destroys a LALCache structure */
void XLALDestroyCache( LALCache *cache );

/** Duplicates a LALCache structure */
LALCache * XLALCacheDuplicate( LALCache *cache );

/** Returns a new LALCache structure that is the merge of two */
LALCache * XLALCacheMerge( LALCache *cache1, LALCache *cache2 );

/** Reads a LAL cache file and produces a LALCache structure */
LALCache * XLALCacheFileRead( LALFILE *fp );

/** Globs a directory and construct LALCache from matching entries */
LALCache * XLALCacheGlob(
                const char *dirstr, /**< colon-delimited list of directories */
                const char *fnptrn  /**< glob pattern for matching files */
                );

/** Writes a LALCache structure to an output LAL cache file */
int XLALCacheFileWrite( LALFILE *fp, LALCache *cache );

/** Sorts entries in a LALCache structure */
int XLALCacheSort( LALCache *cache );

/** Prunes duplicate entries keeping the second one; cache is reduced in
 * length if there are.  Entries are duplicates if their metadata are
 * the same (even if the urls are different */
int XLALCacheUniq( LALCache *cache );

/** Selects only matching entries in a LALCache structure -- other entries
 * are deleted from the LALCache structure */
int XLALCacheSieve(
                LALCache *cache, /**< The LALCache structure - modified */
                INT4 t0, /**< Remove entries ending before t0 (0 to disable) */
                INT4 t1, /**< Remove entries ending after t1 (0 to disable) */
                const char *srcregex, /**< Regular expression to match src field (NULL to disable) */
                const char *dscregex, /**< Regular expression to match dsc field (NULL to disable) */
                const char *urlregex /**< Regular expression to match url field (NULL to disable) */
                );

/** Open a file identified by an entry in a LALCache structure */
LALFILE * XLALCacheEntryOpen( LALCacheEntry *entry );

#ifdef __cplusplus
}
#endif

#endif /* _LALCACHE_H_ */
