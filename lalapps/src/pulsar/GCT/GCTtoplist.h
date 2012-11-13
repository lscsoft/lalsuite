/*
*  Copyright (C) 2007, 2009 Bernd Machenschalk, Reinhard Prix, Holger Pletsch
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

#ifndef GCTFSTATTOPLIST_H
#define GCTFSTATTOPLIST_H

#include "HeapToplist.h"
#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif
extern char**global_argv;
extern int global_argc;
#ifdef __cplusplus
}
#endif

#ifndef GCTTOP_MAX_IFOS
#define GCTTOP_MAX_IFOS 10
#endif
/** Type to hold the fields that will be kept in a "toplist"  */
typedef struct {
  REAL8 Freq;  /**< frequency */
  REAL8 F1dot;/**< spindown value f1dot = df/dt */
  REAL8 F2dot;/**< spindown value f2dot = d2f/dt2 */
  REAL8 Alpha; /**< skyposition: longitude in equatorial coords, radians */
  REAL8 Delta;/**< skyposition: latitude */
  REAL4 sumTwoF;  /**< sum of 2F-values */
  UINT4 nc;       /**< number count */
  REAL4 LV;       /**< Line Veto statistic */
  UINT4 numDetectors; /**< number of detectors for optional sumTwoFX arrays */
  REAL4 sumTwoFX[GCTTOP_MAX_IFOS]; /**< fixed-size array of single-detector 2F-values */
  REAL4 sumTwoFrecalc;  /**< sum of 2F-values as recomputed by recalcToplistStats */
  REAL4 sumTwoFXrecalc[GCTTOP_MAX_IFOS];  /**< fixed-size array of single-detector 2F-values as recomputed by recalcToplistStats */
} GCTtopOutputEntry;

/// enumerate all toplist-sorting options: by F (0), number-count (1), LV-stat (2), "dual" toplists F + LV (3)
typedef enum
  {
    SORTBY_F 		= 0,	//< sort by multi-IFO F-stat (averaged over segments)
    SORTBY_NC 		= 1,	//< sort by number-count 'nc'
    SORTBY_LV 		= 2,	//< sort by line-veto statistic 'LV'
    SORTBY_DUAL_F_LV 	= 3,	//< dual toplists: one sorted by F, one by LV
    SORTBY_LAST			//< end-marker
  } SortBy_t;

/* This has by now been reduced to an interface to the HeapToplist functions */

/** creates a toplist with length elements,
   returns -1 on error (usually out of memory), else 0 */
extern int create_gctFStat_toplist(toplist_t**list, UINT8 length, SortBy_t whatToSortBy);

/** frees the space occupied by the toplist */
extern void free_gctFStat_toplist(toplist_t**list);

/** Inserts an element in to the toplist either if there is space left
   or the element is larger than the smallest element in the toplist.
   In the latter case, remove the smallest element from the toplist
   Returns 1 if the element was actually inserted, 0 if not. */
extern int insert_into_gctFStat_toplist(toplist_t*list, GCTtopOutputEntry line);


/** Writes the toplist to an (already open) filepointer
   Returns the number of written charactes
   sets the checksum if non-NULL
   Returns something <0 on error */
extern int write_gctFStat_toplist_to_fp(toplist_t*list, FILE*fp, UINT4*checksum);


/** sorts the toplist with an internal sorting function,
   used before finally writing it */
extern void sort_gctFStat_toplist(toplist_t*list);


/** sorts the toplist with an internal sorting function,
 used before doing the follow-up analysis */
extern void sort_gctFStat_toplist_strongest(toplist_t*list);



/** File IO */

/** new, simpler checkpointing for HierarchicalSearch */

/** writes a checkpoint:
    - constructs temporary filename (by appending .TMP)
    - writes number of elements ("elems") in toplist to tempfile
    - dumps data to tempfile
    - appends counter
    - appends checksum (of elems, data and counter)
    - renames tempfile to final name
    returns
    -1 in case of an I/O error,
    -2 if out of memory,
     0 otherwise (successful)
*/
extern int write_hfs_checkpoint(const char*filename, toplist_t*tl, UINT4 counter, BOOLEAN do_sync);

/** tries to read a checkpoint
    - tries to open the file, returns 1 if no file found
    - reads elems, data, counter and checksum
    - verifies checksum
    - restores the heap by sorting
    returns
     0 if successfully read a checkpoint
     1 if no checkpoint was found
    -1 in case of an I/O error
    -2 if the checksum was wrong or elems was unreasonable
*/
extern int read_hfs_checkpoint(const char*filename, toplist_t*tl, UINT4*counter);

/** write the final output file:
    - re-sort the toplist into freq/alpha/delta/fdot order
    - write out the toplist in ASCII format with end marker to a temporary file
    - rename the file to the final name
*/
extern int write_hfs_oputput(const char*filename, toplist_t*tl);

#endif /* GCTFSTATTOPLIST_H - double inclusion protection */
