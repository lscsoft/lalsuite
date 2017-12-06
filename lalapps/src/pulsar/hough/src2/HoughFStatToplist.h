/*
*  Copyright (C) 2007 Bernd Machenschalk, Reinhard Prix
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

#ifndef HOUGHFSTATTOPLIST_H
#define HOUGHFSTATTOPLIST_H

#include "HeapToplist.h"
#include <lal/LALDatatypes.h>


/** Type to hold the fields that will be kept in a "toplist"  */
typedef struct {
  REAL8 Freq;		/**< Frequency at maximum (?) of the cluster */
  REAL8 f1dot;		/**< spindown value f1dot = df/dt */
  REAL8 Alpha; 		/**< Skyposition: longitude in equatorial coords, radians */
  REAL8 Delta;		/**< skyposition: latitude */
  REAL8 HoughFStat;	/**< Hough significance */
  REAL8 AlphaBest;      /**< skyposition of best candidate: longitude */
  REAL8 DeltaBest;      /**< skyposition of best candidate: latitude */
  REAL8 MeanSig;        /**< mean of significance values in hough map*/
  REAL8 VarianceSig;    /**< variance of significance values in hough map*/
  REAL4 sumTwoF;        /**< sum of 2F-values as recomputed in LV postprocessing */
  REAL4Vector *sumTwoFX; /**< sum of 2F-values per detector, computed in LV postprocessing */
} HoughFStatOutputEntry;


/* This has by now been reduced to an interface to the HeapToplist functions */

/**
 * creates a toplist with length elements,
 * returns -1 on error (usually out of memory), else 0
 */
extern int create_houghFStat_toplist(toplist_t**list, UINT8 length);

/** frees the space occupied by the toplist */
extern void free_houghFStat_toplist(toplist_t**list);

/**
 * Inserts an element in to the toplist either if there is space left
 * or the element is larger than the smallest element in the toplist.
 * In the latter case, remove the smallest element from the toplist
 * Returns 1 if the element was actually inserted, 0 if not.
 */
extern int insert_into_houghFStat_toplist(toplist_t*list, HoughFStatOutputEntry line);

/**
 * Writes the toplist to an (already open) filepointer
 * Returns the number of written charactes
 * sets the checksum if non-NULL
 * Returns something <0 on error
 */
extern int write_houghFStat_toplist_to_fp(toplist_t*list, FILE*fp, UINT4*checksum);

/**
 * reads a (created!) toplist from an open filepointer
 * sets the checksum if non-NULL
 * reads maximum maxbytes, all that is there if maxbytes is 0
 * returns the number of bytes read,
 * 0 if we found a %DONE marker at the end,
 * -1 if the file contained a syntax error,
 * -2 if given an improper toplist
 */
extern int read_houghFStat_toplist_from_fp(toplist_t*list, FILE*fp, UINT4*checksum, UINT4 maxbytes);

/**
 * sorts the toplist with an internal sorting function,
 * used before finally writing it
 */
extern void sort_houghFStat_toplist(toplist_t*list);




/** File IO */

/**
 * writes an HoughFStatOutputEntry line to an open filepointer.
 * Returns the number of chars written, -1 if in error
 * Updates checksum if given (i.e. not NULL)
 */
extern int write_houghFStat_toplist_item_to_fp(HoughFStatOutputEntry line, FILE*fp, UINT4*checksum);

/**
 * writes the given toplitst to a temporary file, then renames the
 * temporary file to filename. The name of the temporary file is
 * derived from the filename by appending ".tmp". Returns the number
 * of chars written or -1 if the temp file could not be opened.
 */
extern int atomic_write_houghFStat_toplist_to_file(toplist_t*list, const char*filename, UINT4*checksum);

/**
 * meant for the final writing of the toplist
 * - reduces toplist precision
 * - sorts the toplist
 * - finally calls atomic_write_houghFStat_toplist_to_file()
 */
extern int final_write_houghFStat_toplist_to_file(toplist_t*list, const char*filename, UINT4*checksum);



/** new, simpler checkpointing for HierarchicalSearch */

/**
 * writes a checkpoint:
 * - constructs temporary filename (by appending .TMP)
 * - writes number of elements ("elems") in toplist to tempfile
 * - dumps data to tempfile
 * - appends counter
 * - appends checksum (of elems, data and counter)
 * - renames tempfile to final name
 * returns
 * -1 in case of an I/O error,
 * -2 if out of memory,
 * 0 otherwise (successful)
 */
extern int write_hfs_checkpoint(const char*filename, toplist_t*tl, UINT4 counter, BOOLEAN do_sync);

/**
 * tries to read a checkpoint
 * - tries to open the file, returns 1 if no file found
 * - reads elems, data, counter and checksum
 * - verifies checksum
 * - restores the heap by sorting
 * returns
 * 0 if successfully read a checkpoint
 * 1 if no checkpoint was found
 * -1 in case of an I/O error
 * -2 if the checksum was wrong or elems was unreasonable
 */
extern int read_hfs_checkpoint(const char*filename, toplist_t*tl, UINT4*counter);

/**
 * write the final output file:
 * - re-sort the toplist into freq/alpha/delta/fdot order
 * - write out the toplist in ASCII format with end marker to a temporary file
 * - rename the file to the final name
 */
extern int write_hfs_oputput(const char*filename, toplist_t*tl);

#endif /* HOUGHFSTATTOPLIST_H - double inclusion protection */
