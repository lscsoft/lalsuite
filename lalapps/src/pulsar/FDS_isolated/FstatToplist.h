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

#ifndef FSTATTOPLIST_H
#define FSTATTOPLIST_H

#include "ComputeFStatistic.h"
#include "HeapToplist.h"



/* This has by now been reduced to an interface to the HeapToplist functions */

/**
 * creates a toplist with length elements,
 * returns -1 on error (usually out of memory), else 0
 */
extern int create_fstat_toplist(toplist_t**list, UINT8 length);

/** frees the space occupied by the toplist */
extern void free_fstat_toplist(toplist_t**list);

/**
 * Inserts an element in to the toplist either if there is space left
 * or the element is larger than the smallest element in the toplist.
 * In the latter case, remove the smallest element from the toplist
 * Returns 1 if the element was actually inserted, 0 if not.
 */
extern int insert_into_fstat_toplist(toplist_t*list, FstatOutputEntry line);

/**
 * Writes the toplist to an (already open) filepointer
 * Returns the number of written charactes
 * sets the checksum if non-NULL
 * Returns something <0 on error
 */
extern int write_fstat_toplist_to_fp(toplist_t*list, FILE*fp, UINT4*checksum);

/**
 * reads a (created!) toplist from an open filepointer
 * sets the checksum if non-NULL
 * reads maximum maxbytes, all that is there if maxbytes is 0
 * returns the number of bytes read,
 * 0 if we found a %DONE marker at the end,
 * -1 if the file contained a syntax error,
 * -2 if given an improper toplist
 */
extern int read_fstat_toplist_from_fp(toplist_t*list, FILE*fp, UINT4*checksum, UINT4 maxbytes);

/**
 * sorts the toplist with an internal sorting function,
 * used before finally writing it
 */
extern void sort_fstat_toplist(toplist_t*list);




/** File IO */

/**
 * writes an FstatOutputEntry line to an open filepointer.
 * Returns the number of chars written, -1 if in error
 * Updates checksum if given (i.e. not NULL)
 */
extern int write_fstat_toplist_item_to_fp(FstatOutputEntry line, FILE*fp, UINT4*checksum);

/**
 * writes the given toplitst to a temporary file, then renames the
 * temporary file to filename. The name of the temporary file is
 * derived from the filename by appending ".tmp". Returns the number
 * of chars written or -1 if the temp file could not be opened.
 */
extern int atomic_write_fstat_toplist_to_file(toplist_t*list, const char*filename, UINT4*checksum);

/**
 * meant for the final writing of the toplist
 * - reduces toplist precision
 * - sorts the toplist
 * - finally calls atomic_write_fstat_toplist_to_file()
 */
extern int final_write_fstat_toplist_to_file(toplist_t*list, const char*filename, UINT4*checksum);




/** a toplist as a checkpointed file */

typedef struct {
  CHAR* filename;    /**< name of the toplist file */
  CHAR* buffer;      /**< write buffer if needed */
  UINT4 bufsize;     /**< buffer size if needed */
  UINT4 bytes;       /**< counts the bytes in the file */
  UINT4 maxsize;     /**< the file must not grow larger than that */
  UINT4 checksum;    /**< keeps the checksum */
  FILE* fp;          /**< FILE* currently associated */
  toplist_t*list;    /**< toplist this file reflects */
} FStatCheckpointFile;

/** creates a FStatCheckpointFile */
extern int fstat_cpt_file_create (FStatCheckpointFile **cptf,
				  CHAR  *filename,
				  UINT4 bufsize,
				  UINT4 maxsize,
				  toplist_t*tl);

/** destroys a FStatCheckpointFile */
extern int fstat_cpt_file_destroy (FStatCheckpointFile **cptf);

/** opens a file for checkpointing the desired toplist */
extern int fstat_cpt_file_open (FStatCheckpointFile *cptf);

/** flushes the checkpoint file (only useful if buffered) */
extern int fstat_cpt_file_flush (FStatCheckpointFile *cptf);

/** returns information for checkpointing */
extern int fstat_cpt_file_info (FStatCheckpointFile *cptf, CHAR**filename, UINT4*bytes, UINT4*checksum);

/**
 * adds an item to the toplist and keeps the file consistent, i.e.
 * adds the entry to the file if it was really inserted
 * and compacts the file if necessary
 */
extern int fstat_cpt_file_add  (FStatCheckpointFile*cptf, FstatOutputEntry line);

/**
 * closes the file, reduces the precision, sorts the toplist,
 * finally rewrites the file (sorted and compact) with end marker
 */
extern int fstat_cpt_file_close(FStatCheckpointFile*cptf);

/** reads a written checkpointed toplist back into memory */
extern int fstat_cpt_file_read (FStatCheckpointFile*cptf, UINT4 checksum, UINT4 maxbytes);

/** compact a toplist file if the length has reached maxbytes */
extern int fstat_cpt_file_compact(FStatCheckpointFile*cptf);

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
extern int write_hs_checkpoint(const char*filename, toplist_t*tl, UINT4 counter, BOOLEAN do_sync);

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
extern int read_hs_checkpoint(const char*filename, toplist_t*tl, UINT4*counter);

/**
 * write the final output file:
 * - re-sort the toplist into freq/alpha/delta/fdot order
 * - write out the toplist in ASCII format with end marker to a temporary file
 * - rename the file to the final name
 */
extern int write_hs_oputput(const char*filename, toplist_t*tl);

#endif /* FSTATTOPLIST_H - double inclusion protection */
