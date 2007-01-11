#ifndef FSTATTOPLIST_H
#define FSTATTOPLIST_H

#include "ComputeFStatistic.h"
#include "HeapToplist.h"



/* This has by now been reduced to an interface to the HeapToplist functions */

/* creates a toplist with length elements,
   returns -1 on error (usually out of memory), else 0 */
extern int create_fstat_toplist(toplist_t**list, UINT8 length);

/* frees the space occupied by the toplist */
extern void free_fstat_toplist(toplist_t**list);

/* Inserts an element in to the toplist either if there is space left
   or the element is larger than the smallest element in the toplist.
   In the latter case, remove the smallest element from the toplist
   Returns 1 if the element was actually inserted, 0 if not. */
extern int insert_into_fstat_toplist(toplist_t*list, FstatOutputEntry line);

/* Writes the toplist to an (already open) filepointer
   Returns the number of written charactes
   sets the checksum if non-NULL
   Returns something <0 on error */
extern int write_fstat_toplist_to_fp(toplist_t*list, FILE*fp, UINT4*checksum);

/* reads a (created!) toplist from an open filepointer
   sets the checksum if non-NULL
   reads maximum maxbytes, all that is there if maxbytes is 0
   returns the number of bytes read,
   -1 if the file contained a syntax error,
   -2 if given an improper toplist */
extern int read_fstat_toplist_from_fp(toplist_t*list, FILE*fp, UINT4*checksum, UINT4 maxbytes);

/* sorts the toplist with an internal sorting function,
   used before finally writing it */
extern void sort_fstat_toplist(toplist_t*list);




/* File IO */

/* writes an FstatOutputEntry line to an open filepointer.
   Returns the number of chars written, -1 if in error
   Updates checksum if given (i.e. not NULL) */
extern int write_fstat_toplist_item_to_fp(FstatOutputEntry line, FILE*fp, UINT4*checksum);

/* writes the given toplitst to a temporary file, then renames the
   temporary file to filename. The name of the temporary file is
   derived from the filename by appending ".tmp". Returns the number
   of chars written or -1 if the temp file could not be opened. */
extern int atomic_write_fstat_toplist_to_file(toplist_t*list, char*filename, UINT4*checksum);

/* meant for the final writing of the toplist
   - reduces toplist precision
   - sorts the toplist
   - the calls atomic_write_fstat_toplist_to_file() */
extern int final_write_fstat_toplist_to_file(toplist_t*list, char*filename, UINT4*checksum);




/* next level: a toplist as a checkpointed file */

typedef struct {
  CHAR* filename;    /* name of the toplist file */
  CHAR* buffer;      /* write buffer if needed */
  UINT4 bufsize;     /* buffer size if needed */
  UINT4 bytes;       /* counts the bytes in the file */
  UINT4 maxsize;     /* the file must not grow larger than that */
  UINT4 checksum;    /* keeps the checksum */
  FILE* fp;          /* FILE* currently associated */
  toplist_t*list;    /* toplist this file reflects */
} FStatCheckpointFile;


extern int fstat_cpt_file_create (FStatCheckpointFile **cptf,
				  CHAR  *filename,
				  UINT4 bufsize,
				  UINT4 maxsize,
				  toplist_t*tl);

extern int fstat_cpt_file_destroy (FStatCheckpointFile **cptf);

/* opens a file for checkpointing the desired toplist */
extern int fstat_cpt_file_open (FStatCheckpointFile *cptf);

/* flushes the checkpoint file (only useful if buffered) */
extern int fstat_cpt_file_flush (FStatCheckpointFile *cptf);

/* returns information for checkpointing */
extern int fstat_cpt_file_info (FStatCheckpointFile *cptf, CHAR**filename, UINT4*bytes, UINT4*checksum);

/* adds an item to the toplist and keeps the file consistent, i.e.
   adds the entry to the file if it was really inserted
   and compacts the file if necessary */
extern int fstat_cpt_file_add  (FStatCheckpointFile*cptf, FstatOutputEntry line);

/* closes and compacts the file */
extern int fstat_cpt_file_close(FStatCheckpointFile*cptf);

/* reads a written checkpointed toplist back into memory */
extern int fstat_cpt_file_read (FStatCheckpointFile*cptf, UINT4 checksum, UINT4 maxbytes);

#endif /* FSTATTOPLIST_H */
