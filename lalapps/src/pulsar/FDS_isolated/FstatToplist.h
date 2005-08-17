#ifndef FSTATTOPLIST_H
#define FSTATTOPLIST_H

#include "ComputeFStatistic.h"

/* toplist structure based on FstatsClusterOutput */
typedef struct {
    UINT8 length;   /* the length (maximal number of entries) of the toplist */
    UINT8 elems;    /* number of elements currently in the toplist */
    UINT8 smallest; /* index of the smallest element in the toplist */
    FstatsClusterOutput *data;    /* points to the actual data */
    FstatsClusterOutput **sorted; /* an array of sorted pointers to data */
} toplist_t;

/* creates a toplist with length elements,
   returns -1 on error (usually out of memory), else 0 */
extern int create_toplist(toplist_t**, UINT8);

/* frees the space occupied by the toplist */
extern void free_toplist(toplist_t**);

/* Inserts an element in to the toplist either if there is space left
   or the element is larger than the smallest element in the toplist.
   In the latter case, remove the smallest element from the toplist and
   look for the now smallest one.
   Returns 1 if the element was actually inserted, 0 if not. */
extern int insert_into_toplist(toplist_t*, FstatsClusterOutput);

/* Writes the toplist to an (already open) filepointer
   Returns the number of written charactes
   sets the checksum if non-NULL
   Returns something <0 on error */
extern int write_toplist_to_fp(toplist_t*, FILE*, UINT4*);

/* reads a (created!) toplist from an open filepointer
   sets the checksum if non-NULL
   returns -1 if the file contained a syntax error,
   -2 if given an improper toplist */
extern int read_toplist_from_fp(toplist_t*, FILE*, UINT4*);

/* sorts the toplist with an internal sorting function,
   used before finally writing it */
extern void sort_toplist(toplist_t*);

/* writes an FstatsClusterOutput line to an open filepointer.
   Returns the number of chars written, -1 if in error
   Updates checksum if given */
extern int write_toplist_item_to_fp(FstatsClusterOutput, FILE*, UINT4*);

/* writes the given toplitst to a temporary file, then renames the
   temporary file to filename. The name of the temporary file is
   derived from the filename by appending ".tmp". Returns the number
   of chars written or -1 if the temp file could not be opened. */
extern int atomic_write_toplist_to_file(toplist_t*, char*, UINT4*);

#endif /* FSTATTOPLIST_H */
