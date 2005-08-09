#ifndef FSTATTOPLIST_H
#define FSTATTOPLIST_H

#include "ComputeFStatistic.h"

/* toplist structure based on FstatsClusterOutput */
typedef struct {
    UINT8 length;   /* the length (maximal number of entries) of the toplist */
    UINT8 elems;    /* number of elements currently in the toplist */
    UINT8 smallest; /* index of the smallest element in the toplist */
    FstatsClusterOutput *data; /* points to the actual data */
} toplist;

/* creates a toplist with length elements,
   returns -1 on error (usually out of memory), else 0 */
extern int create_toplist(toplist**, UINT8);

/* frees the space occupied by the toplist */
extern void free_toplist(toplist**);

/* Inserts an element in to the toplist either if there is space left
   or the element is larger than the smallest element in the toplist.
   In the latter case, remove the smallest element from the toplist and
   look for the now smallest one.
   Returns 1 if the element was actually inserted, 0 if not. */
extern int insert_into_toplist(toplist*, FstatsClusterOutput);

/* Writes the toplist to an (already open) filepointer
   Returns the number of written charactes
   Returns something <0 on error */
extern int write_toplist_to_fp(toplist*, FILE*);

/* reads a (created!) toplist from an open filepointer
   returns -1 if the file contained a syntax error, -2 if given an improper toplist */
extern int read_toplist_from_fp(toplist*, FILE*);

#endif /* FSTATTOPLIST_H */

