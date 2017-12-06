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

/* This module keeps a "toplist", i.e. a list of the top n elements (accoding to an externally
   supplied comparison function) in a standard heap structure.

   Author (of this implementation): Bernd Machenschalk
 */

#ifndef HEAPTOPLIST_H
#define HEAPTOPLIST_H

#include <stddef.h>

/* toplist structure */
typedef struct {
  size_t length;   /* the length (maximal number of entries) of the toplist */
  size_t elems;    /* number of elements currently in the toplist */
  size_t size;     /* size of an element */
  char   *data;    /* the actual data array of 'length'*'size' chars */
  char   **heap;   /* array of 'length' pointers into data */
  int    (*smaller)(const void *, const void *); /* comparison function */
} toplist_t;


/* creates a toplist with 'length' elements of size 'size', with
   odering based on comparison function 'smaller'.
   returns -1 on error (out of memory), else 0 */
extern int create_toplist(toplist_t**list, size_t length, size_t size,
			  int (*smaller)(const void *, const void *));


/* frees the space occupied by the toplist */
extern void free_toplist(toplist_t**list);


/* Inserts an element in to the toplist either if there is space left
   or the element is larger than the smallest element in the toplist.
   In the latter case, remove the smallest element from the toplist.
   Returns 1 if the element was actually inserted, 0 if not. */
extern int insert_into_toplist(toplist_t*list, void *element);


/* clears an existing toplist of all elements inserted so far */
extern void clear_toplist(toplist_t*list);


/* apply a function to all elements of the list in the current order
   (possibly after calling qsort_toplist(), e.g. for writing out */
extern void go_through_toplist(toplist_t*list, void (*handle)(void *));


/* sort the toplist with an arbitrary sorting function
   (potentially) destroying the heap property.

   note that a (q-)sorted list is a heap, but due to the interface of qsort
   the same comparison function will give the reverse order than the heap.
   in order to restore a heap with qsort_toplist() (e.g. to add more elements) you must
   qsort_toplist() with the inverse function of the "smaller" function of the heap. */
extern void qsort_toplist(toplist_t*list, int (*compare)(const void *, const void *));


/* therefore we provide a qsort_toplist_r() function that gives the reverse ordering of
   qsort_toplist(), which restores the heap property with the same comparison function */
extern void qsort_toplist_r(toplist_t*list, int (*compare)(const void *, const void *));


/* access a certain element of the toplist (shouldn't be done by fiddling in toplist_t)
   returns a NULL pointer on error (i.e. index out of bounds) */
extern void* toplist_elem(toplist_t*list, size_t idx);


/* compare two toplists
   returns -1 if list1 is "smaller", 1 if list2 is "smaller", 0 if they are equal,
   2 if they are uncomparable (different data types or "smaller" functions */
extern int compare_toplists(toplist_t*list1, toplist_t*list2);


#endif /* HEAPTOPLIST_H - double inclusion protection */
