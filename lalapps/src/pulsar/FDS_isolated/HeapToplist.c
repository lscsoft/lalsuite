/*
*  Copyright (C) 2007 Bernd Machenschalk
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

   A Heap is a partially sorted structure that is typically represented as a (binary) tree.
   A tree is a heap if for all nodes the value of the node is smaller than that of all its
   successors according to a given comparison function.
   Footnote: in most other implementation of a heap the order is different from here, i.e.
   replace "smaller" with "greater" in the above.

   The nice thing about such a heap for our application (and others) is that this heap property
   doesn't imply any special relation between two nodes in different branches of the tree, so no
   effort is necessary to keep such a relation. This allows to perform all operations on the heap
   (including removal and insertion of elements) with O(log n), i.e. at most the depth of the tree.

   Here the tree is a binary tree and stored in an array where the successors succ(n) of a node
   with index n have indices 2*n+1 and 2*n+2:

            0
       1          2
     3   4     5     6
    7 8 9 10 11 12 13 14

   There's nothing special about that - look for "Heapsort" in the WWW or an algorithms book.

   Bernd Machenschalk
*/


#include <stdlib.h>
#include <string.h>
#include "HeapToplist.h"

/* this function gets a "partial heap", i.e. a heap where only the top
   element (potentially) violates the heap property. It "bubbles
   down" this element so that the heap property is restored */
static void down_heap(toplist_t*list) {
  size_t node = 0;
  size_t succ;
  char *exch;
  while ((succ = node+node+1) < list->elems) {
    if (succ+1 < list->elems)
      if ((list->smaller)((list->heap)[succ+1], (list->heap)[succ]) > 0)
	succ++;
    if ((list->smaller)((list->heap)[succ], (list->heap)[node]) > 0) {
      exch = (list->heap)[node];
      (list->heap)[node] = (list->heap)[succ];
      (list->heap)[succ] = exch;
      node = succ;
    } else
      break;
  }
}


/* this function gets a "partial heap", i.e. a heap where only an element on
   the lowest level (potentially) violates the heap property. "node" is the 
   index of this element. The function "bubbles up" this element so that the
   heap property is restored */
static void up_heap(toplist_t*list, size_t node) {
  size_t pred;
  char *exch;
  while (node > 0) {
    pred = (node-1)/2;
    if ((list->smaller)((list->heap)[node], (list->heap)[pred]) > 0) {
      exch = (list->heap)[node];
      (list->heap)[node] = (list->heap)[pred];
      (list->heap)[pred] = exch;
      node = pred;
    } else
      break;
  }
}


/* creates a toplist with length elements,
   returns -1 on error (out of memory), else 0 */
int create_toplist(toplist_t**list,
		   size_t length,
		   size_t size,
		   int (*smaller)(const void *, const void *)) {
  toplist_t *listp;

  if (!(listp = malloc(sizeof(toplist_t))))
    return(-1);
  if (!(listp->data = malloc(size * length))) {
    free(listp);
    return(-1);
  }
  if (!(listp->heap = malloc(sizeof(char*) * length))) {
    free(listp->data);
    free(listp);
    return(-1);
  }
  
  listp->length  = length;
  listp->elems   = 0;
  listp->size    = size;
  listp->smaller = smaller;

  *list = listp;
  return(0);
}


/* clears an existing toplist of all elements inserted so far */
void clear_toplist(toplist_t*list) {
  list->elems = 0;
}


/* frees the space occupied by the toplist */
void free_toplist(toplist_t**list) {
    free((*list)->heap);
    free((*list)->data);
    free(*list);
}


/* Inserts an element in to the toplist either if there is space left
   or the element is larger than the smallest element in the toplist.
   In the latter case, remove the smallest element from the toplist.
   Returns 1 if the element was actually inserted, 0 if not. */
int insert_into_toplist(toplist_t*list, void *element) {

  /* if there is room left, add it at the end (and update the heap) */
  if (list->elems < list->length) {
    list->heap[list->elems] = list->data + list->elems * list->size;
    memcpy(list->heap[list->elems], element, list->size);
    list->elems++;
    up_heap(list,list->elems-1);
    return(1);

  /* if it is smaller than the smallest element, simply drop it.
     if it is bigger, replace the smallest element (root of the heap)
     and update the heap */
  } else if ((list->smaller)(element, (list->heap)[0]) < 0) {
    memcpy(list->heap[0], element, list->size);
    down_heap(list);
    return(1);

  } else
    return(0);
}


/* apply the function "handle" to all elements of the list in the current order */
void go_through_toplist(toplist_t*list, void (*handle)(void *)) {
  size_t i;
  for(i=0;i<list->elems;i++)
    handle(list->heap[i]);
}

void* toplist_elem(toplist_t*list, size_t ind) {
  if (list == NULL)
    return(NULL);
  if (ind >= list->elems)
    return(NULL);
  return(list->heap[ind]);
}

int compare_toplists(toplist_t*list1, toplist_t*list2) {
  size_t i=0;
  int res=0;
  if ((list1 == NULL) || (list2 == NULL))
    return(2);
  if ((list1->smaller != list2->smaller) ||
      (list1->size    != list2->size))
    return(2);
  while((i < list1->elems) &&
	(i < list2->elems) &&
	(res == 0)) {
    res = list1->smaller(toplist_elem(list1,i),toplist_elem(list2,i));
    i++;
  }
  if (res == 0) {
    if (list1->elems < list2->elems)
      return(1);
    else if (list1->elems > list2->elems)
      return(-1);
  }
  return(res);
}


/* using qsort requires some "global" help */

/* global function pointer for qsort */
static int (*_qsort_compare1)(const void*, const void*);
/* wrapper function for qsort */
static int _qsort_compare2(const void*a, const void*b){
  return (_qsort_compare1(*(void*const*)a,*(void*const*)b));
}
/* inverse wrapper */
static int _qsort_compare3(const void*b, const void*a){
  return (_qsort_compare1(*(void*const*)a,*(void*const*)b));
}

/* sorts the toplist with an arbitrary sorting function
   (potentially) destroying the heap property */
void qsort_toplist(toplist_t*list, int (*compare)(const void*, const void*)) {
  /* point the global function pointer to compare, then call qsort with the wrapper */
  _qsort_compare1 = compare;
  qsort(list->heap,list->elems,sizeof(char*),_qsort_compare2);
}

/* qsort function that gives the reverse ordering of the previous */
void qsort_toplist_r(toplist_t*list, int (*compare)(const void*, const void*)) {
  /* point the global function pointer to compare, then call qsort with the wrapper */
  _qsort_compare1 = compare;
  qsort(list->heap,list->elems,sizeof(char*),_qsort_compare3);
}
