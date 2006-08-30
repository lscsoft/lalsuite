#include <stdlib.h>
#include <string.h>
#include "HeapToplist.h"

static volatile const char *HEAPTOPLISTCID  = "$Id$";

/* successors of node n: 2*n+1 and 2*n+2
        0
   1          2
 3   4     5     6
7 8 9 10 11 12 13 14
*/


/* this function gets a "partial heap", i.e. a heap where only the top
   element (potentially) violates the heap property. It "bubbles
   down" this element so that the heap property is restored */
static void down_heap(toplist_t*list) {
  UINT8 node = 0;
  UINT8 succ;
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
   the lowest level (potentially) violates the heap property. It "bubbles up"
   this element so that the heap property is restored */
static void up_heap(toplist_t*list, UINT8 node) {
  UINT8 pred;
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
		   UINT8 length,
		   size_t size,
		   int (*smaller)(const void *, const void *)) {
  toplist_t *listp;

  if (!(listp = malloc(sizeof(toplist_t))))
    return(-1);
  if (!(listp->data = malloc(size * length))) {
    free(listp);
    return(-1);
  }
  if (!(listp->heap = malloc(sizeof(void*) * length))) {
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
void go_through_toplist(toplist_t*list, void (*handle)(const void *)) {
  UINT8 i;
  for(i=0;i<list->elems;i++)
    handle(list->heap[i]);
}

void* toplist_elem(toplist_t*list, UINT8 index) {
  if (list == NULL)
    return(NULL);
  if (index >= list->elems)
    return(NULL);
  return(list->heap[index]);
}

int compare_toplists(toplist_t*list1, toplist_t*list2) {
  UINT8 i=0, res=0;
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
  return (_qsort_compare1(*(void**)a,*(void**)b));
}
/* inverse wrapper */
static int _qsort_compare3(const void*b, const void*a){
  return (_qsort_compare1(*(void**)a,*(void**)b));
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
