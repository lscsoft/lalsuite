#include <stdlib.h>
#include <string.h>
#include "HeapToplist.h"

/* successors of node n: 2*n+1 and 2*n+2
          0
    1          2
 3     4    5     6
7 8  9 10 11 12 13 14
*/

/* this function gets a "partial heap", i.e. a heap where
   only the top element (potentially) violates the heap property. It "bubbles
   down" this element so that the heap property is restored */
static void down_heap(toplist_t*list) {
  UINT8 node = 0;
  UINT8 succ;
  char *exch;
  while ((succ = node+node+1) < list->elems) {
    if (succ <= list->elems)
      if ((list->smaller)((*(list->heap)[succ]), (*(list->heap)[succ+1])) > 0)
	succ++;
    if ((list->smaller)((*(list->heap)[node]), (*(list->heap)[succ])) > 0) {
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
  while ((pred = (node-1)/2) > 0) {
    if ((list->smaller)((*(list->heap)[node]), (*(list->heap)[pred])) > 0) {
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
  toplist_p *listp;

  if (!(listp = malloc(sizeof(toplist_p))))
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
  listp->lemes   = 0;
  listp->isheap  = !0;
  listp->smaller = smaller;
  *list = listp;
  return(0);
}

/* frees the space occupied by the toplist */
void free_toplist(toplist_t**list) {
    free((*list)->heap);
    free((*list)->data);
    free(list);
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
  } else if ((list->smaller)(element, (*(list->heap)[0])) < 0) {
    memcpy(list->heap[0], element, list->size);
    down_heap();
    return(1);
  } else
    return(0);
}

/* sorts the toplist with a sorting function different from "smaller",
   (potentially) destroying the heap property */
void sort_toplist_f(toplist_t*list, int (*compare)(const void *, const void *)) {
}

/* apply the function "handle" to all elements of the list in the current sorting order */
void go_through_toplist(toplist_t*list, int (*handle)(const void *)) {
  UINT8 i;
  for(i=0;i<list->elems;i++)
    handle(list->heap[i]);
}

