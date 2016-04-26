/*
 *  Copyright (C) 2016 Karl Wette
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

#include <lal/LALHeap.h>

#define LEFT(i)     (2*(i) + 1)         /* Left child of binary heap element 'i' */
#define RIGHT(i)    (2*(i) + 2)         /* Right child of binary heap element 'i' */
#define PARENT(i)   (((i) - 1)/2)       /* Parent of binary heap element 'i' */

/* Swap elements x and y */
#define SWAP(x, y)   do { void *z = (x); (x) = (y); (y) = z; } while (0)

/* Evaluates true if the pair (x, y) is NOT ordered, according to the heap compare function */
#define UNORDERED(h, x, y)   ((h)->cmp((h)->cmp_param, (x), (y)) * (h)->min_or_max_heap > 0)

struct tagLALHeap {
  void **data;                  /* Binary heap data */
  int data_len;                 /* Size of the memory block 'data', in number of elements */
  int n;                        /* Number of valid elements in the heap */
  LALHeapDtorFcn dtor;          /* Function to free memory of elements of heap, if required */
  int max_size;                 /* Maximum size of the heap; if zero, heap has unlimited size */
  int min_or_max_heap;          /* -1|+1 if root of heap is minimum|maximum element */
  LALHeapCmpParamFcn cmp;       /* Parameterised heap element comparison function */
  void *cmp_param;              /* Parameter to pass to comparison function */
};

/* Call a non-parameterised compare function, which is passed in 'param' */
static int heap_no_param_cmp( void *param, const void *x, const void *y )
{
  LALHeapCmpFcn cmp = ( LALHeapCmpFcn ) param;
  return cmp( x, y );
}

/* Resize the binary heap to twice the current number of heap elements */
static int heap_resize( LALHeap *h )
{
  int new_len = ( h->n > 0 ) ? 2*h->n : 1;
  h->data = XLALRealloc( h->data, new_len * sizeof( h->data[0] ) );
  XLAL_CHECK( h->data != NULL, XLAL_ENOMEM );
  h->data_len = new_len;
  return XLAL_SUCCESS;
}

/* Swap element 'i' with its parent until the heap property is satisfied */
static void heap_bubble_up( LALHeap *h, int i )
{
  int p = PARENT( i );
  while ( i > 0 && UNORDERED( h, h->data[i], h->data[p] ) ) {
    SWAP( h->data[i], h->data[p] );
    i = p;
    p = PARENT( i );
  }
}

/* Swap element 'i' with either of its children until the heap property is satisfied */
static void heap_trickle_down( LALHeap *h, int i )
{
  do {
    int j = -1, r = RIGHT( i );
    if ( r < h->n && UNORDERED( h, h->data[r], h->data[i] ) ) {
      int l = LEFT( i );
      if ( UNORDERED( h, h->data[l], h->data[r] ) ) {
        j = l;
      } else {
        j = r;
      }
    } else {
      int l = LEFT( i );
      if ( l < h->n && UNORDERED( h, h->data[l], h->data[i] ) ) {
        j = l;
      }
    }
    if ( j >= 0 ) {
      SWAP( h->data[i], h->data[j] );
    }
    i = j;
  } while ( i >= 0 );
}

LALHeap *XLALHeapCreate(
  LALHeapDtorFcn dtor,
  int max_size,
  int min_or_max_heap,
  LALHeapCmpFcn cmp
  )
{

  /* Create a heap using heap_no_param_cmp as the comparison function */
  LALHeap *h = XLALHeapCreate2( dtor, max_size, min_or_max_heap, heap_no_param_cmp, cmp );
  XLAL_CHECK_NULL( h != NULL, XLAL_EFUNC );

  return h;

}

LALHeap *XLALHeapCreate2(
  LALHeapDtorFcn dtor,
  int max_size,
  int min_or_max_heap,
  LALHeapCmpParamFcn cmp,
  void *cmp_param
  )
{

  /* Check input */
  XLAL_CHECK_NULL( max_size >= 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( abs( min_or_max_heap ) == 1, XLAL_EINVAL );
  XLAL_CHECK_NULL( cmp != NULL, XLAL_EFAULT );

  /* Allocate memory for heap struct */
  LALHeap *h = XLALCalloc( 1, sizeof( *h ) );
  XLAL_CHECK_NULL( h != NULL, XLAL_ENOMEM );

  /* Set heap struct parameters */
  h->dtor = dtor;
  h->max_size = max_size;
  h->min_or_max_heap = min_or_max_heap;
  h->cmp = cmp;
  h->cmp_param = cmp_param;

  return h;

}

void XLALHeapDestroy(
  LALHeap *h
  )
{
  if ( h != NULL ) {
    if ( h->data != NULL ) {
      if ( h->dtor != NULL ) {
        for ( int i = 0; i < h->n; ++i ) {
          h->dtor( h->data[i] );
        }
      }
      XLALFree( h->data );
    }
    XLALFree( h );
  }
}

int XLALHeapSize(
  const LALHeap *h
  )
{
  XLAL_CHECK( h != NULL, XLAL_EFAULT );
  return h->n;
}

const void *XLALHeapRoot(
  const LALHeap *h
  )
{
  XLAL_CHECK_NULL( h != NULL, XLAL_EFAULT );
  return ( h->n > 0 ) ? h->data[0] : NULL;
}

int XLALHeapResize(
  LALHeap *h,
  int max_size
  )
{

  /* Check input */
  XLAL_CHECK( h != NULL, XLAL_EFAULT );
  XLAL_CHECK( max_size >= 0, XLAL_EINVAL );

  /* If too many heap elements for new (fixed) maximum size, remove some */
  while ( max_size > 0 && h->n > max_size ) {
    XLAL_CHECK( XLALHeapRemoveRoot( h ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  /* Set new maximum size */
  h->max_size = max_size;

  return XLAL_SUCCESS;

}

int XLALHeapAdd(
  LALHeap *h,
  void **x
  )
{

  /* Check input */
  XLAL_CHECK( h != NULL, XLAL_EFAULT );
  XLAL_CHECK( x != NULL, XLAL_EFAULT );
  XLAL_CHECK( *x != NULL, XLAL_EINVAL );

  if ( h->max_size == 0 || h->n < h->max_size ) { /* If heap is unlimited, or not full yet */

    /* Resize binary heap; designed so that resizing costs amortized constant time */
    if ( h->n + 1 > h->data_len ) {
      XLAL_CHECK( heap_resize( h ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

    /* Add new element to end of binary heap, and bubble up to restore heap property */
    h->data[h->n++] = *x;
    *x = NULL;
    heap_bubble_up( h, h->n - 1 );

  } else if ( UNORDERED( h, h->data[0], *x ) ) { /* If new element should replace root */
    XLAL_CHECK( XLALHeapExchangeRoot( h, x ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}

void *XLALHeapExtractRoot(
  LALHeap *h
  )
{

  /* Check input */
  XLAL_CHECK_NULL( h != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( h->n > 0, XLAL_ESIZE );

  /* Save root element */
  void *x = h->data[0];

  /* Replace root with last element in binary heap, and trickle down to restore heap property */
  h->data[0] = h->data[--h->n];
  heap_trickle_down( h, 0 );

  /* Resize binary heap; designed so that resizing costs amortized constant time */
  if ( 3*h->n < h->data_len ) {
    XLAL_CHECK_NULL( heap_resize( h ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return x;

}

int XLALHeapRemoveRoot(
  LALHeap *h
  )
{

  /* Check input */
  XLAL_CHECK( h != NULL, XLAL_EFAULT );
  XLAL_CHECK( h->n > 0, XLAL_ESIZE );

  /* Extract root */
  void *x = XLALHeapExtractRoot( h );
  XLAL_CHECK( x != NULL, XLAL_EFUNC );

  /* Free memory associated with root element, if required */
  if ( h->dtor != NULL ) {
    h->dtor( x );
  }

  return XLAL_SUCCESS;

}

int XLALHeapExchangeRoot(
  LALHeap *h,
  void **x
  )
{

  /* Check input */
  XLAL_CHECK( h != NULL, XLAL_EFAULT );
  XLAL_CHECK( h->n > 0, XLAL_ESIZE );
  XLAL_CHECK( x != NULL, XLAL_EFAULT );

  /* Swap root with *x, and trickle down new root to restore heap property */
  SWAP( h->data[0], *x );
  heap_trickle_down( h, 0 );

  return XLAL_SUCCESS;

}

int XLALHeapVisit(
  LALHeap *h,
  LALHeapVisitFcn visit,
  void *visit_param
  )
{

  /* Check input */
  XLAL_CHECK( h != NULL, XLAL_EFAULT );
  XLAL_CHECK( visit != NULL, XLAL_EFAULT );

  /* Create internal min-heap (without destructor) to get elements in order */
  LALHeap *h2 = XLALHeapCreate2( NULL, 0, -1, h->cmp, h->cmp_param );
  XLAL_CHECK( h2 != NULL, XLAL_EFUNC );

  /* Add all elements in heap to internal min-heap */
  for ( int i = 0; i < h->n; ++i ) {
    void *x = h->data[i];
    XLAL_CHECK( XLALHeapAdd( h2, &x ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  /* Remove all elements from original heap */
  h->n = 0;

  /* Extract roots element of internal min-heap, until empty */
  while ( h2->n > 0 ) {
    void *x = XLALHeapExtractRoot( h2 );
    XLAL_CHECK( x != NULL, XLAL_EFUNC );

    /* Visit element, possibly modifying it */
    XLAL_CHECK( visit( visit_param, x ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* Add element back to original heap */
    XLAL_CHECK( XLALHeapAdd( h, &x ) == XLAL_SUCCESS, XLAL_EFUNC );

  }

  /* Cleanup */
  XLALHeapDestroy( h2 );

  return XLAL_SUCCESS;

}
