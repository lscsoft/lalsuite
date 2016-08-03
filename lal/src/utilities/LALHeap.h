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

#ifndef _LALHEAP_H
#define _LALHEAP_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALHeap_h Header LALHeap.h
 * \ingroup lal_utilities
 * \author Karl Wette
 * \brief Implementation of a generic heap, following Chapter 10.1 of \cite open-data-structs .
 */
/*@{*/

/**
 * Generic heap with elements of type <tt>void *</tt>
 */
typedef struct tagLALHeap LALHeap;

/**
 * Function which free memory associated with heap element <tt>x</tt>
 */
typedef void ( *LALHeapDtorFcn )( void *x );

/**
 * Function which compares heap elements <tt>x</tt> and <tt>y</tt>
 */
typedef int ( *LALHeapCmpFcn )( const void *x, const void *y );

/**
 * Function which compares heap elements <tt>x</tt> and <tt>y</tt>, with a parameter \c param
 */
typedef int ( *LALHeapCmpParamFcn )( void *param, const void *x, const void *y );

/**
 * Function to call when visiting heap element <tt>x</tt>, with a parameter \c param.
 * Return XLAL_SUCCESS if successful, or XLAL_FAILURE otherwise.
 */
typedef int ( *LALHeapVisitFcn )( void *param, const void *x );

/**
 * Function to call when visiting (and possibly modify) heap element <tt>x</tt>, with a parameter \c param.
 * Return XLAL_SUCCESS if successful, or XLAL_FAILURE otherwise.
 */
typedef int ( *LALHeapModifyFcn )( void *param, void *x );

/**
 * Create a heap
 */
LALHeap *XLALHeapCreate(
  LALHeapDtorFcn dtor,          /**< [in] Heap element destructor function, if required */
  int max_size,                 /**< [in] Maximum size of the heap; if zero, heap has unlimited size */
  int min_or_max_heap,          /**< [in] -1|+1 if root of heap is minimum|maximum element */
  LALHeapCmpFcn cmp             /**< [in[ Heap element comparison function */
  );

/**
 * Create a heap with a parameterised comparison function
 */
LALHeap *XLALHeapCreate2(
  LALHeapDtorFcn dtor,          /**< [in] Heap element destructor function, if required */
  int max_size,                 /**< [in] Maximum size of the heap; if zero, heap has unlimited size */
  int min_or_max_heap,          /**< [in] -1|+1 if root of heap is minimum|maximum element */
  LALHeapCmpParamFcn cmp,       /**< [in] Parameterised heap element comparison function */
  void *cmp_param               /**< [in] Parameter to pass to comparison function */
  );

/**
 * Destroy a heap and its elements
 */
void XLALHeapDestroy(
  LALHeap *h                    /**< [in] Pointer to heap */
  );

/**
 * Return the size of a heap
 */
int XLALHeapSize(
  const LALHeap *h              /**< [in] Pointer to heap */
  );

/**
 * Return the maximum size of a heap
 */
int XLALHeapMaxSize(
  const LALHeap *h              /**< [in] Pointer to heap */
  );

/**
 * Return the root element of a heap
 */
const void *XLALHeapRoot(
  const LALHeap *h              /**< [in] Pointer to heap */
  );

/**
 * Change the maximum size of a heap; if the heap is contracted, excess elements are destroyed
 */
int XLALHeapResize(
  LALHeap *h,                   /**< [in] Pointer to heap */
  int max_size                  /**< [in] New maximum size of the heap; if zero, heap has unlimited size */
  );

/**
 * Add a new element to a heap; if the heap is of fixed size and full, the root element is removed
 */
int XLALHeapAdd(
  LALHeap *h,                   /**< [in] Pointer to heap */
  void **x                      /**< [in/out] Pointer to new element. If the root element is removed, it
                                   is returned in <tt>*x</tt>; otherwise <tt>*x</tt> is set to \c NULL */
  );

/**
 * Remove the root element of a heap
 */
void *XLALHeapExtractRoot(
  LALHeap *h                    /**< [in] Pointer to heap */
  );

/**
 * Remove and destroy the root element of a heap
 */
int XLALHeapRemoveRoot(
  LALHeap *h                    /**< [in] Pointer to heap */
  );

/**
 * Exchange the root element of a non-empty heap with the new element in <tt>*x</tt>
 */
int XLALHeapExchangeRoot(
  LALHeap *h,                   /**< [in] Pointer to heap */
  void **x                      /**< [in/out] Pointer to new element. On return, contains root element */
  );

/**
 * Visit each element in the heap in the order given by the comparison function
 */
int XLALHeapVisit(
  const LALHeap *h,             /**< [in] Pointer to heap */
  LALHeapVisitFcn visit,        /**< [in] Visitor function to call for each heap element */
  void *visit_param             /**< [in] Parameter to pass to visitor function */
  );

/**
 * Visit (and possibly modify) each element in the heap in the order given by the comparison function
 */
int XLALHeapModify(
  LALHeap *h,                   /**< [in] Pointer to heap */
  LALHeapModifyFcn modify,      /**< [in] Modifier function to call for each heap element */
  void *modify_param            /**< [in] Parameter to pass to modifier function */
  );

/**
 * Allocate and return an array containing all elements in the heap
 */
const void **XLALHeapElements(
  const LALHeap *h              /**< [in] Pointer to heap */
  );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif // _LALHEAP_H
