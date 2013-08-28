/*
*  Copyright (C) 2007 Jolien Creighton
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

/*-----------------------------------------------------------------------
 *
 * File Name: Sort.h
 *
 * Author: Creighton, T. D.
 *
 *
 *-----------------------------------------------------------------------*/

#ifndef _SORT_H
#define _SORT_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif


/**
 * \addtogroup Sort_h
 *
 * \brief Provides routines for sorting, indexing, and ranking real vector elements.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/Sort.h>
 * \endcode
 *
 * \heading{Description}
 *
 * These routines sort a vector <tt>*data</tt> (of type \c REAL4Vector
 * or \c REAL8Vector) into ascending order using the in-place
 * heapsort algorithm, or construct an index vector <tt>*index</tt> that
 * indexes <tt>*data</tt> in increasing order (leaving <tt>*data</tt>
 * unchanged), or construct a rank vector <tt>*rank</tt> that gives the
 * rank order of the corresponding <tt>*data</tt> element.
 *
 * The relationship between sorting, indexing, and ranking can be a bit
 * confusing.  One way of looking at it is that the original array is
 * ordered by index, while the sorted array is ordered by rank.  The
 * index array gives the index as a function of rank; i.e.\ if you're
 * looking for a given rank (say the 0th, or smallest element), the index
 * array tells you where to look it up in the unsorted array:
 * \code
 * unsorted_array[index[i]] = sorted_array[i]
 * \endcode
 * The rank array gives the rank as a function of index; i.e.\ it tells
 * you where a given element in the unsorted array will appear in the
 * sorted array:
 * \code
 * unsorted_array[j] = sorted_array[rank[j]]
 * \endcode
 * Clearly these imply the following relationships, which can be used to
 * construct the index array from the rank array or vice-versa:
 * \code
 * index[rank[j]] = j
 * rank[index[i]] = i
 * \endcode
 *
 * The XLAL versions of these routines, \c XLALHeapSort(), \c XLALHeapIndex(),
 * and \c XLALHeapRank(), perform the same operations but on arrays of
 * \c nobj generic objects of size \c size pointed to by \c base and
 * using the comparison function \c compar.  The function \c compar() has
 * the prototype
 * \code
 * int compar( void *p, const void *x, const void *y )
 * \endcode
 * and returns \f$-1\f$ if \f${\mathtt{x}}<{\mathtt{y}}\f$,
 * \f$0\f$ if \f${\mathtt{x}}={\mathtt{y}}\f$,
 * and \f$+1\f$ if \f${\mathtt{x}}>{\mathtt{y}}\f$.  Here \c p (which may be NULL)
 * is a pointer to additional data that may be used in the comparison function.
 * This pointer is passed to the comparison function unchanged from the argument
 * \c params of \c XLALHeapSort(), \c XLALHeapIndex(), and
 * \c XLALHeapRank().
 *
 * \heading{Algorithm}
 *
 * These routines use the standard heap sort algorithm described in
 * Sec. 8.3 of Ref. [\ref ptvf1992].
 *
 * The <tt>LALSHeapSort()</tt> and <tt>LALDHeapSort()</tt> routines are entirely
 * in-place, with no auxiliary storage vector.  The <tt>LALSHeapIndex()</tt>
 * and <tt>LALDHeapIndex()</tt> routines are also technically in-place, but
 * they require two input vectors (the data vector and the index vector),
 * and leave the data vector unchanged.  The <tt>LALSHeapRank()</tt> and
 * <tt>LALDHeapRank()</tt> routines require two input vectors (the data and
 * rank vectors), and also allocate a temporary index vector internally;
 * these routines are therefore the most memory-intensive.  All of these
 * algorithms are \f$N\log_2(N)\f$ algorithms, regardless of the ordering of
 * the initial dataset.
 *
 * Note: if you can use \c qsort(), you should.
 *
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define SORTH_ENUL 1		/**< Null pointer */
#define SORTH_ELEN 2		/**< Length mismatch */
#define SORTH_EMEM 3		/**< Memory allocation error */
/*@}*/

/** \cond DONT_DOXYGEN */
#define SORTH_MSGENUL "Null pointer"
#define SORTH_MSGELEN "Length mismatch"
#define SORTH_MSGEMEM "Memory allocation error"
/** \endcond */

/* Function prototypes. */

/* ----- HeapSort.c ----- */

/** \see See \ref Sort_h for documentation */
int XLALHeapSort( void *base, UINT4 nobj, UINT4 size, void *params,
    int (*compar)(void *, const void *, const void *) );

/** \see See \ref Sort_h for documentation */
int XLALHeapIndex( INT4 *indx, void *base, UINT4 nobj, UINT4 size, void *params,
    int (*compar)(void *, const void *, const void *) );

/** \see See \ref Sort_h for documentation */
int XLALHeapRank( INT4 *rank, void *base, UINT4 nobj, UINT4 size, void *params,
    int (*compar)(void *, const void *, const void *) );

/** \see See \ref Sort_h for documentation */
void LALSHeapSort(LALStatus      *status,
	       REAL4Vector *vector);

/** \see See \ref Sort_h for documentation */
void LALSHeapIndex(LALStatus      *status,
		INT4Vector  *indx,
		REAL4Vector *vector);

/** \see See \ref Sort_h for documentation */
void LALSHeapRank(LALStatus      *status,
	       INT4Vector  *rank,
	       REAL4Vector *vector);

/** \see See \ref Sort_h for documentation */
void LALDHeapSort(LALStatus      *status,
	       REAL8Vector *vector);

/** \see See \ref Sort_h for documentation */
void LALDHeapIndex(LALStatus      *status,
		INT4Vector  *indx,
		REAL8Vector *vector);

/** \see See \ref Sort_h for documentation */
void LALDHeapRank(LALStatus      *status,
	       INT4Vector  *rank,
	       REAL8Vector *vector);


/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _SORT_H */
