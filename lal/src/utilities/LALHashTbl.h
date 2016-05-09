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

#ifndef _LALHASHTBL_H
#define _LALHASHTBL_H

#include <lal/LALStdlib.h>
#include <lal/LALHashFunc.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALHashTbl_h Header LALHashTbl.h
 * \ingroup lal_utilities
 * \author Karl Wette
 * \brief Implementation of a generic hash table, following Chapter 5.2 of \cite open-data-structs .
 */
/*@{*/

/**
 * Generic hash table with elements of type <tt>void *</tt>
 */
typedef struct tagLALHashTbl LALHashTbl;

/**
 * Function which free memory associated with hash table element <tt>x</tt>
 */
typedef void ( *LALHashTblDtorFcn )( void *x );

/**
 * Hash function for the hash table element <tt>x</tt>
 */
typedef UINT8( *LALHashTblHashFcn )( const void *x );

/**
 * Hash function for the hash table element <tt>x</tt>, with a parameter \c param
 */
typedef UINT8( *LALHashTblHashParamFcn )( void *param, const void *x );

/**
 * Function which compares hash table elements <tt>x</tt> and <tt>y</tt>
 */
typedef int ( *LALHashTblCmpFcn )( const void *x, const void *y );

/**
 * Function which compares hash table elements <tt>x</tt> and <tt>y</tt>, with a parameter \c param
 */
typedef int ( *LALHashTblCmpParamFcn )( void *param, const void *x, const void *y );

/**
 * Create a hash table
 */
LALHashTbl *XLALHashTblCreate(
  LALHashTblDtorFcn dtor,       /**< [in] Function to free memory of elements of hash, if required */
  LALHashTblHashFcn hash,       /**< [in] Hash function for hash table elements */
  LALHashTblCmpFcn cmp          /**< [in] Hash table element comparison function */
  );

/**
 * Create a hash table with parameterised hash and comparison functions
 */
LALHashTbl *XLALHashTblCreate2(
  LALHashTblDtorFcn dtor,       /**< [in] Function to free memory of elements of hash, if required */
  LALHashTblHashParamFcn hash,  /**< [in] Parameterised hash function for hash table elements */
  void *hash_param,             /**< [in] Parameter to pass to hash function */
  LALHashTblCmpParamFcn cmp,    /**< [in] Parameterised hash table element comparison function */
  void *cmp_param               /**< [in] Parameter to pass to comparison function */
  );

/**
 * Destroy a hash table
 */
void XLALHashTblDestroy(
  LALHashTbl *ht
  );

/**
 * Return the size of a hash table
 */
int XLALHashTblSize(
  const LALHashTbl *ht          /**< [in] Pointer to hash table */
  );

/**
 * Find the element matching <tt>x</tt> in a hash table; if found, return in <tt>*y</tt>
 */
int XLALHashTblFind(
  const LALHashTbl *ht,         /**< [in] Pointer to hash table */
  const void *x,                /**< [in] Hash element to match */
  const void **y                /**< [out] Pointer to matched hash element, or NULL if not found */
  );

/**
 * Add an element to a hash table
 */
int XLALHashTblAdd(
  LALHashTbl *ht,               /**< [in] Pointer to hash table */
  void *x                       /**< [in] Hash element to add */
  );

/**
 * Find the element matching <tt>x</tt> in a hash table; if found, remove it and return in <tt>*y</tt>
 */
int XLALHashTblExtract(
  LALHashTbl *ht,               /**< [in] Pointer to hash table */
  const void *x,                /**< [in] Hash element to match */
  void **y                      /**< [out] Pointer to matched hash element, which has been removed from
                                   the hash table, or NULL if not found */
  );

/**
 * Find the element matching <tt>x</tt> in a hash table; if found, remove and destroy it
 */
int XLALHashTblRemove(
  LALHashTbl *ht,               /**< [in] Pointer to hash table */
  const void *x                 /**< [in] Hash element to match */
  );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif // _LALHASHTBL_H
