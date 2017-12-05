/*
 *  Copyright (C) 2017 Karl Wette
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

#ifndef _LALBITSET_H
#define _LALBITSET_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALBitset_h Header LALBitset.h
 * \ingroup lal_utilities
 * \author Karl Wette
 * \brief Implementation of an arbitrary-size bitset
 */
/*@{*/

/**
 * Arbitrary-size bitset
 */
typedef struct tagLALBitset LALBitset;

/**
 * Create a bitset
 */
LALBitset *XLALBitsetCreate(
  void
  );

/**
 * Destroy a bitset and its elements
 */
void XLALBitsetDestroy(
  LALBitset *bs                   /**< [in] Pointer to bitset */
  );

/**
 * Clear a bitset
 */
int XLALBitsetClear(
  LALBitset *bs                   /**< [in] Pointer to bitset */
  );

/**
 * Set/unset a bit in the bitset
 */
int XLALBitsetSet(
  LALBitset *bs,                  /**< [in] Pointer to bitset */
  const UINT8 idx,                /**< [in] Index of bit in set */
  const BOOLEAN is_set            /**< [in] Whether bit is set */
  );

/**
 * Get whether a bit in the bitset is set
 */
int XLALBitsetGet(
  const LALBitset *bs,            /**< [in] Pointer to bitset */
  const UINT8 idx,                /**< [in] Index of bit in set */
  BOOLEAN *is_set                 /**< [out] Whether bit is set */
  );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif // _LALBITSET_H
