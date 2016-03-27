/*
 *  Copyright (C) 2016 John Veitch
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

#ifndef _LALHASHFUNC_H
#define _LALHASHFUNC_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALHashFunc_h Header LALHashFunc.h
 * \ingroup lal_utilities
 * \authors John Veitch, Karl Wette
 * \brief Implementations of various hash functions
 */
/*@{*/

/**
 * \name Pearson hash function
 */
/*@{*/

/**
 * Compute a arbitrary-sized Pearson hash value for the given arbitrary data
 */
int XLALPearsonHash(
  void *hval,			/**< [in/out] Hash value; should be either zero, or the result
                                   of a previous XLALPearsonHash() call when hashing multiple data */
  const size_t hval_len,	/**< [in] Length of hash value */
  const void *data,		/**< [in] Arbitrary data to hash */
  const size_t data_len		/**< [in] Length of arbitrary data */
);

/*@}*/

/**
 * \name CityHash hash functions
 */
/*@{*/

/** Hash function for a byte array. */
UINT8 XLALCityHash64(const char *buf, size_t len);

/** Hash function for a byte array.  For convenience, a 64-bit seed is also
 * hashed into the result.
 */
UINT8 XLALCityHash64WithSeed(const char *buf, size_t len, UINT8 seed);

/** Hash function for a byte array.  For convenience, two seeds are also
 * hashed into the result.
 */
UINT8 XLALCityHash64WithSeeds(const char *buf, size_t len,
                           UINT8 seed0, UINT8 seed1);

/** Hash function for a byte array.  Most useful in 32-bit binaries. */
UINT4 XLALCityHash32(const char *buf, size_t len);

/*@}*/

/*@}*/

#ifdef __cplusplus
}
#endif

#endif // _LALHASHFUNC_H
