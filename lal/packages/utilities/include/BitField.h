/*
 *  Copyright (C) 2007, 2008 Karl Wette
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

#ifndef _BITFIELD_H
#define _BITFIELD_H

/**
 * \addtogroup BitField_h
 * \author Karl Wette
 * \brief Macros for manipulating integers as bit fields
 */
/*@{*/

/**
 * Return a mask where the (zero-based) ith bit is set
 */
#define ONE_BIT(T, i) (((T)1) << (i))

/**
 * Return a mask where all bits from 0 to n-1 are set
 */
#define ALL_BITS(T, n) ((((T)1) << (n)) - ((T)1))

/**
 * Get the value of the (zero-based) ith bit of x
 */
#define GET_BIT(T, x, i) (((x) & ONE_BIT(T, i)) == ONE_BIT(T, i) ? 1 : 0)

/**
 * Sets the (zero-based) ith bit of x to the truth of v
 */
#define SET_BIT(T, x, i, v) ((v) ? ((x) |= ONE_BIT(T, i)) : ((x) &= ~ONE_BIT(T, i)))

/**
 * Get if all bits from 0 to n-1 of x are set
 */
#define GET_ALL(T, x, n) (((x) & ALL_BITS(T, n)) == ALL_BITS(T, n) ? 1 : 0)

/**
 * Sets all bits from 0 to n of x to the truth of v
 */
#define SET_ALL(T, x, n, v) ((v) ? ((x) |= ALL_BITS(T, n)) : ((x) &= ~ALL_BITS(T, n)))

/*@}*/

#endif
