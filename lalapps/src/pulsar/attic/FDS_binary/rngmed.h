/*
*  Copyright (C) 2007 Reinhard Prix
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

#ifndef RNGMED_H
#define RNGMED_H

int rngmed(const double *, unsigned int, unsigned int, double *);
int rngmed_sortindex(const void *, const void *); 

/*----------------------------
Structure for storing elements of
an array along with position
of the element in the array
---------------------------*/
    struct rngmed_val_index{
        double data;
        unsigned int index;
    };

#else
#endif

/**
 * \struct node
 * This structure is used to make a linked list. The list holds the samples
 * in one block in the running median algorithm.
 *
 * \param data Holds a single number.
 * \param next_sorted Points to the next node in the sorted list.
 * \param prev_sorted Points to the previous node in the sorted list.
 * \param next_sequence point to the next node in the sequential list.
 */

/**
 * \struct rngmed_val_index
 * A structure to store values and indices
 * of elements in an array
 *
 * \param data Stores a single number
 * \param index Stores the original position of the number
 *
 * This structure is used to track the indices
 * of elements after sorting by qsort. An array of
 * rngmed_val_index is passed to qsort which
 * rearranges the array according to the values in
 * the data member. The
 * indices of these elements in the original unsorted
 * array can then be read off from the index member.
 */
