/*
*  Copyright (C) 2007 Jolien Creighton
*  Copyright (C) 2016 Kipp Cannon
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
 * File Name: InsertionSort.c
 *
 * Author: Cannon, Kipp
 *
 *
 *-----------------------------------------------------------------------*/

/* ---------- see Sort.h for doxygen documentation ---------- */

#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/Sort.h>


int XLALInsertionSort(void *base, size_t nobj, size_t size, void *params, int (*compar)(void *, const void *, const void *))
{
	char *i;
	char *end = (char *) base + nobj * size;
	void *temp;

	/* 0 or 1 objects are already sorted. */
	if(nobj < 2)
		return 0;

	temp = XLALMalloc(size);
	if(!temp)
		XLAL_ERROR(XLAL_ENOMEM);

	for(i = (char *) base + size; i < end; i += size) {
		char *j;
		size_t len;
		for(j = i; j > (char *) base && compar(params, j - size, i) > 0; j -= size);
		len = i - j;
		if(len) {
			memcpy(temp, i, size);
			memmove(j + size, j, len);
			memcpy(j, temp, size);
		}
	}

	XLALFree(temp);
	return 0;
}
