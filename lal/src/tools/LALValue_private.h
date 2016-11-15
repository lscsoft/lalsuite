/*
*  Copyright (C) 2016 Jolien Creighton
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

#ifndef _LAL_VALUE_PRIVATE_H
#define _LAL_VALUE_PRIVATE_H
#include <lal/LALDatatypes.h>
struct tagLALValue {
	LALTYPECODE type;
	size_t size;
	union { /* align for any type of data */
		CHAR i1;
		INT2 i2;
		INT4 i4;
		INT8 i8;
		UCHAR u1;
		UINT2 u2;
		UINT4 u4;
		UINT8 u8;
		REAL4 s;
		REAL8 d;
		COMPLEX8 c;
		COMPLEX16 z;
	} data[];
};
#endif /* _LAL_VALUE_PRIVATE_H */
