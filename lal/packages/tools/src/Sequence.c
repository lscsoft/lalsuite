/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#include <string.h>
#include <lal/Date.h>
#include <lal/LALComplex.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALError.h>
#include <lal/Sequence.h>
#include <lal/XLALError.h>

/*
 * Shift the bytes in the buffer buff, whose length is length, count bytes to
 * higher addresses.  If the magnitude of count is greater than or equal to
 * that of length, then nothing is done.
 */


static void memshift(void *buff, size_t length, int count)
{
	if(count >= 0) {
		if(length > (size_t) count)
			memmove((char *) buff + count, buff, length - count);
	} else {
		if(length > (size_t) -count)
			memmove(buff, (char *) buff - count, length + count);
	}
}


/*
 * Begin library functions...
 */

#define DATATYPE REAL4
#define SQUAREDATATYPE REAL4
#define ZERO 0
#define SUM sum += *data;
#define SUMSQ *data * *data
#include "Sequence_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ

#define DATATYPE REAL8
#define SQUAREDATATYPE REAL8
#define ZERO 0
#define SUM sum += *data;
#define SUMSQ *data * *data
#include "Sequence_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ

#define DATATYPE COMPLEX8
#define SQUAREDATATYPE REAL4
#define ZERO LAL_COMPLEX8_ZERO
#define SUM sum = XLALCOMPLEX8Add(sum, *data);
#define SUMSQ XLALCOMPLEX8Abs2(*data)
#include "Sequence_source.c"
#include "SequenceComplex_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ

#define DATATYPE COMPLEX16
#define SQUAREDATATYPE REAL8
#define ZERO LAL_COMPLEX16_ZERO
#define SUM sum = XLALCOMPLEX16Add(sum, *data);
#define SUMSQ XLALCOMPLEX16Abs2(*data)
#include "Sequence_source.c"
#include "SequenceComplex_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ

#define DATATYPE INT2
#define SQUAREDATATYPE UINT2
#define ZERO 0
#define SUM sum += *data;
#define SUMSQ *data * *data
#include "Sequence_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ


#define DATATYPE UINT2
#define SQUAREDATATYPE UINT2
#define ZERO 0
#define SUM sum += *data;
#define SUMSQ *data * *data
#include "Sequence_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ

#define DATATYPE INT4
#define SQUAREDATATYPE UINT4
#define ZERO 0
#define SUM sum += *data;
#define SUMSQ *data * *data
#include "Sequence_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ

#define DATATYPE UINT4
#define SQUAREDATATYPE UINT4
#define ZERO 0
#define SUM sum += *data;
#define SUMSQ *data * *data
#include "Sequence_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ

#define DATATYPE INT8
#define SQUAREDATATYPE UINT8
#define ZERO 0
#define SUM sum += *data;
#define SUMSQ *data * *data
#include "Sequence_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ

#define DATATYPE UINT8
#define SQUAREDATATYPE UINT8
#define ZERO 0
#define SUM sum += *data;
#define SUMSQ *data * *data
#include "Sequence_source.c"
#undef DATATYPE
#undef SQUAREDATATYPE
#undef ZERO
#undef SUM
#undef SUMSQ
