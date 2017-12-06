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

/**
 * \addtogroup LALStdio_h
 *
 * \brief Provides LAL functions similar to the non-file functions in <tt><stdio.h></tt>.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LALStdio.h>
 * #include <lal/FileIO.h>
 * \endcode
 *
 */

#ifndef _LALSTDIO_H
#define _LALSTDIO_H

#include <stdio.h>
#include <stdarg.h>
#include <inttypes.h>
#include <lal/LALConfig.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

#if LAL_BOINC_ENABLED
    extern FILE *boinc_fopen(const char *path, const char *mode);
#define LALFopen boinc_fopen
#else
#define LALFopen fopen
#endif
#define LALFclose fclose


#define LAL_INT2_PRId PRId16
#define LAL_INT2_PRIi PRIi16
#define LAL_INT2_PRIo PRIo16
#define LAL_INT2_PRIu PRIu16
#define LAL_INT2_PRIx PRIx16
#define LAL_INT2_PRIX PRIX16

#define LAL_INT4_PRId PRId32
#define LAL_INT4_PRIi PRIi32
#define LAL_INT4_PRIo PRIo32
#define LAL_INT4_PRIu PRIu32
#define LAL_INT4_PRIx PRIx32
#define LAL_INT4_PRIX PRIX32

#define LAL_INT8_PRId PRId64
#define LAL_INT8_PRIi PRIi64
#define LAL_INT8_PRIo PRIo64
#define LAL_INT8_PRIu PRIu64
#define LAL_INT8_PRIx PRIx64
#define LAL_INT8_PRIX PRIX64

#define LAL_INT2_SCNd SCNd16
#define LAL_INT2_SCNi SCNi16
#define LAL_INT2_SCNo SCNo16
#define LAL_INT2_SCNu SCNu16
#define LAL_INT2_SCNx SCNx16

#define LAL_INT4_SCNd SCNd32
#define LAL_INT4_SCNi SCNi32
#define LAL_INT4_SCNo SCNo32
#define LAL_INT4_SCNu SCNu32
#define LAL_INT4_SCNx SCNx32

#define LAL_INT8_SCNd SCNd64
#define LAL_INT8_SCNi SCNi64
#define LAL_INT8_SCNo SCNo64
#define LAL_INT8_SCNu SCNu64
#define LAL_INT8_SCNx SCNx64

/* convenient versions of above that can be used in
 * either scanf or printf (decimal integers only) */
#define LAL_INT2_FORMAT  LAL_INT2_SCNd
#define LAL_INT4_FORMAT  LAL_INT4_SCNd
#define LAL_INT8_FORMAT  LAL_INT8_SCNd
#define LAL_UINT2_FORMAT LAL_INT2_SCNu
#define LAL_UINT4_FORMAT LAL_INT4_SCNu
#define LAL_UINT8_FORMAT LAL_INT8_SCNu
#define LAL_REAL4_FORMAT "g"
#define LAL_REAL8_FORMAT "lg"

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LALSTDIO_H */
