/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 * Copyright (C) 2011 Adam Mercer, Leo Singer
 */

#ifndef _LOWLATENCYDATA_H
#define _LOWLATENCYDATA_H

#if defined(__cplusplus)
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

/**
   \defgroup LowLatencyData_h Header LowLatencyData.h
   \ingroup pkg_framedata
*/
/*@{*/

/** Structure representing reader's state and resources for monitoring /dev/shm. */
struct tagLowLatencyData;
typedef struct tagLowLatencyData LowLatencyData;

/**
 * Create a new instance of LowLatencyData.
 * It will be ready to read the next available frame file.
 */
LowLatencyData *XLALLowLatencyDataOpen(const char *data_path,   /* for example:  "/dev/shm/H1"    */
    const char *observatory,    /* for example:  "H"              */
    const char *frame_type      /* for example:  "H1_DMT_C00_L0"  */
    );

/**
 * Destroy an instance of LowLatencyData.
 * Release resources associated with monitoring /dev/shm.
 */
void XLALLowLatencyDataClose(LowLatencyData *);

/**
 * Retrieve the next available frame file.
 * On success, a pointer to a new in-memory frame file is returned.  If size is
 * a non-NULL pointer, then the size of the file in bytes is written to (*size).
 * Any prevoiusly returned pointer from this instance of LowLatencyData is
 * invalidated.
 *
 * On failure, NULL is returned and (*size) is not modified.  It is unspecified
 * whether any previously returned pointer is invalidated on failure.
 */
void *XLALLowLatencyDataNextBuffer(LowLatencyData *, size_t * size);

/*@}*/

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}       /* extern "C" */
#endif

#endif /* _LOWLATENCYDATA_h */
