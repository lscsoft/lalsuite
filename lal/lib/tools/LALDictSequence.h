/*
*  Copyright (C) 2023 Jolien Creighton
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef _LAL_DICTSEQUENCE_H
#define _LAL_DICTSEQUENCE_H

#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDict.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

struct tagLALDictSequence {
    size_t length;
    LALDict **data;
};
typedef struct tagLALDictSequence LALDictSequence;

void XLALDestroyDictSequence(LALDictSequence *sequence);
LALDictSequence * XLALCreateDictSequence(size_t length);
LALDictSequence * XLALCutDictSequence(const LALDictSequence *sequence, size_t first, size_t length);
LALDictSequence * XLALCopyDictSequence(const LALDictSequence *sequence);
void XLALShiftDictSequence(LALDictSequence *sequence, int count);
LALDictSequence * XLALResizeDictSequence(LALDictSequence *sequence, int first, size_t length);

size_t XLALDictSequenceLength(LALDictSequence *sequence);
LALDict * XLALDictSequenceGet(LALDictSequence *sequence, int pos);
int XLALDictSequenceSet(LALDictSequence *sequence, LALDict *inparmas, int pos);

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LAL_DICTSEQUENCE_H */
