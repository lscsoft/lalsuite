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

#include <lal/LALStdlib.h>
#include <lal/LALDict.h>
#include <lal/LALDictSequence.h>

void XLALDestroyDictSequence(LALDictSequence *sequence)
{
    if (sequence) {
        if (sequence->data) {
            for (size_t i = 0; i < sequence->length; ++i)
                XLALDestroyDict(sequence->data[i]);
            LALFree(sequence->data);
        }
        LALFree(sequence);
    }
    return;
}

LALDictSequence * XLALCreateDictSequence(size_t length)
{
    LALDictSequence *sequence;
    sequence = LALMalloc(sizeof(*sequence));
    XLAL_CHECK_NULL(sequence, XLAL_ENOMEM);
    sequence->length = length;
    sequence->data = LALCalloc(length, sizeof(*sequence->data));
    if (sequence->data == NULL) {
        LALFree(sequence);
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    }
    for (size_t i = 0; i < length; ++i) {
        sequence->data[i] = XLALCreateDict();
        if (sequence->data[i] == NULL) {
            XLALDestroyDictSequence(sequence);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }
    }
    return sequence;
}

LALDictSequence * XLALCutDictSequence(const LALDictSequence *sequence, size_t first, size_t length)
{
    LALDictSequence *new;
    XLAL_CHECK_NULL(sequence, XLAL_EFAULT);
    XLAL_CHECK_NULL(first + length <= sequence->length, XLAL_EBADLEN);
    new = XLALCreateDictSequence(length);
    XLAL_CHECK_NULL(new, XLAL_EFUNC);
    for (size_t i = 0; i < length; ++i) {
        int retval = XLALDictUpdate(new->data[i], sequence->data[i + first]);
        if (retval < 0) {
            XLALDestroyDictSequence(new);
            XLAL_ERROR_NULL(XLAL_EFUNC);
        }
    }
    return new;
}

LALDictSequence * XLALCopyDictSequence(const LALDictSequence *sequence)
{
    XLAL_CHECK_NULL(sequence, XLAL_EFAULT);
    return XLALCutDictSequence(sequence, 0, sequence->length);
}

/* helper routine to reverse elements in sequence in range [start, end) */
static void reverse(LALDict **a, int start, int end)
{
    int left = start;
    int right = end - 1;
    while (left < right) {
        LALDict *tmp = a[left];
        a[left] = a[right];
        a[right] = tmp;
        ++left;
        --right;
    }
}

void XLALShiftDictSequence(LALDictSequence *sequence, int count)
{
    size_t abscount = count > 0 ? count : -count;

    XLAL_CHECK_VOID(sequence, XLAL_EFAULT);
    XLAL_CHECK_VOID(sequence->data || sequence->length == 0, XLAL_EINVAL);

    if (!count)
        return;

    /* shift memory if abs(count) < sequence->length */
    if (sequence->length > abscount) {
        int shift = count > 0 ? sequence->length - count : -count;
        reverse(sequence->data, 0, shift);
        reverse(sequence->data, shift, sequence->length);
        reverse(sequence->data, 0, sequence->length);
    }

    /* clear unshifted memory */
    if (abscount >= sequence->length) {
        for (size_t i = 0; i < sequence->length; ++i)
            XLALClearDict(sequence->data[i]);
    } else if (count > 0) {
        for (size_t i = 0; i < abscount; ++i)
            XLALClearDict(sequence->data[i]);
    } else {
        for (size_t i = sequence->length - abscount; i < sequence->length; ++i)
            XLALClearDict(sequence->data[i]);
    }

    return;
}

LALDictSequence * XLALResizeDictSequence(LALDictSequence *sequence, int first, size_t length)
{
    LALDict **data;
    XLAL_CHECK_NULL(sequence, XLAL_EFAULT);
    XLAL_CHECK_NULL(sequence->data || sequence->length == 0, XLAL_EINVAL);

    if (length > sequence->length) { /* need to increase memory */
        data = XLALRealloc(sequence->data, length * sizeof(*sequence->data));
        XLAL_CHECK_NULL(data, XLAL_EFUNC);
        for (size_t i = sequence->length; i < length; ++i)
            data[i] = XLALCreateDict();
        sequence->data = data;
        sequence->length = length;
        XLALShiftDictSequence(sequence, -first);
    } else if (length > 0) { /* need to decrease memory */
        XLALShiftDictSequence(sequence, -first);
        for (size_t i = length; i < sequence->length; ++i)
            XLALDestroyDict(sequence->data[i]);
        data = XLALRealloc(sequence->data, length * sizeof(*sequence->data));
        XLAL_CHECK_NULL(data, XLAL_EFUNC);
        sequence->data = data;
        sequence->length = length;
    } else { /* length == 0: need to release all memory */
        for (size_t i = 0; i < sequence->length; ++i)
            XLALDestroyDict(sequence->data[i]);
        XLALFree(sequence->data);
        sequence->data = NULL;
        sequence->length = 0;
    }

    return sequence;
}

size_t XLALDictSequenceLength(LALDictSequence *sequence)
{
    XLAL_CHECK(sequence, XLAL_EFAULT);
    return sequence->length;
}

LALDict * XLALDictSequenceGet(LALDictSequence *sequence, int pos)
{
    XLAL_CHECK_NULL(sequence, XLAL_EFAULT);
    XLAL_CHECK_NULL(sequence->data || sequence->length == 0, XLAL_EINVAL);
    if (pos < 0) /* negative positions are from end */
        pos += sequence->length;
    XLAL_CHECK_NULL(pos >= 0 && (unsigned)pos < sequence->length, XLAL_EBADLEN);
    return XLALDictDuplicate(sequence->data[pos]);
}

int XLALDictSequenceSet(LALDictSequence *sequence, LALDict *dict, int pos)
{
    XLAL_CHECK(sequence, XLAL_EFAULT);
    XLAL_CHECK(sequence->data || sequence->length == 0, XLAL_EINVAL);
    if (pos < 0) /* negative positions are from end */
        pos += sequence->length;
    XLAL_CHECK(pos >= 0 && (unsigned)pos < sequence->length, XLAL_EBADLEN);
    XLALDestroyDict(sequence->data[pos]);
    sequence->data[pos] = XLALDictDuplicate(dict);
    return 0;
}
