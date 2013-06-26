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

/** \file
 * \addtogroup LALString_h
 * \author Creighton, J. D. E.
 * \brief XLAL string manipulation routines.
 *
 */

#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>

/** Like strcat but dynamically reallocates string with LALRealloc. */
char *XLALStringAppend(char *s, const char *append)
{
    size_t curlen;
    size_t newlen;
    if (!append)
        return s;
    curlen = s ? strlen(s) : 0;
    newlen = curlen + strlen(append);
    s = LALRealloc(s, newlen + 1);
    if (!s)
        XLAL_ERROR_NULL(XLAL_ENOMEM);
    strcpy(s + curlen, append);
    return s;
}

/** Like strdup but uses LAL allocation routines (free with LALFree). */
char *XLALStringDuplicate(const char *s)
{
    char *dup;
    dup = XLALStringAppend(NULL, s);
    return dup;
}

/** Copy sources string src to destination string dst.
 * Up to size - 1 characters are copied and destination string dst is
 * guaranteed to be NUL-terminated.
 * Return value is the length of source string src.  If this is greater than
 * or equal to the size of the destination string buffer, size, then truncation
 * has occurred. Should be nearly equivalent to strlcpy. */
size_t XLALStringCopy(char *dst, const char *src, size_t size)
{
    size_t srclen;
    if (!src)
        src = "";
    srclen = strlen(src);
    if (!dst || size < 1)       /* no copy */
        return srclen;
    if (size == 1) {    /* NUL terminate and exit */
        dst[0] = 0;
        return srclen;
    }
    strncpy(dst, src, size - 1);
    dst[size - 1] = 0;
    return srclen;
}

/** Concatenate sources string src to the end of destination string dst.
 * Characters are added to destination string dst until the source string src
 * is exhausted or the length of destination string dst is size - 1 characters.
 * Destination string dst is guaranteed to be NUL-terminated.
 * Return value is the initial length of destination string dst plus the
 * length of source string src.  If this is greater than
 * or equal to the size of the destination string buffer, size, then truncation
 * has occurred. Should be nearly equivalent to strlcat. */
size_t XLALStringConcatenate(char *dst, const char *src, size_t size)
{
    size_t srclen;
    size_t dstlen;
    if (!src)
        src = "";
    srclen = strlen(src);
    if (!dst || size < 1)       /* no copy */
        return srclen;
    if (size == 1) {    /* NUL terminate and exit */
        dst[0] = 0;
        return srclen;
    }
    dstlen = strlen(dst);
    strncat(dst, src, size - dstlen - 1);
    return srclen + dstlen;
}


/**
 * Helper function:  turn a string in-place into lowercase without
 * using locale-dependent functions.
 */
int XLALStringToLowerCase(CHAR * string)
{
/**< [in/out] string to convert */
    XLAL_CHECK(string != NULL, XLAL_EINVAL);

    /* ctype replacements w/o locale */
    const char upper_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const char lower_chars[] = "abcdefghijklmnopqrstuvwxyz";

    for (UINT4 i = 0; i < strlen(string); i++) {
        int c = string[i];
        if (c) {
            char *p = strchr(upper_chars, c);
            if (p) {
                int offset = p - upper_chars;
                c = lower_chars[offset];
            }
        }
        string[i] = c;

    }   // for i < len(string)

    return XLAL_SUCCESS;

}       /* XLALStringToLowerCase() */


/**
 * Helper function:  turn a string in-place into uppercase without
 * using locale-dependent functions.
 */
int XLALStringToUpperCase(CHAR * string)
{
/**< [in/out] string to convert */
    XLAL_CHECK(string != NULL, XLAL_EINVAL);

    /* ctype replacements w/o locale */
    const char upper_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const char lower_chars[] = "abcdefghijklmnopqrstuvwxyz";

    for (UINT4 i = 0; i < strlen(string); i++) {
        int c = string[i];
        if (c) {
            char *p = strchr(lower_chars, c);
            if (p) {
                int offset = p - lower_chars;
                c = upper_chars[offset];
            }
        }
        string[i] = c;

    }   // for i < len(string)

    return XLAL_SUCCESS;

}       /* XLALStringToUpperCase() */
