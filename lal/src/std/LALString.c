/*
*  Copyright (C) 2015, 2016 Karl Wette
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

#include <stdint.h>
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

/**
 * Append the formatted string 'fmt' to the string 's', which
 * is reallocated with XLALRealloc() to the required size.
 */
char *XLALStringAppendFmt(char *s, const char *fmt, ...)
{
  XLAL_CHECK_NULL(fmt != NULL, XLAL_EFAULT);
  const size_t n = (s == NULL) ? 0 : strlen(s);
  va_list ap;
  va_start(ap, fmt);
  char tmp[1];
  const int m = vsnprintf(tmp, sizeof(tmp), fmt, ap);
  va_end(ap);
  XLAL_CHECK_NULL(m >= 0, XLAL_ESYS, "vsnprintf('%s', ...) failed", fmt);
  const size_t l = (n + m + 1) * sizeof(*s);
  s = XLALRealloc(s, l);
  XLAL_CHECK_NULL(s != NULL, XLAL_ENOMEM, "XLALRealloc(n=%zu) failed", l);
  va_start(ap, fmt);
  XLAL_CHECK_NULL(vsnprintf(s + n, m + 1, fmt, ap) >= 0, XLAL_ESYS, "vsnprintf('%s', ...) failed", fmt);
  va_end(ap);
  return s;
}

/** Like strdup but uses LAL allocation routines (free with LALFree). */
char *XLALStringDuplicate(const char *s)
{
    char *dup;
    dup = XLALStringAppend(NULL, s);
    return dup;
}

/**
 * Copy sources string src to destination string dst.
 * Up to size - 1 characters are copied and destination string dst is
 * guaranteed to be NUL-terminated.
 * Return value is the length of source string src.  If this is greater than
 * or equal to the size of the destination string buffer, size, then truncation
 * has occurred. Should be nearly equivalent to strlcpy.
 */
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

/**
 * Concatenate sources string src to the end of destination string dst.
 * Characters are added to destination string dst until the source string src
 * is exhausted or the length of destination string dst is size - 1 characters.
 * Destination string dst is guaranteed to be NUL-terminated.
 * Return value is the initial length of destination string dst plus the
 * length of source string src.  If this is greater than
 * or equal to the size of the destination string buffer, size, then truncation
 * has occurred. Should be nearly equivalent to strlcat.
 */
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
 * Turn a string in-place into lowercase without using locale-dependent functions.
 */
int XLALStringToLowerCase(char * string)
{
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
 * Turn a string in-place into uppercase without using locale-dependent functions.
 */
int XLALStringToUpperCase(char * string)
{
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

/**
 * Compare two strings, ignoring case and without using locale-dependent functions.
 */
int XLALStringCaseCompare(const char *s1, const char *s2)
{
    size_t n = ( (s1 == NULL) ? 0 : strlen(s1) ) + ( (s2 == NULL) ? 0 : strlen(s2) );
    return XLALStringNCaseCompare(s1, s2, n);
}

/**
 * Compare the first N characters of two strings, ignoring case and without using locale-dependent functions.
 */
int XLALStringNCaseCompare(const char *s1, const char *s2, size_t n)
{

    /* ctype replacements w/o locale */
    const char upper_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const char lower_chars[] = "abcdefghijklmnopqrstuvwxyz";

    int c1 = 0, c2 = 0;

    /* based on implementation of strncmp() in glibc */
    while (n-- > 0) {

        c1 = (s1 == NULL) ? 0 : *s1++;
        c2 = (s2 == NULL) ? 0 : *s2++;

        /* convert c1 to lower case */
        if (c1) {
            char *p = strchr(upper_chars, c1);
            if (p) {
                int offset = p - upper_chars;
                c1 = lower_chars[offset];
            }
        }

        /* convert c2 to lower case */
        if (c2) {
            char *p = strchr(upper_chars, c2);
            if (p) {
                int offset = p - upper_chars;
                c2 = lower_chars[offset];
            }
        }

        /* compare characters */
        if (c1 == '\0' || c1 != c2) {
            return c1 - c2;
        }

    }

    return c1 - c2;

}

/**
 * Locates substring needle in string haystack, ignoring case and without
 * using locale-dependent functions.
 */
char * XLALStringCaseSubstring(const char *haystack, const char *needle)
{
    size_t haystack_length;
    size_t needle_length = strlen(needle);

    /* return haystack if needle is empty */
    if (needle_length == 0)
        return (char *)(intptr_t)(haystack);

    haystack_length = strlen(haystack);
    while (needle_length <= haystack_length) {
        if (XLALStringNCaseCompare(haystack, needle, needle_length) == 0)
            return (char *)(intptr_t)(haystack);
        --haystack_length;
        ++haystack;
    }

    /* needle not found in haystack */
    return NULL;
}

/**
 * Return the next token delimited by any character in 'delim' from the string 's', which is updated
 * to point just pass the returned token. If 'empty' is true, empty tokens are accepted.
 */
char *XLALStringToken(char **s, const char *delim, int empty)
{

    if (*s == NULL) {
        return NULL;
    }

    /* based on implementations of strtok_r() and strsep() in glibc */
    char *begin = *s;
    if (!empty) {
        begin += strspn(begin, delim);
        if (*begin == '\0') {
            *s = NULL;
            return NULL;
        }
    }
    char *end = strpbrk(begin, delim);
    if (end != NULL) {
        *end++ = '\0';
        *s = end;
    } else {
        *s = NULL;
    }

    return begin;

}

/**
 * Return the string 's' with all characters 'from' replaced with 'to'
 */
char *XLALStringReplaceChar(char *s, const int from, const int to)
{
    for (char *c = s; c != NULL && *c != '\0'; ++c) {
        if (*c == from) {
            *c = to;
        }
    }
    return s;
}
