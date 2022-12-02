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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef _LALSTRING_H
#define _LALSTRING_H

#include <stddef.h>
#include <lal/LALAtomicDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

/**
 * \defgroup LALString_h Header LALString.h
 * \ingroup lal_std
 * \author Creighton, J. D. E.
 * \brief XLAL string manipulation routines.
 *//** @{ */

int XLALStringPrint(char *s, size_t n, const char *fmt, ...) _LAL_GCC_PRINTF_FORMAT_(3,4);
char *XLALStringAppend(char *s, const char *append);
char *XLALStringAppendFmt(char *s, const char *fmt, ...) _LAL_GCC_PRINTF_FORMAT_(2,3);
char *XLALStringDuplicate(const char *s);
size_t XLALStringCopy(char *dst, const char *src, size_t size);
size_t XLALStringConcatenate(char *dst, const char *src, size_t size);
int XLALStringToLowerCase(char * string);
int XLALStringToUpperCase(char * string);
int XLALStringCaseCompare(const char *s1, const char *s2);
int XLALStringNCaseCompare(const char *s1, const char *s2, size_t n);
char *XLALStringCaseSubstring(const char *haystack, const char *needle);
char *XLALStringToken(char **s, const char *delim, int empty);
char *XLALStringTranslate(char *s, int (*f)(int, void*), void *param);
char *XLALStringStripChars(char *s, int (*f)(int));
char *XLALStringKeepChars(char *s, int (*f)(int));
char *XLALStringReplaceChar(char *s, const int from, const int to);

/** @} */

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LALSTRING_H */
