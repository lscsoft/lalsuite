/*
*  Copyright (C) 2007 Jolien Creighton, Peter Shawhan
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

#ifndef _STRINGINPUT_H
#define _STRINGINPUT_H

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif


/**
 * \addtogroup StringInput_h
 * \author Creighton, T. D.
 *
 * \brief Provides routines to parse \c CHARVectors into other LAL datatypes.
 *
 * \heading{Synopsis}
 * \code
 * #include "StringInput.h"
 * \endcode
 *
 * This header provides prototypes for routines that construct
 * LAL data structures using the data from a character string.  As in
 * standard C, a \e string is a block of non-null bytes of arbitrary
 * length, terminated by a null byte <tt>'\0'</tt>, and referred to by a
 * value of type <tt>CHAR *</tt> pointing to the first byte in the string.
 * It is not to be confused with a \c CHARVector, a LAL structure
 * referring to a block of data of a specified length, which may or may
 * not contain one or more instances of <tt>'\0'</tt>.
 *
 * In general, the routines under this header will have string inputs of
 * type <tt>const CHAR *</tt> (in order to allow, for instance, string
 * literals to be used as inputs), but will allocate \c CHARVector
 * structures to store string outputs.  Unless otherwise specified, these
 * outputs are guaranteed to contain at least one <tt>'\0'</tt> character,
 * so their \c data fields are valid strings.  It is the
 * responsibility of the calling routine to ensure that the string input
 * contains a terminating <tt>'\0'</tt> within the memory segment pointed
 * to by the <tt>CHAR *</tt> input, in order to avoid segmentation
 * violation.
 *
 * These routines are intended to work in conjunction with the functions
 * in <tt>StreamInput.h</tt> to add LAL robustness to otherwise ad-hoc data
 * input routines.  However, the functions in \ref StringInput.h are
 * fully LAL-compliant and use only LAL types, so they are included in
 * \c liblal proper.
 *
 * \heading{Constants}
 *
 * The following constants are format strings that can be used by the
 * various C <tt>scanf()</tt> or <tt>printf()</tt> functions to parse or
 * write sequences of characters corresponding to base LAL datatypes.
 * Since the C datatypes (\c short, \c int, \c long,
 * <tt>long long</tt>, \c float, \c double, etc.) do not have fixed
 * mappings to LAL base datatypes (\c INT2, \c INT4, \c INT8,
 * \c REAL4, \c REAL8, etc.), the appropriate format strings for
 * each LAL datatype must be determined at configuration time and set at
 * compile time.
 *
 * These format strings give only the conversion character preceded by
 * any length modifier according to the type (\c short, \c long,
 * etc.).  In particular they do \e not contain the initial
 * <tt>'%'</tt> character that initiates the conversion specification.
 * However, being <tt>\#define</tt>d string literals, they can be combined
 * with <tt>"%"</tt> string literals or more complicated format strings
 * through implicit concatenation.  Thus to scan \c string for a
 * \c UINT4 number \c n one would write:
 * \code
 * sscanf( string, "%" LAL_UINT4_FORMAT, &n );
 * \endcode
 * Similarly, to print a \c REAL8 number \c x with 12 digits
 * following the decimal place, one could use the following:
 * \code
 * printf( "%.12" LAL_REAL8_FORMAT, x );
 * \endcode
 * Of course, floating-point numbers are more commonly printed using the
 * <tt>"%e"</tt> conversion specifier, which does not generally require
 * type-dependent length modifiers.
 *
 * <table>
 * <tr><th>Name</th><th>Usual value</th></tr>
 * <tr><td>#LAL_INT2_FORMAT</td><td><tt>"hd"</tt></td></tr>
 * <tr><td>#LAL_INT4_FORMAT</td><td><tt>"d"</tt>  or <tt>"ld"</tt></td></tr>
 * <tr><td>#LAL_INT8_FORMAT</td><td><tt>"ld"</tt> or <tt>"lld"</tt></td></tr>
 * <tr><td>#LAL_UINT2_FORMAT</td><td><tt>"hu"</tt></td></tr>
 * <tr><td>#LAL_UINT4_FORMAT</td><td><tt>"u"</tt>  or <tt>"lu"</tt></td></tr>
 * <tr><td>#LAL_UINT8_FORMAT</td><td><tt>"lu"</tt> or <tt>"llu"</tt></td></tr>
 * <tr><td>#LAL_REAL4_FORMAT</td><td><tt>"f"</tt></td></tr>
 * <tr><td>#LAL_REAL8_FORMAT</td><td><tt>"lf"</tt></td></tr>
 * </table>
 *
 */
/*@{*/

/** \name Error Codes *//*@{ */
#define STRINGINPUTH_ENUL 1     /**< Unexpected null pointer in arguments */
#define STRINGINPUTH_EOUT 2     /**< Output handle points to a non-null pointer */
#define STRINGINPUTH_EMEM 3     /**< Memory allocation error */
/*@}*/
/** \cond DONT_DOXYGEN */
#define STRINGINPUTH_MSGENUL "Unexpected null pointer in arguments"
#define STRINGINPUTH_MSGEOUT "Output handle points to a non-null pointer"
#define STRINGINPUTH_MSGEMEM "Memory allocation error"
/** \endcond */

/**
 * This structure stores a number of null-terminated strings of arbitrary
 * length.  The entire list is stored flattened in a \c CHARVector,
 * and individual tokens are pointed to by a <tt>CHAR *[]</tt> handle.
 */
typedef struct tagTokenList {
    UINT4 nTokens;  /**< The number of tokens in the list */
    CHAR **tokens;  /**< A list of pointers to the individual tokens;
                     * the elements <tt>tokens[0..nTokens-1]</tt> point to tokens, and
                     * the element <tt>tokens[nTokens]</tt> is explicitly \c NULL (as is
                     the convention for an \c argv argument list */
    CHARVector *list;
                    /**< The flattened list of tokens, separated by (and terminated with) <tt>'\0'</tt> characters */
} TokenList;

/*@}*/

/* Function prototypes. */

void
LALCreateTokenList(LALStatus * status,
                   TokenList ** list,
                   const CHAR * string, const CHAR * delimiters);

void LALDestroyTokenList(LALStatus * status, TokenList ** list);


int
XLALCreateTokenList(TokenList ** list,
                    const CHAR * string, const CHAR * delimiters);

void XLALDestroyTokenList(TokenList * list);

void
LALStringToU2(LALStatus * status, UINT2 * value, const CHAR * string,
              CHAR ** endptr);

void
LALStringToU4(LALStatus * status, UINT4 * value, const CHAR * string,
              CHAR ** endptr);

void
LALStringToU8(LALStatus * status, UINT8 * value, const CHAR * string,
              CHAR ** endptr);

void
LALStringToI2(LALStatus * status, INT2 * value, const CHAR * string,
              CHAR ** endptr);

void
LALStringToI4(LALStatus * status, INT4 * value, const CHAR * string,
              CHAR ** endptr);

void
LALStringToI8(LALStatus * status, INT8 * value, const CHAR * string,
              CHAR ** endptr);

void
LALStringToS(LALStatus * status, REAL4 * value, const CHAR * string,
             CHAR ** endptr);

void
LALStringToD(LALStatus * status, REAL8 * value, const CHAR * string,
             CHAR ** endptr);

void
LALStringToC(LALStatus * status, COMPLEX8 * value, const CHAR * string,
             CHAR ** endptr);

void
LALStringToZ(LALStatus * status, COMPLEX16 * value, const CHAR * string,
             CHAR ** endptr);

void
LALStringToGPS(LALStatus * status, LIGOTimeGPS * value,
               const CHAR * string, CHAR ** endptr);

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _STRINGINPUT_H */
