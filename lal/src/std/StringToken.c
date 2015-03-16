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
 * \file
 * \ingroup StringInput_h
 * \author Creighton, T. D.
 *
 * \brief Converts a string into a series of tokens, for use by other routines.
 *
 * ### Description ###
 *
 * The routine <tt>XLALCreateTokenList()</tt> parses <tt>*string</tt> as a
 * sequence of tokens (substrings of non-null characters that do not
 * appear in \c delimiters), separated by delimiters (substrings
 * consisting only of characters that appear in \c delimiters), and
 * terminated by the null character <tt>'\0'</tt>.  The structure
 * <tt>**list</tt> is created, storing the sequence of tokens as a list
 * null-terminated character strings.
 *
 * The output \c list should be a non-\c NULL handle that points
 * to the value \c NULL (i.e.\ \c list\f$\neq\f$\c NULL but
 * <tt>*list</tt>=\c NULL).  Even if no tokens were found, <tt>*list</tt>
 * will be created, but will have <tt>(*list)->nTokens</tt>=0,
 * <tt>(*list)->tokens[0]</tt>=\c NULL, and
 * <tt>(*list)->list</tt>=\c NULL.  Note that this is \e not an
 * error, so the calling routine need not guarantee in advance that
 * \c string contain any non-delimiter characters.
 *
 * The routine <tt>XLALDestroyTokenList()</tt> destroys a list of tokens as
 * created by <tt>XLALCreateTokenList()</tt>, setting <tt>*list</tt> to
 * \c NULL.
 *
 * ### Algorithm ###
 *
 * The <tt>XLALCreateTokenList()</tt> function is not particularly
 * memory-efficient, requiring internal storage up to twice the length of
 * <tt>*string</tt>.  It first creates a working copy of
 * <tt>string->data</tt>, and replaces all occurences of characters
 * appearing in <tt>*delimiters</tt> with <tt>'\0'</tt>, while at the same
 * time keeping track of the number and total length of all tokens.  It
 * then allocates a contiguous block of memory to store all the tokens
 * (separated by and terminated with single <tt>'\0'</tt> characters), and
 * a set of <tt>CHAR *</tt> pointers to point to the individual tokens in
 * this block.  Then the routine proceeds through the working copy one
 * last time, copying tokens into the token list and setting the token
 * pointers accordingly, before destroying the working copy.
 *
 */

#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/StringInput.h>

/**
 * \deprecated Use XLALCreateTokenList() instead
 */
void
LALCreateTokenList(LALStatus * stat,
                   TokenList ** list,
                   const CHAR * string, const CHAR * delimiters)
{
    BOOLEAN delimiter = 1;      /* whether current character is a delimiter */
    UINT4 i = 0, j = 0; /* indecies */
    UINT4 nTokens = 0;  /* number of tokens */
    UINT4 sLength;      /* length of string */
    UINT4 tLength = 0;  /* length of token list */
    CHAR *copy; /* working copy of token list */

    INITSTATUS(stat);
    ATTATCHSTATUSPTR(stat);

    /* Check for valid input arguments. */
    ASSERT(list, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(string, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(delimiters, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(!*list, stat, STRINGINPUTH_EOUT, STRINGINPUTH_MSGEOUT);

    /* Create working copy of token list. */
    sLength = strlen(string) + 1;
    if (!(copy = (CHAR *) LALMalloc(sLength * sizeof(CHAR)))) {
        ABORT(stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    }
    for (i = 0; i < sLength; i++) {
        CHAR c = string[i];
        if (strchr(delimiters, c)) {
            copy[i] = '\0';
            delimiter = 1;
        } else {
            copy[i] = c;
            tLength++;
            if (delimiter) {
                delimiter = 0;
                nTokens++;
            }
        }
    }

    /* Create the token list. */
    if (!(*list = (TokenList *) LALMalloc(sizeof(TokenList)))) {
        LALFree(copy);
        ABORT(stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    }
    if (!((*list)->tokens =
          (CHAR **) LALMalloc((nTokens + 1) * sizeof(CHAR *)))) {
        LALFree(*list);
        *list = NULL;
        LALFree(copy);
        ABORT(stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    }
    (*list)->nTokens = nTokens;
    (*list)->list = NULL;


    /* If tokens were found, copy them over and set up pointers. */
    if (nTokens) {
        CHAR *listData; /* pointer to token list data */
        LALCHARCreateVector(stat->statusPtr, &((*list)->list),
                            nTokens + tLength);
        BEGINFAIL(stat) {
            LALFree((*list)->tokens);
            LALFree(*list);
            *list = NULL;
            LALFree(copy);
        }
        ENDFAIL(stat);
        listData = (*list)->list->data;
        i = 0;
        while (i < sLength) {
            if (copy[i]) {
                tLength = strlen(copy + i) + 1;
                memcpy(listData, copy + i, tLength * sizeof(CHAR));
                (*list)->tokens[j++] = listData;
                i += tLength;
                listData += tLength;
            } else
                i++;
        }
    }
    (*list)->tokens[j] = NULL;

    /* Clean up and exit. */
    LALFree(copy);
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}

/** Split given input string into a list of 'tokens' separated by any
 * of the characters given in 'delimiters'
 */
int XLALCreateTokenList(TokenList ** list,         //!< [out] list of tokens
                        const CHAR * string,       //!< [in] string to split into tokens
                        const CHAR * delimiters    //!< [in] set of token-delimiter characters
                        )
{
    XLAL_CHECK((list != NULL) && ((*list) == NULL), XLAL_EINVAL);
    XLAL_CHECK(string != NULL, XLAL_EINVAL);
    XLAL_CHECK(delimiters != NULL, XLAL_EINVAL);

    // prepare output TokenList structure
    TokenList *ret;
    XLAL_CHECK((ret = XLALCalloc(1, sizeof(*ret))) != NULL, XLAL_ENOMEM);

    size_t stringLen = strlen(string);
    if ((ret->list = XLALCreateCHARVector(stringLen + 1)) == NULL) {
        XLALFree(ret);
        XLAL_ERROR(XLAL_ENOMEM);
    }
    strcpy(ret->list->data, string);

    // initialize pointers to walk along local copy of input string
    char *ptr = ret->list->data;
    const char *endPtr = ptr + stringLen;

    UINT4 nTokens = 0;
    UINT4 nTokensAlloc = 0;

    while ((ptr != NULL) && (ptr < endPtr)) {
        // skip and nuke delimiter
        size_t skip = strspn(ptr, delimiters);
        memset(ptr, 0, skip);   // fill with '0'
        ptr += skip;

        if (ptr >= endPtr) {
            break;
        }
        // 'ptr' points at next token
        nTokens++;

        // allocate next batch of token-pointers if required
        if (nTokens > nTokensAlloc) {
            nTokensAlloc = 2 * nTokens; // proceed by doubling current space
            if ((ret->tokens = XLALRealloc(ret->tokens, nTokensAlloc * sizeof(char *))) == NULL) {
                XLALDestroyTokenList(ret);
                XLAL_ERROR(XLAL_ENOMEM);
            }
        }       // if nTokens > nTokensAlloc

        // enter new token-pointer into list
        ret->tokens[nTokens - 1] = ptr;

        // advance to next delimiter
        ptr = strpbrk(ptr, delimiters);

    }   // while ptr < endPtr

    // reduce tokens-array to actual size
    if ((ret->tokens = XLALRealloc(ret->tokens, nTokens * sizeof(char *))) == NULL) {
        XLALDestroyTokenList(ret);
        XLAL_ERROR(XLAL_ENOMEM);
    }
    ret->nTokens = nTokens;

    // return result
    (*list) = ret;

    return XLAL_SUCCESS;

}       // XLALCreateTokenList()

/**
 * \deprecated Use XLALDestroyTokenList() instead
 */
void LALDestroyTokenList(LALStatus * stat, TokenList ** list)
{
    INITSTATUS(stat);
    ATTATCHSTATUSPTR(stat);

    /* Check for valid input arguments. */
    ASSERT(list, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);
    ASSERT(*list, stat, STRINGINPUTH_ENUL, STRINGINPUTH_MSGENUL);

    /* Free everything and exit. */
    if ((*list)->list) {
        TRY(LALCHARDestroyVector(stat->statusPtr, &((*list)->list)), stat);
    }
    LALFree((*list)->tokens);
    LALFree(*list);
    *list = NULL;
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}

/** See StringToken.c for documentation */
void XLALDestroyTokenList(TokenList * list)
{
    /* Free everything and exit. */
    if (list) {
        if (list->list)
            XLALDestroyCHARVector(list->list);
        XLALFree(list->tokens);
        XLALFree(list);
    }
}
