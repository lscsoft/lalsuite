/*
 * Copyright (C) 2008 Reinhard Prix
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

#ifndef _STRINGVECTOR_H  /* Double-include protection. */
#define _STRINGVECTOR_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/*---------- INCLUDES ----------*/
#include <stdarg.h>

#include <lal/LALDatatypes.h>

/**
 * \author Reinhard Prix
 * \date 2008
 * \addtogroup StringVector_h
 *
 * \brief Creation/destruction/manipulation API for ::LALStringVector objects,
 *  which are just LAL-type vectors of CHAR * pointers.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/StringVector.h>
 * \endcode
 *
 */
/*@{*/

/*---------- DEFINES ----------*/

/*----- Error-codes -----*/
/*---------- exported types ----------*/
/*---------- Global variables ----------*/
/*---------- exported prototypes [API] ----------*/

#ifdef SWIG /* SWIG interface directives */
/* disable keywords arguments for XLALCreateStringVector */
%feature("kwargs", 0) XLALCreateStringVector;
/* ensure that SWIG generates "compact" default arguments,
   i.e. it only generates one wrapping function where all
   missing arguments are assigned NULL (it should do this
   anyway, since we're inside an extern "C" block and SWIG
   doesn't do any overloading for C-linkage functions) */
%feature("compactdefaultargs") XLALCreateStringVector;
/* add 20 optional arguments to XLALCreateStringVector,
   so that it can be called with up to 21 arguments */
%varargs(20, CHAR *arg = NULL) XLALCreateStringVector;
/* but since XLALCreateStringVector will crash if none
   of its arguments is NULL, add a contract to require
   that at least the final 21st argument is NULL */
%contract XLALCreateStringVector {
require:
  arg21 == NULL;
}
#endif /* SWIG */

LALStringVector *XLALCreateStringVector ( const CHAR *str1, ... );
LALStringVector *XLALAppendString2Vector (LALStringVector *vect, const CHAR *string );
void XLALDestroyStringVector ( LALStringVector *vect );

int XLALSortStringVector (LALStringVector *strings);
LALStringVector *XLALParseCSV2StringVector ( const CHAR *CSVlist );
INT4 XLALFindStringInVector ( const char *needle, const LALStringVector *haystack );

/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
