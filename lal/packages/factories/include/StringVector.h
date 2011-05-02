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

/**
 * \author Reinhard Prix
 * \date 2008
 * \file
 * \ingroup factories
 * \brief Creation/destruction/manipulation API for 'StringVector' type objects,
 *  which are just LAL 'vectors' of CHAR * pointers.
 *
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
#include <lal/LALRCSID.h>

/*---------- DEFINES ----------*/

/*----- Error-codes -----*/
/*---------- exported types ----------*/
/*---------- Global variables ----------*/
/*---------- exported prototypes [API] ----------*/

LALStringVector *XLALCreateStringVector ( const CHAR *str1, ... );
LALStringVector *XLALAppendString2Vector (LALStringVector *vect, const CHAR *string );
void XLALDestroyStringVector ( LALStringVector *vect );

int XLALSortStringVector (LALStringVector *strings);
LALStringVector *XLALParseCSV2StringVector ( const CHAR *CSVlist );
INT4 XLALFindStringInVector ( const char *needle, const LALStringVector *haystack );

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
