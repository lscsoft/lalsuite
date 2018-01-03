/*
 * Copyright (C) 2011 Reinhard Prix
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

/*---------- INCLUDES ----------*/
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <lal/LALStdlib.h>
#include <lal/XLALError.h>
#include <lal/StringVector.h>
#include <lal/LALMalloc.h>

/*---------- DEFINES ----------*/
/*---------- internal types ----------*/
/*---------- empty initializers ---------- */
/*---------- Global variables ----------*/
/*---------- internal prototypes ----------*/
int XLALStringVector_TEST ( void );

/**
 * \author Reinhard Prix
 * \file
 * \ingroup StringVector_h
 */

/*==================== FUNCTION DEFINITIONS ====================*/

int main( int argc, char *argv[] )
{

  /* sanity checks */
  if ( argc != 1 )
    XLAL_ERROR ( XLAL_EINVAL, "The executable '%s' doesn't support any input arguments right now.\n", argv[0] );

  if ( XLALStringVector_TEST() != XLAL_SUCCESS )
    XLAL_ERROR ( XLAL_EFUNC, "StringVector TEST failed.\n" );

  return XLAL_SUCCESS;

} /* main() */

/**
 * Test various StringVector functions
 */
int
XLALStringVector_TEST ( void )
{
  return XLAL_SUCCESS;
}
