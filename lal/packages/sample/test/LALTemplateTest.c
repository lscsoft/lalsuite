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

/****************************** <lalVerbatim file="LALTemplateTestCV">
Author: Creighton, T. D.
$Id$
******************************* </lalVerbatim> */

/* Usage and description */

#include <lal/LALTemplate.h>  /* Include any required headers */

/* Define RCS ID string */
NRCSID(LALTEMPLATETESTC,"$Id$");

/* Define local constants and macros */

/* Declare and set global lalDebugLevel */
int lalDebugLevel = 0;

/* Declare local (static) functions (definitions can go here or at the
   end of the file) */

/* Define main function */
/* int main( int argc, char **argv ) */
int main( void )
{
  /* Variable declarations */
  static LALStatus stat;

  /* Code */
  LALTemplate( &stat );
  return 0;
}
