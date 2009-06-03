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


/*
 * Prints Errors in the search to specified file
 */

#include "LALDoc.h"

int
LALDocErr(const char *msg, const char *inFileName, int laldocLine ,
          const char *laldocSourceFile, int fatal )
{

fprintf(ptrLALErrorFile,"\n     %s    \n \
   Input file: %s.                    \n \
   Input line Number: %i              \n \
   LALDoc Source-file line number: %i \n \
   LALDoc Source-file Name: %s        \n",
   msg, inFileName , lineNum, laldocLine , laldocSourceFile
);

/*
fprintf(stderr,"\n   %s            \n \
   Input file: %s.                   \n \
   Input line Number: %i             \n \
   LALDoc Source-file line number: %i\n \
   LALDoc Source-file Name: %s       \n",
   msg, inFileName , lineNum, laldocLine , laldocSourceFile
);
*/

if (fatal){
        CloseAFile( ptrLALErrorFile ,  1 );
        fprintf(ptrLALErrorFile,"Determined to be fatal. Exiting.");
        exit(BADEXIT);
}

return 0;
}

