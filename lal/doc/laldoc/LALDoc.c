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

/* static char *rcsid="$Id$"; */
/* Authors: Patrick Brady, Teviet Creighton and Alan Wiseman */
#include "LALDoc.h"

FILE *ptrLALErrorFile;   /* Pointer to file where errors are stored  */
int lineNum =0;          /* and current line number are gobal.       */

int main(int argc, char *argv[])
{
        char line[MAXSTR], sourceFileName[MAXSTR] ;
        FILE *ptrInputFile;
        LALEnvironment LaTeXEnv, VerbatimEnv , ErrTableEnv ;

        ParseCommandLine( argc , argv , sourceFileName , &ptrInputFile );

       /*
        *   Set the Enviroment flags:
       */
        SetLaTeXFlags   ( &LaTeXEnv    , sourceFileName , ptrInputFile );
        SetVerbatimFlags( &VerbatimEnv , sourceFileName , ptrInputFile );
        SetErrTableFlags( &ErrTableEnv , sourceFileName , ptrInputFile );



       /*
        *  Loop over the lines of code in the source file
       */
        while( LalDocGetLine(line,MAXSTR, ptrInputFile)>0  ){

               lineNum++;

               CheckEnvironment( line   ,  &LaTeXEnv     );
               CheckEnvironment( line   ,  &VerbatimEnv  );
               CheckEnvironment( line   ,  &ErrTableEnv  );
               CheckForEnvConflicts( &LaTeXEnv , &VerbatimEnv , &ErrTableEnv );

               if( LaTeXEnv.On    ) WriteLaTeX    ( line , &LaTeXEnv    ) ;
               if( VerbatimEnv.On ) WriteVerbatim ( line , &VerbatimEnv ) ;
               if( ErrTableEnv.On ) WriteErrTable ( line , &ErrTableEnv ) ;

        } /* end loop over lines */


        CleanExit( ptrInputFile , ptrLALErrorFile ,
                   &LaTeXEnv , &VerbatimEnv , &ErrTableEnv );

return 0;
}
