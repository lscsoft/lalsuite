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
