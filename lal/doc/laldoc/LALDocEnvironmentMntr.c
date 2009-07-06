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

#include "LALDoc.h"

/*
 *
 * Check to see if enviroment switch is on this line,
 * and set the switch accordingly.
 *
 */
int CheckEnvironment(char *line, LALEnvironment *Env)
{
        int i;
        size_t szchr,szStr;
        szchr = sizeof(char);
        szStr = (size_t)(MAXSTR);
        /* If there is an ON FLAG on the line ...  */
        if(strstr(line, Env->OnFlag) ){
           /* we check to see it is already on .... */
           if(Env->On){
                 LALDocErr("Opening environment that is already open.",
                            Env->sourceFile , __LINE__ , __FILE__ , 1 );
           }

           else   /* if not on, then we turn it on and open an output file */
           {

                   Env->On = 1;

                   Env->fileName   = calloc( szStr , szchr );
                   Env->allCaps    = calloc( szStr , szchr );
                   Env->dfltFile   = calloc( szStr , szchr );
                   Env->errCodePrfx= calloc( szStr , szchr );
                   /* and initialize the strings */
                   for(i=0;i<MAXSTR;i++)  {
                           *(Env->fileName     +i ) = '\0';
                           *(Env->allCaps      +i ) = '\0';
                           *(Env->dfltFile     +i ) = '\0';
                           *(Env->errCodePrfx  +i ) = '\0';
                   }

                   FindNameOfOutPutFile( line , Env );
                   FindDefaultFileName(  Env );
                   /* if no file name specified, use the default */
                   if( *(Env->fileName) == '\0' )
                          for(i=0;i<MAXSTR;i++){
                                  *(Env->fileName +i) = *(Env->dfltFile +i);
                          }
                   Env->OutFilePtr = OpenAFile( Env->fileName,"a+",1 );
                   Env->Preamble ( Env );
           }

        }

        /* If there is an OFF FLAG on the line, turn it off the environment */
        if(strstr(line, Env->OffFlag) ){
           if(!(Env->On)){
                   LALDocErr("Closing environment that is already closed.",
                            Env->sourceFile , __LINE__ , __FILE__ , 1 );
           }
           else
           {
                   /*
                    * print any closing remarks to the file
                    * eg \end{verbatim}
                   */
                   Env->PostFix ( Env );
                   Env->On = 0; /* ie turn off the environment */

                   /* free the space and release the files */
                   CloseAFile(Env->OutFilePtr,1);
                   Env->lineCount = 0 ;

                   free( Env->fileName    );
                   free( Env->allCaps     );
                   free( Env->dfltFile    );
                   free( Env->errCodePrfx );


           }
        }
        return 0;
}

/*
 *       Check For Environment Conflicts
 *
 * Routine checks to see if there are any environment conflicts.
 * It is very forgiving.  It allows everything, but writing
 * different environments to the same file at the same time.
 *
*/
int
CheckForEnvConflicts( LALEnvironment *LaTeXEnv    ,
                      LALEnvironment *VerbatimEnv ,
                      LALEnvironment *ErrTableEnv )
{

    /* Note: The nested 'if' statements are necessary because strcmp()
     * will give a SEGV if it has not been intialized.
    */
    if( LaTeXEnv->On ){
            if( VerbatimEnv->On ){
                    if( (strcmp(LaTeXEnv->fileName,VerbatimEnv->fileName)!=0  )){
                            LALDocErr("LaTeX/Verbatim output conflict.",
                            LaTeXEnv->sourceFile , __LINE__ , __FILE__ , 1 );
                    }
            }
            if( ErrTableEnv->On ){
                    if( ( strcmp( LaTeXEnv->fileName, ErrTableEnv->fileName)!=0  )){
                            LALDocErr("LaTeX/ErrTable output conflict.",
                            LaTeXEnv->sourceFile , __LINE__ , __FILE__ , 1 );
                    }
            }
    }
    if( VerbatimEnv->On ){
            if( ErrTableEnv->On ){
                    if((strcmp( ErrTableEnv->fileName, VerbatimEnv->fileName)==0)){
                            LALDocErr("Verbatim/ErrTable output conflict.",
                            LaTeXEnv->sourceFile , __LINE__ , __FILE__ , 1 );
                    }
            }
    }

return 1;
}


int
CleanExit( FILE *ptrInputFile , FILE *ptrErrorFile,
           LALEnvironment *LaTeXEnv    ,
           LALEnvironment *VerbatimEnv ,
           LALEnvironment *ErrTableEnv )
{

    /* Check to see that all environments terminated */
    if( LaTeXEnv->On != 0) {
        LALDocErr("LaTeX environment never terminated.",
                   LaTeXEnv->sourceFile , __LINE__ , __FILE__ , 1 );
    }
    if( VerbatimEnv->On != 0) {
        LALDocErr("Verbatim environment never terminated.",
                   VerbatimEnv->sourceFile , __LINE__ , __FILE__ , 1 );
    }
    if( ErrTableEnv->On != 0) {
        LALDocErr("Error Table environment never terminated.",
                   ErrTableEnv->sourceFile , __LINE__ , __FILE__ , 1 );
    }


    /* close the input files */
    CloseAFile( ptrInputFile ,  0 );
    CloseAFile( ptrErrorFile ,  1 );
return 1;
}



/* K&R read a line,  return length */
int LalDocGetLine(char *line, int max, FILE *fpin)
{
   int i;
   for (i=0;i<MAXSTR;i++) { *(line+i) = '\0' ; }
   if (fgets(line, max, fpin) == NULL)
      return 0;
   else
      return strlen(line);
}

