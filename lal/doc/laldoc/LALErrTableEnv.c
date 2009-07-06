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
 * Function that sets flag for Error Table environment.
 *
 */
int SetErrTableFlags( LALEnvironment *ErrTableEnv, const char *sourcefile , FILE *inSrc )
{
        ErrTableEnv->OnFlag       = "<lalErrTable"   ;
        ErrTableEnv->OffFlag      = "</lalErrTable>" ;
        ErrTableEnv->closer       = ">"              ;
        ErrTableEnv->On           = 0                ;
        ErrTableEnv->InFilePtr    = inSrc            ;
        ErrTableEnv->OutFilePtr   = NULL             ;
        ErrTableEnv->fileName     = NULL             ;
        ErrTableEnv->sourceFile   = sourcefile       ;
        ErrTableEnv->allCaps      = '\0'             ;
        ErrTableEnv->errCodePrfx  = '\0'             ;
        ErrTableEnv->dfltFile     = '\0'             ;
        ErrTableEnv->suffix       = ".tex\0"         ;
        ErrTableEnv->cnvtnVioltn  = 0                ;
        ErrTableEnv->lineCount    = 0                ;
        ErrTableEnv->Preamble     = TopOfErrTableEnv ;
        ErrTableEnv->PostFix      = EndOfErrTableEnv ;


        return 0;
}

/*
 *
 * Routine that writes  the Error Code Table:
 *
 */
int WriteErrTable( char *line , LALEnvironment *Env )
{
        fpos_t  linePtr;
        char *onFlagOnCurrentLine;
        onFlagOnCurrentLine = strstr(line, Env->OnFlag );

        if( onFlagOnCurrentLine  ){

           /*
            * This routine determines that the Error Table On Flag
            * has been set. It then searches for the file and
            * writes the table.  Then it repositions the file
            * pointer back to the line just below the
            * Error Table On flag as if nothing has happened.
           */
           Env->cnvtnVioltn = 0;  /* initialize the convention violation flag*/
           /* Get the file position at the top of the environment, ... */
           fgetpos( Env->InFilePtr , &linePtr);
           /* ... find and print the entire table, ... */
           LALDocConstructErrorTable( Env );
           /* and reset the file  position as if nothing has happened. */
           fsetpos( Env->InFilePtr , &linePtr);

        }
        else
        {
                /* if we are in the in body table, skip it we have
                 * already done the work above
                */
        }


        return 0;
}

void *
TopOfErrTableEnv(void *recastToEnv )
{
    LALEnvironment *Env;
    Env = (LALEnvironment *)recastToEnv;

    /* put the header in for the table  */
    fprintf(Env->OutFilePtr,"\\begin{center}");
    fprintf(Env->OutFilePtr,"\\begin{tabular}{|ccp{4.5in}|}");
    fprintf(Env->OutFilePtr,"\\hline");
    fprintf(Env->OutFilePtr,"\n \\texttt{<}\\textit{name}\\texttt{>} & code & description \\\\ ");
    fprintf(Env->OutFilePtr,"\\hline");

    return recastToEnv;
}

/*
 *
 * Routine that closes the table  environment and writes the caption:
 *
*/
void *
EndOfErrTableEnv(void  *recastToEnv )
{
    LALEnvironment *Env;
    Env = (LALEnvironment *)recastToEnv;

    /* close out the table */
    fprintf(Env->OutFilePtr,"\\hline");
    fprintf(Env->OutFilePtr,"\n\\end{tabular}");
    fprintf(Env->OutFilePtr,"\n\\end{center}");

    fprintf(Env->OutFilePtr,"\n\\vspace{-.1in}");

    if( Env->cnvtnVioltn ) {
    fprintf(Env->OutFilePtr,"The status codes in the table above did not obey  ");
    fprintf(Env->OutFilePtr, "the LAL naming convention, i.e. \n");
    fprintf(Env->OutFilePtr, "the code does not use the file name in all caps ");
    fprintf(Env->OutFilePtr, "as the prefix for the error names.\n");
    fprintf(Env->OutFilePtr, "Consult the source code for the full names.\n");
    fprintf(Env->OutFilePtr, "Better yet, fix the code.\n");
    }
    else{
    fprintf(Env->OutFilePtr,"The status codes in the table above are stored in ");
    fprintf(Env->OutFilePtr, "the constants\n");
    fprintf(Env->OutFilePtr, "\\verb@%s_E@\\texttt{<}\\textit{name}\\texttt{>},\n",
          Env->errCodePrfx );
    fprintf(Env->OutFilePtr, "and the status descriptions in\n");
    fprintf(Env->OutFilePtr, "\\verb@%s_MSGE@\\texttt{<}\\textit{name}\\texttt{>}. ",
          Env->errCodePrfx );

    }
    fprintf(Env->OutFilePtr,"The source code ");
    fprintf(Env->OutFilePtr,"with these messages is in \\verb@%s@",
          Env->sourceFile );
    fprintf(Env->OutFilePtr, " on line \\verb@l.%i@.", lineNum );
    return recastToEnv;
}
