#include "LALDoc.h"

/*
 *
 * Function that sets flag for Error Table environment.
 *
 */
int SetErrTableFlags( LALEnvironment *ErrTableEnv, char *sourcefile , FILE *inSrc )
{
        ErrTableEnv->OnFlag     = "<lalErrTable"   ;
        ErrTableEnv->OffFlag    = "</lalErrTable" ;
        ErrTableEnv->closer     = ">"              ;
        ErrTableEnv->On         = 0                ;
        ErrTableEnv->InFilePtr  = inSrc            ;
        ErrTableEnv->OutFilePtr = NULL             ;
        ErrTableEnv->sourceFile = sourcefile       ;
        ErrTableEnv->suffix     = ".tex\0"         ;
        ErrTableEnv->Preamble   = TopOfErrTableEnv ;
        ErrTableEnv->PostFix    = EndOfErrTableEnv ;

      
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

           fgetpos( Env->InFilePtr , &linePtr);
           /* Build and print the entire table ... */
           LALDocConstructErrorTable( Env );
           /* and reset the file  position. */
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
    fprintf(Env->OutFilePtr,"\\begin{tabular}{|lcp{4.5in}|}");
    fprintf(Env->OutFilePtr,"\\hline");
    fprintf(Env->OutFilePtr,"\n\\it name & code & description \\\\ ");
    fprintf(Env->OutFilePtr,"\\hline");

    return recastToEnv;
}

/*
 *
 * Routine that closes the verbatim environment:
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
    fprintf(Env->OutFilePtr,"The status codes in the table above are stored in ");
    fprintf(Env->OutFilePtr, "the constants\n");
    fprintf(Env->OutFilePtr, "\\verb@%s_E@\\texttt{<}\\textit{name}\\texttt{>},\n",
          Env->errCodePrfx );
    fprintf(Env->OutFilePtr, "and the status descriptions in\n");
    fprintf(Env->OutFilePtr, "\\verb@%s_MSGE@\\texttt{<}\\textit{name}\\texttt{>}.",
          Env->errCodePrfx );
    if(strcmp(Env->errCodePrfx,Env->allCaps)){
    fprintf(Env->OutFilePtr,"By convention the status codes are usually stored in ");
    fprintf(Env->OutFilePtr, "the constants\n");
    fprintf(Env->OutFilePtr, "\\verb@%s_E@\\texttt{<}\\textit{name}\\texttt{>},\n",
          Env->allCaps );
    fprintf(Env->OutFilePtr, "and the status descriptions in\n");
    fprintf(Env->OutFilePtr, "\\verb@%s_MSGE@\\texttt{<}\\textit{name}\\texttt{>}. ",
          Env->allCaps );
    fprintf(Env->OutFilePtr,"This file violates the naming convention.");
    }


    return recastToEnv;
}

