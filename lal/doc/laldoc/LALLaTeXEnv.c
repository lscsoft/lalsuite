#include "LALDoc.h"
/*
 *
 * Function that sets flag for LaTeX Verbatime environment.
 *
 */

int SetLaTeXFlags( LALEnvironment *LaTeXEnv, char *sourcefile , FILE *inSrc )
{
        LaTeXEnv->OnFlag     = "<lalLaTeX"   ;
        LaTeXEnv->OffFlag    = "</lalLaTeX"  ;
        LaTeXEnv->closer     = ">"           ;
        LaTeXEnv->On         = 0             ;
        LaTeXEnv->InFilePtr  = inSrc         ;
        LaTeXEnv->OutFilePtr = NULL          ;
        LaTeXEnv->sourceFile = sourcefile    ;
        LaTeXEnv->suffix     = ".tex\0"      ;
        LaTeXEnv->Preamble   = TopOfLaTeXEnv ;
        LaTeXEnv->PostFix    = EndOfLaTeXEnv ;
    
        return 0;
}

/*
 * 
 * Routine that writes to the  LaTeX file 
 *
 */
int WriteLaTeX( char *line , LALEnvironment *Env )
{  
        char *onFlagOnCurrentLine;
        onFlagOnCurrentLine = strstr(line, Env->OnFlag );

        if( onFlagOnCurrentLine  ){ 
             /* We do not print anything from line with On-Flag.*/
             return 0;
        }
        else
        {
             /* if we are in the in body of the LaTeX,  we print the line */
             fputs(line,Env->OutFilePtr );
        }
        return 0;
}


void *
TopOfLaTeXEnv(void *recastToEnv)
{
    LALEnvironment *Env;
    Env = (LALEnvironment *)recastToEnv;
    /* nothing is needed in the latex file, ie no \begin{...} */
    return recastToEnv;
}

void *
EndOfLaTeXEnv(void *recastToEnv)
{
    LALEnvironment *Env;
    Env = (LALEnvironment *)recastToEnv;
    /* nothing is needed in the latex file, ie no \end{...} */
    return recastToEnv;
}

