#include "LALDoc.h"
/*
 *
 * Function that sets flag for LaTeX Verbatim environment.
 *
*/
int SetVerbatimFlags( LALEnvironment *VerbatimEnv, char *sourcefile , FILE *inSrc )
{
        VerbatimEnv->OnFlag     = "<lalVerbatim"   ;
        VerbatimEnv->OffFlag    = "</lalVerbatim"  ;
        VerbatimEnv->closer     = ">"              ;
        VerbatimEnv->On         = 0                ;
        VerbatimEnv->InFilePtr  = inSrc            ;
        VerbatimEnv->OutFilePtr = NULL             ;
        VerbatimEnv->fileName   = NULL             ;
        VerbatimEnv->sourceFile = sourcefile       ;
        VerbatimEnv->allCaps    = '\0'             ;
        VerbatimEnv->errCodePrfx= '\0'             ;
        VerbatimEnv->dfltFile   = '\0'             ;
        VerbatimEnv->suffix     = ".tex\0"         ;
        VerbatimEnv->cnvtnVioltn= 0                ;
        VerbatimEnv->lineCount  = 0                ;
        VerbatimEnv->Preamble   = TopOfVerbatimEnv ;
        VerbatimEnv->PostFix    = EndOfVerbatimEnv ;

    
        return 0;
}

/*
 * 
 * Routine that writes each line to the  LaTeX Verbatim  file 
 *
 */
int WriteVerbatim( char *line , LALEnvironment *Env )
{  
        char *onFlagOnCurrentLine;
        onFlagOnCurrentLine = strstr(line, Env->OnFlag );

        if( onFlagOnCurrentLine  ){ 

             return 0;
        }
        else
        {
             /* if we are in the in body of the Verbatim,  we print the line */
	     /* deal with tabs */
	     char tabout[MAXSTR];                        
	     TabAlign(tabout, line, MAXSTR, TABSIZE);
             line = tabout ;
             /* fputs(line,stderr); */
             fputs(line,Env->OutFilePtr );
        }
        return 0;
}

/*
 *
 * Routine that puts the margin par on the verbatim environment:
 *
 */
void *
TopOfVerbatimEnv(void *recastToEnv )
{
    LALEnvironment *Env;
    Env = (LALEnvironment *)recastToEnv;
    /* write the margin par */
    fprintf(Env->OutFilePtr,"\n\\vspace{-1ex}\n");
    fprintf(Env->OutFilePtr,"\\mbox{}\\marginpar{\\tiny\\texttt{l.%i}\\\\\\texttt{%s}}\n",
            lineNum,Env->sourceFile);
    fprintf(Env->OutFilePtr,"\\vspace{-3ex}\n");
    fprintf(Env->OutFilePtr,"\\begin{verbatim}\n");
    return recastToEnv;
}

/*
 *
 * Routine that closes the verbatim environment:
 *
 */
void *
EndOfVerbatimEnv(void  *recastToEnv )
{
    LALEnvironment *Env;
    Env = (LALEnvironment *)recastToEnv;
    fprintf(Env->OutFilePtr,"\\end{verbatim}\n");
    return recastToEnv;
}


/* 
 * Replaces tabs in verbatim code so that it is better formatted:  
 * Teviet Creighton 
 *
*/
void 
TabAlign(char *out, char *in, int n, int tabsize)
     /* This routine writes an input string *in to string *out,
        replacing tab characters with spaces.  Each tab is replaced by
        enough spaces (at least 1) so that the next character's column
        is an integral multiple of tabsize (counting from column 0).
        The parameter n is the maximum length of either array
        (including the '\0' character denoting the end of a string).

	If the tab replacement stretches *out past this length, then
	it will be truncated.  Memory allocation and the determination
	of n are the responsibility of the calling routine.  This
	routine is leak-proof provided it given the right n. */
{
  int i=0,j=0; /* Indecies. */
  char c='x';  /* Initialized to anything nonzero. */

  while((i<n)&&(j<n)&&(c!='\0'))
    if((c=in[i++])=='\t'){
      if(j<n)
	out[j++]=' ';
      while((j%tabsize)&&(j<n))
	out[j++]=' ';
    }
    else
      out[j++]=c;
  return;
}


