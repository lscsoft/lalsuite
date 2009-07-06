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
        LaTeXEnv->fileName   = NULL           ;
        LaTeXEnv->sourceFile = sourcefile    ;
        LaTeXEnv->allCaps    = '\0'          ;
        LaTeXEnv->errCodePrfx= '\0'          ;
        LaTeXEnv->dfltFile   = '\0'          ;
        LaTeXEnv->suffix     = ".tex\0"      ;
        LaTeXEnv->cnvtnVioltn= 0             ;
        LaTeXEnv->lineCount  = 0             ;
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
        int i;
        onFlagOnCurrentLine = strstr(line, Env->OnFlag );

        if( onFlagOnCurrentLine  ){
             /* We do not print anything from line with On-Flag.*/
             return 0;
        }
        else
        {
             /*
              * Elliminate leading *'s from the latex, ie like in
              * this comment block
              *
             */
             i = 0;
             while(*(line+i) == ' ' && *(line+i) != '\0') { i++; }
             if(*(line+i) == '*') { *(line+i)=' '; }
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

