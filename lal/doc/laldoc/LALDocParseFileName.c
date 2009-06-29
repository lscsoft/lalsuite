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
 *   Routine for parsing a line for the file name or assigning a
 *   default.
 *
 *   Exits via LALDocErr() if there is bad syntax on the line
 *   Returns the output file in Env->fileName if one is found.
 *   Returns Env->fileName = NULL if there is no file name.
 *
 * TO DO: The original parser that patrick wrote was better
 * than this one.
*/

int
FindNameOfOutPutFile( char *line , LALEnvironment *Env )
{
     int sizech,lOfFileName,i;
     char *shift,*ptrCloser,*ptrPHILE,*mark,*ptrFileStr;
     char li[MAXSTR];
     sizech = sizeof(char);

     for(i=0;i<MAXSTR;i++)  (*(li+i)) = '\0'; ;
     for(i=0;i<MAXSTR;i++)  (*(li+i)) = (*(line+i)) ;

     /* Shift string pointer over to the end of the environment OnFlag */
     shift = strstr(li,Env->OnFlag) + sizech*strlen( Env->OnFlag ) ;

     /* If there no closer on line to the right of the flag ... eject   */
     ptrCloser = strstr(shift , Env->closer);
     if ( !ptrCloser ) {
              LALDocErr("No parsing closer on the line.",
                        Env->sourceFile , __LINE__ , __FILE__ , 1 );
     }

     /* pointer to PHILE (ie "file=") string the input line */
     ptrPHILE  = strstr(shift , PHILE );

     /* possibly no file requested */
     if ( ptrPHILE==NULL || (ptrPHILE > ptrCloser)  ) {
             mark = shift ;
             /* check for no extraneous junk on the line */
             while(  (*(mark) == ' ')  ) mark =mark + sizech ;
             if ( mark != ptrCloser ){
                 LALDocErr("Bad syntax in source. Extraneous junk on the line.",
                           Env->sourceFile, __LINE__ , __FILE__ , 1);
             }
             /* line is good, but no file specified, so return a blank string. */
             *(Env->fileName) = '\0'; ;
             return 1;
     }
     else /* dig for the file name */
     {
             /* Squeeze in on the file name ... */

             /* ... from the left (take out leading spaces and quotes) */
             ptrFileStr = ptrPHILE + sizech*strlen( PHILE );
             while((*(ptrFileStr) == ' ')||(*(ptrFileStr)== '\"'))
                    { ptrFileStr = ptrFileStr + sizech ; }

             /* from the right  (take out trailing spaces and quotes) */
             mark = ptrCloser;
             while( (*(mark-sizech)==' ')||(*(mark-sizech)=='\"')||
                    (*(mark-sizech)=='\n') )mark=mark-sizech;
             if(ptrFileStr > (mark-sizech) ){
                LALDocErr("Bad syntax in source. Output file needs to specified.",
                          Env->sourceFile , __LINE__ , __FILE__ , 1 );
             }

             /* assign the file name on the line to the Env */
             lOfFileName = strlen(ptrFileStr) - strlen(mark);
             *(ptrFileStr+ lOfFileName) = '\0';
             for(i=0;i<lOfFileName;i++) {
                     Env->fileName[i] = ptrFileStr[i];
             }
             strcat(Env->fileName,Env->suffix);
             if (  strstr(Env->fileName," ")  ) {
               LALDocErr("Spaces in the file name.",
                         Env->sourceFile , __LINE__ , __FILE__  , 1 );
             }
             return 1;
     }
}


/*
 *
 * Returns default file name, ie crap.H is the input -> crapH.tex
 * is the default output file result.
 *
*/


int
FindDefaultFileName( LALEnvironment *Env )
{
        char *period;
        char c;
        int  i, lensrc , lendot , lensuffix, lenBaseNameStr ;

        period = strstr( Env->sourceFile , "." ) ;
        if( !period ){
                LALDocErr("No dot in in input source file name.",
                        Env->sourceFile , __LINE__ , __FILE__ , 1 );
        }
        lensrc = strlen(Env->sourceFile);
        lendot = strlen(period);
        lenBaseNameStr = lensrc-lendot;

        /* put the base (input) file name into the default (output) file name */
        for (i=0;i < lenBaseNameStr ; i++) {
                c = *(Env->sourceFile +i);
                *(Env->dfltFile + i ) = c;
                *(Env->allCaps  + i ) = (char)toupper( (int)c );
        }

        /* tack on the file Extention (.c or .h) shift to upper case */
        *(Env->dfltFile + lenBaseNameStr)   =  toupper(*(period+1));
        *(Env->allCaps  + lenBaseNameStr)   =  toupper(*(period+1));

        /* and put on the output file extention (.tex) */
        lensuffix = strlen(Env->suffix);
        i= 0 ;
        while(i < lensuffix ) {
              *(Env->dfltFile +lenBaseNameStr + 1  +i  )  = *(Env->suffix+i);
              i++;
        }


        return 1;
}


