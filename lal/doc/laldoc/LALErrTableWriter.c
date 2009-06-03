/*
*  Copyright (C) 2007 Alan Wiseman, Jolien Creighton, Reinhard Prix
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

int
ParseErrLine(char *Ptr , LALEnvironment *Env,
             char *caps , char *errName , char *errStr , char *errNum );


/*
 * Module Contains the routine that sorts the error code information
 * and generates a nice LaTeX table.
 *
 * The table information in the code should be bracketed
 * by delimiters.  The entries must obey the file naming
 * convention Crap1.H --> CRAP1H_EMYMESSAGE, etc.
 * The entries may be grouped any which way.   The routine will
 * abort if it determines that for some Error Code there
 * does not exist a corresponding error string.
 *
 * An Example:
*/

/*<lalErrTable file="MyCodes">

   The program will step over extraneous crap int he middle
   of the environment.

   #define CRAP1C_ENUL  11
   #define CRAP1C_EOUT  234

   You can rearragne these entries to your hearts content.

   #define CRAP1C_MSGENUL  "Null pointer"
   #define CRAP1C_MSGEOUT  "Output already exists"


  </lalErrTable> */
/*
 *
 *
*/

int
LALDocConstructErrorTable(  LALEnvironment *Env  )
{
    char *linePtr , *errNumPrint , *errStrPrint;
    char caps1[MAXSTR],errName1[MAXSTR],errStr1[MAXSTR],errNum1[MAXSTR];
    char caps2[MAXSTR],errName2[MAXSTR],errStr2[MAXSTR],errNum2[MAXSTR],line[MAXSTR];
    int  j,szHASH,numberOfHashDefines,numLinesInTable,numPairsFound,FoundTheEnd;
    fpos_t  filePositionX,filePositionY,filePositionAtArrival;

    fgetpos( Env->InFilePtr , &filePositionAtArrival);
    szHASH = strlen(HASHDEFINE);


    /* Determine the number of lines of text in this table environment */
    fgetpos( Env->InFilePtr , &filePositionX );
    FoundTheEnd = 0;
    Env->lineCount = 1;
    while (!FoundTheEnd){
        if( !fgets(line,MAXSTR,Env->InFilePtr) ){
                LALDocErr("Table Environment has no end.",
                           Env->sourceFile , __LINE__ , __FILE__ , 1 );
        }
        Env->lineCount ++ ;
        if(strstr(line,Env->OffFlag)) FoundTheEnd=1;
    }
    fsetpos( Env->InFilePtr , &filePositionX);

    /* Determine the number of entries in the table, ... */
    fgetpos( Env->InFilePtr , &filePositionX );
    numberOfHashDefines=0;
    j = 1;
    while( fgets(line,MAXSTR,Env->InFilePtr) && (j<Env->lineCount) ){
            if(strstr(line,HASHDEFINE)) numberOfHashDefines++;
            j++;
    }
    /* and make sure it is even number Hash Defines. */
    numLinesInTable = numberOfHashDefines/2;
    if( numberOfHashDefines % 2 ){
        LALDocErr("Table has odd number of errors and messages.",
        Env->sourceFile , __LINE__ , __FILE__ , 1 );
    }
    fsetpos( Env->InFilePtr , &filePositionX);



    /*
     *
     * Loop over the lines in the Table environment
    */
    fgetpos( Env->InFilePtr , &filePositionY );
    numPairsFound = 0;
    while (numPairsFound < numLinesInTable){

        filePositionX   =   filePositionY;
        fsetpos( Env->InFilePtr , &filePositionX);

        /* Find one of the entries that make up a line in the table */
        for(j=0;j<MAXSTR;j++) line[j] = '\0';
        fgets( line , MAXSTR , Env->InFilePtr ) ;
        while ( !strstr(line , HASHDEFINE ) ){
            fgets( line , MAXSTR , Env->InFilePtr ) ;
            if( strstr(line , Env->OffFlag) ){
                LALDocErr("Incomplete Error Table.",
                Env->sourceFile , __LINE__ , __FILE__  , 1 );
            }
        }
        fgetpos( Env->InFilePtr , &filePositionY );

        /* inner loop over remaining entries in the Table env
         * to find the other piece of the pair. */
        linePtr=strstr(line,HASHDEFINE) + szHASH;
        ParseErrLine(linePtr,Env,caps1,errName1,errStr1,errNum1);
        while( (strcmp(errName1,errName2)!= 0) && !strstr(line,Env->OffFlag) ){
             fgets( line , MAXSTR , Env->InFilePtr )  ;
             linePtr=strstr(line,HASHDEFINE) ;
             if(linePtr){
                 linePtr += szHASH ;
                 ParseErrLine(linePtr,Env,caps2,errName2,errStr2,errNum2);
             }
             if( strcmp(errName1,errName2)==0 ){
                 numPairsFound++;
                 errNumPrint = errNum2;
                 errStrPrint = errStr1;
                 if ( strlen(errNum2) == 0 ) {
                     errNumPrint = errNum1 ;
                     errStrPrint = errStr2 ;
                 }
                 fprintf(Env->OutFilePtr,
                        "\\tt %s & %s & \\vspace{-1.4ex}\\tt %s  \\\\ \n",
                         errName1 , errNumPrint, errStrPrint);
             }
        } /* end of inner loop over remaining entries */
    } /* end loop over lines in the Table environment */


    /* restore to same position with in the input file */
    fsetpos( Env->InFilePtr , &filePositionAtArrival);
    return 0;
}

int
ParseErrLine( char *Ptr          , /* Ptr to string after #define)      */
              LALEnvironment *Env, /* pass the Env variables            */
              char *caps         , /* returns MYFILEH                   */
              char *errName      , /* returns EYOUSUCK, or MSGEYOUSUCK  */
              char *errStr       , /* returns string "you really suck"  */
              char *errNum       ) /* Error Number                      */
{
    char *position, *savePosition;
    char line[MAXSTR];
    int i,ni;

    /* Initialize the strings */
    for(i=0;i<MAXSTR;i++){
            *(errName+i) = '\0';
            *(errStr +i) = '\0';
            *(errNum +i) = '\0';
            *(caps   +i) = '\0';
    }


    /*
     * Find the all caps file name
    */
    sscanf(Ptr,"%s", caps );
    position = strstr(caps,"_");
    *position = '\0';
    sscanf(caps,"%s", Env->errCodePrfx);

    /*
     * If there is ever an error code code convention violation flag it.
     * [Flag it the first time, but ignore it after that.]
    */
    if( !Env->cnvtnVioltn &&  strcmp(Env->errCodePrfx, Env->allCaps) ) {
            Env->cnvtnVioltn = 1;
            LALDocErr("Violation of LAL Error Table Naming Convention.",
                       Env->sourceFile , __LINE__ , __FILE__ , 0 );
    }

    /*
     * Find the Error Name and Descripiton
     * Look for E's and M's, stop if they aren't there.
    */
    position = strstr(Ptr,"_");
    switch ( *(position+1)  ){
            case 'E' :  { sscanf(position+2,"%s",errName);
                          break;
                        }
            case 'M' :  { sscanf(position+5,"%s",errName);
                          break;
                        }
            default  :  { LALDocErr("Bad Error Table: Mind your E's and MSGE's.",
                          Env->sourceFile , __LINE__ , __FILE__ , 1 );
                          break;
                        }
    }

    position += strcspn(position, " \t");
    savePosition = position ;
    sscanf(position,"%s",errStr);

    /*
     * Read the Error Number, ...
    */
    if ( isdigit( (int)errStr[0] ) ) {
         sscanf(position,"%s",errNum);
         *errStr = '\0' ;
    }
    /*
     * ... or read the error description string.
    */
    else
    {
        *errNum = '\0' ;
        position = strstr(Ptr,"\"" );
        if( errStr[0] == '"' ){
            i=1;
            ni=0;
            while( position[i-ni] != '"' ){
                /* if there is a line break in an error explanation */
                if( (position[i-ni] =='\n')  ){
                     LalDocGetLine(line,MAXSTR, Env->InFilePtr );
                     position = line;
                     errStr[i] = ' ';
                     i++;
                     ni = i;
                }
                errStr[i] = *(position+i-ni) ;
                i++;
            }
            errStr[i]='"';
        }
    }
    return 0;
}
