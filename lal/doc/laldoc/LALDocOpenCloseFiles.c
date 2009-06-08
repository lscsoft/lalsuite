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
 * Open a file
 *
*/

FILE
*OpenAFile(const char *file, const char *rwa , int timeStamp  )
{
        FILE *filePtr;
        time_t tp;
        filePtr = NULL ;
        filePtr= fopen((const char *)file,rwa);
        if (timeStamp){
             time( &tp );
             fprintf(filePtr ,"%%Laldoc Opened at: %s",
                     asctime(localtime(&tp)));
        }
        return filePtr;
}

/*
 *
 * Close  a file
 *
*/

int
CloseAFile(FILE *filePtr, int timeStamp )
{
        time_t tp;
        if (timeStamp){
             time( &tp );
             fprintf(filePtr ,"%%Laldoc Closed at: %s\n",
                     asctime(localtime(&tp)));
        }
        fclose((FILE *)filePtr);
   return 0;
}

