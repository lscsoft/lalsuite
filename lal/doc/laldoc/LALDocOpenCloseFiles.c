#include "LALDoc.h"

/*
 *
 * Open a file
 * 
*/

FILE 
*OpenAFile(char *file, char *rwa , int timeStamp  )
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

