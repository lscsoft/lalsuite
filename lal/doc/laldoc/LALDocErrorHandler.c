
/*
 * Prints Errors in the search to specified file 
 */

#include "LALDoc.h"

int
LALDocErr(const char *msg, const char *inFileName, int laldocLine ,
          const char *laldocSourceFile, int fatal ) 
{

fprintf(ptrLALErrorFile,"\n     %s    \n \
   Input file: %s.                    \n \
   Input line Number: %i              \n \
   LALDoc Source-file line number: %i \n \
   LALDoc Source-file Name: %s        \n",
   msg, inFileName , lineNum, laldocLine , laldocSourceFile 
);

/*
fprintf(stderr,"\n   %s            \n \
   Input file: %s.                   \n \
   Input line Number: %i             \n \
   LALDoc Source-file line number: %i\n \
   LALDoc Source-file Name: %s       \n",
   msg, inFileName , lineNum, laldocLine , laldocSourceFile 
);
*/

if (fatal){
        CloseAFile( ptrLALErrorFile ,  1 );
        fprintf(ptrLALErrorFile,"Determined to be fatal. Exiting.");
        exit(BADEXIT);
}

return 0;
}

