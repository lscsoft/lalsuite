#include <lal/LALInspiralBank.h>

NRCSID (LALADDVECTORSC, "$Id$");

void 
LALAddVectors(
   LALStatus *status, 
   REAL4Vector *vector, 
   AddVectorsIn in)
{
   INT4 i;

   INITSTATUS (status, "LALAddVectors", LALADDVECTORSC);

   i=vector->length;
   if ( (in.v1->length != in.v2->length) ||
        (in.v1->length != vector->length) ) {
      fprintf(stderr, "LALAddVectors: Incompatible lengths\n");
      exit(1);
   }
   while (i--)
      vector->data[i] = in.a1 * in.v1->data[i] + in.a2 * in.v2->data[i];
   RETURN(status);
}

