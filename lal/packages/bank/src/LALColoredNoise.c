#include <lal/LALInspiralBank.h>

NRCSID (LALCOLOREDNOISEC, "$Id$");
void 
LALColoredNoise (
   LALStatus   *status,
   REAL4Vector *noisy, 
   REAL8Vector  psd) 
{

   int i, j;
   REAL8 x, length;

   INITSTATUS (status, "LALColoredNoise", LALCOLOREDNOISEC);
   if (psd.length != noisy->length/2+1) {
     fprintf(stderr, "LALColoredNoise: Incompatible lengths\n");
     exit(1);
   }
   length = noisy->length;
   *(noisy->data+0) = 0;
   *(noisy->data+1) = 0;

   for (i=2; i<(int)psd.length; i++) {
      j = 2*i-2;
      x = sqrt(*(psd.data+i) / length);
      *(noisy->data+j) *= x;
      *(noisy->data+j+1) *= x;
   }
   RETURN(status); 
}
