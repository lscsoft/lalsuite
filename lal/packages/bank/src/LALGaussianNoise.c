#include <stdlib.h>
#include <lal/LALInspiralBank.h>

#define random() rand()
#define srandom( seed ) srand( seed )

NRCSID (LALGAUSSIANNOISEC, "$Id$");
void 
LALGaussianNoise (
   LALStatus   *status,
   REAL4Vector *noisy, 
   UINT8       *seed) 
{
   int i;

   INITSTATUS (status, "LALGaussianNoise", LALGAUSSIANNOISEC);
   ATTATCHSTATUSPTR(status);
   srandom(*seed);
   *seed = random();
   for (i=0; i<(int)noisy->length; i++){
      LALGaussianRandomNumber(status->statusPtr,noisy->data+i);
   }
   CHECKSTATUSPTR(status);
   DETATCHSTATUSPTR(status);
   RETURN(status); 
}

