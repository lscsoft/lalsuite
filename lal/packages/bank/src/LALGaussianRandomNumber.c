#include <stdlib.h>
#include <lal/LALInspiralBank.h>

long int random(void);
void srandom(unsigned int seed);

NRCSID (LALGAUSSIANRANDOMNUMBERC, "$Id$");
void
LALGaussianRandomNumber(
   LALStatus *status, 
   REAL4 *randnum)
{
   static int iset=0;
   static REAL8 gset;
   REAL8 fac,rsq,v1,v2;

   INITSTATUS (status, "LALGaussianRandomNumber", LALGAUSSIANRANDOMNUMBERC);
   if  (iset == 0) {
      do {
         v1 = 2.*(random()/(REAL8)RAND_MAX) - 1.0;
         v2 = 2.*(random()/(REAL8)RAND_MAX) - 1.0;
         rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      *randnum = (REAL4) (v2*fac);
      RETURN(status); 
   } else {
      iset=0;
      *randnum = (REAL4) (gset);
      RETURN(status); 
   }
}
