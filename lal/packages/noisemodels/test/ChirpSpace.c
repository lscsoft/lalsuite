/* 
   This code generates the chirp parameter space for a given
   minimum companion mass mMin and maximum total mass MMax. One can use
   xmgr to plot the resulting file to plot as in the demo script.
   August 8 , 00.
*/
#include <stdio.h>
#include <lal/LALInspiral.h>

INT4 lalDebugLevel=0;

int main ( void ) {
   static InspiralTemplate p;
   static LALStatus status;
   double mmin, Mmax, totalMmax, compmmin, m2;

/**************************************************/
/* Change the parameters of the search space here */
/**************************************************/
   mmin = 3.0;
   Mmax = 20.;
   mmin = log10(mmin);
   Mmax = log10(Mmax);

   p.ieta=1; 
   p.fLower=40.0; 
   p.order = twoPN;

   totalMmax = pow(10.,Mmax);
   compmmin = pow(10.,mmin);
   fprintf(stderr, "mmin=%e Mmax=%e\n", pow(10., mmin), pow(10., Mmax));

   p.massChoice=m1Andm2;
   p.mass1 = compmmin;
   for (m2=mmin; m2<=Mmax; m2+=0.01) {
      p.mass2 = pow(10.,m2);
      LALInspiralParameterCalc (&status, &p);
      if (p.totalMass > totalMmax) break;
      printf("%e %e %e %e %e %e %e %e %e %e %e\n", 
         p.t0,
         p.t3,
         p.t2,
         p.mass2,
         p.mass1, 
         p.t4,
         p.totalMass,
         p.eta,
         p.mu,
         p.chirpMass,
         p.tC);
   }
   printf("&\n");

   p.totalMass = totalMmax;
   for (m2=log10(totalMmax-compmmin); m2>=mmin; m2-=0.01) {
      p.mass2 = pow(10.,m2);
      if ((p.mass1=p.totalMass - p.mass2) > p.totalMass/2.) break;
      LALInspiralParameterCalc (&status, &p);
      printf("%e %e %e %e %e %e %e %e %e %e %e\n", 
         p.t0,
         p.t3,
         p.t2,
         p.mass2,
         p.mass1, 
         p.t4,
         p.totalMass,
         p.eta,
         p.mu,
         p.chirpMass,
         p.tC);
   }
   printf("&\n");
   p.massChoice=totalMassAndEta;
   p.eta = 0.25;
   for (m2=log10(totalMmax); m2>=mmin; m2-=0.01) {
      if ((p.totalMass = pow(10.,m2)) < 2.*compmmin) break;
      LALInspiralParameterCalc (&status, &p);
      printf("%e %e %e %e %e %e %e %e %e %e %e\n", 
         p.t0,
         p.t3,
         p.t2,
         p.mass2,
         p.mass1, 
         p.t4,
         p.totalMass,
         p.eta,
         p.mu,
         p.chirpMass,
         p.tC);
   }
   return 0;
}
