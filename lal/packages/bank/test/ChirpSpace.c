/* 
   This code generates the chirp parameter space for a given
   mMin and MMax as command line arguments. One can use
   xmgr to plot the resulting file to plot as in the demo script.

   August 8 , 00.
*/
#include <stdio.h>
#include <lal/LALInspiral.h>

INT4 lalDebugLevel=1;

int main (int argc, char *argv[]) {
#if FIXME
   static InspiralTemplate p;
   static LALStatus status;
   double mmin, Mmax, totalMmax, compmmin, m2;

   mmin = log10(atof(argv[1]));
   Mmax = log10(atof(argv[2]));
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
         p.t2,
         p.t3,
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
         p.t2,
         p.t3,
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
         p.t2,
         p.t3,
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
#else
   argc = 0; argv = NULL; return 77;
#endif
}
