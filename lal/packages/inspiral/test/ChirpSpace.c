/* <lalVerbatim file="ChirpSpaceCV">
Author: Sathyaprakash, B. S.
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Program \texttt{ChirpSpace.c}}
\label{ss:bank:ChirpSpace.c}

Test code for \texttt{LALInspiralParameterCalc} module.
If the variable \texttt{type} is set to 1 the code works 
out the boundary of the region enclosed
by the parameter space specified by {\em maximum total
mass} and {\em minimum companion masses,} as given in \texttt{mmin} 
and \texttt{Mmax}. If the variable \texttt{type} is set to 0  
it computes the boundary of the region corresponding to the companion
masses in the range defined by \texttt{mmin} and \texttt{mmax}.

\subsubsection*{Uses}
\begin{verbatim}
LALInspiralParameterCalc
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ChirpSpaceCV}}
</lalLaTeX> */
/* 
   This code generates the chirp parameter space for a given
   minimum companion mass mMin and maximum total mass MMax. One can use
   xmgr to plot the resulting file to plot as in the demo script.
   August 8 , 00.
*/
#include <stdio.h>
#include <lal/LALInspiral.h>

INT4 lalDebugLevel=4;

int main ( void ) {
   static InspiralTemplate p;
   static LALStatus status;
   double mmin, mmax, Mmax, totalMmax, compmmin, m1, m2, finalmass;
   UINT2 type;

   /*
    Change the parameters of the search space here 
   
    type=0 creates a region defined by mMin and mMax
    i.e. maximum mass of the companion given by mMax
    type=1 creates a region defined by mMin and MMax
    i.e. total mass of the body given by MMax
   */
   type = 0;
   mmin = 3.0;
   mmax = 30.;
   Mmax = mmax*2.;
   /*
    Don't change anything below: 
   */
   mmin = log10(mmin);
   mmax = log10(mmax);
   Mmax = log10(Mmax);

   p.ieta=1; 
   p.fLower=20.0; 
   p.order = twoPN;

   totalMmax = pow(10.,Mmax);
   compmmin = pow(10.,mmin);
   fprintf(stderr, "mmin=%e Mmax=%e\n", pow(10., mmin), pow(10., Mmax));

   p.massChoice=m1Andm2;
   p.mass1 = compmmin;

   if (type) finalmass=Mmax; else finalmass=mmax;

   for (m2=mmin; m2<=finalmass; m2+=0.01) {
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

   if (type)
   {
	   p.totalMass = totalMmax;
	   for (m2=log10(totalMmax-compmmin); m2>=mmin; m2-=0.01) 
	   {
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
   }
   else
   {
	   p.totalMass = totalMmax;
	   p.mass2 = p.totalMass/2.;
	   for (m1=mmin; m1<=mmax; m1+=0.01) 
	   {
		   p.mass1 = pow(10.L,m1);
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
