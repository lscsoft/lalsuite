/*
*  Copyright (C) 2007 Bernd Machenschalk, David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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

/* <lalVerbatim file="ChirpSpaceCV">
Author: Sathyaprakash, B. S., Cokelaer, T.
$Id$
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

#include <lal/LALRCSID.h>
NRCSID (CHIRPSPACEC,"$Id$");

INT4 lalDebugLevel=0;

int main ( void ) {
   static InspiralTemplate p;
   static LALStatus status;
   double mmin, mmax, Mmax, totalMmax, compmmin, m1, m2, finalmass;
   UINT2 type;
   FILE *fpr;

   fpr = fopen("ChirSpace.out", "w");
   /*
    Change the parameters of the search space here

    type=0 creates a region defined by mMin and mMax
    i.e. maximum mass of the companion given by mMax
    type=1 creates a region defined by mMin and MMax
    i.e. total mass of the body given by MMax
   */
   type = 0;
   mmin = 1.;
   mmax = 20.;
   Mmax = mmax*2.;
   p.ieta=1;
   p.fLower=40.0;
   /*
    Don't change anything below:
   */
   mmin = log10(mmin);
   mmax = log10(mmax);
   Mmax = log10(Mmax);

   p.order = LAL_PNORDER_TWO;

   totalMmax = pow(10.,Mmax);
   compmmin = pow(10.,mmin);
   fprintf(fpr, "#mmin=%e Mmax=%e\n", pow(10., mmin), pow(10., Mmax));

   p.massChoice=m1Andm2;
   p.mass1 = compmmin;

   if (type) finalmass=Mmax; else finalmass=mmax;

   for (m2=mmin; m2<=finalmass; m2+=0.01) {
      p.mass2 = pow(10.,m2);
      LALInspiralParameterCalc (&status, &p);
      if (p.totalMass > totalMmax) break;
      fprintf(fpr, "%e %e %e %e %e %e %e %e %e %e %e\n",
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
   fprintf(fpr, "&\n");

   if (type)
   {
	   p.totalMass = totalMmax;
	   for (m2=log10(totalMmax-compmmin); m2>=mmin; m2-=0.01)
	   {
		   p.mass2 = pow(10.,m2);
		   if ((p.mass1=p.totalMass - p.mass2) > p.totalMass/2.) break;
		   LALInspiralParameterCalc (&status, &p);
		   fprintf(fpr, "%e %e %e %e %e %e %e %e %e %e %e\n",
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
		   fprintf(fpr, "%e %e %e %e %e %e %e %e %e %e %e\n",
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
   fprintf(fpr, "&\n");
   p.massChoice=totalMassAndEta;
   p.eta = 0.25;
   for (m2=log10(totalMmax); m2>=mmin; m2-=0.01) {
      if ((p.totalMass = pow(10.,m2)) < 2.*compmmin) break;
      LALInspiralParameterCalc (&status, &p);
      fprintf(fpr, "%e %e %e %e %e %e %e %e %e %e %e\n",
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


   fclose(fpr);
   return 0;
}
