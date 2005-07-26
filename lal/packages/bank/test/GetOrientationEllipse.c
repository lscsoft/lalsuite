/* <lalVerbatim file="SpaceCoveringCV">
Author: Thomas Cokelaer 
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Program \texttt{SpaceCovering.c}}
\label{ss:SpaceCovering.c}

Test code for the \texttt{bank} modules.

\subsubsection*{Usage}
\begin{verbatim}
SpaceCovering --template [TaylorT1, EOB ....] --grid-type [square, hexagonal, squareOriented, hexagonalOriented]
\end{verbatim}

\subsubsection*{Description}

This test code gives an example of how to generate a template bank and
generates vertices of the ambiguity 'rectangle' around each lattice point
suitable for plotting with xmgr or xgrace.

\subsubsection*{Exit codes}
\input{CoarseTest2CE}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALRectangleVertices
LALInspiralCreateCoarseBank
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CoarseTest2CV}}
</lalLaTeX> */

/* <lalErrTable file="CoarseTest2CE"> */
/* </lalErrTable> */

#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>


void LALInspiralCreateBoundarySpace(InspiralCoarseBankIn coarseIn);

typedef struct {
  int calque;
   int gridSpacing;
}
UserParams; 

INT4 lalDebugLevel=33;

int
main(int argc, char **argv)
{
  INT4 arg;
  /* top-level status structure */
  static LALStatus status;     
  /* Structure specifying the nature of the bank needed */
  static InspiralCoarseBankIn coarseIn;



  /*  void *noisemodel = LALLIGOIPsd;*/
  void *noisemodel = LALVIRGOPsd;
  UINT4   j, numPSDpts=262144/4/4;
  InspiralTemplate tempPars;
  InspiralMetric metric;
  InspiralMomentsEtc moments;


  coarseIn.LowGM        = -2;
  coarseIn.HighGM       = 6;
  coarseIn.fLower       = 20.L;
  coarseIn.fUpper       = 2047.L;
  coarseIn.tSampling    = 4096.L;
  coarseIn.order        = twoPN;
  coarseIn.space        = Tau0Tau3;
  coarseIn.mmCoarse     = 0.97;
  coarseIn.mmFine       = 0.97;
  coarseIn.iflso        = 0.0L;
  coarseIn.mMin         = 3;
  coarseIn.mMax         = 30.0;
  coarseIn.MMax         = coarseIn.mMax * 2.;
  coarseIn.massRange    = MinMaxComponentMass; 
  /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/
  /* minimum value of eta */
  coarseIn.etamin       = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);
  coarseIn.psi0Min      = 1.e0;
  coarseIn.psi0Max      = 2.5e4;
  coarseIn.psi3Min      = -1e4;
  coarseIn.psi3Max      = -10;
  coarseIn.alpha        = 0.L;
  coarseIn.numFcutTemplates = 4;

  memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
  coarseIn.shf.f0 = 0;
  LALDCreateVector( &status, &(coarseIn.shf.data), numPSDpts );
  coarseIn.shf.deltaF = coarseIn.tSampling / (2.*(REAL8) coarseIn.shf.data->length + 1.L);
  LALNoiseSpectralDensity (&status,
			   coarseIn.shf.data,
			   noisemodel, coarseIn.shf.deltaF );



  tempPars.totalMass = coarseIn.MMax;
  tempPars.eta = 0.25;
  tempPars.ieta = 1.L;
  tempPars.fLower = coarseIn.fLower;
  tempPars.massChoice = totalMassAndEta; 

  metric.space = coarseIn.space;
  LALInspiralSetParams( &status, &tempPars, coarseIn );


  {
    double tau0, tau3;
    int validPars;
    InspiralBankParams bankPars, bankParsOld;

    LALInspiralSetSearchLimits( &status, &bankPars, coarseIn );

    for (tau0=.1; tau0<40; tau0+=2){
      for (tau3=.1; tau3<2; tau3+=0.1){
	tempPars.t0= tau0;
	tempPars.t3=tau3;
	LALInspiralParameterCalc( &status, &tempPars );
	LALGetInspiralMoments( &status, &moments, &coarseIn.shf, &tempPars );
	LALInspiralComputeMetric( &status, &metric, &tempPars, &moments );
	LALInspiralValidTemplate( &status,
				  &validPars, bankPars, coarseIn );

	validPars=1;
	if (validPars) 
	  printf("%f  %f %f\n", tempPars.t0, tempPars.t3, metric.theta*180/3.14159);
	else
	  printf("%f  %f %f\n", tempPars.t0, tempPars.t3, -1.);
	fflush(stdout);

      }
    }
  }


  LALInspiralCreateBoundarySpace(coarseIn);

  /* Free the list, and exit. */

  LALDDestroyVector( &status, &(coarseIn.shf.data) );
  LALCheckMemoryLeaks();
  return(0);
}








void LALInspiralCreateBoundarySpace(InspiralCoarseBankIn coarseIn)
{

  static InspiralTemplate p;
  static LALStatus status;
  
  
  
  double mmin, mmax, Mmax, totalMmax, compmmin, m1, m2, finalmass;
  
  UINT2 type;
  FILE *fpr;
  
  fpr = fopen("ChirpSpace.out", "w");
  
  /*
    Change the parameters of the search space here 
   
    type=0 creates a region defined by mMin and mMax
    i.e. maximum mass of the companion given by mMax
    type=1 creates a region defined by mMin and MMax
    i.e. total mass of the body given by MMax
   */
   type = 0;
   mmin = coarseIn.mMin;
   mmax = coarseIn.mMax;
   Mmax = mmax*2.;
   p.ieta=1; 
   p.fLower=coarseIn.fLower; 
  
   /*
    Don't change anything below: 
   */
   mmin = log10(mmin);
   mmax = log10(mmax);
   Mmax = log10(Mmax);

   p.order = twoPN;

   totalMmax = pow(10.,Mmax);
   compmmin = pow(10.,mmin);
   

   p.massChoice=m1Andm2;
   p.mass1 = compmmin;

   if (type)
     finalmass=Mmax;
   else
     finalmass=mmax;

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

}
