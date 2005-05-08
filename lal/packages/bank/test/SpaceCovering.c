/* <lalVerbatim file="CoarseTest2CV">
Author: Churches, D. K. and Sathyaprakash, B. S.
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Program \texttt{CoarseTest2.c}}
\label{ss:CoarseTest2.c}

Test code for the \texttt{bank} modules.

\subsubsection*{Usage}
\begin{verbatim}
CoarseTest2
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
   int gridType;
}
UserParams; 

INT4 lalDebugLevel=0;

int
main(int argc, char **argv)
{
  INT4 arg;
  /* top-level status structure */
  static LALStatus status;     
  /* Structure specifying the nature of the bank needed */
  static InspiralCoarseBankIn coarseIn;
/* Template bank lists */
  static InspiralTemplateList *list1, *list2;
/* Number of templates in list1 and list2 */
  INT4 nlist1=0, nlist2=0;

  void *noisemodel = LALLIGOIPsd;
  UINT4   j, numPSDpts=262144;
  FILE *fpr;
  UserParams userParams;
  INT4 i=1;
  double scaling = sqrt(2.);


  while(i < argc)
    {
      if ( strcmp(argv[i],	"--template") 	== 0 ) {
	if ( strcmp(argv[++i],	"TaylorT1")	==0)	userParams.calque = TaylorT1;
	else if ( strcmp(argv[i],	"TaylorT2")	==0)	userParams.calque = TaylorT2;
	else if ( strcmp(argv[i],	"TaylorT3")	==0)	userParams.calque = TaylorT3;
	else if ( strcmp(argv[i],	"TaylorF1")	==0)	userParams.calque = TaylorF1;
	else if ( strcmp(argv[i],	"TaylorF2")	==0)	userParams.calque = TaylorF2;
	else if ( strcmp(argv[i],	"PadeT1")	==0)	userParams.calque = PadeT1;
	else if ( strcmp(argv[i],	"PadeF1")	==0)	userParams.calque = PadeF1;
	else if ( strcmp(argv[i],	"EOB")		==0)	userParams.calque = EOB;
	else if ( strcmp(argv[i],	"BCV")		==0)    userParams.calque = BCV;
	else if ( strcmp(argv[i],	"SpinTaylorT3")	==0)	userParams.calque = SpinTaylorT3;

      }
      else if  ( strcmp(argv[i],	"--grid-type") 	== 0 ) {
	i++;
	
	if (strcmp(argv[i], "square") == 0)   coarseIn.gridType = BANKGRID_SQUARE;
	else if (strcmp(argv[i], "hexagonal")         == 0)   coarseIn.gridType = BANKGRID_HEXAGONAL;
	else if (strcmp(argv[i], "squareOriented")         == 0)   coarseIn.gridType = BANKGRID_SQUARE_AND_ORIENTED;
	else if (strcmp(argv[i], "hexagonalOriented")         == 0)   coarseIn.gridType = BANKGRID_HEXAGONAL_AND_ORIENTED;
	else {fprintf(stderr, "bank-grid-type is either square or hexagonal\n"); exit(0);}
	
      }
      else if  ( strcmp(argv[i],	"--scaling") 	== 0 ) {
	      i++;
      scaling = sqrt(atof(argv[i]));
      printf("%f", scaling);
      }
      i++;
    }
  fprintf(stderr, "save results in SpaceCovering.out\n");
  fpr = fopen("SpaceCovering.out", "w");

  coarseIn.LowGM = 3;
  coarseIn.HighGM= 6;
  coarseIn.fLower = 40.L;
  coarseIn.fUpper = 2000.L;
  coarseIn.tSampling = 4096.L;
  coarseIn.order = twoPN;
  coarseIn.space = Tau0Tau3;
  coarseIn.mmCoarse = 0.95;
  coarseIn.mmFine = 0.97;
  coarseIn.iflso = 0.0L;
  coarseIn.mMin = 3.0;
  coarseIn.mMax = 20.0;
  coarseIn.MMax = coarseIn.mMax * 2.;
  coarseIn.massRange = MinMaxComponentMass; 
  /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/
  /* minimum value of eta */
  coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);
  coarseIn.psi0Min = 1.e4;
  coarseIn.psi0Max = 2.5e5;
  coarseIn.psi3Min = -3e4;
  coarseIn.psi3Max = -10;
  coarseIn.alpha = 0.L;
  coarseIn.numFcutTemplates = 4;

  memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
  coarseIn.shf.f0 = 0;
  LALDCreateVector( &status, &(coarseIn.shf.data), numPSDpts );
  coarseIn.shf.deltaF = coarseIn.tSampling / (2.*(REAL8) coarseIn.shf.data->length + 1.L);
  LALNoiseSpectralDensity (&status,
			   coarseIn.shf.data,
			   noisemodel, coarseIn.shf.deltaF );

  if (userParams.calque == BCV)
  {
	coarseIn.approximant = BCV;
	coarseIn.space 	= Psi0Psi3;
  }
  else
  { 
	coarseIn.approximant = TaylorT3;
	coarseIn.space 	= Tau0Tau3;
  } 
      	
  LALInspiralCreateCoarseBank(&status, &list1, &nlist1, coarseIn);
  
  for (j=0; j<nlist1; j++)
  {
	  switch(userParams.calque){
	  case TaylorT1:
	  case TaylorT2:
	  case TaylorT3:
	  case EOB:

 	  fprintf(fpr, "%e %e %e %e\n", 
	  		  list1[j].params.t0, 
			  list1[j].params.t3, 
			  list1[j].params.totalMass, 
			  list1[j].params.fFinal);
		  break;
	  case BCV:
	  fprintf(fpr, "%e %e %e %e\n", 
	  		  list1[j].params.psi0, 
			  list1[j].params.psi3, 
			  list1[j].params.totalMass, 
			  list1[j].params.fFinal);
	  break;
	  }
  }
	

	
  fprintf(fpr, "&\n");
  
  {
    UINT4 j;
    UINT4 valid;
  
    static RectangleIn RectIn;
    static RectangleOut RectOut;

  
    RectIn.dx = sqrt(2.0 * (1. - coarseIn.mmCoarse)/list1[0].metric.g00 );
    RectIn.dy = sqrt(2.0 * (1. - coarseIn.mmCoarse)/list1[0].metric.g11 );
    RectIn.theta = list1[0].metric.theta ;
    
    /* Print out the template parameters */
    for (j=0; j<nlist1; j++)
    {
	/*
	Retain only those templates that have meaningful masses:
	*/
	RectIn.x0 = (REAL8) list1[j].params.psi0;
	RectIn.y0 = (REAL8) list1[j].params.psi3;
	/*
	LALInspiralValidParams(&status, &valid, bankParams, coarseIn);
	*/
	valid = 1;
        if (valid) 
	{
		LALRectangleVertices(&status, &RectOut, &RectIn);
		fprintf(fpr, "%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n", 
				RectOut.x1, RectOut.y1, 
				RectOut.x2, RectOut.y2, 
				RectOut.x3, RectOut.y3, 
				RectOut.x4, RectOut.y4, 
				RectOut.x5, RectOut.y5);
		fprintf(fpr, "&\n");
	}

	{
	  int Nth =100;
	  double th;
	  double x,y,phi,a,b, theta;
	  for (theta=0; theta<2*3.14; theta+=2*3.14/(double)Nth)
	    {
	      a = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/list1[0].metric.g00 );
	      b = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/list1[0].metric.g11 );
	      x = a * cos(theta)/scaling;
	      y = b * sin(theta)/scaling;
	      phi=list1[0].metric.theta ;

	      th = x*cos(phi)-y*sin(phi)+list1[j].params.psi0;
	      y  = x*sin(phi)+y*cos(phi)+list1[j].params.psi3;
	      x  = th;
	      fprintf(fpr, "%f %f\n", x, y);	    
	    }
		theta = 0;
	      a = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/list1[0].metric.g00 );
	      b = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/list1[0].metric.g11 );
	      x = a * cos(theta) /scaling;
	      y = b * sin(theta)/scaling;
	      phi=list1[0].metric.theta ;

	      th = x*cos(phi)-y*sin(phi)+list1[j].params.psi0;
	      y  = x*sin(phi)+y*cos(phi)+list1[j].params.psi3;
	      x  = th;
	      fprintf(fpr, "%f %f\n", x, y);	    

	  
	      fprintf(fpr, "&\n");
	}
    }
  }

  LALInspiralCreateBoundarySpace(coarseIn);


  fclose(fpr);
  /* Free the list, and exit. */
  if (list1 != NULL) LALFree (list1);
  if (list2 != NULL) LALFree (list2);
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
