/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Thomas Cokelaer
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
SpaceCovering
--template [TaylorT1, EOB ....]
--grid-spacing [square, hexagonal, squareOriented, hexagonalOriented]
\end{verbatim}

--template BCV --grid-spacing Hexagonal gives 178 templates
and --template EOB grid-spacing Hexagonal gives 532 templates

\subsubsection*{Description}
This test code gives an example of how to generate a template bank and
generates vertices of the ambiguity 'rectangle' around each lattice point
suitable for plotting with xmgr or xgrace.

This is valid for various template bank placement (Hexagonal, SquareNorOriented).
This is valid for BCV and physical template families.


\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALRectangleVertices
LALInspiralCreateCoarseBank
\end{verbatim}

\subsubsection*{Notes}
\vfill{\footnotesize\input{SpaceCoveringCV}}

</lalLaTeX> */

#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>

#include <lal/LALRCSID.h>
NRCSID (SPACECOVERINGC,"$Id$");


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
  /* Template bank lists */
  static InspiralTemplateList *list1, *list2;
  /* Number of templates in list1 and list2 */
  INT4 nlist1=0, nlist2=0;

  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;
  INT4   j, numPSDpts=16384/2;
  FILE *fpr;
  UserParams userParams;
  INT4 i=1;
  double scaling = sqrt(2.);


  while(i < argc)
    {
      if (strcmp(argv[i], "--help") == 0) {
        fprintf(stderr, "--template [TaylorT1, TaylorT3, EOB, ...]\n--grid-spacing [Hexagonal, SquareNotOriented, ...]\n`--noise-model [LIGO, VIRGO...]\n");
        exit(0);
      }
      else if ( strcmp(argv[i],	"--noise-model") 	== 0 ) {
	if ( strcmp(argv[++i],	"VIRGO")	==0)	{
          (noisemodel) = &LALVIRGOPsd;
        }
        else if ( strcmp(argv[i],	"LIGOI")	==0){
          (noisemodel) = &LALLIGOIPsd;
        }
        else if ( strcmp(argv[i],	"LIGOA")	==0){
          (noisemodel) = &LALAdvLIGOPsd;
        }
        else if ( strcmp(argv[i],	"EGO")	==0){
          (noisemodel) = &LALEGOPsd;
        }
      }
      else if ( strcmp(argv[i],	"--template") 	== 0 ) {
	if ( strcmp(argv[++i],	"TaylorT1")	==0)	{
          userParams.calque = TaylorT1;
        }
	else if ( strcmp(argv[i],	"TaylorT2")	==0)	{
          userParams.calque = TaylorT2;
        }
	else if ( strcmp(argv[i],	"TaylorT3")	==0)	{
          userParams.calque = TaylorT3;
        }
	else if ( strcmp(argv[i],	"TaylorF1")	==0)	{
          userParams.calque = TaylorF1;
        }
	else if ( strcmp(argv[i],	"TaylorF2")	==0)	{
          userParams.calque = TaylorF2;
        }
	else if ( strcmp(argv[i],	"PadeT1")	==0)	{
          userParams.calque = PadeT1;
        }
	else if ( strcmp(argv[i],	"PadeF1")	==0)	{
          userParams.calque = PadeF1;
        }
	else if ( strcmp(argv[i],	"EOB")		==0)	{
          userParams.calque = EOB;
        }
	else if ( strcmp(argv[i],	"BCV")		==0)    {
          userParams.calque = BCV;
        }
	else if ( strcmp(argv[i],	"SpinTaylorT3")	==0)	{
          userParams.calque = SpinTaylorT3;
        }
      } /*end of --template if*/
      else if  ( strcmp(argv[i],	"--grid-spacing") 	== 0 ) {
	i++;
	if (strcmp(argv[i], "Square") == 0)
          coarseIn.gridSpacing = Square;
	else if (strcmp(argv[i], "Hexagonal")         == 0)
          coarseIn.gridSpacing = Hexagonal;
	else if (strcmp(argv[i], "HybridHexagonal")         == 0)
          coarseIn.gridSpacing = HybridHexagonal;
	else if (strcmp(argv[i], "SquareNotOriented")         == 0)
          coarseIn.gridSpacing = SquareNotOriented;
	else if (strcmp(argv[i], "HexagonalNotOriented")         == 0)
          coarseIn.gridSpacing = HexagonalNotOriented;
	else {fprintf(stderr, "grid-spacing is either square or hexagonal\n"); exit(0);}
      }
      i++;
    }

  coarseIn.LowGM        = -2;
  coarseIn.HighGM       = 6;
  coarseIn.fLower       = 40.L;
  coarseIn.fUpper       = 2047L;
  coarseIn.tSampling    = 4096.L;
  coarseIn.order        = LAL_PNORDER_TWO;
  coarseIn.space        = Tau0Tau3;
  coarseIn.mmCoarse     = 0.98;
  coarseIn.mmFine       = 0.98;
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


  fprintf(stderr, "save %d template in results in SpaceCovering.out\n", nlist1);
  fflush(stdout);
  fflush(stderr);

  fpr = fopen("SpaceCovering.out", "w");
  for (j=0; j<nlist1; j++)
  {
	  switch(userParams.calque){
	  case TaylorT1:
	  case TaylorT2:
	  case TaylorT3:
	  case EOB:

 	  fprintf(fpr, "%e %e %e %e %e %e %e\n",
	  		  list1[j].params.t0,
			  list1[j].params.t3,
			  list1[j].metric.g00,
			  list1[j].metric.g11,
                          2*sqrt((1-coarseIn.mmCoarse)/list1[j].metric.g00),
                          2*sqrt((1-coarseIn.mmCoarse)/list1[j].metric.g11),
			  list1[j].metric.theta*180/LAL_PI);
		  break;
	  case BCV:
	  fprintf(fpr, "%e %e %e %e %e %e %e\n",
	  		  list1[j].params.psi0,
			  list1[j].params.psi3,
			  list1[j].metric.g00,
			  list1[j].metric.g11,
                          2*sqrt((1-coarseIn.mmCoarse)/list1[j].metric.g00),
                          2*sqrt((1-coarseIn.mmCoarse)/list1[j].metric.g11),
			  list1[j].metric.theta);
	  break;
	  }
  }



  fprintf(fpr, "&\n");

  {
    UINT4 j;
    UINT4 valid;
    static RectangleIn RectIn;
    static RectangleOut RectOut;
    static HexagonOut HexaOut;

    /* Print out the template parameters */
    for (j=0; j<nlist1; j++)
    {
      RectIn.dx = sqrt(2.0 * (1. - coarseIn.mmCoarse)/list1[j].metric.g00 );
      RectIn.dy = sqrt(2.0 * (1. - coarseIn.mmCoarse)/list1[j].metric.g11 );
      RectIn.theta = list1[j].metric.theta  ;

      if (userParams.calque == BCV)
      {
        RectIn.x0 = (REAL8) list1[j].params.psi0;
        RectIn.y0 = (REAL8) list1[j].params.psi3;
      }
      else
      {
        RectIn.x0 = (REAL8) list1[j].params.t0;
        RectIn.y0 = (REAL8) list1[j].params.t3;
      }
      /*
       * LALInspiralValidParams(&status, &valid, bankParams, coarseIn);
       * */
      valid = 1;
      if (valid)
      {
        if (coarseIn.gridSpacing == SquareNotOriented)
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
        else if (coarseIn.gridSpacing == Hexagonal || coarseIn.gridSpacing == HybridHexagonal)
        {
          RectIn.dx = sqrt(3.0 * (1. - coarseIn.mmCoarse)/list1[j].metric.g00 );
          RectIn.dy = sqrt(3.0 * (1. - coarseIn.mmCoarse)/list1[j].metric.g11 );
          /*RectIn.theta = list1[j].metric.theta + LAL_PI/6;*/
          LALHexagonVertices(&status, &HexaOut, &RectIn);
          fprintf(fpr, "%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n",
              HexaOut.x1, HexaOut.y1,
              HexaOut.x2, HexaOut.y2,
              HexaOut.x3, HexaOut.y3,
              HexaOut.x4, HexaOut.y4,
              HexaOut.x5, HexaOut.y5,
              HexaOut.x6, HexaOut.y6,
              HexaOut.x7, HexaOut.y7);
          fprintf(fpr, "&\n");

        }
      }

      /*plots the ellipses*/
      {
        int Nth =100;
        double th;
        double x,y,phi,a,b, theta;
        for (theta=0; theta<2*3.14; theta+=2*3.14/(double)Nth)
        {
          a = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/list1[j].metric.g00 );
          b = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/list1[j].metric.g11 );
          x = a * cos(theta)/scaling;
          y = b * sin(theta)/scaling;
          phi=list1[j].metric.theta ;
          /*              phi *=1.1;*/
          if (userParams.calque == BCV)
          {
            th = x*cos(phi)-y*sin(phi)+list1[j].params.psi0;
            y  = x*sin(phi)+y*cos(phi)+list1[j].params.psi3;
          }
          else
          {
            th = x*cos(phi)-y*sin(phi)+list1[j].params.t0;
            y  = x*sin(phi)+y*cos(phi)+list1[j].params.t3;
          }
          x  = th;
          fprintf(fpr, "%f %f\n", x, y);
        }
        theta = 0;
        a = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/list1[j].metric.g00 );
        b = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/list1[j].metric.g11 );
        x = a * cos(theta) /scaling;
        y = b * sin(theta)/scaling;
        phi=list1[j].metric.theta ;

        if (userParams.calque == BCV)
        {
          th = x*cos(phi)-y*sin(phi)+list1[j].params.psi0;
          y  = x*sin(phi)+y*cos(phi)+list1[j].params.psi3;
        }
        else
        {
          th = x*cos(phi)-y*sin(phi)+list1[j].params.t0;
          y  = x*sin(phi)+y*cos(phi)+list1[j].params.t3;
        }
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






/*  we should move this function into lal */

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

   p.order = LAL_PNORDER_TWO;

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
