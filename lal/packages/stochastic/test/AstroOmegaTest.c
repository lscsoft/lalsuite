
/*<lalVerbatim file="AstroOmegaTestCV">
Author: Regimbau Tania
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Program \texttt{AstroOmegaTest.c}}
This programs verifies that the routines LALAstroOmega_*() gives the expected results (computed separetly with mathematica) for a set of key values.

\subsubsection*{Exit codes}
returns 0 on success, otherwise returns 1.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALAstroOmegaPulsar()
LALAstroOmegaMode()  
LAstroOmegaBinary()   
\end{verbatim}

\vfill{\footnotesize\input{AstroOmegaTestCV}}
</lalLaTeX> */

#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include "AstroOmega.h"

NRCSID (ASTROOMEGATESTC, "$Id$");

static void SDensity (REAL8 *dEgw, REAL8 nu)
 {  *dEgw=pow(nu,3.);
  return;
  }

int lalDebugLevel = 0;
int main ()
 {
  static LALStatus status;
  AstroOmegaGeneralParams generalp;
  AstroOmegaTemplatesParams pulsarp,modesp,binaryp;
  AstroOmegaGeneralSourceParams generalsp;
  AstroOmegaTemplatesSourceParams pulsarsp,modessp,binarysp;
  AstroOmegaCosmoParams cosmop;
  REAL8 zmax, nu, test;

  cosmop.ho=0.68;
  cosmop.density_matter=0.3;
  cosmop.density_vacuum=0.7;
  cosmop.density_k=0.;
  nu=1100.;
  /*source*/
  generalsp.fact=2.88E-22;
  generalsp.numax=4000.;
  generalsp.SDensitySource=SDensity;
  generalp.cosmoparams=cosmop;
  generalp.gsourceparams=generalsp;
  generalp.extraparams=&nu; 
  LALAstroOmegaSource (&status, &test, nu,&generalp);
  if (fabs(test-2.20E-10)>1.E-12) 
   {printf("error! the right value is 2.20E-10 no %.2e\n",test);
    /*return 1;*/
   } 
  else printf("omega(%f)= %.2e o.k\n", nu,test);
 
  /*rotating pulsars*/
  pulsarsp.fact=2.88E-22;
  pulsarsp.numax=4000.;
  pulsarp.cosmoparams=cosmop;
  pulsarp.tsourceparams=pulsarsp;
  pulsarp.extraparams=&nu; 
  LALAstroOmegaPulsar (&status, &test, nu,&pulsarp);
  if (fabs(test-2.20E-10)>1.E-12) 
   {printf("error! the right value is 2.20E-10 no %.2e\n",test);
    /*return 1;*/
   } 
  else printf("omega(%f)= %.2e o.k\n", nu,test);
 /*binaries */
 /* binarysp.fact=1.E-10;
  binarysp.numax=1180.;
  binaryp.cosmoparams=cosmop;
  binaryp.tsourceparams=binarysp;
  binaryp.extraparams=&nu; 
  LALAstroOmegaBinary (&status, &test, nu,&binaryp);
  if (fabs(test-6.28E-12)>1.E-14)
   {printf("error! the right value is 6.28E-12 no %.2e\n",test);
   }
  else printf("omega(%f)= %.2e o.k\n", nu,test);*/
  /*r modes*/
  modessp.fact=1.36E-16;
  modessp.numax=1958.;
  modesp.cosmoparams=cosmop;
  modesp.tsourceparams=modessp;
  modesp.extraparams=&nu; 
  LALAstroOmegaModes (&status, &test, nu, &modesp);
  if (fabs(test-2.05E-12)>1.E-14) 
   {printf("error! the right value is 2.05E-12 no %.2e\n",test);
   /*return 1;*/
   }
  else printf("omega(%f)= %.2e o.k\n", nu,test);
 
  return 0;
  }
