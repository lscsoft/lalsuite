
/*<lalVerbatim file="AstroOmegaTestCV">
Author: Regimbau Tania
$Id$
</lalVerbatim> */

/* 
<lalLaTeX>
\subsection{Program \texttt{AstroOmegaTest.c}}
This programs verifies that the routine LALAstroOmega() gives the expected results (computed separetly with mathematica) for a set of input parameters.

\subsubsection*{Exit codes}
returns 0 on success, otherwise returns 1.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALAstroOmega() 
\end{verbatim}

\vfill{\footnotesize\input{AstroOmegaTestCV}}
</lalLaTeX> 
*/

#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include <lal/AstroOmega.h>

NRCSID (ASTROOMEGATESTC, "$Id$");

static void SDensity (REAL8 *dEgw, REAL8 nu)
 {  
  *dEgw=pow(nu,3.);
  return;
 }

int lalDebugLevel = 0;
int main ( void )
 {
  static LALStatus status;
  AstroOmegaParams p;
  AstroOmegaSourceParams sourcep;
  AstroOmegaCosmoParams cosmop;
  REAL8 nu, test;

  cosmop.ho = 0.68;
  cosmop.density_matter = 0.3;
  cosmop.density_vacuum = 0.7;
  cosmop.density_k = 0.;
  nu = 1100.;
  /*source*/
  sourcep.fact = 2.88E-22;
  sourcep.numax = 4000.;
  sourcep.SDensitySource = SDensity;
  p.cosmoparams = cosmop;
  p.sourceparams = sourcep;

  
   p.extraparams = &nu; 
   LALAstroOmega (&status, &test, nu,&p);
   if (fabs(test-6.74E-11) > 1.E-12)
    {
     printf("error! the right value is 6.74E-11 not %.2e\n",test);
    }
   else printf("omega(%f)= %.2e o.k\n", nu,test);


   
   /* 
  for (nu = 0.;nu < sourcep.numax;nu = nu + 10.)
   {
    p.extraparams = &nu;
    LALAstroOmega (&status, &test, nu, &p);
    printf("%f\t %.2e\n", nu,test);}
   */
  return 0;
  }
