
/* <lalVerbatim file="AstroOmegaGeneralCV">
Author: Regimbau T.
$Id$  
</lalVerbatim> */

/*<lalLaTeX>

\subsection*{Module \texttt{AstroOmegaGeneral.c}}

[compute the energy density spectrum of stochastic backgrounds produced
by cosmological population of astrophysical sources]

\subsubsection*{Prototypes}
%\input{AstroOmegaGeneralCP}
\index{\texttt{LALAstroOmegaSource()}}
\subsubsection*{Description}

The function of this module compute the parameter $\Omega_{gw}(\nu_{o})$
for given cosmological and source models in a general case where the user
define himself the single spectral energy density of his source model.


\subsubsection*{Operating Instructions}
the following program shows how to use the function LALAstroOmegaSource

\begin{verbatim}


#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include <lal/AstroOmega.h>

NRCSID (ASTROOMEGATESTC, "$Id$");

//define here the spectral energy density of a single source 
static void SDensity (REAL8 *dEgw, REAL8 nu)
 {  *dEgw=pow(nu,3.);
  return;
  }

int lalDebugLevel = 0;
int main ()
 {
  static LALStatus status;
  DIntegrateIn  zint;
  AstroOmegaGeneralParams params;
  AstroOmegaGeneralSourceParams gsourcep;
  AstroOmegaCosmoParams cosmop;
  REAL8 omegaz, zmax, nmaxSource, factSource, nu, test;
  //define here the model parameters
  //cosmological model parameters
  cosmop.ho=0.68;
  cosmop.density_matter=0.3;
  cosmop.density_vacuum=0.7;
  cosmop.density_k=0.;
  //source model parameters
  gsourcep.fact=2.88E-22;
  gsourcep.numax=4000.;
  gsourcep.SDensitySource=SDensity;
  params.cosmoparams=cosmop;
  params.gsourceparams=gsourcep;
  for (nu=0.;nu,numax;nu=nu+10.)
   {params.extraparams=&nu;
    LALAstroOmegaSource (&status, &test, nu, &params); 
    printf("omega(%f)= %.2e\n", nu,test);}
   return 0;
  }


\end{verbatim}

\subsubsection*{Notes}
The cosmic star formation rate used in these functions corresponds to Madau model. To obtain the results for Hopkins model just multiplicate the parameter fact by 6.67.
\vfill{\footnotesize\input{AstroOmegaGeneralCV}}

</lalLaTeX> */ 
#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include "AstroOmega.h"

NRCSID (ASTROOMEGAGENERALC, "$Id$");


static void SFR (REAL8 *result, REAL8 z);
static void Ez (REAL8 *result, REAL8 z, AstroOmegaCosmoParams cosmop);
static void dAstroOmegaSource (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p);

void LALAstroOmegaSource (LALStatus *s, REAL8 *omeganu, REAL8 nu, void *p)
 {   
   DIntegrateIn  zint;
   AstroOmegaGeneralParams params;
   AstroOmegaGeneralSourceParams gsourcep;
   AstroOmegaCosmoParams cosmop;
   REAL8 omegaz, zmax, numaxSource, factSource;   
   INITSTATUS (s, "LALAstroOmegaSource", ASTROOMEGAGENERALC);
   ATTATCHSTATUSPTR (s);
  
   params = *((AstroOmegaGeneralParams *)p);
   gsourcep=params.gsourceparams;
   cosmop=params.cosmoparams;
   numaxSource=gsourcep.numax;
   factSource=gsourcep.fact;

   if((nu>=numaxSource)||(nu<=0.)){*omeganu=0.;}
   else
    {
     if (nu<(numaxSource/6.)) {zmax=5.;}
     else {zmax=(numaxSource/nu)-1.;}
 
     zint.function = dAstroOmegaSource;
     zint.xmin     = 0;
     zint.xmax     = zmax;
     zint.type     = ClosedInterval;

     LALDRombergIntegrate (s->statusPtr, &omegaz, &zint, &params); 
     *omeganu = factSource*nu*omegaz;
     }
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }


/*event rate*/
/*source formation rate: lambda*SFR(z), lambda is the mass fraction of NS progenitor*/
static void SFR (REAL8 *result, REAL8 z)
 {  
  /*cosmic star formation rate*/ 
  /*Madau model*/
  *result = 0.23*exp(3.4*z)/(44.7+exp(3.8*z));
  /*Hopkins model*/
  /**result = 1.207*exp(3.836*z)/(39.97+exp(4.163*z));*/
  return;
 }


static void Ez (REAL8 *result, REAL8 z, AstroOmegaCosmoParams cosmop)
 {  
  *result=sqrt(cosmop.density_matter*(1.+z)*(1.+z)*(1.+z)+cosmop.density_k*(1.+z)*(1.+z)+cosmop.density_vacuum);
  return;
 }


static void dAstroOmegaSource (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p)
 {
  
  AstroOmegaGeneralParams params;
  AstroOmegaGeneralSourceParams gsourcep;
  AstroOmegaCosmoParams cosmop; 
  REAL8LALSDensity *SDensitySource;
  REAL8 vRz, vEz, dEgw, nu, nuz, numaxSource, factSource;
  
  INITSTATUS (s, "dAstroOmegaSource", ASTROOMEGAGENERALC);
  ATTATCHSTATUSPTR (s); 
  params = *((AstroOmegaGeneralParams *)p);
  cosmop=params.cosmoparams;
  gsourcep=params.gsourceparams;
  SDensitySource=gsourcep.SDensitySource;
  numaxSource=gsourcep.numax;
  factSource=gsourcep.fact;
  nu = *((REAL8 *)params.extraparams);
  
  /*frequency in the source frame*/
  nuz=(1.+z)*nu;
  /*single spectral energy density in the source frame*/
  SDensitySource(&dEgw, nuz);
  /*event formation rate*/
  SFR(&vRz, z); 
  Ez(&vEz, z, cosmop);
  *domegaz = dEgw*vRz/((1.+z)*vEz);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }
 
static void SDensity (REAL8 *dEgw, REAL8 nu)
 {  
  *dEgw=pow(nu,3.);
  return;
  }


