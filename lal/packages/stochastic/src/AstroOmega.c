
/* <lalVerbatim file="AstroOmegaCV">
Author: Regimbau T.
$Id$  
</lalVerbatim> */

/*<lalLaTeX>

\subsection*{Module \texttt{AstroOmega.c}}

[compute the energy density spectrum of stochastic backgrounds produced
by cosmological population of astrophysical sources]

\subsubsection*{Prototypes}
%\input{AstroOmegaCP}
\index{\texttt{LALAstroOmega()}}
\subsubsection*{Description}
The function of this module computes the energy density parameter $\Omega_{gw}(\nu_{o})$ for a given source and a given cosmological model.
The spectral properties of the stochastic background are characterized by the dimensionless parameter
\begin{equation}
\Omega _{gw}=\frac{1}{\rho _{c}}\frac{d\rho _{gw}}{d\log \nu _{o}}
\end{equation}
where $\rho _{gw}$ is the gravitational energy density, $\nu _{o}$ the frequency in the
observer frame and
\begin{equation}
\rho _{c}=\frac{3H_{o}^{2}}{8\pi G}
\end{equation}
the critical energy density to close the Universe today.

For a given cosmological model, a stochastic background of
astrophysical origin is fully determined by the source formation rate as a function of redshift,
$dR(z)$ and the gravitational spectral density of a single source, ${dE_{gw}}\over{d\nu}$.
Following Ferrari et al. (1999):
\begin{equation}
\Omega _{gw}(\nu _{o})=\frac{1}{c^{3} \rho _{c}}{\nu _{o}}F_{\nu _{o}}
\end{equation}
where
\begin{equation}
F_{\nu _{o}}=\int_{0}^{z_{\max }}f_{\nu _{o}}dR(z)
\end{equation}
is the gravitational wave flux at the
frequency $\nu _{o}$ integrated over all the sources.
The gravitational flux of a single source located between z, z+dz is:
\begin{equation}
f_{\nu _{o}}=\frac{1}{4\pi d_{L}^{2}}\frac{dE_{gw}}{d\nu }(1+z)
\end{equation}
where $d_{L}=(1+z)r$ is the distance luminosity and  $\nu=(1+z)\nu _{o}$
the frequency in the source frame.
The event rate between z, z+dz as observed in our frame is given by:
\begin{equation}
dR(z)=\lambda _{p}\frac{R_{c}(z)}{1+z}{\frac{{dV}}{{dz}}}dz
\end{equation}
where $R_{c}(z)$ is the cosmic star formation rate and  $\lambda_{p}$
is the mass fraction of the source progenitors. The
term (1 +z) in the denominator accounts for the time dilatation by cosmic expansion of the
observed rate.
The element of the comoving volume is
\begin{equation}
dV = 4\pi r^2{{c}\over{H_o}}{{dz}\over{E(\Omega_i,z)}}
\end{equation}
where the function E$(\Omega_i,z)$ is
defined by the equation:
\begin{equation}
E(\Omega_i,z) =[(1+z)^2(1+z\Omega_m)-z(2+z)\Omega_v]^{1/2}
\end{equation}
where $\Omega_m$ and $\Omega_v$ are respectively the density parameters
due to matter (baryonic and non-baryonic) and the vacuum, corresponding
to a non-zero cosmological constant. The equivalent
density parameter due to the spatial curvature satisfies $\Omega_k =
1-\Omega_m -\Omega_v$. 
The cosmic star formation rate is computed as:
\begin{equation}
R_{c}(z) = 
R_{SF2}(z)
\frac{E(\Omega_i,z)}{(1+z)^{3/2}}\; 
\mathrm{M}_{\odot}\,\mathrm{yr}^{-1}\,\mathrm{Mpc}^{-3}
\end{equation}
where 
\begin{equation}
R_{SF2}(z) = \frac{0.15\,e^{3.4z}}{(22 + e^{3.4z})} \;
\mathrm{M}_{\odot}\,\mathrm{yr}^{-1}\,\mathrm{Mpc}^{-3}
\end{equation}
is the cosmic star formation rate (Madau \&Porciani, 2001) in a matter-dominateduniverse ($\Omega_{m}=1$).
We consider an Einstein - de Sitter universe,
with non-zero cosmological constant and with $\Omega_m$ = 0.30
 and $\Omega_v$ = 0.70.
These values are suggested by
recent BOOMERANG and MAXIMA results (de Bernardis et al. 2000; Hanany et al. 2000) on
the power spectra of the cosmic microwave background and by observations of distant type 
Ia supernovas (Perlmutter et al. 1999; Schmidt et al. 1998).
The Hubble parameter H$_o$ is taken to be $73\,\mathrm{km\,s}^{-1}\mathrm{Mpc}^{-1}$ (Freedman et al., 2001).
Combining the previous equations one obtains:
\begin{equation}
\Omega _{gw}(\nu _{o})=\frac{8\pi G}{3c^{2}H_{o}^{3}}\lambda_{p}\, \nu_{o}\int_{0}^{z_{\sup }}\frac{dE_{gw}}{d\nu }\frac{R_{c}(z)}{E(z)(1+z)^2}dz
\end{equation} 
The upper limit of the integral is determined by the cutoff frequency in the
source frame, $\nu _{\sup }$, as:
\begin{equation}
z_{\sup }=\frac{\nu _{\sup }}{\nu_{o}} - 1
\end{equation}
Note that we are limited to $z<5$ since the cosmic star formation rate is not
modeled above. This restriction is without significant consequence
since objects with $z>5 $ contribute very little to the integrated signal.



\subsubsection*{Operating Instructions}
the following program shows how to use the function LALAstroOmega

\begin{verbatim}


#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include "AstroOmega.h"

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
  //cosmological parameters
  cosmop.ho=0.73;
  cosmop.density_matter=0.3;
  cosmop.density_vacuum=0.7;
  cosmop.density_k=0.;

  //source parameters for rotating pulsars
  gsourcep.fact=2.88E-22;
  gsourcep.numax=4000.;
  gsourcep.SDensitySource=SDensity;
  params.cosmoparams=cosmop;
  params.gsourceparams=gsourcep;
  for (nu=0.; nu < numax ; nu=nu+10.)
   {params.extraparams=&nu;

    LALAstroOmegaSource (&status, &test, nu, &params); 
    printf("omega(%f)= %.2e\n", nu,test);}
   return 0;
  }


\end{verbatim}

\subsubsection*{Notes}
\vfill{\footnotesize\input{AstroOmegaGeneralCV}}

</lalLaTeX> */ 

#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include "AstroOmega.h"

NRCSID (ASTROOMEGAC, "$Id$");


static void SFR (REAL8 *result, REAL8 z);
static void dAstroOmega (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p);

void LALAstroOmega (LALStatus *s, REAL8 *omeganu, REAL8 nu, void *p)
 {   
   DIntegrateIn  zint;
   AstroOmegaParams params;
   AstroOmegaSourceParams sourcep;
   AstroOmegaCosmoParams cosmop;
   REAL8 omegaz, zmax, numaxSource, factSource;   
   INITSTATUS (s, "LALAstroOmega", ASTROOMEGAC);
   ATTATCHSTATUSPTR (s);
  
   params = *((AstroOmegaParams *)p);
   cosmop=params.cosmoparams;
   sourcep = params.sourceparams;
   numaxSource= sourcep.numax;
   factSource = sourcep.fact;

   if((nu>=numaxSource)||(nu<=0.)){*omeganu=0.;}
   else
    {
     if (nu<(numaxSource/6.)) {zmax=5.;}
     else {zmax=(numaxSource/nu)-1.;}
 
     zint.function = dAstroOmega;
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



/*cosmic star formation rate */
static void SFR (REAL8 *result, REAL8 z)
 {  
  /*cosmic star formation rate*/ 
  /*Madau & Pozetti, 2001 */
  *result = 0.15 * exp(3.4*z) / (22. + exp(3.4*z));
  
  return;
 }



static void dAstroOmega (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p)
 {
  
  AstroOmegaParams params;
  AstroOmegaSourceParams sourcep; 
  AstroOmegaCosmoParams cosmop;
  REAL8LALSDensity *SDensitySource;
  REAL8 vRz, dEgw, nu, nuz, numaxSource, factSource;
  
  INITSTATUS (s, "dAstroOmega", ASTROOMEGAC);
  ATTATCHSTATUSPTR (s); 
  params = *((AstroOmegaParams *)p);
  cosmop = params.cosmoparams;
  sourcep = params.sourceparams;
  SDensitySource = sourcep.SDensitySource;
  numaxSource = sourcep.numax;
  factSource = sourcep.fact;
  nu = *((REAL8 *)params.extraparams);
  
  /*frequency in the source frame*/
  nuz = (1.+z) * nu;
  /*single spectral energy density in the source frame*/
  SDensitySource(&dEgw, nuz);
  /*cosmic formation rate*/
  SFR(&vRz, z); 
  *domegaz = dEgw * vRz / pow((1.+z),3.5);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }

