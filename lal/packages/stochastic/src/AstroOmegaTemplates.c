/* <lalVerbatim file="AstroOmegaTemplatesCV">
Author: Regimbau T.
$Id$  
</lalVerbatim> */

/*<lalLaTeX>

\subsection*{Module \texttt{AstroOmegaTermplates.c}}

[compute the energy density spectrum of stochastic backgrounds produced
by cosmological population of astrophysical sources]

\subsubsection*{Prototypes}
%\input{AstroOmegaTemplatesCP}
\index{\texttt{LALAstroOmegaPulsar()}}
\index{\texttt{LALAstroOmegaModes()}}
\index{\texttt{LALAstroOmegaBinary()}}

\subsubsection*{Description}
The stochastic background is characterized by the unitless parameter 
$$
\Omega _{gw}=\frac{1}{\rho _{c}}\frac{d\rho _{gw}}{d\log \nu _{o}}
$$
where $\rho _{gw}$ is the gravitational energy density, $\nu _{o}$ the
frequency in the observer frame and $\rho _{c}=\frac{3H_{o}^{2}}{8\pi G}$
the energy density to close the Universe. For a given cosmological model, a
gravitational background of astrophysical origin is then fully determined by
the gravitational energy spectral density of a single source $\frac{dE_{gw}}{d\nu }$ and the source formation rate (or the coalescence rate for binary
systems) as a function of redshift $dR(z)$. We have: 
$$
\Omega _{gw}(\nu _{o})=\frac{1}{\rho _{c}c^{3}}{\nu _{o}}F_{\nu _{o}}
$$
where 
$$
F_{\nu _{o}}=\int_{0}^{z_{\sup }}f_{\nu _{o}}dR(z)
$$
is the total flux spectral density. The contribution of a single source
located between z and z+dz is given by: 
$$
f_{\nu _{o}}=\frac{1}{4\pi d_{L}^{2}}\frac{dE_{gw}}{d\nu _{o}}
$$
where $d_{L}=(1+z)r$ is the distance luminosity. The observed frequency
being redshifted with respect to the source frequency by, 
$$
\frac{dE_{gw}}{d\nu _{o}}=\frac{dE_{gw}}{d\nu }(1+z)
$$
{\large {Warning! }}The upper limit of the integral upper depends on the
frequency cutoff $\nu _{\sup }$ in the star frame. 
$$\left\lbrace 
            \begin{array}{l}
              5$ if $\nu_{o}<\frac{\nu_{sup}}{6} \\ 
              \\ 
              \frac{\nu_{sup}}{\nu_{o}}-1$ if  $\nu_{o}\geq \frac{\nu_{sup}}{6}
            \end{array}
          \right.$$
The functions of this module compute the parameter $\Omega_{gw}(\nu_{o})$
for given cosmological and source model. LALAstroOmegaPulsar(),
LALAstroOmegaMode() and LALAstroOmegaBinary() correspond respectively to
contributions from rotating pulsars, bar and r mode from fast newborn
neutron stars and binary system coalescence. The function
LALAstroOmegaSource()corresponds to a more general case where the user
define himself the single spectral energy density of his source model.
\paragraph*{LALAstroOmegaPulsar()}
\subparagraph*{single spectral energy density}
$$
\frac{dE_{gw}}{d\nu }=\frac{256G\pi^{6}}{5c^{5}}\varepsilon^{2}I^{2}\frac{\tau_{0}}{P_{0}^{2}}\nu^{3}=K\nu^{3} 
$$
where $\nu \leq \nu_{K}$ (the frequency cutoff corresponds to the Keplerian
rotational angular velocity $\Omega_{K}$)
$P_{0}$is the initial rotational period and $\tau_{0}$ the
magnetic breaking timescale.
$K$ can be written in term of $<\frac{1}{B_{0}^{2}}>$ (where $B_{0}$ is the
initial magnetic field) using the relation $B_{0}^{2}=\frac{3Ic^{3}}{4\pi
^{2}R^{6}}\frac{P_{0}^{2}}{\tau_{0}}$.

\subparagraph*{source formation rate}
The neutron star formation rate between z and z+dz is: 
$$
dR(z)=\lambda _{NS}R_{c}(z){\frac{{dV}}{{dz}}}dz
$$
where $R_{c}(z)$ is the cosmic star formation rate. An estimate of
this parameter is given by Madau and Pozzetti (1999) 
$$
R_{c}(z)=\frac{0.23e^{3.4z}}{44.7+e^{3.8z}} \mathrm{M}_{\odot}\mathrm{an}^{-1}\mathrm{Mpc}^{-3}
$$
and by Hopkins et al. (2001) 
$$
R_{c}(z)=\frac{1.207e^{3.84z}}{39.97+e^{4.16z}} \mathrm{M}_{\odot}\mathrm{an}^{-1}\mathrm{Mpc}^{-3}
$$
The factor $\lambda_{NS}$ is the mass fraction of neutron star progenitors. 
$$
\lambda_{NS}=\int_{10M_{\odot }}^{40M_{\odot }}\xi (m)dm
$$
where the initial mass function is a Salpeter law, 
$$
\xi (m)\propto m^{-(1+x)}
$$
with $x=1.35$ and normalizing between $0.1-80$ M$_{\odot }$ through $\int_{0.1M_{\odot }}^{80M_{\odot }}m\xi (m)dm$ we obtain: $\lambda_{p}=4.84\times 10^{-3}$ M$_{\odot }^{-1}$

The comoving volume element is 
$$
dV=4\pi r^{2}{\frac{{c}}{{H_{o}}}}{\frac{{dz}}{E{(\Omega _{i},z)}}}
$$
where $E(\Omega _{i},z)$ is defined by 
$$
E(z)=[\Omega _{m}(1+z)^{3}+\Omega _{k}(1+z)^{2}+\Omega _{v}]^{1/2}
$$
where 
$$
\Omega _{m}+\Omega _{v}+\Omega _{k}=1
$$
$\Omega _{m},\Omega _{v}$ and $\Omega _{k}$ are respectively the density
parameters due to the matter (both baryonic and non-baryonic), to the vacuum
and to the space curvature, for a non zero cosmological constant.
\subparagraph{stochastic background}
In order to save computer time, we groupe the various factors which appear
in the expression of $\Omega_{gw}$.
Defining 
$$
\frac{dE_{gw}}{d\nu }=K\frac{d\bar{E}_{gw}}{d\nu }
$$
we write 
$$
\Omega _{gw}(\nu _{o})=fact\nu _{o}\int_{0}^{z_{\max }}\frac{d\bar{E}_{gw}}{d\nu }\frac{R_{c}(z)}{E(z)(1+z)}dz
$$
where $fact=\frac{8\pi G}{3c^{2}H_{0}^{3}}K\lambda_{NS}$
\paragraph{LALAstroOmegaMode()}
\subparagraph{single spectral energy density}
$$
\frac{dE_{gw}}{d\nu }=\frac{2E_{o}}{\nu _{sup }^{2}}\nu
$$
For r modes, $E_{o}$ corresponds to the rotational energy dissipated during the
instability phase and $\nu_{sup}=\frac{4}{3}\nu_{K}$
For bar modes,$E_{o}$corresponds to the difference between the
equilibrium energies of Maclaurin spheroid and Dedekind ellipsoid.
\subparagraph{formation rate}
$$
dR(z)=\xi \lambda _{NS}R_{c}(z){\frac{{dV}}{{dz}}}dz
$$
where $\xi $ is the fraction of objects fast and hot enough to cross
the instability window.
\subparagraph{stochastic background}
$$
\Omega _{gw}(\nu _{o})=fact\nu _{o}\int_{0}^{z_{\max }}\frac{d\bar{E}_{gw}}{d\nu }\frac{R_{c}(z)}{E(z)(1+z)}dz
$$
where $fact=\frac{8\pi G}{3c^{2}H_{0}^{3}}K \xi \lambda_{NS}$
\paragraph{LALAstroOmegaBinary()}
\subparagraph{single spectral energy density}
$$
\frac{dE_{gw}}{d\nu }=\frac{(G\pi )^{2/3}}{3}\frac{m_{1}m_{2}}{M^{1/3}}\nu ^{-1/3}=K\nu ^{-1/3}
$$
where $\nu \leq \nu_{sup}$. The frequency cutoff corresponds to the orbital
frequency when the two stars start to merge. For NS/NS systems, $m_{1}=m_{2}=1.4$ $M_{\odot }$, $K=5.2\times 10^{50}$ cgs and $\nu _{\sup}=1180$ Hz.
\subparagraph{coalescence rate}
{\Large Warning!}:In this case, $dR(z)$ is the coalescence rate and not the
formation rate.

If we call $R_{b}$ the NS/NS formation rate and $P_{z}(z)$ the probability
for such systems to collapse within z , the cosmic coalescence rate is given
by: 
$$
dR(z)=\frac{{dV}}{{dz}}dz\int_{z}^{z_{\sup }}R_{b}(z^{\prime
}-z)P_{z}(z^{\prime })dz^{\prime }
$$
The binary system formation rate at z is: 
$$
R_{b}(z)=\lambda _{b}R_{c}(z)
$$
where $R_{c}(z)$ is the cosmic star formation rate and $\lambda_{b}$ the
mass fraction of binary system progenitors.
Under the assumption that $\lambda _{b}$ is the same in all galaxies and
that it doesn't depend on their evolution stage: 
$$
\lambda _{b}=\frac{R_{b}}{R}
$$
where $R$ is the present galactic star formation rate and $R_{b}$ the actual
binary system birth rate in our galaxy. $R$ is computed under the assumption
that the star formation rate at a given time is proportional to the
available mass of gas $M_{g}(t)$: 
$$
R(t)=kM_{g}(t)
$$
Under the so called ''Instantaneous Recycling Approximation'', the ejected
mass is: 
$$
\frac{dM_{e}(t)}{dt}=xkM_{g}(t)
$$
and the variation rate of the available gas mass is: 
$$
\frac{dM_{g}(t)}{dt}=-R(t)+\frac{dM_{e}(t)}{dt}=-(1-x)kM_{g}(t)
$$
Integrating this expression we obtain: 
$$
M_{g}(t)=M_{0}e^{-(1-x)kt}
$$
where $M_{0}$ is the initial gas mass.
Writting the previous equation today ($t=T$) , 
$$
M_{HI}=M_{0}e^{-(1-x)kT} 
$$
where $T$ is the galaxy age and $M_{HI}$ the present gas mass, we deduce $k$: 
$$
k=\frac{1}{T(1-x)}\ln(\frac{M_{o}}{M_{HI}})
$$
At $T$ we have 
$$
R=kM_{HI}
$$
Considering $T=15\mathrm{Gyr}$, $x=0.1$, $\frac{M_{HI}}{M_{0}}=8$ and
$M_{HI}=6.5\mathrm{GM}_{\odot}$, we obtain
$R=1.2\mathrm{M}_{\odot}\mathrm{yr}^{-1}$.

Binary systems are born with a set of orbital parameters which determine
their evolution and in particular the time before their collapse. The
probability $P_{z}(z)$ can be deduced from the probability $P_{\tau }(\tau )$
to form a system with a coalescence time $\tau $: 
$$
P_{z}(z)=P_{\tau }(\tau )\frac{d\tau }{dz}
$$
According to de Freitas Pacheco (1997) [27] we have: 
$$
P_{\tau }(\tau )=\frac{B}{\tau }
$$
Assuming that $\tau $ is in the range $0.1-20$ Gyr, the normalization
constante is $B=\frac{1}{\ln 200}=0.19$. Using instead of $\tau $ the
variable $t=\tau -t_{o}$ where $t_{o}=0.1$Gyr is the minimal coalescence
time: 
$$
P_{t}(t)=\frac{B}{t_{o}+t}
$$
and: 
$$
P_{z}(z)=P_{t}(t)\frac{dt}{dz}
$$
With 
$$
\frac{dt^{\prime }}{dz^{\prime }}=\frac{1}{H_{o}}\frac{1}{(1+z^{\prime
})E(z^{\prime })}
$$
we obtain 
$$
P_{z}(z^{\prime})=\frac{B}{H_{o}(t_{o}+t^{\prime }(z^{\prime }))}\frac{1}{(1+z^{\prime })E(z^{\prime })}
$$
where 
$$
H_{o}t^{\prime }(z^{\prime })=\int_{z}^{z^{\prime }}\frac{dz^{\prime \prime }}{(1+z^{\prime \prime })E(z^{\prime \prime })}
$$
Introducing 
$$
\bar{P}(z)=\frac{P(z)}{B}
$$
and 
$$
F(z)=\int_{z}^{z_{\sup }}R_{c}(z^{\prime }-z)\bar{P}(z^{\prime })dz^{\prime }
$$
the final expression of the coalescence rate is given by: 
$$
dR(z)=\lambda _{b}BF(z)\frac{{dV}}{{dz}}dz
$$
\subparagraph{stochastic background}
$$
\Omega _{gw}(\nu _{o})=fact\nu _{o}\int_{0}^{z_{\max }}\frac{d\bar{E}%
_{gw}(\nu _{o})}{d\nu }\frac{F(z)}{E(z)(1+z)}dz
$$
where $fact=\frac{8\pi G}{3c^{2}H_{0}^{3}}K\lambda _{b}B$

\subsubsection*{Operating Instructions}

the following program shows how to use the predefined functions: here LALAstroOmegaPulsar
\begin {verbatim}

#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include "AstroOmega.h"

NRCSID (ASTROOMEGATESTC, "$Id$");

int lalDebugLevel = 0;
int main ()
 {
  static LALStatus status;
  DIntegrateIn  zint;
  AstroOmegaTemplatesParams params;
  AstroOmegaTemplatesSourceParams pulsarp;
  AstroOmegaCosmoParams cosmop;
  REAL8 omegaz, zmax, nu, omega;
  //define here the model parameters
  //cosmological model parameters
  cosmop.ho=0.68;
  cosmop.density_matter=0.3;
  cosmop.density_vacuum=0.7;
  cosmop.density_k=0.;
  //source model parameters
  pulsarp.fact=2.88E-22;
  pulsarp.numax=4000.;
  params.cosmoparams=cosmop;
  params.gsourceparams=gsourcep;
  for (nu=0.;nu,gsourcep.numax;nu=nu+10.)
   {params.extraparams=&nu;
    LALAstroOmegaPulsar (&status, &omega, nu, &params); 
    printf("omega(%f)= %.2e\n", nu,omega);}
   return 0;
  }
\end{verbatim}
\subsubsection*{Uses}
\begin{verbatim}
These routines use the function \index{\texttt{LALDRombergIntegrate()}}
\end{verbatim}

\subsubsection*{Notes}
The cosmic star formation rate used in these functions corresponds to Madau model. To obtain the results for Hopkins model just multiplicate the parameter fact by 6.67.
\vfill{\footnotesize\input{AstroOmegaTemplatesCV}}

</lalLaTeX> */ 

#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include <lal/AstroOmega.h>
NRCSID (ASTROOMEGATEMPLATESC, "$Id$");


static void SFR (REAL8 *result, REAL8 z);
static void Ez (REAL8 *result, REAL8 z, AstroOmegaCosmoParams cosmop);
static void Hodt (LALStatus *s, REAL8 *result, REAL8 z2, void *p);
static void dSCR (LALStatus *s, REAL8 *result, REAL8 z1, void *p);
static void SCR (LALStatus *s, REAL8 *result, REAL8 z, void *p);
static void SDensityPulsar (REAL8 *dEgw, REAL8 nu);
static void SDensityModes (REAL8 *dEgw, REAL8 nu);
static void SDensityBinary (REAL8 *dEgw, REAL8 nu);
static void dAstroOmegaPulsar (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p);
static void dAstroOmegaModes (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p);
static void dAstroOmegaBinary(LALStatus *s, REAL8 *domegaz, REAL8 z, void *p);




void LALAstroOmegaPulsar (LALStatus *s, REAL8 *omeganu, REAL8 nu, void *p)
 {   
  DIntegrateIn  zint;
  AstroOmegaTemplatesParams params;
  AstroOmegaTemplatesSourceParams tsourcep;
  AstroOmegaCosmoParams cosmop;
  REAL8 omegaz, zmax, numaxPulsar,factPulsar; 
  INITSTATUS (s, "LALAstroOmegaPulsar", ASTROOMEGATEMPLATESC);
  ATTATCHSTATUSPTR (s);
  params = *((AstroOmegaTemplatesParams *)p);
  tsourcep=params.tsourceparams;
  cosmop=params.cosmoparams;
  numaxPulsar=tsourcep.numax;
  factPulsar=tsourcep.fact;
  if((nu>=numaxPulsar)||(nu<=0.)){*omeganu=0.;}
   else
     {
     if (nu<(numaxPulsar/6.)) {zmax=5.;}
     else {zmax=(numaxPulsar/nu)-1.;}
 
     zint.function = dAstroOmegaPulsar;
     zint.xmin     = 0;
     zint.xmax     = zmax;
     zint.type     = ClosedInterval;

     LALDRombergIntegrate (s->statusPtr, &omegaz, &zint, &params); 
     *omeganu = factPulsar*nu*omegaz;
     }
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }

void LALAstroOmegaModes (LALStatus *s, REAL8 *omeganu, REAL8 nu, void *p)
 {
  DIntegrateIn  zint;
  AstroOmegaTemplatesParams params;
  AstroOmegaTemplatesSourceParams tsourcep;
  AstroOmegaCosmoParams cosmop;
  REAL8 omegaz, zmax, numaxModes,factModes; 
  INITSTATUS (s, "LALAstroOmegaModes", ASTROOMEGATEMPLATESC);
  ATTATCHSTATUSPTR (s);
  params = *((AstroOmegaTemplatesParams *)p);
  tsourcep=params.tsourceparams;
  cosmop=params.cosmoparams;
  numaxModes=tsourcep.numax;
  factModes=tsourcep.fact;   
 
  if((nu>=numaxModes)||(nu<=0.)){*omeganu=0.;}
   else
     {
     if (nu<(numaxModes/6.)) {zmax=5.;}
     else {zmax=(numaxModes/nu)-1.;}
 
     zint.function = dAstroOmegaModes;
     zint.xmin     = 0;
     zint.xmax     = zmax;
     zint.type     = ClosedInterval;

     LALDRombergIntegrate (s->statusPtr, &omegaz, &zint, &params); 
     *omeganu = factModes*nu*omegaz;
     }
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }

void LALAstroOmegaBinary (LALStatus *s, REAL8 *omeganu, REAL8 nu,void *p)
 {
  DIntegrateIn  zint;
  AstroOmegaTemplatesParams params;
  AstroOmegaTemplatesSourceParams tsourcep;
  AstroOmegaCosmoParams cosmop;
  REAL8 omegaz, zmax, numaxBinary,factBinary; 
  INITSTATUS (s, "LALAstroOmegaBinary", ASTROOMEGATEMPLATESC);
  ATTATCHSTATUSPTR (s);
  params = *((AstroOmegaTemplatesParams *)p);
  tsourcep=params.tsourceparams;
  cosmop=params.cosmoparams;
  numaxBinary=tsourcep.numax;
  factBinary=tsourcep.fact;   
  if((nu>=numaxBinary)||(nu<=0.)){*omeganu=0.;}
   else
     {
      if (nu<(numaxBinary/6.)) {zmax=5.;}
       else {zmax=(numaxBinary/nu)-1.;}

       zint.function = dAstroOmegaBinary;
       zint.xmin     = 0;
       zint.xmax     = zmax;
       zint.type     = ClosedInterval;

       LALDRombergIntegrate (s->statusPtr, &omegaz, &zint, &params);
       *omeganu = factBinary*nu*omegaz;
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
static void Hodt (LALStatus *s, REAL8 *result, REAL8 z2, void *p)
 {
  AstroOmegaCosmoParams cosmop;
  REAL8 vEz;
  INITSTATUS (s, "Hodt", ASTROOMEGATEMPLATESC);
  ATTATCHSTATUSPTR (s);
  cosmop = *((AstroOmegaCosmoParams *)p); 
  Ez(&vEz,z2,cosmop);
  *result=1./((1.+z2)*vEz);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }

static void dSCR (LALStatus *s, REAL8 *result, REAL8 z1, void *p)
 {
  DIntegrateIn  z2int; 
  AstroOmegaTemplatesParams params;
  AstroOmegaTemplatesSourceParams tsourcep;
  AstroOmegaCosmoParams cosmop;  
  REAL8 vHoto, vHot, vEz, vRc,vP,z,deltaz; 
  INITSTATUS (s, "dSCR", ASTROOMEGATEMPLATESC);
  ATTATCHSTATUSPTR (s);
  params = *((AstroOmegaTemplatesParams *)p);
  cosmop=params.cosmoparams;
  z = *((REAL8 *)params.extraparams);
  vHoto=0.007;
  deltaz=z1-z;
  if (deltaz<=0.){vHot=0.;}
   else  
    {  
     z2int.function = Hodt;
     z2int.xmin     = z;
     z2int.xmax     = z1;
     z2int.type     = ClosedInterval;
     LALDRombergIntegrate (s->statusPtr, &vHot, &z2int,&cosmop);
     }
  SFR(&vRc, deltaz); 
  Ez(&vEz,z1,cosmop);
  vP=1./((vHoto+vHot)*(1.+z1)*vEz);  
  *result=vRc*vP;
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }


static void SCR (LALStatus *s, REAL8 *result, REAL8 z, void *p)
 {  

  DIntegrateIn  z1int; 
  AstroOmegaTemplatesParams params, paramsd;
  AstroOmegaTemplatesSourceParams tsourcep;
  REAL8 zmax,nu,deltaz,vSCR;
  INITSTATUS (s, "dSCR", ASTROOMEGATEMPLATESC);
  ATTATCHSTATUSPTR (s);
  params = *((AstroOmegaTemplatesParams *)p);
  tsourcep=params.tsourceparams;
  nu = *((REAL8 *)params.extraparams);
  paramsd.cosmoparams=params.cosmoparams;
  paramsd.tsourceparams=params.tsourceparams;
  paramsd.extraparams=&z;
  if (nu<(tsourcep.numax/6.)) {zmax=5.;}
   else {zmax=(tsourcep.numax/nu)-1.;}
  deltaz=zmax-z; 
  if (deltaz<=0.){vSCR=0.;}
   else  
    {   
     z1int.function = dSCR;
     z1int.xmin     = z;
     z1int.xmax     = zmax;
     z1int.type     = ClosedInterval;
     LALDRombergIntegrate (s->statusPtr, &vSCR, &z1int, &paramsd);
    } 
  *result=vSCR;
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }



/*single spectral density*/
/*single spectral density: K*dEgw*/
static void SDensityPulsar (REAL8 *dEgw, REAL8 nu)
 {
  *dEgw=pow(nu,3.);
  return;
 }
static void SDensityModes (REAL8 *dEgw, REAL8 nu)
 {
  *dEgw=nu;
  return;
 }
static void SDensityBinary (REAL8 *dEgw, REAL8 nu)
 {
  *dEgw=pow(nu,-(1./3.));
  return;
 }


static void dAstroOmegaPulsar (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p)
 {

  AstroOmegaTemplatesParams params;
  AstroOmegaTemplatesSourceParams tsourcep;
  AstroOmegaCosmoParams cosmop; 
  REAL8 vRz, vEz, dEgw,nu,nuz, numaxPulsar, factPulsar;
  INITSTATUS (s, "dAstroOmegaPulsar", ASTROOMEGATEMPLATESC);  
  ATTATCHSTATUSPTR (s); 
  params = *((AstroOmegaTemplatesParams *)p);
  cosmop=params.cosmoparams;
  tsourcep=params.tsourceparams;
  numaxPulsar=tsourcep.numax;
  factPulsar=tsourcep.fact;
  nu = *((REAL8 *)params.extraparams);
  /*frequency in the NS frame*/
  nuz=(1.+z)*nu;
  /*single spectral energy density in the source frame*/
  SDensityPulsar(&dEgw, nuz);
  /*event formation rate*/
  SFR(&vRz, z); 
  Ez(&vEz, z, cosmop);
  *domegaz = dEgw*vRz/((1.+z)*vEz);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }

static void dAstroOmegaModes (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p)
 {
  AstroOmegaTemplatesParams params;
  AstroOmegaTemplatesSourceParams tsourcep;
  AstroOmegaCosmoParams cosmop; 
  REAL8 vRz, vEz, dEgw,nu,nuz, numaxModes, factModes;
  INITSTATUS (s, "dAstroOmegaModes", ASTROOMEGATEMPLATESC);  
  ATTATCHSTATUSPTR (s); 
  params = *((AstroOmegaTemplatesParams *)p);
  cosmop=params.cosmoparams;
  tsourcep=params.tsourceparams;
  numaxModes=tsourcep.numax;
  factModes=tsourcep.fact;
  nu = *((REAL8 *)params.extraparams);
  /*frequency in the NS frame*/
  nuz=(1.+z)*nu;
  /*single spectral energy density in the source frame*/
  SDensityModes(&dEgw, nuz);
  /*event formation rate*/
  SFR(&vRz, z); 
  Ez(&vEz, z, cosmop);
  *domegaz = dEgw*vRz/((1.+z)*vEz);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }

static void dAstroOmegaBinary (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p)
 {
  AstroOmegaTemplatesParams params;
  AstroOmegaTemplatesSourceParams tsourcep;
  AstroOmegaCosmoParams cosmop; 
  REAL8 vRz, vEz, dEgw,nu,nuz, numaxBinary, factBinary;
  INITSTATUS (s, "dAstroOmegaBinary", ASTROOMEGATEMPLATESC);
  ATTATCHSTATUSPTR (s);
  params = *((AstroOmegaTemplatesParams *)p);
  cosmop=params.cosmoparams;
  tsourcep=params.tsourceparams;
  numaxBinary=tsourcep.numax;
  factBinary=tsourcep.fact;
  nu = *((REAL8 *)params.extraparams);
  /*frequency in the system frame*/
  nuz=(1.+z)*nu;
  /*single spectral energy density in the system frame*/
  SDensityBinary(&dEgw, nuz);
  /*coalescence rate*/
  SCR (s->statusPtr,&vRz, z,&params);
  Ez(&vEz, z, cosmop);
  *domegaz = dEgw*vRz/((1.+z)*vEz);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }
