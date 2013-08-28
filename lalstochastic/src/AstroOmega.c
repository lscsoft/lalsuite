/*
*  Copyright (C) 2007 Robert Adam Mercer, Tania Regimbau
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


/**
 * \file
 * \author Regimbau T.
 * \ingroup AstroOmega_h
 *
 * \brief Compute the energy density spectrum of stochastic backgrounds produced
 * by cosmological population of astrophysical sources.
 *
 * \heading{Prototypes}
 *
 * <tt>LALAstroOmega()</tt>
 * \heading{Description}
 * The function of this module computes the energy density parameter \f$\Omega_{gw}(\nu_{o})\f$ for a given source and a given cosmological model.
 * The spectral properties of the stochastic background are characterized by the dimensionless parameter
 * \f{equation}{
 * \Omega _{gw}=\frac{1}{\rho _{c}}\frac{d\rho _{gw}}{d\log \nu _{o}}
 * \f}
 * where \f$\rho _{gw}\f$ is the gravitational energy density, \f$\nu _{o}\f$ the frequency in the
 * observer frame and
 * \f{equation}{
 * \rho _{c}=\frac{3H_{o}^{2}}{8\pi G}
 * \f}
 * the critical energy density to close the Universe today.
 *
 * For a given cosmological model, a stochastic background of
 * astrophysical origin is fully determined by the source formation rate as a function of redshift,
 * \f$dR(z)\f$ and the gravitational spectral density of a single source, \f$\frac{dE_{gw}}{d\nu}\f$.
 * Following Ferrari et al. (1999):
 * \f{equation}{
 * \Omega _{gw}(\nu _{o})=\frac{1}{c^{3} \rho _{c}}{\nu _{o}}F_{\nu _{o}}
 * \f}
 * where
 * \f{equation}{
 * F_{\nu _{o}}=\int_{0}^{z_{\max }}f_{\nu _{o}}dR(z)
 * \f}
 * is the gravitational wave flux at the
 * frequency \f$\nu _{o}\f$ integrated over all the sources.
 * The gravitational flux of a single source located between z, z+dz is:
 * \f{equation}{
 * f_{\nu _{o}}=\frac{1}{4\pi d_{L}^{2}}\frac{dE_{gw}}{d\nu }(1+z)
 * \f}
 * where \f$d_{L}=(1+z)r\f$ is the distance luminosity and  \f$\nu=(1+z)\nu _{o}\f$
 * the frequency in the source frame.
 * The event rate between z, z+dz as observed in our frame is given by:
 * \f{equation}{
 * dR(z)=\lambda _{p}\frac{R_{c}(z)}{1+z}{\frac{{dV}}{{dz}}}dz
 * \f}
 * where \f$R_{c}(z)\f$ is the cosmic star formation rate and  \f$\lambda_{p}\f$
 * is the mass fraction of the source progenitors. The
 * term (1 +z) in the denominator accounts for the time dilatation by cosmic expansion of the
 * observed rate.
 * The element of the comoving volume is
 * \f{equation}{
 * dV = 4\pi r^2 \frac{{c}}{{H_o}} \frac{{dz}}{{E(\Omega_i,z)}}
 * \f}
 * where the function E\f$(\Omega_i,z)\f$ is
 * defined by the equation:
 * \f{equation}{
 * E(\Omega_i,z) =[\Omega_m(1+z)^3 + \Omega_v]^{1/2}
 * \f}
 * where \f$\Omega_m\f$ and \f$\Omega_v\f$ are respectively the density parameters
 * due to matter (baryonic and non-baryonic) and the vacuum.
 * The cosmic star formation rate is computed as:
 * \f{equation}{
 * R_{c}(z) =
 * R_{SF2}(z)
 * \frac{h_0}{0.65}\frac{E(\Omega_i,z)}{(1+z)^{3/2}}\;
 * \mathrm{M}_{\odot}\,\mathrm{yr}^{-1}\,\mathrm{Mpc}^{-3}
 * \f}
 * where
 * \f{equation}{
 * R_{SF2}(z) = \frac{0.15\,e^{3.4z}}{(22 + e^{3.4z})} \;
 * \mathrm{M}_{\odot}\,\mathrm{yr}^{-1}\,\mathrm{Mpc}^{-3}
 * \f}
 * is the cosmic star formation rate (Madau \&Porciani, 2001) in a matter-dominateduniverse (\f$\Omega_{m}=1\f$) with \f$H_o = 65\,\mathrm{km\,s}^{-1}\mathrm{Mpc}^{-1}\f$.
 * Combining the previous equations one obtains:
 * \f{equation}{
 * \Omega _{gw}(\nu _{o})=\frac{8\pi G}{3c^{2}H_{o}^{3}}\lambda_{p}\, \nu_{o}\int_{0}^{z_{\sup }}\frac{dE_{gw}}{d\nu }\frac{R_{c}(z)}{E(z)(1+z)^2}dz
 * \f}
 * The upper limit of the integral is determined by the cutoff frequency in the
 * source frame, \f$\nu _{\sup }\f$, as:
 * \f{equation}{
 * z_{\sup }=\frac{\nu _{\sup }}{\nu_{o}} - 1
 * \f}
 * Note that we are limited to \f$z<5\f$ since the cosmic star formation rate is not
 * modeled above. This restriction is without significant consequence
 * since objects with \f$z>5 \f$ contribute very little to the integrated signal.
 *
 * \heading{Operating Instructions}
 * the following program shows how to use the function LALAstroOmega
 *
 * \code
 *
 * #include <stdio.h>
 * #include <math.h>
 * #include <lal/LALConfig.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/Integrate.h>
 * #include <lal/AstroOmega.h>
 *
 * //define here the spectral energy density of a single source
 *
 * static void SDensity (REAL8 *dEgw, REAL8 nu)
 * {
 *   *dEgw=pow(nu,3.);
 *   return;
 * }
 *
 * int main ()
 * {
 *   static LALStatus status;
 *   DIntegrateIn  zint;
 *   AstroOmegaGeneralParams params;
 *   AstroOmegaSourceParams sourcep;
 *   AstroOmegaCosmoParams cosmop;
 *   REAL8 omegaz, nu, test;
 *
 *   //define here the model parameters
 *   //cosmological parameters
 *   cosmop.ho=0.65;
 *   cosmop.density_matter=0.3;
 *   cosmop.density_vacuum=0.7;
 *   cosmop.density_k=0.;
 *
 *   //source parameters for rotating pulsars
 *   sourcep.fact = 2.88E-22;
 *   sourcep.numax = 4000.;
 *   sourcep.SDensitySource = SDensity;
 *   params.cosmoparams = cosmop;
 *   params.sourceparams = sourcep;
 *   for (nu = 0.; nu < sourcep.numax ; nu = nu + 10.)
 *   {
 *     params.extraparams = &nu;
 *     LALAstroOmegaSource (&status, &test, nu, &params);
 *     printf("omega(%f)= %.2e\n", nu, test);
 *   }
 *   return 0;
 * }
 *
 * \endcode
 *
 */

#include <stdio.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/Integrate.h>
#include <lal/AstroOmega.h>

static void SFR (REAL8 *result, REAL8 z);
static void dAstroOmega (LALStatus *s, REAL8 *domegaz, REAL8 z, void *p);

void LALAstroOmega (LALStatus *s, REAL8 *omeganu, REAL8 nu, void *p)
 {
   DIntegrateIn  zint;
   AstroOmegaParams params;
   AstroOmegaSourceParams sourcep;
   AstroOmegaCosmoParams cosmop;
   REAL8 omegaz, zmax, numax, lambda;
   INITSTATUS(s);
   ATTATCHSTATUSPTR (s);

   params = *((AstroOmegaParams *)p);
   cosmop=params.cosmoparams;
   sourcep = params.sourceparams;
   numax= sourcep.numax;
   lambda = sourcep.lambda;

   if((nu >= numax)||(nu <= 0.)){*omeganu = 0.;}
   else
    {
     if (nu < (numax / 6.)) {zmax = 5.;}
     else {zmax = (numax / nu) - 1.;}

     zint.function = dAstroOmega;
     zint.xmin     = 0;
     zint.xmax     = zmax;
     zint.type     = ClosedInterval;

     LALDRombergIntegrate (s->statusPtr, &omegaz, &zint, &params);
     *omeganu = 4.66e-56 * lambda / (cosmop.ho * cosmop.ho) * nu * omegaz;
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
  /*AstroOmegaCosmoParams cosmop;*/
  REAL8LALSDensity *SDensitySource;
  REAL8 Rc, dEgw, nu, nuz;

  INITSTATUS(s);
  ATTATCHSTATUSPTR (s);
  params = *((AstroOmegaParams *)p);
  sourcep = params.sourceparams;
  SDensitySource = sourcep.SDensitySource;
  nu = *((REAL8 *)params.extraparams);

  /*frequency in the source frame*/
  nuz = (1. + z) * nu;
  /*single spectral energy density in the source frame*/
  SDensitySource(&dEgw, nuz);
  /*cosmic formation rate*/
  SFR(&Rc, z);
  *domegaz = dEgw * Rc / pow((1.+z),3.5);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
 }
