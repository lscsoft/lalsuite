/*
*  Copyright (C) 2007 David Churches, Jolien Creighton, David McKechan, B.S. Sathyaprakash, Thomas Cokelaer, Duncan Brown, Riccardo Sturani, Laszlo Vereb, Drew Keppel
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
\author Sathyaprakash, B. S.
\file
\ingroup LALSimInspiral_h

\brief Module containing the energy and flux functions for waveform generation.

\heading{Prototypes}

<tt>REAL8 SAMPLEFUNCTION()</tt>
<ul>
<li> \c v: PN parameter.<\li>
<li> \c ak: Input containing the PN expnasion coefficients.</li>
</ul>

\heading{Description}
This module gives the post-Newtonian expansions and/or P-approximants
to the energy, its derivative and gravitational-wave flux functions. More
specifically, the <tt>REAL8</tt> functions below give Taylor expansions
of \f$dE/dv,\f$ and \f${\cal F}(v),\f$ P-approximants of \f$e(v),\f$ \f$dE/dv\f$
(derived from \f$e(v)\f$) and \f${\cal F}(v).\f$

\heading{Notes}
<ul>
<li> See Damour, Iyer and Sathyaprakash, PRD 57, 885, 1998 for further details.
Damour, Iyer and Sathyaprakash, PRD 63, 044023, 2001 is a resource paper that
summarizes how to generate waveforms in different approximations to the dynamics
of a compact binary under radiation reaction.</li>
<li> The Pade Approximant for the 1PN expansion is undefined as also
EOB at orders less than 2PN. BCV is independent of the PN order.
Spinning waveforms are only defined at the highest PN order.</li>
</ul>

*/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiraldEnergyFlux.h>



REAL8 Et0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy;
   dEnergy = ak->ETaN * v;
   return (dEnergy);
}

REAL8 Et2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 Energy;
   REAL8 x = v*v;
   Energy = ak->ETaN * x * (1. + ak->ETa1*x);
   return Energy;
}

REAL8 Et4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 Energy;
   REAL8 x = v*v;
   REAL8 x2 = x*x;
   Energy = ak->ETaN * x * (1. + ak->ETa1*x + ak->ETa2*x2);
   return Energy;
}

REAL8 Et6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 Energy;
   REAL8 x = v*v;
   REAL8 x2 = x*x;
   REAL8 x3 = x*x2;
   Energy = ak->ETaN * x * (1. + ak->ETa1*x + ak->ETa2*x2 + ak->ETa3*x3);
   return Energy;
}

REAL8 dEt0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy;
   dEnergy = ak->dETaN * v;
   return (dEnergy);
}


REAL8 dEt2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x);
   return (dEnergy);
}


REAL8 dEt4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa2*x*x);
   return (dEnergy);
}


REAL8 dEt6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa2*x*x + ak->dETa3*x*x*x);
   return (dEnergy);
}


REAL8 Ft0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10;
   return (flux);
}


REAL8 Ft2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2);
   return (flux);
}


REAL8 Ft3(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2 + ak->FTa3*v2*v);
   return (flux);
}


REAL8 Ft4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4);
   return (flux);
}


REAL8 Ft5(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v);
   return (flux);
}


REAL8 Ft6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(v))*v6);
   return (flux);
}


REAL8 Ft7(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(v))*v6 + ak->FTa7*v6*v);
   return (flux);
}


/*
REAL8 ep0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 x, energy;
   ak = NULL;
   x = v*v;
   energy = -x;
   return (energy);
}
*/


REAL8 ep2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1 * x);
   return (energy);
}


REAL8 ep4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1*x /(1. + ak->ePa2*x));
   return (energy);
}


REAL8 ep6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1*x /(1. + ak->ePa2*x /(1. + ak->ePa3*x)));
   return (energy);
}


/*
REAL8 dEp0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, y;
   energy = ep0(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   denergy = -1;
   dEnergy = v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}
*/


REAL8 dEp2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, x, y;
   x = v*v;
   energy = ep2(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   denergy = -1. / ((1. + ak->ePa1*x)*(1. + ak->ePa1*x));
   dEnergy = v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}


REAL8 dEp4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, denom, x, y;
   x = v*v;
   energy = ep4(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;

   denom = 1. + (ak->ePa1 + ak->ePa2) * x;
   denom = denom * denom;
   denergy = (1. + 2.*ak->ePa2*x + ((ak->ePa1 + ak->ePa2) * ak->ePa2 * x*x))/denom;
   dEnergy = - v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}


REAL8 dEp6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, denom, x, y;
   x = v*v;
   energy = ep6(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   
   denom = 1. + (ak->ePa1 + ak->ePa2 + ak->ePa3) * x + ak->ePa1*ak->ePa3*x*x;
   denom = denom * denom;

   denergy = (1. + 2.*(ak->ePa2+ak->ePa3)*x + (ak->ePa1*ak->ePa2
           + ak->ePa2*ak->ePa2 + 2.* ak->ePa2*ak->ePa3
           + ak->ePa3*ak->ePa3) * x*x)
           /denom;
   dEnergy = - v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}


/*
REAL8 Fp0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10;
   return (flux);
}
*/


/*
REAL8 Fp1(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v) * (1.-v/ak->vpoleP4));
   return (flux);
}
*/


/*
REAL8 Fp2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v / (1.+ak->fPa2*v)) * (1.-v/ak->vpoleP4));
   return (flux);
}
*/


REAL8 Fp3(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v)))
	* (1.-v/ak->vpoleP4));
   return (flux);
}


REAL8 Fp4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
	/ (1.+ak->fPa4*v)))) * (1.-v/ak->vpoleP4));
   return (flux);
}


REAL8 Fp5(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
	/ (1.+ak->fPa4*v / (1.+ak->fPa5*v))))) * (1.-v/ak->vpoleP4));
   return (flux);
}


REAL8 Fp6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v))))))
        * (1.-v/ak->vpoleP6));
   /* */
   flux *= (1.+  log(v/ak->vlsoP4) * ak->FTl6*v6) ;
   return (flux);
}


REAL8 Fp7(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v / (1.+ak->fPa7*v)))))))
        * (1.-v/ak->vpoleP6));
   flux *= (1.+  log(v/ak->vlsoP4) * ak->FTl6*v6) ;
   return (flux);
}


/* Flux for the EOBNRv2 model */
REAL8 Fp8PP(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10, l6, l8;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   l6 = ak->FTl6;
   l8 = ak->FTl8 - ak->FTa2*ak->FTl6;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v / (1.+ak->fPa7*v
        / (1.+ak->fPa8*v))))))))
        * (1.-v/ak->vpolePP));
   flux *= (1.+  log(v/ak->vlsoPP) * (l6*v6 + l8*v8) ) ;
   return (flux);
}


/*
REAL8 Fp8(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10, l6, l8;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   l6 = ak->FTl6;
   l8 = ak->FTl8 - ak->FTa2*ak->FTl6;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v / (1.+ak->fPa7*v
	/ (1.+ak->fPa8*v))))))))
        * (1.-v/ak->vpoleP6));
   flux *= (1.+  log(v/ak->vlsoP4) * (l6*v6 + l8*v8) ) ;
   return (flux);
}
*/
