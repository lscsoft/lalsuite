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

#ifndef LALSIMINSPIRALDENERGYFLUC_C
#define LALSIMINSPIRALDENERGYFLUC_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <math.h>
#include <lal/LALStdlib.h>

typedef struct
{
   /* coefficients in the Pade expression of new energy function */
   REAL8 ePaN, ePa1, ePa2, ePa3;
   /* coefficients in the Taylor expansion of usual energy function */
   REAL8 ETaN, ETa1, ETa2, ETa3, ETa5, ETa6;
   /* coefficients in the Taylor expansion of the derivative of the
    usual energy function */
   REAL8 dETaN, dETa1, dETa2, dETa3, dETa5, dETa6;

   /* Taylor expansion coefficients of energy flux*/
   REAL8 FTaN, FTa1, FTa2, FTa3, FTa4, FTa5, FTa6, FTa7, FTa8, FTl6, FTl8, FTa10, FTa12;
   /* Coefficients of the corresponding P-approximant */
   REAL8 fPaN, fPa1, fPa2, fPa3, fPa4, fPa5, fPa6, fPa7, fPa8;

   /* symmetric mass ratio, total mass, component masses */
   REAL8 eta, totalmass;

   /* initial and final values of frequency, time, velocity; lso
    values of velocity and frequency; final phase. */
   REAL8 vlso;

   /* last stable orbit and pole defined by various Taylor and P-approximants */
   REAL8 vlsoT0, vlsoT2, vlsoT4, vlsoT6;
   REAL8 vlsoP0, vlsoP2, vlsoP4, vlsoP6;
   REAL8 vlsoPP;
   REAL8 vpoleP4, vpoleP6;
   REAL8 vpolePP;
}  expnCoeffsdEnergyFlux;


typedef REAL8 (*EnergyFunction)(
   REAL8 v,
   expnCoeffsdEnergyFlux *ak);


typedef REAL8 (*FluxFunction)(
   REAL8 v,
   expnCoeffsdEnergyFlux *ak);


typedef struct
{
   EnergyFunction dEnergy;
   FluxFunction flux;
   expnCoeffsdEnergyFlux *coeffs;
} TofVIntegrandIn;


typedef struct
{
   REAL8 t;
   REAL8 v0;
   REAL8 t0;
   REAL8 vlso;
   REAL8 totalmass;
   EnergyFunction dEnergy;
   FluxFunction flux;
   expnCoeffsdEnergyFlux *coeffs;
} TofVIn;


static REAL8 UNUSED XLALSimInspiralEt0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy;
   REAL8 x = v*v;
   dEnergy = ak->ETaN * x;
   return (dEnergy);
}

static REAL8 UNUSED XLALSimInspiralEt2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 Energy;
   REAL8 x = v*v;
   REAL8 x5 = x*x*x*x*x;
   REAL8 x6 = x*x5;
   Energy = ak->ETaN * x * (1. + ak->ETa1*x + ak->ETa5*x5 + ak->ETa6*x6);
   return Energy;
}

static REAL8 UNUSED XLALSimInspiralEt4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 Energy;
   REAL8 x = v*v;
   REAL8 x2 = x*x;
   REAL8 x5 = x*x2*x2;
   REAL8 x6 = x*x5;
   Energy = ak->ETaN * x * (1. + ak->ETa1*x + ak->ETa2*x2 + ak->ETa5*x5
          + ak->ETa6*x6);
   return Energy;
}

static REAL8 UNUSED XLALSimInspiralEt6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 Energy;
   REAL8 x = v*v;
   REAL8 x2 = x*x;
   REAL8 x3 = x*x2;
   REAL8 x5 = x3*x2;
   REAL8 x6 = x*x5;
   Energy = ak->ETaN * x * (1. + ak->ETa1*x + ak->ETa2*x2 + ak->ETa3*x3
          + ak->ETa5*x5 + ak->ETa6*x6);
   return Energy;
}

static REAL8 UNUSED XLALSimInspiraldEt0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy;
   dEnergy = ak->dETaN * v;
   return (dEnergy);
}


static REAL8 UNUSED XLALSimInspiraldEt2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy;
   REAL8 x = v*v;
   REAL8 x5 = x*x*x*x*x;
   REAL8 x6 = x*x5;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa5*x5 + ak->dETa6*x6);
   return (dEnergy);
}


static REAL8 UNUSED XLALSimInspiraldEt4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy;
   REAL8 x = v*v;
   REAL8 x2 = x*x;
   REAL8 x5 = x*x2*x2;
   REAL8 x6 = x*x5;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa2*x2 + ak->dETa5*x5
           + ak->dETa6*x6);
   return (dEnergy);
}


static REAL8 UNUSED XLALSimInspiraldEt6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 dEnergy;
   REAL8 x = v*v;
   REAL8 x2 = x*x;
   REAL8 x3 = x*x2;
   REAL8 x5 = x3*x2;
   REAL8 x6 = x*x5;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa2*x2 + ak->dETa3*x3
           + ak->dETa5*x5 + ak->dETa6*x6);
   return (dEnergy);
}


static REAL8 UNUSED XLALSimInspiralFt0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10;
   return (flux);
}


static REAL8 UNUSED XLALSimInspiralFt2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10,v12;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   v12 = v10*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2 + ak->FTa10*v10 + ak->FTa12*v12);
   return (flux);
}


static REAL8 UNUSED XLALSimInspiralFt3(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10,v12;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   v12 = v10*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa10*v10
        + ak->FTa12*v12);
   return (flux);
}


static REAL8 UNUSED XLALSimInspiralFt4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10,v12;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   v12 = v10*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa10*v10 + ak->FTa12*v12);
   return (flux);
}


static REAL8 UNUSED XLALSimInspiralFt5(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v8,v10,v12;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   v12 = v10*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v + ak->FTa10*v10 + ak->FTa12*v12);
   return (flux);
}


static REAL8 UNUSED XLALSimInspiralFt6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10,v12;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   v12 = v10*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(16.0*v2))*v6 + ak->FTa10*v10
        + ak->FTa12*v12);
   return (flux);
}


static REAL8 UNUSED XLALSimInspiralFt7(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10,v12;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   v12 = v10*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(16.0*v2))*v6 + ak->FTa7*v6*v
        + ak->FTa10*v10 + ak->FTa12*v12);
   return (flux);
}


/*
static REAL8 UNUSED XLALSimInspiralep0(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 x, energy;
   ak = NULL;
   x = v*v;
   energy = -x;
   return (energy);
}
*/


static REAL8 UNUSED XLALSimInspiralep2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1 * x);
   return (energy);
}


static REAL8 UNUSED XLALSimInspiralep4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1*x /(1. + ak->ePa2*x));
   return (energy);
}


static REAL8 UNUSED XLALSimInspiralep6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1*x /(1. + ak->ePa2*x /(1. + ak->ePa3*x)));
   return (energy);
}


/*
static REAL8 UNUSED XLALSimInspiraldEp0(REAL8 v, expnCoeffsdEnergyFlux *ak)
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


static REAL8 UNUSED XLALSimInspiraldEp2(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, x, y;
   x = v*v;
   energy = XLALSimInspiralep2(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   denergy = -1. / ((1. + ak->ePa1*x)*(1. + ak->ePa1*x));
   dEnergy = v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}


static REAL8 UNUSED XLALSimInspiraldEp4(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, denom, x, y;
   x = v*v;
   energy = XLALSimInspiralep4(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;

   denom = 1. + (ak->ePa1 + ak->ePa2) * x;
   denom = denom * denom;
   denergy = (1. + 2.*ak->ePa2*x + ((ak->ePa1 + ak->ePa2) * ak->ePa2 * x*x))/denom;
   dEnergy = - v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}


static REAL8 UNUSED XLALSimInspiraldEp6(REAL8 v, expnCoeffsdEnergyFlux *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, denom, x, y;
   x = v*v;
   energy = XLALSimInspiralep6(v, ak);
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
static REAL8 UNUSED XLALSimInspiralFp0(REAL8 v, expnCoeffsdEnergyFlux *ak)
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
static REAL8 UNUSED XLALSimInspiralFp1(REAL8 v, expnCoeffsdEnergyFlux *ak)
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
static REAL8 UNUSED XLALSimInspiralFp2(REAL8 v, expnCoeffsdEnergyFlux *ak)
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


static REAL8 UNUSED XLALSimInspiralFp3(REAL8 v, expnCoeffsdEnergyFlux *ak)
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


static REAL8 UNUSED XLALSimInspiralFp4(REAL8 v, expnCoeffsdEnergyFlux *ak)
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


static REAL8 UNUSED XLALSimInspiralFp5(REAL8 v, expnCoeffsdEnergyFlux *ak)
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


static REAL8 UNUSED XLALSimInspiralFp6(REAL8 v, expnCoeffsdEnergyFlux *ak)
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


static REAL8 UNUSED XLALSimInspiralFp7(REAL8 v, expnCoeffsdEnergyFlux *ak)
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
static REAL8 UNUSED XLALSimInspiralFp8PP(REAL8 v, expnCoeffsdEnergyFlux *ak)
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
static REAL8 UNUSED XLALSimInspiralFp8(REAL8 v, expnCoeffsdEnergyFlux *ak)
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

#endif /* LALSIMINSPIRALDENERGYFLUC_C */
