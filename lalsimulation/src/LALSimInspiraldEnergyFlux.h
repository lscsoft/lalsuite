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

#ifndef _LALSIMINSPIRALDENERGYFLUX_H
#define _LALSIMINSPIRALDENERGYFLUX_H

/* remove SWIG interface directives */
#if !defined(SWIG) && !defined(SWIGLAL_STRUCT)
#define SWIGLAL_STRUCT(...)
#endif

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

typedef struct
tagexpnCoeffsdEnergyFlux {
  SWIGLAL_STRUCT(expnCoeffsdEnergyFlux);

   /* coefficients in the Pade expression of new energy function */
   REAL8 ePaN, ePa1, ePa2, ePa3;
   /* coefficients in the Taylor expansion of usual energy function */
   REAL8 ETaN, ETa1, ETa2, ETa3;
   /* coefficients in the Taylor expansion of the derivative of the
    usual energy function */
   REAL8 dETaN, dETa1, dETa2, dETa3;

   /* Taylor expansion coefficients of energy flux*/
   REAL8 FTaN, FTa1, FTa2, FTa3, FTa4, FTa5, FTa6, FTa7, FTa8, FTl6, FTl8;
   /* Coefficients of the corresponding P-approximant */
   REAL8 fPaN, fPa1, fPa2, fPa3, fPa4, fPa5, fPa6, fPa7, fPa8;

   /* symmetric mass ratio, total mass, component masses */
   REAL8 eta;

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


typedef REAL8 EnergyFunction(
   REAL8 v,
   expnCoeffsdEnergyFlux *ak);


typedef REAL8 FluxFunction(
   REAL8 v,
   expnCoeffsdEnergyFlux *ak);


typedef struct
tagTofVIntegrandIn
{
  SWIGLAL_STRUCT(TofVIntegrandIn);

   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffsdEnergyFlux *coeffs;
} TofVIntegrandIn;


typedef struct
tagTofVIn
{
  SWIGLAL_STRUCT(TofVIn);

   REAL8 t;
   REAL8 v0;
   REAL8 t0;
   REAL8 vlso;
   REAL8 totalmass;
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffsdEnergyFlux *coeffs;
} TofVIn;


REAL8
XLALSimInspiralTofV (
   REAL8 v,
   void *params
   );


REAL8
XLALSimInspiralTofVIntegrand(
   REAL8      v,
   void      *params
   );


REAL8 Et0(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Et2(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Et4(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Et6(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 dEt0(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 dEt2(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 dEt4(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 dEt6(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Ft0(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Ft2(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Ft3(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Ft4(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Ft5(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Ft6(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Ft7(REAL8 v, expnCoeffsdEnergyFlux *ak);


/*
REAL8 ep0(REAL8 v, expnCoeffsdEnergyFlux *ak);
*/


REAL8 ep2(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 ep4(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 ep6(REAL8 v, expnCoeffsdEnergyFlux *ak);


/*
REAL8 dEp0(REAL8 v, expnCoeffsdEnergyFlux *ak);
*/


REAL8 dEp2(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 dEp4(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 dEp6(REAL8 v, expnCoeffsdEnergyFlux *ak);


/*
REAL8 Fp0(REAL8 v, expnCoeffsdEnergyFlux *ak);
*/

/*
REAL8 Fp1(REAL8 v, expnCoeffsdEnergyFlux *ak);
*/

/*
REAL8 Fp2(REAL8 v, expnCoeffsdEnergyFlux *ak);
*/


REAL8 Fp3(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Fp4(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Fp5(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Fp6(REAL8 v, expnCoeffsdEnergyFlux *ak);


REAL8 Fp7(REAL8 v, expnCoeffsdEnergyFlux *ak);


/* Flux for the EOBNRv2 model */
REAL8 Fp8PP(REAL8 v, expnCoeffsdEnergyFlux *ak);

/*
REAL8 Fp8(REAL8 v, expnCoeffsdEnergyFlux *ak);
*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRALDENERGYFLUX_H */
