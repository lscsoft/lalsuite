/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   include/mcmc_lal.h:        LAL-interface header file
   
   
   Copyright 2007, 2008, 2009 Christian Roever, Marc van der Sluys, Vivien Raymond, Ilya Mandel
   
   
   This file is part of SPINspiral.
   
   SPINspiral is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   SPINspiral is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with SPINspiral.  If not, see <http://www.gnu.org/licenses/>.
   
*/



#ifndef mcmc_lal_h
#define mcmc_lal_h


#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>

#include <mcmc.h>


// Declare function prototypes:
//void LALHpHc(CoherentGW *waveform, double *hplus, double *hcross, int *l, int length, struct parset *par, struct interferometer *ifo, int ifonr);
void LALHpHc12(LALStatus *status, CoherentGW *waveform, SimInspiralTable *injParams, PPNParamStruc *ppnParams, int *l, struct parset *par, struct interferometer *ifo);
void LALHpHc15(LALStatus *status, CoherentGW *waveform, SimInspiralTable *injParams, PPNParamStruc *ppnParams, int *l, struct parset *par, struct interferometer *ifo);
//double LALFpFc(CoherentGW *waveform, double *wave, int *l, int length, struct parset *par, int ifonr);
double LALFpFc(LALStatus *status, CoherentGW *waveform, SimInspiralTable *injParams, PPNParamStruc *ppnParams, double *wave, int length, struct parset *par, struct interferometer *ifo, int ifonr);
void LALfreedom(CoherentGW *waveform);


#endif
