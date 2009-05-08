/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   mcmc_nolal.c:              dummy routines to compile and run SPINspiral without LAL
   
   
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



#include <mcmc.h>


//Use the LAL 3.5/2.5 PN spinning waveform, with 1 spinning object (12 parameters)
void templateLAL12(struct parset *par, struct interferometer *ifo[], int ifonr)
{
  //Get rid of warnings at compile time (x set but never used)
  double x = 0.0;
  x = par->par[0];
  x = ifo[0]->lati; 
  x = (double)ifonr;
  x = x;
  
  printf("\n\n");
  printf("\n    *** ERROR ***");
  printf("\n    You have compiled SPINspiral without LAL support, but are trying to run it using a LAL waveform.");
  printf("\n    Please change mcmcWaveform or injectionWaveform to use a non-LAL waveform  *or*  recompile with LAL support.");
  printf("\n\n\n");
  exit(1);
}

//Use the LAL 3.5/2.5 PN spinning waveform, with 2 spinning objects (15 parameters)
void templateLAL15(struct parset *par, struct interferometer *ifo[], int ifonr)
{
  //Get rid of warnings at compile time (x set but never used)
  double x = 0.0;
  x = par->par[0];
  x = ifo[0]->lati; 
  x = (double)ifonr;
  x = x;
  
  printf("\n\n");
  printf("\n    *** ERROR ***");
  printf("\n    You have compiled SPINspiral without LAL support, but are trying to run it using a LAL waveform.");
  printf("\n    Please change mcmcWaveform or injectionWaveform to use a non-LAL waveform  *or*  recompile with LAL support.");
  printf("\n\n\n");
  exit(1);
}

