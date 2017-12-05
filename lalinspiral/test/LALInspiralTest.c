/*
*  Copyright (C) 2007 Stas Babak, Bernd Machenschalk, David Churches, Duncan Brown, David Chin, Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief Test routine for codes that generate inspiral waveform from non-spinning black
 * hole binaries.
 *
 * Time domain signals are returned when the \c approximant is
 * one of <tt>TaylorT1, TaylorT2, TaylorT3, PadeT1, EOB, SpinTaylorT3</tt>
 * and frequency domain signals are returned when the \c approximant is
 * one of <tt>TaylorF1, TaylorF2, BCV.</tt> This code checks every available approximant
 * at every order and reports whether or not there was any problem with the
 * generation codes.
 *
 * To generate a waveform first set the \c InspiralTemplate structure (see
 * below for an example).  Next, to measure the length of the array required
 * call the function
 * \code
 * LALInspiralWaveLength (&status, &n, params)
 * \endcode
 * The length will be returned in \c n. Finally, call the function
 * \code
 * LALInspiralWave(&status, signal1, params);
 * \endcode
 * to generate the wave, which will be returned in \c signal.
 *
 * Example values of the parameters that can be set (with options in brackets) is:
 * \code
 * params.OmegaS = 0.;     (Unknown 3PN parameter in energy; shown to be 0 by DJS)
 * params.Theta = 0.;      (Unknown 3PN parameter in flux; arbitrarily set to 0)
 * params.ieta=1;          (1 for comparable masses model, 0 for test mass model)
 * params.mass1=1.4;       (masses of the component stars in solar masses)
 * params.mass2=1.4;
 * params.startTime=0.0;   (defined so that the instantaneous GW frequency
 * is params.fLower at params.startTime)
 * params.startPhase=0.0;  (0 to LAL_PI_2)
 * params.fLower=40.0;     (in Hz)
 * params.fCutoff=1000.0;  (in Hz)
 * params.tSampling=4000.; (in Hz; should be larger than 2 fCutoff or 2 flso,
 * whichever is smaller)
 * params.signalAmplitude=1.0;
 * params.nStartPad=0;     (number of leading zero bins)
 * params.nEndPad=0;       (number of trailing zero bins)
 * params.approximant=TaylorF2; (TaylorT1, PadeT1=ODE solver;
 * TaylorT2=implicit phasing formula solved in quadrature;
 * TaylorT3=explicit time-domain phasing;
 * TaylorF1=stationary phase approx. using ODEs;
 * TaylorF2=usual stationary phase approx.;
 * EOB=effective-one-body approach)
 * params.order=twoPN;     (also newtonian, onePN, oneAndHalfPN, twoPN,
 * twoAndHalfPN, threePN, threeAndHalfPN)
 * params.massChoice=m1Andm2; (also t0t2, t0t3, t0t4, totalMassAndEta,totalMassAndMu)
 * params.psi0 = 132250.;   (parameter required to generate BCV detection templates)
 * params.psi3 = -1014.2;   (parameter required to generate BCV detection templates)
 * params.alpha = 0.528;    (amplitude correction used in BCV templates)
 * params.fFinal = 868.7;  (frequency at which the BCV template is terminated)
 *
 * \endcode
 *
 */

#include <stdio.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>


void printf_timeseries (int n, float *sig, double delta, double t0) ;
void printf_timeseries (int n, float *sig, double delta, double t0)
{
  int i=0;
  FILE *outfile1;

  outfile1=fopen("wave1.dat","a");

  do
     fprintf (outfile1,"%e %e\n", i*delta+t0, *(sig+i));
  while (n-++i);

  fprintf(outfile1,"&\n");
  fclose(outfile1);
}


int main (void) {
   static REAL4Vector *signal1;
   static LALStatus status;
   static InspiralTemplate params;
   static REAL8 dt;
   UINT4 n;

   params.OmegaS = 0.;
   params.Theta = 0.;
   params.ieta=1;
   params.mass1=5.1;
   params.mass2=5.;
   params.startTime=0.0;
   params.startPhase=0.0;
   params.fLower=40.0;
   params.fCutoff=2048.00;
   params.tSampling=8192.0;
   params.order=6;
   params.approximant=IMRPhenomB;
/* SpinTaylorT3 */
   params.signalAmplitude=1.0;
   params.nStartPad=0;
   params.nEndPad=0;
   params.massChoice=m1Andm2;
   params.distance = 100.;
   dt = 1./params.tSampling;
   params.spin1[2]=1.0;
   params.spin2[2]=0.99;

   LALInspiralWaveLength(&status, &n, params);
   LALInspiralParameterCalc(&status, &params);
   fprintf(stderr, "Testing Inspiral Signal Generation Codes:\n");
   fprintf(stderr, "Signal length=%d, m1=%e, m2=%e, fLower=%e, fUpper=%e\n",
		   n, params.mass1, params.mass2, params.fLower, params.fCutoff);
   LALCreateVector(&status, &signal1, n);

   LALInspiralWave(&status, signal1, &params);
   fprintf(stderr,"fFinal = %e\n", params.fFinal);
	      printf_timeseries(signal1->length, signal1->data, dt, params.startTime);
REPORTSTATUS(&status);

   printf("duration is %f \n", params.tC);
   printf("final frequency is %f \n", params.fFinal);

   LALDestroyVector(&status, &signal1);
   LALCheckMemoryLeaks();
   return 0;
}
