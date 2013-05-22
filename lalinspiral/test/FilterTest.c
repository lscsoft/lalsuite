/*
*  Copyright (C) 2007 Bernd Machenschalk, David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Anand Sengupta, Thomas Cokelaer
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
\author Cokelaer, T. and Sathyaprakash, B. S.
\file

\brief Test code for the inspiral bank modules.

\heading{Usage}
\code
BankEfficiency [options]

The options are :
   -alpha : BCV amplitude correction parameter
   -approximant : Post-Newtonian model such as TaylorT1, PadeT1, EOB, BCV ...
-fl : lower frequency cutoff
-mMin : minimal mass of component stars
-mMax : maximal mass of component stars
-mm : minimal match for template bank
-n : number of trials
-order : order of PN model
-quiet : if this flag is present, the output is restricted to the min
-seed : seed for random generation
-sigAmp : amplitude of the signal
-simType : type of simulation, 0, 1 or 2
-x0Max : Max value of psi0
-x1Min : Min value of psi


\endcode

\heading{Description}

This test code gives an example of how one might generate inspiral
waveforms and use them to compute the overlap of a random signal
(with or witnout simulated noise) of a certain strength. The parameter
Note that one must calculate the length of the waveform and allocate memory for it
\e before calling
\c InspiralWave. The length of the waveform can be calculated by calling the function
\c InspiralWaveLength beforehand, as shown.

There are only two functions which one can call to generate waveforms. These are
\c InspiralWave,
which will return a \e single waveform, and \c InspiralWaveTemplates, which
returns a \e pair
of waveforms which have phases which differ by \f$\pi/2\f$.

\heading{Uses}
This code directly uses the following functions (see those functions
to find out what they call in turn):
\code
LALInspiralWaveLength
LALInspiralCreateCoarseBank
LALRandomInspiralSignal
LALInspiralParameterCalc
LALNoiseSpectralDensity
LALCreateForwardRealFFTPlan
LALCreateReverseRealFFTPlan
LALForwardRealFFT
LALReverseRealFFT
LALDestroyRealFFTPlan
LALInspiralWaveOverlap
LALInspiralParameterCalc
\endcode

\heading{Notes}


*/

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

void printf_timeseries (INT4 n, REAL4 *signal1, REAL8 delta, REAL8 t0);


int
main (  int argc, char **argv )
{
   static LALStatus status;
   static INT4 i, approx, tmplt;
   static UINT4 j;
   static UINT4 psdLength, quietFlag = 0;
   static REAL8 df, norm;
   static REAL4Vector signalvec, correlation;
   void   (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;
   static InspiralWaveOverlapIn overlapin;
   static InspiralWaveOverlapOut overlapout;
   static REAL8FrequencySeries shf;
   static RealFFTPlan *fwdp=NULL,*revp=NULL;
   static InspiralTemplate tmpltParam, param;
   static InspiralWaveNormaliseIn normin;


   quietFlag = 0;
   approx = PadeT1;
   tmplt = BCV;

   param.approximant = approx;
   param.massChoice = m1Andm2;
   param.ieta = 1;
   param.fLower = 40;
   param.order = 5;
   param.mass1 = 5.0;
   param.mass2 = 5.0;
   param.startTime=0.0;
   param.startPhase=0.0;
   param.nStartPad=0;
   param.fCutoff = 1000.;
   param.tSampling = 4096.;
   param.signalAmplitude = 1.0;
   param.nEndPad = 0;
   param.alpha = 0.L;

   tmpltParam = param;
   tmpltParam.massChoice = psi0Andpsi3;
   tmpltParam.approximant = tmplt;
   tmpltParam.alpha = 0.669L;
   tmpltParam.psi0 = 224810.L;
   tmpltParam.psi3 = -867.58L;
   tmpltParam.fFinal = 764.5L;
   param.alpha = 0.669L;
   param.psi0 = 224810.L;
   param.psi3 = -867.58L;
   param.fFinal = 764.5L;
   param.fFinal = tmpltParam.fFinal;


   i=1;
   while(i < argc)
   {
	   if (strcmp(argv[i],"-fl")==0)
	   {
		   param.fLower = atof(argv[++i]);
		   tmpltParam.fLower = param.fLower;
	   }
	   else if (strcmp(argv[i],"-order")==0)
	   {
		   param.order = atoi(argv[++i]);
		   tmpltParam.order = param.order;
	   }
	   else if (strcmp(argv[i],"-m1")==0)
		   param.mass1 = atof(argv[++i]);
	   else if (strcmp(argv[i],"-m2")==0)
		   param.mass2 = atof(argv[++i]);
	   else if (strcmp(argv[i],"-quiet")==0)
		   quietFlag = 1;
	   else if (strcmp(argv[i],"-alpha")==0)
		   tmpltParam.alpha = atof(argv[++i]);
	   else if (strcmp(argv[i],"-psi0")==0)
		   tmpltParam.psi0 = atof(argv[++i]);
	   else if (strcmp(argv[i],"-psi3")==0)
		   tmpltParam.psi3 = atof(argv[++i]);
	   else if (strcmp(argv[i],"-fcut")==0)
		   tmpltParam.fFinal = atof(argv[++i]);
	   else if (strcmp(argv[i], "-approximant")==0 || (strcmp(argv[i], "-signal")==0))
	   {
		   if (strcmp(argv[++i],"TaylorT1")==0)
			   approx = 0;
		   else if (strcmp(argv[i],"TaylorT2")==0)
			   approx = 1;
		   else if (strcmp(argv[i],"TaylorT3")==0)
			   approx = 2;
		   else if (strcmp(argv[i],"TaylorF1")==0)
			   approx = 3;
		   else if (strcmp(argv[i],"TaylorF2")==0)
			   approx = 4;
		   else if (strcmp(argv[i],"PadeT1")==0)
			   approx = 5;
		   else if (strcmp(argv[i],"PadeF1")==0)
			   approx = 6;
		   else if (strcmp(argv[i],"EOB")==0)
			   approx = 7;
		   else if (strcmp(argv[i],"BCV")==0)
			   approx = 8;
		   else if (strcmp(argv[i],"SpinTaylorT3")==0)
			   approx = 9;
		   param.approximant = approx;
	   }
	   else if (strcmp(argv[i], "-template")==0)
	   {
		   if (strcmp(argv[++i],"TaylorT1")==0)
			   tmplt = 0;
		   else if (strcmp(argv[i],"TaylorT2")==0)
			   tmplt = 1;
		   else if (strcmp(argv[i],"TaylorT3")==0)
			   tmplt = 2;
		   else if (strcmp(argv[i],"TaylorF1")==0)
			   tmplt = 3;
		   else if (strcmp(argv[i],"TaylorF2")==0)
			   tmplt = 4;
		   else if (strcmp(argv[i],"PadeT1")==0)
			   tmplt = 5;
		   else if (strcmp(argv[i],"PadeF1")==0)
			   tmplt = 6;
		   else if (strcmp(argv[i],"EOB")==0)
			   tmplt = 7;
		   else if (strcmp(argv[i],"BCV")==0)
			   tmplt = 8;
		   else if (strcmp(argv[i],"SpinTaylorT3")==0)
			   tmplt = 9;
		   tmpltParam.approximant = tmplt;
	   }
	   else
	   {
		   fprintf(stderr,"\nUSAGE: %s [options]\n", argv[0]);
		   fprintf(stderr,"The options are (with default values in brackets)\n");
		   fprintf(stderr,"      -quiet : if this flag is present, the output is restricted to the minimum\n");
		   fprintf(stderr,"      -alpha : BCV amplitude correction parameter (%7.2f)\n", param.alpha);
		   fprintf(stderr,"-approximant : Post-Newtonian model such as TaylorT1, PadeT1, EOB, BCV ...  (PadeT1)\n");
		   fprintf(stderr,"         -fl : lower frequency cutoff (%7.2f) Hz\n", param.fLower);
		   fprintf(stderr,"         -m1 : mass of primary (%7.2f) SolarMass\n", param.mass1);
		   fprintf(stderr,"         -m2 : mass of companion (%7.2f) SolarMass\n", param.mass2);
		   fprintf(stderr,"      -order : order of PN model (%7.2d)\n", param.order);
		   fprintf(stderr,"       -psi0 : Max value of psi0 (%7.2f)\n", param.psi0);
		   fprintf(stderr,"       -psi3 : Min value of psi  (%7.2f)\n", param.psi3);
		   fprintf(stderr,"       -fcut : Cutoff frequency for BCV (%7.2f)\n\n", param.fFinal);
		   return 1;

	   }
	   i++;
   }

/*---------------------------------------------------------------------------*/
/* User can choose allowed values of the various parameters below this line  */
/*---------------------------------------------------------------------------*/
   overlapin.nBegin = 0;
   overlapin.nEnd = 0;
   signalvec.length = 0.;
   param.approximant = EOB;
   LALInspiralWaveLength (&status, &signalvec.length, param);
   if (!quietFlag) fprintf(stdout, "signal length = %d\n", signalvec.length);

   correlation.length = signalvec.length;
   psdLength = signalvec.length/2 + 1;
   signalvec.data = (REAL4*) LALMalloc(sizeof(REAL4)*signalvec.length);
   correlation.data = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   memset( &(shf), 0, sizeof(REAL8FrequencySeries) );
   shf.f0 = 0;
   LALDCreateVector( &status, &(shf.data), psdLength );
   shf.deltaF = param.tSampling / signalvec.length;
   df = param.tSampling/(float) signalvec.length;
   LALNoiseSpectralDensity (&status, shf.data, noisemodel, df);

/*
 * Estimate the plans
 */
   LALCreateForwardRealFFTPlan(&status, &fwdp, signalvec.length, 0);
   LALCreateReverseRealFFTPlan(&status, &revp, signalvec.length, 0);
/* REPORTSTATUS(&status); */

   param.approximant = approx;
   if (param.approximant == BCV)
	   param.massChoice = psi0Andpsi3;
   else
	   param.massChoice = m1Andm2;

   LALInspiralParameterCalc (&status, &param);
   if (param.approximant != BCV && param.approximant != TaylorF1 && param.approximant != TaylorF2 )
   {
	   LALInspiralWave(&status, &correlation, &param);
	   LALREAL4VectorFFT(&status, &signalvec, &correlation, fwdp);
   }
   else
   {
	   LALInspiralWave(&status, &signalvec, &param);
   }

   normin.psd = shf.data;
   normin.df = param.tSampling / (REAL8) signalvec.length;
   normin.fCutoff = param.fFinal;
   normin.samplingRate = param.tSampling;
   LALInspiralWaveNormaliseLSO(&status, &signalvec, &norm, &normin);
   /*
   LALREAL4VectorFFT(&status, &correlation, &signalvec, revp);
   if (!quietFlag) for (i=0; i<correlation.length; i++) printf("%e\n", correlation.data[i]);
   */

   overlapin.param = tmpltParam;
   overlapin.psd = *shf.data;
   overlapin.signal = signalvec;
   overlapin.fwdp = fwdp;
   overlapin.revp = revp;

   overlapin.ifExtOutput = 0;
   LALInspiralWaveOverlap (&status,&correlation,&overlapout,&overlapin);
   if (!quietFlag) for (j=0; j<correlation.length; j++) printf("%e\n", correlation.data[j]);
   fprintf (stdout, "%e %e %e %e %e %e %e %e %e\n",
		   tmpltParam.psi0,
		   tmpltParam.psi3,
		   param.mass1,
		   param.mass2,
		   tmpltParam.totalMass,
		   param.totalMass,
		   overlapin.param.fFinal,
		   param.fFinal,
		   overlapout.max
		   );
   /* destroy the plans, correlation and signal */

   /*
      LALDDestroyVector( &status, &(shf.data) );
      if (signalvec.data != NULL) LALFree(signalvec.data);
      if (correlation.data != NULL) LALFree(correlation.data);
      LALDestroyRealFFTPlan(&status,&fwdp);
      LALDestroyRealFFTPlan(&status,&revp);
      LALCheckMemoryLeaks();
    */

   return(0);

}

void printf_timeseries (INT4 n, REAL4 *signal1, double delta, double t0)
{
  int i=0;

  do
     printf ("%e %e\n", i*delta+t0, *(signal1+i));
  while (n-++i);
  printf("&\n");
}
