/*
*  Copyright (C) 2007 B.S. Sathyaprakash, Thomas Cokelaer
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
 * \author Cokelaer, T. and Sathyaprakash, B. S.
 * \file
 * \ingroup inspiral
 *
 * \brief Code for computing overalps of two waveforms.
 *
 * ### Usage ###
 *
 * \code
 * InspiralOverlap [options]
 *
 * The options are :
 * -quiet : if this flag is present, the output is minimal
 * -alpha : BCV amplitude correction parameter
 * -approximant : Post-Newtonian model of signal
 * -signal : same as -approximant
 * -template : Post-Newtonian model of template
 * -fl : lower frequency cutoff
 * -m1 : mass of primary
 * -m2 : mass of companion
 * -m3 : mass of primary
 * -m4 : mass of companion
 * -order : order of PN model
 * -psi0 : Max value of psi0
 * -psi3 : Min value of psi3
 * -fcut : Cutoff frequency for BCV
 * \endcode
 *
 * ### Description ###
 *
 * This code helps to compute the overlap of two waveforms with default
 * (or user specified) parameters.
 *
 * ### Exit codes ###
 *
 *
 * ### Uses ###
 *
 * This code directly uses the following functions (see those functions
 * to find out what they call in turn):
 * \code
 * LALInspiralWaveLength
 * LALInspiralCreateCoarseBank
 * LALInspiralParameterCalc
 * LALInspiralWave
 * LALInspiralWaveNormaliseLSO
 * LALInspiralWaveOverlap
 * LALNoiseSpectralDensity
 * LALCreateForwardRealFFTPlan
 * LALCreateReverseRealFFTPlan
 * LALForwardRealFFT
 * LALReverseRealFFT
 * LALDestroyRealFFTPlan
 * \endcode
 *
 * ### Notes ###
 *
 */


#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALNoiseModelsInspiral.h>
#include <lalapps.h>



void printf_timeseries (INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0);
void PrintParams(InspiralTemplate params, InspiralTemplate parans);


int
main (  int argc, char **argv )
{
   static LALStatus status;
   static INT4 i, approx, tmplt;
   static UINT4 psdLength, quietFlag = 0;
   static REAL8 df, norm;
   static REAL4Vector signalvec, correlation;
   void (*noisemodel)(LALStatus*, REAL8*, REAL8) = NULL;
   static InspiralWaveOverlapIn overlapin;
   static InspiralWaveOverlapOut overlapout;
   static REAL8FrequencySeries shf;
   static RealFFTPlan *fwdp=NULL,*revp=NULL;
   static InspiralTemplate tmpltParams, params;
   static InspiralWaveNormaliseIn normin;

   noisemodel = &LALVIRGOPsd;

   quietFlag = 0;

	/* ALL default Values here */
   approx = TaylorF2;
   tmplt = TaylorF2;

   params.approximant = approx;
   params.massChoice = m1Andm2;
   params.ieta = 1;
   params.fLower = 40;
   params.order = 5;
   params.mass1 = 5.0;
   params.mass2 = 5.0;
   params.startTime=0.0;
   params.startPhase=0.8819;
   params.nStartPad=1000;
   params.fCutoff = 1000.;
   params.tSampling = 2048.;
   params.signalAmplitude = 1.0;
   params.nEndPad = 0;
   params.alpha = 0.L;

   tmpltParams = params;

   tmpltParams.massChoice = psi0Andpsi3;
   tmpltParams.approximant = tmplt;
   tmpltParams.alpha = 0.669L;
   tmpltParams.psi0 = 224810.L;
   tmpltParams.psi3 = -867.58L;
   tmpltParams.fFinal = 764.5L;
   params.alpha = tmpltParams.alpha;
   params.psi0 = tmpltParams.psi0;
   params.psi3 = tmpltParams.psi3;
   params.fFinal = tmpltParams.fFinal;

   /* User parameters here */
   i=1;
   while(i < argc)
   {
	   if (strcmp(argv[i],"--fl1")==0)
		   params.fLower      = atof(argv[++i]);
	   else if (strcmp(argv[i],"--fl2")==0)
		   tmpltParams.fLower = atof(argv[++i]);
	   else if (strcmp(argv[i],"--order1")==0)
		   params.order = atoi(argv[++i]);
	   else if (strcmp(argv[i],"--order2")==0)
		   tmpltParams.order = atoi(argv[++i]);
	   else if (strcmp(argv[i],"--zeta2")==0)
	     params.Zeta2 = atof(argv[++i]);
	   else if  (strcmp(argv[i],"--sampling")==0)
	     {
	       params.tSampling = atof(argv[++i]);
	       tmpltParams.tSampling = params.tSampling;
	     }
	   else if (strcmp(argv[i],"--m1")==0)
		   params.mass1 = atof(argv[++i]);
	   else if (strcmp(argv[i],"--m2")==0)
		   params.mass2 = atof(argv[++i]);
	   else if (strcmp(argv[i],"--m3")==0)
		   tmpltParams.mass1 = atof(argv[++i]);
	   else if (strcmp(argv[i],"--m4")==0)
		   tmpltParams.mass2 = atof(argv[++i]);
	   else if (strcmp(argv[i],"--quiet")==0)
		   quietFlag = 1;
	   else if (strcmp(argv[i],"--alpha")==0)
		   tmpltParams.alpha = atof(argv[++i]);
	   else if (strcmp(argv[i],"--psi0")==0)
		   tmpltParams.psi0 = atof(argv[++i]);
	   else if (strcmp(argv[i],"--psi3")==0)
		   tmpltParams.psi3 = atof(argv[++i]);
	   else if (strcmp(argv[i],"--fcut")==0)
		   tmpltParams.fFinal = atof(argv[++i]);
	   else if ((strcmp(argv[i], "--signal")==0))
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
		   params.approximant = approx;
	   }
	   else if (strcmp(argv[i], "--template")==0)
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
		   tmpltParams.approximant = tmplt;
	   }
	   else
	   {
		   fprintf(stderr,"\nUSAGE: %s [options]\n", argv[0]);
		   fprintf(stderr,"The options are (with default values in brackets)\n");
		   fprintf(stderr,"      --quiet : if this flag is present, the output is restricted to the minimum\n");
		   fprintf(stderr,"      --alpha : BCV amplitude correction paramseter (%7.2f)\n", params.alpha);
		   fprintf(stderr,"     --signal : same as -approximant\n");
		   fprintf(stderr,"   --template : Post-Newtonian model of template (BCV)\n");
		   fprintf(stderr,"        --fl1 : lower frequency cutoff (%7.2f) Hz of signal\n", params.fLower);
		   fprintf(stderr,"        --fl2 : lower frequency cutoff (%7.2f) Hz of template\n", params.fLower);
		   fprintf(stderr,"      --zeta2 : zeta2 for EOB 3pn (%7.2f)\n", params.Zeta2);
		   fprintf(stderr,"         --m1 : mass of primary (%7.2f) SolarMass\n", params.mass1);
		   fprintf(stderr,"         --m2 : mass of companion (%7.2f) SolarMass\n", params.mass2);
		   fprintf(stderr,"         --m3 : mass of primary (%7.2f) in template\n", tmpltParams.mass1);
		   fprintf(stderr,"         --m4 : mass of companion (%7.2f) in template\n", tmpltParams.mass2);
		   fprintf(stderr,"     --order1 : order of PN model (%7.2d) of signal\n", params.order);
		   fprintf(stderr,"     --order2 : order of PN model (%7.2d) of template\n", params.order);
		   fprintf(stderr,"       --psi0 : Max value of psi0 (%7.2f)\n", params.psi0);
		   fprintf(stderr,"       --psi3 : Min value of psi3 (%7.2f)\n", params.psi3);
		   fprintf(stderr,"       --fcut : Cutoff frequency for BCV (%7.2f)\n\n", params.fFinal);
		   return 1;

	   }
	   i++;
   }

/*---------------------------------------------------------------------------*/
/* User can choose allowed values of the various parameters below this line  */
/*---------------------------------------------------------------------------*/
   overlapin.nBegin 	= 0;
   overlapin.nEnd 		= 0;
   signalvec.length 		= 0.;
   params.approximant 	= EOB;  /*just for allocating memory */
   LAL_CALL(LALInspiralWaveLength (&status, &signalvec.length, params),&status);
   if (!quietFlag) fprintf(stdout, "#signal length = %d\n", signalvec.length);


   correlation.length = signalvec.length;
   psdLength = signalvec.length/2 + 1;
   signalvec.data = (REAL4*) LALMalloc(sizeof(REAL4)*signalvec.length);
   correlation.data = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   memset( &(shf), 0, sizeof(REAL8FrequencySeries) );
   shf.f0 = 0;
   LAL_CALL(LALDCreateVector( &status, &(shf.data), psdLength ),&status);
   shf.deltaF = params.tSampling / signalvec.length;
   df = params.tSampling/(float) signalvec.length;
   LAL_CALL(LALNoiseSpectralDensity (&status, shf.data, noisemodel, df),&status);

/*
 * Estimate the plans
 */
   LAL_CALL(LALCreateForwardRealFFTPlan(&status, &fwdp, signalvec.length, 0),&status);
   LAL_CALL(LALCreateReverseRealFFTPlan(&status, &revp, signalvec.length, 0),&status);
/* REPORTSTATUS(&status); */


   params.approximant = approx;
   if (params.approximant == BCV)
	   params.massChoice = psi0Andpsi3;
   else
	   params.massChoice = m1Andm2;

   LAL_CALL(LALInspiralParameterCalc (&status, &params),&status);
   if (params.approximant != BCV && params.approximant != TaylorF1 && params.approximant != TaylorF2 )
   {
	   LAL_CALL(LALInspiralWave(&status, &correlation, &params),&status);
	   LAL_CALL(LALREAL4VectorFFT(&status, &signalvec, &correlation, fwdp),&status);
   }
   else
   {
	   LAL_CALL(LALInspiralWave(&status, &signalvec, &params),&status);
   }

   normin.psd = shf.data;
   normin.df = params.tSampling / (REAL8)signalvec.length;
   normin.fCutoff = params.fFinal;

  /* normin.fLower = 0;*/
/*      normin.fCutoff = params.tSampling / 2. -1;*/
   normin.samplingRate = params.tSampling;


   LAL_CALL(LALInspiralWaveNormaliseLSO(&status, &signalvec, &norm, &normin),&status);


   tmpltParams.approximant = tmplt;
   if (tmpltParams.approximant == BCV)
   {
	   tmpltParams.massChoice = psi0Andpsi3;
   }
   else
   {
	   tmpltParams.massChoice = m1Andm2;
   }
   LAL_CALL(LALInspiralParameterCalc(&status, &tmpltParams),&status);

   overlapin.param = tmpltParams;
   overlapin.psd = *shf.data;
   overlapin.signal = signalvec;
   overlapin.fwdp = fwdp;
   overlapin.revp = revp;
   for (i=0; (UINT4) i<correlation.length; i++) correlation.data[i] = 0.;

   LAL_CALL(LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin),&status);


   if (!quietFlag) for (i=0; (UINT4) i<correlation.length; i++)
     printf("%e\n", correlation.data[i]);
   if (tmplt == BCV)
	   fprintf (stdout, "## %d %d %e %e %e %e %e %e %e %e %e\n",
			   tmpltParams.approximant,
			   params.approximant,
			   tmpltParams.psi0,
			   tmpltParams.psi3,
			   params.mass1,
			   params.mass2,
			   tmpltParams.totalMass,
			   params.totalMass,
			   overlapin.param.fFinal,
			   params.fFinal,
			   overlapout.max
		   );
   else
	   fprintf (stdout, "## %d %d %e %e %e %e %e %e %e %e %e %e\n",
			   tmpltParams.approximant,
			   params.approximant,
			   tmpltParams.mass1,
			   params.mass1,
			   tmpltParams.mass2,
			   params.mass2,
			   tmpltParams.totalMass,
			   params.totalMass,
			   overlapin.param.fFinal,
			   params.fFinal,tmpltParams.fLower,
			   overlapout.max
		   );
   /* destroy the plans, correlation and signal */

      LAL_CALL(LALDDestroyVector( &status, &(shf.data) ),&status);
      if (signalvec.data != NULL) LAL_CALL(LALFree(signalvec.data),&status);
      if (correlation.data != NULL) LAL_CALL(LALFree(correlation.data),&status);
      LALDestroyRealFFTPlan(&status,&fwdp);
      LALDestroyRealFFTPlan(&status,&revp);
      LALCheckMemoryLeaks();
      if (!quietFlag)
      	PrintParams(params, tmpltParams);
   return(0);

}

void printf_timeseries (INT4 n, REAL4 *signalvec, double delta, double t0)
{
  int i=0;

  do
     printf ("%e %e\n", i*delta+t0, *(signalvec+i));
  while (n-++i);
  printf("&\n");
}




void
PrintParams(InspiralTemplate params1, InspiralTemplate params2)
{
  printf("\n#ieta     = %7.2d %7.2d\n",params1.ieta,params2.ieta);
  printf("#level      = %7.2d %7.2d\n",params1.level,params2.level);
  printf("#nStartPad  = %7.2d %7.2d\n",params1.nStartPad,params2.nStartPad);
  printf("#nEndPad    = %7.2d %7.2d\n",params1.nEndPad,params2.nEndPad);
  printf("#minMatch   = %7.2f %7.2f\n",params1.minMatch,params2.minMatch);
  printf("#mass1      = %7.2f %7.2f\n",params1.mass1,params2.mass1);
  printf("#mass2      = %7.2f %7.2f\n",params1.mass2,params2.mass2);
  printf("#totalMass  = %7.2f %7.2f\n",params1.totalMass,params2.totalMass);
  printf("#chirpmass  = %7.2f %7.2f\n",params1.chirpMass,params2.chirpMass);
  printf("#psi0       = %7.2f %7.2f\n",params1.psi0,params2.psi0);
  printf("#psi3       = %7.2f %7.2f\n",params1.psi3,params2.psi3);
  printf("#fendBCV    = %7.2f %7.2f\n",params1.fFinal,params2.fFinal);
  printf("#alpha      = %7.2f %7.2f\n",params1.alpha,params2.alpha);
  printf("#alpha1     = %7.2f %7.2f\n",params1.alpha1,params2.alpha1);
  printf("#alpha2     = %7.2f %7.2f\n",params1.alpha2,params2.alpha1);
  printf("#beta       = %7.2f %7.2f\n",params1.beta,params2.beta);
  printf("#tc         = %7.2f %7.2f\n",params1.tC,params2.tC);
  printf("#eta        = %7.2f %7.2f\n",params1.eta,params2.eta);
  printf("#fLower     = %7.2f %7.2f\n",params1.fLower,params2.fLower);
  printf("#fcutoff    = %7.2f %7.2f\n",params1.fCutoff,params2.fCutoff);
  printf("#tsampling  = %7.2f %7.2f\n",params1.tSampling,params2.tSampling);
  printf("#phase0     = %7.2f %7.2f\n",params1.startPhase,params2.startPhase);
  printf("#t0         = %7.2f %7.2f\n",params1.startTime,params2.startTime);
  printf("#fFinal     = %7.2f %7.2f\n",params1.fFinal,params2.fFinal);
  printf("#zeta2      = %7.2f %7.2f\n",params1.Zeta2,params2.Zeta2);
  printf("#order      = %7.2d %7.2d\n",params1.order,params2.order);
  printf("#approximant= %7.2d %7.2d\n",params1.approximant,params2.approximant);
}
