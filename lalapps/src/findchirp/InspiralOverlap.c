/******************************** <lalVerbatim file="InspiralOverlap">
Author: Cokelaer, T. and Sathyaprakash, B. S.
< $Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{InspiralOverlap.c}}
\label{ss:InspiralOverlap.c}

Code for computing overalps of two waveforms.

\subsubsection*{Usage}
\begin{verbatim}
InspiralOverlap [options]

The options are :
      -quiet : if this flag is present, the output is minimal
      -alpha : BCV amplitude correction parameter 
-approximant : Post-Newtonian model of signal 
     -signal : same as -approximant
   -template : Post-Newtonian model of template 
         -fl : lower frequency cutoff 
         -m1 : mass of primary 
         -m2 : mass of companion 
         -m3 : mass of primary 
         -m4 : mass of companion
      -order : order of PN model 
       -psi0 : Max value of psi0 
       -psi3 : Min value of psi3 
       -fcut : Cutoff frequency for BCV 
\end{verbatim}

\subsubsection*{Description}

This code helps to compute the overlap of two waveforms with default
(or user specified) parameters.

\subsubsection*{Exit codes}
\input{InspiralOverlapCE}

\subsubsection*{Uses}
This code directly uses the following functions (see those functions
to find out what they call in turn):
\begin{verbatim}
LALInspiralWaveLength
LALInspiralCreateCoarseBank
LALInspiralParameterCalc
LALInspiralWave
LALInspiralWaveNormaliseLSO
LALInspiralWaveOverlap
LALNoiseSpectralDensity 
LALCreateForwardRealFFTPlan
LALCreateReverseRealFFTPlan
LALForwardRealFFT
LALReverseRealFFT
LALDestroyRealFFTPlan
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{InspiralOverlapCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="InspiralOverlapCE"> */
/***************************** </lalErrTable> */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
 
void printf_timeseries (INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0);
void PrintParams(InspiralTemplate params, InspiralTemplate parans);
INT4 lalDebugLevel=1;
     
int 
main (  int argc, char **argv ) 
{
   static LALStatus status;
   static INT4 i, approx, tmplt;
   static UINT4 psdLength, quietFlag = 0, optFlag=0;
   static REAL8 df, norm;
   static REAL4Vector signal, correlation;
   void   *noisemodel = LALLIGOIPsd;
/*   void   *noisemodel = LALVIRGOPsd;*/
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
   param.startPhase=0.88189; 
   param.nStartPad=1000;
   param.fCutoff = 1000.;
   param.tSampling = 2048.;
   param.signalAmplitude = 1.0;
   param.nEndPad = 0;
   param.alpha = 0.L;
   
   tmpltParam = param;

   tmpltParam.massChoice = psi0Andpsi3;
   tmpltParam.approximant = tmplt;
   tmpltParam.alpha = 0.669L;
   tmpltParam.psi0 = 224810.L;
   tmpltParam.psi3 = -867.58L;
   tmpltParam.fendBCV = 764.5L;
   param.alpha = tmpltParam.alpha;
   param.psi0 = tmpltParam.psi0;
   param.psi3 = tmpltParam.psi3;
   param.fendBCV = tmpltParam.fendBCV;


   i=1;
   while(i < argc)
   {
	   if (strcmp(argv[i],"-fl1")==0)
		   param.fLower      = atof(argv[++i]);
	   else if (strcmp(argv[i],"-fl2")==0)
		   tmpltParam.fLower = atof(argv[++i]);
	   else if (strcmp(argv[i],"-order1")==0)
		   param.order = atoi(argv[++i]); 
	   else if (strcmp(argv[i],"-order2")==0)
		   tmpltParam.order = atoi(argv[++i]);
	   else if (strcmp(argv[i],"-zeta2")==0)
	     param.Zeta2 = atof(argv[++i]);
	   else if  (strcmp(argv[i],"-sampling")==0)
	     {
	       param.tSampling = atof(argv[++i]);
	       tmpltParam.tSampling = param.tSampling;
	     }
	   else if (strcmp(argv[i],"-m1")==0)
		   param.mass1 = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-m2")==0)
		   param.mass2 = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-m3")==0)
		   tmpltParam.mass1 = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-m4")==0)
		   tmpltParam.mass2 = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-quiet")==0)
		   quietFlag = 1;
	   else if (strcmp(argv[i],"-opt")==0)
		   optFlag = 1;
	   else if (strcmp(argv[i],"-alpha")==0)
		   tmpltParam.alpha = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-psi0")==0)
		   tmpltParam.psi0 = atof(argv[++i]);
	   else if (strcmp(argv[i],"-psi3")==0)
		   tmpltParam.psi3 = atof(argv[++i]);
	   else if (strcmp(argv[i],"-fcut")==0)
		   tmpltParam.fendBCV = atof(argv[++i]);
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
		   fprintf(stderr,"-approximant : Post-Newtonian model of signal (PadeT1)\n");
		   fprintf(stderr,"     -signal : same as -approximant\n");
		   fprintf(stderr,"   -template : Post-Newtonian model of template (BCV)\n");
		   fprintf(stderr,"        -fl1 : lower frequency cutoff (%7.2f) Hz of signal\n", param.fLower);
		   fprintf(stderr,"        -fl2 : lower frequency cutoff (%7.2f) Hz of template\n", param.fLower);
		   fprintf(stderr,"      -Zeta2 : zeta2 for EOB 3pn (%7.2f)\n", param.Zeta2);
		   fprintf(stderr,"         -m1 : mass of primary (%7.2f) SolarMass\n", param.mass1);
		   fprintf(stderr,"         -m2 : mass of companion (%7.2f) SolarMass\n", param.mass2);
		   fprintf(stderr,"         -m3 : mass of primary (%7.2f) in template\n", tmpltParam.mass1);
		   fprintf(stderr,"         -m4 : mass of companion (%7.2f) in template\n", tmpltParam.mass2);
		   fprintf(stderr,"      -order : order of PN model (%7.2d) of signal\n", param.order);
		   fprintf(stderr,"      -order : order of PN model (%7.2d) of template\n", param.order);
		   fprintf(stderr,"       -psi0 : Max value of psi0 (%7.2f)\n", param.psi0);
		   fprintf(stderr,"       -psi3 : Min value of psi3 (%7.2f)\n", param.psi3);
		   fprintf(stderr,"       -fcut : Cutoff frequency for BCV (%7.2f)\n\n", param.fendBCV);
		   return 1;	

	   }
	   i++;       
   }

/*---------------------------------------------------------------------------*/
/* User can choose allowed values of the various parameters below this line  */
/*---------------------------------------------------------------------------*/
   overlapin.nBegin = 0;
   overlapin.nEnd = 0;
   signal.length = 0.;
   param.approximant = EOB;
   LALInspiralWaveLength (&status, &signal.length, param);
   if (!quietFlag) fprintf(stdout, "signal length = %d\n", signal.length);

   correlation.length = signal.length;
   psdLength = signal.length/2 + 1;
   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*signal.length);
   correlation.data = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   memset( &(shf), 0, sizeof(REAL8FrequencySeries) );
   shf.f0 = 0;
   LALDCreateVector( &status, &(shf.data), psdLength );
   shf.deltaF = param.tSampling / signal.length;
   df = param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, shf.data, noisemodel, df);

/* 
 * Estimate the plans 
 */
   LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0);
   LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0);
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
	   LALREAL4VectorFFT(&status, &signal, &correlation, fwdp);
   }
   else
   {
	   LALInspiralWave(&status, &signal, &param);
   }

   normin.psd = shf.data;
   normin.df = param.tSampling / (REAL8)signal.length;
   normin.fCutoff = param.fFinal;
  /* normin.fLower = 0;*/
/*      normin.fCutoff = param.tSampling / 2. -1;*/
   normin.samplingRate = param.tSampling;


   LALInspiralWaveNormaliseLSO(&status, &signal, &norm, &normin);
   

   tmpltParam.approximant = tmplt;
   if (tmpltParam.approximant == BCV)
   {
	   tmpltParam.massChoice = psi0Andpsi3;
   }
   else
   {
	   tmpltParam.massChoice = m1Andm2;
   }
   LALInspiralParameterCalc (&status, &tmpltParam);

   overlapin.param = tmpltParam;
   overlapin.psd = *shf.data;
   overlapin.signal = signal;
   overlapin.fwdp = fwdp;
   overlapin.revp = revp;
   for (i=0; (UINT4) i<correlation.length; i++) correlation.data[i] = 0.;
   LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);


   if (!quietFlag) for (i=0; (UINT4) i<correlation.length; i++) printf("%e\n", correlation.data[i]);
   if (tmplt == BCV)
	   fprintf (stdout, "%d %d %e %e %e %e %e %e %e %e %e\n", 
			   tmpltParam.approximant,
			   param.approximant,
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
   else
	   fprintf (stdout, "%d %d %e %e %e %e %e %e %e %e %e\n", 
			   tmpltParam.approximant,
			   param.approximant,
			   tmpltParam.mass1, 
			   param.mass1,
			   tmpltParam.mass2, 
			   param.mass2,
			   tmpltParam.totalMass, 
			   param.totalMass,
			   overlapin.param.fFinal, 
			   param.fFinal,
			   overlapout.max
		   );
   /* destroy the plans, correlation and signal */

      LALDDestroyVector( &status, &(shf.data) );
      if (signal.data != NULL) LALFree(signal.data);
      if (correlation.data != NULL) LALFree(correlation.data);
      LALDestroyRealFFTPlan(&status,&fwdp);   
      LALDestroyRealFFTPlan(&status,&revp);
      LALCheckMemoryLeaks();
      if (!quietFlag) 
      	PrintParams(param, tmpltParam);
   return(0);

}

void printf_timeseries (INT4 n, REAL4 *signal, double delta, double t0) 
{
  int i=0;

  do 
     printf ("%e %e\n", i*delta+t0, *(signal+i));
  while (n-++i); 
  printf("&\n");
}




void
PrintParams(InspiralTemplate params1, InspiralTemplate param2)
{
  printf("\n#ieta     = %7.2d %7.2d\n",params1.ieta,param2.ieta);
  printf("#level      = %7.2d %7.2d\n",params1.level,param2.level);
  printf("#nStartPad  = %7.2d %7.2d\n",params1.nStartPad,param2.nStartPad);
  printf("#nEndPad    = %7.2d %7.2d\n",params1.nEndPad,param2.nEndPad);
  printf("#minMatch   = %7.2f %7.2f\n",params1.minMatch,param2.minMatch);
  printf("#mass1      = %7.2f %7.2f\n",params1.mass1,param2.mass1); 
  printf("#mass2      = %7.2f %7.2f\n",params1.mass2,param2.mass2);
  printf("#totalMass  = %7.2f %7.2f\n",params1.totalMass,param2.totalMass); 
  printf("#chirpmass  = %7.2f %7.2f\n",params1.chirpMass,param2.chirpMass); 
  printf("#psi0       = %7.2f %7.2f\n",params1.psi0,param2.psi0);
  printf("#psi3       = %7.2f %7.2f\n",params1.psi3,param2.psi3);
  printf("#fendBCV    = %7.2f %7.2f\n",params1.fendBCV,param2.fendBCV);
  printf("#alpha      = %7.2f %7.2f\n",params1.alpha,param2.alpha);
  printf("#alpha1     = %7.2f %7.2f\n",params1.alpha1,param2.alpha1);
  printf("#alpha2     = %7.2f %7.2f\n",params1.alpha2,param2.alpha1);
  printf("#beta       = %7.2f %7.2f\n",params1.beta,param2.beta);
  printf("#tc         = %7.2f %7.2f\n",params1.tC,param2.tC); 
  printf("#eta        = %7.2f %7.2f\n",params1.eta,param2.eta);
  printf("#fLower     = %7.2f %7.2f\n",params1.fLower,param2.fLower);
  printf("#fcutoff    = %7.2f %7.2f\n",params1.fCutoff,param2.fCutoff);
  printf("#tsampling  = %7.2f %7.2f\n",params1.tSampling,param2.tSampling);
  printf("#phase0     = %7.2f %7.2f\n",params1.startPhase,param2.startPhase);
  printf("#t0         = %7.2f %7.2f\n",params1.startTime,param2.startTime);
  printf("#fFinal     = %7.2f %7.2f\n",params1.fFinal,param2.fFinal);
  printf("#zeta2      = %7.2f %7.2f\n",params1.Zeta2,param2.Zeta2);
  printf("#order      = %7.2d %7.2d\n",params1.order,param2.order);
  printf("#approximant= %7.2d %7.2d\n",params1.approximant,param2.approximant);
}
