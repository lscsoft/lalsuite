#include <lal/LALInspiralBank.h>
#include <lal/LALInspiral.h>

INT4 lalDebugLevel = 0;
void printf_timeseries (INT4, REAL4 *, REAL4, REAL4); 

int main () 
{
#if FIXME
   static LALStatus status;
   REAL4 randnum, dt, df, norm, overlap;
   INT4 n,i=10,pad;
   REAL4Vector signal;
   RandomInspiralSignalIn randIn;
   Detector choice;
   OverlapIn overlapin;
   OverlapOut overlapout;

   randIn.useed = 29184045;
   randIn.mMin = 3.0;
   randIn.MMax = 20.0;
   randIn.SignalAmp = 10.0;
   randIn.NoiseAmp = 1.0;
   randIn.type = 2;

   randIn.param.mass1=randIn.mMin;
   randIn.param.mass2=randIn.mMin;
   randIn.param.ieta=1; 
   randIn.param.startTime=0.0; 
   randIn.param.startPhase=0.0; 
   randIn.param.fLower=40.0; 
   randIn.param.fCutoff=1000.0;
   randIn.param.tSampling=4000.0;
   randIn.param.signalAmplitude=1.0;
   pad = randIn.param.nStartPad=2000;
   randIn.param.nEndPad=0;
   randIn.param.method=one;
   randIn.param.order=twoPointFivePN;
   randIn.param.domain=TimeDomain;
   randIn.param.approximant=pade;
   randIn.param.massChoice=m1Andm2;
   dt = 1./randIn.param.tSampling;

   LALInspiralWaveLength (&status, &signal.length, randIn.param);
   randIn.psd.length = signal.length/2 + 1;
   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*signal.length);
   randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length);
   choice = ligo;
   df = randIn.param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, &randIn.psd, choice, df);
   overlapin.psd = randIn.psd;
   while (i--) {
      randIn.param.nStartPad = (i+1) * pad;
      LALRandomInspiralSignal(&status, &signal, &randIn);
/*
      printf_timeseries (signal.length, signal.data, dt, 0.0);
*/
      overlapin.signal = signal;
      overlapin.param = randIn.param;
      LALInspiralWaveOverlap(&status, &overlap, &overlapout, &overlapin);
      fprintf(stderr,"%e, %e %e\n", randIn.param.mass1, randIn.param.mass2, overlap);
   }

/*
   n=16384;
   while (n--) {
      LALGaussianRandomNumber(&status, &randnum);
      printf("%e\n",randnum);
   }
*/
   return(1);
#else
   return 77;
#endif
}

void printf_timeseries (INT4 n, REAL4 *signal, REAL4 delta, REAL4 t_0) 
{
  int i=0;
  do 
     printf ("%e %e\n", i*delta+t_0, *(signal+i));

  while (n-++i); 
  printf("&\n");
}

