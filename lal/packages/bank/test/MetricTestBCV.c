#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/AVFactories.h>

INT4 lalDebugLevel=1;

int 
main ( void )
{  
  InspiralMetric metric;
  LALStatus            *status = LALCalloc(1, sizeof(*status));
  InspiralTemplate     *params;
  REAL8FrequencySeries    psd;
  void *(noisemodel) = LALLIGOIPsd;
  UINT4 numPSDpts = 262144;    
  REAL8 tSampling;
  REAL8 mismatch;

  mismatch = 0.03;
  params = (InspiralTemplate *)LALMalloc(sizeof(InspiralTemplate));

  params->alpha = 0.L;
  params->fLower = 30;    
  params->fCutoff = 400;

  tSampling = 4096.L;

  memset( &(psd), 0, sizeof(REAL8FrequencySeries) );
  psd.f0 = 0;
  LALDCreateVector( status, &(psd.data), numPSDpts );
  psd.deltaF = tSampling / (2.L*(REAL8) psd.data->length + 1.L);
  LALNoiseSpectralDensity (status, psd.data, noisemodel, psd.deltaF );
   
  LALInspiralComputeMetricBCV(status->statusPtr, &metric, &psd, params);

  fprintf(stderr, "%e %e %e\n", metric.G00, metric.G01, metric.G11);
  fprintf(stderr, "%e %e %e\n", metric.g00, metric.g11, metric.theta);
  fprintf(stderr, "dp0=%e dp1=%e\n", sqrt (mismatch/metric.G00), sqrt (mismatch/metric.G11));
  fprintf(stderr, "dP0=%e dP1=%e\n", sqrt (mismatch/metric.g00), sqrt (mismatch/metric.g11));

  
  {
  double MM;
  double dp0, dp1;
  long n=100;
  double dp0min=-5750;
  double dp0max=5750;
  double dp1min=-220;
  double dp1max=220;
  double d0=(dp0max-dp0min)/(double)n;
  double d1=(dp1max-dp1min)/(double)n;
  for ( dp0= dp0min; dp0<=dp0max ; dp0+=d0)
    {
      for ( dp1= dp1min; dp1<=dp1max ; dp1+=d1)
	{  
	  
	  MM = 1. - (metric.G00 * dp0 * dp0
		  +  metric.G01 * dp0 * dp1
		  +  metric.G01 * dp1 * dp0
		  +  metric.G11 * dp1 * dp1);
	  printf("%f %f %f\n", dp0, dp1, MM);	  
	}
              printf("\n");
    }
  }



  return 0;

}
