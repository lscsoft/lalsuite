#include "LALInference.h"
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Random.h>

int main(void) {
  const UINT4 N = 8192;
  const UINT4 NPSD = N/2+1;
  const REAL8 T = 8.38;
  const REAL8 dT = T/(N-1);
  const REAL8 fNy = 1.0/(2.0*dT);
  const REAL8 dF = fNy/(NPSD-1);

  const REAL8 Af0 = 15.0;
  const REAL8 Afdot0 = 10.0;

  const REAL8 fMin = 5;
  const REAL8 fMax = 99;

  UINT4 i;

  LIGOTimeGPS zero = {0,0};

  REAL8FrequencySeries *PSD;
  COMPLEX16FrequencySeries *Af, *noiseF;
  REAL8TimeSeries *TDW, *A, *noiseT;

  REAL8FFTPlan *fwd, *rev;

  REAL8 sum = 0.0, overlap;

  RandomParams *params;

  params = XLALCreateRandomParams(0); /* Seeded off the current time. */
  PSD=XLALCreateREAL8FrequencySeries("PSD", &zero, 0.0, dF, &lalDimensionlessUnit, NPSD);
  Af=XLALCreateCOMPLEX16FrequencySeries("Af", &zero, 0.0, dF, &lalDimensionlessUnit, NPSD);
  noiseF=XLALCreateCOMPLEX16FrequencySeries("noise", &zero, 0.0, dF, &lalDimensionlessUnit, NPSD);

  TDW=XLALCreateREAL8TimeSeries("TDW", &zero, 0.0, dT, &lalDimensionlessUnit, N);
  A=XLALCreateREAL8TimeSeries("A", &zero, 0.0, dT, &lalDimensionlessUnit, N);
  noiseT=XLALCreateREAL8TimeSeries("noiseT", &zero, 0.0, dT, &lalDimensionlessUnit, N);

  fwd = XLALCreateForwardREAL8FFTPlan(N, 1);
  rev = XLALCreateReverseREAL8FFTPlan(N, 1);

  for (i = 0; i < N; i++) {
    REAL8 t = i*dT;
    const REAL8 sigAmp = 0.0;

    A->data->data[i] = sigAmp*t*cos(2.0*M_PI*(Af0 + Afdot0*t)*t);
  }

  PSD->data->data[0] = 1.0;
  for (i = 1; i < NPSD; i++) {
    REAL8 f = i*dF;

    PSD->data->data[i] = fabs(1.0/f); 

    /* Make some noise. */
    noiseF->data->data[i].re = 0.5*sqrt(PSD->data->data[i]/dF)*XLALNormalDeviate(params);
    noiseF->data->data[i].im = 0.5*sqrt(PSD->data->data[i]/dF)*XLALNormalDeviate(params);
  }

  XLALREAL8FreqTimeFFT(noiseT, noiseF, rev);

  for (i = 0; i < N; i++) {
    A->data->data[i] += noiseT->data->data[i];
  }

  PSDToTDW(TDW, PSD, rev, fMin, fMax);

  XLALREAL8TimeFreqFFT(Af, A, fwd);

  overlap = timeDomainOverlap(TDW, A, A);

  sum = 0.0;
  for (i = 0; i < NPSD; i++) {
    REAL8 re = Af->data->data[i].re;
    REAL8 im = Af->data->data[i].im;
    REAL8 dFLocal = Af->deltaF;
    REAL8 f = dFLocal*i;

    if (fMin <= f && f <= fMax) {
      sum += 4.0*dFLocal*(re*re + im*im)/PSD->data->data[i];
    } else {
      /* Do Nothing. */
    }
  }

  fprintf(stderr, "Time domian overlap = %g, freq domain sum = %g, ratio = %g\n",
          overlap, sum, fabs(overlap/sum));

  return 0;
}
