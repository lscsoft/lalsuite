#include "LALInference.h"
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>

int main() {
  const UINT4 N = 256;
  const UINT4 NPSD = N/2+1;
  const REAL8 T = 1.0;
  const REAL8 dT = T/(N-1);
  const REAL8 fNy = 1.0/(2.0*dT);
  const REAL8 dF = fNy/(NPSD-1);

  const REAL8 Af0 = 15.0;
  const REAL8 Afdot0 = 10.0;

  UINT4 i;

  LIGOTimeGPS zero = {0,0};

  REAL8FrequencySeries *PSD;
  COMPLEX16FrequencySeries *Af;
  REAL8TimeSeries *TDW, *A;

  REAL8FFTPlan *fwd, *rev;

  REAL8 sum = 0.0, overlap;

  PSD=XLALCreateREAL8FrequencySeries("PSD", &zero, 0.0, dF, &lalDimensionlessUnit, NPSD);
  Af=XLALCreateCOMPLEX16FrequencySeries("Af", &zero, 0.0, dF, &lalDimensionlessUnit, NPSD);

  TDW=XLALCreateREAL8TimeSeries("TDW", &zero, 0.0, dT, &lalDimensionlessUnit, N);
  A=XLALCreateREAL8TimeSeries("A", &zero, 0.0, dT, &lalDimensionlessUnit, N);

  fwd = XLALCreateForwardREAL8FFTPlan(N, 1);
  rev = XLALCreateReverseREAL8FFTPlan(N, 1);

  for (i = 0; i < N; i++) {
    REAL8 t = i*dT;

    A->data->data[i] = t*cos(2.0*M_PI*(Af0 + Afdot0*t)*t);
  }

  PSD->data->data[0] = 1.0;
  for (i = 1; i < NPSD; i++) {
    REAL8 f = i*dF;

    PSD->data->data[i] = 1.0/f; /* 1/f noise. */
  }

  PSDToTDW(TDW, PSD, rev);

  XLALREAL8TimeFreqFFT(Af, A, fwd);

  overlap = timeDomainOverlap(TDW, A, A);

  for (i = 0; i < NPSD; i++) {
    REAL8 re = Af->data->data[i].re;
    REAL8 im = Af->data->data[i].im;

    sum += 2.0*(re*re + im*im)/PSD->data->data[i];
  }

  fprintf(stderr, "Time domian overlap = %g, freq domain sum = %g, ratio = %g\n",
          overlap, sum, fabs(overlap/sum));

  return 0;
}
