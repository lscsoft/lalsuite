#include "LALInference.h"
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/RealFFT.h>

int main(void) {
  REAL8FrequencySeries *PSD;
  REAL8TimeSeries *TDW;
  REAL8FFTPlan *plan;
  UINT4 i;
  const UINT4 NPSD = 10000;
  const REAL8 PSDFMax = 100;
  const REAL8 PSDDF = PSDFMax / (NPSD - 1);
  LIGOTimeGPS epochZero = {0, 0};

  PSD = XLALCreateREAL8FrequencySeries("Flat PSD", &epochZero, 0.0, PSDDF, &lalDimensionlessUnit, NPSD);
  TDW = XLALCreateREAL8TimeSeries("TDW for flat PSD", &epochZero, 0.0, 1.0, &lalDimensionlessUnit, 2*(NPSD - 1));

  plan = XLALCreateReverseREAL8FFTPlan(TDW->data->length, 1); /* Let's measure carefully, even though we're only doing one of these. */

  for (i = 0; i < NPSD; i++) {
    /* 1/f noise. */
    PSD->data->data[i] = 2.0/(i*PSDDF); /* Flat two-sided PSD between -100 and 100 would be 1.0 */
  }

  PSDToTDW(TDW, PSD, plan, 0.0, 1.0/0.0);

  for (i = 0; i < TDW->data->length; i++) {
    double t = XLALGPSGetREAL8(&TDW->epoch) + i*TDW->deltaT;
    printf("%g %g\n", t, TDW->data->data[i]);
  }

  XLALDestroyREAL8FFTPlan(plan);
  XLALDestroyREAL8FrequencySeries(PSD);
  XLALDestroyREAL8TimeSeries(TDW);

  return 0;
}
