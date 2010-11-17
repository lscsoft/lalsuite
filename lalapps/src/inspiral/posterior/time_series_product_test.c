#include "LALInference.h"
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Date.h>

int main() {
  REAL8TimeSeries *xSamp, *x2Samp;
  LIGOTimeGPS xStart, x2Start;
  REAL8 dtx = 0.05, dtx2 = 0.0033;
  REAL8 T0x = 0.5, T0x2 = 0.333;
  REAL8 Tx = 1.0, Tx2 = 1.0;
  size_t Nx = (size_t)(round(Tx/dtx))+1, Nx2 = (size_t)(round(Tx2/dtx2))+1; /* Approx one-second duration for both. */
  UINT4 i;
  REAL8 lowerLimit, upperLimit, trueIntegral, approxIntegral;

  /* Re-calibrate the time lengths due to rounding in Nx, Nx2. */
  Tx = dtx*(Nx-1);
  Tx2 = dtx2*(Nx2-1);

  lowerLimit = T0x;
  upperLimit = T0x2 + Tx2;
  
  trueIntegral = (1.0/4.0)*(upperLimit*upperLimit*upperLimit*upperLimit - lowerLimit*lowerLimit*lowerLimit*lowerLimit); /* 1/4 (thigh^4 - tlow^4) */

  XLALGPSSet(&xStart, 0, (INT4)(round(T0x*1e9)) ); /* 0.5 seconds */
  XLALGPSSet(&x2Start, 0, (INT4)(round(T0x2*1e9))); /* 0.333 seconds */

  xSamp = XLALCreateREAL8TimeSeries("f(x) = x", &xStart, 0.0, dtx, &lalDimensionlessUnit, Nx);
  x2Samp = XLALCreateREAL8TimeSeries("f(x) = x^2", &x2Start, 0.0, dtx2, &lalDimensionlessUnit, Nx2);

  for (i = 0; i < Nx; i++) {
    xSamp->data->data[i] = T0x + dtx*i; /* f(x) = x */
  }

  for (i = 0; i < Nx2; i++) {
    double t = T0x2 + dtx2*i;
    x2Samp->data->data[i] = t*t; /* f(x) = x^2 */
  }

  approxIntegral = integrateSeriesProduct(xSamp, x2Samp);

  printf("Exact integral = %g, approx integral = %g, delta = %g\n", trueIntegral, approxIntegral, fabs(approxIntegral-trueIntegral));
  
  return 0;
}
