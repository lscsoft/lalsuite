#include "LALInference.h"
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

/* Bad idea if a,b are anything but variable references---i.e. MIN_INT(a++, b) will *break*. */
#define MIN_INT(a, b) ((a) < (b) ? (a) : (b))
#define MAX_INT(a, b) ((a) > (b) ? (a) : (b))

int main(void) {
  REAL8TimeSeries *conv, *data, *response, *exactConv;
  const UINT4 N = 101;
  const LIGOTimeGPS zero = {0, 0};
  UINT4 i;

  conv = XLALCreateREAL8TimeSeries("FFT Convolution", &zero, 0.0, 0.01, &lalDimensionlessUnit, N);
  data = XLALCreateREAL8TimeSeries("data", &zero, 0.0, 0.01, &lalDimensionlessUnit, N);
  response = XLALCreateREAL8TimeSeries("response function", &zero, 0.0, 0.01, &lalDimensionlessUnit, N);
  exactConv = XLALCreateREAL8TimeSeries("exact convolution", &zero, 0.0, 0.01, &lalDimensionlessUnit, N);

  data->data->data[0] = 1.0;
  for (i = 1; i < N; i++) {
    /* Data is 1/i */
    data->data->data[i] = 1.0 / i;
  }

  response->data->data[0] = 1.0;
  for (i = 1; i <= (N-1)/2; i++) {
    response->data->data[i] = exp(-(i/10.0)); /* Decaying exponential with time constant 10. */
    response->data->data[N-i] = exp(-(i/10.0)); /* Same in the negative, wrapped dimension. */
  }
  if (N % 2 == 0) {
    response->data->data[N/2] = exp(-(N/20.0));
  }

  /* Exact, O(N^2) convolution */
  for (i = 0; i < N; i++) {
    REAL8 sum = 0.0;
    UINT4 j;

    for (j = MAX_INT(0, ((int)i)-((int) N/2)); j <= MIN_INT(i+N/2, N-1); j++) {
      UINT4 rIndex = (j > i ? i + (N - j) : i - j);

      sum += data->data->data[j]*response->data->data[rIndex];
    }

    exactConv->data->data[i] = sum;
  }

  convolveTimeSeries(conv, data, response);

  for (i = 0; i < N; i++) {
    if (fabs(exactConv->data->data[i] - conv->data->data[i]) > 1e-8) {
      fprintf(stderr, "Index %d differs (exact = %g, FFT = %g)\n", 
              i, exactConv->data->data[i], conv->data->data[i]);
    }
  }

  XLALDestroyREAL8TimeSeries(conv);
  XLALDestroyREAL8TimeSeries(data);
  XLALDestroyREAL8TimeSeries(response);
  XLALDestroyREAL8TimeSeries(exactConv);

  return 0;
}
