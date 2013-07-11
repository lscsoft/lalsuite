
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_spline.h>

#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Interpolate.h>

int main(void) {
	size_t len=100;
	double dt=0.02;
	double drift=0.0;
	REAL8Sequence *sample_time = XLALCreateREAL8Sequence(len);
	size_t i=0;
	sample_time->data[0] = dt*i;
	for(i=1; i<len; i++) {
		sample_time->data[i] = dt*i+drift;
		drift += 1e-3;
	}

	double frequency = 20.0 / 2 / 3.14159;

	FILE* fref = fopen( "fref.txt", "w" );
	REAL8Sequence *fcn = XLALCreateREAL8Sequence(len);
	for(i=0; i<len; i++) {
		fcn->data[i] = sin( sample_time->data[i] * frequency );
		fprintf( fref, "%f %f\n", sample_time->data[i], fcn->data[i] );
	}
	fclose(fref);

	LIGOTimeGPS ep = {0, 0};
	
	REAL8TimeSeries *ts = XLALCreateREAL8TimeSeries("intrp test", &ep, 0.0, dt*0.9, &lalDimensionlessUnit, len);

	int ret = XLALREAL8TimeSeriesInterpolation(ts, fcn, sample_time, NULL, len, gsl_interp_cspline);

	REAL8 start = XLALGPSGetREAL8(&ts->epoch);
	REAL8 tolerance = 1e-6;
	FILE* fout = fopen( "fout.txt", "w" );
	for(i=0; i<ts->data->length; i++) {
		REAL8 t = ts->deltaT * i + start;
		REAL8 fcnval = sin(frequency*t);
		REAL8 diff = fabs(ts->data->data[i] - fcnval);
		fprintf( fout, "%f %f %f\n", t, fcnval, ts->data->data[i] );
		if (diff > tolerance) {
			fprintf(stderr, "%f %g\n", t, diff);
			return -1;
		}
	}
	fclose(fout);

	return ret;
}
