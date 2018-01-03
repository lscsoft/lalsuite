#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/Date.h>
#include <lal/LALSimulation.h>
#include <lal/LALSimBurst.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#define DELTA_T		(1.0 / 16384)	/* seconds */
#define SIMLENGTH	(16384 * 16)	/* samples */
#define DSTLENGTH	(16384 * 128)	/* samples */
#define OFFSET		86.332874431	/* seconds */
#define REAL4THRESH	.5e-6
#define REAL8THRESH	1e-12


static int TestXLALSimAddInjectionREAL4TimeSeries(void)
{
	LIGOTimeGPS epoch = {0, 0};
	REAL8TimeSeries *x;
	REAL4TimeSeries *source = XLALCreateREAL4TimeSeries(NULL, &epoch, 0.0, DELTA_T, &lalDimensionlessUnit, SIMLENGTH);
	REAL4TimeSeries *target = XLALCreateREAL4TimeSeries(NULL, &epoch, 0.0, DELTA_T, &lalDimensionlessUnit, DSTLENGTH);
	double abs_before, abs_after;
	unsigned i;

	/* zero target time series */
	memset(target->data->data, 0, target->data->length * sizeof(*target->data->data));
	/* set "injection" to all 1s */
	for(i = 0; i < source->data->length; i++)
		source->data->data[i] = 1.0;

	x = XLALConvertREAL4TimeSeriesToREAL8(source);
	abs_before = XLALMeasureIntS1S2DT(x, x);
	XLALDestroyREAL8TimeSeries(x);

	/* set injection epoch to something non-trivial */
	XLALGPSAdd(&source->epoch, OFFSET);

	/* add injection to target */
	XLALSimAddInjectionREAL4TimeSeries(target, source, NULL);

	x = XLALConvertREAL4TimeSeriesToREAL8(target);
	abs_after = XLALMeasureIntS1S2DT(x, x);
	XLALDestroyREAL8TimeSeries(x);

#if 0
	/* dump result */
	for(i = 0; i < target->data->length; i++) {
		epoch = target->epoch;
		XLALGPSAdd(&epoch, i * target->deltaT);
		printf("%d.%09d %g\n", epoch.gpsSeconds, epoch.gpsNanoSeconds, target->data->data[i]);
	}
#endif

	fprintf(stderr, "%s(): square integral before = %.17g, square integral after = %.17g, fractional difference = %g\n", __func__, abs_before, abs_after, fabs(abs_after - abs_before) / abs_before);
	return fabs(abs_after - abs_before) / abs_before > REAL4THRESH;
}


static int TestXLALSimAddInjectionREAL8TimeSeries(void)
{
	LIGOTimeGPS epoch = {0, 0};
	REAL8TimeSeries *source = XLALCreateREAL8TimeSeries(NULL, &epoch, 0.0, DELTA_T, &lalDimensionlessUnit, SIMLENGTH);
	REAL8TimeSeries *target = XLALCreateREAL8TimeSeries(NULL, &epoch, 0.0, DELTA_T, &lalDimensionlessUnit, DSTLENGTH);
	double abs_before, abs_after;
	unsigned i;

	/* zero target time series */
	memset(target->data->data, 0, target->data->length * sizeof(*target->data->data));
	/* set "injection" to all 1s */
	for(i = 0; i < source->data->length; i++)
		source->data->data[i] = 1.0;

	abs_before = XLALMeasureIntS1S2DT(source, source);

	/* set injection epoch to something non-trivial */
	XLALGPSAdd(&source->epoch, OFFSET);

	/* add injection to target */
	XLALSimAddInjectionREAL8TimeSeries(target, source, NULL);

	abs_after = XLALMeasureIntS1S2DT(target, target);

#if 0
	/* dump result */
	for(i = 0; i < target->data->length; i++) {
		epoch = target->epoch;
		XLALGPSAdd(&epoch, i * target->deltaT);
		printf("%d.%09d %g\n", epoch.gpsSeconds, epoch.gpsNanoSeconds, target->data->data[i]);
	}
#endif

	fprintf(stderr, "%s(): square integral before = %.17g, square integral after = %.17g, fractional difference = %g\n", __func__, abs_before, abs_after, fabs(abs_after - abs_before) / abs_before);
	return fabs(abs_after - abs_before) / abs_before > REAL8THRESH;
}


int main(int argc, char *argv[])
{
	(void) argc;	/* silence unused parameter warning */
	(void) argv;	/* silence unused parameter warning */
	return TestXLALSimAddInjectionREAL4TimeSeries() || TestXLALSimAddInjectionREAL8TimeSeries();
}
