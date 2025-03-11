#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <lal/Date.h>
#include <lal/LALSimulation.h>
#include <lal/LALSimBurst.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#define DELTA_T		(1.0 / 16384)	/* seconds */
#define SOURCE_LENGTH	642.	/* metre */

/* In some error, the hplus wave only takes 0, and this is the test to check hplus surely has some value.  */
static int check_nonzero(REAL8TimeSeries **target){
	unsigned i;

	for(i=0; i<(*target)->data->length; i++){
		if((*target)->data->data[i] != 0){
			return 1;
		}
	}
	return 0;
}

/* hcross must only get the value of 0, and this test return 0 if hcross has nonzero value. */
static int check_zero(REAL8TimeSeries **target){
	unsigned i;

	for(i=0; i<(*target)->data->length; i++){
		if((*target)->data->data[i] != 0){
			return 0;
		}
	}
	return 1;
}


/* If the spectrum of low freauency is too influential to the waveform, the waveform becomes only taking
 * positve value although it unlikes the wave form of Cherenkov Burst. This test return 0 if the wave
 * only gets positive value. */
static int check_positive(REAL8TimeSeries **target){
	unsigned i;

	for(i=0; i<(*target)->data->length; i++){
		if((*target)->data->data[i] < 0){
			return 1;
		}
	}
	return 0;
}

/* Perforiming the tests above. */
static int TestWaveform(void)
{
	REAL8TimeSeries *hplus = NULL, *hcross = NULL;
	double dE_over_dA = 1.;
	int ret = XLALSimBurstCherenkovRadiation(&hplus,&hcross,SOURCE_LENGTH,dE_over_dA,DELTA_T);
	if(ret != XLAL_SUCCESS){
		XLALDestroyREAL8TimeSeries(hplus);
		XLALDestroyREAL8TimeSeries(hcross);
	}
	else if(check_nonzero(&hplus) || check_zero(&hcross) || check_positive(&hplus)){
		fprintf(stderr, "error waveform");
		XLALDestroyREAL8TimeSeries(hplus);
		XLALDestroyREAL8TimeSeries(hcross);

		return 0;
	}
	return 1;

}

int main(void)
{
	return TestWaveform();
}
