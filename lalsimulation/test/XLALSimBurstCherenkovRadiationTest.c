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
#define REAL8THRESH	1e-12
#define LOWESTFREQ	10.
#define HIGHESTFREQ	1.e3
#define LOWESTVELOCITY	0.5
#define HIGHESTVELOCITY	2.0
#define SAMPLEFREQ 300 /* Hz */
#define SAMPLEVELOCITY 1.0

static int check_nonzero(REAL8TimeSeries **target){
	unsigned i;

	for(i=0; i<(*target)->data->length; i++){
		if((*target)->data->data[i] != 0){
			return 1;
		}
	}
	return 0;
}

static int check_zero(REAL8TimeSeries **target){
	unsigned i;

	for(i=0; i<(*target)->data->length; i++){
		if((*target)->data->data[i] != 0){
			return 0;
		}
	}
	return 1;
}

static int check_positive(REAL8TimeSeries **target){
	unsigned i;

	for(i=0; i<(*target)->data->length; i++){
		if((*target)->data->data[i] < 0){
			return 1;
		}
	}
	return 0;
}

static int TestFrequencyAndVelocity(void)
{
	double f;
	double beta;

	for(f = LOWESTFREQ; f < HIGHESTFREQ; f += 1.){
		for(beta = LOWESTVELOCITY; beta < HIGHESTVELOCITY; beta += 0.05){
			REAL8TimeSeries *hplus = NULL, *hcross = NULL;
			double Eover_Rsquared = 1.;
			int ret = XLALSimBurstCherenkovRadiation(&hplus,&hcross,f,beta,Eover_Rsquared,DELTA_T);
			if(ret != XLAL_SUCCESS){
				XLALDestroyREAL8TimeSeries(hplus);
				XLALDestroyREAL8TimeSeries(hcross);
			}
			else if(check_nonzero(&hplus) || check_zero(&hcross) || check_positive(&hplus)){
				fprintf(stderr, "error waveform at f = %f, beta = %f\n",f, beta);

				XLALDestroyREAL8TimeSeries(hplus);
				XLALDestroyREAL8TimeSeries(hcross);

				return 0;
			}
		}
	}
	return 1;

}

static int TestEoverRsquared(void)
{

	double Eover_Rsquared = (double) rand()/RAND_MAX;
	REAL8TimeSeries *hplus = NULL, *hcross = NULL;
	XLALSimBurstCherenkovRadiation(&hplus, &hcross, SAMPLEFREQ, SAMPLEVELOCITY, Eover_Rsquared,DELTA_T);
	double Eafter = XLALMeasureEoverRsquared(hplus, hcross);
	if(fabs(Eafter - Eover_Rsquared) / Eover_Rsquared > REAL8THRESH){
		fprintf(stderr, "Expected energy is not acquired.");
	}

	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);

	return fabs(Eafter - Eover_Rsquared) / Eover_Rsquared > REAL8THRESH;

}

int main(void)
{
	return TestFrequencyAndVelocity() || TestEoverRsquared(); 
}
