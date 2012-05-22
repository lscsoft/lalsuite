/* 
 *  LALInferenceLikelihoodTest.c:  Unit tests for LALInference.c library code
 *
 *  Copyright (C) 2011 Ben Aylott, Ilya Mandel, Chiara Mingarelli Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch, Will Vousden
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <lal/LALInference.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Sequence.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInference.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

int LALInferenceComputeFrequencyDomainOverlapNullTest(void);
int LALInferenceComputeFrequencyDomainOverlapTest(void);
int LALInferenceNullLogLikelihoodNullTest(void);
int LALInferenceTimeDomainNullLogLikelihoodNullTest(void);
int LALInferenceWhitenedTimeDomainOverlapNullTest(void);
int LALInferenceTimeDomainNullLogLikelihoodNullTest(void);
int LALInferenceIntegrateSeriesProductNullTest(void);
int LALInferenceConvolveTimeSeriesNullTest(void);

int LALInferenceComputeFrequencyDomainOverlapNullTest(void){
	REAL8 answer;
	fprintf(stdout, " Testing LALInferenceComputeFrequencyDomainOverlap...\n");
	fprintf(stdout, "...with NULL \n");
	answer=LALInferenceComputeFrequencyDomainOverlap(NULL,NULL, NULL);
	//return 0;
        return answer;
}

int LALInferenceComputeFrequencyDomainOverlapTest(void){
/*
Parameters used in this function are:
	dataPtr->timeData->deltaT
	dataPtr->timeData->data->length
	dataPtr->fLow
	dataPtr->fHigh
	freqData1->data[i].re
	freqData2->data[i].re
	freqData1->data[i].im
	freqData2->data[i].im
	dataPtr->oneSidedNoisePowerSpectrum->data->data[i]
*/
	REAL8 answer;
	LALInferenceIFOData * test_data = XLALMalloc(sizeof(LALInferenceIFOData));
	COMPLEX16Vector * test_freq1 = XLALMalloc(sizeof(COMPLEX16Vector));
	COMPLEX16Vector * test_freq2 = XLALMalloc(sizeof(COMPLEX16Vector));
	COMPLEX16 * test_freq_data1 = XLALMalloc(sizeof(COMPLEX16));
	COMPLEX16 * test_freq_data2 = XLALMalloc(sizeof(COMPLEX16));
	/*REAL8TimeSeries * test_data_time_data = XLALMalloc(sizeof(REAL8TimeSeries));*/
	/*REAL8Sequence * test_data_time_data_data = XLALMalloc(sizeof(REAL8Sequence));*/
	/*REAL8FrequencySeries *test_data_noise = XLALMalloc(sizeof(REAL8FrequencySeries));*/
	/*REAL8Sequence * test_data_noise_data = XLALMalloc(sizeof(REAL8Sequence));*/
	
	REAL8 deltaT = 1;
	UINT4 length = 2;
	
	test_freq1->length=length;
	test_freq2->length=length;
	
	for (UINT4 i=1; i<length; i++){
	    test_freq_data1[i] = CX16rect(1, 1);
	}
	
	for (UINT4 i=1; i<length; i++){
	    test_freq_data2[i] = CX16rect(1, 1);
	}	

	test_freq1->data=test_freq_data1;
	test_freq2->data=test_freq_data2;

	test_data->timeData=XLALMalloc(sizeof(REAL8TimeSeries));
	test_data->timeData->deltaT= deltaT;
	test_data->timeData->data=XLALMalloc(sizeof(REAL8Sequence));
	test_data->timeData->data->length=2;
	test_data->fLow=10;
	test_data->fHigh=1000;	

	fprintf(stdout, " Testing LALInferenceComputeFrequencyDomainOverlap...\n");
	fprintf(stdout, "...with test variables \n");
	answer=LALInferenceComputeFrequencyDomainOverlap(NULL,NULL, NULL);
	//return 0;
        return answer;
}

int LALInferenceNullLogLikelihoodNullTest(void){
	REAL8 answer;
	fprintf(stdout, " Testing LALInferenceNullLogLikelihood...\n");
	fprintf(stdout, "...with Null \n");
	answer=LALInferenceNullLogLikelihood(NULL);
	//return 0;
        return answer;
}


	
	

int LALInferenceWhitenedTimeDomainOverlapNullTest(void){
	REAL8 answer;
	fprintf(stdout, " Testing LALInferenceWhitenedTimeDomainOverlap...\n");
	fprintf(stdout, "...with NULL \n");
	answer=LALInferenceWhitenedTimeDomainOverlap(NULL, NULL);
	//return 0;
        return answer;
}

int LALInferenceTimeDomainNullLogLikelihoodNullTest(void){
	REAL8 answer;
	fprintf(stdout, " Testing LALInferenceTimeDomainNullLogLikelihood...\n");
	fprintf(stdout, "...with NULL \n");
	answer = LALInferenceTimeDomainNullLogLikelihood(NULL);
	//return 0;
        return answer;
}

int LALInferenceIntegrateSeriesProductNullTest(void){
	REAL8 answer;
	fprintf(stdout, " Testing LALInferenceIntegrateSeriesProduct...\n");
	fprintf(stdout, "...with NULL \n");
	answer = LALInferenceIntegrateSeriesProduct(NULL, NULL);
	//return 0;
        return answer;
}

int LALInferenceConvolveTimeSeriesNullTest(){
	fprintf(stdout, " Testing LALInferenceConvolveTimeSeries...\n");
	fprintf(stdout, "...with NULL \n");
 	LALInferenceConvolveTimeSeries(NULL, NULL, NULL);
	return 0;
}

int main(void){

	LALInferenceComputeFrequencyDomainOverlapNullTest();
	LALInferenceComputeFrequencyDomainOverlapTest();
	LALInferenceNullLogLikelihoodNullTest();
	LALInferenceWhitenedTimeDomainOverlapNullTest();
	LALInferenceTimeDomainNullLogLikelihoodNullTest();
	LALInferenceIntegrateSeriesProductNullTest();
	LALInferenceConvolveTimeSeriesNullTest();
	return 0;                                
}                                    

REAL8 FreqDomainNullLogLikelihood(LALInferenceIFOData *data);

REAL8 FreqDomainNullLogLikelihood(LALInferenceIFOData *data)
/* calls the `FreqDomainLogLikelihood()' function in conjunction   */
/* with the `templateNullFreqdomain()' template in order to return */
/* the "Null likelihood" without having to bother specifying       */
/* parameters or template while ensuring computations are exactly  */
/* the same as in usual likelihood calculations.                   */
{
	LALInferenceVariables dummyParams;
	double dummyValue;
	double loglikeli;
	/* set some (basically arbitrary) dummy values for intrinsic parameters */
	/* (these shouldn't make a difference, but need to be present):         */
	dummyParams.head      = NULL;
	dummyParams.dimension = 0;
	dummyValue = 0.5;
	LALInferenceAddVariable(&dummyParams, "rightascension", &dummyValue, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceAddVariable(&dummyParams, "declination",    &dummyValue, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddVariable(&dummyParams, "polarisation",   &dummyValue, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddVariable(&dummyParams, "distance",       &dummyValue, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
	dummyValue = XLALGPSGetREAL8(&data->timeData->epoch) 
	+ (((double) data->timeData->data->length) / 2.0) * data->timeData->deltaT;
	LALInferenceAddVariable(&dummyParams, "time",           &dummyValue, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
	loglikeli = LALInferenceFreqDomainLogLikelihood(&dummyParams, data, &LALInferenceTemplateNullFreqdomain);
	LALInferenceDestroyVariables(&dummyParams);
	return(loglikeli);
}
