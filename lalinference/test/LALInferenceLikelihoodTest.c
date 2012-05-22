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

/* functions to check that likelihood tests return cleanly for NULL inputs */
int LALInferenceComputeFrequencyDomainOverlapNullTest(void);
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
