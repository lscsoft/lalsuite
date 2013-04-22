/* 
 *  LALInferenceTest.c:  Unit tests for LALInference.c library code
 *
 *  Copyright (C) 2011 Ben Aylott, Ilya Mandel, Chiara Mingarelli Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch, Will Vousden
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

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALInference.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/Date.h>
#include <lal/XLALError.h>
#include <lal/TimeSeries.h>
#include <lal/Sequence.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/LALInferenceReadData.h> 
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferencePrior.h>

#include "LALInferenceTest.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

size_t LALInferenceTypeSize[] = {sizeof(INT4), 
                                 sizeof(INT8),
                                 sizeof(UINT4),
                                 sizeof(REAL4), 
                                 sizeof(REAL8), 
                                 sizeof(COMPLEX8), 
                                 sizeof(COMPLEX16), 
                                 sizeof(gsl_matrix *),
                                 sizeof(REAL8Vector *),
                                 sizeof(UINT4Vector *),
                                 sizeof(CHAR *),
                                 sizeof(void *)
};

/* for LALInferenceParseCommandLine tests*/
int LALInferenceParseCommandLineTEST_NODBLDASH(void);
int LALInferenceParseCommandLineTEST_STDINPUT(void);
int LALInferenceParseCommandLineTEST_DASHINPUT(void);

/*  LALInferenceProcessParamLine tests */
int LALInferenceProcessParamLine_TEST(void);
int LALInferenceProcessParamLine_TEST_EMPTYFILE(void);
int LALInferenceProcessParamLine_TEST_CHARFILE(void);

/*  LALInferenceExecuteFT tests */
int LALInferenceExecuteFTTEST_NULLPLAN(void);

int main(void){
    lalDebugLevel |= LALERROR;
    
	int failureCount = 0;

	failureCount += LALInferenceParseCommandLineTEST_NODBLDASH();
	printf("\n");
	failureCount += LALInferenceParseCommandLineTEST_STDINPUT();
	printf("\n");
	failureCount += LALInferenceParseCommandLineTEST_DASHINPUT();
	printf("\n");
	failureCount += LALInferenceProcessParamLine_TEST();
	printf("\n");
	failureCount += LALInferenceProcessParamLine_TEST_EMPTYFILE();
	printf("\n");
	failureCount += LALInferenceProcessParamLine_TEST_CHARFILE();
	printf("\n");
	failureCount += LALInferenceExecuteFTTEST_NULLPLAN();
	printf("\n");
	printf("Test results: %i failure(s).\n", failureCount);

	return failureCount;

}

/*****************     TEST CODE for LALInferenceParseCommandLine     *****************/

/*this function tests to see if LALInferenceParseCommandLine catches the double dash condition. Expect fail. */
int LALInferenceParseCommandLineTEST_NODBLDASH(void){
    TEST_HEADER();
    int errnum,i;
    const int number=3;
    ProcessParamsTable *answer;
    
    char** list=(char**)XLALCalloc(3,sizeof(char*));
    for(i=0;i<number;i++){
        list[i]=(char*)XLALCalloc(10,sizeof(char));
    }
    
    strcpy(list[0],"foo");
    strcpy(list[1],"bar");
    strcpy(list[2],"baz");
   
    XLAL_TRY(answer=LALInferenceParseCommandLine(number,list), errnum);
    (void)number;
    
    if (errnum == XLAL_SUCCESS||answer!=NULL)
    {
        TEST_FAIL("Did not use double dash on fist input, should fail; XLAL error code: %i.", errnum);
    } 

    TEST_FOOTER();

}

/*this function tests to see if LALInferenceParseCommandLine tests standard input. Expect pass. */
int LALInferenceParseCommandLineTEST_STDINPUT(void){
    TEST_HEADER();
    int errnum,i;
    ProcessParamsTable *answer;
    const int number=3;
    
    char** list=(char**)XLALCalloc(3,sizeof(char*));
    for(i=0;i<number;i++){
        list[i]=(char*)XLALCalloc(10,sizeof(char));
    }
    
    strcpy(list[0],"filename");
    strcpy(list[1],"--bar");
    strcpy(list[2],"--baz");
   
    XLAL_TRY(answer=LALInferenceParseCommandLine(number,list), errnum);
    if (errnum == XLAL_SUCCESS||answer!=NULL)
    {
        printf("Standard input woking ... \n ");
    } 
    else{
      TEST_FAIL("Input error; XLAL error code: %i.", errnum);
    }

    TEST_FOOTER();

}

/*this function tests to see if LALInferenceParseCommandLine tests inputs. We expect this to fail bc !dbldash. */
int LALInferenceParseCommandLineTEST_DASHINPUT(void){
    TEST_HEADER();
    int errnum,i;
    ProcessParamsTable *answer;
    const int number=3;
    
    char** list=(char**)XLALCalloc(3,sizeof(char*));
    for(i=0;i<number;i++){
        list[i]=(char*)XLALCalloc(10,sizeof(char));
    }
    
    strcpy(list[0],"foo");
    strcpy(list[1],"-bar");
    strcpy(list[2],"-baz");
   
    XLAL_TRY(answer=LALInferenceParseCommandLine(number,list), errnum);
    if (errnum == XLAL_SUCCESS||answer!=NULL)
    {
        TEST_FAIL("Did not use double dash on fist input, should fail; XLAL error: %s.", XLALErrorString(errnum));
    }

    TEST_FOOTER();

}




/*****************     TEST CODE for LALInferenceProcessParamLine     *****************/

//this should work, and it does
int LALInferenceProcessParamLine_TEST(void){
    TEST_HEADER();
    int errnum,i;

    FILE *file;
    LALInferenceVariables *vars;

    const int argc=3;
    const int argLength=10;
    char **headers=XLALCalloc(argc + 1,sizeof(char*));

    for(i=0;i<argc;i++){
        headers[i]=XLALCalloc(argLength,sizeof(char));
    }

    strcpy(headers[0],"foo");
    strcpy(headers[1],"barbie");
    strcpy(headers[2],"baz");
    headers[3]=NULL;
    vars=XLALCalloc(1,sizeof(LALInferenceVariables));

    // Create a temporary file, populate it, and rewind to the beginning before using it.
    file=tmpfile();
    fprintf(file, "-1.01 3.01 4.01");
    rewind(file);

    XLAL_TRY(LALInferenceProcessParamLine(file, headers,vars), errnum);
    if (errnum != XLAL_SUCCESS)
    {
        TEST_FAIL("Could not read file; XLAL error: %s.", XLALErrorString(errnum));
    }

    fclose(file);
    TEST_FOOTER();
}

//what happens with an empty file? expect failure and fails.
int LALInferenceProcessParamLine_TEST_EMPTYFILE(void){
    TEST_HEADER();
    int errnum,i;

    FILE *file;
    LALInferenceVariables *vars;

    const int argc=3;
    const int argLength=10;
    char **headers=XLALCalloc(argc + 1,sizeof(char*));

    for(i=0;i<argc;i++){
        headers[i]=XLALCalloc(argLength,sizeof(char));
    }

    strcpy(headers[0],"foo");
    strcpy(headers[1],"barbie");
    strcpy(headers[2],"baz");
    headers[3]=NULL;
    vars=XLALCalloc(1,sizeof(LALInferenceVariables));

    // Create a temporary file and leave it empty.
    file=tmpfile();

    XLAL_TRY(LALInferenceProcessParamLine(file, headers,vars), errnum);
    if (errnum == XLAL_SUCCESS)
    {
        TEST_FAIL("Reading empty file succeeded but should have failed!.");
    }

    fclose(file);
    TEST_FOOTER();
}

/*what happens if it gets chars instad of doubles? */
int LALInferenceProcessParamLine_TEST_CHARFILE(void){
    TEST_HEADER();
    int errnum,i;

    FILE *file;
    LALInferenceVariables *vars;
    
    char** headers=(char**)XLALCalloc(4,sizeof(char*));
    
    for(i=0;i<3;i++){
        headers[i]=(char*)XLALCalloc(10,sizeof(char));
    }
    
    strcpy(headers[0],"foo");
    strcpy(headers[1],"barbie");
    strcpy(headers[2],"baz");
    headers[3]=NULL;
    vars=XLALCalloc(1,sizeof(LALInferenceVariables));
    
    // Create a temporary file, populate it, and rewind to the beginning before using it.
    file=tmpfile();
    fprintf(file, "2.99 3.01 b");
    rewind(file);

    XLAL_TRY(LALInferenceProcessParamLine(file, headers,vars), errnum);
    if (errnum == XLAL_SUCCESS)
    {
        TEST_FAIL("Should not pass; non-numeric characters in file!.");
    }

    TEST_FOOTER();
}




/*****************     TEST CODE for LALInferenceExecuteFT     *****************/
/* Test that LALInferenceExecuteFT fails if the FFT plan is NULL .*/

int LALInferenceExecuteFTTEST_NULLPLAN(void){
    
    TEST_HEADER();
    
    UINT4 i, length;
    REAL8 deltaF;
    LIGOTimeGPS epoch;
    int errnum;
    
    length = 1;
    
    deltaF=0.1;
    
    LALInferenceIFOData *testIFOData=(LALInferenceIFOData*)LALCalloc(1, sizeof(LALInferenceIFOData));
    //LALInferenceIFOData  *testNULLIFOData = NULL;
    
    REAL8TimeSeries *timeModelhPlus=(REAL8TimeSeries*)XLALCalloc(1, sizeof(REAL8TimeSeries));
    REAL8TimeSeries *timeModelhCross=(REAL8TimeSeries*)XLALCalloc(1, sizeof(REAL8TimeSeries));
    REAL8TimeSeries *TData=(REAL8TimeSeries*)XLALCreateREAL8TimeSeries("timeData",&epoch,0.0,0.1,&lalDimensionlessUnit,length);
    
    REAL8Sequence *TimedataPlus=XLALCreateREAL8Sequence(length);
    REAL8Sequence *TimedataCross=XLALCreateREAL8Sequence(length);
    
    REAL8Window *window=XLALCalloc(1, sizeof(REAL8Window));
    REAL8Sequence *Windowdata=XLALCreateREAL8Sequence(length);
    
    COMPLEX16Sequence *Freqdata=XLALCreateCOMPLEX16Sequence(length);
    
    COMPLEX16FrequencySeries *freqData=XLALCalloc(1, sizeof(COMPLEX16FrequencySeries));
    
    testIFOData->freqData=freqData;
    testIFOData->freqData->deltaF=deltaF;
    testIFOData->freqData->data=Freqdata;
    
    testIFOData->freqData->data->length=length;
    
    testIFOData->timeModelhPlus=timeModelhPlus;
    testIFOData->timeModelhPlus->data=TimedataPlus;
    testIFOData->timeModelhPlus->data->length=length;
    
    testIFOData->timeModelhCross=timeModelhCross;
    testIFOData->timeModelhCross->data=TimedataCross;
    testIFOData->timeModelhCross->data->length=length;
    
    testIFOData->timeData=TData;
    testIFOData->timeData->epoch=epoch;
    
    testIFOData->window=window;
    
    testIFOData->window->data=Windowdata;
    
    testIFOData->timeToFreqFFTPlan=NULL;

    for (i=0; i<length; ++i){
        testIFOData->timeModelhPlus->data->data[i]  = 0.0; 
        testIFOData->timeModelhCross->data->data[i] = 0.0;
    }

    if (testIFOData!=NULL){            
        XLAL_TRY(LALInferenceExecuteFT(testIFOData), errnum);
        
        if( errnum == XLAL_SUCCESS ){
          TEST_FAIL("Should not pass; timeToFreqFFTPlan is NULL!");
        }
    }
    
    TEST_FOOTER();

}


/******************************************
 * 
 * Old tests
 * 
 ******************************************/

LALInferenceVariables variables;
LALInferenceVariables variables2;
LALInferenceVariables currentParams;
LALInferenceIFOData *IfoPtr;
REAL4 number,five;
REAL8 numberR8;
INT4 numberI4;
INT8 numberI8;
COMPLEX8 numberC8;
COMPLEX16 numberC16;
REAL8 likelihood, nulllikelihood;

LALStatus status;	
ProcessParamsTable *ppt, *ptr;
LALInferenceRunState *runstate=NULL;
int i, j, k;

LALInferenceRunState *initialize(ProcessParamsTable *commandLine);

//Test LALEvolveOneStepFunction
void BasicMCMCOneStep(LALInferenceRunState *runState);

//Test LALAlgorithm
void MCMCAlgorithm (struct tagLALInferenceRunState *runState);
void NelderMeadEval(struct tagLALInferenceRunState *runState,
                    char **names, REAL8 *values, int dim,
                    REAL8 *logprior, REAL8 *loglikelihood);
void NelderMeadAlgorithm(struct tagLALInferenceRunState *runState, LALInferenceVariables *subset);


void LALVariablesTest(void);
void DataTest(void);
void TemplateStatPhaseTest(void);
void SingleIFOLikelihoodTest(void);
void BasicMCMCTest(void);
void TemplateDumpTest(void);
void PTMCMCTest(void);

// gsl_rng * InitializeRandomSeed(void);
// unsigned long int random_seed();

//REAL8 NullLogLikelihood(LALInferenceIFOData *data)
///*Idential to FreqDomainNullLogLikelihood                        */
//{
//	REAL8 loglikeli, totalChiSquared=0.0;
//	LALInferenceIFOData *ifoPtr=data;
//	
//	/* loop over data (different interferometers): */
//	while (ifoPtr != NULL) {
//		totalChiSquared+=ComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data);
//		ifoPtr = ifoPtr->next;
//	}
//	loglikeli = -0.5 * totalChiSquared; // note (again): the log-likelihood is unnormalised!
//	return(loglikeli);
//}


LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
	LALInferenceRunState *irs=NULL;
	LALInferenceIFOData *ifoPtr, *ifoListStart;
	//ProcessParamsTable *ppt=NULL;
	unsigned long int randomseed;
	struct timeval tv;
	FILE *devrandom;
	
	irs = XLALCalloc(1, sizeof(LALInferenceRunState));
	/* read data from files: */
	fprintf(stdout, " LALInferenceReadData(): started.\n");
	irs->data = LALInferenceReadData(commandLine);
	/* (this will already initialise each LALInferenceIFOData's following elements:  */
	/*     fLow, fHigh, detector, timeToFreqFFTPlan, freqToTimeFFTPlan,     */
	/*     window, oneSidedNoisePowerSpectrum, timeDate, freqData         ) */
	fprintf(stdout, " LALInferenceReadData(): finished.\n");
	if (irs->data != NULL) {
		fprintf(stdout, " initialize(): successfully read data.\n");
		
		fprintf(stdout, " LALInferenceInjectInspiralSignal(): started.\n");
		LALInferenceInjectInspiralSignal(irs->data,commandLine);
		fprintf(stdout, " LALInferenceInjectInspiralSignal(): finished.\n");
		
		ifoPtr = irs->data;
		ifoListStart = irs->data;
		while (ifoPtr != NULL) {
			/*If two IFOs have the same sampling rate, they should have the same timeModelh*,
			 freqModelh*, and modelParams variables to avoid excess computation 
			 in model waveform generation in the future*/
			LALInferenceIFOData * ifoPtrCompare=ifoListStart;
			int foundIFOwithSameSampleRate=0;
			while(ifoPtrCompare != NULL && ifoPtrCompare!=ifoPtr) {
				if(ifoPtrCompare->timeData->deltaT == ifoPtr->timeData->deltaT){
					ifoPtr->timeModelhPlus=ifoPtrCompare->timeModelhPlus;
					ifoPtr->freqModelhPlus=ifoPtrCompare->freqModelhPlus;
					ifoPtr->timeModelhCross=ifoPtrCompare->timeModelhCross;				
					ifoPtr->freqModelhCross=ifoPtrCompare->freqModelhCross;				
					ifoPtr->modelParams=ifoPtrCompare->modelParams;	
					foundIFOwithSameSampleRate=1;	
					break;
				}
			}
			if(!foundIFOwithSameSampleRate){
				ifoPtr->timeModelhPlus  = XLALCreateREAL8TimeSeries("timeModelhPlus",
																	&(ifoPtr->timeData->epoch),
																	0.0,
																	ifoPtr->timeData->deltaT,
																	&lalDimensionlessUnit,
																	ifoPtr->timeData->data->length);
				ifoPtr->timeModelhCross = XLALCreateREAL8TimeSeries("timeModelhCross",
																	&(ifoPtr->timeData->epoch),
																	0.0,
																	ifoPtr->timeData->deltaT,
																	&lalDimensionlessUnit,
																	ifoPtr->timeData->data->length);
				ifoPtr->freqModelhPlus = XLALCreateCOMPLEX16FrequencySeries("freqModelhPlus",
																			&(ifoPtr->freqData->epoch),
																			0.0,
																			ifoPtr->freqData->deltaF,
																			&lalDimensionlessUnit,
																			ifoPtr->freqData->data->length);
				ifoPtr->freqModelhCross = XLALCreateCOMPLEX16FrequencySeries("freqModelhCross",
																			 &(ifoPtr->freqData->epoch),
																			 0.0,
																			 ifoPtr->freqData->deltaF,
																			 &lalDimensionlessUnit,
																			 ifoPtr->freqData->data->length);
				ifoPtr->modelParams = XLALCalloc(1, sizeof(LALInferenceVariables));
			}
			ifoPtr = ifoPtr->next;
		}
		irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
		printf("Injection Null Log Likelihood: %g\n", irs->currentLikelihood);
	}
	else
		fprintf(stdout, " initialize(): no data read.\n");
	
	/* set up GSL random number generator: */
	gsl_rng_env_setup();
	irs->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
	/* (try to) get random seed from command line: */
	ppt = LALInferenceGetProcParamVal(commandLine, "--randomseed");
	if (ppt != NULL)
		randomseed = atoi(ppt->value);
	else { /* otherwise generate "random" random seed: */
		if ((devrandom = fopen("/dev/random","r")) == NULL) {
			gettimeofday(&tv, 0);
			randomseed = tv.tv_sec + tv.tv_usec;
		} 
		else {
			if(!fread(&randomseed, sizeof(randomseed), 1, devrandom))
                exit(1);
			fclose(devrandom);
		}
	}
	fprintf(stdout, " initialize(): random seed: %lu\n", randomseed);
	gsl_rng_set(irs->GSLrandom, randomseed);
	
	return(irs);
}




//Test LALEvolveOneStepFunction
void BasicMCMCOneStep(LALInferenceRunState *runState)
// Metropolis-Hastings sampler.
{
  REAL8 logPriorCurrent, logPriorProposed;
  REAL8 logLikelihoodCurrent, logLikelihoodProposed;
  LALInferenceVariables proposedParams;
  REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
  REAL8 logAcceptanceProbability;

  // current values:
  logPriorCurrent      = runState->prior(runState, runState->currentParams);
  logLikelihoodCurrent = runState->currentLikelihood;

  // generate proposal:
  proposedParams.head = NULL;
  proposedParams.dimension = 0;
  runState->proposal(runState, &proposedParams);
  if (LALInferenceCheckVariable(runstate->proposalArgs, "logProposalRatio"))
    logProposalRatio = *(REAL8*) LALInferenceGetVariable(runstate->proposalArgs, "logProposalRatio");

  // compute prior & likelihood:
  logPriorProposed = runState->prior(runState, &proposedParams);
  if (logPriorProposed > -HUGE_VAL)
    logLikelihoodProposed = runState->likelihood(&proposedParams, runState->data, runState->templt);
  else
    logLikelihoodProposed = -HUGE_VAL;

  // determine acceptance probability:
  logAcceptanceProbability = (logLikelihoodProposed - logLikelihoodCurrent) 
                             + (logPriorProposed - logPriorCurrent)
                             + logProposalRatio;

  // accept/reject:
  if ((logAcceptanceProbability > 0) 
      || (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
    LALInferenceCopyVariables(&proposedParams, runState->currentParams);
    runState->currentLikelihood = logLikelihoodProposed;
  }

  LALInferenceClearVariables(&proposedParams);	
}



//Test LALAlgorithm
void MCMCAlgorithm(struct tagLALInferenceRunState *runState)
{
  //int i;
  REAL8 dummyR8;
  
  printf(" MCMCAlgorithm(); starting parameter values:\n");
  LALInferencePrintVariables(runState->currentParams);
  // initialize starting likelihood value:
  runState->currentLikelihood = runState->likelihood(runstate->currentParams, runState->data, runState->templt);
  // iterate:
  for(i=0; i<100; i++) {
    printf(" MCMC iteration: %d\n", i+1);
    dummyR8 = runState->currentLikelihood;
    runState->evolve(runState);
    if (runState->currentLikelihood != dummyR8) {
      printf(" accepted! new parameter values:\n");
      LALInferencePrintVariables(runState->currentParams);
    }
  }
}


void NelderMeadEval(struct tagLALInferenceRunState *runState,
                    char **names, REAL8 *values, int dim,
                    REAL8 *logprior, REAL8 *loglikelihood)
// Auxiliary function for "NelderMeadAlgorithm()" (see below).
// Evaluates Prior & Likelihood for a given (sub-) set of parameters.
//  /!\  Side effect: alters value of "runState->currentParams" !
{
  //int i;
  // copy over (subset of) values from "value" argument
  // (other parameter values, if any, remain as they are):
  for (i=0; i<dim; ++i)
    LALInferenceSetVariable(runState->currentParams, names[i], &values[i]);
  // evaluate prior & likelihood:
  *logprior = runstate->prior(runstate, runstate->currentParams);
  if (*logprior > -HUGE_VAL)
    *loglikelihood = runState->likelihood(runstate->currentParams, runState->data, runState->templt);
  else
    *loglikelihood = -HUGE_VAL;
  runState->currentLikelihood = *loglikelihood;
  // printf(" x");
  return;
}


void NelderMeadAlgorithm(struct tagLALInferenceRunState *runState, LALInferenceVariables *subset)
/************************************************************************************/
/*  Nelder-Mead (flexible polyhedron search) algorithm                              */
/*  following D. M. Himmelblau (1972): Applied nonlinear programming. McGraw-Hill.  */
/************************************************************************************/
/* Starting values are generated from the "runState->currentParams" value, by       */
/* using repeated draws from "runState->proposal()" function while ensuring         */
/* non-zero prior density for these.                                                */
/* Depending on the "ML" setting (still hard-coded, see below), it will either      */
/* aim for Maximum-Likelihood (ML) or Maximum-A-Posteriori (MAP) values.            */
/* In future, one should be able to specify the subset of parameters to be          */
/* optimized over (since e.g. one wouldn't want to optimize over PN order, which    */
/* may also be part of the parameters. Or one may want to keep sky location fixed). */
/* By now the algorithm can only handle REAL8 parameters.                           */
/************************************************************************************/
/* TO DO:                                                                           */
/*  - use (named) "subset" argument to determine subset of parameters over which to */
/*    optimize. By now simply all "REAL8" values are used.                          */
/*  - allow to specify "ML" option from outside function.                           */
/*  - allow to supply starting simplex?                                             */
/*  - allow to specify stop criteria from outside.                                  */
/*  - get rid of text output.                                                       */
/*  - somehow allow parameters like phase or rightascension to wrap around          */ 
/*    (i.e. let the simplex move across parameter space bounds)                     */
/************************************************************************************/
{
  int ML = 1; // ML==1 --> Maximum-Likelihood (ML);  ML==0 --> Maximum-A-Posteriori (MAP).
  //REAL8 e = sqrt(LAL_REAL8_EPS); // stop criterion
  REAL8 epsilon = 0.001;  // stop criterion 
  int maxiter = 500;      // also stop criterion

  //int i, j;
  LALInferenceVariables param;
  LALInferenceVariables startval;
  char str[VARNAME_MAX];
  int nmDim;            // dimension of (sub-) space to be optimized over.
  char **nameVec=NULL;  // vector of dimensions' names.
  REAL8 *R8Vec=NULL;
  REAL8 *simplex=NULL;  // (d x (d+1)) - matrix containing simplex vertices.
  REAL8 *val_simplex;   // corresponding function values (likelihood or posterior).
  REAL8 logprior, loglikelihood; // dummy variables.
  REAL8 *centroid, *reflected, *expanded, *contracted; // proposed new vertices...
  REAL8 val_reflected, val_expanded, val_contracted; // ...and corresponding function values
  int iteration;
  int terminate=0;
  int mini, maxi;
  
  printf(" NelderMeadAlgorithm(); current parameter values:\n");
  LALInferencePrintVariables(runState->currentParams);
  startval.head=NULL;
  startval.dimension=0;
  LALInferenceCopyVariables(runState->currentParams, &startval);

  // initialize "param":
  param.head=NULL;
  param.dimension=0;
  // "subset" specified? If not, simply gather all REAL8 elements of "currentParams" to optimize over:
  if (subset==NULL) {
    if (runstate->currentParams == NULL) {
      fprintf(stderr," ERROR in NelderMeadAlgorithm(): no \"runstate->currentParams\" vector provided.\n");
      exit(1);
    }
    i = LALInferenceGetVariableDimension(runstate->currentParams);
    if (i==0) {
      fprintf(stderr," ERROR in NelderMeadAlgorithm(): empty \"runstate->currentParams\" vector provided.\n");
      exit(1);
    }
    for (j=1; j<=i; ++j) {  // check "currentParams" entries and copy all REAL( values:
      if (LALInferenceGetVariableTypeByIndex(runstate->currentParams, j) == LALINFERENCE_REAL8_t){
	strcpy(str, LALInferenceGetVariableName(runstate->currentParams, j));
        LALInferenceAddVariable(&param, str, LALInferenceGetVariable(runstate->currentParams, str), LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
      }
    }
  }
  else {
    fprintf(stderr," ERROR in NelderMeadAlgorithm(): \"subset\" feature not yet implemented.\n"); exit(0);
    // TODO: take a (named) "subset" vector of zeroes/ones indicating which variables optimize over and which to keep fixed.
  }

  // Figure out parameter space dimension, names, &c.:
  nmDim   = LALInferenceGetVariableDimension(&param);
  R8Vec   = (REAL8*) XLALMalloc(sizeof(REAL8) * nmDim);
  nameVec = (char**) XLALMalloc(sizeof(char*) * nmDim);
  for (i=0; i<nmDim; ++i) {
    nameVec[i] = (char*) XLALMalloc(sizeof(char) * VARNAME_MAX);
    strcpy(nameVec[i], LALInferenceGetVariableName(&param, i+1));
    R8Vec[i] = *(REAL8*) LALInferenceGetVariable(&param, nameVec[i]);
  }
  simplex       = (REAL8*) XLALMalloc(sizeof(REAL8) * nmDim * (nmDim+1));
  val_simplex   = (REAL8*) XLALMalloc(sizeof(REAL8) * (nmDim+1));
  centroid   = (REAL8*) XLALMalloc(sizeof(REAL8) * nmDim);
  reflected  = (REAL8*) XLALMalloc(sizeof(REAL8) * nmDim);
  expanded   = (REAL8*) XLALMalloc(sizeof(REAL8) * nmDim);
  contracted = (REAL8*) XLALMalloc(sizeof(REAL8) * nmDim);

  // populate simplex;
  // first corner is starting value:
  for (j=0; j<nmDim; ++j)
    simplex[j] = *(REAL8*) LALInferenceGetVariable(&param, nameVec[j]);
  NelderMeadEval(runState, nameVec, &simplex[0], nmDim, &logprior, &loglikelihood);
  if (!(loglikelihood>-HUGE_VAL)) {
    fprintf(stderr," ERROR in NelderMeadAlgorithm(): invalid starting value provided.\n");
    exit(1);
  }
  val_simplex[0] = ML ? loglikelihood : logprior+loglikelihood;
  // remaining corners are drawn from "runState->proposal()" function:
  for (i=1; i<(nmDim+1); ++i) {  // (loop over vertices (except 1st))
    logprior = -HUGE_VAL;
    while (!(logprior > -HUGE_VAL)) {
      // draw a proposal & copy over:
      LALInferenceCopyVariables(&startval, runState->currentParams);
      runState->proposal(runState, &param);
      for (j=0; j<nmDim; ++j)
        simplex[i*nmDim+j] = *(REAL8*) LALInferenceGetVariable(&param, nameVec[j]);
      // compute prior & likelihood:
      NelderMeadEval(runState, nameVec, &simplex[i*nmDim], nmDim, &logprior, &loglikelihood);
      val_simplex[i] = ML ? loglikelihood : logprior+loglikelihood;
    }    
  }
  // determine minimum & maximum in simplex:
  mini = maxi = 0;
  for (i=1; i<(nmDim+1); ++i) {
    if (val_simplex[i] < val_simplex[mini]) mini = i;
    if (val_simplex[i] > val_simplex[maxi]) maxi = i;
  }

  // start actual Nelder-Mead iterations:
  iteration = 0;
  while (!terminate) {
    ++iteration;
    // determine centroid of simplex, excluding the worst (minimal) point:
    for (i=0; i<nmDim; ++i) {      // (loop over parameter dimensions)
      centroid[i] = 0.0;
      for (j=0; j<(nmDim+1); ++j)  // (loop over simplex vertices)
        centroid[i] += (j==mini) ? 0.0 : (1.0/((double)nmDim)) * simplex[j*nmDim+i];
    }
    NelderMeadEval(runState, nameVec, centroid, nmDim, &logprior, &loglikelihood);
    // UNUSED!!: REAL8 val_centroid = ML ? loglikelihood : logprior+loglikelihood;

    // REFLECT:
    for (i=0; i<nmDim; ++i)
      reflected[i] = centroid[i] + 1.0*(centroid[i] - simplex[mini*nmDim + i]);
    NelderMeadEval(runState, nameVec, reflected, nmDim, &logprior, &loglikelihood);
    val_reflected = ML ? loglikelihood : logprior+loglikelihood;
    if (val_reflected > val_simplex[maxi]) { // reflected better than best so far?
      // EXPAND:
      for (i=0; i<nmDim; ++i)
        expanded[i] = centroid[i] + 2.9*(reflected[i] - centroid[i]);
      NelderMeadEval(runState, nameVec, expanded, nmDim, &logprior, &loglikelihood);
      val_expanded = ML ? loglikelihood : logprior+loglikelihood;
      if (val_expanded > val_simplex[maxi]) { // expanded better than best so far?
        for (i=0; i<nmDim; ++i) // adopt expanded
          simplex[mini*nmDim+i] = expanded[i];
        val_simplex[mini] = val_expanded;
      }
      else {
        for (i=0; i<nmDim; ++i) // adopt reflected
          simplex[mini*nmDim+i] = reflected[i];
        val_simplex[mini] = val_reflected;
      }
    }
    else { // (reflected is worse that best so far)
      // check: reflected better than any of current (except worst)?
      j=0;
      for (i=0; i<(nmDim+1); ++i)
        j += ((i!=mini) && (val_reflected > val_simplex[i]));
      if (j>0) {
        for (i=0; i<nmDim; ++i) // adopt reflected
          simplex[mini*nmDim+i] = reflected[i];
        val_simplex[mini] = val_reflected;
      }
      else { // (reflected is either worst or 2nd worst)
        if (val_reflected > val_simplex[mini]) { // if 2nd worst, adopt
          for (i=0; i<nmDim; ++i) // adopt reflected
            simplex[mini*nmDim+i] = reflected[i];
          val_simplex[mini] = val_reflected;
        }
        // either way: CONTRACT:
        for (i=0; i<nmDim; ++i)
          contracted[i] = centroid[i] + 0.5*(simplex[mini*nmDim+i] - centroid[i]);
        NelderMeadEval(runState, nameVec, contracted, nmDim, &logprior, &loglikelihood);
        val_contracted = ML ? loglikelihood : logprior+loglikelihood;
        if (val_contracted > val_simplex[mini]) { // adopt contracted
          for (i=0; i<nmDim; ++i)
            simplex[mini*nmDim+i] = contracted[i];
          val_simplex[mini] = val_contracted;
        }
        else { // contraction didn't help, REDUCE:
          for (i=0; i<(nmDim+1); ++i)  // loop over vertices
            if (i!=maxi) {
              for (j=0; j<nmDim; ++j)  // loop over parameters
                simplex[i*nmDim+j] += 0.5 * (simplex[maxi*nmDim+j] - simplex[i*nmDim+j]);
              NelderMeadEval(runState, nameVec, &simplex[i*nmDim], nmDim, &logprior, &loglikelihood);
              val_simplex[i] = ML ? loglikelihood : logprior+loglikelihood;
	    }
        }
      }
    }
    // re-determine minimum & maximum:
    mini = maxi = 0;
    for (i=1; i<(nmDim+1); ++i) {
      if (val_simplex[i] < val_simplex[mini]) mini = i;
      if (val_simplex[i] > val_simplex[maxi]) maxi = i;
    }
    printf(" iter=%d,  maxi=%f,  range=%f\n", 
           iteration, val_simplex[maxi], val_simplex[maxi]-val_simplex[mini]);
    // termination condition:
    terminate = ((val_simplex[maxi]-val_simplex[mini]<epsilon) || (iteration>=maxiter));
  }
  // copy optimized value over to "runState->currentParams":
  for (j=0; j<nmDim; ++j)
    LALInferenceSetVariable(runState->currentParams, nameVec[j], &simplex[maxi*nmDim+j]);
  runState->currentLikelihood = ML ? val_simplex[maxi] : runState->likelihood(runstate->currentParams, runState->data, runState->templt);

  printf(" NelderMeadAlgorithm(); done.\n");
  LALInferencePrintVariables(runState->currentParams);

  LALInferenceClearVariables(&startval);
  LALInferenceClearVariables(&param);
  XLALFree(R8Vec);
  for (i=0; i<nmDim; ++i) XLALFree(nameVec[i]);
  XLALFree(nameVec);
  XLALFree(simplex);
  XLALFree(val_simplex);
  XLALFree(centroid);
  XLALFree(reflected);
  XLALFree(expanded);
  XLALFree(contracted);
}


void LALVariablesTest(void)
{
  number = 10.0;
  five=5.0;
  variables.head=NULL;
  variables.dimension=0;
	
  memset(&status,0,sizeof(status));
  LALInferenceAddVariable(&variables, "number", &number, LALINFERENCE_REAL4_t,LALINFERENCE_PARAM_FIXED);
  numberR8 = 7.0;
  LALInferenceAddVariable(&variables, "seven", &numberR8, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  numberR8 = LAL_PI;
  LALInferenceAddVariable(&variables, "pi", &numberR8, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  numberI4 = 123;
  LALInferenceAddVariable(&variables, "small", &numberI4, LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
  numberI8 = 256*256*256*64;
  LALInferenceAddVariable(&variables, "large", &numberI8, LALINFERENCE_INT8_t,LALINFERENCE_PARAM_FIXED);
  numberC8 = crectf( 2.0, 3.0 );
  LALInferenceAddVariable(&variables, "complex1", &numberC8, LALINFERENCE_COMPLEX8_t,LALINFERENCE_PARAM_FIXED);
  numberC16 = crect( 1.23, -3.45 );
  LALInferenceAddVariable(&variables, "complex2", &numberC16, LALINFERENCE_COMPLEX16_t,LALINFERENCE_PARAM_FIXED);

  number=*(REAL4 *)LALInferenceGetVariable(&variables,"number");
  fprintf(stdout,"Got %lf\n",number);
  LALInferenceSetVariable(&variables,"number",&five);
  number=*(REAL4 *)LALInferenceGetVariable(&variables,"number");
  fprintf(stdout,"Got %lf\n",number);
  fprintf(stdout,"Checkvariable?: %i\n",LALInferenceCheckVariable(&variables,"number"));
  LALInferencePrintVariables(&variables);
  LALInferenceCopyVariables(&variables, &variables2);
  LALInferencePrintVariables(&variables2);
  fprintf(stdout,"LALInferenceCompareVariables?: %i\n",
          LALInferenceCompareVariables(&variables,&variables2));
  numberC16 = crect( creal(numberC16), 4.56 );
  LALInferenceSetVariable(&variables2,"complex2",&numberC16);
  fprintf(stdout,"LALInferenceCompareVariables?: %i\n",
          LALInferenceCompareVariables(&variables,&variables2));
  numberC16 = crect( creal(numberC16), -3.45 );
  LALInferenceSetVariable(&variables2,"complex2",&numberC16);
  fprintf(stdout,"LALInferenceCompareVariables?: %i\n",
          LALInferenceCompareVariables(&variables,&variables2));

  LALInferenceRemoveVariable(&variables,"number");
  fprintf(stdout,"Removed, Checkvariable?: %i\n",LALInferenceCheckVariable(&variables,"number"));
  
  fprintf(stdout,"LALInferenceCompareVariables?: %i\n",
          LALInferenceCompareVariables(&variables,&variables2));
  LALInferenceClearVariables(&variables);
  LALInferenceClearVariables(&variables2);
  LALInferencePrintVariables(&variables2);

  fprintf(stdout," ----------\n");
}

void DataTest(void)
{
	fprintf(stdout," data found --> trying some template computations etc.\n");
    
    /* print some information on individual "runstate->data" elements: */
    IfoPtr = runstate->data;  i = 1;
    while (IfoPtr != NULL) {
      if (IfoPtr->timeData)
        fprintf(stdout, " [%d] timeData (\"%s\"): length=%d, deltaT=%f, epoch=%.3f\n", 
                i, IfoPtr->timeData->name, IfoPtr->timeData->data->length, IfoPtr->timeData->deltaT, 
                XLALGPSGetREAL8(&IfoPtr->timeData->epoch));
      if (IfoPtr->freqData)
        fprintf(stdout, "     freqData (\"%s\"): length=%d, deltaF=%f\n", 
                IfoPtr->freqData->name, IfoPtr->freqData->data->length, IfoPtr->freqData->deltaF);
      fprintf(stdout, "     fLow=%.1f Hz,  fHigh=%.1f Hz  (%d freq bins w/in range)\n", 
              IfoPtr->fLow, IfoPtr->fHigh, 
              ((int) (floor(IfoPtr->fHigh / IfoPtr->freqData->deltaF) - ceil(IfoPtr->fLow / IfoPtr->freqData->deltaF)))+1);
      fprintf(stdout, "     detector location: (%.1f, %.1f, %.1f)\n",
              IfoPtr->detector->location[0], IfoPtr->detector->location[1], IfoPtr->detector->location[2]);
      fprintf(stdout, "     detector response matrix:\n");
      for (j=0; j<3; ++j){
        fprintf(stdout, "     ");
        for (k=0; k<3; ++k)
          fprintf(stdout, "%f  ", IfoPtr->detector->response[j][k]);
        fprintf(stdout, "\n");
      }
      IfoPtr = IfoPtr->next;
    }

    IfoPtr=runstate->data;
	SimInspiralTable *injTable=NULL;
	printf("Ninj: %d\n", SimInspiralTableFromLIGOLw(&injTable,LALInferenceGetProcParamVal(ppt,"--injXML")->value,0,0));
	//REAL4 m1 = 10.0;
    //LALInferenceAddVariable(runstate->data->modelParams,"m1",&m1,LALINFERENCE_REAL4_t);
    //REAL4 m2 = 1.4;
    //LALInferenceAddVariable(runstate->data->modelParams,"m2",&m2,LALINFERENCE_REAL4_t);
	//REAL4 inc = 0.0;
    //LALInferenceAddVariable(runstate->data->modelParams,"inc",&inc,LALINFERENCE_REAL4_t);
    //REAL4 phii = 0.0;
    //LALInferenceAddVariable(runstate->data->modelParams,"phii",&phii,LALINFERENCE_REAL4_t);
	//ProcessParamsTable *procparam=LALInferenceGetProcParamVal(ppt,"--trigtime");
	//LIGOTimeGPS trigger_time;
	//char * chartmp;
	//LALStringToGPS(&status,&trigger_time,procparam->value,&chartmp);
	//REAL8 tc = XLALGPSGetREAL8(&trigger_time);
	//LALInferenceAddVariable(runstate->data->modelParams,"time",&tc,LALINFERENCE_REAL8_t);
	
	REAL8 mc = injTable->mchirp;
	REAL8 eta = injTable->eta;
    REAL8 iota = injTable->inclination;
    REAL8 phi = injTable->coa_phase;
	LIGOTimeGPS trigger_time=injTable->geocent_end_time;
	REAL8 tc = XLALGPSGetREAL8(&trigger_time);
	REAL8 ra_current = injTable->longitude;
	REAL8 dec_current = injTable->latitude;
	REAL8 psi_current = injTable->polarization;
	REAL8 distMpc_current = injTable->distance;
	
	
	LALInferenceAddVariable(&currentParams, "chirpmass",       &mc,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "massratio",       &eta,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "inclination",     &iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "phase",           &phi,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(&currentParams, "time",            &tc   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
    LALInferenceAddVariable(&currentParams, "rightascension",  &ra_current,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(&currentParams, "declination",     &dec_current,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "polarisation",    &psi_current,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "distance",        &distMpc_current, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
   /* fprintf(stdout, " trying 'templateLAL' likelihood...\n");
    numberI4 = TaylorF2;
    LALInferenceAddVariable(&currentParams, "LAL_APPROXIMANT", &numberI4,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    numberI4 = LAL_PNORDER_TWO;
    LALInferenceAddVariable(&currentParams, "LAL_PNORDER",     &numberI4,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);*/
	 fprintf(stdout, " trying 'LALTemplateGeneratePPN' likelihood..\n");
    likelihood = LALInferenceFreqDomainLogLikelihood(&currentParams, runstate->data, LALInferenceLALTemplateGeneratePPN);
    nulllikelihood = LALInferenceNullLogLikelihood(runstate->data);
printf("Likelihood %g NullLikelihood %g RelativeLikelihood %g\n", likelihood, nulllikelihood, likelihood-nulllikelihood);

/*
    LALTemplateGeneratePPN(runstate->data);
	  executeFT(runstate->data);
	  
	  FILE *testout=fopen("test_FD.txt","w");
	  for (i=0;i<runstate->data->freqModelhPlus->data->length;i++){
		  fprintf(testout,"%g %g %g %g %g\n",i*runstate->data->freqModelhPlus->deltaF,
				  runstate->data->freqModelhPlus->data->data[i].re,
				  runstate->data->freqModelhPlus->data->data[i].im,
				  runstate->data->freqModelhCross->data->data[i].re,
				  runstate->data->freqModelhCross->data->data[i].im);
	  }
	  fclose(testout);
	  testout=fopen("test_TD.txt","w");
	  for (i=0;i<runstate->data->timeModelhPlus->data->length;i++){
		  fprintf(testout,"%10.10lf %g %g\n",runstate->data->timeData->epoch.gpsSeconds
					+(1e-9*runstate->data->timeData->epoch.gpsNanoSeconds)+i*runstate->data->timeModelhPlus->deltaT,
				  runstate->data->timeModelhPlus->data->data[i],
				  runstate->data->timeModelhCross->data->data[i]);
	  }
	  fclose(testout);
	  testout=fopen("PSD.txt","w");
	  for (i=0;i<runstate->data->oneSidedNoisePowerSpectrum->data->length;i++){
		  fprintf(testout,"%g %g\n",i*runstate->data->oneSidedNoisePowerSpectrum->deltaF,
				  runstate->data->oneSidedNoisePowerSpectrum->data->data[i]);
	  }
	  fclose(testout);
	  testout=fopen("noise_TD.txt","w");
	  for (i=0;i<runstate->data->timeData->data->length;i++){
		  fprintf(testout,"%10.10lf %g\n",runstate->data->timeData->epoch.gpsSeconds
					+(1e-9*runstate->data->timeData->epoch.gpsNanoSeconds)+i*runstate->data->timeData->deltaT,
				  runstate->data->timeData->data->data[i]);
	  }
	  fclose(testout);
	  testout=fopen("noise_FD.txt","w");
	  for (i=0;i<runstate->data->freqData->data->length;i++){
	          //fprintf(testout,"%g %g %g %g %g\n",i*runstate->data->freqData->deltaF,
		  fprintf(testout,"%g %g %g\n",i*runstate->data->freqData->deltaF,
				  runstate->data->freqData->data->data[i].re,
			          //runstate->data->freqData->data->data[i].im,
				  //runstate->data->freqData->data->data[i].re,
				  runstate->data->freqData->data->data[i].im);
	  }
	  
*/	  
    fprintf(stdout," ----------\n");
}


void TemplateStatPhaseTest(void)
{
    fprintf(stdout, " trying out 'templateStatPhase()'...\n");
    REAL8 mc   = 4.0;
    REAL8 eta  = 0.24;
    REAL8 iota = 0.4;
    REAL8 phi  = 2.0;
    REAL8 tcoal   = XLALGPSGetREAL8(&(runstate->data->timeData->epoch)) + 
		(((double)runstate->data->timeData->data->length) * runstate->data->timeData->deltaT) - 1.0;
    printf("TCOAL: %f\n",tcoal);
	REAL8 tc=*((REAL8 *) LALInferenceGetVariable(runstate->data->modelParams,"time"));
	printf("t_c: %f\n", tc);
    LALInferenceClearVariables(runstate->data->modelParams);
    LALInferenceAddVariable(runstate->data->modelParams, "chirpmass",   &mc,    LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(runstate->data->modelParams, "massratio",   &eta,   LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(runstate->data->modelParams, "inclination", &iota,  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(runstate->data->modelParams, "phase",       &phi,   LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(runstate->data->modelParams, "time",        &tcoal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferencePrintVariables(runstate->data->modelParams);
    LALInferenceTemplateStatPhase(runstate->data);
    fprintf(stdout, " ...done.\n");

	  
	  // Parameters for which I am going to compute the likelihood
	  
	  REAL8 ra_current        = 0.0;	/* radian      */
	  REAL8 dec_current       = 0.0;	/* radian      */
	  REAL8 psi_current       = 0.8;	/* radian      */
	  REAL8 distMpc_current   = 10.0;	/* Mpc         */
	  
    LALInferenceAddVariable(&currentParams, "chirpmass",       &mc,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "massratio",       &eta,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "inclination",     &iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "phase",           &phi,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(&currentParams, "time",            &tc   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
    LALInferenceAddVariable(&currentParams, "rightascension",  &ra_current,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(&currentParams, "declination",     &dec_current,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "polarisation",    &psi_current,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "distance",        &distMpc_current, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    fprintf(stdout, " trying 'templateLAL' likelihood...\n");
    numberI4 = TaylorT1;
    LALInferenceAddVariable(&currentParams, "LAL_APPROXIMANT", &numberI4,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    numberI4 = LAL_PNORDER_TWO;
    LALInferenceAddVariable(&currentParams, "LAL_PNORDER",     &numberI4,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    likelihood = LALInferenceFreqDomainLogLikelihood(&currentParams, runstate->data, LALInferenceTemplateLAL);
    fprintf(stdout, " ...done.\n");
    fprintf(stdout," templateLAL log-likelihood %f\n", likelihood);  
    fprintf(stdout," ----------\n");
}


void SingleIFOLikelihoodTest(void)
{
	fprintf(stdout, "Single IFO likelihood test\n");
	COMPLEX16Vector *freqModel1=XLALCreateCOMPLEX16Vector(runstate->data->freqData->data->length);
	COMPLEX16Vector *freqModel2=XLALCreateCOMPLEX16Vector(runstate->data->freqData->data->length);
	numberI4 = LAL_PNORDER_TWO;
    LALInferenceSetVariable(&currentParams, "LAL_PNORDER",     &numberI4);	
	numberI4 = TaylorT1;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);														  																  
	LALInferenceComputeFreqDomainResponse(&currentParams, runstate->data, LALInferenceTemplateLAL, freqModel1);
	freqModel2=runstate->data->freqData->data;
	//ComputeFreqDomainResponse(&currentParams, runstate->data, templateLAL, freqModel2);
	FILE * freqModelFile=fopen("freqModelFile.dat", "w");
	for(i=0; i<(int)runstate->data->freqData->data->length; i++){
		fprintf(freqModelFile, "%g\t %g\t %g\t %g\t %g\t %g\n", 
		((double)i)*1.0/ (((double)runstate->data->timeData->data->length) * runstate->data->timeData->deltaT),
		creal(freqModel1->data[i]), cimag(freqModel1->data[i]), creal(freqModel2->data[i]), cimag(freqModel2->data[i]),
		runstate->data->oneSidedNoisePowerSpectrum->data->data[i]);
	}
	fprintf(stdout, "overlap=%g\n", 
		LALInferenceComputeFrequencyDomainOverlap(runstate->data, freqModel1, freqModel2));
	fprintf(stdout, "<d|d>=%g, <d|h>=%g, <h|h>=%g, <d|h>-1/2<h|h>=%g\n", 
		LALInferenceComputeFrequencyDomainOverlap(runstate->data, freqModel2, freqModel2),
		LALInferenceComputeFrequencyDomainOverlap(runstate->data, freqModel1, freqModel2),
		LALInferenceComputeFrequencyDomainOverlap(runstate->data, freqModel1, freqModel1),
		LALInferenceComputeFrequencyDomainOverlap(runstate->data, freqModel2, freqModel1)
			-0.5*LALInferenceComputeFrequencyDomainOverlap(runstate->data, freqModel1, freqModel1)
		);				
	fprintf(stdout, "likelihood %g\n",
		LALInferenceFreqDomainLogLikelihood(&currentParams, runstate->data, LALInferenceTemplateLAL));
	fprintf(stdout, "undecomposed likelihood %g \n", 
		LALInferenceUndecomposedFreqDomainLogLikelihood(&currentParams, runstate->data, LALInferenceTemplateLAL));
	fprintf(stdout, "null likelihood %g decomposed null likelihood %g\n",
		LALInferenceNullLogLikelihood(runstate->data),
		LALInferenceNullLogLikelihood(runstate->data));
    XLALDestroyCOMPLEX16Vector(freqModel1);
    //	XLALDestroyCOMPLEX16Vector(freqModel2);
}

//void BasicMCMCTest(void)
//{
	//fprintf(stdout, "Try MCMC basic Sampler test\n");
	//runstate->algorithm=MCMCAlgorithm;
	//runstate->evolve=BasicMCMCOneStep;
	//runstate->prior=BasicUniformLALPrior;
	//runstate->proposal=BasicMCMCLALProposal;
        //runstate->proposalArgs = XLALMalloc(sizeof(LALInferenceVariables));
        //runstate->proposalArgs->head=NULL;
        //runstate->proposalArgs->dimension=0;
	//runstate->likelihood=LALInferenceFreqDomainLogLikelihood;
	////runstate->templt=templateLAL;
	//runstate->templt=LALInferenceTemplateStatPhase;
	//runstate->currentParams=&currentParams;
	//MCMCAlgorithm(runstate);
	//fprintf(stdout, "End of MCMC basic Sampler test\n");
//}


void TemplateDumpTest(void)
{
 /* NOTE: try out the "forceTimeLocation" flag within the "templateLAL()" function */
    /*       for aligning (time domain) templates.                                    */
    fprintf(stdout," generating templates & writing to files...:\n");
    LALInferenceDumptemplateFreqDomain(&currentParams, runstate->data, LALInferenceTemplateStatPhase, "test_FTemplate25SP.csv");
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateStatPhase, "test_TTemplate25SP.csv");
    LALInferenceDumptemplateFreqDomain(&currentParams, runstate->data, LALInferenceTemplate3525TD, "test_FTemplate3525TD.csv");
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplate3525TD, "test_TTemplate3525TD.csv");

    fprintf(stdout," ----------\n");
	 
	 double mass1=10.;
	 double mass2=1.4;
    LALInferenceAddVariable(&currentParams, "m1",       &mass1,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddVariable(&currentParams, "m2",       &mass2,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	  double spin1x = 0.5;
	  double spin1y = 0.1;
	  double spin1z = 0.0;
	  double spin2x = 0.2;
	  double spin2y = 0.0;
	  double spin2z = 0.3;
	  LALInferenceAddVariable(&currentParams, "spin1x",       &spin1x,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	  LALInferenceAddVariable(&currentParams, "spin1y",       &spin1y,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	  LALInferenceAddVariable(&currentParams, "spin1z",       &spin1z,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	  LALInferenceAddVariable(&currentParams, "spin2x",       &spin2x,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	  LALInferenceAddVariable(&currentParams, "spin2y",       &spin2y,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	  LALInferenceAddVariable(&currentParams, "spin2z",       &spin2z,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	  double shift0 = 0.3;
	  LALInferenceAddVariable(&currentParams, "shift0",       &shift0,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	  double coa_phase = 0.1;
	  LALInferenceAddVariable(&currentParams, "coa_phase",    &coa_phase,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);	  
	  double PNorder = 3.5;
	  LALInferenceAddVariable(&currentParams, "PNorder",      &PNorder,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);	  
	  LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLALGenerateInspiral, "test_TTemplateLALSTPN.csv");

	  
    /* These are the LAL templates that (...seem to...) work right now: */
    /* TaylorT1, TaylorT2, TaylorT3, TaylorF2, IMRPhenomA, PadeT1, EOB  */
    numberI4 = LAL_PNORDER_TWO;
    LALInferenceSetVariable(&currentParams, "LAL_PNORDER",     &numberI4);
    numberI4 = TaylorF2;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_TTemplateLAL-TF2.csv");
    numberI4 = TaylorT1;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_TTemplateLAL-TT1.csv");
    numberI4 = TaylorT2;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_TTemplateLAL-TT2.csv");
    numberI4 = TaylorT3;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_TTemplateLAL-TT3.csv");

    numberI4 = IMRPhenomA;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_TTemplateLAL-Phenom.csv");
    LALInferenceDumptemplateFreqDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_FTemplateLAL-Phenom.csv");

    numberI4 = PadeT1;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_TTemplateLAL-PadeT1.csv");

    numberI4 = EOB;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    numberI4 = LAL_PNORDER_PSEUDO_FOUR;
    LALInferenceSetVariable(&currentParams, "LAL_PNORDER", &numberI4);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_TTemplateLAL-EOB.csv");

    numberI4 = BCV;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    numberI4 = LAL_PNORDER_TWO;
    LALInferenceSetVariable(&currentParams, "LAL_PNORDER",     &numberI4);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_TTemplateLAL-BCV.csv");

    numberI4 = EOBNR;
    LALInferenceSetVariable(&currentParams, "LAL_APPROXIMANT", &numberI4);
    numberI4 = LAL_PNORDER_PSEUDO_FOUR;
    LALInferenceSetVariable(&currentParams, "LAL_PNORDER", &numberI4);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateLAL, "test_TTemplateLAL-EOBNR.csv");

    fprintf(stdout," ----------\n");

    numberR8 = 440;
    LALInferenceAddVariable(&currentParams, "frequency", &numberR8, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    numberR8 = 1e-19;
    LALInferenceAddVariable(&currentParams, "amplitude", &numberR8, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateSinc, "test_TTemplateSinc.csv");

    numberR8 = 0.01;
    LALInferenceAddVariable(&currentParams, "sigma", &numberR8, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateSineGaussian, "test_TTemplateSineGauss.csv");

    numberR8 = 0.01;
    LALInferenceAddVariable(&currentParams, "tau", &numberR8, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceDumptemplateTimeDomain(&currentParams, runstate->data, LALInferenceTemplateDampedSinusoid, "test_TTemplateDampedSinus.csv");

    LALInferenceClearVariables(&currentParams);
    fprintf(stdout," ----------\n");
}

void PTMCMCTest(void)
{
	fprintf(stdout, "PTMCMC test\n");
	//runstate->algorithm=PTMCMCAlgorithm;
	//runstate->evolve=PTMCMCOneStep;
	//runstate->prior=PTUniformLALPrior;
	//runstate->prior=PTUniformGaussianPrior;
	//runstate->proposal=PTMCMCLALProposal;
	//runstate->proposal=PTMCMCGaussianProposal;
	runstate->proposalArgs = XLALMalloc(sizeof(LALInferenceVariables));
	runstate->proposalArgs->head=NULL;
	runstate->proposalArgs->dimension=0;
	runstate->likelihood=LALInferenceFreqDomainLogLikelihood;
	//runstate->likelihood=GaussianLikelihood;
	runstate->templt=LALInferenceTemplateLAL;
	
	
	SimInspiralTable *injTable=NULL;
	printf("Ninj: %d\n", SimInspiralTableFromLIGOLw(&injTable,LALInferenceGetProcParamVal(ppt,"--injXML")->value,0,0));
	
	REAL8 mc = injTable->mchirp;
	REAL8 eta = injTable->eta;
    REAL8 iota = injTable->inclination;
    REAL8 phi = injTable->coa_phase;
	LIGOTimeGPS trigger_time=injTable->geocent_end_time;
	REAL8 tc = XLALGPSGetREAL8(&trigger_time);
	REAL8 ra_current = injTable->longitude;
	REAL8 dec_current = injTable->latitude;
	REAL8 psi_current = injTable->polarization;
	REAL8 distMpc_current = injTable->distance;
	
    numberI4 = TaylorF2;
    LALInferenceAddVariable(&currentParams, "LAL_APPROXIMANT", &numberI4,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
    numberI4 = LAL_PNORDER_TWO;
    LALInferenceAddVariable(&currentParams, "LAL_PNORDER",     &numberI4,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
	
	LALInferenceAddVariable(&currentParams, "chirpmass",       &mc,              LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "massratio",       &eta,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(&currentParams, "inclination",     &iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(&currentParams, "phase",           &phi,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(&currentParams, "time",            &tc   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
    LALInferenceAddVariable(&currentParams, "rightascension",  &ra_current,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(&currentParams, "declination",     &dec_current,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(&currentParams, "polarisation",    &psi_current,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceAddVariable(&currentParams, "distance",        &distMpc_current, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	
	
	REAL8 x0 = 0.9;
	LALInferenceAddVariable(&currentParams, "x0", &x0,  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	
	
	
	
	runstate->currentParams=&currentParams;
	//PTMCMCAlgorithm(runstate);
	fprintf(stdout, "End of PTMCMC test\n");
}


