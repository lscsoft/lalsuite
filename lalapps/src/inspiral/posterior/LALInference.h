/*

   LALInference:       		Bayesian Followup        
   include/LALInference.h:      main header file


   Copyright 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys, John Veitch

 Test 123
*/

/**
 * \file LALInference.h
 * \brief Main header file
 */

#ifndef LALInference_h
#define LALInference_h

#include <stdio.h>
//...other includes

//Structure containing inference run state 
typedef struct 
tagLALInferenceRunState
{
  ProcParamTable * commandLine;
  LALAlgorithm * algorithm;
  LALEvolveOneStepFunction * evolve;
  LALPriorFunction * prior;
  LALLikelihoodFunction * likelihood;
  LALProposalFunction * proposal;
  LALTemplateFunction * template;
  LALIFOData * data;
  LALVariables * currentParams, *priorArgs, *proposalArgs;
} LALInferenceRunState;


//Main driver function for a run; will distinguish MCMC from NestedSampling
typedef void (Algorithm) (LALInferenceRunState *);

//Compute next state along chain; replaces currentParams
typedef void (LALEvolveOneStepFunction) (LALVariables * currentParams, 
	LALIFOData * data, LALPriorFunction * prior, 
	LALLikelihoodFunction * likelihood, LALProposalFunction * proposal);

//Jump proposal distribution
//Computes proposedParams based on currentParams and additional variables
//stored as proposalArgs, which could include correlation matrix, etc.,
//as well as forward and reverse proposal probability.
//A jump proposal distribution function could call other jump proposal
//distribution functions with various probabilities to allow for multiple
//jump proposal distributions
typedef REAL8 (LALProposalFunction) (LALVariables *currentParams,
	LALVariables *proposedParams, LALVariables *proposalArgs);

typedef REAL8 (LALPriorFunction) (LALVariables *currentParams, 
	LALVariables *priorArgs);

//Likelihood calculator 
//Should take care to perform expensive evaluation of h+ and hx 
//only once if possible, unless necessary because different IFOs 
//have different data lengths or sampling rates 
typedef REAL8 (LALLikelihoodFunction) (LALVariables *currentParams,
        IFOData * data, LALTemplateFunction * template);

//Wrapper for template computation 
//(relies on LAL libraries for implementation) <- could be a #DEFINE ?
typedef void (LALTemplateFunction) (LALVariables *currentParams,
	IFOData *data);

LALInferenceRunState * Initialize (ProcParamsTable * commandLine);

IFOData * ReadData (ProcParamsTable * commandLine);

/*Data storage type definitions*/

typedef enum VariableType {REAL8, REAL4, gslMatrix};  //, ..., ...

//VariableItem should NEVER be accessed directly, only through special
//access functions as defined below.
//Implementation may change from linked list to hashtable for faster access
typedef struct
tagVariableItem
{
  char * name;
  void * value;
  VariableType type;
  VariableItem * next;
}VariableItem;

typedef struct
tagLALVariables
{
  VariableItem * head;
  INT4 dimension;
}LALVariables;

void getVariable(LALVariables * vars, char * name);
void setVariable(LALVariables * vars, char * name, void * value);
void addVariable(LALVariables * vars, char * name, void * value, 
	VariableType type);

typedef struct
tagIFOData
{
  REAL8TimeSeries *timeData, *timeModelhPlus, *timeModelhCross;
  COMPLEX16FrequencySeries *freqData, *freqModelhPlus, *freqModelhCross;
  REAL8FrequencySeries *oneSidedNoisePowerSpectrum;
  REAL8TimeSeries *window;
  REAL8FFTPlan *timeToFreqFFTPlan, *freqToTimeFFTPlan;
  REAL8 fLow, fHigh;	//integration limits;
  LALDetector *detector;
  IFOData *next;
}IFOData;

#endif

