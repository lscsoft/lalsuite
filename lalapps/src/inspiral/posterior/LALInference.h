/*
 *
 *  LALInference:             Bayesian Followup        
 *  include/LALInference.h:   main header file
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
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

/**
 * \file LALInference.h
 * \brief Main header file
 */

#ifndef LALInference_h
#define LALInference_h

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

#define VARNAME_MAX 128

# include <lal/LALStdlib.h>
# include <lal/LALConstants.h>
# include <lal/SimulateCoherentGW.h>
# include <lal/GeneratePPNInspiral.h>
# include <lal/LIGOMetadataTables.h>
# include <lal/LALDatatypes.h>
# include <lal/FindChirp.h>
# include <lal/Window.h>

#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/LALDetectors.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/BinaryPulsarTiming.h>


#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>

//...other includes

struct tagLALInferenceRunState;
struct tagLALIFOData;

/*Data storage type definitions*/

typedef enum {
  INT4_t, 
  INT8_t,
  UINT4_t,
  REAL4_t, 
  REAL8_t, 
  COMPLEX8_t, 
  COMPLEX16_t, 
  gslMatrix_t,
  REAL8Vector_t,
  UINT4Vector_t,
  string_t
} VariableType;

typedef enum {
  timeDomain, 
  frequencyDomain
} LALDomain;

typedef enum {
	PARAM_LINEAR,
	PARAM_CIRCULAR,
	PARAM_FIXED,    /* Never changes */
	PARAM_OUTPUT    /* Changed by the inner code and passed out */
} ParamVaryType;


typedef enum
{
	INFERENCE_TAPER_NONE,
	INFERENCE_TAPER_START,
	INFERENCE_TAPER_END,
	INFERENCE_TAPER_STARTEND,
	INFERENCE_TAPER_NUM_OPTS,
	INFERENCE_RING,
	INFERENCE_SMOOTH
}  LALInferenceApplyTaper;


extern size_t typeSize[];

//VariableItem should NEVER be accessed directly, only through special
//access functions as defined below.
//Implementation may change from linked list to hashtable for faster access
typedef struct
tagVariableItem
{
  char                    name[VARNAME_MAX];
  void                    *value;
  VariableType            type;
  ParamVaryType			  vary;
  struct tagVariableItem  *next;
} LALVariableItem;

typedef struct
tagLALVariables
{
  LALVariableItem * head;
  INT4 dimension;
} LALVariables;

const char *translateInternalToExternalParamName(const char *inName);
int fprintParameterNonFixedHeaders(FILE *out, LALVariables *params);

void *getVariable(LALVariables * vars, const char * name);
INT4 getVariableDimension(LALVariables *vars);
INT4 getVariableDimensionNonFixed(LALVariables *vars);
VariableType getVariableTypeByIndex(LALVariables *vars, int index);
VariableType getVariableType(LALVariables *vars, const char *name);
ParamVaryType getVariableVaryType(LALVariables *vars, const char *name);
char *getVariableName(LALVariables *vars, int index);
void setVariable(LALVariables * vars, const char * name, void * value);
void addVariable(LALVariables * vars, const char * name, void * value, 
	VariableType type, ParamVaryType vary);
void removeVariable(LALVariables *vars,const char *name);
int  checkVariable(LALVariables *vars,const char *name);
void destroyVariables(LALVariables *vars);
void copyVariables(LALVariables *origin, LALVariables *target);
void printVariables(LALVariables *var);
int compareVariables(LALVariables *var1, LALVariables *var2);


//Wrapper for template computation 
//(relies on LAL libraries for implementation) <- could be a #DEFINE ?
//typedef void (LALTemplateFunction) (LALVariables *currentParams, struct tagLALIFOData *data); //Parameter Set is modelParams of LALIFOData
typedef void (LALTemplateFunction) (struct tagLALIFOData *data);


//Jump proposal distribution
//Computes proposedParams based on currentParams and additional variables
//stored as proposalArgs, which could include correlation matrix, etc.,
//as well as forward and reverse proposal probability.
//A jump proposal distribution function could call other jump proposal
//distribution functions with various probabilities to allow for multiple
//jump proposal distributions
typedef void (LALProposalFunction) (struct tagLALInferenceRunState *runState,
	LALVariables *proposedParams);

typedef REAL8 (LALPriorFunction) (struct tagLALInferenceRunState *runState,
	LALVariables *params);

//Likelihood calculator 
//Should take care to perform expensive evaluation of h+ and hx 
//only once if possible, unless necessary because different IFOs 
//have different data lengths or sampling rates 
typedef REAL8 (LALLikelihoodFunction) (LALVariables *currentParams,
        struct tagLALIFOData * data, LALTemplateFunction *template);

//Compute next state along chain; replaces currentParams
typedef void (LALEvolveOneStepFunction) (struct tagLALInferenceRunState *runState);

//Main driver function for a run; will distinguish MCMC from NestedSampling
typedef void (LALAlgorithm) (struct tagLALInferenceRunState *runState);

//Structure containing inference run state 
typedef struct 
tagLALInferenceRunState
{
  ProcessParamsTable        *commandLine;
  LALAlgorithm              *algorithm;
  LALEvolveOneStepFunction  *evolve;
  LALPriorFunction          *prior;
  LALLikelihoodFunction     *likelihood;
  LALProposalFunction       *proposal;
  LALTemplateFunction       *template;
  struct tagLALIFOData      *data;
  LALVariables              *currentParams,
    *priorArgs,
    *proposalArgs,
    *algorithmParams; /* Parameters which control the running of the algorithm */
  LALVariables				**livePoints; /* Array of live points for Nested Sampling */
  LALVariables **differentialPoints;
  size_t differentialPointsLength;
  REAL8						currentLikelihood;
  REAL8                     currentPrior;
  gsl_rng                   *GSLrandom;
} LALInferenceRunState;


LALInferenceRunState *initialize(ProcessParamsTable *commandLine);

struct tagLALIFOData * readData (ProcessParamsTable * commandLine);

void injectSignal(struct tagLALIFOData *IFOdata, ProcessParamsTable *commandLine);

#define DETNAMELEN 256
typedef struct
tagLALIFOData
{
  char                       name[DETNAMELEN];
  REAL8TimeSeries           *timeData, 
                            *timeModelhPlus, *timeModelhCross,
                            *whiteTimeData, *windowedTimeData; /* white is not really white, but over-white. */
  /* Stores the log(L) for the model in presence of data.  These were
     added to allow for individual-detector log(L) output.  The
     convention is that loglikelihood always stores the log(L) for the
     model in freqModel... or timeModel....  When a jump is accepted,
     that value is copied into acceptedloglikelihood, which is the
     quantity that is actually output in the output files. */
  REAL8                      nullloglikelihood, loglikelihood, acceptedloglikelihood; 
  REAL8                      fPlus, fCross;
  REAL8                      timeshift;
  COMPLEX16FrequencySeries  *freqData, 
                            *freqModelhPlus, *freqModelhCross,
                            *whiteFreqData; /* Over-white. */
  COMPLEX16TimeSeries       *compTimeData, *compModelData;
  LIGOTimeGPSVector         *dataTimes;
  LALVariables              *modelParams;
  LALVariables				*dataParams; /* Optional data parameters */
  LALDomain                 modelDomain;
  REAL8FrequencySeries      *oneSidedNoisePowerSpectrum;
  REAL8TimeSeries           *timeDomainNoiseWeights; /* Roughly, InvFFT(1/Noise PSD). */
  REAL8Window               *window;
  REAL8FFTPlan              *timeToFreqFFTPlan, *freqToTimeFFTPlan;
  REAL8                     fLow, fHigh;	//integration limits;
  LALDetector               *detector;
  BarycenterInput           *bary;
  EphemerisData             *ephem;
  LIGOTimeGPS				epoch;

  struct tagLALIFOData      *next;
} LALIFOData;

/* Returns the element of the process params table with "name" */
ProcessParamsTable *getProcParamVal(ProcessParamsTable *procparams,const char *name);

void parseCharacterOptionString(char *input, char **strings[], UINT4 *n);

ProcessParamsTable *parseCommandLine(int argc, char *argv[]);

void printCommandLine(ProcessParamsTable *procparams, char *str);

REAL8 UndecomposedFreqDomainLogLikelihood(LALVariables *currentParams, LALIFOData *data, 
                              LALTemplateFunction *template);

/* For testing purposes, likelihood that returns 0.0 = log(1) every
   time.  Activated with the --zeroLogLike command flag. */
REAL8 ZeroLogLikelihood(LALVariables *currentParams, LALIFOData *data, LALTemplateFunction *template);

REAL8 FreqDomainLogLikelihood(LALVariables *currentParams, LALIFOData * data, 
                              LALTemplateFunction *template);
REAL8 ChiSquareTest(LALVariables *currentParams, LALIFOData * data, 
                              LALTemplateFunction *template);
void ComputeFreqDomainResponse(LALVariables *currentParams, LALIFOData * dataPtr, 
                              LALTemplateFunction *template, COMPLEX16Vector *freqWaveform);	

REAL8 TimeDomainLogLikelihood(LALVariables *currentParams, LALIFOData * data, 
                              LALTemplateFunction *template);
void ComputeTimeDomainResponse(LALVariables *currentParams, LALIFOData * dataPtr, 
                               LALTemplateFunction *template, REAL8TimeSeries *timeWaveform);						  
REAL8 ComputeFrequencyDomainOverlap(LALIFOData * data, 
	COMPLEX16Vector * freqData1, COMPLEX16Vector * freqData2);
void COMPLEX16VectorSubtract(COMPLEX16Vector * out, const COMPLEX16Vector * in1, const COMPLEX16Vector * in2);

REAL8 NullLogLikelihood(LALIFOData *data);
REAL8 TimeDomainNullLogLikelihood(LALIFOData *data);

void dumptemplateFreqDomain(LALVariables *currentParams, LALIFOData * data, 
                            LALTemplateFunction *template, const char *filename);
void dumptemplateTimeDomain(LALVariables *currentParams, LALIFOData * data, 
                            LALTemplateFunction *template, const char *filename);

void executeFT(LALIFOData *IFOdata);
void executeInvFT(LALIFOData *IFOdata);

void die(const char *message);
void LALTemplateGeneratePPN(LALIFOData *IFOdata);
void templateStatPhase(LALIFOData *IFOdata);
void templateNullFreqdomain(LALIFOData *IFOdata);
void templateNullTimedomain(LALIFOData *IFOdata);
void templateLAL(LALIFOData *IFOdata);
void template3525TD(LALIFOData *IFOdata);
void templateSineGaussian(LALIFOData *IFOdata);
void templateDampedSinusoid(LALIFOData *IFOdata);
void templateSinc(LALIFOData *IFOdata);
void templateLALSTPN(LALIFOData *IFOdata);
void templateASinOmegaT(LALIFOData *IFOdata);
void templateLALGenerateInspiral(LALIFOData *IFOdata);

void addMinMaxPrior(LALVariables *priorArgs, const char *name, void *min, void *max, VariableType type);
void getMinMaxPrior(LALVariables *priorArgs, const char *name, void *min, void *max);

LALVariableItem *getItem(LALVariables *vars,const char *name);
LALVariableItem *getItemNr(LALVariables *vars, int index);
void fprintSample(FILE *fp,LALVariables *sample);
void fprintSampleNonFixed(FILE *fp,LALVariables *sample);

void mc2masses(double mc, double eta, double *m1, double *m2);

/* Time-Domain Likelihood Utility Functions. */
void makeWhiteData(LALIFOData *IFOdata);
REAL8 WhitenedTimeDomainOverlap(const REAL8TimeSeries *whitenedData, const REAL8TimeSeries *data);
void PSDToTDW(REAL8TimeSeries *TDW, const REAL8FrequencySeries *PSD, const REAL8FFTPlan *plan,
              const REAL8 fMin, const REAL8 fMax);
UINT4 nextPowerOfTwo(const UINT4 n);
void padREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data);
void padWrappedREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data);
UINT4 LIGOTimeGPSToNearestIndex(const LIGOTimeGPS *time, const REAL8TimeSeries *series);
REAL8 integrateSeriesProduct(const REAL8TimeSeries *s1, const REAL8TimeSeries *s2);
void convolveTimeSeries(REAL8TimeSeries *conv, const REAL8TimeSeries *data, const REAL8TimeSeries *response);
UINT4 NTDWFromNPSD(const UINT4 NPSD);
void wrappedTimeSeriesToLinearTimeSeries(REAL8TimeSeries *linear, const REAL8TimeSeries *wrapped);
void linearTimeSeriesToWrappedTimeSeries(REAL8TimeSeries *wrapped, const REAL8TimeSeries *linear);

REAL8 timeDomainOverlap(const REAL8TimeSeries *TDW, const REAL8TimeSeries *A, const REAL8TimeSeries *B);

/* Differential Evolution and Common Format Posterior File Parsing
   Utilities. */

/* Returns an array of header strings (terminated by NULL). */
char **getHeaderLine(FILE *inp);

/* Turns common-format column names into our internal parameter
   names. */
const char *colNameToParamName(const char *colName);

/* Reads one line from the given file and stores the values there into
   the variable structure, using the given header array to name the
   columns.  Returns 0 on success. */
int processParamLine(FILE *inp, char **headers, LALVariables *vars);

#endif

